//! Tabix (`.tbi`) index and point queries over a BGZF-compressed, sorted,
//! tab-delimited file (htslib's TBI layout: binning index plus 16 kb linear
//! index). Used for dbSNP lookups, which only ever ask "which records start
//! at this position".

use std::fs::File;
use std::io::{self, BufReader, Read};
use std::os::unix::fs::FileExt;
use std::path::Path;

use crate::bgzf::{BgzfError, decode_block, read_raw_block};

fn bgzf_err(e: BgzfError) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, format!("{e:?}"))
}

fn invalid(msg: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg.to_string())
}

/// Decompress a whole BGZF file into memory.
fn read_bgzf(path: &Path) -> io::Result<Vec<u8>> {
    let mut r = BufReader::new(File::open(path)?);
    let (mut raw, mut out, mut coffset) = (Vec::new(), Vec::new(), 0u64);
    while read_raw_block(&mut r, &mut raw, coffset).map_err(bgzf_err)? {
        decode_block(&raw, coffset, &mut out).map_err(bgzf_err)?;
        coffset += raw.len() as u64;
    }
    Ok(out)
}

struct RefIndex {
    /// (bin, chunks as virtual-offset pairs)
    bins: Vec<(u32, Vec<(u64, u64)>)>,
    linear: Vec<u64>,
}

pub struct Tabix {
    names: Vec<Vec<u8>>,
    refs: Vec<RefIndex>,
    /// 0-based column of the sequence name and start position.
    col_seq: usize,
    col_beg: usize,
    meta: u8,
    /// Read with `read_at`, so queries take `&self` and run concurrently.
    file: File,
}

/// `Read` over a shared file at a private offset.
struct At<'a> {
    file: &'a File,
    pos: u64,
}

impl Read for At<'_> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let n = self.file.read_at(buf, self.pos)?;
        self.pos += n as u64;
        Ok(n)
    }
}

struct Cursor<'a> {
    b: &'a [u8],
    at: usize,
}

impl Cursor<'_> {
    fn take(&mut self, n: usize) -> io::Result<&[u8]> {
        let s = self
            .b
            .get(self.at..self.at + n)
            .ok_or_else(|| invalid("truncated tabix index"))?;
        self.at += n;
        Ok(s)
    }
    fn i32(&mut self) -> io::Result<i32> {
        Ok(i32::from_le_bytes(self.take(4)?.try_into().unwrap()))
    }
    fn u32(&mut self) -> io::Result<u32> {
        Ok(u32::from_le_bytes(self.take(4)?.try_into().unwrap()))
    }
    fn u64(&mut self) -> io::Result<u64> {
        Ok(u64::from_le_bytes(self.take(8)?.try_into().unwrap()))
    }
}

/// Bins that may hold records overlapping the 0-based half-open `[beg, end)`.
fn reg2bins(beg: u64, end: u64) -> Vec<u32> {
    let end = end - 1;
    let mut bins = vec![0];
    for (shift, offset) in [(26, 1), (23, 9), (20, 73), (17, 585), (14, 4681)] {
        bins.extend((offset + (beg >> shift))..=(offset + (end >> shift)));
    }
    bins.into_iter().map(|b| b as u32).collect()
}

impl Tabix {
    /// VCF records overlapping 1-based closed `[beg, end]`, in file order.
    /// REF length and INFO/END determine the span. Unlike start-only queries,
    /// this includes earlier indels/SVs that can mask a SNP in FeatureContext.
    pub fn vcf_lines_overlapping(
        &self,
        contig: &[u8],
        beg: u64,
        end: u64,
    ) -> io::Result<Vec<Vec<u8>>> {
        if beg == 0 || end < beg || end > 1 << 29 {
            return Err(invalid("invalid tabix VCF interval"));
        }
        let Some(tid) = self.names.iter().position(|n| n == contig) else {
            return Ok(Vec::new());
        };
        let wanted = reg2bins(beg - 1, end);
        let mut chunks: Vec<_> = self.refs[tid]
            .bins
            .iter()
            .filter(|(b, _)| wanted.contains(b))
            .flat_map(|(_, cs)| cs.iter().copied())
            .collect();
        chunks.sort_unstable();
        let mut merged: Vec<(u64, u64)> = Vec::new();
        for (b, e) in chunks {
            if let Some(last) = merged.last_mut().filter(|last| b <= last.1) {
                last.1 = last.1.max(e);
            } else {
                merged.push((b, e));
            }
        }
        let mut out = Vec::new();
        for (b, e) in merged {
            for line in self.read_range_exclusive(b, e)?.split(|&b| b == b'\n') {
                if line.is_empty() || line[0] == self.meta {
                    continue;
                }
                let record = crate::vcf::parse_line(line)?;
                if record.chrom != contig {
                    continue;
                }
                if record.pos > end {
                    break;
                }
                let mut stop = record.pos + record.ref_allele.len() as u64 - 1;
                for field in record.info.split(|&b| b == b';') {
                    if let Some(value) = field.strip_prefix(b"END=") {
                        stop = std::str::from_utf8(value)
                            .ok()
                            .and_then(|s| s.parse().ok())
                            .ok_or_else(|| invalid("invalid VCF END"))?;
                    }
                }
                if stop >= beg {
                    out.push(line.to_vec());
                }
            }
        }
        Ok(out)
    }

    /// Open `data` with its index at `data` + `.tbi`.
    pub fn open(data: &Path) -> io::Result<Self> {
        let mut tbi = data.as_os_str().to_owned();
        tbi.push(".tbi");
        let idx = read_bgzf(Path::new(&tbi))?;
        let mut c = Cursor { b: &idx, at: 0 };
        if c.take(4)? != b"TBI\x01" {
            return Err(invalid("not a tabix index"));
        }
        let n_ref = c.i32()? as usize;
        let _format = c.i32()?;
        let col_seq = c.i32()? as usize - 1;
        let col_beg = c.i32()? as usize - 1;
        let _col_end = c.i32()?;
        let meta = c.i32()? as u8;
        let _skip = c.i32()?;
        let l_nm = c.i32()? as usize;
        let names = c
            .take(l_nm)?
            .split(|&b| b == 0)
            .filter(|n| !n.is_empty())
            .map(<[u8]>::to_vec)
            .collect::<Vec<_>>();
        let mut refs = Vec::with_capacity(n_ref);
        for _ in 0..n_ref {
            let n_bin = c.i32()? as usize;
            let mut bins = Vec::with_capacity(n_bin);
            for _ in 0..n_bin {
                let bin = c.u32()?;
                let n_chunk = c.i32()? as usize;
                let chunks = (0..n_chunk)
                    .map(|_| Ok((c.u64()?, c.u64()?)))
                    .collect::<io::Result<Vec<_>>>()?;
                bins.push((bin, chunks));
            }
            let n_intv = c.i32()? as usize;
            let linear = (0..n_intv).map(|_| c.u64()).collect::<io::Result<_>>()?;
            refs.push(RefIndex { bins, linear });
        }
        if names.len() != n_ref {
            return Err(invalid("tabix name count mismatch"));
        }
        Ok(Self {
            names,
            refs,
            col_seq,
            col_beg,
            meta,
            file: File::open(data)?,
        })
    }

    /// Lines (without newline) whose sequence is `contig` and whose start
    /// column equals the 1-based `pos`, in file order.
    pub fn lines_starting_at(&self, contig: &[u8], pos: u64) -> io::Result<Vec<Vec<u8>>> {
        self.lines_starting_in(contig, pos, pos)
    }

    /// Lines whose start column lies in the 1-based closed `[beg, end]`, in
    /// file order.
    pub fn lines_starting_in(&self, contig: &[u8], beg: u64, end: u64) -> io::Result<Vec<Vec<u8>>> {
        let Some(tid) = self.names.iter().position(|n| n == contig) else {
            return Ok(Vec::new());
        };
        let r = &self.refs[tid];
        let beg0 = beg - 1;
        let min_off = r
            .linear
            .get((beg0 >> 14) as usize)
            .or(r.linear.last())
            .copied()
            .unwrap_or(0);
        let wanted = reg2bins(beg0, end);
        let mut chunks: Vec<(u64, u64)> = r
            .bins
            .iter()
            .filter(|(b, _)| wanted.contains(b))
            .flat_map(|(_, cs)| cs.iter().copied())
            .filter(|&(_, e)| e > min_off)
            .collect();
        chunks.sort_unstable();
        let mut out = Vec::new();
        let mut seen_past = false;
        for (cbeg, cend) in chunks {
            if seen_past {
                break;
            }
            let cbeg = cbeg.max(min_off);
            let text = self.read_range(cbeg, cend)?;
            for line in text.split(|&b| b == b'\n') {
                if line.is_empty() || line[0] == self.meta {
                    continue;
                }
                let mut f = line.split(|&b| b == b'\t');
                let cols: Vec<&[u8]> = f
                    .by_ref()
                    .take(self.col_seq.max(self.col_beg) + 1)
                    .collect();
                let (Some(seq), Some(start)) = (cols.get(self.col_seq), cols.get(self.col_beg))
                else {
                    continue;
                };
                if *seq != contig {
                    continue;
                }
                let start: u64 = std::str::from_utf8(start)
                    .ok()
                    .and_then(|s| s.parse().ok())
                    .ok_or_else(|| invalid("bad start column"))?;
                if start > end {
                    seen_past = true;
                    break;
                }
                if start >= beg {
                    out.push(line.to_vec());
                }
            }
        }
        Ok(out)
    }

    // TBI chunks end on record boundaries. Keep the end exclusive, including
    // when both offsets are in the same block; do not read the next record.
    fn read_range_exclusive(&self, beg: u64, end: u64) -> io::Result<Vec<u8>> {
        let mut coff = beg >> 16;
        let mut file = BufReader::new(At {
            file: &self.file,
            pos: coff,
        });
        let (mut raw, mut block, mut out) = (Vec::new(), Vec::new(), Vec::new());
        loop {
            if coff > end >> 16 || (coff == end >> 16 && end & 0xffff == 0) {
                break;
            }
            if !read_raw_block(&mut file, &mut raw, coff).map_err(bgzf_err)? {
                break;
            }
            block.clear();
            decode_block(&raw, coff, &mut block).map_err(bgzf_err)?;
            let start = if coff == beg >> 16 {
                (beg & 0xffff) as usize
            } else {
                0
            };
            let stop = if coff == end >> 16 {
                (end & 0xffff) as usize
            } else {
                block.len()
            };
            out.extend_from_slice(
                block
                    .get(start..stop)
                    .ok_or_else(|| invalid("invalid TBI chunk"))?,
            );
            coff += raw.len() as u64;
        }
        Ok(out)
    }

    /// Decompressed bytes from virtual offset `vbeg` up to `vend`, cut to
    /// whole lines (a line straddling `vend` is completed).
    fn read_range(&self, vbeg: u64, vend: u64) -> io::Result<Vec<u8>> {
        let (mut coff, uoff) = (vbeg >> 16, (vbeg & 0xffff) as usize);
        let (cend, uend) = (vend >> 16, (vend & 0xffff) as usize);
        let mut file = BufReader::new(At {
            file: &self.file,
            pos: coff,
        });
        let (mut raw, mut out) = (Vec::new(), Vec::new());
        let mut first = true;
        // Once the chunk's last block is decoded: cut at the first newline at or after this index.
        let mut cut_from: Option<usize> = None;
        while read_raw_block(&mut file, &mut raw, coff).map_err(bgzf_err)? {
            let at = out.len();
            decode_block(&raw, coff, &mut out).map_err(bgzf_err)?;
            if first {
                out.drain(at..at + uoff.min(out.len() - at));
                first = false;
            }
            if coff == cend {
                cut_from = Some((at + uend).min(out.len()));
            }
            coff += raw.len() as u64;
            if let Some(from) = cut_from
                && let Some(nl) = out[from..].iter().position(|&b| b == b'\n')
            {
                out.truncate(from + nl);
                break;
            }
        }
        Ok(out)
    }
}

#[cfg(test)]
mod overlap_tests {
    use super::*;
    use std::io::Write;
    #[test]
    fn vcf_overlap_keeps_earlier_spans_and_duplicate_records_once_each() {
        let p = std::env::temp_dir().join(format!(
            "nomics-tabix-overlap-{}.vcf.gz",
            std::process::id()
        ));
        let mut w = crate::bgzf::BgzfWriter::new(
            File::create(&p).unwrap(),
            crate::bgzf::Stored,
            crate::bgzf::Profile::HTSJDK,
        );
        writeln!(w, "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100000>\n##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();
        for line in [
            "chr1\t1\tlong\tA\t<DEL>\t.\t.\tEND=40000",
            "chr1\t16380\tdel\tAAA\tA\t.\t.\t.",
            "chr1\t16381\tsnp\tA\tC\t.\t.\t.",
            "chr1\t16381\tsnp\tA\tC\t.\t.\t.",
            "chr1\t40001\tnext\tA\tG\t.\t.\t.",
        ] {
            writeln!(w, "{line}").unwrap();
        }
        w.finish().unwrap();
        let mut tbi = p.as_os_str().to_owned();
        tbi.push(".tbi");
        crate::htsjdk_index::write_tbi_htsjdk_vcf(&p, Path::new(&tbi)).unwrap();
        let q = Tabix::open(&p).unwrap();
        let lines = q.vcf_lines_overlapping(b"chr1", 16381, 16381).unwrap();
        assert_eq!(lines.len(), 4); // END span, REF span, both identical SNP records
        let ids: Vec<_> = lines
            .iter()
            .map(|l| crate::vcf::parse_line(l).unwrap().id)
            .collect();
        assert_eq!(ids, [b"long".as_slice(), b"del", b"snp", b"snp"]);
        assert_eq!(
            q.vcf_lines_overlapping(b"chr1", 40000, 40000)
                .unwrap()
                .len(),
            1
        );
        assert_eq!(
            q.vcf_lines_overlapping(b"chr1", 40001, 40001)
                .unwrap()
                .len(),
            1
        );
        assert!(q.vcf_lines_overlapping(b"chr2", 1, 100).unwrap().is_empty());
        assert!(q.vcf_lines_overlapping(b"chr1", 0, 1).is_err());
        std::fs::remove_file(p).unwrap();
        std::fs::remove_file(tbi).unwrap();
    }
}
