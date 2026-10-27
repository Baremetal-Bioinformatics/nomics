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
