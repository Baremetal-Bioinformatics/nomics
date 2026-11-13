//! `.bai` and `.tbi` indexes laid out exactly as htsjdk 4.2.0 writes them
//! while GATK writes the indexed file (`BAMIndexer` for MarkDuplicates /
//! ApplyBQSR `CREATE_INDEX`, `TabixIndexCreator` for VCF output). These
//! differ from `samtools index` / `tabix`: chunks coalesce only within the
//! same or the next block address, htsjdk's linear-index rules apply, and a
//! BAI carries the metadata pseudo-bin and the no-coordinate count.
//!
//! The index is computed from the finished file: every record's virtual
//! offsets are recovered from the BGZF block layout. htsjdk takes them from
//! `BlockCompressedOutputStream.getFilePointer()` while writing, which points
//! at the *next* block as soon as a block has been deflated (full, or
//! flushed). So a position at the end of a block maps to `(next block, 0)`.
//! At the very end of the data the two writers differ: `BAMIndexer.finish`
//! runs while the BAM's last (partial) block is still held, giving
//! `(last block, length)`, whereas the VCF writer flushes the stream before
//! finalizing the index, giving `(next block, 0)`.
//!
//! Sources (htsjdk 4.2.0): `BAMIndexer`, `BinningIndexBuilder`, `Bin.addChunk`,
//! `BAMIndexMetaData`, `BinaryBAMIndexWriter`, `LinearIndex`,
//! `BlockCompressedFilePointerUtil`, `GenomicIndexUtil.regionToBin`,
//! `SAMRecord.computeIndexingBin`, `TabixIndexCreator`, `TabixIndex.write`.
//!
//! Gate: the golden `md.bai`, `recal.bai`, `.vcf.gz.tbi` and `.g.vcf.gz.tbi`
//! are reproduced byte for byte from their golden data files.

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{self, BufReader, Read, Write};
use std::path::Path;

use crate::bgzf::{BgzfError, BgzfWriter, Profile, decode_block, read_raw_block};
use crate::deflate::HtsjdkDeflater;

fn invalid(msg: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg.into())
}

fn bgzf_err(e: BgzfError) -> io::Error {
    e.into()
}

/// `GenomicIndexUtil.MAX_BINS`: the BAI metadata pseudo-bin.
const MAX_BINS: u32 = 37450;
/// htsjdk's full-block payload (`DEFAULT_UNCOMPRESSED_BLOCK_SIZE`).
const FULL_BLOCK: usize = Profile::HTSJDK.block_payload;

/// `GenomicIndexUtil.regionToBin` (0-based, half-open).
pub fn region_to_bin(beg: i32, end: i32) -> u32 {
    let end = end - 1;
    let b = if beg >> 14 == end >> 14 {
        ((1 << 15) - 1) / 7 + (beg >> 14)
    } else if beg >> 17 == end >> 17 {
        ((1 << 12) - 1) / 7 + (beg >> 17)
    } else if beg >> 20 == end >> 20 {
        ((1 << 9) - 1) / 7 + (beg >> 20)
    } else if beg >> 23 == end >> 23 {
        ((1 << 6) - 1) / 7 + (beg >> 23)
    } else if beg >> 26 == end >> 26 {
        1 + (beg >> 26) // ((1 << 3) - 1) / 7
    } else {
        0
    };
    b as u32
}

/// `LinearIndex.convertToLinearIndexOffset`.
fn linear_window(pos: i32) -> usize {
    let p = if pos <= 0 { 0 } else { pos - 1 };
    (p >> 14) as usize
}

fn block_address(vfp: u64) -> u64 {
    (vfp >> 16) & 0xFFFF_FFFF_FFFF
}

type Bins = BTreeMap<u32, Vec<(u64, u64)>>;

/// `BinningIndexBuilder` for one reference (fill-in of uninitialised windows on).
#[derive(Default)]
struct Builder {
    bins: Bins,
    linear: Vec<i64>,
}

impl Builder {
    /// `processFeature`: 1-based inclusive `start..=end`, `end == 0` = unset.
    fn feature(&mut self, start: i32, end: i32, bin: u32, chunk: (u64, u64)) {
        let chunks = self.bins.entry(bin).or_default();
        // Bin.addChunk: coalesce when in the same or the next block address.
        let coalesce = chunks.last().is_some_and(|last| {
            let (a, b) = (block_address(last.1), block_address(chunk.0));
            a == b || a + 1 == b
        });
        if coalesce {
            chunks.last_mut().unwrap().1 = chunk.1;
        } else {
            chunks.push(chunk);
        }
        let (start_w, end_w) = if end == 0 {
            let w = linear_window(start - 1);
            (w, w)
        } else {
            (linear_window(start), linear_window(end))
        };
        if end_w + 1 > self.linear.len() {
            self.linear.resize(end_w + 1, -1);
        }
        for w in start_w..=end_w {
            let v = &mut self.linear[w];
            if *v == -1 || (chunk.0 as i64) < *v {
                *v = chunk.0 as i64;
            }
        }
    }

    /// `generateIndexContent`: `None` without bins; else the bins and the
    /// linear index (as long as the largest window seen), gaps filled with
    /// the last earlier value (0 before any).
    fn finish(mut self) -> Option<(Bins, Vec<u64>)> {
        if self.bins.is_empty() {
            return None;
        }
        let mut last = 0i64;
        for v in &mut self.linear {
            if *v == -1 {
                *v = last;
            } else {
                last = *v;
            }
        }
        Some((self.bins, self.linear.into_iter().map(|v| v as u64).collect()))
    }
}

/// Sequential reader of a BGZF file's payload that knows each byte's
/// virtual offset as htsjdk's writer reported it.
struct VirtualReader<R> {
    r: R,
    raw: Vec<u8>,
    /// Buffered payload; `blocks` holds, per data block in it,
    /// (start in `buf`, coffset, payload length, member length).
    buf: Vec<u8>,
    blocks: Vec<(usize, u64, usize, u64)>,
    next_coffset: u64,
    eof: bool,
}

impl<R: Read> VirtualReader<R> {
    fn new(r: R) -> Self {
        VirtualReader {
            r,
            raw: Vec::new(),
            buf: Vec::new(),
            blocks: Vec::new(),
            next_coffset: 0,
            eof: false,
        }
    }

    /// Buffer until `buf.len() >= n` or the file ends; empty members (the
    /// EOF marker) carry no data and are skipped.
    fn fill(&mut self, n: usize) -> io::Result<()> {
        while self.buf.len() < n && !self.eof {
            if !read_raw_block(&mut self.r, &mut self.raw, self.next_coffset).map_err(bgzf_err)? {
                self.eof = true;
                break;
            }
            let at = self.buf.len();
            let info = decode_block(&self.raw, self.next_coffset, &mut self.buf).map_err(bgzf_err)?;
            if info.isize > 0 {
                self.blocks
                    .push((at, info.coffset, info.isize as usize, info.csize as u64));
            }
            self.next_coffset += info.csize as u64;
        }
        Ok(())
    }

    /// `n` bytes at `pos`, or `None` at the end of the data.
    fn bytes(&mut self, pos: usize, n: usize) -> io::Result<Option<&[u8]>> {
        if n == 0 {
            return Ok(Some(&[]));
        }
        self.fill(pos + n)?;
        if self.buf.len() < pos + n {
            if self.buf.len() == pos || n == 0 {
                return Ok(None);
            }
            return Err(invalid("BGZF data ends inside a record"));
        }
        Ok(Some(&self.buf[pos..pos + n]))
    }

    /// The virtual offset htsjdk reported after writing the first `pos`
    /// bytes of the buffered stream. At the end of the data, `flushed` says
    /// whether the writer had already deflated its last block when it asked.
    fn vo_at(&mut self, pos: usize, flushed: bool) -> io::Result<u64> {
        self.fill(pos + 1)?;
        for &(s, c, len, csize) in self.blocks.iter().rev() {
            if s > pos {
                continue;
            }
            if pos < s + len {
                return Ok(c << 16 | (pos - s) as u64);
            }
            if pos == s + len {
                // End of the data (no later block holds `pos`): htsjdk had
                // deflated this block only if it was full.
                return Ok(if len == FULL_BLOCK || flushed {
                    (c + csize) << 16
                } else {
                    c << 16 | len as u64
                });
            }
            break;
        }
        if pos == 0 && self.blocks.is_empty() {
            return Ok(0);
        }
        Err(invalid("virtual offset beyond the data"))
    }

    /// [`Self::vo_at`] for a position inside the data, or its end when the
    /// writer had not flushed (`BAMIndexer.finish` before the BAM is closed).
    fn vo(&mut self, pos: usize) -> io::Result<u64> {
        self.vo_at(pos, false)
    }

    /// Forget payload before `pos` (whole blocks only); returns the number
    /// of bytes dropped, which callers subtract from their positions.
    fn compact(&mut self, pos: usize) -> usize {
        let Some(k) = self.blocks.iter().rposition(|&(s, _, _, _)| s <= pos) else {
            return 0;
        };
        let cut = self.blocks[k].0;
        if cut == 0 {
            return 0;
        }
        self.buf.drain(..cut);
        self.blocks.drain(..k);
        for b in &mut self.blocks {
            b.0 -= cut;
        }
        cut
    }
}

fn le_i32(b: &[u8], at: usize) -> i32 {
    i32::from_le_bytes(b[at..at + 4].try_into().unwrap())
}

fn le_u16(b: &[u8], at: usize) -> u16 {
    u16::from_le_bytes(b[at..at + 2].try_into().unwrap())
}

/// htsjdk's `.bai` for a coordinate-sorted BAM, as `BAMIndexer` builds it
/// while writing that BAM.
pub fn bai_htsjdk(bam: &Path) -> io::Result<Vec<u8>> {
    let mut vr = VirtualReader::new(BufReader::with_capacity(1 << 20, File::open(bam)?));
    // Header: magic, l_text, text, n_ref, refs.
    let h = vr
        .bytes(0, 8)?
        .ok_or_else(|| invalid("empty BAM"))?
        .to_vec();
    if &h[0..4] != b"BAM\x01" {
        return Err(invalid("not a BAM file"));
    }
    let l_text = le_i32(&h, 4) as usize;
    let mut pos = 8 + l_text;
    let text = vr
        .bytes(8, l_text)?
        .ok_or_else(|| invalid("truncated BAM header"))?;
    let text = String::from_utf8_lossy(text).into_owned();
    let so = text
        .lines()
        .next()
        .filter(|l| l.starts_with("@HD"))
        .and_then(|l| l.split('\t').find_map(|f| f.strip_prefix("SO:")))
        .unwrap_or("unsorted")
        .to_owned();
    if so != "coordinate" && so != "unsorted" {
        return Err(invalid(format!("indexing requires a coordinate-sorted BAM (SO:{so})")));
    }
    let n_ref = le_i32(vr.bytes(pos, 4)?.ok_or_else(|| invalid("truncated BAM header"))?, 0);
    pos += 4;
    for _ in 0..n_ref {
        let l_name = le_i32(vr.bytes(pos, 4)?.ok_or_else(|| invalid("truncated @SQ"))?, 0);
        pos += 4 + l_name as usize + 4;
    }
    let n_ref = n_ref as usize;

    let mut out = Vec::new();
    out.extend_from_slice(b"BAI\x01");
    out.extend_from_slice(&(n_ref as i32).to_le_bytes());
    let null_content = |out: &mut Vec<u8>| out.extend_from_slice(&0u64.to_le_bytes());

    let mut current = 0usize; // BAMIndexer.currentReference
    let mut builder = Builder::default();
    // BAMIndexMetaData: per reference, and the file-wide no-coordinate count.
    let (mut first_off, mut last_off, mut aligned, mut unaligned): (i64, u64, i32, i32) =
        (-1, 0, 0, 0);
    let mut no_coordinate: u64 = 0;

    // Write the content of `current` (BinaryBAMIndexWriter.writeReference).
    let write_reference = |out: &mut Vec<u8>,
                           builder: Builder,
                           meta: (i64, u64, i32, i32)| {
        let Some((bins, linear)) = builder.finish() else {
            null_content(out);
            return;
        };
        let size = bins.keys().filter(|&&b| b != MAX_BINS).count();
        if size == 0 {
            null_content(out);
            return;
        }
        out.extend_from_slice(&(size as i32 + 1).to_le_bytes());
        for (&bin, chunks) in &bins {
            if bin == MAX_BINS {
                continue;
            }
            out.extend_from_slice(&bin.to_le_bytes());
            out.extend_from_slice(&(chunks.len() as i32).to_le_bytes());
            for &(s, e) in chunks {
                out.extend_from_slice(&s.to_le_bytes());
                out.extend_from_slice(&e.to_le_bytes());
            }
        }
        out.extend_from_slice(&MAX_BINS.to_le_bytes());
        out.extend_from_slice(&2i32.to_le_bytes());
        out.extend_from_slice(&meta.0.to_le_bytes());
        out.extend_from_slice(&meta.1.to_le_bytes());
        out.extend_from_slice(&(meta.2 as i64).to_le_bytes());
        out.extend_from_slice(&(meta.3 as i64).to_le_bytes());
        out.extend_from_slice(&(linear.len() as i32).to_le_bytes());
        for v in linear {
            out.extend_from_slice(&v.to_le_bytes());
        }
    };

    loop {
        let dropped = vr.compact(pos);
        pos -= dropped;
        let Some(len) = vr.bytes(pos, 4)? else { break };
        let block_size = le_i32(len, 0);
        if block_size < 32 {
            return Err(invalid("bad BAM record length"));
        }
        let rec = vr
            .bytes(pos + 4, block_size as usize)?
            .ok_or_else(|| invalid("truncated BAM record"))?
            .to_vec();
        let start_vo = vr.vo(pos)?;
        let end = pos + 4 + block_size as usize;
        let end_vo = vr.vo(end)?;
        pos = end;

        let ref_id = le_i32(&rec, 0);
        let pos0 = le_i32(&rec, 4);
        let l_read_name = rec[8] as usize;
        let n_cigar = le_u16(&rec, 12) as usize;
        let flag = le_u16(&rec, 14);
        let unmapped = flag & 0x4 != 0;
        // BAMIndexer.processAlignment: a new reference flushes the earlier ones.
        if ref_id >= 0 && ref_id as usize != current {
            let r = ref_id as usize;
            if r < current || r >= n_ref {
                return Err(invalid("BAM is not coordinate sorted by reference"));
            }
            while current < r {
                write_reference(
                    &mut out,
                    std::mem::take(&mut builder),
                    (first_off, last_off, aligned, unaligned),
                );
                current += 1;
                (first_off, last_off, aligned, unaligned) = (-1, 0, 0, 0);
            }
        }
        // getAlignmentStart: pos + 1, 0 when unset.
        let start = pos0 + 1;
        if start == 0 {
            no_coordinate += 1;
            continue;
        }
        if ref_id < 0 {
            return Err(invalid("record with a position but no reference"));
        }
        // BAMIndexMetaData.recordMetaData
        if unmapped {
            unaligned += 1;
        } else {
            aligned += 1;
        }
        if first_off == -1 || start_vo <= first_off as u64 {
            first_off = start_vo as i64;
        }
        if last_off <= end_vo {
            last_off = end_vo;
        }
        // getAlignmentEnd: 0 for unmapped reads, else start + reference length - 1.
        let alignment_end = if unmapped {
            0
        } else {
            let cig = 32 + l_read_name;
            let mut ref_len = 0i32;
            for k in 0..n_cigar {
                let op = u32::from_le_bytes(rec[cig + 4 * k..cig + 4 * k + 4].try_into().unwrap());
                if matches!(op & 0xf, 0 | 2 | 3 | 7 | 8) {
                    ref_len += (op >> 4) as i32;
                }
            }
            start + ref_len - 1
        };
        // SAMRecord.computeIndexingBin
        let bin_end = if alignment_end <= 0 {
            start - 1 + 1
        } else {
            alignment_end
        };
        let bin = region_to_bin(start - 1, bin_end);
        builder.feature(start, alignment_end, bin, (start_vo, end_vo));
    }
    // finish(): the remaining references, then the no-coordinate count.
    while current < n_ref {
        write_reference(
            &mut out,
            std::mem::take(&mut builder),
            (first_off, last_off, aligned, unaligned),
        );
        current += 1;
        (first_off, last_off, aligned, unaligned) = (-1, 0, 0, 0);
    }
    out.extend_from_slice(&no_coordinate.to_le_bytes());
    Ok(out)
}

/// htsjdk's `.tbi` for a BGZF VCF (`TabixFormat.VCF`), as GATK's on-the-fly
/// `TabixIndexCreator` builds it, BGZF-compressed like `TabixIndex.write`.
pub fn tbi_htsjdk_vcf(vcf_gz: &Path) -> io::Result<Vec<u8>> {
    let payload = tbi_payload_vcf(vcf_gz)?;
    let mut w = BgzfWriter::new(Vec::new(), HtsjdkDeflater::new(), Profile::HTSJDK);
    w.write_all(&payload)?;
    w.finish()
}

/// The uncompressed `.tbi` content.
pub fn tbi_payload_vcf(vcf_gz: &Path) -> io::Result<Vec<u8>> {
    let mut vr = VirtualReader::new(BufReader::with_capacity(1 << 20, File::open(vcf_gz)?));
    let mut names: Vec<String> = Vec::new();
    let mut contents: Vec<Option<(Bins, Vec<u64>)>> = Vec::new();
    let mut builder: Option<Builder> = None;
    // The previous feature waits for the next one's position (its end).
    let mut previous: Option<(i32, i32, u64)> = None;
    let mut pos = 0usize;
    loop {
        let dropped = vr.compact(pos);
        pos -= dropped;
        let line_start = pos;
        // Find the end of the line.
        let mut end = pos;
        loop {
            vr.fill(end + 1)?;
            if end >= vr.buf.len() {
                break;
            }
            let nl = vr.buf[end..].iter().position(|&b| b == b'\n');
            match nl {
                Some(k) => {
                    end += k + 1;
                    break;
                }
                None => {
                    end = vr.buf.len();
                    vr.fill(end + 1)?;
                    if end >= vr.buf.len() {
                        break;
                    }
                }
            }
        }
        if end == line_start {
            break;
        }
        let line = vr.buf[line_start..end].to_vec();
        pos = end;
        if line.first() == Some(&b'#') {
            continue;
        }
        let vo = vr.vo(line_start)?;
        let text = std::str::from_utf8(&line).map_err(|_| invalid("VCF line is not UTF-8"))?;
        let mut f = text.trim_end_matches(['\n', '\r']).split('\t');
        let chrom = f.next().unwrap_or("");
        let p: i32 = f
            .next()
            .and_then(|s| s.parse().ok())
            .ok_or_else(|| invalid("VCF line without POS"))?;
        let _id = f.next();
        let reference = f.next().ok_or_else(|| invalid("VCF line without REF"))?;
        let info = f.nth(3).unwrap_or(".");
        // VariantContext.getEnd: INFO END when present, else POS + |REF| - 1.
        let end_pos = info
            .split(';')
            .find_map(|kv| kv.strip_prefix("END="))
            .and_then(|v| v.parse::<i32>().ok())
            .unwrap_or(p + reference.len() as i32 - 1);
        if let Some((ps, pe, pvo)) = previous.take() {
            if pvo >= vo {
                return Err(invalid("feature start position >= feature end position"));
            }
            let s = ps - 1;
            let e = if pe <= 0 { s + 1 } else { pe };
            builder
                .as_mut()
                .expect("builder exists with a previous feature")
                .feature(ps, pe, region_to_bin(s, e), (pvo, vo));
        }
        if names.last().map(String::as_str) != Some(chrom) {
            if names.iter().any(|n| n == chrom) {
                return Err(invalid(format!("contig {chrom} out of order")));
            }
            if let Some(b) = builder.take() {
                contents.push(b.finish());
            }
            builder = Some(Builder::default());
            names.push(chrom.to_owned());
        }
        previous = Some((p, end_pos, vo));
    }
    if let Some((ps, pe, pvo)) = previous.take() {
        // IndexingVariantContextWriter.close flushes the BGZF stream before
        // `finalizeIndex(getPosition())`: the last feature ends at the next
        // block address even when the last block was not full.
        let final_vo = vr.vo_at(pos, true)?;
        let s = ps - 1;
        let e = if pe <= 0 { s + 1 } else { pe };
        builder
            .as_mut()
            .expect("builder exists with a previous feature")
            .feature(ps, pe, region_to_bin(s, e), (pvo, final_vo));
    }
    if let Some(b) = builder.take() {
        contents.push(b.finish());
    }

    let mut out = Vec::new();
    out.extend_from_slice(b"TBI\x01");
    out.extend_from_slice(&(names.len() as i32).to_le_bytes());
    // TabixFormat.VCF: flags 2 (VCF), sequence column 1, start 2, end 0, '#', skip 0.
    for v in [2i32, 1, 2, 0, b'#' as i32, 0] {
        out.extend_from_slice(&v.to_le_bytes());
    }
    let name_block: usize = names.iter().map(|n| n.len() + 1).sum();
    out.extend_from_slice(&(name_block as i32).to_le_bytes());
    for n in &names {
        out.extend_from_slice(n.as_bytes());
        out.push(0);
    }
    for c in contents {
        match c {
            None => out.extend_from_slice(&0i32.to_le_bytes()),
            Some((bins, linear)) => {
                out.extend_from_slice(&(bins.len() as i32).to_le_bytes());
                for (bin, chunks) in &bins {
                    out.extend_from_slice(&bin.to_le_bytes());
                    out.extend_from_slice(&(chunks.len() as i32).to_le_bytes());
                    for &(s, e) in chunks {
                        out.extend_from_slice(&s.to_le_bytes());
                        out.extend_from_slice(&e.to_le_bytes());
                    }
                }
                out.extend_from_slice(&(linear.len() as i32).to_le_bytes());
                for v in linear {
                    out.extend_from_slice(&v.to_le_bytes());
                }
            }
        }
    }
    Ok(out)
}

/// Write `bai_htsjdk(bam)` to `out`.
pub fn write_bai_htsjdk(bam: &Path, out: &Path) -> io::Result<()> {
    std::fs::write(out, bai_htsjdk(bam)?)
}

/// Write `tbi_htsjdk_vcf(vcf_gz)` to `out`.
pub fn write_tbi_htsjdk_vcf(vcf_gz: &Path, out: &Path) -> io::Result<()> {
    std::fs::write(out, tbi_htsjdk_vcf(vcf_gz)?)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bins_match_sam_spec() {
        assert_eq!(region_to_bin(0, 1), 4681);
        assert_eq!(region_to_bin(0, 1 << 14), 4681);
        assert_eq!(region_to_bin(0, (1 << 14) + 1), 585);
        assert_eq!(region_to_bin(0, 1 << 29), 0);
    }

    #[test]
    fn linear_index_fills_gaps_and_keeps_minimum() {
        let mut b = Builder::default();
        b.feature(1, 10, 4681, (100, 200));
        b.feature(40000, 40010, 4683, (300, 400));
        b.feature(40005, 40020, 4683, (250, 500));
        let (bins, linear) = b.finish().unwrap();
        assert_eq!(linear, vec![100, 100, 250]);
        assert_eq!(bins[&4683], vec![(300, 500)]); // same block address 0: coalesced
    }

    #[test]
    fn virtual_offsets_follow_htsjdk_flush_rule() {
        // Two full blocks then a partial one, as htsjdk's writer lays them out.
        let data: Vec<u8> = (0..FULL_BLOCK * 2 + 10).map(|i| i as u8).collect();
        let mut w = BgzfWriter::new(Vec::new(), HtsjdkDeflater::new(), Profile::HTSJDK);
        w.write_all(&data).unwrap();
        let file = w.finish().unwrap();
        let blocks = crate::bgzf::split_blocks(&file).unwrap();
        let mut vr = VirtualReader::new(&file[..]);
        assert_eq!(vr.vo(0).unwrap(), 0);
        assert_eq!(vr.vo(FULL_BLOCK).unwrap(), (blocks[1].start as u64) << 16);
        let end = vr.vo(2 * FULL_BLOCK + 10).unwrap();
        assert_eq!(end, (blocks[2].start as u64) << 16 | 10);
    }
}
