//! Read-pair filtering and chunking, equivalent to sarek's FASTP step:
//! fastp 1.1.0 with `--thread 12 --detect_adapter_for_pe
//! --disable_adapter_trimming --split_by_lines 200000000 --length_required 15`.
//!
//! What fastp does with those arguments (verified against fastp v1.1.0 source
//! and the tier0 golden chunks):
//!
//! - Pairs are read in packs of `PACK_SIZE = 256`; pack `p` goes to worker
//!   thread `p % threads`, and with split output each worker writes its own
//!   files, so chunk `t + 1` holds packs `t, t + threads, ...` in input order.
//! - `Filter::passFilter` per read, first failure wins: quality (more than
//!   `40 * len / 100.0` bases with quality char < `'0'`, i.e. Q < 15), then
//!   N bases (more than 5 `N`), then length (< 15). A pair is written only if
//!   both reads pass. Records are written unchanged unless poly-G trimmed.
//! - Adapter trimming is disabled, so overlap analysis, adapter-dimer marking
//!   and adapter detection change nothing.
//! - Poly-G trimming switches on when the first R1 read name starts with a
//!   two-colour instrument prefix (`Evaluator::isTwoColorSystem`, main.cpp
//!   L504-507; `--poly_g_min_len` 10). It runs on both reads before the
//!   filter (`peprocessor.cpp` L428-431, after a no-op `trimAndCut`) and
//!   shortens sequence and qualities in place ([`trim_poly_g`]).
//!
//! Gate: the chunk FASTQ contents (`*.fastp.fastq.gz`, decompressed).

use std::io::{self, BufRead, BufReader, Write};
use std::path::Path;

use nomics_io::bgzf::ParallelBgzfWriter;
use nomics_io::gzip::GzReader;

pub const PACK_SIZE: u64 = 256;
const QUALIFIED_QUAL: u8 = 15 + 33;
const UNQUALIFIED_PERCENT_LIMIT: f64 = 40.0;
const N_BASE_LIMIT: usize = 5;
const LENGTH_REQUIRED: usize = 15;
/// `--split_by_lines 200000000`: lines per output file before fastp rolls to a new one.
const SPLIT_LINES: u64 = 200_000_000;

/// fastp `Evaluator::isTwoColorSystem` read-name prefixes.
/// `--poly_g_min_len` default.
const POLY_G_MIN_LEN: usize = 10;

/// `PolyX::trimPolyG` for one read (fastp 1.1.0 `polyx.cpp` L16-42): scan
/// from the 3' end counting non-G bases; stop when they exceed 5, or exceed
/// one per 8 scanned bases once at least `compare_req` bases are scanned;
/// if the scan covered `compare_req` bases, cut at the leftmost G seen.
/// Returns the new length.
pub fn trim_poly_g(seq: &[u8], compare_req: usize) -> usize {
    const ALLOW_ONE_MISMATCH_FOR_EACH: usize = 8;
    const MAX_MISMATCH: usize = 5;
    let rlen = seq.len();
    let mut mismatch = 0;
    let mut first_g = rlen as i64 - 1;
    let mut i = 0;
    while i < rlen {
        if seq[rlen - i - 1] != b'G' {
            mismatch += 1;
        } else {
            first_g = (rlen - i - 1) as i64;
        }
        let allowed = (i + 1) / ALLOW_ONE_MISMATCH_FOR_EACH;
        if mismatch > MAX_MISMATCH || (mismatch > allowed && i + 1 >= compare_req) {
            break;
        }
        i += 1;
    }
    // Read::resize ignores lengths outside [0, length].
    if i >= compare_req && first_g >= 0 && (first_g as usize) <= rlen {
        first_g as usize
    } else {
        rlen
    }
}

const TWO_COLOR_PREFIXES: [&[u8]; 11] = [
    b"@FS", b"@MN", b"@SH", b"@NS", b"@NB", b"@NDX", b"@VL", b"@VH", b"@A", b"@NA", b"@LH",
];

fn invalid(msg: String) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg)
}

/// `Filter::passFilter` results; discriminants are fastp's codes, whose
/// numeric `max` decides how a failed pair is counted.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Verdict {
    Pass = 0,
    TooManyN = 12,
    TooShort = 16,
    LowQuality = 20,
}

/// fastp `Filter::passFilter` with this pipeline's options.
pub fn filter_read(seq: &[u8], qual: &[u8]) -> Verdict {
    let len = seq.len();
    if len == 0 {
        return Verdict::TooShort;
    }
    let low_qual = qual.iter().filter(|&&q| q < QUALIFIED_QUAL).count();
    let n_bases = seq.iter().filter(|&&b| b == b'N').count();
    if low_qual as f64 > UNQUALIFIED_PERCENT_LIMIT * len as f64 / 100.0 {
        Verdict::LowQuality
    } else if n_bases > N_BASE_LIMIT {
        Verdict::TooManyN
    } else if len < LENGTH_REQUIRED {
        Verdict::TooShort
    } else {
        Verdict::Pass
    }
}

#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct Stats {
    pub pairs_in: u64,
    pub pairs_out: u64,
    /// fastp's `filtering_result`: each pair adds 2 under the worse
    /// (numerically larger) result of its two reads.
    pub passed_reads: u64,
    pub low_quality_reads: u64,
    pub too_many_n_reads: u64,
    pub too_short_reads: u64,
}

/// Surviving pairs, split into fastp's chunks: `r1[c]`/`r2[c]` hold chunk
/// `c + 1` as FASTQ text, records unchanged.
pub struct Chunks {
    pub r1: Vec<Vec<u8>>,
    pub r2: Vec<Vec<u8>>,
}

struct Fastq<R> {
    r: BufReader<R>,
    lines: [Vec<u8>; 4],
}

impl<R: io::Read> Fastq<R> {
    fn next(&mut self) -> io::Result<bool> {
        for (i, line) in self.lines.iter_mut().enumerate() {
            line.clear();
            if self.r.read_until(b'\n', line)? == 0 {
                return if i == 0 {
                    Ok(false)
                } else {
                    Err(invalid("truncated FASTQ record".into()))
                };
            }
        }
        if self.lines[0].first() != Some(&b'@') || self.lines[2].first() != Some(&b'+') {
            return Err(invalid(format!(
                "malformed FASTQ record {:?}",
                String::from_utf8_lossy(&self.lines[0])
            )));
        }
        Ok(true)
    }

    fn seq(&self) -> &[u8] {
        trim_eol(&self.lines[1])
    }

    fn qual(&self) -> &[u8] {
        trim_eol(&self.lines[3])
    }

    /// Poly-G trim this record in place (sequence and qualities).
    fn trim_poly_g(&mut self, compare_req: usize) {
        let len = self.seq().len();
        let new = trim_poly_g(self.seq(), compare_req);
        if new < len {
            for k in [1, 3] {
                let eol: Vec<u8> = self.lines[k][len..].to_vec();
                self.lines[k].truncate(new);
                self.lines[k].extend_from_slice(&eol);
            }
        }
    }

    fn append_to(&self, out: &mut Vec<u8>) {
        for l in &self.lines {
            out.extend_from_slice(l);
        }
    }

    fn write_to(&self, out: &mut impl Write) -> io::Result<()> {
        for l in &self.lines {
            out.write_all(l)?;
        }
        Ok(())
    }
}

fn trim_eol(l: &[u8]) -> &[u8] {
    let l = l.strip_suffix(b"\n").unwrap_or(l);
    l.strip_suffix(b"\r").unwrap_or(l)
}

/// Filter and chunk a FASTQ pair exactly as fastp does with `fastp_threads` workers.
pub fn run(fq1: &Path, fq2: &Path, fastp_threads: usize) -> io::Result<(Chunks, Stats)> {
    let threads = fastp_threads.max(1);
    let mut chunks = Chunks {
        r1: vec![Vec::new(); threads],
        r2: vec![Vec::new(); threads],
    };
    let stats = run_into(fq1, fq2, threads, |chunk, a, b| {
        a.append_to(&mut chunks.r1[chunk]);
        b.append_to(&mut chunks.r2[chunk]);
        Ok(())
    })?;
    Ok((chunks, stats))
}

/// Paths of the uncompressed chunk FASTQs [`run_to_dir`] writes.
pub fn spill_paths(dir: &Path, chunk: usize) -> [std::path::PathBuf; 2] {
    [1, 2].map(|mate| dir.join(format!("{:04}_R{mate}.fastq", chunk + 1)))
}

/// [`run`], but each chunk is written as uncompressed FASTQ text to
/// [`spill_paths`] instead of being kept in memory (whole genomes).
pub fn run_to_dir(fq1: &Path, fq2: &Path, fastp_threads: usize, dir: &Path) -> io::Result<Stats> {
    std::fs::create_dir_all(dir)?;
    let threads = fastp_threads.max(1);
    let mut writers = Vec::with_capacity(threads);
    for chunk in 0..threads {
        let [p1, p2] = spill_paths(dir, chunk);
        writers.push([
            io::BufWriter::with_capacity(1 << 20, std::fs::File::create(p1)?),
            io::BufWriter::with_capacity(1 << 20, std::fs::File::create(p2)?),
        ]);
    }
    let stats = run_into(fq1, fq2, threads, |chunk, a, b| {
        let [w1, w2] = &mut writers[chunk];
        a.write_to(w1)?;
        b.write_to(w2)
    })?;
    for [mut w1, mut w2] in writers {
        w1.flush()?;
        w2.flush()?;
    }
    Ok(stats)
}

/// The fastp filter and pack-to-chunk assignment; `emit(chunk, r1, r2)`
/// receives every surviving pair in input order.
fn run_into(
    fq1: &Path,
    fq2: &Path,
    fastp_threads: usize,
    mut emit: impl FnMut(usize, &Fastq<GzReader>, &Fastq<GzReader>) -> io::Result<()>,
) -> io::Result<Stats> {
    let threads = fastp_threads.max(1) as u64;
    let mut a = Fastq {
        r: BufReader::with_capacity(1 << 20, GzReader::open(fq1)?),
        lines: Default::default(),
    };
    let mut b = Fastq {
        r: BufReader::with_capacity(1 << 20, GzReader::open(fq2)?),
        lines: Default::default(),
    };
    let mut lines_out = vec![0u64; threads as usize];
    let mut poly_g = false;
    let mut stats = Stats::default();
    loop {
        match (a.next()?, b.next()?) {
            (false, false) => break,
            (true, true) => {}
            _ => {
                return Err(invalid(format!(
                    "FASTQ inputs differ in length after {} pairs",
                    stats.pairs_in
                )));
            }
        }
        if stats.pairs_in == 0 {
            poly_g = TWO_COLOR_PREFIXES.iter().any(|p| a.lines[0].starts_with(p));
        }
        if poly_g {
            a.trim_poly_g(POLY_G_MIN_LEN);
            b.trim_poly_g(POLY_G_MIN_LEN);
        }
        let chunk = ((stats.pairs_in / PACK_SIZE) % threads) as usize;
        stats.pairs_in += 1;
        let (v1, v2) = (
            filter_read(a.seq(), a.qual()),
            filter_read(b.seq(), b.qual()),
        );
        match v1.max(v2) {
            Verdict::Pass => stats.passed_reads += 2,
            Verdict::LowQuality => stats.low_quality_reads += 2,
            Verdict::TooManyN => stats.too_many_n_reads += 2,
            Verdict::TooShort => stats.too_short_reads += 2,
        }
        if v1 == Verdict::Pass && v2 == Verdict::Pass {
            stats.pairs_out += 1;
            lines_out[chunk] += 4;
            if lines_out[chunk] > SPLIT_LINES {
                return Err(invalid(
                    "a chunk exceeds --split_by_lines; fastp's file roll-over is not implemented"
                        .into(),
                ));
            }
            emit(chunk, &a, &b)?;
        }
    }
    Ok(stats)
}

/// Write chunks as `{dir}/{NNNN}.{prefix}_R{1,2}.fastp.fastq.gz` (BGZF, valid gzip).
pub fn write_chunks(chunks: &Chunks, dir: &Path, prefix: &str, threads: usize) -> io::Result<()> {
    std::fs::create_dir_all(dir)?;
    for (c, (r1, r2)) in chunks.r1.iter().zip(&chunks.r2).enumerate() {
        write_chunk(r1, r2, dir, prefix, c + 1, threads)?;
    }
    Ok(())
}

/// Write one chunk with the same naming and compression as [`write_chunks`].
pub fn write_chunk(
    r1: &[u8],
    r2: &[u8],
    dir: &Path,
    prefix: &str,
    chunk: usize,
    threads: usize,
) -> io::Result<()> {
    std::fs::create_dir_all(dir)?;
    for (mate, text) in [(1, r1), (2, r2)] {
        let path = dir.join(format!("{chunk:04}.{prefix}_R{mate}.fastp.fastq.gz"));
        let mut w =
            ParallelBgzfWriter::new(io::BufWriter::new(std::fs::File::create(path)?), threads);
        w.write_all(text)?;
        w.finish()?.flush()?;
    }
    Ok(())
}

/// Compress chunk `chunk` (0-based) of a [`run_to_dir`] spill into the
/// files [`write_chunk`] would write, streaming from disk.
pub fn compress_spill(
    spill: &Path,
    chunk: usize,
    dir: &Path,
    prefix: &str,
    threads: usize,
) -> io::Result<()> {
    std::fs::create_dir_all(dir)?;
    for (mate, src) in [1, 2].into_iter().zip(spill_paths(spill, chunk)) {
        let path = dir.join(format!("{:04}.{prefix}_R{mate}.fastp.fastq.gz", chunk + 1));
        let mut w =
            ParallelBgzfWriter::new(io::BufWriter::new(std::fs::File::create(path)?), threads);
        io::copy(&mut std::fs::File::open(src)?, &mut w)?;
        w.finish()?.flush()?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn poly_g_follows_fastp() {
        // Ten trailing G's: cut at the first of them.
        let mut s = b"ACGTACGTAC".to_vec();
        s.extend(std::iter::repeat_n(b'G', 10));
        assert_eq!(trim_poly_g(&s, 10), 10);
        // Nine G's: the first non-G (i = 9) is within one mismatch per 8
        // bases, the scan breaks at i = 10 >= 10, so they are cut too.
        let s9 = [b"ACGTACGTAC".as_slice(), &[b'G'; 9]].concat();
        assert_eq!(trim_poly_g(&s9, 10), 10);
        // Eight G's after two non-G: the scan breaks at i = 9 < 10.
        let s8 = [b"ACGTACGTTT".as_slice(), &[b'G'; 8]].concat();
        assert_eq!(trim_poly_g(&s8, 10), s8.len());
        // One mismatch per 8 bases is tolerated once 10 are scanned.
        let s1 = [b"TTTTTTTTTT".as_slice(), b"GGGGAGGGGGGGGGGG"].concat();
        assert_eq!(trim_poly_g(&s1, 10), 10);
        // All G: cut everything.
        assert_eq!(trim_poly_g(&[b'G'; 20], 10), 0);
        // No G tail: unchanged.
        assert_eq!(trim_poly_g(b"ACGTACGTACGTACGT", 10), 16);
    }

    #[test]
    fn spill_matches_in_memory_chunks() {
        let dir = std::env::temp_dir().join(format!("nomics-prep-spill-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let mut texts = [Vec::new(), Vec::new()];
        for i in 0..2000u32 {
            for (mate, text) in texts.iter_mut().enumerate() {
                let qual = if i % 7 == mate as u32 { b'#' } else { b'F' };
                let len = 20 + (i % 50) as usize;
                text.extend_from_slice(format!("@r{i}/{}\n", mate + 1).as_bytes());
                text.extend(std::iter::repeat_n(b"ACGT"[i as usize % 4], len));
                text.extend_from_slice(b"\n+\n");
                text.extend(std::iter::repeat_n(qual, len));
                text.push(b'\n');
            }
        }
        let fq: Vec<_> = (0..2)
            .map(|m| {
                let path = dir.join(format!("in_{m}.fastq.gz"));
                let mut w = ParallelBgzfWriter::new(std::fs::File::create(&path).unwrap(), 2);
                w.write_all(&texts[m]).unwrap();
                w.finish().unwrap();
                path
            })
            .collect();
        let (chunks, stats) = run(&fq[0], &fq[1], 3).unwrap();
        let spill = dir.join("spill");
        assert_eq!(run_to_dir(&fq[0], &fq[1], 3, &spill).unwrap(), stats);
        for c in 0..3 {
            let [p1, p2] = spill_paths(&spill, c);
            assert_eq!(std::fs::read(p1).unwrap(), chunks.r1[c]);
            assert_eq!(std::fs::read(p2).unwrap(), chunks.r2[c]);
        }
        assert!(stats.pairs_out > 0 && stats.pairs_out < stats.pairs_in);
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn filter_boundaries_follow_fastp() {
        let seq = vec![b'A'; 100];
        // 40 low-quality bases out of 100: 40 > 40.0 is false -> pass; 41 fails.
        let mut q = vec![b'F'; 100];
        q[..40].fill(b'/'); // '/' = Q14
        assert_eq!(filter_read(&seq, &q), Verdict::Pass);
        q[40] = b'/';
        assert_eq!(filter_read(&seq, &q), Verdict::LowQuality);
        // Q15 ('0') is qualified.
        assert_eq!(filter_read(&seq, &[b'0'; 100]), Verdict::Pass);
        // Five Ns pass, six fail; quality is checked first.
        let mut s = seq.clone();
        s[..5].fill(b'N');
        assert_eq!(filter_read(&s, &[b'F'; 100]), Verdict::Pass);
        s[5] = b'N';
        assert_eq!(filter_read(&s, &[b'F'; 100]), Verdict::TooManyN);
        assert_eq!(filter_read(&s, &[b'#'; 100]), Verdict::LowQuality);
        // Length 14 fails, 15 passes.
        assert_eq!(filter_read(&[b'A'; 14], &[b'F'; 14]), Verdict::TooShort);
        assert_eq!(filter_read(&[b'A'; 15], &[b'F'; 15]), Verdict::Pass);
    }
}
