//! Build the tier0 development fixture: a subset of the tier1 FASTQ pairs,
//! chosen from tier1 goldens, small enough for a minutes-long reference run.
//!
//! A read pair (keyed by the run-local spot number `N` in `SRR14724482.N`) is
//! kept when:
//! - **target**: any mapped record of the pair overlaps a target interval on
//!   `--chrom` in the golden BAM, or
//! - **fastp-dropped** (sampled): the pair is absent from the golden BAM, i.e.
//!   fastp filtered it, and `N` falls in the 1-in-`--sample-every` sample, or
//! - **unmapped** (sampled): both mates are unmapped, and `N` is sampled.
//!
//! Pairs are written in their original order, byte-for-byte, as BGZF-compressed
//! FASTQ. tier0 goldens must come from running the reference on these FASTQs:
//! BQSR and duplicate marking depend on the whole input.

use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use nomics_io::bam::{Header, RecordView, read_raw_record};
use nomics_io::bgzf::{BgzfReader, BgzfWriter, Profile};
use nomics_io::deflate::HtsjdkDeflater;
use nomics_io::gzip::GzReader;

pub struct Args {
    pub golden_bam: PathBuf,
    pub targets: PathBuf,
    pub chrom: String,
    pub fastq1: PathBuf,
    pub fastq2: PathBuf,
    pub outdir: PathBuf,
    pub sample_every: u64,
}

#[derive(Default)]
struct Bits(Vec<u64>);

impl Bits {
    fn set(&mut self, n: u64) {
        let w = (n / 64) as usize;
        if w >= self.0.len() {
            self.0.resize(w + 1, 0);
        }
        self.0[w] |= 1 << (n % 64);
    }

    fn get(&self, n: u64) -> bool {
        self.0
            .get((n / 64) as usize)
            .is_some_and(|w| w >> (n % 64) & 1 == 1)
    }
}

fn invalid(msg: String) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg)
}

/// Spot number from `SRR14724482.N` (BAM) or `@SRR14724482.N N/1` (FASTQ).
fn spot(name: &[u8]) -> io::Result<u64> {
    let name = name.strip_prefix(b"@").unwrap_or(name);
    let end = name
        .iter()
        .position(|&c| c == b' ' || c == b'\t' || c == b'\n')
        .unwrap_or(name.len());
    let dot = name[..end].iter().rposition(|&c| c == b'.');
    dot.and_then(|d| std::str::from_utf8(&name[d + 1..end]).ok()?.parse().ok())
        .ok_or_else(|| {
            invalid(format!(
                "no spot number in {:?}",
                String::from_utf8_lossy(name)
            ))
        })
}

/// splitmix64, for a sample that does not correlate with spot order.
fn sampled(n: u64, every: u64) -> bool {
    let mut z = n.wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    (z ^ (z >> 31)).is_multiple_of(every)
}

/// Merged, sorted half-open intervals on `chrom` from a BED file.
fn read_targets(path: &Path, chrom: &str) -> io::Result<Vec<(i64, i64)>> {
    let mut iv = Vec::new();
    for line in BufReader::new(File::open(path)?).lines() {
        let line = line?;
        let mut f = line.split('\t');
        if f.next() != Some(chrom) {
            continue;
        }
        let parse = |s: Option<&str>| s.and_then(|s| s.trim().parse::<i64>().ok());
        match (parse(f.next()), parse(f.next())) {
            (Some(a), Some(b)) if a < b => iv.push((a, b)),
            _ => return Err(invalid(format!("bad BED line: {line}"))),
        }
    }
    iv.sort_unstable();
    let mut merged: Vec<(i64, i64)> = Vec::with_capacity(iv.len());
    for (a, b) in iv {
        match merged.last_mut() {
            Some(last) if a <= last.1 => last.1 = last.1.max(b),
            _ => merged.push((a, b)),
        }
    }
    Ok(merged)
}

fn overlaps(iv: &[(i64, i64)], beg: i64, end: i64) -> bool {
    let i = iv.partition_point(|&(_, e)| e <= beg);
    i < iv.len() && iv[i].0 < end
}

struct Fastq<R> {
    r: BufReader<R>,
    lines: [Vec<u8>; 4],
}

impl<R: io::Read> Fastq<R> {
    /// Next record as its four raw lines (with newlines). `Ok(false)` at EOF.
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

    fn write_to(&self, w: &mut impl Write) -> io::Result<()> {
        self.lines.iter().try_for_each(|l| w.write_all(l))
    }
}

pub fn run(a: &Args) -> io::Result<()> {
    let targets = read_targets(&a.targets, &a.chrom)?;
    eprintln!("{} merged {} target intervals", targets.len(), a.chrom);

    let mut bam = BgzfReader::new(BufReader::with_capacity(
        1 << 20,
        File::open(&a.golden_bam)?,
    ));
    let header = Header::read(&mut bam)?;
    let chrom_id = header
        .refs
        .iter()
        .position(|r| r.name == a.chrom.as_bytes())
        .ok_or_else(|| invalid(format!("{} not in BAM header", a.chrom)))?
        as i32;
    let (mut in_bam, mut on_target, mut unmapped) =
        (Bits::default(), Bits::default(), Bits::default());
    let (mut raw, mut records) = (Vec::new(), 0u64);
    while read_raw_record(&mut bam, &mut raw)? {
        records += 1;
        let v = RecordView::new(&raw)?;
        let n = spot(v.read_name())?;
        in_bam.set(n);
        let flag = v.flag();
        if flag & 4 == 0 && v.ref_id() == chrom_id {
            let beg = v.pos() as i64;
            if overlaps(&targets, beg, beg + (v.reference_length() as i64).max(1)) {
                on_target.set(n);
            }
        }
        if flag & 0xc == 0xc {
            unmapped.set(n);
        }
    }
    eprintln!("scanned {records} golden records");

    std::fs::create_dir_all(&a.outdir)?;
    let open_out = |r: u8| -> io::Result<_> {
        let f = File::create(a.outdir.join(format!("tier0_{r}.fastq.gz")))?;
        Ok(BgzfWriter::new(
            BufWriter::new(f),
            HtsjdkDeflater::new(),
            Profile::HTSJDK,
        ))
    };
    let (mut w1, mut w2) = (open_out(1)?, open_out(2)?);
    let mut f1 = Fastq {
        r: BufReader::with_capacity(1 << 20, GzReader::open(&a.fastq1)?),
        lines: Default::default(),
    };
    let mut f2 = Fastq {
        r: BufReader::with_capacity(1 << 20, GzReader::open(&a.fastq2)?),
        lines: Default::default(),
    };
    let (
        mut pairs,
        mut kept_target,
        mut kept_dropped,
        mut kept_unmapped,
        mut dropped_total,
        mut unmapped_total,
    ) = (0u64, 0u64, 0u64, 0u64, 0u64, 0u64);
    loop {
        match (f1.next()?, f2.next()?) {
            (false, false) => break,
            (true, true) => {}
            _ => {
                return Err(invalid(format!(
                    "FASTQ files differ in length after {pairs} pairs"
                )));
            }
        }
        pairs += 1;
        let n = spot(&f1.lines[0])?;
        if spot(&f2.lines[0])? != n {
            return Err(invalid(format!("mate names disagree at pair {pairs}")));
        }
        let dropped = !in_bam.get(n);
        let both_unmapped = unmapped.get(n);
        dropped_total += dropped as u64;
        unmapped_total += both_unmapped as u64;
        let keep = if on_target.get(n) {
            kept_target += 1;
            true
        } else if dropped && sampled(n, a.sample_every) {
            kept_dropped += 1;
            true
        } else if both_unmapped && sampled(n, a.sample_every) {
            kept_unmapped += 1;
            true
        } else {
            false
        };
        if keep {
            f1.write_to(&mut w1)?;
            f2.write_to(&mut w2)?;
        }
    }
    w1.finish()?.flush()?;
    w2.finish()?.flush()?;

    let kept = kept_target + kept_dropped + kept_unmapped;
    let manifest = format!(
        "key\tvalue\n\
         source_fastq1\t{}\nsource_fastq2\t{}\ngolden_bam\t{}\ntargets\t{}\nchrom\t{}\nsample_every\t{}\n\
         input_pairs\t{pairs}\ngolden_records\t{records}\nfastp_dropped_pairs\t{dropped_total}\n\
         both_unmapped_pairs\t{unmapped_total}\nkept_target\t{kept_target}\nkept_fastp_dropped\t{kept_dropped}\n\
         kept_unmapped\t{kept_unmapped}\nkept_total\t{kept}\n",
        a.fastq1.display(),
        a.fastq2.display(),
        a.golden_bam.display(),
        a.targets.display(),
        a.chrom,
        a.sample_every
    );
    std::fs::write(a.outdir.join("MANIFEST.tsv"), &manifest)?;
    print!("{manifest}");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn spot_numbers() {
        assert_eq!(spot(b"SRR14724482.104665").unwrap(), 104665);
        assert_eq!(spot(b"@SRR14724482.1 1/1\n").unwrap(), 1);
        assert!(spot(b"@nodot 1/1").is_err());
    }

    #[test]
    fn interval_overlap_is_half_open() {
        let iv = [(100, 200), (300, 400)];
        assert!(!overlaps(&iv, 50, 100));
        assert!(overlaps(&iv, 50, 101));
        assert!(overlaps(&iv, 199, 250));
        assert!(!overlaps(&iv, 200, 300));
        assert!(overlaps(&iv, 399, 500));
        assert!(!overlaps(&iv, 400, 500));
    }

    #[test]
    fn bitset() {
        let mut b = Bits::default();
        b.set(0);
        b.set(41_000_123);
        assert!(b.get(0) && b.get(41_000_123) && !b.get(1) && !b.get(u64::MAX / 2));
    }
}
