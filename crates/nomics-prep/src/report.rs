//! fastp v1.1.0's sarek-mode report. See docs/spec/qc-fastp-report.md.
use std::fmt::Write as _;
use std::io;

use nomics_mem::HostBuf;

mod html;
pub use html::timestamp_now;

#[derive(Default)]
struct Cycle {
    count: [u64; 8],
    quality: [u64; 8],
}

struct ReadStats {
    reads: u64,
    bases: u64,
    q: [u64; 3],
    gc: u64,
    cycles: Vec<Cycle>,
    kmers: [u64; 1024],
    qual_hist: [u64; 128],
}

impl Default for ReadStats {
    fn default() -> Self {
        Self {
            reads: 0,
            bases: 0,
            q: [0; 3],
            gc: 0,
            cycles: Vec::new(),
            kmers: [0; 1024],
            qual_hist: [0; 128],
        }
    }
}

impl ReadStats {
    fn add(&mut self, seq: &[u8], qual: &[u8]) {
        self.reads += 1;
        self.bases += seq.len() as u64;
        self.cycles
            .resize_with(self.cycles.len().max(seq.len()), Cycle::default);
        let mut kmer = 0;
        let mut valid = 0;
        for (i, (&b, &q)) in seq.iter().zip(qual).enumerate() {
            self.qual_hist[q as usize] += 1;
            let bucket = (b & 7) as usize;
            self.cycles[i].count[bucket] += 1;
            self.cycles[i].quality[bucket] += u64::from(q - 33);
            self.gc += u64::from(bucket == (b'C' & 7) as usize || bucket == (b'G' & 7) as usize);
            self.q[0] += u64::from(q >= 53);
            self.q[1] += u64::from(q >= 63);
            self.q[2] += u64::from((73..127).contains(&q));
            let val = match b {
                b'A' => 0,
                b'T' => 1,
                b'C' => 2,
                b'G' => 3,
                _ => {
                    valid = 0;
                    continue;
                }
            };
            kmer = ((kmer << 2) & 1023) | val;
            valid += 1;
            if valid >= 5 {
                self.kmers[kmer] += 1;
            }
        }
    }

    fn mean_length(&self) -> u64 {
        self.bases.checked_div(self.reads).unwrap_or(0)
    }

    fn curve(&self, base: &str, quality: bool) -> Vec<f64> {
        self.cycles
            .iter()
            .map(|c| {
                let count: u64 = c.count.iter().sum();
                if quality {
                    let (q, n) = if base == "mean" {
                        (c.quality.iter().sum(), count)
                    } else {
                        let b = (base.as_bytes()[0] & 7) as usize;
                        if c.count[b] == 0 {
                            (c.quality.iter().sum(), count)
                        } else {
                            (c.quality[b], c.count[b])
                        }
                    };
                    q as f64 / n as f64
                } else {
                    let n = if base == "GC" {
                        c.count[(b'G' & 7) as usize] + c.count[(b'C' & 7) as usize]
                    } else {
                        c.count[(base.as_bytes()[0] & 7) as usize]
                    };
                    n as f64 / count as f64
                }
            })
            .collect()
    }

    fn json(&self, s: &mut String) {
        s.push_str("{\n");
        for (name, value) in [
            ("total_reads", self.reads),
            ("total_bases", self.bases),
            ("q20_bases", self.q[0]),
            ("q30_bases", self.q[1]),
            ("q40_bases", self.q[2]),
            ("total_cycles", self.cycles.len() as u64),
        ] {
            writeln!(s, "\t\t\"{name}\": {value},").unwrap();
        }
        for (name, names, quality) in [
            ("quality_curves", &["A", "T", "C", "G", "mean"][..], true),
            (
                "content_curves",
                &["A", "T", "C", "G", "N", "GC"][..],
                false,
            ),
        ] {
            writeln!(s, "\t\t\"{name}\": {{").unwrap();
            for (i, base) in names.iter().enumerate() {
                write!(s, "\t\t\t\"{base}\":[").unwrap();
                list(s, self.curve(base, quality).into_iter().map(float));
                writeln!(s, "]{}", if i + 1 == names.len() { "" } else { "," }).unwrap();
            }
            s.push_str("\t\t},\n");
        }
        s.push_str("\t\t\"kmer_count\": {\n");
        for i in 0..1024 {
            s.push_str("\t\t\t\"");
            for shift in (0..5).rev() {
                s.push(b"ATCG"[(i >> (shift * 2)) & 3] as char);
            }
            write!(s, "\":{}", self.kmers[i]).unwrap();
            if i != 1023 {
                s.push(',');
            }
            if i % 16 == 15 {
                s.push('\n');
            }
        }
        s.push_str("\t\t},\n\t\t\"overrepresented_sequences\": {\n\t\t}\n\t},\n");
    }
}

/// Serial FASTQ-order evaluation of fastp's two-plane Bloom filter.
/// The reference's concurrent two-plane updates are not a transaction.
struct Duplicate {
    bits: HostBuf,
    primes: [u64; 1024],
    total: u64,
    duplicates: u64,
}

impl Duplicate {
    fn new() -> io::Result<Self> {
        let mut primes = [0; 1024];
        let mut n = 10000u64;
        for p in &mut primes {
            loop {
                n += 1;
                let mut divisor = 2;
                while divisor * divisor <= n && !n.is_multiple_of(divisor) {
                    divisor += 1;
                }
                if divisor * divisor > n {
                    break;
                }
            }
            *p = n;
            n += 10000;
        }
        Ok(Self {
            bits: HostBuf::new(1 << 30)?,
            primes,
            total: 0,
            duplicates: 0,
        })
    }

    fn add(&mut self, a: &[u8], b: &[u8]) {
        let mut hash = [0u64; 2];
        for (p, &b) in a.iter().chain(b).enumerate() {
            let base = match b {
                b'A' => 7,
                b'T' => 222,
                b'C' => 74,
                b'G' => 31,
                _ => 13,
            };
            for (plane, h) in hash.iter_mut().enumerate() {
                *h = h.wrapping_add(
                    self.primes[(p * 2 + plane) & 1023].wrapping_mul(base + p as u64),
                );
            }
        }
        let bytes = self.bits.as_mut_slice();
        let mut duplicate = true;
        for (plane, h) in hash.into_iter().enumerate() {
            let pos = h & 0xffff_ffff;
            let byte = &mut bytes[plane * (1 << 29) + (pos >> 3) as usize];
            let bit = 1 << (pos & 7);
            duplicate &= *byte & bit != 0;
            *byte |= bit;
        }
        self.total += 1;
        self.duplicates += u64::from(duplicate);
    }

    fn rate(&self) -> f64 {
        ratio(self.duplicates, self.total)
    }
}

/// Statistics collected optionally in the existing filter pass, before and
/// after trimming/filtering. No FASTQ text or extra input pass is retained.
pub struct Report {
    reads: [ReadStats; 4],
    duplicate: Duplicate,
    insert: [u64; 513],
}

impl Report {
    pub fn new() -> io::Result<Self> {
        Ok(Self {
            reads: std::array::from_fn(|_| ReadStats::default()),
            duplicate: Duplicate::new()?,
            insert: [0; 513],
        })
    }

    pub(crate) fn before(&mut self, a: &[u8], aq: &[u8], b: &[u8], bq: &[u8]) {
        self.reads[0].add(a, aq);
        self.reads[1].add(b, bq);
        self.duplicate.add(a, b);
    }

    pub(crate) fn sample_insert(&mut self, a: &[u8], b: &[u8]) {
        self.insert[insert_size(a, b)] += 1;
    }

    pub(crate) fn after(&mut self, a: &[u8], aq: &[u8], b: &[u8], bq: &[u8]) {
        self.reads[2].add(a, aq);
        self.reads[3].add(b, bq);
    }

    fn peak(&self) -> usize {
        // First maximum, including bin zero, excluding the unknown bin.
        let mut peak = 0;
        for i in 1..512 {
            if self.insert[i] > self.insert[peak] {
                peak = i;
            }
        }
        peak
    }

    /// Literal fastp JSON, including the command's supplied trailing space.
    /// `command` is inserted verbatim, as fastp does (no escaping).
    pub fn json(&self, stats: &super::Stats, command: &str) -> String {
        let mut s = format!(
            "{{\n\t\"summary\": {{\n\t\t\"fastp_version\": \"1.1.0\",\n\t\t\"sequencing\": \"paired end ({} cycles + {} cycles)\",\n",
            self.reads[0].cycles.len(),
            self.reads[1].cycles.len()
        );
        for (part, name) in ["before_filtering", "after_filtering"]
            .into_iter()
            .enumerate()
        {
            let a = &self.reads[part * 2];
            let b = &self.reads[part * 2 + 1];
            writeln!(s, "\t\t\"{name}\": {{").unwrap();
            for (key, value) in [
                ("total_reads", a.reads + b.reads),
                ("total_bases", a.bases + b.bases),
                ("q20_bases", a.q[0] + b.q[0]),
                ("q30_bases", a.q[1] + b.q[1]),
            ] {
                writeln!(s, "\t\t\t\"{key}\":{value},").unwrap();
            }
            for (key, value) in [
                ("q20_rate", ratio(a.q[0] + b.q[0], a.bases + b.bases)),
                ("q30_rate", ratio(a.q[1] + b.q[1], a.bases + b.bases)),
            ] {
                writeln!(s, "\t\t\t\"{key}\":{},", float(value)).unwrap();
            }
            writeln!(s, "\t\t\t\"read1_mean_length\":{},\n\t\t\t\"read2_mean_length\":{},\n\t\t\t\"gc_content\":{}\n\t\t}}{}", a.mean_length(), b.mean_length(), float(ratio(a.gc+b.gc, a.bases+b.bases)), if part == 0 { "," } else { "" }).unwrap();
        }
        s.push_str("\t},\n\t\"filtering_result\": {\n");
        for (key, value) in [
            ("passed_filter_reads", stats.passed_reads),
            ("low_quality_reads", stats.low_quality_reads),
            ("too_many_N_reads", stats.too_many_n_reads),
            ("too_short_reads", stats.too_short_reads),
        ] {
            writeln!(s, "\t\t\"{key}\": {value},").unwrap();
        }
        write!(s, "\t\t\"too_long_reads\": 0\n\t}},\n\t\"duplication\": {{\n\t\t\"rate\": {}\n\t}},\n\t\"insert_size\": {{\n\t\t\"peak\": {},\n\t\t\"unknown\": {},\n\t\t\"histogram\": [", float(self.duplicate.rate()), self.peak(), self.insert[512]).unwrap();
        list(&mut s, self.insert[..512].iter().map(u64::to_string));
        s.push_str("]\n\t},\n");
        for (i, name) in [
            "read1_before_filtering",
            "read2_before_filtering",
            "read1_after_filtering",
            "read2_after_filtering",
        ]
        .into_iter()
        .enumerate()
        {
            write!(s, "\t\"{name}\": ").unwrap();
            self.reads[i].json(&mut s);
        }
        write!(s, "\t\"command\": \"{command}\"\n}}").unwrap();
        s
    }
}

fn ratio(n: u64, d: u64) -> f64 {
    if d == 0 { 0.0 } else { n as f64 / d as f64 }
}

fn list(s: &mut String, items: impl Iterator<Item = String>) {
    for (i, item) in items.enumerate() {
        if i > 0 {
            s.push(',');
        }
        s.push_str(&item);
    }
}

/// C++ defaultfloat, precision 6. Decide notation after decimal rounding.
fn float(v: f64) -> String {
    if !v.is_finite() {
        return v.to_string().to_lowercase();
    }
    if v == 0.0 {
        return if v.is_sign_negative() { "-0" } else { "0" }.into();
    }
    let scientific = format!("{v:.5e}");
    let (mantissa, exponent) = scientific.split_once('e').unwrap();
    let exponent: i32 = exponent.parse().unwrap();
    if !(-4..6).contains(&exponent) {
        return format!(
            "{}e{exponent:+03}",
            mantissa.trim_end_matches('0').trim_end_matches('.')
        );
    }
    let negative = mantissa.starts_with('-');
    let digits = mantissa.trim_start_matches('-').replace('.', "");
    let point = exponent + 1;
    let mut text = if point <= 0 {
        format!("0.{}{digits}", "0".repeat((-point) as usize))
    } else {
        let p = point as usize;
        if p == 6 {
            digits
        } else {
            format!("{}.{}", &digits[..p], &digits[p..])
        }
    };
    if text.contains('.') {
        text = text.trim_end_matches('0').trim_end_matches('.').into();
    }
    if negative {
        text.insert(0, '-');
    }
    text
}

fn insert_size(a: &[u8], b: &[u8]) -> usize {
    let rc: Vec<_> = b
        .iter()
        .rev()
        .map(|b| match b {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            _ => b'N',
        })
        .collect();
    let matches = |a: &[u8], b: &[u8]| {
        let len = a.len().min(b.len());
        let limit = 5.min((len as f64 * 0.2) as usize);
        let mut diff = 0;
        let mut i = 0;
        while i < len {
            if a[i] != b[i] {
                diff += 1;
                if diff > limit && i < 50 {
                    break;
                }
            }
            i += 1;
        }
        // Preserve fastp's acceptance of late mismatches.
        diff <= limit || (diff > limit && i > 50)
    };
    for offset in 0..a.len().saturating_sub(30) {
        let len = (a.len() - offset).min(rc.len());
        if matches(&a[offset..], &rc) {
            return if offset > 0 {
                a.len() + b.len() - len
            } else {
                len
            }
            .min(512);
        }
    }
    for offset in 0..rc.len().saturating_sub(30) {
        let len = a.len().min(rc.len() - offset);
        if matches(a, &rc[offset..]) {
            return len.min(512);
        }
    }
    512
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stream_format_boundaries() {
        for (x, expected) in [
            (0.0, "0"),
            (1.0, "1"),
            (1.234565, "1.23456"),
            (0.00001, "1e-05"),
            (0.0001, "0.0001"),
            (999999.5, "1e+06"),
            (12.34565, "12.3456"),
            (1e20, "1e+20"),
        ] {
            assert_eq!(float(x), expected);
        }
    }

    #[test]
    fn cycle_fallback_thresholds_and_kmer_reset() {
        let mut s = ReadStats::default();
        s.add(b"ATCGGNATCGG", b"5?IH~!55555");
        assert_eq!(s.q, [10, 4, 2]); // Q40 excludes DEL (127), includes '~'.
        assert_eq!(s.kmers[0b00_01_10_11_11], 2);
        assert_eq!(s.kmers.iter().sum::<u64>(), 2);
        assert_eq!(s.curve("C", true)[0], 20.0); // absent C uses mean
        assert_eq!(s.mean_length(), 11);
        s.add(b"A", b"5");
        assert_eq!(s.mean_length(), 6); // integer division
    }

    #[test]
    fn overlap_boundaries_and_late_mismatch_bug() {
        // Exactly 30 bases are never tested; forward and reverse loops are strict.
        assert_eq!(insert_size(&[b'A'; 30], &[b'T'; 30]), 512);
        assert_eq!(insert_size(&[b'A'; 31], &[b'T'; 31]), 31);
        // Excess mismatches past index 50 are accepted, even beyond limit 5.
        let mut rc = [b'A'; 100];
        rc[51..].fill(b'C');
        let b: Vec<_> = rc
            .iter()
            .rev()
            .map(|x| if *x == b'A' { b'T' } else { b'G' })
            .collect();
        assert_eq!(insert_size(&[b'A'; 100], &b), 100);
        // Early mismatches at positions 0..5 cause rejection instead.
        assert_eq!(insert_size(&[b'A'; 100], &[b'G'; 100]), 512);
    }

    #[test]
    fn duplicate_serial_order_and_prime_sequence() {
        let mut d = Duplicate::new().unwrap();
        assert_eq!(&d.primes[..3], &[10007, 20011, 30013]);
        d.add(b"ATCGATCG", b"GCTAGCTA");
        d.add(b"ATCGATCG", b"GCTAGCTA");
        assert_eq!(d.rate(), 0.5);
    }
}
