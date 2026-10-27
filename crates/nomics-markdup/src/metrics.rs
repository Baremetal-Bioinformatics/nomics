//! `picard.sam.DuplicationMetrics` accumulation, library-size estimate and the
//! metrics file with htsjdk's number formatting (spec §6).

use std::collections::BTreeMap;
use std::io::{self, Write};
use std::path::Path;

use nomics_io::bam::{Header, RecordView};

#[derive(Debug, Clone, Default, PartialEq)]
pub struct DuplicationMetrics {
    pub library: String,
    pub unpaired_reads_examined: u64,
    pub read_pairs_examined: u64,
    pub secondary_or_supplementary_rds: u64,
    pub unmapped_reads: u64,
    pub unpaired_read_duplicates: u64,
    pub read_pair_duplicates: u64,
    pub read_pair_optical_duplicates: u64,
    pub percent_duplication: f64,
    pub estimated_library_size: Option<i64>,
}

pub struct Accumulator {
    m: DuplicationMetrics,
    records: u64,
}

impl Accumulator {
    pub fn new(header: &Header) -> Self {
        let library = String::from_utf8_lossy(&header.text)
            .lines()
            .filter(|l| l.starts_with("@RG"))
            .find_map(|l| {
                l.split('\t')
                    .find_map(|f| f.strip_prefix("LB:"))
                    .map(str::to_owned)
            })
            .unwrap_or_else(|| "Unknown Library".into());
        Accumulator {
            m: DuplicationMetrics {
                library,
                ..Default::default()
            },
            records: 0,
        }
    }

    /// `addReadToLibraryMetrics` + `addDuplicateReadToMetrics` (§6.1).
    pub fn add(&mut self, raw: &[u8], duplicate: bool) -> io::Result<()> {
        let flag = RecordView::new(raw)?.flag();
        self.records += 1;
        let m = &mut self.m;
        if flag & 0x4 != 0 {
            m.unmapped_reads += 1;
        } else if flag & 0x900 != 0 {
            m.secondary_or_supplementary_rds += 1;
        } else if flag & 0x1 == 0 || flag & 0x8 != 0 {
            m.unpaired_reads_examined += 1;
        } else {
            m.read_pairs_examined += 1;
        }
        if duplicate && flag & 0x904 == 0 {
            if flag & 0x1 == 0 || flag & 0x8 != 0 {
                m.unpaired_read_duplicates += 1;
            } else {
                m.read_pair_duplicates += 1;
            }
        }
        Ok(())
    }

    /// `finalizeAndWriteMetrics` + `calculateDerivedFields`, with the
    /// library's optical duplicate pair count.
    pub fn finish(mut self, optical_pairs: u64) -> DuplicationMetrics {
        let m = &mut self.m;
        m.read_pairs_examined /= 2;
        m.read_pair_duplicates /= 2;
        m.read_pair_optical_duplicates = optical_pairs;
        m.estimated_library_size = estimate_library_size(
            (m.read_pairs_examined - m.read_pair_optical_duplicates) as i64,
            (m.read_pairs_examined - m.read_pair_duplicates) as i64,
        );
        let denom = m.unpaired_reads_examined + m.read_pairs_examined * 2;
        m.percent_duplication = if denom == 0 {
            0.0
        } else {
            (m.unpaired_read_duplicates + m.read_pair_duplicates * 2) as f64 / denom as f64
        };
        self.m
    }
}

fn f(x: f64, c: f64, n: f64) -> f64 {
    c / x - 1.0 + (-n / x).exp()
}

/// `DuplicationMetrics.estimateLibrarySize`: 40-step bisection, truncated.
pub fn estimate_library_size(read_pairs: i64, unique_read_pairs: i64) -> Option<i64> {
    let dups = read_pairs - unique_read_pairs;
    if read_pairs <= 0 || dups <= 0 {
        return None;
    }
    let (n, c) = (read_pairs as f64, unique_read_pairs as f64);
    let (mut m, mut big_m) = (1.0f64, 100.0f64);
    if unique_read_pairs >= read_pairs || f(m * c, c, n) < 0.0 {
        return None; // Picard throws here; unreachable for valid counts
    }
    while f(big_m * c, c, n) > 0.0 {
        big_m *= 10.0;
    }
    for _ in 0..40 {
        let r = (m + big_m) / 2.0;
        let u = f(r * c, c, n);
        if u == 0.0 {
            break;
        } else if u > 0.0 {
            m = r;
        } else {
            big_m = r;
        }
    }
    Some((c * (m + big_m) / 2.0) as i64)
}

/// Java `NumberFormat` with max 6 fraction digits, HALF_DOWN, no grouping,
/// trailing zeros dropped (htsjdk `FormatUtil` for doubles).
pub fn format_double(v: f64) -> String {
    if !v.is_finite() {
        return "?".into();
    }
    // Round the exact binary value: work in exact decimal via 30 digits, then
    // apply HALF_DOWN at 6 places (ties only occur at exact ...5 decimals).
    let s = format!("{:.30}", v.abs());
    let (int_part, frac) = s.split_once('.').unwrap();
    let keep = &frac[..6];
    let rest = &frac[6..];
    let round_up = match rest.as_bytes()[0] {
        b'6'..=b'9' => true,
        b'5' => rest[1..].bytes().any(|c| c != b'0'),
        _ => false,
    };
    let mut digits: Vec<u8> = int_part.bytes().chain(keep.bytes()).collect();
    if round_up {
        let mut i = digits.len();
        loop {
            if i == 0 {
                digits.insert(0, b'1');
                break;
            }
            i -= 1;
            if digits[i] == b'9' {
                digits[i] = b'0';
            } else {
                digits[i] += 1;
                break;
            }
        }
    }
    let split = digits.len() - 6;
    let int_s = String::from_utf8(digits[..split].to_vec()).unwrap();
    let frac_s = String::from_utf8(digits[split..].to_vec()).unwrap();
    let frac_s = frac_s.trim_end_matches('0');
    let neg = v < 0.0 && (!int_s.trim_start_matches('0').is_empty() || !frac_s.is_empty());
    let body = if frac_s.is_empty() {
        int_s
    } else {
        format!("{int_s}.{frac_s}")
    };
    if neg { format!("-{body}") } else { body }
}

/// Java `Double.toString` for the histogram's small integral keys.
fn java_double_key(k: f64) -> String {
    if k.fract() == 0.0 && k.abs() < 1e7 {
        format!("{k:.1}")
    } else {
        format!("{k}")
    }
}

/// Picard's duplicate-set histograms, keyed by set size (`LibraryIdGenerator`).
#[derive(Default, Debug, Clone, PartialEq)]
pub struct SetHistograms {
    pub all: BTreeMap<u64, u64>,
    pub optical: BTreeMap<u64, u64>,
    pub non_optical: BTreeMap<u64, u64>,
}

pub fn write_file(
    path: &Path,
    command_line: &str,
    m: &DuplicationMetrics,
    hist: &SetHistograms,
) -> io::Result<()> {
    let mut out = io::BufWriter::new(std::fs::File::create(path)?);
    writeln!(out, "## htsjdk.samtools.metrics.StringHeader")?;
    writeln!(out, "# {command_line}")?;
    writeln!(out, "## htsjdk.samtools.metrics.StringHeader")?;
    writeln!(out, "# Started on: (nomics)")?;
    writeln!(out)?;
    writeln!(out, "## METRICS CLASS\tpicard.sam.DuplicationMetrics")?;
    writeln!(
        out,
        "LIBRARY\tUNPAIRED_READS_EXAMINED\tREAD_PAIRS_EXAMINED\tSECONDARY_OR_SUPPLEMENTARY_RDS\tUNMAPPED_READS\tUNPAIRED_READ_DUPLICATES\tREAD_PAIR_DUPLICATES\tREAD_PAIR_OPTICAL_DUPLICATES\tPERCENT_DUPLICATION\tESTIMATED_LIBRARY_SIZE"
    )?;
    writeln!(
        out,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        m.library,
        m.unpaired_reads_examined,
        m.read_pairs_examined,
        m.secondary_or_supplementary_rds,
        m.unmapped_reads,
        m.unpaired_read_duplicates,
        m.read_pair_duplicates,
        m.read_pair_optical_duplicates,
        format_double(m.percent_duplication),
        m.estimated_library_size
            .map(|v| v.to_string())
            .unwrap_or_default()
    )?;
    writeln!(out)?;
    // ROI histogram (single library), then all_sets, optical_sets,
    // non_optical_sets in Picard's order; htsjdk drops empty histograms and
    // prints 0 for a bin a histogram lacks.
    let sets: Vec<(&str, &BTreeMap<u64, u64>)> = [
        ("all_sets", &hist.all),
        ("optical_sets", &hist.optical),
        ("non_optical_sets", &hist.non_optical),
    ]
    .into_iter()
    .filter(|(_, h)| !h.is_empty())
    .collect();
    let mut keys: std::collections::BTreeSet<u64> =
        sets.iter().flat_map(|(_, h)| h.keys().copied()).collect();
    let roi: Option<BTreeMap<u64, f64>> = m.estimated_library_size.map(|l| {
        let l = l as f64;
        let pairs = m.read_pairs_examined as f64;
        let unique = (m.read_pairs_examined - m.read_pair_duplicates) as f64;
        (1..=100u64)
            .map(|x| (x, l * (1.0 - (-(x as f64 * pairs) / l).exp()) / unique))
            .collect()
    });
    if let Some(r) = &roi {
        keys.extend(r.keys());
    }
    if roi.is_some() || !sets.is_empty() {
        writeln!(out, "## HISTOGRAM\tjava.lang.Double")?;
        let mut label = String::from("BIN");
        if roi.is_some() {
            label += "\tCoverageMult";
        }
        for (name, _) in &sets {
            label += &format!("\t{name}");
        }
        writeln!(out, "{label}")?;
        for k in keys {
            let mut row = java_double_key(k as f64);
            if let Some(r) = &roi {
                row += &format!("\t{}", format_double(r.get(&k).copied().unwrap_or(0.0)));
            }
            for (_, h) in &sets {
                row += &format!("\t{}", h.get(&k).copied().unwrap_or(0));
            }
            writeln!(out, "{row}")?;
        }
        writeln!(out)?;
    }
    out.flush()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn library_size_matches_golden() {
        assert_eq!(estimate_library_size(487639, 487639 - 69990), Some(1531991));
    }

    #[test]
    fn java_number_formatting() {
        assert_eq!(format_double(140016.0 / 975378.0), "0.143551");
        assert_eq!(format_double(1.0), "1");
        assert_eq!(format_double(0.99999997), "1");
        assert_eq!(format_double(3.668130), "3.66813");
        assert_eq!(format_double(0.0), "0");
        assert_eq!(format_double(0.0000004), "0");
        assert_eq!(format_double(1.7273825000001), "1.727383");
        assert_eq!(java_double_key(1.0), "1.0");
        assert_eq!(java_double_key(137.0), "137.0");
    }
}
