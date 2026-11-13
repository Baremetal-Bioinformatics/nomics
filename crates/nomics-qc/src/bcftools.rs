//! bcftools 1.23.1 vcfstats.c, default one-file mode (no sample selection).
//! Typing follows htslib 1.23.1 vcf.c; histograms follow bcftools dist.c.
use crate::{Header, Record};
use std::collections::BTreeMap;
use std::fmt::Write;
use std::io;
#[path = "bcftools_headers.rs"]
mod headers;

const SNP: u8 = 1;
const MNP: u8 = 2;
const INDEL: u8 = 4;
const OTHER: u8 = 8;
const BND: u8 = 16;
const OVERLAP: u8 = 32;

/// The legacy variant mask and signed indel length. Only indels use the length.
fn variant_type(reference: &str, alternate: &str) -> (u8, i32) {
    if alternate == "*" {
        return (OVERLAP, 0);
    }
    let r = reference.as_bytes();
    let a = alternate.as_bytes();
    if r.len() == 1 && a.len() == 1 {
        return (
            if a == b"." || r == a || a == b"X" {
                0
            } else {
                SNP
            },
            0,
        );
    }
    if a.starts_with(b"<") {
        return (
            if matches!(alternate, "<X>" | "<*>" | "<NON_REF>") {
                0
            } else {
                OTHER
            },
            0,
        );
    }
    if matches!(a.first(), Some(b'[' | b']')) {
        return (BND, 0);
    }
    let common = r
        .iter()
        .zip(a)
        .take_while(|(x, y)| x.eq_ignore_ascii_case(y))
        .count();
    let rr = &r[common..];
    let aa = &a[common..];
    if matches!(aa.last(), Some(b'[' | b']')) {
        return (BND, 0);
    }
    match (rr.is_empty(), aa.is_empty()) {
        (true, true) => return (0, 0),
        (true, false) | (false, true) => return (INDEL, a.len() as i32 - r.len() as i32),
        _ => (),
    }
    // Trim equal suffixes while leaving at least one base on either side.
    let suffix = rr
        .iter()
        .rev()
        .take(rr.len() - 1)
        .zip(aa.iter().rev().take(aa.len() - 1))
        .take_while(|(x, y)| x.eq_ignore_ascii_case(y))
        .count();
    let nr = rr.len() - suffix;
    let na = aa.len() - suffix;
    if nr == 1 && na == 1 {
        return (SNP, 0);
    }
    if nr == 1 || na == 1 {
        return (
            if rr[nr - 1].eq_ignore_ascii_case(&aa[na - 1]) {
                INDEL
            } else {
                OTHER
            },
            na as i32 - nr as i32,
        );
    }
    (if nr == na { MNP } else { OTHER }, 0)
}

fn base(b: u8) -> Option<usize> {
    b"ACGT".iter().position(|x| *x == b.to_ascii_uppercase())
}

// dist_init(5): exact through 100000; later levels have decimal steps 10,100,...
fn quality_bin(q: f32) -> u32 {
    let value = if q.is_nan() || q < 0.0 {
        0
    } else {
        1 + (q * 10.0) as u32
    };
    if value <= 100000 {
        return value;
    }
    let power = value.ilog10();
    let level = power - 4;
    100000 + 90000 * (level - 1) + (value - 10u32.pow(power)) / 10u32.pow(level)
}
fn quality_label(bin: u32) -> String {
    if bin == 0 {
        return ".".into();
    }
    let beg = if bin < 100000 {
        bin
    } else {
        let level = (bin - 100000) / 90000 + 1;
        10u32.pow(level + 4) + 10u32.pow(level) * (bin - 100000 - 90000 * (level - 1))
    };
    // vcfstats stores the double expression into a float, then printf promotes.
    format!("{:.1}", (0.1 * (beg - 1) as f64) as f32)
}

#[derive(Default)]
pub(crate) struct Stats {
    summary: [u64; 8], // records,noalt,snp,mnp,indel,other,multi,multi-SNP
    af: Vec<[u64; 4]>, // SNP alleles, ts, tv, indel alleles (all repeat n/a)
    first: [u64; 2],
    qual: BTreeMap<u32, [u64; 3]>,
    indels: BTreeMap<i32, u64>,
    substitutions: [u64; 16],
    depth: BTreeMap<i32, u64>,
}

impl Stats {
    pub fn new(samples: usize) -> Self {
        Self {
            af: vec![[0; 4]; 101.max(samples + 1)],
            ..Self::default()
        }
    }
    fn af_indices(
        &self,
        record: &Record<'_>,
        header: &Header,
        nalt: usize,
    ) -> io::Result<Vec<usize>> {
        let unknown = || vec![0; nalt];
        if !header.integer_info.contains("AC") || !header.integer_info.contains("AN") {
            return Ok(unknown());
        }
        let Some(an) = record
            .info("AN")
            .and_then(|s| s.parse::<i32>().ok())
            .filter(|n| *n >= 0)
        else {
            return Ok(unknown());
        };
        let Some(ac) = record.info("AC") else {
            return Ok(unknown());
        };
        let ac: Vec<_> = ac
            .split(',')
            .map(|s| if s == "." { Ok(-128) } else { s.parse::<i32>() })
            .collect();
        if ac.len() != nalt || ac.iter().any(|s| s.is_err()) {
            return Ok(unknown());
        }
        let ac: Vec<_> = ac.into_iter().map(Result::unwrap).collect();
        if ac.iter().map(|n| i64::from(*n)).sum::<i64>() > i64::from(an) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "AC sum exceeds AN",
            ));
        }
        Ok(ac
            .into_iter()
            .map(|n| {
                if n == 1 {
                    0
                } else if an == 0 {
                    1
                } else {
                    let af = (n as f32 / an as f32).clamp(0.0, 1.0);
                    1 + (af * (self.af.len() - 2) as f32) as usize
                }
            })
            .collect())
    }
    pub fn add(&mut self, record: &Record<'_>, header: &Header) -> io::Result<()> {
        let alts = record.alternates();
        // htslib stores ALT '.' as no alternate, including for multibase REF.
        let alts: Vec<_> = if record.fields[4] == "." {
            Vec::new()
        } else {
            alts
        };
        let variants: Vec<_> = alts
            .iter()
            .map(|a| variant_type(record.reference(), a))
            .collect();
        let mask = variants.iter().fold(0, |s, t| s | t.0);
        self.summary[0] += 1;
        self.summary[1] += u64::from(mask == 0);
        for (idx, bit) in [(2, SNP), (3, MNP), (4, INDEL), (5, OTHER)] {
            self.summary[idx] += u64::from(mask & bit != 0);
        }
        if alts.len() > 1 {
            self.summary[6] += 1;
            self.summary[7] += u64::from(mask == SNP);
        }
        let iaf = self.af_indices(record, header, alts.len())?;
        let q = record.fields[5].parse::<f64>().unwrap_or(f64::NAN) as f32;
        if q.is_infinite() || q >= (i32::MAX as f32 / 10.0) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "QUAL exceeds upstream histogram integer range",
            ));
        }
        let iq = quality_bin(q);
        if mask & INDEL != 0 {
            self.qual.entry(iq).or_default()[2] += 1;
        }
        for (i, ((kind, len), alt)) in variants.iter().zip(&alts).enumerate() {
            if *kind == INDEL {
                self.af[iaf[i]][3] += 1;
                *self.indels.entry((*len).clamp(-60, 60)).or_default() += 1;
            }
            if *kind == SNP
                && let (Some(r), Some(a)) = (
                    record.reference().bytes().next().and_then(base),
                    alt.bytes().next().and_then(base),
                )
                && r != a
            {
                self.substitutions[r * 4 + a] += 1;
                self.af[iaf[i]][0] += 1;
                let k = usize::from(r.abs_diff(a) != 2);
                self.af[iaf[i]][k + 1] += 1;
                if i == 0 {
                    self.first[k] += 1;
                    self.qual.entry(iq).or_default()[k] += 1;
                }
            }
        }
        if header.integer_info.contains("DP")
            && let Some(dp) = record.info("DP").and_then(|s| {
                if !s.contains(',') {
                    Some(s.parse::<i32>().unwrap_or(i32::MIN))
                } else {
                    None
                }
            })
        {
            *self
                .depth
                .entry(if dp < 0 {
                    -1
                } else if dp > 500 {
                    501
                } else {
                    dp
                })
                .or_default() += 1;
        }
        Ok(())
    }
    pub fn render(mut self, filename: &str, samples: usize) -> String {
        let mut out = format!(
            "# This file was produced by bcftools stats (1.23.1+htslib-1.23.1) and can be plotted using plot-vcfstats.\n# The command line was:\tbcftools stats  {filename}\n#\n# Definition of sets:\n# ID\t[2]id\t[3]tab-separated file names\nID\t0\t{filename}\n"
        );
        out.push_str(headers::SN);
        writeln!(out, "SN\t0\tnumber of samples:\t{samples}").unwrap();
        for (key, n) in [
            "records",
            "no-ALTs",
            "SNPs",
            "MNPs",
            "indels",
            "others",
            "multiallelic sites",
            "multiallelic SNP sites",
        ]
        .iter()
        .zip(self.summary)
        {
            writeln!(out, "SN\t0\tnumber of {key}:\t{n}").unwrap();
        }
        out.push_str(headers::TSTV);
        let ts: u64 = self.af.iter().map(|a| a[1]).sum();
        let tv: u64 = self.af.iter().map(|a| a[2]).sum();
        let [ts1, tv1] = self.first;
        // Both ratios are computed as float, then promoted for printf %.2f.
        let ratio = |a, b| if b == 0 { 0.0 } else { a as f32 / b as f32 };
        writeln!(
            out,
            "TSTV\t0\t{ts}\t{tv}\t{:.2}\t{ts1}\t{tv1}\t{:.2}",
            ratio(ts, tv),
            ratio(ts1, tv1)
        )
        .unwrap();
        out.push_str(headers::SIS);
        let [sn, ts, tv, ins] = self.af[0];
        writeln!(out, "SiS\t0\t1\t{sn}\t{ts}\t{tv}\t{ins}\t0\t0\t{ins}").unwrap();
        for i in 0..4 {
            self.af[1][i] += self.af[0][i];
        }
        out.push_str(headers::AF);
        for (i, [sn, ts, tv, ins]) in self.af.iter().enumerate().skip(1) {
            if sn + ts + tv + ins == 0 {
                continue;
            }
            writeln!(
                out,
                "AF\t0\t{:.6}\t{sn}\t{ts}\t{tv}\t{ins}\t0\t0\t{ins}",
                (i - 1) as f64 / (self.af.len() - 1) as f64
            )
            .unwrap();
        }
        out.push_str(headers::QUAL);
        for (i, [ts, tv, ins]) in self.qual {
            writeln!(
                out,
                "QUAL\t0\t{}\t{}\t{ts}\t{tv}\t{ins}",
                quality_label(i),
                ts + tv
            )
            .unwrap();
        }
        out.push_str(headers::IDD);
        for (length, n) in self.indels {
            writeln!(out, "IDD\t0\t{length}\t{n}\t0\t.").unwrap();
        }
        out.push_str(headers::ST);
        for (r, x) in b"ACGT".iter().enumerate() {
            for (a, y) in b"ACGT".iter().enumerate() {
                if r == a {
                    continue;
                }
                writeln!(
                    out,
                    "ST\t0\t{}>{}\t{}",
                    *x as char,
                    *y as char,
                    self.substitutions[r * 4 + a]
                )
                .unwrap();
            }
        }
        out.push_str(headers::DP);
        let total: u64 = self.depth.values().sum();
        for (dp, n) in self.depth {
            let label = match dp {
                -1 => "<0".into(),
                501 => ">500".into(),
                _ => dp.to_string(),
            };
            writeln!(
                out,
                "DP\t0\t{label}\t0\t0.000000\t{n}\t{:.6}",
                n as f64 * 100.0 / total as f64
            )
            .unwrap();
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn htslib_variant_classification() {
        // htslib 1.23.1 bcf_set_variant_type / bcf_get_variant_types: prefix
        // and suffix trimming, legacy mask, and reserved reference alleles.
        for (r, a, t, len) in [
            ("A", "G", SNP, 0),
            ("A", "a", SNP, 0),
            ("AA", "aa", 0, 0),
            ("A", "N", SNP, 0),
            ("A", "X", 0, 0),
            ("A", ".", 0, 0),
            ("A", "*", OVERLAP, 0),
            ("A", "<NON_REF>", 0, 0),
            ("A", "<*>", 0, 0),
            ("A", "<X>", 0, 0),
            ("A", "<DEL>", OTHER, 0),
            ("CC", "TT", MNP, 0),
            ("ACT", "TCGA", OTHER, 0),
            ("A", "ATG", INDEL, 2),
            ("ATG", "A", INDEL, -2),
            ("GAC", "GTC", SNP, 0),
            ("GAC", "GTTC", OTHER, 2),
            ("AC", "TAC", INDEL, 1),
            ("TAC", "AC", INDEL, -1),
            ("A", "]chr1:2]A", BND, 0),
            ("A", "A[chr1:2[", BND, 0),
        ] {
            let actual = variant_type(r, a);
            assert_eq!(actual.0, t, "{r}>{a}");
            if t == INDEL {
                assert_eq!(actual.1, len, "{r}>{a}");
            }
        }
    }
    #[test]
    fn decimal_histogram_boundaries() {
        // bcftools dist_insert/dist_get, dist_init(5), and do_snp_stats.
        for (q, bin, label) in [
            (f32::NAN, 0, "."),
            (-1.0, 0, "."),
            (0.0, 1, "0.0"),
            (9999.9, 100000, "9999.9"),
            (10000.0, 100000, "9999.9"),
            (10000.9, 100001, "10000.9"),
            (100000.0, 190000, "99999.9"),
        ] {
            assert_eq!(quality_bin(q), bin);
            assert_eq!(quality_label(bin), label);
        }
    }
}
