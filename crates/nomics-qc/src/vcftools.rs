//! vcftools v0.1.17 variant_file_output.cpp and vcf_entry{,_setters}.cpp.
use crate::{Record, ReferenceArch, default_float};
use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::fmt::Write;
use std::io;

// Only finite keys enter this map. Numeric equality deliberately joins +/-0,
// preserving the first key, as the upstream std::map<double,...> does.
#[derive(Clone, Copy, PartialEq)]
struct Quality(f64);
impl Eq for Quality {}
impl PartialOrd for Quality {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Quality {
    fn cmp(&self, other: &Self) -> Ordering {
        self.0.partial_cmp(&other.0).unwrap()
    }
}

#[derive(Default)]
pub(crate) struct Stats {
    count: Vec<[u64; 2]>,
    quality: BTreeMap<Quality, [u64; 2]>,
    filters: BTreeMap<String, [u64; 3]>,
}

fn substitution(reference: &str, alternate: &str) -> Option<usize> {
    if reference.len() != 1 || alternate.len() != 1 {
        return None;
    }
    let r = reference.as_bytes()[0].to_ascii_uppercase();
    let a = alternate.as_bytes()[0].to_ascii_uppercase();
    if !b"ACGT".contains(&r) || !b"ACGT".contains(&a) || r == a {
        return None;
    }
    Some(usize::from(!matches!(
        (r, a),
        (b'A', b'G') | (b'G', b'A') | (b'C', b'T') | (b'T', b'C')
    )))
}

impl Stats {
    pub fn new(samples: usize) -> Self {
        Self {
            count: vec![[0; 2]; samples * 2],
            ..Self::default()
        }
    }
    pub fn add(&mut self, record: &Record<'_>) -> io::Result<()> {
        // set_ALT drops every '.' entry, preserves duplicates, and uppercases.
        let alts: Vec<_> = record
            .alternates()
            .into_iter()
            .filter(|a| *a != ".")
            .collect();
        let first = alts.first().copied().unwrap_or(".");
        // set_FILTER sorts individual flags; summary ties reverse the whole key.
        let mut flags: Vec<_> = record.fields[6].split(';').collect();
        flags.sort_unstable();
        let filter = self.filters.entry(flags.join(";")).or_default();
        filter[0] += 1;
        if let Some(kind) = substitution(record.reference(), first) {
            filter[kind + 1] += 1;
        }
        if alts.len() != 1 {
            return Ok(());
        }
        if let Some(kind) = substitution(record.reference(), first) {
            let qual = record.fields[5].parse::<f64>().unwrap_or(-1.0);
            if !qual.is_finite() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "non-finite vcftools quality key",
                ));
            }
            self.quality.entry(Quality(qual)).or_default()[kind] += 1;
            let mut ac = 0;
            if let Some(gt_index) = record.fields.get(8).and_then(|f| {
                f.split(':')
                    .enumerate()
                    .filter_map(|(i, s)| (s == "GT").then_some(i))
                    .last()
            }) {
                for sample in record.fields.iter().skip(9) {
                    let gt = sample.split(':').nth(gt_index).unwrap_or(".");
                    let alleles: Vec<_> = gt.split(['/', '|']).collect();
                    if alleles.len() > 2 {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            "vcftools does not support polyploid genotypes",
                        ));
                    }
                    for allele in alleles {
                        if allele == "1" {
                            ac += 1;
                        } else if allele != "0" && allele != "." {
                            return Err(io::Error::new(
                                io::ErrorKind::InvalidData,
                                "invalid biallelic genotype",
                            ));
                        }
                    }
                }
            }
            // output_TsTv_by_count writes bin 2*nsamples out of bounds. The
            // authorized interpretation ignores that write; visible bins stay.
            if let Some(bin) = self.count.get_mut(ac) {
                bin[kind] += 1;
            }
        }
        Ok(())
    }
    pub fn render(self, arch: ReferenceArch) -> [String; 3] {
        let mut count = "ALT_ALLELE_COUNT\tN_Ts\tN_Tv\tTs/Tv\n".to_owned();
        for (i, [ts, tv]) in self.count.iter().enumerate() {
            writeln!(
                count,
                "{i}\t{ts}\t{tv}\t{}",
                arch.ratio(*ts as f64, *tv as f64)
            )
            .unwrap();
        }
        let mut filter = "FILTER\tN_VARIANTS\tN_Ts\tN_Tv\tTs/Tv\n".to_owned();
        let mut rows: Vec<_> = self.filters.into_iter().collect();
        rows.sort_unstable_by(|a, b| b.1[0].cmp(&a.1[0]).then_with(|| b.0.cmp(&a.0)));
        for (name, [n, ts, tv]) in rows {
            writeln!(
                filter,
                "{name}\t{n}\t{ts}\t{tv}\t{}",
                arch.ratio(ts as f64, tv as f64)
            )
            .unwrap();
        }
        let mut qual = "QUAL_THRESHOLD\tN_Ts_LT_QUAL_THRESHOLD\tN_Tv_LT_QUAL_THRESHOLD\tTs/Tv_LT_QUAL_THRESHOLD\tN_Ts_GT_QUAL_THRESHOLD\tN_Tv_GT_QUAL_THRESHOLD\tTs/Tv_GT_QUAL_THRESHOLD\n".to_owned();
        let mut upper = self
            .quality
            .values()
            .fold([0u64; 2], |sum, bin| [sum[0] + bin[0], sum[1] + bin[1]]);
        let mut below = [0u64; 2];
        for (Quality(q), bin) in self.quality {
            // Source upper vector is 1-based. ui+1 excludes the current bin;
            // its otherwise invalid N+1 element is explicitly zero here.
            upper[0] -= bin[0];
            upper[1] -= bin[1];
            writeln!(
                qual,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                default_float(q),
                default_float(below[0] as f64),
                default_float(below[1] as f64),
                arch.ratio(below[0] as f64, below[1] as f64),
                default_float(upper[0] as f64),
                default_float(upper[1] as f64),
                arch.ratio(upper[0] as f64, upper[1] as f64)
            )
            .unwrap();
            below[0] += bin[0];
            below[1] += bin[1];
        }
        [count, qual, filter]
    }
}
