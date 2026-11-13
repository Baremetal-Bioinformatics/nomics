//! GATK 4.6.2.0 walkers/contamination/{CalculateContamination,ContaminationModel,
//! ContaminationRecord,MinorAlleleFractionRecord}.java.
use crate::{
    math::{
        self, binomial::binomial_probability, brent::brent_optimizer, javafmt::double_to_string,
        stream_sum,
    },
    pileup_summary::{self, PileupSummary},
    segmentation,
};
use std::{
    collections::HashSet,
    fs::File,
    io::{self, BufWriter, Write},
    path::Path,
};

pub fn filter_sites_by_coverage(
    all: Vec<PileupSummary>,
    low: f64,
    high: f64,
) -> Vec<PileupSummary> {
    let covered: Vec<_> = all.into_iter().filter(|s| s.total() > 10).collect();
    let cov: Vec<_> = covered.iter().map(|s| f64::from(s.total())).collect();
    let min = math::percentile(&cov, 50.0) * low;
    let max = math::mean(&cov) * high;
    covered
        .into_iter()
        .filter(|s| f64::from(s.total()) > min && f64::from(s.total()) < max)
        .collect()
}
fn error_rate(sites: &[PileupSummary]) -> f64 {
    let total = sites.iter().fold(0i32, |n, s| n.wrapping_add(s.total()));
    let other = sites
        .iter()
        .fold(0i32, |n, s| n.wrapping_add(s.other_alt_count));
    1.5 * (f64::from(other) / f64::from(total))
}
fn genotype_likelihoods(s: &PileupSummary, c: f64, error: f64, maf: f64) -> [f64; 4] {
    let f = s.allele_frequency;
    let priors = [(1.0 - f) * (1.0 - f), f * (1.0 - f), f * (1.0 - f), f * f];
    let afs = [error / 3.0, maf, 1.0 - maf, 1.0 - error];
    std::array::from_fn(|i| {
        priors[i]
            * binomial_probability(
                s.alt_count + s.ref_count,
                s.alt_count,
                (1.0 - c) * afs[i] + c * f,
            )
    })
}
fn segment_log_likelihood(
    sites: &[PileupSummary],
    segment: &[usize],
    c: f64,
    error: f64,
    maf: f64,
) -> f64 {
    stream_sum(segment.iter().map(|&i| {
        math::java_log::log(
            genotype_likelihoods(&sites[i], c, error, maf)
                .iter()
                .sum::<f64>(),
        )
    }))
}
struct Model<'a> {
    sites: &'a [PileupSummary],
    segments: Vec<Vec<usize>>,
    mafs: Vec<f64>,
    contamination: f64,
    error: f64,
}
impl<'a> Model<'a> {
    fn new(sites: &'a [PileupSummary]) -> io::Result<Self> {
        let error = error_rate(sites);
        let segments = segmentation::find_segments(sites);
        let mut mafs = vec![0.5; segments.len()];
        let mut c = 0.0;
        for _ in 0..3 {
            for (segment, maf) in segments.iter().zip(&mut mafs) {
                *maf = brent_optimizer(
                    |m| segment_log_likelihood(sites, segment, c, error, m),
                    0.1,
                    0.5,
                    0.4,
                    0.01,
                    0.01,
                    20,
                )?
                .0;
            }
            let total: usize = segments.iter().map(Vec::len).sum();
            let mut threshold = 0.4;
            let mut selected = (0..segments.len()).collect::<Vec<_>>();
            while threshold > 0.0 {
                let indices: Vec<_> = (0..segments.len())
                    .filter(|&i| mafs[i] > threshold)
                    .collect();
                let num: usize = indices.iter().map(|&i| segments[i].len()).sum();
                if num as f64 / total as f64 > 0.25 {
                    selected = indices;
                    break;
                }
                threshold -= 0.04;
            }
            let mut best = None;
            for initial in [0.02, 0.05, 0.1, 0.2] {
                // IndexRange.sum is plain addition over segments; each segment is a compensated stream sum.
                let optimum = brent_optimizer(
                    |c| {
                        selected
                            .iter()
                            .map(|&i| {
                                segment_log_likelihood(sites, &segments[i], c, error, mafs[i])
                            })
                            .sum()
                    },
                    0.0,
                    0.5,
                    initial,
                    1e-4,
                    1e-4,
                    30,
                )?;
                if best.is_none_or(|b: (f64, f64)| optimum.1 > b.1) {
                    best = Some(optimum);
                }
            }
            c = best.unwrap().0;
        }
        Ok(Self {
            sites,
            segments,
            mafs,
            contamination: c,
            error,
        })
    }
    fn estimate(&self, tumor: &[PileupSummary], threshold: f64, strategy: usize) -> (f64, f64) {
        let loci: HashSet<(&str, i32)> = if strategy == 2 {
            let candidates: Vec<_> = tumor.iter().filter(|s| s.alt_fraction() < 0.15).collect();
            let fractions: Vec<_> = candidates.iter().map(|s| s.alt_fraction()).collect();
            let p = math::percentile(&fractions, 90.0);
            let cutoff = if p.is_nan() { p } else { 0.1_f64.max(p) };
            candidates
                .into_iter()
                .filter(|s| s.alt_fraction() <= cutoff)
                .map(|s| (s.contig.as_str(), s.position))
                .collect()
        } else {
            let mut loci = HashSet::new();
            let genotype = if strategy == 0 { 3 } else { 0 };
            for (segment, &maf) in self.segments.iter().zip(&self.mafs) {
                if maf <= threshold {
                    continue;
                }
                for &i in segment {
                    let s = &self.sites[i];
                    let gl = genotype_likelihoods(s, self.contamination, self.error, maf);
                    if gl[genotype] / gl.iter().sum::<f64>() > 0.5 {
                        loci.insert((s.contig.as_str(), s.position));
                    }
                }
            }
            loci
        };
        let homs: Vec<_> = tumor
            .iter()
            .filter(|s| loci.contains(&(s.contig.as_str(), s.position)))
            .collect();
        let opposite_frequency = |s: &&PileupSummary| {
            if strategy == 0 {
                1.0 - s.allele_frequency
            } else {
                s.allele_frequency
            }
        };
        let total: i64 = homs.iter().map(|s| i64::from(s.total())).sum();
        let opposite: i64 = homs
            .iter()
            .map(|s| {
                i64::from(if strategy == 0 {
                    s.ref_count
                } else {
                    s.alt_count
                })
            })
            .sum();
        let errors = math::java_round(total as f64 * error_rate(tumor) / 3.0);
        let q = stream_sum(
            homs.iter()
                .map(|s| f64::from(s.total()) * opposite_frequency(s)),
        );
        let estimate = (opposite - errors).max(0) as f64 / q;
        let a = stream_sum(
            homs.iter()
                .map(|s| opposite_frequency(s) * f64::from(s.total())),
        );
        let b = stream_sum(homs.iter().map(|s| {
            let f = opposite_frequency(s);
            let d = f64::from(s.total());
            f * (1.0 - f) * (d * d)
        }));
        let error = |c: f64| {
            if homs.is_empty() {
                1.0
            } else {
                (a * c * (1.0 - c) + b * c * c).sqrt() / q
            }
        };
        let upper = math::binary_search_zero(|c| c - error(c) - estimate, estimate, 1.0);
        let lower = math::binary_search_zero(|c| c + error(c) - estimate, 0.0, estimate);
        let se = if let (Some(u), Some(l)) = (upper, lower) {
            (u - estimate).max(estimate - l)
        } else {
            error(estimate)
        };
        (
            if estimate.is_nan() {
                estimate
            } else {
                estimate.min(1.0)
            },
            se,
        )
    }
    fn calculate_from_homs(&self, tumor: &[PileupSummary]) -> (f64, f64) {
        let mut threshold = 0.4;
        while threshold >= 0.0 {
            let strategy = if threshold > 0.25 {
                0
            } else if threshold > 0.20 {
                1
            } else {
                2
            };
            let result = self.estimate(tumor, threshold, strategy);
            if !result.0.is_nan() && result.1 < result.0 * 0.2 + 0.001 {
                return result;
            }
            threshold -= 0.04;
        }
        let result = self.estimate(tumor, 0.0, 2);
        if result.0.is_nan() {
            (0.0, 1.0)
        } else {
            result
        }
    }
    fn write_segmentation(&self, sample: &str, path: &Path) -> io::Result<()> {
        let mut w = BufWriter::new(File::create(path)?);
        writeln!(
            w,
            "#<METADATA>SAMPLE={sample}\ncontig\tstart\tend\tminor_allele_fraction"
        )?;
        for (s, &maf) in self.segments.iter().zip(&self.mafs) {
            let first = &self.sites[s[0]];
            let last = &self.sites[*s.last().unwrap()];
            writeln!(
                w,
                "{}\t{}\t{}\t{}",
                first.contig,
                first.position,
                last.position,
                double_to_string(maf)
            )?;
        }
        w.flush()
    }
}
/// Execute CalculateContamination. Tumor-only mode reuses its genotyping model.
pub fn calculate(
    input: &Path,
    normal: Option<&Path>,
    output: &Path,
    segmentation: Option<&Path>,
    low: f64,
    high: f64,
) -> io::Result<()> {
    let (sample, all) = pileup_summary::read(input)?;
    let tumor = filter_sites_by_coverage(all, low, high);
    let normal = normal
        .map(|p| pileup_summary::read(p).map(|t| filter_sites_by_coverage(t.1, low, high)))
        .transpose()?;
    let model = Model::new(normal.as_deref().unwrap_or(&tumor))?;
    if let Some(path) = segmentation {
        if normal.is_some() {
            Model::new(&tumor)?.write_segmentation(&sample, path)?;
        } else {
            model.write_segmentation(&sample, path)?;
        }
    }
    let (c, e) = model.calculate_from_homs(&tumor);
    let mut w = BufWriter::new(File::create(output)?);
    writeln!(
        w,
        "sample\tcontamination\terror\n{sample}\t{}\t{}",
        double_to_string(c),
        double_to_string(e)
    )?;
    w.flush()
}

#[cfg(test)]
mod tests {
    use super::*;
    fn row(position: i32, reference: i32, alt: i32, other: i32) -> PileupSummary {
        PileupSummary {
            contig: "chr1".into(),
            position,
            ref_count: reference,
            alt_count: alt,
            other_alt_count: other,
            allele_frequency: 0.1,
        }
    }
    #[test]
    fn coverage_thresholds_are_strict_and_error_sum_wraps() {
        // Surviving depths [11,20,40,80], median=30: 0.5*median=15.
        let rows = vec![
            row(1, 10, 0, 0),
            row(2, 11, 0, 0),
            row(3, 20, 0, 0),
            row(4, 40, 0, 0),
            row(5, 80, 0, 0),
        ];
        assert_eq!(
            filter_sites_by_coverage(rows, 0.5, 3.)
                .iter()
                .map(|s| s.position)
                .collect::<Vec<_>>(),
            vec![3, 4, 5]
        );
        assert!(filter_sites_by_coverage(vec![row(1, 20, 0, 0)], 1., 3.).is_empty());
        assert!(filter_sites_by_coverage(vec![row(1, 20, 0, 0)], 0.5, 1.).is_empty());
        assert_eq!(error_rate(&[row(1, 99, 0, 1)]), 0.015);
        let expected = 1.5 * (2.0 / -2147483646_f64);
        assert_eq!(
            error_rate(&[row(1, 1073741824, 0, 1), row(2, 1073741824, 0, 1)]),
            expected
        );
    }
}
