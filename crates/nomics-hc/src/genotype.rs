//! Single-sample, diploid per-event genotyping.
//!
//! This is an operation-order-preserving port of the GATK 4.6.2.0 path from
//! `GenotypeLikelihoodCalculator` through `GenotypingEngine`. The likelihood
//! matrix is indexed by allele, then read, and genotypes use VCF order.

use crate::jmath;
use crate::tables::{
    approximate_log10_sum_log10, java_round, log10_e, log10_sum_log10, max_element_index,
};

const REF_PSEUDOCOUNT: f64 = 10.0;
const SNP_PSEUDOCOUNT: f64 = 0.01;
const INDEL_PSEUDOCOUNT: f64 = 0.00125;
const ALLELE_COUNT_CONVERGENCE: f64 = 0.1;
const CALLING_CONFIDENCE: f64 = 30.0;
const AF_EPSILON: f64 = 1.0e-10;
const SPAN_DEL: &[u8] = b"*";

/// Sentinel used only if GATK's hom-ref/non-informative no-call rule fires.
pub const NO_CALL_ALLELE: usize = usize::MAX;

/// Calculate diploid genotype log10 likelihoods in VCF genotype order.
pub fn genotype_likelihoods(matrix: &[Vec<f64>]) -> Vec<f64> {
    assert!(!matrix.is_empty(), "the likelihood matrix needs an allele");
    let read_count = matrix[0].len();
    assert!(
        matrix.iter().all(|row| row.len() == read_count),
        "likelihood matrix rows have unequal lengths"
    );

    let mut result = Vec::with_capacity(matrix.len() * (matrix.len() + 1) / 2);
    let log10_two = jmath::log10(2.0);
    for high in 0..matrix.len() {
        for low in 0..=high {
            let mut likelihood = 0.0;
            if low == high {
                for &read_likelihood in &matrix[low] {
                    likelihood += read_likelihood;
                }
            } else {
                let log10_one = jmath::log10(1.0);
                for (&low_likelihood, &high_likelihood) in matrix[low].iter().zip(&matrix[high]) {
                    likelihood += approximate_log10_sum_log10(
                        low_likelihood + log10_one,
                        high_likelihood + log10_one,
                    );
                }
                likelihood -= read_count as f64 * log10_two;
            }
            result.push(likelihood);
        }
    }
    result
}

/// Convert log10 genotype likelihoods to normalized integer PLs using htsjdk
/// `GenotypeLikelihoods.GLsToPLs` semantics.
pub fn pls(gls: &[f64]) -> Vec<i32> {
    assert!(!gls.is_empty(), "the likelihood vector must not be empty");
    let mut adjust = f64::NEG_INFINITY;
    for &gl in gls {
        adjust = java_max(adjust, gl);
    }
    gls.iter()
        .map(|&gl| {
            let phred = java_min(-10.0 * (gl - adjust), i32::MAX as f64);
            java_round(phred) as i32
        })
        .collect()
}

/// A genotyped HaplotypeCaller event.
#[derive(Clone, Debug, PartialEq)]
pub struct Call {
    pub log10_p_error: f64,
    /// Indices into the input allele list, reference first.
    pub alleles_kept: Vec<usize>,
    /// Indices into `alleles_kept`; `NO_CALL_ALLELE` represents `.`.
    pub gt: [usize; 2],
    pub pl: Vec<i32>,
    pub gq: i32,
    pub mleac: Vec<i32>,
    pub mleaf: Vec<f64>,
}

#[derive(Clone, Copy, Debug)]
struct Deletion {
    contig: usize,
    start: i64,
    end: i64,
}

/// Stateful per-event genotyper. Emitted deletions are retained so a later
/// `*` allele is accepted only while it is actually spanned by one.
#[derive(Default, Debug)]
pub struct GenotypingEngine {
    upstream_deletions: Vec<Deletion>,
    pub(crate) emit_gvcf: bool,
}

#[derive(Debug)]
struct AfResult {
    allele_counts: Vec<i32>,
    log10_p_no_variant: f64,
    log10_p_absent: Vec<f64>,
}

impl GenotypingEngine {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn gvcf() -> Self {
        Self {
            emit_gvcf: true,
            ..Self::default()
        }
    }

    /// Apply HaplotypeCaller's default single-sample genotyping policy.
    pub fn calculate(
        &mut self,
        contig: usize,
        start: i64,
        _end: i64,
        alleles: &[Vec<u8>],
        input_pls: &[i32],
    ) -> Option<Call> {
        assert!(alleles.len() > 1, "genotyping requires an alternate allele");
        assert_eq!(
            input_pls.len(),
            alleles.len() * (alleles.len() + 1) / 2,
            "PL count does not match diploid allele count"
        );

        // More than --max-alternate-alleles (6) alts: the AF calculation runs
        // on the most likely alleles only (`calculateGenotypes`). `used` maps
        // the reduced allele list to indices into `alleles`.
        let confidence = if self.emit_gvcf {
            -0.0
        } else {
            CALLING_CONFIDENCE
        };
        let (used, reduced_pls) = most_likely_alleles(
            input_pls,
            alleles.len(),
            alleles.iter().position(|a| a == b"<NON_REF>"),
        );
        let used_alleles: Vec<Vec<u8>> = used.iter().map(|&i| alleles[i].clone()).collect();
        let af = calculate_af(&used_alleles, &reduced_pls);
        let mut output_alts = Vec::new();
        let mut mleac = Vec::new();
        let mut site_is_monomorphic = true;

        for (reduced_index, &allele_index) in used.iter().enumerate().skip(1) {
            let allele = &alleles[allele_index];
            let plausible = af.log10_p_absent[reduced_index] + AF_EPSILON < confidence * -0.1;
            let spurious_span_del =
                allele == SPAN_DEL && !self.is_covered_by_deletion(contig, start);
            site_is_monomorphic &= !(plausible && !spurious_span_del);
            if (plausible || self.emit_gvcf) && !spurious_span_del {
                output_alts.push(allele_index);
                mleac.push(af.allele_counts[reduced_index - 1]);
            }
        }

        let log10_confidence = if !site_is_monomorphic || self.emit_gvcf {
            af.log10_p_no_variant + 0.0
        } else {
            log10_one_minus_pow10(af.log10_p_no_variant) + 0.0
        };
        let phred_confidence = (-10.0 * log10_confidence) + 0.0;
        if (site_is_monomorphic || phred_confidence < confidence)
            && output_alts
                .first()
                .is_none_or(|&a| alleles[a] != b"<NON_REF>")
        {
            return None;
        }
        if output_alts.len() == 1 && alleles[output_alts[0]] == SPAN_DEL {
            return None;
        }

        let mut kept = Vec::with_capacity(output_alts.len() + 1);
        kept.push(0);
        kept.extend(output_alts);

        let (gt, call_pl, gq) = subset_and_call(input_pls, &kept);
        let called_alleles = gt.iter().filter(|&&a| a != NO_CALL_ALLELE).count() as f64;
        let mleaf = mleac
            .iter()
            .map(|&count| java_min(1.0, count as f64 / called_alleles))
            .collect();

        self.record_deletions(contig, start, alleles, &kept);

        Some(Call {
            log10_p_error: log10_confidence,
            alleles_kept: kept,
            gt,
            pl: call_pl,
            gq,
            mleac,
            mleaf,
        })
    }

    fn is_covered_by_deletion(&self, contig: usize, start: i64) -> bool {
        self.upstream_deletions.iter().any(|deletion| {
            deletion.contig == contig && deletion.start < start && start <= deletion.end
        })
    }

    fn record_deletions(&mut self, contig: usize, start: i64, alleles: &[Vec<u8>], kept: &[usize]) {
        self.upstream_deletions
            .retain(|deletion| deletion.contig == contig && deletion.end >= start);

        let reference_length = alleles[0].len() as i64;
        for &allele_index in kept {
            let deletion_size = reference_length - allele_length(&alleles[allele_index]) as i64;
            if deletion_size > 0 {
                self.upstream_deletions.push(Deletion {
                    contig,
                    start,
                    end: start + deletion_size,
                });
            }
        }
    }
}

fn allele_length(allele: &[u8]) -> usize {
    if allele.starts_with(b"<") {
        0
    } else {
        allele.len()
    }
}

const MAX_ALTERNATE_ALLELES: usize = 6;

/// `AlleleSubsettingUtils.calculateMostLikelyAlleles` + `subsetAlleles` for one
/// diploid sample. Each alt in the most likely genotype scores
/// `|GL(best) - GL(hom-ref)|` (others 0); the top 6 by a stable descending sort
/// are kept, in their original order, and the PLs are subset, rescaled to a
/// maximum of 0 and rounded back to PLs. Returns all alleles unchanged when
/// there are at most 6 alts.
fn most_likely_alleles(
    input_pls: &[i32],
    allele_count: usize,
    non_ref: Option<usize>,
) -> (Vec<usize>, Vec<i32>) {
    if allele_count - 1 - usize::from(non_ref.is_some()) <= MAX_ALTERNATE_ALLELES {
        return ((0..allele_count).collect(), input_pls.to_vec());
    }
    let gls: Vec<f64> = input_pls.iter().map(|&pl| pl as f64 / -10.0).collect();
    let best = max_element_index(&gls);
    let diff = (gls[best] - gls[0]).abs();
    let (low, high) = diploid_genotypes(allele_count)
        .nth(best)
        .expect("genotype index in range");
    let mut sums = vec![0.0; allele_count];
    for (allele, sum) in sums.iter_mut().enumerate().skip(1) {
        if allele == low || allele == high {
            *sum += diff;
        }
    }
    let mut ranked: Vec<usize> = (1..allele_count).filter(|a| Some(*a) != non_ref).collect();
    // Comparator.comparingDouble(...).reversed() on a stable sort.
    ranked.sort_by(|&a, &b| sums[b].total_cmp(&sums[a]));
    let mut used: Vec<usize> = std::iter::once(0)
        .chain(ranked.into_iter().take(MAX_ALTERNATE_ALLELES))
        .chain(non_ref)
        .collect();
    used.sort_unstable();
    let mut subset = Vec::with_capacity(used.len() * (used.len() + 1) / 2);
    for (low, high) in diploid_genotypes(used.len()) {
        subset.push(gls[genotype_index(used[low], used[high])]);
    }
    // scaleLogSpaceArrayForNumericalStability, then GenotypeBuilder.PL(double[]).
    let maximum = subset[max_element_index(&subset)];
    for value in &mut subset {
        *value -= maximum;
    }
    (used, pls(&subset))
}

fn calculate_af(alleles: &[Vec<u8>], input_pls: &[i32]) -> AfResult {
    let allele_count = alleles.len();
    let gls: Vec<f64> = input_pls.iter().map(|&pl| pl as f64 / -10.0).collect();
    let reference_length = alleles[0].len();
    let prior_pseudocounts: Vec<f64> = alleles
        .iter()
        .enumerate()
        .map(|(index, allele)| {
            if index == 0 {
                REF_PSEUDOCOUNT
            } else if allele_length(allele) == reference_length {
                SNP_PSEUDOCOUNT
            } else {
                INDEL_PSEUDOCOUNT
            }
        })
        .collect();

    let flat_frequency = -jmath::log10(allele_count as f64);
    let mut log10_frequencies = vec![flat_frequency; allele_count];
    let mut allele_counts = vec![0.0; allele_count];

    loop {
        let new_counts = effective_allele_counts(&gls, &log10_frequencies);
        let mut maximum_difference = f64::NEG_INFINITY;
        for index in 0..allele_count {
            let difference = (allele_counts[index] - new_counts[index]).abs();
            if difference > maximum_difference {
                maximum_difference = difference;
            }
        }
        allele_counts = new_counts;

        let mut posterior_pseudocounts = Vec::with_capacity(allele_count);
        for index in 0..allele_count {
            posterior_pseudocounts.push(prior_pseudocounts[index] + allele_counts[index]);
        }
        let mut pseudocount_sum = 0.0;
        for &value in &posterior_pseudocounts {
            pseudocount_sum += value;
        }
        log10_frequencies.clear();
        for &value in &posterior_pseudocounts {
            log10_frequencies.push(jmath::log10(value / pseudocount_sum));
        }

        if maximum_difference <= ALLELE_COUNT_CONVERGENCE {
            break;
        }
    }

    let posteriors = normalized_genotype_posteriors(&gls, &log10_frequencies);
    let span_del_index = alleles.iter().position(|allele| allele == SPAN_DEL);
    let mut log10_p_no_variant = 0.0;
    if let Some(span_del) = span_del_index {
        let non_variant = [
            posteriors[genotype_index(0, 0)],
            posteriors[genotype_index(0, span_del)],
            posteriors[genotype_index(span_del, span_del)],
        ];
        log10_p_no_variant += java_min_zero(log10_sum_log10(&non_variant));
    } else {
        log10_p_no_variant += posteriors[0];
    }

    let mut log10_p_absent = vec![0.0; allele_count];
    if allele_count == 2 && span_del_index.is_none() {
        log10_p_absent[1] = log10_p_no_variant;
    } else {
        let mut absent_posteriors = vec![Vec::new(); allele_count];
        for (index, (low, high)) in diploid_genotypes(allele_count).enumerate() {
            for (allele, absent) in absent_posteriors.iter_mut().enumerate() {
                if allele != low && allele != high {
                    absent.push(posteriors[index]);
                }
            }
        }
        for allele in 0..allele_count {
            log10_p_absent[allele] = java_min_zero(log10_sum_log10(&absent_posteriors[allele]));
        }
    }

    let rounded_counts: Vec<i32> = allele_counts
        .iter()
        .map(|&count| java_round(count) as i32)
        .collect();
    AfResult {
        allele_counts: rounded_counts[1..].to_vec(),
        log10_p_no_variant,
        log10_p_absent,
    }
}

fn effective_allele_counts(gls: &[f64], log10_frequencies: &[f64]) -> Vec<f64> {
    let posteriors = normalized_genotype_posteriors(gls, log10_frequencies);
    let allele_count = log10_frequencies.len();
    let mut log10_result = vec![f64::NEG_INFINITY; allele_count];
    for (index, (low, high)) in diploid_genotypes(allele_count).enumerate() {
        if low == high {
            let term = posteriors[index] + jmath::log10(2.0);
            log10_result[low] = log10_sum_two(log10_result[low], term);
        } else {
            let low_term = posteriors[index] + jmath::log10(1.0);
            log10_result[low] = log10_sum_two(log10_result[low], low_term);
            let high_term = posteriors[index] + jmath::log10(1.0);
            log10_result[high] = log10_sum_two(log10_result[high], high_term);
        }
    }
    for value in &mut log10_result {
        *value = jmath::pow(10.0, *value);
    }
    log10_result
}

fn normalized_genotype_posteriors(gls: &[f64], log10_frequencies: &[f64]) -> Vec<f64> {
    let mut posteriors = Vec::with_capacity(gls.len());
    for (index, (low, high)) in diploid_genotypes(log10_frequencies.len()).enumerate() {
        let combination_count = if low == high {
            0.0
        } else {
            // commons-math 3.5 factorialLog(2) is FastMath.log(2); for this
            // input it has the same fdlibm value as the reference JVM.
            jmath::log(2.0) * log10_e()
        };
        let mut frequency_sum = 0.0;
        if low == high {
            frequency_sum += 2.0 * log10_frequencies[low];
        } else {
            frequency_sum += log10_frequencies[low];
            frequency_sum += log10_frequencies[high];
        }
        posteriors.push(combination_count + gls[index] + frequency_sum);
    }
    let normalization = log10_sum_log10(&posteriors);
    for posterior in &mut posteriors {
        *posterior -= normalization;
    }
    posteriors
}

fn subset_and_call(input_pls: &[i32], kept: &[usize]) -> ([usize; 2], Vec<i32>, i32) {
    let mut likelihoods = Vec::with_capacity(kept.len() * (kept.len() + 1) / 2);
    for high in 0..kept.len() {
        for low in 0..=high {
            likelihoods.push(input_pls[genotype_index(kept[low], kept[high])] as f64 / -10.0);
        }
    }

    let max_index = max_element_index(&likelihoods);
    let maximum = likelihoods[max_index];
    for likelihood in &mut likelihoods {
        *likelihood -= maximum;
    }

    let max_index = max_element_index(&likelihoods);
    let mut second_best = f64::NEG_INFINITY;
    for (index, &likelihood) in likelihoods.iter().enumerate() {
        if index != max_index && likelihood >= second_best {
            second_best = likelihood;
        }
    }
    let difference = likelihoods[max_index] - second_best;
    let log10_gq = -difference;
    let gq = java_round(log10_gq * -10.0) as i32;

    let mut chosen = [0, 0];
    for (index, (low, high)) in diploid_genotypes(kept.len()).enumerate() {
        if index == max_index {
            chosen = [low, high];
            break;
        }
    }
    if max_index == 0 && log10_gq > -0.1 {
        chosen = [NO_CALL_ALLELE, NO_CALL_ALLELE];
    }

    (chosen, pls(&likelihoods), gq)
}

#[inline]
fn genotype_index(a: usize, b: usize) -> usize {
    let (low, high) = if a <= b { (a, b) } else { (b, a) };
    high * (high + 1) / 2 + low
}

fn diploid_genotypes(allele_count: usize) -> impl Iterator<Item = (usize, usize)> {
    (0..allele_count).flat_map(|high| (0..=high).map(move |low| (low, high)))
}

#[inline]
fn log10_sum_two(a: f64, b: f64) -> f64 {
    if a > b {
        a + jmath::log10(1.0 + jmath::pow(10.0, b - a))
    } else {
        b + jmath::log10(1.0 + jmath::pow(10.0, a - b))
    }
}

#[inline]
fn java_min_zero(value: f64) -> f64 {
    java_min(0.0, value)
}

#[inline]
fn java_min(a: f64, b: f64) -> f64 {
    if a.is_nan() {
        a
    } else if b.is_nan() {
        b
    } else if a == 0.0 && b == 0.0 {
        if a.is_sign_negative() || b.is_sign_negative() {
            -0.0
        } else {
            0.0
        }
    } else if a <= b {
        a
    } else {
        b
    }
}

#[inline]
fn java_max(a: f64, b: f64) -> f64 {
    if a.is_nan() {
        a
    } else if b.is_nan() {
        b
    } else if a == 0.0 && b == 0.0 {
        if a.is_sign_positive() || b.is_sign_positive() {
            0.0
        } else {
            -0.0
        }
    } else if a >= b {
        a
    } else {
        b
    }
}

fn log10_one_minus_pow10(a: f64) -> f64 {
    if a > 0.0 {
        return f64::NAN;
    }
    if a == 0.0 {
        return f64::NEG_INFINITY;
    }
    let log_ten = jmath::log(10.0);
    let b = a * log_ten;
    let threshold = jmath::log(0.5);
    let natural = if b < threshold {
        jmath::log1p(-jmath::exp(b))
    } else {
        jmath::log(-jmath::expm1(b))
    };
    natural * (1.0 / log_ten)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn diploid_order_and_likelihoods() {
        let matrix = vec![vec![-1.0], vec![-2.0], vec![-3.0]];
        let gls = genotype_likelihoods(&matrix);
        assert_eq!(gls.len(), 6);
        assert_eq!(gls[0], -1.0);
        assert_eq!(gls[2], -2.0);
        assert_eq!(gls[5], -3.0);
    }

    #[test]
    fn genotype_index_is_vcf_order() {
        assert_eq!(genotype_index(0, 0), 0);
        assert_eq!(genotype_index(0, 1), 1);
        assert_eq!(genotype_index(1, 1), 2);
        assert_eq!(genotype_index(0, 2), 3);
        assert_eq!(genotype_index(1, 2), 4);
        assert_eq!(genotype_index(2, 2), 5);
    }
}
