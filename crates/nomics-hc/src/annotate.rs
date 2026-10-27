//! GATK 4.6.2.0 HC default annotations, ported from `VariantAnnotatorEngine`,
//! `MannWhitneyU`, and commons-math3's Loader saddle-point code.  The gate is
//! `tests/annotate.rs` against `GatkCall.java`'s tier0 annotation matrices.

use std::collections::BTreeMap;

use crate::javafmt::format_fixed;
use crate::jmath;
use crate::jrandom::JavaRandom;
use crate::tables::{java_round, log10_e, log10_sum_log10};
use crate::vcf::InfoValue;

const INFORMATIVE_THRESHOLD: f64 = 0.2;

/// The called genotype and precomputed maximum-likelihood allele fields.
#[derive(Clone, Debug)]
pub struct AnnotationGenotype {
    pub alleles: [usize; 2],
    pub pl: Vec<i32>,
    pub gq: i32,
    pub mle_ac: Vec<i64>,
    pub mle_af: Vec<f64>,
}

/// One row of the realigned read-by-allele annotation matrix.
#[derive(Clone, Debug)]
pub struct AnnotationRead {
    pub flag: u16,
    pub mapq: u8,
    pub reverse: bool,
    /// One-based aligned start.
    pub start: i64,
    /// BAM-packed CIGAR elements.
    pub cigar: Vec<u32>,
    pub bases: Vec<u8>,
    pub quals: Vec<u8>,
    /// Log10 likelihoods in `AnnotationCall::alleles` order.
    pub likelihoods: Vec<f64>,
}

/// Inputs to one `VariantAnnotatorEngine.annotateContext` call.
#[derive(Clone, Debug)]
pub struct AnnotationCall {
    /// REF and emitted ALTs, in VC order.
    pub alleles: Vec<Vec<u8>>,
    /// Alleles in the annotation matrix; this can retain dropped ALTs.
    pub matrix_alleles: Vec<Vec<u8>>,
    /// GATK's exact log10 probability of error (negative; QUAL is `-10*x`).
    pub log10_p_error: f64,
    pub genotype: AnnotationGenotype,
    pub contig: String,
    pub pos: i64,
    pub end: i64,
    pub reads: Vec<AnnotationRead>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct AnnotationResult {
    /// Alphabetical, as htsjdk's INFO `TreeMap` writes it.
    pub info: Vec<(String, InfoValue)>,
    pub ad: Vec<i32>,
    pub dp: i32,
    pub sb: Option<[i32; 4]>,
}

#[derive(Clone, Copy)]
struct BestAllele {
    allele: usize,
    likelihood: f64,
    second_likelihood: f64,
}

impl BestAllele {
    fn informative(self) -> bool {
        let confidence = if self.likelihood == self.second_likelihood {
            0.0
        } else {
            self.likelihood - self.second_likelihood
        };
        confidence > INFORMATIVE_THRESHOLD
    }
}

/// `AlleleLikelihoods.searchBestAllele(..., reference-priority)`.
fn best_allele(values: &[f64]) -> BestAllele {
    assert!(!values.is_empty());
    let mut best = 0;
    let mut second = 0;
    let mut best_lk = values[0];
    let mut second_lk = f64::NEG_INFINITY;
    for (a, &candidate) in values.iter().enumerate().skip(1) {
        if candidate > best_lk {
            second = best;
            best = a;
            second_lk = best_lk;
            best_lk = candidate;
        } else if candidate > second_lk {
            second = a;
            second_lk = candidate;
        }
    }
    if best_lk - second_lk < INFORMATIVE_THRESHOLD {
        let mut best_priority = usize::from(best == 0);
        let mut second_priority = usize::from(second == 0);
        for (a, &candidate) in values.iter().enumerate() {
            if a == best || best_lk - candidate > INFORMATIVE_THRESHOLD {
                continue;
            }
            let priority = usize::from(a == 0);
            if priority > best_priority {
                second = best;
                best = a;
                second_priority = best_priority;
                best_priority = priority;
            } else if priority > second_priority {
                second = a;
                second_priority = priority;
            }
        }
    }
    BestAllele {
        allele: best,
        likelihood: values[best],
        second_likelihood: if second == best {
            f64::NEG_INFINITY
        } else {
            values[second]
        },
    }
}

#[inline]
fn cigar_op(c: u32) -> u32 {
    c & 0xf
}

#[inline]
fn cigar_len(c: u32) -> i64 {
    i64::from(c >> 4)
}

#[inline]
fn consumes_read(op: u32) -> bool {
    matches!(op, 0 | 1 | 4 | 7 | 8)
}

#[inline]
fn consumes_ref(op: u32) -> bool {
    matches!(op, 0 | 2 | 3 | 7 | 8)
}

fn aligned_end(read: &AnnotationRead) -> i64 {
    let n: i64 = read
        .cigar
        .iter()
        .filter(|&&c| consumes_ref(cigar_op(c)))
        .map(|&c| cigar_len(c))
        .sum();
    if n == 0 {
        read.start
    } else {
        read.start + n - 1
    }
}

fn soft_start(read: &AnnotationRead) -> i64 {
    let mut start = read.start;
    for &c in &read.cigar {
        match cigar_op(c) {
            4 => start -= cigar_len(c),
            5 => {}
            _ => break,
        }
    }
    start
}

fn soft_end(read: &AnnotationRead) -> i64 {
    let mut end = aligned_end(read);
    let mut found_aligned = false;
    for &c in read.cigar.iter().rev() {
        match cigar_op(c) {
            4 => end += cigar_len(c),
            5 => {}
            _ => {
                found_aligned = true;
                break;
            }
        }
    }
    if found_aligned {
        end
    } else {
        aligned_end(read)
    }
}

/// Returns `(read index, cigar operator)` with soft clips consuming reference,
/// exactly as `ReadUtils.getReadIndexForReferenceCoordinate`.
fn read_index_at(read: &AnnotationRead, coordinate: i64) -> Option<(usize, u32)> {
    let soft_start = soft_start(read);
    let mut last_read = 0_i64;
    let mut last_ref = soft_start;
    if coordinate < soft_start {
        return None;
    }
    for &c in &read.cigar {
        let op = cigar_op(c);
        let first_read = last_read;
        let first_ref = last_ref;
        if consumes_read(op) {
            last_read += cigar_len(c);
        }
        if consumes_ref(op) || op == 4 {
            last_ref += cigar_len(c);
        }
        if first_ref <= coordinate && coordinate < last_ref {
            let at = first_read
                + if consumes_read(op) {
                    coordinate - first_ref
                } else {
                    0
                };
            return usize::try_from(at).ok().map(|i| (i, op));
        }
    }
    None
}

fn base_quality(read: &AnnotationRead, pos: i64) -> Option<f64> {
    if pos < read.start || pos > aligned_end(read) {
        return None;
    }
    let (i, op) = read_index_at(read, pos)?;
    consumes_read(op).then(|| f64::from(read.quals[i]))
}

fn read_position(read: &AnnotationRead, pos: i64, end: i64) -> Option<f64> {
    let first_core = read
        .cigar
        .iter()
        .map(|&c| cigar_op(c))
        .find(|&op| !matches!(op, 4 | 5));
    if read.start == end + 1 && first_core == Some(1) {
        return Some(0.0);
    }
    if soft_start(read) > end + 1 || soft_end(read) < pos {
        return None;
    }
    let (i, _) = read_index_at(read, pos)?;
    let leading_hard = read
        .cigar
        .first()
        .filter(|&&c| cigar_op(c) == 5)
        .map_or(0, |&c| cigar_len(c) as usize);
    let trailing_hard = read
        .cigar
        .last()
        .filter(|&&c| cigar_op(c) == 5)
        .map_or(0, |&c| cigar_len(c) as usize);
    let left = leading_hard + i;
    let right = read.bases.len() - 1 - i + trailing_hard;
    Some(left.min(right) as f64)
}

fn rank_sum(
    reads: &[AnnotationRead],
    best: &[BestAllele],
    matrix_to_call: &[Option<usize>],
    mut value: impl FnMut(&AnnotationRead) -> Option<f64>,
) -> Option<String> {
    let mut alt = Vec::new();
    let mut reference = Vec::new();
    for (read, &ba) in reads.iter().zip(best) {
        if !ba.informative() || read.mapq == 0 || read.mapq == 255 {
            continue;
        }
        if let Some(v) = value(read) {
            match matrix_to_call[ba.allele] {
                Some(0) => reference.push(v),
                Some(_) => alt.push(v),
                None => {}
            }
        }
    }
    mann_whitney_z(&mut alt, &mut reference).map(|z| format_fixed(z, 3))
}

fn ranked(series1: &mut [f64], series2: &mut [f64]) -> (Vec<(f32, u8)>, f64) {
    series1.sort_by(f64::total_cmp);
    series2.sort_by(f64::total_cmp);
    let mut ranks = Vec::with_capacity(series1.len() + series2.len());
    let (mut i, mut j) = (0, 0);
    while i < series1.len() || j < series2.len() {
        let which = if i >= series1.len() {
            2
        } else if j >= series2.len() || series1[i] <= series2[j] {
            1
        } else {
            2
        };
        let value = if which == 1 {
            let x = series1[i];
            i += 1;
            x
        } else {
            let x = series2[j];
            j += 1;
            x
        };
        ranks.push((ranks.len() as f32 + 1.0, which, value));
    }
    let mut ties = 0.0;
    let mut at = 0;
    while at < ranks.len() {
        let mut end = at + 1;
        let mut rank = ranks[at].0;
        while end < ranks.len() && ranks[end].2 == ranks[at].2 {
            rank += ranks[end].0;
            end += 1;
        }
        let count = end - at;
        if count > 1 {
            rank /= count as f32;
            for item in &mut ranks[at..end] {
                item.0 = rank;
            }
            if count != ranks.len() {
                ties += jmath::pow(count as f64, 3.0) - count as f64;
            }
        }
        at = end;
    }
    (ranks.into_iter().map(|(r, s, _)| (r, s)).collect(), ties)
}

fn mann_whitney_z(series1: &mut [f64], series2: &mut [f64]) -> Option<f64> {
    let (n1, n2) = (series1.len(), series2.len());
    if n1 == 0 || n2 == 0 {
        return None;
    }
    let (ranks, ties) = ranked(series1, series2);
    let mut r1 = 0.0_f32;
    for &(rank, series) in &ranks {
        if series == 1 {
            r1 += rank;
        }
    }
    let u = f64::from(r1) - (n1 * (n1 + 1)) as f64 / 2.0;
    if n1 >= 10 || n2 >= 10 {
        let m = (n1 * n2) as f64 / 2.0;
        let correction = if ties == 0.0 { 0.0 } else { -0.5 };
        let divisor = (n1 + n2) * (n1 + n2 - 1);
        let sigma =
            (((n1 * n2) as f64 / 12.0) * ((n1 + n2 + 1) as f64 - ties / divisor as f64)).sqrt();
        Some((u - m - correction) / sigma)
    } else {
        let target = java_round(2.0 * u);
        let rank_values: Vec<f32> = ranks.iter().map(|x| x.0).collect();
        let mut observed = 0_u64;
        let mut below = 0_u64;
        let mut total = 0_u64;
        enumerate_rank_subsets(
            &rank_values,
            n1,
            0,
            0,
            0.0,
            target,
            &mut observed,
            &mut below,
            &mut total,
        );
        let p = (observed as f64 / 2.0 + below as f64) / total as f64;
        Some(normal_inverse(p))
    }
}

#[allow(clippy::too_many_arguments)]
fn enumerate_rank_subsets(
    ranks: &[f32],
    choose: usize,
    at: usize,
    chosen: usize,
    sum: f32,
    target: i64,
    observed: &mut u64,
    below: &mut u64,
    total: &mut u64,
) {
    if chosen == choose {
        let u = f64::from(sum) - (choose * (choose + 1)) as f64 / 2.0;
        let bin = java_round(2.0 * u);
        *total += 1;
        if bin == target {
            *observed += 1;
        }
        if bin < target {
            *below += 1;
        }
        return;
    }
    if ranks.len() - at < choose - chosen {
        return;
    }
    for i in at..=ranks.len() - (choose - chosen) {
        enumerate_rank_subsets(
            ranks,
            choose,
            i + 1,
            chosen + 1,
            sum + ranks[i],
            target,
            observed,
            below,
            total,
        );
    }
}

/// commons-math3 3.5 `sqrt(2) * Erf.erfInv(2p-1)`.
fn normal_inverse(p0: f64) -> f64 {
    let x = 2.0 * p0 - 1.0;
    let mut w = -jmath::log((1.0 - x) * (1.0 + x));
    let mut p;
    if w < 6.25 {
        w -= 3.125;
        p = -3.6444120640178196996e-21;
        for c in [
            -1.685059138182016589e-19,
            1.2858480715256400167e-18,
            1.115787767802518096e-17,
            -1.333171662854620906e-16,
            2.0972767875968561637e-17,
            6.6376381343583238325e-15,
            -4.0545662729752068639e-14,
            -8.1519341976054721522e-14,
            2.6335093153082322977e-12,
            -1.2975133253453532498e-11,
            -5.4154120542946279317e-11,
            1.051212273321532285e-09,
            -4.1126339803469836976e-09,
            -2.9070369957882005086e-08,
            4.2347877827932403518e-07,
            -1.3654692000834678645e-06,
            -1.3882523362786468719e-05,
            0.0001867342080340571352,
            -0.00074070253416626697512,
            -0.0060336708714301490533,
            0.24015818242558961693,
            1.6536545626831027356,
        ] {
            p = c + p * w;
        }
    } else if w < 16.0 {
        w = w.sqrt() - 3.25;
        p = 2.2137376921775787049e-09;
        for c in [
            9.0756561938885390979e-08,
            -2.7517406297064545428e-07,
            1.8239629214389227755e-08,
            1.5027403968909827627e-06,
            -4.013867526981545969e-06,
            2.9234449089955446044e-06,
            1.2475304481671778723e-05,
            -4.7318229009055733981e-05,
            6.8284851459573175448e-05,
            2.4031110387097893999e-05,
            -0.0003550375203628474796,
            0.00095328937973738049703,
            -0.0016882755560235047313,
            0.0024914420961078508066,
            -0.0037512085075692412107,
            0.005370914553590063617,
            1.0052589676941592334,
            3.0838856104922207635,
        ] {
            p = c + p * w;
        }
    } else if w.is_finite() {
        w = w.sqrt() - 5.0;
        p = -2.7109920616438573243e-11;
        for c in [
            -2.5556418169965252055e-10,
            1.5076572693500548083e-09,
            -3.7894654401267369937e-09,
            7.6157012080783393804e-09,
            -1.4960026627149240478e-08,
            2.9147953450901080826e-08,
            -6.7711997758452339498e-08,
            2.2900482228026654717e-07,
            -9.9298272942317002539e-07,
            4.5260625972231537039e-06,
            -1.9681778105531670567e-05,
            7.5995277030017761139e-05,
            -0.00021503011930044477347,
            -0.00013871931833623122026,
            1.0103004648645343977,
            4.8499064014085844221,
        ] {
            p = c + p * w;
        }
    } else {
        p = f64::INFINITY;
    }
    2.0_f64.sqrt() * (p * x)
}

const STIRLING: [f64; 31] = [
    0.0,
    0.1534264097200273452913848,
    0.0810614667953272582196702,
    0.0548141210519176538961390,
    0.0413406959554092940938221,
    0.03316287351993628748511048,
    0.02767792568499833914878929,
    0.02374616365629749597132920,
    0.02079067210376509311152277,
    0.01848845053267318523077934,
    0.01664469118982119216319487,
    0.01513497322191737887351255,
    0.01387612882307074799874573,
    0.01281046524292022692424986,
    0.01189670994589177009505572,
    0.01110455975820691732662991,
    0.010411265261972096497478567,
    0.009799416126158803298389475,
    0.009255462182712732917728637,
    0.008768700134139385462952823,
    0.008330563433362871256469318,
    0.007934114564314020547248100,
    0.007573675487951840794972024,
    0.007244554301320383179543912,
    0.006942840107209529865664152,
    0.006665247032707682442354394,
    0.006408994188004207068439631,
    0.006171712263039457647532867,
    0.005951370112758847735624416,
    0.005746216513010115682023589,
    0.005554733551962801371038690,
];

fn stirling_error(z: f64) -> f64 {
    if z < 15.0 {
        STIRLING[(2.0 * z) as usize]
    } else {
        let z2 = z * z;
        (0.083333333333333333333
            - (0.00277777777777777777778
                - (0.00079365079365079365079365
                    - (0.000595238095238095238095238 - 0.0008417508417508417508417508 / z2) / z2)
                    / z2)
                / z2)
            / z
    }
}

fn deviance_part(x: f64, mu: f64) -> f64 {
    if (x - mu).abs() < 0.1 * (x + mu) {
        let d = x - mu;
        let mut v = d / (x + mu);
        let mut s1 = v * d;
        let mut s = f64::NAN;
        let mut ej = 2.0 * x * v;
        v *= v;
        let mut j = 1;
        while s1 != s {
            s = s1;
            ej *= v;
            s1 = s + ej / (j * 2 + 1) as f64;
            j += 1;
        }
        s1
    } else {
        x * jmath::log(x / mu) + mu - x
    }
}

fn log_binomial_probability(x: i32, n: i32, p: f64, q: f64) -> f64 {
    if x == 0 {
        if p < 0.1 {
            -deviance_part(n as f64, n as f64 * q) - n as f64 * p
        } else {
            n as f64 * jmath::log(q)
        }
    } else if x == n {
        if q < 0.1 {
            -deviance_part(n as f64, n as f64 * p) - n as f64 * q
        } else {
            n as f64 * jmath::log(p)
        }
    } else {
        let mut r = stirling_error(n as f64)
            - stirling_error(x as f64)
            - stirling_error((n - x) as f64)
            - deviance_part(x as f64, n as f64 * p)
            - deviance_part((n - x) as f64, n as f64 * q);
        let f = (2.0 * std::f64::consts::PI * x as f64 * (n - x) as f64) / n as f64;
        r += -0.5 * jmath::log(f);
        r
    }
}

fn hypergeometric_log_probability(x: i32, population: i32, successes: i32, sample: i32) -> f64 {
    let p = sample as f64 / population as f64;
    let q = (population - sample) as f64 / population as f64;
    log_binomial_probability(x, successes, p, q)
        + log_binomial_probability(sample - x, population - successes, p, q)
        - log_binomial_probability(sample, population, p, q)
}

fn fisher_pvalue(mut table: [[i32; 2]; 2]) -> f64 {
    let total: i32 = table.iter().flatten().sum();
    if total > 400 {
        let norm = total as f64 / 200.0;
        for x in table.iter_mut().flatten() {
            *x = (*x as f64 / norm) as i32;
        }
    }
    let m = table[0][0] + table[0][1];
    let n = table[1][0] + table[1][1];
    let k = table[0][0] + table[1][0];
    let lo = 0.max(k - n);
    let hi = k.min(m);
    if hi - lo + 1 <= 1 {
        return 1.0;
    }
    let logs: Vec<f64> = (lo..=hi)
        .map(|x| hypergeometric_log_probability(x, m + n, m, k))
        .collect();
    let threshold = logs[(table[0][0] - lo) as usize] * (1.0 - 1e-6);
    let selected: Vec<f64> = logs
        .into_iter()
        .filter(|&x| x <= threshold)
        .map(|x| x * log10_e())
        .collect();
    jmath::pow(10.0, log10_sum_log10(&selected)).min(1.0)
}

/// `FisherStrand`: phred-scaled two-sided Fisher p-value of a ref/alt by
/// forward/reverse table (normalized to 200 when larger than 400).
pub fn fisher_strand(table: [[i32; 2]; 2]) -> f64 {
    let p = fisher_pvalue(table).max(1e-320);
    (-10.0 * jmath::log10(p.max(f64::MIN_POSITIVE))).abs()
}

/// `StrandOddsRatio.calculateSOR` with the +1 pseudocount.
pub fn strand_odds_ratio(table: [[i32; 2]; 2]) -> f64 {
    let t00 = f64::from(table[0][0]) + 1.0;
    let t01 = f64::from(table[0][1]) + 1.0;
    let t10 = f64::from(table[1][0]) + 1.0;
    let t11 = f64::from(table[1][1]) + 1.0;
    let ratio = (t00 / t01) * (t11 / t10) + (t01 / t00) * (t10 / t11);
    let ref_ratio = t00.min(t01) / t00.max(t01);
    let alt_ratio = t10.min(t11) / t10.max(t11);
    jmath::log(ratio) + jmath::log(ref_ratio) - jmath::log(alt_ratio)
}

fn contingency(
    reads: &[AnnotationRead],
    best: &[BestAllele],
    matrix_to_call: &[Option<usize>],
    min_count: i32,
) -> [[i32; 2]; 2] {
    let mut flat = [0_i32; 4];
    for (read, &ba) in reads.iter().zip(best) {
        if ba.informative() {
            if let Some(call_allele) = matrix_to_call[ba.allele] {
                let allele = usize::from(call_allele != 0);
                let strand = usize::from(read.reverse);
                flat[allele * 2 + strand] += 1;
            }
        }
    }
    if flat.iter().sum::<i32>() > min_count {
        [[flat[0], flat[1]], [flat[2], flat[3]]]
    } else {
        [[0, 0], [0, 0]]
    }
}

/// Compute the default single-sample VCF-mode annotations.
pub fn annotate(call: &AnnotationCall, rng: &mut JavaRandom) -> AnnotationResult {
    annotate_mode(call, rng, false)
}

pub fn annotate_mode(call: &AnnotationCall, rng: &mut JavaRandom, gvcf: bool) -> AnnotationResult {
    assert!(!call.alleles.is_empty());
    assert!(
        call.reads
            .iter()
            .all(|r| r.likelihoods.len() == call.matrix_alleles.len())
    );
    let matrix_to_call: Vec<Option<usize>> = call
        .matrix_alleles
        .iter()
        .map(|a| call.alleles.iter().position(|b| a == b))
        .collect();
    let selected: Vec<usize> = call
        .alleles
        .iter()
        .map(|a| {
            call.matrix_alleles
                .iter()
                .position(|b| a == b)
                .expect("VC allele absent from annotation matrix")
        })
        .collect();
    let best: Vec<BestAllele> = call
        .reads
        .iter()
        .map(|r| best_allele(&r.likelihoods))
        .collect();
    let subset_best: Vec<BestAllele> = call
        .reads
        .iter()
        .map(|r| {
            let values: Vec<f64> = selected.iter().map(|&a| r.likelihoods[a]).collect();
            best_allele(&values)
        })
        .collect();
    let mut ad = vec![0_i32; call.alleles.len()];
    for ba in &subset_best {
        if ba.informative() {
            ad[ba.allele] += 1;
        }
    }
    let mut dp = ad.iter().sum();
    let mut info = BTreeMap::new();

    let mut ac = vec![0_i64; call.alleles.len() - 1];
    for &a in &call.genotype.alleles {
        if a != 0 && a != crate::genotype::NO_CALL_ALLELE {
            ac[a - 1] += 1;
        }
    }
    info.insert("AC".to_string(), InfoValue::Ints(ac.clone()));
    info.insert(
        "AF".to_string(),
        InfoValue::Doubles(ac.iter().map(|&x| x as f64 / 2.0).collect()),
    );
    info.insert("AN".to_string(), InfoValue::Ints(vec![2]));

    if let Some(x) = rank_sum(&call.reads, &best, &matrix_to_call, |r| {
        base_quality(r, call.pos)
    }) {
        info.insert("BaseQRankSum".to_string(), InfoValue::Str(x));
    }
    info.insert(
        "DP".to_string(),
        InfoValue::Str(call.reads.len().to_string()),
    );
    info.insert(
        "ExcessHet".to_string(),
        InfoValue::Str("0.0000".to_string()),
    );

    let fs_table = contingency(&call.reads, &best, &matrix_to_call, 2);
    info.insert(
        "FS".to_string(),
        InfoValue::Str(format_fixed(fisher_strand(fs_table), 3)),
    );
    info.insert(
        "MLEAC".to_string(),
        InfoValue::Ints(call.genotype.mle_ac.clone()),
    );
    info.insert(
        "MLEAF".to_string(),
        InfoValue::Doubles(call.genotype.mle_af.clone()),
    );

    let mut square_sum = 0_i64;
    let mut mq_n = 0_i64;
    for read in &call.reads {
        if read.mapq != 255 {
            square_sum += i64::from(read.mapq) * i64::from(read.mapq);
            mq_n += 1;
        }
    }
    let mq = (square_sum as f64 / mq_n as f64).sqrt();
    info.insert("MQ".to_string(), InfoValue::Str(format_fixed(mq, 2)));
    if let Some(x) = rank_sum(&call.reads, &best, &matrix_to_call, |r| {
        Some(f64::from(r.mapq))
    }) {
        info.insert("MQRankSum".to_string(), InfoValue::Str(x));
    }

    // QualByDepth.getDepth: het/hom-var genotypes only; their AD sum, or,
    // when AD sums to zero, the sample's evidence count (INFO DP).
    let qd_depth = if !call.genotype.alleles.iter().any(|&a| a != 0) {
        0
    } else if dp != 0 {
        dp
    } else {
        call.reads.len() as i32
    };
    if !gvcf && qd_depth != 0 {
        let mut qd = (-10.0 * call.log10_p_error) / qd_depth as f64;
        if qd >= 35.0 {
            qd = 30.0 + rng.next_gaussian() * 3.0;
        }
        info.insert("QD".to_string(), InfoValue::Str(format_fixed(qd, 2)));
    }
    if let Some(x) = rank_sum(&call.reads, &best, &matrix_to_call, |r| {
        read_position(r, call.pos, call.end)
    }) {
        info.insert("ReadPosRankSum".to_string(), InfoValue::Str(x));
    }

    let table = contingency(&call.reads, &best, &matrix_to_call, 0);
    info.insert(
        "SOR".to_string(),
        InfoValue::Str(format_fixed(strand_odds_ratio(table), 3)),
    );

    if gvcf {
        for key in ["AC", "AF", "AN", "FS", "SOR", "MQ"] {
            info.remove(key);
        }
        info.insert(
            "RAW_MQandDP".into(),
            InfoValue::Ints(vec![square_sum, mq_n]),
        );
        // Coverage/RMSMappingQuality emit nothing without evidence. The
        // genotype annotations and ExcessHet require a called genotype.
        if call.reads.is_empty() {
            info.remove("DP");
            info.remove("RAW_MQandDP");
        }
        if call
            .genotype
            .alleles
            .contains(&crate::genotype::NO_CALL_ALLELE)
        {
            info.remove("ExcessHet");
            ad.clear();
            dp = -1;
        }
    }
    AnnotationResult {
        sb: (gvcf && dp >= 0).then_some([table[0][0], table[0][1], table[1][0], table[1][1]]),
        info: info.into_iter().collect(),
        ad,
        dp,
    }
}
