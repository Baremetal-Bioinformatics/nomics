//! GATK 4.6.2.0 numeric tables and log-space helpers (`QualityUtils`,
//! `MathUtils`), bit-exact: every transcendental goes through `jmath`
//! (fdlibm, which the reference JVM's `Math` matched), and expressions keep
//! Java's evaluation order. Gate: `tests/gatk_tables.rs` against values
//! dumped from the reference jar (`reference/oracle/hc/GatkTables.java`).

use std::sync::OnceLock;

use crate::jmath;

/// `QualityUtils.MAX_QUAL`.
pub const MAX_QUAL: usize = 254;

/// `MathUtils.LOG10_ONE_THIRD = -Math.log10(3.0)`.
pub fn log10_one_third() -> f64 {
    -jmath::log10(3.0)
}

/// `MathUtils.LOG10_E = Math.log10(Math.E)`.
pub fn log10_e() -> f64 {
    jmath::log10(std::f64::consts::E)
}

struct QualTables {
    error_prob: [f64; MAX_QUAL + 1],
    prob_log10: [f64; MAX_QUAL + 1],
}

fn qual_tables() -> &'static QualTables {
    static T: OnceLock<QualTables> = OnceLock::new();
    T.get_or_init(|| {
        let mut t = QualTables {
            error_prob: [0.0; MAX_QUAL + 1],
            prob_log10: [0.0; MAX_QUAL + 1],
        };
        for i in 0..=MAX_QUAL {
            t.error_prob[i] = qual_to_error_prob_f(i as f64);
            t.prob_log10[i] = jmath::log10(1.0 - t.error_prob[i]);
        }
        t
    })
}

/// `QualityUtils.qualToErrorProb(double) = Math.pow(10.0, qual / -10.0)`.
pub fn qual_to_error_prob_f(qual: f64) -> f64 {
    jmath::pow(10.0, qual / -10.0)
}

/// `QualityUtils.qualToErrorProb(byte)` (cached; `qual & 0xff`).
pub fn qual_to_error_prob(qual: u8) -> f64 {
    qual_tables().error_prob[qual as usize]
}

/// `QualityUtils.qualToProb(byte) = 1.0 - qualToErrorProb(qual)`.
pub fn qual_to_prob(qual: u8) -> f64 {
    1.0 - qual_to_error_prob(qual)
}

/// `QualityUtils.qualToProbLog10(byte)` (cached).
pub fn qual_to_prob_log10(qual: u8) -> f64 {
    qual_tables().prob_log10[qual as usize]
}

/// `QualityUtils.qualToErrorProbLog10(byte) = (qual & 0xff) * -0.1`.
pub fn qual_to_error_prob_log10(qual: u8) -> f64 {
    qual as f64 * -0.1
}

/// `MathUtils.fastRound`.
#[inline]
pub fn fast_round(d: f64) -> i32 {
    if d > 0.0 {
        (d + 0.5) as i32
    } else {
        (d - 0.5) as i32
    }
}

const JACOBIAN_MAX_TOLERANCE: f64 = 8.0;
const JACOBIAN_TABLE_STEP: f64 = 0.0001;

fn jacobian_table() -> &'static [f64] {
    static T: OnceLock<Vec<f64>> = OnceLock::new();
    T.get_or_init(|| {
        // new IndexRange(0, (int) (MAX_TOLERANCE / TABLE_STEP) + 1)
        let n = (JACOBIAN_MAX_TOLERANCE / JACOBIAN_TABLE_STEP) as i32 + 1;
        (0..n)
            .map(|k| jmath::log10(1.0 + jmath::pow(10.0, -(k as f64) * JACOBIAN_TABLE_STEP)))
            .collect()
    })
}

/// `MathUtils.approximateLog10SumLog10(a, b)`: Jacobian-table approximation
/// of `log10(10^a + 10^b)`.
pub fn approximate_log10_sum_log10(a: f64, b: f64) -> f64 {
    let (a, b) = if a > b { (b, a) } else { (a, b) };
    if a == f64::NEG_INFINITY {
        return b;
    }
    let diff = b - a;
    b + if diff < JACOBIAN_MAX_TOLERANCE {
        // INV_STEP is the computed double 1.0 / TABLE_STEP, not 10000.
        let inv_step = 1.0 / JACOBIAN_TABLE_STEP;
        jacobian_table()[fast_round(diff * inv_step) as usize]
    } else {
        0.0
    }
}

/// `MathUtils.normalDistribution(mean, sd, x)`.
pub fn normal_distribution(mean: f64, sd: f64, x: f64) -> f64 {
    let root_two_pi = (2.0 * std::f64::consts::PI).sqrt();
    jmath::exp(-(x - mean) * (x - mean) / (2.0 * sd * sd)) / (sd * root_two_pi)
}

/// `MathUtils.maxElementIndex`: first index of the maximum.
pub fn max_element_index(a: &[f64]) -> usize {
    let mut m = 0;
    for i in 1..a.len() {
        if a[i] > a[m] {
            m = i;
        }
    }
    m
}

/// `MathUtils.log10SumLog10(values)`.
pub fn log10_sum_log10(v: &[f64]) -> f64 {
    if v.is_empty() {
        return f64::NEG_INFINITY;
    }
    let mi = max_element_index(v);
    let max = v[mi];
    if max == f64::NEG_INFINITY {
        return max;
    }
    let mut sum = 1.0;
    for (i, &x) in v.iter().enumerate() {
        if i == mi || x == f64::NEG_INFINITY {
            continue;
        }
        sum += jmath::pow(10.0, x - max);
    }
    max + if sum != 1.0 { jmath::log10(sum) } else { 0.0 }
}

/// `MathUtils.normalizeFromLog10ToLinearSpace` (in place).
pub fn normalize_from_log10_to_linear(v: &mut [f64]) {
    let s = log10_sum_log10(v);
    for x in v.iter_mut() {
        *x -= s;
    }
    for x in v.iter_mut() {
        *x = jmath::pow(10.0, *x);
    }
}

/// `Math.round(double)` (Java 7+: `floor(x + 1/2)` computed exactly).
pub fn java_round(x: f64) -> i64 {
    if x.is_nan() {
        return 0;
    }
    let f = x.floor();
    let r = if x - f >= 0.5 { f + 1.0 } else { f };
    r as i64
}

/// htsjdk `GenotypeLikelihoods.GLsToPLs` then `PL / -10.0`: the likelihoods
/// a `GenotypeBuilder.PL(double[])` genotype hands back.
pub fn pl_rounded(gls: &[f64]) -> Vec<f64> {
    let adjust = gls[max_element_index(gls)];
    gls.iter()
        .map(|&g| (java_round((-10.0 * (g - adjust)).min(i32::MAX as f64)) as i32) as f64 / -10.0)
        .collect()
}
