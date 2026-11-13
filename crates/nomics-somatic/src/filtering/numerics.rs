//! GATK NaturalLogUtils, ThresholdCalculator and Commons Math Beta CDF.
use crate::math::{fastmath, java_round, stream_sum};
use crate::orientation::{jmath, math};
pub(super) fn log_sum(x: &[f64]) -> f64 {
    let mut max = 0;
    for i in 1..x.len() {
        if x[i] > x[max] {
            max = i;
        }
    }
    if x.is_empty() || x[max] == f64::NEG_INFINITY {
        return f64::NEG_INFINITY;
    }
    let mut sum = 1.0;
    for (i, &v) in x.iter().enumerate() {
        if i != max && v != f64::NEG_INFINITY {
            sum += jmath::exp(v - x[max]);
        }
    }
    x[max] + if sum != 1.0 { jmath::log(sum) } else { 0.0 }
}
pub(super) fn normalize(x: &[f64]) -> Vec<f64> {
    let s = log_sum(x);
    x.iter().map(|&v| jmath::exp(v - s)).collect()
}
pub(super) fn log1mexp(x: f64) -> f64 {
    if x < -std::f64::consts::LN_2 {
        jmath::log1p(-jmath::exp(x))
    } else {
        jmath::log(-jmath::expm1(x))
    }
}
pub(super) fn posterior(odds: f64, prior: f64) -> f64 {
    normalize(&[odds + prior, log1mexp(prior)])[1]
}
pub(super) fn combined_error(p: [f64; 3], order: [usize; 3]) -> f64 {
    let mut truth = 1.0;
    for t in order {
        truth *= 1.0 - p[t];
    }
    (1.0 - truth).clamp(0.0, 1.0)
}
pub(super) fn threshold(mut p: Vec<f64>) -> f64 {
    p.sort_by(f64::total_cmp);
    let mut fn_ = stream_sum(p.iter().map(|p| 1.0 - p));
    let (mut tp, mut fp, mut best, mut index) = (0.0, 0.0, 0.0, None);
    for (i, &v) in p.iter().enumerate() {
        tp += 1.0 - v;
        fp += v;
        fn_ -= 1.0 - v;
        let f = 2.0 * tp / (2.0 * tp + fn_ + fp);
        if f >= best {
            best = f;
            index = Some(i);
        }
    }
    match index {
        None => 0.0,
        Some(i) if i == p.len() - 1 => 1.0,
        Some(i) => p[i],
    }
}
pub(super) fn round(x: f64, places: u32) -> String {
    let ulp = if x == 0.0 {
        f64::from_bits(1)
    } else {
        let v = x.abs();
        v.next_up() - v
    };
    let scale = jmath::pow(10.0, places as f64);
    crate::math::javafmt::double_to_string(java_round((x + ulp) * scale) as f64 / scale)
}
pub(super) fn qual(p: f64) -> i64 {
    java_round(-10.0 * jmath::log10(p)).clamp(1, 93)
}
pub(super) fn regularized_beta(x: f64, a: f64, b: f64) -> f64 {
    if x > (a + 1.0) / (2.0 + b + a) && 1.0 - x <= (b + 1.0) / (2.0 + b + a) {
        return 1.0 - regularized_beta(1.0 - x, b, a);
    }
    let (mut h, mut d, mut c) = (1.0, 0.0, 1.0);
    let mut converged = false;
    for n in 1..i32::MAX {
        let bn = if n % 2 == 0 {
            let m = n as f64 / 2.0;
            (m * (b - m) * x) / ((a + (2.0 * m) - 1.0) * (a + (2.0 * m)))
        } else {
            let m = (n as f64 - 1.0) / 2.0;
            -((a + m) * (a + b + m) * x) / ((a + (2.0 * m)) * (a + (2.0 * m) + 1.0))
        };
        d = 1.0 + bn * d;
        if d.abs() <= 1e-50 {
            d = 1e-50;
        }
        c = 1.0 + bn / c;
        if c.abs() <= 1e-50 {
            c = 1e-50;
        }
        d = 1.0 / d;
        let delta = c * d;
        h *= delta;
        assert!(h.is_finite(), "Beta continued fraction diverged");
        if (delta - 1.0).abs() < 1e-14 {
            converged = true;
            break;
        }
    }
    assert!(
        converged,
        "unsupported: Beta continued fraction iteration limit reached"
    );
    fastmath::exp((a * math::log(x)) + (b * math::log1p(-x)) - math::log(a) - math::log_beta(a, b))
        * 1.0
        / h
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn java_filtering_numerics() {
        let expected: Vec<_> = include_str!("../../tests/oracle/filtering_numerics.bits")
            .lines()
            .map(|s| u64::from_str_radix(s, 16).unwrap())
            .collect();
        let mut actual = Vec::new();
        for x in [
            0.01, 0.1, 0.5, 0.7, 1.0, 1.2, 1.7, 2.0, 2.5, 3.0, 7.9, 8.0, 10.0, 49.0, 1000.0,
        ] {
            actual.push(crate::math::gamma::log_gamma(x));
            actual.push(crate::math::gamma::digamma(x));
        }
        for (a, b) in [
            (0.01, 0.2),
            (0.3, 20.0),
            (1.79, 3.17),
            (10.02, 0.5),
            (99.0, 10001.0),
            (1.0, 50000.0),
            (0.5, 2.0),
            (2.0, 2.0),
        ] {
            actual.push(math::log_beta(a, b));
            for x in [0.001, 0.1, 0.5, 0.9] {
                actual.push(regularized_beta(x, a, b));
            }
        }
        assert_eq!(actual.len(), expected.len());
        for (i, (a, e)) in actual.into_iter().zip(expected).enumerate() {
            assert_eq!(a.to_bits(), e, "filtering oracle {i}: {a}");
        }
    }
    #[test]
    fn threshold_boundaries() {
        assert_eq!(threshold(vec![]), 0.0);
        assert_eq!(threshold(vec![0.0, 0.0]), 1.0);
        assert_eq!(threshold(vec![0.0, 1.0]), 0.0);
        assert_eq!(qual(0.0), 93);
        assert_eq!(qual(1.0), 1);
    }
    #[test]
    fn error_type_order_can_reach_written_threshold_digits() {
        // Arithmetic-boundary witness, not a claim that a particular VCF
        // realizes these three filter probabilities. Checked independently in
        // FilteringLearning.java, including Java's ulp-adjusted rounding.
        let p = [
            0.04163411606541967,
            0.006318578253251994,
            0.03865375732824039,
        ];
        let a = combined_error(p, [0, 1, 2]);
        let b = combined_error(p, [0, 2, 1]);
        assert_eq!(a.to_bits(), 0x3fb5_a1ca_c083_1270);
        assert_eq!(b.to_bits(), 0x3fb5_a1ca_c083_1268);
        assert_eq!(round(threshold(vec![0.0, a, 1.0]), 3), "0.085");
        assert_eq!(round(threshold(vec![0.0, b, 1.0]), 3), "0.084");
    }
}
