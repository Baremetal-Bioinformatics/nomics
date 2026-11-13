//! Numerical boundary ports. See individual modules for upstream file names.
pub(crate) mod binomial;
pub(crate) mod gamma;
pub(crate) mod brent;
pub(crate) mod fastmath;
mod fastmath_tables;
pub(crate) mod java_log;
pub(crate) mod javafmt;
pub(crate) mod svd;

/// java.util.stream.DoublePipeline.sum / Collectors.sumWithCompensation.
/// Keep the separate simple sum for infinity recovery and the final *minus* correction.
pub(crate) fn stream_sum(values: impl IntoIterator<Item = f64>) -> f64 {
    let (mut sum, mut correction, mut simple) = (0.0, 0.0, 0.0);
    for x in values {
        simple += x;
        let corrected = x - correction;
        let next = sum + corrected;
        correction = (next - sum) - corrected;
        sum = next;
    }
    let result = sum - correction;
    if result.is_nan() && simple.is_infinite() {
        simple
    } else {
        result
    }
}
/// Commons Math stat/descriptive/moment/Mean.evaluate, two-pass correction.
pub(crate) fn mean(x: &[f64]) -> f64 {
    let initial = x.iter().copied().sum::<f64>() / x.len() as f64;
    initial + x.iter().map(|v| v - initial).sum::<f64>() / x.len() as f64
}
/// Commons Math stat/descriptive/rank/Percentile, default LEGACY estimator.
pub(crate) fn percentile(x: &[f64], p: f64) -> f64 {
    let mut x: Vec<_> = x.iter().copied().filter(|v| !v.is_nan()).collect();
    if x.is_empty() {
        return f64::NAN;
    }
    x.sort_by(f64::total_cmp);
    let pos = (p / 100.0) * (x.len() as f64 + 1.0);
    if pos < 1.0 {
        return x[0];
    }
    if pos >= x.len() as f64 {
        return x[x.len() - 1];
    }
    let floor = pos.floor();
    let i = floor as usize;
    x[i - 1] + (pos - floor) * (x[i] - x[i - 1])
}
/// GATK utils/MathUtils.binarySearchFindZero, including endpoint-sign tests.
pub(crate) fn binary_search_zero(
    f: impl Fn(f64) -> f64,
    mut bottom: f64,
    mut top: f64,
) -> Option<f64> {
    fn sign(x: f64) -> f64 {
        if x == 0.0 { x } else { x.signum() }
    }
    while top - bottom > 1e-6 {
        let mid = (bottom + top) / 2.0;
        let bv = f(bottom);
        let tv = f(top);
        let mv = f(mid);
        if sign(bv) == sign(tv) {
            return None;
        }
        if sign(bv) == sign(mv) {
            bottom = mid;
        } else {
            top = mid;
        }
    }
    Some((bottom + top) / 2.0)
}

/// OpenJDK java/lang/Math.round(double), avoiding a rounded x + 0.5.
pub(crate) fn java_round(x: f64) -> i64 {
    let bits = x.to_bits();
    let exponent = ((bits & 0x7ff0000000000000) >> 52) as i64;
    let shift = 1074 - exponent;
    if (0..64).contains(&shift) {
        let mut r = ((bits & 0x000fffffffffffff) | 0x0010000000000000) as i64;
        if x.is_sign_negative() {
            r = -r;
        }
        ((r >> shift) + 1) >> 1
    } else {
        x as i64
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn numerical_java_oracles() {
        // Independently compiled unmodified Commons Math MATH_3_5 Java sources
        // with OpenJDK 17. Numerics.java inputs are repeated here; fixture stores
        // Double.doubleToRawLongBits in call order (not Rust-derived values).
        let expected: Vec<_> = include_str!("java_oracle_bits.txt")
            .lines()
            .map(|s| u64::from_str_radix(s, 16).unwrap())
            .collect();
        let mut actual = Vec::new();
        for x in [
            -745., -710., -709., -100., -4., -0.1, 0., 0.01, 1., 100., 709.,
        ] {
            actual.push(fastmath::exp(x));
        }
        for x in [
            f64::from_bits(1),
            1e-300,
            0.001,
            0.1,
            0.995,
            1.,
            1.005,
            2.,
            100.,
            1e300,
        ] {
            actual.push(fastmath::log(x));
        }
        for (x, y) in [
            (3., 4.),
            (1e-300, 1e-300),
            (1e300, 1e300),
            (0., 0.),
            (1., 0.123),
        ] {
            actual.push(fastmath::hypot(x, y));
        }
        for (n, k) in [(0, 0), (10, 0), (10, 10), (10, 3), (150, 60), (300, 4)] {
            for p in [0., 0.01, 0.1, 0.4, 1.] {
                actual.push(binomial::binomial_probability(n, k, p));
            }
        }
        for guess in [0.02, 0.05, 0.1, 0.2, 0.4] {
            let r = brent::brent_optimizer(
                |x| -(x - 0.137) * (x - 0.137),
                0.,
                0.5,
                guess,
                0.0001,
                0.0001,
                30,
            )
            .unwrap();
            actual.extend([r.0, r.1]);
        }
        for matrix in [
            vec![vec![1., 2.], vec![2., 4.]],
            vec![vec![1., 0.2, 0.3], vec![0.2, 1., 0.7], vec![0.3, 0.7, 1.]],
        ] {
            let (s, u) = svd::singular_value_decomposition(matrix);
            actual.extend(s);
            actual.extend(u.into_iter().flatten());
        }
        for z in [0., 1., 14., 15., 20., 100.] {
            actual.push(binomial::stirling_error(z));
        }
        for (x, mu) in [(10., 9.9), (10., 4.), (1., 1.), (100., 110.)] {
            actual.push(binomial::deviance_part(x, mu));
        }
        assert_eq!(expected.len(), actual.len());
        for (i, (a, e)) in actual.into_iter().zip(expected).enumerate() {
            assert_eq!(a.to_bits(), e, "oracle index {i}: {a}");
        }
    }
    #[test]
    fn additional_java_oracles() {
        let expected: Vec<_> = include_str!("java_more_bits.txt")
            .lines()
            .map(|s| u64::from_str_radix(s, 16).unwrap())
            .collect();
        let mut actual = Vec::new();
        for x in [
            f64::from_bits(1),
            1e-300,
            0.001,
            0.1,
            0.995,
            1.,
            1.005,
            2.,
            100.,
            1e300,
        ] {
            actual.push(java_log::log(x));
        }
        for x in [
            vec![1., 2., 3., 4.],
            vec![1e16, 1., 1., -1e16],
            vec![0.1, 0.2, 0.3, 0.4, 0.5],
            vec![],
        ] {
            actual.extend([
                mean(&x),
                percentile(&x, 50.),
                percentile(&x, 90.),
                stream_sum(x),
            ]);
        }
        for (a, e) in actual.into_iter().zip(expected) {
            if a.is_nan() {
                assert!(f64::from_bits(e).is_nan());
            } else {
                assert_eq!(a.to_bits(), e);
            }
        }
        assert_eq!(java_round(0.5_f64.next_down()), 0);
        assert_eq!(java_round(-0.5), 0);
        assert_eq!(java_round(-1.5), -1);
        assert_eq!(java_round(f64::NAN), 0);
        for (x, s) in [
            (0.0, "0.0"),
            (-0.0, "-0.0"),
            (0.001, "0.001"),
            (0.0001, "1.0E-4"),
            (1e7, "1.0E7"),
            (0.45404420727340083, "0.45404420727340083"),
            (f64::INFINITY, "Infinity"),
        ] {
            assert_eq!(javafmt::double_to_string(x), s);
        }
    }
    #[test]
    fn reductions_and_search() {
        assert_eq!(stream_sum([1e16, 1., 1., -1e16]), 2.0); // Java DoubleStream.sum = 2
        assert_eq!(stream_sum([f64::INFINITY, 1.]), f64::INFINITY);
        assert_eq!(mean(&[1., 2., 3., 4.]), 2.5);
        assert_eq!(percentile(&[1., 2., 3., 4.], 50.), 2.5);
        assert_eq!(percentile(&[1., 2., 3., 4.], 90.), 4.);
        assert_eq!(
            binary_search_zero(|x| x - 0.25, 0., 1.),
            Some(0.2499995231628418)
        );
        assert!(binary_search_zero(|x| x + 1., 0., 1.).is_none());
        assert!(brent::brent_optimizer(|x| x, 0., 1., 0.5, 1e-4, 1e-4, 1).is_err());
    }
}
