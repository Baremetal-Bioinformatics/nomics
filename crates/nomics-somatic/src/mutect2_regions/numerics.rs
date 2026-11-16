//! GATK Mutect2Engine.logLikelihoodRatio, flat beta prior only.
use crate::{
    math::{fastmath, gamma},
    orientation::{jmath, math::binomial_log},
};
use std::sync::OnceLock;

/// Operation-ordered variational activity likelihood (natural logs).
pub fn log_likelihood_ratio(n_ref: usize, alt_quals: &[u8]) -> f64 {
    repeated_log_likelihood_ratio(n_ref, alt_quals, 1)
}

/// AdaptiveChainPruner's constant-error overload: quality rounding precedes LLR.
pub fn constant_error_log_likelihood_ratio(n_ref: i32, n_alt: i32, error: f64) -> f64 {
    let q = ((-10.0 * nomics_hc::jmath::log10(error) + 0.5).floor() as i32).clamp(1, 93) as u8;
    repeated_log_likelihood_ratio(n_ref as usize, &[q], n_alt as usize)
}

fn repeated_log_likelihood_ratio(n_ref: usize, alt_quals: &[u8], repeat: usize) -> f64 {
    static QUALS: OnceLock<[(f64, f64, f64); 255]> = OnceLock::new();
    let quals = QUALS.get_or_init(|| {
        std::array::from_fn(|q| {
            let epsilon = jmath::pow(10.0, -(q as f64) / 10.0);
            let log_epsilon = q as f64 * (-jmath::log(10.0) / 10.0);
            let log_prob = if log_epsilon == 0.0 {
                f64::NEG_INFINITY
            } else if log_epsilon < jmath::log(0.5) {
                jmath::log1p(-jmath::exp(log_epsilon))
            } else {
                jmath::log(-jmath::expm1(log_epsilon))
            };
            (epsilon, log_epsilon, log_prob)
        })
    });
    let n_alt = alt_quals.len() * repeat;
    let n = n_ref + n_alt;
    let ratio =
        fastmath::exp(gamma::digamma((n_ref + 1) as f64) - gamma::digamma((n_alt + 1) as f64));
    let mut read_sum = 0.0;
    for &q in alt_quals {
        let (e, log_e, log_p) = quals[q as usize];
        let z = (1.0 - e) / (1.0 - e + e * ratio);
        let p = z * (1.0 - z);
        read_sum += z * (log_p - log_e) + p * (11.0 + 33.0 * p) / (2.0 + 20.0 * p);
    }
    -jmath::log((n + 1) as f64) - binomial_log(n as u32, n_alt as u32) + read_sum * repeat as f64
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn unmodified_jar_likelihood_bits() {
        for line in include_str!("numerics_oracle.tsv")
            .lines()
            .filter(|l| l.starts_with("LLR\t"))
        {
            let f: Vec<_> = line.split('\t').collect();
            let r = f[1].parse().unwrap();
            let qs: Vec<u8> = f[2]
                .split(',')
                .filter(|s| !s.is_empty())
                .map(|s| s.parse().unwrap())
                .collect();
            let expected = u64::from_str_radix(f[3], 16).unwrap();
            let actual = log_likelihood_ratio(r, &qs);
            assert_eq!(
                actual.to_bits(),
                expected,
                "nRef={r}, quals={qs:?}: {actual}"
            );
        }
    }
    #[test]
    fn unmodified_jar_smoothing_and_cuts() {
        use nomics_hc::activity::BinaryActivityProfile;
        let mut profile = BinaryActivityProfile::new(1000);
        let mut regions = Vec::new();
        for pos in 1..=1000 {
            regions.extend(profile.push(0, pos, if pos == 400 || pos == 600 { 1.0 } else { 0.0 }));
        }
        regions.extend(profile.finish(0));
        let actual: String = regions
            .iter()
            .map(|r| format!("REGION\t{}\t{}\t{}\n", r.start, r.end, r.active))
            .collect();
        let expected: String = include_str!("numerics_oracle.tsv")
            .lines()
            .filter(|l| l.starts_with("REGION\t"))
            .map(|l| format!("{l}\n"))
            .collect();
        assert_eq!(actual, expected);
    }
}
