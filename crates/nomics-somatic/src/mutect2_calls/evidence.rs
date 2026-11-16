//! GATK 4.6.2.0 SomaticLikelihoodsEngine, default unweighted evidence.
//! Matrices are allele-major and in natural-log space. Reduction order is part
//! of the contract: only the entropy term uses DoubleStream's compensated sum.
use crate::math::{gamma, stream_sum};
use nomics_hc::jmath::{exp, log};

fn sum(xs: impl IntoIterator<Item = f64>) -> f64 {
    xs.into_iter().fold(0.0, |a, b| a + b)
}

/// NaturalLogUtils.logSumExp (first maximum; ordered, uncompensated loop).
pub fn log_sum_exp(xs: &[f64]) -> f64 {
    assert!(!xs.is_empty());
    let mut best = 0;
    for i in 1..xs.len() {
        if xs[i] > xs[best] {
            best = i;
        }
    }
    let max = xs[best];
    if max == f64::NEG_INFINITY {
        return max;
    }
    let mut total = 1.0;
    for (i, &v) in xs.iter().enumerate() {
        if i != best && v != f64::NEG_INFINITY {
            total += exp(v - max);
        }
    }
    assert!(total.is_finite(), "invalid likelihoods in logSumExp");
    max + if total != 1.0 { log(total) } else { 0.0 }
}

pub fn normalize(xs: &[f64]) -> Vec<f64> {
    let total = log_sum_exp(xs);
    xs.iter().map(|&x| exp(x - total)).collect()
}

fn weights(alpha: &[f64]) -> Vec<f64> {
    let psi_sum = gamma::digamma(sum(alpha.iter().copied()));
    alpha.iter().map(|&a| gamma::digamma(a) - psi_sum).collect()
}

fn responsibilities(matrix: &[Vec<f64>], column: usize, weights: &[f64]) -> Vec<f64> {
    normalize(
        &matrix
            .iter()
            .zip(weights)
            .map(|(r, &w)| w + r[column])
            .collect::<Vec<_>>(),
    )
}

fn validate(matrix: &[Vec<f64>], prior: &[f64]) -> usize {
    assert!(!matrix.is_empty());
    assert_eq!(matrix.len(), prior.len());
    let n = matrix[0].len();
    assert!(matrix.iter().all(|row| row.len() == n));
    assert!(prior.iter().all(|&p| p > 0.0 && p.is_finite()));
    n
}

/// Dirichlet posterior parameters. The caller handles the empty-matrix case,
/// just as SomaticGenotypingEngine does before invoking the Java engine.
pub fn posterior(matrix: &[Vec<f64>], prior: &[f64]) -> Vec<f64> {
    let n = validate(matrix, prior);
    assert!(n > 0, "SomaticLikelihoodsEngine requires nonempty evidence");
    let mut old = vec![1.0; prior.len()];
    loop {
        let w = weights(&old);
        // sumArrayFunction initializes from column zero, rather than zeros.
        let mut counts = responsibilities(matrix, 0, &w);
        for column in 1..n {
            for (a, r) in counts.iter_mut().zip(responsibilities(matrix, column, &w)) {
                *a += r;
            }
        }
        let next: Vec<_> = counts.iter().zip(prior).map(|(&a, &b)| a + b).collect();
        let distance = sum(old.iter().zip(&next).map(|(&a, &b)| (a - b).abs()));
        let converged = distance / sum(next.iter().copied()) < 0.001;
        old = next;
        if converged {
            return old;
        }
    }
}

pub fn log_dirichlet_normalization(alpha: &[f64]) -> f64 {
    gamma::log_gamma(sum(alpha.iter().copied())) - sum(alpha.iter().map(|&a| gamma::log_gamma(a)))
}

/// Variational log evidence; this is not a diploid likelihood.
pub fn log_evidence(matrix: &[Vec<f64>], prior: &[f64]) -> f64 {
    let alpha = posterior(matrix, prior);
    let prior_term = log_dirichlet_normalization(prior);
    let posterior_term = -log_dirichlet_normalization(&alpha);
    let w = weights(&alpha);
    let contribution = sum((0..matrix[0].len()).map(|column| {
        let r = responsibilities(matrix, column, &w);
        let entropy = stream_sum(r.iter().map(|&x| if x < 1e-8 { 0.0 } else { x * log(x) }));
        let likelihood = sum(r.iter().enumerate().map(|(a, &x)| {
            if x < 1e-10 {
                0.0
            } else {
                matrix[a][column] * x
            }
        }));
        likelihood - entropy
    }));
    prior_term + posterior_term + contribution
}

/// SomaticGenotypingEngine.somaticLogOdds, reference at index zero, minAF=0.
pub fn somatic_log_odds(matrix: &[Vec<f64>]) -> Vec<f64> {
    validate(matrix, &vec![1.0; matrix.len()]);
    if matrix[0].is_empty() {
        return vec![0.0; matrix.len() - 1];
    }
    let all = log_evidence(matrix, &vec![1.0; matrix.len()]);
    (1..matrix.len())
        .map(|alt| {
            let without: Vec<_> = matrix
                .iter()
                .enumerate()
                .filter(|(i, _)| *i != alt)
                .map(|(_, r)| r.clone())
                .collect();
            all - log_evidence(&without, &vec![1.0; without.len()])
        })
        .collect()
}

/// Hom-ref minus heterozygous normal likelihood. There is no hom-alt term.
pub fn diploid_alt_log_odds(matrix: &[Vec<f64>]) -> Vec<f64> {
    validate(matrix, &vec![1.0; matrix.len()]);
    let hom_ref = sum(matrix[0].iter().copied());
    matrix[1..]
        .iter()
        .map(|alt| {
            hom_ref
                - sum(matrix[0]
                    .iter()
                    .zip(alt)
                    .map(|(&r, &a)| log_sum_exp(&[r, a]) + log(0.5)))
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn gatk_independent_bit_patterns() {
        // Generated by MutectSomaticNumerics.java against the unmodified jar.
        // Input and expected output bits share each row; none are Rust-derived.
        for line in include_str!("evidence.tsv").lines() {
            let f: Vec<_> = line.split_whitespace().collect();
            let a: usize = f[0].parse().unwrap();
            let n: usize = f[1].parse().unwrap();
            let bits: Vec<_> = f[2..]
                .iter()
                .map(|s| u64::from_str_radix(s, 16).unwrap())
                .collect();
            let m: Vec<Vec<_>> = (0..a)
                .map(|i| {
                    bits[i * n..(i + 1) * n]
                        .iter()
                        .map(|&b| f64::from_bits(b))
                        .collect()
                })
                .collect();
            let prior: Vec<_> = bits[a * n..a * n + a]
                .iter()
                .map(|&b| f64::from_bits(b))
                .collect();
            let mut actual = posterior(&m, &prior);
            actual.push(log_evidence(&m, &prior));
            actual.push(log_dirichlet_normalization(&prior));
            actual.extend(somatic_log_odds(&m));
            actual.extend(diploid_alt_log_odds(&m));
            for (i, (v, &e)) in actual.iter().zip(&bits[a * n + a..]).enumerate() {
                assert_eq!(v.to_bits(), e, "{a}x{n}, result {i}: {v}");
            }
            assert_eq!(actual.len(), bits.len() - a * n - a);
        }
        assert_eq!(somatic_log_odds(&[vec![], vec![]]), vec![0.0]);
        assert_eq!(diploid_alt_log_odds(&[vec![], vec![]]), vec![0.0]);
    }
}
