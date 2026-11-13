//! Frozen fragment-length model; ascending integer-count observations followed
//! by the same two normalization boundaries as BAM -> RAD -> inference.
use statrs::distribution::{Binomial, ContinuousCDF, Discrete, Normal};

pub(crate) const MAX: usize = 1000;

pub(crate) fn log_add(a: f64, b: f64) -> f64 {
    if a == f64::NEG_INFINITY {
        return b;
    }
    if b == f64::NEG_INFINITY {
        return a;
    }
    a.max(b) + (a.min(b) - a.max(b)).exp().ln_1p()
}

pub(crate) struct Distribution {
    pub pmf: Vec<f64>,
    pub cmf: Vec<f64>,
    means: Vec<f64>,
}

impl Distribution {
    pub fn new(counts: &[u64]) -> Self {
        let normal = Normal::new(250.0, 25.0).unwrap();
        let mut hist: Vec<f64> = (0..=MAX)
            .map(|i| {
                let p = normal.cdf(i as f64 + 0.5) - normal.cdf(i as f64 - 0.5);
                if p == 0.0 { -1e10 } else { p.ln() }
            })
            .collect();
        let kernel = Binomial::new(0.5, 4).unwrap();
        let kernel: Vec<f64> = (0..=4).map(|i| kernel.pmf(i).ln()).collect();
        let mut total = hist.iter().copied().fold(f64::NEG_INFINITY, log_add);
        for (length, &count) in counts.iter().enumerate() {
            if count == 0 {
                continue;
            }
            for (k, &mass) in kernel.iter().enumerate() {
                let bin = length as isize + k as isize - 2;
                if bin > 0 && bin <= MAX as isize {
                    let contribution = (count as f64).ln() + mass;
                    hist[bin as usize] = log_add(hist[bin as usize], contribution);
                    total = log_add(total, contribution);
                }
            }
        }
        let mut pmf: Vec<f64> = hist.iter().map(|p| p - total).collect();
        normalize(&mut pmf);
        // from_log_pmf reconstructs hist/tot_mass, then cache normalizes again.
        let total = pmf.iter().copied().fold(f64::NEG_INFINITY, log_add);
        for p in &mut pmf {
            *p -= total;
        }
        normalize(&mut pmf);
        let mut sum = f64::NEG_INFINITY;
        let cmf = pmf
            .iter()
            .map(|&p| {
                sum = log_add(sum, p);
                sum
            })
            .collect();
        let (mut numerator, mut denominator) = (0.0, 0.0);
        let means = pmf
            .iter()
            .enumerate()
            .map(|(i, p)| {
                let p = p.exp();
                numerator += i as f64 * p;
                denominator += p;
                if denominator > 0.0 {
                    numerator / denominator
                } else {
                    0.0
                }
            })
            .collect();
        Self { pmf, cmf, means }
    }

    pub fn effective_length(&self, length: usize) -> f64 {
        let effective = length as f64 - self.means[length.min(MAX)];
        if effective < 1.0 {
            length as f64
        } else {
            effective
        }
    }
}

fn normalize(values: &mut [f64]) {
    let sum = values.iter().copied().fold(f64::NEG_INFINITY, log_add);
    for p in values {
        *p -= sum;
    }
}
