//! GATK 4.6.2.0 mutect/clustering, including the initialization split bug.
use super::numerics::{log_sum, normalize, posterior};
use crate::math::{
    binomial::binomial_probability,
    gamma::{digamma, log_gamma},
    stream_sum,
};
use crate::orientation::{jmath, math::beta_binomial};
use std::collections::BTreeMap;
#[derive(Clone, Copy)]
pub(super) struct Datum {
    pub lod: f64,
    pub artifact: f64,
    pub non_seq: f64,
    pub k: u32,
    pub n: u32,
    pub len: i32,
}
#[derive(Clone)]
struct Cluster {
    a: f64,
    b: f64,
    discrete: bool,
}
impl Cluster {
    fn fuzzy(mean: f64) -> Self {
        let m = mean.min(0.99);
        let t = (1.0 - m) / (m * (0.01 * 0.01)) - 1.0;
        let a = m * t;
        Self {
            a,
            b: t - a,
            discrete: true,
        }
    }
    fn log_likelihood(&self, n: u32, k: u32) -> f64 {
        beta_binomial(n, k, self.a, self.b)
    }
    fn corrected(&self, d: &Datum) -> f64 {
        fn g(a: f64, b: f64) -> f64 {
            log_gamma(a + b) - (log_gamma(a) + log_gamma(b))
        }
        let correction =
            g(self.a, self.b) - g(self.a + d.k as f64, self.b + (d.n - d.k) as f64) - g(1.0, 1.0)
                + g(1.0 + d.k as f64, 1.0 + (d.n - d.k) as f64);
        d.lod + correction
    }
    fn learn(&mut self, data: &[Datum], r: &[f64]) {
        if self.discrete {
            let alt = stream_sum(data.iter().zip(r).map(|(d, r)| d.k as f64 * r)) + 0.0001;
            let total = stream_sum(data.iter().zip(r).map(|(d, r)| d.n as f64 * r)) + 0.0001;
            *self = Self::fuzzy(alt / total);
            return;
        }
        for _ in 0..10 {
            for (d, r) in data.iter().zip(r) {
                let dt = digamma(d.n as f64 + self.a + self.b);
                let ds = digamma(self.a + self.b);
                let ag = digamma(self.a + d.k as f64) - dt - digamma(self.a) + ds;
                let bg = digamma(self.b + (d.n - d.k) as f64) - dt - digamma(self.b) + ds;
                self.a = (self.a + 0.01 * ag * r).max(1.0);
                self.b = (self.b + 0.01 * bg * r).max(0.5);
            }
        }
    }
}
pub(super) struct Model {
    clusters: Vec<Cluster>,
    weights: Vec<f64>,
    priors: BTreeMap<i32, f64>,
    pub artifact_prior: f64,
    pub data: Vec<Datum>,
    pub obvious: u32,
    callable: f64,
    initialized: bool,
}
impl Model {
    pub fn new(callable: f64) -> Self {
        let mut priors: BTreeMap<_, _> = (-10..=10)
            .map(|n| (n, -7.0 * std::f64::consts::LN_10))
            .collect();
        priors.insert(0, -6.0 * std::f64::consts::LN_10);
        Self {
            clusters: vec![
                Cluster {
                    a: 1.0,
                    b: 1.0,
                    discrete: false,
                },
                Cluster {
                    a: 10.0,
                    b: 1.0,
                    discrete: false,
                },
            ],
            weights: vec![jmath::log1p(0.01), jmath::log(0.01)],
            priors,
            artifact_prior: -std::f64::consts::LN_10,
            data: Vec::new(),
            obvious: 0,
            callable,
            initialized: false,
        }
    }
    pub fn register_length(&mut self, len: i32) {
        if !self.priors.contains_key(&len) {
            let minimum = self.priors.values().copied().fold(f64::INFINITY, f64::min);
            self.priors.insert(len, minimum);
        }
    }
    pub fn prior(&self, len: i32) -> f64 {
        self.priors
            .get(&len)
            .copied()
            .unwrap_or_else(|| self.priors.values().copied().fold(f64::INFINITY, f64::min))
            // MathUtils.LOG_ONE_THIRD is -Math.log(3.0), not Math.log(1/3).
            + if len == 0 { -jmath::log(3.0) } else { 0.0 }
    }
    pub fn sequencing(&self, d: &Datum) -> f64 {
        posterior(
            log_sum(
                &self
                    .clusters
                    .iter()
                    .zip(&self.weights)
                    .map(|(c, w)| w + c.corrected(d))
                    .collect::<Vec<_>>(),
            ),
            self.prior(d.len),
        )
    }
    fn somatic(&self, d: &Datum) -> f64 {
        (1.0 - d.artifact) * (1.0 - d.non_seq) * (1.0 - self.sequencing(d))
    }
    fn logs(&self, n: u32, k: u32) -> Vec<f64> {
        self.clusters
            .iter()
            .zip(&self.weights)
            .map(|(c, w)| w + c.log_likelihood(n, k))
            .collect()
    }
    pub fn likelihood(&self, n: u32, k: u32) -> f64 {
        log_sum(&self.logs(n, k))
    }
    fn em(&mut self, update: bool) {
        let mut counts: BTreeMap<i32, f64> = (-10..=10).map(|i| (i, 0.0)).collect();
        let mut insertion: Vec<_> = (-10..=10).collect();
        let mut responsibilities = Vec::new();
        let mut totals = vec![0.0; self.clusters.len()];
        for d in &self.data {
            let s = self.somatic(d);
            if !counts.contains_key(&d.len) {
                insertion.push(d.len);
            }
            *counts.entry(d.len).or_default() += s;
            let r: Vec<_> = normalize(&self.logs(d.n, d.k))
                .into_iter()
                .map(|p| s * p)
                .collect();
            for (t, r) in totals.iter_mut().zip(&r) {
                *t += r;
            }
            responsibilities.push(r);
        }
        for t in &mut totals {
            *t += 1.0;
        }
        let sum: f64 = totals.iter().sum();
        self.weights = totals.iter().map(|v| jmath::log(v / sum)).collect();
        let artifacts = self.obvious as f64 + stream_sum(self.data.iter().map(|d| d.artifact));
        // Collectors.toMap's two-argument overload uses putIfAbsent (JDK 17),
        // as do the added lengths. Resize splitting retains chain order.
        // This supported key set does not treeify.
        let mut capacity = 16;
        while counts.len() > capacity * 3 / 4 {
            capacity *= 2;
        }
        insertion.sort_by_key(|n| {
            let h = *n as u32;
            (h ^ (h >> 16)) as usize & (capacity - 1)
        });
        let count = stream_sum(insertion.iter().map(|n| counts[n]));
        if update {
            self.artifact_prior = jmath::log((count + 1.0) / (count + artifacts + 2.0));
            for n in -10..=10 {
                self.priors.insert(
                    n,
                    jmath::log((counts[&n] / self.callable).max(if n == 0 { 1e-8 } else { 1e-9 })),
                );
            }
        }
        for (i, c) in self.clusters.iter_mut().enumerate() {
            c.learn(
                &self.data,
                &responsibilities.iter().map(|r| r[i]).collect::<Vec<_>>(),
            );
        }
    }
    fn initialize(&mut self) {
        let somatic: Vec<_> = self.data.iter().map(|d| self.somatic(d)).collect();
        let mut previous = f64::NEG_INFINITY;
        for _ in 0..5 {
            let old = self.weights.clone();
            let bg: Vec<_> = self
                .data
                .iter()
                .zip(&somatic)
                .map(|(d, s)| s * normalize(&self.logs(d.n, d.k))[0])
                .collect();
            let mut af: Vec<_> = self
                .data
                .iter()
                .map(|d| (d.k as f64 / d.n as f64, self.somatic(d)))
                .collect();
            af.sort_by(|a, b| a.0.total_cmp(&b.0));
            let step = stream_sum(af.iter().map(|p| p.1)) / 50.0;
            if step == 0.0 {
                break;
            }
            let mut next = step;
            let mut cumulative = 0.0;
            let mut quantiles = Vec::new();
            for (f, s) in af {
                cumulative += s;
                if cumulative > next {
                    if quantiles.last() != Some(&f) {
                        quantiles.push(f);
                    }
                    while cumulative > next {
                        next += step;
                    }
                }
            }
            let mut density = vec![0.0; quantiles.len()];
            for (d, bg) in self.data.iter().zip(bg) {
                for (v, &q) in density.iter_mut().zip(&quantiles) {
                    *v += binomial_probability(d.n as i32, d.k as i32, q) * bg * (d.n + 1) as f64;
                }
            }
            let (mut mass, mut peak, mut peak_r) = (0.0, 0.0, 0.0);
            let mut peaks = Vec::new();
            for q in 0..quantiles.len() {
                let left = if q == 0 { 0.0 } else { density[q - 1] };
                let right = if q + 1 == quantiles.len() {
                    0.0
                } else {
                    density[q + 1]
                };
                let left_af = if q == 0 { 0.0 } else { quantiles[q - 1] };
                let r = density[q];
                mass += (quantiles[q] - left_af) * (left + r) / 2.0;
                if r > peak_r {
                    peak = quantiles[q];
                    peak_r = r;
                }
                let min = (r < left && r <= right) || (r <= left && r < right);
                if (min && q > 0) || q + 1 == quantiles.len() {
                    peaks.push((peak, mass));
                    mass = 0.0;
                    peak = quantiles[q];
                    peak_r = r;
                }
            }
            let total_mass = stream_sum(peaks.iter().map(|p| p.1));
            peaks.sort_by(|a, b| b.1.total_cmp(&a.1));
            let Some(&(peak, mass)) = peaks.first() else {
                break;
            };
            if peak < quantiles[5.min(quantiles.len() - 1)] {
                break;
            }
            let fraction = (mass / total_mass).min(0.9);
            self.clusters.push(Cluster::fuzzy(peak));
            self.weights.push(jmath::log(fraction) + self.weights[0]);
            self.weights[0] += jmath::log1p(fraction);
            for _ in 0..5 {
                self.em(false);
            }
            let ll: f64 = self
                .data
                .iter()
                .zip(&somatic)
                .map(|(d, s)| s * self.likelihood(d.n, d.k))
                .sum();
            let effective: f64 = somatic.iter().sum();
            let bic = ll - (2 * self.clusters.len()) as f64 * jmath::log(effective);
            if bic < previous {
                self.clusters.pop();
                self.weights = old;
                break;
            }
            previous = bic;
        }
        self.initialized = true;
    }
    pub fn learn(&mut self) {
        if !self.initialized {
            self.initialize();
        }
        for _ in 0..5 {
            self.em(true);
        }
        self.data.clear();
        self.obvious = 0;
    }
    pub fn metadata(&self) -> String {
        let mut out = String::new();
        for n in -10i32..=10 {
            let name = if n == 0 {
                "SNV".to_owned()
            } else {
                format!(
                    "{} of length {}",
                    if n < 0 { "deletion" } else { "insertion" },
                    n.abs()
                )
            };
            out += &format!(
                "#<METADATA>Ln prior of {name}={}\n",
                crate::math::javafmt::double_to_string(self.priors[&n])
            );
        }
        for i in 0..2 {
            let c = &self.clusters[i];
            out += &format!(
                "#<METADATA>{} beta-binomial cluster=weight = {:.4}, alpha = {:.2}, beta = {:.2}\n",
                if i == 0 { "Background" } else { "High-AF" },
                jmath::exp(self.weights[i]),
                c.a,
                c.b
            );
        }
        let mut indices: Vec<_> = (2..self.clusters.len()).collect();
        indices.sort_by(|&a, &b| self.weights[b].total_cmp(&self.weights[a]));
        for i in indices {
            let c = &self.clusters[i];
            out += &format!(
                "#<METADATA>Binomial cluster=weight = {:.4}, mean = {:.3}\n",
                jmath::exp(self.weights[i]),
                c.a / (c.a + c.b)
            );
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn snv_offset_uses_negative_log_three() {
        let mut model = Model::new(100.0);
        model.priors.insert(0, 0.0);
        // MathUtils.LOG_ONE_THIRD, from FilteringLearning.java. log(1/3)
        // differs by one bit, even though the initial prior can hide it.
        assert_eq!(model.prior(0).to_bits(), 0xbff1_93ea_7aad_030a);
        assert_ne!(model.prior(0).to_bits(), jmath::log(1.0 / 3.0).to_bits());
    }
}
