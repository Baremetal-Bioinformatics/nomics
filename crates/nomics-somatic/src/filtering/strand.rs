//! GATK StrandArtifactFilter E/M steps (log10 normalization is intentional).
use crate::math::{brent::brent_optimizer, stream_sum};
use crate::orientation::{
    jmath,
    math::{beta_binomial as bb, binomial_log},
};
#[derive(Clone, Copy)]
pub(super) struct Step {
    pub f: f64,
    pub r: f64,
    pub nf: u32,
    pub nr: u32,
    pub kf: u32,
    pub kr: u32,
}
impl Step {
    pub fn p(&self) -> f64 {
        self.f + self.r
    }
}
pub(super) struct Strand {
    pub a: f64,
    pub b: f64,
    pub prior: f64,
    pub steps: Vec<Step>,
}
impl Strand {
    pub fn new() -> Self {
        Self {
            a: 1.0,
            b: 20.0,
            prior: 0.001,
            steps: Vec::new(),
        }
    }
    pub fn step(&self, nf: u32, nr: u32, kf: u32, kr: u32, len: u32) -> Step {
        let mut s = Step {
            f: 0.0,
            r: 0.0,
            nf,
            nr,
            kf,
            kr,
        };
        if kf + kr == 0 || len > 4 {
            return s;
        }
        let beta = if len == 0 {
            1000.0
        } else if len < 3 {
            5000.0
        } else {
            50000.0
        };
        let f = bb(nf, kf, self.a, self.b) + bb(nr, kr, 1.0, beta);
        let r = bb(nr, kr, self.a, self.b) + bb(nf, kf, 1.0, beta);
        let none = binomial_log(nf, kf) + binomial_log(nr, kr) - binomial_log(nf + nr, kf + kr)
            + bb(nf + nr, kf + kr, 1.0, 1.0);
        let x = [
            (f + jmath::log(self.prior / 2.0)) * std::f64::consts::LOG10_E,
            (r + jmath::log(self.prior / 2.0)) * std::f64::consts::LOG10_E,
            (none + jmath::log(1.0 - self.prior)) * std::f64::consts::LOG10_E,
        ];
        let mut max = 0;
        for i in 1..3 {
            if x[i] > x[max] {
                max = i;
            }
        }
        let mut sum = 1.0;
        for i in 0..3 {
            if i != max && x[i] != f64::NEG_INFINITY {
                sum += jmath::pow(10.0, x[i] - x[max]);
            }
        }
        let logsum = x[max] + if sum != 1.0 { jmath::log10(sum) } else { 0.0 };
        s.f = jmath::pow(10.0, x[0] - logsum);
        s.r = jmath::pow(10.0, x[1] - logsum);
        s
    }
    pub fn learn(&mut self) -> std::io::Result<()> {
        let potential: Vec<_> = self.steps.iter().filter(|s| s.p() > 0.1).collect();
        let artifacts = stream_sum(potential.iter().map(|s| s.p()));
        let non = stream_sum(self.steps.iter().map(|s| 1.0 - s.p()));
        self.prior = (artifacts + 1.0) / (artifacts + 1.0 + non + 1000.0);
        let alt = stream_sum(
            potential
                .iter()
                .map(|s| s.f * s.kf as f64 + s.r * s.kr as f64),
        );
        let depth = stream_sum(
            potential
                .iter()
                .map(|s| s.f * s.nf as f64 + s.r * s.nr as f64),
        );
        let mean = (alt + 1.0) / (depth + 1.0 + 20.0);
        let objective = |a| {
            let b = (1.0 / mean - 1.0) * a;
            stream_sum(
                potential
                    .iter()
                    .map(|s| s.f * bb(s.nf, s.kf, a, b) + s.r * bb(s.nr, s.kr, a, b)),
            )
        };
        self.a = brent_optimizer(objective, 0.01, 100.0, 1.0, 0.01, 0.01, 100)
            .map_err(std::io::Error::other)?
            .0;
        self.b = (1.0 / mean - 1.0) * self.a;
        self.steps.clear();
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::filtering::numerics::{combined_error, round};

    #[test]
    fn java_scalar_and_strand_learning_bits() {
        let expected: Vec<_> = include_str!("../../tests/oracle/filtering_learning.bits")
            .lines()
            .map(|s| u64::from_str_radix(s, 16).unwrap())
            .collect();
        let mut actual = vec![-jmath::log(3.0), jmath::log(1.0 / 3.0)];
        for x in [-20.0, -0.7, -0.1, -1e-8, 0.0, 0.01, 0.9] {
            actual.extend([jmath::log1p(x), jmath::expm1(x)]);
        }
        for x in [-300.0, -17.61, -5.186226408608968, -0.1, 0.0, 2.0, 3.0, 8.0] {
            actual.push(jmath::pow(10.0, x));
        }
        for x in [1e-300, 0.001, 0.1, 1.0, 1.5, 3.0, 1000.0] {
            actual.push(jmath::log10(x));
        }
        let p = [
            0.04163411606541967,
            0.006318578253251994,
            0.03865375732824039,
        ];
        for order in [[0, 1, 2], [0, 2, 1]] {
            let c = combined_error(p, order);
            actual.extend([c, round(c, 3).parse().unwrap()]);
        }
        let mut strand = Strand::new();
        for pass in 0..3 {
            for [nf, nr, kf, kr, len] in [
                [1, 37, 0, 15, 0],
                [6, 23, 1, 1, 0],
                [70, 45, 28, 0, 0],
                [50, 50, 0, 21, 1],
                [90, 70, 35, 2, 2],
                [100, 100, 1, 40, 3],
                [20, 30, 10, 0, 4],
            ] {
                let s = strand.step(nf, nr, kf, kr, len);
                actual.extend([s.f, s.r]);
                strand.steps.push(s);
            }
            if pass < 2 {
                strand.learn().unwrap();
                actual.extend([strand.a, strand.b, strand.prior]);
            }
        }
        assert_eq!(actual.len(), expected.len());
        for (i, (a, e)) in actual.into_iter().zip(expected).enumerate() {
            if f64::from_bits(e).is_nan() {
                assert!(a.is_nan(), "Java learning oracle {i}");
            } else {
                assert_eq!(a.to_bits(), e, "Java learning oracle {i}: {a}");
            }
        }
    }
}
