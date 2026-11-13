//! Quality arithmetic for BQSR: empirical quality (GATK `RecalDatum`,
//! commons-math 3.5 `SaddlePointExpansion`) and the quality quantizer
//! (`QualQuantizer`). Section numbers (§) refer to `docs/spec/bqsr.md`.
//!
//! Java semantics replicated on purpose: `(int)` casts truncate toward zero
//! and saturate (Rust `as` does the same), `Math.round` is half-up, and the
//! argmax keeps the first maximum. Java's `Math.pow`/`log10` and commons-math
//! `FastMath.log` are replaced by Rust's `powf`/`log10`/`ln`; all are within
//! 1 ulp, which only matters at exact near-ties (gated against goldens).

pub const MAX_RECALIBRATED_Q_SCORE: i32 = 93;
const MAX_REASONABLE_Q_SCORE: i32 = 60;
const MAX_GATK_USABLE_Q_SCORE: i32 = 40;
const MAX_NUMBER_OF_OBSERVATIONS: i64 = i32::MAX as i64 - 1;

/// `QualityUtils.qualToErrorProb(double)`.
pub fn qual_to_error_prob(q: f64) -> f64 {
    10f64.powf(q / -10.0)
}

/// commons-math 3.5 `SaddlePointExpansion.EXACT_STIRLING_ERRORS` (z = 0, 0.5, ..., 15).
// Copied verbatim from the source; the decimal literals parse to the same doubles as in Java.
#[allow(clippy::excessive_precision)]
const EXACT_STIRLING_ERRORS: [f64; 31] = [
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

// Series coefficients copied verbatim from commons-math 3.5.
#[allow(clippy::excessive_precision)]
fn stirling_error(z: f64) -> f64 {
    if z < 15.0 {
        let z2 = 2.0 * z;
        assert!(z2.floor() == z2, "non-half-integer Stirling argument {z}");
        EXACT_STIRLING_ERRORS[z2 as usize]
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
        let mut j = 1i32;
        while s1 != s {
            s = s1;
            ej *= v;
            s1 = s + ej / ((j * 2) + 1) as f64;
            j += 1;
        }
        s1
    } else {
        x * (x / mu).ln() + mu - x
    }
}

/// `SaddlePointExpansion.logBinomialProbability(x, n, p, q)`.
fn saddle_log_binomial(x: i32, n: i32, p: f64, q: f64) -> f64 {
    let (xf, nf) = (x as f64, n as f64);
    if x == 0 {
        if p < 0.1 {
            -deviance_part(nf, nf * q) - nf * p
        } else {
            nf * q.ln()
        }
    } else if x == n {
        if q < 0.1 {
            -deviance_part(nf, nf * p) - nf * q
        } else {
            nf * p.ln()
        }
    } else {
        let ret = stirling_error(nf)
            - stirling_error(xf)
            - stirling_error(nf - xf)
            - deviance_part(xf, nf * p)
            - deviance_part(nf - xf, nf * q);
        let f = (2.0 * std::f64::consts::PI * xf * (nf - xf)) / nf;
        -0.5 * f.ln() + ret
    }
}

/// `BinomialDistribution(n, p).logProbability(k)`.
fn log_binomial_probability(n: i32, k: i32, p: f64) -> f64 {
    if n == 0 {
        return if k == 0 { 0.0 } else { f64::NEG_INFINITY };
    }
    if k < 0 || k > n {
        f64::NEG_INFINITY
    } else {
        saddle_log_binomial(k, n, p, 1.0 - p)
    }
}

/// `RecalDatum.getLogBinomialLikelihood`, including the > 2^31 scaling.
fn log_binomial_likelihood(q: f64, mut n: i64, mut k: i64) -> f64 {
    if n == 0 {
        return 0.0;
    }
    if n > MAX_NUMBER_OF_OBSERVATIONS {
        let fraction = MAX_NUMBER_OF_OBSERVATIONS as f64 / n as f64;
        k = java_round(k as f64 * fraction);
        n = MAX_NUMBER_OF_OBSERVATIONS;
    }
    let ll = log_binomial_probability(n as i32, k as i32, qual_to_error_prob(q));
    if ll.is_finite() { ll } else { -f64::MAX }
}

/// Java `Math.round(double)` for finite values: round half up, to `long`.
fn java_round(x: f64) -> i64 {
    (x + 0.5).floor() as i64
}

fn log_prior_cache() -> &'static [f64; (MAX_GATK_USABLE_Q_SCORE + 1) as usize] {
    static CACHE: std::sync::OnceLock<[f64; 41]> = std::sync::OnceLock::new();
    CACHE.get_or_init(|| {
        // NormalDistribution(mean 0, sd 0.5).logDensity(i)
        let sd = 0.5f64;
        let log_sd_plus_half_log_2pi = sd.ln() + 0.5 * (2.0 * std::f64::consts::PI).ln();
        std::array::from_fn(|i| {
            let x1 = (i as f64 - 0.0) / sd;
            -0.5 * x1 * x1 - log_sd_plus_half_log_2pi
        })
    })
}

fn log_prior(q: f64, prior: f64) -> f64 {
    let diff = ((q - prior) as i32)
        .unsigned_abs()
        .min(MAX_GATK_USABLE_Q_SCORE as u32);
    log_prior_cache()[diff as usize]
}

/// `RecalDatum.bayesianEstimateOfEmpiricalQuality`: argmax over Q0..Q60, first max wins.
pub fn bayesian_empirical_quality(n: i64, k: i64, prior: f64) -> i32 {
    let mut best = 0usize;
    let mut best_v = f64::NEG_INFINITY;
    for q in 0..=MAX_REASONABLE_Q_SCORE {
        let v = log_prior(q as f64, prior) + log_binomial_likelihood(q as f64, n, k);
        if q == 0 || v > best_v {
            best_v = v;
            best = q as usize;
        }
    }
    best as i32
}

/// `RecalDatum.calcEmpiricalQuality`: one smoothing error and non-error, capped at 93.
pub fn empirical_quality(observations: u64, errors: f64, prior: f64) -> i32 {
    let k = (errors + 0.5) as i64 + 1;
    let n = observations as i64 + 2;
    bayesian_empirical_quality(n, k, prior).min(MAX_RECALIBRATED_Q_SCORE)
}

/// `MathUtils.roundToNDecimalPlaces(x, 4)`: `Math.round((x + ulp(x)) * 1e4) / 1e4`.
pub fn round4(x: f64) -> f64 {
    let ulp = x.abs().next_up() - x.abs();
    java_round((x + ulp) * 10000.0) as f64 / 10000.0
}

/// `QualityUtils.errorProbToQual(rate)` (max 93), incl. the +Inf saturation to 93.
fn error_prob_to_qual(rate: f64) -> u8 {
    let d = (-10.0 * rate.log10()).round();
    // Java: Math.round -> long (saturating), then (int) of that double (saturating).
    let as_long = if d.is_nan() { 0 } else { d as i64 };
    let q = (as_long as f64) as i32;
    q.clamp(1, 93) as u8
}

// ---------------------------------------------------------------------------
// QualQuantizer (§6)

struct Interval {
    q_start: usize,
    q_end: usize,
    n_obs: i64,
    n_errors: i64,
    fixed_qual: Option<u8>,
    children: Option<(usize, usize)>,
}

struct Quantizer {
    nodes: Vec<Interval>,
    min_interesting_qual: usize,
}

impl Quantizer {
    fn error_rate(&self, i: usize) -> f64 {
        let n = &self.nodes[i];
        match n.fixed_qual {
            Some(q) => qual_to_error_prob(q as f64),
            None if n.n_obs == 0 => 0.0,
            None => (n.n_errors + 1) as f64 / (1.0 * (n.n_obs + 1) as f64),
        }
    }

    fn calc_penalty(&self, i: usize, global_rate: f64) -> f64 {
        if global_rate == 0.0 {
            return 0.0;
        }
        let n = &self.nodes[i];
        match n.children {
            None => {
                if n.q_end <= self.min_interesting_qual {
                    0.0
                } else {
                    (self.error_rate(i).log10() - global_rate.log10()).abs() * n.n_obs as f64
                }
            }
            Some((l, r)) => self.calc_penalty(l, global_rate) + self.calc_penalty(r, global_rate),
        }
    }

    fn qual(&self, i: usize) -> u8 {
        match self.nodes[i].fixed_qual {
            Some(q) => q,
            None => error_prob_to_qual(self.error_rate(i)),
        }
    }
}

/// `QualQuantizer(hist, nLevels, minInterestingQual).getOriginalToQuantizedMap()`.
pub fn quantize(hist: &[u64], n_levels: usize, min_interesting_qual: usize) -> Vec<u8> {
    let mut qz = Quantizer {
        nodes: Vec::new(),
        min_interesting_qual,
    };
    let mut current: Vec<usize> = Vec::new();
    for (q, &n_obs) in hist.iter().enumerate() {
        let n_errors = (n_obs as f64 * qual_to_error_prob(q as f64)).floor() as i32 as i64;
        qz.nodes.push(Interval {
            q_start: q,
            q_end: q,
            n_obs: n_obs as i64,
            n_errors,
            fixed_qual: Some(q as u8),
            children: None,
        });
        current.push(qz.nodes.len() - 1);
    }
    while current.len() > n_levels {
        let mut best: Option<(usize, f64)> = None;
        for pair in 0..current.len() - 1 {
            let (l, r) = (current[pair], current[pair + 1]);
            let merged = Interval {
                q_start: qz.nodes[l].q_start,
                q_end: qz.nodes[r].q_end,
                n_obs: qz.nodes[l].n_obs + qz.nodes[r].n_obs,
                n_errors: qz.nodes[l].n_errors + qz.nodes[r].n_errors,
                fixed_qual: None,
                children: Some((l, r)),
            };
            qz.nodes.push(merged);
            let m = qz.nodes.len() - 1;
            let penalty = qz.calc_penalty(m, qz.error_rate(m));
            qz.nodes.pop();
            if best.is_none_or(|(_, p)| penalty < p) {
                best = Some((pair, penalty));
            }
        }
        let (pair, _) = best.unwrap();
        let (l, r) = (current[pair], current[pair + 1]);
        qz.nodes.push(Interval {
            q_start: qz.nodes[l].q_start,
            q_end: qz.nodes[r].q_end,
            n_obs: qz.nodes[l].n_obs + qz.nodes[r].n_obs,
            n_errors: qz.nodes[l].n_errors + qz.nodes[r].n_errors,
            fixed_qual: None,
            children: Some((l, r)),
        });
        current.splice(pair..pair + 2, [qz.nodes.len() - 1]);
    }
    let mut map = vec![0u8; hist.len()];
    for &i in &current {
        let qual = qz.qual(i);
        for slot in &mut map[qz.nodes[i].q_start..=qz.nodes[i].q_end] {
            *slot = qual;
        }
    }
    map
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empirical_qualities_match_tier0_golden_rows() {
        // RecalTable1 / RecalTable0 rows of the tier0 golden (prior = reported quality).
        assert_eq!(empirical_quality(1_320_847, 100_924.0, 11.0), 11);
        assert_eq!(empirical_quality(2_851_349, 5_427.0, 25.0), 27);
        assert_eq!(empirical_quality(87_437_550, 14_884.0, 37.0), 38);
        assert_eq!(empirical_quality(91_609_746, 121_235.0, 28.4341), 29);
    }

    #[test]
    fn large_counts_take_the_scaling_branch() {
        // tier1 RecalTable0: 4,022,916,294 observations (> 2^31 - 2).
        assert_eq!(empirical_quality(4_022_916_294, 5_286_607.0, 28.3912), 29);
    }

    #[test]
    fn rounding_and_saturation() {
        assert_eq!(round4(28.43414999), 28.4341);
        assert_eq!(round4(28.43415), 28.4342);
        assert_eq!(error_prob_to_qual(0.0), 93);
        assert_eq!(error_prob_to_qual(0.001), 30);
    }

    #[test]
    fn quantizer_reproduces_tier0_pattern() {
        let mut hist = vec![0u64; 94];
        hist[11] = 1_320_847;
        hist[27] = 2_851_349;
        hist[38] = 87_437_550;
        let map = quantize(&hist, 16, 6);
        assert!(map[..=10].iter().all(|&q| q == 93));
        assert_eq!((map[11], map[27], map[38]), (11, 27, 38));
        let distinct_intervals = 1 + map.windows(2).filter(|w| w[0] != w[1]).count();
        assert!(distinct_intervals <= 16);
    }
}
