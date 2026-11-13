//! Commons Math 3.5 distribution/{BinomialDistribution,SaddlePointExpansion}.java.
//! Only integer trial/count inputs are called; the non-half-integer Gamma branch is unreachable.
#![allow(clippy::excessive_precision)]
use super::fastmath::{exp, log};
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
pub fn stirling_error(z: f64) -> f64 {
    if z < 15.0 {
        return STIRLING[(2.0 * z) as usize];
    }
    let z2 = z * z;
    (0.083333333333333333333
        - (0.00277777777777777777778
            - (0.00079365079365079365079365
                - (0.000595238095238095238095238 - 0.0008417508417508417508417508 / z2) / z2)
                / z2)
            / z2)
        / z
}
pub fn deviance_part(x: f64, mu: f64) -> f64 {
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
            s1 = s + ej / f64::from(j * 2 + 1);
            j += 1;
        }
        s1
    } else {
        x * log(x / mu) + mu - x
    }
}
pub(crate) fn log_binomial_probability(n: i32, k: i32, p: f64) -> f64 {
    if k < 0 || k > n {
        return f64::NEG_INFINITY;
    }
    if n == 0 {
        return 0.0;
    }
    let q = 1.0 - p;
    let nf = f64::from(n);
    let x = f64::from(k);
    if k == 0 {
        if p < 0.1 {
            -deviance_part(nf, nf * q) - nf * p
        } else {
            nf * log(q)
        }
    } else if k == n {
        if q < 0.1 {
            -deviance_part(nf, nf * p) - nf * q
        } else {
            nf * log(p)
        }
    } else {
        let ret = stirling_error(nf)
            - stirling_error(x)
            - stirling_error(nf - x)
            - deviance_part(x, nf * p)
            - deviance_part(nf - x, nf * q);
        let f = (std::f64::consts::TAU * x * (nf - x)) / nf;
        -0.5 * log(f) + ret
    }
}
pub fn binomial_probability(n: i32, k: i32, p: f64) -> f64 {
    let ret = log_binomial_probability(n, k, p);
    if ret == f64::NEG_INFINITY { 0.0 } else { exp(ret) }
}
