//! Operation-ordered Commons Math 3.5 Beta.logBeta / FastMath.log / log1p.
//! Shared with filtering: general positive shapes call math::gamma; the
//! orientation learner still reaches only Gamma(1) and Gamma(2).
#![allow(clippy::excessive_precision, clippy::approx_constant)]
use super::log_table::LN_MANT;
const SPLIT: f64 = 1073741824.0;
const QUICK: [[f64; 2]; 9] = [
    [1.0, 5.669184079525E-24],
    [-0.25, -0.25],
    [0.3333333134651184, 1.986821492305628E-8],
    [-0.25, -6.663542893624021E-14],
    [0.19999998807907104, 1.1921056801463227E-8],
    [-0.1666666567325592, -7.800414592973399E-9],
    [0.1428571343421936, 5.650007086920087E-9],
    [-0.12502530217170715, -7.44321345601866E-11],
    [0.11113807559013367, 9.219544613762692E-9],
];
const HIGH: [[f64; 2]; 6] = [
    [1.0, -6.032174644509064E-23],
    [-0.25, -0.25],
    [0.3333333134651184, 1.9868161777724352E-8],
    [-0.2499999701976776, -2.957007209750105E-8],
    [0.19999954104423523, 1.5830993332061267E-10],
    [-0.16624879837036133, -2.6033824355191673E-8],
];
fn split(x: f64) -> (f64, f64) {
    let tmp = x * SPLIT;
    let a = x + tmp - tmp;
    (a, x - a)
}
fn polynomial(xa: f64, xb: f64, coefficients: &[[f64; 2]]) -> (f64, f64) {
    let [mut ya, mut yb] = coefficients[coefficients.len() - 1];
    for &[ca, cb] in coefficients[..coefficients.len() - 1].iter().rev() {
        let aa = ya * xa;
        let ab = ya * xb + yb * xa + yb * xb;
        (ya, yb) = split(aa);
        yb += ab;
        let aa = ya + ca;
        let ab = yb + cb;
        (ya, yb) = split(aa);
        yb += ab;
    }
    let aa = ya * xa;
    let ab = ya * xb + yb * xa + yb * xb;
    (aa, ab)
}
pub(crate) fn log(x: f64) -> f64 {
    log_parts(x, false).0
}
fn log_parts(x: f64, high: bool) -> (f64, [f64; 2]) {
    if x == 0.0 {
        return (f64::NEG_INFINITY, [f64::NEG_INFINITY, 0.0]);
    }
    if x < 0.0 || x.is_nan() {
        return (f64::NAN, [f64::NAN, 0.0]);
    }
    if x.is_infinite() {
        return (x, [x, 0.0]);
    }
    let mut bits = x.to_bits();
    let mut exp = (bits >> 52) as i32 - 1023;
    if bits & 0x7ff0000000000000 == 0 {
        bits <<= 1;
        while bits & 0x0010000000000000 == 0 {
            exp -= 1;
            bits <<= 1;
        }
    }
    if (exp == -1 || exp == 0) && x < 1.01 && x > 0.99 && !high {
        let (xa, xb) = split(x - 1.0);
        let (aa, ab) = polynomial(xa, xb, &QUICK);
        let (ya, yb) = split(aa);
        return (ya + (yb + ab), [0.0; 2]);
    }
    let lnm = LN_MANT[((bits & 0x000ffc0000000000) >> 42) as usize];
    let numer = (bits & 0x3ffffffffff) as f64;
    let denom = 4503599627370496.0 + (bits & 0x000ffc0000000000) as f64;
    let epsilon = numer / denom;
    let (lnza, lnzb) = if high {
        let (xa, mut xb) = split(epsilon);
        let aa = numer - xa * denom - xb * denom;
        xb += aa / denom;
        let (aa, ab) = polynomial(xa, xb, &HIGH);
        let a = aa + ab;
        (a, -(a - aa - ab))
    } else {
        let mut z = -0.16624882440418567;
        z = z * epsilon + 0.19999954120254515;
        z = z * epsilon + -0.2499999997677497;
        z = z * epsilon + 0.3333333333332802;
        z = z * epsilon + -0.5;
        z = z * epsilon + 1.0;
        (z * epsilon, 0.0)
    };
    let mut a = 0.693147063255310059 * exp as f64;
    let mut b = 0.0;
    for z in [
        lnm[0],
        lnza,
        1.17304635250823482e-7 * exp as f64,
        lnm[1],
        lnzb,
    ] {
        let c = a + z;
        let d = -(c - a - z);
        a = c;
        b += d;
    }
    (a + b, [a, b])
}
pub(crate) fn log1p(x: f64) -> f64 {
    if x == -1.0 {
        return f64::NEG_INFINITY;
    }
    if x == f64::INFINITY {
        return x;
    }
    if !(-1e-6..=1e-6).contains(&x) {
        let xpa = 1.0 + x;
        let xpb = -(xpa - 1.0 - x);
        let (lores, hi) = log_parts(xpa, true);
        if lores.is_infinite() {
            return lores;
        }
        let fx1 = xpb / xpa;
        let epsilon = 0.5 * fx1 + 1.0;
        epsilon * fx1 + hi[1] + hi[0]
    } else {
        ((x * (1.0 / 3.0) - 0.5) * x + 1.0) * x
    }
}
const DELTA: [f64; 15] = [
    0.833333333333333333333333333333E-01,
    -0.277777777777777777777777752282E-04,
    0.793650793650793650791732130419E-07,
    -0.595238095238095232389839236182E-09,
    0.841750841750832853294451671990E-11,
    -0.191752691751854612334149171243E-12,
    0.641025640510325475730918472625E-14,
    -0.295506514125338232839867823991E-15,
    0.179643716359402238723287696452E-16,
    -0.139228964661627791231203060395E-17,
    0.133802855014020915603275339093E-18,
    -0.154246009867966094273710216533E-19,
    0.197701992980957427278370133333E-20,
    -0.234065664793997056856992426667E-21,
    0.171348014966398575409015466667E-22,
];
fn delta(a: f64, b: f64) -> f64 {
    let h = a / b;
    let p = h / (1.0 + h);
    let q = 1.0 / (1.0 + h);
    let q2 = q * q;
    let mut s = [1.0; 15];
    for i in 1..15 {
        s[i] = 1.0 + (q + q2 * s[i - 1]);
    }
    let sqrt_t = 10.0 / b;
    let t = sqrt_t * sqrt_t;
    let mut w = DELTA[14] * s[14];
    for i in (0..14).rev() {
        w = t * w + DELTA[i] * s[i];
    }
    w * p / b
}
fn gamma_difference(a: f64, b: f64) -> f64 {
    let (d, w) = if a <= b {
        (b + (a - 0.5), delta(a, b))
    } else {
        (a + (b - 0.5), delta(b, a))
    };
    let u = d * log1p(a / b);
    let v = a * (log(b) - 1.0);
    if u <= v { (w - u) - v } else { (w - v) - u }
}
fn gamma_reduced(x: f64) -> f64 {
    crate::math::gamma::log_gamma(x)
}
fn gamma_sum(a: f64, b: f64) -> f64 {
    let x = (a - 1.0) + (b - 1.0);
    if x <= 0.5 {
        crate::math::gamma::log_gamma1p(1.0 + x)
    } else if x <= 1.5 {
        crate::math::gamma::log_gamma1p(x) + log1p(x)
    } else {
        crate::math::gamma::log_gamma1p(x - 1.0) + log(x * (1.0 + x))
    }
}
pub(crate) fn log_beta(p: f64, q: f64) -> f64 {
    if p <= 0.0 || q <= 0.0 || p.is_nan() || q.is_nan() {
        return f64::NAN;
    }
    let a = p.min(q);
    let b = p.max(q);
    if a >= 10.0 {
        let sqrt_t = 10.0 / a;
        let t = sqrt_t * sqrt_t;
        let mut z = DELTA[14];
        for i in (0..14).rev() {
            z = t * z + DELTA[i];
        }
        let w = z / a + delta(a, b);
        let h = a / b;
        let c = h / (1.0 + h);
        let u = -(a - 0.5) * log(c);
        let v = b * log1p(h);
        let r = (-0.5 * log(b) + 0.9189385332046727) + w;
        if u <= v { (r - u) - v } else { (r - v) - u }
    } else if a > 2.0 {
        if b > 1000.0 {
            let n = (a - 1.0).floor() as usize;
            let mut prod = 1.0;
            let mut ar = a;
            for _ in 0..n {
                ar -= 1.0;
                prod *= ar / (1.0 + ar / b);
            }
            (log(prod) - n as f64 * log(b)) + (gamma_reduced(ar) + gamma_difference(ar, b))
        } else {
            let mut prod1 = 1.0;
            let mut ar = a;
            while ar > 2.0 {
                ar -= 1.0;
                let h = ar / b;
                prod1 *= h / (1.0 + h);
            }
            if b < 10.0 {
                let mut prod2 = 1.0;
                let mut br = b;
                while br > 2.0 {
                    br -= 1.0;
                    prod2 *= br / (ar + br);
                }
                log(prod1)
                    + log(prod2)
                    + (gamma_reduced(ar) + (gamma_reduced(br) - gamma_sum(ar, br)))
            } else {
                log(prod1) + gamma_reduced(ar) + gamma_difference(ar, b)
            }
        }
    } else if a < 1.0 {
        if b >= 10.0 {
            gamma_reduced(a) + gamma_difference(a, b)
        } else {
            log(
                crate::math::gamma::gamma_small(a) * crate::math::gamma::gamma_small(b)
                    / crate::math::gamma::gamma_small(a + b),
            )
        }
    } else if b > 2.0 {
        if b < 10.0 {
            let mut prod = 1.0;
            let mut br = b;
            while br > 2.0 {
                br -= 1.0;
                prod *= br / (a + br);
            }
            log(prod) + (gamma_reduced(a) + (gamma_reduced(br) - gamma_sum(a, br)))
        } else {
            gamma_reduced(a) + gamma_difference(a, b)
        }
    } else {
        gamma_reduced(a) + gamma_reduced(b) - gamma_sum(a, b)
    }
}
pub(crate) fn binomial_log(n: u32, k: u32) -> f64 {
    if n == k || k == 0 {
        return 0.0;
    }
    if k == 1 || k == n - 1 {
        return log(n as f64);
    }
    let k = k.min(n - k);
    if n < 67 {
        // Exact integer result; u128 avoids Java's gcd overflow workaround.
        let mut r = 1u128;
        for j in 1..=k {
            r = r * (n - k + j) as u128 / j as u128;
        }
        log(r as f64)
    } else if n < 1030 {
        let mut r = 1.0;
        for i in 1..=k {
            r *= (n - k + i) as f64 / i as f64;
        }
        log((r + 0.5).floor())
    } else {
        let mut r = 0.0;
        for i in n - k + 1..=n {
            r += log(i as f64);
        }
        for i in 2..=k {
            r -= log(i as f64);
        }
        r
    }
}
pub(crate) fn beta_binomial(n: u32, k: u32, a: f64, b: f64) -> f64 {
    if k > n {
        f64::NEG_INFINITY
    } else {
        binomial_log(n, k) + log_beta(k as f64 + a, (n - k) as f64 + b) - log_beta(a, b)
    }
}
