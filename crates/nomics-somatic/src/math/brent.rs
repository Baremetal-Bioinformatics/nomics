//! Commons Math 3.5 optim/univariate/BrentOptimizer.java (maximizing, no checker).
use std::io;
fn best(a: (f64, f64), b: (f64, f64)) -> (f64, f64) {
    if a.1 >= b.1 { a } else { b }
}
fn equal(a: f64, b: f64) -> bool {
    a == b
        || (a.is_sign_positive() == b.is_sign_positive() && a.to_bits().abs_diff(b.to_bits()) <= 1)
}
pub fn brent_optimizer(
    f: impl Fn(f64) -> f64,
    mut a: f64,
    mut b: f64,
    initial: f64,
    rel: f64,
    abs: f64,
    max_eval: usize,
) -> io::Result<(f64, f64)> {
    let (mut x, mut v, mut w) = (initial, initial, initial);
    let (mut d, mut e) = (0.0_f64, 0.0_f64);
    let mut fx = -f(x);
    let (mut fv, mut fw) = (fx, fx);
    let mut current = (x, -fx);
    let mut optimum = current;
    let mut evaluations = 1;
    loop {
        let m = 0.5 * (a + b);
        let tol1 = rel * x.abs() + abs;
        let tol2 = 2.0 * tol1;
        if (x - m).abs() <= tol2 - 0.5 * (b - a) {
            return Ok(optimum);
        }
        if e.abs() > tol1 {
            let mut r = (x - w) * (fx - fv);
            let mut q = (x - v) * (fx - fw);
            let mut p = (x - v) * q - (x - w) * r;
            q = 2.0 * (q - r);
            if q > 0.0 {
                p = -p;
            } else {
                q = -q;
            }
            r = e;
            e = d;
            if p > q * (a - x) && p < q * (b - x) && p.abs() < (0.5 * q * r).abs() {
                d = p / q;
                let u = x + d;
                if u - a < tol2 || b - u < tol2 {
                    d = if x <= m { tol1 } else { -tol1 };
                }
            } else {
                e = if x < m { b - x } else { a - x };
                d = 0.5 * (3.0 - 5.0_f64.sqrt()) * e;
            }
        } else {
            e = if x < m { b - x } else { a - x };
            d = 0.5 * (3.0 - 5.0_f64.sqrt()) * e;
        }
        let u = if d.abs() < tol1 {
            if d >= 0.0 { x + tol1 } else { x - tol1 }
        } else {
            x + d
        };
        evaluations += 1;
        if evaluations > max_eval {
            return Err(io::Error::other(
                "BrentOptimizer maximum evaluations exceeded",
            ));
        }
        let fu = -f(u);
        let previous = current;
        current = (u, -fu);
        optimum = best(optimum, best(previous, current));
        if fu <= fx {
            if u < x {
                b = x;
            } else {
                a = x;
            }
            v = w;
            fv = fw;
            w = x;
            fw = fx;
            x = u;
            fx = fu;
        } else {
            if u < x {
                a = u;
            } else {
                b = u;
            }
            if fu <= fw || equal(w, x) {
                v = w;
                fv = fw;
                w = u;
                fw = fu;
            } else if fu <= fv || equal(v, x) || equal(v, w) {
                v = u;
                fv = fu;
            }
        }
    }
}
