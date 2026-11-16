//! Square-matrix path of Commons Math 3.5 linear/SingularValueDecomposition.java.
//! Bidiagonalization and four QR cases retain upstream operation and loop order.
//! Apache-2.0; algorithm attribution: ASF, JAMA.
#![allow(clippy::needless_range_loop, clippy::assign_op_pattern)]
use super::fastmath::hypot;
pub fn singular_value_decomposition(mut a: Vec<Vec<f64>>) -> (Vec<f64>, Vec<Vec<f64>>) {
    let n = a.len();
    assert!(n > 0 && a.iter().all(|r| r.len() == n));
    let mut s = vec![0.0; n];
    let mut u = vec![vec![0.0; n]; n];
    let mut v = u.clone();
    let mut e = vec![0.0; n];
    let mut work = vec![0.0; n];
    let nct = n - 1;
    let nrt = n.saturating_sub(2);
    for k in 0..nct.max(nrt) {
        if k < nct {
            s[k] = 0.0;
            for i in k..n {
                s[k] = hypot(s[k], a[i][k]);
            }
            if s[k] != 0.0 {
                if a[k][k] < 0.0 {
                    s[k] = -s[k];
                }
                for i in k..n {
                    a[i][k] /= s[k];
                }
                a[k][k] += 1.0;
            }
            s[k] = -s[k];
        }
        for j in k + 1..n {
            if k < nct && s[k] != 0.0 {
                let mut t = 0.0;
                for i in k..n {
                    t += a[i][k] * a[i][j];
                }
                t = -t / a[k][k];
                for i in k..n {
                    a[i][j] += t * a[i][k];
                }
            }
            e[j] = a[k][j];
        }
        if k < nct {
            for i in k..n {
                u[i][k] = a[i][k];
            }
        }
        if k < nrt {
            e[k] = 0.0;
            for i in k + 1..n {
                e[k] = hypot(e[k], e[i]);
            }
            if e[k] != 0.0 {
                if e[k + 1] < 0.0 {
                    e[k] = -e[k];
                }
                for i in k + 1..n {
                    e[i] /= e[k];
                }
                e[k + 1] += 1.0;
            }
            e[k] = -e[k];
            if k + 1 < n && e[k] != 0.0 {
                for i in k + 1..n {
                    work[i] = 0.0;
                }
                for j in k + 1..n {
                    for i in k + 1..n {
                        work[i] += e[j] * a[i][j];
                    }
                }
                for j in k + 1..n {
                    let t = -e[j] / e[k + 1];
                    for i in k + 1..n {
                        a[i][j] += t * work[i];
                    }
                }
            }
            for i in k + 1..n {
                v[i][k] = e[i];
            }
        }
    }
    let mut p = n;
    s[nct] = a[nct][nct];
    if nrt + 1 < p {
        e[nrt] = a[nrt][p - 1];
    }
    e[p - 1] = 0.0;
    for j in nct..n {
        for i in 0..n {
            u[i][j] = 0.0;
        }
        u[j][j] = 1.0;
    }
    for k in (0..nct).rev() {
        if s[k] != 0.0 {
            for j in k + 1..n {
                let mut t = 0.0;
                for i in k..n {
                    t += u[i][k] * u[i][j];
                }
                t = -t / u[k][k];
                for i in k..n {
                    u[i][j] += t * u[i][k];
                }
            }
            for i in k..n {
                u[i][k] = -u[i][k];
            }
            u[k][k] = 1.0 + u[k][k];
            // Upstream deliberately uses i < k - 1 (not i < k).
            for i in 0..k.saturating_sub(1) {
                u[i][k] = 0.0;
            }
        } else {
            for i in 0..n {
                u[i][k] = 0.0;
            }
            u[k][k] = 1.0;
        }
    }
    for k in (0..n).rev() {
        if k < nrt && e[k] != 0.0 {
            for j in k + 1..n {
                let mut t = 0.0;
                for i in k + 1..n {
                    t += v[i][k] * v[i][j];
                }
                t = -t / v[k + 1][k];
                for i in k + 1..n {
                    v[i][j] += t * v[i][k];
                }
            }
        }
        for i in 0..n {
            v[i][k] = 0.0;
        }
        v[k][k] = 1.0;
    }
    let pp = n - 1;
    while p > 0 {
        let mut k = p as isize - 2;
        while k >= 0 {
            let i = k as usize;
            let threshold =
                f64::from_bits(57u64 << 52) + f64::EPSILON * (s[i].abs() + s[i + 1].abs());
            if e[i].abs().partial_cmp(&threshold) != Some(std::cmp::Ordering::Greater) {
                e[i] = 0.0;
                break;
            }
            k -= 1;
        }
        let case;
        if k == p as isize - 2 {
            case = 4;
        } else {
            let mut ks = p as isize - 1;
            while ks >= k {
                if ks == k {
                    break;
                }
                let i = ks as usize;
                let t = (if i != p { e[i].abs() } else { 0.0 })
                    + (if ks != k + 1 { e[i - 1].abs() } else { 0.0 });
                if s[i].abs() <= f64::from_bits(57u64 << 52) + f64::EPSILON * t {
                    s[i] = 0.0;
                    break;
                }
                ks -= 1;
            }
            if ks == k {
                case = 3;
            } else if ks == p as isize - 1 {
                case = 1;
            } else {
                case = 2;
                k = ks;
            }
        }
        let mut k = (k + 1) as usize;
        match case {
            1 => {
                let mut f = e[p - 2];
                e[p - 2] = 0.0;
                for j in (k..=p - 2).rev() {
                    let mut t = hypot(s[j], f);
                    let cs = s[j] / t;
                    let sn = f / t;
                    s[j] = t;
                    if j != k {
                        f = -sn * e[j - 1];
                        e[j - 1] = cs * e[j - 1];
                    }
                    for i in 0..n {
                        t = cs * v[i][j] + sn * v[i][p - 1];
                        v[i][p - 1] = -sn * v[i][j] + cs * v[i][p - 1];
                        v[i][j] = t;
                    }
                }
            }
            2 => {
                let mut f = e[k - 1];
                e[k - 1] = 0.0;
                for j in k..p {
                    let mut t = hypot(s[j], f);
                    let cs = s[j] / t;
                    let sn = f / t;
                    s[j] = t;
                    f = -sn * e[j];
                    e[j] = cs * e[j];
                    for i in 0..n {
                        t = cs * u[i][j] + sn * u[i][k - 1];
                        u[i][k - 1] = -sn * u[i][j] + cs * u[i][k - 1];
                        u[i][j] = t;
                    }
                }
            }
            3 => {
                let scale = s[p - 1]
                    .abs()
                    .max(s[p - 2].abs())
                    .max(e[p - 2].abs())
                    .max(s[k].abs())
                    .max(e[k].abs());
                let sp = s[p - 1] / scale;
                let spm1 = s[p - 2] / scale;
                let epm1 = e[p - 2] / scale;
                let sk = s[k] / scale;
                let ek = e[k] / scale;
                let b = ((spm1 + sp) * (spm1 - sp) + epm1 * epm1) / 2.0;
                let c = (sp * epm1) * (sp * epm1);
                let mut shift = 0.0;
                if b != 0.0 || c != 0.0 {
                    shift = (b * b + c).sqrt();
                    if b < 0.0 {
                        shift = -shift;
                    }
                    shift = c / (b + shift);
                }
                let mut f = (sk + sp) * (sk - sp) + shift;
                let mut g = sk * ek;
                for j in k..p - 1 {
                    let mut t = hypot(f, g);
                    let mut cs = f / t;
                    let mut sn = g / t;
                    if j != k {
                        e[j - 1] = t;
                    }
                    f = cs * s[j] + sn * e[j];
                    e[j] = cs * e[j] - sn * s[j];
                    g = sn * s[j + 1];
                    s[j + 1] = cs * s[j + 1];
                    for i in 0..n {
                        t = cs * v[i][j] + sn * v[i][j + 1];
                        v[i][j + 1] = -sn * v[i][j] + cs * v[i][j + 1];
                        v[i][j] = t;
                    }
                    t = hypot(f, g);
                    cs = f / t;
                    sn = g / t;
                    s[j] = t;
                    f = cs * e[j] + sn * s[j + 1];
                    s[j + 1] = -sn * e[j] + cs * s[j + 1];
                    g = sn * e[j + 1];
                    e[j + 1] = cs * e[j + 1];
                    if j < n - 1 {
                        for i in 0..n {
                            t = cs * u[i][j] + sn * u[i][j + 1];
                            u[i][j + 1] = -sn * u[i][j] + cs * u[i][j + 1];
                            u[i][j] = t;
                        }
                    }
                }
                e[p - 2] = f;
            }
            _ => {
                if s[k] <= 0.0 {
                    s[k] = if s[k] < 0.0 { -s[k] } else { 0.0 };
                    for i in 0..=pp {
                        v[i][k] = -v[i][k];
                    }
                }
                while k < pp {
                    if s[k] >= s[k + 1] {
                        break;
                    }
                    s.swap(k, k + 1);
                    if k < n - 1 {
                        for i in 0..n {
                            v[i].swap(k, k + 1);
                        }
                        for i in 0..n {
                            u[i].swap(k, k + 1);
                        }
                    }
                    k += 1;
                }
                p -= 1;
            }
        }
    }
    (s, u)
}
