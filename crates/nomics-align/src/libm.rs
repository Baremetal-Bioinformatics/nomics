//! `erfc` exactly as glibc 2.39 computes it on aarch64
//! (`sysdeps/ieee754/dbl-64/s_erf.c`, fdlibm-derived).
//!
//! glibc's aarch64 libm is built with floating-point contraction, so most
//! `a + b * c` in the C source are single fused multiply-adds (39 `fmadd`/
//! `fmsub` in the object code), mirrored here with `mul_add`; the two that
//! glibc's build left unfused are marked.
//! `exp` is `glibc_math::exp`, a bit-exact port of glibc 2.39's.
//!
//! Verified bit-for-bit against the Spark's libm by
//! `examples/erfc_vs_libm.rs`.

#![allow(clippy::excessive_precision)] // constants copied verbatim from glibc

const ERX: f64 = 8.45062911510467529297e-01;
const PP: [f64; 5] = [
    1.28379167095512558561e-01,
    -3.25042107247001499370e-01,
    -2.84817495755985104766e-02,
    -5.77027029648944159157e-03,
    -2.37630166566501626084e-05,
];
const QQ: [f64; 6] = [
    0.0,
    3.97917223959155352819e-01,
    6.50222499887672944485e-02,
    5.08130628187576562776e-03,
    1.32494738004321644526e-04,
    -3.96022827877536812320e-06,
];
const PA: [f64; 7] = [
    -2.36211856075265944077e-03,
    4.14856118683748331666e-01,
    -3.72207876035701323847e-01,
    3.18346619901161753674e-01,
    -1.10894694282396677476e-01,
    3.54783043256182359371e-02,
    -2.16637559486879084300e-03,
];
const QA: [f64; 7] = [
    0.0,
    1.06420880400844228286e-01,
    5.40397917702171048937e-01,
    7.18286544141962662868e-02,
    1.26171219808761642112e-01,
    1.36370839120290507362e-02,
    1.19844998467991074170e-02,
];
const RA: [f64; 8] = [
    -9.86494403484714822705e-03,
    -6.93858572707181764372e-01,
    -1.05586262253232909814e+01,
    -6.23753324503260060396e+01,
    -1.62396669462573470355e+02,
    -1.84605092906711035994e+02,
    -8.12874355063065934246e+01,
    -9.81432934416914548592e+00,
];
const SA: [f64; 9] = [
    0.0,
    1.96512716674392571292e+01,
    1.37657754143519042600e+02,
    4.34565877475229228821e+02,
    6.45387271733267880336e+02,
    4.29008140027567833386e+02,
    1.08635005541779435134e+02,
    6.57024977031928170135e+00,
    -6.04244152148580987438e-02,
];
const RB: [f64; 7] = [
    -9.86494292470009928597e-03,
    -7.99283237680523006574e-01,
    -1.77579549177547519889e+01,
    -1.60636384855821916062e+02,
    -6.37566443368389627722e+02,
    -1.02509513161107724954e+03,
    -4.83519191608651397019e+02,
];
const SB: [f64; 8] = [
    0.0,
    3.03380607434824582924e+01,
    3.25792512996573918826e+02,
    1.53672958608443695994e+03,
    3.19985821950859553908e+03,
    2.55305040643316442583e+03,
    4.74528541206955367215e+02,
    -2.24409524465858183362e+01,
];
const TINY: f64 = 1e-300;

/// `a + b * c`, contracted.
#[inline]
fn fma(a: f64, b: f64, c: f64) -> f64 {
    b.mul_add(c, a)
}

/// glibc `__erfc`.
pub fn erfc(x: f64) -> f64 {
    let hx = (x.to_bits() >> 32) as u32 as i32;
    let ix = hx & 0x7fff_ffff;
    if ix >= 0x7ff0_0000 {
        return (((hx as u32) >> 31) << 1) as f64 + 1.0 / x;
    }
    if ix < 0x3feb_0000 {
        // |x| < 0.84375
        if ix < 0x3c70_0000 {
            return 1.0 - x;
        }
        let z = x * x;
        let r1 = fma(PP[0], z, PP[1]);
        let z2 = z * z;
        let r2 = fma(PP[2], z, PP[3]);
        let z4 = z2 * z2;
        let s1 = fma(1.0, z, QQ[1]);
        let s2 = fma(QQ[2], z, QQ[3]);
        let s3 = fma(QQ[4], z, QQ[5]);
        let r = fma(fma(r1, z2, r2), z4, PP[4]);
        let s = fma(fma(s1, z2, s2), z4, s3);
        let y = r / s;
        // Not contracted in glibc's object code: x * y is a separate fmul.
        let xy = x * y;
        return if hx < 0x3fd0_0000 {
            1.0 - (x + xy)
        } else {
            0.5 - (xy + (x - 0.5))
        };
    }
    if ix < 0x3ff4_0000 {
        // 0.84375 <= |x| < 1.25
        let s = x.abs() - 1.0;
        let p1 = fma(PA[0], s, PA[1]);
        let s2 = s * s;
        let q1 = fma(1.0, s, QA[1]);
        let s4 = s2 * s2;
        let p2 = fma(PA[2], s, PA[3]);
        let s6 = s4 * s2;
        let q2 = fma(QA[2], s, QA[3]);
        let p3 = fma(PA[4], s, PA[5]);
        let q3 = fma(QA[4], s, QA[5]);
        let p = fma(fma(fma(p1, s2, p2), s4, p3), s6, PA[6]);
        let q = fma(fma(fma(q1, s2, q2), s4, q3), s6, QA[6]);
        return if hx >= 0 {
            (1.0 - ERX) - p / q
        } else {
            1.0 + (ERX + p / q)
        };
    }
    if ix < 0x403c_0000 {
        // |x| < 28
        let x = x.abs();
        let s = 1.0 / (x * x);
        let (r, ss);
        if ix < 0x4006_DB6D {
            let r1 = fma(RA[0], s, RA[1]);
            let s2 = s * s;
            let s1 = fma(1.0, s, SA[1]);
            let s4 = s2 * s2;
            let r2 = fma(RA[2], s, RA[3]);
            let s6 = s4 * s2;
            let ss2 = fma(SA[2], s, SA[3]);
            let s8 = s4 * s4;
            let r3 = fma(RA[4], s, RA[5]);
            let ss3 = fma(SA[4], s, SA[5]);
            let r4 = fma(RA[6], s, RA[7]);
            let ss4 = fma(SA[6], s, SA[7]);
            r = fma(fma(fma(r1, s2, r2), s4, r3), s6, r4);
            ss = fma(fma(fma(fma(s1, s2, ss2), s4, ss3), s6, ss4), s8, SA[8]);
        } else {
            if hx < 0 && ix >= 0x4018_0000 {
                return 2.0 - TINY;
            }
            let r1 = fma(RB[0], s, RB[1]);
            let s2 = s * s;
            let s1 = fma(1.0, s, SB[1]);
            let s4 = s2 * s2;
            let r2 = fma(RB[2], s, RB[3]);
            let s6 = s4 * s2;
            let ss2 = fma(SB[2], s, SB[3]);
            let r3 = fma(RB[4], s, RB[5]);
            let ss3 = fma(SB[4], s, SB[5]);
            let ss4 = fma(SB[6], s, SB[7]);
            r = fma(fma(fma(r1, s2, r2), s4, r3), s6, RB[6]);
            ss = fma(fma(fma(s1, s2, ss2), s4, ss3), s6, ss4);
        }
        let z = f64::from_bits(x.to_bits() & 0xffff_ffff_0000_0000);
        let e = crate::glibc_math::exp((-z).mul_add(z, -0.5625))
            * crate::glibc_math::exp((z - x).mul_add(z + x, r / ss));
        return if hx > 0 { e / x } else { 2.0 - e / x };
    }
    if hx > 0 { TINY * TINY } else { 2.0 - TINY }
}
