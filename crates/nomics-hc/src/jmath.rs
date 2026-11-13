#![allow(
    clippy::approx_constant,
    clippy::eq_op,
    clippy::excessive_precision,
    clippy::needless_late_init,
    clippy::collapsible_if
)]

//! OpenJDK 17 fdlibm routines used by `StrictMath`.
//!
//! This is a literal, operation-order-preserving port of the JDK sources.
//! In particular, do not replace polynomial expressions with `mul_add`.

const TWO54: f64 = 1.801_439_850_948_198_4e16;
const LN2_HI: f64 = 6.931_471_803_691_238_164_90e-1;
const LN2_LO: f64 = 1.908_214_929_270_587_700_02e-10;

#[inline]
fn hi(x: f64) -> i32 {
    (x.to_bits() >> 32) as u32 as i32
}
#[inline]
fn lo(x: f64) -> i32 {
    x.to_bits() as u32 as i32
}
#[inline]
fn with_hi(x: f64, h: i32) -> f64 {
    f64::from_bits(((h as u32 as u64) << 32) | (x.to_bits() & 0xffff_ffff))
}
#[inline]
fn with_lo(x: f64, l: i32) -> f64 {
    f64::from_bits((x.to_bits() & !0xffff_ffff) | l as u32 as u64)
}

/// OpenJDK native fdlibm `__ieee754_log`.
pub fn log(mut x: f64) -> f64 {
    const L1: f64 = 6.666_666_666_666_735_130e-1;
    const L2: f64 = 3.999_999_999_940_941_908e-1;
    const L3: f64 = 2.857_142_874_366_239_149e-1;
    const L4: f64 = 2.222_219_843_214_978_396e-1;
    const L5: f64 = 1.818_357_216_161_805_012e-1;
    const L6: f64 = 1.531_383_769_920_937_332e-1;
    const L7: f64 = 1.479_819_860_511_658_591e-1;
    let mut hx = hi(x);
    let lx = lo(x) as u32;
    let mut k = 0i32;
    if hx < 0x0010_0000 {
        if ((hx & 0x7fff_ffff) as u32 | lx) == 0 {
            return -TWO54 / 0.0;
        }
        if hx < 0 {
            return (x - x) / 0.0;
        }
        k -= 54;
        x *= TWO54;
        hx = hi(x);
    }
    if hx >= 0x7ff0_0000 {
        return x + x;
    }
    k += (hx >> 20) - 1023;
    hx &= 0x000f_ffff;
    let i = (hx + 0x95f64) & 0x10_0000;
    x = with_hi(x, hx | (i ^ 0x3ff0_0000));
    k += i >> 20;
    let f = x - 1.0;
    if (0x000f_ffff & (2 + hx)) < 3 {
        if f == 0.0 {
            if k == 0 {
                return 0.0;
            }
            let dk = k as f64;
            return dk * LN2_HI + dk * LN2_LO;
        }
        let r = f * f * (0.5 - 0.333_333_333_333_333_33 * f);
        if k == 0 {
            return f - r;
        }
        let dk = k as f64;
        return dk * LN2_HI - ((r - dk * LN2_LO) - f);
    }
    let s = f / (2.0 + f);
    let dk = k as f64;
    let z = s * s;
    let ii = hx - 0x6147a;
    let w = z * z;
    let j = 0x6b851 - hx;
    let t1 = w * (L2 + w * (L4 + w * L6));
    let t2 = z * (L1 + w * (L3 + w * (L5 + w * L7)));
    let r = t2 + t1;
    if (ii | j) > 0 {
        let hfsq = 0.5 * f * f;
        if k == 0 {
            f - (hfsq - s * (hfsq + r))
        } else {
            dk * LN2_HI - ((hfsq - (s * (hfsq + r) + dk * LN2_LO)) - f)
        }
    } else if k == 0 {
        f - s * (f - r)
    } else {
        dk * LN2_HI - ((s * (f - r) - dk * LN2_LO) - f)
    }
}

/// OpenJDK native fdlibm `__ieee754_log10`.
pub fn log10(mut x: f64) -> f64 {
    const IVLN10: f64 = 4.342_944_819_032_518_166_68e-1;
    const HI: f64 = 3.010_299_956_636_117_713_06e-1;
    const LO: f64 = 3.694_239_077_158_930_786_16e-13;
    let mut hx = hi(x);
    let lx = lo(x) as u32;
    let mut k = 0i32;
    if hx < 0x0010_0000 {
        if ((hx & 0x7fff_ffff) as u32 | lx) == 0 {
            return -TWO54 / 0.0;
        }
        if hx < 0 {
            return (x - x) / 0.0;
        }
        k -= 54;
        x *= TWO54;
        hx = hi(x);
    }
    if hx >= 0x7ff0_0000 {
        return x + x;
    }
    k += (hx >> 20) - 1023;
    let i = ((k as u32 & 0x8000_0000) >> 31) as i32;
    hx = (hx & 0x000f_ffff) | ((0x3ff - i) << 20);
    let y = (k + i) as f64;
    x = with_hi(x, hx);
    let z = y * LO + IVLN10 * log(x);
    z + y * HI
}

/// OpenJDK native fdlibm `log1p`.
pub fn log1p(x: f64) -> f64 {
    const L1: f64 = 6.666_666_666_666_735_130e-1;
    const L2: f64 = 3.999_999_999_940_941_908e-1;
    const L3: f64 = 2.857_142_874_366_239_149e-1;
    const L4: f64 = 2.222_219_843_214_978_396e-1;
    const L5: f64 = 1.818_357_216_161_805_012e-1;
    const L6: f64 = 1.531_383_769_920_937_332e-1;
    const L7: f64 = 1.479_819_860_511_658_591e-1;
    let hx = hi(x);
    let ax = hx & 0x7fff_ffff;
    let mut k = 1i32;
    let mut f = 0.0;
    let mut c = 0.0;
    let mut hu = 0i32;
    if hx < 0x3fda_827a {
        if ax >= 0x3ff0_0000 {
            if x == -1.0 && hx == 0xbff0_0000u32 as i32 {
                return -TWO54 / 0.0;
            }
            return (x - x) / (x - x);
        }
        if ax < 0x3e20_0000 {
            if TWO54 + x > 0.0 && ax < 0x3c90_0000 {
                return x;
            }
            return x - x * x * 0.5;
        }
        if hx > 0 || hx <= 0xbfd2_bec3u32 as i32 {
            k = 0;
            f = x;
            hu = 1;
        }
    }
    if hx >= 0x7ff0_0000 {
        return x + x;
    }
    if k != 0 {
        let mut u;
        if hx < 0x4340_0000 {
            u = 1.0 + x;
            hu = hi(u);
            k = (hu >> 20) - 1023;
            c = if k > 0 { 1.0 - (u - x) } else { x - (u - 1.0) };
            c /= u;
        } else {
            u = x;
            hu = hi(u);
            k = (hu >> 20) - 1023;
            c = 0.0;
        }
        hu &= 0x000f_ffff;
        if hu < 0x6a09e {
            u = with_hi(u, hu | 0x3ff0_0000);
        } else {
            k += 1;
            u = with_hi(u, hu | 0x3fe0_0000);
            hu = (0x0010_0000 - hu) >> 2;
        }
        f = u - 1.0;
    }
    let hfsq = 0.5 * f * f;
    if hu == 0 {
        if f == 0.0 {
            if k == 0 {
                return 0.0;
            }
            c += (k as f64) * LN2_LO;
            return (k as f64) * LN2_HI + c;
        }
        let r = hfsq * (1.0 - 0.666_666_666_666_666_66 * f);
        return if k == 0 {
            f - r
        } else {
            (k as f64) * LN2_HI - ((r - ((k as f64) * LN2_LO + c)) - f)
        };
    }
    let s = f / (2.0 + f);
    let z = s * s;
    let r = z * (L1 + z * (L2 + z * (L3 + z * (L4 + z * (L5 + z * (L6 + z * L7))))));
    if k == 0 {
        f - (hfsq - s * (hfsq + r))
    } else {
        (k as f64) * LN2_HI - ((hfsq - (s * (hfsq + r) + ((k as f64) * LN2_LO + c))) - f)
    }
}

/// OpenJDK native fdlibm `expm1`.
pub fn expm1(mut x: f64) -> f64 {
    const O: f64 = 7.097_827_128_933_839_730_96e2;
    const Q1: f64 = -3.333_333_333_333_313_164_28e-2;
    const Q2: f64 = 1.587_301_587_254_814_601_65e-3;
    const Q3: f64 = -7.936_507_578_674_879_424_73e-5;
    const Q4: f64 = 4.008_217_827_329_362_395_52e-6;
    const Q5: f64 = -2.010_992_181_836_243_713_26e-7;
    let mut hx = hi(x) as u32;
    let xsb = hx & 0x8000_0000;
    hx &= 0x7fff_ffff;
    let mut k = 0i32;
    let mut c = 0.0;
    if hx >= 0x4043_687a {
        if hx >= 0x4086_2e42 {
            if hx >= 0x7ff0_0000 {
                return if ((hx & 0xfffff) | lo(x) as u32) != 0 {
                    x + x
                } else if xsb == 0 {
                    x
                } else {
                    -1.0
                };
            }
            if x > O {
                return 1.0e300 * 1.0e300;
            }
        }
        if xsb != 0 {
            return 1.0e-300 - 1.0;
        }
    }
    if hx > 0x3fd6_2e42 {
        let (h, l);
        if hx < 0x3ff0_a2b2 {
            if xsb == 0 {
                h = x - LN2_HI;
                l = LN2_LO;
                k = 1;
            } else {
                h = x + LN2_HI;
                l = -LN2_LO;
                k = -1;
            }
        } else {
            k = (1.442_695_040_888_963_387 * x + if xsb == 0 { 0.5 } else { -0.5 }) as i32;
            let t = k as f64;
            h = x - t * LN2_HI;
            l = t * LN2_LO;
        }
        x = h - l;
        c = (h - x) - l;
    } else if hx < 0x3c90_0000 {
        let t = 1.0e300 + x;
        return x - (t - (1.0e300 + x));
    }
    let hfx = 0.5 * x;
    let hxs = x * hfx;
    let r1 = 1.0 + hxs * (Q1 + hxs * (Q2 + hxs * (Q3 + hxs * (Q4 + hxs * Q5))));
    let t = 3.0 - r1 * hfx;
    let mut e = hxs * ((r1 - t) / (6.0 - x * t));
    if k == 0 {
        return x - (x * e - hxs);
    }
    e = x * (e - c) - c;
    e -= hxs;
    if k == -1 {
        return 0.5 * (x - e) - 0.5;
    }
    if k == 1 {
        return if x < -0.25 {
            -2.0 * (e - (x + 0.5))
        } else {
            1.0 + 2.0 * (x - e)
        };
    }
    if k <= -2 || k > 56 {
        let y = with_hi(1.0 - (e - x), hi(1.0 - (e - x)).wrapping_add(k << 20));
        return y - 1.0;
    }
    let mut t = 1.0;
    let y;
    if k < 20 {
        t = with_hi(t, 0x3ff0_0000 - (0x20_0000 >> k));
        let q = t - (e - x);
        y = with_hi(q, hi(q).wrapping_add(k << 20));
    } else {
        t = with_hi(t, (0x3ff - k) << 20);
        let mut q = x - (e + t);
        q += 1.0;
        y = with_hi(q, hi(q).wrapping_add(k << 20));
    }
    y
}

/// Java `FdLibm.Exp.compute`.
pub fn exp(mut x: f64) -> f64 {
    const O: f64 = f64::from_bits(0x4086_2e42_fefa_39ef);
    const U: f64 = f64::from_bits(0xc087_4910_d52d_3051);
    const P1: f64 = f64::from_bits(0x3fc5_5555_5555_553e);
    const P2: f64 = f64::from_bits(0xbf66_c16c_16be_bd93);
    const P3: f64 = f64::from_bits(0x3f11_566a_af25_de2c);
    const P4: f64 = f64::from_bits(0xbebb_bd41_c5d2_6bf1);
    const P5: f64 = f64::from_bits(0x3e66_3769_72be_a4d0);
    let mut hx = hi(x) as u32;
    let xsb = (hx >> 31) & 1;
    hx &= 0x7fff_ffff;
    let mut k = 0i32;
    let mut low = 0.0;
    let mut h = 0.0;
    if hx >= 0x4086_2e42 {
        if hx >= 0x7ff0_0000 {
            if ((hx & 0xfffff) | lo(x) as u32) != 0 {
                return x + x;
            }
            return if xsb == 0 { x } else { 0.0 };
        }
        if x > O {
            return 1.0e300 * 1.0e300;
        }
        if x < U {
            return f64::from_bits(0x0170_0000_0000_0000) * f64::from_bits(0x0170_0000_0000_0000);
        }
    }
    if hx > 0x3fd6_2e42 {
        if hx < 0x3ff0_a2b2 {
            h = x - if xsb == 0 { LN2_HI } else { -LN2_HI };
            low = if xsb == 0 { LN2_LO } else { -LN2_LO };
            k = 1 - (xsb as i32) - (xsb as i32);
        } else {
            k = (f64::from_bits(0x3ff7_1547_652b_82fe) * x + if xsb == 0 { 0.5 } else { -0.5 })
                as i32;
            let t = k as f64;
            h = x - t * LN2_HI;
            low = t * LN2_LO;
        }
        x = h - low;
    } else if hx < 0x3e30_0000 {
        if 1.0e300 + x > 1.0 {
            return 1.0 + x;
        }
    }
    let t = x * x;
    let c = x - t * (P1 + t * (P2 + t * (P3 + t * (P4 + t * P5))));
    let mut y = if k == 0 {
        1.0 - ((x * c) / (c - 2.0) - x)
    } else {
        1.0 - ((low - (x * c) / (2.0 - c)) - h)
    };
    if k == 0 {
        return y;
    }
    if k >= -1021 {
        y = with_hi(y, hi(y).wrapping_add(k << 20));
        y
    } else {
        y = with_hi(y, hi(y).wrapping_add((k + 1000) << 20));
        y * f64::from_bits(0x0170_0000_0000_0000)
    }
}

/// Java `FdLibm.Pow.compute` (not the platform `pow`).
pub fn pow(x: f64, y: f64) -> f64 {
    if y == 0.0 {
        return 1.0;
    }
    if x.is_nan() || y.is_nan() {
        return x + y;
    }
    let ya = y.abs();
    let mut xa = x.abs();
    if y == 2.0 {
        return x * x;
    }
    if y == 0.5 && x >= -f64::MAX {
        return (x + 0.0).sqrt();
    }
    if ya == 1.0 {
        return if y == 1.0 { x } else { 1.0 / x };
    }
    if ya == f64::INFINITY {
        if xa == 1.0 {
            return y - y;
        }
        if xa > 1.0 {
            return if y >= 0.0 { y } else { 0.0 };
        }
        return if y < 0.0 { -y } else { 0.0 };
    }
    let hx = hi(x);
    let mut ix = hx & 0x7fff_ffff;
    let mut yi = 0i32;
    if hx < 0 {
        if ya >= 9_007_199_254_740_992.0 {
            yi = 2;
        } else if ya >= 1.0 {
            let l = ya as i64;
            if l as f64 == ya {
                yi = 2 - (l as i32 & 1);
            }
        }
    }
    if xa == 0.0 || xa == f64::INFINITY || xa == 1.0 {
        let mut z = xa;
        if y < 0.0 {
            z = 1.0 / z;
        }
        if hx < 0 {
            if ((ix - 0x3ff0_0000) | yi) == 0 {
                z = (z - z) / (z - z);
            } else if yi == 1 {
                z = -z;
            }
        }
        return z;
    }
    let mut n = (hx >> 31) + 1;
    if (n | yi) == 0 {
        return (x - x) / (x - x);
    }
    let mut sign = 1.0;
    if (n | (yi - 1)) == 0 {
        sign = -1.0;
    }
    let (mut ph, mut pl, t1, t2): (f64, f64, f64, f64);
    if ya > f64::from_bits(0x41e0_0000_ffff_ffff) {
        const IL2: f64 = f64::from_bits(0x3ff7_1547_652b_82fe);
        const IL2H: f64 = f64::from_bits(0x3ff7_1547_6000_0000);
        const IL2L: f64 = f64::from_bits(0x3e54_ae0b_f85d_df44);
        if xa < f64::from_bits(0x3fef_ffff_0000_0000) {
            return if y < 0.0 {
                sign * f64::INFINITY
            } else {
                sign * 0.0
            };
        }
        if xa > f64::from_bits(0x3ff0_0000_ffff_ffff) {
            return if y > 0.0 {
                sign * f64::INFINITY
            } else {
                sign * 0.0
            };
        }
        let t = xa - 1.0;
        let w = (t * t) * (0.5 - t * (0.333_333_333_333_333_3 - t * 0.25));
        let u = IL2H * t;
        let v = t * IL2L - w * IL2;
        t1 = with_lo(u + v, 0);
        t2 = v - (t1 - u);
    } else {
        const CP: f64 = f64::from_bits(0x3fee_c709_dc3a_03fd);
        const CPH: f64 = f64::from_bits(0x3fee_c709_e000_0000);
        const CPL: f64 = f64::from_bits(0xbe3e_2fe0_145b_01f5);
        const BP: [f64; 2] = [1.0, 1.5];
        const DPH: [f64; 2] = [0.0, f64::from_bits(0x3fe2_b803_4000_0000)];
        const DPL: [f64; 2] = [0.0, f64::from_bits(0x3e4c_fdeb_43cf_d006)];
        const L1: f64 = f64::from_bits(0x3fe3_3333_3333_3303);
        const L2: f64 = f64::from_bits(0x3fdb_6db6_db6f_abff);
        const L3: f64 = f64::from_bits(0x3fd5_5555_518f_264d);
        const L4: f64 = f64::from_bits(0x3fd1_7460_a91d_4101);
        const L5: f64 = f64::from_bits(0x3fcd_864a_93c9_db65);
        const L6: f64 = f64::from_bits(0x3fca_7e28_4a45_4eef);
        n = 0;
        if ix < 0x0010_0000 {
            xa *= 9_007_199_254_740_992.0;
            n -= 53;
            ix = hi(xa);
        }
        n += (ix >> 20) - 0x3ff;
        let j0 = ix & 0x000f_ffff;
        ix = j0 | 0x3ff0_0000;
        let k: usize;
        if j0 <= 0x3988e {
            k = 0;
        } else if j0 < 0xbb67a {
            k = 1;
        } else {
            k = 0;
            n += 1;
            ix -= 0x0010_0000;
        }
        xa = with_hi(xa, ix);
        let u = xa - BP[k];
        let v = 1.0 / (xa + BP[k]);
        let ss = u * v;
        let sh = with_lo(ss, 0);
        let th = with_hi(
            0.0,
            ((ix >> 1) | 0x2000_0000) + 0x0008_0000 + ((k as i32) << 18),
        );
        let tl = xa - (th - BP[k]);
        let sl = v * ((u - sh * th) - sh * tl);
        let s2 = ss * ss;
        let mut r = s2 * s2 * (L1 + s2 * (L2 + s2 * (L3 + s2 * (L4 + s2 * (L5 + s2 * L6)))));
        r += sl * (sh + ss);
        let sq = sh * sh;
        let th2 = with_lo(3.0 + sq + r, 0);
        let tl2 = r - ((th2 - 3.0) - sq);
        let u2 = sh * th2;
        let v2 = sl * th2 + tl2 * ss;
        ph = with_lo(u2 + v2, 0);
        pl = v2 - (ph - u2);
        let zh = CPH * ph;
        let zl = CPL * ph + pl * CP + DPL[k];
        let q = n as f64;
        t1 = with_lo((zh + zl) + DPH[k] + q, 0);
        t2 = zl - (((t1 - q) - DPH[k]) - zh);
    }
    let y1 = with_lo(y, 0);
    pl = (y - y1) * t1 + y * t2;
    ph = y1 * t1;
    let mut z = pl + ph;
    let mut j = hi(z);
    let i = lo(z);
    if j >= 0x4090_0000 {
        if ((j - 0x4090_0000) | i) != 0 {
            return sign * f64::INFINITY;
        }
        if pl + 8.008_566_259_537_294_437_2e-17 > z - ph {
            return sign * f64::INFINITY;
        }
    } else if (j & 0x7fff_ffff) >= 0x4090_cc00 {
        if ((j - 0xc090_cc00u32 as i32) | i) != 0 {
            return sign * 0.0;
        }
        if pl <= z - ph {
            return sign * 0.0;
        }
    }
    const P1: f64 = f64::from_bits(0x3fc5_5555_5555_553e);
    const P2: f64 = f64::from_bits(0xbf66_c16c_16be_bd93);
    const P3: f64 = f64::from_bits(0x3f11_566a_af25_de2c);
    const P4: f64 = f64::from_bits(0xbebb_bd41_c5d2_6bf1);
    const P5: f64 = f64::from_bits(0x3e66_3769_72be_a4d0);
    let ii = j & 0x7fff_ffff;
    let mut kk = (ii >> 20) - 0x3ff;
    n = 0;
    if ii > 0x3fe0_0000 {
        n = j.wrapping_add(0x0010_0000 >> (kk + 1));
        kk = ((n & 0x7fff_ffff) >> 20) - 0x3ff;
        let tt = with_hi(0.0, n & !(0x000f_ffff >> kk));
        n = ((n & 0x000f_ffff) | 0x0010_0000) >> (20 - kk);
        if j < 0 {
            n = -n;
        }
        ph -= tt;
    }
    let t = with_lo(pl + ph, 0);
    let u = t * f64::from_bits(0x3fe6_2e43_0000_0000);
    let v = (pl - (t - ph)) * f64::from_bits(0x3fe6_2e42_fefa_39ef)
        + t * f64::from_bits(0xbe20_5c61_0ca8_6c39);
    z = u + v;
    let w = v - (z - u);
    let ts = z * z;
    let q = z - ts * (P1 + ts * (P2 + ts * (P3 + ts * (P4 + ts * P5))));
    let r = z * q / (q - 2.0) - (w + z * w);
    z = 1.0 - (r - z);
    j = hi(z).wrapping_add(n << 20);
    if (j >> 20) <= 0 {
        z *= 2.0_f64.powi(n);
    } else {
        z = with_hi(z, hi(z).wrapping_add(n << 20));
    }
    sign * z
}
