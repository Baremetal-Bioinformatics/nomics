//! OpenJDK 17 fdlibm e_log.c (__ieee754_log), used by Java Math.log on
//! the reference path. Preserve the argument reduction and polynomial ordering.
//! Source: openjdk/jdk17u, jdk-17.0.17-ga, src/java.base/share/native/libfdlibm/e_log.c.
//! Copyright (c) 1998, 2001 Oracle and/or its affiliates; GPL-2.0 with Classpath exception.
#![allow(clippy::excessive_precision)]
pub fn log(mut x: f64) -> f64 {
    const HI: f64 = 6.93147180369123816490e-01;
    const LO: f64 = 1.90821492927058770002e-10;
    let mut hx = (x.to_bits() >> 32) as i32;
    let lx = x.to_bits() as u32;
    let mut k = 0;
    if hx < 0x00100000 {
        if ((hx as u32 & 0x7fffffff) | lx) == 0 {
            return f64::NEG_INFINITY;
        }
        if hx < 0 {
            return f64::NAN;
        }
        k -= 54;
        x *= 1.80143985094819840000e16;
        hx = (x.to_bits() >> 32) as i32;
    }
    if hx >= 0x7ff00000 {
        return x + x;
    }
    k += (hx >> 20) - 1023;
    hx &= 0x000fffff;
    let i = (hx + 0x95f64) & 0x100000;
    x = f64::from_bits(((hx | (i ^ 0x3ff00000)) as u64) << 32 | (x.to_bits() & 0xffffffff));
    k += i >> 20;
    let f = x - 1.0;
    if (0x000fffff & (2 + hx)) < 3 {
        if f == 0.0 {
            if k == 0 {
                return 0.0;
            }
            let dk = f64::from(k);
            return dk * HI + dk * LO;
        }
        let r = f * f * (0.5 - 0.33333333333333333 * f);
        if k == 0 {
            return f - r;
        }
        let dk = f64::from(k);
        return dk * HI - ((r - dk * LO) - f);
    }
    let s = f / (2.0 + f);
    let dk = f64::from(k);
    let z = s * s;
    let i = hx - 0x6147a;
    let w = z * z;
    let j = 0x6b851 - hx;
    let t1 = w
        * (3.999999999940941908e-01
            + w * (2.222219843214978396e-01 + w * 1.531383769920937332e-01));
    let t2 = z
        * (6.666666666666735130e-01
            + w * (2.857142874366239149e-01
                + w * (1.818357216161805012e-01 + w * 1.479819860511658591e-01)));
    let r = t2 + t1;
    if (i | j) > 0 {
        let hfsq = 0.5 * f * f;
        if k == 0 {
            f - (hfsq - s * (hfsq + r))
        } else {
            dk * HI - ((hfsq - (s * (hfsq + r) + dk * LO)) - f)
        }
    } else if k == 0 {
        f - s * (f - r)
    } else {
        dk * HI - ((s * (f - r) - dk * LO) - f)
    }
}
