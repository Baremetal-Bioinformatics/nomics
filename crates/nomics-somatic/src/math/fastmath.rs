//! Operation-ordered real arithmetic from Commons Math 3.5 util/FastMath.java.
//! Only the overloads used by the contamination chain are implemented.
#![allow(clippy::excessive_precision)]
use super::fastmath_tables::*;

pub fn exp(x: f64) -> f64 {
    if x.is_nan() {
        return x;
    }
    let mut int_val = x as i32;
    if x < 0.0 {
        if x < -746.0 {
            return 0.0;
        }
        if int_val < -709 {
            return exp(x + 40.19140625) / 285040095144011776.0;
        }
        if int_val == -709 {
            return exp(x + 1.494140625) / 4.455505956692756620;
        }
        int_val -= 1;
    } else if int_val > 709 {
        return f64::INFINITY;
    }
    let a = EXP_INT_A[(750 + int_val) as usize];
    let b = EXP_INT_B[(750 + int_val) as usize];
    let frac = ((x - f64::from(int_val)) * 1024.0) as usize;
    let c = EXP_FRAC_A[frac];
    let d = EXP_FRAC_B[frac];
    let epsilon = x - (f64::from(int_val) + frac as f64 / 1024.0);
    let mut z = 0.04168701738764507;
    z = z * epsilon + 0.1666666505023083;
    z = z * epsilon + 0.5000000000042687;
    z = z * epsilon + 1.0;
    z = z * epsilon + -3.940510424527919E-20;
    let ta = a * c;
    let tb = a * d + b * c + b * d;
    let tc = tb + ta;
    tc * z + tb + ta
}

pub fn log(x: f64) -> f64 {
    if x == 0.0 {
        return f64::NEG_INFINITY;
    }
    if x < 0.0 || x.is_nan() {
        return f64::NAN;
    }
    if x == f64::INFINITY {
        return x;
    }
    let mut bits = x.to_bits();
    let mut exponent = (bits >> 52) as i32 - 1023;
    if bits & 0x7ff0000000000000 == 0 {
        bits <<= 1;
        while bits & 0x0010000000000000 == 0 {
            exponent -= 1;
            bits <<= 1;
        }
    }
    if (exponent == -1 || exponent == 0) && x < 1.01 && x > 0.99 {
        const COEF: [[f64; 2]; 9] = [
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
        let xa = x - 1.0;
        let tmp = xa * 1073741824.0;
        let aa = xa + tmp - tmp;
        let xb = xa - aa;
        let xa = aa;
        let mut ya = COEF[8][0];
        let mut yb = COEF[8][1];
        for coef in COEF[..8].iter().rev() {
            let aa = ya * xa;
            let ab = ya * xb + yb * xa + yb * xb;
            let tmp = aa * 1073741824.0;
            ya = aa + tmp - tmp;
            yb = aa - ya + ab;
            let aa = ya + coef[0];
            let ab = yb + coef[1];
            let tmp = aa * 1073741824.0;
            ya = aa + tmp - tmp;
            yb = aa - ya + ab;
        }
        let aa = ya * xa;
        let ab = ya * xb + yb * xa + yb * xb;
        let tmp = aa * 1073741824.0;
        ya = aa + tmp - tmp;
        yb = aa - ya + ab;
        return ya + yb;
    }
    let lnm = LN_MANT[((bits & 0x000ffc0000000000) >> 42) as usize];
    let epsilon =
        (bits & 0x3ffffffffff) as f64 / (4503599627370496.0 + (bits & 0x000ffc0000000000) as f64);
    let mut lnz = -0.16624882440418567;
    lnz = lnz * epsilon + 0.19999954120254515;
    lnz = lnz * epsilon + -0.2499999997677497;
    lnz = lnz * epsilon + 0.3333333333332802;
    lnz = lnz * epsilon + -0.5;
    lnz = lnz * epsilon + 1.0;
    lnz *= epsilon;
    let exponent = f64::from(exponent);
    let mut a = 0.693147063255310059 * exponent;
    let mut b = 0.0;
    for term in [lnm[0], lnz, 1.17304635250823482e-7 * exponent, lnm[1], 0.0] {
        let c = a + term;
        let d = -(c - a - term);
        a = c;
        b += d;
    }
    a + b
}

// FastMath.hypot uses exponent scaling and sqrt, *not* the platform hypot.
pub fn hypot(x: f64, y: f64) -> f64 {
    if x.is_infinite() || y.is_infinite() {
        return f64::INFINITY;
    }
    if x.is_nan() || y.is_nan() {
        return f64::NAN;
    }
    let ex = ((x.to_bits() >> 52) & 2047) as i32 - 1023;
    let ey = ((y.to_bits() >> 52) & 2047) as i32 - 1023;
    if ex > ey + 27 {
        return x.abs();
    }
    if ey > ex + 27 {
        return y.abs();
    }
    let middle = (ex + ey) / 2;
    let sx = scalb(x, -middle);
    let sy = scalb(y, -middle);
    scalb((sx * sx + sy * sy).sqrt(), middle)
}
// Exact power-of-two scaling, including subnormal rounding; FastMath.scalb.
fn scalb(d: f64, n: i32) -> f64 {
    if n > -1023 && n < 1024 {
        return d * f64::from_bits(((n + 1023) as u64) << 52);
    }
    if d.is_nan() || d.is_infinite() || d == 0.0 {
        return d;
    }
    if n < -2098 {
        return 0.0_f64.copysign(d);
    }
    if n > 2097 {
        return f64::INFINITY.copysign(d);
    }
    let bits = d.to_bits();
    let sign = bits & 0x8000000000000000;
    let exponent = ((bits >> 52) & 0x7ff) as i32;
    let mut mantissa = bits & 0x000fffffffffffff;
    let mut scaled = exponent + n;
    if n < 0 {
        if scaled > 0 {
            return f64::from_bits(sign | ((scaled as u64) << 52) | mantissa);
        }
        if scaled > -53 {
            mantissa |= 1u64 << 52;
            let lost = mantissa & (1u64 << (-scaled));
            mantissa >>= 1 - scaled;
            if lost != 0 {
                mantissa += 1;
            }
            return f64::from_bits(sign | mantissa);
        }
        0.0_f64.copysign(d)
    } else {
        if exponent == 0 {
            while mantissa >> 52 != 1 {
                mantissa <<= 1;
                scaled -= 1;
            }
            scaled += 1;
            mantissa &= 0x000fffffffffffff;
        }
        if scaled < 2047 {
            f64::from_bits(sign | ((scaled as u64) << 52) | mantissa)
        } else {
            f64::INFINITY.copysign(d)
        }
    }
}
