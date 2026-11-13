//! Fresh implementation of OpenJDK 17 FloatingDecimal's binary-to-decimal
//! decisions, confirmed against the reference text. This is intentionally not Rust's
//! Display: Java's symmetric power-of-two bounds and small-integer paths differ.
use std::cmp::Ordering;

#[allow(clippy::approx_constant)] // The truncated upstream coefficient is observable.
pub(super) fn java(v: f64) -> String {
    if v.is_nan() {
        return "NaN".into();
    }
    if v.is_infinite() {
        return if v.is_sign_negative() {
            "-Infinity"
        } else {
            "Infinity"
        }
        .into();
    }
    if v == 0.0 {
        return if v.is_sign_negative() { "-0.0" } else { "0.0" }.into();
    }
    let bits = v.abs().to_bits();
    let mut fraction = bits & ((1 << 52) - 1);
    let mut exponent = ((bits >> 52) & 2047) as i32;
    let significant;
    if exponent == 0 {
        let shift = fraction.leading_zeros() as i32 - 11;
        significant = 64 - fraction.leading_zeros() as i32;
        fraction <<= shift;
        exponent = 1 - shift;
    } else {
        fraction |= 1 << 52;
        significant = 53;
    }
    exponent -= 1023;
    let zeros = fraction.trailing_zeros() as i32;
    let nbits = 53 - zeros;
    let tiny = (nbits - exponent - 1).max(0);
    let (digits, exp) = if (-21..=62).contains(&exponent) && tiny == 0 {
        let mut n = if exponent >= 52 {
            fraction << (exponent - 52)
        } else {
            fraction >> (52 - exponent)
        };
        let power = exponent - significant - 1;
        let insignificant = if exponent > significant && power > 1 {
            (1u64 << power).ilog10()
        } else {
            0
        };
        let mut exp = 0;
        if insignificant > 0 {
            let scale = 10u64.pow(insignificant);
            let r = n % scale;
            n /= scale;
            if r >= scale / 2 {
                n += 1;
            }
            exp += insignificant as i32;
        }
        let s = n.to_string();
        exp += s.len() as i32;
        (s.trim_end_matches('0').as_bytes().to_vec(), exp)
    } else {
        let d2 = f64::from_bits((1023u64 << 52) | (fraction & ((1 << 52) - 1)));
        let mut decimal =
            ((d2 - 1.5) * 0.289529654 + 0.176091259 + exponent as f64 * 0.301029995663981).floor()
                as i32;
        let b5 = (-decimal).max(0);
        let s5 = decimal.max(0);
        let mut b2 = b5 + tiny + exponent - (nbits - 1);
        let mut s2 = s5 + tiny;
        let mut m2 = b5 + tiny + exponent - significant;
        let common = b2.min(s2);
        b2 -= common;
        s2 -= common;
        m2 -= common;
        if nbits == 1 {
            m2 -= 1;
        }
        if m2 < 0 {
            b2 -= m2;
            s2 -= m2;
            m2 = 0;
        }
        let pow_bits = |p: i32| {
            if p == 0 {
                0
            } else if p < 27 {
                64 - 5u64.pow(p as u32).leading_zeros() as i32
            } else {
                p * 3
            }
        };
        let bbits = nbits + b2 + pow_bits(b5);
        let sbits = s2 + 1 + pow_bits(s5 + 1);
        let width = if bbits < 32 && sbits < 32 {
            32
        } else if bbits < 64 && sbits < 64 {
            64
        } else {
            0
        };
        let mut digits = Vec::new();
        let (low, high, diff) = if width > 0 {
            let wrap = |x: i64| if width == 32 { x as i32 as i64 } else { x };
            let mut b = wrap(
                (fraction >> zeros)
                    .wrapping_mul(5u64.pow(b5 as u32))
                    .wrapping_shl(b2 as u32) as i64,
            );
            let s = wrap(5u64.pow(s5 as u32).wrapping_shl(s2 as u32) as i64);
            let mut m = wrap(5u64.pow(b5 as u32).wrapping_shl(m2 as u32) as i64);
            let ten = wrap(s.wrapping_mul(10));
            let mut first = true;
            loop {
                let q = b / s;
                b = wrap((b % s).wrapping_mul(10));
                m = wrap(m.wrapping_mul(10));
                let mut low = b < m;
                let mut high = wrap(b.wrapping_add(m)) > ten;
                if !first && m <= 0 {
                    low = true;
                    high = true;
                }
                if first && q == 0 && !high {
                    decimal -= 1;
                } else {
                    digits.push(b'0' + q as u8);
                }
                if first && !(-3..8).contains(&decimal) {
                    low = false;
                    high = false;
                }
                first = false;
                if low || high {
                    break (
                        low,
                        high,
                        wrap(wrap(b.wrapping_shl(1)).wrapping_sub(ten)).cmp(&0),
                    );
                }
            }
        } else {
            let mut b = Big::scaled(fraction >> zeros, b5, b2);
            let s = Big::scaled(1, s5, s2);
            let mut m = Big::scaled(1, b5 + 1, m2 + 1);
            let ten = Big::scaled(1, s5 + 1, s2 + 1);
            let mut first = true;
            loop {
                let mut q = 0;
                while b.cmp(&s) != Ordering::Less {
                    b.sub(&s);
                    q += 1;
                }
                b.mul(10);
                if !first {
                    m.mul(10);
                }
                let mut low = b.cmp(&m) == Ordering::Less;
                let mut sum = b.clone();
                sum.add(&m);
                let mut high = ten.cmp(&sum) != Ordering::Greater;
                if first && q == 0 && !high {
                    decimal -= 1;
                } else {
                    digits.push(b'0' + q);
                }
                if first && !(-3..8).contains(&decimal) {
                    low = false;
                    high = false;
                }
                first = false;
                if low || high {
                    b.mul(2);
                    break (low, high, b.cmp(&ten));
                }
            }
        };
        if high
            && (!low
                || diff == Ordering::Greater
                || (diff == Ordering::Equal && digits.last().unwrap() % 2 != 0))
        {
            let mut i = digits.len() - 1;
            while i > 0 && digits[i] == b'9' {
                digits[i] = b'0';
                i -= 1;
            }
            if digits[i] == b'9' {
                digits[0] = b'1';
                decimal += 1;
            } else {
                digits[i] += 1;
            }
        }
        (digits, decimal + 1)
    };
    let mut out = if v.is_sign_negative() {
        String::from("-")
    } else {
        String::new()
    };
    let s = std::str::from_utf8(&digits).unwrap();
    if exp > 0 && exp < 8 {
        let n = digits.len().min(exp as usize);
        out.push_str(&s[..n]);
        if n < exp as usize {
            out.push_str(&"0".repeat(exp as usize - n));
            out.push_str(".0");
        } else {
            out.push('.');
            out.push_str(if n < digits.len() { &s[n..] } else { "0" });
        }
    } else if exp <= 0 && exp > -3 {
        out.push_str("0.");
        out.push_str(&"0".repeat(-exp as usize));
        out.push_str(s);
    } else {
        out.push(digits[0] as char);
        out.push('.');
        out.push_str(if s.len() > 1 { &s[1..] } else { "0" });
        out.push('E');
        out.push_str(&(exp - 1).to_string());
    }
    out
}
// Only used while rendering, never per read. Base 2^32, little endian.
#[derive(Clone)]
struct Big(Vec<u32>);
impl Big {
    fn scaled(n: u64, fives: i32, twos: i32) -> Self {
        let mut b = Self(vec![n as u32, (n >> 32) as u32]);
        b.trim();
        for _ in 0..fives {
            b.mul(5);
        }
        for _ in 0..twos {
            b.mul(2);
        }
        b
    }
    fn trim(&mut self) {
        while self.0.len() > 1 && self.0.last() == Some(&0) {
            self.0.pop();
        }
    }
    fn cmp(&self, other: &Self) -> Ordering {
        self.0
            .len()
            .cmp(&other.0.len())
            .then_with(|| self.0.iter().rev().cmp(other.0.iter().rev()))
    }
    fn mul(&mut self, n: u32) {
        let mut carry = 0u64;
        for x in &mut self.0 {
            let v = *x as u64 * n as u64 + carry;
            *x = v as u32;
            carry = v >> 32;
        }
        if carry > 0 {
            self.0.push(carry as u32);
        }
    }
    fn sub(&mut self, other: &Self) {
        let mut borrow = 0i64;
        for (i, x) in self.0.iter_mut().enumerate() {
            let v = *x as i64 - other.0.get(i).copied().unwrap_or(0) as i64 - borrow;
            *x = v as u32;
            borrow = i64::from(v < 0);
        }
        self.trim();
    }
    fn add(&mut self, other: &Self) {
        self.0.resize(self.0.len().max(other.0.len()), 0);
        let mut carry = 0u64;
        for (i, x) in self.0.iter_mut().enumerate() {
            let v = *x as u64 + other.0.get(i).copied().unwrap_or(0) as u64 + carry;
            *x = v as u32;
            carry = v >> 32;
        }
        if carry > 0 {
            self.0.push(carry as u32);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::java;
    #[test]
    fn legacy_jvm_rounding_and_layout() {
        // OpenJDK 17 FloatOracle; includes values where modern shortest output
        // differs, scientific cutoffs, signed zero and the smallest subnormal.
        for (v, text) in [
            (1e23, "9.999999999999999E22"),
            (1e18, "1.0E18"),
            (0.001, "0.001"),
            (0.0001, "1.0E-4"),
            (10_000_000.0, "1.0E7"),
            (f64::from_bits(1), "4.9E-324"),
            (36.407130129813, "36.407130129813"),
            (-0.0, "-0.0"),
            (f64::NAN, "NaN"),
            (f64::INFINITY, "Infinity"),
        ] {
            assert_eq!(java(v), text);
        }
    }
}
