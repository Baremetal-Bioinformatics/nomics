//! OpenJDK 17-compatible formatting of binary64 values.
//!
//! The digit generator is a Rust port of the binary-to-ASCII half of
//! `jdk.internal.math.FloatingDecimal`. `format_fixed` and `format_sci` then
//! apply `FormattedFloatingDecimal`'s decimal-digit rounding and the padding
//! rules used by `java.util.Formatter`.

use std::cell::RefCell;
use std::cmp::Ordering;
use std::sync::OnceLock;

const EXP_SHIFT: i32 = 52;
const EXP_BIAS: i32 = 1023;
const FRACT_HOB: u64 = 1_u64 << EXP_SHIFT;
const SIGNIF_MASK: u64 = FRACT_HOB - 1;
const EXP_MASK: u64 = 0x7ff0_0000_0000_0000;
const SIGN_MASK: u64 = 1_u64 << 63;

#[derive(Clone)]
struct Decimal {
    negative: bool,
    dec_exp: i32,
    digits: Vec<u8>,
}

#[derive(Clone, Default, Eq, PartialEq)]
struct Big {
    // Little-endian base-2^32 limbs. An empty vector is zero.
    words: Vec<u32>,
}

impl Big {
    fn from_u64(value: u64) -> Self {
        if value == 0 {
            Self::default()
        } else {
            let mut words = vec![value as u32];
            let hi = (value >> 32) as u32;
            if hi != 0 {
                words.push(hi);
            }
            Self { words }
        }
    }

    fn normalize(&mut self) {
        while self.words.last() == Some(&0) {
            self.words.pop();
        }
    }

    fn mul_small(&mut self, factor: u32) {
        if factor == 0 || self.words.is_empty() {
            self.words.clear();
            return;
        }
        let mut carry = 0_u64;
        for word in &mut self.words {
            let product = u64::from(*word) * u64::from(factor) + carry;
            *word = product as u32;
            carry = product >> 32;
        }
        if carry != 0 {
            self.words.push(carry as u32);
        }
    }

    fn shl(&mut self, bits: usize) {
        if bits == 0 || self.words.is_empty() {
            return;
        }
        let word_shift = bits / 32;
        let bit_shift = bits % 32;
        let old = std::mem::take(&mut self.words);
        self.words = Vec::with_capacity(old.len() + word_shift + usize::from(bit_shift != 0));
        self.words.resize(word_shift, 0);
        if bit_shift == 0 {
            self.words.extend(old);
            return;
        }
        let mut carry = 0_u64;
        for word in old {
            let shifted = (u64::from(word) << bit_shift) | carry;
            self.words.push(shifted as u32);
            carry = shifted >> 32;
        }
        if carry != 0 {
            self.words.push(carry as u32);
        }
    }

    fn add_assign(&mut self, other: &Self) {
        let len = self.words.len().max(other.words.len());
        self.words.resize(len, 0);
        let mut carry = 0_u64;
        for i in 0..len {
            let sum = u64::from(self.words[i])
                + u64::from(other.words.get(i).copied().unwrap_or(0))
                + carry;
            self.words[i] = sum as u32;
            carry = sum >> 32;
        }
        if carry != 0 {
            self.words.push(carry as u32);
        }
    }

    fn sub_assign(&mut self, other: &Self) {
        debug_assert!((*self).cmp(other) != Ordering::Less);
        let mut borrow = 0_i64;
        for i in 0..self.words.len() {
            let difference = i64::from(self.words[i])
                - i64::from(other.words.get(i).copied().unwrap_or(0))
                - borrow;
            self.words[i] = difference as u32;
            borrow = i64::from(difference < 0);
        }
        debug_assert_eq!(borrow, 0);
        self.normalize();
    }

    fn cmp_sum(&self, a: &Self, b: &Self) -> Ordering {
        let mut sum = a.clone();
        sum.add_assign(b);
        self.cmp(&sum)
    }

    /// Return the next decimal digit and replace `self` with ten times the
    /// remainder. This is `FDBigInteger.quoRemIteration` for the range used by
    /// FloatingDecimal, where the quotient is a single digit.
    fn quo_rem_iteration(&mut self, divisor: &Self) -> u8 {
        let mut quotient = 0_u8;
        while (*self).cmp(divisor) != Ordering::Less {
            self.sub_assign(divisor);
            quotient += 1;
            debug_assert!(quotient < 10);
        }
        self.mul_small(10);
        quotient
    }
}

impl Ord for Big {
    fn cmp(&self, other: &Self) -> Ordering {
        match self.words.len().cmp(&other.words.len()) {
            Ordering::Equal => self.words.iter().rev().cmp(other.words.iter().rev()),
            order => order,
        }
    }
}

impl PartialOrd for Big {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

static POW5: OnceLock<Vec<Big>> = OnceLock::new();

fn pow5(power: usize) -> Big {
    let powers = POW5.get_or_init(|| {
        let mut powers = Vec::with_capacity(326);
        powers.push(Big::from_u64(1));
        for i in 1..326 {
            let mut next = powers[i - 1].clone();
            next.mul_small(5);
            powers.push(next);
        }
        powers
    });
    powers[power].clone()
}

fn pow52(power5: i32, power2: i32) -> Big {
    debug_assert!(power5 >= 0 && power2 >= 0);
    let mut value = pow5(power5 as usize);
    value.shl(power2 as usize);
    value
}

fn mul_pow52(value: u64, power5: i32, power2: i32) -> Big {
    // `value` is at most 53 bits. Shift-and-add keeps this helper independent
    // of machine-width multiplication.
    let mut result = Big::default();
    let mut addend = pow5(power5 as usize);
    let mut bits = value;
    while bits != 0 {
        if bits & 1 != 0 {
            result.add_assign(&addend);
        }
        bits >>= 1;
        if bits != 0 {
            addend.shl(1);
        }
    }
    result.shl(power2 as usize);
    result
}

#[allow(clippy::approx_constant)] // The deliberately truncated JDK 17 literal is the specification.
fn estimate_dec_exp(fract_bits: u64, bin_exp: i32) -> i32 {
    let d2 = f64::from_bits((u64::from(EXP_BIAS as u32) << EXP_SHIFT) | (fract_bits & SIGNIF_MASK));
    // Keep the evaluation order used by the Java expression; an occasional
    // one-off estimate is intentionally corrected by dtoa below.
    let estimate = (d2 - 1.5_f64) * 0.289_529_654_f64;
    let estimate = estimate + 0.176_091_259_f64;
    let estimate = estimate + f64::from(bin_exp) * 0.301_029_995_663_981_f64;
    let bits = estimate.to_bits();
    let exponent = ((bits & EXP_MASK) >> EXP_SHIFT) as i32 - EXP_BIAS;
    let negative = bits & SIGN_MASK != 0;
    if (0..52).contains(&exponent) {
        let mask = SIGNIF_MASK >> exponent;
        let integer = (((bits & SIGNIF_MASK) | FRACT_HOB) >> (EXP_SHIFT - exponent)) as i32;
        if negative {
            if mask & bits == 0 {
                -integer
            } else {
                -integer - 1
            }
        } else {
            integer
        }
    } else if exponent < 0 {
        if bits & !SIGN_MASK == 0 {
            0
        } else if negative {
            -1
        } else {
            0
        }
    } else {
        estimate as i32
    }
}

const INSIGNIFICANT_DIGITS: [u8; 64] = [
    0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 8, 8, 8, 9, 9,
    9, 9, 10, 10, 10, 11, 11, 11, 12, 12, 12, 12, 13, 13, 13, 14, 14, 14, 15, 15, 15, 15, 16, 16,
    16, 17, 17, 17, 18, 18, 18, 19,
];

fn develop_long_digits(mut value: u64, mut dec_exp: i32, insignificant: usize) -> Decimal {
    if insignificant != 0 {
        let power10 = 10_u64.pow(insignificant as u32);
        let residue = value % power10;
        value /= power10;
        dec_exp += insignificant as i32;
        if residue >= power10 >> 1 {
            value += 1;
        }
    }

    let mut reversed = Vec::with_capacity(20);
    let mut digit = (value % 10) as u8;
    value /= 10;
    while digit == 0 {
        dec_exp += 1;
        digit = (value % 10) as u8;
        value /= 10;
    }
    while value != 0 {
        reversed.push(digit);
        dec_exp += 1;
        digit = (value % 10) as u8;
        value /= 10;
    }
    reversed.push(digit);
    reversed.reverse();
    Decimal {
        negative: false,
        dec_exp: dec_exp + 1,
        digits: reversed,
    }
}

fn dtoa(bin_exp: i32, mut fract_bits: u64, significant_bits: i32, compatible: bool) -> Decimal {
    let trailing_zeros = fract_bits.trailing_zeros() as i32;
    let fract_bit_count = EXP_SHIFT + 1 - trailing_zeros;
    let tiny_bits = 0.max(fract_bit_count - bin_exp - 1);

    if (-21..=62).contains(&bin_exp) && tiny_bits == 0 {
        let insignificant = if bin_exp > significant_bits {
            let p2 = bin_exp - significant_bits - 1;
            if p2 > 1 && (p2 as usize) < INSIGNIFICANT_DIGITS.len() {
                usize::from(INSIGNIFICANT_DIGITS[p2 as usize])
            } else {
                0
            }
        } else {
            0
        };
        if bin_exp >= EXP_SHIFT {
            fract_bits <<= bin_exp - EXP_SHIFT;
        } else {
            fract_bits >>= EXP_SHIFT - bin_exp;
        }
        return develop_long_digits(fract_bits, 0, insignificant);
    }

    let mut dec_exp = estimate_dec_exp(fract_bits, bin_exp);
    let b5 = 0.max(-dec_exp);
    let mut b2 = b5 + tiny_bits + bin_exp;
    let s5 = 0.max(dec_exp);
    let mut s2 = s5 + tiny_bits;
    let m5 = b5;
    let mut m2 = b2 - significant_bits;

    fract_bits >>= trailing_zeros;
    b2 -= fract_bit_count - 1;
    let common2 = b2.min(s2);
    b2 -= common2;
    s2 -= common2;
    m2 -= common2;
    if fract_bit_count == 1 {
        m2 -= 1;
    }
    if m2 < 0 {
        b2 -= m2;
        s2 -= m2;
        m2 = 0;
    }

    let s = pow52(s5, s2);
    let mut b = mul_pow52(fract_bits, b5, b2);
    // The first quoRemIteration has multiplied the remainder by ten when the
    // JDK compares it with M, hence the extra powers here.
    let mut m = pow52(m5 + 1, m2 + 1);
    let ten_s = {
        let mut value = s.clone();
        value.mul_small(10);
        value
    };

    let mut digits = Vec::with_capacity(20);
    let digit = b.quo_rem_iteration(&s);
    let mut low = b < m;
    let mut high = ten_s.cmp_sum(&b, &m) != Ordering::Greater;
    debug_assert!(digit < 10);
    if digit == 0 && !high {
        dec_exp -= 1;
    } else {
        digits.push(digit);
    }
    if !compatible || !(-3..8).contains(&dec_exp) {
        low = false;
        high = false;
    }
    while !low && !high {
        let digit = b.quo_rem_iteration(&s);
        debug_assert!(digit < 10);
        m.mul_small(10);
        low = b < m;
        high = ten_s.cmp_sum(&b, &m) != Ordering::Greater;
        digits.push(digit);
    }

    if high {
        let round_up = if low {
            let mut twice_b = b.clone();
            twice_b.shl(1);
            match twice_b.cmp(&ten_s) {
                Ordering::Equal => digits.last().is_some_and(|digit| digit & 1 != 0),
                Ordering::Greater => true,
                Ordering::Less => false,
            }
        } else {
            true
        };
        if round_up {
            let mut i = digits.len() - 1;
            while digits[i] == 9 && i > 0 {
                digits[i] = 0;
                i -= 1;
            }
            if digits[i] == 9 {
                digits[0] = 1;
                dec_exp += 1;
            } else {
                digits[i] += 1;
            }
        }
    }

    Decimal {
        negative: false,
        dec_exp: dec_exp + 1,
        digits,
    }
}

thread_local! {
    static CONVERSION_CACHE: RefCell<Option<(u64, bool, Decimal)>> = const { RefCell::new(None) };
}

fn finite_decimal(d: f64, compatible: bool) -> Decimal {
    let bits = d.to_bits();
    CONVERSION_CACHE.with(|slot| {
        if let Some((cached_bits, cached_mode, decimal)) = &*slot.borrow()
            && *cached_bits == bits
            && *cached_mode == compatible
        {
            return decimal.clone();
        }

        let negative = bits & SIGN_MASK != 0;
        let mut fract_bits = bits & SIGNIF_MASK;
        let mut bin_exp = ((bits & EXP_MASK) >> EXP_SHIFT) as i32;
        let result = if bin_exp == 0 && fract_bits == 0 {
            Decimal {
                negative,
                dec_exp: 0,
                digits: vec![0],
            }
        } else {
            let significant_bits;
            if bin_exp == 0 {
                let leading_zeros = fract_bits.leading_zeros() as i32;
                let shift = leading_zeros - (63 - EXP_SHIFT);
                fract_bits <<= shift;
                bin_exp = 1 - shift;
                significant_bits = 64 - leading_zeros;
            } else {
                fract_bits |= FRACT_HOB;
                significant_bits = EXP_SHIFT + 1;
            }
            bin_exp -= EXP_BIAS;
            let mut decimal = dtoa(bin_exp, fract_bits, significant_bits, compatible);
            decimal.negative = negative;
            decimal
        };
        *slot.borrow_mut() = Some((bits, compatible, result.clone()));
        result
    })
}

fn apply_precision(decimal: &mut Decimal, precision: i32) {
    let digit_count = decimal.digits.len() as i32;
    if precision >= digit_count || precision < 0 {
        return;
    }
    if precision == 0 {
        if decimal.digits[0] >= 5 {
            decimal.digits[0] = 1;
            decimal.digits[1..].fill(0);
            decimal.dec_exp += 1;
        } else {
            decimal.digits.fill(0);
        }
        return;
    }

    let precision = precision as usize;
    if decimal.digits[precision] >= 5 {
        let mut i = precision - 1;
        while decimal.digits[i] == 9 && i > 0 {
            i -= 1;
        }
        if decimal.digits[i] == 9 {
            decimal.digits[0] = 1;
            decimal.digits[1..].fill(0);
            decimal.dec_exp += 1;
        } else {
            decimal.digits[i] += 1;
            decimal.digits[i + 1..].fill(0);
        }
    } else {
        decimal.digits[precision..].fill(0);
    }
}

fn push_sign(output: &mut String, negative: bool) {
    if negative {
        output.push('-');
    }
}

fn push_digits(output: &mut String, digits: &[u8]) {
    output.extend(digits.iter().map(|digit| char::from(b'0' + digit)));
}

fn add_zeros(output: &mut String, precision: usize) {
    let dot = output.find('.');
    let current = dot.map_or(0, |at| output.len() - at - 1);
    if current == precision {
        return;
    }
    if dot.is_none() {
        output.push('.');
    }
    output.extend(std::iter::repeat_n('0', precision - current));
}

fn exceptional(d: f64) -> Option<String> {
    if d.is_nan() {
        Some("NaN".to_owned())
    } else if d == f64::INFINITY {
        Some("Infinity".to_owned())
    } else if d == f64::NEG_INFINITY {
        Some("-Infinity".to_owned())
    } else {
        None
    }
}

/// Equivalent to `String.format(Locale.US, "%.<prec>f", d)`.
pub fn format_fixed(d: f64, prec: usize) -> String {
    if let Some(value) = exceptional(d) {
        return value;
    }
    let precision = i32::try_from(prec).expect("precision exceeds Java Formatter's integer range");
    let mut decimal = finite_decimal(d, false);
    let significant_precision = decimal.dec_exp.saturating_add(precision);
    apply_precision(&mut decimal, significant_precision);

    let mut output = String::new();
    push_sign(&mut output, decimal.negative);
    if decimal.dec_exp > 0 {
        let exp = decimal.dec_exp as usize;
        if decimal.digits.len() < exp {
            push_digits(&mut output, &decimal.digits);
            output.extend(std::iter::repeat_n('0', exp - decimal.digits.len()));
        } else {
            push_digits(&mut output, &decimal.digits[..exp]);
            let count = (decimal.digits.len() - exp).min(prec);
            if count != 0 {
                output.push('.');
                push_digits(&mut output, &decimal.digits[exp..exp + count]);
            }
        }
    } else {
        let zeros = ((-decimal.dec_exp).max(0) as usize).min(prec);
        let available = precision.saturating_add(decimal.dec_exp).max(0) as usize;
        let count = decimal.digits.len().min(available);
        output.push('0');
        if zeros != 0 || count != 0 {
            output.push('.');
            output.extend(std::iter::repeat_n('0', zeros));
            push_digits(&mut output, &decimal.digits[..count]);
        }
    }
    add_zeros(&mut output, prec);
    output
}

/// Equivalent to `String.format(Locale.US, "%.<prec>e", d)`.
pub fn format_sci(d: f64, prec: usize) -> String {
    if let Some(value) = exceptional(d) {
        return value;
    }
    let precision = i32::try_from(prec).expect("precision exceeds Java Formatter's integer range");
    let mut decimal = finite_decimal(d, false);
    apply_precision(&mut decimal, precision.saturating_add(1));

    let mut mantissa = String::new();
    push_sign(&mut mantissa, decimal.negative);
    mantissa.push(char::from(b'0' + decimal.digits[0]));
    let count = decimal.digits.len().saturating_sub(1).min(prec);
    if count != 0 {
        mantissa.push('.');
        push_digits(&mut mantissa, &decimal.digits[1..=count]);
    }
    add_zeros(&mut mantissa, prec);

    let exponent = if d == 0.0 { 0 } else { decimal.dec_exp - 1 };
    let sign = if exponent < 0 { '-' } else { '+' };
    let magnitude = exponent.unsigned_abs();
    mantissa.push('e');
    mantissa.push(sign);
    if magnitude < 10 {
        mantissa.push('0');
    }
    mantissa.push_str(&magnitude.to_string());
    mantissa
}

/// Equivalent to OpenJDK 17's `Double.toString(d)`.
pub fn double_to_string(d: f64) -> String {
    if let Some(value) = exceptional(d) {
        return value;
    }
    let decimal = finite_decimal(d, true);
    let mut output = String::new();
    push_sign(&mut output, decimal.negative);
    if decimal.dec_exp > 0 && decimal.dec_exp < 8 {
        let exp = decimal.dec_exp as usize;
        let count = decimal.digits.len().min(exp);
        push_digits(&mut output, &decimal.digits[..count]);
        if count < exp {
            output.extend(std::iter::repeat_n('0', exp - count));
            output.push_str(".0");
        } else {
            output.push('.');
            if count < decimal.digits.len() {
                push_digits(&mut output, &decimal.digits[count..]);
            } else {
                output.push('0');
            }
        }
    } else if decimal.dec_exp <= 0 && decimal.dec_exp > -3 {
        output.push_str("0.");
        output.extend(std::iter::repeat_n('0', (-decimal.dec_exp) as usize));
        push_digits(&mut output, &decimal.digits);
    } else {
        output.push(char::from(b'0' + decimal.digits[0]));
        output.push('.');
        if decimal.digits.len() > 1 {
            push_digits(&mut output, &decimal.digits[1..]);
        } else {
            output.push('0');
        }
        output.push('E');
        let exponent = decimal.dec_exp - 1;
        if exponent < 0 {
            output.push('-');
        }
        output.push_str(&exponent.unsigned_abs().to_string());
    }
    output
}
