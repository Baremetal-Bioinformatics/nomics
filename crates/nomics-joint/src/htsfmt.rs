//! The text path after GenotypeGVCFs: sarek pipes each interval's VCF
//! through `bcftools sort` (htslib 1.23.1), which parses every Float field
//! into a 32-bit float and prints it back with `kputd`, and then through
//! GATK `MergeVcfs`, which keeps INFO strings verbatim but re-prints QUAL
//! with htsjdk's `formatQualValue`.

use nomics_hc::vcf::format_qual;

/// htslib `hts_str2dbl` (textutils_internal.h): a quick decimal parser that
/// falls back to `strtod` for exponents, long inputs and non-numbers.
pub fn hts_str2dbl(s: &str) -> f64 {
    let b = s.trim_start().as_bytes();
    let (neg, mut v) = match b.first() {
        Some(b'-') => (true, &b[1..]),
        Some(b'+') => (false, &b[1..]),
        _ => (false, b),
    };
    match v.first() {
        Some(b'1'..=b'9') => {}
        Some(b'0') if !matches!(v.get(1), Some(b'x' | b'X')) => {}
        _ => return strtod(s),
    }
    while v.first() == Some(&b'0') {
        v = &v[1..];
    }
    let mut n: u64 = 0;
    let mut max_len = 15;
    let mut i = 0;
    let mut point: Option<usize> = None;
    loop {
        max_len -= 1;
        if max_len == 0 || !matches!(v.get(i), Some(b'0'..=b'9')) {
            break;
        }
        n = n * 10 + u64::from(v[i] - b'0');
        i += 1;
    }
    if max_len != 0 && v.get(i) == Some(&b'.') {
        point = Some(i);
        i += 1;
        loop {
            max_len -= 1;
            if max_len == 0 || !matches!(v.get(i), Some(b'0'..=b'9')) {
                break;
            }
            n = n * 10 + u64::from(v[i] - b'0');
            i += 1;
        }
    }
    let point = point.unwrap_or(i);
    if max_len == 0 || matches!(v.get(i), Some(b'e' | b'E')) {
        return strtod(s);
    }
    const D: [f64; 22] = [
        1.0, 1.0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12, 1e13, 1e14,
        1e15, 1e16, 1e17, 1e18, 1e19, 1e20,
    ];
    let d = n as f64 / D[i - point];
    if neg { -d } else { d }
}

/// `strtod` on the longest numeric prefix (correctly rounded, as glibc is).
fn strtod(s: &str) -> f64 {
    let t = s.trim_start();
    let end = t
        .char_indices()
        .take_while(|&(_, c)| c.is_ascii_digit() || matches!(c, '+' | '-' | '.' | 'e' | 'E'))
        .map(|(i, c)| i + c.len_utf8())
        .last()
        .unwrap_or(0);
    (0..=end)
        .rev()
        .find_map(|e| t[..e].parse::<f64>().ok())
        .unwrap_or(0.0)
}

/// htslib `kputd`: at most six significant digits, trailing zeros culled;
/// stdio `%g` outside [0.0001, 999999].
pub fn kputd(mut d: f64) -> String {
    if d == 0.0 {
        return if d.is_sign_negative() { "-0" } else { "0" }.to_string();
    }
    let mut out = String::new();
    if d < 0.0 {
        out.push('-');
        d = -d;
    }
    if !(0.0001..=999999.0).contains(&d) {
        out.push_str(&c_format_g(d));
        return out;
    }
    let (mut i, k): (u32, usize) = if d < 0.001 {
        ((d * 1e9).round_ties_even() as u32, 1)
    } else if d < 0.01 {
        ((d * 1e8).round_ties_even() as u32, 2)
    } else if d < 0.1 {
        ((d * 1e7).round_ties_even() as u32, 3)
    } else if d < 1.0 {
        ((d * 1e6).round_ties_even() as u32, 4)
    } else if d < 10.0 {
        ((d * 1e5).round_ties_even() as u32, 5)
    } else if d < 100.0 {
        ((d * 1e4).round_ties_even() as u32, 6)
    } else if d < 1000.0 {
        ((d * 1e3).round_ties_even() as u32, 7)
    } else if d < 10000.0 {
        ((d * 1e2).round_ties_even() as u32, 8)
    } else if d < 100000.0 {
        ((d * 1e1).round_ties_even() as u32, 9)
    } else {
        (d.round_ties_even() as u32, 10)
    };
    let mut buf = [0u8; 21];
    let mut cp = 20 - k;
    for last in [false, false, true] {
        cp -= 2;
        let pair = i % 100;
        buf[cp] = b'0' + (pair / 10) as u8;
        buf[cp + 1] = b'0' + (pair % 10) as u8;
        if !last {
            i /= 100;
        }
    }
    if i >= 100 {
        cp -= 1;
        buf[cp] = b'0' + (i / 100) as u8;
    }
    let mut p = 20 - cp;
    let mut ep;
    if p <= 10 {
        ep = cp + 5;
        while p < 10 {
            cp -= 1;
            buf[cp] = b'0';
            p += 1;
        }
        cp -= 1;
        buf[cp] = b'.';
        cp -= 1;
        buf[cp] = b'0';
    } else {
        cp -= 1;
        let mut xp = cp;
        ep = cp + 6;
        while p > 10 {
            buf[xp] = buf[xp + 1];
            xp += 1;
            p -= 1;
        }
        buf[xp] = b'.';
    }
    while buf[ep] == b'0' && ep > cp {
        ep -= 1;
    }
    if buf[ep] != 0 && buf[ep] != b'.' {
        ep += 1;
    }
    out.push_str(std::str::from_utf8(&buf[cp..ep]).expect("ASCII digits"));
    out
}

/// C `printf("%g", d)` for finite positive `d`.
fn c_format_g(d: f64) -> String {
    if d.is_infinite() {
        return "inf".into();
    }
    if d.is_nan() {
        return "nan".into();
    }
    let sci = format!("{d:.5e}");
    let (mantissa, exp) = sci.split_once('e').expect("exponent");
    let exp: i32 = exp.parse().expect("exponent digits");
    if exp < -4 || exp >= 6 {
        let m = trim_fraction(mantissa);
        let sign = if exp < 0 { '-' } else { '+' };
        format!("{m}e{sign}{:02}", exp.abs())
    } else {
        trim_fraction(&format!("{d:.*}", (5 - exp) as usize)).to_string()
    }
}

fn trim_fraction(s: &str) -> &str {
    if s.contains('.') {
        s.trim_end_matches('0').trim_end_matches('.')
    } else {
        s
    }
}

/// A VCF Float field after a `bcftools` round trip.
pub fn bcf_float(s: &str) -> String {
    if s == "." {
        return s.to_string();
    }
    kputd(f64::from(hts_str2dbl(s) as f32))
}

/// QUAL after `bcftools` (`atof`, stored as float, `kputd`) and then
/// htsjdk (`Double.valueOf`, `log10PError = qual / -10`, `formatQualValue`).
pub fn merged_qual(gatk: &str) -> String {
    let bcf = kputd(f64::from(strtod(gatk) as f32));
    let qual: f64 = bcf.parse().expect("kputd output parses");
    let log10_p_error = qual / -10.0;
    format_qual(-10.0 * log10_p_error)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn kputd_matches_htslib() {
        for (d, s) in [
            (3.2_f32, "3.2"),
            (0.5, "0.5"),
            (22.0, "22"),
            (25.36, "25.36"),
            (-0.187, "-0.187"),
            (121.84, "121.84"),
            (1492.02, "1492.02"),
            (123456.7, "123457"),
            (0.00012, "0.00012"),
            (0.00001, "1e-05"),
            (2.0e6, "2e+06"),
            (0.9999999, "1"),
        ] {
            assert_eq!(kputd(f64::from(d)), s, "{d}");
        }
    }

    #[test]
    fn bcf_float_round_trips_gatk_strings() {
        assert_eq!(bcf_float("3.20"), "3.2");
        assert_eq!(bcf_float("0.000"), "0");
        assert_eq!(bcf_float("-3.816e+00"), "-3.816");
        assert_eq!(bcf_float("0.500"), "0.5");
        assert_eq!(bcf_float("2.833"), "2.833");
    }

    #[test]
    fn merged_qual_formats() {
        assert_eq!(merged_qual("121.84"), "121.84");
        assert_eq!(merged_qual("30"), "30");
        assert_eq!(merged_qual("12345.67"), "12345.70");
    }
}
