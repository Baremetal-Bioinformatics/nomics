//! ITF-8 and LTF-8 variable-length integers (CRAMv3 §2.3).

use std::io;

/// Append `v` as ITF-8 (1-5 bytes; negative values take 5).
#[inline]
pub fn put_itf8(out: &mut Vec<u8>, v: i32) {
    let u = v as u32;
    if u < 0x80 {
        out.push(u as u8);
    } else if u < 0x4000 {
        out.extend_from_slice(&[0x80 | (u >> 8) as u8, u as u8]);
    } else if u < 0x20_0000 {
        out.extend_from_slice(&[0xC0 | (u >> 16) as u8, (u >> 8) as u8, u as u8]);
    } else if u < 0x1000_0000 {
        out.extend_from_slice(&[
            0xE0 | (u >> 24) as u8,
            (u >> 16) as u8,
            (u >> 8) as u8,
            u as u8,
        ]);
    } else {
        out.extend_from_slice(&[
            0xF0 | ((u >> 28) & 0x0f) as u8,
            (u >> 20) as u8,
            (u >> 12) as u8,
            (u >> 4) as u8,
            (u & 0x0f) as u8,
        ]);
    }
}

/// Append `v` as LTF-8 (1-9 bytes).
pub fn put_ltf8(out: &mut Vec<u8>, v: i64) {
    let u = v as u64;
    let n = if u < 1 << 7 {
        1
    } else if u < 1 << 14 {
        2
    } else if u < 1 << 21 {
        3
    } else if u < 1 << 28 {
        4
    } else if u < 1 << 35 {
        5
    } else if u < 1 << 42 {
        6
    } else if u < 1 << 49 {
        7
    } else if u < 1 << 56 {
        8
    } else {
        9
    };
    if n == 9 {
        out.push(0xFF);
        out.extend_from_slice(&u.to_be_bytes());
        return;
    }
    if n == 8 {
        out.push(0xFE);
        out.extend_from_slice(&u.to_be_bytes()[1..]);
        return;
    }
    // First byte: (n-1) leading ones, then the top bits of the value.
    let prefix: u8 = !(0xFFu8 >> (n - 1));
    let bytes = u.to_be_bytes();
    let tail = &bytes[8 - (n - 1)..];
    let top = (u >> (8 * (n - 1))) as u8;
    out.push(prefix | top);
    out.extend_from_slice(tail);
}

fn eof() -> io::Error {
    io::Error::new(io::ErrorKind::UnexpectedEof, "truncated ITF-8/LTF-8")
}

/// Read an ITF-8 integer from `b[*at..]`, advancing `at`.
pub fn get_itf8(b: &[u8], at: &mut usize) -> io::Result<i32> {
    let b0 = *b.get(*at).ok_or_else(eof)? as u32;
    let n = if b0 < 0x80 {
        0
    } else if b0 < 0xC0 {
        1
    } else if b0 < 0xE0 {
        2
    } else if b0 < 0xF0 {
        3
    } else {
        4
    };
    let rest = b.get(*at + 1..*at + 1 + n).ok_or_else(eof)?;
    *at += 1 + n;
    let r = |i: usize| rest[i] as u32;
    Ok(match n {
        0 => b0,
        1 => ((b0 & 0x3f) << 8) | r(0),
        2 => ((b0 & 0x1f) << 16) | (r(0) << 8) | r(1),
        3 => ((b0 & 0x0f) << 24) | (r(0) << 16) | (r(1) << 8) | r(2),
        _ => ((b0 & 0x0f) << 28) | (r(0) << 20) | (r(1) << 12) | (r(2) << 4) | (r(3) & 0x0f),
    } as i32)
}

/// Read an LTF-8 integer from `b[*at..]`, advancing `at`.
pub fn get_ltf8(b: &[u8], at: &mut usize) -> io::Result<i64> {
    let b0 = *b.get(*at).ok_or_else(eof)?;
    let n = b0.leading_ones() as usize; // bytes that follow
    let rest = b.get(*at + 1..*at + 1 + n).ok_or_else(eof)?;
    *at += 1 + n;
    let mut v: u64 = if n >= 7 {
        0
    } else {
        (b0 & (0xFFu8 >> (n + 1))) as u64
    };
    for &x in rest {
        v = (v << 8) | x as u64;
    }
    Ok(v as i64)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn itf8_known_bytes() {
        let cases: &[(i32, &[u8])] = &[
            (0, &[0x00]),
            (0x7f, &[0x7f]),
            (0x80, &[0x80, 0x80]),
            (0x3fff, &[0xbf, 0xff]),
            (0x4000, &[0xc0, 0x40, 0x00]),
            (4542278, &[0xe0, 0x45, 0x4f, 0x46]), // EOF container alignment start
            (-1, &[0xff, 0xff, 0xff, 0xff, 0x0f]),
            (0x4F515A, &[0xe0, 0x4f, 0x51, 0x5a]), // spec: OQ:Z tag key
        ];
        for &(v, bytes) in cases {
            let mut out = Vec::new();
            put_itf8(&mut out, v);
            assert_eq!(out, bytes, "{v}");
            let mut at = 0;
            assert_eq!(get_itf8(&out, &mut at).unwrap(), v);
            assert_eq!(at, bytes.len());
        }
    }

    #[test]
    fn itf8_round_trips() {
        let mut v: Vec<i32> = vec![i32::MIN, i32::MAX, -2, 1 << 28, (1 << 28) - 1];
        for s in 0..32 {
            v.push(1 << s);
            v.push(((1i64 << s) - 1) as i32);
            v.push((-(1i64 << s)) as i32);
        }
        for x in v {
            let mut out = Vec::new();
            put_itf8(&mut out, x);
            let mut at = 0;
            assert_eq!(get_itf8(&out, &mut at).unwrap(), x);
            assert_eq!(at, out.len());
        }
    }

    #[test]
    fn ltf8_round_trips() {
        let mut v: Vec<i64> = vec![0, -1, i64::MIN, i64::MAX];
        for s in 0..64 {
            v.push(1i64 << s);
            v.push(((1u64 << s) - 1) as i64);
        }
        for x in v {
            let mut out = Vec::new();
            put_ltf8(&mut out, x);
            let mut at = 0;
            assert_eq!(get_ltf8(&out, &mut at).unwrap(), x, "{x:#x}");
            assert_eq!(at, out.len());
        }
        let mut out = Vec::new();
        put_ltf8(&mut out, 0x80);
        assert_eq!(out, [0x80, 0x80]);
        out.clear();
        put_ltf8(&mut out, 1 << 35);
        assert_eq!(out, [0xF8, 0x08, 0, 0, 0, 0]);
    }
}
