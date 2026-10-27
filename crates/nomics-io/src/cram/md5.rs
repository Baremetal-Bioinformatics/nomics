//! MD5 (RFC 1321), for `@SQ M5:` and CRAM slice reference checksums.

const S: [u32; 64] = [
    7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17, 22, 5, 9, 14, 20, 5, 9, 14, 20, 5, 9,
    14, 20, 5, 9, 14, 20, 4, 11, 16, 23, 4, 11, 16, 23, 4, 11, 16, 23, 4, 11, 16, 23, 6, 10, 15,
    21, 6, 10, 15, 21, 6, 10, 15, 21, 6, 10, 15, 21,
];

const K: [u32; 64] = {
    // floor(abs(sin(i + 1)) * 2^32), precomputed.
    [
        0xd76aa478, 0xe8c7b756, 0x242070db, 0xc1bdceee, 0xf57c0faf, 0x4787c62a, 0xa8304613,
        0xfd469501, 0x698098d8, 0x8b44f7af, 0xffff5bb1, 0x895cd7be, 0x6b901122, 0xfd987193,
        0xa679438e, 0x49b40821, 0xf61e2562, 0xc040b340, 0x265e5a51, 0xe9b6c7aa, 0xd62f105d,
        0x02441453, 0xd8a1e681, 0xe7d3fbc8, 0x21e1cde6, 0xc33707d6, 0xf4d50d87, 0x455a14ed,
        0xa9e3e905, 0xfcefa3f8, 0x676f02d9, 0x8d2a4c8a, 0xfffa3942, 0x8771f681, 0x6d9d6122,
        0xfde5380c, 0xa4beea44, 0x4bdecfa9, 0xf6bb4b60, 0xbebfbc70, 0x289b7ec6, 0xeaa127fa,
        0xd4ef3085, 0x04881d05, 0xd9d4d039, 0xe6db99e5, 0x1fa27cf8, 0xc4ac5665, 0xf4292244,
        0x432aff97, 0xab9423a7, 0xfc93a039, 0x655b59c3, 0x8f0ccc92, 0xffeff47d, 0x85845dd1,
        0x6fa87e4f, 0xfe2ce6e0, 0xa3014314, 0x4e0811a1, 0xf7537e82, 0xbd3af235, 0x2ad7d2bb,
        0xeb86d391,
    ]
};

pub struct Md5 {
    state: [u32; 4],
    buf: [u8; 64],
    buf_len: usize,
    total: u64,
}

impl Default for Md5 {
    fn default() -> Self {
        Self::new()
    }
}

impl Md5 {
    pub fn new() -> Md5 {
        Md5 {
            state: [0x67452301, 0xefcdab89, 0x98badcfe, 0x10325476],
            buf: [0; 64],
            buf_len: 0,
            total: 0,
        }
    }

    fn block(state: &mut [u32; 4], chunk: &[u8]) {
        let mut m = [0u32; 16];
        for (i, w) in m.iter_mut().enumerate() {
            *w = u32::from_le_bytes(chunk[4 * i..4 * i + 4].try_into().unwrap());
        }
        let [mut a, mut b, mut c, mut d] = *state;
        // Fully unrolled: round i updates the register that is `a` in RFC
        // 1321's notation; the four registers rotate roles instead of values.
        macro_rules! step {
            ($f:expr, $a:ident, $b:ident, $c:ident, $d:ident, $i:expr, $g:expr) => {
                $a = $b.wrapping_add(
                    $a.wrapping_add($f($b, $c, $d))
                        .wrapping_add(K[$i])
                        .wrapping_add(m[$g])
                        .rotate_left(S[$i]),
                );
            };
        }
        let f = |b: u32, c: u32, d: u32| (b & c) | (!b & d);
        let g = |b: u32, c: u32, d: u32| (d & b) | (!d & c);
        let h = |b: u32, c: u32, d: u32| b ^ c ^ d;
        let ii = |b: u32, c: u32, d: u32| c ^ (b | !d);
        macro_rules! round {
            ($fun:expr, $base:expr, $gi:expr) => {
                step!($fun, a, b, c, d, $base, $gi($base));
                step!($fun, d, a, b, c, $base + 1, $gi($base + 1));
                step!($fun, c, d, a, b, $base + 2, $gi($base + 2));
                step!($fun, b, c, d, a, $base + 3, $gi($base + 3));
            };
        }
        const fn g0(i: usize) -> usize {
            i
        }
        const fn g1(i: usize) -> usize {
            (5 * i + 1) % 16
        }
        const fn g2(i: usize) -> usize {
            (3 * i + 5) % 16
        }
        const fn g3(i: usize) -> usize {
            (7 * i) % 16
        }
        round!(f, 0, g0);
        round!(f, 4, g0);
        round!(f, 8, g0);
        round!(f, 12, g0);
        round!(g, 16, g1);
        round!(g, 20, g1);
        round!(g, 24, g1);
        round!(g, 28, g1);
        round!(h, 32, g2);
        round!(h, 36, g2);
        round!(h, 40, g2);
        round!(h, 44, g2);
        round!(ii, 48, g3);
        round!(ii, 52, g3);
        round!(ii, 56, g3);
        round!(ii, 60, g3);
        state[0] = state[0].wrapping_add(a);
        state[1] = state[1].wrapping_add(b);
        state[2] = state[2].wrapping_add(c);
        state[3] = state[3].wrapping_add(d);
    }

    pub fn update(&mut self, mut data: &[u8]) {
        self.total += data.len() as u64;
        if self.buf_len > 0 {
            let take = (64 - self.buf_len).min(data.len());
            self.buf[self.buf_len..self.buf_len + take].copy_from_slice(&data[..take]);
            self.buf_len += take;
            data = &data[take..];
            if self.buf_len < 64 {
                return;
            }
            let buf = self.buf;
            Self::block(&mut self.state, &buf);
            self.buf_len = 0;
        }
        let (blocks, rest) = data.as_chunks::<64>();
        for c in blocks {
            Self::block(&mut self.state, c);
        }
        self.buf[..rest.len()].copy_from_slice(rest);
        self.buf_len = rest.len();
    }

    /// Feed `data` upper-cased (CRAM reference checksums are over upper case).
    pub fn update_upper(&mut self, data: &[u8]) {
        let mut tmp = [0u8; 4096];
        for c in data.chunks(tmp.len()) {
            for (t, &b) in tmp.iter_mut().zip(c) {
                *t = b.to_ascii_uppercase();
            }
            self.update(&tmp[..c.len()]);
        }
    }

    pub fn finish(mut self) -> [u8; 16] {
        let bits = self.total.wrapping_mul(8);
        self.update(&[0x80]);
        while self.buf_len != 56 {
            self.update(&[0]);
        }
        self.update(&bits.to_le_bytes());
        let mut out = [0u8; 16];
        for (i, w) in self.state.iter().enumerate() {
            out[4 * i..4 * i + 4].copy_from_slice(&w.to_le_bytes());
        }
        out
    }
}

pub fn hex(d: &[u8; 16]) -> String {
    d.iter().map(|b| format!("{b:02x}")).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn md5(s: &[u8]) -> String {
        let mut m = Md5::new();
        m.update(s);
        hex(&m.finish())
    }

    #[test]
    fn rfc1321_vectors() {
        assert_eq!(md5(b""), "d41d8cd98f00b204e9800998ecf8427e");
        assert_eq!(md5(b"a"), "0cc175b9c0f1b6a831c399e269772661");
        assert_eq!(md5(b"abc"), "900150983cd24fb0d6963f7d28e17f72");
        assert_eq!(md5(b"message digest"), "f96b697d7cb7938d525a2f31aaf161d0");
        assert_eq!(
            md5(
                b"12345678901234567890123456789012345678901234567890123456789012345678901234567890"
            ),
            "57edf4a22be3c955ac49da2e2107b67a"
        );
    }

    #[test]
    fn chunking_and_upper() {
        let data: Vec<u8> = (0..10_000u32).map(|i| b"acgtN"[i as usize % 5]).collect();
        let mut a = Md5::new();
        for c in data.chunks(37) {
            a.update_upper(c);
        }
        let upper: Vec<u8> = data.iter().map(|b| b.to_ascii_uppercase()).collect();
        let mut b = Md5::new();
        b.update(&upper);
        assert_eq!(a.finish(), b.finish());
    }
}
