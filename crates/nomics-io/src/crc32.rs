//! CRC-32 (IEEE 802.3, reflected polynomial 0xEDB88320) as used by gzip/BGZF.
//! ARMv8 CRC32 instructions where available, else slice-by-8 tables.

const POLY: u32 = 0xEDB8_8320;

const TABLES: [[u32; 256]; 8] = build_tables();

const fn build_tables() -> [[u32; 256]; 8] {
    let mut t = [[0u32; 256]; 8];
    let mut i = 0;
    while i < 256 {
        let mut c = i as u32;
        let mut k = 0;
        while k < 8 {
            c = if c & 1 != 0 { POLY ^ (c >> 1) } else { c >> 1 };
            k += 1;
        }
        t[0][i] = c;
        i += 1;
    }
    let mut s = 1;
    while s < 8 {
        let mut i = 0;
        while i < 256 {
            let prev = t[s - 1][i];
            t[s][i] = (prev >> 8) ^ t[0][(prev & 0xff) as usize];
            i += 1;
        }
        s += 1;
    }
    t
}

/// Incremental CRC-32. `Crc32::new().update(a).update(b).finish()` equals the
/// CRC of `a ++ b`.
#[derive(Clone, Copy, Debug)]
pub struct Crc32(u32);

impl Default for Crc32 {
    fn default() -> Self {
        Self::new()
    }
}

impl Crc32 {
    pub const fn new() -> Self {
        Crc32(!0)
    }

    /// Uses the ARMv8 CRC32 instructions when the CPU has them
    /// (`nomics_simd::crc32_update`), else [`Self::update_portable`].
    pub fn update(self, data: &[u8]) -> Self {
        match nomics_simd::crc32_update(self.0, data) {
            Some(c) => Crc32(c),
            None => self.update_portable(data),
        }
    }

    /// Slice-by-8 tables; the reference for the hardware path.
    pub fn update_portable(mut self, data: &[u8]) -> Self {
        let mut c = self.0;
        let (chunks, rest) = data.as_chunks::<8>();
        for b in chunks {
            let lo = c ^ u32::from_le_bytes([b[0], b[1], b[2], b[3]]);
            let hi = u32::from_le_bytes([b[4], b[5], b[6], b[7]]);
            c = TABLES[7][(lo & 0xff) as usize]
                ^ TABLES[6][((lo >> 8) & 0xff) as usize]
                ^ TABLES[5][((lo >> 16) & 0xff) as usize]
                ^ TABLES[4][(lo >> 24) as usize]
                ^ TABLES[3][(hi & 0xff) as usize]
                ^ TABLES[2][((hi >> 8) & 0xff) as usize]
                ^ TABLES[1][((hi >> 16) & 0xff) as usize]
                ^ TABLES[0][(hi >> 24) as usize];
        }
        for &b in rest {
            c = (c >> 8) ^ TABLES[0][((c ^ b as u32) & 0xff) as usize];
        }
        self.0 = c;
        self
    }

    pub fn finish(self) -> u32 {
        !self.0
    }
}

pub fn crc32(data: &[u8]) -> u32 {
    Crc32::new().update(data).finish()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn known_vectors() {
        assert_eq!(crc32(b""), 0);
        assert_eq!(crc32(b"123456789"), 0xCBF4_3926);
        assert_eq!(
            crc32(b"The quick brown fox jumps over the lazy dog"),
            0x414F_A339
        );
    }

    #[test]
    fn hardware_matches_portable() {
        let data: Vec<u8> = (0..70_000u32)
            .map(|i| (i.wrapping_mul(2654435761) >> 11) as u8)
            .collect();
        for (a, b) in [(0, 0), (0, 1), (1, 8), (3, 65501), (0, 70_000), (9, 69_999)] {
            let d = &data[a..b];
            assert_eq!(
                Crc32::new().update(d).finish(),
                Crc32::new().update_portable(d).finish(),
                "{a}..{b}"
            );
        }
    }

    #[test]
    fn incremental_matches_one_shot() {
        let data: Vec<u8> = (0..10_000u32)
            .map(|i| (i.wrapping_mul(2654435761) >> 13) as u8)
            .collect();
        for split in [0, 1, 7, 8, 9, 4096, 9999, 10_000] {
            let (a, b) = data.split_at(split);
            assert_eq!(Crc32::new().update(a).update(b).finish(), crc32(&data));
        }
    }
}
