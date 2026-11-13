/// gzip CRC-32 (reflected 0xEDB88320) with the ARMv8 CRC32 instructions.
///
/// `state` is the running register as in `nomics_io::crc32::Crc32`: the
/// caller starts from `!0` and inverts at the end. Returns `None` when the
/// CPU has no CRC32 extension (or is not aarch64); the caller then uses its
/// portable table code, which gives the same value.
pub fn crc32_update(state: u32, data: &[u8]) -> Option<u32> {
    #[cfg(all(target_arch = "aarch64", target_endian = "little"))]
    if std::arch::is_aarch64_feature_detected!("crc") {
        // SAFETY: the CRC extension was detected above; the kernel only
        // calls CRC32 intrinsics on values read through checked slices.
        return Some(unsafe { crc32_arm(state, data) });
    }
    let _ = (state, data);
    None
}

#[cfg(all(target_arch = "aarch64", target_endian = "little"))]
#[target_feature(enable = "crc")]
fn crc32_arm(mut c: u32, data: &[u8]) -> u32 {
    use core::arch::aarch64::{__crc32b, __crc32d};
    let (words, rest) = data.as_chunks::<8>();
    for w in words {
        c = __crc32d(c, u64::from_le_bytes(*w));
    }
    for &b in rest {
        c = __crc32b(c, b);
    }
    c
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Bitwise reference.
    fn portable(mut c: u32, data: &[u8]) -> u32 {
        for &b in data {
            c ^= b as u32;
            for _ in 0..8 {
                c = if c & 1 != 0 {
                    0xEDB8_8320 ^ (c >> 1)
                } else {
                    c >> 1
                };
            }
        }
        c
    }

    #[test]
    fn matches_portable() {
        let data: Vec<u8> = (0..5000u32)
            .map(|i| (i.wrapping_mul(2_654_435_761) >> 13) as u8)
            .collect();
        for (start, end) in [(0, 0), (0, 1), (3, 10), (1, 4096), (0, 5000), (7, 4999)] {
            let d = &data[start..end];
            for init in [!0u32, 0, 0x1234_5678] {
                if let Some(c) = crc32_update(init, d) {
                    assert_eq!(c, portable(init, d), "{start}..{end}");
                }
            }
        }
        assert_eq!(
            crc32_update(!0, b"123456789")
                .map(|c| !c)
                .unwrap_or(0xCBF4_3926),
            0xCBF4_3926
        );
    }
}
