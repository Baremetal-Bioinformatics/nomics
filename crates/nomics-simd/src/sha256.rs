const K: [u32; 64] = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
];

/// Portable SHA-256 compression, also available for differential tests and benchmarks.
///
/// # Panics
/// Panics if `blocks.len()` is not a multiple of 64. No padding is added.
pub fn sha256_compress_portable(state: &mut [u32; 8], blocks: &[u8]) {
    assert!(
        blocks.len().is_multiple_of(64),
        "whole SHA-256 blocks required"
    );
    for block in blocks.as_chunks::<64>().0 {
        let mut w = [0u32; 64];
        for (word, bytes) in w[..16].iter_mut().zip(block.as_chunks::<4>().0) {
            *word = u32::from_be_bytes(*bytes);
        }
        for i in 16..64 {
            let s0 = w[i - 15].rotate_right(7) ^ w[i - 15].rotate_right(18) ^ (w[i - 15] >> 3);
            let s1 = w[i - 2].rotate_right(17) ^ w[i - 2].rotate_right(19) ^ (w[i - 2] >> 10);
            w[i] = w[i - 16]
                .wrapping_add(s0)
                .wrapping_add(w[i - 7])
                .wrapping_add(s1);
        }
        let [mut a, mut b, mut c, mut d, mut e, mut f, mut g, mut h] = *state;
        for i in 0..64 {
            let s1 = e.rotate_right(6) ^ e.rotate_right(11) ^ e.rotate_right(25);
            let ch = (e & f) ^ (!e & g);
            let t1 = h
                .wrapping_add(s1)
                .wrapping_add(ch)
                .wrapping_add(K[i])
                .wrapping_add(w[i]);
            let s0 = a.rotate_right(2) ^ a.rotate_right(13) ^ a.rotate_right(22);
            let maj = (a & b) ^ (a & c) ^ (b & c);
            let t2 = s0.wrapping_add(maj);
            h = g;
            g = f;
            f = e;
            e = d.wrapping_add(t1);
            d = c;
            c = b;
            b = a;
            a = t1.wrapping_add(t2);
        }
        for (state, value) in state.iter_mut().zip([a, b, c, d, e, f, g, h]) {
            *state = state.wrapping_add(value);
        }
    }
}

/// Compress whole SHA-256 blocks, using ARMv8 SHA2 when available.
/// `state` contains the eight chaining words in host byte order; input words
/// are big-endian. Empty input leaves `state` unchanged. No padding is added.
///
/// # Panics
/// Panics if `blocks.len()` is not a multiple of 64.
pub fn sha256_compress(state: &mut [u32; 8], blocks: &[u8]) {
    assert!(
        blocks.len().is_multiple_of(64),
        "whole SHA-256 blocks required"
    );
    #[cfg(all(target_arch = "aarch64", target_endian = "little"))]
    if std::arch::is_aarch64_feature_detected!("sha2") {
        // SAFETY: SHA2 support was detected above. All loads/stores in the
        // kernel use whole arrays obtained from checked slices.
        unsafe { compress_sha2(state, blocks) };
        return;
    }
    sha256_compress_portable(state, blocks);
}

#[cfg(all(target_arch = "aarch64", target_endian = "little"))]
#[target_feature(enable = "sha2")]
fn compress_sha2(state: &mut [u32; 8], blocks: &[u8]) {
    use core::arch::aarch64::*;

    // SAFETY: each pointer refers to an entire array of the loaded/stored size;
    // NEON loads/stores do not require alignment. SHA2 is enabled for this fn.
    let load = |v: &[u32; 4]| unsafe { vld1q_u32(v.as_ptr()) };
    let load_be = |v: &[u8; 16]| unsafe { vreinterpretq_u32_u8(vrev32q_u8(vld1q_u8(v.as_ptr()))) };
    let store = |v: &mut [u32; 4], x| unsafe { vst1q_u32(v.as_mut_ptr(), x) };
    let mut abcd = load(state[..4].try_into().unwrap());
    let mut efgh = load(state[4..].try_into().unwrap());
    let constants = K.as_chunks::<4>().0;
    for block in blocks.as_chunks::<64>().0 {
        let saved_abcd = abcd;
        let saved_efgh = efgh;
        let words = block.as_chunks::<16>().0;
        let mut w0 = load_be(&words[0]);
        let mut w1 = load_be(&words[1]);
        let mut w2 = load_be(&words[2]);
        let mut w3 = load_be(&words[3]);

        // Each hash instruction computes four rounds. H2 needs the OLD abcd.
        // The schedule rotates through four vectors; update a vector only
        // after its four words have been consumed. Unroll to keep the schedule
        // and constants in registers, with no indexing in the round loop.
        macro_rules! rounds {
            ($w:ident, $next:ident, $third:ident, $last:ident, $i:expr) => {{
                let wk = vaddq_u32($w, load(&constants[$i]));
                let old_abcd = abcd;
                abcd = vsha256hq_u32(abcd, efgh, wk);
                efgh = vsha256h2q_u32(efgh, old_abcd, wk);
                if $i < 12 {
                    $w = vsha256su1q_u32(vsha256su0q_u32($w, $next), $third, $last);
                }
            }};
        }
        rounds!(w0, w1, w2, w3, 0);
        rounds!(w1, w2, w3, w0, 1);
        rounds!(w2, w3, w0, w1, 2);
        rounds!(w3, w0, w1, w2, 3);
        rounds!(w0, w1, w2, w3, 4);
        rounds!(w1, w2, w3, w0, 5);
        rounds!(w2, w3, w0, w1, 6);
        rounds!(w3, w0, w1, w2, 7);
        rounds!(w0, w1, w2, w3, 8);
        rounds!(w1, w2, w3, w0, 9);
        rounds!(w2, w3, w0, w1, 10);
        rounds!(w3, w0, w1, w2, 11);
        rounds!(w0, w1, w2, w3, 12);
        rounds!(w1, w2, w3, w0, 13);
        rounds!(w2, w3, w0, w1, 14);
        rounds!(w3, w0, w1, w2, 15);
        abcd = vaddq_u32(abcd, saved_abcd);
        efgh = vaddq_u32(efgh, saved_efgh);
    }
    store((&mut state[..4]).try_into().unwrap(), abcd);
    store((&mut state[4..]).try_into().unwrap(), efgh);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dispatch_matches_portable_chaining_states() {
        let mut seed = 0x1234_5678_9abc_def0u64;
        let mut random = || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed as u32
        };
        for blocks in (0..=64).chain([127, 128, 129, 1024, 4096]) {
            // Deliberately unaligned, with an arbitrary initial chaining state.
            let mut input = vec![0; blocks * 64 + 1];
            for byte in &mut input {
                *byte = random() as u8;
            }
            let initial = std::array::from_fn(|_| random());
            let mut expected = initial;
            sha256_compress_portable(&mut expected, &input[1..]);
            let mut actual = initial;
            sha256_compress(&mut actual, &input[1..]);
            assert_eq!(actual, expected, "{blocks} blocks");
            let mut split = initial;
            for block in input[1..].chunks(64) {
                sha256_compress(&mut split, block);
            }
            assert_eq!(split, expected, "split {blocks} blocks");
            if blocks == 0 {
                assert_eq!(actual, initial);
            }
        }
    }

    #[test]
    #[should_panic(expected = "whole SHA-256 blocks required")]
    fn dispatch_rejects_partial_block() {
        sha256_compress(&mut [0; 8], &[0; 65]);
    }

    #[test]
    #[should_panic(expected = "whole SHA-256 blocks required")]
    fn portable_rejects_partial_block() {
        sha256_compress_portable(&mut [0; 8], &[0; 63]);
    }
}
