//! `ZlibFast` against classic zlib 1.3 output. Fixtures come from
//! `gen_zlib_fixtures.py`, run on the Spark (the zlib the GATK JVM loaded);
//! `synth` below must stay identical to `gen` there.

use std::path::PathBuf;

use nomics_io::deflate::{Finish, HTSJDK_OUTPUT_CAPACITY, ZlibFast};
use nomics_io::inflate::inflate;

fn synth(seed: u64, n: usize, alphabet: &[u8], copies: bool) -> Vec<u8> {
    let mut x = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15).wrapping_add(1);
    let mut next = || {
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        x.wrapping_mul(0x2545_F491_4F6C_DD1D)
    };
    let mut out = Vec::with_capacity(n);
    while out.len() < n {
        let r = next();
        if copies && out.len() > 16 && r & 7 == 0 {
            let dist = 1 + (next() % out.len().min(40000) as u64) as usize;
            let length = 3 + (next() % 300) as usize;
            for _ in 0..length {
                if out.len() >= n {
                    break;
                }
                out.push(out[out.len() - dist]);
            }
        } else {
            out.push(alphabet[((r >> 8) % alphabet.len() as u64) as usize]);
        }
    }
    out
}

fn cases() -> Vec<(&'static str, Vec<u8>)> {
    let dna = b"ACGT";
    let qual: Vec<u8> = (33..75).collect();
    let all: Vec<u8> = (0..=255).collect();
    let mut dna_n = dna.to_vec();
    dna_n.push(b'N');
    let mut mixed = synth(6, 30000, &qual, true);
    mixed.extend(synth(7, 35498, dna, true));
    vec![
        ("empty", vec![]),
        ("one", b"x".to_vec()),
        ("hello", b"hello".to_vec()),
        ("dna_1k", synth(1, 1000, dna, true)),
        ("qual_20k", synth(2, 20000, &qual, true)),
        ("random_5k", synth(3, 5000, &all, true)),
        ("random_65498", synth(4, 65498, &all, true)),
        ("dna_65498", synth(5, 65498, &dna_n, true)),
        ("mixed_65498", mixed),
        ("zeros_65498", vec![0; 65498]),
        ("noise_65498", synth(8, 65498, &all, false)),
        ("noise_40000", synth(9, 40000, &all, false)),
    ]
}

fn fixture(name: &str) -> Vec<u8> {
    let p = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/zlib")
        .join(name);
    std::fs::read(&p).unwrap_or_else(|e| panic!("{}: {e}", p.display()))
}

/// `deflate_slow` (lazy matching). Fixtures from `gen_zlib_fixtures.py
/// --lazy` (zlib 1.2.13; deflate_slow is unchanged through 1.3.1).
#[test]
fn fresh_stream_matches_zlib_lazy_levels() {
    for (name, input) in cases() {
        for level in [4, 5, 6, 9] {
            let expected = fixture(&format!("{name}.l{level}.deflate"));
            let mut out = Vec::new();
            nomics_io::deflate::ZlibDeflater::new(level).compress(&input, &mut out);
            assert!(
                out == expected,
                "{name} level {level}: {} bytes vs zlib {}",
                out.len(),
                expected.len()
            );
            let mut back = Vec::new();
            inflate(&out, &mut back, usize::MAX).unwrap();
            assert_eq!(back, input, "{name} l{level} round trip");
        }
    }
}

#[test]
fn fresh_stream_matches_zlib_levels_1_to_3() {
    for (name, input) in cases() {
        for level in 1..=3 {
            let expected = fixture(&format!("{name}.l{level}.deflate"));
            let mut z = ZlibFast::new(level);
            let mut out = Vec::new();
            assert_eq!(
                z.deflate_finish(&input, &mut out, usize::MAX),
                Finish::Done,
                "{name} l{level}"
            );
            assert!(
                out == expected,
                "{name} level {level}: {} bytes vs zlib {}",
                out.len(),
                expected.len()
            );
            let mut back = Vec::new();
            inflate(&out, &mut back, usize::MAX).unwrap();
            assert_eq!(back, input, "{name} l{level} round trip");
        }
    }
}

#[test]
fn reused_stream_matches_zlib_deflate_reset_sequence() {
    // Mirrors htsjdk: one deflater, deflateReset per block, 65,518-byte output.
    let c = cases();
    let get = |n: &str| c.iter().find(|(k, _)| *k == n).unwrap().1.clone();
    let seq = [
        "dna_65498",
        "mixed_65498",
        "dna_1k",
        "noise_65498",
        "random_65498",
        "qual_20k",
        "noise_40000",
        "dna_65498",
    ];
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/zlib");
    let mut z = ZlibFast::new(2);
    for (i, name) in seq.iter().enumerate() {
        let done = dir.join(format!("seq{i}.done.deflate"));
        let (expected_finish, expected) = if done.exists() {
            (Finish::Done, std::fs::read(done).unwrap())
        } else {
            (
                Finish::OutOfSpace,
                fixture(&format!("seq{i}.outofspace.deflate")),
            )
        };
        let mut out = Vec::new();
        let finish = z.deflate_finish(&get(name), &mut out, HTSJDK_OUTPUT_CAPACITY);
        assert_eq!(finish, expected_finish, "block {i} ({name})");
        // zlib hands back at most `capacity` bytes when it runs out of space.
        let got = &out[..out.len().min(HTSJDK_OUTPUT_CAPACITY)];
        assert!(
            got == expected,
            "block {i} ({name}): {} bytes vs zlib {}",
            got.len(),
            expected.len()
        );
    }
}
