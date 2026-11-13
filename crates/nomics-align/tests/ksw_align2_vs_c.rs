//! Differential tests for the striped `ksw_align2` port.

use std::{
    env, fs,
    path::{Path, PathBuf},
    process::Command,
    sync::atomic::{AtomicUsize, Ordering},
};

use nomics_align::ksw::{KSW_XBYTE, KSW_XSTART, KSW_XSTOP, KSW_XSUBO, Kswr, Profile, align2};

static NEXT: AtomicUsize = AtomicUsize::new(0);

#[derive(Debug)]
struct Case {
    query: Vec<u8>,
    target: Vec<u8>,
    mat: Vec<i8>,
    o_del: i32,
    e_del: i32,
    o_ins: i32,
    e_ins: i32,
    xtra: i32,
}

#[derive(Default)]
struct Stats {
    u8_path: usize,
    i16_path: usize,
    byte_overflow: usize,
    reverse_pass: usize,
    score2_found: usize,
    score2_missing: usize,
}

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u32 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0 as u32
    }

    fn below(&mut self, n: u32) -> u32 {
        self.next() % n
    }
}

fn bwa_mat() -> Vec<i8> {
    (0..5)
        .flat_map(|i| {
            (0..5).map(move |j| {
                if i == 4 || j == 4 {
                    -1
                } else if i == j {
                    1
                } else {
                    -4
                }
            })
        })
        .collect()
}

fn random_mat(rng: &mut Rng) -> Vec<i8> {
    let mut mat = Vec::with_capacity(25);
    for i in 0..5 {
        for j in 0..5 {
            mat.push(if i == j {
                1 + rng.below(5) as i8
            } else {
                -(1 + rng.below(8) as i8)
            });
        }
    }
    mat
}

fn random_sequence(rng: &mut Rng, len: usize) -> Vec<u8> {
    (0..len)
        .map(|_| {
            if rng.below(29) == 0 {
                4
            } else {
                rng.below(4) as u8
            }
        })
        .collect()
}

fn add_flanks(rng: &mut Rng, middle: &[u8], long_window: bool) -> Vec<u8> {
    let flank_limit = if long_window { 601 } else { 25 };
    let left = rng.below(flank_limit) as usize;
    let right = rng.below(flank_limit) as usize;
    let mut target = random_sequence(rng, left);
    target.extend_from_slice(middle);
    target.extend(random_sequence(rng, right));
    target.truncate(1500);
    target
}

fn make_case(rng: &mut Rng, number: usize) -> Case {
    let qlen = if number.is_multiple_of(127) {
        15 + rng.below(286) as usize
    } else {
        15 + rng.below(21) as usize
    };
    let mut query = random_sequence(rng, qlen);
    if number.is_multiple_of(41) {
        let motif = [0, 1, 2, 3, 4, 2, 1];
        for (i, base) in query.iter_mut().enumerate() {
            *base = motif[i % motif.len()];
        }
    }

    let long_window = number.is_multiple_of(997);
    let target = match number % 8 {
        0 => add_flanks(rng, &query, long_window),
        1 => {
            let mut copies = query.clone();
            let separation = rng.below(24) as usize;
            copies.extend(random_sequence(rng, separation));
            copies.extend_from_slice(&query);
            add_flanks(rng, &copies, long_window)
        }
        2 => {
            let mut copy = query.clone();
            for _ in 0..1 + rng.below(7) {
                let at = rng.below(copy.len() as u32) as usize;
                copy[at] = rng.below(5) as u8;
            }
            add_flanks(rng, &copy, long_window)
        }
        3 => {
            let mut copy = query.clone();
            if rng.below(2) == 0 {
                let at = rng.below((copy.len() + 1) as u32) as usize;
                let insertion_len = 1 + rng.below(12) as usize;
                let insertion = random_sequence(rng, insertion_len);
                copy.splice(at..at, insertion);
            } else {
                let count = (1 + rng.below(12) as usize).min(copy.len() - 1);
                let at = rng.below((copy.len() - count + 1) as u32) as usize;
                copy.drain(at..at + count);
            }
            add_flanks(rng, &copy, long_window)
        }
        4 => {
            let mut copy = query.clone();
            for _ in 0..1 + rng.below(5) {
                let at = rng.below(copy.len() as u32) as usize;
                copy[at] = 4;
            }
            add_flanks(rng, &copy, long_window)
        }
        5 => {
            let len = if long_window {
                1500
            } else {
                10 + rng.below(90) as usize
            };
            random_sequence(rng, len)
        }
        6 => {
            let len = rng.below(qlen as u32) as usize;
            random_sequence(rng, len)
        }
        _ => {
            let len = rng.below(9) as usize;
            random_sequence(rng, len)
        }
    };

    if number.is_multiple_of(9) {
        let mat = random_mat(rng);
        let byte = rng.below(2) == 0;
        let mut xtra = rng.below(31) as i32;
        if byte {
            xtra |= KSW_XBYTE;
        }
        if rng.below(2) == 0 {
            xtra |= KSW_XSUBO;
        }
        if rng.below(3) == 0 {
            xtra |= KSW_XSTOP;
        }
        // Avoid BWA's documented byte-overflow/XSTART invalid combination.
        if !byte && rng.below(2) == 0 {
            xtra |= KSW_XSTART;
        }
        Case {
            query,
            target,
            mat,
            o_del: 1 + rng.below(20) as i32,
            e_del: 1 + rng.below(6) as i32,
            o_ins: 1 + rng.below(20) as i32,
            e_ins: 1 + rng.below(6) as i32,
            xtra,
        }
    } else {
        // This is exactly mem_matesw's xtra construction for a=1 and
        // min_seed_len=19.
        let xtra = KSW_XSUBO | KSW_XSTART | if qlen < 250 { KSW_XBYTE } else { 0 } | 19;
        Case {
            query,
            target,
            mat: bwa_mat(),
            o_del: 6,
            e_del: 1,
            o_ins: 6,
            e_ins: 1,
            xtra,
        }
    }
}

fn special_cases() -> Vec<Case> {
    let mat = bwa_mat();
    let tandem_query: Vec<u8> = (0..63).map(|i| (i % 4) as u8).collect();
    let mut tandem_target = vec![4, 4, 3, 2];
    tandem_target.extend_from_slice(&tandem_query);
    tandem_target.extend_from_slice(&[4; 21]);
    tandem_target.extend_from_slice(&tandem_query);
    vec![
        Case {
            query: vec![0; 15],
            target: vec![0; 15],
            mat: mat.clone(),
            o_del: 6,
            e_del: 1,
            o_ins: 6,
            e_ins: 1,
            xtra: KSW_XBYTE | KSW_XSUBO | KSW_XSTART | 15,
        },
        Case {
            query: tandem_query,
            target: tandem_target,
            mat: mat.clone(),
            o_del: 6,
            e_del: 1,
            o_ins: 6,
            e_ins: 1,
            xtra: KSW_XBYTE | KSW_XSUBO | KSW_XSTART | 19,
        },
        Case {
            query: vec![0; 300],
            target: vec![0; 300],
            mat: mat.clone(),
            o_del: 6,
            e_del: 1,
            o_ins: 6,
            e_ins: 1,
            xtra: KSW_XBYTE | KSW_XSUBO,
        },
        Case {
            query: vec![0, 1, 2, 3, 4, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1],
            target: Vec::new(),
            mat,
            o_del: 6,
            e_del: 1,
            o_ins: 6,
            e_ins: 1,
            xtra: KSW_XSTART,
        },
    ]
}

fn i32_bytes(out: &mut Vec<u8>, value: i32) {
    out.extend_from_slice(&value.to_ne_bytes());
}

fn read_i32(input: &[u8], cursor: &mut usize) -> i32 {
    let bytes: [u8; 4] = input[*cursor..*cursor + 4].try_into().unwrap();
    *cursor += 4;
    i32::from_ne_bytes(bytes)
}

fn compile_oracle(dir: &Path) -> PathBuf {
    let manifest = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let c = manifest.join("tests/c");
    let executable = dir.join("ksw_align2_oracle");
    let status = Command::new(env::var("CC").unwrap_or_else(|_| "cc".into()))
        .args(["-std=c99", "-O3", "-I"])
        .arg(&c)
        .arg(c.join("ksw.c"))
        .arg(c.join("harness_ksw_align2/ksw_align2_oracle.c"))
        .arg("-o")
        .arg(&executable)
        .status()
        .expect("run C compiler");
    assert!(status.success(), "C oracle compiler failed: {status}");
    executable
}

fn oracle_results(cases: &[Case], root: &Path) -> Vec<Kswr> {
    let input = root.join("input.bin");
    let output = root.join("output.bin");
    let mut bytes = Vec::new();
    bytes.extend_from_slice(&(cases.len() as u32).to_ne_bytes());
    for case in cases {
        for value in [
            case.query.len() as i32,
            case.target.len() as i32,
            5,
            case.o_del,
            case.e_del,
            case.o_ins,
            case.e_ins,
            case.xtra,
        ] {
            i32_bytes(&mut bytes, value);
        }
        bytes.extend(case.mat.iter().map(|&score| score as u8));
        bytes.extend_from_slice(&case.query);
        bytes.extend_from_slice(&case.target);
    }
    fs::write(&input, bytes).unwrap();
    let oracle = compile_oracle(root);
    assert!(
        Command::new(oracle)
            .arg(&input)
            .arg(&output)
            .status()
            .unwrap()
            .success()
    );
    let actual = fs::read(output).unwrap();
    assert_eq!(actual.len(), cases.len() * 7 * size_of::<i32>());
    let mut cursor = 0;
    (0..cases.len())
        .map(|_| Kswr {
            score: read_i32(&actual, &mut cursor),
            te: read_i32(&actual, &mut cursor),
            qe: read_i32(&actual, &mut cursor),
            score2: read_i32(&actual, &mut cursor),
            te2: read_i32(&actual, &mut cursor),
            tb: read_i32(&actual, &mut cursor),
            qb: read_i32(&actual, &mut cursor),
        })
        .collect()
}

fn run_differential(count: usize) {
    let root = env::temp_dir().join(format!(
        "nomics-ksw-align2-{}-{}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ));
    fs::create_dir_all(&root).unwrap();
    let mut rng = Rng(0x82a2_b175_229d_6a5b);
    let mut cases = special_cases();
    cases.extend((cases.len()..count).map(|number| make_case(&mut rng, number)));
    let expected = oracle_results(&cases, &root);
    let mut stats = Stats::default();
    for (case_no, (case, expected)) in cases.iter().zip(expected).enumerate() {
        let got = align2(
            case.query.len() as i32,
            &case.query,
            case.target.len() as i32,
            &case.target,
            5,
            &case.mat,
            case.o_del,
            case.e_del,
            case.o_ins,
            case.e_ins,
            case.xtra,
            None,
        );
        assert_eq!(got, expected, "case {case_no}: {case:#?}");
        if case.xtra & KSW_XBYTE != 0 {
            stats.u8_path += 1;
            stats.byte_overflow += usize::from(got.score == 255);
        } else {
            stats.i16_path += 1;
        }
        let reverse_skipped_for_threshold =
            case.xtra & KSW_XSUBO != 0 && got.score < case.xtra & 0xffff;
        if case.xtra & KSW_XSTART != 0 && !reverse_skipped_for_threshold {
            stats.reverse_pass += 1;
        }
        if got.score2 >= 0 {
            stats.score2_found += 1;
        } else {
            stats.score2_missing += 1;
        }
    }
    eprintln!(
        "ksw_align2 cases={count} u8={} i16={} byte_overflow={} reverse={} score2_found={} score2_minus_one={}",
        stats.u8_path,
        stats.i16_path,
        stats.byte_overflow,
        stats.reverse_pass,
        stats.score2_found,
        stats.score2_missing
    );
    let _ = fs::remove_dir_all(root);
}

#[test]
fn align2_matches_c_debug_suite() {
    run_differential(4_096);
}

#[test]
fn reusable_profile_matches_fresh_profile() {
    let query: Vec<u8> = (0..73).map(|i| (i * 3 % 5) as u8).collect();
    let mat = bwa_mat();
    let mut profile = Profile::new(1, query.len() as i32, &query, 5, &mat);
    for offset in 0..12 {
        let mut target = vec![4; offset];
        target.extend_from_slice(&query);
        target.extend_from_slice(&query[offset..offset + 20]);
        let xtra = KSW_XBYTE | KSW_XSUBO | KSW_XSTART | 19;
        let fresh = align2(
            query.len() as i32,
            &query,
            target.len() as i32,
            &target,
            5,
            &mat,
            6,
            1,
            6,
            1,
            xtra,
            None,
        );
        let reused = align2(
            query.len() as i32,
            &query,
            target.len() as i32,
            &target,
            5,
            &mat,
            6,
            1,
            6,
            1,
            xtra,
            Some(&mut profile),
        );
        assert_eq!(reused, fresh);
    }
}

#[test]
#[ignore = "200k-case C differential suite; run with --release --ignored"]
fn align2_matches_c_heavy_suite() {
    run_differential(200_000);
}
