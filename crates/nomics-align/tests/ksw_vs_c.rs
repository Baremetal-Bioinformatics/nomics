//! Differential testing against the unmodified BWA 0.7.19 ksw.c.

use std::{
    env, fs,
    path::{Path, PathBuf},
    process::Command,
    sync::atomic::{AtomicUsize, Ordering},
};

use nomics_align::ksw::{Extend, extend2, global2};

static NEXT: AtomicUsize = AtomicUsize::new(0);

#[derive(Clone)]
struct Case {
    extend: bool,
    query: Vec<u8>,
    target: Vec<u8>,
    mat: Vec<i8>,
    o_del: i32,
    e_del: i32,
    o_ins: i32,
    e_ins: i32,
    w: i32,
    end_bonus: i32,
    zdrop: i32,
    h0: i32,
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

fn random_base(rng: &mut Rng) -> u8 {
    rng.below(5) as u8
}

// The distribution intentionally makes the large required case count cheap,
// but injects 1..300bp independent, repeat, substitution and indel examples.
fn make_case(rng: &mut Rng, extend: bool, number: usize) -> Case {
    let long = number.is_multiple_of(257);
    let qlen = if long {
        1 + rng.below(300) as usize
    } else {
        rng.below(17) as usize
    };
    let mut query: Vec<u8> = (0..qlen).map(|_| random_base(rng)).collect();
    if number.is_multiple_of(31) && !query.is_empty() {
        // tandem repeat
        let motif = [0, 1, 2, 3];
        for (i, base) in query.iter_mut().enumerate() {
            *base = motif[i % motif.len()];
        }
    }
    let mut target = query.clone();
    match number % 7 {
        0 => {
            target = (0..if long {
                1 + rng.below(300) as usize
            } else {
                rng.below(17) as usize
            })
                .map(|_| random_base(rng))
                .collect()
        }
        1 if !target.is_empty() => {
            for _ in 0..1 + rng.below(5) {
                let p = rng.below(target.len() as u32) as usize;
                target[p] = random_base(rng);
            }
        }
        2 => {
            let n = 1 + rng.below(30) as usize;
            let at = rng.below((target.len() + 1) as u32) as usize;
            target.splice(at..at, (0..n).map(|_| random_base(rng)));
        }
        3 if !target.is_empty() => {
            let n = (1 + rng.below(30) as usize).min(target.len());
            let at = rng.below((target.len() - n + 1) as u32) as usize;
            target.drain(at..at + n);
        }
        4 if !target.is_empty() => {
            let p = rng.below(target.len() as u32) as usize;
            target[p] = 4;
        }
        5 if !target.is_empty() => {
            let at = rng.below(target.len() as u32) as usize;
            let n = (1 + rng.below(30) as usize).min(target.len() - at);
            let copy = target[at..at + n].to_vec();
            target.splice(at + n..at + n, copy);
        }
        _ => {}
    }
    if target.len() > 300 {
        target.truncate(300);
    }
    let (o_del, e_del, o_ins, e_ins, w, end_bonus, zdrop, h0) = if number.is_multiple_of(3) {
        (
            6,
            1,
            6,
            1,
            [1, 5, 50, 100, 200][number % 5],
            5,
            100,
            [1, 19, 30, 151][number % 4],
        )
    } else {
        (
            1 + rng.below(12) as i32,
            1 + rng.below(4) as i32,
            1 + rng.below(12) as i32,
            1 + rng.below(4) as i32,
            1 + rng.below(200) as i32,
            rng.below(12) as i32,
            rng.below(150) as i32,
            1 + rng.below(200) as i32,
        )
    };
    // ksw_global2 has no no-path guard: when the terminal cell is outside its
    // band its C backtracker indexes outside z. Keep global oracle calls in
    // the defined portion of its interface.
    let w = if extend {
        w
    } else {
        w.max((query.len() as i32 - target.len() as i32).abs())
    };
    Case {
        extend,
        query,
        target,
        mat: bwa_mat(),
        o_del,
        e_del,
        o_ins,
        e_ins,
        w,
        end_bonus,
        zdrop,
        h0,
    }
}

fn i32_bytes(out: &mut Vec<u8>, value: i32) {
    out.extend_from_slice(&value.to_ne_bytes());
}
fn u32_bytes(out: &mut Vec<u8>, value: u32) {
    out.extend_from_slice(&value.to_ne_bytes());
}
fn read_i32(input: &[u8], cursor: &mut usize) -> i32 {
    let b: [u8; 4] = input[*cursor..*cursor + 4].try_into().unwrap();
    *cursor += 4;
    i32::from_ne_bytes(b)
}
fn read_u32(input: &[u8], cursor: &mut usize) -> u32 {
    let b: [u8; 4] = input[*cursor..*cursor + 4].try_into().unwrap();
    *cursor += 4;
    u32::from_ne_bytes(b)
}

fn compile_oracle(dir: &Path) -> PathBuf {
    let manifest = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let c = manifest.join("tests/c");
    let exe = dir.join("ksw_oracle");
    let status = Command::new(env::var("CC").unwrap_or_else(|_| "cc".into()))
        .args(["-std=c99", "-O3", "-I"])
        .arg(&c)
        .arg(c.join("ksw.c"))
        .arg(c.join("harness_ksw/ksw_oracle.c"))
        .arg("-o")
        .arg(&exe)
        .status()
        .expect("run C compiler");
    assert!(status.success(), "C oracle compiler failed: {status}");
    exe
}

fn run_cases(cases: Vec<Case>) {
    let root = env::temp_dir().join(format!(
        "nomics-ksw-{}-{}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ));
    fs::create_dir_all(&root).unwrap();
    let input = root.join("input.bin");
    let output = root.join("output.bin");
    let mut bytes = Vec::new();
    u32_bytes(&mut bytes, cases.len() as u32);
    for case in &cases {
        bytes.push(if case.extend { 0 } else { 1 });
        for n in [
            case.query.len() as i32,
            case.target.len() as i32,
            5,
            case.o_del,
            case.e_del,
            case.o_ins,
            case.e_ins,
            case.w,
            case.end_bonus,
            case.zdrop,
            case.h0,
        ] {
            i32_bytes(&mut bytes, n);
        }
        bytes.extend(case.mat.iter().map(|&x| x as u8));
        bytes.extend_from_slice(&case.query);
        bytes.extend_from_slice(&case.target);
    }
    fs::write(&input, bytes).unwrap();
    let oracle = compile_oracle(&root);
    assert!(
        Command::new(oracle)
            .arg(&input)
            .arg(&output)
            .status()
            .unwrap()
            .success()
    );
    let actual = fs::read(&output).unwrap();
    let mut cursor = 0;
    for (case_no, case) in cases.iter().enumerate() {
        if case.extend {
            let expected = Extend {
                score: read_i32(&actual, &mut cursor),
                qle: read_i32(&actual, &mut cursor),
                tle: read_i32(&actual, &mut cursor),
                gtle: read_i32(&actual, &mut cursor),
                gscore: read_i32(&actual, &mut cursor),
                max_off: read_i32(&actual, &mut cursor),
            };
            let got = extend2(
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
                case.w,
                case.end_bonus,
                case.zdrop,
                case.h0,
            );
            assert_eq!(
                got, expected,
                "extend case {case_no}: q={:?}, t={:?}",
                case.query, case.target
            );
        } else {
            let expected_score = read_i32(&actual, &mut cursor);
            let n = read_u32(&actual, &mut cursor) as usize;
            let expected_cigar = (0..n)
                .map(|_| read_u32(&actual, &mut cursor))
                .collect::<Vec<_>>();
            let got = global2(
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
                case.w,
            );
            assert_eq!(
                got,
                (expected_score, expected_cigar),
                "global case {case_no}: q={:?}, t={:?}",
                case.query,
                case.target
            );
        }
    }
    assert_eq!(cursor, actual.len());
    fs::remove_dir_all(root).unwrap();
}

#[test]
fn ksw_smoke_vs_bwa_c() {
    let mut rng = Rng(0x8d1a_7a13_5eed_0123);
    let mut cases = Vec::new();
    for n in 0..1_000 {
        let mut c = make_case(&mut rng, true, n);
        if c.query.is_empty() {
            c.query.push(0);
        }
        cases.push(c);
    }
    for n in 0..500 {
        let mut c = make_case(&mut rng, false, n);
        if c.query.is_empty() || c.target.is_empty() {
            c.query = vec![0];
            c.target = vec![0];
        }
        cases.push(c);
    }
    // ksw_global2's empty sides are defined by the C implementation (unlike
    // extend2's empty query allocation bug) and belong in the oracle run.
    cases.push(Case {
        extend: false,
        query: Vec::new(),
        target: vec![0, 1, 2],
        mat: bwa_mat(),
        o_del: 6,
        e_del: 1,
        o_ins: 6,
        e_ins: 1,
        w: 100,
        end_bonus: 5,
        zdrop: 100,
        h0: 19,
    });
    cases.push(Case {
        extend: false,
        query: vec![0, 1, 2],
        target: Vec::new(),
        mat: bwa_mat(),
        o_del: 6,
        e_del: 1,
        o_ins: 6,
        e_ins: 1,
        w: 100,
        end_bonus: 5,
        zdrop: 100,
        h0: 19,
    });
    run_cases(cases);
}

#[test]
#[ignore = "200k/50k C differential campaign; run with --release -- --ignored"]
fn ksw_large_vs_bwa_c() {
    let mut rng = Rng(0x8d1a_7a13_5eed_0123);
    let mut cases = Vec::with_capacity(250_000);
    for n in 0..200_000 {
        let mut c = make_case(&mut rng, true, n);
        if c.query.is_empty() {
            c.query.push(0);
        }
        cases.push(c);
    }
    for n in 0..50_000 {
        let mut c = make_case(&mut rng, false, n);
        if c.query.is_empty() || c.target.is_empty() {
            c.query = vec![0];
            c.target = vec![0];
        }
        cases.push(c);
    }
    run_cases(cases);
}

#[test]
fn empty_inputs_have_defined_rust_results() {
    let mat = bwa_mat();
    assert_eq!(
        extend2(0, &[], 0, &[], 5, &mat, 6, 1, 6, 1, 100, 5, 100, 19),
        Extend {
            score: 19,
            qle: 0,
            tle: 0,
            gtle: 0,
            gscore: -1,
            max_off: 0
        }
    );
    assert_eq!(
        global2(0, &[], 3, &[0, 1, 2], 5, &mat, 6, 1, 6, 1, 100),
        (-9, vec![3 << 4 | 2])
    );
    assert_eq!(
        global2(3, &[0, 1, 2], 0, &[], 5, &mat, 6, 1, 6, 1, 100),
        (-9, vec![3 << 4 | 1])
    );
}

#[test]
#[should_panic(expected = "h0 > 0")]
fn extend_h0_zero_matches_bwa_assertion_contract() {
    let mat = bwa_mat();
    let _ = extend2(1, &[0], 1, &[0], 5, &mat, 6, 1, 6, 1, 100, 5, 100, 0);
}
