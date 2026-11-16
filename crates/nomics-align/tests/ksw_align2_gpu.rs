//! The GPU `ksw_align2` (mate rescue) against the CPU port on synthetic
//! rescue-shaped calls, including windows up to `GPU_RESCUE_MAX_TARGET`
//! (4,000 bases). Those long targets shrink how many alignments share a
//! threadgroup on Metal (32 KiB of threadgroup memory), a path the recorded
//! dumps (`ksw_bench`, targets up to ~1.3 kb) never reach. Read lengths up to
//! 500 bp cover every kernel width (4, 5, 6, 8 and 16 columns per lane); the
//! dumps' 151 bp reads use only the first two.
#![cfg(any(feature = "cuda", feature = "metal"))]

use nomics_align::ksw::{KSW_XBYTE, KSW_XSTART, KSW_XSUBO, Kswr, align2};
use nomics_gpu::{GpuKswAlign2, GpuKswAlign2Input, GpuKswResult};

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u32 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 >> 32) as u32
    }

    fn below(&mut self, n: u32) -> u32 {
        self.next() % n
    }
}

struct Case {
    query: Vec<u8>,
    target: Vec<u8>,
    xtra: i32,
}

/// BWA's default scoring: match 1, mismatch -4, N -1.
fn matrix() -> Vec<i8> {
    let mut mat = vec![0i8; 25];
    for i in 0..5 {
        for j in 0..5 {
            mat[i * 5 + j] = if i == 4 || j == 4 {
                -1
            } else if i == j {
                1
            } else {
                -4
            };
        }
    }
    mat
}

/// A random window holding up to two mutated copies of the read (two exercise
/// the second-best score), as a rescue window around a mate does.
fn case(rng: &mut Rng, qlens: std::ops::Range<usize>, tlen: usize) -> Case {
    let qlen = qlens.start + rng.below((qlens.end - qlens.start) as u32) as usize;
    // Half the cases are short tandem repeats: many cells tie for a row or
    // column maximum there, so endpoint tie-breaks decide the result.
    let period = if rng.below(2) == 0 {
        0
    } else {
        1 + rng.below(4) as usize
    };
    let motif: Vec<u8> = (0..period.max(1)).map(|_| rng.below(4) as u8).collect();
    let mut base = |rng: &mut Rng, k: usize| {
        if period == 0 || rng.below(25) == 0 {
            rng.below(4) as u8
        } else {
            motif[k % period]
        }
    };
    let query: Vec<u8> = (0..qlen).map(|k| base(rng, k)).collect();
    let mut target: Vec<u8> = (0..tlen).map(|k| base(rng, k)).collect();
    // No copy (a failed rescue: weak hits, many ties), one, or two.
    let copies = rng.below(3);
    for _ in 0..copies {
        if tlen <= qlen {
            break;
        }
        let at = rng.below((tlen - qlen) as u32) as usize;
        for (k, &base) in query.iter().enumerate() {
            target[at + k] = match rng.below(40) {
                0 => rng.below(4) as u8,
                1 => 4,
                _ => base,
            };
        }
    }
    // mem_matesw: XSUBO | XSTART | (l_ms * a < 250 ? XBYTE : 0) | min_seed_len * a
    // (a = 1). Short reads also take the 16-bit path at random: always valid,
    // and it covers that path at every width.
    let byte = if qlen < 250 && rng.below(2) == 0 {
        KSW_XBYTE
    } else {
        0
    };
    Case {
        query,
        target,
        xtra: KSW_XSUBO | KSW_XSTART | byte | 19,
    }
}

fn cpu(case: &Case, mat: &[i8], g: [i32; 4]) -> Kswr {
    align2(
        case.query.len() as i32,
        &case.query,
        case.target.len() as i32,
        &case.target,
        5,
        mat,
        g[0],
        g[1],
        g[2],
        g[3],
        case.xtra,
        None,
    )
}

fn gpu(cases: &[Case], mat: &[i8], g: [i32; 4]) -> Vec<GpuKswResult> {
    let inputs: Vec<GpuKswAlign2Input<'_>> = cases
        .iter()
        .map(|c| GpuKswAlign2Input {
            query: &c.query,
            target: &c.target,
            m: 5,
            mat,
            o_del: g[0],
            e_del: g[1],
            o_ins: g[2],
            e_ins: g[3],
            xtra: c.xtra,
            size: if c.xtra & KSW_XBYTE != 0 { 1 } else { 2 },
        })
        .collect();
    let mut out = vec![GpuKswResult::default(); inputs.len()];
    GpuKswAlign2::new()
        .compute(&inputs, &mut out)
        .expect("GPU ksw_align2");
    out
}

fn check(cases: &[Case]) {
    check_with(cases, &matrix(), [6, 1, 6, 1]);
}

fn check_with(cases: &[Case], mat: &[i8], gaps: [i32; 4]) {
    let got = gpu(cases, mat, gaps);
    for (i, (case, g)) in cases.iter().zip(&got).enumerate() {
        let want = cpu(case, mat, gaps);
        let g = Kswr {
            score: g.score,
            te: g.te,
            qe: g.qe,
            score2: g.score2,
            te2: g.te2,
            tb: g.tb,
            qb: g.qb,
        };
        assert_eq!(
            (g.score, g.te, g.qe, g.score2, g.te2, g.tb, g.qb),
            (
                want.score,
                want.te,
                want.qe,
                want.score2,
                want.te2,
                want.tb,
                want.qb
            ),
            "call {i}: qlen {} tlen {} xtra {:#x}",
            case.query.len(),
            case.target.len(),
            case.xtra
        );
    }
}

#[test]
fn long_rescue_windows_match_cpu() {
    let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
    // Every target near the 4,000-base cap: the fewest alignments per group.
    let long: Vec<Case> = (0..96)
        .map(|_| {
            let tlen = 3500 + rng.below(501) as usize;
            case(&mut rng, 100..152, tlen)
        })
        .collect();
    check(&long);
    // Mixed lengths in one launch: the longest sets the group shape.
    let mixed: Vec<Case> = (0..400)
        .map(|_| {
            let tlen = [200, 800, 1500, 2600, 4000][rng.below(5) as usize];
            case(&mut rng, 100..152, tlen)
        })
        .collect();
    check(&mixed);
}

/// Bounded sub-launches (several per batch, each with its own width class
/// and group shape) give the one-launch results in input order.
#[test]
fn bounded_launches_match_one_launch() {
    let mut rng = Rng(0x51_7cc1_b727_220a);
    let cases: Vec<Case> = (0..600)
        .map(|_| {
            let tlen = [150, 700, 1500, 2600, 4000][rng.below(5) as usize];
            case(&mut rng, 100..300, tlen)
        })
        .collect();
    let mat = matrix();
    let inputs: Vec<GpuKswAlign2Input<'_>> = cases
        .iter()
        .map(|c| GpuKswAlign2Input {
            query: &c.query,
            target: &c.target,
            m: 5,
            mat: &mat,
            o_del: 6,
            e_del: 1,
            o_ins: 6,
            e_ins: 1,
            xtra: c.xtra,
            size: if c.xtra & KSW_XBYTE != 0 { 1 } else { 2 },
        })
        .collect();
    let fields = |r: &GpuKswResult| (r.score, r.te, r.qe, r.score2, r.te2, r.tb, r.qb);
    let mut engine = GpuKswAlign2::new();
    let mut whole = vec![GpuKswResult::default(); inputs.len()];
    engine
        .compute_bounded(&inputs, &mut whole, usize::MAX)
        .expect("one launch");
    // 1 byte: one call per launch; then a few calls per launch, reusing the
    // same engine (its buffers are larger than each launch needs).
    for max_bytes in [1, 5_000, 40_000] {
        let mut sliced = vec![GpuKswResult::default(); inputs.len()];
        engine
            .compute_bounded(&inputs, &mut sliced, max_bytes)
            .expect("bounded launches");
        for (i, (a, b)) in whole.iter().zip(&sliced).enumerate() {
            assert_eq!(fields(a), fields(b), "call {i}, max_bytes {max_bytes}");
        }
    }
    check(&cases);
}

#[test]
fn every_kernel_width_matches_cpu() {
    let mut rng = Rng(0x2545_f491_4f6c_dd1d);
    // The widest padded query of a launch picks the kernel: one launch per
    // width class, so each width runs (padding is to 16 or 8 lanes).
    for qlens in [100..129, 140..161, 170..193, 230..257, 300..513] {
        let cases: Vec<Case> = (0..160)
            .map(|_| {
                let tlen = [qlens.end, 900, 2000, 4000][rng.below(4) as usize];
                case(&mut rng, qlens.clone(), tlen)
            })
            .collect();
        check(&cases);
    }
}

/// CUDA runs byte-sized calls with P = 160 (145-151 bp reads at C = 5) two
/// per warp. Pairs of equal and of different target lengths, most targets
/// free of N (a pair needs bases 0-3) and some not, an odd count (the last
/// call runs alone), and a match score of 2, with which reads whose copy is
/// in the window reach the byte ceiling (score 255; without KSW_XSTART: a
/// saturated result has qe = -1, an empty query for the reverse pass).
#[test]
fn paired_calls_match_cpu() {
    let mut rng = Rng(0x7f4a_7c15_9e37_79b9);
    let mut mat2 = matrix();
    for i in 0..4 {
        mat2[i * 5 + i] = 2;
    }
    for (n, mat, start) in [(1201, matrix(), KSW_XSTART), (301, mat2, 0)] {
        let cases: Vec<Case> = (0..n)
            .map(|i| {
                let tlen = if i % 3 == 0 {
                    1334
                } else {
                    [160, 700, 1334, 2600][rng.below(4) as usize] - rng.below(30) as usize
                };
                let mut c = case(&mut rng, 145..152, tlen);
                c.xtra = (c.xtra & !KSW_XSTART) | KSW_XBYTE | start;
                if rng.below(8) != 0 {
                    for b in &mut c.target {
                        if *b == 4 {
                            *b = rng.below(4) as u8;
                        }
                    }
                }
                c
            })
            .collect();
        check_with(&cases, &mat, [6, 1, 6, 1]);
        // Unequal deletion and insertion gaps (the kernels share one
        // gap term when they are equal).
        check_with(&cases, &mat, [5, 2, 7, 1]);
    }
}

/// Scoring the GPU kernels (CUDA and Metal) do not reproduce is refused at launch, so the
/// caller computes those calls on the CPU: gaps that BWA's byte path
/// truncates to uint8_t, and 16-bit calls whose scores could pass BWA's
/// saturation at 32767. The CPU results are BWA's (from the review's
/// counterexamples).
#[test]
fn out_of_range_scoring_goes_to_the_cpu() {
    let run = |q: Vec<u8>, t: Vec<u8>, mat: &[i8], g: [i32; 4], xtra: i32| {
        let case = Case {
            query: q,
            target: t,
            xtra,
        };
        let want = cpu(&case, mat, g);
        let input = GpuKswAlign2Input {
            query: &case.query,
            target: &case.target,
            m: 5,
            mat,
            o_del: g[0],
            e_del: g[1],
            o_ins: g[2],
            e_ins: g[3],
            xtra,
            size: if xtra & KSW_XBYTE != 0 { 1 } else { 2 },
        };
        let mut out = [GpuKswResult::default()];
        let got = GpuKswAlign2::new().compute(&[input], &mut out);
        (want, got.is_err(), out[0])
    };
    let rep = |b: u8, n: usize| std::iter::repeat_n(b, n);
    // Byte path, o_del + e_del = 256: BWA's uint8_t gap terms wrap.
    let q: Vec<u8> = rep(0, 60).chain(rep(1, 60)).collect();
    let t: Vec<u8> = rep(0, 60).chain(rep(2, 20)).chain(rep(1, 60)).collect();
    let (want, refused, got) = run(q, t, &matrix(), [255, 1, 6, 1], 0xd0013);
    let fields = (want.score, want.te, want.qe, want.score2, want.te2, want.tb, want.qb);
    assert_eq!(fields, (120, 139, 119, -1, -1, 0, 0));
    assert!(refused, "the byte call is refused, not {got:?}");
    // 16-bit path, 259 matches of score 127: BWA saturates at 32767.
    let mut mat = matrix();
    for i in 0..4 {
        mat[i * 5 + i] = 127;
    }
    let (want, refused, got) = run(rep(0, 259).collect(), rep(0, 259).collect(), &mat, [6, 1, 6, 1], 0xc0013);
    assert_eq!(want.score, 32767);
    assert!(refused, "the 16-bit call is refused, not {got:?}");
}
