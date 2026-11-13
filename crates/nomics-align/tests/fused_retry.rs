//! The fused-extension CPU stand-in against a scalar transcription of
//! `chain2aln`'s two band-retry loops, on synthetic seeds that really retry
//! (75-100 bp indels inside a 256 bp query), plus edge cases. The jobs and
//! the stand-in's results are committed as `fixtures/fused-retry.bin` (the
//! `NOMICS_DUMP_FUSED` format) for `fused_bench` and the CUDA kernel gate;
//! `fixture_is_current` regenerates it and byte-compares, so it cannot drift.
//! `NOMICS_BLESS_FUSED=1` rewrites it.

use nomics_align::extend::{
    E_DEL, E_INS, ExtKey, MAT, O_DEL, O_INS, PEN_CLIP3, PEN_CLIP5, W, ZDROP, fused_sides,
};
use nomics_align::fused::{self, dump};
use nomics_align::ksw::{self, Extend};
use nomics_gpu::{GpuFusedInput, GpuFusedResult};

const FIXTURE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/tests/fixtures/fused-retry.bin"
);

struct Rng(u64);

impl Rng {
    fn below(&mut self, n: u32) -> u32 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 as u32) % n
    }

    fn bases(&mut self, n: usize) -> Vec<u8> {
        (0..n).map(|_| self.below(4) as u8).collect()
    }
}

/// One synthetic job; `expect` is `(n_left, n_right)` where the case was
/// built to pin it.
struct Case {
    name: &'static str,
    lq: Vec<u8>,
    lt: Vec<u8>,
    rq: Vec<u8>,
    rt: Vec<u8>,
    h0: i32,
    expect: Option<(i32, i32)>,
}

impl Case {
    fn input(&self) -> GpuFusedInput<'_> {
        GpuFusedInput {
            lq: &self.lq,
            lt: &self.lt,
            rq: &self.rq,
            rt: &self.rt,
            h0: self.h0,
        }
    }
}

/// Query and target that agree on `segs[0]` bases, then for each later
/// segment an indel of `gaps[k]` (positive: deletion, extra target bases;
/// negative: insertion, extra query bases) followed by `segs[k + 1]` matching
/// bases, then `tail` random target bases.
fn indels(rng: &mut Rng, segs: &[usize], gaps: &[i32], tail: usize) -> (Vec<u8>, Vec<u8>) {
    assert_eq!(segs.len(), gaps.len() + 1);
    let (mut q, mut t) = (Vec::new(), Vec::new());
    for (k, &n) in segs.iter().enumerate() {
        if k > 0 {
            let gap = rng.bases(gaps[k - 1].unsigned_abs() as usize);
            if gaps[k - 1] > 0 {
                t.extend(&gap);
            } else {
                q.extend(&gap);
            }
        }
        let same = rng.bases(n);
        q.extend(&same);
        t.extend(&same);
    }
    t.extend(rng.bases(tail));
    (q, t)
}

/// A read-like pair: `qlen` query bases, the target a copy with a few
/// substitutions (some N) and small indels, plus a random tail.
fn readlike(rng: &mut Rng, qlen: usize) -> (Vec<u8>, Vec<u8>) {
    let q: Vec<u8> = (0..qlen)
        .map(|_| {
            if rng.below(50) == 0 {
                4
            } else {
                rng.below(4) as u8
            }
        })
        .collect();
    let mut t = q.clone();
    for _ in 0..rng.below(5) {
        if t.is_empty() {
            break;
        }
        let p = rng.below(t.len() as u32) as usize;
        match rng.below(3) {
            0 => t[p] = rng.below(5) as u8,
            1 => t.insert(p, rng.below(4) as u8),
            _ => {
                t.remove(p);
            }
        }
    }
    let tail = rng.below(120) as usize;
    t.extend(rng.bases(tail));
    (q, t)
}

fn cases() -> Vec<Case> {
    let mut rng = Rng(0x0f05_ed00_5eed_0001);
    let mut v = Vec::new();
    let none = || (Vec::new(), Vec::new());
    let mut push = |name, l: (Vec<u8>, Vec<u8>), r: (Vec<u8>, Vec<u8>), h0, expect| {
        // Only jobs the planner could emit: the GPU packer rejects the rest.
        assert!(l.0.len() <= 256 && r.0.len() <= 256, "{name}: query over 256");
        v.push(Case {
            name,
            lq: l.0,
            lt: l.1,
            rq: r.0,
            rt: r.1,
            h0,
            expect,
        })
    };
    // Band 0 (w = 100) ends 80-95 off the diagonal (max_off >= 75) with a
    // new score: both sides call band 1. h0 plus the first segment outscores
    // the gap, so the local extension does not die inside it.
    let l = indels(&mut rng, &[40, 120], &[80], 30);
    push("left retry (80 bp deletion)", l, none(), 60, Some((2, 0)));
    let r = indels(&mut rng, &[30, 130], &[-80], 40);
    push("right retry (80 bp insertion)", none(), r, 70, Some((0, 2)));
    let l = indels(&mut rng, &[10, 140], &[-90], 50);
    let r = indels(&mut rng, &[25, 120], &[85], 10);
    push("both retry", l, r, 100, Some((2, 2)));
    let l = indels(&mut rng, &[5, 150], &[95], 0);
    let r = indels(&mut rng, &[5, 150], &[-95], 0);
    push(
        "both retry (95 bp, no target tail)",
        l,
        r,
        100,
        Some((2, 2)),
    );
    // Two deletions (80 then 40): band 0 is stuck 80 off, band 1 reaches the
    // third segment and scores higher, so right's h0 is band 1's score.
    let l = indels(&mut rng, &[20, 100, 60], &[80, 40], 20);
    let r = readlike(&mut rng, 60);
    push("left retry, band 1 rescores", l, r, 80, Some((2, 1)));
    let l = readlike(&mut rng, 70);
    let r = indels(&mut rng, &[20, 100, 60], &[80, 40], 20);
    push("right retry, band 1 rescores", l, r, 21, Some((1, 2)));
    // No retries.
    let l = readlike(&mut rng, 40);
    push("left only", l, none(), 60, Some((1, 0)));
    let r = readlike(&mut rng, 90);
    push("right only", none(), r, 61, Some((0, 1)));
    let (l, r) = (readlike(&mut rng, 30), readlike(&mut rng, 80));
    push("both", l, r, 41, Some((1, 1)));
    // A right query over 256: the planner sends a left-only job.
    let l = readlike(&mut rng, 12);
    push(
        "right query > 256 (left-only job)",
        l,
        none(),
        19,
        Some((1, 0)),
    );
    // Empty targets.
    let l = (rng.bases(30), Vec::new());
    let r = readlike(&mut rng, 30);
    push("left tlen 0", l, r, 50, Some((1, 1)));
    let l = readlike(&mut rng, 30);
    let r = (rng.bases(30), Vec::new());
    push("right tlen 0", l, r, 50, Some((1, 1)));
    // Query lengths 1 and 256.
    let (l, r) = (readlike(&mut rng, 1), readlike(&mut rng, 1));
    push("qlen 1", l, r, 148, Some((1, 1)));
    let (l, r) = (readlike(&mut rng, 256), readlike(&mut rng, 256));
    push("qlen 256", l, r, 19, None);
    let l = indels(&mut rng, &[40, 136], &[-80], 64);
    let r = indels(&mut rng, &[40, 136], &[80], 64);
    push("qlen 256 with retries", l, r, 60, Some((2, 2)));
    // Read-like jobs of every shape.
    for k in 0..240 {
        let side = |rng: &mut Rng, on: bool| {
            if on {
                let n = 1 + rng.below(256) as usize;
                readlike(rng, n)
            } else {
                (Vec::new(), Vec::new())
            }
        };
        let (l, r) = match k % 4 {
            0 => (side(&mut rng, true), side(&mut rng, false)),
            1 => (side(&mut rng, false), side(&mut rng, true)),
            _ => (side(&mut rng, true), side(&mut rng, true)),
        };
        let h0 = 19 + rng.below(130) as i32;
        push("read-like", l, r, h0, None);
    }
    // Kernel v3 bins jobs by max(lq, rq) onto <W, C> instantiations: every
    // multiple of 8 and the length after it, the longer side either left or
    // right, the other side short (1 base) or random.
    for q in (8..=256).step_by(8).flat_map(|q| [q, q + 1]).filter(|&q| q <= 256) {
        for k in 0..2 {
            let long = readlike(&mut rng, q);
            let n = if k == 0 { 1 } else { 1 + rng.below(q as u32) as usize };
            let short = readlike(&mut rng, n);
            let h0 = 19 + rng.below(130) as i32;
            if (q + k) % 2 == 0 {
                push("bin edge, long left", long, short, h0, Some((1, 1)));
            } else {
                push("bin edge, long right", short, long, h0, Some((1, 1)));
            }
        }
    }
    // Large h0 on each side, around 2^22 - qlen (the bound a packed
    // score/column reduction tried for kernel v3 needed; kept as large-score
    // coverage), and a left whose sc0 carries the right past it.
    const BIG: i32 = 1 << 22;
    for (q, h0) in [(150, BIG - 151), (150, BIG - 150), (40, BIG - 41), (40, BIG)] {
        let l = readlike(&mut rng, q);
        push("large h0, left", l, none(), h0, Some((1, 0)));
        let r = readlike(&mut rng, q);
        push("large h0, right", none(), r, h0, Some((0, 1)));
    }
    let (l, r) = (readlike(&mut rng, 100), readlike(&mut rng, 150));
    push("large h0, left then right", l, r, BIG - 200, Some((1, 1)));
    v
}

/// One `ksw_extend2` call as `chain2aln` makes it: band, h0 and result.
type Call = (i32, i32, Extend);

/// `chain2aln`'s two loops for one seed, transcribed as BWA writes them
/// (`a.score` starts at -1; `prev` is sampled at the top of each try).
fn chain2aln_loops(c: &GpuFusedInput) -> (Vec<Call>, Vec<Call>) {
    let call = |q: &[u8], t: &[u8], w, end_bonus, h0| {
        ksw::extend2(
            q.len() as i32,
            q,
            t.len() as i32,
            t,
            5,
            &MAT,
            O_DEL,
            E_DEL,
            O_INS,
            E_INS,
            w,
            end_bonus,
            ZDROP,
            h0,
        )
    };
    let (mut left, mut right) = (Vec::new(), Vec::new());
    let mut a_score = -1;
    if !c.lq.is_empty() {
        for i in 0..2 {
            let prev = a_score;
            let aw = W << i;
            let x = call(c.lq, c.lt, aw, PEN_CLIP5, c.h0);
            left.push((aw, c.h0, x));
            a_score = x.score;
            if a_score == prev || x.max_off < (aw >> 1) + (aw >> 2) {
                break;
            }
        }
    } else {
        a_score = c.h0;
    }
    if !c.rq.is_empty() {
        let sc0 = a_score;
        for i in 0..2 {
            let prev = a_score;
            let aw = W << i;
            let x = call(c.rq, c.rt, aw, PEN_CLIP3, sc0);
            right.push((aw, sc0, x));
            a_score = x.score;
            if a_score == prev || x.max_off < (aw >> 1) + (aw >> 2) {
                break;
            }
        }
    }
    (left, right)
}

fn check(c: &Case, r: &GpuFusedResult) {
    let (left, right) = chain2aln_loops(&c.input());
    let name = c.name;
    let mut keys = [Vec::new(), Vec::new()];
    let mut results = [Vec::new(), Vec::new()];
    fused::push_keys(&fused::params(), 7, 3, c.h0, r, &mut keys, &mut results);
    for (side, calls) in [left, right].iter().enumerate() {
        assert_eq!(keys[side].len(), calls.len(), "{name}: side {side} calls");
        for (i, &(w, h0, x)) in calls.iter().enumerate() {
            let key = ExtKey {
                chain: 7,
                seed: 3,
                right: side == 1,
                w,
                h0,
            };
            assert_eq!(keys[side][i], key, "{name}: side {side} call {i} key");
            assert_eq!(results[side][i], x, "{name}: side {side} call {i} result");
        }
    }
    if let Some((n_left, n_right)) = c.expect {
        assert_eq!((r.n_left, r.n_right), (n_left, n_right), "{name}");
    }
}

#[test]
fn stand_in_matches_chain2aln_loops() {
    let p = fused::params();
    let cases = cases();
    let mut tries = [[0; 3]; 2];
    for c in &cases {
        let r = fused::extend_cpu(&p, &c.input());
        check(c, &r);
        tries[0][r.n_left as usize] += 1;
        tries[1][r.n_right as usize] += 1;
    }
    // The rescoring cases: band 1 changes the left score, and right's h0
    // follows it.
    let c = cases
        .iter()
        .find(|c| c.name == "left retry, band 1 rescores")
        .unwrap();
    let r = fused::extend_cpu(&p, &c.input());
    assert!(r.left[1].score > r.left[0].score, "{:?}", r);
    assert_eq!(r.right_h0, r.left[1].score);
    let c = cases
        .iter()
        .find(|c| c.name == "right retry, band 1 rescores")
        .unwrap();
    let r = fused::extend_cpu(&p, &c.input());
    assert!(r.right[1].score > r.right[0].score, "{:?}", r);
    assert!(tries[0][2] >= 5 && tries[1][2] >= 5, "{tries:?}");
}

#[test]
fn engine_matches_per_job_stand_in() {
    // The engine's parallel stand-in, over uneven groups, in order.
    let p = fused::params();
    let cases = cases();
    let bounds = [0, 3, 3, 40, 41, cases.len()];
    let groups: Vec<_> = bounds
        .windows(2)
        .map(|b| cases[b[0]..b[1]].iter().map(Case::input))
        .collect();
    let mut engine = fused::FusedEngine::Cpu(Vec::new());
    for threads in [1, 3, 8] {
        let (out, _) = engine.compute(&p, &groups, threads).unwrap();
        assert_eq!(out.len(), cases.len());
        for (c, r) in cases.iter().zip(out) {
            assert_eq!(*r, fused::extend_cpu(&p, &c.input()), "{}", c.name);
        }
    }
}

#[test]
fn sides_follow_the_query_limit() {
    // (qbeg, len, l_query) -> (left, right)
    assert_eq!(fused_sides(10, 20, 150), Some((true, true)));
    assert_eq!(fused_sides(0, 20, 150), Some((false, true)));
    assert_eq!(fused_sides(130, 20, 150), Some((true, false)));
    assert_eq!(fused_sides(0, 150, 150), None);
    assert_eq!(fused_sides(256, 20, 300), Some((true, true)));
    assert_eq!(fused_sides(257, 20, 300), None);
    assert_eq!(fused_sides(10, 20, 286), Some((true, true)));
    assert_eq!(fused_sides(10, 20, 287), Some((true, false)));
    assert_eq!(fused_sides(0, 20, 287), None);
}

fn fixture_bytes() -> Vec<u8> {
    let p = fused::params();
    let mut b = Vec::new();
    dump::encode_header(&mut b, &p);
    for c in cases() {
        dump::encode_job(&mut b, &c.input(), &fused::extend_cpu(&p, &c.input()));
    }
    b
}

#[test]
fn fixture_is_current() {
    let bytes = fixture_bytes();
    if std::env::var_os("NOMICS_BLESS_FUSED").is_some() {
        std::fs::write(FIXTURE, &bytes).unwrap();
    }
    assert!(bytes.len() < 1 << 20);
    let committed =
        std::fs::read(FIXTURE).expect("fused-retry.bin (NOMICS_BLESS_FUSED=1 writes it)");
    assert!(
        committed == bytes,
        "fused-retry.bin is stale; NOMICS_BLESS_FUSED=1 rewrites it"
    );
    // And it parses back to the same jobs and results.
    let d = dump::Dump::parse(committed);
    assert_eq!(d.params, fused::params());
    let cases = cases();
    assert_eq!(d.jobs.len(), cases.len());
    for (j, c) in d.jobs.iter().zip(&cases) {
        let x = d.input(j);
        assert_eq!(
            (x.lq, x.lt, x.rq, x.rt, x.h0),
            (&c.lq[..], &c.lt[..], &c.rq[..], &c.rt[..], c.h0)
        );
        let mut want = fused::extend_cpu(&d.params, &x);
        // Unused slots are not stored.
        for i in want.n_left as usize..2 {
            want.left[i] = Default::default();
        }
        for i in want.n_right as usize..2 {
            want.right[i] = Default::default();
        }
        assert_eq!(j.result, want, "{}", c.name);
    }
}
