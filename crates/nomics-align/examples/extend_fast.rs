//! Compare `ksw::extend2` (with exact fast paths) to the literal
//! `ksw::extend2_reference` on a dump of real calls.
//!
//! Usage: `extend_fast <extend2.bin> [threads] [reps]`. Verifies every call of
//! both functions against the recorded BWA result, prints the share of calls
//! and of reference time that are perfect-match extensions, then times both
//! functions at 1 thread (one pass) and at `threads` (eight passes per
//! timing), `reps` alternating runs each.

use std::hint::black_box;
use std::ops::Range;
use std::path::Path;
use std::time::Instant;

use nomics_align::ksw::{self, Extend};
use rayon::prelude::*;

struct Call {
    query: Range<usize>,
    target: Range<usize>,
    m: i32,
    mat: usize,
    args: [i32; 8],
    expected: Extend,
}

struct Dump {
    bytes: Vec<u8>,
    mats: Vec<Vec<i8>>,
    calls: Vec<Call>,
}

fn take(bytes: &[u8], at: &mut usize, len: usize) -> Range<usize> {
    let r = *at..*at + len;
    assert!(r.end <= bytes.len(), "truncated dump");
    *at = r.end;
    r
}

fn int(bytes: &[u8], at: &mut usize) -> i32 {
    let r = take(bytes, at, 4);
    i32::from_le_bytes(bytes[r].try_into().unwrap())
}

impl Dump {
    fn read(path: &Path) -> Self {
        let bytes = std::fs::read(path).expect("read dump");
        let (mut at, mut mats, mut calls) = (0usize, Vec::<Vec<i8>>::new(), Vec::new());
        while at < bytes.len() {
            let ql = int(&bytes, &mut at) as usize;
            let query = take(&bytes, &mut at, ql);
            let tl = int(&bytes, &mut at) as usize;
            let target = take(&bytes, &mut at, tl);
            let m = int(&bytes, &mut at);
            let r = take(&bytes, &mut at, (m * m) as usize);
            let v: Vec<i8> = bytes[r].iter().map(|&x| x as i8).collect();
            let mat = mats.iter().position(|x| *x == v).unwrap_or_else(|| {
                mats.push(v);
                mats.len() - 1
            });
            let args = std::array::from_fn(|_| int(&bytes, &mut at));
            let r: [i32; 6] = std::array::from_fn(|_| int(&bytes, &mut at));
            calls.push(Call {
                query,
                target,
                m,
                mat,
                args,
                expected: Extend {
                    score: r[0],
                    qle: r[1],
                    tle: r[2],
                    gtle: r[3],
                    gscore: r[4],
                    max_off: r[5],
                },
            });
        }
        Self { bytes, mats, calls }
    }

    fn run(&self, c: &Call, reference: bool) -> Extend {
        let q = &self.bytes[c.query.clone()];
        let t = &self.bytes[c.target.clone()];
        let a = c.args;
        let f = if reference {
            ksw::extend2_reference
        } else {
            ksw::extend2
        };
        f(
            q.len() as i32,
            q,
            t.len() as i32,
            t,
            c.m,
            &self.mats[c.mat],
            a[0],
            a[1],
            a[2],
            a[3],
            a[4],
            a[5],
            a[6],
            a[7],
        )
    }

    /// Query equals the start of the target, and every query base scores the
    /// matrix maximum against itself.
    fn perfect(&self, c: &Call) -> bool {
        let q = &self.bytes[c.query.clone()];
        let t = &self.bytes[c.target.clone()];
        let mat = &self.mats[c.mat];
        let m = c.m as usize;
        let best = *mat.iter().max().unwrap();
        !q.is_empty()
            && t.len() >= q.len()
            && q == &t[..q.len()]
            && q.iter().all(|&b| mat[b as usize * m + b as usize] == best)
    }
}

fn time_serial(d: &Dump, idx: &[usize], reference: bool) -> f64 {
    let s = Instant::now();
    for &i in idx {
        black_box(d.run(black_box(&d.calls[i]), reference));
    }
    s.elapsed().as_secs_f64()
}

const PAR_PASSES: usize = 8;

fn time_par(pool: &rayon::ThreadPool, d: &Dump, reference: bool) -> f64 {
    let s = Instant::now();
    // Several passes so that one timing is long enough to be stable.
    for _ in 0..PAR_PASSES {
        pool.install(|| {
            d.calls.par_iter().with_min_len(256).for_each(|c| {
                black_box(d.run(black_box(c), reference));
            })
        });
    }
    s.elapsed().as_secs_f64()
}

fn stats(v: &[f64]) -> (f64, f64, f64) {
    let mean = v.iter().sum::<f64>() / v.len() as f64;
    let min = v.iter().cloned().fold(f64::INFINITY, f64::min);
    let max = v.iter().cloned().fold(0.0, f64::max);
    (mean, min, max)
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let d = Dump::read(Path::new(&args[1]));
    let threads: usize = args.get(2).map_or_else(
        || std::thread::available_parallelism().map_or(1, usize::from),
        |s| s.parse().unwrap(),
    );
    let reps: usize = args.get(3).map_or(5, |s| s.parse().unwrap());
    let n = d.calls.len();

    // Exactness of both entry points against the recorded BWA output.
    for reference in [true, false] {
        let bad = d
            .calls
            .par_iter()
            .filter(|c| d.run(c, reference) != c.expected)
            .count();
        println!(
            "{}: {bad} mismatches / {n} calls",
            if reference { "reference" } else { "extend2  " }
        );
        assert_eq!(bad, 0);
    }

    // Workload shape.
    let all: Vec<usize> = (0..n).collect();
    let perfect: Vec<usize> = all
        .iter()
        .copied()
        .filter(|&i| d.perfect(&d.calls[i]))
        .collect();
    let other: Vec<usize> = all
        .iter()
        .copied()
        .filter(|&i| !d.perfect(&d.calls[i]))
        .collect();
    let _ = time_serial(&d, &all[..n.min(20000)], true);
    let t_perf = time_serial(&d, &perfect, true);
    let t_other = time_serial(&d, &other, true);
    let qsum: usize = d.calls.iter().map(|c| c.query.len()).sum();
    let tsum: usize = d.calls.iter().map(|c| c.target.len()).sum();
    println!(
        "calls {n}, mean qlen {:.1}, mean tlen {:.1}; perfect-match calls {} ({:.1}%), reference time {:.3}s of {:.3}s ({:.1}%)",
        qsum as f64 / n as f64,
        tsum as f64 / n as f64,
        perfect.len(),
        100.0 * perfect.len() as f64 / n as f64,
        t_perf,
        t_perf + t_other,
        100.0 * t_perf / (t_perf + t_other),
    );

    // Alternating timed runs.
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .unwrap();
    for (label, par) in [("1 thread", false), ("N thr x8", true)] {
        let (mut r, mut f) = (Vec::new(), Vec::new());
        for _ in 0..reps {
            for reference in [true, false] {
                let t = if par {
                    time_par(&pool, &d, reference)
                } else {
                    time_serial(&d, &all, reference)
                };
                if reference { &mut r } else { &mut f }.push(t);
            }
        }
        let (rm, rmin, rmax) = stats(&r);
        let (fm, fmin, fmax) = stats(&f);
        let ratios: Vec<f64> = r.iter().zip(&f).map(|(a, b)| a / b).collect();
        let (sm, smin, smax) = stats(&ratios);
        println!(
            "{label:<9} ({}): reference {rm:.3}s [{rmin:.3}-{rmax:.3}], extend2 {fm:.3}s [{fmin:.3}-{fmax:.3}], speedup {sm:.3}x [{smin:.3}-{smax:.3}] over {reps} paired runs",
            if par { threads } else { 1 },
        );
    }
}
