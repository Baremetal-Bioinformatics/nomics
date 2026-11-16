//! Replay a fused-extension dump (`NOMICS_DUMP_FUSED`, format beside
//! `nomics_align::fused::dump`) through the CPU stand-in and compare every
//! result field with the recording.
//!
//! Usage: `fused_bench <fused.bin> [threads]`. With `--features cuda` and
//! `EXTEND_GPU=1`, also runs `GpuFusedExtend` and compares it with the dump.
//! Exits non-zero on any mismatch.

use std::path::Path;
use std::time::Instant;

use nomics_align::fused::{self, dump::Dump};
use nomics_gpu::{GpuExtendResult, GpuFusedResult};
use rayon::prelude::*;

const FIELDS: [&str; 6] = ["score", "qle", "tle", "gtle", "gscore", "max_off"];

fn fields(x: &GpuExtendResult) -> [i32; 6] {
    [x.score, x.qle, x.tle, x.gtle, x.gscore, x.max_off]
}

/// Mismatch counts per field: n_left, n_right, right_h0, then left and
/// right calls by band and field. Calls are compared up to the smaller of
/// the two call counts.
#[derive(Default)]
struct Mismatches {
    n: [u64; 3],
    calls: [[[u64; 6]; 2]; 2],
    jobs: u64,
    first: Option<usize>,
}

impl Mismatches {
    fn compare(want: &[GpuFusedResult], got: &[GpuFusedResult]) -> Self {
        let mut m = Mismatches::default();
        for (j, (w, g)) in want.iter().zip(got).enumerate() {
            let before = m.n.iter().sum::<u64>() + m.calls.iter().flatten().flatten().sum::<u64>();
            m.n[0] += u64::from(w.n_left != g.n_left);
            m.n[1] += u64::from(w.n_right != g.n_right);
            m.n[2] += u64::from(w.right_h0 != g.right_h0);
            for (side, (wc, gc, wn, gn)) in [
                (&w.left, &g.left, w.n_left, g.n_left),
                (&w.right, &g.right, w.n_right, g.n_right),
            ]
            .into_iter()
            .enumerate()
            {
                for i in 0..wn.min(gn).clamp(0, 2) as usize {
                    for (f, (a, b)) in fields(&wc[i]).iter().zip(fields(&gc[i])).enumerate() {
                        m.calls[side][i][f] += u64::from(*a != b);
                    }
                }
            }
            let after = m.n.iter().sum::<u64>() + m.calls.iter().flatten().flatten().sum::<u64>();
            if after > before {
                m.jobs += 1;
                m.first.get_or_insert(j);
            }
        }
        m
    }

    fn report(&self, what: &str) -> bool {
        let mut parts = vec![
            format!("n_left {}", self.n[0]),
            format!("n_right {}", self.n[1]),
            format!("right_h0 {}", self.n[2]),
        ];
        for (side, name) in ["left", "right"].iter().enumerate() {
            for band in 0..2 {
                for (f, field) in FIELDS.iter().enumerate() {
                    parts.push(format!(
                        "{name}[{band}].{field} {}",
                        self.calls[side][band][f]
                    ));
                }
            }
        }
        println!(
            "{what}: {} mismatching jobs{}; {}",
            self.jobs,
            self.first
                .map_or(String::new(), |j| format!(" (first: job {j})")),
            parts.join(", ")
        );
        self.jobs == 0
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let path = args
        .get(1)
        .expect("usage: fused_bench <fused.bin> [threads]");
    let threads: usize = args.get(2).map_or(8, |s| s.parse().unwrap());
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .unwrap();
    let t = Instant::now();
    let bytes = std::fs::read(Path::new(path)).unwrap_or_else(|e| panic!("read {path}: {e}"));
    let d = Dump::parse(bytes);
    let want: Vec<GpuFusedResult> = d.jobs.iter().map(|j| j.result).collect();
    let mut tries = [[0u64; 3]; 2];
    for r in &want {
        tries[0][r.n_left as usize] += 1;
        tries[1][r.n_right as usize] += 1;
    }
    let calls: u64 = want.iter().map(|r| (r.n_left + r.n_right) as u64).sum();
    let cells: u64 = d
        .jobs
        .iter()
        .map(|j| {
            (j.lq.len() * j.lt.len()) as u64 * j.result.n_left as u64
                + (j.rq.len() * j.rt.len()) as u64 * j.result.n_right as u64
        })
        .sum();
    println!(
        "{path}: {} jobs, {calls} calls, {:.2} Gcells (full rectangles), read in {:.2}s; n_left 0/1/2 = {:?}, n_right 0/1/2 = {:?}",
        d.jobs.len(),
        cells as f64 * 1e-9,
        t.elapsed().as_secs_f64(),
        tries[0],
        tries[1]
    );

    let p = d.params;
    let t = Instant::now();
    let cpu: Vec<GpuFusedResult> = d
        .jobs
        .par_iter()
        .map(|j| fused::extend_cpu(&p, &d.input(j)))
        .collect();
    let cpu_s = t.elapsed().as_secs_f64();
    println!(
        "stand-in: {cpu_s:.3}s at {threads} threads ({:.1} M jobs/s)",
        d.jobs.len() as f64 / cpu_s * 1e-6
    );
    let mut ok = Mismatches::compare(&want, &cpu).report("stand-in vs dump");

    if std::env::var_os("EXTEND_GPU").is_some() {
        let mut gpu = nomics_gpu::GpuFusedExtend::new();
        let groups: Vec<_> = d
            .jobs
            .chunks(4096)
            .map(|c| c.iter().map(|j| d.input(j)))
            .collect();
        // Run 0 includes first-touch and module load; later runs are steady.
        let reps: usize = std::env::var("REPS").map_or(3, |s| s.parse().unwrap());
        for rep in 0..reps.max(1) {
            let t = Instant::now();
            match gpu.compute_groups(&p, &groups, threads) {
                Ok((out, timings)) => {
                    let wall = t.elapsed().as_secs_f64();
                    println!(
                        "gpu run {rep}: pack {:.3}s, kernel {:.2} ms, wall {wall:.3}s ({:.2}x the stand-in)",
                        timings.prep.as_secs_f64(),
                        timings.kernel.as_secs_f64() * 1e3,
                        cpu_s / wall
                    );
                    if rep == 0 {
                        ok &= Mismatches::compare(&want, out).report("gpu vs dump");
                    }
                }
                Err(e) => {
                    println!("gpu: {e}");
                    ok = false;
                    break;
                }
            }
        }
    }
    if !ok {
        std::process::exit(1);
    }
}
