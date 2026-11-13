//! Offline Smith-Waterman oracle check and benchmark.
//!
//! Usage: `sw_bench <dump> [reps=3] [threads=all]`. The dump format is
//! documented beside `NOMICS_DUMP_SW` in `assembly/sw.rs`. Every variant's
//! output is checked against the recorded one for every call, then each
//! variant is timed on one thread and on a rayon pool.

use std::io::Read;
use std::time::Instant;

use nomics_hc::assembly::sw_bench::{Variant, run};
use rayon::prelude::*;

struct Call {
    leading_indel: bool,
    w: [i32; 4],
    reference: Vec<u8>,
    alternate: Vec<u8>,
    cigar: Vec<(u32, u8)>,
    offset: i32,
}

struct Cursor<'a>(&'a [u8]);

impl Cursor<'_> {
    fn take(&mut self, n: usize) -> &[u8] {
        let (a, b) = self.0.split_at(n);
        self.0 = b;
        a
    }
    fn u32(&mut self) -> u32 {
        u32::from_le_bytes(self.take(4).try_into().unwrap())
    }
    fn i32(&mut self) -> i32 {
        i32::from_le_bytes(self.take(4).try_into().unwrap())
    }
}

fn load(path: &str) -> Vec<Call> {
    let mut bytes = Vec::new();
    std::fs::File::open(path)
        .expect("open dump")
        .read_to_end(&mut bytes)
        .expect("read dump");
    let mut c = Cursor(&bytes);
    let mut calls = Vec::new();
    while !c.0.is_empty() {
        let leading_indel = c.take(1)[0] == 1;
        let w = [c.i32(), c.i32(), c.i32(), c.i32()];
        let n = c.u32() as usize;
        let reference = c.take(n).to_vec();
        let n = c.u32() as usize;
        let alternate = c.take(n).to_vec();
        let n = c.u32() as usize;
        let cigar = (0..n)
            .map(|_| {
                let len = c.u32();
                (len, c.take(1)[0])
            })
            .collect();
        let offset = c.i32();
        calls.push(Call {
            leading_indel,
            w,
            reference,
            alternate,
            cigar,
            offset,
        });
    }
    calls
}

fn call(v: Variant, c: &Call) -> (Vec<(u32, u8)>, i32) {
    run(v, &c.reference, &c.alternate, c.w, c.leading_indel)
}

fn percentile(sorted: &[u64], p: f64) -> u64 {
    sorted[((sorted.len() - 1) as f64 * p).round() as usize]
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let path = args
        .get(1)
        .expect("usage: sw_bench <dump> [reps] [threads]");
    let reps: usize = args.get(2).map_or(3, |s| s.parse().unwrap());
    let threads: usize = args.get(3).map_or_else(
        || std::thread::available_parallelism().unwrap().get(),
        |s| s.parse().unwrap(),
    );
    let calls = load(path);
    println!("{} calls", calls.len());

    // Call-size distribution, by weight set.
    let mut groups: std::collections::BTreeMap<String, Vec<&Call>> = Default::default();
    for c in &calls {
        let key = format!(
            "{:?} {}",
            c.w,
            if c.leading_indel {
                "LeadingIndel"
            } else {
                "Softclip"
            }
        );
        groups.entry(key).or_default().push(c);
    }
    for (key, g) in &groups {
        let exact = g
            .iter()
            .filter(|c| {
                !c.leading_indel
                    && c.cigar.len() == 1
                    && c.cigar[0].1 == 0
                    && c.cigar[0].0 as usize == c.alternate.len()
                    && c.reference
                        .windows(c.alternate.len())
                        .any(|x| x == c.alternate)
            })
            .count();
        let mut r: Vec<u64> = g.iter().map(|c| c.reference.len() as u64 + 1).collect();
        let mut a: Vec<u64> = g.iter().map(|c| c.alternate.len() as u64 + 1).collect();
        let mut cells: Vec<u64> = g
            .iter()
            .map(|c| (c.reference.len() as u64 + 1) * (c.alternate.len() as u64 + 1))
            .collect();
        let total: u64 = cells.iter().sum();
        r.sort_unstable();
        a.sort_unstable();
        cells.sort_unstable();
        let q = |v: &[u64]| {
            format!(
                "p10 {} p50 {} p90 {} p99 {} max {}",
                percentile(v, 0.1),
                percentile(v, 0.5),
                percentile(v, 0.9),
                percentile(v, 0.99),
                v[v.len() - 1]
            )
        };
        println!("{key}: {} calls ({exact} substring hits, no DP)", g.len());
        println!("  rows  {}", q(&r));
        println!("  cols  {}", q(&a));
        println!("  cells {}  total {:.3} G", q(&cells), total as f64 / 1e9);
    }

    let variants = [
        Variant::Reference,
        Variant::WavefrontPortable,
        Variant::Wavefront,
    ];
    for v in variants {
        let bad = calls
            .par_iter()
            .filter(|c| call(v, c) != (c.cigar.clone(), c.offset))
            .count();
        println!("{v:?}: {bad} mismatches / {}", calls.len());
    }

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .unwrap();
    let fmt = |x: &[f64]| {
        x.iter()
            .map(|t| format!("{t:.3}"))
            .collect::<Vec<_>>()
            .join(" ")
    };
    for (key, g) in &groups {
        println!("timing {key} ({} calls)", g.len());
        for v in variants {
            let mut one = Vec::new();
            let mut many = Vec::new();
            for _ in 0..reps {
                let t = Instant::now();
                let mut sink = 0_usize;
                for c in g {
                    sink = sink.wrapping_add(call(v, c).0.len());
                }
                std::hint::black_box(sink);
                one.push(t.elapsed().as_secs_f64());
                let t = Instant::now();
                let sink: usize = pool.install(|| g.par_iter().map(|c| call(v, c).0.len()).sum());
                std::hint::black_box(sink);
                many.push(t.elapsed().as_secs_f64());
            }
            println!(
                "  {v:?}: 1 thread [{}] s; {threads} threads [{}] s",
                fmt(&one),
                fmt(&many)
            );
        }
    }
}
