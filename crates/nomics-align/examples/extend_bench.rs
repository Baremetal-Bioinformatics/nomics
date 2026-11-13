//! Benchmark scalar and across-call lane-parallel `ksw_extend2` on a dump.
//!
//! Usage: `extend_bench <extend2.bin> [threads]`. With `--features cuda` and
//! `EXTEND_GPU=1`, compares the GPU kernel with the dump and the scalar CPU.

use std::hint::black_box;
use std::ops::Range;
use std::path::Path;
use std::time::{Duration, Instant};

use nomics_align::ksw::{self, Extend};
use nomics_align::ksw_lanes::{Extend2Call, LaneStats, extend2_batch};
use rayon::prelude::*;

const ZERO: Extend = Extend {
    score: 0,
    qle: 0,
    tle: 0,
    gtle: 0,
    gscore: 0,
    max_off: 0,
};

#[derive(Clone)]
struct Record {
    query: Range<usize>,
    target: Range<usize>,
    m: i32,
    matrix: usize,
    args: [i32; 8],
    expected: Extend,
}

struct Dump {
    bytes: Vec<u8>,
    matrices: Vec<Box<[i8]>>,
    records: Vec<Record>,
}

impl Dump {
    fn read(path: &Path) -> Self {
        let bytes = std::fs::read(path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
        let mut at = 0usize;
        let mut matrices = Vec::<Box<[i8]>>::new();
        let mut records = Vec::new();
        while at < bytes.len() {
            let qlen = take_u32(&bytes, &mut at) as usize;
            let query = take_range(&bytes, &mut at, qlen);
            let tlen = take_u32(&bytes, &mut at) as usize;
            let target = take_range(&bytes, &mut at, tlen);
            let m = take_i32(&bytes, &mut at);
            assert!(m > 0, "record {} has m={m}", records.len());
            let mat_len = (m as usize)
                .checked_mul(m as usize)
                .expect("matrix length overflow");
            let mat_range = take_range(&bytes, &mut at, mat_len);
            let matrix_values: Vec<i8> = bytes[mat_range].iter().map(|&x| x as i8).collect();
            let matrix = matrices
                .iter()
                .position(|known| known.as_ref() == matrix_values)
                .unwrap_or_else(|| {
                    matrices.push(matrix_values.into_boxed_slice());
                    matrices.len() - 1
                });
            let args = std::array::from_fn(|_| take_i32(&bytes, &mut at));
            let result: [i32; 6] = std::array::from_fn(|_| take_i32(&bytes, &mut at));
            records.push(Record {
                query,
                target,
                m,
                matrix,
                args,
                expected: Extend {
                    score: result[0],
                    qle: result[1],
                    tle: result[2],
                    gtle: result[3],
                    gscore: result[4],
                    max_off: result[5],
                },
            });
        }
        assert_eq!(at, bytes.len());
        Self {
            bytes,
            matrices,
            records,
        }
    }

    fn calls<'a>(&'a self, order: &[usize]) -> Vec<Extend2Call<'a>> {
        order
            .iter()
            .map(|&index| {
                let record = &self.records[index];
                Extend2Call {
                    query: &self.bytes[record.query.clone()],
                    target: &self.bytes[record.target.clone()],
                    m: record.m,
                    mat: &self.matrices[record.matrix],
                    o_del: record.args[0],
                    e_del: record.args[1],
                    o_ins: record.args[2],
                    e_ins: record.args[3],
                    w: record.args[4],
                    end_bonus: record.args[5],
                    zdrop: record.args[6],
                    h0: record.args[7],
                }
            })
            .collect()
    }

    fn expected(&self, order: &[usize]) -> Vec<Extend> {
        order
            .iter()
            .map(|&index| self.records[index].expected)
            .collect()
    }
}

fn take_range(bytes: &[u8], at: &mut usize, len: usize) -> Range<usize> {
    let end = at.checked_add(len).expect("dump offset overflow");
    assert!(end <= bytes.len(), "truncated dump at byte {at}");
    let range = *at..end;
    *at = end;
    range
}

fn take_u32(bytes: &[u8], at: &mut usize) -> u32 {
    let range = take_range(bytes, at, 4);
    u32::from_le_bytes(bytes[range].try_into().unwrap())
}

fn take_i32(bytes: &[u8], at: &mut usize) -> i32 {
    take_u32(bytes, at) as i32
}

fn scalar(call: &Extend2Call<'_>) -> Extend {
    ksw::extend2(
        call.query.len() as i32,
        call.query,
        call.target.len() as i32,
        call.target,
        call.m,
        call.mat,
        call.o_del,
        call.e_del,
        call.o_ins,
        call.e_ins,
        call.w,
        call.end_bonus,
        call.zdrop,
        call.h0,
    )
}

fn run_scalar(calls: &[Extend2Call<'_>]) -> (Duration, Vec<Extend>) {
    let mut out = vec![ZERO; calls.len()];
    let start = Instant::now();
    for (call, result) in calls.iter().zip(&mut out) {
        *result = black_box(scalar(black_box(call)));
    }
    (start.elapsed(), out)
}

fn run_lanes(calls: &[Extend2Call<'_>]) -> (Duration, Vec<Extend>, LaneStats) {
    let mut out = vec![ZERO; calls.len()];
    let start = Instant::now();
    let stats = black_box(extend2_batch(black_box(calls), &mut out));
    (start.elapsed(), out, stats)
}

fn run_scalar_threads(
    pool: &rayon::ThreadPool,
    calls: &[Extend2Call<'_>],
) -> (Duration, Vec<Extend>) {
    let mut out = vec![ZERO; calls.len()];
    let start = Instant::now();
    pool.install(|| {
        calls
            .par_iter()
            .zip(out.par_iter_mut())
            .for_each(|(call, result)| *result = black_box(scalar(black_box(call))));
    });
    (start.elapsed(), out)
}

fn run_lanes_threads(
    pool: &rayon::ThreadPool,
    calls: &[Extend2Call<'_>],
    threads: usize,
) -> (Duration, Vec<Extend>, LaneStats) {
    let mut out = vec![ZERO; calls.len()];
    let chunk = calls
        .len()
        .div_ceil(threads * 4)
        .div_ceil(nomics_align::ksw_lanes::LANES)
        * nomics_align::ksw_lanes::LANES;
    let start = Instant::now();
    let stats = pool.install(|| {
        calls
            .par_chunks(chunk)
            .zip(out.par_chunks_mut(chunk))
            .map(|(call_chunk, out_chunk)| {
                black_box(extend2_batch(black_box(call_chunk), out_chunk))
            })
            .reduce(LaneStats::default, |mut total, part| {
                total.add(part);
                total
            })
    });
    let elapsed = start.elapsed();
    (elapsed, out, stats)
}

fn verify(label: &str, actual: &[Extend], expected: &[Extend]) {
    assert_eq!(actual.len(), expected.len());
    for (index, (&actual, &expected)) in actual.iter().zip(expected).enumerate() {
        assert_eq!(
            actual, expected,
            "{label} differs from dump at call {index}"
        );
    }
    eprintln!(
        "verified {label}: {} / {} calls",
        actual.len(),
        expected.len()
    );
}

fn report(label: &str, elapsed: Duration, calls: usize, cells: u64, baseline: Duration) {
    let seconds = elapsed.as_secs_f64();
    println!(
        "{label:<18} calls {calls:>9}  cells {cells:>12}  wall {seconds:>8.3}s  {ns:>9.1} ns/call  {mcells:>8.1} Mcell/s  {speedup:>6.2}x vs scalar-1",
        ns = elapsed.as_nanos() as f64 / calls as f64,
        mcells = cells as f64 / seconds / 1e6,
        speedup = baseline.as_secs_f64() / seconds,
    );
}

fn main() {
    let args: Vec<_> = std::env::args_os().collect();
    assert!(
        (2..=3).contains(&args.len()),
        "usage: extend_bench <extend2.bin> [threads]"
    );
    let threads = args
        .get(2)
        .map_or_else(
            || std::thread::available_parallelism().map_or(1, usize::from),
            |s| s.to_string_lossy().parse().expect("thread count"),
        )
        .max(1);
    let load_start = Instant::now();
    let dump = Dump::read(Path::new(&args[1]));
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .expect("build rayon pool");
    let order: Vec<_> = (0..dump.records.len()).collect();
    let calls = dump.calls(&order);
    let expected = dump.expected(&order);
    eprintln!(
        "loaded {} calls, {} distinct matrices, {:.1} MiB in {:.3}s",
        calls.len(),
        dump.matrices.len(),
        dump.bytes.len() as f64 / (1 << 20) as f64,
        load_start.elapsed().as_secs_f64()
    );

    #[cfg(feature = "cuda")]
    if std::env::var_os("EXTEND_GPU").is_some() {
        // Empty queries return before any DP (the CPU handles them).
        let keep: Vec<usize> = (0..calls.len())
            .filter(|&i| !calls[i].query.is_empty())
            .collect();
        let calls: Vec<_> = keep.iter().map(|&i| calls[i]).collect();
        let expected: Vec<_> = keep.iter().map(|&i| expected[i]).collect();
        gpu(&calls, &expected, &pool);
        return;
    }

    let warm = calls.len().min(1024);
    let _ = run_scalar(&calls[..warm]);
    let _ = run_lanes(&calls[..warm]);

    let (scalar_one_time, scalar_one) = run_scalar(&calls);
    verify("scalar-1", &scalar_one, &expected);
    let (scalar_n_time, scalar_n) = run_scalar_threads(&pool, &calls);
    verify("scalar-N", &scalar_n, &expected);
    let (lanes_one_time, lanes_one, lanes_one_stats) = run_lanes(&calls);
    verify("lanes-1", &lanes_one, &expected);
    let (lanes_n_time, lanes_n, lanes_n_stats) = run_lanes_threads(&pool, &calls, threads);
    verify("lanes-N", &lanes_n, &expected);

    let cells = lanes_one_stats.useful_cells;
    println!("\noriginal dump order ({threads} worker threads):");
    report(
        "scalar-1",
        scalar_one_time,
        calls.len(),
        cells,
        scalar_one_time,
    );
    report(
        "scalar-N",
        scalar_n_time,
        calls.len(),
        cells,
        scalar_one_time,
    );
    report(
        "lanes-1",
        lanes_one_time,
        calls.len(),
        cells,
        scalar_one_time,
    );
    report("lanes-N", lanes_n_time, calls.len(), cells, scalar_one_time);
    println!(
        "lane utilization: {:.2}% ({} useful cells / {} lane slots); same-thread lane speedups: 1T {:.2}x, {threads}T {:.2}x",
        lanes_one_stats.utilization() * 100.0,
        cells,
        lanes_one_stats.lane_slots(),
        scalar_one_time.as_secs_f64() / lanes_one_time.as_secs_f64(),
        scalar_n_time.as_secs_f64() / lanes_n_time.as_secs_f64(),
    );
    assert_eq!(lanes_one_stats.useful_cells, lanes_n_stats.useful_cells);

    // Stable sorting does not alter calls; it only puts similarly shaped
    // problems in the same lane group to isolate divergence overhead.
    let mut sorted_order = order;
    sorted_order.sort_by_key(|&index| {
        let r = &dump.records[index];
        (r.query.len(), r.target.len(), r.args[4], r.args[6])
    });
    let sorted_calls = dump.calls(&sorted_order);
    let sorted_expected = dump.expected(&sorted_order);
    let (sorted_one_time, sorted_one, sorted_one_stats) = run_lanes(&sorted_calls);
    verify("lanes-1 sorted", &sorted_one, &sorted_expected);
    let (sorted_n_time, sorted_n, sorted_n_stats) =
        run_lanes_threads(&pool, &sorted_calls, threads);
    verify("lanes-N sorted", &sorted_n, &sorted_expected);

    println!("\nlength/band/z-drop sorted order:");
    report(
        "lanes-1 sorted",
        sorted_one_time,
        calls.len(),
        cells,
        scalar_one_time,
    );
    report(
        "lanes-N sorted",
        sorted_n_time,
        calls.len(),
        cells,
        scalar_one_time,
    );
    println!(
        "lane utilization: {:.2}% ({} lane slots, {:+.2} points vs original); lane-slot reduction {:.2}%",
        sorted_one_stats.utilization() * 100.0,
        sorted_one_stats.lane_slots(),
        (sorted_one_stats.utilization() - lanes_one_stats.utilization()) * 100.0,
        (1.0 - sorted_one_stats.lane_slots() as f64 / lanes_one_stats.lane_slots() as f64) * 100.0,
    );
    assert_eq!(sorted_one_stats.useful_cells, cells);
    assert_eq!(sorted_one_stats.useful_cells, sorted_n_stats.useful_cells);
    black_box((
        scalar_one, scalar_n, lanes_one, lanes_n, sorted_one, sorted_n,
    ));
}

#[cfg(feature = "cuda")]
fn gpu(calls: &[Extend2Call<'_>], expected: &[Extend], pool: &rayon::ThreadPool) {
    use nomics_gpu::{GpuExtendInput, GpuExtendResult, GpuKswExtend2};
    let (cpu_time, cpu) = run_scalar_threads(pool, calls);
    verify("scalar-N", &cpu, expected);
    let cells: u64 = calls
        .iter()
        .map(|c| c.query.len() as u64 * c.target.len() as u64)
        .sum();
    let inputs: Vec<GpuExtendInput> = calls
        .iter()
        .map(|c| GpuExtendInput {
            query: c.query,
            target: c.target,
            m: c.m,
            mat: c.mat,
            o_del: c.o_del,
            e_del: c.e_del,
            o_ins: c.o_ins,
            e_ins: c.e_ins,
            w: c.w,
            end_bonus: c.end_bonus,
            zdrop: c.zdrop,
            h0: c.h0,
        })
        .collect();
    let mut engine = GpuKswExtend2::new();
    let mut out = vec![GpuExtendResult::default(); inputs.len()];
    let _ = engine.compute(&inputs[..inputs.len().min(4096)], &mut out);
    let started = Instant::now();
    let t = engine.compute(&inputs, &mut out).expect("GPU extend2");
    let total = started.elapsed();
    let mut bad = 0usize;
    for (i, (g, e)) in out.iter().zip(expected).enumerate() {
        let same = g.score == e.score
            && g.qle == e.qle
            && g.tle == e.tle
            && g.gtle == e.gtle
            && g.gscore == e.gscore
            && g.max_off == e.max_off;
        if !same {
            bad += 1;
            if bad <= 5 {
                println!("  mismatch {i}: gpu {g:?} cpu {e:?}");
            }
        }
    }
    println!(
        "CPU: {:.3}s, {:.3e} full-matrix cells/s; GPU: total {:.3}s (prep {:.3}s, kernel {:.3}s), {:.3e} cells/s kernel; {bad} mismatches of {}",
        cpu_time.as_secs_f64(),
        cells as f64 / cpu_time.as_secs_f64(),
        total.as_secs_f64(),
        t.prep.as_secs_f64(),
        t.kernel.as_secs_f64(),
        cells as f64 / t.kernel.as_secs_f64(),
        calls.len()
    );
}
