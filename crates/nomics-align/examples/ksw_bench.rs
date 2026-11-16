//! Offline CPU/GPU experiment for BWA 0.7.19 `ksw_align2`.
//!
//! Usage: `ksw_bench <align2.bin> [threads=20] [gpu_batch_calls=8192]`.
//! The dump format is documented beside `NOMICS_DUMP_KSW` in `ksw.rs`.

use std::hint::black_box;
use std::ops::Range;
use std::path::Path;
use std::time::{Duration, Instant};

use nomics_align::ksw::{self, KSW_XBYTE, Kswr, Profile};
use rayon::prelude::*;

#[cfg(any(feature = "cuda", feature = "metal"))]
use nomics_gpu::{GpuKswAlign2, GpuKswAlign2Input, GpuKswResult, KswTimings};

#[derive(Clone)]
struct Record {
    query: Range<usize>,
    target: Range<usize>,
    m: i32,
    matrix: usize,
    args: [i32; 5],
    size: u8,
    expected: Kswr,
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
            let size = take_i32(&bytes, &mut at);
            assert!(size == 1 || size == 2, "record has score size {size}");
            let result: [i32; 7] = std::array::from_fn(|_| take_i32(&bytes, &mut at));
            records.push(Record {
                query,
                target,
                m,
                matrix,
                args,
                size: size as u8,
                expected: Kswr {
                    score: result[0],
                    te: result[1],
                    qe: result[2],
                    score2: result[3],
                    te2: result[4],
                    tb: result[5],
                    qb: result[6],
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

    fn query(&self, record: &Record) -> &[u8] {
        &self.bytes[record.query.clone()]
    }

    fn target(&self, record: &Record) -> &[u8] {
        &self.bytes[record.target.clone()]
    }

    #[cfg(any(feature = "cuda", feature = "metal"))]
    fn gpu_input(&self, record: &Record) -> GpuKswAlign2Input<'_> {
        GpuKswAlign2Input {
            query: self.query(record),
            target: self.target(record),
            m: record.m,
            mat: &self.matrices[record.matrix],
            o_del: record.args[0],
            e_del: record.args[1],
            o_ins: record.args[2],
            e_ins: record.args[3],
            xtra: record.args[4],
            size: record.size,
        }
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

fn scalar(dump: &Dump, record: &Record) -> Kswr {
    let query = dump.query(record);
    let target = dump.target(record);
    let mat = &dump.matrices[record.matrix];
    let xtra_size = if record.args[4] & KSW_XBYTE != 0 {
        1
    } else {
        2
    };
    if record.size == xtra_size {
        ksw::align2(
            query.len() as i32,
            query,
            target.len() as i32,
            target,
            record.m,
            mat,
            record.args[0],
            record.args[1],
            record.args[2],
            record.args[3],
            record.args[4],
            None,
        )
    } else {
        let mut profile = Profile::new(record.size, query.len() as i32, query, record.m, mat);
        ksw::align2(
            query.len() as i32,
            query,
            target.len() as i32,
            target,
            record.m,
            mat,
            record.args[0],
            record.args[1],
            record.args[2],
            record.args[3],
            record.args[4],
            Some(&mut profile),
        )
    }
}

fn cpu(dump: &Dump, indices: &[usize], threads: usize) -> (Vec<Kswr>, Duration) {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .expect("build rayon pool");
    let start = Instant::now();
    let results = pool.install(|| {
        indices
            .par_iter()
            .map(|&index| black_box(scalar(black_box(dump), &dump.records[index])))
            .collect()
    });
    (results, start.elapsed())
}

fn mismatch_examples(
    label: &str,
    indices: &[usize],
    actual: &[Kswr],
    expected: &[Kswr],
) -> (usize, Vec<String>) {
    let mut mismatches = 0usize;
    let mut examples = Vec::new();
    for ((&index, &actual), &expected) in indices.iter().zip(actual).zip(expected) {
        if actual != expected {
            mismatches += 1;
            if examples.len() < 8 {
                examples.push(format!(
                    "{label} call {index}: expected {expected:?}, actual {actual:?}"
                ));
            }
        }
    }
    (mismatches, examples)
}

fn cells(dump: &Dump, indices: &[usize]) -> u64 {
    indices
        .iter()
        .map(|&index| {
            let record = &dump.records[index];
            record.query.len() as u64 * record.target.len() as u64
        })
        .sum()
}

fn rate(cells: u64, elapsed: Duration) -> f64 {
    cells as f64 / elapsed.as_secs_f64()
}

#[cfg(any(feature = "cuda", feature = "metal"))]
fn from_gpu(result: GpuKswResult) -> Kswr {
    Kswr {
        score: result.score,
        te: result.te,
        qe: result.qe,
        score2: result.score2,
        te2: result.te2,
        tb: result.tb,
        qb: result.qb,
    }
}

#[cfg(any(feature = "cuda", feature = "metal"))]
fn gpu(
    dump: &Dump,
    indices: &[usize],
    expected_all: &[Kswr],
    batch_calls: usize,
    show_device: bool,
) -> Result<(KswTimings, usize, Vec<String>), nomics_gpu::GpuError> {
    if show_device {
        let info = nomics_gpu::device_info(0)?;
        // CUDA: sm_XY and SMs; Metal: AppleN family and GPU cores.
        let (arch, units) = if cfg!(feature = "metal") {
            (format!("Apple{}", info.compute_capability.0), "GPU cores")
        } else {
            let (major, minor) = info.compute_capability;
            (format!("sm_{major}{minor}"), "SMs")
        };
        println!(
            "GPU: {} {arch} ({} {units}, pageable host access: {})",
            info.name, info.multiprocessors, info.pageable_memory_access
        );
        assert!(
            info.pageable_memory_access,
            "KSW experiment requires pageable host memory access"
        );
    }
    let mut engine = GpuKswAlign2::new();
    let mut total = KswTimings::default();
    let mut mismatches = 0usize;
    let mut examples = Vec::new();
    for batch in indices.chunks(batch_calls) {
        let build_start = Instant::now();
        let inputs: Vec<_> = batch
            .iter()
            .map(|&index| dump.gpu_input(&dump.records[index]))
            .collect();
        let mut output = vec![GpuKswResult::default(); batch.len()];
        total.prep += build_start.elapsed();
        let call_start = Instant::now();
        let timing = engine.compute(&inputs, &mut output)?;
        let call_wall = call_start.elapsed();
        let accounted = timing.prep + timing.kernel + timing.post;
        total.prep += timing.prep + call_wall.saturating_sub(accounted);
        total.kernel += timing.kernel;
        total.post += timing.post;
        for (&index, actual) in batch.iter().zip(output.into_iter().map(from_gpu)) {
            let expected = expected_all[index];
            if actual != expected {
                mismatches += 1;
                if examples.len() < 8 {
                    examples.push(format!(
                        "GPU call {index}: CPU {expected:?}, GPU {actual:?}"
                    ));
                }
            }
        }
    }
    Ok((total, mismatches, examples))
}

#[cfg(any(feature = "cuda", feature = "metal"))]
fn hybrid(
    dump: &Dump,
    expected: &[Kswr],
    frac: f64,
    threads: usize,
    batch_calls: usize,
) -> (Duration, Duration, Duration, usize) {
    assert!((0.0..=1.0).contains(&frac));
    let (mut accumulated, mut gpu_cells) = (0u64, 0u64);
    let (mut gpu_indices, mut cpu_indices) = (Vec::new(), Vec::new());
    for index in 0..dump.records.len() {
        let call_cells = cells(dump, std::slice::from_ref(&index));
        accumulated += call_cells;
        if (gpu_cells as f64) < frac * accumulated as f64 {
            gpu_cells += call_cells;
            gpu_indices.push(index);
        } else {
            cpu_indices.push(index);
        }
    }
    let start = Instant::now();
    let (cpu_elapsed, gpu_elapsed, mismatches) = std::thread::scope(|scope| {
        let gpu_thread = scope.spawn(|| {
            let gpu_start = Instant::now();
            let (_, mismatches, _) =
                gpu(dump, &gpu_indices, expected, batch_calls, false).expect("run hybrid CUDA KSW");
            (gpu_start.elapsed(), mismatches)
        });
        let (cpu_results, cpu_elapsed) = cpu(dump, &cpu_indices, threads);
        let cpu_expected: Vec<_> = cpu_indices.iter().map(|&index| expected[index]).collect();
        let (cpu_mismatches, _) =
            mismatch_examples("hybrid CPU", &cpu_indices, &cpu_results, &cpu_expected);
        let (gpu_elapsed, gpu_mismatches) = gpu_thread.join().unwrap();
        (cpu_elapsed, gpu_elapsed, cpu_mismatches + gpu_mismatches)
    });
    (start.elapsed(), cpu_elapsed, gpu_elapsed, mismatches)
}

fn main() {
    let args: Vec<_> = std::env::args_os().collect();
    if args.len() < 2 || args.len() > 4 {
        eprintln!("usage: ksw_bench <align2.bin> [threads=20] [gpu_batch_calls=8192]");
        std::process::exit(2);
    }
    let threads = args
        .get(2)
        .map(|s| {
            s.to_string_lossy()
                .parse()
                .expect("threads must be an integer")
        })
        .unwrap_or(20usize);
    let batch_calls = args
        .get(3)
        .map(|s| {
            s.to_string_lossy()
                .parse()
                .expect("gpu_batch_calls must be an integer")
        })
        .unwrap_or(8192usize);
    assert!(threads > 0);
    assert!(batch_calls > 0);

    let load_start = Instant::now();
    let dump = Dump::read(Path::new(&args[1]));
    let load_elapsed = load_start.elapsed();
    let indices: Vec<_> = (0..dump.records.len()).collect();
    let total_cells = cells(&dump, &indices);
    let max_query = dump
        .records
        .iter()
        .map(|r| r.query.len())
        .max()
        .unwrap_or(0);
    let max_target = dump
        .records
        .iter()
        .map(|r| r.target.len())
        .max()
        .unwrap_or(0);
    println!(
        "loaded {} calls in {:.3}s; {total_cells} cells; max query {max_query}, max target {max_target}",
        dump.records.len(),
        load_elapsed.as_secs_f64()
    );

    let (cpu_results, cpu_elapsed) = cpu(&dump, &indices, threads);
    let dumped: Vec<_> = dump.records.iter().map(|record| record.expected).collect();
    let (cpu_mismatches, cpu_examples) =
        mismatch_examples("CPU/dump", &indices, &cpu_results, &dumped);
    println!(
        "CPU ({threads} threads): {:.6}s, {:.3e} cells/s; {cpu_mismatches} mismatches vs dump",
        cpu_elapsed.as_secs_f64(),
        rate(total_cells, cpu_elapsed)
    );
    for example in cpu_examples {
        println!("  {example}");
    }
    assert_eq!(cpu_mismatches, 0, "CPU rerun differs from dumped results");

    #[cfg(any(feature = "cuda", feature = "metal"))]
    {
        let (timing, mismatches, examples) =
            gpu(&dump, &indices, &cpu_results, batch_calls, true).expect("run CUDA KSW");
        let total = timing.prep + timing.kernel + timing.post;
        println!(
            "GPU: {:.6}s total (prep {:.6}s, kernel {:.6}s, post {:.6}s), {:.3e} cells/s total, {:.3e} cells/s kernel",
            total.as_secs_f64(),
            timing.prep.as_secs_f64(),
            timing.kernel.as_secs_f64(),
            timing.post.as_secs_f64(),
            rate(total_cells, total),
            rate(total_cells, timing.kernel)
        );
        println!(
            "exactness: {mismatches} field-complete mismatches out of {} calls",
            dump.records.len()
        );
        for example in examples {
            println!("  {example}");
        }
        // NOMICS_HYBRID="frac:threads,...", e.g. "0.3:20,0.4:19".
        if let Ok(spec) = std::env::var("NOMICS_HYBRID") {
            for item in spec.split(',') {
                let (fraction, hybrid_threads) = item.split_once(':').unwrap_or((item, "20"));
                let fraction: f64 = fraction.parse().expect("hybrid GPU fraction");
                let hybrid_threads: usize =
                    hybrid_threads.parse().expect("hybrid CPU thread count");
                let (wall, cpu_part, gpu_part, hybrid_mismatches) =
                    hybrid(&dump, &cpu_results, fraction, hybrid_threads, batch_calls);
                println!(
                    "hybrid gpu {:.0}% cpu {hybrid_threads} threads: {:.3}s wall (cpu part {:.3}s, gpu part {:.3}s), {:.3e} cells/s, {hybrid_mismatches} mismatches",
                    fraction * 100.0,
                    wall.as_secs_f64(),
                    cpu_part.as_secs_f64(),
                    gpu_part.as_secs_f64(),
                    rate(total_cells, wall)
                );
            }
        }
    }

    #[cfg(not(any(feature = "cuda", feature = "metal")))]
    {
        let _ = (cpu_results, batch_calls);
        println!("GPU: unavailable (build with --features cuda on the Spark)");
    }
}
