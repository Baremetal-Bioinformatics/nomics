//! Offline PairHMM CPU/GPU experiment.
//!
//! Usage: `pairhmm_bench <dump.bin> [threads=20] [gpu_batch_pairs=65536]`.
//! The dump format is documented beside `NOMICS_DUMP_PAIRHMM` in
//! `likelihood.rs`.

use std::fs::File;
use std::io::{self, BufReader};
use std::path::Path;
use std::time::Instant;

use nomics_hc::pairhmm::{LANES, PairHmmLanes};
use rayon::prelude::*;

#[cfg(feature = "cuda")]
use std::time::Duration;

#[cfg(feature = "cuda")]
use nomics_gpu::{GpuPairHmm, GpuPairHmmPair};
#[cfg(feature = "cuda")]
use nomics_hc::pairhmm::{initial_condition, row_params};

#[derive(Clone)]
struct Read {
    bases: Vec<u8>,
    quals: Vec<u8>,
    ins: Vec<u8>,
    del: Vec<u8>,
}

#[derive(Clone)]
struct Record {
    haps: Vec<Vec<u8>>,
    reads: Vec<Read>,
}

fn next_u32(reader: &mut impl io::Read) -> io::Result<Option<u32>> {
    let mut bytes = [0u8; 4];
    if reader.read(&mut bytes[..1])? == 0 {
        return Ok(None);
    }
    reader.read_exact(&mut bytes[1..])?;
    Ok(Some(u32::from_le_bytes(bytes)))
}

fn required_u32(reader: &mut impl io::Read) -> io::Result<u32> {
    next_u32(reader)?.ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))
}

fn bytes(reader: &mut impl io::Read, len: usize) -> io::Result<Vec<u8>> {
    let mut value = vec![0; len];
    reader.read_exact(&mut value)?;
    Ok(value)
}

fn load(path: &Path) -> io::Result<Vec<Record>> {
    let mut reader = BufReader::with_capacity(8 << 20, File::open(path)?);
    let mut records = Vec::new();
    while let Some(n_haps) = next_u32(&mut reader)? {
        let mut haps = Vec::with_capacity(n_haps as usize);
        for _ in 0..n_haps {
            let len = required_u32(&mut reader)? as usize;
            haps.push(bytes(&mut reader, len)?);
        }
        let n_reads = required_u32(&mut reader)? as usize;
        let mut reads = Vec::with_capacity(n_reads);
        for _ in 0..n_reads {
            let len = required_u32(&mut reader)? as usize;
            reads.push(Read {
                bases: bytes(&mut reader, len)?,
                quals: bytes(&mut reader, len)?,
                ins: bytes(&mut reader, len)?,
                del: bytes(&mut reader, len)?,
            });
        }
        records.push(Record { haps, reads });
    }
    Ok(records)
}

fn stats(records: &[Record]) -> (u64, u64, usize, usize) {
    let mut pairs = 0u64;
    let mut cells = 0u64;
    let mut max_read = 0usize;
    let mut max_hap = 0usize;
    for record in records {
        pairs += (record.reads.len() * record.haps.len()) as u64;
        for read in &record.reads {
            max_read = max_read.max(read.bases.len());
            for hap in &record.haps {
                max_hap = max_hap.max(hap.len());
                cells += (read.bases.len() * hap.len()) as u64;
            }
        }
    }
    (pairs, cells, max_read, max_hap)
}

/// The prefix-sharing traversal, `values[h][r]` flattened read-major per
/// record like `cpu`.
fn cpu_trie(records: &[Record], threads: usize) -> (Vec<Vec<f64>>, f64) {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .expect("build rayon pool");
    let start = Instant::now();
    let results = pool.install(|| {
        records
            .par_iter()
            .map(|record| {
                let gcps: Vec<Vec<u8>> = record
                    .reads
                    .iter()
                    .map(|r| vec![10; r.bases.len()])
                    .collect();
                let reads: Vec<nomics_hc::pairhmm_trie::HmmRead> = record
                    .reads
                    .iter()
                    .zip(&gcps)
                    .map(|(r, gcp)| nomics_hc::pairhmm_trie::HmmRead {
                        bases: &r.bases,
                        quals: &r.quals,
                        ins: &r.ins,
                        del: &r.del,
                        gcp,
                    })
                    .collect();
                let haps: Vec<&[u8]> = record.haps.iter().map(|h| h.as_slice()).collect();
                let v = nomics_hc::pairhmm_trie::likelihoods(&haps, &reads);
                let mut flat = Vec::with_capacity(reads.len() * haps.len());
                for r in 0..reads.len() {
                    for h in 0..haps.len() {
                        flat.push(v[h][r]);
                    }
                }
                flat
            })
            .collect()
    });
    (results, start.elapsed().as_secs_f64())
}

fn cpu(records: &[Record], threads: usize, portable: bool) -> (Vec<Vec<f64>>, f64) {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .expect("build rayon pool");
    let start = Instant::now();
    let results = pool.install(|| {
        records
            .par_iter()
            .map_init(PairHmmLanes::default, |hmm, record| {
                let mut values = Vec::with_capacity(record.reads.len() * record.haps.len());
                for read in &record.reads {
                    let gcp = vec![10; read.bases.len()];
                    for chunk in record.haps.chunks(LANES) {
                        let mut hap_refs: [&[u8]; LANES] = [&[]; LANES];
                        for (slot, hap) in hap_refs.iter_mut().zip(chunk) {
                            *slot = hap;
                        }
                        let mut out = [0.0; LANES];
                        let f = if portable {
                            PairHmmLanes::log10_likelihoods_portable
                        } else {
                            PairHmmLanes::log10_likelihoods
                        };
                        f(
                            hmm,
                            &hap_refs[..chunk.len()],
                            &read.bases,
                            &read.quals,
                            &read.ins,
                            &read.del,
                            &gcp,
                            &mut out[..chunk.len()],
                        );
                        values.extend_from_slice(&out[..chunk.len()]);
                    }
                }
                values
            })
            .collect()
    });
    (results, start.elapsed().as_secs_f64())
}

#[cfg(feature = "cuda")]
struct Packed {
    read_bases: Vec<u8>,
    hap_bases: Vec<u8>,
    params: Vec<[f64; 7]>,
    pairs: Vec<GpuPairHmmPair>,
}

#[cfg(feature = "cuda")]
fn pack(records: &[Record], init: f64) -> Packed {
    let mut packed = Packed {
        read_bases: Vec::new(),
        hap_bases: Vec::new(),
        params: Vec::new(),
        pairs: Vec::new(),
    };
    for record in records {
        let mut hap_offsets = Vec::with_capacity(record.haps.len());
        for hap in &record.haps {
            hap_offsets.push(packed.hap_bases.len() as u64);
            packed.hap_bases.extend_from_slice(hap);
        }
        for read in &record.reads {
            let read_offset = packed.read_bases.len() as u64;
            let params_offset = packed.params.len() as u64;
            packed.read_bases.extend_from_slice(&read.bases);
            packed.params.extend(
                read.ins
                    .iter()
                    .zip(&read.del)
                    .zip(&read.quals)
                    .map(|((&ins, &del), &q)| row_params(ins, del, q, 10)),
            );
            for (hap, &hap_offset) in record.haps.iter().zip(&hap_offsets) {
                packed.pairs.push(GpuPairHmmPair {
                    read_offset,
                    hap_offset,
                    params_offset,
                    read_len: read.bases.len() as u32,
                    hap_len: hap.len() as u32,
                    initial: init / hap.len() as f64,
                });
            }
        }
    }
    packed
}

#[cfg(feature = "cuda")]
fn gpu(
    records: &[Record],
    cpu_results: &[Vec<f64>],
    batch_pairs: usize,
) -> Result<(Duration, Duration, Duration, u64, Vec<String>), nomics_gpu::GpuError> {
    let info = nomics_gpu::device_info(0)?;
    println!(
        "GPU: {} sm_{}{} ({} SMs, pageable host access: {})",
        info.name,
        info.compute_capability.0,
        info.compute_capability.1,
        info.multiprocessors,
        info.pageable_memory_access
    );
    assert!(
        info.pageable_memory_access,
        "PairHMM experiment requires pageable host memory access"
    );

    let (init, init_log10) = initial_condition();
    let mut engine = GpuPairHmm::new();
    let mut prep = Duration::ZERO;
    let mut kernel = Duration::ZERO;
    let mut post = Duration::ZERO;
    let mut mismatches = 0u64;
    let mut examples = Vec::new();
    let mut begin = 0usize;
    while begin < records.len() {
        let mut end = begin;
        let mut count = 0usize;
        while end < records.len() {
            let record_pairs = records[end].reads.len() * records[end].haps.len();
            if end > begin && count + record_pairs > batch_pairs {
                break;
            }
            count += record_pairs;
            end += 1;
        }

        let pack_start = Instant::now();
        let packed = pack(&records[begin..end], init);
        let mut sums = vec![0.0; packed.pairs.len()];
        prep += pack_start.elapsed();
        let call_start = Instant::now();
        let timing = engine.compute(
            &packed.read_bases,
            &packed.hap_bases,
            &packed.params,
            &packed.pairs,
            &mut sums,
        )?;
        let call_wall = call_start.elapsed();
        prep += call_wall.saturating_sub(timing.kernel + timing.post);
        kernel += timing.kernel;
        post += timing.post;

        let log_start = Instant::now();
        for value in &mut sums {
            *value = nomics_hc::jmath::log10(*value) - init_log10;
        }
        post += log_start.elapsed();

        let mut k = 0usize;
        for (local_record, record) in records[begin..end].iter().enumerate() {
            let record_index = begin + local_record;
            for read_index in 0..record.reads.len() {
                for hap_index in 0..record.haps.len() {
                    let expected =
                        cpu_results[record_index][read_index * record.haps.len() + hap_index];
                    let actual = sums[k];
                    if actual.to_bits() != expected.to_bits() {
                        mismatches += 1;
                        if examples.len() < 8 {
                            examples.push(format!(
                                "record {record_index}, read {read_index}, hap {hap_index}: CPU {expected:.17e} (0x{:016x}), GPU {actual:.17e} (0x{:016x})",
                                expected.to_bits(),
                                actual.to_bits()
                            ));
                        }
                    }
                    k += 1;
                }
            }
        }
        assert_eq!(k, sums.len());
        begin = end;
    }
    Ok((prep, kernel, post, mismatches, examples))
}

/// CPU pool and GPU at the same time, the GPU taking `frac` of the cells
/// (records interleaved). Returns wall seconds and GPU mismatches.
#[cfg(feature = "cuda")]
fn hybrid(
    records: &[Record],
    cpu_results: &[Vec<f64>],
    frac: f64,
    threads: usize,
    batch_pairs: usize,
) -> (f64, f64, f64, u64) {
    let cells = |r: &Record| -> u64 {
        let reads: u64 = r.reads.iter().map(|x| x.bases.len() as u64).sum();
        let haps: u64 = r.haps.iter().map(|h| h.len() as u64).sum();
        reads * haps
    };
    let (mut acc, mut on_gpu) = (0u64, 0u64);
    let (mut gpu_records, mut gpu_expected, mut cpu_records) = (Vec::new(), Vec::new(), Vec::new());
    for (r, expected) in records.iter().zip(cpu_results) {
        let c = cells(r);
        acc += c;
        if (on_gpu as f64) < frac * acc as f64 {
            on_gpu += c;
            gpu_records.push(r.clone());
            gpu_expected.push(expected.clone());
        } else {
            cpu_records.push(r.clone());
        }
    }
    let start = Instant::now();
    let (cpu_secs, gpu_secs, mismatches) = std::thread::scope(|s| {
        let g = s.spawn(|| {
            let t = Instant::now();
            let r = gpu(&gpu_records, &gpu_expected, batch_pairs).expect("run CUDA PairHMM");
            (t.elapsed().as_secs_f64(), r.3)
        });
        let (_, cpu_secs) = cpu(&cpu_records, threads, false);
        let (gpu_secs, mismatches) = g.join().unwrap();
        (cpu_secs, gpu_secs, mismatches)
    });
    (
        start.elapsed().as_secs_f64(),
        cpu_secs,
        gpu_secs,
        mismatches,
    )
}

fn rate(cells: u64, seconds: f64) -> f64 {
    cells as f64 / seconds
}

fn main() {
    let args: Vec<_> = std::env::args_os().collect();
    if args.len() < 2 || args.len() > 4 {
        eprintln!("usage: pairhmm_bench <dump.bin> [threads=20] [gpu_batch_pairs=65536]");
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
    let batch_pairs = args
        .get(3)
        .map(|s| {
            s.to_string_lossy()
                .parse()
                .expect("gpu_batch_pairs must be an integer")
        })
        .unwrap_or(65_536usize);
    assert!(threads > 0);
    assert!(batch_pairs > 0);

    let load_start = Instant::now();
    let records = load(Path::new(&args[1])).expect("load PairHMM dump");
    let load_seconds = load_start.elapsed().as_secs_f64();
    let (pairs, cells, max_read, max_hap) = stats(&records);
    println!(
        "loaded {} records in {load_seconds:.3}s; {pairs} pairs, {cells} cells; max read {max_read}, max hap {max_hap}",
        records.len()
    );

    let (portable_results, portable_seconds) = cpu(&records, threads, true);
    println!(
        "CPU portable lanes ({threads} threads): {portable_seconds:.6}s, {:.3e} cells/s",
        rate(cells, portable_seconds)
    );
    let (trie_results, trie_seconds) = cpu_trie(&records, threads);
    let trie_differing = trie_results
        .iter()
        .flatten()
        .zip(portable_results.iter().flatten())
        .filter(|(a, b)| a.to_bits() != b.to_bits())
        .count();
    println!(
        "CPU prefix trie ({threads} threads): {trie_seconds:.6}s, {:.3e} pair-cells/s; {trie_differing} results differ from the portable lanes (bitwise)",
        rate(cells, trie_seconds)
    );
    let (cpu_results, cpu_seconds) = cpu(&records, threads, false);
    let differing = cpu_results
        .iter()
        .flatten()
        .zip(portable_results.iter().flatten())
        .filter(|(a, b)| a.to_bits() != b.to_bits())
        .count();
    println!(
        "CPU ({threads} threads): {cpu_seconds:.6}s, {:.3e} cells/s; {differing} results differ from the portable lanes (bitwise)",
        rate(cells, cpu_seconds)
    );

    #[cfg(feature = "cuda")]
    {
        let (prep, kernel, post, mismatches, examples) =
            gpu(&records, &cpu_results, batch_pairs).expect("run CUDA PairHMM");
        let total = prep + kernel + post;
        println!(
            "GPU: {:.6}s total (prep/transfer {:.6}s, kernel {:.6}s, post {:.6}s), {:.3e} cells/s total, {:.3e} cells/s kernel",
            total.as_secs_f64(),
            prep.as_secs_f64(),
            kernel.as_secs_f64(),
            post.as_secs_f64(),
            rate(cells, total.as_secs_f64()),
            rate(cells, kernel.as_secs_f64())
        );
        println!("exactness: {mismatches} bit mismatches out of {pairs} pairs");
        for example in examples {
            println!("  {example}");
        }
        // NOMICS_HYBRID="frac:threads,...", e.g. "0.3:20,0.4:19".
        if let Ok(spec) = std::env::var("NOMICS_HYBRID") {
            for item in spec.split(',') {
                let (f, t) = item.split_once(':').unwrap_or((item, "20"));
                let (frac, t): (f64, usize) = (f.parse().unwrap(), t.parse().unwrap());
                let (wall, cpu_s, gpu_s, mism) =
                    hybrid(&records, &cpu_results, frac, t, batch_pairs);
                println!(
                    "hybrid gpu {:.0}% cpu {t} threads: {wall:.3}s wall (cpu part {cpu_s:.3}s, gpu part {gpu_s:.3}s), {:.3e} cells/s, {mism} mismatches",
                    frac * 100.0,
                    rate(cells, wall)
                );
            }
        }
    }

    #[cfg(not(feature = "cuda"))]
    {
        let _ = (cpu_results, batch_pairs);
        println!("GPU: unavailable (build with --features cuda on the Spark)");
    }
}
