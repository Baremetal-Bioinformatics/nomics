//! Exactness gate for the k-mer table and configurable GPU seed kernel.
//!
//! Usage: seed_gate <index> <R1.fastq.gz> <R2.fastq.gz> <table-cache> [pairs] [k]

use std::io::BufReader;
use std::path::Path;
use std::time::Instant;

use nomics_align::index::Index;
use nomics_align::kmer_table::KmerTable;
use nomics_align::mem::read_fastq;
use nomics_align::seed::{Intv, SmemAux, collect_intv, collect_intv_with_table};
use nomics_gpu::{GpuBwt, GpuSeedBuffers, SeedWork};
use nomics_io::gzip::GzReader;
use nomics_mem::HostBuf;

struct PackedReads {
    seqs: HostBuf,
    // HostBufs like every GPU-read buffer (docs/METHODOLOGY.md); heap pages here once
    // coincided with wrong seeds in the pipeline.
    offs: HostBuf,
    lens: HostBuf,
    n: usize,
}

impl PackedReads {
    fn offs(&self) -> &[u64] {
        &self.offs.as_u64s()[..self.n]
    }
    fn lens(&self) -> &[u32] {
        &self.lens.as_u32s()[..self.n]
    }
}

fn read_pairs(r1_path: &Path, r2_path: &Path, max_pairs: usize) -> Vec<nomics_align::sam::Read> {
    let mut r1 = BufReader::new(GzReader::open(r1_path).unwrap());
    let mut r2 = BufReader::new(GzReader::open(r2_path).unwrap());
    let mut line = Vec::new();
    let mut reads = Vec::new();
    for _ in 0..max_pairs {
        match (
            read_fastq(&mut r1, &mut line).unwrap(),
            read_fastq(&mut r2, &mut line).unwrap(),
        ) {
            (Some(a), Some(b)) => reads.extend([a, b]),
            (None, None) => break,
            _ => panic!("R1 and R2 differ in read count"),
        }
    }
    reads
}

fn pack(reads: &[nomics_align::sam::Read]) -> PackedReads {
    let bytes: usize = reads.iter().map(|read| read.seq.len()).sum();
    let mut seqs = HostBuf::new(bytes).unwrap();
    let mut offs = HostBuf::new(reads.len().max(1) * 8).unwrap();
    let mut lens = HostBuf::new(reads.len().max(1) * 4).unwrap();
    let mut at = 0usize;
    for (i, read) in reads.iter().enumerate() {
        offs.as_u64s_mut()[i] = at as u64;
        lens.as_u32s_mut()[i] = read.seq.len() as u32;
        seqs.as_mut_slice()[at..at + read.seq.len()].copy_from_slice(&read.seq);
        at += read.seq.len();
    }
    PackedReads { seqs, offs, lens, n: reads.len() }
}

fn as_intvs(gpu: &[nomics_gpu::GpuIntv]) -> Vec<Intv> {
    gpu.iter()
        .map(|interval| Intv {
            x: interval.x,
            info: interval.info,
        })
        .collect()
}

/// GPU intervals vs the CPU reference for every read the GPU finished
/// (threaded). Returns the overflow read indices (those go to the CPU in the
/// pipeline); `SEED_GATE_OVERFLOW_OUT=<prefix>` also writes them to
/// `<prefix>.<label>.txt` so the set can be compared between kernels.
fn check_gpu(
    label: &str,
    idx: &Index,
    table: Option<&nomics_gpu::GpuKmerTable<'_>>,
    reads: &[nomics_align::sam::Read],
    packed: &PackedReads,
    cap: u32,
    work: SeedWork,
) -> Result<Vec<usize>, nomics_gpu::GpuError> {
    let bwt = GpuBwt {
        primary: idx.bwt.primary,
        l2: idx.bwt.l2,
        seq_len: idx.bwt.seq_len,
        words: idx.bwt.words(),
    };
    let mut output = GpuSeedBuffers::new(reads.len(), cap)?;
    output.touch_and_collapse();
    let ms = output.collect(
        &bwt,
        table,
        packed.seqs.as_slice(),
        packed.offs(),
        packed.lens(),
        work,
    )?;
    let threads = std::thread::available_parallelism().map_or(8, |n| n.get());
    let chunk = reads.len().div_ceil(threads).max(1);
    let output = &output;
    let checked: usize = std::thread::scope(|scope| {
        let handles: Vec<_> = reads
            .chunks(chunk)
            .enumerate()
            .map(|(c, part)| {
                scope.spawn(move || {
                    let mut aux = SmemAux::default();
                    let mut checked = 0usize;
                    for (o, read) in part.iter().enumerate() {
                        let read_index = c * chunk + o;
                        let Some(gpu) = output.get(read_index) else {
                            continue;
                        };
                        collect_intv(&idx.bwt, &read.seq, &mut aux);
                        assert_eq!(
                            aux.mem,
                            as_intvs(gpu),
                            "{label}: GPU differs at read {read_index} ({})",
                            String::from_utf8_lossy(&read.name)
                        );
                        checked += 1;
                    }
                    checked
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join().unwrap()).sum()
    });
    let overflow: Vec<usize> = (0..reads.len())
        .filter(|&r| output.get(r).is_none())
        .collect();
    let with_n = overflow
        .iter()
        .filter(|&&r| reads[r].seq.iter().any(|&b| b > 3))
        .count();
    println!(
        "{label}: GPU {ms:.3} ms, {} reads, {checked} identical to the CPU, {} overflow to CPU ({with_n} of them with N)",
        reads.len(),
        overflow.len()
    );
    Ok(overflow)
}

/// Gate `reads` in batches of `SEED_GATE_BATCH` reads (default 662,000, a
/// 30x pipeline batch) so the result buffers stay batch-sized.
fn check_all(
    label: &str,
    idx: &Index,
    table: Option<&nomics_gpu::GpuKmerTable<'_>>,
    reads: &[nomics_align::sam::Read],
    cap: u32,
    work: SeedWork,
) -> Result<(), nomics_gpu::GpuError> {
    let batch = std::env::var("SEED_GATE_BATCH").map_or(662_000, |v| v.parse().unwrap());
    let mut overflow = Vec::new();
    for (b, chunk) in reads.chunks(batch).enumerate() {
        let packed = pack(chunk);
        let local = check_gpu(label, idx, table, chunk, &packed, cap, work)?;
        overflow.extend(local.into_iter().map(|r| b * batch + r));
    }
    println!(
        "{label} total: {} reads, {} overflow to CPU, every other read identical to the CPU",
        reads.len(),
        overflow.len()
    );
    if let Some(prefix) = std::env::var_os("SEED_GATE_OVERFLOW_OUT") {
        let mut path = prefix;
        path.push(format!(".{}.txt", label.replace(' ', "_")));
        let text: String = overflow.iter().map(|r| format!("{r}\n")).collect();
        std::fs::write(path, text).unwrap();
    }
    Ok(())
}

/// Gates the pipeline's first pass (cap 64, no table), fused two-pass at
/// stride 64 and 128, and cap 128 on every read (the pipeline's two-pass
/// retries cap-64 overflow reads with cap 128; `NOMICS_SEED_VERIFY=1` gates
/// that in a full run);
/// `SEED_GATE_TABLE=1` also gates the CPU and GPU k-mer table paths (loads or
/// builds the table).
fn main() {
    let args: Vec<String> = std::env::args().collect();
    assert!(
        (5..=7).contains(&args.len()),
        "usage: seed_gate <index> <R1.fastq.gz> <R2.fastq.gz> <table-cache> [pairs] [k]"
    );
    let pairs = args.get(5).map_or(100_000, |value| value.parse().unwrap());
    let k = args.get(6).map_or(12, |value| value.parse().unwrap());
    let idx = Index::load(Path::new(&args[1])).unwrap();
    if let Err(e) = idx.lock_memory() {
        println!("lock_memory: {e}");
    }
    let reads = read_pairs(Path::new(&args[2]), Path::new(&args[3]), pairs);
    let with_n = reads
        .iter()
        .filter(|r| r.seq.iter().any(|&b| b > 3))
        .count();
    let long = reads.iter().filter(|r| r.seq.len() > 256).count();
    println!(
        "{} reads, {with_n} with N, {long} longer than 256",
        reads.len()
    );
    for (label, cap, work) in [
        ("cap64", 64, SeedWork::Cap64),
        ("two-pass-stride64", 64, SeedWork::TwoPass64Then128),
        ("cap128-stride128", 128, SeedWork::Cap128),
        ("two-pass", 128, SeedWork::TwoPass64Then128),
    ] {
        match check_all(label, &idx, None, &reads, cap, work) {
            Ok(_) => {}
            Err(nomics_gpu::GpuError::Unavailable) => {
                println!("GPU unavailable; rerun with --features cuda on the Spark");
                return;
            }
            Err(error) => panic!("GPU gate failed: {error}"),
        }
    }
    if std::env::var_os("SEED_GATE_TABLE").is_none() {
        return;
    }
    let started = Instant::now();
    let table = KmerTable::load_or_build(Path::new(&args[4]), &idx.bwt, k).unwrap();
    println!(
        "k={k} table: {:.3} GiB, ready in {:.3}s",
        table.bytes() as f64 / (1u64 << 30) as f64,
        started.elapsed().as_secs_f64()
    );
    let (mut baseline, mut accelerated) = (SmemAux::default(), SmemAux::default());
    for (read_index, read) in reads.iter().enumerate() {
        collect_intv(&idx.bwt, &read.seq, &mut baseline);
        collect_intv_with_table(&idx.bwt, &read.seq, &mut accelerated, &table);
        assert_eq!(
            baseline.mem,
            accelerated.mem,
            "table differs at read {read_index} ({})",
            String::from_utf8_lossy(&read.name)
        );
    }
    println!("CPU baseline == CPU table for {} reads", reads.len());
    let gpu_table = table.gpu_table();
    check_all(
        "table two-pass",
        &idx,
        Some(&gpu_table),
        &reads,
        128,
        SeedWork::TwoPass64Then128,
    )
    .unwrap();
}
