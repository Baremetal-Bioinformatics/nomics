//! Benchmark exact CPU/GPU seeding variants and GPU-to-CPU contention.
//!
//! Usage: seed_bench <index> <R1.fastq.gz> <R2.fastq.gz> <table-cache> [pairs] [k] [work-reads]

use std::hint::black_box;
use std::io::BufReader;
use std::path::Path;
use std::sync::{Arc, Barrier};
use std::time::Instant;

use nomics_align::extend::{AlignAux, align1_core};
use nomics_align::index::Index;
use nomics_align::kmer_table::KmerTable;
use nomics_align::mem::read_fastq;
use nomics_align::seed::{SmemAux, collect_intv, collect_intv_with_table};
use nomics_gpu::{GpuBwt, GpuKmerTable, GpuSeedBuffers, SeedWork};
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

fn cpu_seed(idx: &Index, reads: &[nomics_align::sam::Read], table: Option<&KmerTable>) -> f64 {
    let mut aux = SmemAux::default();
    let mut intervals = 0usize;
    let started = Instant::now();
    for read in reads {
        match table {
            Some(table) => collect_intv_with_table(&idx.bwt, &read.seq, &mut aux, table),
            None => collect_intv(&idx.bwt, &read.seq, &mut aux),
        }
        intervals = intervals.wrapping_add(aux.mem.len());
    }
    black_box(intervals);
    started.elapsed().as_secs_f64()
}

fn gpu_seed(
    bwt: &GpuBwt<'_>,
    table: Option<&GpuKmerTable<'_>>,
    packed: &PackedReads,
    cap: u32,
    work: SeedWork,
) -> Result<(f32, usize), nomics_gpu::GpuError> {
    let mut output = GpuSeedBuffers::new(packed.n, cap)?;
    let ms = output.collect(
        bwt,
        table,
        packed.seqs.as_slice(),
        packed.offs(),
        packed.lens(),
        work,
    )?;
    Ok((ms, output.overflow_count()))
}

fn alignment_workload(idx: &Index, reads: &[nomics_align::sam::Read], count: usize) -> f64 {
    let mut aux = AlignAux::default();
    let mut regions = 0usize;
    let started = Instant::now();
    for read in reads.iter().take(count) {
        regions = regions.wrapping_add(align1_core(idx, &read.seq, &mut aux).len());
    }
    black_box(regions);
    started.elapsed().as_secs_f64()
}

fn contended_alignment(
    idx: &Index,
    reads: &[nomics_align::sam::Read],
    work_reads: usize,
    bwt: GpuBwt<'_>,
    table: Option<GpuKmerTable<'_>>,
    packed: &PackedReads,
) -> Result<(f64, f32, usize), nomics_gpu::GpuError> {
    let mut output = GpuSeedBuffers::new(packed.n, 128)?;
    let barrier = Arc::new(Barrier::new(2));
    std::thread::scope(|scope| {
        let gpu_barrier = Arc::clone(&barrier);
        let gpu = scope.spawn(move || {
            gpu_barrier.wait();
            let ms = output.collect(
                &bwt,
                table.as_ref(),
                packed.seqs.as_slice(),
                packed.offs(),
                packed.lens(),
                SeedWork::TwoPass64Then128,
            )?;
            Ok::<_, nomics_gpu::GpuError>((ms, output.overflow_count()))
        });
        barrier.wait();
        let cpu_seconds = alignment_workload(idx, reads, work_reads);
        let (gpu_ms, overflows) = gpu.join().unwrap()?;
        Ok((cpu_seconds, gpu_ms, overflows))
    })
}

/// FNV-1a over every read's count and intervals, in original read order
/// (`pos[r]` is where read `r` was packed): a quick identity check between
/// kernel variants (the exactness gate is `seed_gate`).
fn digest(output: &GpuSeedBuffers, pos: &[usize]) -> u64 {
    let mut h = 0xcbf29ce484222325u64;
    let mut eat = |v: u64| {
        for b in v.to_le_bytes() {
            h = (h ^ b as u64).wrapping_mul(0x100000001b3);
        }
    };
    for &p in pos {
        eat(output.count(p) as u64);
        for interval in output.get(p).unwrap_or(&[]) {
            interval.x.iter().for_each(|&x| eat(x));
            eat(interval.info);
        }
    }
    h
}

/// Bytes of `[ptr, ptr + len)` on transparent huge pages, from smaps.
fn huge_bytes(ptr: usize, len: usize) -> Option<usize> {
    let smaps = std::fs::read_to_string("/proc/self/smaps").ok()?;
    let mut total = 0;
    let mut in_range = false;
    for line in smaps.lines() {
        if let Some((range, _)) = line.split_once(' ')
            && let Some((lo, hi)) = range.split_once('-')
            && let (Ok(lo), Ok(hi)) = (usize::from_str_radix(lo, 16), usize::from_str_radix(hi, 16))
        {
            in_range = lo < ptr + len && hi > ptr;
            continue;
        }
        if in_range && let Some(rest) = line.strip_prefix("AnonHugePages:") {
            total += rest
                .trim()
                .trim_end_matches("kB")
                .trim()
                .parse::<usize>()
                .ok()?
                * 1024;
        }
    }
    Some(total)
}

/// `SEED_BENCH_GPU_REPS=n`: only the pipeline's configuration (cap 64, no
/// table), `n` timed reps, plus one two-pass run; prints a digest of the
/// output. The index is collapsed to 2 MiB pages first, as the pipeline does.
/// `SEED_BENCH_VEC=1` passes fresh `Vec` buffers each rep like
/// `mem::gpu_seeds`. `SEED_BENCH_ORDER=kmer|reverse` packs the reads in that
/// order (outputs are compared in the original order). `SEED_BENCH_TABLE=k`
/// passes the k-mer table (the `<table-cache>` argument, built if missing).
fn gpu_only(idx: &Index, reads: &[nomics_align::sam::Read], reps: usize, table_path: &Path) {
    if let Err(e) = idx.lock_memory() {
        println!("lock_memory: {e}");
    }
    let (total, huge) = idx.memory();
    let words = idx.bwt.words();
    let bwt_huge = huge_bytes(words.as_ptr() as usize, std::mem::size_of_val(words));
    println!(
        "index {:.1} GB, {:.0}% on 2 MiB pages; BWT {:.2} GB, {:.1}% on 2 MiB pages",
        total as f64 / 1e9,
        100.0 * huge.unwrap_or(0) as f64 / total as f64,
        std::mem::size_of_val(words) as f64 / 1e9,
        100.0 * bwt_huge.unwrap_or(0) as f64 / std::mem::size_of_val(words) as f64
    );
    let mut order: Vec<usize> = (0..reads.len()).collect();
    match std::env::var("SEED_BENCH_ORDER").as_deref() {
        Ok("reverse") => order.reverse(),
        Ok("kmer") => order.sort_by_key(|&r| {
            reads[r]
                .seq
                .iter()
                .take(12)
                .fold(0u32, |a, &c| a << 2 | (c.min(3) as u32))
        }),
        _ => {}
    }
    let mut pos = vec![0; reads.len()];
    for (p, &r) in order.iter().enumerate() {
        pos[r] = p;
    }
    let ordered: Vec<_> = order.iter().map(|&r| reads[r].clone()).collect();
    let packed = pack(&ordered);
    let _ = packed.seqs.collapse_huge();
    let bwt = GpuBwt {
        primary: idx.bwt.primary,
        l2: idx.bwt.l2,
        seq_len: idx.bwt.seq_len,
        words,
    };
    let table = std::env::var("SEED_BENCH_TABLE").ok().map(|k| {
        let table = KmerTable::load_or_build(table_path, &idx.bwt, k.parse().unwrap()).unwrap();
        println!("k-mer table k={k}, {:.3} GiB", table.bytes() as f64 / (1u64 << 30) as f64);
        table
    });
    let gpu_table = table.as_ref().map(|t| t.gpu_table());
    if std::env::var_os("SEED_BENCH_VEC").is_some() {
        let mut times = Vec::new();
        for _ in 0..reps {
            let seqs: Vec<u8> = packed.seqs.as_slice().to_vec();
            let mut out = vec![nomics_gpu::GpuIntv::default(); reads.len() * 64];
            let mut n = vec![0u32; reads.len()];
            let ms = nomics_gpu::seed_collect(
                &bwt,
                &seqs,
                packed.offs(),
                packed.lens(),
                64,
                &mut out,
                &mut n,
            )
            .unwrap();
            times.push(format!("{ms:.1}"));
        }
        println!("GPU cap64 (fresh Vec buffers): ms {times:?}");
        return;
    }
    for (cap, work, label) in [
        (64, SeedWork::Cap64, "cap64"),
        (128, SeedWork::TwoPass64Then128, "two-pass"),
    ] {
        if cap == 128 && std::env::var_os("SEED_BENCH_NO_TWOPASS").is_some() {
            continue;
        }
        let mut output = GpuSeedBuffers::new(reads.len(), cap).unwrap();
        output.touch_and_collapse();
        let n = if cap == 64 { reps } else { 1 };
        let mut times = Vec::new();
        for _ in 0..n {
            let ms = output
                .collect(
                    &bwt,
                    gpu_table.as_ref(),
                    packed.seqs.as_slice(),
                    packed.offs(),
                    packed.lens(),
                    work,
                )
                .unwrap();
            times.push(ms);
        }
        let intervals: u64 = (0..reads.len())
            .filter_map(|r| output.get(r))
            .map(|s| s.len() as u64)
            .sum();
        println!(
            "GPU {label}: {} reads, ms {:?}, {} overflow, {intervals} intervals, digest {:016x}",
            reads.len(),
            times.iter().map(|t| format!("{t:.1}")).collect::<Vec<_>>(),
            output.overflow_count(),
            digest(&output, &pos)
        );
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    assert!(
        (5..=8).contains(&args.len()),
        "usage: seed_bench <index> <R1.fastq.gz> <R2.fastq.gz> <table-cache> [pairs] [k] [work-reads]"
    );
    let pairs = args.get(5).map_or(100_000, |value| value.parse().unwrap());
    let k = args.get(6).map_or(12, |value| value.parse().unwrap());
    let work_reads = args.get(7).map_or(20_000, |value| value.parse().unwrap());
    let idx = Index::load(Path::new(&args[1])).unwrap();
    let reads = read_pairs(Path::new(&args[2]), Path::new(&args[3]), pairs);
    if let Some(reps) = std::env::var_os("SEED_BENCH_GPU_REPS") {
        gpu_only(
            &idx,
            &reads,
            reps.to_str().unwrap().parse().unwrap(),
            Path::new(&args[4]),
        );
        return;
    }

    let ready = Instant::now();
    let table = KmerTable::load_or_build(Path::new(&args[4]), &idx.bwt, k).unwrap();
    println!(
        "{} reads, k={k}, table {:.3} GiB (ready {:.3}s)",
        reads.len(),
        table.bytes() as f64 / (1u64 << 30) as f64,
        ready.elapsed().as_secs_f64()
    );
    for (label, selected) in [("baseline", None), ("table", Some(&table))] {
        let seconds = cpu_seed(&idx, &reads, selected);
        println!(
            "CPU {label}: {:.3}s, {:.3} us/read",
            seconds,
            seconds * 1e6 / reads.len() as f64
        );
    }

    let packed = pack(&reads);
    let bwt = GpuBwt {
        primary: idx.bwt.primary,
        l2: idx.bwt.l2,
        seq_len: idx.bwt.seq_len,
        words: idx.bwt.words(),
    };
    let gpu_table = table.gpu_table();
    let variants = [
        (64, SeedWork::Cap64, "cap64"),
        (96, SeedWork::Cap96, "cap96"),
        (128, SeedWork::Cap128, "cap128"),
        (128, SeedWork::TwoPass64Then128, "two-pass-64-128"),
    ];
    let mut gpu_available = true;
    for (table_label, selected) in [("baseline", None), ("table", Some(&gpu_table))] {
        for (cap, work, cap_label) in variants {
            match gpu_seed(&bwt, selected, &packed, cap, work) {
                Ok((ms, overflows)) => println!(
                    "GPU {table_label} {cap_label}: {ms:.3} ms/batch, {overflows} overflow reads"
                ),
                Err(nomics_gpu::GpuError::Unavailable) => {
                    gpu_available = false;
                    break;
                }
                Err(error) => panic!("GPU benchmark failed: {error}"),
            }
        }
        if !gpu_available {
            break;
        }
    }
    if !gpu_available {
        println!("GPU unavailable; rerun with --features cuda on the Spark");
        return;
    }

    let work_reads = work_reads.min(reads.len());
    let alone = alignment_workload(&idx, &reads, work_reads);
    println!(
        "CPU align alone: {:.3}s for {work_reads} reads ({:.3} us/read)",
        alone,
        alone * 1e6 / work_reads as f64
    );
    for (label, selected) in [("GPU baseline", None), ("GPU table", Some(gpu_table))] {
        let (cpu_seconds, gpu_ms, overflows) =
            contended_alignment(&idx, &reads, work_reads, bwt, selected, &packed).unwrap();
        println!(
            "CPU align + {label}: {:.3}s ({:.3}x alone); GPU {gpu_ms:.3} ms, {overflows} overflows",
            cpu_seconds,
            cpu_seconds / alone
        );
    }
}
