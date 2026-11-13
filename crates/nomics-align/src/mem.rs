//! `bwa mem -K 100000000 -Y` driver for paired FASTQ (`fastmap.c` `process`,
//! `bwa.c` `bseq_read`, `bwamem.c` `mem_process_seqs`).
//!
//! Gate: the SAM records bwa writes to stdout, before sorting.

use std::io::{self, BufRead};
use std::sync::atomic::{AtomicUsize, Ordering};

use crate::extend::{AlignAux, Region, align1_core, align1_core_from};
use crate::index::Index;
use crate::pair::{PeStat, RescueKey, RescuePlan, Rescued, bam_pe, pestat, plan_rescue, sam_pe};
use crate::sam::{AlnSink, BamRecords, ChunkRecords, Read};
use crate::samsort::sort_chunk;
use crate::seed::Intv;
use crate::seed::nt4;
use nomics_io::bam::{Header, RecordSet};

/// `-K 100000000`: bases per batch (R1 + R2).
pub const CHUNK_SIZE: i64 = 100_000_000;

/// One FASTQ record as bwa's `kseq` reads it: name up to the first
/// whitespace, minus a trailing `/[0-9]` (`trim_readno`); comment dropped.
pub fn read_fastq(r: &mut impl BufRead, line: &mut Vec<u8>) -> io::Result<Option<Read>> {
    line.clear();
    if r.read_until(b'\n', line)? == 0 {
        return Ok(None);
    }
    if line.first() != Some(&b'@') {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "FASTQ header without '@'",
        ));
    }
    let end = line[1..]
        .iter()
        .position(|c| c.is_ascii_whitespace())
        .map_or(line.len(), |p| p + 1);
    let mut name = line[1..end].to_vec();
    let l = name.len();
    if l > 2 && name[l - 2] == b'/' && name[l - 1].is_ascii_digit() {
        name.truncate(l - 2);
    }
    let mut take = |r: &mut dyn BufRead| -> io::Result<Vec<u8>> {
        line.clear();
        r.read_until(b'\n', line)?;
        while matches!(line.last(), Some(b'\n' | b'\r')) {
            line.pop();
        }
        Ok(line.clone())
    };
    let seq: Vec<u8> = take(r)?.iter().map(|&b| nt4(b)).collect();
    take(r)?; // '+'
    let qual = take(r)?;
    if qual.len() != seq.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "FASTQ quality length differs",
        ));
    }
    Ok(Some(Read { name, seq, qual }))
}

/// One `-K` batch of pairs, reads interleaved R1, R2.
pub fn read_batch(r1: &mut impl BufRead, r2: &mut impl BufRead) -> io::Result<Vec<Read>> {
    read_batch_limited(r1, r2, usize::MAX)
}

/// As [`read_batch`], stopping after at most `max_pairs` pairs.
pub fn read_batch_limited(
    r1: &mut impl BufRead,
    r2: &mut impl BufRead,
    max_pairs: usize,
) -> io::Result<Vec<Read>> {
    let mut reads = Vec::new();
    let mut size = 0i64;
    let mut line = Vec::new();
    loop {
        let a = read_fastq(r1, &mut line)?;
        let b = read_fastq(r2, &mut line)?;
        match (a, b) {
            (Some(a), Some(b)) => {
                size += (a.seq.len() + b.seq.len()) as i64;
                reads.push(a);
                reads.push(b);
                if size >= CHUNK_SIZE || reads.len() / 2 >= max_pairs {
                    break;
                }
            }
            (None, None) => break,
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "R1 and R2 differ in read count",
                ));
            }
        }
    }
    Ok(reads)
}

/// Pairs per unit of parallel work.
const BLOCK_PAIRS: usize = 256;

/// Run `f` over work items with at most `threads` workers, on the
/// process-wide host pool ([`nomics_gpu::pool`]: persistent threads, not
/// threads spawned per call). Items are independent, so the result does not
/// depend on scheduling. Returns the summed busy time of the workers.
pub(crate) fn for_each_parallel<T: Send>(
    threads: usize,
    items: Vec<T>,
    f: impl Fn(T) + Sync,
) -> f64 {
    nomics_gpu::pool::for_each(threads, items, f)
}

trait BatchOutput: AlnSink + Default + Send {
    /// A finished batch's output as a chunk takes it ([`Self::finish`] runs
    /// on the batch's own thread).
    type Batch: Send;
    /// A chunk's batches, accumulated in batch order.
    type Chunk: Default + Send;
    fn append(&mut self, other: Self);
    fn finish(self) -> Self::Batch;
    fn absorb(chunk: &mut Self::Chunk, batch: Self::Batch);
    /// A chunk's bytes held: (in use, allocated). Diagnostics only.
    fn held(chunk: &Self::Chunk) -> (usize, usize);
}

impl BatchOutput for Vec<u8> {
    type Batch = Vec<u8>;
    type Chunk = Vec<u8>;
    fn append(&mut self, mut other: Self) {
        Vec::append(self, &mut other);
    }
    fn finish(self) -> Self::Batch {
        self
    }
    fn absorb(chunk: &mut Self::Chunk, mut batch: Self::Batch) {
        chunk.append(&mut batch);
    }
    fn held(chunk: &Self::Chunk) -> (usize, usize) {
        (chunk.len(), chunk.capacity())
    }
}

impl BatchOutput for BamRecords {
    type Batch = crate::sam::BinnedRecords;
    type Chunk = ChunkRecords;
    fn append(&mut self, other: Self) {
        BamRecords::append(self, other);
    }
    fn finish(self) -> Self::Batch {
        self.into_bins()
    }
    fn absorb(chunk: &mut Self::Chunk, batch: Self::Batch) {
        chunk.absorb(batch);
    }
    fn held(chunk: &Self::Chunk) -> (usize, usize) {
        chunk.held()
    }
}

/// The time origin of `[timing]` and `[mem]` lines: the first batch's start.
fn timing_epoch() -> std::time::Instant {
    static EPOCH: std::sync::OnceLock<std::time::Instant> = std::sync::OnceLock::new();
    *EPOCH.get_or_init(std::time::Instant::now)
}

/// Seed buffer sets currently allocated (live or pooled). Diagnostics only.
static SEED_SETS: AtomicUsize = AtomicUsize::new(0);

/// One `[mem]` line under `NOMICS_TIMING`: anonymous RSS, HostBuf bytes, seed
/// sets (allocated / pooled), the engines idle in the rescue and fused pools
/// (count and bytes; engines in use are not counted), and what the caller
/// adds.
pub(crate) fn mem_line(what: &str) {
    let gb = |b: u64| b as f64 / 1e9;
    let pooled = SEED_POOL.lock().map_or(0, |p| p.len());
    let (rescue_n, rescue_bytes) = RESCUE_ENGINES
        .lock()
        .map_or((0, 0), |p| (p.len(), p.iter().map(|e| e.bytes()).sum::<usize>()));
    let (fused_n, fused_bytes) = FUSED_ENGINES
        .lock()
        .map_or((0, 0), |p| (p.len(), p.iter().map(|e| e.bytes()).sum::<usize>()));
    eprintln!(
        "[mem] {:.3} {what}: rss_anon {:.2} GB, hostbuf {:.2} GB, seed sets {} ({pooled} pooled), idle rescue engines {rescue_n} {:.2} GB, idle fused engines {fused_n} {:.2} GB",
        timing_epoch().elapsed().as_secs_f64(),
        gb(nomics_mem::rss_anon_bytes().unwrap_or(0)),
        gb(nomics_mem::hostbuf_bytes() as u64),
        SEED_SETS.load(Ordering::Relaxed),
        gb(rescue_bytes as u64),
        gb(fused_bytes as u64),
    );
}

/// Align one batch whose first pair has global index `first_pair` on
/// `threads` threads; appends SAM records to `out` (in input order) and
/// returns the batch's insert-size stats and log. Like bwa, every read is
/// aligned independently (`worker1`), insert sizes are estimated from the
/// whole batch, then every pair is resolved independently (`worker2`).
pub fn process_batch(
    idx: &Index,
    rg_id: &[u8],
    reads: &[Read],
    first_pair: u64,
    threads: usize,
    gpu: bool,
    out: &mut Vec<u8>,
) -> ([PeStat; 4], Vec<String>) {
    let seeds = gpu.then(|| try_gpu_seeds(idx, reads)).flatten();
    process_batch_seeded(idx, rg_id, reads, seeds, first_pair, threads, out)
}

/// [`process_batch`] with seed intervals already computed (e.g. on the GPU,
/// overlapped with the previous batch); `None` seeds on the CPU.
pub fn process_batch_seeded(
    idx: &Index,
    rg_id: &[u8],
    reads: &[Read],
    seeds: Option<GpuSeeds>,
    first_pair: u64,
    threads: usize,
    out: &mut Vec<u8>,
) -> ([PeStat; 4], Vec<String>) {
    process_batch_seeded_to(idx, rg_id, reads, seeds, first_pair, threads, out, sam_pe)
}

pub fn process_batch_bam(
    idx: &Index,
    rg_id: &[u8],
    reads: &[Read],
    seeds: Option<GpuSeeds>,
    first_pair: u64,
    threads: usize,
    out: &mut BamRecords,
) -> ([PeStat; 4], Vec<String>) {
    process_batch_seeded_to(idx, rg_id, reads, seeds, first_pair, threads, out, bam_pe)
}

#[allow(clippy::too_many_arguments)]
fn process_batch_seeded_to<O: BatchOutput>(
    idx: &Index,
    rg_id: &[u8],
    reads: &[Read],
    seeds: Option<GpuSeeds>,
    first_pair: u64,
    threads: usize,
    out: &mut O,
    pair_fn: PairFn<O>,
) -> ([PeStat; 4], Vec<String>) {
    let timing = std::env::var_os("NOMICS_TIMING").is_some();
    let t0 = std::time::Instant::now();
    let epoch = timing_epoch();
    let started = epoch.elapsed().as_secs_f64();
    let _diag_batch = crate::diag::batch(reads, first_pair);
    let mut regs: Vec<Vec<Region>> = vec![Vec::new(); reads.len()];
    let t_seed = t0.elapsed().as_secs_f64();
    let work: Vec<_> = regs
        .chunks_mut(2 * BLOCK_PAIRS)
        .zip(reads.chunks(2 * BLOCK_PAIRS))
        .enumerate()
        .collect();
    let mode = gpu_extend_mode();
    // Only align1 reads the seeds; rescue and pairing just need to know
    // whether the batch was GPU-seeded.
    let gpu_seeded = seeds.is_some();
    let busy1 = {
        let seeds = seeds.as_ref();
        if mode == ExtendMode::Fused || (seeds.is_some() && mode == ExtendMode::TwoPass) {
            align1_gpu_extend(idx, reads, seeds, threads, work, mode)
        } else {
            for_each_parallel(threads, work, |(b, (slots, rs))| {
                let mut aux = AlignAux::default();
                for (j, (slot, r)) in slots.iter_mut().zip(rs).enumerate() {
                    let i = b * 2 * BLOCK_PAIRS + j;
                    crate::diag::at_read(r);
                    *slot = match seeds.and_then(|s| s.get(i)) {
                        Some(mem) => align1_core_from(idx, &r.seq, mem, &mut aux),
                        None => align1_core(idx, &r.seq, &mut aux),
                    };
                }
            })
        }
    };
    drop(seeds);

    let t1 = t0.elapsed().as_secs_f64();
    if timing {
        let vec = std::mem::size_of::<Vec<u8>>();
        let reads_bytes: usize = reads
            .iter()
            .map(|r| 3 * vec + r.name.capacity() + r.seq.capacity() + r.qual.capacity())
            .sum();
        let regs_bytes: usize = regs
            .iter()
            .map(|r| vec + r.capacity() * std::mem::size_of::<Region>())
            .sum();
        mem_line(&format!(
            "batch {first_pair} aligned, reads {:.2} GB, regions {:.2} GB",
            reads_bytes as f64 / 1e9,
            regs_bytes as f64 / 1e9
        ));
    }
    crate::diag::check_batch(
        "after align1, before pestat and GPU mate rescue",
        idx,
        reads,
        &regs,
        first_pair,
    );
    let (pes, log) = pestat(idx.bns.l_pac, &regs);
    let t2 = t0.elapsed().as_secs_f64();
    let rescue = (gpu_seeded && gpu_rescue_enabled())
        .then(|| gpu_rescue(idx, &pes, reads, &regs, threads))
        .flatten();
    crate::diag::check_batch(
        "after the GPU mate-rescue kernel, before pairing",
        idx,
        reads,
        &regs,
        first_pair,
    );
    let t_rescue = t0.elapsed().as_secs_f64() - t2;

    let mut outs: Vec<O> = (0..regs.len().div_ceil(2 * BLOCK_PAIRS))
        .map(|_| O::default())
        .collect();
    let work: Vec<_> = regs
        .chunks_mut(2 * BLOCK_PAIRS)
        .zip(reads.chunks(2 * BLOCK_PAIRS))
        .zip(outs.iter_mut())
        .enumerate()
        .collect();
    let busy2 = for_each_parallel(threads, work, |(b, ((slots, rs), o))| {
        for (i, (pair, s)) in slots
            .chunks_exact_mut(2)
            .zip(rs.chunks_exact(2))
            .enumerate()
        {
            let s = [s[0].clone(), s[1].clone()];
            let mut a = [std::mem::take(&mut pair[0]), std::mem::take(&mut pair[1])];
            let p = b * BLOCK_PAIRS + i;
            let id = first_pair + p as u64;
            let rescued = rescue.as_ref().map_or(Rescued::default(), |r| r.pair(p));
            pair_fn(idx, rg_id, &pes, id, &s, &mut a, &rescued, o);
        }
    });
    let t3 = t0.elapsed().as_secs_f64();
    for o in outs {
        out.append(o);
    }
    if timing {
        let ended = epoch.elapsed().as_secs_f64();
        eprintln!(
            "[timing] batch {first_pair}..{} [{started:.3},{ended:.3}] workers {threads}: {} reads, gpu seed {t_seed:.2}s, align {:.2}s (busy {busy1:.1}s), pestat {:.3}s, gpu rescue {t_rescue:.2}s, pair {:.2}s (busy {busy2:.1}s), concat {:.3}s",
            first_pair + (reads.len() / 2) as u64,
            reads.len(),
            t1 - t_seed,
            t2 - t1,
            t3 - t2 - t_rescue,
            t0.elapsed().as_secs_f64() - t3
        );
    }
    (pes, log)
}

/// The pairing step: `sam_pe` or `bam_pe`.
type PairFn<O> =
    fn(&Index, &[u8], &[PeStat; 4], u64, &[Read; 2], &mut [Vec<Region>; 2], &Rescued, &mut O);

/// Mate-rescue alignments run on the GPU whenever seeding does;
/// `NOMICS_GPU_RESCUE=0` keeps them on the CPU (identical output).
fn gpu_rescue_enabled() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| std::env::var("NOMICS_GPU_RESCUE").map_or(true, |v| v != "0"))
}

/// Longest rescue window sent to the GPU (the kernel stages it in shared
/// memory); longer ones are aligned on the CPU.
const GPU_RESCUE_MAX_TARGET: usize = 4000;

/// Engines kept for reuse across batches (and lanes): their huge-page
/// buffers stay on pages the GPU has translated, and each keeps its CUDA
/// stream. At most this many of each kind are kept.
// 12 covers 3 batches x up to 4 lanes; idle warmed engines hold only a stream.
const ENGINES_KEEP: usize = 12;
static RESCUE_ENGINES: std::sync::Mutex<Vec<nomics_gpu::GpuKswAlign2>> =
    std::sync::Mutex::new(Vec::new());
static FUSED_ENGINES: std::sync::Mutex<Vec<crate::fused::FusedEngine>> =
    std::sync::Mutex::new(Vec::new());
static TWO_PASS_ENGINES: std::sync::Mutex<Vec<nomics_gpu::GpuKswExtend2>> =
    std::sync::Mutex::new(Vec::new());

/// GPU calls that failed and ran on the CPU instead (same output, slower):
/// seed batches, extension tiles, rescue batches.
static GPU_FALLBACKS: [AtomicUsize; 3] = [AtomicUsize::new(0), AtomicUsize::new(0), AtomicUsize::new(0)];

/// One line for the end of a run when any GPU call fell back to the CPU, so
/// a slow run says why (output is unaffected), else `None`.
pub fn gpu_fallback_summary() -> Option<String> {
    let [seed, extend, rescue] = GPU_FALLBACKS.each_ref().map(|c| c.load(Ordering::Relaxed));
    (seed + extend + rescue > 0).then(|| {
        format!(
            "GPU FALLBACK: {seed} seed batches, {extend} extension tiles and {rescue} rescue batches ran on the CPU because GPU calls failed (output identical, run slower; see the first GPU warning above)"
        )
    })
}

/// Before loading the index or starting prep: create the GPU seeding stream
/// and the pooled fused/two-pass and rescue engines with their streams, and
/// run an empty kernel on each, so the driver allocates its per-stream
/// context buffers while host memory is still contiguous. On 2026-10-22 a
/// 30x run started after memory had fragmented, the driver could not
/// allocate them (NV_ERR_NO_MEMORY), and every batch ran on the CPU.
/// One rescue engine per batch in flight and two fused lanes per batch, for
/// the run's batch count (`whole_genome` as for [`align_chunks_bam_files`]).
/// Returns the number of streams warmed.
pub fn gpu_warm_up(threads: usize, whole_genome: bool) -> Result<usize, nomics_gpu::GpuError> {
    nomics_gpu::warm_seeding()?;
    let mut warmed = 1;
    let inflight = batch_inflight(threads, whole_genome);
    let lanes = std::env::var("NOMICS_GPU_EXTEND_LANES")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .filter(|&n| n > 0)
        .unwrap_or(2);
    if gpu_rescue_enabled() {
        for _ in 0..(inflight * RescueSplit::from_env().engines_per_batch()).min(ENGINES_KEEP) {
            let mut engine = nomics_gpu::GpuKswAlign2::new();
            engine.warm()?;
            keep(&RESCUE_ENGINES, engine);
            warmed += 1;
        }
    }
    if gpu_extend_mode() == ExtendMode::Fused {
        for _ in 0..(inflight * lanes).min(ENGINES_KEEP) {
            let mut engine = crate::fused::FusedEngine::from_env();
            engine.warm()?;
            keep(&FUSED_ENGINES, engine);
            warmed += 1;
        }
    }
    Ok(warmed)
}

/// Unmap the alignment GPU buffers kept for reuse (seed buffers and the
/// fused, two-pass and rescue engines, several GB at 30x) once alignment is
/// done, so later stages do not carry them. Pooled buffers are idle: every
/// launch has been waited on before its engine returns to a pool.
pub fn release_gpu_buffers() {
    fn clear<T>(pool: &std::sync::Mutex<Vec<T>>) -> usize {
        let mut pool = pool.lock().unwrap_or_else(|e| e.into_inner());
        let n = pool.len();
        pool.clear();
        n
    }
    let seeds = clear(&SEED_POOL);
    SEED_SETS.fetch_sub(seeds, Ordering::Relaxed);
    clear(&FUSED_ENGINES);
    clear(&TWO_PASS_ENGINES);
    clear(&RESCUE_ENGINES);
}

/// Return an engine to its pool.
fn keep<T>(pool: &std::sync::Mutex<Vec<T>>, engine: T) {
    let mut pool = pool.lock().unwrap_or_else(|e| e.into_inner());
    if pool.len() < ENGINES_KEEP {
        pool.push(engine);
    }
}

/// An engine from its pool, or a new one.
fn take<T>(pool: &std::sync::Mutex<Vec<T>>, new: impl FnOnce() -> T) -> T {
    let engine = pool.lock().unwrap_or_else(|e| e.into_inner()).pop();
    engine.unwrap_or_else(new)
}

/// A batch's mate-rescue alignments, computed ahead of pairing, with each
/// pair's range of them.
struct BatchRescue {
    keys: Vec<RescueKey>,
    results: Vec<crate::ksw::Kswr>,
    pairs: Vec<(u32, u32)>,
}

impl BatchRescue {
    fn pair(&self, p: usize) -> Rescued<'_> {
        let (a, b) = self.pairs[p];
        Rescued {
            keys: &self.keys[a as usize..b as usize],
            results: &self.results[a as usize..b as usize],
        }
    }
}

/// Fraction of read pairs whose mate-rescue calls go to the GPU
/// (`NOMICS_GPU_RESCUE_SHARE`, 0..=1, default 1). The others are computed on
/// the CPU during pairing, as windows too long for the kernel already are, so
/// the output is the same. At 30x on Metal the GPU, not the CPU, is the bound
/// (docs/strategy/metal-sprint.md): this moves rescue work to idle cores.
fn gpu_rescue_share() -> f64 {
    static SHARE: std::sync::OnceLock<f64> = std::sync::OnceLock::new();
    *SHARE.get_or_init(|| {
        std::env::var("NOMICS_GPU_RESCUE_SHARE")
            .ok()
            .and_then(|v| v.parse::<f64>().ok())
            .filter(|f| (0.0..=1.0).contains(f))
            .unwrap_or(1.0)
    })
}

/// Plan every pair's rescue alignments (in parallel), compute them (by
/// default in one GPU launch; [`RescueSplit`] shares them with CPU workers),
/// and return them by pair. A failed GPU launch is computed on the CPU, with
/// identical results; calls not planned here are aligned during pairing.
fn gpu_rescue(
    idx: &Index,
    pes: &[PeStat; 4],
    reads: &[Read],
    regs: &[Vec<Region>],
    threads: usize,
) -> Option<BatchRescue> {
    struct Block {
        plans: Vec<RescuePlan>,
        bases: Vec<u8>,
        spans: Vec<(usize, usize, usize)>,
        pairs: Vec<(u32, u32)>,
    }
    let started = std::time::Instant::now();
    let n_blocks = reads.len().div_ceil(2 * BLOCK_PAIRS);
    let mut blocks: Vec<Block> = (0..n_blocks)
        .map(|_| Block {
            plans: Vec::new(),
            bases: Vec::new(),
            spans: Vec::new(),
            pairs: Vec::new(),
        })
        .collect();
    let work: Vec<_> = blocks
        .iter_mut()
        .zip(
            regs.chunks(2 * BLOCK_PAIRS)
                .zip(reads.chunks(2 * BLOCK_PAIRS)),
        )
        .collect();
    let share = gpu_rescue_share();
    for_each_parallel(threads, work, |(block, (rs, ss))| {
        for (p, (a, s)) in rs.chunks_exact(2).zip(ss.chunks_exact(2)).enumerate() {
            let first = block.plans.len() as u32;
            // Pair p of the block goes to the GPU when floor((p+1)f) >
            // floor(pf): the CPU's pairs spread evenly through the block.
            if share < 1.0 && ((p + 1) as f64 * share).floor() == (p as f64 * share).floor() {
                block.pairs.push((first, first));
                continue;
            }
            plan_rescue(
                idx,
                pes,
                s,
                a,
                &mut block.plans,
                &mut block.bases,
                &mut block.spans,
            );
            // Windows too long for the kernel stay with the CPU.
            let mut k = first as usize;
            while k < block.plans.len() {
                if block.spans[k].2 > GPU_RESCUE_MAX_TARGET {
                    block.plans.remove(k);
                    block.spans.remove(k);
                } else {
                    k += 1;
                }
            }
            block.pairs.push((first, block.plans.len() as u32));
        }
    });
    let planned = started.elapsed().as_secs_f64();
    let mut inputs = Vec::new();
    let mut keys = Vec::new();
    let mut pairs = Vec::with_capacity(reads.len() / 2);
    for block in &blocks {
        let base = keys.len() as u32;
        for (plan, &(at, ql, tl)) in block.plans.iter().zip(&block.spans) {
            keys.push(plan.key);
            inputs.push(nomics_gpu::GpuKswAlign2Input {
                query: &block.bases[at..at + ql],
                target: &block.bases[at + ql..at + ql + tl],
                m: 5,
                mat: &crate::extend::MAT,
                o_del: crate::extend::O_DEL,
                e_del: crate::extend::E_DEL,
                o_ins: crate::extend::O_INS,
                e_ins: crate::extend::E_INS,
                xtra: plan.xtra,
                size: if plan.xtra & crate::ksw::KSW_XBYTE != 0 {
                    1
                } else {
                    2
                },
            });
        }
        pairs.extend(block.pairs.iter().map(|&(a, b)| (base + a, base + b)));
    }
    let (results, stats) = rescue_split(&inputs, threads, RescueSplit::from_env());
    if std::env::var_os("NOMICS_TIMING").is_some() {
        eprintln!(
            "[timing] gpu rescue: {} calls ({} on the GPU in {} launches, {} on the CPU), plan {planned:.2}s, pack {:.2}s, kernel {:.2}s, total {:.2}s",
            inputs.len(),
            stats.gpu_jobs,
            stats.launches,
            inputs.len() - stats.gpu_jobs,
            stats.pack,
            stats.kernel,
            started.elapsed().as_secs_f64()
        );
    }
    dump_rescue_batch(&inputs, &results);
    Some(BatchRescue {
        keys,
        results,
        pairs,
    })
}

/// How a batch's planned rescue calls are shared between the GPU and CPU
/// workers. Every call is independent and `ksw_align2` gives the same
/// `Kswr` on either side, so the split never changes a result.
#[derive(Clone, Copy, Debug)]
struct RescueSplit {
    /// CPU workers computing calls next to the GPU
    /// (`NOMICS_RESCUE_CPU_THREADS`, at most the batch's workers; default 0,
    /// or all the batch's workers with a CPU share).
    cpu_threads: Option<usize>,
    /// With `Some(f)` (`NOMICS_RESCUE_CPU_SHARE`, 0..=1) the longest `f` of
    /// the calls are fixed to the CPU workers and the rest to the GPU, with
    /// no stealing; `1` puts every call on the CPU.
    cpu_share: Option<f64>,
    /// Calls per GPU launch (`NOMICS_RESCUE_GPU_TILE`; default the GPU's
    /// whole part without stealing, at least [`RESCUE_GPU_TILE`] with it).
    /// GPU lanes take tiles from the short end while CPU workers take
    /// [`RESCUE_CPU_CHUNK`] calls at a time from the long end, until they meet.
    gpu_tile: Option<usize>,
    /// With stealing, a GPU launch takes this fraction of the calls still
    /// unclaimed, if more than `gpu_tile` (`NOMICS_RESCUE_GPU_GUIDED`, 0..=1,
    /// default [`RESCUE_GPU_GUIDED`]; 0 gives fixed tiles).
    gpu_guided: Option<f64>,
    /// GPU launches in flight per batch (`NOMICS_RESCUE_GPU_LANES`), each on
    /// its own pooled engine and stream: one packs while another runs.
    gpu_lanes: usize,
}

/// Calls a CPU rescue worker takes at a time.
const RESCUE_CPU_CHUNK: usize = 256;

/// Smallest GPU launch when CPU workers steal from the batch.
const RESCUE_GPU_TILE: usize = 16_384;

/// Fraction of the unclaimed calls a GPU launch takes when CPU workers
/// steal: launches of 3/4, 3/16, ... of the batch.
const RESCUE_GPU_GUIDED: f64 = 0.75;

impl RescueSplit {
    fn from_env() -> Self {
        static SPLIT: std::sync::OnceLock<RescueSplit> = std::sync::OnceLock::new();
        *SPLIT.get_or_init(|| {
            let num = |name: &str| {
                std::env::var(name)
                    .ok()
                    .and_then(|s| s.parse::<usize>().ok())
            };
            RescueSplit {
                cpu_threads: num("NOMICS_RESCUE_CPU_THREADS"),
                cpu_share: std::env::var("NOMICS_RESCUE_CPU_SHARE")
                    .ok()
                    .and_then(|s| s.parse::<f64>().ok())
                    .filter(|f| (0.0..=1.0).contains(f)),
                gpu_tile: num("NOMICS_RESCUE_GPU_TILE").filter(|&n| n > 0),
                gpu_guided: std::env::var("NOMICS_RESCUE_GPU_GUIDED")
                    .ok()
                    .and_then(|s| s.parse::<f64>().ok())
                    .filter(|f| (0.0..=1.0).contains(f)),
                gpu_lanes: num("NOMICS_RESCUE_GPU_LANES")
                    .filter(|&n| n > 0)
                    .unwrap_or(2),
            }
        })
    }

    /// Rescue engines one batch may hold at once (for the warm-up).
    fn engines_per_batch(&self) -> usize {
        let tiled = self.gpu_tile.is_some()
            || self.cpu_share.is_some()
            || self.cpu_threads.is_some_and(|n| n > 0);
        if tiled { self.gpu_lanes } else { 1 }
    }
}

#[derive(Default)]
struct RescueStats {
    gpu_jobs: usize,
    launches: usize,
    pack: f64,
    kernel: f64,
}

/// One rescue call on the CPU: BWA's `ksw_align2`, as `matesw` makes it.
fn rescue_cpu(input: &nomics_gpu::GpuKswAlign2Input<'_>) -> crate::ksw::Kswr {
    crate::ksw::align2(
        input.query.len() as i32,
        input.query,
        input.target.len() as i32,
        input.target,
        input.m,
        input.mat,
        input.o_del,
        input.e_del,
        input.o_ins,
        input.e_ins,
        input.xtra,
        None,
    )
}

/// Every call of `inputs`, in input order, computed on the GPU and/or up to
/// `split.cpu_threads` CPU workers. Calls are ordered by target length; GPU
/// lanes take tiles from the short end and CPU workers chunks from the long
/// end until they meet (with a fixed share, each side drains its own part).
/// Each result is stored at its call's index, so the output does not depend
/// on who computed what or in which order. A failed GPU launch is computed
/// on the CPU by its lane, which stays on the CPU for the rest of the batch.
fn rescue_split(
    inputs: &[nomics_gpu::GpuKswAlign2Input<'_>],
    threads: usize,
    split: RescueSplit,
) -> (Vec<crate::ksw::Kswr>, RescueStats) {
    use std::sync::Mutex;
    type Piece = (usize, Vec<crate::ksw::Kswr>);
    let n = inputs.len();
    let mut order: Vec<u32> = (0..n as u32).collect();
    order.sort_by_key(|&i| inputs[i as usize].target.len());
    let order = &order;
    let all_cpu = split.cpu_share == Some(1.0);
    let default_cpu = if split.cpu_share.is_some() {
        threads
    } else {
        0
    };
    let cpu_threads = split.cpu_threads.unwrap_or(default_cpu).min(threads.max(1));
    let cpu_threads = if all_cpu {
        cpu_threads.max(1)
    } else {
        cpu_threads
    };
    // Each queue is the unclaimed range [lo, hi) of `order`: GPU lanes take
    // from lo, CPU workers from hi. A fixed share gives the CPU its own.
    let fixed = split.cpu_share.is_some() && cpu_threads > 0;
    let cut = if fixed {
        n - ((n as f64 * split.cpu_share.unwrap()).round() as usize).min(n)
    } else {
        n
    };
    // GPU claims: the whole GPU part in one launch without stealing; with
    // stealing, `guided` of the unclaimed range at a time (at least `tile`),
    // so the GPU makes a few large launches and the CPU workers take what
    // the last ones leave.
    let stealing = cpu_threads > 0 && !fixed;
    let tile = split
        .gpu_tile
        .unwrap_or(if stealing { RESCUE_GPU_TILE } else { cut })
        .max(1);
    let guided = split
        .gpu_guided
        .unwrap_or(if stealing { RESCUE_GPU_GUIDED } else { 0.0 });
    let gpu_lanes = if all_cpu || cut == 0 {
        0
    } else {
        split.gpu_lanes.min(cut.div_ceil(tile))
    };
    let queues = [Mutex::new((0, cut)), Mutex::new((cut, n))];
    let cpu_queue = &queues[usize::from(fixed)];
    let stats = Mutex::new(RescueStats::default());
    let failed = std::sync::atomic::AtomicBool::new(false);
    let pieces: Mutex<Vec<Piece>> = Mutex::new(Vec::new());
    let cpu_range = |range: std::ops::Range<usize>| -> Piece {
        let start = range.start;
        (
            start,
            range
                .map(|k| rescue_cpu(&inputs[order[k] as usize]))
                .collect(),
        )
    };
    std::thread::scope(|scope| {
        for _ in 0..gpu_lanes {
            scope.spawn(|| {
                // Each lane takes its own engine (buffers and stream) from the
                // pool, so lanes and batches in flight launch concurrently.
                let mut gpu = Some(take(&RESCUE_ENGINES, nomics_gpu::GpuKswAlign2::new));
                let mut tile_inputs = Vec::new();
                let mut out = Vec::new();
                loop {
                    let range = {
                        let mut q = queues[0].lock().unwrap();
                        if q.0 >= q.1 {
                            break;
                        }
                        let left = q.1 - q.0;
                        let want = tile.max((left as f64 * guided) as usize);
                        let hi = q.0 + want.min(left);
                        let range = q.0..hi;
                        q.0 = hi;
                        range
                    };
                    let Some(engine) = gpu.as_mut() else {
                        let piece = cpu_range(range);
                        pieces.lock().unwrap().push(piece);
                        continue;
                    };
                    tile_inputs.clear();
                    tile_inputs.extend(range.clone().map(|k| inputs[order[k] as usize]));
                    out.clear();
                    out.resize(tile_inputs.len(), nomics_gpu::GpuKswResult::default());
                    match engine.compute(&tile_inputs, &mut out) {
                        Ok(t) => {
                            let mut s = stats.lock().unwrap();
                            s.gpu_jobs += tile_inputs.len();
                            s.launches += 1;
                            s.pack += t.prep.as_secs_f64();
                            s.kernel += t.kernel.as_secs_f64();
                            drop(s);
                            let results = out.iter().map(gpu_kswr).collect();
                            pieces.lock().unwrap().push((range.start, results));
                        }
                        Err(e) => {
                            nomics_gpu::diag_gpu_error("GPU mate rescue", &e);
                            static WARNED: std::sync::Once = std::sync::Once::new();
                            WARNED.call_once(|| {
                                eprintln!("GPU mate rescue failed ({e:?}); rescuing on the CPU")
                            });
                            failed.store(true, Ordering::Relaxed);
                            keep(&RESCUE_ENGINES, gpu.take().unwrap());
                            let piece = cpu_range(range);
                            pieces.lock().unwrap().push(piece);
                        }
                    }
                }
                if let Some(engine) = gpu {
                    keep(&RESCUE_ENGINES, engine);
                }
            });
        }
        for _ in 0..cpu_threads {
            scope.spawn(|| {
                loop {
                    let range = {
                        let mut q = cpu_queue.lock().unwrap();
                        if q.0 >= q.1 {
                            break;
                        }
                        let lo = q.1.saturating_sub(RESCUE_CPU_CHUNK).max(q.0);
                        let range = lo..q.1;
                        q.1 = lo;
                        range
                    };
                    let piece = cpu_range(range);
                    pieces.lock().unwrap().push(piece);
                }
            });
        }
    });
    if failed.into_inner() {
        GPU_FALLBACKS[2].fetch_add(1, Ordering::Relaxed);
    }
    // A side with no worker (e.g. a fixed share with no GPU tile to run)
    // leaves its range unclaimed; compute it here.
    for queue in &queues {
        let (lo, hi) = *queue.lock().unwrap();
        if lo < hi {
            let piece = cpu_range(lo..hi);
            pieces.lock().unwrap().push(piece);
        }
    }
    let mut results = vec![crate::ksw::Kswr::default(); n];
    let mut filled = 0;
    for (start, out) in pieces.into_inner().unwrap() {
        filled += out.len();
        for (k, r) in out.into_iter().enumerate() {
            results[order[start + k] as usize] = r;
        }
    }
    assert_eq!(filled, n, "every rescue call is computed exactly once");
    (results, stats.into_inner().unwrap())
}

fn gpu_kswr(r: &nomics_gpu::GpuKswResult) -> crate::ksw::Kswr {
    crate::ksw::Kswr {
        score: r.score,
        te: r.te,
        qe: r.qe,
        score2: r.score2,
        te2: r.te2,
        tb: r.tb,
        qb: r.qb,
    }
}

/// `NOMICS_DUMP_RESCUE=<file>`: write the first rescue batch's calls and
/// results in the `align2.bin` format of `NOMICS_DUMP_KSW` (`ksw.rs`), a
/// whole-batch replay for `ksw_bench` with real launch shapes.
fn dump_rescue_batch(inputs: &[nomics_gpu::GpuKswAlign2Input<'_>], out: &[crate::ksw::Kswr]) {
    static ONCE: std::sync::Once = std::sync::Once::new();
    let Some(path) = std::env::var_os("NOMICS_DUMP_RESCUE") else {
        return;
    };
    ONCE.call_once(|| {
        let mut b = Vec::new();
        let ints = |b: &mut Vec<u8>, v: &[i32]| {
            v.iter().for_each(|x| b.extend_from_slice(&x.to_le_bytes()))
        };
        for (input, r) in inputs.iter().zip(out) {
            for s in [input.query, input.target] {
                b.extend_from_slice(&(s.len() as u32).to_le_bytes());
                b.extend_from_slice(s);
            }
            ints(&mut b, &[input.m]);
            b.extend(input.mat.iter().map(|&x| x as u8));
            ints(
                &mut b,
                &[
                    input.o_del,
                    input.e_del,
                    input.o_ins,
                    input.e_ins,
                    input.xtra,
                    i32::from(input.size),
                ],
            );
            ints(&mut b, &[r.score, r.te, r.qe, r.score2, r.te2, r.tb, r.qb]);
        }
        std::fs::write(&path, b).expect("write rescue dump");
        eprintln!(
            "[dump] {} rescue calls to {}",
            inputs.len(),
            path.to_string_lossy()
        );
    });
}

/// How seed extensions (`ksw_extend2`) are computed ahead of `chain2aln`.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum ExtendMode {
    /// On the CPU, inside `chain2aln`: `NOMICS_GPU_EXTEND=0`, and the
    /// default without `cuda` (where fused would run its CPU stand-in).
    Off,
    /// `NOMICS_GPU_EXTEND=1` (any value but `0` and `fused`): first-band left
    /// then right launches on the GPU, with GPU seeding.
    TwoPass,
    /// `NOMICS_GPU_EXTEND=fused`, and the default with `cuda` (30x on the
    /// Spark: 3,818 s vs 4,804 s for `Off`, both 41/41): one launch per tile
    /// of both sides with band retries ([`crate::fused`]), on the GPU or
    /// (`NOMICS_EXTEND_ENGINE=cpu`, or without `cuda`) the CPU stand-in.
    /// Runs with CPU seeds too.
    Fused,
}

/// Without `NOMICS_GPU_EXTEND`: fused where its kernel exists.
const DEFAULT_EXTEND_MODE: ExtendMode = if cfg!(feature = "cuda") {
    ExtendMode::Fused
} else {
    ExtendMode::Off
};

fn gpu_extend_mode() -> ExtendMode {
    static MODE: std::sync::OnceLock<ExtendMode> = std::sync::OnceLock::new();
    *MODE.get_or_init(|| match std::env::var("NOMICS_GPU_EXTEND").as_deref() {
        Err(_) => DEFAULT_EXTEND_MODE,
        Ok("0") => ExtendMode::Off,
        Ok("fused") => ExtendMode::Fused,
        Ok(_) => ExtendMode::TwoPass,
    })
}

/// Results stay in their planning blocks: no serial flattening of millions
/// of keys or per-read offsets. Blocks and reads retain input order.
struct ExtendBlock {
    keys: Vec<crate::extend::ExtKey>,
    results: Vec<crate::ksw::Extend>,
    reads: Vec<(u32, u32)>,
}

struct BatchExtend {
    blocks: Vec<ExtendBlock>,
}

impl BatchExtend {
    fn read(&self, i: usize) -> (&[crate::extend::ExtKey], &[crate::ksw::Extend]) {
        let block = &self.blocks[i / (2 * BLOCK_PAIRS)];
        let (a, b) = block.reads[i % (2 * BLOCK_PAIRS)];
        (
            &block.keys[a as usize..b as usize],
            &block.results[a as usize..b as usize],
        )
    }
}

/// Plan and pack independent blocks in parallel. The private engine lets
/// other tiles pack/execute concurrently, without a global extension lock.
/// `None` on failure leaves the original CPU extension fallback in charge.
fn gpu_extend_side(
    gpu: &mut nomics_gpu::GpuKswExtend2,
    reads: usize,
    threads: usize,
    plan: impl Fn(usize, &mut Vec<crate::extend::ExtPlan>, &mut Vec<u8>, &mut Vec<u8>) + Sync,
) -> Option<BatchExtend> {
    use crate::extend::{E_DEL, E_INS, ExtPlan, MAT, O_DEL, O_INS, ZDROP};
    let started = std::time::Instant::now();
    struct Block {
        plans: Vec<ExtPlan>,
        bases: Vec<u8>,
        reads: Vec<(u32, u32)>,
    }
    let per = 2 * BLOCK_PAIRS;
    let mut blocks: Vec<Block> = (0..reads.div_ceil(per))
        .map(|_| Block {
            plans: Vec::new(),
            bases: Vec::new(),
            reads: Vec::new(),
        })
        .collect();
    let work: Vec<_> = blocks.iter_mut().enumerate().collect();
    for_each_parallel(threads, work, |(b, block)| {
        let mut rseq = Vec::new();
        for i in b * per..((b + 1) * per).min(reads) {
            let first = block.plans.len() as u32;
            plan(i, &mut block.plans, &mut block.bases, &mut rseq);
            block.reads.push((first, block.plans.len() as u32));
        }
    });
    let planned = started.elapsed();
    // Iterators borrow the planners' arenas; no GpuExtendInput vector or
    // serial key/by_read concatenation. The engine packs each block directly.
    let groups: Vec<_> = blocks
        .iter()
        .map(|block| {
            block.plans.iter().map(|p| nomics_gpu::GpuExtendInput {
                query: &block.bases[p.at..p.at + p.qlen],
                target: &block.bases[p.at + p.qlen..p.at + p.qlen + p.tlen],
                m: 5,
                mat: &MAT,
                o_del: O_DEL,
                e_del: E_DEL,
                o_ins: O_INS,
                e_ins: E_INS,
                w: p.key.w,
                end_bonus: p.end_bonus,
                zdrop: ZDROP,
                h0: p.key.h0,
            })
        })
        .collect();
    let computed = std::time::Instant::now();
    let (mut out, timings) = match gpu.compute_groups(&groups, threads) {
        Ok(t) => t,
        Err(e) => {
            nomics_gpu::diag_gpu_error("GPU two-pass extension", &e);
            static WARNED: std::sync::Once = std::sync::Once::new();
            WARNED.call_once(|| eprintln!("GPU extension failed ({e:?}); extending on the CPU"));
            GPU_FALLBACKS[1].fetch_add(1, Ordering::Relaxed);
            return None;
        }
    };
    let compute_wall = computed.elapsed();
    let calls = out.len();
    drop(groups);
    let post = std::time::Instant::now();
    let mut result: Vec<_> = (0..blocks.len()).map(|_| None).collect();
    let work: Vec<_> = blocks
        .into_iter()
        .zip(&mut result)
        .map(|(block, slot)| {
            let (values, rest) = out.split_at(block.plans.len());
            out = rest;
            (block, slot, values)
        })
        .collect();
    for_each_parallel(threads, work, |(block, slot, values)| {
        *slot = Some(ExtendBlock {
            keys: block.plans.into_iter().map(|p| p.key).collect(),
            results: values
                .iter()
                .map(|r| crate::ksw::Extend {
                    score: r.score,
                    qle: r.qle,
                    tle: r.tle,
                    gtle: r.gtle,
                    gscore: r.gscore,
                    max_off: r.max_off,
                })
                .collect(),
            reads: block.reads,
        });
        // Drop the plan/base arenas on the workers, included in post time.
    });
    let result = BatchExtend {
        blocks: result.into_iter().map(Option::unwrap).collect(),
    };
    if std::env::var_os("NOMICS_TIMING").is_some() {
        eprintln!(
            "[timing] gpu extend: {calls} calls, plan {:.3}s, pack {:.3}s, kernel {:.3}s, execute {:.3}s, post {:.3}s, lock-wait none (private engine), total {:.3}s",
            planned.as_secs_f64(),
            timings.prep.as_secs_f64(),
            timings.kernel.as_secs_f64(),
            compute_wall.saturating_sub(timings.prep).as_secs_f64(),
            post.elapsed().as_secs_f64(),
            started.elapsed().as_secs_f64()
        );
    }
    Some(result)
}

type ExtendWork<'a> = Vec<(usize, (&'a mut [Vec<Region>], &'a [Read]))>;

/// Two bounded tiles progress independently through chains -> left -> right
/// -> chain2aln (two-pass), or chains -> fused -> chain2aln. In the two-pass
/// path a tile's right h0 still comes only from its left result. The
/// original -K batch and output slots stay intact for pestat and pairing.
fn align1_gpu_extend(
    idx: &Index,
    reads: &[Read],
    seeds: Option<&GpuSeeds>,
    threads: usize,
    work: ExtendWork<'_>,
    mode: ExtendMode,
) -> f64 {
    let positive_env = |name: &str, default: usize| {
        std::env::var(name)
            .ok()
            .and_then(|s| s.parse::<usize>().ok())
            .filter(|&n| n > 0)
            .unwrap_or(default)
    };
    let per = 2 * BLOCK_PAIRS;
    let tile_blocks = positive_env("NOMICS_GPU_EXTEND_READS", 65_536)
        .min(reads.len().max(1))
        .div_ceil(per);
    let split = Split {
        cpu_lanes: std::env::var("NOMICS_GPU_EXTEND_CPU_LANES")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(0),
        cpu_share: std::env::var("NOMICS_GPU_EXTEND_CPU_SHARE")
            .ok()
            .and_then(|s| s.parse::<f64>().ok())
            .filter(|f| (0.0..=1.0).contains(f)),
    };
    let busy = align1_gpu_extend_tiled(
        idx,
        reads,
        seeds,
        threads,
        work,
        tile_blocks,
        positive_env("NOMICS_GPU_EXTEND_LANES", DEFAULT_GPU_EXTEND_LANES),
        mode,
        split,
    );
    if std::env::var_os("NOMICS_TIMING").is_some() {
        use crate::extend::ext_stats;
        use std::sync::atomic::Ordering::Relaxed;
        ext_stats::ON.store(true, Relaxed);
        eprintln!(
            "[timing] gpu extend lookups so far: {} hits, {} misses, {} CPU fallback seeds (query > {}), {} band retries in chain2aln, {} fused jobs, {} fused band retries",
            ext_stats::HITS.load(Relaxed),
            ext_stats::MISSES.load(Relaxed),
            ext_stats::FALLBACKS.load(Relaxed),
            nomics_gpu::EXTEND_MAX_QUERY,
            ext_stats::RETRIES.load(Relaxed),
            ext_stats::FUSED_JOBS.load(Relaxed),
            ext_stats::FUSED_RETRIES.load(Relaxed),
        );
    }
    busy
}

/// GPU extension lanes (tiles in flight) when `NOMICS_GPU_EXTEND_LANES` is
/// unset. Metal has no queue priorities, so a lane's tile can wait behind
/// seeding and rescue work; more tiles in flight keep its single GPU queue
/// full. Tier1 on an M4 Max with fused extension (docs/strategy/metal-sprint.md):
/// align 314 s (2 lanes), 250 (4), 239 (5), 238/240 (6), 245 (7), 261 (8).
/// CUDA keeps the Spark's 2 (its extension stream has the highest priority).
const DEFAULT_GPU_EXTEND_LANES: usize = if cfg!(feature = "metal") { 6 } else { 2 };

/// A lane's private extension engine. `Cpu` lanes align their tiles
/// entirely on the CPU, as the default path does.
enum LaneEngine {
    TwoPass(nomics_gpu::GpuKswExtend2),
    Fused(crate::fused::FusedEngine),
    Cpu,
}

/// Lanes that extend on the CPU next to the GPU lanes
/// (`NOMICS_GPU_EXTEND_CPU_LANES`). With `cpu_share`
/// (`NOMICS_GPU_EXTEND_CPU_SHARE`, 0..=1) that fraction of the tiles is fixed
/// to the CPU lanes up front; otherwise every lane takes the next tile from
/// one queue (work stealing). Each read's regions are the same either way.
#[derive(Clone, Copy, Debug, Default)]
struct Split {
    cpu_lanes: usize,
    cpu_share: Option<f64>,
}

#[allow(clippy::too_many_arguments)]
fn align1_gpu_extend_tiled(
    idx: &Index,
    reads: &[Read],
    seeds: Option<&GpuSeeds>,
    threads: usize,
    work: ExtendWork<'_>,
    tile_blocks: usize,
    lanes: usize,
    mode: ExtendMode,
    split: Split,
) -> f64 {
    assert!(tile_blocks > 0 && lanes > 0);
    let per = 2 * BLOCK_PAIRS;
    let mut items = work.into_iter();
    let mut tiles = Vec::new();
    loop {
        let tile: Vec<_> = items.by_ref().take(tile_blocks).collect();
        if tile.is_empty() {
            break;
        }
        tiles.push(tile);
    }
    let gpu_lanes = lanes.min(threads.max(1)).min(tiles.len());
    let cpu_lanes = split.cpu_lanes.min(tiles.len());
    let lanes = gpu_lanes + cpu_lanes;
    if gpu_lanes == 0 {
        return 0.0;
    }
    // Fixed share: tile t goes to the CPU when floor((t+1)f) > floor(tf),
    // which spreads the CPU's tiles evenly through the batch.
    let (gpu_tiles, cpu_tiles): (Vec<_>, Vec<_>) = match split.cpu_share {
        Some(f) if cpu_lanes > 0 => {
            let mut parts = (Vec::new(), Vec::new());
            for (t, tile) in tiles.into_iter().enumerate() {
                if ((t + 1) as f64 * f).floor() > (t as f64 * f).floor() {
                    parts.1.push(tile);
                } else {
                    parts.0.push(tile);
                }
            }
            parts
        }
        _ => (tiles, Vec::new()),
    };
    let stealing = cpu_tiles.is_empty();
    let queues = [
        std::sync::Mutex::new(gpu_tiles.into_iter()),
        std::sync::Mutex::new(cpu_tiles.into_iter()),
    ];
    std::thread::scope(|scope| {
        let handles: Vec<_> = (0..lanes)
            .map(|lane| {
                let on_cpu = lane >= gpu_lanes;
                let queue = &queues[usize::from(on_cpu && !stealing)];
                // Split NOMICS_ALIGN_WORKERS rather than multiplying its CPU budget.
                let workers = threads.max(1) / lanes + usize::from(lane < threads.max(1) % lanes);
                let workers = workers.max(1);
                scope.spawn(move || {
                    // Engines outlive the batch: taken from and returned to
                    // the process-wide pools.
                    let mut engine = match mode {
                        _ if on_cpu => LaneEngine::Cpu,
                        ExtendMode::Fused => LaneEngine::Fused(take(
                            &FUSED_ENGINES,
                            crate::fused::FusedEngine::from_env,
                        )),
                        _ => LaneEngine::TwoPass(take(
                            &TWO_PASS_ENGINES,
                            nomics_gpu::GpuKswExtend2::new,
                        )),
                    };
                    let mut busy = 0.0;
                    loop {
                        let tile = queue.lock().unwrap().next();
                        let Some(mut tile) = tile else {
                            break;
                        };
                        let first_block = tile[0].0;
                        let first = first_block * per;
                        let n: usize = tile.iter().map(|(_, (_, rs))| rs.len()).sum();
                        for (b, _) in &mut tile {
                            *b -= first_block;
                        }
                        let reads = &reads[first..first + n];
                        busy += match &mut engine {
                            LaneEngine::TwoPass(gpu) => {
                                align1_gpu_extend_tile(idx, reads, seeds, first, workers, tile, gpu)
                            }
                            LaneEngine::Fused(engine) => align1_fused_extend_tile(
                                idx, reads, seeds, first, workers, tile, engine,
                            ),
                            LaneEngine::Cpu => {
                                align1_cpu_tile(idx, reads, seeds, first, workers, tile)
                            }
                        };
                    }
                    match engine {
                        LaneEngine::Fused(engine) => keep(&FUSED_ENGINES, engine),
                        LaneEngine::TwoPass(engine) => keep(&TWO_PASS_ENGINES, engine),
                        LaneEngine::Cpu => {}
                    }
                    busy
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join().unwrap()).sum()
    })
}

/// One tile aligned entirely on the CPU, as the default path aligns a batch.
fn align1_cpu_tile(
    idx: &Index,
    reads: &[Read],
    seeds: Option<&GpuSeeds>,
    first: usize,
    threads: usize,
    work: ExtendWork<'_>,
) -> f64 {
    let started = std::time::Instant::now();
    let busy = for_each_parallel(threads, work, |(b, (slots, rs))| {
        let mut aux = AlignAux::default();
        for (j, (slot, r)) in slots.iter_mut().zip(rs).enumerate() {
            let i = first + b * 2 * BLOCK_PAIRS + j;
            crate::diag::at_read(r);
            *slot = match seeds.and_then(|s| s.get(i)) {
                Some(mem) => align1_core_from(idx, &r.seq, mem, &mut aux),
                None => align1_core(idx, &r.seq, &mut aux),
            };
        }
    });
    if std::env::var_os("NOMICS_TIMING").is_some() {
        eprintln!(
            "[timing] cpu extend tile: {} reads, total {:.3}s",
            reads.len(),
            started.elapsed().as_secs_f64()
        );
    }
    busy
}

/// The filtered chains of every read of a tile, in parallel.
fn tile_chains(
    idx: &Index,
    reads: &[Read],
    seeds: Option<&GpuSeeds>,
    first: usize,
    threads: usize,
) -> (Vec<Vec<crate::chain::Chain>>, f64) {
    let mut chains: Vec<Vec<crate::chain::Chain>> = vec![Vec::new(); reads.len()];
    let chain_work: Vec<_> = chains
        .chunks_mut(2 * BLOCK_PAIRS)
        .zip(reads.chunks(2 * BLOCK_PAIRS))
        .enumerate()
        .collect();
    let busy = for_each_parallel(threads, chain_work, |(b, (slots, rs))| {
        let mut aux = AlignAux::default();
        for (j, (slot, r)) in slots.iter_mut().zip(rs).enumerate() {
            let mem = seeds.and_then(|s| s.get(first + b * 2 * BLOCK_PAIRS + j));
            *slot = crate::extend::align1_chains(idx, &r.seq, mem, &mut aux);
        }
    });
    (chains, busy)
}

/// One tile; seed indices are global within the original -K batch, while
/// chains and extension caches use tile-local indices.
#[allow(clippy::too_many_arguments)]
fn align1_gpu_extend_tile(
    idx: &Index,
    reads: &[Read],
    seeds: Option<&GpuSeeds>,
    first: usize,
    threads: usize,
    work: ExtendWork<'_>,
    gpu: &mut nomics_gpu::GpuKswExtend2,
) -> f64 {
    use crate::extend::{Extended, align1_extend, plan_left, plan_right};
    let (chains, mut busy) = tile_chains(idx, reads, seeds, first, threads);
    let left = gpu_extend_side(gpu, reads.len(), threads, |i, plans, bases, rseq| {
        plan_left(idx, &reads[i].seq, &chains[i], plans, bases, rseq)
    });
    let left_of = |i: usize| left.as_ref().map_or((&[][..], &[][..]), |l| l.read(i));
    let right = left.as_ref().and_then(|_| {
        gpu_extend_side(gpu, reads.len(), threads, |i, plans, bases, rseq| {
            let (k, r) = left_of(i);
            let pre = Extended::new(&chains[i], [k, &[]], [r, &[]]);
            plan_right(idx, &reads[i].seq, &chains[i], &pre, plans, bases, rseq)
        })
    });
    let right_of = |i: usize| right.as_ref().map_or((&[][..], &[][..]), |r| r.read(i));
    busy += for_each_parallel(threads, work, |(b, (slots, rs))| {
        let mut aux = AlignAux::default();
        for (j, (slot, r)) in slots.iter_mut().zip(rs).enumerate() {
            let i = b * 2 * BLOCK_PAIRS + j;
            let ((lk, lr), (rk, rr)) = (left_of(i), right_of(i));
            let pre = Extended::new(&chains[i], [lk, rk], [lr, rr]);
            crate::diag::at_read(r);
            *slot = align1_extend(idx, &r.seq, &chains[i], &mut aux, &pre);
        }
    });

    busy
}

/// A planning block's fused results as replay keys, per side, with each
/// read's `(first, end)` into them.
struct FusedBlock {
    keys: [Vec<crate::extend::ExtKey>; 2],
    results: [Vec<crate::ksw::Extend>; 2],
    reads: Vec<[(u32, u32); 2]>,
}

type FusedRead<'a> = (
    [&'a [crate::extend::ExtKey]; 2],
    [&'a [crate::ksw::Extend]; 2],
);

impl FusedBlock {
    fn read(&self, j: usize) -> FusedRead<'_> {
        let r = self.reads[j];
        let span = |s: usize| r[s].0 as usize..r[s].1 as usize;
        (
            [&self.keys[0][span(0)], &self.keys[1][span(1)]],
            [&self.results[0][span(0)], &self.results[1][span(1)]],
        )
    }
}

/// One tile through the fused path: chains -> plan_fused -> one engine launch
/// -> replay through `align1_extend`, with each read's `Extended` built from
/// its fused results. On engine failure the replay extends on the CPU.
#[allow(clippy::too_many_arguments)]
fn align1_fused_extend_tile(
    idx: &Index,
    reads: &[Read],
    seeds: Option<&GpuSeeds>,
    first: usize,
    threads: usize,
    work: ExtendWork<'_>,
    engine: &mut crate::fused::FusedEngine,
) -> f64 {
    use crate::extend::{Extended, FusedPlan, align1_extend, ext_stats, plan_fused};
    let started = std::time::Instant::now();
    let (chains, mut busy) = tile_chains(idx, reads, seeds, first, threads);
    let chained = started.elapsed();

    struct Plan {
        plans: Vec<FusedPlan>,
        bases: Vec<u8>,
        reads: Vec<(u32, u32)>,
    }
    let per = 2 * BLOCK_PAIRS;
    let mut blocks: Vec<Plan> = (0..reads.len().div_ceil(per))
        .map(|_| Plan {
            plans: Vec::new(),
            bases: Vec::new(),
            reads: Vec::new(),
        })
        .collect();
    let plan_work: Vec<_> = blocks.iter_mut().enumerate().collect();
    busy += for_each_parallel(threads, plan_work, |(b, block)| {
        let mut rseq = Vec::new();
        for i in b * per..((b + 1) * per).min(reads.len()) {
            let at = block.plans.len() as u32;
            plan_fused(
                idx,
                &reads[i].seq,
                &chains[i],
                &mut block.plans,
                &mut block.bases,
                &mut rseq,
            );
            block.reads.push((at, block.plans.len() as u32));
        }
    });
    let planned = started.elapsed();

    let params = crate::fused::params();
    let groups: Vec<_> = blocks
        .iter()
        .map(|block| {
            block.plans.iter().map(|p| {
                let (lq, lt, rq) = (p.lq as usize, p.lt as usize, p.rq as usize);
                let b = &block.bases[p.at..];
                nomics_gpu::GpuFusedInput {
                    lq: &b[..lq],
                    lt: &b[lq..lq + lt],
                    rq: &b[lq + lt..lq + lt + rq],
                    rt: &b[lq + lt + rq..lq + lt + rq + p.rt as usize],
                    h0: p.h0,
                }
            })
        })
        .collect();
    let computed = std::time::Instant::now();
    let (out, timings) = match engine.compute(&params, &groups, threads) {
        Ok(t) => t,
        Err(e) => {
            nomics_gpu::diag_gpu_error("GPU fused extension", &e);
            static WARNED: std::sync::Once = std::sync::Once::new();
            WARNED.call_once(|| eprintln!("fused extension failed ({e:?}); extending on the CPU"));
            GPU_FALLBACKS[1].fetch_add(1, Ordering::Relaxed);
            busy += for_each_parallel(threads, work, |(b, (slots, rs))| {
                let mut aux = AlignAux::default();
                for (j, (slot, r)) in slots.iter_mut().zip(rs).enumerate() {
                    let i = b * per + j;
                    crate::diag::at_read(r);
                    *slot = align1_extend(idx, &r.seq, &chains[i], &mut aux, &Extended::default());
                }
            });
            return busy;
        }
    };
    let compute_wall = computed.elapsed();
    crate::fused::dump::launch(&params, &groups, out);
    drop(groups);
    let jobs = out.len();
    let post = std::time::Instant::now();

    // Per block: fused results -> replay keys, in (chain, seed) order.
    let mut fused: Vec<Option<FusedBlock>> = (0..blocks.len()).map(|_| None).collect();
    let mut rest = out;
    let post_work: Vec<_> = blocks
        .into_iter()
        .zip(&mut fused)
        .map(|(block, slot)| {
            let (values, tail) = rest.split_at(block.plans.len());
            rest = tail;
            (block, slot, values)
        })
        .collect();
    busy += for_each_parallel(threads, post_work, |(block, slot, values)| {
        let mut keys = [Vec::new(), Vec::new()];
        let mut results = [Vec::new(), Vec::new()];
        let mut spans = Vec::with_capacity(block.reads.len());
        let mut retries = 0;
        for &(a, b) in &block.reads {
            let start = [keys[0].len() as u32, keys[1].len() as u32];
            let (a, b) = (a as usize, b as usize);
            for (p, r) in block.plans[a..b].iter().zip(&values[a..b]) {
                retries += u64::from(r.n_left == 2) + u64::from(r.n_right == 2);
                crate::fused::push_keys(&params, p.chain, p.seed, p.h0, r, &mut keys, &mut results);
            }
            spans.push([
                (start[0], keys[0].len() as u32),
                (start[1], keys[1].len() as u32),
            ]);
        }
        ext_stats::add(&ext_stats::FUSED_JOBS, values.len() as u64);
        ext_stats::add(&ext_stats::FUSED_RETRIES, retries);
        *slot = Some(FusedBlock {
            keys,
            results,
            reads: spans,
        });
        // Drop the plan/base arenas on the workers, included in post time.
    });
    let fused: Vec<FusedBlock> = fused.into_iter().map(Option::unwrap).collect();
    let post = post.elapsed();

    let replay = std::time::Instant::now();
    busy += for_each_parallel(threads, work, |(b, (slots, rs))| {
        let mut aux = AlignAux::default();
        for (j, (slot, r)) in slots.iter_mut().zip(rs).enumerate() {
            let i = b * per + j;
            let (keys, results) = fused[i / per].read(i % per);
            let pre = Extended::new(&chains[i], keys, results);
            crate::diag::at_read(r);
            *slot = align1_extend(idx, &r.seq, &chains[i], &mut aux, &pre);
        }
    });
    if std::env::var_os("NOMICS_TIMING").is_some() {
        eprintln!(
            "[timing] fused extend ({}): {jobs} jobs, chains {:.3}s, plan {:.3}s, pack {:.3}s, kernel {:.3}s, execute {:.3}s, post {:.3}s, replay {:.3}s, total {:.3}s",
            engine.name(),
            chained.as_secs_f64(),
            (planned - chained).as_secs_f64(),
            timings.prep.as_secs_f64(),
            timings.kernel.as_secs_f64(),
            compute_wall.saturating_sub(timings.prep).as_secs_f64(),
            post.as_secs_f64(),
            replay.elapsed().as_secs_f64(),
            started.elapsed().as_secs_f64()
        );
    }
    busy
}

/// Interval cap per read of the GPU's first seeding pass (the result stride
/// of the batch-sized buffer).
const GPU_SEED_CAP: u32 = 64;
/// Interval cap of the retry pass, which re-runs only the first pass's
/// overflow reads; reads that overflow it too fall back to the CPU.
const GPU_SEED_RETRY_CAP: u32 = 128;

/// GPU seeding work per batch, `NOMICS_SEED_WORK`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SeedMode {
    /// `cap64`: one 64-entry pass; every overflow read re-seeds on the CPU.
    Cap64,
    /// `two-pass` (the default, see [`DEFAULT_SEED_MODE`]): the 64-entry pass on the whole batch, then the
    /// 128-entry kernel on its overflow reads only, into a buffer sized for
    /// them (the work of `SeedWork::TwoPass64Then128` at stride 128, without
    /// a stride-128 buffer for the whole batch). Per-read results are pure
    /// functions of the read, so both modes give the same seeds; two-pass
    /// leaves fewer reads to the CPU (tier0 on Metal: 473
    /// instead of 3,421 of 991,118).
    TwoPass,
}

/// Default when `NOMICS_SEED_WORK` is unset, on every backend. Metal used
/// cap64 while each launch re-bound the 3.2 GB BWT (~85 ms); with the BWT
/// registered once (`nomics_gpu::GpuHostBuf`) two-pass is faster there too
/// (tier1 on an M4 Max: 372.0 s vs 380 s; docs/strategy/metal-sprint.md).
const DEFAULT_SEED_MODE: SeedMode = SeedMode::TwoPass;

fn parse_seed_mode(v: Option<&str>, default: SeedMode) -> SeedMode {
    match v {
        None => default,
        Some("cap64") => SeedMode::Cap64,
        Some("two-pass") => SeedMode::TwoPass,
        Some(other) => {
            eprintln!(
                "WARNING: NOMICS_SEED_WORK={other:?} is not cap64 or two-pass; using {default:?}"
            );
            default
        }
    }
}

fn seed_mode() -> SeedMode {
    static MODE: std::sync::OnceLock<SeedMode> = std::sync::OnceLock::new();
    *MODE.get_or_init(|| {
        parse_seed_mode(
            std::env::var("NOMICS_SEED_WORK").ok().as_deref(),
            DEFAULT_SEED_MODE,
        )
    })
}

/// `NOMICS_SEED_TABLE=<cache path>` (k = `NOMICS_SEED_TABLE_K`, default 12):
/// the GPU seed kernel starts SMEM forward walks and LAST seeds from this
/// exact k-mer table instead of the first k BWT steps (same intervals, fewer
/// BWT lookups; built and cached at that path if missing). `None` when unset
/// or when the table cannot be used for `idx` (with a warning: seeding is the
/// same without it).
fn seed_table(idx: &Index) -> Option<&'static crate::kmer_table::KmerTable> {
    let path = std::env::var_os("NOMICS_SEED_TABLE")?;
    let k = std::env::var("NOMICS_SEED_TABLE_K")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(12);
    seed_table_at(idx, std::path::Path::new(&path), k)
}

/// Identity of a loaded BWT: where its words live, its header, and a sample
/// of its words (so an index loaded later at a reused address with the same
/// header is still told apart).
#[derive(Clone, Copy, PartialEq, Eq)]
struct BwtKey {
    words: (usize, usize),
    primary: u64,
    seq_len: u64,
    l2: [u64; 5],
    sample: [u32; 16],
}

impl BwtKey {
    fn of(bwt: &crate::index::Bwt) -> Self {
        let words = bwt.words();
        let mut sample = [0u32; 16];
        if !words.is_empty() {
            for (i, s) in sample.iter_mut().enumerate() {
                *s = words[i * (words.len() - 1) / 15];
            }
        }
        Self {
            words: (words.as_ptr() as usize, words.len()),
            primary: bwt.primary,
            seq_len: bwt.seq_len,
            l2: bwt.l2,
            sample,
        }
    }
}

type SeedTableEntry = (
    BwtKey,
    std::path::PathBuf,
    usize,
    Option<&'static crate::kmer_table::KmerTable>,
);

/// The table at `path` for `idx`'s BWT, loaded (or built) once per index,
/// path and k. Tables are cached per BWT identity, never shared between
/// indexes; `KmerTable::load` also checks the file against the BWT (a
/// fingerprint of every word), so a cache built for another reference is
/// refused and the result is `None`. Loaded tables live for the process (one
/// per index in practice).
fn seed_table_at(
    idx: &Index,
    path: &std::path::Path,
    k: usize,
) -> Option<&'static crate::kmer_table::KmerTable> {
    static TABLES: std::sync::Mutex<Vec<SeedTableEntry>> = std::sync::Mutex::new(Vec::new());
    let key = BwtKey::of(&idx.bwt);
    let mut tables = TABLES.lock().unwrap_or_else(|e| e.into_inner());
    if let Some(e) = tables
        .iter()
        .find(|e| e.0 == key && e.1 == path && e.2 == k)
    {
        return e.3;
    }
    let table = match crate::kmer_table::KmerTable::load_or_build(path, &idx.bwt, k) {
        Ok(table) => Some(&*Box::leak(Box::new(table))),
        Err(e) => {
            eprintln!(
                "WARNING: GPU seed k-mer table {} not used ({e}); seeding is unchanged",
                path.display()
            );
            None
        }
    };
    tables.push((key, path.to_path_buf(), k, table));
    table
}

/// Retry rows a seed set holds for `overflow` first-pass overflow reads in a
/// set of `capacity` reads: 25% headroom, at least 1/64 of the set (fewer
/// regrowths), at most 1/8 (bounded memory; overflow reads past that re-seed
/// on the CPU, which is still exact). 30x: ~4% of reads overflow.
fn seed_retry_rows(overflow: usize, capacity: usize) -> usize {
    (overflow + overflow / 4)
        .max(capacity / 64)
        .min(capacity / 8)
        .max(1)
}

/// A finished batch's GPU seeds: the first pass, and the retry of its
/// overflow reads (`retried[j]` is the read of retry row `j`, ascending).
struct SeedResults<'a> {
    first: &'a nomics_gpu::GpuSeedBuffers,
    retried: &'a [u32],
    retry: Option<&'a nomics_gpu::GpuSeedBuffers>,
}

impl SeedResults<'_> {
    fn get(&self, read: usize) -> Option<&[nomics_gpu::GpuIntv]> {
        self.first.get(read).or_else(|| {
            let row = self.retried.binary_search(&(read as u32)).ok()?;
            self.retry?.get(row)
        })
    }
}

/// Seed intervals of every read, computed on the GPU. `get(i)` is `None`
/// for reads the GPU could not finish (the CPU recomputes those).
///
/// Ordinary CPU memory, compacted from the GPU's fixed-stride result buffer
/// (cap 64 intervals per read, ~1.7 GB per 30x batch) once the kernel has
/// finished, so that buffer goes back to [`SEED_POOL`] at once and a batch
/// holds only the intervals it has (0.25 GB per 30x batch; the copy takes
/// ~25 ms on the seeder thread).
pub struct GpuSeeds {
    intvs: Vec<Intv>,
    /// Each read's range in `intvs`; `start == NO_SEEDS` for a CPU fallback.
    spans: Vec<(u64, u64)>,
}

const NO_SEEDS: u64 = u64::MAX;

impl GpuSeeds {
    pub fn get(&self, i: usize) -> Option<&[Intv]> {
        let (a, b) = self.spans[i];
        (a != NO_SEEDS).then(|| &self.intvs[a as usize..b as usize])
    }

    /// Bytes held (diagnostics).
    fn bytes(&self) -> usize {
        self.intvs.capacity() * std::mem::size_of::<Intv>()
            + self.spans.capacity() * std::mem::size_of::<(u64, u64)>()
    }

    /// The first `n` reads' intervals from a batch's finished GPU results.
    fn compact(out: &SeedResults<'_>, n: usize) -> Self {
        let total: usize = (0..n).filter_map(|i| out.get(i)).map(<[_]>::len).sum();
        let mut intvs = Vec::with_capacity(total);
        let spans = (0..n)
            .map(|i| match out.get(i) {
                Some(g) => {
                    let a = intvs.len() as u64;
                    intvs.extend(g.iter().map(|g| Intv {
                        x: g.x,
                        info: g.info,
                    }));
                    (a, intvs.len() as u64)
                }
                None => (NO_SEEDS, 0),
            })
            .collect();
        Self { intvs, spans }
    }

    /// Test seeds, read by read (`None`: falls back to the CPU).
    #[cfg(test)]
    fn from_lists(lists: &[Option<Vec<Intv>>]) -> Self {
        let mut intvs = Vec::new();
        let spans = lists
            .iter()
            .map(|l| match l {
                Some(l) => {
                    let a = intvs.len() as u64;
                    intvs.extend_from_slice(l);
                    (a, intvs.len() as u64)
                }
                None => (NO_SEEDS, 0),
            })
            .collect();
        Self { intvs, spans }
    }
}

/// Return a seed buffer set to [`SEED_POOL`] (or unmap it if the pool is
/// full). Only sets whose GPU work has been waited on come here.
fn release_seed_bufs(bufs: SeedBufs) {
    let mut pool = SEED_POOL.lock().unwrap_or_else(|e| e.into_inner());
    if pool.len() < SEED_POOL_KEEP {
        pool.push(bufs);
    } else {
        SEED_SETS.fetch_sub(1, Ordering::Relaxed);
    }
}

/// One batch's GPU seeding buffers: packed read codes, each read's offset
/// and length, and the results, all huge-page host memory that the kernel
/// reads and writes in place. Batches reuse them through [`SEED_POOL`],
/// because the GPU pays first-touch translation work on fresh pages (~1 s
/// per 30x batch, more than the kernel itself). Nothing the kernel reads
/// lives in allocator memory: recycled 4 KiB heap pages are what the HostBuf
/// rule (docs/METHODOLOGY.md) keeps away from the GPU.
struct SeedBufs {
    // Registered: the GPU binds these pooled buffers without re-wrapping them.
    seqs: nomics_gpu::GpuHostBuf,
    offs: nomics_gpu::GpuHostBuf,
    lens: nomics_gpu::GpuHostBuf,
    out: nomics_gpu::GpuSeedBuffers,
    /// Two-pass retry rows, allocated on first need and grown on demand
    /// (see [`seed_retry_rows`]).
    retry: Option<SeedRetry>,
    /// Reads of the current batch in the retry, ascending (CPU only).
    retried: Vec<u32>,
}

/// The retry pass's inputs and results for up to `out.capacity()` reads.
/// `offs`/`lens` point into the batch's `seqs`, so no sequence is copied.
struct SeedRetry {
    offs: nomics_mem::HostBuf,
    lens: nomics_mem::HostBuf,
    out: nomics_gpu::GpuSeedBuffers,
}

impl SeedRetry {
    fn new(rows: usize) -> Result<Self, nomics_gpu::GpuError> {
        let host = |e: std::io::Error| nomics_gpu::GpuError::HostMemory(e.to_string());
        let mut offs = nomics_mem::HostBuf::new(rows * 8).map_err(host)?;
        let mut lens = nomics_mem::HostBuf::new(rows * 4).map_err(host)?;
        let mut out = nomics_gpu::GpuSeedBuffers::new(rows, GPU_SEED_RETRY_CAP)?;
        offs.touch_huge("GPU seed retry offsets");
        lens.touch_huge("GPU seed retry lengths");
        out.touch_and_collapse();
        Ok(Self { offs, lens, out })
    }
}

/// Buffers kept for reuse. A set is live only while its batch seeds (the
/// results are then compacted into a [`GpuSeeds`]), and the seeder runs one
/// launch at a time, so one set is in use; the second slot keeps a larger
/// batch's set too. A dropped set has always been waited on, so unmapping it
/// is safe.
const SEED_POOL_KEEP: usize = 2;
static SEED_POOL: std::sync::Mutex<Vec<SeedBufs>> = std::sync::Mutex::new(Vec::new());

/// Pooled buffers for `reads` reads of `bases` bases in total (new ones,
/// with 25% headroom, if no pooled set is big enough).
fn seed_bufs(reads: usize, bases: usize) -> Result<SeedBufs, nomics_gpu::GpuError> {
    {
        let mut pool = SEED_POOL.lock().unwrap_or_else(|e| e.into_inner());
        if let Some(at) = pool
            .iter()
            .position(|b| b.out.capacity() >= reads && b.seqs.len() >= bases)
        {
            return Ok(pool.swap_remove(at));
        }
    }
    let host = |e: std::io::Error| nomics_gpu::GpuError::HostMemory(e.to_string());
    let capacity = (reads + reads / 4).max(1);
    let mut seqs = nomics_mem::HostBuf::new((bases + bases / 4).max(1)).map_err(host)?;
    let mut offs = nomics_mem::HostBuf::new(capacity * 8).map_err(host)?;
    let mut lens = nomics_mem::HostBuf::new(capacity * 4).map_err(host)?;
    let mut out = nomics_gpu::GpuSeedBuffers::new(capacity, GPU_SEED_CAP)?;
    // Touch on the CPU and promote to 2 MiB pages; warns if that fails.
    seqs.touch_huge("GPU seed reads");
    offs.touch_huge("GPU seed offsets");
    lens.touch_huge("GPU seed lengths");
    out.touch_and_collapse();
    SEED_SETS.fetch_add(1, Ordering::Relaxed);
    Ok(SeedBufs {
        seqs: nomics_gpu::GpuHostBuf::register(seqs),
        offs: nomics_gpu::GpuHostBuf::register(offs),
        lens: nomics_gpu::GpuHostBuf::register(lens),
        out,
        retry: None,
        retried: Vec::new(),
    })
}

/// Two-pass mode's second pass: re-run the first pass's overflow reads
/// (`out` holds the first `offs.len()` reads' results) with the 128-entry
/// kernel into `retry`, growing it if needed. Returns the kernel time in
/// milliseconds; on error the overflow reads stay on the CPU (same output).
#[allow(clippy::too_many_arguments)]
fn seed_retry(
    bwt: &nomics_gpu::GpuBwt,
    table: Option<&nomics_gpu::GpuKmerTable<'_>>,
    seqs: &[u8],
    offs: &[u64],
    lens: &[u32],
    out: &nomics_gpu::GpuSeedBuffers,
    retry: &mut Option<SeedRetry>,
    retried: &mut Vec<u32>,
) -> Result<f32, nomics_gpu::GpuError> {
    retried.clear();
    retried.extend(
        (0..offs.len() as u32).filter(|&r| out.count(r as usize) == nomics_gpu::SEED_OVERFLOW),
    );
    if retried.is_empty() {
        return Ok(0.0);
    }
    let rows = seed_retry_rows(retried.len(), out.capacity());
    if retry
        .as_ref()
        .is_none_or(|r| r.out.capacity() < rows.min(retried.len()))
    {
        // The old set (if any) has been waited on, so dropping it is safe.
        *retry = None;
        *retry = Some(SeedRetry::new(rows).inspect_err(|_| retried.clear())?);
    }
    let r = retry.as_mut().expect("retry set allocated above");
    // Past the set's rows, overflow reads stay on the CPU.
    retried.truncate(r.out.capacity());
    let m = retried.len();
    let (roffs, rlens) = (
        &mut r.offs.as_u64s_mut()[..m],
        &mut r.lens.as_u32s_mut()[..m],
    );
    for (j, &read) in retried.iter().enumerate() {
        roffs[j] = offs[read as usize];
        rlens[j] = lens[read as usize];
    }
    let ms = r
        .out
        .collect(bwt, table, seqs, roffs, rlens, nomics_gpu::SeedWork::Cap128);
    if ms.is_err() {
        retried.clear();
    }
    ms
}

/// [`gpu_seeds`], or `None` (with a warning) when the GPU call fails, e.g.
/// when the GB10 cannot allocate because the page cache holds the free
/// memory. Reads then seed on the CPU, which gives the same intervals.
pub fn try_gpu_seeds(idx: &Index, reads: &[Read]) -> Option<GpuSeeds> {
    match gpu_seeds_result(idx, reads) {
        Ok(seeds) => Some(seeds),
        Err(e) => {
            nomics_gpu::diag_gpu_error("GPU seed search", &e);
            eprintln!(
                "WARNING: GPU seed search failed ({e:?}); this batch seeds on the CPU (same output, slower)"
            );
            GPU_FALLBACKS[0].fetch_add(1, Ordering::Relaxed);
            None
        }
    }
}

/// Run `mem_collect_intv` for all reads on the GPU. Panics if the GPU call
/// fails; [`try_gpu_seeds`] falls back to the CPU instead.
pub fn gpu_seeds(idx: &Index, reads: &[Read]) -> GpuSeeds {
    gpu_seeds_result(idx, reads).expect("GPU seed search failed")
}

fn gpu_seeds_result(idx: &Index, reads: &[Read]) -> Result<GpuSeeds, nomics_gpu::GpuError> {
    let bases: usize = reads.iter().map(|r| r.seq.len()).sum();
    let mut bufs = seed_bufs(reads.len(), bases)?;
    let n = reads.len();
    let SeedBufs {
        seqs,
        offs,
        lens,
        out,
        retry,
        retried,
    } = &mut bufs;
    let (seqs, offs, lens) = (
        seqs.as_mut_slice(),
        &mut offs.as_u64s_mut()[..n],
        &mut lens.as_u32s_mut()[..n],
    );
    let mut at = 0usize;
    for (r, (off, len)) in reads.iter().zip(offs.iter_mut().zip(lens.iter_mut())) {
        *off = at as u64;
        *len = r.seq.len() as u32;
        seqs[at..at + r.seq.len()].copy_from_slice(&r.seq);
        at += r.seq.len();
    }
    let bwt = nomics_gpu::GpuBwt {
        primary: idx.bwt.primary,
        l2: idx.bwt.l2,
        seq_len: idx.bwt.seq_len,
        words: idx.bwt.words(),
    };
    let table = seed_table(idx).map(|t| t.gpu_table());
    let ms = out.collect(
        &bwt,
        table.as_ref(),
        &seqs[..at],
        offs,
        lens,
        nomics_gpu::SeedWork::Cap64,
    );
    let ms = match ms {
        Ok(ms) => ms,
        Err(e) => {
            // The buffers are fine; keep them for the next batch.
            release_seed_bufs(bufs);
            return Err(e);
        }
    };
    let first_overflow = out.overflow_count();
    retried.clear();
    let retry_ms = match seed_mode() {
        SeedMode::Cap64 => 0.0,
        SeedMode::TwoPass => match seed_retry(
            &bwt,
            table.as_ref(),
            &seqs[..at],
            offs,
            lens,
            out,
            retry,
            retried,
        ) {
            Ok(ms) => ms,
            Err(e) => {
                nomics_gpu::diag_gpu_error("GPU seed retry pass", &e);
                // Once per process: a failing allocation would repeat per batch.
                static WARNED: std::sync::atomic::AtomicBool =
                    std::sync::atomic::AtomicBool::new(false);
                if !WARNED.swap(true, Ordering::Relaxed) {
                    eprintln!(
                        "WARNING: GPU seed retry pass failed ({e:?}); overflow reads of such \
                         batches seed on the CPU (same output, slower)"
                    );
                }
                0.0
            }
        },
    };
    let results = SeedResults {
        first: out,
        retried,
        retry: retry.as_ref().map(|r| &r.out),
    };
    if std::env::var_os("NOMICS_SEED_VERIFY").is_some() {
        verify_seeds(idx, reads, &results, offs);
    }
    let compacting = std::time::Instant::now();
    let seeds = GpuSeeds::compact(&results, n);
    if std::env::var_os("NOMICS_TIMING").is_some() {
        let fallback = seeds.spans.iter().filter(|s| s.0 == NO_SEEDS).count();
        eprintln!(
            "[timing] gpu seed kernel {:.0} ms ({ms:.0} + {retry_ms:.0} retry), {fallback} reads fall back to the CPU ({first_overflow} after cap {GPU_SEED_CAP}, {} retried at cap {GPU_SEED_RETRY_CAP}), compact {:.3}s to {:.2} GB",
            ms + retry_ms,
            retried.len(),
            compacting.elapsed().as_secs_f64(),
            seeds.bytes() as f64 / 1e9
        );
    }
    release_seed_bufs(bufs);
    Ok(seeds)
}

/// Diagnostic (`NOMICS_SEED_VERIFY=1`): compare every read's GPU intervals
/// with the CPU's and report each mismatch (read, seed block of 128, and
/// whether the intervals are another read's of the batch). Output unchanged.
fn verify_seeds(idx: &Index, reads: &[Read], out: &SeedResults<'_>, offs: &[u64]) {
    use crate::seed::{SmemAux, collect_intv};
    static BATCH: AtomicUsize = AtomicUsize::new(0);
    let batch = BATCH.fetch_add(1, Ordering::Relaxed);
    let as_intvs = |g: &[nomics_gpu::GpuIntv]| -> Vec<Intv> {
        g.iter().map(|g| Intv { x: g.x, info: g.info }).collect()
    };
    let cpu = |r: usize, aux: &mut SmemAux| {
        collect_intv(&idx.bwt, &reads[r].seq, aux);
        aux.mem.clone()
    };
    let blocks: Vec<usize> = (0..reads.len().div_ceil(512)).collect();
    let bad = std::sync::Mutex::new(Vec::new());
    let threads = std::thread::available_parallelism().map_or(8, |n| n.get());
    for_each_parallel(threads, blocks, |b| {
        let mut aux = SmemAux::default();
        for r in b * 512..((b + 1) * 512).min(reads.len()) {
            let Some(gpu) = out.get(r) else { continue };
            if as_intvs(gpu) != cpu(r, &mut aux) {
                bad.lock().unwrap().push(r);
            }
        }
    });
    let mut bad = bad.into_inner().unwrap();
    bad.sort_unstable();
    let mut aux = SmemAux::default();
    for &r in bad.iter().take(300) {
        let gpu = as_intvs(out.get(r).unwrap());
        // Whose intervals are these? Look at reads a whole number of seed
        // blocks away.
        let mut from = None;
        for d in (-(64 * 128) as isize..=64 * 128).step_by(128) {
            let s = r as isize + d;
            if d != 0 && s >= 0 && (s as usize) < reads.len() && cpu(s as usize, &mut aux) == gpu {
                from = Some(d);
                break;
            }
        }
        eprintln!(
            "SEED VERIFY batch {batch}: read {r} (block {}) differs: gpu {} intervals, cpu {}; offs {}; intervals of read {:?}",
            r / 128,
            gpu.len(),
            cpu(r, &mut aux).len(),
            offs[r],
            from.map(|d| r as isize + d)
        );
    }
    eprintln!(
        "SEED VERIFY batch {batch}: {} reads, {} differ",
        reads.len(),
        bad.len()
    );
}

/// Align a whole paired FASTQ stream in `-K` batches, appending SAM records
/// to `out`. With `gpu`, a producer thread reads batch N+1 and runs its seed
/// search on the GPU while the CPU threads align batch N; batch boundaries
/// and pair ids are unchanged, so the output is identical to the CPU path.
/// Returns the number of pairs and each batch's pestat log.
pub fn align_stream<R: BufRead + Send>(
    idx: &Index,
    rg_id: &[u8],
    r1: &mut R,
    r2: &mut R,
    threads: usize,
    gpu: bool,
    out: &mut Vec<u8>,
) -> io::Result<(u64, Vec<Vec<String>>)> {
    align_stream_to(idx, rg_id, r1, r2, threads, gpu, out, process_batch_seeded)
}

pub fn align_stream_bam<R: BufRead + Send>(
    idx: &Index,
    rg_id: &[u8],
    r1: &mut R,
    r2: &mut R,
    threads: usize,
    gpu: bool,
    out: &mut BamRecords,
) -> io::Result<(u64, Vec<Vec<String>>)> {
    align_stream_to(idx, rg_id, r1, r2, threads, gpu, out, process_batch_bam)
}

pub fn align_chunks<R: BufRead + Send>(
    idx: &Index,
    rg_id: &[u8],
    streams: Vec<(R, R)>,
    threads: usize,
    gpu: bool,
) -> io::Result<Vec<(u64, Vec<u8>, Vec<Vec<String>>)>> {
    let rg_ids = vec![rg_id.to_vec(); streams.len()];
    align_chunks_to(idx, &rg_ids, streams, threads, gpu, process_batch_seeded)
}

/// Align BAM-direct chunks and deliver each coordinate-sorted chunk as soon
/// as all of its `-K` batches have completed. Alignment remains limited to
/// [`batch_inflight`] FIFO batch jobs; one additional worker sorts completed
/// chunks in input order while later chunks continue aligning.
///
/// The return value retains the per-chunk pair counts and pestat logs without
/// retaining a second copy of the records. `rg_ids[chunk]` is the chunk's
/// read-group ID (lanes differ).
#[allow(clippy::type_complexity)]
pub fn align_chunks_bam<R, F>(
    idx: &Index,
    rg_ids: &[Vec<u8>],
    streams: Vec<(R, R)>,
    headers: Vec<Header>,
    threads: usize,
    gpu: bool,
    mut completed: F,
) -> io::Result<Vec<(u64, Vec<Vec<String>>)>>
where
    R: BufRead + Send,
    F: FnMut(usize, u64, RecordSet) -> io::Result<()> + Send,
{
    align_chunks_bam_with(
        idx,
        rg_ids,
        streams,
        headers,
        threads,
        gpu,
        false,
        move |chunk, pairs, header, records| {
            let set = sort_chunk(header, records)?;
            completed(chunk, pairs, set)
        },
    )
}

/// [`align_chunks_bam`] for whole genomes: each sorted chunk is written
/// straight to its BAM file (`output(chunk)` gives the path and the BGZF
/// threads) without ever holding the sorted chunk in memory, then
/// `completed(chunk, pairs, records)` runs. The file's bytes are those
/// [`RecordSet::write`] gives for the chunk [`align_chunks_bam`] delivers.
#[allow(clippy::type_complexity, clippy::too_many_arguments)]
pub fn align_chunks_bam_files<R, P, F>(
    idx: &Index,
    rg_ids: &[Vec<u8>],
    streams: Vec<(R, R)>,
    headers: Vec<Header>,
    threads: usize,
    gpu: bool,
    output: P,
    mut completed: F,
) -> io::Result<Vec<(u64, Vec<Vec<String>>)>>
where
    R: BufRead + Send,
    P: Fn(usize) -> (std::path::PathBuf, usize) + Send,
    F: FnMut(usize, u64, usize) -> io::Result<()> + Send,
{
    align_chunks_bam_with(
        idx,
        rg_ids,
        streams,
        headers,
        threads,
        gpu,
        true,
        move |chunk, pairs, header, records| {
            let (path, write_threads) = output(chunk);
            let n = crate::samsort::write_sorted_chunk(header, records, &path, write_threads)?;
            completed(chunk, pairs, n)
        },
    )
}

/// The chunk loop behind [`align_chunks_bam`] and [`align_chunks_bam_files`]:
/// `deliver(chunk, pairs, header, records)` runs on the sorter thread, in
/// chunk order, as each chunk completes. `whole_genome` sets the default
/// batch count ([`batch_inflight`]).
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
fn align_chunks_bam_with<R, F>(
    idx: &Index,
    rg_ids: &[Vec<u8>],
    streams: Vec<(R, R)>,
    headers: Vec<Header>,
    threads: usize,
    gpu: bool,
    whole_genome: bool,
    mut deliver: F,
) -> io::Result<Vec<(u64, Vec<Vec<String>>)>>
where
    R: BufRead + Send,
    F: FnMut(usize, u64, Header, ChunkRecords) -> io::Result<()> + Send,
{
    let count = streams.len();
    if headers.len() != count {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "alignment chunk/header count differs: {count} chunks, {} headers",
                headers.len()
            ),
        ));
    }
    let mut results: Vec<Option<(u64, Vec<Vec<String>>)>> = (0..count).map(|_| None).collect();
    std::thread::scope(|scope| -> io::Result<()> {
        let (sort_tx, sort_rx) = std::sync::mpsc::channel::<(usize, u64, ChunkRecords)>();
        let sorter = scope.spawn(move || -> io::Result<()> {
            let mut headers: Vec<Option<Header>> = headers.into_iter().map(Some).collect();
            for (chunk, pairs, records) in sort_rx {
                let header = headers[chunk]
                    .take()
                    .ok_or_else(|| io::Error::other("alignment chunk delivered twice"))?;
                let timing = std::env::var_os("NOMICS_TIMING").is_some();
                if timing {
                    mem_line(&format!("chunk {chunk} sort start"));
                }
                deliver(chunk, pairs, header, records)?;
                if timing {
                    mem_line(&format!("chunk {chunk} delivered"));
                }
            }
            Ok(())
        });
        let aligned = align_chunks_to_each(
            idx,
            rg_ids,
            streams,
            threads,
            gpu,
            whole_genome,
            process_batch_bam,
            |chunk, pairs, records, logs| {
                results[chunk] = Some((pairs, logs));
                sort_tx
                    .send((chunk, pairs, records))
                    .map_err(|_| io::Error::other("alignment chunk sorter stopped"))
            },
        );
        drop(sort_tx);
        let sorted = sorter
            .join()
            .map_err(|_| io::Error::other("alignment chunk sorter panicked"))?;
        // A sorter error (e.g. no space left writing a chunk BAM) is the
        // cause: the aligner only sees its send fail ("sorter stopped").
        sorted?;
        aligned
    })?;
    results
        .into_iter()
        .map(|result| result.ok_or_else(|| io::Error::other("missing alignment chunk result")))
        .collect()
}

/// What the reader hands the seeder: a batch as cut (chunk, first pair id,
/// reads), or the end of a chunk.
enum Feed {
    /// (chunk, first pair id, batch ordinal in the chunk, reads)
    Batch(io::Result<(usize, u64, usize, Vec<Read>)>),
    ChunkRead(usize),
}

type SeededBatch = io::Result<(usize, u64, usize, Vec<Read>, Option<GpuSeeds>)>;

/// Spawn the seeder between the reader and the aligners: it takes one batch
/// at a time from the reader (a rendezvous, so the reader reads the next
/// batch while this one seeds, and never blocks on the GPU itself), seeds it
/// on the GPU (with `gpu`; one seeding launch at a time, on the seeding
/// stream) and queues it, in order, for the aligners. After queuing a
/// chunk's last batch it marks the chunk read (`read_complete`), which the
/// consumer relies on: a chunk is complete only once all of its batches are
/// queued. Returns the reader's end of the rendezvous.
fn spawn_seeder<'scope, 'env>(
    scope: &'scope std::thread::Scope<'scope, 'env>,
    idx: &'env Index,
    gpu: bool,
    out: std::sync::mpsc::SyncSender<SeededBatch>,
    read_complete: Option<&'env AtomicUsize>,
) -> std::sync::mpsc::SyncSender<Feed> {
    let (feed, batches) = std::sync::mpsc::sync_channel::<Feed>(0);
    scope.spawn(move || {
        for item in batches {
            match item {
                Feed::Batch(batch) => {
                    let batch = batch.map(|(chunk, first_pair, ordinal, reads)| {
                        let seeds = gpu.then(|| try_gpu_seeds(idx, &reads)).flatten();
                        (chunk, first_pair, ordinal, reads, seeds)
                    });
                    let stop = batch.is_err();
                    if out.send(batch).is_err() || stop {
                        return;
                    }
                }
                Feed::ChunkRead(chunk) => {
                    if let Some(done) = read_complete {
                        done.store(chunk + 1, Ordering::Release);
                    }
                }
            }
        }
    });
    feed
}

/// `-K` batches aligned at once (`NOMICS_BATCH_INFLIGHT`, 1..=4, overrides).
/// Default 3 for the fused extension on whole genomes (`whole_genome`: chunk
/// BAMs streamed to files, [`align_chunks_bam_files`]), else 2. On the fused
/// path a third batch is faster (30x chunk 1: 234 -> 212 s; with 2, chunk
/// time is set by each batch's serial tail of pestat, rescue and pairing)
/// and fits the 90 GB budget since the rescue engines stopped mapping
/// buffers the wavefront kernel never uses (docs/strategy/
/// align-pipeline-design.md). Exomes do not gain from it (tier1 fused 306 s
/// at 2, 313-319 s at 3), and the CPU-extension path is CPU-bound.
fn batch_inflight(threads: usize, whole_genome: bool) -> usize {
    let default = if whole_genome && gpu_extend_mode() == ExtendMode::Fused {
        3
    } else {
        2
    };
    std::env::var("NOMICS_BATCH_INFLIGHT")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or(default)
        .clamp(1, 4)
        .min(threads.max(1))
}

/// The `-K` batches aligned at once ([`batch_inflight`]); `whole_genome` as
/// for [`gpu_warm_up`].
pub fn batches_in_flight(threads: usize, whole_genome: bool) -> usize {
    batch_inflight(threads, whole_genome)
}

/// Seeded batches queued for the aligners ahead of those in flight. One is
/// enough: the seeder holds the next batch (blocked on the queue) and the
/// reader the one after, so a finishing batch is replaced at once, and each
/// further queued batch held ~0.5 GB at 30x (reads and seeds); depth 3 was
/// no faster on the 30x chunk at 3 in flight (2026-10-23). Batches are cut,
/// numbered and consumed in order whatever the depth.
const SEEDED_QUEUE: usize = 1;

#[allow(clippy::type_complexity)]
fn align_chunks_to<R: BufRead + Send, O: BatchOutput>(
    idx: &Index,
    rg_ids: &[Vec<u8>],
    streams: Vec<(R, R)>,
    threads: usize,
    gpu: bool,
    process: fn(
        &Index,
        &[u8],
        &[Read],
        Option<GpuSeeds>,
        u64,
        usize,
        &mut O,
    ) -> ([PeStat; 4], Vec<String>),
) -> io::Result<Vec<(u64, O::Chunk, Vec<Vec<String>>)>> {
    let count = streams.len();
    let mut results: Vec<Option<(u64, O::Chunk, Vec<Vec<String>>)>> = (0..count).map(|_| None).collect();
    align_chunks_to_each(
        idx,
        rg_ids,
        streams,
        threads,
        gpu,
        false,
        process,
        |chunk, pairs, out, logs| {
            results[chunk] = Some((pairs, out, logs));
            Ok(())
        },
    )?;
    results
        .into_iter()
        .map(|result| result.ok_or_else(|| io::Error::other("missing alignment chunk result")))
        .collect()
}

#[allow(clippy::type_complexity, clippy::too_many_arguments)]
fn align_chunks_to_each<R: BufRead + Send, O: BatchOutput>(
    idx: &Index,
    rg_ids: &[Vec<u8>],
    streams: Vec<(R, R)>,
    threads: usize,
    gpu: bool,
    whole_genome: bool,
    process: fn(
        &Index,
        &[u8],
        &[Read],
        Option<GpuSeeds>,
        u64,
        usize,
        &mut O,
    ) -> ([PeStat; 4], Vec<String>),
    mut completed: impl FnMut(usize, u64, O::Chunk, Vec<Vec<String>>) -> io::Result<()>,
) -> io::Result<()> {
    if rg_ids.len() != streams.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "one read-group ID per alignment chunk",
        ));
    }
    let inflight = batch_inflight(threads, whole_genome);
    if inflight == 1 {
        for (chunk, (mut r1, mut r2)) in streams.into_iter().enumerate() {
            let mut out = O::Chunk::default();
            let mut logs = Vec::new();
            let pairs = align_stream_each(
                idx,
                &rg_ids[chunk],
                &mut r1,
                &mut r2,
                threads,
                gpu,
                process,
                |batch: O, log| {
                    O::absorb(&mut out, batch.finish());
                    logs.push(log);
                    Ok(())
                },
            )?;
            completed(chunk, pairs, out, logs)?;
        }
        return Ok(());
    }

    let count = streams.len();
    let mut outputs: Vec<Option<O::Chunk>> = (0..count).map(|_| Some(O::Chunk::default())).collect();
    let mut pairs = vec![0u64; count];
    let mut logs = vec![Vec::new(); count];
    let read_complete = AtomicUsize::new(0);
    std::thread::scope(|scope| -> io::Result<()> {
        let (tx, rx) = std::sync::mpsc::sync_channel::<
            io::Result<(usize, u64, usize, Vec<Read>, Option<GpuSeeds>)>,
        >(SEEDED_QUEUE);
        let read_complete = &read_complete;
        let feed = spawn_seeder(scope, idx, gpu, tx, Some(read_complete));
        scope.spawn(move || {
            // Streams are consumed in order and dropped at EOF, so a chunk's
            // input is freed as soon as it has been read. Pair ids count from
            // 0 in each chunk and are fixed here, where batches are cut.
            for (chunk, (mut r1, mut r2)) in streams.into_iter().enumerate() {
                // NOMICS_ALIGN_SKIP_BATCHES / NOMICS_ALIGN_MAX_BATCHES.
                let mut cutter = crate::diag::BatchCutter::new(crate::diag::batch_range());
                loop {
                    if !cutter.wants_more() {
                        break;
                    }
                    let batch = match read_batch(&mut r1, &mut r2) {
                        Ok(reads) if reads.is_empty() => break,
                        Ok(reads) => match cutter.cut((reads.len() / 2) as u64) {
                            crate::diag::Cut::Skip => continue,
                            crate::diag::Cut::Align {
                                first_pair,
                                ordinal,
                            } => Ok((chunk, first_pair, ordinal, reads)),
                        },
                        Err(error) => Err(error),
                    };
                    let stop = batch.is_err();
                    if feed.send(Feed::Batch(batch)).is_err() || stop {
                        return;
                    }
                }
                if feed.send(Feed::ChunkRead(chunk)).is_err() {
                    return;
                }
            }
        });

        // Every in-flight batch gets all threads (oversubscribed): while one
        // waits on the GPU, the other uses its cores. NOMICS_ALIGN_WORKERS
        // overrides.
        let workers = std::env::var("NOMICS_ALIGN_WORKERS")
            .ok()
            .and_then(|v| v.parse::<usize>().ok())
            .unwrap_or(threads.max(1));
        let mut pending = std::collections::VecDeque::new();
        let mut positions = crate::diag::PairCursor::new(count, crate::diag::batch_range().0);
        let mut finish =
            |handle: std::thread::ScopedJoinHandle<'_, (usize, u64, u64, O::Batch, Vec<String>)>,
             outputs: &mut [Option<O::Chunk>],
             pairs: &mut [u64],
             logs: &mut [Vec<Vec<String>>]|
             -> io::Result<()> {
                let (chunk, first_pair, count, out, log) = handle
                    .join()
                    .map_err(|_| io::Error::other("alignment batch thread panicked"))?;
                let in_order = positions.accept(chunk, first_pair, count);
                debug_assert!(in_order, "chunk {chunk}: batch at pair {first_pair} out of order");
                let output = outputs[chunk].as_mut().unwrap();
                O::absorb(output, out);
                if std::env::var_os("NOMICS_TIMING").is_some() {
                    let (used, held) = O::held(output);
                    mem_line(&format!(
                        "chunk {chunk} batch {first_pair} done, records {:.2} GB ({:.2} GB allocated)",
                        used as f64 / 1e9,
                        held as f64 / 1e9
                    ));
                }
                pairs[chunk] += count;
                logs[chunk].push(log);
                Ok(())
            };
        let mut next_complete = 0usize;
        let mut emit_ready = |pending: &std::collections::VecDeque<(
            usize,
            std::thread::ScopedJoinHandle<'_, (usize, u64, u64, O::Batch, Vec<String>)>,
        )>,
                              outputs: &mut [Option<O::Chunk>],
                              pairs: &mut [u64],
                              logs: &mut [Vec<Vec<String>>]|
         -> io::Result<()> {
            let read_through = read_complete.load(Ordering::Acquire);
            while next_complete < read_through
                && pending
                    .front()
                    .is_none_or(|(chunk, _)| *chunk > next_complete)
            {
                let out = outputs[next_complete].take().unwrap();
                let chunk_pairs = pairs[next_complete];
                let chunk_logs = std::mem::take(&mut logs[next_complete]);
                completed(next_complete, chunk_pairs, out, chunk_logs)?;
                next_complete += 1;
            }
            Ok(())
        };
        for batch in rx {
            let (chunk, first_pair, ordinal, reads, seeds) = batch?;
            pending.push_back((
                chunk,
                scope.spawn(move || {
                    crate::diag::set_chunk(chunk, ordinal, &rg_ids[chunk]);
                    let count = (reads.len() / 2) as u64;
                    let mut out = O::default();
                    let (_, log) = process(
                        idx,
                        &rg_ids[chunk],
                        &reads,
                        seeds,
                        first_pair,
                        workers,
                        &mut out,
                    );
                    drop(reads);
                    (chunk, first_pair, count, out.finish(), log)
                }),
            ));
            if pending.len() == inflight {
                finish(
                    pending.pop_front().unwrap().1,
                    &mut outputs,
                    &mut pairs,
                    &mut logs,
                )?;
                emit_ready(&pending, &mut outputs, &mut pairs, &mut logs)?;
            }
        }
        while let Some((_, handle)) = pending.pop_front() {
            finish(handle, &mut outputs, &mut pairs, &mut logs)?;
            emit_ready(&pending, &mut outputs, &mut pairs, &mut logs)?;
        }
        emit_ready(&pending, &mut outputs, &mut pairs, &mut logs)?;
        debug_assert_eq!(next_complete, count);
        Ok(())
    })
}

#[allow(clippy::type_complexity)]
fn align_stream_to<R: BufRead + Send, O: BatchOutput>(
    idx: &Index,
    rg_id: &[u8],
    r1: &mut R,
    r2: &mut R,
    threads: usize,
    gpu: bool,
    out: &mut O,
    process: fn(
        &Index,
        &[u8],
        &[Read],
        Option<GpuSeeds>,
        u64,
        usize,
        &mut O,
    ) -> ([PeStat; 4], Vec<String>),
) -> io::Result<(u64, Vec<Vec<String>>)> {
    let mut logs = Vec::new();
    let pairs = align_stream_each(idx, rg_id, r1, r2, threads, gpu, process, |batch, log| {
        out.append(batch);
        logs.push(log);
        Ok(())
    })?;
    Ok((pairs, logs))
}

/// [`align_stream`], but each `-K` batch's SAM records are handed to `emit`
/// (in input order, with the batch's pestat log) as soon as the batch is
/// aligned, so a caller can stream them out like `bwa mem` does.
pub fn align_stream_sam_each<R: BufRead + Send>(
    idx: &Index,
    rg_id: &[u8],
    r1: &mut R,
    r2: &mut R,
    threads: usize,
    emit: impl FnMut(Vec<u8>, Vec<String>) -> io::Result<()>,
) -> io::Result<u64> {
    align_stream_each(
        idx,
        rg_id,
        r1,
        r2,
        threads,
        false,
        process_batch_seeded,
        emit,
    )
}

#[allow(clippy::type_complexity, clippy::too_many_arguments)]
fn align_stream_each<R: BufRead + Send, O: BatchOutput>(
    idx: &Index,
    rg_id: &[u8],
    r1: &mut R,
    r2: &mut R,
    threads: usize,
    gpu: bool,
    process: fn(
        &Index,
        &[u8],
        &[Read],
        Option<GpuSeeds>,
        u64,
        usize,
        &mut O,
    ) -> ([PeStat; 4], Vec<String>),
    mut emit: impl FnMut(O, Vec<String>) -> io::Result<()>,
) -> io::Result<u64> {
    let inflight = batch_inflight(threads, false);
    let mut first_pair = 0u64;
    if inflight == 1 && !gpu {
        loop {
            let reads = read_batch(r1, r2)?;
            if reads.is_empty() {
                break;
            }
            let mut batch_out = O::default();
            let (_, log) = process(idx, rg_id, &reads, None, first_pair, threads, &mut batch_out);
            emit(batch_out, log)?;
            first_pair += (reads.len() / 2) as u64;
        }
        return Ok(first_pair);
    }
    std::thread::scope(|s| {
        let (tx, rx) = std::sync::mpsc::sync_channel::<
            io::Result<(usize, u64, usize, Vec<Read>, Option<GpuSeeds>)>,
        >(SEEDED_QUEUE);
        let feed = spawn_seeder(s, idx, gpu, tx, None);
        s.spawn(move || {
            let mut pair = 0u64;
            loop {
                let batch = read_batch(r1, r2);
                let batch = match batch {
                    Ok(reads) if reads.is_empty() => break,
                    Ok(reads) => {
                        let count = (reads.len() / 2) as u64;
                        let result = Ok((0, pair, 0, reads));
                        pair += count;
                        result
                    }
                    Err(e) => Err(e),
                };
                let stop = batch.is_err();
                if feed.send(Feed::Batch(batch)).is_err() || stop {
                    break;
                }
            }
        });
        // Every in-flight batch gets all threads (oversubscribed): while one
        // waits on the GPU, the other uses its cores. NOMICS_ALIGN_WORKERS
        // overrides.
        let workers = std::env::var("NOMICS_ALIGN_WORKERS")
            .ok()
            .and_then(|v| v.parse::<usize>().ok())
            .unwrap_or(threads.max(1));
        let mut pending = std::collections::VecDeque::new();
        for batch in rx {
            let (_, pair, _, reads, seeds) = batch?;
            pending.push_back(s.spawn(move || {
                let count = (reads.len() / 2) as u64;
                let mut batch_out = O::default();
                let (_, log) = process(
                    idx,
                    rg_id,
                    &reads,
                    seeds,
                    pair,
                    workers,
                    &mut batch_out,
                );
                (pair, count, batch_out, log)
            }));
            if pending.len() == inflight {
                let (pair, count, batch_out, log) = pending
                    .pop_front()
                    .unwrap()
                    .join()
                    .map_err(|_| io::Error::other("alignment batch thread panicked"))?;
                debug_assert_eq!(pair, first_pair);
                emit(batch_out, log)?;
                first_pair += count;
            }
        }
        while let Some(handle) = pending.pop_front() {
            let (pair, count, batch_out, log) = handle
                .join()
                .map_err(|_| io::Error::other("alignment batch thread panicked"))?;
            debug_assert_eq!(pair, first_pair);
            emit(batch_out, log)?;
            first_pair += count;
        }
        Ok(first_pair)
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::samsort::{header, sort_bam, sort_sam};
    use nomics_io::gzip::GzReader;
    use std::io::{BufReader, Cursor};

    #[test]
    fn seed_mode_parses() {
        for d in [SeedMode::Cap64, SeedMode::TwoPass] {
            assert_eq!(parse_seed_mode(None, d), d);
            assert_eq!(parse_seed_mode(Some("two-pass"), d), SeedMode::TwoPass);
            assert_eq!(parse_seed_mode(Some("cap64"), d), SeedMode::Cap64);
            assert_eq!(parse_seed_mode(Some("cap128"), d), d);
        }
        assert_eq!(DEFAULT_SEED_MODE, SeedMode::TwoPass);
    }

    #[test]
    fn seed_retry_rows_are_bounded() {
        // A 30x set (662k reads + 25%): ~4% overflow gets 25% headroom.
        let cap = 827_500;
        assert_eq!(seed_retry_rows(26_480, cap), 33_100);
        // Few overflows: the 1/64 floor; many: the 1/8 ceiling.
        assert_eq!(seed_retry_rows(10, cap), cap / 64);
        assert_eq!(seed_retry_rows(cap, cap), cap / 8);
        assert_eq!(seed_retry_rows(0, 1), 1);
    }

    /// Every split of a batch's rescue calls between GPU lanes and CPU
    /// workers (stealing or a fixed share, any tile) gives each call the CPU
    /// `ksw_align2` result at its own index. Without a GPU the lanes fall
    /// back to the CPU, which exercises the fallback path too.
    #[test]
    fn rescue_split_matches_cpu_whatever_the_split() {
        let mut state = 7u64;
        let mut next = |n: u64| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state % n
        };
        let mut seqs = Vec::new();
        for i in 0..1500 {
            let qlen = 100 + next(52) as usize;
            let tlen = [300, 800, 1334, 2000][i % 4] - next(40) as usize;
            let q: Vec<u8> = (0..qlen).map(|_| next(5) as u8).collect();
            let mut t: Vec<u8> = (0..tlen).map(|_| next(5) as u8).collect();
            let at = next((tlen - qlen) as u64) as usize;
            t[at..at + qlen].copy_from_slice(&q);
            seqs.push((q, t));
        }
        let inputs: Vec<_> = seqs
            .iter()
            .map(|(q, t)| nomics_gpu::GpuKswAlign2Input {
                query: q,
                target: t,
                m: 5,
                mat: &crate::extend::MAT,
                o_del: crate::extend::O_DEL,
                e_del: crate::extend::E_DEL,
                o_ins: crate::extend::O_INS,
                e_ins: crate::extend::E_INS,
                xtra: crate::ksw::KSW_XSUBO | crate::ksw::KSW_XSTART | crate::ksw::KSW_XBYTE | 19,
                size: 1,
            })
            .collect();
        let want: Vec<_> = inputs.iter().map(rescue_cpu).collect();
        let splits = [
            (None, None, None, None, 2),
            (Some(4), None, Some(100), Some(0.0), 2),
            (Some(4), None, None, None, 2),
            (Some(3), None, Some(7), Some(0.5), 3),
            (None, Some(0.0), Some(64), None, 1),
            (Some(2), Some(0.3), Some(200), None, 2),
            (None, Some(0.5), None, None, 2),
            (Some(5), Some(1.0), None, None, 2),
            (Some(8), None, Some(1), Some(1.0), 1),
        ];
        for (cpu_threads, cpu_share, gpu_tile, gpu_guided, gpu_lanes) in splits {
            let split = RescueSplit {
                cpu_threads,
                cpu_share,
                gpu_tile,
                gpu_guided,
                gpu_lanes,
            };
            let (got, stats) = rescue_split(&inputs, 6, split);
            assert!(got == want, "{split:?}");
            assert!(stats.gpu_jobs <= inputs.len());
        }
        assert!(rescue_split(&[], 4, RescueSplit::from_env()).0.is_empty());
    }

    /// Calls the GPU kernels refuse (byte-path gaps BWA truncates to
    /// uint8_t; 16-bit scores past BWA's saturation) come back from the
    /// split with the CPU's results.
    #[test]
    fn refused_rescue_scoring_falls_back_to_the_cpu() {
        let rep = |b: u8, n: usize| std::iter::repeat_n(b, n);
        let q1: Vec<u8> = rep(0, 60).chain(rep(1, 60)).collect();
        let t1: Vec<u8> = rep(0, 60).chain(rep(2, 20)).chain(rep(1, 60)).collect();
        let q2: Vec<u8> = rep(0, 259).collect();
        let mut mat2 = crate::extend::MAT;
        for i in 0..4 {
            mat2[i * 5 + i] = 127;
        }
        let call = |query, target, mat, o_del, xtra: i32| nomics_gpu::GpuKswAlign2Input {
            query,
            target,
            m: 5,
            mat,
            o_del,
            e_del: 1,
            o_ins: 6,
            e_ins: 1,
            xtra,
            size: if xtra & crate::ksw::KSW_XBYTE != 0 { 1 } else { 2 },
        };
        let inputs = [
            call(&q1, &t1, &crate::extend::MAT, 255, 0xd0013),
            call(&q2, &q2, &mat2, 6, 0xc0013),
        ];
        let want: Vec<_> = inputs.iter().map(rescue_cpu).collect();
        assert_eq!(want[1].score, 32767);
        assert!(rescue_split(&inputs, 2, RescueSplit::from_env()).0 == want);
    }
    use std::path::Path;

    /// Tiny real suffix-array/BWT index, built by sorting suffixes in the test.
    /// No external index or CUDA is needed to exercise tile offsets/fallback.
    fn extension_test_index() -> (Index, Vec<u8>, Vec<usize>) {
        test_index(17)
    }

    /// [`extension_test_index`] over a random reference drawn from `seed`.
    fn test_index(seed: u64) -> (Index, Vec<u8>, Vec<usize>) {
        let mut state = seed;
        let reference: Vec<_> = (0..1024)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                (state & 3) as u8
            })
            .collect();
        let text: Vec<_> = reference
            .iter()
            .copied()
            .chain(reference.iter().rev().map(|&b| 3 - b))
            .collect();
        let mut sa: Vec<_> = (0..=text.len()).collect();
        sa.sort_by(|&a, &b| text[a..].cmp(&text[b..]));
        let primary = sa.iter().position(|&s| s == 0).unwrap() as u64;
        let bwt: Vec<_> = sa
            .iter()
            .filter(|&&s| s > 0)
            .map(|&s| text[s - 1])
            .collect();
        let mut l2 = [0u64; 4];
        for &b in &text {
            for n in &mut l2[b as usize..] {
                *n += 1;
            }
        }
        let dir = std::env::temp_dir().join(format!(
            "nomics-extension-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir(&dir).unwrap();
        let words64 = |words: &[u64]| {
            words
                .iter()
                .flat_map(|n| n.to_le_bytes())
                .collect::<Vec<_>>()
        };
        let mut bwt_bytes = words64(&[primary, l2[0], l2[1], l2[2], l2[3]]);
        let mut counts = [0u64; 4];
        for block in bwt.chunks(128) {
            bwt_bytes.extend(words64(&counts));
            for bases in block.chunks(16) {
                let mut word = 0u32;
                for (i, &b) in bases.iter().enumerate() {
                    word |= (b as u32) << (30 - 2 * i);
                    counts[b as usize] += 1;
                }
                bwt_bytes.extend(word.to_le_bytes());
            }
        }
        bwt_bytes.extend(words64(&counts));
        std::fs::write(dir.join("ref.bwt"), bwt_bytes).unwrap();
        let mut sa_bytes = words64(&[primary, l2[0], l2[1], l2[2], l2[3], 1, text.len() as u64]);
        for &s in &sa[1..] {
            sa_bytes.extend((s as u64).to_le_bytes());
        }
        std::fs::write(dir.join("ref.sa"), sa_bytes).unwrap();
        let mut pac = vec![0u8; reference.len() / 4 + 1];
        for (i, &b) in reference.iter().enumerate() {
            pac[i / 4] |= b << (6 - 2 * (i % 4));
        }
        std::fs::write(dir.join("ref.pac"), pac).unwrap();
        std::fs::write(dir.join("ref.ann"), "1024 1 11\n0 chrTest\n0 1024 0\n").unwrap();
        std::fs::write(dir.join("ref.amb"), "1024 1 0\n").unwrap();
        let idx = Index::load(&dir.join("ref")).unwrap();
        std::fs::remove_dir_all(dir).unwrap();
        (idx, reference, sa)
    }

    /// Two indexes in one process: each gets its own k-mer table, and a cache
    /// file built for one reference is never used for the other.
    #[test]
    fn seed_table_is_per_index() {
        let (a, _, _) = test_index(17);
        let (b, _, _) = test_index(99);
        let dir = std::env::temp_dir().join(format!("nomics-seed-table-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let (path_a, path_b) = (dir.join("a.kmer"), dir.join("b.kmer"));
        let ta = seed_table_at(&a, &path_a, 6).expect("table for a");
        assert!(std::ptr::eq(seed_table_at(&a, &path_a, 6).unwrap(), ta));
        // The file at path_a is a's cache: b must not get a's table from it.
        assert!(seed_table_at(&b, &path_a, 6).is_none());
        let tb = seed_table_at(&b, &path_b, 6).expect("table for b");
        assert!(!std::ptr::eq(ta, tb));
        let words = |t: &crate::kmer_table::KmerTable| t.gpu_table().intervals.to_vec();
        assert_ne!(words(ta), words(tb));
        // Each table is the one built from its own BWT.
        let fresh_a = crate::kmer_table::KmerTable::load_or_build(&dir.join("a2.kmer"), &a.bwt, 6);
        let fresh_b = crate::kmer_table::KmerTable::load_or_build(&dir.join("b2.kmer"), &b.bwt, 6);
        assert_eq!(words(ta), words(&fresh_a.unwrap()));
        assert_eq!(words(tb), words(&fresh_b.unwrap()));
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn tiled_extension_matches_cpu_in_input_order() {
        let (idx, reference, sa) = extension_test_index();
        let n = 2 * (2 * BLOCK_PAIRS) + 7; // Two full blocks and a partial tail.
        let mut reads = Vec::new();
        let mut lists: Vec<Option<Vec<Intv>>> = vec![Some(Vec::new()); n];
        for i in 0..n {
            let start = (37 * i) % 400;
            let len = [80, 300, 400][i % 3];
            let mut seq = reference[start..start + len].to_vec();
            seq[60] = (seq[60] + 1) % 4;
            if i % 5 == 0 {
                seq[61] = 4;
            }
            reads.push(Read {
                name: i.to_string().into_bytes(),
                seq,
                qual: vec![b'I'; len],
            });
            // Empty reads' seed lists exercise zero-call blocks/ranges. Other
            // reads include no-left seeds and >256-base left/right fallbacks.
            if i % 11 != 0 {
                for qbeg in [0, 23, len - 25] {
                    let row = sa.iter().position(|&s| s == start + qbeg).unwrap();
                    lists[i].as_mut().unwrap().push(Intv {
                        x: [row as u64, 0, 1],
                        info: ((qbeg as u64) << 32) | (qbeg + 19) as u64,
                    });
                }
            }
        }
        let seeds = GpuSeeds::from_lists(&lists);
        for (i, l) in lists.iter().enumerate() {
            assert_eq!(seeds.get(i), l.as_deref());
        }
        let mut aux = AlignAux::default();
        let expected: Vec<_> = reads
            .iter()
            .enumerate()
            .map(|(i, r)| align1_core_from(&idx, &r.seq, seeds.get(i).unwrap(), &mut aux))
            .collect();
        assert!(expected.iter().filter(|r| !r.is_empty()).count() > n / 2);
        // The fused path also runs on CPU seeds.
        let expected_cpu: Vec<_> = reads
            .iter()
            .map(|r| align1_core(&idx, &r.seq, &mut aux))
            .collect();
        for (mode, gpu_seeds) in [
            (ExtendMode::TwoPass, true),
            (ExtendMode::Fused, true),
            (ExtendMode::Fused, false),
        ] {
            let (seeds, expected) = if gpu_seeds {
                (Some(&seeds), &expected)
            } else {
                (None, &expected_cpu)
            };
            let splits = [
                Split::default(),
                Split {
                    cpu_lanes: 1,
                    cpu_share: None,
                },
                Split {
                    cpu_lanes: 2,
                    cpu_share: Some(0.5),
                },
                Split {
                    cpu_lanes: 1,
                    cpu_share: Some(1.0),
                },
            ];
            for split in splits {
                for (threads, tile_blocks, lanes) in [(1, 1, 1), (3, 1, 2), (4, 2, 2)] {
                    let mut got = vec![Vec::new(); n];
                    let work = got
                        .chunks_mut(2 * BLOCK_PAIRS)
                        .zip(reads.chunks(2 * BLOCK_PAIRS))
                        .enumerate()
                        .collect();
                    align1_gpu_extend_tiled(
                        &idx,
                        &reads,
                        seeds,
                        threads,
                        work,
                        tile_blocks,
                        lanes,
                        mode,
                        split,
                    );
                    assert_eq!(
                        &got, expected,
                        "{mode:?}, GPU seeds {gpu_seeds}, threads={threads}, tile_blocks={tile_blocks}, lanes={lanes}, {split:?}"
                    );
                }
            }
            assert_eq!(
                align1_gpu_extend_tiled(&idx, &[], None, 0, vec![], 1, 2, mode, Split::default()),
                0.0
            );
        }
    }

    /// Reads with 80-90 bp indels make `chain2aln` retry with band 200. The
    /// fused stand-in's results replay with no misses, the band-1 keys are
    /// consumed, and the regions (and the tiled fused path) match the CPU.
    #[test]
    fn fused_band_retries_replay_exactly() {
        use crate::extend::{Extended, align1_chains, align1_extend, plan_fused};
        let (idx, reference, _) = extension_test_index();
        let mut state = 99u64;
        let mut random = |n: usize| -> Vec<u8> {
            (0..n)
                .map(|_| {
                    state ^= state << 13;
                    state ^= state >> 7;
                    state ^= state << 17;
                    (state & 3) as u8
                })
                .collect()
        };
        let mut reads = Vec::new();
        for i in 0..24usize {
            let start = 29 * i % 600;
            let gap = 80 + i % 11;
            let mut seq = reference[start..start + 100].to_vec();
            if i % 2 == 0 {
                // Deletion: the read skips `gap` reference bases.
                seq.extend(&reference[start + 100 + gap..start + 250 + gap]);
            } else {
                seq.extend(random(gap));
                seq.extend(&reference[start + 100..start + 250]);
            }
            reads.push(Read {
                name: i.to_string().into_bytes(),
                qual: vec![b'I'; seq.len()],
                seq,
            });
        }
        let p = crate::fused::params();
        let mut engine = crate::fused::FusedEngine::Cpu(Vec::new());
        let mut aux = AlignAux::default();
        let (mut lookups, mut engine_retries) = ([0u32; 3], 0);
        let mut expected = Vec::new();
        for r in &reads {
            let chains = align1_chains(&idx, &r.seq, None, &mut aux);
            let (mut plans, mut bases, mut rseq) = (Vec::new(), Vec::new(), Vec::new());
            plan_fused(&idx, &r.seq, &chains, &mut plans, &mut bases, &mut rseq);
            let group = plans.iter().map(|p| {
                let (lq, lt, rq) = (p.lq as usize, p.lt as usize, p.rq as usize);
                let b = &bases[p.at..];
                nomics_gpu::GpuFusedInput {
                    lq: &b[..lq],
                    lt: &b[lq..lq + lt],
                    rq: &b[lq + lt..lq + lt + rq],
                    rt: &b[lq + lt + rq..lq + lt + rq + p.rt as usize],
                    h0: p.h0,
                }
            });
            let (out, _) = engine.compute(&p, &[group], 2).unwrap();
            let mut keys = [Vec::new(), Vec::new()];
            let mut results = [Vec::new(), Vec::new()];
            for (plan, x) in plans.iter().zip(out) {
                engine_retries += u32::from(x.n_left == 2) + u32::from(x.n_right == 2);
                crate::fused::push_keys(
                    &p,
                    plan.chain,
                    plan.seed,
                    plan.h0,
                    x,
                    &mut keys,
                    &mut results,
                );
            }
            let pre = Extended::new(&chains, [&keys[0], &keys[1]], [&results[0], &results[1]]);
            let got = align1_extend(&idx, &r.seq, &chains, &mut aux, &pre);
            let want = align1_core(&idx, &r.seq, &mut aux);
            assert_eq!(got, want, "read {}", String::from_utf8_lossy(&r.name));
            expected.push(want);
            for (n, k) in lookups.iter_mut().zip(pre.lookups()) {
                *n += k;
            }
        }
        assert_eq!(lookups[1], 0, "misses: {lookups:?}");
        assert!(
            lookups[2] >= 12 && engine_retries >= 12,
            "{lookups:?} {engine_retries}"
        );
        for (threads, lanes) in [(1, 1), (4, 2)] {
            let mut got = vec![Vec::new(); reads.len()];
            let work = got
                .chunks_mut(2 * BLOCK_PAIRS)
                .zip(reads.chunks(2 * BLOCK_PAIRS))
                .enumerate()
                .collect();
            align1_gpu_extend_tiled(
                &idx,
                &reads,
                None,
                threads,
                work,
                1,
                lanes,
                ExtendMode::Fused,
                Split::default(),
            );
            assert_eq!(got, expected);
        }
    }

    fn fastq_prefix(path: &str, records: usize) -> Vec<u8> {
        let mut reader =
            BufReader::with_capacity(1 << 20, GzReader::open(Path::new(path)).unwrap());
        let mut bytes = Vec::new();
        let mut line = Vec::new();
        for _ in 0..records * 4 {
            line.clear();
            assert_ne!(reader.read_until(b'\n', &mut line).unwrap(), 0);
            bytes.extend_from_slice(&line);
        }
        bytes
    }

    /// A chunk accumulated by reference and sorted per reference
    /// ([`sort_chunk`]) gives exactly [`sort_bam`]'s records over the same
    /// batches appended whole, including ties broken by batch order.
    #[test]
    fn chunk_sort_matches_whole_sort() {
        let (idx, reference, _) = extension_test_index();
        let mut state = 7u64;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        let revcomp = |s: &[u8]| s.iter().rev().map(|&b| 3 - b).collect::<Vec<u8>>();
        let mut batch = |n: usize| -> Vec<Read> {
            let mut reads = Vec::new();
            for p in 0..n {
                let (a, b) = if p % 7 == 3 {
                    // Unmapped: random bases.
                    let r: Vec<u8> = (0..100).map(|_| (next() & 3) as u8).collect();
                    (r.clone(), revcomp(&r))
                } else {
                    // Few distinct positions, so keys tie within and across batches.
                    let at = (next() % 8) as usize * 60;
                    let s = reference[at..at + 100].to_vec();
                    let m = revcomp(&reference[at + 150..at + 250]);
                    (s, m)
                };
                for seq in [a, b] {
                    reads.push(Read {
                        name: format!("p{p}").into_bytes(),
                        qual: vec![b'I'; seq.len()],
                        seq,
                    });
                }
            }
            reads
        };
        let batches: Vec<Vec<Read>> = (0..4).map(|b| batch(1500 + 110 * b)).collect();
        let run = |reads: &[Read], first: u64| {
            let mut out = BamRecords::default();
            process_batch_bam(&idx, b"t", reads, None, first, 2, &mut out);
            out
        };
        let mut whole = BamRecords::default();
        let mut chunks = [ChunkRecords::default(), ChunkRecords::default()];
        let mut first = 0u64;
        for reads in &batches {
            whole.append(run(reads, first));
            for chunk in &mut chunks {
                chunk.absorb(run(reads, first).into_bins());
            }
            first += (reads.len() / 2) as u64;
        }
        let [chunk, streamed] = chunks;
        let hdr = header(&idx, "@RG\tID:t\tSM:t", "t");
        let expected = sort_bam(hdr.clone(), whole).unwrap();
        let got = sort_chunk(hdr.clone(), chunk).unwrap();
        assert!(expected.len() > 10_000);
        // Several BGZF blocks and writer batches (16 blocks at 2 threads).
        assert!(expected.bytes() > 40 * 65_536, "{} bytes", expected.bytes());
        assert_eq!(got.header, expected.header);
        assert_eq!(got.len(), expected.len());
        for i in 0..got.len() {
            assert_eq!(got.record(i), expected.record(i), "record {i}");
        }
        // Streamed straight to a file: the bytes RecordSet::write gives.
        let dir = std::env::temp_dir().join(format!("nomics-chunk-sort-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let (a, b) = (dir.join("set.bam"), dir.join("streamed.bam"));
        expected.write(&a, 3).unwrap();
        let n = crate::samsort::write_sorted_chunk(hdr, streamed, &b, 2).unwrap();
        assert_eq!(n, expected.len());
        assert_eq!(std::fs::read(&a).unwrap(), std::fs::read(&b).unwrap());
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    #[ignore = "requires NOMICS_BWA_INDEX, NOMICS_DENSE_SA, NOMICS_FASTQ1 and NOMICS_FASTQ2"]
    fn bam_chunk_callback_matches_collected_output() {
        let index = std::env::var("NOMICS_BWA_INDEX").unwrap();
        let dense = std::env::var("NOMICS_DENSE_SA").unwrap();
        let fq1 = fastq_prefix(&std::env::var("NOMICS_FASTQ1").unwrap(), 32);
        let fq2 = fastq_prefix(&std::env::var("NOMICS_FASTQ2").unwrap(), 32);
        let mut idx = Index::load(Path::new(&index)).unwrap();
        idx.attach_dense_sa(Path::new(&dense), 2).unwrap();
        let streams = || {
            vec![
                (Cursor::new(fq1.clone()), Cursor::new(fq2.clone())),
                (Cursor::new(fq1.clone()), Cursor::new(fq2.clone())),
            ]
        };
        let headers: Vec<_> = (0..2)
            .map(|chunk| {
                header(
                    &idx,
                    "@RG\tID:equivalence.L1\tSM:equivalence",
                    &format!("equivalence chunk {chunk}"),
                )
            })
            .collect();

        let rg_ids = vec![b"equivalence.L1".to_vec(); 2];
        let mut expected = Vec::new();
        for (chunk, (mut r1, mut r2)) in streams().into_iter().enumerate() {
            let mut records = BamRecords::default();
            let (pairs, logs) =
                align_stream_bam(&idx, &rg_ids[chunk], &mut r1, &mut r2, 2, false, &mut records)
                    .unwrap();
            expected.push((
                pairs,
                logs,
                sort_bam(headers[chunk].clone(), records).unwrap(),
            ));
        }

        let mut delivered = Vec::new();
        let metadata = align_chunks_bam(
            &idx,
            &rg_ids,
            streams(),
            headers,
            2,
            false,
            |chunk, pairs, set| {
                delivered.push((chunk, pairs, set));
                Ok(())
            },
        )
        .unwrap();

        assert_eq!(delivered.len(), expected.len());
        for (chunk, (delivered_chunk, delivered_pairs, delivered_set)) in
            delivered.into_iter().enumerate()
        {
            let (expected_pairs, expected_logs, expected_set) = &expected[chunk];
            assert_eq!(delivered_chunk, chunk);
            assert_eq!(delivered_pairs, *expected_pairs);
            assert_eq!(metadata[chunk], (*expected_pairs, expected_logs.clone()));
            assert_eq!(delivered_set.header, expected_set.header);
            assert_eq!(delivered_set.len(), expected_set.len());
            for record in 0..delivered_set.len() {
                assert_eq!(
                    delivered_set.record(record),
                    expected_set.record(record),
                    "chunk {chunk}, record {record}"
                );
            }
        }
    }

    #[test]
    #[ignore = "requires NOMICS_BWA_INDEX, NOMICS_DENSE_SA, NOMICS_FASTQ1 and NOMICS_FASTQ2"]
    fn bam_direct_matches_sam_text_on_fastq() {
        let index = std::env::var("NOMICS_BWA_INDEX").unwrap();
        let dense = std::env::var("NOMICS_DENSE_SA").unwrap();
        let fq1 = std::env::var("NOMICS_FASTQ1").unwrap();
        let fq2 = std::env::var("NOMICS_FASTQ2").unwrap();
        let mut idx = Index::load(Path::new(&index)).unwrap();
        idx.attach_dense_sa(Path::new(&dense), 20).unwrap();
        let open = |path: &str| {
            BufReader::with_capacity(1 << 20, GzReader::open(Path::new(path)).unwrap())
        };

        let (mut text_r1, mut text_r2) = (open(&fq1), open(&fq2));
        let mut sam = Vec::new();
        let (text_pairs, _) = align_stream(
            &idx,
            b"equivalence.L1",
            &mut text_r1,
            &mut text_r2,
            20,
            false,
            &mut sam,
        )
        .unwrap();
        let text = sort_sam(
            header(
                &idx,
                "@RG\tID:equivalence.L1\tSM:equivalence",
                "equivalence",
            ),
            &sam,
            20,
        )
        .unwrap();

        let (mut bam_r1, mut bam_r2) = (open(&fq1), open(&fq2));
        let mut records = BamRecords::default();
        let (bam_pairs, _) = align_stream_bam(
            &idx,
            b"equivalence.L1",
            &mut bam_r1,
            &mut bam_r2,
            20,
            false,
            &mut records,
        )
        .unwrap();
        let bam = sort_bam(
            header(
                &idx,
                "@RG\tID:equivalence.L1\tSM:equivalence",
                "equivalence",
            ),
            records,
        )
        .unwrap();

        assert_eq!(text_pairs, bam_pairs);
        assert_eq!(text.header, bam.header);
        assert_eq!(text.len(), bam.len());
        for i in 0..text.len() {
            assert_eq!(text.record(i), bam.record(i), "record {i}");
        }
    }
}
