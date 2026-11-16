//! `HaplotypeCallerEngine.callRegion` for one active region, VCF mode
//! (spec §3–§8): finalize → assemble → trim → read likelihoods → realign →
//! genotype each event → annotate (+ dbSNP) → VCF records.
//! Each stage has its own gate; this module only threads their boundaries.
//! End-to-end gate: `examples/hc_call.rs` against the golden VCF.
//!
//! Parallelism: `prepare` (finalize through realignment) is a pure function
//! of the region and runs on a worker pool; `Finisher::finish` runs on the
//! traversal thread in region order, because the genotyper keeps state and
//! QD's jitter draws from the RNG the downsampler also uses. GVCF mode
//! (`run_gvcf`) draws no QD jitter, so it also computes reference confidence
//! on the pool and keeps only genotyping, variant-span reference confidence
//! and block combining in order.

use std::collections::BTreeMap;
use std::io;
use std::sync::mpsc;

use rayon::prelude::*;
use std::time::{Duration, Instant};

use nomics_io::bam::{RecordSet, RecordView};

use crate::activity::{Region, RegionSink, Traversal, traverse};
use crate::annotate::{AnnotationCall, AnnotationGenotype, AnnotationRead, annotate_mode};
use crate::assembly::{AsmRead, assemble};
use crate::call::{CallRegionInput, EvidenceIndex, call_region};
use crate::dbsnp::{Dbsnp, DbsnpSlice};
use crate::finalize::{HcRead, finalize_region};
use crate::genotype::{GenotypingEngine, NO_CALL_ALLELE};
use crate::jrandom::JavaRandom;
use crate::likelihood::{self, LkRead};
use crate::reads::{Interval, ReadStream};
use crate::refconfidence::{RefSite, Site};
use crate::realign::realign_reads_to_best_haplotype;
use crate::trim::{TrimHap, TrimResult, trim_region};
use crate::vcf::{Genotype, InfoValue, Record};

/// `getAssemblyRegionReference(…, 500)` padding around the padded span.
const REFERENCE_PADDING: i64 = 500;
/// Reads per likelihood task (`NOMICS_LK_CHUNK` overrides, for tuning).
fn likelihood_chunk() -> usize {
    static C: std::sync::OnceLock<usize> = std::sync::OnceLock::new();
    *C.get_or_init(|| {
        std::env::var("NOMICS_LK_CHUNK")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(32)
    })
}

/// Run-wide, read-only inputs shared by the worker threads.
pub struct Context<'a> {
    pub contig_names: Vec<String>,
    pub contig_lens: Vec<i64>,
    /// Raw FASTA bases of a contig (any case).
    pub reference: &'a (dyn Fn(usize) -> &'a [u8] + Sync),
    pub dbsnp: Option<&'a Dbsnp>,
}

/// `CachingIndexedFastaSequenceFile` without case or IUPAC preservation.
fn reference_base(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b @ (b'A' | b'C' | b'G' | b'T' | b'N') => b,
        _ => b'N',
    }
}

/// A region through realignment: everything `callRegion` computes before
/// `assignGenotypeLikelihoods`. Free of RNG and cross-region state.
pub struct Prepared {
    contig: usize,
    span: (i64, i64),
    padded_span: (i64, i64),
    full_ref: Vec<u8>,
    haplotypes: Vec<TrimHap>,
    realigned: Vec<HcRead>,
    likelihoods: Vec<Vec<f64>>,
    filtered: Vec<HcRead>,
    /// dbSNP records starting in the genotyping span (all event starts).
    dbsnp: DbsnpSlice,
}

/// Wall time per stage, summed over threads: finalize, assemble, trim,
/// likelihoods, realign, call loop, annotate + dbSNP.
pub type Timing = [Duration; 7];

fn lap(timing: &mut Timing, stage: usize, t: &mut Instant) {
    let now = Instant::now();
    timing[stage] += now - *t;
    *t = now;
}

/// A region through trimming and the dbSNP prefetch: the inputs of its read
/// likelihoods. Owned, so the likelihoods can be computed on the GPU thread.
pub struct Staged {
    contig: usize,
    span: (i64, i64),
    padded_span: (i64, i64),
    full_ref: Vec<u8>,
    haplotypes: Vec<TrimHap>,
    reads: Vec<HcRead>,
    filtered: Vec<HcRead>,
    dbsnp: DbsnpSlice,
}

impl Staged {
    /// PairHMM cells (read bases x haplotype bases).
    pub fn cells(&self) -> u64 {
        let r: u64 = self.reads.iter().map(|r| r.bases.len() as u64).sum();
        let h: u64 = self.haplotypes.iter().map(|h| h.bases.len() as u64).sum();
        r * h
    }
}

/// Finalize → assemble → trim → dbSNP prefetch (`None`: no calls).
pub fn stage(
    ctx: &Context,
    bam: &RecordSet,
    region: &Region,
    members: &[u32],
    timing: &mut Timing,
) -> io::Result<Option<Staged>> {
    if !region.active || members.is_empty() {
        return Ok(None);
    }
    let c = region.contig;
    let (name, len) = (&ctx.contig_names[c], ctx.contig_lens[c]);
    let (start, end) = (region.start, region.end);

    let mut t = Instant::now();
    let mut reads = Vec::with_capacity(members.len());
    for &i in members {
        reads.push(HcRead::from_record(
            &RecordView::new(bam.record(i as usize))?,
            i,
        ));
    }
    let finalized = finalize_region(len, start, end, &reads);
    lap(timing, 0, &mut t);

    let padded = ((start - 100).max(1), (end + 100).min(len));
    let full_start = (padded.0 - REFERENCE_PADDING).max(1);
    let full_end = (padded.1 + REFERENCE_PADDING).min(len);
    let full_ref: Vec<u8> = (ctx.reference)(c)[(full_start - 1) as usize..full_end as usize]
        .iter()
        .map(|&b| reference_base(b))
        .collect();
    let asm_reads: Vec<AsmRead> = finalized
        .iter()
        .map(|r| AsmRead {
            name: r.name.clone(),
            flag: r.flag,
            start: r.start,
            cigar: r.cigar.clone(),
            bases: r.bases.clone(),
            quals: r.quals.clone(),
        })
        .collect();
    let assembled: Vec<TrimHap> = assemble(name, start, end, full_start, &full_ref, &asm_reads)
        .into_iter()
        .map(|h| TrimHap {
            bases: h.bases,
            is_ref: h.is_ref,
            alignment_start_hap_wrt_ref: i64::from(h.alignment_start_hap_wrt_ref),
            cigar: h.cigar,
            contig: h.contig,
            start: h.start,
            end: h.end,
        })
        .collect();
    lap(timing, 1, &mut t);

    // NOMICS_DEBUG_REGION=<contig>:<start> describes one region on stderr.
    let debug = std::env::var("NOMICS_DEBUG_REGION").is_ok_and(|v| v == format!("{name}:{start}"));
    if debug {
        eprintln!(
            "[region {name}:{start}-{end}] {} reads, {} finalized, {} haplotypes",
            members.len(),
            finalized.len(),
            assembled.len()
        );
        for h in &assembled {
            let cigar: String = h
                .cigar
                .iter()
                .map(|&c| format!("{}{}", c >> 4, b"MIDNSHP=X"[(c & 0xf) as usize] as char))
                .collect();
            eprintln!(
                "  hap ref={} {}-{} {cigar} {}",
                h.is_ref,
                h.start,
                h.end,
                String::from_utf8_lossy(&h.bases)
            );
        }
    }
    let trimmed = trim_region(
        name, start, end, len, full_start, &full_ref, &assembled, &finalized,
    );
    if debug {
        match &trimmed {
            TrimResult::Trimmed {
                span,
                padded_span,
                variation_present,
                haplotypes,
                reads,
                ..
            } => eprintln!(
                "  trimmed {span:?} padded {padded_span:?} variation={variation_present} {} haplotypes {} reads",
                haplotypes.len(),
                reads.len()
            ),
            _ => eprintln!("  not trimmed"),
        }
    }
    let TrimResult::Trimmed {
        span,
        padded_span,
        variation_present,
        haplotypes,
        reads: kept_reads,
        filtered,
    } = trimmed
    else {
        lap(timing, 2, &mut t);
        return Ok(None);
    };
    lap(timing, 2, &mut t);
    if !variation_present || kept_reads.is_empty() {
        return Ok(None);
    }
    let dbsnp = match ctx.dbsnp {
        Some(db) => db.slice(name.as_bytes(), span.0, span.1)?,
        None => DbsnpSlice::default(),
    };
    Ok(Some(Staged {
        contig: c,
        span,
        padded_span,
        full_ref,
        haplotypes,
        reads: kept_reads,
        filtered,
        dbsnp,
    }))
}

/// Read likelihoods on the CPU. Read chunks are rayon tasks: idle threads
/// steal chunks of whatever region is still running, which keeps drains short.
fn cpu_likelihoods(st: &Staged) -> likelihood::Likelihoods {
    let hap_bases: Vec<&[u8]> = st.haplotypes.iter().map(|h| h.bases.as_slice()).collect();
    let lk_reads: Vec<LkRead> = st
        .reads
        .iter()
        .map(|r| LkRead {
            bases: &r.bases,
            quals: &r.quals,
            mapq: r.mapq,
        })
        .collect();
    let parts: Vec<likelihood::Likelihoods> = lk_reads
        .par_chunks(likelihood_chunk())
        .map(|c| likelihood::compute(&hap_bases, c))
        .collect();
    likelihood::merge(parts, likelihood_chunk(), hap_bases.len())
}

/// Realignment, completing a region once its likelihoods are known.
pub fn complete(st: Staged, lk: likelihood::Likelihoods, timing: &mut Timing) -> Prepared {
    let mut t = Instant::now();
    let evidence: Vec<HcRead> = lk.kept.iter().map(|&i| st.reads[i].clone()).collect();
    let realigned =
        realign_reads_to_best_haplotype(&evidence, &st.haplotypes, st.padded_span.0, &lk.values);
    lap(timing, 4, &mut t);
    Prepared {
        contig: st.contig,
        span: st.span,
        padded_span: st.padded_span,
        full_ref: st.full_ref,
        haplotypes: st.haplotypes,
        realigned,
        likelihoods: lk.values,
        filtered: st.filtered,
        dbsnp: st.dbsnp,
    }
}

/// Finalize → assemble → trim → likelihoods (CPU) → realign for one region
/// (`None`: the region yields no calls).
pub fn prepare(
    ctx: &Context,
    bam: &RecordSet,
    region: &Region,
    members: &[u32],
    timing: &mut Timing,
) -> io::Result<Option<Prepared>> {
    let Some(st) = stage(ctx, bam, region, members, timing)? else {
        return Ok(None);
    };
    let t = Instant::now();
    let lk = cpu_likelihoods(&st);
    timing[3] += t.elapsed();
    Ok(Some(complete(st, lk, timing)))
}

/// A region whose likelihoods go to the GPU thread.
#[cfg_attr(not(feature = "cuda"), allow(dead_code))]
pub struct GpuJob {
    seq: u64,
    staged: Staged,
    cells: u64,
    timing: Timing,
}

/// Where regions send likelihood work for the GPU, with a backlog budget:
/// past it, a region computes on the CPU instead, which balances the two.
struct GpuRoute<'a> {
    tx: mpsc::Sender<GpuJob>,
    pending: &'a std::sync::atomic::AtomicU64,
    budget: u64,
}

impl GpuRoute<'_> {
    fn try_submit(&self, seq: u64, staged: Staged, timing: Timing) -> Result<(), Staged> {
        use std::sync::atomic::Ordering::AcqRel;
        let cells = staged.cells();
        #[cfg(feature = "cuda")]
        let fits = staged
            .reads
            .iter()
            .all(|r| (1..=GPU_MAX_READ_LEN).contains(&r.bases.len()));
        // The stub's maximum is zero: only the empty set of reads fits.
        #[cfg(not(feature = "cuda"))]
        let fits = staged.reads.is_empty();
        if !fits
            || self
                .pending
                .fetch_update(AcqRel, std::sync::atomic::Ordering::Acquire, |p| {
                    (p + cells <= self.budget).then_some(p + cells)
                })
                .is_err()
        {
            return Err(staged);
        }
        match self.tx.send(GpuJob {
            seq,
            staged,
            cells,
            timing,
        }) {
            Ok(()) => Ok(()),
            Err(e) => {
                self.pending.fetch_sub(cells, AcqRel);
                Err(e.0.staged)
            }
        }
    }
}

#[cfg(feature = "cuda")]
const GPU_MAX_READ_LEN: usize = nomics_gpu::PAIRHMM_MAX_READ_LEN;

/// The GPU thread: batches likelihood jobs across regions, computes the
/// PairHMM sums on the GPU (bit-identical to the CPU lanes; see
/// `examples/pairhmm_bench.rs`), finishes them exactly as the CPU path does,
/// and hands each region's realignment back to the pool.
#[cfg(feature = "cuda")]
fn gpu_loop(
    rx: mpsc::Receiver<GpuJob>,
    pending: &std::sync::atomic::AtomicU64,
    pool: &rayon::ThreadPool,
    done: mpsc::Sender<Msg>,
    batch_pairs: usize,
) {
    use crate::pairhmm::{initial_condition, row_params};
    use crate::readprep::hmm_qualities;
    use nomics_gpu::{GpuPairHmm, GpuPairHmmPair};
    let (init, init_log10) = initial_condition();
    let mut engine = GpuPairHmm::new();
    let (mut read_bases, mut hap_bases, mut params, mut pairs, mut sums) =
        (Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new());
    let pairs_of = |j: &GpuJob| j.staged.reads.len() * j.staged.haplotypes.len();
    while let Ok(first) = rx.recv() {
        let mut n = pairs_of(&first);
        let mut jobs = vec![first];
        while n < batch_pairs {
            match rx.try_recv() {
                Ok(j) => {
                    n += pairs_of(&j);
                    jobs.push(j);
                }
                Err(_) => break,
            }
        }
        let t = Instant::now();
        read_bases.clear();
        hap_bases.clear();
        params.clear();
        pairs.clear();
        for j in &jobs {
            let hap_offsets: Vec<u64> = j
                .staged
                .haplotypes
                .iter()
                .map(|h| {
                    let o = hap_bases.len() as u64;
                    hap_bases.extend_from_slice(&h.bases);
                    o
                })
                .collect();
            for r in &j.staged.reads {
                let (q, ins, del) = hmm_qualities(&r.bases, &r.quals, r.mapq);
                let (read_offset, params_offset) = (read_bases.len() as u64, params.len() as u64);
                read_bases.extend_from_slice(&r.bases);
                params.extend((0..r.bases.len()).map(|i| row_params(ins[i], del[i], q[i], 10)));
                for (h, &hap_offset) in j.staged.haplotypes.iter().zip(&hap_offsets) {
                    pairs.push(GpuPairHmmPair {
                        read_offset,
                        hap_offset,
                        params_offset,
                        read_len: r.bases.len() as u32,
                        hap_len: h.bases.len() as u32,
                        initial: init / h.bases.len() as f64,
                    });
                }
            }
        }
        sums.clear();
        sums.resize(pairs.len(), 0.0);
        let ok = engine
            .compute(&read_bases, &hap_bases, &params, &pairs, &mut sums)
            .is_ok();
        let gpu_time = t.elapsed();
        let mut k = 0usize;
        for mut j in jobs {
            let (nr, nh) = (j.staged.reads.len(), j.staged.haplotypes.len());
            let share = gpu_time.mul_f64((nr * nh) as f64 / pairs.len().max(1) as f64);
            j.timing[3] += share;
            let lk = ok.then(|| {
                let mut values = vec![Vec::with_capacity(nr); nh];
                for _ in 0..nr {
                    for v in values.iter_mut() {
                        v.push(crate::jmath::log10(sums[k]) - init_log10);
                        k += 1;
                    }
                }
                let lens: Vec<usize> = j.staged.reads.iter().map(|r| r.bases.len()).collect();
                likelihood::finish(values, &lens)
            });
            pending.fetch_sub(j.cells, std::sync::atomic::Ordering::AcqRel);
            let done = done.clone();
            pool.spawn(move || {
                let GpuJob {
                    seq,
                    staged,
                    mut timing,
                    ..
                } = j;
                // A failed launch falls back to the CPU for the whole batch.
                let lk = lk.unwrap_or_else(|| cpu_likelihoods(&staged));
                let p = complete(staged, lk, &mut timing);
                progress::PREPARED.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                let _ = done.send(Msg::Done(seq, Ok(Some(p)), timing));
            });
        }
    }
}

/// The ordered tail of `callRegion`: the genotyper (whose span-deletion
/// memory carries across regions), annotations (QD draws from the shared
/// RNG) and dbSNP.
pub struct Finisher {
    engine: GenotypingEngine,
}

impl Default for Finisher {
    fn default() -> Self {
        Self::new()
    }
}

impl Finisher {
    pub fn new() -> Self {
        Self {
            engine: GenotypingEngine::new(),
        }
    }

    pub fn gvcf() -> Self {
        Self {
            engine: GenotypingEngine::gvcf(),
        }
    }

    pub fn finish(
        &mut self,
        ctx: &Context,
        p: Prepared,
        rng: &mut JavaRandom,
        timing: &mut Timing,
    ) -> io::Result<Vec<Record>> {
        let mut t = Instant::now();
        let c = p.contig;
        let name = &ctx.contig_names[c];
        let reference_hap = p
            .haplotypes
            .iter()
            .find(|h| h.is_ref)
            .expect("reference haplotype");
        // Haplotype scores only break ties when > 44 merged alleles must be
        // pruned; trimming does not carry them, so they are left at 0.
        let scores = vec![0.0; p.haplotypes.len()];
        let events = call_region(
            &mut self.engine,
            CallRegionInput {
                contig: c,
                active_span: p.span,
                padded_reference_start: p.padded_span.0 - reference_hap.alignment_start_hap_wrt_ref,
                padded_reference: &p.full_ref,
                haplotypes: &p.haplotypes,
                haplotype_scores: &scores,
                reads: &p.realigned,
                likelihoods: &p.likelihoods,
                filtered_reads: &p.filtered,
            },
        );
        lap(timing, 5, &mut t);

        let called_indices: std::collections::BTreeSet<usize> = events
            .iter()
            .flat_map(|e| e.called_haplotypes.iter().copied())
            .collect();
        let mut out = Vec::with_capacity(events.len());
        let debug_site = std::env::var("NOMICS_DEBUG_SITE").ok();
        for ev in events {
            let call = &ev.call;
            // NOMICS_DEBUG_SITE=<contig>:<pos> lists the reads annotating that call.
            if debug_site.as_deref() == Some(&format!("{name}:{}", ev.start)) {
                for row in &ev.annotation_rows {
                    let (kind, r) = match row.evidence {
                        EvidenceIndex::Read(i) => ("read", &p.realigned[i]),
                        EvidenceIndex::Filtered(i) => ("filtered", &p.filtered[i]),
                    };
                    let cigar: String = r
                        .cigar
                        .iter()
                        .map(|&c| format!("{}{}", c >> 4, b"MIDNSHP=X"[(c & 0xf) as usize] as char))
                        .collect();
                    eprintln!(
                        "[site {name}:{}] {kind} {} flag={} start={} end={} cigar={cigar} mapq={} len={}",
                        ev.start,
                        String::from_utf8_lossy(&r.name),
                        r.flag,
                        r.start,
                        r.end(),
                        r.mapq,
                        r.bases.len()
                    );
                }
            }
            assert!(
                self.engine.emit_gvcf || !call.gt.contains(&NO_CALL_ALLELE),
                "emitted call without a genotype"
            );
            let rows = ev
                .annotation_rows
                .iter()
                .map(|row| {
                    let r = match row.evidence {
                        EvidenceIndex::Read(i) => &p.realigned[i],
                        EvidenceIndex::Filtered(i) => &p.filtered[i],
                    };
                    AnnotationRead {
                        flag: r.flag,
                        mapq: r.mapq,
                        reverse: r.flag & 0x10 != 0,
                        start: r.start,
                        cigar: r.cigar.clone(),
                        bases: r.bases.clone(),
                        quals: r.quals.clone(),
                        likelihoods: row.likelihoods.clone(),
                    }
                })
                .collect();
            let ann = annotate_mode(
                &AnnotationCall {
                    alleles: ev.call_alleles.clone(),
                    matrix_alleles: ev.annotation_alleles.clone(),
                    log10_p_error: call.log10_p_error,
                    genotype: AnnotationGenotype {
                        alleles: call.gt,
                        pl: call.pl.clone(),
                        gq: call.gq,
                        mle_ac: call.mleac.iter().map(|&x| i64::from(x)).collect(),
                        mle_af: call.mleaf.clone(),
                    },
                    contig: name.clone(),
                    pos: ev.start,
                    end: ev.end,
                    reads: rows,
                },
                rng,
                self.engine.emit_gvcf,
            );
            let id = match ctx.dbsnp {
                Some(_) => {
                    let alleles: Vec<&[u8]> = ev.call_alleles.iter().map(Vec::as_slice).collect();
                    p.dbsnp.id_for(ev.start, &alleles)
                }
                None => None,
            };
            let mut info = ann.info;
            if id.is_some() {
                info.push(("DB".to_string(), InfoValue::Flag));
            }
            let text = |a: &Vec<u8>| String::from_utf8_lossy(a).into_owned();
            out.push(Record {
                contig: c,
                pos: ev.final_start,
                id,
                ref_allele: text(&ev.final_alleles[0]),
                alts: ev.final_alleles[1..].iter().map(text).collect(),
                qual: -10.0 * call.log10_p_error + 0.0,
                info,
                genotype: Genotype {
                    alleles: call.gt.to_vec(),
                    ad: ann.ad,
                    dp: ann.dp,
                    gq: call.gq,
                    pl: call.pl.clone(),
                    gvcf: self.engine.emit_gvcf.then_some(crate::vcf::GvcfGenotype {
                        sb: ann.sb,
                        ..Default::default()
                    }),
                },
            });
        }
        if self.engine.emit_gvcf {
            crate::gvcf::phase_calls(
                &mut out,
                &p.haplotypes
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| called_indices.contains(i))
                    .map(|(_, h)| h.clone())
                    .collect::<Vec<_>>(),
                &p.full_ref,
                p.padded_span.0 - reference_hap.alignment_start_hap_wrt_ref,
            );
        }
        lap(timing, 6, &mut t);
        Ok(out)
    }
}

/// Progress counters for the `NOMICS_TIMING` utilization timeline.
pub(crate) mod progress {
    use std::sync::atomic::{AtomicBool, AtomicU64};
    pub static SENT: AtomicU64 = AtomicU64::new(0);
    pub static PREPARED: AtomicU64 = AtomicU64::new(0);
    pub static FINISHED: AtomicU64 = AtomicU64::new(0);
    pub static IN_DRAIN: AtomicBool = AtomicBool::new(false);
    pub static IN_PRECOMPUTE: AtomicBool = AtomicBool::new(false);
}

/// Process CPU seconds (user + system) from `/proc/self/stat`.
fn cpu_seconds() -> Option<f64> {
    let stat = std::fs::read_to_string("/proc/self/stat").ok()?;
    let rest = &stat[stat.rfind(')')? + 2..];
    let f: Vec<&str> = rest.split(' ').collect();
    let ticks: f64 = f.get(11)?.parse::<f64>().ok()? + f.get(12)?.parse::<f64>().ok()?;
    Some(ticks / 100.0)
}

/// Messages to the finisher thread, from the workers and the traversal.
enum Msg {
    /// A prepared region, by dispatch order.
    Done(u64, io::Result<Option<Prepared>>, Timing),
    /// A drain: finish every region below this sequence number, then send
    /// the RNG back.
    Sync(u64),
    /// The RNG after the downsampler's draws; resume finishing.
    Rng(JavaRandom),
}

/// The traversal side: dispatches regions to the rayon pool and hands the
/// RNG to the finisher thread between downsampler draws.
struct Pipeline<'s, 'c, 'a> {
    ctx: &'c Context<'a>,
    bam: &'c RecordSet,
    scope: &'s rayon::Scope<'c>,
    to_finisher: mpsc::Sender<Msg>,
    rng_back: mpsc::Receiver<JavaRandom>,
    sent: u64,
    /// Time the traversal spent blocked in `before_draw`.
    drain_wait: Duration,
    gpu: Option<GpuRoute<'c>>,
}

impl<'c> RegionSink for Pipeline<'_, 'c, '_> {
    fn on_call(&mut self, region: &Region, reads: &[u32], _rng: &mut JavaRandom) {
        if !region.active || reads.is_empty() {
            return;
        }
        let (ctx, bam, tx) = (self.ctx, self.bam, self.to_finisher.clone());
        let (seq, region, members) = (self.sent, *region, reads.to_vec());
        let gpu = self.gpu.as_ref().map(|g| GpuRoute {
            tx: g.tx.clone(),
            pending: g.pending,
            budget: g.budget,
        });
        self.scope.spawn(move |_| {
            let mut timing = Timing::default();
            let r = match stage(ctx, bam, &region, &members, &mut timing) {
                Ok(Some(st)) => {
                    let st = match &gpu {
                        Some(g) => match g.try_submit(seq, st, timing) {
                            Ok(()) => return, // the GPU thread completes it
                            Err(st) => st,
                        },
                        None => st,
                    };
                    let t = Instant::now();
                    let lk = cpu_likelihoods(&st);
                    timing[3] += t.elapsed();
                    Ok(Some(complete(st, lk, &mut timing)))
                }
                other => other.map(|_| None),
            };
            progress::PREPARED.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            let _ = tx.send(Msg::Done(seq, r, timing));
        });
        self.sent += 1;
        progress::SENT.store(self.sent, std::sync::atomic::Ordering::Relaxed);
    }

    fn before_draw(&mut self, rng: &mut JavaRandom) {
        let t = Instant::now();
        progress::IN_DRAIN.store(true, std::sync::atomic::Ordering::Relaxed);
        self.to_finisher
            .send(Msg::Sync(self.sent))
            .expect("finisher alive");
        *rng = self.rng_back.recv().expect("finisher alive");
        progress::IN_DRAIN.store(false, std::sync::atomic::Ordering::Relaxed);
        self.drain_wait += t.elapsed();
    }

    fn after_draw(&mut self, rng: &mut JavaRandom) {
        self.to_finisher
            .send(Msg::Rng(rng.clone()))
            .expect("finisher alive");
    }
}

/// The finisher thread: finishes prepared regions strictly in dispatch
/// order, holding the RNG except while the downsampler draws.
fn finish_loop(
    ctx: &Context,
    mut finisher: Finisher,
    rx: mpsc::Receiver<Msg>,
    rng_back: mpsc::Sender<JavaRandom>,
) -> (io::Result<Vec<Record>>, Timing) {
    let mut rng = Some(JavaRandom::gatk());
    let mut ready: BTreeMap<u64, io::Result<Option<Prepared>>> = BTreeMap::new();
    let (mut next, mut sync): (u64, Option<u64>) = (0, None);
    let (mut records, mut timing, mut error) = (Vec::new(), Timing::default(), None);
    loop {
        // Finish everything that is next in order.
        if let Some(r) = rng.as_mut() {
            while let Some(p) = ready.remove(&next) {
                next += 1;
                progress::FINISHED.store(next, std::sync::atomic::Ordering::Relaxed);
                let out = p.and_then(|p| match p {
                    Some(p) => finisher.finish(ctx, p, r, &mut timing),
                    None => Ok(Vec::new()),
                });
                match out {
                    Ok(recs) => records.extend(recs),
                    Err(e) => {
                        error.get_or_insert(e);
                    }
                }
            }
        }
        if let Some(upto) = sync
            && next >= upto
        {
            sync = None;
            let _ = rng_back.send(rng.take().expect("finisher holds the RNG"));
        }
        match rx.recv() {
            Ok(Msg::Done(seq, r, t)) => {
                for (a, b) in timing.iter_mut().zip(t) {
                    *a += b;
                }
                ready.insert(seq, r);
            }
            Ok(Msg::Sync(upto)) => sync = Some(upto),
            Ok(Msg::Rng(r)) => rng = Some(r),
            Err(_) => break,
        }
    }
    let result = match error {
        Some(e) => Err(e),
        None => Ok(records),
    };
    (result, timing)
}

/// The whole HaplotypeCaller run: `activity::traverse` on this thread,
/// region preparation on a `threads`-wide rayon pool, the ordered tail on a
/// finisher thread; VCF records in output order.
#[allow(clippy::too_many_arguments)]
pub fn run(
    ctx: &Context,
    bam: &RecordSet,
    stream: &ReadStream,
    unpadded: &[Vec<Interval>],
    reference: &(dyn Fn(usize, i64) -> u8 + Sync),
    threads: usize,
) -> io::Result<(Traversal, Vec<Record>, Timing)> {
    // NOMICS_SEQ_ACTIVITY=1 computes every pileup on the traversal thread.
    let activity_threads = if std::env::var_os("NOMICS_SEQ_ACTIVITY").is_some() {
        1
    } else {
        threads
    };
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads.max(1))
        .stack_size(64 << 20)
        .build()
        .map_err(io::Error::other)?;
    let (to_finisher, rx) = mpsc::channel::<Msg>();
    let (rng_tx, rng_back) = mpsc::channel::<JavaRandom>();
    let stop = std::sync::atomic::AtomicBool::new(false);
    let gpu_pending = std::sync::atomic::AtomicU64::new(0);
    let env_num = |k: &str, d: u64| {
        std::env::var(k)
            .ok()
            .and_then(|v| v.parse::<f64>().ok())
            .map_or(d, |v| v as u64)
    };
    // NOMICS_GPU=0 disables the hybrid PairHMM; NOMICS_GPU_BUDGET is the GPU
    // backlog in PairHMM cells beyond which regions stay on the CPU.
    let gpu_enabled = cfg!(feature = "cuda") && std::env::var("NOMICS_GPU").as_deref() != Ok("0");
    let gpu_budget = env_num("NOMICS_GPU_BUDGET", 2_000_000_000);
    let gpu_batch = env_num("NOMICS_GPU_BATCH", 65_536) as usize;
    std::thread::scope(|threads_scope| {
        // Deep recursion in the call loop/annotations wants more than the
        // 2 MiB default (the main thread has 8 MiB).
        let finisher = std::thread::Builder::new()
            .name("hc-finisher".into())
            .stack_size(256 << 20)
            .spawn_scoped(threads_scope, move || {
                finish_loop(ctx, Finisher::new(), rx, rng_tx)
            })?;
        if std::env::var_os("NOMICS_TIMING").is_some() {
            let stop = &stop;
            threads_scope.spawn(move || {
                use std::sync::atomic::Ordering::Relaxed;
                let t0 = Instant::now();
                let mut last = cpu_seconds().unwrap_or(0.0);
                while !stop.load(Relaxed) {
                    std::thread::sleep(Duration::from_secs(1));
                    let now = cpu_seconds().unwrap_or(0.0);
                    eprintln!(
                        "t={:>3.0}s cores={:>5.1} sent={} prepared={} finished={}{}{}",
                        t0.elapsed().as_secs_f64(),
                        now - last,
                        progress::SENT.load(Relaxed),
                        progress::PREPARED.load(Relaxed),
                        progress::FINISHED.load(Relaxed),
                        if progress::IN_DRAIN.load(Relaxed) {
                            " DRAIN"
                        } else {
                            ""
                        },
                        if progress::IN_PRECOMPUTE.load(Relaxed) {
                            " PRECOMPUTE"
                        } else {
                            ""
                        },
                    );
                    last = now;
                }
            });
        }
        let (gpu_tx, gpu_rx) = mpsc::channel::<GpuJob>();
        #[cfg(feature = "cuda")]
        if gpu_enabled {
            let (pool, pending, done) = (&pool, &gpu_pending, to_finisher.clone());
            std::thread::Builder::new()
                .name("hc-gpu".into())
                .spawn_scoped(threads_scope, move || {
                    gpu_loop(gpu_rx, pending, pool, done, gpu_batch)
                })?;
        }
        #[cfg(not(feature = "cuda"))]
        let _ = (gpu_rx, gpu_batch);
        let gpu = gpu_enabled.then(|| GpuRoute {
            tx: gpu_tx,
            pending: &gpu_pending,
            budget: gpu_budget,
        });
        let (tr, drain_wait) = pool.in_place_scope(|scope| {
            let mut pipe = Pipeline {
                gpu,
                ctx,
                bam,
                scope,
                to_finisher,
                rng_back,
                sent: 0,
                drain_wait: Duration::ZERO,
            };
            let mut rng = JavaRandom::gatk();
            let tr = traverse(
                bam,
                stream,
                unpadded,
                reference,
                &ctx.contig_lens,
                &mut rng,
                &mut pipe,
                activity_threads,
            );
            (tr, pipe.drain_wait)
        });
        // The pool scope has joined every task, and the Pipeline (with its
        // sender) is gone, so the finisher drains and exits.
        let (records, timing) = finisher.join().expect("finisher panicked");
        stop.store(true, std::sync::atomic::Ordering::Relaxed);
        if std::env::var_os("NOMICS_TIMING").is_some() {
            eprintln!(
                "traversal blocked in drains: {:.1}s",
                drain_wait.as_secs_f64()
            );
        }
        Ok((tr?, records?, timing))
    })
}

/// The header of a HaplotypeCaller VCF-mode run (default annotations,
/// dbSNP), minus the run-specific `##GATKCommandLine`, which the caller
/// supplies (normalized by `docs/GATES.md`).
pub fn header_meta(command_line: &str) -> Vec<String> {
    let mut m: Vec<String> = [
        r#"FILTER=<ID=LowQual,Description="Low quality">"#,
        r#"FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">"#,
        r#"FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">"#,
        r#"FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">"#,
        r#"FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        r#"FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">"#,
        r#"INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">"#,
        r#"INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">"#,
        r#"INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">"#,
        r#"INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities">"#,
        r#"INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP Membership">"#,
        r#"INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">"#,
        r#"INFO=<ID=ExcessHet,Number=1,Type=Float,Description="Phred-scaled p-value for exact test of excess heterozygosity">"#,
        r#"INFO=<ID=FS,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher's exact test to detect strand bias">"#,
        r#"INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description="Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation">"#,
        r#"INFO=<ID=MLEAC,Number=A,Type=Integer,Description="Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed">"#,
        r#"INFO=<ID=MLEAF,Number=A,Type=Float,Description="Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed">"#,
        r#"INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">"#,
        r#"INFO=<ID=MQRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities">"#,
        r#"INFO=<ID=QD,Number=1,Type=Float,Description="Variant Confidence/Quality by Depth">"#,
        r#"INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias">"#,
        r#"INFO=<ID=SOR,Number=1,Type=Float,Description="Symmetric Odds Ratio of 2x2 contingency table to detect strand bias">"#,
        "source=HaplotypeCaller",
    ]
    .iter()
    .map(|s| s.to_string())
    .collect();
    m.push(command_line.to_string());
    m
}

/// The ordered tail of GVCF mode: the stateful genotyper, reference
/// confidence for each region's spans, and the block combiner, fed regions in
/// traversal order.
struct GvcfTail<'a, 'b> {
    ctx: &'a Context<'b>,
    finisher: Finisher,
    blocks: crate::gvcf::BlockCombiner,
    timing: Timing,
    /// Never drawn from: GVCF annotations exclude QD, the only RNG consumer
    /// in `Finisher::finish`; the downsampler's draws stay in the traversal.
    rng: JavaRandom,
    /// Tail time in reference confidence (with block combining) and overall.
    refconf: Duration,
    busy: Duration,
}

/// Reference confidence for one span, with `variants` passed through at
/// their positions.
fn reference_sites(
    ctx: &Context,
    c: usize,
    span: (i64, i64),
    padded: (i64, i64),
    reads: &[HcRead],
    variants: Vec<Record>,
    submit: impl FnMut(Site),
) {
    if span.0 > span.1 {
        return;
    }
    let bases: Vec<u8> = (ctx.reference)(c)[(padded.0 - 1) as usize..padded.1 as usize]
        .iter()
        .map(|&b| reference_base(b))
        .collect();
    crate::refconfidence::emit_sites(span, padded.0, &bases, reads, variants, submit);
}


impl GvcfTail<'_, '_> {
    fn reference_region(
        &mut self,
        c: usize,
        span: (i64, i64),
        padded: (i64, i64),
        reads: &[HcRead],
        variants: Vec<Record>,
    ) {
        let t = Instant::now();
        let blocks = &mut self.blocks;
        reference_sites(self.ctx, c, span, padded, reads, variants, |s| match s {
            Site::Ref(r) => blocks.submit_ref(c, &r),
            Site::Variant(v) => blocks.submit(v),
        });
        self.refconf += t.elapsed();
    }

    /// Everything after `prepare` for one region, in order.
    fn submit_sites(&mut self, c: usize, sites: &[RefSite]) {
        let t = Instant::now();
        for r in sites {
            self.blocks.submit_ref(c, r);
        }
        self.refconf += t.elapsed();
    }

    /// Everything after `prepare` for one region, in order: genotype, then
    /// the precomputed sites of whichever spans the calls select.
    fn finish_region(&mut self, region: &Region, staged: GvcfStaged) -> io::Result<()> {
        let c = region.contig;
        let GvcfStaged {
            prepared,
            whole,
            flanks,
        } = staged;
        if let (Some(p), Some((left, right))) = (prepared, flanks) {
            let variant_span = p.span;
            let variant_padding = p.padded_span;
            let realigned = p.realigned.clone();
            let calls = self
                .finisher
                .finish(self.ctx, p, &mut self.rng, &mut self.timing)?;
            if calls
                .iter()
                .any(|r| r.genotype.alleles.iter().all(|&a| a != NO_CALL_ALLELE))
            {
                self.submit_sites(c, &left);
                self.reference_region(c, variant_span, variant_padding, &realigned, calls);
                self.submit_sites(c, &right);
                return Ok(());
            }
        }
        self.submit_sites(c, &whole);
        Ok(())
    }
}

/// What the pool hands the tail for one region.
struct GvcfStaged {
    prepared: Option<Prepared>,
    /// Reference confidence over the whole region: used unless a call is made.
    whole: Vec<RefSite>,
    /// With a prepared region, reference confidence over the spans left and
    /// right of its variant span: used if a call is made.
    flanks: Option<(Vec<RefSite>, Vec<RefSite>)>,
}

/// The pure, per-region part of GVCF mode: finalized reads, `prepare`, and
/// reference confidence for every span the tail might need (the whole region,
/// and with a prepared region its two flanks), so the ordered tail only
/// genotypes and runs reference confidence over variant spans.
fn gvcf_stage(
    ctx: &Context,
    bam: &RecordSet,
    region: &Region,
    members: &[u32],
    timing: &mut Timing,
) -> io::Result<GvcfStaged> {
    let c = region.contig;
    let len = ctx.contig_lens[c];
    let raw: Vec<_> = members
        .iter()
        .map(|&i| RecordView::new(bam.record(i as usize)).map(|r| HcRead::from_record(&r, i)))
        .collect::<io::Result<_>>()?;
    let finalized = finalize_region(len, region.start, region.end, &raw);
    let prepared = prepare(ctx, bam, region, members, timing)?;
    let sites = |span: (i64, i64)| {
        let mut out = Vec::new();
        if span.0 > span.1 {
            return out;
        }
        let pad = ((span.0 - 100).max(1), (span.1 + 100).min(len));
        let (reads, _) = crate::trim::trim_reads(&finalized, pad.0, pad.1);
        reference_sites(ctx, c, span, pad, &reads, Vec::new(), |s| match s {
            Site::Ref(r) => out.push(r),
            Site::Variant(_) => unreachable!("no variants were passed"),
        });
        out
    };
    let whole = sites((region.start, region.end));
    let flanks = prepared.as_ref().map(|p| {
        (
            sites((region.start, p.span.0 - 1)),
            sites((p.span.1 + 1, region.end)),
        )
    });
    Ok(GvcfStaged {
        prepared,
        whole,
        flanks,
    })
}

/// Single-sample reference-confidence traversal. GVCF annotations exclude QD,
/// so the shared RNG is consumed only by traversal downsampling and region
/// work needs no RNG hand-off: `gvcf_stage` runs on the `threads`-wide pool
/// and [`GvcfTail`] consumes the results in traversal order on its own
/// thread. `NOMICS_GVCF_SERIAL=1` runs every region inline instead.
pub fn run_gvcf(
    ctx: &Context,
    bam: &RecordSet,
    stream: &ReadStream,
    unpadded: &[Vec<Interval>],
    reference: &(dyn Fn(usize, i64) -> u8 + Sync),
    threads: usize,
) -> io::Result<(Traversal, Vec<Record>, Timing)> {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads.max(1))
        .stack_size(64 << 20)
        .build()
        .map_err(io::Error::other)?;
    let tail = || GvcfTail {
        ctx,
        finisher: Finisher::gvcf(),
        blocks: Default::default(),
        timing: Timing::default(),
        rng: JavaRandom::gatk(),
        refconf: Duration::ZERO,
        busy: Duration::ZERO,
    };
    let report = |tail: &GvcfTail| {
        if std::env::var_os("NOMICS_TIMING").is_some() {
            eprintln!(
                "gvcf tail: busy {:.1}s, reference confidence {:.1}s",
                tail.busy.as_secs_f64(),
                tail.refconf.as_secs_f64()
            );
        }
    };

    if std::env::var_os("NOMICS_GVCF_SERIAL").is_some() {
        struct Serial<'a, 'b> {
            bam: &'a RecordSet,
            tail: GvcfTail<'a, 'b>,
            error: Option<io::Error>,
        }
        impl RegionSink for Serial<'_, '_> {
            fn on_call(&mut self, region: &Region, members: &[u32], _rng: &mut JavaRandom) {
                if self.error.is_some() {
                    return;
                }
                let r = gvcf_stage(self.tail.ctx, self.bam, region, members, &mut self.tail.timing)
                    .and_then(|st| self.tail.finish_region(region, st));
                self.error = r.err();
            }
        }
        return pool.install(|| {
            let mut sink = Serial {
                bam,
                tail: tail(),
                error: None,
            };
            let tr = traverse(
                bam,
                stream,
                unpadded,
                reference,
                &ctx.contig_lens,
                &mut JavaRandom::gatk(),
                &mut sink,
                threads,
            )?;
            if let Some(e) = sink.error {
                return Err(e);
            }
            report(&sink.tail);
            Ok((tr, sink.tail.blocks.finish(), sink.tail.timing))
        });
    }

    type Staged = io::Result<GvcfStaged>;
    struct Dispatch<'s, 'c, 'a> {
        ctx: &'c Context<'a>,
        bam: &'c RecordSet,
        scope: &'s rayon::Scope<'c>,
        tx: mpsc::Sender<(u64, Region, Staged, Timing)>,
        sent: u64,
    }
    impl<'c> RegionSink for Dispatch<'_, 'c, '_> {
        fn on_call(&mut self, region: &Region, members: &[u32], _rng: &mut JavaRandom) {
            let (ctx, bam, tx) = (self.ctx, self.bam, self.tx.clone());
            let (seq, region, members) = (self.sent, *region, members.to_vec());
            self.scope.spawn(move |_| {
                let mut timing = Timing::default();
                let r = gvcf_stage(ctx, bam, &region, &members, &mut timing);
                let _ = tx.send((seq, region, r, timing));
            });
            self.sent += 1;
        }
    }

    let (tx, rx) = mpsc::channel::<(u64, Region, Staged, Timing)>();
    std::thread::scope(|threads_scope| {
        let mut tail = tail();
        let tail_thread = std::thread::Builder::new()
            .name("hc-gvcf-tail".into())
            .stack_size(256 << 20)
            .spawn_scoped(threads_scope, move || -> io::Result<(Vec<Record>, Timing)> {
                let mut ready = BTreeMap::new();
                let mut next = 0u64;
                let mut error = None;
                for (seq, region, r, t) in rx {
                    for (a, b) in tail.timing.iter_mut().zip(t) {
                        *a += b;
                    }
                    ready.insert(seq, (region, r));
                    while let Some((region, r)) = ready.remove(&next) {
                        next += 1;
                        if error.is_none() {
                            let t = Instant::now();
                            error = r
                                .and_then(|st| tail.finish_region(&region, st))
                                .err();
                            tail.busy += t.elapsed();
                        }
                    }
                }
                report(&tail);
                match error {
                    Some(e) => Err(e),
                    None => Ok((tail.blocks.finish(), tail.timing)),
                }
            })?;
        let tr = pool.in_place_scope(|scope| {
            let mut dispatch = Dispatch {
                ctx,
                bam,
                scope,
                tx,
                sent: 0,
            };
            traverse(
                bam,
                stream,
                unpadded,
                reference,
                &ctx.contig_lens,
                &mut JavaRandom::gatk(),
                &mut dispatch,
                threads,
            )
        });
        // Every task has joined and the dispatcher's sender is gone.
        let (records, timing) = tail_thread.join().expect("GVCF tail panicked")?;
        Ok((tr?, records, timing))
    })
}
