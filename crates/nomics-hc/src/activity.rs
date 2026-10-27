//! Active-region determination (spec §2): per-locus pileups over the
//! unpadded intervals, the ref-vs-any activity probability, the band-pass
//! activity profile and region cutting, as GATK 4.6.2.0's
//! `AssemblyRegionIterator` + `BandPassActivityProfile` produce them.
//!
//! Gate: the regions of `--assembly-region-out` (start, end, active).

use nomics_io::bam::{RecordSet, RecordView};

use crate::reads::Interval;
use crate::tables::{
    approximate_log10_sum_log10, log10_one_third, max_element_index, normal_distribution,
    normalize_from_log10_to_linear, pl_rounded, qual_to_error_prob_log10, qual_to_prob_log10,
};

/// `--active-probability-threshold`.
const ACTIVE_PROB_THRESHOLD: f64 = 0.002;
/// `--max-prob-propagation-distance`.
const MAX_PROB_PROPAGATION: i64 = 50;
const MIN_REGION_SIZE: usize = 50;
const MAX_REGION_SIZE: usize = 300;
const MAX_FILTER_SIZE: usize = 50;
const SIGMA: f64 = 17.0;
const MIN_PROB_TO_KEEP_IN_FILTER: f64 = 1e-5;
/// `--min-base-quality-score`.
const MIN_BASE_QUAL: u8 = 10;
/// `--reference-model-deletion-quality`.
const DELETION_QUAL_FOR_MODEL: u8 = 30;
/// `ReadUtils.DEFAULT_ADAPTOR_SIZE`.
const ADAPTOR_SIZE: i32 = 100;
const HQ_SOFTCLIP_QUAL: u8 = 28;
const AVERAGE_HQ_SOFTCLIPS_THRESHOLD: f64 = 6.0;

/// The biallelic single-sample prior terms
/// `logToLog10(binomialCoefficientLog(2, n) + logGamma(n + snp) + logGamma(2 - n + ref))`
/// with GATK's default pseudocounts (ref 10.0, snp 0.01), exactly as the
/// reference jar computes them (commons-math 3.5 `FastMath`; dumped by
/// `reference/oracle/hc/GatkTables.java`).
const BIALLELIC_PRIOR_TERMS: [u64; 3] = [
    0x4023_3286_c49f_86cc,
    0x401b_6eeb_e634_d3df,
    0x4016_3f17_83a7_242a,
];

/// Receives regions at GATK's call time. Work that draws from the shared
/// RNG must be done in region order, and every region already handed to
/// `on_call` must have finished its draws when `before_draw` returns: the
/// downsampler draws next.
pub trait RegionSink {
    fn on_call(&mut self, region: &Region, reads: &[u32], rng: &mut crate::jrandom::JavaRandom);
    fn before_draw(&mut self, _rng: &mut crate::jrandom::JavaRandom) {}
    /// The downsampler has just drawn; `rng` is the state later regions see.
    fn after_draw(&mut self, _rng: &mut crate::jrandom::JavaRandom) {}
}

impl<F: FnMut(&Region, &[u32], &mut crate::jrandom::JavaRandom)> RegionSink for F {
    fn on_call(&mut self, region: &Region, reads: &[u32], rng: &mut crate::jrandom::JavaRandom) {
        self(region, reads, rng)
    }
}

/// One assembly region (1-based closed span on a contig).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Region {
    pub contig: usize,
    pub start: i64,
    pub end: i64,
    pub active: bool,
}

/// What a read contributes at one locus (`PileupElement`).
struct Element {
    is_deletion: bool,
    base: u8,
    qual: u8,
    is_alt_flags: bool,
    next_to_soft_clip: bool,
}

/// A read with the fields pileups need.
struct Read<'a> {
    view: RecordView<'a>,
    cigar: Vec<u32>,
    start: i64,
    end: i64,
    /// `isBaseInsideAdaptor` as (reverse, boundary), when it can be true.
    adaptor: Option<(bool, i64)>,
    /// CIGAR cursor: (element index, reference position and read offset at
    /// its start). Loci only move forward, so it only advances.
    cursor: (usize, i64, i64),
}

impl<'a> Read<'a> {
    fn new(view: RecordView<'a>) -> Self {
        let cigar: Vec<u32> = view.cigar().collect();
        let start = view.pos() as i64 + 1;
        let end = start + view.reference_length() as i64 - 1;
        let mut r = Self {
            view,
            cigar,
            start,
            end,
            adaptor: None,
            cursor: (0, start, 0),
        };
        r.adaptor = r.adaptor_boundary();
        r
    }

    fn inside_adaptor(&self, pos: i64) -> bool {
        match self.adaptor {
            Some((true, boundary)) => pos <= boundary,
            Some((false, boundary)) => pos >= boundary,
            None => false,
        }
    }

    /// Move the cursor to the first reference-consuming element that
    /// contains `pos` (or past the end), as `element`'s walk from the start
    /// would.
    fn advance(&mut self, pos: i64) {
        let mut cur = self.cursor;
        self.advance_cursor(&mut cur, pos);
        self.cursor = cur;
    }

    fn advance_cursor(&self, cursor: &mut (usize, i64, i64), pos: i64) {
        let (mut i, mut p, mut ro) = *cursor;
        while i < self.cigar.len() {
            let (op, len) = (self.cigar[i] & 0xf, (self.cigar[i] >> 4) as i64);
            match op {
                0 | 7 | 8 => {
                    if pos < p + len {
                        break;
                    }
                    p += len;
                    ro += len;
                }
                2 | 3 => {
                    if pos < p + len {
                        break;
                    }
                    p += len;
                }
                1 | 4 => ro += len,
                _ => {}
            }
            i += 1;
        }
        *cursor = (i, p, ro);
    }

    /// `element(pos)` from the cursor (after `advance(pos)`).
    fn element_at(&self, pos: i64) -> Option<Element> {
        self.element_with(self.cursor, pos)
    }

    fn element_with(&self, cursor: (usize, i64, i64), pos: i64) -> Option<Element> {
        let (i, p, ro) = cursor;
        let c = *self.cigar.get(i)?;
        let (op, len) = (c & 0xf, (c >> 4) as i64);
        match op {
            0 | 7 | 8 => Some(self.make(i, pos - p, len, Some((ro + pos - p) as usize))),
            2 => Some(self.make(i, pos - p, len, None)),
            _ => None,
        }
    }

    /// `ReadUtils.isBaseInsideAdaptor`, reduced to its boundary.
    fn adaptor_boundary(&self) -> Option<(bool, i64)> {
        let flag = self.view.flag();
        let tlen = self.view.tlen();
        let rev = flag & 0x10 != 0;
        let mate_start = self.view.mate_pos() as i64 + 1;
        // hasWellDefinedFragmentSize
        let well_defined = tlen != 0
            && flag & 0x1 != 0
            && flag & 0x4 == 0
            && flag & 0x8 == 0
            && rev != (flag & 0x20 != 0)
            && if rev {
                self.end > mate_start
            } else {
                self.start <= mate_start + tlen as i64
            };
        if !well_defined || tlen > ADAPTOR_SIZE {
            return None;
        }
        let boundary = if rev {
            mate_start - 1
        } else {
            self.start + tlen.unsigned_abs() as i64
        };
        Some((rev, boundary))
    }

    /// The pileup element at reference position `pos`, if any (M/=/X or D;
    /// never I, S, N).
    fn element(&self, pos: i64) -> Option<Element> {
        let ops = &self.cigar;
        let (mut p, mut ro) = (self.start, 0i64);
        for (i, &c) in ops.iter().enumerate() {
            let (op, len) = (c & 0xf, (c >> 4) as i64);
            match op {
                0 | 7 | 8 => {
                    if pos < p + len {
                        let off = pos - p;
                        let r = (ro + off) as usize;
                        return Some(self.make(i, off, len, Some(r)));
                    }
                    p += len;
                    ro += len;
                }
                2 => {
                    if pos < p + len {
                        return Some(self.make(i, pos - p, len, None));
                    }
                    p += len;
                }
                3 => {
                    if pos < p + len {
                        return None;
                    }
                    p += len;
                }
                1 | 4 => ro += len,
                _ => {}
            }
        }
        None
    }

    fn make(&self, i: usize, off: i64, len: i64, read_offset: Option<usize>) -> Element {
        let ops = &self.cigar;
        let op_at = |j: isize| -> Option<u32> {
            (j >= 0 && (j as usize) < ops.len()).then(|| ops[j as usize] & 0xf)
        };
        let at_start = off == 0;
        let at_end = off == len - 1;
        let on_genome = |o: u32| matches!(o, 0 | 7 | 8 | 2);
        let nearest = |dir: isize| -> Option<u32> {
            let mut j = i as isize + dir;
            while j >= 0 && (j as usize) < ops.len() {
                let o = ops[j as usize] & 0xf;
                if on_genome(o) {
                    return Some(o);
                }
                j += dir;
            }
            None
        };
        let is_deletion = ops[i] & 0xf == 2;
        let before_del = !is_deletion && at_end && nearest(1) == Some(2);
        let after_del = !is_deletion && at_start && nearest(-1) == Some(2);
        let before_ins = at_end && op_at(i as isize + 1) == Some(1);
        let after_ins = at_start && op_at(i as isize - 1) == Some(1);
        let after_soft = at_start && op_at(i as isize - 1) == Some(4);
        let before_soft = at_end && op_at(i as isize + 1) == Some(4);
        let (base, qual) = match read_offset {
            Some(r) => (self.view.base(r), self.view.quals().unwrap()[r]),
            None => (b'D', 16),
        };
        Element {
            is_deletion,
            base,
            qual,
            is_alt_flags: before_del
                || after_del
                || before_ins
                || after_ins
                || after_soft
                || before_soft,
            next_to_soft_clip: after_soft || before_soft,
        }
    }

    /// `AlignmentUtils.countHighQualitySoftClips(read, 28)`.
    fn hq_soft_clips(&self) -> u32 {
        let q = self.view.quals().unwrap();
        let (mut n, mut at) = (0, 0usize);
        for &c in &self.cigar {
            let (op, len) = (c & 0xf, (c >> 4) as usize);
            if op == 4 {
                n += q[at..at + len]
                    .iter()
                    .filter(|&&x| x > HQ_SOFTCLIP_QUAL)
                    .count() as u32;
                at += len;
            } else if matches!(op, 0 | 1 | 7 | 8) {
                at += len;
            }
        }
        n
    }
}

/// `ActivityProfileState`: probability, plus the HQ-soft-clip expansion.
#[derive(Clone, Copy)]
struct State {
    pos: i64,
    prob: f64,
    hq_soft_clips: Option<f64>,
}

/// `HaplotypeCallerEngine.isActive` on one locus's pileup.
fn is_active(
    elements: &[(Element, &Read)],
    ref_base: u8,
    log10_ploidy: f64,
    one_third: f64,
) -> (f64, Option<f64>) {
    if elements.is_empty() {
        return (0.0, None);
    }
    // ReferenceConfidenceModel.calcGenotypeLikelihoodsOfRefVsAny, ploidy 2.
    let mut gl = [0.0f64; 3];
    let mut read_count = 0i64;
    let (mut hq_mean, mut hq_n) = (0.0f64, 0i64);
    for (p, read) in elements {
        let qual = if p.is_deletion {
            DELETION_QUAL_FOR_MODEL
        } else {
            p.qual
        };
        if qual <= MIN_BASE_QUAL && !p.is_deletion {
            continue;
        }
        read_count += 1;
        let is_alt = p.base != ref_base || p.is_deletion || p.is_alt_flags;
        let (r, nr) = if is_alt {
            (
                qual_to_error_prob_log10(qual) + one_third,
                qual_to_prob_log10(qual),
            )
        } else {
            (
                qual_to_prob_log10(qual),
                qual_to_error_prob_log10(qual) + one_third,
            )
        };
        gl[0] += r + log10_ploidy;
        gl[2] += nr + log10_ploidy;
        // i = 1, j = 1: log10(1) = 0.0
        gl[1] += approximate_log10_sum_log10(r + 0.0, nr + 0.0);
        if is_alt && p.next_to_soft_clip {
            let obs = read.hq_soft_clips() as f64;
            hq_n += 1;
            hq_mean += (obs - hq_mean) / hq_n as f64;
        }
    }
    let denominator = read_count as f64 * log10_ploidy;
    for g in gl.iter_mut() {
        *g -= denominator;
    }
    // GenotypeBuilder.PL(double[]) then getLikelihoods().getAsVector().
    let gl = pl_rounded(&gl);
    // AlleleFrequencyCalculator.calculateSingleSampleBiallelicNonRefPosterior(gl, true)
    let prob = if max_element_index(&gl) == 0 {
        0.0
    } else {
        let mut un: Vec<f64> = (0..3)
            .map(|n| gl[n] + f64::from_bits(BIALLELIC_PRIOR_TERMS[n]))
            .collect();
        if max_element_index(&un) == 0 {
            0.0
        } else {
            normalize_from_log10_to_linear(&mut un);
            1.0 - un[0]
        }
    };
    let hq = (hq_mean > AVERAGE_HQ_SOFTCLIPS_THRESHOLD).then_some(hq_mean);
    (prob, hq)
}

/// `BandPassActivityProfile` + `ActivityProfile` for one contig.
struct Profile {
    states: Vec<(i64, f64)>,
    region_start: Option<i64>,
    region_stop: i64,
    contig_len: i64,
    kernel: Vec<f64>,
    filter_size: i64,
}

fn make_kernel(filter_size: usize) -> Vec<f64> {
    let band = 2 * filter_size + 1;
    let k: Vec<f64> = (0..band)
        .map(|i| normal_distribution(filter_size as f64, SIGMA, i as f64))
        .collect();
    let mut sum = 0.0;
    for &v in &k {
        sum += v;
    }
    k.iter().map(|&v| v / sum).collect()
}

impl Profile {
    fn new(contig_len: i64) -> Self {
        let full = make_kernel(MAX_FILTER_SIZE);
        // determineFilterSize
        let middle = (full.len() - 1) / 2;
        let mut filter_end = middle;
        while filter_end > 0 {
            if full[filter_end - 1] < MIN_PROB_TO_KEEP_IN_FILTER {
                break;
            }
            filter_end -= 1;
        }
        let filter_size = middle - filter_end;
        Self {
            states: Vec::new(),
            region_start: None,
            region_stop: 0,
            contig_len,
            kernel: make_kernel(filter_size),
            filter_size: filter_size as i64,
        }
    }

    fn max_prob_propagation(&self) -> i64 {
        MAX_PROB_PROPAGATION + self.filter_size
    }

    fn is_empty(&self) -> bool {
        self.states.is_empty()
    }

    fn end(&self) -> i64 {
        self.region_stop
    }

    fn loc_ok(&self, pos: i64) -> bool {
        pos >= 1 && pos <= self.contig_len
    }

    fn incorporate(&mut self, pos: i64, prob: f64) {
        let position = pos - self.region_start.unwrap();
        assert!(position <= self.states.len() as i64);
        if position >= 0 {
            let p = position as usize;
            if p < self.states.len() {
                self.states[p].1 += prob;
            } else {
                self.states.push((pos, prob));
            }
        }
    }

    fn add(&mut self, s: State) {
        if self.region_start.is_none() {
            self.region_start = Some(s.pos);
        } else {
            assert_eq!(self.region_stop, s.pos - 1, "non-contiguous profile add");
        }
        self.region_stop = s.pos;
        // Common case (no soft-clip spread, zero probability): the state
        // passes through unconvolved, as below, without the temporaries.
        if s.hq_soft_clips.is_none() && s.prob <= 0.0 {
            self.incorporate(s.pos, s.prob);
            return;
        }
        // ActivityProfile.processState: HQ soft clips spread the state.
        let supers: Vec<f64> = match s.hq_soft_clips {
            Some(v) => {
                let n = (v as i64).min(self.max_prob_propagation());
                (-n..=n)
                    .filter(|&i| self.loc_ok(s.pos + i))
                    .map(|_| s.prob)
                    .collect()
            }
            None => vec![s.prob],
        };
        // BandPassActivityProfile.processState: convolve each around the
        // *just-added* locus; zero-probability states pass the original through.
        let mut out = Vec::new();
        for p in supers {
            if p > 0.0 {
                for i in -self.filter_size..=self.filter_size {
                    let pos = s.pos + i;
                    if self.loc_ok(pos) {
                        out.push((pos, p * self.kernel[(i + self.filter_size) as usize]));
                    }
                }
            } else {
                out.push((s.pos, s.prob));
            }
        }
        for (pos, prob) in out {
            self.incorporate(pos, prob);
        }
    }

    /// Every queued probability is at most the active threshold, so the
    /// next regions are inactive and cut only by length.
    fn all_inactive(&self) -> bool {
        self.states.iter().all(|s| s.1 <= ACTIVE_PROB_THRESHOLD)
    }

    /// `add` of zero-probability states without soft clips for every locus
    /// up to `to` (exclusive): each passes through unconvolved, appending
    /// where the profile ends and adding nothing to states already queued.
    fn extend_zeros(&mut self, to: i64) {
        let start = self.region_start.expect("non-empty profile");
        while start + (self.states.len() as i64) < to {
            let pos = start + self.states.len() as i64;
            self.states.push((pos, 0.0));
        }
        self.region_stop = to - 1;
    }

    /// First locus after `pos` whose (non-forced) pop_ready can cut a
    /// region: there the queue, `max(len, locus - start)`, reaches
    /// MAX_REGION_SIZE + propagation.
    fn next_pop_locus(&self, pos: i64) -> i64 {
        let need = MAX_REGION_SIZE as i64 + self.max_prob_propagation();
        if self.states.len() as i64 >= need {
            pos + 1
        } else {
            (self.region_start.expect("non-empty profile") + need).max(pos + 1)
        }
    }

    fn find_end(&self, is_active: bool, force: bool) -> Option<usize> {
        if !force
            && (self.states.len() as i64) < MAX_REGION_SIZE as i64 + self.max_prob_propagation()
        {
            return None;
        }
        let mut end = 0usize;
        while end < self.states.len() && end < MAX_REGION_SIZE {
            if (self.states[end].1 > ACTIVE_PROB_THRESHOLD) != is_active {
                break;
            }
            end += 1;
        }
        if is_active && end == MAX_REGION_SIZE {
            // findBestCutSite
            let mut min_i = end - 1;
            let mut min_p = f64::MAX;
            let mut i = end as isize - 1;
            while i >= MIN_REGION_SIZE as isize - 1 {
                let iu = i as usize;
                let cur = self.states[iu].1;
                let is_min = iu != self.states.len() - 1
                    && iu >= 1
                    && cur <= self.states[iu + 1].1
                    && cur < self.states[iu - 1].1;
                if cur < min_p && is_min {
                    min_p = cur;
                    min_i = iu;
                }
                i -= 1;
            }
            end = min_i + 1;
        }
        Some(end - 1)
    }

    fn pop_ready(&mut self, contig: usize, force: bool, out: &mut Vec<Region>) {
        loop {
            if self.states.is_empty() {
                return;
            }
            if force {
                let span = (self.region_stop - self.region_start.unwrap() + 1) as usize;
                self.states.truncate(span);
            }
            let first = self.states[0];
            let is_active = first.1 > ACTIVE_PROB_THRESHOLD;
            let Some(off) = self.find_end(is_active, force) else {
                return;
            };
            self.states.drain(..=off);
            if self.states.is_empty() {
                self.region_start = None;
            } else {
                self.region_start = Some(self.states[0].0);
            }
            out.push(Region {
                contig,
                start: first.0,
                end: first.0 + off as i64,
                active: is_active,
            });
        }
    }
}

/// Result of a traversal.
pub struct Traversal {
    /// Every assembly region (active and inactive), in traversal order.
    pub regions: Vec<Region>,
    /// `nextInt` draws made by the downsampler.
    pub draws: u64,
}

/// Per-locus `isActive` results computed off the traversal thread.
struct Precomputed {
    /// (probability, HQ soft clips) by locus ordinal along the intervals.
    states: Vec<(f64, Option<f64>)>,
    /// Loci a downsampled start position's reads may cover: their pileups
    /// depend on the RNG, so the traversal computes them itself.
    dirty: Vec<(i64, i64)>,
}

/// `isActive` for every locus of a contig's intervals, from all filtered
/// reads (no downsampling), in `threads` chunks. Away from start positions
/// with more than `MAX_READS_PER_START` reads the downsampler keeps every
/// read in stream order, so these pileups equal the traversal's exactly.
#[allow(clippy::too_many_arguments)]
fn precompute_activity(
    all: &[Read],
    groups: &[crate::reads::Group],
    ivs: &[Interval],
    c: usize,
    reference: &(dyn Fn(usize, i64) -> u8 + Sync),
    max_span: i64,
    threads: usize,
    log10_ploidy: f64,
    one_third: f64,
) -> Precomputed {
    let mut dirty = Vec::new();
    let mut at = 0usize;
    for g in groups {
        let n = g.reads.len();
        if n > crate::reads::MAX_READS_PER_START {
            let end = all[at..at + n].iter().map(|r| r.end).max().unwrap();
            // Starts ascend; merge overlaps so a forward cursor can test loci.
            match dirty.last_mut() {
                Some((_, b)) if g.pos <= *b => *b = (*b).max(end),
                _ => dirty.push((g.pos, end)),
            }
        }
        at += n;
    }
    // Locus ordinals: interval k starts at ordinal offsets[k].
    let mut offsets = Vec::with_capacity(ivs.len() + 1);
    let mut total = 0usize;
    for iv in ivs {
        offsets.push(total);
        total += (iv.end - iv.start + 1) as usize;
    }
    offsets.push(total);
    let locus = |o: usize| -> i64 {
        let k = offsets.partition_point(|&x| x <= o) - 1;
        ivs[k].start + (o - offsets[k]) as i64
    };
    let chunks = (threads * 8).max(1);
    let per = total.div_ceil(chunks).max(1);
    let ranges: Vec<(usize, usize)> = (0..total)
        .step_by(per)
        .map(|o| (o, (o + per).min(total)))
        .collect();
    let next = std::sync::atomic::AtomicUsize::new(0);
    let mut parts: Vec<Vec<(f64, Option<f64>)>> = vec![Vec::new(); ranges.len()];
    let slots: Vec<std::sync::Mutex<&mut Vec<(f64, Option<f64>)>>> =
        parts.iter_mut().map(std::sync::Mutex::new).collect();
    std::thread::scope(|s| {
        for _ in 0..threads.max(1) {
            s.spawn(|| {
                let mut elements: Vec<(Element, &Read)> = Vec::new();
                // (read, its CIGAR cursor), in stream order.
                let mut active: Vec<(usize, (usize, i64, i64))> = Vec::new();
                loop {
                    let k = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    let Some(&(o0, o1)) = ranges.get(k) else {
                        return;
                    };
                    let first = locus(o0);
                    let mut admitted = all.partition_point(|r| r.start < first - max_span);
                    active.clear();
                    let mut out = Vec::with_capacity(o1 - o0);
                    for o in o0..o1 {
                        let pos = locus(o);
                        while admitted < all.len() && all[admitted].start <= pos {
                            active.push((admitted, (0, all[admitted].start, 0)));
                            admitted += 1;
                        }
                        active.retain(|&(i, _)| all[i].end >= pos);
                        elements.clear();
                        for (i, cur) in active.iter_mut() {
                            let r = &all[*i];
                            r.advance_cursor(cur, pos);
                            if !r.inside_adaptor(pos)
                                && let Some(e) = r.element_with(*cur, pos)
                            {
                                elements.push((e, r));
                            }
                        }
                        out.push(is_active(
                            &elements,
                            reference(c, pos),
                            log10_ploidy,
                            one_third,
                        ));
                    }
                    **slots[k].lock().unwrap() = out;
                }
            });
        }
    });
    drop(slots);
    Precomputed {
        states: parts.into_iter().flatten().collect(),
        dirty,
    }
}

/// First position `>= x` with a non-empty LIBS context (some read with an
/// M/=/X/D element there, outside its adaptor), or `i64::MAX`. `reads` are
/// the contig's filtered reads sorted by start; `max_span` bounds how far
/// back a covering read can start.
fn next_context(reads: &[Read], x: i64, max_span: i64) -> i64 {
    let mut y = x;
    loop {
        let hi = reads.partition_point(|r| r.start <= y);
        let mut covered = false;
        let mut i = hi;
        while i > 0 {
            i -= 1;
            let r = &reads[i];
            if r.start < y - max_span {
                break;
            }
            if r.end >= y {
                covered = true;
                if r.element(y).is_some() && !r.inside_adaptor(y) {
                    return y;
                }
            }
        }
        if covered {
            y += 1;
        } else if hi < reads.len() {
            y = reads[hi].start; // nothing covers y: jump to the next read
        } else {
            return i64::MAX;
        }
    }
}

/// `HaplotypeCaller`'s traversal (spec §1.5, §2): per contig, interval loci in
/// order, pileups from the downsampled stream, the activity profile and
/// region cutting; `on_call(region, reads, rng)` runs once per region (active
/// or not) at the point GATK calls it, with the region's reads (BAM record
/// indices, in `AssemblyRegion` order: the previous region's reads that overlap
/// the padded span, then newly cached reads in stream order), so draws from the shared `rng` interleave
/// with the downsampler's exactly as in GATK:
///
/// * Start position `P_k`'s reservoir draws happen when the last read of the
///   previous position leaves the push-to-pull downsampler: during the LIBS
///   step at `P_{k-1}`, or at `P_{k-2}` if `P_{k-1}` holds one read (the
///   LIBS end-of-step peek pulls it).
/// * `AssemblyRegionIterator.next()` loads region `R+1` before returning `R`,
///   and the interval locus iterator advances LIBS to the first non-empty
///   context at or after the next interval locus; region `R` is called after
///   that (or after LIBS is drained at the end of the contig).
pub fn traverse(
    bam: &RecordSet,
    stream: &crate::reads::ReadStream,
    unpadded: &[Vec<Interval>],
    reference: &(dyn Fn(usize, i64) -> u8 + Sync),
    contig_lens: &[i64],
    rng: &mut crate::jrandom::JavaRandom,
    sink: &mut dyn RegionSink,
    threads: usize,
) -> std::io::Result<Traversal> {
    let log10_ploidy = crate::jmath::log10(2.0);
    let timing_on = std::env::var_os("NOMICS_TIMING").is_some();
    // ingest, pileup, isActive, profile + regions
    let mut tm = [std::time::Duration::ZERO; 4];
    let mut tick = |k: usize, t: &mut std::time::Instant| {
        if timing_on {
            let now = std::time::Instant::now();
            tm[k] += now - *t;
            *t = now;
        }
    };
    let one_third = log10_one_third();
    let mut out = Traversal {
        regions: Vec::new(),
        draws: 0,
    };
    for (c, groups) in stream.groups.iter().enumerate() {
        let ivs = &unpadded[c];
        // All filtered reads of the contig, by start (for context lookahead).
        let all: Vec<Read> = groups
            .iter()
            .flat_map(|g| g.reads.iter())
            .map(|&i| RecordView::new(bam.record(i as usize)).map(Read::new))
            .collect::<Result<_, _>>()?;
        let max_span = all.iter().map(|r| r.end - r.start + 1).max().unwrap_or(0);
        // Draw time of each group (LIBS location; i64::MIN = at construction).
        let t: Vec<i64> = (0..groups.len())
            .map(|k| match k {
                0 => i64::MIN,
                _ if groups[k - 1].reads.len() >= 2 => groups[k - 1].pos,
                1 => i64::MIN,
                _ => groups[k - 2].pos,
            })
            .collect();
        // Region calls waiting for their time: (LIBS location, region).
        type Call = (i64, Region, Vec<u32>);
        let mut calls: std::collections::VecDeque<Call> = Default::default();
        let fire_until = |limit: i64,
                          sink: &mut dyn RegionSink,
                          calls: &mut std::collections::VecDeque<Call>,
                          rng: &mut crate::jrandom::JavaRandom| {
            while calls.front().is_some_and(|(y, _, _)| *y < limit) {
                let (_, r, m) = calls.pop_front().unwrap();
                sink.on_call(&r, &m, rng);
            }
        };
        if ivs.is_empty() {
            // No loci, so no shard and no traversal: these reads are never
            // queried (padded intervals exist only around loci).
            continue;
        }
        // isActive off this thread, except where downsampling can matter.
        let pre = (threads > 1).then(|| {
            crate::caller::progress::IN_PRECOMPUTE
                .store(true, std::sync::atomic::Ordering::Relaxed);
            let p = precompute_activity(
                &all,
                groups,
                ivs,
                c,
                reference,
                max_span,
                threads,
                log10_ploidy,
                one_third,
            );
            crate::caller::progress::IN_PRECOMPUTE
                .store(false, std::sync::atomic::Ordering::Relaxed);
            p
        });
        let mut li = 0usize; // locus ordinal along the intervals
        let mut dirty_at = 0usize; // first merged dirty span not ending before `pos`
        let mut active_valid = true;
        let mut reads: Vec<Read> = Vec::new(); // downsampled, stream order
        let mut ids: Vec<u32> = Vec::new(); // BAM record index of each of `reads`
        // fillNextAssemblyRegionWithReads state: reads[cache_at..] are cached,
        // not yet offered to a region; `prev` are the last region's members.
        let mut cache_at = 0usize;
        let mut prev: Vec<usize> = Vec::new();
        let padded_start = |r: &Region| (r.start - crate::reads::REGION_PADDING).max(1);
        let fill = |r: &Region,
                    reads: &[Read],
                    prev: &mut Vec<usize>,
                    cache_at: &mut usize|
         -> Vec<usize> {
            let (ps, pe) = (
                padded_start(r),
                (r.end + crate::reads::REGION_PADDING).min(contig_lens[c]),
            );
            let mut m: Vec<usize> = prev
                .iter()
                .copied()
                .filter(|&i| reads[i].start <= pe && reads[i].end >= ps)
                .collect();
            while *cache_at < reads.len() && reads[*cache_at].start <= pe {
                if reads[*cache_at].end >= ps {
                    m.push(*cache_at);
                }
                *cache_at += 1;
            }
            *prev = m.clone();
            m
        };
        let mut next_group = 0usize;
        let mut admitted = 0usize; // reads[..admitted] have joined `active`
        let mut active: Vec<usize> = Vec::new();
        let mut profile = Profile::new(contig_lens[c]);
        let mut pending: std::collections::VecDeque<Region> = Default::default();
        let mut emitted: Option<(Region, Vec<u32>)> = None; // last emitted, not yet called
        let padded_end = |r: &Region| (r.end + crate::reads::REGION_PADDING).min(contig_lens[c]);
        let mut iv_i = 0usize;
        let mut pos = ivs[0].start;
        let mut clk = std::time::Instant::now();
        let mut nz_at = 0usize; // next locus ordinal with a nonzero precomputed state
        loop {
            // Read-free stretches: every locus has probability 0 and no soft
            // clips, no group enters, nothing is downsampled and the profile
            // holds only inactive states. Then the only events are length
            // cuts of inactive regions and region emissions, so jump between
            // them instead of stepping every locus (exact: same regions,
            // members and call times as the per-locus loop below).
            if let Some(p) = pre.as_ref()
                && !profile.is_empty()
                && pos == profile.end() + 1
            {
                while nz_at < p.states.len()
                    && (nz_at < li || (p.states[nz_at].0 == 0.0 && p.states[nz_at].1.is_none()))
                {
                    nz_at += 1;
                }
                while dirty_at < p.dirty.len() && p.dirty[dirty_at].1 < pos {
                    dirty_at += 1;
                }
                let mut run_end = ivs[iv_i].end; // the interval's last locus runs normally
                run_end = run_end.min(pos + (nz_at - li) as i64);
                if next_group < groups.len() {
                    run_end = run_end.min(groups[next_group].pos);
                }
                if dirty_at < p.dirty.len() {
                    run_end = run_end.min(p.dirty[dirty_at].0.max(pos));
                }
                if run_end - pos >= 1024 && profile.all_inactive() {
                    active_valid = false;
                    admitted = reads.len();
                    while pos < run_end {
                        let mut popped = Vec::new();
                        profile.pop_ready(c, false, &mut popped);
                        pending.extend(popped);
                        profile.extend_zeros(pos + 1);
                        if pending.front().is_some_and(|r| pos > padded_end(r)) {
                            let r = pending.pop_front().unwrap();
                            let members: Vec<u32> = fill(&r, &reads, &mut prev, &mut cache_at)
                                .into_iter()
                                .map(|i| ids[i])
                                .collect();
                            if let Some((pr, pm)) = emitted.replace((r, members)) {
                                let y = next_context(&all, pos + 1, max_span);
                                calls.push_back((y, pr, pm));
                            }
                            out.regions.push(r);
                        }
                        // Skip to the next locus where a pop or an emission can occur.
                        let mut next = profile.next_pop_locus(pos).min(run_end);
                        if let Some(r) = pending.front() {
                            next = next.min((padded_end(r) + 1).max(pos + 1));
                        }
                        profile.extend_zeros(next);
                        li += (next - pos) as usize;
                        pos = next;
                    }
                    continue;
                }
            }
            // Groups starting at or before this locus enter the stream.
            while next_group < groups.len() && groups[next_group].pos <= pos {
                fire_until(t[next_group], &mut *sink, &mut calls, rng);
                if groups[next_group].reads.len() > crate::reads::MAX_READS_PER_START {
                    sink.before_draw(rng);
                }
                let (kept, d) = crate::reads::downsample(&groups[next_group], rng);
                if d > 0 {
                    sink.after_draw(rng);
                }
                out.draws += d;
                for i in kept {
                    reads.push(Read::new(RecordView::new(bam.record(i as usize))?));
                    ids.push(i);
                }
                next_group += 1;
            }
            let dirty = pre.as_ref().is_none_or(|p| {
                while dirty_at < p.dirty.len() && p.dirty[dirty_at].1 < pos {
                    dirty_at += 1;
                }
                dirty_at < p.dirty.len() && p.dirty[dirty_at].0 <= pos
            });
            if dirty {
                if active_valid {
                    active.extend(admitted..reads.len());
                    active.retain(|&i| reads[i].end >= pos);
                } else {
                    // Admitted reads still covering `pos`, in admission order.
                    let lo = reads.partition_point(|r| r.start < pos - max_span);
                    active = (lo..reads.len()).filter(|&i| reads[i].end >= pos).collect();
                    active_valid = true;
                }
            } else {
                active_valid = false;
            }
            admitted = reads.len();
            tick(0, &mut clk);
            // Pop ready regions before this locus is added.
            if !profile.is_empty() {
                let force = pos != profile.end() + 1;
                let mut popped = Vec::new();
                profile.pop_ready(c, force, &mut popped);
                pending.extend(popped);
            }
            let (prob, hq) = if dirty {
                for &i in &active {
                    reads[i].advance(pos);
                }
                let mut elements: Vec<(Element, &Read)> = Vec::with_capacity(active.len());
                for &i in &active {
                    let r = &reads[i];
                    if !r.inside_adaptor(pos)
                        && let Some(e) = r.element_at(pos)
                    {
                        elements.push((e, r));
                    }
                }
                tick(1, &mut clk);
                is_active(&elements, reference(c, pos), log10_ploidy, one_third)
            } else {
                pre.as_ref().unwrap().states[li]
            };
            li += 1;
            tick(2, &mut clk);
            profile.add(State {
                pos,
                prob,
                hq_soft_clips: hq,
            });
            let next_pos = if pos < ivs[iv_i].end {
                Some(pos + 1)
            } else if iv_i + 1 < ivs.len() {
                Some(ivs[iv_i + 1].start)
            } else {
                None
            };
            // At most one region is emitted per locus.
            if pending.front().is_some_and(|r| pos > padded_end(r)) {
                let r = pending.pop_front().unwrap();
                let members: Vec<u32> = fill(&r, &reads, &mut prev, &mut cache_at)
                    .into_iter()
                    .map(|i| ids[i])
                    .collect();
                if let Some((p, pm)) = emitted.replace((r, members)) {
                    // `p` is called after this region is loaded, i.e. once
                    // LIBS reached the first context at/after the next locus.
                    let y = match next_pos {
                        Some(x) => next_context(&all, x, max_span),
                        None => i64::MAX,
                    };
                    calls.push_back((y, p, pm));
                }
                out.regions.push(r);
            }
            tick(3, &mut clk);
            match next_pos {
                Some(x) => {
                    if x != pos + 1 {
                        iv_i += 1;
                    }
                    pos = x;
                }
                None => break,
            }
        }
        // End of the shard: LIBS is drained first (the remaining groups are
        // read and downsampled), then the profile is forced out; each
        // remaining region is loaded by a next() call after that.
        while next_group < groups.len() {
            fire_until(t[next_group], &mut *sink, &mut calls, rng);
            if groups[next_group].reads.len() > crate::reads::MAX_READS_PER_START {
                sink.before_draw(rng);
            }
            let (kept, d) = crate::reads::downsample(&groups[next_group], rng);
            if d > 0 {
                sink.after_draw(rng);
            }
            out.draws += d;
            for i in kept {
                reads.push(Read::new(RecordView::new(bam.record(i as usize))?));
                ids.push(i);
            }
            next_group += 1;
        }
        let mut popped = Vec::new();
        profile.pop_ready(c, true, &mut popped);
        pending.extend(popped);
        while let Some(r) = pending.pop_front() {
            let members: Vec<u32> = fill(&r, &reads, &mut prev, &mut cache_at)
                .into_iter()
                .map(|i| ids[i])
                .collect();
            if let Some((p, pm)) = emitted.replace((r, members)) {
                calls.push_back((i64::MAX, p, pm));
            }
            out.regions.push(r);
        }
        if let Some((last, lm)) = emitted.take() {
            calls.push_back((i64::MAX, last, lm));
        }
        fire_until(i64::MAX, &mut *sink, &mut calls, rng);
        while let Some((_, r, m)) = calls.pop_front() {
            sink.on_call(&r, &m, rng);
        }
    }
    if timing_on {
        eprintln!(
            "traverse: ingest {:.1}s, pileup {:.1}s, isActive {:.1}s, profile+regions {:.1}s",
            tm[0].as_secs_f64(),
            tm[1].as_secs_f64(),
            tm[2].as_secs_f64(),
            tm[3].as_secs_f64()
        );
    }
    Ok(out)
}
