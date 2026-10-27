//! Region -> alignment -> SAM text (`bwamem.c` `mem_approx_mapq_se`,
//! `mem_reg2aln`, `mem_aln2sam`, `mem_reg2sam`, `infer_bw`, `get_rlen`;
//! `bwamem_extra.c` `mem_gen_alt`).
//!
//! Options are bwa-mem's defaults plus `-Y` (soft-clip supplementaries);
//! `-a`, `-C`, `-M`, `-5`, `-q`, `-V` and XB are not supported.

use crate::extend::{A, B, E_DEL, E_INS, O_DEL, O_INS, Region, W, gen_cigar2};
use crate::index::Index;
use nomics_io::bam::{BinRule, IntType, Record, Value};
use std::fmt::Write as _;
use std::io;

/// `opt->T`: minimum score to output.
pub const T: i32 = 30;
const MIN_SEED_LEN: i32 = 19;
const XA_DROP_RATIO: f32 = 0.80;
const MAX_XA_HITS: i32 = 5;
const MAX_XA_HITS_ALT: i32 = 200;
/// `opt->mapQ_coef_len` (float) and `opt->mapQ_coef_fac = log(50)` stored
/// in an int, i.e. 3.
const MAPQ_COEF_LEN: f32 = 50.0;
const MAPQ_COEF_FAC: i32 = 3;

/// One read as bwa holds it: name, 2-bit codes (0..4), raw qualities.
#[derive(Debug, Clone, Default)]
pub struct Read {
    pub name: Vec<u8>,
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
}

/// `mem_aln_t`; `md` holds the MD string bwa appends after the CIGAR.
#[derive(Debug, Clone, Default)]
pub struct Aln {
    pub pos: i64,
    pub rid: i32,
    pub flag: i32,
    pub is_rev: bool,
    pub is_alt: bool,
    pub mapq: u8,
    pub nm: i32,
    pub cigar: Vec<u32>,
    pub md: String,
    pub xa: Option<String>,
    pub score: i32,
    pub sub: i32,
    pub alt_sc: i32,
}

/// `mem_approx_mapq_se`. Fused multiply-adds follow the golden binary's
/// object code (bwamem.c:996, 1001, 1004).
pub fn approx_mapq_se(a: &Region) -> i32 {
    let mut sub = if a.sub != 0 { a.sub } else { MIN_SEED_LEN * A };
    sub = sub.max(a.csub);
    if sub >= a.score {
        return 0;
    }
    let l = ((a.qe - a.qb) as i64).max(a.re - a.rb) as i32;
    let identity = 1.0 - (l * A - a.score) as f64 / (A + B) as f64 / l as f64;
    let mut mapq = if a.score == 0 {
        0
    } else {
        let mut tmp = if (l as f32) < MAPQ_COEF_LEN {
            1.0
        } else {
            MAPQ_COEF_FAC as f64 / crate::glibc_math::log(l as f64)
        };
        tmp *= identity * identity;
        let x = (a.score - sub) as f64 * 6.02 / A as f64 * tmp;
        x.mul_add(tmp, 0.499) as i32
    };
    if a.sub_n > 0 {
        mapq -= crate::glibc_math::log((a.sub_n + 1) as f64).mul_add(4.343, 0.499) as i32;
    }
    mapq = mapq.clamp(0, 60);
    (mapq as f64).mul_add(1.0 - a.frac_rep as f64, 0.499) as i32
}

/// `infer_bw`.
fn infer_bw(l1: i32, l2: i32, score: i32, a: i32, q: i32, r: i32) -> i32 {
    if l1 == l2 && l1 * a - score < (q + r - a) << 1 {
        return 0;
    }
    let w = ((l1.min(l2) * a - score - q) as f64 / r as f64 + 2.0) as i32;
    w.max((l1 - l2).abs())
}

/// `get_rlen`: reference length of a CIGAR (M and D).
pub fn get_rlen(cigar: &[u32]) -> i64 {
    cigar
        .iter()
        .filter(|&&c| matches!(c & 0xf, 0 | 2))
        .map(|&c| (c >> 4) as i64)
        .sum()
}

/// `mem_reg2aln`: the final alignment (CIGAR, NM, MD, position) of a region,
/// or an unmapped record for `None`.
pub fn reg2aln(idx: &Index, query: &[u8], ar: Option<&Region>) -> Aln {
    let mut a = Aln::default();
    let ar = match ar {
        Some(r) if r.rb >= 0 && r.re >= 0 => r,
        _ => {
            a.rid = -1;
            a.pos = -1;
            a.flag |= 0x4;
            return a;
        }
    };
    let l_query = query.len() as i32;
    let (qb, qe, rb, re) = (ar.qb, ar.qe, ar.rb, ar.re);
    a.mapq = if ar.secondary < 0 {
        approx_mapq_se(ar) as u8
    } else {
        0
    };
    if ar.secondary >= 0 {
        a.flag |= 0x100;
    }
    let tmp = infer_bw(qe - qb, (re - rb) as i32, ar.truesc, A, O_DEL, E_DEL);
    let mut w2 = infer_bw(qe - qb, (re - rb) as i32, ar.truesc, A, O_INS, E_INS);
    w2 = w2.max(tmp);
    if w2 > W {
        w2 = w2.min(ar.w);
    }
    let mut last_sc = -(1 << 30);
    let mut i = 0;
    let mut g;
    loop {
        w2 = w2.min(W << 2);
        g = gen_cigar2(idx, w2, &query[qb as usize..qe as usize], rb, re, true)
            .expect("mem_reg2aln: bwa_gen_cigar2 failed");
        if g.score == last_sc || w2 == W << 2 {
            break;
        }
        last_sc = g.score;
        w2 <<= 1;
        i += 1;
        if !(i < 3 && g.score < ar.truesc - A) {
            break;
        }
    }
    a.nm = g.nm;
    let mut cigar = g.cigar;
    let (mut pos, is_rev) = idx.bns.depos(if rb < idx.bns.l_pac { rb } else { re - 1 });
    a.is_rev = is_rev;
    if !cigar.is_empty() {
        // Squeeze out a leading or trailing deletion.
        if cigar[0] & 0xf == 2 {
            pos += (cigar[0] >> 4) as i64;
            cigar.remove(0);
        } else if cigar[cigar.len() - 1] & 0xf == 2 {
            cigar.pop();
        }
    }
    if qb != 0 || qe != l_query {
        let (clip5, clip3) = if is_rev {
            (l_query - qe, qb)
        } else {
            (qb, l_query - qe)
        };
        if clip5 != 0 {
            cigar.insert(0, (clip5 as u32) << 4 | 3);
        }
        if clip3 != 0 {
            cigar.push((clip3 as u32) << 4 | 3);
        }
    }
    a.cigar = cigar;
    a.md = g.md;
    a.rid = idx.bns.pos2rid(pos);
    assert_eq!(a.rid, ar.rid);
    a.pos = pos - idx.bns.anns[a.rid as usize].offset;
    a.score = ar.score;
    a.sub = ar.sub.max(ar.csub);
    a.is_alt = ar.is_alt;
    a.alt_sc = ar.alt_sc;
    a
}

/// `get_pri_idx`.
fn pri_idx(a: &[Region], i: usize) -> Option<usize> {
    let k = a[i].secondary_all;
    (k >= 0 && a[i].score as f64 >= a[k as usize].score as f64 * XA_DROP_RATIO as f64)
        .then_some(k as usize)
}

/// `mem_gen_alt`: per-region XA strings (only after `mark_primary_se`).
pub fn gen_alt(idx: &Index, a: &[Region], query: &[u8]) -> Option<Vec<Option<String>>> {
    let n = a.len();
    let mut cnt = vec![0i32; n];
    let mut has_alt = vec![false; n];
    let mut tot = 0;
    for i in 0..n {
        if let Some(r) = pri_idx(a, i) {
            cnt[r] += 1;
            tot += 1;
            if a[i].is_alt {
                has_alt[r] = true;
            }
        }
    }
    if tot == 0 {
        return None;
    }
    let mut aln: Vec<Option<String>> = vec![None; n];
    for i in 0..n {
        let Some(r) = pri_idx(a, i) else { continue };
        if cnt[r] > MAX_XA_HITS_ALT || (!has_alt[r] && cnt[r] > MAX_XA_HITS) {
            continue;
        }
        let t = reg2aln(idx, query, Some(&a[i]));
        let s = aln[r].get_or_insert_with(String::new);
        s.push_str(&idx.bns.anns[t.rid as usize].name);
        s.push(',');
        s.push(if t.is_rev { '-' } else { '+' });
        write!(s, "{},", t.pos + 1).unwrap();
        for &c in &t.cigar {
            write!(s, "{}{}", c >> 4, b"MIDSHN"[(c & 0xf) as usize] as char).unwrap();
        }
        write!(s, ",{};", t.nm).unwrap();
    }
    Some(aln)
}

fn push_cigar(out: &mut Vec<u8>, cigar: &[u32]) {
    // -Y: clipping is always soft, so ops print as stored.
    for &c in cigar {
        write!(VecWriter(out), "{}", c >> 4).unwrap();
        out.push(b"MIDSH"[(c & 0xf) as usize]);
    }
}

/// Lets `write!` append to a `Vec<u8>`.
struct VecWriter<'a>(&'a mut Vec<u8>);
impl std::fmt::Write for VecWriter<'_> {
    fn write_str(&mut self, s: &str) -> std::fmt::Result {
        self.0.extend_from_slice(s.as_bytes());
        Ok(())
    }
}

/// One reference's records within a batch's mapping: at `start` (page
/// aligned), `count` little-endian `u64` sort keys, then the `count`
/// length-prefixed BAM records in emission order; `len` bytes in all.
#[derive(Clone, Copy)]
pub(crate) struct BinSpan {
    pub(crate) start: usize,
    pub(crate) count: usize,
    pub(crate) len: usize,
}

/// A batch's records grouped by reference id (ascending), all in one
/// mapping of their own.
///
/// Batches are built on their own threads and freed on the sorter thread;
/// with heap bins the allocator kept such cross-thread frees resident (a
/// 30x chunk sort grew by 16 GB), so record memory bypasses it: each sorted
/// reference's pages go back to
/// the system at once ([`nomics_mem::PageBuf::release`]), and the mapping
/// is unmapped when the chunk is done.
pub(crate) struct BinnedRecords {
    pub(crate) buf: nomics_mem::PageBuf,
    pub(crate) bins: Vec<(usize, BinSpan)>,
    pub(crate) error: Option<io::Error>,
}

/// A chunk's records as its batches complete: every batch's mapping, and
/// for each reference id its bins in batch order. Appending a batch moves
/// it without copying records, and [`crate::samsort::sort_chunk`] frees
/// each reference's records once they are in the sorted stream, so a chunk
/// is never held twice in full.
#[derive(Default)]
pub struct ChunkRecords {
    pub(crate) batches: Vec<nomics_mem::PageBuf>,
    /// By reference id: (batch, span), in batch order.
    pub(crate) refs: Vec<Vec<(usize, BinSpan)>>,
    pub(crate) error: Option<io::Error>,
}

impl ChunkRecords {
    pub(crate) fn absorb(&mut self, batch: BinnedRecords) {
        if self.error.is_none() {
            self.error = batch.error;
        }
        let b = self.batches.len();
        self.batches.push(batch.buf);
        for (tid, span) in batch.bins {
            if self.refs.len() <= tid {
                self.refs.resize_with(tid + 1, Vec::new);
            }
            self.refs[tid].push((b, span));
        }
    }

    /// Bytes held: (records plus sort keys, mapped).
    pub(crate) fn held(&self) -> (usize, usize) {
        let used = self.refs.iter().flatten().map(|(_, s)| s.len).sum();
        (used, self.batches.iter().map(|b| b.len()).sum())
    }
}

#[derive(Default)]
pub struct BamRecords {
    pub(crate) bytes: Vec<u8>,
    pub(crate) keys: Vec<(u64, u64)>,
    error: Option<io::Error>,
}

impl BamRecords {
    fn push(&mut self, idx: &Index, record: Record) {
        if self.error.is_some() {
            return;
        }
        let tid = if record.ref_id < 0 {
            idx.bns.anns.len() as u64
        } else {
            record.ref_id as u64
        };
        let key = tid << 33 | ((record.pos as i64 + 1) as u64) << 1 | (record.flag >> 4 & 1) as u64;
        let off = self.bytes.len() as u64;
        self.bytes.extend_from_slice(&0u32.to_le_bytes());
        let start = self.bytes.len();
        if let Err(e) = record.encode(BinRule::Sam, &mut self.bytes) {
            self.error = Some(e);
            self.bytes.truncate(off as usize);
            return;
        }
        let len = (self.bytes.len() - start) as u32;
        self.bytes[off as usize..off as usize + 4].copy_from_slice(&len.to_le_bytes());
        self.keys.push((key, off));
    }

    pub(crate) fn append(&mut self, mut other: Self) {
        if self.error.is_none() {
            self.error = other.error.take();
        }
        let base = self.bytes.len() as u64;
        self.keys
            .extend(other.keys.into_iter().map(|(key, off)| (key, base + off)));
        self.bytes.append(&mut other.bytes);
    }

    /// Regroup by reference id (the sort key's top bits; unmapped is last)
    /// into one mapping ([`BinnedRecords`]), each reference's keys and
    /// records kept in emission order. Sizes are counted first, and every
    /// reference starts on a page, so the sort can free each one alone.
    pub(crate) fn into_bins(mut self) -> BinnedRecords {
        let tid = |key: u64| (key >> 33) as usize;
        let len = |off: usize| {
            4 + u32::from_le_bytes(self.bytes[off..off + 4].try_into().unwrap()) as usize
        };
        let n_tids = self.keys.iter().map(|&(k, _)| tid(k) + 1).max().unwrap_or(0);
        let mut sizes = vec![(0usize, 0usize); n_tids];
        for &(key, off) in &self.keys {
            let s = &mut sizes[tid(key)];
            s.0 += len(off as usize);
            s.1 += 1;
        }
        let page = nomics_mem::page_size();
        let mut slot = vec![usize::MAX; n_tids];
        let mut bins = Vec::new();
        let mut total = 0usize;
        for (t, &(bytes, count)) in sizes.iter().enumerate() {
            if count > 0 {
                slot[t] = bins.len();
                let span = BinSpan {
                    start: total,
                    count,
                    len: 8 * count + bytes,
                };
                total = (total + span.len).next_multiple_of(page);
                bins.push((t, span));
            }
        }
        let mut buf = match nomics_mem::PageBuf::new(total) {
            Ok(buf) => buf,
            Err(e) => {
                return BinnedRecords {
                    buf: nomics_mem::PageBuf::new(0).expect("empty mapping"),
                    bins: Vec::new(),
                    error: Some(self.error.take().unwrap_or(e)),
                };
            }
        };
        // Fill: (next key, next record) position of each reference.
        let mut fill: Vec<(usize, usize)> =
            bins.iter().map(|(_, s)| (s.start, s.start + 8 * s.count)).collect();
        let out = buf.as_mut_slice();
        for &(key, off) in &self.keys {
            let f = &mut fill[slot[tid(key)]];
            let (off, n) = (off as usize, len(off as usize));
            out[f.0..f.0 + 8].copy_from_slice(&key.to_le_bytes());
            out[f.1..f.1 + n].copy_from_slice(&self.bytes[off..off + n]);
            f.0 += 8;
            f.1 += n;
        }
        BinnedRecords {
            buf,
            bins,
            error: self.error.take(),
        }
    }

    pub(crate) fn take_error(&mut self) -> Option<io::Error> {
        self.error.take()
    }
}

pub trait AlnSink {
    fn emit(
        &mut self,
        idx: &Index,
        rg_id: &[u8],
        s: &Read,
        list: &[Aln],
        which: usize,
        mate: Option<&Aln>,
    );
}

impl AlnSink for Vec<u8> {
    fn emit(
        &mut self,
        idx: &Index,
        rg_id: &[u8],
        s: &Read,
        list: &[Aln],
        which: usize,
        mate: Option<&Aln>,
    ) {
        aln2sam(idx, rg_id, self, s, list, which, mate);
    }
}

impl AlnSink for BamRecords {
    fn emit(
        &mut self,
        idx: &Index,
        rg_id: &[u8],
        s: &Read,
        list: &[Aln],
        which: usize,
        mate: Option<&Aln>,
    ) {
        aln2bam(idx, rg_id, self, s, list, which, mate);
    }
}

fn prepare_aln(list: &[Aln], which: usize, m_: Option<&Aln>) -> (Aln, Option<Aln>, u16) {
    let mut p = list[which].clone();
    let mut m = m_.cloned();
    if m.is_some() {
        p.flag |= 0x1;
    }
    if p.rid < 0 {
        p.flag |= 0x4;
    }
    if let Some(m) = &m
        && m.rid < 0
    {
        p.flag |= 0x8;
    }
    if let Some(m) = &m
        && p.rid < 0
        && m.rid >= 0
    {
        p.rid = m.rid;
        p.pos = m.pos;
        p.is_rev = m.is_rev;
        p.cigar.clear();
    }
    if let Some(m) = &mut m
        && m.rid < 0
        && p.rid >= 0
    {
        m.rid = p.rid;
        m.pos = p.pos;
        m.is_rev = p.is_rev;
        m.cigar.clear();
    }
    if p.is_rev {
        p.flag |= 0x10;
    }
    if m.as_ref().is_some_and(|m| m.is_rev) {
        p.flag |= 0x20;
    }
    let flag = ((p.flag & 0xffff) | if p.flag & 0x10000 != 0 { 0x100 } else { 0 }) as u16;
    (p, m, flag)
}

/// `mem_aln2sam`: append the SAM line of `list[which]` to `out`.
pub fn aln2sam(
    idx: &Index,
    rg_id: &[u8],
    out: &mut Vec<u8>,
    s: &Read,
    list: &[Aln],
    which: usize,
    m_: Option<&Aln>,
) {
    let (p, m, flag) = prepare_aln(list, which, m_);

    let mut w = VecWriter(out);
    w.0.extend_from_slice(&s.name);
    write!(w, "\t{flag}\t").unwrap();
    if p.rid >= 0 {
        write!(
            w,
            "{}\t{}\t{}\t",
            idx.bns.anns[p.rid as usize].name,
            p.pos + 1,
            p.mapq
        )
        .unwrap();
        if p.cigar.is_empty() {
            w.0.push(b'*');
        } else {
            push_cigar(w.0, &p.cigar);
        }
    } else {
        w.0.extend_from_slice(b"*\t0\t0\t*");
    }
    w.0.push(b'\t');

    let mut w = VecWriter(out);
    match &m {
        Some(m) if m.rid >= 0 => {
            if p.rid == m.rid {
                w.0.push(b'=');
            } else {
                w.0.extend_from_slice(idx.bns.anns[m.rid as usize].name.as_bytes());
            }
            write!(w, "\t{}\t", m.pos + 1).unwrap();
            if p.rid == m.rid {
                let p0 = p.pos + if p.is_rev { get_rlen(&p.cigar) - 1 } else { 0 };
                let p1 = m.pos + if m.is_rev { get_rlen(&m.cigar) - 1 } else { 0 };
                if m.cigar.is_empty() || p.cigar.is_empty() {
                    w.0.push(b'0');
                } else {
                    let d = p0 - p1 + (p0 > p1) as i64 - (p0 < p1) as i64;
                    write!(w, "{}", -d).unwrap();
                }
            } else {
                w.0.push(b'0');
            }
        }
        _ => w.0.extend_from_slice(b"*\t0\t0"),
    }
    out.push(b'\t');

    // SEQ and QUAL (no secondaries are printed, and -Y keeps full SEQ on
    // supplementaries).
    if p.flag & 0x100 != 0 {
        out.extend_from_slice(b"*\t*");
    } else if !p.is_rev {
        out.extend(s.seq.iter().map(|&c| b"ACGTN"[c as usize]));
        out.push(b'\t');
        out.extend_from_slice(&s.qual);
    } else {
        out.extend(s.seq.iter().rev().map(|&c| b"TGCAN"[c as usize]));
        out.push(b'\t');
        out.extend(s.qual.iter().rev());
    }

    let mut w = VecWriter(out);
    if !p.cigar.is_empty() {
        write!(w, "\tNM:i:{}\tMD:Z:{}", p.nm, p.md).unwrap();
    }
    if let Some(m) = &m
        && !m.cigar.is_empty()
    {
        w.0.extend_from_slice(b"\tMC:Z:");
        push_cigar(w.0, &m.cigar);
    }
    let mut w = VecWriter(out);
    if let Some(m) = &m {
        write!(w, "\tMQ:i:{}", m.mapq).unwrap();
    }
    if p.score >= 0 {
        write!(w, "\tAS:i:{}", p.score).unwrap();
    }
    if p.sub >= 0 {
        write!(w, "\tXS:i:{}", p.sub).unwrap();
    }
    if !rg_id.is_empty() {
        w.0.extend_from_slice(b"\tRG:Z:");
        w.0.extend_from_slice(rg_id);
    }
    if p.flag & 0x100 == 0 {
        let others = (0..list.len()).any(|i| i != which && list[i].flag & 0x100 == 0);
        if others {
            w.0.extend_from_slice(b"\tSA:Z:");
            for (i, r) in list.iter().enumerate() {
                if i == which || r.flag & 0x100 != 0 {
                    continue;
                }
                let mut w = VecWriter(out);
                write!(
                    w,
                    "{},{},{},",
                    idx.bns.anns[r.rid as usize].name,
                    r.pos + 1,
                    if r.is_rev { '-' } else { '+' }
                )
                .unwrap();
                push_cigar(out, &r.cigar);
                write!(VecWriter(out), ",{},{};", r.mapq, r.nm).unwrap();
            }
        }
        if p.alt_sc > 0 {
            write!(
                VecWriter(out),
                "\tpa:f:{:.3}",
                p.score as f64 / p.alt_sc as f64
            )
            .unwrap();
        }
    }
    if let Some(xa) = &p.xa {
        out.extend_from_slice(b"\tXA:Z:");
        out.extend_from_slice(xa.as_bytes());
    }
    out.push(b'\n');
}

fn int_value(v: i64) -> Value {
    let ty = if v < 0 {
        if v >= -128 {
            IntType::I8
        } else if v >= -32768 {
            IntType::I16
        } else {
            IntType::I32
        }
    } else if v <= 255 {
        IntType::U8
    } else if v <= 65535 {
        IntType::U16
    } else {
        IntType::U32
    };
    Value::Int(ty, v)
}

fn cigar_text(cigar: &[u32]) -> Vec<u8> {
    let mut text = Vec::new();
    push_cigar(&mut text, cigar);
    text
}

fn aln2bam(
    idx: &Index,
    rg_id: &[u8],
    out: &mut BamRecords,
    s: &Read,
    list: &[Aln],
    which: usize,
    m_: Option<&Aln>,
) {
    let (p, m, mut flag) = prepare_aln(list, which, m_);
    if p.rid < 0 || p.cigar.is_empty() {
        flag |= 0x4;
    }
    let (next_ref_id, next_pos, tlen) = match &m {
        Some(m) if m.rid >= 0 => {
            let tlen = if p.rid == m.rid && !m.cigar.is_empty() && !p.cigar.is_empty() {
                let p0 = p.pos + if p.is_rev { get_rlen(&p.cigar) - 1 } else { 0 };
                let p1 = m.pos + if m.is_rev { get_rlen(&m.cigar) - 1 } else { 0 };
                let d = p0 - p1 + (p0 > p1) as i64 - (p0 < p1) as i64;
                -d as i32
            } else {
                0
            };
            (m.rid, m.pos as i32, tlen)
        }
        _ => (-1, -1, 0),
    };
    let (seq, qual) = if p.flag & 0x100 != 0 {
        (b"*".to_vec(), vec![0xff])
    } else if !p.is_rev {
        (
            s.seq.iter().map(|&c| b"ACGTN"[c as usize]).collect(),
            s.qual.iter().map(|&q| q - 33).collect(),
        )
    } else {
        (
            s.seq.iter().rev().map(|&c| b"TGCAN"[c as usize]).collect(),
            s.qual.iter().rev().map(|&q| q - 33).collect(),
        )
    };
    let mut tags = Vec::new();
    if !p.cigar.is_empty() {
        tags.push((*b"NM", int_value(p.nm as i64)));
        tags.push((*b"MD", Value::String(p.md.as_bytes().to_vec())));
    }
    if let Some(m) = &m
        && !m.cigar.is_empty()
    {
        tags.push((*b"MC", Value::String(cigar_text(&m.cigar))));
    }
    if let Some(m) = &m {
        tags.push((*b"MQ", int_value(m.mapq as i64)));
    }
    if p.score >= 0 {
        tags.push((*b"AS", int_value(p.score as i64)));
    }
    if p.sub >= 0 {
        tags.push((*b"XS", int_value(p.sub as i64)));
    }
    if !rg_id.is_empty() {
        tags.push((*b"RG", Value::String(rg_id.to_vec())));
    }
    if p.flag & 0x100 == 0 {
        let others = (0..list.len()).any(|i| i != which && list[i].flag & 0x100 == 0);
        if others {
            let mut sa = Vec::new();
            for (i, r) in list.iter().enumerate() {
                if i == which || r.flag & 0x100 != 0 {
                    continue;
                }
                write!(
                    VecWriter(&mut sa),
                    "{},{},{},",
                    idx.bns.anns[r.rid as usize].name,
                    r.pos + 1,
                    if r.is_rev { '-' } else { '+' }
                )
                .unwrap();
                push_cigar(&mut sa, &r.cigar);
                write!(VecWriter(&mut sa), ",{},{};", r.mapq, r.nm).unwrap();
            }
            tags.push((*b"SA", Value::String(sa)));
        }
        if p.alt_sc > 0 {
            let text = format!("{:.3}", p.score as f64 / p.alt_sc as f64);
            tags.push((*b"pa", Value::Float(text.parse::<f32>().unwrap().to_bits())));
        }
    }
    if let Some(xa) = &p.xa {
        tags.push((*b"XA", Value::String(xa.as_bytes().to_vec())));
    }
    out.push(
        idx,
        Record {
            ref_id: p.rid,
            pos: p.pos as i32,
            mapq: p.mapq,
            flag,
            next_ref_id,
            next_pos,
            tlen,
            read_name: s.name.clone(),
            // bwa's internal ops are MIDSH; BAM's are MIDNSHP=X.
            cigar: p
                .cigar
                .iter()
                .map(|&c| (c & !0xf) | [0, 1, 2, 4, 5][(c & 0xf) as usize])
                .collect(),
            seq,
            qual,
            tags,
        },
    );
}

/// `mem_reg2sam`: all lines of one read on the unpaired path.
pub fn reg2out<S: AlnSink>(
    idx: &Index,
    rg_id: &[u8],
    out: &mut S,
    s: &Read,
    a: &[Region],
    extra_flag: i32,
    m: Option<&Aln>,
) {
    let xa = gen_alt(idx, a, &s.seq);
    let mut aa: Vec<Aln> = Vec::new();
    let mut l = 0;
    for (k, p) in a.iter().enumerate() {
        if p.score < T {
            continue;
        }
        if p.secondary >= 0 {
            continue; // not -a: secondaries are never printed
        }
        let mut q = reg2aln(idx, &s.seq, Some(p));
        assert!(q.rid >= 0);
        q.xa = xa.as_ref().and_then(|x| x[k].clone());
        q.flag |= extra_flag;
        if l > 0 && p.secondary < 0 {
            q.flag |= 0x800;
        }
        if l > 0 && !p.is_alt && q.mapq > aa[0].mapq {
            q.mapq = aa[0].mapq;
        }
        aa.push(q);
        l += 1;
    }
    if aa.is_empty() {
        let mut t = reg2aln(idx, &s.seq, None);
        t.flag |= extra_flag;
        out.emit(idx, rg_id, s, std::slice::from_ref(&t), 0, m);
    } else {
        for k in 0..aa.len() {
            out.emit(idx, rg_id, s, &aa, k, m);
        }
    }
}

pub fn reg2sam(
    idx: &Index,
    rg_id: &[u8],
    out: &mut Vec<u8>,
    s: &Read,
    a: &[Region],
    extra_flag: i32,
    m: Option<&Aln>,
) {
    reg2out(idx, rg_id, out, s, a, extra_flag, m);
}
