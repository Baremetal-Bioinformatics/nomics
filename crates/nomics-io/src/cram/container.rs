//! One data container (one slice) from a run of BAM records.
//!
//! Layout choices (all within CRAM 3.0; see `mod.rs` for why each is lossless):
//!
//! - One slice per container, one reference per slice (`-1` for the unplaced
//!   tail). No multi-reference slices, no embedded reference.
//! - Every data series is EXTERNAL in its own block, except integer series
//!   whose value is constant in the slice, which are a zero-bit HUFFMAN code
//!   (no bytes at all). The core block is empty.
//! - Tags are stored verbatim with their BAM type codes and in their BAM
//!   order, including MD, NM and RG. The RG data series is constant -1.
//! - Read names are always kept. Mates are "attached" (CF 0x4 + NF) only when
//!   htslib's, htsjdk 4's and htsjdk 5's decoders would all rebuild the mate fields and
//!   TLEN exactly; every other paired record is detached with MF/NS/NP/TS.
//! - Quality scores are always stored as an array (CF 0x1).
//!
//! Records the format cannot carry losslessly in this layout (unmapped reads
//! with a CIGAR or MAPQ, `=`/`X` CIGAR ops, mapped reads without SEQ, CIGARs
//! the decoders would rewrite, unpaired records with mate fields, flags of
//! 0x1000 and up) are rejected with an error rather than
//! silently changed.

use std::collections::HashMap;
use std::io;

use super::itf8::{put_itf8, put_ltf8};
use super::md5::Md5;
use super::rans;
use crate::crc32::crc32;
use crate::deflate::FastDeflater;

fn invalid(msg: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg.into())
}

/// A value CRAM stores as a signed 32-bit integer (ITF-8, int32), checked.
fn i32c<T: TryInto<i32> + Copy + std::fmt::Display>(v: T, what: &str) -> io::Result<i32> {
    v.try_into()
        .map_err(|_| invalid(format!("{what} {v} does not fit CRAM's 32-bit field")))
}

/// Block content types.
const FILE_HEADER: u8 = 0;
const COMPRESSION_HEADER: u8 = 1;
const SLICE_HEADER: u8 = 2;
const EXTERNAL_DATA: u8 = 4;
const CORE_DATA: u8 = 5;

/// Block compression methods.
const RAW: u8 = 0;
const GZIP: u8 = 1;
const RANS: u8 = 4;

/// Codec ids.
const E_EXTERNAL: i32 = 1;
const E_HUFFMAN: i32 = 3;
const E_BYTE_ARRAY_LEN: i32 = 4;
const E_BYTE_ARRAY_STOP: i32 = 5;

/// CRAM record flags (CF).
const CF_QS_ARRAY: i32 = 0x1;
const CF_DETACHED: i32 = 0x2;
const CF_MATE_DOWNSTREAM: i32 = 0x4;

// BAM flags.
const PAIRED: u16 = 0x1;
const UNMAPPED: u16 = 0x4;
const MATE_UNMAPPED: u16 = 0x8;
const REVERSE: u16 = 0x10;
const MATE_REVERSE: u16 = 0x20;
const READ1: u16 = 0x40;
const SECONDARY: u16 = 0x100;
const SUPPLEMENTARY: u16 = 0x800;

/// Integer data series, in our block content id order (id = index + 1).
#[derive(Clone, Copy)]
enum IntDs {
    Bf,
    Cf,
    Rl,
    Ap,
    Rg,
    Mf,
    Ns,
    Np,
    Ts,
    Nf,
    Tl,
    Fn,
    Fp,
    Dl,
    Rs,
    Pd,
    Hc,
    Mq,
}
const INT_DS: [(&[u8; 2], IntDs); 18] = [
    (b"BF", IntDs::Bf),
    (b"CF", IntDs::Cf),
    (b"RL", IntDs::Rl),
    (b"AP", IntDs::Ap),
    (b"RG", IntDs::Rg),
    (b"MF", IntDs::Mf),
    (b"NS", IntDs::Ns),
    (b"NP", IntDs::Np),
    (b"TS", IntDs::Ts),
    (b"NF", IntDs::Nf),
    (b"TL", IntDs::Tl),
    (b"FN", IntDs::Fn),
    (b"FP", IntDs::Fp),
    (b"DL", IntDs::Dl),
    (b"RS", IntDs::Rs),
    (b"PD", IntDs::Pd),
    (b"HC", IntDs::Hc),
    (b"MQ", IntDs::Mq),
];

/// Byte series (one byte per value).
#[derive(Clone, Copy)]
enum ByteDs {
    Fc,
    Bs,
    Ba,
    Qs,
}
const BYTE_DS: [(&[u8; 2], ByteDs); 4] = [
    (b"FC", ByteDs::Fc),
    (b"BS", ByteDs::Bs),
    (b"BA", ByteDs::Ba),
    (b"QS", ByteDs::Qs),
];

/// Byte-array series, BYTE_ARRAY_STOP with a NUL stop byte.
#[derive(Clone, Copy)]
enum ArrDs {
    Rn,
    In,
    Sc,
}
const ARR_DS: [(&[u8; 2], ArrDs); 3] = [(b"RN", ArrDs::Rn), (b"IN", ArrDs::In), (b"SC", ArrDs::Sc)];

const INT_ID0: i32 = 1;
const BYTE_ID0: i32 = INT_ID0 + INT_DS.len() as i32;
const ARR_ID0: i32 = BYTE_ID0 + BYTE_DS.len() as i32;

/// Upper-case base → substitution matrix index (A C G T N), or `None`.
#[inline]
fn acgtn(b: u8) -> Option<usize> {
    match b {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        b'N' => Some(4),
        _ => None,
    }
}

const SEQ_CODES: &[u8; 16] = b"=ACMGRSVTWYHKDBN";

/// How tags are stored, which decides what CRAM readers decode.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TagMode {
    /// Every tag verbatim, in BAM order; RG series -1; `cF:C` on mapped
    /// reads lacking MD/NM. Decodes to the BAM's records.
    Lossless,
    /// As htslib's `samtools view -C` (sarek's md.cram): RG:Z of a header
    /// read group becomes the RG series; MD and NM are dropped when htslib
    /// would regenerate the same values; other tags keep BAM order.
    Htslib,
    /// As htsjdk 4's CRAM writer fed by an htsjdk read of an htslib CRAM
    /// (sarek's recal.cram, ApplyBQSR writing CRAM from md.cram): the same
    /// drops as [`TagMode::Htslib`], remaining tags sorted by the 3-byte
    /// tag key as htsjdk's `CompressionHeaderFactory` does.
    HtsjdkSorted,
    /// [`TagMode::HtsjdkSorted`], then as `samtools merge` (1.24) of such
    /// CRAMs leaves them (sarek's multi-group recal.cram): `bam_translate`
    /// deletes and re-appends RG and PG on every record, so a stored PG:Z
    /// (and an RG:Z not in the RG series) follows the other tags.
    HtsjdkSortedMerged,
}

/// Per-file encoding parameters shared by every container.
#[derive(Debug, Clone)]
pub struct EncodeParams {
    pub tags: TagMode,
    /// `@RG` IDs in header order (the RG data series indexes them).
    pub read_groups: Vec<Vec<u8>>,
}

/// Integer value of a BAM integer tag.
fn int_tag(ty: u8, v: &[u8]) -> Option<i64> {
    Some(match (ty, v.len()) {
        (b'c', 1) => v[0] as i8 as i64,
        (b'C', 1) => v[0] as i64,
        (b's', 2) => i16::from_le_bytes([v[0], v[1]]) as i64,
        (b'S', 2) => u16::from_le_bytes([v[0], v[1]]) as i64,
        (b'i', 4) => i32::from_le_bytes(v.try_into().unwrap()) as i64,
        (b'I', 4) => u32::from_le_bytes(v.try_into().unwrap()) as i64,
        _ => return None,
    })
}

/// MD and NM as htslib's encoder computes them (`process_one_read`, htslib
/// 1.21-1.24) against the upper-cased reference, and whether it would keep
/// both verbatim regardless (an `N` read base on an `N` reference base, or
/// an alignment running off the contig).
fn htslib_md_nm(r: &Rec, reference: &[u8]) -> (Vec<u8>, i64, bool) {
    use std::io::Write as _;
    let ref_end = reference.len() as i64;
    let rb = |p: i64| reference[p as usize].to_ascii_uppercase();
    let mut md = Vec::new();
    let (mut nm, mut verbatim) = (0i64, false);
    let mut apos = r.apos - 1;
    let mut spos = 0usize;
    let mut md_last = apos;
    for c in r.cigar() {
        let (op, len) = (c & 0xf, (c >> 4) as i64);
        match op {
            0 | 7 | 8 => {
                let end = if len + apos < ref_end {
                    len
                } else {
                    ref_end - apos
                };
                let mut l = 0i64;
                while l < end {
                    let (rp, sp) = (rb(apos + l), r.base(spos + l as usize));
                    if rp == b'N' && sp == b'N' {
                        verbatim = true;
                    }
                    if rp != sp {
                        let _ = write!(md, "{}", apos + l - md_last);
                        md.push(rp);
                        md_last = apos + l + 1;
                        nm += 1;
                    }
                    l += 1;
                }
                let l = l.max(0);
                spos += l as usize;
                apos += l;
                if l < len {
                    verbatim = true;
                    spos += (len - l) as usize;
                    apos += len; // as htslib (it adds the whole length here)
                }
            }
            2 => {
                let _ = write!(md, "{}", apos - md_last);
                if apos < ref_end {
                    md.push(b'^');
                    for p in apos..(apos + len).min(ref_end) {
                        md.push(rb(p));
                    }
                }
                nm += len;
                apos += len;
                md_last = apos;
            }
            3 => {
                apos += len;
                md_last += len;
            }
            1 => {
                nm += len;
                spos += len as usize;
            }
            4 => spos += len as usize,
            _ => {}
        }
    }
    let _ = write!(md, "{}", apos - md_last);
    (md, nm, verbatim)
}

/// The reference a slice's records are aligned to.
pub struct SliceRef<'a> {
    /// Contig bases (any case), or `None` for the unplaced (`-1`) slice.
    pub bases: Option<&'a [u8]>,
}

/// One container's worth of raw BAM records (without `block_size`), all on
/// reference `ref_id`.
pub struct ContainerInput {
    pub data: Vec<u8>,
    pub ends: Vec<usize>,
    pub ref_id: i32,
    /// File-wide index of the first record.
    pub record_counter: u64,
}

impl ContainerInput {
    pub fn len(&self) -> usize {
        self.ends.len()
    }

    pub fn is_empty(&self) -> bool {
        self.ends.is_empty()
    }

    fn record(&self, i: usize) -> &[u8] {
        let s = if i == 0 { 0 } else { self.ends[i - 1] };
        &self.data[s..self.ends[i]]
    }
}

/// A block compression a trial can choose.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Codec {
    Raw,
    Rans0,
    Rans1,
    /// gzip with our fast deflater.
    Gzip,
    /// gzip with zlib-exact level 5 (htslib's default CRAM level).
    Zlib5,
}

/// The codec each external block (by content id) is compressed with,
/// decided by a trial container from its own data ([`Encoder::encode`]
/// with no plan). Blocks with an id missing from the plan get a full trial.
pub type CodecPlan = HashMap<i32, Codec>;

/// An encoded container plus what the `.crai` needs.
pub struct Encoded {
    pub bytes: Vec<u8>,
    pub ref_id: i32,
    pub start: i64,
    pub span: i64,
    /// Slice header offset from the end of the container header.
    pub landmark: usize,
    /// Slice header block plus all its data blocks.
    pub slice_size: usize,
    pub records: usize,
    /// Each external block's content id, codec and uncompressed length.
    pub codecs: Vec<(i32, Codec, usize)>,
}

/// Fixed fields of a parsed record (offsets into the raw bytes).
struct Rec<'a> {
    raw: &'a [u8],
    ref_id: i32,
    /// 1-based alignment start (0 for unplaced).
    apos: i64,
    /// 1-based inclusive alignment end (apos - 1 + reference span; apos for unmapped).
    aend: i64,
    flag: u16,
    mapq: u8,
    next_ref: i32,
    next_pos: i32,
    tlen: i32,
    name: &'a [u8],
    /// Offset of the CIGAR ops: in the record, or in its `CG:B:I` tag.
    cigar_at: usize,
    n_cigar: usize,
    /// Offset in `aux()` of a `CG` tag that holds the real CIGAR (dropped).
    cg_at: Option<usize>,
    seq_at: usize,
    l_seq: usize,
    qual_at: usize,
    aux_at: usize,
}

impl<'a> Rec<'a> {
    fn parse(raw: &'a [u8]) -> io::Result<Rec<'a>> {
        if raw.len() < 32 {
            return Err(invalid("BAM record shorter than its fixed fields"));
        }
        let i32_at = |at: usize| i32::from_le_bytes(raw[at..at + 4].try_into().unwrap());
        let l_read_name = raw[8] as usize;
        let n_cigar = u16::from_le_bytes([raw[12], raw[13]]) as usize;
        let flag = u16::from_le_bytes([raw[14], raw[15]]);
        let l_seq = u32::from_le_bytes(raw[16..20].try_into().unwrap()) as usize;
        let cigar_at = 32 + l_read_name;
        let seq_at = cigar_at + 4 * n_cigar;
        let qual_at = seq_at + l_seq.div_ceil(2);
        let aux_at = qual_at + l_seq;
        if aux_at > raw.len() || l_read_name == 0 {
            return Err(invalid("truncated BAM record"));
        }
        let pos = i32_at(4);
        let mut r = Rec {
            raw,
            ref_id: i32_at(0),
            apos: pos as i64 + 1,
            aend: 0,
            flag,
            mapq: raw[9],
            next_ref: i32_at(20),
            next_pos: i32_at(24),
            tlen: i32_at(28),
            name: &raw[32..32 + l_read_name - 1],
            cigar_at,
            n_cigar,
            cg_at: None,
            seq_at,
            l_seq,
            qual_at,
            aux_at,
        };
        r.take_cg_cigar()?;
        r.aend = if r.mapped() {
            r.apos - 1
                + r.cigar()
                    .filter(|c| matches!(c & 0xf, 0 | 2 | 3 | 7 | 8))
                    .map(|c| (c >> 4) as i64)
                    .sum::<i64>()
        } else {
            r.apos
        };
        Ok(r)
    }

    fn mapped(&self) -> bool {
        self.flag & UNMAPPED == 0
    }

    /// A BAM record with more than 65,535 CIGAR ops stores `<l_seq>S<span>N`
    /// and the real CIGAR in `CG:B:I`. Ported from htslib 1.24
    /// `bam_tag2cigar` (applied when samtools reads a BAM, so `samtools view
    /// -C` writes the real CIGAR): on a placed record (tid >= 0, pos >= 0)
    /// with a CIGAR whose first op soft-clips the whole read, the first `CG`
    /// tag, if `B:I`/`B:i` with n_cigar <= length < 2^29, is the CIGAR and
    /// is dropped. Anything else leaves the record as it is.
    fn take_cg_cigar(&mut self) -> io::Result<()> {
        if self.n_cigar == 0 || self.ref_id < 0 || self.apos < 1 {
            return Ok(());
        }
        let first = self.cigar().next().unwrap();
        if first & 0xf != 4 || (first >> 4) as usize != self.l_seq {
            return Ok(());
        }
        let aux = self.aux();
        let mut at = 0;
        while aux.len() - at >= 3 {
            let (name, ty, body) = ([aux[at], aux[at + 1]], aux[at + 2], &aux[at + 3..]);
            let len = crate::bam::aux_value_len(ty, body)
                .filter(|&l| l <= body.len())
                .ok_or_else(|| invalid("malformed tag"))?;
            if &name == b"CG" {
                // bam_aux_get: the first CG tag decides. It must be B:I or
                // B:i, and hold at least as many ops as the placeholder (and
                // fewer than 2^29), else the record is left untouched.
                if ty != b'B' || !(body[0] == b'I' || body[0] == b'i') {
                    return Ok(());
                }
                let cg_len = u32::from_le_bytes(body[1..5].try_into().unwrap());
                if (cg_len as usize) < self.n_cigar || cg_len >= 1 << 29 {
                    return Ok(());
                }
                self.cigar_at = self.aux_at + at + 3 + 5;
                self.n_cigar = cg_len as usize;
                self.cg_at = Some(at);
                return Ok(());
            }
            at += 3 + len;
        }
        Ok(())
    }

    /// What every later step relies on: a placed mapped record whose CIGAR
    /// consumes exactly its bases, within i32 coordinates.
    fn validate_mapped(&self) -> io::Result<()> {
        let name = String::from_utf8_lossy(self.name);
        if self.l_seq == 0 {
            return Err(invalid(format!(
                "{name}: mapped record without SEQ is not supported"
            )));
        }
        if self.apos < 1 {
            return Err(invalid(format!("{name}: mapped record without a position")));
        }
        let query: u64 = self
            .cigar()
            .filter(|c| matches!(c & 0xf, 0 | 1 | 4 | 7 | 8))
            .map(|c| (c >> 4) as u64)
            .sum();
        if query != self.l_seq as u64 {
            return Err(invalid(format!(
                "{name}: CIGAR consumes {query} bases but SEQ has {}",
                self.l_seq
            )));
        }
        if self.aend > i32::MAX as i64 {
            return Err(invalid(format!("{name}: alignment end past 2^31")));
        }
        Ok(())
    }

    fn cigar(&self) -> impl Iterator<Item = u32> + '_ {
        (0..self.n_cigar).map(move |i| {
            let at = self.cigar_at + 4 * i;
            u32::from_le_bytes(self.raw[at..at + 4].try_into().unwrap())
        })
    }

    #[inline]
    fn base(&self, i: usize) -> u8 {
        let b = self.raw[self.seq_at + i / 2];
        SEQ_CODES[(if i.is_multiple_of(2) { b >> 4 } else { b & 0xf }) as usize]
    }

    fn quals(&self) -> &'a [u8] {
        &self.raw[self.qual_at..self.qual_at + self.l_seq]
    }

    fn aux(&self) -> &'a [u8] {
        &self.raw[self.aux_at..]
    }

    fn primary(&self) -> bool {
        self.flag & (SECONDARY | SUPPLEMENTARY) == 0
    }
}

/// TLEN as htslib's `cram_decode_slice_xref` rebuilds it for an attached
/// pair (`a` upstream, `b` downstream, same reference, both mapped).
fn htslib_tlen(a: &Rec, b: &Rec) -> (i64, i64) {
    let (mut aleft, mut aright) = (a.apos, a.aend);
    let (mut left_cnt, mut right_cnt) = (0, 0);
    for r in [a, b] {
        if aleft > r.apos {
            aleft = r.apos;
            left_cnt = 1;
        } else if aleft == r.apos {
            left_cnt += 1;
        }
        if aright < r.aend {
            aright = r.aend;
            right_cnt = 1;
        } else if aright == r.aend {
            right_cnt += 1;
        }
    }
    let tlen = aright - aleft + 1;
    if a.apos == aleft && (a.aend < aright || left_cnt <= 1) {
        (tlen, -tlen)
    } else if a.apos == aleft && a.aend == aright && left_cnt > 1 && right_cnt > 1 {
        if a.flag & READ1 != 0 {
            (tlen, -tlen)
        } else {
            (-tlen, tlen)
        }
    } else {
        (-tlen, tlen)
    }
}

/// Upstream TLEN as htsjdk 2.x-4.x rebuild it (`CRAMCompressionRecord.
/// computeInsertSize`, 5'-to-5'); the downstream record gets its negation.
/// GATK 4.6 bundles htsjdk 4.
fn htsjdk4_tlen(first: &Rec, second: &Rec) -> i64 {
    let five = |r: &Rec| {
        if r.flag & REVERSE != 0 {
            r.aend
        } else {
            r.apos
        }
    };
    let (f, s) = (five(first), five(second));
    s - f + if s >= f { 1 } else { -1 }
}

/// Upstream TLEN as htsjdk 5.0 rebuilds it: leftmost/rightmost span, signed
/// by start, ties broken by the upstream record's first-of-pair flag (unlike
/// htslib, which breaks start ties by end). Downstream gets the negation.
fn htsjdk5_tlen(first: &Rec, second: &Rec) -> i64 {
    let m = first.aend.max(second.aend) - first.apos.min(second.apos) + 1;
    if first.apos < second.apos || (first.apos == second.apos && first.flag & READ1 != 0) {
        m
    } else {
        -m
    }
}

/// May `a` (upstream) and `b` (downstream) be stored as an attached pair?
fn attachable(a: &Rec, b: &Rec) -> bool {
    let consistent = |x: &Rec, y: &Rec| {
        x.flag & PAIRED != 0
            && x.next_ref == y.ref_id
            && x.next_pos as i64 + 1 == y.apos
            && (x.flag & MATE_REVERSE != 0) == (y.flag & REVERSE != 0)
            && x.flag & MATE_UNMAPPED == 0
    };
    if !(a.mapped() && b.mapped() && a.ref_id == b.ref_id && a.ref_id >= 0) {
        return false;
    }
    if !(consistent(a, b) && consistent(b, a)) {
        return false;
    }
    // Every decoder must rebuild both TLENs; they disagree on ties and on
    // reads whose 5' ends are not their outer ends, so those stay detached.
    let (ta, tb) = htslib_tlen(a, b);
    let (j4, j5) = (htsjdk4_tlen(a, b), htsjdk5_tlen(a, b));
    let (at, bt) = (a.tlen as i64, b.tlen as i64);
    ta == at && tb == bt && j4 == at && -j4 == bt && j5 == at && -j5 == bt
}

/// Values of one tag (key = name and BAM type) across the slice.
struct TagSeries {
    key: i32,
    ty: u8,
    data: Vec<u8>,
    /// A `Z`/`H` value contains a tab, so the tab stop byte cannot be used.
    has_tab: bool,
}

/// Per-thread reusable state.
pub struct Encoder {
    deflater: FastDeflater,
    zlib: crate::deflate::ZlibDeflater,
}

impl Default for Encoder {
    fn default() -> Self {
        Self::new()
    }
}

/// Encoded values of every series in the slice.
struct Series {
    ints: Vec<Vec<i32>>,
    bytes: Vec<Vec<u8>>,
    arrs: Vec<Vec<u8>>,
    /// Substitutions as (reference index, read index); codes assigned later.
    subs: Vec<(u8, u8)>,
    /// Positions in `bytes[Bs]` matching `subs`.
    tags: Vec<TagSeries>,
    tag_index: HashMap<i32, usize>,
    td: Vec<Vec<u8>>,
    td_index: HashMap<Vec<u8>, i32>,
}

impl Series {
    fn new() -> Series {
        Series {
            ints: vec![Vec::new(); INT_DS.len()],
            bytes: vec![Vec::new(); BYTE_DS.len()],
            arrs: vec![Vec::new(); ARR_DS.len()],
            subs: Vec::new(),
            tags: Vec::new(),
            tag_index: HashMap::new(),
            td: Vec::new(),
            td_index: HashMap::new(),
        }
    }

    #[inline]
    fn int(&mut self, ds: IntDs, v: i32) {
        self.ints[ds as usize].push(v);
    }

    #[inline]
    fn byte(&mut self, ds: ByteDs, v: u8) {
        self.bytes[ds as usize].push(v);
    }

    fn arr(&mut self, ds: ArrDs, v: &[u8]) {
        let a = &mut self.arrs[ds as usize];
        a.extend_from_slice(v);
        a.push(0);
    }

    fn tag(&mut self, name: [u8; 2], ty: u8) -> &mut TagSeries {
        let key = (name[0] as i32) << 16 | (name[1] as i32) << 8 | ty as i32;
        let n = self.tags.len();
        let i = *self.tag_index.entry(key).or_insert(n);
        if i == n {
            self.tags.push(TagSeries {
                key,
                ty,
                data: Vec::new(),
                has_tab: false,
            });
        }
        &mut self.tags[i]
    }
}

/// A read feature being collected for one record.
enum Feature {
    Sub(u8, u8),
    Base(u8, u8),
    Ins(Vec<u8>),
    Del(i32),
    Skip(i32),
    Soft(Vec<u8>),
    Hard(i32),
    Pad(i32),
}

impl Encoder {
    pub fn new() -> Encoder {
        Encoder {
            deflater: FastDeflater::new(),
            zlib: crate::deflate::ZlibDeflater::new(5),
        }
    }

    /// Compress one block's data with the smallest of raw, gzip and rANS
    /// order 0/1 (a full trial). Returns (method, bytes, codec).
    fn compress(&mut self, data: &[u8]) -> (u8, Vec<u8>, Codec) {
        if data.len() < 16 {
            return (RAW, data.to_vec(), Codec::Raw);
        }
        let mut best = (RAW, data.to_vec(), Codec::Raw);
        let r0 = rans::encode_o0(data);
        if r0.len() < best.1.len() {
            best = (RANS, r0, Codec::Rans0);
        }
        let r1 = rans::encode_o1(data);
        if r1.len() < best.1.len() {
            best = (RANS, r1, Codec::Rans1);
        }
        let gz = gzip(&mut self.deflater, data);
        // Our greedy deflater is a quick probe; when gzip is in contention,
        // zlib level 5 (htslib's default CRAM level) usually beats it by
        // 5-15% on qualities and names.
        if gz.len() * 10 < best.1.len() * 11 {
            let z = gzip_zlib(&mut self.zlib, data);
            if z.len() < gz.len() && z.len() < best.1.len() {
                return (GZIP, z, Codec::Zlib5);
            }
        }
        if gz.len() < best.1.len() {
            best = (GZIP, gz, Codec::Gzip);
        }
        best
    }

    /// Compress one block's data with `codec` alone (raw if that is not
    /// smaller). Returns (method, bytes, codec).
    fn compress_with(&mut self, data: &[u8], codec: Codec) -> (u8, Vec<u8>, Codec) {
        let (method, bytes) = match codec {
            _ if data.len() < 16 => return (RAW, data.to_vec(), Codec::Raw),
            Codec::Raw => return (RAW, data.to_vec(), Codec::Raw),
            Codec::Rans0 => (RANS, rans::encode_o0(data)),
            Codec::Rans1 => (RANS, rans::encode_o1(data)),
            Codec::Gzip => (GZIP, gzip(&mut self.deflater, data)),
            Codec::Zlib5 => (GZIP, gzip_zlib(&mut self.zlib, data)),
        };
        if bytes.len() < data.len() {
            (method, bytes, codec)
        } else {
            (RAW, data.to_vec(), Codec::Raw)
        }
    }

    /// Encode one container. `slice_ref.bases` must be the contig for
    /// `input.ref_id` (or `None` when `ref_id == -1`).
    ///
    /// `plan`: the codec per external block from an earlier trial container;
    /// `None` tries every codec on every block (a trial container). The
    /// output depends only on the input, `params` and `plan`.
    pub fn encode(
        &mut self,
        input: &ContainerInput,
        slice_ref: &SliceRef,
        params: &EncodeParams,
        plan: Option<&CodecPlan>,
    ) -> io::Result<Encoded> {
        let n = input.len();
        let mut recs = Vec::with_capacity(n);
        for i in 0..n {
            let r = Rec::parse(input.record(i))?;
            if r.ref_id != input.ref_id {
                return Err(invalid("container mixes references"));
            }
            if r.mapped() {
                r.validate_mapped()?;
            }
            recs.push(r);
        }
        let ref_bases = slice_ref.bases;
        if input.ref_id >= 0 && ref_bases.is_none() && recs.iter().any(|r| r.mapped()) {
            return Err(invalid(format!(
                "reference sequence {} is not in the FASTA",
                input.ref_id
            )));
        }

        // Slice extent.
        let (start, end) = if input.ref_id >= 0 {
            let start = recs.iter().map(|r| r.apos).min().unwrap_or(0);
            let mut end = recs.iter().map(|r| r.aend.max(r.apos)).max().unwrap_or(0);
            // As htslib: a slice ends at the contig end even if an alignment
            // runs past it (those bases are stored verbatim).
            if let Some(b) = ref_bases {
                end = end.min(b.len() as i64).max(start);
            }
            (start, end)
        } else {
            (0, 0)
        };
        // Unplaced slices: start 0, span 1, as htslib's CRAM 3.0 writer
        // (its first/last base over the slice's records, apos = aend = 0).
        // Not span 0: htslib's index reader takes a slice's end as
        // start + span - 1, and `cram_index_query` skips slices with
        // end < pos, so with span 0 a `*` query (pos 0) jumps past every
        // unplaced slice but the last one.
        let span = if input.ref_id >= 0 {
            end - start + 1
        } else {
            1
        };
        let ap_delta = recs.windows(2).all(|w| w[0].apos <= w[1].apos);

        // Mate pairing: attach a primary record to the next primary record of
        // the same name when the decoders would rebuild both exactly.
        let mut cf = vec![CF_QS_ARRAY; n];
        let mut nf = vec![0i32; n];
        {
            let mut first: HashMap<&[u8], usize> = HashMap::new();
            let mut taken = vec![false; n];
            for j in 0..n {
                let r = &recs[j];
                if r.flag & PAIRED == 0 || !r.primary() {
                    continue;
                }
                match first.remove(r.name) {
                    Some(i) if !taken[i] && attachable(&recs[i], r) => {
                        cf[i] |= CF_MATE_DOWNSTREAM;
                        nf[i] = i32c(j - i - 1, "NF")?;
                        taken[i] = true;
                        taken[j] = true;
                    }
                    _ => {
                        first.insert(r.name, j);
                    }
                }
            }
            for j in 0..n {
                if taken[j] {
                    continue;
                }
                let r = &recs[j];
                // A record with no mate information at all needs no mate data.
                let bare = r.flag & (PAIRED | MATE_REVERSE | MATE_UNMAPPED) == 0
                    && r.next_ref == -1
                    && r.next_pos == -1
                    && r.tlen == 0;
                if !bare {
                    if r.flag & PAIRED == 0 {
                        return Err(invalid(format!(
                            "{}: unpaired record with mate fields cannot round-trip through CRAM",
                            String::from_utf8_lossy(r.name)
                        )));
                    }
                    cf[j] |= CF_DETACHED;
                }
            }
        }

        let mut s = Series::new();
        let mut last_apos = start;
        let mut feats: Vec<(usize, Feature)> = Vec::new();
        let mut taglist = Vec::new();
        let mut tags: Vec<([u8; 2], u8, &[u8])> = Vec::new();
        let mut num_bases: i64 = 0;
        for (i, r) in recs.iter().enumerate() {
            if r.flag >= 0x1000 {
                return Err(invalid(format!(
                    "{}: flag {:#x} has bits CRAM cannot store",
                    String::from_utf8_lossy(r.name),
                    r.flag
                )));
            }
            s.int(IntDs::Bf, r.flag as i32);
            s.int(IntDs::Cf, cf[i]);
            s.int(IntDs::Rl, i32c(r.l_seq, "read length")?);
            num_bases += r.l_seq as i64;
            if ap_delta {
                s.int(IntDs::Ap, i32c(r.apos - last_apos, "AP delta")?);
                last_apos = r.apos;
            } else {
                s.int(IntDs::Ap, i32c(r.apos, "alignment start")?);
            }
            s.arr(ArrDs::Rn, r.name);
            if cf[i] & CF_DETACHED != 0 {
                let mf = (r.flag & MATE_REVERSE != 0) as i32
                    | ((r.flag & MATE_UNMAPPED != 0) as i32) << 1;
                s.int(IntDs::Mf, mf);
                s.int(IntDs::Ns, r.next_ref);
                s.int(
                    IntDs::Np,
                    r.next_pos.checked_add(1).ok_or_else(|| {
                        invalid(format!(
                            "{}: mate position {} overflows",
                            String::from_utf8_lossy(r.name),
                            r.next_pos
                        ))
                    })?,
                );
                s.int(IntDs::Ts, r.tlen);
            } else if cf[i] & CF_MATE_DOWNSTREAM != 0 {
                s.int(IntDs::Nf, nf[i]);
            }

            // Tags: parsed, then stored as the tag mode says.
            taglist.clear();
            tags.clear();
            let aux = r.aux();
            let mut at = 0;
            let (mut has_md, mut has_nm) = (false, false);
            while at < aux.len() {
                if aux.len() - at < 3 {
                    return Err(invalid("truncated tag"));
                }
                let name = [aux[at], aux[at + 1]];
                let ty = aux[at + 2];
                let body = &aux[at + 3..];
                let len = crate::bam::aux_value_len(ty, body)
                    .filter(|&l| l <= body.len())
                    .ok_or_else(|| invalid("malformed tag"))?;
                has_md |= &name == b"MD";
                has_nm |= &name == b"NM";
                // A CG tag that held the CIGAR is not stored (htslib drops it).
                if r.cg_at != Some(at) {
                    tags.push((name, ty, &body[..len]));
                }
                at += 3 + len;
            }
            let mut rg = -1;
            if params.tags != TagMode::Lossless {
                // RG:Z naming a header read group moves to the RG series
                // (htslib cram_encode_aux; htsjdk always does this).
                if let Some(k) = tags.iter().position(|t| &t.0 == b"RG" && t.1 == b'Z') {
                    let id = tags[k].2.strip_suffix(&[0]).unwrap_or(tags[k].2);
                    if let Some(g) = params.read_groups.iter().position(|g| g[..] == *id) {
                        rg = g as i32;
                        tags.remove(k);
                    }
                }
                // MD/NM go when htslib would regenerate them identically.
                if r.mapped() && r.l_seq > 0 {
                    if let Some(b) = ref_bases {
                        let (md, nm, verbatim) = htslib_md_nm(r, b);
                        if !verbatim {
                            tags.retain(|t| {
                                let value = t.2;
                                if &t.0 == b"MD" && t.1 == b'Z' {
                                    let v = value.strip_suffix(&[0]).unwrap_or(value);
                                    !v.eq_ignore_ascii_case(&md)
                                } else if &t.0 == b"NM" {
                                    int_tag(t.1, value) != Some(nm)
                                } else {
                                    true
                                }
                            });
                        }
                    }
                }
                if matches!(params.tags, TagMode::HtsjdkSorted | TagMode::HtsjdkSortedMerged) {
                    // htsjdk 4 CompressionHeaderFactory sorts by the 3-byte key.
                    tags.sort_by_key(|t| (t.0, t.1));
                }
                if params.tags == TagMode::HtsjdkSortedMerged {
                    // samtools merge: RG, then PG, re-appended last.
                    let last = |t: &(_, u8, _)| match (t.0, t.1) {
                        (n, b'Z') if &n == b"RG" => 1,
                        (n, b'Z') if &n == b"PG" => 2,
                        _ => 0,
                    };
                    tags.sort_by_key(last);
                }
            }
            s.int(IntDs::Rg, rg);
            for &(name, ty, value) in &tags {
                let len = value.len();
                taglist.extend_from_slice(&[name[0], name[1], ty]);
                let t = s.tag(name, ty);
                match ty {
                    b'Z' | b'H' => {
                        t.has_tab |= value.contains(&b'\t');
                        t.data.extend_from_slice(value); // includes the NUL
                        t.data.push(b'\t');
                    }
                    b'B' => {
                        put_itf8(&mut t.data, i32c(len, "B tag length")?);
                        t.data.extend_from_slice(value);
                    }
                    _ => t.data.extend_from_slice(value),
                }
            }
            // Lossless: htslib generates MD/NM on decode for mapped reads that
            // lack them; cF:C (stripped by htslib's decoder) says not to.
            if params.tags == TagMode::Lossless && r.mapped() && (!has_md || !has_nm) {
                let cf_tag = (!has_md) as u8 | ((!has_nm) as u8) << 1;
                taglist.extend_from_slice(b"cFC");
                s.tag(*b"cF", b'C').data.push(cf_tag);
            }
            let tl = match s.td_index.get(&taglist) {
                Some(&t) => t,
                None => {
                    let t = i32c(s.td.len(), "tag lines")?;
                    s.td.push(taglist.clone());
                    s.td_index.insert(taglist.clone(), t);
                    t
                }
            };
            s.int(IntDs::Tl, tl);

            if r.mapped() {
                self.features(r, ref_bases.unwrap_or(&[]), &mut feats)?;
                s.int(IntDs::Fn, i32c(feats.len(), "feature count")?);
                let mut last = 0usize;
                for (pos, f) in feats.drain(..) {
                    let (code, ds_pos) = match &f {
                        Feature::Sub(..) => (b'X', pos),
                        Feature::Base(..) => (b'B', pos),
                        Feature::Ins(_) => (b'I', pos),
                        Feature::Del(_) => (b'D', pos),
                        Feature::Skip(_) => (b'N', pos),
                        Feature::Soft(_) => (b'S', pos),
                        Feature::Hard(_) => (b'H', pos),
                        Feature::Pad(_) => (b'P', pos),
                    };
                    s.byte(ByteDs::Fc, code);
                    s.int(IntDs::Fp, i32c(ds_pos - last, "feature position")?);
                    last = ds_pos;
                    match f {
                        Feature::Sub(rb, sb) => {
                            s.subs.push((rb, sb));
                            s.byte(ByteDs::Bs, 0); // patched once codes are known
                        }
                        Feature::Base(b, q) => {
                            s.byte(ByteDs::Ba, b);
                            s.byte(ByteDs::Qs, q);
                        }
                        Feature::Ins(b) => s.arr(ArrDs::In, &b),
                        Feature::Soft(b) => s.arr(ArrDs::Sc, &b),
                        Feature::Del(l) => s.int(IntDs::Dl, l),
                        Feature::Skip(l) => s.int(IntDs::Rs, l),
                        Feature::Hard(l) => s.int(IntDs::Hc, l),
                        Feature::Pad(l) => s.int(IntDs::Pd, l),
                    }
                }
                s.int(IntDs::Mq, r.mapq as i32);
            } else {
                if r.n_cigar != 0 || r.mapq != 0 {
                    return Err(invalid(format!(
                        "{}: unmapped record with CIGAR or MAPQ cannot round-trip through CRAM",
                        String::from_utf8_lossy(r.name)
                    )));
                }
                for k in 0..r.l_seq {
                    s.byte(ByteDs::Ba, r.base(k));
                }
            }
            s.bytes[ByteDs::Qs as usize].extend_from_slice(r.quals());
        }

        // Substitution matrix: per reference base, codes by descending count.
        let mut counts = [[0u32; 5]; 5];
        for &(rb, sb) in &s.subs {
            counts[rb as usize][sb as usize] += 1;
        }
        let mut code = [[0u8; 5]; 5];
        let mut sm = [0u8; 5];
        for rb in 0..5 {
            let mut others: Vec<usize> = (0..5).filter(|&b| b != rb).collect();
            others.sort_by_key(|&b| std::cmp::Reverse(counts[rb][b])); // stable: ties keep ACGTN order
            for (c, &b) in others.iter().enumerate() {
                code[rb][b] = c as u8;
            }
            let mut byte = 0u8;
            for b in (0..5).filter(|&b| b != rb) {
                byte = (byte << 2) | code[rb][b];
            }
            sm[rb] = byte;
        }
        {
            let bs = &mut s.bytes[ByteDs::Bs as usize];
            for (v, &(rb, sb)) in bs.iter_mut().zip(&s.subs) {
                *v = code[rb as usize][sb as usize];
            }
        }

        // Reference MD5 over the slice span.
        let mut md5 = [0u8; 16];
        if let Some(b) = ref_bases.filter(|_| input.ref_id >= 0 && span > 0) {
            // Bases past the contig end (if any) are `N` to decoders and not hashed.
            let lo = (start - 1).clamp(0, b.len() as i64) as usize;
            let hi = end.clamp(0, b.len() as i64) as usize;
            let mut m = Md5::new();
            m.update_upper(&b[lo..hi.max(lo)]);
            md5 = m.finish();
        }

        // Compression header.
        let mut ds_map = Vec::new();
        let mut ds_n = 0;
        let mut blocks: Vec<(i32, &[u8])> = Vec::new();
        for (k, (key, _)) in INT_DS.iter().enumerate() {
            let v = &s.ints[k];
            if v.is_empty() {
                continue;
            }
            ds_n += 1;
            ds_map.extend_from_slice(*key);
            if v.iter().all(|&x| x == v[0]) {
                put_huffman1(&mut ds_map, v[0]);
            } else {
                put_external(&mut ds_map, INT_ID0 + k as i32);
            }
        }
        let mut int_blocks: Vec<(i32, Vec<u8>)> = Vec::new();
        for (k, v) in s.ints.iter().enumerate() {
            if v.is_empty() || v.iter().all(|&x| x == v[0]) {
                continue;
            }
            let mut b = Vec::with_capacity(v.len() * 2);
            for &x in v {
                put_itf8(&mut b, x);
            }
            int_blocks.push((INT_ID0 + k as i32, b));
        }
        for (k, (key, ds)) in BYTE_DS.iter().enumerate() {
            // QS is always declared: every record has CF 0x1 (quality array),
            // and samtools looks its codec up even when no bases follow.
            if s.bytes[k].is_empty() && !matches!(ds, ByteDs::Qs) {
                continue;
            }
            ds_n += 1;
            ds_map.extend_from_slice(*key);
            put_external(&mut ds_map, BYTE_ID0 + k as i32);
            blocks.push((BYTE_ID0 + k as i32, &s.bytes[k]));
        }
        for (k, (key, _)) in ARR_DS.iter().enumerate() {
            if s.arrs[k].is_empty() {
                continue;
            }
            ds_n += 1;
            ds_map.extend_from_slice(*key);
            put_byte_array_stop(&mut ds_map, 0, ARR_ID0 + k as i32);
            blocks.push((ARR_ID0 + k as i32, &s.arrs[k]));
        }

        // Tag encodings. Z/H values that contain a tab are re-encoded with
        // an explicit length instead of the tab stop byte.
        let mut tag_blocks: Vec<(i32, Vec<u8>)> = Vec::new();
        let mut tag_map = Vec::new();
        for t in &s.tags {
            let mut key_bytes = Vec::new();
            put_itf8(&mut key_bytes, t.key);
            tag_map.extend_from_slice(&key_bytes);
            match t.ty {
                b'Z' | b'H' if !t.has_tab => {
                    put_byte_array_stop(&mut tag_map, b'\t', t.key);
                    tag_blocks.push((t.key, t.data.clone()));
                }
                b'Z' | b'H' => {
                    // Each value was appended as "bytes NUL TAB" and holds no
                    // NUL of its own, so split on NUL. Length and value share
                    // the block, interleaved in the order the decoder reads.
                    let d = &t.data;
                    let mut blk = Vec::with_capacity(d.len() + d.len() / 8);
                    let mut start = 0;
                    while start < d.len() {
                        let nul = start + d[start..].iter().position(|&c| c == 0).unwrap();
                        let v = &d[start..=nul];
                        put_itf8(&mut blk, i32c(v.len(), "Z tag length")?);
                        blk.extend_from_slice(v);
                        start = nul + 2;
                    }
                    put_byte_array_len_external(&mut tag_map, t.key, t.key);
                    tag_blocks.push((t.key, blk));
                }
                b'B' => {
                    put_byte_array_len_external(&mut tag_map, t.key, t.key);
                    tag_blocks.push((t.key, t.data.clone()));
                }
                ty => {
                    let width = match ty {
                        b'A' | b'c' | b'C' => 1,
                        b's' | b'S' => 2,
                        b'i' | b'I' | b'f' => 4,
                        _ => return Err(invalid(format!("unknown tag type {:?}", ty as char))),
                    };
                    put_byte_array_len_const(&mut tag_map, width, t.key);
                    tag_blocks.push((t.key, t.data.clone()));
                }
            }
        }

        let mut comp = Vec::new();
        // Preservation map.
        {
            let mut td = Vec::new();
            for l in &s.td {
                td.extend_from_slice(l);
                td.push(0);
            }
            let mut body = Vec::new();
            put_itf8(&mut body, 5);
            body.extend_from_slice(b"RN");
            body.push(1);
            body.extend_from_slice(b"AP");
            body.push(ap_delta as u8);
            body.extend_from_slice(b"RR");
            body.push(1);
            body.extend_from_slice(b"SM");
            body.extend_from_slice(&sm);
            body.extend_from_slice(b"TD");
            put_itf8(&mut body, i32c(td.len(), "TD size")?);
            body.extend_from_slice(&td);
            put_itf8(&mut comp, i32c(body.len(), "map size")?);
            comp.extend_from_slice(&body);
        }
        {
            let mut body = Vec::new();
            put_itf8(&mut body, ds_n);
            body.extend_from_slice(&ds_map);
            put_itf8(&mut comp, i32c(body.len(), "map size")?);
            comp.extend_from_slice(&body);
        }
        {
            let mut body = Vec::new();
            put_itf8(&mut body, s.tags.len() as i32);
            body.extend_from_slice(&tag_map);
            put_itf8(&mut comp, i32c(body.len(), "map size")?);
            comp.extend_from_slice(&body);
        }

        // Data blocks, compressed.
        let mut ext: Vec<(i32, u8, usize, Vec<u8>)> = Vec::new(); // id, method, raw len, bytes
        let mut codecs = Vec::new();
        for (id, data) in int_blocks
            .iter()
            .map(|(i, d)| (*i, d.as_slice()))
            .chain(blocks.iter().copied())
            .chain(tag_blocks.iter().map(|(i, d)| (*i, d.as_slice())))
        {
            let (method, bytes, codec) = match plan.and_then(|p| p.get(&id)) {
                Some(&codec) => self.compress_with(data, codec),
                None => self.compress(data),
            };
            codecs.push((id, codec, data.len()));
            ext.push((id, method, data.len(), bytes));
        }

        // Slice header.
        let mut sh = Vec::new();
        put_itf8(&mut sh, input.ref_id);
        put_itf8(&mut sh, i32c(start, "slice start")?);
        put_itf8(&mut sh, i32c(span, "slice span")?);
        put_itf8(&mut sh, i32c(n, "slice records")?);
        put_ltf8(&mut sh, input.record_counter as i64);
        put_itf8(&mut sh, i32c(1 + ext.len(), "block count")?);
        put_itf8(&mut sh, i32c(ext.len(), "block count")?);
        for e in &ext {
            put_itf8(&mut sh, e.0);
        }
        put_itf8(&mut sh, -1);
        sh.extend_from_slice(&md5);

        let mut body = Vec::new();
        put_block(&mut body, RAW, COMPRESSION_HEADER, 0, comp.len(), &comp)?;
        let landmark = body.len();
        put_block(&mut body, RAW, SLICE_HEADER, 0, sh.len(), &sh)?;
        put_block(&mut body, RAW, CORE_DATA, 0, 0, &[])?;
        for (id, method, raw_len, bytes) in &ext {
            put_block(&mut body, *method, EXTERNAL_DATA, *id, *raw_len, bytes)?;
        }
        let slice_size = body.len() - landmark;

        let mut bytes = Vec::with_capacity(body.len() + 64);
        put_container_header(
            &mut bytes,
            body.len(),
            input.ref_id,
            start,
            span,
            n,
            input.record_counter,
            num_bases,
            3 + ext.len(),
            &[landmark],
        )?;
        bytes.extend_from_slice(&body);
        Ok(Encoded {
            bytes,
            ref_id: input.ref_id,
            start,
            span,
            landmark,
            slice_size,
            records: n,
            codecs,
        })
    }

    /// Read features of a mapped record against `reference` (any case).
    fn features(
        &self,
        r: &Rec,
        reference: &[u8],
        out: &mut Vec<(usize, Feature)>,
    ) -> io::Result<()> {
        out.clear();
        let name = || String::from_utf8_lossy(r.name).into_owned();
        if r.l_seq == 0 {
            return Err(invalid(format!(
                "{}: mapped record without SEQ is not supported",
                name()
            )));
        }
        let quals = r.quals();
        let mut spos = 0usize;
        let mut rpos = (r.apos - 1) as usize;
        for c in r.cigar() {
            let (op, len) = (c & 0xf, (c >> 4) as usize);
            match op {
                0 => {
                    if spos + len > r.l_seq {
                        return Err(invalid(format!("{}: CIGAR longer than SEQ", name())));
                    }
                    for l in 0..len {
                        let sb = r.base(spos + l);
                        let rb = reference.get(rpos + l).map(|b| b.to_ascii_uppercase());
                        if rb == Some(sb) {
                            continue;
                        }
                        match (rb.and_then(acgtn), acgtn(sb)) {
                            (Some(ri), Some(si)) => {
                                out.push((spos + l + 1, Feature::Sub(ri as u8, si as u8)))
                            }
                            _ => out.push((spos + l + 1, Feature::Base(sb, quals[spos + l]))),
                        }
                    }
                    spos += len;
                    rpos += len;
                }
                1 => {
                    if spos + len > r.l_seq {
                        return Err(invalid(format!("{}: CIGAR longer than SEQ", name())));
                    }
                    out.push((
                        spos + 1,
                        Feature::Ins((spos..spos + len).map(|k| r.base(k)).collect()),
                    ));
                    spos += len;
                }
                2 => {
                    out.push((spos + 1, Feature::Del(i32c(len, "CIGAR op length")?)));
                    rpos += len;
                }
                3 => {
                    out.push((spos + 1, Feature::Skip(i32c(len, "CIGAR op length")?)));
                    rpos += len;
                }
                4 => {
                    if spos + len > r.l_seq {
                        return Err(invalid(format!("{}: CIGAR longer than SEQ", name())));
                    }
                    out.push((
                        spos + 1,
                        Feature::Soft((spos..spos + len).map(|k| r.base(k)).collect()),
                    ));
                    spos += len;
                }
                5 => out.push((spos + 1, Feature::Hard(i32c(len, "CIGAR op length")?))),
                6 => out.push((spos + 1, Feature::Pad(i32c(len, "CIGAR op length")?))),
                _ => {
                    return Err(invalid(format!(
                        "{}: CIGAR op {} is decoded as M by CRAM readers",
                        name(),
                        crate::bam::CIGAR_OPS
                            .get(op as usize)
                            .map(|&b| b as char)
                            .unwrap_or('?')
                    )));
                }
            }
        }
        if spos != r.l_seq {
            return Err(invalid(format!("{}: CIGAR and SEQ lengths differ", name())));
        }
        // The decoders rebuild the CIGAR from the features, merging runs; make
        // sure that gives back exactly this one.
        let rebuilt = rebuild_cigar(out, r.l_seq);
        if !rebuilt.iter().copied().eq(r.cigar()) {
            return Err(invalid(format!(
                "{}: CIGAR would not survive CRAM decoding (zero-length or adjacent same-type ops)",
                name()
            )));
        }
        Ok(())
    }
}

/// The CIGAR a CRAM decoder rebuilds from read features (htslib
/// `cram_decode_seq` without `USE_X`): matches and substitutions are `M`,
/// consecutive ops of one type merge.
fn rebuild_cigar(features: &[(usize, Feature)], l_seq: usize) -> Vec<u32> {
    let mut ops: Vec<u32> = Vec::new();
    let push = |ops: &mut Vec<u32>, op: u32, len: usize| {
        if len == 0 {
            return;
        }
        if let Some(last) = ops.last_mut().filter(|l| **l & 0xf == op) {
            *last += (len as u32) << 4;
            return;
        }
        ops.push((len as u32) << 4 | op);
    };
    let mut seq = 0usize; // bases consumed
    for (pos, f) in features {
        let p = pos - 1;
        if p > seq {
            push(&mut ops, 0, p - seq);
            seq = p;
        }
        match f {
            Feature::Sub(..) | Feature::Base(..) => {
                push(&mut ops, 0, 1);
                seq += 1;
            }
            Feature::Ins(b) => {
                push(&mut ops, 1, b.len());
                seq += b.len();
            }
            Feature::Soft(b) => {
                push(&mut ops, 4, b.len());
                seq += b.len();
            }
            Feature::Del(l) => push(&mut ops, 2, *l as usize),
            Feature::Skip(l) => push(&mut ops, 3, *l as usize),
            Feature::Hard(l) => push(&mut ops, 5, *l as usize),
            Feature::Pad(l) => push(&mut ops, 6, *l as usize),
        }
    }
    if seq < l_seq {
        push(&mut ops, 0, l_seq - seq);
    }
    ops
}

fn put_encoding(out: &mut Vec<u8>, codec: i32, params: &[u8]) {
    put_itf8(out, codec);
    put_itf8(out, params.len() as i32);
    out.extend_from_slice(params);
}

fn put_external(out: &mut Vec<u8>, id: i32) {
    let mut p = Vec::new();
    put_itf8(&mut p, id);
    put_encoding(out, E_EXTERNAL, &p);
}

/// HUFFMAN with a single symbol of code length 0: the value costs no bits.
fn put_huffman1(out: &mut Vec<u8>, value: i32) {
    let mut p = Vec::new();
    put_itf8(&mut p, 1);
    put_itf8(&mut p, value);
    put_itf8(&mut p, 1);
    put_itf8(&mut p, 0);
    put_encoding(out, E_HUFFMAN, &p);
}

fn put_byte_array_stop(out: &mut Vec<u8>, stop: u8, id: i32) {
    let mut p = vec![stop];
    put_itf8(&mut p, id);
    put_encoding(out, E_BYTE_ARRAY_STOP, &p);
}

/// BYTE_ARRAY_LEN with a constant length (zero-bit HUFFMAN) and EXTERNAL values.
fn put_byte_array_len_const(out: &mut Vec<u8>, len: i32, id: i32) {
    let mut p = Vec::new();
    put_huffman1(&mut p, len);
    put_external(&mut p, id);
    put_encoding(out, E_BYTE_ARRAY_LEN, &p);
}

/// BYTE_ARRAY_LEN with EXTERNAL lengths and values.
fn put_byte_array_len_external(out: &mut Vec<u8>, len_id: i32, id: i32) {
    let mut p = Vec::new();
    put_external(&mut p, len_id);
    put_external(&mut p, id);
    put_encoding(out, E_BYTE_ARRAY_LEN, &p);
}

/// Block: method, content type, content id, sizes, data, CRC32.
pub fn put_block(
    out: &mut Vec<u8>,
    method: u8,
    content_type: u8,
    id: i32,
    raw_len: usize,
    data: &[u8],
) -> io::Result<()> {
    let start = out.len();
    out.push(method);
    out.push(content_type);
    put_itf8(out, id);
    put_itf8(out, i32c(data.len(), "block size")?);
    put_itf8(out, i32c(raw_len, "block size")?);
    out.extend_from_slice(data);
    let crc = crc32(&out[start..]);
    out.extend_from_slice(&crc.to_le_bytes());
    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub fn put_container_header(
    out: &mut Vec<u8>,
    length: usize,
    ref_id: i32,
    start: i64,
    span: i64,
    records: usize,
    record_counter: u64,
    bases: i64,
    blocks: usize,
    landmarks: &[usize],
) -> io::Result<()> {
    let at = out.len();
    out.extend_from_slice(&i32c(length, "container size")?.to_le_bytes());
    put_itf8(out, ref_id);
    put_itf8(out, i32c(start, "container start")?);
    put_itf8(out, i32c(span, "container span")?);
    put_itf8(out, i32c(records, "container records")?);
    put_ltf8(out, record_counter as i64);
    put_ltf8(out, bases);
    put_itf8(out, i32c(blocks, "block count")?);
    put_itf8(out, i32c(landmarks.len(), "landmarks")?);
    for &l in landmarks {
        put_itf8(out, i32c(l, "landmark")?);
    }
    let crc = crc32(&out[at..]);
    out.extend_from_slice(&crc.to_le_bytes());
    Ok(())
}

/// A gzip member (RFC 1952) around our fast raw deflate.
pub fn gzip(deflater: &mut FastDeflater, data: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(data.len() / 3 + 32);
    out.extend_from_slice(&[0x1f, 0x8b, 8, 0, 0, 0, 0, 0, 0, 0xff]);
    deflater.compress(data, &mut out);
    out.extend_from_slice(&crc32(data).to_le_bytes());
    out.extend_from_slice(&(data.len() as u32).to_le_bytes());
    out
}

/// A gzip member around zlib-exact deflate (fresh state per call).
fn gzip_zlib(zlib: &mut crate::deflate::ZlibDeflater, data: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(data.len() / 3 + 32);
    out.extend_from_slice(&[0x1f, 0x8b, 8, 0, 0, 0, 0, 0, 0, 0xff]);
    zlib.compress(data, &mut out);
    out.extend_from_slice(&crc32(data).to_le_bytes());
    out.extend_from_slice(&(data.len() as u32).to_le_bytes());
    out
}

/// The SAM header container (block content type FILE_HEADER), gzip.
pub fn header_container(text: &[u8], deflater: &mut FastDeflater) -> io::Result<Vec<u8>> {
    let mut data = Vec::with_capacity(text.len() + 4);
    data.extend_from_slice(&i32c(text.len(), "header size")?.to_le_bytes());
    data.extend_from_slice(text);
    let gz = gzip(deflater, &data);
    let mut body = Vec::new();
    put_block(&mut body, GZIP, FILE_HEADER, 0, data.len(), &gz)?;
    let mut out = Vec::new();
    put_container_header(&mut out, body.len(), 0, 0, 0, 0, 0, 0, 1, &[])?;
    out.extend_from_slice(&body);
    Ok(out)
}

/// The 38-byte CRAM 3 end-of-file container.
pub const EOF_CONTAINER: [u8; 38] = [
    0x0f, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0x0f, 0xe0, 0x45, 0x4f, 0x46, 0x00, 0x00, 0x00,
    0x00, 0x01, 0x00, 0x05, 0xbd, 0xd9, 0x4f, 0x00, 0x01, 0x00, 0x06, 0x06, 0x01, 0x00, 0x01, 0x00,
    0x01, 0x00, 0xee, 0x63, 0x01, 0x4b,
];

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn eof_container_matches_our_writer() {
        let mut comp = Vec::new();
        for _ in 0..3 {
            put_itf8(&mut comp, 1);
            put_itf8(&mut comp, 0);
        }
        let mut body = Vec::new();
        put_block(&mut body, RAW, COMPRESSION_HEADER, 0, comp.len(), &comp).unwrap();
        let mut out = Vec::new();
        put_container_header(&mut out, body.len(), -1, 4542278, 0, 0, 0, 0, 1, &[]).unwrap();
        out.extend_from_slice(&body);
        assert_eq!(out, EOF_CONTAINER);
    }

    #[test]
    fn rebuild_merges_runs() {
        let f = vec![
            (3usize, Feature::Sub(0, 1)),
            (6, Feature::Del(2)),
            (6, Feature::Ins(vec![b'A'])),
        ];
        let c = rebuild_cigar(&f, 10);
        assert_eq!(c, vec![5 << 4, 2 << 4 | 2, 1 << 4 | 1, 4 << 4]);
    }
}
