//! Per-read logic shared by BaseRecalibrator and ApplyBQSR: read decoding,
//! filters, clipping, reference comparison, known-site masks and covariates.
//! Section numbers (§) refer to `docs/spec/bqsr.md`.

use nomics_io::bam::RecordView;

pub const OP_M: u8 = 0;
pub const OP_I: u8 = 1;
pub const OP_D: u8 = 2;
pub const OP_N: u8 = 3;
pub const OP_S: u8 = 4;
pub const OP_H: u8 = 5;
pub const OP_EQ: u8 = 7;
pub const OP_X: u8 = 8;

pub fn consumes_read(op: u8) -> bool {
    matches!(op, OP_M | OP_I | OP_S | OP_EQ | OP_X)
}

pub fn consumes_ref(op: u8) -> bool {
    matches!(op, OP_M | OP_D | OP_N | OP_EQ | OP_X)
}

const SEQ_CODES: &[u8; 16] = b"=ACMGRSVTWYHKDBN";

/// `BaseUtils.simpleBaseToBaseIndex`: A/a/* -> 0, C/c -> 1, G/g -> 2, T/t -> 3, else -1.
pub fn base_index(b: u8) -> i32 {
    match b {
        b'A' | b'a' | b'*' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => -1,
    }
}

/// `BaseUtils.isRegularBase`.
pub fn is_regular_base(b: u8) -> bool {
    base_index(b) != -1
}

fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'a' => b't',
        b'c' => b'g',
        b'g' => b'c',
        b't' => b'a',
        x => x,
    }
}

/// The fields of a BAM record that BQSR needs, decoded once.
pub struct Read<'a> {
    pub flag: u16,
    pub ref_id: i32,
    /// 1-based POS; 0 when unplaced.
    pub pos1: i64,
    pub mapq: u8,
    pub cigar: Vec<(u32, u8)>,
    pub bases: Vec<u8>,
    pub quals: &'a [u8],
    pub tlen: i32,
    pub mate_ref: i32,
    /// 1-based mate POS; 0 when unplaced.
    pub mate_pos1: i64,
    pub rg: Option<&'a [u8]>,
    /// Byte offset of the QUAL field inside `raw`.
    pub qual_offset: usize,
}

impl<'a> Read<'a> {
    pub fn decode(raw: &'a [u8]) -> Option<Read<'a>> {
        let v = RecordView::new(raw).ok()?;
        let cigar: Vec<(u32, u8)> = v.cigar().map(|c| (c >> 4, (c & 0xf) as u8)).collect();
        let l_seq = u32::from_le_bytes(raw[16..20].try_into().ok()?) as usize;
        let seq_start = 32 + raw[8] as usize + 4 * cigar.len();
        let packed = raw.get(seq_start..seq_start + l_seq.div_ceil(2))?;
        let bases = (0..l_seq)
            .map(|i| SEQ_CODES[((packed[i / 2] >> (4 * (1 - i % 2))) & 0xf) as usize])
            .collect();
        let qual_offset = seq_start + l_seq.div_ceil(2);
        Some(Read {
            flag: v.flag(),
            ref_id: v.ref_id(),
            pos1: v.pos() as i64 + 1,
            mapq: raw[9],
            cigar,
            bases,
            quals: raw.get(qual_offset..qual_offset + l_seq)?,
            tlen: i32::from_le_bytes(raw[28..32].try_into().ok()?),
            mate_ref: v.mate_ref_id(),
            mate_pos1: i32::from_le_bytes(raw[24..28].try_into().ok()?) as i64 + 1,
            rg: v.string_tag(b"RG"),
            qual_offset,
        })
    }

    pub fn reverse(&self) -> bool {
        self.flag & 0x10 != 0
    }

    /// GATK `isUnmapped`: flag 0x4, no contig, or POS 0.
    pub fn is_unmapped(&self) -> bool {
        self.flag & 0x4 != 0 || self.ref_id < 0 || self.pos1 == 0
    }

    fn mate_is_unmapped(&self) -> bool {
        self.flag & 0x8 != 0 || self.mate_ref < 0 || self.mate_pos1 == 0
    }

    pub fn reference_length(&self) -> i64 {
        self.cigar
            .iter()
            .filter(|(_, op)| consumes_ref(*op))
            .map(|(l, _)| *l as i64)
            .sum()
    }

    /// 1-based inclusive alignment end (`POS + refLen - 1`).
    pub fn alignment_end(&self) -> i64 {
        self.pos1 + self.reference_length() - 1
    }

    fn cigar_read_length(&self) -> usize {
        self.cigar
            .iter()
            .filter(|(_, op)| consumes_read(*op))
            .map(|(l, _)| *l as usize)
            .sum()
    }

    /// `WellformedReadFilter` (§1.3).
    pub fn is_wellformed(&self, contig_lens: &[u32]) -> bool {
        let unmapped = self.is_unmapped();
        if !unmapped {
            if self.pos1 <= 0 {
                return false;
            }
            match contig_lens.get(self.ref_id as usize) {
                Some(&len) if self.pos1 <= len as i64 => {}
                _ => return false,
            }
        }
        self.rg.is_some()
            && self.bases.len() == self.quals.len()
            && (unmapped || self.bases.len() == self.cigar_read_length())
            && !self.bases.is_empty()
            && !self.cigar.iter().any(|(_, op)| *op == OP_N)
    }

    /// BaseRecalibrator's default read filters (§1.3).
    pub fn passes_recalibrator_filters(&self, contig_lens: &[u32]) -> bool {
        self.mapq != 0
            && self.mapq != 255
            && !self.is_unmapped()
            && self.flag & 0x100 == 0
            && self.flag & 0x400 == 0
            && self.flag & 0x200 == 0
            && self.is_wellformed(contig_lens)
    }

    /// The span used to query intervals (htsjdk overlap rule, §1.2):
    /// mate-placed unmapped reads use `[POS, POS]`.
    pub fn query_span(&self) -> Option<(i64, i64)> {
        if self.ref_id < 0 || self.pos1 == 0 {
            return None;
        }
        let end = if self.flag & 0x4 != 0 {
            self.pos1
        } else {
            self.alignment_end()
        };
        Some((self.pos1, end))
    }
}

/// `ReadUtils.getReadIndexForReferenceCoordinate(alignmentStart, cigar, refCoord)`:
/// the read index of the base at `coord` (for a coordinate inside a deletion,
/// the index of the base after it), treating `S` as consuming reference.
pub fn read_index_for_ref(
    alignment_start: i64,
    cigar: &[(u32, u8)],
    coord: i64,
) -> Option<(usize, u8)> {
    if coord < alignment_start {
        return None;
    }
    let (mut last_read, mut last_ref) = (0i64, alignment_start);
    for &(len, op) in cigar {
        let (first_read, first_ref) = (last_read, last_ref);
        last_read += if consumes_read(op) { len as i64 } else { 0 };
        last_ref += if consumes_ref(op) || op == OP_S {
            len as i64
        } else {
            0
        };
        if first_ref <= coord && coord < last_ref {
            let pos = first_read
                + if consumes_read(op) {
                    coord - first_ref
                } else {
                    0
                };
            return Some((pos as usize, op));
        }
    }
    None
}

/// A read after BaseRecalibrator's transforms (§2.1): bases `[lo, hi)` of the
/// original read, its alignment start and a CIGAR without clips.
pub struct Clipped {
    pub lo: usize,
    pub hi: usize,
    pub start: i64,
    pub cigar: Vec<(u32, u8)>,
}

impl Clipped {
    pub fn len(&self) -> usize {
        self.hi - self.lo
    }

    pub fn end(&self) -> i64 {
        self.start
            + self
                .cigar
                .iter()
                .filter(|(_, op)| consumes_ref(*op))
                .map(|(l, _)| *l as i64)
                .sum::<i64>()
            - 1
    }
}

/// Adaptor boundary (`ReadUtils.getAdaptorBoundary`), if well defined.
fn adaptor_boundary(r: &Read) -> Option<i64> {
    if r.tlen == 0 || r.flag & 0x1 == 0 || r.is_unmapped() || r.mate_is_unmapped() {
        return None;
    }
    let mate_reverse = r.flag & 0x20 != 0;
    if r.reverse() == mate_reverse {
        return None;
    }
    let ok = if r.reverse() {
        r.alignment_end() > r.mate_pos1
    } else {
        r.pos1 <= r.mate_pos1 + r.tlen as i64
    };
    if !ok {
        return None;
    }
    Some(if r.reverse() {
        r.mate_pos1 - 1
    } else {
        r.pos1 + (r.tlen as i64).abs()
    })
}

/// Consolidate the CIGAR like `CigarBuilder` for the cases bwa can produce:
/// merge equal neighbours and drop zero-length elements. Returns `None` for
/// shapes this port does not handle (leading/trailing deletions, `I D`).
fn consolidate(cigar: &[(u32, u8)]) -> Option<Vec<(u32, u8)>> {
    let mut out: Vec<(u32, u8)> = Vec::with_capacity(cigar.len());
    for &(len, op) in cigar {
        if len == 0 {
            continue;
        }
        match out.last_mut() {
            Some((l, o)) if *o == op => *l += len,
            _ => out.push((len, op)),
        }
    }
    let core: Vec<u8> = out
        .iter()
        .map(|(_, op)| *op)
        .filter(|op| !matches!(*op, OP_S | OP_H))
        .collect();
    if core.first() == Some(&OP_D)
        || core.last() == Some(&OP_D)
        || core.windows(2).any(|w| w == [OP_I, OP_D])
    {
        return None;
    }
    Some(out)
}

/// BaseRecalibrator's read transform (§2.1): consolidate CIGAR, hard-clip
/// adaptor sequence (soft clips count as aligned), then hard-clip soft clips.
/// `None` if the read ends up empty or has a CIGAR shape not handled here.
pub fn clip_for_recalibration(r: &Read) -> Option<Clipped> {
    let cigar = consolidate(&r.cigar)?;
    let len = r.bases.len();
    let (mut lo, mut hi) = (0usize, len);
    if let Some(boundary) = adaptor_boundary(r) {
        let (start, end) = (r.pos1, r.alignment_end());
        if start <= boundary && boundary <= end {
            let lead_s: i64 = cigar
                .iter()
                .take_while(|(_, op)| matches!(*op, OP_S | OP_H))
                .filter(|(_, op)| *op == OP_S)
                .map(|(l, _)| *l as i64)
                .sum();
            let soft_start = r.pos1 - lead_s;
            // Hard clips consume nothing; skip them when walking.
            let walk: Vec<(u32, u8)> = cigar
                .iter()
                .copied()
                .filter(|(_, op)| *op != OP_H)
                .collect();
            if let Some((idx, op)) = read_index_for_ref(soft_start, &walk, boundary) {
                if r.reverse() {
                    let stop = idx as i64 - if consumes_read(op) { 0 } else { 1 };
                    lo = (stop + 1).max(0) as usize;
                } else {
                    hi = idx;
                }
            }
        }
    }
    // Soft clips that survive are removed too.
    let lead_s: usize = cigar
        .iter()
        .take_while(|(_, op)| matches!(*op, OP_S | OP_H))
        .filter(|(_, op)| *op == OP_S)
        .map(|(l, _)| *l as usize)
        .sum();
    let trail_s: usize = cigar
        .iter()
        .rev()
        .take_while(|(_, op)| matches!(*op, OP_S | OP_H))
        .filter(|(_, op)| *op == OP_S)
        .map(|(l, _)| *l as usize)
        .sum();
    lo = lo.max(lead_s);
    hi = hi.min(len - trail_s);
    if lo >= hi {
        return None;
    }
    // Alignment start: POS shifted by the reference consumed before read index
    // `lo` (deletions right after the clip count, soft clips do not).
    let (mut start, mut read_i) = (r.pos1, 0usize);
    let mut new_cigar = Vec::new();
    for &(l, op) in &cigar {
        if op == OP_H {
            continue;
        }
        let l = l as usize;
        if consumes_read(op) {
            let (a, b) = (read_i, read_i + l);
            if b <= lo {
                start += if consumes_ref(op) { l as i64 } else { 0 };
            } else if a < lo {
                start += if consumes_ref(op) { (lo - a) as i64 } else { 0 };
            }
            let (ka, kb) = (a.max(lo), b.min(hi));
            if ka < kb && op != OP_S {
                new_cigar.push(((kb - ka) as u32, op));
            }
            read_i = b;
        } else if read_i <= lo {
            start += l as i64; // deletion before or immediately after the clip
        } else if read_i < hi {
            new_cigar.push((l as u32, op)); // internal deletion
        }
    }
    Some(Clipped {
        lo,
        hi,
        start,
        cigar: new_cigar,
    })
}

/// Per-base mismatch flags against the reference (§2.3): M/=/X bases only.
pub fn snp_errors(bases: &[u8], c: &Clipped, reference: &[u8]) -> Vec<u8> {
    let mut snp = vec![0u8; c.len()];
    let (mut read_i, mut ref_i) = (0usize, c.start - 1);
    for &(l, op) in &c.cigar {
        let l = l as usize;
        match op {
            OP_M | OP_EQ | OP_X => {
                for k in 0..l {
                    let rb = reference
                        .get((ref_i + k as i64) as usize)
                        .copied()
                        .unwrap_or(b'N');
                    if base_index(bases[read_i + k]) != base_index(rb) {
                        snp[read_i + k] = 1;
                    }
                }
                read_i += l;
                ref_i += l as i64;
            }
            OP_I => read_i += l,
            OP_D | OP_N => ref_i += l as i64,
            _ => {}
        }
    }
    snp
}

/// Known-site mask on the clipped read (§2.2, `calculateKnownSites`).
pub fn known_site_mask(c: &Clipped, sites: impl Iterator<Item = (i64, i64)>) -> Vec<bool> {
    let len = c.len();
    let mut mask = vec![false; len];
    let (soft_start, soft_end) = (c.start, c.end());
    for (s, e) in sites {
        if e < soft_start || s > soft_end {
            continue;
        }
        let mut fs: i64 = match read_index_for_ref(c.start, &c.cigar, s) {
            Some((i, op)) => i as i64 - if op == OP_D { 1 } else { 0 },
            None => 0,
        };
        let mut fe: i64 =
            read_index_for_ref(c.start, &c.cigar, e).map_or(len as i64, |(i, _)| i as i64);
        if fs > len as i64 {
            fs = len as i64;
            fe = len as i64;
        }
        let (a, b) = (fs.max(0) as usize, (fe + 1).min(len as i64).max(0) as usize);
        for m in mask.iter_mut().take(b).skip(a) {
            *m = true;
        }
    }
    mask
}

/// Context covariate keys for every offset of a read (§3.3): mismatch
/// context size 2 with low-quality-tail N masking. -1 = no context.
pub fn context_keys(bases: &[u8], quals: &[u8], reverse: bool, low_qual_tail: u8) -> Vec<i32> {
    const CONTEXT_SIZE: usize = 2;
    let n = bases.len();
    // clipLowQualEnds(WRITE_NS)
    let mut right = n as isize - 1;
    while right >= 0 && quals[right as usize] <= low_qual_tail {
        right -= 1;
    }
    let mut left = 0usize;
    while left < n && quals[left] <= low_qual_tail {
        left += 1;
    }
    if left as isize > right {
        return vec![0; n]; // emptied read: covariate values zeroed
    }
    let mut stranded: Vec<u8> = bases.to_vec();
    for b in &mut stranded[(right + 1) as usize..] {
        *b = b'N';
    }
    for b in &mut stranded[..left] {
        *b = b'N';
    }
    if reverse {
        stranded.reverse();
        for b in &mut stranded {
            *b = complement(*b);
        }
    }
    let keys = read_context_at_each_position(&stranded, CONTEXT_SIZE);
    let mut out = vec![0i32; n];
    for (i, &k) in keys.iter().enumerate() {
        let offset = if reverse { n - i - 1 } else { i };
        out[offset] = k;
    }
    out
}

pub(crate) fn key_from_context(dna: &[u8]) -> i32 {
    let mut key = dna.len() as i32;
    let mut bit_offset = 4;
    for &b in dna {
        let bi = base_index(b);
        if bi == -1 {
            return -1;
        }
        key |= bi << bit_offset;
        bit_offset += 2;
    }
    key
}

/// `ContextCovariate.getReadContextAtEachPosition`, ported literally.
fn read_context_at_each_position(bases: &[u8], context_size: usize) -> Vec<i32> {
    let read_length = bases.len();
    let mut keys = Vec::with_capacity(read_length);
    for _ in (1..context_size).take_while(|&i| i <= read_length) {
        keys.push(-1);
    }
    if read_length < context_size {
        return keys;
    }
    let mask = {
        let mut m = 0i32;
        for _ in 0..context_size {
            m = (m << 2) | 3;
        }
        m << 4
    };
    let new_base_offset = 2 * (context_size as i32 - 1) + 4;
    let mut current_key = key_from_context(&bases[..context_size]);
    keys.push(current_key);
    let mut current_n_penalty: i32 = 0;
    if current_key == -1 {
        current_key = 0;
        current_n_penalty = context_size as i32 - 1;
        let mut offset = new_base_offset;
        loop {
            let bi = base_index(bases[current_n_penalty as usize]);
            if bi == -1 {
                break;
            }
            current_key |= bi << offset;
            offset -= 2;
            current_n_penalty -= 1;
        }
    }
    for &b in &bases[context_size..] {
        let bi = base_index(b);
        if bi == -1 {
            current_n_penalty = context_size as i32;
            current_key = 0;
        } else {
            current_key = (current_key >> 2) & mask;
            current_key |= bi << new_base_offset;
            current_key |= context_size as i32;
        }
        if current_n_penalty == 0 {
            keys.push(current_key);
        } else {
            current_n_penalty -= 1;
            keys.push(-1);
        }
    }
    keys
}

/// Cycle covariate key for offset `i` of a read of length `len` (§3.4).
pub fn cycle_key(i: usize, len: usize, flag: u16) -> i32 {
    let f: i64 = if flag & 0x1 != 0 && flag & 0x80 != 0 {
        -1
    } else {
        1
    };
    let cycle = if flag & 0x10 != 0 {
        len as i64 * f + i as i64 * (-f)
    } else {
        f + i as i64 * f
    };
    assert!(
        cycle.abs() <= 500,
        "cycle {cycle} exceeds maximum_cycle_value 500"
    );
    ((cycle.abs() << 1) | (cycle < 0) as i64) as i32
}

/// Printed form of a cycle key.
pub fn cycle_from_key(key: usize) -> i64 {
    let c = (key >> 1) as i64;
    if key & 1 != 0 { -c } else { c }
}

/// Printed form of a size-2 context key (oldest base first).
pub fn context_from_key(key: usize) -> String {
    const B: &[u8; 4] = b"ACGT";
    let size = key & 15;
    (0..size)
        .map(|j| B[(key >> (4 + 2 * j)) & 3] as char)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn context_keys_match_spec_examples() {
        assert_eq!(key_from_context(b"AC"), 66);
        assert_eq!(key_from_context(b"TT"), 242);
        assert_eq!(context_from_key(66), "AC");
        // Forward read: first base has no context; N blocks its own and the next position.
        let k = read_context_at_each_position(b"ACNGTA", 2);
        assert_eq!(
            k,
            vec![
                -1,
                66,
                -1,
                -1,
                key_from_context(b"GT"),
                key_from_context(b"TA")
            ]
        );
    }

    #[test]
    fn reverse_strand_context_is_reverse_complemented() {
        // read AC (reverse): stranded = GT; offset 1 (first sequenced) has no context,
        // offset 0 gets context "GT".
        let keys = context_keys(b"AC", &[30, 30], true, 2);
        assert_eq!(keys, vec![key_from_context(b"GT"), -1]);
    }

    #[test]
    fn low_quality_tail_is_masked() {
        let keys = context_keys(b"ACGTA", &[30, 30, 30, 2, 2], false, 2);
        assert_eq!(keys[..3], [-1, 66, key_from_context(b"CG")]);
        assert_eq!(keys[3..], [-1, -1]);
    }

    #[test]
    fn cycles() {
        // first of pair forward: 1..len; second of pair forward: -1..-len; reverse counts down.
        assert_eq!(cycle_from_key(cycle_key(0, 151, 0x41) as usize), 1);
        assert_eq!(cycle_from_key(cycle_key(150, 151, 0x41) as usize), 151);
        assert_eq!(cycle_from_key(cycle_key(0, 151, 0x81) as usize), -1);
        assert_eq!(cycle_from_key(cycle_key(0, 151, 0x51) as usize), 151);
        assert_eq!(cycle_from_key(cycle_key(150, 151, 0x91) as usize), -1);
    }

    #[test]
    fn adaptor_and_soft_clip_examples_from_spec() {
        // Forward POS 1000 5S146M TLEN 140, mate reverse at 1100 -> boundary 1140:
        // bases 145..150 removed, then the 5 soft-clipped bases.
        let mut raw_cigar = vec![(5u32, OP_S), (146, OP_M)];
        let read = |cigar: Vec<(u32, u8)>, flag: u16, pos1: i64, mate: i64, tlen: i32| Read {
            flag,
            ref_id: 0,
            pos1,
            mapq: 60,
            cigar,
            bases: vec![b'A'; 151],
            quals: &[30; 151],
            tlen,
            mate_ref: 0,
            mate_pos1: mate,
            rg: Some(b"x"),
            qual_offset: 0,
        };
        let r = read(raw_cigar.clone(), 0x1 | 0x20, 1000, 1100, 140);
        let c = clip_for_recalibration(&r).unwrap();
        assert_eq!((c.lo, c.hi, c.start), (5, 145, 1000));
        assert_eq!(c.cigar, vec![(140, OP_M)]);
        // Reverse POS 1000 146M5S, mate forward at 1003 -> boundary 1002 -> first 3 bases clipped.
        raw_cigar = vec![(146, OP_M), (5, OP_S)];
        let r = read(raw_cigar, 0x1 | 0x10, 1000, 1003, -150);
        let c = clip_for_recalibration(&r).unwrap();
        assert_eq!((c.lo, c.hi, c.start), (3, 146, 1003));
    }

    #[test]
    fn known_site_in_deletion_masks_both_neighbours() {
        // 10M2D10M at 100: bases 0..9 on 100..109, deletion 110..111, bases 10..19 on 112..121.
        let c = Clipped {
            lo: 0,
            hi: 20,
            start: 100,
            cigar: vec![(10, OP_M), (2, OP_D), (10, OP_M)],
        };
        let mask = known_site_mask(&c, [(110, 110)].into_iter());
        let masked: Vec<usize> = (0..20).filter(|&i| mask[i]).collect();
        assert_eq!(masked, vec![9, 10]);
        // A deletion record REF=ACG at 104 covers 104..106 (padding base included).
        let mask = known_site_mask(&c, [(104, 106)].into_iter());
        assert_eq!(
            (0..20).filter(|&i| mask[i]).collect::<Vec<_>>(),
            vec![4, 5, 6]
        );
    }
}
