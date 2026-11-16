//! A fast raw-DEFLATE compressor for BGZF blocks whose compressed bytes are
//! not part of any gate.
//!
//! md.bam, recal.bam and the per-chunk BAMs are gated on their records and
//! uncompressed stream (`docs/GATES.md`), so they need not carry htsjdk's
//! exact zlib level-2 bytes. This encoder trades a little size for speed:
//!
//! - **Stateless per call.** Every BGZF member is an independent DEFLATE
//!   stream; nothing carries over between blocks, so blocks can be compressed
//!   on any thread in any order.
//! - **Greedy single-probe matching.** A hash of the next 4 bytes indexes a
//!   table holding the last position with that hash. One candidate is checked
//!   (no chains, no lazy evaluation); matches are extended 8 bytes at a time.
//!   Positions inside a match are inserted too, which keeps the table fresh
//!   for the repetitive record layout of BAM.
//! - **Exact block-type choice.** Symbols are collected per DEFLATE block of
//!   at most `block_symbols`; the dynamic Huffman trees are built with zlib's
//!   length-limited `build_tree` (shared with [`super::ZlibFast`]), and the
//!   encoding among stored, static and dynamic is chosen by zlib's
//!   `_tr_flush_block` rule (stored when its length + 4 bytes is no more than
//!   the coded size), which can be a byte above the true minimum.
//! - **64-bit bit writer**, whole 32-bit words at a time.
//!
//! The output is valid DEFLATE for any input; correctness is checked against
//! two independent decoders (ours and the test oracle) and Python's zlib.

use super::{
    BL_CODES, BL_ORDER, BitSink, Ct, D_CODES, DYN_TREES, END_BLOCK, EXTRA_BLBITS, EXTRA_DBITS,
    EXTRA_LBITS, HEAP_SIZE, L_CODES, LITERALS, MAX_BITS, MAX_BL_BITS, MAX_MATCH, MIN_MATCH,
    STATIC_TREES, STORED_BLOCK, StaticDesc, scan_tree, send_tree, tables,
};
use crate::bgzf::Deflate;

/// DEFLATE's window.
const MAX_DISTANCE: usize = 32768;
/// Hash and minimum match length: 4 bytes suits BAM's binary fields.
const MIN_LEN: usize = 4;
/// A stored block holds at most 65,535 bytes; cut DEFLATE blocks before that.
const MAX_BLOCK_SPAN: usize = 65535 - MAX_MATCH;
/// log2 of the hash table size (128 KiB of `u32`).
const HASH_BITS: u32 = 15;
const HASH_SIZE: usize = 1 << HASH_BITS;
/// Symbol buffer; `block_symbols` stays below it.
const SYM_CAP: usize = 1 << 16;

/// Tuning knobs (benchmarked in `examples/deflate_bench.rs`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FastConfig {
    /// Matches up to this length have every covered position inserted into
    /// the hash table; longer ones only their last two.
    pub max_insert: usize,
    /// Symbols per DEFLATE block (Huffman trees are rebuilt per block). At
    /// most 65,000 so that the 16-bit frequency counts cannot overflow.
    pub block_symbols: usize,
}

impl Default for FastConfig {
    fn default() -> Self {
        FastConfig {
            max_insert: MAX_MATCH,
            block_symbols: 16383,
        }
    }
}

/// See the module docs.
pub struct FastDeflater {
    cfg: FastConfig,
    /// `base + position` of the last occurrence of each hash. Entries from
    /// earlier calls are at least `MAX_DISTANCE + 1` below the current base,
    /// so they fail the distance check and the table never needs clearing
    /// (except when `base` would wrap).
    table: Box<[u32; HASH_SIZE]>,
    base: u32,
    /// Literal (`< 256`) or `256 + (length - 3) | distance << 9`.
    syms: Box<[u32; SYM_CAP]>,
    nsyms: usize,
    ltree: Box<[Ct; HEAP_SIZE]>,
    dtree: [Ct; 2 * D_CODES + 1],
    bl_tree: [Ct; 2 * BL_CODES + 1],
}

impl Default for FastDeflater {
    fn default() -> Self {
        Self::new()
    }
}

impl FastDeflater {
    pub fn new() -> Self {
        Self::with_config(FastConfig::default())
    }

    pub fn with_config(cfg: FastConfig) -> Self {
        assert!(
            (1..=65000).contains(&cfg.block_symbols),
            "block_symbols must be 1..=65000"
        );
        FastDeflater {
            cfg,
            table: Box::new([0; HASH_SIZE]),
            base: MAX_DISTANCE as u32 + 1,
            syms: Box::new([0; SYM_CAP]),
            nsyms: 0,
            ltree: Box::new([Ct::default(); HEAP_SIZE]),
            dtree: [Ct::default(); 2 * D_CODES + 1],
            bl_tree: [Ct::default(); 2 * BL_CODES + 1],
        }
    }

    /// Blocks are flushed at `block_symbols` (< `SYM_CAP`), so the mask
    /// never wraps; it only spares the bounds check.
    #[inline]
    fn push(&mut self, sym: u32) {
        self.syms[self.nsyms & (SYM_CAP - 1)] = sym;
        self.nsyms += 1;
    }

    /// Always below `HASH_SIZE`, so indexing the fixed-size table needs no
    /// bounds check.
    #[inline]
    fn hash(v: u32) -> usize {
        (v.wrapping_mul(0x9E37_79B1) >> (32 - HASH_BITS)) as usize
    }

    fn clear_freqs(&mut self) {
        for ct in self.ltree[..L_CODES].iter_mut() {
            ct.fc = 0;
        }
        for ct in self.dtree[..D_CODES].iter_mut() {
            ct.fc = 0;
        }
        for ct in self.bl_tree[..BL_CODES].iter_mut() {
            ct.fc = 0;
        }
        self.ltree[END_BLOCK].fc = 1;
        self.nsyms = 0;
    }

    /// Compress `input` as one complete raw DEFLATE stream onto `out`.
    pub fn compress(&mut self, input: &[u8], out: &mut Vec<u8>) {
        let n = input.len();
        // Keep `base + position` (plus the gap to the next call) in range.
        let span = n as u64 + MAX_DISTANCE as u64 + 1;
        if self.base as u64 + span > u32::MAX as u64 {
            self.table.fill(0);
            self.base = MAX_DISTANCE as u32 + 1;
        }
        let base = self.base;
        // Saturate rather than wrap: a saturated base forces a reset next call.
        self.base = (base as u64 + span).min(u32::MAX as u64) as u32;

        let start = out.len();
        let mut bits = Bits64 {
            out: &mut *out,
            buf: 0,
            n: 0,
        };
        bits.out.reserve(n + n / 8 + 64);
        self.clear_freqs();
        let t = tables();
        let mut block_start = 0usize;
        let mut i = 0usize;
        while i + MIN_LEN <= n {
            let cur = load32(input, i);
            let h = Self::hash(cur);
            let cand = self.table[h];
            let pos = base + i as u32;
            self.table[h] = pos;
            let dist = pos.wrapping_sub(cand) as usize;
            // dist >= 1 and dist <= i hold by construction for BGZF-sized
            // calls (the base gap); check them anyway so no input can match
            // its own position or before the start (review 2026-10-23).
            let len = if (1..=MAX_DISTANCE.min(i)).contains(&dist) && load32(input, i - dist) == cur {
                match_len(input, i - dist, i, (n - i).min(MAX_MATCH))
            } else {
                0
            };
            if len >= MIN_LEN {
                self.push((dist as u32) << 9 | (256 + (len - MIN_MATCH)) as u32);
                self.ltree[t.length_code[len - MIN_MATCH] as usize + LITERALS + 1].fc += 1;
                self.dtree[t.d_code(dist - 1)].fc += 1;
                // Index the positions a match covers: all of them for short
                // matches, only the last two for long ones.
                let end = (i + len).min(n + 1 - MIN_LEN);
                let from = if len <= self.cfg.max_insert {
                    i + 1
                } else {
                    (i + len).saturating_sub(2).max(i + 1)
                };
                for p in from..end {
                    let h = Self::hash(load32(input, p));
                    self.table[h] = base + p as u32;
                }
                i += len;
            } else {
                self.push(input[i] as u32);
                self.ltree[input[i] as usize].fc += 1;
                i += 1;
            }
            if self.nsyms >= self.cfg.block_symbols || i - block_start >= MAX_BLOCK_SPAN {
                self.flush_block(input, block_start, i, false, &mut bits);
                block_start = i;
            }
        }
        while i < n {
            self.push(input[i] as u32);
            self.ltree[input[i] as usize].fc += 1;
            i += 1;
            if self.nsyms >= self.cfg.block_symbols || i - block_start >= MAX_BLOCK_SPAN {
                self.flush_block(input, block_start, i, false, &mut bits);
                block_start = i;
            }
        }
        self.flush_block(input, block_start, n, true, &mut bits);
        bits.windup();
        // Several stored sub-blocks cost 5 bytes each; one stored stream of
        // the whole input is never larger than that, and keeps a 65,498-byte
        // BGZF payload within one member.
        if out.len() - start > n + 5 && n <= 65535 {
            out.truncate(start);
            crate::bgzf::Stored.deflate(input, out);
        }
    }

    /// zlib `_tr_flush_block`'s choice among stored, static and dynamic
    /// encodings of `syms`, which cover `input[start..end]`.
    fn flush_block(
        &mut self,
        input: &[u8],
        start: usize,
        end: usize,
        last: bool,
        bits: &mut Bits64,
    ) {
        let t = tables();
        let l = StaticDesc {
            stree: Some(&t.static_ltree[..]),
            extra: &EXTRA_LBITS,
            extra_base: LITERALS + 1,
            elems: L_CODES,
            max_length: MAX_BITS,
        };
        let d = StaticDesc {
            stree: Some(&t.static_dtree[..]),
            extra: &EXTRA_DBITS,
            extra_base: 0,
            elems: D_CODES,
            max_length: MAX_BITS,
        };
        let bl = StaticDesc {
            stree: None,
            extra: &EXTRA_BLBITS,
            extra_base: 0,
            elems: BL_CODES,
            max_length: MAX_BL_BITS,
        };
        let (l_max, l_dyn, l_static) = build_code(&mut self.ltree[..], &l);
        let (d_max, d_dyn, d_static) = build_code(&mut self.dtree, &d);
        scan_tree(&mut self.ltree[..], l_max, &mut self.bl_tree);
        scan_tree(&mut self.dtree, d_max, &mut self.bl_tree);
        let (_, bl_dyn, _) = build_code(&mut self.bl_tree, &bl);
        let mut max_blindex = BL_CODES - 1;
        while max_blindex >= 3 && self.bl_tree[BL_ORDER[max_blindex]].dl == 0 {
            max_blindex -= 1;
        }
        let opt_len = l_dyn + d_dyn + bl_dyn + 3 * (max_blindex as u64 + 1) + 5 + 5 + 4;
        let opt_lenb = (opt_len + 3 + 7) >> 3;
        let static_lenb = (l_static + d_static + 3 + 7) >> 3;
        let stored_len = (end - start) as u64;
        let last_bit = last as u64;
        if stored_len + 4 <= opt_lenb.min(static_lenb) {
            bits.put((STORED_BLOCK as u64) << 1 | last_bit, 3);
            bits.align();
            let len = stored_len as u16;
            bits.out.extend_from_slice(&len.to_le_bytes());
            bits.out.extend_from_slice(&(!len).to_le_bytes());
            bits.out.extend_from_slice(&input[start..end]);
        } else if static_lenb <= opt_lenb {
            bits.put((STATIC_TREES as u64) << 1 | last_bit, 3);
            emit_symbols(
                &self.syms[..self.nsyms],
                &t.static_ltree,
                &t.static_dtree,
                bits,
            );
        } else {
            bits.put((DYN_TREES as u64) << 1 | last_bit, 3);
            let (lcodes, dcodes, blcodes) = (l_max + 1, d_max + 1, max_blindex as i32 + 1);
            bits.put((lcodes - 257) as u64, 5);
            bits.put((dcodes - 1) as u64, 5);
            bits.put((blcodes - 4) as u64, 4);
            for &o in &BL_ORDER[..blcodes as usize] {
                bits.put(self.bl_tree[o].dl as u64, 3);
            }
            send_tree(bits, &self.ltree[..], lcodes - 1, &self.bl_tree);
            send_tree(bits, &self.dtree, dcodes - 1, &self.bl_tree);
            emit_symbols(&self.syms[..self.nsyms], &self.ltree[..], &self.dtree, bits);
        }
        self.clear_freqs();
    }
}

impl Deflate for FastDeflater {
    fn deflate(&mut self, input: &[u8], out: &mut Vec<u8>) {
        self.compress(input, out);
    }
}

#[inline]
fn load32(b: &[u8], at: usize) -> u32 {
    u32::from_le_bytes(b[at..at + 4].try_into().unwrap())
}

#[inline]
fn load64(b: &[u8], at: usize) -> u64 {
    u64::from_le_bytes(b[at..at + 8].try_into().unwrap())
}

/// Length of the common prefix of `b[a..]` and `b[c..]` (`a < c`), at most
/// `max` (`c + max <= b.len()`). The first `MIN_LEN` bytes are known equal.
#[inline]
fn match_len(b: &[u8], a: usize, c: usize, max: usize) -> usize {
    let mut l = MIN_LEN;
    while l + 8 <= max {
        let x = load64(b, a + l) ^ load64(b, c + l);
        if x != 0 {
            return l + (x.trailing_zeros() / 8) as usize;
        }
        l += 8;
    }
    while l < max && b[a + l] == b[c + l] {
        l += 1;
    }
    l.min(max)
}

/// A length-limited canonical Huffman code for `tree[..desc.elems]`: `fc`
/// holds frequencies on entry and codes on return, `dl` the lengths.
/// Returns `(max_code, dynamic bits, static bits)`, extra bits included, for
/// the symbols counted (the static figure only when `desc.stree` is set).
///
/// Lengths come from Moffat and Katajainen's in-place minimum-redundancy
/// algorithm on the frequency-sorted symbols. Lengths above the limit are
/// clamped and the length counts repaired to a complete code (Kraft sum
/// exactly 1), then dealt out longest-first to the least frequent symbols. Fewer than two used symbols get a two-symbol code of
/// length-1 codes, as zlib forces too (single-code trees trip some decoders).
fn build_code(tree: &mut [Ct], desc: &StaticDesc) -> (i32, u64, u64) {
    let (elems, max_len) = (desc.elems, desc.max_length);
    // Sort keys: frequency above the symbol, so ties keep symbol order.
    let mut order = [0u32; L_CODES];
    let mut m = 0;
    let mut max_code: i32 = -1;
    for (n, ct) in tree[..elems].iter_mut().enumerate() {
        ct.dl = 0;
        if ct.fc != 0 {
            order[m] = (ct.fc as u32) << 9 | n as u32;
            m += 1;
            max_code = n as i32;
        }
    }
    let mut bl_count = [0u16; MAX_BITS + 1];
    if m < 2 {
        let present = (m == 1).then(|| (order[0] & 0x1ff) as usize);
        let other = if present == Some(0) { 1 } else { 0 };
        for s in [present.unwrap_or(1), other] {
            tree[s].dl = 1;
            max_code = max_code.max(s as i32);
        }
        bl_count[1] = 2;
    } else {
        let order = &mut order[..m];
        order.sort_unstable();
        let mut a = [0u32; L_CODES];
        for (w, &k) in a.iter_mut().zip(order.iter()) {
            *w = k >> 9;
        }
        minimum_redundancy(&mut a[..m]);
        for &len in &a[..m] {
            bl_count[(len as usize).min(max_len)] += 1;
        }
        // Clamping can oversubscribe the code. Kraft sum in units of
        // 2^-max_len: lengthen the longest codes below the limit until it
        // fits, then shorten the longest codes until it is complete again
        // (the deficit is always a multiple of the longest code's unit).
        let full = 1u64 << max_len;
        let mut kraft: u64 = (1..=max_len)
            .map(|b| (bl_count[b] as u64) << (max_len - b))
            .sum();
        while kraft > full {
            let b = (1..max_len).rev().find(|&b| bl_count[b] > 0).unwrap();
            bl_count[b] -= 1;
            bl_count[b + 1] += 1;
            kraft -= 1 << (max_len - b - 1);
        }
        while kraft < full {
            let l = (2..=max_len).rev().find(|&l| bl_count[l] > 0).unwrap();
            bl_count[l] -= 1;
            bl_count[l - 1] += 1;
            kraft += 1 << (max_len - l);
        }
        let mut k = 0;
        for bits in (1..=max_len).rev() {
            for _ in 0..bl_count[bits] {
                tree[(order[k] & 0x1ff) as usize].dl = bits as u16;
                k += 1;
            }
        }
    }
    let (mut dynamic, mut fixed) = (0u64, 0u64);
    for n in 0..=max_code as usize {
        let (f, len) = (tree[n].fc as u64, tree[n].dl as u64);
        let xbits = if n >= desc.extra_base {
            desc.extra[n - desc.extra_base] as u64
        } else {
            0
        };
        dynamic += f * (len + xbits);
        if let Some(stree) = desc.stree {
            fixed += f * (stree[n].dl as u64 + xbits);
        }
    }
    super::gen_codes(tree, max_code, &bl_count);
    (max_code, dynamic, fixed)
}

/// Moffat and Katajainen (1995), "In-place calculation of minimum-redundancy
/// codes": `a` holds at least two weights in non-decreasing order and is
/// replaced by their code lengths (non-increasing).
fn minimum_redundancy(a: &mut [u32]) {
    let n = a.len();
    // Phase 1: build the tree, leaving parent pointers in a[0..n-1].
    a[0] += a[1];
    let (mut root, mut leaf) = (0usize, 2usize);
    for next in 1..n - 1 {
        // First child.
        if leaf >= n || a[root] < a[leaf] {
            a[next] = a[root];
            a[root] = next as u32;
            root += 1;
        } else {
            a[next] = a[leaf];
            leaf += 1;
        }
        // Second child.
        if leaf >= n || (root < next && a[root] < a[leaf]) {
            a[next] += a[root];
            a[root] = next as u32;
            root += 1;
        } else {
            a[next] += a[leaf];
            leaf += 1;
        }
    }
    // Phase 2: internal node depths.
    a[n - 2] = 0;
    for next in (0..n.saturating_sub(2)).rev() {
        a[next] = a[a[next] as usize] + 1;
    }
    // Phase 3: leaf depths.
    let (mut avail, mut used, mut depth) = (1i64, 0i64, 0u32);
    let (mut root, mut next) = (n as isize - 2, n as isize - 1);
    while avail > 0 {
        while root >= 0 && a[root as usize] == depth {
            used += 1;
            root -= 1;
        }
        while avail > used {
            a[next as usize] = depth;
            next -= 1;
            avail -= 1;
        }
        avail = 2 * used;
        depth += 1;
        used = 0;
    }
}

/// Literal/length and distance codes with their extra bits.
fn emit_symbols(syms: &[u32], ltree: &[Ct], dtree: &[Ct], bits: &mut Bits64) {
    let t = tables();
    // Length symbol + extra bits per match length, merged into one put.
    let mut len_code = [(0u32, 0u32); MAX_MATCH - MIN_MATCH + 1];
    for (l, slot) in len_code.iter_mut().enumerate() {
        let code = t.length_code[l] as usize;
        let ct = ltree[code + LITERALS + 1];
        // A code this block never uses can hold scan_tree's 0xffff length
        // sentinel; it is never emitted, so leave its entry empty rather
        // than shift by it (which overflows; review 2026-10-23).
        if ct.dl as u32 >= 16 {
            *slot = (0, 0);
            continue;
        }
        let extra = EXTRA_LBITS[code] as u32;
        // Code 28 (length 258) has no extra bits (and zlib's base_length[28] is 0).
        let extra_value = if extra == 0 {
            0
        } else {
            l as u32 - t.base_length[code] as u32
        };
        let v = ct.fc as u32 | (extra_value << ct.dl);
        *slot = (v, ct.dl as u32 + extra);
    }
    for &s in syms {
        let lc = s & 0x1ff;
        if lc < 256 {
            let ct = ltree[lc as usize];
            bits.put(ct.fc as u64, ct.dl as u32);
        } else {
            let (v, n) = len_code[(lc - 256) as usize];
            bits.put(v as u64, n);
            let d = (s >> 9) as usize - 1;
            let code = t.d_code(d);
            let ct = dtree[code];
            let extra = EXTRA_DBITS[code] as u32;
            let v = ct.fc as u64 | (((d - t.base_dist[code] as usize) as u64) << ct.dl);
            bits.put(v, ct.dl as u32 + extra);
        }
    }
    let ct = ltree[END_BLOCK];
    bits.put(ct.fc as u64, ct.dl as u32);
}

/// LSB-first bit writer with a 64-bit accumulator; whole 32-bit words go out.
struct Bits64<'a> {
    out: &'a mut Vec<u8>,
    buf: u64,
    n: u32,
}

impl Bits64<'_> {
    /// Append the low `len` bits of `v` (`len <= 32`, `v < 2^len`).
    #[inline]
    fn put(&mut self, v: u64, len: u32) {
        self.buf |= v << self.n;
        self.n += len;
        if self.n >= 32 {
            self.out.extend_from_slice(&(self.buf as u32).to_le_bytes());
            self.buf >>= 32;
            self.n -= 32;
        }
    }

    /// Flush every pending bit, padding the last byte with zeros.
    fn align(&mut self) {
        while self.n > 0 {
            self.out.push(self.buf as u8);
            self.buf >>= 8;
            self.n = self.n.saturating_sub(8);
        }
        self.buf = 0;
    }

    fn windup(&mut self) {
        self.align();
    }
}

impl BitSink for Bits64<'_> {
    fn send_bits(&mut self, value: i32, length: i32) {
        self.put(value as u32 as u64 & ((1u64 << length) - 1), length as u32);
    }
}

#[cfg(test)]
mod tests {
    /// Review case 2026-10-23: no matches, dynamic wins, the highest symbol
    /// is EOB, so length code 257 holds scan_tree's sentinel. Must not panic
    /// (debug builds check shifts) and must round-trip.
    #[test]
    fn sentinel_length_code_round_trips() {
        let hex = "000000000100000002000000030000010100000102000001030000020100000202";
        let input: Vec<u8> = (0..hex.len())
            .step_by(2)
            .map(|k| u8::from_str_radix(&hex[k..k + 2], 16).unwrap())
            .collect();
        let mut out = Vec::new();
        FastDeflater::new().compress(&input, &mut out);
        let mut back = Vec::new();
        crate::inflate::inflate(&out, &mut back, input.len()).unwrap();
        assert_eq!(back, input);
    }

    use super::*;
    use crate::inflate::inflate;

    fn next(seed: &mut u64) -> u64 {
        *seed ^= *seed << 13;
        *seed ^= *seed >> 7;
        *seed ^= *seed << 17;
        *seed
    }

    fn roundtrip(d: &mut FastDeflater, data: &[u8]) -> usize {
        let mut out = vec![0xAA]; // compress appends
        d.compress(data, &mut out);
        let mut back = Vec::new();
        let used = inflate(&out[1..], &mut back, data.len()).unwrap();
        assert_eq!(used, out.len() - 1, "trailing bytes after the final block");
        assert!(back == data, "roundtrip mismatch, len {}", data.len());
        out.len() - 1
    }

    #[test]
    fn edge_inputs_roundtrip() {
        let mut d = FastDeflater::new();
        let mut seed = 0x1234_5678_9abc_def1u64;
        let random: Vec<u8> = (0..70_000).map(|_| next(&mut seed) as u8).collect();
        for len in [
            0, 1, 2, 3, 4, 5, 7, 8, 9, 15, 16, 17, 258, 259, 65497, 65498, 65499,
        ] {
            roundtrip(&mut d, &vec![b'A'; len]);
            roundtrip(&mut d, &random[..len]);
        }
        // Incompressible input must come out as stored blocks (~5 bytes overhead).
        let n = roundtrip(&mut d, &random[..65498]);
        assert_eq!(n, 65498 + 5, "random 65,498 bytes");
        // Long runs: 258-byte matches at distance 1.
        let n = roundtrip(&mut d, &vec![0u8; 65498]);
        assert!(n < 300, "zeros took {n}");
    }

    #[test]
    fn randomized_roundtrip_all_configs() {
        let mut seed = 0x9e37_79b9_7f4a_7c15u64;
        let configs = [
            FastConfig::default(),
            FastConfig {
                max_insert: 0,
                block_symbols: 1,
            },
            FastConfig {
                max_insert: 16,
                block_symbols: 65000,
            },
            FastConfig {
                max_insert: 0,
                block_symbols: 300,
            },
        ];
        for cfg in configs {
            let mut d = FastDeflater::with_config(cfg);
            for case in 0..120 {
                let n = match case {
                    0 => 0,
                    1 => 65536,
                    2 => 200_000,
                    _ => next(&mut seed) as usize % 70_000,
                };
                let mut data = Vec::with_capacity(n);
                let alphabet = 1 + next(&mut seed) as usize % 256;
                while data.len() < n {
                    let r = next(&mut seed);
                    if r & 3 == 0 && !data.is_empty() {
                        let dist = 1 + next(&mut seed) as usize % data.len().min(40_000);
                        for _ in 0..(r >> 8) as usize % 300 {
                            if data.len() == n {
                                break;
                            }
                            data.push(data[data.len() - dist]);
                        }
                    } else {
                        data.push(((r >> 20) as usize % alphabet) as u8);
                    }
                }
                roundtrip(&mut d, &data);
            }
        }
    }

    /// Kraft sum of the lengths in units of 2^-limit, and the longest length.
    fn kraft(tree: &[Ct], elems: usize, limit: usize) -> (u64, usize) {
        let lens = tree[..elems]
            .iter()
            .map(|c| c.dl as usize)
            .filter(|&l| l > 0);
        (
            lens.clone().map(|l| 1u64 << (limit - l)).sum(),
            lens.max().unwrap_or(0),
        )
    }

    #[test]
    fn huffman_codes_are_complete_and_length_limited() {
        let t = tables();
        let mut seed = 0x5eed_u64;
        for case in 0..2000 {
            let (elems, limit) = if case % 2 == 0 {
                (L_CODES, MAX_BITS)
            } else {
                (BL_CODES, MAX_BL_BITS)
            };
            let mut tree = vec![Ct::default(); 2 * L_CODES + 1];
            let used = 1 + next(&mut seed) as usize % elems;
            let (mut a, mut b) = (1u32, 1u32);
            for _ in 0..used {
                let s = next(&mut seed) as usize % elems;
                let f = match case % 3 {
                    // Fibonacci weights force lengths past the limit.
                    0 => {
                        (a, b) = (b, (a + b).min(60_000));
                        a
                    }
                    1 => 1 + next(&mut seed) as u32 % 3,
                    _ => 1 + next(&mut seed) as u32 % 60_000,
                };
                tree[s].fc = f.min(u16::MAX as u32) as u16;
            }
            let freqs: Vec<u16> = tree[..elems].iter().map(|c| c.fc).collect();
            let desc = StaticDesc {
                stree: (elems == L_CODES).then_some(&t.static_ltree[..]),
                extra: &EXTRA_LBITS,
                extra_base: if elems == L_CODES {
                    LITERALS + 1
                } else {
                    usize::MAX
                },
                elems,
                max_length: limit,
            };
            build_code(&mut tree, &desc);
            let (sum, longest) = kraft(&tree, elems, limit);
            assert_eq!(
                sum,
                1 << limit,
                "case {case}: incomplete or oversubscribed code"
            );
            assert!(longest <= limit, "case {case}: length {longest} > {limit}");
            for (n, &f) in freqs.iter().enumerate() {
                assert!(
                    f == 0 || tree[n].dl > 0,
                    "case {case}: used symbol {n} has no code"
                );
            }
        }
    }

    #[test]
    fn table_reuse_across_calls_and_base_wrap() {
        // Stale entries from earlier calls must never produce matches.
        let mut d = FastDeflater::new();
        let block: Vec<u8> = (0..65498u32).map(|i| (i % 251) as u8).collect();
        for _ in 0..4 {
            roundtrip(&mut d, &block);
        }
        d.base = u32::MAX - 70_000;
        for _ in 0..4 {
            roundtrip(&mut d, &block);
            roundtrip(&mut d, &block[..1000]);
        }
    }
}
