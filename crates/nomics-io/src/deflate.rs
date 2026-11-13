//! Raw DEFLATE compressor reproducing zlib 1.3.x `deflate_fast` (levels 1-3)
//! and `deflate_slow` (levels 4-9, lazy matching) byte for byte.
//!
//! Spec: zlib 1.3.1 `deflate.c` (`fill_window`, `longest_match`,
//! `deflate_fast`, `deflate_slow`) and `trees.c`, default build (no `FASTEST`, no
//! `UNALIGNED_OK`), `windowBits = -15`, `memLevel = 8`, default strategy.
//! The structure mirrors zlib deliberately, including:
//!
//! - **State that survives `deflateReset`.** htsjdk reuses one Deflater per
//!   file and resets it per BGZF block. Reset clears the hash heads and the
//!   tree state but not the window contents or `high_water`, and stale window
//!   bytes past the input can decide ties at the end of a block. One
//!   [`ZlibFast`] must therefore compress all blocks of a file, in order.
//! - **Output-capacity aborts.** htsjdk gives zlib a 65,518-byte output buffer
//!   ([`HTSJDK_OUTPUT_CAPACITY`]). If a block flush leaves no room, zlib stops
//!   mid-stream: the window holds only what was read so far, and htsjdk
//!   re-deflates the block at level 0.
//! - `ct_data` unions (`fc` = Freq/Code, `dl` = Dad/Len) and `ulg` wrapping
//!   arithmetic in the tree cost accounting.

use crate::bgzf::{Deflate, Stored};

pub mod fast;
pub use fast::FastDeflater;

const W_SIZE: usize = 1 << 15;
const W_MASK: usize = W_SIZE - 1;
const WINDOW_SIZE: usize = 2 * W_SIZE;
const HASH_BITS: u32 = 8 + 7; // memLevel + 7
const HASH_SIZE: usize = 1 << HASH_BITS;
const HASH_MASK: u32 = HASH_SIZE as u32 - 1;
const HASH_SHIFT: u32 = HASH_BITS.div_ceil(3);
const MIN_MATCH: usize = 3;
const MAX_MATCH: usize = 258;
const MIN_LOOKAHEAD: usize = MAX_MATCH + MIN_MATCH + 1;
const MAX_DIST: usize = W_SIZE - MIN_LOOKAHEAD;
const WIN_INIT: usize = MAX_MATCH;
/// `lit_bufsize - 1` symbols per block (memLevel 8).
const SYM_END: usize = (1 << (8 + 6)) - 1;

const LITERALS: usize = 256;
const LENGTH_CODES: usize = 29;
const L_CODES: usize = LITERALS + 1 + LENGTH_CODES;
const D_CODES: usize = 30;
const BL_CODES: usize = 19;
const HEAP_SIZE: usize = 2 * L_CODES + 1;
const MAX_BITS: usize = 15;
const MAX_BL_BITS: usize = 7;
const END_BLOCK: usize = 256;
const REP_3_6: usize = 16;
const REPZ_3_10: usize = 17;
const REPZ_11_138: usize = 18;
const STORED_BLOCK: i32 = 0;
const STATIC_TREES: i32 = 1;
const DYN_TREES: i32 = 2;

const EXTRA_LBITS: [u8; LENGTH_CODES] = [
    0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0,
];
const EXTRA_DBITS: [u8; D_CODES] = [
    0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13,
    13,
];
const EXTRA_BLBITS: [u8; BL_CODES] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 3, 7];
const BL_ORDER: [usize; BL_CODES] = [
    16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15,
];

/// htsjdk `BlockCompressedOutputStream` compressed buffer: 65536 - 18.
pub const HTSJDK_OUTPUT_CAPACITY: usize = 65518;

/// zlib `ct_data`: `fc` is Freq (while counting) or Code (after gen_codes);
/// `dl` is Dad (while building) or Len.
#[derive(Clone, Copy, Default)]
struct Ct {
    fc: u16,
    dl: u16,
}

struct StaticDesc {
    stree: Option<&'static [Ct]>,
    extra: &'static [u8],
    extra_base: usize,
    elems: usize,
    max_length: usize,
}

struct Tables {
    static_ltree: [Ct; L_CODES + 2],
    static_dtree: [Ct; D_CODES],
    length_code: [u8; MAX_MATCH - MIN_MATCH + 1],
    dist_code: [u8; 512],
    base_length: [i32; LENGTH_CODES],
    base_dist: [i32; D_CODES],
}

fn bi_reverse(mut code: u32, mut len: i32) -> u16 {
    let mut res = 0u32;
    loop {
        res |= code & 1;
        code >>= 1;
        res <<= 1;
        len -= 1;
        if len <= 0 {
            break;
        }
    }
    (res >> 1) as u16
}

fn gen_codes(tree: &mut [Ct], max_code: i32, bl_count: &[u16; MAX_BITS + 1]) {
    let mut next_code = [0u16; MAX_BITS + 1];
    let mut code: u32 = 0;
    for bits in 1..=MAX_BITS {
        code = (code + bl_count[bits - 1] as u32) << 1;
        next_code[bits] = code as u16;
    }
    for n in 0..=max_code.max(-1) {
        let n = n as usize;
        let len = tree[n].dl as i32;
        if len == 0 {
            continue;
        }
        tree[n].fc = bi_reverse(next_code[len as usize] as u32, len);
        next_code[len as usize] = next_code[len as usize].wrapping_add(1);
    }
}

impl Tables {
    /// zlib `tr_static_init`.
    fn new() -> Box<Tables> {
        let mut t = Box::new(Tables {
            static_ltree: [Ct::default(); L_CODES + 2],
            static_dtree: [Ct::default(); D_CODES],
            length_code: [0; MAX_MATCH - MIN_MATCH + 1],
            dist_code: [0; 512],
            base_length: [0; LENGTH_CODES],
            base_dist: [0; D_CODES],
        });
        let mut length = 0usize;
        let mut code = 0usize;
        while code < LENGTH_CODES - 1 {
            t.base_length[code] = length as i32;
            for _ in 0..(1 << EXTRA_LBITS[code]) {
                t.length_code[length] = code as u8;
                length += 1;
            }
            code += 1;
        }
        t.length_code[length - 1] = code as u8;
        let mut dist = 0usize;
        code = 0;
        while code < 16 {
            t.base_dist[code] = dist as i32;
            for _ in 0..(1 << EXTRA_DBITS[code]) {
                t.dist_code[dist] = code as u8;
                dist += 1;
            }
            code += 1;
        }
        dist >>= 7;
        while code < D_CODES {
            t.base_dist[code] = (dist << 7) as i32;
            for _ in 0..(1 << (EXTRA_DBITS[code] - 7)) {
                t.dist_code[256 + dist] = code as u8;
                dist += 1;
            }
            code += 1;
        }
        let mut bl_count = [0u16; MAX_BITS + 1];
        for (n, ct) in t.static_ltree.iter_mut().enumerate() {
            let len = match n {
                0..=143 => 8,
                144..=255 => 9,
                256..=279 => 7,
                _ => 8,
            };
            ct.dl = len;
            bl_count[len as usize] += 1;
        }
        gen_codes(&mut t.static_ltree, L_CODES as i32 + 1, &bl_count);
        for (n, ct) in t.static_dtree.iter_mut().enumerate() {
            ct.dl = 5;
            ct.fc = bi_reverse(n as u32, 5);
        }
        t
    }

    fn d_code(&self, dist: usize) -> usize {
        if dist < 256 {
            self.dist_code[dist] as usize
        } else {
            self.dist_code[256 + (dist >> 7)] as usize
        }
    }
}

/// Heap, depth and cost accumulators shared by the three trees (zlib keeps
/// these in `deflate_state`).
struct Builder {
    heap: [i32; HEAP_SIZE],
    heap_len: i32,
    heap_max: i32,
    depth: [u8; HEAP_SIZE],
    bl_count: [u16; MAX_BITS + 1],
    opt_len: u64,
    static_len: u64,
}

fn smaller(tree: &[Ct], n: i32, m: i32, depth: &[u8]) -> bool {
    let (n, m) = (n as usize, m as usize);
    tree[n].fc < tree[m].fc || (tree[n].fc == tree[m].fc && depth[n] <= depth[m])
}

impl Builder {
    fn pqdownheap(&mut self, tree: &[Ct], mut k: i32) {
        let v = self.heap[k as usize];
        let mut j = k << 1;
        while j <= self.heap_len {
            if j < self.heap_len
                && smaller(
                    tree,
                    self.heap[j as usize + 1],
                    self.heap[j as usize],
                    &self.depth,
                )
            {
                j += 1;
            }
            if smaller(tree, v, self.heap[j as usize], &self.depth) {
                break;
            }
            self.heap[k as usize] = self.heap[j as usize];
            k = j;
            j <<= 1;
        }
        self.heap[k as usize] = v;
    }

    fn gen_bitlen(&mut self, tree: &mut [Ct], max_code: i32, desc: &StaticDesc) {
        let max_length = desc.max_length as i32;
        let mut overflow = 0;
        self.bl_count = [0; MAX_BITS + 1];
        let root = self.heap[self.heap_max as usize] as usize;
        tree[root].dl = 0;
        let mut h = self.heap_max as usize + 1;
        while h < HEAP_SIZE {
            let n = self.heap[h] as usize;
            let mut bits = tree[tree[n].dl as usize].dl as i32 + 1;
            if bits > max_length {
                bits = max_length;
                overflow += 1;
            }
            tree[n].dl = bits as u16;
            h += 1;
            if n as i32 > max_code {
                continue;
            }
            self.bl_count[bits as usize] += 1;
            let xbits = if n >= desc.extra_base {
                desc.extra[n - desc.extra_base] as u64
            } else {
                0
            };
            let f = tree[n].fc as u64;
            self.opt_len = self.opt_len.wrapping_add(f * (bits as u64 + xbits));
            if let Some(stree) = desc.stree {
                self.static_len = self
                    .static_len
                    .wrapping_add(f * (stree[n].dl as u64 + xbits));
            }
        }
        if overflow == 0 {
            return;
        }
        loop {
            let mut bits = max_length - 1;
            while self.bl_count[bits as usize] == 0 {
                bits -= 1;
            }
            self.bl_count[bits as usize] -= 1;
            self.bl_count[bits as usize + 1] += 2;
            self.bl_count[max_length as usize] -= 1;
            overflow -= 2;
            if overflow <= 0 {
                break;
            }
        }
        let mut h = HEAP_SIZE;
        let mut bits = max_length;
        while bits != 0 {
            let mut n = self.bl_count[bits as usize];
            while n != 0 {
                h -= 1;
                let m = self.heap[h];
                if m > max_code {
                    continue;
                }
                let m = m as usize;
                if tree[m].dl as i32 != bits {
                    let delta = (bits as u64)
                        .wrapping_sub(tree[m].dl as u64)
                        .wrapping_mul(tree[m].fc as u64);
                    self.opt_len = self.opt_len.wrapping_add(delta);
                    tree[m].dl = bits as u16;
                }
                n -= 1;
            }
            bits -= 1;
        }
    }

    /// zlib `build_tree`; returns `max_code`.
    // Index loops are kept in zlib's shape on purpose: this is a byte-exact port.
    #[allow(clippy::needless_range_loop)]
    fn build_tree(&mut self, tree: &mut [Ct], desc: &StaticDesc) -> i32 {
        let elems = desc.elems;
        let mut max_code: i32 = -1;
        self.heap_len = 0;
        self.heap_max = HEAP_SIZE as i32;
        for n in 0..elems {
            if tree[n].fc != 0 {
                self.heap_len += 1;
                self.heap[self.heap_len as usize] = n as i32;
                max_code = n as i32;
                self.depth[n] = 0;
            } else {
                tree[n].dl = 0;
            }
        }
        while self.heap_len < 2 {
            let node = if max_code < 2 {
                max_code += 1;
                max_code
            } else {
                0
            };
            self.heap_len += 1;
            self.heap[self.heap_len as usize] = node;
            let node = node as usize;
            tree[node].fc = 1;
            self.depth[node] = 0;
            self.opt_len = self.opt_len.wrapping_sub(1);
            if let Some(stree) = desc.stree {
                self.static_len = self.static_len.wrapping_sub(stree[node].dl as u64);
            }
        }
        let mut n = self.heap_len / 2;
        while n >= 1 {
            self.pqdownheap(tree, n);
            n -= 1;
        }
        let mut node = elems;
        loop {
            // pqremove
            let n = self.heap[1];
            self.heap[1] = self.heap[self.heap_len as usize];
            self.heap_len -= 1;
            self.pqdownheap(tree, 1);
            let m = self.heap[1];
            self.heap_max -= 1;
            self.heap[self.heap_max as usize] = n;
            self.heap_max -= 1;
            self.heap[self.heap_max as usize] = m;
            let (nu, mu) = (n as usize, m as usize);
            tree[node].fc = tree[nu].fc.wrapping_add(tree[mu].fc);
            self.depth[node] = self.depth[nu].max(self.depth[mu]) + 1;
            tree[nu].dl = node as u16;
            tree[mu].dl = node as u16;
            self.heap[1] = node as i32;
            node += 1;
            self.pqdownheap(tree, 1);
            if self.heap_len < 2 {
                break;
            }
        }
        self.heap_max -= 1;
        self.heap[self.heap_max as usize] = self.heap[1];
        self.gen_bitlen(tree, max_code, desc);
        gen_codes(tree, max_code, &self.bl_count);
        max_code
    }
}

/// zlib `scan_tree`: count bit-length code frequencies into `bl_tree`.
fn scan_tree(tree: &mut [Ct], max_code: i32, bl_tree: &mut [Ct]) {
    let mut prevlen: i32 = -1;
    let mut nextlen = tree[0].dl as i32;
    let mut count = 0;
    let (mut max_count, mut min_count) = if nextlen == 0 { (138, 3) } else { (7, 4) };
    tree[(max_code + 1) as usize].dl = 0xffff; // guard
    for n in 0..=max_code {
        let curlen = nextlen;
        nextlen = tree[(n + 1) as usize].dl as i32;
        count += 1;
        if count < max_count && curlen == nextlen {
            continue;
        } else if count < min_count {
            bl_tree[curlen as usize].fc = bl_tree[curlen as usize].fc.wrapping_add(count as u16);
        } else if curlen != 0 {
            if curlen != prevlen {
                bl_tree[curlen as usize].fc = bl_tree[curlen as usize].fc.wrapping_add(1);
            }
            bl_tree[REP_3_6].fc = bl_tree[REP_3_6].fc.wrapping_add(1);
        } else if count <= 10 {
            bl_tree[REPZ_3_10].fc = bl_tree[REPZ_3_10].fc.wrapping_add(1);
        } else {
            bl_tree[REPZ_11_138].fc = bl_tree[REPZ_11_138].fc.wrapping_add(1);
        }
        count = 0;
        prevlen = curlen;
        (max_count, min_count) = if nextlen == 0 {
            (138, 3)
        } else if curlen == nextlen {
            (6, 3)
        } else {
            (7, 4)
        };
    }
}

/// zlib's bit writer (`bi_buf`/`bi_valid`, 16-bit buffer) over a byte vector.
struct Bits<'a> {
    out: &'a mut Vec<u8>,
    buf: u16,
    valid: i32,
}

/// Where `send_tree` writes its bits: zlib's 16-bit writer here, a 64-bit
/// one in [`fast`].
trait BitSink {
    fn send_bits(&mut self, value: i32, length: i32);
    fn send_code(&mut self, c: usize, tree: &[Ct]) {
        self.send_bits(tree[c].fc as i32, tree[c].dl as i32);
    }
}

impl BitSink for Bits<'_> {
    fn send_bits(&mut self, value: i32, length: i32) {
        Bits::send_bits(self, value, length);
    }
}

impl Bits<'_> {
    fn put_short(&mut self, w: u16) {
        self.out.extend_from_slice(&w.to_le_bytes());
    }

    fn send_bits(&mut self, value: i32, length: i32) {
        let v = value as u16 as u32;
        if self.valid > 16 - length {
            self.buf |= (v << self.valid) as u16;
            let b = self.buf;
            self.put_short(b);
            self.buf = (v >> (16 - self.valid)) as u16;
            self.valid += length - 16;
        } else {
            self.buf |= (v << self.valid) as u16;
            self.valid += length;
        }
    }

    fn send_code(&mut self, c: usize, tree: &[Ct]) {
        self.send_bits(tree[c].fc as i32, tree[c].dl as i32);
    }

    fn windup(&mut self) {
        if self.valid > 8 {
            let b = self.buf;
            self.put_short(b);
        } else if self.valid > 0 {
            self.out.push(self.buf as u8);
        }
        self.buf = 0;
        self.valid = 0;
    }
}

/// zlib `configuration_table` rows for levels 1-9.
#[derive(Clone, Copy)]
struct Config {
    /// `good_length`: quarter the chain once a match this long is in hand
    /// (only reachable in `deflate_slow`).
    good_length: usize,
    /// `max_lazy`: in `deflate_fast` the longest match whose positions are
    /// all inserted (`max_insert_length`); in `deflate_slow` the longest
    /// previous match for which a better one is still sought.
    max_insert_length: usize,
    nice_length: usize,
    max_chain: usize,
}

const fn cfg(good_length: usize, max_lazy: usize, nice_length: usize, max_chain: usize) -> Config {
    Config {
        good_length,
        max_insert_length: max_lazy,
        nice_length,
        max_chain,
    }
}

const CONFIGS: [Config; 9] = [
    cfg(4, 4, 8, 4),
    cfg(4, 5, 16, 8),
    cfg(4, 6, 32, 32),
    cfg(4, 4, 16, 16),
    cfg(8, 16, 32, 32),
    cfg(8, 16, 128, 128),
    cfg(8, 32, 128, 256),
    cfg(32, 128, 258, 1024),
    cfg(32, 258, 258, 4096),
];

/// zlib `TOO_FAR`: a length-3 match farther than this is not worth coding
/// (`deflate_slow` only).
const TOO_FAR: usize = 4096;

/// Outcome of one `deflate(..., Z_FINISH)` call into a bounded buffer.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Finish {
    /// The stream is complete.
    Done,
    /// The output reached capacity at a block flush; zlib stopped there.
    OutOfSpace,
}

/// Stateful zlib-exact compressor: `deflate_fast` for levels 1-3,
/// `deflate_slow` for 4-9. See the module docs.
pub struct ZlibFast {
    cfg: Config,
    lazy: bool,
    // deflate_slow state.
    prev_length: usize,
    prev_match: usize,
    match_available: bool,
    window: Box<[u8]>,
    high_water: usize,
    prev: Box<[u16]>,
    head: Box<[u16]>,
    ins_h: u32,
    strstart: usize,
    block_start: isize,
    lookahead: usize,
    insert: usize,
    match_start: usize,
    match_length: usize,
    // Current input.
    next_in: usize,
    avail_in: usize,
    // Trees for the current block.
    dyn_ltree: Box<[Ct; HEAP_SIZE]>,
    dyn_dtree: [Ct; 2 * D_CODES + 1],
    bl_tree: [Ct; 2 * BL_CODES + 1],
    builder: Box<Builder>,
    l_max_code: i32,
    d_max_code: i32,
    syms: Vec<(u16, u8)>,
    // Bit writer state persists across flushes within one stream.
    bi_buf: u16,
    bi_valid: i32,
}

impl ZlibFast {
    /// A fresh compressor (zlib `deflateInit2`) for `level` 1-9.
    pub fn new(level: u32) -> ZlibFast {
        assert!((1..=9).contains(&level), "zlib levels are 1-9");
        ZlibFast {
            cfg: CONFIGS[level as usize - 1],
            lazy: level >= 4,
            prev_length: MIN_MATCH - 1,
            prev_match: 0,
            match_available: false,
            // zlib mallocs the window; bytes are zeroed via high_water before any
            // read can observe them, so starting from zeros is equivalent.
            window: vec![0u8; WINDOW_SIZE].into_boxed_slice(),
            high_water: 0,
            prev: vec![0u16; W_SIZE].into_boxed_slice(),
            head: vec![0u16; HASH_SIZE].into_boxed_slice(),
            ins_h: 0,
            strstart: 0,
            block_start: 0,
            lookahead: 0,
            insert: 0,
            match_start: 0,
            match_length: MIN_MATCH - 1,
            next_in: 0,
            avail_in: 0,
            dyn_ltree: Box::new([Ct::default(); HEAP_SIZE]),
            dyn_dtree: [Ct::default(); 2 * D_CODES + 1],
            bl_tree: [Ct::default(); 2 * BL_CODES + 1],
            builder: Box::new(Builder {
                heap: [0; HEAP_SIZE],
                heap_len: 0,
                heap_max: 0,
                depth: [0; HEAP_SIZE],
                bl_count: [0; MAX_BITS + 1],
                opt_len: 0,
                static_len: 0,
            }),
            l_max_code: 0,
            d_max_code: 0,
            syms: Vec::with_capacity(SYM_END),
            bi_buf: 0,
            bi_valid: 0,
        }
    }

    /// Back to `deflateInit2` state: zeroed window and `high_water` 0.
    /// See [`crate::bgzf::ParallelBgzfWriter`] for when this reproduces a
    /// reused deflater exactly.
    pub fn reinit(&mut self) {
        self.window.fill(0);
        self.high_water = 0;
    }

    /// zlib `deflateReset`: lm_init + _tr_init. Window and high_water persist.
    fn reset(&mut self) {
        self.head.fill(0);
        self.strstart = 0;
        self.block_start = 0;
        self.lookahead = 0;
        self.insert = 0;
        self.match_length = MIN_MATCH - 1;
        self.prev_length = MIN_MATCH - 1;
        self.match_available = false;
        self.ins_h = 0;
        self.bi_buf = 0;
        self.bi_valid = 0;
        self.init_block();
    }

    fn init_block(&mut self) {
        for ct in self.dyn_ltree[..L_CODES].iter_mut() {
            ct.fc = 0;
        }
        for ct in self.dyn_dtree[..D_CODES].iter_mut() {
            ct.fc = 0;
        }
        for ct in self.bl_tree[..BL_CODES].iter_mut() {
            ct.fc = 0;
        }
        self.dyn_ltree[END_BLOCK].fc = 1;
        self.builder.opt_len = 0;
        self.builder.static_len = 0;
        self.syms.clear();
    }

    #[inline]
    fn update_hash(&mut self, c: u8) {
        self.ins_h = ((self.ins_h << HASH_SHIFT) ^ c as u32) & HASH_MASK;
    }

    /// INSERT_STRING; returns the previous head of the chain.
    #[inline]
    fn insert_string(&mut self, s: usize) -> usize {
        self.update_hash(self.window[s + MIN_MATCH - 1]);
        let h = self.ins_h as usize;
        let head = self.head[h];
        self.prev[s & W_MASK] = head;
        self.head[h] = s as u16;
        head as usize
    }

    fn slide_hash(&mut self) {
        for p in self.head.iter_mut().chain(self.prev.iter_mut()) {
            let m = *p as usize;
            *p = if m >= W_SIZE { (m - W_SIZE) as u16 } else { 0 };
        }
    }

    fn fill_window(&mut self, input: &[u8]) {
        loop {
            let mut more = WINDOW_SIZE - self.lookahead - self.strstart;
            if self.strstart >= W_SIZE + MAX_DIST {
                self.window.copy_within(W_SIZE..W_SIZE + (W_SIZE - more), 0);
                self.match_start = self.match_start.wrapping_sub(W_SIZE);
                self.strstart -= W_SIZE;
                self.block_start -= W_SIZE as isize;
                if self.insert > self.strstart {
                    self.insert = self.strstart;
                }
                self.slide_hash();
                more += W_SIZE;
            }
            if self.avail_in == 0 {
                break;
            }
            let n = self.avail_in.min(more);
            let dst = self.strstart + self.lookahead;
            self.window[dst..dst + n].copy_from_slice(&input[self.next_in..self.next_in + n]);
            self.next_in += n;
            self.avail_in -= n;
            self.lookahead += n;
            if self.lookahead + self.insert >= MIN_MATCH {
                let mut s = self.strstart - self.insert;
                self.ins_h = self.window[s] as u32;
                self.update_hash(self.window[s + 1]);
                while self.insert != 0 {
                    self.update_hash(self.window[s + MIN_MATCH - 1]);
                    let h = self.ins_h as usize;
                    self.prev[s & W_MASK] = self.head[h];
                    self.head[h] = s as u16;
                    s += 1;
                    self.insert -= 1;
                    if self.lookahead + self.insert < MIN_MATCH {
                        break;
                    }
                }
            }
            if !(self.lookahead < MIN_LOOKAHEAD && self.avail_in != 0) {
                break;
            }
        }
        if self.high_water < WINDOW_SIZE {
            let curr = self.strstart + self.lookahead;
            if self.high_water < curr {
                let init = (WINDOW_SIZE - curr).min(WIN_INIT);
                self.window[curr..curr + init].fill(0);
                self.high_water = curr + init;
            } else if self.high_water < curr + WIN_INIT {
                let init = (curr + WIN_INIT - self.high_water).min(WINDOW_SIZE - self.high_water);
                self.window[self.high_water..self.high_water + init].fill(0);
                self.high_water += init;
            }
        }
    }

    /// `prev_length` is `MIN_MATCH - 1` in `deflate_fast` (so `good_match`
    /// never applies there) and the previous match length in `deflate_slow`.
    fn longest_match(&mut self, mut cur_match: usize, prev_length: usize) -> usize {
        let w = &self.window;
        let mut chain_length = self.cfg.max_chain;
        if prev_length >= self.cfg.good_length {
            chain_length >>= 2;
        }
        let scan = self.strstart;
        let mut best_len = prev_length;
        let nice_match = self.cfg.nice_length.min(self.lookahead);
        let limit = self.strstart.saturating_sub(MAX_DIST);
        let strend = scan + MAX_MATCH;
        let mut scan_end1 = w[scan + best_len - 1];
        let mut scan_end = w[scan + best_len];
        loop {
            let m = cur_match;
            // zlib skips comparing scan[2]/match[2]: equal hashes imply equal bytes.
            if w[m + best_len] == scan_end
                && w[m + best_len - 1] == scan_end1
                && w[m] == w[scan]
                && w[m + 1] == w[scan + 1]
            {
                let (mut s, mut mm) = (scan + 2, m + 2);
                'scan: loop {
                    for _ in 0..8 {
                        s += 1;
                        mm += 1;
                        if w[s] != w[mm] {
                            break 'scan;
                        }
                    }
                    if s >= strend {
                        break;
                    }
                }
                let len = MAX_MATCH - (strend - s);
                if len > best_len {
                    self.match_start = cur_match;
                    best_len = len;
                    if len >= nice_match {
                        break;
                    }
                    scan_end1 = w[scan + best_len - 1];
                    scan_end = w[scan + best_len];
                }
            }
            cur_match = self.prev[cur_match & W_MASK] as usize;
            if cur_match <= limit {
                break;
            }
            chain_length -= 1;
            if chain_length == 0 {
                break;
            }
        }
        best_len.min(self.lookahead)
    }

    /// Record a literal (dist 0) or match; true when the symbol buffer is full.
    #[inline]
    fn tally(&mut self, dist: usize, lc: usize) -> bool {
        self.syms.push((dist as u16, lc as u8));
        if dist == 0 {
            self.dyn_ltree[lc].fc += 1;
        } else {
            let d = dist - 1;
            let lcode = tables().length_code[lc] as usize;
            self.dyn_ltree[lcode + LITERALS + 1].fc += 1;
            let dc = tables().d_code(d);
            self.dyn_dtree[dc].fc += 1;
        }
        self.syms.len() == SYM_END
    }

    fn build_bl_tree(&mut self) -> usize {
        scan_tree(&mut self.dyn_ltree[..], self.l_max_code, &mut self.bl_tree);
        scan_tree(&mut self.dyn_dtree, self.d_max_code, &mut self.bl_tree);
        let bl_desc = StaticDesc {
            stree: None,
            extra: &EXTRA_BLBITS,
            extra_base: 0,
            elems: BL_CODES,
            max_length: MAX_BL_BITS,
        };
        self.builder.build_tree(&mut self.bl_tree, &bl_desc);
        let mut max_blindex = BL_CODES - 1;
        while max_blindex >= 3 {
            if self.bl_tree[BL_ORDER[max_blindex]].dl != 0 {
                break;
            }
            max_blindex -= 1;
        }
        self.builder.opt_len = self
            .builder
            .opt_len
            .wrapping_add(3 * (max_blindex as u64 + 1) + 5 + 5 + 4);
        max_blindex
    }

    fn compress_block(&self, bits: &mut Bits, ltree: &[Ct], dtree: &[Ct]) {
        let t = tables();
        for &(dist, lc) in &self.syms {
            let (mut dist, mut lc) = (dist as usize, lc as i32);
            if dist == 0 {
                bits.send_code(lc as usize, ltree);
            } else {
                let code = t.length_code[lc as usize] as usize;
                bits.send_code(code + LITERALS + 1, ltree);
                let extra = EXTRA_LBITS[code] as i32;
                if extra != 0 {
                    lc -= t.base_length[code];
                    bits.send_bits(lc, extra);
                }
                dist -= 1;
                let code = t.d_code(dist);
                bits.send_code(code, dtree);
                let extra = EXTRA_DBITS[code] as i32;
                if extra != 0 {
                    bits.send_bits(dist as i32 - t.base_dist[code], extra);
                }
            }
        }
        bits.send_code(END_BLOCK, ltree);
    }

    /// zlib `_tr_flush_block` (level > 0).
    fn tr_flush_block(
        &mut self,
        buf: Option<usize>,
        stored_len: u64,
        last: bool,
        out: &mut Vec<u8>,
    ) {
        let tables = tables();
        let l_desc = StaticDesc {
            stree: Some(&tables.static_ltree[..]),
            extra: &EXTRA_LBITS,
            extra_base: LITERALS + 1,
            elems: L_CODES,
            max_length: MAX_BITS,
        };
        let d_desc = StaticDesc {
            stree: Some(&tables.static_dtree[..]),
            extra: &EXTRA_DBITS,
            extra_base: 0,
            elems: D_CODES,
            max_length: MAX_BITS,
        };
        self.l_max_code = self.builder.build_tree(&mut self.dyn_ltree[..], &l_desc);
        self.d_max_code = self.builder.build_tree(&mut self.dyn_dtree, &d_desc);
        let max_blindex = self.build_bl_tree();
        let mut opt_lenb = self.builder.opt_len.wrapping_add(3 + 7) >> 3;
        let static_lenb = self.builder.static_len.wrapping_add(3 + 7) >> 3;
        if static_lenb <= opt_lenb {
            opt_lenb = static_lenb;
        }
        let mut bits = Bits {
            out,
            buf: self.bi_buf,
            valid: self.bi_valid,
        };
        let last_bit = last as i32;
        if let Some(start) = buf.filter(|_| stored_len + 4 <= opt_lenb) {
            bits.send_bits((STORED_BLOCK << 1) + last_bit, 3);
            bits.windup();
            bits.put_short(stored_len as u16);
            bits.put_short(!(stored_len as u16));
            bits.out
                .extend_from_slice(&self.window[start..start + stored_len as usize]);
        } else if static_lenb == opt_lenb {
            bits.send_bits((STATIC_TREES << 1) + last_bit, 3);
            self.compress_block(&mut bits, &tables.static_ltree, &tables.static_dtree);
        } else {
            bits.send_bits((DYN_TREES << 1) + last_bit, 3);
            let (lcodes, dcodes, blcodes) = (
                self.l_max_code + 1,
                self.d_max_code + 1,
                max_blindex as i32 + 1,
            );
            bits.send_bits(lcodes - 257, 5);
            bits.send_bits(dcodes - 1, 5);
            bits.send_bits(blcodes - 4, 4);
            for &o in &BL_ORDER[..blcodes as usize] {
                bits.send_bits(self.bl_tree[o].dl as i32, 3);
            }
            send_tree(&mut bits, &self.dyn_ltree[..], lcodes - 1, &self.bl_tree);
            send_tree(&mut bits, &self.dyn_dtree, dcodes - 1, &self.bl_tree);
            self.compress_block(&mut bits, &self.dyn_ltree[..], &self.dyn_dtree);
        }
        if last {
            bits.windup();
        }
        self.bi_buf = bits.buf;
        self.bi_valid = bits.valid;
        self.init_block();
    }

    /// FLUSH_BLOCK; returns false if the output has reached `capacity`.
    fn flush_block(&mut self, last: bool, out: &mut Vec<u8>, base: usize, capacity: usize) -> bool {
        let buf = (self.block_start >= 0).then_some(self.block_start as usize);
        let stored_len = (self.strstart as isize - self.block_start) as u64;
        self.tr_flush_block(buf, stored_len, last, out);
        self.block_start = self.strstart as isize;
        out.len() - base < capacity
    }

    /// One zlib `deflateReset` + `deflate(Z_FINISH)` call over `input` with an
    /// output buffer of `capacity` bytes. On [`Finish::OutOfSpace`] the bytes
    /// appended to `out` are a partial stream that the caller must discard.
    pub fn deflate_finish(&mut self, input: &[u8], out: &mut Vec<u8>, capacity: usize) -> Finish {
        self.reset();
        self.next_in = 0;
        self.avail_in = input.len();
        if self.lazy {
            return self.deflate_slow(input, out, capacity);
        }
        let base = out.len();
        loop {
            if self.lookahead < MIN_LOOKAHEAD {
                self.fill_window(input);
                if self.lookahead == 0 {
                    break;
                }
            }
            let mut hash_head = 0;
            if self.lookahead >= MIN_MATCH {
                hash_head = self.insert_string(self.strstart);
            }
            if hash_head != 0 && self.strstart - hash_head <= MAX_DIST {
                self.match_length = self.longest_match(hash_head, MIN_MATCH - 1);
            }
            let bflush;
            if self.match_length >= MIN_MATCH {
                bflush = self.tally(
                    self.strstart - self.match_start,
                    self.match_length - MIN_MATCH,
                );
                self.lookahead -= self.match_length;
                if self.match_length <= self.cfg.max_insert_length && self.lookahead >= MIN_MATCH {
                    self.match_length -= 1;
                    loop {
                        self.strstart += 1;
                        self.insert_string(self.strstart);
                        self.match_length -= 1;
                        if self.match_length == 0 {
                            break;
                        }
                    }
                    self.strstart += 1;
                } else {
                    self.strstart += self.match_length;
                    self.match_length = 0;
                    self.ins_h = self.window[self.strstart] as u32;
                    self.update_hash(self.window[self.strstart + 1]);
                }
            } else {
                bflush = self.tally(0, self.window[self.strstart] as usize);
                self.lookahead -= 1;
                self.strstart += 1;
            }
            if bflush && !self.flush_block(false, out, base, capacity) {
                return Finish::OutOfSpace;
            }
        }
        self.insert = self.strstart.min(MIN_MATCH - 1);
        if self.flush_block(true, out, base, capacity) {
            Finish::Done
        } else {
            Finish::OutOfSpace
        }
    }

    /// zlib `deflate_slow` with `flush == Z_FINISH` (default strategy).
    fn deflate_slow(&mut self, input: &[u8], out: &mut Vec<u8>, capacity: usize) -> Finish {
        let base = out.len();
        loop {
            if self.lookahead < MIN_LOOKAHEAD {
                self.fill_window(input);
                if self.lookahead == 0 {
                    break;
                }
            }
            let mut hash_head = 0;
            if self.lookahead >= MIN_MATCH {
                hash_head = self.insert_string(self.strstart);
            }
            self.prev_length = self.match_length;
            self.prev_match = self.match_start;
            self.match_length = MIN_MATCH - 1;
            if hash_head != 0
                && self.prev_length < self.cfg.max_insert_length
                && self.strstart - hash_head <= MAX_DIST
            {
                self.match_length = self.longest_match(hash_head, self.prev_length);
                if self.match_length <= 5
                    && self.match_length == MIN_MATCH
                    && self.strstart - self.match_start > TOO_FAR
                {
                    self.match_length = MIN_MATCH - 1;
                }
            }
            if self.prev_length >= MIN_MATCH && self.match_length <= self.prev_length {
                let max_insert = self.strstart + self.lookahead - MIN_MATCH;
                let bflush = self.tally(
                    self.strstart - 1 - self.prev_match,
                    self.prev_length - MIN_MATCH,
                );
                self.lookahead -= self.prev_length - 1;
                self.prev_length -= 2;
                loop {
                    self.strstart += 1;
                    if self.strstart <= max_insert {
                        self.insert_string(self.strstart);
                    }
                    self.prev_length -= 1;
                    if self.prev_length == 0 {
                        break;
                    }
                }
                self.match_available = false;
                self.match_length = MIN_MATCH - 1;
                self.strstart += 1;
                if bflush && !self.flush_block(false, out, base, capacity) {
                    return Finish::OutOfSpace;
                }
            } else if self.match_available {
                let bflush = self.tally(0, self.window[self.strstart - 1] as usize);
                if bflush && !self.flush_block(false, out, base, capacity) {
                    return Finish::OutOfSpace;
                }
                self.strstart += 1;
                self.lookahead -= 1;
            } else {
                self.match_available = true;
                self.strstart += 1;
                self.lookahead -= 1;
            }
        }
        if self.match_available {
            self.tally(0, self.window[self.strstart - 1] as usize);
            self.match_available = false;
        }
        self.insert = self.strstart.min(MIN_MATCH - 1);
        if self.flush_block(true, out, base, capacity) {
            Finish::Done
        } else {
            Finish::OutOfSpace
        }
    }
}

/// A whole-buffer zlib deflater for any level, stateless across calls
/// (fresh `deflateInit2` per stream), so its output depends only on the input.
pub struct ZlibDeflater {
    zlib: ZlibFast,
}

impl ZlibDeflater {
    pub fn new(level: u32) -> ZlibDeflater {
        ZlibDeflater {
            zlib: ZlibFast::new(level),
        }
    }

    /// Compress `input` as one complete raw DEFLATE stream onto `out`.
    pub fn compress(&mut self, input: &[u8], out: &mut Vec<u8>) {
        self.zlib.reinit();
        let done = self.zlib.deflate_finish(input, out, usize::MAX);
        debug_assert_eq!(done, Finish::Done);
    }
}

/// The static trees and code tables (zlib builds them once per process too).
fn tables() -> &'static Tables {
    static TABLES: std::sync::OnceLock<Box<Tables>> = std::sync::OnceLock::new();
    TABLES.get_or_init(Tables::new)
}

fn send_tree(bits: &mut impl BitSink, tree: &[Ct], max_code: i32, bl_tree: &[Ct]) {
    let mut prevlen: i32 = -1;
    let mut nextlen = tree[0].dl as i32;
    let mut count = 0;
    let (mut max_count, mut min_count) = if nextlen == 0 { (138, 3) } else { (7, 4) };
    for n in 0..=max_code {
        let curlen = nextlen;
        nextlen = tree[(n + 1) as usize].dl as i32;
        count += 1;
        if count < max_count && curlen == nextlen {
            continue;
        } else if count < min_count {
            loop {
                bits.send_code(curlen as usize, bl_tree);
                count -= 1;
                if count == 0 {
                    break;
                }
            }
        } else if curlen != 0 {
            if curlen != prevlen {
                bits.send_code(curlen as usize, bl_tree);
                count -= 1;
            }
            bits.send_code(REP_3_6, bl_tree);
            bits.send_bits(count - 3, 2);
        } else if count <= 10 {
            bits.send_code(REPZ_3_10, bl_tree);
            bits.send_bits(count - 3, 3);
        } else {
            bits.send_code(REPZ_11_138, bl_tree);
            bits.send_bits(count - 11, 7);
        }
        count = 0;
        prevlen = curlen;
        (max_count, min_count) = if nextlen == 0 {
            (138, 3)
        } else if curlen == nextlen {
            (6, 3)
        } else {
            (7, 4)
        };
    }
}

/// htsjdk's BGZF block compressor: zlib level 2 into a 65,518-byte buffer,
/// falling back to level 0 (one stored block) when it does not fit.
pub struct HtsjdkDeflater {
    zlib: ZlibFast,
    scratch: Vec<u8>,
}

impl HtsjdkDeflater {
    pub fn new() -> Self {
        HtsjdkDeflater {
            zlib: ZlibFast::new(2),
            scratch: Vec::with_capacity(1 << 17),
        }
    }
}

impl HtsjdkDeflater {
    /// Fresh `deflateInit2` state (see [`ZlibFast::reinit`]).
    pub fn reinit(&mut self) {
        self.zlib.reinit();
    }
}

impl Default for HtsjdkDeflater {
    fn default() -> Self {
        Self::new()
    }
}

impl Deflate for HtsjdkDeflater {
    fn deflate(&mut self, input: &[u8], out: &mut Vec<u8>) {
        self.scratch.clear();
        match self
            .zlib
            .deflate_finish(input, &mut self.scratch, HTSJDK_OUTPUT_CAPACITY)
        {
            Finish::Done => out.extend_from_slice(&self.scratch),
            Finish::OutOfSpace => Stored.deflate(input, out),
        }
    }
}
