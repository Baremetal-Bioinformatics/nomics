//! The BWA FM-index, read from `bwa index` output (`<prefix>.64.{bwt,sa,pac,ann,amb}`).
//!
//! Layouts follow bwa 0.7.19 `bwt.c` / `bntseq.c`, little-endian:
//!
//! - `.bwt`: `primary: u64`, `L2[1..=4]: u64`, then `u32` words. Every
//!   `OCC_INTERVAL = 128` BWT positions take 16 words: four `u64` occurrence
//!   counts (A, C, G, T before the block) followed by 8 words of 2-bit bases,
//!   most significant pair first. The `$` row (`primary`) is not stored.
//! - `.sa`: `primary`, 4 skipped `u64`, `sa_intv`, `seq_len`, then
//!   `SA[intv], SA[2*intv], ...`. `SA[0]` is implicitly `-1` (it would be `seq_len`).
//! - `.pac`: the forward strand, 4 bases per byte, first base in the high bits;
//!   `l_pac / 4 + 1` bytes are read. Ambiguous bases were replaced by random ones
//!   at index time; their runs are listed in `.amb`.
//! - `.ann` / `.amb`: text, read with the same `fscanf` patterns as bwa.
//!
//! The BWT is of the forward strand followed by its reverse complement, so
//! `seq_len = 2 * l_pac` and SA values index that concatenation.
//!
//! bwa also reads `<prefix>.alt` when present to mark ALT contigs; we refuse to
//! load an index that has one until ALT handling is implemented.

use std::fs::File;
use std::io::{self, Read};
use std::path::{Path, PathBuf};

use nomics_mem::HostBuf;

const OCC_INTV_SHIFT: u32 = 7;
const OCC_INTERVAL: u64 = 1 << OCC_INTV_SHIFT;
const OCC_INTV_MASK: u64 = OCC_INTERVAL - 1;
/// `u32` words per 128-base block: 4 `u64` counts + 8 words of bases.
const BLOCK_WORDS: u64 = 16;

fn invalid(msg: String) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg)
}

fn read_u64s<const N: usize>(r: &mut impl Read) -> io::Result<[u64; N]> {
    let mut out = [0u64; N];
    let mut b = [0u8; 8];
    for v in &mut out {
        r.read_exact(&mut b)?;
        *v = u64::from_le_bytes(b);
    }
    Ok(out)
}

/// Read `len` bytes from `r` into a fresh huge-page buffer.
fn read_into_hostbuf(r: &mut impl Read, len: usize, offset: usize, label: &str) -> io::Result<HostBuf> {
    let mut buf = HostBuf::new(len)?;
    buf.set_label(label);
    r.read_exact(&mut buf.as_mut_slice()[offset..])?;
    Ok(buf)
}

/// bwa's `cnt_table`: byte of 4 bases -> per-base counts packed one per byte.
const fn cnt_table() -> [u32; 256] {
    let mut t = [0u32; 256];
    let mut i = 0;
    while i < 256 {
        let mut x = 0u32;
        let mut j = 0;
        while j < 4 {
            let n = ((i & 3) == j) as u32
                + ((i >> 2 & 3) == j) as u32
                + ((i >> 4 & 3) == j) as u32
                + ((i >> 6) == j) as u32;
            x |= n << (j << 3);
            j += 1;
        }
        t[i] = x;
        i += 1;
    }
    t
}
static CNT_TABLE: [u32; 256] = cnt_table();

#[inline]
fn occ_aux4(b: u32) -> u64 {
    (CNT_TABLE[(b & 0xff) as usize]
        + CNT_TABLE[(b >> 8 & 0xff) as usize]
        + CNT_TABLE[(b >> 16 & 0xff) as usize]
        + CNT_TABLE[(b >> 24) as usize]) as u64
}

/// bwa `__occ_aux`: occurrences of `c` among the 32 bases of `y`.
#[inline]
fn occ_aux(y: u64, c: u8) -> u64 {
    let y = (if c & 2 != 0 { y } else { !y }) >> 1
        & (if c & 1 != 0 { y } else { !y })
        & 0x5555_5555_5555_5555;
    y.count_ones() as u64
}

/// BWT string with interleaved occurrence counts.
pub struct Bwt {
    pub primary: u64,
    /// `L2[c]`: number of bases smaller than `c`; `L2[4] = seq_len`.
    pub l2: [u64; 5],
    pub seq_len: u64,
    /// Registered with the GPU backend for the index's life (Metal binds the
    /// 3.2 GB once instead of on every seed launch; free elsewhere).
    words: nomics_gpu::GpuHostBuf,
}

/// bwa's `(bwtint_t)-1`, used as an empty-interval sentinel.
pub const NONE: u64 = u64::MAX;

impl Bwt {
    pub fn load(path: &Path) -> io::Result<Self> {
        let mut f = File::open(path)?;
        let size = f.metadata()?.len();
        if size < 40 || (size - 40) % 4 != 0 {
            return Err(invalid(format!("{}: bad size {size}", path.display())));
        }
        let [primary, l1, l2, l3, l4] = read_u64s::<5>(&mut f)?;
        let l2 = [0, l1, l2, l3, l4];
        let seq_len = l4;
        let n_words = (size - 40) / 4;
        // bwt_bwtupdate_core: one count block per 128 positions plus a trailing
        // one, and the 2-bit bases packed 16 per word.
        let expect = (seq_len.div_ceil(OCC_INTERVAL) + 1) * 8 + seq_len.div_ceil(16);
        if primary > seq_len || !l2.windows(2).all(|w| w[0] <= w[1]) || n_words != expect {
            return Err(invalid(format!(
                "{}: inconsistent header (primary {primary}, L2 {l2:?}, {n_words} words)",
                path.display()
            )));
        }
        let words = read_into_hostbuf(&mut f, (n_words * 4) as usize, 0, "index BWT")?;
        Ok(Self {
            primary,
            l2,
            seq_len,
            words: nomics_gpu::GpuHostBuf::register(words),
        })
    }

    /// The interleaved occurrence/base words (64-byte-aligned blocks).
    pub fn words(&self) -> &[u32] {
        self.words.as_u32s()
    }

    #[inline]
    fn w(&self) -> &[u32] {
        self.words.as_u32s()
    }

    #[inline]
    fn block(&self, k: u64) -> usize {
        ((k >> OCC_INTV_SHIFT) * BLOCK_WORDS) as usize
    }

    #[inline]
    fn block_count(&self, base: usize, c: usize) -> u64 {
        let w = self.w();
        w[base + 2 * c] as u64 | (w[base + 2 * c + 1] as u64) << 32
    }

    /// `bwt_B0`: the base at stored position `k` (after removing `$`).
    #[inline]
    pub fn base(&self, k: u64) -> u8 {
        let word = self.w()[self.block(k) + 8 + ((k & 0x7f) >> 4) as usize];
        (word >> ((!k & 0xf) << 1) & 3) as u8
    }

    /// `bwt_occ`: occurrences of `c` in BWT[0..=k].
    pub fn occ(&self, k: u64, c: u8) -> u64 {
        if k == self.seq_len {
            return self.l2[c as usize + 1] - self.l2[c as usize];
        }
        if k == NONE {
            return 0;
        }
        let k = k - (k >= self.primary) as u64;
        let base = self.block(k);
        let mut n = self.block_count(base, c as usize);
        let w = self.w();
        let mut p = base + 8;
        let end = p + (((k >> 5) - ((k & !OCC_INTV_MASK) >> 5)) << 1) as usize;
        while p < end {
            n += occ_aux((w[p] as u64) << 32 | w[p + 1] as u64, c);
            p += 2;
        }
        let y = ((w[p] as u64) << 32 | w[p + 1] as u64) & !((1u64 << ((!k & 31) << 1)) - 1);
        n += occ_aux(y, c);
        if c == 0 {
            n -= !k & 31;
        }
        n
    }

    /// `bwt_occ4`: occurrences of each base in BWT[0..=k].
    pub fn occ4(&self, k: u64) -> [u64; 4] {
        if k == NONE {
            return [0; 4];
        }
        let k = k - (k >= self.primary) as u64;
        let base = self.block(k);
        let mut cnt = [0u64; 4];
        for (c, v) in cnt.iter_mut().enumerate() {
            *v = self.block_count(base, c);
        }
        let w = self.w();
        let mut p = base + 8;
        let end = p + ((k >> 4) - ((k & !OCC_INTV_MASK) >> 4)) as usize;
        let mut x = 0u64;
        while p < end {
            x += occ_aux4(w[p]);
            p += 1;
        }
        let tmp = w[p] & !((1u32 << ((!k & 15) << 1)).wrapping_sub(1));
        x += occ_aux4(tmp) - (!k & 15);
        cnt[0] += x & 0xff;
        cnt[1] += x >> 8 & 0xff;
        cnt[2] += x >> 16 & 0xff;
        cnt[3] += x >> 24;
        cnt
    }

    /// `bwt_2occ4`, requiring `k <= l`.
    pub fn occ4_pair(&self, k: u64, l: u64) -> ([u64; 4], [u64; 4]) {
        // Same results as two `occ4` calls; bwa's shared-block shortcut is
        // purely an optimisation.
        (self.occ4(k), self.occ4(l))
    }

    /// `bwt_invPsi`: the row of the suffix one position earlier in the text.
    #[inline]
    pub fn inv_psi(&self, k: u64) -> u64 {
        if k == self.primary {
            return 0;
        }
        let x = k - (k > self.primary) as u64;
        let c = self.base(x);
        self.l2[c as usize] + self.occ(k, c)
    }
}

/// Sampled suffix array.
pub struct Sa {
    pub intv: u64,
    /// `entries[i] = SA[i * intv]`, with `entries[0] = -1`.
    entries: HostBuf,
}

impl Sa {
    pub fn load(path: &Path, bwt: &Bwt) -> io::Result<Self> {
        let mut f = File::open(path)?;
        let [primary, _, _, _, _, intv, seq_len] = read_u64s::<7>(&mut f)?;
        if primary != bwt.primary || seq_len != bwt.seq_len {
            return Err(invalid(format!(
                "{}: SA-BWT inconsistency (primary {primary} vs {}, seq_len {seq_len} vs {})",
                path.display(),
                bwt.primary,
                bwt.seq_len
            )));
        }
        if intv == 0 || !intv.is_power_of_two() {
            return Err(invalid(format!(
                "{}: bad SA interval {intv}",
                path.display()
            )));
        }
        let n_sa = (seq_len + intv) / intv;
        let expect = 56 + (n_sa - 1) * 8;
        let size = f.metadata()?.len();
        if size != expect {
            return Err(invalid(format!(
                "{}: size {size}, expected {expect}",
                path.display()
            )));
        }
        let mut entries = read_into_hostbuf(&mut f, (n_sa * 8) as usize, 8, "index SA")?;
        entries.as_u64s_mut()[0] = NONE;
        Ok(Self { intv, entries })
    }

    /// `bwt_sa`: the text position of BWT row `k`.
    pub fn sa(&self, bwt: &Bwt, mut k: u64) -> u64 {
        let mask = self.intv - 1;
        let mut sa = 0u64;
        while k & mask != 0 {
            sa += 1;
            k = bwt.inv_psi(k);
        }
        sa.wrapping_add(self.entries.as_u64s()[(k / self.intv) as usize])
    }

    /// bwa's samples, `entries()[i] = SA[i * intv]` (`SA[0]` = -1).
    pub fn entries(&self) -> &[u64] {
        self.entries.as_u64s()
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Ann {
    pub offset: i64,
    pub len: i32,
    pub n_ambs: i32,
    pub gi: u32,
    pub is_alt: bool,
    pub name: String,
    pub anno: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Amb {
    pub offset: i64,
    pub len: i32,
    pub amb: u8,
}

/// `bntseq_t` without the file handle.
#[derive(Debug, Clone)]
pub struct Bns {
    pub l_pac: i64,
    pub seed: u32,
    pub anns: Vec<Ann>,
    pub ambs: Vec<Amb>,
}

/// Whitespace-token reader mimicking `fscanf` `%d`/`%s` on a byte slice.
struct Scan<'a> {
    s: &'a [u8],
    i: usize,
}

impl<'a> Scan<'a> {
    fn skip_ws(&mut self) {
        while self.i < self.s.len() && self.s[self.i].is_ascii_whitespace() {
            self.i += 1;
        }
    }

    fn token(&mut self) -> Option<&'a [u8]> {
        self.skip_ws();
        let start = self.i;
        while self.i < self.s.len() && !self.s[self.i].is_ascii_whitespace() {
            self.i += 1;
        }
        (self.i > start).then(|| &self.s[start..self.i])
    }

    fn num<T: std::str::FromStr>(&mut self) -> Option<T> {
        std::str::from_utf8(self.token()?).ok()?.parse().ok()
    }

    /// Rest of the current line, without the newline, consuming it.
    fn rest_of_line(&mut self) -> Option<&'a [u8]> {
        let start = self.i;
        while self.i < self.s.len() && self.s[self.i] != b'\n' {
            self.i += 1;
        }
        if self.i == self.s.len() {
            return None; // bwa treats EOF before the newline as an error
        }
        self.i += 1;
        Some(&self.s[start..self.i - 1])
    }
}

impl Bns {
    pub fn load(ann: &Path, amb: &Path) -> io::Result<Self> {
        let bad = |p: &Path| invalid(format!("parse error reading {}", p.display()));
        let text = std::fs::read(ann)?;
        let mut s = Scan { s: &text, i: 0 };
        let (l_pac, n_seqs, seed): (i64, i32, u32) = (
            s.num().ok_or_else(|| bad(ann))?,
            s.num().ok_or_else(|| bad(ann))?,
            s.num().ok_or_else(|| bad(ann))?,
        );
        let mut anns = Vec::with_capacity(n_seqs.max(0) as usize);
        for _ in 0..n_seqs {
            let gi: u32 = s.num().ok_or_else(|| bad(ann))?;
            let name = String::from_utf8_lossy(s.token().ok_or_else(|| bad(ann))?).into_owned();
            let line = s.rest_of_line().ok_or_else(|| bad(ann))?;
            let line = &line[..line.len().min(8191)];
            let anno = if line.len() > 1 && line != b" (null)" {
                String::from_utf8_lossy(&line[1..]).into_owned()
            } else {
                String::new()
            };
            anns.push(Ann {
                gi,
                name,
                anno,
                offset: s.num().ok_or_else(|| bad(ann))?,
                len: s.num().ok_or_else(|| bad(ann))?,
                n_ambs: s.num().ok_or_else(|| bad(ann))?,
                is_alt: false,
            });
        }
        let text = std::fs::read(amb)?;
        let mut s = Scan { s: &text, i: 0 };
        let (amb_l_pac, amb_n_seqs, n_holes): (i64, i32, i32) = (
            s.num().ok_or_else(|| bad(amb))?,
            s.num().ok_or_else(|| bad(amb))?,
            s.num().ok_or_else(|| bad(amb))?,
        );
        if amb_l_pac != l_pac || amb_n_seqs != n_seqs {
            return Err(invalid("inconsistent .ann and .amb files".into()));
        }
        let mut ambs = Vec::with_capacity(n_holes.max(0) as usize);
        for _ in 0..n_holes {
            ambs.push(Amb {
                offset: s.num().ok_or_else(|| bad(amb))?,
                len: s.num().ok_or_else(|| bad(amb))?,
                amb: s.token().ok_or_else(|| bad(amb))?[0],
            });
        }
        Ok(Self {
            l_pac,
            seed,
            anns,
            ambs,
        })
    }

    /// `bns_depos`: forward coordinate of a concatenated-space position, and
    /// whether it was on the reverse strand.
    #[inline]
    pub fn depos(&self, pos: i64) -> (i64, bool) {
        if pos >= self.l_pac {
            ((self.l_pac << 1) - 1 - pos, true)
        } else {
            (pos, false)
        }
    }

    /// `bns_pos2rid`: contig of a forward position, or -1 past the end.
    pub fn pos2rid(&self, pos_f: i64) -> i32 {
        if pos_f >= self.l_pac {
            return -1;
        }
        let (mut left, mut mid, mut right) = (0usize, 0usize, self.anns.len());
        while left < right {
            mid = (left + right) >> 1;
            if pos_f >= self.anns[mid].offset {
                if mid == self.anns.len() - 1 || pos_f < self.anns[mid + 1].offset {
                    break;
                }
                left = mid + 1;
            } else {
                right = mid;
            }
        }
        mid as i32
    }

    /// `bns_intv2rid`: contig containing `[rb, re)`, -1 if it bridges
    /// contigs, -2 if it crosses the strand boundary.
    pub fn intv2rid(&self, rb: i64, re: i64) -> i32 {
        if rb < self.l_pac && re > self.l_pac {
            return -2;
        }
        debug_assert!(rb <= re);
        let rid_b = self.pos2rid(self.depos(rb).0);
        let rid_e = if rb < re {
            self.pos2rid(self.depos(re - 1).0)
        } else {
            rid_b
        };
        if rid_b == rid_e { rid_b } else { -1 }
    }
}

/// 2-bit packed forward strand.
pub struct Pac {
    bytes: HostBuf,
}

impl Pac {
    pub fn load(path: &Path, l_pac: i64) -> io::Result<Self> {
        let len = (l_pac / 4 + 1) as usize;
        let mut f = File::open(path)?;
        let size = f.metadata()?.len() as usize;
        if size < len {
            return Err(invalid(format!(
                "{}: {size} bytes, need {len}",
                path.display()
            )));
        }
        Ok(Self {
            bytes: read_into_hostbuf(&mut f, len, 0, "index pac")?,
        })
    }

    /// `_get_pac`: base `k` of the forward strand (0..4 = A, C, G, T).
    #[inline]
    pub fn get(&self, k: u64) -> u8 {
        self.bytes.as_slice()[(k >> 2) as usize] >> ((!k & 3) << 1) & 3
    }

    /// The packed bytes (`l_pac / 4 + 1`).
    pub fn bytes(&self) -> &[u8] {
        self.bytes.as_slice()
    }
}

/// The four bases of each packed byte, first base in the high bits.
const UNPACK: [[u8; 4]; 256] = {
    let mut t = [[0u8; 4]; 256];
    let mut b = 0;
    while b < 256 {
        t[b] = [
            (b >> 6) as u8 & 3,
            (b >> 4) as u8 & 3,
            (b >> 2) as u8 & 3,
            b as u8 & 3,
        ];
        b += 1;
    }
    t
};

/// `_get_pac(pac, k)` for `k` in `[beg, end)`, appended to `out`.
pub(crate) fn unpack_2bit(pac: &[u8], beg: u64, end: u64, out: &mut Vec<u8>) {
    if beg >= end {
        return;
    }
    let get = |k: u64| pac[(k >> 2) as usize] >> ((!k & 3) << 1) & 3;
    let start = out.len();
    out.reserve((end - beg) as usize);
    let mut k = beg;
    while k < end && k & 3 != 0 {
        out.push(get(k));
        k += 1;
    }
    let full = end & !3;
    if k < full {
        let src = &pac[(k >> 2) as usize..(full >> 2) as usize];
        let at = out.len();
        out.resize(at + src.len() * 4, 0);
        for (dst, &b) in out[at..].as_chunks_mut::<4>().0.iter_mut().zip(src) {
            dst.copy_from_slice(&UNPACK[b as usize]);
        }
        k = full;
    }
    while k < end {
        out.push(get(k));
        k += 1;
    }
    debug_assert_eq!(out.len() - start, (end - beg) as usize);
}

/// Everything `bwa_idx_load(prefix, BWA_IDX_ALL)` loads.
pub struct Index {
    pub bwt: Bwt,
    pub sa: Sa,
    pub bns: Bns,
    pub pac: Pac,
    /// Optional denser SA (see `dense_sa`); same results, fewer LF steps.
    pub dense: Option<crate::dense_sa::DenseSa>,
}

impl Index {
    /// Load from a FASTA path or index prefix, preferring `<hint>.64.*` like
    /// `bwa_idx_infer_prefix`.
    pub fn load(hint: &Path) -> io::Result<Self> {
        let with = |p: &Path, ext: &str| {
            let mut s = p.as_os_str().to_owned();
            s.push(ext);
            PathBuf::from(s)
        };
        let prefix = if with(hint, ".64.bwt").exists() {
            with(hint, ".64")
        } else if with(hint, ".bwt").exists() {
            hint.to_path_buf()
        } else {
            return Err(io::Error::new(
                io::ErrorKind::NotFound,
                format!("no BWA index at {}", hint.display()),
            ));
        };
        if with(&prefix, ".alt").exists() {
            return Err(invalid(format!(
                "{}.alt present: ALT-aware alignment is not implemented",
                prefix.display()
            )));
        }
        let bwt = Bwt::load(&with(&prefix, ".bwt"))?;
        let sa = Sa::load(&with(&prefix, ".sa"), &bwt)?;
        let bns = Bns::load(&with(&prefix, ".ann"), &with(&prefix, ".amb"))?;
        if bwt.seq_len != 2 * bns.l_pac as u64 {
            return Err(invalid(format!(
                "BWT length {} is not twice l_pac {}",
                bwt.seq_len, bns.l_pac
            )));
        }
        let pac = Pac::load(&with(&prefix, ".pac"), bns.l_pac)?;
        Ok(Self {
            bwt,
            sa,
            bns,
            pac,
            dense: None,
        })
    }

    /// Load (or build and cache) the legacy interval-4 dense SA at `path`.
    ///
    /// New callers that select an interval should use
    /// [`Index::attach_dense_sa_with_interval`].
    pub fn attach_dense_sa(&mut self, path: &Path, threads: usize) -> io::Result<()> {
        self.attach_dense_sa_with_interval(path, 4, threads)
    }

    /// Load (or build and cache) an interval-1, -2, or -4 dense SA at `path`.
    pub fn attach_dense_sa_with_interval(
        &mut self,
        path: &Path,
        interval: u64,
        threads: usize,
    ) -> io::Result<()> {
        self.dense = Some(crate::dense_sa::DenseSa::load_or_build(
            path, &self.bwt, &self.sa, interval, threads,
        )?);
        Ok(())
    }

    /// `bwt_sa(bwt, k)`, through the dense SA when attached.
    #[inline]
    pub fn sa_at(&self, k: u64) -> u64 {
        match &self.dense {
            Some(d) => d.sa(&self.bwt, k),
            None => self.sa.sa(&self.bwt, k),
        }
    }

    /// Promote the large buffers to 2 MiB pages where the kernel left them
    /// on 4 KiB ones (see [`nomics_mem::HostBuf::collapse_huge`]; failure is
    /// reported, not fatal), then lock them in RAM
    /// (see [`nomics_mem::HostBuf::lock`]).
    pub fn lock_memory(&self) -> io::Result<()> {
        let (total, before) = self.memory();
        if before.is_some_and(|h| h < total) {
            let started = std::time::Instant::now();
            let mut failed = None;
            for buf in [&self.bwt.words, &self.sa.entries, &self.pac.bytes] {
                failed = failed.or(buf.collapse_huge().err());
            }
            if let Some(dense) = &self.dense {
                failed = failed.or(dense.collapse_huge().err());
            }
            let after = self.memory().1.unwrap_or(0);
            eprintln!(
                "index: collapsed to 2 MiB pages: {:.0}% -> {:.0}% in {:.1} s{}",
                100.0 * before.unwrap_or(0) as f64 / total as f64,
                100.0 * after as f64 / total as f64,
                started.elapsed().as_secs_f64(),
                failed.map_or(String::new(), |e| format!(" ({e})"))
            );
        }
        for buf in [&self.bwt.words, &self.sa.entries, &self.pac.bytes] {
            buf.lock()?;
        }
        if let Some(dense) = &self.dense {
            dense.lock()?;
        }
        Ok(())
    }

    /// Total bytes of the large buffers, and how many are on huge pages
    /// (`None` where smaps is unavailable).
    pub fn memory(&self) -> (usize, Option<usize>) {
        let bufs = [&self.bwt.words, &self.sa.entries, &self.pac.bytes];
        let mut total = bufs.iter().map(|b| b.len()).sum();
        let mut huge: Option<usize> = bufs.iter().map(|b| b.huge_page_bytes()).sum();
        if let Some(dense) = &self.dense {
            total += dense.bytes();
            huge = huge.and_then(|n| dense.huge_page_bytes().map(|d| n + d));
        }
        (total, huge)
    }

    /// Base `i` of the BWT text: the forward strand, then its reverse complement.
    #[inline]
    pub fn text(&self, i: u64) -> u8 {
        let l = self.bns.l_pac as u64;
        if i < l {
            self.pac.get(i)
        } else {
            3 - self.pac.get(2 * l - 1 - i)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn tiny_bwt_and_sa(text: &[u8], sa_intv: u64) -> (Bwt, Sa, Vec<u64>) {
        assert!(!text.is_empty() && text.len() < 16);
        let n = text.len() as u64;
        let mut suffixes: Vec<u64> = (0..=n).collect();
        suffixes.sort_by(|&a, &b| text[a as usize..].cmp(&text[b as usize..]));
        assert_eq!(suffixes[0], n);
        let primary = suffixes.iter().position(|&p| p == 0).unwrap() as u64;

        let mut counts = [0u64; 4];
        for &base in text {
            counts[base as usize] += 1;
        }
        let l2 = [
            0,
            counts[0],
            counts[0] + counts[1],
            counts[0] + counts[1] + counts[2],
            n,
        ];
        let mut words = HostBuf::new(17 * 4).unwrap();
        for (row, &pos) in suffixes.iter().enumerate() {
            if row as u64 == primary {
                continue;
            }
            let stored = row as u64 - (row as u64 > primary) as u64;
            let base = if pos == 0 {
                unreachable!()
            } else {
                text[(pos - 1) as usize]
            };
            words.as_u32s_mut()[8] |= (base as u32) << ((!stored & 15) << 1);
        }
        let bwt = Bwt {
            primary,
            l2,
            seq_len: n,
            words: nomics_gpu::GpuHostBuf::register(words),
        };

        let sample_count = (n + sa_intv) / sa_intv;
        let mut entries = HostBuf::new(sample_count as usize * 8).unwrap();
        for (i, value) in entries.as_u64s_mut().iter_mut().enumerate() {
            let pos = suffixes[i * sa_intv as usize];
            *value = if pos == n { NONE } else { pos };
        }
        let expected = suffixes
            .into_iter()
            .map(|pos| if pos == n { NONE } else { pos })
            .collect();
        (
            bwt,
            Sa {
                intv: sa_intv,
                entries,
            },
            expected,
        )
    }

    #[test]
    fn occ_aux_counts_each_base() {
        // 32 bases: A x8, C x8, G x8, T x8.
        let mut y = 0u64;
        for i in 0..32 {
            y = y << 2 | (i / 8) as u64;
        }
        for c in 0..4 {
            assert_eq!(occ_aux(y, c), 8);
        }
    }

    #[test]
    fn cnt_table_matches_bwa() {
        // 0b11_10_01_00 = T G C A: one of each.
        assert_eq!(CNT_TABLE[0b1110_0100], 0x0101_0101);
        assert_eq!(CNT_TABLE[0], 4);
        assert_eq!(CNT_TABLE[0xff], 4 << 24);
    }

    #[test]
    fn dense_sa_intervals_invert_the_same_tiny_index() {
        let (bwt, sa, expected) = tiny_bwt_and_sa(&[0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 3], 4);
        for interval in [1, 2, 4] {
            let dense = crate::dense_sa::DenseSa::build(&bwt, &sa, interval, 3).unwrap();
            assert_eq!(dense.interval(), interval);
            for row in 0..=bwt.seq_len {
                assert_eq!(dense.sa(&bwt, row), expected[row as usize], "row {row}");
                assert_eq!(dense.sa(&bwt, row), sa.sa(&bwt, row), "row {row}");
            }
        }
    }

    #[test]
    fn dense_sa_reads_the_interval_four_v1_cache() {
        let (bwt, sa, expected) = tiny_bwt_and_sa(&[0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 3], 4);
        let unique = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let path = std::env::temp_dir().join(format!(
            "nomics-dense-sa-v1-{}-{unique}.nsa",
            std::process::id()
        ));
        let mut hash = 0xcbf2_9ce4_8422_2325u64;
        for &value in sa.entries() {
            hash = (hash ^ value).wrapping_mul(0x0000_0100_0000_01b3);
        }
        let n = bwt.seq_len / 4 + 1;
        {
            let mut file = File::create(&path).unwrap();
            file.write_all(b"NOMSA4v1").unwrap();
            for word in [bwt.primary, bwt.seq_len, hash, n] {
                file.write_all(&word.to_le_bytes()).unwrap();
            }
            for row in (0..=bwt.seq_len).step_by(4) {
                file.write_all(&expected[row as usize].to_le_bytes()[..5])
                    .unwrap();
            }
        }
        let dense = crate::dense_sa::DenseSa::load_or_build(&path, &bwt, &sa, 4, 2).unwrap();
        for row in 0..=bwt.seq_len {
            assert_eq!(dense.sa(&bwt, row), expected[row as usize]);
        }
        std::fs::remove_file(path).unwrap();
    }
}
