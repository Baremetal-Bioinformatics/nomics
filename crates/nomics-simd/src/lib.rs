//! Hand-written NEON kernels. Calling a `#[target_feature]` intrinsic needs
//! `unsafe` at the boundary even though NEON is part of the aarch64 baseline,
//! so the kernels live in this crate, the only other place besides
//! `nomics-mem` and `nomics-gpu` where `unsafe` is allowed. Every `unsafe`
//! here is a NEON load/store of a whole `[T; N]` array or a call into a
//! `#[target_feature]` function after checking CPU support.
//!
//! Every kernel has a portable counterpart here or in its caller; both must
//! give identical results on all architectures.

#![allow(unsafe_code)]

mod crc32;
mod sha256;
pub mod sw_diag;
pub use crc32::crc32_update;
pub use sha256::{sha256_compress, sha256_compress_portable};

/// BWA 0.7.19 `ksw_u8` / `ksw_i16` (striped Smith-Waterman, SSE2 in BWA):
/// one target row, i.e. the column loop over the `slen` query segments
/// followed by the lazy-F loop.
#[cfg(target_arch = "aarch64")]
pub mod ksw {
    use std::arch::aarch64::*;

    /// Gap penalties as BWA passes them to the kernels.
    #[derive(Clone, Copy)]
    pub struct Gaps<T> {
        pub e_del: T,
        pub oe_del: T,
        pub e_ins: T,
        pub oe_ins: T,
    }

    /// One row of `ksw_u8`: reads `h0` (previous row) and `e`, writes `h1`
    /// and `e`, and returns the row maximum. `scores` is the query profile
    /// row for this target base.
    pub fn row_u8(
        scores: &[[u8; 16]],
        h0: &[[u8; 16]],
        h1: &mut [[u8; 16]],
        e: &mut [[u8; 16]],
        shift: u8,
        gaps: Gaps<u8>,
    ) -> u8 {
        let slen = scores.len();
        assert!(slen > 0 && h0.len() >= slen && h1.len() >= slen && e.len() >= slen);
        // SAFETY: NEON is part of the aarch64 baseline; the slices are
        // checked above.
        unsafe { row_u8_neon(scores, h0, &mut h1[..slen], &mut e[..slen], shift, gaps) }
    }

    #[target_feature(enable = "neon")]
    fn row_u8_neon(
        scores: &[[u8; 16]],
        h0: &[[u8; 16]],
        h1: &mut [[u8; 16]],
        e: &mut [[u8; 16]],
        shift: u8,
        gaps: Gaps<u8>,
    ) -> u8 {
        let slen = scores.len();
        let load = |v: &[u8; 16]| unsafe { vld1q_u8(v.as_ptr()) };
        let store = |dst: &mut [u8; 16], v: uint8x16_t| unsafe { vst1q_u8(dst.as_mut_ptr(), v) };
        let zero = vdupq_n_u8(0);
        let shift = vdupq_n_u8(shift);
        let (e_del, oe_del) = (vdupq_n_u8(gaps.e_del), vdupq_n_u8(gaps.oe_del));
        let (e_ins, oe_ins) = (vdupq_n_u8(gaps.e_ins), vdupq_n_u8(gaps.oe_ins));
        let mut f = zero;
        let mut row_max = zero;
        // _mm_slli_si128(h, 1): lanes move up by one, lane 0 becomes 0.
        let mut h = vextq_u8::<15>(zero, load(&h0[slen - 1]));
        for j in 0..slen {
            let mut v = vqsubq_u8(vqaddq_u8(h, load(&scores[j])), shift);
            let ev = load(&e[j]);
            v = vmaxq_u8(vmaxq_u8(v, ev), f);
            row_max = vmaxq_u8(row_max, v);
            store(
                &mut e[j],
                vmaxq_u8(vqsubq_u8(ev, e_del), vqsubq_u8(v, oe_del)),
            );
            f = vmaxq_u8(vqsubq_u8(f, e_ins), vqsubq_u8(v, oe_ins));
            store(&mut h1[j], v);
            h = load(&h0[j]);
        }
        'lazy_f: for _ in 0..16 {
            f = vextq_u8::<15>(zero, f);
            for cell in h1.iter_mut() {
                let v = vmaxq_u8(load(cell), f);
                store(cell, v);
                let h_gap = vqsubq_u8(v, oe_ins);
                f = vqsubq_u8(f, e_ins);
                if vmaxvq_u8(vcgtq_u8(f, h_gap)) == 0 {
                    break 'lazy_f;
                }
            }
        }
        vmaxvq_u8(row_max)
    }

    /// One row of `ksw_i16`, as [`row_u8`] (no score shift).
    pub fn row_i16(
        scores: &[[i16; 8]],
        h0: &[[i16; 8]],
        h1: &mut [[i16; 8]],
        e: &mut [[i16; 8]],
        gaps: Gaps<u16>,
    ) -> i16 {
        let slen = scores.len();
        assert!(slen > 0 && h0.len() >= slen && h1.len() >= slen && e.len() >= slen);
        // SAFETY: as in `row_u8`.
        unsafe { row_i16_neon(scores, h0, &mut h1[..slen], &mut e[..slen], gaps) }
    }

    #[target_feature(enable = "neon")]
    fn row_i16_neon(
        scores: &[[i16; 8]],
        h0: &[[i16; 8]],
        h1: &mut [[i16; 8]],
        e: &mut [[i16; 8]],
        gaps: Gaps<u16>,
    ) -> i16 {
        let slen = scores.len();
        let load = |v: &[i16; 8]| unsafe { vld1q_s16(v.as_ptr()) };
        let store = |dst: &mut [i16; 8], v: int16x8_t| unsafe { vst1q_s16(dst.as_mut_ptr(), v) };
        // _mm_subs_epu16: unsigned saturating subtraction on signed lanes.
        let subs_u = |a: int16x8_t, b: uint16x8_t| {
            vreinterpretq_s16_u16(vqsubq_u16(vreinterpretq_u16_s16(a), b))
        };
        let zero = vdupq_n_s16(0);
        let (e_del, oe_del) = (vdupq_n_u16(gaps.e_del), vdupq_n_u16(gaps.oe_del));
        let (e_ins, oe_ins) = (vdupq_n_u16(gaps.e_ins), vdupq_n_u16(gaps.oe_ins));
        let mut f = zero;
        let mut row_max = zero;
        let mut h = vextq_s16::<7>(zero, load(&h0[slen - 1]));
        for j in 0..slen {
            let mut v = vqaddq_s16(h, load(&scores[j]));
            let ev = load(&e[j]);
            v = vmaxq_s16(vmaxq_s16(v, ev), f);
            row_max = vmaxq_s16(row_max, v);
            store(&mut e[j], vmaxq_s16(subs_u(ev, e_del), subs_u(v, oe_del)));
            f = vmaxq_s16(subs_u(f, e_ins), subs_u(v, oe_ins));
            store(&mut h1[j], v);
            h = load(&h0[j]);
        }
        'lazy_f: for _ in 0..16 {
            f = vextq_s16::<7>(zero, f);
            for cell in h1.iter_mut() {
                let v = vmaxq_s16(load(cell), f);
                store(cell, v);
                let h_gap = subs_u(v, oe_ins);
                f = subs_u(f, e_ins);
                if vmaxvq_u16(vcgtq_s16(f, h_gap)) == 0 {
                    break 'lazy_f;
                }
            }
        }
        vmaxvq_s16(row_max)
    }
}

/// GATK's Java `LoglessPairHMM` inner loop, up to eight haplotypes per pass
/// (lanes in groups of four, two `float64x2_t` per group and state; with
/// two groups their dependency chains interleave). Each lane performs
/// exactly the scalar recurrence with Java's operand order and separate
/// multiply and add roundings (no fused multiply-add), so results are
/// bit-identical to the portable lanes in `nomics-hc`.
#[cfg(target_arch = "aarch64")]
pub mod pairhmm {
    use std::arch::aarch64::*;

    /// Lanes per pass.
    pub const LANES: usize = 8;

    /// Per-read-base transition and prior values, as `nomics-hc` computes them.
    #[derive(Clone, Copy)]
    pub struct RowParams {
        pub m2m: f64,
        pub m2i: f64,
        pub m2d: f64,
        pub i2m: f64,
        pub i2i: f64,
        pub p_match: f64,
        pub p_mismatch: f64,
    }

    type Row<'a> = (&'a [[f64; LANES]], &'a [[f64; LANES]], &'a [[f64; LANES]]);
    type RowMut<'a> = (
        &'a mut [[f64; LANES]],
        &'a mut [[f64; LANES]],
        &'a mut [[f64; LANES]],
    );

    /// One read row over columns `1..hap.len() + 1`: `prev` = (m, ins, del)
    /// of the previous row, `cur` = this row, both `hap.len() + 1` long;
    /// column 0 of `cur` must already be zero. `hap[j]` holds column
    /// `j + 1`'s haplotype base per lane, widened to `u64`. Only lanes 0..4
    /// are computed when `lanes <= 4`.
    pub fn row(x: u8, p: RowParams, hap: &[[u64; LANES]], prev: Row, cur: RowMut, lanes: usize) {
        let n = hap.len();
        assert!(
            prev.0.len() > n
                && prev.1.len() > n
                && prev.2.len() > n
                && cur.0.len() > n
                && cur.1.len() > n
                && cur.2.len() > n
        );
        // SAFETY: NEON is part of the aarch64 baseline; lengths checked above.
        unsafe {
            if lanes <= 4 {
                row_neon::<1>(x, p, hap, prev, cur)
            } else {
                row_neon::<2>(x, p, hap, prev, cur)
            }
        }
    }

    /// `G` groups of four lanes.
    #[target_feature(enable = "neon")]
    fn row_neon<const G: usize>(
        x: u8,
        p: RowParams,
        hap: &[[u64; LANES]],
        (mp, ip, dp): Row,
        (mc, ic, dc): RowMut,
    ) {
        const H: usize = LANES / 2;
        // Lane pair k of a row: lanes 2k, 2k + 1.
        let ld = |v: &[f64; LANES], k: usize| unsafe { vld1q_f64(v.as_ptr().add(2 * k)) };
        let st = |d: &mut [f64; LANES], k: usize, a: float64x2_t| unsafe {
            vst1q_f64(d.as_mut_ptr().add(2 * k), a)
        };
        let ldh = |v: &[u64; LANES], k: usize| unsafe { vld1q_u64(v.as_ptr().add(2 * k)) };
        let dup = vdupq_n_f64;
        let (m2m, m2i, m2d) = (dup(p.m2m), dup(p.m2i), dup(p.m2d));
        let (i2m, i2i) = (dup(p.i2m), dup(p.i2i));
        let (pm, pmm) = (dup(p.p_match), dup(p.p_mismatch));
        let xv = vdupq_n_u64(x as u64);
        let nv = vdupq_n_u64(b'N' as u64);
        // x == 'N' matches every lane.
        let x_is_n = vdupq_n_u64(if x == b'N' { u64::MAX } else { 0 });
        let pairs = 2 * G;
        debug_assert!(pairs <= H);
        let mut mcl = [dup(0.0); H];
        let mut dcl = [dup(0.0); H];
        for j in 1..=hap.len() {
            for k in 0..pairs {
                let y = ldh(&hap[j - 1], k);
                let hit = vorrq_u64(vorrq_u64(vceqq_u64(y, xv), vceqq_u64(y, nv)), x_is_n);
                let prior = vbslq_f64(hit, pm, pmm);
                // prior * (mp * m2m + ip * i2m + dp * i2m), left to right.
                let mv = vmulq_f64(
                    prior,
                    vaddq_f64(
                        vaddq_f64(
                            vmulq_f64(ld(&mp[j - 1], k), m2m),
                            vmulq_f64(ld(&ip[j - 1], k), i2m),
                        ),
                        vmulq_f64(ld(&dp[j - 1], k), i2m),
                    ),
                );
                let iv = vaddq_f64(vmulq_f64(ld(&mp[j], k), m2i), vmulq_f64(ld(&ip[j], k), i2i));
                // d2d = i2i; uses this row's previous column.
                let dv = vaddq_f64(vmulq_f64(mcl[k], m2d), vmulq_f64(dcl[k], i2i));
                st(&mut mc[j], k, mv);
                st(&mut ic[j], k, iv);
                st(&mut dc[j], k, dv);
                mcl[k] = mv;
                dcl[k] = dv;
            }
        }
    }
}

/// The PairHMM row of `nomics-hc`'s prefix-sharing traversal: lanes are
/// eight reads, each with its own row parameters, against one haplotype base
/// per column. Same per-lane arithmetic as [`pairhmm`] (Java operand order,
/// no FMA), so results are bit-identical to the portable row.
#[cfg(target_arch = "aarch64")]
pub mod pairhmm_reads {
    use std::arch::aarch64::*;

    pub const LANES: usize = 8;
    type V = [f64; LANES];

    /// Columns `s + 1..=hap.len()` of one row; `cur[s]` must hold column `s`.
    /// `x` is each lane's read base, `p` the row parameters (m2m, m2i, m2d,
    /// i2m, i2i, p_match, p_mismatch) per lane.
    #[inline]
    pub fn row(
        x: &[u8; LANES],
        p: &[V; 7],
        hap: &[u8],
        s: usize,
        prev: (&[V], &[V], &[V]),
        cur: (&mut [V], &mut [V], &mut [V]),
    ) {
        let n = hap.len();
        assert!(
            s <= n
                && prev.0.len() > n
                && prev.1.len() > n
                && prev.2.len() > n
                && cur.0.len() > n
                && cur.1.len() > n
                && cur.2.len() > n
        );
        // SAFETY: NEON is part of the aarch64 baseline; lengths checked above.
        unsafe { row_neon(x, p, hap, s, prev, cur) }
    }

    #[target_feature(enable = "neon")]
    fn row_neon(
        x: &[u8; LANES],
        p: &[V; 7],
        hap: &[u8],
        s: usize,
        (mp, ip, dp): (&[V], &[V], &[V]),
        (mc, ic, dc): (&mut [V], &mut [V], &mut [V]),
    ) {
        const H: usize = LANES / 2;
        let ld = |v: &V, k: usize| unsafe { vld1q_f64(v.as_ptr().add(2 * k)) };
        let st = |d: &mut V, k: usize, a: float64x2_t| unsafe {
            vst1q_f64(d.as_mut_ptr().add(2 * k), a)
        };
        let params: [[float64x2_t; H]; 7] =
            std::array::from_fn(|q| std::array::from_fn(|k| ld(&p[q], k)));
        let [m2m, m2i, m2d, i2m, i2i, pm, pmm] = params;
        let xs: [u64; LANES] = x.map(|b| b as u64);
        let xv: [uint64x2_t; H] =
            std::array::from_fn(|k| unsafe { vld1q_u64(xs.as_ptr().add(2 * k)) });
        let nv = vdupq_n_u64(b'N' as u64);
        // Lanes whose read base is N match every haplotype base.
        let x_is_n: [uint64x2_t; H] = std::array::from_fn(|k| vceqq_u64(xv[k], nv));
        let mut mcl: [float64x2_t; H] = std::array::from_fn(|k| ld(&mc[s], k));
        let mut dcl: [float64x2_t; H] = std::array::from_fn(|k| ld(&dc[s], k));
        for j in s + 1..=hap.len() {
            let y = hap[j - 1];
            let yv = vdupq_n_u64(y as u64);
            let y_is_n = vdupq_n_u64(if y == b'N' { u64::MAX } else { 0 });
            for k in 0..H {
                let hit = vorrq_u64(vorrq_u64(vceqq_u64(xv[k], yv), x_is_n[k]), y_is_n);
                let prior = vbslq_f64(hit, pm[k], pmm[k]);
                // prior * (mp * m2m + ip * i2m + dp * i2m), left to right.
                let mv = vmulq_f64(
                    prior,
                    vaddq_f64(
                        vaddq_f64(
                            vmulq_f64(ld(&mp[j - 1], k), m2m[k]),
                            vmulq_f64(ld(&ip[j - 1], k), i2m[k]),
                        ),
                        vmulq_f64(ld(&dp[j - 1], k), i2m[k]),
                    ),
                );
                let iv = vaddq_f64(
                    vmulq_f64(ld(&mp[j], k), m2i[k]),
                    vmulq_f64(ld(&ip[j], k), i2i[k]),
                );
                // d2d = i2i; uses this row's previous column.
                let dv = vaddq_f64(vmulq_f64(mcl[k], m2d[k]), vmulq_f64(dcl[k], i2i[k]));
                st(&mut mc[j], k, mv);
                st(&mut ic[j], k, iv);
                st(&mut dc[j], k, dv);
                mcl[k] = mv;
                dcl[k] = dv;
            }
        }
    }
}
