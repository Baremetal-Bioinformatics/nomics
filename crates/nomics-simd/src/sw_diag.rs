//! One anti-diagonal of GATK's `SmithWatermanJavaAligner` matrix fill.
//!
//! Cell `k` of the diagonal reads its up neighbour at index `k` and its left
//! neighbour at index `k + 1` of the previous diagonal (`h1`, `e1`, `f1`),
//! and its diagonal neighbour at index `k` of the diagonal before that
//! (`h2`). Per cell, in Java's order:
//!
//! ```text
//! diag = h2 + (ref == alt ? match : mismatch)
//! e    = h1[up] + open  >  e1[up] + extend   ? open (bit 2) : extend
//! f    = h1[left] + open > f1[left] + extend ? open (bit 3) : extend
//! diag >= e && diag >= f  -> h = diag, dir 0
//! else f >= e             -> h = f,    dir 1 (insertion, Java's -size_h)
//! else                    -> h = e,    dir 2 (deletion,  Java's  size_v)
//! h = max(h, -100_000_000)
//! ```
//!
//! `back[k]` is `dir | e_open << 2 | f_open << 3`. Additions wrap in both
//! versions; the caller keeps every value far from `i32` overflow.

/// GATK's `MATRIX_MIN_CUTOFF`.
pub const MATRIX_MIN_CUTOFF: i32 = -100_000_000;

/// Back-pointer bits.
pub const DIR_MASK: u8 = 3;
pub const DIR_INS: u8 = 1;
pub const DIR_DEL: u8 = 2;
pub const E_OPEN: u8 = 4;
pub const F_OPEN: u8 = 8;

#[derive(Clone, Copy, Debug)]
pub struct Weights {
    pub match_value: i32,
    pub mismatch: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
}

/// Inputs of one diagonal of `n = back.len()` cells: `refs[k]`/`alts[k]`
/// are the two bases (widened) of cell `k`; `h2` and `e1` need `n` values,
/// `h1` and `f1` need `n + 1`.
pub struct Prev<'a> {
    pub refs: &'a [i32],
    pub alts: &'a [i32],
    pub h2: &'a [i32],
    pub h1: &'a [i32],
    pub e1: &'a [i32],
    pub f1: &'a [i32],
}

/// Outputs of one diagonal, each `n` long.
pub struct Cur<'a> {
    pub h: &'a mut [i32],
    pub e: &'a mut [i32],
    pub f: &'a mut [i32],
    pub back: &'a mut [u8],
}

fn check(p: &Prev, c: &Cur) -> usize {
    let n = c.back.len();
    assert!(
        p.refs.len() >= n
            && p.alts.len() >= n
            && p.h2.len() >= n
            && p.e1.len() >= n
            && p.h1.len() > n
            && p.f1.len() > n
            && c.h.len() >= n
            && c.e.len() >= n
            && c.f.len() >= n
    );
    n
}

#[inline(always)]
fn cell(w: Weights, p: &Prev, c: &mut Cur, k: usize) {
    let s = if p.refs[k] == p.alts[k] {
        w.match_value
    } else {
        w.mismatch
    };
    let diag = p.h2[k].wrapping_add(s);
    let e_open = p.h1[k].wrapping_add(w.gap_open);
    let e_ext = p.e1[k].wrapping_add(w.gap_extend);
    let eo = e_open > e_ext;
    let e = if eo { e_open } else { e_ext };
    let f_open = p.h1[k + 1].wrapping_add(w.gap_open);
    let f_ext = p.f1[k + 1].wrapping_add(w.gap_extend);
    let fo = f_open > f_ext;
    let f = if fo { f_open } else { f_ext };
    let (h, dir) = if diag >= e && diag >= f {
        (diag, 0)
    } else if f >= e {
        (f, DIR_INS)
    } else {
        (e, DIR_DEL)
    };
    c.h[k] = h.max(MATRIX_MIN_CUTOFF);
    c.e[k] = e;
    c.f[k] = f;
    c.back[k] = dir | (u8::from(eo) << 2) | (u8::from(fo) << 3);
}

/// Portable version; [`diag`] must match it exactly.
pub fn diag_portable(w: Weights, p: &Prev, c: &mut Cur) {
    let n = check(p, c);
    for k in 0..n {
        cell(w, p, c, k);
    }
}

/// One diagonal: NEON on aarch64, [`diag_portable`] elsewhere.
#[inline]
pub fn diag(w: Weights, p: &Prev, c: &mut Cur) {
    #[cfg(target_arch = "aarch64")]
    {
        let n = check(p, c);
        // SAFETY: NEON is part of the aarch64 baseline; lengths checked.
        unsafe { neon::diag(w, p, c, n) }
    }
    #[cfg(not(target_arch = "aarch64"))]
    diag_portable(w, p, c)
}

#[cfg(target_arch = "aarch64")]
mod neon {
    use super::{Cur, MATRIX_MIN_CUTOFF, Prev, Weights, cell};
    use std::arch::aarch64::*;

    #[target_feature(enable = "neon")]
    pub(super) fn diag(w: Weights, p: &Prev, c: &mut Cur, n: usize) {
        // Every access below is to cells k..k + 4 (or k + 1..k + 5 for the
        // left neighbours) with k + 8 <= n, inside the lengths `check`
        // asserted: whole `[i32; 4]` loads and stores.
        let ld = |s: &[i32], k: usize| -> int32x4_t {
            debug_assert!(k + 4 <= s.len());
            unsafe { vld1q_s32(s.as_ptr().add(k)) }
        };
        let st = |s: &mut [i32], k: usize, v: int32x4_t| {
            debug_assert!(k + 4 <= s.len());
            unsafe { vst1q_s32(s.as_mut_ptr().add(k), v) }
        };
        let mv = vdupq_n_s32(w.match_value);
        let mmv = vdupq_n_s32(w.mismatch);
        let go = vdupq_n_s32(w.gap_open);
        let ge = vdupq_n_s32(w.gap_extend);
        let floor = vdupq_n_s32(MATRIX_MIN_CUTOFF);
        let one = vdupq_n_u32(1);
        let two = vdupq_n_u32(2);
        let four = vdupq_n_u32(4);
        let eight = vdupq_n_u32(8);
        // One group of four cells; returns the four back codes as u32 lanes.
        let quad = |c: &mut Cur, k: usize| -> uint32x4_t {
            let s = vbslq_s32(vceqq_s32(ld(p.refs, k), ld(p.alts, k)), mv, mmv);
            let diag = vaddq_s32(ld(p.h2, k), s);
            let h_up = ld(p.h1, k);
            let h_left = ld(p.h1, k + 1);
            let e_open = vaddq_s32(h_up, go);
            let e_ext = vaddq_s32(ld(p.e1, k), ge);
            let eo = vcgtq_s32(e_open, e_ext);
            let e = vbslq_s32(eo, e_open, e_ext);
            let f_open = vaddq_s32(h_left, go);
            let f_ext = vaddq_s32(ld(p.f1, k + 1), ge);
            let fo = vcgtq_s32(f_open, f_ext);
            let f = vbslq_s32(fo, f_open, f_ext);
            let take_diag = vandq_u32(vcgeq_s32(diag, e), vcgeq_s32(diag, f));
            let take_f = vcgeq_s32(f, e);
            let gap = vbslq_s32(take_f, f, e);
            let h = vbslq_s32(take_diag, diag, gap);
            st(c.h, k, vmaxq_s32(h, floor));
            st(c.e, k, e);
            st(c.f, k, f);
            let dir = vbicq_u32(vbslq_u32(take_f, one, two), take_diag);
            vorrq_u32(vorrq_u32(dir, vandq_u32(eo, four)), vandq_u32(fo, eight))
        };
        let mut k = 0;
        while k + 8 <= n {
            let lo = quad(c, k);
            let hi = quad(c, k + 4);
            let codes = vmovn_u16(vcombine_u16(vmovn_u32(lo), vmovn_u32(hi)));
            unsafe { vst1_u8(c.back.as_mut_ptr().add(k), codes) };
            k += 8;
        }
        while k < n {
            cell(w, p, c, k);
            k += 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lcg(state: &mut u64) -> u32 {
        *state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (*state >> 33) as u32
    }

    #[test]
    fn neon_matches_portable() {
        let mut seed = 7_u64;
        let w = Weights {
            match_value: 200,
            mismatch: -150,
            gap_open: -260,
            gap_extend: -11,
        };
        for n in 0..70 {
            for _ in 0..200 {
                let mut v = |m: usize, spread: u32| -> Vec<i32> {
                    (0..m)
                        .map(|_| {
                            let x = lcg(&mut seed);
                            // Mix narrow values (many ties) and wide ones.
                            if x & 1 == 0 {
                                ((x >> 1) % 8) as i32 * 11 - 40
                            } else {
                                ((x >> 1) % spread) as i32 - (spread / 2) as i32
                            }
                        })
                        .collect()
                };
                let refs: Vec<i32> = v(n, 4);
                let alts: Vec<i32> = v(n, 4);
                let h2 = v(n, 200_000_000);
                let h1 = v(n + 1, 2000);
                let e1 = v(n, 2000);
                let f1 = v(n + 1, 2000);
                let p = Prev {
                    refs: &refs,
                    alts: &alts,
                    h2: &h2,
                    h1: &h1,
                    e1: &e1,
                    f1: &f1,
                };
                let mut out = [vec![0; n], vec![0; n], vec![0; n]];
                let mut out2 = out.clone();
                let mut b = vec![0_u8; n];
                let mut b2 = b.clone();
                {
                    let [h, e, f] = &mut out;
                    diag(
                        w,
                        &p,
                        &mut Cur {
                            h,
                            e,
                            f,
                            back: &mut b,
                        },
                    );
                    let [h, e, f] = &mut out2;
                    diag_portable(
                        w,
                        &p,
                        &mut Cur {
                            h,
                            e,
                            f,
                            back: &mut b2,
                        },
                    );
                }
                assert_eq!(out, out2);
                assert_eq!(b, b2);
            }
        }
    }
}
