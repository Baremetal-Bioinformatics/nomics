//! Exact fast path for `SmithWatermanJavaAligner` (oracle:
//! [`super::cigar::smith_waterman_reference`]), plus the `NOMICS_DUMP_SW`
//! dump.
//!
//! Why it is exact by construction:
//!
//! - **Order.** Every cell's scalar computation (Java's comparison order and
//!   tie rules, see [`nomics_simd::sw_diag`]) reads only `(i-1, j-1)`,
//!   `(i-1, j)` and `(i, j-1)`: `best_v[j]` at `(i, j)` is a function of
//!   `best_v`/score at `(i-1, j)`, `best_h` of `best_h`/score at `(i, j-1)`.
//!   Evaluating whole anti-diagonals therefore performs exactly the same
//!   per-cell operations on the same operands as the row loop.
//! - **Back pointers.** One byte per cell instead of Java's `i32`: the
//!   chosen direction plus whether `best_v` / `best_h` *opened* at that
//!   cell. Java's `size_v[j]` at `(i, j)` (1 on open, +1 on extend) equals
//!   the number of rows from `i` up column `j` to the nearest open cell,
//!   inclusive, and `size_h` likewise along the row; the backtrack walks
//!   those bits. A walk always ends: in row 1, `best_v = LOW + extend`
//!   (`LOW = i32::MIN / 2 = -2^30`) while the open candidate is a row-0 or
//!   clamped score, `> -2^28` under [`in_range`], so it opens; column 1
//!   likewise for `best_h`.
//! - **Overflow.** [`in_range`] bounds every score by
//!   `(rows + cols) * max|weight| < 2^28`, so no addition wraps in either
//!   implementation; calls outside it (and empty inputs) go to the oracle.
//! - **Final scans / backtrack.** Same code as the oracle, over the stored
//!   last column and last row.

use super::cigar::{D, Elem, I, M, S, Strategy, java_last_index_of, smith_waterman_reference};
use nomics_simd::sw_diag::{self, DIR_DEL, DIR_INS, DIR_MASK, E_OPEN, F_OPEN, Weights};
use std::cell::RefCell;

const LOW: i32 = i32::MIN / 2;

/// Diagonals are evaluated in whole blocks of this many cells.
const PAD: usize = 8;

/// The fast paths keep every intermediate far from `i32` overflow (and from
/// `LOW`): `|score| <= (rows + cols) * max|weight|` must stay below 2^28.
fn in_range(reference: &[u8], alternate: &[u8], w: [i32; 4]) -> bool {
    let max_w = w.iter().map(|x| i64::from(*x).abs()).max().unwrap_or(0);
    let span = (reference.len() + alternate.len() + 2) as i64;
    !reference.is_empty() && !alternate.is_empty() && span * (max_w + 1) < (1 << 28)
}

pub(crate) fn smith_waterman_fast(
    reference: &[u8],
    alternate: &[u8],
    w: [i32; 4],
    strategy: Strategy,
) -> (Vec<Elem>, i32) {
    wavefront(reference, alternate, w, strategy)
}

fn reference_call(r: &[u8], a: &[u8], w: [i32; 4], strategy: Strategy) -> (Vec<Elem>, i32) {
    smith_waterman_reference(r, a, w[0], w[1], w[2], w[3], strategy)
}

fn exact_hit(reference: &[u8], alternate: &[u8], strategy: Strategy) -> Option<(Vec<Elem>, i32)> {
    if matches!(strategy, Strategy::Softclip)
        && let Some(offset) = java_last_index_of(reference, alternate)
    {
        return Some((
            vec![Elem {
                len: alternate.len(),
                op: M,
            }],
            i32::try_from(offset).expect("alignment offset too large"),
        ));
    }
    None
}

/// Row 0 / column 0 scores.
fn edge(k: usize, w: [i32; 4], strategy: Strategy) -> i32 {
    match strategy {
        Strategy::LeadingIndel if k > 0 => w[2] + w[3] * (k as i32 - 1),
        _ => 0,
    }
}

/// Java's final scans and backtrack. `last_col[i]` = score at
/// `(i, alt.len())` and `last_row[j]` = score at `(ref.len(), j)` for
/// `i, j >= 1`; `code(i, j)` is the byte back code of cell `(i, j)`.
fn finish(
    reference: &[u8],
    alternate: &[u8],
    strategy: Strategy,
    last_col: &[i32],
    last_row: &[i32],
    code: impl Fn(usize, usize) -> u8,
) -> (Vec<Elem>, i32) {
    let rows = reference.len() + 1;
    let cols = alternate.len() + 1;
    let mut p1 = 0_usize;
    let mut max_score = i32::MIN;
    for (i, &score) in last_col.iter().enumerate().take(rows).skip(1) {
        if score >= max_score {
            p1 = i;
            max_score = score;
        }
    }
    let mut p2 = alternate.len();
    let mut trailing = 0_usize;
    if matches!(strategy, Strategy::Softclip) {
        for (j, &score) in last_row.iter().enumerate().take(cols).skip(1) {
            if score > max_score
                || (score == max_score && reference.len().abs_diff(j) < p1.abs_diff(p2))
            {
                p1 = reference.len();
                p2 = j;
                max_score = score;
                trailing = alternate.len() - j;
            }
        }
    }

    let mut reverse = Vec::<Elem>::new();
    if trailing > 0 && matches!(strategy, Strategy::Softclip) {
        reverse.push(Elem {
            len: trailing,
            op: S,
        });
    }
    let mut state = M;
    let mut segment = 0_usize;
    while p1 > 0 && p2 > 0 {
        let c = code(p1, p2);
        let (new_state, step) = match c & DIR_MASK {
            DIR_DEL => {
                let mut r = p1;
                while code(r, p2) & E_OPEN == 0 {
                    r -= 1;
                }
                (D, p1 - r + 1)
            }
            DIR_INS => {
                let mut q = p2;
                while code(p1, q) & F_OPEN == 0 {
                    q -= 1;
                }
                (I, p2 - q + 1)
            }
            _ => (M, 1),
        };
        match new_state {
            M => {
                p1 -= 1;
                p2 -= 1;
            }
            I => p2 -= step,
            _ => p1 -= step,
        }
        if new_state == state {
            segment += step;
        } else {
            if segment != 0 {
                reverse.push(Elem {
                    len: segment,
                    op: state,
                });
            }
            segment = step;
            state = new_state;
        }
    }
    let offset = if matches!(strategy, Strategy::Softclip) {
        reverse.push(Elem {
            len: segment,
            op: state,
        });
        if p2 > 0 {
            reverse.push(Elem { len: p2, op: S });
        }
        i32::try_from(p1).expect("offset too large")
    } else {
        reverse.push(Elem {
            len: segment,
            op: state,
        });
        if p1 > 0 {
            reverse.push(Elem { len: p1, op: D });
        } else if p2 > 0 {
            reverse.push(Elem { len: p2, op: I });
        }
        0
    };
    reverse.reverse();
    (reverse, offset)
}

/// Grow-only buffer: contents are garbage between calls; every value read
/// is written first in the same call.
fn grow<T: Copy + Default>(v: &mut Vec<T>, n: usize) -> &mut [T] {
    if v.len() < n {
        v.resize(n, T::default());
    }
    &mut v[..n]
}

thread_local! {
    static WAVE: RefCell<WaveScratch> = RefCell::default();
}

#[derive(Default)]
struct WaveScratch {
    back: Vec<u8>,
    off: Vec<usize>,
    refs: Vec<i32>,
    alts: Vec<i32>,
    h: [Vec<i32>; 3],
    e: [Vec<i32>; 2],
    f: [Vec<i32>; 2],
    last_col: Vec<i32>,
    last_row: Vec<i32>,
}

/// Anti-diagonal fill ([`nomics_simd::sw_diag`]).
pub(crate) fn wavefront(
    reference: &[u8],
    alternate: &[u8],
    w: [i32; 4],
    strategy: Strategy,
) -> (Vec<Elem>, i32) {
    wavefront_with(reference, alternate, w, strategy, sw_diag::diag)
}

/// [`wavefront`] with the portable diagonal kernel.
pub(crate) fn wavefront_portable(
    reference: &[u8],
    alternate: &[u8],
    w: [i32; 4],
    strategy: Strategy,
) -> (Vec<Elem>, i32) {
    wavefront_with(reference, alternate, w, strategy, sw_diag::diag_portable)
}

type DiagFn = fn(Weights, &sw_diag::Prev, &mut sw_diag::Cur);

fn wavefront_with(
    reference: &[u8],
    alternate: &[u8],
    w: [i32; 4],
    strategy: Strategy,
    kernel: DiagFn,
) -> (Vec<Elem>, i32) {
    if let Some(hit) = exact_hit(reference, alternate, strategy) {
        return hit;
    }
    if !in_range(reference, alternate, w) {
        return reference_call(reference, alternate, w, strategy);
    }
    WAVE.with_borrow_mut(|s| wave_fill(reference, alternate, w, strategy, s, kernel))
}

fn wave_fill(
    reference: &[u8],
    alternate: &[u8],
    w: [i32; 4],
    strategy: Strategy,
    s: &mut WaveScratch,
    kernel: DiagFn,
) -> (Vec<Elem>, i32) {
    let weights = Weights {
        match_value: w[0],
        mismatch: w[1],
        gap_open: w[2],
        gap_extend: w[3],
    };
    let rows = reference.len() + 1;
    let cols = alternate.len() + 1;
    // Cell (i, j) lies on diagonal d = i + j, 2 <= d <= rows + cols - 2,
    // rows lo(d)..=hi(d); its back code is at off[d] + i - lo(d).
    let lo = |d: usize| 1.max((d + 1).saturating_sub(cols));
    let hi = |d: usize| (rows - 1).min(d - 1);
    let last_d = rows + cols - 2;
    let back = grow(&mut s.back, (rows - 1) * (cols - 1) + PAD);
    let off = grow(&mut s.off, last_d + 1);
    let mut total = 0;
    for (d, o) in off.iter_mut().enumerate().skip(2) {
        *o = total;
        total += hi(d) + 1 - lo(d);
    }
    debug_assert_eq!(total + PAD, back.len());
    let refs = grow(&mut s.refs, rows - 1 + PAD);
    for (dst, &b) in refs.iter_mut().zip(reference) {
        *dst = i32::from(b);
    }
    // alts[t] = alternate[cols - 2 - t]: cell (i, d - i) reads
    // alternate[d - i - 1] = alts[cols - 1 - d + i].
    let alts = grow(&mut s.alts, cols - 1 + PAD);
    for (dst, &b) in alts.iter_mut().zip(alternate.iter().rev()) {
        *dst = i32::from(b);
    }
    let last_col = grow(&mut s.last_col, rows);
    let last_row = grow(&mut s.last_row, cols);
    let [h_a, h_b, h_c] = &mut s.h;
    let [e_a, e_b] = &mut s.e;
    let [f_a, f_b] = &mut s.f;
    // h2 = diagonal d - 2, h1 = d - 1, h0 = d (being written).
    let len = rows + PAD + 1;
    let (mut h2, mut h1, mut h0) = (grow(h_a, len), grow(h_b, len), grow(h_c, len));
    let (mut e1, mut e0) = (grow(e_a, len), grow(e_b, len));
    let (mut f1, mut f0) = (grow(f_a, len), grow(f_b, len));
    // Boundary cells of a finished diagonal t: (0, t) and (t, 0).
    let boundaries = |t: usize, h: &mut [i32], e: &mut [i32], f: &mut [i32]| {
        if t < cols {
            h[0] = edge(t, w, strategy);
            e[0] = LOW;
        }
        if t < rows {
            h[t] = edge(t, w, strategy);
            f[t] = LOW;
        }
    };
    // Diagonal 0 is (0, 0); diagonal 1 is (0, 1) and (1, 0).
    boundaries(0, h2, e0, f0);
    boundaries(1, h1, e1, f1);
    for d in 2..=last_d {
        let (l, u) = (lo(d), hi(d));
        // Cells l..=u, rounded up to whole kernel blocks: the extra cells
        // compute garbage into slots that are rewritten before being read
        // (h0/e0/f0 beyond u are only read after this diagonal's boundaries
        // are set or a later diagonal rewrites them; back codes beyond n
        // belong to diagonal d + 1, written next).
        let n = (u + 1 - l).next_multiple_of(PAD);
        let a0 = cols - 1 + l - d;
        let prev = sw_diag::Prev {
            refs: &refs[l - 1..l - 1 + n],
            alts: &alts[a0..a0 + n],
            h2: &h2[l - 1..l - 1 + n],
            h1: &h1[l - 1..l + n],
            e1: &e1[l - 1..l - 1 + n],
            f1: &f1[l - 1..l + n],
        };
        let mut cur = sw_diag::Cur {
            h: &mut h0[l..l + n],
            e: &mut e0[l..l + n],
            f: &mut f0[l..l + n],
            back: &mut back[off[d]..off[d] + n],
        };
        kernel(weights, &prev, &mut cur);
        if d >= cols && d + 1 - cols <= u {
            last_col[d + 1 - cols] = h0[d + 1 - cols];
        }
        if u == rows - 1 {
            last_row[d + 1 - rows] = h0[rows - 1];
        }
        boundaries(d, h0, e0, f0);
        std::mem::swap(&mut h2, &mut h1);
        std::mem::swap(&mut h1, &mut h0);
        std::mem::swap(&mut e1, &mut e0);
        std::mem::swap(&mut f1, &mut f0);
    }
    let back = &*back;
    let off = &*off;
    finish(
        reference,
        alternate,
        strategy,
        last_col,
        last_row,
        |i, j| {
            let d = i + j;
            back[off[d] + i - lo(d)]
        },
    )
}

/// `NOMICS_DUMP_SW=<path>`: record every [`super::cigar::smith_waterman`]
/// call (at most `NOMICS_DUMP_SW_MAX`, default 2,000,000). Little-endian
/// records: `u8 strategy (0 Softclip, 1 LeadingIndel), i32 match, mismatch,
/// gap_open, gap_extend, u32 len, reference, u32 len, alternate,
/// u32 n_elems, (u32 len, u8 op)*, i32 offset`.
pub(crate) fn dump_sink() -> Option<&'static std::sync::Mutex<(std::fs::File, u64)>> {
    static SINK: std::sync::OnceLock<Option<std::sync::Mutex<(std::fs::File, u64)>>> =
        std::sync::OnceLock::new();
    SINK.get_or_init(|| {
        let path = std::env::var_os("NOMICS_DUMP_SW")?;
        let max = std::env::var("NOMICS_DUMP_SW_MAX")
            .ok()
            .map_or(2_000_000, |v| v.parse().expect("NOMICS_DUMP_SW_MAX"));
        let f = std::fs::File::create(path).expect("create Smith-Waterman dump");
        Some(std::sync::Mutex::new((f, max)))
    })
    .as_ref()
}

pub(crate) fn dump_call(
    sink: &std::sync::Mutex<(std::fs::File, u64)>,
    reference: &[u8],
    alternate: &[u8],
    w: [i32; 4],
    strategy: Strategy,
    (cigar, offset): &(Vec<Elem>, i32),
) {
    use std::io::Write as _;
    let mut b = Vec::with_capacity(reference.len() + alternate.len() + 64);
    let u32le = |b: &mut Vec<u8>, x: usize| b.extend_from_slice(&(x as u32).to_le_bytes());
    b.push(u8::from(matches!(strategy, Strategy::LeadingIndel)));
    for x in w {
        b.extend_from_slice(&x.to_le_bytes());
    }
    u32le(&mut b, reference.len());
    b.extend_from_slice(reference);
    u32le(&mut b, alternate.len());
    b.extend_from_slice(alternate);
    u32le(&mut b, cigar.len());
    for e in cigar {
        u32le(&mut b, e.len);
        b.push(e.op);
    }
    b.extend_from_slice(&offset.to_le_bytes());
    let mut guard = sink.lock().unwrap();
    if guard.1 > 0 {
        guard.1 -= 1;
        guard.0.write_all(&b).expect("write Smith-Waterman dump");
    }
}

/// Benchmark entry points (`examples/sw_bench.rs`); not a stable API.
#[doc(hidden)]
pub mod bench {
    use super::super::cigar::{Elem, Strategy, smith_waterman_reference};

    #[derive(Clone, Copy, Debug)]
    pub enum Variant {
        Reference,
        WavefrontPortable,
        Wavefront,
    }

    /// One call; returns `(len, op)` elements and the offset.
    pub fn run(
        variant: Variant,
        reference: &[u8],
        alternate: &[u8],
        w: [i32; 4],
        leading_indel: bool,
    ) -> (Vec<(u32, u8)>, i32) {
        let strategy = if leading_indel {
            Strategy::LeadingIndel
        } else {
            Strategy::Softclip
        };
        let (cigar, offset) = match variant {
            Variant::Reference => {
                smith_waterman_reference(reference, alternate, w[0], w[1], w[2], w[3], strategy)
            }
            Variant::WavefrontPortable => {
                super::wavefront_portable(reference, alternate, w, strategy)
            }
            Variant::Wavefront => super::wavefront(reference, alternate, w, strategy),
        };
        (
            cigar.iter().map(|e: &Elem| (e.len as u32, e.op)).collect(),
            offset,
        )
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

    /// Mutated copies of a random sequence (the shape HC aligns), plus
    /// unrelated pairs and tiny sizes, over all three weight sets.
    #[test]
    fn fast_paths_match_reference() {
        let weights = [
            [200, -150, -260, -11],
            [25, -50, -110, -6],
            [10, -15, -30, -5],
            [1, -1, -1, -1],
        ];
        let mut seed = 1_u64;
        for case in 0..6000 {
            let alphabet: &[u8] = if case % 3 == 0 { b"AC" } else { b"ACGTN" };
            let len = 1 + (lcg(&mut seed) % if case % 4 == 0 { 8 } else { 160 }) as usize;
            let base: Vec<u8> = (0..len)
                .map(|_| alphabet[(lcg(&mut seed) as usize) % alphabet.len()])
                .collect();
            let other: Vec<u8> = if case % 5 == 0 {
                let l = 1 + (lcg(&mut seed) % 160) as usize;
                (0..l)
                    .map(|_| alphabet[(lcg(&mut seed) as usize) % alphabet.len()])
                    .collect()
            } else {
                let mut v = Vec::new();
                for &b in &base {
                    match lcg(&mut seed) % 40 {
                        0 => {}
                        1 => {
                            v.push(b);
                            for _ in 0..lcg(&mut seed) % 12 {
                                v.push(alphabet[(lcg(&mut seed) as usize) % alphabet.len()]);
                            }
                        }
                        2 => v.push(alphabet[(lcg(&mut seed) as usize) % alphabet.len()]),
                        _ => v.push(b),
                    }
                }
                if v.is_empty() {
                    v.push(b'A');
                }
                v
            };
            let (r, a) = if case % 2 == 0 {
                (&base, &other)
            } else {
                (&other, &base)
            };
            for w in weights {
                for strategy in [Strategy::Softclip, Strategy::LeadingIndel] {
                    let want = reference_call(r, a, w, strategy);
                    assert_eq!(wavefront(r, a, w, strategy), want, "wave case {case}");
                    assert_eq!(
                        wavefront_portable(r, a, w, strategy),
                        want,
                        "wave portable case {case}"
                    );
                }
            }
        }
    }
}
