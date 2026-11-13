//! Lane-parallel experiment for BWA's `ksw_extend2`.
//!
//! A vector contains independent extension problems, not adjacent cells of
//! one problem.  Every live lane therefore follows the scalar recurrence and
//! its control flow independently.  Fixed-size arrays make the arithmetic
//! amenable to compiler auto-vectorization without requiring intrinsics.

use std::cell::RefCell;

use crate::ksw::Extend;

/// Independent extension problems evaluated in one vector pass.
pub const LANES: usize = 8;

type I32xL = [i32; LANES];
type I8xL = [i8; LANES];

/// Borrowed arguments for one `ksw_extend2` call.
#[derive(Clone, Copy)]
pub struct Extend2Call<'a> {
    pub query: &'a [u8],
    pub target: &'a [u8],
    pub m: i32,
    pub mat: &'a [i8],
    pub o_del: i32,
    pub e_del: i32,
    pub o_ins: i32,
    pub e_ins: i32,
    pub w: i32,
    pub end_bonus: i32,
    pub zdrop: i32,
    pub h0: i32,
}

/// Work accounting for a lane-parallel batch.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct LaneStats {
    /// Calls in the batch.
    pub calls: usize,
    /// Cells the corresponding scalar calls evaluate.
    pub useful_cells: u64,
    /// Vector iterations, each of which reserves [`LANES`] lane slots.
    pub vector_steps: u64,
}

impl LaneStats {
    pub fn lane_slots(self) -> u64 {
        self.vector_steps * LANES as u64
    }

    pub fn utilization(self) -> f64 {
        if self.vector_steps == 0 {
            1.0
        } else {
            self.useful_cells as f64 / self.lane_slots() as f64
        }
    }

    pub fn add(&mut self, other: Self) {
        self.calls += other.calls;
        self.useful_cells += other.useful_cells;
        self.vector_steps += other.vector_steps;
    }
}

#[derive(Default)]
struct LaneScratch {
    qp: Vec<I8xL>,
    h: Vec<I32xL>,
    e: Vec<I32xL>,
}

thread_local! {
    static SCRATCH: RefCell<LaneScratch> = RefCell::default();
}

/// Evaluate a batch of independent `ksw_extend2` calls in groups of
/// [`LANES`]. Results have the same order as `calls`.
///
/// The final partial group leaves its unused lanes masked. Grouping calls of
/// similar lengths generally improves lane utilization.
pub fn extend2_batch(calls: &[Extend2Call<'_>], out: &mut [Extend]) -> LaneStats {
    assert_eq!(calls.len(), out.len());
    let mut stats = LaneStats {
        calls: calls.len(),
        ..LaneStats::default()
    };
    SCRATCH.with_borrow_mut(|scratch| {
        for (group, group_out) in calls.chunks(LANES).zip(out.chunks_mut(LANES)) {
            let group_stats = extend2_group(scratch, group, group_out);
            stats.useful_cells += group_stats.useful_cells;
            stats.vector_steps += group_stats.vector_steps;
        }
    });
    stats
}

fn check_call(call: &Extend2Call<'_>) {
    assert!(call.query.len() <= i32::MAX as usize);
    assert!(call.target.len() <= i32::MAX as usize);
    assert!(call.m > 0);
    assert_eq!(call.mat.len(), call.m as usize * call.m as usize);
    assert!(
        call.query
            .iter()
            .chain(call.target)
            .all(|&base| base < call.m as u8)
    );
    assert!(call.h0 > 0, "BWA ksw_extend2 asserts h0 > 0");
}

fn extend2_group(
    scratch: &mut LaneScratch,
    calls: &[Extend2Call<'_>],
    out: &mut [Extend],
) -> LaneStats {
    debug_assert!(!calls.is_empty() && calls.len() <= LANES && calls.len() == out.len());
    for call in calls {
        check_call(call);
    }

    let mut stats = LaneStats::default();
    let mut running = [false; LANES];
    let mut qlen = [0i32; LANES];
    let mut tlen = [0i32; LANES];
    let mut m = [0usize; LANES];
    let mut oe_del = [0i32; LANES];
    let mut oe_ins = [0i32; LANES];
    let mut e_del = [0i32; LANES];
    let mut e_ins = [0i32; LANES];
    let mut w = [0i32; LANES];
    let mut max = [0i32; LANES];
    let mut max_i = [-1i32; LANES];
    let mut max_j = [-1i32; LANES];
    let mut max_ie = [-1i32; LANES];
    let mut gscore = [-1i32; LANES];
    let mut max_off = [0i32; LANES];
    let mut beg = [0i32; LANES];
    let mut end = [0i32; LANES];

    let max_qlen = calls.iter().map(|call| call.query.len()).max().unwrap();
    let max_m = calls.iter().map(|call| call.m as usize).max().unwrap();
    scratch.qp.clear();
    scratch.qp.resize(max_qlen * max_m, [0; LANES]);
    scratch.h.clear();
    scratch.h.resize(max_qlen + 1, [0; LANES]);
    scratch.e.clear();
    scratch.e.resize(max_qlen + 1, [0; LANES]);

    for (lane, call) in calls.iter().enumerate() {
        qlen[lane] = call.query.len() as i32;
        tlen[lane] = call.target.len() as i32;
        m[lane] = call.m as usize;
        oe_del[lane] = call.o_del.wrapping_add(call.e_del);
        oe_ins[lane] = call.o_ins.wrapping_add(call.e_ins);
        e_del[lane] = call.e_del;
        e_ins[lane] = call.e_ins;
        max[lane] = call.h0;
        end[lane] = qlen[lane];

        if qlen[lane] == 0 {
            out[lane] = Extend {
                score: call.h0,
                qle: 0,
                tle: 0,
                gtle: 0,
                gscore: -1,
                max_off: 0,
            };
            continue;
        }
        running[lane] = tlen[lane] != 0;

        for (base, row) in call.mat.chunks_exact(m[lane]).enumerate() {
            for (j, &query_base) in call.query.iter().enumerate() {
                scratch.qp[base * max_qlen + j][lane] = row[query_base as usize];
            }
        }

        scratch.h[0][lane] = call.h0;
        scratch.h[1][lane] = if call.h0 > oe_ins[lane] {
            call.h0.wrapping_sub(oe_ins[lane])
        } else {
            0
        };
        let mut j = 2usize;
        while j <= call.query.len() && scratch.h[j - 1][lane] > call.e_ins {
            scratch.h[j][lane] = scratch.h[j - 1][lane].wrapping_sub(call.e_ins);
            j += 1;
        }

        let mut max_score = 0i32;
        for &score in call.mat {
            max_score = max_score.max(score as i32);
        }
        let max_ins = (((qlen[lane]
            .wrapping_mul(max_score)
            .wrapping_add(call.end_bonus)
            .wrapping_sub(call.o_ins)) as f64
            / call.e_ins as f64
            + 1.0) as i32)
            .max(1);
        w[lane] = call.w.min(max_ins);
        let max_del = (((qlen[lane]
            .wrapping_mul(max_score)
            .wrapping_add(call.end_bonus)
            .wrapping_sub(call.o_del)) as f64
            / call.e_del as f64
            + 1.0) as i32)
            .max(1);
        w[lane] = w[lane].min(max_del);
    }

    let max_tlen = tlen.iter().copied().max().unwrap_or(0);
    for i in 0..max_tlen {
        let mut row_live = [false; LANES];
        let mut h1 = [0i32; LANES];
        let mut f = [0i32; LANES];
        let mut row_max = [0i32; LANES];
        let mut mj = [-1i32; LANES];
        let mut target_base = [0usize; LANES];
        let mut global_beg = i32::MAX;
        let mut global_end = 0i32;

        for lane in 0..LANES {
            if !running[lane] || i >= tlen[lane] {
                continue;
            }
            row_live[lane] = true;
            target_base[lane] = calls[lane].target[i as usize] as usize;
            if beg[lane] < i.wrapping_sub(w[lane]) {
                beg[lane] = i.wrapping_sub(w[lane]);
            }
            if end[lane] > i.wrapping_add(w[lane]).wrapping_add(1) {
                end[lane] = i.wrapping_add(w[lane]).wrapping_add(1);
            }
            if end[lane] > qlen[lane] {
                end[lane] = qlen[lane];
            }
            h1[lane] = if beg[lane] == 0 {
                let call = &calls[lane];
                let gap = call
                    .o_del
                    .wrapping_add(call.e_del.wrapping_mul(i.wrapping_add(1)));
                call.h0.wrapping_sub(gap).max(0)
            } else {
                0
            };
            global_beg = global_beg.min(beg[lane]);
            global_end = global_end.max(end[lane]);
            stats.useful_cells += (end[lane] - beg[lane]) as u64;
        }

        if global_beg == i32::MAX {
            break;
        }
        stats.vector_steps += (global_end - global_beg) as u64;
        for j in global_beg..global_end {
            let j_usize = j as usize;
            let diagonal = scratch.h[j_usize];
            let ev = scratch.e[j_usize];
            let mut score = [0i32; LANES];
            let mut cell_live = [false; LANES];
            for lane in 0..LANES {
                cell_live[lane] = row_live[lane] && j >= beg[lane] && j < end[lane];
                score[lane] = scratch.qp[target_base[lane] * max_qlen + j_usize][lane] as i32;
            }

            // Keep this as fixed-size array arithmetic: LLVM can lower the
            // independent integer lanes to SIMD while the masks preserve each
            // call's scalar band and early-termination state.
            let mut diagonal_score = [0i32; LANES];
            let mut hv = [0i32; LANES];
            let mut next_e = [0i32; LANES];
            let mut next_f = [0i32; LANES];
            for lane in 0..LANES {
                diagonal_score[lane] = if diagonal[lane] != 0 {
                    diagonal[lane].wrapping_add(score[lane])
                } else {
                    0
                };
                hv[lane] = diagonal_score[lane].max(ev[lane]).max(f[lane]);
                next_e[lane] = ev[lane]
                    .wrapping_sub(e_del[lane])
                    .max(diagonal_score[lane].wrapping_sub(oe_del[lane]).max(0));
                next_f[lane] = f[lane]
                    .wrapping_sub(e_ins[lane])
                    .max(diagonal_score[lane].wrapping_sub(oe_ins[lane]).max(0));
            }
            for lane in 0..LANES {
                let live = cell_live[lane];
                scratch.h[j_usize][lane] = if live { h1[lane] } else { diagonal[lane] };
                scratch.e[j_usize][lane] = if live { next_e[lane] } else { ev[lane] };
                f[lane] = if live { next_f[lane] } else { f[lane] };
                h1[lane] = if live { hv[lane] } else { h1[lane] };
                let next_mj = if row_max[lane] > hv[lane] {
                    mj[lane]
                } else {
                    j
                };
                mj[lane] = if live { next_mj } else { mj[lane] };
                row_max[lane] = if live {
                    row_max[lane].max(hv[lane])
                } else {
                    row_max[lane]
                };
            }
        }

        for lane in 0..LANES {
            if !row_live[lane] {
                continue;
            }
            let end_usize = end[lane] as usize;
            scratch.h[end_usize][lane] = h1[lane];
            scratch.e[end_usize][lane] = 0;
            if end[lane] == qlen[lane] {
                max_ie[lane] = if gscore[lane] > h1[lane] {
                    max_ie[lane]
                } else {
                    i
                };
                gscore[lane] = gscore[lane].max(h1[lane]);
            }
            if row_max[lane] == 0 {
                running[lane] = false;
                continue;
            }
            if row_max[lane] > max[lane] {
                max[lane] = row_max[lane];
                max_i[lane] = i;
                max_j[lane] = mj[lane];
                max_off[lane] = max_off[lane].max((mj[lane] - i).abs());
            } else if calls[lane].zdrop > 0 {
                let di = i - max_i[lane];
                let dj = mj[lane] - max_j[lane];
                let dropped = if di > dj {
                    max[lane]
                        .wrapping_sub(row_max[lane])
                        .wrapping_sub((di - dj).wrapping_mul(calls[lane].e_del))
                        > calls[lane].zdrop
                } else {
                    max[lane]
                        .wrapping_sub(row_max[lane])
                        .wrapping_sub((dj - di).wrapping_mul(calls[lane].e_ins))
                        > calls[lane].zdrop
                };
                if dropped {
                    running[lane] = false;
                    continue;
                }
            }

            let begin = beg[lane] as usize;
            let leading_zeroes = (begin..end_usize)
                .position(|j| scratch.h[j][lane] != 0 || scratch.e[j][lane] != 0)
                .unwrap_or(end_usize - begin);
            beg[lane] += leading_zeroes as i32;
            let begin = beg[lane] as usize;
            let last_nonzero = (begin..=end_usize)
                .rev()
                .find(|&j| scratch.h[j][lane] != 0 || scratch.e[j][lane] != 0)
                .map_or(beg[lane] - 1, |j| j as i32);
            end[lane] = (last_nonzero + 2).min(qlen[lane]);
        }
    }

    for lane in 0..calls.len() {
        if qlen[lane] != 0 {
            out[lane] = Extend {
                score: max[lane],
                qle: max_j[lane] + 1,
                tle: max_i[lane] + 1,
                gtle: max_ie[lane] + 1,
                gscore: gscore[lane],
                max_off: max_off[lane],
            };
        }
    }
    stats
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ksw;

    const MAT: [i8; 25] = [
        1, -4, -4, -4, -1, -4, 1, -4, -4, -1, -4, -4, 1, -4, -1, -4, -4, -4, 1, -1, -1, -1, -1, -1,
        -1,
    ];

    #[test]
    fn lanes_match_scalar_with_divergent_calls() {
        let queries: Vec<Vec<u8>> = (0..19)
            .map(|n| (0..n).map(|j| ((j * 3 + n) % 5) as u8).collect())
            .collect();
        let targets: Vec<Vec<u8>> = (0..19)
            .map(|n| (0..(n * 2 + 3)).map(|j| ((j * 7 + n) % 5) as u8).collect())
            .collect();
        let calls: Vec<_> = (0..19)
            .map(|n| Extend2Call {
                query: &queries[n],
                target: &targets[18 - n],
                m: 5,
                mat: &MAT,
                o_del: 6,
                e_del: 1,
                o_ins: 5 + (n % 3) as i32,
                e_ins: 1,
                w: 2 + (n % 9) as i32,
                end_bonus: 5,
                zdrop: if n % 4 == 0 { 0 } else { 3 + n as i32 },
                h0: 1 + n as i32,
            })
            .collect();
        let expected: Vec<_> = calls
            .iter()
            .map(|c| {
                ksw::extend2(
                    c.query.len() as i32,
                    c.query,
                    c.target.len() as i32,
                    c.target,
                    c.m,
                    c.mat,
                    c.o_del,
                    c.e_del,
                    c.o_ins,
                    c.e_ins,
                    c.w,
                    c.end_bonus,
                    c.zdrop,
                    c.h0,
                )
            })
            .collect();
        let mut actual = vec![expected[0]; calls.len()];
        let stats = extend2_batch(&calls, &mut actual);
        assert_eq!(actual, expected);
        assert_eq!(stats.calls, calls.len());
        assert!(stats.useful_cells <= stats.lane_slots());
    }
}
