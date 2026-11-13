//! Scalar, literal ports of BWA 0.7.19's `ksw_extend2` and `ksw_global2`.
//!
//! The order of the comparisons here is part of the interface: it determines
//! both the chosen end point and the global-alignment CIGAR on ties.

use std::cell::RefCell;

/// Result of [`extend2`]. End positions are lengths (and thus zero-based
/// positions in BWA plus one).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Extend {
    pub score: i32,
    pub qle: i32,
    pub tle: i32,
    pub gtle: i32,
    pub gscore: i32,
    pub max_off: i32,
}

#[derive(Clone, Copy, Default)]
struct Eh {
    h: i32,
    e: i32,
}

#[derive(Default)]
struct ExtendScratch {
    qp: Vec<i8>,
    h: Vec<i32>,
    e: Vec<i32>,
    eh: Vec<Eh>,
}

thread_local! {
    static EXTEND_SCRATCH: RefCell<ExtendScratch> = const {
        RefCell::new(ExtendScratch {
            qp: Vec::new(),
            h: Vec::new(),
            e: Vec::new(),
            eh: Vec::new(),
        })
    };
}

const MINUS_INF: i32 = -0x4000_0000;

fn checked_lengths(qlen: i32, query: &[u8], tlen: i32, target: &[u8], m: i32, mat: &[i8]) {
    assert!(qlen >= 0 && tlen >= 0 && m > 0);
    assert_eq!(query.len(), qlen as usize);
    assert_eq!(target.len(), tlen as usize);
    assert_eq!(mat.len(), (m as usize) * (m as usize));
    assert!(query.iter().chain(target).all(|&base| base < m as u8));
}

/// BWA's `ksw_extend2`.
///
/// `query` and `target` contain alphabet codes in `0..m`. As in BWA, the
/// first base of a gap costs `open + extend`.
#[allow(clippy::too_many_arguments)]
pub fn extend2(
    qlen: i32,
    query: &[u8],
    tlen: i32,
    target: &[u8],
    m: i32,
    mat: &[i8],
    o_del: i32,
    e_del: i32,
    o_ins: i32,
    e_ins: i32,
    w: i32,
    end_bonus: i32,
    zdrop: i32,
    h0: i32,
) -> Extend {
    let r = extend2_impl(
        qlen, query, tlen, target, m, mat, o_del, e_del, o_ins, e_ins, w, end_bonus, zdrop, h0,
    );
    if let Some(sink) = dump::sink(dump::EXTEND2) {
        let mut b = Vec::with_capacity(64 + query.len() + target.len());
        dump::seq(&mut b, &query[..qlen as usize]);
        dump::seq(&mut b, &target[..tlen as usize]);
        dump::ints(&mut b, &[m]);
        b.extend(mat[..(m * m) as usize].iter().map(|&x| x as u8));
        dump::ints(
            &mut b,
            &[o_del, e_del, o_ins, e_ins, w, end_bonus, zdrop, h0],
        );
        dump::ints(
            &mut b,
            &[r.score, r.qle, r.tle, r.gtle, r.gscore, r.max_off],
        );
        dump::write(sink, &b);
    }
    r
}

/// The literal row-by-row port of `ksw_extend2`, without any fast path. This
/// is the spec that [`extend2`] must match exactly; kept for tests and
/// benchmarks.
#[allow(clippy::too_many_arguments)]
pub fn extend2_reference(
    qlen: i32,
    query: &[u8],
    tlen: i32,
    target: &[u8],
    m: i32,
    mat: &[i8],
    o_del: i32,
    e_del: i32,
    o_ins: i32,
    e_ins: i32,
    w: i32,
    end_bonus: i32,
    zdrop: i32,
    h0: i32,
) -> Extend {
    checked_lengths(qlen, query, tlen, target, m, mat);
    assert!(h0 > 0, "BWA ksw_extend2 asserts h0 > 0");
    if qlen == 0 {
        return Extend {
            score: h0,
            qle: 0,
            tle: 0,
            gtle: 0,
            gscore: -1,
            max_off: 0,
        };
    }
    EXTEND_SCRATCH.with_borrow_mut(|scratch| {
        extend2_inner(
            scratch, qlen, query, tlen, target, m, mat, o_del, e_del, o_ins, e_ins, w, end_bonus,
            zdrop, h0,
        )
    })
}

#[allow(clippy::too_many_arguments)]
fn extend2_impl(
    qlen: i32,
    query: &[u8],
    tlen: i32,
    target: &[u8],
    m: i32,
    mat: &[i8],
    o_del: i32,
    e_del: i32,
    o_ins: i32,
    e_ins: i32,
    w: i32,
    end_bonus: i32,
    zdrop: i32,
    h0: i32,
) -> Extend {
    checked_lengths(qlen, query, tlen, target, m, mat);
    assert!(h0 > 0, "BWA ksw_extend2 asserts h0 > 0");

    // The C implementation has undefined behaviour for qlen == 0 (it stores
    // eh[1] after allocating only eh[0]). This makes the Rust API total for
    // the documented input domain, without claiming a C-oracle equivalence.
    if qlen == 0 {
        return Extend {
            score: h0,
            qle: 0,
            tle: 0,
            gtle: 0,
            gscore: -1,
            max_off: 0,
        };
    }

    EXTEND_SCRATCH.with_borrow_mut(|scratch| {
        if o_del == o_ins && e_del == e_ins {
            extend2_fast::<true>(
                scratch, qlen, query, tlen, target, m, mat, o_del, e_del, o_ins, e_ins, w,
                end_bonus, zdrop, h0,
            )
        } else {
            extend2_fast::<false>(
                scratch, qlen, query, tlen, target, m, mat, o_del, e_del, o_ins, e_ins, w,
                end_bonus, zdrop, h0,
            )
        }
    })
}

/// Same recurrence, comparisons and tie rules as [`extend2_inner`], with a
/// cheaper inner loop: H and E interleaved (one paired load per cell), a
/// single index counter, one compare for the row maximum and its column, a
/// branchless `M = H ? H + q : 0` (a branch here mispredicts at band edges),
/// and, when the deletion and insertion penalties are equal, one shared
/// `max(M - (o + e), 0)` for E and F. Kept exact by `extend_fast` (dumped
/// BWA calls), `tests::extend2_matches_reference` and `tests/ksw_vs_c.rs`.
#[allow(clippy::too_many_arguments)]
fn extend2_fast<const SAME_GAPS: bool>(
    scratch: &mut ExtendScratch,
    qlen: i32,
    query: &[u8],
    tlen: i32,
    target: &[u8],
    m: i32,
    mat: &[i8],
    o_del: i32,
    e_del: i32,
    o_ins: i32,
    e_ins: i32,
    mut w: i32,
    end_bonus: i32,
    zdrop: i32,
    h0: i32,
) -> Extend {
    let qlen_usize = qlen as usize;
    let tlen_usize = tlen as usize;
    let m_usize = m as usize;
    let oe_del = o_del.wrapping_add(e_del);
    let oe_ins = o_ins.wrapping_add(e_ins);
    scratch.qp.resize(qlen_usize * m_usize, 0);
    for (out, row) in scratch
        .qp
        .chunks_exact_mut(qlen_usize)
        .zip(mat.chunks_exact(m_usize))
    {
        for (score, &base) in out.iter_mut().zip(query) {
            *score = row[base as usize];
        }
    }
    scratch.eh.clear();
    scratch.eh.resize(qlen_usize + 1, Eh::default());
    let qp = &scratch.qp;
    let eh = &mut scratch.eh;

    eh[0].h = h0;
    eh[1].h = if h0 > oe_ins {
        h0.wrapping_sub(oe_ins)
    } else {
        0
    };
    let mut j = 2usize;
    while j <= qlen_usize && eh[j - 1].h > e_ins {
        eh[j].h = eh[j - 1].h.wrapping_sub(e_ins);
        j += 1;
    }

    let mut max_score = 0i32;
    for &score in mat {
        max_score = max_score.max(score as i32);
    }
    let max_ins = (((qlen
        .wrapping_mul(max_score)
        .wrapping_add(end_bonus)
        .wrapping_sub(o_ins)) as f64
        / e_ins as f64
        + 1.0) as i32)
        .max(1);
    w = w.min(max_ins);
    let max_del = (((qlen
        .wrapping_mul(max_score)
        .wrapping_add(end_bonus)
        .wrapping_sub(o_del)) as f64
        / e_del as f64
        + 1.0) as i32)
        .max(1);
    w = w.min(max_del);

    let mut max = h0;
    let mut max_i = -1i32;
    let mut max_j = -1i32;
    let mut max_ie = -1i32;
    let mut gscore = -1i32;
    let mut max_off = 0i32;
    let mut beg = 0i32;
    let mut end = qlen;

    for i_usize in 0..tlen_usize {
        let i = i_usize as i32;
        let q = &qp[target[i_usize] as usize * qlen_usize..][..qlen_usize];
        if beg < i.wrapping_sub(w) {
            beg = i.wrapping_sub(w);
        }
        if end > i.wrapping_add(w).wrapping_add(1) {
            end = i.wrapping_add(w).wrapping_add(1);
        }
        if end > qlen {
            end = qlen;
        }
        let mut h1 = if beg == 0 {
            let gap = o_del.wrapping_add(e_del.wrapping_mul(i.wrapping_add(1)));
            let value = h0.wrapping_sub(gap);
            if value < 0 { 0 } else { value }
        } else {
            0
        };
        let begin = beg as usize;
        let end_usize = end as usize;
        let row = &mut eh[begin..end_usize];
        let q = &q[begin..end_usize];
        let n = row.len().min(q.len());
        let mut f = 0i32;
        let mut row_max = 0i32;
        // Column offset of the row maximum; -1 - beg reproduces `mj = -1`.
        let mut mj_off = -1i32 - beg;
        for k in 0..n {
            let cell = row[k];
            let score = q[k] as i32;
            let diagonal =
                std::hint::select_unpredictable(cell.h != 0, cell.h.wrapping_add(score), 0);
            let h = diagonal.max(cell.e).max(f);
            // `mj = row_max > h ? mj : j; row_max = max(row_max, h)`.
            let take = h >= row_max;
            mj_off = if take { k as i32 } else { mj_off };
            row_max = if take { h } else { row_max };
            let t_del = diagonal.wrapping_sub(oe_del).max(0);
            let t_ins = if SAME_GAPS {
                t_del
            } else {
                diagonal.wrapping_sub(oe_ins).max(0)
            };
            row[k] = Eh {
                h: h1,
                e: cell.e.wrapping_sub(e_del).max(t_del),
            };
            f = f.wrapping_sub(e_ins).max(t_ins);
            h1 = h;
        }
        let mj = beg.wrapping_add(mj_off);
        eh[end_usize] = Eh { h: h1, e: 0 };
        if end == qlen {
            max_ie = if gscore > h1 { max_ie } else { i };
            gscore = gscore.max(h1);
        }
        if row_max == 0 {
            break;
        }
        if row_max > max {
            max = row_max;
            max_i = i;
            max_j = mj;
            max_off = max_off.max((mj - i).abs());
        } else if zdrop > 0 {
            let di = i - max_i;
            let dj = mj - max_j;
            if di > dj {
                if max
                    .wrapping_sub(row_max)
                    .wrapping_sub((di - dj).wrapping_mul(e_del))
                    > zdrop
                {
                    break;
                }
            } else if max
                .wrapping_sub(row_max)
                .wrapping_sub((dj - di).wrapping_mul(e_ins))
                > zdrop
            {
                break;
            }
        }
        let leading_zeroes = eh[begin..end_usize]
            .iter()
            .position(|c| c.h != 0 || c.e != 0)
            .unwrap_or(end_usize - begin);
        beg += leading_zeroes as i32;
        let last_nonzero = eh[beg as usize..=end_usize]
            .iter()
            .rposition(|c| c.h != 0 || c.e != 0)
            .map_or(beg - 1, |offset| beg + offset as i32);
        end = (last_nonzero + 2).min(qlen);
    }

    Extend {
        score: max,
        qle: max_j + 1,
        tle: max_i + 1,
        gtle: max_ie + 1,
        gscore,
        max_off,
    }
}

#[allow(clippy::too_many_arguments)]
fn extend2_inner(
    scratch: &mut ExtendScratch,
    qlen: i32,
    query: &[u8],
    tlen: i32,
    target: &[u8],
    m: i32,
    mat: &[i8],
    o_del: i32,
    e_del: i32,
    o_ins: i32,
    e_ins: i32,
    mut w: i32,
    end_bonus: i32,
    zdrop: i32,
    h0: i32,
) -> Extend {
    let qlen_usize = qlen as usize;
    let tlen_usize = tlen as usize;
    let m_usize = m as usize;
    let oe_del = o_del.wrapping_add(e_del);
    let oe_ins = o_ins.wrapping_add(e_ins);
    scratch.qp.resize(qlen_usize * m_usize, 0);
    for (out, row) in scratch
        .qp
        .chunks_exact_mut(qlen_usize)
        .zip(mat.chunks_exact(m_usize))
    {
        for (score, &base) in out.iter_mut().zip(query) {
            *score = row[base as usize];
        }
    }
    scratch.h.resize(qlen_usize + 1, 0);
    scratch.e.resize(qlen_usize + 1, 0);
    scratch.h.fill(0);
    scratch.e.fill(0);
    let qp = &scratch.qp;
    let h_buf = &mut scratch.h;
    let e_buf = &mut scratch.e;

    h_buf[0] = h0;
    h_buf[1] = if h0 > oe_ins {
        h0.wrapping_sub(oe_ins)
    } else {
        0
    };
    let mut j = 2usize;
    while j <= qlen_usize && h_buf[j - 1] > e_ins {
        h_buf[j] = h_buf[j - 1].wrapping_sub(e_ins);
        j += 1;
    }

    let mut max_score = 0i32;
    for &score in mat {
        max_score = max_score.max(score as i32);
    }
    let max_ins = (((qlen
        .wrapping_mul(max_score)
        .wrapping_add(end_bonus)
        .wrapping_sub(o_ins)) as f64
        / e_ins as f64
        + 1.0) as i32)
        .max(1);
    w = w.min(max_ins);
    let max_del = (((qlen
        .wrapping_mul(max_score)
        .wrapping_add(end_bonus)
        .wrapping_sub(o_del)) as f64
        / e_del as f64
        + 1.0) as i32)
        .max(1);
    w = w.min(max_del);

    let mut max = h0;
    let mut max_i = -1i32;
    let mut max_j = -1i32;
    let mut max_ie = -1i32;
    let mut gscore = -1i32;
    let mut max_off = 0i32;
    let mut beg = 0i32;
    let mut end = qlen;

    for i_usize in 0..tlen_usize {
        let i = i_usize as i32;
        let mut f = 0i32;
        let q = &qp[target[i_usize] as usize * qlen_usize..][..qlen_usize];
        if beg < i.wrapping_sub(w) {
            beg = i.wrapping_sub(w);
        }
        if end > i.wrapping_add(w).wrapping_add(1) {
            end = i.wrapping_add(w).wrapping_add(1);
        }
        if end > qlen {
            end = qlen;
        }
        let mut h1 = if beg == 0 {
            let gap = o_del.wrapping_add(e_del.wrapping_mul(i.wrapping_add(1)));
            let value = h0.wrapping_sub(gap);
            if value < 0 { 0 } else { value }
        } else {
            0
        };
        let mut row_max = 0i32;
        let mut mj = -1i32;
        let begin = beg as usize;
        let end_usize = end as usize;
        for (offset, ((h_slot, e_slot), &score)) in h_buf[begin..end_usize]
            .iter_mut()
            .zip(&mut e_buf[begin..end_usize])
            .zip(&q[begin..end_usize])
            .enumerate()
        {
            let jj = beg + offset as i32;
            let mut diagonal = *h_slot;
            let mut e = *e_slot;
            *h_slot = h1;
            diagonal = if diagonal != 0 {
                diagonal.wrapping_add(score as i32)
            } else {
                0
            };
            let mut h = diagonal.max(e);
            h = h.max(f);
            h1 = h;
            mj = if row_max > h { mj } else { jj };
            row_max = row_max.max(h);
            let mut value = diagonal.wrapping_sub(oe_del);
            value = value.max(0);
            e = e.wrapping_sub(e_del).max(value);
            *e_slot = e;
            value = diagonal.wrapping_sub(oe_ins);
            value = value.max(0);
            f = f.wrapping_sub(e_ins).max(value);
        }
        h_buf[end_usize] = h1;
        e_buf[end_usize] = 0;
        if end == qlen {
            max_ie = if gscore > h1 { max_ie } else { i };
            gscore = gscore.max(h1);
        }
        if row_max == 0 {
            break;
        }
        if row_max > max {
            max = row_max;
            max_i = i;
            max_j = mj;
            max_off = max_off.max((mj - i).abs());
        } else if zdrop > 0 {
            let di = i - max_i;
            let dj = mj - max_j;
            if di > dj {
                if max
                    .wrapping_sub(row_max)
                    .wrapping_sub((di - dj).wrapping_mul(e_del))
                    > zdrop
                {
                    break;
                }
            } else if max
                .wrapping_sub(row_max)
                .wrapping_sub((dj - di).wrapping_mul(e_ins))
                > zdrop
            {
                break;
            }
        }
        let leading_zeroes = h_buf[begin..end_usize]
            .iter()
            .zip(&e_buf[begin..end_usize])
            .position(|(&h, &e)| h != 0 || e != 0)
            .unwrap_or(end_usize - begin);
        beg += leading_zeroes as i32;
        let last_nonzero = h_buf[beg as usize..=end_usize]
            .iter()
            .zip(&e_buf[beg as usize..=end_usize])
            .rposition(|(&h, &e)| h != 0 || e != 0)
            .map_or(beg - 1, |offset| beg + offset as i32);
        end = (last_nonzero + 2).min(qlen);
    }

    Extend {
        score: max,
        qle: max_j + 1,
        tle: max_i + 1,
        gtle: max_ie + 1,
        gscore,
        max_off,
    }
}

/// BWA's symmetric-gap `ksw_extend` wrapper.
#[allow(clippy::too_many_arguments)]
pub fn extend(
    qlen: i32,
    query: &[u8],
    tlen: i32,
    target: &[u8],
    m: i32,
    mat: &[i8],
    gapo: i32,
    gape: i32,
    w: i32,
    end_bonus: i32,
    zdrop: i32,
    h0: i32,
) -> Extend {
    extend2(
        qlen, query, tlen, target, m, mat, gapo, gape, gapo, gape, w, end_bonus, zdrop, h0,
    )
}

/// BWA's `ksw_global2`; CIGAR words use BAM's `len << 4 | op` encoding.
#[allow(clippy::too_many_arguments)]
pub fn global2(
    qlen: i32,
    query: &[u8],
    tlen: i32,
    target: &[u8],
    m: i32,
    mat: &[i8],
    o_del: i32,
    e_del: i32,
    o_ins: i32,
    e_ins: i32,
    w: i32,
) -> (i32, Vec<u32>) {
    checked_lengths(qlen, query, tlen, target, m, mat);
    assert!(w >= 0);
    let qlen_usize = qlen as usize;
    let tlen_usize = tlen as usize;
    let m_usize = m as usize;
    let oe_del = o_del.wrapping_add(e_del);
    let oe_ins = o_ins.wrapping_add(e_ins);
    let n_col = qlen.min(w.saturating_mul(2).saturating_add(1)) as usize;
    let mut sc = GLOBAL_SCRATCH.take();
    let GlobalScratch { qp, eh, z } = &mut sc;
    global_init(qp, eh, qlen, query, m_usize, mat, o_ins, e_ins, w);
    // Zeroed as a fresh allocation would be.
    z.clear();
    z.resize(n_col * tlen_usize, 0);
    for i_usize in 0..tlen_usize {
        let i = i_usize as i32;
        let mut f = MINUS_INF;
        let q = &qp[target[i_usize] as usize * qlen_usize..][..qlen_usize];
        let beg = if i > w { i - w } else { 0 };
        let end = if i.wrapping_add(w).wrapping_add(1) < qlen {
            i + w + 1
        } else {
            qlen
        };
        let mut h1 = if beg == 0 {
            -(o_del.wrapping_add(e_del.wrapping_mul(i.wrapping_add(1))))
        } else {
            MINUS_INF
        };
        // Columns [beg, end) of row i; their directions go to z[i][j - beg]
        // (end - beg <= n_col, and both are 0 when qlen is).
        // When i - w passes qlen the band is empty (beg > end): ksw.c's loop
        // runs zero times and only eh[end] is written below.
        let (b, e_) = (beg as usize, end as usize);
        let b = b.min(e_);
        let zi = &mut z[i_usize * n_col..][..e_ - b];
        for ((slot, &qv), zd) in eh[b..e_].iter_mut().zip(&q[b..e_]).zip(zi) {
            let mut diagonal = slot.h;
            let mut e = slot.e;
            slot.h = h1;
            diagonal = diagonal.wrapping_add(qv as i32);
            let mut direction = if diagonal >= e { 0u8 } else { 1u8 };
            let mut h = diagonal.max(e);
            direction = if h >= f { direction } else { 2 };
            h = h.max(f);
            h1 = h;
            let candidate = diagonal.wrapping_sub(oe_del);
            e = e.wrapping_sub(e_del);
            direction |= if e > candidate { 1 << 2 } else { 0 };
            e = e.max(candidate);
            slot.e = e;
            let candidate = diagonal.wrapping_sub(oe_ins);
            f = f.wrapping_sub(e_ins);
            direction |= if f > candidate { 2 << 4 } else { 0 };
            f = f.max(candidate);
            *zd = direction;
        }
        eh[e_] = Eh {
            h: h1,
            e: MINUS_INF,
        };
    }
    let score = eh[qlen_usize].h;
    let z = &sc.z;

    let mut cigar = Vec::<u32>::new();
    if qlen != 0 && tlen != 0 {
        let mut i = tlen - 1;
        let mut k = (if i + w + 1 < qlen { i + w + 1 } else { qlen }) - 1;
        let mut which = 0u8;
        while i >= 0 && k >= 0 {
            // Negative past an empty band row; C indexes `z` with signed
            // arithmetic, landing in the previous row, so do the same.
            let offset = (k - if i > w { i - w } else { 0 }) as isize;
            which = (z[(i as isize * n_col as isize + offset) as usize] >> (which << 1)) & 3;
            if which == 0 {
                push_cigar(&mut cigar, 0, 1);
                i -= 1;
                k -= 1;
            } else if which == 1 {
                push_cigar(&mut cigar, 2, 1);
                i -= 1;
            } else {
                push_cigar(&mut cigar, 1, 1);
                k -= 1;
            }
        }
        if i >= 0 {
            push_cigar(&mut cigar, 2, (i + 1) as u32);
        }
        if k >= 0 {
            push_cigar(&mut cigar, 1, (k + 1) as u32);
        }
        cigar.reverse();
    } else if tlen != 0 {
        push_cigar(&mut cigar, 2, tlen as u32);
    } else if qlen != 0 {
        push_cigar(&mut cigar, 1, qlen as u32);
    }
    GLOBAL_SCRATCH.set(sc);
    (score, cigar)
}

/// `ksw_global2`'s query profile and first DP row (ksw.c:550-561), into
/// reused buffers.
#[allow(clippy::too_many_arguments)]
fn global_init(
    qp: &mut Vec<i8>,
    eh: &mut Vec<Eh>,
    qlen: i32,
    query: &[u8],
    m: usize,
    mat: &[i8],
    o_ins: i32,
    e_ins: i32,
    w: i32,
) {
    qp.clear();
    for k in 0..m {
        let row = &mat[k * m..(k + 1) * m];
        qp.extend(query.iter().map(|&base| row[base as usize]));
    }
    eh.clear();
    eh.resize(qlen as usize + 1, Eh::default());
    eh[0] = Eh { h: 0, e: MINUS_INF };
    let mut j = 1i32;
    while j <= qlen && j <= w {
        eh[j as usize] = Eh {
            h: -(o_ins.wrapping_add(e_ins.wrapping_mul(j))),
            e: MINUS_INF,
        };
        j += 1;
    }
    while j <= qlen {
        eh[j as usize] = Eh {
            h: MINUS_INF,
            e: MINUS_INF,
        };
        j += 1;
    }
}

/// Per-thread buffers of `global2` and `global2_score`.
#[derive(Default)]
struct GlobalScratch {
    qp: Vec<i8>,
    eh: Vec<Eh>,
    z: Vec<u8>,
}

thread_local! {
    static GLOBAL_SCRATCH: std::cell::Cell<GlobalScratch> = const {
        std::cell::Cell::new(GlobalScratch {
            qp: Vec::new(),
            eh: Vec::new(),
            z: Vec::new(),
        })
    };
}

/// `ksw_global2` called without a CIGAR (`n_cigar_ == 0`, as
/// `mem_patch_reg` does through `bwa_gen_cigar2`): BWA then skips the
/// backtrack matrix and runs the same recurrence in its second inner loop.
/// Returns the score, which `global2` also computes before backtracking.
#[allow(clippy::too_many_arguments)]
pub fn global2_score(
    qlen: i32,
    query: &[u8],
    tlen: i32,
    target: &[u8],
    m: i32,
    mat: &[i8],
    o_del: i32,
    e_del: i32,
    o_ins: i32,
    e_ins: i32,
    w: i32,
) -> i32 {
    checked_lengths(qlen, query, tlen, target, m, mat);
    assert!(w >= 0);
    let qlen_usize = qlen as usize;
    let m_usize = m as usize;
    let oe_del = o_del.wrapping_add(e_del);
    let oe_ins = o_ins.wrapping_add(e_ins);
    let mut sc = GLOBAL_SCRATCH.take();
    let score = {
        let GlobalScratch { qp, eh, .. } = &mut sc;
        global_init(qp, eh, qlen, query, m_usize, mat, o_ins, e_ins, w);
        for (i_usize, &base) in target.iter().enumerate() {
            let i = i_usize as i32;
            let mut f = MINUS_INF;
            let q = &qp[base as usize * qlen_usize..][..qlen_usize];
            let beg = if i > w { i - w } else { 0 };
            let end = if i.wrapping_add(w).wrapping_add(1) < qlen {
                i + w + 1
            } else {
                qlen
            };
            let mut h1 = if beg == 0 {
                -(o_del.wrapping_add(e_del.wrapping_mul(i.wrapping_add(1))))
            } else {
                MINUS_INF
            };
            for jj in beg as usize..end as usize {
                // ksw.c:604-617, the same recurrence as `global2` minus `d`.
                let slot = &mut eh[jj];
                let (m, mut e) = (slot.h, slot.e);
                slot.h = h1;
                let m = m.wrapping_add(q[jj] as i32);
                let h = m.max(e).max(f);
                h1 = h;
                let t = m.wrapping_sub(oe_del);
                e = e.wrapping_sub(e_del).max(t);
                slot.e = e;
                let t = m.wrapping_sub(oe_ins);
                f = f.wrapping_sub(e_ins).max(t);
            }
            eh[end as usize] = Eh {
                h: h1,
                e: MINUS_INF,
            };
        }
        eh[qlen_usize].h
    };
    GLOBAL_SCRATCH.set(sc);
    score
}

fn push_cigar(cigar: &mut Vec<u32>, op: u32, len: u32) {
    if cigar.last().is_none_or(|word| word & 0xf != op) {
        cigar.push((len << 4) | op);
    } else {
        *cigar.last_mut().expect("non-empty after last check") += len << 4;
    }
}

/// BWA's symmetric-gap `ksw_global` wrapper.
#[allow(clippy::too_many_arguments)]
pub fn global(
    qlen: i32,
    query: &[u8],
    tlen: i32,
    target: &[u8],
    m: i32,
    mat: &[i8],
    gapo: i32,
    gape: i32,
    w: i32,
) -> (i32, Vec<u32>) {
    global2(qlen, query, tlen, target, m, mat, gapo, gape, gapo, gape, w)
}

/// Use byte-wide scores in [`align2`].
pub const KSW_XBYTE: i32 = 0x1_0000;
/// Stop once the best score reaches the low 16 bits of `xtra`.
pub const KSW_XSTOP: i32 = 0x2_0000;
/// Find the second-best score above the low 16 bits of `xtra`.
pub const KSW_XSUBO: i32 = 0x4_0000;
/// Find the alignment start positions with BWA's reverse pass.
pub const KSW_XSTART: i32 = 0x8_0000;

/// Result of BWA 0.7.19's `ksw_align2`.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Kswr {
    pub score: i32,
    pub te: i32,
    pub qe: i32,
    pub score2: i32,
    pub te2: i32,
    pub tb: i32,
    pub qb: i32,
}

impl Default for Kswr {
    fn default() -> Self {
        Self {
            score: 0,
            te: -1,
            qe: -1,
            score2: -1,
            te2: -1,
            tb: -1,
            qb: -1,
        }
    }
}

type U8x16 = [u8; 16];
type I16x8 = [i16; 8];

enum ProfileData {
    U8 {
        qp: Vec<U8x16>,
        h0: Vec<U8x16>,
        h1: Vec<U8x16>,
        e: Vec<U8x16>,
        hmax: Vec<U8x16>,
    },
    I16 {
        qp: Vec<I16x8>,
        h0: Vec<I16x8>,
        h1: Vec<I16x8>,
        e: Vec<I16x8>,
        hmax: Vec<I16x8>,
    },
}

/// Reusable equivalent of BWA's `kswq_t` query profile.
///
/// `score_size` is normalized in the same way as `ksw_qinit`: values greater
/// than one select signed 16-bit lanes; every other value selects byte lanes.
pub struct Profile {
    qlen: usize,
    slen: usize,
    size: u8,
    shift: u8,
    max: u8,
    data: ProfileData,
    rows: Vec<(i32, i32)>,
}

impl Profile {
    pub fn new(score_size: u8, qlen: i32, query: &[u8], m: i32, mat: &[i8]) -> Self {
        assert!(qlen > 0, "BWA ksw_qinit requires a non-empty query");
        assert!(m > 0);
        assert_eq!(query.len(), qlen as usize);
        assert_eq!(mat.len(), (m as usize) * (m as usize));
        assert!(query.iter().all(|&base| base < m as u8));

        let size = if score_size > 1 { 2 } else { 1 };
        let lanes = 8 * (3 - size as usize);
        let qlen = qlen as usize;
        let slen = qlen.div_ceil(lanes);
        let mut minimum = 127i8;
        let mut maximum = 0i8;
        for &score in mat {
            if score < minimum {
                minimum = score;
            }
            if score > maximum {
                maximum = score;
            }
        }
        let shift = 0u8.wrapping_sub(minimum as u8);
        let max = maximum as u8;

        let data = if size == 1 {
            let mut qp = vec![[0; 16]; slen * m as usize];
            for a in 0..m as usize {
                for i in 0..slen {
                    for (lane, value) in qp[a * slen + i].iter_mut().enumerate() {
                        let k = i + lane * slen;
                        let score = if k >= qlen {
                            0
                        } else {
                            mat[a * m as usize + query[k] as usize]
                        };
                        *value = (score as u8).wrapping_add(shift);
                    }
                }
            }
            ProfileData::U8 {
                qp,
                h0: vec![[0; 16]; slen],
                h1: vec![[0; 16]; slen],
                e: vec![[0; 16]; slen],
                hmax: vec![[0; 16]; slen],
            }
        } else {
            let mut qp = vec![[0; 8]; slen * m as usize];
            for a in 0..m as usize {
                for i in 0..slen {
                    for (lane, value) in qp[a * slen + i].iter_mut().enumerate() {
                        let k = i + lane * slen;
                        *value = if k >= qlen {
                            0
                        } else {
                            mat[a * m as usize + query[k] as usize] as i16
                        };
                    }
                }
            }
            ProfileData::I16 {
                qp,
                h0: vec![[0; 8]; slen],
                h1: vec![[0; 8]; slen],
                e: vec![[0; 8]; slen],
                hmax: vec![[0; 8]; slen],
            }
        };
        Self {
            qlen,
            slen,
            size,
            shift,
            max,
            data,
            rows: Vec::new(),
        }
    }

    fn reset(&mut self, score_size: u8, qlen: i32, query: &[u8], m: i32, mat: &[i8]) {
        let size = if score_size > 1 { 2 } else { 1 };
        if size != self.size {
            *self = Self::new(score_size, qlen, query, m, mat);
            return;
        }
        let qlen = qlen as usize;
        let m = m as usize;
        let lanes = 8 * (3 - size as usize);
        let slen = qlen.div_ceil(lanes);
        let mut minimum = 127i8;
        let mut maximum = 0i8;
        for &score in mat {
            minimum = minimum.min(score);
            maximum = maximum.max(score);
        }
        let shift = 0u8.wrapping_sub(minimum as u8);

        match &mut self.data {
            ProfileData::U8 {
                qp,
                h0,
                h1,
                e,
                hmax,
            } => {
                qp.resize(slen * m, [0; 16]);
                for a in 0..m {
                    for i in 0..slen {
                        for (lane, value) in qp[a * slen + i].iter_mut().enumerate() {
                            let k = i + lane * slen;
                            let score = if k >= qlen {
                                0
                            } else {
                                mat[a * m + query[k] as usize]
                            };
                            *value = (score as u8).wrapping_add(shift);
                        }
                    }
                }
                h0.resize(slen, [0; 16]);
                h1.resize(slen, [0; 16]);
                e.resize(slen, [0; 16]);
                hmax.resize(slen, [0; 16]);
            }
            ProfileData::I16 {
                qp,
                h0,
                h1,
                e,
                hmax,
            } => {
                qp.resize(slen * m, [0; 8]);
                for a in 0..m {
                    for i in 0..slen {
                        for (lane, value) in qp[a * slen + i].iter_mut().enumerate() {
                            let k = i + lane * slen;
                            *value = if k >= qlen {
                                0
                            } else {
                                mat[a * m + query[k] as usize] as i16
                            };
                        }
                    }
                }
                h0.resize(slen, [0; 8]);
                h1.resize(slen, [0; 8]);
                e.resize(slen, [0; 8]);
                hmax.resize(slen, [0; 8]);
            }
        }
        self.qlen = qlen;
        self.slen = slen;
        self.shift = shift;
        self.max = maximum as u8;
        self.rows.clear();
    }

    pub fn score_size(&self) -> u8 {
        self.size
    }

    pub fn query_len(&self) -> usize {
        self.qlen
    }
}

#[derive(Default)]
struct AlignScratch {
    forward_u8: Option<Profile>,
    forward_i16: Option<Profile>,
    reverse_u8: Option<Profile>,
    reverse_i16: Option<Profile>,
    reverse_query: Vec<u8>,
    reverse_target: Vec<u8>,
}

thread_local! {
    static ALIGN_SCRATCH: RefCell<AlignScratch> = RefCell::default();
}

fn scratch_profile<'a>(
    slot: &'a mut Option<Profile>,
    score_size: u8,
    query: &[u8],
    m: i32,
    mat: &[i8],
) -> &'a mut Profile {
    if let Some(profile) = slot {
        profile.reset(score_size, query.len() as i32, query, m, mat);
        profile
    } else {
        slot.insert(Profile::new(score_size, query.len() as i32, query, m, mat))
    }
}

#[cfg(not(target_arch = "aarch64"))]
#[inline(always)]
fn shift_u8(a: U8x16) -> U8x16 {
    let mut out = [0; 16];
    out[1..].copy_from_slice(&a[..15]);
    out
}

#[cfg(not(target_arch = "aarch64"))]
#[inline(always)]
fn shift_i16(a: I16x8) -> I16x8 {
    let mut out = [0; 8];
    out[1..].copy_from_slice(&a[..7]);
    out
}

#[cfg(not(target_arch = "aarch64"))]
#[inline(never)]
#[allow(clippy::too_many_arguments)]
fn step_u8(
    h: &mut U8x16,
    scores: &U8x16,
    ev: &mut U8x16,
    f: &mut U8x16,
    row_max: &mut U8x16,
    shift: u8,
    e_del: u8,
    oe_del: u8,
    e_ins: u8,
    oe_ins: u8,
) {
    macro_rules! lane {
        ($lane:expr) => {{
            let value = h[$lane]
                .saturating_add(scores[$lane])
                .saturating_sub(shift)
                .max(ev[$lane])
                .max(f[$lane]);
            row_max[$lane] = row_max[$lane].max(value);
            ev[$lane] = ev[$lane]
                .saturating_sub(e_del)
                .max(value.saturating_sub(oe_del));
            f[$lane] = f[$lane]
                .saturating_sub(e_ins)
                .max(value.saturating_sub(oe_ins));
            h[$lane] = value;
        }};
    }
    lane!(0);
    lane!(1);
    lane!(2);
    lane!(3);
    lane!(4);
    lane!(5);
    lane!(6);
    lane!(7);
    lane!(8);
    lane!(9);
    lane!(10);
    lane!(11);
    lane!(12);
    lane!(13);
    lane!(14);
    lane!(15);
}

#[cfg(not(target_arch = "aarch64"))]
#[inline(never)]
fn lazy_f_u8(h: &mut U8x16, f: &mut U8x16, e_ins: u8, oe_ins: u8) -> bool {
    let mut active = false;
    for (h, f) in h.iter_mut().zip(f) {
        *h = (*h).max(*f);
        let h_gap = h.saturating_sub(oe_ins);
        *f = f.saturating_sub(e_ins);
        active |= *f > h_gap;
    }
    active
}

#[cfg(not(target_arch = "aarch64"))]
#[inline(never)]
#[allow(clippy::too_many_arguments)]
fn step_i16(
    h: &mut I16x8,
    scores: &I16x8,
    ev: &mut I16x8,
    f: &mut I16x8,
    row_max: &mut I16x8,
    e_del: u16,
    oe_del: u16,
    e_ins: u16,
    oe_ins: u16,
) {
    for ((((h, &score), e), f), row_max) in h.iter_mut().zip(scores).zip(ev).zip(f).zip(row_max) {
        let value = h.saturating_add(score).max(*e).max(*f);
        *row_max = (*row_max).max(value);
        *e = ((*e as u16).saturating_sub(e_del) as i16)
            .max((value as u16).saturating_sub(oe_del) as i16);
        *f = ((*f as u16).saturating_sub(e_ins) as i16)
            .max((value as u16).saturating_sub(oe_ins) as i16);
        *h = value;
    }
}

#[cfg(not(target_arch = "aarch64"))]
#[inline(never)]
fn lazy_f_i16(h: &mut I16x8, f: &mut I16x8, e_ins: u16, oe_ins: u16) -> bool {
    let mut active = false;
    for (h, f) in h.iter_mut().zip(f) {
        *h = (*h).max(*f);
        let h_gap = h.cast_unsigned().saturating_sub(oe_ins) as i16;
        *f = f.cast_unsigned().saturating_sub(e_ins) as i16;
        active |= *f > h_gap;
    }
    active
}

fn record_row(rows: &mut Vec<(i32, i32)>, score: i32, target_end: i32) {
    if let Some(last) = rows.last_mut()
        && last.1 + 1 == target_end
    {
        if last.0 < score {
            *last = (score, target_end);
        }
    } else {
        rows.push((score, target_end));
    }
}

fn second_best(r: &mut Kswr, rows: &[(i32, i32)], max_score: u8) {
    if rows.is_empty() {
        return;
    }
    let radius = (r.score + max_score as i32 - 1) / max_score as i32;
    let low = r.te - radius;
    let high = r.te + radius;
    for &(score, end) in rows {
        if (end < low || end > high) && score > r.score2 {
            r.score2 = score;
            r.te2 = end;
        }
    }
}

fn ksw_u8(
    profile: &mut Profile,
    target: &[u8],
    o_del: i32,
    e_del: i32,
    o_ins: i32,
    e_ins: i32,
    xtra: i32,
) -> Kswr {
    let slen = profile.slen;
    let shift = profile.shift;
    let profile_max = profile.max;
    let minsc = if xtra & KSW_XSUBO != 0 {
        xtra & 0xffff
    } else {
        0x1_0000
    };
    let endsc = if xtra & KSW_XSTOP != 0 {
        xtra & 0xffff
    } else {
        0x1_0000
    };
    let oe_del = o_del.wrapping_add(e_del) as u8;
    let e_del = e_del as u8;
    let oe_ins = o_ins.wrapping_add(e_ins) as u8;
    let e_ins = e_ins as u8;
    profile.rows.clear();
    let mut result = Kswr::default();
    let mut gmax = 0i32;

    let ProfileData::U8 {
        qp,
        h0,
        h1,
        e,
        hmax,
    } = &mut profile.data
    else {
        unreachable!()
    };
    e.fill([0; 16]);
    h0.fill([0; 16]);
    hmax.fill([0; 16]);

    for (i, &target_base) in target.iter().enumerate() {
        let scores = &qp[target_base as usize * slen..][..slen];
        #[cfg(target_arch = "aarch64")]
        let imax = nomics_simd::ksw::row_u8(
            scores,
            h0,
            h1,
            e,
            shift,
            nomics_simd::ksw::Gaps {
                e_del,
                oe_del,
                e_ins,
                oe_ins,
            },
        ) as i32;
        #[cfg(not(target_arch = "aarch64"))]
        let imax = {
            let mut f = [0; 16];
            let mut row_max = [0; 16];
            let mut h = shift_u8(h0[slen - 1]);
            for j in 0..slen {
                let mut ev = e[j];
                step_u8(
                    &mut h,
                    &scores[j],
                    &mut ev,
                    &mut f,
                    &mut row_max,
                    shift,
                    e_del,
                    oe_del,
                    e_ins,
                    oe_ins,
                );
                h1[j] = h;
                e[j] = ev;
                h = h0[j];
            }
            'lazy_f: for _ in 0..16 {
                f = shift_u8(f);
                for cell in h1.iter_mut().take(slen) {
                    if !lazy_f_u8(cell, &mut f, e_ins, oe_ins) {
                        break 'lazy_f;
                    }
                }
            }
            *row_max.iter().max().expect("non-empty SIMD vector") as i32
        };
        if imax >= minsc {
            record_row(&mut profile.rows, imax, i as i32);
        }
        if imax > gmax {
            gmax = imax;
            result.te = i as i32;
            hmax.copy_from_slice(h1);
            if gmax + shift as i32 >= 255 || gmax >= endsc {
                break;
            }
        }
        std::mem::swap(h0, h1);
    }
    result.score = if gmax + (shift as i32) < 255 {
        gmax
    } else {
        255
    };
    if result.score != 255 {
        let mut maximum = -1i32;
        for i in 0..slen * 16 {
            let value = hmax[i / 16][i % 16] as i32;
            let query_end = (i / 16 + i % 16 * slen) as i32;
            if value > maximum {
                maximum = value;
                result.qe = query_end;
            } else if value == maximum && query_end < result.qe {
                result.qe = query_end;
            }
        }
        second_best(&mut result, &profile.rows, profile_max);
    }
    result
}

fn ksw_i16(
    profile: &mut Profile,
    target: &[u8],
    o_del: i32,
    e_del: i32,
    o_ins: i32,
    e_ins: i32,
    xtra: i32,
) -> Kswr {
    let slen = profile.slen;
    let profile_max = profile.max;
    let minsc = if xtra & KSW_XSUBO != 0 {
        xtra & 0xffff
    } else {
        0x1_0000
    };
    let endsc = if xtra & KSW_XSTOP != 0 {
        xtra & 0xffff
    } else {
        0x1_0000
    };
    let oe_del = o_del.wrapping_add(e_del) as u16;
    let e_del = e_del as u16;
    let oe_ins = o_ins.wrapping_add(e_ins) as u16;
    let e_ins = e_ins as u16;
    profile.rows.clear();
    let mut result = Kswr::default();
    let mut gmax = 0i32;

    let ProfileData::I16 {
        qp,
        h0,
        h1,
        e,
        hmax,
    } = &mut profile.data
    else {
        unreachable!()
    };
    e.fill([0; 8]);
    h0.fill([0; 8]);
    hmax.fill([0; 8]);

    for (i, &target_base) in target.iter().enumerate() {
        let scores = &qp[target_base as usize * slen..][..slen];
        #[cfg(target_arch = "aarch64")]
        let imax = nomics_simd::ksw::row_i16(
            scores,
            h0,
            h1,
            e,
            nomics_simd::ksw::Gaps {
                e_del,
                oe_del,
                e_ins,
                oe_ins,
            },
        ) as i32;
        #[cfg(not(target_arch = "aarch64"))]
        let imax = {
            let mut f = [0; 8];
            let mut row_max = [0; 8];
            let mut h = shift_i16(h0[slen - 1]);
            for j in 0..slen {
                let mut ev = e[j];
                step_i16(
                    &mut h,
                    &scores[j],
                    &mut ev,
                    &mut f,
                    &mut row_max,
                    e_del,
                    oe_del,
                    e_ins,
                    oe_ins,
                );
                h1[j] = h;
                e[j] = ev;
                h = h0[j];
            }
            'lazy_f: for _ in 0..16 {
                f = shift_i16(f);
                for cell in h1.iter_mut().take(slen) {
                    if !lazy_f_i16(cell, &mut f, e_ins, oe_ins) {
                        break 'lazy_f;
                    }
                }
            }
            *row_max.iter().max().expect("non-empty SIMD vector") as i32
        };
        if imax >= minsc {
            record_row(&mut profile.rows, imax, i as i32);
        }
        if imax > gmax {
            gmax = imax;
            result.te = i as i32;
            hmax.copy_from_slice(h1);
            if gmax >= endsc {
                break;
            }
        }
        std::mem::swap(h0, h1);
    }
    result.score = gmax;
    let mut maximum = -1i32;
    result.qe = -1;
    for i in 0..slen * 8 {
        let value = hmax[i / 8][i % 8] as u16 as i32;
        let query_end = (i / 8 + i % 8 * slen) as i32;
        if value > maximum {
            maximum = value;
            result.qe = query_end;
        } else if value == maximum && query_end < result.qe {
            result.qe = query_end;
        }
    }
    second_best(&mut result, &profile.rows, profile_max);
    result
}

fn run_profile(
    profile: &mut Profile,
    target: &[u8],
    o_del: i32,
    e_del: i32,
    o_ins: i32,
    e_ins: i32,
    xtra: i32,
) -> Kswr {
    if profile.size == 2 {
        ksw_i16(profile, target, o_del, e_del, o_ins, e_ins, xtra)
    } else {
        ksw_u8(profile, target, o_del, e_del, o_ins, e_ins, xtra)
    }
}

/// BWA 0.7.19's striped `ksw_align2`.
///
/// Passing a profile reuses it exactly like a non-null `kswq_t **`: its score
/// width and query profile take precedence over `query`, `mat`, and `KSW_XBYTE`.
#[allow(clippy::too_many_arguments)]
pub fn align2(
    qlen: i32,
    query: &[u8],
    tlen: i32,
    target: &[u8],
    m: i32,
    mat: &[i8],
    o_del: i32,
    e_del: i32,
    o_ins: i32,
    e_ins: i32,
    xtra: i32,
    qry: Option<&mut Profile>,
) -> Kswr {
    let Some(sink) = dump::sink(dump::ALIGN2) else {
        return align2_impl(
            qlen, query, tlen, target, m, mat, o_del, e_del, o_ins, e_ins, xtra, qry,
        );
    };
    // Effective score width: a passed profile's width overrides KSW_XBYTE.
    let size = match &qry {
        Some(p) => p.size,
        None if xtra & KSW_XBYTE != 0 => 1,
        None => 2,
    };
    let r = align2_impl(
        qlen, query, tlen, target, m, mat, o_del, e_del, o_ins, e_ins, xtra, qry,
    );
    let mut b = Vec::with_capacity(80 + query.len() + target.len());
    dump::seq(&mut b, &query[..qlen as usize]);
    dump::seq(&mut b, &target[..tlen as usize]);
    dump::ints(&mut b, &[m]);
    b.extend(mat[..(m * m) as usize].iter().map(|&x| x as u8));
    dump::ints(&mut b, &[o_del, e_del, o_ins, e_ins, xtra, i32::from(size)]);
    dump::ints(&mut b, &[r.score, r.te, r.qe, r.score2, r.te2, r.tb, r.qb]);
    dump::write(sink, &b);
    r
}

#[allow(clippy::too_many_arguments)]
fn align2_impl(
    qlen: i32,
    query: &[u8],
    tlen: i32,
    target: &[u8],
    m: i32,
    mat: &[i8],
    o_del: i32,
    e_del: i32,
    o_ins: i32,
    e_ins: i32,
    xtra: i32,
    qry: Option<&mut Profile>,
) -> Kswr {
    checked_lengths(qlen, query, tlen, target, m, mat);
    assert!(qlen > 0, "BWA ksw_align2 requires a non-empty query");
    let (mut result, size) = if let Some(profile) = qry {
        (
            run_profile(profile, target, o_del, e_del, o_ins, e_ins, xtra),
            profile.size,
        )
    } else {
        let size = if xtra & KSW_XBYTE != 0 { 1 } else { 2 };
        let result = ALIGN_SCRATCH.with_borrow_mut(|scratch| {
            let slot = if size == 1 {
                &mut scratch.forward_u8
            } else {
                &mut scratch.forward_i16
            };
            let profile = scratch_profile(slot, size, query, m, mat);
            run_profile(profile, target, o_del, e_del, o_ins, e_ins, xtra)
        });
        (result, size)
    };
    if xtra & KSW_XSTART == 0 || (xtra & KSW_XSUBO != 0 && result.score < xtra & 0xffff) {
        return result;
    }

    let reverse_qlen = result.qe + 1;
    let reverse_tlen = result.te + 1;
    let reverse = ALIGN_SCRATCH.with_borrow_mut(|scratch| {
        scratch.reverse_query.clear();
        scratch
            .reverse_query
            .extend(query[..reverse_qlen as usize].iter().rev());
        scratch.reverse_target.clear();
        scratch
            .reverse_target
            .extend(target[..reverse_tlen as usize].iter().rev());
        let reverse_query = scratch.reverse_query.as_slice();
        let slot = if size == 1 {
            &mut scratch.reverse_u8
        } else {
            &mut scratch.reverse_i16
        };
        let profile = scratch_profile(slot, size, reverse_query, m, mat);
        run_profile(
            profile,
            &scratch.reverse_target,
            o_del,
            e_del,
            o_ins,
            e_ins,
            KSW_XSTOP | result.score,
        )
    });
    if result.score == reverse.score {
        result.tb = result.te - reverse.te;
        result.qb = result.qe - reverse.qe;
    }
    result
}

/// `NOMICS_DUMP_KSW=<dir>`: append every `extend2` / `align2` call (inputs and
/// result) to `<dir>/extend2.bin` / `<dir>/align2.bin`, for offline kernel
/// benchmarks. Little-endian records:
/// - extend2: `u32 qlen, query, u32 tlen, target, i32 m, mat (m*m i8),
///   i32 o_del e_del o_ins e_ins w end_bonus zdrop h0,
///   i32 score qle tle gtle gscore max_off`
/// - align2: `u32 qlen, query, u32 tlen, target, i32 m, mat (m*m i8),
///   i32 o_del e_del o_ins e_ins xtra size(1|2),
///   i32 score te qe score2 te2 tb qb`
mod dump {
    use std::fs::File;
    use std::io::Write as _;
    use std::sync::{Mutex, OnceLock};

    pub const EXTEND2: usize = 0;
    pub const ALIGN2: usize = 1;

    pub fn sink(which: usize) -> Option<&'static Mutex<File>> {
        static SINKS: OnceLock<Option<[Mutex<File>; 2]>> = OnceLock::new();
        SINKS
            .get_or_init(|| {
                let dir = std::path::PathBuf::from(std::env::var_os("NOMICS_DUMP_KSW")?);
                let open =
                    |n: &str| Mutex::new(File::create(dir.join(n)).expect("create ksw dump"));
                Some([open("extend2.bin"), open("align2.bin")])
            })
            .as_ref()
            .map(|s| &s[which])
    }

    pub fn seq(b: &mut Vec<u8>, s: &[u8]) {
        b.extend_from_slice(&(s.len() as u32).to_le_bytes());
        b.extend_from_slice(s);
    }

    pub fn ints(b: &mut Vec<u8>, v: &[i32]) {
        for x in v {
            b.extend_from_slice(&x.to_le_bytes());
        }
    }

    pub fn write(sink: &Mutex<File>, b: &[u8]) {
        sink.lock().unwrap().write_all(b).expect("write ksw dump");
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    struct Rng(u64);
    impl Rng {
        fn below(&mut self, n: u32) -> i32 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            ((self.0 as u32) % n) as i32
        }
    }

    /// `global2_score` is `global2`'s score, on read-like pairs with equal
    /// and unequal gap penalties and bands from the minimum to wider than
    /// the query.
    /// An empty query with a band narrower than the target leaves rows with
    /// beg > end; ksw.c skips them (review finding 2026-10-21). Expected
    /// values are the pre-optimisation port's output.
    #[test]
    fn global2_empty_band_rows() {
        let mat: Vec<i8> = (0..25).map(|k| if k % 6 == 0 { 1 } else { -4 }).collect();
        let (score, cigar) = global2(0, &[], 2, &[0, 0], 5, &mat, 6, 1, 6, 1, 0);
        assert_eq!((score, cigar), (-1_073_741_824, vec![(2 << 4) | 2]));
        assert_eq!(global2_score(0, &[], 2, &[0, 0], 5, &mat, 6, 1, 6, 1, 0), score);
        let (score, cigar) = global2(1, &[0], 4, &[0, 0, 0, 0], 5, &mat, 6, 1, 6, 1, 0);
        assert_eq!(global2_score(1, &[0], 4, &[0, 0, 0, 0], 5, &mat, 6, 1, 6, 1, 0), score);
        assert!(!cigar.is_empty());
    }

    #[test]
    fn global2_score_matches_global2() {
        let mat: Vec<i8> = (0..25)
            .map(|k| match (k / 5, k % 5) {
                (4, _) | (_, 4) => -1,
                (i, j) if i == j => 1,
                _ => -4,
            })
            .collect();
        let mut rng = Rng(0x0123_4567_89ab_cdef);
        for n in 0..20_000 {
            let qlen = 1 + rng.below(200) as usize;
            let query: Vec<u8> = (0..qlen).map(|_| rng.below(5) as u8).collect();
            let mut target = query.clone();
            for _ in 0..rng.below(8) {
                let p = rng.below(target.len() as u32) as usize;
                match rng.below(3) {
                    0 => target[p] = rng.below(5) as u8,
                    1 => target.insert(p, rng.below(4) as u8),
                    _ if target.len() > 1 => {
                        target.remove(p);
                    }
                    _ => {}
                }
            }
            if n % 3 == 0 {
                target.extend((0..rng.below(60)).map(|_| rng.below(4) as u8));
            }
            let (o_del, e_del, o_ins, e_ins) = if n % 2 == 0 {
                (6, 1, 6, 1)
            } else {
                (
                    rng.below(12),
                    1 + rng.below(4),
                    rng.below(12),
                    1 + rng.below(4),
                )
            };
            let (ql, tl) = (qlen as i32, target.len() as i32);
            // bwa_gen_cigar2 always bands at least |tlen - qlen| + 3, which
            // global2's backtrack needs.
            let w = [0, 1, 3, 20, 100, 400][n % 6].max((tl - ql).abs() + 3);
            let full = global2(
                ql, &query, tl, &target, 5, &mat, o_del, e_del, o_ins, e_ins, w,
            );
            let score = global2_score(
                ql, &query, tl, &target, 5, &mat, o_del, e_del, o_ins, e_ins, w,
            );
            assert_eq!(score, full.0, "case {n}: w={w} q={query:?} t={target:?}");
        }
    }

    /// The fast `extend2` path against the literal `extend2_reference`, on
    /// read-like pairs (substitutions, indels, N) with equal and unequal gap
    /// penalties, narrow and wide bands, and z-drop on and off.
    #[test]
    fn extend2_matches_reference() {
        let mat: Vec<i8> = (0..25)
            .map(|k| match (k / 5, k % 5) {
                (4, _) | (_, 4) => -1,
                (i, j) if i == j => 1,
                _ => -4,
            })
            .collect();
        let mut rng = Rng(0x5eed_1234_abcd_0001);
        for n in 0..20_000 {
            let qlen = 1 + rng.below(160) as usize;
            let query: Vec<u8> = (0..qlen).map(|_| rng.below(5) as u8).collect();
            let mut target = query.clone();
            for _ in 0..rng.below(6) {
                let p = rng.below(target.len() as u32) as usize;
                match rng.below(3) {
                    0 => target[p] = rng.below(5) as u8,
                    1 => target.insert(p, rng.below(4) as u8),
                    _ if target.len() > 1 => {
                        target.remove(p);
                    }
                    _ => {}
                }
            }
            target.extend((0..rng.below(120)).map(|_| rng.below(4) as u8));
            let (o_del, e_del, o_ins, e_ins) = if n % 2 == 0 {
                (6, 1, 6, 1)
            } else {
                (
                    rng.below(12),
                    1 + rng.below(4),
                    rng.below(12),
                    1 + rng.below(4),
                )
            };
            let w = [0, 1, 5, 30, 100][n % 5];
            let zdrop = [0, 5, 100][n % 3];
            let h0 = 1 + rng.below(200);
            let end_bonus = rng.below(12);
            let args = (
                qlen as i32,
                target.len() as i32,
                o_del,
                e_del,
                o_ins,
                e_ins,
                w,
                end_bonus,
                zdrop,
                h0,
            );
            let fast = extend2(
                args.0, &query, args.1, &target, 5, &mat, o_del, e_del, o_ins, e_ins, w, end_bonus,
                zdrop, h0,
            );
            let reference = extend2_reference(
                args.0, &query, args.1, &target, 5, &mat, o_del, e_del, o_ins, e_ins, w, end_bonus,
                zdrop, h0,
            );
            assert_eq!(
                fast, reference,
                "case {n}: {args:?} q={query:?} t={target:?}"
            );
        }
    }
}
