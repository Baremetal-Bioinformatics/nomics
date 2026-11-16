// BWA 0.7.19 ksw_extend2 for one group of W lanes: kernel v3's body, used by
// fused.cu. It computes exactly what extend_group (extend_group.cuh, the body
// behind nomics_ksw_extend2, which stays as it is) computes, with the same
// lane/column layout, scan, shuffles and band trimming, and fewer
// instructions and registers per cell:
//
//  1. One profile word per column. Five symbols times a 6-bit signed field
//     (bits 6k..6k+5 hold mat[k*m + base]) replace extend_group's two words of
//     8-bit fields. Requires m <= 5 and every matrix entry in [-32, 31]; the
//     launcher checks both (fused.cu) and refuses the launch otherwise.
//     Symbols are 0..m-1 as in BWA, where a larger one is an out-of-bounds
//     read; here a query symbol >= m and a target symbol above 4 score 0
//     (memory-safe, but outside the exactness claim). The launcher also bounds
//     h0 (<= 2^24), gap penalties and end bonuses (<= 2^16) and w (<= 2^20),
//     so no score, F term or band bound can overflow int.
//  2. Branch-free columns: out-of-band columns compute and discard (selects),
//     so the compiler emits no per-column branch and reconvergence.
//  3. F without multiplications. extend_group keeps running = prefix max of
//     gg_k = max(M_k - oe_ins, 0) + (k+1)*e_ins and reads F_j = running - j*e_ins.
//     Here f = running - j*e_ins is carried directly: f_{j+1} = max(f_j - e_ins,
//     t_j) for an in-band column (t_j = max(M_j - oe_ins, 0)), f_j - e_ins
//     otherwise (running unchanged). Integer identities, no rounding. The
//     lane's scan input max over in-band k of gg_k is formed the same way:
//     g = max over in-band k of t_k - (j0 + C - 1 - k)*e_ins, local_g = g +
//     (j0 + C)*e_ins; INT_MIN without an in-band column, as before. g starts
//     at -2^30; the first in-band column replaces it (t >= 0), which needs
//     C*e_ins < 2^30 (launcher: gap penalties in [0, 2^16]).
//
// Tried and dropped: one reduction for the row's best score and its column
// (key = (row_max << 9) | (mj + 1), guarded by h0 + qlen * max_score < 2^22)
// measured 1% slower than the two reductions on the GB10.
//
// Everything else is extend_group line for line; see its comments and the
// notes at the top of extend.cu.
#pragma once

#include "nomics_gpu.h"

#include <cuda_runtime.h>

#include <limits.h>

template <int W, int C>
__device__ __forceinline__ nomics_extend_result extend_group_v3(
    const uint8_t* query, const uint8_t* target, const int qlen, const int tlen,
    const int m, const int8_t* mat, const int o_del, const int e_del,
    const int o_ins, const int e_ins, const int band_w, const int end_bonus,
    const int zdrop, const int h0) {
    const int lane = threadIdx.x & (W - 1);
    const unsigned mask = (0xffffffffu >> (32 - W)) << ((threadIdx.x & 31) & ~(W - 1));
    const int j0 = lane * C;
    const int oe_del = o_del + e_del, oe_ins = o_ins + e_ins;

    int h[C], e[C];
    unsigned profile[C];
#pragma unroll
    for (int c = 0; c < C; ++c) {
        const int j = j0 + c;
        int initial = 0;
        if (j == 0) initial = h0;
        else if (h0 > oe_ins) initial = max(h0 - oe_ins - (j - 1) * e_ins, 0);
        h[c] = initial;
        e[c] = 0;
        unsigned p = 0;
        if (j < qlen) {
            // Query symbols are 0..m-1 (BWA's contract; nt4 with m = 5). A larger
            // one gets an all-zero profile instead of reading past the matrix:
            // memory-safe, but outside the exactness claim, as in BWA itself.
            const int base = query[j];
            if (base < m) {
#pragma unroll
                for (int k = 0; k < 5; ++k)
                    if (k < m) p |= ((unsigned)mat[k * m + base] & 63u) << (6 * k);
            }
        }
        profile[c] = p;
    }
    int max_score = 0;
    for (int k = lane; k < m * m; k += W) max_score = max(max_score, (int)mat[k]);
    max_score = __reduce_max_sync(mask, max_score);
    int w = band_w;
    const int max_ins = max((int)((double)(qlen * max_score + end_bonus - o_ins) / e_ins + 1.0), 1);
    const int max_del = max((int)((double)(qlen * max_score + end_bonus - o_del) / e_del + 1.0), 1);
    w = min(w, min(max_ins, max_del));
    const int je0 = j0 * e_ins, jeC = (j0 + C) * e_ins;

    int mx = h0, max_i = -1, max_j = -1, max_ie = -1, gscore = -1, max_off = 0;
    int beg = 0, end = qlen;
    for (int i = 0; i < tlen; ++i) {
        beg = max(beg, i - w);
        end = min(end, min(i + w + 1, qlen));
        const int h_start = beg == 0 ? max(h0 - (o_del + e_del * (i + 1)), 0) : 0;
        const int shift = 6 * min((int)target[i], 5);
        int mh[C], t[C];
        int g = -(1 << 30);
#pragma unroll
        for (int c = 0; c < C; ++c) {
            const int j = j0 + c;
            const bool in = j >= beg && j < end;
            const int score = (int)((profile[c] >> shift) << 26) >> 26;
            const int diagonal = h[c] ? h[c] + score : 0;
            mh[c] = max(diagonal, e[c]);
            t[c] = max(diagonal - oe_ins, 0);
            const int e_next = max(e[c] - e_del, max(diagonal - oe_del, 0));
            e[c] = in ? e_next : e[c];
            g = in ? max(g - e_ins, t[c]) : g - e_ins;
        }
        const bool lane_in = beg < end && j0 < end && j0 + C > beg;
        const int seed = beg * e_ins;
        int incl = lane_in ? g + jeC : INT_MIN;
#pragma unroll
        for (int offset = 1; offset < W; offset <<= 1) {
            const int up = __shfl_up_sync(mask, incl, offset, W);
            if (lane >= offset) incl = max(incl, up);
        }
        const int excl = __shfl_up_sync(mask, incl, 1, W);
        int f = (lane == 0 ? seed : max(seed, excl)) - je0;
        int row_max = 0, mj = -1, last_h = 0;
        int first = INT_MAX, last_any = -1;
#pragma unroll
        for (int c = 0; c < C; ++c) {
            const int j = j0 + c;
            const bool in = j >= beg && j < end;
            const int value = max(mh[c], f);
            f = in ? max(f - e_ins, t[c]) : f - e_ins;
            const bool better = in && value >= row_max;
            row_max = better ? value : row_max;
            mj = better ? j : mj;
            // c == 0 needs the preceding lane's last H, obtained below. At
            // j == end the shifted H is (beg == end ? h_start : last_h),
            // which is the in-band rule since j == beg exactly when beg == end.
            if (c != 0) {
                h[c] = in || j == end ? (j == beg ? h_start : last_h) : h[c];
                const bool live = in && (h[c] != 0 || e[c] != 0);
                first = live ? min(first, j) : first;
                last_any = live ? j : last_any;
            }
            last_h = in ? value : last_h;
            e[c] = j == end ? 0 : e[c];
        }
        const int prev_h = __shfl_up_sync(mask, last_h, 1, W);
        if (j0 >= beg && j0 <= end) {
            h[0] = j0 == beg ? h_start : prev_h;
            if (j0 < end && (h[0] != 0 || e[0] != 0)) {
                first = min(first, j0);
                last_any = max(last_any, j0);
            }
        }
        const int end_h = __shfl_sync(mask, last_h, max(end - 1, 0) / C, W);
        const int h1 = beg == end ? h_start : end_h;
        if (end == qlen) {
            max_ie = gscore > h1 ? max_ie : i;
            gscore = max(gscore, h1);
        }
        const int best = __reduce_max_sync(mask, row_max);
        mj = __reduce_max_sync(mask, row_max == best ? mj : -1);
        row_max = best;
        if (row_max == 0) break;
        if (row_max > mx) {
            mx = row_max, max_i = i, max_j = mj;
            max_off = max(max_off, abs(mj - i));
        } else if (zdrop > 0) {
            const int di = i - max_i, dj = mj - max_j;
            if (di > dj) {
                if (mx - row_max - (di - dj) * e_del > zdrop) break;
            } else if (mx - row_max - (dj - di) * e_ins > zdrop) {
                break;
            }
        }
        first = __reduce_min_sync(mask, first);
        const int nb = first == INT_MAX ? end : first;
        int last = __reduce_max_sync(mask, last_any);
        if (h1 != 0) last = end;
        if (last < nb) last = nb - 1;
        beg = nb;
        end = min(last + 2, qlen);
    }
    nomics_extend_result r;
    r.score = mx;
    r.qle = max_j + 1;
    r.tle = max_i + 1;
    r.gtle = max_ie + 1;
    r.gscore = gscore;
    r.max_off = max_off;
    return r;
}
