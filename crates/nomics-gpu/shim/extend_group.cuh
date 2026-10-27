// BWA 0.7.19 ksw_extend2 for one group of W lanes (8, 16 or 32), shared by
// extend.cu (one call per job) and fused.cu (a seed's left and right band
// retries per job). The algorithm notes are at the top of extend.cu.
#pragma once

#include "nomics_gpu.h"

#include <cuda_runtime.h>

#include <limits.h>

#define EXT_MAX_C 8  // 32 lanes * 8 columns = 256 query bases

// Keeping C static scalarizes every array below. In particular, no dynamic
// indexing of h/e/profile is allowed: it would put the DP in local memory.
//
// Arguments follow ksw::extend2 / BWA's ksw_extend2. `query` and `target`
// point at this call's first base (callers add their offsets), `mat` at an
// m*m matrix. Requires 0 < qlen <= W * C. Every lane of the group returns the
// same result: the row loop and every scalar it updates are group-uniform
// (reduced, or shuffled from a uniform source lane). Callers do the store.
template <int W, int C>
__device__ __forceinline__ nomics_extend_result extend_group(
    const uint8_t* query, const uint8_t* target, const int qlen, const int tlen,
    const int m, const int8_t* mat, const int o_del, const int e_del,
    const int o_ins, const int e_ins, const int band_w, const int end_bonus,
    const int zdrop, const int h0) {
    const int lane = threadIdx.x & (W - 1);
    const unsigned mask = (0xffffffffu >> (32 - W)) << ((threadIdx.x & 31) & ~(W - 1));
    const int j0 = lane * C;
    const int oe_del = o_del + e_del, oe_ins = o_ins + e_ins;

    // Two words per column hold the profile for an alphabet of up to eight
    // symbols. Extract as int8_t, not uint8_t: negative substitution scores
    // and ambiguous bases must have exactly the CPU's sign extension.
    int h[C], e[C];
    unsigned profile_lo[C], profile_hi[C];
#pragma unroll
    for (int c = 0; c < C; ++c) {
        const int j = j0 + c;
        int initial = 0;
        if (j == 0) initial = h0;
        else if (h0 > oe_ins) initial = max(h0 - oe_ins - (j - 1) * e_ins, 0);
        h[c] = initial;
        e[c] = 0;
        unsigned lo = 0, hi = 0;
        if (j < qlen) {
            const int base = query[j];
#pragma unroll
            for (int k = 0; k < 8; ++k) {
                if (k < m) {
                    const unsigned score = (uint8_t)mat[k * m + base];
                    if (k < 4) lo |= score << (8 * k);
                    else hi |= score << (8 * (k - 4));
                }
            }
        }
        profile_lo[c] = lo;
        profile_hi[c] = hi;
    }
    int max_score = 0;
    for (int k = lane; k < m * m; k += W) max_score = max(max_score, (int)mat[k]);
    max_score = __reduce_max_sync(mask, max_score);
    int w = band_w;
    const int max_ins = max((int)((double)(qlen * max_score + end_bonus - o_ins) / e_ins + 1.0), 1);
    const int max_del = max((int)((double)(qlen * max_score + end_bonus - o_del) / e_del + 1.0), 1);
    w = min(w, min(max_ins, max_del));

    int mx = h0, max_i = -1, max_j = -1, max_ie = -1, gscore = -1, max_off = 0;
    int beg = 0, end = qlen;
    for (int i = 0; i < tlen; ++i) {
        beg = max(beg, i - w);
        end = min(end, min(i + w + 1, qlen));
        const int h_start = beg == 0 ? max(h0 - (o_del + e_del * (i + 1)), 0) : 0;
        const int symbol = target[i];
        const int shift = (symbol & 3) * 8;
        int mh[C], gg[C];
        const int seed = beg * e_ins;
        int local_g = INT_MIN;
#pragma unroll
        for (int c = 0; c < C; ++c) {
            const int j = j0 + c;
            if (j >= beg && j < end) {
                const unsigned profile = symbol < 4 ? profile_lo[c] : profile_hi[c];
                const int score = (int)(int8_t)(profile >> shift);
                const int diagonal = h[c] ? h[c] + score : 0;
                mh[c] = max(diagonal, e[c]);
                gg[c] = max(diagonal - oe_ins, 0) + (j + 1) * e_ins;
                local_g = max(local_g, gg[c]);
                e[c] = max(e[c] - e_del, max(diagonal - oe_del, 0));
            }
        }
        int incl = local_g;
#pragma unroll
        for (int offset = 1; offset < W; offset <<= 1) {
            const int up = __shfl_up_sync(mask, incl, offset, W);
            if (lane >= offset) incl = max(incl, up);
        }
        const int excl = __shfl_up_sync(mask, incl, 1, W);
        int running = lane == 0 ? seed : max(seed, excl);
        int row_max = 0, mj = -1, last_h = 0;
        int first = INT_MAX, last_any = -1;
#pragma unroll
        for (int c = 0; c < C; ++c) {
            const int j = j0 + c;
            if (j >= beg && j < end) {
                const int value = max(mh[c], running - j * e_ins);
                if (value >= row_max) row_max = value, mj = j;
                running = max(running, gg[c]);
                // c == 0 needs the preceding lane's last H, obtained below.
                // For other cells, shift and trim while H is already live.
                if (c != 0) {
                    h[c] = j == beg ? h_start : last_h;
                    if (h[c] != 0 || e[c] != 0) {
                        first = min(first, j);
                        last_any = j;
                    }
                }
                last_h = value;
            } else if (j == end) {
                e[c] = 0;
                if (c != 0) h[c] = beg == end ? h_start : last_h;
            }
        }
        const int prev_h = __shfl_up_sync(mask, last_h, 1, W);
        if (j0 >= beg && j0 <= end) {
            h[0] = j0 == beg ? h_start : prev_h;
            if (j0 < end && (h[0] != 0 || e[0] != 0)) {
                first = min(first, j0);
                last_any = max(last_any, j0);
            }
        }
        // eh[qlen] need not occupy a register: it is never a diagonal or E
        // input, and is overwritten before every gscore / trim observation.
        // Clamp the source for an empty band; all group lanes still shuffle.
        const int end_h = __shfl_sync(mask, last_h, max(end - 1, 0) / C, W);
        const int h1 = beg == end ? h_start : end_h;
        if (end == qlen) {
            max_ie = gscore > h1 ? max_ie : i;
            gscore = max(gscore, h1);
        }
        // Separate score/column reductions avoid the old score<<10 packing
        // limit. Ties still select the last column, even when its score is 0.
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
        // Leading scan excludes end; trailing scan includes it. Preserve
        // the CPU's all-zero fallback beg - 1 and its extra end + 1 cell.
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
