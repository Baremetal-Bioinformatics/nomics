// A seed's whole BWA extension in one GPU job (chain2aln's left band-retry loop,
// then its right one): a port of shim/fused.cu and its body extend_group_v3
// (shim/extend_v3.cuh). The algorithm, its exactness argument and the layout
// are those files'; read their comments first. Only the primitives change:
//
//  - A job is a group of W lanes (8, 16 or 32) inside a 32-wide simdgroup.
//    CUDA's width-W shuffles and reductions become group-scoped helpers below:
//      __shfl_up_sync(mask, x, d, W) -> simd_shuffle_up, keeping x for group
//                                       lanes below d (CUDA's defined result)
//      __shfl_sync(mask, x, s, W)    -> simd_shuffle from the group's lane s
//      __reduce_{max,min}_sync       -> an XOR butterfly over offsets W/2..1,
//                                       which never leaves an aligned group
//    Every helper reads only lanes of the caller's own group, so a group that
//    finished early (inactive lanes) or runs a different band never feeds its
//    neighbours. The launcher checks threadExecutionWidth == 32.
//  - max_ins / max_del. CUDA computes (int)((double)N / e + 1.0) with int N and
//    e >= 1; Metal has no double. With |N| < 2^31 the quotient N / e is either
//    an integer or at least 1/e from one, far above double rounding, so adding
//    1.0 and truncating toward zero equals C's integer (N + e) / e, which also
//    truncates toward zero. The launcher's bounds (scores, gaps, bonuses
//    <= 2^24 / 2^16) keep N + e in int. Then max(.., 1) as before.
//  - __maxnreg__ hints are dropped; Metal allocates registers itself.
#include <metal_stdlib>
using namespace metal;

#define FUSED_MAX_BAND_TRY 2
#define INT_MIN_ (-2147483647 - 1)
#define INT_MAX_ 2147483647

// nomics_extend_result / nomics_fused_params / nomics_fused_job /
// nomics_fused_result (nomics_gpu.h).
struct ExtendResult {
    int score, qle, tle, gtle, gscore, max_off;
};
struct FusedParams {
    int m, o_del, e_del, o_ins, e_ins, w, max_band_try, end_bonus_left, end_bonus_right, zdrop;
    char mat[64];
};
struct FusedJob {
    ulong lq_off, lt_off, rq_off, rt_off;
    uint lq_len, lt_len, rq_len, rt_len;
    int h0, pad;
};
struct FusedResult {
    ExtendResult left[2], right[2];
    int n_left, n_right, right_h0, pad;
};
static_assert(sizeof(ExtendResult) == 24, "C ABI layout changed");
static_assert(sizeof(FusedParams) == 104, "C ABI layout changed");
static_assert(sizeof(FusedJob) == 56, "C ABI layout changed");
static_assert(sizeof(FusedResult) == 112, "C ABI layout changed");

// ---- group-scoped collectives (W lanes, aligned within the simdgroup) ----

template <int W>
static inline int group_shfl_up(int x, ushort d, int lane) {
    const int up = simd_shuffle_up(x, d);
    return lane >= (int)d ? up : x;
}

template <int W>
static inline int group_shfl(int x, int src, uint simd_lane) {
    return simd_shuffle(x, (ushort)((simd_lane & ~(uint)(W - 1)) + (uint)src));
}

template <int W>
static inline int group_max(int x) {
    for (ushort off = W / 2; off > 0; off >>= 1) x = max(x, simd_shuffle_xor(x, off));
    return x;
}

template <int W>
static inline int group_min(int x) {
    for (ushort off = W / 2; off > 0; off >>= 1) x = min(x, simd_shuffle_xor(x, off));
    return x;
}

// ---- extend_group_v3 (shim/extend_v3.cuh), line for line ----

template <int W, int C>
static ExtendResult extend_group_v3(device const uchar* query, device const uchar* target,
                                    const int qlen, const int tlen, const int m,
                                    threadgroup const char* mat, const int o_del, const int e_del,
                                    const int o_ins, const int e_ins, const int band_w,
                                    const int end_bonus, const int zdrop, const int h0,
                                    const uint simd_lane) {
    const int lane = (int)(simd_lane & (W - 1));
    const int j0 = lane * C;
    const int oe_del = o_del + e_del, oe_ins = o_ins + e_ins;

    int h[C], e[C];
    uint profile[C];
    for (int c = 0; c < C; ++c) {
        const int j = j0 + c;
        int initial = 0;
        if (j == 0) initial = h0;
        else if (h0 > oe_ins) initial = max(h0 - oe_ins - (j - 1) * e_ins, 0);
        h[c] = initial;
        e[c] = 0;
        uint p = 0;
        if (j < qlen) {
            const int base = query[j];
            if (base < m) {
                for (int k = 0; k < 5; ++k)
                    if (k < m) p |= ((uint)(uchar)mat[k * m + base] & 63u) << (6 * k);
            }
        }
        profile[c] = p;
    }
    int max_score = 0;
    for (int k = lane; k < m * m; k += W) max_score = max(max_score, (int)mat[k]);
    max_score = group_max<W>(max_score);
    int w = band_w;
    // (int)((double)N / e + 1.0) == (N + e) / e (see the header comment).
    const int n_ins = qlen * max_score + end_bonus - o_ins;
    const int n_del = qlen * max_score + end_bonus - o_del;
    const int max_ins = max((n_ins + e_ins) / e_ins, 1);
    const int max_del = max((n_del + e_del) / e_del, 1);
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
        int incl = lane_in ? g + jeC : INT_MIN_;
        for (int offset = 1; offset < W; offset <<= 1) {
            const int up = group_shfl_up<W>(incl, (ushort)offset, lane);
            if (lane >= offset) incl = max(incl, up);
        }
        const int excl = group_shfl_up<W>(incl, 1, lane);
        int f = (lane == 0 ? seed : max(seed, excl)) - je0;
        int row_max = 0, mj = -1, last_h = 0;
        int first = INT_MAX_, last_any = -1;
        for (int c = 0; c < C; ++c) {
            const int j = j0 + c;
            const bool in = j >= beg && j < end;
            const int value = max(mh[c], f);
            f = in ? max(f - e_ins, t[c]) : f - e_ins;
            const bool better = in && value >= row_max;
            row_max = better ? value : row_max;
            mj = better ? j : mj;
            if (c != 0) {
                h[c] = in || j == end ? (j == beg ? h_start : last_h) : h[c];
                const bool live = in && (h[c] != 0 || e[c] != 0);
                first = live ? min(first, j) : first;
                last_any = live ? j : last_any;
            }
            last_h = in ? value : last_h;
            e[c] = j == end ? 0 : e[c];
        }
        const int prev_h = group_shfl_up<W>(last_h, 1, lane);
        if (j0 >= beg && j0 <= end) {
            h[0] = j0 == beg ? h_start : prev_h;
            if (j0 < end && (h[0] != 0 || e[0] != 0)) {
                first = min(first, j0);
                last_any = max(last_any, j0);
            }
        }
        const int end_h = group_shfl<W>(last_h, max(end - 1, 0) / C, simd_lane);
        const int h1 = beg == end ? h_start : end_h;
        if (end == qlen) {
            max_ie = gscore > h1 ? max_ie : i;
            gscore = max(gscore, h1);
        }
        const int best = group_max<W>(row_max);
        mj = group_max<W>(row_max == best ? mj : -1);
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
        first = group_min<W>(first);
        const int nb = first == INT_MAX_ ? end : first;
        int last = group_max<W>(last_any);
        if (h1 != 0) last = end;
        if (last < nb) last = nb - 1;
        beg = nb;
        end = min(last + 2, qlen);
    }
    ExtendResult r;
    r.score = mx;
    r.qle = max_j + 1;
    r.tle = max_i + 1;
    r.gtle = max_ie + 1;
    r.gscore = gscore;
    r.max_off = max_off;
    return r;
}

// ---- fused_group / fused_bin_kernel (shim/fused.cu) ----

template <int W, int C>
static void fused_group(device const uchar* queries, device const uchar* targets,
                        threadgroup const char* mat, constant FusedParams& p,
                        device const FusedJob& job, device FusedResult* out, uint simd_lane) {
    const int lane = (int)(simd_lane & (W - 1));
    const int lq = (int)job.lq_len, rq = (int)job.rq_len;
    int side = lq > 0 ? 0 : (rq > 0 ? 1 : 2);
    int band = 0, prev = -1, sc0 = job.h0, n_left = 0, n_right = 0;
    if (side == 1) prev = sc0;
    while (side < 2) {
        const bool right = side == 1;
        const int w = p.w << band;
        const ExtendResult x = extend_group_v3<W, C>(
            queries + (right ? job.rq_off : job.lq_off), targets + (right ? job.rt_off : job.lt_off),
            right ? rq : lq, right ? (int)job.rt_len : (int)job.lt_len, p.m, mat, p.o_del,
            p.e_del, p.o_ins, p.e_ins, w, right ? p.end_bonus_right : p.end_bonus_left, p.zdrop,
            right ? sc0 : job.h0, simd_lane);
        if (lane == 0) {
            if (right) out->right[band] = x;
            else out->left[band] = x;
        }
        const bool stop = x.score == prev || x.max_off < (w >> 1) + (w >> 2);
        if (!stop && band + 1 < p.max_band_try) {
            prev = x.score;
            ++band;
            continue;
        }
        if (!right) {
            n_left = band + 1;
            sc0 = x.score;
            prev = sc0;
            band = 0;
            side = rq > 0 ? 1 : 2;
        } else {
            n_right = band + 1;
            side = 2;
        }
    }
    if (lane == 0) {
        ExtendResult none = {0, 0, 0, 0, 0, 0};
        for (int k = n_left; k < FUSED_MAX_BAND_TRY; ++k) out->left[k] = none;
        for (int k = n_right; k < FUSED_MAX_BAND_TRY; ++k) out->right[k] = none;
        out->n_left = n_left;
        out->n_right = n_right;
        out->right_h0 = sc0;
        out->pad = 0;
    }
}

// One bin: group g runs job order[g]; every job of the bin has
// 0 < max(lq_len, rq_len) <= W * C (checked on the host).
template <int W, int C>
kernel void fused_bin(device const uchar* queries [[buffer(0)]],
                      device const uchar* targets [[buffer(1)]],
                      constant FusedParams& p [[buffer(2)]],
                      device const FusedJob* jobs [[buffer(3)]],
                      device const uint* order [[buffer(4)]],
                      constant uint& count [[buffer(5)]],
                      device FusedResult* out [[buffer(6)]],
                      threadgroup char* mat [[threadgroup(0)]],
                      uint tid [[thread_position_in_grid]],
                      uint lid [[thread_index_in_threadgroup]],
                      uint width [[threads_per_threadgroup]],
                      uint simd_lane [[thread_index_in_simdgroup]]) {
    // Before any return: every thread of the threadgroup reaches the barrier.
    for (uint k = lid; k < 64; k += width) mat[k] = p.mat[k];
    threadgroup_barrier(mem_flags::mem_threadgroup);
    const uint g = tid / W;
    if (g >= count) return;  // whole groups: the grid is a multiple of W
    const uint id = order[g];
    fused_group<W, C>(queries, targets, mat, p, jobs[id], out + id, simd_lane);
}

#define FUSED_KERNEL(W, C)                                                                   \
    template [[host_name("fused_" #W "x" #C)]] kernel void fused_bin<W, C>(                  \
        device const uchar*, device const uchar*, constant FusedParams&, device const FusedJob*, \
        device const uint*, constant uint&, device FusedResult*, threadgroup char*, uint, uint, \
        uint, uint);
#define FUSED_W(W) \
    FUSED_KERNEL(W, 1) FUSED_KERNEL(W, 2) FUSED_KERNEL(W, 3) FUSED_KERNEL(W, 4) \
    FUSED_KERNEL(W, 5) FUSED_KERNEL(W, 6) FUSED_KERNEL(W, 7) FUSED_KERNEL(W, 8)
FUSED_W(8)
FUSED_W(16)
FUSED_W(32)
