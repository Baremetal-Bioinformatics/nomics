// BWA 0.7.19 ksw_align2 (mate rescue) on Metal: a port of ksw_wave_kernel in
// shim/ksw.cu, one 32-wide simdgroup per alignment. The recurrences, their
// evaluation order and every tie-break are those of the CUDA kernel (see the
// derivation above ksw.cu's wave_run); only the warp primitives change:
//   __shfl_up_sync(x, 1) -> simd_shuffle_up(x, 1)
//   __shfl_sync(x, l)    -> simd_shuffle(x, l)
//   __syncwarp()         -> simdgroup_barrier(mem_threadgroup)
// The launcher checks threadExecutionWidth == 32.
#include <metal_stdlib>
using namespace metal;

#define KSW_XSTOP 0x20000
#define KSW_XSUBO 0x40000
#define KSW_XSTART 0x80000

// nomics_ksw_job / nomics_ksw_result (nomics_gpu.h).
struct ksw_job {
    ulong query_offset;
    ulong target_offset;
    ulong mat_offset;
    ulong scratch_offset;
    ulong row_offset;
    uint query_len;
    uint target_len;
    int m;
    int o_del;
    int e_del;
    int o_ins;
    int e_ins;
    int xtra;
    int size;
};
static_assert(sizeof(ksw_job) == 80, "C ABI layout changed");

struct ksw_result {
    int score;
    int te;
    int qe;
    int score2;
    int te2;
    int tb;
    int qb;
};
static_assert(sizeof(ksw_result) == 28, "C ABI layout changed");

struct wave_pass {
    int gmax;
    int te;
    int qe;
    uint processed;
};

static inline uchar query_base(device const uchar* query, uint qlen, uint k, bool reverse) {
    return reverse ? query[qlen - 1 - k] : query[k];
}

static inline uchar target_base(threadgroup const uchar* target, uint tlen, uint i, bool reverse) {
    return reverse ? target[tlen - 1 - i] : target[i];
}

template <int C>
static wave_pass wave_run(device const uchar* query, uint qlen, threadgroup const uchar* target,
                          uint tlen, device const char* mat, int m, int o_del, int e_del, int o_ins,
                          int e_ins, int xtra, int size, bool reverse, threadgroup ushort* rows,
                          bool write_rows, int lane) {
    const int lanes = size == 1 ? 16 : 8;
    const int slen = ((int)qlen + lanes - 1) / lanes;
    const int P = slen * lanes;
    int minimum = 127;
    for (int k = 0; k < m * m; ++k) minimum = minimum < (int)mat[k] ? minimum : (int)mat[k];
    const int shift = (int)(uchar)(0u - (uchar)(char)minimum);
    const int endsc = (xtra & KSW_XSTOP) ? (xtra & 0xffff) : 0x10000;
    const int oe_del = o_del + e_del, oe_ins = o_ins + e_ins;
    const int k0 = lane * C;
    const int last_lane = (P - 1) / C;

    // Per column: its scores against each target base, one byte each (m <= 8;
    // padding columns score 0), so the inner loop does no memory access.
    ulong tab[C];
    int hprev[C], e[C];
    uint stripe_start = 0;  // bit c: column k0+c starts a stripe (Fp resets)
    for (int c = 0; c < C; ++c) {
        const int k = k0 + c;
        if (k % slen == 0) stripe_start |= 1u << c;
        tab[c] = 0;
        if (k < (int)qlen) {
            const int qb = (int)query_base(query, qlen, (uint)k, reverse);
            for (int b = 0; b < m; ++b) tab[c] |= (ulong)(uchar)mat[b * m + qb] << (8 * b);
        }
        hprev[c] = 0;
        e[c] = 0;
    }
    // Values lane l sends to lane l+1 for the row it just finished.
    int send_h = 0, send_f = 0, send_fp = 0, send_max = -1, send_arg = -1;
    int diag_left = 0;  // H[r-1][k0-1] for this lane's next row
    int gmax = 0, te = -1, qe = 0;
    uint processed = tlen;
    bool stop = false;
    const int steps = (int)tlen + last_lane;
    for (int s = 0; s < steps; ++s) {
        int in_h = simd_shuffle_up(send_h, 1);
        int in_f = simd_shuffle_up(send_f, 1);
        int in_fp = simd_shuffle_up(send_fp, 1);
        int in_max = simd_shuffle_up(send_max, 1);
        int in_arg = simd_shuffle_up(send_arg, 1);
        if (lane == 0) in_h = in_f = in_fp = 0, in_max = -1, in_arg = -1;
        const int r = s - lane;
        if (r >= 0 && r < (int)tlen && lane <= last_lane) {
            const int tshift = 8 * (int)target_base(target, tlen, (uint)r, reverse);
            int diag = diag_left, f = in_f, fp = in_fp, mx = in_max, arg = in_arg, h_last = in_h;
            for (int c = 0; c < C; ++c) {
                const int k = k0 + c;
                if (k < P) {
                    if (stripe_start >> c & 1) fp = 0;
                    const int sc = (int)(char)(uchar)(tab[c] >> tshift);
                    int vp = diag + sc;
                    vp = vp > 0 ? vp : 0;
                    vp = vp > e[c] ? vp : e[c];
                    vp = vp > fp ? vp : fp;
                    const int h = vp > f ? vp : f;
                    int t1 = e[c] - e_del, t2 = vp - oe_del;
                    t1 = t1 > t2 ? t1 : t2;
                    e[c] = t1 > 0 ? t1 : 0;
                    t1 = fp - e_ins;
                    t2 = vp - oe_ins;
                    t1 = t1 > t2 ? t1 : t2;
                    fp = t1 > 0 ? t1 : 0;
                    t1 = f - e_ins;
                    t2 = h - oe_ins;
                    t1 = t1 > t2 ? t1 : t2;
                    f = t1 > 0 ? t1 : 0;
                    if (vp > mx) mx = vp, arg = k;
                    diag = hprev[c];
                    hprev[c] = h;
                    h_last = h;
                }
            }
            diag_left = in_h;
            send_h = h_last, send_f = f, send_fp = fp, send_max = mx, send_arg = arg;
            if (lane == last_lane) {
                if (write_rows) rows[r] = (ushort)(mx > 0xffff ? 0xffff : mx);
                if (mx > gmax) {
                    gmax = mx, te = r, qe = arg;
                    if ((size == 1 && gmax + shift >= 255) || gmax >= endsc) {
                        stop = true;
                        processed = (uint)r + 1;
                    }
                }
            }
        }
        if (simd_shuffle((int)stop, (ushort)last_lane)) break;
    }
    wave_pass out;
    out.gmax = simd_shuffle(gmax, (ushort)last_lane);
    out.te = simd_shuffle(te, (ushort)last_lane);
    out.qe = simd_shuffle(qe, (ushort)last_lane);
    out.processed = simd_shuffle(processed, (ushort)last_lane);
    return out;
}

static void consider_second(thread ksw_result& result, int low, int high, int score, int end) {
    if ((end < low || end > high) && score > result.score2) {
        result.score2 = score;
        result.te2 = end;
    }
}

template <int C>
kernel void ksw_wave(device const uchar* queries [[buffer(0)]],
                     device const uchar* targets [[buffer(1)]],
                     device const char* matrices [[buffer(2)]],
                     device const ksw_job* jobs [[buffer(3)]],
                     device ksw_result* out [[buffer(4)]],
                     constant uint& n [[buffer(5)]],
                     constant uint& shared_per_warp [[buffer(6)]],
                     threadgroup uchar* smem [[threadgroup(0)]],
                     uint group [[threadgroup_position_in_grid]],
                     uint warp_in_group [[simdgroup_index_in_threadgroup]],
                     uint warps_per_group [[simdgroups_per_threadgroup]],
                     uint lane_id [[thread_index_in_simdgroup]]) {
    const uint id = group * warps_per_group + warp_in_group;
    const int lane = (int)lane_id;
    if (id >= n) return;  // whole simdgroups only
    const ksw_job job = jobs[id];
    device const uchar* query = queries + job.query_offset;
    device const char* mat = matrices + job.mat_offset;
    // The target and the row maxima live in this simdgroup's threadgroup
    // memory: the inner loop then touches no host-mapped memory.
    threadgroup uchar* target = smem + warp_in_group * shared_per_warp;
    threadgroup ushort* rows = (threadgroup ushort*)(target + ((job.target_len + 15) & ~15u));
    for (uint i = lane_id; i < job.target_len; i += 32) target[i] = targets[job.target_offset + i];
    simdgroup_barrier(mem_flags::mem_threadgroup);

    const wave_pass fwd = wave_run<C>(query, job.query_len, target, job.target_len, mat, job.m,
                                      job.o_del, job.e_del, job.o_ins, job.e_ins, job.xtra,
                                      job.size, false, rows, true, lane);
    int minimum = 127, profile_max = 0;
    for (int k = 0; k < job.m * job.m; ++k) {
        minimum = minimum < (int)mat[k] ? minimum : (int)mat[k];
        profile_max = profile_max > (int)mat[k] ? profile_max : (int)mat[k];
    }
    const int shift = (int)(uchar)(0u - (uchar)(char)minimum);
    ksw_result result = {0, -1, -1, -1, -1, -1, -1};
    result.te = fwd.te;
    if (job.size == 1 && fwd.gmax + shift >= 255) {
        result.score = 255;
    } else {
        result.score = fwd.gmax;
        result.qe = fwd.qe;
        simdgroup_barrier(mem_flags::mem_threadgroup);  // rows[] was written by another lane
        // second_best over record_row's merged row list (lane 0 replays it).
        if (lane == 0 && profile_max > 0) {
            const int minsc = (job.xtra & KSW_XSUBO) ? (job.xtra & 0xffff) : 0x10000;
            const int radius = (result.score + profile_max - 1) / profile_max;
            const int low = result.te - radius, high = result.te + radius;
            bool have = false;
            int last_score = 0, last_end = -1;
            for (uint i = 0; i < fwd.processed; ++i) {
                const int score = rows[i];
                if (score < minsc) continue;
                if (have && last_end + 1 == (int)i) {
                    if (last_score < score) last_score = score, last_end = (int)i;
                } else {
                    if (have) consider_second(result, low, high, last_score, last_end);
                    have = true, last_score = score, last_end = (int)i;
                }
            }
            if (have) consider_second(result, low, high, last_score, last_end);
        }
        result.score2 = simd_shuffle(result.score2, (ushort)0);
        result.te2 = simd_shuffle(result.te2, (ushort)0);
    }
    if ((job.xtra & KSW_XSTART) && !((job.xtra & KSW_XSUBO) && result.score < (job.xtra & 0xffff)) &&
        result.qe >= 0) {
        const wave_pass rev = wave_run<C>(query, (uint)(result.qe + 1), target, (uint)(result.te + 1),
                                          mat, job.m, job.o_del, job.e_del, job.o_ins, job.e_ins,
                                          KSW_XSTOP | result.score, job.size, true, rows, false, lane);
        if (rev.gmax == result.score) {
            result.tb = result.te - rev.te;
            result.qb = result.qe - rev.qe;
        }
    }
    if (lane == 0) out[id] = result;
}

#define KSW_WAVE(C)                                                                              \
    template [[host_name("ksw_wave_" #C)]] kernel void ksw_wave<C>(                              \
        device const uchar*, device const uchar*, device const char*, device const ksw_job*,     \
        device ksw_result*, constant uint&, constant uint&, threadgroup uchar*, uint, uint, uint, \
        uint);
KSW_WAVE(4)
KSW_WAVE(5)
KSW_WAVE(6)
KSW_WAVE(8)
KSW_WAVE(16)
