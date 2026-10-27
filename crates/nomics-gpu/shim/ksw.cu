// Experimental BWA 0.7.19 ksw_align2 kernel. One CUDA thread owns an
// alignment and evaluates the striped vectors, lanes, and lazy-F loop in the
// same order as ksw_u8/ksw_i16. This deliberately forgoes intra-alignment
// parallel reductions: their evaluation order would change endpoint ties.
#include "nomics_gpu.h"
#include "ticket.cuh"

#include <cuda_runtime.h>

#include <limits.h>

#define KSW_XBYTE 0x10000
#define KSW_XSTOP 0x20000
#define KSW_XSUBO 0x40000
#define KSW_XSTART 0x80000

static_assert(sizeof(nomics_ksw_job) == 80, "C ABI layout changed");
static_assert(sizeof(nomics_ksw_result) == 28, "C ABI layout changed");

struct pass_result {
    int32_t score;
    int32_t te;
    int32_t qe;
    uint32_t processed;
};

__device__ __forceinline__ uint8_t sub_u8(uint8_t a, uint8_t b) { return a > b ? a - b : 0; }

__device__ __forceinline__ uint8_t add_u8(uint8_t a, uint8_t b) {
    const unsigned value = (unsigned)a + (unsigned)b;
    return value > UCHAR_MAX ? UCHAR_MAX : (uint8_t)value;
}

__device__ __forceinline__ uint16_t sub_u16(uint16_t a, uint16_t b) { return a > b ? a - b : 0; }

__device__ __forceinline__ int16_t add_i16(int16_t a, int16_t b) {
    const int value = (int)a + (int)b;
    return value > SHRT_MAX ? SHRT_MAX : (value < SHRT_MIN ? SHRT_MIN : (int16_t)value);
}

__device__ __forceinline__ uint8_t query_base(const uint8_t* query, uint32_t qlen, uint32_t k,
                                               bool reverse) {
    return reverse ? query[qlen - 1 - k] : query[k];
}

__device__ __forceinline__ uint8_t target_base(const uint8_t* target, uint32_t tlen, uint32_t i,
                                                bool reverse) {
    return reverse ? target[tlen - 1 - i] : target[i];
}

__device__ pass_result run_u8(const uint8_t* query, uint32_t qlen, const uint8_t* target, uint32_t tlen,
                              int32_t m, const int8_t* mat, int32_t o_del, int32_t e_del, int32_t o_ins,
                              int32_t e_ins, int32_t xtra, bool reverse, uint16_t* workspace,
                              uint32_t workspace_stride, uint16_t* row_maxima) {
    const uint32_t lanes = 16;
    const uint32_t slen = (qlen + lanes - 1) / lanes;
    const uint32_t padded = slen * lanes;
    uint16_t* h0 = workspace;
    uint16_t* h1 = workspace + workspace_stride;
    uint16_t* e = workspace + 2 * workspace_stride;
    uint16_t* hmax = workspace + 3 * workspace_stride;
    for (uint32_t k = 0; k < padded; ++k) h0[k] = e[k] = hmax[k] = 0;

    int minimum = 127;
    int maximum = 0;
    for (int k = 0; k < m * m; ++k) {
        minimum = minimum < (int)mat[k] ? minimum : (int)mat[k];
        maximum = maximum > (int)mat[k] ? maximum : (int)mat[k];
    }
    const uint8_t shift = (uint8_t)(0u - (uint8_t)(int8_t)minimum);
    const int32_t endsc = (xtra & KSW_XSTOP) ? (xtra & 0xffff) : 0x10000;
    const uint8_t gap_del = (uint8_t)e_del;
    const uint8_t open_del = (uint8_t)((uint32_t)o_del + (uint32_t)e_del);
    const uint8_t gap_ins = (uint8_t)e_ins;
    const uint8_t open_ins = (uint8_t)((uint32_t)o_ins + (uint32_t)e_ins);

    int32_t gmax = 0;
    int32_t te = -1;
    uint32_t processed = 0;
    for (uint32_t i = 0; i < tlen; ++i) {
        const uint8_t tb = target_base(target, tlen, i, reverse);
        uint8_t f[16] = {};
        uint8_t vmax[16] = {};
        uint8_t h[16];
        h[0] = 0;
        for (uint32_t lane = 1; lane < lanes; ++lane)
            h[lane] = (uint8_t)h0[(slen - 1) * lanes + lane - 1];

        for (uint32_t j = 0; j < slen; ++j) {
            for (uint32_t lane = 0; lane < lanes; ++lane) {
                const uint32_t index = j * lanes + lane;
                const uint32_t qindex = j + lane * slen;
                const int score = qindex < qlen ? (int)mat[(int)tb * m +
                                                            (int)query_base(query, qlen, qindex, reverse)]
                                                : 0;
                const uint8_t profile = (uint8_t)((uint8_t)(int8_t)score + shift);
                const uint8_t ev = (uint8_t)e[index];
                uint8_t value = sub_u8(add_u8(h[lane], profile), shift);
                value = value > ev ? value : ev;
                value = value > f[lane] ? value : f[lane];
                vmax[lane] = vmax[lane] > value ? vmax[lane] : value;
                const uint8_t next_e = sub_u8(ev, gap_del) > sub_u8(value, open_del)
                                               ? sub_u8(ev, gap_del)
                                               : sub_u8(value, open_del);
                const uint8_t next_f = sub_u8(f[lane], gap_ins) > sub_u8(value, open_ins)
                                               ? sub_u8(f[lane], gap_ins)
                                               : sub_u8(value, open_ins);
                const uint8_t old_h = (uint8_t)h0[index];
                h1[index] = value;
                e[index] = next_e;
                f[lane] = next_f;
                h[lane] = old_h;
            }
        }

        bool stop = false;
        for (uint32_t lazy = 0; lazy < lanes && !stop; ++lazy) {
            for (uint32_t lane = lanes - 1; lane > 0; --lane) f[lane] = f[lane - 1];
            f[0] = 0;
            for (uint32_t j = 0; j < slen; ++j) {
                bool active = false;
                for (uint32_t lane = 0; lane < lanes; ++lane) {
                    const uint32_t index = j * lanes + lane;
                    uint8_t value = (uint8_t)h1[index];
                    value = value > f[lane] ? value : f[lane];
                    h1[index] = value;
                    const uint8_t h_gap = sub_u8(value, open_ins);
                    f[lane] = sub_u8(f[lane], gap_ins);
                    active |= f[lane] > h_gap;
                }
                if (!active) {
                    stop = true;
                    break;
                }
            }
        }

        int32_t imax = 0;
        for (uint32_t lane = 0; lane < lanes; ++lane) imax = imax > vmax[lane] ? imax : vmax[lane];
        if (row_maxima) row_maxima[i] = (uint16_t)imax;
        processed = i + 1;
        if (imax > gmax) {
            gmax = imax;
            te = (int32_t)i;
            for (uint32_t k = 0; k < padded; ++k) hmax[k] = h1[k];
            if (gmax + (int32_t)shift >= 255 || gmax >= endsc) break;
        }
        uint16_t* swap = h0;
        h0 = h1;
        h1 = swap;
    }

    pass_result result = {};
    result.score = gmax + (int32_t)shift < 255 ? gmax : 255;
    result.te = te;
    result.qe = -1;
    result.processed = processed;
    if (result.score != 255) {
        int32_t best = -1;
        for (uint32_t index = 0; index < padded; ++index) {
            const int32_t value = (uint8_t)hmax[index];
            const int32_t query_end = (int32_t)(index / lanes + (index % lanes) * slen);
            if (value > best) {
                best = value;
                result.qe = query_end;
            } else if (value == best && query_end < result.qe) {
                result.qe = query_end;
            }
        }
    }
    return result;
}

__device__ pass_result run_i16(const uint8_t* query, uint32_t qlen, const uint8_t* target, uint32_t tlen,
                               int32_t m, const int8_t* mat, int32_t o_del, int32_t e_del, int32_t o_ins,
                               int32_t e_ins, int32_t xtra, bool reverse, uint16_t* workspace,
                               uint32_t workspace_stride, uint16_t* row_maxima) {
    const uint32_t lanes = 8;
    const uint32_t slen = (qlen + lanes - 1) / lanes;
    const uint32_t padded = slen * lanes;
    uint16_t* h0 = workspace;
    uint16_t* h1 = workspace + workspace_stride;
    uint16_t* e = workspace + 2 * workspace_stride;
    uint16_t* hmax = workspace + 3 * workspace_stride;
    for (uint32_t k = 0; k < padded; ++k) h0[k] = e[k] = hmax[k] = 0;

    const int32_t endsc = (xtra & KSW_XSTOP) ? (xtra & 0xffff) : 0x10000;
    const uint16_t gap_del = (uint16_t)e_del;
    const uint16_t open_del = (uint16_t)((uint32_t)o_del + (uint32_t)e_del);
    const uint16_t gap_ins = (uint16_t)e_ins;
    const uint16_t open_ins = (uint16_t)((uint32_t)o_ins + (uint32_t)e_ins);

    int32_t gmax = 0;
    int32_t te = -1;
    uint32_t processed = 0;
    for (uint32_t i = 0; i < tlen; ++i) {
        const uint8_t tb = target_base(target, tlen, i, reverse);
        int16_t f[8] = {};
        int16_t vmax[8] = {};
        int16_t h[8];
        h[0] = 0;
        for (uint32_t lane = 1; lane < lanes; ++lane)
            h[lane] = (int16_t)h0[(slen - 1) * lanes + lane - 1];

        for (uint32_t j = 0; j < slen; ++j) {
            for (uint32_t lane = 0; lane < lanes; ++lane) {
                const uint32_t index = j * lanes + lane;
                const uint32_t qindex = j + lane * slen;
                const int16_t profile = qindex < qlen
                                            ? (int16_t)mat[(int)tb * m +
                                                           (int)query_base(query, qlen, qindex, reverse)]
                                            : 0;
                const int16_t ev = (int16_t)e[index];
                int16_t value = add_i16(h[lane], profile);
                value = value > ev ? value : ev;
                value = value > f[lane] ? value : f[lane];
                vmax[lane] = vmax[lane] > value ? vmax[lane] : value;
                const int16_t next_e_a = (int16_t)sub_u16((uint16_t)ev, gap_del);
                const int16_t next_e_b = (int16_t)sub_u16((uint16_t)value, open_del);
                const int16_t next_f_a = (int16_t)sub_u16((uint16_t)f[lane], gap_ins);
                const int16_t next_f_b = (int16_t)sub_u16((uint16_t)value, open_ins);
                const int16_t old_h = (int16_t)h0[index];
                h1[index] = (uint16_t)value;
                e[index] = (uint16_t)(next_e_a > next_e_b ? next_e_a : next_e_b);
                f[lane] = next_f_a > next_f_b ? next_f_a : next_f_b;
                h[lane] = old_h;
            }
        }

        bool stop = false;
        for (uint32_t lazy = 0; lazy < 16 && !stop; ++lazy) {
            for (uint32_t lane = lanes - 1; lane > 0; --lane) f[lane] = f[lane - 1];
            f[0] = 0;
            for (uint32_t j = 0; j < slen; ++j) {
                bool active = false;
                for (uint32_t lane = 0; lane < lanes; ++lane) {
                    const uint32_t index = j * lanes + lane;
                    int16_t value = (int16_t)h1[index];
                    value = value > f[lane] ? value : f[lane];
                    h1[index] = (uint16_t)value;
                    const int16_t h_gap = (int16_t)sub_u16((uint16_t)value, open_ins);
                    f[lane] = (int16_t)sub_u16((uint16_t)f[lane], gap_ins);
                    active |= f[lane] > h_gap;
                }
                if (!active) {
                    stop = true;
                    break;
                }
            }
        }

        int32_t imax = 0;
        for (uint32_t lane = 0; lane < lanes; ++lane)
            imax = imax > (int32_t)vmax[lane] ? imax : (int32_t)vmax[lane];
        if (row_maxima) row_maxima[i] = (uint16_t)imax;
        processed = i + 1;
        if (imax > gmax) {
            gmax = imax;
            te = (int32_t)i;
            for (uint32_t k = 0; k < padded; ++k) hmax[k] = h1[k];
            if (gmax >= endsc) break;
        }
        uint16_t* swap = h0;
        h0 = h1;
        h1 = swap;
    }

    pass_result result = {};
    result.score = gmax;
    result.te = te;
    result.qe = -1;
    result.processed = processed;
    int32_t best = -1;
    for (uint32_t index = 0; index < padded; ++index) {
        const int32_t value = (uint16_t)hmax[index];
        const int32_t query_end = (int32_t)(index / lanes + (index % lanes) * slen);
        if (value > best) {
            best = value;
            result.qe = query_end;
        } else if (value == best && query_end < result.qe) {
            result.qe = query_end;
        }
    }
    return result;
}

__device__ void consider_second(nomics_ksw_result* result, int32_t low, int32_t high, int32_t score,
                                int32_t end) {
    if ((end < low || end > high) && score > result->score2) {
        result->score2 = score;
        result->te2 = end;
    }
}

__global__ void ksw_kernel(const uint8_t* queries, const uint8_t* targets, const int8_t* matrices,
                           const nomics_ksw_job* jobs, uint32_t n, uint16_t* scratch,
                           uint16_t* row_maxima, nomics_ksw_result* out) {
    const uint32_t id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id >= n) return;
    const nomics_ksw_job job = jobs[id];
    const uint8_t* query = queries + job.query_offset;
    const uint8_t* target = targets + job.target_offset;
    const int8_t* mat = matrices + job.mat_offset;
    uint16_t* workspace = scratch + job.scratch_offset;
    uint16_t* rows = row_maxima + job.row_offset;
    const uint32_t lanes = job.size == 1 ? 16 : 8;
    const uint32_t stride = ((job.query_len + lanes - 1) / lanes) * lanes;

    pass_result forward = job.size == 1
                              ? run_u8(query, job.query_len, target, job.target_len, job.m, mat, job.o_del,
                                       job.e_del, job.o_ins, job.e_ins, job.xtra, false, workspace, stride, rows)
                              : run_i16(query, job.query_len, target, job.target_len, job.m, mat, job.o_del,
                                        job.e_del, job.o_ins, job.e_ins, job.xtra, false, workspace, stride, rows);
    nomics_ksw_result result = {forward.score, forward.te, forward.qe, -1, -1, -1, -1};

    int profile_max = 0;
    for (int k = 0; k < job.m * job.m; ++k) profile_max = profile_max > (int)mat[k] ? profile_max : (int)mat[k];
    const int32_t minsc = (job.xtra & KSW_XSUBO) ? (job.xtra & 0xffff) : 0x10000;
    if (!(job.size == 1 && result.score == 255) && forward.processed && profile_max > 0) {
        const int32_t radius = (result.score + profile_max - 1) / profile_max;
        const int32_t low = result.te - radius;
        const int32_t high = result.te + radius;
        bool have = false;
        int32_t last_score = 0;
        int32_t last_end = -1;
        for (uint32_t i = 0; i < forward.processed; ++i) {
            const int32_t score = rows[i];
            if (score < minsc) continue;
            if (have && last_end + 1 == (int32_t)i) {
                if (last_score < score) {
                    last_score = score;
                    last_end = (int32_t)i;
                }
            } else {
                if (have) consider_second(&result, low, high, last_score, last_end);
                have = true;
                last_score = score;
                last_end = (int32_t)i;
            }
        }
        if (have) consider_second(&result, low, high, last_score, last_end);
    }

    if ((job.xtra & KSW_XSTART) &&
        !((job.xtra & KSW_XSUBO) && result.score < (job.xtra & 0xffff))) {
        const uint32_t reverse_qlen = (uint32_t)(result.qe + 1);
        const uint32_t reverse_tlen = (uint32_t)(result.te + 1);
        const int32_t reverse_xtra = KSW_XSTOP | result.score;
        pass_result reverse = job.size == 1
                                  ? run_u8(query, reverse_qlen, target, reverse_tlen, job.m, mat, job.o_del,
                                           job.e_del, job.o_ins, job.e_ins, reverse_xtra, true, workspace, stride,
                                           nullptr)
                                  : run_i16(query, reverse_qlen, target, reverse_tlen, job.m, mat, job.o_del,
                                            job.e_del, job.o_ins, job.e_ins, reverse_xtra, true, workspace, stride,
                                            nullptr);
        if (result.score == reverse.score) {
            result.tb = result.te - reverse.te;
            result.qb = result.qe - reverse.qe;
        }
    }
    out[id] = result;
}

extern "C" int nomics_ksw_align2_launch(nomics_stream stream, const uint8_t* queries, const uint8_t* targets,
                                         const int8_t* matrices, const nomics_ksw_job* jobs, uint32_t n,
                                         uint16_t* scratch, uint16_t* row_maxima, nomics_ksw_result* out,
                                         nomics_ticket* ticket) {
    if (ticket == nullptr) return cudaErrorInvalidValue;
    *ticket = nullptr;
    if (n == 0) return cudaSuccess;
    nomics_ticket_s* t = nullptr;
    cudaError_t error = ticket_begin((cudaStream_t)stream, &t);
    if (error != cudaSuccess) return error;
    const int threads = 128;
    ksw_kernel<<<(n + threads - 1) / threads, threads, 0, t->stream>>>(queries, targets, matrices, jobs, n, scratch,
                                                                       row_maxima, out);
    error = cudaGetLastError();
    return ticket_finish(t, error, ticket);
}

// ---------------------------------------------------------------------------
// Wavefront ksw_align2: one warp per alignment. Lane `l` owns query columns
// [l*C, l*C+C) of the padded query (BWA pads to the striped width with
// score-0 columns) and processes target row `s - l` at step `s`, receiving
// the row's left boundary from lane l-1 by shuffle. The result equals the
// striped kernel's because only its observable recurrences are reproduced:
//   Vp = max(H[i-1][k-1] + s, E, Fp, 0)  (Fp: insertions within the stripe,
//        reset at every stripe start k % slen == 0, as in the striped loop)
//   E' = max(E - e_del, Vp - oe_del, 0)   (from Vp: BWA's lazy-F pass never
//        updates E)
//   H  = max(Vp, F)                        (F: all insertions; = lazy-F result)
// Row maxima are taken over Vp (equal to those over H: a cell raised by F is
// below the cell its gap opened from). qe is the smallest column holding the
// row maximum at te, as BWA's scan of Hmax picks it.
// ---------------------------------------------------------------------------

struct wave_pass {
    int32_t gmax;
    int32_t te;
    int32_t qe;
    uint32_t processed;
};

// The score of a column against target base t: byte t of the column's
// 8-byte table (lo: bases 0-3, hi: 4-7), sign-extended, in one PRMT. `sel`
// is wave_sel(t), computed once per row. (__byte_perm masks the selector's
// sign-replicate bits, hence the inline PTX.)
__device__ __forceinline__ int32_t wave_score(uint32_t lo, uint32_t hi, uint32_t sel) {
    uint32_t d;
    asm("prmt.b32 %0, %1, %2, %3;" : "=r"(d) : "r"(lo), "r"(hi), "r"(sel));
    return (int32_t)d;
}

// A per-column constant kept in a register: without this the compiler
// recomputes it from the lane index inside the step loop.
__device__ __forceinline__ uint32_t wave_opaque(uint32_t x) {
    asm volatile("mov.b32 %0, %0;" : "+r"(x));
    return x;
}

// Shared-memory accesses by 32-bit shared address, held in a register: the
// step loops otherwise rebuild the warp's pointer from threadIdx each step.
__device__ __forceinline__ uint32_t smem_addr(const void* p) {
    return wave_opaque((uint32_t)__cvta_generic_to_shared(p));
}
__device__ __forceinline__ uint32_t lds_u8(uint32_t at) {
    uint32_t v;
    asm volatile("ld.shared.u8 %0, [%1];" : "=r"(v) : "r"(at) : "memory");
    return v;
}
__device__ __forceinline__ uint32_t lds_u16(uint32_t at) {
    uint32_t v;
    asm volatile("ld.shared.u16 %0, [%1];" : "=r"(v) : "r"(at) : "memory");
    return v;
}
__device__ __forceinline__ void sts_u16(uint32_t at, uint32_t v) {
    asm volatile("st.shared.u16 [%0], %1;" ::"r"(at), "r"(v) : "memory");
}
__device__ __forceinline__ void sts_u32(uint32_t at, uint32_t v) {
    asm volatile("st.shared.u32 [%0], %1;" ::"r"(at), "r"(v) : "memory");
}

// Byte 0 = byte t; bytes 1-3 replicate its sign (selector nibble t | 8).
__device__ __forceinline__ uint32_t wave_sel(uint32_t t) { return t * 0x1111u + 0x8880u; }

// Steps between checks of the last lane's stop flag. The last lane freezes
// its result at the stopping row, so later steps only compute rows nobody
// reads (rows[] past `processed` is never scanned).
constexpr int WAVE_STOP_EVERY = 8;

// FULL: every lane owns C real columns (P == 32 * C, so last_lane == 31 and
// no column is padding past P), the 151-bp rescue shape at C = 5; the
// per-column bound checks then drop out.
//
// The row maximum and its column travel as one key, vp * 1024 + (1023 - k):
// its maximum has the largest vp and, among equal vp, the smallest k, which
// is what the strict `vp > mx` scan in increasing k selects. vp < 2^16
// (at most 512 columns of score <= 127), so the key fits. H and F travel
// packed as two 16-bit halves, both in [0, 2^16).
template <int C, bool FULL>
__device__ wave_pass wave_run(const uint8_t* query, uint32_t qlen, const uint8_t* target, uint32_t tlen,
                              const int8_t* mat, int32_t m, int32_t o_del, int32_t e_del, int32_t o_ins,
                              int32_t e_ins, int32_t xtra, int32_t size, bool reverse, uint16_t* rows) {
    // Held in registers through the loop (else rebuilt from threadIdx each step).
    const int lane = (int)wave_opaque(threadIdx.x & 31);
    const uint32_t target_at = smem_addr(target);
    const uint32_t rows_at = rows ? smem_addr(rows) : 0u;
    const unsigned full = 0xffffffffu;
    const int32_t lanes = size == 1 ? 16 : 8;
    const int32_t slen = ((int32_t)qlen + lanes - 1) / lanes;
    const int32_t P = slen * lanes;
    int32_t minimum = 127;
    for (int k = 0; k < m * m; ++k) minimum = minimum < (int32_t)mat[k] ? minimum : (int32_t)mat[k];
    const int32_t shift = (int32_t)(uint8_t)(0u - (uint8_t)(int8_t)minimum);
    const int32_t endsc = (xtra & KSW_XSTOP) ? (xtra & 0xffff) : 0x10000;
    const int32_t oe_del = o_del + e_del, oe_ins = o_ins + e_ins;
    const int32_t k0 = lane * C;
    const int last_lane = FULL ? 31 : (P - 1) / C;

    // Per column: its scores against each target base, one byte each (m <= 8;
    // padding columns score 0), so the inner loop does no memory access.
    uint32_t tab_lo[C], tab_hi[C];
    int32_t keep[C];  // 0 where column k starts a stripe (Fp resets), else ~0
    int32_t kinv[C];  // 1023 - k, the key's column part
    int32_t hprev[C], e[C];
#pragma unroll
    for (int c = 0; c < C; ++c) {
        const int32_t k = k0 + c;
        keep[c] = (int32_t)wave_opaque(k % slen == 0 ? 0u : ~0u);
        kinv[c] = (int32_t)wave_opaque((uint32_t)(1023 - k));
        uint64_t tab = 0;
        if (k < (int32_t)qlen) {
            const int32_t qb = (int32_t)query_base(query, qlen, (uint32_t)k, reverse);
            for (int b = 0; b < m; ++b) tab |= (uint64_t)(uint8_t)mat[b * m + qb] << (8 * b);
        }
        tab_lo[c] = (uint32_t)tab;
        tab_hi[c] = (uint32_t)(tab >> 32);
        hprev[c] = 0;
        e[c] = 0;
    }
    // Values lane l sends to lane l+1 for the row it just finished.
    int32_t send_hf = 0, send_fp = 0, send_key = -1;
    int32_t diag_left = 0;  // H[r-1][k0-1] for this lane's next row
    int32_t gmax = 0, te = -1, qe = 0;
    uint32_t processed = tlen;
    bool stop = false;
    const int32_t steps = (int32_t)tlen + last_lane;
    for (int32_t s = 0; s < steps; ++s) {
        int32_t in_hf = __shfl_up_sync(full, send_hf, 1);
        int32_t in_fp = __shfl_up_sync(full, send_fp, 1);
        int32_t in_key = __shfl_up_sync(full, send_key, 1);
        if (lane == 0) in_hf = in_fp = 0, in_key = -1;
        const int32_t r = s - lane;
        if (r >= 0 && r < (int32_t)tlen && (FULL || lane <= last_lane)) {
            const uint32_t sel = wave_sel(lds_u8(target_at + (reverse ? tlen - 1 - (uint32_t)r : (uint32_t)r)));
            int32_t diag = diag_left, f = (int32_t)((uint32_t)in_hf >> 16), fp = in_fp, key = in_key;
            int32_t h_last = in_hf & 0xffff;
#pragma unroll
            for (int c = 0; c < C; ++c) {
                if (FULL || k0 + c < P) {
                    // FULL: slen = 32 * C / lanes is a multiple of C, so only a
                    // lane's first column can start a stripe.
                    if (!FULL || c == 0) fp &= keep[c];
                    const int32_t sc = wave_score(tab_lo[c], tab_hi[c], sel);
                    int32_t vp = __vimax_s32_relu(diag + sc, e[c]);
                    vp = max(vp, fp);
                    const int32_t h = max(vp, f);
                    e[c] = __vimax_s32_relu(e[c] - e_del, vp - oe_del);
                    fp = __vimax_s32_relu(fp - e_ins, vp - oe_ins);
                    // max(f - e_ins, h - oe_ins) with h = max(vp, f) is
                    // max(f - e_ins, vp - oe_ins), as f - oe_ins <= f - e_ins
                    // (o_ins >= 0, checked at launch); vp - oe_ins is fp's.
                    f = __vimax_s32_relu(f - e_ins, vp - oe_ins);
                    key = max(key, vp * 1024 + kinv[c]);
                    diag = hprev[c];
                    hprev[c] = h;
                    h_last = h;
                }
            }
            diag_left = in_hf & 0xffff;
            send_hf = (int32_t)((uint32_t)h_last | (uint32_t)f << 16), send_fp = fp, send_key = key;
            if (lane == last_lane && !stop) {
                const int32_t mx = key >> 10;
                if (rows) sts_u16(rows_at + 2 * (uint32_t)r, (uint32_t)(mx > 0xffff ? 0xffff : mx));
                if (mx > gmax) {
                    gmax = mx, te = r, qe = 1023 - (key & 1023);
                    if ((size == 1 && gmax + shift >= 255) || gmax >= endsc) {
                        stop = true;
                        processed = (uint32_t)r + 1;
                    }
                }
            }
        }
        if (s % WAVE_STOP_EVERY == WAVE_STOP_EVERY - 1 && __shfl_sync(full, (int)stop, last_lane)) break;
    }
    wave_pass out;
    out.gmax = __shfl_sync(full, gmax, last_lane);
    out.te = __shfl_sync(full, te, last_lane);
    out.qe = __shfl_sync(full, qe, last_lane);
    out.processed = __shfl_sync(full, processed, last_lane);
    return out;
}

template <int C>
__device__ wave_pass wave_run_any(const uint8_t* query, uint32_t qlen, const uint8_t* target, uint32_t tlen,
                                  const int8_t* mat, int32_t m, int32_t o_del, int32_t e_del, int32_t o_ins,
                                  int32_t e_ins, int32_t xtra, int32_t size, bool reverse, uint16_t* rows) {
    const int32_t lanes = size == 1 ? 16 : 8;
    const int32_t P = ((int32_t)qlen + lanes - 1) / lanes * lanes;
    // Warp-uniform: every lane sees the same job.
    if (P == 32 * C)
        return wave_run<C, true>(query, qlen, target, tlen, mat, m, o_del, e_del, o_ins, e_ins, xtra, size,
                                 reverse, rows);
    return wave_run<C, false>(query, qlen, target, tlen, mat, m, o_del, e_del, o_ins, e_ins, xtra, size,
                              reverse, rows);
}

// second_best over record_row's merged row list for one job, by the whole
// warp; `mx(i)` is row i's maximum. The sequential replay is, over rows
// i < processed with s(i) = mx(i) >= minsc:
//   if (have && E + 1 == i) { if (S < s(i)) S = s(i), E = i; }
//   else { if (have) consider(S, E); have = true, S = s(i), E = i; }
// and consider(S, E) at the end: score2/te2 become the S/E of the first
// considered pair with E outside [low, high] and the largest S.
//
// Call row i active when it set (S, E). Then E = i - 1 exactly when row
// i - 1 is active, so active(i) = q(i) && (!active(i - 1) || s(i) > s(i - 1)),
// with q(i) = s(i) >= minsc, and the pairs considered are (s(e), e) for the
// last row e of each run of consecutive active rows. Each lane takes a
// contiguous chunk of rows: the map from the state entering a chunk to the
// state leaving it composes across lanes by a shuffle scan, and each lane
// then finds the run ends in its chunk.
template <typename Row>
__device__ void wave_second_best(nomics_ksw_result* result, int32_t xtra, int32_t profile_max,
                                 uint32_t processed, Row mx) {
    const int lane = threadIdx.x & 31;
    const unsigned full = 0xffffffffu;
    const int32_t minsc = (xtra & KSW_XSUBO) ? (xtra & 0xffff) : 0x10000;
    const int32_t radius = (result->score + profile_max - 1) / profile_max;
    const int32_t low = result->te - radius, high = result->te + radius;
    const uint32_t per = (processed + 31) / 32;
    const uint32_t lo = min((uint32_t)lane * per, processed), hi = min(lo + per, processed);
    // This chunk's map as (f(false), f(true)).
    bool f0 = false, f1 = true;
    int32_t prev = lo > 0 ? mx(lo - 1) : 0;
    for (uint32_t i = lo; i < hi; ++i) {
        const int32_t s = mx(i);
        const bool q = s >= minsc, g1 = q && s > prev;  // g(false) = q, g(true) = g1
        f0 = f0 ? g1 : q;
        f1 = f1 ? g1 : q;
        prev = s;
    }
    // Inclusive scan: lane l gets F_l o ... o F_0, then the state entering
    // its chunk is lane l - 1's map applied to false (no row before row 0).
    for (int d = 1; d < 32; d <<= 1) {
        const bool g0 = __shfl_up_sync(full, (int)f0, d), g1 = __shfl_up_sync(full, (int)f1, d);
        if (lane >= d) {
            const bool n0 = g0 ? f1 : f0, n1 = g1 ? f1 : f0;
            f0 = n0, f1 = n1;
        }
    }
    bool active = __shfl_up_sync(full, (int)f0, 1);
    if (lane == 0) active = false;
    // Run ends in this chunk: active(i) && !active(i + 1).
    prev = lo > 0 ? mx(lo - 1) : 0;
    int32_t s = lo < hi ? mx(lo) : 0;
    uint32_t best_s = 0, best_e = 0xffffffffu;  // best_s = S + 1 of the best eligible end
    for (uint32_t i = lo; i < hi; ++i) {
        const bool now = s >= minsc && (!active || s > prev);
        const int32_t next = i + 1 < processed ? mx(i + 1) : 0;
        const bool after = i + 1 < processed && next >= minsc && (!now || next > s);
        if (now && !after && ((int32_t)i < low || (int32_t)i > high) && (uint32_t)s + 1 > best_s)
            best_s = (uint32_t)s + 1, best_e = i;
        active = now;
        prev = s;
        s = next;
    }
    const uint32_t top = __reduce_max_sync(full, best_s);
    const uint32_t end = __reduce_min_sync(full, best_s == top ? best_e : 0xffffffffu);
    if (top > 0 && (int32_t)top - 1 > result->score2) {
        result->score2 = (int32_t)top - 1;
        result->te2 = (int32_t)end;
    }
}

__device__ void mat_range(const int8_t* mat, int32_t m, int32_t* shift, int32_t* profile_max) {
    int32_t minimum = 127, maximum = 0;
    for (int k = 0; k < m * m; ++k) {
        minimum = minimum < (int32_t)mat[k] ? minimum : (int32_t)mat[k];
        maximum = maximum > (int32_t)mat[k] ? maximum : (int32_t)mat[k];
    }
    *shift = (int32_t)(uint8_t)(0u - (uint8_t)(int8_t)minimum);
    *profile_max = maximum;
}

// One job on one warp. `scratch` is the warp's shared memory, at least
// 3 * round_up(target_len, 16) bytes: the target, then the row maxima.
template <int C>
__device__ void wave_job(const nomics_ksw_job& job, const uint8_t* queries, const uint8_t* targets,
                         const int8_t* matrices, uint8_t* scratch, nomics_ksw_result* out) {
    const int lane = threadIdx.x & 31;
    const uint8_t* query = queries + job.query_offset;
    const int8_t* mat = matrices + job.mat_offset;
    // The target and the row maxima live in shared memory: the inner loop
    // then touches no host-mapped memory.
    uint8_t* target = scratch;
    uint16_t* rows = (uint16_t*)(target + ((job.target_len + 15) & ~15u));
    __syncwarp();  // a previous job on this warp may still read scratch
    for (uint32_t i = lane; i < job.target_len; i += 32) target[i] = targets[job.target_offset + i];
    __syncwarp();

    const wave_pass fwd = wave_run_any<C>(query, job.query_len, target, job.target_len, mat, job.m, job.o_del,
                                          job.e_del, job.o_ins, job.e_ins, job.xtra, job.size, false, rows);
    int32_t shift, profile_max;
    mat_range(mat, job.m, &shift, &profile_max);
    nomics_ksw_result result = {0, -1, -1, -1, -1, -1, -1};
    result.te = fwd.te;
    if (job.size == 1 && fwd.gmax + shift >= 255) {
        result.score = 255;
    } else {
        result.score = fwd.gmax;
        result.qe = fwd.qe;
        __syncwarp();  // rows[] was written by another lane
        if (profile_max > 0)
            wave_second_best(&result, job.xtra, profile_max, fwd.processed,
                             [&](uint32_t i) { return (int32_t)rows[i]; });
    }
    if ((job.xtra & KSW_XSTART) && !((job.xtra & KSW_XSUBO) && result.score < (job.xtra & 0xffff)) &&
        result.qe >= 0) {
        const wave_pass rev =
            wave_run_any<C>(query, (uint32_t)(result.qe + 1), target, (uint32_t)(result.te + 1), mat, job.m,
                            job.o_del, job.e_del, job.o_ins, job.e_ins, KSW_XSTOP | result.score, job.size, true,
                            nullptr);
        if (rev.gmax == result.score) {
            result.tb = result.te - rev.te;
            result.qb = result.qe - rev.qe;
        }
    }
    if (lane == 0) *out = result;
}

template <int C>
__global__ void ksw_wave_kernel(const uint8_t* queries, const uint8_t* targets, const int8_t* matrices,
                                const nomics_ksw_job* jobs, uint32_t n, nomics_ksw_result* out,
                                uint32_t shared_per_warp) {
    const uint32_t id = (blockIdx.x * blockDim.x + threadIdx.x) >> 5;
    if (id >= n) return;
    extern __shared__ uint8_t smem[];
    wave_job<C>(jobs[id], queries, targets, matrices, smem + (threadIdx.x >> 5) * shared_per_warp, out + id);
}

// ---------------------------------------------------------------------------
// Two jobs per warp, one in each 16-bit half of every value (VIADD.16x2 and
// VIMNMX.S16x2 are single instructions on sm_121). The recurrences are
// wave_run's, so each half gives the scalar result. A pair qualifies when
// both jobs are byte-sized (size 1: every score stays below 2^15, see
// pair_ok), have the FULL shape P == 32 * C, and their targets hold only
// bases 0-3; otherwise the warp runs the two jobs one after the other with
// wave_job.
//
// The score of a cell pair is one PRMT: operand a holds job A's scores for
// target bases 0-3 at this column, b job B's, and the selector picks byte tA
// of a and byte 4 + tB of b, each sign-extended to 16 bits.
//
// The row maximum and its column travel as a 16-bit key per half,
// min(vp, 255) * 256 + (255 - k), maximised unsigned: the largest vp, then
// the smallest k, as in wave_run. Saturating vp at 255 changes no result:
// a byte-sized job stops at the first row whose maximum reaches
// thr = min(255 - shift, endsc) <= 255, so rows before it are below 255 and
// the stopping row's maximum is only compared against thr. The last lane
// stores each row's keys; the maxima, the first row reaching each maximum
// and the stopping row are then found by a scan of the stored keys, which
// is the strict `mx > gmax` scan of wave_run.
// ---------------------------------------------------------------------------

__device__ __forceinline__ uint32_t pack16(int32_t lo, int32_t hi) {
    return ((uint32_t)lo & 0xffffu) | (uint32_t)hi << 16;
}

struct pair_side {
    const uint8_t* query;
    const int8_t* mat;
    uint32_t qlen, tlen;  // this pass's query and target lengths
    int32_t m, te0;       // te0: the reverse pass reads target[te0 - r]
    int32_t thr;          // stop at the first row whose maximum reaches thr
};

struct pair_pass {
    int32_t gmax, te, qe;
    uint32_t processed;
};

// Run both sides' passes; keys[r] gets row r's two keys. Returns the rows
// stored (every row of a side that did not stop is stored).
template <int C, bool REV>
__device__ uint32_t pair_run(const pair_side& a, const pair_side& b, const uint16_t* t2, uint32_t gaps_del,
                             uint32_t gaps_ode, uint32_t gaps_ins, uint32_t gaps_oie, uint32_t* keys) {
    // Held in registers through the loop (else rebuilt from threadIdx each step).
    const int lane = (int)wave_opaque(threadIdx.x & 31);
    const uint32_t t2_at = smem_addr(t2), keys_at = smem_addr(keys);
    const unsigned full = 0xffffffffu;
    const int32_t slen = 2 * C;  // P = 32 * C = 16 * slen
    const int32_t k0 = lane * C;
    uint32_t tab_a[C], tab_b[C], keep[C], kinv[C], hprev[C], e[C];
#pragma unroll
    for (int c = 0; c < C; ++c) {
        const int32_t k = k0 + c;
        keep[c] = wave_opaque(k % slen == 0 ? 0u : ~0u);
        kinv[c] = wave_opaque((uint32_t)(255 - k) * 0x10001u);
        tab_a[c] = tab_b[c] = 0;
        if (k < (int32_t)a.qlen) {
            const int32_t q = (int32_t)query_base(a.query, a.qlen, (uint32_t)k, REV);
            for (int t = 0; t < 4 && t < a.m; ++t) tab_a[c] |= (uint32_t)(uint8_t)a.mat[t * a.m + q] << (8 * t);
        }
        if (k < (int32_t)b.qlen) {
            const int32_t q = (int32_t)query_base(b.query, b.qlen, (uint32_t)k, REV);
            for (int t = 0; t < 4 && t < b.m; ++t) tab_b[c] |= (uint32_t)(uint8_t)b.mat[t * b.m + q] << (8 * t);
        }
        hprev[c] = 0;
        e[c] = 0;
    }
    const uint32_t tmax = a.tlen > b.tlen ? a.tlen : b.tlen;
    const uint32_t thr_a = (uint32_t)a.thr << 8, thr_b = (uint32_t)b.thr << 8;
    uint32_t send_h = 0, send_f = 0, send_fp = 0, send_key = 0, diag_left = 0;
    uint32_t kmax = 0;  // last lane: running maximum key of each side
    const int32_t steps = (int32_t)tmax + 31;
    int32_t s = 0;
    for (; s < steps; ++s) {
        uint32_t in_h = __shfl_up_sync(full, send_h, 1);
        uint32_t in_f = __shfl_up_sync(full, send_f, 1);
        uint32_t in_fp = __shfl_up_sync(full, send_fp, 1);
        uint32_t in_key = __shfl_up_sync(full, send_key, 1);
        if (lane == 0) in_h = in_f = in_fp = in_key = 0;
        const int32_t r = s - lane;
        if (r >= 0 && r < (int32_t)tmax) {
            uint32_t t;
            if (REV) {
                // Rows past a side's te0 are past its target: any base will do.
                const uint32_t ta = r <= a.te0 ? lds_u16(t2_at + 2 * (uint32_t)(a.te0 - r)) & 0xffu : 0u;
                const uint32_t tb = r <= b.te0 ? lds_u16(t2_at + 2 * (uint32_t)(b.te0 - r)) >> 8 : 0u;
                t = ta | tb << 8;
            } else {
                t = lds_u16(t2_at + 2 * (uint32_t)r);
            }
            const uint32_t sel = t * 0x11u + 0xc480u;
            uint32_t diag = diag_left, f = in_f, fp = in_fp, key = in_key, h = in_h;
#pragma unroll
            for (int c = 0; c < C; ++c) {
                if (c == 0) fp &= keep[c];  // slen = 2 * C: only column 0 starts a stripe
                uint32_t sc;
                asm("prmt.b32 %0, %1, %2, %3;" : "=r"(sc) : "r"(tab_a[c]), "r"(tab_b[c]), "r"(sel));
                uint32_t vp = __vimax_s16x2_relu(__vadd2(diag, sc), e[c]);
                vp = __vmaxs2(vp, fp);
                h = __vmaxs2(vp, f);
                // f: max(f - e_ins, h - oe_ins) = max(f - e_ins, vp - oe_ins),
                // as in wave_run.
                const uint32_t t_ins = __vsub2(vp, gaps_oie);
                e[c] = __vimax_s16x2_relu(__vsub2(e[c], gaps_del), __vsub2(vp, gaps_ode));
                fp = __vimax_s16x2_relu(__vsub2(fp, gaps_ins), t_ins);
                f = __vimax_s16x2_relu(__vsub2(f, gaps_ins), t_ins);
                key = __vmaxu2(key, __vminu2(vp, 0x00ff00ffu) * 256u + kinv[c]);
                diag = hprev[c];
                hprev[c] = h;
            }
            diag_left = in_h;
            send_h = h, send_f = f, send_fp = fp, send_key = key;
            if (lane == 31) {
                sts_u32(keys_at + 4 * (uint32_t)r, key);
                kmax = __vmaxu2(kmax, key);
            }
        }
        if ((s & 7) == 7) {
            // Lane 31 has stored rows [0, s - 30): a side is done once it
            // reached its threshold or all its rows are stored.
            const int32_t stored = s - 30;
            const bool done_a = (kmax & 0xffffu) >= thr_a || stored >= (int32_t)a.tlen;
            const bool done_b = (kmax >> 16) >= thr_b || stored >= (int32_t)b.tlen;
            if (__shfl_sync(full, (int)(done_a && done_b), 31)) {
                ++s;
                break;
            }
        }
    }
    const int32_t stored = s - 31;
    return stored < 0 ? 0u : (uint32_t)stored < tmax ? (uint32_t)stored : tmax;
}

// One side's result from the stored keys (`half` 0 or 1): the stopping row,
// the maximum over the rows before it with the first row reaching it, and
// that row's column. All lanes take part and get the result.
__device__ pair_pass pair_scan(const uint32_t* keys, uint32_t stored, uint32_t tlen, int32_t thr, int half) {
    const int lane = threadIdx.x & 31;
    const unsigned full = 0xffffffffu;
    const uint32_t rows = stored < tlen ? stored : tlen;
    uint32_t first = rows;
    for (uint32_t i = lane; i < rows; i += 32) {
        const int32_t mx = (int32_t)((keys[i] >> (16 * half + 8)) & 0xffu);
        if (mx >= thr) {
            first = i;
            break;
        }
    }
    pair_pass p;
    first = __reduce_min_sync(full, first);
    p.processed = first < rows ? first + 1 : rows;
    // (mx, first row) as mx * 2^16 + (2^16 - 1 - row): the maximum is the
    // largest mx at its first row. Rows are < 2^16 (targets <= 4000).
    uint32_t best = 0;
    for (uint32_t i = lane; i < p.processed; i += 32) {
        const uint32_t mx = (keys[i] >> (16 * half + 8)) & 0xffu;
        const uint32_t v = mx << 16 | (0xffffu - i);
        best = best > v ? best : v;
    }
    best = __reduce_max_sync(full, best);
    p.gmax = (int32_t)(best >> 16);
    if (p.gmax > 0) {
        p.te = (int32_t)(0xffffu - (best & 0xffffu));
        p.qe = 255 - (int32_t)((keys[p.te] >> (16 * half)) & 0xffu);
    } else {
        p.te = -1;  // wave_run's start: no row beat gmax = 0
        p.qe = 0;
    }
    return p;
}

// Whether job j can take a half of a pair with the FULL shape at C: size 1,
// P == 32 * C, m <= 8, and every value below 2^15 (H <= qlen * max score;
// E, F and Fp are at most H) with gaps small enough not to wrap.
template <int C>
__device__ bool pair_ok(const nomics_ksw_job& j, const int8_t* matrices) {
    if (j.size != 1 || j.m > 8 || j.m < 1) return false;
    if ((j.query_len + 15) / 16 * 16 != 32 * C) return false;
    int32_t shift, profile_max;
    mat_range(matrices + j.mat_offset, j.m, &shift, &profile_max);
    (void)shift;
    const auto small = [](int32_t x) { return x >= 0 && x < 8192; };
    return (int32_t)j.query_len * profile_max < 30000 && small(j.o_del) && small(j.e_del) && small(j.o_ins) &&
           small(j.e_ins) && j.target_len <= 0xfff0u;
}

__device__ int32_t pair_thr(const int8_t* mat, int32_t m, int32_t xtra) {
    int32_t shift, profile_max;
    mat_range(mat, m, &shift, &profile_max);
    const int32_t endsc = (xtra & KSW_XSTOP) ? (xtra & 0xffff) : 0x10000;
    int32_t thr = 255 - shift < endsc ? 255 - shift : endsc;
    // wave_run stops only on a row beating gmax >= 0, so at mx >= 1.
    return thr > 1 ? thr : 1;
}

template <int C>
__global__ void ksw_wave_pair_kernel(const uint8_t* queries, const uint8_t* targets, const int8_t* matrices,
                                     const nomics_ksw_job* jobs, uint32_t n, nomics_ksw_result* out,
                                     uint32_t shared_per_warp) {
    const uint32_t pair = (blockIdx.x * blockDim.x + threadIdx.x) >> 5;
    const int lane = threadIdx.x & 31;
    const unsigned full = 0xffffffffu;
    const uint32_t ia = 2 * pair, ib = ia + 1;
    if (ia >= n) return;
    extern __shared__ uint8_t smem[];
    uint8_t* scratch = smem + (threadIdx.x >> 5) * shared_per_warp;
    if (ib >= n) {
        wave_job<C>(jobs[ia], queries, targets, matrices, scratch, out + ia);
        return;
    }
    const nomics_ksw_job ja = jobs[ia], jb = jobs[ib];
    bool ok = pair_ok<C>(ja, matrices) && pair_ok<C>(jb, matrices);
    const uint32_t tmax = ja.target_len > jb.target_len ? ja.target_len : jb.target_len;
    const uint32_t span = (tmax + 15) & ~15u;
    uint16_t* t2 = (uint16_t*)scratch;               // 2 * span bytes: both targets interleaved
    uint32_t* keys = (uint32_t*)(scratch + 2 * span);  // 4 * span bytes: both rows' keys
    if (ok) {
        bool bad = false;
        for (uint32_t i = lane; i < span; i += 32) {
            const uint32_t x = i < ja.target_len ? targets[ja.target_offset + i] : 0u;
            const uint32_t y = i < jb.target_len ? targets[jb.target_offset + i] : 0u;
            bad |= (x | y) >= 4u;
            t2[i] = (uint16_t)(x | y << 8);
        }
        ok = !__any_sync(full, bad);
    }
    if (!ok) {
        wave_job<C>(ja, queries, targets, matrices, scratch, out + ia);
        wave_job<C>(jb, queries, targets, matrices, scratch, out + ib);
        return;
    }
    __syncwarp();
    const nomics_ksw_job* job[2] = {&ja, &jb};
    pair_side side[2];
    for (int h = 0; h < 2; ++h) {
        const nomics_ksw_job& j = *job[h];
        side[h].query = queries + j.query_offset;
        side[h].mat = matrices + j.mat_offset;
        side[h].qlen = j.query_len;
        side[h].tlen = j.target_len;
        side[h].m = j.m;
        side[h].te0 = 0;
        side[h].thr = pair_thr(side[h].mat, j.m, j.xtra);
    }
    const uint32_t del = pack16(ja.e_del, jb.e_del), ode = pack16(ja.o_del + ja.e_del, jb.o_del + jb.e_del);
    const uint32_t ins = pack16(ja.e_ins, jb.e_ins), oie = pack16(ja.o_ins + ja.e_ins, jb.o_ins + jb.e_ins);
    const uint32_t stored = pair_run<C, false>(side[0], side[1], t2, del, ode, ins, oie, keys);
    __syncwarp();  // keys[] was written by lane 31
    nomics_ksw_result result[2];
    bool rev[2];
    for (int h = 0; h < 2; ++h) {
        const nomics_ksw_job& j = *job[h];
        const pair_pass fwd = pair_scan(keys, stored, j.target_len, side[h].thr, h);
        int32_t shift, profile_max;
        mat_range(side[h].mat, j.m, &shift, &profile_max);
        nomics_ksw_result& res = result[h];
        res = {0, -1, -1, -1, -1, -1, -1};
        res.te = fwd.te;
        if (fwd.gmax + shift >= 255) {
            res.score = 255;
        } else {
            res.score = fwd.gmax;
            res.qe = fwd.qe;
            if (profile_max > 0)
                wave_second_best(&res, j.xtra, profile_max, fwd.processed,
                                 [&](uint32_t i) { return (int32_t)((keys[i] >> (16 * h + 8)) & 0xffu); });
        }
        rev[h] = (j.xtra & KSW_XSTART) && !((j.xtra & KSW_XSUBO) && res.score < (j.xtra & 0xffff)) && res.qe >= 0;
    }
    // Reverse passes: paired when both have the FULL shape, else one by one.
    const bool both = rev[0] && rev[1] && (uint32_t)(result[0].qe + 16) / 16 * 16 == 32 * C &&
                      (uint32_t)(result[1].qe + 16) / 16 * 16 == 32 * C;
    __syncwarp();  // the scans above read keys[]
    if (both) {
        for (int h = 0; h < 2; ++h) {
            side[h].qlen = (uint32_t)(result[h].qe + 1);
            side[h].tlen = (uint32_t)(result[h].te + 1);
            side[h].te0 = result[h].te;
            side[h].thr = pair_thr(side[h].mat, side[h].m, KSW_XSTOP | result[h].score);
        }
        const uint32_t stored_rev = pair_run<C, true>(side[0], side[1], t2, del, ode, ins, oie, keys);
        __syncwarp();
        for (int h = 0; h < 2; ++h) {
            const pair_pass p = pair_scan(keys, stored_rev, side[h].tlen, side[h].thr, h);
            if (p.gmax == result[h].score) {
                result[h].tb = result[h].te - p.te;
                result[h].qb = result[h].qe - p.qe;
            }
        }
    } else {
        // The scalar pass reads a byte target: unpack this side's into keys[].
        uint8_t* target = (uint8_t*)keys;
        for (int h = 0; h < 2; ++h) {
            if (!rev[h]) continue;
            const nomics_ksw_job& j = *job[h];
            __syncwarp();
            for (uint32_t i = lane; i < j.target_len; i += 32) target[i] = (uint8_t)(t2[i] >> (8 * h));
            __syncwarp();
            const wave_pass p = wave_run_any<C>(side[h].query, (uint32_t)(result[h].qe + 1), target,
                                                (uint32_t)(result[h].te + 1), side[h].mat, j.m, j.o_del, j.e_del,
                                                j.o_ins, j.e_ins, KSW_XSTOP | result[h].score, j.size, true, nullptr);
            if (p.gmax == result[h].score) {
                result[h].tb = result[h].te - p.te;
                result[h].qb = result[h].qe - p.qe;
            }
        }
    }
    if (lane == 0) out[ia] = result[0];
    if (lane == 1) out[ib] = result[1];
}

template <int C>
static void wave_launch(bool pair, uint32_t n, uint32_t longest, cudaStream_t s, const uint8_t* queries,
                        const uint8_t* targets, const int8_t* matrices, const nomics_ksw_job* jobs,
                        nomics_ksw_result* out) {
    const uint32_t span = (longest + 15) & ~15u;
    if constexpr (C <= 8) {  // the pair keys hold columns < 256
        if (pair) {
            const int threads = 64;  // two warps, two alignments each
            const uint32_t per_warp = 6 * span;  // both targets (2 bytes a row) and keys (4)
            const unsigned blocks = (unsigned)(((uint64_t)(n + 1) / 2 * 32 + threads - 1) / threads);
            ksw_wave_pair_kernel<C><<<blocks, threads, 2 * per_warp, s>>>(queries, targets, matrices, jobs, n,
                                                                        out, per_warp);
            return;
        }
    }
    {
        const int threads = 128;  // four warps, one alignment each
        const uint32_t per_warp = 3 * span;  // target, then 16-bit row maxima
        const unsigned blocks = (unsigned)(((uint64_t)n * 32 + threads - 1) / threads);
        ksw_wave_kernel<C><<<blocks, threads, 4 * per_warp, s>>>(queries, targets, matrices, jobs, n, out,
                                                               per_warp);
    }
}

extern "C" int nomics_ksw_align2_wave_launch(nomics_stream stream, const uint8_t* queries,
                                              const uint8_t* targets, const int8_t* matrices,
                                              const nomics_ksw_job* jobs, uint32_t n, uint16_t* row_maxima,
                                              nomics_ksw_result* out, nomics_ticket* ticket) {
    (void)row_maxima;
    if (ticket == nullptr) return cudaErrorInvalidValue;
    *ticket = nullptr;
    if (n == 0) return cudaSuccess;
    uint32_t widest = 0;  // padded query columns
    uint32_t longest = 0;  // target bases
    for (uint32_t i = 0; i < n; ++i) {
        if (jobs[i].m > 8) return cudaErrorInvalidValue;  // score table holds 8 bases
        // The F update relies on o_ins >= 0 (wave_run); BWA never passes negative gaps.
        if (jobs[i].o_del < 0 || jobs[i].e_del < 0 || jobs[i].o_ins < 0 || jobs[i].e_ins < 0) return cudaErrorInvalidValue;
        // Scoring the kernels do not reproduce goes to the CPU (the caller's
        // exact fallback): BWA's byte path casts e and o + e to uint8_t, and
        // its 16-bit path saturates at 32767 (_mm_adds_epi16) where the
        // kernels add in 32 bits. H <= qlen * max score bounds every cell.
        if (jobs[i].size == 1 && (jobs[i].e_del > 255 || jobs[i].e_ins > 255 ||
                                  jobs[i].o_del + jobs[i].e_del > 255 || jobs[i].o_ins + jobs[i].e_ins > 255))
            return cudaErrorInvalidValue;
        if (jobs[i].size != 1) {
            int32_t top = 0;
            for (int32_t k = 0; k < jobs[i].m * jobs[i].m; ++k)
                top = top > matrices[jobs[i].mat_offset + k] ? top : matrices[jobs[i].mat_offset + k];
            if ((int64_t)jobs[i].query_len * top > 32767) return cudaErrorInvalidValue;
        }
        longest = longest > jobs[i].target_len ? longest : jobs[i].target_len;
        const uint32_t lanes = jobs[i].size == 1 ? 16 : 8;
        const uint32_t padded = (jobs[i].query_len + lanes - 1) / lanes * lanes;
        widest = widest > padded ? widest : padded;
    }
    // Both kernels stage 12 * round_up(longest, 16) bytes per block.
    if (12 * ((longest + 15) & ~15u) > 48 * 1024 || widest > 32 * 16) return cudaErrorInvalidValue;
    // NOMICS_KSW_PAIR=0: one job per warp everywhere.
    static const bool pair = [] {
        const char* v = getenv("NOMICS_KSW_PAIR");
        return !(v && v[0] == '0' && v[1] == 0);
    }();
    nomics_ticket_s* t = nullptr;
    cudaError_t error = ticket_begin((cudaStream_t)stream, &t);
    if (error != cudaSuccess) return error;
    const cudaStream_t s = t->stream;
    if (widest <= 32 * 4)
        wave_launch<4>(pair, n, longest, s, queries, targets, matrices, jobs, out);
    else if (widest <= 32 * 5)
        wave_launch<5>(pair, n, longest, s, queries, targets, matrices, jobs, out);
    else if (widest <= 32 * 6)
        wave_launch<6>(pair, n, longest, s, queries, targets, matrices, jobs, out);
    else if (widest <= 32 * 8)
        wave_launch<8>(pair, n, longest, s, queries, targets, matrices, jobs, out);
    else
        wave_launch<16>(false, n, longest, s, queries, targets, matrices, jobs, out);
    error = cudaGetLastError();
    return ticket_finish(t, error, ticket);
}
