// A seed's whole BWA extension in one GPU job: chain2aln's left band-retry
// loop, then its right one (crates/nomics-align/src/extend.rs), up to four
// ksw_extend2 calls. Design: docs/strategy/gpu-extend-fused-design.md.
//
//   left (lq_len > 0):  prev = -1; for i in 0..max_band_try:
//                         w = p.w << i; x = ksw_extend2(qL, tL, w, h0, end_bonus_left)
//                         if x.score == prev || x.max_off < (w>>1)+(w>>2): break
//                         prev = x.score
//                       sc0 = x.score            (no left: sc0 = h0)
//   right (rq_len > 0): the same loop with prev = sc0, h0 = sc0, end_bonus_right.
//
// Every call writes its result; n_left / n_right count them and right_h0 = sc0
// (also when the right is skipped). Unused slots are zeroed so dumps of the
// result buffer are deterministic.
//
// Kernel v3 (v2 was one launch per group width with C = the warp's widest
// need; measurements in the design doc):
//  - Length bins. The launcher sorts job indices into bins by
//    qmax = max(lq_len, rq_len); each bin is one launch of
//    fused_bin_kernel<W, C> with a fixed group width W and C columns per lane,
//    W * C >= the bin's upper bound, so each instantiation carries only its own
//    C's registers. Within a bin, jobs are ordered by descending target rows (a
//    stable counting sort), so groups sharing a warp run similar row counts. A
//    group reads its job index from `order` and writes out[index]: results stay
//    in input order.
//  - The call body is extend_group_v3 (extend_v3.cuh): extend_group's
//    arithmetic with one profile word per column, branch-free columns, F
//    carried without multiplications.
//  - Scoring is copied into scalars and the job is read at each call site, so
//    nothing lives on the stack and the job's 56 bytes are not held in
//    registers through the row loop.
//
// Exactness. Each call is extend_group_v3<W, C>, which computes what
// extend_group (the body behind nomics_ksw_extend2) computes (argued in
// extend_v3.cuh), with the same argument values a nomics_extend_job for that
// call would carry: pointers to the call's first base, ints from the job, the
// unclipped band w (the body clips its own copy, and the break test uses the
// unclipped w as BWA's aw[] does). The body is exact for any W and any C with
// 0 < qlen <= W*C: lanes and columns past qlen stay inert. Its six results are
// group-uniform, so every lane evaluates the break test identically and the
// group never splits. sc0 >= h0 > 0, since score starts at h0 and only grows:
// the right call keeps BWA's h0 > 0 precondition. Collectives use the group
// mask, so a group that finishes early (or takes a retry its neighbours
// don't), or a partial last warp of a bin, never stalls or corrupts its
// neighbours. Scoring outside extend_group_v3's preconditions (m <= 5, matrix
// entries in [-32, 31], gap penalties in range) is refused with
// cudaErrorInvalidValue; the caller then extends on the CPU.
//
// One call site. The four possible calls go through a single extend_group_v3
// instance per <W, C>, driven by a side/band loop.
//
// Parameters are a by-value kernel argument, not __constant__: callers may
// run concurrently on private streams (the nomics_ksw_extend2 contract), and a
// __constant__ symbol is one per module. The 64-byte matrix is copied once per
// block into shared memory, so the profile build's per-lane (divergent)
// matrix reads hit shared memory, not the parameter bank.
//
// Tuning knobs (development only; defaults are the measured best on the GB10):
//   NOMICS_FUSED_BINS="WxC:hi,..."  bins in ascending hi, the last hi = 256,
//                                   hi <= W*C, W in {8,16,32},
//                                   C in 1..8, W*C <= 256
//   NOMICS_FUSED_SORT=0             bins only, input order within a bin
//   NOMICS_FUSED_THREADS=n          block size (multiple of 32, <= 256; default 32)
//   NOMICS_FUSED_PRIORITY=0         default-priority stream (the default is the
//                                   highest priority, so pending extension blocks
//                                   go ahead of pending seeding blocks)
// A malformed value makes every call fail with cudaErrorInvalidValue.
#include "nomics_gpu.h"
#include "ticket.cuh"

#include "extend_group.cuh"
#include "extend_v3.cuh"

#include <cuda_runtime.h>

#include <stdlib.h>
#include <string.h>

#include <vector>

static_assert(sizeof(nomics_fused_params) == 104, "C ABI layout changed");
static_assert(sizeof(nomics_fused_job) == 56, "C ABI layout changed");
static_assert(sizeof(nomics_fused_result) == 112, "C ABI layout changed");

#define FUSED_MAX_BAND_TRY 2  // slots per side in nomics_fused_result

// nomics_fused_params minus the matrix, as registers.
struct FusedScoring {
    int m, o_del, e_del, o_ins, e_ins, w, max_band_try, end_bonus_left, end_bonus_right, zdrop;
};

template <int W, int C>
__device__ __forceinline__ void fused_group(
    const uint8_t* queries, const uint8_t* targets, const int8_t* mat,
    const FusedScoring p, const nomics_fused_job& job, nomics_fused_result* out) {
    const int lane = threadIdx.x & (W - 1);
    const int lq = (int)job.lq_len, rq = (int)job.rq_len;
    // side 0 = left, 1 = right, 2 = done. All of these are group-uniform.
    int side = lq > 0 ? 0 : (rq > 0 ? 1 : 2);
    int band = 0, prev = -1, sc0 = job.h0, n_left = 0, n_right = 0;
    if (side == 1) prev = sc0;
    while (side < 2) {
        const bool right = side == 1;
        const int w = p.w << band;
        const nomics_extend_result x = extend_group_v3<W, C>(
            queries + (right ? job.rq_off : job.lq_off),
            targets + (right ? job.rt_off : job.lt_off),
            right ? rq : lq, right ? (int)job.rt_len : (int)job.lt_len,
            p.m, mat, p.o_del, p.e_del, p.o_ins, p.e_ins, w,
            right ? p.end_bonus_right : p.end_bonus_left, p.zdrop,
            right ? sc0 : job.h0);
        nomics_extend_result* const slots = right ? out->right : out->left;
        if (lane == 0) slots[band] = x;
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
        nomics_extend_result none;
        none.score = none.qle = none.tle = none.gtle = none.gscore = none.max_off = 0;
        for (int k = n_left; k < FUSED_MAX_BAND_TRY; ++k) out->left[k] = none;
        for (int k = n_right; k < FUSED_MAX_BAND_TRY; ++k) out->right[k] = none;
        out->n_left = n_left;
        out->n_right = n_right;
        out->right_h0 = sc0;
        out->pad = 0;
    }
}

// One bin: group g runs job order[g]. Every job of the bin has
// 0 < max(lq_len, rq_len) <= W * C (checked on the host).
// Register cap: up to five columns per lane the uncapped allocation is
// 75-91 registers; capping at 80 trades a few rematerializations for
// occupancy (measured faster on the GB10, no spills). Wider columns keep
// what the compiler picks.
template <int W, int C>
__global__ void __maxnreg__(C <= 5 ? 80 : 255) fused_bin_kernel(const uint8_t* queries, const uint8_t* targets,
                                 const __grid_constant__ nomics_fused_params p,
                                 const nomics_fused_job* jobs,
                                 const uint32_t* order, uint32_t count, nomics_fused_result* out) {
    // Before any return: every thread of the block reaches the barrier.
    __shared__ int8_t mat[64];
    for (int k = threadIdx.x; k < 64; k += blockDim.x) mat[k] = p.mat[k];
    __syncthreads();
    const uint64_t g = ((uint64_t)blockIdx.x * blockDim.x + threadIdx.x) / W;
    if (g >= count) return;
    const uint32_t id = order[g];
    // A reference, not a copy: the call site reloads the side it needs
    // instead of holding all 56 bytes in registers through the row loop.
    const nomics_fused_job& job = jobs[id];
    const FusedScoring s = {p.m, p.o_del, p.e_del, p.o_ins, p.e_ins, p.w, p.max_band_try,
                            p.end_bonus_left, p.end_bonus_right, p.zdrop};
    fused_group<W, C>(queries, targets, mat, s, job, out + id);
}

namespace {

struct Bin {
    int w, c, hi;
};

bool supported(int w, int c) {
    const bool cw = w == 8 || w == 16 || w == 32;
    const bool cc = c >= 1 && c <= 8;
    return cw && cc && w * c <= 32 * EXT_MAX_C;
}

// Measured on the GB10 (docs/strategy/gpu-extend-fused-design.md, kernel v3).
const Bin kDefaultBins[] = {
    {8, 1, 8},   {8, 2, 16},  {8, 3, 24},  {8, 4, 32},  {16, 3, 48},  {16, 4, 64},
    {16, 5, 80}, {32, 3, 96},  {32, 4, 128}, {32, 5, 160}, {32, 6, 192}, {32, 7, 224}, {32, 8, 256},
};

struct Config {
    bool ok = true;
    std::vector<Bin> bins;
    bool sort = true;
    int threads = 32;
};

Config read_config() {
    Config cfg;
    cfg.bins.assign(kDefaultBins, kDefaultBins + sizeof(kDefaultBins) / sizeof(kDefaultBins[0]));
    if (const char* s = getenv("NOMICS_FUSED_BINS")) {
        cfg.bins.clear();
        const char* at = s;
        int lo = 0;
        while (*at) {
            char* end = nullptr;
            const long w = strtol(at, &end, 10);
            if (*end != 'x') { cfg.ok = false; break; }
            const long c = strtol(end + 1, &end, 10);
            if (*end != ':') { cfg.ok = false; break; }
            const long hi = strtol(end + 1, &end, 10);
            // Range-check the longs before narrowing them.
            if (w < 1 || w > 32 || c < 1 || c > EXT_MAX_C || hi <= lo || hi > w * c ||
                !supported((int)w, (int)c)) {
                cfg.ok = false;
                break;
            }
            cfg.bins.push_back({(int)w, (int)c, (int)hi});
            lo = (int)hi;
            if (*end == ',') ++end;
            else if (*end) { cfg.ok = false; break; }
            at = end;
        }
        if (cfg.bins.empty() || lo != 32 * EXT_MAX_C) cfg.ok = false;
    }
    if (const char* s = getenv("NOMICS_FUSED_SORT")) cfg.sort = strcmp(s, "0") != 0;
    if (const char* s = getenv("NOMICS_FUSED_THREADS")) {
        cfg.threads = atoi(s);
        if (cfg.threads < 32 || cfg.threads > 256 || cfg.threads % 32) cfg.ok = false;
    }
    return cfg;
}

const Config& config() {
    static const Config cfg = read_config();
    return cfg;
}

template <int W, int C>
cudaError_t launch(int threads, cudaStream_t stream, const uint8_t* queries, const uint8_t* targets,
                   const nomics_fused_params& p, const nomics_fused_job* jobs, const uint32_t* order,
                   uint32_t count, nomics_fused_result* out) {
    const unsigned blocks = (unsigned)(((uint64_t)count * W + threads - 1) / threads);
    fused_bin_kernel<W, C><<<blocks, threads, 0, stream>>>(queries, targets, p, jobs, order, count, out);
    return cudaGetLastError();
}

#define FUSED_COLUMNS(X) X(1) X(2) X(3) X(4) X(5) X(6) X(7) X(8)

template <int W>
cudaError_t launch_w(int c, int threads, cudaStream_t stream, const uint8_t* queries,
                     const uint8_t* targets, const nomics_fused_params& p,
                     const nomics_fused_job* jobs, const uint32_t* order, uint32_t count,
                     nomics_fused_result* out) {
    switch (c) {
#define FUSED_CASE(C)                                                                         \
    case C:                                                                                   \
        if constexpr (W * C <= 32 * EXT_MAX_C)                                                \
            return launch<W, C>(threads, stream, queries, targets, p, jobs, order, count, out); \
        break;
        FUSED_COLUMNS(FUSED_CASE)
#undef FUSED_CASE
    }
    return cudaErrorInvalidValue;
}

cudaError_t launch_bin(const Bin& b, int threads, cudaStream_t stream, const uint8_t* queries,
                       const uint8_t* targets, const nomics_fused_params& p,
                       const nomics_fused_job* jobs, const uint32_t* order, uint32_t count,
                       nomics_fused_result* out) {
    switch (b.w) {
        case 8: return launch_w<8>(b.c, threads, stream, queries, targets, p, jobs, order, count, out);
        case 16: return launch_w<16>(b.c, threads, stream, queries, targets, p, jobs, order, count, out);
        case 32: return launch_w<32>(b.c, threads, stream, queries, targets, p, jobs, order, count, out);
    }
    return cudaErrorInvalidValue;
}

constexpr int kRowBuckets = 64;  // 16 target rows each; the last is open-ended

}  // namespace

extern "C" int nomics_fused_extend_launch(nomics_stream stream, const uint8_t* queries, const uint8_t* targets,
                                          const nomics_fused_params* params, const nomics_fused_job* jobs,
                                          uint32_t n, uint32_t* order, nomics_fused_result* out,
                                          nomics_ticket* ticket) {
    if (ticket == nullptr) return cudaErrorInvalidValue;
    *ticket = nullptr;
    if (n == 0) return cudaSuccess;
    const Config& cfg = config();
    if (!cfg.ok || params == nullptr || order == nullptr || params->max_band_try < 1 ||
        params->max_band_try > FUSED_MAX_BAND_TRY)
        return cudaErrorInvalidValue;
    const nomics_fused_params p = *params;
    // extend_group_v3's preconditions (BWA mem's scoring meets them; callers
    // extend on the CPU when a launch is refused).
    if (p.m < 1 || p.m > 5 || p.o_del < 0 || p.o_ins < 0 || p.o_del > 1 << 16 || p.o_ins > 1 << 16 ||
        p.e_del < 1 || p.e_ins < 1 || p.e_del > 1 << 16 || p.e_ins > 1 << 16 ||
        p.end_bonus_left < 0 || p.end_bonus_left > 1 << 16 || p.end_bonus_right < 0 ||
        p.end_bonus_right > 1 << 16 || p.w < 1 || p.w > 1 << 20)
        return cudaErrorInvalidValue;
    for (int k = 0; k < p.m * p.m; ++k)
        if (p.mat[k] < -32 || p.mat[k] > 31) return cudaErrorInvalidValue;
    const int nbins = (int)cfg.bins.size();
    // Bin of each qmax (1..256), then a counting sort by (bin, row bucket).
    int bin_of[32 * EXT_MAX_C + 1];
    for (int q = 1, b = 0; q <= 32 * EXT_MAX_C; ++q) {
        while (q > cfg.bins[b].hi) ++b;
        bin_of[q] = b;
    }
    std::vector<uint32_t> start((size_t)nbins * kRowBuckets + 1, 0);
    auto key = [&](const nomics_fused_job& j) {
        const uint32_t lq = j.lq_len, rq = j.rq_len;
        const uint32_t q = lq > rq ? lq : rq;
        int bucket = 0;
        if (cfg.sort) {
            const uint64_t rows = (lq ? (uint64_t)j.lt_len : 0) + (rq ? (uint64_t)j.rt_len : 0);
            bucket = kRowBuckets - 1 - (int)(rows / 16 < kRowBuckets - 1 ? rows / 16 : kRowBuckets - 1);
        }
        return (size_t)bin_of[q] * kRowBuckets + bucket;
    };
    for (uint32_t i = 0; i < n; ++i) {
        const uint32_t lq = jobs[i].lq_len, rq = jobs[i].rq_len;
        // Target lengths are cast to int in the kernel; reject what would wrap.
        // h0 <= 2^24 keeps every score, F term (+ (j+1)*e_ins) and band bound far
        // from int overflow: scores stay below h0 + 2 * 256 * 31.
        if ((lq == 0 && rq == 0) || lq > 32 * EXT_MAX_C || rq > 32 * EXT_MAX_C ||
            jobs[i].lt_len > (uint32_t)INT_MAX || jobs[i].rt_len > (uint32_t)INT_MAX ||
            jobs[i].h0 < 1 || jobs[i].h0 > 1 << 24)
            return cudaErrorInvalidValue;
        ++start[key(jobs[i]) + 1];
    }
    for (size_t k = 1; k < start.size(); ++k) start[k] += start[k - 1];
    {
        std::vector<uint32_t> at(start.begin(), start.end() - 1);
        for (uint32_t i = 0; i < n; ++i) order[at[key(jobs[i])]++] = i;
    }

    // The stream is the caller's (one per engine, created once): nonblocking,
    // so independent tiles and seeding/rescue do not serialize, and at the
    // greatest priority by default: a lane waits on its own tile's launch,
    // and at default priority that launch queued behind the pending blocks of
    // a ~1.8 s seeding kernel (p90 2.2 s). Scheduling only; results are
    // unaffected.
    nomics_ticket_s* t = nullptr;
    cudaError_t error = ticket_begin((cudaStream_t)stream, &t);
    if (error != cudaSuccess) return error;
    for (int b = 0; b < nbins && error == cudaSuccess; ++b) {
        const uint32_t first = start[(size_t)b * kRowBuckets];
        const uint32_t count = start[(size_t)(b + 1) * kRowBuckets] - first;
        if (count == 0) continue;
        error = launch_bin(cfg.bins[b], cfg.threads, t->stream, queries, targets, p, jobs,
                           order + first, count, out);
    }
    // Even if a later launch failed, this finishes earlier launches before
    // the Rust owner may reuse or free the host buffers.
    return ticket_finish(t, error, ticket);
}
