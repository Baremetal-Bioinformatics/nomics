// BWA 0.7.19 ksw_extend2, one target row at a time. Short queries use
// independent 8/16-lane groups; longer queries use a whole warp. Three launch
// filters keep the original job indices, without sorting or scratch buffers.
//
// Each lane permanently owns C consecutive query columns (C is a template
// constant). H/E and the signed-byte query profile stay in registers. At a
// row boundary h[c] is BWA's eh[j].h, NOT H_j: the shifted H_{j-1}. Only
// [beg,end) and the end sentinel are overwritten; stale cells outside that
// interval must survive band shrinkage and subsequent expansion.
//
// M = diag ? diag + score : 0
// F_j = max(F_{j-1} - e_ins, max(M_{j-1} - oe_ins, 0)), F_beg = 0
//     = exclusive prefix max of max(M_k - oe_ins, 0) + (k+1)*e_ins,
//       seeded with beg*e_ins, minus j*e_ins.
// One group scan supplies F at each lane's first column. A single shuffle
// supplies the shifted H at lane boundaries; there is no shared DP storage
// or row barrier. All collectives use the group's mask, including after a
// neighboring group has exited on z-drop or row_max == 0.
//
// The group body, extend_group, lives in extend_group.cuh so that fused.cu
// runs the same arithmetic.
#include "nomics_gpu.h"
#include "ticket.cuh"

#include "extend_group.cuh"

#include <cuda_runtime.h>

static_assert(sizeof(nomics_extend_job) == 72, "C ABI layout changed");
static_assert(sizeof(nomics_extend_result) == 24, "C ABI layout changed");

// One job per group: unpack it, run extend_group, and store from lane 0.
template <int W, int C>
__device__ __forceinline__ void extend_job(
    const uint8_t* queries, const uint8_t* targets, const int8_t* matrices,
    const nomics_extend_job& job, nomics_extend_result* result) {
    const nomics_extend_result r = extend_group<W, C>(
        queries + job.query_offset, targets + job.target_offset, (int)job.query_len,
        (int)job.target_len, job.m, matrices + job.mat_offset, job.o_del, job.e_del,
        job.o_ins, job.e_ins, job.w, job.end_bonus, job.zdrop, job.h0);
    if ((threadIdx.x & (W - 1)) == 0) *result = r;
}

// Each pass reads just query_len for rejected jobs. No compaction allocation,
// host sorting, or changed result order. C is uniform across a physical warp,
// so neighboring groups cannot serialize different template switch arms.
template <int W>
__global__ void extend_kernel(const uint8_t* queries, const uint8_t* targets, const int8_t* matrices,
                              const nomics_extend_job* jobs, uint32_t n, nomics_extend_result* out) {
    const uint64_t id = ((uint64_t)blockIdx.x * blockDim.x + threadIdx.x) / W;
    const unsigned qlen = id < n ? jobs[id].query_len : 0;
    const unsigned lower = W == 8 ? 0 : (W == 16 ? 32 : 64);
    const unsigned upper = W == 8 ? 32 : (W == 16 ? 64 : 256);
    const bool eligible = qlen > lower && qlen <= upper;
    const unsigned group_q = eligible ? qlen : 0;
    const unsigned warp_q = __reduce_max_sync(0xffffffffu, group_q);
    if (!eligible) return;
    const nomics_extend_job job = jobs[id];
    const int columns = (warp_q + W - 1) / W;
#define EXT_CASE(C) case C: extend_job<W, C>(queries, targets, matrices, job, out + id); break
    switch (columns) {
        EXT_CASE(1);
        EXT_CASE(2);
        EXT_CASE(3);
        EXT_CASE(4);
        // Only full-warp jobs need more than four columns per lane.
        case 5: if constexpr (W == 32) extend_job<W, 5>(queries, targets, matrices, job, out + id); break;
        case 6: if constexpr (W == 32) extend_job<W, 6>(queries, targets, matrices, job, out + id); break;
        case 7: if constexpr (W == 32) extend_job<W, 7>(queries, targets, matrices, job, out + id); break;
        case 8: if constexpr (W == 32) extend_job<W, 8>(queries, targets, matrices, job, out + id); break;
    }
#undef EXT_CASE
}

extern "C" int nomics_ksw_extend2_launch(nomics_stream stream, const uint8_t* queries, const uint8_t* targets,
                                         const int8_t* matrices, const nomics_extend_job* jobs, uint32_t n,
                                         nomics_extend_result* out, nomics_ticket* ticket) {
    if (ticket == nullptr) return cudaErrorInvalidValue;
    *ticket = nullptr;
    if (n == 0) return cudaSuccess;
    bool have_short = false, have_medium = false, have_long = false;
    for (uint32_t i = 0; i < n; ++i) {
        if (jobs[i].query_len == 0 || jobs[i].query_len > 32 * EXT_MAX_C || jobs[i].m * jobs[i].m > 64)
            return cudaErrorInvalidValue;
        have_short |= jobs[i].query_len <= 32;
        have_medium |= jobs[i].query_len > 32 && jobs[i].query_len <= 64;
        have_long |= jobs[i].query_len > 64;
    }
    const int threads = 128;
    // Each caller owns its buffers and stream (nonblocking: no serialization
    // with other launches or the legacy default stream).
    nomics_ticket_s* t = nullptr;
    cudaError_t error = ticket_begin((cudaStream_t)stream, &t);
    if (error != cudaSuccess) return error;
    const cudaStream_t s = t->stream;
    // The ticket times all launch filters, including their skipped jobs.
    if (error == cudaSuccess && have_short) {
        const unsigned blocks = (unsigned)(((uint64_t)n * 8 + threads - 1) / threads);
        extend_kernel<8><<<blocks, threads, 0, s>>>(queries, targets, matrices, jobs, n, out);
        error = cudaGetLastError();
    }
    if (error == cudaSuccess && have_medium) {
        const unsigned blocks = (unsigned)(((uint64_t)n * 16 + threads - 1) / threads);
        extend_kernel<16><<<blocks, threads, 0, s>>>(queries, targets, matrices, jobs, n, out);
        error = cudaGetLastError();
    }
    if (error == cudaSuccess && have_long) {
        const unsigned blocks = (unsigned)(((uint64_t)n * 32 + threads - 1) / threads);
        extend_kernel<32><<<blocks, threads, 0, s>>>(queries, targets, matrices, jobs, n, out);
        error = cudaGetLastError();
    }
    // On error this drains the earlier launches before the Rust owner may
    // reuse or free the host buffers.
    return ticket_finish(t, error, ticket);
}
