// Experimental GATK LoglessPairHMM kernel. One warp owns a (read, haplotype)
// pair and walks its matrix by anti-diagonals. Read row i belongs permanently
// to lane (i - 1) % 32, with up to eight rows per lane.
//
// Every floating-point recurrence operation uses an explicit IEEE round-to-
// nearest intrinsic. nvcc also compiles the whole shim with --fmad=false. The
// final read row is encountered in ascending haplotype-column order, so its
// owner accumulates the result with the CPU's exact left-to-right order.
#include "nomics_gpu.h"

#include <cuda_runtime.h>

#define WARP 32
#define ROWS_PER_LANE ((NOMICS_PAIRHMM_MAX_READ_LEN + WARP - 1) / WARP)

static_assert(ROWS_PER_LANE == 5, "update PairHMM row storage");
static_assert(sizeof(nomics_pairhmm_pair) == 40, "C ABI layout changed");

__device__ __forceinline__ double add_rn(double a, double b) { return __dadd_rn(a, b); }
__device__ __forceinline__ double mul_rn(double a, double b) { return __dmul_rn(a, b); }

__global__ void pairhmm_kernel(const uint8_t* read_bases, const uint8_t* hap_bases, const double* row_params,
                               const nomics_pairhmm_pair* pairs, uint32_t n, double* out_sum) {
    const uint32_t warp_in_grid = (blockIdx.x * blockDim.x + threadIdx.x) / WARP;
    const int lane = threadIdx.x & (WARP - 1);
    if (warp_in_grid >= n) return;

    const nomics_pairhmm_pair pair = pairs[warp_in_grid];
    if (pair.read_len == 0 || pair.hap_len == 0 || pair.read_len > NOMICS_PAIRHMM_MAX_READ_LEN) {
        if (lane == 0) out_sum[warp_in_grid] = __longlong_as_double(0x7ff8000000000000ull);
        return;
    }

    // d1 and d2 are the preceding one and two anti-diagonals. Invalid cells
    // are explicitly zeroed, which also establishes the row/column boundaries.
    double m_d1[ROWS_PER_LANE] = {}, m_d2[ROWS_PER_LANE] = {};
    double i_d1[ROWS_PER_LANE] = {}, i_d2[ROWS_PER_LANE] = {};
    double d_d1[ROWS_PER_LANE] = {}, d_d2[ROWS_PER_LANE] = {};
    double final_sum = 0.0;
    const unsigned mask = 0xffffffffu;

    // Cell (i,j) is on diagonal i+j. (1,1) is diagonal 2.
    const uint32_t last_diagonal = pair.read_len + pair.hap_len;
    for (uint32_t diagonal = 2; diagonal <= last_diagonal; ++diagonal) {
#pragma unroll
        for (int k = ROWS_PER_LANE - 1; k >= 0; --k) {
            const uint32_t row = (uint32_t)(lane + 1 + k * WARP);
            const int64_t col_signed = (int64_t)diagonal - (int64_t)row;
            const bool valid = row <= pair.read_len && col_signed >= 1 && col_signed <= (int64_t)pair.hap_len;

            // Values for row i-1. Lane zero wraps to row group k-1. All lanes
            // execute both shuffles so every named source lane is active.
            const double m2_same = __shfl_up_sync(mask, m_d2[k], 1);
            const double i2_same = __shfl_up_sync(mask, i_d2[k], 1);
            const double d2_same = __shfl_up_sync(mask, d_d2[k], 1);
            const double m1_same = __shfl_up_sync(mask, m_d1[k], 1);
            const double i1_same = __shfl_up_sync(mask, i_d1[k], 1);
            const double m2_wrap = __shfl_sync(mask, k == 0 ? 0.0 : m_d2[k - 1], WARP - 1);
            const double i2_wrap = __shfl_sync(mask, k == 0 ? 0.0 : i_d2[k - 1], WARP - 1);
            const double d2_wrap = __shfl_sync(mask, k == 0 ? 0.0 : d_d2[k - 1], WARP - 1);
            const double m1_wrap = __shfl_sync(mask, k == 0 ? 0.0 : m_d1[k - 1], WARP - 1);
            const double i1_wrap = __shfl_sync(mask, k == 0 ? 0.0 : i_d1[k - 1], WARP - 1);

            double mp2 = lane == 0 ? m2_wrap : m2_same;
            double ip2 = lane == 0 ? i2_wrap : i2_same;
            double dp2 = lane == 0 ? d2_wrap : d2_same;
            const double mp1 = lane == 0 ? m1_wrap : m1_same;
            const double ip1 = lane == 0 ? i1_wrap : i1_same;

            // CPU row zero has M=I=0 and D=initial in every column, including
            // column zero. It is not represented in the lane arrays.
            if (row == 1) {
                mp2 = 0.0;
                ip2 = 0.0;
                dp2 = pair.initial;
            }

            double mv = 0.0, iv = 0.0, dv = 0.0;
            if (valid) {
                const uint32_t col = (uint32_t)col_signed;
                const double* p = row_params + (pair.params_offset + row - 1) * NOMICS_PAIRHMM_ROW_PARAMS;
                const uint8_t x = read_bases[pair.read_offset + row - 1];
                const uint8_t y = hap_bases[pair.hap_offset + col - 1];
                const double prior = (x == y || x == 'N' || y == 'N') ? p[5] : p[6];

                const double mt = add_rn(add_rn(mul_rn(mp2, p[0]), mul_rn(ip2, p[3])), mul_rn(dp2, p[3]));
                mv = mul_rn(prior, mt);
                iv = add_rn(mul_rn(mp1, p[1]), mul_rn(ip1, p[4]));
                dv = add_rn(mul_rn(m_d1[k], p[2]), mul_rn(d_d1[k], p[4]));

                if (row == pair.read_len) final_sum = add_rn(final_sum, add_rn(mv, iv));
            }
            // Descending row-group order leaves group k-1 untouched until its
            // wraparound values have been shuffled, avoiding three more local
            // arrays solely for the new diagonal.
            m_d2[k] = m_d1[k];
            i_d2[k] = i_d1[k];
            d_d2[k] = d_d1[k];
            m_d1[k] = mv;
            i_d1[k] = iv;
            d_d1[k] = dv;
        }
    }

    const int final_lane = (pair.read_len - 1) & (WARP - 1);
    if (lane == final_lane) out_sum[warp_in_grid] = final_sum;
}

extern "C" int nomics_pairhmm(const uint8_t* read_bases, const uint8_t* hap_bases, const double* row_params,
                               const nomics_pairhmm_pair* pairs, uint32_t n, double* out_sum, float* out_ms) {
    if (n == 0) {
        if (out_ms) *out_ms = 0.0f;
        return cudaSuccess;
    }

    cudaEvent_t t0 = nullptr, t1 = nullptr;
    cudaError_t error = cudaEventCreate(&t0);
    if (error == cudaSuccess) error = cudaEventCreate(&t1);
    if (error == cudaSuccess) error = cudaEventRecord(t0);
    if (error == cudaSuccess) {
        const int threads = 128;
        const int warps = threads / WARP;
        pairhmm_kernel<<<(n + warps - 1) / warps, threads>>>(read_bases, hap_bases, row_params, pairs, n, out_sum);
        error = cudaGetLastError();
    }
    if (error == cudaSuccess) error = cudaEventRecord(t1);
    if (error == cudaSuccess) error = cudaEventSynchronize(t1);
    if (error == cudaSuccess && out_ms) error = cudaEventElapsedTime(out_ms, t0, t1);
    // On any error after the launch, drain before returning: the caller falls
    // back to the CPU and may reuse or regrow these buffers right away.
    if (error != cudaSuccess) cudaDeviceSynchronize();
    if (t0) cudaEventDestroy(t0);
    if (t1) cudaEventDestroy(t1);
    return error;
}
