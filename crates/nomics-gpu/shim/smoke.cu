// Smoke kernels proving the unified-memory premise: the GPU reads a 2 MiB-page
// host buffer in place. Not part of any pipeline stage.
#include "nomics_gpu.h"

#include <cub/cub.cuh>
#include <cuda_runtime.h>
#include <string.h>

extern "C" int nomics_device_info_get(int device, nomics_device_info* out) {
    cudaDeviceProp p;
    cudaError_t e = cudaGetDeviceProperties(&p, device);
    if (e != cudaSuccess) return e;
    memset(out, 0, sizeof(*out));
    strncpy(out->name, p.name, sizeof(out->name) - 1);
    out->major = p.major;
    out->minor = p.minor;
    out->multiprocessors = p.multiProcessorCount;
    out->pageable_memory_access = p.pageableMemoryAccess;
    out->total_global_mem = p.totalGlobalMem;
    out->max_buffer_bytes = SIZE_MAX;
    return cudaSuccess;
}

extern "C" int nomics_sum_u64_in_place(const uint64_t* data, size_t n, uint64_t* out_sum, float* out_ms) {
    // Small device-side scalars only; the input stays in host memory.
    uint64_t* d_sum = nullptr;
    void* d_tmp = nullptr;
    size_t tmp_bytes = 0;
    cudaEvent_t t0 = nullptr, t1 = nullptr;
    cudaError_t e = cub::DeviceReduce::Sum(nullptr, tmp_bytes, data, (uint64_t*)nullptr, n);
    if (e == cudaSuccess) e = cudaMallocAsync(&d_tmp, tmp_bytes, 0);
    if (e == cudaSuccess) e = cudaMallocAsync((void**)&d_sum, sizeof(uint64_t), 0);
    if (e == cudaSuccess) { cudaEventCreate(&t0); cudaEventCreate(&t1); cudaEventRecord(t0); }
    if (e == cudaSuccess) e = cub::DeviceReduce::Sum(d_tmp, tmp_bytes, data, d_sum, n);
    if (e == cudaSuccess) { cudaEventRecord(t1); e = cudaEventSynchronize(t1); }
    if (e == cudaSuccess) { cudaEventElapsedTime(out_ms, t0, t1); e = cudaMemcpy(out_sum, d_sum, sizeof(uint64_t), cudaMemcpyDeviceToHost); }
    if (d_tmp) cudaFreeAsync(d_tmp, 0);
    if (d_sum) cudaFreeAsync(d_sum, 0);
    if (t0) cudaEventDestroy(t0);
    if (t1) cudaEventDestroy(t1);
    cudaStreamSynchronize(0);
    return e;
}
