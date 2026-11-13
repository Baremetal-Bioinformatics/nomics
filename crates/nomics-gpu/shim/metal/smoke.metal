// Smoke kernel proving the unified-memory premise on Apple Silicon: the GPU
// reads a host buffer in place (wrapped with newBufferWithBytesNoCopy). Not
// part of any pipeline stage. Mirrors shim/smoke.cu.
#include <metal_stdlib>
using namespace metal;

// Each threadgroup writes one partial sum; the host adds the partials.
// Integer addition wraps, so the result equals the CPU's wrapping sum
// regardless of the order the partials are combined in.
kernel void sum_u64(device const ulong* data [[buffer(0)]],
                    device ulong* partial [[buffer(1)]],
                    constant ulong& n [[buffer(2)]],
                    threadgroup ulong* scratch [[threadgroup(0)]],
                    uint gid [[thread_position_in_grid]],
                    uint grid [[threads_per_grid]],
                    uint lid [[thread_index_in_threadgroup]],
                    uint group [[threadgroup_position_in_grid]],
                    uint width [[threads_per_threadgroup]]) {
    ulong sum = 0;
    for (ulong i = gid; i < n; i += grid) {
        sum += data[i];
    }
    scratch[lid] = sum;
    threadgroup_barrier(mem_flags::mem_threadgroup);
    for (uint stride = width / 2; stride > 0; stride /= 2) {
        if (lid < stride) {
            scratch[lid] += scratch[lid + stride];
        }
        threadgroup_barrier(mem_flags::mem_threadgroup);
    }
    if (lid == 0) {
        partial[group] = scratch[0];
    }
}
