// Streams and tickets: the asynchronous launch ABI (nomics_gpu.h).
#include "ticket.cuh"

#include <new>

extern "C" int nomics_stream_create(int priority, nomics_stream* out) {
    if (out == nullptr) return cudaErrorInvalidValue;
    *out = nullptr;
    int least = 0, greatest = 0;
    cudaError_t e = cudaDeviceGetStreamPriorityRange(&least, &greatest);
    if (e != cudaSuccess) return e;
    cudaStream_t stream = nullptr;
    // Nonblocking: no implicit synchronization with the legacy default stream.
    e = cudaStreamCreateWithPriority(&stream, cudaStreamNonBlocking, priority ? greatest : least);
    if (e != cudaSuccess) return e;
    *out = (nomics_stream)stream;
    return cudaSuccess;
}

__global__ void nomics_warm_kernel() {}

extern "C" int nomics_stream_warm(nomics_stream stream) {
    if (stream == nullptr) return cudaErrorInvalidValue;
    cudaStream_t s = (cudaStream_t)stream;
    nomics_warm_kernel<<<1, 1, 0, s>>>();
    cudaError_t e = cudaGetLastError();
    const cudaError_t synced = cudaStreamSynchronize(s);
    return e != cudaSuccess ? e : synced;
}

extern "C" int nomics_stream_destroy(nomics_stream stream) {
    if (stream == nullptr) return cudaSuccess;
    cudaStream_t s = (cudaStream_t)stream;
    const cudaError_t synced = cudaStreamSynchronize(s);
    const cudaError_t destroyed = cudaStreamDestroy(s);
    return synced != cudaSuccess ? synced : destroyed;
}

cudaError_t ticket_begin(cudaStream_t stream, nomics_ticket_s** out) {
    *out = nullptr;
    nomics_ticket_s* t = new (std::nothrow) nomics_ticket_s{stream, nullptr, nullptr};
    if (t == nullptr) return cudaErrorMemoryAllocation;
    cudaError_t e = cudaEventCreate(&t->t0);
    if (e == cudaSuccess) e = cudaEventCreate(&t->t1);
    if (e == cudaSuccess) e = cudaEventRecord(t->t0, stream);
    if (e != cudaSuccess) {
        // Nothing of this launch is enqueued (recording t0 failed or never ran).
        if (t->t0) cudaEventDestroy(t->t0);
        if (t->t1) cudaEventDestroy(t->t1);
        delete t;
        return e;
    }
    *out = t;
    return cudaSuccess;
}

static void ticket_free(nomics_ticket_s* t) {
    cudaEventDestroy(t->t0);
    cudaEventDestroy(t->t1);
    delete t;
}

cudaError_t ticket_finish(nomics_ticket_s* t, cudaError_t pending, nomics_ticket* out) {
    *out = nullptr;
    cudaError_t e = pending;
    if (e == cudaSuccess) e = cudaEventRecord(t->t1, t->stream);
    if (e != cudaSuccess) {
        // Earlier enqueues of this launch may be running: finish them before
        // the caller reuses or frees the buffers they touch.
        cudaStreamSynchronize(t->stream);
        ticket_free(t);
        return e;
    }
    *out = (nomics_ticket)t;
    return cudaSuccess;
}

extern "C" int nomics_ticket_wait(nomics_ticket ticket, float* out_ms) {
    if (out_ms) *out_ms = 0.0f;
    if (ticket == nullptr) return cudaSuccess;  // nothing was launched
    nomics_ticket_s* t = (nomics_ticket_s*)ticket;
    // t1 follows every enqueue of the launch on its stream. Waiting on the
    // event, not the stream, leaves later launches on the stream running.
    cudaError_t e = cudaEventSynchronize(t->t1);
    float ms = 0.0f;
    if (e == cudaSuccess) e = cudaEventElapsedTime(&ms, t->t0, t->t1);
    // On failure, drain the stream before the caller touches the buffers.
    if (e != cudaSuccess) cudaStreamSynchronize(t->stream);
    if (e == cudaSuccess && out_ms) *out_ms = ms;
    ticket_free(t);
    return e;
}

// CUDA dereferences host memory in place: registration is a no-op.
extern "C" int nomics_host_register(const void*, size_t) { return cudaSuccess; }
extern "C" int nomics_host_unregister(const void*) { return cudaSuccess; }
