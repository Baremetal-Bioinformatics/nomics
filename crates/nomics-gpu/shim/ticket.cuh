// Internal helpers for the asynchronous launch ABI (nomics_gpu.h, "Streams
// and tickets"). Defined in stream.cu.
#pragma once
#include "nomics_gpu.h"

#include <cuda_runtime.h>

// One launch: its stream (owned by the caller) and its start/end events.
struct nomics_ticket_s {
    cudaStream_t stream;
    cudaEvent_t t0, t1;
};

// Creates the ticket's events and records t0 on `stream`. On error nothing
// has been enqueued and *out is null.
cudaError_t ticket_begin(cudaStream_t stream, nomics_ticket_s** out);

// Ends a launch whose enqueues returned `pending`. On success records t1 and
// hands the ticket out. On any error (pending, or recording t1) drains the
// stream so the caller may reuse its buffers, frees the ticket, sets *out to
// null and returns the error.
cudaError_t ticket_finish(nomics_ticket_s* t, cudaError_t pending, nomics_ticket* out);
