// C ABI shared by the CUDA shim (shim/*.cu), the Metal shim (shim/metal/) and
// the Rust bindings in src/lib.rs. Every function returns 0 on success; on
// failure a cudaError_t value (CUDA) or a 1001+ status (Metal, see
// shim/metal/nomics_metal.m).
#pragma once
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    char name[256];
    int major;
    int minor;
    int multiprocessors;
    int pageable_memory_access;  // GPU can dereference ordinary host pointers
    size_t total_global_mem;
    size_t max_buffer_bytes;  // largest single buffer the GPU can map (SIZE_MAX: no limit)
} nomics_device_info;

int nomics_device_info_get(int device, nomics_device_info* out);

// Sum `n` u64 words read in place from host memory `data` (no staging copy).
int nomics_sum_u64_in_place(const uint64_t* data, size_t n, uint64_t* out_sum, float* out_ms);

// ---- Streams and tickets (shim/stream.cu) ----
// Every `*_launch` entry point validates, enqueues its kernels on `stream`
// and returns a ticket; nomics_ticket_wait completes it. The buffers a launch
// reads or writes stay in use until its ticket is waited: the caller must not
// touch, reuse or free them before that. If a launch fails after enqueuing
// anything, it drains the stream before returning (no ticket), so on any
// error return the buffers are free again. A null ticket (nothing to launch,
// e.g. n == 0) is valid and waits as an immediate success.

// A cudaStream_t created with cudaStreamNonBlocking.
typedef void* nomics_stream;
// One launch's start and end events.
typedef struct nomics_ticket_s* nomics_ticket;

// priority 0: the least (default) priority; otherwise the greatest.
int nomics_stream_create(int priority, nomics_stream* out);
// Waits for the stream's work, then destroys it.
int nomics_stream_destroy(nomics_stream stream);
// Runs an empty kernel on the stream and waits, so the driver allocates the
// stream's context buffers now (from contiguous host memory on the GB10),
// before the index and prep fragment memory. A no-op where not needed.
int nomics_stream_warm(nomics_stream stream);
// Waits for the launch (its end event, so later launches on the same stream
// keep running), sets *out_ms (may be null) to the launch's GPU time, and
// frees the ticket. On error it first drains the launch's stream. Wait every
// ticket exactly once.
int nomics_ticket_wait(nomics_ticket ticket, float* out_ms);

// ---- Host buffer registration ----
// Long-lived host buffers that many launches read or write. The Metal shim
// wraps a registered range once (wrapping costs ~17 ms per GiB on an M4 Pro)
// and binds that wrapper for any launch pointer inside it; unregistered
// pointers are wrapped per launch. CUDA reads host memory in place, so both
// are no-ops there. [ptr, ptr + len) must be page aligned (HostBuf is).
// Unregister before freeing, and only after every launch that touched the
// range has been waited.
int nomics_host_register(const void* ptr, size_t len);
int nomics_host_unregister(const void* ptr);

// ---- Seed search (shim/seed.cu) ----

// bwa's BWT: header fields and the interleaved occurrence/base words,
// read in place from host memory.
typedef struct {
    uint64_t primary;
    uint64_t l2[5];
    uint64_t seq_len;
    const uint32_t* words;
} nomics_bwt;

// bwtintv_t: bi-interval x[0..3) and info = qbeg << 32 | qend.
typedef struct {
    uint64_t x[3];
    uint64_t info;
} nomics_intv;

#define NOMICS_KMER_MAX_K 13
typedef struct {
    uint32_t k;
    uint64_t level_offsets[NOMICS_KMER_MAX_K + 1];
    const uint64_t* intervals;  // three u64s per endpoint
} nomics_kmer_table;

#define NOMICS_SEED_MAX_LEN 256
#define NOMICS_SEED_OVERFLOW 0xffffffffu
#define NOMICS_SEED_WORK_64 64u
#define NOMICS_SEED_WORK_96 96u
#define NOMICS_SEED_WORK_128 128u
#define NOMICS_SEED_WORK_TWO_PASS 192u

// mem_collect_intv for `n` reads (2-bit codes, read r at seqs[offs[r]..+lens[r]]).
// Read r's sorted intervals go to out[r*cap ..], their count to out_n[r];
// NOMICS_SEED_OVERFLOW means "more than cap, or longer than
// NOMICS_SEED_MAX_LEN": recompute on the CPU. All pointers are host memory.
// `table` may be null. Work modes select a 64-, 96-, or 128-entry local list;
// TWO_PASS first runs 64 and retries only overflow reads with 128.
//
// `scratch` (may be null) holds the kernel's per-thread interval lists, one
// buffer per list size (index 0: 64, 1: 96, 2: 128), each at least
// nomics_seed_scratch_bytes(size) bytes, in host memory (huge pages) that no
// other code touches during the call; concurrent calls may share it. A null
// buffer for a size the mode uses selects the local-memory kernel instead
// (same output, about half the speed).
typedef struct {
    nomics_intv* lists[3];
    size_t bytes[3];
} nomics_seed_scratch;
int nomics_seed_scratch_bytes(uint32_t list_size, size_t* bytes);
//
// `bwt` and `table` are copied at launch, but the BWT words, the k-mer
// intervals, seqs/offs/lens, out/out_n and the scratch are read or written
// until the ticket is waited.
int nomics_seed_collect_launch(nomics_stream stream, const nomics_bwt* bwt,
                               const nomics_kmer_table* table, const uint8_t* seqs,
                               const uint64_t* offs, const uint32_t* lens, uint32_t n, uint32_t cap,
                               uint32_t work_mode, nomics_intv* out, uint32_t* out_n,
                               const nomics_seed_scratch* scratch, nomics_ticket* ticket);

// ---- PairHMM experiment (shim/pairhmm.cu) ----

#define NOMICS_PAIRHMM_ROW_PARAMS 7
#define NOMICS_PAIRHMM_MAX_READ_LEN 160

// Base offsets count bytes/bases; params_offset counts seven-double rows.
// `initial` is 2^1020 / hap_len, computed on the CPU.
typedef struct {
    uint64_t read_offset;
    uint64_t hap_offset;
    uint64_t params_offset;
    uint32_t read_len;
    uint32_t hap_len;
    double initial;
} nomics_pairhmm_pair;

// One warp computes each pair. Row parameters are seven doubles per read row:
// m2m, m2i, m2d, i2m, i2i (= d2d), p_match, p_mismatch. `out_sum` receives
// the ordered sum of M+I on the final row; log10 remains a CPU operation.
int nomics_pairhmm(const uint8_t* read_bases, const uint8_t* hap_bases, const double* row_params,
                   const nomics_pairhmm_pair* pairs, uint32_t n, double* out_sum, float* out_ms);

// ---- BWA ksw_align2 experiment (shim/ksw.cu) ----

// Offsets into the packed host buffers. `scratch_offset` counts uint16_t
// elements and reserves four striped-query arrays; `row_offset` counts
// uint16_t target-row maxima. Both are assigned by the Rust wrapper.
typedef struct {
    uint64_t query_offset;
    uint64_t target_offset;
    uint64_t mat_offset;
    uint64_t scratch_offset;
    uint64_t row_offset;
    uint32_t query_len;
    uint32_t target_len;
    int32_t m;
    int32_t o_del;
    int32_t e_del;
    int32_t o_ins;
    int32_t e_ins;
    int32_t xtra;
    int32_t size;
} nomics_ksw_job;

typedef struct {
    int32_t score;
    int32_t te;
    int32_t qe;
    int32_t score2;
    int32_t te2;
    int32_t tb;
    int32_t qb;
} nomics_ksw_result;

// One CUDA thread evaluates one alignment in the same striped order as
// BWA 0.7.19. All pointers refer to huge-page-backed host memory. `scratch`
// and `row_maxima` are workspace owned by the caller; there is no device-side
// allocation or registration.
int nomics_ksw_align2_launch(nomics_stream stream, const uint8_t* queries, const uint8_t* targets,
                             const int8_t* matrices, const nomics_ksw_job* jobs, uint32_t n,
                             uint16_t* scratch, uint16_t* row_maxima, nomics_ksw_result* out,
                             nomics_ticket* ticket);

// The same results from a wavefront over the padded query (see ksw.cu): one
// warp per alignment, or (CUDA) two byte-sized alignments per warp in 16-bit
// halves (NOMICS_KSW_PAIR=0 turns that off); ~10x the striped kernel above on
// the GB10 before pairing. Targets and row maxima are staged in shared memory
// (targets up to ~4 kb); `row_maxima` is unused and kept for the ABI (callers
// pass NULL; implementations must not read or write it). Negative gap
// penalties are refused (cudaErrorInvalidValue), as are m > 8.
int nomics_ksw_align2_wave_launch(nomics_stream stream, const uint8_t* queries,
                                  const uint8_t* targets, const int8_t* matrices,
                                  const nomics_ksw_job* jobs, uint32_t n, uint16_t* row_maxima,
                                  nomics_ksw_result* out, nomics_ticket* ticket);

// One BWA ksw_extend2 call (see extend.cu). Offsets index the packed host
// buffers; all pointers refer to huge-page-backed host memory.
typedef struct {
    uint64_t query_offset;
    uint64_t target_offset;
    uint64_t mat_offset;
    uint32_t query_len;
    uint32_t target_len;
    int32_t m;
    int32_t o_del;
    int32_t e_del;
    int32_t o_ins;
    int32_t e_ins;
    int32_t w;
    int32_t end_bonus;
    int32_t zdrop;
    int32_t h0;
} nomics_extend_job;

typedef struct {
    int32_t score;
    int32_t qle;
    int32_t tle;
    int32_t gtle;
    int32_t gscore;
    int32_t max_off;
} nomics_extend_result;

// Queries up to 256 bases; every job needs query_len > 0 and m*m <= 64.
// Callers may launch concurrently with disjoint buffers and streams.
int nomics_ksw_extend2_launch(nomics_stream stream, const uint8_t* queries, const uint8_t* targets,
                              const int8_t* matrices, const nomics_extend_job* jobs, uint32_t n,
                              nomics_extend_result* out, nomics_ticket* ticket);

// ---- Fused left+right seed extension with band retries (shim/fused.cu) ----
// Design: docs/strategy/gpu-extend-fused-design.md.

// Scoring shared by every job of a launch (BWA mem defaults today).
typedef struct {
    int32_t m;                 // alphabet size, m*m <= 64
    int32_t o_del, e_del, o_ins, e_ins;
    int32_t w;                 // first band, 100; band i is w << i
    int32_t max_band_try;      // 2
    int32_t end_bonus_left;    // PEN_CLIP5
    int32_t end_bonus_right;   // PEN_CLIP3
    int32_t zdrop;
    int8_t  mat[64];           // m*m used
} nomics_fused_params;         // 104 bytes

typedef struct {
    uint64_t lq_off, lt_off;   // left query (reversed prefix) and target (reversed window)
    uint64_t rq_off, rt_off;   // right query (suffix) and target
    uint32_t lq_len, lt_len;   // lq_len == 0: no left (qbeg == 0)
    uint32_t rq_len, rt_len;   // rq_len == 0: no right on the GPU
    int32_t  h0;               // seed len * A
    int32_t  pad;
} nomics_fused_job;            // 56 bytes

typedef struct {
    nomics_extend_result left[2], right[2];
    int32_t n_left, n_right;   // calls made per side, 0..2
    int32_t right_h0;          // sc0
    int32_t pad;
} nomics_fused_result;         // 112 bytes

// Queries up to 256 bases per side; lq_len + rq_len > 0. Same contract as
// nomics_ksw_extend2_launch; `params` is copied at launch. Offsets index `queries`
// (lq, rq) and `targets` (lt, rt). Call i on a side used band w << i; left
// calls used h0 = job.h0, right calls right_h0. `order` is caller-owned
// workspace of n entries (huge-page host memory, read by the GPU in place):
// the launcher writes the jobs' length-binned order there. Results are in
// input order. The stream's priority is the caller's choice (nomics-gpu
// gives fused engines the greatest unless NOMICS_FUSED_PRIORITY=0).
int nomics_fused_extend_launch(nomics_stream stream, const uint8_t* queries, const uint8_t* targets,
                               const nomics_fused_params* params, const nomics_fused_job* jobs,
                               uint32_t n, uint32_t* order, nomics_fused_result* out,
                               nomics_ticket* ticket);

#ifdef __cplusplus
}
#endif
