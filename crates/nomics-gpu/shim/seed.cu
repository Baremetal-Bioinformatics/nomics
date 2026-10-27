#include <cstdio>
#include <cstdlib>
#include <cstring>
// Seed search on the GPU: bwa 0.7.19's mem_collect_intv (SMEMs, re-seeding,
// LAST-like seeds) for a batch of reads, one thread per read.
//
// A line-for-line port of crates/nomics-align/src/seed.rs (itself gated
// against bwa on every read); pure integer code, so results are identical by
// construction and gated again against the CPU path. The BWT, the reads, the
// outputs and the per-thread interval lists (collect_slot_kernel) are read and
// written in place in host memory.
#include "nomics_gpu.h"
#include "ticket.cuh"

#include <cuda_runtime.h>

#include <atomic>

#define OCC_INTV_SHIFT 7
#define OCC_INTV_MASK 127ull
#define MIN_SEED_LEN 19
#define SPLIT_WIDTH 10
#define MAX_MEM_INTV 20
#define MAXL NOMICS_SEED_MAX_LEN
#define SMEM_OVERFLOW (-2)

typedef nomics_intv Intv;


struct Ctx {
    nomics_bwt b;
};

// Occurrences of each base among the 16 2-bit bases of `w`, packed one count
// per byte (bwa's cnt_table sum), via popcount instead of a table: divergent
// indices into __constant__ memory would serialize the warp.
__device__ __forceinline__ uint32_t occ_aux4(uint32_t w) {
    uint32_t hi = w >> 1 & 0x55555555u, lo = w & 0x55555555u;
    uint32_t a = __popc(~hi & ~lo & 0x55555555u);
    uint32_t c = __popc(~hi & lo & 0x55555555u);
    uint32_t g = __popc(hi & ~lo & 0x55555555u);
    uint32_t t = __popc(hi & lo);
    return a | c << 8 | g << 16 | t << 24;
}

// bwt_occ4
__device__ void occ4(const nomics_bwt& b, uint64_t k, uint64_t cnt[4]) {
    if (k == ~0ull) {
        cnt[0] = cnt[1] = cnt[2] = cnt[3] = 0;
        return;
    }
    k -= (k >= b.primary);
    const uint32_t* p = b.words + ((k >> OCC_INTV_SHIFT) << 4);
    const uint64_t* c = (const uint64_t*)p;
    cnt[0] = c[0];
    cnt[1] = c[1];
    cnt[2] = c[2];
    cnt[3] = c[3];
    p += 8;
    const uint32_t* end = p + ((k >> 4) - ((k & ~OCC_INTV_MASK) >> 4));
    uint64_t x = 0;
    for (; p < end; ++p) x += occ_aux4(*p);
    uint32_t tmp = *p & ~((1u << ((~k & 15) << 1)) - 1);
    x += occ_aux4(tmp) - (~k & 15);
    cnt[0] += x & 0xff;
    cnt[1] += x >> 8 & 0xff;
    cnt[2] += x >> 16 & 0xff;
    cnt[3] += x >> 24;
}

// bwt_set_intv
__device__ __forceinline__ Intv set_intv(const nomics_bwt& b, int c) {
    Intv r;
    r.x[0] = b.l2[c] + 1;
    r.x[1] = b.l2[3 - c] + 1;
    r.x[2] = b.l2[c + 1] - b.l2[c];
    r.info = 0;
    return r;
}

__device__ __forceinline__ Intv table_intv(const nomics_kmer_table& table, int len, uint64_t code) {
    const uint64_t* p = table.intervals + 3 * (table.level_offsets[len] + code);
    Intv r;
    r.x[0] = p[0];
    r.x[1] = p[1];
    r.x[2] = p[2];
    r.info = 0;
    return r;
}

// bwt_extend: is_back selects backward (prepend) or forward extension.
__device__ void extend(const nomics_bwt& b, const Intv& ik, Intv ok[4], bool is_back) {
    int bi = is_back ? 1 : 0, nb = is_back ? 0 : 1;
    uint64_t k = ik.x[nb] - 1;
    uint64_t tk[4], tl[4];
    occ4(b, k, tk);
    occ4(b, k + ik.x[2], tl);
    for (int i = 0; i < 4; ++i) {
        ok[i].x[nb] = b.l2[i] + 1 + tk[i];
        ok[i].x[2] = tl[i] - tk[i];
        ok[i].info = 0;
    }
    uint64_t spans = (ik.x[nb] <= b.primary && ik.x[nb] + ik.x[2] - 1 >= b.primary) ? 1 : 0;
    ok[3].x[bi] = ik.x[bi] + spans;
    ok[2].x[bi] = ok[3].x[bi] + ok[3].x[2];
    ok[1].x[bi] = ok[2].x[bi] + ok[2].x[2];
    ok[0].x[bi] = ok[1].x[bi] + ok[1].x[2];
}

// Interval lists for smem1. PlainList stores bwtintv_t as is. PackedList
// stores one in 16 bytes: x[0], x[1], x[2] < 2^36 and qbeg, qend < 2^10
// (reads are at most NOMICS_SEED_MAX_LEN = 256 bases). put() refuses a value
// outside those bounds, and smem1 then reports SMEM_OVERFLOW, so the read is
// recomputed on the CPU: the packing can never change a result.
struct PlainList {
    Intv* p;
    __device__ __forceinline__ Intv get(int i) const { return p[i]; }
    __device__ __forceinline__ bool put(int i, const Intv& v) const {
        p[i] = v;
        return true;
    }
};

struct PackedList {
    uint4* p;
    __device__ __forceinline__ Intv get(int i) const {
        uint4 w = p[i];
        Intv v;
        v.x[0] = (uint64_t)w.x | (uint64_t)(w.w & 0xf) << 32;
        v.x[1] = (uint64_t)w.y | (uint64_t)(w.w >> 4 & 0xf) << 32;
        v.x[2] = (uint64_t)w.z | (uint64_t)(w.w >> 8 & 0xf) << 32;
        v.info = (uint64_t)(w.w >> 22) << 32 | (uint64_t)(w.w >> 12 & 0x3ff);
        return v;
    }
    __device__ __forceinline__ bool put(int i, const Intv& v) const {
        uint64_t qbeg = v.info >> 32, qend = (uint32_t)v.info;
        if (((v.x[0] | v.x[1] | v.x[2]) >> 36) != 0 || qbeg >= 1024 || qend >= 1024) return false;
        uint4 w;
        w.x = (uint32_t)v.x[0];
        w.y = (uint32_t)v.x[1];
        w.z = (uint32_t)v.x[2];
        w.w = (uint32_t)(v.x[0] >> 32) | (uint32_t)(v.x[1] >> 32) << 4 | (uint32_t)(v.x[2] >> 32) << 8 |
              (uint32_t)qend << 12 | (uint32_t)qbeg << 22;
        p[i] = w;
        return true;
    }
};

// bwt_smem1 (bwt_smem1a with max_intv = 0). `mem` receives the SMEMs covering
// x, in reverse order of start (see the end); returns where the next search
// starts, or SMEM_OVERFLOW if a list outgrew MAXI.
template<int MAXI, class L>
__device__ int smem1(const nomics_bwt& b, const nomics_kmer_table& table, bool use_table,
                     const uint8_t* q, int len, int x, uint64_t min_intv, L mem, int* n_mem,
                     L prev, L curr) {
    *n_mem = 0;
    if (q[x] > 3) return x + 1;
    if (min_intv < 1) min_intv = 1;
    Intv ik = set_intv(b, q[x]), ok[4];
    ik.info = x + 1;
    int n_curr = 0, n_prev;
    int i = x + 1;
    bool stopped = false;
    if (use_table && x + (int)table.k <= len) {
        bool canonical = true;
        for (int j = 0; j < (int)table.k; ++j) canonical &= q[x + j] < 4;
        if (canonical) {
            uint64_t code = q[x];
            for (int prefix_len = 2; prefix_len <= (int)table.k; ++prefix_len) {
                code = code << 2 | q[x + prefix_len - 1];
                Intv next = table_intv(table, prefix_len, code);
                if (next.x[2] != ik.x[2]) {
                    if (n_curr == MAXI || !curr.put(n_curr++, ik)) return SMEM_OVERFLOW;
                    if (next.x[2] < min_intv) {
                        i = x + prefix_len - 1;
                        stopped = true;
                        break;
                    }
                }
                ik = next;
                ik.info = x + prefix_len;
            }
            if (!stopped) i = x + table.k;
        }
    }
    for (; !stopped && i < len; ++i) {
        if (q[i] < 4) {
            int c = 3 - q[i];
            extend(b, ik, ok, false);
            if (ok[c].x[2] != ik.x[2]) {
                if (n_curr == MAXI || !curr.put(n_curr++, ik)) return SMEM_OVERFLOW;
                if (ok[c].x[2] < min_intv) break;
            }
            ik = ok[c];
            ik.info = i + 1;
        } else {
            if (n_curr == MAXI || !curr.put(n_curr++, ik)) return SMEM_OVERFLOW;
            break;
        }
    }
    if (i == len) {
        if (n_curr == MAXI || !curr.put(n_curr++, ik)) return SMEM_OVERFLOW;
    }
    // bwa reverses curr here (longest match first). Instead, the first
    // backward step below reads it from the end.
    int ret = (int)(uint32_t)curr.get(n_curr - 1).info;
    L t = prev;
    prev = curr;
    curr = t;
    n_prev = n_curr;

    // The list reads bwa makes for the last curr and mem entries are served
    // from registers: a successful put() stores its value exactly, so these
    // are the values the lists hold.
    uint64_t last_mem_beg = 0;
    bool from_end = true;
    for (i = x - 1; i >= -1; --i) {
        int c = (i < 0 || q[i] > 3) ? -1 : q[i];
        n_curr = 0;
        uint64_t last_curr_size = 0;
        for (int j = 0; j < n_prev; ++j) {
            const Intv p = prev.get(from_end ? n_prev - 1 - j : j);
            bool keep_ext = false;
            if (c >= 0) {
                extend(b, p, ok, true);
                keep_ext = ok[c].x[2] >= min_intv;
            }
            if (keep_ext) {
                if (n_curr == 0 || ok[c].x[2] != last_curr_size) {
                    ok[c].info = p.info;
                    if (n_curr == MAXI || !curr.put(n_curr++, ok[c])) return SMEM_OVERFLOW;
                    last_curr_size = ok[c].x[2];
                }
            } else if (n_curr == 0 && (*n_mem == 0 || (uint64_t)(i + 1) < last_mem_beg)) {
                Intv m = p;
                m.info |= (uint64_t)(i + 1) << 32;
                if (*n_mem == MAXI || !mem.put((*n_mem)++, m)) return SMEM_OVERFLOW;
                last_mem_beg = m.info >> 32;
            }
        }
        if (n_curr == 0) break;
        t = prev;
        prev = curr;
        curr = t;
        n_prev = n_curr;
        from_end = false;
    }
    // bwa reverses mem here to sort it by start. The only caller copies mem
    // into the read's intervals in order, and those are sorted by info at
    // the end (equal keys are identical intervals), so the order is not
    // observable and the reversal is skipped.
    return ret;
}

// bwt_seed_strategy1
__device__ int seed_strategy1(const nomics_bwt& b, const nomics_kmer_table& table, bool use_table,
                              const uint8_t* q, int len, int x, int min_len, uint64_t max_intv, Intv* m) {
    m->x[0] = m->x[1] = m->x[2] = m->info = 0;
    if (q[x] > 3) return x + 1;
    Intv ik = set_intv(b, q[x]), ok[4];
    int start = x + 1;
    if (use_table && table.k <= (uint32_t)min_len && x + (int)table.k <= len) {
        uint64_t code = 0;
        bool canonical = true;
        for (int j = 0; j < (int)table.k; ++j) {
            canonical &= q[x + j] < 4;
            code = code << 2 | q[x + j];
        }
        if (canonical) {
            ik = table_intv(table, table.k, code);
            start = x + table.k;
        }
    }
    for (int i = start; i < len; ++i) {
        if (q[i] < 4) {
            int c = 3 - q[i];
            extend(b, ik, ok, false);
            if (ok[c].x[2] < max_intv && i - x >= min_len) {
                *m = ok[c];
                m->info = (uint64_t)x << 32 | (uint64_t)(i + 1);
                return i + 1;
            }
            ik = ok[c];
        } else
            return i + 1;
    }
    return len;
}

// mem_collect_intv for read r; mem1, prev and curr are MAXI-entry lists.
template<int MAXI, class L>
__device__ __forceinline__ void collect_read(const nomics_bwt& b, const nomics_kmer_table& table, bool use_table,
                                             const uint8_t* seqs, const uint64_t* offs, const uint32_t* lens,
                                             uint32_t cap, bool retry_only, Intv* out, uint32_t* out_n, uint32_t r,
                                             L mem1, L prev, L curr) {
    if (retry_only && out_n[r] != NOMICS_SEED_OVERFLOW) return;
    // Compare before narrowing: a length above INT_MAX must overflow too.
    if (lens[r] > (uint32_t)MAXL) {
        out_n[r] = NOMICS_SEED_OVERFLOW;
        return;
    }
    int len = (int)lens[r];
    if (len > MAXL) {
        out_n[r] = NOMICS_SEED_OVERFLOW;
        return;
    }
    const uint8_t* q = seqs + offs[r];
    Intv* mem = out + (uint64_t)r * cap;
    uint32_t n_mem = 0;
    int n1;
    bool overflow = false;
    // Pass 1: SMEMs.
    int x = 0;
    while (x < len && !overflow) {
        if (q[x] < 4) {
            x = smem1<MAXI, L>(b, table, use_table, q, len, x, 1, mem1, &n1, prev, curr);
            if (x == SMEM_OVERFLOW) { overflow = true; break; }
            for (int i = 0; i < n1; ++i) {
                Intv m = mem1.get(i);
                uint32_t sl = (uint32_t)m.info - (uint32_t)(m.info >> 32);
                if (sl >= MIN_SEED_LEN) {
                    if (n_mem == cap) { overflow = true; break; }
                    mem[n_mem++] = m;
                }
            }
        } else
            ++x;
    }
    // Pass 2: re-seed long SMEMs with few occurrences from their midpoint.
    // This is the seed kernel's only floating-point arithmetic. It evaluates
    // to 28 with or without contraction, so the shim-wide --fmad=false used by
    // PairHMM cannot change a seed boundary.
    int split_len = (int)((float)MIN_SEED_LEN * 1.5f + .499);
    uint32_t old_n = n_mem;
    for (uint32_t k = 0; k < old_n && !overflow; ++k) {
        Intv p = mem[k];
        int start = (int)(p.info >> 32), end = (int)(uint32_t)p.info;
        if (end - start < split_len || p.x[2] > SPLIT_WIDTH) continue;
        if (smem1<MAXI, L>(b, table, use_table, q, len, (start + end) >> 1, p.x[2] + 1,
                        mem1, &n1, prev, curr) == SMEM_OVERFLOW) {
            overflow = true;
            break;
        }
        for (int i = 0; i < n1; ++i) {
            Intv m = mem1.get(i);
            uint32_t sl = (uint32_t)m.info - (uint32_t)(m.info >> 32);
            if (sl >= MIN_SEED_LEN) {
                if (n_mem == cap) { overflow = true; break; }
                mem[n_mem++] = m;
            }
        }
    }
    // Pass 3: LAST-like seeds.
    x = 0;
    while (x < len && !overflow) {
        if (q[x] < 4) {
            Intv m;
            x = seed_strategy1(b, table, use_table, q, len, x, MIN_SEED_LEN, MAX_MEM_INTV, &m);
            if (m.x[2] > 0) {
                if (n_mem == cap) { overflow = true; break; }
                mem[n_mem++] = m;
            }
        } else
            ++x;
    }
    if (overflow) {
        out_n[r] = NOMICS_SEED_OVERFLOW;
        return;
    }
    // Sort by info. Equal keys are identical intervals, so any sort is bwa's.
    for (uint32_t i = 1; i < n_mem; ++i) {
        Intv t = mem[i];
        uint32_t j = i;
        while (j > 0 && mem[j - 1].info > t.info) {
            mem[j] = mem[j - 1];
            --j;
        }
        mem[j] = t;
    }
    out_n[r] = n_mem;
}

#define SEED_THREADS 128
// Minimum resident blocks per SM for collect_slot_kernel (caps registers;
// 0 leaves the register count to the compiler).
#define SEED_MIN_BLOCKS 0
#if SEED_MIN_BLOCKS
#define SEED_SLOT_BOUNDS __launch_bounds__(SEED_THREADS, SEED_MIN_BLOCKS)
#else
#define SEED_SLOT_BOUNDS __launch_bounds__(SEED_THREADS)
#endif

// Lists in local memory (the fallback when no scratch is given). Local memory
// is interleaved across a warp at word granularity, so when lanes touch
// different list entries, which they nearly always do, each 32-byte interval
// costs eight separate sectors. That traffic, not the FM-index lookups, was
// the kernel's bound (see collect_slot_kernel).
template<int MAXI>
__global__ void __launch_bounds__(SEED_THREADS)
    collect_kernel(nomics_bwt b, nomics_kmer_table table, bool use_table, const uint8_t* seqs, const uint64_t* offs,
                   const uint32_t* lens, uint32_t n, uint32_t cap, bool retry_only, Intv* out, uint32_t* out_n) {
    uint32_t r = blockIdx.x * blockDim.x + threadIdx.x;
    if (r >= n) return;
    Intv mem1[MAXI], prev[MAXI], curr[MAXI];
    collect_read<MAXI, PlainList>(b, table, use_table, seqs, offs, lens, cap, retry_only, out, out_n, r,
                                  PlainList{mem1}, PlainList{prev}, PlainList{curr});
}

// Slot pool for the scratch lists: one flag per resident-block slot. A block
// claims a free slot at start and releases it at exit, so the scratch is
// sized by how many blocks can be resident at once, not by the grid.
#define SEED_MAX_SLOTS 4096
// The slot flags are per device while the Rust scratch is per process: one
// CUDA device per process (as everywhere in nomics-gpu).
__device__ unsigned int g_seed_slots[3][SEED_MAX_SLOTS];

// The same per-read code with each thread's three lists contiguous in a
// scratch area (host memory, read in place like everything else), so an
// interval is one sector (PlainList, 32 bytes) or half of one (PackedList).
// On a 30x batch (662,000 reads): collect_kernel 1,525 ms, plain lists
// 629 ms, packed 588 ms, identical output. Blocks stay one read per thread,
// so a higher-priority stream still gets SMs as blocks retire.
template<int MAXI, class L>
__global__ void SEED_SLOT_BOUNDS
    collect_slot_kernel(nomics_bwt b, nomics_kmer_table table, bool use_table, const uint8_t* seqs,
                        const uint64_t* offs, const uint32_t* lens, uint32_t n, uint32_t cap, bool retry_only,
                        Intv* out, uint32_t* out_n, void* scratch, unsigned int* slots, uint32_t n_slots) {
    __shared__ uint32_t slot;
    if (threadIdx.x == 0) {
        // At most n_slots blocks of this kernel are resident, so a slot is
        // free or about to be freed by a running block.
        uint32_t s = blockIdx.x % n_slots;
        while (atomicCAS(&slots[s], 0u, 1u) != 0u) {
            s = s + 1 == n_slots ? 0 : s + 1;
            if (s == 0) __nanosleep(1000);
        }
        slot = s;
    }
    __syncthreads();
    uint32_t r = blockIdx.x * blockDim.x + threadIdx.x;
    if (r < n) {
        // Regions are laid out at the plain (32-byte) stride whatever the list
        // type, so a slot number owns the same bytes for packed and plain
        // launches; packed lists use the first half of their region.
        typedef decltype(L::p) P;
        P s = (P)((char*)scratch +
                  ((uint64_t)slot * SEED_THREADS + threadIdx.x) * (3 * MAXI) * sizeof(Intv));
        collect_read<MAXI, L>(b, table, use_table, seqs, offs, lens, cap, retry_only, out, out_n, r, L{s},
                              L{s + MAXI}, L{s + 2 * MAXI});
    }
    // Every list write is performed before the slot can be reused. A thread
    // only ever reads list entries it wrote itself for the same read, so the
    // next owner never sees this block's data.
    __threadfence();
    __syncthreads();
    if (threadIdx.x == 0) atomicExch(&slots[slot], 0u);
}

static int maxi_index(int maxi) { return maxi == 64 ? 0 : maxi == 96 ? 1 : 2; }

// Resident-block slots for collect_slot_kernel<MAXI>: SMs x blocks per SM.
template<int MAXI>
static cudaError_t seed_slots(uint32_t* n_slots) {
    // Atomic: C ABI callers may initialise from several threads at once.
    static std::atomic<uint32_t> cached{0};
    if (uint32_t c = cached.load(std::memory_order_acquire)) {
        *n_slots = c;
        return cudaSuccess;
    }
    int dev = 0, sms = 0, per_sm = 0, packed_per_sm = 0;
    cudaError_t e = cudaGetDevice(&dev);
    if (e == cudaSuccess) e = cudaDeviceGetAttribute(&sms, cudaDevAttrMultiProcessorCount, dev);
    // The lists stay in L1 when they can: prefer all of it over shared
    // memory (packed: 588 vs 625 ms on a 30x batch). Set before the
    // occupancy queries so they see the launch configuration.
    if (e == cudaSuccess)
        e = cudaFuncSetAttribute(collect_slot_kernel<MAXI, PlainList>, cudaFuncAttributePreferredSharedMemoryCarveout, 0);
    if (e == cudaSuccess)
        e = cudaFuncSetAttribute(collect_slot_kernel<MAXI, PackedList>, cudaFuncAttributePreferredSharedMemoryCarveout,
                                 0);
    if (e == cudaSuccess)
        e = cudaOccupancyMaxActiveBlocksPerMultiprocessor(&per_sm, collect_slot_kernel<MAXI, PlainList>, SEED_THREADS, 0);
    if (e == cudaSuccess)
        e = cudaOccupancyMaxActiveBlocksPerMultiprocessor(&packed_per_sm, collect_slot_kernel<MAXI, PackedList>,
                                                          SEED_THREADS, 0);
    if (e != cudaSuccess) return e;
    // Enough slots for whichever variant is resident.
    if (packed_per_sm > per_sm) per_sm = packed_per_sm;
    uint64_t total = (uint64_t)sms * (uint64_t)per_sm;
    if (total == 0 || total > SEED_MAX_SLOTS) return cudaErrorInvalidConfiguration;
    cached.store((uint32_t)total, std::memory_order_release);
    *n_slots = (uint32_t)total;
    return cudaSuccess;
}

static cudaError_t seed_slots_for(int maxi, uint32_t* n_slots) {
    return maxi == 64 ? seed_slots<64>(n_slots) : maxi == 96 ? seed_slots<96>(n_slots) : seed_slots<128>(n_slots);
}

extern "C" int nomics_seed_scratch_bytes(uint32_t maxi, size_t* bytes) {
    if (maxi != 64 && maxi != 96 && maxi != 128) return cudaErrorInvalidValue;
    uint32_t n_slots = 0;
    cudaError_t e = seed_slots_for((int)maxi, &n_slots);
    if (e != cudaSuccess) return e;
    *bytes = (size_t)n_slots * SEED_THREADS * 3 * maxi * sizeof(Intv);
    return cudaSuccess;
}

template<int MAXI>
static cudaError_t launch_collect(cudaStream_t stream, const nomics_bwt& b, const nomics_kmer_table& table, bool use_table,
                                  const uint8_t* seqs, const uint64_t* offs, const uint32_t* lens, uint32_t n,
                                  uint32_t cap, bool retry_only, Intv* out, uint32_t* out_n,
                                  const nomics_seed_scratch* scratch) {
    int blocks = (int)((n + SEED_THREADS - 1) / SEED_THREADS);
    int at = maxi_index(MAXI);
    void* lists = scratch ? (void*)scratch->lists[at] : nullptr;
    if (lists) {
        uint32_t n_slots = 0;
        cudaError_t e = seed_slots<MAXI>(&n_slots);
        if (e != cudaSuccess) return e;
        if (scratch->bytes[at] < (size_t)n_slots * SEED_THREADS * 3 * MAXI * sizeof(Intv)) return cudaErrorInvalidValue;
        unsigned int* slots = nullptr;
        e = cudaGetSymbolAddress((void**)&slots, g_seed_slots);
        if (e != cudaSuccess) return e;
        // Packed lists hold x < 2^36, i.e. every interval of a BWT shorter
        // than that (a larger one would still be exact, with every read
        // going to the CPU; plain lists avoid that). NOMICS_SEED_PACK=0
        // selects plain lists (measurement only; same output).
        bool packed = b.seq_len + 2 < (1ull << 36);
        if (const char* s = getenv("NOMICS_SEED_PACK")) packed = packed && strcmp(s, "0") != 0;
        slots += at * SEED_MAX_SLOTS;
        if (packed)
            collect_slot_kernel<MAXI, PackedList><<<blocks, SEED_THREADS, 0, stream>>>(
                b, table, use_table, seqs, offs, lens, n, cap, retry_only, out, out_n, lists, slots, n_slots);
        else
            collect_slot_kernel<MAXI, PlainList><<<blocks, SEED_THREADS, 0, stream>>>(
                b, table, use_table, seqs, offs, lens, n, cap, retry_only, out, out_n, lists, slots, n_slots);
    } else {
        collect_kernel<MAXI><<<blocks, SEED_THREADS, 0, stream>>>(b, table, use_table, seqs, offs, lens, n, cap, retry_only, out,
                                                       out_n);
    }
    return cudaGetLastError();
}

extern "C" int nomics_seed_collect_launch(nomics_stream stream_handle, const nomics_bwt* bwt,
                                          const nomics_kmer_table* table_ptr, const uint8_t* seqs,
                                          const uint64_t* offs, const uint32_t* lens, uint32_t n, uint32_t cap,
                                          uint32_t work_mode, nomics_intv* out, uint32_t* out_n,
                                          const nomics_seed_scratch* scratch, nomics_ticket* ticket) {
    if (ticket == nullptr) return cudaErrorInvalidValue;
    *ticket = nullptr;
    if (n == 0) return cudaSuccess;
    if (work_mode != NOMICS_SEED_WORK_64 && work_mode != NOMICS_SEED_WORK_96 && work_mode != NOMICS_SEED_WORK_128 &&
        work_mode != NOMICS_SEED_WORK_TWO_PASS)
        return cudaErrorInvalidValue;
    cudaError_t e = cudaGetLastError();
    if (e != cudaSuccess) {
        fprintf(stderr, "nomics_seed_collect: error pending on entry: %d (%s)\n", (int)e, cudaGetErrorString(e));
        return e;
    }
    cudaStream_t stream = (cudaStream_t)stream_handle;
    nomics_ticket_s* t = nullptr;
    if ((e = ticket_begin(stream, &t)) != cudaSuccess) {
        fprintf(stderr, "nomics_seed_collect: event setup failed: %d (%s)\n", (int)e, cudaGetErrorString(e));
        return e;
    }
    nomics_kmer_table table = {};
    bool use_table = table_ptr != nullptr;
    if (use_table) table = *table_ptr;
    switch (work_mode) {
        case NOMICS_SEED_WORK_64:
            e = launch_collect<64>(stream, *bwt, table, use_table, seqs, offs, lens, n, cap, false, out, out_n, scratch);
            break;
        case NOMICS_SEED_WORK_96:
            e = launch_collect<96>(stream, *bwt, table, use_table, seqs, offs, lens, n, cap, false, out, out_n, scratch);
            break;
        case NOMICS_SEED_WORK_128:
            e = launch_collect<128>(stream, *bwt, table, use_table, seqs, offs, lens, n, cap, false, out, out_n,
                                    scratch);
            break;
        case NOMICS_SEED_WORK_TWO_PASS:
            e = launch_collect<64>(stream, *bwt, table, use_table, seqs, offs, lens, n, cap, false, out, out_n, scratch);
            // The retry pass follows the first on the same stream.
            if (e == cudaSuccess)
                e = launch_collect<128>(stream, *bwt, table, use_table, seqs, offs, lens, n, cap, true, out, out_n,
                                        scratch);
            break;
    }
    if (e != cudaSuccess)
        fprintf(stderr, "nomics_seed_collect: launch failed: %d (%s)\n", (int)e, cudaGetErrorString(e));
    // On error this drains the stream (an earlier pass may be queued) before
    // the caller reuses its buffers.
    return ticket_finish(t, e, ticket);
}
