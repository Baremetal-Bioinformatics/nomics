// Seed search on Metal: bwa 0.7.19's mem_collect_intv (SMEMs, re-seeding,
// LAST-like seeds) for a batch of reads, one thread per read. A port of
// collect_slot_kernel in shim/seed.cu (itself a line-for-line port of
// crates/nomics-align/src/seed.rs); pure integer code. Differences:
//   - The per-thread interval lists live in scratch indexed by the thread's
//     position in its dispatch. The launcher splits a batch into dispatches of
//     at most the scratch's thread count, run serially, so CUDA's slot pool
//     (an atomicCAS spin with __nanosleep, which Metal cannot express safely)
//     is not needed. There is no local-memory variant: scratch is required.
//   - __popc -> popcount; the re-seeding split length is the constant 28 (see
//     collect_read).
#include <metal_stdlib>
using namespace metal;

#define OCC_INTV_SHIFT 7
#define OCC_INTV_MASK 127ul
#define MIN_SEED_LEN 19
#define SPLIT_WIDTH 10
#define MAX_MEM_INTV 20
#define MAXL 256  // NOMICS_SEED_MAX_LEN
#define SEED_OVERFLOW 0xffffffffu  // NOMICS_SEED_OVERFLOW
#define SMEM_OVERFLOW (-2)
#define KMER_MAX_K 13  // NOMICS_KMER_MAX_K

// nomics_intv
struct Intv {
    ulong x[3];
    ulong info;
};
static_assert(sizeof(Intv) == 32, "C ABI layout changed");

// Everything a launch passes by value (the pointer fields of nomics_bwt and
// nomics_kmer_table are bound as buffers instead).
struct SeedParams {
    ulong primary;
    ulong l2[5];
    ulong seq_len;
    ulong level_offsets[KMER_MAX_K + 1];
    uint table_k;
    uint use_table;
    uint n;
    uint cap;
    uint retry_only;
    uint base;  // first read of this dispatch
};

struct Bwt {
    constant SeedParams& p;
    device const uint* words;
};

struct Table {
    constant SeedParams& p;
    device const ulong* intervals;
};

// Occurrences of each base among the 16 2-bit bases of `w`, one count per byte.
static inline uint occ_aux4(uint w) {
    uint hi = w >> 1 & 0x55555555u, lo = w & 0x55555555u;
    uint a = popcount(~hi & ~lo & 0x55555555u);
    uint c = popcount(~hi & lo & 0x55555555u);
    uint g = popcount(hi & ~lo & 0x55555555u);
    uint t = popcount(hi & lo);
    return a | c << 8 | g << 16 | t << 24;
}

// bwt_occ4
static void occ4(Bwt b, ulong k, thread ulong* cnt) {
    if (k == ~0ul) {
        cnt[0] = cnt[1] = cnt[2] = cnt[3] = 0;
        return;
    }
    k -= (k >= b.p.primary);
    device const uint* p = b.words + ((k >> OCC_INTV_SHIFT) << 4);
    device const ulong* c = (device const ulong*)p;
    cnt[0] = c[0];
    cnt[1] = c[1];
    cnt[2] = c[2];
    cnt[3] = c[3];
    p += 8;
    device const uint* end = p + ((k >> 4) - ((k & ~OCC_INTV_MASK) >> 4));
    ulong x = 0;
    for (; p < end; ++p) x += occ_aux4(*p);
    uint tmp = *p & ~((1u << ((~k & 15) << 1)) - 1);
    x += occ_aux4(tmp) - (~k & 15);
    cnt[0] += x & 0xff;
    cnt[1] += x >> 8 & 0xff;
    cnt[2] += x >> 16 & 0xff;
    cnt[3] += x >> 24;
}

// bwt_set_intv
static inline Intv set_intv(Bwt b, int c) {
    Intv r;
    r.x[0] = b.p.l2[c] + 1;
    r.x[1] = b.p.l2[3 - c] + 1;
    r.x[2] = b.p.l2[c + 1] - b.p.l2[c];
    r.info = 0;
    return r;
}

static inline Intv table_intv(Table table, int len, ulong code) {
    device const ulong* p = table.intervals + 3 * (table.p.level_offsets[len] + code);
    Intv r;
    r.x[0] = p[0];
    r.x[1] = p[1];
    r.x[2] = p[2];
    r.info = 0;
    return r;
}

// bwt_extend: is_back selects backward (prepend) or forward extension.
static void extend(Bwt b, thread const Intv& ik, thread Intv* ok, bool is_back) {
    int bi = is_back ? 1 : 0, nb = is_back ? 0 : 1;
    ulong k = ik.x[nb] - 1;
    ulong tk[4], tl[4];
    occ4(b, k, tk);
    occ4(b, k + ik.x[2], tl);
    for (int i = 0; i < 4; ++i) {
        ok[i].x[nb] = b.p.l2[i] + 1 + tk[i];
        ok[i].x[2] = tl[i] - tk[i];
        ok[i].info = 0;
    }
    ulong spans = (ik.x[nb] <= b.p.primary && ik.x[nb] + ik.x[2] - 1 >= b.p.primary) ? 1 : 0;
    ok[3].x[bi] = ik.x[bi] + spans;
    ok[2].x[bi] = ok[3].x[bi] + ok[3].x[2];
    ok[1].x[bi] = ok[2].x[bi] + ok[2].x[2];
    ok[0].x[bi] = ok[1].x[bi] + ok[1].x[2];
}

// Interval lists for smem1, in device scratch. PlainList stores bwtintv_t as
// is. PackedList stores one in 16 bytes: x[0], x[1], x[2] < 2^36 and qbeg,
// qend < 2^10 (reads are at most MAXL bases). put() refuses a value outside
// those bounds, and smem1 then reports SMEM_OVERFLOW, so the read is
// recomputed on the CPU: the packing can never change a result.
struct PlainList {
    device Intv* p;
    Intv get(int i) const { return p[i]; }
    bool put(int i, Intv v) const {
        p[i] = v;
        return true;
    }
};

struct PackedList {
    device uint4* p;
    Intv get(int i) const {
        uint4 w = p[i];
        Intv v;
        v.x[0] = (ulong)w.x | (ulong)(w.w & 0xf) << 32;
        v.x[1] = (ulong)w.y | (ulong)(w.w >> 4 & 0xf) << 32;
        v.x[2] = (ulong)w.z | (ulong)(w.w >> 8 & 0xf) << 32;
        v.info = (ulong)(w.w >> 22) << 32 | (ulong)(w.w >> 12 & 0x3ff);
        return v;
    }
    bool put(int i, Intv v) const {
        ulong qbeg = v.info >> 32, qend = (uint)v.info;
        if (((v.x[0] | v.x[1] | v.x[2]) >> 36) != 0 || qbeg >= 1024 || qend >= 1024) return false;
        uint4 w;
        w.x = (uint)v.x[0];
        w.y = (uint)v.x[1];
        w.z = (uint)v.x[2];
        w.w = (uint)(v.x[0] >> 32) | (uint)(v.x[1] >> 32) << 4 | (uint)(v.x[2] >> 32) << 8 |
              (uint)qend << 12 | (uint)qbeg << 22;
        p[i] = w;
        return true;
    }
};

template <class L>
static void reverse_list(L a, int n) {
    for (int i = 0, j = n - 1; i < j; ++i, --j) {
        Intv t = a.get(i), u = a.get(j);
        a.put(i, u);  // both came out of this list, so both fit
        a.put(j, t);
    }
}

// bwt_smem1 (bwt_smem1a with max_intv = 0). `mem` receives the SMEMs covering
// x sorted by start; returns where the next search starts, or SMEM_OVERFLOW
// if a list outgrew MAXI.
template <int MAXI, class L>
static int smem1(Bwt b, Table table, bool use_table, device const uchar* q, int len, int x,
                 ulong min_intv, L mem, thread int* n_mem, L prev, L curr) {
    *n_mem = 0;
    if (q[x] > 3) return x + 1;
    if (min_intv < 1) min_intv = 1;
    Intv ik = set_intv(b, q[x]), ok[4];
    ik.info = x + 1;
    int n_curr = 0, n_prev;
    int i = x + 1;
    bool stopped = false;
    if (use_table && x + (int)table.p.table_k <= len) {
        bool canonical = true;
        for (int j = 0; j < (int)table.p.table_k; ++j) canonical &= q[x + j] < 4;
        if (canonical) {
            ulong code = q[x];
            for (int prefix_len = 2; prefix_len <= (int)table.p.table_k; ++prefix_len) {
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
            if (!stopped) i = x + table.p.table_k;
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
    reverse_list(curr, n_curr);
    int ret = (int)(uint)curr.get(0).info;
    L t = prev;
    prev = curr;
    curr = t;
    n_prev = n_curr;

    for (i = x - 1; i >= -1; --i) {
        int c = (i < 0 || q[i] > 3) ? -1 : q[i];
        n_curr = 0;
        for (int j = 0; j < n_prev; ++j) {
            const Intv p = prev.get(j);
            bool keep_ext = false;
            if (c >= 0) {
                extend(b, p, ok, true);
                keep_ext = ok[c].x[2] >= min_intv;
            }
            if (keep_ext) {
                if (n_curr == 0 || ok[c].x[2] != curr.get(n_curr - 1).x[2]) {
                    ok[c].info = p.info;
                    if (n_curr == MAXI || !curr.put(n_curr++, ok[c])) return SMEM_OVERFLOW;
                }
            } else if (n_curr == 0 &&
                       (*n_mem == 0 || (ulong)(i + 1) < mem.get(*n_mem - 1).info >> 32)) {
                Intv m = p;
                m.info |= (ulong)(i + 1) << 32;
                if (*n_mem == MAXI || !mem.put((*n_mem)++, m)) return SMEM_OVERFLOW;
            }
        }
        if (n_curr == 0) break;
        t = prev;
        prev = curr;
        curr = t;
        n_prev = n_curr;
    }
    reverse_list(mem, *n_mem);
    return ret;
}

// bwt_seed_strategy1
static int seed_strategy1(Bwt b, Table table, bool use_table, device const uchar* q, int len, int x,
                          int min_len, ulong max_intv, thread Intv* m) {
    m->x[0] = m->x[1] = m->x[2] = m->info = 0;
    if (q[x] > 3) return x + 1;
    Intv ik = set_intv(b, q[x]), ok[4];
    int start = x + 1;
    if (use_table && table.p.table_k <= (uint)min_len && x + (int)table.p.table_k <= len) {
        ulong code = 0;
        bool canonical = true;
        for (int j = 0; j < (int)table.p.table_k; ++j) {
            canonical &= q[x + j] < 4;
            code = code << 2 | q[x + j];
        }
        if (canonical) {
            ik = table_intv(table, table.p.table_k, code);
            start = x + table.p.table_k;
        }
    }
    for (int i = start; i < len; ++i) {
        if (q[i] < 4) {
            int c = 3 - q[i];
            extend(b, ik, ok, false);
            if (ok[c].x[2] < max_intv && i - x >= min_len) {
                *m = ok[c];
                m->info = (ulong)x << 32 | (ulong)(i + 1);
                return i + 1;
            }
            ik = ok[c];
        } else
            return i + 1;
    }
    return len;
}

// mem_collect_intv for read r; mem1, prev and curr are MAXI-entry lists.
template <int MAXI, class L>
static void collect_read(Bwt b, Table table, bool use_table, device const uchar* seqs,
                         device const ulong* offs, device const uint* lens, uint cap, bool retry_only,
                         device Intv* out, device uint* out_n, uint r, L mem1, L prev, L curr) {
    if (retry_only && out_n[r] != SEED_OVERFLOW) return;
    // Compare before narrowing: a length above INT_MAX must overflow too.
    if (lens[r] > (uint)MAXL) {
        out_n[r] = SEED_OVERFLOW;
        return;
    }
    int len = (int)lens[r];
    device const uchar* q = seqs + offs[r];
    device Intv* mem = out + (ulong)r * cap;
    uint n_mem = 0;
    int n1;
    bool overflow = false;
    // Pass 1: SMEMs.
    int x = 0;
    while (x < len && !overflow) {
        if (q[x] < 4) {
            x = smem1<MAXI, L>(b, table, use_table, q, len, x, 1, mem1, &n1, prev, curr);
            if (x == SMEM_OVERFLOW) {
                overflow = true;
                break;
            }
            for (int i = 0; i < n1; ++i) {
                Intv m = mem1.get(i);
                uint sl = (uint)m.info - (uint)(m.info >> 32);
                if (sl >= MIN_SEED_LEN) {
                    if (n_mem == cap) {
                        overflow = true;
                        break;
                    }
                    mem[n_mem++] = m;
                }
            }
        } else
            ++x;
    }
    // Pass 2: re-seed long SMEMs with few occurrences from their midpoint.
    // bwa computes split_len = (int)(min_seed_len * split_factor + .499) in
    // double: (int)(19 * 1.5 + .499) = 28. Metal has no double; the value is
    // a constant of the fixed options (seed.cu notes it is 28 either way).
    const int split_len = 28;
    uint old_n = n_mem;
    for (uint k = 0; k < old_n && !overflow; ++k) {
        Intv p = mem[k];
        int start = (int)(p.info >> 32), end = (int)(uint)p.info;
        if (end - start < split_len || p.x[2] > SPLIT_WIDTH) continue;
        if (smem1<MAXI, L>(b, table, use_table, q, len, (start + end) >> 1, p.x[2] + 1, mem1, &n1,
                           prev, curr) == SMEM_OVERFLOW) {
            overflow = true;
            break;
        }
        for (int i = 0; i < n1; ++i) {
            Intv m = mem1.get(i);
            uint sl = (uint)m.info - (uint)(m.info >> 32);
            if (sl >= MIN_SEED_LEN) {
                if (n_mem == cap) {
                    overflow = true;
                    break;
                }
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
                if (n_mem == cap) {
                    overflow = true;
                    break;
                }
                mem[n_mem++] = m;
            }
        } else
            ++x;
    }
    if (overflow) {
        out_n[r] = SEED_OVERFLOW;
        return;
    }
    // Sort by info. Equal keys are identical intervals, so any sort is bwa's.
    for (uint i = 1; i < n_mem; ++i) {
        Intv t = mem[i];
        uint j = i;
        while (j > 0 && mem[j - 1].info > t.info) {
            mem[j] = mem[j - 1];
            --j;
        }
        mem[j] = t;
    }
    out_n[r] = n_mem;
}

// One thread per read r = base + tid. Its three lists sit at the plain
// (32-byte) stride in `scratch`, whatever the list type, so one scratch serves
// packed and plain launches; packed lists use the first half of their region.
template <int MAXI, class L>
kernel void collect(constant SeedParams& p [[buffer(0)]],
                    device const uint* words [[buffer(1)]],
                    device const ulong* table_intervals [[buffer(2)]],
                    device const uchar* seqs [[buffer(3)]],
                    device const ulong* offs [[buffer(4)]],
                    device const uint* lens [[buffer(5)]],
                    device Intv* out [[buffer(6)]],
                    device uint* out_n [[buffer(7)]],
                    device uchar* scratch [[buffer(8)]],
                    uint tid [[thread_position_in_grid]]) {
    const uint r = p.base + tid;
    if (r >= p.n) return;
    Bwt b = {p, words};
    Table table = {p, table_intervals};
    device uchar* lists = scratch + (ulong)tid * (3 * MAXI) * sizeof(Intv);
    typedef decltype(L::p) P;
    P s = (P)lists;
    collect_read<MAXI, L>(b, table, p.use_table != 0, seqs, offs, lens, p.cap, p.retry_only != 0,
                          out, out_n, r, L{s}, L{s + MAXI}, L{s + 2 * MAXI});
}

#define SEED_KERNEL(MAXI, L, NAME)                                                           \
    template [[host_name(NAME)]] kernel void collect<MAXI, L>(                               \
        constant SeedParams&, device const uint*, device const ulong*, device const uchar*,  \
        device const ulong*, device const uint*, device Intv*, device uint*, device uchar*,  \
        uint);
SEED_KERNEL(64, PlainList, "seed_64_plain")
SEED_KERNEL(96, PlainList, "seed_96_plain")
SEED_KERNEL(128, PlainList, "seed_128_plain")
SEED_KERNEL(64, PackedList, "seed_64_packed")
SEED_KERNEL(96, PackedList, "seed_96_packed")
SEED_KERNEL(128, PackedList, "seed_128_packed")
