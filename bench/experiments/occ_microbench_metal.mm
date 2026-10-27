// Dependent random occ4 lookups over bwa's .bwt on Apple Silicon: Metal
// (reading host memory in place via newBufferWithBytesNoCopy) vs all CPU
// threads. The Metal twin of occ_microbench.cu (same chains, same occ4, same
// seeds), run before porting seed search (Metal sprint A4, docs/strategy/metal-sprint.md).
//
//   xcrun clang++ -std=c++17 -O3 -fobjc-arc -framework Metal -framework Foundation \
//       occ_microbench_metal.mm -o occ_microbench_metal
//   ./occ_microbench_metal <prefix>.64.bwt [chains_gpu] [steps]
//
// Each "chain" mimics an LF walk: occ4 at row k (one 64-byte block read),
// then k' = L2[c] + occ[c] for a base c derived from the counts, so every
// lookup depends on the previous one, as in backward search. The GPU's first
// sample of chains is checked against the CPU's walk of the same chains.
#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <sys/sysctl.h>
#include <time.h>

typedef struct {
    uint64_t primary, l2[5], seq_len;
} BwtHdr;

static inline uint64_t occ_aux4(const uint32_t* cnt, uint32_t b) {
    return cnt[b & 0xff] + cnt[b >> 8 & 0xff] + cnt[b >> 16 & 0xff] + cnt[b >> 24];
}

static inline void occ4(const BwtHdr* b, const uint32_t* w, const uint32_t* cnt, uint64_t k,
                        uint64_t out[4]) {
    k -= (k >= b->primary);
    const uint32_t* p = w + ((k >> 7) << 4);
    const uint64_t* c = (const uint64_t*)p;
    out[0] = c[0]; out[1] = c[1]; out[2] = c[2]; out[3] = c[3];
    p += 8;
    const uint32_t* end = p + ((k >> 4) - ((k & ~127ull) >> 4));
    uint64_t x = 0;
    for (; p < end; ++p) x += occ_aux4(cnt, *p);
    uint32_t tmp = *p & ~((1u << ((~k & 15) << 1)) - 1);
    x += occ_aux4(cnt, tmp) - (~k & 15);
    out[0] += x & 0xff; out[1] += x >> 8 & 0xff; out[2] += x >> 16 & 0xff; out[3] += x >> 24;
}

static inline uint64_t walk(const BwtHdr* b, const uint32_t* w, const uint32_t* cnt, uint64_t k,
                            int steps) {
    uint64_t o[4], acc = 0;
    for (int i = 0; i < steps; ++i) {
        occ4(b, w, cnt, k, o);
        int c = (int)((o[0] ^ o[1] ^ o[2] ^ o[3] ^ k) & 3);
        k = (b->l2[c] + o[c] + (k * 0x9e3779b97f4a7c15ull >> 40)) % b->seq_len;
        acc += k;
    }
    return acc;
}

static const char* const kSource = R"MSL(
#include <metal_stdlib>
using namespace metal;
struct BwtHdr { ulong primary; ulong l2[5]; ulong seq_len; };

static inline ulong occ_aux4(constant uint* cnt, uint b) {
    return cnt[b & 0xff] + cnt[b >> 8 & 0xff] + cnt[b >> 16 & 0xff] + cnt[b >> 24];
}

static inline void occ4(constant BwtHdr& b, device const uint* w, constant uint* cnt, ulong k,
                        thread ulong* out) {
    k -= (k >= b.primary);
    device const uint* p = w + ((k >> 7) << 4);
    device const ulong* c = (device const ulong*)p;
    out[0] = c[0]; out[1] = c[1]; out[2] = c[2]; out[3] = c[3];
    p += 8;
    device const uint* end = p + ((k >> 4) - ((k & ~127ul) >> 4));
    ulong x = 0;
    for (; p < end; ++p) x += occ_aux4(cnt, *p);
    uint tmp = *p & ~((1u << ((~k & 15) << 1)) - 1);
    x += occ_aux4(cnt, tmp) - (~k & 15);
    out[0] += x & 0xff; out[1] += x >> 8 & 0xff; out[2] += x >> 16 & 0xff; out[3] += x >> 24;
}

kernel void chains(constant BwtHdr& b [[buffer(0)]], device const uint* w [[buffer(1)]],
                   constant uint* cnt [[buffer(2)]], constant ulong& n [[buffer(3)]],
                   constant int& steps [[buffer(4)]], device ulong* sink [[buffer(5)]],
                   uint t [[thread_position_in_grid]]) {
    if (t >= n) return;
    ulong k = ((ulong)t * 0x9e3779b97f4a7c15ul) % b.seq_len;
    ulong o[4], acc = 0;
    for (int i = 0; i < steps; ++i) {
        occ4(b, w, cnt, k, o);
        int c = (int)((o[0] ^ o[1] ^ o[2] ^ o[3] ^ k) & 3);
        k = (b.l2[c] + o[c] + (k * 0x9e3779b97f4a7c15ul >> 40)) % b.seq_len;
        acc += k;
    }
    sink[t] = acc;
}
)MSL";

typedef struct {
    const BwtHdr* b;
    const uint32_t* w;
    const uint32_t* cnt;
    uint64_t first, stride, end;
    int steps;
    uint64_t acc;
} CpuJob;

static void* cpu_thread(void* arg) {
    CpuJob* j = (CpuJob*)arg;
    uint64_t acc = 0;
    for (uint64_t c = j->first; c < j->end; c += j->stride)
        acc += walk(j->b, j->w, j->cnt, (c * 0x9e3779b97f4a7c15ull) % j->b->seq_len, j->steps);
    j->acc = acc;
    return NULL;
}

static double now(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

int main(int argc, char** argv) {
    @autoreleasepool {
        if (argc < 2) {
            fprintf(stderr, "usage: %s <prefix>.64.bwt [chains_gpu] [steps]\n", argv[0]);
            return 2;
        }
        const char* path = argv[1];
        uint64_t chains_gpu = argc > 2 ? strtoull(argv[2], 0, 10) : (1u << 20);
        int steps = argc > 3 ? atoi(argv[3]) : 200;
        FILE* f = fopen(path, "rb");
        if (!f) { perror(path); return 1; }
        fseek(f, 0, SEEK_END);
        size_t size = ftell(f);
        fseek(f, 0, SEEK_SET);
        BwtHdr b;
        fread(&b.primary, 8, 1, f);
        fread(b.l2 + 1, 8, 4, f);
        b.l2[0] = 0;
        b.seq_len = b.l2[4];
        size_t words = size - 40;
        // 2 MiB aligned and rounded, as HostBuf allocates (whole 16 KiB pages).
        size_t map = (words + (2 << 20)) & ~((size_t)(2 << 20) - 1);
        void* mem = mmap(0, map, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
        if (mem == MAP_FAILED) { perror("mmap"); return 1; }
        double t_read = now();
        if (fread(mem, 1, words, f) != words) { fprintf(stderr, "short read\n"); return 1; }
        fclose(f);
        printf("bwt: %.2f GB, seq_len %llu, read in %.1f s\n", words / 1e9,
               (unsigned long long)b.seq_len, now() - t_read);
        const uint32_t* w = (const uint32_t*)mem;
        uint32_t cnt[256];
        for (int i = 0; i < 256; ++i) {
            uint32_t x = 0;
            for (int j = 0; j < 4; ++j)
                x |= (((i & 3) == j) + ((i >> 2 & 3) == j) + ((i >> 4 & 3) == j) + (i >> 6 == j))
                     << (j << 3);
            cnt[i] = x;
        }

        // CPU: every hardware thread, each running chains sequentially.
        int nt = 0;
        size_t len = sizeof nt;
        sysctlbyname("hw.ncpu", &nt, &len, NULL, 0);
        uint64_t chains_cpu = (uint64_t)nt * 4096;
        CpuJob jobs[256];
        pthread_t th[256];
        double t0 = now();
        for (int t = 0; t < nt; ++t) {
            jobs[t] = CpuJob{&b, w, cnt, (uint64_t)t, (uint64_t)nt, chains_cpu, steps, 0};
            pthread_create(&th[t], NULL, cpu_thread, &jobs[t]);
        }
        for (int t = 0; t < nt; ++t) pthread_join(th[t], NULL);
        double cpu_s = now() - t0;
        printf("CPU %d threads: %.1f M lookups/s (%.2f s)\n", nt,
               chains_cpu * (double)steps / cpu_s / 1e6, cpu_s);

        // GPU: one thread per chain, BWT read in place from host memory.
        id<MTLDevice> dev = MTLCreateSystemDefaultDevice();
        NSError* err = nil;
        MTLCompileOptions* opts = [MTLCompileOptions new];
        opts.mathMode = MTLMathModeSafe;
        id<MTLLibrary> lib = [dev newLibraryWithSource:[NSString stringWithUTF8String:kSource] options:opts error:&err];
        if (!lib) { NSLog(@"compile: %@", err); return 1; }
        id<MTLComputePipelineState> pso =
            [dev newComputePipelineStateWithFunction:[lib newFunctionWithName:@"chains"] error:&err];
        if (!pso) { NSLog(@"pipeline: %@", err); return 1; }
        double t_wrap = now();
        id<MTLBuffer> bw = [dev newBufferWithBytesNoCopy:mem
                                                  length:map
                                                 options:MTLResourceStorageModeShared
                                             deallocator:nil];
        if (!bw) { fprintf(stderr, "no-copy wrap of %.2f GB failed\n", map / 1e9); return 1; }
        printf("GPU: %s, wrapped %.2f GB in %.0f ms\n", dev.name.UTF8String, map / 1e9,
               (now() - t_wrap) * 1e3);
        id<MTLBuffer> sink = [dev newBufferWithLength:chains_gpu * 8
                                              options:MTLResourceStorageModeShared];
        id<MTLCommandQueue> q = [dev newCommandQueue];
        for (int rep = 0; rep < 3; ++rep) {  // the first run warms the GPU's page mappings
            id<MTLCommandBuffer> cmd = [q commandBuffer];
            id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
            [enc setComputePipelineState:pso];
            [enc setBytes:&b length:sizeof b atIndex:0];
            [enc setBuffer:bw offset:0 atIndex:1];
            [enc setBytes:cnt length:sizeof cnt atIndex:2];
            [enc setBytes:&chains_gpu length:8 atIndex:3];
            [enc setBytes:&steps length:4 atIndex:4];
            [enc setBuffer:sink offset:0 atIndex:5];
            [enc dispatchThreads:MTLSizeMake(chains_gpu, 1, 1)
                threadsPerThreadgroup:MTLSizeMake(256, 1, 1)];
            [enc endEncoding];
            [cmd commit];
            [cmd waitUntilCompleted];
            double s = cmd.GPUEndTime - cmd.GPUStartTime;
            printf("GPU run %d: %.1f M lookups/s (%.2f s)%s\n", rep,
                   chains_gpu * (double)steps / s / 1e6, s,
                   cmd.status == MTLCommandBufferStatusCompleted ? "" : " FAILED");
        }
        // Spot-check: the GPU's chains equal the CPU's walk of the same chains.
        const uint64_t* got = (const uint64_t*)sink.contents;
        int bad = 0;
        const uint64_t sample = chains_gpu < 4096 ? chains_gpu : 4096;
        for (uint64_t t = 0; t < sample; ++t)
            bad += got[t] != walk(&b, w, cnt, (t * 0x9e3779b97f4a7c15ull) % b.seq_len, steps);
        printf("check: %d of %llu sampled chains differ from the CPU\n", bad,
               (unsigned long long)sample);
        return bad != 0;
    }
}
