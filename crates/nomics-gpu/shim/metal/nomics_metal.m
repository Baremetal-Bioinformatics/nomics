// Metal backend for the C ABI in ../nomics_gpu.h (built with `--features metal`).
// Kernels live in shim/metal/*.metal; build.rs compiles them into one
// metallib and embeds its bytes as nomics_metallib / nomics_metallib_len.
//
// Host memory is read in place: a launch pointer inside a registered range
// (nomics_host_register) binds that range's wrapper at its offset; any other
// pointer is wrapped with newBufferWithBytesNoCopy over the whole pages that
// contain it, for that launch only (Apple Silicon shares one physical memory
// between CPU and GPU). HostBuf allocations are 2 MiB aligned, so they cover
// whole 16 KiB pages and a wrap never reaches outside the mapping.
//
// Streams are MTLCommandQueues (no priorities on Metal; the flag is ignored)
// and a ticket is one committed MTLCommandBuffer. Launches encode and commit
// at the very end, so a failed launch has enqueued nothing.
#import <Foundation/Foundation.h>
#import <IOKit/IOKitLib.h>
#import <Metal/Metal.h>
#include <os/lock.h>
#include <pthread.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "../nomics_gpu.h"

extern const unsigned char nomics_metallib[];
extern const size_t nomics_metallib_len;

// Status codes returned in place of a cudaError_t (0 = success).
enum {
    NM_OK = 0,
    NM_UNSUPPORTED = 1000,  // kernel not ported to Metal: callers use the CPU
    NM_NO_DEVICE = 1001,
    NM_LIBRARY = 1002,
    NM_PIPELINE = 1003,
    NM_WRAP = 1004,
    NM_ALLOC = 1005,
    NM_COMMAND = 1006,
    NM_INVALID = 1007,
};

typedef struct {
    id<MTLDevice> device;
    id<MTLLibrary> library;
    id<MTLCommandQueue> queue;  // synchronous smoke calls
    int status;
} nm_context;

static nm_context* nm_get(void) {
    static nm_context ctx;
    static dispatch_once_t once;
    dispatch_once(&once, ^{
        ctx.device = MTLCreateSystemDefaultDevice();
        if (!ctx.device) {
            ctx.status = NM_NO_DEVICE;
            return;
        }
        dispatch_data_t data = dispatch_data_create(nomics_metallib, nomics_metallib_len, NULL,
                                                    DISPATCH_DATA_DESTRUCTOR_DEFAULT);
        NSError* error = nil;
        ctx.library = [ctx.device newLibraryWithData:data error:&error];
        if (!ctx.library) {
            NSLog(@"nomics-gpu: metallib load failed: %@", error);
            ctx.status = NM_LIBRARY;
            return;
        }
        ctx.queue = [ctx.device newCommandQueue];
        ctx.status = ctx.queue ? NM_OK : NM_NO_DEVICE;
    });
    return &ctx;
}

static id<MTLComputePipelineState> nm_pipeline(nm_context* ctx, NSString* name, int* status) {
    id<MTLFunction> fn = [ctx->library newFunctionWithName:name];
    if (!fn) {
        *status = NM_PIPELINE;
        return nil;
    }
    NSError* error = nil;
    id<MTLComputePipelineState> pso = [ctx->device newComputePipelineStateWithFunction:fn
                                                                                 error:&error];
    if (!pso) {
        NSLog(@"nomics-gpu: pipeline %@ failed: %@", name, error);
        *status = NM_PIPELINE;
    }
    return pso;
}

// ---- Host buffer registry ----

@interface NMRegion : NSObject
@property(nonatomic) uintptr_t start;
@property(nonatomic) uintptr_t end;
@property(nonatomic, strong) id<MTLBuffer> buffer;
@end
@implementation NMRegion
@end

static os_unfair_lock nm_regions_lock = OS_UNFAIR_LOCK_INIT;
static NSMutableArray<NMRegion*>* nm_regions;

int nomics_host_register(const void* ptr, size_t len) {
    @autoreleasepool {
        nm_context* ctx = nm_get();
        if (ctx->status) return ctx->status;
        const uintptr_t page = (uintptr_t)getpagesize();
        const uintptr_t start = (uintptr_t)ptr;
        if (!ptr || !len || (start & (page - 1)) || (len & (page - 1))) return NM_INVALID;
        if (len > ctx->device.maxBufferLength) return NM_WRAP;
        id<MTLBuffer> buffer = [ctx->device newBufferWithBytesNoCopy:(void*)ptr
                                                              length:len
                                                             options:MTLResourceStorageModeShared
                                                         deallocator:nil];
        if (!buffer) return NM_WRAP;
        NMRegion* region = [NMRegion new];
        region.start = start;
        region.end = start + len;
        region.buffer = buffer;
        os_unfair_lock_lock(&nm_regions_lock);
        if (!nm_regions) nm_regions = [NSMutableArray new];
        [nm_regions addObject:region];
        os_unfair_lock_unlock(&nm_regions_lock);
        return NM_OK;
    }
}

int nomics_host_unregister(const void* ptr) {
    @autoreleasepool {
        int status = NM_INVALID;
        os_unfair_lock_lock(&nm_regions_lock);
        for (NSUInteger i = 0; i < nm_regions.count; i++) {
            if (nm_regions[i].start == (uintptr_t)ptr) {
                [nm_regions removeObjectAtIndex:i];
                status = NM_OK;
                break;
            }
        }
        os_unfair_lock_unlock(&nm_regions_lock);
        return status;
    }
}

// The buffer to bind for [p, p + n): a registered range containing it, else a
// per-launch wrap of its pages. *offset is p's offset into the buffer. The
// command buffer retains whatever is bound, so a per-launch wrap lives until
// the launch completes.
static id<MTLBuffer> nm_bind(nm_context* ctx, const void* p, size_t n, NSUInteger* offset) {
    const uintptr_t at = (uintptr_t)p;
    id<MTLBuffer> found = nil;
    os_unfair_lock_lock(&nm_regions_lock);
    for (NMRegion* region in nm_regions) {
        if (at >= region.start && at + n <= region.end) {
            found = region.buffer;
            *offset = (NSUInteger)(at - region.start);
            break;
        }
    }
    os_unfair_lock_unlock(&nm_regions_lock);
    if (found) return found;
    const uintptr_t page = (uintptr_t)getpagesize();
    const uintptr_t start = at & ~(page - 1);
    const uintptr_t end = (at + (n ? n : 1) + page - 1) & ~(page - 1);
    if (end - start > ctx->device.maxBufferLength) return nil;
    *offset = (NSUInteger)(at - start);
    return [ctx->device newBufferWithBytesNoCopy:(void*)start
                                          length:end - start
                                         options:MTLResourceStorageModeShared
                                     deallocator:nil];
}

// ---- Device query and smoke test ----

// GPU core count from the IORegistry (Metal has no API for it); 0 if unknown.
static int nm_gpu_cores(void) {
    int cores = 0;
    io_iterator_t it = 0;
    if (IOServiceGetMatchingServices(kIOMainPortDefault, IOServiceMatching("AGXAccelerator"),
                                     &it) != KERN_SUCCESS) {
        return 0;
    }
    io_object_t entry;
    while (!cores && (entry = IOIteratorNext(it))) {
        CFTypeRef value = IORegistryEntryCreateCFProperty(entry, CFSTR("gpu-core-count"),
                                                          kCFAllocatorDefault, 0);
        if (value && CFGetTypeID(value) == CFNumberGetTypeID()) {
            CFNumberGetValue((CFNumberRef)value, kCFNumberIntType, &cores);
        }
        if (value) CFRelease(value);
        IOObjectRelease(entry);
    }
    IOObjectRelease(it);
    return cores;
}

int nomics_device_info_get(int device, nomics_device_info* out) {
    @autoreleasepool {
        if (device != 0) return NM_INVALID;
        nm_context* ctx = nm_get();
        if (ctx->status) return ctx->status;
        memset(out, 0, sizeof(*out));
        strncpy(out->name, ctx->device.name.UTF8String, sizeof(out->name) - 1);
        // Highest MTLGPUFamilyAppleN the device supports, reported as (N, 0).
        for (int family = 12; family >= 1; family--) {
            if ([ctx->device supportsFamily:(MTLGPUFamily)(MTLGPUFamilyApple1 - 1 + family)]) {
                out->major = family;
                break;
            }
        }
        out->multiprocessors = nm_gpu_cores();
        out->pageable_memory_access = 1;  // via newBufferWithBytesNoCopy
        out->total_global_mem = (size_t)ctx->device.recommendedMaxWorkingSetSize;
        out->max_buffer_bytes = (size_t)ctx->device.maxBufferLength;
        return NM_OK;
    }
}

int nomics_sum_u64_in_place(const uint64_t* data, size_t n, uint64_t* out_sum, float* out_ms) {
    @autoreleasepool {
        nm_context* ctx = nm_get();
        if (ctx->status) return ctx->status;
        static id<MTLComputePipelineState> pso;
        static dispatch_once_t once;
        static int pso_status;
        dispatch_once(&once, ^{
            pso = nm_pipeline(ctx, @"sum_u64", &pso_status);
        });
        if (pso_status) return pso_status;

        NSUInteger offset = 0;
        id<MTLBuffer> input = nm_bind(ctx, data, n * sizeof(uint64_t), &offset);
        if (!input) return NM_WRAP;
        const NSUInteger width = MIN((NSUInteger)256, pso.maxTotalThreadsPerThreadgroup);
        const NSUInteger groups = 4096;
        id<MTLBuffer> partial = [ctx->device newBufferWithLength:groups * sizeof(uint64_t)
                                                         options:MTLResourceStorageModeShared];
        if (!partial) return NM_ALLOC;
        const uint64_t count = n;

        id<MTLCommandBuffer> cmd = [ctx->queue commandBuffer];
        id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
        if (!enc) return NM_COMMAND;  // a nil encoder would commit an empty buffer
        [enc setComputePipelineState:pso];
        [enc setBuffer:input offset:offset atIndex:0];
        [enc setBuffer:partial offset:0 atIndex:1];
        [enc setBytes:&count length:sizeof(count) atIndex:2];
        [enc setThreadgroupMemoryLength:width * sizeof(uint64_t) atIndex:0];
        [enc dispatchThreadgroups:MTLSizeMake(groups, 1, 1)
            threadsPerThreadgroup:MTLSizeMake(width, 1, 1)];
        [enc endEncoding];
        [cmd commit];
        [cmd waitUntilCompleted];
        if (cmd.status != MTLCommandBufferStatusCompleted) {
            NSLog(@"nomics-gpu: sum_u64 failed: %@", cmd.error);
            return NM_COMMAND;
        }
        const uint64_t* parts = (const uint64_t*)partial.contents;
        uint64_t sum = 0;
        for (NSUInteger i = 0; i < groups; i++) sum += parts[i];
        *out_sum = sum;
        if (out_ms) *out_ms = (float)((cmd.GPUEndTime - cmd.GPUStartTime) * 1e3);
        return NM_OK;
    }
}

// ---- Streams and tickets ----

// Wait for every command buffer committed to `queue` so far. A queue starts
// its command buffers in order, so an empty one completing after them all
// means they have finished.
static int nm_drain(id<MTLCommandQueue> queue) {
    id<MTLCommandBuffer> drain = [queue commandBuffer];
    if (!drain) return NM_COMMAND;
    [drain commit];
    [drain waitUntilCompleted];
    return drain.status == MTLCommandBufferStatusCompleted ? NM_OK : NM_COMMAND;
}

int nomics_stream_create(int priority, nomics_stream* out) {
    (void)priority;  // Metal queues have no priorities
    if (!out) return NM_INVALID;
    *out = NULL;
    @autoreleasepool {
        nm_context* ctx = nm_get();
        if (ctx->status) return ctx->status;
        id<MTLCommandQueue> queue = [ctx->device newCommandQueue];
        if (!queue) return NM_NO_DEVICE;
        *out = (__bridge_retained void*)queue;
        return NM_OK;
    }
}

int nomics_stream_warm(nomics_stream stream) {
    // Metal command queues need no up-front context allocation.
    return stream ? NM_OK : NM_INVALID;
}

int nomics_stream_destroy(nomics_stream stream) {
    if (!stream) return NM_OK;
    @autoreleasepool {
        id<MTLCommandQueue> queue = (__bridge_transfer id<MTLCommandQueue>)stream;
        return nm_drain(queue);
    }
}

int nomics_ticket_wait(nomics_ticket ticket, float* out_ms) {
    if (out_ms) *out_ms = 0.0f;
    if (!ticket) return NM_OK;  // nothing was launched
    @autoreleasepool {
        id<MTLCommandBuffer> cmd = (__bridge_transfer id<MTLCommandBuffer>)(void*)ticket;
        [cmd waitUntilCompleted];
        if (cmd.status != MTLCommandBufferStatusCompleted) {
            NSLog(@"nomics-gpu: launch failed: %@", cmd.error);
            // As the ABI requires: drain the launch's stream before the caller
            // touches buffers that later launches on it may still use.
            nm_drain(cmd.commandQueue);
            return NM_COMMAND;
        }
        if (out_ms) *out_ms = (float)((cmd.GPUEndTime - cmd.GPUStartTime) * 1e3);
        return NM_OK;
    }
}

static nomics_ticket nm_commit(id<MTLCommandBuffer> cmd) {
    [cmd commit];
    return (nomics_ticket)(__bridge_retained void*)cmd;
}

// ---- BWA ksw_align2: wavefront mate rescue (ksw.metal) ----

int nomics_ksw_align2_wave_launch(nomics_stream stream, const uint8_t* queries,
                                  const uint8_t* targets, const int8_t* matrices,
                                  const nomics_ksw_job* jobs, uint32_t n, uint16_t* row_maxima,
                                  nomics_ksw_result* out, nomics_ticket* ticket) {
    (void)row_maxima;  // rows live in threadgroup memory, as in ksw.cu
    if (!ticket || !stream) return NM_INVALID;
    *ticket = NULL;
    if (n == 0) return NM_OK;
    @autoreleasepool {
        nm_context* ctx = nm_get();
        if (ctx->status) return ctx->status;
        static const int widths[5] = {4, 5, 6, 8, 16};
        static id<MTLComputePipelineState> psos[5];
        static dispatch_once_t once;
        static int pso_status;
        dispatch_once(&once, ^{
            for (int i = 0; i < 5 && !pso_status; i++) {
                psos[i] = nm_pipeline(ctx, [NSString stringWithFormat:@"ksw_wave_%d", widths[i]],
                                      &pso_status);
                // The kernel is written for 32-wide simdgroups.
                if (!pso_status && psos[i].threadExecutionWidth != 32) pso_status = NM_PIPELINE;
                if (!pso_status && getenv("NOMICS_METAL_INFO")) {
                    fprintf(stderr, "[metal] ksw_wave_%d: %lu threads/threadgroup max, %lu B static\n",
                            widths[i], (unsigned long)psos[i].maxTotalThreadsPerThreadgroup,
                            (unsigned long)psos[i].staticThreadgroupMemoryLength);
                }
            }
        });
        if (pso_status) return pso_status;

        uint32_t widest = 0;   // padded query columns
        uint32_t longest = 0;  // target bases
        size_t query_end = 0, target_end = 0, mat_end = 0;
        for (uint32_t i = 0; i < n; ++i) {
            const nomics_ksw_job* job = &jobs[i];
            if (job->m > 8) return NM_INVALID;  // score table holds 8 bases
            // As in ksw.cu: scoring the kernel does not reproduce goes to the
            // CPU. BWA's byte path casts e and o + e to uint8_t; its 16-bit
            // path saturates at 32767, and H <= qlen * max score.
            if (job->size == 1 && (job->e_del > 255 || job->e_ins > 255 || job->o_del + job->e_del > 255 ||
                                   job->o_ins + job->e_ins > 255))
                return NM_INVALID;
            if (job->size != 1) {
                int32_t top = 0;
                for (int32_t k = 0; k < job->m * job->m; ++k)
                    top = MAX(top, (int32_t)matrices[job->mat_offset + k]);
                if ((int64_t)job->query_len * top > 32767) return NM_INVALID;
            }
            longest = longest > job->target_len ? longest : job->target_len;
            const uint32_t lanes = job->size == 1 ? 16 : 8;
            const uint32_t padded = (job->query_len + lanes - 1) / lanes * lanes;
            widest = widest > padded ? widest : padded;
            query_end = MAX(query_end, (size_t)(job->query_offset + job->query_len));
            target_end = MAX(target_end, (size_t)(job->target_offset + job->target_len));
            mat_end = MAX(mat_end, (size_t)(job->mat_offset + (size_t)job->m * (size_t)job->m));
        }
        if (widest > 32 * 16) return NM_INVALID;
        // In size_t, and bounded before narrowing: a huge target must not wrap
        // to a small threadgroup allocation that staging then overruns.
        const size_t per_warp_bytes = 3 * (((size_t)longest + 15) & ~(size_t)15);
        if (per_warp_bytes > ctx->device.maxThreadgroupMemoryLength) return NM_INVALID;
        const uint32_t per_warp = (uint32_t)per_warp_bytes;
        // Up to four simdgroups per threadgroup, as many as fit in its memory
        // (32 KiB on Apple GPUs, against 48 KiB of CUDA shared memory).
        // Wide variants hold more per-thread state, which can lower the
        // pipeline's thread limit; dispatching past it fails the command
        // buffer (and the whole batch would fall back to the CPU).
        const int c = widest <= 32 * 4 ? 0 : widest <= 32 * 5 ? 1 : widest <= 32 * 6 ? 2
                    : widest <= 32 * 8 ? 3 : 4;
        id<MTLComputePipelineState> pso = psos[c];
        const NSUInteger limit = ctx->device.maxThreadgroupMemoryLength;
        NSUInteger warps = MIN((NSUInteger)4, per_warp ? limit / per_warp : 4);
        warps = MIN(warps, pso.maxTotalThreadsPerThreadgroup / 32);
        if (warps == 0) return NM_INVALID;

        NSUInteger off_q = 0, off_t = 0, off_m = 0, off_j = 0, off_o = 0;
        id<MTLBuffer> bq = nm_bind(ctx, queries, query_end, &off_q);
        id<MTLBuffer> bt = nm_bind(ctx, targets, target_end, &off_t);
        id<MTLBuffer> bm = nm_bind(ctx, matrices, mat_end, &off_m);
        id<MTLBuffer> bj = nm_bind(ctx, jobs, (size_t)n * sizeof(nomics_ksw_job), &off_j);
        id<MTLBuffer> bo = nm_bind(ctx, out, (size_t)n * sizeof(nomics_ksw_result), &off_o);
        if (!bq || !bt || !bm || !bj || !bo) return NM_WRAP;

        id<MTLCommandQueue> queue = (__bridge id<MTLCommandQueue>)stream;
        id<MTLCommandBuffer> cmd = [queue commandBuffer];
        id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
        // Nothing is committed on these paths: a nil encoder would commit an
        // empty buffer whose wait "succeeds" with stale results.
        if (!cmd || !enc) return NM_COMMAND;
        [enc setComputePipelineState:pso];
        [enc setBuffer:bq offset:off_q atIndex:0];
        [enc setBuffer:bt offset:off_t atIndex:1];
        [enc setBuffer:bm offset:off_m atIndex:2];
        [enc setBuffer:bj offset:off_j atIndex:3];
        [enc setBuffer:bo offset:off_o atIndex:4];
        [enc setBytes:&n length:sizeof(n) atIndex:5];
        [enc setBytes:&per_warp length:sizeof(per_warp) atIndex:6];
        [enc setThreadgroupMemoryLength:MAX((NSUInteger)16, warps * per_warp) atIndex:0];
        [enc dispatchThreadgroups:MTLSizeMake((n + warps - 1) / warps, 1, 1)
            threadsPerThreadgroup:MTLSizeMake(32 * warps, 1, 1)];
        [enc endEncoding];
        *ticket = nm_commit(cmd);
        return NM_OK;
    }
}

// ---- Not ported to Metal (yet): callers fall back to the CPU ----

int nomics_ksw_align2_launch(nomics_stream stream, const uint8_t* queries, const uint8_t* targets,
                             const int8_t* matrices, const nomics_ksw_job* jobs, uint32_t n,
                             uint16_t* scratch, uint16_t* row_maxima, nomics_ksw_result* out,
                             nomics_ticket* ticket) {
    (void)stream, (void)queries, (void)targets, (void)matrices, (void)jobs, (void)n;
    (void)scratch, (void)row_maxima, (void)out;
    if (ticket) *ticket = NULL;
    return NM_UNSUPPORTED;
}

// ---- Seed search: mem_collect_intv (seed.metal) ----
//
// The kernel's per-thread interval lists live in caller scratch indexed by
// the thread's position in its dispatch, so a batch runs as serial dispatches
// of at most SEED_DISPATCH reads and scratch is sized by that, not by the
// batch (CUDA sizes it by resident blocks and claims slots with a spin lock,
// which Metal cannot express safely). Seed launches share the scratch, so
// they are serialized: each waits for the previous seed launch to complete
// before it commits (on the CPU, so a failed launch cannot deadlock the next).

#define SEED_DISPATCH 131072u

// seed.metal's SeedParams.
typedef struct {
    uint64_t primary;
    uint64_t l2[5];
    uint64_t seq_len;
    uint64_t level_offsets[NOMICS_KMER_MAX_K + 1];
    uint32_t table_k;
    uint32_t use_table;
    uint32_t n;
    uint32_t cap;
    uint32_t retry_only;
    uint32_t base;
} nm_seed_params;
_Static_assert(sizeof(nm_seed_params) == 192, "seed.metal SeedParams layout");
_Static_assert(sizeof(nomics_intv) == 32, "seed.metal Intv layout");

static pthread_mutex_t nm_seed_lock = PTHREAD_MUTEX_INITIALIZER;
static id<MTLCommandBuffer> nm_seed_last;  // guarded by nm_seed_lock

int nomics_seed_scratch_bytes(uint32_t list_size, size_t* bytes) {
    if (!bytes) return NM_INVALID;
    *bytes = 0;
    if (list_size != 64 && list_size != 96 && list_size != 128) return NM_INVALID;
    *bytes = (size_t)SEED_DISPATCH * 3 * list_size * sizeof(nomics_intv);
    return NM_OK;
}

int nomics_seed_collect_launch(nomics_stream stream, const nomics_bwt* bwt,
                               const nomics_kmer_table* table, const uint8_t* seqs,
                               const uint64_t* offs, const uint32_t* lens, uint32_t n,
                               uint32_t cap, uint32_t work, nomics_intv* out, uint32_t* out_n,
                               const nomics_seed_scratch* scratch, nomics_ticket* ticket) {
    if (!ticket || !stream || !bwt) return NM_INVALID;
    *ticket = NULL;
    if (n == 0) return NM_OK;
    if (work != NOMICS_SEED_WORK_64 && work != NOMICS_SEED_WORK_96 &&
        work != NOMICS_SEED_WORK_128 && work != NOMICS_SEED_WORK_TWO_PASS)
        return NM_INVALID;
    @autoreleasepool {
        nm_context* ctx = nm_get();
        if (ctx->status) return ctx->status;
        static NSString* const names[2][3] = {
            {@"seed_64_plain", @"seed_96_plain", @"seed_128_plain"},
            {@"seed_64_packed", @"seed_96_packed", @"seed_128_packed"}};
        static id<MTLComputePipelineState> psos[2][3];
        static dispatch_once_t once;
        static int pso_status;
        dispatch_once(&once, ^{
            for (int p = 0; p < 2 && !pso_status; p++)
                for (int i = 0; i < 3 && !pso_status; i++)
                    psos[p][i] = nm_pipeline(ctx, names[p][i], &pso_status);
        });
        if (pso_status) return pso_status;

        // The passes this mode runs: list size index (0: 64, 1: 96, 2: 128)
        // and whether it only retries the previous pass's overflow reads.
        int passes[2][2], n_passes = 0;
        switch (work) {
            case NOMICS_SEED_WORK_64: passes[n_passes][0] = 0, passes[n_passes++][1] = 0; break;
            case NOMICS_SEED_WORK_96: passes[n_passes][0] = 1, passes[n_passes++][1] = 0; break;
            case NOMICS_SEED_WORK_128: passes[n_passes][0] = 2, passes[n_passes++][1] = 0; break;
            default:
                passes[n_passes][0] = 0, passes[n_passes++][1] = 0;
                passes[n_passes][0] = 2, passes[n_passes++][1] = 1;
        }
        static const uint32_t sizes[3] = {64, 96, 128};
        for (int i = 0; i < n_passes; i++) {
            const int at = passes[i][0];
            // No scratch for a size: the caller runs this batch on the CPU.
            if (!scratch || !scratch->lists[at]) return NM_UNSUPPORTED;
            size_t need = 0;
            nomics_seed_scratch_bytes(sizes[at], &need);
            if (scratch->bytes[at] < need) return NM_INVALID;
        }

        // Packed lists hold x < 2^36: every interval of a BWT shorter than
        // that. NOMICS_SEED_PACK=0 selects plain lists (same output).
        bool packed = bwt->seq_len + 2 < (1ull << 36);
        const char* pack_env = getenv("NOMICS_SEED_PACK");
        if (pack_env && strcmp(pack_env, "0") == 0) packed = false;

        nm_seed_params params;
        memset(&params, 0, sizeof params);
        params.primary = bwt->primary;
        memcpy(params.l2, bwt->l2, sizeof params.l2);
        params.seq_len = bwt->seq_len;
        params.use_table = table != NULL;
        size_t table_bytes = 0;
        if (table) {
            if (table->k < 2 || table->k > NOMICS_KMER_MAX_K) return NM_INVALID;
            params.table_k = table->k;
            memcpy(params.level_offsets, table->level_offsets, sizeof params.level_offsets);
            table_bytes = 3 * sizeof(uint64_t) *
                          (size_t)(table->level_offsets[table->k] + (1ull << (2 * table->k)));
        }
        params.n = n;
        params.cap = cap;

        // Extents of every buffer the kernel reads or writes.
        size_t seq_end = 0;
        for (uint32_t r = 0; r < n; r++) seq_end = MAX(seq_end, (size_t)(offs[r] + lens[r]));
        const size_t words_bytes = (size_t)((bwt->seq_len + 127) >> 7) * 16 * sizeof(uint32_t);
        NSUInteger o_words = 0, o_table = 0, o_seqs = 0, o_offs = 0, o_lens = 0, o_out = 0,
                   o_out_n = 0;
        id<MTLBuffer> b_words = nm_bind(ctx, bwt->words, words_bytes, &o_words);
        id<MTLBuffer> b_seqs = nm_bind(ctx, seqs, seq_end, &o_seqs);
        id<MTLBuffer> b_offs = nm_bind(ctx, offs, (size_t)n * 8, &o_offs);
        id<MTLBuffer> b_lens = nm_bind(ctx, lens, (size_t)n * 4, &o_lens);
        id<MTLBuffer> b_out = nm_bind(ctx, out, (size_t)n * cap * sizeof(nomics_intv), &o_out);
        id<MTLBuffer> b_out_n = nm_bind(ctx, out_n, (size_t)n * 4, &o_out_n);
        // Without a table the kernel never reads binding 2; bind anything.
        id<MTLBuffer> b_table = table ? nm_bind(ctx, table->intervals, table_bytes, &o_table) : b_lens;
        if (!b_words || !b_seqs || !b_offs || !b_lens || !b_out || !b_out_n || !b_table) return NM_WRAP;
        id<MTLBuffer> b_scratch[3] = {nil, nil, nil};
        NSUInteger o_scratch[3] = {0, 0, 0};
        for (int i = 0; i < n_passes; i++) {
            const int at = passes[i][0];
            b_scratch[at] = nm_bind(ctx, scratch->lists[at], scratch->bytes[at], &o_scratch[at]);
            if (!b_scratch[at]) return NM_WRAP;
        }

        id<MTLCommandQueue> queue = (__bridge id<MTLCommandQueue>)stream;
        id<MTLCommandBuffer> cmd = [queue commandBuffer];
        id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
        if (!cmd || !enc) return NM_COMMAND;
        [enc setBuffer:b_words offset:o_words atIndex:1];
        [enc setBuffer:b_table offset:o_table atIndex:2];
        [enc setBuffer:b_seqs offset:o_seqs atIndex:3];
        [enc setBuffer:b_offs offset:o_offs atIndex:4];
        [enc setBuffer:b_lens offset:o_lens atIndex:5];
        [enc setBuffer:b_out offset:o_out atIndex:6];
        [enc setBuffer:b_out_n offset:o_out_n atIndex:7];
        // A serial encoder: each dispatch (and the retry pass) starts after the
        // previous one has finished, so the scratch is reused safely.
        for (int i = 0; i < n_passes; i++) {
            const int at = passes[i][0];
            id<MTLComputePipelineState> pso = psos[packed ? 1 : 0][at];
            [enc setComputePipelineState:pso];
            [enc setBuffer:b_scratch[at] offset:o_scratch[at] atIndex:8];
            params.retry_only = (uint32_t)passes[i][1];
            const NSUInteger width = MIN((NSUInteger)128, pso.maxTotalThreadsPerThreadgroup);
            for (uint32_t base = 0; base < n; base += SEED_DISPATCH) {
                params.base = base;
                [enc setBytes:&params length:sizeof params atIndex:0];
                const uint32_t count = MIN(SEED_DISPATCH, n - base);
                [enc dispatchThreads:MTLSizeMake(count, 1, 1)
                    threadsPerThreadgroup:MTLSizeMake(width, 1, 1)];
            }
        }
        [enc endEncoding];

        pthread_mutex_lock(&nm_seed_lock);
        if (nm_seed_last) [nm_seed_last waitUntilCompleted];
        nm_seed_last = cmd;
        *ticket = nm_commit(cmd);
        pthread_mutex_unlock(&nm_seed_lock);
        return NM_OK;
    }
}

// PairHMM is FP64; Metal has no double. It stays on the CPU (NEON).
int nomics_pairhmm(const uint8_t* read_bases, const uint8_t* hap_bases, const double* row_params,
                   const nomics_pairhmm_pair* pairs, uint32_t n, double* out_sum, float* out_ms) {
    (void)read_bases, (void)hap_bases, (void)row_params, (void)pairs, (void)n, (void)out_sum;
    if (out_ms) *out_ms = 0.0f;
    return NM_UNSUPPORTED;
}

int nomics_ksw_extend2_launch(nomics_stream stream, const uint8_t* queries, const uint8_t* targets,
                              const int8_t* matrices, const nomics_extend_job* jobs, uint32_t n,
                              nomics_extend_result* out, nomics_ticket* ticket) {
    (void)stream, (void)queries, (void)targets, (void)matrices, (void)jobs, (void)n, (void)out;
    if (ticket) *ticket = NULL;
    return NM_UNSUPPORTED;
}

// ---- Fused left+right seed extension (fused.metal; launcher ported from fused.cu) ----

#define FUSED_MAX_BAND_TRY 2
#define EXT_MAX_C 8
#define FUSED_ROW_BUCKETS 64  // 16 target rows each; the last is open-ended
_Static_assert(sizeof(nomics_fused_params) == 104, "fused.metal layout");
_Static_assert(sizeof(nomics_fused_job) == 56, "fused.metal layout");
_Static_assert(sizeof(nomics_fused_result) == 112, "fused.metal layout");

typedef struct {
    int w, c, hi;
} nm_bin;

// fused.cu's defaults (measured on the GB10; NOMICS_FUSED_BINS overrides).
static const nm_bin nm_default_bins[] = {
    {8, 1, 8},   {8, 2, 16},  {8, 3, 24},   {8, 4, 32},   {16, 3, 48},  {16, 4, 64}, {16, 5, 80},
    {32, 3, 96}, {32, 4, 128}, {32, 5, 160}, {32, 6, 192}, {32, 7, 224}, {32, 8, 256},
};

typedef struct {
    int ok, sort, threads, nbins;
    nm_bin bins[32 * EXT_MAX_C];
} nm_fused_config;

static int nm_fused_supported(long w, long c) {
    return (w == 8 || w == 16 || w == 32) && c >= 1 && c <= EXT_MAX_C && w * c <= 32 * EXT_MAX_C;
}

static const nm_fused_config* nm_fused_cfg(void) {
    static nm_fused_config cfg;
    static dispatch_once_t once;
    dispatch_once(&once, ^{
        cfg.ok = 1;
        cfg.sort = 1;
        cfg.threads = 32;
        cfg.nbins = (int)(sizeof(nm_default_bins) / sizeof(nm_default_bins[0]));
        memcpy(cfg.bins, nm_default_bins, sizeof(nm_default_bins));
        const char* s = getenv("NOMICS_FUSED_BINS");
        if (s) {
            cfg.nbins = 0;
            const char* at = s;
            long lo = 0;
            while (*at) {
                char* e = NULL;
                const long w = strtol(at, &e, 10);
                if (*e != 'x') { cfg.ok = 0; break; }
                const long c = strtol(e + 1, &e, 10);
                if (*e != ':') { cfg.ok = 0; break; }
                const long hi = strtol(e + 1, &e, 10);
                if (w < 1 || w > 32 || c < 1 || c > EXT_MAX_C || hi <= lo || hi > w * c ||
                    !nm_fused_supported(w, c) || cfg.nbins == 32 * EXT_MAX_C) {
                    cfg.ok = 0;
                    break;
                }
                cfg.bins[cfg.nbins++] = (nm_bin){(int)w, (int)c, (int)hi};
                lo = hi;
                if (*e == ',') ++e;
                else if (*e) { cfg.ok = 0; break; }
                at = e;
            }
            if (cfg.nbins == 0 || lo != 32 * EXT_MAX_C) cfg.ok = 0;
        }
        if ((s = getenv("NOMICS_FUSED_SORT"))) cfg.sort = strcmp(s, "0") != 0;
        if ((s = getenv("NOMICS_FUSED_THREADS"))) {
            cfg.threads = atoi(s);
            if (cfg.threads < 32 || cfg.threads > 256 || cfg.threads % 32) cfg.ok = 0;
        }
    });
    return &cfg;
}

// fused_<W>x<C>, built on first use; nil (status set) on failure.
static id<MTLComputePipelineState> nm_fused_pso(nm_context* ctx, int w, int c, int* status) {
    static id<MTLComputePipelineState> cache[33][EXT_MAX_C + 1];
    static pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_lock(&lock);
    id<MTLComputePipelineState> pso = cache[w][c];
    if (!pso) {
        pso = nm_pipeline(ctx, [NSString stringWithFormat:@"fused_%dx%d", w, c], status);
        // Group collectives assume 32-wide simdgroups.
        if (pso && pso.threadExecutionWidth != 32) {
            pso = nil;
            *status = NM_PIPELINE;
        }
        cache[w][c] = pso;
    }
    pthread_mutex_unlock(&lock);
    return pso;
}

int nomics_fused_extend_launch(nomics_stream stream, const uint8_t* queries, const uint8_t* targets,
                               const nomics_fused_params* params, const nomics_fused_job* jobs,
                               uint32_t n, uint32_t* order, nomics_fused_result* out,
                               nomics_ticket* ticket) {
    if (!ticket || !stream) return NM_INVALID;
    *ticket = NULL;
    if (n == 0) return NM_OK;
    const nm_fused_config* cfg = nm_fused_cfg();
    if (!cfg->ok || !params || !order || params->max_band_try < 1 ||
        params->max_band_try > FUSED_MAX_BAND_TRY)
        return NM_INVALID;
    const nomics_fused_params p = *params;
    // extend_group_v3's preconditions, as fused.cu (callers extend on the CPU
    // when a launch is refused).
    if (p.m < 1 || p.m > 5 || p.o_del < 0 || p.o_ins < 0 || p.o_del > 1 << 16 || p.o_ins > 1 << 16 ||
        p.e_del < 1 || p.e_ins < 1 || p.e_del > 1 << 16 || p.e_ins > 1 << 16 ||
        p.end_bonus_left < 0 || p.end_bonus_left > 1 << 16 || p.end_bonus_right < 0 ||
        p.end_bonus_right > 1 << 16 || p.w < 1 || p.w > 1 << 20)
        return NM_INVALID;
    for (int k = 0; k < p.m * p.m; ++k)
        if (p.mat[k] < -32 || p.mat[k] > 31) return NM_INVALID;
    @autoreleasepool {
        nm_context* ctx = nm_get();
        if (ctx->status) return ctx->status;
        const int nbins = cfg->nbins;
        int bin_of[32 * EXT_MAX_C + 1];
        for (int q = 1, b = 0; q <= 32 * EXT_MAX_C; ++q) {
            while (q > cfg->bins[b].hi) ++b;
            bin_of[q] = b;
        }
        // Counting sort by (bin, descending target-row bucket), as fused.cu.
        const size_t nkeys = (size_t)nbins * FUSED_ROW_BUCKETS;
        uint32_t* start = calloc(nkeys + 1, sizeof(uint32_t));
        uint32_t* at = malloc(nkeys * sizeof(uint32_t));
        if (!start || !at) {
            free(start);
            free(at);
            return NM_ALLOC;
        }
        size_t q_end = 0, t_end = 0;
#define NM_FUSED_KEY(j)                                                                     \
    ({                                                                                     \
        const uint32_t lq_ = (j).lq_len, rq_ = (j).rq_len, q_ = lq_ > rq_ ? lq_ : rq_;      \
        int bucket_ = 0;                                                                   \
        if (cfg->sort) {                                                                   \
            const uint64_t rows_ = (lq_ ? (uint64_t)(j).lt_len : 0) + (rq_ ? (uint64_t)(j).rt_len : 0); \
            bucket_ = FUSED_ROW_BUCKETS - 1 -                                              \
                      (int)(rows_ / 16 < FUSED_ROW_BUCKETS - 1 ? rows_ / 16 : FUSED_ROW_BUCKETS - 1); \
        }                                                                                  \
        (size_t)bin_of[q_] * FUSED_ROW_BUCKETS + (size_t)bucket_;                          \
    })
        for (uint32_t i = 0; i < n; ++i) {
            const nomics_fused_job* j = &jobs[i];
            // Target lengths are cast to int in the kernel; h0 <= 2^24 keeps
            // every score far from int overflow (fused.cu).
            if ((j->lq_len == 0 && j->rq_len == 0) || j->lq_len > 32 * EXT_MAX_C ||
                j->rq_len > 32 * EXT_MAX_C || j->lt_len > (uint32_t)INT_MAX ||
                j->rt_len > (uint32_t)INT_MAX || j->h0 < 1 || j->h0 > 1 << 24) {
                free(start);
                free(at);
                return NM_INVALID;
            }
            if (j->lq_len) {
                q_end = MAX(q_end, (size_t)(j->lq_off + j->lq_len));
                t_end = MAX(t_end, (size_t)(j->lt_off + j->lt_len));
            }
            if (j->rq_len) {
                q_end = MAX(q_end, (size_t)(j->rq_off + j->rq_len));
                t_end = MAX(t_end, (size_t)(j->rt_off + j->rt_len));
            }
            ++start[NM_FUSED_KEY(*j) + 1];
        }
        for (size_t k = 1; k <= nkeys; ++k) start[k] += start[k - 1];
        memcpy(at, start, nkeys * sizeof(uint32_t));
        for (uint32_t i = 0; i < n; ++i) order[at[NM_FUSED_KEY(jobs[i])]++] = i;
#undef NM_FUSED_KEY
        free(at);

        int status = NM_OK;
        NSUInteger o_q = 0, o_t = 0, o_j = 0, o_o = 0, o_out = 0;
        id<MTLBuffer> bq = nm_bind(ctx, queries, q_end, &o_q);
        id<MTLBuffer> bt = nm_bind(ctx, targets, t_end, &o_t);
        id<MTLBuffer> bj = nm_bind(ctx, jobs, (size_t)n * sizeof(nomics_fused_job), &o_j);
        id<MTLBuffer> bo = nm_bind(ctx, order, (size_t)n * sizeof(uint32_t), &o_o);
        id<MTLBuffer> bout = nm_bind(ctx, out, (size_t)n * sizeof(nomics_fused_result), &o_out);
        id<MTLCommandQueue> queue = (__bridge id<MTLCommandQueue>)stream;
        id<MTLCommandBuffer> cmd = [queue commandBuffer];
        id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
        if (!bq || !bt || !bj || !bo || !bout) status = NM_WRAP;
        else if (!cmd || !enc) status = NM_COMMAND;
        for (int b = 0; b < nbins && status == NM_OK; ++b) {
            const uint32_t first = start[(size_t)b * FUSED_ROW_BUCKETS];
            const uint32_t count = start[(size_t)(b + 1) * FUSED_ROW_BUCKETS] - first;
            if (count == 0) continue;
            const nm_bin bin = cfg->bins[b];
            id<MTLComputePipelineState> pso = nm_fused_pso(ctx, bin.w, bin.c, &status);
            if (!pso) break;
            const NSUInteger threads = MIN((NSUInteger)cfg->threads, pso.maxTotalThreadsPerThreadgroup);
            [enc setComputePipelineState:pso];
            [enc setBuffer:bq offset:o_q atIndex:0];
            [enc setBuffer:bt offset:o_t atIndex:1];
            [enc setBytes:&p length:sizeof p atIndex:2];
            [enc setBuffer:bj offset:o_j atIndex:3];
            [enc setBuffer:bo offset:o_o + (NSUInteger)first * sizeof(uint32_t) atIndex:4];
            [enc setBytes:&count length:sizeof count atIndex:5];
            [enc setBuffer:bout offset:o_out atIndex:6];
            [enc setThreadgroupMemoryLength:64 atIndex:0];
            [enc dispatchThreads:MTLSizeMake((NSUInteger)count * (NSUInteger)bin.w, 1, 1)
                threadsPerThreadgroup:MTLSizeMake(threads, 1, 1)];
        }
        free(start);
        // Nothing is committed on failure: the caller's buffers are free again.
        if (status != NM_OK) return status;
        [enc endEncoding];
        *ticket = nm_commit(cmd);
        return NM_OK;
    }
}
