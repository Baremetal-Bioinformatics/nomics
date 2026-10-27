// Dependent random occ4 lookups over bwa's .bwt: GPU (reading host memory in
// place) vs 20 CPU threads. Decides whether seed search moves to the GPU.
//
//   nvcc -O3 -arch=sm_121 -std=c++17 occ_microbench.cu -o occ_microbench
//   ./occ_microbench <prefix>.64.bwt [chains_gpu] [steps]
//
// Each "chain" mimics an LF walk: occ4 at row k (one 64-byte block read),
// then k' = L2[c] + occ[c] for a base c derived from the counts, so every
// lookup depends on the previous one, as in backward search.
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <chrono>
#include <thread>
#include <vector>
#include <sys/mman.h>
#include <cuda_runtime.h>

struct Bwt {
	uint64_t primary, l2[5], seq_len;
	const uint32_t *w;
};

__host__ __device__ static inline uint64_t occ_aux4(const uint32_t *cnt, uint32_t b) {
	return cnt[b & 0xff] + cnt[b >> 8 & 0xff] + cnt[b >> 16 & 0xff] + cnt[b >> 24];
}

__host__ __device__ static inline void occ4(const Bwt &b, const uint32_t *cnt, uint64_t k, uint64_t out[4]) {
	k -= (k >= b.primary);
	const uint32_t *p = b.w + ((k >> 7) << 4);
	const uint64_t *c = (const uint64_t *)p;
	out[0] = c[0]; out[1] = c[1]; out[2] = c[2]; out[3] = c[3];
	p += 8;
	const uint32_t *end = p + ((k >> 4) - ((k & ~127ull) >> 4));
	uint64_t x = 0;
	for (; p < end; ++p) x += occ_aux4(cnt, *p);
	uint32_t tmp = *p & ~((1u << ((~k & 15) << 1)) - 1);
	x += occ_aux4(cnt, tmp) - (~k & 15);
	out[0] += x & 0xff; out[1] += x >> 8 & 0xff; out[2] += x >> 16 & 0xff; out[3] += x >> 24;
}

__host__ __device__ static inline uint64_t walk(const Bwt &b, const uint32_t *cnt, uint64_t k, int steps) {
	uint64_t o[4], acc = 0;
	for (int i = 0; i < steps; ++i) {
		occ4(b, cnt, k, o);
		int c = (int)((o[0] ^ o[1] ^ o[2] ^ o[3] ^ k) & 3);
		k = (b.l2[c] + o[c] + (k * 0x9e3779b97f4a7c15ull >> 40)) % b.seq_len;
		acc += k;
	}
	return acc;
}

__constant__ uint32_t d_cnt[256];

__global__ void kernel(Bwt b, uint64_t n, int steps, uint64_t *sink) {
	uint64_t t = blockIdx.x * (uint64_t)blockDim.x + threadIdx.x;
	if (t >= n) return;
	uint64_t k = (t * 0x9e3779b97f4a7c15ull) % b.seq_len;
	sink[t] = walk(b, d_cnt, k, steps);
}

int main(int argc, char **argv) {
	const char *path = argv[1];
	uint64_t chains_gpu = argc > 2 ? strtoull(argv[2], 0, 10) : (1u << 20);
	int steps = argc > 3 ? atoi(argv[3]) : 200;
	FILE *f = fopen(path, "rb");
	fseek(f, 0, SEEK_END);
	size_t size = ftell(f);
	fseek(f, 0, SEEK_SET);
	Bwt b;
	fread(&b.primary, 8, 1, f);
	fread(b.l2 + 1, 8, 4, f);
	b.l2[0] = 0;
	b.seq_len = b.l2[4];
	size_t words = size - 40;
	size_t map = (words + (2 << 20)) & ~((size_t)(2 << 20) - 1);
	void *mem = mmap(0, map, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
	madvise(mem, map, MADV_HUGEPAGE);
	fread(mem, 1, words, f);
	fclose(f);
	b.w = (const uint32_t *)mem;
	uint32_t cnt[256];
	for (int i = 0; i < 256; ++i) {
		uint32_t x = 0;
		for (int j = 0; j < 4; ++j)
			x |= (((i & 3) == j) + ((i >> 2 & 3) == j) + ((i >> 4 & 3) == j) + (i >> 6 == j)) << (j << 3);
		cnt[i] = x;
	}

	// CPU: 20 threads, each running chains sequentially.
	unsigned nt = 20;
	uint64_t chains_cpu = 20 * 4096;
	std::vector<uint64_t> sinks(nt);
	auto t0 = std::chrono::steady_clock::now();
	std::vector<std::thread> th;
	for (unsigned t = 0; t < nt; ++t)
		th.emplace_back([&, t] {
			uint64_t acc = 0;
			for (uint64_t c = t; c < chains_cpu; c += nt)
				acc += walk(b, cnt, (c * 0x9e3779b97f4a7c15ull) % b.seq_len, steps);
			sinks[t] = acc;
		});
	for (auto &x : th) x.join();
	double cpu_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
	double cpu_rate = chains_cpu * (double)steps / cpu_s;
	printf("CPU 20 threads: %.1f M lookups/s (%.2f s)\n", cpu_rate / 1e6, cpu_s);

	// GPU: one thread per chain, BWT read in place from host memory.
	cudaMemcpyToSymbol(d_cnt, cnt, sizeof cnt);
	uint64_t *sink;
	cudaMalloc(&sink, chains_gpu * 8);
	for (int rep = 0; rep < 2; ++rep) { // first run warms the GPU's page mappings
		cudaEvent_t e0, e1;
		cudaEventCreate(&e0);
		cudaEventCreate(&e1);
		cudaEventRecord(e0);
		kernel<<<(chains_gpu + 255) / 256, 256>>>(b, chains_gpu, steps, sink);
		cudaEventRecord(e1);
		cudaEventSynchronize(e1);
		float ms;
		cudaEventElapsedTime(&ms, e0, e1);
		cudaError_t err = cudaGetLastError();
		printf("GPU run %d: %.1f M lookups/s (%.2f s) %s\n", rep, chains_gpu * (double)steps / (ms / 1e3) / 1e6, ms / 1e3,
		       err == cudaSuccess ? "" : cudaGetErrorString(err));
	}
	return 0;
}
