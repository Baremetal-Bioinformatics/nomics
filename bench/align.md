# nomics-align (BWA-MEM 0.7.19 + samtools sort 1.22.1)

Gate: per-chunk `preprocessing/mapped/*.NNNN.bam` records vs the golden (`nomics dev gate-bam`).
Stage gates against an instrumented oracle bwa (docs/GATES.md "Test oracles"): seeds, chains, regions, raw SAM.

## tier0 (12 chunks, 489k pairs)

| | reference (bwa -t 20 \| samtools sort) | nomics align_chunk, 20 threads |
|---|---|---|
| per chunk (40.7k pairs) | 2.1 s realtime at ~530% CPU | load 0.3 s (page cache) + align 0.5 s + sort/write 0.2 s |
| gate | | **IDENTICAL**, all 12 chunks (978,738 records) |

Single-threaded nomics: 6.3 s per chunk (~13k reads/s/core; bwa ~7.4k reads per CPU-second at tier0).

## tier1 (12 chunks, 40.9M pairs, ~10 `-K` batches per chunk)

**IDENTICAL on all 12 chunks (81,580,207 records).** First version: 105-107 s wall per chunk (align ~92 s, sort+write ~15 s),
19.7 GB peak RSS; reference 90 s per chunk (18 min total). Not yet faster: see profiling below.

### Profiling and first optimizations (tier1 chunk 0001)

Per read, single thread, before: seed 56 us, chain (SA lookups + tree) 85 us, extend 23 us; pairing 39 us, of which mate rescue 36 us.
Tier1 reads are far more repetitive than tier0's (20 chains/read vs 4), so per-read cost is ~2.2x tier0's; bwa shows the same ratio.

| change | chunk 0001 wall | gate |
|---|---|---|
| first threaded version | 107 s (align 92 s, sort+write 15 s) | IDENTICAL |
| dense SA cache (every 4th row; 8.0 GB, built once in 252 s) | 79 s (chain 85 -> 16.5 us/read) | IDENTICAL |
| parallel SAM->BAM encode + keyed stable sort | **68.5 s** (sort 2.0 s, write 1.6 s) | IDENTICAL |

Reference: 90 s per chunk. Next: seed search (62 us/read, memory-latency bound; GPU candidate) and mate-rescue SW (36 us/read; Codex vectorizing).

All 12 tier1 chunks re-gated on the dense-SA build (committed HEAD d67073f): **IDENTICAL**, 68.5-69.4 s each,
13.8 min total vs the reference's 18 min; 23.9 GB peak RSS.

### GPU vs CPU for seed search (bench/experiments/occ_microbench.cu)

Dependent random `occ4` lookups over the real 3.2 GB `.bwt` in a huge-page host buffer (LF-walk-like chains):

| | lookups/s |
|---|---|
| CPU, 20 threads | 149 M |
| GPU (GB10, reading host memory in place), 1M-4M chains | 2.36 G (**~16x**) |

Seed search is a dependent chain of ~800 such lookups per read; the GPU hides the latency with parallelism.
Decision: port `collect_intv` to the GPU (CPU path stays the reference; GPU path ships only after the 12-chunk tier1 gate).

| ksw.rs allocation/bounds-check removal (Codex, gpt-5.6-sol; align2 1.79x, extend2 1.18x per call) | chunk 0001 align 63.7 -> 55.3 s (~60 s wall) | IDENTICAL (tier0 x12 + tier1 0001) |

### GPU seed search (`crates/nomics-gpu/shim/seed.cu`)

Line-for-line CUDA port of `seed.rs` (one thread per read; BWT and reads read in place from host memory; up to 64 intervals
per read, overflowing reads recomputed on the CPU). Gate: GPU intervals equal the oracle's `M` dump for all 978,738 tier0 reads.

| step | GPU kernel per 662k-read tier1 batch | chunk 0001 wall | gate |
|---|---|---|---|
| first port (24.8 KB stack/thread) | 2.0-2.8 s, 1.1% fallback | ~58 s (serial: GPU then CPU) | IDENTICAL |
| working lists capped at 64 (6.3 KB stack) | ~2.2 s, 3% fallback | | |
| popcount occ instead of `__constant__` table | ~2.0 s | | |
| **overlap: producer thread reads + GPU-seeds batch N+1 while CPU aligns batch N** | ~1.9 s, hidden | **49.6 s** | IDENTICAL |

Findings:
- The kernel reaches ~3 us/read on tier1 (0.73 us on tier0), ~15x below the dependent-lookup microbench: per-read control flow
  diverges within a warp (SMEM lengths, re-seeding), not memory or constant-cache limits.
- Unified memory is shared bandwidth: while the GPU kernel runs, CPU busy time for the same alignment work rises from ~36 s to
  44-62 s per batch. GPU and CPU stages that both do random access contend; overlap still wins net.

All 12 tier1 chunks on the GPU path (committed HEAD 0fcf205): **IDENTICAL**; 48.9-51.0 s for 10 chunks, 74.6 s and 77.0 s for
0007 and 0012 (concurrent load on the Spark during the run). ~10 min total vs the reference's 18 min.
