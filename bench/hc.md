# HaplotypeCaller

GATK 4.6.2.0 HaplotypeCaller as sarek 3.10.0 runs it (VCF mode, dbSNP, PCR
indel model CONSERVATIVE). Port: `crates/nomics-hc`, driver `examples/hc_call.rs`.

## Exactness (2026-10-09)

| tier | records | result |
|---|---:|---|
| tier0 | 1,047 | header (excl. `##GATKCommandLine`) and all records identical, in order |
| tier1 | 42,358 | header (excl. `##GATKCommandLine`) and all records identical, in order |

Per-stage gates (tier0 oracles from the unmodified jar, `reference/oracle/hc/`):
regions and read stream (tier1 too), finalize, assembly, trimming, read
likelihoods, realignment, calling loop, genotyping (synthetic), annotations,
dbSNP (both tiers).

## Timing (Spark, tier1, 2026-10-09)

| | wall (excl. 10 s load) | notes |
|---|---:|---|
| GATK (sarek trace) | 26 min 27 s | ~1 core (102 % CPU) |
| nomics, first exact port | 22 min 0 s | 1 thread |
| + region pool, ordered tail | 114 s | 20 threads; rayon, 4-read likelihood tasks |
| + finisher thread | 108.7 s | RNG handed over around downsampler draws |
| + 4-lane rolling-row PairHMM | 65.4 s | likelihoods 1,046 s -> 265 s single-thread, bit-exact |
| + parallel isActive, dbSNP prefetch | 46.7 s | 1421 % CPU |
| + 32-read likelihood tasks | 41.9 s | 812 CPU-s of 838 available: 97 % saturated |
| + assembly 2.2x faster (packed kmers, two-row SW) | 35.8 s | 678 CPU-s |
| + hybrid CPU+GPU PairHMM (`--features cuda`) | **29.2 s** | GPU backlog budget 2e9 cells |

Every row produced the golden's 42,358 records, identical and in order.

### How it stays exact while parallel

- **Downsampling and QD share one random-number stream.** The downsampler and QD's jitter both draw from GATK's shared `java.util.Random`, so their relative order is part of the output.
  - Only start positions with more than 50 reads draw (25 on tier1).
  - `traverse` calls `RegionSink::before_draw` before each such draw. Every region dispatched so far is finished first, including its QD draws.
  - The RNG is handed back afterwards through `after_draw`.
- **Region preparation is RNG-free and runs in the pool.** The stages are finalize, assembly, trimming, likelihoods, realignment and the dbSNP prefetch.
- **The rest runs in region order on a finisher thread.** That is the call loop, with its spanning-deletion memory, plus annotations and the QD draw.
- **isActive is precomputed in parallel from all reads.** This is exact wherever no downsampled start position's reads reach. The few loci they do reach are piled up on the traversal thread.

### Where tier1's time goes now

- **The pool is saturated.** The run uses 812 CPU-s against 600 single-threaded. The difference is hardware: 10 of the 20 cores are Cortex-A725s, and all-core clocks drop. Drains no longer waste much, so further gains have to cut work.
- **Single-threaded work by stage:** likelihoods (PairHMM) about 265 s, assembly about 176 s, activity about 45 s, the rest small.
- **Parallel assembly across k-mer sizes was tried and reverted.** It took 52.8 s against 46.5 s at the time. No region's assembly exceeds 0.4 s, so assembly never drives the tail. Nested rayon joins also let a waiting worker steal whole other regions, which delays in-order completion.
- **Likelihood task size:** 4 reads → 46.4 s, 32 → 42.1 s, 64 → 42.4 s, 512 → 56.2 s, unsplit → 66.3 s.

Next levers reduce work: PairHMM on the GPU (it must beat this saturated CPU path) and single-thread assembly efficiency.

### Experiments (2026-10-09)

**Assembly efficiency.** Kept: 2.2x on `asm_bench`, end to end 41.9 s → 35.8 s. Remaining hotspots are Smith-Waterman CIGARs and graph threading.

**GPU PairHMM.** Not shipped: bit-exact, but slower than the CPU path. The run used `pairhmm_bench` on the full tier1 dump: 22.06M pairs, 2.80e11 cells.
- CPU, 20 threads: 21.1 s (1.33e10 cells/s).
- GPU total: 30.9 s (0.91e10 cells/s). The kernel alone takes 26.3 s (1.06e10 cells/s).
- 0 bit mismatches.

The kernel needs f64 in Java's operand order, without FMA. It uses one warp per pair and an anti-diagonal wavefront. A CPU+GPU split could still add throughput, but memory-bandwidth contention has to be measured first.

**Hybrid CPU+GPU PairHMM: shipped** (`--features cuda`; `NOMICS_GPU=0` disables).
- `pairhmm_bench` split mode (GPU share, all 20 CPU threads): 25 % → 16.7 s, 35 % → 14.6 s, 45 % → 14.1 s, against 21.1 s CPU-only. All runs bit-exact. Contention slows the CPU's share by about 22 %. Giving the GPU feeder one of the 20 threads was worse.
- In the pipeline, a region sends its likelihood batch to the GPU thread while the GPU backlog is under a cell budget, and otherwise computes on the CPU. This balances the two dynamically.
- Tier1 by budget: 0.5e9 → 30.0 s, 1e9 → 29.4 s, 2e9 → 29.3 s, 4e9 → 29.9 s, 8e9 → 31.6 s, 16e9 → 34.5 s. CPU-only takes 35.7 s.
- Large backlogs hurt for two reasons: regions wait on the GPU while CPU cores could have computed them, and in-order completion is delayed.

### GVCF mode in parallel (sprint 2, 2026-10-16)

`nomics hc --emit-gvcf` on the tier1 recal.bam (Spark, 20 threads; harness: g.vcf gate + `joint-genotype` + joint VCF gate, all IDENTICAL at every step):

| change | HC wall | peak RSS | ordered tail busy (ref. confidence) |
|---|---:|---:|---:|
| before: every region in order on one thread (`NOMICS_GVCF_SERIAL=1`) | 392.5 s | 31.6 GB | n/a (ref. confidence 111.7 s) |
| S-01 regions prepared on the pool, genotyping + ref. confidence + blocks in order | 134.0 s | 46.3 GB | 116.5 s (91.0 s) |
| S-05 whole-region ref. confidence on the pool when nothing is prepared | 80.2 s | 65.4 GB | 60.7 s (31.4 s) |
| S-06 compact `RefSite` into the block combiner (no per-locus `Record`) | 65.6 s | 43.7 GB | 47.6 s (18.5 s) |
| S-07 flank ref. confidence on the pool too; the tail keeps only variant spans | 45.3 s | 31.6 GB | 23.0 s (6.4 s) |

All rows ran with `NOMICS_TIMING=1`, which on the wgs30s VCF-mode harness costs ~25% of HC wall (286 vs 220 s, same binary; board S-11), so absolute times here are pessimistic; the ratios are like for like.

GVCF annotations draw no QD jitter, so region work needs no RNG hand-off (the downsampler's draws stay in the traversal). Rejected: dropping ended reads from the per-locus scan (S-03), exact but 91.0 → 89.4 s.

### Sprint 2 kernel and allocator changes (wgs30s VCF-mode HC harness, Spark, `NOMICS_TIMING` off)

`nomics run --resume-from hc` on the wgs30s recal.bam (whole-genome intervals), VCF gated against the x86 Java-mode golden: IDENTICAL on every run.

| build | runs | HC wall |
|---|---|---:|
| glibc malloc | 2 | 220.5, 223.4 s |
| + mimalloc (S-10) | 3 | 200.5, 199.8, 198.7 s |
| + NEON wavefront Smith-Waterman for haplotype CIGARs (S-09) | 2 | 157.0, 157.5 s |
| + assembly graph: reusable k-mer table, array sets for vertex IDs, fewer allocations (S-13) | 2 | 150.6, 145.3 s |

Before S-09, `assembly::cigar::smith_waterman` was 28% of HC CPU (perf, wgs30s), almost all via `calculate_cigar`. The wavefront evaluates anti-diagonals (each cell reads only (i-1,j-1), (i-1,j), (i,j-1)) with Java's operation order and tie rules; on the 600k-call dump it is 3.75x on the `calculate_cigar` workload on the Spark (7.05 -> 1.88 s, 1 thread) with 0 mismatches against the reference port.

S-13 on `asm_bench` (8,398 dumped exome regions, `~/omics-fixtures/asm-dump`, checksum identical): graph 6.10 -> 2.92 s, all of assemble 12.94 -> 9.93 s per pass. 30x HC tier 1,015.6 -> 983.1 s (4,773,797 records identical); GATK assembly oracle test passes on both fixture sets.

### 30x whole genome, end of sprint 2 (2026-10-16, binary 4681d4d: GVCF parallel, extend2 loop, mimalloc, SW wavefront; before S-13)

41/41 digests identical. prep 534 s (alongside align), align 4,744, markdup 299, BQSR 393, HC 1,015 s; FASTQ -> VCF 6,467 s (last window's final run 7,135 s); peak RSS 74.5 GB (was 77).
