# What we learned

Lessons from reimplementing sarek's germline pipeline exactly and making it fast on one
DGX Spark. Each comes with the measurement behind it and the document that records it. The
rules these lessons produced are in [METHODOLOGY.md](METHODOLOGY.md).

## Exactness

### The reference does not always reproduce itself

- **Merged BAMs differ from run to run even with identical reads.** sarek's `samtools
  merge` gives clashing `@RG`/`@PG` IDs random hex suffixes. The "random" generator is seeded
  with the clock (`hts_srand48(time(NULL))`), so we could recover the seed of a golden run
  by search (the task's start second) and reproduce the suffixes.
  ([REFERENCE-VARIABILITY.md](REFERENCE-VARIABILITY.md), [GATES.md](GATES.md))
- **A header difference moves every compressed byte after it.** Two reference runs order
  their 12 `@PG` lines differently (task-completion order). Because BGZF blocks hold a fixed
  65,498 bytes, the whole BAM shifts. Records are identical, so the gate is on records and
  the uncompressed stream. ([GATES.md](GATES.md))
- **Whole-genome outputs depend on which parallel task finished first.** The order of
  `GatherBQSRReports` inputs changes a floating-point sum; the order of `samtools merge`
  inputs changes tie order. A golden therefore includes its recorded scatter order.
  ([REFERENCE-VARIABILITY.md](REFERENCE-VARIABILITY.md))
- **The same GATK calls different variants on x86 and arm64.** On a chromosome-20 test,
  125,121 of 125,153 VCF records agreed. Forcing Java PairHMM and Smith-Waterman on both left
  7 differences; those came from HotSpot's x86 `Math.log/exp/pow` intrinsics. Turning them
  off (`-XX:-UseLibmIntrinsic`) gave 125,153 of 125,153. At 30×, sarek's default x86 run
  (Intel GKL) differs from the Java-mode run in 93,882 of 4,773,796 records.
  ([REFERENCE-VARIABILITY.md](REFERENCE-VARIABILITY.md))
- **Test data can manufacture non-determinism.** A smoke test that reused one FASTQ pair as
  two lanes made MarkDuplicates look unstable (3,858 recalibration-table lines differed).
  Every read existed twice, and which copy survived followed lane-merge order. Real data
  showed no content variability. ([REFERENCE-VARIABILITY.md](REFERENCE-VARIABILITY.md))
- **A harness can pass when both sides fail.** An early comparison hashed two empty outputs
  and printed PASS. Gates now abort on any tool error and print record counts.
  ([REFERENCE-VARIABILITY.md](REFERENCE-VARIABILITY.md))

### Some behavior is only visible in source

- **BWA aligns against random bases where the reference has N.** At index build, every
  non-ACGT base is written as `lrand48() & 3` with `srand48(11)`, and `bwa mem` never reads
  the N mask. Reads can align across N runs, and `MD` prints the random bases. A port must
  use the built index, not the FASTA. ([spec/bwa.md](spec/bwa.md))
- **BWA's output depends on an unstable sort and a B-tree's shape.** Chains are filtered
  with an introsort under comparators that are not total, and duplicate chains depend on a
  B-tree with nodes of up to nine keys. Both had to be ported exactly, not replaced by
  "equivalent" library structures. Score ties are broken by a hash of the read's index
  within its chunk, so chunk boundaries are part of the spec. ([spec/bwa.md](spec/bwa.md))
- **A QD above 35 is a random number.** GATK replaces it with `30 + 3·nextGaussian()` from
  a fixed-seed `java.util.Random`, and the read downsampler draws from the same stream. The
  port reproduces both streams interleaved; the first exome record where they meet pinned
  the first downsampling draw. ([spec/haplotypecaller.md](spec/haplotypecaller.md))
- **Joint genotyping depends on GenomicsDB's contig order.** GenomicsDB visits contigs as
  `chr1, chr10…chr19, chr2…`, and the QD jitter stream is consumed in that order; `bcftools
  sort` later restores dictionary order. The first chr2 record matched draw 3,627: exactly
  chr1 plus chr10–19. ([spec/joint.md](spec/joint.md))
- **MarkDuplicates silently loses optical-duplicate coordinates on SRA read names.** Its
  name parser throws on `SRR14724482.104665`, the exception is caught, and tile/x/y stay −1,
  which changes tie-breaking. ([spec/markdup.md](spec/markdup.md))
- **Chunking is part of the specification.** BWA's batches and insert-size statistics
  restart at every fastp chunk, and fastp splits per worker thread: the exome becomes 12
  chunks, not the 1 its size suggests. ([REFERENCE-SPEC.md](REFERENCE-SPEC.md))

### Numbers and formatting decide bytes

- **The compiler that built the reference matters.** The golden `bwa` binary (GCC 13.3,
  default `-ffp-contract=fast` on arm64) fuses multiply-adds at 58 sites. Rust never fuses,
  so each site uses an explicit `mul_add`, read from the disassembly.
  ([spec/bwa.md](spec/bwa.md))
- **Different tools in one pipeline need different math libraries.** BWA's pairing score
  needs glibc's `log` and `erfc` to the last bit; GATK needs Java's fdlibm (`StrictMath`),
  which glibc does not match. ([spec/bwa.md](spec/bwa.md),
  [spec/haplotypecaller.md](spec/haplotypecaller.md))
- **Java's arithmetic and formatting are part of the output.** `Math.round(+∞)` cast to
  `int` saturates, which is why an empty quantization bin reports quality 93. GATK's
  quality prior truncates toward zero, so a prior of 36.6 gives the same value for q=36 and
  q=37. Picard prints `1.0` as `1`.
  ([spec/bqsr.md](spec/bqsr.md), [spec/markdup.md](spec/markdup.md))
- **Two tools in a chain can reformat the same number.** After GATK writes a VCF,
  `bcftools sort` re-parses every float INFO value as a 32-bit float and prints six
  significant digits (`3.20` becomes `3.2`); `MergeVcfs` then re-prints QUAL. We reproduce
  both instead of masking them. ([spec/joint.md](spec/joint.md))
- **Byte-exact compression is possible but optional.** Our deflate reproduces htsjdk's
  exome md.bam (5,047,709,513 bytes) and recal.bam exactly. It is the opt-in mode; the
  default deflater keeps the same block boundaries and records, and is faster. ([REFERENCE-SPEC.md](REFERENCE-SPEC.md), [strategy/structural-gap.md](strategy/structural-gap.md))

### Gating

- **Small inputs hide real bugs.** A reduced interval file reproduced the loci under study
  but skipped earlier random draws, changing downsampling elsewhere; only the full run
  decides. A CRAM index edge case (unplaced slices with span 0, which htslib's query then
  skips) appeared only at 30×. ([spec/gvcf.md](spec/gvcf.md), [GATES.md](GATES.md))
- **"First exact match" is not the end.** The exome gVCF matched 1,713,849 of 1,713,854
  records; the last five came from three separate source-level behaviors (a signed
  fragment-length filter re-applied after realignment, annotations skipped for uncalled
  genotypes, a likelihood refreshed at annotation time). ([spec/gvcf.md](spec/gvcf.md))

## Performance on unified memory

### Memory layout matters more than kernel speed

- **Huge pages are worth more than any single kernel change.** GPU seeding is about 7×
  slower when the index sits on 4 KiB pages. On the same 30× genome and build, a run right
  after boot took 3,619 s; the next three, on fragmented memory, took 3,930–4,024 s. The
  last kernel round saved about 200 s. ([research/alignment.md](research/alignment.md),
  [strategy/gpu-alignment-scope.md](strategy/gpu-alignment-scope.md))
- **Evicting caches and compacting memory is not enough.** After a 220 GB copy, only 85% of
  the index reached huge pages despite both. The harness now flags any run whose buffers
  fall back to small pages, instead of reporting its time.
  ([strategy/align-pipeline-design.md](strategy/align-pipeline-design.md))
- **The GPU reads host memory in place at near-device speed, given the right allocation.**
  A host buffer that is 2 MiB-aligned with `MADV_HUGEPAGE` sums at 154 GiB/s from the GB10
  GPU, with no copy and no registration. ([strategy/metal-sprint.md](strategy/metal-sprint.md),
  [crates/nomics-mem](../crates/nomics-mem/src/lib.rs))
- **Unified memory is shared, even when nothing is bandwidth-bound.** The aligner moves
  9–14 GB/s, a tenth of the ceiling, yet a CPU-path chunk took 60% longer next to a
  saturating GPU copy: the GPU wins arbitration and CPU latency rises.
  ([strategy/gpu-extend-fused-design.md](strategy/gpu-extend-fused-design.md))
- **Footprint is a speed lever.** Compacting GPU seed sets and storing records by reference
  cut system peak from 90.9 to 64.5 GB and one 30× chunk's alignment from 324.8 to 294.9 s,
  and made room for a third batch in flight. ([strategy/align-pipeline-design.md](strategy/align-pipeline-design.md))

### Balancing CPU and GPU

- **Moving work to the GPU can make things slower.** The first fused GPU extension was
  exact and cut CPU work per read from ~130 to ~45 µs, but chunk alignment went from 379 to
  423 s: seeding, extension and rescue kernels ran one after another, and the GPU became the
  critical path. ([strategy/gpu-extend-fused-design.md](strategy/gpu-extend-fused-design.md))
- **A faster kernel can change nothing.** Extension kernel v3 was 46% faster, and the chunk
  time did not move until the kernel got a high-priority stream (425 → 368 s): the saved
  time had gone into queueing behind seeding.
  ([strategy/gpu-extend-fused-design.md](strategy/gpu-extend-fused-design.md))
- **The bottleneck moves, so profiles expire.** In nine days alignment went from CPU-bound
  with the GPU half idle, to balanced, to GPU-bound at 90%, and back to CPU-bound after the
  rescue kernel halved its GPU time. A planned GPU chaining kernel would have made the
  GPU-bound phase slower. ([strategy/gpu-alignment-scope.md](strategy/gpu-alignment-scope.md))
- **Splitting work between CPU and GPU only helps when one side is saturated.** On the
  Spark, the best extension split tied the CPU-only path (369.9 vs 370.1 s). On an M4 Max,
  where the GPU was saturated, the same split saved 8%.
  ([strategy/gpu-extend-fused-design.md](strategy/gpu-extend-fused-design.md),
  [bench/metal.md](../bench/metal.md))
- **A GPU PairHMM that loses alone wins as a hybrid.** Bit-exact in double precision, it
  was slower than 20 CPU threads (30.9 vs 21.1 s); a dynamically balanced split gave
  14.1 s. ([bench/hc.md](../bench/hc.md))

### Kernels under exactness

- **Copying the CPU algorithm onto the GPU can be exact and too slow.** The first mate-rescue
  kernel emulated BWA's striped vectors and ran at an eighth of 20 CPU threads; a warp
  wavefront rewrite, exact on the first try, ran at 4.4× the CPU.
  ([research/alignment.md](research/alignment.md))
- **The GPU computes, the CPU decides.** GPU results are keyed and checked against the
  CPU's replay of BWA's control flow; a miss falls back to the CPU. Hundreds of millions of
  lookups per 30× chunk, zero misses, and every kernel is gated against oracle dumps.
  ([strategy/gpu-extend-fused-design.md](strategy/gpu-extend-fused-design.md))
- **Some kernels can only get faster by doing less work.** FM-index seeding is bound by
  random-memory requests: about 2.6 G random requests per second, the memory fabric's
  random-miss ceiling, whatever the request size. Occupancy,
  prefetching and read reordering were all flat; to stay exact, nomics must make the same
  lookups BWA does. ([strategy/gpu-alignment-scope.md](strategy/gpu-alignment-scope.md))
- **A silent wrong-memory read was caught only by the gates.** One run gave 64 consecutive
  reads the seeds of the read 2,048 positions earlier. The inputs involved were the only
  GPU-read buffers on ordinary heap pages; every GPU-read buffer is now a pooled huge-page
  buffer, and a diagnostics mode adds guard pages and canaries.
  ([strategy/align-pipeline-design.md](strategy/align-pipeline-design.md), [DIAG.md](DIAG.md))

### I/O and measurement

- **The gate decides what compression costs.** Requiring htsjdk's exact compressed bytes
  cost a third of the markdup and BQSR stages; gating records instead let a faster
  deflater cut those stages from 675 to 528 s. ([strategy/structural-gap.md](strategy/structural-gap.md))
- **Estimates of I/O savings over-promised four- to five-fold.** A faster intermediate
  format decoded 4× faster per core and changed markdup by 5 s, because the stage was bound
  by rewriting and compression. ([strategy/structural-gap.md](strategy/structural-gap.md))
- **Trimming is nearly free.** fastp-equivalent trimming overlaps alignment and costs about
  0–60 s of a 30× run. ([strategy/structural-gap.md](strategy/structural-gap.md))
- **Noise is bigger than you think, and instruments distort.** One exome binary varied by
  39 s in alignment across runs; per-lookup timing counters cost 15% of Metal alignment;
  nsys billed the same rescue kernel 85.5 and 62.9 GPU-seconds on two runs. Claims under
  ~30 s need repeats. ([research/alignment.md](research/alignment.md), [bench/metal.md](../bench/metal.md))
- **The reference is tail-bound.** sarek's 30× run kept 9.4 of 20 cores busy on average;
  its longest HaplotypeCaller scatter task ran 4 h 02 min against a 51-minute median.
  (sarek's Nextflow trace, 2026-10-17)
