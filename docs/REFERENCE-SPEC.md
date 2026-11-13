# The reference as run: what sarek actually does

Our toolchain must reproduce these behaviours, which are easy to miss from the pipeline's headline description. Every item is to be confirmed from tool source + goldens before it becomes a gate.

## Preprocessing: fastp (runs even with trimming off)

- Sarek runs FASTP whenever `split_fastq > 0` (default 50,000,000). With `trim_fastq=false` the args are `--disable_adapter_trimming --split_by_lines 200000000`.
- **fastp's default filters stay on**: quality filtering (`-q 15`, `-u 40`), N-base limit (`-n 5`), `--length_required 15`. Pairs failing these are dropped before alignment. Exact count TBD from `*.fastp.json`.
- Tier1 produced **12 chunks** (`0001…0012.HG002_WES-L1_R{1,2}.fastp.fastq.gz`), not the 1 chunk that 41M pairs / 50M-read split would suggest. fastp's split is per worker thread. **Chunk membership and order are part of the spec** because BWA's `-K` batches restart at each chunk.
  - TODO: confirm whether chunk contents are deterministic across runs.

## Alignment: BWA-MEM 0.7.19-r1273

- `bwa mem -K 100000000 -Y -R <RG> -t <cpus>` per chunk.
  - `-K` fixes the insert-size batch at 100 Mbp of input, so output is independent of `-t`.
  - Batches start fresh at every fastp chunk.
- Each chunk is piped to `samtools sort` (1.22.1). Chunks are merged in **task-completion order**: see `REFERENCE-VARIABILITY.md`.

## Duplicates and BQSR: GATK 4.6.2.0

- Sarek's `gatk4_gcnvkernel` conda env is replaced by `gatk4=4.6.2.0 htslib=1.21 samtools=1.21` (Wave arm64 build `gatk4_htslib_samtools:6fd46722be032bdd`). gcnvkernel needs MKL.
- On aarch64, Intel GKL native code is not available, so deflate and PairHMM use Java implementations.
- BQSR and HaplotypeCaller are scattered over interval groups derived from `Agilent_v7.bed`. Record the group boundaries with the goldens.

## BAM writing (htsjdk, GATK 4.6.2.0): verified on tier1 r1 goldens

- BGZF: blocks hold exactly 65,498 uncompressed bytes (no flush after the header or at record boundaries), then one partial block, then the 28-byte EOF member.
- Every probed block is **zlib level 2, raw deflate, independent per block** (the Deflater is reset). Recompressing with zlib 1.3 reproduces the bytes. The GATK image's JVM links classic zlib 1.3.x.
- `bin` is 0 for reads without an alignment start (pos -1), not the SAMv1 4680 (`BAMRecordCodec.encode`).
- `nomics-io` decode→re-encode reproduces `md.bam` (81,578,107 records) and `recal.bam` (40,074,075) byte-for-byte in the uncompressed stream, with identical block layout.
- **Byte-identical rewrite (milestone 2):** decoding each golden and re-encoding it with `nomics-io` (BAM codec, BGZF `Profile::HTSJDK`, `HtsjdkDeflater` = zlib-exact level 2 with one reused deflater and the level-0 fallback) reproduces `md.bam` (5,047,709,513 bytes) and `recal.bam` (3,516,306,871 bytes) exactly. Single-threaded: 219 s and 123 s, including decode. `nomics` BAM outputs do this only with `NOMICS_BAM_COMPRESSION=htsjdk`; the default fast deflater keeps the block cuts but not the compressed bytes (`docs/GATES.md`).
- **Fresh vs reused deflater:** in theory, htsjdk's reused deflater can differ from a fresh one when stale window bytes break a match tie at a block's end. Measured on tier1 goldens: 0 of 668,293 blocks differ (`examples/fresh_vs_stateful`). Parallel BGZF compression with fresh per-block state is therefore viable, guarded by a detector for matches that scan past end-of-input (sequential replay for those blocks).
- **ApplyBQSR is run with `--intervals`**, so `recal.bam` contains only reads overlapping the targets (40.1M of 81.6M); unmapped reads are dropped.

## Calling

- HaplotypeCaller, single sample, per-interval-group VCFs merged. `haplotypecaller_filter` is skipped (CNNScoreVariants), so the reference output is the raw merged HC VCF.

- With `--joint_germline`, sarek adds `-ERC GVCF` to HC. The tier0 joint fixture uses one BED group spanning the full Agilent v7 BED, with reads concentrated on chr20. GVCF emits uncovered targets on every contig. Default bands, raw annotations, physical phasing, and RNG differences are specified and gated in `docs/spec/gvcf.md`.
