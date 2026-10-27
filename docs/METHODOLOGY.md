# Methodology

nomics makes one promise: for the same input, it publishes what nf-core/sarek 3.10.0
publishes. This document describes how that promise is defined, checked and kept while the
code gets faster. The rules came out of the work; [LEARNINGS.md](LEARNINGS.md) collects the
evidence behind them.

## 1. Pin the reference, then measure it

The reference is one pinned pipeline:

- **sarek 3.10.0** (`8ccac7ad`), Nextflow 26.04.6, iGenomes `GATK.GRCh38`. The exome,
  tier0 and chromosome-20 goldens ran on arm64 (`-profile docker,arm64`); the whole-genome
  goldens ran on x86 with HaplotypeCaller in Java mode, which was shown to match arm64
  record for record on chromosome 20. Tool versions are whatever that release runs: fastp 1.1.0, BWA-MEM
  0.7.19-r1273, samtools 1.22.1 (sort) and 1.24 (merge), GATK 4.6.2.0 with htsjdk 4.2.0.
- The command lines, configs, params and samplesheets that produced every golden are in
  [`reference/`](../reference/) (`run_sarek.sh`, `params.*.yaml`, `*.config`).
- **The reference is measured before it is gated against.** The exome and tier0 goldens
  were each run twice, and the r1-vs-r2 diff is the gate's specification: whatever differs
  between two reference runs cannot be a byte gate, and whatever is stable must be
  reproduced exactly. Whole-genome goldens are single runs (12 hours each); their
  order-dependent state is recorded instead (below). The measured differences are in
  [REFERENCE-VARIABILITY.md](REFERENCE-VARIABILITY.md); what sarek actually runs, as opposed
  to what its documentation suggests, is in [REFERENCE-SPEC.md](REFERENCE-SPEC.md).

Two findings shaped everything after: the reference is not always stable (sarek's merge
appends time-seeded random suffixes; task-completion order changes headers and float sums),
and it is not platform-neutral (GATK on x86 and arm64 call different variants). So "the
reference" always means a specific run, recorded with its scatter and merge order
(`reference/*.scatter-order.txt`) and a digest manifest (`reference/*.digest.tsv`).

## 2. Gate every stage boundary

A **gate** compares one published artifact against the golden. The full list, with the
exact comparison for each, is [GATES.md](GATES.md). The policy:

- **Bytes** where the reference is byte-stable (recalibration table, chunk FASTQ contents,
  metrics bodies).
- **Records** where the bytes depend on something outside the pipeline's logic: BAM
  compression by libdeflate, `@PG` header order from task completion, CRAM codec choices.
  The uncompressed record stream must still be identical.
- **Normalizations are written down** with their reason, in GATES.md, and nothing else is
  masked. Typical ones drop a command-line or date header line. No sorting, rounding or
  tolerance is applied to data.
- **Indexes are gated two-sided**: our indexer applied to the golden BAM must reproduce the
  golden index byte for byte, and our index must be that same layout applied to our BAM.
- **Gates check themselves.** A decoder that exits non-zero or decodes nothing fails the gate;
  record counts are printed next to every PASS. (An early harness compared two empty outputs
  and reported PASS.)

Each stage is first gated on the golden's own upstream input, so a defect shows up in the
stage that caused it rather than compounding downstream. End-to-end runs (`nomics run`, FASTQ
to VCF) are then gated at every boundary at once: 18 artifacts on an exome, 41 on a 30×
genome.

## 3. Climb a ladder of inputs

Small inputs find most bugs quickly; only large inputs find the rest. Every change is gated
up a ladder of HG002 inputs, each with its own sarek goldens:

| Tier | Input | Size | Gate |
|---|---|---|---|
| tier0 | chr20 exome targets plus a sample of filtered and unmapped pairs | 0.5M pairs | every stage; runs on a laptop |
| tier0nv | tier0 as two NovaSeq lanes with Illumina read names | 0.5M pairs | 31 artifacts (optical duplicates, per-lane read groups) |
| tier1 | Agilent v7 whole exome, SRR14724482 | 40.9M pairs | 18 artifacts |
| wgs20 | chromosome 20 of the 30× genome | 8.3M pairs | 18, 21 interval groups |
| wgs30s | the first 1/12 of the 30× read pairs, all of the genome | 29.6M pairs | 18 |
| wgs30 | HG002 30× whole genome, SRR14724532 | 350M pairs | 41 |

Large inputs were not optional: the 30× genome exposed a CRAM index edge case and memory
effects that no smaller tier triggered.

## 4. The source is the spec

When a reference output can't be explained, the answer is in the tool's source, not its
documentation. Every port starts with a written specification derived from source and
checked against goldens before implementation: [docs/spec/](spec/) holds them for BWA,
MarkDuplicates, BQSR, HaplotypeCaller, GVCF mode and joint genotyping. They record the
behaviors a reimplementation must copy, including ones that look like bugs:

- unstable sort tie-breaks and data-structure shapes that decide output order;
- random-number streams (seeds, draw order, which code paths consume draws);
- float evaluation (fused multiply-add sites in the compiled reference binary, the libm that
  computed a logarithm, float vs double, Java vs C formatting);
- silent fallbacks (a read-name parser that throws and leaves optical coordinates unset).

Where the reference never writes an intermediate, an **oracle** supplies it: the reference
tool built from its pinned source, shown to reproduce the goldens, and instrumented to dump
internal state (`reference/oracle/`). Oracles are test fixtures only; nothing in nomics
links or calls them.

## 5. Write it fresh, keep two implementations of every kernel

- **No bioinformatics or compression libraries.** nomics has its own inflate, a deflate that
  reproduces zlib's output exactly, and its own BAM, CRAM, VCF and index code. General-purpose
  crates (clap, rayon, libc, mimalloc) are allowed.
- **`unsafe` is confined** to three crates: memory management, the CUDA/Metal FFI, and NEON
  kernels.
- **Every GPU kernel and every SIMD kernel has a portable CPU version**, and the two must
  produce identical results. The CPU version is the reference for the kernel's tests and the
  fallback on machines without the accelerator.
- **The HostBuf rule.** Every large buffer, and every buffer the GPU reads, comes from
  `nomics_mem::HostBuf`: anonymous memory, 2 MiB aligned, `MADV_HUGEPAGE` before first touch.
  The GPU reads these buffers in place. No `cudaHostRegister` and no large `cudaMalloc` (on
  GB10 it is system memory and fragments into small pages).

## 6. Measure honestly

- **Speedups are quoted against our own saturated CPU path**, not against GATK. A GPU path
  ships as the default only when it beats that path on the gated tiers.
- **One heavy job at a time**, under a lock, with page caches evicted and memory compacted
  before each run. Runs record their peak memory per stage, and check how much of the index
  and of the GPU buffers landed on 2 MiB pages; if any fell back to small pages, the run
  writes `timing_suspect.txt` and its time is not quoted.
- **Interleave A/B runs** (A, B, B, A) and repeat chunk-level measurements; only differences
  above the measured noise are claimed. On a long-running machine, memory fragmentation
  moves a 30× run by 300–400 s, more than most single optimizations.
- **Comparisons with other tools run on the same box**, the same day where possible, with the
  same inputs, reference bundle, known sites and flags (`scripts/bakeoff/`).

## 7. Keep only what passes

Every performance change follows the same loop: a hypothesis, a measurement on a gated tier,
and a decision. A change is kept only if it is faster and every gate stays identical; the
rejected ones are recorded with the reason (`bench/`, `docs/strategy/`). The exact path is
never traded for speed: when an optimization would need a different algorithm than the
reference's, it is not an optimization nomics can make.
