# nomics

[![ci](https://github.com/Baremetal-Bioinformatics/nomics/actions/workflows/ci.yml/badge.svg)](https://github.com/Baremetal-Bioinformatics/nomics/actions/workflows/ci.yml)

**The nf-core/sarek germline pipeline, reimplemented from scratch: the same records at every
stage, 11× faster on one NVIDIA DGX Spark.**

If your lab validated its germline calls on sarek 3.10.0, nomics gives the same reads, the
same duplicate marks, the same recalibration and the same variant records from the same
FASTQs, in an hour instead of half a day, on one desktop-sized machine. It is a single
binary: no Nextflow, JVM or containers.

| HG002 30× genome, FASTQ → VCF, one DGX Spark | Wall time | Same records as sarek? |
|---|---:|---|
| sarek 3.10.0 (Nextflow, Docker) | 11 h 24 min | — |
| **nomics** | **1 h 00 min** | **Yes: 41/41 artifacts, all 4,773,797 VCF records** |

Same box, one job at a time, page caches evicted before each run. sarek was measured on
2026-10-17/18 and 2026-10-23, nomics on 2026-10-27. The nomics time is its first
run after a reboot. On a machine that had been up for eleven days, thirteen consecutive
default runs took 64.5–65.8 min, each with identical output. For the exome (Agilent v7, 40.9M
pairs) the times are 1 h 08 min for sarek and **4 min 58 s** for
nomics, 18/18 artifacts identical. Details and history: [RESULTS.md](RESULTS.md).

**Not only HG002.** Six more exomes and three more genomes pass the same gates, each against
its own pair of sarek runs: seven individuals in all, three exome kits and four instruments
([RESULTS.md](RESULTS.md#exactness-across-samples)).

**Which sarek.** The match is to sarek as run on arm64, and to sarek on x86 with GATK's Java
PairHMM and Smith-Waterman and HotSpot's libm intrinsics turned off; those two agree record
for record. sarek's default x86 configuration (Intel GKL) calls about 2% of 30× records
differently from both ([docs/REFERENCE-VARIABILITY.md](docs/REFERENCE-VARIABILITY.md)).

## Why it matters

- **Speed without revalidation.** GPU pipelines are fast, but their output differs from
  the reference, so moving to one means validating it again. nomics produces the reference's
  output, so the validation you already have still applies.
- **One box, one binary.** The reference index and the GPU's working buffers live in the
  DGX Spark's 128 GB of unified memory, where CPU and GPU read them in place. Exome runs keep
  every intermediate in memory; whole genomes stream between stages through local disk.
- **Proof at every stage, not just at the end.** Each stage boundary (trimmed reads, aligned
  chunks, duplicate-marked BAM, recalibration table, recalibrated BAM, VCF) has a written
  gate. Tables and metrics are byte-identical; BAMs and VCFs are record-identical, because
  their compression and run-provenance header lines differ even between two sarek runs.
  Every normalization is written down in [docs/GATES.md](docs/GATES.md).

## How exactness is achieved

The short version (the long one is [docs/METHODOLOGY.md](docs/METHODOLOGY.md)):

1. **The reference tools' source code is the spec**, bugs included. Tie-breaks, rounding,
   random-number use and float formatting are read from BWA, GATK, htsjdk and samtools source
   and confirmed against real outputs, never assumed. The per-tool specifications are in
   [docs/spec/](docs/spec/).
2. **Measure the reference before gating against it.** Several reference outputs are not
   stable from run to run (`@PG` header order, `samtools merge`'s random ID suffixes, the
   VCF's command-line header). The exome goldens were run twice; the gate is whatever is
   stable, and the rest is documented
   ([docs/REFERENCE-VARIABILITY.md](docs/REFERENCE-VARIABILITY.md)).
3. **Gate every boundary on a ladder of inputs:** a small exome slice, the full exome, a
   chromosome-20 slice, 1/12 of a genome, then the whole 30× genome. Then widen it: other
   individuals, exome kits, instruments and depths, each with its own reference runs.
4. **Write it fresh.** No bioinformatics or compression libraries; nomics has its own
   inflate, a zlib-exact deflate, BAM/CRAM/VCF codecs and indexers. Every GPU kernel has a
   CPU version, and the two must agree exactly.
5. **GPU paths are judged against our own saturated CPU path** and ship only when they beat
   it; the end-to-end numbers above are against sarek.

What we learned along the way, from reference tools that do not reproduce themselves to
memory fragmentation that costs more than a kernel rewrite saves:
[docs/LEARNINGS.md](docs/LEARNINGS.md).

## Status

| Pipeline | Status |
|---|---|
| sarek germline, single sample (VCF, GVCF, CRAM or BAM output) | **Exact** on all test tiers up to a 30× genome, and on eleven datasets from seven individuals |
| sarek `--joint_germline` (GenotypeGVCFs) | **Exact**, exome |
| Multi-lane input with optical duplicates | **Exact** |
| sarek's default QC reports (`nomics run --reports`): fastp, FastQC data, samtools stats, mosdepth, bcftools, vcftools | **Exact** but for FastQC's rendered pages; BAM output only, not yet tuned for speed |
| sarek tumour-normal, tool by tool (`nomics somatic …`): Mutect2, the pileup and contamination tools, LearnReadOrientationModel, FilterMutectCalls, the two merges | Experimental: each tool **exact** against GATK 4.6.2.0 on one exome pair; Mutect2 on nine chromosome slices of it (two held out while it was written). Not yet one run, and Mutect2 holds its reads in memory (about 12 GiB per exome chromosome) |
| DGX Spark (GB10, CUDA) | **Primary platform**; all numbers above |
| Apple Silicon (Metal) | Experimental: exact (41/41 at 30×), 1 h 17 min on an M4 Max |
| x86-64, CPU only | Exact at 30× (4 h 46 min, [docs/research/x86-lane-r1](docs/research/x86-lane-r1/)); no GPU path |
| nf-core/rnaseq: Salmon, featureCounts | Experimental: exact on test samples, not yet faster than the reference |
| nf-core/scrnaseq: STARsolo counting | Experimental: exact counting from an oracle BAM; no aligner yet |

Not supported: somatic calling from `nomics run` (the tools above are separate commands),
Mutect2 modes other than one tumour with one matched normal, other references than GRCh38 (the iGenomes GATK bundle),
other versions of the tools above. Design notes, including what could come next, are in
[docs/strategy/](docs/strategy/).

## Quick start

Requirements: Rust (pinned in `rust-toolchain.toml`), and for the GPU path the CUDA toolkit
(`/usr/local/cuda/bin/nvcc`, or set `NOMICS_NVCC`; kernels build for `sm_121`, the GB10, or
set `NOMICS_CUDA_ARCH`). Without `--features cuda` the workspace builds and tests as a
CPU-only build.

```bash
cargo build --release -p nomics --features cuda
```

```bash
target/release/nomics selftest
```

```bash
NOMICS_GPU_SEED=1 target/release/nomics run \
  --fastq1 sample_R1.fastq.gz --fastq2 sample_R2.fastq.gz \
  --reference /path/to/GATK.GRCh38 \
  --targets wgs_calling_regions_noseconds.hg38.bed \
  --dense-sa-cache /path/to/sa2.GRCh38.nsa --dense-sa-interval 2 \
  --sample HG002 --prefix HG002 --outdir out/
```

`--reference` is a directory laid out like iGenomes `GATK.GRCh38`; `--targets` is the BED
sarek gets as `--intervals`. `NOMICS_GPU_SEED=1` turns on GPU seeding (the published runs
use it, with a dense suffix array sampled every 2 positions); the cache is built on first use.

`nomics run --help` lists every option; several mirror a sarek parameter:
`--save-mapped`, `--save-output-as-bam`, and `--emit-gvcf` for sarek's `--joint_germline`. Reproducing the published
numbers, including generating your own sarek goldens and gating against them, is described
in [REPRODUCE.md](REPRODUCE.md).

## Check it yourself

`nomics verify` runs pinned nf-core/sarek 3.10.0 and nomics side by side on a slice of your
own reads, on your own machine, and reports whether every artifact is identical. Docker and
the iGenomes `GATK.GRCh38` directory are all it needs.

```bash
target/release/nomics verify fetch
```

```bash
target/release/nomics verify --fastq1 R1.fastq.gz --fastq2 R2.fastq.gz \
  --reference /path/to/GATK.GRCh38 --targets targets.bed --outdir verify-out/
```

It exits 0 on PASS, 1 on FAIL (with the first differing record of each artifact) and 2 when
no real comparison was possible. [docs/VERIFY.md](docs/VERIFY.md) covers the slice, what is
compared, the report and the current limitations.

## Repository map

| Path | What it holds |
|---|---|
| `crates/nomics` | the `nomics` binary: `run` (FASTQ → VCF) and per-stage commands |
| `crates/nomics-prep`, `-align`, `-markdup`, `-bqsr`, `-hc`, `-joint` | one crate per sarek stage |
| `crates/nomics-io` | BGZF, inflate/deflate, BAM, CRAM, VCF, indexes |
| `crates/nomics-gpu`, `-mem`, `-simd` | CUDA/Metal kernels behind a C ABI, huge-page buffers, NEON kernels |
| `crates/nomics-gate` | the comparison tools behind every gate |
| `crates/nomics-quant`, `-count`, `-solo` | experimental RNA-seq and single-cell ports |
| `docs/` | methodology, gates, per-tool specs, design notes, research |
| `bench/` | measured results per stage and the bake-off |
| `reference/` | how the sarek goldens were produced, oracles, samplesheets |
| `scripts/` | gate and benchmark drivers |

## License

GPL-3.0-only (see [LICENSE](LICENSE)). nomics contains no code copied from BWA, GATK or
htslib, but nomics-align reproduces BWA's algorithms closely enough that we license it, and
the whole project, under BWA's license. Third-party notices: [THIRD_PARTY.md](THIRD_PARTY.md).

nomics is not affiliated with nf-core, the Broad Institute or NVIDIA.
