# nf-core drop-in: sarek with nomics inside it (design, 2026-10-17)

Goal: `nextflow run nf-core/sarek -r 3.10.0 -profile docker,arm64,nomics …` produces the
same published outputs as stock sarek, faster, with no change to the user's samplesheet,
params, reports or downstream tooling.

This is the **bridge**, not the destination. It meets nf-core users inside the workflow they
already run, and it proves concordance there. People who could not get nf-core running at
all are better served by the standalone CLI writing the nf-core output layout (brainstorm
B3): one binary, no Nextflow, no containers, no Wave.

## Mechanism: argv-compatible shims, "exact or delegate"

- `nomics` is a multi-call binary. When it is invoked as `fastp`, `bwa`, `samtools` or
  `gatk` (symlinks in a shim `bin/`), it parses argv exactly as that tool would.
- **Exact or delegate.** If every flag and subcommand is in the supported set, nomics runs
  its own implementation, whose outputs are gated identical to the real tool's. Otherwise it
  `exec`s the real tool with the untouched argv. It never approximates.
- Every invocation logs one line to stderr: `nomics-shim: handled gatk HaplotypeCaller`, or
  `nomics-shim: delegated gatk VariantFiltration (unsupported: --filter-expression)`.
- The real tools stay in the image (the profile extends the stock Wave image), so delegation
  always works.
- Precedent: the RNA-seq gate replayed nf-core's Salmon and featureCounts `.command.sh`
  unchanged with nomics first on `PATH`; `nomics salmon-quant` accepts Salmon's flags.

### sarek germline tasks (tier0/tier1 command lines, from the golden work dirs)

| Process | Command | nomics | Notes |
|---|---|---|---|
| FASTP | `fastp --in1 --in2 --out1 --out2 --json --html --thread --detect_adapter_for_pe --disable_adapter_trimming --split_by_lines N --length_required 15` | `prep` | Chunk files are gated on decompressed content; the gzip bytes differ (fastp's own compressor). JSON/HTML are QC side outputs: delegate or emit compatible JSON. |
| BWAMEM1_MEM | `bwa mem -K 100000000 -Y -R RG -t N $INDEX r1 r2 \| samtools sort --threads N -o out.bam -` | align | v1: the `bwa` shim emits bwa 0.7.19's exact SAM on stdout and `samtools sort` is delegated. v2: the samtools shim recognizes the pattern. The mapped-BAM gate is on records (libdeflate bytes are not reproduced). |
| GATK4_MARKDUPLICATES | `gatk MarkDuplicates --INPUT … --OUTPUT … --METRICS_FILE …` | `markdup` | md.bam rewrites byte-for-byte; metrics identical except lines 2 and 4. |
| GATK4_BASERECALIBRATOR | `gatk BaseRecalibrator --input --output --reference --intervals --known-sites ×3` | `bqsr` | Table byte-identical. |
| GATK4_APPLYBQSR | `gatk ApplyBQSR --input --output --reference --bqsr-recal-file --intervals` | `bqsr` | Records identical. |
| GATK4_HAPLOTYPECALLER | `gatk HaplotypeCaller --input --output --reference --native-pair-hmm-threads --dbsnp --intervals --pcr-indel-model CONSERVATIVE [-ERC GVCF]` | `hc` | VCF/g.vcf identical. One RNG reset per invocation, as GATK does. |
| GATK4_GENOMICSDBIMPORT + GATK4_GENOTYPEGVCFS | workspace dir, then `--variant gendb://…` | `joint-genotype` | Both are ours: the import shim writes a workspace marker naming the g.vcf; the genotype shim reads it. |
| samtools index/merge/stats, mosdepth, FastQC, MultiQC, bcftools, vcftools | | delegate | QC and bookkeeping. |

Per-task shims keep nf-core's orchestration: separate tasks, intermediate files, and
container start-ups. So the drop-in speedup is the kernel share only. That split is useful
evidence in itself. Stock sarek vs sarek + nomics profile vs `nomics run` separates "faster
kernels" from "no orchestration".

## Gates

1. **Task replay harness (no Nextflow, no Docker).** For each task in a golden sarek work
   dir: copy it to scratch, remap the reference/input symlinks to local paths, run its
   `.command.sh` with the shim `bin/` first on `PATH`, and compare every output with the
   golden task's output using the existing comparators (`gate-chunks`, `gate-bam`,
   `gate-report`, `gate-vcf`, metrics diff). Report "handled/delegated" per task. Runs on
   the Mac against `~/omics-fixtures/sarek-tier0/r1`.
2. **Full drop-in run.** `nextflow run sarek -profile docker,arm64,nomics` on tier0, then
   tier1, gated against the goldens with `full_exome.sh`'s gate list.

## Hackability (northstar, unproven)

The sprint's working method is the product's second value proposition:
- exact kernels with oracle dumps and replay benches (ksw, SW, assembly, PairHMM);
- gates at every stage boundary;
- digest kits that verify a run without shipping goldens;
- workers that change core code with the harness deciding what is kept.

If the harness ships as part of the product, the pitch becomes "exact results and optimal
performance for *your* pipeline and *your* hardware": tune or extend a stage, and the gates
tell you in minutes whether you are still identical to the reference. Before claiming it
publicly we need someone other than us to make a change through the harness successfully, and
the harness to be runnable on a user's own data (brainstorm A2).

## Setup tax (motivation)

Getting nf-core running on new infrastructure costs real time: Nextflow and Java versions,
container engine, Wave or registry access, cloud executor config, reference staging. One
data point: half a day to set up sarek benchmarks on AWS. On the Mac today: no Java runtime,
container daemon not running. The standalone CLI removes all of it; the drop-in removes none
of it. That is why both exist.

## Phase 1 status (2026-10-17): shims + task replay, no Nextflow

`nomics` is multi-call: `crates/nomics/src/shim/` dispatches on `argv[0]` (`fastp`, `bwa`,
`gatk`, `samtools`). Each tool has an explicit supported-set table; anything outside it
(flag, value, subcommand, repeated scalar, abbreviation, non-BAM/non-BED input) is delegated
with the untouched argv to the real tool (`NOMICS_SHIM_REAL_PATH`, else `PATH` minus every
directory whose tool name resolves to `nomics`).

| Process | Shim | Notes |
|---|---|---|
| FASTP | handled (`nomics_prep`) | `--split_by_lines 200000000` only; chunk names `<dir of out1>/NNNN.<basename>`. JSON/HTML come from the real fastp run alongside in report-only mode (not reproduced). A chunk over the line limit delegates. |
| BWAMEM1_MEM | `bwa mem` handled, `samtools sort` delegated | bwa 0.7.19's header (`@HD`, `@SQ`, unescaped `-R`, `@PG ID:bwa ... CL:<argv>`) and SAM, streamed per `-K` batch. Plain `Index::load` (no dense-SA cache). |
| GATK4_MARKDUPLICATES | handled | `.md.bai` written with htsjdk's layout. |
| GATK4_BASERECALIBRATOR | handled | |
| GATK4_APPLYBQSR | handled | Table parsed back as `RecalibrationReport` does (`nomics_bqsr::apply_from_report`); a table with other `Arguments` delegates. `.recal.bai` htsjdk layout. |
| GATK4_HAPLOTYPECALLER | handled, VCF and `-ERC GVCF` | `.tbi` htsjdk layout. |
| GENOMICSDBIMPORT / GENOTYPEGVCFS | delegated (`tool ...`) | Not yet: no golden work dir with these tasks is local. |
| samtools *, mosdepth, bcftools, vcftools, FastQC, MultiQC | delegated / untouched | |

Replay: `scripts/nfcore/replay_tasks.py` (tier0 r1 + a synthesized `-ERC GVCF` task from
`gvcf-tier0`): 17 tasks with a handled call, all PASS under GATES.md (12 bwa chunks, MD, BQSR,
ApplyBQSR, HC VCF, HC GVCF). FASTP, FastQC and MultiQC inputs are not on the Mac
(`/home/<user>/omics/data/tier0/tier0_{1,2}.fastq.gz`, sarek's `multiqc_config.yml`).

## Phase 2 status (2026-11-12): sarek under Nextflow with the shims

`scripts/nfcore/dropin_run.sh tier0|tier1` runs the pinned sarek 3.10.0 with
`reference/nomics-dropin.config`, which does one thing: it puts the shim directory first on
`PATH` inside every task container (Nextflow's `env` scope; the directory is under the volume
Nextflow already mounts, so no image is built). The binary is a CPU-only nomics; it runs in
all six task images as they are. Three gates, all in the script:

1. **Shim gate.** From the trace and each task's stderr: every task of FASTP, BWAMEM1_MEM,
   MARKDUPLICATES, BASERECALIBRATOR, APPLYBQSR and HAPLOTYPECALLER must have logged
   `nomics-shim: handled`. A run where tasks quietly delegate would pass the other two gates
   and prove nothing.
2. **Output gate.** `reference/compare_runs.sh` against the stock run: records and headers
   for BAMs and VCFs, table bodies; `mapped/*.sorted.bam` is not gated, as in `docs/GATES.md`.
3. **Reports gate.** Every file of the stock run's `reports/` tree.

| | tier0 | tier1 exome |
|---|---|---|
| Expected tasks handled by nomics | 17 of 17 | 17 of 17 (HaplotypeCaller twice: see below) |
| Main outputs | 17 of 17 | 17 of 17 |
| Report files | 28 of 28 | 28 of 28 |
| Wall time, stock sarek (two runs) | 2 min 14 s, 2 min 16 s | 68 min 00 s, 68 min 26 s |
| Wall time, sarek with the shims | 2 min 04 s | 25 min 54 s |
| `nomics run`, same input, GPU | 37 s, QC reports included | 4 min 38 s |

At tier1 the three arms are 4,080 to 4,106 s, 1,554 s and 278 s: the kernels alone take sarek from 68 to
26 minutes, and removing the orchestration and adding the GPU takes it from 26 to under 5. The
middle arm is not a controlled kernel measurement: it is CPU-only, loads the BWA index once per
chunk task, and the GATK shims use every core of the machine whatever the task was given (only
`bwa -t` and, since this phase, HaplotypeCaller's `--native-pair-hmm-threads` bound them).

Task time by process at tier1, stock against drop-in (seconds, completed attempts):

| Process | Stock | Drop-in | |
|---|---|---|---|
| GATK4_HAPLOTYPECALLER | 1,587 | 171 | handled |
| BWAMEM1_MEM (12 tasks) | 1,075 | 898 | `bwa mem` handled, `samtools sort` real |
| GATK4_BASERECALIBRATOR | 518 | 57 | handled |
| GATK4_MARKDUPLICATES | 358 | 57 | handled |
| GATK4_APPLYBQSR | 355 | 34 | handled |
| FASTP | 45 | 175 | handled, and slower: the report statistics are collected on one thread |
| MERGE_BAM, FASTQC, SAMTOOLS_STATS, MOSDEPTH, INDEX_MERGE_BAM | 621 | 622 | real tools |

What the first tier1 runs found, all from sarek's per-task memory limits, which `nomics run`
never meets:

- **MarkDuplicates** held the whole record set and was killed at 60 GB (exit 137). The shim now
  uses the file-based path (`mark_files`, read ends only in memory): 13.7 GiB peak.
- **fastp** held every chunk and was killed at 4 GB. Nextflow's retry asked for 20 threads,
  which the shim does not support, so the real fastp ran and wrote 20 chunks instead of 12:
  the run finished but could not be compared with the stock run. The shim now writes chunks
  through disk and fits in 4 GB.
- **HaplotypeCaller** is killed at 12 GB on its first attempt and passes on sarek's automatic
  retry at 24 GB (about two minutes lost). nomics loads one interval group's reads at a time,
  and an exome is one group. Bounding that is a change to the caller, not the shim; it would
  also lower `nomics run`'s own peak. Until then the drop-in costs one retry here.
- Two shims in one pipeline (`bwa mem | samtools sort`) share stderr, and their log lines
  interleaved; the shim now writes its line in one piece.

Not done: the QC tools (samtools stats, mosdepth, bcftools, vcftools, FastQC) are still the
real ones under the drop-in, although nomics reproduces them; a whole genome; a GPU arm.
