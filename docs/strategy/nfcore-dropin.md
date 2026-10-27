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
