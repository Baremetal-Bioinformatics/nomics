# Single-cell decision run: do STARsolo's count matrices repeat? (2026-10-15)

**Verdict: go.** Two runs of nf-core/scrnaseq 4.2.0's STARsolo command on 10x PBMC 1k v3 produced byte-identical count matrices.

## Setup, and one deviation

- nf-core/scrnaseq 4.2.0 `--aligner star --protocol 10XV3` on the Spark, with GRCh38 primary and GENCODE v49.
- Genome preparation ran in the pipeline: the GTF gene filter, then `STAR_GENOMEGENERATE` in the pipeline's own container.
- **`STAR_ALIGN` cannot run on arm64 as shipped.**
  - The `arm64` profile builds containers from conda through Wave, with strategy `conda,container`. Bioconda has no `linux-aarch64` build of `star=2.7.10b`, so the build fails ([Wave build report](https://wave.seqera.io/view/builds/bd-2575d215d4a24f44_7)).
  - The pipeline's `containers_docker_arm64.config` is empty.
  - The module's x86 container, `quay.io/biocontainers/star:2.7.10b--h9ee0642_0`, has no arm64 image.
- **What we ran instead:** STAR 2.7.10b built from the tagged GitHub source for arm64 (`local/star:2.7.10b-arm64`, built with `make STAR CXXFLAGS_SIMD=""`), executing nf-core's exact logged `STAR_ALIGN` command.
  - Inputs were staged from run 1's work directories: the index, the filtered GTF, the pipeline's 10x v3 whitelist and the two-lane FASTQs.
  - The command used `--runThreadN 12`, `--soloType CB_UMI_Simple`, `--soloCellFilter CellRanger2.2 1000 0.99 10`, `--twopassMode Basic`, and a sorted BAM.
  - STAR reported "Genome version is compatible".
- Two runs, s1 and s2, took 507 s and 510 s.
- Scripts: `~/omics/logs/solo_decision.sh`. Outputs: `~/omics/runs/expr-decision/solo/`.

## Result

| Output | s1 vs s2 |
|---|---|
| `Solo.out/Gene/raw/{matrix.mtx,barcodes.tsv,features.tsv}` | **identical** |
| `Solo.out/Gene/filtered/{matrix.mtx,barcodes.tsv,features.tsv}` | **identical** |
| `Summary.csv`, `Barcodes.stats`, `Features.stats`, `UMIperCellSorted.txt` | **identical** |
| `SJ.out.tab`, all four `Signal.*.bg` bedGraphs | **identical** |
| `Aligned.sortedByCoord.out.bam` | differs (thread-dependent order, as for bulk STAR) |
| `Log.out`, `Log.final.out`, `Log.progress.out` | differ (timings) |

The pipeline gzips the `.tsv` and `.mtx` files afterwards. Those `.gz` files are gated after decompression, because gzip headers carry mtimes.

## What it means

- The roadmap's first single-cell milestone, STARsolo counting with identical matrices on PBMC 1k, has a stable target.
- The golden for that milestone is a source build of the pinned STAR version, because no packaged arm64 binary exists. Two options:
  - accept the source build as the arm64 reference, as the bwa golden was an ARM build;
  - or confirm the matrices against the x86 bioconda container on the x86 lane.

  The second is cheap and worth doing once, before Phase 1 gates are written.
- The BAM needs the same treatment as bulk STAR: a normalized gate, or a single-thread golden.
