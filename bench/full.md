# Full exome, FASTQ → VCF (tier1, HG002 WES, 40.9M read pairs)

`scripts/full_exome.sh tier1` on the Spark (2026-10-09, run `nomics-full-tier1-20261009-1747`). Every stage consumes our own previous output, and every boundary is then gated against the sarek 3.10.0 reference run (`tier1-wes-20261008-r1`).

| Stage | nomics wall | peak RSS | Gate |
|---|---:|---:|---|
| prep (fastp 1.1.0) | 47.6 s | 30.7 GB | 24/24 chunk FASTQs IDENTICAL |
| align, 12 chunks (BWA-MEM 0.7.19 + samtools sort; GPU seeds, dense SA) | 614.3 s (49–65 s each) | 23.5 GB | 12/12 chunk BAMs IDENTICAL |
| preprocess (MarkDuplicates → BaseRecalibrator → ApplyBQSR) | 81.9 s | 54.0 GB | md.bam, metrics, recal table (bytes), recal.bam IDENTICAL |
| HaplotypeCaller (20 threads + GPU PairHMM share) | 41.1 s (29.3 s after a 9.7 s load) | 19.7 GB | 42,358/42,358 records + header IDENTICAL |
| **Total** | **784.8 s = 13 min 5 s** | | **all 19 gates PASS** |

The reference sarek run took **68 min 26 s** wall (`time.txt`), so the end-to-end ratio is 5.2×.

Alignment is now 78 % of the total: 12 chunks run one after another, with the 5.6 GB index and the 8 GB dense SA reloaded from page cache each time. Next targets:
- one process that keeps the index resident;
- feeding prep's chunks straight into alignment, which also drops prep's chunk writes;
- HaplotypeCaller reading `recal.bam` from memory instead of reloading it.

## One process: `nomics run` (2026-10-09, run `nomics-fused-tier1-20261009-r1`)

`scripts/full_exome.sh tier1` now calls `nomics run --write-intermediates` once:
- the index and dense SA are loaded once;
- prep chunks are aligned from memory;
- chunk BAMs, md.bam and recal.bam are written in the background;
- markdup reads the chunk record sets, and HaplotypeCaller reads BQSR's records, from memory.

| Stage | separate commands | `nomics run` |
|---|---:|---:|
| prep | 47.6 s | 29.7 s |
| align (12 chunks) | 614.3 s | 575.4 s |
| markdup + BQSR | 81.9 s | 42.3 + 23.3 s |
| HaplotypeCaller | 41.1 s | 33.7 s |
| **Total** | 784.8 s | **704.4 s = 11 min 44 s** (5.8× sarek) |

Peak RSS is 57.6 GB. All 18 gates are IDENTICAL; `gate-vcf` now covers the VCF header and records.

Loading the index once saved little (about 3 s per chunk), because it came from page cache. Alignment compute itself is now 82 % of the run.

## Alignment ideas evaluated (2026-10-09)

Details are in `docs/research/alignment.md`. All tier1 runs below are gated 18/18.

| Run | prep | align | markdup | bqsr | HC | total |
|---|---:|---:|---:|---:|---:|---:|
| baseline (SAM text, one batch) | 26.0 | 610.6 | 41.7 | 23.2 | 35.7 | 737.2 |
| BAM-direct | 29.1 | 573.8 | 42.4 | 23.4 | 34.7 | 703.5 |
| two batches | 26.4 | 594.7 | 42.8 | 23.3 | 34.8 | 722.0 |
| both | 26.4 | 556.5 | 42.6 | 23.4 | 34.7 | 683.6 |
| both + markdup overlap | 26.2 | 553.1 | 41.9 | 23.5 | 34.7 | 679.3 |
| **+ index loaded during prep, writer widened after alignment** (2 runs) | 32.4 / 35.2 | 528.5 / 524.1 | 42.0 / 42.5 | 23.5 / 23.6 | 35.0 / 35.4 | **660.9 / 660.6 = 11 min 1 s** |
