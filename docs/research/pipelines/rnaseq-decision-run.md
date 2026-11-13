# RNA-seq decision run: does nf-core/rnaseq 3.27.0 repeat itself? (2026-10-15)

**Verdict: go.** The roadmap's condition for RNA-seq was a byte-identical Salmon `quant.sf` across two stock runs, despite STAR's thread-dependent record order. It holds.

## Setup

- nf-core/rnaseq 3.27.0, `-profile arm64`, `--aligner star_salmon`, on the Spark.
- Inputs:
  - ERR188140 subsampled to 20M pairs (the predecessor project tier1);
  - GRCh38 primary assembly with GENCODE v49;
  - the predecessor project's STAR 2.7.11b full index.
- Parameters: `gencode = true`, `skip_bbsplit = true` (set through a config `params` block, because bare flags fail nf-schema validation as strings).
- Two runs, r1 and r2, each taking 44–45 min wall time (2,723 s and 2,605 s).
- Comparison: every published file, r1 against r2, byte for byte, with gzip files compared after decompression.
- Scripts: `~/omics/logs/expr_decision.sh` on the Spark. Outputs: `~/omics/runs/expr-decision/`.

## Result: 222 files same, 195 differ, 4 missing in r2

Most differences are reports: MultiQC data and plots, PDFs, HTML, and `pipeline_info` files. The table classifies the differences in pipeline outputs.

| Output | r1 vs r2 | Cause |
|---|---|---|
| Salmon `quant.sf`, `ambig_info.tsv`, `expected_bias.gz`, `fld.gz`, observed bias | **identical** | Salmon 2.7.0 is deterministic, as its 2.6.0 release notes claim, even on STAR transcriptome BAMs that differ in order |
| featureCounts table, summary and biotype counts | **identical** | |
| STAR `SJ.out.tab` | **identical** | |
| bigWig, dupRadar `dupMatrix`, intercept/slope, DESeq2 `.dds` | **identical** | |
| Salmon `meta_info.json` | differs | run time and peak RSS only |
| STAR `Log.out`, `Log.final.out`, `Log.progress.out` | differ | timings (expected) |
| Picard `markdup.sorted.metrics.txt` | differs on line 4 only | the "Started on" timestamp, already normalized in GATES.md |
| `markdup.sorted.bam` and `.bai` | differ | STAR's thread-dependent record order and transcriptome primary choice, as measured in the predecessor project |
| `samtools stats` | 1 line differs ("bases duplicated" 658,314,215 vs 658,314,228) | Picard flags the same number of duplicates, but different reads, because ties follow STAR's order |
| RSeQC `read_distribution` totals, junction files, inner distance | differ | depend on the BAM above; junction BEDs are identical after sorting (encounter order) |
| **StringTie 3.0.3** `gene.abundance.txt`, `transcripts.gtf`, ballgown | differ | **New finding:** TPM values differ in the last digits (e.g. 20.971356 vs 20.971352); sorted, 37,962 lines differ |
| fq-lint logs, FastQC zips | differ | timestamps |

## What it means for the plan

- **Gateable now as byte gates** against a 16-thread golden: Salmon quantification, featureCounts, `SJ.out.tab`, the bigWig, and DESeq2 inputs.
- **Need a single-thread golden, or a normalized gate:**
  - the STAR BAM, the markdup BAM and samtools stats (STAR order, as known);
  - RSeQC totals.
- **StringTie 3.0.3 is not deterministic at default threads.** Gate it against a `-p 1` golden, or treat it as a documented residual. The brief's Phase 0 already planned a StringTie `-p 1` run.
- The roadmap's first milestone (the post-alignment core on golden STAR BAMs) stands as written.
