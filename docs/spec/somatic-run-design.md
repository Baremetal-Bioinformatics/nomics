# Tumour-normal FASTQ integration into `nomics run`

Implemented on `somatic/run-pair`, based on `426bdc8308cd2a219ed786e5bf4990ff3b98ced4`.
**Experimental: local orchestration tests pass; the FASTQ-to-sarek reference gate
has not been run.** The existing BAM chain's reference results do not establish
this new boundary. See [implementation report](somatic-run-report.md).

## CLI and identity

```sh
nomics run --samplesheet samplesheet.csv --reference "$REF" --targets targets.bed \
  --dense-sa-cache sa4.GATK.GRCh38.nsa --outdir new-pair-output \
  --germline-resource af-only-gnomad.hg38.vcf.gz \
  --panel-of-normals 1000g_pon.hg38.vcf.gz \
  --pileup-resource af-only-gnomad.hg38.vcf.gz \
  --orientation-order chr17_38877931-38878296,chr1_11981-12351
```

The CSV has exactly the columns `patient,sex,status,sample,lane,fastq_1,fastq_2`
(in any column order). Status is `0` normal, `1` tumour; sex is `XX`, `XY` or
`NA`, consistent across the patient. Multiple lanes retain samplesheet order.
One patient, one normal and one tumour are required. SM is `patient_sample`;
publication prefix and LB are `sample`. This follows sarek's identity convention;
the [3.10.0 input schema](https://github.com/nf-core/sarek/blob/3.10.0/assets/schema_input.json)
defines numeric status and sex values.

Explicit conservative restrictions: no extra columns, omitted statuses, multiline
CSV fields, ambiguous pairs, conflicting sample statuses, duplicate sample/lane
combinations, same-path mates, or unsafe names. Identifiers use ASCII letters,
digits, `_`, `-`, `.`; they cannot start with `-` or `.`. `+` is refused in pair
mode because the existing chain's label validator does not accept it. Quoted
commas and doubled quotes are supported. Relative FASTQ paths resolve against
the process working directory. Errors identify the CSV row; all sheet validation
and both samples' lane/read-group preflights happen before any output is created.
The parser does not silently rename or trim identifiers.

`--samplesheet` excludes `--fastq1/2`, `--lane`, `--sample`, `--prefix`,
`--resume-from`, `--hc-only`, `--emit-gvcf`, `--reports` and `--keep-bams`.
BAM output is required (`--save-output-as-bam false` is refused). The three
somatic resources are required with a sheet and refused without one, as are
orientation and per-sample scatter-order options. `NOMICS_STOP_AFTER` is refused
in pair mode. No HaplotypeCaller runs, matching sarek `tools: mutect2`.

`--orientation-order` has the chain's comma-separated/repeated BED-stem syntax.
If absent, the log explicitly records the scatter-group order chosen. This is
not asserted to reproduce a reference's task-completion order. The standalone
`somatic chain` still requires an explicit order. See sarek's unsorted
[`groupTuple()` f1r2 reduction](https://github.com/nf-core/sarek/blob/3.10.0/subworkflows/local/bam_variant_calling_somatic_mutect2/main.nf).
`--normal-scatter-order` and `--tumour-scatter-order` override the shared
`--scatter-order` independently for BQSR; unspecified orders use the shared
option, then BED order. Each file is checked as a complete permutation.

## Execution and publication

Normal preprocessing returns before tumour preprocessing starts. Both call the
existing pipeline with an internal `stop_after_bqsr` flag, using its existing
BQSR stop boundary to write recal.bam. Every sample's stack, buffers and joined
workers leave scope before the next call. Files, not read buffers, bridge the
samples. The existing htsjdk index writer produces recal.bai. The library
`somatic_chain::Chain::run` then consumes the recalibrated BAMs and the same
200,000 nt/s scatter algorithm. There are no subprocesses in pair execution.
Neither caller nor preprocessing/compression algorithms were modified.

The existing single-sample layout remains unchanged. Pair mode renames its
`mapped/`, `<sample>.md.bam`, `<sample>.recal.bam` and `<sample>.recal.table`
into the publication tree below; metrics are renamed from `.md.metrics` to
`.md.bam.metrics`. Chain publication is reused without modification.

```text
OUT/
  results/
    preprocessing/
      mapped/S/                         # per-chunk BAMs; --save-mapped true
      markduplicates/S/S.md.bam
      recal_table/S/S.recal.table
      recalibrated/S/S.recal.bam
      recalibrated/S/S.recal.bai
    reports/markduplicates/S/S.md.bam.metrics
    variant_calling/mutect2/
      T_vs_N/T_vs_N.mutect2.vcf.gz
      T_vs_N/T_vs_N.mutect2.vcf.gz.stats
      T_vs_N/T_vs_N.mutect2.artifactprior.tar.gz
      T_vs_N/T_vs_N.mutect2.contamination.table
      T_vs_N/T_vs_N.mutect2.segmentation.table
      T_vs_N/T_vs_N.mutect2.filtered.vcf.gz
      T_vs_N/T_vs_N.mutect2.filtered.vcf.gz.filteringStats.tsv
      T/T.mutect2.pileups.table
      N/N.mutect2.pileups.table
  intervals/*.bed
  timings/S.tsv                        # per-sample existing stage timings
  stages.tsv                           # sample/chain/total time, process peak RSS
  work/preprocessing/S/                              # only with --write-intermediates
  work/chain/                          # per-group caller/pileup outputs for gates
```

`S` occurs once for each sample. `mapped/` includes the nomics ordered chunk
manifest; it is not a sarek artifact. Preprocessing work is nested under `work/preprocessing/` so a sample named `chain`
cannot collide with retained chain intermediates. The final output must not exist and its
parent must exist. A sibling `.nomics-pair-PID.incomplete` staging directory is
renamed to OUT only after successful filtering. Failed staging is retained for
diagnosis, with its path printed; it is not resumable. Logs go to stderr. Timing
peaks are cumulative process high-water marks, not individual-stage allocations.

Not produced: mapped merged `S.sorted.bam` and its index; per-chunk mapped indexes;
`S.md.bai`; raw and filtered Mutect2 `.vcf.gz.tbi`; CRAM/CRAI; FastQC HTML/ZIP/text
reports, fastp JSON/HTML/log reports, samtools statistics, mosdepth coverage and
index files, variant QC reports (bcftools/vcftools), MultiQC HTML/data, or Nextflow
`pipeline_info` provenance/trace/report/DAG/version files. MarkDuplicates metrics
are the only QC report published here. Per-group f1r2 archives, VCFs and pileups
are work products retained for the gate, not additional sarek publications.
The complete reference publication inventory is unavailable on this Mac; the
lead should audit it on the gate machine rather than interpret this as full
sarek artifact coverage.

## Resume (not implemented)

A future pair manifest must bind both patient/sample/SM identities, sex/status,
ordered lanes and FASTQ checksums, reference and resource checksums/index versions,
all preprocessing options, per-sample gather/merge permutations and merge seeds,
ordered interval contents and orientation reduction order, plus output checksums
and generation/parent relationships. A caller resume requires two complete,
validated recalibrated BAM/index generations. Filtering requires the complete
merged VCF/stats, priors and contamination/segmentation set. Changing either
sample invalidates every pair-derived stage. Never combine one fresh sample with
one stale sample.

## Gate-machine results (Spark, 2026-11-15, tree `a6c33a0`, CUDA build)

`NOMICS_REQUIRE_FIXTURES=1 scripts/somatic_run_gate.sh ~/omics/runs/somatic/S1-r1 OUT`, with
the single-sample exome gate's settings (20 threads, GPU seeding, dense SA interval 2, the
three known-sites files, `--spill-prep` because a sample's FASTQs exceed 8 GB gzipped) and
sarek's recorded orientation order.

| Boundary | Result |
|---|---|
| HCC1395BL_WES (normal, penalty 4): fastp chunks, 12 mapped BAMs, md.bam, metrics, recal table, recal.bam | 17 of 17 identical to sarek |
| HCC1395_WES (tumour, penalty 3): the same boundaries | 17 of 17 identical to sarek |
| Scatter groups from the target BED | two groups, names and bytes identical |
| Chain: 16 comparisons | 16 of 16 identical (6,351 records, 2,024 PASS) |

FASTQ to filtered VCF in one process: 58 min 23 s, whole-process peak 59.5 GiB.

| Stage | Normal | Tumour |
|---|---:|---:|
| prep | 57.5 s | 69.0 s |
| align | 247.4 s | 319.5 s (seed extension and mate rescue on the CPU) |
| markdup | 62.1 s | 73.4 s |
| bqsr | 40.9 s | 49.7 s |

The chain took 2,512 s, serial. sarek's reference run of the same pair took 4 h 58 min on the
same machine; that is context, not a like-for-like claim.

Later the same day, on an integration tree (`23a4272`: the chain calling the contigs of a
group in parallel, and GPU scoring at penalty 3), the same gate passed in 21 min 10 s at a
53.5 GiB peak: normal 444 s, tumour 499 s (alignment 275 s), chain 326 s. About 0.041 kWh
at the wall against 0.064 kWh for the run above.

Two earlier attempts are part of the record. The first was killed at 115 GiB because the gate
script ran the 88M-pair normal in memory with six threads; the script now uses the exome
gate's settings. The second completed and failed on the tumour's 12 mapped BAMs: sarek aligns
a sample whose status is tumour with `bwa mem -B 3`, and the aligner then implemented only
the default of 4 (next section).

### Per-sample BWA mismatch penalty

Like sarek 3.10.0 `conf/modules/aligner.config:62-63`, preprocessing sets B=4
for status 0 (normal) and B=3 for status 1 (tumour). Each sample's run log
records its status and selected penalty. Pair mode refuses an explicit
`--bwa-mismatch-penalty`; single-sample mode accepts 3/4, default 4. B=3 uses
CPU extension and CPU mate rescue even in CUDA/Metal builds; GPU seeding can
still run. See `docs/spec/bwa.md` for the source audit and `docs/GATES.md` for
the two-chunk fixture gate. Full tumour/pair equivalence must be rerun on the
gate machine; local fixture equality alone does not establish it.
