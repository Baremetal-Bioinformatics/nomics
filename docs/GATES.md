# Gates

What "identical to the reference" means at each stage boundary. A normalization may only appear here with its justification. Nothing else is masked.

| Boundary (golden) | Written by | Gate | Normalization and why |
|---|---|---|---|
| RNA-seq `star_salmon/featurecounts/*.featureCounts.tsv`, `.tsv.summary` | Subread featureCounts 2.1.1 | **bytes** with matching staged filenames; **bytes after path-only provenance normalization** otherwise | The comment embeds quoted argv, and count/summary column headers embed input BAM paths. Only the executable spelling (to `featureCounts`) and the supplied annotation/input/output path tokens may be replaced by agreed staging names; input BAM column headers use those same names. Preserve the version, argument order, all flags and values, quotes, trailing space, annotation columns, counts and every summary row. Do not discard the comment line. No normalization is needed in the five-sample CLI gate, which stages reference basenames. See `docs/FEATURECOUNTS.md`: nf-core's `-p` without `--countReadPairs` counts alignments, and two fixture samples are SE (no `-p`). |
| RNA-seq `star_salmon/featurecounts/*.biotype_counts*_mqc.tsv` | nf-core/rnaseq 3.27.0 custom biotype script | **bytes** | None. All ten fixture files match. |
| `preprocessing/markduplicates/*.md.bam` | GATK MarkDuplicates (htsjdk) | **records + uncompressed stream identical after `@PG` normalization** | The reference itself is not byte-stable. r1 and r2 order the 12 BWA/samtools-sort `@PG` lines differently (task-completion order), and that changes header length. Header blocks span fixed 65,498-byte offsets, so every downstream BGZF block shifts. *Corollary:* our htsjdk-mode writer is byte-exact (both tier1 goldens rewritten identically), so with the reference's own header our compressed bytes equal the golden's. That mode is opt-in (`NOMICS_BAM_COMPRESSION=htsjdk`). By default md.bam, recal.bam and the per-chunk BAMs use the fast deflater (`nomics_io::deflate::FastDeflater`): same 65,498-byte block cuts, records and uncompressed stream, different (smaller) compressed bytes. |
| `preprocessing/recalibrated/*.recal.bam` | GATK ApplyBQSR (htsjdk) | same as above | Same `@PG` ordering effect. |
| `preprocessing/markduplicates/*.md.cram`, `preprocessing/recalibrated/*.recal.cram` + `.crai` (default sarek, `save_output_as_bam false`) | md: `samtools view -Ch -T <fasta>` (samtools 1.21, CRAM 3.0) in GATK4_MARKDUPLICATES. recal: GATK 4.6.2.0 ApplyBQSR writing CRAM from md.cram (htsjdk 4.2.0, CRAM 3.0); one interval group at tier0, so no `samtools merge` (which would write `version=3.0`) | **decoded records and header equal the reference CRAM's, under each decoder**: htslib (samtools 1.24), htsjdk 5.0.0 (Picard 3.5.0) and htsjdk 4.2.0 (GATK 4.6.2.0's ViewSam in sarek's container). Header normalization: `@SQ UR:` (a path) and all `@PG` lines dropped on both sides, as in the md.bam gate (provenance: ours records nomics, sarek's the tools it ran). Our `.crai` equals `samtools index`'s; htslib region and htsjdk interval queries equal the reference's. Every decoder's exit status is checked and an empty decode fails the gate; each region query's count must also equal the source BAM's. Unplaced slices carry start 0, span 1 (as htslib's CRAM 3.0 writer; with span 0 htslib's `cram_index_query` skips all but the last unplaced slice on a `*` query, which the 30x gate caught). References rebuilt on the tier0 golden by `scripts/cram_reference.sh` (sarek's own containers); gate `scripts/cram_gate.sh <bam> <fasta> <reference.cram> sarek-md\|sarek-recal`. tier0 md (977,738 records) and recal (873,049): all pass, both from `nomics dev bam-to-cram` on the golden BAMs and from `nomics run --save-output-as-bam false` (`SAVE_OUTPUT_AS_BAM=false scripts/mac_tier0.sh`, in-memory and `--spill-prep` streaming paths). | Bytes are not gated (htslib's/htsjdk's codec choices and layout are their own). What we reproduce, from their source: **md** (`nomics_io::cram::TagMode::Htslib`): `RG:Z` of a header read group moves to the RG series; `MD`/`NM` are dropped exactly when htslib's `process_one_read` recomputes the same value (case-insensitive MD, integer NM, kept verbatim on an N-on-N base or an alignment off the contig); other tags keep BAM order; so htslib decodes them as `..., MD, NM, RG` and htsjdk 4 (which does not regenerate) shows no MD/NM. **recal** (`TagMode::HtsjdkSorted`): the same drops (ApplyBQSR read them from md.cram), remaining tags sorted by the 3-byte key as htsjdk's `CompressionHeaderFactory` does. Headers: `M5`/`UR` on every `@SQ` (data, as htslib adds them), plus our own `@PG ID:nomics-cram` line; we write no `@PG` for tools we did not run. Mates: we attach only pairs every decoder rebuilds exactly, htslib attaches by its own rule within its own slices; both decode to the true mate fields whenever TLEN follows BWA's 5'-to-5' convention, which tier0 confirms. **Known divergences where the reference is lossy and we are not:** a start-tie pair whose TLEN breaks the 5' convention (htslib's encoder and decoder break the tie differently), and a `Z` value containing a TAB (invalid SAM; htslib truncates it); both shown in `tests/cram_roundtrip.rs`. `TagMode::Lossless` (decodes to the BAM's own records) stays available via `--profile lossless`. **Multi-group recal** (whole genomes): sarek merges the per-group htsjdk CRAMs with `samtools merge --output-fmt cram,version=3.0`. In CRAM mode each group's header is md.cram's, i.e. md.bam's plus the `@PG` lines `samtools view -Ch` adds, then ApplyBQSR's. Those extra IDs take more of merge's seeded `lrand48` draws, so the `@RG`/`@PG` suffixes (and every record's RG:Z/PG:Z) differ from recal.bam's from the third input on. `bam_translate` also re-appends RG and PG on every record. We reproduce both (`nomics_bqsr::CramMerge`: the seed is recovered from recal.bam's suffixes, and the BAM-mode merge is re-derived and checked; `TagMode::HtsjdkSortedMerged`: PG:Z after the other stored tags). This reproduces how sarek renames IDs; it is not a provenance claim. The header carries the groups' `samtools view` `@PG` lines, which the gate drops. Gated on tier0 split into 3 groups against sarek's chain built from the same run's md.bam and table (`scripts/cram_merge_gate.sh`): streaming, in memory, a resume at HC, and dev mode `--groups 3`. All pass, and the pre-fix writer fails. Its IDs reproduce the 30x cloud reference's `@RG`/`@PG` lists exactly; **Whole genome (30×, 21 interval groups; built on an x86 cloud instance):** the reference is sarek's CRAM-mode chain rebuilt from our gated md.bam (md.cram, ApplyBQSR per group writing CRAM, then `samtools merge --output-fmt cram,version=3.0` in the golden's merge order, seeded with the golden merge's `time(NULL)`). **md.cram passes everything** (703,470,853 records and the header, under htslib and htsjdk 4; M5, `.crai`, regions). recal.cram failed before the multi-group fix above, on exactly those two causes (merge ID suffixes; PG:Z position), with M5, `.crai` and region counts passing; with the fix (608843e, run on the lane 2026-10-25 at 2a03e00) **recal.cram passes everything at 30×**: all 676,643,121 records and the header identical under htslib 1.24 and htsjdk 4.2.0, M5, `.crai` (67,678 slices), 14 region queries and their counts against recal.bam. Every artifact sarek publishes by default is now gated at 30×. |
| fastp chunk FASTQs (`work/*/NNNN.*_R{1,2}.fastp.fastq.gz`, not published) | fastp 1.1.0 | **decompressed content, per chunk** | Intermediate files that decide what BWA sees (chunk membership and order). Compression is fastp's own and not reproduced. `nomics prep`: 24/24 identical at both tiers. fastp's JSON, HTML and log reports have their own rows below (`reports/fastp/**`). |
| `preprocessing/mapped/*.sorted.bam` | samtools merge 1.24 (htslib + libdeflate) | not gated | `samtools merge` appends random hex suffixes to clashing `@RG`/`@PG` IDs (`HG002_WES.L1-128D133E`), so it differs every run. It is only a side output (`--save_mapped`); MarkDuplicates reads the per-chunk BAMs directly. |
| `preprocessing/mapped/*.NNNN.bam` (per chunk) | samtools sort 1.22.1 (libdeflate) | **records** (decompressed stream) | Byte-identical across r1/r2. We do not reproduce libdeflate, so the aligner is gated on records. |
| `preprocessing/recal_table/*.recal.table` | GATK BaseRecalibrator | **bytes** | Byte-identical across r1/r2 at tier1 and tier0. `nomics bqsr` output is byte-identical at both tiers (`bench/bqsr.md`). |
| `variant_calling/haplotypecaller/*.vcf.gz` | GATK HaplotypeCaller (htsjdk) | **records identical; header identical except `##GATKCommandLine`** | `##GATKCommandLine` carries the run date and absolute work-directory paths. HaplotypeCaller writes no `##fileDate`, and `##source` is constant. Records were identical across r1/r2 at both tiers. |
| `variant_calling/haplotypecaller/*.haplotypecaller.g.vcf.gz` | GATK HaplotypeCaller `-ERC GVCF` + GVCFBlockCombiner | **records identical; header identical except `##GATKCommandLine`**, same normalization as VCF above | Joint-germline fixtures: tier0 **251,885/251,885** and full-exome tier1 **1,713,854/1,713,854** records identical, with normalized headers identical and zero mismatches. Tier0 includes 250,733 hom-ref blocks and 1,152 variant records (105 physically phased). Independent reference runs agree except gzip header and `##GATKCommandLine`. `scripts/gvcf_tier0_gate.sh` / `nomics dev gate-gvcf` compare all records and header fields; no block expansion, allele reordering, or annotation masking. See `docs/spec/gvcf.md`. |
| `variant_calling/haplotypecaller/joint_variant_calling/joint_germline.vcf.gz` | GATK GenomicsDBImport + GenotypeGVCFs → `bcftools sort` → GATK MergeVcfs | **records identical; header identical except `##GATKCommandLine`** | Tier1 x86 golden pair (r1 == r2): **42,362/42,362** records identical and the normalized header identical, from the golden gVCF and from nomics' own HC gVCF. `nomics run --emit-gvcf` writes it as `joint_variant_calling/joint_germline.vcf.gz` from the in-memory gVCF (2.0 s of an 866 s FASTQ → joint VCF run). `scripts/full_exome.sh tier1joint` gates it with `nomics dev gate-vcf`; `nomics dev digest` lists it as `joint_vcf`. Standalone: `nomics joint-genotype`. The command line embeds work-directory paths. Float fields carry bcftools' float32/`kputd` re-encoding and are reproduced, not masked. See `docs/spec/joint.md`. |
| `reports/vcftools/**/*.TsTv.qual` | vcftools 0.1.17 `output_TsTv_by_quality` | **bytes except explicitly exempt fields**; tier0 and tier1 pass | Column 5 (`N_Ts_GT_QUAL_THRESHOLD`) is exempt only where nomics's own column 5 is `0`; column 6 (`N_Tv_GT_QUAL_THRESHOLD`) is exempt only where nomics's own column 6 is `0`; column 7 is exempt only where nomics's own column 5 or 6 is `0`. All other fields, all of columns 1–4, the header, row count, separators and newlines remain byte-identical. nomics computes with the invalid boundary element taken as `0.0`, including ordinary double division and six-significant-digit default stream formatting (`0/0` prints `nan`). Why: upstream allocates N+1 elements at [variant_file_output.cpp:3324](https://github.com/vcftools/vcftools/blob/v0.1.17/src/cpp/variant_file_output.cpp#L3324), then reads element N+1 at [3331–3332](https://github.com/vcftools/vcftools/blob/v0.1.17/src/cpp/variant_file_output.cpp#L3331) and [3344](https://github.com/vcftools/vcftools/blob/v0.1.17/src/cpp/variant_file_output.cpp#L3344). These undefined reads expose heap contents only in zero-count fields on the supplied goldens. Five fields exempted on each tier: all five differ on tier0; all five already equal on tier1. `nomics_qc::compare_tstv_qual` implements this asymmetric comparison and returns every exempt field, including equal ones. Evidence and exact coordinates: `docs/spec/qc-vcf.md`. |
| `reports/vcftools/**/*.TsTv.count`, `reports/vcftools/**/*.FILTER.summary`, `reports/bcftools/**/*.bcftools_stats.txt` | vcftools 0.1.17; bcftools 1.23.1 + htslib 1.23.1 | **bytes**; six tier0/tier1 comparisons pass | None, but the vcftools files depend on the architecture sarek ran on: its Ts/Tv ratios go through glibc, which prints a NaN's sign, and 0/0 is `-nan` on x86 and `nan` on arm64 (the four x86 goldens against the arm64 ones). nomics writes the spelling of `--reference-arch`, by default the machine it runs on, so it matches sarek on that machine. For `.TsTv.count`, implement only the source's visible indices 0 through 2×samples−1, treating its invalid maximum-AC write ([3195/3197](https://github.com/vcftools/vcftools/blob/v0.1.17/src/cpp/variant_file_output.cpp#L3195)) as having no effect; this changes no golden byte. bcftools provenance retains its pinned version and the input VCF basename (independent of the output prefix), including the command's two spaces. |
| `reports/mosdepth/**/*.mosdepth.global.dist.txt`, `*.mosdepth.region.dist.txt`, `*.mosdepth.summary.txt` | mosdepth 0.3.14 | **bytes**; all 18 comparisons pass. Wider gate 2026-11-12 (`scripts/qc_bam_gate.sh` on the Spark): every sarek BAM that still has its reports, 130 BAMs in 64 runs, 130 pass (129 have mosdepth reports: 110 in the exome form over five target BEDs, 19 in the genome form up to 30x) | None. Same BAMs as sarek; BED targets for exomes, windows with `--no-per-base --fast-mode` for genomes. Retain all rounding, histogram, contig-inclusion and summary quirks documented in `docs/spec/qc-mosdepth.md`. |
| `reports/mosdepth/**/*.per-base.bed.gz`, `*.regions.bed.gz` | mosdepth 0.3.14 via hts-nim/htslib | **decompressed bytes + identical uncompressed BGZF block-size sequence**; all 10 comparisons pass; and in every BAM of the wider gate (see the row above) | Decode BGZF only; no text normalization. Each full block has **65,280** uncompressed bytes; lines/newlines may cross blocks; retain the partial final block and standard empty EOF member. Compressed bytes are not reproduced: all 100 sampled tier0 md per-base blocks and all 95 nonempty region blocks match libdeflate level 1 exactly, versus zero matches at zlib levels 1–9. nomics uses its existing `FastDeflater` and the same block cuts; no new compressor or `nomics-io` change. |
| `reports/mosdepth/**/*.bed.gz.csi` | hts-nim BGZI/CSI wrappers + htslib index builder | **hts-nim/htslib layout, two-sided**, as for `.bai`; all 10 gates pass, and both sides in every BAM of the wider gate; CSI payload compared after BGZF decoding | (1) `nomics qc csi` applied to the **golden** BED reproduces the **golden CSI payload**. (2) our CSI is the same builder applied to **our BED** (complete compressed bytes also compared). Report index BGZF block-size sequences separately. Min shift 14; levels from largest BAM-header contig (5 on GRCh38); columns 1/2/3, zero-based BED preset; push after each line; flush before finish; retain metadata bins, khash bin order and compressed-distance compaction. Shared BED block cuts preserve virtual offsets within corresponding boundaries, but do **not** ensure identical index shapes: htslib compacts a bin into an existing parent when its **compressed byte span <65536**. WGS20 md/recal first cross-file payload difference: byte **79402**, chr1 `n_bin`, ours **32**, golden **34**; payload sizes **380010 / 380634**. Thus cross-file CSI block-size and within-block-field sequences can differ, even though (1) and (2) are exact. This is the stated two-sided gate, not a field exemption; see `docs/spec/qc-mosdepth.md`. |
| `reports/fastqc/**/*_fastqc.zip` members `fastqc_data.txt`, `summary.txt` | FastQC 0.12.1, default raw FASTQ per mate | **bytes**. Wider gate 2026-11-12 (`scripts/qc_reads_gate.sh` on the Spark): every sarek samplesheet row whose FASTQs are still on disk, 60 of 95 rows, in 56 runs (exomes of seven samples, 14 synthetic sets, a tumour-normal pair, one 20x genome), both mates, both files: 240 of 240 identical. The 35 rows not run: 31 whose FASTQs are gone or were staged on another machine (the three matrix genomes, the x86 and CI runs) and 4 with no reports (two failed runs and a directory that holds only a samplesheet). Two of the x86 runs were then gated from the Spark's own copies of their FASTQs: the HG002 30x genome (355.8 M pairs) and its slice, both mates, both files identical, and all three fastp files identical under the rows below | None. Retain the staged input basename (`--input-name`), module order/status, all numbers, tabs and final newline. Images, HTML, `.fo` and ZIP framing are outside this boundary. See `docs/spec/qc-fastqc.md`. |
| `reports/fastp/**/*.fastp.json` | fastp 1.1.0, `--thread 12` | **bytes except the number at `duplication.rate`**. Wider gate 2026-11-12 (`scripts/qc_reads_gate.sh` on the Spark, the 60 rows of the FastQC row above): 47 of 60 identical with nothing exempted, 13 identical but for that number; the 30x genome and its slice, against x86 goldens, identical but for that number (0.0495891 against 0.049589) | The number is replaced by a fixed token on both sides; every other byte counts, including the command line, whitespace, key order, every histogram and curve, and the final brace with no newline. No tolerance. Why: fastp's duplicate filter sets one bit in each of two Bloom planes with two separate atomic updates (`duplicate.cpp:154-166`), so two of its 12 workers that meet the same unseen pair can each win one plane and neither counts a duplicate. The reference does not reproduce the number itself: the two runs of the same exome report 0.0902664 and 0.0902662 (HG001), 0.0977926 and 0.0977925 (HG002). nomics counts in FASTQ pair order, which is what fastp does with one worker: two `--thread 1` runs of fastp 1.1.0 on the HG002 exome give the same JSON as each other and as `nomics prep --fastp-threads 1`, byte for byte, rate 0.0977929 (`scripts/fastp_thread1_check.sh`). In all 13 rows that differ, the nomics rate is the higher one, by at most 0.000001, as lost duplicates predict. The command line is not exempt: the gate passes fastp's own. See `docs/spec/qc-fastp-report.md`. |
| `reports/fastp/**/*.fastp.html` | fastp 1.1.0 | **bytes except the two timestamps and the duplication rate**. Wider gate 2026-11-12 (`scripts/qc_reads_gate.sh` on the Spark, the 60 rows of the FastQC row above): 60 of 60 | Exempt: the text between `<title>fastp report at ` and ` </title>`, the text between the footer's `fastp 1.1.0, at ` and ` </div>`, and the number in the summary row `duplication rate:` (the `%` stays). The timestamps are wall clock; the rate is the JSON row's race. Every other byte counts, plot data included. No tolerance. |
| `reports/fastp/**/*.fastp.log` | fastp 1.1.0 | **bytes except the duplication rate and the elapsed seconds**. Wider gate 2026-11-12 (`scripts/qc_reads_gate.sh` on the Spark, the 60 rows of the FastQC row above): 14 of 60 identical with nothing exempted, 46 identical but for those two numbers | Exempt: the number on `Duplication rate: N%` and the integer in `fastp v1.1.0, time used: N seconds`. Every other byte counts, including blank lines, report paths, the command and its trailing space. The file gated is the task's log; what the shim prints to its own stderr is not part of it. No tolerance. |
| `reports/samtools/**/*.{md,recal}.bam.stats` | samtools 1.24 + htslib 1.24 `stats --threads 1 --reference REF INPUT.bam` | **bytes**; all six tier0/wgs20/tier1 md/recal comparisons pass. Wider gate 2026-11-12 (`scripts/qc_bam_gate.sh` on the Spark): every sarek BAM that still has its reports, 130 BAMs in 64 runs, 130 pass (79 distinct report contents) | None. Version, command basenames (overridable), comments, all sections and floating-point formatting are reproduced. The global `--reference` decoder option does not activate stats’ `-r/--ref-seq` mode: for BAM it opens no FASTA, emits no MPC, and estimates GC-depth from reads. Supplementary records, per-bin insert halving/trimming, the coverage ring’s shortened full flush and byte-count resize copies, and GC-depth’s final-bin round/sort/output mismatch follow the pinned source. CRAM is not gated. See `docs/spec/qc-bam-stats.md`; `tests/bam_stats_goldens.rs` invokes the release CLI and literal `cmp`. |
| `reports/markduplicates/*.md.bam.metrics` | Picard `MetricsFile` | **bytes except lines 2 and 4** | Line 2 is `# <command line>` (lists `--INPUT` in task-completion order); line 4 is `# Started on: <wall clock>`. Everything else (metrics row, all histogram rows) must match byte for byte. |
| `*.md.bai`, `*.recal.bai` (Picard naming: `x.bam` -> `x.bai`) | htsjdk `BAMIndexer` (MarkDuplicates `--CREATE_INDEX true`, ApplyBQSR default) | **htsjdk layout, two-sided** | The index holds virtual offsets into its own BAM, so it can only equal the golden byte for byte when the BAM does, and the BAM differs in `@PG` (above). The gate is therefore: (1) `nomics_io::htsjdk_index` applied to the *golden* BAM reproduces the *golden* index byte for byte, and (2) our index is that layout applied to *our* BAM. (1) holds for tier0 `md.bai` (1,620,120 B) and `recal.bai` (1,130,784 B). The layout is htsjdk 4.2.0's, not samtools': chunks coalesce only within the same or the next block address, linear-index gaps take the last earlier value, a metadata pseudo-bin 37450 and the no-coordinate count are written, and the last record of a BAM ends at `(last block, length)` because `BAMIndexer.finish` runs before the BAM is closed. The layout depends on the BAM's block payload sizes, which both compression modes share; only block addresses differ. `samtools index` of the same BAM differs (1,604,512 B). `nomics dev htsjdk-index`. |
| `*.vcf.gz.tbi`, `*.g.vcf.gz.tbi` | htsjdk `TabixIndexCreator` (GATK `--create-output-variant-index`, default true) | **htsjdk layout, two-sided**, as for `.bai` | Same two-sided gate. (1) holds for tier0 `haplotypecaller.vcf.gz.tbi`, tier0 and tier1 `g.vcf.gz.tbi` (548,349 B and 585,150 B, the latter over all contigs). TBI payload BGZF-compressed with htsjdk's level-2 deflater. Differs from `BAMIndexer` in one place: GATK's VCF writer flushes the stream before `finalizeIndex`, so the last feature ends at the next block address. `tabix -p vcf` of the same file differs. |
| `*.artifactprior.tar.gz` (somatic orientation model, S1-r1 tumour sample) | GATK 4.6.2.0 LearnReadOrientationModel | **ordered archive member names + exact member bytes**; both original scatter inputs, chr17 then chr1, pass | Decode tar/gzip envelope only: timestamps, ownership and compressor bytes are incidental and excluded. No prior-table normalization. One `./HCC1395_HCC1395_WES.orientation_priors` leaf, 15,642 bytes / 64 rows. Reversed input order changes floating low bits and is tested separately, never sorted away. Single-output-sample collector layout; unsupported filesystem-dependent multi-sample/member orders fail explicitly. See `docs/spec/mutect2-filtering.md` §1.5. |

## The `reports/` tree of a run

`nomics run --reports` writes the reports above at sarek's paths, from the run's own reads,
BAMs and VCF. `scripts/qc_run_gate.sh GOLDEN_RUN OUTDIR` accounts for every file under the
golden's `results/reports`: it is compared under its row above, or it is one of FastQC's
rendered pages (`*_fastqc.html`, and the images inside the zip), which are not reproduced, or
it is reported missing. A file of ours with no golden counterpart fails the gate too. FastQC's
two text files are written in the directory FastQC makes with `--extract`
(`reports/fastqc/<sample>-<lane>/<sample>-<lane>_N_fastqc/`), not as a zip. The fastp command
line, the staged FASTQ names and the 12-thread insert-size sample are sarek's, written by
nomics; nothing is read from the golden.

Spark, 2026-11-12, each with the tier's own gates passing with the flag on:

| Tier | Report files compared | Out of scope |
|---|---|---|
| tier0 exome | 26 of 26 | 2 FastQC pages |
| tier0nv, two lanes | 31 of 31 | 4 |
| tier1 exome | 26 of 26 | 2 |
| tier1 exome in GVCF mode (`--emit-gvcf`, against a `--joint_germline` run) | 22 of 22: no VCF reports on either side | 2 |
| wgs20 (chr20 genome, streamed stages) | 22 of 22 | 2 |
| wgs30s (30x slice, against an x86 run: `--reference-arch x86_64`) | 22 of 22 | 2 |
| wgs30 (HG002 30x genome, 41 of 41 main gates with the flag on) | 18 of 18 against the x86 run; the 4 VCF reports against the reference tools, see below | 2 |

The 30x golden's four VCF reports were made from its default-x86 (GKL) VCF, which is not the
VCF nomics matches (`docs/REFERENCE-VARIABILITY.md`), so they cannot be compared. Instead
sarek's own vcftools 0.1.17 and bcftools 1.23.1 containers were run on the Java-mode reference
VCF: the reports of the nomics run's VCF are identical to theirs in every byte, all four.

Not covered: CRAM output (`--reports` refuses it: the samtools and mosdepth ports read BAM
only) and `--resume-from` (refused: the read reports come from the FASTQ pass). In GVCF mode
no bcftools or vcftools report is written, as in sarek.

## Somatic pileup gathering and contamination (sarek 3.10.0, GATK 4.6.2.0)

| Boundary (S1-r1 golden) | Written by | Gate | Normalization and evidence |
|---|---|---|---|
| `gps-chr21/{HCC1395_WES,HCC1395BL_WES}.chr21.pileups.table` | GATK 4.6.2.0 GetPileupSummaries from BAM + resource VCF + BED | **bytes; both cmp pass** | None. Tumour 3,084 rows; normal 3,115. All 11 ordered read-filter counts and processed totals match GATK logs. Full recalibrated BAMs, 2026-11-12 (`scripts/somatic_boundary1_gate.sh` on the Spark): the four per-group tables (tumour 50,034 and 194,263 rows, normal 50,384 and 195,763) are identical and the read-filter counts equal GATK's logs; the gathered, contamination, segmentation and orientation-prior files computed from those tables are identical too. See `docs/spec/mutect2-filtering.md` §2.4. |
| `GATHERPILEUPSUMMARIES_TUMOR-HCC1395_WES/HCC1395_WES.mutect2.pileups.table` | GatherPileupSummaries | **bytes; cmp passes** | None. 244,297 rows. Both per-group inputs in recorded argument order; actual GRCh38 sequence dictionary. |
| `GATHERPILEUPSUMMARIES_NORMAL-HCC1395BL_WES/HCC1395BL_WES.mutect2.pileups.table` | GatherPileupSummaries | **bytes; cmp passes** | None. 246,147 rows; same gather contract. |
| `CALCULATECONTAMINATION-HCC1395_WES_vs_HCC1395BL_WES/*.contamination.table` | CalculateContamination, matched normal | **bytes; cmp passes** | None. Entire 86-byte table, including Java-formatted estimate and standard error. |
| `CALCULATECONTAMINATION-HCC1395_WES_vs_HCC1395BL_WES/*.segmentation.table` | CalculateContamination `--tumor-segmentation` | **bytes; cmp passes** | None. Entire 1,855-byte table, 42 rows, inclusive boundaries, Java HashMap contig order and exact MAF doubles. |

`crates/nomics-somatic/tests/contamination_goldens.rs` runs literal `cmp` on all
four files and skips with a message if the golden directory is absent. Both
published reference runs agreed. Synthetic no-normal and empty-input cases have
no golden. Numerical Java oracles and remaining validation limits are in
`docs/spec/mutect2-filtering.md` §3.6. That test starts from the golden per-group
tables; BAM counting is covered by the first row. None of this implies performance
relative to any reference tool.

Beyond the one golden (2026-11-12): CalculateContamination and LearnReadOrientationModel fit a
model to their input, so both are also compared with GATK 4.6.2.0 on inputs the ports were
not developed on (`scripts/somatic_boundary1_cases.sh make`, then `gate`). Contamination:
eight row subsets of the gathered tables, estimates from 0.001 to 0.062 and 2 to 42 segments;
contamination and segmentation tables identical in all eight. Orientation: each scatter
archive alone and both orders; all four identical, and GATK's own two orders differ from
each other, as the port's do.

## RNA-seq alignment quantification (nf-core/rnaseq 3.27.0)

| Boundary | Reference | Gate | Normalization and why |
|---|---|---|---|
| `star_salmon/<sample>/quant.sf`, `quant.genes.sf` | Salmon 2.7.0 deterministic alignment mode | **bytes: 6/6 samples, both files** | None. Five original ISR/SR samples plus the human IU fixture; flags follow each sample's recorded command. A on the human fixture also resolves to IU and matches the IU golden. |
| `lib_format_counts.json`, `aux_info/ambig_info.tsv`, `aux_info/{expected_bias,observed_bias,observed_bias_3p}.gz` | Same | **bytes: 6/6 samples** | None, including compressed bytes of all three bias files. |
| `aux_info/fld.gz` | Same | **decompressed bytes: 6/6 samples** | Decode gzip only. Both files contain the identical little-endian i32 expected fragment-length histogram. The reference uses zlib-rs level 6, while nomics-io's fresh compressor implements zlib-fast levels 1–3 (level 3 here). Compressed differences and first byte offsets are retained in `docs/research/pipelines/nomics-quant.md`; no histogram bins are changed or rounded by the gate. |
| `aux_info/meta_info.json` | Same | **not emitted / not gated** | Timings/RSS vary; explicitly outside this task's requested gate. |

Source ledger, all twelve explicit library types plus A, incompatible priors, mixed paired/single-end fixture commands, and human IU validation: [nomics-quant](research/pipelines/nomics-quant.md).

## Single-cell oracle counting (nf-core/scrnaseq 4.2.0)

| Boundary | Reference | Gate | Normalization and why |
|---|---|---|---|
| `Solo.out/Barcodes.stats`, `Gene/{Features.stats,Summary.csv,UMIperCellSorted.txt}` | STAR 2.7.10b, PBMC 1k v3 | **bytes: all four identical** | None. |
| `Solo.out/Gene/{raw,filtered}/{matrix.mtx,barcodes.tsv,features.tsv}` | Same | **bytes: all six identical** | For nf-core's gzipped published leaves, decompress gzip only: gzip headers contain mtimes. No sorting, rounding, or content masking. |

The counting milestone consumes an oracle BAM with raw CR/CY/UR/UY and lowercase
gx gene evidence. Corrected CB/UB are not counting inputs. The original ordered
features.tsv is an explicitly authorized input until GTF derivation is implemented.
Source ledger, exact input contract, measurements, and test command: [SOLO-COUNT.md](SOLO-COUNT.md).

## Test oracles (verification only)

Some gates compare intermediate state that the reference never writes. For those we run an **oracle**: the reference tool built from its pinned
source, reproducing the golden outputs, optionally instrumented to dump internal state. Oracles are test fixtures in the same category as the
samtools calls in `reference/compare_runs.sh`; nothing in `nomics` links, calls or ships them.

- **bwa 0.7.19 oracle** (`reference/oracle/`, built at `~/omics/oracle/bwa-0.7.19` with the golden binary's flags). Reproduces all tier0 chunk
  BAMs record-identically after the golden's `samtools sort`. Stage dumps (SMEMs, chains, regions, pairing) gate `nomics-align` per read.

## FilterMutectCalls (sarek 3.10.0, GATK 4.6.2.0)

| Boundary (golden) | Written by | Gate | Normalization and evidence |
|---|---|---|---|
| `*.mutect2.filtered.vcf.gz.filteringStats.tsv` | GATK 4.6.2.0 FilterMutectCalls | **bytes** | None. S1-r1 (HCC1395 tumour-normal exomes): all 1,900 bytes. The file prints the learned log priors to 17 digits, so it is the sharpest check on the learning arithmetic. |
| `*.mutect2.filtered.vcf.gz` | GATK 4.6.2.0 FilterMutectCalls (htsjdk) | **records identical; header identical except `##GATKCommandLine`** | The existing rule for GATK VCFs. S1-r1: 6,351 of 6,351 records, 2,024 PASS, every INFO and FORMAT field compared. Compressed bytes and the `.tbi` are not gated; nomics writes valid BGZF and no index yet. |

One golden is one learned model, so the port is also gated on 17 inputs it was not developed
on: subsets of the golden's raw calls (alternate records, thirds, contig ranges, random 70%
samples, the first 500) and altered contamination estimates (0.02, 0.1, 0.3), each filtered
by GATK 4.6.2.0 in sarek's container (`scripts/somatic_filter_cases.sh`), then compared
under the two rows above (`scripts/somatic_filter_gate.sh`): 17 of 17 pass, 2026-11-12. The
first five of these found a real difference: the port matched the golden in every byte, but
on two subsets ten and two of the 21 log priors differed from GATK's in the last digits,
because `pow`, `log10`, `log1p` and `expm1` came from the platform's libm where GATK uses
Java's. GATK itself gave the same bytes on four runs of each. The port now uses the
workspace's Java-exact routines throughout, and the twelve later cases, made after the fix,
passed first time. See `docs/spec/mutect2-filtering.md` §9.3 to §9.5.

Limits: one tumour-normal pair and sarek's default options. Tumour-only and mitochondrial
modes, several tumour samples and custom thresholds are not supported: the command has no
options for them and they have no golden.
GATK multiplies three error terms in the order of a HashMap keyed by enum identity hashes;
all six orders give the golden's bytes, but an order that changed a printed digit can be
constructed, so this is a possible, unobserved source of reference variability
(`docs/REFERENCE-VARIABILITY.md`).

## Mutect2 caller, first stage: reads, activity and assembly regions (GATK 4.6.2.0)

`nomics somatic mutect2-regions`. The goldens are GATK's own Mutect2 on one-contig slices of
the HCC1395 tumour-normal exomes with sarek's arguments (`scripts/somatic_mutect2_slice.sh`);
on chr21 and chr22 the slice run's calls are the full sarek run's calls there (38 and 134
records). All comparisons are exact: no tolerance, sorting or deduplication.

| Boundary (chr21 slice) | Reference | Gate | Evidence |
|---|---|---|---|
| Assembly regions | `m2-chr21/oracle/regions.igv`, GATK 4.6.2.0 | Entire file bytes | 7,167 regions in traversal order, 699 active, including markers, sizes and formatting. |
| Callable sites | `m2-chr21/r1/*.vcf.gz.stats` | Entire file bytes | 35 bytes, `callable\t1230643.0`. |
| Read-filter summary | `m2-chr21/r1/gatk.log` | All twelve count/name lines, their trailing spaces and order, then the filtered/processed totals | Extract only the log prefix on the first line; no count/content normalization. 1,025,272 filtered of 3,023,196 processed. |
| Raw activity | `MutectRegions.java`, unmodified jar evaluator | Entire `dump.activity.tsv` bytes | All 1,493,217 traversed loci, including zeros and gaps. |
| Region read assignment | Same harness, unmodified region iterator | Entire `dump.reads.tsv` bytes | 4,180,613 ordered assignments, including inactive regions. Original SAM flags, sample, name, start and CIGAR; no flag reconstruction. |
| Read payload at next-stage callback | Same harness after the source's `cleanOverlappingReadPairs` | Entire `dump.prepared.bin` bytes | In region/read order: big-endian u32 length, then bases, BQ, BI, BD (length bytes each). No numeric tolerance. This is before assembly or unmarked-duplicate removal. |


The first three rows are gated on five contigs (`scripts/somatic_mutect2_regions_gate.sh`,
2026-11-12): chr21, on which the port was developed (7,167 regions, 699 active), and chr22
(13,815 regions, 1,170 active; 2,689,447 callable) and chr19 (34,338 regions, 3,386 active;
6,845,249 callable), made afterwards, and chr20 (15,433; 1,367 active) and chr13 (11,531;
1,166 active), made after the second stage; all identical. With the final stage they were
run on the held-out chr11 (35,007 regions, 3,219 active) and chr17 (an intact GATK file: see
the final stage), also identical. The last three rows are internal state
dumped from the unmodified GATK jar by `reference/oracle/mutect2/MutectRegions.java` and are
gated on chr21 only (`reference/oracle/mutect2/VERIFICATION.md`); the integration test
`mutect2_regions_goldens.rs` compares them when `NOMICS_M2_ORACLE_PREFIX` is set. Only sarek's
tumour-normal form is supported.

## Mutect2 caller, second stage: assembly and haplotypes (GATK 4.6.2.0)

No normalization or tolerance. The command is
`nomics somatic mutect2-haplotypes`; source/limits are in `docs/spec/mutect2.md` §12.

| Boundary | Oracle | Gate | Scope |
|---|---|---|---|
| Finalized assembly reads/reference | `MutectAssembly.java`, qualified local JDK 17 and unmodified GATK 4.6.2.0 jar | Entire ordered `R`, `REF`, `D` rows by bytes, including BQ/BI/BD | Every assembled chr21 active region. Shared-read edits remain stage-one behavior. |
| Untrimmed and trimmed assembly results | Same jar harness | Entire internal `.assembly.tsv` by bytes; reference state, bases, CIGAR, genomic span, alignment offset, Java score text, k-mer size and list order | chr21; no-variation exits included. |
| Artificial haplotype BAM records | Original sarek Mutect2 `oracle/haplotypes.bam` | Contig, position, CIGAR, bases, `HC:i`, `RG:Z`, `CR:Z`, in file order; no field tolerance | chr21 **3,160**, chr22 **6,883**, chr19 **19,219** records, all exact; and on two contigs made after the port was finished, chr20 **7,802** and chr13 **5,244**, all exact (2026-11-12). `ArtificialHaplotypes.java` rejects any unexpected tag. Biological realigned reads depend on later stages and are not projected. |

`reference/oracle/mutect2/run_assembly_oracle.sh OUTDIR` qualifies the JVM on
chr21 VCF/stats/regions first, then generates the assembly dump and independent
adaptive-pruner fixtures. `NOMICS_M2_ASSEMBLY_ORACLE=.../assembly.assembly.tsv bash
scripts/somatic_mutect2_haplotypes_gate.sh OUTDIR` performs the CLI comparisons.
The integration test is `mutect2_haplotypes_goldens.rs`, with
`NOMICS_M2_ASSEMBLY_GOLDENS` and `NOMICS_M2_REFERENCE` overrides. Missing contig
goldens skip with a message; absent internal oracle skips only that comparison.
The results and remaining limits are recorded in the appended stage-two
verification record. The first-stage gate passes on the same five contigs. The
HaplotypeCaller code this stage shares was parametrized, not changed: the tier0 germline
run still passes 41 of 41 digests on this branch.

## Mutect2 caller, third stage: likelihoods and realigned reads (GATK 4.6.2.0)

Exact gates. They cover the
fixed sarek 3.10.0 tumour-normal/default GATK 4.6.2.0 Java
`LOGLESS_CACHING` path. No numerical tolerance or rounding. The record gate compares a
projection of the BAM records, not BAM bytes: tags are taken by name and integer tags by
value, because the port writes no BAM at this stage. When this stage was accepted the XA tag was
left out, because genotyping sets it; since the final stage the projection includes it. Beyond the five contigs the port was developed on (1,337,133
records), it was run on two made afterwards, 2026-11-12: chr16 (526,456 records) and chr12
(537,627), all identical, with the first two stages passing there too (14,892 and 19,179
haplotype records).

| Boundary | Gate | Evidence / scope |
|---|---|---|
| PairHMM inputs and all pre-genotyping likelihood adjustments | Byte-identical ordered dump: region traversal ordinal, trimmed padded span, stub count; ordered haplotypes; sorted header sample list (including empty samples); reads in matrix order with name/flag/start/CIGAR/bases/BQ/BI/BD/GCP; complete raw log10, capped log10, retained log10 and ln matrices, each value as its raw binary64 hexadecimal bits; realigned read coordinates and optional HC tag | `reference/oracle/mutect2/MutectLikelihoods.java` against unmodified jar, qualified local JDK 17; public-engine consistency assertion; chr21. `NOMICS_M2_LIKELIHOODS_ORACLE` supplies the dump to the integration test and gate script. |
| Realigned biological records, file order | Every non-artificial BAM record: name, flags, reference, 1-based position, MAPQ, CIGAR, mate reference/position, TLEN, bases, base qualities and every tag (XA included since the final stage). Tags are sorted by name, integers compared by SAM integer value rather than BAM storage width; floats remain raw bits. No coordinate or record sorting is applied to the golden. | `scripts/somatic_mutect2_likelihoods_gate.sh OUTDIR`; five contigs chr21/22/19/20/13. AR/CR/HC/BI/BD/PG/RG/AS/XS/MC/MD/MQ/NM/SA included wherever present. PG/RG IDs are translated by the source-defined header merge, not normalized. AR/CR are attached inside `callMutations`, but their span values are already fixed here and included. XA depends on the called alleles; it is compared since the final stage. |
| Shared HaplotypeCaller regression | Existing HC tests and 41/41 tier0 stage digests unchanged | `cargo test --release -p nomics-somatic -p nomics-io -p nomics-hc`; `scripts/mac_tier0.sh cpu OUTDIR`, with `NOMICS_RG_DS` set to the reference path of the golden run. |

`nomics somatic mutect2-likelihoods` produces boundary projections, not a somatic VCF or
BAM; the VCF is the final stage's. Verification
commands, counts, timings and limitations are appended to
`reference/oracle/mutect2/VERIFICATION.md`; no stage-three gate is claimed on a
missing/skipped fixture. The integration test checks all five available
contigs, printing an explicit skip for absent golden directories; override the
root with `NOMICS_M2_LIKELIHOODS_GOLDENS` and reference with `NOMICS_M2_REFERENCE`.

## Mutect2 caller, final stage: raw VCF, stats and f1r2 archive (GATK 4.6.2.0)

`nomics somatic mutect2` takes the tumour and normal BAMs, the reference, the intervals, the
panel of normals and the germline resource, as sarek 3.10.0 runs GATK4_MUTECT2 for a
tumour-normal pair, and writes the raw VCF, `<output>.stats` and, with `--f1r2-tar-gz`, the
f1r2 archive, all from one pass over the reads. Goldens: GATK's own Mutect2 on one-contig
slices of the HCC1395 exome pair (`scripts/somatic_mutect2_slice.sh`).

| Output | Reference | Rule | Normalization |
|---|---|---|---|
| `*.mutect2.vcf.gz` records | GATK 4.6.2.0 Mutect2 | **every record, every byte** of the decompressed text | None |
| `*.mutect2.vcf.gz` header | same | **every line, in order, except `##GATKCommandLine`** | The existing rule for GATK VCFs. Compressed bytes are not gated; no `.tbi` is written yet. |
| `*.mutect2.vcf.gz.stats` | same | **bytes** | None |
| `*.f1r2.tar.gz` member names | same | **the set of three names**, `./` prefix included | Member order: GATK writes the members in the order its temporary directory lists them, which differs between machines (`docs/REFERENCE-VARIABILITY.md`). nomics always writes `ref_histogram`, `alt_histogram`, `alt_table`. The tar and gzip envelopes are not gated. |
| `ref_histogram`, `alt_table` in the archive | same | **bytes** | None |
| `alt_histogram` in the archive | same | **every value under its column label; every other line and the line endings by bytes** | Column order: GATK's 384 data columns come in a hash-map order that changes with the command line (`docs/REFERENCE-VARIABILITY.md`). The label set must be equal; rows are not sorted. |
| Orientation priors learned from the archive | `nomics somatic learn-read-orientation-model` on GATK's archive | **bytes** of the priors archive | None. A check that the archive rule above loses nothing that LearnReadOrientationModel reads. |
| Realigned reads with their `XA` tags | GATK's `--bam-output` | The third stage's record projection, now with every tag | As in the third stage; `XA` is no longer left out. |

`scripts/somatic_mutect2_gate.sh [CONTIG ...]` runs the command once per contig and applies
the first seven rows; a contig without a golden fails. The last row is
`scripts/somatic_mutect2_likelihoods_gate.sh`.

Results, 2026-11-13, from a build of the merged branch on the Spark. The port was developed
against seven contigs; chr17 and chr11 were made beforehand, kept off the machine the port
was written on, and run once it was finished.

| Contig | VCF records | Header, stats | f1r2 archive, priors | Reads with tags |
|---|---|---|---|---|
| chr17 (held out) | 305 of 305 | identical | identical | 534,667 identical |
| chr11 (held out) | 306 of 306 | identical | identical | 496,871 identical |
| chr21 | 38 of 38 | identical | identical | 95,695 identical |
| chr22 | 134 of 134 | identical | identical | 262,031 identical |
| chr19 | 335 of 335 | identical | identical | 563,178 identical |
| chr20 | 144 of 144 | identical | identical | 250,879 identical |
| chr13 | 97 of 97 | identical | identical | 165,350 identical |
| chr16 | 449 of 449 | identical | identical | 526,456 identical |
| chr12 | 337 of 337 | identical | identical | 537,627 identical |

2,145 records in all. The held-out pair passed on its first run, and the earlier stages'
gates pass on all nine from the same build (assembly regions, stats and filter counts;
haplotype records, 16,719 and 15,612 on the held-out pair).
Three runs of chr17 give the same VCF (decompressed bytes). The archive from the combined command and from
`nomics somatic mutect2-f1r2` alone are the same bytes on both held-out contigs. The HaplotypeCaller code the caller shares is unchanged in output: tier0 passes 41 of
41 digests on this branch.

One golden was wrong, not the port: GATK's `--assembly-region-out` file for chr17 was missing
18,881 bytes in the middle. A second GATK run wrote the whole file, which equals the port's
(`docs/REFERENCE-VARIABILITY.md`); the intact file is the golden now.

Limits. One tumour and one named normal with sarek's default arguments; no other mode
is supported, and flow-order read groups are refused. Where GATK's output order comes from the
iteration order of a Java hash map, the port reproduces that order for the cases the source
makes deterministic and refuses the one it cannot: a hash bin that Java would turn into a
tree (nine or more colliding allele keys), whose order can depend on identity hashes. No
golden reaches it; the most alleles in one record is five (`docs/spec/mutect2.md` §16). A
run holds all reads of its intervals in memory: 10.9 GiB (chr11) to 12.5 GiB (chr19) at peak
for one chromosome of this exome pair, so the two scatter groups of a whole exome need a machine with far more than a
laptop's memory until the traversal is streamed. The scatter-and-merge of sarek's two
interval groups has its own commands, below; the whole pair run through them is still to do.

### Merging the scatter groups (GATK 4.6.2.0 MergeVcfs, MergeMutectStats)

sarek runs Mutect2 once per interval group and merges: `nomics somatic merge-vcfs` and
`nomics somatic merge-mutect-stats`.

| Output | Reference | Rule | Normalization |
|---|---|---|---|
| `*.mutect2.vcf.gz` (merged) | GATK 4.6.2.0 MergeVcfs (Picard) | **every decompressed byte**, header included | None: given the inputs in sarek's order even `##GATKCommandLine` is GATK's, because the first input's line is kept. Compressed bytes are not gated; no `.tbi` is written. |
| `*.mutect2.vcf.gz.stats` (merged) | GATK 4.6.2.0 MergeMutectStats | **bytes** | None |

Evidence, 2026-11-13. The golden (test `merge_goldens`): the two groups of the HCC1395 pair,
6,351 records, and `callable` 1.20404485E8. Beyond it, GATK was run on three sets of inputs
the commands were not written against (`scripts/somatic_merge_cases.sh`, on the Spark): the
nine slice VCFs in a scrambled order (2,145 records), two of them, and three with a
records-free VCF among them; VCF and stats identical in all three.

Scope: the commands are exact for Mutect2's scatter outputs and refuse the rest. The inputs'
headers must be equal but for `##GATKCommandLine` (Picard merges differing headers by rules
that are not ported), the sequence dictionary must list the header's contigs, and two inputs
may not have a record at the same position (Picard keeps both; in which order is not
confirmed, and disjoint interval groups cannot produce it). Only the `callable` statistic is
summed.

Internal boundaries used while porting (fragment grouping, the somatic evidence matrices and
LODs, unrounded annotation values, phase sets) have their own exact dumps against the
unmodified GATK jar on chr21, recorded in `reference/oracle/mutect2/VERIFICATION.md`; the
diagnostic command `nomics somatic mutect2-genotyping` writes them.
