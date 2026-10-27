# Bulk RNA-seq: extending nomics to nf-core/rnaseq (research brief, 2026-10-15)

*"The predecessor project" is an earlier, unpublished RNA-seq reimplementation on the same DGX Spark; citations to its files are kept for provenance.*

Research brief. Implementation updates (2026-10-15): [nomics-quant](nomics-quant.md) gates Salmon 2.7.0 alignment-mode transcript/gene tables byte-identically, and `nomics-count` matches featureCounts (20/20 files; `docs/FEATURECOUNTS.md`), both on all five rnaseq-test samples; the other RNA-seq stages below remain research. Tags: **[measured]** = measured by us or by the predecessor project `the predecessor project` on the same Spark (cited by file), **[source]** = read in the tool's or pipeline's code/config, **[claimed]** = vendor/author claim, unverified, **[estimate]** = ours.

## Summary

- **Recommendation: go, conditionally, scoped to the default `star_salmon` path.** The condition is one cheap experiment (Phase 0): run nf-core/rnaseq 3.27.0 twice on the same input and check whether `quant.sf` comes out byte-identical r1 vs r2. If it does, then for the first time this pipeline has a byte-gateable quantification output. If it doesn't, Salmon falls back to a statistical gate, and the "exact" pitch loses its best asset. In that case the answer becomes "not yet".
- **The fact that changes the picture since the predecessor project:** the current release (3.27.0, 2026-09-23) moved from Salmon 1.10.3 (C++) to **Salmon 2.7.0, a from-scratch Rust rewrite that has been deterministic by default since 2.6.0**, and from StringTie 2.2.3 to 3.0.3. Everything the predecessor project measured about Salmon non-determinism was against 1.10.3: the run-to-run TPM envelope, the CAS contention, the proposed statistical-equivalence policy. **None of that applies to the current pin.** The STAR findings still apply, because STAR 2.7.11b is unchanged.
- **STAR is the long pole and cannot be byte-exact against a multithreaded reference.** By construction, its unsorted output order and the primary flag on transcriptome multimappers depend on which thread got which read. The predecessor already measured and documented the exact normalization that works [measured].
- **The QC chain is no longer a differentiator.** Seqera's RustQC already ships inside nf-core/rnaseq as `--use_rustqc`. It replaces 15 QC tools in one pass, claims ">60x", and admits "near-identical" outputs. NVIDIA Parabricks (`--use_parabricks_star`) and Sentieon STAR are already alternative STAR branches. Our differentiators are exactness plus covering alignment and quantification in one resident process.
- **Effort [estimate]:** about 13.5–18 engineer-weeks for FASTQ → published `star_salmon` outputs at the pace nomics has shown. That is roughly 3× the germline effort, mostly STAR (5–8 weeks). The first gateable milestone is about 2 weeks after Phase 0: the post-alignment core on golden STAR BAMs (sort, Picard markdup, Salmon `-a`, featureCounts biotype, samtools stats).

## 1. Pipeline anatomy (nf-core/rnaseq 3.27.0)

The latest release is **3.27.0 "Vanadium Vulture", published 2026-09-23**. Before it: 3.26.0 (2026-05-07) and 3.25.0 (2026-04-24) ([releases API](https://api.github.com/repos/nf-core/rnaseq/releases?per_page=5)). Five releases (3.23.0 to 3.27.0) in seven months, so the pin moves quickly. Versions below come from each module's `environment.yml` at tag 3.27.0 ([tree](https://github.com/nf-core/rnaseq/tree/3.27.0/modules/nf-core)) [source].

**Default path** (`aligner = 'star_salmon'`, `trimmer = 'trimgalore'`, `skip_markduplicates = false`, `skip_stringtie = false`, `skip_qc = false`, `skip_preseq = true`, `use_rustqc = false`; [nextflow.config](https://github.com/nf-core/rnaseq/blob/3.27.0/nextflow.config)):

| Step | Tool, pinned version | Notes |
|---|---|---|
| FASTQ lint / QC | fq 0.x lint, FastQC 0.12.1 | |
| Trimming (default) | **Trim Galore 2.3.0** (Rust rewrite; the 2.x line claims byte compatibility with 0.6.11 on CI-covered paths, [Zenodo](https://zenodo.org/records/7579519/latest), [MultiQC docs](https://docs.seqera.io/multiqc/modules/trim_galore)) | `--trimmer fastp` selects **fastp 1.3.6** instead |
| Optional filters | SortMeRNA 4.3.7, BBSplit (bbmap 39.18), ribodetector, Kraken2/Bracken, sylph | off by default (`skip_bbsplit = true`, etc.) |
| Strandedness (when `auto`) | fq subsample + Salmon | |
| Alignment | **STAR 2.7.11b** (+ samtools 1.21 in the module) | two-pass `--twopassMode Basic`, `--quantMode TranscriptomeSAM`, `--outSAMtype BAM Unsorted`, `--outSAMattributes NH HI AS NM MD`, `--sjdbGTFfile`. The predecessor's argv template comes from a real run [source: `predecessor:bench/PHASE2E-production-gate.md`]; the 3.27.0 config is [align_star.config](https://github.com/nf-core/rnaseq/blob/3.27.0/conf/modules/align_star.config) |
| Sort / index / stats | samtools 1.24: sort, index, stats, flagstat, idxstats | |
| Duplicate marking | **Picard 3.5.0 MarkDuplicates** `--ASSUME_SORTED true --REMOVE_DUPLICATES false` ([markduplicates.config](https://github.com/nf-core/rnaseq/blob/3.27.0/conf/modules/markduplicates.config)) | the published genome BAM is this output |
| Quantification | **Salmon 2.7.0** alignment mode on STAR's transcriptome BAM, `--libType` from samplesheet strandedness, no `--gcBias`/`--seqBias` by default ([quantify_pseudo_alignment.config](https://github.com/nf-core/rnaseq/blob/3.27.0/conf/modules/quantify_pseudo_alignment.config)) | |
| Merge / R objects | tximeta 1.20.1 / tximport, SummarizedExperiment, DESeq2 QC (R) | cohort-level tables and `.rds` |
| Transcript assembly | **StringTie 3.0.3** (`-e` with `--fr/--rf`) | |
| Coverage | bedtools 2.31.1 genomecov → UCSC bedClip/bedGraphToBigWig | |
| Biotype QC | Subread 2.1.1 featureCounts `-B -C -g gene_biotype` ([featurecounts.config](https://github.com/nf-core/rnaseq/blob/3.27.0/conf/modules/featurecounts.config)) | not the count matrix; counts come from Salmon |
| QC | RSeQC 5.0.4 (bam_stat, inner_distance, infer_experiment, junction_annotation, junction_saturation, read_distribution, read_duplication; TIN optional), Qualimap 2.3 rnaseq, dupRadar 1.38.0, preseq 3.2.0 (off) | `--use_rustqc` swaps these for RustQC 0.2.1 |
| Report | MultiQC 1.35 | |

**Other branches** (`aligner` enum: `star_salmon`, `star_rsem`, `hisat2`, `bowtie2_salmon`; `pseudo_aligner`: `salmon`, `kallisto`) [source: [nextflow_schema.json](https://github.com/nf-core/rnaseq/blob/3.27.0/nextflow_schema.json)]: RSEM 1.3.3, HISAT2 2.2.3, Bowtie2, kallisto 0.52.0. UMI: umi_tools 1.1.6 extract/dedup, or UMICollapse 1.1.0. Alternative STAR implementations: Parabricks `rna_fq2bam` 4.7.1 and Sentieon 202503.02 STAR/RSEM.

**Gate boundaries.** These are the published leaves under `results/star_salmon/`, plus `fastqc/`, `trimgalore/` and `multiqc/`. The predecessor's two-run screen of a 5M-pair input under 3.26-era pins counted **414 published files per run**, 102 non-PDF leaves under `star_salmon/` alone [measured: `predecessor:bench/PHASE3B-reference-variability.md`]. Ranked by scientific weight:

1. Salmon `quant.sf` / `quant.genes.sf` per sample, plus merged `salmon.merged.gene_{counts,tpm,…}.tsv` and transcript tables (the product most users consume).
2. The markdup genome BAM + `.bai` (input to everything downstream, and to rnavar-like use).
3. STAR `SJ.out.tab`, `Log.final.out`.
4. featureCounts biotype counts, samtools stats/flagstat/idxstats, Picard metrics.
5. StringTie GTF/ctab, bigWigs.
6. RSeQC / Qualimap / dupRadar / preseq text, then MultiQC data.

R objects (`.rds`), PDFs and HTML belong to a later report tier, as in the predecessor.

## 2. Where the time goes

The best data we have is the predecessor's run on this exact box: nf-core/rnaseq **3.26.0** stock, 6 Geuvadis LCLs at 78M pairs each, Spark GB10 [measured: `predecessor:bench/PHASE2-tier2a-throughput.md`, `bench/COST-CURVE.md`]:

| Bucket | CPU-min / sample | Notes |
|---|---:|---|
| Whole pipeline (39 processes) | **223** | 105 min serial wall for one sample; 34.7 min amortized when 6 overlap; only 32% core utilisation |
| STAR | ~65 | 5.2× parallel efficiency; ~40 GB RSS |
| Salmon 1.10.3 | ~51 | 6.3× parallel efficiency; **obsolete: 3.27.0 runs Salmon 2.7.0, unmeasured** |
| Post-alignment BAM/QC chain (Qualimap, Picard, samtools ×5, dupRadar, RSeQC ×7, bedtools, featureCounts) | 88 | 66.6 min task-wall; almost all single-threaded; 58% of task-wall |
| Trimming + FastQC | ~19 | |

Task-wall leaders over 6 samples: Qualimap 128 min, STAR 75, Picard 64, Salmon 49, dupRadar 45, RSeQC read_duplication 28 [measured: same file]. Independent corroboration: Seqera reports the 15 QC tools taking **>15 h sequential wall** on a ~186M-read GRCh38 sample, against ~15 min for RustQC [claimed: [Seqera blog, 2026-04-02](https://seqera.io/blog/rustqc/)].

The predecessor already showed what a resident CPU pass does to the QC chain: `umbam` took it from 66.6 → **2.8 min wall** (CPU) and **1.7 min** (GPU), byte-identical except for documented residuals [measured: `predecessor:bench/COST-CURVE.md`]. After that, STAR + Salmon are ~90% of CPU.

**GPU vs CPU suitability:**

| Stage | Shape | GPU verdict |
|---|---|---|
| STAR seed search (suffix-array binary search, ~30 GB uncompressed SA) | latency-bound random access, batchable | **Plausible but unproven end to end.** Kernel 6.6× on 1M real STAR requests, results identical [measured: `predecessor:bench/PHASE2C-real-requests.md`]. Integrated into stock STAR, it saved only −5.6% CPU vs the huge-page-advised CPU arm, under the 8% target [measured: the predecessor project commit `c8e46b1`]. **THP alone gave −9.4% STAR CPU** under the production two-pass argv [measured: `PHASE2E-production-gate.md` §3]. Seed search is ~39% of STAR CPU [measured: `docs/review-seed-triage.md`]. A ground-up STAR could batch seeds the way `nomics-align` batches BWA seeding; the prize is a fraction of 65 CPU-min. |
| STAR stitching / scoring / two-pass junction insertion | branchy, per-read, plus a mid-run SA rebuild | CPU |
| Coordinate sort, markdup | sort-shaped | GPU 3.7–12× in umbam stages [measured: COST-CURVE fig 4], but small in absolute terms |
| BGZF deflate | | CPU. nvCOMP was 7.5× slower at equal ratio [measured: COST-CURVE] |
| Salmon `-a` | BAM parse, fragment probabilities, VBEM over equivalence classes | CPU. Parse/IO-bound; the EM is small. Deterministic reductions matter more than raw speed. |
| QC chain | one resident pass | CPU-sufficient (2.8 min); GPU bought 1.1 min |
| StringTie | per-bundle graph and flow | CPU |

Net: like germline, the GPU's natural job is index lookups (seeding). The rest is a CPU engineering problem, solved by doing everything in one resident pass.

## 3. Exactness feasibility, tool by tool

| Tool (pin) | Determinism of the reference | Byte-identical reimplementation? |
|---|---|---|
| **STAR 2.7.11b** | **Not byte-deterministic at 16 threads** [measured: PHASE2E §1]. `Aligned.out.bam` differs between two stock runs, but the record multiset is identical after header drop + whole-line sort; `SJ.out.tab` (both passes) is byte-identical. In `Aligned.toTranscriptome.out.bam` the primary flag (0x100) and `HI:i` move between equivalent multimapper records. Cause [source]: `ReadAlign.cpp:11` seeds `rngMultOrder.seed(P.runRNGseed*(iChunk+1))` per thread, and `ReadAlign_quantTranscriptome.cpp:69` picks the primary with `alignT[int(rngUniformReal0to1(rngMultOrder)*nAlignT)]` ([STAR 2.7.11b source](https://github.com/alexdobin/STAR/tree/2.7.11b/source)). Read-to-thread assignment is dynamic. NVIDIA independently documents the same issue and patches the baseline with `rngMultOrder.seed(777)` to match Parabricks ([rna_fq2bam docs](https://docs.nvidia.com/clara/parabricks/tool-reference/tools/rna_fq2bam.md)). Alignment content is identical across aarch64 and x86 [measured: `PHASE2F-x86-transfer.md`]. | **Yes, against two references:** (a) the normalized gate vs the 16-thread golden (the predecessor's comparator, `run_production_gate.py::normalize_bam`); (b) byte-identical vs a `--runThreadN 1` golden, where one RNG stream is consumed in input order. nf-core passes `--runRNGseed 0` ([align_star.config](https://github.com/nf-core/rnaseq/blob/3.27.0/conf/modules/align_star.config)), so every thread seeds `mt19937(0)`. The 16-thread variance therefore comes from how many draws a thread made before a given read, not from different seeds. (777 is STAR's default and the Parabricks patch value; it is not the reference's seed.) Under the default `--outMultimapperOrder Old_2.4`, the shuffles in `ReadAlign_multMapSelect.cpp:70-79` are skipped, so the transcriptome primary pick is the only draw that matters [source]. To hit (b) with a parallel aligner, align in parallel and make the RNG draws (primary choice) in a serial read-order pass at output time [estimate]. Hard part: two-pass `sjdbInsertJunctions` rebuilds the SA and SAindex mid-run (`STAR.cpp:148-149`, `twoPassRunPass1.cpp:92`, `sjdbBuildIndex.cpp:289-290` [source, via the predecessor project commit `496cdb5`]), and it must be exact. |
| **Picard 3.5.0 MarkDuplicates** | Deterministic given input order. **Run-to-run dup-flag swaps in the reference** (50 gained / 50 lost of 13.4M records) are inherited from STAR's order [measured: `predecessor:docs/E2E-BIO-VARIABILITY.md`]. Picard ties follow input order; after `samtools sort`, input order among equal keys traces back to STAR's thread order. | Yes, against a golden with fixed upstream order (1-thread STAR golden, or the golden BAM fed in). `nomics-markdup` targets GATK 4.6.2.0's bundled Picard; a source diff against Picard 3.5.0 is needed (rnavar pins Picard 3.4.0). |
| **Salmon 2.7.0** (`-a`) | **[claimed]** Since 2.6.0, `quant.sf` is "byte-identical across runs and thread counts — in every input mode", including `-a`. Alignment mode scores by the BAM `AS` tag by default, removing 1.x's sampled online error model ([v2.6.0 notes](https://github.com/COMBINE-lab/salmon/releases/tag/v2.6.0)). 2.7.0 output is byte-identical to 2.6.0 on the authors' validation sets ([v2.7.0 notes](https://github.com/COMBINE-lab/salmon/releases/tag/v2.7.0)). **Unknown:** whether output depends on BAM record order or on 0x100/HI, which is exactly what STAR varies. | **Yes, if Phase 0 confirms.** The spec is readable Rust (BSD-3, [repo](https://github.com/COMBINE-lab/salmon)). The EM must reproduce the same reduction order. Note the trap: 2.8.0 already exists (2026-09-25); the pin is 2.7.0. |
| Salmon 1.10.3 (pre-3.27 pins) | **Not deterministic, even single-threaded**: workers seed from `random_device` for error-model sampling [source: `predecessor:docs/review-salmon-next.md`]. r1 vs r2: 12–25% of TPM rows changed [measured: PHASE3B]. | No. Out of scope; superseded by the 3.27.0 pin. |
| featureCounts 2.1.1 | Deterministic | Yes. umbam matched `-p --countReadPairs` CORE exactly [measured: `predecessor:crates/umbam/COMPAT.md`]. nf-core's own snapshots treat it as stable ([Seqera blog](https://seqera.io/blog/rustqc/)). |
| samtools sort/stats/flagstat/idxstats 1.24 | Deterministic given input | Yes. umbam did it; `.bai` is already a planned nomics deliverable (`docs/GATES.md`). |
| RSeQC 5.0.4 | Mostly deterministic. `junction_saturation` uses an unseeded `random.shuffle` (`SAM.py:3973`); junction files come out in encounter order [source: E2E-BIO-VARIABILITY]. | Yes, except junction_saturation's 5–95% points, which are gated at endpoints only, as in umbam [measured: COMPAT.md]. |
| Qualimap 2.3 | nf-core excludes Qualimap's MultiQC files from snapshots as non-deterministic across x86/ARM ([CHANGELOG PR #1770](https://github.com/nf-core/rnaseq/blob/3.27.0/CHANGELOG.md)). umbam did not reach parity on the gene/ambiguous split [measured: COMPAT.md]. | Partial. Plan on documented residuals. |
| dupRadar 1.38.0 | Deterministic | Yes, except a known ±1 on 17/78,900 genes in the multimapper columns (Rsubread tie-break) [measured: `PHASE2A-tier1-validation.md`]. |
| StringTie 3.0.3 | **Unverified.** Under 2.2.3 the predecessor saw numeric differences and 26 extra zero-expression records r1 vs r2, but its input BAM also differed between runs, so the cause was not isolated [measured: E2E-BIO-VARIABILITY]. | Unknown until Phase 0. |
| Trim Galore 2.3.0 / fastp 1.3.6 | Unverified for these versions | Likely. Trimming is per-read. `nomics-prep` currently implements fastp 1.1.0 filtering without adapter trimming. |
| RSEM 1.3.3, HISAT2 2.2.3, umi_tools 1.1.6 | **Unverified.** HISAT2 runs without `--reorder` ([hisat2/align main.nf](https://github.com/nf-core/rnaseq/blob/3.27.0/modules/nf-core/hisat2/align/main.nf)), so multithreaded output order is unordered. | Not planned. Non-default enum branches. |
| tximport / DESeq2 / SummarizedExperiment (R) | Deterministic given inputs. nf-core 3.27.0 now sorts the inputs to keep cache keys stable (PR #1885). | Merged TSVs: feasible (R `write.table` formatting must be matched). `.rds`: leave to stock R. |

**Gate policy proposal.** For each stage: byte gate against a fixed-input golden (feed the golden upstream output, the convention in docs/METHODOLOGY.md), plus a normalized gate against the stock 16-thread golden. Generate a second golden family with STAR `--runThreadN 1` (and StringTie `-p 1` if Phase 0 shows thread sensitivity). Then end-to-end FASTQ → `quant.sf` can be byte-gated, while every deviation from "the reference as run" is written down in `GATES.md`. Nothing is masked except the measured 0x100/HI normalization on the transcriptome BAM.

## 4. Reuse from nomics vs new work

**Reusable as is or with small changes:**

| Component | Use in RNA-seq |
|---|---|
| `nomics-io`: FASTQ, BGZF, own inflate, zlib-exact deflate, BAM read/write | all BAM outputs; needs a `.bai`/`.csi` writer matching samtools' layout (already a planned germline gap) |
| `nomics-mem` (2 MiB huge-page buffers) | STAR SA/genome resident in THP. The predecessor measured −9.4% STAR CPU from THP alone [measured] |
| `nomics-gate` (digests, normalizations) | add the STAR normalization (header drop, whole-line sort, transcriptome HI strip / 0x100 mask) |
| `nomics-markdup` (Picard semantics) | Picard 3.5.0 diff; spliced reads (`N` in CIGAR) in unclipped-5′ math |
| `nomics-gpu` batching pattern + seeding kernels | the *pattern* transfers (batched random lookups in a resident index). The BWA FM-index kernels do not. |
| `nomics-prep` | only if goldens pin `--trimmer fastp`; still needs fastp 1.3.6 adapter trimming |
| `nomics-bqsr`, `nomics-hc` (GATK 4.6.2.0) | **rnavar only**: it pins GATK 4.6.2.0 and STAR 2.7.11b, the same as our germline pins ([rnavar modules](https://github.com/nf-core/rnavar/tree/1.3.0/modules/nf-core)) |

**Does not transfer:** BWA SMEM/chaining/ksw extension and mate rescue. STAR uses maximum-mappable-prefix search in an uncompressed suffix array plus its own seed clustering and stitching ([Dobin et al. 2013](https://doi.org/10.1093/bioinformatics/bts635)). The predecessor project supplies lessons and goldens, not crates.

**New code [estimate, Rust lines; the reference is `nomics-align` at 13.2k for BWA]:**

| Component | Size | Notes |
|---|---:|---|
| STAR aligner: index load, SA/SAindex search, seed windows/stitching/scoring, PE logic, multimapper ordering, SAM/BAM emit, `SJ.out.tab`, `Log.final.out` | 15–22k | STAR's own source is ~1.2 MB of C++ ([languages](https://api.github.com/repos/alexdobin/STAR/languages)). Load STAR's `genomeGenerate` index rather than building it. |
| STAR two-pass `sjdbInsertJunctions` (SA insertion, SAindex rebuild) | 2–3k | runs twice per sample in production; must be exact |
| STAR `quantTranscriptome` (genome → transcriptome projection) | 1–2k | |
| Salmon 2.7.0 `-a` path: fragment-length model, AS-based scoring, eq-classes, VBEM, `quant.sf`, `quant.genes.sf`, `meta_info.json`, `lib_format_counts.json` | 4–7k | readable Rust spec |
| Trimming (Trim Galore 2.3.0 or fastp 1.3.6 adapters) | 1.5–3k | |
| featureCounts biotype + samtools stats/flagstat/idxstats + `.bai` | 2–3k | |
| RSeQC ×7, dupRadar, Qualimap, preseq | 6–9k | umbam (7.2k lines incl. GPU) is the size reference |
| bedtools genomecov + bigWig writer | 1–1.5k | |
| StringTie 3.0.3 `-e` | 6–10k | `-e` mode (quantify reference transcripts only) is smaller than full assembly |
| tximport-equivalent merged TSVs | 0.5–1k | optional |

## 5. Test data and goldens

| Dataset | What | Why | Access / licence |
|---|---|---|---|
| **Geuvadis LCLs** (E-GEUV-1; ERR188140, ERR188217/NA11832_F, NA12399…) | 1000 Genomes LCL PE RNA-seq; the predecessor used 6 × 78M pairs, plus 20M and 5M-pair slices and a chr22 tier0 fixture [measured: the predecessor project bench] | Already on the Spark/NAS with 3.26-era goldens, STAR/x86 parity evidence and timings. Salmon 2.0 was itself validated on Geuvadis ERR188044 ([v2.0.0 notes](https://github.com/COMBINE-lab/salmon/releases/tag/v2.0.0)). | ENA/ArrayExpress, open access. **New 3.27.0 goldens are required.** |
| **nf-core `test_full`** | 8 ENCODE samples (GM12878, K562, MCF7, H1 × 2 reps, reverse-stranded) on `s3://ngi-igenomes/test-data/rnaseq/`, GRCh37 ([samplesheet](https://raw.githubusercontent.com/nf-core/test-datasets/626c8fab639062eade4b10747e919341cbf9b41a/samplesheet/v3.10/samplesheet_full.csv), [test_full.config](https://github.com/nf-core/rnaseq/blob/3.27.0/conf/test_full.config)) | the pipeline's own full-size test; multiple cell types | ENCODE data is unrestricted; sizes not verified here |
| **SEQC/MAQC-III** | UHRR (A), HBRR (B), 3:1 and 1:3 mixes (C, D) with ERCC spike-ins; GEO GSE47774 / BioProject PRJNA208369 ([Sci Data 2014](https://www.nature.com/articles/sdata201420), [BioProject](https://ncbi.nlm.nih.gov/bioproject/208369)) | known mixing ratios and ERCCs: a sanity check on the reference (like GIAB for germline), not the gate | public SRA; run sizes not verified here |
| **GIAB HG002 RNA** | draft Illumina + long-read RNA-seq for HG002 LCL/iPSC, HG004/HG005 LCLs, in `s3://giab` ([NIST GIAB](https://www.nist.gov/programs-projects/genome-bottle)) | pairs with nomics' HG002 DNA goldens; the natural rnavar test (GIAB truth VCF for sanity) | public; NIST marks it draft |
| nf-core test-datasets (rnaseq branch) | tiny yeast/human subsets | CI tier | MIT |

**Fixture ladder:** tier0 chr22 subset (the predecessor's fixture recipe) → tier1 5M pairs (ERR188140) → tier2 full depth (78M pairs, plus one ENCODE `test_full` sample on GRCh38). Phase 0 per tier: two stock 16-thread runs, plus one STAR/StringTie single-thread run.

## 6. Demand

- nf-core/rnaseq has **1,389 GitHub stars and 904 forks**, vs **614 / 561 for sarek** ([rnaseq](https://api.github.com/repos/nf-core/rnaseq), [sarek](https://api.github.com/repos/nf-core/sarek); fetched 2026-10-15). That is more than twice sarek's stars. Download/launch counts were not available (the nf-core stats page did not load).
- Citations (OpenAlex, 2026-10-15): STAR **51,241**; featureCounts 27,019; RSEM 22,074; HISAT2 14,881; StringTie 14,306; Salmon 12,576; kallisto 9,690; nf-core framework paper 4,053 ([OpenAlex API](https://api.openalex.org/works/doi:10.1093/bioinformatics/bts635)).
- Who would value a faster *exact* version: core facilities and CROs running hundreds of samples a month (the predecessor's cost model put one Spark at ~41 samples/day stock [measured]); pharma and regulated groups whose validated SOP names a pinned nf-core version (an exact drop-in needs no revalidation); and consortium reprocessing efforts.
- Caveat on the "exact" pitch: RNA-seq users tolerate version drift. nf-core's own 3.27.0 notes say the Salmon/StringTie bump causes "double-digit-percent TPM/FPKM swings" for low-to-moderate genes, and it shipped anyway ([CHANGELOG](https://github.com/nf-core/rnaseq/blob/3.27.0/CHANGELOG.md)). So exactness sells as "certified identical to 3.27.0, ~5× faster", not as a scientific necessity the way it is for clinical variants.

## 7. Competition

| Product | Scope | Claims | Exact vs stock? |
|---|---|---|---|
| **RustQC 0.2.1** (Seqera, GPL-3; `--use_rustqc`, experimental, off by default) | 15 QC tools in one BAM pass: dupRadar, featureCounts, 8 RSeQC, preseq, Qualimap, samtools ×3 | ">60x faster", >15 h → ~15 min on a 186M-read sample [claimed: [blog](https://seqera.io/blog/rustqc/), [docs.rs](https://docs.rs/rustqc)] | "functionally identical" / "near-identical"; byte-identical only for deterministic tools like featureCounts [claimed] |
| **Parabricks `rna_fq2bam` 4.7.1** (`--use_parabricks_star`) | STAR alignment, sorted BAM, optional markdup | NVIDIA publishes only absolute times: ~5.7–10 min per melanoma sample on 2–4 L4/RTX PRO GPUs, no CPU baseline ([overview](https://docs.nvidia.com/clara/parabricks/latest/overview.html)). One user reports 8–9× (45 → 5–7 min on 4 GPUs), another reports no speedup ([forum](https://forums.developer.nvidia.com/t/clara-parabricks-4-0-0-1-and-4-2-0-1-star-equivalent-parameters-no-speed-up-observed-in-rna-fq2bam-compared-to-equivalent-star/277323)). | Compatible with **STAR 2.7.2a**, not 2.7.11b. Deterministic transcriptome primary differs from stock unless the baseline is patched ([docs](https://docs.nvidia.com/clara/parabricks/tool-reference/tools/rna_fq2bam.md)). GB10 is supported from 4.7 ([4.7 notes](https://docs.nvidia.com/clara/parabricks/latest/whatsnew/WhatsNew4.7.0-1.html)), but a user reports fq2bam exhausting a single GB10's 128 GB ([forum](https://forums.developer.nvidia.com/t/can-i-use-clara-parabricks-on-dgx-spark/341589)). A junction-rich-sample crash has been reported on 4.7.0 ([forum](https://forums.developer.nvidia.com/t/rna-fq2bam-gpu-star-illegal-memory-access-at-1-000-000-collapsed-junctions-independent-of-max-collapsed-junctions/375099)). |
| **Sentieon STAR/RSEM 202503.02** (licensed) | STAR- and RSEM-compatible CPU implementation | none found | nf-core skips Sentieon tests on conda because "license-server-driven output drifts" (PR #1819) |
| Upstream Rust rewrites (Salmon 2.x, Trim Galore 2.x) | the tools themselves | faster index build, lower RSS (Salmon 2.7.0) | they *are* the reference, which makes the reference more deterministic: good for us |
| kallisto / Salmon selective alignment (`pseudo_aligner`) | skip alignment entirely | minutes per sample on CPU | a different output set; not a competitor to the aligned path, but the cheap option users choose when they don't need BAMs |

No GPU spliced aligner other than Parabricks turned up in a search for 2025–26 work. Our angle: the only path that is identical to the stock pinned branch, single-process, on one desk-side box.

## 8. Recommendation

**Go, conditional on Phase 0. Scope: the `star_salmon` default path, GRCh38 + GENCODE, paired-end, no UMI. Out of scope until asked: `star_rsem`, `hisat2`, `bowtie2_salmon`, kallisto, UMI dedup, contaminant screening, R `.rds`/PDF/HTML.**

| Phase | Deliverable / gate | Effort [estimate] |
|---|---|---|
| **0. Reference characterisation** | 3.27.0 goldens on tier0/tier1: two stock runs plus one with STAR `--runThreadN 1` (and StringTie `-p 1`). Classify every leaf r1 vs r2 as in PHASE3B. **Decision point: is `quant.sf` byte-identical r1 vs r2? Probe order/flag sensitivity by feeding fixed Salmon 2.7.0 the r1 and r2 stock transcriptome BAMs, which already differ in order and 0x100/HI. Do not whole-line sort: `-a` needs each read's records to stay contiguous, so a sort would change results for the wrong reason.** Pick `--trimmer`. | 1 wk |
| **1. First gateable milestone: post-alignment core on golden STAR BAMs** | `samtools sort` → Picard 3.5.0 markdup BAM byte-identical; Salmon 2.7.0 `-a` `quant.sf`/`quant.genes.sf` byte-identical from the golden transcriptome BAM; featureCounts biotype; samtools stats/flagstat/idxstats; `.bai` | 2 wk |
| 2. STAR 2.7.11b | normalized gate vs the 16-thread golden; byte gate vs the 1-thread golden; `SJ.out.tab` both passes byte-exact; then GPU seed batching vs the THP CPU arm | 5–8 wk |
| 3. Trimming + FASTQ QC | Trim Galore 2.3.0 (or fastp 1.3.6) outputs byte-identical | 1 wk |
| 4. QC chain + coverage | RSeQC ×7, dupRadar, Qualimap (documented residuals), bigWig, MultiQC-ingestible files | 2–3 wk |
| 5. StringTie 3.0.3 `-e` | byte or documented | 1.5–2 wk |
| 6. Fusion: one process FASTQ → published outputs | end-to-end tier2, timings vs our own saturated CPU path | 1 wk |
| **Total** | | **~13.5–18 wk** |
| Optional: rnavar 1.3.0 | add GATK SplitNCigarReads; reuse BQSR/HC (same GATK 4.6.2.0 pins) and the Phase 2 STAR | +2–3 wk |
| Not recommended now: rnafusion 4.1.3 | Arriba, STAR-Fusion, FusionCatcher; chimeric STAR mode; small user base (176 stars) | — |

**Speed target [estimate]:** today stock takes 105 min serial for one 78M-pair sample [measured, 3.26]. With STAR in a resident, huge-paged, seed-batched implementation, and Salmon/QC/StringTie in the same process, ~15–25 min per sample is plausible, roughly where germline WES landed (68 → 9.8 min). The THP-only STAR result (−9.4%) and GPU-seed result (−5.6%) both say: measure the CPU arm first.

**Top risks:**
1. **Salmon 2.x order sensitivity (decides go vs not-yet).** If `quant.sf` depends on transcriptome-BAM order or on 0x100/HI, the 16-thread reference is not reproducible at the main product, and only the 1-thread golden family is gateable.
2. **STAR size and two-pass exactness.** At 5–8 weeks it is the long pole. `sjdbInsertJunctions` rebuilds the SA mid-run and must match bit for bit. STAR keeps a ~30–40 GB resident index (40 GB RSS measured) alongside our buffers in 121 GB.
3. **Moving reference.** 3.23 → 3.27 took seven months. Salmon 2.x ships about weekly (2.8.0 is already out). Pin 3.27.0 and treat bumps as new projects.
4. **Commoditised QC.** RustQC already serves the QC-speed buyer. Spend QC effort on exactness, not on headline speedups.
5. **GPU may not pay for STAR.** The predecessor's integrated GPU seed search missed its 8% target. Keep the GPU confined to seeding, and only ship it if it beats the THP CPU arm.
6. **Licensing hygiene.** STAR (MIT) and Salmon (BSD-3) are permissive. RSeQC, Qualimap, featureCounts and RustQC are GPL. Reading GPL source as the spec for a fresh reimplementation is how nomics already works with bwa, but it should be a conscious choice.
