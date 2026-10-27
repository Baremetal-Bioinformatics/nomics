# Reference pipeline run-to-run variability

Gates are defined from measured r1-vs-r2 differences, never assumed (a lesson from an earlier, unpublished RNA-seq project on the same machine, where STAR, Salmon and Picard all varied run to run).

## Smoke test (sarek `test` profile, chr22 40 kb), runs 02 vs 03, 2026-10-08

| Boundary | Result | Cause |
|---|---|---|
| HaplotypeCaller VCF | Normalized header identical; **body has 0 records in both runs**, so this boundary says nothing | The 40 kb test region yields no calls |
| MarkDuplicates CRAM | QNAME/FLAG/RNAME/POS identical row-by-row; other record fields and `@SQ UR:` differ | (a) `UR:` embeds the Nextflow work-dir path, which is path noise. (b) The test data uses one FASTQ pair as two "lanes" (L1, L2), so every read exists twice with the same QNAME. Which lane's copy sorts first and is kept by MarkDuplicates depends on the order sarek merges lane BAMs, i.e. Nextflow task completion order. RG tags therefore differ. |
| BQSR recal table | 3,858 lines differ: per-read-group rows for L1/L2 swap values between runs | Same cause as (b): the per-RG covariates follow the arbitrary lane/duplicate assignment |
| Recalibrated CRAM | records differ | Follows from the table |

**Implications**
- Normalize `@SQ UR:` (and `@PG`) in every BAM/CRAM gate.
- Sarek's merge of per-chunk/per-lane BAMs is ordered by task completion. Where records tie on the sort key, their relative order (and duplicate choice) can vary run to run. The tier1 r1/r2 diff will show whether this matters with unique QNAMEs. The 12 fastp chunks in tier1 are merged the same way.
- The smoke data's duplicated lanes exaggerate the effect. It is not representative of real data.
- Correction (2026-10-08): the first version of `compare_runs.sh` ran samtools with relative paths that did not exist inside the container. When both sides failed, both empty outputs hashed the same and reported PASS. The script now uses absolute paths, aborts on any tool error, and reports record counts so empty matches are visible.

## Tier1 (HG002 WES, SRR14724482) and tier0, r1 vs r2 (2026-10-08)

Measured with the corrected `compare_runs.sh`. Same result at both scales:

| Boundary | tier1 | tier0 | Notes |
|---|---|---|---|
| `mapped/*.NNNN.bam` (12 BWA chunks) | bytes identical | bytes identical | fastp chunking and BWA (`-K`) are deterministic |
| `mapped/*.sorted.bam` | **differs** | **differs** | `samtools merge` adds random hex suffixes to clashing `@RG`/`@PG` IDs, so every record's RG tag differs. Rows also differ in tie order (tier1 19,418,742; tier0 228,797). Side output only. |
| `markduplicates/*.md.bam` | records + header (excl. `@PG`) identical | same | Bytes differ by 1-2: `@PG` lines are ordered by task completion, which changes header length and shifts all blocks |
| `recal_table/*.recal.table` | **bytes identical** | **bytes identical** | |
| `recalibrated/*.recal.bam` | records + header (excl. `@PG`) identical | same | as md.bam |
| `haplotypecaller/*.vcf.gz` | 42,358 records identical; normalized header identical | 1,047 records identical | bytes differ (command line/date header lines) |

**Conclusion:** the reference is reproducible in *content* at every boundary we gate, with no residual run-to-run variability. The only nondeterminism is header metadata (`@PG` order, merge ID suffixes, VCF command/date lines). "Perfect concordance" is therefore achievable as a target: any difference between our output and a golden is ours to explain.

Runtimes: tier1 68 min per run; tier0 2 min 15 s per run.

## WGS mode: interval scatter (HG002 chr20 slice, run `wgs20-20261010-r1`, 2026-10-10)

Not yet measured with an r2. Three things in WGS mode depend on how the run happened, not only on its inputs. They come from reading sarek's commands and the samtools 1.24 source, and they are recorded per reference run in `reference/<run>.scatter-order.txt`:

| What | Where it comes from | Effect |
|---|---|---|
| GatherBQSRReports `--input` order | Nextflow task-completion order | Read-group `EstimatedQReported` is a floating-point combine of the per-group (rounded) estimates, in this order |
| `samtools merge` input order for `recal.bam` | Nextflow task-completion order | Records that tie on (contig, position, strand) across groups come out in this order |
| `samtools merge` ID suffixes | `hts_srand48(time(NULL))`: the second the task started | Inputs 2..n get `@RG`/`@PG` IDs `ID-%08lX`, and every record's `RG`/`PG` tag takes its input's ID (moved to the end of the tags). Recovered for r1 by searching seeds: 1791621718, the task's start second |

A second sarek run will differ in all three, so its recal BAM records (RG tags) and possibly its recal table will differ from r1. A gate against a different run needs that run's own record of these three.

## Architecture: arm64 (Spark) vs x86-64 (cloud), HG002 chr20 slice (2026-10-10)

sarek 3.10.0 run `wgs20-x86-r1` on an x86 box (2x Xeon, AVX-512; amd64 containers, same per-task resource caps as the Spark) against the arm64 golden `wgs20-20261010-r1`:

| Boundary | x86 vs arm64 | Why |
|---|---|---|
| fastp chunks, 12 mapped chunk BAMs | records identical (16.5M) | bwa's 58 FMA sites (arm64 contracts, x86 builds do not) flip nothing on this data |
| md.bam, md metrics, recal table, recal BAM | identical | the RG `DS:` holds the FASTA path as passed to sarek (`/data/omics/...` vs `/home/...`): header metadata, set with `nomics run --read-group-ds` |
| HaplotypeCaller VCF | 125,121 of 125,153 identical; 1 call only on arm64 | x86 GATK loads Intel GKL: AVX-512 PairHMM (FTZ, float with double fallback) and AVX Smith-Waterman. arm64 has no GKL and uses the Java implementations |
| same, with `--pair-hmm-implementation LOGLESS_CACHING --smith-waterman JAVA` | 125,146 of 125,153 | HotSpot's x86 `UseLibmIntrinsic` stubs compute `Math.log/exp/pow` differently from fdlibm, which arm64 JDK 17 uses (and `jmath` ports) |
| same, plus `-XX:+UnlockDiagnosticVMOptions -XX:-UseLibmIntrinsic` | **125,153 of 125,153 identical** | |

Consequences:
- An x86 reference is equivalent to the arm64 one up to the recal BAM, and for the VCF when HaplotypeCaller runs with the Java PairHMM, Java Smith-Waterman and LIBM intrinsics off. The full 30x WGS golden (`wgs30-x86-r1`) is built that way.
- Our own pipeline is architecture-independent: the x86 build of nomics (CPU-only) gives the same 125,153 records as the arm64 build and as the arm64 golden.

### Whole genome, 30x (`wgs30-x86-r1`, 2026-10-11)

The same check at full scale, on the 21-group HaplotypeCaller scatter of the 30x WGS run:

| VCF | Records | vs Java-mode x86 |
|---|---:|---|
| x86, Java PairHMM + Java SW + LIBM intrinsics off (all 21 groups) | 4,773,797 | (golden) |
| arm64 GATK, the two heaviest groups (chr1_122026460…, chr15_20729747…) | 217,059 + 270,669 | **identical** |
| x86 as sarek ran it (Intel GKL) | 4,773,796 | 4,679,906 identical; 93,882 differ at shared sites; 8 only native, 9 only Java |

The native-vs-Java differences are not spread evenly. 19 groups differ in ~0.01% of records, the chr20 rate. Two groups, `chr4_10001-1429358` (65,966 of 335,092) and `chr9_41229379-41237752` (27,388 of 214,786), differ in 20% and 13%: on chr4, 65,821 of the 65,973 differences are QD alone, starting at one point in the group and persisting to its end. That is the shape of a desync in a group-wide random stream rather than of independent arithmetic noise. It was not investigated further, because the golden is the Java-mode VCF.

On 2026-10-12, seven more groups were recomputed with arm64 GATK (Java PairHMM) to adjudicate nomics' first 30x run: `chr9_41229379`, `chr10_39590436`, `chr17_60001`, `chr20_36314720`, `chr2_238904048`, `chr4_190173122` and `chr5_139453660` (1,684,665 records) are byte-identical to their x86 Java-mode versions. Nine of the 21 groups are now verified on arm64, with no difference found.

The 30x golden VCF is therefore the Java-mode merge: 19 group VCFs computed on x86 and the chr1/chr15 groups recomputed on arm64 (byte-identical to their x86 Java-mode versions), merged with MergeVcfs.
