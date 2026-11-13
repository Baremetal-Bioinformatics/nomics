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

Measured with an r2 on three genomes in November (below). Three things in WGS mode depend on how the run happened, not only on its inputs. They come from reading sarek's commands and the samtools 1.24 source, and they are recorded per reference run in `reference/<run>.scatter-order.txt`:

| What | Where it comes from | Effect |
|---|---|---|
| GatherBQSRReports `--input` order | Nextflow task-completion order | Read-group `EstimatedQReported` is a floating-point combine of the per-group (rounded) estimates, in this order |
| `samtools merge` input order for `recal.bam` | Nextflow task-completion order | Records that tie on (contig, position, strand) across groups come out in this order |
| `samtools merge` ID suffixes | `hts_srand48(time(NULL))`: the second the task started | Inputs 2..n get `@RG`/`@PG` IDs `ID-%08lX`, and every record's `RG`/`PG` tag takes its input's ID (moved to the end of the tags). Recovered for r1 by searching seeds: 1791621718, the task's start second |

**Gather order, measured 2026-11-11.** Re-running GatherBQSRReports on one run's per-group
tables in three input orders (as run, reversed, sorted by name) gave byte-identical gathered
tables, for a run with 18 scatter groups that all hold reads and for the 21-group chr20 slice.
The order changes the arithmetic but did not reach the printed table in any case tried.

A second sarek run will differ in all three, so its recal BAM records (RG tags) and possibly its recal table will differ from r1. A gate against a different run needs that run's own record of these three.

**Measured, 2026-11-07 to 11-09** (sample-matrix genomes G1 HG002 10×, G2 HG001 20×, G3 HG003 20×; `scripts/goldens/genome_golden.sh`). Each genome ran sarek twice. In all three, r1 and r2 produced identical digests for 40 of the 41 artifacts and differed in one, `recal_bam` (same record count, different RG tags: the merge ID suffixes). The recal table and the VCF did not differ in these runs. nomics, given each run's own scatter record, matched r1 41/41 and r2 41/41 on every genome. The scattered exome (E5, 88.5M pairs) shows the same: r1 and r2 match in 17 of 17 judged files, and the recal BAM is not judged between them for the same reason; nomics is gated on it against r1 with r1's recovered seed.

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

## Somatic: Mutect2 tumour-normal (decision run, 2026-11-11)

sarek 3.10.0 with `tools: mutect2`, twice, on the SEQC2 pair HCC1395 (tumour, SRR7890918) and
HCC1395BL (normal) exomes, on the Spark (`scripts/goldens/somatic_golden.sh`; r1 4 h 58 min,
r2 5 h 13 min). The question was whether the Mutect2 chain repeats, so that a port has
something exact to match. It does:

| Output | r1 vs r2 |
|---|---|
| Mutect2 VCF, per scatter group and merged | records identical (6,351), header identical after the usual normalization |
| `.vcf.gz.stats` | bytes identical |
| f1r2 tars, per group | contents identical |
| LearnReadOrientationModel artifact priors (`.artifactprior.tar.gz`) | contents identical (the tar's bytes differ) |
| GetPileupSummaries tables, tumour and normal | bytes identical |
| CalculateContamination contamination and segmentation tables | bytes identical |
| FilterMutectCalls VCF | records identical (6,351; 2,024 PASS), header identical; `.filteringStats.tsv` bytes identical |

Upstream of Mutect2 the two runs differ only where germline runs do: the merged sorted BAM's
tie order, `@PG` order in md.bam, and the recalibrated BAM's `RG` tags (the merge ID suffixes;
no record differs in name, flag or position). Mutect2's output does not depend on those tags.
The recalibration tables are byte-identical. Mutect2 itself is the cost: its two scatter
groups took 45 min and 3 h 9 min, each on one core.

### FilterMutectCalls: a possible order dependence, not observed (2026-11-12)

`ErrorProbabilities` combines the sequencing, artifact and non-somatic error terms as
`1 - (1-a)(1-b)(1-c)`, iterating a `HashMap` keyed by an enum. Enum hash codes are JVM
identity hashes, so the iteration order is not fixed by the source. Floating-point products
are not associative: a worked example gives 0.08450000000000002 in one order and
0.08449999999999991 in another, which print as 0.085 and 0.084. No run has shown it: the two
sarek runs agree in every byte of the filtered VCF and the stats file, GATK gave identical
output on four repeats of two other inputs, and all six orders reproduce the golden through
the port. nomics uses one fixed order (sequencing, artifact, non-somatic). If a future golden
differs from the port in one such rounded digit, this is the first place to look.

### Mutect2 f1r2 archive: column order follows JVM hashes (observed, 2026-11-12)

On the chr21 slices of the HCC1395 pair, two plain Mutect2 runs give the same VCF, stats file
and f1r2 archive members. A third run with only debug outputs added (`--assembly-region-out`,
`--bam-output`) gives the same VCF and stats, the same `ref_histogram` and `alt_table`, but an
`alt_histogram` whose 384 data columns come in another order (`ATT_A_F2R1, ATT_A_F1R2, …`
against `ATT_A_F1R2, ATT_A_F2R1, …`). The values under each column name are the same, and
LearnReadOrientationModel learns identical priors from either archive. The columns are written
in the iteration order of a hash map whose keys include enums, whose hash codes are JVM
identity hashes: stable for one command line, not across options. So an f1r2 gate has to
compare `alt_histogram` by column name, not position: that is the rule in `docs/GATES.md`. Files: `~/omics/runs/somatic/m2-chr21/{r1,r2,oracle}`.


### Mutect2 f1r2 archive: member order follows the directory listing (observed, 2026-11-12)

The f1r2 archive holds three files for a single tumour. Their order in the tar differs
between GATK runs of the same input on different machines:

| Run of unmodified GATK 4.6.2.0 on the chr21 slice | Member order |
|---|---|
| sarek's container on the Spark (the golden) | `.alt_histogram`, `.ref_histogram`, `.alt_table` |
| Homebrew JDK 17.0.21 on the Mac (2026-11-13; `~/omics/runs/somatic/m2-chr21/mac-jdk`) | `.ref_histogram`, `.alt_table`, `.alt_histogram` |

The members themselves agree: `ref_histogram` (43,305 bytes) and `alt_table` (134,856 bytes)
are the same bytes in both archives, and the `alt_histogram`s hold the same values under the
same column labels (their column order is the previous section's subject). The Mac run's VCF
(but for the command line) and stats equal the golden's. On the Spark the order has not
varied: all 14 archives there (nine slice goldens, three repeats, two tasks of the full sarek
run) have the golden's order.

Source: `F1R2CountsCollector.closeAndArchiveFiles` (line 188) passes `tmpDir.listFiles()` to
`IOUtils.writeTarGz`, which writes the array as given (lines 541 to 543). `File.listFiles()`
has no defined order, so the archive inherits the file system's directory enumeration. What
is measured is a difference between two file systems; whether the order can change between
runs on one file system is not known; it did not in those 14.

It does not reach the results: `LearnReadOrientationModel.doWork` extracts each input archive
to its own directory and finds the three kinds of file by suffix, so with one tumour the
order of members inside an archive cannot change what is learned. nomics learns the same
priors (15,362 bytes) from the golden archive and from a copy with its members permuted. The
order of the input archives themselves does matter, in GATK and in the port.

Gate rule (`docs/GATES.md`): the archive's member names are compared as a set. nomics writes
`ref_histogram`, `alt_histogram`, `alt_table`, always in that order.

### Mutect2 `--assembly-region-out`: one run wrote a file with a hole (observed once, 2026-11-13)

`--assembly-region-out` is a debugging output, not one sarek asks for; the slice goldens add it
to compare the port's regions. Of ten GATK 4.6.2.0 runs with it (nine contigs, chr17 twice),
one wrote a damaged file: the first chr17 run's `regions.igv` lacks 18,881 bytes starting at
offset 2,445,312 (a multiple of 4,096), so one line stops after `size=300\t-` and continues
with a record 448 lines later. GATK logged no error and exited 0. The same run's VCF, stats
and f1r2 archive equal the second run's. The second run's file has no hole and is, byte for
byte, the file nomics writes; the first run's file equals it before and after the hole.

Cause not established. The disk was 98% full, with tens of gigabytes free. GATK writes this file
through a `PrintStream` (`AssemblyRegionWalker`), which does not raise write errors, so a failed write would look
like this. The lesson for goldens: a reference file that fails a gate is checked against a
second reference run before the port is changed. The damaged file is kept beside the golden
as `~/omics/runs/somatic/m2-chr17/oracle/regions.igv.damaged-20261112`.
