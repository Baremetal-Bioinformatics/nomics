# Mutect2 first-stage verification

Note (lead, 2026-11-13): this file is the workers' record, stage by stage, in the order the
work was done. Each part describes the port as it stood then: the earlier parts say the
caller is incomplete, that `--f1r2-tar-gz` is not wired, and that gate rows are proposals.
None of that holds now. The current state and the adopted gates are in `docs/GATES.md`; the
lead's acceptance runs (nine contigs, two held out) are recorded there, not here.

Note (lead, 2026-11-12): `scripts/somatic_mutect2_regions_gate.sh` was generalized after this
record was written: it now takes contig names, gates the three external boundaries on every
slice golden present, and no longer compares the internal dumps (the integration test does).
The commands below are as the worker ran them.

2026-11-12, local arm64 Mac. Source revision GATK tag 4.6.2.0,
`76edc75c26504da94bbaee66584e107e76ee15de`; htsjdk 4.2.0
`4cc010022ac038fb30f26e6f9717fabff3e808c1`. Both sources are outside the repo.
No reference code is linked into or called by nomics. Stack map was not edited.

## Exact commands

From the repository root (the fixture paths default to `$HOME/omics`):

```sh
cargo build --release -p nomics
NOMICS_M2_ORACLE_PREFIX=/tmp/nomics-m2-oracle/dump \
  cargo test --release -p nomics-somatic -p nomics-io -p nomics-hc
cargo clippy --release --no-deps -p nomics-somatic -p nomics --all-targets -- -D warnings
NOMICS_M2_ORACLE_PREFIX=/tmp/nomics-m2-oracle/dump \
  scripts/somatic_mutect2_regions_gate.sh /tmp/nomics-m2-port
grep -rnE '\.(powf|powi|ln|exp|log10|log2|ln_1p|exp_m1)\(' \
  crates/nomics-somatic/src/mutect2_regions \
  reference/oracle/mutect2/MutectRegions.java \
  reference/oracle/mutect2/MutectNumerics.java
```

The CLI gate script executes this command, with no stdout filtering:

```sh
target/release/nomics somatic mutect2-regions \
  --input "$HOME/omics/runs/somatic/m2-chr21/HCC1395BL_WES.chr21.bam" \
  --input "$HOME/omics/runs/somatic/m2-chr21/HCC1395_WES.chr21.bam" \
  --normal-sample HCC1395_HCC1395BL_WES \
  --reference "$HOME/omics/ref/GATK.GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta" \
  --panel-of-normals "$HOME/omics/runs/somatic/m2-chr21/1000g_pon.chr21.vcf.gz" \
  --germline-resource "$HOME/omics/runs/somatic/m2-chr21/af-only-gnomad.chr21.vcf.gz" \
  --intervals "$HOME/omics/runs/somatic/m2-chr21/chr21.bed" \
  --assembly-region-out /tmp/nomics-m2-port/regions.igv \
  --stats /tmp/nomics-m2-port/X.stats \
  --dump-prefix /tmp/nomics-m2-port/dump \
  2> /tmp/nomics-m2-port/filters.txt
cmp /tmp/nomics-m2-port/regions.igv "$HOME/omics/runs/somatic/m2-chr21/oracle/regions.igv"
cmp /tmp/nomics-m2-port/X.stats "$HOME"/omics/runs/somatic/m2-chr21/r1/*.vcf.gz.stats
diff -u /tmp/nomics-m2-port/reference-filters.txt /tmp/nomics-m2-port/filters.txt
cmp /tmp/nomics-m2-port/dump.activity.tsv /tmp/nomics-m2-oracle/dump.activity.tsv
cmp /tmp/nomics-m2-port/dump.reads.tsv /tmp/nomics-m2-oracle/dump.reads.tsv
cmp /tmp/nomics-m2-port/dump.prepared.bin /tmp/nomics-m2-oracle/dump.prepared.bin
```

`reference-filters.txt` is extracted by the script from the first numeric count
through `reads processed`, preserving all intervening characters, including
trailing spaces. Only GATK's timestamp/logger prefix is removed. `diff` is empty.
The integration test performs its own extraction and exact comparison.

## Results

All commands succeed; `grep` has zero hits (exit 1, as expected). The Rust suite
reports **157 passed, zero failed, 10 existing ignored**. None were newly ignored.
The new Mutect2 integration test executes, rather than skipping, and compares
all six boundaries when `NOMICS_M2_ORACLE_PREFIX` is set. The 42 Java-produced LLR
bit patterns, smoothing/cut fixture and joint-sample reservoir fixture all pass.
The only build/lint warning is the preexisting HC `Atomic::fetch_update`
deprecation; the test build additionally warns about existing `unused_mut` in
`nomics-hc/tests/dbsnp.rs`. No unrelated code was changed to silence them.

```text
regions cmp: PASS (7167, 699 active)
stats cmp: PASS (35 bytes)
filters diff: PASS
308376 read(s) filtered by: MappingQualityReadFilter
0 read(s) filtered by: MappingQualityAvailableReadFilter
0 read(s) filtered by: MappingQualityNotZeroReadFilter
0 read(s) filtered by: MappedReadFilter
0 read(s) filtered by: NotSecondaryAlignmentReadFilter
716896 read(s) filtered by: NotDuplicateReadFilter
0 read(s) filtered by: PassesVendorQualityCheckReadFilter
0 read(s) filtered by: NonChimericOriginalAlignmentReadFilter
0 read(s) filtered by: NonZeroReferenceLengthAlignmentReadFilter
0 read(s) filtered by: ReadLengthReadFilter
0 read(s) filtered by: GoodCigarReadFilter
0 read(s) filtered by: WellformedReadFilter
1025272 total reads filtered out of 3023196 reads processed
activity.tsv cmp: PASS
reads.tsv cmp: PASS
prepared.bin cmp: PASS
```

The display above omits the otherwise invisible trailing space on each filter
line; compared files retain it. Regions are **7,167/7,167 identical**, **699 active**.
Stats are exactly `statistic\tvalue\ncallable\t1230643.0\n` (35 bytes). Activities:
**1,493,217/1,493,217 loci**, including zeros. Ordered read assignments:
**4,180,613/4,180,613**. Prepared bases/BQ/BI/BD: **2,123,746,776 bytes**, exact.

Local logs from the required commands are `/tmp/m2-build.log`,
`/tmp/m2-tests.log`, `/tmp/m2-clippy.log`, `/tmp/m2-gates.log`,
`/tmp/m2-gates-time.log` and `/tmp/m2-transcendentals.log` (empty).
The final CLI wall time with all dumps enabled was **8.628 seconds**. The
first full scripted gate took 8.836 seconds; timings are in `m2-gates-time.log`.
The first CLI run without the large prepared-read dump took **7.69 seconds**.
These are wall times, not speedup claims; no GPU path is involved.

## Oracle qualification and reproduction

Homebrew OpenJDK **17.0.21.0.3+0**, explicit executables:
`/opt/homebrew/opt/openjdk@17/bin/{java,javac}`. Sarek's recorded container was
17.0.17-internal; it was not used for the local dumps. Before compiling a harness,
the local JVM ran the unmodified jar on the exact BAMs/resources/BED and normal
name, with VCF, stats, f1r2 and `--assembly-region-out`. It matched the golden
regions and stats byte for byte, all filter counts, and all 38 VCF records and
header after excluding only `##GATKCommandLine`. Full-run traversal took 86.0 s.

The f1r2 reference histogram and alt table match byte for byte. The alt-histogram
column order differs; matching columns by label shows identical values. That is
an observation about JVM hash iteration, **not** an approved normalization or a
new f1r2 gate. No f1r2 implementation is part of this change.

Then `MutectRegions.java` ran the unmodified engine, filters, post-transformer,
downsampler and region iterator; its apply method only dumps and runs the
original pre-assembly overlap cleaner. It reproduces the same region and stats
bytes as the full unmodified run. The final harness took 26.4 s to traverse with
all dumps. Read sets use raw SAM flags; reconstructing flags with ReadUtils
would clear mate-strand bits for two rows, which was corrected in the harness.
No masking or port adjustment was used. All final comparisons are exact.

To reproduce this qualification and the dumps:

```sh
JAVA=/opt/homebrew/opt/openjdk@17/bin/java \
JAVAC=/opt/homebrew/opt/openjdk@17/bin/javac \
  reference/oracle/mutect2/run_oracle.sh /tmp/nomics-m2-oracle
```

This script records the same commands as the original runs; its baseline
compares regions, stats, counts and VCF (excluding provenance) before using the
harness. Numeric fixture generation uses `MutectNumerics.java` against that same
jar. `numerics_oracle.tsv` is the only committed JVM output (small deterministic
unit fixtures); BAM/activity/read dumps remain outside the repository.

## Reuse, source details and limits

HC reuse is limited to its interval parser, HcRead, JavaRandom, hard clipper and
unchanged smoothing/cutting implementation through `BinaryActivityProfile`.
The HC caller never invokes that wrapper; hard clipping only gained visibility.
Its enabled tests pass. The ten existing ignored tests were not run, so this
is not a claim to have rerun every HC BAM-level golden.

The source/prose clarifications and all unsupported modes are recorded in
`docs/spec/mutect2.md` §11: start-only header agreement, full htsjdk CIGAR validity,
soft-start mate ordering/reversed ties, signed Java quality bytes, and the
correct raw flag serialization. No default argument discrepancy was found.
The old spec's statement that BAM traversal could not be checked locally is
superseded by these chr21 measurements. Only this default tumour-normal scatter
is validated. Later stages and full-exome/arbitrary-input parity remain out of
scope. The three required gates and three additional internal gates pass;
there are no residual differences in those six comparisons.

Final checks completed at 2026-11-12 22:59:38 UTC, approximately 16 minutes
after branch setup (22:44:00 UTC). All six gates pass on the final code.

# Stage two: region finalization, adaptive assembly and trimming

2026-11-12, local arm64 Mac, branch `somatic/mutect2-assembly` based on
`origin/somatic/caller-a`. No subagents, push, PR or merge. Stack map unchanged.
The pinned source was fetched outside the repository with
`git -C /tmp/mutect2-spec-gatk fetch origin tag 4.6.2.0`; HEAD is
`76edc75c26504da94bbaee66584e107e76ee15de`. No oracle code is linked into nomics.

## Oracle qualification, before trusting assembly dumps

Executed:

```sh
JAVA=/opt/homebrew/opt/openjdk@17/bin/java \
JAVAC=/opt/homebrew/opt/openjdk@17/bin/javac \
  reference/oracle/mutect2/run_oracle.sh /tmp/nomics-m2-stage2-oracle \
  > /tmp/m2-stage2-qualification.log 2>&1
```

The local Homebrew JDK is OpenJDK **17.0.21.0.3+0**. Before the assembly
harness ran, the unmodified sarek jar reproduced chr21's **7,167 regions / 699
active**, stats, filter counts, and **38 VCF records plus header** after removing
only `##GATKCommandLine`. The full Mutect2 run reported 1.57 minutes. The script
compares against `r1/`; a subsequent direct comparison against the requested
`oracle/` VCF/stats/regions also passed exactly under that existing provenance
rule. No JVM-dependent VCF difference was accepted. f1r2 remains out of scope.

`MutectAssembly.java` uses the unmodified walker/engine, filters, transform,
downsampler, overlap cleaner, private unmarked-duplicate remover (reflection),
`assembleReads(...false...)`, EventMap(1), trimmer and `trimTo`. It stops before
likelihoods. Dumps contain traversal ordinal, original/padded spans, full padded
reference, ordered finalized reads (sample/name/raw flags/start/CIGAR/bases and
hex BQ/BI/BD), and both ordered haplotype lists with all requested metadata.
No read/haplotype sorting is applied to the **internal** gate. Artificial BAM
projection follows the original BAM's file order.

The assembly harness was compiled with:

```sh
JAR=$HOME/omics/oracle/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar
/opt/homebrew/opt/openjdk@17/bin/javac -cp "$JAR" -d /tmp/nomics-m2-stage2-oracle \
  reference/oracle/mutect2/MutectAssembly.java \
  reference/oracle/mutect2/AdaptivePruner.java \
  reference/oracle/mutect2/ArtificialHaplotypes.java
D=$HOME/omics/runs/somatic/m2-chr21
/usr/bin/time -p /opt/homebrew/opt/openjdk@17/bin/java -Xmx4g -XX:-UsePerfData \
  -cp "/tmp/nomics-m2-stage2-oracle:$JAR" MutectAssembly \
  --input "$D/HCC1395BL_WES.chr21.bam" --input "$D/HCC1395_WES.chr21.bam" \
  --reference "$HOME/omics/ref/GATK.GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta" \
  --panel-of-normals "$D/1000g_pon.chr21.vcf.gz" \
  --germline-resource "$D/af-only-gnomad.chr21.vcf.gz" --intervals "$D/chr21.bed" \
  --normal-sample HCC1395_HCC1395BL_WES \
  --dump-prefix /tmp/nomics-m2-stage2-oracle/final-assembly \
  --assembly-region-out /tmp/nomics-m2-stage2-oracle/final-assembly.igv \
  > /tmp/m2-stage2-final-oracle.log 2>&1
```

Final harness wall time **31.69 s**. Its region bytes and stats match the full
qualified run. The earlier harness dumped BQ only; extending D rows with clipped
BI/BD and removing unused first-stage dump writers did not alter assembly. Both
versions' comparisons passed; the final required internal dump includes BI/BD.
`run_assembly_oracle.sh OUTDIR` packages qualification and this harness invocation
for reproduction, writing the same format under `assembly.assembly.tsv`.

## Gate commands and outputs

Required build/test/lint and stage-one commands:

```sh
cargo build --release -p nomics > /tmp/m2-stage2-final-build.log 2>&1
NOMICS_M2_ASSEMBLY_ORACLE=/tmp/nomics-m2-stage2-oracle/final-assembly.assembly.tsv \
  cargo test --release -p nomics-somatic -p nomics-io -p nomics-hc \
  > /tmp/m2-stage2-final-tests.log 2>&1
cargo clippy --release --no-deps -p nomics-somatic -p nomics-hc -p nomics --all-targets -- -D warnings \
  > /tmp/m2-stage2-final-clippy.log 2>&1
bash scripts/somatic_mutect2_regions_gate.sh > /tmp/m2-stage2-regions-gate.log 2>&1
```

The three-contig CLI/external gate executed:

```sh
NOMICS_M2_ASSEMBLY_ORACLE=/tmp/nomics-m2-stage2-oracle/assembly.assembly.tsv \
  bash scripts/somatic_mutect2_haplotypes_gate.sh /tmp/nomics-m2-stage2-gate \
  > /tmp/m2-stage2-haplotypes-gate.log 2>&1
```

That script spells out the CLI for each contig and uses the independent Java
`ArtificialHaplotypes` reader on the original BAM. After extending D rows, chr21
was rerun with the same script CLI arguments, using output prefix
`/tmp/nomics-m2-stage2-port/chr21`, and these final comparisons succeeded:

```sh
cmp /tmp/nomics-m2-stage2-port/chr21.assembly.tsv \
    /tmp/nomics-m2-stage2-oracle/final-assembly.assembly.tsv
cmp /tmp/nomics-m2-stage2-port/chr21.artificial.tsv \
    /tmp/nomics-m2-stage2-oracle/chr21.artificial.tsv
cmp /tmp/nomics-m2-stage2-oracle/final-assembly.igv \
    "$HOME/omics/runs/somatic/m2-chr21/oracle/regions.igv"
cmp /tmp/nomics-m2-stage2-oracle/final-assembly.stats \
    /tmp/nomics-m2-stage2-oracle/X.vcf.gz.stats
```

All cmp commands exit zero with no output. Final internal file: 202,832,436 bytes. The integration test repeats the
three external comparisons with nomics-io's own BAM reader and, when the env
variable above is set, compares the final chr21 internal dump too. It never
calls Java. Missing golden directories skip with a message.

Recorded CLI gate stdout:

```text
PASS chr21 artificial
assembly: AssemblyStats { assembled: 699, non_reference: 479, before: 5028, after: 3160, artificial: 3160 }
real 9.15
chr21 golden artificial records: 3160
PASS chr22 artificial
assembly: AssemblyStats { assembled: 1170, non_reference: 850, before: 9555, after: 6883, artificial: 6883 }
real 30.51
chr22 golden artificial records: 6883
PASS chr19 artificial
assembly: AssemblyStats { assembled: 3386, non_reference: 2376, before: 26253, after: 19219, artificial: 19219 }
real 197.31
chr19 golden artificial records: 19219
PASS chr21 internal
```

The final chr21 CLI with BI/BD dumping took **10.80 s**. The preceding timings
were under concurrent local validation load, not benchmarks or speedup claims.
All artificial fields and tags (`HC:i` base hash, `RG:Z`, `CR:Z` callable span)
match. No artificial tag depends on a later stage in this fixed mode; none was
excluded. Realigned biological records are deliberately outside this gate.

The stage-one command reports:

```text
PASS chr21 (7167 regions, 699 active; 1230643.0 callable; 1025272 / 3023196 filtered/processed; 7 s)
PASS chr22 (13815 regions, 1170 active; 2689447.0 callable; 2445097 / 8935845 filtered/processed; 25 s)
PASS chr19 (34338 regions, 3386 active; 6845249.0 callable; 5481759 / 20756074 filtered/processed; 202 s)
mutect2 regions gate: 3 pass, 0 fail, 0 skipped
```

Each line abbreviates the script's full filter-summary wording; its comparisons
are the unmodified stage-one script's region/stats bytes and all filter counts.

## Required checks, reuse and limits

Build: **PASS**, final release build 27.76 s; existing HC `fetch_update`
deprecation warning only. The earlier complete three-crate test run reports
**159 passed, zero failed, 10 existing ignored**, including all three external
assembly goldens. The final run with BI/BD and the final internal oracle enabled also reports
**159 passed, zero failed, 10 existing ignored**; the three-contig integration
test passes in **136.55 s**. All comparisons ran (no golden skips).

**Required full Clippy command: FAIL**, 112 existing errors in unchanged HC
files: activity, annotate, assembly/cigar, call, caller, dbsnp, likelihood,
pairhmm_trie, reads and refconfidence. The first is `annotate.rs:625`
`clippy::int_plus_one`; `caller.rs:340` also fails on the existing deprecated
`Atomic::fetch_update`. No lint finding references a modified HC file.
These unrelated files were not edited or suppressed to manufacture a pass.
The supplemental command
`cargo clippy --release --no-deps -p nomics-somatic -p nomics --all-targets -- -D warnings`
**passes** (2.13 s), checking the new implementation, CLI and integration test.

The requested forbidden-transcendental grep was executed over every touched or
added Rust implementation/test and Java harness source (graph.rs, assembly/mod.rs,
finalize.rs, trim.rs, somatic/lib.rs, mutect2_regions/numerics.rs,
mutect2_assembly/, mutect2_haplotypes_goldens.rs, nomics/main.rs, and the three
new .java files). Pattern:
`\.(powf|powi|ln|exp|log10|log2|ln_1p|exp_m1)\(`. **Zero hits**, exit 1 as expected;
`/tmp/m2-stage2-transcendentals.log` is empty. The compact final command was:

```sh
git diff --name-only origin/somatic/caller-a -- '*.rs' '*.java' | \
  xargs grep -rnE '\.(powf|powi|ln|exp|log10|log2|ln_1p|exp_m1)\(' \
  > /tmp/m2-stage2-transcendentals.log
```

`git diff --check` and shell syntax
checks on both added scripts pass.

Shared HC reuse/unchanged defaults and all source/prose refinements are in spec
§12. The HC entry points still use overlap correction true, fixed pruning and
dangling cutoff 2, one sample, and MNP distance 0. New entry points enable the
somatic settings without changing HC callers. HC tests passed (32 enabled tests); nine existing HC
extended fixture tests remained ignored, alongside one ignored IO test, so no claim is made to have rerun those
unavailable oracles. The assembler's SW/CIGAR, simplification, k-best and retry
code is shared, not copied.

Independent pruner fixtures: **25 Java graphs / 232 chains and 96 exact LLR bit
patterns**. The large graph has 112 chains (110 branches plus two flanks), of
which eight are removed; it checks the <=100 loop and already-added 101st
variant behavior. Tests also cover equal multiplicities, reference protection
and Q93 at zero estimated error. Generate with `AdaptivePruner.java`; committed
output is only the small `pruner_oracle.tsv`, never a full BAM/dump.

First differing region during development: **2477, chr21:33296904-33297000**.
All read/reference rows matched; the first haplotype score/order differed
(-0.8293374103867619 vs -0.8845183800306204, with a different leading path).
Tracing raw graph edge counts to source identified **per-k-mer usable-run sample
insertion**, not an LLR threshold. Before fixing this there were three differing
assembly regions (2477, 2901, 6034); after fixing it, all chr21 internal bytes
match and the external held-out contigs pass. No threshold tuning, sorting-away
of internal differences, normalization or residual numerical tolerance was used.

Not confirmed: arbitrary fully comparator-equal distinct chains and JVM identity
hash ordering, nondefault assembly modes, and extended ignored HC fixtures.
Unsupported CLI modes and exact boundaries are listed in spec §12.1. The passing
contigs establish these gates, not arbitrary-graph/JVM or full-exome parity.
The GATES additions are explicitly **proposals awaiting the lead**.

## Stage three: read likelihoods and realignment (2026-11-12)

Implementation branch: `somatic/mutect2-likelihoods`, based on
`origin/somatic/caller-b`. No push, PR or merge. Source tag
4.6.2.0 was refreshed with `git -C /tmp/mutect2-spec-gatk fetch origin tag
4.6.2.0`; its HEAD is `76edc75c26504da94bbaee66584e107e76ee15de`.
htsjdk 4.2.0 was read in `/tmp/mutect2-spec-htsjdk`.

### Oracle qualification and commands

Homebrew JDK **17.0.21.0.3+0**, explicit Java/javac under
`/opt/homebrew/opt/openjdk@17/bin`. Unmodified sarek jar SHA256:
`40b494b1e356931c9143be09c57ba6dfbb8d150167361aac2ceec24089b75667`.
Golden `oracle/gatk.log` reports native OpenMP AVX unsupported and fallback to
**LOGLESS_CACHING**. The same fallback appears in the local logs. This is a
Java-mode gate, not Intel GKL/NEON/GPU parity.

Before trusting likelihood dumps, ran:

```sh
JAVA=/opt/homebrew/opt/openjdk@17/bin/java \
JAVAC=/opt/homebrew/opt/openjdk@17/bin/javac \
  bash reference/oracle/mutect2/run_oracle.sh /tmp/m2-likelihood-oracle \
  > /tmp/m2-likelihood-qualification.log 2>&1
```

Output:

```text
Unmodified JVM: regions, stats, filter counts exact; VCF exact except command provenance
Harness: stage goldens exact. Dumps and numerical fixtures ready.
```

The full unmodified Mutect2 run took **87 seconds / 1.45 minutes** according
to its log. Its VCF check removes only `##GATKCommandLine=`. The qualification
also compares the region file and stats bytes and all read-filter counts.
This JVM is not the container's 17.0.17-internal build. No f1r2 equality claim
is made; its known histogram ordering issue is outside this task.

Then compiled and ran the stage-three walker (the public-engine consistency
assertion was subsequently added and this same command rerun successfully):

```sh
jar=$HOME/omics/oracle/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar
out=/tmp/m2-likelihood-oracle
root=$HOME/omics/runs/somatic/m2-chr21
ref=$HOME/omics/ref/GATK.GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta
/opt/homebrew/opt/openjdk@17/bin/javac -cp "$jar" -d "$out" \
  reference/oracle/mutect2/MutectLikelihoods.java
/usr/bin/time -p /opt/homebrew/opt/openjdk@17/bin/java -Xmx4g -XX:-UsePerfData \
  -cp "$out:$jar" MutectLikelihoods \
  --input "$root/HCC1395BL_WES.chr21.bam" --input "$root/HCC1395_WES.chr21.bam" \
  --reference "$ref" --panel-of-normals "$root/1000g_pon.chr21.vcf.gz" \
  --germline-resource "$root/af-only-gnomad.chr21.vcf.gz" --intervals "$root/chr21.bed" \
  --normal-sample HCC1395_HCC1395BL_WES --dump-prefix "$out/likelihood" \
  --assembly-region-out "$out/likelihood.igv" > "$out/likelihood.log" 2>&1
cmp "$out/likelihood.igv" "$root/oracle/regions.igv"
cmp "$out/likelihood.stats" "$out/X.vcf.gz.stats"
```

Initial harness wall time **83.19 s**; final run with the independent public
`computeReadLikelihoods` assertion **153.96 s**. All assertions and both `cmp`s
succeeded. The harness invokes unmodified jar code, including private raw HMM
and preparation methods through reflection; it does not instrument/rebuild
the jar. `reference/oracle/mutect2/run_likelihood_oracle.sh` packages these
qualification/dump/fixture steps for reuse. The individual commands above are
the ones executed during development.

Independent small fixtures were produced by compiling/running
`MutectLikelihoodNumerics.java` and `MutectLikelihoodTies.java` with the same
explicit javac/Java and classpath. The five numerical rows contain input and
prepared BQ/BI/BD, raw PairHMM bits and ln bits. The eight tie rows contain
input difference bits, chosen allele, informative flag and confidence bits;
the exact ln-window endpoint and adjacent doubles are included. Fixtures are
committed under `mutect2_likelihoods/numerics.tsv` and HC's
`tests/fixtures/mutect2-ties.tsv`. JVM outputs compare exactly to those files.
The fixture generation calls Java `Math` functions; nomics uses its existing
Java-exact math ports, never platform Rust transcendental functions.

### Rust commands and exact gate commands

```sh
cargo build --release -p nomics > /tmp/m2-lk-build-final.log 2>&1
NOMICS_M2_LIKELIHOODS_ORACLE=/tmp/m2-likelihood-oracle/likelihood.likelihoods.tsv \
  cargo test --release -p nomics-somatic -p nomics-io -p nomics-hc -- --nocapture \
  > /tmp/m2-lk-tests-final.log 2>&1
cargo clippy --release --no-deps -p nomics-somatic -p nomics --all-targets -- -D warnings \
  > /tmp/m2-lk-clippy-final.log 2>&1
bash scripts/somatic_mutect2_regions_gate.sh > /tmp/m2-lk-regions-gate.log 2>&1
bash scripts/somatic_mutect2_haplotypes_gate.sh /tmp/m2-lk-haplotypes \
  > /tmp/m2-lk-haplotypes-gate.log 2>&1
NOMICS_BUILD_DIR="$PWD/target" \
NOMICS_RG_DS=<reference path of the golden run> \
  scripts/mac_tier0.sh cpu /tmp/m2-lk-tier0 > /tmp/m2-lk-tier0.log 2>&1
NOMICS_M2_LIKELIHOODS_ORACLE=/tmp/m2-likelihood-oracle/likelihood.likelihoods.tsv \
  bash scripts/somatic_mutect2_likelihoods_gate.sh /tmp/m2-lk-port-final \
  > /tmp/m2-lk-gates-final.log 2>&1
cmp /tmp/m2-lk-port-final/chr21.likelihoods.tsv \
  /tmp/m2-likelihood-oracle/likelihood.likelihoods.tsv
```

The likelihood gate script runs this exact CLI for each `c` in
`chr21 chr22 chr19 chr20 chr13`, then independently projects the golden BAM
using `BiologicalReads.java` and compares bytes:

```sh
d=$HOME/omics/runs/somatic/m2-$c
out=/tmp/m2-lk-port-final
/usr/bin/time -p target/release/nomics somatic mutect2-likelihoods \
  --input "$d/HCC1395BL_WES.$c.bam" --input "$d/HCC1395_WES.$c.bam" \
  --normal-sample HCC1395_HCC1395BL_WES --reference "$ref" --intervals "$d/$c.bed" \
  --panel-of-normals "$d/1000g_pon.$c.vcf.gz" --germline-resource "$d/af-only-gnomad.$c.vcf.gz" \
  --assembly-region-out "$out/$c.igv" --stats "$out/$c.stats" \
  --output "$out/$c.likelihoods.tsv" --realigned-records "$out/$c.realigned.tsv" \
  2> "$out/$c.log"
/opt/homebrew/opt/openjdk@17/bin/java -cp "$out:$jar" BiologicalReads \
  "$d/oracle/haplotypes.bam" > "$out/$c.golden.tsv"
cmp "$out/$c.realigned.tsv" "$out/$c.golden.tsv"
```

For the requested transcendental grep, the added/changed Rust paths were
collected with `git diff --name-only` plus `git ls-files --others
--exclude-standard`, then:

```sh
grep -rnE '\.(powf|powi|ln|exp|log10|log2|ln_1p|exp_m1)\(' \
  crates/nomics-io/src/bam_query.rs crates/nomics-hc/src/{finalize,readprep,realign}.rs \
  crates/nomics-somatic/src/{lib.rs,mutect2_assembly/mod.rs,mutect2_regions/mod.rs,mutect2_regions/reads.rs} \
  crates/nomics-somatic/src/mutect2_likelihoods \
  crates/nomics-somatic/tests/mutect2_likelihoods_goldens.rs crates/nomics/src/main.rs
```

No hits; exit 1 as expected. `/tmp/m2-lk-transcendentals.log` is empty.

### Source discrepancies, exclusions and implementation choices

No default likelihood rule in §5 disagreed with source. §5/§10's uncertainty
about the actual PairHMM fallback is superseded by the golden log. Newly
confirmed details are written in spec §13: actual-length `<10` stubs, signed
quality bytes, overlap-tag presence even after clipping, sample/header order,
PG/RG merge translations and the BAM writer's late position in `callRegion`.
The first external development difference was chr21's first biological record
`SRR7890918.102799831`, PG `MarkDuplicates` versus `MarkDuplicates.P`, affecting
44,726 tumour records. Tracing htsjdk's ancestry/counter algorithm fixed it;
no constant or tag exclusion was fitted to conceal it. There was no internal
numerical development difference on chr21.

All stage-determined biological fields and tags are compared. AR/CR are attached inside genotyping but their span values are fixed here and
are included. Only XA is excluded because it depends on called alleles.
The exact tag list and source functions are in spec §13 and proposed gates.
No later genotyper filtering removes matrix reads in this default mode. No
record is deduplicated across regions, and golden BAM order is preserved.

Reused HC's portable prefix-sharing PairHMM, PCR cache/repeat finder,
`likelihood::finish`, hard clipping, coordinate comparator and SW/CIGAR read
realignment. New quality preparation accepts supplied BI/BD; its Java signed
byte semantics do not change HC's existing entry point. Realignment adds an
ln-aware entry point using a parameterized tie search; HC still supplies 0.2.
Assembly adds only a synchronous varying-region callback. Stage one carries
a Boolean for overlap-materialized BI/BD tags and an exact input-file/BGZF
end-offset identity, without retaining full BAM payloads; only stage three performs the second tag-recovery I/O pass.
No SIMD/GPU code or new dependencies were introduced.

Limits: default two-indexed-BAM/BED sarek mode only. No fragments, genotyping,
annotations, VCF, f1r2, native GKL, GPU, flow/collapse, custom likelihood settings,
DRAGstr, forced alleles, mitochondria, gVCF or other sample configurations.
Flow read groups are explicit errors. Input-file/virtual-offset identity
distinguishes even different records with identical sample/name/flag. Exotic auxiliary tag types are supported
by the Rust projection but unconfirmed by the simple JVM projection (which
fails explicitly instead of ignoring one). SW offset -1 and arbitrary graph
identity-hash ordering remain unconfirmed. See spec §§12.1/13.3. No speedup
claim is made; timings below include local verification contention.

An initial chr20 development run exposed why sample/name/flag is not an input
record identity: normal read `SRR7890919.48708291`, flag 2209, has separate
supplementary records at chr20:23627891 (85S41M) and chr20:23627934 (32M94S).
In particular, multiple source records can share all three. That run stopped
explicitly instead of selecting one. The final port retains the source file
ordinal and the BGZF virtual offset after each input record; tag recovery
uses that exact key. `BamQuery::virtual_offset` is a read-only accessor, with no
change to query traversal. The final runs below use this corrected identity.
Earlier redundant test runs were stopped and replaced by the final full run.
AR/CR were also added to the final comparison after auditing that their values,
although assigned inside genotyping, are entirely fixed by this boundary.

### Stage-three results

All five final biological-record comparisons and the chr21 internal byte
comparison pass. The latter is **166,025,968 bytes**, containing **1,253,491**
raw and capped binary64 likelihoods each, and **1,213,461** retained log10 and
ln values each. No tolerance or last-bit rounding is applied.

| Contig | Likelihood regions | Stubs dropped | Reads evaluated | Poorly modelled dropped | Realigned records compared | CLI wall s |
|---|---:|---:|---:|---:|---:|---:|
| chr21 | 364 | 7,288 | 97,758 | 2,063 | 95,695 | 20.02 |
| chr22 | 702 | 22,355 | 267,334 | 5,303 | 262,031 | 59.99 |
| chr19 | 1,909 | 48,127 | 576,806 | 13,628 | 563,178 | 164.60 |
| chr20 | 801 | 20,740 | 256,440 | 5,561 | 250,879 | 60.37 |
| chr13 | 596 | 13,299 | 168,802 | 3,452 | 165,350 | 35.29 |

Counts below are **reads affected** per rule, except the last column, which
counts matrix cells capped. These categories overlap; they are not summed as
disjoint reads. Quality caps/floors do not drop a read. Counts include reads
subsequently disqualified. “Realigned records” counts every retained output
record, including alignments whose coordinates do not change.

| Contig | MQ BQ cap | PCR BI/BD cap | BQ floor | BI/BD floor | Likelihood cap reads | Likelihood cap cells |
|---|---:|---:|---:|---:|---:|---:|
| chr21 | 1,222 | 63,610 | 24,480 | 0 | 29,283 | 374,370 |
| chr22 | 3,994 | 166,747 | 64,598 | 0 | 75,285 | 1,127,987 |
| chr19 | 3,920 | 356,299 | 150,655 | 0 | 169,610 | 2,810,859 |
| chr20 | 1,004 | 164,133 | 65,239 | 0 | 66,959 | 1,173,622 |
| chr13 | 1,291 | 110,055 | 42,056 | 0 | 51,670 | 780,463 |

Exact final gate output (`/tmp/m2-lk-gates-final.log`):

```text
PASS chr21 biological records
likelihoods: LikelihoodStats { regions: 364, stubs: 7288, evaluated: 97758, mapq_capped_reads: 1222, base_quality_floored_reads: 24480, pcr_capped_reads: 63610, indel_quality_floored_reads: 0, capped_reads: 29283, capped_values: 374370, poorly_modelled: 2063, realigned: 95695 }
real 20.02
chr21 golden records:    95695
PASS chr22 biological records
likelihoods: LikelihoodStats { regions: 702, stubs: 22355, evaluated: 267334, mapq_capped_reads: 3994, base_quality_floored_reads: 64598, pcr_capped_reads: 166747, indel_quality_floored_reads: 0, capped_reads: 75285, capped_values: 1127987, poorly_modelled: 5303, realigned: 262031 }
real 59.99
chr22 golden records:   262031
PASS chr19 biological records
likelihoods: LikelihoodStats { regions: 1909, stubs: 48127, evaluated: 576806, mapq_capped_reads: 3920, base_quality_floored_reads: 150655, pcr_capped_reads: 356299, indel_quality_floored_reads: 0, capped_reads: 169610, capped_values: 2810859, poorly_modelled: 13628, realigned: 563178 }
real 164.60
chr19 golden records:   563178
PASS chr20 biological records
likelihoods: LikelihoodStats { regions: 801, stubs: 20740, evaluated: 256440, mapq_capped_reads: 1004, base_quality_floored_reads: 65239, pcr_capped_reads: 164133, indel_quality_floored_reads: 0, capped_reads: 66959, capped_values: 1173622, poorly_modelled: 5561, realigned: 250879 }
real 60.37
chr20 golden records:   250879
PASS chr13 biological records
likelihoods: LikelihoodStats { regions: 596, stubs: 13299, evaluated: 168802, mapq_capped_reads: 1291, base_quality_floored_reads: 42056, pcr_capped_reads: 110055, indel_quality_floored_reads: 0, capped_reads: 51670, capped_values: 780463, poorly_modelled: 3452, realigned: 165350 }
real 35.29
chr13 golden records:   165350
PASS chr21 internal
```

### Upstream and shared-code regression results

Final build: `Finished release profile [optimized + debuginfo]` in **31.55 s**.
Strict requested Clippy command succeeds in **3.31 s**. The pre-existing HC
`Atomic::fetch_update` deprecation remains a dependency warning; the full test
build also reports the existing unused-mut warning in `tests/dbsnp.rs`. No
unrelated warning was suppressed.

The requested stage-one script output:

```text
PASS chr21  (7167 regions, 699 active; 1230643.0 callable; 1025272 total reads filtered out of 3023196 reads processed; 9 s)
PASS chr22  (13815 regions, 1170 active; 2689447.0 callable; 2445097 total reads filtered out of 8935845 reads processed; 62 s)
PASS chr19  (34338 regions, 3386 active; 6845249.0 callable; 5481759 total reads filtered out of 20756074 reads processed; 206 s)
PASS chr20  (15433 regions, 1367 active; 3161042.0 callable; 2283986 total reads filtered out of 8908067 reads processed; 27 s)
PASS chr13  (11531 regions, 1166 active; 2157790.0 callable; 1458591 total reads filtered out of 5253705 reads processed; 16 s)
mutect2 regions gate: 5 pass, 0 fail, 0 skipped
```

The requested stage-two script output:

```text
PASS chr21 artificial
assembly: AssemblyStats { assembled: 699, non_reference: 479, before: 5028, after: 3160, artificial: 3160 }
real 10.11
chr21 golden artificial records:     3160
PASS chr22 artificial
assembly: AssemblyStats { assembled: 1170, non_reference: 850, before: 9555, after: 6883, artificial: 6883 }
real 32.69
chr22 golden artificial records:     6883
PASS chr19 artificial
assembly: AssemblyStats { assembled: 3386, non_reference: 2376, before: 26253, after: 19219, artificial: 19219 }
real 94.76
chr19 golden artificial records:    19219
PASS chr20 artificial
assembly: AssemblyStats { assembled: 1367, non_reference: 962, before: 11368, after: 7802, artificial: 7802 }
real 32.37
chr20 golden artificial records:     7802
PASS chr13 artificial
assembly: AssemblyStats { assembled: 1166, non_reference: 777, before: 7810, after: 5244, artificial: 5244 }
real 18.46
chr13 golden artificial records:     5244
SKIP chr21 internal: set NOMICS_M2_ASSEMBLY_ORACLE to qualified Java dump
```

The script skipped its optional internal dump check because that environment
variable was not set. It was then checked explicitly, successfully:

```sh
cmp /tmp/m2-lk-haplotypes/chr21.assembly.tsv \
  /tmp/nomics-m2-stage2-oracle/final-assembly.assembly.tsv
```

The final stage-three CLI regions/stats/filter summaries were also compared
to their golden files on all five contigs after the source-identity change
(`/tmp/m2-lk-upstream-from-final.log`):

```text
PASS chr21: final CLI regions, stats, all filter counts
PASS chr22: final CLI regions, stats, all filter counts
PASS chr19: final CLI regions, stats, all filter counts
PASS chr20: final CLI regions, stats, all filter counts
PASS chr13: final CLI regions, stats, all filter counts
```

Tier0 result: **41/41 digests PASS**, wall **32 s**, peak RSS **14.4 GB**.
The exact output tail was:

```text

== cpu: tier0 FASTQ -> VCF, 14 threads -> /tmp/m2-lk-tier0
== cpu: wall 32s, peak RSS 14.4 GB
PASS	recal_table	count=1090
PASS	vcf	count=1047
verdict: PASS
```

Across all five contigs, the complete compared tag set is:
`AR, AS, BD, BI, CR, HC, MC, MD, MQ, NM, PG, RG, SA, XS`. Only `XA` is excluded
for its dependency on somatic calls. Every projected field and tag value is
identical; there is no residual stage-three difference on these gates.

### Final Rust suite output

The final required Cargo test command completed with exit 0: **162 passed,
0 failed, 10 existing ignored**, across 31 test groups (including doc tests).
No tests were newly ignored. The two new independent numerical unit tests
pass. The new integration test executes all five contigs and the qualified
chr21 internal dump check; it does not skip them.

```text
test realign::tests::mutect2_natural_log_ties_match_independent_jar ... ok
test mutect2_likelihoods::tests::independently_generated_gatk_quality_and_likelihood_bits ... ok
chr21: LikelihoodStats { regions: 364, stubs: 7288, evaluated: 97758, mapq_capped_reads: 1222, base_quality_floored_reads: 24480, pcr_capped_reads: 63610, indel_quality_floored_reads: 0, capped_reads: 29283, capped_values: 374370, poorly_modelled: 2063, realigned: 95695 }; golden 95695; 21.049s
chr22: LikelihoodStats { regions: 702, stubs: 22355, evaluated: 267334, mapq_capped_reads: 3994, base_quality_floored_reads: 64598, pcr_capped_reads: 166747, indel_quality_floored_reads: 0, capped_reads: 75285, capped_values: 1127987, poorly_modelled: 5303, realigned: 262031 }; golden 262031; 63.179s
chr19: LikelihoodStats { regions: 1909, stubs: 48127, evaluated: 576806, mapq_capped_reads: 3920, base_quality_floored_reads: 150655, pcr_capped_reads: 356299, indel_quality_floored_reads: 0, capped_reads: 169610, capped_values: 2810859, poorly_modelled: 13628, realigned: 563178 }; golden 563178; 164.304s
chr20: LikelihoodStats { regions: 801, stubs: 20740, evaluated: 256440, mapq_capped_reads: 1004, base_quality_floored_reads: 65239, pcr_capped_reads: 164133, indel_quality_floored_reads: 0, capped_reads: 66959, capped_values: 1173622, poorly_modelled: 5561, realigned: 250879 }; golden 250879; 65.518s
chr13: LikelihoodStats { regions: 596, stubs: 13299, evaluated: 168802, mapq_capped_reads: 1291, base_quality_floored_reads: 42056, pcr_capped_reads: 110055, indel_quality_floored_reads: 0, capped_reads: 51670, capped_values: 780463, poorly_modelled: 3452, realigned: 165350 }; golden 165350; 38.020s
test mutect2_stage_three_available_goldens ... ok
```

All test-group result lines from `/tmp/m2-lk-tests-final.log`:

```text
test result: ok. 18 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 1.52s
test result: ok. 0 passed; 0 failed; 1 ignored; 0 measured; 0 filtered out; finished in 0.00s
test result: ok. 0 passed; 0 failed; 1 ignored; 0 measured; 0 filtered out; finished in 0.00s
test result: ok. 0 passed; 0 failed; 1 ignored; 0 measured; 0 filtered out; finished in 0.00s
test result: ok. 1 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.00s
test result: ok. 0 passed; 0 failed; 1 ignored; 0 measured; 0 filtered out; finished in 0.00s
test result: ok. 1 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.07s
test result: ok. 2 passed; 0 failed; 1 ignored; 0 measured; 0 filtered out; finished in 0.00s
test result: ok. 3 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.00s
test result: ok. 3 passed; 0 failed; 1 ignored; 0 measured; 0 filtered out; finished in 0.34s
test result: ok. 1 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.02s
test result: ok. 0 passed; 0 failed; 1 ignored; 0 measured; 0 filtered out; finished in 0.00s
test result: ok. 1 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 1.00s
test result: ok. 1 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.08s
test result: ok. 0 passed; 0 failed; 1 ignored; 0 measured; 0 filtered out; finished in 0.00s
test result: ok. 2 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.05s
test result: ok. 0 passed; 0 failed; 1 ignored; 0 measured; 0 filtered out; finished in 0.00s
test result: ok. 79 passed; 0 failed; 1 ignored; 0 measured; 0 filtered out; finished in 0.54s
test result: ok. 10 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.20s
test result: ok. 3 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.02s
test result: ok. 27 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 4.39s
test result: ok. 2 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 3.49s
test result: ok. 2 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 2.48s
test result: ok. 1 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 143.47s
test result: ok. 1 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 352.75s
test result: ok. 1 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 7.09s
test result: ok. 1 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 2.42s
test result: ok. 2 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 4.02s
test result: ok. 0 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.00s
test result: ok. 0 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.00s
test result: ok. 0 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.00s
```

The new integration test uses the same Rust projection codec for the BAM
input, while the CLI gate independently uses htsjdk/Java for that projection.
Both pass. No new code changes were made after these successful checks; final
documentation records the commands and outputs.

# Stage four, partial: fragments and pre-annotation somatic evidence

2026-11-12, local arm64 Mac, branch `somatic/mutect2-calls` based on
`origin/somatic/caller-c`. **The requested final caller is not complete.**
There is no `nomics somatic mutect2`, VCF writer, f1r2 collector/archive writer,
index or XA implementation. The new command is explicitly named
`mutect2-genotyping`, and its output is a diagnostic dump, not a VCF. No final
output parity or full internal annotation gate is claimed. No push or PR was made. No reference tool is invoked by nomics.

## Oracle qualification, fixtures and exact internal commands

Pinned source checkout reused outside the repo: `/tmp/mutect2-spec-gatk`,
`git rev-parse HEAD` = `76edc75c26504da94bbaee66584e107e76ee15de`,
`git describe --tags --exact-match` = `4.6.2.0`; htsjdk 4.2.0 and Commons Math
3.5 sources are also outside the repo. JDK executables are explicitly
`/opt/homebrew/opt/openjdk@17/bin/{java,javac}` (17.0.21).

```sh
JAVA=/opt/homebrew/opt/openjdk@17/bin/java \
JAVAC=/opt/homebrew/opt/openjdk@17/bin/javac \
  bash reference/oracle/mutect2/run_oracle.sh /tmp/m2-stage4-oracle \
  > /tmp/m2-stage4-qualification.log 2>&1
bash reference/oracle/mutect2/run_calls_oracle.sh /tmp/m2-stage4-final-oracle \
  > /tmp/m2-stage4-final-oracle.log 2>&1
```

Both qualification runs report:

```text
Unmodified JVM: regions, stats, filter counts exact; VCF exact except command provenance
Harness: stage goldens exact. Dumps and numerical fixtures ready.
```

The fresh unmodified jar reproduced chr21's header and **38/38** VCF records,
excluding only `##GATKCommandLine`, the **35 stats bytes**, regions and all
filter counts, before stage-four dumps were trusted. The final script also
runs the complete unmodified engine with an instrumented SomaticGenotypingEngine:

```text
PASS instrumented chr21: VCF except command provenance, stats bytes; actual/projected LOD assertions executed
```

The test-only walker delegates assembly, likelihoods, annotation and all call
production to the jar. Its genotyper subclass independently projects fragment
and event matrices, and asserts that every **actual** call to the protected
somaticLogOdds method returns exactly the predicted raw bits in the predicted
order. Normal diploid LOD comes from the jar's actual private method by
reflection. The retained original engine is called to emit the real VCF/stats.
A separate `.annotations.tsv` file serializes actual returned INFO/genotype
attributes before formatting, with doubles as raw bits; Rust does not yet
produce its counterpart. The diagnostic projection copies likelihood values
and does not substitute its results into the real caller.

`MutectSomaticNumerics.java` supplies **80** independent matrices (1–5 alleles,
1–70 evidence columns, flat and nonflat priors, weak/strong evidence and
-Infinity values). Expected bits cover posterior parameters, variational
logEvidence, Dirichlet normalization and actual somatic/normal diploid LOD
methods. `MutectFragments.java` supplies **145 reads**, exercising String hash
collisions/tree bins/resizes and the one/two, >2, single-clean and all-flagged
representative-read cases. Its warning says “randomly”, but the chosen reads
are source-defined first survivors; no RNG draw occurs there. Both fixtures
are committed as TSVs and byte-compared when regenerated by the final script.

Exact Rust diagnostic command and comparison:

```sh
mkdir -p /tmp/m2-stage4-port
/usr/bin/time -p target/release/nomics somatic mutect2-genotyping \
 --input "$HOME/omics/runs/somatic/m2-chr21/HCC1395BL_WES.chr21.bam" \
 --input "$HOME/omics/runs/somatic/m2-chr21/HCC1395_WES.chr21.bam" \
 --normal-sample HCC1395_HCC1395BL_WES \
 --reference "$HOME/omics/ref/GATK.GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta" \
 --panel-of-normals "$HOME/omics/runs/somatic/m2-chr21/1000g_pon.chr21.vcf.gz" \
 --germline-resource "$HOME/omics/runs/somatic/m2-chr21/af-only-gnomad.chr21.vcf.gz" \
 --intervals "$HOME/omics/runs/somatic/m2-chr21/chr21.bed" \
 --assembly-region-out /tmp/m2-stage4-port/chr21.igv \
 --stats /tmp/m2-stage4-port/chr21.stats \
 --output /tmp/m2-stage4-port/chr21.calls.tsv --likelihoods /dev/null \
 --realigned-records /tmp/m2-stage4-port/chr21.reads.tsv \
 > /tmp/m2-stage4-port/chr21.stdout 2> /tmp/m2-stage4-port/chr21.log
cmp /tmp/m2-stage4-port/chr21.calls.tsv /tmp/m2-stage4-final-oracle/chr21.calls.tsv
```

`cmp` exits **0**, no output. **38,757,793 bytes** identical. The dump covers
**364** varying regions, **58,307** ordered fragments, **862** events, **1,724**
per-sample marginal matrices, every TLOD/NLOD/NALOD in raw binary64 bits and
all emission decisions. **38** events pass. This is the **pre-annotation
subset**, not the full requested internal gate. The measured Rust command wall
time is **19.86 s**, user 19.09 s, sys 0.61 s. The instrumented JVM run with
BAM/f1r2 output took **95.44 s**, user 107.88 s, sys 1.91 s. These are wall times,
not speedup claims; no GPU path was added.

## Reference archive-order discrepancy (no new normalization)

The final gate's Python comparator was also executed directly on the fresh
**JVM** output (not Rust output), with these arguments:

```sh
python3 - <<'PY'
from pathlib import Path
import sys
source = Path('scripts/somatic_mutect2_gate.sh').read_text()
code = source.split("<<'PY'", 1)[1].split('\nPY', 1)[0]
sys.argv = ['comparator', str(Path.home()/'omics/runs/somatic/m2-chr21'),
            '/tmp/m2-stage4-final-oracle', 'chr21']
exec(compile(code, 'somatic_mutect2_gate.sh comparator', 'exec'))
PY
```

It correctly fails, exit **1**, after:

```text
PASS chr21: records 38/38, header except GATKCommandLine
PASS chr21: stats 35 bytes
AssertionError: archive member names/order: first difference at 0, lengths 3/3
```

The golden member suffix order is `.alt_histogram, .ref_histogram, .alt_table`.
Both the fresh **unmodified** qualification archive and the instrumented archive
have `.ref_histogram, .alt_table, .alt_histogram`. Names themselves are identical,
including the `./HCC1395_HCC1395_WES` basename. Source
`F1R2CountsCollector.java:188` passes unsorted `tmpDir.listFiles()` to
`IOUtils.writeTarGz`; `IOUtils.java:541–543` preserves the array order.
This is directory enumeration, not the genotyper RNG or an identity hash. No
within-filesystem nondeterminism or cross-filesystem order emulation is proved.

By member name in **both** local archives, ref_histogram's **43,305 bytes** and
alt_table's **134,856 bytes** equal the golden exactly. SHA256 respectively:
`5ee4eeff6f99db42285dc741331d3c554e17b3690fd3bc003fd13db4e1fd6107` and
`554aaa283ffd874337e997202fa7cf392dab6ab415c607a94f9a722fa5eb4cfe`.
All alt_histograms are 158,850 bytes; values by column label and non-column
lines match. Their byte hashes differ, as already expected from enum/pair-key
column ordering. The **member-order check is unchanged**. No fixed golden order
or additional exemption was implemented. This extra reference discrepancy is
recorded in `docs/REFERENCE-VARIABILITY.md` and reported with the partial port.

## Required verification commands

```sh
cargo build --release -p nomics > /tmp/m2-stage4-build-final.log 2>&1
NOMICS_M2_CALLS_ORACLE=/tmp/m2-stage4-final-oracle \
NOMICS_M2_LIKELIHOODS_ORACLE=/tmp/m2-likelihood-oracle/likelihood.likelihoods.tsv \
NOMICS_M2_ASSEMBLY_ORACLE=/tmp/nomics-m2-stage2-oracle/final-assembly.assembly.tsv \
NOMICS_M2_ORACLE_PREFIX=/tmp/nomics-m2-oracle/dump \
 RUST_TEST_THREADS=1 cargo test --release -p nomics-somatic -p nomics-io -p nomics-hc -- --nocapture \
 > /tmp/m2-stage4-tests.log 2>&1
cargo clippy --release --no-deps -p nomics-somatic -p nomics --all-targets -- -D warnings \
 > /tmp/m2-stage4-clippy-final.log 2>&1
bash scripts/somatic_mutect2_gate.sh
bash scripts/somatic_mutect2_regions_gate.sh > /tmp/m2-stage4-regions-gate.log 2>&1
NOMICS_M2_ASSEMBLY_ORACLE=/tmp/nomics-m2-stage2-oracle/assembly.assembly.tsv \
 bash scripts/somatic_mutect2_haplotypes_gate.sh /tmp/m2-stage4-haplotypes \
 > /tmp/m2-stage4-haplotypes-gate.log 2>&1
NOMICS_M2_ASSEMBLY_ORACLE=/tmp/nomics-m2-stage2-oracle/final-assembly.assembly.tsv \
 bash scripts/somatic_mutect2_haplotypes_gate.sh /tmp/m2-stage4-haplotypes chr21 \
 > /tmp/m2-stage4-haplotypes-corrected.log 2>&1
cmp /tmp/m2-stage4-haplotypes/chr21.assembly.tsv \
 /tmp/nomics-m2-stage2-oracle/final-assembly.assembly.tsv
NOMICS_M2_LIKELIHOODS_ORACLE=/tmp/m2-likelihood-oracle/likelihood.likelihoods.tsv \
 bash scripts/somatic_mutect2_likelihoods_gate.sh /tmp/m2-stage4-likelihoods \
 > /tmp/m2-stage4-likelihoods-gate.log 2>&1
NOMICS_RG_DS=<reference path of the golden run> \
 bash scripts/mac_tier0.sh cpu /tmp/m2-stage4-tier0 > /tmp/m2-stage4-tier0.log 2>&1
grep -nE '\.(powf|powi|ln|exp|log10|log2|ln_1p|exp_m1)\(' \
 crates/nomics-hc/src/call.rs crates/nomics-hc/src/trim.rs \
 crates/nomics-somatic/src/lib.rs crates/nomics-somatic/src/mutect2_assembly/mod.rs \
 crates/nomics-somatic/src/mutect2_likelihoods/mod.rs \
 crates/nomics-somatic/src/mutect2_calls/*.rs \
 crates/nomics-somatic/tests/mutect2_calls_goldens.rs crates/nomics/src/main.rs \
 > /tmp/m2-stage4-transcendentals.log
```

Build and clippy exit **0**; grep exits **1** with **zero hits**, as expected.
Existing HC Atomic::fetch_update deprecation and unused_mut warnings are left
unchanged. No new ignores, numerical tolerances or output masks were added.
The tier0 script first failed in the sandbox at sysctl with “Operation not
permitted”; rerunning the same authorized command with sandbox disabled passed:

```text
== cpu: wall 27s, peak RSS 14.5 GB
PASS recal_table count=1090
PASS vcf count=1047
verdict: PASS
```

All **41 digests** match. Shared HC changes only expose Event/MNP event extraction
and add the spanning-event parameterization around the existing merger/mapper.
Its existing call path still uses spanning events and maxMnpDistance zero. Its
Java-exact math is reused, not replaced; the somatic evidence model is separate.
Assembly exposes its original padded reference and stage three gains a
synchronous callback; original stage commands preserve their output contracts.

The required **final caller gate fails on all seven contigs**, exit **1**:

```text
FAIL chr21: caller failed
FAIL chr22: caller failed
FAIL chr19: caller failed
FAIL chr20: caller failed
FAIL chr13: caller failed
FAIL chr16: caller failed
FAIL chr12: caller failed
```

First error for every contig: `error: unrecognized subcommand 'mutect2'`.
There is no Rust VCF record to compare (0 compared out of 38/134/335/144/97/449/337,
not 0 matching records), no generated archive and **zero XA tags compared**.
No BGZF VCF or tabix index is written. Stats comparisons below exercise the
existing traversal through the diagnostic, not a completed final caller.

Local resource note: an initial concurrent full-suite attempt was intentionally
terminated (SIGTERM to its own test child) during the new integration test to
relieve memory/swap pressure. It exited 101; it is not counted as a test pass
or as a numerical mismatch. Its log is `/tmp/m2-stage4-tests-interrupted.log`.
The full command above was subsequently rerun after the other gates. Large
already-compared temporary files were pruned with SHA256 values retained in
`/tmp/m2-stage4-pruned-files.sha256`; qualification dumps and all logs were kept.
Contig wall times below were measured during these local checks, with competing
processes for part of the run; they are not saturated-path benchmarks.

Assembly fixture correction: the first all-seven stage-two command above passed
all seven artificial-record comparisons, then exited **1** solely because its
optional chr21 internal comparison used the obsolete `assembly.assembly.tsv`.
The first difference was byte 1542, line 3: the obsolete `D` row has 8 fields,
whereas the current dump has 10, with BI/BD appended. The earlier stage-two
record already identifies `final-assembly.assembly.tsv` as the final qualified
BI/BD fixture; using that file, `cmp` exits **0**, no output. The corrected
chr21 script run exits **0**, with `PASS chr21 artificial` and `PASS chr21 internal`
(wall 9.46 s). No source or expected values were changed to repair this path
mistake. The full-suite rerun uses the final fixture.

## Completed pre-annotation count/stats and first-stage results

Standalone integration command:

```sh
cargo test --release -p nomics-somatic --test mutect2_calls_goldens -- --nocapture \
 > /tmp/m2-stage4-integration.log 2>&1
```

It exits **0**, one passed test, **958.41 s** including the seven runs. No golden
directory is skipped. In this standalone run the optional internal dump paths
were unset and explicitly skipped; the explicit chr21 cmp above and the later
full-suite run supply the qualified dump. Passing-event counts are not VCF
record comparisons; no archive or XA tag is compared by this test.

| Contig | Varying regions | Events scored | Passing-event count / golden VCF count | Stats bytes, exact | Diagnostic wall s |
|---|---:|---:|---:|---:|---:|
| chr21 | 364 | 862 | 38 / 38 | 35 | 20.917 |
| chr22 | 702 | 1,984 | 134 / 134 | 35 | 64.195 |
| chr19 | 1,909 | 5,522 | 335 / 335 | 35 | 244.348 |
| chr20 | 801 | 2,294 | 144 / 144 | 35 | 122.449 |
| chr13 | 596 | 1,599 | 97 / 97 | 35 | 43.818 |
| chr16 | 1,457 | 4,110 | 449 / 449 | 35 | 303.951 |
| chr12 | 1,858 | 5,373 | 337 / 337 | 35 | 157.708 |

First-stage gate exits **0**, `7 pass, 0 fail, 0 skipped`. Its complete summary:

```text
PASS chr21 (7167 regions, 699 active; 1230643.0 callable; 1025272 total reads filtered out of 3023196 reads processed; 8 s)
PASS chr22 (13815 regions, 1170 active; 2689447.0 callable; 2445097 total reads filtered out of 8935845 reads processed; 29 s)
PASS chr19 (34338 regions, 3386 active; 6845249.0 callable; 5481759 total reads filtered out of 20756074 reads processed; 232 s)
PASS chr20 (15433 regions, 1367 active; 3161042.0 callable; 2283986 total reads filtered out of 8908067 reads processed; 29 s)
PASS chr13 (11531 regions, 1166 active; 2157790.0 callable; 1458591 total reads filtered out of 5253705 reads processed; 15 s)
PASS chr16 (26083 regions, 2352 active; 4859121.0 callable; 5195319 total reads filtered out of 15802761 reads processed; 222 s)
PASS chr12 (35375 regions, 3551 active; 6840733.0 callable; 4524389 total reads filtered out of 16361829 reads processed; 71 s)
mutect2 regions gate: 7 pass, 0 fail, 0 skipped
```

Second-stage external artificial records are exact on all seven: **3,160**,
**6,883**, **19,219**, **7,802**, **5,244**, **14,892**, **19,179**, in contig order
as above. CLI wall seconds: 10.15, 32.78, 354.51, 37.50, 20.19, 71.23, 76.18.
The optional internal fixture path correction and its successful recheck are
recorded above rather than hidden as a passing initial script invocation.

## Source corrections, hash order and remaining limits

Confirmed refinements/corrections to the prose spec, now reflected in §14 and
in the relevant §§5–6: (1) groups of at most two reads are not flag-filtered;
larger groups use first clean survivors or the first original read if none
survive; (2) createAlleleMapper can append `*` even when the merged VC excludes
spanning events, and that row affects evidence/LODs; (3) the Q45 second cap is
`45*(-Math.log(10)/10)` in that operation order; (4) sumArrayFunction starts from
column zero, not a zero vector; (5) logSumExp still adds positive zero when its
scaled sum equals one, retaining the signed-zero behavior omitted in the
former wording; (6) the archive order is unsorted directory enumeration.
No numerical constants were adjusted to match the goldens.

The genotyper RNG is unseeded, but its only nextDouble call is short-circuited
by default-false genotypeGermlineSites. No unseeded value was observed to reach
an output in this configuration. The oracle rejects that non-default switch.

String name-map order is the existing OpenJDK HashMap model (computeIfAbsent
head insertion, resize splitting and tree-root movement), now independently
exercised for fragment grouping. Group members retain original read order.
Header samples are sorted; event starts are sorted; Event HashSet is used for
membership rather than allele ordering; merged alleles retain first-seen
order; allele mapping is LinkedHashMap order, with an appended `*` last.
Haplotype hashes are Arrays.hashCode(bases), SimpleAllele hashes use bases and
reference status, and Event hashes use start and alleles, rather than JVM
object identity. The chr21 raw dump proves the exercised orders. It is not an
arbitrary-graph identity-order proof for the prior assembly stage, nor proof
of every later HashMap/HashSet path in annotations/phasing. Existing assembly
identity-order limitations remain unchanged. F1R2 enum/pair-key column order
is already known to vary and has only the requested column-label comparison
proposal. The newly observed directory/member order is not normalized.

Incomplete/unconfirmed in Rust: population AF lookup/POPAF annotation, initial
genotype construction and AF/AD emission, call trimming/remapping, all standard
Mutect annotations and their independent strand/other-statistics fixtures,
physical phasing, ECNT/ECNTH, genotype-dependent XA tags, VCF header/formatting,
BGZF VCF/index writing, f1r2 collection and archive writing. The Java annotation
dump has 38 calls and 41,524 bytes; no Rust annotation dump is compared. Thus
the required full internal gate is **incomplete**, despite the exact
pre-annotation subset. No VCF/archive output is borrowed from the reference.
Only the inherited two-indexed-BAM, one-normal/one-tumour BED/default contract
is exposed for the diagnostic. Tumour-only, multiple tumour/normal samples,
forced alleles, germline-site genotyping, mitochondrial/flow/gVCF modes,
non-default minAF/likelihood/quality/filter/region settings, independent mates,
CRAM and SAM are not exposed/supported. No GPU or speedup claim is made.

Third-stage script exits **0**, including its qualified chr21 internal dump.
The biological projection still excludes XA and therefore is **not** the
requested completed read-tag gate. All remaining projected fields/tags match:

| Contig | Biological records exact / total | Artificial records exact / total | Stage-three wall s |
|---|---:|---:|---:|
| chr21 | 95,695 / 95,695 | 3,160 / 3,160 | 22.04 |
| chr22 | 262,031 / 262,031 | 6,883 / 6,883 | 102.74 |
| chr19 | 563,178 / 563,178 | 19,219 / 19,219 | 327.00 |
| chr20 | 250,879 / 250,879 | 7,802 / 7,802 | 63.92 |
| chr13 | 165,350 / 165,350 | 5,244 / 5,244 | 37.53 |
| chr16 | 526,456 / 526,456 | 14,892 / 14,892 | 115.82 |
| chr12 | 537,627 / 537,627 | 19,179 / 19,179 | 135.67 |

All seven external assembly/biological boundaries pass with their original
normalizations unchanged. Complete stdout and wall times are retained in the
stage-two/stage-three gate logs named above. XA tags compared: **0 on every
contig**. The new full-caller gate remains a failure on every contig.

## Final full-suite rerun

The full release test command recorded above completed with exit **0**: **165 passed, zero failed, 10 existing ignored** across 32 test binaries/doc-test summaries. No new ignored tests were added. The two new independent numerical/fragment unit tests execute and pass. The new integration test runs all seven available goldens and the qualified chr21 internal comparison (`PASS chr21: internal pre-annotation dump bytes`); its other six optional internal dumps are absent and explicitly skipped, not claimed. The existing assembly, likelihood and region integration tests also pass with their qualified internal dumps supplied. All 41 independent tier0 germline digests remain exact.

Seven-contig pre-annotation wall times in this serialized full-suite run (same exact passing counts and 35-byte stats as the standalone table):

| Contig | Wall s |
|---|---:|
| chr21 | 21.413 |
| chr22 | 65.417 |
| chr19 | 156.067 |
| chr20 | 63.013 |
| chr13 | 36.803 |
| chr16 | 118.536 |
| chr12 | 140.831 |

Its new integration test took **602.29 s**. These are verification wall times, not a benchmark or speedup. The final raw-VCF/archive/XA implementation and full unrounded annotation gate remain **incomplete**; the final CLI gate still fails on all seven contigs. None of the passing regression checks changes that status.

# Continuation: read/fragment annotations and physical phasing (partial)

Branch `somatic/mutect2-vcf`, based on `somatic/mutect2-calls` at 81658dc.
**This continuation is incomplete.** There is still no final `somatic mutect2`
command, raw VCF writer, f1r2 collector/archive writer, or XA implementation.
No output is borrowed from a golden/reference tool. No push, PR or merge was performed.

The existing diagnostic now adds three sidecars to its --output path:

* `.annotations.tsv`: raw subset, INFO AS_SB_TABLE/DP/MBQ/MFRL/MMQ/MPOS and
  genotype AD/DP/AF/F1R2/F2R1/FAD/SB. AF is raw binary64 bits.
* `.annotation-fields.tsv`: the same subset formatted by the existing
  Java-exact HC VCF number formatter, for comparisons needing no Java.
* `.phasing.tsv`: final GT/PGT/PID/PS per emitted call and sample.

These are diagnostic projections, **not VCFs**. The integration test compares
all represented values and coordinates/alleles/sample order without tolerances.
The final VCF gate and full internal post-annotation gate remain incomplete.
Existing pre-annotation dump bytes remain an independent comparison.

## Source audit and independent fixtures

Source checkouts outside the repo were verified at GATK
`76edc75c26504da94bbaee66584e107e76ee15de` and htsjdk
`4cc010022ac038fb30f26e6f9717fabff3e808c1`. The qualified unmodified-jar
chr21 pre/post-annotation dumps from the prior worker were reused; this
continuation did not requalify another JVM or replace that evidence.

`MutectAnnotationNumerics.java` generates 89 matrix cases (80 from the earlier
independent evidence fixture, five empty-evidence cases, a tie-boundary matrix,
and 13/20/49-allele resize cases), nine median cases, and 119 CIGAR/coordinate
cases. The wrapper invokes the actual private getEffectiveCounts/addGenotypes
methods, and public best-allele, coordinate and median helpers. Raw output bits
are committed in `mutect2_calls/annotations.tsv`. `MutectPhasing.java` invokes
the private constructPhaseSetMapping on 120 deterministic set graphs; output
group IDs and slots are in `mutect2_calls/phasing.tsv`. These test oracles are
never invoked or linked by nomics. Regeneration commands:

```sh
jar=$HOME/omics/oracle/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar
jdk=/opt/homebrew/opt/openjdk@17/bin
out=/tmp/m2-stage4-final-oracle
"$jdk/javac" -cp "$jar" -d "$out" reference/oracle/mutect2/MutectAnnotationNumerics.java reference/oracle/mutect2/MutectPhasing.java
"$jdk/java" -cp "$out:$jar" MutectAnnotationNumerics crates/nomics-somatic/src/mutect2_calls/evidence.tsv > /tmp/m2-vcf-annotation-fixture.tsv
cmp /tmp/m2-vcf-annotation-fixture.tsv crates/nomics-somatic/src/mutect2_calls/annotations.tsv
"$jdk/java" -cp "$out:$jar" MutectPhasing > /tmp/m2-vcf-phasing-fixture.tsv
cmp /tmp/m2-vcf-phasing-fixture.tsv crates/nomics-somatic/src/mutect2_calls/phasing.tsv
```

Source refinements, not normalizations: matrix-row order after HashMap
marginalization is independent of emitted allele order; the first maximum and
reference-priority ln tie window precede the literal .2 informative test;
ReadPosition's outer guard rejects a leading-insertion position before the
rank-sum helper can supply zero; deletion read coordinates use the first base
after the deletion, despite a ReadUtils overload's prose saying "preceding".
All are described in spec §15. AS_SB_TABLE has a strict per-sample >2 threshold;
orientation counts require BQ>=20 and the first selected fragment read.
No formula or constant was tuned to a golden.

The first chr21 serialized-field comparison failed at record 5, coordinate
chr21:10716938: the diagnostic had emitted `ATCATTTCCATTCCC,A,*` whereas the
VCF has `ATCATTTCCATTCCC,A`. The existing merged/matrix API appends the spanning
row, which must remain in evidence but is not a merged-VC emission candidate.
Filtering it only when selecting emitted rows fixed that implementation error.
The next chr21 comparison passed all represented fields and raw AF bits.
An initial tiny-fixture attempt also exposed a test-wrapper error: leaving
AlleleLikelihoods at its constructor's log10 scale tested the wrong tie window.
The generator now switches the scale flag **before** setting ln values.
The corrected independent bits pass; the Rust constants were unchanged.

Ordering audit: samples stay sorted header order; alleles stay merged first-seen
order with threshold filtering; annotation matrices use modeled SimpleAllele
HashMap list-bin order, whose hash is bases/reference status rather than object
identity. Unsupported tree-bin order is an explicit error. No identity-hash
value was observed to reach these tested fields. Fragment String HashMap order
and member order remain the already-gated implementation. Phase connection
logic uses call order plus set membership/size/containment; no set traversal
chooses a phase. Existing assembly identity-order limits are unchanged. The
only unseeded genotyper nextDouble remains short-circuited by the default-false
genotypeGermlineSites flag. No RNG was added.

## Archive-order proposal and reader check

The continuation's requested normalization changes the proposal to exact
member-name **sets**, independently retaining byte comparison for ref_histogram
and alt_table and label-keyed comparison for alt_histogram. Duplicate member
names remain errors. Header normalization is still only GATKCommandLine.
No new VCF exemption or numerical tolerance was introduced.

Source LearnReadOrientationModel.doWork extracts each archive separately and
independently discovers the three member types by suffix. With one tumour and
one member per type, tar order cannot alter the within-type lists. Input
archive order remains significant and is not normalized. Intended fixed writer
order for a future collector is ref_histogram, alt_histogram, alt_table; **no
Rust archive is written by this continuation**.

Reader check: the nomics learner was run on the chr21 golden archive and a tar
permutation containing exactly the same member bytes. Their prior-member sets
and payload bytes are identical: `./HCC1395_HCC1395_WES.orientation_priors`,
15,362 bytes, SHA256
`67087bb44890f3f59c6deeb4ecb4b5035e7a588a60ad8a136db9026bb7924e38`.
Artifacts and logs: `/tmp/m2-vcf-orientation/`. This is an archive-order reader
check using reference counts, **not** the requested golden-prior check on a
Rust-collected archive, which cannot yet run.

The final-gate comparator was independently exercised using copied reference
VCF/stats in a test directory, first permuting archive members, then swapping
two alt-histogram columns together with their labels: both comparisons passed.
Changing one histogram count was rejected with
`.alt_histogram by column label: first difference at 1`. This is comparator
validation only, not Rust VCF/archive parity. Log: `/tmp/m2-vcf-comparator.log`.

## Seven-contig annotation/phasing results

The final-source integration test in the full release suite passed, **595.21 s**:

```text
PASS chr21: six INFO annotation fields and seven genotype fields on 38 calls
PASS chr21: GT/PGT/PID/PS on 38 calls (all samples)
PASS chr21: internal pre-annotation dump bytes
PASS chr21: raw annotation subset bytes, AF binary64 bits (not the full internal gate)
PASS chr22: annotation subset and GT/PGT/PID/PS on 134 calls
PASS chr19: annotation subset and GT/PGT/PID/PS on 335 calls
PASS chr20: annotation subset and GT/PGT/PID/PS on 144 calls
PASS chr13: annotation subset and GT/PGT/PID/PS on 97 calls
PASS chr16: annotation subset and GT/PGT/PID/PS on 449 calls
PASS chr12: annotation subset and GT/PGT/PID/PS on 337 calls
```

No contig was skipped. The other six raw internal dumps are absent and are
explicitly reported as skipped; only chr21 has a qualified raw-bit comparison.

| Contig | Calls with every represented field exact / total | Stats bytes exact | Diagnostic wall s |
|---|---:|---:|---:|
| chr21 | 38 / 38 | 35 | 20.717 |
| chr22 | 134 / 134 | 35 | 62.608 |
| chr19 | 335 / 335 | 35 | 154.429 |
| chr20 | 144 / 144 | 35 | 62.003 |
| chr13 | 97 / 97 | 35 | 36.377 |
| chr16 | 449 / 449 | 35 | 117.756 |
| chr12 | 337 / 337 | 35 | 141.009 |

These are annotation-subset/genotype comparisons, **not entire VCF records**.
Stats are written by the diagnostic traversal, not an unavailable final caller.
Walls are verification times, not saturated-path benchmarks or speedups.

A direct full post-annotation comparison fails at byte **0**: Rust's subset
starts `V chr21:9248963 G,A`; the jar's complete dump starts
`R chr21:5039034-5039053`. Lengths are 24,681 versus 41,524 bytes. The complete
region/call structure, remaining INFO fields and combined genotype dump are
absent. No normalization is proposed to make this full comparison pass.
Log: `/tmp/m2-vcf-full-post-dump.log`.

The separate standalone annotation-only integration run (before the phasing
sidecar was added) also passed all seven, 604.44 s. Its intermediate chr21
failure and correction are recorded above; it is not used as evidence of the
later physical-phasing code. The focused final-source chr21 phase run passed
in 20.88 s, explicitly skipping the six other directories in its one-contig
fixture root; the full-suite run above supplies all seven.

## Required build/test/CLI checks

Commands (all from this worktree; default seven-contig lists):

```sh
cargo build --release -p nomics
NOMICS_M2_CALLS_ORACLE=/tmp/m2-stage4-final-oracle \
NOMICS_M2_LIKELIHOODS_ORACLE=/tmp/m2-likelihood-oracle/likelihood.likelihoods.tsv \
NOMICS_M2_ASSEMBLY_ORACLE=/tmp/nomics-m2-stage2-oracle/final-assembly.assembly.tsv \
NOMICS_M2_ORACLE_PREFIX=/tmp/nomics-m2-oracle/dump \
RUST_TEST_THREADS=1 cargo test --release -p nomics-somatic -p nomics-io -p nomics-hc -- --nocapture
cargo clippy --release --no-deps -p nomics-somatic -p nomics --all-targets -- -D warnings
bash scripts/somatic_mutect2_gate.sh
bash scripts/somatic_mutect2_regions_gate.sh
NOMICS_M2_ASSEMBLY_ORACLE=/tmp/nomics-m2-stage2-oracle/final-assembly.assembly.tsv \
 bash scripts/somatic_mutect2_haplotypes_gate.sh /tmp/m2-vcf-haplotypes
NOMICS_M2_LIKELIHOODS_ORACLE=/tmp/m2-likelihood-oracle/likelihood.likelihoods.tsv \
 bash scripts/somatic_mutect2_likelihoods_gate.sh /tmp/m2-vcf-likelihoods
```

| Check | Exit | Wall s |
|---|---:|---:|
| build | 0 | 21.005 |
| tests | 0 | 1116.457 |
| clippy | 0 | 1.275 |
| final-gate | 1 | 0.050 |
| regions-gate | 0 | 242.515 |
| haplotypes-gate | 0 | 294.422 |
| likelihoods-gate | 0 | 559.753 |

Build and clippy pass. Full release suite: **167 passed, zero failed, 10 existing
ignored**, 32 summaries including doc tests. No new ignored test. The existing
assembly integration test covers three contigs, the existing likelihood test
five; the **CLI scripts** below supply seven-contig checks. New annotation and
phasing integration checks cover all seven. Existing HC deprecation/unused_mut
warnings remain unchanged. Full logs are `/tmp/m2-vcf-{build,tests,clippy}.log`.

Region CLI stdout (`/tmp/m2-vcf-regions-gate.log`):

```text
PASS chr21  (7167 regions, 699 active; 1230643.0 callable; 1025272 total reads filtered out of 3023196 reads processed; 7 s)
PASS chr22  (13815 regions, 1170 active; 2689447.0 callable; 2445097 total reads filtered out of 8935845 reads processed; 21 s)
PASS chr19  (34338 regions, 3386 active; 6845249.0 callable; 5481759 total reads filtered out of 20756074 reads processed; 71 s)
PASS chr20  (15433 regions, 1367 active; 3161042.0 callable; 2283986 total reads filtered out of 8908067 reads processed; 23 s)
PASS chr13  (11531 regions, 1166 active; 2157790.0 callable; 1458591 total reads filtered out of 5253705 reads processed; 14 s)
PASS chr16  (26083 regions, 2352 active; 4859121.0 callable; 5195319 total reads filtered out of 15802761 reads processed; 44 s)
PASS chr12  (35375 regions, 3551 active; 6840733.0 callable; 4524389 total reads filtered out of 16361829 reads processed; 62 s)
mutect2 regions gate: 7 pass, 0 fail, 0 skipped
```

Assembly and likelihood CLI projections all pass, with **no XA comparison**.
Qualified chr21 internal comparisons both print `PASS chr21 internal`.
Logs: `/tmp/m2-vcf-{haplotypes,likelihoods}-gate.log`.

| Contig | Artificial exact / total | Stage-two wall s | Biological exact / total, excluding XA | Stage-three wall s |
|---|---:|---:|---:|---:|
| chr21 | 3,160 / 3,160 | 9.21 | 95,695 / 95,695 | 18.64 |
| chr22 | 6,883 / 6,883 | 27.96 | 262,031 / 262,031 | 55.09 |
| chr19 | 19,219 / 19,219 | 84.11 | 563,178 / 563,178 | 143.64 |
| chr20 | 7,802 / 7,802 | 28.40 | 250,879 / 250,879 | 56.45 |
| chr13 | 5,244 / 5,244 | 17.16 | 165,350 / 165,350 | 33.62 |
| chr16 | 14,892 / 14,892 | 54.24 | 526,456 / 526,456 | 108.80 |
| chr12 | 19,179 / 19,179 | 68.65 | 537,627 / 537,627 | 131.95 |

After successful comparisons, only this run's bulky assembly/likelihood dumps
were pruned to conserve local disk; SHA256, size and paths are retained in
`/tmp/m2-vcf-pruned.sha256`. Golden projections, read projections, small
annotation sidecars (where retained), logs and qualified independent fixtures
were not replaced. These CLI walls are verification timings, not benchmarks.

The **final caller gate fails on all seven**, exit **1**. Actual stdout:

```text
Outputs: /var/folders/26/rpv6m28d2_g1xf44f1f58h5r0000gp/T//nomics-m2-calls.qsLR67
FAIL chr21: caller failed (see /var/folders/26/rpv6m28d2_g1xf44f1f58h5r0000gp/T//nomics-m2-calls.qsLR67/chr21.log)
FAIL chr22: caller failed (see /var/folders/26/rpv6m28d2_g1xf44f1f58h5r0000gp/T//nomics-m2-calls.qsLR67/chr22.log)
FAIL chr19: caller failed (see /var/folders/26/rpv6m28d2_g1xf44f1f58h5r0000gp/T//nomics-m2-calls.qsLR67/chr19.log)
FAIL chr20: caller failed (see /var/folders/26/rpv6m28d2_g1xf44f1f58h5r0000gp/T//nomics-m2-calls.qsLR67/chr20.log)
FAIL chr13: caller failed (see /var/folders/26/rpv6m28d2_g1xf44f1f58h5r0000gp/T//nomics-m2-calls.qsLR67/chr13.log)
FAIL chr16: caller failed (see /var/folders/26/rpv6m28d2_g1xf44f1f58h5r0000gp/T//nomics-m2-calls.qsLR67/chr16.log)
FAIL chr12: caller failed (see /var/folders/26/rpv6m28d2_g1xf44f1f58h5r0000gp/T//nomics-m2-calls.qsLR67/chr12.log)
```

Every per-contig first error is `error: unrecognized subcommand 'mutect2'`.
No port VCF record exists to compare. The following are **uncompared**, not
zero matching records:

| Contig | Whole raw VCF records compared / total | Header | Final stats file | Rust f1r2 members compared / 3 | Biological records compared with XA / total |
|---|---:|---|---|---:|---:|
| chr21 | 0 / 38 | not compared | absent | 0 / 3 | 0 / 95,695 |
| chr22 | 0 / 134 | not compared | absent | 0 / 3 | 0 / 262,031 |
| chr19 | 0 / 335 | not compared | absent | 0 / 3 | 0 / 563,178 |
| chr20 | 0 / 144 | not compared | absent | 0 / 3 | 0 / 250,879 |
| chr13 | 0 / 97 | not compared | absent | 0 / 3 | 0 / 165,350 |
| chr16 | 0 / 449 | not compared | absent | 0 / 3 | 0 / 526,456 |
| chr12 | 0 / 337 | not compared | absent | 0 / 3 | 0 / 537,627 |

Diagnostic traversal stats are still exact (35 bytes each), as reported above;
this does not create a final `.vcf.gz.stats`. No BGZF VCF or `.tbi` is written.
No Rust-collected archive exists, so its orientation-prior gate cannot run.
The reference-count reader/permutation check above is reported separately.

Forbidden-function scan of **all changed/added Rust and Java source files**:
`grep -rnE "\.(powf|powi|ln|exp|log10|log2|ln_1p|exp_m1)\(" ...` exits **1**,
zero hits. Log `/tmp/m2-vcf-source-transcendentals.log` is empty. Extending the
scan literally to every changed text file produces two **pre-existing** Java
math expressions in `docs/spec/mutect2.md`, lines 242 and 1007; these are source
specification prose, not new calls. They were not removed to hide grep hits.
That broader scan is `/tmp/m2-vcf-transcendentals.log`.

Still not implemented/confirmed: final CLI and VCF header/record/BGZF/index
serialization, POPAF/PON, stored INFO LOD arrays, ECNT/ECNTH, RPA/RU/STR,
f1r2 collection and archive output, XA, and the complete raw post-annotation
boundary. Other options and assembly identity-order limitations remain as in
spec §§11–15. Tree-bin annotation allele-map order explicitly errors rather
than guessing an unproved order. No GPU path or speedup claim is made.

## HaplotypeCaller unchanged: tier0 CPU check

No files under `crates/nomics-hc` differ from base 81658dc:
`git diff --exit-code somatic/mutect2-calls -- crates/nomics-hc` exits **0**.
Its release tests pass as part of the full suite above. The requested tier0
command, with the existing workspace build directory reused, exits **0**:

```sh
NOMICS_BUILD_DIR=$PWD/target \
NOMICS_RG_DS=<reference path of the golden run> \
 bash scripts/mac_tier0.sh cpu /tmp/m2-vcf-tier0
```

```text
== cpu: tier0 FASTQ -> VCF, 14 threads -> /tmp/m2-vcf-tier0
== cpu: wall 17s, peak RSS 15.0 GB
PASS recal_table count=1090
PASS vcf count=1047
verdict: PASS
```

All **41 digests** match the independent tier0 manifest. Logs/artifacts:
`/tmp/m2-vcf-tier0.log`, `/tmp/m2-vcf-tier0/digest-compare.txt`,
`/tmp/m2-vcf-tier0/digests.tsv`. The command ran with sandbox disabled because
sysctl and existing cache access lie outside the sandbox, as established by
the prior verification. No approval rejection occurred. This is a regression
check, not a speedup benchmark.

No push, pull request or merge was performed.
The final caller remains incomplete even though the implemented annotation
subset, phasing, earlier stage gates and germline regression pass. A successor
must finish the remaining outputs and full raw annotation gate before claiming
Mutect2 parity; no final-output gate failure has been converted to a skip.


## Complete raw-VCF record gate continuation (2026-11-13)

Implemented `nomics somatic mutect2` on base 1ba2d7b: complete default raw VCF
schema, resource annotations, stored LODs, ECNT/ECNTH, tandem repeats, BGZF and
caller stats. The stage-three projection now includes genotype-dependent XA,
including its existing input prefix. No reference tool is used by nomics.

| Contig | Records identical | Header* | Stats bytes | Biological records incl. XA | Caller wall s |
|---|---:|---|---|---:|---:|
| chr21 | 38/38 | PASS | PASS (35) | 95695/95695 | 19.17 |
| chr22 | 134/134 | PASS | PASS (35) | 262031/262031 | 61.91 |
| chr19 | 335/335 | PASS | PASS (35) | 563178/563178 | 192.57 |
| chr20 | 144/144 | PASS | PASS (35) | 250879/250879 | 60.73 |
| chr13 | 97/97 | PASS | PASS (35) | 165350/165350 | 35.92 |
| chr16 | 449/449 | PASS | PASS (35) | 526456/526456 | 117.67 |
| chr12 | 337/337 | PASS | PASS (35) | 537627/537627 | 139.51 |

*Header excludes only `##GATKCommandLine`. No record normalization. All 1,534
records and 2,401,216 biological records match. No differing-field rows remain.
Full final gate stdout is [records-gate-20261113.txt](records-gate-20261113.txt).
The gate preserves record line endings and reports per-key differences and
header order. BGZF framing/EOF was checked; no `.tbi` is generated.

The final caller walls above were measured in a serial run under
`lockf -k /tmp/nomics-heavy.lock`. They are verification times, not speedup or
saturated-CPU benchmarks. `/usr/bin/time -l` on **chr12**, the longest golden
contig, reported peak RSS **6,518,816,768 bytes (6.071 GiB)**. The timed
caller's VCF is itself the chr12 output passing the final gate. After the lead's
memory constraint, existing commands were allowed to finish, the outer stage
runner was stopped before launching another unlocked gate, and remaining heavy
commands ran sequentially under the shared lock. The wrapper also waited for
the previously observed other-worker unlocked test to exit.

Verification:

* `cargo build --release -p nomics`: PASS.
* `RUST_TEST_THREADS=1 cargo test --release -p nomics-somatic -p nomics-io -p nomics-hc`:
  167 passed, zero failed, 10 existing ignored, 32 summaries. The two subsequently
  added focused VCF tests also pass in release mode (mixed-site eligibility and
  missing versus malformed AF). No new ignored tests. The golden tests now
  remove each contig's artifacts only after its checks pass, retaining failures.
* `cargo clippy --release --no-deps -p nomics-somatic -p nomics --all-targets -- -D warnings`:
  PASS; inherited HC dependency warnings are unchanged.
* `scripts/somatic_mutect2_gate.sh`: all seven PASS, final source, under the lock.
* Region gate: 7 pass, 0 fail, 0 skipped. Assembly and likelihood CLI gates:
  all seven pass, with XA included. Their initially skipped optional chr21
  internal checks were resolved by fresh qualified jar dumps: assembly and
  likelihood dumps compare byte-for-byte, and the region integration test
  passes with `NOMICS_M2_ORACLE_PREFIX` set to the regenerated dump.
* HC source has no diff from `somatic/mutect2-vcf`. Requested tier0 CPU command
  with the recorded `NOMICS_RG_DS` ends `verdict: PASS`: all 41 digests match.
* Requested transcendental scan: zero changed Rust/Java source hits. Scanning
  every changed text file also finds two unchanged mathematical formulas in
  the spec (FastMath.exp and Math.log10); neither is executable port code.
* `--f1r2-tar-gz` returns `not yet wired` before touching inputs or outputs.

Source corrections/refinements are in spec §16: start-only resource queries;
no marginalization-only `*` in ECNT; actual reference-haplotype priority and
stable support-count ties for ECNTH; MIXED sites excluded from TandemRepeat;
signed repeat counts; XA appends rather than replaces an existing tag; AF `.`
defaults while malformed numbers fail. Common-REF extension follows the source's
length rule, including its lack of an extra prefix check for unequal lengths.
The first chr21 run exposed ECNT/ECNTH/POPAF differences, and chr13 exposed two
mixed-site repeat annotations. These were fixed from source, not normalized.

Ordering: emitted ALT order and sorted sample order are preserved; GT precedes
lexically sorted FORMAT keys and INFO keys sort lexically. Resource allele maps
and event/support maps supply lookup/membership, not iteration-derived output.
Stable haplotype-list order breaks support ties. Existing OpenJDK String-key
fragment map emulation and SimpleAllele list-bin annotation order remain in
force. Annotation tree bins (at least nine colliding keys at capacity >=64)
explicitly fail; equal-hash non-comparable keys can require JVM identity order.
No golden hits this: maximum total alleles is five. Inherited assembly identity
limits remain documented; no new identity-derived output order is guessed.

Unsupported/unconfirmed: f1r2 collection wiring, index output, non-default modes
listed in §§11–15, arbitrary assembly identity-order cases and native/GPU paths.
The complete internal post-annotation raw-bit boundary is not claimed by these
record gates. No new output exemption, push, PR or merge.

Logs and outputs are retained in this worktree's `target/m2-verification/`;
qualified Java dumps are `/tmp/m2-records-oracle-{assembly,likelihood}`.

The raw chr12 time report is [records-rss-chr12-20261113.txt](records-rss-chr12-20261113.txt).
It additionally reports 10,559,381,424 bytes of peak memory footprint (distinct
from RSS). Java 17's log-conversion constant and Rust's `LOG10_E` were checked
independently: both have binary64 bits `3fdbcb7b1526e50e`.

# F1R2 collector implementation (2026-11-13)

Branch `somatic/mutect2-f1r2`, based on local `somatic/mutect2-vcf` at
`1ba2d7b0bdf699c4755b88ae5be42f738de395fe`. No delegates, push, PR or merge.
Caller/VCF modules and stack map are unchanged. Authority: the source fetched
outside the repo at `/tmp/mutect2-spec-gatk`, verified tag 4.6.2.0 / HEAD
`76edc75c26504da94bbaee66584e107e76ee15de`, and htsjdk 4.2.0. Audited the
collector and its histogram/table/filter helpers, Mutect2Engine.isActive,
ReadPileup/PileupElement, ReadUtils/GATKRead, MathUtils.median and htsjdk
MetricsFile/FormatUtil. No reference code is linked or executed by nomics.
The prior JDK17 qualification above remains applicable. No collector mismatch
required a new Java run or internal per-locus dump.

Final-caller integration (one BAM traversal):
`Collector::default()` -> `mutect2_regions::run_with_sinks(inputs, region_sink,
&mut collector)` -> `collector.write_archive(output)`. The standalone CLI uses
that same path. The hook receives every nonempty joint activity pileup in the
interval union before LOD/normal/PoN rejection. It selects the tumour, reuses
`activity::covered`/`element` for LIBS adapter/CIGAR rules, and consumes the exact
shared read arena and seeded-downsampled IDs after the twelve filters and
palindrome clipping. Header-only sample discovery handles empty tumour pileups.

Source refinements to spec §8.2: "before region mate-quality adjustment" means
the current region; earlier callRegion edits can already be present in later
pileups. Deletion Q16 is removed by the **signed-byte** BQ >20 filter before
the indel heuristic. Retained before-deletion / **after-insertion** elements
still count against integer depth/100. Non-ACGT elements enter median MQ but
not depth; even-count median averages the middle two (an exact integer
threshold comparison suffices). Ties choose first A/C/G/T ALT; table depth is
ref+dominant ALT, histogram depth includes minor ALTs. First-of-pair is false
for unpaired reads, regardless of raw flag64. '*' counts as A for depth but not
orientation. Reference histogram String-key HashMap order is reproduced from
64 lexicographic insertions and 128 buckets; integer counts serialize exactly
as FormatUtil's ungrouped integer-valued doubles. No transcendentals are needed.

Outside the new module: lib.rs export, main.rs CLI/dispatch; traversal gains
LocusSink/run_with_sinks and visibility for its existing pileup/sample helpers;
orientation exposes kmer/hash/url encoding and its existing tar read/write
transport (no algorithm changes). The learner's reduced Record(n,k,x,alt)
cannot represent ref_f1r2/original context, so a lossless AltSiteRecord is used.
Histograms/records use small allocations; serialization reuses the orientation
tar/native-deflate path. No dependencies, unsafe or numerical libraries added.
Fixed member order: **ref_histogram, alt_histogram, alt_table**. The two order
normalizations are proposals awaiting the lead in docs/GATES.md. The member
finding already exists in docs/REFERENCE-VARIABILITY.md and was not duplicated.

Commands (logs under `/tmp/m2-f1r2-*`):

```sh
cargo build --release -p nomics > /tmp/m2-f1r2-build.log 2>&1
cargo test --release -p nomics-somatic -p nomics-io > /tmp/m2-f1r2-tests.log 2>&1
cargo clippy --release --no-deps -p nomics-somatic -p nomics --all-targets -- -D warnings > /tmp/m2-f1r2-clippy.log 2>&1
M2_F1R2_OUT=/tmp/nomics-m2-f1r2-gate \
  bash scripts/somatic_mutect2_f1r2_gate.sh > /tmp/m2-f1r2-gates.log 2>&1
bash scripts/somatic_mutect2_regions_gate.sh > /tmp/m2-f1r2-regions-gates.log 2>&1
mkdir -p /tmp/nomics-m2-f1r2-priors
target/release/nomics somatic learn-read-orientation-model \
  --input /tmp/nomics-m2-f1r2-gate/chr21.f1r2.tar.gz \
  --output /tmp/nomics-m2-f1r2-priors/port.tar.gz
target/release/nomics somatic learn-read-orientation-model \
  --input "$HOME/omics/runs/somatic/m2-chr21/oracle/HCC1395_WES_vs_HCC1395BL_WES.mutect2.chr21.f1r2.tar.gz" \
  --output /tmp/nomics-m2-f1r2-priors/golden.tar.gz
cmp /tmp/nomics-m2-f1r2-priors/port.tar.gz /tmp/nomics-m2-f1r2-priors/golden.tar.gz
```

The integration test repeats chr21 archive/prior comparisons without Java,
explicitly skipping with a message if the golden directory is absent. Three
unit tests cover signed BQ, median/sample/depth membership, ALT ties, unpaired
flags, table depth, insertion side and integer indel thresholds. Scope is the
inherited default two-BAM, one-normal/one-tumour traversal, BQ >20, median MQ
>=50, depth cap200, and the shared tar writer's 100-ASCII-byte member-name
limit (after URL encoding). Other collector/caller modes are not exposed or
claimed. No GPU or speedup claim. This completes a collector boundary, not the
final Mutect2/VCF command. A mid-session environment restart removed temporary
logs; the full required checks were rerun for the evidence below.

Collector gate: **7 pass, 0 fail, 0 skipped**. Every member-name set is exact;
every reference histogram and ALT table byte is exact. On each contig the
singleton histogram compares **384 unique data columns × 200 ordered depth
rows = 76,800 exact numeric-text values**, plus the depth column and every
non-column byte. No internal-gate fallback was necessary.

| Contig | ref_histogram exact bytes | alt_table exact bytes | Collector wall s |
|---|---:|---:|---:|
| chr21 | 43,305 | 134,856 | 11 |
| chr22 | 47,323 | 556,295 | 48 |
| chr19 | 52,512 | 1,129,143 | 280 |
| chr20 | 49,384 | 522,844 | 34 |
| chr13 | 46,519 | 262,134 | 19 |
| chr16 | 50,996 | 946,564 | 60 |
| chr12 | 52,794 | 754,747 | 73 |

These are script wall times including comparison, under variable shared-Mac
memory pressure, not benchmarks. The initial overlapping test/regions runs
were stopped by exact verified process IDs to release memory; only the
collector gate remained active. The final test/regions checks were restarted
serially. No other worktree's process was stopped.

The chr21 prior CLI check passes on the final code, including complete
nomics-produced tar.gz bytes. Both contain the 15,362-byte
`./HCC1395_HCC1395_WES.orientation_priors`, SHA256
`67087bb44890f3f59c6deeb4ecb4b5035e7a588a60ad8a136db9026bb7924e38`.

Stage-one CLI regression: **7 pass, 0 fail, 0 skipped**, regions and stats
by bytes and all ordered filter counts/totals exact. Wall seconds in contig
order chr21/22/19/20/13/16/12: **7 / 29 / 98 / 27 / 16 / 54 / 125**.
Final release build passes (25.72 s), as does the required Clippy command
(32.18 s). The sole build/lint warning is the existing dependency HC
Atomic::fetch_update deprecation, not a warning in the modified packages.
`git diff --check`, rustfmt checks on the new module/test and shell syntax
checks pass. The forbidden-transcendental scan on the added implementation
and integration test has zero hits. An additional javap audit of the unmodified
jar's TableWriter.close confirms header emission even for an empty ALT table.

Required release tests: **PASS**, exit 0; **138 passed, 0 failed, 1 existing
ignored**, across 15 test/doc-test summaries. No new ignores. The new chr21
collector/prior integration test executes and passes (8.78 s); all seven
caller goldens pass (685.66 s), three-contig assembly passes (152.69 s), and
five-contig likelihood/read-realignment passes (348.70 s). Optional internal
oracle dumps were not supplied; no new internal-dump verification is claimed.

The lead's machine-wide-lock directive arrived while this final `cargo test`
command was already running. It was allowed to finish as explicitly directed;
no second heavy job was started. The CLI gate runs above also predate that
directive. All subsequent heavy commands use `lockf -k /tmp/nomics-heavy.lock`;
future reproductions of the test/gate commands above must use that prefix too.

Peak memory, `/usr/bin/time -l`, with one traversal at a time under the global
lock: chr19 is the largest slice by processed read count (20,756,074), chr12
is the longest chromosome among the seven. Both repeated complete archives
also compare byte for byte with the earlier gated outputs.

| Slice | Maximum resident bytes | RSS GiB | Peak memory footprint bytes | Wall s |
|---|---:|---:|---:|---:|
| chr19 | 10,734,518,272 | 9.9973 | 12,953,429,520 | 80.12 |
| chr12 | 9,942,745,088 | 9.2599 | 10,168,901,168 | 67.07 |

The first locked timing attempt completed chr19's traversal but `/usr/bin/time`
failed at `sysctl kern.clockrate: Operation not permitted` inside the sandbox,
without reporting RSS. The successful measurements above reran with sandbox
escalation approved, still holding the same machine-wide lock. No memory
optimization or source change was made. Actual command: `lockf -k
/tmp/nomics-heavy.lock bash /tmp/nomics-m2-f1r2-memory.sh
"$PWD/target/release/nomics"`; the helper loops over chr19/chr12 serially, using
`/usr/bin/time -l` on the exact collector CLI from the gate script, writing
`/tmp/nomics-m2-f1r2-memory/CONTIG.{f1r2.tar.gz,time.txt}` and comparing each
archive with `/tmp/nomics-m2-f1r2-gate/CONTIG.f1r2.tar.gz`. Its exit status is 0.
The lock was released when the helper finished. No heavy process from this
worker remains running. No required gate failed; there are no residual
collector payload or downstream-prior differences under the two proposed order
normalizations. No push, pull request or merge was performed.
