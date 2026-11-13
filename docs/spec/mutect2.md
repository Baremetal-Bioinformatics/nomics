# Mutect2 tumour-normal specification (GATK 4.6.2.0)

Target: the two-scatter SEQC2 HCC1395 / HCC1395BL exome run of sarek 3.10.0,
`tools: mutect2`. This is a specification for a fresh implementation, not a claim
that nomics implements somatic calling. Filtering and auxiliary tools are in
[mutect2-filtering.md](mutect2-filtering.md).

## 0. Sources, evidence and boundaries

**[source]** means code inspected at the following pinned revisions. GATK paths
below are relative to `src/main/java/org/broadinstitute/hellbender/`.
A citation naming a class abbreviates its package only after that package is
introduced. Method references are deliberate: they identify the complete rule,
including branches which have no example in the published files.

* [GATK 4.6.2.0](https://github.com/broadinstitute/gatk/tree/4.6.2.0), commit
  `76edc75c26504da94bbaee66584e107e76ee15de`.
* [htsjdk 4.2.0](https://github.com/samtools/htsjdk/tree/4.2.0), commit
  `4cc010022ac038fb30f26e6f9717fabff3e808c1`; paths relative to `src/main/java/htsjdk/`.
* [Picard 3.4.0](https://github.com/broadinstitute/picard/tree/3.4.0), commit
  `6c3f23bc2e0d229d75e9f9b04200396bcd067526`; paths relative to `src/main/java/picard/`.
* [Commons Math MATH_3_5](https://github.com/apache/commons-math/tree/MATH_3_5),
  commit `b3c5dae8f253fcb4484e5cd3cc5662587803efc2`, fetched separately under
  `/tmp/mutect2-spec-math3` for the numerical checks in the filtering spec.
* GATK `build.gradle`, `htsjdkVersion`, `picardVersion`, `gklVersion`, and the
  `commons-math3` dependency pin: 4.2.0, 3.4.0, 0.8.11, and strictly 3.5.
  Source trees were fetched under `/tmp/mutect2-spec-{gatk,htsjdk,picard}`, outside
  this repository. No upstream implementation is vendored. **[source]**

**[golden]** means inspected locally under
`~/omics/runs/somatic/S1-r1/results/variant_calling/mutect2/`.
Let `X = HCC1395_WES_vs_HCC1395BL_WES`; pair files are in `X/`.
Normal and tumour pileups are in the respective `HCC1395BL_WES/` and
`HCC1395_WES/` directories. BAM sample names are
`HCC1395_HCC1395BL_WES` and `HCC1395_HCC1395_WES` (normal first in both VCFs).
`params.yaml` and `samplesheet.csv` at the run root retain the run inputs.
The BED is `S07604624_Padded_Agilent_SureSelectXT_allexons_V6_UTR.bed`,
not the Agilent v7 BED used by the germline HC spec. **[golden: params.yaml]**

The raw and filtered VCFs each contain **6,351** records; **2,024** filtered
records are PASS. Raw FILTER and QUAL are `.` throughout. ALT counts are
5,867 biallelic, 343 two-alt, 117 three-alt, 22 four-alt and 2 five-alt records.
The merged raw header has 3,412 lines. **[golden: X.mutect2{,.filtered}.vcf.gz]**
The final section of `docs/REFERENCE-VARIABILITY.md` reports the independent
r1/r2 comparison: records and auxiliary table bytes repeat, tar contents repeat
but artifact-prior archive bytes do not. That comparison was not rerun here.

The prescribed chain is Mutect2 per BED group → Picard MergeVcfs plus
MergeMutectStats; LearnReadOrientationModel consumes every group's f1r2 tar;
GetPileupSummaries per sample/group → GatherPileupSummaries;
CalculateContamination with matched normal; FilterMutectCalls on the merged VCF.
The raw golden's retained Mutect2 command names the group
`chr17_38877931-38878296`; it is not evidence that only that group was run.
**[golden: raw ##GATKCommandLine; task's recorded commands; REFERENCE-VARIABILITY]**

### 0.1 Recorded command contract

Placeholders below stand for the staged basenames and exact scatter group BED;
all omitted arguments are defaults in §1. These reproduce the supplied golden
command contract, not newly executed commands. **[task record; golden caller and
filter command headers]**

```text
gatk Mutect2 --input NORMAL.recal.bam --input TUMOR.recal.bam --output X.mutect2.<group>.vcf.gz --reference Homo_sapiens_assembly38.fasta --panel-of-normals 1000g_pon.hg38.vcf.gz --germline-resource af-only-gnomad.hg38.vcf.gz --intervals <group>.bed --tmp-dir . --f1r2-tar-gz X.mutect2.<group>.f1r2.tar.gz --normal-sample HCC1395_HCC1395BL_WES
gatk MergeVcfs --INPUT <each group VCF> --OUTPUT X.mutect2.vcf.gz --SEQUENCE_DICTIONARY Homo_sapiens_assembly38.dict
gatk MergeMutectStats --stats <each group stats> --output X.mutect2.vcf.gz.stats
gatk LearnReadOrientationModel --input <each f1r2 tar> --output X.mutect2.artifactprior.tar.gz
gatk GetPileupSummaries --input S.recal.bam --variant af-only-gnomad.hg38.vcf.gz --reference Homo_sapiens_assembly38.fasta --intervals <group>.bed --output S.mutect2.<group>.pileups.table
gatk GatherPileupSummaries --I <each sample group table> --O S.mutect2.pileups.table --sequence-dictionary Homo_sapiens_assembly38.dict
gatk CalculateContamination --input TUMOR.pileups.table --matched-normal NORMAL.pileups.table --output X.mutect2.contamination.table -tumor-segmentation X.mutect2.segmentation.table
gatk FilterMutectCalls --variant X.mutect2.vcf.gz --output X.mutect2.filtered.vcf.gz --reference Homo_sapiens_assembly38.fasta --orientation-bias-artifact-priors X.mutect2.artifactprior.tar.gz --tumor-segmentation X.mutect2.segmentation.table --contamination-table X.mutect2.contamination.table
```

## 1. Effective arguments

All following values are in the raw golden's Mutect2 command header and agree
with the argument classes. The header prints `--af-of-alleles-not-in-resource -1.0`;
the effective value is **1e-6**, selected because a normal was supplied.
`M2ArgumentCollection.getDefaultAlleleFrequency()` uses 5e-8 without a normal.
Internal emission/activity thresholds are converted from log10 to natural logs.
**[source: tools/walkers/mutect/M2ArgumentCollection.java, getters] [golden]**

| Area | Effective settings | Source at pinned tag |
|---|---|---|
| Calling | initial tumour LOD 2.0; emit LOD 3.0; normal LOD 2.2; minimum AF 0.0 | `M2ArgumentCollection.java L53–64, L263–282, L333, L367`; threshold getters |
| Resource shortcuts | genotype PoN sites false; genotype germline sites false; germline fraction 1.0; maximum population AF 0.01 | `M2ArgumentCollection.java L132–148, L287–348` |
| PCR/activity | SNV quality 40; indel quality 40; substitution correction 5; ignore ITR artifacts false | `M2ArgumentCollection.java L132–148, L287–348` |
| Modes | mitochondria false; flow NONE; reference confidence NONE; independent mates false; no forced alleles; no Permutect output | `M2ArgumentCollection.java L95–115, L186–245, L367–374`; inherited `AssemblyBasedCallerArgumentCollection.java` |
| Downsampling | max reads/start 50; stride 1; max suspicious reads/start 0 (disabled); callable depth 10 | `M2ArgumentCollection.java L314–324`; `engine/spark/AssemblyRegionArgumentCollection.defaultMaxReadsPerAlignmentStart` |
| Regions | min 50; max 300; padding 100; activity threshold 0.002; propagation distance 50; force active false | `AssemblyRegionArgumentCollection.java L70–158`, default getters |
| Trimming | SNP padding 20; indel padding 75; STR padding 75; legacy trimming false; legacy max extension 25 | `AssemblyRegionArgumentCollection`; `AssemblyRegionTrimmer.trim` |
| Assembly | k 10,25; increase k on failure/cycles allowed; nonunique reference kmers disallowed; pruning samples 1; adaptive pruning ON | `haplotypecaller/MutectReadThreadingAssemblerArgumentCollection.makeReadThreadingAssembler` |
| Pruning | initial error 0.001; extension ln threshold 2.302585092994046; seed ln threshold 9.210340371976184; max unpruned variants 100 | `ReadThreadingAssemblerArgumentCollection.java L133–154` |
| Graph | min dangling length 4; recover all dangling false; recovery enabled; max haplotypes 128; linked de Bruijn false; legacy cycle detection false; matching bases -1 | `MutectReadThreadingAssemblerArgumentCollection.makeReadThreadingAssembler` |
| Correction | error-correct reads false; error-correction LOD -Infinity; assembly minimum BQ 10 | `ReadThreadingAssemblerArgumentCollection`; `AssemblyBasedCallerArgumentCollection.createReadThreadingAssembler` |
| Events | maximum MNP distance **1**, rather than HC's 0 | `M2ArgumentCollection.getDefaultMaxMnpDistance` |
| Likelihoods | PairHMM; FASTEST_AVAILABLE; native threads 4; native double precision false; GCP 10; global mismapping Q45 | `LikelihoodEngineArgumentCollection.java L61–109`; native arguments |
| Likelihood qualities | PCR indel CONSERVATIVE; BQ threshold 18; cap BQ to MAPQ; symmetric normalization true; dynamic disqualification false; expected mismatch rate 0.02 | `LikelihoodEngineArgumentCollection.java L49–113` |
| Clipping | use soft-clipped bases when fragment permits; override fragment check false; soft-clip low-quality ends false | `AssemblyBasedCallerArgumentCollection.java L112–205` |
| SW | FASTEST_AVAILABLE; dangling 25,-50,-110,-6; hap→ref 200,-150,-260,-11; read→hap 10,-15,-30,-5 | `AssemblyBasedCallerArgumentCollection` SW parameter getters |
| Pileup/flow extras | pileup detection false; PDHMM false; flow homopolymer collapse 0; allele prefilter false | inherited caller/likelihood argument collections |
| f1r2 | min BQ 20; median MQ 50; depth cap 200 | `readorientation/CollectF1R2CountsArgumentCollection.java L13–23` |
| Intervals/output | UNION; merge ALL; padding 0; validation SILENT; index true; md5 false; command line true; sites-only false | `IntervalArgumentCollection`; `GATKTool`; golden command |

`min-pruning 2` is printed but the Mutect assembler passes prune factor **0**
when adaptive pruning is enabled. `disable-artificial-haplotype-recovery false`
is passed directly to `setArtificialHaplotypeRecoveryMode`, despite the negative
option name. Preserve the executable code rather than interpreting the name.
**[source: MutectReadThreadingAssemblerArgumentCollection.makeReadThreadingAssembler]**

## 2. Traversal, filters and read transformations

Mutect2 is an `AssemblyRegionWalker`. BED decoding, dictionary sorting,
merging overlapping/abutting intervals, padded queries, activity smoothing,
read assignment and readiness of regions use the shared implementation described
in `haplotypecaller.md` §§1–2. `AssemblyRegionWalker.makeReadShards` makes a shard
per contig; each shard queries its union of intervals padded by 100, but activity
is evaluated at the unpadded traversal intervals. The chr21 implementation and
end-to-end query checks are recorded in §11. **[source: engine/AssemblyRegionWalker.java, makeReadShards,
traverse, processReadShard; engine/AssemblyRegionIterator.java]**

The stream is pre-transformer → read filters → post-transformer → downsampler.
Mutect2 does **not** override the pre-filter transformer: it inherits identity;
HC installs `IUPACReadTransformer`. Post-filter Mutect2 installs palindrome
artifact clipping by default. **[source: Mutect2.java, makePostReadFilterTransformer;
engine/GATKTool.java, makePreReadFilterTransformer;
engine/MultiIntervalLocalReadShard.java, iterator]**

With two BAM inputs, `engine/ReadsPathDataSource.prepareIteratorsForTraversal`
queries each reader and wraps them in htsjdk MergingSamRecordIterator using
SamFileHeaderMerger. Coordinate comparator ties use SAMRecordCoordinateComparator's
full comparison (strand/name/flags/MQ/mate coordinates/TLEN after position),
not blanket “normal before tumour”. RG/PG collisions are translated by the
header merger. Confirm exact tie streams on the two original recal BAMs before
using reservoir results as a caller oracle. **[source: ReadsPathDataSource.
prepareIteratorsForTraversal/createHeaderMerger; htsjdk samtools/MergingSamRecordIterator.next,
ComparableSamRecordIterator.compareTo; SAMRecordCoordinateComparator.compare]**

### 2.1 Filters, in execution order

All must pass; short circuit at the first failure. **[source:
Mutect2Engine.makeStandardMutect2ReadFilters; engine/filters/ReadFilterLibrary.java;
MappingQualityReadFilter.test; ReadLengthReadFilter.test; WellformedReadFilter.test]**

1. MAPQ ≥20.
2. MAPQ available (not 255).
3. MAPQ nonzero (redundant with 1 at the default).
4. Mapped.
5. Not secondary; supplementary reads remain eligible.
6. Not duplicate.
7. Passes vendor quality check.
8. Non-chimeric original alignment: **only if both OA and the original mate-contig
   XM tag exist**, require their contigs equal. This is not the current mate-contig
   test used by HC after assembly (`ReadFilterLibrary.NonChimericOriginalAlignmentReadFilter.test`).
9. Nonzero reference-consuming alignment length.
10. Read length 30 through `Integer.MAX_VALUE`, inclusive.
11. Good CIGAR: htsjdk validity, no adjacent indel operators, no terminal deletion
    after ignoring clipping (`utils/read/CigarUtils.isGood`).
12. Well formed: alignment coordinates/header agreement, read group, matching
    SEQ/QUAL and CIGAR lengths, stored bases, no `N` CIGAR operator.

Compared with HC's initial filter list, add 3, 8 and 10. Do not add a
not-supplementary or proper-pair filter. **[source: Mutect2Engine.makeStandardMutect2ReadFilters;
haplotypecaller/HaplotypeCallerEngine.makeStandardHCReadFilters]**

### 2.2 Palindrome artifact clipping

`transformers/PalindromeArtifactClipReadTransformer.apply` requires proper
pairing and a computable adaptor boundary. With positive TLEN it examines the
first CIGAR element, otherwise the last; only `S` or `I` qualifies. Let its length
be a and n=min(a+5, read length). For positive TLEN compare reference
[boundary−n,boundary−1] with the reverse complemented first n read bases;
otherwise compare [boundary+1,boundary+n] with the reverse complemented last n.
Out-of-contig reference, or reference crossing the indicated alignment boundary,
returns the read unchanged. Match fraction **≥0.9** hard-clips the candidate a
bases; negative-TLEN clipping passes stop=`read.length` to `ClippingOp`.
This runs before activity and f1r2 collection. **[source: that method;
Mutect2Engine.MIN_PALINDROME_SIZE=5]**

### 2.3 MutectDownsampler

`utils/downsampling/MutectDownsampler` pools starts from the first start in a
stride up to but excluding first+stride; a contig change also closes the pool.
Default stride 1 means one start, **across both samples**, not 50 per sample.
Up to 50 pending reads: keep all in arrival order. Above 50: discard every read
with MAPQ ≤50, then reservoir-sample the remaining reads to 50; do not refill
from low-MAPQ reads if fewer than 50 remain. Finalized lists are stably sorted
by assigned start. **[source: submit, handlePositionalChange,
finalizePendingReads, consumeFinalizedItems]**

The suspicious-read switch 0 maps to `Integer.MAX_VALUE`; otherwise MAPQ ≤50
increments a counter and reaching the limit rejects the entire stride.
Reservoir replacement consumes `Utils.getRandomGenerator()` (shared seed
47382911), preserving reservoir-slot order among equal starts. Each scattered
process starts its own RNG; changing scatter boundaries can change retained
reads. HC's `PositionalDownsampler` does not do this preliminary MAPQ selection.
**[source: MutectDownsampler constructor/submit; ReservoirDownsampler.submit;
utils/Utils.java, randomGenerator; HC spec §1.5]**

## 3. Activity and callable sites

`Mutect2Engine.isActive` applies the following rules **in order**. Outputs are
binary activity probabilities 0 or 1, not HC's diploid posterior. **[source]**

1. A passing forced allele gives activity 1 before all other work (no forced
   alleles in this chain).
2. Null/empty pileup gives 0.
3. If combined pileup size ≥10, increment callable-sites **before** testing
   tumour evidence, normal evidence or PoN. Size is the pileup's size, not
   usable-base depth and not tumour-only depth.
4. Extract tumour pileup; feed it to F1R2CountsCollector before activity rejection.
5. Accumulate six quality buffers: A,C,G,T, other substitution, and all indels.
   Choose the largest quality sum, first index on ties (`>` update, starting A).
6. Compute tumour ln likelihood ratio with nRef=`tumourPileup.size−nAlt`.
   If `< ln(10)*2`, activity 0; equality survives.
7. With a normal and germline genotyping off, accumulate normal buffers. Only
   when its winning category equals the tumour winner, reject if normal alt
   count **>0.3×normal pileup size AND** normal quality sum **>100**. This is a
   count/quality shortcut, not a normal LOD calculation.
8. Without a normal and with germline genotyping off, use resource AF ≥0.01
   shortcuts. The substitution branch compares the tumour substitution base
   with **germline REF[0]**, despite the comment saying dominant alt. Any indel
   category matches any resource indel. This branch is not used in this pair.
9. If genotype-PoN-sites is false and any PoN record overlaps the single locus,
   return 0; allele identity is irrelevant.
10. A surviving site with depth <10 also increments callable-sites; return 1.

### 3.1 Activity qualities and approximate likelihood

`Mutect2Engine.PileupQualBuffer.accumulateQuals` checks current/following indel
first: use Q=min(30+10×(length−1),127). Next to a useful soft clip use indel Q30:
current and immediately adjacent clipped base must both have BQ **>6**.
Otherwise a nonreference substitution with BQ **>6** goes to its base buffer.
If proper-paired with mapped mate and
`mateStart ≤ position < mateStart+read.length`, cap its BQ at PCR-SNV-Q/2=20.
Then add substitution correction 5, capped at `QualityUtils.MAX_QUAL` (254).
An unrecognized substitution goes to the other buffer **without** the +5.
Ref and ignored bases still contribute to nRef above. **[source: buffer methods;
isNextToUsefulSoftClip; getCurrentOrFollowingIndelLength; indelQual]**

The ln LLR is a variational approximation. Put A=nAlt, R=nRef, N=A+R,
`t=FastMath.exp(MathUtils.digamma(R+1)−MathUtils.digamma(A+1))`.
For each alt Q, e=10^(−Q/10), z=(1−e)/(1−e+e×t), add
`z×(ln(1−e)−ln(e)) + Hfast(z)`.
`Hfast(z)=p×(11+33p)/(2+20p)`, p=z(1−z), **not exact Bernoulli entropy**.
Add `−ln(N+1)−binomialCoefficientLog(N,A)` for the flat beta prior.
The repeated-alt overload multiplies the read sum by repeatFactor and sets
A=repeatFactor×list length; adaptive pruning uses this overload after
`QualityUtils.errorProbToQual` quantizes its error estimate.
**[source: Mutect2Engine.logLikelihoodRatio overloads;
utils/MathUtils.fastBernoulliEntropy, digamma; utils/NaturalLogUtils quality methods]**

### 3.2 Smoothing and cuts

The shared `utils/activityprofile/BandPassActivityProfile` uses maximum filter
size 50, sigma 17.0, normalized width-101 Gaussian kernel and adaptive truncation
at 1e-5. `ActivityProfile` compares smoothed probability **>0.002**. An active
run reaching 300 is cut at the lowest qualifying local minimum searched from
index 299 down to 49; output end is index+1. Gaps force conversion and discard
propagated tail outside the observed span. Padding is 100 clipped to the contig.
There are no HIGH_QUALITY_SOFT_CLIPS activity states from Mutect2, so HC's
extra soft-clip propagation is not taken. **[source: engine/AssemblyRegionIterator;
BandPassActivityProfile constructor/processState; ActivityProfile.findEndOfRegion,
findBestCutSite; Mutect2Engine.isActive; detailed shared algorithm in HC §2.4]**

## 4. Region preparation and assembly

`Mutect2Engine.callRegion` first adjusts overlapping mates on the original
region, **even before the inactive/empty-region return**. Agreeing bases get
BQ=min(BQ,20); conflicting bases keep their qualities (`setConflictingToZero=false`).
Insertion/deletion qualities in overlaps are each capped at 20, materializing
BI/BD arrays when absent. The index walk is base-for-base across overlaps and
does not align indels between mates; preserve that limitation.
HC instead zeros conflicting BQs and supplies no indel-quality cap.
**[source: callRegion; utils/fragments/FragmentUtils.adjustQualsOfOverlappingPairedFragments;
AssemblyBasedCallerUtils.cleanOverlappingReadPairs]**

After the early return, remove unmarked duplicates among paired reads with a
mapped mate on another contig or |TLEN|>1,000,000. Group by sample and signed
unclipped start (first-of-pair positive, otherwise negative), then mate contig.
Keep only the first read of a mate-contig subgroup if its size is **>half** the
whole group; discard all others, including every subgroup without a strict
majority. Do not apply HC's later blanket mate-on-other-contig exclusion.
**[source: Mutect2Engine.removeUnmarkedDuplicates, callRegion]**

Shared finalization hard-clips or reverts soft clips according to fragment
geometry, trims tails at Q≤9, adaptor-clips, clips to padded region and stably
sorts with ReadCoordinateComparator. `assembleReads(...,false,...)` disables its
own overlapping-base correction because Mutect2 did it already. It fetches the
same extra 500 reference bases and reference haplotype as HC.
**[source: haplotypecaller/AssemblyBasedCallerUtils.assembleReads/finalizeRegion;
Mutect2Engine.callRegion; shared clipping details HC §3]**

ReadThreadingAssembler otherwise follows the shared k=10,25 graph construction,
cycle retries, graph simplification, dangling-end recovery, k-best extraction,
haplotype/reference SW, event construction and trimming. Graph input sequences
are grouped **by sample**, with ref first and sample groups in insertion order;
MultiSampleEdge flushes per sample. Do not collapse both samples into an
unlabelled HC read list. Bases below assembly BQ10 or N split a read's usable
runs. Candidate kmers retry at +10 up to six attempts if initial sizes fail.
**[source: haplotypecaller/readthreading/ReadThreadingAssembler.runLocalAssembly,
assemble, createGraph; AbstractReadThreadingGraph.addRead/buildGraphIfNecessary;
HC §4 for shared graph rules]**

### 4.1 Adaptive pruning is the principal graph difference

`haplotypecaller/graphs/AdaptiveChainPruner.chainsToRemove` runs
`likelyErrorChains` at error 0.001; estimates error rate as the sum of last-edge
multiplicities of those error chains divided by total multiplicities over every
edge of every chain; then reruns at that estimate. Final removal never removes
a chain containing a reference edge. HC's default is LowWeightChainPruner(2).
**[source: those methods; MutectReadThreadingAssemblerArgumentCollection]**

The same pruneFactor **0** also reaches `recoverDanglingTails` and
`recoverDanglingHeads`; HC supplies 2 there. Adapting only chain removal while
leaving the dangling-path cutoff at 2 is therefore insufficient. Pruning runs
before cycle counting with legacy cycle detection disabled. **[source:
readthreading/ReadThreadingAssembler.createGraph, calls at L654–689]**

For each chain end, compare its terminal edge count against total multiplicity
of the incident outgoing/incoming edges with Mutect2's approximate LLR; a
source/sink side has LLR 0. Build good incoming/outgoing adjacency when that
side LLR ≥2.302585092994046 or the first edge is reference. Seed-worthy chains
must pass **both** 9.210340371976184 side thresholds. Seed vertices require
more than two incident seedable chains; also seed the max-weight chain.
**[source: AdaptiveChainPruner.chainLogOdds, likelyErrorChains]**

Max-weight ties compare largest edge multiplicity, then chain length, then
first-vertex bases. The growth priority queue sorts descending LLR, then
first-vertex bases, then complete chain bases. Track visited vertices and good
chains in insertion-ordered sets. Adding another outgoing chain at a previously
used source counts as a variant. The loop permits entry while count **≤100**;
the chain that increments beyond 100 has already been added to the good set
before expansion is skipped. Preserve that boundary behavior.
**[source: AdaptiveChainPruner.getMaxWeightChain, likelyErrorChains]**

### 4.2 Trimming and event representation

`getVariationEvents(maxMnpDistance=1)` coalesces adjacent substitution events
according to EventMap; this affects MNP alleles, region trimming and filters.
Use the shared nonlegacy `AssemblyRegionTrimmer.trim`: variant span intersected
with original active span, SNP/indel/STR expansion, then limit by original padded
span. `AssemblyResultSet.trimTo` trims haplotypes and reads; identical trimmed
haplotypes collapse with reference winning, then sort by length and base order.
Mutect2 removes stubs of **read length <10** before PairHMM.
**[source: Mutect2Engine.callRegion/removeReadStubs; utils/haplotype/EventMap;
haplotypecaller/AssemblyRegionTrimmer.trim; AssemblyResultSet.trimTo]**

## 5. PairHMM, realignment and fragments

`utils/pairhmm/PairHMM.Implementation.FASTEST_AVAILABLE` tries native OpenMP AVX,
then single-thread AVX, then Java LoglessPairHMM, catching HardwareFeatureException
at each native attempt. This is capability/load dependent, not an unconditional
x86/arm switch. Linux arm64 without Intel GKL takes Java double-precision
LoglessPairHMM; supported x86 with loaded GKL normally takes OpenMP native
(default 4 native threads, float mode requested). The local final VCF proves
FASTEST_AVAILABLE was requested, not which library loaded; somatic task logs
are needed to confirm the fallback for this run. SmithWaterman similarly
selects the fastest loadable implementation. **[source: PairHMM enum;
utils/smithwaterman/SmithWatermanAligner.getAligner; golden command]**

The same `PairHMMLikelihoodCalculationEngine` is used as HC: original read bases,
cloned BQs and BI/BD, default indel Q45 **except arrays capped at overlaps**;
CONSERVATIVE repeat-dependent PCR indel cap; base Q capped to MAPQ then Q<18→6;
indel Q<6→6; GCP10; tristate mismatch division by 3. Java recurrence and initial
2^1020 scaling are those in HC §6. No DRAGstr input is supplied.
**[source: PairHMMLikelihoodCalculationEngine.modifyReadQualities,
applyPCRErrorModel, capMinimumReadQualities; LoglessPairHMM]**

The engine caps each likelihood at best−4.5 (symmetric). Poorly modeled reads
are discarded with the shared static mismatch rule, using HMM-quality length,
expected rate .02 and cap of two errors: reject when best log10 likelihood
< `−4×min(2,ceil(HMM-quality-array length×.02))`. Then `Mutect2Engine.callRegion` converts
the entire matrix **from log10 to ln before read realignment**. Shared SW
realignment uses the read/haplotype scores above. Best-haplotype priority is
reference bonus plus 1−CIGAR-element-count; close-tie window is ln(10)×0.2
because `getInformativeThreshold()` checks matrix scale.
**[source: computeReadLikelihoods; ReadLikelihoodCalculationEngine.filterPoorlyModeledEvidence;
Mutect2Engine.callRegion; AssemblyBasedCallerUtils.realignReadsToTheirBestHaplotype;
utils/genotyper/AlleleLikelihoods.switchToNaturalLog/searchBestAllele]**

**Separate threshold trap:** `AlleleLikelihoods.BestAllele.isInformative()` tests
`confidence > LOG_10_INFORMATIVE_THRESHOLD` (literal **0.2**) even on a natural-log
matrix. It does not use getInformativeThreshold. Thus annotation informative
counts must use **ln difference >0.2**, not >0.4605170185988092. This is source
behavior, not yet independently isolated with a BAM-level golden.
**[source: AlleleLikelihoods.java L38–50, BestAllele constructor/isInformative]**

`SomaticGenotypingEngine.callMutations` caps the ln matrix again with
`NaturalLogUtils.qualToLogErrorProb(45)` before grouping. `groupEvidence` groups
**per sample by read name** through Java HashMap (`Collectors.groupingBy`),
then sums all group members' log likelihoods for each haplotype, in group-member
order. It is not max-over-mates. `Fragment.createAndAvoidFailure` uses one/two input reads unchanged when the
input group has at most two reads (including any supplementary flags). For a
larger group, it filters duplicate/secondary/supplementary reads, takes the
**first** two survivors when more than two remain (despite the warning saying
“randomly”), uses the one/two survivors otherwise, or falls back to the first
original read when none survive. Group likelihoods still include **all** group
members, including those omitted from the representative interval. Fragment overlap uses the whole spanning interval,
including the gap between mates. **[source: SomaticGenotypingEngine.callMutations;
AlleleLikelihoods.groupEvidence; utils/read/Fragment constructors/createAndAvoidFailure]**

## 6. Somatic genotyping, in source order

`tools/walkers/mutect/SomaticGenotypingEngine.callMutations` rebuilds EventMaps,
uses sorted event starts in the genotyping active span, and for each start:

1. Collect haplotype events with `includeSpanningEvents=false` (HC uses true).
   Merge through `AssemblyBasedCallerUtils.makeMergedVariantContext` /
   `GATKVariantContextUtils.simpleMerge`: longest REF, right-padding shorter
   alleles, reference then first-seen distinct ALT order. Do not sort ALT strings.
2. Make allele→haplotype mapper; marginalize **fragment** likelihoods by maximum
   over supporting haplotypes. The LinkedHashMap mapper can append a `*` matrix
   allele for a preceding overlapping event even though the merged VC excluded
   spanning events. This additional allele participates in variational evidence
   and absence LODs, but is not one of the merged VC's emission candidates.
   Retain fragments overlapping VC±2 bases.
3. Concatenate tumour matrices in sample order; compute per-alt somatic ln LOD.
   Compute normal diploid absence ln LOD and normal somatic artifact ln LOD.
4. Keep tumour alts whose LOD **>3×ln(10)** (strict). Make a further list whose
   normal absence LOD **>2.2×ln(10)** (strict), unless forced/mode overrides.
   If this further list is empty, skip the record. If any survives, emit **all
   tumour-threshold alts**, including the germline-like alts at that site.
5. Independently accumulate potential somatic Event objects in a HashSet using
   both strict thresholds; this later determines ECNT/ECNTH.
6. Obtain population frequencies for emitted alleles, set POPAF/TLOD/NLOD/NALOD;
   annotate PON if any PoN entry starts at that position (no allele matching).
7. Build genotypes, trim common suffix and prefix (`trimAlleles(true,true)`),
   remap fragment and read allele likelihoods to trimmed alleles; annotate.
8. After all loci, phase calls; then compute event-count annotations.

**[source: callMutations; AssemblyBasedCallerUtils.getVariantsFromActiveHaplotypes,
createAlleleMapper; utils/variant/GATKVariantContextUtils.simpleMerge/trimAlleles]**

### 6.1 SomaticLikelihoodsEngine numerics

`somaticLogOdds` computes variational logEvidence(all alleles) minus
logEvidence(excluding this alt), not a diploid genotype comparison. Empty
matrices yield zero evidence. Priors are one pseudocount per allele at minAF=0;
nonzero minAF would set alt prior `1−ln(2)/ln(minAF)` while ref stays 1.
**[source: SomaticGenotypingEngine constructor/makePriorPseudocounts/somaticLogOdds]**

Posterior iteration starts with ones. For each fragment, compute normalized
responsibilities from ln likelihoods plus `digamma(alpha_i)−digamma(sum alpha)`;
sum responsibilities and add prior. Stop when L1(old,new)/sum(new) **<0.001**;
there is no maximum-iteration parameter in this loop. Posterior means are
alpha/sum(alpha), not responsibilities normalized after a single iteration.
**[source: SomaticLikelihoodsEngine.alleleFractionsPosterior/getEffectiveCounts;
utils/Dirichlet.effectiveLogMultinomialWeights]**

Let C(alpha)=logGamma(sum alpha)−sum logGamma(alpha_i). Evidence is
C(prior)−C(posterior) plus, fragment by fragment,
`sum_i responsibility_i×lnLikelihood_i − sum_i xLogx(responsibility_i)`.
The likelihood term is zeroed for responsibility **<1e-10**; xLogx is zero for
x **<1e-8**. Preserve both cutoffs and the different accumulation operations:
IndexRange.sum, MathUtils.sumArrayFunction, and Java double-stream sum are not
interchangeable reductions. **[source: SomaticLikelihoodsEngine.logEvidence,
likelihoodsContribution/xLogx/logDirichletNormalization]**

`NaturalLogUtils.logSumExp` finds the first maximum, starts sum at 1, skips that
index and −Infinity terms, adds exp(term−max) in index order, and returns `max+(sum==1 ? 0 : ln(sum))`. The addition of +0 is retained,
including its effect on a negative-zero maximum. This is not HC's Jacobian log10 sum.
Commons Math Gamma/digamma and FastMath remain numeric dependencies to match.
**[source: utils/NaturalLogUtils.logSumExp/normalizeLog; utils/Dirichlet]**

Normal absence LOD is `sum lnL(ref) − sum(logSumExp(lnL(ref),lnL(alt))+ln(.5))`.
There is **no hom-alt likelihood term**, despite the NLOD header mentioning
“diploid het or hom alt”. NALOD uses the normal's somatic evidence difference
with the same sign as TLOD; filtering negates it when modeling absence of a
normal artifact. Stored TLOD/NLOD/NALOD are converted back to log10.
**[source: SomaticGenotypingEngine.diploidAltLogOdds/callMutations]**

### 6.2 Genotypes and population resource

Initial genotype AD is rounded soft responsibility count from a flat per-read
normalization of the **fragment** allele matrix; AF is the Dirichlet posterior
mean excluding reference. Annotation subsequently **overwrites AD with hard
informative read counts**. AF remains the fragment posterior, so AF≠ADalt/DP
is expected. A no-evidence genotype gets uniform AF prior and zero counts.
Normal alleles are always two copies of REF. Tumour GT is the entire allele
list: e.g. `0/1/2`; it is not a ploidy estimate. GQ/PL are never assigned here.
**[source: SomaticGenotypingEngine.addGenotypes/getEffectiveCounts;
annotator/DepthPerAlleleBySample.annotate]**

Population resource lookup uses alleleIndices to match equivalent right-padded
representations at the current locus. Start with 1e-6 for every emitted allele;
only accept a resource AF list of length Nalleles−1. A null or `.` element
becomes the default via Mutect2Engine's fixed double-list reader; an absent
attribute yields no usable values. Malformed numeric text raises an error. Later resource entries
can overwrite an earlier match. POPAF=−log10(AF). No germline posterior is
written by the caller; GERMQ belongs to filtering.
**[source: SomaticGenotypingEngine.getGermlineAltAlleleFrequencies,
getNegativeLogPopulationAFAnnotation; Mutect2Engine.getAttributeAsDoubleList]**

## 7. Annotation contract

The standard group is `annotator/StandardMutectAnnotation`. All annotations see
the trimmed, realigned **original-quality** reads, not the quality-modified
PairHMM copies. Jumbo annotations additionally receive fragment-allele and
fragment-haplotype likelihoods. Count informativeness uses §5's literal .2.
**[source: Mutect2Engine.getStandardMutect2AnnotationGroups;
SomaticGenotypingEngine.callMutations; annotator/VariantAnnotatorEngine.annotateContext]**

| Field | Computation and ordering | Source (annotator package unless stated) |
|---|---|---|
| INFO DP | total evidenceCount in read allele likelihoods, including uninformative reads | `Coverage.annotate` |
| FORMAT AD | informative best-allele read counts per sample, ref then ALT; overwrites initial soft fragment AD | `DepthPerAlleleBySample.annotateWithLikelihoods` |
| FORMAT DP | informative read count in emitted-allele subset; equals sum AD | `DepthPerSampleHC.annotate` |
| FAD | same hard-count routine as AD but using fragments | `FragmentDepthPerAlleleBySample.annotate` |
| AF | fragment posterior means, ref omitted; no hard-count recomputation | `SomaticGenotypingEngine.addGenotypes` |
| F1R2/F2R1 | informative best fragment, inspect **first read**: usable MQ≠0/255 and base quality at VC ≥20; bin by pair orientation; ref plus each alt | `OrientationBiasReadCounts.annotate` |
| SB | per-sample informative read counts `[refF,refR,allAltF,allAltR]`; minimum sample count >0 | `StrandBiasBySample.annotate`; `StrandBiasTest.getContingencyTable` |
| AS_SB_TABLE | per-allele read strand counts across samples; include a sample only with total informative count **>2**; encode refF,refR then `|altF,altR` per alt | `AS_StrandBiasMutectAnnotation.annotate`; `allelespecific/StrandBiasUtils` |
| MBQ | median read base quality at VC by allele; usable MQ and informative; empty=0 | `BaseQuality`; `BaseQualityRankSumTest.getReadBaseQuality` |
| MMQ | median read MAPQ by allele, including REF; empty=60 | `MappingQuality` |
| MFRL | median absolute TLEN by allele, including REF; empty=0 | `FragmentLength` |
| MPOS | median min distance to read end, ALT only; includes terminal hard clips as original bases; empty=50 | `ReadPosition`; `ReadPosRankSumTest.getReadPosition` |
| ECNT | number of distinct threshold-passing somatic Event objects in region, not simply output record count | `SomaticGenotypingEngine.callMutations` |
| ECNTH | per ALT: choose supporting haplotype with most best-read assignments, stable first on ties; count its events present in potential-somatic set | `SomaticGenotypingEngine.callMutations` |
| TLOD,NLOD,NALOD | §6 models, per ALT in emitted order, log10 | `SomaticGenotypingEngine.callMutations` |
| PON | flag for any PoN entry starting at the call position, not an exact REF/ALT match | same method |
| POPAF | §6.2, per ALT | `getNegativeLogPopulationAFAnnotation` |
| RPA,RU,STR | indels only, tandem-repeat routine succeeds: repeat counts ref then ALT, repeat-unit string, flag | `TandemRepeat.annotate/getNumTandemRepeatUnits`; `GATKVariantContextUtils.getNumTandemRepeatUnits` |
| PGT,PID,PS | shared physical phasing of called haplotypes, described below | `AssemblyBasedCallerUtils.phaseCalls` |

`PerAlleleAnnotation.annotate` pools samples for the medians. ReadUtils.isF1R2 is `reverseStrand != firstOfPair`; F2R1 uses equality,
without an additional proper-pair requirement. **[source: utils/read/ReadUtils.isF1R2/isF2R1]**
Integer median is
`FastMath.round(Commons Math Median.evaluate(...))`, including rounding the mean
of the two middle values. Base-quality lookup uses aligned read index and the
shared rank-sum helper; MPOS uses leading-H+offset and length−1−offset+trailing-H,
then takes the smaller. It is not the distance to the indel's far endpoint.
**[source: PerAlleleAnnotation; utils/MathUtils.median(int[]);
ReadPosRankSumTest.getReadPosition; BaseQualityRankSumTest.getReadBaseQuality]**

Tandem-repeat annotation excludes the common anchor base, examines reference
context to the right, chooses the repeated unit from the longer non-anchor
allele and counts repetitions in alleles plus reference continuation. Failure
for any alt suppresses the whole annotation. Use the full shared routine for
mixed multiallelic indels rather than independently normalizing each ALT.
**[source: TandemRepeat; utils/variant/GATKVariantContextUtils.getNumTandemRepeatUnits overloads]**

Physical phasing is enabled here, unlike plain-VCF HC. Shared phasing maps
variants to haplotypes, connects compatible same/opposite haplotype sets, and
assigns PGT `0|1` or `1|0`, PID from the first call's position/REF/ALT, PS from
its position; GT is marked phased too. Preserve eligibility and conflicts in
`AssemblyBasedCallerUtils.constructHaplotypeMapping`, `constructPhaseSetMapping`,
`constructPhaseGroups` and `phaseCalls`. **[source: these methods]**
Only calls with exactly one site-specific alternate allele enter the haplotype
map (ignore NON_REF and `*`). Compare calls in output order. Equal haplotype sets,
or the appropriate all-alt-haplotypes containment, give same phase; disjoint sets
whose sizes sum to the total alt-haplotype union give opposite phase. An all-alt
call narrows its available haplotypes to its first partner. If an unassigned
call tries to create a group with an already-assigned partner, clear **all**
phase mappings and return unphased. Existing partner mappings are otherwise
left as assigned. A group must contain at least two calls. PID has exact form
`<position>_<REF>_<firstALT>`; all samples receive PGT/PID/PS, including hom-ref
normal. Heterozygous GT alleles are reversed if necessary to put ALT in the
phase group's indicated slot; hom-ref GT still becomes phased. **[source:
constructHaplotypeMapping/constructPhaseSetMapping/constructPhaseGroups,
createUniqueID/phaseVC, in AssemblyBasedCallerUtils]**

There are 1,010 phased golden records; normal GT `0|0`, tumour 977 `0|1` and 33
`1|0`. Unphased FORMAT is always `GT:AD:AF:DP:F1R2:F2R1:FAD:SB` (5,341 records);
phased adds `PGT:PID:PS` before SB. **[golden]**

### 7.1 Header declaration cardinality and type

The following is the complete caller INFO/FORMAT declaration inventory in
serialized order. GQ/PL and filtering quality declarations remain even where
no caller record emits them. Descriptions come from standard htsjdk lines or
GATKVCFHeaderLines; retain their literal golden wording, including NLOD and
AF wording discussed above. **[golden: raw VCF header; source:
Mutect2Engine.writeHeader; utils/variant/GATKVCFHeaderLines.java static definitions;
htsjdk variant/vcf/VCFStandardHeaderLines.java]**

| Kind | ID | Number | Type |
|---|---|---|---|
| FORMAT | AD | R | Integer |
| FORMAT | AF | A | Float |
| FORMAT | DP | 1 | Integer |
| FORMAT | F1R2 | R | Integer |
| FORMAT | F2R1 | R | Integer |
| FORMAT | FAD | R | Integer |
| FORMAT | GQ | 1 | Integer |
| FORMAT | GT | 1 | String |
| FORMAT | PGT | 1 | String |
| FORMAT | PID | 1 | String |
| FORMAT | PL | G | Integer |
| FORMAT | PS | 1 | Integer |
| FORMAT | SB | 4 | Integer |
| INFO | AS_SB_TABLE | 1 | String |
| INFO | AS_UNIQ_ALT_READ_COUNT | A | Integer |
| INFO | CONTQ | 1 | Float |
| INFO | DP | 1 | Integer |
| INFO | ECNT | 1 | Integer |
| INFO | ECNTH | A | Integer |
| INFO | GERMQ | 1 | Integer |
| INFO | MBQ | R | Integer |
| INFO | MFRL | R | Integer |
| INFO | MMQ | R | Integer |
| INFO | MPOS | A | Integer |
| INFO | NALOD | A | Float |
| INFO | NCount | 1 | Integer |
| INFO | NLOD | A | Float |
| INFO | OCM | 1 | Integer |
| INFO | PON | 0 | Flag |
| INFO | POPAF | A | Float |
| INFO | ROQ | 1 | Float |
| INFO | RPA | R | Integer |
| INFO | RU | 1 | String |
| INFO | SEQQ | 1 | Integer |
| INFO | STR | 0 | Flag |
| INFO | STRANDQ | 1 | Integer |
| INFO | STRQ | 1 | Integer |
| INFO | TLOD | A | Float |

Raw INFO fields present on every record: AS_SB_TABLE, DP, ECNT, ECNTH, MBQ,
MFRL, MMQ, MPOS, NALOD, NLOD, POPAF, TLOD. RPA/RU/STR occur together on 1,602;
PON occurs on 637. GERMQ, ROQ and STRQ are **not** raw annotations. Other declared
but absent fields include CONTQ, SEQQ, STRANDQ, NCount, OCM,
AS_UNIQ_ALT_READ_COUNT, GQ and PL. Declaration is not emission. The filtered
records add AS_FilterStatus/GERMQ/ROQ everywhere and STRQ on those 1,602 STRs.
**[golden: both VCFs; source: Mutect2Engine.writeHeader;
filtering/Mutect2FilteringEngine.applyFiltersAndAccumulateOutputStats]**

## 8. Extra outputs and serialization

### 8.1 Stats

`Mutect2.onTraversalSuccess` writes `<output>.stats` even when no variants were
emitted. `MutectStats` uses a TSV header `statistic\tvalue`; callable is the
counter in §3 converted to double, rendered with Java Double.toString, newline
terminated. `MergeMutectStats.doWork` sums callable across its insertion-ordered
set of input files (double stream sum); unknown statistic names fail. No union
or overlap correction is applied to callable bases across scatter groups.
**[source: these methods; utils/tsv/DataLine.set(double)]**
The published merged table is exactly `statistic\tvalue\ncallable\t1.20404485E8\n`
(38 bytes; SHA256 `934c6dd019d4b5e6b2186290195544cd3b5b2561491e88e88fa6b6f586bf4a46`).
It is a filtering input controlling learned somatic priors. **[golden; filtering spec §4]**

### 8.2 F1R2CountsCollector

Collection occurs during activity evaluation, from tumour reads after filters,
ITR clipping and downsampling but **before** region mate-quality adjustment.
Require a complete 3-base reference context without N. Per sample, retain
pileup elements with BQ **>20** (different from the annotation's ≥20).
BaseCounts counts A/C/G/T. Reject zero depth, or indel evidence count
(deletion/after insertion/before deletion start) **> integer depth/100**, or
median MQ<50. Use all eligible pileup mapping qualities in that median.
**[source: readorientation/F1R2CountsCollector.process/isPileupGood]**

Replace ref count by −1 in a temporary A,C,G,T array and choose first max to
select the dominant ALT. Zero alt count increments the reference-context
histogram at min(depth,200). Exactly one ALT read increments the context,
alt-base, F1R2/F2R1 histogram at capped depth. Larger ALT count writes an
AltSiteRecord with context, ref/alt counts and ref/alt F1R2 counts; other ALT
bases do not enter those two counts. Early `return` branches in process leave
the entire method, not just the current sample; relevant if adding multi-tumour
support later. **[source: F1R2CountsCollector.process; DepthOneHistograms.increment]**

Each tumour contributes URL-encoded sample basename plus `.ref_histogram`,
`.alt_histogram`, `.alt_table`. Ref histograms cover 64 contexts, depth bins
1..200; alt histograms cover 6 substitutions/orientations per context.
Picard MetricsFile emits sample StringHeader and histograms; AltSiteRecord
uses TSV metadata SAMPLE and its explicit column enum. Temporary-directory
`listFiles()` is archived by IOUtils.writeTarGz. File/member order and archive
metadata need a content-based comparison proposal; do not assume tar bytes.
**[source: F1R2CountsCollector constructor/writeHistograms/closeAndArchiveFiles;
F1R2FilterConstants; F1R2FilterUtils; AltSiteRecord; utils/io/IOUtils.writeTarGz]**
The alt-table column order is `context,ref_count,alt_count,ref_f1r2,alt_f1r2,depth,alt`
(tab separated); depth=ref+dominant-alt, excluding other ALT bases. Histogram
labels are context or `<context>_<ALT>_<F1R2-or-F2R1>`, with key label `depth`.
All depth bins are prefilled, including zeros. MetricsFile prints
`## HISTOGRAM\tjava.lang.Integer`, column order from supplied histogram lists,
keys in numeric order. FormatUtil uses at most 6 fractional digits,
HALF_DOWN, grouping disabled; whole-valued counts have no decimal suffix.
Reference columns follow context HashMap order; alt columns additionally inherit
enum/pair-key HashMap order. **[source: AltSiteRecord column enum/constructor;
F1R2FilterUtils; DepthOneHistograms.getHistograms; htsjdk samtools/metrics/MetricsFile.printHistogram;
samtools/util/FormatUtil constructor/format]**

Per-group f1r2 archives are absent locally. The final learned prior tar contains
one member `./HCC1395_HCC1395_WES.orientation_priors`; its 15,642-byte payload
contains 64 context rows, with 12 state probabilities and example counts.
**[golden; learning details in filtering spec §1]**

### 8.3 VCF writing and merge

`Mutect2Engine.writeHeader` adds MutectVersion=2.2, unfiltered warning,
annotation declarations and reserved filtering fields, standard GT/AD/GQ/DP/PL,
AF and phasing definitions, tumour_sample and normal_sample, and the BAM
sequence dictionary. VCFHeader/VCFWriter sort metadata; FORMAT declarations
precede GATKCommandLine then INFO, with contigs in dictionary order. Samples
come from ReadUtils.getSamplesFromHeader, sorted by sample name. Do not derive
normal/tumour order from input option order. **[source: writeHeader;
utils/read/ReadUtils.getSamplesFromHeader; htsjdk variant/vcf/VCFHeader;
variant/variantcontext/writer/VCFWriter.writeHeader]**

The htsjdk encoder sorts INFO keys with TreeMap. Genotype keys are unioned and
sorted, GT first. Numeric integer arrays are comma-joined; true flag has no
`=value`. Missing QUAL is `.`; raw unfiltered status is `.`, filters-applied
empty set becomes PASS. Filter names are sorted at site serialization.
**[source: htsjdk variant/vcf/VCFEncoder.encode/formatVCFField/getFilterString;
variant/variantcontext/VariantContext.calcVCFGenotypeKeys]**

For a Double d, htsjdk `VCFEncoder.formatVCFDouble` uses Locale.US:
if d≥1, `%.2f`; else if d≥.01, `%.3f`; else if |d|≥1e-20, `%.3e`;
otherwise `0.00`. **Every negative number** takes the scientific branch,
including NALOD=−1.301e+00. Already-decoded VCF strings need not be reformatted:
FilterMutectCalls preserves input attributes and genotypes it does not replace.
`String.format` rounding must follow Java formatting, not Rust's default float
Display or a blanket precision. **[source: VCFEncoder.formatVCFField/formatVCFDouble]**

Per-group traversal is sequential: contigs, assembly regions, sorted event
starts; Mutect2.apply writes callRegion's returned order. Picard
`vcf/MergeVcfs.doWork` uses smartMergeHeaders and a MergingIterator comparing
dictionary contig then start; it **does not deduplicate** variants. Equal-position
records need the actual source iterator/heap tie behavior, not an added allele
sort. SEQUENCE_DICTIONARY supplies the writer dictionary; an input's nonempty
contig list remains the comparator basis and must be compatible with others.
**[source: Mutect2.apply; SomaticGenotypingEngine.callMutations;
Picard MergeVcfs.doWork; htsjdk samtools/util/MergingIterator]**

smartMergeHeaders identifies compound lines by key+ID, otherwise by key; keeps
the first conflicting ordinary line. Hence different scatter Mutect2 command
lines collapse to the first input header's line, consistent with the one
retained golden command. MergingIterator inserts input iterators into a Java
PriorityQueue in input-list order, polls the minimum, consumes one item and
reinserts the iterator. Its comparator adds no secondary key for equal-position
VCFs: preserve heap ties if such records occur. **[source: htsjdk
variant/vcf/VCFUtils.smartMergeHeaders; samtools/util/MergingIterator constructor/next]**

Compression/index writing shares htsjdk's BGZF and Tabix layout specified in
`haplotypecaller.md` §9 and `docs/GATES.md`: compressed offsets depend on header
bytes, block cuts and deflate implementation. Compare normalized decompressed
VCF first; a same-input rewrite is needed to gate compressed bytes and .tbi.
This document does not authorize new normalization. **[source: htsjdk
VCFWriter/TabixIndexCreator; existing gates]**

## 9. Port plan

This is a proposed implementation/gating order, not a new accepted gate.
Somatic rows must be added to `docs/GATES.md` with justification before any
normalization is used. Today its HC command-line exception is not automatically
an exception for arbitrary somatic fields or tar contents. **[project rule]**

| Boundary | File/oracle and proposed comparison | Current golden |
|---|---|---|
| Pileup production | exact gathered TSV bytes, including sample metadata and row order | both sample `.pileups.table`, 246,147 normal / 244,297 tumour rows |
| Contamination | exact TSV bytes; preserve segmentation's Java map contig order | `.contamination.table`, `.segmentation.table` (42 rows) |
| Caller | all raw VCF records and header, proposing only GATKCommandLine provenance exception | merged `.vcf.gz`; per-group files must be retrieved |
| Callable bases | exact TSV bytes, not rounded callable total | merged `.vcf.gz.stats` |
| F1R2 collection | member names and member payload bytes; archive wrapper normalization requires GATES entry | per-group tar unavailable here |
| Orientation learning | prior member bytes in original row order; no probability tolerance; wrapper normalization requires GATES entry | `.artifactprior.tar.gz` payload |
| Filtering | exact records, header with justified command-line exception; exact stats TSV bytes | `.filtered.vcf.gz`, `.filteringStats.tsv` |
| Compression/index | two-sided writer test: rewrite golden with original header; index our own stream equivalently | raw/filtered `.vcf.gz.tbi` |
| Internal caller stages | instrument pinned GATK: filtered read stream, activities/regions, graph/pruning, haplotypes, likelihoods, fragment groups, LODs | needs recal BAMs and group BEDs from Spark |

### 9.1 Reuse assessment from this checkout

This assessment is based on inspected Rust files, not just package names.

* Reuse arithmetic kernels in `crates/nomics-hc/src/pairhmm.rs`,
  `pairhmm_trie.rs`, `assembly/sw.rs`, `assembly/cigar.rs`; Java-compatible
  `jmath.rs`, `javafmt.rs`, `jrandom.rs`, `tables.rs` provide existing numerical
  machinery. Reuse BGZF/VCF numeric formatting and index primitives in nomics-io
  and `vcf.rs`, with a generalized record schema. **[local source]**
* Reuse smoothing/cutting from `activity.rs`, region/read clipping operations in
  `finalize.rs`, haplotype trimming in `trim.rs`, haplotype alignment and event
  construction in `assembly/mod.rs`, and SW composition in `realign.rs`.
  Each needs an explicit Mutect mode at the orchestration boundary: different
  activity evidence; no IUPAC pre-transform; overlap policy and indel arrays;
  max-MNP-distance 1; ln-aware realignment tie threshold. **[local source]**
* `assembly/graph.rs` contains fixed prune-at-2 logic (`prune`, dangling-path
  checks) and reference/single-sample flush state; add adaptive pruning and correct multi-sample
  threading. Do not merely change a global cutoff. `readprep.rs.hmm_qualities`
  currently synthesizes Q45 indel arrays: accept preexisting overlap-capped
  arrays. `likelihood.rs` already implements mismatch rejection and cap−4.5;
  retain its CPU kernel but add scale/sample orchestration. **[local source]**
* `call.rs` has shared event merging/marginalization ideas but assumes HC
  spanning-event behavior and diploid genotype reduction. `genotype.rs` is HC's
  germline model and is not the somatic evidence engine. `annotate.rs` uses
  single-sample/diploid annotation structures and log10 thresholds; reuse its
  coordinate helpers, not its whole pipeline. `vcf.rs.Record` has one genotype,
  mandatory germline fields and one sample in Header: extend for multiple
  samples, AF/FAD/orientation fields, missing QUAL and filter state.
  `caller.rs` owns single-sample traversal/downsampling and must gain sample
  routing and scattered job state. **[local source: these types/functions]**
* New components: Mutect activity/LLR, palindrome clipper, MutectDownsampler,
  fragment likelihood container with Java name-map ordering, variational somatic
  genotyper, population/PoN resource queries, ECNT/ECNTH, f1r2 collection,
  orientation EM, pileup/contamination tools, clustering and every filtering pass.
  These are proposed fresh modules, not existing crates. **[port recommendation]**

Port in four steps. First build TSV/VCF readers and auxiliary pileup,
contamination and orientation-learning tools against fixed published inputs;
retrieve f1r2 tar inputs before claiming learning parity. Second implement
filtering on the golden raw VCF, gating its four passes and final stats, so
caller differences cannot hide behind filtering. Third port read/region/assembly
mode differences and gate internal states on each exact group. Fourth add
somatic likelihoods, annotations, phasing and serialization, then gate the full
scatter/merge/filter chain. **[recommendation based on the dependency boundaries above]**

### 9.2 Likely costs and variation

Floating-point reduction order is part of the target: fragment HashMap values,
allele maps, per-sample concatenation, stream compensated sums versus explicit
loops, Gamma/digamma/FastMath, polynomial entropy and Java formatting. Do not
replace approximate LLR or variational evidence with a mathematically cleaner
integral. `SomaticLikelihoodsEngine` and filtering learn from those differences.
**[source: §§3,5,6; SomaticLikelihoodsEngine; filtering spec]**

A Java-compatible map traversal matters even if the set is unchanged: name
HashMap ordering defines fragment columns and hence iterative sums. Reference
and haplotype order controls equal-likelihood choices. The adaptive-pruning
queue explicitly supplies base-sequence tie breaks; do not erase them.
Somatic genotype RNG is unseeded JDKRandomGenerator, but its nextDouble call
short-circuits when genotypeGermlineSites=false, and Permutect is disabled.
The default downsampler uses a seeded shared RNG. **[source:
AlleleLikelihoods.groupEvidence; AdaptiveChainPruner; SomaticGenotypingEngine.callMutations]**

Do not assume region preparation is independent merely because existing HC
preparation is parallel. `AssemblyRegionIterator.next` loads the next ready
region before handing the current region to apply; `fillNextAssemblyRegionWithReads`
and `AssemblyRegion.add` share read object references. Mutect2 overlap correction
mutates those reads before finalizeRegion makes defensive copies. Preserve
traversal/call timing until a BAM-level oracle settles effects across overlapping
regions and future pileups. **[source: these methods; Mutect2Engine.callRegion;
AssemblyBasedCallerUtils.HardClipAndPossiblyAddToCollection]**

Native x86 PairHMM can change low bits versus Java arm64; this target must be
confirmed with the somatic task logs before optimizing for native/GPU parity.
A matching Java path is the first target. Changing scatter membership changes
RNG streams, query/padding context, assembly boundaries and callable totals.
Archive metadata and command-line dates/paths vary independently of scientific
content. No performance claim or GPU shipping decision follows from this spec.
**[source: PairHMM.Implementation; AssemblyRegionWalker; MutectDownsampler path
in §2; IOUtils.writeTarGz; project rules]**

### 9.3 Checks performed for this specification

Read-only checks on all records confirmed raw/filtered coordinate, REF/ALT and
genotype-column identity; 6,351 total and 2,024 PASS; DP=sum(AD) for every sample
genotype; 64 prior contexts with probability sums within 1e-12 for the inspection
check (this tolerance is **not a port gate**). The tables, headers, field presence,
phasing and hashes recorded above were inspected directly. No reference tool
was rerun and no executable nomics code was changed. **[golden inspection]**

## 10. Not confirmed

These are source-derived rules or proposed gates, not independently validated
internal calculations. Published output observations above are confirmed.

* Filtered read order/counts, ITR clips, overlap edits (including BI/BD), duplicate
  removal, downsampling RNG draws, exact active states/regions, adaptive-pruning
  choices, haplotypes and likelihood matrices: settle with **both recal BAMs,
  each scatter BED and an instrumented 4.6.2.0 caller**. The chr21 recal BAMs now
  provide exact filtered-count, activity, region,
  read-assignment and pre-assembly quality gates (§11). Assembly internals
  remain unconfirmed.
* Exact somatic evidence/AF iteration trajectories, fragment map order, .2-ln
  informative-boundary effects, ECNTH ties and phasing groups: same inputs plus
  stage dumps. Final VCF values alone cannot prove their independent derivation.
* Per-group VCF/header merge provenance, equal-coordinate merge ties and
  callable partials: retrieve **per-group VCFs and .stats, MergeVcfs command/log**.
  Merged values and final ordering are locally observable.
* Collector histogram column order, AltSiteRecord bytes and f1r2 tar member
  names/order: retrieve **per-group f1r2 tars**. Final priors do not establish the
  individual input counts or histogram layouts.
* Actual PairHMM/SW/deflater library selected for this somatic run and precise
  JDK build: retrieve **Mutect2 .command.log**. Header arguments alone prove
  requested settings only. Intel GKL native arithmetic was not audited here;
  no x86/native parity claim is made.
* Full Commons Math 3.5 and JDK math/formatting implementation parity, Java
  HashMap tree-bin behavior, and extreme/unsupported-mode branches were not
  executed. Use small JVM oracles plus BAM-level boundary dumps; alternate
  modes (tumour-only, mito, flow, forced alleles, gVCF) need separate fixtures.
* Compression/index byte parity and any tar/header normalization: requires
  same-header rewrite tests plus an explicit `docs/GATES.md` addition. This
  documentation change does not silently create those exceptions.

## 11. Implemented first-stage boundary (2026-11-12)

`nomics somatic mutect2-regions` implements the default sarek tumour-normal
traversal described in §§1–3, plus callable-stat serialization from §8.1.
`--input` is repeated twice, with `--normal-sample`, `--reference`, `--intervals`
(BED), `--panel-of-normals`, `--germline-resource`, `--assembly-region-out`
and `--stats`. Germline resource existence is checked, but matched-normal
activity does not query it (source's `hasNormal` branch). No VCF or f1r2 is
produced. Proposed gates, awaiting lead acceptance, are in `docs/GATES.md`.

The library's `run_with_sink` visits every region in order with its assigned
reads, after the original-region overlapping-mate adjustment, before assembly
or unmarked-duplicate removal. Reads carry sample identity, alignment, bases,
BQ, BI and BD. The arena remains shared across regions. The callback consumes
references immediately; retaining/replaying independent region copies would
lose the source's mutation timing. One ready-region lookahead is retained:
region R is called only after R+1 has been filled. Inactive regions are cleaned
as well. Optional `--dump-prefix` emits all raw activity states, ordered read
assignments, and binary read payloads after adjustment.

### 11.1 Source confirmations and corrections

Pinned GATK tag 4.6.2.0 was fetched into `/tmp/mutect2-spec-gatk` and checked at
`76edc75c26504da94bbaee66584e107e76ee15de`. The already fetched htsjdk 4.2.0
source was inspected alongside it. The source rules in §§1–3 were confirmed
on the chr21 recalibrated BAMs, including cross-sample MQ selection before
reservoir draws, input-head merging using the full coordinate comparator
(without re-sorting each BAM's equal-start reads), adaptor exclusion from
pileups, indel/soft-clip evidence, and PoN overlap independent of alleles.

No effective default-argument discrepancy was found in §§1–3. Earlier claims
that the BAMs/internal states were unavailable are superseded by the checks
below. Additional source details needed beyond the prose spec (the source wins):

* Wellformed's header-agreement check limits **alignment start**, not end, to
  contig length (`ReadUtils.alignmentAgreesWithHeader`). ValidAlignmentEnd
  checks `end-start+1 >= 0`. Do not copy HC's stricter local end-bound check
  into this caller. GoodCigar also includes htsjdk padding validity and its
  repeated-indel checks, beyond adjacent indel rejection.
* FragmentUtils orders paired reads by **soft start**, reversing on equal
  soft starts. Its overlap offset treats S as reference-consuming and D as
  returning the following read offset; the subsequent base-array loop ignores
  internal indels. This existing source behavior is preserved, not repaired.
* Activity quality sums add **signed Java bytes**, including corrected
  substitutions capped at 254; quality lookup tables interpret the same bytes
  unsigned. Default chr21 reads do not exercise the high-byte edge, but the
  independent numerical fixture includes Q254.
* The Java dump must use `convertToSAMRecord(...).getFlags()` for original
  flags. `ReadUtils.getSAMFlagsForRead` clears mate-strand flags for unmapped
  mates: it changed two dump rows at regions 6206/6207. This was a harness
  serialization mistake, not a traversal difference; no port normalization
  was added.

Reuse: HC's `Intervals`, `HcRead`, `JavaRandom`, hard clipping and the exact
existing band-pass profile/cuts. `BinaryActivityProfile` is a new public wrapper
around the unchanged private HC implementation, never invoked by HC;
`hard_clip_bases` only changes visibility. No HC numerical/filter behavior is
changed. Large reference-byte buffers use `nomics_mem::HostBuf`; only contigs in the
traversal are loaded. The somatic filter helpers are local, derived from the same small HC
helpers but corrected to the upstream rules above. Traversal orchestration,
activity evidence and Mutect downsampling are separate somatic code because
sample handling and shared-read edits differ from HC. All transcendental math
uses existing Java-exact ports, including Commons Math digamma/FastMath.exp,
Commons binomialCoefficientLog and the NaturalLogUtils quality-cache formula.

### 11.2 Evidence and remaining limits

A local Homebrew OpenJDK **17.0.21.0.3+0** first ran the unmodified sarek GATK
jar with the exact chr21 scientific arguments. It reproduced golden regions,
stats, all filter counts, and VCF bytes after removing the command-line
provenance line. This is not the container's 17.0.17-internal build. Its f1r2
reference histogram and alt table match, but alt-histogram column order differs
(the enum/pair-key hash iteration discussed in §8.2); no f1r2 parity is claimed.
Only after the requested first-stage boundaries reproduced was the same JVM
used for `reference/oracle/mutect2/MutectRegions.java`. That walker uses the
unmodified engine/filters/downsampler/iterator and source overlap cleaner,
skipping assembly. Its regions and stats also match the unmodified full run.
`MutectNumerics.java` independently supplies 42 exact LLR bit fixtures, a
143-read cross-sample downsampling fixture (63 retained, including the
low-MQ no-refill case), and a 1,000-locus two-impulse smoothing/cutting fixture. No oracle is linked into or
called by nomics.

The chr21 results supersede the earlier §2 and §10 statements that recal BAM
traversal/activities could not be checked locally: 7,167 regions (699 active),
callable 1,230,643, all twelve filter counts and totals are confirmed. Internal
activity and read-assignment dumps are also compared exactly. Detailed commands,
payload results, timings and required Rust checks are in the verification
record. No new accepted normalization follows from these observations.

**Not supported / not confirmed:** only the fixed default two-BAM, one-normal,
one-tumour, coordinate-sorted indexed-BAM/BED contract is exposed. Different BAM
dictionaries and sample cardinalities fail explicitly. Tumour-only, multiple
tumours/normals, CRAM/SAM, forced alleles, mitochondrial/flow/gVCF modes, custom
filters, interval/downsampling/region/quality thresholds, disabled ITR clipping
and genotype-PoN/germline switches are not exposed. Non-default inputs need
independent fixtures before parity is claimed; this is one chr21 scatter, not
full-exome or arbitrary-input validation. No assembly, likelihoods, genotyping,
annotations, VCF, f1r2, or post-assembly preparation is implemented here. The
future assembly stage may consume the region callback without re-traversing
reads; it must preserve callback order. No speedup or GPU claim is made.

## 12. Stage-two implementation and source confirmations (2026-11-12)

`nomics somatic mutect2-haplotypes` consumes the ordered stage-one callbacks.
It emits an internal assembly dump (`--output`) and the stage-determined
artificial-BAM projection (`--artificial-records`). These are proposed gates,
awaiting the lead's acceptance, not a somatic variant caller. No PairHMM,
genotyping, annotations, VCF or f1r2 work is invoked by the Rust command.

The source remains tag 4.6.2.0 at
`76edc75c26504da94bbaee66584e107e76ee15de`, fetched outside the repository.
The following details refine §§4 and 9.2; no threshold was fitted to a golden:

* **Sample insertion order is per k-mer attempt and per usable run.**
  `AbstractReadThreadingGraph.addRead` adds a sample to its pending LinkedHashMap
  only when a BQ10/non-N run has length **at least k**. An earlier finalized
  short read does not establish sample order for that attempt. Each attempt
  groups runs in that order, reference first, and flushes MultiSampleEdge once
  per sample. Merely grouping finalized reads before assembly is incorrect.
  This explained the first development difference, chr21 region 2477:
  a six-base tumour read precedes a normal 41-base read, so normal threads first
  at k25. The first differing quantity had been path score/order, not reads or
  bases; fixing the pending insertion rule fixed all three differing regions.
* **Adaptive error rates are quantized.**
  `Mutect2Engine.logLikelihoodRatio(refCount, altCount, errorProbability)` calls
  `QualityUtils.errorProbToQual` (`Math.round(-10*Math.log10(error))`, bounded
  to [1,93]), then evaluates **one** quality with repeatFactor=altCount. It
  multiplies that one read contribution at the end, rather than repeatedly
  accumulating it. The digamma arguments still use the full counts. Zero
  estimated error becomes Q93. This is more specific than §4.1's continuous
  error-rate prose. All transcendental operations use existing Java-exact ports.
* The shared finalizer has `correctOverlappingBaseQualities=false` here;
  stage one has already edited BQ/BI/BD on shared reads, including inactive
  regions. Copies are finalized after unmarked-duplicate removal. Hard clipping
  retains the corresponding contiguous slices of BQ/BI/BD. BI/BD do not enter
  the assembly graph, but the internal dump includes them for the next worker.
* Trimming uses EventMap distance 1 **within each match CIGAR element**.
  Trimmed non-reference haplotypes retain their original score and k-mer size.
  Reconstructing the reference in `AssemblyResultSet.trimDownHaplotypes` retains
  the score but omits k-mer size, leaving zero; this source behavior is retained.
* Artificial BAM records have `RG:Z:ArtificialHaplotypeRG`, names
  `HC_<global counter starting at 1>_<region index starting at 0>`, and flag 16.
  `HC:i` is Java `Arrays.hashCode(bases)`; `CR:Z` is the trimmed callable span.
  These three tags are all the tags observed on the artificial records and all
  are determined by this stage. No later-stage tag is excluded. The default
  writer does not enable its optional special/collapsed tags. The gate compares
  contig, 1-based position, CIGAR, bases and all three tags **in BAM file order**.
  The coordinate writer orders equal-position artificial records by read name;
  the port retains the same global names for this ordering (names themselves
  are not a projected field).
* `Mutect2Engine.callRegion` actually calls `writeBamOutput` **after**
  `genotypingEngine.callMutations`, using the uncollapsed likelihood allele
  list. In this fixed non-flow mode, likelihood calculation/realignment does
  not change that list. It is the trimmed list, written whenever the region
  passes both no-variation exits, regardless of whether genotyping emits a
  variant. Realigned biological reads remain outside this task.

Reuse from HC: default finalization gained an explicit overlap-correction
parameter with its existing entry point still passing true. Assembly gained
an optional chain-removal callback and sample IDs; its existing entry point
passes neither, retaining fixed pruning/dangling cutoff 2 and one read sample.
The callback receives graph-owned maximal chains in traversal order, so the
adaptive LLR implementation stays in nomics-somatic without a dependency cycle.
The new mode passes dangling cutoff 0. The SW, CIGAR, threading, simplification,
k-best code and k-mer retry schedule are shared. The existing trimming entry
point still uses MNP distance 0; a new entry point accepts distance 1. Mutect2
uses only its haplotypes, avoiding HC's later genotyping read filters.

Ordering audit: `BaseGraph.getSources` returns LinkedHashSet; chain discovery
uses FIFO starts and outgoing-edge order. Graph vertices/edges use jgrapht's
insertion-preserving storage. MultiDeBruijnVertex equality is object identity;
its cached hash is System.identityHashCode. Seedable multimap key iteration is
therefore **not assumed insertion order**: all seed queues are populated before
extraction, whose comparator is descending LLR, then ascending first-vertex and
complete-chain bases. The Rust implementation preserves that comparator and
chain discovery order; graph integer IDs are never used as a queue tie-break.
Full assembly dumps, including path scores and ordering, are the parity check.
See the verification record for the tested scope and exact commands.

### 12.1 Not supported / not confirmed

Only the existing sarek two-indexed-BAM/BED contract is exposed. Custom k-mer,
pruning/quality/region thresholds, error correction, linked de Bruijn graphs,
legacy cycle/trimming modes, forced/pileup-injected alleles, flow/homopolymer
collapse, custom SW scores and alternate ploidy/sample/mode switches are not
exposed. No optional user argument silently selects one of those modes.

Object-identity hash iteration and **fully comparator-equal distinct chains**
are not proved for arbitrary graphs/JVMs. The small oracle graphs cover equal
weights, reference protection, zero estimated error, and the 100-variant limit;
the contig gates cover the actual path order on these inputs. They are not a
proof of JVM identity hashes on arbitrary data or every HC golden. No new
normalization or tolerance is proposed to conceal this limitation. Extended
HC fixture tests remain ignored unless the independent fixtures are supplied.

## 13. Stage-three likelihood boundary (2026-11-12)

`nomics somatic mutect2-likelihoods` runs the existing ordered traversal and
assembly followed by Java-mode PairHMM, likelihood capping/disqualification,
natural-log conversion and best-haplotype realignment. `--output` is an ordered
internal diagnostic dump; `--realigned-records` is the biological BAM projection,
not a BAM. It stops before `SomaticGenotypingEngine.callMutations`: fragment
grouping (§5), the genotyper's second likelihood cap, annotations, VCF and f1r2
are outside this command. Proposed gates await lead acceptance in `docs/GATES.md`.

### 13.1 Source audit and ordering

Source: GATK 4.6.2.0, commit `76edc75c26504da94bbaee66584e107e76ee15de`,
and htsjdk 4.2.0, inspected outside the repository. The actual golden
`oracle/gatk.log` explicitly reports fallback to **LOGLESS_CACHING**, superseding
§5/§10's uncertainty about the implementation selected on this arm64 run.

The rules in §5 were confirmed, including CONSERVATIVE PCR caps on the existing
BI/BD arrays, symmetric best-minus-4.5 capping, static poorly-modelled-read
filtering, multiplication by `Math.log(10)`, the ln-aware tie window and the
literal-0.2 informative test. No default numerical constant disagreed with §5.
The following additional details matter at this boundary:

* `AssemblyRegion.trim` hard-clips finalized reads to the trimmed **padded**
  span and sorts with `ReadCoordinateComparator`. `removeReadStubs` then removes
  reads whose actual remaining base-array length is **less than 10**. HC's
  later MQ/mate-contig filtering is not used by Mutect2.
* `ReadUtils.getSamplesFromHeader` returns a TreeSet, so sample order is sorted
  sample-name order, including samples with no surviving reads. Within each
  sample, splitting retains the trimmed coordinate order. Haplotype order is
  the trimmed `AssemblyResultSet` order. The dump preserves all three orders;
  no order-dependent reduction is replaced by a map traversal.
* PCR capping leaves the last indel-quality element untouched; minimum-quality
  flooring still visits it. Java compares quality bytes as **signed** after
  the unsigned BQ-to-MQ cap. The new somatic preparation entry point preserves
  this, including BQ128/254 when MQ255 is supplied in the independent fixture.
  Valid BI/BD FASTQ qualities are at most 93. No platform transcendental
  function is used.
* `AlignmentUtils.createReadAlignedToRef` assigns `HC:i` only when the chosen
  priority-adjusted best allele is informative. Choosing a higher-priority
  lower-likelihood haplotype may give negative confidence and hence no HC tag.
  Mate positions, insert size, flags, MAPQ and original alignment tags (e.g.
  NM/MD/AS/XS/MC/MQ) are **not recomputed** after read realignment.
* BI/BD tags are materialized by overlapping-pair adjustment, even when all
  their changed bases are later clipped away. A Boolean alongside stage-one
  read arrays preserves that fact; merely testing the final arrays for Q20
  would lose tags. The arrays written on aligned reads are the clipped,
  overlap-edited arrays, **not** the PCR-adjusted temporary HMM arrays.
* GATK's multi-input header merger rewrites PG/RG tag IDs. htsjdk
  `SamFileHeaderMerger` merges PG records by ancestry layers, grouping equal
  IDs/attributes in insertion order, remapping PP before the next layer. Its
  global counter increments even for IDs that do not collide; colliding IDs
  receive `.` plus uppercase base36(counter). On these inputs the tumour PG
  `MarkDuplicates` becomes `MarkDuplicates.P`. This is implemented from the
  source, not normalized away or hard-coded to this suffix.
* `writeBamOutput` occurs **after** genotyping. With
  `ALL_POSSIBLE_HAPLOTYPES`, the writer emits every surviving matrix read,
  sample order then read order, without deduplication across regions. The
  coordinate BAM writer sorts those records using the full SAM coordinate
  comparator; equal-comparator records retain their insertion order. The port
  projects them in this same file order.
* `AR:Z`, `CR:Z` and `XA:Z` are later-stage tags: `callMutations` calls
  `annotateReadLikelihoodsWithRegions` (AR/CR), then supported-allele annotation
  adds XA. AR/CR values are already fixed by the trimmed padded/callable spans, so
  both are included in the projection without invoking genotyping. Only XA
  depends on called alleles and is excluded. AR, CR, HC, BI, BD, PG, RG, AS, XS,
  MC, MD, MQ, NM and SA are compared wherever present. These tags are not cleared
  merely because they
  describe a former alignment. This refines §12's note that biological reads
  were outside stage two, without importing annotations into stage three.

### 13.2 Reuse and validation contract

HC's portable Rust prefix-sharing PairHMM is called directly with prepared
read arrays. It contains no nomics-simd/NEON call. The complete chr21 raw-matrix
bit dump is the evidence that this path matches Mutect2's Java settings.
HC's original preparation entry point keeps its allocation count and behavior;
a new entry point accepts BI/BD and uses Java signed-byte flooring. The PCR
cache and repeat finder are shared. HC's `likelihood::finish` supplies the
unchanged default capping/disqualification. Read clipping and coordinate
comparison became public; realignment gained a separate ln-aware entry point,
sharing SW/CIGAR composition/left alignment while leaving HC's log10 entry point
unchanged. The assembly callback is synchronous and preserves traversal timing.

The JVM oracle uses the unmodified jar's private preparation/raw-computation
methods by reflection, then its public normalization, filtering, log conversion
and realignment methods. It independently asserts that its intermediate-state
orchestration gives the same retained evidence and binary64 bits as the public
`computeReadLikelihoods` entry point. Raw, capped, filtered and ln matrices are
serialized as 16 lowercase hexadecimal raw binary64 bits in allele/read order.
Prepared bases/BQ/BI/BD/GCP and realigned names/flags/positions/CIGAR/HC are also
dumped. The independent tiny fixtures cover supplied indel qualities, PCR
repeats, floors/MQ caps, raw PairHMM/ln bits and ln tie-window endpoints.

### 13.3 Not supported / not confirmed

The existing two-indexed-BAM, one-normal/one-tumour, BED/default-settings
contract remains. No non-default likelihood/PCR/GCP/disqualification settings,
DRAGstr, native GKL, GPU/NEON, flow, homopolymer collapse, allele prefiltering,
forced alleles, mitochondria, gVCF, CRAM or multi-tumour modes are exposed.
Flow-order read groups are rejected explicitly. Source records are recovered
for tag projection by input-file ordinal and exact BGZF virtual offset, so
differing records with the same sample/name/flag remain distinct. Arbitrary JVM/graph identity-order limitations from §12.1
remain. No arbitrary-input, x86/native or somatic genotyping parity follows
from the five-contig gates. SW failure with offset -1 and exotic BAM auxiliary
array/hex types are not independently exercised by these goldens. The Java BAM
projection oracle fails explicitly on an unhandled tag value type. Rust's BAM
projection supports those value types, but independent evidence is still needed.


## 14. Partial stage-four implementation (2026-11-12)

Note (lead, 2026-11-13): this section records the port as it stood when the section was
written. The caller is complete: `nomics somatic mutect2` writes the VCF, stats and f1r2
archive (§16), and the gates are adopted in `docs/GATES.md`. Statements below that the
caller, the VCF or the f1r2 collector is missing, or that a gate is a proposal, are history.

**This is not the final Mutect2 caller.** `nomics somatic mutect2-genotyping`
is a diagnostic command for capped read/fragment likelihoods, event merging,
fragment marginalization, default tumour/normal/artifact LODs and the strict
emission tests. It writes no VCF, f1r2 archive or genotyping-dependent XA tags.
`nomics somatic mutect2` is not implemented. Consequently the required final
caller gate fails; counts of passing events are not record parity.

### 14.1 Confirmed source details and corrections

Pinned sources reused outside the repository: GATK 4.6.2.0 commit
`76edc75c26504da94bbaee66584e107e76ee15de`, htsjdk 4.2.0 and Commons Math 3.5.
The unmodified local jar first reproduced the chr21 VCF (header and all 38
records, excluding only GATKCommandLine), stats, regions and filter counts on
the Homebrew JDK 17. Only then was the instrumented genotyper trusted.

* The fragment fallback rules and the appended spanning-deletion matrix allele
  refine §§5/6 above. They are not output normalizations. The independent
  fragment fixture exercises <=2 flagged reads, >2 unflagged reads, a single
  clean read, all-flagged groups, collisions, treeification and resizing.
* The genotyper's second cap is calculated as
  `45 * (-Math.log(10)/10)`, then added to the best likelihood. Preserve that
  operation order; do not algebraically reduce it to `-4.5*Math.log(10)`.
* `MathUtils.sumArrayFunction` initializes its result from column zero, then
  adds following columns in order. It is undefined for a zero-column Java
  RealMatrix. The caller's evidence path returns 0 before invoking the engine
  for an empty matrix. The Rust numerical engine requires nonempty matrices;
  the caller wrapper implements the same empty-evidence branch.
* The variational prior, posterior and per-fragment contributions accumulate
  as `(priorContribution + posteriorContribution) + likelihoodContribution`.
  Only the per-fragment entropy reduction uses DoubleStream's compensated sum.
  Dirichlet and L1 sums, responsibilities and per-fragment outer sums are plain
  loops. The port reuses Java-exact Math.log/exp and Commons Gamma/digamma.
* `NaturalLogUtils.logSumExp` still adds positive zero when its scaled sum is
  one. The wording in §6.1 formerly omitted this signed-zero detail.
* The RNG expression is literally `MTAC.genotypeGermlineSites &&
  rng.nextDouble() < MTAC.genotypeGermlineSitesFraction`. The fixed default
  false short-circuits the only genotyper draw; no random value is consumed in
  this configuration. The Java diagnostic rejects a true germline-site switch.

Ordering: sample order remains the sorted header TreeSet order. Grouping uses
existing OpenJDK-17 String-key HashMap emulation, including computeIfAbsent head
insertion, stable resize splits and tree-bin root movement. Members are summed
in evidence order from positive zero, not stream order or max-over-mates.
Fragment intervals use the selected reads after realignment and span mate gaps.
Event starts are a sorted set, merging retains first-seen ALT order, and the
allele mapper uses LinkedHashMap order with the optional `*` appended last.
The existing HC event merging/mapping implementation is exposed through a
parameterized entry point; HC still includes spanning events and keeps MNP
max distance zero. No HC germline genotyper or diploid annotation structures
are used by the somatic evidence engine.

The code used here does not depend on JVM identity hashes: read names are
Strings; Haplotype.hashCode is Arrays.hashCode(bases); Event.hashCode uses
start and its alleles; SimpleAllele.hashCode uses bases and reference status.
Hash sets used to deduplicate events are not traversed to establish allele
order. Later phasing, annotation and f1r2 order are not proved by this partial
port. Existing assembly-stage identity-hash limitations remain in §12.1.

### 14.2 Validation scope and incomplete work

The pre-annotation internal dump contains all visited varying regions, including
non-emitting events: capped read matrices, fragment member and representative
indices, intervals, fragment matrices, per-event alleles and retained fragment
indices, per-sample allele matrices, TLOD/NLOD/NALOD as raw binary64 bits and
strict threshold decisions. The Java wrapper uses the actual unmodified
SomaticGenotypingEngine, asserts its intercepted somaticLogOdds calls equal
the independently projected ones, and writes a real VCF and stats for
qualification. The separate annotation dump serializes the actual returned
call attributes and genotype attributes before VCF formatting, with doubles
as raw bits. There is no Rust annotation dump or annotation parity claim yet.

`mutect2_calls_goldens.rs` runs all seven available contigs, compares stats by
bytes and checks passing-event counts. Qualified internal dumps are compared
by bytes when NOMICS_M2_CALLS_ORACLE points to their directory. Missing goldens
or missing internal dumps are explicitly reported. Neither passing counts nor
stats imply VCF/annotation/f1r2/XA parity. The pure-shell/Python final gate
`scripts/somatic_mutect2_gate.sh [CONTIG ...]` is provided now and deliberately
fails when the full caller cannot execute. No final output is substituted
from a reference VCF or reference tool.

**Not implemented / not confirmed in Rust:** population-resource annotation,
initial somatic genotype construction, common-prefix/suffix call trimming and
remapping, the full StandardMutectAnnotation contract, physical phasing,
ECNT/ECNTH, XA tags, VCF headers/formatting, BGZF output, index generation, f1r2
collection and archive output. There are no final-output byte comparisons or
independent annotation-statistic unit fixtures. No archive or index is written.
Only the inherited two-BAM/default tumour-normal contract is exposed for the
diagnostic; all non-default modes enumerated in §§11–13 remain unsupported.
Detailed commands, results, timings and explicit final-gate failures are
appended to reference/oracle/mutect2/VERIFICATION.md.

### 14.3 Additional archive-order discrepancy

The qualification run also exposed an output-order limit: its **unmodified**
jar writes `.ref_histogram, .alt_table, .alt_histogram`, whereas the chr21
container golden writes `.alt_histogram, .ref_histogram, .alt_table`. File names
are otherwise identical; by name, reference histogram and alt table bytes
match and alt histogram values/non-column lines match under the already
specified column-label rule. Source `F1R2CountsCollector.closeAndArchiveFiles`
passes unsorted `tmpDir.listFiles()` to `IOUtils.writeTarGz`, which preserves
that array's order. The source caution in §8.2 was therefore warranted: member
order depends on directory enumeration, beyond the known column hash order.
This is a local/container discrepancy, not proof of within-filesystem run
nondeterminism. It is documented in `docs/REFERENCE-VARIABILITY.md`.

No member-order normalization is proposed or silently applied. The required
exact-order gate is retained and rejects this local oracle archive at member
zero. No filesystem-order emulation or fixed golden order is substituted.
This reference-order discrepancy is reported to the lead along with the
incomplete Rust output implementation.

## 15. Read/fragment annotation continuation (partial)

Note (lead, 2026-11-13): this section records the port as it stood when the section was
written. The caller is complete: `nomics somatic mutect2` writes the VCF, stats and f1r2
archive (§16), and the gates are adopted in `docs/GATES.md`. Statements below that the
caller, the VCF or the f1r2 collector is missing, or that a gate is a proposal, are history.

This continuation does **not** finish the caller. The diagnostic
`mutect2-genotyping` retains its byte-identical pre-annotation output and adds
two sidecars by appending to `--output`: `.annotations.tsv` contains raw
attributes for a limited subset, with AF as binary64 bits;
`.annotation-fields.tsv` contains those same fields in VCF text form. They
are not VCFs. A third sidecar, `.phasing.tsv`, contains GT/PGT/PID/PS after
region-wide phasing. The annotation subset is exactly INFO
AS_SB_TABLE/DP/MBQ/MFRL/MMQ/MPOS and genotype AD/DP/AF/F1R2/F2R1/FAD/SB.
The integration test compares coordinates, emitted alleles, sample order and
every value of those fields against available raw goldens, and compares the
raw subset with the qualified jar dump when NOMICS_M2_CALLS_ORACLE is supplied.
The complete post-annotation gate remains incomplete.

### 15.1 Source refinements and independent fixtures

Pinned GATK/htsjdk source revisions remain those in §0, inspected outside the
repository. No source or reference tool is linked or called by nomics.

* `trimmedToUntrimmedAlleleMap` and the AD/DP marginalization maps are ordinary
  HashMaps collected in call-allele order. Matrix rows follow key iteration,
  not VCF allele order. `SimpleAllele.hashCode` is
  `31*Arrays.hashCode(bases)+Boolean.hashCode(isRef)`. The subset models list
  bins, putIfAbsent tail insertion and resize splitting; a bin needing
  tree-order handling is rejected explicitly rather than guessed. No JVM
  identity order was observed to reach the tested annotation values.
* Best-allele search chooses the first strict maximum. With a gap strictly
  below `Math.log(10)*0.2`, reference priority can replace the maximum; its
  resulting confidence may be negative. Counting still uses **confidence
  >0.2**, as §5 specifies. Dropped alts and an appended matrix `*` are removed
  before the annotation matrices are formed. The spanning row remains in
  upstream evidence, but is never an emitted allele in this configuration.
* `getEffectiveCounts` initializes from the first column and then adds following
  normalized columns. Empty evidence returns zero counts. Initial AF uses
  the all-ones prior for empty evidence, normalized in emitted allele order;
  hard AD/FAD do not replace AF. No new transcendental implementation is used.
* Per-allele medians exclude MQ0/255, but AD/DP/SB/FAD do not. AS_SB_TABLE pools
  a sample only when its informative read count is strictly >2. Orientation
  annotations test the first selected fragment read with BQ **>=20**, whereas
  collection in §8.2 uses **>20**. The two predicates remain separate.
* Base-quality lookup rejects deletion/skipped-reference operators. Read-position
  lookup uses the first read index after a deletion, even though a ReadUtils
  overload's prose says “preceding”. ReadPosition's outer alignment-overlap
  guard runs before the rank-sum helper's leading-insertion exception, so
  that exception does not supply an MPOS outside the aligned span. Terminal
  hard clips count as bases; soft clips participate in the index walk.
* Coverage stores INFO DP as a Java String internally, despite the Integer
  VCF declaration; genotype DP is an Integer. The raw subset preserves this
  distinction rather than treating internal attribute types as VCF types.
* `MathUtils.median(int[])` applies Commons Median then FastMath.round; even
  medians therefore round upward at .5, including negative values. Empty
  defaults remain MBQ/MFRL=0, MMQ=60, MPOS=50.

`MutectAnnotationNumerics.java` invokes the unmodified jar's private initial
genotype/count methods and public best-allele, coordinate and median helpers.
Its committed fixture contains 89 matrix cases (including five empty matrices,
explicit tie/informative endpoints and 13/20/49-allele hash-map resizes), nine
median cases, and 119 CIGAR-coordinate cases. Expected AF/soft-count/confidence
values are raw binary64 bit patterns produced by Java, not by Rust. The Java
fixture explicitly switches the likelihood object's scale flag before filling
its ln matrix; leaving that flag at the constructor's log10 default would test
the wrong tie window. Unit tests consume the fixture without requiring Java.

Physical phasing uses the source's original call order and set membership,
size and containment operations, so set iteration does not establish output
order. The called-haplotype union includes support for emitted REF as well as
ALTs, as in callMutations. Only biallelic concrete calls enter the alternate
haplotype map. `MutectPhasing.java` independently invokes the jar's private
constructPhaseSetMapping on 120 deterministic set graphs; expected group IDs
and slots are in `phasing.tsv`. The integration test compares the resulting
GT/PGT/PID/PS fields with the VCFs, including phased hom-ref normals. This
sidecar is separate from the raw annotation subset; the full raw internal
dump is still incomplete.

The default-false germline-site flag still short-circuits the only unseeded
genotyper draw in SomaticGenotypingEngine.callMutations (§14.1). The subset adds
no RNG. No annotation constant was changed to match a golden.

### 15.2 Archive proposals and remaining work

At the lead's explicit request, the stage-four proposal now compares archive
member **sets**, while retaining exact member content and the alt-histogram
column-label rule. Source LearnReadOrientationModel.doWork extracts inputs
separately and discovers each type by suffix; one member per type cannot
inherit tar enumeration order. A permutation of the chr21 golden gives the
same nomics prior bytes. This does not reorder the input archive list and
does not qualify a Rust collector. Intended fixed write order for a future
collector is ref_histogram, alt_histogram, alt_table. `docs/GATES.md` records
both proposed normalizations awaiting the lead.

Still absent: final `mutect2` CLI, VCF header/records/BGZF/index, population/PoN
annotations, ECNT/ECNTH, tandem repeats, f1r2
collection/archive writing, XA, and a full post-annotation raw-bit comparison.
Unrounded LODs remain gated at the existing pre-annotation boundary, not as
stored INFO arrays. All unsupported modes in §§11–14 remain unsupported.
No full-output gate or arbitrary-input parity is claimed for this subset.

## 16. Raw caller and supported-allele tags (record-gate continuation)

`nomics somatic mutect2` accepts the two-BAM sarek tumour-normal contract and
writes BGZF VCF plus `<output>.stats`. With `--f1r2-tar-gz` the collector of §8.2
(`mutect2_f1r2`) is the traversal's locus sink, so the archive comes from the same pass
over the reads as the calls. No index is written. Other modes remain unsupported as listed above. The default header
schema includes declarations even for annotations that are not emitted in this
mode. Contigs come from the reference index; samples come from the sorted BAM
sample-name set. INFO sorts by key; FORMAT puts GT first, then sorts remaining
keys. Numeric serialization reuses HC's Java-compatible VCF formatting; LODs
multiply by LOG10_E, and POPAF uses the existing Java-exact log10 port.

The caller reuses the existing annotation/phasing computations in memory; no
oracle files or external programs are read or invoked by the implementation.
The separate evidence diagnostic retains its output boundaries. The biological
read projection now includes XA, which requires genotyping even when running
`mutect2-likelihoods`.

Source details (GATK 4.6.2.0, `SomaticGenotypingEngine.callMutations`,
`getGermlineAltAlleleFrequencies`, `FeatureContext.getValues(input, int)`,
`TandemRepeat`, `GATKVariantContextUtils.getNumTandemRepeatUnits`, and
`AssemblyBasedCallerUtils.annotateReadLikelihoodsWithSupportedAlleles`):

* POPAF and PON query records **starting** at the original event position, not
  all resource records overlapping it. Alleles are padded to the longer REF
  before matching AF, and later matching resource records overwrite earlier AF.
  As in `determineReferenceAllele`, equal-length conflicting REFs fail; for
  unequal lengths, extension uses the longer REF suffix without an additional
  prefix-consistency check. Numeric AF text other than `.` must parse.
* ECNT counts the distinct minimal potential somatic events passing both LOD
  thresholds. The marginalization-only spanning `*` row is not a merged ALT
  and does not enter this set. It is not simply the number of VCF records.
* ECNTH counts potential somatic events on the ALT's best-supported haplotype.
  Every best-read assignment counts, including uninformative ones. Reference
  priority uses the actual reference haplotype index. Equal support preserves
  haplotype order; map iteration does not decide this tie.
* TandemRepeat first requires `VariantContext.isIndel()`: a mixed indel/MNP
  or indel/SNP record receives no STR/RU/RPA, even if its bases repeat. The
  computation excludes the leading anchor base, keeps the REF
  count from the first ALT and the repeat unit from the last ALT. A zero count
  aborts the whole annotation. Counts use signed arithmetic, as Java does.
* XA uses informative fragment assignments, the fragment's selected reads,
  and emitted allele indices at the trimmed call position. Existing XA text
  (including BWA alternative-alignments text) is retained and followed by
  `, ` before each new `contig:position=allele` entry. Tags are not recomputed
  from AD or from per-read assignments.

Ordering limits: existing String-key fragment HashMap emulation and assembly
limitations remain. Annotation SimpleAllele map list-bin order is unchanged.
It explicitly rejects a bin requiring treeification (at least nine colliding
keys after capacity reaches 64). OpenJDK tree tie-breaking for distinct,
non-comparable keys with equal full hashes can involve identity hashes and
change traversal order. No such order is guessed. The seven raw goldens have
at most five total alleles per record, so none reaches this rejection. Event
membership and support-count lookup do not expose hash iteration in the new
annotations; output INFO and FORMAT ordering is lexical.

Full-record, header, stats and XA measurements for this continuation are in the
appended verification report; the f1r2 archive's are in `docs/GATES.md`. Passing these
gates does not qualify non-default modes, arbitrary assembly identity-order cases, or GPU execution.
