# Porting nomics from GATK 4.6.2.0 to 4.0.12.0: source diff and effort

**Scope:** what changes when our exact GATK 4.6.2.0 port (MarkDuplicates, BaseRecalibrator + ApplyBQSR, `HaplotypeCaller -ERC GVCF`, GenotypeGVCFs) is re-pinned to **GATK 4.0.12.0**, the version used for the UK Biobank 500k WGS GATK release. This page covers GATK, Picard, htsjdk and GKL internals only. The pipeline side (BWA, samtools, CRAM, goldens) is in `port-ukb-pipeline-and-goldens.md`.

**Method.** I diffed sparse checkouts by class name, since `git diff` across the two tags is meaningless. The checkouts were: GATK `4.0.12.0` (commit `02682d54`) and `4.6.2.0` (`76edc75c`); Picard `2.18.16`, `2.18.26`, `2.20.3` and `3.4.0`; htsjdk `2.18.1` and `4.2.0`; GKL `0.8.6` and `0.8.11`. To date PRs I used a treeless clone of the GATK history (`git tag --contains`). It covers 1,077 commits between the two tags.

**Citations.** `old:` means `src/main/java/org/broadinstitute/hellbender/` at 4.0.12.0. `new:` means the same path at 4.6.2.0. "Lx" is a line number. **[unverified]** marks inference that I did not confirm in source or by running code.

## Bundled versions (build.gradle)

| | GATK 4.0.12.0 | GATK 4.6.2.0 |
|---|---|---|
| htsjdk | 2.18.1 (L59, forced at L151-152) | 4.2.0 (L62, `strictly`) |
| Picard | 2.18.16 (L60) | 3.4.0 (L63) |
| GKL | 0.8.6 (L313) | 0.8.11 (L75) |
| Java | 8 (L58) | 17 |

**UK Biobank did not use the bundled Picard.** deCODE ran Picard MarkDuplicates **2.20.3** (built against htsjdk 2.19.0) in queryname mode, per read group. Sanger ran **2.18.26** (htsjdk 2.18.2) in queryname mode, per sample, with `OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500`. samblaster and sambamba only added mate tags or sorted. BQSR, ApplyBQSR (with SQQ 10/20/30) and `HaplotypeCaller -ERC GVCF` were GATK 4.0.12.0. (Coordinator's reading of the Nature 2025 SI; I did not re-verify it.) So the MarkDuplicates port to size is **Picard 2.18.26/2.20.3 in queryname mode**, not the bundled 2.18.16.

## Headline findings

1. **The QUAL model is the old exact AF calculator.** `use-new-qual-calculator` defaults to `false` at old: `tools/walkers/genotyper/GenotypeCalculationArgumentCollection.java` L51, and to `true` (deprecated flag) at new: L66-67. The switch is commit `6faaa03` "Make new qual the default in HaplotypeCaller and GenotypeGVCFs (#5484)", dated 2019-01-28, first in **4.1.0.0** and not an ancestor of 4.0.12.0.
   - The whole old calculator family was deleted by 4.6.2.0: `AFCalculator`, `ExactAFCalculator`, `IndependentAllelesDiploidExactAFCalculator`, `ReferenceDiploidExactAFCalculator`, `StateTracker`, `AFCalculatorProvider` and the rest, about 1,500 lines. HC uses `FixedAFCalculatorProvider` (old: `HaplotypeCallerEngine.java` L203). The default is `EXACT_INDEPENDENT` (old: `afcalc/AFCalculatorImplementation.java` L20, L56).
   - The calculator is also used inside **isActive** (old: `HaplotypeCallerEngine.java` L303-312). PR #6343 "single sample fast mode for qual in HC isActive" (4.1.5.0) replaced that, so active-region boundaries depend on it too.
   - It uses `MathUtils.approximateLog10SumLog10`, the Jacobian table (old: `IndependentAllelesDiploidExactAFCalculator.java` L388-389; `ReferenceDiploidExactAFCalculator.java` L155-156). That code is byte-identical between tags, and we already have the table.
2. **4.0.12.0 HaplotypeCaller is not fully deterministic, and this is the biggest risk.** `SeqVertex.hashCode()` is `System.identityHashCode(this)` (old: `graphs/SeqVertex.java` L59-64). `BaseGraph.getSources/getSinks/outgoingVerticesOf/incomingVerticesOf` collect into `Collectors.toSet()`, a HashSet ordered by identity hash (old: `graphs/BaseGraph.java` L81-90, L334, L344).
   - GATK fixed this in #6195 "Fix non-determinism in KbestHaplotypeFinding/SeqGraphZippingCode" (4.1.4.0), which switched to `LinkedHashSet`. A related priority-queue tie fix is #6104 (4.1.4.0).
   - Wherever the iteration order reaches graph zipping or path enumeration, 4.0.12.0 output depends on JVM identity hashes, not on the inputs. Our rule ("source is the spec") cannot be met by construction there.
   - Options: (a) measure the run-to-run variation of the 4.0.12.0 reference (as in `REFERENCE-VARIABILITY.md`) and document the residual; or (b) emulate HotSpot's identity-hash generator and its exact call sequence. Option (b) is **[unverified]** and very likely impractical.
3. **4.0.12.0 HaplotypeCaller is mostly the code that 4.6.2.0 rewrote.** The rewrites that change default `-ERC GVCF` output, all after 4.0.12.0:
   - KBestHaplotypeFinder: recursion replaced by Dijkstra (#5462, 4.1.0.0), plus the lost-bestPaths fix (#5952).
   - Symmetric likelihood capping (#6196, 4.1.4.1).
   - Unfilled likelihoods abolished and the annotator revamped (#6172).
   - Read selection for PL/AD/DP (#6055).
   - Window and padding refactor (#6358).
   - leftAlignIndels rewrite (#6427).
   - Realignment to the best haplotype (#6461, #6498).
   - CigarBuilder indel order (#6510).
   - Overlapping-pair quality fixes (#5926, #6815, #6943, #6886).
   - Spanning-deletion QUAL fix (#6859).
   - Phasing fixes (#7019, #6937).
   - Minimal-representation and no-data hom-ref fixes (#8567, #8715, #7644).
4. **Smith-Waterman defaults to JAVA in 4.0.12.0** (old: `AssemblyBasedCallerArgumentCollection.java` L128) and to `FASTEST_AVAILABLE`, which means GKL AVX on x86, in 4.6.2.0 (new: L169). So **GKL Smith-Waterman is out of scope for 4.0.12.0**. GKL matters only for PairHMM and the deflater.
5. **BQSR and ApplyBQSR have the same defaults, table format, covariates and quantizer, but two changes alter the output:**
   - the empirical-quality posterior in `RecalDatum` was rewritten;
   - the known-site masking at deletions differs.

   SQQ handling is identical (section 3).
6. **On x86, GATK uses GKL by default in both versions:** the PairHMM is `FASTEST_AVAILABLE`, which picks AVX/OMP, and the deflater is Intel. If UKB ran with defaults on x86 (likely, **[unverified]**), its gVCFs are **GKL-mode** output, and a Java-mode port alone will not match them.

## 1. HC traversal, assembly, Smith-Waterman, PairHMM, downsampling

**Argument defaults that are unchanged** (diffed every `@Argument` in the HC, assembly, likelihood, native-PairHMM and region collections):
- min/max region 50/300, padding 100, `max-reads-per-alignment-start` 50, active-prob 0.002, propagation 50 (old: `HaplotypeCaller.java` L142-147; new: `engine/spark/AssemblyRegionArgumentCollection.java` L22-27).
- kmers [10,25], `min-pruning` 2, `max-num-haplotypes-in-population` 128, `min-dangling-branch-length` 4, adaptive pruning off, recover-dangling-heads off.
- `pcr-indel-model` CONSERVATIVE (old: `LikelihoodEngineArgumentCollection.java` L50; new: L83), `pair-hmm-gap-continuation-penalty` 10, `phred-scaled-global-read-mismapping-rate` 45, `pair-hmm-implementation` FASTEST_AVAILABLE.
- `--native-pair-hmm-threads` exists in both, default 4, with use-double-precision false (old: `PairHMMNativeArgumentCollection.java` L11-16).
- `GVCFGQBands` 1..60, 70, 80, 90, 99 (old: `HaplotypeCallerArgumentCollection.java` L80-86; new: L117-123). `indelSizeToEliminateInRefModel` 10.
- The HC read filter list is identical: MQ≥20, MQ available, mapped, primary, not duplicate, vendor-QC pass, non-zero reference length, good CIGAR, wellformed (old: `HaplotypeCallerEngine.java` L319-331; new: L463-475). The filter classes have 0 diff lines.
- Overlapping-pair correction is **on** in both. Old L131 names the field `doNotCorrectOverlappingBaseQualities=false` and the call site negates it (old: `HaplotypeCallerEngine.java` L681). #6528 only renamed the flag.
- Smith-Waterman parameters are identical: haplotype-to-reference `NEW_SW_PARAMETERS` (200,-150,-260,-11), read-to-haplotype (10,-15,-30,-5), dangling ends `STANDARD_NGS` (25,-50,-110,-6) (old: `utils/read/CigarUtils.java` L19, L27, `readthreading/ReadThreadingGraph.java` L738, L769; new: `AssemblyBasedCallerArgumentCollection.java` L61-69).
- RNG seed 47382911 (old: `utils/Utils.java` L57-58; new: L52-53).

**Behaviour changes that alter output:**

| Component | Change | Old | New | Port size |
|---|---|---|---|---|
| Graph iteration order | Identity-hash HashSets (see finding 2) | `graphs/BaseGraph.java` L81-90, L334, L344; `SeqVertex.java` L59-64 | LinkedHashSet (#6195) | Exact port impossible in places; measure the residual |
| KBestHaplotypeFinder | Recursive finder (`RecursiveSubHaplotypeFinder`, `AggregatedSubHaplotypeFinder` with PriorityQueue, `KBestHaplotype.compareTo` score only) vs Dijkstra with a path tie-break | `graphs/KBestHaplotypeFinder.java` (424 lines) + 4 deleted files | `graphs/KBestHaplotypeFinder.java` | L (5-10 d): new code; tie order is undocumented |
| Read-threading graph, dangling ends, SeqGraph cleanup | Restructured into `AbstractReadThreadingGraph`; dangling-head indels added (#6113), cycle and merge fixes (#5960, #7086) | `readthreading/ReadThreadingGraph.java` (1,138 diff lines), `ReadThreadingAssembler.java` L495-508 | `AbstractReadThreadingGraph.java`, `ReadThreadingAssembler.java` | L (5-8 d) |
| Region trimming and padding | #6358 refactor (4.1.5.0). 4.6.2.0 keeps `trimLegacy` behind `--enable-legacy-assembly-region-trimming` (default false) | `AssemblyRegionTrimmer.java`, `AssemblyRegionTrimmerArgumentCollection.java` (indel padding 150, SNP 20, max extension 25) | `AssemblyRegionTrimmer.java` L169-170, L218; `engine/spark/AssemblyRegionArgumentCollection.java` L101-128 (indel 75) | M (3-6 d) if `trimLegacy` reproduces 4.0.12.0 **[unverified]** |
| Likelihood capping | Cap relative to the best *alt* allele, vs symmetric (#6196) | `utils/genotyper/ReadLikelihoods.java` L367-383 | `utils/genotyper/AlleleLikelihoods.java` L441-458; flag `LikelihoodEngineArgumentCollection.java` L105 | S (<1 d): the flag already exists |
| Isactive | Uses the old AF-calc genotyper; #6343 fast path absent | `HaplotypeCallerEngine.java` L303-312, L464-479 | `HaplotypeCallerEngine.java` L440-453, L674-693 | M (2-3 d), on top of the AF calculator |
| CIGAR, left-align, realign | leftAlignIndels rewrite (#6427, #6792), realign to best haplotype (#6461, #6498), CigarBuilder order (#6510), hardClipCigar (#6280) | `utils/read/AlignmentUtils.java` (872 diff lines), `CigarUtils.java` (576) | same files + `CigarBuilder.java` (new) | L (5-10 d) |
| Overlapping-pair quality adjustment | `FragmentUtils` fixes (#5926 disjoint pairs, #6815 insertion vs deletion quals, #6943, #6886) | `utils/fragments/FragmentUtils.java` | same file (128 diff lines) | M (2-4 d) |
| Downsampling | `ReservoirDownsampler` end-of-stream semantics (#5594, 4.1.0.0). The `nextInt` draw site is the same | `utils/downsampling/ReservoirDownsampler.java` L103, L133-171 | L119-123, L133-218 | S-M (1-2 d); it shifts RNG use if release timing changes **[unverified]** |
| Java Smith-Waterman | `consolidateCigar` vs skipping zero-length segments | `utils/smithwaterman/SmithWatermanJavaAligner.java` L79-81, L374-375 | L83, L338-340, L378 | none; equivalent **[unverified edge cases]** |
| Java PairHMM | `LoglessPairHMM`, `PairHMMModel`: 0 diff lines | | | none |

## 2. Genotyping, QUAL, annotations, ReferenceConfidenceModel, GVCF writer, GenotypeGVCFs

- **Calling threshold.** In `-ERC GVCF` both versions force the call threshold to `-0.0` and `annotateAllSitesWithPLs=true` (old: `HaplotypeCallerEngine.java` L245-256; new: L387-393). The 10 (old: `GenotypeCalculationArgumentCollection.java` L114) vs 30 (new: L28) default therefore does not matter for GVCF. It does matter for VCF mode and GenotypeGVCFs.
- **Old QUAL path** (old: `genotyper/GenotypingEngine.java` L225-317): `afCalculatorProvider.getInstance(...).getLog10PNonRef(reducedVC, ploidy, maxAlt, priors)`, then `calculateOutputAlleleSubset` keeps alleles by `isPolymorphicPhredScaledQual(allele, STANDARD_CONFIDENCE_FOR_CALLING)` (L389), MLEAC/MLEAF from the AF result, and `log10Confidence` from `getLog10PosteriorOfAFEq0`. **Size: L (5-8 d)**: about 1,500 lines of FP-order-sensitive Java plus `AFCalculationResult`.
- **Genotype likelihoods.**
  - Old: each read value for a homozygous genotype is `lk + log10(ploidy)`, summed, then `- readCount*log10(ploidy)` (old: `genotyper/GenotypeLikelihoodCalculator.java` L276-285, L392-430).
  - New: `MathUtils.sum(lk)` directly for homozygous genotypes (new: L80-86).
  - These differ at ulp level. #6351 (4.5.0.0) says "no change in output", but ulp changes can flip PL rounding. **S (1 d)**, plus verification.
- **Read sets for genotyping and annotation.** Old marginalizes with an `ALLELE_EXTENSION` interval (old: `HaplotypeCallerGenotypingEngine.java` L171, L185). New uses `retainEvidence` with an overlap predicate, and appends filtered reads with likelihood 0 (new: L192-197, L586-618). The changes come from #6055 and #6172. **M (2-4 d)**. The annotation read set feeds DP, RAW_MQandDP and the rank sums.
- **`qualToErrorProbLog10`.** Old is `qual / -10.0`, new is `qual * -0.1` (old: `utils/QualityUtils.java` L190-192; new: L188-190). They differ by 1 ulp for 91 of the 256 integer quals (q = 3, 6, 7, 12, …).
  - It is used by the `ReferenceConfidenceModel` ref-vs-any likelihoods (new: L479, L483) and the global mismapping rate.
  - The PairHMM uses the unchanged `qualToErrorProb` (`Math.pow(10, q/-10.0)`).
  - For HC with UKB's SQQ 10/20/30 binned quals the change is moot: those values and 2 are not affected, and the Q45/Q30 constants are equal.
  - It is **not** moot for BaseRecalibrator, whose old posterior evaluates it at every bin from 0 to 93 (section 3). **S** to port; it must be used in both places.
- **ReferenceConfidenceModel.**
  - The constants are equal: deletion qual 30, base-qual threshold 6, HQ soft-clip threshold 28, at most 40 indel-informative reads, indel error rate -4.5 (old L57-104; new L81-128).
  - What changed: the informative-read computation (#5469, #5607), `getDeletionQuality` with a realigned-reads flag (new L416), and soft-clip "os"/"oe" tags. **M (2-4 d).**
- **GVCF writer.**
  - Same bands.
  - New also folds no-call genotypes with `PL[0]==0` into reference blocks; old does not (old: `utils/variant/writers/GVCFBlockCombiner.java` L192; new: L195-197).
  - `HomRefBlock` GQ: old uses `calculateGQFromPLs(minPLs)`; new uses `getMinGQ()`, which is the same value when PLs exist (new L54-68, L150-156).
  - Block ploidy comes from the default ploidy in old and from the genotype in new; these are equal for diploid samples. **S (1-2 d).**
- **Default annotations in GVCF mode** are the same at both tags: `StandardAnnotation` + `StandardHCAnnotation`, then `filterReferenceConfidenceAnnotations` adds StrandBiasBySample and removes ChromosomeCounts, FisherStrand, StrandOddsRatio and QualByDepth (old: `HaplotypeCallerEngine.java` L168-182; new: L264-278). The emitted set is DP, RAW_MQandDP, BaseQRankSum, MQRankSum, ReadPosRankSum, ExcessHet, InbreedingCoeff (header only), and FORMAT AD/DP/GQ/PL/SB. Consequence: no QD jitter in GVCF mode, so in GVCF mode the RNG is consumed only by downsampling.
- **RAW_MQandDP.** Both write `RAW_MQandDP=%d,%d` (Σ MAPQ², n over reads with MAPQ≠255) (old: `annotator/RMSMappingQuality.java` L60, L115-117, L182-198; new: L58, L139, L213-229). The legacy `RAW_MQ` is only read for compatibility. Unchanged, except for which reads the matrix holds (above).
- **Rank sums change.**
  - ReadPosRankSum. Old: `getReadCoordinateForReferenceCoordinate(..., RIGHT_TAIL)`; a deletion becomes INVALID and is excluded; position computed by `calcAlignmentByteArrayOffset`, folded at `numOriginalBases/2`; usable if `softEnd ≥ refLoc` (old: `annotator/ReadPosRankSumTest.java` L44-81). New: `getReadIndexForReferenceCoordinate`, `min(left, right)` distance, a special case for an insertion at `vc.end+1`, and a two-sided overlap test (new: L42-80).
  - BaseQRankSum. Old: `getReadCoordinateForReferenceCoordinateUpToEndOfRead`, which returns a base even inside deletions. New: empty when the position falls in a deletion (old: `BaseQualityRankSumTest.java` L41-43; new: L41-43).
  - `MannWhitneyU` is the same algorithm.
  - **M (2-4 d)**: the old `ReadUtils` coordinate helpers are needed.
- **ExcessHet changes on every record.** 4.0.12.0 uses a mid-p exact test (`rareCopies<=1 → 0.5`, divide by `2·sum`, no clamp). For one diploid sample it gives `ExcessHet=3.0103` on every GVCF variant record, where 4.6.2.0 gives `0.0000` (old: `annotator/ExcessHet.java` L83-198; new: L137-207; #7394, 4.2.4.0). **S.** (Parallel read; I confirmed the old mid-p branch at old L133-145.)
- **GT, PL and GQ on GVCF variant records.**
  - 4.6.2.0 turns a best genotype containing `<NON_REF>` into `0/0` with all-zero PLs (#6437, 4.1.5.0), and some low-confidence hom-refs into `./.` (#8741). 4.0.12.0 emits, for example, `0/2` (NON_REF) as is (old: `utils/variant/GATKVariantContextUtils.java` L248-290; new: L310-370).
  - When likelihoods are flat, 4.0.12.0 also drops PL and GQ (old: `genotyper/AlleleSubsettingUtils.java` L82-93).
  - **S-M (1-3 d).** (Parallel read, not re-checked.)
- **Old allele plausibility and `*` handling.** `isPolymorphicPhredScaledQual`, and the upstream-deletion bookkeeping, which counts `<NON_REF>` (old: `GenotypingEngine.java` L389-449). This is folded into the AF-calculator item.
- **Other annotations.** Coverage, DepthPerAlleleBySample, DepthPerSampleHC, StrandBiasBySample and ExcessHet differ only in refactoring for HC, because HC always has filled likelihoods.
- **GenotypeGVCFs.**
  - 4.0.12.0 is a single class: `tools/walkers/GenotypeGVCFs.java`, with `GeneralPloidyFailOverAFCalculatorProvider` (old L186).
  - New is `GenotypeGVCFs` + `GenotypeGVCFsEngine`, plus #6437 (no-call NON_REF, PLs zeroed), #5835, #6401 and #8318.
  - UKB joint-called with GraphTyper, so this is **optional**: L (5-10 d).

## 3. BQSR, ApplyBQSR, GATKReport

- **Arguments and covariates.** `RecalibrationArgumentCollection` has 0 diff lines, so every BaseRecalibrator default is unchanged, as is the Arguments table that the report writes. The default covariates are ReadGroup, QualityScore, Context (mismatch 2, indel 3, low-quality tail 2) and Cycle (max 500), in the same order (`StandardCovariateList`: 15 lines of refactor).
- **Read filters** are identical (old: `bqsr/BaseRecalibrator.java` L150-161; new: L154-165).
- **The empirical-quality posterior was rewritten. It affects output: M (3-5 d).**
  - 4.0.12.0 works in log10 (old: `utils/recalibration/RecalDatum.java` L256-316).
    - Prior: `Math.log10(new Gaussian(0.9, 0, 0.5).value(d))`, cached for `d = min(|(int)(Qemp−Qrep)|, 93)`. A value of −∞ becomes `−Double.MAX_VALUE`.
    - Likelihood: `MathUtils.log10BinomialProbability(n, k, qualToErrorProbLog10(Qemp))`. This uses the **old `q/-10.0`** form and `log10Factorial` (a `Math.log10` cumulative cache below 10,000, commons-math `Gamma.logGamma` above).
  - 4.6.2.0 uses `NormalDistribution.logDensity` and a natural-log saddle-point binomial (new L288-367).
  - The argmax is mathematically the same. It differs in two cases:
    - large deltas, where `exp(−2d²)` underflows to 0 around d ≥ 20, so those bins are forbidden in 4.0.12.0 (`FastMath` behaviour **[unverified]**);
    - last-ulp ties.
  - Smoothing, the cap of 93, the 61-bin range, the first-max tie-break and the >2³¹ scaling of n are unchanged.
  - (BQSR finding by a parallel read, spot-checked here.)
- **Known-site mask: S.** 4.0.12.0 maps both ends of a site with `getReadCoordinateForReferenceCoordinate(…, LEFT_TAIL, true)` (old: `BaseRecalibrationEngine.java` L329-356). When a site's end falls inside a read deletion, 4.0.12.0 does not mask the base after the deletion; 4.6.2.0 does (new L353-377). Our `bqsr.md` §2.2 text is 4.6-only. This is common at known indels in WES.
- **CIGAR consolidation.** 4.0.12.0 uses `AlignmentUtils.consolidateCigar`, not `CigarBuilder`. That only matters for CIGARs BWA does not emit. **S.**
- **Rest of the recal table arithmetic is unchanged.**
  - Renames such as `estimatedQReported` becoming `reportedQuality`.
  - `roundToNDecimalPlaces` is identical (old: `MathUtils.java` L1375-1379; new: L927-931).
  - The table collapse is split into two methods but does the same arithmetic (old: `BaseRecalibrationEngine.java` L160-205; new: L161-213).
- **GATKReport.** `GATKReport`, `GATKReportColumnFormat` and `GATKReportVersion` have 0 diff lines (v1.1). `GATKReportTable`/`GATKReportColumn` now print non-finite doubles with `toString`, which is unreachable here (new: `GATKReportTable.java` L562-568).
- **ApplyBQSR** (`transformers/BQSRReadTransformer.java`):
  - Recalibration math and `getRecalibratedQual` bounding are the same.
  - A read group missing from the table: old returns the read unchanged (old L159-160); new throws unless `--allow-missing-read-groups` (new L171-185, #9020). It cannot occur when one table covers the BAM.
- **SQQ.** `constructStaticQuantizedMapping` is identical apart from a comment (old L232-…; new L292-…). `--static-quantized-quals`, `--round-down-quantized` (false), `--preserve-qscores-less-than` (6), `--emit-original-quals` (false) and `--global-qscore-prior` (-1) have the same defaults. ApplyBQSR binds `ApplyBQSRArgumentCollection` directly (old: `bqsr/ApplyBQSR.java` L94, L103).
  - Old `ApplyBQSRUniqueArgumentCollection.toApplyBQSRArgumentCollection` drops `staticQuantizationQuals` and `roundDown` (old L61-68; fixed by #5949 in 4.1.3.0). That method is used only by pipelines that convert collections (Spark), **not by ApplyBQSR**, so SQQ works in 4.0.12.0 ApplyBQSR. **[unverified]**: I did not grep all callers.
- **ContextCovariate.** Non-ACGT, non-N bases (IUPAC codes) break the context in new (#6625, 4.1.9.0) but not in old (old: `ContextCovariate.java`, `while (bases[..] != 'N')`). BWA output has only ACGTN, so this is S.
- **Output BAM.** Only the @PG line changes (version and CL). BAM encoding is in section 4.
- **Size:** M (4-8 d), mostly the posterior and the known-site mask, plus verification.

## 4. Picard MarkDuplicates, htsjdk, GKL

### MarkDuplicates (Picard 2.18.26 / 2.20.3 vs our 3.4.0)

- **Dispatch.** At both GATK tags, `gatk MarkDuplicates` runs Picard's class from the bundled Picard jar **[unverified: dispatch path not traced]**. There is no hellbender MarkDuplicates in either tree.
- **Duplicate flagging in queryname mode** is the same code in 2.18.16, 2.18.26, 2.20.3 and 3.4.0: `isDuplicate = index == nextDuplicateIndex || (queryname && index > next && name.equals(duplicateQueryName))` (2.20.3 L349-353, L434-436; 3.4.0 L345-349, L410-412). In this mode secondary, supplementary and unmapped-mate records inherit their template's flag. **Our port is coordinate-only, so queryname input is new work:**
  - building ReadEnds in queryname order;
  - `MemoryBasedReadEndsForMarkDuplicatesMap` for queryname in 3.4.0 versus the disk-based map in 2.x, which changes memory use, not output **[unverified]**;
  - representative-read and DT/DS tagging (3.4.0 changed the `TAG_DUPLICATE_SET_MEMBERS` logic; it is off by default).

  **M-L (4-8 d).**
- **Mate-pairing key.** 2.x uses the literal `"RG" + readName`; 3.4.0 uses `readGroupId + readName` (2.20.3 L576; 3.4.0 L558). This only matters when read names collide across read groups, which does not happen with Illumina names. **S.**
- **Optical duplicate finder.** It is identical between 2.18.26 and 3.4.0 (0 diff lines; graph-based for sets of 3-4 or more). The bundled 2.18.16 is the older pairwise version. Optical calls affect only the metrics file, and DT tags if tagging is on. Sanger's distance of 2500 is just an argument. **S.**
- **Scoring.** `DuplicateScoringStrategy` (htsjdk) is identical between htsjdk 2.18.1 and 4.2.0 (0 diff lines). htsjdk 2.18.2 and 2.19.0 were not checked **[unverified]**.
- **Metrics file** (from a parallel read of 2.18.16 vs 3.4.0; I did not re-check it):
  - 2.18.16 lacks the `all_sets`/`optical_sets`/`non_optical_sets` histograms and the `CoverageMult` value label.
  - It has no zero-denominator guard on PERCENT_DUPLICATION.
  - It orders libraries by HashMap rather than TreeMap.
  - `@PG CL` comes from Barclay 2.1.0 vs 5.0.0 and is **unknown**.
  - Check the same points for 2.18.26 and 2.20.3. **S-M (1-2 d).**

### htsjdk 2.18.1 vs 4.2.0 (BAM and VCF bytes)

- **BGZF.** `deflateBlock` is equivalent: same reset, level-0 fallback and block layout. Old keeps a dead partial-block branch (old: `util/BlockCompressedOutputStream.java` L366-405; new: L382-420). `DeflaterFactory`, `BlockCompressedStreamConstants` and `BinaryCodec` have 0 diff lines.
- **BAM bin.** Old reuses a record's existing bin (`computeIndexingBinIfAbsent`); new always recomputes (old: `BAMRecordCodec.java` L140-144, `SAMRecord.java` L1601-1607; new: L143-148). Both write 0 for an unplaced read. This differs only if an input bin was non-canonical. **S.**
- **VCF.** `VCFEncoder` differs only in the `Locale.US` formatting and braces. `formatVCFDouble` is the same thresholds and formats (2.18.1 L244-262; 4.2.0 L261-279). Header-line *text* (GATK `GATKVCFHeaderLines`, `##source`, `##GATKCommandLine`) differs by version and needs a GATES normalization or a per-version header table. A parallel read found two more htsjdk differences (I did not re-check them):
  - 4.2.0 quotes `Version="…"` in `##GATKCommandLine`; 2.18.1 does not (`VCFHeaderLine.java` L160-164).
  - 2.18.1 expands a missing FORMAT value for a Number>1 key to `.,.` (`VCFEncoder.java` L341-357), which does not matter for single-sample gVCFs.

  **S (1 d).**
- **Launcher.** Both scripts set `-Dsamjdk.compression_level=2` and async-write on (launcher L37-46).

### GKL (0.8.6 in 4.0.12.0; 0.8.11 in 4.6.2.0)

**What is on by default on x86:**
- **PairHMM.** `FASTEST_AVAILABLE` → `VectorLoglessPairHMM` OMP, then AVX, with 4 native threads (old: `utils/pairhmm/PairHMM.java` L81-106).
- **Deflater and inflater.** `IntelDeflaterFactory`/`IntelInflaterFactory` unless `--use-jdk-deflater/inflater` (old: `cmdline/CommandLineProgram.java` L77-81, L176-180; new: L78-82, L187-190).
- **GKL Smith-Waterman** is used only by 4.6.2.0 (finding 4).

**The PairHMM kernel:**
- **Float first, double fallback.** Float is computed first; if `result_float < MIN_ACCEPTED = 1e-28f`, it recomputes in double. The result is `log10(result_double) - LOG10_INITIAL_CONSTANT` or `log10f(result_float) - LOG10_INITIAL_CONSTANT` (`IntelPairHmm.cc` L124-134 at 0.8.6, L156-164 at 0.8.11; `pairhmm_common.h` L16).
- **Scaling and FTZ.** The initial constant is 2^120 in float and 2^1020 in double (`Context.h` L109-110, L149-150). FTZ is turned on with `_MM_SET_FLUSH_ZERO_MODE` (`IntelPairHmm.cc` L66-70); DAZ is not set.
- **Reduction order.** It keeps separate running sums ΣM and ΣX over the last row, then adds them (`avx-pairhmm-template.h` L308-343). Java instead adds `M+X` per column.
- **Compiler flags.** `-O3 -mavx` (avx512: `-mavx512f/dq/vl`); no `-mfma` or `-ffast-math` (`CMakeLists.txt` L10; `pairhmm/CMakeLists.txt` L7-8). Whether `-mavx512f` lets GCC contract into FMA is **[unverified]**.
- **0.8.6 → 0.8.11 changes output.**
  - `Context.h` replaced `LN10 = log(10)` / `INV_LN10 = 1/LN10` with truncated literals `2.302585` / `0.434294`. These enter the `matchToMatchProb` table (`log1p(...) * INV_LN10`, `Context.h` L58-59).
  - The `isinf(...) == -1` check also changed.
  - So 4.0.12.0 on x86 and 4.6.2.0 on x86 produce different PairHMM likelihoods.
  - The kernel files themselves differ only by casts (`(NUMBER)(tc->haplen)`).

**Deflater.** Levels 1 and 2 go to **ISA-L igzip stateless deflate**, a patched vendored copy; other levels go to Intel-patched zlib (`IntelDeflater.cc` L87-103 at 0.8.6; L81, L247 at 0.8.11). On x86, both the BAMs and `.g.vcf.gz` at level 2 are ISA-L bytes, not zlib bytes.

**GKL 0.8.6 deflater bug.** On every reuse, `resetNative` calls `isal_deflate_stateless_init` without restoring `level` and `level_buf` (`IntelDeflater.cc` L108-111). That function sets `level=0` and `level_buf=NULL` (vendored `igzip.c` L1004-1006). htsjdk resets the deflater before every block. So under 4.0.12.0 on x86, **block 1 of each BGZF stream is ISA-L level 2 and every later block is ISA-L level 0** (source-read; confirm on one real file). 0.8.11 restores the level, and the vendored ISA-L version also differs.

**Feasibility of GKL-exact output:**
- **PairHMM is exact by construction in principle.** All operations are IEEE single or double add, mul and div with no FMA (AVX path). A scalar Rust emulation is feasible if it does three things:
  - flush subnormal *results* to zero with x86 semantics (tininess detected after rounding **[unverified]**), with no DAZ;
  - reproduce `Context` table construction, which uses double `pow`/`log10`/`log1p` from glibc and then casts to float;
  - use the ΣM, ΣX order and the 1e-28f switch.
- **Unknowns.**
  - Whether glibc's `log10f`/`log10`/`pow` on the UKB hosts are correctly rounded; GKL links the host libm dynamically.
  - Whether AVX and AVX-512 kernels give identical per-cell values. They should, because only the lane width differs, **[unverified]**.
  - ARM NEON FZ also flushes inputs (DAZ-like), so the Spark needs software FTZ, not hardware FTZ.
  - GKL sets FTZ on the calling JVM thread and never restores it, so later *Java* FP work on that thread may also run with FTZ on **[unverified]**. That would be a GKL-mode effect outside the PairHMM.
  - The shipped AVX-512 `.so` must be checked for `vfmadd` with `objdump`.
- **Verification** needs an x86 box running GKL 0.8.6 to dump per-read×haplotype likelihoods, using the `--pair-hmm-results-file` idea, which exists only in 4.6.2.0, or a JNI harness.
- **ISA-L level-2 deflate** is deterministic but is a different algorithm from zlib (hashing, match finder, Huffman). Byte-exact BAM or gVCF.gz output would mean porting ISA-L's level-2 path, which is large. I recommend gating on decompressed content instead, because the UKB boundary is the gVCF records.

## Unknowns and risks

- **Non-determinism in 4.0.12.0** (finding 2): how often it changes records is unmeasured. Run 4.0.12.0 twice on the tier0 slice first.
- **JDK 8 vs 17 transcendentals.** Java-mode equivalence currently rests on `-XX:-UseLibmIntrinsic` under JDK 17 (`port-ukb-pipeline-and-goldens.md`). JDK 8 x86 intrinsics are unmeasured.
- **Whether UKB ran GKL PairHMM**, i.e. x86 with defaults. This decides whether the Java-mode port has any UKB golden to match.
- **Diff coverage.** Large diffs that I could not read line by line and only sized by PR history: `GATKVariantContextUtils` (1,294 diff lines), `AssemblyBasedCallerUtils` (913), `ClippingOp` (503), `EventMap` (321), `AlleleSubsettingUtils` (324), `ReadUtils` (429).

## Effort (engineer-days, given working 4.6.2.0 code and specs)

Assumptions:
- One engineer who knows our codebase.
- x86 4.0.12.0 goldens exist (Java mode and GKL mode), per `port-ukb-pipeline-and-goldens.md`.
- "Done" means record-identical, or a documented residual for the identity-hash non-determinism.
- Each item includes its spec and gate updates.

### (a) Java-mode 4.0.12.0 port

| Component | Low | High | Notes |
|---|---:|---:|---|
| Version switch, argument plumbing, header-text tables | 3 | 5 | One binary with a `--gatk-version` profile |
| HC traversal, filters, downsampler, isActive on the old AF path | 3 | 6 | Filters and region defaults are unchanged |
| Read-threading graph, dangling ends, SeqGraph cleanup | 5 | 8 | |
| KBestHaplotypeFinder (recursive) | 5 | 10 | Tie order is undocumented |
| Identity-hash non-determinism: measure and document | 2 | 5 | Plus an **unknown** residual. Exact emulation is not costed and likely infeasible |
| Region trimming and padding (`trimLegacy` check) | 3 | 6 | |
| CIGAR, left-align, realign, hard-clip | 5 | 10 | |
| Overlapping-pair quality (`FragmentUtils`) | 2 | 4 | |
| Asymmetric likelihood cap | 0.5 | 1 | |
| Old exact AF calculator and old `GenotypingEngine` QUAL/MLEAC | 5 | 8 | |
| GL ulp path, PL/AD read sets, events, phasing, spanning deletions | 6 | 12 | |
| ReferenceConfidenceModel | 2 | 4 | |
| GVCF writer and `HomRefBlock` | 1 | 2 | |
| Rank-sum annotations (old `ReadUtils` helpers), ExcessHet, NON_REF GT/PL rules | 3 | 7 | |
| BQSR / ApplyBQSR / GATKReport (posterior, known-site mask) | 4 | 8 | |
| Picard 2.18.26/2.20.3 MarkDuplicates, queryname mode | 4 | 8 | Plus 1-2 d for metrics if needed |
| htsjdk 2.18.1 BAM/VCF details | 1 | 2 | |
| JDK 8 transcendental check and fixes | 2 | 10 | **Unknown** until measured |
| Scale validation and residual triage (tier0 → tier1) | 10 | 20 | |
| **Total (a)** | **66** | **136** | GenotypeGVCFs 4.0.12.0 is optional: **+5-10** |

### (b) Additional cost of GKL-mode exactness (x86 reference)

| Component | Low | High | Notes |
|---|---:|---:|---|
| PairHMM float kernel emulation (software FTZ, 0.8.6 `Context` tables, ΣM/ΣX order) | 8 | 15 | Exact by construction if no FMA **[unverified for AVX-512]** |
| Double fallback (1e-28f) and `log10f`/`log10` libm matching | 2 | 5 | **Unknown:** host glibc |
| x86 oracle harness (JNI dump of GKL 0.8.6 likelihoods), AVX vs AVX-512 check | 3 | 5 | Needs x86 hardware |
| GPU/NEON path of the float kernel matching the scalar emulation | 3 | 8 | Only if we accelerate it |
| GKL Smith-Waterman | 0 | 0 | Not used by 4.0.12.0 (JAVA default) |
| **Subtotal (b), required** | **16** | **33** | |
| ISA-L level-2 deflate, byte-exact BAM/gVCF.gz (optional) | 15 | 30 | Recommend a decompressed-content gate instead |

**Bottom line.** A Java-mode 4.0.12.0 port is roughly **3-6 engineer-months**. Most of that is HaplotypeCaller code that 4.6.2.0 rewrote, and our existing code carries over less than the first estimate assumed. It has a structural ceiling: 4.0.12.0's own identity-hash non-determinism. GKL-mode PairHMM adds **3-7 weeks**, and is probably required to match UKB's x86 gVCFs.
