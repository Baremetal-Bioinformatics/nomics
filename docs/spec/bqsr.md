# BaseRecalibrator / ApplyBQSR byte-exactness specification (GATK 4.6.2.0 / htsjdk 4.2.0)

Target: reproduce byte-for-byte (a) the recalibration table written by `gatk BaseRecalibrator` and (b) the records of the BAM
written by `gatk ApplyBQSR` (modulo BGZF framing), as run by nf-core/sarek 3.10.0 with exactly the command lines in the task
context (one interval group = the full Agilent v7 BED, three `--known-sites` VCFs, all other arguments default, no
`GatherBQSRReports`).

Every claim is tagged **[source]** (repo/file/method/lines at the pinned tag) or **[inferred]**. Line numbers refer to the
files as checked out at the tag (not to comment-stripped views).

Pinned sources:

* `gatk` = github.com/broadinstitute/gatk tag `4.6.2.0`, paths under `src/main/java/org/broadinstitute/hellbender/` unless noted.
* `htsjdk` = github.com/samtools/htsjdk tag `4.2.0`, paths under `src/main/java/htsjdk/`.
* `cm3` = Apache commons-math3 **3.5** (`gatk/build.gradle` L279-281 pins `strictly '3.5'`), paths under `src/main/java/org/apache/commons/math3/`.
* `barclay` = github.com/broadinstitute/barclay tag `5.0.0`.
* Version verification **[source]**: `gatk/build.gradle` L62-64 (`htsjdk 4.2.0`, `picard 3.4.0`, `barclay 5.0.0`), L541-543
  (manifest `Implementation-Version` = archive version `4.6.2.0`, `Toolkit-Short-Name` = `GATK`).
* The `gatk` launcher script passes `-Dsamjdk.compression_level=2` and `-Dsamjdk.use_async_io_write_samtools=true`
  (`gatk/gatk` L43-46) - BGZF level 2 only; out of scope. **[source]**

htsjdk record re-encoding facts are reused from `markdup.md` §5.3-5.4 (tag key order, smallest-int retyping, `bin` rule) and are
not re-derived here.

---

## 0. Execution path and the big-picture findings

* Both tools are `ReadWalker`s (`tools/walkers/bqsr/BaseRecalibrator.java` L85, `tools/walkers/bqsr/ApplyBQSR.java` L76). The
  traversal is `engine/ReadWalker.java` `traverse()` L92-108: `getTransformedReadStream(makeReadFilter())` then for each read
  `apply(read, new ReferenceContext(reference, readInterval), new FeatureContext(features, readInterval))` with
  `readInterval = getReadInterval(read)` (L115-117: `null` for unmapped reads, else `SimpleInterval(contig, start, end)`). **[source]**
* `GATKTool.getTransformedReadStream` (`engine/GATKTool.java` L382-394): `reads` iterator → `map(preTransformer)` →
  `filter(countingFilter)` → `map(postTransformer)`. Pre-transformer is identity for both tools (L355-357); the post-transformer is
  identity for BaseRecalibrator and `BQSRReadTransformer` for ApplyBQSR (`ApplyBQSR.java` L106-109). **[source]**
* Findings that drive everything below (each proven in its section):
  1. **BAQ is off** (`utils/recalibration/RecalibrationArgumentCollection.java` L102-104 `enableBAQ = false`, hidden argument). With BAQ
     off every per-base error weight is exactly `0.0` or `1.0`, so all `Errors` are integers and summation order is irrelevant (§4).
  2. **Only the `M` (base substitution) event type is accumulated** (`computeIndelBQSRTables = false`, L106-108; engine L101). No `I`/`D` rows.
  3. **Empirical qualities are integers** (`utils/recalibration/RecalDatum.java` L38 `private int empiricalQuality`, Feb-2025 change);
     the only non-integer number in the whole report is `EstimatedQReported` of `RecalTable0` (§5.2).
  4. **ApplyBQSR does not quantize by default**: `--quantize-quals` defaults to `0` which calls `QuantizationInfo.noQuantization()`
     (identity map for Q0..Q92) (§8.3). The `Quantized` table is written but effectively unused.
  5. ApplyBQSR **recomputes** empirical qualities from the table's counts with hierarchical priors; it never reads the
     `EmpiricalQuality` column (§8.4).
  6. BaseRecalibrator uses **supplementary** reads (no `NotSupplementary` filter) but not duplicates/secondary/QC-fail/MAPQ 0/MAPQ 255;
     ApplyBQSR emits and recalibrates **every** well-formed read overlapping the intervals, including duplicates, supplementary,
     secondary, QC-fail, MAPQ 0 and mate-placed unmapped reads (§1.3, §8.1).

---

## 1. Traversal: which reads are visited, in which order

### 1.1 Intervals

* `--intervals chr1_12146-12310.bed`: `IntervalUtils.parseIntervalArguments` (`utils/IntervalUtils.java` L306-320) detects a Feature
  file and calls `featureFileToIntervals` (L354-364), which decodes the BED with htsjdk `BEDCodec` (`tribble/bed/BEDCodec.java`
  `decode` L127-199, L140: `start = parseInt(col2) + 1`, `end = parseInt(col3)`, i.e. 0-based half-open → 1-based closed). **[source]**
* `IntervalArgumentCollection` defaults (`cmdline/argumentcollections/IntervalArgumentCollection.java`): `intervalSetRule = UNION`
  L75, `intervalPadding = 0` L82, `intervalExclusionPadding = 0` L89, `intervalMergingRule = ALL` L97. `loadIntervals`
  (`IntervalUtils.java` L229-249) → `sortAndMergeIntervals`; `mergeIntervalLocations` (L964-985) merges overlapping **and abutting**
  intervals. `parseIntervals` (`IntervalArgumentCollection.java` L186-241) sets `traverseUnmapped` only if the literal interval
  `unmapped` was given (L234-237) - it was not, so **unplaced unmapped reads are never visited**. **[source]**
* Merging/sorting never changes the set of reads selected (selection is "overlaps any interval"), so a Rust port may use the raw BED.

### 1.2 Query semantics (htsjdk)

* `ReadWalker.setReadTraversalBounds` (L59-66) → `ReadsPathDataSource.setTraversalBounds` (`engine/ReadsPathDataSource.java` L296-304);
  iteration (L398-433) wraps the single reader in a `SamReaderQueryingIterator` (`utils/iterators/SamReaderQueryingIterator.java`):
  intervals are converted to `QueryInterval`s and `QueryInterval.optimizeIntervals` (htsjdk `samtools/QueryInterval.java` L74-96:
  sort, merge overlapping/abutting) (L90-108), then `reader.queryOverlapping(queryIntervals)` (L128-131). **[source]**
* `BAMFileReader.createIndexIterator` (htsjdk `samtools/BAMFileReader.java` L942-955) iterates the index chunks of all intervals in
  file order and filters with `BAMQueryMultipleIntervalsIteratorFilter` (`samtools/BAMQueryMultipleIntervalsIteratorFilter.java`
  `compareToFilter` L21-38, `compareIntervalToRecord` L40-58). A record is emitted iff it overlaps an interval (`contained=false`),
  with: `alignmentEnd = POS` for a read with flag 0x4 **and** POS != 0 (mate-placed unmapped) and `alignmentEnd = getAlignmentEnd()`
  (POS + reference length of CIGAR − 1) otherwise; overlap = `interval.start <= alignmentEnd && POS <= interval.end`, same reference index.
  Because the record stream is coordinate sorted and the interval cursor only advances when `interval.end < POS`, each record is
  emitted **at most once** even if it spans several intervals. **[source]**
* Consequences: (a) a mate-placed unmapped read is visited iff its POS (1-based) lies inside a BED interval; (b) the 71 such reads in
  the golden are those; (c) unplaced unmapped reads (refID −1) are never visited; (d) record order = input (coordinate) order, so
  **ApplyBQSR output order = md.bam order restricted to the visited set**. **[source + inferred(c,d)]**

### 1.3 Read filters

`GATKTool.makeReadFilter` (L280-286) → `GATKReadFilterPluginDescriptor.getResolvedInstances`
(`cmdline/GATKPlugin/GATKReadFilterPluginDescriptor.java` L198-219): the tool's default filters in declared order, then user filters
(none). All are pure predicates AND-ed together, so order only affects log counts, never output. **[source]**

*BaseRecalibrator* (`BaseRecalibrator.java` L153-181), in order:
`MAPPING_QUALITY_NOT_ZERO` (MAPQ != 0, `engine/filters/ReadFilterLibrary.java` L107-110), `MAPPING_QUALITY_AVAILABLE` (MAPQ != 255, L100-103),
`MAPPED` (`!isUnmapped()`, L93-96), `NOT_SECONDARY_ALIGNMENT` (0x100 clear, L166-169), `NOT_DUPLICATE` (0x400 clear, L152-155),
`PASSES_VENDOR_QUALITY_CHECK` (0x200 clear, L187-190), then `WellformedReadFilter`. **Supplementary (0x800) reads are kept.** **[source]**

*ApplyBQSR*: only `WellformedReadFilter` (`ReadWalker.getDefaultReadFilters` L130-132; `ApplyBQSR` does not override). **[source]**

`WellformedReadFilter.createFilter` (`engine/filters/WellformedReadFilter.java` L62-73) = `VALID_ALIGNMENT_START` (unmapped or POS > 0,
L250-253) ∧ `VALID_ALIGNMENT_END` (unmapped or end−start+1 ≥ 0, L273-276) ∧ `AlignmentAgreesWithHeader` (`utils/read/ReadUtils.java`
L894-903: mapped reads must have a contig in the dictionary and POS ≤ contig length) ∧ `HAS_READ_GROUP` (`RG` tag present, L76-79) ∧
`HAS_MATCHING_BASES_AND_QUALS` (SEQ length == QUAL length, L69-72) ∧ `READLENGTH_EQUALS_CIGARLENGTH` (unmapped or SEQ length == CIGAR
read length, L224-228) ∧ `SEQ_IS_STORED` (SEQ length > 0, L239-242) ∧ `CIGAR_CONTAINS_NO_N_OPERATOR` (L31-34). **[source]**

`isUnmapped()` for GATK reads is `flag 0x4 || RNAME == "*" || POS == 0` (`utils/read/SAMRecordToGATKReadAdapter.java` L425-429). **[source]**
For bwa/MarkDuplicates output every visited read passes Wellformed **[inferred]**; the golden count 873,049 is "records overlapping
the BED under §1.2 that pass Wellformed".

### 1.4 Whole-read processing

A visited read is processed in its entirety (all bases), not only the bases inside the interval(s). **[source: no interval clipping anywhere in §2/§8]**

---

## 2. BaseRecalibrator per-read pipeline (`utils/recalibration/BaseRecalibrationEngine.java`)

`BaseRecalibrator.apply` (L188-190) calls `recalibrationEngine.processRead(read, referenceDataSource, featureContext.getValues(knownSites))`.
`processRead` (L115-149):

1. `read = makeReadTransform().apply(originalRead)` (§2.1). If `read.isEmpty()` (no bases left) the read is skipped entirely and not counted (L119-121).
2. `RecalUtils.updatePlatformForRead` (`utils/recalibration/RecalUtils.java` L490-511): `FORCE_PLATFORM`/`DEFAULT_PLATFORM` are null; a read group with no `PL` throws. `PL:ILLUMINA` is present, nothing changes. **[source]**
3. `calculateIsSNPOrIndel` (§2.3) → `nErrors`.
4. `baqArray = flatBAQArray(read)` because `enableBAQ == false` (L134; `flatBAQArray` L507-511 fills `'@'` = 64). **[source]**
5. `covariates = RecalUtils.computeCovariates(read, header, covariates, true, keyCache)` (§3).
6. `skip = calculateSkipArray(read, knownSites)` (§2.2).
7. `snpErrors = calculateFractionalErrorArray(isSNP, baqArray)` (§2.4) (insertion/deletion arrays computed but unused).
8. `updateRecalTablesForRead(info)` (§4).

### 2.1 Read transformers (L290-299), applied in this order

1. **`consolidateCigar`** (L301-306): `read.setCigar(new CigarBuilder().addAll(cigar).make())` on the **original** read.
   `utils/read/CigarBuilder.java` `add` L53-115 / `make` L128-144 (default `removeDeletionsAtEnds = true`): drops zero-length elements,
   merges consecutive equal operators, drops a deletion that directly follows leading clips (or a leading insertion), drops a trailing
   deletion (also one just before trailing clips / before a trailing insertion), and reorders `I D` → `D I`. The alignment start is
   **not** shifted when a leading deletion is removed. bwa mem does not emit such CIGARs **[inferred]**; for well-formed bwa CIGARs this
   is the identity. **[source]**
2. `setDefaultBaseQualities`: no-op because `defaultBaseQualities = -1` (L315-321). **[source]**
3. `resetOriginalBaseQualities`: no-op because `useOriginalBaseQualities = false` (L308-313). `OQ` is never consulted. **[source]**
4. **`ReadClipper::hardClipAdaptorSequence`** (`utils/clipping/ReadClipper.java` L367-378) **[source]**:
   * `adaptorBoundary = read.getAdaptorBoundary()` (`ReadUtils.getAdaptorBoundary` L527-536):
     `CANNOT_COMPUTE` (Integer.MIN_VALUE) unless `hasWellDefinedFragmentSize` (L547-580): TLEN != 0, paired (0x1), neither read nor
     mate unmapped (GATK `mateIsUnmapped` = 0x8 || RNEXT "*" || PNEXT 0, adapter L456-460), strands differ, and (reverse: `end > mateStart`;
     forward: `start <= mateStart + TLEN`). Then boundary = reverse ? `mateStart − 1` : `start + |TLEN|`.
   * No clip if boundary cannot be computed or `!isInsideRead` (`ReadUtils.isInsideRead` L743-745: `start <= boundary <= end`).
   * Reverse strand: `hardClipByReferenceCoordinatesLeftTail(boundary)` = `clipByReferenceCoordinates(-1, boundary, HARDCLIP)`;
     forward: `hardClipByReferenceCoordinatesRightTail(boundary)` = `clipByReferenceCoordinates(boundary, -1, HARDCLIP)` (L167-187).
   * `clipByReferenceCoordinates` (L419-472): for the left tail, `stop = idx(refStop) − (op consumes read bases ? 0 : 1)` and `start = 0`;
     for the right tail, `start = idx(refStart)`, `stop = len − 1`; where `idx = ReadUtils.getReadIndexForReferenceCoordinate(read, coord)`
     (L713-715 → L683-706) walks the CIGAR from `read.getSoftStart()` treating **`S` as consuming reference**, returns the read index of the
     base at that coordinate, or for a coordinate inside a deletion the index of the first base **after** it (operator `D`); returns
     `READ_INDEX_NOT_FOUND` (−1) if before soft start or past the end → then no clipping (L452-454). Soft clips still exist at this
     stage, so soft-clipped bases are treated as aligned for adaptor clipping.
   * Clipping (`utils/clipping/ClippingOp.java` `applyHardClipBases` L182-224): new bases/quals = the unclipped range; CIGAR =
     `CigarUtils.clipCigar(cigar, start, stop+1, H)` (`utils/read/CigarUtils.java` L342-384: converts the clipped read bases to `H`,
     drops deletions that end up at the clip edge, keeps existing hard clips); if clipping the left side, POS += `alignmentStartShift`
     (L389-411: number of reference bases consumed by the first `stop+1` read bases, counting deletions lying wholly inside the clipped
     span, not counting `S`). If the whole read is clipped, `ReadUtils.emptyRead` (L784-800: unmapped, no bases) is returned → read skipped.
   * Example: forward read POS 1000, CIGAR `5S146M`, TLEN 140 → boundary 1140, `softStart` 995, `start = 1140−995 = 145` → bases 145..150
     removed → `5S140M6H`, POS unchanged. Reverse read POS 1000, `146M5S`, mate at 1003 → boundary 1002 → `stop = 2` → `3H143M5S`, POS 1003.
5. **`ReadClipper::hardClipSoftClippedBases`** (L281-324): converts leading/trailing `S` to `H`, removing those bases from
   SEQ/QUAL (right op first, then left); POS unchanged (`alignmentStartShift` over `S` is 0). After this step the read has no `S`. **[source]**

Each clipping step returns a **copy** (`ClippingOp.apply` L43-84 copies the SAMRecord); only step 1 mutates the original, which is
irrelevant because BaseRecalibrator writes no reads. Everything below operates on the clipped read: `read.getLength()` is the clipped
length, `read.getStart()/getEnd()` the clipped span. **[source]**

### 2.2 Skip array (L339-347) and known sites (L353-379)

`skip[i] = !isRegularBase(base_i) || qual_i < PRESERVE_QSCORES_LESS_THAN (6) || knownSites[i]`. `isRegularBase` (`utils/BaseUtils.java`
L173-175, map L44-56): `A C G T` (either case) and `*` are regular; `N` and IUPAC codes are not. **[source]**

Known-site features: `featureContext.getValues(knownSites)` (`engine/FeatureContext.java` L256-266 → L169-174) concatenates, for each
`--known-sites` input in command-line order, `FeatureManager.getFeatures` (`engine/FeatureManager.java` L357-364) =
`FeatureDataSource.queryAndPrefetch(readInterval)` (`engine/FeatureDataSource.java` L577-596, refill L608-627 with lookahead 1000,
`ReadWalker` L42) which, through `FeatureCache` (`engine/FeatureCache.java` `trimToNewStartPosition` L173-204,
`getCachedFeaturesUpToStopPosition` L213-224), returns exactly the records with `feature.start <= read.end && feature.end >= read.start`,
where the query interval is the **original (unclipped-cigar, pre-transform) read's** `[POS, alignmentEnd]`. **[source]**

VCF record extent (htsjdk `variant/vcf/AbstractVCFCodec.java` L432-440): `start = POS`; `end = INFO/END` if present, else
`POS + len(REF) − 1`. So a SNP covers 1 base; an MNP `len(REF)` bases; a deletion `REF=ACGT ALT=A` covers POS..POS+3 (**the padding base
is included**); an insertion `REF=A` covers only POS; multi-allelic records use `len(REF)`; symbolic records use `END`. **FILTER is ignored**
(filtered records count). ID/QUAL/genotypes are irrelevant. **[source]**

`calculateKnownSites` maps each feature onto the clipped read (L356-377) **[source]**:

```
softStart = read.getSoftStart(); softEnd = read.getSoftEnd()         // == start/end: no S left
skip feature if feature.end < softStart || feature.start > softEnd
(s, opS) = idx(feature.start); featureStart = (s == -1) ? 0 : s; if opS == D: featureStart -= 1
(e, opE) = idx(feature.end);   featureEnd   = (e == -1) ? readLength : e
if featureStart > readLength: featureStart = featureEnd = readLength
fill knownSites[max(0, featureStart) .. min(readLength, featureEnd + 1)) = true
```

with `idx` as in §2.1 (for a coordinate inside a deletion: index of the base after it). Hence a known site falling in a read deletion
masks the base before and the base after the deletion; a site starting before the read masks from base 0; a site ending after the read
masks to the last base; insertions in the read are never separately masked (they inherit the mask through the index range).

### 2.3 Reference comparison (`calculateIsSNPOrIndel` L392-445)

* `refBases = ref.queryAndPrefetch(contig, read.getStart(), read.getEnd()).getBases()` from a second `ReferenceDataSource`
  (`BaseRecalibrator.java` L150) = `ReferenceFileSource` (`engine/ReferenceFileSource.java` L77-79) over
  `CachingIndexedFastaSequenceFile` with `preserveCase=false, preserveIUPAC=false` (`utils/fasta/CachingIndexedFastaSequenceFile.java`
  L324-325, L343-344): reference bases are **upper-cased and IUPAC codes become `N`**. The `MD` tag is never used. **[source]**
* Walk the CIGAR: `M/=/X`: `snp[readPos] = basesAreEqual(readBase, refBase) ? 0 : 1` (`BaseUtils.basesAreEqual` L118-120 compares
  `simpleBaseToBaseIndex` values, L145-149: case-insensitive, `A/a/*→0, C/c→1, G/g→2, T/t→3`, everything else −1; so `N` vs `N` is
  "equal", `N` vs `A` is a mismatch). `D`: marks `isDel` at `readPos−1` (forward) or `readPos` (reverse). `N`: skip ref. `I`: marks
  `isIns` at `readPos−1` (forward) or `readPos+len` (reverse). `S`: skip read bases (none left). `H/P`: nothing. **[source]**
* Only `snp[]` matters for the output (`M` event only); `nErrors` only decides whether to run BAQ, which is off.

### 2.4 Fractional errors (`calculateFractionalErrorArray` L454-486)

With a flat BAQ array every entry equals `NO_BAQ_UNCERTAINTY`, so `fractionalErrors[i] = (double) errorArray[i]` ∈ {0.0, 1.0}
(L466-468); the block-averaging branch (L470, L488-499) is never entered. **[source]** Therefore every `isError` passed to the tables is
exactly 0 or 1 and all `Errors` values are integers (§4.2).

---

## 3. Covariates (`utils/recalibration/covariates/`)

`StandardCovariateList` (`StandardCovariateList.java` L37-50): exactly `[ReadGroupCovariate, QualityScoreCovariate, ContextCovariate,
CycleCovariate]`, indices 0..3; `covariateNames()` = `"ReadGroupCovariate,QualityScoreCovariate,ContextCovariate,CycleCovariate"`
(L110-112); report names are the class name before `Covariate` (`Covariate.java` L56-58): `ReadGroup`, `QualityScore`, `Context`, `Cycle`.
`computeCovariates` (`RecalUtils.java` L528-532) allocates a `PerReadCovariateMatrix` `[3 events][readLength][4]` (reused per read
length through `CovariateKeyCache`, LRU of 500, `CovariateKeyCache.java` L16-18) and each covariate writes its key for every offset
(`populatePerReadCovariateMatrix` L134-140). Every covariate writes all `readLength` offsets for the `M` event, so the cache never leaks
stale keys. **[source]**

### 3.1 ReadGroupCovariate (`ReadGroupCovariate.java`)

* Identifier = `PU` of the read's `@RG` if present, else `ID` (`getReadGroupIdentifier` L70-73). Keys are assigned in header `@RG` order
  starting at 0 (L32-47, `getReadGroupIDs` L109-111). Same key for every base (L50-61). `formatKey` returns the identifier (L76-79). **[source]**
* **The `ReadGroup` column of the tables is the PU string, not the ID** (sarek sets `PU`) **[inferred - verify on golden]**.

### 3.2 QualityScoreCovariate (`QualityScoreCovariate.java` L20-35)

Key for the `M` event = `read.getBaseQuality(i)` (the signed byte, 0..93 in practice). `formatKey` = `%d` (L38-40); `maximumKeyValue` = 93
(`QualityUtils.MAX_SAM_QUAL_SCORE` = htsjdk `SAMUtils.MAX_PHRED_SCORE` = 93, `utils/QualityUtils.java` L16, htsjdk `samtools/SAMUtils.java`
L119). In BaseRecalibrator the quality is that of the clipped read (so `BI/BD` insertion/deletion qualities and the defaults of 45 are only
used for the never-emitted `I`/`D` events). **[source]**

### 3.3 ContextCovariate (`ContextCovariate.java`)

Settings: `mismatchesContextSize = 2`, `indelsContextSize = 3`, `lowQualTail = 2` (L41-61). `M`-event keys only (L100-103). **[source]**

1. `strandedClippedBases = getStrandedClippedBytes(read, 2)` (L126-137): `ReadClipper.clipLowQualEnds(read, 2, WRITE_NS)`
   (`ReadClipper.java` L229-258 + `ClippingOp.applyWriteNs` L113-117): on a **copy**, every base from the 3' end backwards while
   `qual <= 2` and every base from the 5' end forwards while `qual <= 2` is overwritten with `N` (`Nucleotide.N.encodeAsByte()` = `'N'`,
   `utils/Nucleotide.java` L241-243). If **all** bases are `<= 2` the read becomes empty (L247-249) - then every context key is 0
   (L81-87; irrelevant because such bases are never recalibrated, §2.2/§8.3). Then, for reverse-strand reads, the bases are
   reverse-complemented (`BaseUtils.simpleReverseComplement` L247-255 / `simpleComplement` L222-245: `ACGT`↔`TGCA`, anything else
   unchanged, so `N` stays `N`). Length is unchanged by `WRITE_NS`.
2. `getReadContextAtEachPosition(bases, 2, mask)` (L176-233) yields one key per position `i` of the stranded sequence: the key of the
   2-mer `bases[i−1..i]` (preceding base first, the base itself last), or `−1` if `i < 1` or if any base of the 2-mer is not `ACGT`
   (`simpleBaseToBaseIndex == −1`); after a non-ACGT base at position `p`, positions `p .. p+contextSize−1` are `−1` (penalty logic
   L212-230); the initial window handles the same rule (L195-210).
3. Key encoding (`keyFromContext` L247-260, `createMask` L152-160, shift L219-221): `key = contextSize | Σ_j baseIndex(b_j) << (4 + 2j)`
   with `j = 0` the **oldest** base. For size 2: `key = 2 | older<<4 | current<<6`, e.g. `"AC"` = `2 | 0 | 1<<6 = 66`, `"TT"` = `2 | 48 | 192 = 242`.
   `formatKey` (L140-145, `contextFromKey` L268-286) prints oldest base first; `maximumKeyValue` (L289-299) for `max(2,3)=3` = 1011 (table dim 1012).
4. Position mapping (L100-103, `getStrandedOffset` L114-116): the key computed at stranded index `i` is stored at read offset
   `i` (forward) or `readLength − i − 1` (reverse). So for a reverse-strand read, the context of read offset `o` is the reverse complement
   of read bases `o+1, o` (i.e. the two bases in sequencing order ending at `o`).

### 3.4 CycleCovariate (`CycleCovariate.java`)

`cycleKey(i, read, indel=false, 500)` (L65-92) with `readLength = read.getLength()` of the read **as seen by the covariate** (clipped
read in BaseRecalibrator, full read in ApplyBQSR) **[source]**:

```
f = (paired && secondOfPair) ? -1 : 1
cycle = reverse ? readLength*f + i*(-f) : f + i*f      // forward: f*(i+1); reverse: f*(readLength - i)
key = |cycle| << 1 | (cycle < 0 ? 1 : 0)                // keyFromCycle L108-120; |cycle| > 500 throws
```

`formatKey` prints the signed cycle (`cycleFromKey` L97-103, `%d`). `maximumKeyValue` = 1001 (L52-54). Hard clipping therefore changes
the cycle numbering of a reverse-strand read's bases (quirk, replicate as is). With 151 bp reads the possible values are ±1..±151.

---

## 4. Accumulation (`updateRecalTablesForRead` L255-288, `RecalUtils.increment*` L595-634)

### 4.1 Tables (`utils/recalibration/RecalibrationTables.java` L43-72, `utils/collections/NestedIntegerArray.java`)

* `readGroupTable[numRG][3]`, `qualityScoreTable[numRG][94][3]`, `contextTable[numRG][94][1012][3]`, `cycleTable[numRG][94][1002][3]`
  (`allTables` order: RG, QS, Context, Cycle, L55/L60/L65-71). Event index = `EventType.ordinal()`: `M=0, I=1, D=2` (`EventType.java` L3-6). **[source]**
* For each offset with `!skip[offset]`, for the single cached event `M` (L266-267): `qual = baseQuals[offset]`
  (`ReadRecalibrationInfo.getQual` L72-80), `isError = snpErrors[offset]`, `rg = cov[0]`, `q = cov[1]`:
  `incrementDatum3keys(qualityScoreTable, qual, isError, rg, q, 0)`; then for `i = 2, 3`: if `cov[i] >= 0`:
  `incrementDatum4keys(table_i, qual, isError, rg, q, cov[i], 0)`. A missing leaf is created as `new RecalDatum(1, isError, qual)`
  (`createDatumObject` L643-645), an existing one gets `increment(1L, isError)`. Context key −1 is skipped (never stored); cycle keys are
  always ≥ 0 for `M`. **[source]**
* Single-threaded; no partial tables are merged. Traversal order is the BAM order of §1.2. **[source: ReadWalker is sequential]**

### 4.2 `RecalDatum` arithmetic (`RecalDatum.java`)

* Fields: `reportedQuality` (double), `empiricalQuality` (int, −1 = unset), `numObservations` (long), `numMismatches` (double, stored
  **×100000.0**, `MULTIPLIER` L24). Constructor L75-84: `numMismatches = numMismatches*MULTIPLIER`, `reportedQuality = (double) byte`.
  `increment(long, double)` L238-242: `numObservations += n; numMismatches += inc*MULTIPLIER`. `getNumMismatches()` = `numMismatches/MULTIPLIER` (L218-220). **[source]**
* With `inc ∈ {0.0, 1.0}` all stored values are integer multiples of 100000.0 (exact in binary64 up to 2^53 ≈ 9·10^15, i.e. ≤ 9·10^10
  errors), and `k·100000/100000` is exactly `k`. Hence `Errors` is always an exact integer and **accumulation order does not matter**. **[source + inferred]**

---

## 5. Finalisation (`BaseRecalibrator.onTraversalSuccess` L193-205 → `finalizeData` L159-164)

### 5.1 Read-group table (`collapseQualityScoreTableToReadGroupTable` L179-202)

Iterate the quality-score table's leaves in `getAllLeaves` order (`NestedIntegerArray.java` L202-221: nested ascending key order =
by rg, then q ascending, then event). For the first leaf of an (rg, event): `rgDatum = new RecalDatum(qualDatum)` (copy L90-95, so
`reportedQuality = q_1` exactly, `empiricalQuality` copied = −1). For each subsequent leaf, `rgDatum.combine(leaf)` (L106-117):

```
expected = N_acc * 10^(r_acc / -10.0) + N_k * 10^(q_k / -10.0)   // calcExpectedErrors L260-262 → QualityUtils.qualToErrorProb(double) L138-141: Math.pow(10.0, qual / -10.0)
N_acc += N_k; E_acc += E_k                                        // increment
r_acc = -10 * Math.log10(expected / N_acc)                        // Java Math.log10 (HotSpot intrinsic)
empiricalQuality = -1
```

Note the expected errors use the **running** `r_acc` (not the sum over original qualities), evaluated before the counts are added. **[source]**

### 5.2 Rounding (`roundTableValues` L208-216, `MathUtils.roundToNDecimalPlaces` `utils/MathUtils.java` L927-931)

For every leaf of every table: `numMismatches = round(getNumMismatches(), 2)`, `reportedQuality = round(reportedQuality, 4)` with
`round(x, n) = Math.round((x + Math.ulp(x)) * 10^n) / 10^n` (`Math.pow(10, n)` = exactly 100.0 / 10000.0; `Math.round(double)` = floor(x+0.5)
as a long, exact). Integers and integer qualities are unchanged; only the read-group row's `EstimatedQReported` is actually rounded.
The rounded value is the nearest double to `k/10000`, so printing it with 4 decimals is exact (§7.4). **[source]**

### 5.3 Empirical quality (`RecalDatum.getEmpiricalQuality` L169-181 → `calcEmpiricalQuality` L267-275)

Computed lazily and cached per datum (first call wins). Report-time calls (`RecalUtils.generateReportTables` L256) use
`getEmpiricalQuality()` = prior `reportedQuality` (L169-171); the quantizer (§6) calls it first for the quality table, with the same prior,
so caching does not change results. **[source]**

```
k' = (long)(getNumMismatches() + 0.5) + 1          // SMOOTHING_CONSTANT = 1 (L60); k integer here → k' = k + 1
n' = numObservations + 2
Q  = bayesianEstimateOfEmpiricalQuality(n', k', prior)   (L303-309)
empiricalQuality = min(Q, 93)                            // MAX_RECALIBRATED_Q_SCORE = 93 (L21)
```

`bayesianEstimateOfEmpiricalQuality(n, k, prior)`: `logPost[q] = getLogPrior(q, prior) + getLogBinomialLikelihood(q, n, k)` for
`q = 0..60` (`MAX_REASONABLE_Q_SCORE = 60`, `QualityUtils.java` L31; 61 bins), result = `MathUtils.maxElementIndex(logPost)`
(`MathUtils.java` L786-797: **first** index of the maximum, strict `>`; ties → smallest q). **[source]**

* `getLogPrior(q, prior)` (L330-333): `d = min(|(int)(q − prior)|, 40)` - **Java `(int)` truncates toward zero**, so e.g. prior 36.6 gives
  `d = 0` for both q = 36 and q = 37; `logPriorCache[d] = NormalDistribution(mean 0, sd 0.5).logDensity(d)` (L315-327) =
  `-0.5*(d/0.5)^2 − (FastMath.log(0.5) + 0.5*FastMath.log(2*FastMath.PI))` (cm3 `distribution/NormalDistribution.java` L147, L176-179) =
  `−2d² − 0.22579135264472744` **[source; constant inferred, equals the correctly rounded value]**. The constant is common to all q and
  cannot change the argmax except through last-bit rounding of the sums.
* `getLogBinomialLikelihood(q, n, k)` (L350-367): `n` ≤ 2^31−2 here (otherwise scaled); `p = QualityUtils.qualToErrorProb((double) q)` =
  `Math.pow(10.0, q / -10.0)`; `ll = MathUtils.logBinomialProbability((int) n, (int) k, p)` (`MathUtils.java` L610-612) =
  `new BinomialDistribution(null, n, p).logProbability(k)` (cm3 `distribution/BinomialDistribution.java` L114-127) =
  `SaddlePointExpansion.logBinomialProbability(k, n, p, 1.0 − p)` (cm3 `distribution/SaddlePointExpansion.java` L177-199); if the result
  is infinite or NaN it becomes `−Double.MAX_VALUE`.

`SaddlePointExpansion` (cm3 3.5, all logs are `FastMath.log`) **[source]**:

```
logBinomialProbability(x, n, p, q):
  if x == 0:      p < 0.1 ? -getDeviancePart(n, n*q) - n*p : n*log(q)
  elif x == n:    q < 0.1 ? -getDeviancePart(n, n*p) - n*q : n*log(p)
  else:           ret = getStirlingError(n) - getStirlingError(x) - getStirlingError(n-x) - getDeviancePart(x, n*p) - getDeviancePart(n-x, n*q)
                  f = (2π * x * (n-x)) / n;  ret = -0.5*log(f) + ret
getStirlingError(z) (L106-127): z < 15 → EXACT_STIRLING_ERRORS[(int)(2z)] (table L51-81; 2z is always integral here) ;
                  else ((0.083333333333333333333 - (0.00277777777777777777778 - (0.00079365079365079365079365 - (0.000595238095238095238095238 - 0.0008417508417508417508417508/z2)/z2)/z2)/z2)/z) with z2 = z*z
getDeviancePart(x, mu) (L144-165): if |x-mu| < 0.1*(x+mu): d = x-mu; v = d/(x+mu); s1 = v*d; s = NaN; ej = 2*x*v; v *= v; j = 1;
                  while (s1 != s) { s = s1; ej *= v; s1 = s + ej/(2j+1); ++j }  return s1
                  else x*log(x/mu) + mu - x
```

Here `x = k' ≥ 1` and `x = n'` is impossible, so only the general branch is used; `n' < 15` uses the table. The 31-entry table must be
copied verbatim (values are decimal literals, parsed to the nearest double).

### 5.4 `EstimatedQReported` and the RG row's empirical quality

`RecalTable0.EstimatedQReported` = the rounded `r_acc` of §5.1-5.2; `RecalTable0.EmpiricalQuality` uses that rounded value as prior
with the summed counts. `RecalTable1`/`RecalTable2` rows use prior = their integer `QualityScore`. **[source]**

---

## 6. Quantization (`utils/recalibration/QuantizationInfo.java`, `QualQuantizer.java`)

`QuantizationInfo(tables, 16)` (L32-48): `hist[0..93] = 0`; for every leaf of the **quality-score table** (`getAllValues`, any order):
`hist[fastRound(datum.getEmpiricalQuality())] += numObservations` (`MathUtils.fastRound` L428-430; the value is already an integer).
So `Count[q]` = total observations of `RecalTable1` rows whose `EmpiricalQuality == q`. Then `new QualQuantizer(hist, 16, 6)`
(`minInterestingQual = MIN_USABLE_Q_SCORE = 6`, L51-54). **[source]**

`QualQuantizer` (L82-97, `quantize` L260-280) **[source]**:

1. Leaves `I_q = QualInterval(q, q, nObs = hist[q], nErrors = (int) floor(hist[q] * qualToErrorProb((byte) q)), level 0, fixedQual = q)`
   for q = 0..93, in a `TreeSet` ordered by `qStart` (`compareTo` L177-179).
2. While `size > 16`: `mergeLowestPenaltyIntervals` (L286-313): walk adjacent pairs `(left, right)` in ascending order; `merged = left.merge(right)`
   (L192-207: `nObs/nErrors` summed, `fixedQual = −1`, `subIntervals = [left, right]` in that order); keep the **first** pair whose
   `merged.getPenalty()` is strictly smaller than the current minimum; replace the pair by `merged`.
3. `getErrorRate()` (L150-157): fixed-qual leaf → `qualToErrorProb((byte) fixedQual)` (cache `QualityUtils.java` L48-57 of
   `Math.pow(10.0, q / -10.0)`); merged with `nObs == 0` → `0.0`; else `(nErrors + 1) / (1.0 * (nObs + 1))`.
4. `getPenalty()` = `calcPenalty(getErrorRate())` (L209-245): `0.0` if the rate is 0; for a leaf: `0` if `qEnd <= 6`, else
   `|Math.log10(leafRate) − Math.log10(rate)| * nObs`; for a merged node: sum over `subIntervals` in insertion order (left subtree first,
   recursively).
5. `getQual()` (L162-167): fixed leaf → `fixedQual`; merged → `QualityUtils.errorProbToQual(rate)` (L208-228) =
   `boundQual((int) Math.round(−10 * Math.log10(rate)), 93)` with `boundQual = max(min(x, 93), 1)` (L367-369). **For an all-empty merged
   interval `rate = 0` → `log10(0) = −∞` → `Math.round(+∞) = Long.MAX_VALUE` → `(int)` of that double saturates to `Integer.MAX_VALUE` → `93`.**
   This is why the golden shows `Q0-10 … QuantizedScore 93`.
6. `intervalsToMap` (L325-338): `QuantizedScore[q] = interval(q).getQual()`.

Behaviour on this data (3 non-empty bins e1 < e2 < e3, everything else 0) **[inferred from the algorithm]**: every merge of two
all-empty intervals (or of empties with bins ≤ 6) has penalty 0 and is chosen leftmost-first; merging an empty interval with a non-empty
bin has a small positive penalty (the leaf's own rate differs from the smoothed merged rate) and is never chosen while a zero-penalty
pair exists. 78 merges are needed (94 → 16) and 87 zero-penalty merges are available, so the final 16 intervals are
`[0..e1−1]`, `[e1]`, `[e1+1..e2−1]`, `[e2]`, `[e2+1..e3−1]`, `[e3]`, then the empties above `e3` merged left-to-right until exactly 16
intervals remain (one long merged interval followed by singletons up to 93). Quantized values: merged empties → 93; non-empty singletons →
their own q; high singletons → their own q.

---

## 7. GATKReport v1.1 text format (`utils/report/`)

### 7.1 Tables and their construction

`RecalUtils.outputRecalibrationReport` (`RecalUtils.java` L300-303) → `createRecalibrationGATKReport` (L339-345) adds, to a
`GATKReport` whose tables live in a `TreeMap<String, GATKReportTable>` (`GATKReport.java` L33, `addTable` L148-150), the tables
`Arguments`, `Quantized`, `RecalTable0`, `RecalTable1`, `RecalTable2`. Output order is therefore **alphabetical by table name**. **[source]**

* **Arguments** (`RecalibrationArgumentCollection.generateReportTable` L148-188): `new GATKReportTable("Arguments",
  "Recalibration argument collection values used in this run", 2, SORT_BY_COLUMN)`; columns `Argument` (format `%s`) and `Value`
  (format `""` → data type `Unknown`, printed format `%s`, `GATKReportColumn.java` L28-39, L106-108). Rows are added with
  `addRowID(name, true)` (row id and first column = name) and `set(name, "Value", value)`; values and their `%s` renderings:
  `covariate` → `ReadGroupCovariate,QualityScoreCovariate,ContextCovariate,CycleCovariate`; `no_standard_covs` → `false`;
  `run_without_dbsnp` → `false`; `solid_recal_mode` → `SET_Q_ZERO`; `solid_nocall_strategy` → `THROW_EXCEPTION`;
  `mismatches_context_size` → `2`; `indels_context_size` → `3`; `mismatches_default_quality` → `-1`; `deletions_default_quality` → `45`;
  `insertions_default_quality` → `45`; `maximum_cycle_value` → `500`; `low_quality_tail` → `2`; `default_platform` → `null`
  (a Java `null` is stored as the string `"null"`, `GATKReportTable.set` L360-361); `force_platform` → `null`; `quantizing_levels` → `16`;
  `recalibration_report` → `null`; `binary_tag_name` → `null`. **[source]**
* **Quantized** (`QuantizationInfo.generateReportTable` L70-83): description `Quality quantization map`, columns `QualityScore %d`,
  `Count %d`, `QuantizedScore %d`; 94 rows q = 0..93 with `(Integer q, Long count, Byte quantized)`. **[source]**
* **RecalTable0/1/2** (`RecalUtils.generateReportTables` L190-271): description `""`; columns, in order:
  * RecalTable0: `ReadGroup %s`, `EventType %s`, `EmpiricalQuality %.4f`, `EstimatedQReported %.4f`, `Observations %d`, `Errors %.2f`;
  * RecalTable1: `ReadGroup %s`, `QualityScore %d`, `EventType %s`, `EmpiricalQuality %.4f`, `Observations %d`, `Errors %.2f`;
  * RecalTable2: `ReadGroup %s`, `QualityScore %d`, `CovariateValue %s`, `CovariateName %s`, `EventType %s`, `EmpiricalQuality %.4f`,
    `Observations %d`, `Errors %.2f`; the Context and Cycle tables are **both** written into this one table (L227-235), context rows first
    (table order) - irrelevant after sorting.
  * Cell values: `ReadGroup` = identifier string; `QualityScore` = `formatKey` string `"37"` which `fixType` (`GATKReportTable.java`
    L396-424) converts to `Long` because the column type is Integer; `CovariateValue` = `formatKey` string (`"AC"`, `"-12"`);
    `CovariateName` = `Context`/`Cycle`; `EventType` = `"M"` (`EventType.toString` L55-57); `EmpiricalQuality` = `Double` (int widened);
    `EstimatedQReported` = `Double`; `Observations` = `Long`; `Errors` = `Double` (`getNumMismatches()`). **[source]**

### 7.2 Row order (`GATKReportTable.write` L485-544, `ROW_COMPARATOR` L24-54)

All five tables use `SORT_BY_COLUMN`: rows are sorted (stable `Collections.sort`) lexicographically over the columns left to right,
comparing cells of equal class natively (`Long`, `Double`, …) and anything else by `toString().compareTo`. Hence:
`Arguments` by argument name (String order, as listed in the task); `Quantized` by q (0..93); `RecalTable1` by (RG, q numeric);
`RecalTable2` by (RG, q numeric, **CovariateValue as a String**, CovariateName). String order puts `-` (0x2D) before digits before
letters, so within one quality the rows come as cycles `-1, -10, -100, -101, …, -109, -11, -110, …, -151, -16, …, -19, -2, -20, …, -99`,
then `1, 10, 100, …, 151, 16, …, 99`, then contexts `AA, AC, AG, AT, CA, …, TT`. **[source]**

### 7.3 Layout (`GATKReport.print` L175-180, `GATKReportTable.write` L495-543, `writeRow` L546-574, `GATKReportColumn`, `GATKReportColumnFormat`)

```
#:GATKReport.v1.1:5                                        println: prefix + version + ":" + number of tables
#:GATKTable:<ncol>:<nrow>:<fmt1>:<fmt2>:...:;              printf; one ":fmt" per column (getFormat(): "%s" for Unknown)
#:GATKTable:<name>:<description>                           printf "...%s:%s\n"  (empty description → line ends with ':')
<col1 name padded>  <col2 name padded>  ...                 names with "%-Ws", joined by two spaces
<row>...                                                    cells with "%Ws" (right) or "%-Ws" (left), joined by two spaces
<empty line>                                                println() after the last row of every table
```

* Column width `W` = max(length of column name, length of every formatted value set in that column) (`updateMaxWidth` L120-122;
  `formatValue` L80-96 formats with the column format, Doubles/Floats in an Unknown column with `%.8f`). **[source]**
* Alignment (L18, L124-127, `isRightAlign` L69-71): starts RIGHT and becomes LEFT (permanently) at the first value whose formatted string
  is not `null`/`NA`/`Infinity`/`-Infinity`/`NaN` and not `NumberUtils.isCreatable` (commons-lang3 numeric literal). Consequences:
  `Arguments.Argument` and `Arguments.Value` (first value is the covariate list) are LEFT; `ReadGroup` is LEFT (unless the identifier is
  numeric); `EventType` (`M`) LEFT; `CovariateName` LEFT; `CovariateValue` becomes LEFT as soon as any context (`AA`) is set - all RecalTable2
  rows are set before writing, so it is LEFT for the whole table (cycle values included); all numeric columns RIGHT. **[source]**
* Padding is applied to the **last** column too, so LEFT-aligned last columns and every header line whose last column name is shorter
  than `W` carry trailing spaces (e.g. `Errors` header when values like `123456.00` are wider than 6). The header line is always
  `%-Ws` for every column. **[source]**
* Separator between columns: exactly two spaces (`out.printf("  ")`). Line terminator: `println` (platform `line.separator`, `\n` on Linux)
  for the report header, column-name line, rows and the blank line; literal `\n` for the name/description line. The file ends with
  the blank line after the last `RecalTable2` row, i.e. `...<last row>\n\n`. **[source]**

Concrete widths for this run **[inferred from the rules]**: `Arguments`: `Argument` W=26 (`mismatches_default_quality`), `Value` W=72
(the covariate list) → every data line is 100 characters, left-aligned, space-padded; header `Argument` + 18 spaces + 2 + `Value` + 67
spaces. `Quantized`: `QualityScore` 12, `Count` max(5, digits of the largest count), `QuantizedScore` 14, all RIGHT. Header line
`#:GATKTable:3:94:%d:%d:%d:;` then `#:GATKTable:Quantized:Quality quantization map`. RecalTable0 header `#:GATKTable:6:1:%s:%s:%.4f:%.4f:%d:%.2f:;`
then `#:GATKTable:RecalTable0:`; RecalTable1 `#:GATKTable:6:3:%s:%d:%s:%.4f:%d:%.2f:;`; RecalTable2 `#:GATKTable:8:954:%s:%d:%s:%s:%s:%.4f:%d:%.2f:;`.

### 7.4 Number formatting (`writeRow` L558-570, `formatValue`)

* `%d` on `Long`/`Integer`/`Byte` → plain decimal, no grouping. `%s` on `Boolean` → `true`/`false`; on `Byte -1` → `-1`.
* `%.4f`/`%.2f` on a finite `Double` → `String.format(fmt, obj)`. Java's `Formatter` rounds **HALF_UP on the shortest-repr decimal digits**
  of the double (`FormattedFloatingDecimal.applyPrecision`), not on the exact binary value, so Java and C/Rust can differ for ties
  (`1.005 → "1.01"` in Java, `"1.00"` in C) **[inferred, JDK behaviour]**. This cannot bite here: `EmpiricalQuality` is an integer
  (`37.0000`), `Errors` is an integer (`12345.00`), and `EstimatedQReported` was pre-rounded to 4 places (§5.2) so the shortest
  representation has ≤ 4 decimals and any correct formatter prints the same digits. Rust `format!("{:.4}", x)` of the pre-rounded value is
  therefore byte-identical. Non-finite doubles would print via `toString()` (`Infinity`, `NaN`) - impossible here.
* Java `(int)`/`(long)` casts and `Math.round` semantics above must be replicated exactly where noted (§5.3, §6.5).

---

## 8. ApplyBQSR

### 8.1 Reads emitted

Same traversal as §1.1-1.2, only the `WellformedReadFilter` (§1.3). Every visited read is passed through `BQSRReadTransformer.apply`
and then `outputWriter.addRead(read)` (`ApplyBQSR.java` L117-120); nothing is dropped or reordered. Duplicates, supplementary, secondary,
QC-fail, MAPQ-0/255 and mate-placed unmapped reads are all **emitted and recalibrated** (the transformer never inspects flags or
alignment; the covariates of §3 only need bases, qualities, strand/pair flags). Unplaced unmapped reads are never visited. **[source]**

### 8.2 Table parsing (`utils/recalibration/RecalibrationReport.java` L55-74)

The table is re-read as text (`GATKReportTable(BufferedReader)` L106-185: fixed-width split on the column-name word starts; `%d` → `Long`,
`%.Nf` → `Double.parseDouble`, `%s` → String, `GATKReportDataType.Parse` L124-143). Read groups = the set of `ReadGroup` strings of
`RecalTable0` (`GATKReport.getReadGroups` L77-84), keyed in sorted order (`initializeReadGroupCovariates` L190-194). For every row a
`RecalDatum(nObservations, nErrors, (byte) 1)` is created and `setReportedQuality(EstimatedQReported)` (RecalTable0) or
`setReportedQuality(QualityScore)` (others) (`getRecalDatum` L283-296). **`EmpiricalQuality` is not read** (commented out, L294).
Because every printed number round-trips exactly (§7.4), the parsed datums equal the in-memory ones of BaseRecalibrator. The `Quantized`
table yields `quantizedQuals[q]` and counts (`initializeQuantizationTable` L304-317). `RecalibrationArgumentCollection` values are
re-read from `Arguments` (L325-397). **[source]**

### 8.3 Per-read transform (`transformers/BQSRReadTransformer.java` L141-234)

Arguments (`tools/ApplyBQSRArgumentCollection.java` L21-31, `ApplyBQSRUniqueArgumentCollection.java` L30-74): `quantizationLevels = 0`,
`staticQuantizationQuals = []`, `roundDown = false`, `emitOriginalQuals = false`, `globalQScorePrior = −1.0`, `allowMissingReadGroups = false`,
`PRESERVE_QSCORES_LESS_THAN = 6`, `useOriginalBaseQualities = false`. **[source]**

Constructor (L79-109): `quantizationLevels == 0` → `quantizationInfo.noQuantization()` (`QuantizationInfo.java` L56-60): `quantizedQuals[i] = i`
for `i = 0..92`; **index 93 keeps the value parsed from the `Quantized` table**. No static mapping. **[source]**

`apply(read)` **[source]**:

1. No OQ handling (`useOriginalBaseQualities`, `emitOriginalQuals` false).
2. `computeCovariates(read, header, covariates, false, keyCache)` on the **unclipped, untransformed** read (soft-clipped bases included;
   adaptor bases included; cycle computed with the full length).
3. `read.clearAttribute("BI")`, `read.clearAttribute("BD")` (L159-160) - absent in the input; `SAMRecord.setAttribute(tag, null)` on an
   absent tag is a no-op on the attribute list but still marks the BAM record's binary block stale (htsjdk `samtools/BAMRecord.java`
   L214-219), so every record is re-encoded from decoded fields (`markdup.md` §5.3-5.4 apply).
4. `rgKey` from the read's `@RG` identifier (PU/ID); unknown → exception. `readGroupDatum = readGroupTable[rgKey][M]`.
5. `recalibratedQuals = copy of QUAL`. For each offset: if `QUAL[offset] < 6` → unchanged (L201-203). Else, with `q = QUAL[offset]`
   (= the quality-score covariate), `ctx = cov[2]`, `cyc = cov[3]`:
   `qualityScoreDatum = qualityScoreTable[rg][q][M]` (may be null if no row); `ctxDatum = contextTable[rg][q][ctx][M]` if `ctx >= 0`
   (null if absent); `cycDatum = cycleTable[rg][q][cyc][M]` (null if absent).
   `prior = readGroupDatum.getReportedQuality()` (= parsed `EstimatedQReported`; `globalQScorePrior <= 0`).
   `raw = hierarchicalBayesianQualityEstimate(prior, readGroupDatum, qualityScoreDatum, ctxDatum, cycDatum)` (L259-277):
   ```
   Qrg   = readGroupDatum.getEmpiricalQuality(prior)                       // §5.3 with n = Observations(RG) + 2, k = Errors(RG) + 1
   Qq    = qualityScoreDatum == null ? Qrg : qualityScoreDatum.getEmpiricalQuality(Qrg)    // prior = Qrg (an integer), NOT q
   delta = Σ over non-null {ctxDatum, cycDatum}: datum.getEmpiricalQuality(Qq) − Qq        // prior = Qq
   raw   = Qq + delta
   ```
   `newQual = quantizedQuals[boundQual(fastRound(raw), 93)]` with `boundQual = max(min(x, 93), 1)` (L237-239, `QualityUtils.boundQual` L367-369).
   All terms are integers, so `raw` is an exact integer-valued double, `fastRound` is exact, and with the identity map
   **`newQual = clamp(Qq + (Qctx − Qq) + (Qcyc − Qq), 1, 93)`** (missing covariate rows / `ctx = −1` simply drop their term; a quality
   with no `RecalTable1` row gives `Qq = Qrg` and no covariate terms). Index 93 (reached only if `raw ≥ 92.5`) maps through the parsed
   `Quantized` table value for Q93.
6. `read.setBaseQualities(recalibratedQuals)` (adapter L313-323 → htsjdk `BAMRecord.setBaseQualities` L208-211, marks stale).

Caching note (L271 TODO): `getEmpiricalQuality(prior)` caches the first result per datum. The RG datum always sees the same prior; the
quality datum for `q` always sees `Qrg` (one read group); each covariate datum `(rg, q, key)` always sees the same `Qq`. So the cache
never changes a result here. The empirical qualities used by ApplyBQSR for `RecalTable1`/`RecalTable2` datums are computed with
**different priors than the ones printed in the table** and must be recomputed by a port. **[source]**

### 8.4 Output record and header

* Record changes relative to md.bam: **only the QUAL bytes** of bases with original `QUAL >= 6` that resolve to a different value.
  No tag is added or removed (`OQ` not written, `BI`/`BD` absent, `PG:Z:MarkDuplicates` kept), flags/CIGAR/POS unchanged. Re-encoding
  (§8.3 step 3) rewrites tags in ascending binary-key order with smallest integer types and recomputes `bin`; because md.bam was itself
  written by htsjdk 4.2.0 with the same rules (`markdup.md` §5.3-5.4), this is the identity on the tag block and `bin`. **[source + inferred]**
* Header: `GATKTool.createSAMWriter` (L855-870) → `getHeaderForSAMWriter` (L962-972): the input header object plus one program record
  (`addOutputSAMProgramRecord = true`, L111) appended **after** the existing `@PG` lines (htsjdk `samtools/SAMFileHeader.java`
  `addProgramRecord` L200-207): `ID` = `createProgramGroupID` (L992-1002) = `getToolName()` = `"GATK ApplyBQSR"` (L1041-1043:
  `Toolkit-Short-Name` + simple class name; collisions would append `.1`, `.2`, …); attributes set in the order `VN`
  (`getVersion()` = manifest `Implementation-Version` = `4.6.2.0`, `cmdline/CommandLineProgram.java` L390-392 → `utils/runtime/RuntimeUtils.java`
  L90-93), `CL` (`getCommandLine()` = Barclay `CommandLineArgumentParser.getCommandLine` L220-246: `ApplyBQSR` + explicitly set arguments
  in parser definition order + non-null defaults; **copy the literal from the golden**), `PN` (`GATK ApplyBQSR`). No `PP`. The resulting
  line is `@PG\tID:GATK ApplyBQSR\tVN:4.6.2.0\tCL:ApplyBQSR …\tPN:GATK ApplyBQSR` **[source; CL literal inferred]**.
* `ReadUtils.createCommonSAMWriterFromFactory` (L1001-1017) writes `header.clone()` (text round trip, attribute order preserved,
  `@HD` re-emitted as `VN:1.6 GO:none SO:coordinate`, `markdup.md` §1.6) through `SAMFileWriterFactory.makeWriter(..., presorted=true)`;
  `createOutputBamIndex = true` (L84) writes a `.bai` (out of scope); `createOutputBamMD5 = false`. **[source]**

---

## 9. Other factors affecting exact output

* `--reference`: used only for §2.3 (and CRAM); the FASTA is upper-cased/IUPAC→N on read. `--tmp-dir`: no effect. **[source]**
* `--known-sites` order: only the union of masks matters. **[source]**
* `low_quality_tail 2`: only affects the context covariate (N-overwrite at both read ends, §3.3); it does **not** clip or skip bases. **[source]**
* `PRESERVE_QSCORES_LESS_THAN 6`: skips bases in BaseRecalibrator (§2.2) and leaves them untouched in ApplyBQSR (§8.3). A Q2 base stays Q2. **[source]**
* `maximum_cycle_value 500`: throws if `|cycle| > 500`; 151 bp reads never reach it. **[source]**
* Read-group `PL`: must be present (else `MalformedRead` exception); value otherwise unused (`Utils.warnOnNonIlluminaReadGroups` L1097 only logs). **[source]**
* Run-to-run nondeterminism: none in the table (confirmed byte-identical twice); in the BAM only the `CL:` string of the new `@PG` can vary
  with argument spelling/order. **[source + inferred]**
* JDK: the GATK 4.6.2.0 image runs Java 17 **[inferred]**; `Math.log10`/`Math.pow` are HotSpot intrinsics (≤ 1 ulp from the exact
  result, not necessarily correctly rounded); `FastMath.log` (cm3) is a table-driven implementation also ≤ 1 ulp. See §11 for where this can matter.

---

## 10. Rust implementation notes (algorithm summary)

1. Load BED → 1-based closed intervals (merge optional). Stream md.bam in order; select records by the §1.2 rule; apply the Wellformed
   predicates; for BaseRecalibrator additionally MAPQ ∉ {0, 255}, mapped, not 0x100/0x400/0x200.
2. BaseRecalibrator per read: consolidate CIGAR (§2.1.1), adaptor hard-clip (§2.1.4), soft→hard clip (§2.1.5); skip empty reads;
   `snp[]` vs upper-cased reference; `skip[]` from non-ACGT, `qual < 6`, known-site mask (VCF POS..END overlap, §2.2); covariates on the
   clipped read (§3); increment `(rg, q)` and `(rg, q, ctx)` if `ctx >= 0` and `(rg, q, cyc)` with `obs += 1`, `err += snp`.
3. Finalise: RG row by iterative combine over q ascending (§5.1) with `pow(10, x/-10.0)` and `log10`; round to 4 decimals (§5.2);
   empirical qualities by the 61-bin argmax with the truncating prior and the SaddlePoint likelihood (§5.3); quantizer (§6).
4. Write the five tables with the exact layout of §7 (sorting, widths, alignment, trailing spaces, two blank-line conventions).
5. ApplyBQSR: parse nothing - reuse the in-memory tables (identical by construction) - or parse the text; recompute `Qrg`, `Qq`, `Qctx`,
   `Qcyc` with hierarchical priors; `newQual = clamp(Qq + Σ(Qcov − Qq), 1, 93)` for bases with `qual >= 6`; covariates on the full read;
   write every visited record unchanged except QUAL; header = input header + the `GATK ApplyBQSR` `@PG` line.

---

## 11. Checks against the tier0 golden

Files: `HG002_WES.recal.table`, `HG002_WES.recal.bam`, `HG002_WES.md.bam`, the BED and the three VCFs under
`/home/<user>/omics/runs/tier0-wes-20261008-r1/results/`.

Table structure:

1. Line 1 is `#:GATKReport.v1.1:5`; tables appear in the order Arguments, Quantized, RecalTable0, RecalTable1, RecalTable2; the
   description lines are exactly `#:GATKTable:Arguments:Recalibration argument collection values used in this run`,
   `#:GATKTable:Quantized:Quality quantization map`, `#:GATKTable:RecalTable0:`, `#:GATKTable:RecalTable1:`, `#:GATKTable:RecalTable2:`
   (trailing colon, nothing after). Every table is followed by exactly one empty line; the file ends with `\n\n`.
2. Arguments: 17 data lines each exactly 100 bytes (`%-26s` + 2 spaces + `%-72s`), header line = `Argument` + 18 spaces + 2 spaces +
   `Value` + 67 spaces (also 100 bytes). Row order and values as in §7.1.
3. Quantized header line is `QualityScore  Count` + (W_count − 5 + 2) spaces + `QuantizedScore` where `W_count` = number of digits of the
   largest `Count` (7 if the largest is 1320847); rows are right-aligned to widths 12/W_count/14.
4. `Σ_q Count[q]` = `RecalTable0.Observations` = `Σ RecalTable1.Observations`; for each q, `Count[q]` = Σ of `Observations` over
   RecalTable1 rows with `EmpiricalQuality == q.0000` (so at most 3 bins are non-zero).
5. QuantizedScore pattern (§6): with non-empty bins e1 < e2 < e3: q in `[0, e1−1]` → 93; `e1` → e1; `[e1+1, e2−1]` → 93; `e2` → e2;
   `[e2+1, e3−1]` → 93; `e3` → e3; above e3 one merged block → 93 followed by singletons printing their own q, with exactly 16 distinct
   intervals in total. Golden: `Q0-10 → 93`, `Q11 → 11` matches with e1 = 11.
6. RecalTable0: one row; `EventType` `M`; `Observations`/`Errors` equal the column sums of RecalTable1; `Errors` ends in `.00`;
   `EstimatedQReported` recomputed by §5.1-5.2 from the RecalTable1 `(QualityScore, Observations)` pairs in ascending q matches to 4 decimals;
   `EmpiricalQuality` = `argmax` of §5.3 with that prior, printed as an integer with `.0000`.
7. `ReadGroup` cell = the `PU` value of the md.bam `@RG` line (not the `ID`) if `PU` is present (§3.1).
8. RecalTable1: 3 rows = the distinct base qualities ≥ 6 present in the visited reads (NovaSeq bins, e.g. 12/23/37); each
   `EmpiricalQuality` = `bayesian(Observations + 2, Errors + 1, q)`; `Errors` integer-valued.
9. RecalTable2 has 954 = 3 × (302 + 16) rows: for every q all 302 cycle values `-151..-1, 1..151` and all 16 dinucleotides; row order per §7.2
   (`-1, -10, -100, …`); `CovariateValue` left-aligned (cycle numbers padded on the right); `CovariateName` is `Context`/`Cycle`.
   For each q: `Σ Observations over Cycle rows = RecalTable1.Observations(q)` and the same for `Errors`; `Σ over Context rows ≤` the
   RecalTable1 value (contexts `-1` at the first sequenced base and after `N`/low-quality-tail bases are excluded).
10. Trailing whitespace: the RecalTable1/2 header lines end with spaces iff the widest `Errors` value is wider than 6 characters;
    RecalTable0 header likewise; data lines of RecalTable0/1/2 never end with spaces (last column numeric, right-aligned).

Semantics checks (computable from md.bam + BED + VCFs):

11. `samtools view -c` of recal.bam = number of md.bam records with refID ≥ 0 that overlap a BED interval under §1.2 (mate-placed
    unmapped: POS inside an interval) = 873,049; the 71 unmapped records all have POS inside a BED interval; no record with refID −1.
12. recal.bam record order = md.bam order filtered; every record is identical to its md.bam counterpart except the QUAL field
    (compare `samtools view` without the QUAL column); `PG:Z:MarkDuplicates` retained; no `OQ`/`BI`/`BD`.
13. Header = md.bam header + exactly one appended `@PG` with `ID:GATK ApplyBQSR`, `VN:4.6.2.0`, `PN:GATK ApplyBQSR`, no `PP`,
    tag order `ID VN CL PN`; `@HD` unchanged.
14. For every base with original QUAL < 6 (e.g. Q2) the output QUAL is unchanged; for every other base
    `QUAL' = clamp(Qq + (Qctx − Qq) + (Qcyc − Qq), 1, 93)` with `Qq/Qctx/Qcyc` recomputed per §8.3 (priors `Qrg`, `Qq`) from the table counts,
    using the context of the **full** read (soft clips included, low-quality-tail bases N-masked) and the cycle from the full read length.
    Duplicates, supplementary and mate-placed unmapped reads follow the same rule (they are **not** emitted unchanged).
15. Bases that are `N` in the read keep the cycle term and drop the context term; the first sequenced base of every read (offset 0 forward,
    offset len−1 reverse) and the base right after a `Q≤2` tail or an `N` drop the context term.
16. RecalTable1 `Observations(q)` = number of (read, base) pairs over BaseRecalibrator-visited reads (MAPQ ∉ {0,255}, mapped, not
    0x100/0x200/0x400, Wellformed) after adaptor+soft-clip hard clipping, with `qual == q`, base ∈ ACGT, and not under a known-site mask
    (VCF `POS..POS+len(REF)−1` or `INFO/END`, FILTER ignored, deletion-adjacent bases masked); `Errors(q)` = those that mismatch the
    upper-cased reference. A spot check on a few reads with soft clips / TLEN shorter than the read is the cheapest validation of §2.1.
17. Cycle rows exist for ±151 at each q only if unclipped 151 bp forward reads (cycle 151) and unclipped reverse reads (cycle 151 at
    offset 0) occur - consistent with 302 cycle rows per q.

---

## 12. Riskiest exactness traps (summary)

* Floating point: `EstimatedQReported` chain (`pow`, `log10`, iterative combine order, `+ulp` then round-half-up to 4 places);
  the 61-bin posterior argmax (SaddlePoint expansion with `FastMath.log`, Stirling table, deviance series, truncating `(int)` prior
  distance, first-max tie-break); quantizer `floor(nObs * 10^(-q/10))` (exact integer products such as `10 * 0.1` flip with a 1-ulp `pow`
  error) and `log10` penalties with leftmost-first zero-penalty ties and the `93` saturation quirk.
* Integer/semantic: supplementary reads counted; cycle computed on the clipped length; adaptor clipping treats `S` as aligned; known-site
  masks include the VCF padding base and deletion-adjacent read bases; `PU` as read-group key; ApplyBQSR priors differ from the printed
  `EmpiricalQuality`; ApplyBQSR does not quantize; mate-placed unmapped reads are recalibrated.
* Text: alphabetical table order, string-sorted cycle values, left/right alignment inference, padded last columns and header lines with
  trailing spaces, `\n\n` at EOF, `null` strings in Arguments.
