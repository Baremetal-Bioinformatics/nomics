# MarkDuplicates byte-exactness specification (GATK 4.6.2.0 / Picard 3.4.0 / htsjdk 4.2.0)

Target: reproduce, record-for-record and byte-for-byte (modulo BGZF framing, which is out of scope), the
output BAM and the metrics file produced by `gatk MarkDuplicates` as run by nf-core/sarek 3.10.0
(process `GATK4_MARKDUPLICATES`) with the exact command line given in the task context.

Every claim is tagged **[source]** (with repo/file/method/lines at the pinned tag) or **[inferred]**.

Pinned sources:

* `picard` = github.com/broadinstitute/picard tag `3.4.0`, paths under `src/main/java/picard/`.
* `htsjdk` = github.com/samtools/htsjdk tag `4.2.0`, paths under `src/main/java/htsjdk/`.
* `gatk` = github.com/broadinstitute/gatk tag `4.6.2.0`.
* `barclay` = github.com/broadinstitute/barclay tag `5.0.0` (GATK's command-line parser).

Version verification **[source]**: `gatk/build.gradle` L62 `htsjdkVersion = ... '4.2.0'`, L63
`picardVersion = ... '3.4.0'`, L64 `barclayVersion = ... '5.0.0'`; `picard/build.gradle` L62 pins htsjdk `4.2.0` too.

---

## 0. Execution path and global settings

* GATK's `Main` finds the Picard class and wraps it in `PicardCommandLineProgramExecutor`
  (`gatk/.../hellbender/Main.java` L366-368); `runCommandLineProgram` strips the tool name and calls
  `program.instanceMain(mainArgs)` (L160-167); the executor forwards to `picard.cmdline.CommandLineProgram.instanceMain`
  (`gatk/.../cmdline/PicardCommandLineProgramExecutor.java` L36-42). GATK sets `picard.useLegacyParser=false`
  (`Main.java` L49), so the Barclay parser is used. **[source]**
* `picard/cmdline/CommandLineProgram.instanceMain` (L184-300) **[source]**:
  * builds the two metrics-file headers: `new StringHeader(commandLine)` and `new StringHeader("Started on: " + new Date())` (L213-216);
  * `BlockCompressedOutputStream.setDefaultCompressionLevel(COMPRESSION_LEVEL=2)` (L221-222);
  * `SAMFileWriterImpl.setDefaultMaxRecordsInRam(500000)` (L226-228) - irrelevant here because output is written presorted;
  * `SAMFileWriterFactory.setDefaultCreateIndexWhileWriting(true)` because `CREATE_INDEX=true` (L231-233);
  * Intel deflater/inflater factories (L240-246) - compression only;
  * `SamReaderFactory.setDefaultValidationStringency(LENIENT)` (L250). LENIENT only logs validation errors; it never mutates records. **[source]**
* `MarkDuplicates()` constructor forces `DUPLICATE_SCORING_STRATEGY = SUM_OF_BASE_QUALITIES` (`MarkDuplicates.java` L246-248). **[source]**
* `OpticalDuplicateFinder` is created in `customCommandLineValidation` with
  `(READ_NAME_REGEX = DEFAULT_READ_NAME_REGEX, OPTICAL_DUPLICATE_PIXEL_DISTANCE = 100, MAX_OPTICAL_DUPLICATE_SET_SIZE = 300000, LOG)`
  (`picard/sam/markduplicates/util/AbstractOpticalDuplicateFinderCommandLineProgram.java` L56, L73-80). **[source]**
* `doWork` (`MarkDuplicates.java` L256-451): `FLOW_MODE=false` so `calcHelper = this`; `useBarcodes = false`
  (no BARCODE_TAG / READ_ONE_BARCODE_TAG / READ_TWO_BARCODE_TAG) (L261-266). Steps: `buildSortedReadEndLists`
  (pass 1), `generateDuplicateIndexes`, re-open inputs, write output (pass 2), `finalizeAndWriteMetrics`. **[source]**

---

## 1. Input merging

### 1.1 Opening

`AbstractMarkDuplicatesCommandLineProgram.openInputs(eagerlyDecode)` (`picard/sam/markduplicates/util/AbstractMarkDuplicatesCommandLineProgram.java` L237-268) **[source]**:

* One `SamReader` per `--INPUT`, in `--INPUT` order (L241-250). Pass 1 uses `EAGERLY_DECODE`; pass 2 does not (L243-245). This has no effect on output content (see 5.4).
* `ASSUME_SORT_ORDER == null` and `ASSUME_SORTED == false`, so no header sort order is overridden (L251-259).
* With more than one input: `new SamFileHeaderMerger(headers.get(0).getSortOrder() /* coordinate */, headers, mergeDictionaries=false)`
  and `new MergingSamRecordIterator(headerMerger, readers, assumeSorted=false)` (L264-266). The merged header is `headerMerger.getMergedHeader()`.
  With exactly one input the input header and reader iterator are used directly (L261-262).

`openInputs` is called twice (pass 1 and pass 2); both calls build an identical merged header for the same `--INPUT` order. **[source]** (L283, L509)

### 1.2 Header merge (`htsjdk/samtools/SamFileHeaderMerger.java`)

Constructor L158-194 **[source]**:

1. `mergedHeader = new SAMFileHeader()` -> attributes `{VN:1.6}` (`SAMFileHeader.java` L107-109, `CURRENT_VERSION="1.6"` L44).
2. Sequence dictionary: `getSequenceDictionary(headers)` (L522-535) takes the **first** header's `SAMSequenceDictionary` object
   and asserts every other one equal via `SequenceUtil.assertSequenceDictionariesEqual` (htsjdk `SequenceUtil.java` L327); differing
   dictionaries throw because `mergeDictionaries=false` (L166-172). `hasMergedSequenceDictionary` stays `false`.
   Consequence: `@SQ` lines (names, lengths, and any extra attributes in their original order) are those of the first `--INPUT` file.
3. `@PG`: `mergeProgramGroups` (L248-335), see 1.3. Added to the merged header in the returned (ID-sorted) order (L178-180).
4. `@RG`: `mergeReadGroups` (L213-238), see 1.4. `setReadGroups(result)` (L184).
5. `setGroupOrder(GroupOrder.none)` then `setSortOrder(coordinate)` (L185-187). `SAMFileHeader.setGroupOrder` / `setSortOrder`
   store `GO`/`SO` via `AbstractSAMHeaderRecord.setAttribute` (`SAMFileHeader.java` L257-260, L279-282) into a `LinkedHashMap`
   (`AbstractSAMHeaderRecord.java` L40, L75-81). Insertion order is therefore `VN, GO, SO`, which is exactly the golden
   `@HD VN:1.6 GO:none SO:coordinate`. **[source]**
6. `@CO`: comments of all inputs, in input order (L189-193). bwa/samtools write none. **[inferred: none expected]**

### 1.3 `@PG` merge and record-ID renaming

`mergeProgramGroups` (L248-335) **[source]**:

* Collects every `@PG` of every input in `(input order, line order)`; records with no `PP` are processed first ("roots"),
  then records whose `PP` points at a just-processed record **from the same input file** (L279-319), level by level.
* Each level is merged by `mergeHeaderRecords` (L415-490): records are binned by ID in a `LinkedHashMap`, then within an ID by
  **record equality** (`SAMProgramRecord.equals` = ID and *all attributes* equal; `SAMProgramRecord.java` L111-121).
  Identical records collapse into one. For each distinct record of a given ID, in first-seen order:
  * if the ID is not yet taken: keep the ID; `++recordCounter` (L458-463);
  * else: `newId = id + "." + positiveFourDigitBase36Str(recordCounter++)` until free (L464-474). `positiveFourDigitBase36Str`
    (L498-512) is plain base-36 with digits `0-9A-Z` (no zero padding; `0 -> "0"`).
  * `recordCounter` is reset to 0 once before the whole `@PG` merge (L268) and keeps counting across levels.
  * A per-input translation map `old ID -> new ID` is recorded (L476-483), and `PP` values of later levels are translated through it (L350-393).
* Final list is sorted by ID with `Collections.sort` (stable) using `String.compareTo` (L116-121, L332).

Worked example for sarek's 12 inputs, each carrying `@PG ID:bwa ... CL:bwa mem ... <chunk-specific fastq names>` and
`@PG ID:samtools PN:samtools PP:bwa ... CL:samtools sort ...<chunk-specific>` (all CLs distinct) **[inferred from sarek; mechanics are source]**:

* bwa level: first input keeps `bwa` (counter -> 1); the remaining 11 get `bwa.1, bwa.2, ..., bwa.9, bwa.A, bwa.B`
  (counter values 1..11), in `--INPUT` order. Counter is now 12.
* samtools level: first input keeps `samtools` (counter -> 13); the rest get `samtools.D, samtools.E, ..., samtools.N`
  (counter values 13..23). Each `PP` is rewritten to that input's bwa ID.
* Output order after sort: `bwa, bwa.1, ..., bwa.9, bwa.A, bwa.B, samtools, samtools.D, ..., samtools.N`.
* Which input's `CL` ends up under which suffix depends on `--INPUT` order. This is the only header nondeterminism and matches the
  observation that the two reference runs differed only in `@PG` lines.
* `hasProgramGroupCollisions()` becomes true, so `MergingSamRecordIterator.next` would rewrite a record's `PG:Z` tag through the
  translation table (`MergingSamRecordIterator.java` L144-150) **if the record has a PG tag**. bwa mem / samtools sort do not add
  `PG` tags to records **[inferred; consistent with the observation that records were identical across runs - if input records
  carried `PG:Z:bwa`, output tags would be `PG:Z:MarkDuplicates`/`MarkDuplicates.<n>` chosen per input and vary with input order]**.

If two inputs had byte-identical `@PG` lines they collapse to one line with the bare ID. **[source]**

### 1.4 `@RG` merge

`mergeReadGroups` (L213-238): same `mergeHeaderRecords` machinery, keyed by `SAMReadGroupRecord.equals` (ID + all attributes,
`SAMReadGroupRecord.java` L223-238); result sorted by ID. All 12 inputs carry the identical `@RG ID:HG002_WES.L1 ...` line
(sarek passes the same `-R` string to every chunk) so they collapse to a single `@RG`, `hasReadGroupCollisions()==false`, and record
`RG:Z` tags are untouched (`MergingSamRecordIterator.java` L134-141). The `@RG` attribute order is the input's order (parsed into a
`LinkedHashMap`; ID is always emitted first, see 1.6). **[source]** If attributes differed between chunks (e.g. a per-chunk `PU`),
IDs would be renamed `HG002_WES.L1.<base36>` and record `RG` tags rewritten; the golden has one `@RG` so this did not happen. **[inferred]**

### 1.5 Record merge order (`htsjdk/samtools/MergingSamRecordIterator.java`)

* Comparator: `hasMergedSequenceDictionary()==false` and sort order `coordinate`, so `getComparator()` returns
  `SortOrder.coordinate.getComparatorInstance()` = `SAMRecordCoordinateComparator` (L178-199). **[source]**
* Mechanism (L98-104, L125-153, L159-165): a `java.util.PriorityQueue<ComparableSamRecordIterator>` holding one peekable iterator per
  input (added in `--INPUT` order). `next()` polls the iterator whose *head record* is smallest, emits its head, re-offers the iterator
  if non-empty. `ComparableSamRecordIterator.compareTo` (`htsjdk/samtools/ComparableSamRecordIterator.java` L67-75) compares the two
  head records with `comparator.compare` (the **full** comparator, not `fileOrderCompare`). **[source]**
* This is a k-way merge of each file's existing sequence; records within one file are never reordered relative to each other.
  Reproduce it as a k-way merge, not as a global sort (per-file order from `samtools sort` need not agree with htsjdk's tie-breaks). **[source + inferred]**
  If every chunk happens to be internally ordered consistently with this comparator, the merge coincides with a global sort, so sorting
  the golden and comparing cannot validate the merge claim; only replaying the merge from the chunk BAMs (check 12) can.
* After polling, `record.setHeader(mergedHeader)` re-resolves reference/mate reference names against the merged dictionary (L132);
  since all dictionaries are identical, indices are unchanged. **[source]**

`SAMRecordCoordinateComparator.compare` (`htsjdk/samtools/SAMRecordCoordinateComparator.java` L47-71, L88-106) **[source]**, in order:

1. `fileOrderCompare`: `refIndex1 == -1` -> (refIndex2 == -1 ? 0 : 1); `refIndex2 == -1` -> -1; else `refIndex1 - refIndex2`;
   then `alignmentStart1 - alignmentStart2` (POS, 1-based; 0 for unplaced). `getReferenceIndex()` is the record's own RNAME index,
   so mate-placed unmapped reads (flag 0x4 with RNAME/POS set) sort at their placed position.
2. If strand flags differ: forward (0x10 clear) first.
3. Read name: `String.compareTo` (UTF-16 code-unit order = byte order for ASCII).
4. `getFlags()` ascending (int compare).
5. `getMappingQuality()` ascending.
6. `getMateReferenceIndex()` ascending (-1 for `*`).
7. `getMateAlignmentStart()` ascending (0 if none).
8. `getInferredInsertSize()` (TLEN) ascending.

Ties (all eight equal) between heads of different files are resolved by `PriorityQueue` heap structure (depends on `--INPUT` order).
They require the same read name in two inputs, which fastp chunking rules out, so output record order is independent of `--INPUT`
order **[inferred; matches observation]**. In particular, the unmapped-unplaced tail (refIndex -1 in every file) is merged by
steps 2-8 across the 12 file tails (each tail kept in its file's internal order).

### 1.6 Output header construction (pass 2)

`doWork` L283-308 **[source]**:

1. `outputHeader = header.clone()` (L287). `SAMFileHeader.clone()` (`SAMFileHeader.java` L396-400) round-trips through text:
   `getSAMString()` encodes with `keepExistingVersionNumber=false`, then decodes with validation SILENT. Encoding
   (`SAMTextHeaderCodec.encode` L411-433): `@HD` (L492-512: a fresh header gets `VN:1.6` first, then all other `@HD` attributes
   except `VN` in stored order), then `@SQ` (L518-525: `SN`, `LN`, then other attributes in stored order), `@RG` (L484-490: `ID`
   first, then attributes in stored order), `@PG` (L472-478: `ID` first, then attributes in stored order), `@CO`. Lines are
   tab-joined, `\n`-terminated. Decoding (L142-148 `transferAttributes`, L150-158 `parsePGLine`, L219-243 `parseHDLine`,
   `ParsedHeaderLine.mKeyValuePairs` is a `LinkedHashMap`) preserves attribute order (with the ID/SN/LN removed from the map and
   stored separately). Net effect: clone preserves order; `@HD` becomes `VN:1.6 GO:none SO:coordinate`.
2. `COMMENT` is empty -> nothing added (L305).
3. `getChainedPgIds(outputHeader)` (`AbstractMarkDuplicatesCommandLineProgram.java` L132-158): for every distinct value in
   `pgIdsSeen` (the set of `PG:Z` values seen on records in pass 1, including `null` for records without a PG tag;
   `MarkDuplicates.java` L529-536) create one `@PG`:
   * `newPgId = PgIdGenerator.getNonCollidingId("MarkDuplicates")` (`SAMFileHeader.java` L410-442): `"MarkDuplicates"` if no existing
     `@PG` has that ID, else `"MarkDuplicates." + base36(counter++)` where `counter` starts at the number of existing `@PG` records.
   * Attributes set in this order: `VN` = `PROGRAM_GROUP_VERSION` (defaults to `getVersion()`), `CL` = `PROGRAM_GROUP_COMMAND_LINE`
     (defaults to `getCommandLine()`), `PN` = `PROGRAM_GROUP_NAME` = `"MarkDuplicates"`, `PP` = the existing PG id; `setPreviousProgramGroupId(null)`
     **removes** the `PP` key (`AbstractSAMHeaderRecord.setAttribute(key, null)` L75-81). So for `pgIdsSeen == {null}` the line is
     `@PG\tID:MarkDuplicates\tVN:<version>\tCL:<command line>\tPN:MarkDuplicates` with no `PP`, appended after all merged `@PG`s (L152). **[source]**
   * `chainedPgIds` maps old PG id (or `null`) -> new id; it is a `HashMap`, so `null` is a valid key.
   * `pgIdsSeen` is a `HashSet`; with more than one distinct value the *assignment* of `MarkDuplicates`, `MarkDuplicates.<n>` to old ids
     follows `HashSet` iteration order (nondeterministic in principle). Not the case here (one value). **[source]**
4. `getVersion()` = Barclay `"Version:" + MarkDuplicates.class.getPackage().getImplementationVersion()` (`barclay/.../CommandLineArgumentParser.java` L197-199).
   Inside the GATK shadow jar the manifest `Implementation-Version` is `archiveVersion` = `4.6.2.0` (`gatk/build.gradle` L541-549), so
   `VN:Version:4.6.2.0` **[inferred - copy the literal from the golden]**.
5. `getCommandLine()` (Barclay L220-246) = `"MarkDuplicates"` + explicitly set arguments (in the parser's argument-definition order,
   repeated `--INPUT x` per value, in `--INPUT` order) + all non-null defaults. This string is used verbatim in `CL:` and in the
   metrics `# ` header line; it depends on `--INPUT` order (the second source of run-to-run nondeterminism). **[source]**

The BAM header block is written by `SAMFileWriterImpl.setHeader` (`htsjdk/samtools/SAMFileWriterImpl.java` L152-175, L267-271):
`header.setSortOrder(coordinate)` (re-put, order unchanged), text = `SAMTextHeaderCodec.encode(writer, header)` with
`keepExistingVersionNumber=false` (same encoder rules as above), then `BAMFileWriter.writeHeader` (`htsjdk/samtools/BAMFileWriter.java`
L176-186): magic `BAM\1`, `l_text` + text (no NUL terminator: `writeString(text, true, false)`), `n_ref`, then per `@SQ`:
`l_name` + name + NUL, `l_ref`. **[source]**

---

## 2. Pass 1: `buildSortedReadEndLists` (`MarkDuplicates.java` L475-634)

### 2.1 Collections

* `pairSort`, `fragSort`: `SortingCollection<ReadEndsForMarkDuplicates>` with `ReadEndsMDComparator(useBarcodes=false)` (L497-507);
  `maxInMemory = min(maxMemory * 0.25 / 69, Integer.MAX_VALUE - 32)` (L483; `ReadEndsForMarkDuplicates.SIZE_OF` = (1*1)+(5*2)+(5*4)+(8*2)+1+8+13 = 69, `ReadEndsForMarkDuplicates.java` L43-45).
  `MAX_RECORDS_IN_RAM` is overwritten at L482 and otherwise unused here; `MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP` is documented obsolete (L165-167). **[source]**
* `tmp` (unmatched first mates) = `DiskBasedReadEndsForMarkDuplicatesMap(MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000, codec)` because
  the sort order is coordinate (L512-513); semantically a map keyed by `(int partition, String key)` built on
  `htsjdk/samtools/CoordinateSortedPairInfoMap` (`put` L~120, `remove` L~64). Putting an existing key throws. **[source]**
* `libraryIdGenerator = new LibraryIdGenerator(header, false)` (L518-520).

### 2.2 Per-record loop (L524-626)

`index` starts at 0 and is incremented at L622 for **every** record delivered by the merged iterator (mapped, unmapped-placed,
secondary, supplementary, QC-fail), except that the loop `break`s (L544-548) at the first record with `flag & 0x4` **and**
`getReferenceIndex() == -1` (the first unmapped-unplaced record; in a coordinate-sorted merge these all come last). Records with
`0x4` set but a reference index (mate-placed unmapped) are skipped but still counted in `index`. **[source]**

For every record with `0x4` clear and `!(0x100 || 0x800)` (`isSecondaryOrSupplementary`, `SAMRecord.java` L1122-1124) (L551):

1. `fragmentEnd = buildReadEnds(header, index, rec, false)` (L653-706) **[source]**:
   * `read1ReferenceIndex = rec.getReferenceIndex()`.
   * `read1Coordinate = neg ? rec.getUnclippedEnd() : rec.getUnclippedStart()` where
     `getUnclippedStart = POS - (sum of leading S and H lengths, stopping at the first non-S/H op)`
     (`SAMUtils.java` L775-786) and
     `getUnclippedEnd = alignmentEnd + (sum of trailing S and H lengths, scanning from the last op backwards until a non-S/H op)`
     (L799-813) with `alignmentEnd = POS + Cigar.getReferenceLength() - 1` (`SAMRecord.java` L610-619; reference length sums
     `M, D, N, =, X`). Both soft and hard clips count. Example: `POS=1000, CIGAR=5S146M` forward -> 995; `POS=1000, CIGAR=3H140M8S` reverse -> 1000+140-1+8 = 1147.
   * `orientation = neg ? R(1) : F(0)` (`ReadEnds.java` L31).
   * `read1IndexInFile = index` (the merged-stream index above).
   * `score = computeDuplicateScore(rec, SUM_OF_BASE_QUALITIES)` (see 2.3).
   * `read2ReferenceIndex = rec.getMateReferenceIndex()` **only if** `0x1` set and `0x8` clear; otherwise stays -1 (L668-670).
     `isPaired()` is defined as `read2ReferenceIndex != -1` (`ReadEnds.java` L49).
   * `libraryId = libraryIdGenerator.getLibraryId(rec)` (see 2.4).
   * Optical location: `opticalDuplicateFinder.addLocationInformation(readName, ends)` returns false for our names (see 4.3), so
     `readGroup` stays -1, `tile/x/y` stay -1 (`PhysicalLocationInt.java` L7-8, `ReadEnds.java` L41).
   * `read2Coordinate`, `read2IndexInFile` stay -1; `orientationForOpticalDuplicates` stays -1.
2. `fragSort.add(fragmentEnd)` (L554) - **every** primary mapped record contributes exactly one fragment entry.
3. If `0x1` set and `0x8` clear (L556):
   * `key = rec.getReadGroup().getReadGroupId() + rec.getReadName()` (string concatenation, no separator; `getReadGroup()` resolves
     the `RG:Z` value against the merged header, `SAMRecord.java` L871-878; a missing RG would NPE) (L557-559).
   * `pairedEnds = tmp.remove(rec.getReferenceIndex(), key)` (L560).
   * **First mate seen** (null): `pairedEnds = fragmentEnd.clone()` (field copy, `ReadEndsForMarkDuplicates.java` L58-77, L85-87) and
     `tmp.put(pairedEnds.read2ReferenceIndex /* = mate ref index */, key, pairedEnds)` (L563-567).
   * **Second mate seen** (L568-617), with `A` = stored first mate, `B` = current record, `b.ref = fragmentEnd.read1ReferenceIndex`, `b.coord = fragmentEnd.read1Coordinate`:
     * `orientationForOpticalDuplicates` (L574-584): if `B` has `0x40`: `getOrientationByte(B.neg, A.neg)` else `getOrientationByte(A.neg, B.neg)`
       (i.e. always (strand of read 1 of the pair, strand of read 2 of the pair)). `A.neg` is read as `pairedEnds.orientation == R`.
     * If `b.ref > A.read1ReferenceIndex || (b.ref == A.read1ReferenceIndex && b.coord >= A.read1Coordinate)` (L587-603):
       `read2ReferenceIndex = b.ref; read2Coordinate = b.coord; read2IndexInFile = index; orientation = getOrientationByte(A.neg, B.neg)`;
       then, if `read2ReferenceIndex == read1ReferenceIndex && read2Coordinate == read1Coordinate && orientation == RF` set `orientation = FR` (L599-603).
     * Else (L604-613, only reachable in coordinate order when same contig and `b.coord < A.read1Coordinate`): read2 <- A's
       `(ref, coord, index)`, read1 <- B's `(ref, coord, index)`, `orientation = getOrientationByte(B.neg, A.neg)`. (No RF->FR fix in this branch; it cannot apply since coords differ.)
     * `getOrientationByte(r1neg, r2neg)`: FF=2, FR=3, RR=4, RF=5 (`ReadEnds.java` L31, L66-74).
     * `pairedEnds.score += getReadDuplicateScore(B)` (short arithmetic; L615, L646-648): pair score = scoreA + scoreB.
     * `pairSort.add(pairedEnds)` (L616).
   * "read1"/"read2" of a pair entry therefore means **lower / higher unclipped-5' position** (ties: first-seen is read1), not first/second of pair.
   * A first mate whose partner never arrives stays in `tmp`; it then has only its fragment entry (with `isPaired()==true`) and can never be flagged (see 4.2). **[source]**

### 2.3 Duplicate score (`htsjdk/samtools/DuplicateScoringStrategy.java`)

`computeDuplicateScore(rec, SUM_OF_BASE_QUALITIES)` (L61-63, L73-116) **[source]**:

* `sum = Σ q over all bytes of rec.getBaseQualities() with q >= 15` (L49-56). All bases of SEQ count (soft-clipped bases included; hard-clipped bases do not exist in SEQ). Qualities are signed bytes; a missing-quality record (`0xFF`, decoded as empty `NULL_QUALS`) contributes 0.
* `score = (short) min(sum, 16383)` (`Short.MAX_VALUE/2 = 16383`, L82).
* If `0x200` (fails vendor QC) set: `score += (short)(Short.MIN_VALUE/2) = -16384` (L109).
* Cached per record in a transient attribute (irrelevant for output).
* `MQ`/`ms` tags are **not** used by this strategy (only `TOTAL_MAPPED_REFERENCE_LENGTH` with `assumeMateCigar` looks at MC; never here). **[source]**
* Fragment entries carry the single read's score; pair entries carry `scoreA + scoreB` evaluated in Java `short` (no overflow possible: max 2*16383, min 2*(-16384) = -32768). **[source]**

### 2.4 Library ids (`picard/sam/markduplicates/util/LibraryIdGenerator.java`)

* Library name = `LB` of the `@RG` referenced by the record's `RG:Z`; `"Unknown Library"` if no RG tag, unknown RG, or no LB (L110-122).
* `getLibraryId` assigns `short` ids starting at 1 in order of first encounter (L125-135). One library here -> id 1 for all. **[source]**
* `metricsByLibrary` is a `TreeMap<String, DuplicationMetrics>` pre-populated with one entry per `@RG` library (L62-74). **[source]**

### 2.5 Spilling and determinism

* `SortingCollection.add` spills when `maxRecordsInRam` is reached (`htsjdk/samtools/util/SortingCollection.java` L166-197):
  `Arrays.parallelSort(ramRecords, comparator)` then codec-encode to a temp file (L243-267). Iteration merges spill files with a
  `TreeSet` ordered by `(comparator on peeked records, then file serial number)` (L502-592, L652-660); in-memory-only data is
  `parallelSort`ed (L450-458). **[source]**
* The comparator (3.1) ends with `read1IndexInFile`, `read2IndexInFile`; entries in one collection never share both (fragment entries
  have unique `read1IndexInFile`; pair entries have unique index pairs), so it is a **total order** and spill boundaries, `parallelSort`
  instability, and JVM heap size cannot change the iteration order. **[source]**
* The codec (`ReadEndsForMarkDuplicatesCodec.java` L52-112) round-trips `score, libraryId, orientation, read1ReferenceIndex, read1Coordinate,
  read1IndexInFile, read2ReferenceIndex, [read2Coordinate, read2IndexInFile only if orientation > R], readGroup, tile, (short)x, (short)y,
  orientationForOpticalDuplicates, duplicateSetSize`. Entries written while `orientation ∈ {F,R}` (the `tmp` map) have `read2Coordinate ==
  read2IndexInFile == -1`, which equals the decode defaults, so nothing is lost. `x`/`y` are already shorts (`PhysicalLocationShort.java` L38-41). **[source]**
* `tmp` (`CoordinateSortedPairInfoMap`) keeps the map of the reference index most recently queried in RAM and the others as per-index append files;
  semantics are those of an exact `(int, String) -> entry` map, so file-handle limits do not change results. **[source]**
* `SortingLongCollection` (duplicate indexes) sorts `long`s ascending; merge via `PriorityQueue` of peeked file iterators (L132-170, L238-255). Order of equal values is irrelevant. **[source]**

---

## 3. Sort orders

### 3.1 `ReadEndsMDComparator` (`MarkDuplicates.java` L982-1046) **[source]**

Used for both `pairSort` and `fragSort`. Compare in order (all via integer subtraction / `Long.compare`; all values are small so no overflow):

1. `libraryId`
2. `read1ReferenceIndex`
3. `read1Coordinate` (unclipped 5' position of the lower end; may be <= 0)
4. `orientation` (F=0 < R=1 < FF=2 < FR=3 < RR=4 < RF=5)
5. `read2ReferenceIndex` (-1 for unpaired fragments; mate's ref index for a paired read's fragment entry; higher end's ref for pair entries)
6. `read2Coordinate` (-1 for every fragment entry)
7. `tile`, 8. `x`, 9. `y` (-1 when the read name does not parse, e.g. SRA names; Illumina names make them decide score ties)
10. `read1IndexInFile`, 11. `read2IndexInFile` (-1 for fragment entries)

### 3.2 Chunk equality (`areComparableForDuplicates`, L810-833) **[source]**

* Pairs (`compareRead2=true`): `libraryId, read1ReferenceIndex, read1Coordinate, orientation, read2ReferenceIndex, read2Coordinate` all equal.
* Fragments (`compareRead2=false`): `libraryId, read1ReferenceIndex, read1Coordinate, orientation` equal. Note `read2ReferenceIndex`
  is **not** part of the fragment key, so a paired read's fragment entry and an unpaired fragment at the same 5' position/strand share a chunk.
* Chunks are formed by comparing each entry with the **first entry of the current chunk** (L748, L770); since the stream is sorted by a
  comparator whose leading keys are exactly the chunk keys, this equals grouping by key equality.

---

## 4. Duplicate selection (`generateDuplicateIndexes`, L739-797)

### 4.1 Pairs (L745-760, `handleChunk` L799-808, `markDuplicatePairs` L879-917) **[source]**

For each chunk of `pairSort` (in sorted order):

* size 1: `addSingletonToCount` -> `all_sets[1] += 1`, `non_optical_sets[1] += 1` (`AbstractMarkDuplicatesCommandLineProgram.java` L327-330).
* size > 1:
  1. `best` = the **first entry in sorted order whose score is strictly greater than every earlier entry's score**, with the first entry
     taken unconditionally (`end.score > maxScore || best == null`, `maxScore` initialised to 0; L880-889). Equivalently: the entry with
     the maximum score; ties broken by the comparator order of 3.1, which inside a chunk reduces to `(tile,x,y)` (all -1) then
     `read1IndexInFile`, then `read2IndexInFile` -> the pair whose lower-coordinate end appears **earliest in the merged stream** wins.
     Negative scores (QC-fail) are handled correctly by the `best == null` guard.
  2. `READ_NAME_REGEX != null`, so `trackOpticalDuplicates(list, best, finder, libraryIdGenerator)` runs (L891-893; see 4.3). With no
     parseable locations it flags nothing but still does `all_sets[n] += 1; non_optical_sets[n] += 1` (`trackDuplicateCounts` L362-377;
     `optical_sets` untouched because `optDupCnt == 0`).
  3. Every entry `!= best` (object identity): `addIndexAsDuplicate(read1IndexInFile)` and, if different, `addIndexAsDuplicate(read2IndexInFile)` (L895-903).
     Both mates of a non-best pair are always flagged.
  4. `TAG_DUPLICATE_SET_MEMBERS=false`: no DI/DS bookkeeping.

### 4.2 Fragments (L762-787, `markDuplicateFragments` L954-977) **[source]**

For each chunk of `fragSort` (key: library, ref, 5' coordinate, strand), with `containsPairs = any entry.isPaired()`,
`containsFrags = any !entry.isPaired()`:

* The in-loop call is skipped when `size <= 1` or `!containsFrags`; the trailing call after the loop is unconditional, but it is a no-op
  in those cases, so semantics are uniform:
* If `containsPairs`: **every unpaired entry** in the chunk is flagged (`addIndexAsDuplicate(read1IndexInFile)`), regardless of score;
  paired entries are never flagged here (L955-960). This is why a mate-unmapped read sharing a 5' position/strand with any paired primary end is a duplicate.
* Else (only unpaired entries): `best` chosen exactly as in 4.1 (max score, ties -> smallest `read1IndexInFile`, i.e. earliest in the
  merged stream); all others flagged (L961-975).
* Fragment chunks never touch the set-size histograms.

### 4.3 Optical duplicates with names like `SRR14724482.104665` **[source]**

* `READ_NAME_REGEX` defaults to the sentinel `DEFAULT_READ_NAME_REGEX` (`ReadNameParser.java` L28), so the optimized parser is used
  (`useOptimizedDefaultParsing`, L69-73, L86).
* `getLastThreeFields(name, ':', tmp)` (L160-186): scans from the end for `':'` or index 0; at `i == 0` it calls
  `rapidParseInt(name.substring(0, len))` = `rapidParseInt("SRR14724482.104665")`, which throws `NumberFormatException` because the first
  char `'S'` is not a digit (L194-219). `readLocationInformation` catches it, logs one warning, returns `false` (L125-132).
  (For a name with ':' but fewer than 5/7 fields it also returns false, L88-97.)
* `addLocationInformation` returns false (L135-152) -> `buildReadEnds` leaves `readGroup = -1`, `tile = x = y = -1` (`MarkDuplicates.java` L676-691).
* `findOpticalDuplicates` (`OpticalDuplicateFinder.java` L121-155): `keeperOrNull` returns null because `keeper.hasLocation()` is false
  (`PhysicalLocation.java` L34-36: `tile != -1`) (L318-325). For `n >= 3` the graph path (L206-288) adds no edges because no entry has a
  location (L222); for `n == 2` the fast path's `closeEnough` is false (`hasLocation` false, L328-335). Result: no optical duplicates,
  `READ_PAIR_OPTICAL_DUPLICATES = 0`, `optical_sets` histogram empty.
* `trackOpticalDuplicates` partitions by `orientationForOpticalDuplicates` only when a chunk contains both FR and RF (L280-321); any FF/RR
  value in such a chunk would throw, but a chunk's position-based `orientation` key makes this impossible (FR/RF chunks only contain
  `{F,R}` strand pairs). **[source]**
* `MAX_OPTICAL_DUPLICATE_SET_SIZE=300000` and `OPTICAL_DUPLICATE_PIXEL_DISTANCE=100` therefore have no effect on this data. **[source]**

### 4.4 Duplicate index collections

`duplicateIndexes` (and `opticalDuplicateIndexes`, allocated because `indexOpticalDuplicates = REMOVE_SEQUENCING_DUPLICATES || TAGGING_POLICY != DontTag`
= false -> **not** allocated; L272, L724-727) are `SortingLongCollection`s; `doneAddingStartIteration()` sorts ascending (L790-793). A given
record index is added at most once (a paired primary record is only ever flagged via its pair entry; an unpaired one only via its fragment entry). **[source]**

This at-most-once invariant is what makes "flag record i iff i is in the set of duplicate indices" equivalent to the Java pass-2 cursor,
which pulls **at most one** index per record (`nextIndexIfNeeded`, L935-946). If the invariant were ever violated (it is not in this
configuration; flow mode or barcodes are different code paths) the Java cursor would lag and mis-flag later records; a set-based Rust
implementation must not try to emulate that lag.

---

## 5. Pass 2: output (`doWork` L283-451)

### 5.1 Iteration and flagging (L338-428) **[source]**

`recordInFileIndex` counts every record of the merged stream (same numbering as pass-1 `index`, and it continues through the
unmapped-unplaced tail). For each record:

1. `addReadToLibraryMetrics(rec)` (see 6.1).
2. `nextDuplicateIndex` is advanced while `recordInFileIndex > nextDuplicateIndex` (coordinate branch of `nextIndexIfNeeded`, L935-946).
   `isDuplicate = (recordInFileIndex == nextDuplicateIndex)`.
3. `rec.setDuplicateReadFlag(isDuplicate)` (L351-357): 0x400 is **set** for duplicates and **cleared** for everything else (pre-existing 0x400 on input would be removed).
   Secondary (0x100), supplementary (0x800), unmapped (0x4, placed or not) and QC-fail records can never be in `duplicateIndexes`
   (pass 1 skips 0x4 and 0x100/0x800; QC-fail records are eligible but score-penalised). This is why 0 of the 2,102 supplementary
   records are flagged. `metrics.addDuplicateReadToMetrics(rec)` for duplicates (see 6.1).
4. `nextOpticalDuplicateIndex`: collection is null -> stays `NO_SUCH_INDEX`; `isOpticalDuplicate` false.
5. `CLEAR_DT=true`: `rec.setAttribute("DT", null)` (L365-367) - removes a `DT` tag if present; also forces attribute decoding and marks the
   binary block stale (see 5.4) on every record, even those without `DT`.
6. `TAGGING_POLICY=DontTag`: no `DT` written. `TAG_DUPLICATE_SET_MEMBERS=false`: no `DI`/`DS`. `BARCODE_TAG=null`: no `MI`.
7. `recordInFileIndex++`; `REMOVE_DUPLICATES=false`, `REMOVE_SEQUENCING_DUPLICATES=false`: every record is written.
8. `PROGRAM_RECORD_ID != null && ADD_PG_TAG_TO_READS=true` (`PGTagArgumentCollection.java` L10):
   `rec.setAttribute("PG", chainedPgIds.get(rec.getStringAttribute("PG")))` (L423-425) -> every record gets `PG:Z:MarkDuplicates`
   (records with no PG map through the `null` key). An input `PG` value is **replaced** (same tag -> value replaced in place, 5.3).
9. `out.addAlignment(rec)`.

### 5.2 Record order

Identical to the merged input order of 1.5 (writer is `presorted=true`; `SAMFileWriterFactory.makeWriter(outputHeader, true, OUTPUT, REFERENCE_SEQUENCE)`
L310 -> `makeBAMWriter` L251-253/L287-307 with `createIndex=true` and the default compression level; `initializeBAMWriter` L310-320).
Unmapped-unplaced reads are emitted at the end in the k-way-merged order of 1.5 (not "input order" across files). **[source]**

### 5.3 Tag list representation (`htsjdk/samtools/SAMBinaryTagAndValue.java`, `BinaryTagCodec.java`) **[source]**

* Binary tag key: `SAMTag.makeBinaryTag(tag) = (short)(tag.charAt(1) << 8 | tag.charAt(0))` (`SAMTag.java` L123-128), i.e. the
  little-endian 16-bit value of the two ASCII bytes; sort key is **(second char, first char)**.
* `BinaryTagCodec.readTags` (L269-306) builds the attribute list **sorted ascending by this key**: a tag greater than the current tail is
  appended, otherwise `head.insert(tmp)`; `insert` (L206-231) places the element before the first element with a larger key, or
  **replaces** an element with an equal key.
* `SAMRecord.setAttribute(tag, value)` (`SAMRecord.java` L1530-1532, L1582-1601): `value == null` -> `remove(tag)` (L233-240); otherwise
  `insert` as above. `find` (L242-246) stops early once it passes the key, relying on sorted order.
* Consequence: **output tags are always emitted in ascending key order, regardless of input order**, with `PG` inserted at its sorted slot.
  Keys for common bwa/MarkDuplicates tags: `SA`=0x4153 < `XA`=0x4158 < `MC`=0x434D < `MD`=0x444D < `PG`=0x4750 < `RG`=0x4752 <
  `NM`=0x4D4E < `AS`=0x5341 < `XS`=0x5358 < `DT`=0x5444 (removed) ... < `ms`=0x736D. So e.g. an input `NM MD MC AS XS RG` becomes
  `MC MD PG RG NM AS XS` on output.

### 5.4 Record re-encoding (`htsjdk/samtools/BAMRecordCodec.encode` L106-205, `BAMRecord.java`) **[source]**

Because every record has `setAttribute` called (step 5 and 8), `BAMRecord.setAttribute` (L214-218) decodes all attributes and sets
`mBinaryDataStale=true`, so `getVariableBinaryRepresentation()` returns null (L154-160) and the variable-length block is **re-encoded
from decoded fields** for every record (never copied from input bytes). Field by field:

* `block_size` recomputed (L121-139).
* `refID`, `pos = POS-1`, `l_read_name = len+1`, `mapq`, `bin`, `n_cigar_op`, `flag` (with 0x400 as decided), `l_seq`, `next_refID`, `next_pos = PNEXT-1`, `tlen` (L150-163).
* `bin` (L143-148, `SAMRecord.computeIndexingBin` L1695-1709): if `POS == 0` (unplaced) -> `bin = 0`; else `reg2bin(POS-1, end)` with
  `end = alignmentEnd` (1-based inclusive = 0-based exclusive) or `POS` (one base) when `alignmentEnd <= 0` (unmapped-placed, or empty
  reference span). htslib writes 4680 for unplaced reads; htsjdk writes **0** - the output differs from the input here. **[source; verify on golden]**
* Read name + NUL; CIGAR ops re-encoded (identical); bases re-packed; qualities as decoded (an input quality block whose first byte is
  `0xFF` decodes to the empty `NULL_QUALS`, `BAMRecord.decodeBaseQualities` L421-434, and is re-written as `0xFF` x `l_seq`, `BAMRecordCodec` L186-191).
* Tags written in list order (5.3) via `writeTag` (L186-227). Integer values are re-typed to the **smallest fitting type**
  (`getIntegerType` L152-181): `v > 2^32-1` error; `> 2^31-1` -> `I`; `> 65535` -> `i`; `> 32767` -> `S`; `> 255` -> `s`; `> 127` -> `C`;
  `>= -128` -> `c`; `>= -32768` -> `s`; else `i`. The input type byte is irrelevant: the output type is chosen by value alone
  (`0..127 -> c`, `128..255 -> C`, `256..32767 -> s`, `32768..65535 -> S`, larger -> `i`/`I`, negatives by range). So `NM:i:0` stored by
  htslib as `C` is written as `c`, and `AS:i:300` is written as `s` whatever its input type.
  Strings `Z`, chars `A`, floats `f` unchanged; `B` arrays keep element type and signedness (L229-260); `H` input would be re-written as `B:c`.
* Long-CIGAR `CG` handling (L113-123) is only for > 65535 ops - not applicable.

### 5.5 Index

`CREATE_INDEX=true` -> `BAMIndexer` writes `<OUTPUT minus .bam>.bai` (`BAMFileWriter.java` L101-125). Out of scope for byte-exactness.

---

## 6. Metrics file

### 6.1 Accumulation (`picard/sam/DuplicationMetrics.java`, `AbstractMarkDuplicatesCommandLineProgram.java`) **[source]**

Per record in pass 2, `addReadToLibraryMetrics(rec, header, libraryIdGenerator, false)` (L207-218) selects the library's row (creating
one for an unknown library) and calls `DuplicationMetrics.addReadToLibraryMetrics` (L264-275), first match wins:

1. `0x4` set -> `UNMAPPED_READS++` (placed and unplaced; secondary/supplementary flags are not checked first).
2. `0x100 || 0x800` -> `SECONDARY_OR_SUPPLEMENTARY_RDS++`.
3. `0x1` clear or `0x8` set -> `UNPAIRED_READS_EXAMINED++`.
4. else `READ_PAIRS_EXAMINED++` (record count; halved later).

For flagged records, `addDuplicateReadToMetrics` (L248-259): only if primary and mapped: unpaired/mate-unmapped -> `UNPAIRED_READ_DUPLICATES++`,
else `READ_PAIR_DUPLICATES++` (record count; halved later).

`finalizeAndWriteMetrics` (L170-205), per library row in `TreeMap` (name) order:

* `READ_PAIRS_EXAMINED /= 2`, `READ_PAIR_DUPLICATES /= 2` (integer division).
* `READ_PAIR_OPTICAL_DUPLICATES = (long) opticalDuplicatesByLibraryId.get(libraryId).getValue()` if that bin exists, else left at 0.
* `calculateDerivedFields()` (`DuplicationMetrics.java` L110-119):
  * `ESTIMATED_LIBRARY_SIZE = estimateLibrarySize(READ_PAIRS_EXAMINED - READ_PAIR_OPTICAL_DUPLICATES, READ_PAIRS_EXAMINED - READ_PAIR_DUPLICATES)`;
  * `PERCENT_DUPLICATION = (UNPAIRED_READ_DUPLICATES + READ_PAIR_DUPLICATES*2) / (double)(UNPAIRED_READS_EXAMINED + READ_PAIRS_EXAMINED*2)`,
    or `0.0` if the denominator is 0. Numerator/denominator are `long`; division in double.
* `metricsFile.addMetric(row)`; if exactly one library, `setHistogram(row.calculateRoiHistogram())`; then `addHistogram` of
  `all_sets`, `optical_sets`, `non_optical_sets` **in that order** (L195-202) - note `optical_sets` is added second.

`estimateLibrarySize(readPairs, uniqueReadPairs)` (L143-178), `f(x,c,n) = c/x - 1 + exp(-n/x)` (L183-185), all IEEE double **[source]**:

```
dups = readPairs - uniqueReadPairs
if readPairs <= 0 or dups <= 0: return null
m = 1.0; M = 100.0
if unique >= readPairs or f(m*unique, unique, readPairs) < 0: throw IllegalStateException
while f(M*unique, unique, readPairs) > 0: M *= 10.0
repeat 40 times: r = (m+M)/2; u = f(r*unique, unique, readPairs); if u == 0 break; if u > 0 m = r else M = r
return (long)(unique * (m+M) / 2.0)        // truncation toward zero
```

With the golden inputs (487639, 417649) this yields 1531991.303... -> 1531991 (verified by re-implementation). **[inferred: identical
results require `exp` to agree with Java's `Math.exp` for the ~45 evaluated arguments; see traps.]**

ROI histogram (`calculateRoiHistogram` L206-226, `estimateRoi` L198-200): for `x = 1.0, 2.0, ..., 100.0` (loop bound is 100, not the javadoc's 10):
`value = L * (1 - exp(-(x*pairs)/L)) / uniquePairs` with `L = ESTIMATED_LIBRARY_SIZE` (long -> double), `pairs = READ_PAIRS_EXAMINED`,
`uniquePairs = READ_PAIRS_EXAMINED - READ_PAIR_DUPLICATES`; bins keyed by `Double x`; `valueLabel = "CoverageMult"`, `binLabel` default `"BIN"`
(`Histogram.java` L42-43). Returns null (no histogram) if `ESTIMATED_LIBRARY_SIZE` is null.

Set-size histograms (`LibraryIdGenerator.java` L54-56): `Histogram<Double>("set_size", "all_sets")`, `("set_size", "optical_sets")`,
`("set_size", "non_optical_sets")`; incremented only from pair chunks (4.1), keys are `(double) chunkSize`. `optical_sets` stays empty.

### 6.2 File format (`htsjdk/samtools/metrics/MetricsFile.java`, `FormatUtil.java`) **[source]**

`write(File)` uses a `FileWriter` (platform charset) and `BufferedWriter.newLine()` (`\n` on Linux). Layout (L148-165):

```
## htsjdk.samtools.metrics.StringHeader
# <command line string, see 1.6 step 5>
## htsjdk.samtools.metrics.StringHeader
# Started on: <java.util.Date.toString(), e.g. "Wed Oct 08 10:11:12 UTC 2026">
<empty line>
## METRICS CLASS	picard.sam.DuplicationMetrics
LIBRARY	UNPAIRED_READS_EXAMINED	READ_PAIRS_EXAMINED	SECONDARY_OR_SUPPLEMENTARY_RDS	UNMAPPED_READS	UNPAIRED_READ_DUPLICATES	READ_PAIR_DUPLICATES	READ_PAIR_OPTICAL_DUPLICATES	PERCENT_DUPLICATION	ESTIMATED_LIBRARY_SIZE
HG002_WES	100	487639	2102	258	36	69990	0	0.143551	1531991
<empty line>
## HISTOGRAM	java.lang.Double
BIN	CoverageMult	all_sets	non_optical_sets
1.0	1	<all_sets[1]>	<all_sets[1]>
2.0	1.727382	<...>	<...>
...
100.0	3.66813	<...>	<...>
<rows for any set size > 100, e.g. "137.0	0	1	1">
<empty line>
```

Details:

* Headers (L168-177): `"## " + h.getClass().getName()` then `"# " + h.toString()`; `StringHeader.toString()` is the raw string (`StringHeader.java`).
  The second header is nondeterministic (wall clock). The first depends on `--INPUT` order (1.6 step 5).
* Metrics (L180-227): `"## METRICS CLASS\t" + class name`; column names are the public fields of `DuplicationMetrics` via
  `Class.getFields()` - declaration order (`DuplicationMetrics.java` L47-103; `MergeableMetricBase`/`MetricBase` declare no public fields)
  **[source for the fields; field order relies on HotSpot's reflection order, which is declaration order - inferred, matches golden]**.
  Values via `FormatUtil.format(Object)` (`FormatUtil.java` L54-70, L86-104): `long` -> `NumberFormat.getIntegerInstance()` with grouping
  off; `Double` -> `NumberFormat.getNumberInstance()` with grouping off, `maximumFractionDigits=6`, `RoundingMode.HALF_DOWN`, default
  minimum fraction digits 0 (so `1.0` prints `1`, `3.668130` prints `3.66813`), NaN/Infinity print `?`; `null` -> empty string.
  Locale-dependent decimal separator (`.` under the usual C/en locales). Rounding is applied to the exact binary value (JDK >= 8).
* Histogram (L230-273): empty histograms are dropped (`optical_sets` here), so there are exactly three value columns in the order
  `CoverageMult, all_sets, non_optical_sets`. The key set is the union of all non-empty histograms' keys (a `TreeSet<Double>` ascending).
  Header `"## HISTOGRAM\t" + key class name` = `java.lang.Double`; label row = first histogram's bin label (`BIN`) + each value label.
  Each row: `key.toString()` (**`Double.toString`**, so `1.0`, `2.0`, ..., `100.0`, `137.0`), then for each histogram
  `formatter.format(bin == null ? 0 : bin.getValue())` (so missing bins print `0`, counts print as integers).
* The file ends with the trailing empty line after the histogram (`out.newLine()` L160).

Golden-derived values (computed with the same algorithm; **verify**): `CoverageMult(1.0) = 0.99999997... -> "1"`, `(2.0) -> "1.727382"`,
`(3.0) -> "2.256465"`, `(10.0) -> "3.516052"`, `(100.0) -> "3.66813"`.

---

## 7. Other factors affecting exact output

* `--TMP_DIR`, `--MAX_FILE_HANDLES_FOR_READ_ENDS_MAP`, `--SORTING_COLLECTION_SIZE_RATIO`, `--MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP`,
  `--MAX_RECORDS_IN_RAM`: no effect on content (2.5). **[source]**
* `--VALIDATION_STRINGENCY LENIENT`: logging only. **[source]**
* `--REFERENCE_SEQUENCE`: only used for CRAM; BAM in/out ignores it. **[source]** (`SamReaderFactory.referenceSequence`, `makeWriter` reference arg)
* `--DUPLEX_UMI`, `--MOLECULAR_IDENTIFIER_TAG`, barcode tags: unused. **[source]**
* `--PROGRAM_GROUP_VERSION`/`--PROGRAM_GROUP_COMMAND_LINE` not given -> derived as in 1.6. **[source]**
* `--COMPRESSION_LEVEL 2`, `--USE_JDK_DEFLATER false`: BGZF only (out of scope). **[source]**
* Run-to-run nondeterminism, exhaustively: (a) `@PG` suffix assignment for colliding input `@PG`s (1.3); (b) the `CL:` string in the
  new `@PG` and the metrics `# ` line (`--INPUT` order); (c) `# Started on:` line. Record content/order and all metric numbers are deterministic. **[source + inferred]**
* `--ASSUME_SORTED false` + header `SO:coordinate` from `samtools sort` -> coordinate path everywhere; a queryname path exists in the code
  (`MemoryBasedReadEndsForMarkDuplicatesMap`, per-name indexing) but is not exercised. **[source]**

---

## 8. Rust implementation notes (summary of the algorithm to reproduce)

1. Open the 12 BAMs; build the merged header per 1.2-1.4 and 1.6 (IDs, suffixes, sorted `@PG`, appended `MarkDuplicates` `@PG` with
   `VN`/`CL` copied verbatim from the golden or regenerated per Barclay rules).
2. K-way merge the 12 record streams with the comparator of 1.5 (heap over heads; keep per-file order).
3. Pass 1 over the merged stream with a running `index`; stop at the first unplaced unmapped record; build fragment entries for all
   primary mapped records and pair entries per 2.2 using unclipped 5' coordinates, F/R orientation, capped Q>=15 sums, QC-fail penalty.
4. Sort both entry lists by 3.1; chunk by 3.2; select per 4.1/4.2; collect flagged indices.
5. Pass 2: same merged stream; set/clear 0x400 by index; drop `DT`; set/replace `PG:Z:MarkDuplicates`; emit tags in ascending
   `(c2<<8|c1)` key order with smallest-fitting integer types; recompute `bin` (0 for unplaced); accumulate metrics per 6.1.
6. Metrics per 6.1/6.2 with `f64` and the exact formatting rules.

---

## 9. Checks against the tier0 golden

Check first (cheapest discriminator of the biggest assumption): the `@PG ID:MarkDuplicates` line must have **no `PP`**. That proves
`pgIdsSeen == {null}`, i.e. no input record carried a `PG` tag. If a `PP` is present, section 1.3's tag rewriting applies, there is one
`MarkDuplicates`/`MarkDuplicates.<n>` `@PG` per distinct input `PG` value, and per-record `PG` values depend on `--INPUT` order.

Header:

1. `samtools view -H` line 1 is exactly `@HD\tVN:1.6\tGO:none\tSO:coordinate`.
2. `@SQ` lines are byte-identical (including attribute order) to those of the first `--INPUT` BAM (the first path in the `CL:` string).
3. Exactly one `@RG` line, byte-identical to the inputs' `@RG` line.
4. `@PG` IDs are exactly `bwa, bwa.1..bwa.9, bwa.A, bwa.B, samtools, samtools.D..samtools.N, MarkDuplicates` in that order (assuming each
   input has one bwa and one samtools `@PG` with distinct `CL`s); every `samtools*` line's `PP` names the `bwa*` ID whose `CL` lists the
   same fastq chunk; the MarkDuplicates line has tags in the order `ID VN CL PN` and **no `PP`**; `VN` is `Version:4.6.2.0`.
5. The two reference runs' `CL:` strings differ only in the order of the `--INPUT` arguments.

Records (use `samtools view` with tag-order-preserving output, e.g. check raw BAM or `samtools view` text which preserves order):

6. Every record carries `PG:Z:MarkDuplicates`; no record carries `DT`.
7. For every record, the tag list is strictly ascending in `(second char << 8 | first char)`; concretely `MC` < `MD` < `PG` < `RG` < `NM` <
   `AS` < `XS`, and `SA`/`XA` (if present) precede `MC`. Input chunk BAMs show a different order (e.g. `NM` first).
8. Integer tag type bytes in the raw BAM: `NM`/`AS`/`XS`/`MC`... values 0..127 use type `c`, 128..255 `C`, 256..32767 `s`, 32768..65535 `S`.
9. The `bin` field of every unplaced unmapped record (refID = -1, pos = -1) is 0 (not 4680); for mate-placed unmapped records it is
   `reg2bin(pos, pos+1)`; for mapped records `reg2bin(pos, pos + reference length)`.
10. No record with flag bit 0x4, 0x100 or 0x800 has 0x400; count of 0x400 records = 140016 = 36 + 2*69990.
11. Every flagged record with 0x1 set and 0x8 clear has its mate (same name, same RG) flagged too.
12. The output record sequence equals the k-way merge of the 12 input chunk streams under the comparator of 1.5 (replay it); in
    particular the unplaced-unmapped tail is **not** simply concatenated input order: at each step the smallest head by
    (strand, name, flags, mapq, mate ref, mate pos, tlen) is taken.
13. `samtools view -c -f 1 -F 0x90c` (primary, mapped, paired, mate mapped) = 2 * 487639 = 975278; `-F 0x904 -f 8` plus unpaired = 100;
    `-f 0x900 -F 4` = 2102; `-f 4` = 258 (includes mate-placed unmapped); total 977738.
14. Unclipped coordinate rule: pick a flagged/unflagged pair of reads where one has a leading soft clip, e.g. forward read `POS=p, CIGAR=5S...`
    groups with a forward read at `POS=p-5` with no clip (same strand, same mate end) - they must be in the same duplicate set.
    For reverse reads use `POS + reference length - 1 + trailing S/H`.
15. In every duplicate set (group of pairs sharing library, lower 5' ref/pos, orientation code, higher 5' ref/pos, where the orientation
    is computed from the lower end's strand then the higher end's strand with the same-position RF->FR rewrite), the unflagged pair has the
    maximum `Σ_{both mates} Σ_{q>=15} q` (capped per read at 16383, minus 16384 per QC-fail mate); among equal scores the kept pair is the one
    whose lower-coordinate mate appears first in the merged stream.
16. Every mate-unmapped primary record (`-f 9 -F 0x904`) that shares (ref, unclipped 5', strand) with any primary record of a mapped pair
    is flagged, even when its score is higher; among mate-unmapped-only groups the max score (earliest on ties) is kept. 36 such flags in total.
17. A pair whose mates are on different contigs forms a pair entry with read1 = the lower contig index; it participates in pair sets normally.

Metrics file:

18. Line 1 `## htsjdk.samtools.metrics.StringHeader`, line 2 `# MarkDuplicates --INPUT ...` identical to the `@PG CL:` string, line 3
    `## htsjdk.samtools.metrics.StringHeader`, line 4 `# Started on: ...` (nondeterministic), line 5 empty, line 6
    `## METRICS CLASS\tpicard.sam.DuplicationMetrics`, line 7 the column header, line 8
    `HG002_WES\t100\t487639\t2102\t258\t36\t69990\t0\t0.143551\t1531991`, line 9 empty, line 10 `## HISTOGRAM\tjava.lang.Double`,
    line 11 `BIN\tCoverageMult\tall_sets\tnon_optical_sets` (no `optical_sets` column), then rows keyed `1.0` ... `100.0` plus any
    set-size keys > 100, and a final empty line (file ends with `\n\n`).
19. Row `1.0` has `CoverageMult` = `1`; row `2.0` = `1.727382`; row `100.0` = `3.66813`; `all_sets` equals `non_optical_sets` on every row;
    `Σ_k k * all_sets[k] = 487639` and `Σ_k (k-1) * all_sets[k] = 69990`; rows with key > 100 have `CoverageMult` = `0`.
20. `ESTIMATED_LIBRARY_SIZE` recomputed by the bisection of 6.1 from `(487639, 417649)` is 1531991; `PERCENT_DUPLICATION` =
    `140016/975378` formatted to 6 decimals = `0.143551`.

Input assumptions to confirm on the chunk BAMs:

21. No input record carries a `PG` tag (otherwise 1.3's tag rewriting and multiple `MarkDuplicates.<n>` `@PG`s would apply).
22. All 12 `@RG` lines are byte-identical; all 12 `@SQ` sets are identical; each input has `@HD ... SO:coordinate`.
23. Reads carry an `RG:Z:HG002_WES.L1` tag (required: the pair key uses it, and a missing RG would crash Picard).
