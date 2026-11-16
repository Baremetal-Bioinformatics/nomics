# Mutect2 auxiliary models and filtering specification (GATK 4.6.2.0)

This specifies the rest of sarek 3.10.0's tumour-normal Mutect2 chain. Caller,
shared HC components, source commit identities and proposed port gates are in
[mutect2.md](mutect2.md). All GATK source references here use **tag 4.6.2.0**,
relative to `src/main/java/org/broadinstitute/hellbender/`. Class/method citations
cover executable behavior; comments are not an alternative specification.
**[golden]** uses the S1-r1 files and X basename defined in that companion spec.
The two S1-r1 scatter f1r2 archives and auxiliary-model goldens are available
under `~/omics/runs/somatic/S1-r1/golden/`; BAMs are not required for §1.

## 1. LearnReadOrientationModel

### 1.1 Inputs, defaults and canonicalization

`tools/walkers/readorientation/LearnReadOrientationModel.doWork` extracts each
`--input` tar to its own temporary directory, discovers `.ref_histogram`,
`.alt_histogram` and `.alt_table` leaves, and groups histograms by their first
MetricsFile StringHeader (sample). Ref and alt histogram sample sets and file
counts must agree when alt histograms exist. Tables are gathered by SAMPLE
metadata, preserving file order and row order within each list. Extraction is
Commons Compress `GzipCompressorInputStream` → `TarArchiveInputStream`, in member
order; discovery then uses **top-level `File.listFiles()`**, not a recursive walk
or sorted member list. Input directories are concatenated in **`--input` order**.
The collector writes one leaf of each kind per sample per archive, so filesystem
order within an archive cannot reorder records of the same sample in this layout.
With multiple leaves of a kind for one sample, `listFiles()` ordering is
filesystem-dependent and unspecified by Java. Defaults:
convergence threshold **1e-4**, maximum EM iterations **20**, max depth **200**.
The output name must end `.tar.gz`. **[source: that method, gatherAltSiteRecords,
sumHistogramsFromFiles; argument fields]**

Reference context is a 3-mer with one base of padding. There are 64 A/C/G/T
contexts and 32 reverse-complement pairs. Canonical representative is the
lexicographically smaller context. For each canonical context, sum its ref
histogram with its reverse complement. For alt-depth-one histograms, complement
the ALT base and swap orientation F1R2↔F2R1 before summing depth-bin counts.
For full alt records, append reverse-complemented records to the canonical list:
new refF1R2=refCount−oldRefF1R2; new altF1R2=altCount−oldAltF1R2. Preserve this
concatenation order. **[source: F1R2FilterConstants.ALL_KMERS/CANONICAL_KMERS;
LearnReadOrientationModel.combineRefHistogramWithRC,
combineAltDepthOneHistogramWithRC, mergeDesignMatrices;
AltSiteRecord.getReverseComplementOfRecord]**

Skip learning a context if its combined ref histogram has zero observations
**or its full alt design matrix is empty**. Having only singleton-alt histograms
does not avoid the second condition. ArtifactPriorCollection initially gives
every context a flat prior; skipped contexts retain that default and zero
example counts. **[source: LearnReadOrientationModel.doWork;
ArtifactPriorCollection constructor]**

### 1.2 Mixture states and beta priors

The ordinal order is F1R2_A,C,G,T; F2R1_A,C,G,T; HOM_REF; GERMLINE_HET;
SOMATIC_HET; HOM_VAR. Reverse complementation swaps artifact orientation and
complements the associated base; nonartifact states map to themselves.
**[source: readorientation/ArtifactState.java enum/getRevCompState]**

For a site with total depth n, alt count k, alt F1R2 count x, state s has
likelihood `BB(k;n,afAlpha_s,afBeta_s) × BB(x;k,oriAlpha_s,oriBeta_s)`.
The state probability multiplies this likelihood. **[source:
LearnReadOrientationModelEngine.computeLogPosterior/computeResponsibilities]**

| State | Allele-fraction Beta(alpha,beta) | Alt-orientation Beta(alpha,beta) |
|---|---|---|
| F1R2 artifact, matching ALT | 1,9 | 100,1 |
| F2R1 artifact, matching ALT | 1,9 | 1,100 |
| HOM_REF | 3,10000 | 10,10 |
| GERMLINE_HET | 5,5 | 10,10 |
| SOMATIC_HET | 2,5 | 10,10 |
| HOM_VAR | 10000,3 | 10,10 |

**[source: LearnReadOrientationModelEngine.getPseudoCountsForAlleleFraction,
getPseudoCountsForAltF1R2Fraction and constants]**

Ref→ref artifact states always have ln responsibility −Infinity. Artifact
states with a different ALT base also have −Infinity. `givenNotHomRef=true`
additionally excludes HOM_REF; learning passes **false**, application during
VCF filtering passes **true**. For reference-only sites ALT=REF and k=x=0,
so every artifact state is excluded. Normalize using NaturalLogUtils's exact
logSumExp order, not HC's approximate log10 Jacobian table.
**[source: computeResponsibilities; utils/NaturalLogUtils.normalizeFromLogToLinearSpace]**

### 1.3 EM and reduction order

`Histogram` uses a **TreeMap<Integer,Bin>**, so its bin iteration and
`getSumOfValues` reductions are ascending-depth ordinary additions. Histogram
files are merged in input order, by **value label**, into the first MetricsFile's
histogram list; counts are doubles (`Histogram.addHistogram` calls increment).
Integer-valued golden totals are exactly representable. No averaging occurs.
The context grouping and sample maps are HashMaps, but neither contributes a
floating reduction across contexts/samples. The engine's state-shape maps and
singleton-responsibility HashMap are queried by key, never summed by iteration.

Initialize mixture probabilities and regularizing pseudocounts with .1 for
each of the ten allowed states, zero for the two ref→ref artifact states.
The nominal pseudocount total is one; ordinary binary64 addition of ten .1
entries gives `0.9999999999999999`. Do not renormalize the initializer.
**[source: LearnReadOrientationModelEngine.getFlatPrior/learnPriorForArtifactStates]**

Each iteration takes an E step, an M step, then Euclidean distance(old,new).
Always execute at least one iteration; continue while distance **>1e-4** and
iteration count **<20**. This is parameter-distance convergence, not the
likelihood-change test suggested by an introductory comment.
**[source: learnPriorForArtifactStates]**

E step computes responsibilities for depths 1..200 for reference-only sites;
for each full alt record using its uncapped ref+alt depth; and for every
(depth,nonref base,orientation) singleton combination. M step separately sums:
full alt matrix in record order; each singleton histogram in **A,C,G,T excluding
REF, then F1R2,F2R1 order**, each in ascending depth; then ref histogram in
ascending depth. `sumArrayFunction` starts with the first vector and adds later
vectors component by component (ordinary additions, no compensated stream sum).
Add **(design + singleton) + reference**, add the .1 pseudocount vector, sum the
12 components in enum order with ordinary additions, and divide by that sum. No new beta shapes
are learned; EM learns only mixture weights. **[source: takeEstep/takeMstep]**

The beta-binomial implementation is
`tools/walkers/validation/basicshortmutpileup/BetaBinomialDistribution.logProbability`;
it returns `binomialCoefficientLog(n,k)+Beta.logBeta(k+alpha,n−k+beta)
−Beta.logBeta(alpha,beta)`, with k>n returning −Infinity. Negative k is an
argument error. Its beta/gamma and combinatorial arithmetic is a fresh-port
numerical boundary.
The default integer beta shapes mean all `Gamma.logGamma` calls reached through
`Beta.logBeta` reduce to arguments **1 or 2**, with `logGamma1p(0) = -0.0` and `logGamma1p(1) = +0.0`.
`logGammaSum` still uses the source's ordered `log1p`/`log` branches. Other Gamma
branches are not needed here. `FastMath.log` uses its 1024 split-mantissa table,
near-one polynomial or high-precision polynomial; `log1p` uses high-precision
`log(1+x)` and a residual, or the small-x cubic. Preserve every split and sum.
`binomialCoefficientLog` uses exact integer coefficients for n<67, rounded
multiplicative doubles for n<1030, otherwise ascending numerator log additions
then ascending denominator log subtractions (k is reduced by symmetry).

`computeLogPosterior` uses **(Math.log(prior) + logBB_AF) + logBB_orientation**;
do not pre-add the two likelihoods. `NaturalLogUtils.logSumExp` finds the first
maximum, starts sum at 1, then adds `Math.exp(value−maximum)` in state order
excluding that maximum and negative infinities. It returns maximum+Math.log(sum)
(or maximum+0 when sum==1); each responsibility is `Math.exp(logPosterior−logSum)`.
The port privately reuses nomics' fresh fdlibm log/exp and OpenJDK 17 decimal
digit generator; the complete golden and independent Java tests verify their
observed binary64/decimal behavior on this Mac. No math/compression/archive crate
is introduced. Large payload and likelihood buffers use nomics-mem.
**[source: computeLogPosterior; see §10 for dependency-level unconfirmed work]**

### 1.4 Output format

ArtifactPriorCollection.set stores both canonical and reverse-complement rows.
Writer iterates Java `HashMap.values()`, **not sorted contexts**. The table
starts `#<METADATA>SAMPLE=HCC1395_HCC1395_WES`, then tab-delimited columns:
`context, rev_comp, f1r2_a, f1r2_c, f1r2_g, f1r2_t, f2r1_a, f2r1_c, f2r1_g,
f2r1_t, hom_ref, germline_het, somatic_het, hom_var, num_examples, num_alt_examples`.
Commas here separate column names in prose; the file uses tabs. Doubles use
DataLine's Java Double.toString; example counts are integers. Output tar leaves
are URL-encoded sample + `.orientation_priors`; archive input is temp-dir
listFiles. `SequenceUtil.generateAllKmers` is a breadth-first A,C,G,T generator
which **removes the first complete AAA then appends it before returning**.
Thus ALL_KMERS starts AAC and ends AAA. Canonical contexts retain their first
encounter in that list; AAA is learned last. The prior HashMap starts with
capacity 64 and resizes to 128 on insertion 49. Its final bucket index is
`(String.hashCode() ^ (hash >>> 16)) & 127`, with insertion order within a bucket;
updates do not move entries. These String keys never treeify. The port models
this order explicitly rather than relying on a Rust map's iteration. **[source: ArtifactPriorCollection.set/writeArtifactPriors;
ArtifactPrior.ArtifactPriorTableWriter/column enum; utils/tsv/DataLine;
LearnReadOrientationModel.doWork; utils/io/IOUtils.writeTarGz]**

The golden archive has one `./HCC1395_HCC1395_WES.orientation_priors` member,
15,642 bytes and 64 rows, beginning ATT, CTT, GTT, TAT, AAA, CAA, AAC.
ATT has HOM_REF `0.999181460062392`, SOMATIC_HET `2.712306031466204E-4`,
4,434,773 examples and 86,357 alt examples. Its reverse-complement AAT has the
same example counts and reversed artifact-state probabilities. The port now
reproduces these values by learning from both original inputs (§1.5). **[golden: artifactprior tar payload]**

### 1.5 Port gate, archive contract and independent oracles

`nomics somatic learn-read-orientation-model --input CHR17 --input CHR1
--output priors.tar.gz` implements the defaults above. S1-r1 has one output
sample and one top-level regular ustar leaf of each input kind per sample per
archive. The port explicitly rejects multiple output samples, multiple leaves
of the same kind/sample within one input, nonregular/nested/extended tar entries,
and output names over 100 ASCII bytes. These regimes need separate gates:
GATK's multi-sample output `tmpPriorDir.listFiles()` order is filesystem-dependent;
the port does **not** invent or silently normalize an archive member order.
Custom max-depth/convergence/iteration CLI parameters are outside this default
command's scope. URL encoding is Java URLEncoder UTF-8 form encoding (space +,
uppercase percent hex; A–Z/a–z/0–9/`.`/`-`/`*`/`_` stay literal).

The gate compares **member names in order and each complete member's bytes**.
No table normalization is allowed. Gzip/tar envelopes are outside this gate:
GATK writes source-file timestamps/ownership and JVM compressor bytes. Our writer
emits deterministic POSIX ustar: 100-byte NUL-padded `./<encoded sample>.orientation_priors` name;
mode 0644; uid/gid 0; octal payload size; mtime 0; checksum computed with spaces
in its field, then six octal digits/NUL/space; regular type `0`; empty linkname;
magic `ustar\0`, version `00`; empty uname/gname/prefix; devmajor/devminor 0;
remaining reserved bytes zero. Payloads pad with zeros to 512 bytes; two zero
blocks terminate the archive. Gzip has mtime 0, flags/XFL 0, OS 255, nomics-io
raw DEFLATE level 2, CRC32 and ISIZE. These fields identify a portable regular
file without machine/user/time provenance; they are not attempts to match
GATK's incidental headers. The reader checks ustar checksums and gzip CRC/ISIZE.

`Mutect2FilteringEngine` extracts each prior archive using IOUtils, passes **all
extracted top-level files** to ReadOrientationFilter, and ArtifactPriorCollection
reads SAMPLE metadata plus canonical rows and reconstructs reverse complements.
It does not require archive mtimes, ownership or a particular compression level.
Keep exactly the prior leaf and no extra manifest. Compatibility is checked with
system tar, GNU tar 1.35 and GATK's pinned Commons Compress 1.26.0 stream classes;
this is not a full FilterMutectCalls rerun.

Independent numerical fixtures live beside the private implementation:
`orientation/Oracle.java` calls the pinned Commons Math 3.5 source's Beta and
CombinatoricsUtils directly and independently writes the 12-state responsibility
formula with Java Math operations. `EmOracle.java` separately performs E/M steps
for REF C; full rows (depth,alt,altF1R2,ALT)=(30,2,2,A),(30,2,0,A),(100,50,23,A);
1000 ref-only observations at depth 30 and 7 A/F1R2 singletons at depth 20.
The fixture has 1010 total / 10 alt examples and converges in ten iterations.
Build with javac 17 `-sourcepath /tmp/mutect2-spec-math3/src/main/java -d OUT
Oracle.java`, then compile EmOracle with `-cp OUT -d OUT`; run each with
`java -cp OUT Oracle` / `EmOracle`. Hex fixtures are Java doubleToLongBits,
not rounded decimal values or numbers generated by Rust. Unit tests compare
56 beta-binomial values (including n=0,66,1029,1200), four responsibility vectors,
20 FastMath log/log1p values and four small-integer logBeta values, and
**every intermediate prior** in that EM trajectory bit for bit. Additional
tests cover default skips, reverse complements, HashMap order, URL encoding and
archive integrity. The golden integration test skips with a message if the
S1-r1 directory is absent (`NOMICS_SOMATIC_GOLDEN` overrides it).

S1-r1 ordered inputs (chr17 before chr1) reproduce the complete 15,642-byte
`./HCC1395_HCC1395_WES.orientation_priors` member. Reversing those inputs
**can change bytes**: the full-alt rows remain canonical rows first, reverse
complements second, but groups swap within both lists, changing ordinary
floating sums in the M step. S1-r1 reversal changes 548 prior fields across all
64 rows. Counts and row order remain identical; neither sorting nor rounding is
permitted. The integration test runs both orders and
asserts this distinction. The release CLI took **1.226 s wall time** on this Mac
(including input/output);
system-tar and GNU-tar extraction followed by literal `cmp` both returned 0.
No GATK-relative speedup is claimed.

## 2. GetPileupSummaries and GatherPileupSummaries

### 2.1 Traversal and site selection

`tools/walkers/contamination/GetPileupSummaries` is a LocusWalker. Requires
reads, intervals and features; requiresReference=false even though this chain
supplies one. Its defaults are minimum population AF **0.01**, maximum **0.2**,
minimum mapping quality **50**. AF declaration must exist in the resource VCF
header. **[source: argument fields, requires* methods, onTraversalStart]**

At a locus take **only the first** resource VariantContext returned by the
FeatureContext. Require biallelic SNP and AF present with **.01 < AF < .2**,
both strict; resource FILTER is not consulted. Documentation describing
inclusive limits disagrees with `alleleFrequencyInRange`. Missing AF warns and
skips; seeing missing AF but no record with AF produces a traversal-success
error. **[source: apply/alleleFrequencyInRange/onTraversalSuccess]**

AF is the scalar INFO value via htsjdk 4.2.0
`CommonInfo.getAttributeAsDouble`, not the first item of an AF list.
`AbstractVCFCodec.parseInfo` stores a comma-separated value as a List;
the scalar getter then throws on its String cast. `AF=.` / empty / bare AF
also fail numeric parsing when visited, whereas an absent AF key warns/skips.
A malformed AF on an unvisited or non-biallelic-SNP record is not evaluated.
NaN and infinities parse but fail the strict range. Duplicate INFO keys use the
last value. htsjdk ignores no-call ALT `.` entries before testing biallelicity;
`*` is **nonsymbolic**, length one, and thus can qualify as a SNP. The later
BaseUtils wildcard mapping counts `*` as A. N alleles pass SNP selection but
cannot index the four-element count array: they fail if AF qualifies. No extra
ACGT-only site-selection rule is added. **[source: htsjdk
variant/variantcontext/{VariantContext.typeOfBiallelicVariant,CommonInfo.getAttributeAsDouble,
Allele.wouldBeSymbolicAllele}; variant/vcf/{AbstractVCFCodec.parseInfo/parseSingleAltAllele,
VCFUtils.parseVcfDouble}; GATK BaseUtils.simpleBaseToBaseIndex]**

FeatureContext uses overlapping records, including REF-length/INFO-END spans
starting before the locus, in resource order. The first record can therefore
hide a later SNP even if that first record fails the SNP/AF filters. **[source:
engine/FeatureContext.getValues; FeatureDataSource.queryAndPrefetch;
htsjdk AbstractVCFCodec.decodeLine]**

Default LocusWalker max depth/sample=**0** (no downsampling), include deletions
true, include N CIGAR skips false, emitEmptyLoci=false. There is no Mutect2
ITR transformer, assembly or f1r2 collection. Only covered loci are traversed;
a deletion-only covered locus can still have zero A/C/G/T counts.
**[source: engine/LocusWalker.defaultMaxDepthPerSample/includeDeletions/includeNs,
emitEmptyLoci; GetPileupSummaries inherits these defaults]**

There is **no minimum A/C/G/T count** for emission: at least one pileup element
is enough to traverse the locus. A deletion-only or N-base-only pileup emits a
zero-count row when its resource site qualifies; zero-coverage loci are absent.
`AlignmentStateMachine.stepForwardOnGenome` consumes I/S without advancing the
reference and ignores H/P; neither creates an independent pileup element.
M/=/X bases and D elements do. CIGAR N is excluded by Wellformed before traversal.
Mate overlap is not collapsed; both ends count independently after adaptor
filtering, regardless of base quality. **[source: utils/locusiterator/
{AlignmentStateMachine.stepForwardOnGenome,LocusIteratorByState.lazyLoadNextAlignmentContext};
ReadPileup.getBaseCounts; GetPileupSummaries.apply]**

The locus iterator additionally excludes `ReadUtils.isBaseInsideAdaptor` bases.
A well-defined fragment requires nonzero TLEN, paired mapped ends on opposite
strands, and reverse end > mate start or forward start ≤ mate start + TLEN.
The boundary is mate start−1 for reverse reads and start+abs(TLEN) for forward;
exclude reverse positions ≤ boundary or forward positions ≥ boundary.
The bypass uses **signed TLEN >100**, not abs(TLEN)>100. Proper-pair is not
required. **[source: LocusIteratorByState.dontIncludeReadInPileup;
utils/read/ReadUtils.{hasWellDefinedFragmentSize,getAdaptorBoundary,isBaseInsideAdaptor}]**

BED start is zero-based, end exclusive (`htsjdk BEDCodec.decode`: start+1,
end unchanged). Intervals are sorted in dictionary order, overlaps and adjacent
intervals unioned once; default padding is zero and merging rule ALL.
**[source: cmdline/argumentcollections/IntervalArgumentCollection fields;
utils/IntervalUtils.mergeIntervalLocations]**

### 2.2 Read filters, in order

1. MAPQ≥50.
2. MAPQ≠255.
3. MAPQ≠0.
4. Mapped.
5. Primary line: neither secondary **nor supplementary**.
6. Not duplicate.
7. Passes vendor check.
8. Nonzero reference-consuming alignment.
9. Mate on same contig **or** no mapped mate.
10. Good CIGAR.
11. WellformedReadFilter.

**[source: GetPileupSummaries.getDefaultReadFilters;
engine/filters/ReadFilterLibrary.PrimaryLineReadFilter and mate test;
MappingQualityReadFilter/WellformedReadFilter; shared good-CIGAR rules in caller spec]**

Wellformed also rejects missing BAM qualities: htsjdk
`BAMRecord.decodeBaseQualities` turns a first quality byte of 0xFF into an
empty array, so bases/quals lengths differ. This is distinct from a base-quality
threshold: quality zero is retained. HC's raw-slice quality-length test does
not account for this encoding; the new path does.

Wellformed's header agreement checks **start ≤ contig length**, not end.
A read extending beyond a contig end still passes this test; only requested
valid loci contribute. GoodCigar delegates to htsjdk `Cigar.isValid`, additionally
rejecting adjacent I/D and leading/trailing D ignoring clips. Padding may lead
a CIGAR; internal padding must be between real operators, and trailing padding
is invalid. HC's private helpers currently use an end-bound check and omit
some padding validity rules, so this port keeps source-correct local helpers
and leaves the existing gated HC behavior unchanged. **[source:
WellformedReadFilter.createFilter; ReadUtils.alignmentAgreesWithHeader;
utils/read/CigarUtils.isGood; htsjdk samtools/Cigar.isValid]**

PileupSummary calls ReadPileup.getBaseCounts: count A/C/G/T bases with no BQ
threshold, skip deletions and non-simple bases; ref/alt counts select the resource
bases, other_alt_count=sum(all four)−ref−alt. It does not collapse overlapping
mates or use posterior/fragment counts. Population AF is copied from resource.
Position is the resource VC start. **[source: PileupSummary(VariantContext,ReadPileup);
utils/pileup/ReadPileup.getBaseCounts; GetPileupSummaries.apply]**

### 2.3 File and gather semantics

TSV begins SAMPLE metadata followed by
`contig\tposition\tref_count\talt_count\tother_alt_count\tallele_frequency`.
Positions/counts are decimal integers; AF is Double.toString. Sample is the
first name from the sorted BAM sample set. **[source: GetPileupSummaries.onTraversalStart;
PileupSummary.PileupSummaryTableWriter; utils/read/ReadUtils.getSamplesFromHeader;
utils/tsv/DataLine]**

GatherPileupSummaries removes files with zero data rows, stably sorts remaining
**files by their first row** (dictionary contig index then start), then concatenates
whole tables. It assumes sorted, non-overlapping input files; it neither merges
interleaved rows nor deduplicates equal loci. First nonempty file supplies SAMPLE;
a differing sample is fatal. Equal first loci preserve argument order. Unknown
contigs compare with index **−1** (Java `List.indexOf`), before known contigs.
All-empty input writes just the column header, without SAMPLE metadata.
**[source: GatherPileupSummaries.doWork/removeEmptyFiles,
PileupSummaryFileComparator; PileupSummary.writeToFile(List<File>,File)]**

Golden normal has **246,147** rows, tumour **244,297**. First normal row is
`chr1 13273 6 0 0 0.115`; first tumour row is `chr1 17407 14 0 0 0.046`
(spaces here stand for file tabs). SHA256 of complete files:
normal `394fb78350724596cf04b2dec000378e4f3afc734d9a815cf1a652044757e738`;
tumour `d14d6e4bc76c1e82352d9b38d02dac70e235dadd5257c73643e1a0dd51a0edf0`.
**[golden: sample pileups]**

### 2.4 GetPileupSummaries port and chr21 gate

`nomics somatic get-pileup-summaries --input BAM --variant VCF.gz
--reference FASTA --intervals BED --output TABLE` implements the sarek defaults.
The reference bases are unused by this tool; `nomics_io::fasta::read_fai`
validates common dictionary lengths without loading the genome. BAM records
use the existing Header/RecordView/BGZF decoder; `nomics-io::bam_query::BamQuery`
adds BAI chunk union and virtual seeking. Exact read/BED overlap is checked
before the filter chain. Selected chunks are decoded once in file order,
without retaining BAM records. Tabix adds overlap queries with exclusive chunk
ends and chunk union, retaining earlier spanning features and duplicate records
in file order. The existing VCF fixed-column parser is exposed and reused;
existing start-only tabix queries are unchanged. Occupied 1-Mb VCF windows are
queried separately, with sparse counters retained for one contig at a time.
The existing Java Double.toString/table writer is shared with Gather.
There is no threading or depth cap in this traversal.

GATK source checked at tag 4.6.2.0, commit
`76edc75c26504da94bbaee66584e107e76ee15de`; its pinned htsjdk 4.2.0 at
`4cc010022ac038fb30f26e6f9717fabff3e808c1`. The read-only fixture is
`~/omics/runs/somatic/gps-chr21`; both complete output tables pass literal
`cmp` with **no normalization**. Tumour: **3,084** rows, **1,327,905** reads
processed; normal: **3,115** rows, **1,695,291** processed. Every ordered filter
count matches the corresponding GATK log. Nonzero counts, in filter order:
tumour MAPQ 155,946 / primary 704 / duplicate 268,210 / mate-contig 2,589;
normal 209,062 / 2,306 / 434,260 / 1,992. All seven other counts are zero.
Final release CLI wall time on this Mac: **1.77 s tumour**, **2.07 s normal**,
including indexed I/O and table writing.
No GATK-relative speedup is claimed.

`tests/pileup_goldens.rs` runs both byte and filter-count gates and skips with a
message when the BAM/resource/reference fixtures are absent
(`NOMICS_PILEUP_GOLDENS` / `NOMICS_PILEUP_REFERENCE` override paths).
The source-derived hand-built tests cover AF edges/missing/list values, feature
precedence, adjacent/overlapping BEDs, deletion-only/N-only/other-alt-only
coverage, zero-BQ bases, insertion/softclip offsets, overlapping mates,
adaptor boundaries, supplementary filters, contig ends and uncovered contigs.
Neither real chr21 golden contains a zero-total or other-alt-only row; those
edge expectations are established by source and the synthetic tests, not by
claiming evidence absent from the real tables.
The IO tests compare BAI chunk-union traversal to an independent positional
oracle across BGZF seams and check VCF earlier spans and identical duplicates.
These expectations are derived from the Java methods cited above; no GATK
execution of the synthetic cases is claimed. Full 13–14 GB sample gates remain
for the lead; only the supplied chr21 data were run here. Nondefault GATK
arguments and other input formats are outside this command's scope.

Verification (2026-11-12, pinned Rust 1.99.1):

* Before edits: `cargo test -p nomics-somatic -p nomics` and
  `cargo test -p nomics-io` passed, including existing somatic goldens.
* After edits: `cargo build --release -p nomics` passed;
  `cargo test -p nomics-somatic -p nomics-io -p nomics` passed, including
  existing contamination/orientation goldens and both new pileup tests.
* `cargo clippy --no-deps -p nomics-somatic -p nomics -p nomics-io -- -D warnings`
  fails on three **pre-existing** nomics-io diagnostics (also confirmed with
  `CARGO_TARGET_DIR=/tmp/nomics-gps-clippy-strict` to avoid cached checks):
  `bam.rs` manual_is_multiple_of, `bgzf.rs` type_complexity (writer worker type),
  `cram/container.rs` collapsible_if. No unrelated code or lint policy was changed.
  The same command with `-A clippy::manual_is_multiple_of
  -A clippy::type_complexity -A clippy::collapsible_if` appended passes.
  `cargo clippy --no-deps -p nomics-somatic -p nomics -- -D warnings` also
  passes. An existing dependency warning for HC's deprecated fetch_update
  appears in build/test/clippy output; HC is unchanged.
* `NOMICS_PILEUP_GOLDENS=/tmp/nomics-gps-absent-fixture cargo test
  -p nomics-somatic --test pileup_goldens pileup_goldens -- --nocapture`
  passes with the expected SKIP message.
* Final release CLI commands (each `cmp` exits **0**):

```sh
mkdir /tmp/nomics-gps-verified-20261112
for S in HCC1395_WES HCC1395BL_WES; do
  /usr/bin/time -p target/release/nomics somatic get-pileup-summaries --input "$HOME/omics/runs/somatic/gps-chr21/$S.chr21.bam" --variant "$HOME/omics/runs/somatic/gps-chr21/af-only-gnomad.chr21.vcf.gz" --reference "$HOME/omics/ref/GATK.GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta" --intervals "$HOME/omics/runs/somatic/gps-chr21/chr21.bed" --output "/tmp/nomics-gps-verified-20261112/$S.table" 2> "/tmp/nomics-gps-verified-20261112/$S.log" || exit
  cmp "/tmp/nomics-gps-verified-20261112/$S.table" "$HOME/omics/runs/somatic/gps-chr21/$S.chr21.pileups.table" || exit
  echo "$S cmp: PASS"
done
```

The filter logs were independently parsed and compared against each GATK log:
all 11 `(count, filter-name)` pairs and `(total filtered, total processed)`
match, in order. No normalization, new bioinformatics/compression dependency,
unsafe code, external publishing or changes to existing filter paths were added.

## 3. CalculateContamination

### 3.1 Coverage and matched-normal routing

`CalculateContamination.filterSitesByCoverage` first removes totalCount **≤10**.
From survivors calculate Commons Math Median and Mean coverage. Retain sites
with totalCount **>0.5×median AND <3×mean** (strict). Other ALT bases are part of
totalCount. Apply independently to tumour and normal; do not intersect first.
**[source: tools/walkers/contamination/CalculateContamination argument fields,
filterSitesByCoverage]**

Build the genotyping ContaminationModel from the matched normal. If segmentation
is requested, build a **separate tumour model** and write its segments. Compute
final contamination by applying normal-model homozygous-locus selection to the
tumour sites. The normal model's own contamination estimate is not the published
tumour estimate. **[source: CalculateContamination.doWork]**

### 3.2 Segmentation

`ContaminationSegmenter.findSegments` groups input sites with
Collectors.groupingBy(contig), then iterates HashMap **values**. Within each
contig select sites with observed ALT fraction in **[.1,.9]**, inclusive;
AF denominator is ref+alt+other. No het sites means no segments for that contig.
These hets locate boundaries; every coverage-filtered site inside each resulting
interval is brought back via OverlapDetector and sorted by position.
**[source: ContaminationSegmenter.findSegments/findContigSegments;
PileupSummary.getAltFraction]**

Kernel uses m=min(observed AF,1−AF),
`FastMath.exp(−(m1−m2)^2/(2×.025))`. Call KernelSegmenter with dimension **100**,
window sizes `[50]` (50 **sites**, not bases), max changepoints/chromosome **10**,
linear and log-linear penalties both **1.0**, sort order INDEX.
Boundary index belongs to the left segment. Prepend −1, append last het index;
segment starts at het[previous+1] and ends at het[current], both inclusive.
These are spans between observed hets, not whole-contig partitions.
**[source: ContaminationSegmenter constants/findContigSegments]**

`tools/copynumber/utils/segmentation/KernelSegmenter.findChangepoints` restarts
Java Random seed **1216** for each call. Low-rank construction samples min(100,N)
points **with replacement**, unless that equals N, in which case it uses the
original list directly. Compute symmetric sampled kernel; Commons Math SVD;
scale U columns by `1/(sqrt(singularValue)+1e-10)`; multiply full-vs-sampled
kernel by scaled U to form Z. Diagonal is **the square of the rounded row-vector norm**: plain sequential
sum of squared components, then sqrt, then square. It is not simply that sum.
**[source: calculateReducedObservationMatrix/calculateKernelApproximationDiagonal]**

Segment cost follows sequential D,W,V recurrence: diagonal sum D, row-vector
sum W, V updated by `2×dot(Znew,Wold)+diagonalNew`; C=D−V/size. Sliding window
costs compare left and right width-50 segments with their union, with **circular
wraparound at contig ends**. The loop updates and overwrites position 0 after
one full cycle; a direct independent cost recomputation need not match its low
bits. **[source: KernelSegmenter.calculateSegmentCost/calculateWindowCosts]**

Candidate minima use PersistenceOptimizer, remove endpoints, take up to ten
per window. Backward selection sorts/distincts candidates, merges the adjacent
pair maximizing oldCostLeft+oldCostRight−mergedCost (**first index on ties**),
and records total costs at every model size. Choose first minimum of cost plus
`k + k×ln(N/(k+1e-10))`, considering k=0..min(10,candidate count); return selected
indices sorted. The persistence ordering below and SVD math require numeric
oracles before port parity can be claimed. **[source: KernelSegmenter.findChangepointCandidates,
selectChangepoints/calculateChangepointPenalty; unconfirmed dependency detail §10]**

PersistenceOptimizer is in `tools/copynumber/utils/optimization/`, not the
segmentation package. Stably sort point indices by ascending cost; watershed
components are born at minima, extended into uncolored neighbors, and merged
when both neighbors are colored. Pair the higher-minimum component with the
current maximum; on equal minima the pairing branch selects the right component,
while component retention breaks ties by lower creation color. Recolor merged
boundaries and assign the maximum the left neighbor's color. Sort resulting
extrema pairs by decreasing persistence using a stable sort, prepend the global
minimum (first in sorted indices). The getter's comment promises increasing-index
ties, but the executable final sort compares **persistence only**; preserve its
inherited pair order rather than adding an index comparator. **[source:
PersistenceOptimizer constructor/findExtremaPairs/mergeComponents/getMinimaIndices]**

### 3.3 Minor-allele-fraction mixture model

`ContaminationModel` error rate e=`1.5×(sum(otherAlt)/sum(total))`; those reductions
use `mapToInt().sum()` before widening to long. Initial MAF=.5 per segment,
contamination c=0. Repeat **3** times: optimize each segment's MAF, select
non-LOH segments, then optimize contamination. **[source: constructor/calculateErrorRate]**

For population ALT frequency f and n=ref+alt, k=alt, genotype weights are
`[(1−f)^2, f(1−f), f(1−f), f^2]` in hom-ref, ALT-minor, ALT-major, hom-alt order.
Their sample ALT fractions are `[e/3, maf, 1−maf, 1−e]`.
Each weighted likelihood is weight×Binomial(k;n,(1−c)×sampleAF+c×f).
Site log likelihood uses **Java Math.log** of a plain sequential sum of the
four weighted likelihoods; segment likelihood sums sites with compensated
Java DoubleStream.sum. Model likelihood then adds segment likelihoods with
plain sequential `IndexRange.sum`, not another DoubleStream.sum. Other alts
affect filtering/error rate but not this binomial's n. **[source: genotypeLikelihoods/segmentLogLikelihood;
utils/MathUtils.binomialProbability]**

MAF maximization: bounds [.1,.5], initial .4, relative and absolute tolerance
.01, maximum evaluations **20**. Contamination model maximization: bounds
[0,.5], initial guesses **.02,.05,.1,.2**, relative/absolute tolerance 1e-4,
max evaluations **30** per guess; pick maximum objective, preserving list tie
behavior. These are OptimizationUtils/Commons Math Brent operations, not a
uniform grid or closed-form MAF estimator. **[source: calculateMinorAlleleFraction,
calculateContamination(double,...); utils/OptimizationUtils.max]**

Select non-LOH segments by decreasing minMaf from .40 in **double subtraction
steps .04**, requiring maf>threshold, until retained sites are **>25%** of all
segmented sites; otherwise use all segments. Preserve subtraction rather than
substituting a decimal literal table. **[source: getNonLOHSegments]**

### 3.4 Final estimate and standard error

Try minMaf=.40,.36,... while minMaf≥0 (same floating subtraction). Above .25
use model hom-alts; above .20 use hom-refs; otherwise heuristic hom-refs.
Model selection requires segment MAF>minMaf and posterior probability of the
requested homozygous genotype **>0.5**. Heuristic candidates have observed
ALT fraction **<.15**; keep those ≤max(.1,90th percentile of candidate AFs).
Intersect chosen positions with coverage-filtered tumour sites through
OverlapDetector. Stop when estimate is not NaN and error **<.2×estimate+.001**.
Final fallback uses heuristic hom-ref at minMaf=0; NaN gives estimate 0,error 1.
**[source: calculateContaminationFromHoms/getType/calculateContamination(Strategy,...)]**

Let d_i=total depth, o_i=ref count and q_i=1−populationAF for the hom-alt
strategy; for hom-ref use alt count and populationAF. Estimate tumour error e_T
from all retained tumour sites. Let D=sum d_i, O=sum o_i,
E=`Math.round(D×e_T/3)`, K=max(O−E,0), Q=sum d_i×q_i.
Contamination estimate is K/Q, capped at 1 only on return.
**[source: calculateContamination(Strategy,...)]**

Define a=sum q_i d_i and b=sum q_i(1−q_i)d_i².
Standard-error function s(c)=sqrt(a c(1−c)+b c²)/Q, or 1 if no sites.
Binary-search c−s(c)−estimate=0 on [estimate,1], and
c+s(c)−estimate=0 on [0,estimate], precision **1e-6**. Return the larger distance
from estimate to either root; if either root unavailable, use s(estimate).
`MathUtils.binarySearchFindZero` compares signs at both endpoints every step,
bisects until width≤precision, and returns the midpoint. **[source: same method;
utils/MathUtils.binarySearchFindZero]**

### 3.5 Tables and observed target

Contamination TSV columns: `sample\tcontamination\terror`, no SAMPLE metadata.
Golden single data row:
`HCC1395_HCC1395_WES\t0.0391300479306773\t0.00100753451566718`.
Whole file 86 bytes; SHA256
`dd545ad9aeef2dea786d32f8f4378b9aa387a3eb1660442d466c55644257cf86`.
**[source: ContaminationRecord.ContaminationTableWriter; golden]**

Segmentation TSV has SAMPLE metadata then `contig\tstart\tend\tminor_allele_fraction`.
Coordinates are inclusive observed-site spans; doubles use Double.toString.
Writer preserves ContaminationModel.segmentationRecords order, including the
HashMap contig order. Golden has **42** rows, 1,855 bytes; SHA256
`2a703e0520bbbf9603eba4b41578f0439d8a96d653a0c8ed2f3d953c85843ff3`.
It begins chr9 (three rows), chr7 (two), chrX, chr8, chr5, chr6; first MAF is
`0.45404420727340083`. Sorting this table is a byte mismatch.
**[source: MinorAlleleFractionRecord writer;
ContaminationModel.segmentationRecords; golden]**

### 3.6 Fresh-port verification and dependency arithmetic

`nomics somatic gather-pileup-summaries` and `calculate-contamination` reproduce
all four S1-r1 tables with literal `cmp`, no normalization. The recorded gather
argument order is chr17 group then chr1 group; stable whole-file sorting reverses
those two files. The actual GRCh38 `.dict` was used, not the `.fai` fallback.
Coverage filtering retains **147,029 normal** and **144,612 tumour** sites.
A Java oracle compiled from the actual four GATK model/segmentation classes and
Commons Math 3.5 independently reproduced all 42 tumour segment rows and the
published contamination/error pair. See
`crates/nomics-somatic/tests/oracle/README.md` for its scaffolding and numerical
fixture generation.

The fresh numerical boundary includes Commons Math FastMath exp/log lookup tables
and arithmetic, scaled hypot/scalb, the square-matrix SVD bidiagonalization and
four QR cases, saddle-point binomial PMF, and Brent maximization. SVD has **no
iteration cap** in upstream; Brent uses the stated evaluation budgets and its
default Integer.MAX_VALUE iteration budget, without a convergence checker.
Precision.equals in Brent allows adjacent representable points, not just ==.
Brent returns the best encountered point and preserves the older objective on
ties; the four-start contamination selection preserves the first optimum on a
tie. No Gamma or Beta special function is called by this chain's integer-input
PMF path: the Stirling-error table is used for small integer n/k, its asymptotic
polynomial for n/k≥15. Other users of Gamma/Beta in this document remain separate
port boundaries.

Java Math.log is distinct from Commons Math FastMath.log. Comparing all 144,612
filtered tumour rows at c=0, maf=.4 matched every genotype-likelihood bit but
found **5,953 platform-ln differences** in site log likelihood. Math.log and
StrictMath.log agreed in every such comparison. The reference log identifies
OpenJDK **17.0.17-internal** on aarch64. Porting OpenJDK 17.0.17 fdlibm
`__ieee754_log` removed a three-ulp error in the first segment's final MAF;
substituting either platform ln or FastMath.log is not permitted by the gate.
Java DoubleStream.sum keeps a compensated sum AND a simple sum for infinity
recovery; preserve the final subtraction of the compensation. Java Math.round
uses exponent/significand arithmetic, avoiding a rounded `x+0.5`. Output uses
the existing nomics OpenJDK-17 FloatingDecimal port, locally reused in this crate.
Commons Math LEGACY percentile divides p by 100 **before** multiplying by N+1.

Synthetic tests exercise no-matched-normal execution, no het segments, empty
inputs, strict coverage boundaries, sample mismatches, empty-file removal,
stable file ties, and retained overlapping loci. Those paths have no published
reference golden. Full tumour-only execution on the published tumour pileups
also ran (0.59 s): contamination `0.0019024981353576948`, error
`2.0702953973346697E-4`, exactly matching the independently compiled Java model;
its segmentation matched the matched-normal run's tumour segments. This is an
oracle comparison, **not a new published golden**. The source implementation
was inspected; final golden parity does not establish every intermediate SVD factor/optimizer trajectory for all
inputs. Java HashMap tree bins, insertion ties and resize splitting also have an
independent Java ordering oracle using deliberately colliding String keys.
Equal-position overlap-detector iteration ties, TSV quoting in unusual contig
names, and malformed numeric input are not validated by these fixtures.

## 4. FilterMutectCalls traversal and defaults

`tools/walkers/mutect/filtering/FilterMutectCalls` requires reference and reads
the merged VCF, inferred sidecar `<drivingVariantFile>.stats`, orientation tar,
contamination and segmentation tables. Missing .stats is fatal even though the
filtering-engine constructor has a no-stats path. onTraversalStart resets the
shared Utils RNG. No reads are supplied; all evidence is VCF annotations and
reference context. **[source: onTraversalStart/requiresReference]**

There are **four VCF traversals**, not two:

| Zero-based pass | Per record | At end |
|---|---|---|
| 0 | accumulate filter/model/threshold data | learn filter parameters, somatic model and threshold; clear accumulators |
| 1 | same with learned parameters | learn again and clear |
| 2 | accumulate again | learn **threshold only**; keep model/filter parameters fixed |
| 3 | calculate probabilities, write every record, collect output stats | write filteringStats TSV |

**[source: FilterMutectCalls.NUMBER_OF_LEARNING_PASSES=2,
numberOfPasses/nthPassApply/afterNthPass;
Mutect2FilteringEngine.accumulateData/learnParameters/learnThreshold]**

The golden FilterMutectCalls command confirms defaults below; the -1 mapping
quality sentinel is resolved to 30 in `getMinMedianMappingQuality`.

| Parameter | Value |
|---|---|
| threshold strategy / beta | OPTIMAL_F_SCORE / 1.0 |
| initial threshold / alternative FDR limit | .1 / .05 |
| max region / haplotype events | 3 / 2 |
| max ALT allele count | 1 |
| unique ALT read count | 0 |
| median BQ / MQ / position | 20 / 30 effective / 1 |
| max median fragment length difference | 10000 |
| max N ratio / min reads per strand / min AF | Infinity / 0 / 0.0 |
| default contamination estimate | 0.0, superseded by sample table |
| initial ln SNV / indel / real-vs-artifact prior | -13.815510557964275 / -16.11809565095832 / -2.302585092994046 |
| normal p-value / slippage length / slippage rate | .001 / 8 / .1 |
| haplotype distance / long indel length | 100 / 5 |
| modes | mitochondria=false, microbial=false |

**[source: M2FiltersArgumentCollection.java L20–268, getMinMedianMappingQuality,
getLogSnvPrior/getLogIndelPrior; golden filtered header]**

### 4.1 Filter list in construction order

The list below is exactly `Mutect2FilteringEngine.buildFiltersList`. A denotes
allele filter, V denotes variant filter (one probability copied to all alts).
Hard tests yield probability 0/1. All are included for this chain, including
ones disabled by default thresholds or absent required annotations.

| # | Class → FILTER | Kind, error type | Executable rule / required INFO |
|---|---|---|---|
| 1 | TumorEvidenceFilter → weak_evidence | A, SEQUENCING | learned sequencing-error posterior from TLOD and tumour AD; requires TLOD |
| 2 | BaseQualityFilter → base_qual | A, ARTIFACT | ALT MBQ<20; requires MBQ |
| 3 | MappingQualityFilter → map_qual | A, ARTIFACT | MQ<30; for signed indel length≥5 replace ALT MQ with REF MQ; requires MMQ |
| 4 | DuplicatedAltReadFilter → duplicate | A, ARTIFACT | unique count≤0; requires AS_UNIQ_ALT_READ_COUNT (absent in this golden) |
| 5 | StrandArtifactFilter → strand_bias | A, ARTIFACT | learned strand mixture, §6; requires AS_SB_TABLE |
| 6 | ContaminationFilter → contamination | A, NON_SOMATIC | §7; requires POPAF |
| 7 | StrictStrandBiasFilter → strict_strand | A, ARTIFACT | default 0 returns no probabilities; when enabled flags either strand count==0 (does not compare with the supplied positive magnitude); requires AS_SB_TABLE |
| 8 | ReadPositionFilter → position | A, ARTIFACT | ALT MPOS>-1 and <1; negative positions explicitly exempt; requires MPOS |
| 9 | MinAlleleFractionFilter → low_allele_frac | A, ARTIFACT | maximum tumour sample AF for allele<0; none missing defaults to 1; no required INFO |
| 10 | NormalArtifactFilter → normal_artifact | V, ARTIFACT | §7.2; requires NALOD,TLOD |
| 11 | NRatioFilter → n_ratio | V, ARTIFACT | NCount/all-sample ALT AD≥Infinity; zero ALT count skips; requires NCount (absent) |
| 12 | PanelOfNormalsFilter → panel_of_normals | V, ARTIFACT | PON attribute present; no required INFO |
| 13 | ReadOrientationFilter → orientation | V, ARTIFACT | added because prior tar provided; §6.2 |
| 14 | ClusteredEventsFilter → clustered_events | V, ARTIFACT | max ECNTH>2 OR ECNT>3; requires both |
| 15 | MultiallelicFilter → multiallelic | V, ARTIFACT | more than one TLOD value exceeds **5.0 in natural logs**; requires TLOD |
| 16 | FragmentLengthFilter → fragment | V, ARTIFACT | abs(MFRL[1]−MFRL[0])>10000; only first ALT; requires MFRL |
| 17 | PolymeraseSlippageFilter → slippage | V, ARTIFACT | §7.3; requires RPA,RU |
| 18 | FilteredHaplotypeFilter → haplotype | V, ARTIFACT | previous-pass neighboring phase probabilities, §8; no required INFO |
| 19 | GermlineFilter → germline | V, NON_SOMATIC | §7.1; requires TLOD,POPAF |

**[source: each class's errorType/filterName/requiredInfoAnnotations and
calculateErrorProbabilityForAlleles, calculateErrorProbability, areAllelesArtifacts
or isArtifact; buildFiltersList]**

htsjdk `VariantContext.getIndelLengths` returns **ALT.length−REF.length**,
so only insertions of length≥5 use reference MAPQ; deletions do not, despite
comments describing long indels generally. MultiallelicFilter obtains LODs through getTumorLogOdds,
which multiplies decoded log10 TLOD by ln(10), then compares with literal 5.0;
it does not compare VCF TLOD with 5. **[source: MappingQualityFilter.areAllelesArtifacts;
MultiallelicFilter.isArtifact; htsjdk variant/variantcontext/VariantContext.getIndelLengths]**

Required-annotation absence differs by filter base class: allele filter returns
an **empty list**, variant filter returns zero replicated to all alts. Both
clamp finite precision errors with max(min(p,1),0); NaN is not specially repaired.
All probabilities are calculated before accumulation/learning for that record.
**[source: Mutect2AlleleFilter.errorProbabilities;
Mutect2VariantFilter.errorProbabilities; ErrorProbabilities constructor]**

## 5. Somatic clustering and priors

`tools/walkers/mutect/clustering/SomaticClusteringModel` creates prior-map keys
−10..10, zero for SNV and signed length for indels. Initial values come from
§4's ln priors. A queried length outside the map is inserted with the current
minimum prior. That inserted value persists: EM updates only keys −10..10,
so a long-indel key first visited during pass 0 keeps its initial value even
after shorter-indel priors change. Recomputing the minimum on every lookup
is incorrect (the first visible S1-r1 mismatch was GERMQ at chr1:1670309). Per-specific-SNV prior adds **−Math.log(3.0)** (`MathUtils.LOG_ONE_THIRD`);
indel priors do not. Do not evaluate Math.log(1.0/3.0): the division rounds
first, and its logarithm has bits `bff193ea7aad030b`, whereas the source
constant has bits `bff193ea7aad030a`. This correction is separate from the
strand scalar-function gap described in §9.3.
The callable statistic supplies the denominator for empirical site priors.
**[source: constructor/getLogPriorOfSomaticVariant/indelLength]**

Initialize beta-binomial background Beta(1,1), high-AF Beta(10,1), log weights
**ln(1+.01), ln(.01)**. The first is `Math.log1p(.01)`, not ln(.99); preserve the
code even though these are not normalized initial weights. Splitting background
later similarly uses `Math.log1p(fraction)` for the new background weight.
**[source: constructor/initializeClusters]**

Each nonsymbolic ALT contributes Datum(TLODln,artifact probability,non-somatic
probability,ALT AD,total tumour AD,signed indel length). Symbolic ALT AD is zeroed
before total. Artifact probability>.9 increments obvious-artifact count and
skips datum; non-somatic probability>.9 just skips it. Datum's stored
non-sequencing-error probability is `1−(1−artifact)×(1−nonSomatic)`.
**[source: SomaticClusteringModel.record; clustering/Datum constructor]**

For datum d, corrected cluster likelihood is TLODln plus the beta-prior evidence
correction. With C(a,b)=logGamma(a+b)−(logGamma(a)+logGamma(b)), correction from
flat prior to Beta(a,b) is
`C(a,b)−C(a+k,b+n−k)−C(1,1)+C(1+k,1+n−k)`.
Sequencing-error probability uses lnSumExp(logWeight+correctedLikelihood) and
prior of somatic variant (§5.2). Probability used as learning weight is
`(1−artifact)×(1−storedNonSequencingError)×(1−sequencingError)`; this means the
artifact factor enters twice through Datum. Do not simplify away that behavior.
**[source: BetaBinomialCluster.correctedLogLikelihood/logOddsCorrection;
SomaticClusteringModel.probabilityOfSequencingError/probabilityOfSomaticVariant]**

### 5.1 Initialization and learning

Only the first learn initializes extra clusters, adding at most **5** discrete
AF clusters. Compute somatic-probability-weighted AF quantiles, target **50**
quantiles; after crossing a target, append that AF once and advance past all
crossed thresholds, then distinct. Calculate background-weighted binomial
posterior densities `(n+1)×Binomial(k;n,f)` at those AFs, and integrate between
local minima with trapezoids. Choose first largest-mass peak from a stable
sort. Stop if no peaks or its AF is below quantile index min(5,last).
**[source: initializeClusters/calculateAlleleFractionQuantiles,
calculateQuantileBackgroundResponsibilities/calculatePeaksAndMasses]**

Split fraction=min(.9,peakMass/totalMass) from background, add fuzzy binomial
cluster and run **5 EM iterations without prior updates**. Score
weighted log likelihood−(2×clusterCount)×ln(effectiveSomaticCount). This is the
source's BIC-like score, without substituting another factor-of-two convention.
If score decreases, remove the new cluster and restore old log weights; other
cluster parameters learned during the trial are not restored. First comparison
starts at −Infinity. **[source: initializeClusters]**

Every learn pass then does **5 EM iterations with prior updates**. Per datum,
calculate current somatic probability and cluster responsibilities from actual
beta-binomial count likelihoods. Add pseudocount **1** to each cluster's total,
normalize weights, and learn shapes. Update real-vs-artifact prior to
ln((variantCount+1)/(variantCount+technicalArtifactCount+2)). Technical artifacts
include obvious-artifact count plus retained Datum artifact probabilities.
With callable sites, each −10..10 prior becomes
ln(max(expectedCount/callable,1e-8 for SNV or 1e-9 for indels)). Clear data and
obvious-artifact count after five iterations. **[source: learnAndClearAccumulatedData,
performEMIteration]**

Exact reduction contract (JDK 17): the per-iteration signed-length counts are
fresh `HashMap<Integer,MutableDouble>` values, initially −10..10 from the
**two-argument** Collectors.toMap. Its accumulator uses putIfAbsent, appending
collision chains; newly encountered lengths also use putIfAbsent. Ordinary
`MutableDouble.add` accumulates somatic weights in retained VCF/ALT order,
not a compensated reduction. The initial capacity grows from 16 to 32; further
resizes split chains while preserving order. For the initial keys, iteration
starts `[-1,0,-2,1,-3,2,...,-10,9,10]`. Count-map values are then reduced with
compensated DoubleStream.sum to get variantCount; that total affects the
real-vs-artifact prior, **not** the empirical per-length denominators.
Technical artifacts use a separate compensated stream sum plus the obvious
count. Cluster responsibility totals accumulate in the same datum order with
ordinary array additions; add 1 per cluster, then normalize with the ordinary
sequential MathUtils.sum, divide, and Math.log. Both normalizations use the
whole array, never pairwise log-sum-exp. **[source: performEMIteration;
MathUtils.normalizeSumToOne; NaturalLogUtils.logSumExp/normalizeLog;
JDK 17 Collectors.uniqKeysMapAccumulator and HashMap.putVal/resize]**

Quantile total weights, peak masses, fuzzy-binomial alt/total counts and all
strand M-step reductions use compensated DoubleStream.sum. Initialization's
BIC likelihood and effective count use ordinary sequential MathUtils.sum.
Beta-binomial shapes update alpha and beta after computing **both** gradients
from the old pair, for ten epochs in datum order. EM updates weights and priors
before learning shapes from its saved responsibilities. Initialization runs once;
its frozen somatic-probability array is reused for background responsibilities
and BIC, but quantile weights are recomputed at each split. Learned clusters,
weights, standard priors, artifact prior and first-access long-indel priors
persist into the second learning pass; only data and obvious count clear.
The third traversal learns only the threshold; the fourth writes output.
**[source: initializeClusters/learnAndClearAccumulatedData;
BetaBinomialCluster.learn; BinomialCluster.learn; Mutect2FilteringEngine.learnParameters;
FilterMutectCalls.afterNthPass]**


BetaBinomialCluster learning makes **10 epochs**, in original datum order,
gradient rate **.01**. Gradients are the beta-binomial digamma derivatives;
alpha=max(alpha+(.01×gradientAlpha)×responsibility,1),
beta=max(beta+(.01×gradientBeta)×responsibility,.5), both gradients calculated
before updating the two values for that datum. The declared MAX_RATE=.1 is not
used to adapt that update. **[source: BetaBinomialCluster.learn]**

A BinomialCluster is actually a narrow beta-binomial: weighted ALT and total
counts each add **.0001**, mean=ALT/total then cap mean≤.99.
alpha+beta=(1−mean)/(mean×.01²)−1; alpha=mean×(alpha+beta).
This fixes SD/mean=.01. **[source: BinomialCluster.learn/getFuzzyBinomial]**

### 5.2 Common posterior calculation

Given ln odds L of real vs error and ln prior s of real, normalize the pair
`[L+s, log1mexp(s)]` and return its error component. Normal-artifact calculation
uses the learned real-vs-artifact prior instead of the genome-wide somatic prior.
NaturalLogUtils.log1mexp switches at ln(.5): log1p(−exp(x)) below, otherwise
log(−expm1(x)); x=0 returns −Infinity. Preserve Math vs FastMath call sites and
Java stream/loop reduction order. **[source: filtering/Mutect2FilteringEngine.
posteriorProbabilityOfError/posteriorProbabilityOfNormalArtifact;
utils/NaturalLogUtils.log1mexp/normalizeLog]**

## 6. Strand and read-orientation artifacts

### 6.1 StrandArtifactFilter

Use AS_SB_TABLE across all samples. For each ALT sum all allele forward and
reverse depths as nF/nR and select ALT forward/reverse kF/kR. Zero ALT support
or |indelLength|>**4** gives probability 0. Sequencing beta prior is (1,1000)
for substitutions, (1,5000) for indels shorter than **3**, otherwise (1,50000).
Artifact strand prior starts Beta(1,20), mixture probability **.001**.
**[source: StrandArtifactFilter constants/calculateArtifactProbabilities]**

Forward artifact likelihood is BB(kF;nF,alphaStrand,betaStrand) ×
BB(kR;nR,sequencing prior); reverse swaps directions. Nonartifact likelihood is
`choose(nF,kF)×choose(nR,kR)/choose(nF+nR,kF+kR)` times
BB(kF+kR;nF+nR,1,1). Priors are artifactPrior/2 each, 1−artifactPrior for none.
The code converts ln terms to log10 and calls **MathUtils.normalizeLog10**;
return forward+reverse responsibilities. Preserve Java **Math.pow** for every
base-10 exponent and **Math.log10** for the array log-sum; platform pow/log10
are not bit-equivalent. The port reuses its own OpenJDK scalar implementations,
not a numerical library. **[source: strandArtifactProbability,
artifactStrandLogLikelihood/nonArtifactStrandLogLikelihood]**

Accumulate ESteps for every annotated allele. Learning uses potential artifacts
with total artifact responsibility **>.1**. Prior becomes
(sumPotentialArtifact+1)/(sumPotentialArtifact+1+sumAllNonartifact+1000).
Empirical beta mean uses responsibility-weighted strand ALT/depth, with initial
pseudocounts 1,20. Fix beta=(1/mean−1)×alpha and optimize alpha over [.01,100],
initial **1**, tolerances .01/.01, max evaluations **100**, with weighted
artifact-strand log likelihood objective. Clear ESteps after learning.
**[source: accumulateDataForLearning/learnParameters/clearAccumulatedData]**

### 6.2 ReadOrientationFilter

Load extracted `.orientation_priors` tables by sample. SNP or MNP only;
no prior for sample, hom-ref genotype, missing F1R2/F2R1, or null/N-containing
3-mer gives 0. Select ALT with largest TLOD (first max). For an MNP, test each
base independently and take the **maximum** probability, using the same ALT
orientation counts for every base. Ref count=F1R2[0]+F2R1[0]; selected ALT count
and F1R2 count use its allele index. Call §1's responsibilities with
**givenNotHomRef=true** and return the **maximum**, not sum, of that ALT's
F1R2 and F2R1 state probabilities. The previous version of this specification
incorrectly said to add them; `ReadOrientationFilter.artifactProbability` calls
`Math.max(posteriorOfF1R2, posteriorOfF2R1)`.
**[source: ReadOrientationFilter.calculateErrorProbability/artifactProbability]**

Across tumour samples use ALT-depth-weighted median posterior: sort ascending
posterior, return first whose cumulative ALT depth×2≥total. Equal-half selects
the lower posterior; zero total returns the first available posterior through
the same test, empty input returns 0. Weight is total ALT AD, not orientation
count depth. The single tumour makes aggregation trivial here.
**[source: Mutect2Filter.weightedMedianPosteriorProbability;
ReadOrientationFilter.calculateErrorProbability]**

## 7. Contamination, germline, normal-artifact and slippage probabilities

### 7.1 GermlineFilter

Select ALT with max TLOD, first tie. POPAF is converted back by 10^(−POPAF).
Population AF<1e-10 gives 0, >1−1e-10 gives 1. Sum tumour AD; zero total gives 0.
MAF is depth-weighted over tumour sample segments overlapping the VC, default
.5 where no segment exists; if multiple overlap use the first detector result.
Tumour AF is similarly weighted from FORMAT AF by total sample AD.
**[source: GermlineFilter.calculateErrorProbability/computeMinorAlleleFraction;
Mutect2FilteringEngine.weightedAverageOfTumorAFs]**

Het germline count likelihood is .5×Binomial(k;n,maf)+.5×Binomial(k;n,1−maf).
Subtract learned somatic count log likelihood. Hom-alt-vs-somatic log odds is
−Infinity for tumour AF<**.9**, otherwise 0. Negate selected NLODln to get
normal-presence log odds; missing NLOD gives 0.
Germline het unnormalized log probability =
ln(2f(1−f))+hetVsSomatic+normalPresence+log1mexp(logSomaticPrior).
Hom-alt uses ln(f²) and homAltVsSomatic instead. Sum those two in log space;
somatic probability uses ln((1−f)²)+logSomaticPrior. Normalize germline vs
somatic, return germline component. **[source: calculateErrorProbability/germlineProbability]**

### 7.2 ContaminationFilter and NormalArtifactFilter

ContaminationFilter reads only the **first row** of each contamination table,
keyed by sample; absent sample uses default 0. Clamp c to [0,1−1e-10].
For each tumour ALT with population frequency f, count k and total AD n:
single-contaminant likelihood = 2f(1−f)Binomial(k;n,c/2)+f²Binomial(k;n,c);
many-contaminant likelihood = Binomial(k;n,c×f). Take **maximum**, not sum,
of these likelihoods. Compare with learned somatic count likelihood and §5.2's
somatic prior. Aggregate per-ALT sample posteriors by ALT-AD weighted median.
The contamination table's error column is not used here.
**[source: ContaminationFilter constructor/calculateErrorProbabilityForAlleles]**

NormalArtifactFilter selects max-TLOD ALT. Sum normal and tumour AD separately.
If normal ALT fraction <**.1×tumour ALT fraction**, probability 0, regardless
of NALOD. Preserve this early `<` return: zero tumour AD makes tumour AF NaN,
so the comparison is false and calculation continues. Replacing the fallthrough
with `normalAF >= .1*tumourAF` loses `normal_artifact` at chr4:1394834 in S1-r1. Otherwise take **negative** NALODln and real-vs-artifact prior.
Additionally compute p=1−BinomialCDF(normalAlt−1;normalDepth,10^(−refMBQ/10));
if p<**.001** return 1 instead. Ref MBQ is the **pooled INFO reference MBQ**,
not a separately emitted normal MBQ, despite the explanatory comment. Missing
MBQ uses 30. **[source: NormalArtifactFilter.calculateErrorProbability]**

### 7.3 PolymeraseSlippageFilter

Requires at least two RPA entries. Reference STR length=RU.length×RPA[0].
Test only reference length **≥8** and |RPA[0]−RPA[1]|==**1**. First ALT controls
eligibility even on multiallelics. Sum tumour AD; somatic likelihood uses total
ALT count, but artifact likelihood is
`Beta.regularizedBeta(.1, AD[1]+1, AD[0]+1)`.
On MaxCountExceededException substitute Binomial(AD[1];depth,.1).
Compare log somatic vs log artifact with the first ALT's somatic prior.
Other contexts return 0. STRQ is written only where required RPA,RU exist.
**[source: PolymeraseSlippageFilter.calculateErrorProbability/requiredInfoAnnotations]**

## 8. Haplotype propagation, error composition and decisions

FilteredHaplotypeFilter uses the **previous learning pass** map keyed by direct
concatenation `PGT+PID` (no separator, no contig/sample field). Choose tumour
genotype with highest maximum FORMAT AF; absent phase key gives 0. Maximum
stored probability within **≤100** bases is returned, including the current
position if stored. Not limited to strictly neighboring records.
**[source: FilteredHaplotypeFilter.calculateErrorProbability/makePhasingString]**

During accumulation, obtain max germline probability. If >**.25**, ignore
normal-artifact contributions. Take maximum of all other probabilities except
NON_SOMATIC types and this haplotype filter itself, and store for each phased
tumour genotype. The code comment says “non-sequencing” but the predicate
excludes only NON_SOMATIC: **sequencing error is included**. `learnParameters`
assigns accumulating map to learned map; clearing replaces the accumulating map
without clearing the learned one. **[source: accumulateDataForLearning,
learnParameters/clearAccumulatedData]**

### 8.1 Error probabilities used for training/thresholding

`ErrorProbabilities` calculates filters in list order, preserving a LinkedHashMap;
removes empty vectors and symbolic ALT entries; groups by ErrorType. For each
allele take **maximum within each type**, then independent union across types:
pCombined=`1−product_type(1−pMaxType)`, clipped to [0,1]. The type-to-vector map
is a HashMap of enum keys; its traversal order is a floating-point concern,
not permission to reorder per-filter AS statuses. It is not the product of
all nineteen individual filters. **[source: ErrorProbabilities constructor]**

Training ignores records whose only ALT is symbolic NON_REF. All regular
records contribute probabilities to the threshold learner, while clustering
may omit obvious artifacts/non-somatic alleles as §5 specifies. No original
FILTER state is used to exempt a record from recalculation.
**[source: Mutect2FilteringEngine.accumulateData]**

### 8.2 ThresholdCalculator

Default OPTIMAL_F_SCORE sorts all allele pCombined ascending. Expected true
variants=sum(1−p) using stream sum. Start TP=FP=0, FN=expectedTrue, best F=0,
best index=−1. Include each allele in sequence: TP+=1−p, FP+=p, FN−=1−p.
`F=(1+beta²)TP/((1+beta²)TP+beta²FN+FP)`, beta=1. Update best on **≥**,
so later equal maxima win. Threshold=0 when none, 1 when best is last,
otherwise that best allele's probability. **[source: ThresholdCalculator.
calculateThresholdBasedOnOptimalFScore]**

Alternative FALSE_DISCOVERY_RATE (not selected here) returns the previous
posterior at the first prefix whose expected FP fraction exceeds requested
.05, 0 if the first already exceeds it, or 1 if no prefix does.
CONSTANT leaves initial .1 unchanged. Threshold learning clears collected
probabilities each time. The .05 parameter is **not** the default strategy's
objective. **[source: relearnThresholdAndClearAcumulatedProbabilities;
calculateThresholdBasedOnFalseDiscoveryRate]**

### 8.3 Final VCF decisions and allele statuses

`applyFiltersAndAccumulateOutputStats` begins by clearing input FILTER set.
Stats record combined probabilities using **threshold−1e-10** (§9).
For decisions clamp threshold to [1e-10,1−1e-10]. Individual filter applies
when its probability **>** that threshold, not ≥. The combined probability
chooses the threshold but is **not itself a named output filter**.
**[source: Mutect2FilteringEngine.applyFiltersAndAccumulateOutputStats]**

For each allele-specific filter, produce name or `SITE` per ALT. Transpose,
distinct filter names preserving construction order, remove SITE if any name
exists, or use SITE alone. Join names by comma and ALT groups by `|` into
**AS_FilterStatus**. Symbolic alts use SITE. Names here are not alphabetical.
If every ALT has the same non-SITE outcome for one allele filter, add that name
as site filter with probability 1. **[source: addFilterStrings,
getDistinctFiltersForAllele/applyFiltersAndAccumulateOutputStats;
annotator/AnnotationUtils.encodeStringList/encodeAnyASListWithRawDelim]**

Variant-filter probabilities apply site-wide. Where required annotations exist,
variant filters that expose a quality annotation add `QualityUtils.errorProbToQual(p)`.
If all alts fail but for different allele-specific reasons and no site filter
was added, add FAIL. If any nonsymbolic ALT's allele status includes SITE,
that fallback is not applied. Site filter reporting keeps only entries with
p≥min(maxSiteErrorProbability,**.1**). This suppresses minor reported filters
without changing the site's failed status. htsjdk sorts resulting site FILTER
names alphabetically; an empty applied set writes PASS.
**[source: applyFiltersAndAccumulateOutputStats;
htsjdk variant/vcf/VCFEncoder.getFilterString]**

Quality-annotation asymmetry: the loop writing qualities handles
**Mutect2VariantFilter only**. Thus GermlineFilter writes GERMQ,
ReadOrientationFilter ROQ, PolymeraseSlippageFilter STRQ; the allele-filter
quality hooks SEQQ, CONTQ, STRANDQ are never written by this loop.
QualityUtils.errorProbToQual evaluates `Math.round(−10×Math.log10(p))`,
converts through double to int, and bounds the byte-quality to **[1,93]**;
zero probability saturates at 93 and p=1 yields **1**, not 0. GERMQ/ROQ/STRQ are serialized as
integers even though ROQ's header Type is Float.
**[source: applyFiltersAndAccumulateOutputStats;
utils/QualityUtils.errorProbToQual/boundQual; golden field values]**

## 9. Filtering output, headers and statistics

FilterMutectCalls preserves record order and emits **every** input record,
including failed variants. onTraversalStart removes the unfiltered warning
`filtering_status` and adds filtered status text; adds STRQ, the complete
GATKVCFConstants.MUTECT_FILTER_NAMES set, MUTECT_AS_FILTER_NAMES INFO definitions,
and the FilterMutectCalls command line. Existing caller declarations and sample
order persist. Header definition AS_FilterStatus still describes ApplyVQSR;
do not rewrite that legacy wording. **[source: FilterMutectCalls.onTraversalStart/
nthPassApply; utils/variant/GATKVCFConstants and GATKVCFHeaderLines; golden]**

Filtered golden header FILTER IDs, in serialized order:
FAIL, PASS, base_qual, clustered_events, contamination, duplicate, fragment,
germline, haplotype, low_allele_frac, map_qual, multiallelic, n_ratio,
normal_artifact, orientation, panel_of_normals, position, possible_numt,
slippage, strand_bias, strict_strand, weak_evidence. possible_numt is declared
without a filter instance in this engine. There are 6,351 records, 2,024 PASS;
raw/filtered CHROM,POS,REF,ALT and all genotype columns agree. Qualities remain `.`.
**[golden; source: header-constant sets and FilterMutectCalls]**

### 9.1 FilteringOutputStats

Counts are **per nonsymbolic ALT**, not per VCF record. For combined p greater
than t=threshold−1e-10 add FN+=1−p; otherwise increment pass and add FP+=p,
TP+=1−p. Per filter, if its allele probability>1e-10 and >t−1e-10, add
1−pCombined to that filter's FN; otherwise, if the allele passes combined
threshold, add its individual probability to that filter's FP. These expected
counts are not the number of records whose FILTER contains the name.
**[source: FilteringOutputStats.recordCall;
Mutect2FilteringEngine.applyFiltersAndAccumulateOutputStats]**

Keep rows only when unrounded expected FP>0 or FN>0. Row order is the filter
construction order. Columns `filter\tFP\tFDR\tFN\tFNR`; FDR=filterFP/pass,
FNR=filterFN/(TP+FN). DataLine.set(value,2) calls roundToNDecimalPlaces then
Double.toString, so `0.0`, `21.1` and `5.88` have different visible widths;
there is **no fixed two-decimal padding**. Rounding formula is
`Math.round((x+Math.ulp(x))*pow(10,n))/pow(10,n)`.
**[source: FilteringOutputStats.writeFilteringStats;
FilterStats.Mutect2FilterStatsWriter.composeLine;
utils/tsv/DataLine.set(double,int); utils/MathUtils.roundToNDecimalPlaces]**

Metadata rows precede header, each `#<METADATA>key=value`, in this order:
ln priors for signed lengths −10..10 (deletion length 10 through 1, SNV,
insertion 1 through 10); background and high-AF cluster parameters; discrete
clusters in descending weight; threshold, fdr, sensitivity. Prior values use
Double.toString; beta cluster weight `%.4f`, alpha/beta `%.2f`; fuzzy binomial
mean `%.3f`. Threshold/FDR/sensitivity use the same ulp-adjusted rounding to
**three** places then Double.toString. **[source: SomaticClusteringModel.clusteringMetadata;
BetaBinomialCluster.toString/BinomialCluster.toString; FilterStats.writeM2FilterSummary]**

Golden metadata: threshold **0.511**, FDR **0.026**, sensitivity **0.963**;
these rounded values do not reveal the exact internal decision threshold.
Background weight .8032, alpha 1.79, beta 3.17; high-AF weight .0027,
alpha 10.02, beta .50; one fuzzy binomial weight .1941, mean .990.
Prior of SNV `-11.08973727984901`; insertion length 5 is the 1e-9 floor
`-20.72326583694641`. **[golden: filteringStats.tsv]**

The eight emitted rows (FP,FDR,FN,FNR) are:

| Filter | FP | FDR | FN | FNR |
|---|---|---|---|---|
| weak_evidence | 5.88 | 0.0 | 19.11 | 0.01 |
| strand_bias | 11.22 | 0.01 | 3.91 | 0.0 |
| contamination | 21.1 | 0.01 | 28.22 | 0.01 |
| normal_artifact | 0.19 | 0.0 | 0.46 | 0.0 |
| orientation | 9.5 | 0.0 | 9.09 | 0.0 |
| slippage | 0.0 | 0.0 | 0.0 | 0.0 |
| haplotype | 2.77 | 0.0 | 2.04 | 0.0 |
| germline | 7.14 | 0.0 | 7.2 | 0.0 |

The all-zero displayed slippage row survives because row selection occurs
**before rounding**. File is 1,900 bytes; SHA256
`f71cc635a02d6533bb30794f41e34132e6b3d8983a6fa61518b46fae4e0de025`.
**[golden; source: FilteringOutputStats.writeFilteringStats]**

### 9.2 Numeric dependency contract inspected beyond GATK

Commons Math references in this subsection use tag **MATH_3_5**, commit
`b3c5dae8f253fcb4484e5cd3cc5662587803efc2`, paths under
`src/main/java/org/apache/commons/math3/`. The tree is outside the repository;
these are rules for a fresh implementation, not imported library code.

`distribution/BinomialDistribution.logProbability` uses
SaddlePointExpansion.logBinomialProbability, not a direct logGamma combination;
n=0 returns log probability 0 only at k=0, otherwise −Infinity. `probability`
returns zero for −Infinity or **FastMath.exp** otherwise. Its CDF returns 0
below k=0, 1 at k≥n, otherwise
`1−Beta.regularizedBeta(p,k+1,n−k)`. Thus NormalArtifactFilter subtracts an
already-subtracted incomplete beta from one; replacing the whole expression
with the original beta value can change cancellation. **[source: those methods]**

`special/Beta.regularizedBeta` defaults epsilon **1e-14** and maximum iterations
**Integer.MAX_VALUE**. It uses symmetry when
`x>(a+1)/(2+b+a) AND 1−x≤(b+1)/(2+b+a)`, then the continued-fraction path.
The source-specific Beta.logBeta routine is also what beta-binomial calls;
logGamma subtraction is not assumed equivalent. **[source: Beta.java,
regularizedBeta overloads and logBeta]**

`special/Gamma.logGamma` branches at x=.5,2.5,8 and uses FastMath.log; small x
uses logGamma1p, larger x uses finite products or Lanczos with that source's
coefficients. `digamma` for positive x≤1e-5 returns −EulerGamma−1/x;
x≥49 uses FastMath.log and its asymptotic polynomial, otherwise recurse on x+1
and subtract 1/x on return. Mutect activity's integer values are cached by
GATK `utils/DigammaCache` (max cache size 100,000). Somatic posterior digammas
have noninteger arguments and cannot be replaced with that integer cache.
**[source: Gamma.logGamma/digamma; GATK DigammaCache.compute/maxSize]**

`optim/univariate/BrentOptimizer.doOptimize` uses
`tol1=relativeTolerance×abs(x)+absoluteTolerance`, tol2=2tol1; stops when
`abs(x−(a+b)/2)≤tol2−(b−a)/2`. It tries an interior parabolic step only within
its strict bound/step checks, otherwise golden section constant
`.5×(3−FastMath.sqrt(5))`; moves at least tol1. Maximization negates objective
values during search and returns the **best encountered** point, not merely
last x. Internal `best` keeps the older point on equal objective. No extra
convergence checker is supplied by OptimizationUtils.max. Interpolation and exhausted-budget behavior have the independent numerical
oracles in §3.6; §9.4 additionally checks two full strand M steps. **[source: BrentOptimizer.
doOptimize/best; GATK utils/OptimizationUtils.max]**

The TSV decimal-place rounding in §9.1 is floating arithmetic plus one ulp,
not BigDecimal decimal rounding. VCF Double formatting uses Java Formatter's
specified precision, while histogram counts use htsjdk FormatUtil HALF_DOWN.
These are distinct contracts; replacing all three with one “round to n” helper
would lose parity. **[source: GATK MathUtils.roundToNDecimalPlaces;
htsjdk VCFEncoder.formatVCFDouble and FormatUtil constructor]**

### 9.3 FilterMutectCalls fresh port and S1-r1 verification

`nomics somatic filter-mutect-calls` implements the four traversals in
`crates/nomics-somatic/src/filtering{.rs,/}`. It accepts `--variant`, `--output`,
`--reference`, `--orientation-bias-artifact-priors`, `--tumor-segmentation`,
`--contamination-table` and optional `--stats` (otherwise `<variant>.stats`).
The latter is explicit in the golden test because the archived stats live in
MERGEMUTECTSTATS rather than beside the MERGE_MUTECT2 VCF. No reference executable
is called or linked. No numerical, bioinformatics or compression dependency was
added; unsafe remains denied. VCF fields retain their input textual values,
new INFO is merged in sorted-key order, and sample/FORMAT columns are preserved.
Contig declarations retain their dictionary order (htsjdk contig comparison
uses the contig index, not lexical ID sorting). New header declarations were
checked against `GATKVCFHeaderLines`, including the legacy ApplyVQSR description.

The port reuses the orientation model's E step and checked tar reader, its
beta-binomial/logBeta/FastMath arithmetic, the earlier Brent, binomial PMF,
Java log/exp, compensated stream sum, Java round and Double.toString ports.
General positive-shape Gamma/logGamma/digamma and the regularized beta continued
fraction extend the numerical boundary. `tests/oracle/FilteringNumerics.java`
independently calls unmodified Commons Math **MATH_3_5** (commit
`b3c5dae8f253fcb4484e5cd3cc5662587803efc2`). Its **70 binary64 fixtures** cover
Gamma branches, noninteger digammas, noninteger/small/large logBeta and CDF
symmetry. All 70 match exactly; existing orientation/contamination numerical
and golden tests still pass. Regenerate with an installed JDK 17:

```sh
mkdir -p /tmp/nomics-filter-oracle
/opt/homebrew/opt/openjdk@17/bin/javac \
  -sourcepath /tmp/mutect2-spec-math3/src/main/java \
  -d /tmp/nomics-filter-oracle \
  crates/nomics-somatic/tests/oracle/FilteringNumerics.java
/opt/homebrew/opt/openjdk@17/bin/java -cp /tmp/nomics-filter-oracle \
  FilteringNumerics > crates/nomics-somatic/tests/oracle/filtering_numerics.bits
```

GATK tag 4.6.2.0 was fetched outside the repository to
`/tmp/nomics-filter-gatk`, commit `76edc75c26504da94bbaee66584e107e76ee15de`.
No third-party source is committed. Source/spec discrepancies corrected here:
orientation **max rather than sum** (§6.2) and multiplication order of cluster
shape gradients (§5.1). Additional source-confirmed clarifications are persistent
out-of-range priors (§5), the NaN-sensitive normal-artifact early return (§7.2),
and contig header order above. No new output exemption was introduced.

The integration test `tests/filtering_goldens.rs` runs the actual port, compares
stats bytes first, every decompressed record next, then all header lines except
`##GATKCommandLine`. It skips with an explicit message if
`NOMICS_SOMATIC_GOLDENS` (default `~/omics/runs/somatic/S1-r1/golden`) is absent.
`NOMICS_SOMATIC_FASTA` optionally relocates the reference. A separate unit gate
runs all six error-type multiplication orders through the entire traversal;
all six produce identical boundary files. A unit test checks long-indel prior
persistence across learning.

Exact release CLI, literal `cmp`, record-by-record comparison, PASS count,
normalized header comparison and BGZF member validation are reproducible with:

```sh
cargo build --release -p nomics
cargo test --release -p nomics-somatic -p nomics-io
cargo clippy --release --no-deps -p nomics-somatic -p nomics --all-targets -- -D warnings
sh crates/nomics-somatic/tests/oracle/filtering_cli_gate.sh /tmp/nomics-filter-verified.vcf.gz
```

The script contains the full CLI arguments and invokes:
`cmp "$E.filteringStats.tsv" "$O.filteringStats.tsv"`, followed by a Python
comparison of every complete VCF record and the normalized header. The only
header normalization removes lines beginning `##GATKCommandLine=`, exactly as
in the existing HaplotypeCaller/joint-genotyping rows. **Compressed bytes and
`.tbi` are outside this task's gate.** The writer emits valid BGZF with its EOF
member and **does not write a tabix index**; it does not copy the reference index.

### 9.4 Perturbation exactness repair (2026-11-12)

Five independent-input reference outputs under `~/omics/runs/somatic/fmc-perturb`
extend the S1-r1 gate. Before this repair, chr1to7 had ten differing prior lines
and rest had two, despite identical records, rounded clusters and per-filter
statistics. The reference's four repeated runs of each were deterministic
(as reported by the lead). No tolerance or new normalization is introduced.

Cause: `StrandArtifactFilter.strandArtifactProbability` calls
`MathUtils.normalizeLog10`, whose log-sum and conversion use Java Math.pow and
Math.log10. The port used platform pow/log10. A snapshot trace of chr1to7 found
2,116 differing forward/reverse responsibilities across four traversals with
identical supplied strand parameters. The first was nF=1,nR=37,kF=0,kR=15,
length=0, alpha=1,beta=20,prior=.001: reverse responsibility Java
`6.51288763118399E-6`, platform `6.5128876311839904E-6`. These are learning
inputs, so their changes propagate through strand learning and the retained
Datum artifact probabilities into the somatic EM priors. Snapshot checks of
somatic sequencing posteriors, somatic weights, count likelihoods and subsequent
cluster-shape/weight updates against Java agreed exactly when supplied the same
state and data; changing summation order was not the remedy.

The repair reuses nomics-hc's existing fresh OpenJDK log10/log1p/expm1/Pow ports
inside nomics-somatic, preserving the existing log/exp implementations. All
filtering scalar call sites now use those routines: strand E-step pow/log10,
initial and split-cluster log1p, posterior log1mexp's log1p/expm1, POPAF inversion,
normal-artifact quality conversion, output quality log10 and TSV decimal scale
Math.pow. The separate source-expression correction is −Math.log(3), §5.
No reduction order, persistent-prior rule or optimizer budget was changed.
Review of cluster initialization/EM, strand M steps, threshold scan and output
statistic accumulation found their reduction orders already matched the source;
§5.1 records those distinctions. The count-map comment was corrected from
HashMap.merge to Collectors.toMap's JDK 17 putIfAbsent; execution is unchanged.

`tests/oracle/FilteringLearning.java` is an independent Java harness using
Commons Math 3.5 directly and the GATK expressions. Its **83 binary64 fixtures**
cover the SNV constants, scalar functions, both product orders below, seven
strand E steps (SNV and indels 1–4), two complete strand M steps with Commons
Math Brent, and the next E steps after each update. No runtime reference tool
or library is called by nomics. Regenerate outside the repository:

```sh
mkdir -p /tmp/nomics-filter-learning
/opt/homebrew/opt/openjdk@17/bin/javac \
  -sourcepath /tmp/mutect2-spec-math3/src/main/java \
  -d /tmp/nomics-filter-learning \
  crates/nomics-somatic/tests/oracle/FilteringLearning.java
/opt/homebrew/opt/openjdk@17/bin/java -cp /tmp/nomics-filter-learning \
  FilteringLearning > crates/nomics-somatic/tests/oracle/filtering_learning.bits
```

The existing golden integration test remains. `filtering_perturbed_inputs` adds
all five cases, stats bytes first, every complete record, then the normalized
header. `NOMICS_SOMATIC_FMC_PERTURB` overrides the default perturbation directory;
an absent directory prints an explicit SKIP. `NOMICS_SOMATIC_FASTA` relocates the
reference for both tests. No fixture output is written beside the inputs by
the integration test.

Final verification on this Mac, 2026-11-12 (all commands exit 0):

```text
$ cargo build --release -p nomics
Finished `release` profile [optimized + debuginfo] target(s) in 26.43s
$ cargo test --release -p nomics-somatic -p nomics-io
121 passed; 0 failed; 1 pre-existing ignored test (across all test binaries)
$ cargo clippy --release --no-deps -p nomics-somatic -p nomics --all-targets -- -D warnings
Finished `release` profile [optimized + debuginfo] target(s)
```

Build/Clippy retain the unrelated nomics-hc deprecated fetch_update warning;
strict targeted Clippy now passes after merging main's two contamination-test
lint fixes, without lint exemptions. Regenerating the new Java oracle into a
fresh /tmp directory and comparing its 83 lines also returns cmp exit 0.
The absent-perturbation integration check prints
`SKIP FilterMutectCalls perturbations: /tmp/nomics-filter-priors-absent absent`.

```text
$ bash ~/omics/runs/somatic/fmc-perturb/check.sh target/release/nomics ~/omics/runs/somatic/fmc-perturb
odd: records 3176/3176 identical; stats identical; header identical; FilterMutectCalls: 3176 records, 956 PASS
chr1to7: records 2666/2666 identical; stats identical; header identical; FilterMutectCalls: 2666 records, 857 PASS
rest: records 3685/3685 identical; stats identical; header identical; FilterMutectCalls: 3685 records, 1105 PASS
contam10: records 6351/6351 identical; stats identical; header identical; FilterMutectCalls: 6351 records, 1788 PASS
same: records 6351/6351 identical; stats identical; header identical; FilterMutectCalls: 6351 records, 2024 PASS
$ sh crates/nomics-somatic/tests/oracle/filtering_cli_gate.sh /tmp/nomics-filter-recheck.vcf.gz
FilterMutectCalls pass 0: threshold 0.5062301244631351
FilterMutectCalls pass 1: threshold 0.5147533777536284
FilterMutectCalls pass 2: threshold 0.5106759443380069
FilterMutectCalls: 6351 records, 2024 PASS
real 0.67
user 0.66
sys 0.01
stats cmp: PASS (1900 bytes)
records: 6351/6351 identical; PASS: 2024
header: PASS (only ##GATKCommandLine excluded)
BGZF: PASS (41 members); tabix index: not emitted
```

All six normalized headers contain 3,435 identical lines. Stats lengths for
odd/chr1to7/rest/contam10/same are 1,951/1,891/1,893/1,955/1,900 bytes.
GATK was not rerun; the stored reference outputs are the boundary oracle.
GATK tag 4.6.2.0 was freshly fetched to `/tmp/nomics-gatk-4.6.2.0`, at the same
commit recorded in §9.3. No third-party source is committed. The gate rows
are adopted in `docs/GATES.md`.


### 9.5 ErrorProbabilities enum-map order can reach output

`ErrorProbabilities` collects per-type maxima into a HashMap keyed by ErrorType,
then traverses its values and reduces `(1−p)` with multiplication, starting at
1.0. Enum identity hashes can change the order. Multiplication of binary64
numbers is commutative here but **not associative**; the final `1−product` and
clamp do not generally erase the difference. The combined probabilities feed
ThresholdCalculator and FilteringOutputStats. There is no universal arithmetic
proof that written bytes are independent of this order.

A concrete arithmetic-boundary witness, checked with Java Math.round and in the
Rust threshold/rounding path: probabilities
`[0.04163411606541967,0.006318578253251994,0.03865375732824039]` produce
`0.08450000000000002` for order [0,1,2], but `0.08449999999999991` for [0,2,1].
Put either combined value between 0 and 1 in the threshold's sorted input:
OPTIMAL_F_SCORE selects that interior value, and the **actual ulp-adjusted
three-place writer** yields threshold text `0.085` versus `0.084`. Thus low bits
can survive all the way to a written metadata byte; clamping and rounding are
not sufficient safeguards. This proves the arithmetic pipeline's sensitivity,
**not** the existence of a VCF realizing those particular per-type probabilities.
No such end-to-end GATK witness is claimed. The six permutations still reproduce
the S1-r1 golden boundaries; the five new default-order fixture gates do not
prove independence for all other inputs/JVMs. No enum order is masked or
silently normalized, and the port retains [SEQUENCING,ARTIFACT,NON_SOMATIC].

## 10. Not confirmed

* Orientation learning is now gated from **both original scatter f1r2 archives**:
  ordered member name and all prior bytes match S1-r1, with no normalization.
  Independent Java numerical/EM fixtures settle default arithmetic and reduction
  order; the golden tests also reverse the inputs. Multi-sample archive ordering,
  extended tar layouts and custom EM/max-depth settings are explicitly unsupported
  rather than claimed equivalent. A full instrumented GATK trajectory for every
  golden context and full FilterMutectCalls execution remain unperformed.
* GetPileupSummaries BAM-derived counts now pass both complete chr21 byte gates
  and all ordered filter counts (§2.4). Source and hand-built tests settle
  supplementary-read rejection and zero/deletion-only coverage. The full
  13–14 GB tumour/normal BAM gates have not been run on this Mac; no claim of
  a full-run recomputation is made. Requested default semantics are confirmed
  in the source; synthetic edge cases have not been rerun through a JVM.
* Contamination and segmentation now pass the fresh-port four-table gate and
  an independent Java-source oracle (§3.6). This confirms filtered-site counts
  and final outputs, with numerical unit oracles; exhaustive intermediate SVD
  factors, optimizer trajectories and fallback strategies still lack stage goldens.
* Commons Math **3.5** Gamma/Beta and FastMath overloads outside the contamination
  PMF/segmentation path remain separate numerical dependencies. SVD, full Brent
  interpolation, integer-input SaddlePointExpansion, and the called FastMath
  exp/log/hypot/scalb paths are now source-ported with Java bit oracles (§3.6).
  The integer-shape Beta.logBeta and FastMath.log/log1p paths needed by the
  orientation model are ported separately and independently tested (§1.5).
  No BigDecimal-based rounding was found in the inspected TSV/VCF writer path;
  unrelated paths were not exhaustively searched.
* FilterMutectCalls now passes the S1-r1 and five perturbation gates (§9.4),
  including statistics,
  all records and the normalized header. GATK itself was not rerun or instrumented
  during this port; the stored reference is the oracle. Per-pass internal cluster
  trajectories, strand parameters and all unrounded probabilities have no full
  instrumented JVM comparison. Final threshold on the port is
  `0.5106759443380069`; the reference file only publishes its rounded `0.511`.
* Enum-key HashMap multiplication order in ErrorProbabilities depends on JVM
  enum identity hashes. All **six permutations**, each through the complete
  four-pass port, reproduce the exact S1-r1 VCF and stats boundaries. The port
  uses SEQUENCING, ARTIFACT, NON_SOMATIC; this is not a claim to recover the
  original JVM's iteration order or identical internal low bits. §9.5 gives a Java-checked arithmetic witness in which the order changes
  the written threshold digits. A realizable full-VCF witness and the exact
  enum iteration order of any particular reference JVM remain unconfirmed.
* Filtering supports only the default two-sample chain with one tumour, one
  normal, one prior archive and the supplied contamination/segmentation tables.
  Tumour-only, mitochondrial/microbial modes, multiple tumour/normal samples,
  symbolic/gVCF alleles, missing callable statistics, no-callable-site fallback,
  overlapping tumour segments, and custom filtering parameters are unsupported.
  The CLI does not expose these modes; sample count, symbolic alleles, missing
  stats and visited overlapping segments return errors. Priors reuse the existing
  checked regular-ustar transport; extended/nested tar entries are rejected.
  Generic htsjdk re-encoding for arbitrary VCF inputs, malformed annotations,
  integer overflow, Integer HashMap tree bins, Beta nonconvergence/fallback, and
  exact Java Formatter midpoint behavior outside this fixture are not validated.
  Filtering Math log1p/expm1/pow/log10 call sites now use the fresh OpenJDK
  scalar ports with independent Java fixtures (§9.4). Exhaustive agreement with
  every JVM intrinsic/platform remains outside the finite fixture gate.
* Alternate parameter regimes, empty/all-error inputs, overlapping segments,
  malformed resource AFs, multiple tumours/normals and symbolic/gVCF alleles
  need synthetic fixtures. This target is the default two-sample VCF chain;
  source branches above are not claims of validation for other modes.
* The orientation tar/gzip envelope exclusion is now explicit in `docs/GATES.md`
  and §1.5, with exact names/member bytes still required. Other caller tar
  metadata and compressed VCF/.tbi bytes still require the explicit normalization
  and two-sided writer gates proposed in `mutect2.md` §9.
