# FastQC 0.12.1 text reports

`nomics qc fastqc --fastq INPUT --outdir OUTPUT [--input-name STAGED_NAME]`
writes `OUTPUT/fastqc_data.txt` and `OUTPUT/summary.txt` for **one raw mate**.
The gate is byte identity to these members of the reference `*_fastqc.zip`,
including module status, order, tabs, number formatting and final newline.
There is no normalization. The default input name is the supplied path's basename;
`--input-name` supplies sarek's symlink basename without renaming the input.
For tier0, invoke this twice with `HG002_WES-L1_1.gz` and
`HG002_WES-L1_2.gz`. These are the inputs **before fastp**. Sarek's confirmed command:

```sh
fastqc --quiet --threads 4 --memory 1024 HG002_WES-L1_1.gz HG002_WES-L1_2.gz
```

The implementation is fresh Rust, without new dependencies, unsafe, FFI,
bioinformatics/compression libraries or production invocation of Java/FastQC.
Gzip (including multiple members) goes through `nomics_io::gzip::GzReader`,
with its existing CRC/ISIZE checks. Plain input is also accepted (gzip magic
is detected independently of the extension). The four reusable line buffers
follow nomics-prep's private `Fastq` reader; that parser cannot be imported
without an out-of-scope change to nomics-prep. No transient per-record allocation
occurs in the hot loop: line/position buffers are reused, and duplication keys are stored inline in the bounded table (at most
100,000 keys), with no per-key heap object. Allocation is limited to aggregate
buffer/table growth on new sizes, tiles or capacity thresholds.

## Source ledger

Source fetched outside the worktree at tag **v0.12.1**, commit
`e7ef390bf10382f60786bdd0cf28abd4f8683ffd`, from
[s-andrews/FastQC](https://github.com/s-andrews/FastQC/tree/v0.12.1).
Paths below are relative to `uk/ac/babraham/FastQC/` in that tree.

| Behavior | Source |
|---|---|
| Standard module order | `Modules/ModuleFactory.java` |
| FASTQ four-line parsing, blank lines, uppercase, colorspace; CASAVA filtering requires a nondefault option | `Sequence/FastQFile.java`, `Sequence/Sequence.java`, `Analysis/AnalysisRunner.java` |
| Counts, base-call denominator, integer %GC, truncated displayed total bases | `Modules/BasicStats.java` |
| Encoding detection from minimum observed character | `Sequence/QualityEncoding/PhredEncoding.java` |
| Per-position histogram, mean, integer percentile rank | `Utilities/QualityCount.java` |
| Position-mean/percentile grouping and status | `Modules/PerBaseQualityScores.java` |
| Tile-name parsing, sampling, cross-tile normalization | `Modules/PerTileQualityScores.java` |
| Integer mean quality per read, lowest modal tie | `Modules/PerSequenceQualityScores.java` |
| Called-base and N denominators | `Modules/PerBaseSequenceContent.java`, `Modules/NContent.java` |
| Linear long-read bins | `Graphs/BaseGroup.java` |
| GC truncation, fractional bins, theoretical distribution | `Modules/PerSequenceGCContent.java`, `Modules/GCModel/{GCModel,GCModelValue}.java`, `Statistics/NormalDistribution.java` |
| Length bins and removal of empty endpoint padding | `Modules/SequenceLengthDistribution.java` |
| Distinct-key admission, stable count sort | `Modules/OverRepresentedSeqs.java` |
| Duplication extrapolation, category boundaries | `Modules/DuplicationLevel.java` |
| Contaminant matching and tie-breaks | `Sequence/Contaminant/{Contaminant,ContaminantHit,ContaminentFinder}.java` |
| Cumulative first-hit adapter counts, expansion, grouped means | `Modules/AdapterContent.java` |
| Module framing, summary spelling, LF | `Report/HTMLReportArchive.java` |
| Limits, enabled modules, adapter/contaminant data | `Configuration/{limits,adapter_list,contaminant_list}.txt` |

The adapter and contaminant values are embedded in Rust with their source-file
names. All status constants below come from `Configuration/limits.txt`, rather
than comments or the website. Kmer is disabled by that file. Module order is
Basic Statistics, Per base sequence quality, Per tile sequence quality (when
available), Per sequence quality scores, Per base sequence content, Per sequence
GC content, Per base N content, Sequence Length Distribution, Sequence
Duplication Levels, Overrepresented sequences, Adapter Content.

## Numeric and collection behavior

FastQC's string append delegates doubles to Java. `fastqc/float.rs` implements
the legacy OpenJDK 17 `jdk.internal.math.FloatingDecimal` conversion decisions:
integer shortcut, exponent estimate (including its truncated constant), symmetric
half-ULP stopping bounds, special treatment of exact powers of two, signed
32/64-bit intermediate overflow, arbitrary-precision path and ties-to-even.
Formatting retains `.0`, `NaN`, signed zero, and uppercase `E` without `+` or
exponent zero-padding. Scientific notation is used below 0.001 and at 10^7.
This is deliberately different from simply applying Rust's shortest formatter:
for example 10^23 becomes `9.999999999999999E22` on the legacy JVM.
The JDK source was read from the locally installed OpenJDK 17 `src.zip`, outside
the worktree. A Java oracle matched 20,010 selected/random bit patterns.
The original golden JVM build is not recorded here; the four golden text files
confirm the selected behavior on their actual values.

`fastqc/java_map.rs` replays the OpenJDK 17 `java.util.HashMap` insertion/resize
and red-black tree-bin ordering for comparable keys. Java String hashes are signed
32-bit polynomial hashes with the normal spread. This preserves stable-sort
count ties in overrepresented output and integer-key iteration order in both
maps used by duplication estimation. Linked-list insertion order alone is
insufficient when a bin treeifies. Standalone Java oracle checks covered 8, 16,
64 and 256 deliberately colliding strings and 100,000 random DNA keys.

## Source details that are easy to miss

- Basic `%GC` excludes N/other characters and uses integer division. Total Bases
  chooses bp/kbp/Mbp/Gbp and **truncates** to at most one fractional character.
  Encoding is Sanger/Illumina 1.9 below ASCII 64; exactly 65 selects Illumina 1.3;
  other 64..126 selects Illumina 1.5. Flags in IDs are not filtered without
  `--casava`, which sarek does not pass.
- Position percentiles require **more than 100** observations. Their rank is
  `floor(total * percentile / 100)`, then the first cumulative bin meeting that
  rank. A grouped percentile averages per-position percentiles, rather than
  pooling histograms. Mean quality also averages position means. Insufficient
  percentile data prints `NaN` and is excluded from quality status decisions.
- Linear base bins are single positions through length 75. Above that, positions
  1..9 remain individual; the first of 2/5/10 times a power of ten giving **fewer
  than 75** groups is used. When interval >10, the next bin is 10..interval-1.
- Per-read quality averages raw characters with integer division before
  subtracting the offset. Modal quality ties choose the lowest score.
- Tile quality samples the first 10,000 nonempty quality strings and then every
  tenth. The first sampled name selects colon field 4 (7+ fields), otherwise 2
  (5+ fields). An unparseable sampled name disables the entire module; more than
  2500 existing tiles also disables it on the next new tile. Means are normalized
  against the unweighted average of tile means, in sorted tile order. Status uses
  the largest **absolute** deviation, including better-than-average tiles.
- Base composition excludes N/other characters; N percentages include every
  character in the denominator. Grouping pools counts before division.
- GC counts only a prefix rounded down to a multiple of 100 for lengths >100,
  or 1000 for lengths >1000. The GC model assigns a count to every rounded
  percentage claimed by its +/-0.5-base interval, weighted by the inverse number
  of claims for that percentage. Additions retain read order. The theoretical
  Gaussian uses the first mode and its contiguous neighbors above **90%**, even
  though the source comment says 95%. If either side reaches a boundary, it uses
  the first mode alone. Variance divides by total weight minus one; deviation is
  absolute observed/theoretical distance divided by total weight. Theoretical
  points affect status only and are not written to the text report.
- Length bins choose 1/2/5 times a power of ten, targeting integer span/interval
  <=50. One-position padding is added on each side; only zero-count **endpoint**
  bins are omitted. Interior zeros remain.
- Duplication comments still mention 75 bp, but **v0.12.1 truncates every read
  longer than 50 bp to 50 bp**. Admission freezes at **100,000 distinct keys**,
  not 100,000 reads. Previously admitted keys continue accumulating. The last
  read count before freezing is the extrapolation limit. The hypergeometric
  non-observation probability is multiplied term by term, with the exact source
  early exits and 0.01-observation cutoff; no log approximation or closed-form
  replacement is used. Counts 10..50 occupy the slot labeled `>10`; similarly
  51..100, 101..500, 501..1000, 1001..5000, 5001..10000, >10000 occupy the
  remaining slots. Percentages are normalized to the inferred raw total.
- Overrepresented keys are retained at strictly >0.1%. Stable descending-count
  sorting retains Java map order for ties. No hits prints `No Hit`; an empty
  module has no column header. Short contaminant queries 8..19 allow exact
  substrings. Longer matching uses both orientations and offset scanning, with
  a segment longer than 20 bp and the source's one-mismatch/reset behavior.
  Within an offset length then identity wins; between offsets/orientations only
  strictly greater length wins. Across contaminants length then identity wins,
  retaining the first list entry on an exact tie.
- Adapters count the first exact substring occurrence and increment all later
  possible positions. Extending the longest read copies the previous cumulative
  tail. Each position is converted to `(count * 100.0) / reads` **before** averaging
  within its group; replacing this with one pooled fraction changes bytes. The
  longest adapter is 12; a maximum read length <=12 leaves the module warning
  with no table in the upstream implementation.

| Module | Warn / fail, with source comparison |
|---|---|
| Basic | Always pass |
| Base quality | Lower quartile <10 / <5 **or** median <25 / <20 |
| Tile quality | Absolute deviation >5 / >10 |
| Sequence quality | Modal score <=27 / <=20 |
| Base composition | A/T or C/G percentage difference >10 / >20 |
| GC distribution | Deviation >15 / >30 percent |
| N content | Any grouped percentage >5 / >20 |
| Length | Warn for multiple lengths; fail if zero-length reads exist |
| Duplication | Deduplicated percentage <70 / <50 |
| Overrepresented | Any key >0.1 / >1 percent |
| Adapter | Any grouped percentage >5 / >10; warn if analysis too short |

## Verification and limits

`tests/fastqc_goldens.rs` reads stored/deflated ZIP members directly through
nomics-io's raw inflate; it checks central-directory size and CRC and handles
FastQC's data-descriptor entries. It compares **both mates and both files**.
It passes with a skip message if either input/golden pair is absent.
`NOMICS_RUNS` defaults to `~/omics/runs`; `NOMICS_FASTQC_INPUTS` defaults to
`~/omics/data/tier0`. Output is generated afresh from FASTQ; no golden statistics
are inputs to the implementation.

Tier0: all four byte comparisons pass. Release elapsed time is recorded below
with the final verification commands. The wider gate ran on the Spark on
2026-11-12 (`scripts/qc_reads_gate.sh`): 60 samplesheet rows in 56 sarek runs,
both mates and both files, 240 of 240 identical; and the HG002 30x genome and
its slice against the x86 runs' reports, 8 of 8. Five additional 400-read FASTQs were compared with the source
compiled under OpenJDK 17: tiles with unequal quality, 1001..1003 bp reads,
9 bp reads, mixed 20/50/51/75/76/100/151/251 bp reads, and colorspace. All ten
text comparisons passed. These exercise nonempty contaminant tables and warn/fail
states missing from tier0. Unit tests retain percentile sufficiency, bin,
encoding, tile-sampling, truncation/freezing, formatter and tree-bin cases.

The CLI contract is default ASCII FASTQ (gzip or plain), with no custom limits,
adapter/contaminant files, grouping overrides, CASAVA filtering, kmer opt-in,
BAM/SAM input or bzip2 input. Read/base totals use 64-bit storage for genome-scale
inputs; per-key duplication counts retain Java's signed 32-bit wrapping. Inputs requiring billions of copies
of a single sequence have not been run end to end. Host f64 square root
and power functions implement the Gaussian status calculation; theoretical
floating-point values are not text output. No known residual difference remains
on a tested input, and no normalization has been introduced.

## Deliberately excluded artifacts

- PNG: would require reproducing Java2D graph drawing, fonts, antialiasing,
  dimensions, legends, heatmap/boxplot rendering and image encoding (including
  environment/JVM-dependent rasterization), then image-byte gating.
- SVG: would require reproducing FastQC's `Utilities/ImageSaver/SVGGenerator.java` and
  graph-specific Java2D-to-SVG output ordering, coordinate/number formatting and
  serialization. No SVG/PNG production renderer is included here.
- `fastqc_report.html`: would require `HTMLReportArchive`'s XML writer behavior,
  full template/CSS/JavaScript/icon assets, generated images/base64, navigation,
  timestamp/provenance and serializer escaping/attribute ordering, plus separate
  investigation of reference repeatability.
- `fastqc.fo`: is an XSLT transformation of the XHTML through
  `Templates/fastqc2fo.xsl`. Matching it needs the same HTML and XSLT/DOM serializer
  behavior, including embedded images and JVM implementation details.
- ZIP: would need the full entry set, names/order, Java ZIP metadata and data
  descriptors, timestamps and compression stream. The two text members alone
  do not reproduce the container. This task only decodes reference ZIPs in tests.

## Recorded local verification (2026-11-12)

```sh
cargo build --release -p nomics
cargo test -p nomics-qc
cargo test -p nomics-qc --lib
cargo clippy -p nomics-qc -p nomics -- -D warnings
cargo clippy --no-deps -p nomics-qc -p nomics -- -D warnings
NOMICS_RUNS=/tmp/fastqc-src/absent NOMICS_FASTQC_INPUTS=/tmp/fastqc-src/absent \
  cargo test -p nomics-qc --test fastqc_goldens -- --nocapture

/usr/bin/time -p target/release/nomics qc fastqc \
  --fastq ~/omics/data/tier0/tier0_1.fastq.gz \
  --input-name HG002_WES-L1_1.gz --outdir /tmp/nomics-fastqc-1
/usr/bin/time -p target/release/nomics qc fastqc \
  --fastq ~/omics/data/tier0/tier0_2.fastq.gz \
  --input-name HG002_WES-L1_2.gz --outdir /tmp/nomics-fastqc-2
cmp /tmp/nomics-fastqc-1/fastqc_data.txt /tmp/fastqc-goldens/tier0-wes-20261008-r1/HG002_WES-L1_1_fastqc/fastqc_data.txt
cmp /tmp/nomics-fastqc-1/summary.txt /tmp/fastqc-goldens/tier0-wes-20261008-r1/HG002_WES-L1_1_fastqc/summary.txt
cmp /tmp/nomics-fastqc-2/fastqc_data.txt /tmp/fastqc-goldens/tier0-wes-20261008-r1/HG002_WES-L1_2_fastqc/fastqc_data.txt
cmp /tmp/nomics-fastqc-2/summary.txt /tmp/fastqc-goldens/tier0-wes-20261008-r1/HG002_WES-L1_2_fastqc/summary.txt
```

Release build passed (the existing nomics-hc `fetch_update` deprecation warning
remains). All four literal `cmp` commands exited **0**. Final warm-cache release
elapsed time was **0.37 / 0.36 s**, 0.73 s total for the two 495,559-read inputs;
user CPU time was 0.67 / 0.67 s, including the concurrent inflate worker. This is
a local wall-time measurement, not a GPU/CPU speedup claim.
The missing-fixture golden test passed with its explicit skip message.

The strict Clippy command exited 101 on existing dependency lints:
`nomics-gpu` `too_many_arguments`, and `nomics-io` `manual_is_multiple_of`,
`type_complexity`, `collapsible_if`. No unrelated files were changed to suppress
them. The same command with `--no-deps` exited **0**, checking both requested
packages with warnings denied; the existing nomics-hc dependency deprecation
remains a dependency warning. Unit tests: **32 passed**, including six new
FastQC tests. The FastQC integration test passed all four byte comparisons.
The full `cargo test -p nomics-qc` run passed **31 unit tests and four integration
tests**, with zero failures and zero ignored tests (plus zero doc tests), including
the existing six-BAM stats, VCF and six-BAM mosdepth integration gates. After the
final boundary test and inline-key storage change, `--lib` passed all **32 unit
tests**, and the FastQC integration test passed again (8.496 / 8.712 s in the
unoptimized test build). The final release binary also passed all four literal
`cmp` commands and all ten synthetic source-oracle comparisons. Rustfmt on the
new Rust modules/test and `git diff --check` passed.
