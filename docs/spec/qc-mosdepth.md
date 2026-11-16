# mosdepth QC: source port and release gates

`nomics qc mosdepth --bam X.bam --prefix DIR/PREFIX --by targets.bed`
implements sarek's exome invocation of mosdepth 0.3.14.
`--by 500 --no-per-base --fast-mode` implements its genome invocation.
Window sizes and BED contents are inputs, not fixture-specific constants.
`--threads` defaults to four and controls BAM decompression only.
`nomics qc csi --bed-gz X.bed.gz --out X.csi [--levels 5]` rebuilds the
on-the-fly index layout independently from the BED's own virtual offsets.

The production path is fresh safe Rust. It reads BAM with `nomics-io`, uses
its existing BGZF writer and `FastDeflater`, and allocates the chromosome
array through `nomics-mem::HostArray<u32>` (HostBuf). Wrapping u32 differences
represent the source's int32 difference array without unsafe casts. The only
added dependency edge is to the existing workspace `nomics-mem`; its lockfile
entry is updated. No external dependency or `nomics-io` change is needed.
No upstream source is vendored, linked, or executed by production code.

## Source ledger

Sources were fetched outside the repository under `/tmp/nomics-*`.
Line references below are to these tags, not moving branches.

| Entry point | Source and lines | Confirmed behavior |
|---|---|---|
| Coverage | [mosdepth v0.3.14 mosdepth.nim:250](https://github.com/brentp/mosdepth/blob/v0.3.14/mosdepth.nim#L250)–360 | chromosome query, filtering, mate correction, fast path |
| CIGAR and runs | [mosdepth.nim:58](https://github.com/brentp/mosdepth/blob/v0.3.14/mosdepth.nim#L58)–96, [146](https://github.com/brentp/mosdepth/blob/v0.3.14/mosdepth.nim#L146)–172 | run boundaries, reference/query consumption |
| BED and windows | [mosdepth.nim:195](https://github.com/brentp/mosdepth/blob/v0.3.14/mosdepth.nim#L195)–209, [362](https://github.com/brentp/mosdepth/blob/v0.3.14/mosdepth.nim#L362)–413 | parsing, stable start sort, interval generation, mean accumulation |
| Distributions and summary | [mosdepth.nim:415](https://github.com/brentp/mosdepth/blob/v0.3.14/mosdepth.nim#L415)–480; [depthstat.nim:38](https://github.com/brentp/mosdepth/blob/v0.3.14/depthstat.nim#L38)–59 | histogram cap, cumulative order, sum/min/max |
| Main output loop | [mosdepth.nim:589](https://github.com/brentp/mosdepth/blob/v0.3.14/mosdepth.nim#L589)–839 | files, contig inclusion, sentinel handling, totals |
| Decimal coordinates | [int2str.nim:60](https://github.com/brentp/mosdepth/blob/v0.3.14/int2str.nim#L60)–88 | ordinary nonnegative decimal integers, no padding |
| Index setup | [mosdepth.nim:565](https://github.com/brentp/mosdepth/blob/v0.3.14/mosdepth.nim#L565)–577, 650–678 | levels from largest BAM-header contig; columns 1/2/3; zero-based; BED preset |
| BGZI wrapper | [hts-nim v0.3.31 bgzi.nim:13](https://github.com/brentp/hts-nim/blob/v0.3.31/src/hts/bgzf/bgzi.nim#L13)–19, 82–106; [bgzf.nim:18](https://github.com/brentp/hts-nim/blob/v0.3.31/src/hts/bgzf.nim#L18)–26 | level 1 for both BED files, line and newline written separately, push after writing; flush before index finish |
| CSI wrapper | [csi.nim:34](https://github.com/brentp/hts-nim/blob/v0.3.31/src/hts/csi.nim#L34)–89 | min_shift=14, supplied levels, mapped count, metadata layout |
| Index internals | [htslib 1.22.1 hts.c:2304](https://github.com/samtools/htslib/blob/1.22.1/hts.c#L2304)–2600, [2709](https://github.com/samtools/htslib/blob/1.22.1/hts.c#L2709)–2828; [khash.h:248](https://github.com/samtools/htslib/blob/1.22.1/htslib/khash.h#L248)–355 | bins, linear offsets, compaction, hash iteration order, serialization |
| BGZF layout | [bgzf.c:550](https://github.com/samtools/htslib/blob/1.22.1/bgzf.c#L550)–668, `bgzf_write`, `bgzf_flush`, `bgzf_close`; [bgzf.h](https://github.com/samtools/htslib/blob/1.22.1/htslib/bgzf.h) | 65,280-byte payload, byte-stream block cuts, standard EOF, libdeflate level mapping |
| BED preset | [tbx.c:46](https://github.com/samtools/htslib/blob/1.22.1/tbx.c#L46), `tbx_set_meta` | preset 0x10000, columns 1/2/3, comment 35, skip 0, NUL-separated names |
| CIGAR/flag wrappers | [hts-nim cigar.nim](https://github.com/brentp/hts-nim/blob/v0.3.31/src/hts/bam/cigar.nim), [flag.nim](https://github.com/brentp/hts-nim/blob/v0.3.31/src/hts/bam/flag.nim) | BAM consumption masks and proper/supplementary flag tests |
| BAM support | [hts-nim bam.nim:184](https://github.com/brentp/hts-nim/blob/v0.3.31/src/hts/bam.nim#L184), 364–384; htslib `sam.c:bam_endpos`, `bam_tag2cigar`; `hts.c:hts_set_fai_filename` | stop from reference-consuming CIGAR (minimum one); CG long-CIGAR expansion; FASTA only configures decoding for CRAM |

mosdepth commit: `821fddb12860d024fef4cf0bfe86918f2413d4e4`.
hts-nim v0.3.31: `93ad65442b16d0191aa0c44803d10a6fb7cde9fa`.
htslib 1.22.1: `c1f35d67dd5ff1e226d94abe7b850f28c60f5910`.
The relevant index routines were also read at htslib 1.17
(`6143086502567c5c4bb5cacb2951f664ba28ed6e`).

**Version provenance limitation:** mosdepth's `mosdepth.nimble:10` specifies
`hts >= 0.3.31`, not an exact pin. Its release-tag build workflow names
htslib 1.22.1 and Nim 1.6.18/2.0.2. The source tree does not establish the exact
hts-nim/htslib/Nim revision installed in the golden container. The audited
wrappers and index algorithm are confirmed against the actual golden payloads;
this is not a claim to have recovered that container's dependency lock.

The imported Nim 1.6.18 runtime was also inspected: `system.nim:1470`
(`toInt`: half away from zero), `pure/algorithm.nim:365` (stable mergesort),
`pure/parseutils.nim:432` (signed decimal with underscores after digits), and
`pure/strutils.nim:2356`–2418 (`ffDecimal`: C `%#.*f`, decimal separator `.`).
Rust's fixed decimal formatting preserves the two fractional digits and
round-to-nearest/even behavior of the reference's default C rounding mode.

## Numeric and output rules

* Filter mask **1796 = unmapped | secondary | QC-fail | duplicate**. Default
  mapq is **0**: MAPQ 0 is retained. Supplementary reads are retained, but do
  not participate in mate correction. No insert-length or read-group filter
  is active in either requested form.
* `M`, `=`, `X` add depth; `D` and `N` advance the reference without adding
  depth; `I`, `S`, `H`, `P` add none. Adjacent match spans separated only by
  query-only operations remain continuous. `--fast-mode` counts the entire
  start–stop span, including deletions/skips, and does no mate correction.
* In default mode only proper, non-supplementary pairs are candidates.
  Save a lower overlapping mate by QNAME when the reference IDs match and
  `stop > matepos`; equal-start pairs save only the first occurrence. Otherwise
  take a previously saved QNAME. General CIGARs subtract the intersection of
  covered spans. **Single-CIGAR shortcut bug:** subtract `[rec.start,mate.stop)`
  without checking either opcode or clipping the end to `rec.stop`. A contained
  second mate can therefore erase coverage past its own end. The port retains
  this behavior, with a unit case at first `[2,10)`, second `[5,7)` producing
  coverage only on `[2,7)`. No base qualities or reference bases are examined.
* The chromosome array has **length+1** entries; the last is a sentinel.
  Global stats/histograms and per-base runs exclude it. BED slices clip against
  the array length, **including the sentinel** if a target extends that far.
  Means divide each visited base by the original BED interval length and add
  those float64 terms in order; replacing this with integer sum/length is not
  equivalent at rounding boundaries. Empty slices have mean zero.
* BED is tab-delimited, comments `#` and `track ` lines are ignored, fewer than
  three columns are skipped, and a blank line terminates parsing (the upstream
  `hts_getline(...) > 0` loop). Column four is stripped and becomes the name;
  later columns are ignored. Stable sort is by start only, per chromosome.
  Regions remain separate when duplicated or overlapping: their bases are
  counted repeatedly in region summaries and BED-mode distributions. Unknown
  contigs warn and produce no output. Plain/gzip BED input is supported.
* Windows cover `[0,L)` consecutively, with a short final window and no empty
  final record. They are generated lazily, including for window size one.
  Window region distributions count **one rounded mean per window**, half away
  from zero, saturated at histogram bin 511. BED region distributions instead
  count **individual covered bases**, including overlaps; they grow the
  histogram as needed.
* Global/BED histograms start at 512 bins, grow to `depth+10` when depth is at
  least the current last index, clamp depths **strictly greater than 400000** to
  **399990**, and ignore negative signed depths. Summary min/max are not capped.
  Distribution iteration is descending: skip empty bins above 300, accumulate
  `float64(count)/float64(total)` term by term, suppress cumulative values below
  `8e-5`, and print two decimal places. No histogram mass means no rows at all.
* Summary columns are `chrom length bases mean min max`. Header order is BAM
  header order, with each chromosome immediately followed by its `_region` row;
  `total` then `total_region` follow. Region min/max are base min/max over all
  clipped target slices, not min/max of region means. Length counts overlaps.
  Empty statistics print 0 for min/max and 0.00 for mean.
* **No reads versus all filtered:** the source sets `found` before filtering.
  A queried contig with no alignments has no summary or distribution rows,
  contributes nothing to totals, but has zero-valued BED regions/windows and
  (if enabled) a single zero per-base run. A contig with only excluded reads
  has zero-valued full-length statistics and histogram mass. Per-base output
  lists all BAM-header contigs; regions list only contigs having BED intervals
  (or all positive-length contigs in window mode). Summary/distributions list
  only contigs with query records. Totals always print.
* The optional upstream `-n` BED optimization skips contigs without targets
  entirely; the implementation preserves it. The requested exome form does
  not activate that optimization. No-per-base omits both per-base leaves.
* For BAM, `--fasta` sets `fn_aux`; `hts_set_fai_filename` only passes a reference
  to a decoder when the format is CRAM. There is no reference use in the
  coverage arithmetic. The port's BAM-only CLI therefore needs no FASTA flag.

## BGZF and CSI

The BGZF payload is **65,280 bytes per full block**, independent of lines:
hts-nim writes a line and then its newline, and htslib's `bgzf_write` fills the
current block from either call, splitting the supplied bytes when necessary.
The existing `BgzfWriter` with `Profile { block_payload: 65280 }` has this exact
layout. Both BED files and indexes end with the standard 28-byte empty member.
The compressed bytes are intentionally outside the authorized gate.

Compression identification was done outside the repo, for verification only:
recompressing the first 100 nonempty tier0 md per-base blocks and all 95
nonempty region blocks matched **libdeflate level 1, 100/100 and 95/95**,
using the installed Homebrew libdeflate. Python zlib levels 1–9 each matched
**zero** of these blocks. This is direct compressed-payload evidence, stronger
than inferring the codec from gzip headers. It does not identify the golden's
libdeflate version. All nonempty blocks in the tier0 md CSI files also matched **libdeflate level
7** exactly: **7/7** per-base-index blocks and **23/23** region-index blocks;
again zero matches at zlib levels 1–9. htslib's default index compression level
6 maps to libdeflate level 7 (`bgzf.c`'s `lvl_map`), consistent with this result.
No libdeflate dependency was added.

Index details retained from source:

* CSI magic `CSI\1`; min shift 14; smallest level count covering the largest
  BAM-header contig (all six fixtures: 5). The standalone builder defaults to
  five and accepts `--levels`: BED text cannot reveal an absent longer contig.
* Metadata is seven little-endian uint32 values `[65536,1,2,3,35,0,names_bytes]`,
  then the names of chromosomes actually written, in output order, each NUL
  terminated. Names do not come from hash iteration.
* Push after each complete line, with the previous end as the record start;
  all entries count as mapped. Chromosome boundaries are continuous. Linear
  bins retain the first offset, missing entries fill **backward from the next**
  entry; bin `loff` is computed before compaction.
* Final index offset follows a flush of the BED partial block, pointing at the
  EOF block address (not the last line's within-block end). There is a metadata
  pseudo-bin with start/end offsets and mapped/unmapped counts, and a final
  uint64 zero no-coordinate count.
* Bin order is htslib's integer khash table order: identity hash, triangular
  probing, power-of-two buckets, 0.77 load threshold, in-place rehash kick-outs.
  Sorting bin numbers would fail byte identity. No arbitrary hash seed is used.
* Small bins compact into **existing** parents, bottom-up, when their compressed
  address span is less than **65536 bytes**. Parent chunks are sorted at the
  source's level boundaries, then chunks whose compressed blocks overlap are
  merged. This threshold is a compressed-byte distance, not a block count.

Consequently **shared BED block cuts do not guarantee identical CSI shapes**.
The two-sided rule remains exact: builder(golden) equals golden payload;
our index equals builder(our BED). Cross-file index block sizes and within-block
fields are diagnostics, not additional normalizations. On WGS both md/recal,
the first cross-file payload difference is zero-based byte **79402**, chr1's
`n_bin`: ours **32**, golden **34**. At that contig ours retains bin 9;
the golden instead retains bins 30, 33, 37. The payload lengths are
**380010 / 380634** bytes. htslib's compressed-distance compaction explains
this difference; forcing the golden's bin choices onto our file would violate
part (2). Where bin/chunk shapes differ, a position-by-position sequence of
within-block fields also differs, rather than just compressed block addresses.

## Release verification

`tests/mosdepth_goldens.rs` runs the release CLI on the six sarek BAMs, checks
all three text leaves by byte identity, compares complete decompressed BED
streams using independent `gzip -dc` processes (both exits checked), and checks
the entire raw BGZF ISIZE sequence, including EOF. For each BED, it invokes
`qc csi` on the golden and our file, comparing both decompressed index payloads
against their respective published indexes. It separately reports rebuilt-index
block-size equality and cross-file index block sizes/within-block fields.
No provenance, coordinates, bin, statistic or formatting field is masked.

`NOMICS_RUNS` defaults to `~/omics/runs`; `NOMICS_QC_BIN` selects a release binary;
`NOMICS_MOSDEPTH_BED` overrides `~/omics/data/Agilent_v7.bed`.
`NOMICS_MOSDEPTH_ONLY` optionally selects one run for iteration. Missing inputs
or goldens explicitly skip that BAM; present BAM sizes must match the six
specified completed-copy sizes before use. Failed gates retain outputs in the
test's temporary directory and report the first differing byte.

Source-derived unit cases cover CIGAR gaps/insertions, default versus fast mate
handling, the contained-mate shortcut bug, histogram growth/cap, absent versus
filtered contigs, overlapping/out-of-bounds targets (including the sentinel),
BED stable ties/names/comments/blank termination, and empty statistics.

### Measured gate ledger (2026-11-11, local Mac)

Timings are the final release pass and include all own report/index writes,
excluding golden verification. The complete mosdepth integration gate took
172.58 seconds.

Prefixes are literal: WES is `HG002_WES.{md,recal}`, WGS20 is
`HG002_WGS20.{md,recal}`. Paths are the three runs in the task. Each text
column below names a separate `.mosdepth.<column>.txt` leaf.

| Run | Prefix | global.dist | region.dist | summary | Release wall seconds (BAM → all own reports/indexes) |
|---|---|---|---|---|---:|
| tier0-wes-20261008-r1 | HG002_WES.md | bytes identical | bytes identical | bytes identical | 8.585 |
| tier0-wes-20261008-r1 | HG002_WES.recal | bytes identical | bytes identical | bytes identical | 6.541 |
| wgs20-20261010-r1 | HG002_WGS20.md | bytes identical | bytes identical | bytes identical | 13.909 |
| wgs20-20261010-r1 | HG002_WGS20.recal | bytes identical | bytes identical | bytes identical | 14.546 |
| tier1-wes-20261008-r1 | HG002_WES.md | bytes identical | bytes identical | bytes identical | 35.858 |
| tier1-wes-20261008-r1 | HG002_WES.recal | bytes identical | bytes identical | bytes identical | 24.880 |

Every BED below has identical **complete decompressed bytes** and the same
**entire block-size sequence**, not just the same total decoded size. Block
counts include the empty EOF member. WGS intentionally has no per-base leaves.

| Run | Prefix | BED suffix | Rule b | BGZF blocks |
|---|---|---|---|---:|
| tier0 | HG002_WES.md | per-base.bed.gz | pass | 330 |
| tier0 | HG002_WES.md | regions.bed.gz | pass | 96 |
| tier0 | HG002_WES.recal | per-base.bed.gz | pass | 289 |
| tier0 | HG002_WES.recal | regions.bed.gz | pass | 96 |
| wgs20 | HG002_WGS20.md | regions.bed.gz | pass | 2889 |
| wgs20 | HG002_WGS20.recal | regions.bed.gz | pass | 2889 |
| tier1 | HG002_WES.md | per-base.bed.gz | pass | 30770 |
| tier1 | HG002_WES.md | regions.bed.gz | pass | 99 |
| tier1 | HG002_WES.recal | per-base.bed.gz | pass | 12649 |
| tier1 | HG002_WES.recal | regions.bed.gz | pass | 99 |

For **all ten** indexes, rebuilding the golden gives identical decoded payload
**and** identical BGZF block-size sequence; rebuilding ours gives identical
complete compressed bytes. The final two columns compare our index directly
with the golden, a diagnostic outside the two-sided identity comparisons.
`same` means the entire sequence, including the last block and EOF.

| Run | Prefix | CSI suffix | Rule c (1) | Rule c (2) | Cross-file index block sizes | Cross-file within-block field sequence |
|---|---|---|---|---|---|---|
| tier0 | HG002_WES.md | per-base.bed.gz.csi | pass | pass | same | same |
| tier0 | HG002_WES.md | regions.bed.gz.csi | pass | pass | same | same |
| tier0 | HG002_WES.recal | per-base.bed.gz.csi | pass | pass | same | same |
| tier0 | HG002_WES.recal | regions.bed.gz.csi | pass | pass | same | same |
| wgs20 | HG002_WGS20.md | regions.bed.gz.csi | pass | pass | different | different |
| wgs20 | HG002_WGS20.recal | regions.bed.gz.csi | pass | pass | different | different |
| tier1 | HG002_WES.md | per-base.bed.gz.csi | pass | pass | different | different |
| tier1 | HG002_WES.md | regions.bed.gz.csi | pass | pass | same | same |
| tier1 | HG002_WES.recal | per-base.bed.gz.csi | pass | pass | different | different |
| tier1 | HG002_WES.recal | regions.bed.gz.csi | pass | pass | same | same |

No first difference exists in any gated comparison. For completeness, the
**raw cross-file decoded CSI** comparisons (which are not the gate) differ as
follows; all positions are zero-based. Bytes are shown as unsigned decimal.
The equal-shape cases first differ in compressed block addresses. Different
shapes arise from the source's compressed-distance compaction described above.

| Run / prefix | CSI suffix | First differing byte | ours / golden byte | Payload bytes ours / golden |
|---|---|---:|---|---|
| tier0 / HG002_WES.md | per-base.bed.gz.csi | 107728 | 14 / 133 | 419002 / 419002 |
| tier0 / HG002_WES.md | regions.bed.gz.csi | 235 | 199 / 80 | 1474941 / 1474941 |
| tier0 / HG002_WES.recal | per-base.bed.gz.csi | 84060 | 13 / 112 | 368602 / 368602 |
| tier0 / HG002_WES.recal | regions.bed.gz.csi | 235 | 199 / 80 | 1474941 / 1474941 |
| wgs20 / HG002_WGS20.md | regions.bed.gz.csi | 79402 | 32 / 34 | 380010 / 380634 |
| wgs20 / HG002_WGS20.recal | regions.bed.gz.csi | 79402 | 32 / 34 | 380010 / 380634 |
| tier1 / HG002_WES.md | per-base.bed.gz.csi | 79402 | 47 / 50 | 503386 / 504858 |
| tier1 / HG002_WES.md | regions.bed.gz.csi | 223 | 22 / 181 | 1474893 / 1474893 |
| tier1 / HG002_WES.recal | per-base.bed.gz.csi | 79424 | 174 / 233 | 519642 / 520074 |
| tier1 / HG002_WES.recal | regions.bed.gz.csi | 223 | 22 / 181 | 1474893 / 1474893 |

Validation: `cargo build --release -p nomics`, `cargo test -p nomics-qc`
(13 unit tests, the existing eight VCF gates, and all **18 text + 10 BED +
10 two-sided CSI** mosdepth gates), and `cargo clippy -p nomics-qc -p nomics`.
An additional `--all-targets --no-deps` clippy run checks the new test code.
There are no new warnings; existing warnings in other workspace crates remain.
`nomics-io` was untouched. None of the six BAMs or any golden was modified.
