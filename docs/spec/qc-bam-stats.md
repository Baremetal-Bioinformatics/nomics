# BAM QC: samtools stats 1.24

Status: implemented in `nomics-qc::samtools_stats`. The gate is **all bytes**, with
no normalization, including comments, section order, command names, numeric
formatting and final newline. BAM input only. This is the command sarek 3.10.0
publishes as `reports/samtools/<sample>/<sample>.{md,recal}.bam.stats`:

```text
samtools stats --threads 1 --reference Homo_sapiens_assembly38.fasta HG002_WES.md.bam
```

CLI:

```sh
nomics qc bam-stats --bam X.bam --reference ref.fasta --out X.bam.stats
```

`--input-name` and `--reference-name` override only the two command-line names;
the defaults are the supplied paths' basenames. The header says
`1.24+htslib-1.24`, with two spaces before `stats`. Statistics and report rendering
are independent of file framing: public `Stats::add(&Record)` takes a borrowed
record, and `Stats::render` consumes it. The BAM adapter streams through
`nomics_io::bgzf::ParallelBgzfReader` and `bam::{Header,read_raw_record,RecordView}`.
There are no new dependencies, upstream executables, FFI or unsafe code.

## Source pins and entry points

Both source trees were fetched outside this repository at tag **1.24**:
samtools commit `dc71c7274044d1050ccb64901731373ec7e915b6`, htslib commit
`4b705e4fada8ee2b6b15746f725ee8ac51631803`. No upstream implementation was copied
into the workspace. Literal output vocabulary follows the tool's format strings.

| Behavior | Pinned source |
|---|---|
| Option defaults, parser, stream, final flush | [stats.c:2313](https://github.com/samtools/samtools/blob/1.24/stats.c#L2313), [2701](https://github.com/samtools/samtools/blob/1.24/stats.c#L2701), [2960](https://github.com/samtools/samtools/blob/1.24/stats.c#L2960) |
| Initial dimensions and allocation | [stats.c:2354](https://github.com/samtools/samtools/blob/1.24/stats.c#L2354), [2391](https://github.com/samtools/samtools/blob/1.24/stats.c#L2391) |
| Record selection and counter ordering | [stats.c:1212](https://github.com/samtools/samtools/blob/1.24/stats.c#L1212) |
| Original-read GC, bases, qualities, flags | [stats.c:918](https://github.com/samtools/samtools/blob/1.24/stats.c#L918) |
| CHK | [stats.c:783](https://github.com/samtools/samtools/blob/1.24/stats.c#L783) |
| Indels and unclipped length | [stats.c:419](https://github.com/samtools/samtools/blob/1.24/stats.c#L419), [465](https://github.com/samtools/samtools/blob/1.24/stats.c#L465) |
| Circular coverage | [stats.c:311](https://github.com/samtools/samtools/blob/1.24/stats.c#L311)–391; resize [690](https://github.com/samtools/samtools/blob/1.24/stats.c#L690)–780 |
| GC-depth collection | [stats.c:1414](https://github.com/samtools/samtools/blob/1.24/stats.c#L1414)–1456 |
| GC-depth sorting, interpolation | [stats.c:1516](https://github.com/samtools/samtools/blob/1.24/stats.c#L1516), [1526](https://github.com/samtools/samtools/blob/1.24/stats.c#L1526) |
| Every report format and insert trimming | [stats.c:1542](https://github.com/samtools/samtools/blob/1.24/stats.c#L1542)–1952 |
| Insert histogram storage | [stats_isize.c:180](https://github.com/samtools/samtools/blob/1.24/stats_isize.c#L180)–275 |
| BAM parsing, packed sequence, qualities | [htslib sam.c:786](https://github.com/samtools/htslib/blob/1.24/sam.c#L786), [sam.h:300](https://github.com/samtools/htslib/blob/1.24/htslib/sam.h#L300) |
| Long-CIGAR CG promotion | [htslib sam.c:685](https://github.com/samtools/htslib/blob/1.24/sam.c#L685) |

## `--reference` is not `-r`: no FASTA statistics in this command

This distinction is confirmed both by source and by the goldens' lack of MPC
rows. `--reference` goes through
[sam_opts.c:103–124](https://github.com/samtools/samtools/blob/1.24/sam_opts.c#L103)
into the generic `reference=` input option. htslib records the filename in
[hts.c:1253](https://github.com/samtools/htslib/blob/1.24/hts.c#L1253), but
[hts_set_opt, hts.c:1911](https://github.com/samtools/htslib/blob/1.24/hts.c#L1911)
returns without applying CRAM options for BAM. It does not open the FASTA.
Consequently the BAM implementation deliberately does not load the reference or
require its index, and cannot derive mismatches or reference GC from it.

Stats' separate `-r/--ref-seq` branch at
[stats.c:2776](https://github.com/samtools/samtools/blob/1.24/stats.c#L2776)
loads `info->fai`. That branch activates `read_ref_seq`, `fai_gc_content` and
`count_mismatches_per_cycle` (stats.c:477–633), uppercasing A/C/G/T, converting
other reference characters to zero, fetching bounded intervals through
[htslib faidx.c `faidx_fetch_seq64`](https://github.com/samtools/htslib/blob/1.24/faidx.c),
and computing reference GC excluding unknown bases. That mode is outside this
CLI's contract. A future implementation of it must use `nomics_io::fasta` and
port those buffer/skip/quality rules separately. Reading that FASTA now would
change the requested tool's behavior rather than reproduce it.

## Record categories and summary counters

There are no user filters, regions, read-length restrictions or read-group
splits. All RGs contribute together. No statistic scans unused contigs.

* CHK is updated first, even for secondary/supplementary records and zero-length
  sequences. Each of the three columns is a wrapping 32-bit sum of independent
  zlib-compatible IEEE CRC32 values, initialized with zero per record: name up
  to its first NUL (no NUL included), packed BAM sequence including an odd final
  padding nibble, and raw quality bytes. Missing qualities are literal `0xff`,
  not an empty string. Zero-length sequence adds only the name CRC.
* Secondary increments `non-primary alignments` and immediately returns.
  Supplementary increments its counter but continues. Empty sequence then
  returns. Duplicate count and duplicate length therefore include nonempty
  supplementary records, but not secondary ones.
* Original means neither secondary nor supplementary. Original, nonempty
  sequences supply raw totals, first/last/other counts, total lengths, read
  length histograms, QC-fail and paired-technology counts, GC, quality/base
  cycles, and mapped/unmapped summaries. QC-fail and duplicates stay included.
  Single-end reads are first fragments. Paired reads with neither or both of
  READ1/READ2 are 'other', receiving no fragment GC/quality/base-cycle counts.
* Mapped original reads contribute their entire sequence length, including
  soft clips, to `bases mapped`; mapq zero increments MQ0. Mapped-and-paired
  requires PAIRED and neither UNMAP nor MUNMAP. Proper pairs additionally
  require PROPER_PAIR. Different-chromosome count is the mapped-and-paired
  anomalous read count divided by two with integer truncation.
* MAPQ alone excludes all of UNMAP, SECONDARY, SUPPLEMENTARY, QCFAIL and DUP;
  0 through 255 are distinct bins. These exclusions do not affect CHK.
* All mapped nonsecondary nonempty records, including supplementary, contribute
  NM (integer aux conversion, absent/other types zero), indels, CIGAR bases,
  coverage and GC-depth. NM is not recomputed from the FASTA. `bases mapped
  (cigar)` counts M, I, =, X, not D, N or clipping. `bases trimmed` is zero.
* Length histograms and maxima use sequence length **plus hard clips**; total
  lengths and cycle counts use stored sequence length. Maxima also include
  supplementary alignments. HTSlib rejects mapped CIGAR/query-length disagreement;
  the adapter does likewise. Long CIGARs encoded with a placeholder and CG:B:I/i
  are promoted before collection, preserving the original sequence/quality CRC.

## Cycle statistics and GC histograms

`ngc=200`, `nquals=257`, initial read capacity 300; the indel-length output limit
stays **300** even when read buffers grow. For recognized first/last reads, exact
C/G contribute GC. Increment integer bins in
`[gc_count*199/seq_len, min((gc_count+1)*199/seq_len,199))`. This interval can be
empty for long reads. The GC output coalesces consecutive equal counts, printing
the midpoint of each run when its count changes, with two decimal places. It
never emits the final run by itself.

Reverse reads reverse the cycle index but FBC/LBC/GCC do **not** complement their
bases. GCT additionally complements A/C/G/T. N and ambiguous/equals bases go to
N and O respectively, excluded from the ACGT denominator. N/O percentages thus
use the ACGT sum, not the total number of bases. GCT has only four columns.
FTC/LTC sum five raw counts across cycles; O is omitted.

Qualities use the reversed cycle for reverse reads, retain all scores including
255, and average the scores over total original sequence length (including
'other' reads, which supplied no quality sum). FFQ/LFQ output all qualities
through one above the maximum, except 255 has no extra column; there is no
per-cycle suppression. Rows extend to each fragment maximum including hard clips,
so trailing zero rows can occur. Base rows suppress cycles with no ACGT.

All mapped nonsecondary indels count once per CIGAR operation, not per base.
The cycle cursor advances for query-consuming operations, excludes D/N/H/P, and
ignores zero-length operations. Insertion cycle is forward cursor or reverse
`seq_len-cursor-insertion_len`. Deletion cycle is forward `cursor-1` or reverse
`seq_len-cursor-1`; negative cycles are discarded including the length count.
IC's forward/reverse labels actually mean first/last fragment; single-end is
first. ID suppresses empty lengths; IC suppresses empty cycles.

Default barcode pairs BC/QT, CR/CY, OX/BZ and RX/QX are also handled: original
first fragments only, string tags, first observed length fixes capacity, longer
later tags skipped, shorter tags update their existing prefix. A/C/G/T/N are
counted case-sensitively; the first other character fixes the separator position,
and a subsequent inconsistent separator aborts that record's remaining tag
processing. Separator positions are omitted. Phred+33 quality length must equal
the barcode length; absent qualities leave the quality rows without columns.
These optional sections appear between LTC and IS when tags are present.

## Insert sizes and floating point

Default dense insert histogram has indices 0–8000. Both mapped original mates
count independently; absolute TLEN is capped at 8000. TLEN zero on different
contigs is omitted, but nonzero TLEN there still counts. Same-strand mates are
'other'. Opposite strands use READ1, reverse flags and signed mate displacement
to select inward/outward; exact overlap is assumed inward.

At output **each orientation at each size** is halved independently with integer
truncation. Total orientation summaries span every bin. The main bulk includes
bins through the first cumulative fraction **strictly greater than the float32
representation of 0.99**, including that bin; if never reached, output ends at the
last occupied bin. IS includes intervening empty bins, starting at zero. Average
uses that bulk only. Standard deviation divides by the bulk count, but its sum
starts at size **1**, omitting zero's variance contribution. These are preserved
quirks, not statistical corrections.

SN error rate and average lengths cast numerator to float32 and perform the
division in float32 before formatting. Proper-pair percent multiplies the integer
numerator by 100 first, then divides in float32. Average quality and insert
moments use double. Error rate is C `%e` (six fractional digits, signed exponent
with at least two digits); lengths use `%.0f`, other SN real fields `%.1f`.
Rust's fixed decimal formatting retains round-to-nearest-even behavior. Base
percentages and GC midpoints use double with two decimals. No locale or platform
NaN exemption is needed: report divisions guard empty denominators except the
internal insert cut-off, whose NaN merely fails its comparison.

## Coverage ring and GC-depth quirks

Coverage starts with 1500 int32 slots (5 × 300). Before adding a mapped read,
flush positions strictly before its start. Insert only M/=/X; D/N advance the
reference cursor; I/S/H/P do not contribute depth. Mate overlaps and duplicate,
QC-fail and supplementary alignments are retained. Zero-depth positions are
never counted. Default output bins are `[1-1]` through `[1000-1000]`, overflow
`[1000<]` labeled 1000. Only occupied bins print.

The actual ring semantics are reproduced, including two upstream mistakes:

* Full flush sets its local endpoint to `old_pos+size-1`, then flushes only
  through endpoint-minus-one, leaving the last slot. After a large gap, its new
  physical start is computed from that shortened endpoint while the logical
  position is set to the new requested position. A contig change requests -1
  and resets the physical start to zero. Do not replace this with a conventional
  pileup implementation.
* When unclipped read length reaches capacity, read buffers grow to
  `old+2*(1+length-old)`. The coverage ring grows to `length*5`, but its two
  memcpy calls at stats.c:772–774 use **byte counts**, despite calculating them
  in int32 slots (and the second destination is displaced in slots). Safe Rust
  reproduces those byte transfers using native-endian bytes. It does not fix
  the upstream loss of part of the active pileup.

Sort checking compares a mapped record's position with the last mapped position
only when its tid matches the current GC-depth tid. The first detected inversion
permanently disables coverage/GC-depth collection and output. It does not check
header sort order, monotonic contig IDs or unmapped positions.

Without `-r`, GC-depth begins with an empty bin at index zero. New bins begin at
the first mapped record, a contig change, or a start displacement **>20000**
(not >=). The entire read goes to one bin regardless of overlap with its end.
Depth is the count of mapped records, not base coverage. GC is a float32 sum of
per-original-read `gc_count/seq_len`; supplementary records contribute depth but
zero GC. Bins are sparse; reference-only windows/contigs never appear.

Before sorting, only indices **less than** the final index are converted to GC
percent, dividing the accumulated float32 GC by depth in double, then `rint`
(ties to even). The final bin stays an unscaled GC sum. Sort **all** bins by GC,
then depth, but output only **count-minus-one** sorted bins: this may omit a
rounded bin and emit the unscaled final bin. The empty initial bin stays included.
Group adjacent GC values differing from the group's first value by <0.1.
Unique-sequence percent is `(group_end+1)*100/total_bin_count`, printed `%.3f`.

For each group's depth percentile p, compute float32 `p*(N+1)/100`, truncate to k,
clamp k<=0 to first depth and k>=N to last depth; otherwise interpolate depths
k-1 and k with float32 arithmetic (the depth difference is uint32 subtraction).
Multiply by float32 average read length and divide by 20000 in float32. GC prints
`%.1f`; the five depths use `%.3f`. These boundaries, including the off-by-one
round/sort/output ranges, are preserved without normalization.

## Validation

`tests/bam_stats_goldens.rs` follows the VCF release-gate convention:
`NOMICS_RUNS` defaults to `~/omics/runs`, `NOMICS_QC_BIN` selects an alternate
release binary, absent BAM/goldens explicitly skip, and present fixtures must
have the supplied exact sizes before use. It invokes the CLI and literal `cmp`
for each file; failure reports the first differing line. Inputs remain read-only.
Unit cases cover category ordering, single-end/reverse/missing qualities,
ambiguous bases, hard/soft clips, indels, lengths >300, fragment histograms,
barcodes, percentile endpoints/interpolation, float exponent formatting,
strict insert trimming/zero-size variance, and ring wrap/overlap/overflow/full
flush. The six fixture results and wall times are recorded below.

## CRAM follow-up

CRAM input is out of scope. A CRAM adapter must supply this same record boundary,
matching htslib 1.24's reference-backed reconstruction of packed sequence and
qualities, exact read names/flags/coordinates/CIGAR/NM/tags, and odd-length packed
sequence padding for CHK. Use `nomics_io::fasta` for required reference access;
confirm CRAM decoder behavior against htslib before claiming checksum identity.
Do not activate stats' `-r` calculations merely because a CRAM decoder uses a
reference. New CRAM goldens (both stages) and source audit are required; renaming
a BAM stats file to `.cram.stats` is not a gate.

## Six release-CLI comparisons (2026-11-11, arm64 Mac)

Every output is byte-identical (`cmp`, no normalization or differing line):

| Run | File | Wall time |
|---|---|---:|
| `tier0-wes-20261008-r1` | `HG002_WES.md.bam.stats` | 0.569 s |
| `tier0-wes-20261008-r1` | `HG002_WES.recal.bam.stats` | 0.493 s |
| `wgs20-20261010-r1` | `HG002_WGS20.md.bam.stats` | 9.320 s |
| `wgs20-20261010-r1` | `HG002_WGS20.recal.bam.stats` | 12.317 s |
| `tier1-wes-20261008-r1` | `HG002_WES.md.bam.stats` | 54.363 s |
| `tier1-wes-20261008-r1` | `HG002_WES.recal.bam.stats` | 20.481 s |

Each time covers one CLI invocation, streaming the entire BAM, calculating every
section and writing the report. The 5,047,709,513-byte tier1 markduplicates BAM
completed in under one minute. These are local single-pass wall times, not a
speedup comparison or a saturated-CPU benchmark.

Validation passed: `cargo build --release -p nomics`, `cargo test -p nomics-qc`
(21 unit tests, all six BAM comparisons, the existing eight VCF comparisons, and
doc tests), `cargo clippy -p nomics-qc -p nomics` with no new warnings, and
`cargo clippy -p nomics-qc --all-targets --no-deps -- -D warnings`. The broader
clippy command still reports existing warnings in other crates, and the release
build retains the existing `nomics-hc` deprecated `fetch_update` warning.
`git diff --check` passes. No golden, BAM, reference or other file under
`~/omics` was changed.

No residual byte differences or unconfirmed rules remain for the requested BAM
command. The supplied fixtures confirm the six real-data reports; source-derived
unit inputs cover the additional categories and edge conditions listed above.
They are not substitute upstream goldens for other tools/options or CRAM, which
requires its own pass. Port and validation took roughly ten minutes.
