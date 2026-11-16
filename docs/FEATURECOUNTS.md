# featureCounts 2.1.1 / nf-core/rnaseq 3.27.0

`nomics-count` implements the alignment-counting path used for biotype QC by
nf-core/rnaseq 3.27.0. The binary entry point is `nomics featurecounts`; the
library supports reusable annotations, streamed BAM counting, and TSV writers.

## The counting unit matters

The actual pipeline invocation is:

```sh
nomics featurecounts -B -C -g gene_biotype -t exon -p -T 6 \
  -a genome_gfp.gtf -s 2 -o SAMPLE.featureCounts.tsv SAMPLE.markdup.sorted.bam
```

In Subread 2.1.1, **`-p` declares paired input; it does not count fragments**.
Only `--countReadPairs` switches on fragment assignment. The pinned nf-core
module does not supply it. The C oracle explicitly reports “The reads are
assigned on the single-end mode.” Both ends of a pair are independently
counted. Changing this to fragment counting would fail the reference gate.

Consequences confirmed in `readSummary.c::process_line_buffer`:

- Unmapped alignments are `Unassigned_Unmapped` before subsequent filters.
- With `-p -B`, an unpaired record in the mixed library is
  `Unassigned_Singleton`. Paired records with an unmapped mate still count in
  this mode. No paired reads at all with `-p` is an input error.
- The paired-fragment `-B` and `-C` checks require
  `is_paired_end_mode_assign`, which is false here. `-C` is inactive, even for
  interchromosomal or same-strand mates. There is no fragment-length filter.
- NH > 1 rejects each alignment as `Unassigned_MultiMapping`; absent or
  noninteger NH defaults to 1. Subread zero-extends c/s NH values, and treats
  both i/I as signed 32-bit values; that observable quirk is retained.
- Default minimum MAPQ is 0. Duplicates, secondary and supplementary
  alignments are included unless rejected by another filter, notably NH.
- R2 (flag 0x80) reverses the read strand before comparison, in either input
  mode. `-s 2` requires inequality against the annotation strand; unknown
  annotation strand (255 in C) therefore passes reverse-strand comparisons.
- M/= /X bases overlap features; I/D/N separate sections, D/N consume reference,
  and soft/hard clips and padding contribute no overlap. At least one base
  normally suffices. The C parser's empty section at I/D/N can hit an exon
  spanning the section boundary; this zero-length-section quirk is retained.
  Section lengths wrap at 16 bits, as in the C parser. The 2.1.1 default section
  limit is **65,536**, not the older advertised 10.
- Multiple exons of the same grouping value form one meta-feature. A read
  overlapping multiple distinct meta-features is `Unassigned_Ambiguity`; there
  is no largest-overlap preference, fractional counting or `-O` assignment.

The two RAP1_UNINDUCED fixture samples are single-end: their golden commands
omit `-p`. The remaining three samples include `-p`. Tests reproduce these
recorded invocations rather than silently changing input mode.

## Annotation and output

GTF selection defaults to exon; `-t` also accepts comma-separated types.
Grouping supports `gene_biotype`, `gene_id`, or another attribute. Attribute names
are matched exactly; spaces and GFF `=` separators, quoted semicolons, repeated
keys (last wins), and the feature-name length limit are handled. Malformed
annotation or BAM input returns an error; C crashes on some malformed inputs
are not reproduced. Input is uncompressed GTF and BGZF BAM through `nomics-io`.

Group rows follow first encounter in the annotation, chromosomes follow first
encounter, and exons are ordered by start using the reference's tie rule
(selection-sort leaves of <=11, right half first on equal merge keys).
Within a group, output intervals are partitioned by chromosome and strand.
**The lists retain overlaps and duplicate exons; Length is the union length
within each chromosome/strand partition.** It is not a merged interval list. The reference's signed 32-bit Length
wraparound for very large meta-features is retained and oracle-tested.
Chromosome matching tries exact names, removes a BAM `chr` prefix, or adds
`chr` only to BAM names of at most two bytes, matching the C lookup.

All 14 summary categories are emitted in reference order, including zero rows
for filters inactive in this mode. Integer counts and annotation order are
independent of inflate thread count. BAM input streams in batches; no whole-BAM
buffer or name-pairing table is retained. CPU assignment is serial, BGZF inflate
is parallel. No speedup claim or GPU path is made.

The comment is exactly `# Program:featureCounts v2.1.1; Command:` followed by
**every argument double-quoted and followed by a space**, including the final
argument. Subread does not escape quotes inside arguments. The CLI uses the
compatibility executable label `featureCounts`, then preserves its arguments
in order, omitting only the nomics-specific `--mqc-prefix` option and value.
Paths in the comment and count/summary column headers are emitted as supplied.
See `GATES.md` for the narrowly allowed path normalization. Staging the same
basenames gives fully identical files without normalization.

The supported scope deliberately excludes `--countReadPairs`, SAF/SAM input,
read assignment detail files, junction counting, and additional filtering or
fractional-counting options. Unsupported CLI flags are rejected by clap.

## MultiQC and API

Add `--mqc-prefix SAMPLE` for the two nf-core biotype tables (single BAM only).
The sample label is taken from the input basename before `.markdup.sorted.bam`.
The API `write_multiqc` accepts an explicit sample label. The conversion follows
the pipeline's `mqc_features_stat.py`: first and seventh columns, default maximum
100 categories, rRNA percentage over assigned biotype counts, missing rRNA
rendered as integer `0`, no rRNA table for zero total counts.
`src/biotypes_header.txt` is the verbatim MIT-licensed nf-core asset; its
notice is retained in `crates/nomics-count/THIRD_PARTY.md`.

Library entry points:

- `Annotation::read` / `Annotation::from_gtf`: load once, reuse across samples.
- `count_bam(&annotation, path, &CountOptions)`: integer `Counts` by group and
  14-row summary. `CountOptions::default()` is the nf-core PE reverse-strand
  configuration; set `paired_end: false` for SE libraries.
- `write_table`, `write_summary`: deterministic writers with caller-supplied
  displayed input names and provenance argv (including argv[0]).
- `run(&Options)`: count one or more BAMs and write output plus `.summary`.
- `write_multiqc`: optional biotype tables.

## Validation (2026-10-15, Mac fixture)

Fixture: `$NOMICS_FIXTURES/rnaseq-test/`, read-only.

| Sample | Assigned | featureCounts TSV | Summary | Both MultiQC TSVs |
|---|---:|---|---|---|
| RAP1_IAA_30M_REP1 | 83,236 | identical | identical | identical |
| RAP1_UNINDUCED_REP1 | 42,946 | identical | identical | identical |
| RAP1_UNINDUCED_REP2 | 85,983 | identical | identical | identical |
| WT_REP1 | 162,121 | identical | identical | identical |
| WT_REP2 | 81,624 | identical | identical | identical |

**20/20 files byte-identical, including provenance comments. No normalization
needed with the original staging basenames.** Library counts match at 1 and 6
inflate threads for every sample. The CLI test stages read-only fixture inputs
via symlinks in a temporary directory and compares all bytes.

A separate `-g gene_id` run on all five BAMs matches the freshly built 2.1.1 C
oracle: 125 feature rows per sample and all summary rows identical, normalizing
only oracle executable spelling and output paths in the provenance comment.
No differences in counts, annotation metadata, or row order.

Tests:

```sh
cargo test -p nomics-count -- --nocapture
cargo test -p nomics --test featurecounts -- --nocapture
cargo clippy -p nomics-count --all-targets --no-deps -- -D warnings
# Optional source oracle (not linked, called, or shipped by production code):
NOMICS_FEATURECOUNTS_ORACLE=/path/to/subread-2.1.1/featureCounts \
  cargo test -p nomics-count --test oracle -- --nocapture
```

Set `NOMICS_RNASEQ_FIXTURE` to override the fixture directory. An absent default
fixture skips with a message; an explicitly configured absent directory fails.
The optional differential test constructs its own BAM/GTF and checks both
biotype/gene grouping at -s 0, 1, 2 against the C oracle, including equal-start
exons, unsorted input, unknown strand, mixed read types, splicing/indels,
duplicates, secondary/supplementary flags, and signed NH types. Unit tests cover
filter precedence, ambiguity, union lengths, command spacing and summary order.

The focused library, oracle, CLI fixture, and strict clippy checks pass with
**the pinned Rust 1.98.1**. The official distribution was checksum-verified and
installed only under `/tmp` because this Mac's system Cargo uses Rust 1.99.0
and has no rustup. No toolchain configuration was changed.

Broader workspace checks expose existing issues outside this change:
`cargo test --workspace` fails in `nomics-hc/tests/annotate.rs:158` because its
external `t0_calls.tsv.gz` fixture (under `$NOMICS_FIXTURES`) is absent. Workspace clippy
reports an empty `1..=GPU_MAX_READ_LEN` range in the non-CUDA HaplotypeCaller
stub at `nomics-hc/src/caller.rs:328`. The changed crate and CLI pass strict
clippy with `--no-deps`; these unrelated files are unchanged.

## Source provenance

The requested GitHub `2.1.1` tag was absent (the mirror advertises only 2.0.2).
The spec and test-only oracle used the official
[Subread 2.1.1 source release](https://sourceforge.net/projects/subread/files/subread-2.1.1/subread-2.1.1-source.tar.gz/download),
SHA-256 `6392d7c66831cdd767e58251892a79a51b6fab8ed0ba9671ad5e85ff1ab01eaa`.
GPL source was read as a behavioral specification; no C implementation is
vendored, linked or invoked by nomics. Fresh Rust uses only `std` and `nomics-io`.

Relevant release source locations:

| Source | Behavior |
|---|---|
| `readSummary.c:6647-6858` | CLI defaults, argv quoting, -p vs countReadPairs |
| `readSummary.c:985-1410` | GTF loading, group attributes, invalid/missing fields |
| `HelperFunctions.c:505-650` | attribute grammar and trimming |
| `readSummary.c:1446-1660`, `core.c:4718-4755` | annotation sort and equal-start order |
| `readSummary.c:2296-2440` | BAM CIGAR and NH parsing |
| `input-files.c:3477-3497` | integer NH tag conversion quirks |
| `readSummary.c:2888-3100` | filter precedence, PE declaration, strand comparison |
| `readSummary.c:3240-3350`, `:3530-3965` | chromosome aliases, overlaps and assignment |
| `readSummary.c:4615-4895`, `interval_merge.c:29-98` | metadata, union lengths, TSV and summary formatting |

Pinned pipeline sources:
[featureCounts module](https://github.com/nf-core/rnaseq/blob/3.27.0/modules/nf-core/subread/featurecounts/main.nf),
[biotype script](https://github.com/nf-core/rnaseq/blob/3.27.0/modules/nf-core/custom/multiqccustombiotype/templates/mqc_features_stat.py),
[MIT header asset](https://github.com/nf-core/rnaseq/blob/3.27.0/assets/biotypes_header.txt).
