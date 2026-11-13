# STARsolo oracle counting milestone

Implemented 2026-10-16 in `nomics-solo`, exposed by `nomics solo-count`.
The reference is STAR **2.7.10b**, using nf-core/scrnaseq **4.2.0**'s
10x v3 command. This is a counting milestone: no STAR aligner is implemented.

## Measured gate

Fixture: `$NOMICS_FIXTURES/starsolo-pbmc1k/`. The regenerated oracle
includes `NH HI AS nM CR CY UR UY GX GN gx gn CB UB` and unmapped reads
(`--outSAMunmapped Within`). The lead confirmed its Solo outputs equal the
independent decision-run goldens. The original uppercase-only oracle was
insufficient for multi-gene statistics; the regenerated fixture closes that gap.

The release integration test compared every byte, without sorting, numeric
rounding, header deletion, or other content normalization:

| File relative to Solo.out | Result |
|---|---|
| `Barcodes.stats` | identical |
| `Gene/Features.stats` | identical |
| `Gene/Summary.csv` | identical |
| `Gene/UMIperCellSorted.txt` | identical |
| `Gene/raw/matrix.mtx` | identical |
| `Gene/raw/barcodes.tsv` | identical |
| `Gene/raw/features.tsv` | identical |
| `Gene/filtered/matrix.mtx` | identical |
| `Gene/filtered/barcodes.tsv` | identical |
| `Gene/filtered/features.tsv` | identical |

No residual differences. The fixture's files are uncompressed. The test also
accepts gzipped golden files, comparing decompressed bytes as specified in
`GATES.md`.

Measured library runtime: **69.669 seconds**, six BAM inflate threads, on this
Mac, including reading the whitelist, reading the full BAM, grouping reads,
barcode matching, UMI collapse, cell calling, and writing all ten files.
The complete fixture test took **70.01 seconds**, excluding compilation.
After adding the zero-count tentative-barcode knee regression, the final
implementation repeated the full gate in **44.437 seconds** (fixture test
**44.78 seconds**) with warm input caches. All ten files were again identical.
These are individual observations, not saturated CPU benchmarks or speedup claims.

- Reads: **66,601,887**.
- Counted unique-gene reads: **34,520,418**.
- Raw UMIs: **10,473,444** across **225,769** nonempty barcodes.
- Raw matrix: **78,899 × 6,794,880**, **3,628,292** nonzero entries.
- CellRanger2.2 cutoff: **2,304** UMIs; **1,119** called cells.
- UMIs in called cells: **9,389,330**.

Verified output retained at
`/var/folders/26/rpv6m28d2_g1xf44f1f58h5r0000gp/T/nomics-solo-gate-80141`.

## Usage and input contract

```sh
cargo build --release -p nomics
fixture=$NOMICS_FIXTURES/starsolo-pbmc1k
# Run only after this marker exists; never read a partially copied BAM.
test -f "$fixture/ORACLE-READY"
target/release/nomics solo-count \
  --bam "$fixture/oracle.bam" \
  --whitelist "$fixture/10x_V3_barcode_whitelist.txt.gz" \
  --features "$fixture/golden-Solo.out/Gene/raw/features.tsv" \
  --output /private/tmp/nomics-pbmc1k-Solo.out \
  --threads 6
```

The output directory must not exist. Temporary files live beside the output
and are removed on completion or ordinary error. Allow several GB of temporary
disk space. Inputs are read-only. Whitelist and features inputs can be plain
text or gzip; outputs are plain text.

The supported counting configuration is fixed to the milestone:
`CB_UMI_Simple`, 16-base CB, 12-base UMI, `Gene`, `1MM_multi`, `1MM_All`,
`soloUMIfiltering -`, `soloMultiMappers Unique`, and
`CellRanger2.2 1000 0.99 10`. Other protocols and filtering options are not
exposed. The oracle must contain every input read and every reported alignment,
including unmapped reads, unique read names, and ordinary single-end records.
Paired/supplementary records, missing required tags, unknown feature IDs, and
incomplete valid-barcode multimapper groups are errors.

**Authorized feature-list exception:** `--features` takes the reference
`features.tsv` in its original order. No counts, statistics, cell barcodes, or
other golden values are read by the implementation. GTF derivation is deferred
because the filtered GTF is absent on this Mac. STAR assigns gene indices on
first occurrence in accepted exon records in `GTF.cpp`, writes `geneInfo.tab`
in that order, and outputs features in that order. It substitutes the gene ID
for an empty gene name. The default third output column is the literal
`Gene Expression`, **not** the GTF gene biotype. This implementation validates
and preserves the supplied three-column feature list.

## Design and source ledger

Fresh safe Rust; only `nomics-io` and `nomics-mem` dependencies. No external
bioinformatics/compression library or reference executable is invoked.
The large whitelist, prior counts, packed read array, per-cell numeric arrays,
and filtered barcode ranks use `HostArray` backed by `HostBuf`.

1. Stream the BAM through `nomics-io::bamstream::BamStream`. Primary alignment
   flags select one copy of each read's raw barcode/UMI and RNA qualities for
   statistics. `NH` determines genomic uniqueness. Lowercase `gx` supplies
   per-alignment gene evidence; the union across all alignments gives the read's
   gene set. At most two distinct gene IDs need to be retained because Unique
   counting and feature statistics only distinguish zero, one, or multiple.
   Multimappers are hash-partitioned to disk by exact read name and reunited
   before feature counting; the observed record count must match NH. Uppercase
   GX/GN is not sufficient: it writes `-` for both multi-gene and no-gene reads.
2. Recompute barcode matching from CR/CY and UMI validation from UR. UY and CY
   contribute the CB+UMI Q30 statistic; BAM RNA qualities supply the RNA Q30
   statistic. UMIs containing non-ACGT bases and homopolymers are rejected,
   overriding the CB match category. Exact priors include only valid-UMI reads,
   regardless of their gene/mapping status. Whitelist barcodes are sorted and
   deduplicated, then indexed by a 16-bit prefix lookup and binary search.
3. Defer all unique-gene read correction until exact priors are complete.
   A unique 1MM candidate requires a nonzero exact prior. Multiple candidates
   use the source's candidate iteration order, Phred cap 33, and single-precision
   product/sum/maximum and threshold `0.975f32`. Error probability itself uses
   double-precision `pow`, then rounds the weighted probability to float.
   A single N in a barcode permits only one whitelist candidate; multiple Ns
   and ambiguous N matches are rejected. No pseudocounts are used.
4. Pack accepted CB index / gene index / raw UMI into 64 bits (23/17/24), sort,
   and process each cell/gene group. `1MM_All` counts connected components of
   the undirected Hamming-distance-one graph, including transitive chains.
   Abundance does not change this count. Default UMI filtering leaves the same
   UMI in separate genes independent. This is implemented with union-find,
   not STAR's graph-color bookkeeping; the graph and counts are identical.
5. Call cells from the descending UMI counts using zero-based robust-max rank
   `round(1000 * (1 - 0.99))`, cutoff `max(1, round(max / 10))`, and inclusive
   cutoff ties. Tentative barcodes with no surviving reads still occupy zero-count
   positions in the ranking. Means use integer division; UMI median uses descending order,
   while read/gene medians use ascending order, all at index `nCells / 2`.
6. Write matrices in whitelist-column / original-gene-row order, all whitelist
   barcodes in raw output, filtered barcodes in whitelist order, exact stats
   field widths, and six-significant-digit C++-style summary formatting.

The production implementation never reads **CB or UB**. The synthetic test
supplies deliberately wrong CB/UB tags and verifies all ten outputs are
unchanged when the tags are removed.

Pinned source specifications:

- [Barcode extraction, matching, UMI precedence](https://github.com/alexdobin/STAR/blob/2.7.10b/source/SoloReadBarcode_getCBandUMI.cpp).
- [Parameter initialization and defaults](https://github.com/alexdobin/STAR/blob/2.7.10b/source/ParametersSolo.cpp),
  [default values](https://github.com/alexdobin/STAR/blob/2.7.10b/source/parametersDefault),
  [floating-point field types](https://github.com/alexdobin/STAR/blob/2.7.10b/source/ParametersSolo.h).
- [BAM GX/gx output](https://github.com/alexdobin/STAR/blob/2.7.10b/source/ReadAlign_alignBAM.cpp),
  [gene evidence](https://github.com/alexdobin/STAR/blob/2.7.10b/source/Transcriptome_classifyAlign.cpp),
  [read feature accounting](https://github.com/alexdobin/STAR/blob/2.7.10b/source/SoloReadFeature_record.cpp).
- [Posterior correction and rejected-read accounting](https://github.com/alexdobin/STAR/blob/2.7.10b/source/SoloReadFeature_inputRecords.cpp).
- [Per-cell grouping and counts](https://github.com/alexdobin/STAR/blob/2.7.10b/source/SoloFeature_collapseUMIall.cpp),
  [1MM_All graph](https://github.com/alexdobin/STAR/blob/2.7.10b/source/SoloFeature_collapseUMI_Graph.cpp).
- [Cell calling and medians](https://github.com/alexdobin/STAR/blob/2.7.10b/source/SoloFeature_cellFiltering.cpp),
  [matrix output](https://github.com/alexdobin/STAR/blob/2.7.10b/source/SoloFeature_outputResults.cpp),
  [summary output](https://github.com/alexdobin/STAR/blob/2.7.10b/source/SoloFeature_statsOutput.cpp).
- [Gene index assignment](https://github.com/alexdobin/STAR/blob/2.7.10b/source/GTF.cpp),
  [geneInfo output](https://github.com/alexdobin/STAR/blob/2.7.10b/source/GTF_transcriptGeneSJ.cpp).

## Validation

```sh
cargo test -p nomics-solo --lib --test synthetic
cargo test --release -p nomics-solo --test fixture -- --nocapture
cargo clippy -p nomics-solo --all-targets --no-deps -- -D warnings
cargo clippy -p nomics-solo -p nomics --all-targets --no-deps
```

The full fixture test automatically runs when the default fixture has
`ORACLE-READY`; set `NOMICS_SOLO_FIXTURE` to use another location. Missing/not-ready
optional fixtures are reported as skipped. On differences, it prints the first
differing line of **every** requested file and retains the generated outputs.
It never uses golden counts as inputs.

Six unit tests cover invalid-UMI precedence, barcode N ambiguity, exact priors,
posterior float behavior, transitive UMI graph edges in both halves, knee rank
and rounding, and summary precision. The synthetic integration test additionally
covers multi-alignment gene unions (including a gene-free primary alignment),
multi-gene exclusion, single-mismatch correction, unmapped reads, and CB/UB
independence. All pass. The full ten-file fixture gate passes.

Strict Clippy passes for the new crate (`--no-deps`). Clippy for the new crate
and CLI passes with an existing `nomics-hc` deprecated `fetch_update` warning.
Running strict Clippy through dependencies also encounters existing
`nomics-io` `manual_is_multiple_of` and `type_complexity` warnings; those
unrelated source files are unchanged. Release CLI build and `solo-count --help`
were verified.
