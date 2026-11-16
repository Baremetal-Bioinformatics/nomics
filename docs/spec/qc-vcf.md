# VCF QC: implemented and gated

Status: **implemented in `nomics-qc`, eight release-CLI gates pass**. The
previous worker stopped at commit `433d696` because the original byte-only
contract could not define upstream heap contents. The project lead subsequently
authorized exactly the `.TsTv.qual` field exemptions now recorded in
`docs/GATES.md`, with invalid boundary elements taken as `0.0`, and ignoring
`output_TsTv_by_count`'s invalid maximum-count write. No other normalization is
used. The source/golden evidence that established that blocker is retained below.

CLI: `nomics qc vcf-stats --vcf X.vcf.gz --outdir DIR [--prefix P]`.
The default prefix removes `.vcf.gz`. All four reports are written under DIR;
bcftools' command and ID lines use the input basename, including `.vcf.gz`,
regardless of P. No changes to pipeline execution or other QC tools are included.

## Upstream entry points

| Report | Exact source | Function |
|---|---|---|
| `.TsTv.count` | [vcftools v0.1.17, variant_file_output.cpp](https://github.com/vcftools/vcftools/blob/v0.1.17/src/cpp/variant_file_output.cpp#L3145) | `variant_file::output_TsTv_by_count` |
| `.TsTv.qual` | [vcftools v0.1.17, variant_file_output.cpp](https://github.com/vcftools/vcftools/blob/v0.1.17/src/cpp/variant_file_output.cpp#L3229) | `variant_file::output_TsTv_by_quality` |
| `.FILTER.summary` | [vcftools v0.1.17, variant_file_output.cpp](https://github.com/vcftools/vcftools/blob/v0.1.17/src/cpp/variant_file_output.cpp#L2867) | `variant_file::output_FILTER_summary` |
| `.bcftools_stats.txt` | [bcftools 1.23.1, vcfstats.c](https://github.com/samtools/bcftools/blob/1.23.1/vcfstats.c) | `do_vcf_stats`, `print_header`, `print_stats`; variant typing in [htslib 1.23.1, vcf.c](https://github.com/samtools/htslib/blob/1.23.1/vcf.c), `bcf_get_variant_types` |

These sources were retrieved at the named tags. Supporting source audited for
this implementation: vcftools `vcf_entry.cpp`, `vcf_entry_setters.cpp`,
`entry_getters.cpp`, `header.{h,cpp}`; bcftools
[dist.c](https://github.com/samtools/bcftools/blob/1.23.1/dist.c); htslib
[vcfutils.c](https://github.com/samtools/htslib/blob/1.23.1/vcfutils.c),
`bcf_calc_ac`. Only report vocabulary/format strings are reproduced verbatim;
the algorithms are fresh safe Rust, with no added dependencies.

## `.TsTv.qual`: observable undefined behavior

Let `N` be the number of distinct recognized biallelic SNP quality values.
At line 3324, both above-threshold cumulative vectors have **N+1** elements:
their valid indices are **0 through N**. The reverse accumulation begins with
`idx = N`. Lines 3331–3332 read element **N+1** on its first iteration. That
element is outside each vector. The final output iteration at line 3344 also
reads element **N+1**. C++ `vector::operator[]` does not check these accesses.

Consequently the above-threshold counts and ratios depend on memory outside
the initialized vectors. Zero-extending the vectors would make a safe,
deterministic algorithm, but would change the tier0 golden's bytes.

The supplied tier0 report has 800 quality bins (801 lines including its
header). Its final two lines, shown with literal TAB separators, are:

```text
9056.06	674	251	2.68526	3.17042e-320	1	3.17042e-320
9501.06	675	251	2.68924	3.17042e-320	7.5261e-320	0.421256
```

Thus line **800**, column **5**, already exposes a non-count value, and line
**801**, columns **5–7**, exposes two subnormal doubles and their ratio.
Parsing the displayed subnormals as IEEE-754 binary64 gives bit patterns
`0x1911` and `0x3b81`. These are consistent with non-count memory contents;
the golden alone cannot establish which allocator objects occupied that
memory. Their provenance as out-of-range reads is established by the source.
They are not rounding errors in integer cumulative sums.

Tier1 has 8657 quality bins (8658 lines). Its final two lines are:

```text
26678.1	27245	10552	2.58197	0	1	0
30433.1	27246	10552	2.58207	0	0	nan
```

These happen to agree with zero boundary values. They do not make the source's
invalid reads defined. The two goldens do not capture process heap state or
establish repeatability across upstream executions.

Other pinned behavior in this function: quality thresholds are sorted unique
`double` keys, lower sums exclude the current bin, and upper sums select
`ui+1` in the **one-based** upper cumulative array. This excludes the current
bin; the earlier blocker note's claim that it also excluded the next bin was
an indexing error in that note. The final two golden rows above confirm the
actual indexing. Output uses the default C++ stream precision of six significant
digits, trailing fractional zeros removed, scientific notation for rounded
exponents below −4 or at least 6. The goldens spell zero-over-zero as `nan`.
None of these rules explains away the subnormals.

## Additional source findings (original inspection, now implemented)

* `.TsTv.count` allocates **2 × sample count** elements, without the extra
  element needed for the maximum alternate allele count. A fully called,
  all-homozygous-alternate biallelic SNP writes past the vector at lines
  3195/3197. Output only visits indices below that maximum. The lead-authorized implementation ignores
  that bin and prints precisely the visible rows: both goldens match byte for
  byte. No normalization applies to this report; arbitrary effects of the
  invalid write are not reproduced. Its SNP model table recognizes the
  twelve directed substitutions between A, C, G and T.
* `.FILTER.summary` counts every retained site despite its log message saying
  "bi-allelic loci only". Its substitution model uses REF and the first ALT
  only, sorting their concatenation and recognizing six two-base models. Rows
  sort by decreasing site count, then decreasing FILTER string on ties. Ratios
  use default stream formatting.
* The earlier worker deferred the bcftools audit. That audit is now complete for
  the requested single-file, no-options mode; its rules are recorded below and
  both reports pass byte identity.

## Defined implementation and gate state

`nomics-qc/src/vcftools.rs` mirrors the three `variant_file_output.cpp`
entry points. `vcf_entry::parse_basic_entry` uppercases alleles, missing QUAL is
`-1.0` (`header::str2double`), `set_ALT` drops `.`, and `set_FILTER` sorts the
individual flags before `get_FILTER` joins them. TsTv uses only recognized
biallelic SNP substitutions; counts come from GT, with haploid and partially
missing calls handled by `get_allele_counts`. The upstream count operation
rejects polyploid GTs. FILTER summary counts all sites and the first ALT's
substitution, including on multiallelic sites. Numeric keys compare signed zeros
as equal and retain the first spelling, as C++ maps do.

`nomics-qc/src/bcftools.rs` implements these default-mode rules:

* `bcf_set_variant_type` / `bcf_get_variant_types` classify each ALT, then OR
  legacy masks. Prefix and suffix trimming are case-insensitive; the single-base
  fast path is case-sensitive. SNP, MNP, indel and OTHER site counts overlap on
  mixed records. Breakends and `*` have separate masks and do not increment
  OTHER. `.`, identical alleles, `X`, `<X>`, `<*>`, `<NON_REF>` are reference-like
  under the source's respective branches. Multiallelic SNP counts require the
  complete mask to equal SNP. Padded SNPs use the **first** REF/ALT bases for
  substitution stats, even when those bases are equal (then no substitution).
* `init_iaf` / `bcf_calc_ac`: default AF uses INFO/AC and INFO/AN only. GT is not
  a fallback without sample options. AC=1 and unknown counts use bin 0, separately
  printed as SiS, then merged into bin 1 for AF. With `m=max(101,samples+1)`,
  other AF indices are `1+trunc(float(AC)/AN * (m-2))`, clamped to [0,1] before
  multiplication. The displayed frequency is double `(index-1)/(m-1)` with six
  decimals. Missing integer AC values clamp to AF0 rather than singleton counts.
* `do_snp_stats` counts all recognized SNP ALTs for AF/ST/TSTV; only ALT1
  contributes QUAL and first-ALT TSTV. Ratios divide in float32 before promotion
  for two-decimal printing, and zero TV yields `0.00` here.
* `do_indel_stats` counts one site in QUAL but each indel ALT in AF and IDD.
  Lengths are signed, capped at ±60. Without a reference all repeat-context counts
  go to n/a, and without sample selection IDD genotype counts are 0/VAF `.`.
* QUAL is parsed to double then stored in float32 (`vcf_parse`); `do_snp_stats`
  and `do_indel_stats` form `1+trunc(QUAL*10)` in float32, with missing/negative
  QUAL at bin 0. `dist_init(5)` is exact through 100000, then uses decimal
  levels with steps 10, 100, etc. `dist_get` maps even index 100000 into the next
  level. Labels use `float(0.1*(bin_start-1))` with one decimal; bin 0 prints `.`.
* `do_vcf_stats` collects only scalar integer INFO/DP for default depth stats:
  bins `<0`, 0–500, `>500`; percentages divide by the number of sites **with DP**,
  not the total VCF records. No per-sample sections are printed. Section order:
  ID, SN, TSTV, SiS, AF, QUAL, IDD, ST, DP. Every header and separator is retained.

Inputs are streamed through `nomics_io::bgzf::BgzfReader`. QC does not run or
link any upstream tool. Invalid column counts, polyploid biallelic GTs, invalid
biallelic allele indices, non-finite vcftools quality keys, and qualities outside
bcftools' defined histogram integer range return errors, not normalized output.

### Eight release-CLI comparisons

| Tier | TsTv.count | FILTER.summary | bcftools_stats.txt | TsTv.qual |
|---|---|---|---|---|
| tier0 (1,047 records) | `cmp` identical | `cmp` identical | `cmp` identical | Identical under authorized normalization; 5 exempt fields, all different |
| tier1 (42,358 records) | `cmp` identical | `cmp` identical | `cmp` identical | Identical under authorized normalization; 5 exempt fields, all equal |

Exempt field ledger (one-based file lines, including header):

| Tier | Line | Column | nomics | Golden |
|---|---:|---:|---|---|
| tier0 | 800 | 5 | `0` | `3.17042e-320` |
| tier0 | 800 | 7 | `0` | `3.17042e-320` |
| tier0 | 801 | 5 | `0` | `3.17042e-320` |
| tier0 | 801 | 6 | `0` | `7.5261e-320` |
| tier0 | 801 | 7 | `nan` | `0.421256` |
| tier1 | 8657 | 5 | `0` | `0` |
| tier1 | 8657 | 7 | `0` | `0` |
| tier1 | 8658 | 5 | `0` | `0` |
| tier1 | 8658 | 6 | `0` | `0` |
| tier1 | 8658 | 7 | `nan` | `nan` |

There are **five**, not four, exempt fields on each tier. No other residual
byte difference exists on these fixtures. The comparator returns every exemption,
including equal fields; its tests also reject mutations of non-exempt fields,
headers, row counts and newlines. It never selects exemptions from golden values.

### Validation and coverage

Build the release binary before the integration gate:

```sh
cargo build --release -p nomics
cargo test -p nomics-qc -- --nocapture
cargo clippy -p nomics-qc
cargo clippy -p nomics-qc --all-targets --no-deps -- -D warnings
```

`tests/goldens.rs` executes the release CLI for both tiers and performs all eight
comparisons (six also invoke `cmp`). `NOMICS_RUNS` selects the runs root, default
`~/omics/runs`; missing VCF/reports produce an explicit per-tier skip message.
`NOMICS_QC_BIN` can select a release binary in a different target directory.
Seven unit tests cover source-derived classification tables, decimal histogram
boundaries, stream formatting, the normalization's scope, and inline VCFs with
missing/negative QUAL, zero upper counts, multiple samples, haploid/partially
missing GTs, filter sorting/ties, INFO-vs-GT disagreement, missing AC/AN, mixed
alleles, MNPs, symbolic alleles, breakends, padded SNPs, unknown SNP bases,
reference-only sites, missing/overflow DP and capped indels.

Both goldens contain only one diploid sample, FILTER `.`, numeric QUAL and
INFO/AC/AN; they exercise SNPs, indels, mixed multiallelic records and `*`.
Neither golden exercises multiple samples (including the >100-sample AF width),
haploid/missing calls, MNPs, OTHER symbolic/complex alleles, breakends,
reference-only records, lowercase alleles, missing QUAL, absent AC/AN, negative
or missing DP, or non-dot FILTER sorting/ties. Inline unit cases cover these
branches except the >100-sample AF width. Tests do not substitute for new
source/golden gates for other samples or record distributions.

Validation passes: seven unit tests, one integration test (8/8 release gates),
release build, and clippy with no warnings in the new crate. Plain clippy reports
three pre-existing `nomics-io` warnings; the release build reports the existing
`nomics-hc` deprecated `fetch_update` warning. Those crates are outside scope.

## Architecture: `nan` and `-nan`

vcftools writes each Ts/Tv ratio with `operator<<(double)`, which glibc formats with the
NaN's sign. `0.0/0.0` yields the platform's default NaN: sign bit set on x86 (`-nan`), clear
on arm64 (`nan`). The gate over every golden on the Spark (2026-11-11) passed on all 79 arm64
runs and failed on the four x86 goldens (wgs20-x86-r1, wgs30s-x86 r1 and r2, wgs30-x86-r1) at
exactly these fields, in `.TsTv.count` and `.TsTv.qual`. `ReferenceArch` carries the choice;
`nomics qc vcf-stats --reference-arch` defaults to the host, which is what sarek on the same
machine writes. `.FILTER.summary` and the bcftools report had no such field in any golden.

## Gate over every golden (Spark, 2026-11-11)

`scripts/qc_vcf_gate.sh ~/omics/runs` computes the four reports from each sarek run's own
VCF and compares them with the reports that run published: 83 runs, 83 pass (tier0, tier0nv,
tier1, wgs20, the four x86 goldens, the exomes E1 to E6, the fuzz cases, the genomes G1 to G3,
the Mutect2 decision run's germline calls, and the sarek side of the verify and CI runs).
`.TsTv.count`, `.FILTER.summary` and the bcftools report are byte-identical in all 83.

`.TsTv.qual` needed its exemption in 47 of them, in 4 to 13 fields each, always subnormal
values or ratios made from them where the true count is 0. The stray values are not stable:
the same tier0 input gave `7.5261e-320`, `7.54191e-320` and `7.50239e-320` in the last
row's column 6 across the golden, the verify runs and one CUDA-side run. So the reference
does not reproduce these fields itself, and no byte-exact port could match all of its runs.
