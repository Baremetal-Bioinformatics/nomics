# nomics-quant: Salmon 2.7.0 alignment gate (2026-10-15)

`nomics salmon-quant` implements the deterministic, default alignment quantifier used by nf-core/rnaseq 3.27.0. Fresh code; no Salmon crates, bioinformatics libraries, or compression libraries are dependencies. `statrs = 0.19.0` supplies general-purpose probability distributions and digamma; Rayon schedules fixed inference shards. BAM/BGZF/inflate/deflate come from nomics-io. Numeric inference and result arrays use `nomics_mem::HostArray`, a sealed numeric view backed by `HostBuf` (huge-page advice before first touch); unsafe stays inside nomics-mem.

## Interface and scope

```sh
nomics salmon-quant --geneMap genome_gfp.gtf --threads 6 \
  -t transcriptome.fasta -a WT_REP1.Aligned.toTranscriptome.out.bam \
  --libType=ISR -o WT_REP1
```

The library exposes `nomics_quant::{Options, quantify, run, Quantification, LibraryType}`. `quantify(&Options)` returns unrounded arrays in BAM @SQ order without writing outputs. `run(&Options)` also writes the quantification files. Input BAM must group every fragment's records contiguously. Coordinate-sorted headers are rejected unless `GO:query` is present. FASTA is scanned for reference names without retaining the sequence dataset; BAM is streamed twice with bounded inflate batches. Counts and effective lengths come from BAM, never from golden quantification files.

Supported types: **IU, ISF, ISR, OU, OSF, OSR, MU, MSF, MSR, U, SF, SR, and A**. Format parsing is case-insensitive, as upstream accepts. The original Mac fixture is mixed: **RAP1_UNINDUCED_REP1 and REP2 were run with SR**, as recorded in each golden `cmd_info.json`; the other three were run with ISR. `commands.txt` only provides WT_REP1's alignment command. Passing ISR to the single-end BAMs would reject their placements, just as upstream's explicit strand rule does. The regression test checks the recorded command type before running each sample; it does not silently infer or change flags.

This targets the default AS-score/deterministic/VBEM command. Explicit error-model training (`--errorModel`), online inference, mapping reads, RAD input/output, bias correction, bootstraps/Gibbs, alternative fragment distributions, and additional Salmon flags are unsupported. The CLI rejects unknown flags. None of the six fixture BAMs carries AS tags, so these fixtures exercise the specified missing-AS fallback to score zero; the AS sum and best-score penalty are implemented and covered separately by focused geometry/scoring tests. Validation covers the five yeast/GFP samples and a human transcriptome with 509,650 references using the supplied truncated BAM. The full 20M-pair BAM is not present on this Mac, so no full-20M runtime is claimed.

## Source ledger

Spec: [COMBINE-lab/salmon v2.7.0](https://github.com/COMBINE-lab/salmon/tree/v2.7.0), resolved commit `d3be8b31d31242683c01fc5510c48fbd8bc917fb`. Source was read from a separate checkout under `/private/tmp`, not vendored into this workspace.

| Upstream source / function | Reproduced behavior |
|---|---|
| `salmon-align/src/lib.rs`: `canonical_name`, `record_to_frag`, `pair_records` | Strip a trailing slash plus one character from QNAME; skip unmapped records; pair read 1 with reciprocal same-transcript mate coordinates; prefer matching HI when available, falling back to the first reciprocal match; secondary/supplementary records are retained. |
| `salmon-align/src/bam_rad.rs`: `analyze_group`, `tally_fld_naive`, `WriteWorker` | Sum mate AS scores, absent AS = 0. Learn fragment lengths only from a single placement containing both mates and positive absolute TLEN, before strand filtering. Count opposite/same-strand buckets as integers. |
| `salmon-rad/src/lib.rs`: `RadHit::for_placement`, `frag_map_type::fragment_level`; `salmon-align/src/rad.rs`: `rad_frag_format` | Emulate RAD's lossy geometry without writing a RAD: 16-bit saturated insert/read lengths, fragment-level most-complete status replacing each hit's status, orphan mate-forward false, and opposite-strand pairs classified inward regardless of coordinate order. Last orphan wins equal-rank status ties. |
| `salmon-model/src/fld.rs`: `DiscreteFld::finish`, `FragmentLengthDistribution::{new,add_val,cache,from_log_pmf}` | Majority orientation (opposite wins ties); ascending length-bin integer counts become log masses. Gaussian prior mean 250, SD 25, total mass 1, support 0..1000; five-tap binomial(4, 0.5) kernel excluding length zero. Preserve both cache normalization boundaries of the BAM/RAD round trip. |
| Same file: `conditional_means`, `smoothed_effective_length` | `length - E[L | L <= length]`, clamp lookup at 1000; fall back to raw length when result < 1. |
| `salmon-align/src/rad.rs`: `score_weight_basis`, `process_rad_fragment`; `salmon-core/src/libtype.rs` | Default log weight = `-(best_AS - AS)` plus length-conditioned PMF for pairs, bounded conditional CMF for orphans in paired libraries, zero length term for single-end libraries. Discard incompatible placements (`incompat_prior = 0`). Tally observed library formats once per fragment per distinct format, before filtering. |
| `salmon-core/src/libtype.rs`: `compatible_paired`, `compatible_single`, `summarize_lib_format_counts` | Explicit formats compare pair orientation and the expected strand (U accepts both). Single/orphan compatibility preserves the S/A vs SA/AS distinction and reverses left/right orphan expectations for opposed pairs. `incompatPrior <= 0` drops incompatible hits; positive values add `ln(prior)` to each incompatible hit's log weight before transcript grouping. Fragment judgments and raw per-format tallies stay independent of whether the prior retains a hit. JSON summaries count histogram entries, so their sum can exceed the number of fragments. Strand bias is zero if no observed format agrees. |
| `salmon-model/src/libdetect.rs`: `infer_format_from_counts`; `salmon-align/src/bam_rad.rs`: FLD bake; `salmon-align/src/rad.rs`: format resolution | A uses the complete unique-proper-pair tally collected with the FLD, never a read-prefix budget. Orientation plurality ties prefer inward, then outward; forward fractions <0.3 resolve reverse-stranded, [0.3,0.7) unstranded, >=0.7 forward-stranded. No unique proper pairs leaves A unresolved, with no compatibility judgments and the writer's mapped-count/ratio-1 fallback. The deterministic alignment path keeps its paired length model under A even for single-end BAMs. Raw-string mode decisions are preserved for lowercase inputs too (e.g. lowercase `a` remains unresolved in the RAD inference path). |
| Same `process_rad_fragment`; `salmon-eqclass/src/lib.rs` | Stable ascending transcript grouping, log-sum-exp of placements per transcript, softmax, range-factorization bin `floor(weight * (floor(sqrt(number of transcripts)) + 4))`. Accumulate rounded nonnegative weights at scale 2^40 in u128; finish sorted by transcript IDs then bins. Divide accumulated weights by effective lengths and normalize within each class. |
| `salmon-infer/src/{lib,packed}.rs` | VBEM uniform count start, flat transcript prior 0.01, digamma expectation with fixed transcript-order sum, singleton classes receive full counts. At least 50, at most 10,000 M-steps. Convergence: max relative change `abs(out-in)/out < 0.01`, only `in > 0.01`, positive output. Logical shards: ceil(incidences/4096), capped at min(128, number of classes); nearest cumulative-incidence boundary, lower boundary wins ties. Reduce in ascending shard order. Final sequential masked M-step, threshold 1e-8, using exp(digamma(alpha+0.01)) with truncated transcripts excluded; no global rescaling. |
| `salmon-align/src/lib.rs`: `write_outputs`; `salmon-core/src/genemap.rs` | Transcript rows in @SQ order; integer Length, three fixed decimal places for EffectiveLength/NumReads, six for TPM. TPM = counts/effective length normalized over all transcripts. Gene rows sorted lexically; counts/TPM sums of unrounded estimates, lengths weighted by unrounded TPM, arithmetic mean if no TPM. First GTF mapping wins; unmapped transcript names become their own genes and are listed in `aux_info/genemap_unmatched_txps.json`. |
| `salmon-model/src/dumps.rs`; `salmon-core/src/libtype.rs` | Library counts JSON field order/ratios; per-transcript unique/ambiguous counts; fld = round(10000 * PMF) little-endian i32. The three bias files are upstream's own neutral legacy dumps (one i32 zero for each observed file, one f64 1 for expected), since no bias flags ran. |

## Gate and residual

All six `quant.sf` and all six `quant.genes.sf` files were generated and compared directly with their goldens: **12/12 byte-identical**, including the human IU fixture and a recheck of all five original samples. `lib_format_counts.json`, `aux_info/ambig_info.tsv`, `expected_bias.gz`, `observed_bias.gz`, and `observed_bias_3p.gz` are also byte-identical on every sample. No numeric tolerance or text normalization is applied to these files.

Only `aux_info/fld.gz` has a residual: decompressed histograms match exactly, but compressed bytes differ because nomics-io currently provides zlib-fast levels 1–3, whereas Salmon writes gzip using flate2's zlib-rs backend at level 6. nomics uses its own level-3 deflate with the same deterministic gzip header/trailer. The normalization is explicit in `docs/GATES.md`. Byte offsets below are zero-based; there is no differing data line in these binary files.

| Sample | First compressed byte difference | nomics / golden compressed bytes |
|---|---:|---:|
| WT_REP1 | 11 | 440 / 394 |
| WT_REP2 | 11 | 460 / 422 |
| RAP1_UNINDUCED_REP1 | 16 | 329 / 292 |
| RAP1_UNINDUCED_REP2 | 16 | 329 / 292 |
| RAP1_IAA_30M_REP1 | 11 | 453 / 413 |
| Human IU (and A resolving to IU) | 11 | 516 / 448 |

`meta_info.json` is intentionally not emitted or gated (per task scope). No timestamps, path fields, or estimates are substituted from the reference. The optional fixture test reads expected files only in the test process, after running the library.

## Reproduce validation

```sh
cargo build
cargo test -p nomics-quant -- --nocapture
cargo test -p nomics-mem
cargo clippy
cargo clippy -p nomics-quant --all-targets --no-deps -- -D warnings
cargo build --release
```

The fixture test uses `$HOME/omics-fixtures/rnaseq-test`, or `NOMICS_RNASEQ_FIXTURE` when set. Without the default fixture it prints a skip and succeeds; an explicitly configured missing fixture is an error. With a fixture it runs all five samples, compares every strict artifact above, and checks fld decompressed bytes. Failures report first differing lines, retain outputs in the temporary directory, and fail the test. Separate tests cover reciprocal mate pairing, summed AS, orphan/status loss, range-factorized accumulation under reversed fragment order, multiple inference shards at one versus six threads, and mass conservation.

Validation environment: macOS arm64, available Homebrew Rust/Cargo 1.99.0 (the repository's requested rust-toolchain.toml remains 1.98.1; rustup is not installed here). Cargo's shared cache was read-only in the sandbox, so these commands used `CARGO_HOME=/private/tmp/nomics-quant-cargo`. No lockfile is bypassed.

Workspace `cargo test` was run: the pre-existing `nomics-hc/tests/annotate.rs::tier0_gatk_default_annotations` requires an unavailable external `t0_calls.tsv.gz` oracle and fails at gzip. Tests outside nomics-hc pass. Workspace clippy completes with existing warnings after a small equivalent no-CUDA predicate fix in `nomics-hc/src/caller.rs`: the stub's maximum read length is zero, so only an empty read set satisfies the old `all(1..=0)` check. CUDA behavior is unchanged. Strict clippy on nomics-quant passes.

Initial five-sample release validation (`--threads 6`, Mac, one invocation per sample). Implementation time includes both BAM passes, FASTA validation, inference, GTF parsing, and writing all outputs; subprocess wall time also includes loading the executable. These are elapsed CPU implementation timings, not speedups against Salmon, GATK, or a GPU comparison.

| Sample | libType | quant.sf | quant.genes.sf | Implementation seconds | Process wall seconds |
|---|---|---|---|---:|---:|
| WT_REP1 | ISR | byte-identical | byte-identical | 0.073 | 0.264 |
| WT_REP2 | ISR | byte-identical | byte-identical | 0.036 | 0.039 |
| RAP1_UNINDUCED_REP1 | SR | byte-identical | byte-identical | 0.032 | 0.036 |
| RAP1_UNINDUCED_REP2 | SR | byte-identical | byte-identical | 0.043 | 0.047 |
| RAP1_IAA_30M_REP1 | ISR | byte-identical | byte-identical | 0.037 | 0.041 |

All samples converged in 50 VBEM iterations. WT_REP1's first executable launch had 0.191 s of process startup overhead; this is why its process wall and implementation times are reported separately. No GPU path is involved.

Initial validation results: build and release build passed; all 7 quantifier unit tests, the external five-sample integration gate, and all 3 memory tests passed. Workspace `cargo test --workspace --no-fail-fast` completed with exactly one failed target (the unavailable pre-existing HC annotation oracle); every other target passed or skipped its declared optional fixture. Workspace `cargo clippy` passed with existing warnings; strict `-D warnings` passed for all nomics-quant targets and the nomics-mem library. Direct release CLI byte comparisons confirmed all 10 primary files and all 25 strict auxiliary files identical; all 5 fld histograms identical after decompression. `git diff --check` passed.

The initial implementation's final release rebuild and CLI repeat after the numeric-buffer bounds check again passed all five byte comparisons. Implementation elapsed seconds in sample-table order were 0.074, 0.036, 0.030, 0.042, and 0.036; the stack map shows the resulting 0.030–0.074 s range.


## Human IU / complete library-type follow-up

The new read-only fixture is `$HOME/omics-fixtures/salmon-iu` (override: `NOMICS_SALMON_IU_FIXTURE`). `COMMAND.txt` records Salmon 2.7.0 arm64 `--geneMap annotation.gtf --threads 6 -t genome.transcripts.fa -a sub.bam --libType=IU`. Both `golden-r1` and `golden-r2` were byte-compared for all seven strict files and fld.gz; the reference repeats itself. The BAM has 509,650 reference transcripts and 53,262 fragments (the reference reports 1,490,348 alignment records). Both implementations form 26,906 equivalence classes and converge in 204 VBEM iterations.

The final release CLI was run with IU and A, plus all five original sample commands. All seven strict files match golden bytes for each run, and all fld histograms match after decompression. A resolves to IU on the human data and produces the same requested outputs as the IU golden. No primary differing line exists. The sole requested auxiliary residual is the fld.gz compressed encoding, first difference at zero-based byte 11 on the human sample, 516 vs 448 bytes; its i32 histogram is identical. Metadata remains intentionally skipped; invocation logs/cmd_info and libParams side outputs are outside the original quantitative-output gate.

Input parsing now uses 1 MiB buffered readers. FASTA validation retains a borrowed hash set of BAM reference names and skips sequence lines without copying them or retaining the 922 MB sequence dataset. GTF parsing streams the entire 3.3 GB annotation through one reused line buffer. A hash map stores one mapping per transcript; repeated exon rows do not allocate another key/value or reparse the gene attribute. First GTF mapping wins, last TSV mapping wins, and gene output still uses sorted names and transcript-order reductions. Gzip is detected by content and decoded through nomics-io; a compressed annotation can retain a bare `.gtf` suffix. No whole-file annotation buffer or bioinformatics/compression dependency was added.

`Options` now includes `incompat_prior: f64` (use 0.0 for the pipeline default), also exposed as `--incompatPrior`. `Quantification.expected_format` records the canonical format used by filtering, or the unresolved auto request. Synthetic BAM integration tests cover all 12 formats plus A through the public API and complete library-count JSON, including paired, single, and orphan placements; positive, zero, and negative incompatible priors; a multitranscript posterior that changes with the prior; and A without unique proper pairs. Unit tests check auto-detection threshold boundaries, orientation ties, the full compatibility matrix, zero-agreement strand bias, annotation duplicate semantics, and gzip content detection.

Current verification: `cargo build` and release build pass; `cargo clippy` passes with existing workspace warnings; `cargo clippy -p nomics-quant --all-targets --no-deps -- -D warnings` passes. All 11 unit tests, 3 library-type integration tests, and both optional fixture tests pass (16 tests total, fixtures present). `cargo test` reaches the same existing missing HC annotation-oracle failure documented above; `cargo test --workspace --exclude nomics-hc` passes. The status page was not changed during this follow-up; its SHA-256 remained `5e304171d6ed7e540fe149abd547d1b3a3eb4642fca256c4a8f36b20b92b19dd`.

An isolated final release IU run on this Mac, six threads, including both BAM passes, FASTA validation, VBEM, streaming the GTF, gene aggregation, and all output writing:

| Measurement | Value |
|---|---:|
| Implementation elapsed | 3.509 s |
| Process wall time | 3.519 s |
| User / system CPU | 4.394 / 0.421 s |
| Peak RSS | 387.5 MiB |

The measured input is the supplied BAM prefix, not the complete 20M-pair sample. No full-input performance extrapolation or speedup against Salmon is asserted. Earlier release IU invocations took 3.582–3.624 s; A took 3.528 s, with identical gate results.

## Single-pass input follow-up (2026-10-16)

The BAM is now decoded once. That scan tallies FLD and library-detection counts from the original placements, then retains the RAD geometry for equivalence-class construction after the model is frozen. Storage is `8 * fragments + 16 * placements` bytes, rounded to 2 MiB chunks of `nomics_mem::HostArray<u64>` (backed by `HostBuf`); no full BAM records are retained. A placement preserves transcript ID, position, signed AS sum, saturated RAD length, fragment-level status, and both strands. Storage is released before inference. `inference.rs`, placement likelihood calculations, equivalence-class arithmetic, and output formatting are unchanged.

`SliceBamReader` borrows record payloads from decoded BGZF buffers, using one reusable seam buffer only when a length prefix or payload crosses buffers. `ParallelBgzfReader` implements `BufRead` and hands off the separate inflate buffers without copying them into a concatenated segment. Inflate uses `std::thread::available_parallelism()` (14 workers on this Mac), independently of the inference `--threads` option. AS and HI are extracted together in one auxiliary-tag traversal; first-occurrence, integer-type, missing-tag, and malformed-tag behavior is preserved. Alignment, placement, and pairing scratch vectors are reused between fragments.

Release CLI measurements used the recorded IU inputs, including the full uncompressed FASTA/GTF scans, all output writing, and `--threads 6`. The original binary was saved before editing. Three alternating before/after runs were executed in the foreground, after fixture checks had warmed the file cache. Every timed run passed the four requested byte gates and the decompressed fld gate.

| Run | Wall (s) | User (s) | System (s) | CPU % |
|---|---:|---:|---:|---:|
| Before 1 | 3.21 | 4.00 | 0.32 | 134.6 |
| After 1 | 3.46 | 3.79 | 0.33 | 119.1 |
| Before 2 | 3.42 | 4.28 | 0.36 | 135.7 |
| After 2 | 3.17 | 3.59 | 0.32 | 123.3 |
| Before 3 | 3.24 | 4.04 | 0.32 | 134.6 |
| After 3 | 3.09 | 3.48 | 0.31 | 122.7 |

Median wall time: **3.24 s before, 3.17 s after** (2.2% lower). Median total CPU time: **4.36 s before, 3.91 s after** (10.3% lower). Median CPU utilization: **134.6% before, 122.7% after**, calculated per run as `100 * (user + system) / wall`. Wall-time ranges overlap, so the small prefix does not establish a substantial elapsed-time gain or predict full-20M performance. Each invocation used `/usr/bin/time -l`; it printed real/user/system timing, then exited 1 because the sandbox denies `sysctl kern.clockrate`. Extended `-l` counters are unavailable. The quantifier completed and its artifacts matched; the timing-tool error is not a quantifier failure. The first pre-edit baseline was 3.50 s wall, 3.94 s user, 0.45 s system (125.4% CPU).

Validation: `cargo build --release -p nomics`; `cargo clippy --workspace --all-targets --no-deps` (existing warnings); `cargo clippy -p nomics-quant --all-targets --no-deps -- -D warnings`; formatting checks on the four changed Rust files; and `git diff --check` pass. `cargo test --workspace --release` stops at the existing HC `tier0_gatk_default_annotations` test because its external `t0_calls.tsv.gz` oracle is absent. With that test explicitly excluded via `-- --skip tier0_gatk_default_annotations`, the workspace passes **184 tests, 16 ignored, 1 filtered**. This includes both fixture gates (all six samples), all library-type/prior tests, and new checks for borrowed record slices, every small-buffer seam offset, truncated/invalid record framing, parallel BGZF partial batches and errors, compact placement fields across mapping boundaries, and AS/HI types and duplicate/malformed tags. A stricter `-D warnings` check including nomics-io also encounters its two pre-existing lints (`bam.rs` modulo check and BGZF writer type complexity); no new lint was suppressed. No commit was made, and the status page was not edited.
