# BWA mismatch penalty implementation report

Base: `7c0d33a1cbeb1047cd1a2cef75a4971866750106` (printed before any other
command); branch `align/mismatch-penalty`. Local Mac verification, 2026-11-15.
Implementation commit: `6c7690f`.
No full tumour alignment, remote execution, push, or PR.

The complete file/line table of direct `opt->b` uses, indirect matrix callers,
and former constants is in [the BWA source audit](../spec/bwa.md#mismatch-penalty-sarek-tumour-alignment).
BWA source was absent from `~/omics/src` and `~/omics/oracle`, so tag v0.7.19
was cloned to `/tmp/nomics-bwa-v0.7.19` and audited at
`b92993c1161e73167181558856567ef2f367e3f0`.
Sarek's upstream tag places its tumour rule at lines 62–63 (the supplied
finding named 67–68); the rule itself is identical.

## Implementation and decisions

- One validated `Options` carries B and its derived matrix in the loaded
  `Index` before workers share it. No global/env/thread-local scoring state,
  no index format change. Only 3/4 are accepted; other scoring flags remain
  unsupported. B alone leaves A=1, gap 6/1, clipping 5/5, T=30, U=17 and
  zdrop=100 unchanged.
- `nomics-align`: selected matrix in extension, patch/global scoring, CIGAR
  reconstruction and rescue; selected B in MAPQ identity and both secondary
  windows. Long-read seed re-scoring remains explicitly refused as before.
- `nomics-simd`: the existing NEON API already receives the parameter through
  the query-profile scores and byte shift. No redundant scalar B argument or
  unsafe code was added. Portable and NEON kernels are exercised separately
  with `nomics-align/portable`; B=4 remains the normal default dispatch.
- GPU kernels are unchanged. B=3 disables fused/two-pass extension and GPU
  mate rescue before planning, regardless of environment overrides. Additional
  assertions protect their entry points. Warm-up skips scoring kernels too.
  GPU seeding remains available. No CUDA or Metal scoring exactness is claimed.
- CLI: `nomics run` and `nomics dev slice` accept
  `--bwa-mismatch-penalty 3|4`. Pair mode automatically selects normal=4,
  tumour=3, logs each selection, and refuses an explicit global override.
  This tree has no `nomics align` subcommand. The two existing SAM/BAM oracle
  examples accept an optional final penalty argument.
- The bwa shim accepts separate or attached `-B`, validates before alignment,
  and preserves the complete argv in its bwa `@PG CL`, including B's position.
  Unsupported B follows the shim's existing refusal/delegation contract; the
  nomics aligner never runs it. Native pipeline headers use nomics provenance,
  adding `--bwa-mismatch-penalty 3` only for B=3 to preserve B=4 headers.
- `nomics verify` passes B to the nomics run, contig slicer, and reference
  BWAMEM1_MEM overlay, including its optional stock comparison. Existing
  `--against` manifests are B=4 and explicitly refuse B=3. Nextflow execution
  of this new overlay was not run on this Mac; only its construction was tested.
- No new dependencies and no changes to somatic, HC, BQSR, prep, or IO crates.
  No new test-only skips: missing B3 fixtures fail the script. Existing ignored
  large/external-fixture tests retain their existing status.
- Stack-map ledger updated locally. Its external published artifact was not
  changed during this local implementation task.

## Commands and results

All heavy commands ran serially under `lockf -k /tmp/nomics-heavy.lock`.
Each build/test/alignment/gate used `/usr/bin/time -l`; build artifacts are
in this checkout's `target/`. Logs: `target/b3-validation/`.
RSS below is maximum resident set size in decimal GB, not a speed comparison;
index paging on the shared Mac affects RSS.

| Command / path | Result | Peak RSS GB |
|---|---|---:|
| `cargo test --release -p nomics-align -p nomics-simd` | 54 passed, 6 existing ignored; no failures | 0.814 |
| Same, `--features nomics-align/portable` | 54 passed, 6 existing ignored; no failures | 0.880 |
| Fixture build, NEON | release nomics + align_chunk built | 1.821 |
| Fixture build, portable | release nomics + align_chunk built | 1.828 |
| NEON chunk 0001 | 81,773 literal SAM records; BAM zero mismatches | 8.307 |
| NEON chunk 0002 | 81,821 literal SAM records; BAM zero mismatches | 7.741 |
| Portable chunk 0001 | 81,773 literal SAM records; BAM zero mismatches | 6.927 |
| Portable chunk 0002 | 81,821 literal SAM records; BAM zero mismatches | 8.021 |
| `scripts/mac_tier0.sh cpu` | B=4 **41/41 PASS**, 1,047 VCF records | 14.894 |

Reproduce the B3 fixture build, alignments, literal SAM comparisons and existing
`nomics dev gate-bam` comparisons with:

```
lockf -k /tmp/nomics-heavy.lock /usr/bin/time -l scripts/mac_bwa_b3_gate.sh
```

Each fixture alignment forces `NOMICS_GPU_EXTEND=fused` and its log confirms
CPU routing. The script uses the fixture README's read group (including DS),
local hg38 index and existing interval-4 dense SA. SAM headers alone are omitted
for literal emission-order comparison; BAM uses the existing normalization.
No differences were masked or added to the normalization.

B4 command:

```
lockf -k /tmp/nomics-heavy.lock /usr/bin/time -l env \
  NOMICS_BUILD_DIR="$PWD/target" NOMICS_THREADS=8 \
  NOMICS_RG_DS="<the gate machine's absolute path to Homo_sapiens_assembly38.fasta, as in the golden read groups>" \
  bash scripts/mac_tier0.sh cpu "$PWD/target/b3-validation/tier0-b4"
```

Tests include source-derived matrix/N scores, secondary window 7, MAPQ
identity at B=3/4, byte shift/saturation/padding, secondary exclusion radius,
rescue byte/word threshold and invariant clipping/T/U values. C differential
suites retain the original cases and add B3 variants for extension, global
alignment and rescue (including byte/word/reverse/overflow paths).

Targeted CLI tests used `cargo test --release -p nomics --bin nomics FILTER`
under the same lock and timer. Filters and peak RSS: `mismatch_penalty_cli`
(1.746 GB), `pair_cli_constraints` (1.747 GB),
`sample_scoring_follows_status` (1.789 GB), and
`mismatch_option_separate_attached_and_rejected` (1.801 GB).
Each passed one test. `reference_and_nomics_select_same_penalty` also
passed one test (1.803 GB peak RSS).
`cargo check -p nomics --all-targets` and `cargo check -p nomics --features metal`
passed; only the existing HC `fetch_update` deprecation warning remained.
Metal required compiler-cache sandbox access; no Metal runtime test was run.
Initial validation harness failures (a missing test-constructor field and
macOS Bash empty-array expansion) were fixed before the successful results.
The initial sandboxed timer could not read `kern.clockrate`; reported RSS
comes from the subsequent permitted `/usr/bin/time -l` runs.

## Lead's gate-machine work

1. Build this branch with `cargo build --release -p nomics --features cuda`
   under the gate machine's heavy-job lock. CUDA cannot be compiled here.
2. Rerun the full matched-pair preprocessing/caller gate. Normal should log
   status 0 / penalty 4; tumour status 1 / penalty 3. Gate all 12 tumour mapped
   BAMs and downstream boundaries. Also rerun the tumour alone with
   `--bwa-mismatch-penalty 3` to confirm the original independent failure.
3. For B=3, CUDA uses GPU seeding when enabled, but portable CPU seed extension
   and CPU mate rescue (on ARM, NEON rescue rows). It never enters CUDA scoring
   kernels, even if `NOMICS_GPU_EXTEND=fused` or `1` is forced. Confirm this exact
   log line in the tumour run:

   ```
   BWA -B 3: seed extension=CPU, mate rescue=CPU (CUDA/Metal scoring disabled); GPU seeding unaffected
   ```

   With `NOMICS_TIMING=1`, tumour extension statistics should have zero fused
   jobs/precomputed lookups. The dispatch guards and assertions, rather than
   GPU fallback counters, enforce this policy. B=4 keeps its previous routes.
4. Before enabling B3 GPU scoring in a future change, parameterize and gate
   the actual CUDA/Metal kernels; this branch intentionally does not do so.
   Run a B3 `nomics verify` reference comparison on its supported Linux host
   to validate Nextflow overlay integration as well.

## Gate machine, 2026-11-15 (added by the lead)

Tree `a6c33a0`, CUDA build. At the default penalty: tier0nv 31 of 31, tier0 18 of 18, tier1
18 of 18, wgs20 18 of 18, no GPU fallback. At penalty 3: the HCC1395 tumour sample's fastp
chunks, 12 mapped BAMs, md.bam, metrics, recal table and recal.bam are all identical to
sarek's (about 101.8M read pairs), with the run log line `BWA -B 3: seed extension=CPU, mate
rescue=CPU`. Tumour alignment took 319.5 s on that route; the same reads at penalty 4 with
the GPU scoring kernels took 298.3 s (and did not match). The Nextflow overlay that
`nomics verify` builds for penalty 3 has not been run.

## GPU scoring at penalty 3, 2026-11-15 (added by the lead)

Steps 3 and 4 above describe the first implementation. They no longer hold: the scoring
kernels turned out to take the matrix as a launch parameter already (`mat` in the fused
parameters, the two-pass inputs and the rescue inputs; `crates/nomics-gpu` tests them with
an arbitrary matrix), and only the launches in `nomics-align` pinned the B=4 matrix. The
launches now pass the run's matrix and the B=3 routing guard is gone. No kernel changed.

Mac (`scripts/mac_bwa_b3_gate.sh`): both fixture chunks identical to BWA on NEON, portable
and Metal; the Metal log shows fused extension and mate rescue on the GPU, no fallback.
Metal has no two-pass extension.

Gate machine, tree `42fb5d3`, CUDA build, the HCC1395 tumour exome alone to BQSR
(`--bwa-mismatch-penalty 3`, 20 threads, GPU seeding, caches evicted before each run),
each run gated against sarek at the 17 preprocessing comparisons:

| Route | Alignment | Run | Gates |
|---|---|---|---|
| fused extension + GPU rescue | 302.7 s | 445 s | 17 of 17 |
| CPU extension + CPU rescue | 351.3 s | 491 s | 17 of 17 |
| fused extension + GPU rescue (repeat) | 294.5 s | 433 s | 17 of 17 |
| two-pass extension + GPU rescue | 311.3 s | 451 s | 17 of 17 |

3,336 fused launches and 264 rescue launches per fused run, 6,672 two-pass launches, no
fallback in any run. One session, one run per cell except the repeat; the CPU figure is
higher than the 319.5 s measured inside pair mode the night before, so compare within
this table only. Germline at the default penalty on the same build: tier0nv 31 of 31,
tier1 18 of 18. Not run: the whole pair through `nomics run --samplesheet` on this tree
(pair mode calls the same alignment code with penalty 3), tier0, wgs20, 30×.
