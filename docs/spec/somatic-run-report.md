# Pair-mode implementation report — 2026-11-14

Base: `426bdc8308cd2a219ed786e5bf4990ff3b98ced4`; branch `somatic/run-pair`.
Implementation commit: `b614836`; gate/test commit: `1e17f34`.
The first command printed that SHA. The checkout was attached to
`t3/tumour-normal-fastq-mode`, not detached; the new branch was created explicitly
from the required SHA. No push or PR. No prohibited library paths or dependencies
changed.

[CLI, output tree, decisions, omissions and resume design](somatic-run-design.md)
are the implementation contract. `somatic_run.rs` owns CSV validation and
orchestration; the only execution change in `pipeline.rs` is an internal boolean
that takes its existing BQSR stop boundary. Single-sample calls always set false.
The standalone chain is unchanged. Pair preprocessing is normal then tumour,
with full function return between them; the chain consumes their published BAMs.

## Local validation

All heavy work used `lockf -k /tmp/nomics-heavy.lock`; no concurrent/background
heavy jobs, no workspace-wide test invocation. Builds used this checkout's
`target/`. `/usr/bin/time -l` required sandbox escalation because the sandbox
blocked `sysctl kern.clockrate`; the first sandboxed test run passed its tests
but the time command exited 1. All reported measurements below are subsequent
successful measurements with that restriction lifted.

Commands were run from the checkout root:

```sh
lockf -k /tmp/nomics-heavy.lock /usr/bin/time -l cargo build --release -p nomics
lockf -k /tmp/nomics-heavy.lock bash /tmp/nomics-pair-fixture.sh
lockf -k /tmp/nomics-heavy.lock bash /tmp/nomics-pair-plumbing/run.sh
lockf -k /tmp/nomics-heavy.lock bash /tmp/nomics-pair-plumbing/chain.sh
lockf -k /tmp/nomics-heavy.lock bash /tmp/nomics-pair-tests.sh
lockf -k /tmp/nomics-heavy.lock python3 /tmp/nomics-pair-plumbing/gate-selftest.py
```

The fixture and pipeline commands are now preserved in
`scripts/somatic_run_plumbing.sh` and `scripts/somatic_run_plumbing_compare.py`.
To repeat exactly that recipe with fresh artifact paths:

```sh
lockf -k /tmp/nomics-heavy.lock bash scripts/somatic_run_plumbing.sh /tmp/nomics-pair-repeat
```

They use installed samtools (`view`, `collate`, `fastq`) on the existing
`~/omics/runs/somatic/m2-chr21/{HCC1395BL_WES,HCC1395_WES}.chr21.bam`, retain the
first 2,000 complete pairs from chr21:10400001-10500000, and gzip paired FASTQs.
All fixtures and outputs are external to git, under `/tmp/nomics-pair-plumbing`.
The three command scripts and their full stderr/time files remain there.
Reference, dense-SA and known-sites cache are the existing `~/omics` copies;
somatic resources are the slice's indexed chr21 VCFs. Six CPU threads, no GPU.

Results: **90 byte comparisons pass**, including prep chunks, mapped chunks,
recal tables, every published chain file and retained per-group chain file.
Four BAM comparisons pass (normal md/recal: 4,005/3,677 records; tumour md/recal:
4,003/3,617 records), with zero mismatches and identical headers excluding `@PG`.
Both metrics comparisons pass excluding only lines 2 and 4, as already specified
in GATES.md. Different command lines/staging paths prevent raw BAM and metrics
byte identity; no new normalization was introduced. Independent single-sample
runs include HC; the pair does not. The fixture produces two somatic records,
zero PASS. This is **orchestration equivalence, not reference qualification**.

Measured `/usr/bin/time -l` results:

- Release build: 33.31 s; peak RSS 1,854,324,736 bytes.
- Pair: 32.36 s; whole-process peak RSS **13,748,879,360 bytes**.
- Standalone normal: 5.94 s; peak RSS 13,742,391,296 bytes.
- Standalone tumour: 5.97 s; peak RSS 13,742,309,376 bytes.
- Standalone BAM chain: 0.35 s; peak RSS 107,986,944 bytes.

The pair log recorded 7.886 s for normal preprocessing/publication and 24.038 s
for tumour, then 0.363 s for the chain. Cumulative high-water RSS was
13,738,885,120 bytes after normal and 13,748,879,360 after tumour. The peak did not
double: the first sample's owned buffers were dropped before the second started.
This is not a measurement of current RSS after freeing, an allocator-release
guarantee, or a large-pair memory bound. Timings on the shared Mac are not a
performance comparison; the second index load was slower in the pair run.

Targeted test script (each cargo invocation was timed with `/usr/bin/time -l`
inside the outer lock):

```sh
for filter in somatic_run::tests somatic_chain::tests tests::pair_cli_constraints \
  tests::run_requires_sample_and_prefix pipeline::identity_tests; do
  /usr/bin/time -l cargo test -p nomics --bin nomics "$filter"
done
NOMICS_SOMATIC_INTERVALS="$HOME/omics/runs/somatic/S1-r1/golden/intervals" \
  /usr/bin/time -l cargo test -p nomics --bin nomics \
  scatter::tests::somatic_reference_groups -- --ignored
```

All **19 tests** passed (2, 7, 1, 1, 7, 1). Peak RSS per invocation, including
cargo/compiler: 546,422,784; 366,428,160; 367,230,976; 367,378,432; 367,476,736;
367,198,208 bytes. The parser table covers missing partners, extra tumours,
multiple patients, invalid/conflicting statuses, duplicate lanes, unsafe
patient/sample/lane names, sex mismatch/invalid sex, malformed CSV/headers,
empty paths and same-path mates. CLI tests cover mutually exclusive flags,
required resources, pair-only resources/orders and the unchanged mandatory
single-sample identity.

Scatter: concatenating the copied reference groups in chr1 then chr17 order and
running the actual Rust scatter implementation reproduced **296,453 records**,
`chr1_11981-12351.bed` and `chr17_38877931-38878296.bed`, byte for byte. The original
full target BED is absent locally; equivalence from that input remains for the
gate machine. The self-test of `CHAIN_PRODUCED` compared reference output copies
to the goldens: 16/16 pass; removing one produced stats file gave a failure.
Those 16 are a comparator-harness self-test, **not nomics reference results**.

`cargo check -q -p nomics`, shell syntax checks, Python compilation and
`git diff --check` pass. Clippy completes; it reports the pre-existing HC
`fetch_update` deprecation and a large CLI enum warning (the new run flags enlarge
that one startup allocation). The new module's manual-flatten suggestion was fixed.

The final work layout was also rerun through the checked-in script:

```sh
lockf -k /tmp/nomics-heavy.lock /usr/bin/time -l cargo build --release -p nomics
lockf -k /tmp/nomics-heavy.lock bash scripts/somatic_run_plumbing.sh /tmp/nomics-pair-final
```

Again 90 byte comparisons, four BAM gates and two metrics comparisons passed.
This verifies the separate `work/preprocessing/S` and `work/chain` layout;
those paths prevent a sample named `chain` from aliasing chain work by construction. Final pair: 25.66 s, peak RSS **14,208,516,096 bytes**;
normal boundary high-water 13,733,773,312 bytes, tumour boundary 14,208,516,096.
Final independent normal/tumour runs: 5.66/5.49 s and
13,738,033,152/13,738,786,816 bytes; standalone chain: 0.35 s and 108,019,712 bytes.
Final artifacts: `/tmp/nomics-pair-final`; full comparison log:
`/tmp/nomics-pair-final.log`. The work naming fix changes only pair staging paths.
Fixture-policy checks also passed: absent fixtures return 77 by default and 1
with `NOMICS_REQUIRE_FIXTURES=1`, without creating the output.

## Mac tier0 germline regression

Run **once**, with the pinned DS (the subsequent pair-only work-directory naming
fix does not enter this single-sample path):

```sh
lockf -k /tmp/nomics-heavy.lock env NOMICS_BUILD_DIR="$PWD/target" \
  NOMICS_RG_DS="<the gate machine's absolute path to Homo_sapiens_assembly38.fasta, as in the golden read groups>" \
  NOMICS_THREADS=6 bash scripts/mac_tier0.sh cpu /tmp/nomics-pair-tier0
```

**PASS, all 41 digests**, including 1,047 HC records. `/usr/bin/time -l`:
16.54 s; peak RSS **14,929,084,416 bytes**. Full evidence:
`/tmp/nomics-pair-tier0/{time.txt,stages.tsv,digest-compare.txt,run.log}`.
This checks the single-sample path against its existing golden, not the new pair
against the somatic reference.

## Lead's gate-machine command

Build the branch on the gate machine first, then run under its resource lock:

```sh
flock "$NOMICS_LOCK" env \
  NOMICS="$PWD/target/release/nomics" NOMICS_REQUIRE_FIXTURES=1 \
  NOMICS_THREADS=20 TMPDIR=/tmp \
  bash scripts/somatic_run_gate.sh \
  "$HOME/omics/runs/somatic/S1-r1" \
  "$HOME/omics/runs/somatic/nomics-run-pair-426bdc8"
```

Use a fresh output path. The script prints all preprocessing stage timings,
chain/total timing and whole-process GNU-time peak RSS. It must finish with both
samples' complete preprocessing comparisons and 16/16 chain comparisons; missing
comparisons fail. It reads FASTQs from the original samplesheet and orientation
order from the recorded golden command; it derives each sample's separate BQSR
orders and merge seed using `scatter_order.py` (fixed to select that sample's task).
Estimate **60–90 minutes**, not measured: the prior BAM-only HCC1395 chain took
41 min 28 s; preprocessing, indexing and comparison add work. Allow more on a
busy host or without accelerated alignment. No ssh/reference run was attempted.

## Remaining qualifications

- Full FASTQ-to-reference preprocessing/chain gate and full target-BED comparison
  are pending; leave results marked pending in the design note.
- Multilane CSV is unit-tested; a real multilane pair and the streamed whole-genome
  pair path were not run on this Mac. Large-pair peak memory is unmeasured.
- Pair resume, CRAM mode, full QC/report coverage and omitted indexes are explicitly
  unsupported; see the design's complete output/omission list.
- The full sarek publication inventory is not copied locally; exact per-file
  omission audit needs the gate machine. No claim of complete sarek publication.
