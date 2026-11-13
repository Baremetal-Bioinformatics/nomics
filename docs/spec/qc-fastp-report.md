# fastp 1.1.0 reports (sarek 3.10.0)

Source pin: OpenGene/fastp tag `v1.1.0`, commit
`925dd872319d5391c347330ff91d1092e5294fa2`, inspected in `/tmp/fastp-src`.
Source line references below are to `src/` at that commit; the source is the
spec, including the overlap acceptance bug. Scope is the existing paired-end,
no-adapter-trimming, no-correction, default-filter, split-by-lines prep path.

## Statistics and text contract

- `peprocessor.cpp:394-402,428-435,481-485,550-571`: pre statistics and
  duplication use original reads; poly-G trimming precedes insert sampling;
  post statistics include both mates only when both pass. Existing pair filter
  result is `max(result1,result2)` counted twice. No second FASTQ pass is needed.
- `stats.cpp:232-346`: eight base buckets indexed by ASCII `& 7`; Q20 is
  quality byte >=53, Q30 >=63; Q40 sums bytes 73..126 (`:169-171`). Counts and
  quality sums are integers. Cycles stop at first zero coverage (`:148-156`).
  Mean read length is **integer division** (`:225-230`).
- `stats.cpp:174-211`: quality curves divide per-base quality sums by counts,
  falling back to the cycle mean for absent bases. Content divides by all bases;
  GC is the C+G count divided once (not a sum of rounded curves).
- `stats.cpp:268-306,334-346,454-475`: all valid uppercase 5-mers, overlapping,
  A/T/C/G = 0/1/2/3; all 1024 keys emitted in that order, sixteen per line.
  An invalid base resets the rolling calculation. Overrepresentation analysis
  is disabled without `-p` (`main.cpp:136,420-421`); the map is empty but emitted.
- `stats.cpp:909-970`, `peprocessor.cpp:212-232`: merge per-worker integer
  arrays and histograms, then summarize; no floating-point reductions. Thus a
  single accumulator gives identical statistics for this mode.
- `jsonreporter.cpp:24-172`, `filterresult.cpp:234-252`, `stats.cpp:406-495`:
  explicit key order, tabs, spaces, commas, newlines, with **no newline after
  the final brace**. C++ default streams use six significant digits, fixed
  for rounded exponents -4..5, scientific otherwise, trimming fraction zeros.
  Q40 is present only in individual read sections, not the summary.
- `main.cpp:429-433`: command is argv joined with one trailing space per token,
  including executable spelling, and is inserted **without JSON escaping**.
  Callers supply provenance explicitly; comparisons stage the original names.

## Insert sizes: deterministic sample, not scheduling noise

`peprocessor.cpp:481-485,710-723,797-802`: only worker 0 samples, before
filtering, on packs assigned by pack index modulo thread count. Pack size is
256. With 12 workers use input pairs where `(pair_index / 256) % 12 == 0`.
`options.cpp:22`, `main.cpp:123-125`: maximum 512, minimum overlap 30, mismatch
limit 5 and percentage 20. No gap search in this mode.
`overlapanalysis.cpp:16-93` reverse-complements mate 2 (`sequence.cpp:24-57`),
tries forward offsets `0 .. len1-30` (exclusive upper bound), then offsets
`0,-1,...` strictly greater than `-(len2-30)`. Limit is
`min(5, int(overlap_len * 0.2))`. A mismatch above the limit breaks **only
before position 50**. Acceptance also allows excess mismatches if loop index
is **greater than** 50, including a completed long comparison. Reproduce this
bug. Positive offset gives len1+len2-overlap; zero/negative gives overlap.
Unknown and sizes >=512 share bin 512. JSON histogram contains bins 0..511,
unknown separate. Peak excludes unknown and chooses the first maximum
(`peprocessor.cpp:337-348`). The sampled input order is fixed by reader routing;
thread timing changes neither the sample nor integer histogram.

## Duplication: the one exempt number

`duplicate.cpp:8-67,72-91,99-114,136-172`, `duplicate.h:36-37`,
`main.cpp:203-209`: accuracy 1 defaults to two 512 MiB Bloom planes, each
2^32 bits. Hash all original mate-1 then mate-2 sequence bytes, using the
position-dependent prime table (1024 primes starting at 10007, then search
from previous prime +10001), A/T/C/G/other = 7/222/74/31/13. Arithmetic wraps
u64. A pair counts as duplicate iff **both** atomic fetch-or calls see an
already-set bit. Total/duplicate counters are atomic; the two-plane transaction
is not. Two workers visiting the same previously unseen hash pair can each
win a different plane, so neither counts as duplicate. Hash collisions also
make serial order relevant. There is no RNG involved.

A reference run serialized in FASTQ pair order (each pair updates plane 0 then
plane 1 completely before the next pair) has a defined deterministic value.
The implementation uses that order. It is not guaranteed to equal a 12-worker
reference run. The gate therefore exempts **only the numeric token at
`duplication.rate`**, with all other bytes compared and no numeric tolerance
(`docs/GATES.md`, adopted 2026-11-12). Tier1 r1/r2: 0.0977926 vs 0.0977925.

Checked on the Spark, 2026-11-12: two `--thread 1` runs of the pinned fastp
1.1.0 image on the tier1 exome write the same JSON as each other and as
`nomics prep --reports-only --fastp-threads 1`, byte for byte, rate
**0.0977929**. So serial order is fastp's own one-worker result. Over the 60
samplesheet rows of `scripts/qc_reads_gate.sh`, 47 JSON files equal the
12-worker golden outright and 13 differ in this number alone; in all 13 the
serial rate is the higher, by at most 0.000001, as lost duplicates predict.
`scripts/fastp_thread1_check.sh` reruns the one-thread comparison.

Cost, measured on the Spark on the 30x slice (4.9 GB of FASTQ.gz): reports alone
56.8 s wall (89 s of CPU: the statistics are collected in pair order on one
thread), against 39.8 s for chunking without reports on 12 threads. The whole
30x genome took 1,256 s for both FastQC mates and the fastp reports, one after
the other. So writing the reports from `nomics run` needs them off the main
path's critical thread before it is free.

## HTML and log (only after JSON validation)

`htmlreporter.cpp` and `stats.cpp` also write their own literal HTML, curves,
K-mer tables, and plots. Title/footer use local wall-clock time; summary
prints duplication with `to_string` (six fraction digits). `Stats::print`
(`stats.cpp:395-404`), `FilterResult::print`, `peprocessor.cpp:234-271`, and
`main.cpp` write the stderr log, including integer wall-clock elapsed seconds.
The HTML and log gates exempt only timestamp values, the printed duplication
numeric token, and the elapsed integer (log only). Over the same 60 rows: HTML
60 of 60 identical outside those fields; log 14 identical outright (the elapsed
seconds happened to agree) and 46 identical outside them.

## Not confirmed

- Exact worker interleaving that generated any stored reference duplication
  estimate cannot be recovered from reports. Serial order is defined, not
  claimed to be the golden's execution order.
- Empty-input nonfinite spelling is based on arm64 (`nan`); no empty reference run is available.
- Non-sarek modes (single-end, correction, adapters, overrepresentation, gaps,
  dedup, alternative accuracy levels) are outside this implementation's API.

## Implemented interfaces and validation

- `run` and `run_to_dir` remain unchanged. `Report::new` plus `run_with_report`
  collects in the existing pass. Its 1 GiB zero-filled Bloom allocation uses
  `nomics_mem::HostBuf`; all new code is safe Rust and no new external crate
  was added. The file-backed pipeline API does not yet collect reports.
- `nomics prep --json PATH --html PATH --log PATH` writes any requested subset.
  `--report-command TEXT` supplies exact fastp argv provenance, including trailing
  space. Otherwise provenance describes the equivalent fastp operation with
  the supplied paths. A caller controls report destinations; they are not
  created implicitly for existing callers.
- The fastp shim handles `--json/-j` and `--html/-h` natively, retaining original
  argv order/spelling and argv[0]. A JSON request writes a sibling `.log` (or
  appends `.log` if the JSON path has another extension). Its own existing
  stderr diagnostics remain separate from that literal report log. It no longer
  invokes real fastp for reports. Unsupported processing modes still delegate.
- HTML timestamps default to UTC (fastp uses its process-local timezone); the
  API accepts separate literal title/footer timestamps for reproducible output.
  Only the timestamp values are exempt. Elapsed log seconds are the
  difference between integer epoch seconds, as in `main.cpp:435,514-519`.
- HTML: `htmlreporter.cpp:26-55,56-281,344-425,427-630`,
  `filterresult.cpp:357-377`, `stats.cpp:507-529,652-907`. Number units use strict
  `>1000`, fixed six fraction digits; percentage labels truncate their formatted
  strings to five bytes. K-mer shade denominator is `(bases+1)/2048`, not 1024.
  Long reads (>300 cycles) keep positions 1..40 then successive `pos *= 1.05`
  floored to integer, ending at the last cycle; plotted values average the
  intervening positions in increasing order. JSON always retains every cycle.
  Quality histograms print nonzero bins 33..127 as Phred 0..94.
- `head.html` contains only the fixed JS/CSS output literals from
  `htmlreporter.cpp:569-610`, required by the byte contract; no golden-derived
  templates or runtime reference dependency exist.
- Tier0: all three raw reports pass `cmp` when supplied the reference's command,
  two timestamps and elapsed seconds. JSON passes without any exemption,
  duplication rate **0.103521** on both sides. Report on/off chunk vectors and
  filtering counters match exactly. `tests/reports.rs` skips absent fixtures;
  a JSON discrepancy on tier0 fails the test, whichever field it is in; the
  rate exemption is applied only by `scripts/qc_reads_gate.sh`.

## Local verification record (2026-11-12)

Commands (all exit 0 unless explicitly noted):

```sh
# Initial sandboxed fetch could not write shared FETCH_HEAD; retry outside
# the sandbox succeeded, as did the requested branch creation.
git fetch origin && git switch -c qc/fastp-report origin/main
git clone --depth 1 --branch v1.1.0 https://github.com/OpenGene/fastp.git /tmp/fastp-src
cargo test -p nomics-prep -p nomics  # before: 61 tests; after: 68 tests
cargo build --release -p nomics
cargo clippy -p nomics-prep -p nomics --no-deps -- -D warnings
cargo clippy -p nomics-prep -p nomics --all-targets --no-deps -- -D warnings
cargo test --release -p nomics-prep --test reports -- --nocapture
NOMICS_FASTP_FIXTURES=/tmp/nomics-fastp-absent cargo test -p nomics-prep --test reports -- --nocapture
NOMICS_FASTP_FIXTURES=/tmp/nomics-fastp-absent cargo test -p nomics --test fastp_reports -- --nocapture
python3 /tmp/check-fastp-json.py
python3 /tmp/bench-fastp-report.py
cmp /tmp/nomics-fastp-tier0-final/all/HG002_WES-L1.fastp.json "$HOME/omics/runs/tier0-wes-20261008-r1/results/reports/fastp/HG002_WES/HG002_WES-L1.fastp.json"
git diff --check
```

Both absent-fixture tests print an explicit SKIP and pass. The two end-to-end
CLI/shim invocations stage the reference basenames, compare raw JSON with `cmp`,
and compare HTML/log after replacing only timestamps / the elapsed integer.
The API golden test injects the golden's command, two timestamps and elapsed
integer, then runs literal `cmp` on JSON, HTML and log. **No duplication
exemption is exercised on tier0.** If only that field differs, the test proves
all other bytes identical and still fails. Intermediate HTML testing found four whitespace
mismatches in the N curve's color literal; correcting the literal made the
complete comparison pass. Initial clippy findings in new code (checked integer
division and a test-only unnecessary Vec) were fixed. Build/test/clippy retain
one pre-existing dependency warning in `nomics-hc` about deprecated
`AtomicU64::fetch_update`; no unrelated code was changed.

Timing: five interleaved release CLI repetitions per mode, alternating order
`off,json,all` and `all,json,off`, with warm local input cache. Timed the entire
CLI including filtering, compression, output writes, report allocation and
report serialization. Inputs are tier0's 495,559 pairs, output is 24 chunks,
12 fastp workers; compression uses the Mac's available parallelism. No speedup
claim is made. Median wall seconds:

| Reports | Wall seconds | Added vs off |
|---|---:|---:|
| off | 0.763344 | — |
| JSON | 1.520741 | 0.757397 |
| JSON + HTML + log | 1.531299 | 0.767955 |

Raw repetitions (seconds): off `[0.729326,0.763344,0.764898,0.754934,0.777495]`,
JSON `[1.520741,1.513583,1.540386,1.514566,2.424779]`,
all `[1.519196,1.566802,1.531299,1.510693,1.650545]`.
The benchmark runs `target/release/nomics prep --fastq1
$HOME/omics/data/tier0/tier0_1.fastq.gz --fastq2
$HOME/omics/data/tier0/tier0_2.fastq.gz --outdir /tmp/nomics-fastp-tier0-final/MODE
--prefix HG002_WES-L1`, adding `--json .../HG002_WES-L1.fastp.json
--report-command <golden command including trailing space>` in JSON mode and
also `--html .../HG002_WES-L1.fastp.html --log .../HG002_WES-L1.fastp.log`
in all mode. It uses Python `time.perf_counter` around `subprocess.run` and
compares each final JSON directly to the golden. The raw scripts and artifacts
remain under `/tmp`; no input/reference files were changed or deleted.

Tier1 and WGS20 golden shapes were read: both have the same nine top-level
keys, 151/151 cycles, and 512 histogram bins. Their rates are 0.0977926 and
0.0463469 respectively. No wider identity claim is made without inputs.

Deviations/limits: HTML wall-clock values default to UTC; the API can take the
reference's literal local timestamps. Standalone shim logs are inferred from
JSON paths, whereas fastp itself logs to stderr; shim stderr diagnostics remain.
Reports are not yet wired into `nomics run` or its spill path (the requested
prep command/API and fastp shim are wired). The stack map is updated locally;
its private hosted artifact is not republished under this local-only brief.
No push, PR, merge, new delegation, external reference dependency, unsafe code,
bioinformatics crate, compression crate or second input pass was introduced.
