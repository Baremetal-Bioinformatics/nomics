#!/usr/bin/env bash
# The evidence behind the one number docs/GATES.md exempts in fastp's reports.
# fastp counts duplicates with two atomic updates that are not one transaction, so with sarek's
# 12 workers the duplication rate varies from run to run. With one worker there is no race. This
# runs the pinned fastp 1.1.0 container twice with `--thread 1` on one paired FASTQ and
# `nomics prep --reports-only --fastp-threads 1` on the same reads, and requires all three JSON
# reports to be the same bytes. Needs docker. On the Spark, run it under the resource lock.
#
#   scripts/fastp_thread1_check.sh FQ1 FQ2 [OUT_DIR]      NOMICS=path/to/nomics to pick the binary
set -uo pipefail
REPO=$(cd "$(dirname "$0")/.." && pwd)
bin=${NOMICS:-$REPO/target/release/nomics}
fq1=$(realpath "$1"); fq2=$(realpath "$2"); out=${3:-$(mktemp -d "${TMPDIR:-/tmp}/fastp-thread1.XXXXXX")}
img=community.wave.seqera.io/library/fastp:1.1.0--5650575a9fd255b5
p=SAMPLE-L1   # sarek stages the reads as <sample>-<lane>_R{1,2}.fastq.gz; the names go into the command line

for t in a b; do
  mkdir -p "$out/ref-$t"
  docker run --rm -u "$(id -u):$(id -g)" -v "$out/ref-$t":/w -w /w \
    -v "$fq1":/w/${p}_R1.fastq.gz:ro -v "$fq2":/w/${p}_R2.fastq.gz:ro "$img" \
    fastp --in1 ${p}_R1.fastq.gz --in2 ${p}_R2.fastq.gz --out1 ${p}_R1.fastp.fastq.gz --out2 ${p}_R2.fastp.fastq.gz \
      --json $p.fastp.json --html $p.fastp.html --thread 1 \
      --detect_adapter_for_pe --disable_adapter_trimming --split_by_lines 200000000 --length_required 15 \
      2> "$out/ref-$t/$p.fastp.log" || { echo "fastp run $t failed: $(tail -n 1 "$out/ref-$t/$p.fastp.log")"; exit 1; }
  find "$out/ref-$t" -name '*.fastp.fastq.gz' -delete
done

mkdir -p "$out/ours"
cl=$(python3 -I -c 'import json,sys; sys.stdout.write(json.load(open(sys.argv[1]))["command"])' "$out/ref-a/$p.fastp.json")
(cd "$out/ours" && "$bin" prep --reports-only --fastq1 "$fq1" --fastq2 "$fq2" --outdir . --prefix "$p" --fastp-threads 1 \
   --json "$p.fastp.json" --html "$p.fastp.html" --log "$p.fastp.log" --report-command "$cl") || { echo "nomics prep failed"; exit 1; }

rate() { sed -n 's/.*"rate": *\([-+.0-9eE]*\).*/\1/p' "$1" | head -n 1; }
echo "duplication rate: fastp $(rate "$out/ref-a/$p.fastp.json") and $(rate "$out/ref-b/$p.fastp.json"), nomics $(rate "$out/ours/$p.fastp.json")"
ok=0
cmp -s "$out/ref-a/$p.fastp.json" "$out/ref-b/$p.fastp.json" && echo "fastp --thread 1, two runs: json identical" || { echo "fastp --thread 1, two runs: json DIFFERS"; ok=1; }
cmp -s "$out/ours/$p.fastp.json" "$out/ref-a/$p.fastp.json" && echo "nomics against fastp --thread 1: json identical" || { echo "nomics against fastp --thread 1: json DIFFERS"; ok=1; }
echo "outputs in $out"
exit $ok
