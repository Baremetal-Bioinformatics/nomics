#!/usr/bin/env bash
# Re-align chunks of the 30x run's fastp output (`--resume-from align`) and
# compare each mapped chunk BAM byte for byte with the 30x run's own (whose
# downstream artifacts are digest-identical to sarek's). One chunk is 1/12 of
# a 30x genome (~29M pairs, ~10 min): an alignment iteration loop at whole-
# genome scale. Per-batch seed/align/pair timings are logged unless NOMICS_TIMING=0.
#
#   CHUNKS=1,7 scripts/wgs30_align_gate.sh [output-dir]   (Spark, under the lock)
set -euo pipefail
REPO=$(cd "$(dirname "$0")/.." && pwd)
BIN=$REPO/target/release/nomics
REF=$HOME/omics/ref/GATK.GRCh38
P=HG002_WGS30
SRC=${WGS30_RUN:-$HOME/omics/runs/nomics-wgs30-r1}
OUT=${1:-$HOME/omics/runs/wgs30-align-$(date +%Y%m%d-%H%M)}
mkdir -p "$OUT/mapped"
ln -sfn "$SRC/prep" "$OUT/prep"
# Free page cache first so the index and GPU buffers get 2 MiB pages: a
# fragmented run once took 500 s instead of 345 s (NOMICS_NO_EVICT=1 skips).
# Then compact memory where allowed (needs a sudoers rule allowing `tee /proc/sys/vm/compact_memory`;
# skipped otherwise): eviction frees pages but does not defragment them.
[ -n "${NOMICS_NO_EVICT:-}" ] || { python3 "$REPO/scripts/evict_caches.py" "$HOME/omics/data" "$HOME/omics/ref" "$HOME/omics/cache" "$HOME/omics/runs"; sync; echo 1 | sudo -n /usr/bin/tee /proc/sys/vm/compact_memory > /dev/null 2>&1 || true; sleep 10; }
# NOMICS_TIMING=0 runs without timing (any value of the variable turns it on).
if [ "${NOMICS_TIMING:-1}" = 0 ]; then unset NOMICS_TIMING; else export NOMICS_TIMING=1; fi
# NOMICS_WRAP: a command prefix (nsys profile ..., perf record ... --, compute-sanitizer ...).
# shellcheck disable=SC2086
NOMICS_ALIGN_CHUNKS=${CHUNKS:-1} NOMICS_GPU_SEED=1 \
${NOMICS_WRAP:-} "$BIN" run --resume-from align --fastq1 /dev/null --fastq2 /dev/null --reference "$REF" \
  --targets "$REF/Annotation/intervals/wgs_calling_regions_noseconds.hg38.bed" --outdir "$OUT" \
  --dense-sa-cache "$HOME/omics/cache/sa${NOMICS_SA_INTERVAL:-2}.GRCh38.nsa" --dense-sa-interval "${NOMICS_SA_INTERVAL:-2}" \
  --read-group-ds /data/omics/ref/GATK.GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta \
  --sample HG002_HG002_WGS30 --prefix $P --threads 20 > "$OUT/run.log" 2>&1
: > "$OUT/gates.txt"
for bam in "$OUT"/mapped/$P.*.bam; do
  name=$(basename "$bam")
  if cmp -s "$bam" "$SRC/mapped/$name"; then echo "PASS $name" >> "$OUT/gates.txt"
  elif "$BIN" dev gate-bam --ours "$bam" --golden "$SRC/mapped/$name" > "$OUT/gate_$name.log" 2>&1; then
    echo "PASS $name (normalized)" >> "$OUT/gates.txt"
  else echo "FAIL $name (see gate_$name.log)" >> "$OUT/gates.txt"; fi
done
# Timing is only comparable with the index and GPU buffers on 2 MiB pages.
grep -E "index: .* on 2 MiB pages" "$OUT/run.log" | grep -qvE "index: ([0-9.]+) GB, \1 GB on" && echo "TIMING SUSPECT: index not fully on 2 MiB pages ($(grep -m1 'index: ' "$OUT/run.log"))" | tee -a "$OUT/gates.txt"
grep -q "on 2 MiB pages (memory fragmented" "$OUT/run.log" && echo "TIMING SUSPECT: GPU buffers on small pages" | tee -a "$OUT/gates.txt"
grep -qE "GPU FALLBACK|GPU warm-up failed" "$OUT/run.log" && echo "TIMING SUSPECT: $(grep -m1 -oE "GPU FALLBACK:[^(]*|GPU warm-up failed" "$OUT/run.log")" | tee -a "$OUT/gates.txt"
echo "== $OUT"; grep "align-only" "$OUT/run.log" || true; cat "$OUT/gates.txt"
