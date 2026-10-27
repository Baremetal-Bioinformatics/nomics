#!/usr/bin/env bash
# Gate and measure the fused GPU extension (NOMICS_GPU_EXTEND=fused) on the Spark.
# Steps, each logged to $OUT/summary.txt:
#   1. build nomics + benches with cuda
#   2. kernel gates: extend_bench EXTEND_GPU=1 on the ksw dumps (refactor guard),
#      fused_bench EXTEND_GPU=1 on the retry fixture and a tier0 chunk dump
#   3. end-to-end gates with fused: tier1 (18/18) and wgs20
#   4. 30x chunk 1: CPU extension, two-pass, fused (timing off), fused (timing on)
#
#   nohup flock "${NOMICS_LOCK:-$HOME/.cache/uni-rnaseq-resource.lock}" timeout 8h scripts/fused_spark_gate.sh [out] &
# Needs $FUSED_DUMP (default ~/omics/bench/fused/tier0-0001.fused.bin, copied from the Mac) and
# $FUSED_WGS_DUMP (default ~/omics/bench/fused/wgs30-0001.fused.bin: every 110th launch of 30x chunk 1).
set -uo pipefail
REPO=$(cd "$(dirname "$0")/.." && pwd)
OUT=${1:-$HOME/omics/runs/fused-gate-$(date +%Y%m%d-%H%M)}
DUMP=${FUSED_DUMP:-$HOME/omics/bench/fused/tier0-0001.fused.bin}
mkdir -p "$OUT"
S=$OUT/summary.txt
say() { echo "[$(date '+%F %T')] $*" | tee -a "$S"; }
export PATH=$HOME/.cargo/bin:/usr/local/cuda/bin:$PATH
cd "$REPO"

say "build"
cargo build --release --features cuda -p nomics > "$OUT/build.log" 2>&1 &&
cargo build --release --features nomics-align/cuda -p nomics-align --example extend_bench --example fused_bench >> "$OUT/build.log" 2>&1 ||
  { say "BUILD FAILED (see build.log)"; tail -30 "$OUT/build.log" >> "$S"; exit 1; }
EX=$REPO/target/release/examples

say "kernel gates"
fail=0
for d in "$HOME"/omics/bench/ksw*/extend2.bin; do
  [ -f "$d" ] || continue
  log=$OUT/extend_bench.$(basename "$(dirname "$d")").log
  EXTEND_GPU=1 "$EX/extend_bench" "$d" 20 > "$log" 2>&1; rc=$?
  # extend_bench exits 0 on GPU mismatches; its summary line carries the count.
  [ $rc -eq 0 ] && grep -q "; 0 mismatches of" "$log" || fail=1
  say "extend_bench $d: rc=$rc $(grep -iE 'gpu|mismatch' "$log" | tail -2 | tr '\n' ' ')"
done
for f in "$REPO/crates/nomics-align/tests/fixtures/fused-retry.bin" "$DUMP" "${FUSED_WGS_DUMP:-$HOME/omics/bench/fused/wgs30-0001.fused.bin}"; do
  [ -f "$f" ] || { say "MISSING $f"; fail=1; continue; }
  log=$OUT/fused_bench.$(basename "$f").log
  EXTEND_GPU=1 "$EX/fused_bench" "$f" 20 > "$log" 2>&1; rc=$?
  [ $rc -eq 0 ] || fail=1
  say "fused_bench $(basename "$f"): rc=$rc $(tail -4 "$log" | tr '\n' ' ')"
done
[ $fail -eq 0 ] || { say "KERNEL GATES FAILED: stopping before end-to-end runs"; exit 1; }

for tier in tier1 wgs20; do
  say "full_exome $tier fused"
  NOMICS_GPU_EXTEND=fused "$REPO/scripts/full_exome.sh" "$tier" "$OUT/$tier-fused" > "$OUT/$tier-fused.log" 2>&1
  say "$tier fused rc=$?: $(grep -cE '^PASS|PASS ' "$OUT/$tier-fused.log") PASS, $(grep -cE 'FAIL' "$OUT/$tier-fused.log") FAIL; $(grep -iE 'total|wall' "$OUT/$tier-fused.log" | tail -1)"
done

# 30x chunk 1 four ways; nvidia-smi at 1 s alongside each.
for mode in cpu twopass fused fused-timing; do
  case $mode in
    cpu) env=(NOMICS_GPU_EXTEND=0 NOMICS_TIMING=0) ;;
    twopass) env=(NOMICS_GPU_EXTEND=1 NOMICS_TIMING=0) ;;
    fused) env=(NOMICS_GPU_EXTEND=fused NOMICS_TIMING=0) ;;
    fused-timing) env=(NOMICS_GPU_EXTEND=fused NOMICS_TIMING=1) ;;
  esac
  d=$OUT/wgs30-chunk1-$mode
  nvidia-smi --query-gpu=timestamp,utilization.gpu,memory.used --format=csv,noheader -l 1 > "$OUT/smi-$mode.csv" 2>/dev/null &
  smi=$!
  t0=$(date +%s.%N)
  env "${env[@]}" CHUNKS=1 "$REPO/scripts/wgs30_align_gate.sh" "$d" > "$OUT/wgs30-$mode.log" 2>&1; rc=$?
  t1=$(date +%s.%N)
  kill $smi 2>/dev/null
  gpu=$(awk -F', ' '{gsub(/ %/,"",$2); s+=$2; n++} END {if (n) printf "%.0f%%", s/n}' "$OUT/smi-$mode.csv")
  say "wgs30 chunk1 $mode rc=$rc wall $(echo "$t1 - $t0" | bc) s, mean GPU $gpu; $(grep -c PASS "$d/gates.txt" 2>/dev/null) PASS $(grep -c FAIL "$d/gates.txt" 2>/dev/null) FAIL; $(grep -iE 'align-only|huge' "$d/run.log" | tail -2 | tr '\n' ' ')"
  [ "$mode" = fused-timing ] && grep -E '\[timing\] gpu extend|extension per' "$d/run.log" | tail -3 >> "$S"
done
say "done"
