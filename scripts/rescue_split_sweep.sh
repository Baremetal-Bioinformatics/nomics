#!/usr/bin/env bash
# Sweep the CPU/GPU mate-rescue split on 30x chunk 1 (scripts/wgs30_align_gate.sh,
# timing on for the per-batch rescue lines); each config is byte-gated against
# the 30x run's own chunk BAM. Words of a config:
#   tileN / lanesN   NOMICS_RESCUE_GPU_TILE / NOMICS_RESCUE_GPU_LANES
#   cpuN / shareF    NOMICS_RESCUE_CPU_THREADS / NOMICS_RESCUE_CPU_SHARE
#   off              NOMICS_GPU_RESCUE=0 (rescue inside pairing, on the CPU)
#   ext0             NOMICS_GPU_EXTEND=0 (the CPU-extension path)
#   dump             NOMICS_DUMP_RESCUE: the first batch's calls and results, for ksw_bench
#   gpu              no override (the default)
#
#   nohup timeout 6h scripts/rescue_split_sweep.sh [out] &
# Each config takes the resource lock for its own run only, so a box shared
# with another worker interleaves between configs.
set -uo pipefail
LOCK=${NOMICS_LOCK:-$HOME/.cache/uni-rnaseq-resource.lock}
REPO=$(cd "$(dirname "$0")/.." && pwd)
OUT=${1:-$HOME/omics/runs/rescuek-sweep-$(date +%Y%m%d-%H%M)}
mkdir -p "$OUT"
S=$OUT/summary.tsv
echo -e "config\talign_s\tmean_gpu\trescue_gpu_share\tmean_rescue_s\tgates" > "$S"
CONFIGS=${CONFIGS:-"gpu
tile65536
cpu4
cpu8
cpu20
share0.25
share0.5
share1
off"}
while IFS= read -r cfg; do
  env=()
  for w in $cfg; do
    case $w in
      tile*) env+=(NOMICS_RESCUE_GPU_TILE=${w#tile}) ;;
      lanes*) env+=(NOMICS_RESCUE_GPU_LANES=${w#lanes}) ;;
      cpu*) env+=(NOMICS_RESCUE_CPU_THREADS=${w#cpu}) ;;
      share*) env+=(NOMICS_RESCUE_CPU_SHARE=${w#share}) ;;
      off) env+=(NOMICS_GPU_RESCUE=0) ;;
      ext0) env+=(NOMICS_GPU_EXTEND=0) ;;
      dump) env+=(NOMICS_DUMP_RESCUE=$OUT/${cfg// /-}.rescue.bin) ;;
      gpu) ;;
    esac
  done
  name=${cfg// /-}
  d=$OUT/$name
  # shellcheck disable=SC2016
  flock "$LOCK" bash -c '
    nvidia-smi --query-gpu=utilization.gpu --format=csv,noheader,nounits -l 1 > "$1" 2>/dev/null &
    smi=$!
    shift
    "$@"
    kill $smi 2>/dev/null' _ "$OUT/smi-$name.csv" \
    env "${env[@]}" CHUNKS=1 "$REPO/scripts/wgs30_align_gate.sh" "$d" > "$OUT/$name.log" 2>&1
  align=$(grep -oE 'align-only L1 chunk 0001: [0-9.]+' "$d/run.log" | grep -oE '[0-9.]+$')
  gpu=$(awk '{s+=$1; n++} END {if (n) printf "%.0f", s/n}' "$OUT/smi-$name.csv")
  share=$(grep -oE 'gpu rescue: [0-9]+ calls \([0-9]+ on the GPU' "$d/run.log" \
    | awk '{t+=$3; g+=substr($5,2)} END {if (t) printf "%.3f", g/t; else print "-"}')
  tail=$(grep -oE 'gpu rescue [0-9.]+s' "$d/run.log" | awk '{s+=$3; n++} END {if (n) printf "%.2f", s/n; else print "-"}')
  gates="$(grep -c PASS "$d/gates.txt" 2>/dev/null) PASS $(grep -c FAIL "$d/gates.txt" 2>/dev/null) FAIL"
  echo -e "$cfg\t$align\t$gpu\t$share\t$tail\t$gates" | tee -a "$S"
  rm -f "$d"/mapped/*.bam  # gated; keep logs only
done <<< "$CONFIGS"
echo done >> "$S"
