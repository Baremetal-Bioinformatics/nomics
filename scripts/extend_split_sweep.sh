#!/usr/bin/env bash
# Sweep the CPU/GPU extension split on 30x chunk 1 (scripts/wgs30_align_gate.sh),
# timing off; each config is byte-gated against the 30x run's own chunk BAM.
#
#   nohup flock "${NOMICS_LOCK:-$HOME/.cache/uni-rnaseq-resource.lock}" timeout 4h scripts/extend_split_sweep.sh [out] &
set -uo pipefail
REPO=$(cd "$(dirname "$0")/.." && pwd)
OUT=${1:-$HOME/omics/runs/extend-split-$(date +%Y%m%d-%H%M)}
mkdir -p "$OUT"
S=$OUT/summary.tsv
echo -e "config\talign_s\twall_s\tmean_gpu\tgates" > "$S"
CONFIGS=${CONFIGS:-"cpu
fused
fused cpu1
fused cpu2
fused cpu4
fused cpu2 r16k
fused cpu4 r16k
fused cpu2 share0.3
fused cpu2 share0.5"}
while IFS= read -r cfg; do
  env=(NOMICS_TIMING=0)
  for w in $cfg; do
    case $w in
      cpu) env+=(NOMICS_GPU_EXTEND=0) ;;
      fused) env+=(NOMICS_GPU_EXTEND=fused) ;;
      cpu[0-9]*) env+=(NOMICS_GPU_EXTEND_CPU_LANES=${w#cpu}) ;;
      r16k) env+=(NOMICS_GPU_EXTEND_READS=16384) ;;
      share*) env+=(NOMICS_GPU_EXTEND_CPU_SHARE=${w#share}) ;;
    esac
  done
  name=${cfg// /-}
  d=$OUT/$name
  nvidia-smi --query-gpu=utilization.gpu --format=csv,noheader,nounits -l 1 > "$OUT/smi-$name.csv" 2>/dev/null &
  smi=$!
  t0=$(date +%s.%N)
  env "${env[@]}" CHUNKS=1 "$REPO/scripts/wgs30_align_gate.sh" "$d" > "$OUT/$name.log" 2>&1
  t1=$(date +%s.%N)
  kill $smi 2>/dev/null
  align=$(grep -oE 'align-only L1 chunk 0001: [0-9.]+' "$d/run.log" | grep -oE '[0-9.]+$')
  gpu=$(awk '{s+=$1; n++} END {if (n) printf "%.0f", s/n}' "$OUT/smi-$name.csv")
  gates="$(grep -c PASS "$d/gates.txt" 2>/dev/null) PASS $(grep -c FAIL "$d/gates.txt" 2>/dev/null) FAIL"
  echo -e "$cfg\t$align\t$(echo "$t1 - $t0" | bc)\t$gpu\t$gates" | tee -a "$S"
  rm -f "$d"/mapped/*.bam  # gated; keep logs only
done <<< "$CONFIGS"
echo done >> "$S"
