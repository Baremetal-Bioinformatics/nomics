#!/usr/bin/env bash
# Profile 30x chunk(s) on the current build, each run byte-gated like
# wgs30_align_gate.sh: a plain timed run, an nsys trace (CUDA, NVTX, OS
# runtime), and a perf CPU profile. Writes summary.txt with align time, GPU
# busy share and per-kernel GPU time (nsys), and the top CPU symbols (perf).
#
#   nohup flock "${NOMICS_LOCK:-$HOME/.cache/uni-rnaseq-resource.lock}" timeout 3h scripts/profile_chunk.sh [out] &   (Spark)
#   CHUNKS=1 (default) picks the chunk(s), as for wgs30_align_gate.sh.
set -uo pipefail
REPO=$(cd "$(dirname "$0")/.." && pwd)
OUT=${1:-$HOME/omics/runs/profile-chunk-$(date +%Y%m%d-%H%M)}
mkdir -p "$OUT/nsys" "$OUT/perf"
S=$OUT/summary.txt
say() { echo "[$(date '+%F %T')] $*" | tee -a "$S"; }
gates() { tr '\n' ' ' < "$1/gates.txt" 2>/dev/null; }
align() { grep -m1 "align-only" "$1/run.log" | sed 's/.*: //'; }

"$REPO/scripts/wgs30_align_gate.sh" "$OUT/plain" > "$OUT/plain.log" 2>&1
say "plain: align $(align "$OUT/plain"); $(gates "$OUT/plain")"

NOMICS_WRAP="nsys profile --trace=cuda,nvtx,osrt --sample=none --cpuctxsw=none --force-overwrite=true -o $OUT/nsys/chunk" \
  "$REPO/scripts/wgs30_align_gate.sh" "$OUT/nsys-run" > "$OUT/nsys-run.log" 2>&1
say "nsys: align $(align "$OUT/nsys-run"); $(gates "$OUT/nsys-run")"
if [ -f "$OUT/nsys/chunk.nsys-rep" ]; then
  nsys stats --force-export=true --report cuda_gpu_kern_sum --format csv -o "$OUT/nsys/kern" "$OUT/nsys/chunk.nsys-rep" > /dev/null 2>&1
  nsys stats --report cuda_gpu_trace --format csv -o "$OUT/nsys/trace" "$OUT/nsys/chunk.nsys-rep" > /dev/null 2>&1
  python3 - "$OUT/nsys" >> "$S" <<'PY'
import csv, glob, sys
d = sys.argv[1]
kern = glob.glob(f"{d}/kern*.csv")
if kern:
    rows = list(csv.DictReader(open(kern[0])))
    tot = sum(float(r["Total Time (ns)"]) for r in rows)
    print(f"  GPU kernel time {tot / 1e9:.1f} s over {len(rows)} kernels:")
    for r in sorted(rows, key=lambda r: -float(r["Total Time (ns)"]))[:8]:
        print(f"    {float(r['Total Time (ns)']) / 1e9:7.1f} s  {int(r['Instances']):>8} launches  {r['Name'][:70]}")
tr = glob.glob(f"{d}/trace*.csv")
if tr:
    iv = []
    for r in csv.DictReader(open(tr[0])):
        if r.get("Name", "").startswith("[CUDA mem"):
            continue
        s, n = float(r["Start (ns)"]), float(r["Duration (ns)"])
        iv.append((s, s + n))
    iv.sort()
    busy, cur_s, cur_e = 0.0, None, None
    for s, e in iv:
        if cur_e is None or s > cur_e:
            if cur_e is not None:
                busy += cur_e - cur_s
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    if cur_e is not None:
        busy += cur_e - cur_s
        span = iv[-1][1] - iv[0][0]
        print(f"  GPU busy (union of kernels) {busy / 1e9:.1f} s of a {span / 1e9:.1f} s span ({100 * busy / span:.0f}%)")
PY
fi

NOMICS_WRAP="perf record -F 199 -g -o $OUT/perf/perf.data --" \
  "$REPO/scripts/wgs30_align_gate.sh" "$OUT/perf-run" > "$OUT/perf-run.log" 2>&1
say "perf: align $(align "$OUT/perf-run"); $(gates "$OUT/perf-run")"
if [ -f "$OUT/perf/perf.data" ]; then
  perf report -i "$OUT/perf/perf.data" --no-children --sort symbol --stdio 2> /dev/null \
    | grep -E '^ +[0-9.]+%' | head -40 > "$OUT/perf/top.txt"
  echo "  top CPU symbols (self):" >> "$S"; head -20 "$OUT/perf/top.txt" | sed 's/^/  /' >> "$S"
fi
say "PROFILE COMPLETE"
