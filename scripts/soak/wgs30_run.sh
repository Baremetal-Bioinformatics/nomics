#!/usr/bin/env bash
# One soak run on the Spark: a full HG002 30x `nomics run` (FASTQ to VCF) with the region checks
# of NOMICS_DIAG=regions on, then all 41 artifacts digested and compared with the wgs30 golden.
#
#   scripts/soak/wgs30_run.sh LABEL ARTIFACT_DIR        (EXTEND0=1: CPU extension, NOMICS_GPU_EXTEND=0)
#                                                       (LOCK=1: host buffers locked in RAM,
#                                                        NOMICS_HOSTBUF_LOCK=1; docs/DIAG.md)
#
# The run is clean when nomics exits 0, the log has no panic, no NOMICS_DIAG FATAL block and no
# GPU FALLBACK, and digest-compare passes 41/41. NOMICS_DIAG=regions keeps the normal HostBuf
# layout (docs/DIAG.md), so a repeat of the 2026-10-24 rid -1 panic stops at its source with the
# read, the batch (replayable with NOMICS_ALIGN_SKIP_BATCHES) and the bad region in the log.
# The run directory is deleted after a clean run; an unclean one is kept.
set -euo pipefail
label=$1 art=$2
REPO=$(cd "$(dirname "$0")/../.." && pwd)
export PATH="$HOME/.cargo/bin:/usr/local/cuda/bin:$PATH"
OMICS=$HOME/omics
CI=${NOMICS_CI_HOME:-$OMICS/ci}
# Peak disk use: the 30x FASTQs are over the 16 GiB streaming threshold, so prep spills them
# uncompressed (about 255 GB for HG002 30x). Each chunk's spill is deleted once the chunk has
# aligned, so the run directory peaks near 290 GB (288 GB measured on 2026-11-11; about 500 GB
# before that change, when r01 on 2026-11-09 ran out of space at 329 GB free).
MIN_FREE_GB=${MIN_FREE_GB:-350}
mkdir -p "$art"
log() { printf '[%s] %s\n' "$(date -u +%H:%M:%S)" "$*"; }

# A lock needs room under RLIMIT_MEMLOCK, and mlock fails (ENOMEM) without it. The GitHub
# runner's systemd service has 8 MiB on the Spark, an ssh login shell unlimited, so stop
# here, before the build and an 80-minute run that could only come back unclean.
if [[ ${LOCK:-} == 1 ]]; then
  ulimit -l unlimited 2>/dev/null || true
  [[ $(ulimit -l) == unlimited ]] || { echo "::error::LOCK=1 needs an unlimited locked-memory limit, this shell has $(ulimit -l) kB: set LimitMEMLOCK=infinity on the service, or start the run from a login shell"; exit 2; }
fi

# Short on disk: a site helper, where there is one, moves finished runs to other storage first
# (each copied and verified before the local copy is swapped for a link).
room=$REPO/scripts/goldens/make_room.sh
[[ -x $room ]] && "$room" "$MIN_FREE_GB"
free_gb=$(df -BG --output=avail "$OMICS" | tail -n 1 | tr -dc 0-9)
log "free on $OMICS: ${free_gb} GB"
(( free_gb >= MIN_FREE_GB )) || { echo "::error::only ${free_gb} GB free on $OMICS, need ${MIN_FREE_GB}"; exit 2; }

target=$CI/target-cuda
cargo build --release --locked -p nomics --features cuda --target-dir "$target" 2>&1 | tail -n 3
commit=$(git -C "$REPO" rev-parse --short HEAD)

mode=default; extra=()
if [[ ${EXTEND0:-} == 1 ]]; then mode=gpu_extend0; extra=(NOMICS_GPU_EXTEND=0); fi
# The migration test (docs/DIAG.md): locked buffers only stay put when the kernel may not
# compact unevictable pages, so the run records that setting and is unclean without it.
compact=$(cat /proc/sys/vm/compact_unevictable_allowed 2>/dev/null || echo unknown)
if [[ ${LOCK:-} == 1 ]]; then mode=${mode}_lock; extra+=(NOMICS_HOSTBUF_LOCK=1); fi
echo "vm.compact_unevictable_allowed=$compact" > "$art/kernel.txt"
out=$OMICS/runs/soak/$label-$mode-$(date -u +%Y%m%d-%H%M)
mkdir -p "$(dirname "$out")"
log "run $label ($mode) at $commit -> $out"
ulimit -c 0   # an 80 GB core would fill the disk; the diag report carries what we need
t0=$(date +%s)
status=0
env NOMICS_DIAG=regions "${extra[@]}" NOMICS_BIN_DIR=$target/release \
  timeout 5h "$REPO/scripts/full_exome.sh" wgs30 "$out" > "$art/full.log" 2>&1 || status=$?
secs=$(( $(date +%s) - t0 ))

cp "$out/stages.tsv" "$out/gates.txt" "$out/run.log" "$art/" 2>/dev/null || true
problems=()
[[ $status == 0 ]] || problems+=("exit status $status")
if [[ -f $out/run.log ]]; then
  grep -m1 -E "^NOMICS_DIAG" "$out/run.log" > "$art/diag_mode.txt" || true
  grep -q "panicked at" "$out/run.log" && problems+=("panic: $(grep -m1 'panicked at' "$out/run.log")")
  grep -q "NOMICS_DIAG FATAL" "$out/run.log" && problems+=("NOMICS_DIAG FATAL")
  grep -qE "GPU FALLBACK" "$out/run.log" && problems+=("GPU FALLBACK")
  if [[ ${LOCK:-} == 1 ]]; then
    grep -q "hostbuf: NOMICS_HOSTBUF_LOCK=1" "$out/run.log" || problems+=("buffers were not locked")
    grep -q "hostbuf: LOCK FAILED" "$out/run.log" && problems+=("$(grep -m1 'hostbuf: LOCK FAILED' "$out/run.log")")
    grep -q "index: mlock failed" "$out/run.log" && problems+=("$(grep -m1 'index: mlock failed' "$out/run.log")")
    [[ $compact == 0 ]] || problems+=("vm.compact_unevictable_allowed=$compact, the lock test needs 0")
  fi
  if (( ${#problems[@]} )); then
    # The report around the failure, and the end of the log.
    grep -n -A80 -E "panicked at|NOMICS_DIAG FATAL|==== NOMICS_DIAG" "$out/run.log" | head -n 400 > "$art/failure_report.txt" || true
    tail -n 200 "$out/run.log" > "$art/run_tail.log"
  fi
else
  problems+=("no run.log")
fi
n_pass=$(grep -c '^PASS' "$out/gates.txt" 2>/dev/null || true)
grep -q '^verdict: PASS' "$out/gates.txt" 2>/dev/null || problems+=("digests: ${n_pass:-0} PASS, verdict not PASS")

verdict=CLEAN; (( ${#problems[@]} )) && verdict=UNCLEAN
printf 'label\tmode\tcommit\tverdict\tdigests_pass\tseconds\tproblems\n%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
  "$label" "$mode" "$commit" "$verdict" "${n_pass:-0}" "$secs" "$(IFS='; '; echo "${problems[*]:-none}")" | tee "$art/outcome.tsv"
column -t "$out/stages.tsv" 2>/dev/null || true
if [[ $verdict == CLEAN ]]; then rm -rf "$out"; else log "kept $out"; fi
[[ $verdict == CLEAN ]]
