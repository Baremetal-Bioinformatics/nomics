# shellcheck shell=bash
# Shared by scripts/goldens/exome_golden.sh and genome_golden.sh (sourced, not run).
log() { printf '[%s] %s\n' "$(date -u +%H:%M:%S)" "$*"; }

fetch() {  # fetch URL DEST MD5
  # Resumes a partial download; a transfer slower than 200 kB/s for 5 min is dropped and resumed,
  # so a stalled connection can't hang the job. Gives up after 30 tries.
  local url=$1 dest=$2 want=$3 try t0 t1 size
  if [[ -f $dest ]] && [[ $(md5sum "$dest" | cut -c1-32) == "$want" ]]; then log "have $(basename "$dest")"; return; fi
  log "fetch $url"
  # A partial file from another URL (an earlier source) can't be resumed.
  if [[ -f $dest.part && "$(cat "$dest.part.url" 2>/dev/null)" != "$url" ]]; then rm -f "$dest.part"; fi
  printf '%s\n' "$url" > "$dest.part.url"
  t0=$(date +%s)
  for try in $(seq 1 30); do
    if curl -fsSL --connect-timeout 60 --speed-limit 200000 --speed-time 300 -C - -o "$dest.part" "$url"; then break; fi
    log "  try $try stopped at $(stat -c%s "$dest.part" 2>/dev/null || echo 0) bytes; resuming"
    [[ $try == 30 ]] && { echo "::error::could not fetch $url"; exit 2; }
    sleep 20
  done
  t1=$(date +%s); size=$(stat -c%s "$dest.part")
  log "  $size bytes in $((t1 - t0)) s ($(( size / (t1 - t0 + 1) / 1000000 )) MB/s)"
  local got; got=$(md5sum "$dest.part" | cut -c1-32)
  [[ $got == "$want" ]] || { echo "::error::md5 mismatch for $url: $got, want $want"; rm -f "$dest.part" "$dest.part.url"; exit 2; }
  mv "$dest.part" "$dest"; rm -f "$dest.part.url"
}

free_gb() { df -BG --output=avail "$1" | tail -n 1 | tr -dc 0-9; }

# The fastp task directory of a sarek run (one lane): the task that wrote at least two
# *.fastp.fastq.gz chunks and is not the BWA index step.
fastp_dir() {  # fastp_dir RUN_DIR
  local f d
  for f in $(grep -l 'fastp' "$1"/work/*/*/.command.sh); do
    d=$(dirname "$f")
    [ "$(ls "$d" | grep -c fastp.fastq.gz)" -ge 2 ] && ! ls "$d" | grep -q BWAIndex || continue
    echo "$d"; return
  done
  return 1
}
