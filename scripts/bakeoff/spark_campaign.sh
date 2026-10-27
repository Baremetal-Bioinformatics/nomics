#!/usr/bin/env bash
# Like-for-like bake-off on the Spark: sarek 3.10.0, nomics and Parabricks 4.7.1,
# raw FASTQ pair -> VCF, same box, one job at a time, page cache evicted first.
#
#   scripts/bakeoff/spark_campaign.sh CAMPAIGN_DIR [phase ...]   (nohup; takes the lock per run)
#   phases: tier1 wgs30s wgs30 (default: all three, in that order); pb = Parabricks only
#   (tier1 raw x2 + fastp, wgs30 raw x2; fastp chunks from $PB_FASTP_FROM/fastp-tier1)
#
# CAMPAIGN_DIR/runs.tsv gets one line per run (tool, tier, replicate, start, end,
# status, run dir); scripts/bakeoff/sampler.py samples the whole box every 5 s
# into CAMPAIGN_DIR/samples.tsv; scripts/bakeoff/summarize.py turns both into
# wall time, CPU-seconds, peak memory above idle and GPU utilization per run.
# After each run its directory moves to the NAS archive (small results are kept
# in CAMPAIGN_DIR/keep/<run>).
set -uo pipefail
CAMP=${1:?campaign dir}; shift
PHASES=${*:-tier1 wgs30s wgs30}
REPO=$(cd "$(dirname "$0")/../.." && pwd)
LOCK=${NOMICS_LOCK:-$HOME/.cache/uni-rnaseq-resource.lock}
ARCHIVE=${BAKEOFF_ARCHIVE:-${NOMICS_ARCHIVE:-$HOME/omics/runs-archive}/bakeoff-20261017}
REF=$HOME/omics/ref/GATK.GRCh38
FASTA=$REF/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta
BWA=$REF/Sequence/BWAIndex/Homo_sapiens_assembly38.fasta.64
B=$REF/Annotation/GATKBundle
WGS_BED=$REF/Annotation/intervals/wgs_calling_regions_noseconds.hg38.bed
PB_IMG=nvcr.io/nvidia/clara/clara-parabricks:4.7.1-1
mkdir -p "$CAMP/keep" "$ARCHIVE"
RUNS=$CAMP/runs.tsv
[ -s "$RUNS" ] || echo -e "tool\ttier\trep\tstart\tend\tstatus\tdir" > "$RUNS"
log() { echo "[$(date '+%F %T')] $*" | tee -a "$CAMP/campaign.log" >&2; }

# One sampler for the whole campaign.
if ! [ -s "$CAMP/sampler.pid" ] || ! kill -0 "$(cat "$CAMP/sampler.pid")" 2>/dev/null; then
  nohup python3 "$REPO/scripts/bakeoff/sampler.py" "$CAMP/samples.tsv" 5 > /dev/null 2>&1 &
  echo $! > "$CAMP/sampler.pid"
fi

tier_inputs() { # tier -> FQ1 FQ2 BED SAMPLE SHEET PARAMS
  case $1 in
    tier1) FQ1=$HOME/omics/data/fastq/SRR14724482_1.fastq.gz; FQ2=$HOME/omics/data/fastq/SRR14724482_2.fastq.gz
           BED=$HOME/omics/data/Agilent_v7.bed; SAMPLE=HG002_WES
           SHEET=$REPO/reference/samplesheet.tier1.csv; PARAMS=$REPO/reference/params.tier1.yaml ;;
    wgs30s) FQ1=$HOME/omics/data/wgs30_slice/wgs30s_1.fastq.gz; FQ2=$HOME/omics/data/wgs30_slice/wgs30s_2.fastq.gz
           BED=$WGS_BED; SAMPLE=HG002_WGS30S
           SHEET=$REPO/reference/samplesheet.wgs30s.csv; PARAMS=$REPO/reference/params.wgs20.yaml ;;
    wgs30) FQ1=$HOME/omics/data/sra_wgs/fastq/SRR14724532_1.fastq.gz; FQ2=$HOME/omics/data/sra_wgs/fastq/SRR14724532_2.fastq.gz
           BED=$WGS_BED; SAMPLE=HG002_WGS30
           SHEET=$REPO/reference/samplesheet.wgs30.csv; PARAMS=$REPO/reference/params.wgs20.yaml ;;
  esac
}

evict() { # drop every input, reference and cache file of ours from the page cache
  python3 - "$@" <<'PY'
import os, sys
for top in sys.argv[1:]:
    for r, _, fs in os.walk(top):
        for f in fs:
            p = os.path.join(r, f)
            if os.path.islink(p):
                continue
            try:
                fd = os.open(p, os.O_RDONLY); os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED); os.close(fd)
            except OSError:
                pass
PY
}

prepare_run() { # evict, settle, record start
  evict "$HOME/omics/data" "$HOME/omics/ref" "$HOME/omics/cache" "$CAMP/pbref"
  sync; sleep 60
}

record() { # tool tier rep start end status dir
  echo -e "$1\t$2\t$3\t$4\t$5\t$6\t$7" >> "$RUNS"
}

archive_run() { # dir: keep small results, move the rest to the NAS
  local d=$1 k=$CAMP/keep/$(basename "$1")
  mkdir -p "$k"
  (cd "$d" && find . -maxdepth 6 -type f \( -name '*.vcf.gz' -o -name '*.tbi' -o -name '*.tsv' -o -name '*.txt' \
     -o -name '*.log' -o -name 'execution_trace*' -o -name '*.html' -o -name 'time.txt' -o -name '*.metrics' \
     -o -name 'exit_status' -o -name '*_at' -o -name 'params.yaml' \) -size -2G -not -path './work/*' -print0 \
     | xargs -0 -I{} cp --parents {} "$k/" 2>/dev/null)
  if cp -a "$d" "$ARCHIVE/" && diff -rq "$d" "$ARCHIVE/$(basename "$d")" > /dev/null; then
    rm -rf "$d"; evict "$ARCHIVE/$(basename "$d")"; log "archived $(basename "$d"); free $(df -h / | awk 'NR==2{print $4}')"
  else
    log "ARCHIVE FAILED for $d (left in place)"
  fi
}

run_sarek() { # tier rep
  tier_inputs "$1"; local d=$HOME/omics/runs/bakeoff-sarek-$1-r$2
  rm -rf "$d"; mkdir -p "$d"; cd "$d" || return
  envsubst < "$SHEET" > samplesheet.csv; envsubst < "$PARAMS" > params.yaml
  export PATH=$HOME/.local/bin:$PATH JAVA_HOME=$HOME/.local/opt/jdk17 NXF_VER=26.04.6
  prepare_run; local s; s=$(date +%s.%N); log "sarek $1 r$2 start"
  /usr/bin/time -v timeout 72h nextflow run "$HOME/omics/ref/sarek-3.10.0" -profile docker,arm64 \
    -c "$REPO/reference/spark.config" -params-file params.yaml \
    -with-trace -with-report report.html -with-timeline timeline.html > nextflow.out 2> time.txt
  local rc=$? e; e=$(date +%s.%N); echo "EXIT=$rc" > exit_status
  record sarek "$1" "$2" "$s" "$e" "$rc" "$d"; log "sarek $1 r$2 exit $rc ($(echo "$e - $s" | bc) s)"
  cd "$CAMP" && archive_run "$d"
}

run_nomics() { # tier rep
  local d=$HOME/omics/runs/bakeoff-nomics-$1-r$2
  rm -rf "$d"; prepare_run; local s; s=$(date +%s.%N); log "nomics $1 r$2 start"
  # full_exome.sh times the run (stages.tsv, written when the run ends), then gates it.
  timeout 6h "$REPO/scripts/full_exome.sh" "$1" "$d" > "$CAMP/nomics-$1-r$2.log" 2>&1
  local rc=$? e
  e=$(stat -c %.9Y "$d/stages.tsv" 2>/dev/null || date +%s.%N)
  record nomics "$1" "$2" "$s" "$e" "$rc" "$d"; log "nomics $1 r$2 exit $rc ($(echo "$e - $s" | bc) s); gates: $(grep -c PASS "$d/gates.txt" 2>/dev/null) PASS lines, $(grep -c FAIL "$d/gates.txt" 2>/dev/null) FAIL"
  archive_run "$d"
}

pb_ref() { # Parabricks wants <fasta>.{amb,ann,bwt,pac,sa} beside the real FASTA
  mkdir -p "$CAMP/pbref"
  for x in fasta fasta.fai dict; do ln -f "${FASTA%.fasta}.$x" "$CAMP/pbref/Homo_sapiens_assembly38.$x"; done
  for x in amb ann bwt pac sa; do ln -f "$BWA.$x" "$CAMP/pbref/Homo_sapiens_assembly38.fasta.$x"; done
}

run_pb() { # tier rep [fastp]
  tier_inputs "$1"; pb_ref; local cfg=${3:-raw} d=$HOME/omics/runs/bakeoff-pb-$1-$cfg-r$2
  rm -rf "$d"; mkdir -p "$d"
  local in1=$FQ1 in2=$FQ2
  if [ "$cfg" = fastp ]; then # sarek's fastp chunks from this campaign's r1 work dir, joined in order
    local w=$CAMP/fastp-$1; in1=$w/fastp_R1.fastq.gz; in2=$w/fastp_R2.fastq.gz
    [ -s "$in1" ] || { log "pb $1 fastp: no chunks at $w, skipping"; return; }
  fi
  local RG='@RG\tID:'$SAMPLE'.L1\tPU:L1\tSM:HG002_'$SAMPLE'\tLB:'$SAMPLE'\tDS:'$FASTA'\tPL:ILLUMINA'
  local R=$CAMP/pbref/Homo_sapiens_assembly38.fasta
  local PB=(docker run --rm --gpus all -u "$(id -u):$(id -g)" -v "$HOME/omics:$HOME/omics" -v "$CAMP:$CAMP" -w "$d" "$PB_IMG" pbrun)
  prepare_run; local s; s=$(date +%s.%N); log "pb $1 $cfg r$2 start"
  local fq=("${PB[@]}" fq2bam --ref "$R" --in-fq "$in1" "$in2" "$RG" --bwa-options="-K 100000000 -Y"
    --knownSites "$B/dbsnp_146.hg38.vcf.gz" --knownSites "$B/beta/Homo_sapiens_assembly38.known_indels.vcf.gz"
    --knownSites "$B/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" --interval-file "$BED"
    --out-bam "$d/$SAMPLE.md.bam" --out-recal-file "$d/recal.txt" --out-duplicate-metrics "$d/dup_metrics.txt" --tmp-dir "$d")
  local t0 rc; t0=$(date +%s.%N)
  timeout 24h "${fq[@]}" > "$d/fq2bam.log" 2>&1; rc=$?
  if [ $rc -ne 0 ]; then # a whole genome may need the low-memory mode on one GB10
    log "pb $1 $cfg r$2 fq2bam exit $rc; retrying with --low-memory"
    t0=$(date +%s.%N); timeout 24h "${fq[@]}" --low-memory > "$d/fq2bam.lowmem.log" 2>&1; rc=$?
    echo "low-memory" > "$d/fq2bam.mode"
  fi
  local t1; t1=$(date +%s.%N)
  if [ $rc -eq 0 ]; then
    timeout 24h "${PB[@]}" haplotypecaller --ref "$R" --in-bam "$d/$SAMPLE.md.bam" --in-recal-file "$d/recal.txt" \
      --interval-file "$BED" --haplotypecaller-options="-pcr-indel-model CONSERVATIVE" \
      --out-variants "$d/$SAMPLE.vcf.gz" --tmp-dir "$d" > "$d/haplotypecaller.log" 2>&1; rc=$?
  fi
  local e; e=$(date +%s.%N)
  printf "fq2bam\t%.1f\nhaplotypecaller\t%.1f\n" "$(echo "$t1 - $t0" | bc)" "$(echo "$e - $t1" | bc)" > "$d/stages.tsv"
  record "pb-$cfg" "$1" "$2" "$s" "$e" "$rc" "$d"; log "pb $1 $cfg r$2 exit $rc ($(echo "$e - $s" | bc) s)"
  archive_run "$d"
}

save_fastp_chunks() { # tier: join sarek r1's 12 fastp chunk pairs (same reads as the goldens) for pb-fastp
  local w=$CAMP/fastp-$1 src=$ARCHIVE/bakeoff-sarek-$1-r1/work
  mkdir -p "$w"
  local dir; dir=$(grep -l 'fastp' "$src"/*/*/.command.sh 2>/dev/null | while read -r f; do x=$(dirname "$f"); [ "$(ls "$x" | grep -c 'fastp.fastq.gz')" -ge 24 ] && echo "$x" && break; done)
  [ -n "$dir" ] || { log "no fastp chunk dir for $1"; return; }
  for m in 1 2; do cat "$dir"/*_R$m.fastp.fastq.gz > "$w/fastp_R$m.fastq.gz"; done
}

locked() { flock "$LOCK" "$@"; }
export -f run_sarek run_nomics run_pb tier_inputs evict prepare_run record archive_run pb_ref log save_fastp_chunks
export CAMP REPO RUNS ARCHIVE REF FASTA BWA B WGS_BED PB_IMG

for phase in $PHASES; do
  log "== phase $phase"
  case $phase in
    tier1|wgs30s)
      for rep in 1 2; do
        locked bash -c "run_nomics $phase $rep"
        locked bash -c "run_pb $phase $rep"
        locked bash -c "run_sarek $phase $rep"
      done
      save_fastp_chunks "$phase"; locked bash -c "run_pb $phase 1 fastp" ;;
    pb) # Parabricks only, same procedure as its bake-off runs (re-run on a new kernel).
      [ -d "$CAMP/fastp-tier1" ] || ln -s "${PB_FASTP_FROM:-$HOME/omics/bakeoff-20261017}/fastp-tier1" "$CAMP/fastp-tier1"
      for rep in 1 2; do locked bash -c "run_pb tier1 $rep"; done
      locked bash -c "run_pb tier1 1 fastp"
      for rep in 1 2; do locked bash -c "run_pb wgs30 $rep"; done ;;
    wgs30)
      for rep in 1 2; do locked bash -c "run_pb wgs30 $rep"; locked bash -c "run_nomics wgs30 $rep"; done
      # sarek's 30x work dir may need ~1 TB: move the previous final-run reference to the NAS first.
      [ -d "$HOME/omics/runs/sprint2-full-wgs30-G" ] && locked bash -c "archive_run $HOME/omics/runs/sprint2-full-wgs30-G"
      locked bash -c "run_sarek wgs30 1" ;;
  esac
done
log "== CAMPAIGN DONE"
kill "$(cat "$CAMP/sampler.pid")" 2>/dev/null
