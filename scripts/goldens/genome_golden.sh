#!/usr/bin/env bash
# One sample-matrix genome on the Spark: fetch, two sarek 3.10.0 reference runs (r1, r2), and a
# full nomics run gated by digest against each of them.
#
#   scripts/goldens/genome_golden.sh ID ARTIFACT_DIR     (ID is a row of reference/matrix/genomes.tsv)
#
# A whole-genome sarek run needs about 5x its FASTQ size on disk, so nothing is kept long:
#   1. sarek r1, then its digest manifest and scatter order (scripts/scatter_order.py) are written
#      and its BAMs and work tree are deleted;
#   2. nomics with r1's scatter order, digest-compared against r1 (all 41 artifacts), then deleted;
#   3. the same for r2, and nomics again with r2's scatter order against r2.
# Each digest manifest, scatter order and gate result is saved under ~/omics/runs/matrix, so a
# stopped job resumes at the step it was on.
#
# r1 against r2: in WGS mode three things follow how a run happened (the BQSR gather order, the
# recal merge order and samtools merge's time-seeded ID suffixes; docs/REFERENCE-VARIABILITY.md),
# so the recal table, recal BAM and VCF may differ between two sarek runs. Everything up to
# MarkDuplicates must be identical. nomics reproducing each run exactly from that run's scatter
# record shows the later differences are only the recorded scatter. Exit 0 only when all hold.
set -euo pipefail

id=$1 art=$2
REPO=$(cd "$(dirname "$0")/../.." && pwd)
export PATH="$HOME/.cargo/bin:/usr/local/cuda/bin:$PATH"
OMICS=$HOME/omics
DATA=$OMICS/data/matrix/$id
RUNS=$OMICS/runs/matrix
CI=${NOMICS_CI_HOME:-$OMICS/ci}
mkdir -p "$DATA" "$RUNS" "$art"
. "$REPO/scripts/goldens/lib.sh"

row=$(awk -F'\t' -v id="$id" '!/^#/ && $1 == id' "$REPO/reference/matrix/genomes.tsv")
[[ -n $row ]] || { echo "unknown genome: $id" >&2; exit 2; }
IFS=$'\t' read -r _ patient P sex r1 r2 md5_1 md5_2 <<<"$row"
SM=${patient}_${P}
log "$id: $patient $P ($sex), $r1"

# ---- where the disk goes (read only), and our own exome outputs that are no longer needed ----
log "disk survey"
df -h "$OMICS" | sed 's/^/  /'
du -x -d1 -BG "$OMICS" 2>/dev/null | sort -n | tail -n 12 | sed 's/^/  omics: /' || true
du -x -d1 -BG "$OMICS/runs" 2>/dev/null | sort -n | tail -n 15 | sed 's/^/  runs: /' || true
# A matrix exome's nomics output is only our side of a gate that already passed: its BAMs go.
# Its gates, timings and VCF stay. sarek's r1 results (the goldens) are kept.
for d in "$RUNS"/E*-nomics; do
  [[ -s $d/gates.txt ]] && ! grep -q '^FAIL' "$d/gates.txt" || continue
  n=$(find "$d" -maxdepth 1 \( -name '*.bam' -o -name '*.bai' \) | wc -l)
  (( n > 0 )) || continue
  find "$d" -maxdepth 1 \( -name '*.bam' -o -name '*.bai' \) -delete
  log "dropped $n BAM files of $(basename "$d") (its gates passed); $(free_gb "$OMICS") GB free"
done

# ---- inputs, verified by md5; check the disk first ----
FQ1=$DATA/${P}_1.fastq.gz FQ2=$DATA/${P}_2.fastq.gz
fq_bytes=0 missing_bytes=0
for i in 1 2; do
  u=$r1 f=$FQ1; [[ $i == 2 ]] && u=$r2 f=$FQ2
  b=$(curl -fsSIL --retry 3 "$u" | tr -d '\r' | awk 'tolower($1) == "content-length:" {n = $2} END {print n + 0}')
  (( b > 0 )) || { echo "::error::could not read the size of $u"; exit 2; }
  fq_bytes=$((fq_bytes + b)); [[ -f $f ]] || missing_bytes=$((missing_bytes + b))
done
fq_gb=$(( (fq_bytes + 999999999) / 1000000000 ))
# Peak is one sarek run (or one nomics run) next to the FASTQs: about 5x the FASTQ size, plus margin.
need_gb=$(( 5 * fq_gb + 25 + missing_bytes / 1000000000 ))
room=$REPO/scripts/goldens/make_room.sh   # site helper: moves passed results to other storage
[[ -x $room ]] && "$room" "$need_gb"
have_gb=$(free_gb "$OMICS")
log "FASTQs ${fq_gb} GB; need about ${need_gb} GB free, have ${have_gb} GB"
(( have_gb >= need_gb )) || { echo "::error::only ${have_gb} GB free on $OMICS; $id needs about ${need_gb} GB"; exit 2; }
fetch "$r1" "$FQ1" "$md5_1"; fetch "$r2" "$FQ2" "$md5_2"

target=$CI/target-cuda
cargo build --release --locked -p nomics --features cuda --target-dir "$target" 2>&1 | tail -n 5
NOMICS=$target/release/nomics

sheet=$DATA/samplesheet.csv params=$REPO/reference/params.matrix-wgs.yaml
printf 'patient,sex,status,sample,lane,fastq_1,fastq_2\n%s,%s,0,%s,L1,%s,%s\n' "$patient" "$sex" "$P" "$FQ1" "$FQ2" > "$sheet"
export NOMICS_SAREK_CONFIG=$REPO/reference/spark.config

drop_bulk() {  # drop_bulk RUN_DIR: keep reports, VCFs and tables; drop work/ and every BAM/CRAM
  rm -rf "$1/work"
  find "$1/results" \( -name '*.bam' -o -name '*.bai' -o -name '*.cram' -o -name '*.crai' -o -name '*.fastq.gz' \) -delete 2>/dev/null || true
}

# ---- one sarek run, reduced to its digest manifest and scatter order ----
sarek_golden() {  # sarek_golden r1|r2
  local r=$1 run=$RUNS/$id-$1 man=$RUNS/$id-$1.digest.tsv so=$RUNS/$id-$1.scatter-order.txt res chunks
  if [[ -s $man && -s $so ]]; then log "sarek $r digested earlier"; cp "$man" "$so" "$art/"; return; fi
  if grep -qx 'EXIT=0' "$run/exit_status" 2>/dev/null && [[ -d $run/work ]]; then
    log "sarek $r already done"
  else
    rm -rf "$run"
    log "sarek $r start"
    "$REPO/reference/run_sarek.sh" "$run" "$sheet" "$params" > "$RUNS/$id-$r.log" 2>&1 || true
    grep -qx 'EXIT=0' "$run/exit_status" 2>/dev/null || { echo "::error::sarek $r failed"; tail -n 60 "$RUNS/$id-$r.log"; tail -n 40 "$run/.nextflow.log" 2>/dev/null; exit 1; }
    log "sarek $r done: $(cat "$run/started_at") to $(cat "$run/finished_at")"
  fi
  log "disk: $r uses $(du -sh "$run" | cut -f1), $(free_gb "$OMICS") GB free"
  res=$run/results
  chunks=$(fastp_dir "$run") || { echo "::error::no fastp chunk directory in $run/work"; exit 1; }
  log "digest $r"
  "$NOMICS" dev digest --no-cache --prefix "$P" --output "$man.part" \
    --chunks "$chunks" --mapped "$res/preprocessing/mapped/$P" \
    --md-bam "$res/preprocessing/markduplicates/$P/$P.md.bam" \
    --md-metrics "$(find "$res/reports" -name '*.md.bam.metrics' | head -n 1)" \
    --recal-table "$res/preprocessing/recal_table/$P/$P.recal.table" \
    --recal-bam "$res/preprocessing/recalibrated/$P/$P.recal.bam" \
    --vcf "$res/variant_calling/haplotypecaller/$P/$P.haplotypecaller.vcf.gz"
  python3 "$REPO/scripts/scatter_order.py" "$run" "$P" > "$so.part"
  mv "$so.part" "$so"; mv "$man.part" "$man"
  log "sarek $r: $(($(wc -l < "$man") - 1)) artifacts, $(grep -c '^gather' "$so") scatter groups"
  cp "$man" "$so" "$art/"
  drop_bulk "$run"
  log "disk: dropped $r's BAMs and work tree, $(free_gb "$OMICS") GB free"
}

# ---- nomics, gated against one sarek run with that run's scatter order ----
# Always streamed (--spill-prep): a 10x genome is under nomics' 16 GiB in-memory limit, but holding
# it in memory needs more than the Spark's 121 GB, and earlyoom stopped G1 at markdup on 2026-11-07.
nomics_gate() {  # nomics_gate r1|r2
  local r=$1 out=$RUNS/$id-nomics-$1 saved=$RUNS/$id-nomics-$1.gates.txt
  if grep -qx 'verdict: PASS' "$saved" 2>/dev/null; then log "nomics against $r passed earlier"; cp "$saved" "$art/gates_$r.txt"; return; fi
  rm -rf "$out"
  log "nomics run, gated against $r"
  MATRIX_FQ1=$FQ1 MATRIX_FQ2=$FQ2 MATRIX_PREFIX=$P MATRIX_SAMPLE=$SM \
  MATRIX_SCATTER_ORDER=$RUNS/$id-$r.scatter-order.txt MATRIX_DIGEST=$RUNS/$id-$r.digest.tsv \
  NOMICS_BIN_DIR=$target/release NOMICS_RUN_EXTRA="${NOMICS_RUN_EXTRA:-} --spill-prep" \
    "$REPO/scripts/full_exome.sh" matrixwgs "$out" > "$art/nomics_$r.log" 2>&1 \
    || { rc=$?
         # 143 = SIGTERM and 137 = SIGKILL from outside the run (an OOM daemon, a timeout, a kill).
         echo "::error::nomics run failed with exit $rc at $(date -u +%H:%M:%SZ) (kept $out)"
         free -g | sed 's/^/  /'
         tail -n 60 "$art/nomics_$r.log"; tail -n 60 "$out/run.log" 2>/dev/null; exit 1; }
  cp "$out/gates.txt" "$saved"; cp "$out/gates.txt" "$art/gates_$r.txt"
  cp "$out/stages.tsv" "$art/stages_$r.tsv" 2>/dev/null || true
  grep -v '^PASS' "$out/gates.txt" | sed 's/^/  /' || true
  log "disk: nomics output $(du -sh "$out" | cut -f1), $(free_gb "$OMICS") GB free"
  if grep -qx 'verdict: PASS' "$saved"; then
    rm -rf "$out"
  else
    # Keep our side of a failed gate for debugging, without its chunk FASTQs; stop here.
    rm -rf "$out/prep"
    echo "::error::nomics differs from sarek $r (see gates_$r.txt; kept $out without its chunk FASTQs)"
    exit 1
  fi
}

sarek_golden r1
nomics_gate r1
sarek_golden r2
nomics_gate r2

# ---- r1 against r2 ----
"$NOMICS" dev digest-compare --ours "$RUNS/$id-r2.digest.tsv" --golden "$RUNS/$id-r1.digest.tsv" > "$art/compare_r1_r2.txt" || true
cp "$art/compare_r1_r2.txt" "$RUNS/$id-compare_r1_r2.txt"
differs=$(awk '/^FAIL/ {printf "%s ", $2}' "$art/compare_r1_r2.txt")
differs=${differs% }
scatter_same=yes; diff -q <(grep -v '^#' "$RUNS/$id-r1.scatter-order.txt") <(grep -v '^#' "$RUNS/$id-r2.scatter-order.txt") > /dev/null || scatter_same=no
# Up to MarkDuplicates nothing depends on the scatter: those artifacts must be identical. BQSR
# and later may differ only when the two runs' scatter records differ.
repro=PASS
for a in $differs; do
  case $a in
    recal_table|recal_bam|vcf) [[ $scatter_same == no ]] || repro=FAIL ;;
    *) repro=FAIL ;;
  esac
done
log "r1 against r2: ${differs:-nothing} differs; scatter records the same: $scatter_same"

n1=$(grep -c '^PASS' "$art/gates_r1.txt" || true); t1=$(grep -cE '^(PASS|FAIL)' "$art/gates_r1.txt" || true)
n2=$(grep -c '^PASS' "$art/gates_r2.txt" || true); t2=$(grep -cE '^(PASS|FAIL)' "$art/gates_r2.txt" || true)
records=$(awk -F'\t' '$1 == "vcf" {print $3}' "$RUNS/$id-r1.digest.tsv")
printf 'id\tpatient\tsarek_r1_r2\tr1_r2_differ\tnomics_vs_r1\tnomics_vs_r2\tvcf_records\n%s\t%s\t%s\t%s\t%s/%s\t%s/%s\t%s\n' \
  "$id" "$patient" "$repro" "${differs:-none}" "$n1" "$t1" "$n2" "$t2" "$records" | tee "$art/outcome.tsv"
ok=no; [[ $repro == PASS && $n1 == "$t1" && $n2 == "$t2" && $t1 -gt 0 ]] && ok=yes
# The FASTQs are Google's and can be fetched again; after a full pass they make room for the next genome.
[[ $ok == yes ]] && { rm -f "$FQ1" "$FQ2"; log "dropped the FASTQs; $(free_gb "$OMICS") GB free"; }
[[ $ok == yes ]]
