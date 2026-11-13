#!/usr/bin/env bash
# One sample-matrix exome on the Spark: fetch, two sarek 3.10.0 reference runs (r1, r2),
# compare r1 with r2, then a full nomics run gated at every stage boundary against r1.
#
#   scripts/goldens/exome_golden.sh ID ARTIFACT_DIR
#
# ID is a row of reference/matrix/exomes.tsv (E*, real exomes) or reference/matrix/fuzz.tsv (F*,
# synthetic reads from scripts/fuzz/make_fuzz.py, run through the same exome flow over the
# windows the reads came from).
#
# Resumable: verified inputs and finished sarek runs (exit_status EXIT=0) are kept and reused.
# Exit 0 only when r1 == r2 in content and every nomics gate passes.
set -euo pipefail

id=$1 art=$2
REPO=$(cd "$(dirname "$0")/../.." && pwd)
export PATH="$HOME/.cargo/bin:/usr/local/cuda/bin:$PATH"
OMICS=$HOME/omics
REF=$OMICS/ref/GATK.GRCh38
DATA=$OMICS/data/matrix/$id
RUNS=$OMICS/runs/matrix
CI=${NOMICS_CI_HOME:-$OMICS/ci}
mkdir -p "$DATA" "$RUNS" "$art"
. "$REPO/scripts/goldens/lib.sh"

fuzz=false
if [[ $id == F* ]]; then
  fuzz=true
  row=$(awk -F'\t' -v id="$id" '!/^#/ && $1 == id' "$REPO/reference/matrix/fuzz.tsv")
  [[ -n $row ]] || { echo "unknown sample: $id" >&2; exit 2; }
  IFS=$'\t' read -r _ scenario seed pairs <<<"$row"
  patient=FUZZ sample=FUZZ_$scenario sex=XX
  MIN_FREE_GB=${MIN_FREE_GB:-40}
  log "$id: fuzz scenario $scenario, seed $seed, $pairs background pairs"
else
  row=$(awk -F'\t' -v id="$id" '!/^#/ && $1 == id' "$REPO/reference/matrix/exomes.tsv")
  [[ -n $row ]] || { echo "unknown sample: $id" >&2; exit 2; }
  IFS=$'\t' read -r _ patient sample sex source r1 r2 md5_1 md5_2 bedname <<<"$row"
  IFS=$'\t' read -r _ bedfile bedurl bedmd5 < <(awk -F'\t' -v n="$bedname" '!/^#/ && $1 == n' "$REPO/reference/matrix/beds.tsv")
  log "$id: $patient $sample, $source $r1, BED $bedfile"
fi

# A sarek run of this sample that didn't finish runs again from scratch below. Drop it now, so a
# partial work tree (up to a whole run's size) counts as free space instead of failing the check.
for r in r1 r2; do
  run=$RUNS/$id-$r
  if [[ -d $run ]] && ! grep -qx 'EXIT=0' "$run/exit_status" 2>/dev/null; then
    log "dropping unfinished sarek $r ($(du -sh "$run" | cut -f1)); it runs again"
    rm -rf "$run"
  fi
done
# Disk: a whole exome needs about 220 GB (fuzz rows set 40 above). Once
# r1 is done, r2 needs about r1's work tree (r2 publishes by hard link), and nomics runs after r2's
# work tree is dropped.
if [[ -n ${MIN_FREE_GB:-} ]]; then
  need=$MIN_FREE_GB
elif grep -qx 'EXIT=0' "$RUNS/$id-r1/exit_status" 2>/dev/null; then
  need=$(( $(du -s -BG "$RUNS/$id-r1/work" 2>/dev/null | cut -f1 | tr -dc 0-9) + 20 ))
else
  need=220
fi
room=$REPO/scripts/goldens/make_room.sh   # site helper: moves passed results to other storage
[[ -x $room ]] && "$room" "$need"
free_gb=$(df -BG --output=avail "$OMICS" | tail -n 1 | tr -dc 0-9)
log "free on $OMICS: ${free_gb} GB, need ${need}"
(( free_gb >= need )) || { echo "::error::only ${free_gb} GB free on $OMICS, need ${need}"; exit 2; }

# ---- inputs: generated (fuzz) or fetched and verified by md5 ----
FQ1=$DATA/${sample}_1.fastq.gz FQ2=$DATA/${sample}_2.fastq.gz
if $fuzz; then
  # Regenerated every time: the generator is deterministic, so a rerun gives the same bytes, and
  # the md5s in the artifact show which reads this golden was made from.
  gen=$DATA/gen
  rm -rf "$gen"
  python3 -I "$REPO/scripts/fuzz/make_fuzz.py" "$scenario" "$seed" "$pairs" \
    "$REF/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta" "$OMICS/data/Agilent_v7.bed" "$gen" | tee "$art/fuzz_summary.txt"
  mv "$gen/R1.fastq.gz" "$FQ1"; mv "$gen/R2.fastq.gz" "$FQ2"
  BED=$DATA/targets.bed; mv "$gen/targets.bed" "$BED"; rmdir "$gen"
  md5sum "$FQ1" "$FQ2" "$BED" | sed "s|$DATA/||" | tee "$art/fuzz_md5.txt"
  if [[ -f $RUNS/$id-inputs.md5 ]] && ! cmp -s "$RUNS/$id-inputs.md5" "$art/fuzz_md5.txt"; then
    # Finished sarek runs were made from other reads (a changed generator): they can't be reused.
    log "generated inputs differ from the ones the saved sarek runs used; running sarek again"
    rm -rf "$RUNS/$id-r1" "$RUNS/$id-r2" "$RUNS/$id-compare_r1_r2.tsv" "$RUNS/$id-r1.scatter-order.txt"
  fi
  cp "$art/fuzz_md5.txt" "$RUNS/$id-inputs.md5"
elif [[ $source == ena ]]; then
  urls=$(curl -fsS --retry 5 "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$r1&result=read_run&fields=fastq_ftp,fastq_md5&format=tsv" | awk -F'\t' 'NR==2 {print $2"\t"$3}')
  IFS=$'\t' read -r ftps md5s <<<"$urls"
  IFS=';' read -ra U <<<"$ftps"; IFS=';' read -ra M <<<"$md5s"
  u1= u2=
  for i in "${!U[@]}"; do
    [[ ${M[$i]} == "$md5_1" ]] && u1=https://${U[$i]}
    [[ ${M[$i]} == "$md5_2" ]] && u2=https://${U[$i]}
  done
  [[ -n $u1 && -n $u2 ]] || { echo "::error::ENA lists no files with the pinned md5s for $r1: $ftps / $md5s"; exit 2; }
  fetch "$u1" "$FQ1" "$md5_1"; fetch "$u2" "$FQ2" "$md5_2"
else
  fetch "$r1" "$FQ1" "$md5_1"; fetch "$r2" "$FQ2" "$md5_2"
fi
if $fuzz; then
  :
elif [[ $bedurl == - ]]; then
  BED=$OMICS/data/$bedfile
  got=$(md5sum "$BED" | cut -c1-32)
  [[ $got == "$bedmd5" ]] || log "note: $BED md5 is $got (nf-core's copy is $bedmd5); using the Spark's file, as the HG002 goldens did"
else
  BED=$OMICS/data/matrix/beds/$bedfile; mkdir -p "$(dirname "$BED")"
  fetch "$bedurl" "$BED" "$bedmd5"
fi

# ---- sarek r1 and r2 ----
sheet=$DATA/samplesheet.csv params=$DATA/params.yaml
printf 'patient,sex,status,sample,lane,fastq_1,fastq_2\n%s,%s,0,%s,L1,%s,%s\n' "$patient" "$sex" "$sample" "$FQ1" "$FQ2" > "$sheet"
sed "s|\${HOME}/omics/data/Agilent_v7.bed|$BED|" "$REPO/reference/params.tier1.yaml" > "$params"
# Hard links instead of copies: the published files are the same bytes, at about half the disk.
echo 'publish_dir_mode: link' >> "$params"
grep -q "^intervals: $BED\$" "$params" || { echo "::error::params: intervals not set"; exit 2; }
export NOMICS_SAREK_CONFIG=$REPO/reference/spark.config
saved_compare=$RUNS/$id-compare_r1_r2.tsv compared=$RUNS/$id-r2/compared_with_r1
for r in r1 r2; do
  run=$RUNS/$id-$r
  if grep -qx 'EXIT=0' "$run/exit_status" 2>/dev/null; then
    # r2's BAMs are dropped after the comparison; a finished r2 without them and without a saved
    # comparison can't be compared again, so it runs again.
    if [[ $r == r2 && ! -f $compared ]] && [[ -z $(find "$run/results" -name '*.bam' -print -quit 2>/dev/null) ]]; then
      log "sarek r2 finished but its BAMs are gone and no comparison was saved; running r2 again"
    else
      log "sarek $r already done"; continue
    fi
  fi
  rm -rf "$run"
  log "sarek $r start"
  "$REPO/reference/run_sarek.sh" "$run" "$sheet" "$params" > "$RUNS/$id-$r.log" 2>&1 || true
  grep -qx 'EXIT=0' "$run/exit_status" 2>/dev/null || { echo "::error::sarek $r failed"; tail -n 60 "$RUNS/$id-$r.log"; tail -n 40 "$run/.nextflow.log" 2>/dev/null; exit 1; }
  log "sarek $r done: $(cat "$run/started_at") to $(cat "$run/finished_at")"
  # The comparison reads published results only, so r2's work tree goes as soon as r2 is done.
  [[ $r == r2 ]] && rm -rf "$run/work"
  log "disk: $r uses $(du -sh "$run" | cut -f1), $(df -BG --output=avail "$OMICS" | tail -n 1 | tr -d ' ') free"
done

# ---- r1 vs r2 ----
if [[ -f $compared && -s $saved_compare ]]; then
  log "r1 r2 compared earlier ($(cat "$compared")); reusing that result"
  cp "$saved_compare" "$art/compare_r1_r2.tsv"
else
  log "compare r1 r2"
  "$REPO/reference/compare_runs.sh" "$RUNS/$id-r1" "$RUNS/$id-r2" > "$art/compare_r1_r2.tsv"
  cp "$art/compare_r1_r2.tsv" "$saved_compare"
  date -u +%Y-%m-%dT%H:%MZ > "$compared"
fi
# ---- scatter order, when sarek scattered ----
scatter=
if grep -q ':GATK4_GATHERBQSRREPORTS (' "$RUNS/$id-r1"/results/pipeline_info/execution_trace*.txt; then
  scatter=$RUNS/$id-r1.scatter-order.txt
  python3 "$REPO/scripts/scatter_order.py" "$RUNS/$id-r1" "$sample" > "$scatter"
  cp "$scatter" "$art/"
  log "sarek scattered: $(grep -c '^gather' "$scatter") gather lines"
fi

repro=PASS
# A scattered run's recal.bam carries samtools merge's time-seeded @RG/@PG ID suffixes on every
# record (docs/REFERENCE-VARIABILITY.md), so r1 and r2 cannot match there; nomics is gated on it
# against r1 with r1's recovered seed instead.
python3 "$REPO/scripts/goldens/compare_summary.py" ${scatter:+--scattered} "$art/compare_r1_r2.tsv" > "$art/compare_summary.txt" || repro=FAIL
cat "$art/compare_summary.txt"
# r2 has served its purpose; keep its report and VCF, drop its work and BAMs.
rm -rf "$RUNS/$id-r2/work"; find "$RUNS/$id-r2/results" -name '*.bam' -delete

# ---- nomics, gated against r1 ----
target=$CI/target-cuda
cargo build --release --locked -p nomics --features cuda --target-dir "$target" 2>&1 | tail -n 5
out=$RUNS/$id-nomics
rm -rf "$out"
log "nomics run + gates"
# nomics holds every chunk in memory below 16 GiB of gzipped FASTQ and streams above it. In memory,
# markdup's merge needs about 7x the gzipped FASTQ plus the 22 GB index: an 88M-pair exome needs
# more than the Spark's 121 GB, and earlyoom stops the run. Stream anything above 8 GB.
fq_gz=$(( $(stat -c %s "$FQ1") + $(stat -c %s "$FQ2") ))
extra=${NOMICS_RUN_EXTRA:-}
if (( fq_gz > 8000000000 )); then extra="$extra --spill-prep"; log "FASTQs $((fq_gz / 1000000000)) GB gzipped: streaming (--spill-prep)"; fi
NOMICS_RUN_EXTRA=$extra MATRIX_FQ1=$FQ1 MATRIX_FQ2=$FQ2 MATRIX_GOLD=$RUNS/$id-r1 MATRIX_BED=$BED MATRIX_PREFIX=$sample \
MATRIX_SAMPLE=${patient}_${sample} MATRIX_SCATTER_ORDER=$scatter NOMICS_BIN_DIR=$target/release \
  "$REPO/scripts/full_exome.sh" matrix "$out" > "$art/nomics.log" 2>&1 || { echo "::error::nomics run failed"; tail -n 60 "$art/nomics.log"; tail -n 60 "$out/run.log" 2>/dev/null; exit 1; }
cp "$out/gates.txt" "$out/stages.tsv" "$art/" 2>/dev/null || true
log "disk: nomics output $(du -sh "$out" | cut -f1), $(df -BG --output=avail "$OMICS" | tail -n 1 | tr -d ' ') free"
grep -q '^FAIL' "$out/gates.txt" 2>/dev/null && cp "$out"/gate_*.log "$art/" 2>/dev/null || true
# The artifact may not be readable where the result is read, so the job log carries each failed gate's report.
for g in $(awk '/^FAIL/ {print $2}' "$out/gates.txt" 2>/dev/null); do
  echo "---- gate_$g.log (first 40 lines)"; head -n 40 "$out/gate_$g.log" 2>/dev/null || true
done
n_pass=$(grep -c '^PASS' "$out/gates.txt" || true); n_all=$(grep -c . "$out/gates.txt" || true)
records=$(zcat "$out/$sample.haplotypecaller.vcf.gz" 2>/dev/null | grep -vc '^#' || true)
gates=FAIL; [[ $n_all -gt 0 && $n_pass == "$n_all" ]] && gates=PASS
# A VCF gate over a nearly empty VCF proves nothing about HaplotypeCaller: nine of the first
# fuzz goldens passed 18/18 with 0 records (2026-11-10). Fewer than the minimum fails the golden.
min_records=${NOMICS_MIN_VCF_RECORDS:-100}
if (( ${records:-0} < min_records )); then
  echo "::error::$id: the VCF has ${records:-0} records, fewer than $min_records; its gate compared too little"
  gates=FAIL
fi
# Disk: after a full PASS keep r1's published results and the reports, drop the work trees
# (the gates already read r1's fastp chunks from work/). After a FAIL keep everything.
if [[ $gates == PASS && $repro == PASS ]]; then rm -rf "$out/prep" "$out/mapped" "$RUNS/$id-r1/work"; fi
# A passed real exome's sarek results then move to other storage, where a site helper exists, so
# the genomes have room. A failed move leaves them in place and does not fail the golden.
mover=$REPO/scripts/goldens/nas_move.sh
if [[ $gates == PASS && $repro == PASS && $fuzz == false && -x $mover ]]; then
  "$mover" "$id" || log "the move did not finish; $id stays on this machine"
fi

printf 'id\tpatient\tsarek_r1_r2\tnomics_gates\tpassed\ttotal\tvcf_records\tscattered\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
  "$id" "$patient" "$repro" "$gates" "$n_pass" "$n_all" "$records" "$([[ -n $scatter ]] && echo yes || echo no)" | tee "$art/outcome.tsv"
cat "$out/gates.txt"
[[ $repro == PASS && $gates == PASS ]]
