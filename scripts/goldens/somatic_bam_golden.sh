#!/usr/bin/env bash
# BAM-started sarek Mutect2 golden. Usage: somatic_bam_golden.sh ARTIFACT_DIR
# Environment: PAIR_ID, PATIENT, TUMOUR_SAMPLE, NORMAL_SAMPLE (samplesheet names),
# TUMOUR_BAM, NORMAL_BAM, TUMOUR_BAI, NORMAL_BAI, SEX, INTERVALS (.interval_list or BED),
# REF, DICT, BUNDLE, DATA, RUN_ROOT, RUNS (r1 or "r1 r2"), MIN_FREE_GB (150).
# No input downloads or automatic deletion of failed runs/work trees.
set -euo pipefail
REPO=$(cd "$(dirname "$0")/../.." && pwd)
ART=${1:?usage: somatic_bam_golden.sh ARTIFACT_DIR}
OMICS=$HOME/omics
PAIR_ID=${PAIR_ID:-S2}; PATIENT=${PATIENT:-HCC1143}; SEX=${SEX:-XX}
TUMOUR_SAMPLE=${TUMOUR_SAMPLE:-tumor}; NORMAL_SAMPLE=${NORMAL_SAMPLE:-normal}
TUMOUR_BAM=${TUMOUR_BAM:-$OMICS/data/somatic2/hcc1143/hcc1143_T_clean.bam}
NORMAL_BAM=${NORMAL_BAM:-$OMICS/data/somatic2/hcc1143/hcc1143_N_clean.bam}
TUMOUR_BAI=${TUMOUR_BAI:-${TUMOUR_BAM%.bam}.bai}; NORMAL_BAI=${NORMAL_BAI:-${NORMAL_BAM%.bam}.bai}
INTERVALS=${INTERVALS:-$OMICS/data/somatic2/hcc1143/exome_calling_regions.v1.interval_list}
REF=${REF:-$OMICS/ref/GATK.GRCh38}
DICT=${DICT:-$REF/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.dict}
BUNDLE=${BUNDLE:-$REF/Annotation/GATKBundle}
DATA=${DATA:-$OMICS/data/somatic2/$PAIR_ID}; RUN_ROOT=${RUN_ROOT:-$OMICS/runs/somatic}
RUNS=${RUNS:-r1}
[[ $PAIR_ID =~ ^[A-Za-z0-9_-]+$ ]] || { echo '::error::unsafe PAIR_ID'; exit 2; }
[[ $RUNS == r1 || $RUNS == 'r1 r2' ]] || { echo '::error::RUNS must be r1 or "r1 r2"'; exit 2; }
[[ $SEX == XX || $SEX == XY || $SEX == NA ]] || { echo '::error::invalid SEX'; exit 2; }
for value in "$PATIENT" "$TUMOUR_SAMPLE" "$NORMAL_SAMPLE"; do
  [[ $value =~ ^[A-Za-z0-9_-]+$ ]] || { echo '::error::unsafe patient/sample name'; exit 2; }
done
[[ $TUMOUR_SAMPLE != "$NORMAL_SAMPLE" ]] || { echo '::error::samples must differ'; exit 2; }
mkdir -p "$DATA" "$RUN_ROOT" "$ART"
DATA=$(cd "$DATA" && pwd); ART=$(cd "$ART" && pwd); RUN_ROOT=$(cd "$RUN_ROOT" && pwd)
. "$REPO/scripts/goldens/lib.sh"
have=$(free_gb "$RUN_ROOT")
log "$PAIR_ID: BAM-started Mutect2; $have GB free"
(( have >= ${MIN_FREE_GB:-150} )) || { echo "::error::$have GB free, need ${MIN_FREE_GB:-150}"; exit 1; }
export REPO ART DATA REF PATIENT SEX TUMOUR_SAMPLE NORMAL_SAMPLE TUMOUR_BAM NORMAL_BAM TUMOUR_BAI NORMAL_BAI INTERVALS DICT BUNDLE
python3 -I "$REPO/scripts/goldens/somatic_bam_prepare.py"
export NOMICS_SAREK_CONFIG=${NOMICS_SAREK_CONFIG:-$REPO/reference/spark.config}
for r in $RUNS; do
  run=$RUN_ROOT/$PAIR_ID-$r
  if grep -qx 'EXIT=0' "$run/exit_status" 2>/dev/null; then log "sarek $r already done"; continue; fi
  [[ ! -e $run ]] || { echo "::error::$run exists without EXIT=0; inspect and move it before retrying"; exit 2; }
  log "sarek $r start"
  "$REPO/reference/run_sarek.sh" "$run" "$DATA/samplesheet.csv" "$DATA/params.yaml" > "$RUN_ROOT/$PAIR_ID-$r.log" 2>&1
  grep -qx 'EXIT=0' "$run/exit_status" || { echo "::error::sarek $r failed; see $RUN_ROOT/$PAIR_ID-$r.log"; exit 1; }
  cp "$run"/results/pipeline_info/execution_trace*.txt "$ART/trace-$r.txt"
done
if [[ $RUNS == 'r1 r2' ]]; then
  "$REPO/reference/compare_runs.sh" "$RUN_ROOT/$PAIR_ID-r1" "$RUN_ROOT/$PAIR_ID-r2" > "$ART/compare_r1_r2.tsv"
  python3 -I "$REPO/scripts/goldens/somatic_work_compare.py" "$RUN_ROOT/$PAIR_ID-r1" "$RUN_ROOT/$PAIR_ID-r2" > "$ART/compare_work_r1_r2.tsv"
fi
if [[ ! -d $RUN_ROOT/$PAIR_ID-r1/golden ]]; then
  "$REPO/scripts/goldens/somatic_harvest.sh" "$RUN_ROOT/$PAIR_ID-r1"
else
  log 'r1 golden already exists; leaving it untouched (somatic_harvest.sh --force to replace)'
fi
echo 'SOMATIC BAM GOLDEN DONE'
