#!/usr/bin/env bash
# Decision run for the Mutect2 chain, before any port work: does sarek 3.10.0's
# tumour-normal Mutect2 output repeat from run to run?
#
#   scripts/goldens/somatic_golden.sh ARTIFACT_DIR        (on the Spark, under nohup)
#
# The pair is SEQC2's HCC1395 (tumour, SRR7890918) and HCC1395BL (normal, SRR7890919, which is
# sample-matrix exome E5), whole exome, Agilent V6 + UTR targets: sarek's own `test_full` pair,
# with a published truth set. It runs sarek twice with `tools: mutect2` and compares r1 with r2
# at every published boundary (reference/compare_runs.sh). There is no nomics gate yet: this
# decides whether an exact port has a stable target, and r1 becomes its golden.
#
# r1 keeps its work tree (the f1r2 tars and per-group VCFs are only there); r2's is dropped
# after the comparison. Needs about 500 GB free at the start.
set -euo pipefail
art=${1:?usage: somatic_golden.sh ARTIFACT_DIR}
REPO=$(cd "$(dirname "$0")/../.." && pwd)
OMICS=$HOME/omics
REF=$OMICS/ref/GATK.GRCh38
BUNDLE=$REF/Annotation/GATKBundle
DATA=$OMICS/data/somatic
RUNS=$OMICS/runs/somatic
id=S1
mkdir -p "$DATA" "$RUNS" "$art"
. "$REPO/scripts/goldens/lib.sh"

S3=https://ngi-igenomes.s3.amazonaws.com
N1=$OMICS/data/matrix/E5/HCC1395BL_WES_1.fastq.gz
N2=$OMICS/data/matrix/E5/HCC1395BL_WES_2.fastq.gz
T1=$DATA/HCC1395_WES_1.fastq.gz
T2=$DATA/HCC1395_WES_2.fastq.gz
BED=$OMICS/data/matrix/beds/S07604624_Padded_Agilent_SureSelectXT_allexons_V6_UTR.bed

have=$(free_gb "$OMICS")
log "$id: HCC1395 tumour vs HCC1395BL normal, WES; $have GB free"
(( have >= ${MIN_FREE_GB:-500} )) || { echo "::error::$have GB free, need ${MIN_FREE_GB:-500}"; exit 1; }

# ---- inputs: md5s are ENA's for SRR7890918 / SRR7890919 ----
fetch "$S3/test-data/sarek/SRR7890919_WES_HCC1395BL-EA_normal_1.fastq.gz" "$N1" 092cbd4c3a126509a3a6d06f14dca748
fetch "$S3/test-data/sarek/SRR7890919_WES_HCC1395BL-EA_normal_2.fastq.gz" "$N2" 9a5a26a33a62cbcfbaa170150c5417b5
fetch "$S3/test-data/sarek/SRR7890918_WES_HCC1395-EA_tumor_1.fastq.gz" "$T1" 8d8f82e4f35e1b4babc4bddb5150d50b
fetch "$S3/test-data/sarek/SRR7890918_WES_HCC1395-EA_tumor_2.fastq.gz" "$T2" 0bdf725a2e2926b480ddf8135413e795
[[ -f $BED ]] || { echo "::error::missing $BED (fetched by the E5 golden)"; exit 2; }
# Mutect2's germline resource and panel of normals: sarek's GATK.GRCh38 defaults. They carry no
# published md5, so the sizes are pinned (S3 listing, 2026-11-10) and the md5s recorded.
resource() {  # resource NAME BYTES
  local dest=$BUNDLE/$1
  if [[ ! -f $dest || $(stat -c%s "$dest") != "$2" ]]; then
    log "fetch $1"
    curl -fsSL --retry 10 -C - -o "$dest.part" "$S3/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/$1"
    [[ $(stat -c%s "$dest.part") == "$2" ]] || { echo "::error::$1: wrong size"; exit 2; }
    mv "$dest.part" "$dest"
  fi
}
resource af-only-gnomad.hg38.vcf.gz 3184275189
resource af-only-gnomad.hg38.vcf.gz.tbi 2443190
resource 1000g_pon.hg38.vcf.gz 17273497
resource 1000g_pon.hg38.vcf.gz.tbi 1534802
md5sum "$T1" "$T2" "$N1" "$N2" "$BED" "$BUNDLE"/af-only-gnomad.hg38.vcf.gz* "$BUNDLE"/1000g_pon.hg38.vcf.gz* \
  | sed "s|$OMICS/||" | tee "$art/inputs.md5"

# ---- sarek r1 and r2 ----
sheet=$DATA/samplesheet.csv params=$DATA/params.yaml
{
  echo 'patient,sex,status,sample,lane,fastq_1,fastq_2'
  echo "HCC1395,XX,0,HCC1395BL_WES,L1,$N1,$N2"
  echo "HCC1395,XX,1,HCC1395_WES,L1,$T1,$T2"
} > "$sheet"
sed -e "s|\${HOME}/omics/data/Agilent_v7.bed|$BED|" \
    -e 's|^tools: haplotypecaller$|tools: mutect2|' \
    -e '/^skip_tools: /d' "$REPO/reference/params.tier1.yaml" > "$params"
{
  echo "germline_resource: $BUNDLE/af-only-gnomad.hg38.vcf.gz"
  echo "germline_resource_tbi: $BUNDLE/af-only-gnomad.hg38.vcf.gz.tbi"
  echo "pon: $BUNDLE/1000g_pon.hg38.vcf.gz"
  echo "pon_tbi: $BUNDLE/1000g_pon.hg38.vcf.gz.tbi"
  # Hard links instead of copies: the same bytes at about half the disk.
  echo 'publish_dir_mode: link'
} >> "$params"
grep -q "^intervals: $BED\$" "$params" && grep -qx 'tools: mutect2' "$params" || { echo "::error::params not set"; exit 2; }
cp "$sheet" "$params" "$art/"
export NOMICS_SAREK_CONFIG=$REPO/reference/spark.config
for r in r1 r2; do
  run=$RUNS/$id-$r
  if grep -qx 'EXIT=0' "$run/exit_status" 2>/dev/null; then log "sarek $r already done"; continue; fi
  rm -rf "$run"
  log "sarek $r start"
  "$REPO/reference/run_sarek.sh" "$run" "$sheet" "$params" > "$RUNS/$id-$r.log" 2>&1 || true
  grep -qx 'EXIT=0' "$run/exit_status" 2>/dev/null || { echo "::error::sarek $r failed"; tail -n 60 "$RUNS/$id-$r.log"; tail -n 40 "$run/.nextflow.log" 2>/dev/null; exit 1; }
  log "sarek $r done: $(cat "$run/started_at") to $(cat "$run/finished_at"); $(du -sh "$run" | cut -f1), $(free_gb "$OMICS") GB free"
  cp "$run"/results/pipeline_info/execution_trace*.txt "$art/trace-$r.txt" 2>/dev/null || true
done

# ---- r1 vs r2 ----
log "compare r1 r2"
"$REPO/reference/compare_runs.sh" "$RUNS/$id-r1" "$RUNS/$id-r2" > "$art/compare_r1_r2.tsv"
# The unpublished per-group Mutect2 outputs, for the port's debugging and for where a difference starts.
{
  body() { zcat "$1" | { grep -v '^#' || true; } | md5sum | cut -c1-32; }
  for f in $(cd "$RUNS/$id-r1/work" && find . -type f \( -name '*f1r2*.tar.gz' -o -name '*mutect2*.vcf.gz' \) | sort); do
    name=$(basename "$f")
    b=$(find "$RUNS/$id-r2/work" -type f -name "$name" | head -n 1)
    [[ -n $b ]] || { printf 'work\t%s\texists\tFAIL\tmissing in r2\n' "$name"; continue; }
    case $name in
      *.tar.gz) a=$(tar -xzOf "$RUNS/$id-r1/work/$f" | md5sum | cut -c1-32); c=$(tar -xzOf "$b" | md5sum | cut -c1-32) ;;
      *) a=$(body "$RUNS/$id-r1/work/$f"); c=$(body "$b") ;;
    esac
    printf 'work\t%s\tcontents\t%s\t\n' "$name" "$([[ $a == "$c" ]] && echo PASS || echo DIFF)"
  done
} > "$art/compare_work_r1_r2.tsv" || true
rm -rf "$RUNS/$id-r2/work"
log "published files: $(grep -c . "$art/compare_r1_r2.tsv") checks, $(awk -F'\t' '$4 != "PASS"' "$art/compare_r1_r2.tsv" | wc -l) not PASS"
log "work files: $(grep -c . "$art/compare_work_r1_r2.tsv") checks, $(awk -F'\t' '$4 != "PASS"' "$art/compare_work_r1_r2.tsv" | wc -l) not PASS"
awk -F'\t' '$4 != "PASS"' "$art/compare_r1_r2.tsv" "$art/compare_work_r1_r2.tsv" | head -n 60
echo "SOMATIC DECISION RUN DONE"
