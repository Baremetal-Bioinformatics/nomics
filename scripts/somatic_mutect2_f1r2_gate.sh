#!/usr/bin/env bash
# GATK 4.6.2.0 collector payload gate; no Java. Default all seven slice goldens.
# Compared under the f1r2 rule in docs/GATES.md (scripts/somatic_f1r2_cmp.py): member
# names as a set, alt_histogram by column label, every other byte as it is.
# NOMICS, M2_GOLDENS and REF have the same meanings as the regions gate.
# M2_F1R2_OUT optionally keeps generated archives for downstream checks.
set -uo pipefail
REPO=$(cd "$(dirname "$0")/.." && pwd)
bin=${NOMICS:-$REPO/target/release/nomics}
root=${M2_GOLDENS:-$HOME/omics/runs/somatic}
ref=${REF:-$HOME/omics/ref/GATK.GRCh38}/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta
contigs=("$@"); (( ${#contigs[@]} )) || contigs=(chr21 chr22 chr19 chr20 chr13 chr16 chr12)
T=$(mktemp -d "${TMPDIR:-/tmp}/m2-f1r2-gate.XXXXXX"); trap 'rm -rf "$T"' EXIT
out=${M2_F1R2_OUT:-$T}; mkdir -p "$out"
pass=0 fail=0 skip=0
for c in "${contigs[@]}"; do
  D=$root/m2-$c
  golden=$D/oracle/HCC1395_WES_vs_HCC1395BL_WES.mutect2.$c.f1r2.tar.gz
  [[ -f $golden ]] || { skip=$((skip+1)); echo "SKIP $c (no $golden)"; continue; }
  s=$(date +%s)
  if ! "$bin" somatic mutect2-f1r2 --input "$D/HCC1395BL_WES.$c.bam" --input "$D/HCC1395_WES.$c.bam" \
      --normal-sample HCC1395_HCC1395BL_WES --reference "$ref" \
      --panel-of-normals "$D/1000g_pon.$c.vcf.gz" --germline-resource "$D/af-only-gnomad.$c.vcf.gz" \
      --intervals "$D/$c.bed" --f1r2-tar-gz "$out/$c.f1r2.tar.gz" 2> "$T/filters"; then
    fail=$((fail+1)); echo "FAIL $c ($(tail -n 1 "$T/filters"))"; continue
  fi
  if python3 -I "$REPO/scripts/somatic_f1r2_cmp.py" "$golden" "$out/$c.f1r2.tar.gz" "$c"
  then pass=$((pass+1)); else fail=$((fail+1)); echo "FAIL $c"; fi
  echo "$c wall $(( $(date +%s)-s )) s"
done
echo "mutect2 f1r2 gate: $pass pass, $fail fail, $skip skipped"
(( fail==0 && pass>0 ))
