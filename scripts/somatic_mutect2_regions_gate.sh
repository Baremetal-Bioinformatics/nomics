#!/usr/bin/env bash
# Gate `nomics somatic mutect2-regions` (Mutect2's traversal, read filters, activity and
# assembly regions) against the slice goldens made by scripts/somatic_mutect2_slice.sh.
# Per contig (docs/GATES.md): the assembly-region file by bytes, the stats file by bytes, and
# the read-filter summary (every per-filter count, in order, and the totals) against GATK's log.
#
#   scripts/somatic_mutect2_regions_gate.sh [CONTIG ...]     default: chr21 chr22 chr19 chr20 chr13 chr16 chr12
#   NOMICS=path/to/nomics; goldens under ${M2_GOLDENS:-~/omics/runs/somatic}/m2-CONTIG
set -uo pipefail
REPO=$(cd "$(dirname "$0")/.." && pwd)
bin=${NOMICS:-$REPO/target/release/nomics}
root=${M2_GOLDENS:-$HOME/omics/runs/somatic}
ref=${REF:-$HOME/omics/ref/GATK.GRCh38}/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta
contigs=("$@"); (( ${#contigs[@]} )) || contigs=(chr21 chr22 chr19 chr20 chr13 chr16 chr12)
T=$(mktemp -d "${TMPDIR:-/tmp}/m2-regions-gate.XXXXXX"); trap 'rm -rf "$T"' EXIT
pass=0 fail=0 skip=0
for c in "${contigs[@]}"; do
  D=$root/m2-$c; X=HCC1395_WES_vs_HCC1395BL_WES.mutect2.$c
  [[ -f $D/oracle/regions.igv ]] || { skip=$((skip + 1)); echo "SKIP $c  (no $D/oracle/regions.igv)"; continue; }
  s=$(date +%s)
  if ! "$bin" somatic mutect2-regions --input "$D/HCC1395BL_WES.$c.bam" --input "$D/HCC1395_WES.$c.bam" \
       --normal-sample HCC1395_HCC1395BL_WES --reference "$ref" \
       --panel-of-normals "$D/1000g_pon.$c.vcf.gz" --germline-resource "$D/af-only-gnomad.$c.vcf.gz" --intervals "$D/$c.bed" \
       --assembly-region-out "$T/regions.igv" --stats "$T/stats" 2> "$T/filters"; then
    fail=$((fail + 1)); echo "FAIL $c  (the port failed: $(tail -n 1 "$T/filters"))"; continue
  fi
  problems=()
  cmp -s "$T/regions.igv" "$D/oracle/regions.igv" || problems+=("regions differ")
  cmp -s "$T/stats" "$D/oracle/$X.vcf.gz.stats" || problems+=("stats differ")
  diff <(grep -E "filtered by|filtered out of" "$D/oracle/gatk.log" | sed -E 's/^.*INFO  Mutect2 - //; s/ *$//') \
       <(grep -E "filtered by|filtered out of" "$T/filters" | sed 's/ *$//') > /dev/null || problems+=("read-filter summary differs")
  if (( ${#problems[@]} )); then fail=$((fail + 1)); echo "FAIL $c  (${problems[*]})"
  else
    pass=$((pass + 1))
    echo "PASS $c  ($(grep -c size= "$D/oracle/regions.igv") regions, $(awk -F'\t' '$4 ~ /size/ && $5 > 0' "$D/oracle/regions.igv" | wc -l | tr -d ' ') active; $(tail -n 1 "$D/oracle/$X.vcf.gz.stats" | cut -f2) callable; $(grep "filtered out of" "$T/filters" | sed 's/ *$//'); $(( $(date +%s) - s )) s)"
  fi
done
echo "mutect2 regions gate: $pass pass, $fail fail, $skip skipped"
(( fail == 0 && pass > 0 ))
