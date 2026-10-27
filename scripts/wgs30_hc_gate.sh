#!/usr/bin/env bash
# Whole-genome HaplotypeCaller gate (tier "wgs30hc"): nomics run --hc-only on the
# 30x recal.bam (identical to sarek's, digest-gated) and the merged VCF compared
# record by record with the 30x golden. Catches caller bugs in regions the exome
# and chr20 tiers never reach, in ~25 min instead of a 3 h full run.
#
#   scripts/wgs30_hc_gate.sh [output-dir]   (on the Spark, under the resource lock)
set -euo pipefail
REPO=$(cd "$(dirname "$0")/.." && pwd)
BIN=$REPO/target/release/nomics
REF=$HOME/omics/ref/GATK.GRCh38
P=HG002_WGS30
SRC=${WGS30_RUN:-$HOME/omics/runs/nomics-wgs30-r1}
GOLD_VCF=$HOME/omics/runs/wgs30-x86-hcjava/$P.haplotypecaller.vcf.gz
OUT=${1:-$HOME/omics/runs/wgs30hc-$(date +%Y%m%d-%H%M)}
mkdir -p "$OUT"
ln -sf "$SRC/$P.recal.bam" "$OUT/$P.recal.bam"
"$BIN" run --resume-from hc --fastq1 /dev/null --fastq2 /dev/null --reference "$REF" \
  --targets "$REF/Annotation/intervals/wgs_calling_regions_noseconds.hg38.bed" --outdir "$OUT" \
  --dbsnp "$REF/Annotation/GATKBundle/dbsnp_146.hg38.vcf.gz" \
  --dense-sa-cache "$HOME/omics/cache/sa4.GRCh38.nsa" \
  --sample HG002_HG002_WGS30 --prefix $P --threads 20 > "$OUT/run.log" 2>&1
if "$BIN" dev gate-vcf --ours "$OUT/$P.haplotypecaller.vcf.gz" --golden "$GOLD_VCF" --examples 10 > "$OUT/gate_vcf.log" 2>&1; then
  echo "PASS vcf" > "$OUT/gates.txt"
else
  echo "FAIL vcf (see gate_vcf.log)" > "$OUT/gates.txt"
fi
echo "== $OUT"; column -t "$OUT/stages.tsv"; cat "$OUT/gates.txt"; head -3 "$OUT/gate_vcf.log"
