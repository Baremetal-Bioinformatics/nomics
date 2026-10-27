#!/usr/bin/env bash
# Timing experiment: BQSR alone at 30x, resumed from the 30x md.bam, at several
# thread counts (NOMICS_STOP_AFTER=bqsr). Checks the recal table byte for byte
# against the 30x run's; the full digest gate is scripts/wgs30_resume_gate.sh.
#
#   THREADS="10 20" scripts/wgs30_bqsr_sweep.sh   (Spark, under the lock)
set -euo pipefail
REPO=$(cd "$(dirname "$0")/.." && pwd)
BIN=$REPO/target/release/nomics
REF=$HOME/omics/ref/GATK.GRCh38
B=$REF/Annotation/GATKBundle
P=HG002_WGS30
SRC=${WGS30_RUN:-$HOME/omics/runs/nomics-wgs30-r1}
for t in ${THREADS:-10 20}; do
  OUT=$HOME/omics/runs/wgs30-bqsr-t$t-$(date +%Y%m%d-%H%M)
  mkdir -p "$OUT"
  ln -sf "$SRC/$P.md.bam" "$OUT/$P.md.bam"; cp "$SRC/$P.md.metrics" "$OUT/"
  NOMICS_STOP_AFTER=bqsr "$BIN" run --resume-from bqsr --fastq1 /dev/null --fastq2 /dev/null --reference "$REF" \
    --targets "$REF/Annotation/intervals/wgs_calling_regions_noseconds.hg38.bed" --outdir "$OUT" \
    --dbsnp "$B/dbsnp_146.hg38.vcf.gz" --known-sites "$B/dbsnp_146.hg38.vcf.gz" \
    --known-sites "$B/beta/Homo_sapiens_assembly38.known_indels.vcf.gz" \
    --known-sites "$B/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" \
    --known-sites-cache "$HOME/omics/cache/known_sites.GATK.GRCh38.nks" \
    --dense-sa-cache "$HOME/omics/cache/sa4.GRCh38.nsa" \
    --scatter-order "$REPO/reference/wgs30-x86-r1.scatter-order.txt" \
    --read-group-ds /data/omics/ref/GATK.GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta \
    --sample HG002_HG002_WGS30 --prefix $P --threads "$t" > "$OUT/run.log" 2>&1
  table=$(cmp -s "$OUT/$P.recal.table" "$SRC/$P.recal.table" && echo identical || echo DIFFERENT)
  echo "threads $t: recal table $table; $(grep -E 'bqsr stage (recalibrator|apply)' "$OUT/run.log" | tr '\n' ' ')"
  rm -f "$OUT/$P.recal.bam"
done
