#!/usr/bin/env bash
# VCF-mode HC harness on the wgs30s slice: resume from hc, gate the VCF.
# Usage (on the Spark, under the resource lock): scripts/hc_wgs30s_gate.sh OUT [perf]
# SRC: a wgs30s run dir with HG002_WGS30S.recal.bam; BIN: the nomics binary.
set -euo pipefail
SRC=${SRC:-$HOME/omics/runs/nomics-wgs30s-r1}
OUT=$1; MODE=${2:-}
REPO=$(cd "$(dirname "$0")/.." && pwd)
BIN=${BIN:-$REPO/target/release/nomics}
REF=$HOME/omics/ref/GATK.GRCh38
P=HG002_WGS30S
mkdir -p "$OUT"; ln -sf "$SRC/$P.recal.bam" "$OUT/$P.recal.bam"
CMD=("$BIN" run --resume-from hc --fastq1 /dev/null --fastq2 /dev/null --reference "$REF" \
  --targets "$REF/Annotation/intervals/wgs_calling_regions_noseconds.hg38.bed" --outdir "$OUT" \
  --dbsnp "$REF/Annotation/GATKBundle/dbsnp_146.hg38.vcf.gz" --dense-sa-cache "$HOME/omics/cache/sa4.GRCh38.nsa" \
  --scatter-order "$REPO/reference/wgs30s-x86-20261014-r1.scatter-order.txt" \
  --sample HG002_$P --prefix $P --threads 20)
if [ "$MODE" = perf ]; then perf record -F 299 -g -o "$OUT/perf.data" -- "${CMD[@]}" > "$OUT/run.log" 2>&1
else "${CMD[@]}" > "$OUT/run.log" 2>&1; fi
"$BIN" dev gate-vcf --ours "$OUT/$P.haplotypecaller.vcf.gz" --golden "$HOME/omics/runs/wgs30s-x86-20261014-r1/results/variant_calling/haplotypecaller/$P/$P.haplotypecaller.vcf.gz" > "$OUT/gate_vcf.log" 2>&1 || true
echo "vcf: $(tail -1 "$OUT/gate_vcf.log")"; cat "$OUT/stages.tsv"
