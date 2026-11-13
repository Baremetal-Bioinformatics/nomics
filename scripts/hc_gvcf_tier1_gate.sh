#!/usr/bin/env bash
# HC-GVCF harness: nomics hc --emit-gvcf on the tier1joint recal.bam, then joint-genotype; gate both.
# Usage (on the Spark, under the resource lock): scripts/hc_gvcf_tier1_gate.sh OUT [perf]
# SRC: a tier1joint run dir with HG002_WES.recal.bam; BIN: the nomics binary. NOMICS_TIMING=1 is set in
# timed mode and costs ~25% of HC wall (board S-11): compare runs of this harness with each other only.
set -euo pipefail
SRC=${SRC:-$HOME/omics/runs/nomics-full-tier1joint-20261016-1527}
OUT=$1; MODE=${2:-}
REPO=$(cd "$(dirname "$0")/.." && pwd)
BIN=${BIN:-$REPO/target/release/nomics}
REF=$HOME/omics/ref/GATK.GRCh38
DB=$REF/Annotation/GATKBundle/dbsnp_146.hg38.vcf.gz
GJ=$HOME/omics/runs/tier1-joint-x86-20261016/r1
mkdir -p "$OUT"
CMD=("$BIN" hc --input "$SRC/HG002_WES.recal.bam" --intervals "$HOME/omics/data/Agilent_v7.bed" \
  --reference "$REF/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta" --dbsnp "$DB" \
  --output "$OUT/HG002_WES.haplotypecaller.g.vcf.gz" --threads 20 --emit-gvcf)
if [ "$MODE" = perf ]; then perf record -F 499 -g -o "$OUT/perf.data" -- "${CMD[@]}" > "$OUT/run.log" 2>&1
else NOMICS_TIMING=1 /usr/bin/time -f "%e s wall %M KB" "${CMD[@]}" > "$OUT/run.log" 2>&1; fi
"$BIN" dev gate-gvcf --ours "$OUT/HG002_WES.haplotypecaller.g.vcf.gz" --golden "$GJ/HG002_WES/HG002_WES.haplotypecaller.g.vcf.gz" > "$OUT/gate_gvcf.log" 2>&1 || true
"$BIN" joint-genotype --gvcf "$OUT/HG002_WES.haplotypecaller.g.vcf.gz" --dbsnp "$DB" -o "$OUT/joint_germline.vcf.gz" --threads 20 > "$OUT/joint.log" 2>&1 || true
"$BIN" dev gate-vcf --ours "$OUT/joint_germline.vcf.gz" --golden "$GJ/joint_variant_calling/joint_germline.vcf.gz" > "$OUT/gate_joint.log" 2>&1 || true
echo "gvcf: $(tail -1 "$OUT/gate_gvcf.log")  joint: $(tail -1 "$OUT/gate_joint.log")"; grep -v "^$" "$OUT/run.log" | grep -iv "kallsyms\|relocation\|vmlinux\|kernel" | grep -v "^t=" | tail -6
