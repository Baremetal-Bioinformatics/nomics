#!/bin/bash
set -euo pipefail
golden="$HOME/omics/runs/somatic/S1-r1/golden"
out=$(mktemp -d /tmp/nomics-contamination-gate.XXXXXX)
printf '%s\n' "$out" > /tmp/nomics-contamination-gate-path
for role in TUMOR NORMAL; do
 if [ "$role" = TUMOR ]; then sample=HCC1395_WES; else sample=HCC1395BL_WES; fi
 target/release/nomics somatic gather-pileup-summaries \
  --I "$golden/GETPILEUPSUMMARIES_${role}-${sample}/${sample}.mutect2.chr17_38877931-38878296.pileups.table" \
  --I "$golden/GETPILEUPSUMMARIES_${role}-${sample}.2/${sample}.mutect2.chr1_11981-12351.pileups.table" \
  --sequence-dictionary "$HOME/omics/ref/GATK.GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.dict" \
  --output "$out/${sample}.mutect2.pileups.table"
 cmp "$out/${sample}.mutect2.pileups.table" "$golden/GATHERPILEUPSUMMARIES_${role}-${sample}/${sample}.mutect2.pileups.table"
 echo "cmp gathered $role: PASS"
done
target/release/nomics somatic calculate-contamination \
 --input "$out/HCC1395_WES.mutect2.pileups.table" \
 --matched-normal "$out/HCC1395BL_WES.mutect2.pileups.table" \
 --output "$out/contamination.table" --tumor-segmentation "$out/segmentation.table"
for table in contamination segmentation; do
 cmp "$out/${table}.table" "$golden/CALCULATECONTAMINATION-HCC1395_WES_vs_HCC1395BL_WES/HCC1395_WES_vs_HCC1395BL_WES.mutect2.${table}.table"
 echo "cmp $table: PASS"
done
