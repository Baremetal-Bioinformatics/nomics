#!/usr/bin/env bash
# Parabricks 4.7.1 on the x86 cloud box (8x RTX PRO 6000), raw FASTQs -> VCF.
#
#   scripts/bakeoff/parabricks_cloud.sh <outdir> [tier1|wgs20|wgs30 ...]
#
# Same reference bundle, known sites, BWA flags and read group as the sarek
# runs; HaplotypeCaller with -pcr-indel-model CONSERVATIVE. tier1 uses the
# Agilent v7 targets, the WGS inputs sarek's wgs_calling_regions. Times and
# GPU count per stage go to <outdir>/stages.tsv; VCFs are compared with the
# goldens on the Spark afterwards (scripts/bakeoff/concordance.py).
set -euo pipefail
OUT=${1:?outdir}; shift
SAMPLES=${*:-tier1 wgs20 wgs30}
O=/data/omics
IMG=nvcr.io/nvidia/clara/clara-parabricks:4.7.1-1
REF=$O/ref/GATK.GRCh38
FASTA=$REF/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta
BWA=$REF/Sequence/BWAIndex/Homo_sapiens_assembly38.fasta.64
B=$REF/Annotation/GATKBundle
NGPU=$(nvidia-smi -L | wc -l)
mkdir -p "$OUT/ref"
# Parabricks wants <fasta>.{amb,ann,bwt,pac,sa} beside the real FASTA.
for x in fasta fasta.fai dict; do ln -f "${FASTA%.fasta}.$x" "$OUT/ref/Homo_sapiens_assembly38.$x"; done
for x in amb ann bwt pac sa; do ln -f "$BWA.$x" "$OUT/ref/Homo_sapiens_assembly38.fasta.$x"; done
R=$OUT/ref/Homo_sapiens_assembly38.fasta
PB=(docker run --rm --gpus all -v "$O:$O" -w "$OUT" "$IMG" pbrun)
LOG=$OUT/stages.tsv
[ -s "$LOG" ] || echo -e "sample\tstage\tgpus\twall_s" > "$LOG"
stage() { # sample stage cmd...
  local s=$1 st=$2; shift 2
  echo "== $(date +%T) $s $st" >&2
  local t0; t0=$(date +%s.%N)
  "$@" > "$OUT/$s.$st.log" 2>&1
  printf "%s\t%s\t%s\t%.1f\n" "$s" "$st" "$NGPU" "$(echo "$(date +%s.%N) - $t0" | bc -l)" >> "$LOG"
}
for s in $SAMPLES; do
  case $s in
    tier1) FQ1=$O/data/fastq/SRR14724482_1.fastq.gz; FQ2=$O/data/fastq/SRR14724482_2.fastq.gz
           BED=$O/data/Agilent_v7.bed; ID=HG002_WES ;;
    wgs20) FQ1=$O/data/wgs_chr20/HG002_WGS_chr20_1.fastq.gz; FQ2=$O/data/wgs_chr20/HG002_WGS_chr20_2.fastq.gz
           BED=$REF/Annotation/intervals/wgs_calling_regions_noseconds.hg38.bed; ID=HG002_WGS20 ;;
    wgs30) FQ1=$O/data/sra_wgs/fastq/SRR14724532_1.fastq.gz; FQ2=$O/data/sra_wgs/fastq/SRR14724532_2.fastq.gz
           BED=$REF/Annotation/intervals/wgs_calling_regions_noseconds.hg38.bed; ID=HG002_WGS30 ;;
  esac
  [ -s "$OUT/$s.vcf.gz" ] && { echo "== $s done, skipping" >&2; continue; }
  RG='@RG\tID:'$ID'.L1\tPU:L1\tSM:HG002_'$ID'\tLB:'$ID'\tDS:'$FASTA'\tPL:ILLUMINA'
  stage "$s" fq2bam "${PB[@]}" fq2bam --ref "$R" --in-fq "$FQ1" "$FQ2" "$RG" --bwa-options="-K 100000000 -Y" \
    --knownSites "$B/dbsnp_146.hg38.vcf.gz" --knownSites "$B/beta/Homo_sapiens_assembly38.known_indels.vcf.gz" \
    --knownSites "$B/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" --interval-file "$BED" \
    --out-bam "$OUT/$s.md.bam" --out-recal-file "$OUT/$s.recal.txt" \
    --out-duplicate-metrics "$OUT/$s.dup_metrics.txt" --tmp-dir "$OUT"
  stage "$s" haplotypecaller "${PB[@]}" haplotypecaller --ref "$R" --in-bam "$OUT/$s.md.bam" \
    --in-recal-file "$OUT/$s.recal.txt" --interval-file "$BED" \
    --haplotypecaller-options="-pcr-indel-model CONSERVATIVE" \
    --out-variants "$OUT/$s.vcf.gz" --tmp-dir "$OUT"
done
column -t "$LOG"
