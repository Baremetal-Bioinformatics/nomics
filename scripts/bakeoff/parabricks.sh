#!/usr/bin/env bash
# Parabricks 4.7.1 bake-off on the Spark, same inputs as sarek/nomics tier1.
#
#   scripts/bakeoff/parabricks.sh <outdir>      (under the resource lock, nohup)
#
# Two configurations, each fq2bam (BWA-MEM + sort + markdup + BQSR table) then
# haplotypecaller with the table applied:
#   raw    the ENA FASTQ pair (Parabricks' natural input)
#   fastp  the 12 fastp chunk pairs sarek aligned (same reads as the goldens)
# Then scripts/bakeoff/concordance.py against the sarek golden VCF and GIAB.
set -euo pipefail
OUT=${1:?outdir}
IMG=nvcr.io/nvidia/clara/clara-parabricks:4.7.1-1
REF=$HOME/omics/ref/GATK.GRCh38
FASTA=$REF/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta
BWA=$REF/Sequence/BWAIndex/Homo_sapiens_assembly38.fasta.64
B=$REF/Annotation/GATKBundle
BED=$HOME/omics/data/Agilent_v7.bed
GOLD=$HOME/omics/runs/tier1-wes-20261008-r1
RG='@RG\tID:HG002_WES.L1\tPU:L1\tSM:HG002_HG002_WES\tLB:HG002_WES\tDS:'$FASTA'\tPL:ILLUMINA'
mkdir -p "$OUT/ref"
# Parabricks wants <fasta>.{amb,ann,bwt,pac,sa} next to the real FASTA (it resolves
# symlinks); ours use bwa's .64 naming, so hard-link everything into one dir.
for x in fasta fasta.fai dict; do ln -f "${FASTA%.fasta}.$x" "$OUT/ref/Homo_sapiens_assembly38.$x"; done
for x in amb ann bwt pac sa; do ln -f "$BWA.$x" "$OUT/ref/Homo_sapiens_assembly38.fasta.$x"; done
R=$OUT/ref/Homo_sapiens_assembly38.fasta
PB=(docker run --rm --gpus all -u "$(id -u):$(id -g)" -v "$HOME/omics:$HOME/omics" -w "$OUT" "$IMG" pbrun)
LOG=$OUT/stages.tsv
[ -s "$LOG" ] || echo -e "config\tstage\twall_s\tpeak_rss_gb" > "$LOG"
stage() { # config stage cmd...
  local c=$1 s=$2; shift 2
  echo "== $(date +%T) $c $s" >&2
  /usr/bin/time -f "%e %M" -o "$OUT/.time" "$@" > "$OUT/$c.$s.log" 2>&1
  read -r wall rss < "$OUT/.time"
  printf "%s\t%s\t%s\t%.1f\n" "$c" "$s" "$wall" "$(echo "$rss / 1048576" | bc -l)" >> "$LOG"
}
fastp_dir=$(dirname "$(grep -l 'fastp' "$GOLD"/work/*/*/.command.sh | while read -r f; do d=$(dirname "$f"); [ "$(ls "$d" | grep -c 'fastp.fastq.gz')" -eq 24 ] && echo "$f" && break; done)")
# Parabricks rejects repeated read-group IDs across --in-fq pairs, so join the
# 12 chunks in order (gzip members concatenate): same reads, one read group.
for m in 1 2; do
  [ -s "$OUT/fastp_R$m.fastq.gz" ] || cat "$fastp_dir"/*_R$m.fastp.fastq.gz > "$OUT/fastp_R$m.fastq.gz"
done

for cfg in raw fastp; do
  case $cfg in
    raw) inputs=(--in-fq "$HOME/omics/data/fastq/SRR14724482_1.fastq.gz" "$HOME/omics/data/fastq/SRR14724482_2.fastq.gz" "$RG") ;;
    fastp) inputs=(--in-fq "$OUT/fastp_R1.fastq.gz" "$OUT/fastp_R2.fastq.gz" "$RG") ;;
  esac
  [ -s "$OUT/$cfg.vcf.gz" ] && { echo "== $cfg done, skipping" >&2; continue; }
  stage $cfg fq2bam "${PB[@]}" fq2bam --ref "$R" "${inputs[@]}" --bwa-options="-K 100000000 -Y" \
    --knownSites "$B/dbsnp_146.hg38.vcf.gz" --knownSites "$B/beta/Homo_sapiens_assembly38.known_indels.vcf.gz" \
    --knownSites "$B/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" --interval-file "$BED" \
    --out-bam "$OUT/$cfg.md.bam" --out-recal-file "$OUT/$cfg.recal.txt" \
    --out-duplicate-metrics "$OUT/$cfg.dup_metrics.txt" --tmp-dir "$OUT"
  stage $cfg haplotypecaller "${PB[@]}" haplotypecaller --ref "$R" --in-bam "$OUT/$cfg.md.bam" \
    --in-recal-file "$OUT/$cfg.recal.txt" --interval-file "$BED" \
    --haplotypecaller-options="-pcr-indel-model CONSERVATIVE" \
    --out-variants "$OUT/$cfg.vcf.gz" --tmp-dir "$OUT"
done
column -t "$LOG"
python3 "$(dirname "$0")/concordance.py" \
  --golden "$GOLD/results/variant_calling/haplotypecaller/HG002_WES/HG002_WES.haplotypecaller.vcf.gz" \
  --truth "$HOME/omics/data/truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" \
  --truth-bed "$HOME/omics/data/truth/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed" \
  --targets "$BED" --vcf raw="$OUT/raw.vcf.gz" --vcf fastp="$OUT/fastp.vcf.gz" | tee "$OUT/concordance.txt"
