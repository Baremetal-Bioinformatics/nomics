#!/usr/bin/env bash
# Rebuild sarek 3.10.0's default (save_output_as_bam false) md.cram and
# recal.cram for the tier0 golden, with the containers the golden run used,
# so the CRAM gate can compare against what a default sarek user gets.
#
#   scripts/cram_reference.sh [golden-run-dir] [out-dir]
#
# Chain, from sarek 3.10.0 source:
# - modules/nf-core/gatk4/markduplicates/main.nf: MarkDuplicates writes
#   md.bam, then `samtools view -Ch -T <fasta> -o <id>.md.cram <id>.md.bam`
#   and `samtools index <id>.md.cram` (samtools 1.21 in the GATK container).
#   CRAM mode differs only in dropping --CREATE_INDEX, so the golden md.bam
#   is that BAM.
# - subworkflows/local/bam_applybqsr + conf/modules/recalibrate.config:
#   GATK4_APPLYBQSR reads md.cram and writes <id>.recal.cram itself (htsjdk
#   CRAM writer) per interval group. tier0 has one group (the golden task
#   wrote <id>.recal.bam, i.e. num_intervals <= 1), so CRAM_MERGE_INDEX_SAMTOOLS
#   skips MERGE_CRAM (which would be `samtools merge --output-fmt
#   cram,version=3.0`) and only runs `samtools index` (samtools 1.24).
# BaseRecalibrator's table and the interval BED are taken from the golden.
set -euo pipefail
G=${1:-${NOMICS_FIXTURES:-$HOME/omics-fixtures}/sarek-tier0/r1}
OUT=${2:-$HOME/omics/ref-cram/tier0}
ID=HG002_WES
REFDIR=$HOME/omics/ref/GATK.GRCh38/Sequence/WholeGenomeFasta
GATK_IMG=community.wave.seqera.io/library/gatk4_htslib_samtools:6fd46722be032bdd # GATK 4.6.2.0 / htsjdk 4.2.0, samtools 1.21
SAM_IMG=community.wave.seqera.io/library/htslib_samtools:1.24--1523e39b702f9721
mkdir -p "$OUT"
OUT=$(cd "$OUT" && pwd)
work=$G/work
md_bam=$G/results/preprocessing/markduplicates/$ID/$ID.md.bam
table=$(ls "$work"/5d/*/"$ID".recal.table)
bed=$(ls "$work"/af/*/chr1_12146-12310.bed)
# Stage inputs as Nextflow does (symlinks under their task names).
for f in Homo_sapiens_assembly38.fasta Homo_sapiens_assembly38.fasta.fai Homo_sapiens_assembly38.dict; do
  ln -sf "$REFDIR/$f" "$OUT/$f"
done
ln -sf "$md_bam" "$OUT/$ID.md.bam"
ln -sf "$table" "$OUT/$ID.recal.table"
ln -sf "$bed" "$OUT/chr1_12146-12310.bed"
D=(docker run --rm --platform linux/arm64 -u "$(id -u):$(id -g)" -v "$HOME/omics:$HOME/omics"
   -v "$G:$G" -v "$OUT:$OUT" -w "$OUT")
"${D[@]}" "$GATK_IMG" bash -c "
  set -euo pipefail
  samtools view -Ch -T Homo_sapiens_assembly38.fasta -o $ID.md.cram $ID.md.bam
  samtools index $ID.md.cram"
"${D[@]}" "$GATK_IMG" bash -c "
  set -euo pipefail
  gatk --java-options '-Xmx3276M -XX:-UsePerfData' ApplyBQSR \
    --input $ID.md.cram --output $ID.recal.cram --reference Homo_sapiens_assembly38.fasta \
    --bqsr-recal-file $ID.recal.table --intervals chr1_12146-12310.bed --tmp-dir ." 2>"$OUT/applybqsr.log"
"${D[@]}" "$SAM_IMG" samtools index "$ID.recal.cram"
for f in "$ID.md.cram" "$ID.recal.cram"; do
  v=$(od -An -tu1 -j4 -N2 "$OUT/$f" | awk 'NF{print $1"."$2}')
  echo "$f: CRAM $v, $(wc -c <"$OUT/$f" | tr -d ' ') bytes"
done
