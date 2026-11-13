#!/usr/bin/env bash
# Mutect2 slice goldens for the caller port: one contig of the HCC1395 tumour-normal pair's
# recalibrated BAMs, of the panel of normals, the germline resource and the intervals, run
# through GATK 4.6.2.0 Mutect2 in sarek's container with sarek's arguments plus GATK's own
# stage outputs (--assembly-region-out, --bam-output with every haplotype). Needs docker and
# the sarek run; on the Spark, under the resource lock. Output: ~/omics/runs/somatic/m2-CONTIG/
# (inputs, and oracle/ with the VCF, stats, f1r2 archive, regions.igv, haplotypes.bam, gatk.log).
# The gate for the first stage is scripts/somatic_mutect2_regions_gate.sh.
#
#   scripts/somatic_mutect2_slice.sh CONTIG ...
set -uo pipefail
R=$HOME/omics/runs/somatic/S1-r1; REF=$HOME/omics/ref/GATK.GRCh38; U="-u $(id -u):$(id -g)"
ST=community.wave.seqera.io/library/htslib_samtools:1.24--1523e39b702f9721
BT=community.wave.seqera.io/library/bcftools_htslib:1.23.1--583f22aacd68faff
GK=community.wave.seqera.io/library/gatk4_htslib_samtools:6fd46722be032bdd
G=$R/golden/MERGE_MUTECT2-HCC1395_WES_vs_HCC1395BL_WES/HCC1395_WES_vs_HCC1395BL_WES.mutect2.vcf.gz
for c in "$@"; do
  D=$HOME/omics/runs/somatic/m2-$c; mkdir -p "$D/oracle"; cd "$D"
  X=HCC1395_WES_vs_HCC1395BL_WES.mutect2.$c
  cat "$R"/golden/intervals/*.bed | grep -P "^$c\t" | sort -k2,2n > $c.bed
  for s in HCC1395_WES HCC1395BL_WES; do
    [ -f $s.$c.bam.bai ] || docker run --rm $U -v "$R/results/preprocessing/recalibrated/$s":/in:ro -v "$D":/w -w /w $ST \
      sh -c "samtools view -b -o $s.$c.bam /in/$s.recal.bam $c && samtools index $s.$c.bam"
  done
  docker run --rm $U -v "$REF/Annotation/GATKBundle":/in:ro -v "$D":/w -w /w $BT \
    sh -c "bcftools view -r $c -Oz -o af-only-gnomad.$c.vcf.gz /in/af-only-gnomad.hg38.vcf.gz && tabix -f -p vcf af-only-gnomad.$c.vcf.gz && bcftools view -r $c -Oz -o 1000g_pon.$c.vcf.gz /in/1000g_pon.hg38.vcf.gz && tabix -f -p vcf 1000g_pon.$c.vcf.gz"
  s=$(date +%s)
  docker run --rm $U -v "$REF":/ref:ro -v "$D":/w -w /w $GK bash -c "gatk --java-options '-Xmx9830M -XX:-UsePerfData' Mutect2 --input HCC1395BL_WES.$c.bam --input HCC1395_WES.$c.bam \
     --output oracle/$X.vcf.gz --reference /ref/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta \
     --panel-of-normals 1000g_pon.$c.vcf.gz --germline-resource af-only-gnomad.$c.vcf.gz \
     --intervals $c.bed --tmp-dir . --f1r2-tar-gz oracle/$X.f1r2.tar.gz --normal-sample HCC1395_HCC1395BL_WES \
     --assembly-region-out oracle/regions.igv --bam-output oracle/haplotypes.bam --bam-writer-type ALL_POSSIBLE_HAPLOTYPES > oracle/gatk.log 2>&1"
  echo "$c: gatk exit $?, $(( $(date +%s) - s )) s; $(zcat oracle/$X.vcf.gz | grep -vc '^#') records; $(grep -c size= oracle/regions.igv) regions; $(cat oracle/$X.vcf.gz.stats | tail -n 1)"
  echo "$c: full run has $(zcat "$G" | awk -v c=$c '$1==c' | wc -l) records there; identical lines in the slice run: $(comm -12 <(zcat "$G" | awk -v c=$c '$1==c' | sort) <(zcat oracle/$X.vcf.gz | grep -v '^#' | sort) | wc -l)"
done
echo M2SLICES_DONE
