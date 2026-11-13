#!/usr/bin/env bash
# Make FilterMutectCalls test cases beyond the one sarek golden: subsets of the golden's raw
# Mutect2 calls and altered contamination estimates, each filtered by GATK 4.6.2.0 in sarek's
# container. Every case changes what the tool learns (cluster priors, thresholds), so a port
# fitted to the golden alone does not pass them. Needs docker; on the Spark, run it under the
# resource lock. The gate is scripts/somatic_filter_gate.sh.
#
#   scripts/somatic_filter_cases.sh [OUT_DIR] [GOLDEN]
#   defaults: OUT_DIR=~/omics/runs/somatic/fmc-cases, GOLDEN=~/omics/runs/somatic/S1-r1/golden
#
# OUT_DIR/cases.tsv lists each case: name, raw calls, contamination table.
set -uo pipefail
O=${1:-$HOME/omics/runs/somatic/fmc-cases}; G=${2:-$HOME/omics/runs/somatic/S1-r1/golden}
X=HCC1395_WES_vs_HCC1395BL_WES
IMG=${GATK_IMAGE:-community.wave.seqera.io/library/gatk4_htslib_samtools:6fd46722be032bdd}
REF=${REF:-$HOME/omics/ref/GATK.GRCh38}/Sequence/WholeGenomeFasta
mkdir -p "$O"; cd "$O" || exit 2
cp "$G/MERGE_MUTECT2-$X/$X.mutect2.vcf.gz" raw.vcf.gz
cp "$G/MERGEMUTECTSTATS-$X/$X.mutect2.vcf.gz.stats" raw.vcf.gz.stats
cp "$G/LEARNREADORIENTATIONMODEL-$X/$X.mutect2.artifactprior.tar.gz" priors.tar.gz
cp "$G/CALCULATECONTAMINATION-$X/$X.mutect2.segmentation.table" seg.table
cp "$G/CALCULATECONTAMINATION-$X/$X.mutect2.contamination.table" contam.table
for c in 0.02 0.1 0.3; do
  awk -v c=$c 'BEGIN{FS=OFS="\t"} NR==1{print; next} {$2=c; print}' contam.table > contam$c.table
done
run() { docker run --rm -u "$(id -u):$(id -g)" -v "$O":/w -v "$REF":/ref:ro -w /w "$IMG" bash -c "$1"; }

# Subsets of the records: alternate ones, thirds, contig ranges, random 70% samples, the head.
run 'set -e
     sub() { zcat raw.vcf.gz | awk "/^#/ {print; next} $2" | bgzip > $1.vcf.gz; }
     sub odd "{n++; if (n % 2) print}"
     for k in 0 1 2; do zcat raw.vcf.gz | awk -v k=$k "/^#/ {print; next} {n++; if (n % 3 == k) print}" | bgzip > mod3_$k.vcf.gz; done
     sub chr1to7 "\$1 ~ /^chr[1-7]\$/"
     sub chr8to14 "\$1 ~ /^chr(8|9|10|11|12|13|14)\$/"
     sub chr15on "\$1 ~ /^chr(15|16|17|18|19|20|21|22|X|Y)\$/"
     sub rest "\$1 !~ /^chr[1-7]\$/"
     for s in 11 23 47; do zcat raw.vcf.gz | awk -v s=$s "BEGIN{srand(s)} /^#/ {print; next} rand() < 0.7" | bgzip > rand$s.vcf.gz; done
     sub first500 "{n++; if (n <= 500) print}"
     for v in raw odd mod3_0 mod3_1 mod3_2 chr1to7 chr8to14 chr15on rest rand11 rand23 rand47 first500; do
       [ $v = raw ] || cp raw.vcf.gz.stats $v.vcf.gz.stats; tabix -f -p vcf $v.vcf.gz
     done' || { echo "making the subsets failed"; exit 1; }

: > cases.tsv
fail=0
filter() {  # filter NAME VARIANT CONTAM
  run "gatk --java-options '-Xmx8g -XX:-UsePerfData' FilterMutectCalls --variant $2.vcf.gz --output $1.filtered.vcf.gz \
       --reference /ref/Homo_sapiens_assembly38.fasta --orientation-bias-artifact-priors priors.tar.gz \
       --tumor-segmentation seg.table --contamination-table $3 --tmp-dir . > $1.log 2>&1" || { echo "$1: GATK failed (see $O/$1.log)"; fail=1; return; }
  printf '%s\t%s\t%s\n' "$1" "$2" "$3" >> cases.tsv
  echo "$1: $(zcat "$1.filtered.vcf.gz" | grep -vc '^#') records, $(zcat "$1.filtered.vcf.gz" | awk -F'\t' '!/^#/ && $7=="PASS"' | wc -l) PASS"
}
for v in odd mod3_0 mod3_1 mod3_2 chr1to7 chr8to14 chr15on rest rand11 rand23 rand47 first500; do filter "$v" "$v" contam.table; done
for c in 0.02 0.1 0.3; do filter "contam$c" raw "contam$c.table"; done
filter rand11c0.3 rand11 contam0.3.table
filter same raw contam.table
# The container must reproduce the sarek golden on the unaltered input, or the cases mean nothing.
if cmp -s <(zcat same.filtered.vcf.gz | grep -v '^#') <(zcat "$G/FILTERMUTECTCALLS-$X/$X.mutect2.filtered.vcf.gz" | grep -v '^#') \
   && cmp -s same.filtered.vcf.gz.filteringStats.tsv "$G/FILTERMUTECTCALLS-$X/$X.mutect2.filtered.vcf.gz.filteringStats.tsv"
then echo "same: this container reproduces the golden's records and statistics"
else echo "same: this container does NOT reproduce the golden"; fail=1; fi
echo "$(wc -l < cases.tsv) cases in $O"
exit $fail
