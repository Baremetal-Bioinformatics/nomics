#!/usr/bin/env bash
# MergeVcfs and MergeMutectStats against GATK 4.6.2.0 on inputs the port was not written
# against: the Mutect2 slice goldens (scripts/somatic_mutect2_slice.sh) in a scrambled order,
# a pair of them, and a set with a records-free VCF. GATK runs here, in sarek's container, so
# this needs docker; the golden itself is the test `merge_goldens`.
#
#   scripts/somatic_merge_cases.sh [OUT]     NOMICS=path/to/nomics, M2_GOLDENS (~/omics/runs/somatic),
#                                            REF (~/omics/ref/GATK.GRCh38)
set -uo pipefail
REPO=$(cd "$(dirname "$0")/.." && pwd)
bin=${NOMICS:-$REPO/target/release/nomics}
R=${M2_GOLDENS:-$HOME/omics/runs/somatic}; O=${1:-$R/merge-cases}; rm -rf "$O"; mkdir -p "$O"; cd "$O"
GK=community.wave.seqera.io/library/gatk4_htslib_samtools:6fd46722be032bdd; U="-u $(id -u):$(id -g)"
DICT=${REF:-$HOME/omics/ref/GATK.GRCh38}/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.dict
v() { ls $R/m2-$1/oracle/*.mutect2.$1.vcf.gz; }
# a VCF with the header only
zcat $(v chr21) | grep '^#' | docker run --rm -i $U -v $O:/w -w /w $GK bgzip -c > empty.vcf.gz
case1="chr19 chr21 chr12 chr13 chr22 chr20 chr16 chr17 chr11"; case2="chr22 chr13"; case3="chr20 EMPTY chr12"
pass=0; fail=0
for n in 1 2 3; do eval "set -- \$case$n"; ins=(); st=()
  for c in "$@"; do if [ $c = EMPTY ]; then ins+=($O/empty.vcf.gz); else ins+=($(v $c)); st+=($(v $c).stats); fi; done
  docker run --rm $U -v $HOME/omics:$HOME/omics -w $O $GK gatk --java-options "-Xmx3276M -XX:-UsePerfData" MergeVcfs $(printf -- '--INPUT %s ' "${ins[@]}") --OUTPUT $O/gatk$n.vcf.gz --SEQUENCE_DICTIONARY $DICT --TMP_DIR . > gatk$n.log 2>&1 || echo "gatk MergeVcfs $n failed"
  docker run --rm $U -v $HOME/omics:$HOME/omics -w $O $GK gatk --java-options "-Xmx3276M -XX:-UsePerfData" MergeMutectStats $(printf -- '--stats %s ' "${st[@]}") --output $O/gatk$n.stats --tmp-dir . > gatks$n.log 2>&1 || echo "gatk MergeMutectStats $n failed"
  $bin somatic merge-vcfs $(printf -- '--input %s ' "${ins[@]}") --sequence-dictionary $DICT --output $O/ours$n.vcf.gz 2> ours$n.log || echo "ours merge-vcfs $n failed: $(cat ours$n.log)"
  $bin somatic merge-mutect-stats $(printf -- '--stats %s ' "${st[@]}") --output $O/ours$n.stats || echo "ours stats $n failed"
  if cmp <(zcat gatk$n.vcf.gz) <(zcat ours$n.vcf.gz); then echo "PASS case $n VCF: $(zcat ours$n.vcf.gz | grep -vc '^#') records, every decompressed byte ($*)"; pass=$((pass+1)); else echo "FAIL case $n VCF"; fail=$((fail+1)); fi
  if cmp gatk$n.stats ours$n.stats; then echo "PASS case $n stats: $(tail -1 ours$n.stats)"; pass=$((pass+1)); else echo "FAIL case $n stats: $(tail -1 gatk$n.stats) vs $(tail -1 ours$n.stats)"; fail=$((fail+1)); fi
done
# what Picard does where the port refuses: the same file twice
docker run --rm $U -v $HOME/omics:$HOME/omics -w $O $GK gatk MergeVcfs --INPUT $(v chr21) --INPUT $(v chr21) --OUTPUT $O/twice.vcf.gz --TMP_DIR . > twice.log 2>&1; echo "GATK on the same file twice: exit $?, $(zcat twice.vcf.gz 2>/dev/null | grep -vc '^#') records"
$bin somatic merge-vcfs --input $(v chr21) --input $(v chr21) --output $O/twice-ours.vcf.gz; echo "ours on the same file twice: exit $?"
echo "merge cases: $pass pass, $fail fail"
(( fail == 0 && pass > 0 ))
