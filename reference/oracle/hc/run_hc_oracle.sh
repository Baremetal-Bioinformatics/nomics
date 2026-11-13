#!/usr/bin/env bash
# Re-run the reference HaplotypeCaller (unmodified GATK 4.6.2.0 jar, same
# container, same arguments as the golden task) with extra debug outputs, and
# check its VCF records still equal the golden's.
#
#   run_hc_oracle.sh <tier0|tier1> <outdir> [extra gatk args...]
#
# Test oracle only (docs/GATES.md); outputs land in <outdir>.
set -euo pipefail
tier=$1; out=$2; shift 2
run=$HOME/omics/runs/$tier-wes-20261008-r1
wd=$(dirname "$(grep -l "HaplotypeCaller" "$run"/work/*/*/.command.sh | head -1)")
IMG=community.wave.seqera.io/library/gatk4_htslib_samtools:6fd46722be032bdd
mkdir -p "$out"
out=$(readlink -f "$out")
# Stage the task's inputs (symlinks resolved) as in the work dir.
for f in "$wd"/*; do ln -sf "$(readlink -f "$f")" "$out/$(basename "$f")"; done
rm -f "$out/HG002_WES.haplotypecaller.vcf.gz" "$out/HG002_WES.haplotypecaller.vcf.gz.tbi"
mounts=(-v "$out":"$out" -v "$HOME/omics":"$HOME/omics")
docker run --rm --pull never -u "$(id -u):$(id -g)" "${mounts[@]}" -w "$out" "$IMG" \
  gatk --java-options "-Xmx9830M -XX:-UsePerfData" HaplotypeCaller \
    --input HG002_WES.recal.bam --output HG002_WES.haplotypecaller.vcf.gz \
    --reference Homo_sapiens_assembly38.fasta --native-pair-hmm-threads 2 \
    --dbsnp dbsnp_146.hg38.vcf.gz --intervals chr1_12146-12310.bed \
    --tmp-dir . --pcr-indel-model CONSERVATIVE "$@" > "$out/gatk.log" 2>&1
golden=$(ls "$run"/results/variant_calling/haplotypecaller/HG002_WES/HG002_WES.haplotypecaller.vcf.gz)
if cmp -s <(zcat "$golden" | grep -v '^##GATKCommandLine') <(zcat "$out/HG002_WES.haplotypecaller.vcf.gz" | grep -v '^##GATKCommandLine'); then
  echo "oracle VCF: IDENTICAL to golden (excl. ##GATKCommandLine)"
else
  echo "oracle VCF: DIFFERS from golden"
fi
