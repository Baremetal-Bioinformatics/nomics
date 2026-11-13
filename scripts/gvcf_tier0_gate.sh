#!/usr/bin/env bash
# Reproduce sarek's single grouped-BED tier0 HC invocation, then compare every
# record and header line under docs/GATES.md. No reference tools are invoked.
set -euo pipefail
repo=$(cd "$(dirname "$0")/.." && pwd)
fixture=${1:-"$HOME/omics-fixtures/gvcf-tier0"}
output=${2:-"${TMPDIR:-/tmp}/nomics-gvcf-tier0.haplotypecaller.g.vcf.gz"}
if [[ ! -f "$fixture/REF-READY" ]]; then
    echo "Reference copy is incomplete: waiting for $fixture/REF-READY" >&2
    while [[ ! -f "$fixture/REF-READY" ]]; do sleep 5; done
fi
cd "$repo"
cargo build --release -p nomics
binary=${CARGO_TARGET_DIR:-target}/release/nomics
"$binary" hc --emit-gvcf \
    --input "$fixture/recal/HG002_WES.recal.bam" \
    --intervals "$fixture/Agilent_v7.bed" \
    --reference "$fixture/ref/Homo_sapiens_assembly38.fasta" \
    --dbsnp "$fixture/ref/dbsnp_146.hg38.vcf.gz" \
    --output "$output" --threads "${NOMICS_THREADS:-8}"
"$binary" dev gate-gvcf --ours "$output" \
    --golden "$fixture/golden/HG002_WES.haplotypecaller.g.vcf.gz"
