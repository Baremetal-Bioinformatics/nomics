#!/usr/bin/env bash
# Run the oracle bwa (a source build of v0.7.19 with the golden binary's flags,
# optionally instrumented) on one golden chunk, sort with the golden's samtools
# image, and gate against the golden chunk BAM.
#
#   run_bwa_oracle.sh BWA_BINARY CHUNK_WORKDIR GOLDEN_BAM OUT_BAM [extra bwa env...]
#
# Test oracle only: never part of the nomics toolchain (see docs/GATES.md).
set -euo pipefail
BWA=$1 WD=$2 GOLDEN=$3 OUT=$4
IMG=community.wave.seqera.io/library/bwa_htslib_samtools:73cb446d71839fbc
INDEX=$HOME/omics/ref/GATK.GRCh38/Sequence/BWAIndex/Homo_sapiens_assembly38.fasta.64
R1=$(ls "$WD"/*_R1.fastp.fastq.gz) R2=$(ls "$WD"/*_R2.fastp.fastq.gz)
"$BWA" mem -K 100000000 -Y \
  -R "@RG\tID:HG002_WES.L1\tPU:L1\tSM:HG002_HG002_WES\tLB:HG002_WES\tDS:${NOMICS_RG_DS:-${NOMICS_REF:-$HOME/omics/ref/GATK.GRCh38}/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta}\tPL:ILLUMINA" \
  -t 20 "$INDEX" "$R1" "$R2" 2>"${OUT%.bam}.bwa.log" \
  | docker run --rm -i -u "$(id -u):$(id -g)" -v "$(dirname "$(readlink -f "$OUT")")":/out "$IMG" \
      samtools sort --threads 20 -o "/out/$(basename "$OUT")" -
"$HOME/omics/new_omics/target/release/nomics" dev gate-bam --ours "$OUT" --golden "$GOLDEN"
