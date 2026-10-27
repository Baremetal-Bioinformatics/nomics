#!/usr/bin/env bash
# gate_stage.sh EXAMPLE TAGS... : run a nomics-align dump example on every tier0
# chunk and compare each tag's lines with the oracle dumps (~/omics/oracle/dumps).
set -uo pipefail
EX=$1; shift
cd ~/omics
IDX=$HOME/omics/ref/GATK.GRCh38/Sequence/BWAIndex/Homo_sapiens_assembly38.fasta
for o in oracle/dumps/*/; do
  c=$(basename $o); wd=$(cat $o/workdir); out=oracle/ours/$c; mkdir -p $out
  new_omics/target/release/examples/$EX $IDX $wd/$c.*_R1.fastp.fastq.gz $wd/$c.*_R2.fastp.fastq.gz $out/x
  line="$c"
  for t in "$@"; do
    n=$(wc -l < $o/$t.txt)
    if cmp -s $out/x.$t.txt $o/$t.txt; then line="$line $t:IDENTICAL($n)"
    else line="$line $t:DIFF($(diff $out/x.$t.txt $o/$t.txt | grep -c '^>')/$n)"; fi
  done
  echo "$line"
done
