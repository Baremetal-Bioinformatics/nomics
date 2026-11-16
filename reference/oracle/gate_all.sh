#!/usr/bin/env bash
# gate_all.sh BWA TAG : run oracle on all tier0 chunks, print verdicts
cd ~/omics
for w in $(grep -l "bwa mem" runs/tier0-wes-20261008-r1/work/*/*/.command.sh); do
  d=$(dirname $w); c=$(basename $(ls $d/*.bam) .bam | sed "s/.*\.//")
  mkdir -p oracle/out/$2
  v=$(new_omics/reference/oracle/run_bwa_oracle.sh $1 $d runs/tier0-wes-20261008-r1/results/preprocessing/mapped/HG002_WES/HG002_WES.$c.bam oracle/out/$2/$c.bam 2>&1 | grep -E "^records|verdict" | paste -sd" " -)
  echo "$2 $c $v"
done
