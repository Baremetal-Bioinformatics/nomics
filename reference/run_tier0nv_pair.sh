#!/usr/bin/env bash
# Two reference runs of tier0nv: tier0 re-headered to NovaSeq-style names (fastp
# poly-G trimming and Picard optical duplicates switch on) and split into two lanes.
# Input built by reference/make_tier0nv.py. Same pinned sarek and params as tier1.
set -uo pipefail
D=$HOME/omics/runs; S=$HOME/omics/new_omics/reference; TAG=${1:?usage: run_tier0nv_pair.sh <YYYYMMDD>}
for r in r1 r2; do
  $S/run_sarek.sh $D/tier0nv-wes-$TAG-$r $S/samplesheet.tier0nv.csv $S/params.tier1.yaml > $D/tier0nv-wes-$TAG-$r.log 2>&1
done
echo PAIR-DONE
