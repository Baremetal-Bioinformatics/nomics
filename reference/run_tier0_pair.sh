#!/usr/bin/env bash
# Two independent tier0 reference runs (same pinned sarek and params as tier1).
set -uo pipefail
D=$HOME/omics/runs; S=$HOME/omics/new_omics/reference; TAG=${1:?usage: run_tier0_pair.sh <YYYYMMDD>}
for r in r1 r2; do
  $S/run_sarek.sh $D/tier0-wes-$TAG-$r $S/samplesheet.tier0.csv $S/params.tier1.yaml > $D/tier0-wes-$TAG-$r.log 2>&1
done
echo PAIR-DONE
