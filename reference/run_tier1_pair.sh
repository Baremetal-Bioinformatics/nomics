#!/usr/bin/env bash
# Two independent full-WES reference runs (r1, r2) to measure the reference's own run-to-run variability.
set -uo pipefail
D=$HOME/omics/runs; S=$HOME/omics/new_omics/reference
for r in r1 r2; do
  $S/run_sarek.sh $D/tier1-wes-20261008-$r $S/samplesheet.tier1.csv $S/params.tier1.yaml > $D/tier1-wes-20261008-$r.log 2>&1
done
echo PAIR-DONE
