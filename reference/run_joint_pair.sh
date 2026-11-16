#!/usr/bin/env bash
# Two independent --joint_germline reference runs for one tier (sarek 3.10.0, tier1 params + joint_germline).
#   run_joint_pair.sh <tier0|tier1> <YYYYMMDD>
set -uo pipefail
D=$HOME/omics/runs; S=$HOME/omics/new_omics/reference
T=${1:?usage: run_joint_pair.sh <tier0|tier1> <YYYYMMDD>}; TAG=${2:?}
for r in r1 r2; do
  $S/run_sarek.sh $D/$T-joint-$TAG-$r $S/samplesheet.$T.csv $S/params.tier1-joint.yaml > $D/$T-joint-$TAG-$r.log 2>&1
done
echo PAIR-DONE
