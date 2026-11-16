#!/usr/bin/env bash
# FASTQs -> both preprocessing boundaries -> 16 somatic comparisons. Linux/GNU time.
# Run under the gate machine's resource lock.
# Missing fixtures exit 77 (SKIP), or 1 with NOMICS_REQUIRE_FIXTURES=1; never PASS.
set -euo pipefail
REPO=$(cd "$(dirname "$0")/.." && pwd)
source "$REPO/scripts/fixture_policy.sh"
R=${1:?usage: somatic_run_gate.sh RUN_DIR OUT_DIR}
O=${2:?usage: somatic_run_gate.sh RUN_DIR OUT_DIR}
for fixture in "$R/samplesheet.csv" "$R/params.yaml" "$R/work" "$R/golden/intervals"; do
  if [[ ! -e $fixture ]]; then fixture_skip "$fixture" || exit 1; exit 77; fi
done
export NOMICS=${NOMICS:-$REPO/target/release/nomics}
python3 -I "$REPO/scripts/somatic_run_gate.py" "$R" "$O"
