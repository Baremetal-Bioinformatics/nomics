#!/usr/bin/env bash
# Promote completed Mutect2 tasks, preserving staged-input provenance.
# Usage: somatic_harvest.sh RUN_DIR [--force]
set -euo pipefail
python3 -I "$(dirname "$0")/somatic_harvest.py" "$@"
