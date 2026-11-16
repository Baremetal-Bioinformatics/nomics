#!/usr/bin/env bash
# Usage: run_sarek_cloud.sh <run_dir> <samplesheet> <params.yaml>
# Pinned reference run on the x86 cloud box: sarek 3.10.0 @ 8ccac7ad (local copy),
# Nextflow 26.04.6, default (amd64) containers. Records start/finish times: samtools
# merge seeds its @RG/@PG ID suffixes with the wall clock (REFERENCE-VARIABILITY.md).
set -uo pipefail
RUN=$1; SHEET=$2; PARAMS=$3
export NXF_VER=26.04.6 NXF_HOME=/data/omics/.nextflow
mkdir -p "$RUN" && cd "$RUN"
cp "$SHEET" samplesheet.csv
cp "$PARAMS" params.yaml
date -Is > started_at
/usr/bin/time -v nextflow run /data/omics/ref/sarek-3.10.0 -profile docker \
  -c /data/omics/new_omics/reference/cloud.config -params-file params.yaml \
  -with-trace -with-report report.html -with-timeline timeline.html \
  > nextflow.out 2> time.txt
echo EXIT=$? | tee exit_status
date -Is > finished_at
