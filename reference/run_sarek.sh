#!/usr/bin/env bash
# Usage: run_sarek.sh <run_dir> <samplesheet> <params.yaml>
# Pinned reference run: sarek 3.10.0 @ 8ccac7ad, Nextflow 26.04.6, arm64 Wave images.
set -uo pipefail
RUN=$1; SHEET=$2; PARAMS=$3
export PATH=$HOME/.local/bin:$PATH JAVA_HOME=$HOME/.local/opt/jdk17 NXF_VER=26.04.6
mkdir -p "$RUN" && cd "$RUN"
envsubst < "$SHEET" > samplesheet.csv   # the samplesheets use ${HOME}
exec 9>${NOMICS_LOCK:-$HOME/.cache/uni-rnaseq-resource.lock}
flock 9
date -Is > started_at
envsubst < "$PARAMS" > params.yaml
/usr/bin/time -v nextflow run $HOME/omics/ref/sarek-3.10.0 -profile docker,arm64 \
  -c $HOME/omics/new_omics/reference/spark.config -params-file params.yaml \
  -with-trace -with-report report.html -with-timeline timeline.html \
  2> time.txt
echo EXIT=$? | tee exit_status
date -Is > finished_at
