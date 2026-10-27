#!/usr/bin/env bash
# Two sarek 3.10.0 reference runs of the 30x slice (first 1/12 of the raw
# SRR14724532 pairs) on the x86 cloud box, HaplotypeCaller in Java mode.
#   run_wgs30s_cloud_pair.sh <YYYYMMDD>
set -uo pipefail
TAG=${1:?usage: run_wgs30s_cloud_pair.sh <YYYYMMDD>}
S=/data/omics/new_omics/reference
export NXF_VER=26.04.6 NXF_HOME=/data/omics/.nextflow
for r in r1 r2; do
  RUN=/data/omics/runs/wgs30s-x86-$TAG-$r
  mkdir -p "$RUN" && cd "$RUN"
  cp $S/samplesheet.cloud-wgs30s.csv samplesheet.csv
  cp $S/params.cloud-wgs.yaml params.yaml
  date -Is > started_at
  /usr/bin/time -v nextflow run /data/omics/ref/sarek-3.10.0 -profile docker \
    -c $S/cloud.config -c ${EXTRA_CONFIG:-$S/cloud-96.config} -c $S/cloud-hcjava.config -params-file params.yaml \
    -with-trace -with-report report.html -with-timeline timeline.html \
    > nextflow.out 2> time.txt
  echo EXIT=$? | tee exit_status
  date -Is > finished_at
done
echo PAIR-DONE
