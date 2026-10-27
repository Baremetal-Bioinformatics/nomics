#!/usr/bin/env bash
# Two sarek 3.10.0 reference runs of the tier1 exome with --joint_germline (HC in GVCF mode,
# GenomicsDBImport, GenotypeGVCFs) on the x86 cloud box, HaplotypeCaller in Java mode.
#   run_tier1_joint_cloud_pair.sh <YYYYMMDD>
set -uo pipefail
TAG=${1:?usage: run_tier1_joint_cloud_pair.sh <YYYYMMDD>}
S=/data/omics/new_omics/reference
export NXF_VER=26.04.6 NXF_HOME=/data/omics/.nextflow
for r in r1 r2; do
  RUN=/data/omics/runs/tier1-joint-x86-$TAG-$r
  mkdir -p "$RUN" && cd "$RUN"
  cp $S/samplesheet.cloud-tier1.csv samplesheet.csv
  cp $S/params.cloud-tier1-joint.yaml params.yaml
  date -Is > started_at
  # Concurrent first pulls of one image can fail a task (exit 125): retry with -resume.
  for attempt in 1 2 3; do
    RESUME=; [ $attempt -gt 1 ] && RESUME=-resume && mv -f report.html report.$attempt.html 2>/dev/null && mv -f timeline.html timeline.$attempt.html 2>/dev/null
    /usr/bin/time -v nextflow run /data/omics/ref/sarek-3.10.0 -profile docker $RESUME \
      -c $S/cloud.config -c ${EXTRA_CONFIG:-$S/cloud-32.config} -c $S/cloud-hcjava.config -params-file params.yaml \
      -with-trace -with-report report.html -with-timeline timeline.html \
      > nextflow.$attempt.out 2> time.$attempt.txt
    rc=$?; [ $rc -eq 0 ] && break
  done
  echo EXIT=$rc attempts=$attempt | tee exit_status
  date -Is > finished_at
done
echo PAIR-DONE
