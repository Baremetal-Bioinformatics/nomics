#!/usr/bin/env bash
# Gate the Mutect2 chain's auxiliary tools end to end against a sarek tumour-normal run:
# GetPileupSummaries from the run's recalibrated BAMs (one table per interval group, read-filter
# counts against GATK's log), GatherPileupSummaries, CalculateContamination and
# LearnReadOrientationModel, each output compared by bytes with the golden (docs/GATES.md).
# RUN is the sarek run directory with results/ and the promoted golden/ task directories.
# On the Spark, run it under the resource lock.
#
#   scripts/somatic_boundary1_gate.sh [RUN] [OUT_DIR]     NOMICS=path/to/nomics to pick the binary
#   defaults: RUN=~/omics/runs/somatic/S1-r1, reference bundle REF=~/omics/ref/GATK.GRCh38
set -uo pipefail
REPO=$(cd "$(dirname "$0")/.." && pwd)
bin=${NOMICS:-$REPO/target/release/nomics}
R=${1:-$HOME/omics/runs/somatic/S1-r1}; G=$R/golden
O=${2:-$(mktemp -d "${TMPDIR:-/tmp}/somatic-boundary1.XXXXXX")}; mkdir -p "$O"
REF=${REF:-$HOME/omics/ref/GATK.GRCh38}
tumour=HCC1395_WES; normal=HCC1395BL_WES; P=${tumour}_vs_${normal}.mutect2
pass=0; fail=0
chk() { if cmp -s "$1" "$2"; then pass=$((pass+1)); echo "  identical: $3"; else fail=$((fail+1)); echo "  DIFFERS: $3"; fi; }

for p in "TUMOR $tumour" "NORMAL $normal"; do set -- $p
  tables=()
  for bed in "$G"/intervals/*.bed; do
    g=$(basename "$bed" .bed); t=$O/$2.mutect2.$g.pileups.table; s=$(date +%s)
    golden=$(command ls "$G"/GETPILEUPSUMMARIES_$1-$2*/"$2.mutect2.$g.pileups.table")
    "$bin" somatic get-pileup-summaries --input "$R/results/preprocessing/recalibrated/$2/$2.recal.bam" \
      --variant "$REF/Annotation/GATKBundle/af-only-gnomad.hg38.vcf.gz" \
      --reference "$REF/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta" \
      --intervals "$bed" --output "$t" 2> "$O/$2.$g.log" || echo "  get-pileup-summaries FAILED: $(tail -n 1 "$O/$2.$g.log")"
    echo "$2 $g: $(( $(date +%s) - s )) s; $(tail -n 1 "$O/$2.$g.log")"
    chk "$t" "$golden" "$2 $g per-group table"
    if diff <(grep "filtered by" "$O/$2.$g.log" | sed 's/ *$//') \
            <(grep "filtered by" "$(dirname "$golden")/.command.log" | sed 's/^.*INFO  GetPileupSummaries - //; s/ *$//') > /dev/null
    then pass=$((pass+1)); echo "  identical: $2 $g read-filter counts and GATK's log"
    else fail=$((fail+1)); echo "  DIFFERS: $2 $g read-filter counts and GATK's log"; fi
    tables+=(-I "$t")
  done
  # sarek gathers in the order its tasks were numbered: `.command.sh` of the gather task has it.
  mapfile -t order < <(grep -o -- '--I [^ ]*' "$G/GATHERPILEUPSUMMARIES_$1-$2/.command.sh" | sed 's/^--I //')
  if (( ${#order[@]} )); then tables=(); for f in "${order[@]}"; do tables+=(-I "$O/$f"); done; fi
  "$bin" somatic gather-pileup-summaries "${tables[@]}" \
    --sequence-dictionary "$REF/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.dict" -O "$O/$2.mutect2.pileups.table"
  chk "$O/$2.mutect2.pileups.table" "$G/GATHERPILEUPSUMMARIES_$1-$2/$2.mutect2.pileups.table" "$2 gathered table"
done

"$bin" somatic calculate-contamination -I "$O/$tumour.mutect2.pileups.table" --matched-normal "$O/$normal.mutect2.pileups.table" \
  -O "$O/$P.contamination.table" --tumor-segmentation "$O/$P.segmentation.table"
for f in contamination segmentation; do
  chk "$O/$P.$f.table" "$G/CALCULATECONTAMINATION-$tumour"_vs_"$normal/$P.$f.table" "$f table"
done

# The orientation model depends on input order; take it from the reference task's command.
L=$G/LEARNREADORIENTATIONMODEL-${tumour}_vs_${normal}
inputs=()
while IFS= read -r f; do inputs+=(--input "$(command ls "$G"/MUTECT2_PAIRED-*/"$f")"); done \
  < <(grep -o -- '--input [^ ]*f1r2\.tar\.gz' "$L/.command.sh" | sed 's/^--input //')
"$bin" somatic learn-read-orientation-model "${inputs[@]}" --output "$O/$P.artifactprior.tar.gz" 2> /dev/null
mkdir -p "$O/ours" "$O/golden"
tar xzf "$O/$P.artifactprior.tar.gz" -C "$O/ours"; tar xzf "$L/$P.artifactprior.tar.gz" -C "$O/golden"
diff <(tar tzf "$O/$P.artifactprior.tar.gz") <(tar tzf "$L/$P.artifactprior.tar.gz") > /dev/null \
  && { pass=$((pass+1)); echo "  identical: orientation archive member names, in order"; } \
  || { fail=$((fail+1)); echo "  DIFFERS: orientation archive member names"; }
for m in "$O"/golden/*; do chk "$O/ours/$(basename "$m")" "$m" "orientation member $(basename "$m")"; done

echo "somatic boundary 1 gate: $pass identical, $fail differ (outputs in $O)"
(( fail == 0 && pass > 0 ))
