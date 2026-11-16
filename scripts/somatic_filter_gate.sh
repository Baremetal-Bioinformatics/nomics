#!/usr/bin/env bash
# Gate `nomics somatic filter-mutect-calls` against GATK 4.6.2.0 FilterMutectCalls on the cases
# made by scripts/somatic_filter_cases.sh (cases.tsv: name, raw calls, contamination table).
# Per case (docs/GATES.md): the filteringStats.tsv by bytes, every VCF record, and the header
# but for the tool's own ##GATKCommandLine line. No reference tool is run here.
#
#   scripts/somatic_filter_gate.sh [CASES_DIR ...]     NOMICS=path/to/nomics to pick the binary
#   default: ~/omics/runs/somatic/fmc-cases; REF=reference bundle (default ~/omics/ref/GATK.GRCh38)
set -uo pipefail
REPO=$(cd "$(dirname "$0")/.." && pwd)
bin=${NOMICS:-$REPO/target/release/nomics}
REF=${REF:-$HOME/omics/ref/GATK.GRCh38}/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta
dirs=("$@"); (( ${#dirs[@]} )) || dirs=("$HOME/omics/runs/somatic/fmc-cases")
T=$(mktemp -d "${TMPDIR:-/tmp}/somatic-filter-gate.XXXXXX"); trap 'rm -rf "$T"' EXIT
pass=0 fail=0
for D in "${dirs[@]}"; do
  [[ -f $D/cases.tsv ]] || { echo "no $D/cases.tsv"; fail=$((fail + 1)); continue; }
  while IFS=$'\t' read -r name variant contam; do
    g=$D/$name.filtered.vcf.gz; o=$T/$name.vcf.gz
    if ! "$bin" somatic filter-mutect-calls --variant "$D/$variant.vcf.gz" --stats "$D/$variant.vcf.gz.stats" --reference "$REF" \
         --orientation-bias-artifact-priors "$D/priors.tar.gz" --tumor-segmentation "$D/seg.table" \
         --contamination-table "$D/$contam" --output "$o" 2> "$T/log"; then
      fail=$((fail + 1)); echo "FAIL $name  (the port failed: $(tail -n 1 "$T/log"))"; continue
    fi
    n=$(gzip -dc "$g" | grep -vc '^#')
    same=$(paste -d '\n' <(gzip -dc "$g" | grep -v '^#') <(gzip -dc "$o" | grep -v '^#') | awk 'NR%2{a=$0; next} a==$0{s++} END{print s+0}')
    ours=$(gzip -dc "$o" | grep -vc '^#')
    problems=()
    [[ $same == "$n" && $ours == "$n" ]] || problems+=("records $same of $n identical, ours has $ours")
    cmp -s "$g.filteringStats.tsv" "$o.filteringStats.tsv" || problems+=("filteringStats.tsv differs")
    diff <(gzip -dc "$g" | grep '^#' | grep -v '^##GATKCommandLine=<ID=FilterMutectCalls') \
         <(gzip -dc "$o" | grep '^#' | grep -v '^##GATKCommandLine=<ID=FilterMutectCalls') > /dev/null || problems+=("header differs")
    if (( ${#problems[@]} )); then fail=$((fail + 1)); echo "FAIL $name  (${problems[*]})"
    else pass=$((pass + 1)); echo "PASS $name  ($n records, $(gzip -dc "$o" | awk -F'\t' '!/^#/ && $7=="PASS"' | wc -l | tr -d ' ') PASS; stats and header identical)"; fi
  done < "$D/cases.tsv"
done
echo "somatic filter gate: $pass pass, $fail fail"
(( fail == 0 && pass > 0 ))
