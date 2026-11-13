#!/usr/bin/env bash
# Gate `nomics qc vcf-stats` against every sarek run under the given directories
# (default ~/omics/runs): for each published HaplotypeCaller VCF that has vcftools
# and bcftools reports beside it, compute the four reports and compare them with
# sarek's (bytes; .TsTv.qual under the normalization in docs/GATES.md). vcftools
# spells 0/0 differently on x86 and arm64, so the run's architecture is passed.
#
#   scripts/qc_vcf_gate.sh [RUNS_DIR ...]      NOMICS=path/to/nomics to pick the binary
set -uo pipefail
REPO=$(cd "$(dirname "$0")/.." && pwd)
bin=${NOMICS:-$REPO/target/release/nomics}
roots=("$@"); (( ${#roots[@]} )) || roots=("$HOME/omics/runs")
pass=0 fail=0
while IFS= read -r qual; do
  vt=$(dirname "$qual"); results=${vt%/reports/vcftools/*}; rel=${vt#"$results"/reports/vcftools/}
  base=$(basename "$qual" .TsTv.qual)
  vcf=$results/variant_calling/$rel/$base.vcf.gz
  bt=$results/reports/bcftools/$rel
  [[ -f $vcf && -d $bt ]] || { echo "SKIP $qual (no VCF or bcftools report beside it)"; continue; }
  arch=aarch64; [[ $results == *x86* ]] && arch=x86_64   # the goldens made on x86 say so in their names
  out=$("$bin" qc vcf-gate --vcf "$vcf" --vcftools "$vt" --bcftools "$bt" --reference-arch "$arch" 2>&1)
  if [[ $out == *"verdict: PASS"* ]]; then
    pass=$((pass + 1))
    note=$(grep -o '[1-9][0-9]* exempt fields differ.*' <<<"$out" || true)
    echo "PASS ${results%/results}  $base${note:+  ($note)}"
  else
    fail=$((fail + 1)); echo "FAIL ${results%/results}  $base"; sed 's/^/     /' <<<"$out"
  fi
done < <(find -L "${roots[@]}" -path '*/results/reports/vcftools/*' -name '*.TsTv.qual' 2>/dev/null | sort)
echo "qc vcf gate: $pass pass, $fail fail"
(( fail == 0 && pass > 0 ))
