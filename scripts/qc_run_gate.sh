#!/usr/bin/env bash
# Gate the reports/ tree of a `nomics run --reports` output against a sarek run's.
# Every file under the golden's results/reports is accounted for: compared under its rule in
# docs/GATES.md, or listed as out of scope (FastQC's HTML; of its zip, only the two text
# members are compared), or reported MISSING. A file of ours with no golden counterpart is
# reported EXTRA. Nothing is read from the golden to make our side: the run wrote it.
#
#   fastp/**.fastp.{json,html,log}         bytes but for the exempt fields (scripts/fastp_norm.py)
#   fastqc/<p>/<p>_N_fastqc.zip            members fastqc_data.txt, summary.txt against our
#                                          fastqc/<p>/<p>_N_fastqc/ directory, bytes
#   markduplicates/**.metrics              bytes except lines 2 and 4 (command line, start time)
#   samtools/**.stats, mosdepth/**.txt     bytes
#   mosdepth/**.bed.gz                     decompressed bytes
#   mosdepth/**.csi                        two-sided: our builder on the golden file gives the
#                                          golden index, and ours is the builder on our file
#   vcftools/**, bcftools/**               bytes; TsTv.qual under its exempt-field rule
#
#   scripts/qc_run_gate.sh GOLDEN_RUN OUR_OUTDIR     NOMICS=path/to/nomics, ARCH=aarch64|x86_64
set -uo pipefail
REPO=$(cd "$(dirname "$0")/.." && pwd)
bin=${NOMICS:-$REPO/target/release/nomics}
G=${1:?usage: qc_run_gate.sh GOLDEN_RUN OUR_OUTDIR}/results/reports; O=${2:?usage: qc_run_gate.sh GOLDEN_RUN OUR_OUTDIR}/reports
[[ -d $G ]] || { echo "no $G"; exit 2; }
[[ -d $O ]] || { echo "no $O"; exit 2; }
arch=${ARCH:-aarch64}; [[ $1 == *x86* ]] && arch=${ARCH:-x86_64}
tmp=$(mktemp -d "${TMPDIR:-/tmp}/qc-run-gate.XXXXXX"); trap 'rm -rf "$tmp"' EXIT
pass=0 fail=0 missing=0 scope=0 extra=0
ok() { pass=$((pass + 1)); echo "PASS     $1${2:+  ($2)}"; }
bad() { fail=$((fail + 1)); echo "FAIL     $1  ($2)"; }
norm() { python3 -I "$REPO/scripts/fastp_norm.py" "$1" "$2"; }
same_gz() { cmp -s <(gzip -dc "$1") <(gzip -dc "$2"); }

while IFS= read -r g; do
  rel=${g#"$G"/}; o=$O/$rel
  case $rel in
    fastqc/*_fastqc.html) scope=$((scope + 1)); echo "SCOPE    $rel  (FastQC's rendered page is not reproduced)"; continue ;;
    fastqc/*_fastqc.zip)
      d=${o%.zip}
      if [[ ! -d $d ]]; then missing=$((missing + 1)); echo "MISSING  $rel  (no ${d#"$O"/}/)"; continue; fi
      problems=()
      for m in fastqc_data.txt summary.txt; do
        unzip -p "$g" "*/$m" > "$tmp/m" 2> /dev/null && cmp -s "$tmp/m" "$d/$m" || problems+=("$m")
      done
      (( ${#problems[@]} )) && bad "$rel" "differs: ${problems[*]}" || ok "$rel" "two text members; images and page not reproduced"
      continue ;;
  esac
  if [[ ! -f $o ]]; then missing=$((missing + 1)); echo "MISSING  $rel"; continue; fi
  case $rel in
    fastp/*.fastp.json|fastp/*.fastp.html|fastp/*.fastp.log)
      e=${rel##*.}
      if cmp -s "$o" "$g"; then ok "$rel"
      elif cmp -s <(norm "$e" "$o") <(norm "$e" "$g"); then ok "$rel" "outside the exempt fields"
      else bad "$rel" "differs outside the exempt fields"; fi ;;
    markduplicates/*.metrics)
      cmp -s <(sed '2d;4d' "$o") <(sed '2d;4d' "$g") && ok "$rel" "but lines 2 and 4" || bad "$rel" "differs outside lines 2 and 4" ;;
    mosdepth/*.bed.gz)
      same_gz "$o" "$g" && ok "$rel" "decompressed" || bad "$rel" "decompressed bytes differ" ;;
    mosdepth/*.csi)
      if "$bin" qc csi --bed-gz "${g%.csi}" --out "$tmp/g.csi" 2> "$tmp/err" && same_gz "$tmp/g.csi" "$g" \
         && "$bin" qc csi --bed-gz "${o%.csi}" --out "$tmp/o.csi" 2>> "$tmp/err" && cmp -s "$tmp/o.csi" "$o"
      then ok "$rel" "two-sided"; else bad "$rel" "two-sided index rule: $(head -n 1 "$tmp/err")"; fi ;;
    vcftools/*.TsTv.qual)
      if cmp -s "$o" "$g"; then ok "$rel"
      else
        # The exempt fields are decided from our VCF: the reports made from it must be the ones
        # the run wrote, and must pass against the golden under the rule.
        vt=$(dirname "$g"); sub=${vt#"$G"/vcftools/}; base=$(basename "$g" .TsTv.qual)
        vcf=$(dirname "$O")/$base.vcf.gz
        out=$("$bin" qc vcf-gate --vcf "$vcf" --vcftools "$vt" --bcftools "$G/bcftools/$sub" --reference-arch "$arch" 2>&1)
        "$bin" qc vcf-stats --vcf "$vcf" --outdir "$tmp/vs" --reference-arch "$arch" 2> /dev/null
        if [[ $out == *"verdict: PASS"* ]] && cmp -s "$tmp/vs/$base.TsTv.qual" "$o"; then ok "$rel" "exempt fields only"
        else bad "$rel" "differs outside the exempt fields"; fi
      fi ;;
    *) cmp -s "$o" "$g" && ok "$rel" || bad "$rel" "bytes differ" ;;
  esac
done < <(find -L "$G" -type f | sort)

while IFS= read -r o; do
  rel=${o#"$O"/}
  case $rel in
    fastqc/*_fastqc/fastqc_data.txt|fastqc/*_fastqc/summary.txt) [[ -f $G/$(dirname "$rel").zip ]] && continue ;;
  esac
  [[ -f $G/$rel ]] || { extra=$((extra + 1)); echo "EXTRA    $rel"; }
done < <(find -L "$O" -type f | sort)

echo "qc run gate: $pass pass, $fail fail, $missing missing, $extra extra, $scope out of scope"
(( fail == 0 && missing == 0 && extra == 0 && pass > 0 ))
