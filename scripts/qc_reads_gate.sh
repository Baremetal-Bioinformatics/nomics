#!/usr/bin/env bash
# Gate `nomics qc fastqc` and the fastp reports of `nomics prep` against every sarek run under
# the given directories (default ~/omics/runs). Inputs come from each run's samplesheet.csv
# (patient,sex,status,sample,lane,fastq_1,fastq_2); a row whose FASTQs are gone is skipped.
#   FastQC  reports/fastqc/<sample>-<lane>/<sample>-<lane>_{1,2}_fastqc.zip
#             members fastqc_data.txt, summary.txt                         bytes
#   fastp   reports/fastp/<sample>/<sample>-<lane>.fastp.{json,html,log}    bytes but for the
#             fields docs/GATES.md names (duplication rate; clock values in the html and log)
# The FASTQ's staged name and fastp's command line and thread count are read from the golden.
# Each fastp file is first compared as it is, then with the named fields blanked, and the
# output says which of the two it took.
#
#   scripts/qc_reads_gate.sh [RUNS_DIR ...]     NOMICS=path/to/nomics to pick the binary
set -uo pipefail
REPO=$(cd "$(dirname "$0")/.." && pwd)
bin=${NOMICS:-$REPO/target/release/nomics}
roots=("$@"); (( ${#roots[@]} )) || roots=("$HOME/omics/runs")
tmp=$(mktemp -d "${TMPDIR:-/tmp}/qc-reads-gate.XXXXXX"); trap 'rm -rf "$tmp"' EXIT
pass=0 fail=0 skip=0 exact=0 exempt=0

# Blank the fields the gate exempts, and nothing else.
norm() { python3 -I "$REPO/scripts/fastp_norm.py" "$1" "$2"; }   # norm json|html|log FILE

while IFS= read -r sheet; do
  run=$(dirname "$sheet"); reports=$run/results/reports
  [[ -d $reports/fastqc || -d $reports/fastp ]] || continue
  while IFS=, read -r _ _ _ sample lane fq1 fq2 _; do
    [[ $sample == sample || -z ${fq2:-} ]] && continue
    p=$sample-$lane; problems=() notes=() checked=0
    n=$((pass + fail + skip)); out=$tmp/$n; mkdir -p "$out"; [[ -d $tmp/$((n - 1)) ]] && rm -r "$tmp/$((n - 1))"
    if [[ ! -f $fq1 || ! -f $fq2 ]]; then skip=$((skip + 1)); continue; fi

    m=0
    for fq in "$fq1" "$fq2"; do
      m=$((m + 1)); zip=$reports/fastqc/$p/${p}_${m}_fastqc.zip
      [[ -f $zip ]] || continue
      mkdir -p "$out/g$m" "$out/o$m"
      unzip -p "$zip" '*/fastqc_data.txt' > "$out/g$m/fastqc_data.txt" && unzip -p "$zip" '*/summary.txt' > "$out/g$m/summary.txt" \
        || { problems+=("fastqc mate $m: cannot read the golden zip"); continue; }
      name=$(sed -n 's/^Filename\t//p' "$out/g$m/fastqc_data.txt" | head -n 1)
      if "$bin" qc fastqc --fastq "$fq" --outdir "$out/o$m" --input-name "$name" 2> "$out/err"; then
        for f in fastqc_data.txt summary.txt; do
          cmp -s "$out/o$m/$f" "$out/g$m/$f" || problems+=("fastqc mate $m: $f differs")
        done
      else problems+=("fastqc mate $m failed: $(head -n 1 "$out/err")"); fi
      checked=$((checked + 1))
    done

    g=$reports/fastp/$sample/$p.fastp
    if [[ -f $g.json ]]; then
      cl=$(python3 -I -c 'import json,sys; sys.stdout.write(json.load(open(sys.argv[1]))["command"])' "$g.json")
      thr=$(sed -n 's/.*--thread \([0-9][0-9]*\).*/\1/p' <<< "$cl"); mkdir -p "$out/fp"
      # Relative report names, as sarek's task has them: the log prints them as given.
      if (cd "$out/fp" && "$bin" prep --reports-only --fastq1 "$fq1" --fastq2 "$fq2" --outdir . --prefix "$p" \
            --fastp-threads "${thr:-12}" --json "$p.fastp.json" --html "$p.fastp.html" --log "$p.fastp.log" \
            --report-command "$cl") 2> "$out/err"; then
        for e in json html log; do
          [[ -f $g.$e ]] || continue
          if cmp -s "$out/fp/$p.fastp.$e" "$g.$e"; then exact=$((exact + 1)); notes+=("$e exact")
          elif cmp -s <(norm "$e" "$out/fp/$p.fastp.$e") <(norm "$e" "$g.$e"); then exempt=$((exempt + 1)); notes+=("$e exempt-fields")
          else problems+=("fastp $e differs outside the exempt fields"); fi
        done
        r1=$(sed -n 's/.*"rate": *\([-+.0-9eE]*\).*/\1/p' "$out/fp/$p.fastp.json" | head -n 1)
        r2=$(sed -n 's/.*"rate": *\([-+.0-9eE]*\).*/\1/p' "$g.json" | head -n 1)
        notes+=("dup ours $r1 golden $r2")
      else problems+=("prep failed: $(head -n 1 "$out/err")"); fi
      checked=$((checked + 1))
    fi

    if (( checked == 0 )); then skip=$((skip + 1)); continue; fi
    if (( ${#problems[@]} )); then
      fail=$((fail + 1)); echo "FAIL $run  $p"; printf '     %s\n' "${problems[@]}"
    else
      pass=$((pass + 1)); echo "PASS $run  $p  ($checked reports; ${notes[*]:-fastqc only})"
    fi
  done < "$sheet"
done < <(find -L "${roots[@]}" -name samplesheet.csv -not -path '*/work/*' 2>/dev/null | sort)
echo "qc reads gate: $pass pass, $fail fail, $skip rows skipped; fastp files: $exact exact, $exempt identical outside the exempt fields"
(( fail == 0 && pass > 0 ))
