#!/usr/bin/env bash
# Gate `nomics qc bam-stats` and `nomics qc mosdepth` against every sarek run under the
# given directories (default ~/omics/runs). For each published `<sample>.{md,recal}.bam`
# that has reports beside it:
#   samtools  reports/samtools/<sample>/<name>.bam.stats            bytes
#   mosdepth  reports/mosdepth/<sample>/<name>.mosdepth.*.txt       bytes
#             <name>.{per-base,regions}.bed.gz                      decompressed bytes
#             <name>.*.bed.gz.csi                                   two-sided (docs/GATES.md)
# mosdepth's form follows the golden: a per-base file means the exome form (`--by <BED>`),
# none means the genome form (`--by 500 --no-per-base --fast-mode`). The BED is `intervals:`
# in the run's params.yaml, looked up as given and then by base name under ~/omics/data;
# BED=path overrides it. Runs published as CRAM are skipped (CRAM input is not ported yet).
#
#   scripts/qc_bam_gate.sh [RUNS_DIR ...]      NOMICS=path/to/nomics to pick the binary
set -uo pipefail
REPO=$(cd "$(dirname "$0")/.." && pwd)
bin=${NOMICS:-$REPO/target/release/nomics}
roots=("$@"); (( ${#roots[@]} )) || roots=("$HOME/omics/runs")
tmp=$(mktemp -d "${TMPDIR:-/tmp}/qc-bam-gate.XXXXXX"); trap 'rm -rf "$tmp"' EXIT
pass=0 fail=0 skip=0
same() { cmp -s "$1" "$2"; }
same_gz() { cmp -s <(gzip -dc "$1") <(gzip -dc "$2"); }

while IFS= read -r bam; do
  dir=$(dirname "$bam"); sample=$(basename "$dir"); name=$(basename "$bam" .bam)   # <sample>.md | <sample>.recal
  results=${dir%/preprocessing/*}; run=${results%/results}
  stats=$results/reports/samtools/$sample/$name.bam.stats
  md=$results/reports/mosdepth/$sample
  problems=() checked=0
  n=$((pass + fail + skip)); out=$tmp/$n; mkdir -p "$out"; [[ -d $tmp/$((n - 1)) ]] && rm -r "$tmp/$((n - 1))"

  if [[ -f $stats ]]; then
    ref=$(sed -n 's/^# The command line was: .*--reference \([^ ]*\) .*/\1/p' "$stats" | head -n 1)
    if "$bin" qc bam-stats --bam "$bam" --reference "${ref:-reference.fasta}" --out "$out/stats" 2> "$out/err"; then
      same "$out/stats" "$stats" || problems+=("samtools stats differs")
    else problems+=("bam-stats failed: $(head -n 1 "$out/err")"); fi
    checked=$((checked + 1))
  fi

  if [[ -f $md/$name.mosdepth.summary.txt ]]; then
    beds=(regions); args=(--by 500 --no-per-base --fast-mode); ready=1
    if [[ -f $md/$name.per-base.bed.gz ]]; then
      bed=${BED:-$(sed -n 's/^intervals: *//p' "$run/params.yaml" 2>/dev/null | tr -d "\"'" | head -n 1)}
      [[ -f $bed ]] || bed=$HOME/omics/data/$(basename "${bed:-none}")
      beds=(per-base regions); args=(--by "$bed")
      [[ -f $bed ]] || { ready=0; problems+=("no target BED for the exome form (set BED=)"); }
    fi
    if (( ready )) && "$bin" qc mosdepth --bam "$bam" --prefix "$out/$name" "${args[@]}" 2> "$out/err"; then
      for f in mosdepth.global.dist.txt mosdepth.region.dist.txt mosdepth.summary.txt; do
        same "$out/$name.$f" "$md/$name.$f" || problems+=("$f differs")
      done
      for b in "${beds[@]}"; do
        same_gz "$out/$name.$b.bed.gz" "$md/$name.$b.bed.gz" || problems+=("$b.bed.gz differs (decompressed)")
        "$bin" qc csi --bed-gz "$md/$name.$b.bed.gz" --out "$out/golden.csi" 2>> "$out/err" \
          && same_gz "$out/golden.csi" "$md/$name.$b.bed.gz.csi" || problems+=("$b csi: builder on the golden file does not give the golden index")
        "$bin" qc csi --bed-gz "$out/$name.$b.bed.gz" --out "$out/ours.csi" 2>> "$out/err" \
          && same "$out/ours.csi" "$out/$name.$b.bed.gz.csi" || problems+=("$b csi: ours is not the builder on our file")
      done
    elif (( ready )); then problems+=("mosdepth failed: $(head -n 1 "$out/err")"); fi
    checked=$((checked + 1))
  fi

  if (( checked == 0 )); then skip=$((skip + 1)); continue; fi
  if (( ${#problems[@]} )); then
    fail=$((fail + 1)); echo "FAIL $run  $name"; printf '     %s\n' "${problems[@]}"
  else
    pass=$((pass + 1)); echo "PASS $run  $name  ($checked of 2 tools)"
  fi
done < <(find -L "${roots[@]}" \( -path '*/results/preprocessing/markduplicates/*' -o -path '*/results/preprocessing/recalibrated/*' \) \
           \( -name '*.md.bam' -o -name '*.recal.bam' \) 2>/dev/null | sort)
echo "qc bam gate: $pass pass, $fail fail, $skip BAMs without reports"
(( fail == 0 && pass > 0 ))
