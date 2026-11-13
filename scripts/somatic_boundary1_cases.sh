#!/usr/bin/env bash
# Cases beyond the one golden for the two tools of the first Mutect2 boundary that fit a model:
# CalculateContamination on subsets of the gathered pileup tables, and LearnReadOrientationModel
# on each scatter archive alone and in both orders (its output depends on the order).
#
#   scripts/somatic_boundary1_cases.sh make [DIR] [GOLDEN]   run GATK 4.6.2.0 in sarek's container
#                                                            (docker; on the Spark, under the lock)
#   scripts/somatic_boundary1_cases.sh gate [DIR]            the nomics ports against those outputs,
#                                                            bytes; NOMICS=path/to/nomics
#   defaults: DIR=~/omics/runs/somatic/b1-cases, GOLDEN=~/omics/runs/somatic/S1-r1/golden
set -uo pipefail
mode=${1:?usage: somatic_boundary1_cases.sh make|gate [DIR] [GOLDEN]}
O=${2:-$HOME/omics/runs/somatic/b1-cases}
REPO=$(cd "$(dirname "$0")/.." && pwd)

if [[ $mode == make ]]; then
  G=${3:-$HOME/omics/runs/somatic/S1-r1/golden}
  IMG=${GATK_IMAGE:-community.wave.seqera.io/library/gatk4_htslib_samtools:6fd46722be032bdd}
  X=HCC1395_WES_vs_HCC1395BL_WES.mutect2
  mkdir -p "$O"; cd "$O" || exit 2
  cp "$G/GATHERPILEUPSUMMARIES_TUMOR-HCC1395_WES/HCC1395_WES.mutect2.pileups.table" tumour.table
  cp "$G/GATHERPILEUPSUMMARIES_NORMAL-HCC1395BL_WES/HCC1395BL_WES.mutect2.pileups.table" normal.table
  cp "$G/MUTECT2_PAIRED-HCC1395_WES_vs_HCC1395BL_WES/$X.chr17_38877931-38878296.f1r2.tar.gz" a.f1r2.tar.gz
  cp "$G/MUTECT2_PAIRED-HCC1395_WES_vs_HCC1395BL_WES.2/$X.chr1_11981-12351.f1r2.tar.gz" b.f1r2.tar.gz
  sub() {  # sub NAME AWK_CONDITION [SEED]: rows of both tables
    for s in tumour normal; do
      awk -v seed="${3:-1}" 'BEGIN{srand(seed)} /^#/ || $1=="contig" {print; next} '"$2" $s.table > "$1.$s.table"
    done
  }
  sub all '{print}'
  sub odd '{n++; if (n % 2) print}'
  sub chr1to7 '$1 ~ /^chr[1-7]$/'
  sub chr8on '$1 !~ /^chr[1-7]$/'
  sub chr1 '$1 == "chr1"'
  sub rand11 'rand() < 0.7' 11
  sub rand23 'rand() < 0.7' 23
  sub rand47 'rand() < 0.3' 47
  run() { docker run --rm -u "$(id -u):$(id -g)" -v "$O":/w -w /w "$IMG" bash -c "$1"; }
  fail=0
  : > contamination.cases
  for c in all odd chr1to7 chr8on chr1 rand11 rand23 rand47; do
    if run "gatk --java-options '-Xmx8g -XX:-UsePerfData' CalculateContamination --input $c.tumour.table --matched-normal $c.normal.table --output $c.contamination.table --tumor-segmentation $c.segmentation.table --tmp-dir . > $c.contam.log 2>&1"
    then echo "$c" >> contamination.cases; echo "contamination $c: $(sed -n 2p "$c.contamination.table")"
    else echo "contamination $c: GATK failed (see $O/$c.contam.log)"; fail=1; fi
  done
  : > orientation.cases
  lrom() {  # lrom NAME INPUT...
    local n=$1; shift; local args=""; for f in "$@"; do args="$args --input $f"; done
    if run "gatk --java-options '-Xmx8g -XX:-UsePerfData' LearnReadOrientationModel $args --output $n.artifactprior.tar.gz --tmp-dir . > $n.lrom.log 2>&1"
    then printf '%s\t%s\n' "$n" "$*" >> orientation.cases; echo "orientation $n: made"
    else echo "orientation $n: GATK failed (see $O/$n.lrom.log)"; fail=1; fi
  }
  lrom ab a.f1r2.tar.gz b.f1r2.tar.gz
  lrom ba b.f1r2.tar.gz a.f1r2.tar.gz
  lrom a a.f1r2.tar.gz
  lrom b b.f1r2.tar.gz
  # The container must reproduce the sarek golden on the unaltered input.
  if cmp -s all.contamination.table "$G/CALCULATECONTAMINATION-HCC1395_WES_vs_HCC1395BL_WES/$X.contamination.table" \
     && cmp -s all.segmentation.table "$G/CALCULATECONTAMINATION-HCC1395_WES_vs_HCC1395BL_WES/$X.segmentation.table"
  then echo "all: this container reproduces the golden contamination and segmentation tables"
  else echo "all: this container does NOT reproduce the golden"; fail=1; fi
  exit $fail
fi

[[ $mode == gate ]] || { echo "mode must be make or gate"; exit 2; }
bin=${NOMICS:-$REPO/target/release/nomics}
T=$(mktemp -d "${TMPDIR:-/tmp}/somatic-b1-cases.XXXXXX"); trap 'rm -rf "$T"' EXIT
pass=0 fail=0
ok() { pass=$((pass + 1)); echo "PASS $1"; }
bad() { fail=$((fail + 1)); echo "FAIL $1  ($2)"; }
while IFS= read -r c; do
  if ! "$bin" somatic calculate-contamination -I "$O/$c.tumour.table" --matched-normal "$O/$c.normal.table" \
       -O "$T/$c.contamination.table" --tumor-segmentation "$T/$c.segmentation.table" 2> "$T/err"; then
    bad "contamination $c" "the port failed: $(head -n 1 "$T/err")"; continue
  fi
  if cmp -s "$T/$c.contamination.table" "$O/$c.contamination.table" && cmp -s "$T/$c.segmentation.table" "$O/$c.segmentation.table"
  then ok "contamination $c  (estimate $(sed -n 2p "$O/$c.contamination.table" | cut -f2), $(( $(wc -l < "$O/$c.segmentation.table") - 2 )) segments)"
  else bad "contamination $c" "contamination or segmentation table differs"; fi
done < "$O/contamination.cases"
while IFS=$'\t' read -r n inputs; do
  args=(); for f in $inputs; do args+=(--input "$O/$f"); done
  if ! "$bin" somatic learn-read-orientation-model "${args[@]}" --output "$T/$n.tar.gz" 2> "$T/err"; then
    bad "orientation $n" "the port failed: $(head -n 1 "$T/err")"; continue
  fi
  mkdir -p "$T/o.$n" "$T/g.$n"; tar xzf "$T/$n.tar.gz" -C "$T/o.$n"; tar xzf "$O/$n.artifactprior.tar.gz" -C "$T/g.$n"
  if diff <(tar tzf "$T/$n.tar.gz") <(tar tzf "$O/$n.artifactprior.tar.gz") > /dev/null && diff -r "$T/o.$n" "$T/g.$n" > /dev/null
  then ok "orientation $n  (inputs: $inputs)"; else bad "orientation $n" "member names or bytes differ"; fi
done < "$O/orientation.cases"
echo "somatic boundary 1 cases: $pass pass, $fail fail"
(( fail == 0 && pass > 0 ))
