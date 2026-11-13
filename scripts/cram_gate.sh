#!/usr/bin/env bash
# CRAM gate (docs/GATES.md): our CRAM must decode to what the reference file
# decodes to, under the same decoder: htslib (samtools 1.24 on PATH), htsjdk
# 5.0.0 (Picard 3.5.0, local) and htsjdk 4.2.0 (GATK 4.6.2.0 in sarek's
# container). Headers must match after dropping @SQ UR (a path) and @PG
# lines (provenance: ours records nomics, sarek's the tools it ran; the
# md.bam gate drops @PG the same way). Our .crai
# must equal `samtools index`'s, and region queries must match the
# reference's (htslib and htsjdk).
#
#   scripts/cram_gate.sh <in.bam> <ref.fa> <reference> <profile> [workdir]
#
# <profile> is a `nomics dev bam-to-cram --profile`: sarek-md / sarek-recal
# with <reference> = sarek's md.cram / recal.cram (scripts/cram_reference.sh
# builds them), or lossless with <reference> = <in.bam> (then @SQ M5 is
# dropped too, since a BAM has none). The reference needs an index.
# Every decoder's exit status is checked, and a decode with no records fails
# the gate. Env: NOMICS (nomics binary), THREADS, JAVA/PICARD_JAR, GATK_IMAGE,
# BASELINE=0 to skip the samtools view -C size/time baseline, CRAM_GROUPS=<n>
# to encode <in.bam> as the merge of n interval groups (--groups), OURS=<cram>
# to gate an existing CRAM (e.g. from `nomics run --save-output-as-bam
# false`, with its .crai) instead of encoding <in.bam>. STREAM=1 compares
# decodes by record count and SHA-256 instead of keeping the decoded SAM
# (whole genomes).
set -euo pipefail
BAM=${1:?usage: cram_gate.sh <in.bam> <ref.fa> <reference> <profile> [workdir]}
REF=${2:?}
GOLD=${3:?}
PROFILE=${4:?}
REPO=$(cd "$(dirname "$0")/.." && pwd)
W=${5:-$(mktemp -d "${TMPDIR:-/tmp}/cram-gate.XXXXXX")}
mkdir -p "$W"
W=$(cd "$W" && pwd)
abs() { echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"; }
realp() { python3 -c 'import os, sys; print(os.path.realpath(sys.argv[1]))' "$1"; }
BAM=$(abs "$BAM"); REF=$(abs "$REF"); GOLD=$(abs "$GOLD")
THREADS=${THREADS:-$(sysctl -n hw.ncpu 2>/dev/null || nproc)}
JAVA=${JAVA:-$HOME/miniconda3/envs/rnaseq/lib/jvm/bin/java}
PICARD_JAR=${PICARD_JAR:-$HOME/miniconda3/envs/rnaseq/share/picard-3.5.0-0/picard.jar}
GATK_IMAGE=${GATK_IMAGE:-community.wave.seqera.io/library/gatk4_htslib_samtools:6fd46722be032bdd}
if [[ -z ${NOMICS:-} ]]; then
  cargo build --manifest-path "$REPO/Cargo.toml" --release -q -p nomics
  NOMICS=$REPO/target/release/nomics
fi
name=$(basename "$BAM" .bam)
CRAM=${OURS:+$(abs "$OURS")}
CRAM=${CRAM:-$W/$name.cram}
fail=0
pass() { echo "PASS  $*"; }
bad() { echo "FAIL  $*"; fail=1; }
now() { python3 -c 'import time; print(time.time())'; }
secs() { python3 -c "print(f'{$2 - $1:.2f}')"; }
strip_m5=$([[ $PROFILE == lossless ]] && echo 1 || echo 0)
# Header normalization: @SQ UR (path) and @PG always; M5 only against a BAM.
norm() {
  python3 -c '
import sys
m5 = sys.argv[1] == "1"
for l in sys.stdin:
    if l.startswith("@PG"):
        continue
    if l.startswith("@SQ"):
        f = [x for x in l.rstrip("\n").split("\t") if not x.startswith("UR:") and not (m5 and x.startswith("M5:"))]
        l = "\t".join(f) + "\n"
    sys.stdout.write(l)' "$strip_m5"
}

echo "== $name -> $CRAM (profile $PROFILE, reference $GOLD, $THREADS threads)"
command -v samtools >/dev/null || { echo "samtools not on PATH" >&2; exit 2; }
# Run "$@" with stdout to <file>; a non-zero exit fails the gate (and returns 1).
to() { local out=$1; shift; if "$@" >"$out"; then return 0; fi; bad "command failed: $* (see $out)"; return 1; }
recs_of() { awk '!/^@/' "$1" >"$2"; }
hdr_of() { awk '/^@/' "$1" | norm >"$2"; }
nrec() { awk 'END{print NR}' "$1"; }
if [[ -z ${OURS:-} ]]; then
  "$NOMICS" dev bam-to-cram "$BAM" "$REF" "$CRAM" --profile "$PROFILE" ${CRAM_GROUPS:+--groups "$CRAM_GROUPS"} --threads "$THREADS" 2>"$W/encode.log" \
    || { cat "$W/encode.log"; bad "encode failed"; exit 1; }
  cat "$W/encode.log"
  "$NOMICS" dev bam-to-cram "$BAM" "$REF" "$W/$name.t1.cram" --profile "$PROFILE" ${CRAM_GROUPS:+--groups "$CRAM_GROUPS"} --threads 1 2>"$W/encode1.log" \
    || { cat "$W/encode1.log"; bad "1-thread encode failed"; exit 1; }
  cat "$W/encode1.log"
  # The file id and header (whose @PG records the command line) differ;
  # compare everything from the first data container (the .crai's offset).
  # awk reads the whole stream (no SIGPIPE under pipefail); validated below.
  first() { gzip -dc "$1.crai" | awk -F'\t' 'NR == 1 { v = $4 } END { print v }'; }
  o1=$(first "$CRAM"); o2=$(first "$W/$name.t1.cram")
  if [[ $o1 =~ ^[0-9]+$ && $o2 =~ ^[0-9]+$ ]]; then
    tail -c +$((o1 + 1)) "$CRAM" >"$W/data.n"; tail -c +$((o2 + 1)) "$W/$name.t1.cram" >"$W/data.1"
    cmp -s "$W/data.n" "$W/data.1" && pass "1-thread and $THREADS-thread CRAMs: data containers byte-identical" \
      || bad "CRAM bytes depend on the thread count"
  else bad "cannot read the first container offset from the .crai ($o1 / $o2)"; fi
  rm -f "$W/data.n" "$W/data.1"
  if [[ $PROFILE != lossless ]]; then
    # Our own lossless encode must decode to the BAM (records via htslib).
    "$NOMICS" dev bam-to-cram "$BAM" "$REF" "$W/$name.lossless.cram" --profile lossless --threads "$THREADS" 2>/dev/null \
      || bad "lossless encode failed"
    if to "$W/lossless.sam" samtools view --no-PG -T "$REF" "$W/$name.lossless.cram" \
       && to "$W/bam.sam" samtools view --no-PG "$BAM"; then
      [[ $(nrec "$W/bam.sam") -gt 0 ]] && cmp -s "$W/lossless.sam" "$W/bam.sam" \
        && pass "lossless profile: records identical to the BAM's (htslib)" || bad "lossless profile differs from the BAM"
    fi
  fi
fi
[[ -e $GOLD.crai || -e $GOLD.bai ]] || samtools index "$GOLD"

# Decoders: htslib always; htsjdk 5 (Picard) and htsjdk 4 (GATK) when present.
decs=(st)
[[ -x $JAVA && -f $PICARD_JAR ]] && decs+=(hj5) || echo "SKIP  htsjdk 5 (no Picard at $PICARD_JAR)"
# Docker sees host paths only through mounts: mount real paths (symlinks
# resolved) and check up front that the container can read every input.
RREF=$(realp "$REF")
dk() { docker run --rm --log-driver none -u "$(id -u):$(id -g)" "$@"; }
if docker image inspect "$GATK_IMAGE" >/dev/null 2>&1; then
  seen=1
  for f in "$RREF" "$RREF.fai" "${RREF%.*}.dict" "$(realp "$CRAM")" "$(realp "$GOLD")"; do
    d=$(dirname "$f")
    if ! dk -v "$d:$d" "$GATK_IMAGE" test -r "$f" >/dev/null 2>&1; then
      bad "Docker cannot read $f (is $d shared with Docker?)"; seen=0
    fi
  done
  [[ $seen == 1 ]] && decs+=(hj4)
else
  echo "SKIP  htsjdk 4 (GATK image $GATK_IMAGE not present)"
fi
run_dec() { # <dec> <file> <out>
  case $1 in
    st) samtools view --no-PG -h -T "$REF" "$2" ;;
    hj5) "$JAVA" -Xmx4g -jar "$PICARD_JAR" ViewSam I="$2" R="$REF" QUIET=true VALIDATION_STRINGENCY=SILENT 2>>"$W/picard.log" ;;
    hj4) local f d rd; f=$(realp "$2"); d=$(dirname "$f"); rd=$(dirname "$RREF")
      dk -v "$d:$d" -v "$rd:$rd" "$GATK_IMAGE" \
        gatk --java-options -Xmx4g ViewSam -I "$f" -R "$RREF" --QUIET true --VALIDATION_STRINGENCY SILENT \
        --ALIGNMENT_STATUS All --PF_STATUS All 2>>"$W/gatk.log" ;;
  esac >"$3"
}
# STREAM=1 (whole genomes): keep only each decode's header, record count and
# the SHA-256 of its record lines, never the decoded SAM (hundreds of GB at 30x).
# Every stage of the pipe must succeed.
dec_stream() { # <dec> <file> <out prefix>
  run_dec "$1" "$2" /dev/stdout | awk -v H="$3.sam" -v C="$3.n" '/^@/ { print > H; next } { print; n++ } END { print n + 0 > C }' \
    | shasum -a 256 | awk '{ print $1 }' >"$3.sha"
  local st=("${PIPESTATUS[@]}")
  [[ ${st[0]} == 0 && ${st[1]} == 0 && ${st[2]} == 0 ]]
}
# All decodes in parallel; every one's exit status is checked.
pids=(); jobs_=()
for dec in "${decs[@]}"; do
  for who in ours gold; do
    f=$([[ $who == ours ]] && echo "$CRAM" || echo "$GOLD")
    if [[ ${STREAM:-0} == 1 ]]; then dec_stream "$dec" "$f" "$W/$who.$dec" &
    else run_dec "$dec" "$f" "$W/$who.$dec.sam" & fi
    pids+=($!); jobs_+=("$who.$dec")
  done
done
rm -f "$W"/ok.*
for i in "${!pids[@]}"; do
  if wait "${pids[$i]}"; then touch "$W/ok.${jobs_[$i]}"; else bad "decoder ${jobs_[$i]} failed"; fi
done
for dec in "${decs[@]}"; do
  if [[ $dec == st ]]; then label="htslib $(samtools --version | awk 'NR == 1 { print $2 }')"
  elif [[ $dec == hj5 ]]; then
    label="htsjdk $(unzip -p "$PICARD_JAR" META-INF/MANIFEST.MF | awk -F': ' '/htsjdk-Version/{print $2}' | tr -d '\r') (Picard)"
  else label="htsjdk 4.2.0 (GATK 4.6.2.0)"; fi
  [[ -e $W/ok.ours.$dec && -e $W/ok.gold.$dec ]] || continue
  if [[ ${STREAM:-0} == 1 ]]; then
    n_ours=$(cat "$W/ours.$dec.n"); n_gold=$(cat "$W/gold.$dec.n")
    if [[ $n_ours -eq 0 || $n_gold -eq 0 ]]; then bad "$label: decoded $n_ours / $n_gold records"
    elif [[ $n_ours == "$n_gold" ]] && cmp -s "$W/ours.$dec.sha" "$W/gold.$dec.sha"; then
      pass "$label: all $n_ours records identical to the reference's (SHA-256 $(cut -c1-16 "$W/ours.$dec.sha"))"
    else bad "$label: records differ ($n_ours vs $n_gold, $(cut -c1-16 "$W/ours.$dec.sha") vs $(cut -c1-16 "$W/gold.$dec.sha"))"; fi
    for who in ours gold; do norm <"$W/$who.$dec.sam" >"$W/$who.$dec.hdr"; done
    if [[ ! -s $W/gold.$dec.hdr ]]; then bad "$label: empty header"
    elif cmp -s "$W/ours.$dec.hdr" "$W/gold.$dec.hdr"; then pass "$label: header identical (@SQ UR and @PG removed)"
    else bad "$label: header differs"; diff "$W/ours.$dec.hdr" "$W/gold.$dec.hdr" | awk 'NR <= 4' || true; fi
    continue
  fi
  for who in ours gold; do
    recs_of "$W/$who.$dec.sam" "$W/$who.$dec.recs"; hdr_of "$W/$who.$dec.sam" "$W/$who.$dec.hdr"
  done
  n_ours=$(nrec "$W/ours.$dec.recs"); n_gold=$(nrec "$W/gold.$dec.recs")
  if [[ $n_ours -eq 0 || $n_gold -eq 0 ]]; then bad "$label: decoded $n_ours / $n_gold records"; continue; fi
  if cmp -s "$W/ours.$dec.recs" "$W/gold.$dec.recs"; then pass "$label: all $n_ours records identical to the reference's"
  else bad "$label: records differ ($n_ours vs $n_gold)"; cmp "$W/ours.$dec.recs" "$W/gold.$dec.recs" | awk 'NR <= 2' || true; fi
  if [[ ! -s $W/gold.$dec.hdr ]]; then bad "$label: empty header"
  elif cmp -s "$W/ours.$dec.hdr" "$W/gold.$dec.hdr"; then
    pass "$label: header identical (@SQ UR and @PG removed$([[ $strip_m5 == 1 ]] && echo ', M5 removed'))"
  else bad "$label: header differs"; diff "$W/ours.$dec.hdr" "$W/gold.$dec.hdr" | awk 'NR <= 4' || true; fi
  rm -f "$W/ours.$dec.recs" "$W/gold.$dec.recs"
done

# @SQ M5 values equal samtools dict's.
sq_m5() { awk -F'\t' '$1=="@SQ"{s="";m="";for(i=2;i<=NF;i++){if($i~/^SN:/)s=substr($i,4); if($i~/^M5:/)m=substr($i,4)}; print s"\t"m}' "$1" | sort; }
if to "$W/dict.sam" samtools dict "$REF" && to "$W/ours.hdr.sam" samtools view -H "$CRAM"; then
  sq_m5 "$W/dict.sam" >"$W/dict.m5"; sq_m5 "$W/ours.hdr.sam" >"$W/ours.m5"
  join -t $'\t' "$W/ours.m5" "$W/dict.m5" | awk -F'\t' '$2==$3{print $1"\t"$2}' >"$W/m5.agree"
  if [[ -s $W/ours.m5 ]] && cmp -s "$W/ours.m5" "$W/m5.agree"; then
    pass "@SQ M5 = samtools dict for all $(nrec "$W/ours.m5") contigs"
  else bad "@SQ M5 mismatch vs samtools dict"; fi
fi

# Index: our .crai equals samtools index's; region queries match the reference.
cp "$CRAM" "$W/reindex.cram"
if samtools index "$W/reindex.cram"; then
  gzip -dc "$CRAM.crai" >"$W/ours.crai.txt" && gzip -dc "$W/reindex.cram.crai" >"$W/st.crai.txt"
  if [[ -s $W/ours.crai.txt ]] && cmp -s "$W/ours.crai.txt" "$W/st.crai.txt"; then
    pass ".crai identical to samtools index ($(nrec "$W/ours.crai.txt") slices)"
  else bad ".crai differs from samtools index"; fi
else bad "samtools index of our CRAM failed"; fi
awk -F'\t' '$1=="@SQ"{print substr($2,4)}' "$W/ours.hdr.sam" >"$W/names"
awk -v names="$W/names" '
  BEGIN { while ((getline n < names) > 0) name[i++] = n }
  $1 >= 0 { k++; r[k] = name[$1]; s[k] = $2; e[k] = $2 + $3 - 1 }
  END {
    step = int(k / 6) + 1
    for (j = 1; j <= k; j += step) {
      print r[j] ":" s[j] "-" int((s[j] + e[j]) / 2)
      if (j < k && r[j + 1] == r[j]) print r[j] ":" e[j] - 500 "-" s[j + 1] + 500
    }
    print r[1]
  }' "$W/ours.crai.txt" >"$W/regions"
echo '*' >>"$W/regions"
nreg=0; nrecs=0
# STREAM=1: a region's records as "<count> <SHA-256>" instead of the SAM.
region() { # <out> <cram> <region>
  # </dev/null: a containerized samtools (docker -i) would read the region list.
  if [[ ${STREAM:-0} != 1 ]]; then to "$1" samtools view --no-PG -T "$REF" "$2" "$3" </dev/null; return; fi
  samtools view --no-PG -T "$REF" "$2" "$3" </dev/null | awk '{ print; n++ } END { print n + 0 > "/dev/stderr" }' 2>"$1.n" \
    | shasum -a 256 | awk '{ print $1 }' >"$1.sha"
  local st=("${PIPESTATUS[@]}")
  [[ ${st[0]} == 0 && ${st[1]} == 0 && ${st[2]} == 0 ]] || { bad "command failed: samtools view $2 $3"; return 1; }
  echo "$(cat "$1.n") $(cat "$1.sha")" >"$1"
}
nrec_of() { if [[ ${STREAM:-0} == 1 ]]; then awk '{ print $1 }' "$1"; else nrec "$1"; fi; }
# Each region's count must also equal the source BAM's (through its .bai):
# an index bug that drops slices shows as a short count on both sides only
# if the reference shares it, never against the BAM.
bam_idx=0
[[ -e $BAM.bai || -e ${BAM%.bam}.bai || -e $BAM.csi ]] && bam_idx=1
[[ $bam_idx == 1 ]] || echo "SKIP  region counts vs the BAM (no index for $BAM)"
ncount=0
while read -r reg; do
  if region "$W/r.ours" "$CRAM" "$reg" && region "$W/r.gold" "$GOLD" "$reg"; then
    if cmp -s "$W/r.ours" "$W/r.gold"; then nreg=$((nreg + 1)); nrecs=$((nrecs + $(nrec_of "$W/r.ours")))
    else bad "htslib region $reg differs"; fi
    if [[ $bam_idx == 1 ]]; then
      if to "$W/r.bam.n" samtools view -c "$BAM" "$reg" </dev/null; then
        n_bam=$(awk 'NR == 1 { print $1 }' "$W/r.bam.n"); n_ours=$(nrec_of "$W/r.ours")
        if [[ $n_bam == "$n_ours" ]]; then ncount=$((ncount + 1))
        else bad "region $reg: $n_ours records from our CRAM, $n_bam from the BAM"; fi
      fi
    fi
  fi
done <"$W/regions"
[[ $nrecs -gt 0 ]] || bad "htslib region queries returned no records"
nq=$(awk 'END { print NR }' "$W/regions")
[[ $nreg == "$nq" ]] && pass "htslib: $nreg region queries identical to the reference's ($nrecs records; $(paste -sd' ' "$W/regions"))"
[[ $bam_idx == 1 ]] && pass "htslib: $ncount region counts equal the source BAM's"
if [[ " ${decs[*]} " == *" hj5 "* ]]; then
  grep '^@SQ' "$W/ours.hdr.sam" | python3 -c 'import sys; [sys.stdout.write("\t".join(f for f in l.rstrip("\n").split("\t") if not f.startswith("UR:"))+"\n") for l in sys.stdin]' >"$W/regions.interval_list"
  awk 'index($0, ":") { split($0, a, /[:-]/); print a[1] "\t" a[2] "\t" a[3] "\t+\tr" NR }' "$W/regions" \
    >>"$W/regions.interval_list"
  iv_ok=1
  for f in ours:"$CRAM" gold:"$GOLD"; do
    to "$W/${f%%:*}.regions.full" "$JAVA" -Xmx4g -jar "$PICARD_JAR" ViewSam I="${f#*:}" R="$REF" \
      INTERVAL_LIST="$W/regions.interval_list" QUIET=true VALIDATION_STRINGENCY=SILENT 2>>"$W/picard.log" || iv_ok=0
    recs_of "$W/${f%%:*}.regions.full" "$W/${f%%:*}.regions.sam"
  done
  if [[ $iv_ok == 1 ]]; then
    if [[ $(nrec "$W/gold.regions.sam") -gt 0 ]] && cmp -s "$W/ours.regions.sam" "$W/gold.regions.sam"; then
      pass "htsjdk: interval query identical to the reference's ($(nrec "$W/ours.regions.sam") records)"
    else bad "htsjdk interval query differs or is empty"; fi
  fi
fi

echo "size: ours $(wc -c <"$CRAM" | tr -d ' ') bytes, reference $(wc -c <"$GOLD" | tr -d ' ') bytes, BAM $(wc -c <"$BAM" | tr -d ' ') bytes"
if [[ ${BASELINE:-1} != 0 ]]; then
  for t in 1 "$THREADS"; do
    a=$(now); samtools view -@ "$t" -C -T "$REF" --output-fmt-option version=3.0 -o "$W/st30.$t.cram" "$BAM"; b=$(now)
    echo "samtools 1.24 view -C version=3.0 -@ $t: $(secs "$a" "$b") s, $(wc -c <"$W/st30.$t.cram" | tr -d ' ') bytes"
  done
  a=$(now); samtools view -@ "$THREADS" -b -o "$W/st.bam" "$BAM"; b=$(now)
  echo "samtools view -b -@ $THREADS (BAM rewrite): $(secs "$a" "$b") s"
fi
echo "== $([[ $fail == 0 ]] && echo PASS || echo FAIL) ($W)"
exit $fail
