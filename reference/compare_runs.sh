#!/usr/bin/env bash
# Compare two reference runs at every stage boundary.
# Usage: compare_runs.sh <runA_dir> <runB_dir> > report.tsv
# Columns: stage  file  check  result  detail
#
# Checks, from strictest to loosest:
#   bytes    cmp of the published file
#   records  BAM: md5 of `samtools view` body (header excluded; @PG/@CO lines carry paths)
#            VCF: md5 of non-header lines
#   header   md5 of header with @PG/@CO and ##*Command* / ##fileDate lines dropped
#   flags    BAM: count of records whose FLAG differs (same read order assumed)
set -uo pipefail
A=$(cd "$1" && pwd) || exit 2
B=$(cd "$2" && pwd) || exit 2
IMG=community.wave.seqera.io/library/htslib_samtools:1.23.1--aedf2bc6d9d9dddb
st(){ docker run --rm -u "$(id -u):$(id -g)" -v "$HOME/omics:$HOME/omics" -w "$PWD" "$IMG" samtools "$@"; }
out(){ printf '%s\t%s\t%s\t%s\t%s\n' "$@"; }
# A failing tool must abort the comparison, never hash its (empty) output.
fatal(){ echo "FATAL: $*" >&2; exit 3; }
digest(){ local d; d=$(set -o pipefail; "$@" | md5sum | cut -c1-32) || fatal "command failed: $*"; echo "$d"; }
# grep exits 1 when it selects no lines; that is a valid (empty) result, not an error.
g(){ grep "$@"; local rc=$?; [ $rc -le 1 ]; }
bam_header(){ set -o pipefail; st view -H "$1" | g -vE '^@(PG|CO)'; }
vcf_body(){ set -o pipefail; zcat "$1" | g -v '^#'; }
vcf_header(){ set -o pipefail; zcat "$1" | g '^#' | g -vE '^##(GATKCommandLine|fileDate|source|.*[Cc]ommand)'; }

rel_files(){ (cd "$1/results" && find preprocessing variant_calling -type f \( -name '*.bam' -o -name '*.cram' -o -name '*.recal.table' -o -name '*.vcf.gz' -o -name '*.g.vcf.gz' \) | sort); }

for f in $(rel_files "$A"); do
  a=$A/results/$f; b=$B/results/$f
  stage=$(cut -d/ -f2 <<<"$f")
  [ -f "$b" ] || { out "$stage" "$f" exists FAIL "missing in B"; continue; }
  if cmp -s "$a" "$b"; then out "$stage" "$f" bytes PASS ""; continue; fi
  out "$stage" "$f" bytes DIFF "$(stat -c%s "$a") vs $(stat -c%s "$b") bytes"
  case "$f" in
    *.bam|*.cram)
      ha=$(digest bam_header "$a") || exit 3
      hb=$(digest bam_header "$b") || exit 3
      [ "$ha" = "$hb" ] && out "$stage" "$f" header PASS "" || out "$stage" "$f" header DIFF "$ha $hb"
      ra=$(digest st view "$a") || exit 3; rb=$(digest st view "$b") || exit 3
      if [ "$ra" = "$rb" ]; then out "$stage" "$f" records PASS ""; else
        out "$stage" "$f" records DIFF "$ra $rb"
        n=$(paste <(st view "$a" | cut -f1,2,3,4) <(st view "$b" | cut -f1,2,3,4) | awk -F'\t' 'NF!=8||$1!=$5||$2!=$6||$3!=$7||$4!=$8' | wc -l)
        out "$stage" "$f" qname-flag-pos DIFF "$n records differ in QNAME/FLAG/RNAME/POS (paired by row)"
      fi ;;
    *.vcf.gz)
      ba=$(digest vcf_body "$a") || exit 3; bb=$(digest vcf_body "$b") || exit 3
      nrec=$(set -o pipefail; vcf_body "$a" | wc -l) || fatal "vcf_body $a"
      if [ "$ba" = "$bb" ]; then out "$stage" "$f" records PASS "$nrec records"; else
        n=$(diff <(zcat "$a" | grep -v '^#') <(zcat "$b" | grep -v '^#') | grep -c '^[<>]')
        out "$stage" "$f" records DIFF "$n differing body lines"
      fi
      ha=$(digest vcf_header "$a") || exit 3
      hb=$(digest vcf_header "$b") || exit 3
      [ "$ha" = "$hb" ] && out "$stage" "$f" header PASS "" || out "$stage" "$f" header DIFF "" ;;
    *.recal.table)
      n=$(diff <(grep -v '^#:GATKReport.v1.1:' "$a" | grep -v 'Argument' ) <(grep -v '^#:GATKReport.v1.1:' "$b" | grep -v 'Argument') | grep -c '^[<>]')
      out "$stage" "$f" body "$([ "$n" = 0 ] && echo PASS || echo DIFF)" "$n differing lines (excl. argument table)" ;;
  esac
done
