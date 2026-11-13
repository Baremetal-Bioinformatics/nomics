#!/usr/bin/env bash
# The Mutect2 tumour-normal chain end to end against a sarek run, from the run's two
# recalibrated BAMs to the filtered VCF. Every nomics step reads nomics outputs, never the
# golden's: only the BAMs, the interval groups, the resources and the order of inputs that
# sarek's tasks used come from the run.
#
#   Mutect2 per interval group        VCF records, header but ##GATKCommandLine, stats bytes,
#                                     f1r2 archive (scripts/somatic_f1r2_cmp.py)
#   MergeVcfs, MergeMutectStats       records, header but ##GATKCommandLine; stats bytes
#   LearnReadOrientationModel         priors archive members, bytes (from our f1r2 archives)
#   GetPileupSummaries .. CalculateContamination    scripts/somatic_boundary1_gate.sh
#   FilterMutectCalls                 records, header but ##GATKCommandLine; filteringStats bytes
#
#   scripts/somatic_chain_gate.sh [RUN] [OUT_DIR]     NOMICS=path/to/nomics
#   defaults: RUN=~/omics/runs/somatic/S1-r1, REF=~/omics/ref/GATK.GRCh38
#   CHAIN_REUSE=1 keeps the per-group Mutect2 outputs already in OUT_DIR (the long step) and
#   gates again from there.
# Linux (GNU time). On the Spark, run it under the resource lock. Rules: docs/GATES.md.
set -uo pipefail
REPO=$(cd "$(dirname "$0")/.." && pwd)
bin=${NOMICS:-$REPO/target/release/nomics}
R=${1:-$HOME/omics/runs/somatic/S1-r1}; G=$R/golden
O=${2:-$(mktemp -d "${TMPDIR:-/tmp}/somatic-chain.XXXXXX")}; mkdir -p "$O"
REF=${REF:-$HOME/omics/ref/GATK.GRCh38}; FA=$REF/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta
tumour=HCC1395_WES; normal=HCC1395BL_WES; X=${tumour}_vs_${normal}; P=$X.mutect2
pass=0; fail=0
ok() { pass=$((pass+1)); echo "  identical: $1"; }
bad() { fail=$((fail+1)); echo "  DIFFERS: $1"; }
chk() { if cmp -s "$1" "$2"; then ok "$3"; else bad "$3"; fi; }
# A VCF's records, and its header without the command lines.
recs() { gzip -dc "$1" | grep -v '^#'; }
head_() { gzip -dc "$1" | grep '^#' | grep -v '^##GATKCommandLine='; }
vcf_chk() {
  if cmp -s <(recs "$1") <(recs "$2"); then ok "$3 records ($(recs "$1" | wc -l | tr -d ' '))"
  else bad "$3 records ($(comm -12 <(recs "$1" | sort) <(recs "$2" | sort) | wc -l | tr -d ' ') of $(recs "$2" | wc -l | tr -d ' ') lines in common)"; fi
  if cmp -s <(head_ "$1") <(head_ "$2"); then ok "$3 header"; else bad "$3 header"; fi
}
# The inputs of a golden task, in the order its command gave them.
order() { grep -o -- "$2 [^ ]*" "$G/$1/.command.sh" | sed "s/^$2 //"; }

echo "== Mutect2, one run per interval group"
for bed in "$G"/intervals/*.bed; do
  g=$(basename "$bed" .bed); s=$(date +%s)
  golden=$(command ls "$G"/MUTECT2_PAIRED-$X*/"$P.$g.vcf.gz")
  if [[ -n ${CHAIN_REUSE:-} && -f $O/$P.$g.vcf.gz && -f $O/$P.$g.f1r2.tar.gz ]]; then
    echo "$g: reusing $O/$P.$g.vcf.gz"
  elif ! /usr/bin/time -v "$bin" somatic mutect2 \
      --input "$R/results/preprocessing/recalibrated/$normal/$normal.recal.bam" \
      --input "$R/results/preprocessing/recalibrated/$tumour/$tumour.recal.bam" \
      --normal-sample HCC1395_$normal --reference "$FA" --intervals "$bed" \
      --panel-of-normals "$REF/Annotation/GATKBundle/1000g_pon.hg38.vcf.gz" \
      --germline-resource "$REF/Annotation/GATKBundle/af-only-gnomad.hg38.vcf.gz" \
      --output "$O/$P.$g.vcf.gz" --f1r2-tar-gz "$O/$P.$g.f1r2.tar.gz" 2> "$O/$g.log"; then
    bad "Mutect2 $g failed: $(grep -v '^\s' "$O/$g.log" | tail -n 1)"; continue
  else
    echo "$g: $(( $(date +%s) - s )) s; peak $(awk '/Maximum resident/{printf "%.1f GiB", $NF/1048576}' "$O/$g.log")"
  fi
  vcf_chk "$O/$P.$g.vcf.gz" "$golden" "Mutect2 $g"
  chk "$O/$P.$g.vcf.gz.stats" "$golden.stats" "Mutect2 $g stats"
  if python3 -I "$REPO/scripts/somatic_f1r2_cmp.py" "${golden%.vcf.gz}.f1r2.tar.gz" "$O/$P.$g.f1r2.tar.gz" "$g" > "$O/$g.f1r2.txt" 2>&1
  then ok "Mutect2 $g f1r2 archive"; else bad "Mutect2 $g f1r2 archive: $(tail -n 1 "$O/$g.f1r2.txt")"; fi
done

echo "== MergeVcfs, MergeMutectStats"
args=(); while IFS= read -r f; do args+=(--input "$O/$f"); done < <(order "MERGE_MUTECT2-$X" --INPUT)
"$bin" somatic merge-vcfs "${args[@]}" --sequence-dictionary "${FA%.fasta}.dict" --output "$O/$P.vcf.gz" 2> "$O/merge.log" \
  || bad "merge-vcfs failed: $(tail -n 1 "$O/merge.log")"
vcf_chk "$O/$P.vcf.gz" "$G/MERGE_MUTECT2-$X/$P.vcf.gz" "merged VCF"
args=(); while IFS= read -r f; do args+=(--stats "$O/$f"); done < <(order "MERGEMUTECTSTATS-$X" --stats)
"$bin" somatic merge-mutect-stats "${args[@]}" --output "$O/$P.vcf.gz.stats" || bad "merge-mutect-stats failed"
chk "$O/$P.vcf.gz.stats" "$G/MERGEMUTECTSTATS-$X/$P.vcf.gz.stats" "merged stats"

echo "== LearnReadOrientationModel, from our f1r2 archives"
L=$G/LEARNREADORIENTATIONMODEL-$X
args=(); while IFS= read -r f; do args+=(--input "$O/$f"); done < <(order "LEARNREADORIENTATIONMODEL-$X" --input)
"$bin" somatic learn-read-orientation-model "${args[@]}" --output "$O/$P.artifactprior.tar.gz" 2> "$O/lrom.log" \
  || bad "learn-read-orientation-model failed: $(tail -n 1 "$O/lrom.log")"
mkdir -p "$O/priors/ours" "$O/priors/golden"
tar xzf "$O/$P.artifactprior.tar.gz" -C "$O/priors/ours"; tar xzf "$L/$P.artifactprior.tar.gz" -C "$O/priors/golden"
for m in "$O"/priors/golden/*; do chk "$O/priors/ours/$(basename "$m")" "$m" "orientation priors $(basename "$m")"; done

echo "== Pileups and contamination"
if NOMICS=$bin REF=$REF bash "$REPO/scripts/somatic_boundary1_gate.sh" "$R" "$O/b1" > "$O/b1.log" 2>&1
then ok "pileup, contamination and segmentation tables ($(tail -n 1 "$O/b1.log"))"
else bad "pileups and contamination ($(tail -n 1 "$O/b1.log"))"; fi

echo "== FilterMutectCalls, on our own files"
F=$G/FILTERMUTECTCALLS-$X
"$bin" somatic filter-mutect-calls --variant "$O/$P.vcf.gz" --stats "$O/$P.vcf.gz.stats" --reference "$FA" \
  --orientation-bias-artifact-priors "$O/$P.artifactprior.tar.gz" \
  --tumor-segmentation "$O/b1/$P.segmentation.table" --contamination-table "$O/b1/$P.contamination.table" \
  --output "$O/$P.filtered.vcf.gz" 2> "$O/filter.log" || bad "filter-mutect-calls failed: $(tail -n 1 "$O/filter.log")"
vcf_chk "$O/$P.filtered.vcf.gz" "$F/$P.filtered.vcf.gz" "filtered VCF"
chk "$O/$P.filtered.vcf.gz.filteringStats.tsv" "$F/$P.filtered.vcf.gz.filteringStats.tsv" "filtering statistics"
echo "  PASS records: $(recs "$O/$P.filtered.vcf.gz" | awk -F'\t' '$7=="PASS"' | wc -l | tr -d ' ') ours, $(recs "$F/$P.filtered.vcf.gz" | awk -F'\t' '$7=="PASS"' | wc -l | tr -d ' ') golden"

echo "somatic chain gate: $pass identical, $fail differ (outputs in $O)"
(( fail == 0 && pass > 0 ))
