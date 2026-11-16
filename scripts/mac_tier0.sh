#!/usr/bin/env bash
# tier0 FASTQ -> VCF on Apple Silicon, one arm per call, then every stage
# artifact digested and compared with the sarek tier0 golden manifest
# (reference/tier0-wes-20261008-r1.digest.tsv, made on the Spark with
# `nomics dev digest`), so no golden files are needed on the Mac.
#
#   scripts/mac_tier0.sh cpu|metal [output-dir]
#
# SAVE_OUTPUT_AS_BAM=false runs as default sarek (md.cram/recal.cram under
# <out>/preprocessing) and gates both CRAMs with scripts/cram_gate.sh against
# sarek's own (scripts/cram_reference.sh, ~/omics/ref-cram/tier0); the
# digest then has no md/recal BAM entries. RUN_EXTRA adds `nomics run`
# flags (e.g. --spill-prep for the whole-genome streaming path).
#
# Builds into ~/omics/build/<arm> (CPU: default features; Metal: --features
# metal and NOMICS_GPU_SEED=1, i.e. Metal seeding and mate rescue). Inputs
# mirror the Spark layout under ~/omics (docs/strategy/metal-sprint.md, A8).
set -euo pipefail
arm=${1:?usage: mac_tier0.sh cpu|metal [output-dir]}
REPO=$(cd "$(dirname "$0")/.." && pwd)
REF=$HOME/omics/ref/GATK.GRCh38
BUNDLE=$REF/Annotation/GATKBundle
DBSNP=$BUNDLE/dbsnp_146.hg38.vcf.gz
P=HG002_WES
# bwa's @RG DS is the reference path as the golden run saw it (the Spark's).
RG_DS=${NOMICS_RG_DS:-${NOMICS_REF:-$HOME/omics/ref/GATK.GRCh38}/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta}
case $arm in
  cpu) FEATURES=(); GPU_ENV=() ;;
  metal) FEATURES=(--features metal); GPU_ENV=(NOMICS_GPU_SEED=1) ;;
  *) echo "arm must be cpu or metal, not $arm" >&2; exit 2 ;;
esac
OUT=${2:-$HOME/omics/runs/mac-tier0-$arm-$(date +%Y%m%d-%H%M)}
THREADS=${NOMICS_THREADS:-$(sysctl -n hw.ncpu)}
SAVE_OUTPUT_AS_BAM=${SAVE_OUTPUT_AS_BAM:-true}
CACHE=$HOME/omics/cache
mkdir -p "$OUT" "$CACHE"

export CARGO_TARGET_DIR=${NOMICS_BUILD_DIR:-$HOME/omics/build/$arm}
cargo build --manifest-path "$REPO/Cargo.toml" --release -q -p nomics ${FEATURES[@]+"${FEATURES[@]}"}
BIN=$CARGO_TARGET_DIR/release/nomics

echo "== $arm: tier0 FASTQ -> VCF, $THREADS threads -> $OUT" >&2
start=$(date +%s)
/usr/bin/time -l -o "$OUT/time.txt" env ${GPU_ENV[@]+"${GPU_ENV[@]}"} NOMICS_TIMING=1 "$BIN" run \
  --fastq1 "$HOME/omics/data/tier0/tier0_1.fastq.gz" --fastq2 "$HOME/omics/data/tier0/tier0_2.fastq.gz" \
  --reference "$REF" --targets "$HOME/omics/data/Agilent_v7.bed" --outdir "$OUT" \
  --dbsnp "$DBSNP" --known-sites "$DBSNP" \
  --known-sites "$BUNDLE/beta/Homo_sapiens_assembly38.known_indels.vcf.gz" \
  --known-sites "$BUNDLE/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" \
  --known-sites-cache "$CACHE/known_sites.GATK.GRCh38.nks" \
  --dense-sa-cache "$CACHE/sa4.GATK.GRCh38.nsa" --dense-sa-interval 4 \
  --read-group-ds "$RG_DS" --sample HG002_$P --prefix $P --threads "$THREADS" --write-intermediates \
  --save-output-as-bam "$SAVE_OUTPUT_AS_BAM" ${RUN_EXTRA:-} \
  > "$OUT/run.log" 2>&1
wall=$(( $(date +%s) - start ))
peak=$(awk '/maximum resident set size/ {printf "%.1f", $1 / 1e9}' "$OUT/time.txt")
echo "== $arm: wall ${wall}s, peak RSS ${peak} GB" | tee "$OUT/summary.txt" >&2

"$BIN" dev digest --run-dir "$OUT" --prefix $P --output "$OUT/digests.tsv" > "$OUT/digest.log" 2>&1
"$BIN" dev digest-compare --ours "$OUT/digests.tsv" \
  --golden "$REPO/reference/tier0-wes-20261008-r1.digest.tsv" | tee "$OUT/digest-compare.txt" | tail -3 \
  || [[ $SAVE_OUTPUT_AS_BAM == false ]] # md_bam/recal_bam are then missing by design

if [[ $SAVE_OUTPUT_AS_BAM == false ]]; then
  RC=${CRAM_REFERENCE:-$HOME/omics/ref-cram/tier0}
  G=${GOLDEN:-${NOMICS_FIXTURES:-$HOME/omics-fixtures}/sarek-tier0/r1}/results/preprocessing
  for kind in md recal; do
    sub=$([[ $kind == md ]] && echo markduplicates || echo recalibrated)
    OURS=$OUT/preprocessing/$sub/$P/$P.$kind.cram BASELINE=0 NOMICS=$BIN \
      "$REPO/scripts/cram_gate.sh" "$G/$sub/$P/$P.$kind.bam" "$REF/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta" \
      "$RC/$P.$kind.cram" "sarek-$kind" "$OUT/cram-gate-$kind" | tee "$OUT/cram-gate-$kind.txt" | grep -E '^(PASS|FAIL|==)'
  done
fi
