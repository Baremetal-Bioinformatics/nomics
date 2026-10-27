#!/usr/bin/env bash
# Multi-group recal.cram gate on tier0 (the 30x case in small): a `nomics run
# --save-output-as-bam false` over targets that split into several interval
# groups, against sarek 3.10.0's CRAM-mode chain built from that run's own
# md.bam and recal table (built on an x86 cloud instance):
#   samtools 1.21 `view -Ch` md.bam -> md.cram; GATK 4.6.2.0 ApplyBQSR per
#   group from md.cram (htsjdk 4.2.0 CRAM); samtools 1.24 `merge --output-fmt
#   cram,version=3.0 -s <the run's merge seed>` in the run's merge order.
# Then scripts/cram_gate.sh on our recal.cram (OURS) and on `nomics dev
# bam-to-cram --groups N` of our recal.bam (dev mode, 1-vs-N threads).
#
#   scripts/cram_merge_gate.sh [out-dir] [run-extra]
#
# run-extra: e.g. --spill-prep (streaming; default) or "" (in memory).
set -euo pipefail
REPO=$(cd "$(dirname "$0")/.." && pwd)
OUT=${1:-$HOME/omics/runs/cram-merge-gate-$(date +%Y%m%d-%H%M)}
EXTRA=${2---spill-prep}
REF=$HOME/omics/ref/GATK.GRCh38
FA=$REF/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta
BUNDLE=$REF/Annotation/GATKBundle
DBSNP=$BUNDLE/dbsnp_146.hg38.vcf.gz
CACHE=$HOME/omics/cache
P=HG002_WES
RG_DS=${NOMICS_RG_DS:-${NOMICS_REF:-$HOME/omics/ref/GATK.GRCh38}/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta}
GATK_IMG=community.wave.seqera.io/library/gatk4_htslib_samtools:6fd46722be032bdd # GATK 4.6.2.0 / htsjdk 4.2.0, samtools 1.21
SAM_IMG=community.wave.seqera.io/library/htslib_samtools:1.24--1523e39b702f9721
THREADS=${NOMICS_THREADS:-$(sysctl -n hw.ncpu)}
mkdir -p "$OUT"; OUT=$(cd "$OUT" && pwd)
export CARGO_TARGET_DIR=${NOMICS_BUILD_DIR:-$HOME/omics/build/cpu}
cargo build --manifest-path "$REPO/Cargo.toml" --release -q -p nomics
BIN=$CARGO_TARGET_DIR/release/nomics

# Targets in three interval groups (~120 Mb each at sarek's 200,000 nt/s).
T=$OUT/targets.bed
printf 'chr1\t0\t125000000\nchr1\t125000000\t248956422\nchr20\t0\t64444167\n' > "$T"
echo "== nomics run (CRAM mode, --keep-bams $EXTRA) -> $OUT" >&2
"$BIN" run \
  --fastq1 "$HOME/omics/data/tier0/tier0_1.fastq.gz" --fastq2 "$HOME/omics/data/tier0/tier0_2.fastq.gz" \
  --reference "$REF" --targets "$T" --outdir "$OUT" \
  --dbsnp "$DBSNP" --known-sites "$DBSNP" \
  --known-sites "$BUNDLE/beta/Homo_sapiens_assembly38.known_indels.vcf.gz" \
  --known-sites "$BUNDLE/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" \
  --known-sites-cache "$CACHE/known_sites.GATK.GRCh38.nks" \
  --dense-sa-cache "$CACHE/sa4.GATK.GRCh38.nsa" --dense-sa-interval 4 \
  --read-group-ds "$RG_DS" --sample HG002_$P --prefix $P --threads "$THREADS" \
  --save-output-as-bam false --keep-bams $EXTRA > "$OUT/run.log" 2>&1
SEED=$(awk '/^scatter: .* merge seed/ { print $NF }' "$OUT/run.log")
# Merge order: the groups in BED order (no --scatter-order), each named by
# its first interval.
ORDER=()
while IFS=$'\t' read -r c s e; do
  g="${c}_$((s + 1))-$e"; [[ -f $OUT/intervals/$g.bed ]] && ORDER+=("$g")
done < "$T"
echo "== ${#ORDER[@]} groups (${ORDER[*]}), merge seed $SEED" >&2
[[ ${#ORDER[@]} -ge 2 ]] || { echo "targets gave ${#ORDER[@]} group(s)" >&2; exit 1; }

CH=$OUT/chain; mkdir -p "$CH"
for f in Homo_sapiens_assembly38.fasta Homo_sapiens_assembly38.fasta.fai Homo_sapiens_assembly38.dict; do
  ln -sf "$REF/Sequence/WholeGenomeFasta/$f" "$CH/$f"
done
ln -sf "$OUT/$P.md.bam" "$CH/$P.md.bam"; ln -sf "$OUT/$P.recal.table" "$CH/$P.recal.table"
for g in "${ORDER[@]}"; do cp "$OUT/intervals/$g.bed" "$CH/$g.bed"; done
D=(docker run --rm --platform linux/arm64 -u "$(id -u):$(id -g)" -v "$HOME/omics:$HOME/omics" -w "$CH")
"${D[@]}" "$GATK_IMG" bash -c "samtools view -Ch -T Homo_sapiens_assembly38.fasta -o $P.md.cram $P.md.bam && samtools index $P.md.cram"
ins=()
for i in "${!ORDER[@]}"; do
  g=${ORDER[$i]}
  "${D[@]}" "$GATK_IMG" gatk --java-options '-Xmx3276M -XX:-UsePerfData' ApplyBQSR \
    --input $P.md.cram --output ${P}_$g.recal.cram --reference Homo_sapiens_assembly38.fasta \
    --bqsr-recal-file $P.recal.table --intervals $g.bed --tmp-dir . > "$CH/applybqsr.$g.log" 2>&1
  mkdir -p "$CH/$((i + 1))"; mv "$CH/${P}_$g.recal.cram" "$CH/$((i + 1))/"
  ins+=("$((i + 1))/${P}_$g.recal.cram")
done
"${D[@]}" "$SAM_IMG" samtools merge -s "$SEED" --threads 1 --output-fmt cram,version=3.0 \
  --reference Homo_sapiens_assembly38.fasta $P.recal.cram "${ins[@]}"
"${D[@]}" "$SAM_IMG" samtools index $P.recal.cram

OURS_RECAL=$OUT/preprocessing/recalibrated/$P/$P.recal.cram
samtools index "$OUT/$P.recal.bam" # region counts vs the BAM
OURS=$OURS_RECAL BASELINE=0 NOMICS=$BIN "$REPO/scripts/cram_gate.sh" "$OUT/$P.recal.bam" "$FA" \
  "$CH/$P.recal.cram" sarek-recal "$OUT/gate-run" > "$OUT/gate-run.txt" 2>&1 || true
CRAM_GROUPS=${#ORDER[@]} BASELINE=0 NOMICS=$BIN "$REPO/scripts/cram_gate.sh" "$OUT/$P.recal.bam" "$FA" \
  "$CH/$P.recal.cram" sarek-recal "$OUT/gate-dev" > "$OUT/gate-dev.txt" 2>&1 || true
for g in run dev; do
  echo "== gate ($g)"; grep -E "^(PASS|FAIL|==)" "$OUT/gate-$g.txt"
done
