#!/usr/bin/env bash
# End-to-end nomics run on the Spark: raw FASTQ pair -> VCF in one process,
# then every stage boundary gated against the sarek 3.10.0 reference run.
#
#   scripts/full_exome.sh [tier1|tier1joint|tier0|tier0nv|wgs20|wgs30] [output-dir]   (on the Spark, under nohup)
#
# Build first (from ~/omics/new_omics):
#   cargo build --release --features cuda -p nomics
set -euo pipefail
tier=${1:-tier1}
REPO=$(cd "$(dirname "$0")/.." && pwd)
BIN=$REPO/target/release
REF=$HOME/omics/ref/GATK.GRCh38
BUNDLE=$REF/Annotation/GATKBundle
DBSNP=$BUNDLE/dbsnp_146.hg38.vcf.gz
GOLD=$HOME/omics/runs/$tier-wes-20261008-r1
BED=$HOME/omics/data/Agilent_v7.bed
P=HG002_WES
SCATTER_ORDER=
RG_DS=
GOLD_VCF=
GOLD_JOINT=
case $tier in
  tier1) FQ1=$HOME/omics/data/fastq/SRR14724482_1.fastq.gz; FQ2=$HOME/omics/data/fastq/SRR14724482_2.fastq.gz ;;
  # sarek --joint_germline: --emit-gvcf, then GenotypeGVCFs. Upstream stages gate against tier1;
  # the gVCF and joint VCF against the x86 joint golden (its gVCF matches the ARM one record for record).
  tier1joint) FQ1=$HOME/omics/data/fastq/SRR14724482_1.fastq.gz; FQ2=$HOME/omics/data/fastq/SRR14724482_2.fastq.gz
         GOLD=$HOME/omics/runs/tier1-wes-20261008-r1; GOLD_JOINT=$HOME/omics/runs/tier1-joint-x86-20261016/r1 ;;
  tier0) FQ1=$HOME/omics/data/tier0/tier0_1.fastq.gz; FQ2=$HOME/omics/data/tier0/tier0_2.fastq.gz ;;
  # tier0 with NovaSeq-style names in two lanes: poly-G trimming, optical duplicates, per-lane read groups.
  tier0nv) FQ1=($HOME/omics/data/tier0nv/tier0nv_L{1,2}_1.fastq.gz); FQ2=($HOME/omics/data/tier0nv/tier0nv_L{1,2}_2.fastq.gz)
         GOLD=$HOME/omics/runs/tier0nv-wes-20261012-r1 ;;
  wgs20) FQ1=$HOME/omics/data/wgs_chr20/HG002_WGS_chr20_1.fastq.gz; FQ2=$HOME/omics/data/wgs_chr20/HG002_WGS_chr20_2.fastq.gz
         GOLD=$HOME/omics/runs/wgs20-20261010-r1; BED=$REF/Annotation/intervals/wgs_calling_regions_noseconds.hg38.bed; P=HG002_WGS20
         SCATTER_ORDER=$REPO/reference/wgs20-20261010-r1.scatter-order.txt ;;
  wgs20x86) FQ1=$HOME/omics/data/wgs_chr20/HG002_WGS_chr20_1.fastq.gz; FQ2=$HOME/omics/data/wgs_chr20/HG002_WGS_chr20_2.fastq.gz
         GOLD=$HOME/omics/runs/wgs20-x86-r1; BED=$REF/Annotation/intervals/wgs_calling_regions_noseconds.hg38.bed; P=HG002_WGS20
         SCATTER_ORDER=$REPO/reference/wgs20-x86-r1.scatter-order.txt
         RG_DS=/data/omics/ref/GATK.GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta ;;
  # 30x slice: the first 1/12 of the raw 30x pairs, whole-genome intervals; ~12 min. Golden from x86
  # sarek with HaplotypeCaller in ARM-equivalent Java mode (reference/cloud-hcjava.config).
  wgs30s) FQ1=$HOME/omics/data/wgs30_slice/wgs30s_1.fastq.gz; FQ2=$HOME/omics/data/wgs30_slice/wgs30s_2.fastq.gz
         GOLD=$HOME/omics/runs/wgs30s-x86-20261014-r1; BED=$REF/Annotation/intervals/wgs_calling_regions_noseconds.hg38.bed; P=HG002_WGS30S
         SCATTER_ORDER=$REPO/reference/wgs30s-x86-20261014-r1.scatter-order.txt
         RG_DS=/data/omics/ref/GATK.GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta ;;
  # 30x whole genome; golden VCF is the Java-mode HaplotypeCaller merge (docs/REFERENCE-VARIABILITY.md).
  wgs30) FQ1=$HOME/omics/data/sra_wgs/fastq/SRR14724532_1.fastq.gz; FQ2=$HOME/omics/data/sra_wgs/fastq/SRR14724532_2.fastq.gz
         GOLD=$HOME/omics/runs/wgs30-x86-r1; BED=$REF/Annotation/intervals/wgs_calling_regions_noseconds.hg38.bed; P=HG002_WGS30
         SCATTER_ORDER=$REPO/reference/wgs30-x86-r1.scatter-order.txt
         RG_DS=/data/omics/ref/GATK.GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta
         GOLD_VCF=$HOME/omics/runs/wgs30-x86-hcjava/HG002_WGS30.haplotypecaller.vcf.gz ;;
esac
OUT=${2:-$HOME/omics/runs/nomics-full-$tier-$(date +%Y%m%d-%H%M)}
mkdir -p "$OUT"
echo "== run: raw FASTQs -> VCF" >&2
# NOMICS_CPU_ONLY=1: no GPU seeding (builds without --features cuda).
# NOMICS_RUN_EXTRA: extra `nomics run` flags, e.g. --spill-prep.
GPU_ENV=(NOMICS_GPU_SEED=1); [ -n "${NOMICS_CPU_ONLY:-}" ] && GPU_ENV=()
# One --fastq1/--fastq2 pair per lane (FQ1/FQ2 may be arrays).
# The fastp chunk FASTQs are only written for the gate (sarek does not publish them);
# NOMICS_NO_INTERMEDIATES=1 measures a production run, and their gates then fail.
WRITE_INTERMEDIATES=(--write-intermediates); [ -n "${NOMICS_NO_INTERMEDIATES:-}" ] && WRITE_INTERMEDIATES=()
# Free page cache first so the index and GPU buffers get 2 MiB pages: a
# fragmented run once took 500 s instead of 345 s (NOMICS_NO_EVICT=1 skips).
# Then compact memory where allowed (needs a sudoers rule allowing `tee /proc/sys/vm/compact_memory`;
# skipped otherwise): eviction frees pages but does not defragment them.
[ -n "${NOMICS_NO_EVICT:-}" ] || { python3 "$REPO/scripts/evict_caches.py" "$HOME/omics/data" "$HOME/omics/ref" "$HOME/omics/cache" "$HOME/omics/runs"; sync; echo 1 | sudo -n /usr/bin/tee /proc/sys/vm/compact_memory > /dev/null 2>&1 || true; sleep 10; }
FQ_ARGS=(); for i in "${!FQ1[@]}"; do FQ_ARGS+=(--fastq1 "${FQ1[$i]}" --fastq2 "${FQ2[$i]}"); done
env "${GPU_ENV[@]}" "$BIN/nomics" run \
  "${FQ_ARGS[@]}" --reference "$REF" --targets "$BED" --outdir "$OUT" \
  --dbsnp "$DBSNP" --known-sites "$DBSNP" \
  --known-sites "$BUNDLE/beta/Homo_sapiens_assembly38.known_indels.vcf.gz" \
  --known-sites "$BUNDLE/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" \
  --known-sites-cache "$HOME/omics/cache/known_sites.GATK.GRCh38.nks" \
  --dense-sa-cache "$HOME/omics/cache/sa${NOMICS_SA_INTERVAL:-2}.GRCh38.nsa" --dense-sa-interval "${NOMICS_SA_INTERVAL:-2}" \
  ${SCATTER_ORDER:+--scatter-order "$SCATTER_ORDER"} ${RG_DS:+--read-group-ds "$RG_DS"} \
  --sample HG002_$P --prefix $P --threads 20 ${WRITE_INTERMEDIATES[@]+"${WRITE_INTERMEDIATES[@]}"} ${GOLD_JOINT:+--emit-gvcf} ${NOMICS_RUN_EXTRA:-} \
  > "$OUT/run.log" 2>&1
GOLD_VCF=${GOLD_VCF:-}
[ -n "$GOLD_VCF" ] || [ -n "$GOLD_JOINT" ] || GOLD_VCF=$(ls "$GOLD"/results/variant_calling/haplotypecaller/$P/$P.haplotypecaller.vcf.gz)

# Timing is only comparable with the index and GPU buffers on 2 MiB pages.
TIMING_NOTE=
grep -m1 "index: " "$OUT/run.log" | grep -qvE "index: ([0-9.]+) GB, \1 GB on" && TIMING_NOTE="TIMING SUSPECT: $(grep -m1 'index: ' "$OUT/run.log")"
grep -q "on 2 MiB pages (memory fragmented" "$OUT/run.log" && TIMING_NOTE="$TIMING_NOTE; GPU buffers on small pages"
grep -qE "GPU FALLBACK|GPU warm-up failed" "$OUT/run.log" && TIMING_NOTE="$TIMING_NOTE; $(grep -m1 -oE "GPU FALLBACK:[^(]*|GPU warm-up failed" "$OUT/run.log")"
[ -n "$TIMING_NOTE" ] && echo "$TIMING_NOTE" | tee "$OUT/timing_suspect.txt"
# ---- gates (not timed) ----
G=$OUT/gates.txt
: > "$G"
if [ "$tier" = wgs30 ]; then
  # Whole-genome BAMs do not fit in memory for gate-bam: compare streaming
  # digests (same normalizations) against the committed golden manifest.
  "$BIN/nomics" dev digest --run-dir "$OUT" --prefix $P --output "$OUT/digests.tsv" > "$OUT/digest.log" 2>&1
  "$BIN/nomics" dev digest-compare --ours "$OUT/digests.tsv" --golden "$REPO/reference/wgs30-x86-r1.digest.tsv" > "$G" 2>&1
  echo "== $OUT"; column -t "$OUT/stages.tsv"; cat "$G"
  exit 0
fi
gate() { local name=$1; shift; if "$@" > "$OUT/gate_$name.log" 2>&1; then echo "PASS $name" >> "$G"; else echo "FAIL $name (see gate_$name.log)" >> "$G"; fi; }
# One fastp task per lane: gate each task's chunks.
n=0
for f in $(grep -l 'fastp' "$GOLD"/work/*/*/.command.sh); do
  d=$(dirname "$f")
  [ "$(ls "$d" | grep -c fastp.fastq.gz)" -ge 2 ] && ! ls "$d" | grep -q BWAIndex || continue
  n=$((n + 1)); gate "prep_chunks_$n" "$BIN/nomics" dev gate-chunks --ours "$OUT/prep" --golden "$d"
done
for b in "$OUT"/mapped/$P[.-]*.bam; do
  name=$(basename "$b")
  gold=$GOLD/results/preprocessing/mapped/$P/$name
  [ -e "$gold" ] || gold=$(find "$GOLD/work" -name "$name" -type f | head -1)
  gate "mapped_${name%.bam}" "$BIN/nomics" dev gate-bam --ours "$b" --golden "$gold"
done
gate md_bam "$BIN/nomics" dev gate-bam --ours "$OUT/$P.md.bam" --golden "$GOLD/results/preprocessing/markduplicates/$P/$P.md.bam"
gold_metrics=$(find "$GOLD/results/reports" -name '*.md.bam.metrics' | head -1)
gate md_metrics bash -c "diff <(sed '2d;4d' '$OUT/$P.md.metrics') <(sed '2d;4d' '$gold_metrics')"
gate recal_table "$BIN/nomics" dev gate-report --ours "$OUT/$P.recal.table" --golden "$GOLD/results/preprocessing/recal_table/$P/$P.recal.table"
gate recal_bam "$BIN/nomics" dev gate-bam --ours "$OUT/$P.recal.bam" --golden "$GOLD/results/preprocessing/recalibrated/$P/$P.recal.bam"
if [ -n "$GOLD_JOINT" ]; then
  gate gvcf "$BIN/nomics" dev gate-gvcf --ours "$OUT/$P.haplotypecaller.g.vcf.gz" --golden "$GOLD_JOINT/$P/$P.haplotypecaller.g.vcf.gz"
  gate joint_vcf "$BIN/nomics" dev gate-vcf --ours "$OUT/joint_variant_calling/joint_germline.vcf.gz" --golden "$GOLD_JOINT/joint_variant_calling/joint_germline.vcf.gz"
else
  gate vcf "$BIN/nomics" dev gate-vcf --ours "$OUT/$P.haplotypecaller.vcf.gz" --golden "$GOLD_VCF"
fi

echo "== $OUT"; column -t "$OUT/stages.tsv"; cat "$G"
