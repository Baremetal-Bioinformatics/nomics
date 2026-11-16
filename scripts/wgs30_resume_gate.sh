#!/usr/bin/env bash
# Resume the 30x run at a later stage in a fresh directory and gate every
# artifact against the committed 30x golden manifest (41 digests): iterate on
# markdup / BQSR / HaplotypeCaller at whole-genome scale without re-aligning.
#
#   scripts/wgs30_resume_gate.sh <markdup|bqsr|hc> [output-dir]   (Spark, under the lock)
#
# markdup needs the source run's mapped/bam.manifest, written when its
# alignment finished. A run directory from before that has none; rebuild the
# chunk set elsewhere (gated chunk by chunk against the source run, which is
# left untouched) and resume from that:
#   CHUNKS=1,2,3,4,5,6,7,8,9,10,11,12 scripts/wgs30_align_gate.sh DIR
#   WGS30_RUN=DIR scripts/wgs30_resume_gate.sh markdup
#
# NOMICS_EXTRA_ARGS: more `nomics run` flags, word-split (e.g.
# "--save-output-as-bam false --keep-bams" for sarek's default CRAM outputs).
set -euo pipefail
FROM=${1:?usage: wgs30_resume_gate.sh <markdup|bqsr|hc> [output-dir]}
REPO=$(cd "$(dirname "$0")/.." && pwd)
BIN=$REPO/target/release/nomics
REF=$HOME/omics/ref/GATK.GRCh38
B=$REF/Annotation/GATKBundle
P=HG002_WGS30
SRC=${WGS30_RUN:-$HOME/omics/runs/nomics-wgs30-20261115}
OUT=${2:-$HOME/omics/runs/wgs30-resume-$FROM-$(date +%Y%m%d-%H%M)}
mkdir -p "$OUT"
ln -sfn "$SRC/prep" "$OUT/prep"
ln -sfn "$SRC/mapped" "$OUT/mapped"
read -r -a EXTRA <<< "${NOMICS_EXTRA_ARGS:-}"
# A resume checks each BAM's generation marker (<bam>.gen), so it travels with the link.
gen() { [ ! -e "$SRC/$P.$1.gen" ] || ln -sf "$SRC/$P.$1.gen" "$OUT/$P.$1.gen"; }
case $FROM in
  bqsr) ln -sf "$SRC/$P.md.bam" "$OUT/$P.md.bam"; gen md.bam; cp "$SRC/$P.md.metrics" "$OUT/" ;;
  hc) for f in md.bam recal.bam; do ln -sf "$SRC/$P.$f" "$OUT/$P.$f"; gen $f; done; cp "$SRC/$P.md.metrics" "$SRC/$P.recal.table" "$OUT/" ;;
esac
"$BIN" run --resume-from "$FROM" --fastq1 /dev/null --fastq2 /dev/null --reference "$REF" \
  --targets "$REF/Annotation/intervals/wgs_calling_regions_noseconds.hg38.bed" --outdir "$OUT" \
  --dbsnp "$B/dbsnp_146.hg38.vcf.gz" --known-sites "$B/dbsnp_146.hg38.vcf.gz" \
  --known-sites "$B/beta/Homo_sapiens_assembly38.known_indels.vcf.gz" \
  --known-sites "$B/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" \
  --known-sites-cache "$HOME/omics/cache/known_sites.GATK.GRCh38.nks" \
  --dense-sa-cache "$HOME/omics/cache/sa4.GRCh38.nsa" \
  --scatter-order "$REPO/reference/wgs30-x86-r1.scatter-order.txt" \
  --read-group-ds /data/omics/ref/GATK.GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta \
  --sample HG002_HG002_WGS30 --prefix $P --threads ${NOMICS_THREADS:-20} ${EXTRA[@]+"${EXTRA[@]}"} > "$OUT/run.log" 2>&1
"$BIN" dev digest --run-dir "$OUT" --prefix $P --output "$OUT/digests.tsv" > "$OUT/digest.log" 2>&1
"$BIN" dev digest-compare --ours "$OUT/digests.tsv" --golden "$REPO/reference/wgs30-x86-r1.digest.tsv" > "$OUT/gates.txt" 2>&1 || true
echo "== $OUT"; column -t "$OUT/stages.tsv"; grep -c "^PASS" "$OUT/gates.txt" || true; grep -E "^FAIL|verdict" "$OUT/gates.txt"
