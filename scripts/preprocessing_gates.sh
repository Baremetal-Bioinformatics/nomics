# Shared by full_exome.sh and somatic_run_gate.sh. Same comparators and normalizations.
preprocessing_gates() {
gate() { local name=$1; shift; if "$@" > "$OUT/gate_$name.log" 2>&1; then echo "PASS $name" >> "$G"; else echo "FAIL $name (see gate_$name.log)" >> "$G"; fi; }
# One fastp task per lane: gate each task's chunks.
n=0
for f in ${PREP_TASK_COMMANDS:-$(grep -l 'fastp' "$GOLD"/work/*/*/.command.sh)}; do
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
gold_metrics=${GOLD_METRICS:-$(find "$GOLD/results/reports" -name '*.md.bam.metrics' | head -1)}
gate md_metrics bash -c "diff <(sed '2d;4d' '$OUT/$P.md.metrics') <(sed '2d;4d' '$gold_metrics')"
gate recal_table "$BIN/nomics" dev gate-report --ours "$OUT/$P.recal.table" --golden "$GOLD/results/preprocessing/recal_table/$P/$P.recal.table"
gate recal_bam "$BIN/nomics" dev gate-bam --ours "$OUT/$P.recal.bam" --golden "$GOLD/results/preprocessing/recalibrated/$P/$P.recal.bam"
}
