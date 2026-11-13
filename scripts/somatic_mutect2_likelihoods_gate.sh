#!/usr/bin/env bash
# Exact stage-three gate; Java is only the independent BAM projection oracle.
set -euo pipefail
repo=$(cd "$(dirname "$0")/.." && pwd)
bin=${NOMICS:-$repo/target/release/nomics}
root=${NOMICS_M2_LIKELIHOODS_GOLDENS:-$HOME/omics/runs/somatic}
ref=${NOMICS_M2_REFERENCE:-$HOME/omics/ref/GATK.GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta}
jar=${GATK_JAR:-$HOME/omics/oracle/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar}
java=${JAVA:-/opt/homebrew/opt/openjdk@17/bin/java}; javac=${JAVAC:-/opt/homebrew/opt/openjdk@17/bin/javac}
out=${1:?usage: somatic_mutect2_likelihoods_gate.sh OUTDIR [CONTIG ...]}
shift; contigs=("$@"); (( ${#contigs[@]} )) || contigs=(chr21 chr22 chr19 chr20 chr13 chr16 chr12)
mkdir -p "$out"
"$javac" -cp "$jar" -d "$out" "$repo/reference/oracle/mutect2/BiologicalReads.java"
fail=0; tested=0
for c in "${contigs[@]}"; do
 d=$root/m2-$c
 [[ -f $d/oracle/haplotypes.bam ]] || { echo "SKIP $c: no golden"; continue; }
 tested=$((tested+1))
 /usr/bin/time -p "$bin" somatic mutect2-likelihoods --input "$d/HCC1395BL_WES.$c.bam" --input "$d/HCC1395_WES.$c.bam" \
   --normal-sample HCC1395_HCC1395BL_WES --reference "$ref" --intervals "$d/$c.bed" \
   --panel-of-normals "$d/1000g_pon.$c.vcf.gz" --germline-resource "$d/af-only-gnomad.$c.vcf.gz" \
   --assembly-region-out "$out/$c.igv" --stats "$out/$c.stats" --output "$out/$c.likelihoods.tsv" \
   --realigned-records "$out/$c.realigned.tsv" 2> "$out/$c.log" || { fail=1; tail "$out/$c.log"; continue; }
 "$java" -cp "$out:$jar" BiologicalReads "$d/oracle/haplotypes.bam" > "$out/$c.golden.tsv"
 if cmp "$out/$c.realigned.tsv" "$out/$c.golden.tsv"; then echo "PASS $c biological records"; else echo "FAIL $c biological records"; fail=1; fi
 grep -E '^likelihoods:|^real' "$out/$c.log"
 echo "$c golden records: $(wc -l < "$out/$c.golden.tsv")"
done
if [[ -n ${NOMICS_M2_LIKELIHOODS_ORACLE:-} ]]; then
 if cmp "$out/chr21.likelihoods.tsv" "$NOMICS_M2_LIKELIHOODS_ORACLE"; then echo 'PASS chr21 internal'; else echo 'FAIL chr21 internal'; fail=1; fi
else
 echo 'SKIP chr21 internal: set NOMICS_M2_LIKELIHOODS_ORACLE to qualified Java dump'
fi
(( tested > 0 )) || { echo 'No contig goldens found'; exit 1; }
exit "$fail"
