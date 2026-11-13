#!/usr/bin/env bash
# Exact stage-two checks; Java is a test oracle only. Keep outputs for diagnosis.
set -uo pipefail
repo=$(cd "$(dirname "$0")/.." && pwd)
bin=${NOMICS:-$repo/target/release/nomics}
root=${NOMICS_M2_ASSEMBLY_GOLDENS:-$HOME/omics/runs/somatic}
ref=${NOMICS_M2_REFERENCE:-$HOME/omics/ref/GATK.GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta}
jar=${GATK_JAR:-$HOME/omics/oracle/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar}
java=${JAVA:-/opt/homebrew/opt/openjdk@17/bin/java}; javac=${JAVAC:-/opt/homebrew/opt/openjdk@17/bin/javac}
out=${1:?usage: somatic_mutect2_haplotypes_gate.sh OUTSIDE_REPO_OUTPUT_DIRECTORY [CONTIG ...]}
shift; contigs=("$@"); (( ${#contigs[@]} )) || contigs=(chr21 chr22 chr19 chr20 chr13 chr16 chr12)
mkdir -p "$out"
"$javac" -cp "$jar" -d "$out" "$repo/reference/oracle/mutect2/ArtificialHaplotypes.java" || exit 1
fail=0; tested=0
for c in "${contigs[@]}"; do
  d=$root/m2-$c
  [[ -f $d/oracle/haplotypes.bam ]] || { echo "SKIP $c: no $d/oracle/haplotypes.bam"; continue; }
  tested=$((tested + 1))
  /usr/bin/time -p "$bin" somatic mutect2-haplotypes --input "$d/HCC1395BL_WES.$c.bam" --input "$d/HCC1395_WES.$c.bam" \
    --normal-sample HCC1395_HCC1395BL_WES --reference "$ref" --intervals "$d/$c.bed" \
    --panel-of-normals "$d/1000g_pon.$c.vcf.gz" --germline-resource "$d/af-only-gnomad.$c.vcf.gz" \
    --assembly-region-out "$out/$c.igv" --stats "$out/$c.stats" --output "$out/$c.assembly.tsv" \
    --artificial-records "$out/$c.artificial.tsv" 2> "$out/$c.log" || { fail=1; tail "$out/$c.log"; continue; }
  "$java" -cp "$out:$jar" ArtificialHaplotypes "$d/oracle/haplotypes.bam" > "$out/$c.golden.tsv" || { fail=1; continue; }
  if cmp "$out/$c.artificial.tsv" "$out/$c.golden.tsv"; then echo "PASS $c artificial"; else echo "FAIL $c artificial"; fail=1; fi
  grep -E '^assembly:|^real' "$out/$c.log"
  echo "$c golden artificial records: $(wc -l < "$out/$c.golden.tsv")"
done
if [[ -n ${NOMICS_M2_ASSEMBLY_ORACLE:-} ]]; then
  if cmp "$out/chr21.assembly.tsv" "$NOMICS_M2_ASSEMBLY_ORACLE"; then echo 'PASS chr21 internal'; else echo 'FAIL chr21 internal'; fail=1; fi
else
  echo 'SKIP chr21 internal: set NOMICS_M2_ASSEMBLY_ORACLE to qualified Java dump'
fi
(( tested > 0 )) || { echo "No contig goldens found"; exit 1; }
exit "$fail"
