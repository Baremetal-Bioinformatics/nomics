#!/usr/bin/env bash
# Test-only JVM qualification followed by assembly dumps from the unmodified jar.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
out=${1:?usage: run_assembly_oracle.sh OUTSIDE_REPO_OUTPUT_DIRECTORY}
root=${NOMICS_M2_GOLDEN:-$HOME/omics/runs/somatic/m2-chr21}
ref=${NOMICS_M2_REFERENCE:-$HOME/omics/ref/GATK.GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta}
jar=${GATK_JAR:-$HOME/omics/oracle/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar}
java=${JAVA:-/opt/homebrew/opt/openjdk@17/bin/java}; javac=${JAVAC:-/opt/homebrew/opt/openjdk@17/bin/javac}
JAVA="$java" JAVAC="$javac" bash "$here/run_oracle.sh" "$out"
"$java" -version > "$out/java-version.txt" 2>&1
shasum -a 256 "$jar" > "$out/jar.sha256"
"$javac" -cp "$jar" -d "$out" "$here/MutectAssembly.java" "$here/AdaptivePruner.java" "$here/ArtificialHaplotypes.java"
/usr/bin/time -p "$java" -Xmx4g -XX:-UsePerfData -cp "$out:$jar" MutectAssembly \
 --input "$root/HCC1395BL_WES.chr21.bam" --input "$root/HCC1395_WES.chr21.bam" \
 --reference "$ref" --panel-of-normals "$root/1000g_pon.chr21.vcf.gz" \
 --germline-resource "$root/af-only-gnomad.chr21.vcf.gz" --intervals "$root/chr21.bed" \
 --normal-sample HCC1395_HCC1395BL_WES --dump-prefix "$out/assembly" \
 --assembly-region-out "$out/assembly.igv" > "$out/assembly.log" 2>&1
cmp "$out/assembly.igv" "$root/oracle/regions.igv"
cmp "$out/assembly.stats" "$out/X.vcf.gz.stats"
"$java" -cp "$out:$jar" AdaptivePruner > "$out/pruner_oracle.tsv"
cmp "$out/pruner_oracle.tsv" "$here/../../../crates/nomics-somatic/src/mutect2_assembly/pruner_oracle.tsv"
echo 'Qualified chr21 assembly dump and independent pruner fixtures ready.'
