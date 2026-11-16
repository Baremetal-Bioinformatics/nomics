#!/usr/bin/env bash
# Test-only JVM qualification and pre-annotation stage-four dumps.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
out=${1:?usage: run_calls_oracle.sh OUTSIDE_REPO_OUTPUT_DIRECTORY [CONTIG ...]}
shift; contigs=("$@"); (( ${#contigs[@]} )) || contigs=(chr21)
root=${NOMICS_M2_CALLS_GOLDENS:-$HOME/omics/runs/somatic}
ref=${NOMICS_M2_REFERENCE:-$HOME/omics/ref/GATK.GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta}
jar=${GATK_JAR:-$HOME/omics/oracle/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar}
java=${JAVA:-/opt/homebrew/opt/openjdk@17/bin/java}; javac=${JAVAC:-/opt/homebrew/opt/openjdk@17/bin/javac}
JAVA="$java" JAVAC="$javac" bash "$here/run_oracle.sh" "$out/qualification"
"$javac" -cp "$jar" -d "$out" "$here/calls/Mutect2.java" "$here/MutectSomaticNumerics.java" "$here/MutectFragments.java"
"$java" -cp "$out:$jar" MutectSomaticNumerics > "$out/evidence.tsv"
cmp "$out/evidence.tsv" "$here/../../../crates/nomics-somatic/src/mutect2_calls/evidence.tsv"
"$java" -cp "$out:$jar" MutectFragments > "$out/fragments.tsv"
cmp "$out/fragments.tsv" "$here/../../../crates/nomics-somatic/src/mutect2_calls/fragments.tsv"
for c in "${contigs[@]}"; do
 d=$root/m2-$c
 /usr/bin/time -p "$java" -Xmx4g -XX:-UsePerfData -cp "$out:$jar" nomics.oracle.Mutect2 \
  --input "$d/HCC1395BL_WES.$c.bam" --input "$d/HCC1395_WES.$c.bam" \
  --normal-sample HCC1395_HCC1395BL_WES --reference "$ref" --intervals "$d/$c.bed" \
  --panel-of-normals "$d/1000g_pon.$c.vcf.gz" --germline-resource "$d/af-only-gnomad.$c.vcf.gz" \
  --output "$out/$c.vcf.gz" --calls-dump "$out/$c.calls.tsv" \
  --f1r2-tar-gz "$out/$c.f1r2.tar.gz" --bam-output "$out/$c.bam" \
  --bam-writer-type ALL_POSSIBLE_HAPLOTYPES > "$out/$c.log" 2>&1
 python3 - "$d" "$out" "$c" <<'PY'
import gzip,pathlib,sys
root,out=map(pathlib.Path,sys.argv[1:3]); c=sys.argv[3]
golden=next((root/'oracle').glob('*.vcf.gz'))
def vcf(p): return [l for l in gzip.open(p) if not l.startswith(b'##GATKCommandLine=')]
assert vcf(out/f'{c}.vcf.gz')==vcf(golden), 'instrumented VCF differs'
assert (out/f'{c}.vcf.gz.stats').read_bytes()==pathlib.Path(str(golden)+'.stats').read_bytes(), 'instrumented stats differ'
print(f'PASS instrumented {c}: VCF except command provenance, stats bytes; actual/projected LOD assertions executed')
PY
done
