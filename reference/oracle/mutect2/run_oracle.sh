#!/usr/bin/env bash
# Test oracle only. Reproduce the stage goldens before trusting internal dumps.
set -euo pipefail
root=${NOMICS_M2_GOLDEN:-$HOME/omics/runs/somatic/m2-chr21}
ref=${NOMICS_M2_REFERENCE:-$HOME/omics/ref/GATK.GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta}
jar=${GATK_JAR:-$HOME/omics/oracle/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar}
out=${1:?usage: run_oracle.sh OUTSIDE_REPO_OUTPUT_DIRECTORY}
mkdir -p "$out"
java=${JAVA:-java}; javac=${JAVAC:-javac}
args=(--input "$root/HCC1395BL_WES.chr21.bam" --input "$root/HCC1395_WES.chr21.bam"
 --reference "$ref" --panel-of-normals "$root/1000g_pon.chr21.vcf.gz"
 --germline-resource "$root/af-only-gnomad.chr21.vcf.gz" --intervals "$root/chr21.bed"
 --normal-sample HCC1395_HCC1395BL_WES)
"$java" -Xmx4g -XX:-UsePerfData -jar "$jar" Mutect2 "${args[@]}" --output "$out/X.vcf.gz" \
 --tmp-dir "$out" --f1r2-tar-gz "$out/X.f1r2.tar.gz" --assembly-region-out "$out/regions.igv" > "$out/gatk.log" 2>&1
cmp "$out/regions.igv" "$root/oracle/regions.igv"
cmp "$out/X.vcf.gz.stats" "$root"/r1/*.vcf.gz.stats
python3 - "$root" "$out" <<'PY'
import gzip, pathlib, sys
root,out=map(pathlib.Path,sys.argv[1:])
def vcf(p):
    return [l for l in gzip.open(p) if not l.startswith(b'##GATKCommandLine=')]
assert vcf(out/'X.vcf.gz') == vcf(next((root/'r1').glob('*.vcf.gz')))
def counts(p):
    text=p.read_text(); start=text.index('308376 read(s) filtered by:')
    return text[start:].split(' reads processed',1)[0]+' reads processed'
assert counts(out/'gatk.log') == counts(root/'r1/gatk.log')
print('Unmodified JVM: regions, stats, filter counts exact; VCF exact except command provenance')
PY
here=$(cd "$(dirname "$0")" && pwd)
"$javac" -cp "$jar" -d "$out" "$here/MutectRegions.java" "$here/MutectNumerics.java"
"$java" -Xmx4g -XX:-UsePerfData -cp "$out:$jar" MutectRegions "${args[@]}" \
 --dump-prefix "$out/dump" --assembly-region-out "$out/dump.igv" > "$out/dump.log" 2>&1
cmp "$out/dump.igv" "$out/regions.igv"
cmp "$out/dump.stats" "$out/X.vcf.gz.stats"
"$java" -cp "$out:$jar" MutectNumerics > "$out/numerics_oracle.tsv"
echo 'Harness: stage goldens exact. Dumps and numerical fixtures ready.'
