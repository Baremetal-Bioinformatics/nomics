#!/usr/bin/env bash
# Mac-only orchestration test, NOT a sarek reference gate. Existing samtools is
# used solely to create fixtures. Invoke under lockf -k /tmp/nomics-heavy.lock.
# No build, downloads or installs. Every fixture/output stays outside the repo.
set -euo pipefail
REPO=$(cd "$(dirname "$0")/.." && pwd)
D=${1:?usage: somatic_run_plumbing.sh NEW_EXTERNAL_DIRECTORY}
[[ ! -e $D ]] || { echo "refusing existing $D" >&2; exit 1; }
mkdir -p "$D"
D=$(cd "$D" && pwd)
S=$HOME/omics/runs/somatic/m2-chr21
mkdir -p "$D"
printf 'chr21\t10400000\t10500000\n' > "$D/targets.bed"
for sample in HCC1395BL_WES HCC1395_WES; do
 samtools view -b -f 1 -F 2304 "$S/$sample.chr21.bam" chr21:10400001-10500000 > "$D/$sample.bam"
 samtools collate -u -O "$D/$sample.bam" | samtools fastq -n -1 "$D/$sample.all_1.fq" -2 "$D/$sample.all_2.fq" -0 /dev/null -s /dev/null -
 head -n 8000 "$D/$sample.all_1.fq" | gzip -c > "$D/${sample}_1.fq.gz"
 head -n 8000 "$D/$sample.all_2.fq" | gzip -c > "$D/${sample}_2.fq.gz"
done
printf 'patient,sex,status,sample,lane,fastq_1,fastq_2\nHCC1395,XX,0,HCC1395BL_WES,L1,%s/HCC1395BL_WES_1.fq.gz,%s/HCC1395BL_WES_2.fq.gz\nHCC1395,XX,1,HCC1395_WES,L1,%s/HCC1395_WES_1.fq.gz,%s/HCC1395_WES_2.fq.gz\n' "$D" "$D" "$D" "$D" > "$D/samplesheet.csv"

BIN=${NOMICS:-$REPO/target/release/nomics}
REF=$HOME/omics/ref/GATK.GRCh38
S=$HOME/omics/runs/somatic/m2-chr21
ARGS=(--reference "$REF" --targets "$D/targets.bed" --dense-sa-cache "$HOME/omics/cache/sa4.GATK.GRCh38.nsa" --known-sites-cache "$HOME/omics/cache/known_sites.GATK.GRCh38.nks" --threads 6 --write-intermediates)
RES=(--germline-resource "$S/af-only-gnomad.chr21.vcf.gz" --panel-of-normals "$S/1000g_pon.chr21.vcf.gz" --pileup-resource "$S/af-only-gnomad.chr21.vcf.gz")
/usr/bin/time -l -o "$D/pair.time" "$BIN" run "${ARGS[@]}" "${RES[@]}" --samplesheet "$D/samplesheet.csv" --outdir "$D/pair" > "$D/pair.log" 2>&1
for sample in HCC1395BL_WES HCC1395_WES; do
 /usr/bin/time -l -o "$D/$sample.time" "$BIN" run "${ARGS[@]}" --sample HCC1395_$sample --prefix "$sample" --fastq1 "$D/${sample}_1.fq.gz" --fastq2 "$D/${sample}_2.fq.gz" --outdir "$D/$sample" > "$D/$sample.log" 2>&1
done

REF=$HOME/omics/ref/GATK.GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta
S=$HOME/omics/runs/somatic/m2-chr21
for sample in HCC1395BL_WES HCC1395_WES; do
 "$BIN" dev htsjdk-index --input "$D/$sample/$sample.recal.bam" --output "$D/$sample/$sample.recal.bai"
done
/usr/bin/time -l -o "$D/chain.time" "$BIN" somatic chain --normal "$D/HCC1395BL_WES/HCC1395BL_WES.recal.bam" --tumour "$D/HCC1395_WES/HCC1395_WES.recal.bam" --reference "$REF" --intervals-dir "$D/pair/intervals" --orientation-order chr21_10400001-10500000 --normal-label HCC1395BL_WES --tumour-label HCC1395_WES --germline-resource "$S/af-only-gnomad.chr21.vcf.gz" --panel-of-normals "$S/1000g_pon.chr21.vcf.gz" --pileup-resource "$S/af-only-gnomad.chr21.vcf.gz" --outdir "$D/chain" --keep-intermediates > "$D/chain.log" 2>&1
python3 -I "$REPO/scripts/somatic_run_plumbing_compare.py" "$D" "$BIN"
