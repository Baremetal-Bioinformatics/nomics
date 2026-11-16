#!/usr/bin/env bash
# How fast does the Spark download from each source the sample matrix can use?
# Reads up to 300 MB from one large file per source, for at most 90 s each, and prints MB/s.
set -uo pipefail
probe() {  # probe NAME URL
  local out
  out=$(curl -sS -o /dev/null --max-time 90 -r 0-300000000 -w '%{size_download} %{time_total} %{http_code}' "$2" 2>/dev/null)
  [[ $out =~ ^[0-9]+\ [0-9.]+\ [0-9]+$ ]] || out="0 0 err"
  awk -v n="$1" -v o="$out" 'BEGIN { split(o, a, " "); printf "%-20s http %-4s %6.1f MB in %5.1f s  %6.2f MB/s\n", n, a[3], a[1]/1e6, a[2], a[1]/1e6/(a[2]+0.001) }'
}
probe "ENA https"   https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/083/SRR14724483/SRR14724483_1.fastq.gz
probe "ENA ftp"     ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/083/SRR14724483/SRR14724483_1.fastq.gz
probe "Google GCS"  https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/novaseq/wes_agilent/50x/HG001.novaseq.wes_agilent.50x.R1.fastq.gz
probe "AWS S3 (us-east-1)" https://ngi-igenomes.s3.amazonaws.com/test-data/sarek/SRR7890919_WES_HCC1395BL-EA_normal_1.fastq.gz
probe "NCBI SRA (S3)" https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR14724483/SRR14724483
