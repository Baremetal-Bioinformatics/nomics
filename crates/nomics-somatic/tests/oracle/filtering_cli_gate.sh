#!/bin/sh
# From repository root; no reference tool is invoked.
set -eu
G=${NOMICS_SOMATIC_GOLDENS:-"$HOME/omics/runs/somatic/S1-r1/golden"}
R=${NOMICS_SOMATIC_FASTA:-"$HOME/omics/ref/GATK.GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"}
X=HCC1395_WES_vs_HCC1395BL_WES
O=${1:-/tmp/nomics-filter-verified.vcf.gz}
/usr/bin/time -p target/release/nomics somatic filter-mutect-calls \
  --variant "$G/MERGE_MUTECT2-$X/$X.mutect2.vcf.gz" \
  --stats "$G/MERGEMUTECTSTATS-$X/$X.mutect2.vcf.gz.stats" \
  --reference "$R" \
  --orientation-bias-artifact-priors "$G/LEARNREADORIENTATIONMODEL-$X/$X.mutect2.artifactprior.tar.gz" \
  --tumor-segmentation "$G/CALCULATECONTAMINATION-$X/$X.mutect2.segmentation.table" \
  --contamination-table "$G/CALCULATECONTAMINATION-$X/$X.mutect2.contamination.table" \
  --output "$O"
E="$G/FILTERMUTECTCALLS-$X/$X.mutect2.filtered.vcf.gz"
cmp "$E.filteringStats.tsv" "$O.filteringStats.tsv"
printf 'stats cmp: PASS (%s bytes)\n' "$(wc -c < "$O.filteringStats.tsv" | tr -d ' ')"
python3 - "$E" "$O" <<'PY'
import gzip, pathlib, struct, sys
expected, actual = [gzip.open(p, 'rb').readlines() for p in sys.argv[1:]]
records = lambda lines: [s for s in lines if not s.startswith(b'#')]
e, a = records(expected), records(actual)
identical = sum(x == y for x, y in zip(e, a))
passing = sum(s.split(b'\t')[6] == b'PASS' for s in a)
print(f'records: {identical}/{len(e)} identical; PASS: {passing}')
for i, (x, y) in enumerate(zip(e, a)):
    if x != y:
        print(f'first mismatch record {i + 1}\nREF {x}\nOUT {y}')
        break
assert len(a) == len(e) == identical == 6351
assert passing == 2024
header = lambda lines: [s for s in lines if s.startswith(b'#') and not s.startswith(b'##GATKCommandLine=')]
assert header(expected) == header(actual)
print('header: PASS (only ##GATKCommandLine excluded)')
# Check every member uses BGZF's BC/BSIZE field, including the EOF member.
b = pathlib.Path(sys.argv[2]).read_bytes()
p = members = 0
while p < len(b):
    assert b[p:p+4] == b'\x1f\x8b\x08\x04'
    assert b[p+10:p+16] == b'\x06\x00BC\x02\x00'
    size = struct.unpack_from('<H', b, p+16)[0] + 1
    assert 26 <= size <= 65536
    p += size
    members += 1
assert p == len(b)
print(f'BGZF: PASS ({members} members); tabix index: not emitted')
PY
