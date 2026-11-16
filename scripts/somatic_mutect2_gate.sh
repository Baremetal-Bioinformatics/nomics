#!/usr/bin/env bash
# Gate `nomics somatic mutect2` against the slice goldens (scripts/somatic_mutect2_slice.sh).
# Per contig, from one run of the command (rules in docs/GATES.md): every VCF record by
# bytes, the header except ##GATKCommandLine, the stats file by bytes, the f1r2 archive
# under its rule (scripts/somatic_f1r2_cmp.py), and the orientation priors learned from our
# archive against those learned from GATK's, by bytes. A contig with no golden fails.
# Shell and python3 only.
#
#   scripts/somatic_mutect2_gate.sh [CONTIG ...]   default: chr21 chr22 chr19 chr20 chr13 chr16 chr12
#   NOMICS=path/to/nomics; NOMICS_M2_CALLS_GOLDENS (~/omics/runs/somatic); NOMICS_M2_CALLS_OUT keeps outputs
set -euo pipefail
repo=$(cd "$(dirname "$0")/.." && pwd)
bin=${NOMICS:-$repo/target/release/nomics}
root=${NOMICS_M2_CALLS_GOLDENS:-$HOME/omics/runs/somatic}
ref=${NOMICS_M2_REFERENCE:-$HOME/omics/ref/GATK.GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta}
out=${NOMICS_M2_CALLS_OUT:-$(mktemp -d "${TMPDIR:-/tmp}/nomics-m2-calls.XXXXXX")}
contigs=("$@"); (( ${#contigs[@]} )) || contigs=(chr21 chr22 chr19 chr20 chr13 chr16 chr12)
mkdir -p "$out"
echo "Outputs: $out"
fail=0; tested=0
for c in "${contigs[@]}"; do
 d=$root/m2-$c
 [[ -f $d/HCC1395BL_WES.$c.bam ]] || { echo "FAIL $c: no golden directory"; fail=1; continue; }
 tested=$((tested+1))
 if ! /usr/bin/time -p "$bin" somatic mutect2 \
  --input "$d/HCC1395BL_WES.$c.bam" --input "$d/HCC1395_WES.$c.bam" \
  --normal-sample HCC1395_HCC1395BL_WES --reference "$ref" --intervals "$d/$c.bed" \
  --panel-of-normals "$d/1000g_pon.$c.vcf.gz" --germline-resource "$d/af-only-gnomad.$c.vcf.gz" \
  --output "$out/$c.vcf.gz" --f1r2-tar-gz "$out/$c.f1r2.tar.gz" > "$out/$c.stdout" 2> "$out/$c.log"; then
  echo "FAIL $c: caller failed (see $out/$c.log)"; fail=1; continue
 fi
 if ! python3 - "$d" "$out" "$c" <<'PY'
import gzip,pathlib,sys,collections,itertools
root,out=map(pathlib.Path,sys.argv[1:3]); c=sys.argv[3]
def vcf(p):
    lines=[line.decode('utf-8') for line in gzip.open(p,'rb').readlines()]
    return ([l for l in lines if l.startswith('#') and not l.startswith('##GATKCommandLine=')],
            [l for l in lines if not l.startswith('#')])
golden=next((root/'oracle').glob('*.vcf.gz'))
header,records=vcf(out/f'{c}.vcf.gz'); gh,gr=vcf(golden)
exact=sum(a==b for a,b in zip(records,gr))
print(f'{c}: records identical {exact}/{len(gr)}; emitted {len(records)}')
counts=collections.Counter(); examples={}; seen=collections.Counter()
def fields(line):
    f=line.rstrip('\n').split('\t')
    d=dict(zip(['CHROM','POS','ID','REF','ALT','QUAL','FILTER'],f[:7]))
    d['record/newline']='CRLF' if line.endswith('\r\n') else ('LF' if line.endswith('\n') else '<absent>')
    for entry in f[7].split(';'):
        k,_,v=entry.partition('=');d['INFO/'+k]=v
    d['FORMAT/order']=f[8]
    for i,s in enumerate(f[9:]):
        d.update({f'FORMAT/{i}/{k}':v for k,v in zip(f[8].split(':'),s.split(':'))})
    return d
for i,(a,b) in enumerate(itertools.zip_longest(gr,records)):
    ga=fields(a) if a else {}; our=fields(b) if b else {}
    for k in sorted(ga.keys()|our.keys()):
        seen[k]+=1
        x,y=ga.get(k,'<absent>'),our.get(k,'<absent>')
        if x!=y: counts[k]+=1; examples.setdefault(k,(i+1,x,y))
print('field | differences | first record | golden | ours')
for k,n in counts.most_common():
    i,x,y=examples[k]; print(f'{k} | {n} | {i} | {x} | {y}')
print('identical fields: '+', '.join(k for k in sorted(seen) if not counts[k]))
print('header: '+('PASS' if header==gh else 'FAIL'))
if header!=gh and set(header)==set(gh):
    for i,(a,b) in enumerate(zip(gh,header),1):
        if a!=b: print(f'header order line {i}: golden {a!r}; ours {b!r}')
for l in gh:
    if l not in header: print('header missing: '+l.rstrip())
for l in header:
    if l not in gh: print('header extra: '+l.rstrip())
s=(out/f'{c}.vcf.gz.stats').read_bytes();gs=pathlib.Path(str(golden)+'.stats').read_bytes()
print(f'stats: {"PASS" if s==gs else "FAIL"} ({len(s)} bytes)')
sys.exit(0 if records==gr and header==gh and s==gs else 1)
PY
 then fail=1; fi
 g=$d/oracle/HCC1395_WES_vs_HCC1395BL_WES.mutect2.$c.f1r2.tar.gz
 python3 -I "$repo/scripts/somatic_f1r2_cmp.py" "$g" "$out/$c.f1r2.tar.gz" "$c f1r2" || { echo "FAIL $c f1r2 archive"; fail=1; }
 if "$bin" somatic learn-read-orientation-model --input "$out/$c.f1r2.tar.gz" --output "$out/$c.priors.tar.gz" 2>> "$out/$c.log" \
    && "$bin" somatic learn-read-orientation-model --input "$g" --output "$out/$c.priors.golden.tar.gz" 2>> "$out/$c.log" \
    && cmp -s "$out/$c.priors.tar.gz" "$out/$c.priors.golden.tar.gz"
 then echo "$c orientation priors: PASS ($(wc -c < "$out/$c.priors.tar.gz" | tr -d ' ') bytes)"
 else echo "$c orientation priors: FAIL"; fail=1; fi
 grep '^real ' "$out/$c.log"
done
(( tested > 0 )) || { echo 'No contig goldens found'; exit 1; }
exit "$fail"
