#!/usr/bin/env bash
# Deterministic: 40 kb deep-pile windows plus greedy call-rich 4 kb windows
# (target 32 golden records/contig, at most 32), plus richest indel/multiallelic
# windows. Ties: lowest start. Clip source BED to the union; at most three contigs.
# SAMTOOLS may name a docker wrapper.
# Usage: somatic_mutect2_interval_cases.sh <slice-name-or-directory> <out-dir>
set -euo pipefail
slice=${1:?slice required}; out=${2:?output directory required}
[[ -d $slice ]] || slice=${NOMICS_M2_CALLS_GOLDENS:-$HOME/omics/runs/somatic}/m2-${slice#m2-}
mkdir -p "$out"
python3 -I - "$slice" "$out" "${SAMTOOLS:-samtools}" <<'PY'
import bisect, collections, gzip, hashlib, json, os, pathlib, random, subprocess, sys
src,out=map(lambda x:pathlib.Path(x).resolve(),sys.argv[1:3]); sam=sys.argv[3]
tmp=pathlib.Path(os.environ.get('TMPDIR',str(out/'tmp'))).resolve()
if not tmp.is_relative_to(out): tmp=out/'tmp'
tmp.mkdir(parents=True,exist_ok=True); os.environ['TMPDIR']=str(tmp)
# Refuse overwriting an existing set: cached reference results must not change underneath it.
if (out/'manifest.json').exists():
    sys.exit('output already contains a manifest; use a fresh directory')
beds=list(src.glob('*.bed'))
if len(beds)!=1: sys.exit(f'expected one slice BED, found {beds}')
bed=beds[0]; name=bed.stem
paths={k:str(src/f'{prefix}.{name}.{suffix}') for k,prefix,suffix in [
    ('normal','HCC1395BL_WES','bam'),('tumour','HCC1395_WES','bam'),
    ('pon','1000g_pon','vcf.gz'),('germline','af-only-gnomad','vcf.gz')]}
for p in paths.values():
    if not pathlib.Path(p).is_file(): sys.exit(f'missing {p}')
header=subprocess.check_output([sam,'view','-H',paths['tumour']],text=True)
lengths={}
for line in header.splitlines():
    if line.startswith('@SQ\t'):
        d=dict(x.split(':',1) for x in line.split('\t')[1:]); lengths[d['SN']]=int(d['LN'])
rank={c:i for i,c in enumerate(lengths)}
rows=[]
for line in bed.read_text().splitlines():
    if not line or line.startswith(('#','track','browser')): continue
    c,s,e,*_=line.split(); s,e=int(s),int(e)
    if c not in lengths or not 0<=s<e<=lengths[c]: sys.exit(f'invalid source interval {line}')
    rows.append((c,s,e))
rows.sort(key=lambda r:(rank[r[0]],r[1],r[2]))
by=collections.defaultdict(list)
for r in rows: by[r[0]].append(r)
# Merge source spans only for membership queries; retain original records in all transforms.
merged={}
for c,rr in by.items():
    mm=[]
    for _,s,e in rr:
        if mm and s<=mm[-1][1]: mm[-1]=(mm[-1][0],max(e,mm[-1][1]))
        else: mm.append((s,e))
    merged[c]=mm
starts={c:[s for s,e in mm] for c,mm in merged.items()}
peaks={}
for sample in ('normal','tumour'):
    # Indexed slice selection; coordinate-sorted BAM means one start group at a time.
    p=subprocess.Popen([sam,'view','-F','3844','-q','20','-L',str(bed),paths[sample]],stdout=subprocess.PIPE,text=True)
    last=None; count=0
    def record(key,n):
        if key is None: return
        c,pos=key
        if c not in merged: return
        i=bisect.bisect_right(starts[c],pos)-1
        if i<0 or pos>=merged[c][i][1]: return
        old=peaks.get(c)
        if old is None or (-n,pos,sample)<(-old['count'],old['start0'],old['sample']):
            peaks[c]={'start0':pos,'count':n,'sample':sample}
    for line in p.stdout:
        f=line.split('\t'); key=(f[2],int(f[3])-1)
        if key!=last: record(last,count); last=key; count=0
        count+=1
    record(last,count)
    if p.wait(): sys.exit(f'samtools scan failed: {sample}')
if set(peaks)!=set(by): sys.exit('at least one BED contig has no eligible covered alignment start')
# Prefer r1, then the oracle layout used by multi; never use nomics outputs.
golden=None
for directory in ('r1','oracle','.'):
    candidates=sorted((src/directory).glob(f'*.mutect2.{name}.vcf.gz'))
    if len(candidates)>1: sys.exit(f'ambiguous golden VCFs: {candidates}')
    if candidates: golden=candidates[0]; break
if golden is None: sys.exit(f'no GATK golden VCF found: expected r1/, oracle/ or root/*.mutect2.{name}.vcf.gz under {src}')
contigs=sorted(by,key=lambda c:rank[c])[:3]
calls={c:[] for c in contigs}
with gzip.open(golden,'rt') as f:
    for line in f:
        if line.startswith('#'): continue
        fields=line.split('\t'); c=fields[0]; pos=int(fields[1])-1
        if c not in calls: continue
        i=bisect.bisect_right(starts[c],pos)-1
        if i<0 or pos>=merged[c][i][1]: continue
        ref,alt=fields[3:5]; alts=alt.split(',')
        calls[c].append({'pos0':pos,'end0':pos+len(ref),'indel':any(len(x)!=len(ref) for x in alts),'multi':len(alts)>1})
for cc in calls.values(): cc.sort(key=lambda v:(v['pos0'],v['end0']))
if not any(calls.values()): sys.exit(f'golden VCF has no records starting in selected source BED contigs: {golden}')
windows=[]
for c in contigs:
    def window(pos): return (max(0,pos-2000),min(lengths[c],pos+2000))
    pos=peaks[c]['start0']; windows.append((c,max(0,pos-20000),min(lengths[c],pos+20000)))
    remaining=calls[c].copy(); covered=0
    for _ in range(32):
        if covered>=min(32,len(calls[c])) or not remaining: break
        lo,hi=min((window(v['pos0']) for v in remaining),key=lambda w:(-sum(w[0]<=v['pos0']<w[1] for v in remaining),w[0]))
        windows.append((c,lo,hi))
        rest=[v for v in remaining if not lo<=v['pos0']<hi]
        covered+=len(remaining)-len(rest); remaining=rest
    for kind in ('indel','multi'):
        candidates=[window(v['pos0']) for v in calls[c] if v[kind]]
        if candidates:
            lo,hi=min(candidates,key=lambda w:(-sum(w[0]<=v['pos0']<w[1] and v[kind] for v in calls[c]),w[0]))
            windows.append((c,lo,hi))
win_union=[]
for c,s,e in sorted(set(windows),key=lambda r:(rank[r[0]],r[1],r[2])):
    if win_union and win_union[-1][0]==c and s<=win_union[-1][2]:
        a,b,z=win_union[-1]; win_union[-1]=(a,b,max(z,e))
    else: win_union.append((c,s,e))
base=sorted([(c,max(s,lo),min(e,hi)) for c,lo,hi in win_union for _,s,e in by[c] if s<hi and e>lo],key=lambda r:(rank[r[0]],r[1],r[2]))
selection='40 kb deep-pile + greedy 4 kb call-rich windows (target/max 32 records/windows per contig) + richest indel/multiallelic windows; ties lowest start; clipped source BED'
manifest={'schema':2,'slice':str(src),'source_bed':str(bed),'inputs':paths,
          'golden_vcf':str(golden),'golden_records':{c:len(v) for c,v in calls.items()},
          'selection':selection,'windows':win_union,
          'peaks':{c:peaks[c] for c in contigs},'pile_window_radius':20000,'call_window_radius':2000,'cases':[]}
def add(name,description,rr):
    description += '; selection: '+selection
    d=out/name; d.mkdir(exist_ok=True)
    data=''.join(f'{c}\t{s}\t{e}\n' for c,s,e in rr)
    (d/'intervals.bed').write_text(data); (d/'description.txt').write_text(description+'\n')
    manifest['cases'].append({'name':name,'description':description,'sha256':hashlib.sha256(data.encode()).hexdigest(),'records':len(rr),'bases_with_multiplicity':sum(e-s for c,s,e in rr)})
add('baseline','Call-rich and deep-pile baseline; GATK should emit supported calls',base)
add('reversed','Reverse all selected source BED records',list(reversed(base)))
shuffled=base.copy(); random.Random(4620).shuffle(shuffled)
add('shuffled','Shuffle all selected source BED records with seed 4620',shuffled)
for overlap in (0,1,50):
    rr=[]
    for c,s,e in base:
        if e-s<2: rr.append((c,s,e)); continue
        m=(s+e)//2; k=min(overlap,e-s)
        left_end=min(e,m+(k+1)//2)
        rr.extend([(c,s,left_end),(c,left_end-k,e)])
    add(f'split-{overlap}',f'Split every selected record at midpoint; overlap {overlap} bp (clamped on short records)',rr)
# Operate on the union to ensure the requested gaps really are uncovered, even where
# the original BED has overlapping transcript annotations.
union=[]
for c,s,e in base:
    if union and union[-1][0]==c and s<=union[-1][2]:
        a,b,z=union[-1]; union[-1]=(a,b,max(z,e))
    else: union.append((c,s,e))
rr=[]
for c,s,e in union:
    if rr and rr[-1][0]==c and s-rr[-1][2]<=300:
        a,b,z=rr[-1]; rr[-1]=(a,b,e)
    else: rr.append((c,s,e))
add('merge-short-gaps','Merge neighbouring union records across gaps <=300 bp',rr)
for gap in (1,50,100,300):
    rr=[]
    for c in contigs:
        pos=peaks[c]['start0']; s=max(0,pos-600); e=min(lengths[c],pos+600)
        rr.extend([(c,s,pos),(c,pos+gap,e)])
    add(f'gap-{gap}',f'Two 600 bp peak-centred flanks per contig separated by exactly {gap} uncovered bases',rr)
deep=[]
for c in contigs:
    pos=peaks[c]['start0']
    if pos<1 or pos+601>lengths[c]: sys.exit('peak too close to contig edge to construct deep cuts')
    deep.extend([(c,pos-600 if pos>=600 else 0,pos),(c,pos+1,min(lengths[c],pos+601))])
add('deep-cuts','Interval ends just before deepest shared alignment start; next starts one base after it (1 bp gap)',deep)
add('contig-edges','First two and last two bases of each selected contig',[(c,s,e) for c in contigs for s,e in [(0,2),(lengths[c]-2,lengths[c])]])
for n in (1,2):
    add(f'tiny-{n}',f'{n}-base interval starting at a covered alignment-start peak on each contig',[(c,peaks[c]['start0'],peaks[c]['start0']+n) for c in contigs])
add('duplicates','Every selected BED record repeated twice consecutively',[r for r in base for _ in range(2)])
if len(contigs)>1:
    add('contig-change','Deep cuts on both sides of direct contig transitions in dictionary order',deep)
    add('contig-nondictionary','Same deep cuts with contig blocks in reverse dictionary order',[r for c in reversed(contigs) for r in deep if r[0]==c])
# One deterministic call per contig: prefer indels, then multiallelics, then POS.
anchors={c:min(v,key=lambda v:(not v['indel'],not v['multi'],v['pos0'])) for c,v in calls.items() if v}
manifest['call_anchors']=anchors
for offset in (-1,0,1):
    add(f'call-edge-{offset:+d}',f'Rightmost included base is call POS {offset:+d}; GATK may change/drop the edge call through activity, assembly or output clipping',
        [(c,max(0,v['pos0']-600),v['pos0']+offset+1) for c,v in anchors.items()])
for n in (1,2):
    add(f'call-tiny-{n}',f'Call POS in a {n}-base interval; GATK may emit an overlapping call using padded reads or drop it if activity/assembly is insufficient',
        [(c,v['pos0'],v['pos0']+n) for c,v in anchors.items()])
indels={c:next((v for v in calls[c] if v['end0']-v['pos0']>1),None) for c in contigs}
indels={c:v for c,v in indels.items() if v is not None}
if indels:
    add('call-indel-span','Right edge inside deletion REF span, after anchor and before final REF base; GATK may emit the full overlapping deletion or change/drop it after reassembly',
        [(c,max(0,v['pos0']-600),v['pos0']+1) for c,v in indels.items()])
# Mutect2 default assembly-region-padding is 100 bp.
for side in ('inside','outside'):
    for distance in (99,101):
        rr=[]
        for c,v in anchors.items():
            edge=v['pos0']-distance if side=='inside' else v['pos0']+distance
            rr.append((c,max(0,edge),min(lengths[c],edge+600)))
        add(f'call-padding-{side}-{distance}',f'Call {distance} bp {side} the left edge (assembly padding 100 bp); GATK uses padded reads but outside calls are not output targets; inside calls may change with assembly context',rr)
add('call-padding-overlap','Two intervals with a 1 bp gap at call POS; 100 bp padded regions overlap across the call; anchor is outside target union so GATK may omit it or emit an overlapping event',
    [(c,s,e) for c,v in anchors.items() for s,e in [(max(0,v['pos0']-600),v['pos0']),(v['pos0']+1,min(lengths[c],v['pos0']+601))]])
# Explicit malformed controls exercise the rejection branch, independently of whether
# GATK accepts reordered, overlapping, or duplicate records.
c=contigs[0]
add('invalid-empty','Malformed BED: empty interval (end equals start)',[(c,10,10)])
add('invalid-beyond-contig','Malformed BED: end exceeds sequence dictionary length',[(c,lengths[c]-1,lengths[c]+1)])
(out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
for case in manifest['cases']: print(f"{case['name']}: {case['description']}")
print('Peaks (0-based): '+json.dumps(manifest['peaks'],sort_keys=True))
PY
