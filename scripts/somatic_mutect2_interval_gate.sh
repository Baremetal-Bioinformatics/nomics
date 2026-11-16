#!/usr/bin/env bash
# Usage: somatic_mutect2_interval_gate.sh <cases-dir>
# NOMICS, JAVA (an executable, including the acceptance docker wrapper), GATK_JAR,
# NOMICS_M2_REFERENCE as in the oracle scripts. M2_INTERVAL_CASES is an optional
# comma-separated subset. Reference results are reused only with matching provenance.
set -euo pipefail
repo=$(cd "$(dirname "$0")/.." && pwd)
export NOMICS=${NOMICS:-$repo/target/release/nomics}
export JAVA=${JAVA:-/opt/homebrew/opt/openjdk@17/bin/java}
export GATK_JAR=${GATK_JAR:-$HOME/omics/oracle/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar}
export NOMICS_M2_REFERENCE=${NOMICS_M2_REFERENCE:-$HOME/omics/ref/GATK.GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta}
python3 -I -u - "${1:?cases directory required}" "$repo" <<'PY'
import gzip, hashlib, itertools, json, os, pathlib, platform, re, shutil, subprocess, sys
out,repo=map(lambda s:pathlib.Path(s).resolve(),sys.argv[1:]); m=json.loads((out/'manifest.json').read_text())
for key in ('NOMICS','JAVA'):
    value=os.environ[key]
    os.environ[key]=shutil.which(value) or str(pathlib.Path(value).resolve())
for key in ('GATK_JAR','NOMICS_M2_REFERENCE'):
    os.environ[key]=str(pathlib.Path(os.environ[key]).resolve())
os.chdir(out)
tmp=pathlib.Path(os.environ.get('TMPDIR',str(out/'tmp'))).resolve()
# macOS exports a system TMPDIR by default; never write there. Explicit in-tree TMPDIR
# is honoured, otherwise use the output-local tmp directory.
if not tmp.is_relative_to(out): tmp=out/'tmp'
tmp.mkdir(parents=True,exist_ok=True)
os.environ['TMPDIR']=str(tmp); os.environ['TMP']=str(tmp); os.environ['TEMP']=str(tmp)
ref=os.environ['NOMICS_M2_REFERENCE']; inputs=m['inputs']
nomics=os.environ['NOMICS']; java=os.environ['JAVA']; jar=os.environ['GATK_JAR']
for p in [ref,ref+'.fai',str(pathlib.Path(ref).with_suffix('.dict')),*inputs.values(),
          inputs['normal']+'.bai',inputs['tumour']+'.bai',inputs['pon']+'.tbi',inputs['germline']+'.tbi',nomics]:
    if not pathlib.Path(p).is_file(): sys.exit(f'missing prerequisite {p} (no indexes or binaries are built by this gate)')
if platform.system()=='Darwin': lock=['lockf','-k','/tmp/nomics-heavy.lock']
else: lock=['flock',str(out/'heavy.lock')]
def run(cmd,d):
    d.mkdir(parents=True,exist_ok=True)
    with (d/'run.log').open('wb') as log:
        try: r=subprocess.run(lock+cmd,stdout=log,stderr=subprocess.STDOUT).returncode
        except OSError as e: log.write(str(e).encode()); r=127
    (d/'exit.txt').write_text(str(r)+'\n'); return r
version=out/'gatk-version'
vr=run([java,'-XX:-UsePerfData',f'-Djava.io.tmpdir={tmp}','-jar',jar,'--version'],version)
if vr or not re.search(r'4\.6\.2\.0', (version/'run.log').read_text(errors='replace')):
    sys.exit(f'GATK 4.6.2.0 qualification failed: {version}/run.log')
common=['--input',inputs['normal'],'--input',inputs['tumour'],'--normal-sample','HCC1395_HCC1395BL_WES',
        '--reference',ref,'--panel-of-normals',inputs['pon'],'--germline-resource',inputs['germline']]
def oracle(case,d):
    return run([java,'-Xmx4g','-XX:-UsePerfData',f'-Djava.io.tmpdir={tmp}','-jar',jar,'Mutect2',*common,
                '--intervals',str(case/'intervals.bed'),'--tmp-dir',str(tmp),
                '--output',str(d/'calls.vcf.gz'),'--f1r2-tar-gz',str(d/'f1r2.tar.gz')],d)
def compare(a,b,log):
    ok=True
    with log.open('w') as f:
        try:
            def vcf(d):
                with gzip.open(d/'calls.vcf.gz','rb') as v: lines=v.readlines()
                return ([l for l in lines if l.startswith(b'#') and not l.startswith(b'##GATKCommandLine=')],
                        [l for l in lines if not l.startswith(b'#')])
            av,bv=vcf(a),vcf(b)
            for label,x,y in [('header',av[0],bv[0]),('record',av[1],bv[1])]:
                if x!=y:
                    ok=False
                    for i,(u,v) in enumerate(itertools.zip_longest(x,y),1):
                        if u!=v: f.write(f'first differing {label} {i}: reference={u!r}; other={v!r}\n'); break
            x=(a/'calls.vcf.gz.stats').read_bytes(); y=(b/'calls.vcf.gz.stats').read_bytes()
            if x!=y:
                ok=False; i=next(i for i,(u,v) in enumerate(itertools.zip_longest(x,y)) if u!=v)
                f.write(f'first differing stats byte {i} (0-based): reference={x[i:i+32]!r}; other={y[i:i+32]!r}\n')
            r=subprocess.run([sys.executable,'-I',str(repo/'scripts/somatic_f1r2_cmp.py'),str(a/'f1r2.tar.gz'),str(b/'f1r2.tar.gz'),'f1r2'],stdout=f,stderr=f)
            ok=ok and r.returncode==0
        except (OSError,EOFError,ValueError) as e: f.write(f'incomplete/invalid output: {e}\n'); ok=False
    return ok, len(av[1]) if 'av' in locals() else 0, len(bv[1]) if 'bv' in locals() else 0
def signature(case):
    # Input metadata plus BED bytes; changing any source or JVM invalidates the cache.
    files=[ref,ref+'.fai',str(pathlib.Path(ref).with_suffix('.dict')),*inputs.values(),
           inputs['normal']+'.bai',inputs['tumour']+'.bai',inputs['pon']+'.tbi',inputs['germline']+'.tbi',jar,java]
    meta=[]
    for p in files:
        q=pathlib.Path(p)
        meta.append((p,q.stat().st_size,q.stat().st_mtime_ns) if q.exists() else (p,'wrapper/container path'))
    return hashlib.sha256(json.dumps([common,meta,str(tmp),hashlib.sha256((repo/'scripts/somatic_mutect2_interval_gate.sh').read_bytes()).hexdigest(),hashlib.sha256((case/'intervals.bed').read_bytes()).hexdigest()],sort_keys=True).encode()).hexdigest()
def rejected(d,rc):
    # A crashed JVM, missing input, or failed lock is not a reference rejection.
    text=(d/'run.log').read_text(errors='replace')
    error=re.search(r'A USER ERROR has occurred:([^\n]*(?:\n[^\n]+)?)',text)
    return rc==2 and error is not None and bool(re.search(r'interval|bed|contig|coordinate|genome loc|feature (start|end)',error[1],re.I))
selected=set(filter(None,os.environ.get('M2_INTERVAL_CASES','').split(',')))
names={c['name'] for c in m['cases']}
if selected-names: sys.exit(f'unknown cases: {sorted(selected-names)}')
# Four nonempty valid cases prevent a baseline-only records gate. Subsets obey
# the same rule and must include baseline; empty cases still compare fully.
minimum=4
fail=0; count=0; rejects=0; nonempty=0; total=0; empty=0; baseline_records=0
for entry in m['cases']:
    name=entry['name']
    if selected and name not in selected: continue
    count+=1; case=out/name; a=case/'gatk'; b=case/'nomics'; sig=signature(case)
    if hashlib.sha256((case/'intervals.bed').read_bytes()).hexdigest()!=entry['sha256']:
        print(f'FAIL {name}: BED differs from manifest'); fail+=1; continue
    cache=a/'signature.txt'
    if cache.exists() and cache.read_text()==sig and (a/'exit.txt').exists(): ar=int((a/'exit.txt').read_text())
    else:
        ar=oracle(case,a)
        # Do not cache infrastructure failures.
        if ar==0 or rejected(a,ar): cache.write_text(sig)
        elif cache.exists(): cache.unlink()
    br=run([nomics,'somatic','mutect2',*common,'--intervals',str(case/'intervals.bed'),
            '--output',str(b/'calls.vcf.gz'),'--f1r2-tar-gz',str(b/'f1r2.tar.gz')],b)
    good=False; records=other_records=0
    if ar==0 and br==0:
        good,records,other_records=compare(a,b,case/'compare.log')
        total+=records
        if not name.startswith('invalid-') and records and other_records: nonempty+=1
        if name=='baseline': baseline_records=records
    record_note=f' records compared={records} (GATK={records}, nomics={other_records})' if ar==br==0 else ' records compared=0'
    if good and records==0: empty+=1
    if rejected(a,ar) and br in (1,2):
        good=True; rejects+=1
    if good:
        print(f"PASS{' (empty)' if ar==0 and records==0 else ''} {name}: GATK={ar} nomics={br}"+(' (both reject)' if ar else ' (records/header/stats/f1r2)')+record_note)
    else:
        fail+=1
        # Every failure gets a second independent GATK invocation, including failures
        # involving rejection; never silently replace the first oracle with the repeat.
        repeat=case/'gatk-repeat'; rr=oracle(case,repeat)
        if ar==rr==0: stability='identical' if compare(a,repeat,case/'reference-repeat.log')[0] else 'DIFFERENT'
        elif rejected(a,ar) and rejected(repeat,rr): stability='both reject'
        else: stability='failed/inconsistent'
        print(f'FAIL {name}: GATK={ar} nomics={br}; GATK repeat={rr} ({stability}); logs={case}'+record_note)
        if (case/'compare.log').exists():
            print((case/'compare.log').read_text().strip())
if baseline_records==0:
    print('FAIL non-empty rule: baseline must be selected and contain records'); fail+=1
if nonempty<minimum:
    print(f'FAIL non-empty rule: {nonempty} valid cases with records; require at least {minimum}'); fail+=1
print(f'interval gate: {count} cases, {nonempty} cases with records, {total} total records compared, {empty} empty passes, {rejects} paired rejections, {fail} failures; outputs: {out}')
sys.exit(bool(fail) or count==0)
PY
