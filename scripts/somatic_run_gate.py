#!/usr/bin/env python3
"""Driver for the shared full_exome preprocessing gates and CHAIN_PRODUCED."""
import csv
import os
from pathlib import Path
import re
import shlex
import subprocess as sp
import sys

repo = Path(__file__).resolve().parent.parent
run, out = map(lambda p: Path(p).resolve(), sys.argv[1:])
bin = str(Path(os.environ['NOMICS']).resolve())
ref = Path(os.environ.get('REF', str(Path.home() / 'omics/ref/GATK.GRCh38')))
rows = list(csv.DictReader((run / 'samplesheet.csv').open()))
# The Rust parser remains authoritative; reject ambiguity here before using labels in paths.
assert rows and len({r['patient'] for r in rows}) == 1
samples = {}
for r in rows:
    assert r['status'] in ('0', '1') and re.fullmatch(r'[A-Za-z0-9_][A-Za-z0-9_.-]*', r['sample'])
    samples.setdefault(r['status'], set()).add(r['sample'])
assert set(samples) == {'0', '1'} and all(len(s) == 1 for s in samples.values())
normal, tumour = (next(iter(samples[s])) for s in ('0', '1'))
params = dict(line.split(':', 1) for line in (run / 'params.yaml').read_text().splitlines() if ':' in line and not line.startswith('#'))
bed = Path(os.environ.get('TARGETS', params['intervals'].strip().strip('\'"')))
assert bed.is_file(), bed
assert not out.exists(), f'refusing existing output {out}'
# Everything written by this gate is under OUT_DIR, except TMPDIR-owned library temporaries.
out.mkdir()
orders = []
for sample in (normal, tumour):
    order = out / f'{sample}.scatter-order.txt'
    with order.open('w') as f:
        sp.run([sys.executable, str(repo / 'scripts/scatter_order.py'), str(run), sample], stdout=f, check=True)
    orders.append(order)
pair = f'{tumour}_vs_{normal}'
prefix = f'{pair}.mutect2'
command = (run / f'golden/LEARNREADORIENTATIONMODEL-{pair}/.command.sh').read_text()
orientation = re.findall(r'--input\s+(\S+)', command)
orientation = [Path(p).name.removeprefix(prefix + '.').removesuffix('.f1r2.tar.gz') for p in orientation]
assert orientation and set(orientation) == {p.stem for p in (run / 'golden/intervals').glob('*.bed')}
# The run settings of scripts/full_exome.sh, so a pair is preprocessed as a gated single sample is:
# the same known sites, dense SA cache and interval, thread count and GPU seeding.
bundle = ref / 'Annotation/GATKBundle'
sa = os.environ.get('NOMICS_SA_INTERVAL', '2')
args = [bin, 'run', '--samplesheet', str(run / 'samplesheet.csv'), '--outdir', str(out / 'pair'),
        '--reference', str(ref), '--targets', str(bed), '--write-intermediates', '--wes',
        '--known-sites', str(bundle / 'dbsnp_146.hg38.vcf.gz'),
        '--known-sites', str(bundle / 'beta/Homo_sapiens_assembly38.known_indels.vcf.gz'),
        '--known-sites', str(bundle / 'Mills_and_1000G_gold_standard.indels.hg38.vcf.gz'),
        '--dense-sa-cache', os.environ.get('DENSE_SA_CACHE', str(Path.home() / f'omics/cache/sa{sa}.GRCh38.nsa')),
        '--dense-sa-interval', sa,
        '--known-sites-cache', str(Path.home() / 'omics/cache/known_sites.GATK.GRCh38.nks'),
        '--threads', os.environ.get('NOMICS_THREADS', '20'),
        '--read-group-ds', params['fasta'].strip().strip('\'"'),
        '--normal-scatter-order', str(orders[0]), '--tumour-scatter-order', str(orders[1]),
        '--orientation-order', ','.join(orientation)]
for flag, name in [('germline-resource', 'af-only-gnomad.hg38.vcf.gz'), ('pileup-resource', 'af-only-gnomad.hg38.vcf.gz'), ('panel-of-normals', '1000g_pon.hg38.vcf.gz')]:
    args += ['--' + flag, str(ref / 'Annotation/GATKBundle' / name)]
# In memory, duplicate marking needs about 7x the gzipped FASTQ plus the index, so a large exome
# is streamed, at the threshold scripts/goldens/exome_golden.sh uses for single samples.
fq_gz = {}
for r in rows:
    fq_gz[r['sample']] = fq_gz.get(r['sample'], 0) + sum(Path(os.path.expandvars(r[k])).stat().st_size for k in ('fastq_1', 'fastq_2'))
if max(fq_gz.values()) > 8_000_000_000:
    args.append('--spill-prep')
args += shlex.split(os.environ.get('NOMICS_RUN_EXTRA', ''))
run_env = dict(os.environ)
if not os.environ.get('NOMICS_CPU_ONLY'):
    run_env['NOMICS_GPU_SEED'] = '1'
(out / 'command.txt').write_text(shlex.join(args) + '\n')
with (out / 'run.log').open('w') as log:
    sp.run(['/usr/bin/time', '-v', '-o', str(out / 'time.txt'), *args], stdout=log, stderr=sp.STDOUT, check=True, env=run_env)
for sample in (normal, tumour):
    print(f'== {sample} stages', flush=True)
    print((out / f'pair/timings/{sample}.tsv').read_text(), flush=True)
print((out / 'pair/stages.tsv').read_text(), flush=True)
print((out / 'time.txt').read_text(), flush=True)
# Exact scatter bytes, not just group names.
ours = {p.name: p.read_bytes() for p in (out / 'pair/intervals').glob('*.bed')}
golden = {p.name: p.read_bytes() for p in (run / 'golden/intervals').glob('*.bed')}
assert ours == golden, 'scatter group names/contents differ'
for sample in (normal, tumour):
    view = out / f'gate-{sample}'
    view.mkdir()
    pre = out / 'pair/results/preprocessing'
    (view / 'prep').symlink_to(out / f'pair/work/preprocessing/{sample}/prep')
    (view / 'mapped').symlink_to(pre / 'mapped' / sample)
    for stage, suffix in [('markduplicates', 'md.bam'), ('recalibrated', 'recal.bam'), ('recal_table', 'recal.table')]:
        name = f'{sample}.{suffix}'
        (view / name).symlink_to(pre / stage / sample / name)
    (view / f'{sample}.md.metrics').symlink_to(out / f'pair/results/reports/markduplicates/{sample}/{sample}.md.bam.metrics')
    tasks = []
    for p in (run / 'work').glob('*/*/.command.sh'):
        if 'fastp' in p.read_text() and list(p.parent.glob(f'*.{sample}-*_R1.fastp.fastq.gz')) and not (p.parent / 'BWAIndex').exists():
            tasks.append(p)
    lanes = sum(r['sample'] == sample for r in rows)
    assert len(tasks) == lanes, f'{sample}: expected {lanes} fastp tasks, got {len(tasks)}'
    chunk_names = set()
    for task in tasks:
        chunks = {p.name for p in task.parent.glob('*.fastp.fastq.gz')}
        assert len(chunks) == 24, f'{task}: missing fastp chunks'
        chunk_names.update(chunks)
    assert len(chunk_names) == 24 * lanes
    assert chunk_names == {p.name for p in (view / 'prep').glob('*.fastp.fastq.gz')}, 'missing/extra produced prep chunks'
    assert len(list((view / 'mapped').glob('*.bam'))) == 12 * lanes, 'missing mapped chunks'
    metrics = list((run / 'results/reports').rglob(f'{sample}.md.bam.metrics'))
    assert len(metrics) == 1, 'missing or ambiguous reference metrics'
    env = dict(os.environ, REPO=str(repo), OUT=str(view), GOLD=str(run), P=sample,
               G=str(view / 'gates.txt'), BIN=str(Path(bin).parent), GOLD_METRICS=str(metrics[0]),
               PREP_TASK_COMMANDS='\n'.join(map(str, tasks)), NOMICS=bin)
    # BIN/nomics is the legacy helper spelling; use a local link to honor arbitrary NOMICS names.
    bindir = view / 'bin'; bindir.mkdir(); (bindir / 'nomics').symlink_to(bin); env['BIN'] = str(bindir)
    sp.run(['bash', '-euo', 'pipefail', '-c', 'source "$REPO/scripts/preprocessing_gates.sh"; preprocessing_gates'], env=env, check=True)
    verdicts = (view / 'gates.txt').read_text().splitlines()
    print('\n'.join(verdicts), flush=True)
    assert len(verdicts) == lanes + 12 * lanes + 4 and all(v.startswith('PASS ') for v in verdicts), 'missing/failed preprocessing comparison'
env = dict(os.environ, NOMICS=bin, REF=str(ref), TUMOUR_SAMPLE=tumour, NORMAL_SAMPLE=normal,
           NORMAL_SM=f"{rows[0]['patient']}_{normal}", CHAIN_PRODUCED=str(out / 'pair/results'),
           CHAIN_WORK=str(out / 'pair/work/chain'), CHAIN_LOG=str(out / 'run.log'))
sp.run(['bash', str(repo / 'scripts/somatic_chain_gate.sh'), str(run), str(out / 'chain-gate')], env=env, check=True)
print('PASS somatic run: both preprocessing paths, scatter bytes and every chain comparison')
