#!/usr/bin/env python3
"""Small, optional fixture differential test; no third-party packages."""
import gzip
import os
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile

binary = sys.argv[1]
home = Path.home()
fixture = Path(os.environ.get('NOMICS_CHAIN_FIXTURE', home / 'omics/runs/somatic/m2-multi'))
ref = Path(os.environ.get('NOMICS_CHAIN_REFERENCE', home / 'omics/ref/GATK.GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta'))
if not fixture.is_dir() or not ref.is_file():
    print('SKIP chain fixture: m2-multi or reference absent')
    sys.exit(0)

def run(*args):
    subprocess.run([binary, 'somatic', *map(str, args)], check=True)

def content(p):
    if str(p).endswith('.tar.gz'):
        with tarfile.open(p) as t:
            return [(m.name, t.extractfile(m).read()) for m in t if m.isfile()]
    if str(p).endswith('.vcf.gz'):
        with gzip.open(p) as f:
            return [l for l in f if not l.startswith(b'##GATKCommandLine=')]
    return p.read_bytes()

with tempfile.TemporaryDirectory(prefix='nomics-chain-fixture-') as tmp:
    tmp = Path(tmp)
    # Three tiny slices centered on actual oracle calls, one per chromosome.
    oracle = next((fixture / 'oracle').glob('*.vcf.gz'))
    loci = {}
    with gzip.open(oracle, 'rt') as f:
        for line in f:
            if not line.startswith('#'):
                fields = line.split('\t')
                loci.setdefault(fields[0], int(fields[1]))
    assert len(loci) >= 2
    beds = []
    for i, (chrom, pos) in enumerate(loci.items()):
        bed = tmp / f'g{i}.bed'
        bed.write_text(f'{chrom}\t{max(0, pos-500)}\t{pos+500}\n')
        beds.append(bed)
    tumour = fixture / 'HCC1395_WES.multi.bam'
    normal = fixture / 'HCC1395BL_WES.multi.bam'
    pon = fixture / '1000g_pon.multi.vcf.gz'
    germ = fixture / 'af-only-gnomad.multi.vcf.gz'
    common = ['--reference', ref, '--panel-of-normals', pon, '--germline-resource', germ]
    args = ['--tumour', tumour, '--normal', normal, *common, '--pileup-resource', germ,
            '--tumour-label', 'T', '--normal-label', 'N', '--outdir', tmp/'chain', '--keep-intermediates']
    for bed in beds:
        args += ['--intervals', bed]
    order = list(reversed(beds))
    for bed in order:
        args += ['--orientation-order', bed.stem]
    run('chain', *args)
    assert any(not row.startswith(b'#') for row in content(tmp/'chain/variant_calling/mutect2/T_vs_N/T_vs_N.mutect2.vcf.gz')), 'fixture must exercise nonempty calls'
    separate = tmp/'separate'
    separate.mkdir()
    p = 'T_vs_N.mutect2'
    vcfs, stats, archives, tp, np = [], [], {}, [], []
    for bed in beds:
        v = separate/f'{p}.{bed.stem}.vcf.gz'
        ar = separate/f'{p}.{bed.stem}.f1r2.tar.gz'
        run('mutect2', '--input', normal, '--input', tumour, '--normal-sample', 'HCC1395_HCC1395BL_WES', *common, '--intervals', bed, '--output', v, '--f1r2-tar-gz', ar)
        vcfs.append(v); stats.append(Path(str(v)+'.stats')); archives[bed.stem] = ar
        for bam, label, tables in [(tumour,'T',tp),(normal,'N',np)]:
            table = separate/f'{label}.mutect2.{bed.stem}.pileups.table'
            run('get-pileup-summaries','--input',bam,'--reference',ref,'--variant',germ,'--intervals',bed,'--output',table)
            tables.append(table)
    def repeated(flag, values):
        return [v for p in values for v in (flag,p)]
    raw = separate/f'{p}.vcf.gz'
    stat = Path(str(raw)+'.stats')
    priors = separate/f'{p}.artifactprior.tar.gz'
    run('merge-vcfs', *repeated('--input',vcfs), '--sequence-dictionary',ref.with_suffix('.dict'),'--output',raw)
    run('merge-mutect-stats', *repeated('--stats',stats),'--output',stat)
    run('learn-read-orientation-model',*repeated('--input',[archives[b.stem] for b in order]),'--output',priors)
    for label, tables in [('T',tp),('N',np)]:
        run('gather-pileup-summaries',*repeated('--input',tables),'--sequence-dictionary',ref.with_suffix('.dict'),'--output',separate/f'{label}.mutect2.pileups.table')
    contam = separate/f'{p}.contamination.table'; seg = separate/f'{p}.segmentation.table'
    run('calculate-contamination','--input',separate/'T.mutect2.pileups.table','--matched-normal',separate/'N.mutect2.pileups.table','--output',contam,'--tumor-segmentation',seg)
    run('filter-mutect-calls','--variant',raw,'--stats',stat,'--reference',ref,'--orientation-bias-artifact-priors',priors,'--contamination-table',contam,'--tumor-segmentation',seg,'--output',separate/f'{p}.filtered.vcf.gz')
    for actual in (tmp/'chain').rglob('*'):
        # Timing metadata has no counterpart among the standalone tools.
        if actual.is_file() and actual != tmp/'chain/chain-timings.tsv':
            assert content(actual) == content(separate/actual.name), actual
    # Reverse all four reduction inputs: table/VCF/stats invariance in the gated regime;
    # orientation is compared separately and is permitted to differ, never normalized.
    for tool, flag, files, output in [('merge-vcfs','--input',vcfs,raw),('merge-mutect-stats','--stats',stats,stat),('gather-pileup-summaries','--input',tp,separate/'T.mutect2.pileups.table')]:
        before = content(output)
        extra = [] if tool == 'merge-mutect-stats' else ['--sequence-dictionary',ref.with_suffix('.dict')]
        run(tool,*repeated(flag,reversed(files)),*extra,'--output',output)
        assert before == content(output), tool
    before = content(priors)
    run('learn-read-orientation-model',*repeated('--input',[archives[b.stem] for b in beds]),'--output',priors)
    print('orientation reversed-order member bytes differ:', before != content(priors))
    print('PASS chain: every published and intermediate output matches separate subcommands; three merge-order checks passed')
    # The two reference orientation archives are small and independent of full BAMs.
    golden = home/'omics/runs/somatic/S1-r1/golden'
    x = 'HCC1395_WES_vs_HCC1395BL_WES'
    directories = [golden/f'MUTECT2_PAIRED-{x}', golden/f'MUTECT2_PAIRED-{x}.2']
    inputs = [list(d.glob('*.f1r2.tar.gz')) for d in directories]
    if all(len(paths) == 1 for paths in inputs):
        inputs = [paths[0] for paths in inputs]
        run('learn-read-orientation-model',*repeated('--input',inputs),'--output',priors)
        forward = content(priors)
        run('learn-read-orientation-model',*repeated('--input',reversed(inputs)),'--output',priors)
        backward = content(priors)
        assert forward != backward, 'reference orientation order must remain observable'
        changed = sum(a != b for (_, left), (_, right) in zip(forward, backward)
                      for lrow, rrow in zip(left.splitlines(), right.splitlines())
                      for a, b in zip(lrow.split(b'\t'), rrow.split(b'\t')))
        print(f'PASS reference orientation reversal changes {changed} table fields')
    else:
        print('SKIP reference orientation reversal: S1-r1 archives absent')
