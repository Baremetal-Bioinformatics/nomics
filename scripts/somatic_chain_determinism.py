#!/usr/bin/env python3
"""Exact chain regression: baseline plus six runs, including compressed bytes.
Run under lockf -k /tmp/nomics-heavy.lock on macOS. Fresh output root required.
The sole excluded file is the new, deliberately variable chain-timings.tsv.
"""
import argparse
import filecmp
import gzip
import hashlib
import os
import re
from pathlib import Path
import subprocess


def manifest(root):
    return {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(root.rglob('*')) if p.is_file() and p.relative_to(root) != Path('chain-timings.tsv')}


def records(root):
    vcfs = list((root / 'variant_calling').rglob('*.vcf.gz'))
    return sum(sum(not line.startswith(b'#') for line in gzip.open(p, 'rb')) for p in vcfs)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('baseline', type=Path)
    p.add_argument('binary', type=Path)
    p.add_argument('out', type=Path)
    p.add_argument('--threads', type=int, nargs='+', help='thread counts (default: 1, 4, all)')
    p.add_argument('--repeats', type=int, default=2)
    p.add_argument('--one-group', action='store_true')
    p.add_argument('--before', type=Path, help='5f2bc31 binary for before timings')
    p.add_argument('--self-check', action='store_true', help='additional serial-boundary/output check at four threads')
    p.add_argument('--baseline-only', action='store_true')
    p.add_argument('--resume', action='store_true', help='reuse completed baseline; candidate directories must be absent')
    p.add_argument('--serial', action='store_true', help='timings-only commit has no thread flag')
    args = p.parse_args()
    thread_counts = args.threads or [1, 4, os.cpu_count() or 1]
    if args.repeats < 1 or any(t < 1 for t in thread_counts):
        p.error('thread counts and repeats must be positive')
    fixture = Path(os.environ.get('NOMICS_CHAIN_FIXTURE', Path.home() / 'omics/runs/somatic/m2-multi'))
    ref = Path(os.environ.get('NOMICS_CHAIN_REFERENCE', Path.home() / 'omics/ref/GATK.GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta'))
    args.out.mkdir(exist_ok=args.resume)
    beds = args.out / 'beds'
    if not args.resume:
        beds.mkdir()
        # Two nonoverlapping groups; first spans two contigs (shared RNG stream).
        lines = (fixture / 'multi.bed').read_text().splitlines(True)
        contigs = list(dict.fromkeys(line.split()[0] for line in lines))
        assert len(contigs) >= 3
        for name, chroms in ([('g0', contigs)] if args.one_group else [('g0', contigs[:-1]), ('g1', contigs[-1:])]):
            (beds / f'{name}.bed').write_text(''.join(l for l in lines if l.split()[0] in chroms))
    common = ['somatic', 'chain', '--tumour', fixture / 'HCC1395_WES.multi.bam',
              '--normal', fixture / 'HCC1395BL_WES.multi.bam', '--reference', ref,
              '--panel-of-normals', fixture / '1000g_pon.multi.vcf.gz',
              '--germline-resource', fixture / 'af-only-gnomad.multi.vcf.gz',
              '--pileup-resource', fixture / 'af-only-gnomad.multi.vcf.gz',
              '--intervals-dir', beds, '--orientation-order', 'g0' if args.one_group else 'g1,g0',
              '--tumour-label', 'T', '--normal-label', 'N', '--keep-intermediates']

    def run(binary, name, threads=None, self_check=False):
        cmd = ['/usr/bin/time', '-l', binary.resolve(), *common, '--outdir', args.out / name]
        if threads is not None and not args.serial:
            cmd += ['--threads', str(threads)]
        print(f'RUN {name}', flush=True)
        env = dict(os.environ)
        env.pop('NOMICS_CONTIG_SELF_CHECK', None)
        if self_check:
            env['NOMICS_CONTIG_SELF_CHECK'] = '1'
        with (args.out / f'{name}.log').open('w') as log:
            subprocess.run(list(map(str, cmd)), stdout=log, stderr=log, check=True, env=env)

    baseline = args.out / 'baseline'
    if not args.resume:
        run(args.baseline, 'baseline')
    expected = manifest(baseline)
    count = records(baseline)
    assert expected and count > 0, 'zero files or records'
    print(f'BASELINE files={len(expected)} published_records={count}', flush=True)
    if args.baseline_only:
        return
    if args.before:
        for threads in thread_counts:
            run(args.before, f'before-t{threads}', threads)
            assert manifest(args.out / f'before-t{threads}') == expected
    comparisons = 0
    expected_calls = None
    for threads in ([1] if args.serial else thread_counts):
        for repeat in ([1] if args.serial else range(1, args.repeats + 1)):
            name = 'serial-timing' if args.serial else f't{threads}-r{repeat}'
            run(args.binary, name, threads)
            if not args.serial:
                calls = sorted(re.findall(r'^contig_rng\t([^\t]+)\t([^\t]+)\t(\d+)$',
                               (args.out / f'{name}.log').read_text(), re.M))
                assert calls and sum(int(n) for _, _, n in calls) > 0, 'no generator calls in multi-contig run'
                assert len({(g, c) for g, c, _ in calls}) == len(calls), 'duplicate contig plan'
                if expected_calls is None:
                    expected_calls = calls
                assert calls == expected_calls, 'generator call counts changed'
                print(f'RNG {name} ' + ','.join(f'{g}/{c}={n}' for g,c,n in calls), flush=True)
            actual = manifest(args.out / name)
            assert actual.keys() == expected.keys(), f'{name}: file set differs'
            for path, digest in expected.items():
                assert actual[path] == digest and filecmp.cmp(baseline / path, args.out / name / path, shallow=False), f'{name}: byte mismatch {path}'
            assert records(args.out / name) == count
            comparisons += len(expected)
            print(f'PASS {name} files={len(expected)} published_records={count} compressed=byte-exact', flush=True)
    if args.self_check:
        run(args.binary, 'self-check', 4, True)
        assert manifest(args.out / 'self-check') == expected
        print('PASS serial boundary and concatenated compressed output self-check', flush=True)
    runs = 1 if args.serial else len(thread_counts) * args.repeats
    print(f'PASS runs={runs} files_per_run={len(expected)} file_comparisons={comparisons} published_records_per_run={count} record_comparisons={runs * count}', flush=True)


if __name__ == '__main__':
    main()
