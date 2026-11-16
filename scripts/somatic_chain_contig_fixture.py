#!/usr/bin/env python3
"""Build a seven-contig hand-off stress fixture from local chr slices.
Run under lockf -k /tmp/nomics-heavy.lock. Uses reference samtools/bgzip/tabix
only for fixture construction, never as runtime dependencies. No input writes.
The oracle for the resulting merged/sampled input is nomics commit 882d6f0.
"""
import argparse
from collections import Counter
import gzip
from pathlib import Path
import subprocess


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('out', type=Path)
    args = p.parse_args()
    args.out.mkdir()
    base = Path.home() / 'omics/runs/somatic'
    slices = [(i, base / f'm2-chr{i}') for i in range(12, 23) if (base / f'm2-chr{i}').is_dir()]
    assert len(slices) >= 4
    fai = Path.home() / 'omics/ref/GATK.GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta.fai'
    lengths = {f[0]: int(f[1]) for line in fai.read_text().splitlines() if (f := line.split())}
    intervals = []
    for i, directory in slices:
        chrom = f'chr{i}'
        # Actual call sites plus deep same-start pools exercise both outputs and
        # RNG hand-off. Keep whole input BAMs; restrict the traversal by BED.
        oracle = next((directory / 'oracle').glob('*.vcf.gz'))
        positions = []
        with gzip.open(oracle, 'rt') as vcf:
            for line in vcf:
                if not line.startswith('#'):
                    positions.append(int(line.split('\t')[1]))
        positions = positions[:30]
        counts = Counter()
        for sample in ['HCC1395BL', 'HCC1395']:
            cmd = ['samtools', 'view', '-q', '51', '-F', str(0xF04), directory / f'{sample}_WES.{chrom}.bam']
            with subprocess.Popen(list(map(str, cmd)), stdout=subprocess.PIPE, text=True) as proc:
                for line in proc.stdout:
                    fields = line.split('\t', 5)
                    counts[int(fields[3])] += 1
                assert proc.wait() == 0
        deepest = counts.most_common(30)
        positions.extend(pos for pos, _ in deepest)
        print(chrom, 'deepest_joint_raw_pools', deepest, flush=True)
        for pos in sorted(set(positions)):
            intervals.append(f'{chrom}\t{max(0, pos-1000)}\t{min(lengths[chrom], pos+1000)}\n')
    (args.out / 'multi.bed').write_text(''.join(intervals))
    for sample in ['HCC1395BL', 'HCC1395']:
        output = args.out / f'{sample}_WES.multi.bam'
        inputs = [directory / f'{sample}_WES.chr{i}.bam' for i, directory in slices]
        # merge preserves coordinate order and combines identical RG/PG IDs.
        subprocess.run(list(map(str, ['samtools', 'merge', '--no-PG', '-c', '-p', '-o', output, *inputs])), check=True)
        subprocess.run(['samtools', 'index', str(output)], check=True)
    for resource in ['1000g_pon', 'af-only-gnomad']:
        plain = args.out / f'{resource}.multi.vcf'
        with plain.open('wb') as out:
            for index, (i, directory) in enumerate(slices):
                with gzip.open(directory / f'{resource}.chr{i}.vcf.gz', 'rb') as src:
                    for line in src:
                        if index == 0 or not line.startswith(b'#'):
                            out.write(line)
        compressed = plain.with_suffix('.vcf.gz')
        with compressed.open('wb') as out:
            subprocess.run(['bgzip', '-c', str(plain)], stdout=out, check=True)
        subprocess.run(['tabix', '-p', 'vcf', str(compressed)], check=True)
    print(f'FIXTURE contigs={len(slices)} windows={len(intervals)} oracle=882d6f0', flush=True)


if __name__ == '__main__':
    main()
