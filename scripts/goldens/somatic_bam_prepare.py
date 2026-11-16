"""Prepare BAM-started reference inputs; invoked by somatic_bam_golden.sh."""
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import sys


def prepare():
    e = os.environ
    art, data = Path(e['ART']), Path(e['DATA'])
    for key in ('TUMOUR_BAM', 'NORMAL_BAM', 'TUMOUR_BAI', 'NORMAL_BAI', 'DICT', 'INTERVALS', 'BUNDLE'):
        e[key] = str(Path(e[key]).resolve())
    dictionary = Path(e['DICT'])
    contigs = {}
    for line in dictionary.read_text().splitlines():
        if line.startswith('@SQ\t'):
            fields = dict(x.split(':', 1) for x in line.split('\t')[1:])
            contigs[fields['SN']] = (len(contigs), int(fields['LN']))
    if not contigs:
        raise ValueError(f'no @SQ records in {dictionary}')
    source = Path(e['INTERVALS'])
    picard = source.name.endswith('.interval_list')
    rows = []
    previous = None
    for number, line in enumerate(source.read_text().splitlines(), 1):
        if not line.strip() or line.startswith(('@', '#', 'track ', 'browser ')):
            continue
        fields = line.split()
        try:
            chrom, start, end = fields[0], int(fields[1]) - int(picard), int(fields[2])
            rank, length = contigs[chrom]
        except (IndexError, ValueError, KeyError) as exc:
            raise ValueError(f'{source}:{number}: invalid interval or unknown contig') from exc
        if not 0 <= start < end <= length:
            raise ValueError(f'{source}:{number}: interval outside contig bounds')
        if previous and (rank < previous[0] or (rank == previous[0] and start < previous[1])):
            raise ValueError(f'{source}:{number}: unsorted or overlapping intervals in dictionary order')
        previous = (rank, end)
        rows.append(f'{chrom}\t{start}\t{end}\n')
    if not rows:
        raise ValueError('empty intervals')
    bed = data / 'intervals.bed'
    bed.write_text(''.join(rows))
    inputs = [Path(e[k]) for k in ('TUMOUR_BAM', 'NORMAL_BAM', 'TUMOUR_BAI', 'NORMAL_BAI')]
    inputs += [source, bed, dictionary]
    bundle = Path(e['BUNDLE'])
    resources = {k: bundle / v for k, v in {
        'germline_resource': 'af-only-gnomad.hg38.vcf.gz',
        'germline_resource_tbi': 'af-only-gnomad.hg38.vcf.gz.tbi',
        'pon': '1000g_pon.hg38.vcf.gz', 'pon_tbi': '1000g_pon.hg38.vcf.gz.tbi'}.items()}
    inputs += list(resources.values())
    checksums = []
    for path in inputs:
        stat = path.stat()
        identity = [stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns, stat.st_ino]
        cache = Path(str(path) + '.md5.json')
        record = None
        if path.suffix == '.bam' and cache.exists():
            record = json.loads(cache.read_text())
        if not record or record.get('identity') != identity:
            h = hashlib.md5()
            with path.open('rb') as fh:
                for block in iter(lambda: fh.read(8 * 1024 * 1024), b''):
                    h.update(block)
            record = {'identity': identity, 'md5': h.hexdigest()}
            if path.suffix == '.bam':
                tmp = Path(str(cache) + '.tmp')
                tmp.write_text(json.dumps(record) + '\n')
                tmp.replace(cache)
        display = str(path).replace(str(Path.home() / 'omics') + '/', '')
        checksums.append(f"{record['md5']}  {display}\n")
    (art / 'inputs.md5').write_text(''.join(checksums))
    with (data / 'samplesheet.csv').open('w', newline='') as fh:
        writer = csv.writer(fh, lineterminator='\n')
        writer.writerow(['patient', 'sex', 'status', 'sample', 'bam', 'bai'])
        for role, status in [('NORMAL', 0), ('TUMOUR', 1)]:
            writer.writerow([e['PATIENT'], e['SEX'], status, e[role + '_SAMPLE'], e[role + '_BAM'], e[role + '_BAI']])
    # Keep the reference settings; omit alignment, BQSR and FASTQ-only settings.
    omit = {'bwa', 'dbsnp', 'dbsnp_tbi', 'known_indels', 'known_indels_tbi',
            'aligner', 'skip_tools', 'save_mapped', 'intervals', 'tools', 'wes', 'save_output_as_bam'}
    lines = [line for line in Path(e['REPO'] + '/reference/params.tier1.yaml').read_text().splitlines()
             if line.split(':', 1)[0] not in omit]
    settings = {'step': 'variant_calling', 'tools': 'mutect2', 'intervals': str(bed), **{k: str(v) for k, v in resources.items()}}
    reference = Path(e.get('REF', str(Path.home() / 'omics/ref/GATK.GRCh38'))) / 'Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta'
    lines = [line for line in lines if line.split(':', 1)[0] not in {'fasta', 'fasta_fai', 'dict'}]
    settings.update(fasta=str(reference), fasta_fai=str(reference) + '.fai', dict=e['DICT'])
    lines += [f'{k}: {json.dumps(v)}' for k, v in settings.items()]
    lines += ['wes: true', 'save_output_as_bam: true', 'publish_dir_mode: link']
    (data / 'params.yaml').write_text('\n'.join(lines) + '\n')
    for name in ('samplesheet.csv', 'params.yaml'):
        (art / name).write_bytes((data / name).read_bytes())
    print(f'prepared {len(rows)} intervals')


if __name__ == '__main__':
    try:
        prepare()
    except (ValueError, OSError, KeyError) as exc:
        sys.exit(f'::error::{exc}')
