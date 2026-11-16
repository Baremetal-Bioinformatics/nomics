"""Trace-driven golden promotion; no external Python dependencies."""
import argparse
import collections
import csv
import hashlib
import os
from pathlib import Path
import re
import shutil
import tempfile


def harvest(run, force=False):
    run = Path(run).resolve()
    golden = run / 'golden'
    if golden.exists() and not force:
        raise ValueError(f'{golden} exists (use --force to replace)')
    traces = sorted((run / 'results/pipeline_info').glob('execution_trace*.txt'))
    if len(traces) != 1:
        raise ValueError(f'expected one execution trace, found {len(traces)}; select one before harvesting')
    counts = collections.Counter()
    tasks = []
    with traces[0].open() as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            name = row['name']
            if 'BAM_VARIANT_CALLING_SOMATIC_MUTECT2:' not in name or row['status'] != 'COMPLETED':
                continue
            match = re.fullmatch(r'.*:([^:]+) \(([^/]+)\)', name)
            if not match or not re.fullmatch(r'[0-9a-f]+/[0-9a-f]+', row['hash']):
                raise ValueError(f'invalid task name/hash: {row}')
            process, tag = match.groups()
            base = f'{process}-{tag}'
            if tag in ('.', '..') or '\n' in tag:
                raise ValueError(f'unsafe tag: {tag}')
            counts[base] += 1
            label = base + (f'.{counts[base]}' if counts[base] > 1 else '')
            prefix = run / 'work' / row['hash']
            dirs = sorted(p for p in prefix.parent.glob(prefix.name + '*') if p.is_dir())
            if len(dirs) != 1:
                raise ValueError(f'{prefix}: expected one work directory, found {len(dirs)}')
            tasks.append((label, process, name, dirs[0]))
    if not tasks or not any(t[1] == 'MUTECT2_PAIRED' for t in tasks):
        raise ValueError('no completed paired Mutect2 tasks')
    stage = Path(tempfile.mkdtemp(prefix='.golden-', dir=run))
    try:
        intervals = stage / 'intervals'
        intervals.mkdir()
        for label, process, name, work in tasks:
            dest = stage / label
            dest.mkdir()
            links = []
            beds = 0
            for f in sorted(work.iterdir()):
                if f.is_symlink():
                    links.append(f'{f.name} -> {os.readlink(f)}\n')
                    if process == 'MUTECT2_PAIRED' and f.suffix == '.bed':
                        beds += 1
                        target = intervals / f.name
                        if target.exists() and target.read_bytes() != f.read_bytes():
                            raise ValueError(f'conflicting interval contents: {f.name}')
                        shutil.copyfile(f, target)
                elif f.is_file():
                    if f.name in ('TASK', 'INPUTS'):
                        raise ValueError(f'reserved filename: {f}')
                    shutil.copy2(f, dest / f.name)
            if process == 'MUTECT2_PAIRED' and beds != 1:
                raise ValueError(f'{work}: expected one staged BED, found {beds}')
            (dest / 'TASK').write_text(name + '\n')
            (dest / 'INPUTS').write_text(''.join(links))
        with (stage / 'SHA256SUMS').open('w') as out:
            for f in sorted(stage.rglob('*')):
                if f.is_file() and f != stage / 'SHA256SUMS':
                    h = hashlib.sha256()
                    with f.open('rb') as fh:
                        for block in iter(lambda: fh.read(1024 * 1024), b''):
                            h.update(block)
                    out.write(f'{h.hexdigest()}  ./{f.relative_to(stage)}\n')
        if golden.exists():
            if golden.is_symlink():
                raise ValueError('refusing to replace a symlinked golden')
            shutil.rmtree(golden)
        stage.rename(golden)
    finally:
        if stage.exists():
            shutil.rmtree(stage)
    print(f'harvested {len(tasks)} tasks into {golden}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run')
    parser.add_argument('--force', action='store_true')
    args = parser.parse_args()
    try:
        harvest(args.run, args.force)
    except (ValueError, OSError, KeyError) as exc:
        parser.exit(2, f'::error::{exc}\n')
