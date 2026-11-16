"""Compare unpublished scatter outputs, without hiding missing or ambiguous matches."""
import gzip
import hashlib
from pathlib import Path
import sys
import tarfile


def digest(path):
    h = hashlib.md5()
    if path.name.endswith('.tar.gz'):
        with tarfile.open(path, 'r:gz') as archive:
            for member in archive:
                if member.isfile():
                    h.update(archive.extractfile(member).read())
    else:
        with gzip.open(path, 'rb') as fh:
            for line in fh:
                if not line.startswith(b'#'):
                    h.update(line)
    return h.digest()


def files(run):
    result = {}
    for path in sorted((Path(run) / 'work').glob('*/*/*')):
        if path.is_symlink() or not path.is_file():
            continue
        if 'f1r2' in path.name and path.name.endswith('.tar.gz') or 'mutect2' in path.name and path.name.endswith('.vcf.gz'):
            if path.name in result:
                raise ValueError(f'ambiguous work output: {path.name}')
            result[path.name] = path
    if not result:
        raise ValueError(f'no scatter outputs in {run}')
    return result


a, b = map(files, sys.argv[1:])
for name in sorted(a.keys() | b.keys()):
    status = 'FAIL' if name not in a or name not in b else ('PASS' if digest(a[name]) == digest(b[name]) else 'DIFF')
    print(f'work\t{name}\tcontents\t{status}\t')
