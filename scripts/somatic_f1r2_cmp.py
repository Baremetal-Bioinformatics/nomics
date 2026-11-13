#!/usr/bin/env python3
"""Compare a Mutect2 f1r2 archive with GATK's under the rule in docs/GATES.md.

  somatic_f1r2_cmp.py GOLDEN.tar.gz OURS.tar.gz LABEL

Member names as a set (GATK's order is its temporary directory's listing order);
ref_histogram and alt_table by bytes; alt_histogram by column label (GATK's column order is
a hash-map order), every value, every other line and the line endings by bytes.
"""
import sys,tarfile

def members(path):
    with tarfile.open(path,'r:gz') as t:
        leaves=t.getmembers()
        names=[m.name for m in leaves]
        assert len(names)==len(set(names))==3, f'{path}: expected three unique members'
        assert all(m.isfile() for m in leaves), f'{path}: non-file member'
        return {m.name:t.extractfile(m).read() for m in leaves}

def compare_alt(a,b):
    x,y=a.splitlines(keepends=True),b.splitlines(keepends=True)
    assert len(x)==len(y), 'alt histogram line count'
    ix=next(i for i,l in enumerate(x) if l.startswith(b'depth\t'))
    assert next(i for i,l in enumerate(y) if l.startswith(b'depth\t'))==ix
    def fields(line):
        assert line.endswith(b'\n') and not line.endswith(b'\r\n'), 'line ending'
        return line[:-1].split(b'\t')
    hx,hy=fields(x[ix]),fields(y[ix])
    assert len(hx)==len(hy)==385 and len(set(hx))==385 and set(hx)==set(hy), 'alt histogram column labels'
    order=[hy.index(c) for c in hx]
    rows=0
    for i,(lx,ly) in enumerate(zip(x,y)):
        if i==ix: continue
        if i>ix and lx.strip():
            fx,fy=fields(lx),fields(ly)
            assert len(fx)==len(fy)==385, f'alt row {i}: width'
            for j,col in enumerate(hx):
                assert fx[j]==fy[order[j]], f'alt line {i+1}, column {col!r}: {fx[j]!r} != {fy[order[j]]!r}'
            rows+=1
        else: assert lx==ly, f'alt non-column line {i+1}: {lx!r} != {ly!r}'
    assert rows==200, f'alt histogram: {rows} depth rows'
    return rows

a,b=members(sys.argv[1]),members(sys.argv[2])
assert set(a)==set(b), f'member sets: {set(a)^set(b)}'
result=[]
for name in a:
    if name.endswith('.alt_histogram'):
        rows=compare_alt(a[name],b[name]); result.append(f'alt_histogram 384 columns x {rows} rows exact')
    else:
        if a[name]!=b[name]:
            offset=next((i for i,(x,y) in enumerate(zip(a[name],b[name])) if x!=y),min(len(a[name]),len(b[name])))
            raise AssertionError(f'{name}: first difference byte {offset}; golden {a[name][max(0,offset-35):offset+70]!r}, port {b[name][max(0,offset-35):offset+70]!r}; lengths {len(a[name])}/{len(b[name])}')
        result.append(f'{name.rsplit(".",1)[-1]} {len(a[name])} bytes exact')
print(f'PASS {sys.argv[3]}: member set exact; '+ '; '.join(sorted(result)))
