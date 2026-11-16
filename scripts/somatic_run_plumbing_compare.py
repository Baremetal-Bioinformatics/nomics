from pathlib import Path
import subprocess as sp
import sys
D=Path(sys.argv[1]); B=Path(sys.argv[2])
count=0
for s in ['HCC1395BL_WES','HCC1395_WES']:
    a=D/s; p=D/'pair/results/preprocessing'
    for folder, dest in [('prep',D/f'pair/work/preprocessing/{s}/prep'), ('mapped',p/'mapped'/s)]:
        for f in (a/folder).iterdir():
            if f.suffix == '.tsv': continue
            assert f.read_bytes()==(dest/f.name).read_bytes(), f
            count+=1
    assert (a/f'{s}.recal.table').read_bytes()==(p/f'recal_table/{s}/{s}.recal.table').read_bytes()
    count+=1
    for stage,kind in [('markduplicates','md'),('recalibrated','recal')]:
        sp.run([str(B),'dev','gate-bam','--ours',str(p/f'{stage}/{s}/{s}.{kind}.bam'),'--golden',str(a/f'{s}.{kind}.bam')],check=True)
    x=(a/f'{s}.md.metrics').read_text().splitlines();y=(D/f'pair/results/reports/markduplicates/{s}/{s}.md.bam.metrics').read_text().splitlines()
    assert [v for i,v in enumerate(x) if i not in (1,3)]==[v for i,v in enumerate(y) if i not in (1,3)]
for folder,other in [(D/'chain/variant_calling',D/'pair/results/variant_calling'),(D/'chain/work',D/'pair/work/chain')]:
    files={f.relative_to(folder) for f in folder.rglob('*') if f.is_file()}
    assert files=={f.relative_to(other) for f in other.rglob('*') if f.is_file()}
    for f in files:
        assert (folder/f).read_bytes()==(other/f).read_bytes(), f
        count+=1
assert count == 90, f'missing or unexpected deterministic comparisons: {count}'
print(f'PASS {count} byte comparisons; 4 BAM gates; 2 normalized metrics comparisons')
