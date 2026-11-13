# Re-header tier0 to NovaSeq-style names and split into two lanes (deterministic).
import gzip, os, sys
src = os.path.expanduser("~/omics/data/tier0/tier0_%d.fastq.gz")
dst = os.path.expanduser("~/omics/data/tier0nv/tier0nv_L%d_%d.fastq.gz")
os.makedirs(os.path.dirname(dst % (1, 1)), exist_ok=True)
ins = [gzip.open(src % m, "rt") for m in (1, 2)]
outs = {(l, m): gzip.open(dst % (l, m), "wt", compresslevel=1) for l in (1, 2) for m in (1, 2)}
i = 0
while True:
    recs = [[f.readline() for _ in range(4)] for f in ins]
    if not recs[0][0]:
        break
    lane = 1 if i % 2 == 0 else 2
    tile = 1101 + (i // 2) % 8
    # x, y in a narrow field so duplicates can land within Picard's 100-pixel distance.
    h = (i * 2654435761) & 0xffffffff
    x, y = 1000 + h % 3000, 1000 + (h >> 12) % 3000
    name = "@A00744:46:HV3GHDSXX:%d:%d:%d:%d" % (lane, tile, x, y)
    for m, rec in enumerate(recs, start=1):
        rec[0] = "%s %d:N:0:ATCACGAT+TTAGGCAT\n" % (name, m)
        outs[(lane, m)].writelines(rec)
    i += 1
for f in list(ins) + list(outs.values()):
    f.close()
print("pairs", i)
