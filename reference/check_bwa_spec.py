#!/usr/bin/env python3
"""Check docs/spec/bwa.md section 13 claims against golden chunk BAMs.

usage: check_bwa_spec.py CHUNK.bam R1.fastq.gz R2.fastq.gz

Parses BAM directly (gzip handles BGZF members) so raw aux type bytes and the
bin field are visible. Prints one line per check: PASS/FAIL + counts.
"""
import gzip
import re
import struct
import sys
from collections import Counter, defaultdict

bam, fq1, fq2 = sys.argv[1:4]
results = []


def check(name, ok, detail=""):
    results.append((name, ok))
    print(f"{'PASS' if ok else 'FAIL'} {name} {detail}")


def reg2bin(beg, end):
    end -= 1
    if beg >> 14 == end >> 14:
        return ((1 << 15) - 1) // 7 + (beg >> 14)
    if beg >> 17 == end >> 17:
        return ((1 << 12) - 1) // 7 + (beg >> 17)
    if beg >> 20 == end >> 20:
        return ((1 << 9) - 1) // 7 + (beg >> 20)
    if beg >> 23 == end >> 23:
        return ((1 << 6) - 1) // 7 + (beg >> 23)
    if beg >> 26 == end >> 26:
        return ((1 << 3) - 1) // 7 + (beg >> 26)
    return 0


AUX_SIZE = {"A": 1, "c": 1, "C": 1, "s": 2, "S": 2, "i": 4, "I": 4, "f": 4}
AUX_FMT = {"c": "<b", "C": "<B", "s": "<h", "S": "<H", "i": "<i", "I": "<I", "f": "<f", "A": "c"}

data = gzip.open(bam).read()
assert data[:4] == b"BAM\1"
(l_text,) = struct.unpack_from("<i", data, 4)
text = data[8 : 8 + l_text].decode()
off = 8 + l_text
(n_ref,) = struct.unpack_from("<i", data, off)
off += 4
refs = []
for _ in range(n_ref):
    (ln,) = struct.unpack_from("<i", data, off)
    name = data[off + 4 : off + 4 + ln - 1].decode()
    (lr,) = struct.unpack_from("<i", data, off + 4 + ln)
    refs.append((name, lr))
    off += 8 + ln

recs = []
while off < len(data):
    (bs,) = struct.unpack_from("<i", data, off)
    b = data[off + 4 : off + 4 + bs]
    off += 4 + bs
    tid, pos, lrn, mapq, binv, ncig, flag, lseq, mtid, mpos, tlen = struct.unpack_from("<iiBBHHHiiii", b, 0)
    p = 32
    qname = b[p : p + lrn - 1].decode()
    p += lrn
    cig = []
    for i in range(ncig):
        (v,) = struct.unpack_from("<I", b, p + 4 * i)
        cig.append((v >> 4, "MIDNSHP=X"[v & 15]))
    p += 4 * ncig
    seq = "".join("=ACMGRSVTWYHKDBN"[(b[p + i // 2] >> (4 * (1 - i % 2))) & 15] for i in range(lseq))
    p += (lseq + 1) // 2
    qual = bytes(q + 33 for q in b[p : p + lseq]).decode() if lseq and b[p] != 0xFF else "*"
    p += lseq
    aux = []
    while p < len(b):
        tag = b[p : p + 2].decode()
        t = chr(b[p + 2])
        p += 3
        if t == "Z":
            e = b.index(0, p)
            val = b[p:e].decode()
            p = e + 1
        elif t in AUX_SIZE:
            val = struct.unpack_from(AUX_FMT[t], b, p)[0]
            p += AUX_SIZE[t]
        else:
            raise SystemExit(f"unexpected aux type {t}")
        aux.append((tag, t, val))
    recs.append(
        dict(tid=tid, pos=pos, mapq=mapq, bin=binv, flag=flag, mtid=mtid, mpos=mpos, tlen=tlen,
             qname=qname, cig=cig, seq=seq, qual=qual, aux=aux)
    )

print(f"# {bam}: {len(recs)} records, {n_ref} refs")
hdr = text.rstrip("\n").split("\n")

# 1-3 header
check("1 HD line", hdr[0] == "@HD\tVN:1.5\tSO:coordinate", repr(hdr[0]))
sq = [l for l in hdr if l.startswith("@SQ")]
check("2 SQ SN/LN only, 3366", len(sq) == 3366 and all(len(l.split("\t")) == 3 for l in sq))
pg = [l for l in hdr if l.startswith("@PG")]
ok = (len(pg) == 2 and pg[0].startswith("@PG\tID:bwa\tPN:bwa\tVN:0.7.19-r1273\tCL:")
      and pg[1].startswith("@PG\tID:samtools\tPN:samtools\tPP:bwa\tVN:1.22.1\tCL:samtools sort"))
check("2 two PG lines", ok, " | ".join(l[:60] for l in pg))
check("3 CL has literal \\t", "\\t" in pg[0] if pg else False)
order = [l[:3] for l in hdr]
check("2 header order HD,SQ*,RG,PG,PG", order == ["@HD"] + ["@SQ"] * 3366 + ["@RG", "@PG", "@PG"], str(Counter(order)))

# FASTQ
def readfq(p):
    out = {}
    order = []
    with gzip.open(p, "rt") as f:
        while True:
            h = f.readline()
            if not h:
                break
            s = f.readline().rstrip("\n")
            f.readline()
            q = f.readline().rstrip("\n")
            name = h[1:].split()[0]
            out[name] = (s, q, h)
            order.append(name)
    return out, order


r1, order1 = readfq(fq1)
r2, order2 = readfq(fq2)
check("20 R1/R2 names pair up", order1 == order2 and not any(re.search(r"/[0-9]$", n) for n in order1))
check("21 qual alphabet", set("".join(q for _, q, _ in r1.values())) <= set("#,:F"),
      "".join(sorted(set("".join(q for _, q, _ in r1.values())))))
idx = {n: i for i, n in enumerate(order1)}

# 5 names
c = Counter(r["qname"] for r in recs)
supp = Counter(r["qname"] for r in recs if r["flag"] & 0x800)
check("5 names = 2 + supplementaries", all(c[n] == 2 + supp[n] for n in c) and set(c) == set(order1),
      f"{len(c)} names vs {len(order1)} pairs")

# 6/7 aux
TAG_ORDER = ["NM", "MD", "MC", "MQ", "AS", "XS", "RG", "SA", "XA"]
bad_order = bad_type = 0
for r in recs:
    tags = [a[0] for a in r["aux"]]
    if [t for t in TAG_ORDER if t in tags] != tags:
        bad_order += 1
    for tag, t, v in r["aux"]:
        if tag in ("NM", "MQ", "AS", "XS"):
            want = "C" if 0 <= v <= 255 else "S" if v <= 65535 else "I" if v >= 0 else None
            if v < 0:
                want = "c" if v >= -128 else "s" if v >= -32768 else "i"
            if t != want:
                bad_type += 1
check("6 tag order", bad_order == 0, f"{bad_order} bad")
check("6 integer typing", bad_type == 0, f"{bad_type} bad")
types = Counter((a[0], a[1]) for r in recs for a in r["aux"])
print("   aux types:", dict(types))
missing = sum(1 for r in recs if not {"MQ", "AS", "RG", "XS"} <= {a[0] for a in r["aux"]})
check("7 MQ/AS/RG/XS everywhere", missing == 0, f"{missing} missing")
check("7 no secondary", not any(r["flag"] & 0x100 for r in recs))
um_bad = 0
mates = defaultdict(dict)
for r in recs:
    if not r["flag"] & 0x800:
        mates[r["qname"]][1 if r["flag"] & 0x40 else 2] = r
for r in recs:
    if r["flag"] & 4:
        a = {x[0]: x[2] for x in r["aux"]}
        m = mates[r["qname"]][2 if r["flag"] & 0x40 else 1]
        mate_mapped = not m["flag"] & 4
        if a.get("AS") != 0 or a.get("XS") != 0 or "NM" in a or "MD" in a or (("MC" in a) != mate_mapped):
            um_bad += 1
check("7 unmapped aux", um_bad == 0, f"{um_bad} bad")

# 8/9 supplementary
sups = [r for r in recs if r["flag"] & 0x800]
cig_h = sum(1 for r in recs if any(op == "H" for _, op in r["cig"]))
sa_h = sum(1 for r in recs for a in r["aux"] if a[0] == "SA" and "H" in re.sub(r"[^A-Z]", "", a[2].replace("+", "").replace("-", "")))
check("8 supplementaries exist with SA", len(sups) > 0 and all(any(a[0] == "SA" for a in s["aux"]) for s in sups), f"{len(sups)}")
pp_bad = sum(1 for s in sups if (s["flag"] & 2) != (mates[s["qname"]][1 if s["flag"] & 0x40 else 2]["flag"] & 2))
check("8 supp 0x2 equals primary", pp_bad == 0, f"{pp_bad} bad")
check("8 no H in CIGAR/SA", cig_h == 0)
mq_bad = sum(1 for s in sups if s["mapq"] > mates[s["qname"]][1 if s["flag"] & 0x40 else 2]["mapq"])
check("9 supp MAPQ <= primary", mq_bad == 0, f"{mq_bad} bad")

# 10/11 flags
pairs = Counter()
for n, m in mates.items():
    pairs[(m[1]["flag"], m[2]["flag"])] += 1
one_unmapped = {k: v for k, v in pairs.items() if (k[0] ^ k[1]) & 4}
both = {k: v for k, v in pairs.items() if k[0] & 4 and k[1] & 4}
# The spec lists mapped/unmapped flag pairs unordered; pairs here are (R1, R2).
allowed = {(73, 133), (121, 181), (137, 69), (185, 117)}
allowed |= {(b, a) for a, b in allowed}
check("10 one-end-unmapped flag pairs", set(one_unmapped) <= allowed, str(one_unmapped))
check("10 both-unmapped 77/141", set(both) <= {(77, 141)}, str(both))


def rc(s):
    return s[::-1].translate(str.maketrans("ACGTN", "TGCAN"))


def reflen(cig):
    return sum(n for n, op in cig if op in "MDN=X")


u_bad = 0
for n, m in mates.items():
    for i, j in ((1, 2), (2, 1)):
        r, o = m[i], m[j]
        if r["flag"] & 4 and not o["flag"] & 4:
            a = {x[0]: x[2] for x in r["aux"]}
            ocig = "".join(f"{l}{op}" for l, op in o["cig"])
            fs, fq, _ = (r1 if i == 1 else r2)[n]
            exp_seq, exp_q = (rc(fs), fq[::-1]) if o["flag"] & 0x10 else (fs, fq)
            if not (r["tid"] == o["tid"] and r["pos"] == o["pos"] and r["mapq"] == 0 and not r["cig"]
                    and r["mtid"] == r["tid"] and r["mpos"] == r["pos"] and r["tlen"] == 0
                    and a.get("MC") == ocig and a.get("MQ") == o["mapq"] and r["seq"] == exp_seq and r["qual"] == exp_q):
                u_bad += 1
check("11 mate-placed unmapped fields", u_bad == 0, f"{u_bad} bad")

# 12 TLEN
t_bad = t_n = 0
for n, m in mates.items():
    a, b2 = m[1], m[2]
    if a["flag"] & 4 or b2["flag"] & 4 or a["tid"] != b2["tid"]:
        continue
    def p5(r):
        return r["pos"] if not r["flag"] & 0x10 else r["pos"] + reflen(r["cig"]) - 1
    for r, o in ((a, b2), (b2, a)):
        p0, p1 = p5(r), p5(o)
        exp = p1 - p0 + 1 if p0 < p1 else p1 - p0 - 1 if p0 > p1 else 0
        t_n += 1
        if r["tlen"] != exp:
            t_bad += 1
check("12 TLEN formula (same-contig mapped pairs)", t_bad == 0, f"{t_bad}/{t_n} bad")

# 13 MD/NM
md_bad = 0
for r in recs:
    a = {x[0]: x[2] for x in r["aux"]}
    if "MD" not in a:
        continue
    md = a["MD"]
    mism = len(re.findall(r"[A-Z]", re.sub(r"\^[A-Z]+", "", md)))
    dels = sum(len(x) for x in re.findall(r"\^([A-Z]+)", md))
    ins = sum(n for n, op in r["cig"] if op == "I")
    if a["NM"] != mism + dels + ins or md.startswith("^") or re.search(r"\^[A-Z]+$", md):
        md_bad += 1
check("13 MD/NM consistent", md_bad == 0, f"{md_bad} bad")

# 14 SEQ vs FASTQ
s_bad = 0
for r in recs:
    if r["flag"] & 4:
        continue
    fs, fq, _ = (r1 if r["flag"] & 0x40 else r2)[r["qname"]]
    exp = rc(fs) if r["flag"] & 0x10 else fs
    if r["seq"] != exp.upper():
        s_bad += 1
check("14 mapped SEQ = FASTQ (rc if reverse)", s_bad == 0, f"{s_bad} bad")

# 15 XA
xa_bad = 0
xa_n = 0
for r in recs:
    for tag, t, v in r["aux"]:
        if tag == "XA":
            xa_n += 1
            ents = v.split(";")[:-1]
            if not v.endswith(";") or len(ents) > 5 or not all(re.fullmatch(r"[^,]+,[+-]\d+,[0-9MIDS]+,\d+", e) for e in ents):
                xa_bad += 1
check("15 XA format and <=5", xa_bad == 0, f"{xa_bad}/{xa_n} bad")

# 16 MAPQ
check("16 MAPQ <= 60", all(r["mapq"] <= 60 for r in recs))

# 17 sort order
def key(r):
    return (r["tid"] if r["tid"] >= 0 else 1 << 31, r["pos"], (r["flag"] >> 4) & 1)

o_bad = sum(1 for x, y in zip(recs, recs[1:]) if key(x) > key(y))
check("17 (tid,pos,strand) non-decreasing", o_bad == 0, f"{o_bad} bad")
# stability: within equal keys, emission order = (pair index, R1 before R2, primary before supplementary)
def emit(r):
    return (idx[r["qname"]], 0 if r["flag"] & 0x40 else 1)

st_bad = st_n = 0
for x, y in zip(recs, recs[1:]):
    if key(x) == key(y):
        st_n += 1
        if emit(x) > emit(y):
            st_bad += 1
check("17 ties in emission order", st_bad == 0, f"{st_bad}/{st_n} ties bad")

# 18 bin
b_bad = 0
for r in recs:
    if r["tid"] < 0:
        exp = 4680
    else:
        end = r["pos"] + (reflen(r["cig"]) or 1)
        exp = reg2bin(r["pos"], end)
    if r["bin"] != exp:
        b_bad += 1
check("18 bin", b_bad == 0, f"{b_bad} bad")

nf = sum(1 for _, ok in results if not ok)
print(f"# {len(results) - nf}/{len(results)} checks passed")
sys.exit(1 if nf else 0)
