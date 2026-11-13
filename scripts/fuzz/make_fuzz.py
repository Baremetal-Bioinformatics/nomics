#!/usr/bin/env python3
"""Synthetic read pairs aimed at one pipeline edge case, for the fuzz goldens.

    make_fuzz.py SCENARIO SEED PAIRS FASTA TARGETS_BED OUTDIR

Writes OUTDIR/{R1,R2}.fastq.gz and OUTDIR/targets.bed (the windows the reads came from, sorted
in FASTA order and merged), and prints a summary line. Everything is a function of
(SCENARIO, SEED, PAIRS, FASTA, TARGETS_BED): the same inputs give the same bytes.

Every scenario mixes PAIRS background pairs (ordinary 2x151 reads at ~30x over target windows)
with the reads that aim at its edge case. The sample is diploid everywhere: each window carries
one to three planted het or hom variants (Gen.plant), so HaplotypeCaller has calls to make in
every scenario and the edge-case reads sit next to them. Scenarios (see
docs/strategy/fuzz-goldens.md):

  homopolymer    indel errors and true indels in homopolymer runs
  long_indel     het and hom deletions of 15-120 bp and insertions of 10-60 bp
  short_insert   fragments shorter than the reads (adapter read-through, dovetails), some < 15 bp
  n_bases        N runs: leading, trailing, internal, at fastp's 5-N limit and over it, all-N
  duplicates     PCR and optical duplicates around Picard's 100-pixel distance, ties on quality
  chimera        mates on different contigs, split reads, inversions, 10 kb and 1 Mb inserts
  unmapped       random and low-complexity mates, both mates unmapped, adapter dimers
  contig_edges   chrM across its origin, contig starts and ends, N gaps, alt, decoy and EBV contigs
  quality        Q0-Q41 extremes, fastp's 40% low-quality limit, poly-G tails of 5-100 bp
  read_length    lengths 1-151, R1 and R2 of different lengths, around fastp's 15 bp minimum
  pileup         three loci at 3,000-8,000x with het SNPs and indels
  complex_var    SNP clusters, MNPs, two alt alleles at one site, indels next to SNPs, STRs
  overlap        fragment == read length, start-tie pairs, mates fully overlapping
  coord_wrap     tile, x and y above 32767 and 65535 (Picard keeps them as Java shorts)
"""
import gzip
import os
import random
import sys

COMP = str.maketrans("ACGTN", "TGCAN")
ADAPTER1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
ADAPTER2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA"
L = 151


def rc(s):
    return s.translate(COMP)[::-1]


class Fasta:
    """Random access through the .fai index."""

    def __init__(self, path):
        self.f = open(path, "rb")
        self.idx = {}
        self.order = []
        for line in open(path + ".fai"):
            name, length, offset, bases, width = line.split("\t")[:5]
            self.idx[name] = (int(length), int(offset), int(bases), int(width))
            self.order.append(name)

    def length(self, name):
        return self.idx[name][0]

    def seq(self, name, start, end):  # 0-based, half open
        length, offset, bases, width = self.idx[name]
        start, end = max(0, start), min(length, end)
        if end <= start:
            return ""
        a = offset + (start // bases) * width + start % bases
        b = offset + (end // bases) * width + end % bases
        self.f.seek(a)
        return self.f.read(b - a).replace(b"\n", b"").replace(b"\r", b"").decode().upper()


class Gen:
    def __init__(self, scenario, seed, fasta, bed):
        self.r = random.Random(f"{scenario}:{seed}")
        self.fa = fasta
        self.targets = []
        for line in open(bed):
            if line.startswith(("#", "track", "browser")) or not line.strip():
                continue
            c, s, e = line.split("\t")[:3]
            if c in fasta.idx and int(e) - int(s) >= 120:
                self.targets.append((c, int(s), int(e)))
        self.windows = []
        self.pairs = []
        self.n = 0
        self.names = "A00744:46:HV3GHDSXX"
        self.used = set()

    # ---- read names: unique (lane, tile, x, y) per pair unless a scenario sets them ----
    def name(self):
        """A NovaSeq-style name: tiles 1101-2478 (surface, swath, tile), x and y in 1000-31000."""
        self.n += 1
        i = self.n
        tile = (1 + i % 2) * 1000 + (1 + (i // 2) % 4) * 100 + 1 + (i // 8) % 78
        x = 1000 + (i * 7919) % 30000
        y = 1000 + (i * 104729 + i // 30000) % 30000
        return f"{self.names}:{1 + i % 2}:{tile}:{x}:{y}"

    def qual(self, n, kind="novaseq"):
        r = self.r
        if kind == "novaseq":
            return "".join("F" if r.random() < 0.86 else (":" if r.random() < 0.7 else ",") for _ in range(n))
        return "".join(chr(33 + min(41, max(2, int(r.gauss(36, 4))))) for _ in range(n))

    def errors(self, s, rate=0.002):
        r = self.r
        out = list(s)
        for i in range(len(out)):
            if r.random() < rate:
                out[i] = r.choice([b for b in "ACGT" if b != out[i]])
        return "".join(out)

    def window(self, pad=300, clean=True, tries=200):
        """A target interval with padding; with clean, one without N bases."""
        for _ in range(tries):
            c, s, e = self.r.choice(self.targets)
            s, e = max(0, s - pad), min(self.fa.length(c), e + pad)
            seq = self.fa.seq(c, s, e)
            if not clean or "N" not in seq:
                self.windows.append((c, s, e))
                return c, s, e, seq
        raise SystemExit("no clean target window found")

    def add(self, r1, q1, r2, q2, name=None):
        name = name or self.name()
        # Names stay unique: a clash moves the pair to another lane, which Picard does not read
        # (its location is tile, x and y), so optical-duplicate geometry is unchanged.
        while name in self.used:
            f = name.split(":")
            f[3] = str(int(f[3]) + 4)
            name = ":".join(f)
        self.used.add(name)
        self.pairs.append((name, r1, q1, r2, q2))
        return name

    def pair_from(self, hap, frag_mean=320, frag_sd=60, frag=None, rl=L, rate=0.002, name=None, qkind="novaseq"):
        """One FR pair from a haplotype sequence; short fragments read into the adapter."""
        r = self.r
        flen = frag if frag is not None else max(30, int(r.gauss(frag_mean, frag_sd)))
        flen = min(flen, len(hap))
        start = r.randrange(0, len(hap) - flen + 1)
        f = hap[start:start + flen]
        if r.random() < 0.5:
            f = rc(f)
        a = (f + ADAPTER1 + "A" * rl)[:rl]
        b = (rc(f) + ADAPTER2 + "A" * rl)[:rl]
        return self.add(self.errors(a, rate), self.qual(rl, qkind), self.errors(b, rate), self.qual(rl, qkind), name)

    def background(self, pairs):
        per = 120
        while pairs > 0:
            _, _, _, seq = self.window()
            haps = self.plant(seq)
            for _ in range(min(per, pairs)):
                self.pair_from(self.r.choice(haps))
            pairs -= per

    # ---- haplotypes ----
    def mutate(self, seq, variants):
        """Apply (pos, ref, alt) variants (positions in seq, non-overlapping) right to left."""
        out = seq
        for pos, ref, alt in sorted(variants, reverse=True):
            assert out[pos:pos + len(ref)] == ref, (pos, ref, out[pos:pos + len(ref)])
            out = out[:pos] + alt + out[pos + len(ref):]
        return out

    def diploid(self, seq, variants, depth, **kw):
        """variants: (pos, ref, alt, genotype) with genotype 'het' or 'hom' (or 'het2': second alt)."""
        h1 = self.mutate(seq, [(p, rf, a) for p, rf, a, g in variants if g in ("hom", "het")])
        h2 = self.mutate(seq, [(p, rf, a) for p, rf, a, g in variants if g in ("hom", "het2")])
        n = int(depth * len(seq) / (2 * L))
        for _ in range(n):
            self.pair_from(h1 if self.r.random() < 0.5 else h2, **kw)

    def plant(self, seq, indels=True):
        """The sample's two haplotypes over seq: one to three planted variants, het or hom, at
        least 40 bases apart, 60 from the ends and 10 from any N. SNVs, and with indels also
        deletions and insertions of 1-8 bases. Without indels both haplotypes keep seq's length,
        so a scenario can index them as it indexes seq."""
        r = self.r
        ok = [p for p in range(60, len(seq) - 70) if "N" not in seq[p - 10:p + 20]]
        if not ok:
            return seq, seq
        variants, last = [], -1000
        for p in sorted(r.sample(ok, min(len(ok), r.randrange(1, 4)))):
            if p - last < 40:
                continue
            gt = r.choice(["het", "het", "hom"])
            kind = r.randrange(4) if indels else 0
            if kind <= 1:
                variants.append((p, seq[p], r.choice([b for b in "ACGT" if b != seq[p]]), gt))
            elif kind == 2:
                variants.append((p, seq[p:p + 1 + r.randrange(1, 9)], seq[p], gt))
            else:
                variants.append((p, seq[p], seq[p] + "".join(r.choice("ACGT") for _ in range(r.randrange(1, 9))), gt))
            last = p + 10
        h1 = self.mutate(seq, [(p, rf, a) for p, rf, a, g in variants])
        h2 = self.mutate(seq, [(p, rf, a) for p, rf, a, g in variants if g == "hom"])
        return h1, h2

    def cover(self, haps, depth=30, **kw):
        """Ordinary pairs at `depth` over a window's two haplotypes."""
        for _ in range(int(depth * len(haps[0]) / (2 * L))):
            self.pair_from(self.r.choice(haps), **kw)

    def write(self, outdir):
        os.makedirs(outdir, exist_ok=True)
        self.r.shuffle(self.pairs)
        f1 = gzip.GzipFile(os.path.join(outdir, "R1.fastq.gz"), "wb", compresslevel=6, mtime=0)
        f2 = gzip.GzipFile(os.path.join(outdir, "R2.fastq.gz"), "wb", compresslevel=6, mtime=0)
        for name, r1, q1, r2, q2 in self.pairs:
            assert len(r1) == len(q1) and len(r2) == len(q2), name
            f1.write(f"@{name} 1:N:0:ATCACGAT+TTAGGCAT\n{r1}\n+\n{q1}\n".encode())
            f2.write(f"@{name} 2:N:0:ATCACGAT+TTAGGCAT\n{r2}\n+\n{q2}\n".encode())
        f1.close()
        f2.close()
        rank = {c: i for i, c in enumerate(self.fa.order)}
        merged = []
        for c, s, e in sorted(set(self.windows), key=lambda w: (rank[w[0]], w[1], w[2])):
            if merged and merged[-1][0] == c and s <= merged[-1][2]:
                merged[-1][2] = max(merged[-1][2], e)
            else:
                merged.append([c, s, e])
        with open(os.path.join(outdir, "targets.bed"), "w") as f:
            for c, s, e in merged:
                f.write(f"{c}\t{s}\t{e}\n")
        return len(self.pairs), len(merged)


# ---------------------------------------------------------------------------------------------
def homopolymer(g):
    found = 0
    for _ in range(4000):
        if found >= 60:
            break
        c, s, e, seq = g.window(pad=400)
        runs = [(i, j) for i, j in runs_of(seq) if j - i >= 8 and 200 < i < len(seq) - 200]
        if not runs:
            g.windows.pop()
            continue
        found += 1
        i, j = g.r.choice(runs)
        kind = g.r.choice(["del1", "ins1", "del3", "none"])
        var = []
        if kind == "del1":
            var = [(i - 1, seq[i - 1:i + 1], seq[i - 1], "het")]
        elif kind == "ins1":
            var = [(i - 1, seq[i - 1], seq[i - 1] + seq[i], g.r.choice(["het", "hom"]))]
        elif kind == "del3":
            var = [(i - 1, seq[i - 1:i + 3], seq[i - 1], "het")]
        h1 = g.mutate(seq, [(p, rf, a) for p, rf, a, gt in var])
        h2 = g.mutate(seq, [(p, rf, a) for p, rf, a, gt in var if gt == "hom"])
        for _ in range(int(40 * len(seq) / (2 * L))):
            hap = h1 if g.r.random() < 0.5 else h2
            # Sequencing slippage: a read-level +-1..3 base change of the run on 15% of reads.
            if g.r.random() < 0.15:
                k = g.r.choice([-3, -2, -1, 1, 2])
                b = seq[i]
                hap = hap[:i] + (b * k if k > 0 else "") + hap[i - k if k < 0 else i:]
            g.pair_from(hap)


def runs_of(seq):
    i = 0
    while i < len(seq):
        j = i
        while j < len(seq) and seq[j] == seq[i]:
            j += 1
        yield i, j
        i = j


def long_indel(g):
    for _ in range(80):
        c, s, e, seq = g.window(pad=500)
        p = g.r.randrange(250, len(seq) - 250)
        if g.r.random() < 0.55:
            n = g.r.randrange(15, 121)
            var = (p, seq[p:p + n + 1], seq[p], g.r.choice(["het", "het", "hom"]))
        else:
            n = g.r.randrange(10, 61)
            ins = seq[p + 1:p + 1 + n] if g.r.random() < 0.5 else "".join(g.r.choice("ACGT") for _ in range(n))
            var = (p, seq[p], seq[p] + ins, g.r.choice(["het", "het", "hom"]))
        g.diploid(seq, [var], 40)


def short_insert(g):
    for _ in range(60):
        c, s, e, seq = g.window()
        haps = g.plant(seq, indels=False)
        for _ in range(150):
            g.pair_from(g.r.choice(haps), frag=g.r.choice([g.r.randrange(1, 15), g.r.randrange(15, 60), g.r.randrange(60, 152)]))


def n_bases(g):
    for _ in range(50):
        c, s, e, seq = g.window()
        haps = g.plant(seq, indels=False)
        for _ in range(200):
            name = g.name()
            g.pair_from(g.r.choice(haps), name=name)
            nm, r1, q1, r2, q2 = g.pairs.pop()
            kind = g.r.randrange(8)
            if kind == 0:
                r1 = "N" * g.r.randrange(1, 6) + r1[5:]
                r1 = r1[:L].ljust(L, "A")
            elif kind == 1:
                k = g.r.randrange(1, 30)
                r2 = r2[:L - k] + "N" * k
            elif kind == 2:  # exactly fastp's limit (5) and one over
                k = g.r.choice([5, 6])
                pos = g.r.sample(range(L), k)
                r1 = "".join("N" if i in pos else b for i, b in enumerate(r1))
            elif kind == 3:
                p = g.r.randrange(10, 120)
                k = g.r.randrange(1, 21)
                r1 = r1[:p] + "N" * k + r1[p + k:]
            elif kind == 4:
                r1, r2 = "N" * L, "N" * L
            elif kind == 5:
                r2 = "N" * L
            # N bases carry the lowest quality, as instruments write them.
            q1 = "".join("#" if b == "N" else q for b, q in zip(r1, q1))
            q2 = "".join("#" if b == "N" else q for b, q in zip(r2, q2))
            g.add(r1, q1, r2, q2, nm)


def duplicates(g):
    for k in range(400):
        c, s, e, seq = g.window()
        haps = g.plant(seq, indels=False)
        if k % 4 == 0:  # a quarter of the duplicate sets sit in called regions
            g.cover(haps)
        g.pair_from(g.r.choice(haps))
        name, r1, q1, r2, q2 = g.pairs[-1]
        lane, tile, x, y = (int(t) for t in name.split(":")[3:7])
        copies = g.r.choice([1, 1, 2, 3, 5, 10, 50])
        for j in range(copies):
            kind = g.r.randrange(6)
            if kind == 0:  # optical: same tile, distance below, at and above 100 pixels
                d = g.r.choice([1, 50, 99, 100, 101, 150])
                nx, ny, nt = x + d, y, tile
            elif kind == 1:
                nx, ny, nt = x + g.r.randrange(-70, 71), y + g.r.randrange(-70, 71), tile
            elif kind == 2:  # same x, y on another tile: not optical
                nx, ny, nt = x, y, tile + 1
            else:  # PCR duplicate anywhere
                nx, ny, nt = g.r.randrange(1000, 30000), g.r.randrange(1000, 30000), 1101 + g.r.randrange(4) * 100
            # Quality ties and near-ties decide which copy MarkDuplicates keeps.
            qq1 = q1 if g.r.random() < 0.5 else g.qual(L)
            qq2 = q2 if g.r.random() < 0.5 else g.qual(L)
            g.add(r1, qq1, r2, qq2, f"{g.names}:{lane}:{nt}:{max(0, nx)}:{max(0, ny)}")


def chimera(g):
    start = len(g.pairs)
    for k in range(4000):
        a = g.window()
        b = g.window()
        ha, hb = g.plant(a[3], indels=False), g.plant(b[3], indels=False)
        if k % 20 == 0:  # some chimeric reads land in called regions
            g.cover(ha)
        fa, fb = g.r.choice(ha), g.r.choice(hb)
        kind = g.r.randrange(4)
        if kind == 0:  # mates on different contigs (or far apart)
            p, q = g.r.randrange(0, len(fa) - L), g.r.randrange(0, len(fb) - L)
            g.add(g.errors(fa[p:p + L]), g.qual(L), g.errors(rc(fb[q:q + L])), g.qual(L))
        elif kind == 1:  # split read: R1 joins two loci
            p, q = g.r.randrange(0, len(fa) - L), g.r.randrange(0, len(fb) - L)
            cut = g.r.randrange(30, 121)
            r1 = fa[p:p + cut] + fb[q:q + L - cut]
            g.add(g.errors(r1), g.qual(L), g.errors(rc(fa[p + 200:p + 200 + L]).ljust(L, "A")[:L]), g.qual(L))
        elif kind == 2:  # inversion: both mates on the same strand
            p = g.r.randrange(0, len(fa) - 400)
            g.add(g.errors(fa[p:p + L]), g.qual(L), g.errors(fa[p + 200:p + 200 + L].ljust(L, "A")[:L]), g.qual(L))
        else:  # large insert: 10 kb or 1 Mb on the same contig
            c, s, e = a[:3]
            gap = g.r.choice([10_000, 1_000_000])
            far = g.fa.seq(c, e + gap, e + gap + L)
            p = g.r.randrange(0, len(fa) - L)
            if len(far) == L and "N" not in far:
                g.windows.append((c, e + gap - 50, e + gap + L + 50))
                g.add(g.errors(fa[p:p + L]), g.qual(L), g.errors(rc(far)), g.qual(L))
        if len(g.pairs) - start >= 3000:
            break


def unmapped(g):
    g.names = "K00167:276:HJ2KYBBXY"  # HiSeq 4000 names: poly-G off
    for i in range(7500):
        if i % 5 == 0:
            c, s, e, ref = g.window()
            haps = g.plant(ref, indels=False)
            if i % 100 == 0:  # some half-mapped pairs land in called regions
                g.cover(haps, qkind="hiseq")
        seq = g.r.choice(haps)
        p = g.r.randrange(0, len(seq) - L)
        mapped = g.errors(seq[p:p + L])
        rand = "".join(g.r.choice("ACGT") for _ in range(L))
        low = (g.r.choice(["A", "AT", "CA", "GGC"]) * L)[:L]
        dimer = (ADAPTER1 + "A" * L)[:L]
        kind = g.r.randrange(5)
        if kind == 0:
            g.add(mapped, g.qual(L, "hiseq"), rand, g.qual(L, "hiseq"))
        elif kind == 1:
            g.add(rand, g.qual(L, "hiseq"), "".join(g.r.choice("ACGT") for _ in range(L)), g.qual(L, "hiseq"))
        elif kind == 2:
            g.add(low, g.qual(L, "hiseq"), mapped, g.qual(L, "hiseq"))
        elif kind == 3:
            g.add(dimer, g.qual(L, "hiseq"), (ADAPTER2 + "A" * L)[:L], g.qual(L, "hiseq"))
        else:
            g.add(mapped, g.qual(L, "hiseq"), rc(seq[p + 150:p + 150 + L]).ljust(L, "A")[:L], g.qual(L, "hiseq"))


def contig_edges(g):
    fa = g.fa
    regions = []
    if "chrM" in fa.idx:  # circular: fragments across the origin are end + start
        m = fa.length("chrM")
        haps = g.plant(fa.seq("chrM", 0, m), indels=False)
        for _ in range(600):
            whole = g.r.choice(haps)
            g.pair_from(whole[-600:] + whole[:600], frag_mean=300, frag_sd=40)
        for _ in range(800):
            g.pair_from(g.r.choice(haps))
        g.windows += [("chrM", 0, m)]
    for c in [x for x in fa.order if not x.startswith(("chrUn", "HLA")) and "_" not in x][:25]:
        n = fa.length(c)
        # The first and last bases that are not N, and the N gap edge next to them.
        for lo, hi in ((0, 200_000), (max(0, n - 200_000), n)):
            seq = fa.seq(c, lo, hi)
            idx = [i for i in range(1, len(seq)) if (seq[i] == "N") != (seq[i - 1] == "N")]
            if not idx:
                continue
            for i in (idx[0], idx[-1]):
                a, b = max(0, i - 600), min(len(seq), i + 600)
                regions.append((c, lo + a, lo + b, seq[a:b]))
    alts = [x for x in fa.order if x.endswith("_alt")][:6]
    decoys = [x for x in fa.order if x.endswith("_decoy")][:6]
    others = [x for x in ("chrEBV",) if x in fa.idx]
    for c in alts + decoys + others:
        n = fa.length(c)
        a = g.r.randrange(0, max(1, n - 3000))
        regions.append((c, a, min(n, a + 3000), fa.seq(c, a, min(n, a + 3000))))
    for c, a, b, seq in regions:
        g.windows.append((c, a, b))
        # Reads overlapping N stretches are kept as the reference has them (N in the read).
        g.cover(g.plant(seq, indels=False), frag_mean=280, frag_sd=50)


def quality(g):
    for _ in range(60):
        c, s, e, seq = g.window()
        haps = g.plant(seq, indels=False)
        for _ in range(200):
            g.pair_from(g.r.choice(haps))
            name, r1, q1, r2, q2 = g.pairs.pop()
            kind = g.r.randrange(7)
            if kind == 0:
                q1 = "#" * L
            elif kind == 1:  # fastp: unqualified (< Q15) bases at 40% and just over
                k = g.r.choice([60, 61])
                pos = set(g.r.sample(range(L), k))
                q1 = "".join("+" if i in pos else q for i, q in enumerate(q1))
            elif kind == 2:
                q2 = q2[:100] + "#" * 51
            elif kind == 3:
                q1 = "".join(chr(33 + g.r.randrange(0, 42)) for _ in range(L))
            elif kind == 4:  # two-colour no-signal: poly-G tails
                k = g.r.choice([5, 9, 10, 11, 30, 100, L])
                r1 = r1[:L - k] + "G" * k
            elif kind == 5:
                k = g.r.choice([10, 40])
                r2 = r2[:L - k] + "G" * k
            else:
                q1 = "J" * L
            g.add(r1, q1, r2, q2, name)


def read_length(g):
    g.names = "K00167:276:HJ2KYBBXY"
    for _ in range(60):
        c, s, e, seq = g.window()
        haps = g.plant(seq, indels=False)
        for _ in range(200):
            g.pair_from(g.r.choice(haps), qkind="hiseq")
            name, r1, q1, r2, q2 = g.pairs.pop()
            l1 = g.r.choice([1, 14, 15, 16, 30, 75, 100, 150, L])
            l2 = g.r.choice([l1, l1, g.r.randrange(1, L + 1)])
            g.add(r1[:l1], q1[:l1], r2[:l2], q2[:l2], name)


def pileup(g):
    for depth in (3000, 5000, 8000):
        c, s, e, seq = g.window(pad=200)
        p = len(seq) // 2
        var = [(p, seq[p], "ACGT"[("ACGT".index(seq[p]) + 1) % 4], "het"),
               (p + 40, seq[p + 40:p + 43], seq[p + 40], "het")]
        g.diploid(seq, var, depth, frag_mean=300, frag_sd=40)


def complex_var(g):
    for _ in range(120):
        c, s, e, seq = g.window(pad=400)
        p = g.r.randrange(250, len(seq) - 300)
        alt = lambda b: g.r.choice([x for x in "ACGT" if x != b])
        kind = g.r.randrange(5)
        if kind == 0:  # SNP cluster
            var = [(p + d, seq[p + d], alt(seq[p + d]), g.r.choice(["het", "hom"])) for d in (0, 3, 7, 8, 15)]
        elif kind == 1:  # MNP
            var = [(p, seq[p:p + 3], "".join(alt(b) for b in seq[p:p + 3]), "het")]
        elif kind == 2:  # two alt alleles at one site
            a1 = alt(seq[p])
            a2 = g.r.choice([x for x in "ACGT" if x not in (seq[p], a1)])
            var = [(p, seq[p], a1, "het"), (p, seq[p], a2, "het2")]
        elif kind == 3:  # deletion next to a SNP
            var = [(p, seq[p:p + 6], seq[p], "het"), (p + 8, seq[p + 8], alt(seq[p + 8]), "het")]
        else:  # STR expansion
            unit = seq[p:p + g.r.choice([2, 3, 4])]
            var = [(p, seq[p], seq[p] + unit * g.r.randrange(2, 8), g.r.choice(["het", "hom"]))]
        g.diploid(seq, var, 45)


def overlap(g):
    for _ in range(80):
        c, s, e, ref = g.window()
        haps = g.plant(ref, indels=False)
        for _ in range(150):
            seq = g.r.choice(haps)
            kind = g.r.randrange(4)
            if kind == 0:
                g.pair_from(seq, frag=L)
            elif kind == 1:  # start-tie: mates share the 5' position on opposite strands
                p = g.r.randrange(0, len(seq) - L)
                g.add(g.errors(seq[p:p + L]), g.qual(L), g.errors(rc(seq[p:p + L])), g.qual(L))
            elif kind == 2:
                g.pair_from(seq, frag=g.r.randrange(L - 20, L))
            else:
                g.pair_from(seq, frag=g.r.randrange(L + 1, L + 40))


def coord_wrap(g):
    edge = [32766, 32767, 32768, 32769, 65535, 65536, 65537, 98303, 98304, 100000]
    for k in range(600):
        c, s, e, seq = g.window()
        haps = g.plant(seq, indels=False)
        if k % 6 == 0:  # some wrapped duplicate sets sit in called regions
            g.cover(haps)
        g.pair_from(g.r.choice(haps))
        name, r1, q1, r2, q2 = g.pairs.pop()
        lane = int(name.split(":")[3])
        x, y = g.r.choice(edge), g.r.choice(edge)
        tile = g.r.choice([1101, 2678, 32767, 32768, 40000])
        g.add(r1, q1, r2, q2, f"{g.names}:{lane}:{tile}:{x}:{y}")
        # Copies whose distance is small only before or only after the 16-bit wrap.
        for dx, dy in ((1, 0), (0, 1), (65536, 0), (0, 65536), (65500, 0), (40, 40)):
            if g.r.random() < 0.5:
                g.add(r1, g.qual(L), r2, g.qual(L), f"{g.names}:{lane}:{tile}:{x + dx}:{y + dy}")


SCENARIOS = {f.__name__: f for f in (homopolymer, long_indel, short_insert, n_bases, duplicates, chimera,
                                     unmapped, contig_edges, quality, read_length, pileup, complex_var,
                                     overlap, coord_wrap)}


def main():
    if len(sys.argv) != 7 or sys.argv[1] not in SCENARIOS:
        sys.exit(__doc__)
    scenario, seed, pairs, fasta, bed, outdir = sys.argv[1:]
    g = Gen(scenario, seed, Fasta(fasta), bed)
    g.background(int(pairs))
    SCENARIOS[scenario](g)
    n, w = g.write(outdir)
    print(f"{scenario}\tseed={seed}\tpairs={n}\twindows={w}")


if __name__ == "__main__":
    main()
