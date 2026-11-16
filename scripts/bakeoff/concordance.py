#!/usr/bin/env python3
"""VCF concordance for the Parabricks bake-off.

For each candidate VCF:
- vs the sarek golden: identical record lines (ignoring ID and INFO, which
  depend on dbSNP annotation), identical (site, alleles, GT), and sites only
  in one side.
- vs GIAB v4.2.1 inside (truth bed AND targets): TP/FP/FN on exact
  (chrom, pos, ref, alt, GT-class) matches. No allele normalization, so this
  is a coarse, comparable proxy, not a hap.py replacement; it is applied
  identically to the golden (== nomics, which is identical to the golden).
"""
import argparse
import bisect
import gzip
from collections import defaultdict


def read_bed(path):
    iv = defaultdict(list)
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            c, s, e = line.split("\t")[:3]
            iv[c].append((int(s), int(e)))
    out = {}
    for c, v in iv.items():
        v.sort()
        merged = []
        for s, e in v:
            if merged and s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        out[c] = ([m[0] for m in merged], [m[1] for m in merged])
    return out


def inside(bed, c, pos0):
    if c not in bed:
        return False
    starts, ends = bed[c]
    i = bisect.bisect_right(starts, pos0) - 1
    return i >= 0 and pos0 < ends[i]


def gt_class(sample):
    gt = sample.split(":")[0].replace("|", "/")
    a = [x for x in gt.split("/") if x not in (".",)]
    if not a or all(x == "0" for x in a):
        return None
    return "hom" if len(set(a)) == 1 else "het"


def read_vcf(path):
    lines, calls = {}, {}
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            t = line.rstrip("\n").split("\t")
            key = (t[0], int(t[1]), t[3], t[4])
            # ignore ID (dbSNP) and INFO (DB flag and annotation formatting)
            lines[key] = "\t".join(t[:2] + t[3:7] + t[8:])
            calls[key] = gt_class(t[9]) if len(t) > 9 else None
    return lines, calls


def vs_truth(calls, truth, region):
    def keep(k):
        return inside(region[0], k[0], k[1] - 1) and inside(region[1], k[0], k[1] - 1)

    q = {k: g for k, g in calls.items() if g and keep(k)}
    t = {k: g for k, g in truth.items() if g and keep(k)}
    tp = sum(1 for k, g in q.items() if t.get(k) == g)
    fp, fn = len(q) - tp, len(t) - tp
    prec = tp / (tp + fp) if tp + fp else 0.0
    rec = tp / (tp + fn) if tp + fn else 0.0
    f1 = 2 * prec * rec / (prec + rec) if prec + rec else 0.0
    return tp, fp, fn, prec, rec, f1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--golden", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--truth-bed", required=True)
    ap.add_argument("--targets", required=True)
    ap.add_argument("--vcf", action="append", required=True, help="name=path")
    a = ap.parse_args()
    region = (read_bed(a.truth_bed), read_bed(a.targets))
    _, truth = read_vcf(a.truth)
    g_lines, g_calls = read_vcf(a.golden)
    print(f"golden (sarek == nomics): {len(g_lines)} records")
    tp, fp, fn, p, r, f1 = vs_truth(g_calls, truth, region)
    print(f"  vs GIAB (truth bed AND targets): TP {tp} FP {fp} FN {fn}  precision {p:.5f} recall {r:.5f} F1 {f1:.5f}")
    for spec in a.vcf:
        name, path = spec.split("=", 1)
        lines, calls = read_vcf(path)
        both = set(lines) & set(g_lines)
        same_line = sum(1 for k in both if lines[k] == g_lines[k])
        same_gt = sum(1 for k in both if calls[k] == g_calls[k])
        print(f"{name}: {len(lines)} records")
        print(f"  vs golden: identical records (excl. ID/INFO) {same_line} of {len(g_lines)} "
              f"({100 * same_line / len(g_lines):.2f}%); same site+alleles {len(both)}, same genotype {same_gt}; "
              f"only in {name} {len(set(lines) - set(g_lines))}, only in golden {len(set(g_lines) - set(lines))}")
        tp, fp, fn, p, r, f1 = vs_truth(calls, truth, region)
        print(f"  vs GIAB (truth bed AND targets): TP {tp} FP {fp} FN {fn}  precision {p:.5f} recall {r:.5f} F1 {f1:.5f}")


if __name__ == "__main__":
    main()
