#!/usr/bin/env python3
"""Is a sarek r2 the same as r1 in content? Reads reference/compare_runs.sh output.

    scripts/goldens/compare_summary.py [--scattered] compare_r1_r2.tsv

A file matches when its bytes match, or (BAM, VCF) its records and normalized header match,
or (recal table) its body matches. `mapped/*.sorted.bam` is a side output whose @RG suffixes
and tie order vary run to run (docs/REFERENCE-VARIABILITY.md); it is listed, not judged.
With --scattered (sarek ran BQSR per interval group and merged), `recalibrated/*.recal.bam` is
listed, not judged either: samtools merge gives inputs 2..n time-seeded @RG/@PG ID suffixes that
every record carries, so two runs differ there by construction. nomics is gated on it against r1.
Exits 1 when any judged file differs.
"""
import sys
from collections import defaultdict

args = sys.argv[1:]
scattered = "--scattered" in args
args = [a for a in args if a != "--scattered"]
checks = defaultdict(dict)
for line in open(args[0]):
    parts = line.rstrip("\n").split("\t")
    if len(parts) < 4:
        continue
    _stage, path, check, result = parts[:4]
    checks[path][check] = result

if not checks:
    print("no files compared")
    sys.exit(1)

bad = 0
for path in sorted(checks):
    c = checks[path]
    if c.get("bytes") == "PASS":
        verdict = "same bytes"
    elif path.endswith(".recal.table"):
        verdict = "same body" if c.get("body") == "PASS" else "DIFFERS"
    elif c.get("records") == "PASS" and c.get("header") == "PASS":
        verdict = "same records and header"
    else:
        verdict = "DIFFERS"
    if verdict == "DIFFERS" and path.endswith(".sorted.bam"):
        verdict = "differs (side output, not judged)"
    elif verdict == "DIFFERS" and scattered and path.endswith(".recal.bam"):
        verdict = "differs (scattered: merge ID suffixes, not judged)"
    elif verdict == "DIFFERS":
        bad += 1
    print(f"{verdict}\t{path}")

print(f"r1 vs r2: {len(checks) - bad} of {len(checks)} files match" + ("" if not bad else f", {bad} differ"))
sys.exit(1 if bad else 0)
