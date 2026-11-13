#!/usr/bin/env python3
"""Record a sarek WGS run's scatter as nomics `--scatter-order` input.

    scripts/scatter_order.py <sarek run dir> <sample> > reference/<run>.scatter-order.txt

Writes `merge-seed`, then `gather <group>` lines (GatherBQSRReports --input
order) and `merge <group>` lines (samtools merge input order of the recal
BAMs): all three follow how that run happened (docs/REFERENCE-VARIABILITY.md).

The seed is samtools merge's `hts_srand48(time(NULL))`: found by searching
the seconds around the recal MERGE_BAM task's start for the one whose first
lrand48() draw gives the suffix of the second @RG line of the merged BAM.
"""
import calendar
import glob
import gzip
import os
import re
import struct
import sys
import time

run, sample = sys.argv[1], sys.argv[2]
trace = sorted(glob.glob(f"{run}/results/pipeline_info/execution_trace*.txt"))[-1]
rows = [line.rstrip("\n").split("\t") for line in open(trace)]
head = rows[0]
col = {name: i for i, name in enumerate(head)}


def task_dir(h):
    return glob.glob(f"{run}/work/{h}*")[0]


def command(h):
    return open(f"{task_dir(h)}/.command.sh").read()


gather = next(r for r in rows[1:] if ":GATK4_GATHERBQSRREPORTS (" in r[col["name"]])
gather_order = re.findall(rf"--input {re.escape(sample)}_(\S+?)\.recal\.table", command(gather[col["hash"]]))
merge = next(
    r
    for r in rows[1:]
    if ":MERGE_BAM (" in r[col["name"]] and ".recal.bam" in command(r[col["hash"]])
)
merge_order = re.findall(rf"\d+/{re.escape(sample)}_(\S+?)\.recal\.bam", command(merge[col["hash"]]))

bam = f"{run}/results/preprocessing/recalibrated/{sample}/{sample}.recal.bam"
with gzip.open(bam, "rb") as f:  # BGZF is gzip-compatible; the header leads the stream
    if f.read(4) != b"BAM\x01":
        sys.exit(f"{bam}: not a BAM")
    (l_text,) = struct.unpack("<i", f.read(4))
    text = f.read(l_text).decode()
rg_ids = re.findall(r"^@RG\tID:(\S+)", text, re.M)
suffix = int(rg_ids[1].rsplit("-", 1)[1], 16)

submit = merge[col["submit"]]  # e.g. 2026-10-10 04:41:58.022 (local time of the run host)
local = time.strptime(submit.split(".")[0], "%Y-%m-%d %H:%M:%S")
approx = calendar.timegm(local)
A, C, M = 0x5DEECE66D, 0xB, 1 << 48
seed = None
for s in range(approx - 26 * 3600, approx + 26 * 3600):  # host timezone unknown: search ±26 h
    x = (((s & 0xFFFFFFFF) << 16) | 0x330E) % M
    x = (A * x + C) % M
    if x >> 17 == suffix:
        seed = s
        break
if seed is None:
    sys.exit("merge seed not found")

print(f"# sarek {os.path.basename(os.path.normpath(run))}: GatherBQSRReports --input order and samtools merge (recal) input order")
print(f"merge-seed {seed}")
for g in gather_order:
    print(f"gather {g}")
for g in merge_order:
    print(f"merge {g}")
