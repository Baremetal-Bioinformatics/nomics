#!/usr/bin/env python3
"""Summarize a Spark bake-off campaign (scripts/bakeoff/spark_campaign.sh).

For each run in CAMPAIGN/runs.tsv, over its [start, end] window in
CAMPAIGN/samples.tsv (system-wide, 5 s):
  wall      end - start
  cpu_s     busy CPU-seconds (all cores, /proc/stat user+nice+system+irq+softirq+steal)
  cores     cpu_s / wall (mean busy cores of 20)
  peak_gb   peak used memory (MemTotal - MemAvailable) minus the idle level just
            before the run, i.e. what the run added (page cache excluded)
  gpu_pct   mean GPU utilization while the run was active
Usage: summarize.py CAMPAIGN_DIR [--markdown]
"""
import csv
import os
import sys

CLK_TCK = os.sysconf("SC_CLK_TCK") if hasattr(os, "sysconf") else 100


def load_samples(path):
    rows = []
    with open(path) as f:
        r = csv.DictReader(f, delimiter="\t")
        for x in r:
            gpu = x["gpu_util"]
            rows.append((float(x["epoch"]), int(x["cpu_busy"]), int(x["cpu_total"]),
                         int(x["mem_used_kb"]), None if gpu in ("NA", "[N/A]", "") else float(gpu)))
    return rows


def fmt_wall(s):
    s = int(round(s))
    h, m, sec = s // 3600, s % 3600 // 60, s % 60
    return f"{h} h {m:02d} min" if h else f"{m} min {sec:02d} s"


def main():
    camp = sys.argv[1]
    md = "--markdown" in sys.argv
    samples = load_samples(os.path.join(camp, "samples.tsv"))
    out = []
    with open(os.path.join(camp, "runs.tsv")) as f:
        for run in csv.DictReader(f, delimiter="\t"):
            s, e = float(run["start"]), float(run["end"])
            win = [x for x in samples if s <= x[0] <= e]
            before = [x for x in samples if s - 60 <= x[0] < s]
            if len(win) < 2:
                out.append((run, e - s, None, None, None, None))
                continue
            cpu_s = (win[-1][1] - win[0][1]) / CLK_TCK
            span = win[-1][0] - win[0][0]
            idle_mem = min((x[3] for x in before), default=win[0][3])
            peak = (max(x[3] for x in win) - idle_mem) / 1048576
            gpus = [x[4] for x in win if x[4] is not None]
            gpu = sum(gpus) / len(gpus) if gpus else None
            out.append((run, e - s, cpu_s * (e - s) / span if span else cpu_s, peak, gpu, None))
    hdr = ["tool", "tier", "rep", "status", "wall", "wall_s", "cpu_s", "cores", "peak_gb", "gpu_pct"]
    lines = []
    for run, wall, cpu_s, peak, gpu, _ in out:
        lines.append([run["tool"], run["tier"], run["rep"], run["status"], fmt_wall(wall), f"{wall:.0f}",
                      "" if cpu_s is None else f"{cpu_s:.0f}",
                      "" if cpu_s is None else f"{cpu_s / wall:.1f}",
                      "" if peak is None else f"{peak:.1f}",
                      "" if gpu is None else f"{gpu:.0f}"])
    if md:
        print("| " + " | ".join(hdr) + " |")
        print("|" + "---|" * len(hdr))
        for l in lines:
            print("| " + " | ".join(l) + " |")
    else:
        print("\t".join(hdr))
        for l in lines:
            print("\t".join(l))


if __name__ == "__main__":
    main()
