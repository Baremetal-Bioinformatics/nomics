#!/usr/bin/env python3
"""System-wide resource sampler for the Spark bake-off campaign.

Every INTERVAL seconds appends one TSV line:
  epoch  cpu_busy_jiffies  cpu_total_jiffies  mem_used_kb  gpu_util_pct  gpu_mem_mb
The same instrument covers sarek (docker + JVM), nomics and Parabricks (docker),
so CPU-seconds, peak memory and GPU use are comparable across tools.
Usage: sampler.py OUT.tsv [interval_s]
"""
import subprocess
import sys
import time


def cpu():
    with open("/proc/stat") as f:
        v = [int(x) for x in f.readline().split()[1:]]
    idle = v[3] + v[4]  # idle + iowait
    return sum(v) - idle, sum(v)


def mem_used_kb():
    m = {}
    with open("/proc/meminfo") as f:
        for line in f:
            k, v = line.split(":", 1)
            m[k] = int(v.split()[0])
    return m["MemTotal"] - m["MemAvailable"]


def gpu():
    try:
        out = subprocess.run(
            ["nvidia-smi", "--query-gpu=utilization.gpu,memory.used", "--format=csv,noheader,nounits"],
            capture_output=True, text=True, timeout=10,
        ).stdout.strip().splitlines()[0]
        u, m = (x.strip() for x in out.split(","))
        return u, m
    except Exception:
        return "NA", "NA"


def main():
    path = sys.argv[1]
    interval = float(sys.argv[2]) if len(sys.argv) > 2 else 5.0
    with open(path, "a", buffering=1) as f:
        if f.tell() == 0:
            f.write("epoch\tcpu_busy\tcpu_total\tmem_used_kb\tgpu_util\tgpu_mem_mb\n")
        while True:
            busy, total = cpu()
            u, m = gpu()
            f.write(f"{time.time():.1f}\t{busy}\t{total}\t{mem_used_kb()}\t{u}\t{m}\n")
            time.sleep(interval)


if __name__ == "__main__":
    main()
