#!/usr/bin/env python3
"""Replay a golden nf-core/sarek work dir task by task through the nomics shims.

No Nextflow and no containers. For every task (a work dir with .command.sh):

1. Recreate the task's inputs in scratch: every symlink (and subdirectory of
   symlinks) is copied, with targets on the Spark remapped to local fixtures
   (REMAP below). Upstream task outputs resolve to the *golden* task's files,
   so every task is gated independently of the others.
2. Run the task's script (the part before Nextflow's "capture process
   environment" tail, which only records tool versions) with
   `bash -C -e -u -o pipefail` and a shim bin dir (fastp, bwa, gatk, samtools
   -> nomics) first on PATH. The real tools are whatever else is on PATH.
3. Compare every output file the golden task wrote with ours, using the
   docs/GATES.md comparators:
     *.fastp.fastq.gz  nomics dev gate-chunks (decompressed content, per chunk)
     *.bam             nomics dev gate-bam (records + header except @PG)
     *.table           nomics dev gate-report (GATKReport bytes, then cells)
     *.vcf.gz          nomics dev gate-vcf (records + header except command lines)
     *.metrics         bytes except lines 2 and 4
     *.bai / *.tbi     handled tasks: htsjdk layout (nomics dev htsjdk-index
                       reproduces the golden index from the golden data file,
                       and ours from our data file); otherwise bytes
     anything else     bytes
   Outputs that GATES.md does not gate (fastp JSON/HTML/log) are listed as such.

The shim's stderr line (`nomics-shim: handled ...` / `delegated ...`) gives
the mode. Exactness bar: every *handled* task must PASS. Delegated and native
tasks are run when the real tool is installed and reported for information
("needs real tool" otherwise).

Usage:
  scripts/nfcore/replay_tasks.py [--golden DIR] [--scratch DIR] [--nomics BIN]
                                 [--only PATTERN ...] [--skip PATTERN ...]
"""

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

HOME = Path.home()
# Home directory of the machine the golden sarek run happened on.
GOLDEN_HOME = os.environ.get("NOMICS_GOLDEN_HOME", os.path.expanduser("~"))
REPO = Path(__file__).resolve().parents[2]

# Spark path prefix -> local path (first match wins). {golden} is --golden.
REMAP = [
    (f"{GOLDEN_HOME}/omics/runs/tier0-wes-20261008-r1/work/", "{golden}/work/"),
    (f"{GOLDEN_HOME}/omics/ref/GATK.GRCh38/Sequence/WholeGenomeFasta/",
     str(HOME / "omics-fixtures/gvcf-tier0/ref") + "/"),
    (f"{GOLDEN_HOME}/omics/ref/GATK.GRCh38/Annotation/GATKBundle/beta/",
     str(HOME / "omics-fixtures/gvcf-tier0/ref/beta") + "/"),
    (f"{GOLDEN_HOME}/omics/ref/GATK.GRCh38/Annotation/GATKBundle/",
     str(HOME / "omics-fixtures/gvcf-tier0/ref") + "/"),
    (f"{GOLDEN_HOME}/omics/ref/GATK.GRCh38/Sequence/BWAIndex",
     str(HOME / "omics-fixtures/gvcf-tier0/ref/BWAIndex")),
    (f"{GOLDEN_HOME}/omics/data/Agilent_v7.bed",
     str(HOME / "omics-fixtures/gvcf-tier0/Agilent_v7.bed")),
]

SHIM_TOOLS = ["fastp", "bwa", "gatk", "samtools"]
NOT_GATED = (".fastp.json", ".fastp.html", ".fastp.log")
SKIP_FILES = {".exitcode", "versions.yml"}


def remap(target, golden):
    for old, new in REMAP:
        if target.startswith(old):
            return new.format(golden=golden) + target[len(old):]
    return target


def process_name(task):
    for line in (task / ".command.run").read_text().splitlines()[:5]:
        m = re.match(r"### name: '(.*)'", line)
        if m:
            return m.group(1)
    return "?"


def short_process(name):
    # NFCORE_SAREK:SAREK:...:GATK4_MARKDUPLICATES (HG002_WES) -> GATK4_MARKDUPLICATES
    return name.split(" (")[0].split(":")[-1]


def stage_inputs(src, dst, golden, missing):
    """Recreate src's symlinks and subdirectories in dst, remapped."""
    dst.mkdir(parents=True, exist_ok=True)
    for entry in sorted(src.iterdir()):
        out = dst / entry.name
        if entry.is_symlink():
            target = remap(os.readlink(entry), golden)
            os.symlink(target, out)
            resolved = target if os.path.isabs(target) else str(dst / target)
            if not os.path.exists(resolved):
                missing.append(f"{entry.name} -> {target}")
        elif entry.is_dir():
            stage_inputs(entry, out, golden, missing)


def golden_outputs(task):
    """Regular files the golden task wrote (its outputs)."""
    outs = []
    for p in sorted(task.rglob("*")):
        if p.is_symlink() or not p.is_file():
            continue
        rel = p.relative_to(task)
        if rel.name.startswith(".command") or rel.name in SKIP_FILES:
            continue
        outs.append(rel)
    return outs


def run_gate(nomics, sub, ours, golden):
    r = subprocess.run([nomics, "dev", sub, "--ours", str(ours), "--golden", str(golden)],
                       capture_output=True, text=True)
    tail = (r.stdout + r.stderr).strip().splitlines()
    return r.returncode == 0, (tail[-1] if tail else "")


def metrics_equal(a, b):
    la, lb = a.read_text().splitlines(), b.read_text().splitlines()
    keep = lambda ls: [l for i, l in enumerate(ls, 1) if i not in (2, 4)]
    return keep(la) == keep(lb)


def htsjdk_index(nomics, data, out):
    r = subprocess.run([nomics, "dev", "htsjdk-index", "--input", str(data), "--output", str(out)],
                       capture_output=True, text=True)
    return r.returncode == 0


def index_data_file(index_path):
    name = index_path.name
    if name.endswith(".tbi"):
        return index_path.with_name(name[:-4])
    if name.endswith(".bam.bai"):
        return index_path.with_name(name[:-4])
    if name.endswith(".bai"):
        return index_path.with_name(name[:-4] + ".bam")
    return None


def compare(nomics, rel, ours_dir, golden_dir, handled, scratch_tmp):
    """Return (status, note) for one output; status in PASS/FAIL/MISSING/NOT GATED."""
    ours, golden = ours_dir / rel, golden_dir / rel
    name = rel.name
    if name.endswith(NOT_GATED):
        return ("NOT GATED", "QC side output" + ("" if ours.exists() else ", not written"))
    if not ours.exists():
        return ("MISSING", "")
    if name.endswith(".sorted.bam"):
        # GATES.md: samtools merge draws random @RG/@PG ID suffixes; not gated.
        ok, note = run_gate(nomics, "gate-bam", ours, golden)
        return ("NOT GATED", "samtools merge random ID suffixes" + (", records identical" if ok else ""))
    if name.endswith(".bam"):
        ok, note = run_gate(nomics, "gate-bam", ours, golden)
        return ("PASS" if ok else "FAIL", "gate-bam" + ("" if ok else f": {note}"))
    if name.endswith(".table"):
        ok, note = run_gate(nomics, "gate-report", ours, golden)
        return ("PASS" if ok else "FAIL", "gate-report" + ("" if ok else f": {note}"))
    if name.endswith(".vcf.gz"):
        ok, note = run_gate(nomics, "gate-vcf", ours, golden)
        return ("PASS" if ok else "FAIL", "gate-vcf" + ("" if ok else f": {note}"))
    if name.endswith(".metrics"):
        ok = metrics_equal(ours, golden)
        return ("PASS" if ok else "FAIL", "bytes except lines 2, 4")
    if handled and name.endswith((".bai", ".tbi")):
        data_ours, data_golden = index_data_file(ours), index_data_file(golden)
        if data_ours is None or not data_ours.exists() or not data_golden.exists():
            return ("FAIL", "index without its data file")
        if ours.read_bytes() == golden.read_bytes():
            return ("PASS", "bytes")
        g = scratch_tmp / ("golden-layout-" + name)
        o = scratch_tmp / ("ours-layout-" + name)
        layout_ok = htsjdk_index(nomics, data_golden, g) and g.read_bytes() == golden.read_bytes()
        ours_ok = htsjdk_index(nomics, data_ours, o) and o.read_bytes() == ours.read_bytes()
        if layout_ok and ours_ok:
            return ("PASS", "htsjdk layout (golden index reproduced from golden data)")
        return ("FAIL", f"htsjdk layout: golden {'ok' if layout_ok else 'NOT reproduced'}, "
                        f"ours {'consistent' if ours_ok else 'inconsistent'}")
    ok = ours.read_bytes() == golden.read_bytes()
    return ("PASS" if ok else "FAIL", "bytes")


def synthesize_gvcf_task(fixture, root, golden):
    """A GATK4_HAPLOTYPECALLER -ERC GVCF task (sarek --joint_germline) built
    from the gvcf-tier0 fixture: its recorded .command.sh, recal.bam and
    golden g.vcf.gz. The tier0 r1 run itself has no GVCF task."""
    task = root / "work" / "gv" / "cf-tier0-joint-germline-hc"
    if task.exists():
        shutil.rmtree(task)
    task.mkdir(parents=True)
    (task / ".command.sh").write_text((fixture / "hc-command.txt").read_text()
                                      .split("chr1_12146-12310.bed\n")[0])
    (task / ".command.run").write_text(
        "#!/bin/bash\n### ---\n### name: 'NFCORE_SAREK:SAREK:BAM_VARIANT_CALLING_GERMLINE_ALL:"
        "BAM_VARIANT_CALLING_HAPLOTYPECALLER:GATK4_HAPLOTYPECALLER (HG002_WES, -ERC GVCF)'\n")
    ref = fixture / "ref"
    links = {
        "HG002_WES.recal.bam": fixture / "recal/HG002_WES.recal.bam",
        "HG002_WES.recal.bam.bai": fixture / "recal/HG002_WES.recal.bam.bai",
        "Homo_sapiens_assembly38.fasta": ref / "Homo_sapiens_assembly38.fasta",
        "Homo_sapiens_assembly38.fasta.fai": ref / "Homo_sapiens_assembly38.fasta.fai",
        "Homo_sapiens_assembly38.dict": ref / "Homo_sapiens_assembly38.dict",
        "dbsnp_146.hg38.vcf.gz": ref / "dbsnp_146.hg38.vcf.gz",
        "dbsnp_146.hg38.vcf.gz.tbi": ref / "dbsnp_146.hg38.vcf.gz.tbi",
        # The same CREATE_INDEX_BED group as the tier0 r1 run.
        "chr1_12146-12310.bed": golden / "work/af/0060662473c92d42ad990b98e1b8c9/chr1_12146-12310.bed",
    }
    for name, target in links.items():
        os.symlink(target, task / name)
    for name in ["HG002_WES.haplotypecaller.g.vcf.gz", "HG002_WES.haplotypecaller.g.vcf.gz.tbi"]:
        shutil.copyfile(fixture / "golden" / name, task / name)
    return task


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--golden", type=Path, default=HOME / "omics-fixtures/sarek-tier0/r1")
    ap.add_argument("--scratch", type=Path,
                    default=Path(os.environ.get("TMPDIR", "/tmp")) / "nomics-replay")
    ap.add_argument("--nomics", type=Path, default=REPO / "target/release/nomics")
    ap.add_argument("--only", action="append", default=[],
                    help="regex on 'hash process'; repeatable")
    ap.add_argument("--skip", action="append", default=[], help="regex; repeatable")
    ap.add_argument("--json", type=Path, help="write the results here")
    ap.add_argument("--gvcf-fixture", type=Path, default=HOME / "omics-fixtures/gvcf-tier0",
                    help="also replay an -ERC GVCF HaplotypeCaller task built from this "
                         "fixture (hc-command.txt, recal/, golden/); '' to skip")
    args = ap.parse_args()

    golden = args.golden.resolve()
    nomics = str(args.nomics.resolve())
    if not os.access(nomics, os.X_OK):
        sys.exit(f"no nomics binary at {nomics} (cargo build --release -p nomics)")
    shim = args.scratch / "shim-bin"
    if shim.exists():
        shutil.rmtree(shim)
    shim.mkdir(parents=True)
    for t in SHIM_TOOLS:
        os.symlink(nomics, shim / t)
    bash = shutil.which("bash") or "/bin/bash"

    tasks = sorted(p.parent for p in (golden / "work").glob("*/*/.command.sh"))
    if str(args.gvcf_fixture) not in ("", ".") and (args.gvcf_fixture / "hc-command.txt").exists():
        tasks.append(synthesize_gvcf_task(args.gvcf_fixture, args.scratch / "synthetic", golden))
    results = []
    for task in tasks:
        tid = f"{task.parent.name}/{task.name[:6]}"
        proc = process_name(task)
        key = f"{task.parent.name}/{task.name} {proc}"
        if args.only and not any(re.search(p, key) for p in args.only):
            continue
        if any(re.search(p, key) for p in args.skip):
            continue
        work = args.scratch / "tasks" / task.parent.name / task.name
        if work.exists():
            shutil.rmtree(work)
        missing = []
        stage_inputs(task, work, str(golden), missing)
        script = (task / ".command.sh").read_text().split("# capture process environment")[0]
        (work / ".command.main.sh").write_text(script)
        row = {"task": tid, "process": short_process(proc), "mode": "", "outputs": [],
               "verdict": "", "seconds": 0.0, "shim": []}
        outs = golden_outputs(task)
        if missing:
            row.update(mode="-", verdict="INPUTS MISSING",
                       note="; ".join(missing))
            results.append(row)
            print(f"{tid} {row['process']}: inputs missing: {row['note']}", flush=True)
            continue
        env = dict(os.environ)
        env["PATH"] = f"{shim}:{env.get('PATH', '')}"
        env["NXF_TASK_WORKDIR"] = str(work)
        t0 = time.time()
        with open(work / ".command.out", "wb") as so, open(work / ".command.err", "wb") as se:
            rc = subprocess.run([bash, "-C", "-e", "-u", "-o", "pipefail", ".command.main.sh"],
                                cwd=work, env=env, stdout=so, stderr=se).returncode
        row["seconds"] = round(time.time() - t0, 1)
        err = (work / ".command.err").read_text(errors="replace")
        # Concurrent processes in one pipe (bwa | samtools sort) can interleave
        # their stderr mid-line, so find every shim decision anywhere in the text.
        decisions = re.findall(r"nomics-shim: (handled|delegated) ([^\n]*?)(?=nomics-shim:|\n|$)", err)
        row["shim"] = [f"nomics-shim: {m} {rest.strip()}" for m, rest in decisions]
        modes = {m for m, _ in decisions}
        # A task whose output passes through any handled call is held to the
        # exactness bar (e.g. bwa handled | samtools sort delegated).
        handled = "handled" in modes
        row["mode"] = ("handled+deleg" if modes == {"handled", "delegated"} else
                       "handled" if handled else "delegated" if "delegated" in modes
                       else "native")
        no_real = ("cannot delegate" in err) or rc == 127 or "command not found" in err
        tmp = work / ".gate-tmp"
        tmp.mkdir(exist_ok=True)
        chunk_files = [o for o in outs if o.name.endswith(".fastp.fastq.gz")]
        if chunk_files:
            ok, note = run_gate(nomics, "gate-chunks", work, task)
            row["outputs"].append({"file": f"{len(chunk_files)} x *.fastp.fastq.gz",
                                   "status": "PASS" if ok else "FAIL",
                                   "note": "gate-chunks" + ("" if ok else f": {note}")})
        for rel in outs:
            if rel in chunk_files:
                continue
            status, note = compare(nomics, rel, work, task, handled, tmp)
            row["outputs"].append({"file": str(rel), "status": status, "note": note})
        gated = [o for o in row["outputs"] if o["status"] != "NOT GATED"]
        all_pass = gated and all(o["status"] == "PASS" for o in gated)
        if rc != 0 and no_real:
            row["verdict"] = "NEEDS REAL TOOL"
        elif rc != 0:
            row["verdict"] = f"FAIL (exit {rc})"
        else:
            row["verdict"] = ("NOT GATED" if not gated else "PASS" if all_pass else "FAIL")
        results.append(row)
        print(f"{tid} {row['process']}: {row['mode']} {row['verdict']} ({row['seconds']}s)",
              flush=True)

    # Table
    print()
    hdr = f"{'task':<10} {'process':<36} {'mode':<10} {'outputs (status)':<60} verdict"
    print(hdr)
    print("-" * len(hdr))
    for r in results:
        outs = ", ".join(f"{o['file']} {o['status']}" for o in r["outputs"]) or "-"
        print(f"{r['task']:<10} {r['process']:<36} {r['mode']:<10} {outs:<60} {r['verdict']}")
    handled_fail = [r for r in results if r["mode"].startswith("handled") and r["verdict"] != "PASS"]
    print()
    print(f"{len(results)} tasks; handled {sum(r['mode'].startswith('handled') for r in results)}, "
          f"handled not PASS: {len(handled_fail)}")
    if args.json:
        args.json.write_text(json.dumps(results, indent=1))
    sys.exit(1 if handled_fail else 0)


if __name__ == "__main__":
    main()
