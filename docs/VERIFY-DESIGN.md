# `nomics verify`: design note (for review before the build)

Status: built, 2026-10-30. This note is the design as reviewed; [VERIFY.md](VERIFY.md)
describes what shipped. Where they differ, VERIFY.md is right: `--reference-config` became
`--also-sarek-default` alone, the sarek source is a pinned git clone (commit checked) rather
than a tarball, and the Nextflow image is pinned by tag with its digest recorded in the report.

`nomics verify` answers one question: **is nomics identical to sarek on my data, on my
machine?** It takes a deterministic slice of the user's reads, runs pinned sarek 3.10.0 and
nomics on that slice, gates every artifact with the `docs/GATES.md` normalizations, and
writes a report anyone can read. The user does not have to trust us: sarek runs on their
machine, from pinned containers.

## 1. Commands

```bash
# Exome: the common case.
nomics verify --fastq1 R1.fastq.gz --fastq2 R2.fastq.gz \
  --reference /path/to/GATK.GRCh38 --targets targets.bed --outdir verify/

# Whole genome: no --targets; reads are sliced and calling is limited to one contig.
nomics verify --fastq1 R1.fastq.gz --fastq2 R2.fastq.gz \
  --reference /path/to/GATK.GRCh38 --wgs --outdir verify/

# No sarek: check a nomics run on the public fixture against a shipped digest manifest.
nomics verify --against reference/tier0-wes-20261008-r1.digest.tsv \
  --fastq1 tier0_1.fastq.gz --fastq2 tier0_2.fastq.gz \
  --reference /path/to/GATK.GRCh38 --targets Agilent_v7.bed --outdir verify/

# Pull the pinned containers and sarek source once, then run with no network at all.
nomics verify fetch
nomics verify --offline ...
```

### Flags

| Flag | Default | Meaning |
|---|---|---|
| `--fastq1`, `--fastq2` | required | One paired FASTQ (gzip or BGZF). Multi-lane input is a later addition; v1 refuses it. |
| `--reference DIR` | required | iGenomes `GATK.GRCh38` layout. Refused if its `.dict` differs from the pinned GRCh38 dictionary. |
| `--targets BED` | required unless `--wgs` | Exome targets, as sarek's `--wes --intervals`. |
| `--wgs` | off | Whole-genome mode: no targets; calling restricted by `--region`. |
| `--region CONTIG` | `chr20` with `--wgs` | Limits HaplotypeCaller intervals (sarek `--intervals`) so a genome fits the time budget. Upstream stages still see every sliced read. |
| `--slice head:N` | `head:4000000` (genome) | Take the first N read pairs, byte for byte. See section 2. |
| `--slice contig:NAME` | `contig:chr20` (exome) | Pairs whose nomics alignment touches NAME (`nomics dev slice`). Slower and version-dependent; see question 3. |
| `--outdir DIR` | required | New or empty directory. Layout in section 5. |
| `--threads N` | all cores | Both sides get the same N; sarek via Nextflow `cpus`. |
| `--device auto\|cpu\|cuda` | `auto` | nomics side only. `auto` uses CUDA when the build has it and a GPU answers `selftest`. |
| `--cache-dir DIR` | `~/.cache/nomics` | nomics' dense-SA and known-sites caches (built on first use, about 6 min, reported separately). |
| `--reference-config pinned\|sarek-default` | `pinned` | `pinned`: arm64 containers + the local GATK image, or x86 with `cloud-hcjava.config` (Java PairHMM). See question 1. |
| `--also-sarek-default` | off | x86 only: also run sarek's stock configuration (Intel GKL) and report its VCF difference as **EXPECTED**, never as PASS or FAIL. |
| `--against MANIFEST` | off | Second mode: no sarek. Requires the exact FASTQs the manifest names (md5-checked); runs nomics and compares digests. |
| `--offline` | off | Never touch the network; error if a pinned container or the sarek source is missing. |
| `--keep-work` | off | Keep sarek's `work/` and nomics' intermediates (otherwise removed after a PASS; always kept after FAIL/ERROR). |
| `--min-records N` | *provisional* `200` | VCF records required on the reference side for the result to mean anything; below it the verdict is ERROR, not PASS. |

`nomics verify fetch [--reference-config ...]` pulls every pinned container by digest,
builds the arm64 GATK image from the embedded `reference/containers/gatk4/Dockerfile`, and
downloads the sarek source tarball for commit `8ccac7ad37b05dd792447763bf9671b719824587`,
checking its SHA-256. A normal run does this on first use; `--offline` refuses to.

### Verdicts and exit codes

| Verdict | Exit | When |
|---|---|---|
| **PASS** | 0 | Every artifact in the descriptor matches, every count is non-zero, the reference VCF has at least `--min-records` records. |
| **FAIL** | 1 | Both sides ran to completion and at least one artifact differs or is missing. The report names each one and its first differing record. |
| **ERROR** | 2 | Anything that prevents a real comparison: missing Docker, a failed sarek task, a nomics error, refused input, empty outputs, a slice too small. Never reported as PASS. |

## 2. The slice

- **Method.** `head:N` takes the first N pairs of R1 and R2 in file order, copying records
  byte for byte, and writes them BGZF-compressed. It needs no alignment, so the slice is the
  same on every machine and every nomics version. The 1/12 genome tier was made the same way.
  BWA's batching and insert-size estimation see a valid input (a prefix of a real run).
- **Checks while slicing** (each an ERROR with the reason): R1/R2 names disagree at the same
  position, one file ends first, a malformed record, fewer than N pairs available (then the
  whole file is used and the report says so).
- **Slice digest.** SHA-256 over the decompressed R1 records then R2 records, so the digest
  does not depend on our compressor. Recorded in the report; the repeatability gate checks it.
- **Why the defaults are provisional.** 500k exome pairs is tier0 scale (tier0 has 489k
  pairs; sarek takes 2 min 15 s on it on the Spark), but tier0 is concentrated on chr20
  targets, while a head slice spreads over the whole exome at about 4x. The build measures
  record counts and wall times on HG002 and the Brief 1 samples and picks N so a 20-core
  machine finishes in under 30 minutes with enough calls to mean something.

## 3. What runs

1. **Preflight** (seconds): Docker reachable; reference `.dict` matches pinned GRCh38; BED
   parses and its contigs exist; outdir empty; disk space; CPU/GPU inventory.
2. **Slice** (section 2).
3. **sarek**, from the descriptor: pinned Nextflow 26.04.6 running inside its own container
   (`nextflow/nextflow` by digest, Docker socket mounted, paths mounted at identical
   locations), so Docker is the only prerequisite (question 2). Params as the goldens:
   `tools: haplotypecaller`, `skip_tools: haplotypecaller_filter`, `save_output_as_bam: true`,
   `save_mapped: true`, `--wes --intervals` or WGS with `--intervals <region>`. Containers run
   with `--pull never` after fetch. The trace must show every task `COMPLETED`, exit 0.
4. **Capture from the sarek run**: scatter order (gather and merge orders, merge seed: the
   logic of `scripts/scatter_order.py`, ported to Rust so verify needs no Python), the
   read-group `DS` and sample fields from md.bam's header, the fastp chunk files from
   `work/`.
5. **nomics run** on the same slice with `--write-intermediates`, `--scatter-order`,
   `--read-group-ds`, the same sample and prefix.
6. **Gate**: digest both sides with the descriptor's comparator per artifact. On a digest
   mismatch, run that artifact's detailed comparator (`compare_bams`, `compare_vcfs`,
   `compare_gatk_reports`, `compare_fastq_chunks`) to find the first differing record.
   Any artifact with count 0 on either side is ERROR ("empty output"), never PASS.
7. **Report** (section 6).

Privacy: no read data leaves the machine or enters the report. The only network traffic is
`fetch` (containers by digest, the sarek tarball). The report records file sizes, md5s and
basenames, not full paths.

## 4. The pipeline descriptor

verify is written around a descriptor so a later "port your own pipeline" harness can add
pipelines as data. sarek 3.10.0 germline is the only one now. It is a small line-oriented
text file embedded in the binary (`crates/nomics/descriptors/sarek-3.10.0-germline.tsv`), in
the same spirit as the digest manifests; no general descriptor language and no new
dependencies.

```text
# nomics pipeline descriptor v1
id              sarek-3.10.0-germline
pipeline        nf-core/sarek
version         3.10.0
commit          8ccac7ad37b05dd792447763bf9671b719824587
source_sha256   <tarball sha256>
nextflow        26.04.6  nextflow/nextflow@sha256:<digest>
config.x86      docker   reference/cloud-hcjava.config
config.arm64    docker,arm64  reference/containers/gatk4/Dockerfile
param           tools                 haplotypecaller
param           skip_tools            haplotypecaller_filter
param           save_output_as_bam    true
param           save_mapped           true
capture         scatter_order         trace:GATK4_GATHERBQSRREPORTS,MERGE_BAM
capture         read_group_ds         header:preprocessing/markduplicates/{prefix}/{prefix}.md.bam
# artifact      name          location                                                   comparator   normalization (GATES.md row)
artifact        chunk/*       work:*.{prefix}-{lane}_R?.fastp.fastq.gz                    fastq_chunk  decompressed content per chunk
artifact        mapped/*      published:preprocessing/mapped/{prefix}/{prefix}.????.bam   bam          records + header without @PG
artifact        md_bam        published:preprocessing/markduplicates/{prefix}/{prefix}.md.bam   bam   records + header without @PG
artifact        md_metrics    published:reports/markduplicates/{prefix}/{prefix}.md.metrics     md_metrics  bytes except lines 2 and 4
artifact        recal_table   published:preprocessing/recal_table/{prefix}/{prefix}.recal.table raw   bytes
artifact        recal_bam     published:preprocessing/recalibrated/{prefix}/{prefix}.recal.bam  bam   records + header without @PG
artifact        vcf           published:variant_calling/haplotypecaller/{prefix}/{prefix}.haplotypecaller.vcf.gz  vcf  records; header without ##GATKCommandLine
refuse          dict_mismatch   reference .dict differs from pinned GRCh38
refuse          multi_lane      more than one FASTQ pair (v1)
```

The comparator names map to a registry in code (`fastq_chunk`, `bam`, `md_metrics`, `raw`,
`vcf`, later `cram`, `gvcf`); the reference runner is a trait with one implementation
(Nextflow). A second descriptor (say UK Biobank BWA 0.7.17 + GATK 4.0.12) would add a file,
any new comparators, and a runner if it isn't Nextflow, and nothing in verify's core. Exact
paths are confirmed against a real sarek output tree during the build.

## 5. Output directory

```text
verify/
  verify.json        machine-readable report (schema below)
  verify.md          the same, readable
  slice/             R1/R2 slice (removed after PASS unless --keep-work)
  sarek/             results/, trace, logs (work/ removed after PASS unless --keep-work)
  nomics/            run directory, stages.tsv
  digests/           reference.tsv, nomics.tsv (the existing manifest format)
  logs/              verify.log, sarek.log, nomics.log
```

## 6. Report

### `verify.json` (schema `nomics-verify/1`)

```json
{
  "schema": "nomics-verify/1",
  "verdict": "PASS",
  "verdict_reason": "41/41 artifacts identical; 1,047 VCF records",
  "created_at": "2026-11-12T14:03:22Z",
  "mode": "reference-run",
  "descriptor": { "id": "sarek-3.10.0-germline", "sha256": "…" },
  "tools": {
    "nomics": { "version": "0.1.0", "commit": "c7203d3…", "features": ["cuda"], "device": "cuda" },
    "reference": {
      "pipeline": "nf-core/sarek", "version": "3.10.0",
      "commit": "8ccac7ad37b05dd792447763bf9671b719824587",
      "nextflow": "26.04.6", "config": "x86 pinned (Java PairHMM, cloud-hcjava.config)",
      "containers": [ { "process": "GATK4_HAPLOTYPECALLER", "image": "…", "digest": "sha256:…" } ]
    }
  },
  "host": { "arch": "x86_64", "cpu": "AMD EPYC 9R14", "cores": 20, "memory_gib": 160,
            "gpu": null, "os": "Ubuntu 24.04", "kernel": "6.8.0-1015-aws" },
  "inputs": {
    "fastq1": { "basename": "R1.fastq.gz", "bytes": 5123456789, "md5": "…" },
    "fastq2": { "basename": "R2.fastq.gz", "bytes": 5234567890, "md5": "…" },
    "targets": { "basename": "targets.bed", "md5": "…", "intervals": 212610 },
    "reference_dict_sha256": "…"
  },
  "slice": { "method": "head", "pairs_requested": 500000, "pairs_taken": 500000,
             "sha256": "…", "region": null },
  "artifacts": [
    { "name": "vcf", "comparator": "vcf",
      "normalization": "records exact; header without ##GATKCommandLine",
      "reference": { "sha256": "…", "count": 1047 },
      "nomics":    { "sha256": "…", "count": 1047 },
      "status": "PASS", "first_difference": null }
  ],
  "summary": { "artifacts": 41, "pass": 41, "fail": 0, "expected": 0 },
  "timing_s": { "preflight": 4, "slice": 21, "reference": 612, "nomics": 138, "gate": 9,
                "first_use_setup": 0 },
  "report_sha256": "…"
}
```

- `status` per artifact: `PASS`, `FAIL`, `MISSING` (on one side), `EMPTY` (count 0, an
  ERROR), or `EXPECTED` (only under `--also-sarek-default`, the GKL VCF).
- `first_difference` on a FAIL: `{ "index": 48211, "reference": "<record>", "nomics":
  "<record>", "fields": ["QUAL", "PL"] }`. For VCFs this is a variant record, which is not
  read data; for BAMs it gives the record index, name and differing fields only, without
  sequence or qualities.
- `report_sha256`: SHA-256 of the canonical JSON (keys in schema order, the field itself
  empty), so an edit is evident. `verify.md` repeats it. This is tamper-evidence, not
  signing.
- In `--against` mode, `mode` is `manifest`, `tools.reference` is replaced by the manifest's
  name and SHA-256, and timing has no `reference` entry.

### `verify.md` (sample, a FAIL; counts borrowed from tier0 for illustration)

```markdown
# nomics verify: FAIL

1 of 41 artifacts differs from sarek 3.10.0 on this machine.

| | |
|---|---|
| nomics | 0.1.0 (c7203d3), CUDA on NVIDIA GB10 |
| Reference | nf-core/sarek 3.10.0 (8ccac7ad), Nextflow 26.04.6, arm64 pinned config |
| Host | aarch64, 20 cores, 128 GiB, Ubuntu 24.04, kernel 6.11 |
| Input | R1.fastq.gz (md5 1a2b…), R2.fastq.gz (md5 3c4d…), targets.bed (212,610 intervals) |
| Slice | first 500,000 pairs, sha256 9f8e… |
| Time | sarek 10 min 12 s, nomics 2 min 18 s |

| Artifact | Normalization | Records (sarek / nomics) | Result |
|---|---|---|---|
| chunk/0001…R1 (24 files) | decompressed content | 40,797 / 40,797 … | PASS |
| mapped/….0001.bam (12 files) | records, header without @PG | 81,777 / 81,777 … | PASS |
| md_bam | records, header without @PG | 977,738 / 977,738 | PASS |
| md_metrics | bytes except lines 2 and 4 | 110 / 110 | PASS |
| recal_table | bytes | 1,090 / 1,090 | PASS |
| recal_bam | records, header without @PG | 873,049 / 873,049 | PASS |
| vcf | records; header without ##GATKCommandLine | 1,047 / 1,047 | **FAIL** |

## vcf: first difference

Record 512 (chr20:31,022,441): QUAL 812.64 (sarek) vs 812.65 (nomics); PL differs.

Report sha256: 7c1d…
```

(Chunk and mapped rows collapse to one line each in the readable report when all pass; every
file is listed individually in `verify.json`.)

## 7. How the brief's gates will be met

| Gate | Plan |
|---|---|
| Positive | HG002 tier0-scale slice PASS on the Spark (CUDA), Graviton (CPU) and x86 (CPU). |
| Negative | Fault-injection tests, each FAIL or ERROR and never PASS: one-record change in an output; a test-only nomics build with a changed tie-break (cargo feature `verify-fault`, never in release builds); a sarek task killed midway; empty outputs on both sides; a truncated BGZF file. |
| Fresh machine | Clean x86 and arm64 instances with Docker and the bundle only: release binary to verdict in 30 min or less with containers pre-pulled; first run with pulls timed separately. |
| Another sample | PASS on one non-HG002 sample from Brief 1 (or any GIAB sample). |
| Repeatable | Two runs, same verdict, slice digest and artifact digests; reports differ only in timestamps and timings. |
| CI | Unit tests (descriptor parsing, slice, report canonicalization, verdict logic) and a synthetic smoke test with no sarek run in the existing matrix. |

## 8. Files this touches

New: `crates/nomics/src/verify/` (mod, slice, descriptor, runner, gate, report),
`crates/nomics/descriptors/`, `docs/VERIFY.md`, tests under `crates/nomics/tests/`. Changed:
`crates/nomics/src/main.rs` (one `Verify` subcommand), `README.md` ("Check it yourself").
Small refactor: move the manifest read/compare helpers out of `main.rs` into
`nomics-gate` so verify and `dev digest-compare` share them. Not touched: `reference/`,
`docs/GATES.md`, `docs/REFERENCE-VARIABILITY.md`, `scripts/` (the sample-diversity work).

## 9. Decisions for review

1. **Which sarek configuration is the default?** Proposed: `pinned`, the configuration the
   goldens use (BAM outputs, `save_mapped`, Java PairHMM on x86), because it gates every stage
   boundary (41 artifacts at tier0) and is platform-independent. sarek's stock defaults
   publish CRAMs instead and, on x86, differ by design in about 2% of genome records (GKL);
   that stays available as `--also-sarek-default`, reported as EXPECTED. CRAM gating in
   verify would need samtools and htsjdk decoders in containers; proposed as a later
   addition.
2. **How is Nextflow provided?** Proposed: inside its pinned container, so Docker is the
   only prerequisite. The alternative is requiring Java 17 and Nextflow 26.04.6 on the host,
   which is simpler to build but one more thing for a user to install, and Nextflow's
   self-download is a network call outside the pinned containers.
3. **What is the default slice?** Proposed: `head:N`, the first N pairs: version- and
   machine-independent, no alignment pass, but shallow across the whole exome.
   `contig:chr20` gives tier0-like depth but needs a nomics alignment pass over the whole
   input first (slower, and the slice then depends on the nomics version). Both ship; the
   question is only the default, which the build then sizes by measurement.
   **Decided 2026-11-03: `contig:chr20` for exomes.** On a 30x exome on the Spark,
   `head:500000` gave 0.8x mean target depth; `contig:chr20` gave 58x and 1,060 calls and
   passed in 6 min 44 s. Genomes keep `head:4000000`, because a chr20 genome slice took 114 min on the Spark (46 of them slicing).
