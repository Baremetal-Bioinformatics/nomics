# `nomics verify`: check nomics against sarek yourself

`nomics verify` answers one question: **is nomics identical to sarek on my reads, on my
machine?** It takes a deterministic slice of your FASTQs, runs pinned nf-core/sarek 3.10.0
and nomics on that slice, compares every intermediate and final artifact under the
normalizations in [GATES.md](GATES.md), and writes a sealed report. You do not have to trust
our numbers: sarek runs on your machine, from pinned containers, and the report says exactly
what ran.

Status: the Spark (GB10, arm64, CUDA) gate has passed (2026-11-02/03): tier0 twice with
identical digests, `--against` the tier0 golden, and `--offline` with no network (about
2 min 15 s per run once nomics' index caches exist, 3 min 40 s the first time); the full
exome on `head:500000` and on `contig:chr20`; WGS chr20 (`--wgs`, 53 min, almost all of it
sarek's one HaplotypeCaller task); a NovaSeq lane; and `--against` the full 30x WGS x86
golden (41/41 identical, 4,773,797 VCF records). Every input on the Spark is HG002.

Off the Spark (2026-11-11): stock Ubuntu 24.04 EC2 instances with nothing of ours installed,
one x86_64 (m7i.8xlarge) and one Graviton3 (m7g.8xlarge), each building nomics from source.
Both pass tier0 against sarek run on that machine (41/41), tier0 `--against` the golden,
tier0 `--offline`, and the full exome on the default `contig:chr20` slice (41/41, 1,060
calls; 16 min 26 s on x86, 12 min 41 s on Graviton, CPU only). The first attempts failed, and
both failures are fixed in this version: on x86 the pinned Nextflow image's Docker client was
too old for Docker Engine 29, and on Graviton3 rustc 1.98.1 miscompiled one sum for 256-bit
SVE ([LEARNINGS.md](LEARNINGS.md)). The reference, the exome and the tier0 fixture came from
our own bucket, checked against the public checksums, not straight from iGenomes and ENA. The tier0 runs used `--slice head:500000`, which takes all of
tier0. On the full 30x exome (SRR14724482) the default `contig:chr20` slice passed with
1,060 calls at 58x mean target depth in 6 min 44 s with GPU slicing (11 min 19 s on the CPU).

## What you need

- Linux on x86_64 or aarch64, with Docker (the daemon reachable through
  `/var/run/docker.sock`). Docker is the only prerequisite: Nextflow itself runs inside its
  pinned container (`nextflow/nextflow:26.04.6`).
- The iGenomes `GATK.GRCh38` reference directory. verify refuses any other: its `.dict`
  must match the public GRCh38 dictionary contig for contig.
- One lane of paired FASTQs (gzip or BGZF). Multi-lane input is refused in this version.
- An exome target BED, or `--wgs` for whole-genome reads.
- About 20 GB of free disk for the cache and run directory, and the memory nomics needs for
  its index (see [REPRODUCE.md](../REPRODUCE.md)).

## Running it

```bash
# Once, with network: the Nextflow image (a local build of the pinned one: on
# arm64 because Nextflow publishes none, on x86 to give it a Docker client that
# current Docker engines accept), sarek source, Nextflow plugins, and on arm64
# the pinned task containers.
nomics verify fetch
```

```bash
# Exome.
nomics verify --fastq1 R1.fastq.gz --fastq2 R2.fastq.gz \
  --reference /path/to/GATK.GRCh38 --targets targets.bed --outdir verify-out/
```

```bash
# Whole genome: no targets; calling is limited to one contig (chr20 by default).
nomics verify --fastq1 R1.fastq.gz --fastq2 R2.fastq.gz \
  --reference /path/to/GATK.GRCh38 --wgs --outdir verify-out/
```

```bash
# Check the report has not been edited since verify wrote it.
nomics verify check-report verify-out/verify.json
```

On x86 the first run pulls sarek's task containers; on arm64 `fetch` already did. Later runs
can add `--offline`: Nextflow and every task container then run with `--network none`, and
the run fails if something pinned is missing.

`nomics verify --help` lists every flag. The ones you are most likely to change:

| Flag | Default | Meaning |
|---|---|---|
| `--slice contig:NAME` | `contig:chr20` exome | Pairs whose nomics alignment touches NAME: a full-depth exome on one chromosome. |
| `--slice head:N` | `head:4000000` genome | The first N read pairs, copied byte for byte. No alignment pass, but shallow: on a 30x exome, `head:500000` gave 0.8x target depth. |
| `--device auto\|cpu\|cuda` | `auto` | nomics side only; sarek is CPU either way. |
| `--threads N` | all cores | nomics threads. sarek keeps the goldens' per-task caps. |
| `--keep-work` | off | Keep the slice, sarek's `work/` and nomics' intermediates after a PASS. |
| `--also-sarek-default` | off | x86: also run sarek's stock configuration and report its VCF difference as EXPECTED. |
| `--online` | off | Allow network beyond the pinned pulls (for example Seqera Wave builds). Not needed on any platform. |

## Verdicts

| Verdict | Exit | Meaning |
|---|---|---|
| **PASS** | 0 | Every artifact matches, none is empty, and the sarek VCF has at least `--min-records` (200) records. |
| **FAIL** | 1 | Both sides finished and at least one artifact differs, is missing on the nomics side, or cannot be read. The report names each one and its first differing record. |
| **ERROR** | 2 | No real comparison was possible: Docker missing, refused input, a failed sarek task, empty outputs on both sides, too few calls. Never reported as PASS. |

## What runs

1. **Preflight.** Docker answers; the reference matches GRCh38; the BED parses; the output
   directory is new or empty.
2. **Slice.** `contig:NAME` aligns every pair with nomics (on `--device`) and keeps the pairs
   that touch NAME; it reads the reference's `BWAIndex` and builds nomics' index cache on
   first use. Which pairs it keeps can change between nomics versions; both sides then run
   on the same slice, so the comparison holds, and the report names the slice. `head:N`
   copies the first N pairs instead, with no alignment. Either way R1 and R2 names must agree
   record for record and neither file may end early. The slice digest is taken over the
   uncompressed records kept, so it does not depend on the compressor, and two runs on the
   same input report the same digest.
3. **sarek 3.10.0** (commit `8ccac7ad`) with the goldens' parameters: `tools:
   haplotypecaller`, `skip_tools: haplotypecaller_filter`, `save_output_as_bam`,
   `save_mapped`, bwa-mem. On x86 it adds `reference/cloud-hcjava.config` (Java PairHMM, so
   results do not depend on the CPU's vector units). On arm64 it turns Seqera Wave off and runs
   each task in the image the arm64 goldens ran, pinned by digest in the descriptor; every
   task in the trace must have used its pin. Task CPUs are capped at 20 on every machine, because fastp's chunking depends on
   its task CPUs and the goldens ran with that cap. Every task in the trace must end
   COMPLETED with exit 0.
4. **Capture** from the sarek run: the read-group `DS` field and, for multi-interval runs,
   the scatter order and merge seed.
5. **nomics** on the same slice, writing intermediates.
6. **Gate.** Each artifact is digested on both sides with its normalization; on a mismatch,
   the artifact's comparator finds the first differing record.
7. **Report** in `verify.json` and `verify.md`. After a PASS, the slice and both work trees
   are removed unless `--keep-work`; after FAIL or ERROR they stay for inspection.

No read data leaves your machine or enters the report. The report records input sizes and
md5s, not paths.

## What is compared

| Artifact | Normalization |
|---|---|
| `chunk/*` (fastp chunks) | decompressed content of each chunk |
| `mapped/*` (per-chunk BAMs) | records, and the header without `@PG` |
| `md_bam` | records, and the header without `@PG` |
| `md_metrics` | bytes, except lines 2 and 4 (command line and date) |
| `recal_table` | raw bytes |
| `recal_bam` | records, and the header without `@PG` |
| `vcf` | records, and the header without `##GATKCommandLine` |

These come from the pipeline descriptor, `crates/nomics/descriptors/sarek-3.10.0-germline.tsv`,
which is data, not code: the pinned source, Nextflow image, per-architecture profiles,
parameters, and where each artifact lives on each side. A second pipeline is a second
descriptor plus, only if it produces a file type none of the comparators reads, a new
comparator.

## Output directory

```text
verify-out/
  verify.json     machine-readable report, schema nomics-verify/1, sealed with its SHA-256
  verify.md       the same, for people
  digests/        reference.tsv and nomics.tsv, one line per artifact
  inputs/         samplesheet, params and Nextflow overlay as run
  logs/           Nextflow, sarek trace, nomics
  slice/ sarek/ nomics/   removed after a PASS unless --keep-work
```

`verify.json` records the verdict and reason; the nomics version, commit and build features
(the commit ends in `-dirty` when tracked files differed from it at build time); the
descriptor's digest; the host's CPU, memory, GPU and OS; the slice method and digest;
the Nextflow image digest and every task container's image digest; and per artifact, both
digests, record counts and the first difference.

## Checking without sarek

`--against MANIFEST` skips sarek and compares nomics with a digest manifest shipped in
`reference/`. It needs the exact FASTQs the manifest was made from (md5-checked: a
re-download or re-compression of the same run has other bytes and is refused).

| Manifest | Input | Flags | Spark wall |
|---|---|---|---|
| `tier0-wes-20261008-r1.digest.tsv` | tier0 exome fixture | `--targets`, and `--read-group-ds` with the FASTA path its golden run's bwa saw | 10 s |
| `wgs30-x86-r1.digest.tsv` | HG002 30x WGS, SRR14724532 as SRA-toolkit FASTQ (not ENA's gz) | none: whole-genome calling regions, the golden's read group and scatter order are built in; `--targets` is refused | about 80 min on the Spark (78 min measured: 66 min nomics, 10 min gate) |

```bash
nomics verify --against reference/wgs30-x86-r1.digest.tsv \
  --fastq1 SRR14724532_1.fastq.gz --fastq2 SRR14724532_2.fastq.gz \
  --reference /path/to/GATK.GRCh38 --outdir verify-wgs30/
```

## Checking that verify can fail

The unit tests inject five faults into synthetic outputs (a changed record, a changed read,
empty outputs, a truncated BGZF file, a missing output) and check each gives FAIL or ERROR.
For a fault inside nomics itself, build the test-only `verify-fault-tiebreak` feature: among
duplicates with equal scores it keeps the last, where Picard keeps the first. The
report lists the feature, and the binary warns on every run. Run it on real reads; it must
FAIL, first at `md_bam` (which reads carry the duplicate flag). On the Spark, `--against` the
tier0 golden, it does: `md_bam`, `recal_table`, `recal_bam` and `vcf` differ, while the
chunks, mapped BAMs and duplicate metrics (the same counts either way) still match.

```bash
cargo build --release -p nomics --features cuda,verify-fault-tiebreak --target-dir target-fault
target-fault/release/nomics verify --fastq1 R1.fastq.gz --fastq2 R2.fastq.gz \
  --reference /path/to/GATK.GRCh38 --targets targets.bed --outdir verify-fault/
```

## On the Spark runner

The development repository runs these checks on a self-hosted DGX Spark runner on every push
to a `verify/**` branch: `tier0`, `tier0-against` by default, and `exome-default` or
`fault-tiebreak` when the head commit's message carries `[spark: ...]` with their names. Each check states the verdict it
expects; the job summary says whether every one ended that way. First runs, 2026-11-03:
tier0 PASS in 2 min 13 s, `--against` tier0 PASS in 7 s, the full exome on the default
chr20 slice PASS in 6 min 54 s, and the tie-break build FAIL as intended in 8 s.
`wgs-default` runs the genome default slice (`head:4000000`) and `wgs-chr20` the 30x
genome's `contig:chr20` slice. `wgs-default` passed on 2026-11-03 (41/41 identical, 4,462
VCF records, 9 min 20 s).
Each check ends by printing where its time went: verify's steps, then sarek's processes by
task time.

## Limitations

- **arm64 pins come from the goldens' images.** The 11 frozen Wave images the arm64 goldens
  ran (named in their `report.html`) are pinned by registry digest. A workflow change that
  adds a process needs a pin for it, or the run ends in ERROR. Calling regions that form more
  than one scatter group add two such processes (the GATK gather steps), whose stock image
  has no arm64 build.
- **x86 task containers are pinned by tag** in sarek; verify records their digests in the
  report, so a later change shows.
- **One lane, GRCh38 only.**
- **Whole genomes do not fit 30 minutes at depth.** On the Spark, the 30x genome's
  `contig:chr20` slice (8.3M pairs, 124,492 calls) took 114 min: 46 min to slice (the slicer
  reads all of the genome's pairs to pick out chr20's), 62 min in sarek, 3.5 min in nomics.
  That run cut the calling regions into 4 Mb pieces to call them in parallel; one piece ran
  50 min, so the calling time sits in one region, and the option was removed. Slicing alone is
  over the budget. The genome default stays `head:4000000`: 9 min 20 s, 4,462 calls, shallow.
  An earlier 53 min `--wgs` run on the Spark reported 118,062 calls from `head:4000000`, which
  does not match today's 4,462 from the same flag; its input was probably a different FASTQ.
- **Exome wall times come from three machines.** The chr20 default took 6 min 44 s on the
  Spark with the GPU, 12 min 41 s on a 32-core Graviton3 and 16 min 26 s on a 32-vCPU x86
  instance (CPU only).
