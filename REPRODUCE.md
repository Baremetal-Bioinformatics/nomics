# Reproducing the results

There are two levels of checking. The first needs only nomics, the public 30× FASTQs and
the digest manifest shipped in `reference/`. The second regenerates the sarek goldens
themselves, which takes about 12 hours of compute for the 30× genome.

The scripts assume this fixed layout under `$HOME/omics`:

```
~/omics/
  data/     FASTQs and the exome BED
  ref/      GATK.GRCh38 (iGenomes layout) and sarek-3.10.0
  cache/    nomics' index caches (built on first use)
  runs/     sarek goldens and nomics runs
  new_omics/  this repository
```

## 1. Inputs

| Input | Source | Check |
|---|---|---|
| HG002 exome, Agilent v7 (tier1) | ENA [SRR14724482](https://www.ebi.ac.uk/ena/browser/view/SRR14724482), project PRJNA734598 | md5 `191fc6bfd2492b4ccc436b34f0fdbbba` (R1), `c104c7cc0411652c1a5c7d775554f80d` (R2) |
| HG002 30× whole genome | ENA [SRR14724532](https://www.ebi.ac.uk/ena/browser/view/SRR14724532) | |
| Agilent v7 targets | `s3://ngi-igenomes/test-data/sarek/Agilent_v7.bed` | 212,610 targets |
| GRCh38 reference and GATK bundle | `s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38` (FASTA, `.fai`, `.dict`, BWA index, dbSNP 146, Mills, known indels, WGS calling regions) | 12 GB |

Put the FASTQs at the paths in `reference/samplesheet.*.csv` (`data/fastq/` for the exome,
`data/sra_wgs/fastq/` for the genome). Download with several connections (for example
`aria2c -x 16`); a single HTTPS stream from ENA can be very slow. Verify the md5s: one
early download had the right size and the wrong content.

The smaller tiers are derived from these. tier0 needs the tier1 golden first: `nomics dev
tier0` selects its pairs from sarek's tier1 md.bam (see `crates/nomics/src/tier0.rs`).
tier0nv is made from tier0 with `reference/make_tier0nv.py`, and the 1/12 genome slice is the
first 1/12 of the 30× read pairs.

## 2. Build

```bash
cargo build --release -p nomics --features cuda
```

```bash
cargo test --workspace
```

Without `--features cuda` you get the CPU-only build (any Linux or macOS machine); the tests
run in either. Some tests need large fixtures that are not in the repository; they skip
unless `NOMICS_FIXTURES` points at a directory holding them.

## 3. Check a run against the shipped digests (no sarek needed)

`reference/wgs30-x86-r1.digest.tsv` (and, for tier0, `reference/tier0-wes-20261008-r1.digest.tsv`)
holds a normalized SHA-256 digest and record count for every artifact sarek published in the
golden run, computed with the normalizations in [docs/GATES.md](docs/GATES.md). A nomics run whose
digests match has the same content as sarek's, artifact by artifact.

30× whole genome (about an hour on a DGX Spark; peak memory about 77 GB):

```bash
scripts/full_exome.sh wgs30
```

The script evicts page caches (and compacts memory if a sudoers rule allows `tee
/proc/sys/vm/compact_memory`; otherwise that step is skipped), runs `nomics run` with GPU
seeding and the golden run's scatter order and read-group `DS`, then runs `nomics dev digest`
and `nomics dev digest-compare` against the manifest. It prints the stage timings, one PASS
line per artifact (41) and a verdict. If GPU
buffers or the index fall back to small pages because memory is fragmented, it writes
`timing_suspect.txt`: the outputs are still exact, but the time is not comparable. A fresh
boot gives the best times.

tier0, on any machine (a few minutes; `cpu` or `metal`), once you have the tier0 FASTQs
in `~/omics/data/tier0/`:

```bash
scripts/mac_tier0.sh cpu
```

## 4. Regenerate the goldens

Each golden is a sarek 3.10.0 run. On arm64, `reference/run_sarek.sh <run-dir> <samplesheet>
<params.yaml>` runs it with the pinned Nextflow and configuration, recording start and finish
times (sarek's `samtools merge` seeds its random ID suffixes from the clock); on x86,
`reference/run_sarek_cloud.sh` does the same with sarek's default containers:

| Tier | Samplesheet | Params | Run directory the gates expect |
|---|---|---|---|
| tier1 | `samplesheet.tier1.csv` | `params.tier1.yaml` | `runs/tier1-wes-20261008-r1` |
| tier0 | `samplesheet.tier0.csv` | `params.tier1.yaml` | `runs/tier0-wes-20261008-r1` |
| tier0nv | `samplesheet.tier0nv.csv` | `params.tier1.yaml` | `runs/tier0nv-wes-20261012-r1` |
| wgs20 | `samplesheet.wgs20.csv` | `params.wgs20.yaml` | `runs/wgs20-20261010-r1` |
| wgs30 (x86) | `samplesheet.cloud-wgs30.csv` | `params.cloud-wgs.yaml` | `runs/wgs30-x86-r1` |

Notes:

- **Run each golden twice** and compare with `reference/compare_runs.sh` before trusting a
  gate; [docs/REFERENCE-VARIABILITY.md](docs/REFERENCE-VARIABILITY.md) lists what differs.
- **Platform matters.** The arm64 goldens use sarek's arm64 containers plus a locally built
  GATK image (`reference/containers/gatk4/`), since sarek's GATK image has no arm64 build.
  The 30× golden ran on x86. Its VCF was then recomputed with `cloud-hcjava.config` added
  (`-c reference/cloud-hcjava.config`), which switches HaplotypeCaller to the Java PairHMM
  and Smith-Waterman and turns off HotSpot's libm intrinsics, and the 21 group VCFs were
  merged with MergeVcfs. That makes x86 GATK agree with arm64 GATK record for record; sarek's
  default x86 run (Intel GKL) calls about 2% of 30× records differently. Everything upstream
  of HaplotypeCaller is the same on both platforms.
- **Record the scatter order.** Whole-genome outputs depend on the order in which sarek's
  parallel tasks finish. `reference/*.scatter-order.txt` holds it for each golden, and
  `nomics run --scatter-order` reproduces it.
- `scripts/full_exome.sh <tier>` gates every artifact of a nomics run against a golden
  directory (`tier1`, `tier1joint`, `tier0`, `tier0nv`, `wgs20`) or, for `wgs30`, against the
  digest manifest.

## 5. Benchmarks

`scripts/bakeoff/spark_campaign.sh` runs the sarek and nomics timings in [RESULTS.md](RESULTS.md):
one run at a time, page caches evicted between runs, a 5-second CPU/GPU/memory sampler
alongside. It can also time third-party pipelines on the same inputs; those need their own
containers and licenses, and their results are not published here.
