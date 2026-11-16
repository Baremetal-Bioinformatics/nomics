# Structural costs of matching sarek's outputs (2026-10-20)

nomics: 5,240 s FASTQ → VCF at 30×, fused GPU extension, 41/41 identical to sarek.
Part of that time is work a pipeline that only needs the VCF would skip:
- trimming;
- writing the chunk BAMs, md.bam and recal.bam.

This note measures what that work costs us, and what sarek actually requires.

## What each pipeline writes

| Artifact | sarek default | Our golden config (`save_mapped`, `save_output_as_bam`) | nomics today | Parabricks |
|---|---|---|---|---|
| fastp chunk FASTQs | work dir only (`save_split_fastqs = false`) | work dir only | written only with `--write-intermediates` (gate) | no trimming |
| 12 per-chunk sorted BAMs | work dir only (`save_mapped = false`) | **published** | written; in whole-genome mode they are also markdup's spill | none (sorts on the GPU) |
| duplicate-marked alignments | **published as CRAM** (MarkDuplicates writes BAM, samtools converts) | md.bam | md.bam (htsjdk level-2 deflate, byte-exact) | md.bam (its own compression) |
| recalibrated alignments | **published as CRAM** | recal.bam | recal.bam (htsjdk level 2, byte-exact) | **never written**: HaplotypeCaller applies the table on the fly (`--in-recal-file`) |
| BQSR table, VCF | published | published | written | written |

Two consequences:
1. A lab on default sarek gets CRAMs, which nomics cannot produce yet. CRAM output is an
   adoption gap, separate from speed.
2. Without `save_mapped`, the chunk BAMs are never published, so in that configuration we
   write them only as our own spill.

## What it costs us (Spark, 30×)

**fastp chunk FASTQs (gate-only).** A full fused run without `--write-intermediates`
(`NOMICS_NO_INTERMEDIATES=1`) took 5,179 s against 5,240 s: **−61 s (1.2%)**, almost all
inside alignment (3,518 vs 3,578 s). The 17 artifacts that run wrote are gated separately
(see below).

**Trimming (fastp semantics).** Prep runs inside alignment: 330–460 s, overlapped. Chunk 1
aligned from finished prep (`--resume-from align`) took 292–306 s across this week's fused
runs. Assuming the other 11 chunks are similar, × 12 gives ~3,500–3,670 s, against
3,518–3,578 s in full runs. So trimming
**costs ~0–60 s** of wall time; it hides behind alignment.

**Writing the chunk BAMs during alignment.** Deflate is ~1.4% of CPU samples on the fused
chunk (`perf`, `ZlibFast`): **tens of seconds at most**, overlapped. In whole-genome mode,
markdup re-reads these files twice; that inflate is counted below.

**md.bam and recal.bam.** Markdup + BQSR at 30× (`--resume-from markdup`, stop after BQSR):
markdup 308 s + BQSR 395 s = **703 s**, writing md.bam (44.8 GB) and recal.bam (69.7 GB).
CPU samples:

| Work | Share of markdup + BQSR CPU |
|---|---:|
| Deflate (htsjdk-exact level 2: `deflate_finish`, `compress_block`, `encode_member`, Huffman) | **~37%** |
| Inflate (re-reading chunk BAMs twice and md.bam twice) | **~21%** |
| BQSR covariates and apply | ~20% |
| markdup sort and classification | ~5% |

The stage is CPU-bound, so those shares translate roughly to wall time: **~250 s of deflate
and ~150 s of inflate** out of 703 s.

A pipeline that applies the recalibration table on the fly inside HaplotypeCaller never
writes a recal.bam.

## What is structural and what is ours to fix

| Item | Seconds at 30× | Required to match sarek? | Lever (all exact) |
|---|---:|---|---|
| fastp chunk FASTQs | ~60 | No: gate only | Skip in production (`--write-intermediates` off) |
| Trimming | ~0–60 | Yes | none needed; it overlaps |
| Chunk BAM writes in align | < 50 | Only with `save_mapped` | lighter spill format when not published |
| Chunk BAM re-reads in markdup (inflate) | part of ~150 | No (our spill design) | lighter spill format; read once |
| md.bam write | part of ~250 | Yes, as BAM or CRAM | faster deflate: the gate is on records and the uncompressed stream, not bytes |
| md.bam re-read by BQSR (inflate) | part of ~150 | No | build BQSR tables during markdup pass 2 (board I-06) |
| recal.bam write | part of ~250 | Yes (sarek publishes it) | faster deflate, as for md.bam |

Neither BAM gate requires htsjdk's exact compressed bytes. GATES.md gates md.bam and
recal.bam on records plus the uncompressed stream. The `.bai` gate is two-sided, our index
of our own BAM, so it holds at any compression level. We write htsjdk's bytes exactly by
choice. That choice makes files byte-comparable to the goldens, and it costs about a third of
markdup + BQSR.

Estimate for the exact, config-preserving levers:
- a faster deflate for md.bam, recal.bam and the chunk BAMs: −100 to −180 s;
- I-06 (one fewer md.bam read): −50 to −100 s;
- a light spill when chunk BAMs aren't published: −30 to −80 s, sarek-default configuration
  only.

Together, about **−150 to −300 s of duplicate marking + BQSR**. What remains
is recal.bam itself: writing 70 GB that sarek publishes.

## Alignment: the larger part of the gap

The clean nsys trace of the fused chunk (2026-10-20, index fully on 2 MiB pages) gives:
- GPU busy 206 s of a 271 s span (76%): extension 101 s, rescue 84 s, seeding 80 s;
- CPU about 16 GB/s of memory traffic (calibrated counters), with no saturation in sight.

Neither side is full, so alignment waits on its stage barriers. The fused-path CPU profile
shows what the CPU spends its time on:

| CPU work (fused path) | Share of samples |
|---|---:|
| Region merging: `sort_dedup_patch`, `patch_reg` → `gen_cigar2` → `global2`, and most reference fetches (`get_seq`) | ~30–35% |
| Chaining: `mem_chain_from`, `mem_chain_flt` | ~22% |
| CPU re-seeding of reads that overflow the GPU's 64-interval cap (~4% of reads) | ~18% |
| BAM encode + deflate | ~2% |

These add two cheap, exact levers to the list in the alignment-gap scope:
- **two-pass GPU seeding** in the pipeline: 11× fewer CPU fallbacks, for ~160 ms of GPU per
  batch;
- **a faster region merge.** It is a dump-replay candidate, like S-08 was for extension.

## Gates

- The production run (no intermediates) wrote 17 of the 41 artifacts: 12 mapped BAMs,
  md.bam, metrics, recal table, recal.bam and the VCF. **All 17 PASS** against the golden
  digest manifest. The other 24 lines are the fastp chunks, digested from the reference run's
  prep (symlinked), so they don't count for this run.
- Every profile chunk passed its byte gate.

## Update 2026-10-23: fast BAM deflate (G-07, 7a4afab)

A from-scratch single-probe DEFLATE encoder now writes md.bam, recal.bam and the chunk BAMs.
It keeps the same BGZF block cuts as htsjdk and adds an ARMv8 CRC32. The byte-exact htsjdk
mode is kept behind `NOMICS_BAM_COMPRESSION=htsjdk`.

| 30× | htsjdk level 2 (exact) | fast |
|---|---:|---:|
| markdup | 283–284 s | 215–216 s |
| BQSR | 392 s | 311–313 s |
| **markdup + BQSR** | **675 s** | **528 s (−22%)** |
| md.bam | 44.83 GB | 41.26 GB |
| recal.bam | 69.73 GB | 66.67 GB |
| chunk BAMs | 59.59 GB | 56.41 GB |

Markdup + BQSR is now ~165 s shorter. Most of what remains is writing recal.bam,
which sarek publishes.

## Update 2026-10-24: I-06, two-pass GPU seeding, light spill

Three exact levers, measured at 30× on the Spark (fused path, not freshly rebooted; run dirs
archived under `runs-archive/session-20261024`).

| 30× fused | Total | align | markdup | BQSR | HC | Gates |
|---|---:|---:|---:|---:|---:|---|
| Clean run, Oct 23 (after reboot) | 4,003 s | 2,495 | 222 | 319 | 944 | 41/41 |
| A: all three levers, chunk BAMs published | **3,974 s** | 2,535 | 285 | 199 | 937 | 41/41 |
| B: resume from A's chunk BAMs, I-06 off | — | — | 233 | 314 | 936 | 41/41 |
| C: `--save-mapped false`, no gate FASTQs | **3,961 s** | 2,514 | 290 | 198 | 940 | 5/5 published |

- **I-06 (BQSR tables built during markdup's last pass): −63 s.** A vs B on the same chunk
  BAMs: markdup + BQSR 484 s vs 547 s. The recal table is identical, since every accumulator
  is an integer count. Markdup itself grows by ~52 s: the table builder (busy ~185 s, 20
  threads) shares cores with markdup's rewrite and deflate. Standalone, the table pass took
  ~115 s. Peak RSS grows ~5 GB (77 vs 72.5 GB) from the reference and known sites the builder
  holds during markdup.
- **Two-pass GPU seeding: exact, ~0 s end to end.** On the 30× chunk, CPU re-seeds drop 11.3×
  (2.33M → 206k reads) for 66 GPU-s of cap-128 retries over 89 batches. The CPU-extension
  path gains ~10% (align 265 vs ~295 s), but the fused path only ~1.5% (197–199 vs ~202 s),
  because its GPU is the busier side and the retries add to it. Full-run alignment was 2,514–
  2,535 s vs 2,495 s on the clean run; these runs were not after a reboot, so that ~1–2% is
  within the fragmentation spread we have seen. Metal keeps cap 64 by default.
- **Light spill (`--save-mapped false`): ~0 s, more disk.** The `.nspill` chunk files decode
  ~4× faster per core than BGZF, but markdup with them (290 s) was no faster than with fast
  BAM chunks (285 s): markdup is bound by its rewrite, deflate and the table builder, not by
  inflate. They take ~79 GB against ~56 GB of BAM. Skipping the gate-only fastp FASTQs saved
  little here (prep 343 vs 329 s, align −21 s).

The estimate was −150 to −300 s; the measured gain is ~30–40 s, nearly all of it I-06. The
remaining markdup + BQSR time is mostly writing recal.bam (or
recal.cram), which sarek publishes.

## Update 2026-10-25: CRAM mode from the stages' own streams

At 30x, `--save-output-as-bam false` cost ~450 s over BAM output (4,408 vs 3,961 s). md.cram
and recal.cram were re-encoded from md.bam / recal.bam on all cores alongside BQSR (+135 s)
and HaplotypeCaller (+260 s). Two commits address this; each can be measured or reverted on
its own.

1. **Pipeline** (`run: CRAM mode encodes md.cram from markdup's last pass ...`).
   - markdup's pass-2 stream is tee'd to the I-06 table builder and an md.cram encoder.
   - ApplyBQSR's output stream is tee'd to a recal.cram encoder.
   - Without `--keep-bams`, md and recal travel between the stages as internal `.nspill`
     files.
   - md.nspill is removed after BQSR, and recal.nspill after HaplotypeCaller, each once
     its CRAM is complete.
   - HaplotypeCaller then runs with no encoder beside it.
2. **Encoder** (`cram: codec trials ...`). Each block tried rANS o0 and o1 plus fast gzip,
   then zlib-5 when gzip was close.
   - Now a trial container tries every codec, once per 32 containers: one with ≥1,000
     records, so a small contig's few reads never set the codecs.
   - The containers after it reuse that trial's codec for each block id.
   - Roles depend on container indices and data only, so 1-thread and N-thread files are
     byte-identical.

Per record, one Mac P-core, CPU time. Header M5 and slice MD5 are subtracted: they track
reference span (about the genome once per file), not record count. Columns:
- old enc: encoder before the codec commit.
- trial: codec commit at tier0, where few containers follow a trial.
- follower: all containers follow (as at 30x).
- BAM w/r, spill w/r: writing and reading the stream as BAM or as a spill file.

| µs / record | old enc | trial | follower | BAM w / r | spill w / r |
|---|---:|---:|---:|---:|---:|
| md (977,738 records) | 2.70 | 2.45 | 1.61 | 0.64 / 0.46 | 0.38 / 0.24 |
| recal (873,049 records) | 6.14 | 4.45 | 3.41 | 1.03 / 0.69 | 0.73 / 0.21 |

Sizes are unchanged within 0.1%: md.cram 13.19 → 13.21 MB, recal.cram 44.02 → 44.05 MB.
- Dropping zlib-5 entirely would cut recal to 1.44 µs, but recal.cram grows 16% (51.0 MB).
- zlib level 4 gives 2.66 µs, but recal.cram grows 8% (47.6 MB).
- Both rejected. zlib-5 on quality blocks is most of the remaining recal cost.

**30x estimate** (703.5M md and 676.6M recal records), extra CPU over BAM mode:

| | Before | After |
|---|---:|---:|
| md.cram + recal.cram encoding | ~6,050 CPU-s | ~3,440 CPU-s (followers) |
| re-reading md.bam / recal.bam | ~790 CPU-s | none |
| writing and reading md and recal as spill instead of BAM | — | about −860 CPU-s |
| **Total** | **~6,850 CPU-s** | **~2,580 CPU-s** |

The encoders now run paced by the stage producing the records, through bounded channels.
That should make the cost close to zero-sum on 20 cores, ~130 s, against the +450 s seen,
which was super-linear contention. Expected CRAM-mode FASTQ→VCF is ~4,100–4,150 s.

The assumptions are unmeasured: Grace cores were taken as equal to the Mac's, and
contention as zero-sum.

**Disk**: on tier0, spill is 1.43× BAM for md and 1.36× for recal.
- Peak internal disk, during BQSR: md.nspill + recal.nspill, ~1.4× (md.bam + recal.bam).
- During HaplotypeCaller: recal.nspill alone, ~1.36× recal.bam, below the old md.bam + recal.bam.

tier0 (Mac, 14 threads, `--spill-prep`, streaming), user+sys CPU:
- Base: 116.5 and 115.7 s.
- Pipeline commit: 114.4, 114.4 and 113.7 s.
- Pipeline + codec commits: 112.5 and 111.8 s.

Stage walls are 1–5 s and do not discriminate. All CRAM gates pass, and the recal table and
VCF digests pass. BAM mode is record- and header-identical to before.

**Stream integrity (review follow-up).** Streams reach their final names only once complete.
- md and recal, as BAM or spill file, are written under `.partial` names, synced, renamed and
  directory-synced.
- Stage readers refuse a BAM without its BGZF EOF marker.
- Each stream a run writes gets a generation in `<stream>.gen`. For recal it also names the
  md generation it was made from.
- CRAM `.done` markers (v2) name their source generation. Resume takes the latest complete
  stream and regenerates any CRAM whose marker names another stream.
- The other format's copy of a rewritten stream is removed only once its replacement is
  recorded.
- The CRAM writer's reorder buffer is bounded at 2 × threads + 2 containers.
- A stream keeps its `.gen` until a new stream takes its name; the writer drops it just
  before then. A stream with no matching `.gen` (an older binary, or a run stopped before
  recording) is never a CRAM's recorded source.
- On resume at HC it is refused, as is a recal whose md parent is not the md beside it:
  resume from bqsr instead.
- On resume at BQSR an unrecorded md is recorded under a new generation.
- `BamStream::open` accepts any BAM. The EOF/footer check applies only to the pipeline's
  own streams.

tier0 failure arms:
- md.bam cut before its EOF block: resume refused.
- CRAM run A, then BAM-mode run B rewriting the streams, then resume of B at HC in CRAM mode:
  both CRAMs regenerated from B's streams, gates pass.

### Measured at 30× (Spark, 2026-10-25, 712fb4c, fused, sarek-default output)

`--save-output-as-bam false --save-mapped false`, no gate FASTQs. The run is flagged
timing-suspect: after 7 h of uptime only 64% of the index got 2 MiB pages, which mostly inflates
align. Compare the downstream stages:

| CRAM mode | markdup | BQSR | HC | markdup→VCF | over BAM mode |
|---|---:|---:|---:|---:|---:|
| 2026-10-24, re-encoded from md.bam/recal.bam | 300 | 333 | 1,200 | 1,833 s | +405 s |
| **2026-10-25, encoded in-stream** | **376** | **317** | **939** | **1,632 s** | **+204 s** |
| BAM mode (run C, 2026-10-24) | 290 | 198 | 940 | 1,428 s | — |

- **Where the CRAM time now goes:** md.cram is encoded during markdup's last pass (371 s
  inside markdup), and recal.cram during BQSR's apply (317 s). HaplotypeCaller runs with no
  encoder beside it and is back to its BAM-mode time.
- **The outcome is about half the estimate:** CRAM's overhead halved (about −200 s). The
  estimate was a further 270–320 s saved, landing near +130–180 s over BAM mode instead of +204 s.
- **Disk and outputs:** peak disk 323 GB above the start. md.cram is 12.65 GB and recal.cram
  36.55 GB; no `.bam` or `.nspill` is left behind.
- **Gates:**
  - The recal table and the VCF are identical; recal.cram has 676,643,121 records.
  - md.cram passed under htsjdk 4 and on every htslib region query and count.
  - The htslib full decode of the reference failed on a disk filled by Docker's json-file log
    of the decoder's stdout. That is now fixed: `--log-driver none` in `cram_gate.sh` and in
    the Spark's samtools wrapper.
  - The lane's run of the same check (2a03e00) passed htslib in full.
