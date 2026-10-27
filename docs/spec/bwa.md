# BWA-MEM byte-exactness specification (bwa 0.7.19-r1273 + samtools sort 1.22.1 / htslib 1.22.1)

Target: reproduce, record-for-record and byte-for-byte (modulo BGZF framing / libdeflate, which is out of scope), the
per-chunk coordinate-sorted BAM produced by nf-core/sarek 3.10.0's `BWAMEM1_MEM` process:

```
bwa mem -K 100000000 -Y -R "@RG\tID:HG002_WES.L1\tPU:L1\tSM:HG002_HG002_WES\tLB:HG002_WES\tDS:/home/<user>/omics/ref/GATK.GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta\tPL:ILLUMINA" -t 20 \
  Homo_sapiens_assembly38.fasta.64 0001.HG002_WES-L1_R1.fastp.fastq.gz 0001.HG002_WES-L1_R2.fastp.fastq.gz \
  | samtools sort --threads 20 -o HG002_WES.0001.bam -
```

Every claim is tagged **[source]** (repo/file/function/lines at the pinned tag) or **[inferred]**.

Pinned sources:

* `bwa` = github.com/lh3/bwa tag `v0.7.19`. Version verification **[source]**: `main.c` L33 `#define PACKAGE_VERSION "0.7.19-r1273"`,
  so the tag is exactly the binary's version string (`[main] Version: 0.7.19-r1273`, `@PG ... VN:0.7.19-r1273`, `main.c` L94, L122).
* `samtools` = github.com/samtools/samtools tag `1.22.1` (`version.sh` L27 `VERSION=1.22.1`; `bamtk.c` L79-81 `samtools_version()` returns it). **[source]**
* `htslib` = github.com/samtools/htslib tag `1.22.1`.

Notation: `l_pac` = forward reference length (sum of contig lengths), the "concatenated" coordinate space is `[0, 2*l_pac)` with
`[l_pac, 2*l_pac)` being the reverse complement of the whole forward sequence (`bntseq.c` `bns_fasta2bntseq` L306-312). A position `p >= l_pac`
is on the reverse strand; `bns_depos(p) = 2*l_pac - 1 - p` is its forward coordinate (`bntseq.h` L87-90). **[source]**

---

## 0. Execution path, options, and global settings

### 0.1 Option parsing (`fastmap.c` `main_mem` L141-406) **[source]**

* `mem_opt_init()` (`bwamem.c` L74-110) defaults, all of which are in effect here except the three flags set on the command line:
  `a=1 b=4` (L79), `o_del=o_ins=6` (L80), `e_del=e_ins=1` (L81), `w=100` (L82), `T=30` (L83), `zdrop=100` (L84), `pen_unpaired=17` (L85),
  `pen_clip5=pen_clip3=5` (L86), `max_mem_intv=20` (L88), `min_seed_len=19` (L90), `split_width=10` (L91), `max_occ=500` (L92),
  `max_chain_gap=10000` (L93), `max_ins=10000` (L94), `mask_level=0.50f` (L95), `drop_ratio=0.50f` (L96), `XA_drop_ratio=0.80f` (L97),
  `split_factor=1.5f` (L98), `chunk_size=10000000` (L99), `max_XA_hits=5` (L101), `max_XA_hits_alt=200` (L102), `max_matesw=50` (L103),
  `mask_level_redun=0.95f` (L104), `min_chain_weight=0` (L105), `max_chain_extend=1<<30` (L106), `mapQ_coef_len=50` (float) and
  `mapQ_coef_fac = log(50)` **stored into an `int` field** (`bwamem.h` L78-79: `float mapQ_coef_len; int mapQ_coef_fac;`), so `mapQ_coef_fac == 3`
  (L107). `mat[25]` from `bwa_fill_scmat(1, 4, mat)` (L108; `bwa.c` L136-145): `mat[i*5+j] = (i==j ? 1 : -4)` for `i,j < 4`, and `-1` whenever
  either base is `N` (code 4).
* Command-line effects: `-K 100000000` -> `fixed_chunk_size = 100000000` (L193); `-Y` -> `opt->flag |= MEM_F_SOFTCLIP` (L174);
  `-R` -> `rg_line = bwa_set_rg(optarg)` (L222-223); `-t 20` -> `opt->n_threads = 20` (L168). `opt0` records nothing for `a` so
  `update_a` (L125-139, L359) is a no-op; `bwa_fill_scmat` is re-run with the same values (L360).
* No `-M` (`MEM_F_NO_MULTI` clear), no `-a`, no `-p`, no `-S`/`-P`, no `-5`/`-q`, no `-C` (`aux.copy_comment == 0`), no `-j`, no `-x`, no `-I`
  (`aux.pes0 == NULL`, so the insert-size distribution is inferred per batch, see 1.3). **[source]**
* `-R`: `bwa_set_rg` (`bwa.c` L457-488) requires the string to start with `@RG`, to contain no literal tabs, unescapes `\t` -> TAB
  (`bwa_escape` L441-455; also `\n`, `\r`, `\\`), extracts `ID:` into the global `bwa_rg_id` (= `HG002_WES.L1`), and returns the
  unescaped line. `main_mem` L261-264 appends it to `hdr_line` via `bwa_insert_header` (`bwa.c` L490-502). **[source]**

### 0.2 Index loading (`fastmap.c` L362-366, `bwa.c`) **[source]**

* `bwa_idx_load_from_shm` returns 0 unless a `/bwactl` shared-memory segment exists (`bwashm.c` L91-116); assume disk loading.
* `bwa_idx_infer_prefix(hint)` (`bwa.c` L245-269) tries `<hint>.64.bwt` first, then `<hint>.bwt`. With `hint = Homo_sapiens_assembly38.fasta.64`
  the file `Homo_sapiens_assembly38.fasta.64.64.bwt` does not exist, so the prefix is `Homo_sapiens_assembly38.fasta.64` and the files used are
  `<prefix>.bwt .sa .ann .amb .pac` (L280-284, L302, `bntseq.c` L173-176). `.64` in the file names only signals a 64-bit `bwtint_t` build
  (`bwt.h` L46 `typedef uint64_t bwtint_t`); it is not a different format.
* `bns_restore` (`bntseq.c` L168-211) reads `.alt` only if the file `<prefix>.alt` exists (L178). It does not for this reference, so every
  `anns[i].is_alt == 0`. **Consequence**: every code path conditioned on `is_alt`/`n_pri < n` is dead (listed where relevant). **[source]**
* No sequence is read from the FASTA: all reference bases come from `.pac`. See section 2 for what that implies for `N`.

### 0.3 Header output (`bwa.c` `bwa_print_sam_hdr` L407-439, `main.c` L94-96) **[source]**

Printed to stdout once, before any alignment, in this order and byte form (`\n`-terminated lines):

1. `@HD\tVN:1.5\tSO:unsorted\tGO:query` (L426; no `@HD` was given with `-H`).
2. One `@SQ\tSN:<name>\tLN:<len>` per contig in `.ann` order (L430-434); `\tAH:*` only for `is_alt` contigs, i.e. never here.
3. The `-R` line, unescaped: `@RG\tID:HG002_WES.L1\tPU:L1\tSM:HG002_HG002_WES\tLB:HG002_WES\tDS:/home/...fasta\tPL:ILLUMINA` (L437).
4. `@PG\tID:bwa\tPN:bwa\tVN:0.7.19-r1273\tCL:<argv[0]> <argv[1]> ...` - `argv` joined with single spaces exactly as received by `main`
   (`main.c` L94-96). Because the `-R` argument reaches `argv` with literal backslash-t sequences (shell quoting), `CL` contains `\t` as two
   characters, never a TAB. **[inferred: depends on how the container invoked bwa; copy the literal line from the golden]**

### 0.4 Pipeline and thread-independence (`fastmap.c` `process` L64-123, L394-395; `kthread.c` `kt_pipeline`, `kt_for`) **[source]**

* `kt_pipeline(2, process, &aux, 3)`: a 3-step pipeline (read batch / align batch / write batch) with two worker slots. `ktp_worker`
  only lets a worker run step `s` of batch `n` when no worker is at a step `<= s` for a batch `< n`, so step 2 (output, L113-121) executes in
  strictly increasing batch order and, inside a batch, in `seqs[]` order (`err_fputs(seqs[i].sam)` for `i = 0..n_seqs-1`). Output order is
  therefore: batch order, then read order within the batch, then for each read the lines of its own SAM string (see 8.3). **[source]**
* `mem_process_seqs` (`bwamem.c` L1235-1264): `kt_for(n_threads, worker1, ...)` over pair index `i` (L1252) writes `w->regs[2i], regs[2i+1]`;
  then `mem_pestat` over the **whole batch** (L1256-1259); then `kt_for(worker2)` (L1260) which writes `seqs[2i].sam`, `seqs[2i+1].sam`.
  Each work item is a pure function of (opt, index, batch inputs, `n_processed`, pes); the per-thread `smem_aux_t` scratch (L1249-1251,
  `tid`-indexed) carries no state between items (every use starts with `a->mem.n = 0`, L145). Hence output does not depend on `-t`, matching
  the observation that two `-t 20` runs were byte-identical. **[source + inferred]**
* `n_processed` (`fastmap.c` L111) is the number of reads (not pairs) processed before the batch; it seeds the hash tie-breakers (6.1, 7.3).

---

## 1. Input batching and per-batch state

### 1.1 FASTQ parsing and read names (`bwa.c` L54-112, `kseq.h`) **[source]**

* `kseq_read` (`kseq.h` `__KSEQ_READ`): the name is everything after `@` up to the first `isspace()` character (`ks_getuntil(ks, KS_SEP_SPACE=0, ...)`,
  `kseq.h` L39, L186); the rest of the header line is the comment (L187; a trailing `\r` is stripped). Sequence lines up to the `+` line are
  concatenated; the quality string must have the same length. For `@SRR14724482.33 33/1` the name is `SRR14724482.33` and the comment `33/1`.
* `trim_readno` (`bwa.c` L54-58) strips a trailing `/<digit>` **from the name only** (`s->l > 2 && s[l-2]=='/' && isdigit(s[l-1])`).
  `SRR14724482.33` has no such suffix, so QNAME is `SRR14724482.33` unchanged. **[source]** (The `/1` in `33/1` lives in the comment.)
* `kseq2bseq1` (L70-77) copies name, comment, seq, qual (`l_seq = seq.l`). Because `-C` is absent, `process` step 0 frees every comment
  (`fastmap.c` L78-82), so comments never reach the SAM output.
* Read bases are converted to 2-bit codes in place by `mem_align1_core` (`bwamem.c` L1087-1088) through `nst_nt4_table`
  (`bntseq.c` L46-63): `A/a=0 C/c=1 G/g=2 T/t=3`, `-`=5, everything else (including `N` and IUPAC codes) = 4. SAM SEQ is later regenerated from
  these codes (8.3), so lowercase is uppercased and any non-ACGT letter becomes `N`. **[source]**
* Qualities are passed through verbatim (bytes), never reinterpreted.

### 1.2 Batch formation (`bwa.c` `bseq_read` L79-112; `fastmap.c` L73, L394) **[source]**

* `actual_chunk_size = 100000000` (`-K`). `bseq_read` reads one record from R1 and one from R2 per iteration, appends them in the order
  R1, R2 (ids `2k`, `2k+1`), adds **both** `l_seq` values to `size`, and stops after the pair at which `size >= chunk_size` (`&& (n&1)==0`, always
  true here) (L85-105). So a batch is the shortest prefix of remaining pairs whose total base count (R1+R2 lengths, post-fastp lengths as they
  are in the file) reaches 100,000,000; the pair that crosses the threshold is included. The last batch is whatever remains. **[source]**
* If R2 runs out first the batch is truncated with a warning (L86-89); the final `size == 0` check only warns (L106-109). fastp emits equal
  counts, so irrelevant. **[inferred]**
* `-p` is not set, so `bseq_classify` is never used; `MEM_F_PE` is set because a second FASTQ is given (`fastmap.c` L378-391).

### 1.3 Per-batch insert-size estimation (`bwamem_pair.c` `mem_pestat` L72-135) **[source]**

Run once per batch on the `worker1` results (regions sorted as in 5.5, i.e. `a[0]` is the best region of each read), producing `pes[4]` for
orientations indexed `d = 0 FF, 1 FR, 2 RF, 3 RR` (the stderr labels use `"FR"[d>>1&1], "FR"[d&1]`, L97).

Candidate collection (L78-90): for each pair `i`:
1. skip if either read has no region (L84);
2. skip if `cal_sub(r) > 0.8 * r->a[0].score` for either read (L85-86; `MIN_RATIO = 0.8` double, L42). `cal_sub` (L58-70): the score of the first
   region `j >= 1` whose query interval overlaps `a[0]`'s by at least `min_l * mask_level` (`int * float` -> float compare), where
   `min_l` = shorter of the two query spans; if none, `min_seed_len * a = 19`. (So a read whose best hit has no competing hit is "unique"
   iff `19 <= 0.8 * score`, i.e. `score >= 24`.)
3. skip unless `r[0]->a[0].rid == r[1]->a[0].rid` (L87);
4. `dir = mem_infer_dir(l_pac, r0.rb, r1.rb, &is)` (L88; L49-56): `r1 = (b1 >= l_pac)`, `r2 = (b2 >= l_pac)`, `p2 = (r1 == r2) ? b2 : 2*l_pac - 1 - b2`,
   `*dist = |p2 - b1|`, `dir = (r1 == r2 ? 0 : 1) ^ (p2 > b1 ? 0 : 3)`. Note `dist` is the difference of the two region *start* coordinates
   after mapping read 2 onto read 1's strand; for an FR pair this is `(forward coordinate of read 2's last aligned base) - (first aligned base of
   read 1)`, i.e. **outer span minus one**, not the SAM TLEN. **[source]**
5. push `is` into `isize[dir]` if `0 < is <= max_ins (10000)` (L89).

Per orientation (L92-127): if fewer than `MIN_DIR_CNT = 10` candidates -> `failed = 1` (L96-100). Otherwise sort ascending (`ks_introsort_64`,
L102) and, with `n = q->n`:
```
p25 = q[(int)(.25*n + .499)]; p50 = q[(int)(.50*n + .499)]; p75 = q[(int)(.75*n + .499)]     // double arithmetic, truncation
low  = (int)(p25 - 2.0*(p75-p25) + .499); if (low < 1) low = 1
high = (int)(p75 + 2.0*(p75-p25) + .499)
avg  = (sum of q[i] with low <= q[i] <= high, in sorted order, accumulated in double) / x       // x = count
std  = sqrt( (sum over the same i of (q[i]-avg)*(q[i]-avg)) / x )                              // double; q[i]-avg is double
low  = (int)(p25 - 3.0*(p75-p25) + .499); high = (int)(p75 + 3.0*(p75-p25) + .499)
if (low  > avg - 4.0*std) low  = (int)(avg - 4.0*std + .499)
if (high < avg + 4.0*std) high = (int)(avg + 4.0*std + .499)
if (low < 1) low = 1
```
(L103-124; constants L45-47.) `q[i]` are `uint64_t`; `q[i] - r->avg` is computed in double. All casts `(int)` truncate toward zero.
Finally (L128-134) any orientation with `n < max_n * 0.05` (`MIN_DIR_RATIO`, double) is marked `failed` even if it had >= 10 pairs. `avg`,
`std`, `low`, `high`, `failed` of all four orientations are inputs to pairing (7.2), mate rescue (7.1) and the proper-pair flag (7.4). **[source]**

**Small-batch consequence [source]**: a batch with fewer than 10 unique FR candidates (typically the short final batch of a chunk, or any batch of
`< 10` pairs) has every orientation `failed`. Then `mem_matesw` skips all orientations (7.1 step 1), `mem_pair` pushes nothing (7.2 step 2 skips failed
`dir`s) so `o == 0` and every pair takes `no_pairing` (7.3), and the proper-pair test there fails (`!pes[d].failed`): no `0x2`, no mate rescue, SE-only
MAPQ for the whole batch. Reproducing batch boundaries exactly (1.2) is therefore load-bearing even for the pairing fields.

`memset(pes, 0, ...)` first (L76), so a failed orientation has `avg = std = 0`, `low = high = 0` unless set above. `failed` orientations
are never used for distance comparisons in `mem_pair`/`mem_sam_pe` (they are checked first), but **are** used in `mem_matesw`'s "orientation
already found" test only via `skip[r] = 1` (7.1). **[source]**

---

## 2. Index contents and reference representation

### 2.1 `.bwt` / `.sa` (`bwt.c`, `bwt.h`) **[source]**

* `.bwt` (`bwt_restore_bwt` L443-462): `primary` (u64), `L2[1..4]` (4 x u64), then the BWT with interleaved occurrence counts.
  `seq_len = L2[4]` = `2*l_pac`. Layout (`bwt.h` L37-39, L74-75): per block of `OCC_INTERVAL = 128` symbols, 16 `uint32`: 4 x `uint64` cumulative
  counts `Occ(c, block start)` for `c = 0..3`, then 8 `uint32` holding 128 2-bit symbols, 16 per word, most significant pair first
  (`bwt_B0(b,k) = word >> ((~k & 0xf) << 1) & 3`, L80). The `$` symbol is not stored; `bwt_occ`/`bwt_2occ4` apply `k -= (k >= primary)` (L114, L177,
  L192-201). `bwt_occ4` (L169-186) and `bwt_2occ4` (L189-220) count via `cnt_table` (L42-51); `bwt_occ` via `__occ_aux` (L98-105). Semantics are
  the standard `Occ(c, k)` = number of `c` in `B[0..k]` and must be reproduced exactly (including `k == (bwtint_t)-1 -> 0`, L113, L173-176).
* `.sa` (`bwt_restore_sa` L421-441): `primary`, 4 skipped words, `sa_intv` (= 32 for `bwa index`, `bwtindex.c` L316), `seq_len`, then
  `n_sa - 1 = (seq_len + 32)/32 - 1` sampled values; `sa[0]` is set to `-1` in memory (L437). `bwt_sa(k)` (L86-96): walk `k = bwt_invPsi(k)` while
  `k & (sa_intv-1)`, counting steps, return `steps + sa[k/sa_intv]`. `bwt_invPsi` (L53-59): `x = k - (k > primary); x = B0(x); x = L2[x] + Occ(x, k)`;
  returns `0` if `k == primary`. Positions returned are in the concatenated `[0, 2*l_pac)` space (forward then reverse strand). **[source]**
* Bidirectional intervals `bwtintv_t {x[0], x[1], x[2], info}` (`bwt.h` L62-64): `x[0]` = SA interval start of the match, `x[1]` = start of the
  interval of its reverse complement, `x[2]` = size, `info = qbeg<<32 | qend` (end exclusive). `bwt_set_intv(c)` (L82): `x[0] = L2[c]+1`,
  `x[2] = L2[c+1]-L2[c]`, `x[1] = L2[3-c]+1`. `bwt_extend` (`bwt.c` L262-275) is the standard backward/forward extension; port verbatim, including the
  `primary` adjustment on L271.

### 2.2 `.ann` / `.amb` / `.pac` (`bntseq.c`) **[source]**

* `.ann` (L107-137): `l_pac n_seqs seed`, then per contig `gi name [anno]` and `offset len n_ambs`. Contig `i` occupies forward positions
  `[offset_i, offset_i + len_i)`; `bns_pos2rid` (L354-368) is a binary search over offsets. Contig boundaries are enforced by `bns_intv2rid`
  (L370-378: `-2` if the interval straddles `l_pac`, `-1` if start and end-1 fall in different contigs) and by `bns_fetch_seq` (L426-451), which clips
  a requested window to the contig containing `mid` (on the same strand).
* `.pac` (`bwa.c` L307-312): `l_pac/4 + 1` bytes of 2-bit forward-strand bases, 4 per byte, most significant pair first (`_get_pac`, `bntseq.c` L230).
  Reverse-strand bases are derived by complementing and reading backwards (`bns_get_seq` L403-424).
* **Reference `N`s**: when the index was built, every non-ACGT base was written to `.pac` as `lrand48() & 3` with `srand48(11)` (L266, L295-296),
  and the run was recorded in `.amb` (L248-263). `bwa mem` never consults `.amb` (`bns_cnt_ambi` L380-401 is only used by `fastmap`); it aligns
  against the random bases stored in `.pac`. A port must therefore read `.pac` as-is; it must not regenerate the reference from the FASTA and must not
  treat holes specially. Reads overlapping an `N` run can align across it (scoring the random bases as ordinary matches/mismatches), and `MD`
  will print the random base letters. **[source]**
* `.amb` is parsed (L138-155) only to validate `l_pac`/`n_seqs` consistency.

---

## 3. Seeding (`bwamem.c` `mem_collect_intv` L140-188; `bwt.c`)

Input: the read as codes 0-4, `len`. Output `a->mem`: a list of `bwtintv_t` sorted by `info` (`ks_introsort(mem_intv)`, L187; comparator
`(a).info < (b).info`, L116). **[source]**

1. **Pass 1 - SMEMs** (L146-157): `x = 0; while (x < len) { if (seq[x] < 4) { x = bwt_smem1(bwt, len, seq, x, 1, &mem1, tmpv); keep every interval
   with length >= 19 } else ++x; }`. `bwt_smem1 = bwt_smem1a(..., min_intv=1, max_intv=0, ...)` (L353-356). `bwt_smem1a` (L289-351):
   * returns `x+1` immediately if `q[x] > 3` (L296);
   * forward phase (L301-321): start with the single-base interval at `x`; for `i = x+1 .. len-1`: `max_intv = 0` so the first branch never fires;
     if `q[i] < 4`, extend with `c = 3 - q[i]` (`bwt_extend(..., is_back=0)`); if the new size differs from the old, push the old interval; if the new
     size `< min_intv (1)` (i.e. 0) stop; else `ik = ok[c]; ik.info = i+1`. An ambiguous base pushes `ik` and stops (L316-319). If the loop reaches
     `len`, push the last `ik` (L321). Reverse `curr` so the longest match is first (L322); the return value is `curr->a[0].info` = end of the
     longest exact match starting at `x` (L323) - this is where pass 1 resumes.
   * backward phase (L326-345): for `i = x-1 .. -1`: `c = (i < 0 || q[i] > 3) ? -1 : q[i]`; for each interval `p` in `prev` (longest first):
     extend backward if `c >= 0` (the `ik.x[2] >= max_intv` test is always true since `max_intv = 0`); if `c < 0` or the extension size
     `< min_intv`, the interval is a MEM: push to `mem` **only if `curr->n == 0`** (no longer match survived this step) **and** it is not contained
     in the previously pushed one (`i + 1 < mem[last].info>>32`), with `info |= (i+1)<<32` (L331-337); otherwise push the extension to `curr`
     unless its size equals the last pushed size (L338-341). Stop when `curr` is empty. Finally reverse `mem` (sorted by start, L346).
   * Net effect: all SMEMs covering position `x`, each reported once with `info = (qbeg<<32 | qend)`; then pass 1 jumps to the end of the longest
     one. Port the loop verbatim rather than re-deriving SMEM semantics - the `curr->n == 0` and size-dedup conditions matter.
2. **Pass 2 - re-seeding** (L159-168): for each pass-1 SMEM (in the order they were collected, `old_n` fixed before the loop) with
   `end - start >= split_len` and `x[2] <= split_width (10)`, where `split_len = (int)(19 * 1.5f + .499) = 28` (L144: `int * float` -> `28.5f`, then
   `+ .499` in double -> `28.999` -> 28), call `bwt_smem1(..., x = (start+end)>>1, min_intv = x[2] + 1, ...)` and keep intervals of length >= 19.
   These are SMEMs through the midpoint whose occurrence exceeds the parent's.
3. **Pass 3 - LAST-like seeds** (L170-185): `max_mem_intv = 20 > 0`: `x = 0; while (x < len) { if (seq[x] < 4) { x = bwt_seed_strategy1(bwt, len, seq, x,
   19, 20, &m); if (m.x[2] > 0) push m; } else ++x; }`. `bwt_seed_strategy1` (`bwt.c` L358-379): forward-extend from `x`; at the first `i` where the
   interval size drops **below 20 and `i - x >= 19`** (so the match length `i - x + 1 >= 20`), report `info = x<<32 | (i+1)` and return `i+1`; an
   ambiguous base returns `i+1` with no seed (`memset` to zero, L363, so `x[2] == 0`); reaching the end returns `len` with no seed. **[source]**
4. Sort by `info` (L187). Entries with equal `info` are byte-identical intervals (same substring -> same SA interval), so instability is harmless.
   **[inferred]**

Ambiguous read bases (code 4) therefore terminate every exact-match extension and are never part of a seed.

---

## 4. Chaining (`bwamem.c` L194-411)

### 4.1 Seed instantiation and B-tree chaining (`mem_chain` L277-341) **[source]**

* If `len < 19` return no chains (L286) -> read is unmapped.
* `frac_rep` (L291-298): over the sorted SMEM list, take only intervals with `x[2] > max_occ (500)`; `l_rep` = total length of the union of their
  query intervals (greedy merge in sorted order); `frac_rep = (float)l_rep / len` (float division, L336) attached to every chain.
* For each interval `p` (sorted order) (L299-327): `slen = qend - qbeg`; `step = x[2] > 500 ? x[2] / 500 : 1` (integer division);
  for `k = 0; k < x[2] && count < 500; k += step, ++count`: `rbeg = bwt_sa(bwt, x[0] + k)`, seed `s = {rbeg, qbeg, len = slen, score = slen}`;
  `rid = bns_intv2rid(bns, rbeg, rbeg + len)`; **skip the seed if `rid < 0`** (bridging contigs or the strand boundary, L312-313).
* Chains live in a `kbtree` keyed by `pos` = `rbeg` of the chain's first seed (`chain_cmp` L212). For each seed: if the tree is non-empty,
  `kb_intervalp` (`kbtree.h` L151-172) finds `lower` = the chain with the greatest key `<= rbeg` (exact match if present; L117-131 `__kb_getp_aux`);
  if no `lower` or `test_and_merge(lower)` returns 0, create a new chain `{n=1, m=4, seeds=[s], rid, is_alt=0}` and `kb_putp` it (L314-325).
* `test_and_merge` (L216-237), with `last` = the chain's last seed: `0` if `seed_rid != c->rid`; `1` (absorbed, nothing stored) if the seed is contained
  in `[seeds[0].qbeg, last.qend)` x `[seeds[0].rbeg, last.rend)`; `0` if the chain is on the forward strand and the seed on the reverse
  (`(last.rbeg < l_pac || seeds[0].rbeg < l_pac) && p.rbeg >= l_pac`); else `x = p.qbeg - last.qbeg` (non-negative because seeds arrive in
  `qbeg` order), `y = p.rbeg - last.rbeg`; append iff `y >= 0 && x - y <= 100 && y - x <= 100 && x - last.len < 10000 && y - last.len < 10000`.
* After all seeds, `__kb_traverse` (`kbtree.h` L336-360) emits chains in key order (in-order traversal) into `chain.a` (L330-334). **[source]**
* Exactness note **[inferred]**: `kb_intervalp` is a deterministic function of the key set only when keys are distinct. Two chains can share a
  `pos` (same `rbeg`, seed not mergeable because `x > w`); then which duplicate is returned, and the traversal order among duplicates, depend
  on the B-tree shape. Port `kbtree.h` (`__KB_INIT` node sizing with `KB_DEFAULT_SIZE`, `__kb_getp_aux`, `kb_intervalp`, `__kb_split`,
  `__kb_putp_aux`, `__kb_traverse`) verbatim rather than substituting an ordered map. On x86-64 `sizeof(mem_chain_t) = 40`, `sizeof(void*) = 8`,
  so `t = ((512 - 4 - 8) / (8 + 40) + 1) >> 1 = 5` and nodes hold at most 9 keys (`kbtree.h` L58-63).

### 4.2 Chain filtering (`mem_chain_flt` L353-411) **[source]**

1. `w = mem_chain_weight(c)` (L239-258): min over query and reference of the union length of seed intervals, scanned in seed order with a running
   `end`. Chains with `w < min_chain_weight (0)` are dropped (none).
2. `ks_introsort(mem_flt)`: descending `w` (L350-351, L367). **Unstable**; ties must follow ksort's exact introsort (see 12.2).
3. Greedy overlap filter (L369-393), `chains` = kept indices starting with `[0]`, `a[0].kept = 3`. For each `i >= 1`, for each kept `j` in order:
   query overlap `e_min - b_max > 0` (`chn_beg/chn_end` L347-348: first seed's `qbeg`, last seed's `qend`) and (`!a[j].is_alt || a[i].is_alt`, always
   true); significant if `e_min - b_max >= min_l * 0.5f` (`int * float`) and `min_l < 10000`: then `large_ovlp = 1`, `a[j].first = i` if unset, and
   **drop `i`** (break) if `a[i].w < a[j].w * 0.5f && a[j].w - a[i].w >= 38`. If no kept chain dropped it, keep it with `kept = large_ovlp ? 2 : 3`.
4. `a[c->first].kept = 1` for every kept chain's first shadowed chain (L394-397) - these are extended only to improve mapq estimates.
5. `max_chain_extend = 1<<30` so the cap loop (L399-404) never triggers. Chains with `kept == 0` are freed; survivors are compacted in place
   (sorted order preserved) (L405-410).

### 4.3 `mem_flt_chained_seeds` (L624-641) **[source]**

`min_l = 5.5f * log(l_query)` (double); returns immediately if `min_l > 0.05f * l_query` (float product, promoted). For `l_query = 151`:
`27.6 > 7.55`, so **this function is a no-op for all reads shorter than ~700 bp** and seed scores stay `= len`. (It becomes active for reads
where `5.5 ln(l) <= 0.05 l`.)

---

## 5. Extension (`bwamem.c` `mem_chain2aln` L658-812; `ksw.c` `ksw_extend2` L416-515)

### 5.1 Reference window and seed order **[source]**

* `cal_max_gap(qlen)` (L647-654) `= clamp(max((qlen - 6)/1 + 1, (qlen - 6)/1 + 1), 1, 200)` = `min(max(qlen - 5, 1), 200)` (double division then
  truncation; exact for these integers).
* `rmax[0] = min over seeds of rbeg - (qbeg + cal_max_gap(qbeg))`, `rmax[1] = max of rbeg + len + ((l_query - qbeg - len) + cal_max_gap(l_query - qbeg - len))`,
  clamped to `[0, 2*l_pac]`; if the window straddles `l_pac` it is cut at `l_pac` on the chain's strand (L668-683). `rseq = bns_fetch_seq(bns, pac,
  &rmax[0], seeds[0].rbeg, &rmax[1], &rid)` further clips to the contig (L685); `assert(c->rid == rid)`.
* Seeds are processed in descending `(score << 32 | index)` order (`srt`, L688-691, `ks_introsort_64`) - with scores all `= len` here this is
  longest-first, ties by **higher** original index first (`k` runs from `c->n - 1` downwards).

### 5.2 Skip test for seeds contained in an existing region (L697-732) **[source]**

For each seed `s` and each region `p` already in `av` (in insertion order): skip `p` unless `s` is fully inside `p` on both query and reference and
`s->len - p->seedlen0 <= 0.1 * l_query` (double). Compute `qd = s.qbeg - p.qb`, `rd = s.rbeg - p.rb`, `w = min(cal_max_gap(min(qd, rd)), p.w)`;
if `qd - rd < w && rd - qd < w` the seed is "around" `p` (break); else the same test on the trailing side. If some `p` matched: scan later seeds
`t` in `srt` order `i = k+1 .. c->n-1` (skipping `srt[i] == 0` entries) with `t.len >= s.len * .95` (double); if any such `t` overlaps `s` on the query
by at least `s.len >> 2` **and lies on a different diagonal**, extension still happens; otherwise mark `srt[k] = 0` and skip. **[source]**

### 5.3 Left/right extension **[source]**

New region `a` is zeroed, `a->w = aw[0] = aw[1] = 100`, `score = truesc = -1`, `rid = c->rid` (L734-738).

Left (L741-770, only if `s->qbeg > 0`): `qs` = `query[0..qbeg)` reversed, `rs` = `rseq[0..s.rbeg - rmax[0])` reversed. Up to `MAX_BAND_TRY = 2`
attempts with `aw[0] = 100 << i`: `a->score = ksw_extend2(qbeg, qs, tmp, rs, 5, mat, 6, 1, 6, 1, aw[0], pen_clip5 = 5, zdrop = 100, h0 = s->len * 1,
&qle, &tle, &gtle, &gscore, &max_off[0])`; stop retrying if the score did not change or `max_off < (aw>>1) + (aw>>2)` (75 for 100). Then:
if `gscore <= 0 || gscore <= a->score - 5` -> local: `qb = qbeg - qle`, `rb = s.rbeg - tle`, `truesc = score`; else to-end: `qb = 0`,
`rb = s.rbeg - gtle`, `truesc = gscore`. Without a left part: `score = truesc = len`, `qb = 0`, `rb = s.rbeg` (L770).

Right (L772-797, only if `qbeg + len != l_query`): `qe = qbeg + len`, `re = s.rbeg + len - rmax[0]`, `sc0 = a->score`; same retry loop with
`ksw_extend2(l_query - qe, query + qe, rmax[1] - rmax[0] - re, rseq + re, ..., aw[1], pen_clip3 = 5, 100, sc0, ...)`; local: `a->qe = qe + qle`,
`a->re = rmax[0] + re + tle`, `truesc += score - sc0`; to-end: `a->qe = l_query`, `a->re = rmax[0] + re + gtle`, `truesc += gscore - sc0`.
Without a right part: `qe = l_query`, `re = s.rbeg + len`.

Then `seedcov` = total length of chain seeds fully inside `[qb,qe) x [rb,re)` (L801-805), `a->w = max(aw[0], aw[1])`, `seedlen0 = s->len`,
`frac_rep = c->frac_rep` (L806-809). Note `a->score` is the best **local** extension score (this becomes `AS`), while `truesc` is the score of
the chosen end points (used for band inference and in `mem_patch_reg`). **[source]**

### 5.4 `ksw_extend2` exact semantics (`ksw.c` L416-515) **[source]**

Scalar affine-gap extension; port verbatim. Key points:
* Query profile `qp[k*qlen + j] = mat[k*5 + query[j]]` (L426-429); `N` (code 4) scores `-1` against anything.
* First row: `eh[0].h = h0`, `eh[1].h = max(h0 - 7, 0)`, `eh[j].h = eh[j-1].h - 1` while `eh[j-1].h > 1` (L431-433); `e` all 0.
* Band: `w = min(w, max_ins, max_del)` with `max_ins = max((int)((qlen*1 + 5 - 6)/1. + 1.), 1) = max(qlen, 1)` (L435-443).
* Row `i` (target base): `beg = max(beg, i - w)`, `end = min(end, i + w + 1, qlen)`; `h1 = beg == 0 ? max(h0 - (6 + 1*(i+1)), 0) : 0` (L452-459).
  Per cell (L460-484): `M = eh[j].h; e = eh[j].e; eh[j].h = h1; M = M ? M + q[j] : 0` (**an H of 0 cannot be extended by a match**);
  `h = max(M, e, f); h1 = h; mj = (m > h) ? mj : j; m = max(m, h)` (row max, **rightmost** column on ties); `e = max(e - 1, max(M - 7, 0))` stored;
  `f = max(f - 1, max(M - 7, 0))`. After the row: `eh[end].h = h1; eh[end].e = 0`; if `end == qlen`: `max_ie = (gscore > h1) ? max_ie : i; gscore = max(gscore, h1)`
  (ties update `max_ie` to the later row) (L485-489). `if (m == 0) break` (L490). If `m > max`: `max = m, max_i = i, max_j = mj`,
  `max_off = max(max_off, |mj - i|)`; else z-drop (L494-500): with `di = i - max_i`, `dj = mj - max_j`: if `di > dj` break when
  `max - m - (di - dj)*e_del > 100`, else break when `max - m - (dj - di)*e_ins > 100`. Band shrink (L502-505): `beg` = first `j` in `[beg,end)` with
  `eh[j].h != 0 || eh[j].e != 0`; `end` = `min(j + 2, qlen)` where `j` is the last index `>= beg` with a non-zero cell (scanning from `end` downwards).
* Returns `max` (initialised to `h0`, L445), `qle = max_j + 1`, `tle = max_i + 1`, `gtle = max_ie + 1`, `gscore` (initial `-1`), `max_off`.

### 5.5 Redundancy removal and patching (`mem_sort_dedup_patch` L463-515, `mem_patch_reg` L432-461) **[source]**

Called once per read on all regions from all chains (`mem_align1_core` L1103) and after every mate rescue (7.1, with `bns = pac = query = 0`).

1. `ks_introsort(mem_ars2)`: ascending `re` only (L417, L467) - **many ties; unstable; port ksort exactly** (12.2). Set `n_comp = 1` on all.
2. For `i = 1..n-1`, `p = a[i]`; skip unless `p.rid == a[i-1].rid && p.rb < a[i-1].re + 10000`. For `j = i-1` down while `a[j].rid == p.rid && p.rb < a[j].re + 10000`:
   `q = a[j]`; skip if `q.qe == q.qb` (already excluded). `or = q.re - p.rb`, `oq = (q.qb < p.qb) ? q.qe - p.qb : p.qe - q.qb`, `mr = min ref span`,
   `mq = min query span`. If `or > 0.95f * mr && oq > 0.95f * mq` (float products; `int64` compared as float): if `p.score < q.score` exclude `p`
   (`p.qe = p.qb`) and stop the inner loop, else exclude `q` (ties exclude `q`, the one with the smaller `re` / earlier in sorted order).
   Else if `q.rb < p.rb` and `mem_patch_reg(q, p) > 0` (only when `bns != 0`): merge `q` into `p`: `n_comp += q.n_comp + 1`, `seedcov = max`,
   `sub = max`, `csub = max`, `qb = q.qb`, `rb = q.rb`, `truesc = score = patched score`, `w = patched w`, exclude `q` (L486-495).
3. Compact (L498-503); `ks_introsort(mem_ars)`: `score` desc, then `rb` asc, then `qb` asc (L420, L504); drop any hit whose `(score, rb, qb)` equals its
   predecessor's (L505-513).

`mem_patch_reg(a=q, b=p)` (L432-461): returns 0 if on different strands, not colinear (`a.qb >= b.qb || a.qe >= b.qe || a.re >= b.re`), or the
bandwidth tests fail: `w = |(a.re - b.rb) - (a.qe - b.qb)|`, `r = |(a.re - b.rb)/(double)(b.re - a.rb) - (a.qe - b.qb)/(double)(b.qe - a.qb)|`; if the
two hits do not overlap (`a.re < b.rb || a.qe < b.qb`): reject if `w > 200 || r >= 0.05f` (float literal promoted: `0.05000000074505806`); else reject
if `w > 400 || r >= 0.1f` (`0.10000000149011612`). Then `w = min(w + a.w + b.w, 400)` and `bwa_gen_cigar2(..., w, ..., b.qe - a.qb, query + a.qb, a.rb, b.re,
&score, 0, 0)` computes a global score (8.1); `q_s = (int)((double)(b.qe - a.qb) / ((b.qe - b.qb) + (a.qe - a.qb)) * (b.score + a.score) + .499)`,
`r_s` analogously on reference spans; accept iff `(double)score / max(q_s, r_s) >= 0.9f` (`0.8999999761581421`). **[source]**

### 5.6 `mem_align1_core` (L1081-1117) **[source]**

Sequence: convert to codes -> `mem_chain` -> `mem_chain_flt` -> `mem_flt_chained_seeds` -> `mem_chain2aln` for each surviving chain in order,
appending to one `regs` vector -> `mem_sort_dedup_patch` -> `is_alt` tagging (never). The result `regs` is sorted by `(score desc, rb, qb)`.

---

## 6. Single-end post-processing

### 6.1 `mem_mark_primary_se` (L547-584) and `_core` (L519-545) **[source]**

Called with `id = (n_processed>>1 + pair_index) << 1 | read_index` from `mem_sam_pe` (`bwamem_pair.c` L303-304, `bwamem.c` L1230).

1. For every region `i` (`i` = array position at call time, i.e. the `(score desc, rb, qb)` order left by the last `mem_sort_dedup_patch` **after** mate
   rescue has interleaved rescued regions, 7.1 step 5): `sub = alt_sc = 0`, `secondary = secondary_all = -1`, `hash = hash_64(id + i)` (`utils.h` L98-109:
   Thomas Wang 64-bit mix, port verbatim; wrapping u64 arithmetic). `n_pri = n` (nothing is alt).
2. `ks_introsort(mem_ars_hash)`: `score` desc, `is_alt` asc, **`hash` asc** (L423, L556). Hashes are distinct unless `hash_64` collides, so this order
   is effectively total: equal-score hits are ordered by the hash of (`id` + original index).
3. `mem_mark_primary_se_core` (L519-545): `tmp = max(a + b, o_del + e_del, o_ins + e_ins) = 7`. `z = [0]`. For `i = 1..n-1`: for each `j` in `z` (in order):
   query overlap `e_min - b_max > 0` and `>= min_l * 0.5f` (float) -> `if (a[j].sub == 0) a[j].sub = a[i].score`; if `a[j].score - a[i].score <= 7`
   (and the alt condition, always true) `++a[j].sub_n`; **break** (first overlapping primary wins). If no primary overlapped, push `i` to `z`
   (it is primary); else `a[i].secondary = z[k]`.
4. `secondary_all = i` then, because `n_pri == n`, the `else` branch (L578-581) sets `secondary_all = secondary` for all. The ALT branch (L564-577)
   is dead. Returns `n_pri = n`.

Outputs used later: `secondary` (`-1` primary, else index of the shadowing primary), `sub` (score of the first overlapping lower hit, 0 if none),
`sub_n` (count of overlapping hits within 7 of this score), and the array order.

### 6.2 `mem_approx_mapq_se` (L982-1006) **[source]**

```
sub = a->sub ? a->sub : 19;  sub = max(a->csub, sub)
if (sub >= a->score) return 0
l = max(qe - qb, re - rb)
identity = 1. - (double)(l*1 - score) / (1 + 4) / l                 // double: ((l - score) / 5.0) / l
if (score == 0) mapq = 0
else { tmp = (l < 50) ? 1. : 3 / log(l);                            // mapQ_coef_fac is the int 3; log() is libm double
       tmp *= identity * identity;
       mapq = (int)(6.02 * (score - sub) / 1 * tmp * tmp + .499); }  // evaluated left to right: ((6.02*(score-sub))/1)*tmp*tmp + .499
if (sub_n > 0) mapq -= (int)(4.343 * log(sub_n + 1) + .499)
mapq = min(mapq, 60); mapq = max(mapq, 0)
mapq = (int)(mapq * (1. - frac_rep) + .499)                          // frac_rep is float, promoted exactly
```
`mapQ_coef_len = 50 > 0` so the `MEM_MAPQ_COEF` branch (L998-999) is dead. **[source]**

### 6.3 `XS` **[source]**

`mem_reg2aln` L1185: `a.sub = max(ar->sub, ar->csub)`; printed as `XS:i:<sub>` whenever `sub >= 0` (`mem_aln2sam` L937), i.e. for every
non-secondary record including unmapped ones (`sub == 0`, see 8.3).

---

## 7. Pairing (`bwamem_pair.c`)

### 7.1 Mate rescue (`mem_sam_pe` L291-302, `mem_matesw` L137-206) **[source]**

Before anything else (no `-S`): for each read `i`, `b[i]` = copy of regions with `score >= a[i].a[0].score - 17` (L294-297, in array order); then
for `i = 0, 1` and `j < min(b[i].n, 50)`: `mem_matesw(opt, bns, pac, pes, &b[i].a[j], l_seq of the mate, mate codes, &a[!i])` (L298-300). The
anchor list for read 1 is snapshotted before read 0's rescues modify `a[1]`, but each rescue sees the current `a[!i]`.

`mem_matesw` (L137-206):
1. `skip[r] = pes[r].failed`. For every existing region `m` of the mate: `r = mem_infer_dir(l_pac, a->rb, m.rb, &dist)`; if
   `pes[r].low <= dist <= pes[r].high` then `skip[r] = 1` (a consistent mate already exists; **this uses `low/high` even of failed orientations**,
   which are 0, so a failed orientation is skipped anyway). If all four are skipped, return 0 (L142-150).
2. For each non-skipped `r`: `is_rev = (r>>1) != (r&1)`, `is_larger = !(r>>1)`. Window in concatenated coordinates (L163-171):
   forward mate: `rb = is_larger ? a.rb + low : a.rb - high`, `re = (is_larger ? a.rb + high : a.rb - low) + l_ms`;
   reverse-complemented mate: `rb = (is_larger ? a.rb + low : a.rb - high) - l_ms`, `re = is_larger ? a.rb + high : a.rb - low`; clamp to `[0, 2*l_pac]`.
   `ref = bns_fetch_seq(bns, pac, &rb, (rb+re)>>1, &re, &rid)` (clips to the contig of the midpoint). Require `rid == a->rid && re - rb >= 19`.
3. `ksw_align2(l_ms, seq, re - rb, ref, 5, mat, 6, 1, 6, 1, xtra, 0)` with `xtra = KSW_XSUBO | KSW_XSTART | (l_ms * 1 < 250 ? KSW_XBYTE : 0) | 19`
   (L176-177). For reads shorter than 250 bp this selects the **8-bit** SSE2 kernel `ksw_u8` (`ksw.c` L379-401, L122-253); longer reads use `ksw_i16`.
   Semantics (port as scalar Smith-Waterman that reproduces them exactly; **do not reuse the `ksw_extend2`/`ksw_global2` recurrence**):
   * local alignment, scores clamped at 0, `oe = 7`, `e = 1`; `minsc = 19` (XSUBO); profile shift `q->shift = 256 - (-4) = 4` (u8), `q->max = 1`.
   * recurrence (`ksw.c` L181-196 / L306-318): `H(i,j) = max(H(i-1,j-1) + S, E(i,j), F(i,j))` clamped at 0; then **gap opening is taken from `H`**, not from
     the match-only `M` used by `ksw_extend2`/`ksw_global2`: `E(i+1,j) = max(E(i,j) - 1, H(i,j) - 7)`, `F(i,j+1) = max(F(i,j) - 1, H(i,j) - 7)` (unsigned
     saturating, so clamped at 0). The "lazy-F" loop (L201-211 / L321-331) only completes `F` propagation; the converged values are those of this recurrence.
   * striped padding (`ksw_qinit` L95-103, L104-112): the query is laid out in `slen = ceil(qlen / 16)` (u8) or `ceil(qlen / 8)` (i16) stripes of
     `p = 16` or `8` lanes; lanes with `k >= qlen` get profile score `0` (+ shift) for every target base and are **computed like real cells**: they inherit
     `H(i-1, k-1)` along the diagonal with no penalty and take part in the row maximum `imax` (L186, L214), hence in `gmax`/`te` (strict `>` so a padded
     carry-over never moves `te`), in the `b[]` entries that produce `score2`/`csub` (row coalescing is affected because padded lanes keep `imax >= minsc`
     for up to `16*slen - qlen` rows after a real maximum), and in the `Hmax` scan for `qe` (L237-239, index map `i/16 + i%16*slen`). A scalar port must
     either reproduce the striped layout literally or emulate the padded lanes with the same propagation; dropping them changes `csub` (-> `XS` and the
     tandem-repeat MAPQ cap of 7.3) for rescued mates. **[source]**
   * per target row `i`: row maximum `imax`; `b[]` records `(imax, i)` for rows with `imax >= 19`, coalescing consecutive rows into one entry holding
     the max (L215-223); `gmax`/`te` update on **strict** improvement, saving the row (L224-229); u8 stops early if `gmax + shift >= 255` (impossible for
     151-bp reads: `gmax <= 151`).
   * `score = gmax`, `te`; `qe` = the **smallest** query index attaining the maximum in the saved row (L234-239 / L352-357, note the lane->index mapping
     `i/16 + i%16*slen`); `score2`/`te2` = best `b[]` entry with `te2 < te - score || te2 > te + score` (`i = (score + 1 - 1)/1 = score`, L241-249).
   * `KSW_XSTART`: if `score >= 19`, re-run on the reversed prefixes `query[0..qe]`, `target[0..te]` with `KSW_XSTOP | score` (stop at the first row reaching
     `score`); if the reverse score equals `score`, `tb = te - rr.te`, `qb = qe - rr.qe`, else `qb = tb = -1` (L392-400).
4. Accept if `aln.score >= 19 && aln.qb >= 0` (L179): region `b` with `rid`, `qb/qe/rb/re` mapped back to concatenated coordinates (reverse case:
   `qb = l_ms - (aln.qe+1)`, `qe = l_ms - aln.qb`, `rb = 2*l_pac - (rb + aln.te + 1)`, `re = 2*l_pac - (rb + aln.tb)`; L182-185), `score = aln.score`,
   `csub = aln.score2` (`-1` if none), `secondary = -1`, `seedcov = min(re - rb, qe - qb) >> 1`, everything else 0 (`sub`, `sub_n`, `w`, `frac_rep`,
   `truesc`, `seedlen0`, `n_comp`, `hash`). Insert into `ma` after the last element with `score >= b.score` (L191-197). `++n`.
5. After each orientation with `n > 0`: `ma->n = mem_sort_dedup_patch(opt, 0, 0, 0, ma->n, ma->a)` (L201) - dedup only, no patching, re-sorted by
   `(score desc, rb, qb)`.

### 7.2 `mem_mark_primary_se` and `mem_pair` **[source]**

After rescue: `n_pri[i] = mem_mark_primary_se(opt, a[i].n, a[i].a, id<<1|i)` (L303-304); `-5` not set. Then, if both reads have regions:
`o = mem_pair(opt, bns, pac, pes, s, a, id, &subo, &n_sub, z, n_pri)` (L311; `mem_pair` L208-269). **`id` is received as `int`** (truncated
`uint64_t`).

1. Keys (L214-224): for `r = 0,1` and each region `i < n_pri[r]` (all): `x = rid<<32 | (forward coordinate of rb - contig offset)` where the forward
   coordinate is `rb` if `rb < l_pac` else `2*l_pac - 1 - rb` (for a reverse hit this is the forward position of its **last** aligned base);
   `y = score<<32 | i<<2 | is_rev<<1 | r`. Sort by `(x, y)` (`ks_introsort_128`, `utils.c` L45-47: `x` then `y`, all distinct).
2. Scan (L225-253): `y[4] = -1` holds the last index seen per `(is_rev, r)` class (`v.a[i].y & 3`). For each `i`, for `rr = 0,1`:
   `dir = rr<<1 | is_rev(i)` (= strand of the earlier hit, strand of this hit); skip failed orientations; `which = rr<<1 | (r(i) ^ 1)`; walk
   `k = y[which]` downwards over **all** earlier entries with `(y & 3) == which` (other read, strand `rr`): `dist = x[i] - x[k]` (int64, includes
   the rid high bits, so different contigs give huge `dist`); `break` if `dist > high`; `continue` if `dist < low`; otherwise
   `ns = (dist - avg) / std` (double) and
   `q = (int)(score_i + score_k + .721 * log(2. * erfc(fabs(ns) * M_SQRT1_2)) * 1 + .499)`, `q = max(q, 0)`;
   push `u = {x = q<<32 | (hash_64(((uint64_t)k<<32 | i) ^ (uint64_t)(int)(id<<8)) & 0xffffffff), y = k<<32 | i}` (L239-249). `id<<8` is an `int` shift
   (well-defined while `id < 2^23`, i.e. < 8,388,608 pairs per bwa process; sarek chunks are smaller **[inferred]**).
   After both `rr`, `y[v.a[i].y & 3] = i`.
3. If `u` is non-empty (L254-266): sort `u` by `(x, y)`; the best pair is `u[last]`: `z[r] = region index` for both reads (`y<<32>>34`), `ret = q_best`,
   `*sub = u.n > 1 ? q of u[n-2] : 0`, `*n_sub` = number of `u[0..n-2]` with `sub - q_i <= 7`. Else `ret = sub = n_sub = 0`.

Ties in `q` are broken by the 32-bit hash, then by `(k, i)`.

### 7.3 `mem_sam_pe` decision tree (L276-419) **[source]**

```
if n_pri[0] && n_pri[1] && (o = mem_pair(...)) > 0:
    is_multi[i] = exists j in 1..n_pri[i]-1 with a[i].a[j].secondary < 0 && score >= T(30)     // L315-319
    if is_multi[0] || is_multi[1]: goto no_pairing                                            // split reads never take the paired path
    score_un = a[0].a[0].score + a[1].a[0].score - 17
    subo = max(subo, score_un)
    q_pe = raw_mapq(o - subo) = (int)(6.02 * (o - subo) / 1 + .499)                           // L274, L325
    if n_sub > 0: q_pe -= (int)(4.343 * log(n_sub + 1) + .499)
    q_pe = clamp(q_pe, 0, 60)
    q_pe = (int)(q_pe * (1. - .5 * (a[0].a[0].frac_rep + a[1].a[0].frac_rep)) + .499)         // the float sum is computed in float, then promoted
    if o > score_un:                                                                           // paired alignment preferred
        c[i] = &a[i].a[z[i]]
        if c[i]->secondary >= 0: c[i]->sub = a[i].a[c[i]->secondary].score; c[i]->secondary = -2   // promote a secondary to primary
        q_se[i] = mem_approx_mapq_se(opt, c[i])
        q_se[i] = q_se[i] > q_pe ? q_se[i] : min(q_pe, q_se[i] + 40)                          // L339-340
        extra_flag |= 2                                                                        // proper pair
        q_se[i] = min(q_se[i], raw_mapq(c[i]->score - c[i]->csub))                            // tandem-repeat cap; csub may be 0 or -1
    else:
        z[0] = z[1] = 0; q_se[i] = mem_approx_mapq_se(opt, &a[i].a[0])                         // L346-348
    swap secondary_all bookkeeping if z[i] was a secondary (L350-359)
    XA[i] = mem_gen_alt(...)                                                                   // 8.2
    h[i] = mem_reg2aln(a[i].a[z[i]]); h[i].mapq = q_se[i]; h[i].flag |= (0x40 << i) | extra_flag; h[i].XA = XA[i][z[i]]
    (ALT-hit supplementary g[i]: dead, needs is_alt)
    one SAM line per read: mem_aln2sam(opt, bns, str, &s[0], 1, aa[0], 0, &h[1]); likewise read 2 with &h[0]
else:
no_pairing:
    for i in 0,1: which = (a[i].n && a[i].a[0].score >= 30) ? 0 : -1   (ALT fallback dead)
                  h[i] = which >= 0 ? mem_reg2aln(a[i].a[0]) : unmapped record
    if h[0].rid == h[1].rid && h[0].rid >= 0:                                                  // L408-413
        d = mem_infer_dir(l_pac, a[0].a[0].rb, a[1].a[0].rb, &dist)
        if !pes[d].failed && pes[d].low <= dist <= pes[d].high: extra_flag |= 2
    mem_reg2sam(opt, bns, pac, &s[0], &a[0], 0x41 | extra_flag, &h[1])
    mem_reg2sam(opt, bns, pac, &s[1], &a[1], 0x81 | extra_flag, &h[0])
```
`err_fatal` if the two names differ (L386, L416). Note the `no_pairing` branch is taken (a) when a read has no regions, (b) when `mem_pair`
finds no pair within bounds (`o == 0`), or (c) when either read has a second non-secondary hit with score >= 30 (chimeric/split reads). In case
(c) the mate fields of every line of read `i` come from `h[!i]` = the mate's **top** region with its SE mapq. **[source]**

### 7.4 Proper-pair flag, mate fields, TLEN (`bwamem.c` `mem_aln2sam` L851-895) **[source]**

Given `p` (this record, copied) and `m` (the mate record copied from `h[!i]`):
* `flag |= 0x1`; `|= 0x4` if `p->rid < 0`; `|= 0x8` if `m->rid < 0` (L858-860).
* If `p` unmapped and `m` mapped: `p->rid = m->rid, p->pos = m->pos, p->is_rev = m->is_rev, p->n_cigar = 0` (L861-862). If `m` unmapped and `p` mapped:
  `m->rid = p->rid, m->pos = p->pos, m->is_rev = p->is_rev, m->n_cigar = 0` (L863-864).
* `flag |= 0x10` if `p->is_rev` (after the copy above, so an unmapped read placed at a reverse mate gets 0x10); `|= 0x20` if `m->is_rev` (after the copy,
  so a mapped reverse read with an unmapped mate gets 0x20) (L865-866). Resulting flags for one-end-unmapped pairs:
  forward mapped read 1: `1+8+64 = 73`, its mate `1+4+128 = 133`; reverse mapped read 1: `1+8+16+32+64 = 121`, its mate `1+4+16+32+128 = 181`;
  symmetric 137/69 and 185/117 for read 2. Flags 89/165 and 153/101 (reverse read, mate unmapped, without 0x20/0x10) cannot occur. Both unmapped: 77/141. **[source]**
* `0x2` comes only through `extra_flag` (7.3); `0x40/0x80` from `0x41/0x81` or `0x40<<i`; `0x100` for secondary (never printed, 8.2); `0x800` per 8.2.
* RNEXT/PNEXT (L882-893): if `m->rid >= 0`: `=` when `p->rid == m->rid` else the mate's contig name; `m->pos + 1`. TLEN if same contig:
  `p0 = p->pos + (p->is_rev ? get_rlen(p) - 1 : 0)`, `p1` likewise for `m` (`get_rlen` = sum of M and D lengths, L827-836); if either `n_cigar == 0`
  print `0`; else print `-(p0 - p1 + sgn(p0 - p1))`, i.e. `p1 - p0 + 1` if `p0 < p1`, `p1 - p0 - 1` if `p0 > p1`, `0` if equal. This is the distance
  between the two reads' 5'-end coordinates (leftmost base for forward, rightmost base for reverse), **not** the SAM-spec leftmost-to-rightmost
  span. With `m` unmapped-but-placed, `m->n_cigar == 0` -> TLEN `0`. Otherwise (`m->rid < 0`, both unmapped): `*\t0\t0`.

---

## 8. Alignment finalisation and the SAM line

### 8.1 `mem_reg2aln` (L1119-1189) and `bwa_gen_cigar2` (`bwa.c` L148-234) **[source]**

* `ar == 0` -> unmapped `mem_aln_t`: all zero except `rid = -1, pos = -1, flag = 0x4` (so `mapq = 0, score = 0, sub = 0, n_cigar = 0, XA = 0`).
* `mapq = secondary < 0 ? mem_approx_mapq_se(ar) : 0`; `flag |= 0x100` if `secondary >= 0` (L1136-1137).
* Band for the global re-alignment (L1138-1142): `infer_bw(l1, l2, truesc, a, q, r)` (L818-825) = `0` if `l1 == l2 && l1 - truesc < 12`
  (`(q + r - a) << 1`), else `max((int)((double)(min(l1,l2) - truesc - q) / r + 2.), |l1 - l2|)`; `w2 = max(infer_bw(... o_del, e_del), infer_bw(... o_ins, e_ins))`;
  if `w2 > 100` then `w2 = min(w2, ar->w)`.
* Loop (L1143-1152): `w2 = min(w2, 400)`; `cigar = bwa_gen_cigar2(mat, 6, 1, 6, 1, w2, l_pac, pac, qe - qb, query + qb, rb, re, &score, &n_cigar, &NM)`;
  break if `score == last_sc` or `w2 == 400`; `last_sc = score; w2 <<= 1`; repeat while `++i < 3 && score < truesc - 1`. (`w2 == 0` stays 0.)
* `bwa_gen_cigar2`: reject if `rb >= re` or bridging `l_pac` (L159). `rseq = bns_get_seq(l_pac, pac, rb, re)` (strand-aware). If `rb >= l_pac` **reverse
  both the query slice and `rseq`** (L162-167): the alignment is then computed in forward-strand orientation, so the CIGAR comes out left-to-right on the
  forward strand without further reversal; the query is reversed back at the end (L227-229). If `l_query == re - rb && w_ == 0` (reachable when `infer_bw`
  returned 0 for both gap types): no DP - `cigar = [l_query M]`, `score = sum of mat[rseq[i]*5 + query[i]]` (L168-176). Otherwise
  `max_ins = max_del = (int)((double)(((l_query+1)>>1) - 6) / 1 + 1.)`, `max_gap = max(max_ins, 1)`, `w = (max_gap + |rlen - l_query| + 1) >> 1`,
  `w = min(w, w_)`, `w = max(w, |rlen - l_query| + 3)` and `ksw_global2(l_query, query, rlen, rseq, 5, mat, 6, 1, 6, 1, w, &n_cigar, &cigar)` (L178-194).
* `ksw_global2` (`ksw.c` L540-642): banded Needleman-Wunsch with affine gaps and a direction matrix; port verbatim. Tie rules: `d = (m >= e) ? 0 : 1`,
  then `d = (h >= f) ? d : 2` (match preferred over deletion over insertion on equal H); gap-extension bits set only when the extension is **strictly**
  better than opening (`e > t`, `f > t`). First row `eh[j].h = -(6 + j)` for `j <= w`, `-inf` beyond; first column `-(6 + (i+1))` when `beg == 0`. Score
  `= eh[qlen].h`. Backtrack from `(tlen-1, k = min(tlen-1+w+1, qlen)-1)` with a 3-state machine reading bits `(which<<1)` of `z` (L627-633); a remaining
  `i >= 0` becomes `(i+1) D`, remaining `k >= 0` becomes `(k+1) I`; CIGAR reversed. Ops are `0=M 1=I 2=D`.
* `NM`/`MD` (L196-226) over the (forward-oriented) CIGAR: for `M` runs compare `query[x+i] != rseq[y+i]` (a read `N`, code 4, always mismatches
  because `.pac` only holds codes 0-3, see 2.2; it scores `-1` in the DP and counts 1 in `NM`); emit `<u><refbase>` per mismatch
  using `int2base = "ACGTN"` for forward hits and `"TGCAN"` for reverse hits (complement back to forward strand); `D` runs emit `<u>^<bases>` and count
  toward NM **only if the op is neither the first nor the last CIGAR op**; `I` runs add their length to NM and do not break the MD match count.
  `MD` ends with the final `u`. The MD string is stored right after the CIGAR array (`str.s = (char*)cigar`, L198) and read back via
  `(char*)(p->cigar + p->n_cigar)` (L932).
* Position and clipping (L1153-1184): `pos = bns_depos(rb if forward else re - 1)`; squeeze **one** leading `D` (`pos += len`, drop op) **else** one
  trailing `D` (drop op, `MD` unchanged since it never included it) (L1157-1166; if both exist only the leading one is removed). Clips:
  `clip5 = is_rev ? l_query - qe : qb`, `clip3 = is_rev ? qb : l_query - qe`; prepend `clip5 S` and append `clip3 S` (op code 3) when non-zero
  (L1167-1181). `rid = bns_pos2rid(pos)`, `pos -= anns[rid].offset` (0-based). `score = ar->score` (local max, 5.3), `sub = max(ar->sub, ar->csub)`.

### 8.2 Which hits are printed: `mem_reg2sam` (L1033-1079), supplementary/secondary, `XA` (`bwamem_extra.c` L116-172), `SA` **[source]**

`mem_reg2sam` (used for both reads in the `no_pairing` path and for all single-end output):
* `XA = mem_gen_alt(...)` (no `-a`). For `k = 0..n-1` in array order: skip if `score < 30`; skip if `secondary >= 0` (secondary alignments are **never
  printed** without `-a`; the `drop_ratio` test on L1050 is unreachable); otherwise `q = mem_reg2aln(p)`, `q.XA = XA[k]`, `q.flag |= extra_flag`;
  the first printed hit (`l == 0`) is the primary; every later one gets `flag |= 0x800` (no `-M`) (L1057-1058); and, because `-q`/`-5` are absent,
  `q.mapq = min(q.mapq, mapq of the first printed hit)` (L1059-1060). If nothing qualifies, one unmapped record with `extra_flag` (L1063-1067).
  Then `mem_aln2sam` for each `k` with `n = aa.n`, `which = k`, mate `m`.
* Primary per read = `a[0]` in the `(score desc, is_alt, hash)` order of 6.1 (or `a[z[i]]` in the paired path). The `0x800` rule: in `aa` order, every
  non-secondary hit after the first.
* `SA` (`mem_aln2sam` L939-958): for a record without `0x100`, if any other `list` entry lacks `0x100`, append `SA:Z:` followed by, for each other
  such entry in `list` order, `<rname>,<pos+1>,<+|->,<CIGAR>,<mapq>,<NM>;`. The CIGAR here prints the stored ops with `"MIDSH"` (with `-Y` all clips are
  `S`); `mapq` is the value after the supplementary cap. In the paired path `n == 1`, so `SA` only appears in the `no_pairing` path (split reads).
* `XA` (`mem_gen_alt`, L124-172): for each region `i` with `k = secondary_all >= 0` and `a[i].score >= a[k].score * 0.8f` (`int * float` -> float,
  compared as float; `0.8f = 0.800000011920929`, so e.g. a secondary with score 4 against a primary with score 5 does **not** qualify), count `cnt[k]`.
  For each such `i` in array order, unless `cnt[k] > 5` (`max_XA_hits`; the 200 limit needs an alt hit), append
  `<rname>,<+|-><pos+1>,<CIGAR with "MIDSHN">,<NM>;` built from `mem_reg2aln(a[i])` (mapq ignored). The result for primary `k` is printed as
  `XA:Z:...` on `k`'s line (L962-965). A primary with more than 5 qualifying secondaries gets **no** `XA`.
* `pa:f` (L959-960) needs `alt_sc > 0`: never. `XR` needs `-V`: never.

### 8.3 SAM text, field by field (`mem_aln2sam` L851-976) **[source]**

One line per printed hit, appended to the read's `kstring`, terminated by `\n` (L975). Integers are written with `kputw`/`kputl`
(`kstring.h` L62-77, L94-109): plain decimal, `-` sign, no padding.

| Field | Rule |
|---|---|
| QNAME | `s->name` (1.1) |
| FLAG | `(flag & 0xffff) \| (flag & 0x10000 ? 0x100 : 0)` (L872); bits per 7.4 / 8.2 |
| RNAME, POS, MAPQ, CIGAR | if `rid >= 0`: contig name, `pos + 1`, `mapq`, CIGAR via `add_cigar` (L838-849: ops `"MIDSH"[c]`; with `-Y` no S->H rewrite; `n_cigar == 0` prints `*`, which is the case for a mate-placed unmapped read); else `*\t0\t0\t*` (L873-878) |
| RNEXT, PNEXT, TLEN | per 7.4 |
| SEQ, QUAL | secondary (`0x100`): `*\t*`. Forward: `"ACGTN"[code]` for all bases and the raw quality string. Reverse: `"TGCAN"[code]` iterating from the last base to the first, and the quality string reversed (L896-927). The hard-clip trimming (`qb/qe` adjustments) requires `!MEM_F_SOFTCLIP`, so with `-Y` the full read is always printed. If `s->qual == 0` print `*` (never for FASTQ) |
| `NM:i:` `MD:Z:` | only if `n_cigar > 0` (L930-933) |
| `MC:Z:` | if `m && m->n_cigar` (L934): the mate's CIGAR via `add_cigar` (omitted when the mate is unmapped) |
| `MQ:i:` | if `m` (every paired record, L935): `m->mapq` (0 for an unmapped mate) |
| `AS:i:` | if `score >= 0` (L936): always, including unmapped records (`AS:i:0`) |
| `XS:i:` | if `sub >= 0` (L937): always except secondary hits (`sub = -1`, L1056); unmapped records print `XS:i:0` |
| `RG:Z:` | if `bwa_rg_id[0]` (L938): `RG:Z:HG002_WES.L1` |
| `SA:Z:` | per 8.2 (only non-secondary records with other non-secondary hits) |
| `pa:f:` | never |
| `XA:Z:` | if `p->XA` (L962-965) |
| comment | never (`-C` absent) |
| `XR:Z:` | never (`-V` absent) |

So the tag order is exactly `NM MD MC MQ AS XS RG SA XA`, each present under the conditions above - matching the golden
`NM MD MC MQ AS XS RG SA`. A fully unmapped pair member prints `<name>\t77|141\t*\t0\t0\t*\t*\t0\t0\t<SEQ>\t<QUAL>\tMQ:i:0\tAS:i:0\tXS:i:0\tRG:Z:HG002_WES.L1`;
an unmapped read with a mapped mate prints `<name>\t<flag>\t<mate rname>\t<mate pos+1>\t0\t*\t=\t<mate pos+1>\t0\t<SEQ>\t<QUAL>\tMC:Z:<mate cigar>\tMQ:i:<mate mapq>\tAS:i:0\tXS:i:0\tRG:Z:...`
where SEQ/QUAL are reverse-complemented if the mate is reverse (because `is_rev` was copied). **[source]**

The mate record `m` for read `i` is always `h[!i]` (the mate's primary line, with the paired-path mapq `q_se` or the SE mapq in `no_pairing`),
so `MC`/`MQ` are identical on every line of a split read. **[source]**

---

## 9. `samtools sort` (samtools 1.22.1 `bam_sort.c`, htslib 1.22.1)

### 9.1 Reading bwa's SAM text (`htslib/sam.c` `sam_parse1` L2877-3055, `aux_parse` L2739-2875) **[source]**

* QNAME copied, NUL-terminated and padded to a 4-byte boundary (`l_extranul`, L2939; padding is dropped again on write, `bam_write1` L856-858).
* FLAG via `parse_sam_flag` (L2713-2729): decimal. RNAME `*` -> `tid = -1`; POS `hts_str2uint - 1`; `tid < 0` forces `BAM_FUNMAP` (L2965, no-op for bwa's
  flags). MAPQ 8-bit. CIGAR `*` keeps `n_cigar = 0` and forces `BAM_FUNMAP` (already set). `bin = hts_reg2bin(pos, pos + cigreflen, 14, 5)` (L2989) with
  `cigreflen` = reference length of the CIGAR, or `1` for unmapped / zero-length (L2979-2983); for `tid = -1, pos = -1` this gives 4680. RNEXT `=` -> `mtid = tid`,
  `*` -> `-1`; PNEXT `- 1`; TLEN `hts_str2int`. SEQ packed with `seq_nt16_table` (`hts.c` L239-255: `A=1 C=2 G=4 T=8 N=15`); QUAL = ASCII `- 33` (L3038).
* Aux integers (`aux_parse` L2807-2834): bwa writes every integer tag as `i`. htslib stores it as the **smallest type by value**, choosing unsigned for
  non-negative values: value with leading `-`: `>= -128 -> c`, `>= -32768 -> s`, else `i`; otherwise `<= 255 -> C`, `<= 65535 -> S`, else `I`. Hence the
  golden's `NM:C:0`, `MQ:C:60`, `AS:C:69`, `XS:C:19`, and e.g. `AS:S:300`. `Z` tags are copied verbatim with a NUL. Aux order is preserved as written
  by bwa (`NM MD MC MQ AS XS RG SA XA`). **[source]**

### 9.2 Options and buffering (`bam_sort` L3733-3880, `bam_sort_core_ext` L3356-3660) **[source]**

* Defaults: coordinate order (`sam_order = Coordinate`), `max_mem = 768 MiB` per thread (L135, L3735), `modeout = "wb"` (default BGZF level);
  `--threads 20` -> `ga.nthreads = 20` -> `n_threads = 20`, `max_mem = 768 MiB * 20` (L3401). Input `-` is the SAM stream on stdin; `-o` names the output;
  `--no-PG` not given (`no_pg = 0`); `arg_list = stringify_argv(argc+1, argv-1)` (L3838; `htslib/textutils.c` L419-447: all argv including the program
  name joined by single spaces, TABs replaced by spaces) = `samtools sort --threads 20 -o HG002_WES.0001.bam -`.
* Records are appended to `buf[]` in input order. When `max_mem` fills (L3466-3472), the block is sorted and written to a temporary BAM (L3521-3560) and
  the buffer restarts; the final block stays in memory. Temporary files are merged with the same comparator and the same tie-breaking (9.3), and the
  file index precedes the in-memory blocks, so **the result is independent of whether spilling occurred**. **[source]**

### 9.3 Comparator and stability **[source]**

* In-memory block sort: `sort_blocks` (L3131-3175) splits the buffer into `n_threads` **contiguous** slices in input order (one slice if `k < 64 * n_threads`),
  each sorted by `worker` -> `ks_radixsort` (L3109-3128, L2381-2460): an LSD byte-radix sort (counting pass placing elements from `i = n-1` down to 0,
  i.e. **stable**) on the key `(tid' , (pos + 1) << 1 | is_rev)` where `tid' = tid == -1 ? n_ref : tid` (unmapped last), bytes little-endian. So each slice
  ends up sorted by `(tid', pos, strand)` with ties in **input order**.
* Merge (`bam_merge_simple` L1838-1960, `heap_lt` L226-272, `heap_add_read` L1790-1836): a `ks_heap` over one head per slice/file, ordered by `tid`
  (as `uint32_t`, so `-1` sorts last), then `pos + 1`, then `rev`, then `i` (slice/file index, lower first), then `idx` (a global counter incremented each time
  a record is pulled into the heap - within one slice this is monotone). `ks_heapmake`/`ks_heapadjust` with a strict comparator that never returns
  "equal" for two live records (`(i, idx)` is unique), so the heap order is a total order and the pop sequence is deterministic. Result: a **stable** sort
  of the entire input by `(tid with -1 last, POS, strand)`. If `n_files == 0 && num_in_mem < 2`, `write_buffer` (L2326-2372) emits the single slice
  directly - same order.
* `bam1_cmp_core` (L2012-2043) / `bam1_lt` (L2291-2309) implement the same `(tid as uint64, pos+1, is_rev)` key; they are used by `ks_mergesort` for the
  non-coordinate orders and are not reached here (L3114-3118).
* Equivalent specification for the port: output = bwa's record stream, stably sorted by `(tid == -1 ? +inf : tid, POS, FLAG & 0x10)`. Mate-placed
  unmapped reads sort at their placed position (and with their copied strand bit); fully unmapped reads keep bwa's relative order at the end.
  **[source]** Note this differs from htsjdk's coordinate comparator (markdup.md 1.5), which is why the downstream merge in markdup.md must be replayed as a
  k-way merge. **[inferred]**

### 9.4 Output header **[source]**

* `set_sort_order` (L3249-3330) after reading everything: `sam_hdr_update_sort(h, "SO", "coordinate")` (macro L3238-3243 -> `sam_hdr_update_line(h, "HD", NULL, NULL,
  "SO", "coordinate", NULL)`, `header.c` L1548-1600 -> `sam_hrecs_vupdate` L2578-2616: an existing `SO` tag is **replaced in place**, keeping tag order);
  then `sam_hdr_remove_tag_hd(h, "GO")` and `(h, "SS")` (`sam_hdr_remove_tag_id` L1984-2008). bwa's `@HD\tVN:1.5\tSO:unsorted\tGO:query` therefore becomes
  `@HD\tVN:1.5\tSO:coordinate`.
* `sam_hdr_add_pg(hout, "samtools", "VN", "1.22.1", "CL", arg_list, NULL)` (L1918-1922 / L2337-2340; `header.c` L2251-2340): `sam_hdr_link_pg` finds the
  chain ends (the `bwa` record has no `PP`, so it is an end); for each end a new record is added with explicit tags `ID`, `PN`, `PP` first and then the
  caller's `VN`, `CL` (`sam_hrecs_vadd` L549-660: explicit varargs before the `va_list`; empty-string values skipped). `ID` = `sam_hdr_pg_id(h, "samtools")`
  (L2200-2230) = `samtools` because no such ID exists. The new line is placed in the global list **after the last existing `@PG`** (L586, `sam_hrecs_global_list_add(..., h_type->prev)`). Result:
  `@PG\tID:samtools\tPN:samtools\tPP:bwa\tVN:1.22.1\tCL:samtools sort --threads 20 -o HG002_WES.0001.bam -`.
* The BAM header text is regenerated from the parsed records (`bam_hdr_write` L343-400 -> `sam_hrecs_rebuild_text` L2015 -> `sam_hrecs_rebuild_lines`
  L755-770 / `build_header_line` L740-753): `@XX` + `\t` + each tag string (stored verbatim) + `\n`, in the original line order. Since bwa's lines are already
  canonical this reproduces them byte-for-byte; `n_ref` and the name/length table follow from the `@SQ` lines. Final header line order:
  `@HD`, 3366 x `@SQ`, `@RG`, `@PG bwa`, `@PG samtools`. **[source]**

### 9.5 Record encoding (`bam_write1` L856-915) **[source]**

Standard BAM: `block_size`, `refID`, `pos`, `l_read_name`(incl. NUL) / `mapq` / `bin`, `n_cigar_op` / `flag`, `l_seq`, `next_refID`, `next_pos`, `tlen`, then
name+NUL (no padding), CIGAR, packed SEQ, QUAL, aux bytes as parsed in 9.1. `bin` is the value computed by `sam_parse1` (9.1).

---

## 10. Floating-point inventory (every FP operation that can change output)

All of bwa's arithmetic is IEEE double unless the operand types below say `float`. `(int)` casts truncate toward zero. Each line gives the
exact expression; keep operation order.

**Correction (verified 2026-10-08): the golden binary is aarch64 and uses fused multiply-add.** The golden `bwa` (bioconda `bwa-0.7.19-h0cbc5ad_1`,
extracted from the Wave image) was built by conda-forge GCC 13.3.0 with `-O2 -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt` (from
`DW_AT_producer`) and no `-std`, so GCC's default `-ffp-contract=fast` fuses `a*b + c` into `fmadd`/`fmsub` (single rounding). rustc never fuses:
each site below must use `f64::mul_add` with the operands the compiler chose (read the disassembly at that address, not just the C). Sites, from
`objdump -d` + `addr2line` on the golden (excluding timing/hash-table code):

| Function | Source line(s) | Instructions |
|---|---|---|
| `mem_flt_chained_seeds` | bwamem.c:627 | 1 fmadd |
| `mem_patch_reg` | bwamem.c:455, 456 | 2 fmadd |
| `mem_approx_mapq_se` | bwamem.c:996, 998, 999, 1001, 1004 | 5 fmadd |
| `mem_pestat` | bwamem_pair.c:103, 104, 105, 106 (fmsub), 108, 117, 120 (fmsub), 121 | 8 |
| `mem_pair` | bwamem_pair.c:244 | 1 fmadd |
| `mem_sam_pe` | bwamem_pair.c:326, 329 (fmadd + fmsub) | 3 |

A source build on the Spark (Ubuntu GCC 13.3.0, same flags) has the identical 58-site FMA inventory and reproduces the golden chunk BAMs; it is
the stage-by-stage test oracle (`reference/oracle/`). Rebuilt with `-ffp-contract=off` (zero FMAs) it **also** reproduces all 12 tier0
chunks, so no tier0 record depends on contraction; the sites stay mirrored with `mul_add`. The
port now matches the golden directly on all of tier1 (81.58M records), which is the binding check; an FMA/no-FMA oracle comparison on
tier1 has not been run.

| Where | Expression (types) | libm |
|---|---|---|
| `mem_collect_intv` L144 | `split_len = (int)(19 * 1.5f + .499)` -> `(float)28.5 + 0.499` in double -> 28 | - |
| `mem_chain` L336 | `frac_rep = (float)l_rep / (float)len` (float) | - |
| `mem_chain_flt` L381, L384 | `e_min - b_max >= (float)min_l * 0.5f`; `(float)a[i].w < (float)a[j].w * 0.5f` | - |
| `mem_flt_chained_seeds` L626-628 | `min_l = 5.5f * log((double)l_query)`; `min_l > (double)(0.05f * (float)l_query)` | `log` |
| `cal_max_gap` L649-650 | `(int)((double)(qlen - 6) / 1 + 1.)` (integers exact) | - |
| `mem_chain2aln` L702, L722 | `s->len - p->seedlen0 > .1 * l_query`; `t->len < s->len * .95` | - |
| `ksw_extend2` L438, L441 | `(int)((double)(qlen*1 + 5 - 6) / 1 + 1.)` | - |
| `mem_sort_dedup_patch` L481 | `(float)or > 0.95f * (float)mr` (and `oq`/`mq`) - int64 -> float conversions | - |
| `mem_patch_reg` L442-458 | `r` (double), `r >= (double)0.05f` / `(double)0.1f`; `q_s`, `r_s` with `+ .499`; `(double)score / max(q_s, r_s) < (double)0.9f` | - |
| `bwa_gen_cigar2` L180-181 | `(int)((double)(((l+1)>>1) - 6) / 1 + 1.)` | - |
| `infer_bw` L822 | `(int)((double)(min(l1,l2) - score - q) / r + 2.)` | - |
| `mem_mark_primary_se_core` L534, `cal_sub` L66 | `e_min - b_max >= (float)min_l * 0.5f` | - |
| `mem_approx_mapq_se` L989-1004 | see 6.2: `identity`, `3 / log(l)`, `6.02 * (score - sub) / 1 * tmp * tmp + .499`, `4.343 * log(sub_n + 1) + .499`, `mapq * (1. - (double)frac_rep) + .499` | `log` |
| `mem_pestat` L85-86, L103-123 | `cal_sub > 0.8 * score`; percentile indices `(int)(.25 * n + .499)` etc.; `avg` as a double sum of `uint64` in sorted order `/ x`; `std = sqrt(sum((q - avg)^2) / x)`; bounds with `2.0`, `3.0`, `4.0` and `+ .499` | `sqrt` |
| `mem_pair` L243-244 | `ns = (double)(dist - avg) / std`; `q = (int)(score_i + score_k + .721 * log(2. * erfc(fabs(ns) * M_SQRT1_2)) * 1 + .499)` (`dist - avg`: int64 - double) | `log`, `erfc` |
| `mem_sam_pe` L325-329, L343-344 | `raw_mapq = (int)(6.02 * diff / 1 + .499)`; `(int)(4.343 * log(n_sub + 1) + .499)`; `(int)(q_pe * (1. - .5 * (double)(frac_rep0 + frac_rep1)) + .499)` where `frac_rep0 + frac_rep1` is a **float** addition | `log` |
| `mem_gen_alt` / `get_pri_idx` L119 | `(float)a[i].score >= (float)a[k].score * 0.8f` | - |
| `mem_reg2sam` L1050 | unreachable (needs `-a`) | - |

Port guidance **[inferred]**: `f32` for every `float` expression above (Rust `as f32`/`f32` ops match C with SSE); `f64` elsewhere; `sqrt` is
correctly rounded everywhere. `log` and `erfc` must match the libm of the reference run (glibc in the sarek container): Rust's `f64::ln` calls the
platform libm on Linux (same results when run on glibc), but Rust has no `erfc` - bind the system `erfc` through `libc`/`extern "C"` on Linux rather
than using the `libm` crate (musl/FreeBSD msun), whose last-ulp differences can flip the `(int)(... + .499)` truncations in `mem_pair`. Results are
only sensitive where the pre-truncation value is within ~1e-12 of an integer boundary, so mismatches would be rare but real; test against the
golden `MAPQ` distribution of proper pairs.

---

## 11. Other things that affect exact output

* **Sort stability in bwa** (12.2): `ks_introsort` (`ksort.h` L176-225 plus `ks_combsort` L155-175 and `__ks_insertsort` L146-153) is unstable and is used
  with non-total comparators in `mem_chain_flt` (`w` only) and `mem_sort_dedup_patch` (`re` only). Reproduce the exact algorithm: `n == 2` special
  case; depth `d` with `1ul<<d < n`; median-of-three pivot `k = i + ((j-i)>>1) + 1` with the comparisons on L197-199; the partition loop on L201-206;
  push the larger side, recurse into the smaller, with partitions of `<= 16` elements left for one final insertion sort; falling back to combsort when
  depth is exhausted. Both `mem_ars_hash` and `mem_intv`/`mem_ars`/`ks_introsort_64/128` are effectively total here.
* `hash_64` (`utils.h` L98-109) is used for (a) region tie-breaks (6.1) and (b) pair tie-breaks (7.2); both depend on the read's global pair index
  within the bwa process (`n_processed`), so the port must count reads from the start of the chunk exactly as `bseq_read` does.
* `mem_reorder_primary5` is only called with `-5` (dead).
* Memory layout concerns: `mem_alnreg_t` is `memset` to 0 in `mem_chain2aln` (L735) and `mem_matesw` (L178); `hash` of rescued regions is 0 until
  `mem_mark_primary_se` sets it; `w = 0` for rescued regions affects `infer_bw` (`w2 = min(w2, ar->w)` only when `w2 > 100`, so a rescued region with
  `w2 > 100` gets `w2 = 0` -> no-DP or minimal-band CIGAR path) - port faithfully.
* Reads with `l_seq < 19` are unmapped (`mem_chain` L286); reads consisting only of `N` likewise.
* `assert`s (`bwa_gen_cigar2`, `mem_reg2aln` L1183, `mem_chain2aln` L686, L776) do not fire on valid data; a port should mirror them as debug checks.
* stderr output (`[M::mem_pestat] ...`, `[M::process] read ...`) is not part of the deliverable but is a cheap cross-check of batch boundaries and insert-size
  statistics if the pipeline's `.command.log` files are available.

---

## 12. Summary of the algorithm to reproduce (Rust implementation notes)

1. Load `.64.{bwt,sa,ann,amb,pac}`; implement `Occ`, `bwt_sa`, `bwt_extend` bit-exactly (2.1); treat `.pac` as the reference (2.2).
2. Parse the two FASTQs with kseq semantics (1.1); form batches by cumulative base count (1.2); keep `n_processed`.
3. Per read: seeding (3) -> chaining with a verbatim kbtree (4.1) -> `mem_chain_flt` with a verbatim introsort (4.2) -> extension (5) -> dedup/patch (5.5).
4. Per batch: `mem_pestat` (1.3).
5. Per pair: mate rescue with a scalar re-implementation of `ksw_u8`/`ksw_align2` (7.1) -> `mem_mark_primary_se` with `hash_64` (6.1) -> `mem_pair` (7.2)
   -> the decision tree (7.3) -> `mem_reg2aln` with `ksw_global2` (8.1) -> `XA`/`SA` (8.2) -> SAM lines (8.3).
6. Convert each SAM line to a BAM record with htslib's typing rules (9.1); stable-sort by `(tid with -1 last, POS, strand)` (9.3); emit the header of 9.4.

Gating order suggested by the data flow: (i) seeds/chains/regions for single reads (compare `AS`, CIGAR, `NM`, `MD`, POS of reads whose pair took
the paired path and have no `XA`/`SA`); (ii) `XS`/`XA` (needs 6.1 hashing); (iii) `mem_pestat` per batch (compare proper-pair flags and TLEN);
(iv) rescue + pairing + MAPQ; (v) split reads (`SA`, 0x800, mapq cap); (vi) sort + header.

---

## 13. Checks against the tier0 golden

Files: chunk FASTQs under `/home/<user>/omics/runs/tier0-wes-20261008-r1/work/.../0001.HG002_WES-L1_R{1,2}.fastp.fastq.gz` and the chunk BAMs
`results/preprocessing/mapped/HG002_WES/HG002_WES.NNNN.bam`. Run `samtools view -h --no-PG` (tag order preserved) unless stated.

Header:

1. Line 1 is exactly `@HD\tVN:1.5\tSO:coordinate` (no `GO`, no `SS`).
2. `@SQ` lines are `SN:` + `LN:` only (no `AH:*`, no `M5`), 3366 of them, in `.ann` order; `@RG` is the `-R` string with TABs; then exactly two `@PG`
   lines: `ID:bwa PN:bwa VN:0.7.19-r1273 CL:...` and `ID:samtools PN:samtools PP:bwa VN:1.22.1 CL:samtools sort ...` with tags in that order.
3. The bwa `CL` contains the characters `\t` (backslash, t) inside the `-R` argument, not TABs.

Batching / determinism:

4. Summing `len(R1) + len(R2)` over pairs in file order, the first batch ends at the first pair where the running sum `>= 100000000`; every pair's
   proper-pair/TLEN behaviour must be explainable with that batch's `mem_pestat` values. If `.command.log` exists, `[M::process] read N sequences (M bp)`
   lines give `N` (reads) and `M` (bases) per batch: `M >= 1e8` for all but the last, `N` even, and the sum of `N` equals twice the pair count.
4b. If the last batch of a chunk (pairs after the final 1e8-base boundary) contains fewer than 10 unique FR pairs, none of its records carries `0x2`
    and none has an `XA`/MAPQ pattern requiring rescue; conversely every batch with `>= 10` such pairs has proper pairs. (Locate the batch's reads by
    name order in the FASTQ.)
5. Record names: QNAME has no `/1`/`/2` suffix and no comment; `samtools view | cut -f1 | sort | uniq -c` shows every name exactly 2 + (number of
   supplementary lines) times.

Records:

6. For every record the aux tags appear in the order `NM MD MC MQ AS XS RG SA XA` (subset), and the raw BAM type bytes are `C` for 0..255 and `S` for
   256..65535 on `NM/MQ/AS/XS` (never `c`/`s`/`i`).
7. Every record has `MQ` and `AS` and `RG`; `XS` is present on every record (no secondary records exist: `samtools view -c -f 0x100` = 0);
   unmapped records (`-f 4`) have `AS:C:0 XS:C:0`, no `NM`/`MD`, and `MC` iff the mate is mapped.
8. `samtools view -c -f 0x800` > 0 and each such record carries an `SA` tag and the proper-pair bit identical to its primary line; all CIGARs use `S`,
   never `H` (`-Y`), including inside `SA`.
9. For any record with `SA`, every other line of the same read is listed once in `SA` in the order (primary first, then supplementaries in `aa` order),
   and the supplementary's MAPQ `<=` the primary's MAPQ.
10. Flag rules for one-end-unmapped pairs: the only flag pairs are (73,133), (121,181), (137,69), (185,117); i.e. an unmapped read whose mate
    is reverse has `0x10` **and** `0x20` set (181/117), and a mapped reverse read with an unmapped mate has `0x20` set (121/185). Flags 89, 165,
    153, 101 never occur. Both-unmapped pairs are exactly 77/141. No record has `0x100`.
11. For an unmapped read with a mapped mate: RNAME/POS equal the mate's, MAPQ 0, CIGAR `*`, RNEXT `=`, PNEXT = POS, TLEN 0, `MC` = mate CIGAR, `MQ` = mate
    MAPQ, and SEQ/QUAL are the reverse complement/reverse of the FASTQ record iff the mate is reverse.
12. TLEN check on proper pairs (same contig, both mapped): with `p0` = POS (forward) or POS + reflen(CIGAR) - 1 (reverse) and `p1` the mate's, TLEN =
    `p1 - p0 + 1` if `p0 < p1`, `p1 - p0 - 1` if `p0 > p1`, else 0; the two mates' TLENs are negatives of each other except when `p0 == p1` (both 0).
13. `MD`/`NM` consistency: NM = mismatches in MD + inserted bases + deleted bases; MD never starts with `^` and never ends with a `^` run.
14. SEQ contains only `ACGTN` uppercase even if the FASTQ had lowercase/IUPAC; reverse-strand records are the exact reverse complement of the FASTQ
    sequence (compare against the FASTQ by name).
15. XA tags: each entry `chr,[+-]pos,CIGAR,NM;`; count per record `<= 5`; a record never has `XA` with more than 5 entries (records with 6+
    qualifying secondaries have none).
16. MAPQ of reads taking the paired path is `<= 60` and equals `min(q_se', raw_mapq(score - csub))` as computable from 6.2/7.3 given the regions; at
    minimum verify `MAPQ <= 60` everywhere and `MAPQ == 0` whenever `XS >= AS` for unpaired-path reads (`sub >= score` -> 0 in 6.2).
17. Sort order: for consecutive records `(tid, POS, 0x10)` is non-decreasing with `tid = -1` last; among equal keys the order equals bwa's emission order
    (pair index, then R1 lines before R2 lines, then primary before supplementary). Checkable on the unplaced tail: it must be in increasing
    original-FASTQ order with read 1 before read 2.
18. `bin` field: unplaced unmapped records have `bin == 4680`; mate-placed unmapped records `bin == reg2bin(pos, pos+1)`.
19. Byte-identical chunk BAMs across two runs (already observed) - the spec predicts this; any `-t`-dependent difference would contradict 0.4.

Input assumptions to confirm on the FASTQs:

20. Every header line is `@<name> <comment>` with `<name>` free of whitespace and not ending in `/[0-9]`; R1/R2 names match pairwise.
21. Qualities are ASCII-33 Phred restricted to `#`, `,`, `:`, `F` (irrelevant to bwa, but confirms QUAL passthrough).

### 13.1 Verification record (2026-10-08)

`reference/check_bwa_spec.py` (parses BAM directly, so raw aux type bytes and `bin` are visible) implements checks 1-3, 5-18, 20-21.
On all 12 tier0 r1 chunk BAMs: **28/28 pass** per chunk. Notes:

- Check 10's flag pairs are unordered; observed (R1, R2) pairs include (69, 137) and (117, 185), i.e. the listed pairs with the mates swapped.
  No both-unmapped pairs occur in tier0.
- Check 4/4b are **not exercised by tier0**: each tier0 chunk is a single bwa batch (`[M::process] read ~81,600 sequences (~12.3 Mbp)`,
  below `-K 1e8`). Tier1 chunks (~1 Gbp each) span ~10 batches and are where batch-boundary behaviour must be gated.
- Check 16 is verified only as `MAPQ <= 60`; the full MAPQ formula is gated by record comparison once alignment runs.

### 13.2 Stage gate coverage (tier0, 12 chunks, oracle dumps)

| Stage | Lines | Status |
|---|---|---|
| `M` seeds (`mem_collect_intv`) | 979k reads | IDENTICAL (chunk 0001; implied on all by `c`) |
| `c` chains (`mem_chain`, kbtree order) | 4,083,621 | IDENTICAL, all chunks |
| `C` filtered chains (`mem_chain_flt`) | 1,874,824 | IDENTICAL, all chunks |
| `R` regions (`mem_chain2aln` + `mem_sort_dedup_patch`) | 1,875,551 | IDENTICAL, all chunks |

Not exercised by tier0: `mem_patch_reg` merges (no region has `n_comp > 1`), band retries (every `w` is 100),
`mem_flt_chained_seeds` (reads < ~700 bp; nomics refuses longer reads). Re-check on tier1.
