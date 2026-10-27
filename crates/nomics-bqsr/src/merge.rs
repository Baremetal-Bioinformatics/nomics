//! `samtools merge` (1.24, default options) of the per-group ApplyBQSR BAMs
//! (§9.3): header translation and the RG/PG tag rewrite. The interleaving of
//! records is done by the caller.
//!
//! Every input carries the same @RG and @PG IDs, so for inputs 2..n each ID
//! clashes and becomes `ID-%08lX` with `lrand48()` (`gen_unique_id`, redrawn
//! while taken). The generator is seeded with `time(NULL)` at startup
//! (`hts_srand48`), so the suffixes are part of the reference as run.
//! Draws are consumed per input in order: its @RG IDs, then its @PG IDs.

use std::collections::{HashMap, HashSet};

use nomics_io::bam::{RecordView, aux_value_len};

/// POSIX `srand48`/`lrand48`, as htslib's `hts_srand48`/`hts_lrand48`.
struct Rand48(u64);

impl Rand48 {
    fn seeded(seed: i64) -> Rand48 {
        Rand48(((seed as u64 & 0xffff_ffff) << 16) | 0x330e)
    }

    fn next(&mut self) -> u64 {
        self.0 = (0x5_deec_e66d_u64.wrapping_mul(self.0).wrapping_add(0xb)) & ((1 << 48) - 1);
        self.0 >> 17
    }
}

fn line_id(line: &str) -> Option<&str> {
    line.split('\t').find_map(|f| f.strip_prefix("ID:"))
}

fn replace_tag(line: &str, tag: &str, value: &str) -> String {
    line.split('\t')
        .map(|f| match f.strip_prefix(tag) {
            Some(_) => format!("{tag}{value}"),
            None => f.to_owned(),
        })
        .collect::<Vec<_>>()
        .join("\t")
}

/// The merged header text and, per input, the RG and PG ID maps.
pub struct Merged {
    pub header_text: String,
    pub rg: Vec<HashMap<Vec<u8>, Vec<u8>>>,
    pub pg: Vec<HashMap<Vec<u8>, Vec<u8>>>,
}

/// Merge `inputs` identical-but-for-@PG-CL headers. `texts[k]` is the
/// header text of the k-th input in merge order.
pub fn merge_headers(texts: &[String], seed: i64) -> Merged {
    let mut rng = Rand48::seeded(seed);
    let (mut known_rg, mut known_pg) = (HashSet::new(), HashSet::new());
    let (mut out_rg, mut out_pg, mut out_co) = (String::new(), String::new(), String::new());
    let mut head = String::new();
    let (mut rg_maps, mut pg_maps) = (Vec::new(), Vec::new());
    for (k, text) in texts.iter().enumerate() {
        let lines: Vec<&str> = text.trim_end_matches('\0').lines().collect();
        if k == 0 {
            for l in &lines {
                if l.starts_with("@HD") || l.starts_with("@SQ") {
                    head.push_str(l);
                    head.push('\n');
                }
            }
        }
        let mut translate = |kind: &str, known: &mut HashSet<String>| {
            let mut map: HashMap<String, String> = HashMap::new();
            let mut kept = Vec::new();
            for l in lines.iter().filter(|l| l.starts_with(kind)) {
                let Some(id) = line_id(l) else { continue };
                let new = if known.contains(id) {
                    loop {
                        let candidate = format!("{id}-{:08X}", rng.next());
                        if !known.contains(&candidate) {
                            break candidate;
                        }
                    }
                } else {
                    id.to_owned()
                };
                known.insert(new.clone());
                map.insert(id.to_owned(), new.clone());
                kept.push(replace_tag(l, "ID:", &new));
            }
            (map, kept)
        };
        let (rg_map, rg_lines) = translate("@RG", &mut known_rg);
        let (pg_map, pg_lines) = translate("@PG", &mut known_pg);
        // finish_rg_pg: @RG PG: and @PG PP: point at this input's new PG IDs.
        for (lines, tag, out) in [
            (rg_lines, "PG:", &mut out_rg),
            (pg_lines, "PP:", &mut out_pg),
        ] {
            for l in lines {
                let fixed = match l
                    .split('\t')
                    .find_map(|f| f.strip_prefix(tag))
                    .and_then(|v| pg_map.get(v))
                {
                    Some(v) => replace_tag(&l, tag, v),
                    None => l,
                };
                out.push_str(&fixed);
                out.push('\n');
            }
        }
        for l in lines.iter().filter(|l| l.starts_with("@CO")) {
            out_co.push_str(l);
            out_co.push('\n');
        }
        let bytes = |m: HashMap<String, String>| {
            m.into_iter()
                .map(|(a, b)| (a.into_bytes(), b.into_bytes()))
                .collect::<HashMap<_, _>>()
        };
        rg_maps.push(bytes(rg_map));
        pg_maps.push(bytes(pg_map));
    }
    Merged {
        header_text: head + &out_rg + &out_pg + &out_co,
        rg: rg_maps,
        pg: pg_maps,
    }
}

/// `bam_translate`: delete RG and append its mapped value at the end of the
/// tags, then the same for PG. A tag whose ID has no mapping is dropped.
/// `raw` is a record payload; the result is appended to `out` with its
/// `block_size` prefix.
pub fn translate_record(
    raw: &[u8],
    rg: &HashMap<Vec<u8>, Vec<u8>>,
    pg: &HashMap<Vec<u8>, Vec<u8>>,
    out: &mut Vec<u8>,
) {
    let aux_at = RecordView::new(raw)
        .map_or(raw.len(), |v| v.aux_offset())
        .min(raw.len());
    let start = out.len();
    out.extend_from_slice(&[0; 4]);
    out.extend_from_slice(&raw[..aux_at]);
    // Fields other than RG and PG are copied in order; a malformed tail
    // (an unparsable field) is dropped from there on.
    let (mut rg_val, mut pg_val) = (None, None);
    let mut aux = &raw[aux_at..];
    while aux.len() >= 3 {
        let (tag, ty) = ([aux[0], aux[1]], aux[2]);
        let Some(len) = aux_value_len(ty, &aux[3..]) else {
            break;
        };
        let field = &aux[..3 + len];
        match (&tag, ty) {
            (b"RG", b'Z') => rg_val = Some(&field[3..field.len() - 1]),
            (b"PG", b'Z') => pg_val = Some(&field[3..field.len() - 1]),
            _ => out.extend_from_slice(field),
        }
        aux = &aux[3 + len..];
    }
    for (tag, value, map) in [(b"RG", rg_val, rg), (b"PG", pg_val, pg)] {
        if let Some(new) = value.and_then(|v| map.get(v)) {
            out.extend_from_slice(tag);
            out.push(b'Z');
            out.extend_from_slice(new);
            out.push(0);
        }
    }
    let len = (out.len() - start - 4) as i32;
    out[start..start + 4].copy_from_slice(&len.to_le_bytes());
}

/// The `@PG` lines htslib's `sam_hdr_add_pg` adds for `samtools` (as
/// `samtools view -Ch` does making sarek's md.cram): one per chain end, an
/// unreferenced line whose `PP:` names another existing line, in header
/// order (with none, the last unreferenced line), each with `PP:` that end
/// and an ID made unique as `sam_hdr_pg_id` does (`samtools`, then
/// `samtools.1`, `.2`, ...). They go after the last `@PG` line.
pub fn add_samtools_pg(text: &str, version: &str, cl: &str) -> String {
    let text = text.trim_end_matches('\0');
    let lines: Vec<&str> = text.lines().collect();
    let pg: Vec<&str> = lines.iter().copied().filter(|l| l.starts_with("@PG\t")).collect();
    let field = |l: &str, key: &str| -> Option<String> {
        l.split('\t').skip(1).find_map(|f| f.strip_prefix(key)).map(str::to_owned)
    };
    let mut ids: HashSet<String> = pg.iter().filter_map(|l| field(l, "ID:")).collect();
    let linked = |l: &str| {
        field(l, "PP:").is_some_and(|pp| ids.contains(&pp) && Some(&pp) != field(l, "ID:").as_ref())
    };
    let pps: HashSet<String> = pg
        .iter()
        .filter(|l| linked(l))
        .filter_map(|l| field(l, "PP:"))
        .collect();
    let unreferenced: Vec<&str> = pg
        .iter()
        .copied()
        .filter(|l| field(l, "ID:").is_some_and(|id| !pps.contains(&id)))
        .collect();
    let mut ends: Vec<Option<String>> = unreferenced
        .iter()
        .filter(|l| linked(l))
        .filter_map(|l| field(l, "ID:"))
        .map(Some)
        .collect();
    if pg.is_empty() {
        ends.push(None);
    } else if ends.is_empty() {
        ends.extend(unreferenced.last().and_then(|l| field(l, "ID:")).map(Some));
    }
    let mut counter = 1;
    let mut added = Vec::new();
    for end in ends {
        let mut id = "samtools".to_owned();
        while ids.contains(&id) {
            id = format!("samtools.{counter}");
            counter += 1;
        }
        ids.insert(id.clone());
        let pp = end.map(|e| format!("\tPP:{e}")).unwrap_or_default();
        added.push(format!("@PG\tID:{id}\tPN:samtools{pp}\tVN:{version}\tCL:{cl}"));
    }
    let last_pg = lines.iter().rposition(|l| l.starts_with("@PG\t"));
    let mut out = String::with_capacity(text.len() + 200 * added.len());
    for (i, l) in lines.iter().enumerate() {
        out.push_str(l);
        out.push('\n');
        if Some(i) == last_pg {
            for a in &added {
                out.push_str(a);
                out.push('\n');
            }
        }
    }
    if last_pg.is_none() {
        for a in &added {
            out.push_str(a);
            out.push('\n');
        }
    }
    out
}

/// The `srand48` seeds whose first two draws are `first` and `second`
/// (`lrand48` outputs): the 17 bits each draw drops are searched, so this
/// recovers the seed samtools merge ran with from the suffixes it wrote.
fn seeds_for(first: u64, second: u64) -> Vec<i64> {
    const A: u64 = 0x5_deec_e66d;
    const C: u64 = 0xb;
    const MASK: u64 = (1 << 48) - 1;
    // A is odd, so it has an inverse mod 2^48 (Newton's iteration).
    let mut inv: u64 = A;
    for _ in 0..6 {
        inv = inv.wrapping_mul(2u64.wrapping_sub(A.wrapping_mul(inv)));
    }
    let inv = inv & MASK;
    let mut seeds = Vec::new();
    for low in 0..(1u64 << 17) {
        let s1 = (first << 17) | low;
        let s0 = (s1.wrapping_sub(C) & MASK).wrapping_mul(inv) & MASK;
        if s0 & 0xffff != 0x330e {
            continue;
        }
        let seed = (s0 >> 16) as i64;
        let mut r = Rand48::seeded(seed);
        if r.next() == first && r.next() == second {
            seeds.push(seed);
        }
    }
    seeds
}

/// How sarek's MERGE_CRAM renames read group and program IDs when ApplyBQSR
/// ran on md.cram (CRAM output) rather than md.bam: each group's header then
/// carries md.cram's extra `samtools view` `@PG` lines, whose collisions
/// take more `lrand48` draws, so from the third input on the suffixes differ
/// from the BAM-mode merge (recal.bam) that our records carry.
pub struct CramMerge {
    /// The merged header as samtools merge wrote it in CRAM mode.
    pub header_text: String,
    rg: HashMap<Vec<u8>, Vec<u8>>,
    pg: HashMap<Vec<u8>, Vec<u8>>,
}

impl CramMerge {
    /// From the BAM-mode merged header `recal_text` of `groups` identical
    /// ApplyBQSR outputs (BQSR's own output stream): input 0's lines are
    /// the group header; the merge seed is recovered from input 1's first
    /// two suffixes; and the BAM-mode merge is recomputed and checked
    /// against `recal_text` before the CRAM-mode one (the group header plus
    /// `samtools_cl`'s `@PG` lines, before ApplyBQSR's) is used.
    pub fn from_recal_header(recal_text: &str, groups: usize, samtools_cl: &str) -> Result<CramMerge, String> {
        let text = recal_text.trim_end_matches('\0');
        let lines: Vec<&str> = text.lines().collect();
        let of = |kind: &str| -> Vec<&str> { lines.iter().copied().filter(|l| l.starts_with(kind)).collect() };
        let (rg, pg, co) = (of("@RG\t"), of("@PG\t"), of("@CO\t"));
        let bad = |why: &str| format!("recal header is not a samtools merge of {groups} ApplyBQSR groups: {why}");
        if groups < 2 || rg.len() % groups != 0 || pg.len() % groups != 0 || co.len() % groups != 0 {
            return Err(bad("line counts"));
        }
        let (nrg, npg, nco) = (rg.len() / groups, pg.len() / groups, co.len() / groups);
        let head: Vec<&str> = lines
            .iter()
            .copied()
            .filter(|l| l.starts_with("@HD") || l.starts_with("@SQ"))
            .collect();
        let group = |pg_lines: &[String]| -> String {
            let mut t = String::new();
            for l in head.iter().chain(&rg[..nrg]) {
                t.push_str(l);
                t.push('\n');
            }
            for l in pg_lines {
                t.push_str(l);
                t.push('\n');
            }
            for l in &co[..nco] {
                t.push_str(l);
                t.push('\n');
            }
            t
        };
        let pg0: Vec<String> = pg[..npg].iter().map(|l| l.to_string()).collect();
        let bam_group = group(&pg0);
        // Input 1's @RG suffix is the first draw, its first @PG's the second.
        let suffix = |l: &str| -> Option<u64> {
            let id = line_id(l)?;
            u64::from_str_radix(id.rsplit_once('-')?.1, 16).ok()
        };
        let (Some(d1), Some(d2)) = (suffix(rg[nrg]), suffix(pg[npg])) else {
            return Err(bad("no suffixed IDs on the second input"));
        };
        let ids = |t: &str| -> Vec<String> {
            t.lines()
                .filter(|l| l.starts_with("@RG\t") || l.starts_with("@PG\t"))
                .filter_map(line_id)
                .map(str::to_owned)
                .collect()
        };
        let want = ids(text);
        let seed = seeds_for(d1, d2)
            .into_iter()
            .find(|&s| ids(&merge_headers(&vec![bam_group.clone(); groups], s).header_text) == want)
            .ok_or_else(|| bad("no seed reproduces its IDs"))?;
        let bam = merge_headers(&vec![bam_group; groups], seed);
        // md.cram's header: the group's @PG lines before ApplyBQSR's (the
        // last), then samtools view's, then ApplyBQSR's.
        let (bqsr, md_pg) = pg0.split_last().ok_or_else(|| bad("no @PG lines"))?;
        if !bqsr.starts_with("@PG\tID:GATK ApplyBQSR") {
            return Err(bad("the last @PG of a group is not ApplyBQSR's"));
        }
        let md_cram = add_samtools_pg(&group(md_pg), "1.21", samtools_cl);
        let cram_group = format!("{md_cram}{bqsr}\n");
        let cram = merge_headers(&vec![cram_group; groups], seed);
        let pair = |a: &[HashMap<Vec<u8>, Vec<u8>>], b: &[HashMap<Vec<u8>, Vec<u8>>]| {
            let mut m = HashMap::new();
            for (ka, kb) in a.iter().zip(b) {
                for (id, new) in ka {
                    if let Some(to) = kb.get(id) {
                        m.insert(new.clone(), to.clone());
                    }
                }
            }
            m
        };
        Ok(CramMerge {
            rg: pair(&bam.rg, &cram.rg),
            pg: pair(&bam.pg, &cram.pg),
            header_text: cram.header_text,
        })
    }

    /// A record of the BAM-mode merge with its RG:Z and PG:Z renamed to
    /// the CRAM-mode merge's (moved to the end, as `bam_translate` leaves
    /// them), `block_size`-prefixed, appended to `out`.
    pub fn translate(&self, raw: &[u8], out: &mut Vec<u8>) {
        translate_record(raw, &self.rg, &self.pg, out);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rand48_matches_the_reference_suffixes() {
        // sarek wgs20-20261010-r1: MERGE_BAM started at 1791621718; input 2's
        // RG got 3AF61EBD and its first @PG (bwa) 646C4E63.
        let mut r = Rand48::seeded(1_791_621_718);
        assert_eq!(format!("{:08X}", r.next()), "3AF61EBD");
        assert_eq!(format!("{:08X}", r.next()), "646C4E63");
    }

    #[test]
    fn second_input_ids_get_suffixes() {
        let h = "@HD\tVN:1.6\n@SQ\tSN:c\tLN:9\n@RG\tID:g\tSM:s\n@PG\tID:a\tPN:a\n@PG\tID:b\tPP:a\n";
        let m = merge_headers(&[h.into(), h.into()], 1_791_621_718);
        assert!(m.header_text.contains("@RG\tID:g\tSM:s\n@RG\tID:g-3AF61EBD\tSM:s\n"));
        assert!(m.header_text.contains("@PG\tID:b-413FD371\tPP:a-646C4E63\n"));
        assert_eq!(m.rg[1][&b"g"[..]], b"g-3AF61EBD");
    }

    #[test]
    fn seed_recovered_from_two_draws() {
        let mut r = Rand48::seeded(1_791_621_718);
        let (a, b) = (r.next(), r.next());
        assert!(seeds_for(a, b).contains(&1_791_621_718));
    }

    #[test]
    fn samtools_pg_lines_per_chain_end() {
        let h = "@HD\tVN:1.6\n@PG\tID:a\tPN:a\n@PG\tID:a.1\tPN:a\n@PG\tID:s\tPP:a\n@PG\tID:s.1\tPP:a.1\n@PG\tID:MD\tVN:1\n@CO\tx\n";
        let out = add_samtools_pg(h, "1.21", "samtools view");
        assert_eq!(
            out,
            "@HD\tVN:1.6\n@PG\tID:a\tPN:a\n@PG\tID:a.1\tPN:a\n@PG\tID:s\tPP:a\n@PG\tID:s.1\tPP:a.1\n@PG\tID:MD\tVN:1\n\
             @PG\tID:samtools\tPN:samtools\tPP:s\tVN:1.21\tCL:samtools view\n\
             @PG\tID:samtools.1\tPN:samtools\tPP:s.1\tVN:1.21\tCL:samtools view\n@CO\tx\n"
        );
    }

    /// The CRAM-mode merge, from a BAM-mode merged header, gives the IDs a
    /// samtools merge of the CRAM-mode group headers gives, with the
    /// records' IDs renamed to match.
    #[test]
    fn cram_merge_matches_merging_cram_mode_headers() {
        let md = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c\tLN:9\n@RG\tID:g\tSM:s\n\
                  @PG\tID:bwa\tPN:bwa\n@PG\tID:bwa.1\tPN:bwa\n@PG\tID:sort\tPP:bwa\n@PG\tID:sort.1\tPP:bwa.1\n\
                  @PG\tID:MarkDuplicates\tVN:1\n";
        let bqsr = "@PG\tID:GATK ApplyBQSR\tVN:4.6.2.0\tCL:x\n";
        let seed = 1_791_676_187;
        let n = 4;
        let bam = merge_headers(&vec![format!("{md}{bqsr}"); n], seed);
        let cram_group = format!("{}{bqsr}", add_samtools_pg(md, "1.21", "sv"));
        let want = merge_headers(&vec![cram_group; n], seed);
        let m = CramMerge::from_recal_header(&bam.header_text, n, "sv").unwrap();
        let ids = |t: &str| -> Vec<String> {
            t.lines().filter_map(line_id).map(str::to_owned).collect()
        };
        assert_eq!(ids(&m.header_text), ids(&want.header_text));
        assert_ne!(ids(&m.header_text)[..n], ids(&bam.header_text)[..n], "suffixes differ by input 3");
        // A record of input 3 under BAM-mode IDs gets the CRAM-mode ones.
        let (rg, pg) = (&bam.rg[3][&b"g"[..]], &bam.pg[3][&b"MarkDuplicates"[..]]);
        let mut raw = vec![0u8; 32];
        raw[8] = 2; // l_read_name
        raw.extend_from_slice(b"r\0");
        for (tag, v) in [(b"RG", rg), (b"PG", pg)] {
            raw.extend_from_slice(tag);
            raw.push(b'Z');
            raw.extend_from_slice(v);
            raw.push(0);
        }
        let mut out = Vec::new();
        m.translate(&raw, &mut out);
        let tail = |id: &[u8], tag: &[u8]| [tag, b"Z", id, b"\0"].concat();
        let want_rg = tail(&want.rg[3][&b"g"[..]], b"RG");
        let want_pg = tail(&want.pg[3][&b"MarkDuplicates"[..]], b"PG");
        assert!(out.ends_with(&[want_rg, want_pg].concat()));
    }

    /// The 30x cloud gate's headers (when present): ours (BAM-mode merge of
    /// 21 groups), the reference recal.cram's.
    #[test]
    fn cram_merge_reproduces_the_30x_reference_ids() {
        let Ok(dir) = std::env::var("NOMICS_CLOUD_RECAL_GATE") else {
            return;
        };
        let read = |f: &str| std::fs::read_to_string(format!("{dir}/{f}")).unwrap();
        // The dumps were made with `samtools view -H` (which adds its own
        // @PG lines), and samtools merge adds one per chain end after
        // renaming; neither takes part in the renaming.
        let strip = |t: String| -> String {
            t.lines()
                .filter(|l| {
                    !l.starts_with("@PG\tID:nomics-cram")
                        && !l.contains("\tCL:samtools view -H")
                        && !l.contains("\tCL:samtools merge")
                })
                .map(|l| format!("{l}\n"))
                .collect()
        };
        let ours = strip(read("header.ours_recal.txt"));
        let m = CramMerge::from_recal_header(&ours, 21, "samtools view").unwrap();
        let ids = |t: &str, kind: &str| -> Vec<String> {
            t.lines().filter(|l| l.starts_with(kind)).filter_map(line_id).map(str::to_owned).collect()
        };
        let reference = strip(read("header.ref_recal.txt"));
        assert_eq!(ids(&m.header_text, "@RG"), ids(&reference, "@RG"));
        assert_eq!(ids(&m.header_text, "@PG"), ids(&reference, "@PG"));
    }
}
