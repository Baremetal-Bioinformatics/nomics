//! GATK 4.6.2.0 GetPileupSummaries, default sarek command. The sparse
//! traversal accumulates only resource sites, while retaining deletion/N-base
//! coverage separately from A/C/G/T counts (LocusIteratorByState semantics).
use crate::pileup_summary::{self, PileupSummary};
use nomics_io::{
    bam::{Header, RecordView},
    bam_query::BamQuery,
    bgzf::BgzfReader,
    fasta::read_fai,
    tabix::Tabix,
    vcf::{VcfRecord, parse_line},
};
use std::{
    fs::File,
    io::{self, BufRead, BufReader},
    path::Path,
};

const FILTERS: [&str; 11] = [
    "MappingQualityReadFilter",
    "MappingQualityAvailableReadFilter",
    "MappingQualityNotZeroReadFilter",
    "MappedReadFilter",
    "PrimaryLineReadFilter",
    "NotDuplicateReadFilter",
    "PassesVendorQualityCheckReadFilter",
    "NonZeroReferenceLengthAlignmentReadFilter",
    "MateOnSameContigOrNoMappedMateReadFilter",
    "GoodCigarReadFilter",
    "WellformedReadFilter",
];
#[derive(Debug, Default)]
pub struct Stats {
    pub processed: u64,
    pub filtered: [u64; 11],
    pub rows: usize,
}
impl Stats {
    pub fn log(&self) {
        for (name, count) in FILTERS.iter().zip(self.filtered) {
            eprintln!("{count} read(s) filtered by: {name}");
        }
        eprintln!(
            "{} total reads filtered out of {} reads processed; {} rows",
            self.filtered.iter().sum::<u64>(),
            self.processed,
            self.rows
        );
    }
}
#[derive(Clone, Copy, Debug)]
struct Interval {
    beg: u32,
    end: u32,
} // zero-based half-open
fn intervals(path: &Path, header: &Header) -> io::Result<Vec<Vec<Interval>>> {
    let mut out = vec![Vec::<Interval>::new(); header.refs.len()];
    for line in BufReader::new(File::open(path)?).lines() {
        let line = line?;
        if line.is_empty()
            || line.starts_with('#')
            || line.starts_with("track")
            || line.starts_with("browser")
        {
            continue;
        }
        let f: Vec<_> = line.split_whitespace().collect();
        if f.len() < 3 {
            return Err(io::Error::other("short BED row"));
        }
        let tid = header
            .refs
            .iter()
            .position(|r| r.name == f[0].as_bytes())
            .ok_or_else(|| io::Error::other(format!("unknown BED contig {}", f[0])))?;
        let beg = f[1].parse::<u32>().map_err(io::Error::other)?;
        let end = f[2].parse::<u32>().map_err(io::Error::other)?;
        if beg >= end || end > header.refs[tid].len {
            return Err(io::Error::other("BED interval out of bounds or empty"));
        }
        out[tid].push(Interval { beg, end });
    }
    for contig in &mut out {
        contig.sort_unstable_by_key(|i| i.beg);
        let mut merged: Vec<Interval> = Vec::new();
        for i in contig.drain(..) {
            if let Some(last) = merged.last_mut().filter(|last| i.beg <= last.end) {
                last.end = last.end.max(i.end);
            } else {
                merged.push(i);
            }
        }
        *contig = merged;
    }
    Ok(out)
}
fn overlaps(intervals: &[Interval], beg: i64, end: i64) -> bool {
    let at = intervals.partition_point(|i| i64::from(i.end) <= beg);
    intervals.get(at).is_some_and(|i| i64::from(i.beg) < end)
}

#[derive(Debug)]
struct Site {
    locus: u32, // zero based traversal position; output uses VariantContext.start
    position: i32,
    reference: u8,
    alternate: u8,
    af: Option<String>,
    covered: bool,
    counts: [i32; 4],
}
fn info<'a>(vc: &VcfRecord<'a>, key: &[u8]) -> Option<&'a [u8]> {
    // AbstractVCFCodec.parseInfo is a map: the last occurrence wins.
    vc.info
        .split(|&b| b == b';')
        .filter_map(|f| {
            let eq = f.iter().position(|&b| b == b'=');
            match eq {
                Some(eq) if &f[..eq] == key => Some(&f[eq + 1..]),
                None if f == key => Some(b".".as_slice()),
                _ => None,
            }
        })
        .next_back()
}
fn snp_alt(vc: &VcfRecord<'_>) -> Option<u8> {
    // AbstractVCFCodec ignores no-call ALT entries. '*' is a nonsymbolic
    // length-one allele in htsjdk and BaseUtils counts that wildcard as A.
    if vc.ref_allele.len() != 1 {
        return None;
    }
    let mut alleles = vc.alt.split(|&b| b == b',').filter(|&a| a != b".");
    let alt = alleles.next()?;
    (alleles.next().is_none() && alt.len() == 1).then(|| alt[0])
}
fn is_biallelic_snp(vc: &VcfRecord<'_>) -> bool {
    snp_alt(vc).is_some()
}
fn af_value(raw: &str) -> io::Result<f64> {
    // CommonInfo.getAttributeAsDouble casts a multi-valued List to String:
    // malformed scalar/multi-valued AF is fatal, never take its first member.
    let lower = raw.to_ascii_lowercase();
    match lower.as_str() {
        "nan" | "+nan" | "-nan" => Ok(f64::NAN),
        "inf" | "+inf" | "infinity" | "+infinity" => Ok(f64::INFINITY),
        "-inf" | "-infinity" => Ok(f64::NEG_INFINITY),
        _ => raw.parse().map_err(io::Error::other),
    }
}
fn in_af_range(af: f64) -> bool {
    0.01 < af && af < 0.2
}
fn base_index(base: u8) -> Option<usize> {
    match base.to_ascii_uppercase() {
        b'A' | b'*' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}
fn sites(tabix: &Tabix, contig: &[u8], intervals: &[Interval]) -> io::Result<Vec<Site>> {
    let mut sites = Vec::new();
    // Query each occupied 1-Mb window once rather than materializing the
    // whole resource. Coarse TBI chunks may also decode flanking records.
    let mut windows = Vec::<u32>::new();
    for i in intervals {
        for window in i.beg / 1_000_000..=(i.end - 1) / 1_000_000 {
            if windows.last() != Some(&window) {
                windows.push(window);
            }
        }
    }
    for window in windows {
        let beg = window * 1_000_000;
        let end = ((window + 1) * 1_000_000).min(1 << 29);
        let mut previous_end = 0;
        for line in tabix.vcf_lines_overlapping(contig, u64::from(beg) + 1, u64::from(end))? {
            let vc = parse_line(&line)?;
            let stop = info(&vc, b"END")
                .map(|s| {
                    std::str::from_utf8(s)
                        .map_err(io::Error::other)?
                        .parse::<u64>()
                        .map_err(io::Error::other)
                })
                .transpose()?
                .unwrap_or(vc.pos + vc.ref_allele.len() as u64 - 1);
            let start = vc.pos.max(previous_end + 1).max(u64::from(beg) + 1);
            previous_end = previous_end.max(stop);
            if !is_biallelic_snp(&vc) {
                continue;
            }
            for position in start..=stop.min(u64::from(end)) {
                let locus = (position - 1) as u32;
                if overlaps(intervals, i64::from(locus), i64::from(locus) + 1) {
                    sites.push(Site {
                        locus,
                        position: vc.pos as i32,
                        reference: vc.ref_allele[0],
                        alternate: snp_alt(&vc).unwrap(),
                        af: info(&vc, b"AF").map(|s| String::from_utf8_lossy(s).into_owned()),
                        covered: false,
                        counts: [0; 4],
                    });
                }
            }
        }
    }
    Ok(sites)
}

// Fresh local counterparts of HC's private filter helpers. In particular,
// AlignmentAgreesWithHeader checks START, not end; reads crossing the contig
// end are accepted by GATK (ReadUtils.alignmentAgreesWithHeader).
fn good_cigar(cigar: &[u32]) -> bool {
    let n = cigar.len();
    let op = |i: usize| cigar[i] & 15;
    let real = |o| matches!(o, 0..=3 | 7 | 8);
    let mut seen_real = false;
    for (i, &c) in cigar.iter().enumerate() {
        if c >> 4 == 0 || op(i) > 8 {
            return false;
        }
        match op(i) {
            5 if i != 0 && i + 1 != n => return false,
            4 if i != 0
                && i + 1 != n
                && !(i == 1 && op(0) == 5)
                && !(i + 2 == n && op(n - 1) == 5) =>
            {
                return false;
            }
            6 if i > 0 && (i + 1 == n || !real(op(i - 1)) || !real(op(i + 1))) => return false,
            o if real(o) => {
                seen_real = true;
                if matches!(o, 1 | 2) {
                    for &next in &cigar[i + 1..] {
                        let next = next & 15;
                        if (real(next) && !matches!(next, 1 | 2)) || next == 6 {
                            break;
                        }
                        if next == o {
                            return false;
                        }
                    }
                }
            }
            _ => {}
        }
        if i > 0 && matches!(op(i), 1 | 2) && matches!(op(i - 1), 1 | 2) {
            return false;
        }
    }
    let mut core = cigar.iter().map(|c| c & 15).filter(|&o| o != 4 && o != 5);
    let first = core.next();
    let last = core.next_back().or(first);
    seen_real && first != Some(2) && last != Some(2)
}
fn unmapped(r: &RecordView<'_>) -> bool {
    r.flag() & 4 != 0 || r.ref_id() < 0 || r.pos() < 0
}
fn mate_unmapped(r: &RecordView<'_>) -> bool {
    r.flag() & 8 != 0 || r.mate_ref_id() < 0 || r.mate_pos() < 0
}
fn read_filter(r: &RecordView<'_>, cigar: &[u32], header: &Header) -> Option<usize> {
    let flags = r.flag();
    let reference_length = r.reference_length();
    let read_length: u64 = cigar
        .iter()
        .filter(|c| matches!(*c & 15, 0 | 1 | 4 | 7 | 8))
        .map(|c| u64::from(c >> 4))
        .sum();
    let passes = [
        r.mapq() >= 50,
        r.mapq() != 255,
        r.mapq() != 0,
        !unmapped(r),
        flags & (0x100 | 0x800) == 0,
        flags & 0x400 == 0,
        flags & 0x200 == 0,
        reference_length > 0,
        flags & 1 == 0 || mate_unmapped(r) || (!unmapped(r) && r.ref_id() == r.mate_ref_id()),
        good_cigar(cigar),
        r.pos() >= 0
            && header
                .refs
                .get(r.ref_id() as usize)
                .is_some_and(|c| (r.pos() as u32) < c.len)
            && r.string_tag(b"RG").is_some()
            && r.quals()
                .is_some_and(|q| q.len() == r.len() && q.first() != Some(&255))
            && r.len() as u64 == read_length
            && !r.is_empty()
            && !cigar.iter().any(|c| c & 15 == 3),
    ];
    passes.iter().position(|&p| !p)
}
fn inside_adaptor(r: &RecordView<'_>, locus: i64) -> bool {
    // ReadUtils.isBaseInsideAdaptor: the >100 shortcut uses SIGNED TLEN.
    let f = r.flag();
    if r.tlen() == 0
        || r.tlen() > 100
        || f & 1 == 0
        || unmapped(r)
        || mate_unmapped(r)
        || (f & 16 != 0) == (f & 32 != 0)
    {
        return false;
    }
    let start = i64::from(r.pos()) + 1;
    let mate_start = i64::from(r.mate_pos()) + 1;
    let position = locus + 1;
    if f & 16 != 0 {
        let end = start + i64::from(r.reference_length()) - 1;
        end > mate_start && position < mate_start
    } else {
        start <= mate_start + i64::from(r.tlen())
            && position >= start + i64::from(r.tlen().wrapping_abs())
    }
}
fn count_read(r: &RecordView<'_>, cigar: &[u32], sites: &mut [Site]) {
    let mut pos = i64::from(r.pos());
    let mut read_offset = 0;
    for &c in cigar {
        let length = i64::from(c >> 4);
        match c & 15 {
            op @ (0 | 2 | 7 | 8) => {
                let at = sites.partition_point(|s| i64::from(s.locus) < pos);
                for site in sites[at..]
                    .iter_mut()
                    .take_while(|s| i64::from(s.locus) < pos + length)
                {
                    if inside_adaptor(r, i64::from(site.locus)) {
                        continue;
                    }
                    site.covered = true;
                    if op != 2 {
                        let base = r.base((read_offset + i64::from(site.locus) - pos) as usize);
                        if let Some(index) = base_index(base) {
                            site.counts[index] = site.counts[index].wrapping_add(1);
                        }
                    }
                }
                pos += length;
                if op != 2 {
                    read_offset += length;
                }
            }
            1 | 4 => read_offset += length,
            3 => pos += length,
            _ => {}
        }
    }
}

/// Run with the defaults used by sarek. Reference bases are never consulted by
/// this tool; use the FASTA's .fai to validate the traversal dictionary.
pub fn run(
    input: &Path,
    variant: &Path,
    reference: &Path,
    bed: &Path,
    output: &Path,
) -> io::Result<Stats> {
    let mut header_reader = BgzfReader::new(BufReader::new(File::open(input)?));
    let header = Header::read(&mut header_reader)?;
    File::open(reference)?;
    let mut fai_path = reference.as_os_str().to_owned();
    fai_path.push(".fai");
    let fai = read_fai(Path::new(&fai_path))?;
    for r in &header.refs {
        if let Some(c) = fai.iter().find(|c| c.name.as_bytes() == r.name)
            && c.len != u64::from(r.len)
        {
            return Err(io::Error::other("BAM/reference dictionary length mismatch"));
        }
    }
    let targets = intervals(bed, &header)?;
    let sample = String::from_utf8_lossy(&header.text)
        .lines()
        .filter(|l| l.starts_with("@RG\t"))
        .filter_map(|l| l.split('\t').find_map(|f| f.strip_prefix("SM:")))
        .map(str::to_owned)
        .min_by(|a, b| a.encode_utf16().cmp(b.encode_utf16()))
        .ok_or_else(|| io::Error::other("no sample in BAM header"))?;
    let mut header_vcf = BufReader::new(BgzfReader::new(BufReader::new(File::open(variant)?)));
    let mut line = String::new();
    let mut has_af = false;
    while header_vcf.read_line(&mut line)? > 0 {
        if !line.starts_with('#') {
            break;
        }
        if let Some(fields) = line.strip_prefix("##INFO=<") {
            has_af |= fields.split(',').any(|f| f == "ID=AF");
        }
        line.clear();
    }
    if !has_af {
        return Err(io::Error::other(
            "Population vcf does not have an allele frequency (AF) info field in its header.",
        ));
    }
    let tabix = Tabix::open(variant)?;
    let queries: Vec<_> = targets
        .iter()
        .enumerate()
        .flat_map(|(tid, v)| v.iter().map(move |i| (tid, i.beg, i.end)))
        .collect();
    let mut bam = BamQuery::open(input, &queries)?;
    let mut stats = Stats::default();
    let mut writer = std::io::BufWriter::new(File::create(output)?);
    pileup_summary::write_header(&mut writer, Some(&sample))?;
    let mut raw = Vec::new();
    let mut pending = bam.next_record(&mut raw)?;
    let (mut saw_without, mut saw_with) = (false, false);
    for (tid, target) in targets.iter().enumerate() {
        let mut sites = sites(&tabix, &header.refs[tid].name, target)?;
        while pending {
            let r = RecordView::new(&raw)?;
            if r.ref_id() > tid as i32 {
                break;
            }
            if r.ref_id() == tid as i32
                && overlaps(
                    target,
                    i64::from(r.pos()),
                    i64::from(r.pos()) + i64::from(r.reference_length().max(1)),
                )
            {
                stats.processed += 1;
                let cigar: Vec<_> = r.cigar().collect();
                if let Some(filter) = read_filter(&r, &cigar, &header) {
                    stats.filtered[filter] += 1;
                } else {
                    count_read(&r, &cigar, &mut sites);
                }
            }
            pending = bam.next_record(&mut raw)?;
        }
        for site in sites.into_iter().filter(|s| s.covered) {
            let Some(af) = site.af else {
                if !saw_without {
                    eprintln!(
                        "WARNING: Variant context at {}:{} lacks allele frequency (AF) field.",
                        String::from_utf8_lossy(&header.refs[tid].name),
                        site.position
                    );
                }
                saw_without = true;
                continue;
            };
            saw_with = true;
            let af = af_value(&af)?;
            if !in_af_range(af) {
                continue;
            }
            let ref_index = base_index(site.reference)
                .ok_or_else(|| io::Error::other("non-simple resource REF base"))?;
            let alt_index = base_index(site.alternate)
                .ok_or_else(|| io::Error::other("non-simple resource ALT base"))?;
            let ref_count = site.counts[ref_index];
            let alt_count = site.counts[alt_index];
            let total = site.counts.iter().fold(0i32, |a, b| a.wrapping_add(*b));
            pileup_summary::write_record(
                &mut writer,
                &PileupSummary {
                    contig: String::from_utf8_lossy(&header.refs[tid].name).into_owned(),
                    position: site.position,
                    ref_count,
                    alt_count,
                    other_alt_count: total.wrapping_sub(ref_count).wrapping_sub(alt_count),
                    allele_frequency: af,
                },
            )?;
            stats.rows += 1;
        }
    }
    use std::io::Write;
    writer.flush()?;
    if saw_without && !saw_with {
        return Err(io::Error::other(
            "No variants in population vcf had an allele frequency (AF) field.",
        ));
    }
    Ok(stats)
}

#[cfg(test)]
mod tests {
    use super::*;
    use nomics_io::bam::{BinRule, Record, Reference, Value};
    fn record(seq: &[u8], cigar: &[u32]) -> Record {
        Record {
            ref_id: 0,
            pos: 9,
            mapq: 60,
            flag: 0,
            next_ref_id: -1,
            next_pos: -1,
            tlen: 0,
            read_name: b"same".to_vec(),
            cigar: cigar.to_vec(),
            seq: seq.to_vec(),
            qual: vec![0; seq.len()],
            tags: vec![(*b"RG", Value::String(b"rg".to_vec()))],
        }
    }
    fn site(locus: u32) -> Site {
        Site {
            locus,
            position: locus as i32 + 1,
            reference: b'A',
            alternate: b'C',
            af: Some("0.1".into()),
            covered: false,
            counts: [0; 4],
        }
    }
    fn add(r: &Record, sites: &mut [Site]) {
        let mut raw = Vec::new();
        r.encode(BinRule::Htsjdk, &mut raw).unwrap();
        count_read(&RecordView::new(&raw).unwrap(), &r.cigar, sites);
    }
    #[test]
    fn java_site_selection_and_scalar_af_rules() {
        // GetPileupSummaries.apply/alleleFrequencyInRange, htsjdk
        // VariantContext.typeOfBiallelicVariant and CommonInfo.getAttributeAsDouble.
        for (reference, alt, expected) in [
            ("A", "C", true),
            ("A", "C,G", false),
            ("AA", "CC", false),
            ("A", "AC", false),
            ("A", "<DEL>", false),
            ("A", "*", true),
            ("A", ".", false),
            ("N", "C", true),
            ("A", "C,.", true),
        ] {
            let line = format!("c\t1\t.\t{reference}\t{alt}\t.\tFAIL\tAF=0.1");
            assert_eq!(
                is_biallelic_snp(&parse_line(line.as_bytes()).unwrap()),
                expected
            );
        }
        for (value, pass) in [
            ("0.01", false),
            ("0.0100000000001", true),
            ("0.199999999", true),
            ("0.2", false),
            ("NaN", false),
            ("-Inf", false),
            ("Infinity", false),
        ] {
            assert_eq!(in_af_range(af_value(value).unwrap()), pass);
        }
        for value in [".", "", "0.1,0.2"] {
            assert!(af_value(value).is_err());
        }
        let absent = parse_line(b"c\t1\t.\tA\tC\t.\t.\t.").unwrap();
        assert_eq!(info(&absent, b"AF"), None);
        let repeated = parse_line(b"c\t1\t.\tA\tC\t.\t.\tAF=0.1;AF=0.2").unwrap();
        assert_eq!(info(&repeated, b"AF"), Some(b"0.2".as_slice()));
    }
    #[test]
    fn java_cigar_and_base_counts_without_quality_or_mate_collapse() {
        // AlignmentStateMachine.stepForwardOnGenome skips I/S/H/P, emits D.
        // ReadPileup.getBaseCounts skips D and non-simple bases. No BQ filter.
        let mut s: Vec<_> = (8..=14).map(site).collect();
        add(&record(b"AC", &[1 << 4, (2 << 4) | 2, 1 << 4]), &mut s);
        assert!(!s[0].covered);
        assert!(!s[5].covered);
        assert_eq!(s[1].counts, [1, 0, 0, 0]);
        assert!(s[2].covered && s[3].covered);
        assert_eq!(s[2].counts, [0; 4]);
        assert_eq!(s[3].counts, [0; 4]);
        assert_eq!(s[4].counts, [0, 1, 0, 0]);
        let mut s = vec![site(9), site(10), site(11)];
        add(
            &record(
                b"TACGA",
                &[(1 << 4) | 4, 1 << 4, (1 << 4) | 1, 1 << 4, (1 << 4) | 4],
            ),
            &mut s,
        );
        assert_eq!(s[0].counts, [1, 0, 0, 0]);
        assert_eq!(s[1].counts, [0, 0, 1, 0]);
        assert!(!s[2].covered);
        // Two reads with the same name, each a different end: both count.
        let mut a = record(b"C", &[1 << 4]);
        a.flag = 0x41;
        let mut b = a.clone();
        b.flag = 0x81;
        let mut s = vec![site(9)];
        add(&a, &mut s);
        add(&b, &mut s);
        assert_eq!(s[0].counts, [0, 2, 0, 0]);
        let mut s = vec![site(9), site(10)];
        add(&record(b"NG", &[2 << 4]), &mut s);
        assert!(s[0].covered);
        assert_eq!(s[0].counts, [0; 4]);
        assert!(s[1].covered);
        assert_eq!(s[1].counts, [0, 0, 1, 0]); // other-alt only
    }
    #[test]
    fn java_adaptor_signed_tlen_and_filter_order_and_contig_end() {
        // ReadUtils.{hasWellDefinedFragmentSize,getAdaptorBoundary,isBaseInsideAdaptor}.
        let mut r = record(&[b'A'; 200], &[200 << 4]);
        r.flag = 1 | 32;
        r.next_ref_id = 0;
        r.next_pos = 20;
        r.tlen = 50;
        let raw = |r: &Record| {
            let mut b = Vec::new();
            r.encode(BinRule::Htsjdk, &mut b).unwrap();
            b
        };
        let bytes = raw(&r);
        let v = RecordView::new(&bytes).unwrap();
        assert!(!inside_adaptor(&v, 58));
        assert!(inside_adaptor(&v, 59));
        r.tlen = 101;
        let bytes = raw(&r);
        assert!(!inside_adaptor(&RecordView::new(&bytes).unwrap(), 199));
        r.flag = 1 | 16;
        r.next_pos = 20;
        r.tlen = -200;
        let bytes = raw(&r);
        let v = RecordView::new(&bytes).unwrap();
        assert!(inside_adaptor(&v, 19));
        assert!(!inside_adaptor(&v, 20));
        let h = Header {
            text: vec![],
            refs: vec![Reference {
                name: b"c".to_vec(),
                len: 100,
            }],
        };
        assert_eq!(read_filter(&v, &r.cigar, &h), None); // end 209 exceeds 100, start is valid
        r.flag |= 0x800 | 0x400;
        let bytes = raw(&r);
        assert_eq!(
            read_filter(&RecordView::new(&bytes).unwrap(), &r.cigar, &h),
            Some(4)
        );
        r.mapq = 0;
        let bytes = raw(&r);
        assert_eq!(
            read_filter(&RecordView::new(&bytes).unwrap(), &r.cigar, &h),
            Some(0)
        );
        // htsjdk BAMRecord.decodeBaseQualities maps first-byte FF to an
        // empty quality array; Wellformed then rejects the length mismatch.
        r.flag = 0;
        r.mapq = 60;
        r.qual.fill(255);
        let bytes = raw(&r);
        assert_eq!(
            read_filter(&RecordView::new(&bytes).unwrap(), &r.cigar, &h),
            Some(10)
        );
        for cigar in [
            vec![0],
            vec![1 << 4, 1 << 4 | 2, 1 << 4 | 1, 1 << 4],
            vec![1 << 4, 1 << 4 | 6],
            vec![1 << 4, 1 << 4 | 4, 1 << 4],
        ] {
            assert!(!good_cigar(&cigar));
        }
        assert!(good_cigar(&[1 << 4 | 6, 1 << 4])); // leading padding allowed by htsjdk
    }
}
