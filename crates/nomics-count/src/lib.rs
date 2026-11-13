//! Fresh implementation of the Subread 2.1.1 featureCounts path used by
//! nf-core/rnaseq 3.27.0. See `docs/FEATURECOUNTS.md` for source anchors and scope.
//!
//! In 2.1.1 `-p` declares paired input; **only `--countReadPairs` enables fragment
//! counting**. nf-core does not pass that option. Each BAM alignment is counted
//! separately, including duplicate and secondary records (unless NH > 1).
//! This crate deliberately exposes the alignment-counting path, not the other
//! featureCounts modes. Integer accumulation is independent of inflate threads.

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use nomics_io::bam::{RecordView, aux_value_len};
use nomics_io::bamstream::BamStream;

pub const VERSION: &str = "2.1.1";
pub const STATUSES: [&str; 14] = [
    "Assigned",
    "Unassigned_Unmapped",
    "Unassigned_Read_Type",
    "Unassigned_Singleton",
    "Unassigned_MappingQuality",
    "Unassigned_Chimera",
    "Unassigned_FragmentLength",
    "Unassigned_Duplicate",
    "Unassigned_MultiMapping",
    "Unassigned_Secondary",
    "Unassigned_NonSplit",
    "Unassigned_NoFeatures",
    "Unassigned_Overlapping_Length",
    "Unassigned_Ambiguity",
];

fn invalid(message: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message.into())
}

/// Options for the supported alignment-counting mode (no `--countReadPairs`).
#[derive(Clone, Debug)]
pub struct CountOptions {
    /// `-p`: declare paired input. R2 strand inversion also applies in SE mode.
    pub paired_end: bool,
    /// `-B`: in this mode rejects unpaired records in a paired input library.
    /// It does NOT reject paired records with an unmapped mate.
    pub both_ends: bool,
    /// `-C`: retained for command reconstruction; inactive without countReadPairs.
    pub no_chimeric: bool,
    /// `-s`: 0 = unstranded, 1 = forward, 2 = reverse.
    pub strand: u8,
    /// BAM inflate threads; assignment uses a deterministic serial stream.
    pub threads: usize,
}

impl Default for CountOptions {
    fn default() -> Self {
        Self {
            paired_end: true,
            both_ends: true,
            no_chimeric: true,
            strand: 2,
            threads: 6,
        }
    }
}

#[derive(Clone, Debug)]
struct Exon {
    start: u32,
    end: u32,   // one-based, inclusive
    strand: u8, // 0 = +, 1 = -, 255 = . (C unsigned char)
    group: usize,
}

#[derive(Debug)]
struct Chromosome {
    name: String,
    exons: Vec<Exon>,
    prefix_end: Vec<u32>,
}

struct IntervalPartition<'a> {
    chromosome: &'a str,
    strand: u8,
    intervals: Vec<(u32, u32)>,
}

/// Immutable annotation, reusable across samples and grouping by any GTF attribute.
#[derive(Debug)]
pub struct Annotation {
    groups: Vec<String>,
    chromosomes: Vec<Chromosome>,
    by_name: HashMap<String, usize>,
}

impl Annotation {
    pub fn read(path: &Path, attribute: &str, feature_type: &str) -> io::Result<Self> {
        Self::from_gtf(BufReader::new(File::open(path)?), attribute, feature_type)
    }

    /// Read uncompressed GTF/GFF attributes. Features and groups retain Subread's
    /// encounter order; chromosomes are sorted internally by feature start.
    pub fn from_gtf(reader: impl BufRead, attribute: &str, feature_type: &str) -> io::Result<Self> {
        let mut annotation = Self {
            groups: Vec::new(),
            chromosomes: Vec::new(),
            by_name: HashMap::new(),
        };
        let mut groups = HashMap::new();
        for (line_no, line) in reader.lines().enumerate() {
            let line = line?;
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let fields: Vec<_> = line.split('\t').collect();
            if fields.len() < 3 {
                return Err(invalid(format!(
                    "GTF line {}: fewer than 3 fields",
                    line_no + 1
                )));
            }
            if !feature_type.split(',').any(|t| t == fields[2]) {
                continue;
            }
            if fields.len() < 9 {
                return Err(invalid(format!(
                    "GTF line {}: fewer than 9 fields",
                    line_no + 1
                )));
            }
            let name = attribute_value(fields[8], attribute).ok_or_else(|| {
                invalid(format!(
                    "GTF line {}: missing or malformed {attribute}",
                    line_no + 1
                ))
            })?;
            let start = coordinate(fields[3])?;
            let end = coordinate(fields[4])?;
            if start == 0 || end < start {
                return Err(invalid(format!(
                    "GTF line {}: invalid interval",
                    line_no + 1
                )));
            }
            let group = *groups.entry(name.clone()).or_insert_with(|| {
                let id = annotation.groups.len();
                annotation.groups.push(name);
                id
            });
            let next_chr = annotation.chromosomes.len();
            let chr = *annotation
                .by_name
                .entry(fields[0].to_owned())
                .or_insert(next_chr);
            if chr == next_chr {
                annotation.chromosomes.push(Chromosome {
                    name: fields[0].to_owned(),
                    exons: Vec::new(),
                    prefix_end: Vec::new(),
                });
            }
            annotation.chromosomes[chr].exons.push(Exon {
                start,
                end,
                group,
                strand: match fields[6] {
                    "+" => 0,
                    "-" => 1,
                    _ => 255,
                },
            });
        }
        if annotation.groups.is_empty() {
            return Err(invalid("no features loaded from GTF"));
        }
        for chr in &mut annotation.chromosomes {
            feature_order(&mut chr.exons);
            let mut max_end = 0;
            chr.prefix_end = chr
                .exons
                .iter()
                .map(|e| {
                    max_end = max_end.max(e.end);
                    max_end
                })
                .collect();
        }
        Ok(annotation)
    }

    pub fn group_names(&self) -> &[String] {
        &self.groups
    }

    fn chromosome(&self, name: &str) -> Option<usize> {
        self.by_name.get(name).copied().or_else(|| {
            if let Some(short) = name.strip_prefix("chr") {
                self.by_name.get(short).copied()
            } else if name.len() <= 2 {
                self.by_name.get(&format!("chr{name}")).copied()
            } else {
                None
            }
        })
    }

    fn hits(&self, chr: usize, begin: u32, end: u32, strand: u8, mode: u8, hits: &mut Vec<usize>) {
        let chr = &self.chromosomes[chr];
        let lo = chr.prefix_end.partition_point(|&e| e < begin);
        let hi = chr.exons.partition_point(|e| e.start <= end);
        for exon in &chr.exons[lo..hi] {
            if exon.end >= begin
                && match mode {
                    0 => true,
                    1 => exon.strand == strand,
                    _ => exon.strand != strand,
                }
            {
                hits.push(exon.group);
            }
        }
    }

    fn metadata(&self) -> Vec<[String; 5]> {
        let mut grouped = vec![Vec::new(); self.groups.len()];
        for chr in &self.chromosomes {
            for exon in &chr.exons {
                grouped[exon.group].push((chr.name.as_str(), exon));
            }
        }
        grouped
            .into_iter()
            .map(|exons| {
                let mut partitions: Vec<IntervalPartition<'_>> = Vec::new();
                for (chr, exon) in exons {
                    if let Some(partition) = partitions
                        .iter_mut()
                        .find(|p| p.chromosome == chr && p.strand == exon.strand)
                    {
                        partition.intervals.push((exon.start, exon.end));
                    } else {
                        partitions.push(IntervalPartition {
                            chromosome: chr,
                            strand: exon.strand,
                            intervals: vec![(exon.start, exon.end)],
                        });
                    }
                }
                let mut columns: [Vec<String>; 4] = std::array::from_fn(|_| Vec::new());
                // The reference writes a signed 32-bit union length, including
                // wraparound for very large meta-features across chromosomes.
                let mut length = 0u32;
                for IntervalPartition {
                    chromosome: chr,
                    strand,
                    intervals,
                } in partitions
                {
                    // Already start-sorted within each chromosome/strand. Output lists
                    // original intervals (including duplicates), but length uses union.
                    let mut union_end = 0;
                    for (start, end) in intervals {
                        columns[0].push(chr.to_owned());
                        columns[1].push(start.to_string());
                        columns[2].push(end.to_string());
                        columns[3].push(
                            match strand {
                                0 => "+",
                                1 => "-",
                                _ => ".",
                            }
                            .to_owned(),
                        );
                        length = length.wrapping_add(end.saturating_sub(union_end.max(start - 1)));
                        union_end = union_end.max(end);
                    }
                }
                [
                    columns[0].join(";"),
                    columns[1].join(";"),
                    columns[2].join(";"),
                    columns[3].join(";"),
                    (length as i32).to_string(),
                ]
            })
            .collect()
    }
}

fn coordinate(value: &str) -> io::Result<u32> {
    // atoll/atoi consume a decimal prefix in the C loader.
    let digits: String = value.chars().take_while(char::is_ascii_digit).collect();
    let n = digits
        .parse::<u32>()
        .map_err(|_| invalid("invalid GTF coordinate"))?;
    if value.len() > 10 || n > i32::MAX as u32 {
        return Err(invalid("GTF coordinate exceeds 2^31-1"));
    }
    Ok(n)
}

/// Parse attributes independently of semicolons inside quoted values. Subread
/// accepts GTF spaces and GFF '=', trims spaces, and uses the last repeated key.
fn attribute_value(attrs: &str, key: &str) -> Option<String> {
    let mut rest = attrs.trim_end_matches(['\r', '\n']);
    let mut found = None;
    while !rest.is_empty() {
        rest = rest.trim_start_matches(' ');
        if rest.is_empty() {
            break;
        }
        let n = rest.find([' ', '=', '"'])?;
        let name = &rest[..n];
        if name.is_empty()
            || name.len() > 190
            || !name.bytes().all(|c| c.is_ascii_alphanumeric() || c == b'_')
        {
            return None;
        }
        rest = rest[n..].trim_start_matches([' ', '=']);
        let value;
        if let Some(quoted) = rest.strip_prefix('"') {
            let end = quoted.find('"')?;
            value = &quoted[..end];
            rest = quoted[end + 1..].trim_start_matches(' ');
            if !rest.is_empty() && !rest.starts_with(';') {
                return None;
            }
        } else {
            let end = rest.find(';').unwrap_or(rest.len());
            value = &rest[..end];
            rest = &rest[end..];
        }
        if name == key {
            let bytes = value.trim_start_matches(' ').as_bytes();
            let truncated = &bytes[..bytes.len().min(254)];
            found = Some(
                String::from_utf8_lossy(truncated)
                    .trim_end_matches(' ')
                    .to_owned(),
            );
        }
        rest = rest.strip_prefix(';').unwrap_or(rest);
    }
    found.filter(|v| !v.is_empty())
}

// Annotation output exposes Subread's tie ordering: selection-sort leaves of
// <=11 entries, recursively merged with the right half winning equal starts.
fn feature_order(exons: &mut [Exon]) {
    if exons.len() <= 11 {
        for i in 0..exons.len() {
            let mut least = i;
            for j in i + 1..exons.len() {
                if exons[j].start < exons[least].start {
                    least = j;
                }
            }
            exons.swap(i, least);
        }
    } else {
        let half = exons.len() / 2;
        feature_order(&mut exons[..half]);
        feature_order(&mut exons[half..]);
        let mut merged = Vec::with_capacity(exons.len());
        let (mut left, mut right) = (0, half);
        while left < half || right < exons.len() {
            let at = if left == half
                || (right < exons.len() && exons[right].start <= exons[left].start)
            {
                let at = right;
                right += 1;
                at
            } else {
                let at = left;
                left += 1;
                at
            };
            merged.push(exons[at].clone());
        }
        exons.clone_from_slice(&merged);
    }
}

/// Per-input integer counts. All 14 summary rows are always emitted.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Counts {
    pub groups: Vec<u64>,
    pub summary: [u64; 14],
}

/// Stream a BAM with nomics-io; memory is bounded by the annotation and a few BAM
/// batches, with no name pairing or input-order assumptions in this mode.
pub fn count_bam(
    annotation: &Annotation,
    path: &Path,
    options: &CountOptions,
) -> io::Result<Counts> {
    if options.strand > 2 || options.threads == 0 {
        return Err(invalid("strand must be 0..2 and threads must be positive"));
    }
    let mut stream = BamStream::open(path, 4096, 2, options.threads)?;
    let chrs: Vec<_> = stream
        .header
        .refs
        .iter()
        .map(|r| {
            std::str::from_utf8(&r.name)
                .ok()
                .and_then(|n| annotation.chromosome(n))
        })
        .collect();
    let mut counts = Counts {
        groups: vec![0; annotation.groups.len()],
        summary: [0; 14],
    };
    let mut hits = Vec::new();
    let mut saw_pair = false;
    while let Some(batch) = stream.next_batch()? {
        for i in 0..batch.len() {
            let record = RecordView::new(batch.record(i))?;
            saw_pair |= record.flag() & 1 != 0;
            let status = assign(annotation, record, &chrs, options, &mut hits)?;
            counts.summary[status] += 1;
            if status == 0 {
                counts.groups[hits[0]] += 1;
            }
        }
    }
    if options.paired_end && !saw_pair {
        return Err(invalid(
            "no paired-end reads detected in paired input library",
        ));
    }
    Ok(counts)
}

fn assign(
    annotation: &Annotation,
    record: RecordView<'_>,
    chrs: &[Option<usize>],
    options: &CountOptions,
    hits: &mut Vec<usize>,
) -> io::Result<usize> {
    let flag = record.flag();
    if !options.paired_end && flag & 1 != 0 {
        return Err(invalid("paired-end read in single-end input library"));
    }
    if flag & 4 != 0 {
        return Ok(1);
    }
    if options.paired_end && options.both_ends && flag & 1 == 0 {
        return Ok(3);
    }
    if nh(record)? > 1 {
        return Ok(8);
    }
    hits.clear();
    let Some(chr) = usize::try_from(record.ref_id())
        .ok()
        .and_then(|r| chrs.get(r))
        .copied()
        .flatten()
    else {
        return Ok(11);
    };
    if record.pos() < 0 {
        return Err(invalid("mapped BAM record has negative position"));
    }
    let strand = u8::from(flag & 16 != 0) ^ u8::from(flag & 128 != 0);
    let mut cursor = record.pos() as u32 + 1;
    let mut start = cursor;
    let mut length = 0u16;
    let mut sections = 0;
    for c in record.cigar() {
        let len = c >> 4;
        match c & 15 {
            0 | 7 | 8 => {
                cursor = cursor.wrapping_add(len);
                length = length.wrapping_add(len as u16);
            }
            1..=3 => {
                if sections < 65536 {
                    // The C parser emits even empty sections at I/D/N. A zero
                    // section gives [start, start-1]; its search can hit a
                    // feature spanning that boundary. Preserve this quirk.
                    annotation.hits(
                        chr,
                        start,
                        start.wrapping_add(u32::from(length)).wrapping_sub(1),
                        strand,
                        options.strand,
                        hits,
                    );
                    sections += 1;
                }
                if matches!(c & 15, 2 | 3) {
                    cursor = cursor.wrapping_add(len);
                }
                start = cursor;
                length = 0;
            }
            _ => {} // S/H/P consume no reference and contribute no overlap.
        }
    }
    if length > 0 && sections < 65536 {
        annotation.hits(
            chr,
            start,
            start.wrapping_add(u32::from(length)).wrapping_sub(1),
            strand,
            options.strand,
            hits,
        );
    }
    hits.sort_unstable();
    hits.dedup();
    Ok(match hits.len() {
        0 => 11,
        1 => 0,
        _ => 13,
    })
}

fn nh(record: RecordView<'_>) -> io::Result<i64> {
    let mut aux = record
        .aux()
        .ok_or_else(|| invalid("truncated BAM auxiliary fields"))?;
    while !aux.is_empty() {
        if aux.len() < 3 {
            return Err(invalid("truncated BAM tag"));
        }
        let ty = aux[2];
        let body = &aux[3..];
        let len = aux_value_len(ty, body)
            .filter(|&n| n <= body.len())
            .ok_or_else(|| invalid("malformed BAM tag"))?;
        if &aux[..2] == b"NH" {
            return Ok(match ty {
                // Subread zeroes an int then memcpy's 1/2/4 bytes into it:
                // c/s are zero-extended, and I has signed 32-bit semantics.
                b'c' | b'C' => i64::from(body[0]),
                b's' | b'S' => i64::from(u16::from_le_bytes(body[..2].try_into().unwrap())),
                b'i' | b'I' => i64::from(i32::from_le_bytes(body[..4].try_into().unwrap())),
                _ => 1,
            });
        }
        aux = &body[len..];
    }
    Ok(1)
}

/// Emit the exact featureCounts comment format, including one space after every
/// quoted argument. Arguments are not shell-escaped (neither are Subread's).
pub fn write_table(
    mut out: impl Write,
    annotation: &Annotation,
    inputs: &[String],
    counts: &[Counts],
    command: &[String],
) -> io::Result<()> {
    validate_columns(annotation, inputs, counts)?;
    write!(out, "# Program:featureCounts v{VERSION}; Command:")?;
    for arg in command {
        write!(out, "\"{arg}\" ")?;
    }
    writeln!(out)?;
    writeln!(
        out,
        "Geneid\tChr\tStart\tEnd\tStrand\tLength\t{}",
        inputs.join("\t")
    )?;
    for (i, fields) in annotation.metadata().iter().enumerate() {
        write!(out, "{}\t{}", annotation.groups[i], fields.join("\t"))?;
        for count in counts {
            write!(out, "\t{}", count.groups[i])?;
        }
        writeln!(out)?;
    }
    Ok(())
}

fn validate_columns(
    annotation: &Annotation,
    inputs: &[String],
    counts: &[Counts],
) -> io::Result<()> {
    if inputs.is_empty()
        || inputs.len() != counts.len()
        || counts
            .iter()
            .any(|c| c.groups.len() != annotation.groups.len())
    {
        return Err(invalid(
            "count table dimensions do not match annotation/inputs",
        ));
    }
    Ok(())
}

pub fn write_summary(mut out: impl Write, inputs: &[String], counts: &[Counts]) -> io::Result<()> {
    if inputs.is_empty() || inputs.len() != counts.len() {
        return Err(invalid("summary dimensions do not match inputs"));
    }
    writeln!(out, "Status\t{}", inputs.join("\t"))?;
    for (i, status) in STATUSES.iter().enumerate() {
        write!(out, "{status}")?;
        for count in counts {
            write!(out, "\t{}", count.summary[i])?;
        }
        writeln!(out)?;
    }
    Ok(())
}

/// File-oriented API. Paths are rendered exactly as supplied, without resolving
/// symlinks or making paths absolute. `command` includes argv[0]; callers may use
/// the golden's staging names to reproduce its provenance comment verbatim.
#[derive(Debug)]
pub struct Options {
    pub annotation: PathBuf,
    pub inputs: Vec<PathBuf>,
    pub output: PathBuf,
    pub attribute: String,
    pub feature_type: String,
    pub counting: CountOptions,
    pub command: Vec<String>,
}

pub fn run(options: &Options) -> io::Result<Vec<Counts>> {
    let annotation = Annotation::read(
        &options.annotation,
        &options.attribute,
        &options.feature_type,
    )?;
    if options.inputs.is_empty() {
        return Err(invalid("at least one input BAM is required"));
    }
    let counts = options
        .inputs
        .iter()
        .map(|p| count_bam(&annotation, p, &options.counting))
        .collect::<io::Result<Vec<_>>>()?;
    let inputs: Vec<_> = options
        .inputs
        .iter()
        .map(|p| p.to_string_lossy().into_owned())
        .collect();
    let mut table = BufWriter::new(File::create(&options.output)?);
    write_table(&mut table, &annotation, &inputs, &counts, &options.command)?;
    table.flush()?;
    let mut summary_path = options.output.as_os_str().to_os_string();
    summary_path.push(".summary");
    let mut summary = BufWriter::new(File::create(summary_path)?);
    write_summary(&mut summary, &inputs, &counts)?;
    summary.flush()?;
    Ok(counts)
}

/// nf-core's simple biotype MultiQC conversion, for a single sample. The header
/// is the pipeline's MIT-licensed assets/biotypes_header.txt. No rrna table is
/// written for zero total counts, just as in mqc_features_stat.py.
/// Upstream limits the biotype plot to 100 rows by default.
pub fn write_multiqc(
    annotation: &Annotation,
    counts: &Counts,
    sample: &str,
    prefix: &Path,
) -> io::Result<()> {
    if annotation.groups.len() > 100 || counts.groups.len() != annotation.groups.len() {
        return Err(invalid(
            "MultiQC requires <=100 biotype groups and matching counts",
        ));
    }
    let mut path = prefix.as_os_str().to_os_string();
    path.push(".biotype_counts_mqc.tsv");
    let mut out = BufWriter::new(File::create(path)?);
    out.write_all(include_bytes!("biotypes_header.txt"))?;
    for (name, count) in annotation.groups.iter().zip(&counts.groups) {
        writeln!(out, "{name}\t{count}")?;
    }
    out.flush()?;
    let total: u64 = counts.groups.iter().sum();
    if total == 0 {
        return Ok(());
    }
    let percent = annotation
        .groups
        .iter()
        .position(|g| g == "rRNA")
        .map_or("0".to_owned(), |i| {
            let n = (counts.groups[i] as f64 / total as f64) * 100.0;
            python_percentage(n)
        });
    let mut path = prefix.as_os_str().to_os_string();
    path.push(".biotype_counts_rrna_mqc.tsv");
    std::fs::write(
        path,
        format!(
            "#id: 'biotype-gs'\n#plot_type: 'generalstats'\n#pconfig:\n#    percent_rRNA:\n#        title: '% rRNA'\n#        namespace: 'Biotype Counts'\n#        description: '% reads overlapping rRNA features'\n#        max: 100\n#        min: 0\n#        scale: 'RdYlGn-rev'\nSample\tpercent_rRNA\n'{sample}'\t{percent}\n"
        ),
    )
}

// Python str(float) uses scientific notation below 1e-4 and pads the exponent
// to two digits. The percentage range is 0..100, so no upper exponent case.
fn python_percentage(n: f64) -> String {
    if n != 0.0 && n < 1e-4 {
        let scientific = format!("{n:e}");
        let (mantissa, exponent) = scientific.split_once('e').unwrap();
        let exponent: i32 = exponent.parse().unwrap();
        format!("{mantissa}e{exponent:+03}")
    } else if n.fract() == 0.0 {
        format!("{n:.1}")
    } else {
        n.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nomics_io::bam::{BinRule, IntType, Record, Value};

    fn annotation(gtf: &str) -> Annotation {
        Annotation::from_gtf(gtf.as_bytes(), "gene_id", "exon").unwrap()
    }
    fn record(flag: u16, pos: i32, cigar: Vec<u32>, nh: Option<(IntType, i64)>) -> Vec<u8> {
        let record = Record {
            ref_id: 0,
            pos,
            mapq: 0,
            flag,
            next_ref_id: 0,
            next_pos: 99,
            tlen: 100,
            read_name: b"read".to_vec(),
            cigar,
            seq: vec![b'A'; 10],
            qual: vec![30; 10],
            tags: nh
                .into_iter()
                .map(|(ty, n)| (*b"NH", Value::Int(ty, n)))
                .collect(),
        };
        let mut raw = Vec::new();
        record.encode(BinRule::Sam, &mut raw).unwrap();
        raw
    }
    fn status(a: &Annotation, raw: &[u8], strand: u8) -> usize {
        assign(
            a,
            RecordView::new(raw).unwrap(),
            &[Some(0)],
            &CountOptions {
                strand,
                ..CountOptions::default()
            },
            &mut Vec::new(),
        )
        .unwrap()
    }

    #[test]
    fn nfcore_counts_alignments_with_r2_strand_inversion() {
        let a = annotation("I\tx\texon\t100\t200\t.\t+\t.\tgene_id \"a\";\n");
        assert_eq!(status(&a, &record(1 | 16 | 64, 99, vec![160], None), 2), 0);
        assert_eq!(status(&a, &record(1 | 128, 99, vec![160], None), 2), 0);
        assert_eq!(status(&a, &record(1 | 64, 99, vec![160], None), 2), 11);
        assert_eq!(
            status(&a, &record(1 | 16 | 128, 99, vec![160], None), 2),
            11
        );
        assert_eq!(status(&a, &record(1 | 64, 99, vec![160], None), 1), 0);
    }

    #[test]
    fn filter_precedence_and_defaults_are_source_exact() {
        let a = annotation("I\tx\texon\t100\t200\t.\t+\t.\tgene_id \"a\";\n");
        // Unmapped precedes NH; unpaired with -B is singleton before NH.
        assert_eq!(
            status(&a, &record(1 | 4, 99, vec![160], Some((IntType::U8, 2))), 0),
            1
        );
        assert_eq!(
            status(&a, &record(0, 99, vec![160], Some((IntType::U8, 2))), 0),
            3
        );
        assert_eq!(
            status(&a, &record(1, 99, vec![160], Some((IntType::U8, 2))), 0),
            8
        );
        // -B/-C do not reject a paired singleton, same-strand mates or other
        // chromosomes in alignment mode. Zero MAPQ/duplicate/secondary/supplementary count.
        assert_eq!(
            status(
                &a,
                &record(1 | 8 | 256 | 1024 | 2048, 99, vec![160], None),
                0
            ),
            0
        );
        assert_eq!(
            status(&a, &record(1, 99, vec![160], Some((IntType::U8, 1))), 0),
            0
        );
    }

    #[test]
    fn union_meta_features_ambiguity_and_splice_gaps() {
        let a = annotation(
            "I\tx\texon\t100\t120\t.\t+\t.\tgene_id \"a\";\nI\tx\texon\t110\t140\t.\t+\t.\tgene_id \"a\";\nI\tx\texon\t130\t150\t.\t+\t.\tgene_id \"b\";\n",
        );
        assert_eq!(a.metadata()[0][4], "41");
        assert_eq!(a.metadata()[0][1], "100;110");
        assert_eq!(status(&a, &record(1, 109, vec![160], None), 0), 0);
        assert_eq!(status(&a, &record(1, 129, vec![160], None), 0), 13);
        let b = annotation("I\tx\texon\t120\t130\t.\t+\t.\tgene_id \"intron\";\n");
        for op in [2, 3] {
            // deletion/intron bases don't overlap
            assert_eq!(
                status(&b, &record(1, 99, vec![160, 50 << 4 | op, 160], None), 0),
                11
            );
        }
        assert_eq!(
            status(
                &b,
                &record(1, 119, vec![5 << 4 | 4, 5 << 4 | 7, 5 << 4 | 8], None),
                0
            ),
            0
        );
    }

    #[test]
    fn attributes_and_unknown_strands() {
        let a = annotation(
            "# comment\nI\tx\texon\t1\t10\t.\t.\t.\tother \"a;b\"; gene_id=first; gene_id \" last \";\n",
        );
        assert_eq!(a.group_names(), &["last"]);
        assert_eq!(status(&a, &record(1, 0, vec![160], None), 1), 11);
        assert_eq!(status(&a, &record(1, 0, vec![160], None), 2), 0);
        assert_eq!(a.chromosome("chrI"), Some(0));
        assert!(
            Annotation::from_gtf(
                b"I\tx\texon\t0\t10\t.\t+\t.\tgene_id \"a\";".as_slice(),
                "gene_id",
                "exon"
            )
            .is_err()
        );
        assert!(attribute_value("gene_id \"a\"; broken-field \"bad\";", "gene_id").is_none());
    }

    #[test]
    fn command_format_keeps_trailing_space_and_all_summary_rows() {
        let a = annotation("I\tx\texon\t1\t10\t.\t+\t.\tgene_id \"a\";\n");
        let counts = [Counts {
            groups: vec![0],
            summary: [0; 14],
        }];
        let mut out = Vec::new();
        write_table(
            &mut out,
            &a,
            &["in.bam".into()],
            &counts,
            &["featureCounts".into(), "in.bam".into()],
        )
        .unwrap();
        assert!(out.starts_with(
            b"# Program:featureCounts v2.1.1; Command:\"featureCounts\" \"in.bam\" \n"
        ));
        out.clear();
        write_summary(&mut out, &["in.bam".into()], &counts).unwrap();
        assert_eq!(out.iter().filter(|&&b| b == b'\n').count(), 15);
    }

    #[test]
    fn source_integer_tag_and_python_percentage_quirks() {
        for (kind, value, expected) in [
            (IntType::I8, -1, 255),
            (IntType::I16, -1, 65535),
            (IntType::I32, -1, -1),
            (IntType::U32, u32::MAX as i64, -1),
        ] {
            let raw = record(1, 99, vec![160], Some((kind, value)));
            assert_eq!(nh(RecordView::new(&raw).unwrap()).unwrap(), expected);
        }
        assert_eq!(python_percentage(0.0), "0.0");
        assert_eq!(python_percentage(100.0), "100.0");
        assert_eq!(python_percentage(0.00001), "1e-05");
        assert_eq!(python_percentage(0.000000123), "1.23e-07");
        assert_eq!(python_percentage(0.0001), "0.0001");
    }
}
