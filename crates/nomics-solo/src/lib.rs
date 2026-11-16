//! Fresh implementation of STAR 2.7.10b's default CB_UMI_Simple/Gene counting.
//! Alignment evidence is supplied by an oracle BAM; corrected CB/UB are never read.
mod barcode;
mod output;
mod umi;

use barcode::{Match, Whitelist, encode};
use nomics_io::bam::{RecordView, aux_value_len};
use nomics_io::bamstream::BamStream;
use nomics_mem::HostArray;
use std::collections::HashMap;
use std::fs::{self, File};
use std::hash::{DefaultHasher, Hash, Hasher};
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

pub struct Config {
    pub bam: PathBuf,
    pub whitelist: PathBuf,
    /// STAR features.tsv in index order; GTF derivation is deliberately deferred.
    pub features: PathBuf,
    /// Directory that will contain Barcodes.stats and Gene/ (usually Solo.out).
    pub output: PathBuf,
    pub threads: usize,
}

#[derive(Debug)]
pub struct Report {
    pub reads: u64,
    pub counted_reads: u64,
    pub umis: u64,
    pub cells: usize,
    pub elapsed_seconds: f64,
}

pub(crate) fn invalid(message: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message.into())
}
pub(crate) fn open_text(path: &Path) -> io::Result<Box<dyn BufRead>> {
    if path.extension().is_some_and(|e| e == "gz") {
        Ok(Box::new(BufReader::new(nomics_io::gzip::GzReader::open(
            path,
        )?)))
    } else {
        Ok(Box::new(BufReader::new(File::open(path)?)))
    }
}

#[derive(Default)]
struct Stats {
    bar: [u64; 12],
    feat: [u64; 11],
    reads: u64,
    cb_bases: u64,
    cb_q30: u64,
    rna_bases: u64,
    rna_q30: u64,
    mapped: u64,
    unique: u64,
}

// A read with potentially several genomic alignments. Retain at most two distinct
// gene IDs: Unique needs one gene, and the statistics only distinguish 0/1/>1.
struct MultiRead {
    genes: Vec<u32>,
    raw: [u8; 44], // CR(16), CY(16), UR(12)
    seen: u32,
    expected: u32,
    primary: bool,
}

fn integer_tag(view: &RecordView<'_>, name: &[u8; 2]) -> io::Result<Option<u32>> {
    let mut aux = view
        .aux()
        .ok_or_else(|| invalid("missing auxiliary block"))?;
    while aux.len() >= 3 {
        let body = &aux[3..];
        let len = aux_value_len(aux[2], body).ok_or_else(|| invalid("malformed BAM tag"))?;
        if &aux[..2] == name {
            return match aux[2] {
                b'C' => Ok(Some(u32::from(body[0]))),
                b'S' => Ok(Some(u32::from(u16::from_le_bytes(
                    body[..2].try_into().unwrap(),
                )))),
                b'I' | b'i' => Ok(Some(u32::from_le_bytes(body[..4].try_into().unwrap()))),
                b'c' if body[0] < 128 => Ok(Some(u32::from(body[0]))),
                b's' => {
                    let v = i16::from_le_bytes(body[..2].try_into().unwrap());
                    u32::try_from(v)
                        .map(Some)
                        .map_err(|_| invalid("negative NH"))
                }
                _ => Err(invalid("invalid integer tag")),
            };
        }
        aux = &body[len..];
    }
    Ok(None)
}

fn required<'a>(view: &RecordView<'a>, tag: &[u8; 2], len: usize) -> io::Result<&'a [u8]> {
    let value = view
        .string_tag(tag)
        .ok_or_else(|| invalid(format!("missing {} tag", String::from_utf8_lossy(tag))))?;
    if value.len() != len {
        return Err(invalid(format!(
            "{} length {}, expected {len}",
            String::from_utf8_lossy(tag),
            value.len()
        )));
    }
    Ok(value)
}
fn gene_evidence(view: &RecordView<'_>, features: &HashMap<Vec<u8>, u32>) -> io::Result<Vec<u32>> {
    let gx = view.string_tag(b"gx").ok_or_else(|| {
        invalid(
            "mapped oracle records require lowercase gx (uppercase GX loses multi-gene evidence)",
        )
    })?;
    if gx == b"-" {
        return Ok(Vec::new());
    }
    let mut genes = Vec::new();
    for g in gx.split(|b| *b == b';') {
        let id = *features
            .get(g)
            .ok_or_else(|| invalid(format!("unknown gene {}", String::from_utf8_lossy(g))))?;
        if !genes.contains(&id) && genes.len() < 2 {
            genes.push(id);
        }
    }
    Ok(genes)
}

fn record_gene(
    genes: &[u32],
    multi: bool,
    raw: &[u8; 44],
    m: &Match,
    stats: &mut Stats,
    spool: &mut BufWriter<File>,
) -> io::Result<()> {
    if m.code < 0 {
        return Ok(());
    }
    match genes.len() {
        0 => stats.feat[1] += 1,
        1 => {
            // Spool the raw tags, not corrected barcode assignments. All priors
            // must be complete before any non-exact barcode can be resolved.
            spool.write_all(&genes[0].to_le_bytes())?;
            spool.write_all(raw)?;
        }
        _ => {
            stats.feat[2] += 1;
            if multi {
                stats.feat[3] += 1;
            }
        }
    }
    Ok(())
}

struct Scratch(PathBuf);
impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

/// Count an oracle BAM. Output must be absent, preventing accidental replacement.
pub fn count(config: &Config) -> io::Result<Report> {
    let start = Instant::now();
    if config.output.exists() {
        return Err(io::Error::new(
            io::ErrorKind::AlreadyExists,
            "output directory already exists",
        ));
    }
    fs::create_dir_all(config.output.parent().unwrap_or(Path::new(".")))?;
    let scratch_path = config
        .output
        .with_extension(format!("solo-tmp-{}", std::process::id()));
    fs::create_dir(&scratch_path)?;
    let scratch = Scratch(scratch_path);
    let mut features_text = String::new();
    open_text(&config.features)?.read_to_string(&mut features_text)?;
    let mut features = HashMap::new();
    for (i, line) in features_text.lines().enumerate() {
        let fields: Vec<_> = line.split('\t').collect();
        if fields.len() != 3
            || fields.iter().any(|f| f.is_empty())
            || fields[2] != "Gene Expression"
        {
            return Err(invalid(
                "features.tsv must have gene ID, name, and Gene Expression columns",
            ));
        }
        if features
            .insert(fields[0].as_bytes().to_vec(), i as u32)
            .is_some()
        {
            return Err(invalid("duplicate feature ID"));
        }
    }
    if features.is_empty() || features.len() > (1 << 17) || !features_text.ends_with('\n') {
        return Err(invalid(
            "features require 1..=2^17 rows and a final newline",
        ));
    }
    let mut wl = Whitelist::load(&config.whitelist)?;
    let mut stats = Stats::default();
    let spool_path = scratch.0.join("gene-reads.bin");
    let mut spool = BufWriter::with_capacity(1 << 20, File::create(&spool_path)?);
    const PARTITIONS: usize = 64;
    let mut partitions = (0..PARTITIONS)
        .map(|i| {
            File::create(scratch.0.join(format!("multi-{i}")))
                .map(|f| BufWriter::with_capacity(1 << 16, f))
        })
        .collect::<io::Result<Vec<_>>>()?;
    let mut stream = BamStream::open(&config.bam, 8192, 2, config.threads.max(1))?;
    let mut next_progress = 5_000_000;
    while let Some(batch) = stream.next_batch()? {
        for i in 0..batch.len() {
            let v = RecordView::new(batch.record(i))?;
            if v.flag() & (0x800 | 0x1) != 0 {
                return Err(invalid(
                    "expected single-end STAR oracle without supplementary records",
                ));
            }
            let primary = v.flag() & 0x100 == 0;
            let mapped = v.flag() & 4 == 0;
            let nh = integer_tag(&v, b"NH")?.unwrap_or(0);
            if mapped && nh == 0 {
                return Err(invalid("mapped record missing positive NH"));
            }
            let cb = required(&v, b"CR", 16)?;
            let cq = required(&v, b"CY", 16)?;
            let ur = required(&v, b"UR", 12)?;
            let uq = required(&v, b"UY", 12)?;
            if cq.iter().chain(uq).any(|&q| !(33..=126).contains(&q)) {
                return Err(invalid("invalid ASCII barcode/UMI quality"));
            }
            let mut raw = [0; 44];
            raw[..16].copy_from_slice(cb);
            raw[16..32].copy_from_slice(cq);
            raw[32..].copy_from_slice(ur);
            let m = wl.match_read(cb, cq, ur);
            if primary {
                stats.reads += 1;
                stats.bar[m.stat()] += 1;
                if m.code == 0 {
                    wl.exact[m.candidates[0].0 as usize] += 1;
                }
                stats.cb_bases += 28;
                stats.cb_q30 += cq.iter().chain(uq).filter(|&&q| q >= 63).count() as u64;
                let qualities = v.quals().ok_or_else(|| invalid("RNA qualities absent"))?;
                stats.rna_bases += qualities.len() as u64;
                stats.rna_q30 += qualities.iter().filter(|&&q| q >= 30).count() as u64;
                if mapped {
                    stats.mapped += 1;
                    if nh == 1 {
                        stats.unique += 1;
                    }
                } else if m.code >= 0 {
                    stats.feat[0] += 1;
                }
            }
            if !mapped {
                continue;
            }
            let genes = gene_evidence(&v, &features)?;
            if nh == 1 {
                if !primary {
                    return Err(invalid("NH=1 secondary alignment"));
                }
                record_gene(&genes, false, &raw, &m, &mut stats, &mut spool)?;
            } else if m.code >= 0 {
                let mut hash = DefaultHasher::new();
                v.read_name().hash(&mut hash);
                let writer = &mut partitions[hash.finish() as usize % PARTITIONS];
                let len = u16::try_from(v.read_name().len())
                    .map_err(|_| invalid("read name too long"))?;
                writer.write_all(&len.to_le_bytes())?;
                writer.write_all(v.read_name())?;
                writer.write_all(&raw)?;
                writer.write_all(&nh.to_le_bytes())?;
                writer.write_all(&[u8::from(primary), genes.len() as u8])?;
                for g in genes {
                    writer.write_all(&g.to_le_bytes())?;
                }
            }
        }
        if stats.reads >= next_progress {
            eprintln!(
                "solo-count: {} reads, {:.1}s",
                stats.reads,
                start.elapsed().as_secs_f64()
            );
            next_progress += 5_000_000;
        }
    }
    drop(stream);
    for p in &mut partitions {
        p.flush()?;
    }
    drop(partitions);
    for i in 0..PARTITIONS {
        let path = scratch.0.join(format!("multi-{i}"));
        let mut reader = BufReader::new(File::open(&path)?);
        let mut reads: HashMap<Vec<u8>, MultiRead> = HashMap::new();
        loop {
            if reader.fill_buf()?.is_empty() {
                break;
            }
            let mut b = [0; 2];
            reader.read_exact(&mut b)?;
            let mut name = vec![0; u16::from_le_bytes(b) as usize];
            reader.read_exact(&mut name)?;
            let mut raw = [0; 44];
            reader.read_exact(&mut raw)?;
            let mut n = [0; 4];
            reader.read_exact(&mut n)?;
            let nh = u32::from_le_bytes(n);
            reader.read_exact(&mut b)?;
            let primary = b[0] != 0;
            let entry = reads.entry(name).or_insert_with(|| MultiRead {
                genes: vec![],
                raw,
                seen: 0,
                expected: nh,
                primary: false,
            });
            if entry.raw != raw || entry.expected != nh || (primary && entry.primary) {
                return Err(invalid(
                    "inconsistent records or duplicate primary read name",
                ));
            }
            entry.primary |= primary;
            entry.seen += 1;
            for _ in 0..b[1] {
                reader.read_exact(&mut n)?;
                let g = u32::from_le_bytes(n);
                if entry.genes.len() < 2 && !entry.genes.contains(&g) {
                    entry.genes.push(g);
                }
            }
        }
        for r in reads.values() {
            if !r.primary || r.seen != r.expected {
                return Err(invalid(
                    "incomplete multimapper alignment group (NH mismatch)",
                ));
            }
            let m = wl.match_read(&r.raw[..16], &r.raw[16..32], &r.raw[32..]);
            record_gene(&r.genes, true, &r.raw, &m, &mut stats, &mut spool)?;
        }
        fs::remove_file(path)?;
    }
    spool.flush()?;
    drop(spool);
    eprintln!(
        "solo-count: gene evidence collected, {:.1}s",
        start.elapsed().as_secs_f64()
    );
    let capacity = File::open(&spool_path)?.metadata()?.len() as usize / 48;
    let mut accepted = HostArray::<u64>::zeroed(capacity)?;
    let mut accepted_len = 0;
    let mut read_counts = HostArray::<u32>::zeroed(wl.len)?;
    let mut tentative = HostArray::<u64>::zeroed(wl.len.div_ceil(64))?;
    let mut reader = BufReader::with_capacity(1 << 20, File::open(&spool_path)?);
    let mut record = [0; 48];
    for _ in 0..capacity {
        reader.read_exact(&mut record)?;
        let gene = u32::from_le_bytes(record[..4].try_into().unwrap());
        let m = wl.match_cb(&record[4..20], &record[20..36]);
        for &(cb, _) in &m.candidates {
            tentative[cb as usize / 64] |= 1u64 << (cb % 64);
        }
        if m.code == 0 {
            stats.feat[7] += 1;
        }
        match m.resolve(&wl) {
            Ok(cb) => {
                let umi = encode(&record[36..]).0;
                accepted[accepted_len] =
                    (u64::from(cb) << 41) | (u64::from(gene) << 24) | u64::from(umi);
                accepted_len += 1;
                read_counts[cb as usize] += 1;
            }
            Err(stat) => stats.feat[stat] += 1,
        }
    }
    drop(reader);
    fs::remove_file(spool_path)?;
    stats.feat[6] = accepted_len as u64;
    stats.feat[8] = accepted_len as u64;
    accepted[..accepted_len].sort_unstable();
    eprintln!(
        "solo-count: {} counted reads sorted, {:.1}s",
        accepted_len,
        start.elapsed().as_secs_f64()
    );
    let matrix_path = scratch.0.join("matrix.bin");
    let mut matrix = BufWriter::new(File::create(&matrix_path)?);
    let mut umi_counts = HostArray::<u32>::zeroed(wl.len)?;
    let mut gene_counts = HostArray::<u32>::zeroed(wl.len)?;
    let mut at = 0;
    let mut entries = 0u64;
    let mut unique_umis = Vec::new();
    while at < accepted_len {
        let key = accepted[at] >> 24;
        let cb = (key >> 17) as usize;
        let gene = (key & 0x1ffff) as u32;
        unique_umis.clear();
        while at < accepted_len && accepted[at] >> 24 == key {
            let u = (accepted[at] & 0xffffff) as u32;
            if unique_umis.last() != Some(&u) {
                unique_umis.push(u);
            }
            at += 1;
        }
        let count = umi::components(&unique_umis);
        matrix.write_all(&(cb as u32).to_le_bytes())?;
        matrix.write_all(&gene.to_le_bytes())?;
        matrix.write_all(&count.to_le_bytes())?;
        umi_counts[cb] += count;
        gene_counts[cb] += 1;
        entries += 1;
    }
    matrix.flush()?;
    drop(matrix);
    drop(accepted);
    stats.feat[9] = gene_counts.iter().filter(|&&x| x > 0).count() as u64;
    stats.feat[10] = umi_counts.iter().map(|&x| u64::from(x)).sum();
    // Publish only after all parsing and counting has succeeded.
    fs::create_dir(&config.output)?;
    let cells = output::write(
        config,
        &features_text,
        &wl,
        &stats,
        &matrix_path,
        entries,
        &read_counts,
        &umi_counts,
        &gene_counts,
        tentative
            .iter()
            .map(|word| word.count_ones() as usize)
            .sum(),
    )?;
    Ok(Report {
        reads: stats.reads,
        counted_reads: stats.feat[6],
        umis: stats.feat[10],
        cells,
        elapsed_seconds: start.elapsed().as_secs_f64(),
    })
}
