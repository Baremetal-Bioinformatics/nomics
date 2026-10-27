use std::collections::{BTreeMap, HashMap};
use std::io::{self, BufWriter, Write};
use std::path::Path;

use serde::Serialize;

use crate::{Options, Quantification, invalid};

const HEADER: &str = "Name\tLength\tEffectiveLength\tTPM\tNumReads";

fn attribute<'a>(text: &'a str, key: &str) -> Option<&'a str> {
    text.split(';').find_map(|item| {
        let item = item.trim();
        let rest = item.strip_prefix(key)?;
        let first = rest.as_bytes().first()?;
        if *first == b'=' {
            Some(rest[1..].trim().trim_matches('"'))
        } else if first.is_ascii_whitespace() {
            Some(rest.trim().trim_matches('"'))
        } else {
            None
        }
    })
}

const MAP_CHUNK: u64 = 4 << 20;

/// Up to `count` chunks of whole lines (each ends with '\n' unless it holds
/// the file's last line). A read error ends the batch after the lines that
/// were complete before it, as a line-at-a-time reader would have seen them.
fn read_chunks(
    reader: &mut (dyn io::BufRead + Send),
    count: usize,
    size: u64,
) -> (Vec<Vec<u8>>, Option<io::Error>) {
    let mut chunks = Vec::new();
    while chunks.len() < count {
        let mut chunk = Vec::with_capacity(size.min(MAP_CHUNK) as usize + 4096);
        let mut read = io::Read::read_to_end(&mut io::Read::take(&mut *reader, size), &mut chunk);
        if read.is_ok() && chunk.last().is_some_and(|&b| b != b'\n') {
            read = reader.read_until(b'\n', &mut chunk);
        }
        if let Err(e) = read {
            chunk.truncate(chunk.iter().rposition(|&b| b == b'\n').map_or(0, |i| i + 1));
            chunks.push(chunk);
            return (chunks, Some(e));
        }
        if chunk.is_empty() {
            break;
        }
        chunks.push(chunk);
    }
    (chunks, None)
}

/// One chunk's mappings, resolved within the chunk with the same rule as the
/// whole file (GTF: first line with both ids; TSV: last line), or the first
/// UTF-8 error. Lines keep their '\n' so error positions match line reads.
fn parse_chunk(chunk: &[u8], gtf: bool) -> io::Result<HashMap<&str, &str>> {
    let mut map = HashMap::new();
    for line in chunk.split_inclusive(|&b| b == b'\n') {
        let line = std::str::from_utf8(line)
            .map_err(|e| invalid(e.to_string()))?
            .trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        if gtf {
            if let Some(attrs) = line.split('\t').nth(8)
                && let Some(tx) = attribute(attrs, "transcript_id")
                && !map.contains_key(tx)
                && let Some(gene) = attribute(attrs, "gene_id")
            {
                map.insert(tx, gene);
            }
        } else {
            let mut columns = line.split_whitespace();
            if let (Some(tx), Some(gene)) = (columns.next(), columns.next()) {
                map.insert(tx, gene);
            }
        }
    }
    Ok(map)
}

/// Stream the annotation in line-aligned chunks: chunks parse in parallel
/// while the next batch is read, and merge in file order (GTF keeps the
/// first mapping per transcript, TSV the last), equal to a serial line scan.
/// Memory: two batches of `2 x threads` chunks of `MAP_CHUNK` bytes (~320 MB
/// on 20 cores), independent of annotation size.
fn gene_map(path: &Path) -> io::Result<HashMap<String, String>> {
    gene_map_chunked(path, MAP_CHUNK)
}

fn gene_map_chunked(path: &Path, size: u64) -> io::Result<HashMap<String, String>> {
    use rayon::prelude::*;
    let mut name = path.to_string_lossy().to_ascii_lowercase();
    for suffix in [".gz", ".bgz", ".bgzf"] {
        if name.ends_with(suffix) {
            name.truncate(name.len() - suffix.len());
            break;
        }
    }
    let gtf = [".gtf", ".gff", ".gff3"].iter().any(|s| name.ends_with(s));
    let mut map: HashMap<String, String> = HashMap::new();
    let mut reader = crate::input::reader(path)?;
    let count = 2 * rayon::current_num_threads().max(1);
    let (mut chunks, mut error) = read_chunks(&mut *reader, count, size);
    while !chunks.is_empty() || error.is_some() {
        let ((next, next_error), parsed) = if error.is_some() {
            (
                (Vec::new(), None),
                chunks.par_iter().map(|c| parse_chunk(c, gtf)).collect(),
            )
        } else {
            rayon::join(
                || read_chunks(&mut *reader, count, size),
                || {
                    chunks
                        .par_iter()
                        .map(|c| parse_chunk(c, gtf))
                        .collect::<Vec<_>>()
                },
            )
        };
        for local in parsed {
            for (tx, gene) in local? {
                // GTF: the first chunk's mapping wins; TSV: the last.
                if !gtf || !map.contains_key(tx) {
                    map.insert(tx.to_owned(), gene.to_owned());
                }
            }
        }
        if let Some(e) = error {
            return Err(e);
        }
        (chunks, error) = (next, next_error);
    }
    if map.is_empty() {
        return Err(invalid(format!(
            "no transcript-to-gene mappings in {}",
            path.display()
        )));
    }
    Ok(map)
}

fn write_genes(
    path: &Path,
    aux: &Path,
    map: &HashMap<String, String>,
    q: &Quantification,
) -> io::Result<()> {
    let mut genes: BTreeMap<&str, Vec<usize>> = BTreeMap::new();
    let mut unmatched = Vec::new();
    for (i, name) in q.names.iter().enumerate() {
        let gene = map.get(name).map(String::as_str).unwrap_or_else(|| {
            unmatched.push(name.as_str());
            name
        });
        genes.entry(gene).or_default().push(i);
    }
    let unmatched_path = aux.join("genemap_unmatched_txps.json");
    if unmatched.is_empty() {
        match std::fs::remove_file(unmatched_path) {
            Ok(()) => {}
            Err(e) if e.kind() == io::ErrorKind::NotFound => {}
            Err(e) => return Err(e),
        }
    } else {
        let json = serde_json::json!({"unmatched_transcripts": unmatched});
        std::fs::write(
            unmatched_path,
            format!("{}\n", serde_json::to_string_pretty(&json)?),
        )?;
    }
    let mut w = BufWriter::new(std::fs::File::create(path)?);
    writeln!(w, "{HEADER}")?;
    for (gene, indices) in genes {
        let tpm: f64 = indices.iter().map(|&i| q.tpm[i]).sum();
        let reads: f64 = indices.iter().map(|&i| q.counts[i]).sum();
        let (mut length, mut effective) = (0.0, 0.0);
        for &i in &indices {
            let fraction = if tpm > f64::from_bits(1) {
                q.tpm[i] / tpm
            } else {
                1.0 / indices.len() as f64
            };
            length += q.lengths[i] as f64 * fraction;
            effective += q.effective_lengths[i] * fraction;
        }
        writeln!(
            w,
            "{gene}\t{length:.3}\t{effective:.3}\t{tpm:.6}\t{reads:.3}"
        )?;
    }
    w.flush()
}

#[derive(Serialize)]
struct LibraryCounts<'a> {
    read_files: [String; 0],
    expected_format: &'a str,
    compatible_fragment_ratio: f64,
    num_compatible_fragments: u64,
    num_incompatible_fragments: u64,
    num_assigned_fragments: u64,
    num_frags_with_concordant_consistent_mappings: u64,
    num_frags_with_inconsistent_or_orphan_mappings: u64,
    strand_mapping_bias: f64,
    #[serde(flatten)]
    formats: Formats,
}

// Struct field order is part of the output contract.
#[derive(Serialize)]
#[allow(non_snake_case)]
struct Formats {
    IU: u64,
    ISF: u64,
    ISR: u64,
    OU: u64,
    OSF: u64,
    OSR: u64,
    MU: u64,
    MSF: u64,
    MSR: u64,
    U: u64,
    SF: u64,
    SR: u64,
}

pub(crate) fn write(o: &Options, q: &Quantification) -> io::Result<()> {
    // Parse before opening output files so a bad annotation cannot leave a
    // seemingly successful transcript-only run behind.
    let map = o.gene_map.as_deref().map(gene_map).transpose()?;
    let aux = o.output.join("aux_info");
    std::fs::create_dir_all(&aux)?;
    let mut w = BufWriter::new(std::fs::File::create(o.output.join("quant.sf"))?);
    writeln!(w, "{HEADER}")?;
    for i in 0..q.names.len() {
        writeln!(
            w,
            "{}\t{}\t{:.3}\t{:.6}\t{:.3}",
            q.names[i], q.lengths[i], q.effective_lengths[i], q.tpm[i], q.counts[i]
        )?;
    }
    w.flush()?;
    if let Some(map) = &map {
        write_genes(&o.output.join("quant.genes.sf"), &aux, map, q)?;
    }
    let mut w = BufWriter::new(std::fs::File::create(aux.join("ambig_info.tsv"))?);
    writeln!(w, "UniqueCount\tAmbigCount")?;
    for (unique, ambiguous) in q.unique.iter().zip(q.ambiguous.iter()) {
        writeln!(w, "{unique}\t{ambiguous}")?;
    }
    w.flush()?;
    let [iu, isf, isr, ou, osf, osr, mu, msf, msr, u, sf, sr] = q.format_counts;
    let judged = q.num_compatible + q.num_incompatible;
    let expected: crate::LibraryType = q.expected_format.parse()?;
    let (concordant, inconsistent, bias) = expected.summary(&q.format_counts);
    let counts = LibraryCounts {
        read_files: [],
        expected_format: &q.expected_format,
        compatible_fragment_ratio: if judged == 0 {
            1.0
        } else {
            q.num_compatible as f64 / judged as f64
        },
        num_compatible_fragments: if judged == 0 {
            q.num_mapped
        } else {
            q.num_compatible
        },
        num_incompatible_fragments: q.num_incompatible,
        num_assigned_fragments: q.num_mapped,
        num_frags_with_concordant_consistent_mappings: concordant,
        num_frags_with_inconsistent_or_orphan_mappings: inconsistent,
        strand_mapping_bias: bias,
        formats: Formats {
            IU: iu,
            ISF: isf,
            ISR: isr,
            OU: ou,
            OSF: osf,
            OSR: osr,
            MU: mu,
            MSF: msf,
            MSR: msr,
            U: u,
            SF: sf,
            SR: sr,
        },
    };
    std::fs::write(
        o.output.join("lib_format_counts.json"),
        serde_json::to_string_pretty(&counts)?,
    )?;
    let bytes: Vec<u8> = q
        .fragment_pmf
        .iter()
        .flat_map(|p| ((p * 10000.0).round() as i32).to_le_bytes())
        .collect();
    gzip(&aux.join("fld.gz"), &bytes)?;
    // These neutral single-element tables are upstream's own legacy-model
    // outputs when no bias flags were requested.
    gzip(&aux.join("observed_bias.gz"), &0i32.to_le_bytes())?;
    gzip(&aux.join("observed_bias_3p.gz"), &0i32.to_le_bytes())?;
    gzip(&aux.join("expected_bias.gz"), &1f64.to_le_bytes())?;
    Ok(())
}

fn gzip(path: &Path, bytes: &[u8]) -> io::Result<()> {
    use nomics_io::deflate::{Finish, ZlibFast};
    let mut compressed = Vec::new();
    let result = ZlibFast::new(3).deflate_finish(bytes, &mut compressed, bytes.len() * 2 + 1024);
    if result != Finish::Done {
        return Err(io::Error::other("auxiliary DEFLATE buffer exhausted"));
    }
    let mut f = BufWriter::new(std::fs::File::create(path)?);
    f.write_all(&[31, 139, 8, 0, 0, 0, 0, 0, 0, 255])?;
    f.write_all(&compressed)?;
    f.write_all(&nomics_io::crc32::crc32(bytes).to_le_bytes())?;
    f.write_all(&(bytes.len() as u32).to_le_bytes())?;
    f.flush()
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn attributes_require_a_key_boundary() {
        let attrs = "other_gene_id \"wrong\"; gene_id \"G\"; transcript_id=T;";
        assert_eq!(attribute(attrs, "gene_id"), Some("G"));
        assert_eq!(attribute(attrs, "transcript_id"), Some("T"));
        assert_eq!(attribute(attrs, "transcript"), None);
    }

    #[test]
    fn streaming_gene_map_preserves_first_gtf_mapping_and_last_tsv_mapping() {
        let dir = std::env::temp_dir().join(format!("nomics-quant-map-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("test.gtf");
        let text = b"# comment\nchr1\tsrc\texon\t1\t10\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T\";\nchr1\tsrc\texon\t1\t10\t.\t+\t.\tgene_id \"G2\"; transcript_id \"T\";\n";
        std::fs::write(&path, text).unwrap();
        assert_eq!(gene_map(&path).unwrap()["T"], "G1");
        // Content detection works even when gzip has retained the GTF suffix.
        gzip(&path, text).unwrap();
        assert_eq!(gene_map(&path).unwrap()["T"], "G1");
        let path = dir.join("test.tsv");
        std::fs::write(&path, "T G1\nT G2\n").unwrap();
        assert_eq!(gene_map(&path).unwrap()["T"], "G2");
        // Chunk boundaries (any size) must not change first/last-wins or
        // the reported UTF-8 error.
        let gtf = dir.join("many.gtf");
        let mut text = String::from("# header\n");
        for i in 0..40 {
            let attrs = if i % 7 == 3 {
                format!("transcript_id \"T{}\";", i % 5)
            } else {
                format!("gene_id \"G{i}\"; transcript_id \"T{}\";", i % 5)
            };
            text.push_str(&format!("c\ts\texon\t1\t2\t.\t+\t.\t{attrs}\n"));
        }
        std::fs::write(&gtf, &text).unwrap();
        let tsv = dir.join("many.tsv");
        let lines: String = (0..40).map(|i| format!("T{} G{i}\n", i % 5)).collect();
        std::fs::write(&tsv, lines.trim_end()).unwrap();
        let bad = dir.join("bad.tsv");
        std::fs::write(&bad, b"T1 G1\nT2 G\xe2\x82\nT3 G3\n").unwrap();
        let expected = (
            gene_map_chunked(&gtf, u64::MAX).unwrap(),
            gene_map_chunked(&tsv, u64::MAX).unwrap(),
            gene_map_chunked(&bad, u64::MAX).unwrap_err().to_string(),
        );
        assert_eq!(expected.0["T0"], "G0");
        assert_eq!(expected.1["T0"], "G35");
        for size in 1..64 {
            assert_eq!(gene_map_chunked(&gtf, size).unwrap(), expected.0);
            assert_eq!(gene_map_chunked(&tsv, size).unwrap(), expected.1);
            assert_eq!(
                gene_map_chunked(&bad, size).unwrap_err().to_string(),
                expected.2
            );
        }
        std::fs::remove_dir_all(dir).unwrap();
    }
}
