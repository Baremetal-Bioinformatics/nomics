//! VCF output exactly as GATK 4.6.2.0 HaplotypeCaller writes it through
//! htsjdk 4.2.0 (`VCFWriter`, `VCFEncoder`, `BlockCompressedOutputStream`).
//! Spec: `docs/spec/haplotypecaller.md` §9.
//!
//! - Header: `##fileformat=VCFv4.2`, then the metadata lines in htsjdk's
//!   sorted order (by line text; `##contig` lines among themselves by
//!   dictionary index), then the `#CHROM` line.
//! - INFO keys alphabetical; pre-formatted annotation strings as given;
//!   `Double` values through `formatVCFDouble`; flags without `=`.
//! - QUAL `%.2f` with a trailing `.00` stripped; `GQ` capped at 99.
//! - BGZF in 65,498-byte blocks, zlib-exact level 2 (htsjdk profile).

use std::io::{self, Write};

use nomics_io::bgzf::ParallelBgzfWriter;

use crate::javafmt::{format_fixed, format_sci};

/// An INFO value, typed as GATK hands it to htsjdk.
#[derive(Debug, Clone, PartialEq)]
pub enum InfoValue {
    /// `Number=0` flag (e.g. `DB`): written as the bare key.
    Flag,
    /// An annotation that formats itself (e.g. `QD` with `%.2f`).
    Str(String),
    Ints(Vec<i64>),
    /// Written with `formatVCFDouble`, comma-joined.
    Doubles(Vec<f64>),
}

/// The single sample's genotype (`GT:AD:DP:GQ:PL`).
#[derive(Debug, Clone, PartialEq)]
pub struct Genotype {
    /// Allele indices, unphased.
    pub alleles: Vec<usize>,
    /// Empty when the genotype's AD annotation is absent.
    pub ad: Vec<i32>,
    /// -1 when absent, as in htsjdk's Genotype.NO_DP.
    pub dp: i32,
    /// Written as `min(gq, 99)`.
    pub gq: i32,
    pub pl: Vec<i32>,
    /// Reference block statistics or variant annotations in GVCF mode.
    pub gvcf: Option<GvcfGenotype>,
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct GvcfGenotype {
    pub min_dp: Option<i32>,
    pub sb: Option<[i32; 4]>,
    /// (PGT, PID, PS); presence also makes GT phased.
    pub phase: Option<(String, String, i64)>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Record {
    /// Index into the header's contig list.
    pub contig: usize,
    /// 1-based.
    pub pos: i64,
    pub id: Option<String>,
    pub ref_allele: String,
    pub alts: Vec<String>,
    /// Phred-scaled (`-10 * log10PError + 0.0`).
    pub qual: f64,
    pub info: Vec<(String, InfoValue)>,
    pub genotype: Genotype,
}

/// htsjdk `VCFEncoder.formatVCFDouble`.
pub fn format_vcf_double(d: f64) -> String {
    if d < 1.0 {
        if d < 0.01 {
            if d.abs() >= 1e-20 {
                format_sci(d, 3)
            } else {
                "0.00".to_string()
            }
        } else {
            format_fixed(d, 3)
        }
    } else {
        format_fixed(d, 2)
    }
}

/// htsjdk `VCFEncoder.formatQualValue`.
pub fn format_qual(q: f64) -> String {
    if q.is_nan() {
        return ".".into();
    }
    let s = format_fixed(q, 2);
    match s.strip_suffix(".00") {
        Some(t) => t.to_string(),
        None => s,
    }
}

/// `htsjdk.variant.vcf.VCFHeader` content for HaplotypeCaller's output.
pub struct Header {
    /// Metadata lines without the leading `##` (e.g. `source=HaplotypeCaller`),
    /// excluding `contig` lines.
    pub meta: Vec<String>,
    /// (name, length) in sequence-dictionary order.
    pub contigs: Vec<(String, u64)>,
    pub sample: String,
}

impl Header {
    /// The header text, in `VCFWriter.writeHeader` order.
    pub fn text(&self) -> String {
        let mut meta = self.meta.clone();
        meta.sort();
        let mut s = String::from("##fileformat=VCFv4.2\n");
        // htsjdk's TreeSet orders lines by their text; contig lines compare
        // among themselves by dictionary index, so they sit as one block
        // where "contig=" sorts.
        let mut contigs_done = false;
        let put_contigs = |s: &mut String| {
            for (name, len) in &self.contigs {
                s.push_str(&format!("##contig=<ID={name},length={len}>\n"));
            }
        };
        for m in &meta {
            if !contigs_done && m.as_str() > "contig=" {
                put_contigs(&mut s);
                contigs_done = true;
            }
            s.push_str("##");
            s.push_str(m);
            s.push('\n');
        }
        if !contigs_done {
            put_contigs(&mut s);
        }
        s.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t");
        s.push_str(&self.sample);
        s.push('\n');
        s
    }
}

fn join<T: ToString>(v: &[T]) -> String {
    v.iter()
        .map(|x| x.to_string())
        .collect::<Vec<_>>()
        .join(",")
}

/// `VCFEncoder.encode` for one record, with `\n`.
pub fn encode(rec: &Record, contigs: &[(String, u64)], out: &mut Vec<u8>) {
    let mut line = String::with_capacity(256);
    line.push_str(&contigs[rec.contig].0);
    line.push('\t');
    line.push_str(&rec.pos.to_string());
    line.push('\t');
    line.push_str(rec.id.as_deref().unwrap_or("."));
    line.push('\t');
    line.push_str(&rec.ref_allele);
    line.push('\t');
    line.push_str(if rec.alts.is_empty() { "." } else { "" });
    line.push_str(&rec.alts.join(","));
    line.push('\t');
    line.push_str(&format_qual(rec.qual));
    line.push_str("\t.\t");
    let mut info: Vec<&(String, InfoValue)> = rec.info.iter().collect();
    info.sort_by(|a, b| a.0.cmp(&b.0));
    if info.is_empty() {
        line.push('.');
    }
    for (i, (k, v)) in info.iter().enumerate() {
        if i > 0 {
            line.push(';');
        }
        line.push_str(k);
        match v {
            InfoValue::Flag => {}
            InfoValue::Str(s) => {
                line.push('=');
                line.push_str(s);
            }
            InfoValue::Ints(x) => {
                line.push('=');
                line.push_str(&join(x));
            }
            InfoValue::Doubles(x) => {
                line.push('=');
                line.push_str(
                    &x.iter()
                        .map(|&d| format_vcf_double(d))
                        .collect::<Vec<_>>()
                        .join(","),
                );
            }
        }
    }
    let g = &rec.genotype;
    let extra = g.gvcf.as_ref();
    let block = extra.and_then(|x| x.min_dp);
    let phase = extra.and_then(|x| x.phase.as_ref());
    let sb = extra.and_then(|x| x.sb);
    let mut keys = vec!["GT"];
    let mut values = vec![
        g.alleles
            .iter()
            .map(|&a| {
                if a == crate::genotype::NO_CALL_ALLELE {
                    ".".into()
                } else {
                    a.to_string()
                }
            })
            .collect::<Vec<_>>()
            .join(if phase.is_some() { "|" } else { "/" }),
    ];
    if block.is_none() && !g.ad.is_empty() {
        keys.push("AD");
        values.push(join(&g.ad));
    }
    if g.dp >= 0 {
        keys.push("DP");
        values.push(g.dp.to_string());
    }
    keys.push("GQ");
    values.push(g.gq.min(99).to_string());
    if let Some(min_dp) = block {
        keys.push("MIN_DP");
        values.push(min_dp.to_string());
    }
    if let Some((pgt, pid, _)) = phase {
        keys.extend(["PGT", "PID"]);
        values.extend([pgt.clone(), pid.clone()]);
    }
    keys.push("PL");
    values.push(join(&g.pl));
    if let Some((_, _, ps)) = phase {
        keys.push("PS");
        values.push(ps.to_string());
    }
    if let Some(sb) = sb {
        keys.push("SB");
        values.push(join(&sb));
    }
    line.push('\t');
    line.push_str(&keys.join(":"));
    line.push('\t');
    line.push_str(&values.join(":"));
    line.push('\n');
    out.extend_from_slice(line.as_bytes());
}

/// Write a whole `.vcf.gz` (BGZF, htsjdk profile) on `threads` threads.
pub fn write_vcf_gz<W: Write + Send + 'static>(
    out: W,
    header: &Header,
    records: &[Record],
    threads: usize,
) -> io::Result<W> {
    let mut w = ParallelBgzfWriter::new(out, threads);
    w.write_all(header.text().as_bytes())?;
    let mut buf = Vec::with_capacity(1 << 20);
    for r in records {
        encode(r, &header.contigs, &mut buf);
        if buf.len() >= 1 << 20 {
            w.write_all(&buf)?;
            buf.clear();
        }
    }
    w.write_all(&buf)?;
    w.finish()
}

/// Parse a VCF written by this module's rules back into records (for the
/// round-trip gate). INFO values are typed as GATK types them.
pub fn parse_record(
    line: &str,
    contig_index: &std::collections::HashMap<String, usize>,
) -> io::Result<Record> {
    let bad = |m: &str| io::Error::new(io::ErrorKind::InvalidData, format!("{m}: {line}"));
    let f: Vec<&str> = line.split('\t').collect();
    if f.len() != 10 {
        return Err(bad("expected 10 columns"));
    }
    let contig = *contig_index
        .get(f[0])
        .ok_or_else(|| bad("unknown contig"))?;
    let num = |s: &str| s.parse::<i64>().map_err(|_| bad("bad integer"));
    let mut info = Vec::new();
    for kv in f[7].split(';') {
        let (k, v) = match kv.split_once('=') {
            Some((k, v)) => (k, Some(v)),
            None => (kv, None),
        };
        let value = match (k, v) {
            (_, None) => InfoValue::Flag,
            ("AC" | "AN" | "MLEAC", Some(v)) => {
                InfoValue::Ints(v.split(',').map(num).collect::<Result<_, _>>()?)
            }
            ("AF" | "MLEAF", Some(v)) => InfoValue::Doubles(
                v.split(',')
                    .map(|x| x.parse::<f64>().map_err(|_| bad("bad float")))
                    .collect::<Result<_, _>>()?,
            ),
            (_, Some(v)) => InfoValue::Str(v.to_string()),
        };
        info.push((k.to_string(), value));
    }
    if f[8] != "GT:AD:DP:GQ:PL" {
        return Err(bad("unexpected FORMAT"));
    }
    let g: Vec<&str> = f[9].split(':').collect();
    let ints = |s: &str| -> io::Result<Vec<i32>> {
        s.split(',')
            .map(|x| x.parse::<i32>().map_err(|_| bad("bad int")))
            .collect()
    };
    Ok(Record {
        contig,
        pos: num(f[1])?,
        id: (f[2] != ".").then(|| f[2].to_string()),
        ref_allele: f[3].to_string(),
        alts: f[4].split(',').map(str::to_string).collect(),
        qual: f[5].parse().map_err(|_| bad("bad QUAL"))?,
        info,
        genotype: Genotype {
            alleles: g[0]
                .split('/')
                .map(|x| x.parse().map_err(|_| bad("bad GT")))
                .collect::<Result<_, _>>()?,
            ad: ints(g[1])?,
            dp: g[2].parse().map_err(|_| bad("bad DP"))?,
            gq: g[3].parse().map_err(|_| bad("bad GQ"))?,
            pl: ints(g[4])?,
            gvcf: None,
        },
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn header_order_follows_htsjdk() {
        let h = Header {
            meta: vec![
                "source=HaplotypeCaller".into(),
                "INFO=<ID=AC>".into(),
                "FILTER=<ID=LowQual>".into(),
                "GATKCommandLine=<ID=HaplotypeCaller>".into(),
                "FORMAT=<ID=GT>".into(),
            ],
            contigs: vec![("chr2".into(), 2), ("chr10".into(), 10)],
            sample: "S".into(),
        };
        let text = h.text();
        let lines: Vec<&str> = text.lines().map(|l| l.split('=').next().unwrap()).collect();
        assert_eq!(
            lines,
            [
                "##fileformat",
                "##FILTER",
                "##FORMAT",
                "##GATKCommandLine",
                "##INFO",
                "##contig",
                "##contig",
                "##source",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS"
            ]
        );
        // contigs stay in dictionary order, not text order
        assert!(h.text().find("ID=chr2,").unwrap() < h.text().find("ID=chr10,").unwrap());
    }
}
