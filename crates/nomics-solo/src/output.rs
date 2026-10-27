use crate::{
    Config, Stats,
    barcode::{Whitelist, decode},
};
use std::fs::{self, File};
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::Path;

const BAR_NAMES: [&str; 12] = [
    "noNoAdapter",
    "noNoUMI",
    "noNoCB",
    "noNinCB",
    "noNinUMI",
    "noUMIhomopolymer",
    "noNoWLmatch",
    "noTooManyMM",
    "noTooManyWLmatches",
    "yesWLmatchExact",
    "yesOneWLmatchWithMM",
    "yesMultWLmatchWithMM",
];
const FEAT_NAMES: [&str; 11] = [
    "noUnmapped",
    "noNoFeature",
    "MultiFeature",
    "subMultiFeatureMultiGenomic",
    "noTooManyWLmatches",
    "noMMtoWLwithoutExact",
    "yesWLmatch",
    "yessubWLmatchExact",
    "yessubWLmatch_UniqueFeature",
    "yesCellBarcodes",
    "yesUMIs",
];

fn stats_file(path: &Path, names: &[&str], counts: &[u64]) -> io::Result<()> {
    let mut out = BufWriter::new(File::create(path)?);
    for (name, count) in names.iter().zip(counts) {
        writeln!(out, "{name:>50}{count:>15}")?;
    }
    out.flush()
}

// C++ ostream defaultfloat, precision 6. Round once to six significant digits,
// then choose fixed/scientific using the exponent of the rounded result (%g).
fn general(value: f64) -> String {
    if value.is_nan() {
        return "-nan".into();
    }
    if value == f64::INFINITY {
        return "inf".into();
    }
    if value == 0.0 {
        return "0".into();
    }
    let scientific = format!("{value:.5e}");
    let (mantissa, exponent) = scientific.split_once('e').unwrap();
    let exp: i32 = exponent.parse().unwrap();
    if !(-4..6).contains(&exp) {
        let m = mantissa.trim_end_matches('0').trim_end_matches('.');
        return format!("{m}e{exp:+03}");
    }
    let digits = (5 - exp).max(0) as usize;
    let result = format!("{value:.digits$}");
    if result.contains('.') {
        result.trim_end_matches('0').trim_end_matches('.').into()
    } else {
        result
    }
}

pub(crate) fn threshold(sorted: &[u32], tentative_cells: usize) -> u32 {
    if tentative_cells == 0 {
        return 1;
    }
    // C++ round: positive halves away from zero; robust max uses zero-based rank.
    let rank = (1000.0f64 * (1.0 - 0.99)).round() as usize;
    // Zero-count tentative CBs remain in STAR's nCB and in its sorted array.
    let robust_max = sorted
        .get(rank.min(tentative_cells - 1))
        .copied()
        .unwrap_or(0);
    (f64::from(robust_max) / 10.0).round().max(1.0) as u32
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn write(
    config: &Config,
    features: &str,
    wl: &Whitelist,
    s: &Stats,
    matrix_path: &Path,
    entries: u64,
    reads: &[u32],
    umis: &[u32],
    genes: &[u32],
    tentative_cells: usize,
) -> io::Result<usize> {
    let gene_dir = config.output.join("Gene");
    fs::create_dir_all(gene_dir.join("raw"))?;
    fs::create_dir_all(gene_dir.join("filtered"))?;
    stats_file(&config.output.join("Barcodes.stats"), &BAR_NAMES, &s.bar)?;
    stats_file(&gene_dir.join("Features.stats"), &FEAT_NAMES, &s.feat)?;
    let nfeatures = features.lines().count();
    let mut sorted: Vec<u32> = umis.iter().copied().filter(|&u| u > 0).collect();
    sorted.sort_unstable_by(|a, b| b.cmp(a));
    let cutoff = threshold(&sorted, tentative_cells);
    let cells = umis.iter().filter(|&&u| u >= cutoff).count();
    let filtered_entries: u64 = genes
        .iter()
        .zip(umis)
        .filter(|&(_, u)| *u >= cutoff)
        .map(|(&g, _)| u64::from(g))
        .sum();
    let mut rank = nomics_mem::HostArray::<u32>::zeroed(wl.len)?;
    let mut cell_reads = Vec::new();
    let mut cell_genes = Vec::new();
    let (mut reads_in_cells, mut umis_in_cells, mut genes_in_cells) = (0u64, 0u64, 0u64);
    let mut raw_bc = BufWriter::new(File::create(gene_dir.join("raw/barcodes.tsv"))?);
    let mut filtered_bc = BufWriter::new(File::create(gene_dir.join("filtered/barcodes.tsv"))?);
    let mut r = 0;
    for cb in 0..wl.len {
        let barcode = decode(wl.values[cb], 16);
        writeln!(raw_bc, "{barcode}")?;
        if umis[cb] >= cutoff {
            r += 1;
            rank[cb] = r;
            writeln!(filtered_bc, "{barcode}")?;
            reads_in_cells += u64::from(reads[cb]);
            umis_in_cells += u64::from(umis[cb]);
            genes_in_cells += u64::from(genes[cb]);
            cell_reads.push(reads[cb]);
            cell_genes.push(genes[cb]);
        }
    }
    raw_bc.flush()?;
    filtered_bc.flush()?;
    for dir in ["raw", "filtered"] {
        fs::write(gene_dir.join(dir).join("features.tsv"), features)?;
    }
    let mut raw = BufWriter::with_capacity(1 << 20, File::create(gene_dir.join("raw/matrix.mtx"))?);
    let mut filtered = BufWriter::new(File::create(gene_dir.join("filtered/matrix.mtx"))?);
    write!(
        raw,
        "%%MatrixMarket matrix coordinate integer general\n%\n{nfeatures} {} {entries}\n",
        wl.len
    )?;
    write!(
        filtered,
        "%%MatrixMarket matrix coordinate integer general\n%\n{nfeatures} {cells} {filtered_entries}\n"
    )?;
    let mut input = BufReader::new(File::open(matrix_path)?);
    let mut detected = vec![false; nfeatures];
    for _ in 0..entries {
        let mut record = [0; 12];
        input.read_exact(&mut record)?;
        let cb = u32::from_le_bytes(record[..4].try_into().unwrap()) as usize;
        let g = u32::from_le_bytes(record[4..8].try_into().unwrap());
        let count = u32::from_le_bytes(record[8..].try_into().unwrap());
        writeln!(raw, "{} {} {count}", g + 1, cb + 1)?;
        if rank[cb] > 0 {
            writeln!(filtered, "{} {} {count}", g + 1, rank[cb])?;
            detected[g as usize] = true;
        }
    }
    raw.flush()?;
    filtered.flush()?;
    let mut umi_file = BufWriter::new(File::create(gene_dir.join("UMIperCellSorted.txt"))?);
    for u in &sorted {
        writeln!(umi_file, "{u}")?;
    }
    umi_file.flush()?;
    cell_reads.sort_unstable();
    cell_genes.sort_unstable();
    let mut summary = BufWriter::new(File::create(gene_dir.join("Summary.csv"))?);
    writeln!(summary, "Number of Reads,{}", s.reads)?;
    let invalid: u64 = s.bar[..9].iter().sum::<u64>() + s.feat[4] + s.feat[5];
    for (label, value) in [
        (
            "Reads With Valid Barcodes",
            1.0 - invalid as f64 / s.reads as f64,
        ),
        (
            "Sequencing Saturation",
            1.0 - s.feat[10] as f64 / s.feat[8] as f64,
        ),
        ("Q30 Bases in CB+UMI", s.cb_q30 as f64 / s.cb_bases as f64),
        (
            "Q30 Bases in RNA read",
            s.rna_q30 as f64 / s.rna_bases as f64,
        ),
        (
            "Reads Mapped to Genome: Unique+Multiple",
            s.mapped as f64 / s.reads as f64,
        ),
        (
            "Reads Mapped to Genome: Unique",
            s.unique as f64 / s.reads as f64,
        ),
    ] {
        writeln!(summary, "{label},{}", general(value))?;
    }
    writeln!(
        summary,
        "Reads Mapped to Gene: Unique+Multiple Gene,NoMulti"
    )?;
    writeln!(
        summary,
        "Reads Mapped to Gene: Unique Gene,{}",
        general(s.feat[8] as f64 / s.reads as f64)
    )?;
    writeln!(summary, "Estimated Number of Cells,{cells}")?;
    writeln!(
        summary,
        "Unique Reads in Cells Mapped to Gene,{reads_in_cells}"
    )?;
    writeln!(
        summary,
        "Fraction of Unique Reads in Cells,{}",
        general(reads_in_cells as f64 / s.feat[8] as f64)
    )?;
    let mean = |n: u64| if cells == 0 { 0 } else { n / cells as u64 };
    for (label, value) in [
        ("Mean Reads per Cell", mean(reads_in_cells)),
        (
            "Median Reads per Cell",
            u64::from(*cell_reads.get(cells / 2).unwrap_or(&0)),
        ),
        ("UMIs in Cells", umis_in_cells),
        ("Mean UMI per Cell", mean(umis_in_cells)),
        (
            "Median UMI per Cell",
            if cells == 0 {
                0
            } else {
                u64::from(sorted[cells / 2])
            },
        ),
        ("Mean Gene per Cell", mean(genes_in_cells)),
        (
            "Median Gene per Cell",
            u64::from(*cell_genes.get(cells / 2).unwrap_or(&0)),
        ),
        (
            "Total Gene Detected",
            detected.iter().filter(|&&x| x).count() as u64,
        ),
    ] {
        writeln!(summary, "{label},{value}")?;
    }
    summary.flush()?;
    eprintln!("solo-count: CellRanger2.2 cutoff {cutoff}, {cells} cells");
    Ok(cells)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn knee_rank_rounds_half_away_from_zero() {
        assert_eq!(threshold(&[999; 10], 10), 100);
        let mut counts = vec![999; 12];
        counts[10] = 145;
        assert_eq!(threshold(&counts, counts.len()), 15);
        assert_eq!(threshold(&[1], 1), 1);
        assert_eq!(threshold(&[1500], 1), 150);
        assert_eq!(threshold(&[1500], 2), 1); // candidate with no surviving reads
        assert_eq!(threshold(&[], 0), 1);
    }
    #[test]
    fn cpp_default_precision() {
        for (v, s) in [
            (0.9349196, "0.93492"),
            (0.00001, "1e-05"),
            (1_000_000.0, "1e+06"),
            (0.9999999, "1"),
            (0.9757761, "0.975776"),
            (0.00009999999, "0.0001"),
        ] {
            assert_eq!(general(v), s);
        }
    }
}
