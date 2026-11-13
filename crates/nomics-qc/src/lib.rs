//! VCF QC with the pinned tools' default behavior, including their visible quirks.
//! See docs/spec/qc-vcf.md and docs/GATES.md for the two defined UB boundaries.

pub mod fastqc;
pub mod mosdepth;
pub mod bcftools;
pub mod samtools_stats;
pub mod vcftools;

use std::collections::BTreeSet;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader};
use std::path::Path;

use nomics_io::bgzf::BgzfReader;

#[derive(Default)]
pub(crate) struct Header {
    samples: usize,
    integer_info: BTreeSet<String>,
}

pub(crate) struct Record<'a> {
    fields: Vec<&'a str>,
}

impl<'a> Record<'a> {
    fn parse(line: &'a str, samples: usize) -> io::Result<Self> {
        let fields: Vec<_> = line.split('\t').collect();
        if fields.len() < 8 || (samples > 0 && fields.len() != 9 + samples) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid VCF column count",
            ));
        }
        Ok(Self { fields })
    }

    fn reference(&self) -> &str {
        self.fields[3]
    }
    fn alternates(&self) -> Vec<&str> {
        self.fields[4].split(',').collect()
    }
    fn info(&self, key: &str) -> Option<&str> {
        self.fields[7]
            .split(';')
            .filter_map(|s| s.split_once('='))
            .find_map(|(k, v)| (k == key).then_some(v))
    }
}

/// The architecture the reference tools ran on. vcftools prints its Ts/Tv
/// ratios through glibc, which writes the sign of a NaN, and 0/0 gives x86's
/// default NaN (sign bit set, `-nan`) or arm64's (clear, `nan`). sarek's x86
/// and arm64 runs therefore differ in these fields, and nomics writes what
/// sarek would have written on the named architecture.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ReferenceArch {
    Aarch64,
    X86_64,
}

impl ReferenceArch {
    /// The architecture this binary runs on (x86_64, else arm64's spelling).
    pub fn host() -> Self {
        if cfg!(target_arch = "x86_64") { Self::X86_64 } else { Self::Aarch64 }
    }

    /// `x86_64`/`amd64` or `aarch64`/`arm64`.
    pub fn parse(name: &str) -> Option<Self> {
        match name {
            "x86_64" | "amd64" => Some(Self::X86_64),
            "aarch64" | "arm64" => Some(Self::Aarch64),
            _ => None,
        }
    }

    /// `a / b` as vcftools prints it.
    pub(crate) fn ratio(self, a: f64, b: f64) -> String {
        let r = a / b;
        if r.is_nan() && self == Self::X86_64 { "-nan".into() } else { default_float(r) }
    }
}

/// Render all four reports from uncompressed VCF text. The filename is the input
/// basename (with `.vcf.gz`), independently of the output prefix.
pub fn reports_from_reader(reader: impl BufRead, filename: &str, arch: ReferenceArch) -> io::Result<[String; 4]> {
    let mut header = Header::default();
    let mut vc = vcftools::Stats::default();
    let mut bc = bcftools::Stats::default();
    let mut have_header = false;
    for line in reader.lines() {
        let line = line?;
        if let Some(meta) = line.strip_prefix("##INFO=<") {
            // Quoted descriptions may contain commas and attribute-like text.
            let mut quoted = false;
            let mut escaped = false;
            let parts: Vec<_> = meta
                .trim_end_matches('>')
                .split(|c| {
                    if escaped {
                        escaped = false;
                        return false;
                    }
                    if c == '\\' && quoted {
                        escaped = true;
                    }
                    if c == '"' {
                        quoted = !quoted;
                    }
                    c == ',' && !quoted
                })
                .map(str::trim)
                .collect();
            if parts.contains(&"Type=Integer")
                && let Some(id) = parts.iter().find_map(|s| s.strip_prefix("ID="))
            {
                header.integer_info.insert(id.to_owned());
            }
        } else if line.starts_with("#CHROM\t") {
            header.samples = line.split('\t').count().saturating_sub(9);
            vc = vcftools::Stats::new(header.samples);
            bc = bcftools::Stats::new(header.samples);
            have_header = true;
        } else if !line.starts_with('#') && !line.is_empty() {
            if !have_header {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "missing #CHROM header",
                ));
            }
            let record = Record::parse(&line, header.samples)?;
            vc.add(&record)?;
            bc.add(&record, &header)?;
        }
    }
    if !have_header {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "missing #CHROM header",
        ));
    }
    let [count, qual, filter] = vc.render(arch);
    Ok([count, qual, filter, bc.render(filename, header.samples)])
}

/// Read BGZF using nomics-io and publish the four reports under one directory.
pub fn vcf_stats(vcf: &Path, outdir: &Path, prefix: Option<&str>, arch: ReferenceArch) -> io::Result<()> {
    let filename = vcf
        .file_name()
        .and_then(|s| s.to_str())
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "VCF filename must be UTF-8"))?;
    let prefix = prefix.unwrap_or_else(|| filename.strip_suffix(".vcf.gz").unwrap_or(filename));
    if prefix.is_empty()
        || prefix.contains('/')
        || Path::new(prefix).components().count() != 1
        || prefix == "."
        || prefix == ".."
    {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "prefix must be a filename",
        ));
    }
    let input = BufReader::new(BgzfReader::new(BufReader::new(File::open(vcf)?)));
    let reports = reports_from_reader(input, filename, arch)?;
    fs::create_dir_all(outdir)?;
    for (suffix, report) in [
        "TsTv.count",
        "TsTv.qual",
        "FILTER.summary",
        "bcftools_stats.txt",
    ]
    .iter()
    .zip(reports)
    {
        fs::write(outdir.join(format!("{prefix}.{suffix}")), report)?;
    }
    Ok(())
}

/// A field exempted by the sole `.TsTv.qual` normalization. Row and column are
/// one-based file coordinates, including the header. Equal fields are retained.
#[derive(Debug, PartialEq, Eq)]
pub struct ExemptField {
    pub row: usize,
    pub column: usize,
    pub nomics: String,
    pub golden: String,
}

/// Compare every byte except authorized fields, using *nomics* zeros to select
/// exemptions. Requires identical header, row count, separators and newlines.
pub fn compare_tstv_qual(nomics: &str, golden: &str) -> Result<Vec<ExemptField>, String> {
    let ours: Vec<_> = nomics.split_inclusive('\n').collect();
    let theirs: Vec<_> = golden.split_inclusive('\n').collect();
    if ours.len() != theirs.len() {
        return Err(format!("row count: {} != {}", ours.len(), theirs.len()));
    }
    if ours.is_empty() {
        return Err("missing header".into());
    }
    let mut exempt = Vec::new();
    for (row, (a, b)) in ours.iter().zip(&theirs).enumerate() {
        if row == 0 {
            if a != b {
                return Err("line 1: header differs".into());
            }
            continue;
        }
        if a.ends_with('\n') != b.ends_with('\n') {
            return Err(format!("line {}: newline differs", row + 1));
        }
        let aa: Vec<_> = a.trim_end_matches('\n').split('\t').collect();
        let bb: Vec<_> = b.trim_end_matches('\n').split('\t').collect();
        if aa.len() != 7 || bb.len() != 7 {
            return Err(format!("line {}: expected seven columns", row + 1));
        }
        for col in 0..7 {
            let allowed = match col {
                4 | 5 => aa[col] == "0",
                6 => aa[4] == "0" || aa[5] == "0",
                _ => false,
            };
            if allowed {
                exempt.push(ExemptField {
                    row: row + 1,
                    column: col + 1,
                    nomics: aa[col].into(),
                    golden: bb[col].into(),
                });
            } else if aa[col] != bb[col] {
                return Err(format!(
                    "line {}, column {}: {:?} != {:?}",
                    row + 1,
                    col + 1,
                    aa[col],
                    bb[col]
                ));
            }
        }
    }
    Ok(exempt)
}

/// C++ defaultfloat: six significant digits, exponent chosen after rounding,
/// no trailing fraction zeros, exponent sign and at least two exponent digits.
pub(crate) fn default_float(value: f64) -> String {
    if value.is_nan() {
        return "nan".into();
    }
    if value.is_infinite() {
        return if value.is_sign_negative() {
            "-inf"
        } else {
            "inf"
        }
        .into();
    }
    if value == 0.0 {
        return if value.is_sign_negative() { "-0" } else { "0" }.into();
    }
    let scientific = format!("{value:.5e}");
    let (mantissa, exponent) = scientific.split_once('e').unwrap();
    let exponent: i32 = exponent.parse().unwrap();
    if !(-4..6).contains(&exponent) {
        return format!(
            "{}e{exponent:+03}",
            mantissa.trim_end_matches('0').trim_end_matches('.')
        );
    }
    let negative = mantissa.starts_with('-');
    let digits = mantissa.trim_start_matches('-').replace('.', "");
    let point = exponent + 1;
    let mut result = if point <= 0 {
        format!("0.{}{}", "0".repeat((-point) as usize), digits)
    } else {
        let p = point as usize;
        if p >= digits.len() {
            digits + &"0".repeat(p - 6)
        } else {
            format!("{}.{}", &digits[..p], &digits[p..])
        }
    };
    if result.contains('.') {
        result = result.trim_end_matches('0').trim_end_matches('.').into();
    }
    if negative {
        result.insert(0, '-');
    }
    result
}

#[cfg(test)]
mod tests {
    /// vcftools prints `ts / tv` with operator<<(double); glibc writes a NaN's
    /// sign, and 0.0/0.0 is negative on x86 (the x86 goldens' `-nan`).
    #[test]
    fn zero_over_zero_follows_the_reference_architecture() {
        assert_eq!(super::ReferenceArch::Aarch64.ratio(0.0, 0.0), "nan");
        assert_eq!(super::ReferenceArch::X86_64.ratio(0.0, 0.0), "-nan");
        assert_eq!(super::ReferenceArch::X86_64.ratio(1.0, 0.0), "inf");
        assert_eq!(super::ReferenceArch::X86_64.ratio(27.0, 10.0), "2.7");
    }

    use super::*;
    #[test]
    fn stream_formatting() {
        for (v, s) in [
            (2.0 / 3.0, "0.666667"),
            (999999.5, "1e+06"),
            (0.00009999999, "0.0001"),
            (1e-5, "1e-05"),
            (1e20, "1e+20"),
            (-1.0, "-1"),
        ] {
            assert_eq!(default_float(v), s);
        }
    }
    #[test]
    fn normalization_is_field_scoped_and_asymmetric() {
        let header = "a\tb\tc\td\te\tf\tg\n";
        let a = format!("{header}5\t1\t2\t0.5\t0\t1\t0\n6\t2\t2\t1\t0\t0\tnan\n");
        let b = a.replace("\t0\t1\t0\n", "\tx\t1\ty\n");
        assert_eq!(compare_tstv_qual(&a, &b).unwrap().len(), 5);
        assert!(compare_tstv_qual(&b, &a).is_err());
        for bad in [
            a.replacen("5\t", "7\t", 1),
            a.replacen("\t1\t2", "\t9\t2", 1),
            a.trim_end().to_owned(),
            a.clone() + "\n",
            a.replacen("a\t", "z\t", 1),
        ] {
            assert!(compare_tstv_qual(&a, &bad).is_err());
        }
    }
}

#[cfg(test)]
mod inline_vcf_tests {
    use super::*;
    fn reports(body: &str) -> [String; 4] {
        let header = "##fileformat=VCFv4.2\n##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Count\">\n##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Number\">\n##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n";
        reports_from_reader((header.to_owned() + body).as_bytes(), "test.vcf.gz", ReferenceArch::Aarch64).unwrap()
    }
    #[test]
    fn vcftools_thresholds_count_boundary_missing_quality_and_filter_ties() {
        // output_TsTv_by_quality: lower excludes the current bin; upper uses
        // ui+1 in its 1-based cumulative vector. Boundary N+1 is zero.
        // output_TsTv_by_count ignores the authorized invalid maximum-AC write.
        // output_FILTER_summary uses the first ALT even on multiallelic sites;
        // vcf_entry::set_FILTER sorts flags, ties reverse the resulting string.
        let [count, qual, filter, _] = reports(concat!(
            "1\t1\t.\ta\tg\t.\tz;a\tAC=0;AN=4\tGT\t1|1\t1/1\n",
            "1\t2\t.\tA\tG\t10\tPASS\tAC=0;AN=4\tGT\t1\t./.\n",
            "1\t3\t.\tA\tC\t20\tPASS\tAC=0;AN=4\tGT\t0/.\t0/0\n",
            "1\t4\t.\tC\tT,A\t30\tz;a\tAC=1,1;AN=4\tGT\t0/1\t0/2\n",
        ));
        assert_eq!(
            count,
            "ALT_ALLELE_COUNT\tN_Ts\tN_Tv\tTs/Tv\n0\t0\t1\t0\n1\t1\t0\tinf\n2\t0\t0\tnan\n3\t0\t0\tnan\n"
        );
        let rows: Vec<_> = qual.lines().skip(1).collect();
        assert_eq!(
            rows,
            [
                "-1\t0\t0\tnan\t1\t1\t1",
                "10\t1\t0\tinf\t0\t1\t0",
                "20\t2\t0\tinf\t0\t0\tnan"
            ]
        );
        assert_eq!(
            filter,
            "FILTER\tN_VARIANTS\tN_Ts\tN_Tv\tTs/Tv\na;z\t2\t2\t0\tinf\nPASS\t2\t1\t1\t1\n"
        );
    }
    #[test]
    fn bcftools_mixed_alleles_and_info_not_genotypes() {
        // do_vcf_stats counts each site type once, do_snp_stats counts each
        // ALT (but QUAL only ALT1), init_iaf/bcf_calc_ac use INFO only. No
        // FORMAT/AD/DP statistics or sample sections without -s/-S.
        let [_, _, _, stats] = reports(concat!(
            "1\t1\t.\tA\tG,C\t10\tq10\tAC=1,2;AN=4;DP=501\tGT:AD:DP\t0/0:9,0,0:9\t0/0:9,0,0:9\n",
            "1\t2\t.\tA\tT,ATG\t.\t.\tAC=0,1;AN=4;DP=.\tGT\t0/0\t0/0\n",
            "1\t3\t.\tCC\tTT\t20\tPASS\tDP=0\tGT\t0/0\t0/0\n",
            "1\t4\t.\tACT\tTCGA\t20\tPASS\tDP=500\tGT\t0/0\t0/0\n",
            "1\t5\t.\tA\t<DEL>\t20\tPASS\t.\tGT\t0/0\t0/0\n",
            "1\t6\t.\tA\t*\t20\tPASS\t.\tGT\t0/0\t0/0\n",
            "1\t7\t.\tA\t<NON_REF>\t20\tPASS\t.\tGT\t0/0\t0/0\n",
            "1\t8\t.\tA\tA[chr2:9[\t20\tPASS\t.\tGT\t0/0\t0/0\n",
            "1\t9\t.\tGAC\tGTC\t20\tPASS\t.\tGT\t0/0\t0/0\n",
            "1\t10\t.\tAG\t.\t20\tPASS\t.\tGT\t0/0\t0/0\n",
        ));
        for expected in [
            "SN\t0\tnumber of records:\t10",
            "SN\t0\tnumber of no-ALTs:\t2",
            "SN\t0\tnumber of SNPs:\t3",
            "SN\t0\tnumber of MNPs:\t1",
            "SN\t0\tnumber of indels:\t1",
            "SN\t0\tnumber of others:\t2",
            "SN\t0\tnumber of multiallelic sites:\t2",
            "SN\t0\tnumber of multiallelic SNP sites:\t1",
            "TSTV\t0\t1\t2\t0.50\t1\t1\t1.00",
            "SiS\t0\t1\t1\t1\t0\t1\t0\t0\t1",
            "AF\t0\t0.000000\t2\t1\t1\t1\t0\t0\t1",
            "AF\t0\t0.490000\t1\t0\t1\t0\t0\t0\t0",
            "QUAL\t0\t.\t1\t0\t1\t1",
            "QUAL\t0\t10.0\t1\t1\t0\t0",
            "IDD\t0\t2\t1\t0\t.",
            "ST\t0\tA>C\t1",
            "ST\t0\tA>G\t1",
            "DP\t0\t<0\t0\t0.000000\t1\t25.000000",
            "DP\t0\t>500\t0\t0.000000\t1\t25.000000",
        ] {
            assert!(stats.lines().any(|l| l == expected), "missing {expected}");
        }
        assert!(!stats.contains("QUAL\t0\t20.0")); // Padded SNP: first bases identical.
        assert!(!stats.contains("PSC\t"));
        assert!(stats.contains("ID\t0\ttest.vcf.gz\n"));
    }
    #[test]
    fn bcftools_unknown_af_missing_ac_and_length_cap() {
        // init_iaf reserves singleton/unknown bin0; missing AC is an integer
        // sentinel and clamps to AF0 (bin1), not the singleton bin. dist/IDD
        // caps long indels at 60. A SNP to N counts a site but no substitution.
        let long = "A".repeat(70);
        let text = format!(
            "1\t1\t.\tA\tG\t.\t.\t.\tGT\t1/1\t1/1\n1\t2\t.\tA\tC\t0\t.\tAC=.;AN=0\tGT\t0/0\t0/0\n1\t3\t.\tA\tN\t10\t.\t.\tGT\t0/0\t0/0\n1\t4\t.\tA\t{long}\t1\t.\tAC=2;AN=4\tGT\t0/0\t0/0\n"
        );
        let stats = &reports(&text)[3];
        assert!(stats.contains("SiS\t0\t1\t1\t1\t0\t0\t0\t0\t0\n"));
        assert!(stats.contains("IDD\t0\t60\t1\t0\t.\n"));
        assert!(stats.contains("SN\t0\tnumber of SNPs:\t3\n"));
    }
}
