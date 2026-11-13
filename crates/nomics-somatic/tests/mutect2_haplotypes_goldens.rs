//! Exact stage-two gates. Goldens are optional; differences are never normalized.
use nomics_io::{
    bam::{CIGAR_OPS, Header, RecordView, read_raw_record},
    bgzf::BgzfReader,
};
use nomics_somatic::{mutect2_assembly, mutect2_regions::Inputs};
use std::{
    fs::{self, File},
    io::{BufReader, Write},
    path::{Path, PathBuf},
};
fn projection(bam: &Path, output: &Path) {
    let mut stream = BgzfReader::new(BufReader::new(File::open(bam).unwrap()));
    let header = Header::read(&mut stream).unwrap();
    let mut raw = Vec::new();
    let mut out = std::io::BufWriter::new(File::create(output).unwrap());
    while read_raw_record(&mut stream, &mut raw).unwrap() {
        let r = RecordView::new(&raw).unwrap();
        if r.string_tag(b"RG") != Some(b"ArtificialHaplotypeRG") {
            continue;
        }
        assert_eq!(r.flag(), 16);
        assert!(r.read_name().starts_with(b"HC_"));
        let mut aux = r.aux().unwrap();
        let mut hash = None;
        while !aux.is_empty() {
            let n = nomics_io::bam::aux_value_len(aux[2], &aux[3..]).unwrap();
            match &aux[..2] {
                b"HC" => {
                    hash = Some(match aux[2] {
                        b'c' => i32::from(aux[3] as i8),
                        b'C' => i32::from(aux[3]),
                        b's' => i32::from(i16::from_le_bytes(aux[3..5].try_into().unwrap())),
                        b'S' => i32::from(u16::from_le_bytes(aux[3..5].try_into().unwrap())),
                        b'i' | b'I' => i32::from_le_bytes(aux[3..7].try_into().unwrap()),
                        _ => panic!("HC tag type"),
                    })
                }
                b"RG" | b"CR" => (),
                _ => panic!("unexpected artificial tag: {:?}", &aux[..2]),
            }
            aux = &aux[3 + n..];
        }
        let cigar: String = r
            .cigar()
            .map(|c| format!("{}{}", c >> 4, CIGAR_OPS[(c & 15) as usize] as char))
            .collect();
        let bases: Vec<_> = (0..r.len()).map(|i| r.base(i)).collect();
        writeln!(
            out,
            "{}\t{}\t{}\t{}\tHC:i:{}\tRG:Z:ArtificialHaplotypeRG\tCR:Z:{}",
            String::from_utf8_lossy(&header.refs[r.ref_id() as usize].name),
            r.pos() + 1,
            cigar,
            String::from_utf8_lossy(&bases),
            hash.unwrap(),
            String::from_utf8_lossy(r.string_tag(b"CR").unwrap())
        )
        .unwrap();
    }
    out.flush().unwrap();
}
fn difference(a: &Path, b: &Path) -> Option<String> {
    let a = fs::read(a).unwrap();
    let b = fs::read(b).unwrap();
    if a == b {
        return None;
    }
    let i = a
        .iter()
        .zip(&b)
        .position(|(a, b)| a != b)
        .unwrap_or(a.len().min(b.len()));
    Some(format!(
        "byte {i}, line {}; lengths {} / {}",
        a[..i].iter().filter(|&&b| b == b'\n').count() + 1,
        a.len(),
        b.len()
    ))
}
#[test]
fn mutect2_stage_two_all_available_contigs() {
    let home = PathBuf::from(std::env::var_os("HOME").unwrap());
    let root = std::env::var_os("NOMICS_M2_ASSEMBLY_GOLDENS")
        .map(PathBuf::from)
        .unwrap_or_else(|| home.join("omics/runs/somatic"));
    let reference = std::env::var_os("NOMICS_M2_REFERENCE")
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            home.join(
                "omics/ref/GATK.GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta",
            )
        });
    let output =
        std::env::temp_dir().join(format!("nomics-m2-assembly-test-{}", std::process::id()));
    fs::create_dir_all(&output).unwrap();
    let mut failures = Vec::new();
    for c in ["chr21", "chr22", "chr19"] {
        let d = root.join(format!("m2-{c}"));
        if !d.join("oracle/haplotypes.bam").exists() {
            eprintln!(
                "skipping Mutect2 {c}: {} absent (override NOMICS_M2_ASSEMBLY_GOLDENS)",
                d.display()
            );
            continue;
        }
        let dump = output.join(format!("{c}.assembly.tsv"));
        let records = output.join(format!("{c}.artificial.tsv"));
        let start = std::time::Instant::now();
        let (_, s) = mutect2_assembly::run(
            &Inputs {
                input: &[
                    d.join(format!("HCC1395BL_WES.{c}.bam")),
                    d.join(format!("HCC1395_WES.{c}.bam")),
                ],
                normal_sample: "HCC1395_HCC1395BL_WES",
                reference: &reference,
                intervals: &d.join(format!("{c}.bed")),
                panel_of_normals: &d.join(format!("1000g_pon.{c}.vcf.gz")),
                germline_resource: &d.join(format!("af-only-gnomad.{c}.vcf.gz")),
                regions: &output.join(format!("{c}.igv")),
                stats: &output.join(format!("{c}.stats")),
                dump_prefix: None,
            },
            &dump,
            &records,
        )
        .unwrap();
        let golden = output.join(format!("{c}.golden.tsv"));
        projection(&d.join("oracle/haplotypes.bam"), &golden);
        eprintln!(
            "{c}: {s:?}; {:.3}s; {} golden artificial records",
            start.elapsed().as_secs_f64(),
            fs::read(&golden)
                .unwrap()
                .iter()
                .filter(|&&b| b == b'\n')
                .count()
        );
        if let Some(diff) = difference(&records, &golden) {
            failures.push(format!("{c} artificial: {diff}"));
        }
        if c == "chr21" {
            if let Some(oracle) = std::env::var_os("NOMICS_M2_ASSEMBLY_ORACLE") {
                if let Some(diff) = difference(&dump, Path::new(&oracle)) {
                    failures.push(format!("chr21 internal: {diff}"));
                }
            } else {
                eprintln!(
                    "skipping internal dump gate: set NOMICS_M2_ASSEMBLY_ORACLE to the qualified Java .assembly.tsv"
                );
            }
        }
        if failures.is_empty() {
            for entry in fs::read_dir(&output).unwrap() {
                let path = entry.unwrap().path();
                if path
                    .file_name()
                    .unwrap()
                    .to_string_lossy()
                    .starts_with(&format!("{c}."))
                {
                    fs::remove_file(path).unwrap();
                }
            }
        }
    }
    assert!(
        failures.is_empty(),
        "{}; retained {}",
        failures.join("; "),
        output.display()
    );
    fs::remove_dir_all(output).unwrap();
}
