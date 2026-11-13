use nomics_hc::assembly::{self, AsmRead, Haplotype};
use nomics_hc::javafmt;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

fn fixtures_dir() -> String {
    std::env::var("NOMICS_FIXTURES").unwrap_or_else(|_| "/nonexistent-nomics-fixtures".into())
}

fn default_fixtures() -> String {
    format!("{}/hc_finalize", fixtures_dir())
}

fn fixture_dir() -> PathBuf {
    std::env::var_os("NOMICS_HC_ASSEMBLY_FIXTURES")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from(default_fixtures()))
}

fn next_line(reader: &mut BufReader<File>) -> Option<String> {
    let mut line = String::new();
    let count = reader.read_line(&mut line).expect("read fixture line");
    if count == 0 {
        None
    } else {
        while matches!(line.as_bytes().last(), Some(b'\n' | b'\r')) {
            line.pop();
        }
        Some(line)
    }
}

fn parse_cigar(text: &str) -> Vec<u32> {
    let mut result = Vec::new();
    let mut length = 0_u32;
    for byte in text.bytes() {
        if byte.is_ascii_digit() {
            length = length * 10 + u32::from(byte - b'0');
        } else {
            let op = match byte {
                b'M' => 0,
                b'I' => 1,
                b'D' => 2,
                b'N' => 3,
                b'S' => 4,
                b'H' => 5,
                b'P' => 6,
                b'=' => 7,
                b'X' => 8,
                _ => panic!("invalid CIGAR operator {byte}"),
            };
            result.push((length << 4) | op);
            length = 0;
        }
    }
    assert_eq!(length, 0, "CIGAR ends in a length");
    result
}

fn cigar_text(cigar: &[u32]) -> String {
    let mut result = String::new();
    for value in cigar {
        result.push_str(&(value >> 4).to_string());
        result.push(match value & 15 {
            0 => 'M',
            1 => 'I',
            2 => 'D',
            3 => 'N',
            4 => 'S',
            5 => 'H',
            6 => 'P',
            7 => '=',
            8 => 'X',
            op => panic!("invalid BAM CIGAR operator {op}"),
        });
    }
    result
}

fn haplotype_line(haplotype: &Haplotype) -> String {
    format!(
        "H {} {} {} {}:{}-{} {} {} {}",
        haplotype.is_ref,
        haplotype.alignment_start_hap_wrt_ref,
        cigar_text(&haplotype.cigar),
        haplotype.contig,
        haplotype.start,
        haplotype.end,
        javafmt::double_to_string(haplotype.score),
        haplotype.kmer_size,
        String::from_utf8_lossy(&haplotype.bases),
    )
}

fn parse_location(location: &str) -> (&str, i64, i64) {
    let (contig, positions) = location.rsplit_once(':').expect("location colon");
    let (start, end) = positions.split_once('-').expect("location dash");
    (
        contig,
        start.parse().expect("location start"),
        end.parse().expect("location end"),
    )
}

fn open(path: &Path) -> BufReader<File> {
    BufReader::with_capacity(
        1024 * 1024,
        File::open(path).expect("open assembly fixture"),
    )
}

#[test]
#[ignore = "requires the 75 MiB GATK local-assembly oracle fixtures"]
fn matches_gatk_4620_local_assembly() {
    let directory = fixture_dir();
    let finalized_path = directory.join("t0_finalized.tsv");
    let assembled_path = directory.join("t0_assembled.tsv");
    if !finalized_path.exists() || !assembled_path.exists() {
        eprintln!("assembly fixtures absent at {}", directory.display());
        return;
    }
    let mut finalized = open(&finalized_path);
    let mut assembled = open(&assembled_path);
    let mut regions = 0_usize;
    let mut matching = 0_usize;
    let mut first_mismatch = None;
    let mut mismatch_headers = Vec::new();

    while let Some(final_header) = next_line(&mut finalized) {
        let fields: Vec<_> = final_header.split_ascii_whitespace().collect();
        assert_eq!(fields[0], "R");
        let contig = fields[1];
        let region_start: i64 = fields[2].parse().expect("region start");
        let region_end: i64 = fields[3].parse().expect("region end");
        let read_count: usize = fields[4].parse().expect("read count");
        let mut reads = Vec::with_capacity(read_count);
        for _ in 0..read_count {
            let line = next_line(&mut finalized).expect("finalized read");
            let fields: Vec<_> = line.split_ascii_whitespace().collect();
            reads.push(AsmRead {
                name: fields[0].as_bytes().to_vec(),
                flag: fields[1].parse().expect("flag"),
                start: fields[2].parse().expect("read start"),
                cigar: parse_cigar(fields[3]),
                bases: fields[4].as_bytes().to_vec(),
                quals: fields[5]
                    .split(',')
                    .map(|quality| quality.parse().expect("quality"))
                    .collect(),
            });
        }

        let expected_header = next_line(&mut assembled).expect("assembled region header");
        let expected_fields: Vec<_> = expected_header.split_ascii_whitespace().collect();
        assert_eq!(&expected_fields[..4], &fields[..4]);
        let haplotype_count: usize = expected_fields[4].parse().expect("haplotype count");
        let reference_line = next_line(&mut assembled).expect("assembled reference");
        let reference_fields: Vec<_> = reference_line.split_ascii_whitespace().collect();
        assert_eq!(reference_fields[0], "REF");
        let (_, padded_ref_start, _) = parse_location(reference_fields[1]);
        let padded_ref = reference_fields[2].as_bytes();
        let expected: Vec<String> = (0..haplotype_count)
            .map(|_| next_line(&mut assembled).expect("expected haplotype"))
            .collect();

        let actual: Vec<String> = assembly::assemble(
            contig,
            region_start,
            region_end,
            padded_ref_start,
            padded_ref,
            &reads,
        )
        .iter()
        .map(haplotype_line)
        .collect();
        regions += 1;
        if actual == expected {
            matching += 1;
        } else {
            mismatch_headers.push(expected_header.clone());
            if first_mismatch.is_none() {
                first_mismatch = Some((expected_header, expected, actual));
            }
        }
    }
    assert!(
        next_line(&mut assembled).is_none(),
        "extra assembled fixture lines"
    );
    eprintln!("matching assembly regions: {matching}/{regions}");
    eprintln!("mismatching regions: {}", mismatch_headers.join(", "));
    if let Some((header, expected, actual)) = first_mismatch {
        eprintln!("first mismatch: {header}");
        eprintln!("expected:\n{}", expected.join("\n"));
        eprintln!("actual:\n{}", actual.join("\n"));
    }
    assert_eq!(matching, regions);
    assert_eq!(regions, 1_110);
}
