use nomics_hc::assembly::{self, AsmRead, Haplotype};
use nomics_hc::javafmt;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};

/// `~/omics-fixtures/asm-dump`: every 5th active region of the exome run
/// (t0_finalized.tsv, t0_assembled.tsv), written by the 02832d1 code.
const DEFAULT_FIXTURES: &str = "omics-fixtures/asm-dump";

struct Region {
    label: String,
    contig: String,
    start: i64,
    end: i64,
    padded_ref_start: i64,
    padded_ref: Vec<u8>,
    reads: Vec<AsmRead>,
}

fn fixture_dir() -> PathBuf {
    std::env::var_os("NOMICS_HC_FIXTURES")
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            PathBuf::from(std::env::var_os("HOME").unwrap_or_default()).join(DEFAULT_FIXTURES)
        })
}

fn open(path: &Path) -> BufReader<File> {
    BufReader::with_capacity(
        1024 * 1024,
        File::open(path).expect("open assembly fixture"),
    )
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

fn parse_location(location: &str) -> (&str, i64, i64) {
    let (contig, positions) = location.rsplit_once(':').expect("location colon");
    let (start, end) = positions.split_once('-').expect("location dash");
    (
        contig,
        start.parse().expect("location start"),
        end.parse().expect("location end"),
    )
}

fn load_regions(directory: &Path) -> Vec<Region> {
    let mut finalized = open(&directory.join("t0_finalized.tsv"));
    let mut assembled = open(&directory.join("t0_assembled.tsv"));
    let mut regions = Vec::new();
    while let Some(final_header) = next_line(&mut finalized) {
        let fields: Vec<_> = final_header.split_ascii_whitespace().collect();
        assert_eq!(fields[0], "R");
        let contig = fields[1].to_owned();
        let start = fields[2].parse().expect("region start");
        let end = fields[3].parse().expect("region end");
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
        let assembled_header = next_line(&mut assembled).expect("assembled region header");
        let assembled_fields: Vec<_> = assembled_header.split_ascii_whitespace().collect();
        assert_eq!(&assembled_fields[..4], &fields[..4]);
        let haplotype_count: usize = assembled_fields[4].parse().expect("haplotype count");
        let reference_line = next_line(&mut assembled).expect("assembled reference");
        let reference_fields: Vec<_> = reference_line.split_ascii_whitespace().collect();
        assert_eq!(reference_fields[0], "REF");
        let (_, padded_ref_start, _) = parse_location(reference_fields[1]);
        for _ in 0..haplotype_count {
            next_line(&mut assembled).expect("assembled haplotype");
        }
        regions.push(Region {
            label: format!("{contig}:{start}-{end}"),
            contig,
            start,
            end,
            padded_ref_start,
            padded_ref: reference_fields[2].as_bytes().to_vec(),
            reads,
        });
    }
    assert!(
        next_line(&mut assembled).is_none(),
        "extra assembled fixture lines"
    );
    regions
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

fn hash_bytes(hash: &mut u64, bytes: &[u8]) {
    for byte in bytes {
        *hash ^= u64::from(*byte);
        *hash = hash.wrapping_mul(0x100000001b3);
    }
}

fn main() {
    let repeats: usize = std::env::args()
        .nth(1)
        .map(|arg| arg.parse().expect("repeat count must be an integer"))
        .unwrap_or(5);
    assert!(repeats != 0, "repeat count must be positive");
    let directory = fixture_dir();
    let regions = load_regions(&directory);
    assert!(!regions.is_empty());
    eprintln!(
        "loaded {} regions from {}",
        regions.len(),
        directory.display()
    );

    let mut totals = vec![Duration::ZERO; regions.len()];
    let mut expected_checksum = None;
    for pass in 0..repeats {
        let pass_start = Instant::now();
        let mut checksum = 0xcbf29ce484222325_u64;
        let mut haplotypes = 0_usize;
        for (index, region) in regions.iter().enumerate() {
            let start = Instant::now();
            let output = assembly::assemble(
                &region.contig,
                region.start,
                region.end,
                region.padded_ref_start,
                &region.padded_ref,
                &region.reads,
            );
            totals[index] += start.elapsed();
            haplotypes += output.len();
            for haplotype in &output {
                hash_bytes(&mut checksum, haplotype_line(haplotype).as_bytes());
                hash_bytes(&mut checksum, b"\n");
            }
        }
        let elapsed = pass_start.elapsed();
        if let Some(expected) = expected_checksum {
            assert_eq!(checksum, expected, "checksum changed between passes");
        } else {
            expected_checksum = Some(checksum);
        }
        println!(
            "pass {}/{}: {:.3}s, {} haplotypes, checksum fnv1a64:{checksum:016x}",
            pass + 1,
            repeats,
            elapsed.as_secs_f64(),
            haplotypes,
        );
    }

    let mut slowest: Vec<_> = totals.iter().enumerate().collect();
    slowest.sort_by(|(_, left), (_, right)| right.cmp(left));
    println!("slowest regions (mean over {repeats} passes):");
    for (index, elapsed) in slowest.into_iter().take(20) {
        println!(
            "  {:8.3} ms  {}  {} reads",
            elapsed.as_secs_f64() * 1000.0 / repeats as f64,
            regions[index].label,
            regions[index].reads.len(),
        );
    }
}
