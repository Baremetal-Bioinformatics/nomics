use nomics_io::bam::{BinRule, Header, IntType, Record, RecordSet, Reference, Value};
use nomics_solo::{Config, count};
use std::path::Path;

fn record(name: &str, flag: u16, nh: u32, gx: &str, cb: &str, umi: &str) -> Record {
    Record {
        ref_id: if flag & 4 == 0 { 0 } else { -1 },
        pos: if flag & 4 == 0 { 0 } else { -1 },
        mapq: 255,
        flag,
        next_ref_id: -1,
        next_pos: -1,
        tlen: 0,
        read_name: name.as_bytes().to_vec(),
        cigar: if flag & 4 == 0 { vec![160] } else { vec![] },
        seq: vec![b'A'; 10],
        qual: vec![35; 10],
        tags: vec![
            (*b"NH", Value::Int(IntType::I32, i64::from(nh))),
            (*b"CR", Value::String(cb.as_bytes().to_vec())),
            (*b"CY", Value::String(vec![73; 16])),
            (*b"UR", Value::String(umi.as_bytes().to_vec())),
            (*b"UY", Value::String(vec![73; 12])),
            (*b"gx", Value::String(gx.as_bytes().to_vec())),
            // Intentionally incorrect corrected tags must never influence counting.
            (*b"CB", Value::String(b"TTTTTTTTTTTTTTTT".to_vec())),
            (*b"UB", Value::String(b"TTTTTTTTTTTT".to_vec())),
        ],
    }
}
fn write_bam(path: &Path, records: &[Record]) {
    let header = Header {
        text: b"@HD\tVN:1.6\n".to_vec(),
        refs: vec![Reference {
            name: b"chr1".to_vec(),
            len: 1000,
        }],
    };
    let mut payload = Vec::new();
    for r in records {
        let mut raw = Vec::new();
        r.encode(BinRule::Sam, &mut raw).unwrap();
        payload.extend_from_slice(&(raw.len() as u32).to_le_bytes());
        payload.extend_from_slice(&raw);
    }
    RecordSet::from_parts(header, &payload)
        .unwrap()
        .write(path, 1)
        .unwrap();
}

#[test]
fn multimapper_union_and_corrected_tags_are_not_inputs() {
    let dir = std::env::temp_dir().join(format!("solo-synthetic-{}", std::process::id()));
    std::fs::create_dir(&dir).unwrap();
    let cb = "AAAAAAAAAAAAAAAA";
    // AAA...AC -> AAA...CC -> AAA...CG: a connected component across Hamming-1 edges.
    let records = vec![
        record("one", 0, 1, "g1", cb, "AAAAAAAAAAAC"),
        record("two", 0, 1, "g1", cb, "AAAAAAAAAACC"),
        record("three", 0, 1, "g1", cb, "AAAAAAAAAACG"),
        // Unique gene even though the primary alignment has no feature.
        record("multi", 0, 2, "-", cb, "CCCCCCCCCCCA"),
        record("ambiguous", 0, 2, "g1", cb, "GGGGGGGGGGGA"),
        record("multi", 256, 2, "g2", cb, "CCCCCCCCCCCA"),
        record("ambiguous", 256, 2, "g2", cb, "GGGGGGGGGGGA"),
        record("badumi", 0, 1, "g1", "NNNNNNNNNNNNNNNN", "AAAAAAAAAAAA"),
        record("unmapped", 4, 0, "-", cb, "GGGGGGGGGGGA"),
        // Single 1MM correction, with a valid prior in the completed stream.
        record("corrected", 0, 1, "g2", "CAAAAAAAAAAAAAAA", "CCCCCCCCCCCA"),
    ];
    write_bam(&dir.join("input.bam"), &records);
    std::fs::write(dir.join("wl"), format!("{cb}\n")).unwrap();
    std::fs::write(
        dir.join("features"),
        "g1\tGene1\tGene Expression\ng2\tGene2\tGene Expression\n",
    )
    .unwrap();
    let config = Config {
        bam: dir.join("input.bam"),
        whitelist: dir.join("wl"),
        features: dir.join("features"),
        output: dir.join("Solo.out"),
        threads: 1,
    };
    let report = count(&config).unwrap();
    assert_eq!(
        (
            report.reads,
            report.counted_reads,
            report.umis,
            report.cells
        ),
        (8, 5, 2, 1)
    );
    assert_eq!(
        std::fs::read_to_string(config.output.join("Gene/raw/matrix.mtx")).unwrap(),
        "%%MatrixMarket matrix coordinate integer general\n%\n2 1 2\n1 1 1\n2 1 1\n"
    );
    let stats = std::fs::read_to_string(config.output.join("Gene/Features.stats")).unwrap();
    let counts: Vec<u64> = stats
        .lines()
        .map(|l| l.split_whitespace().nth(1).unwrap().parse().unwrap())
        .collect();
    assert_eq!(counts, vec![1, 0, 1, 1, 0, 0, 5, 4, 5, 1, 2]);
    // A second run with absent CB/UB tags must match every requested file.
    let stripped: Vec<_> = records
        .into_iter()
        .map(|mut r| {
            r.tags.retain(|(tag, _)| tag != b"CB" && tag != b"UB");
            r
        })
        .collect();
    write_bam(&dir.join("stripped.bam"), &stripped);
    let second = Config {
        bam: dir.join("stripped.bam"),
        output: dir.join("without-tags"),
        ..config
    };
    count(&second).unwrap();
    for file in [
        "Barcodes.stats",
        "Gene/Features.stats",
        "Gene/Summary.csv",
        "Gene/UMIperCellSorted.txt",
        "Gene/raw/matrix.mtx",
        "Gene/raw/barcodes.tsv",
        "Gene/raw/features.tsv",
        "Gene/filtered/matrix.mtx",
        "Gene/filtered/barcodes.tsv",
        "Gene/filtered/features.tsv",
    ] {
        assert_eq!(
            std::fs::read(dir.join("Solo.out").join(file)).unwrap(),
            std::fs::read(second.output.join(file)).unwrap(),
            "{file}"
        );
    }
    std::fs::remove_dir_all(dir).unwrap();
}
