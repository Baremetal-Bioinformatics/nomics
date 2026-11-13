//! Optional differential test using a separately built, unmodified Subread 2.1.1
//! executable as a test oracle. The production library never invokes it.
use std::path::PathBuf;
use std::process::Command;

use nomics_count::{Annotation, CountOptions, count_bam, write_summary, write_table};
use nomics_io::bam::{
    BinRule, Header, IntType, Record, RecordSet, Reference, Value, write_raw_record,
};

#[test]
fn synthetic_source_oracle() {
    let Some(executable) = std::env::var_os("NOMICS_FEATURECOUNTS_ORACLE") else {
        eprintln!("SKIP: set NOMICS_FEATURECOUNTS_ORACLE to a Subread 2.1.1 executable");
        return;
    };
    let executable = std::fs::canonicalize(executable).unwrap();
    let version = Command::new(&executable).arg("-v").output().unwrap();
    let version = format!(
        "{}{}",
        String::from_utf8_lossy(&version.stdout),
        String::from_utf8_lossy(&version.stderr)
    );
    assert!(
        version.contains("v2.1.1"),
        "requires Subread 2.1.1, got {version}"
    );
    let temp = std::env::temp_dir().join(format!("nomics-count-oracle-{}", std::process::id()));
    std::fs::create_dir_all(&temp).unwrap();
    // Deliberately unsorted, equal-start exons, duplicates, two chromosomes,
    // unknown strand, overlapping genes and attributes with quoted semicolons.
    let mut gtf = String::new();
    for i in [7, 3, 9, 0, 12, 8, 1, 10, 6, 2, 11, 4, 5, 13, 14, 15] {
        let chr = if i < 12 { "I" } else { "II" };
        let strand = ["+", "-", "."][i % 3];
        let start = 100 + i % 4 * 10;
        gtf.push_str(&format!("{chr}\ttest\texon\t{start}\t{}\t.\t{strand}\t.\tother \"a;b\"; gene_id \"g{}\"; gene_biotype \"b{}\";\n", start + 30 + i, i % 4, i % 2));
    }
    // The output Length field is a signed C int, even for meta-features whose
    // union across chromosomes exceeds 2^31 bases.
    for chr in ["L1", "L2", "L3"] {
        gtf.push_str(&format!("{chr}\ttest\texon\t1\t1000000000\t.\t+\t.\tgene_id \"large\"; gene_biotype \"large\";\n"));
    }
    let gtf_path = temp.join("genome.gtf");
    std::fs::write(&gtf_path, gtf).unwrap();
    let mut data = Vec::new();
    let tags = [
        None,
        Some((IntType::U8, 1)),
        Some((IntType::I8, -1)),
        Some((IntType::I16, -1)),
        Some((IntType::I32, -1)),
        Some((IntType::U32, u32::MAX as i64)),
        Some((IntType::U8, 2)),
    ];
    let cigars = [
        vec![10 << 4],
        vec![10 << 4, 60 << 4 | 3, 10 << 4],
        vec![10 << 4, 4 << 4 | 2, 10 << 4],
        vec![5 << 4, 2 << 4 | 1, 5 << 4],
        vec![2 << 4 | 1, 10 << 4],
        vec![2 << 4 | 4, 5 << 4 | 7, 5 << 4 | 8, 2 << 4 | 5],
    ];
    let mut id = 0;
    for flag in [
        0,
        4,
        1,
        1 | 4,
        1 | 8,
        1 | 16 | 64,
        1 | 128,
        1 | 16 | 128,
        1 | 64 | 256 | 1024,
        1 | 16 | 2048,
    ] {
        for pos in [98, 109, 129, 149, 175] {
            for (ty, cigar) in tags.iter().cycle().zip(cigars.iter()).take(cigars.len()) {
                let length: usize = cigar
                    .iter()
                    .filter(|&&c| matches!(c & 15, 0 | 1 | 4 | 7 | 8))
                    .map(|&c| (c >> 4) as usize)
                    .sum();
                let record = Record {
                    ref_id: id % 2,
                    pos,
                    mapq: 0,
                    flag,
                    next_ref_id: (id + 1) % 2,
                    next_pos: pos + 100,
                    tlen: 100000,
                    read_name: format!("read{id}").into_bytes(),
                    cigar: cigar.clone(),
                    seq: vec![b'A'; length],
                    qual: vec![30; length],
                    tags: ty
                        .iter()
                        .map(|&(kind, n)| (*b"NH", Value::Int(kind, n)))
                        .collect(),
                };
                let mut raw = Vec::new();
                record.encode(BinRule::Sam, &mut raw).unwrap();
                write_raw_record(&mut data, &raw).unwrap();
                id += 1;
            }
        }
    }
    let header = Header {
        text: b"@HD\tVN:1.6\tSO:unsorted\n".to_vec(),
        refs: [b"I".to_vec(), b"chrII".to_vec()]
            .map(|name| Reference { name, len: 1000 })
            .to_vec(),
    };
    let bam = temp.join("test.bam");
    RecordSet::from_parts(header, &data)
        .unwrap()
        .write(&bam, 2)
        .unwrap();
    for group in ["gene_id", "gene_biotype"] {
        let annotation = Annotation::read(&gtf_path, group, "exon").unwrap();
        for strand in 0..=2 {
            let output = temp.join(format!("{group}-{strand}.tsv"));
            let mut args: Vec<String> = ["-B", "-C", "-p", "-g", group, "-s"]
                .map(str::to_owned)
                .to_vec();
            args.extend([
                strand.to_string(),
                "-T".into(),
                "6".into(),
                "-a".into(),
                gtf_path.display().to_string(),
                "-o".into(),
                output.display().to_string(),
                bam.display().to_string(),
            ]);
            let result = Command::new(&executable).args(&args).output().unwrap();
            assert!(
                result.status.success(),
                "{}",
                String::from_utf8_lossy(&result.stderr)
            );
            let counts = [count_bam(
                &annotation,
                &bam,
                &CountOptions {
                    strand,
                    ..CountOptions::default()
                },
            )
            .unwrap()];
            let mut command = vec![executable.display().to_string()];
            command.extend(args);
            let inputs = [bam.display().to_string()];
            let mut ours = Vec::new();
            write_table(&mut ours, &annotation, &inputs, &counts, &command).unwrap();
            let golden = std::fs::read(&output).unwrap();
            assert_eq!(
                String::from_utf8_lossy(&ours),
                String::from_utf8_lossy(&golden),
                "{group} -s {strand} table"
            );
            ours.clear();
            write_summary(&mut ours, &inputs, &counts).unwrap();
            assert_eq!(
                String::from_utf8_lossy(&ours),
                std::fs::read_to_string(PathBuf::from(format!("{}.summary", output.display())))
                    .unwrap(),
                "{group} -s {strand} summary"
            );
        }
    }
    std::fs::remove_dir_all(temp).unwrap();
}
