//! Exercise all library types through BAM parsing, weighting, inference, and
//! the public JSON writer. The fixed expectations include orphan encodings.
use std::path::Path;

use nomics_io::bam::{BinRule, Header, Record, RecordSet, Reference};
use nomics_quant::{Options, run};

fn record(name: &str, tid: i32, flags: u16) -> Record {
    Record {
        ref_id: tid,
        pos: 100,
        mapq: 255,
        flag: flags,
        next_ref_id: -1,
        next_pos: -1,
        tlen: 0,
        read_name: name.as_bytes().to_vec(),
        cigar: vec![75 << 4],
        seq: vec![b'A'; 75],
        qual: vec![30; 75],
        tags: vec![],
    }
}

fn fixture(dir: &Path, pairs: bool) {
    std::fs::create_dir_all(dir).unwrap();
    let mut fa = String::new();
    for i in 0..4 {
        fa.push_str(&format!(">t{i} description\n{}\n", "A".repeat(1000)));
    }
    std::fs::write(dir.join("tx.fa"), fa).unwrap();
    let mut records = Vec::new();
    if pairs {
        for (i, fw, mfw) in [
            (0, true, false),
            (1, false, true),
            (2, true, true),
            (3, false, false),
        ] {
            let mut a = record(&format!("p{i}"), i, 1 | 64 | if fw { 0 } else { 16 });
            let mut b = record(&format!("p{i}"), i, 1 | 128 | if mfw { 0 } else { 16 });
            a.next_ref_id = i;
            a.next_pos = 200;
            a.tlen = 175;
            b.pos = 200;
            b.next_ref_id = i;
            b.next_pos = 100;
            b.tlen = -175;
            records.extend([a, b]);
        }
    }
    for (name, flags) in [
        ("sf", 0),
        ("sr", 16),
        ("lf", 1 | 8 | 64),
        ("lr", 1 | 8 | 64 | 16),
        ("rf", 1 | 8 | 128),
        ("rr", 1 | 8 | 128 | 16),
    ] {
        records.push(record(name, 0, flags));
    }
    write_records(dir, records);
}

fn write_records(dir: &Path, records: Vec<Record>) {
    let header = Header {
        text: b"@HD\tVN:1.6\tGO:query\n".to_vec(),
        refs: (0..4)
            .map(|i| Reference {
                name: format!("t{i}").into_bytes(),
                len: 1000,
            })
            .collect(),
    };
    let mut bytes = Vec::new();
    for r in records {
        let mut raw = Vec::new();
        r.encode(BinRule::Sam, &mut raw).unwrap();
        bytes.extend_from_slice(&(raw.len() as i32).to_le_bytes());
        bytes.extend(raw);
    }
    RecordSet::from_parts(header, &bytes)
        .unwrap()
        .write(&dir.join("test.bam"), 1)
        .unwrap();
}

fn options(dir: &Path, library: &str, prior: f64) -> Options {
    Options {
        transcripts: dir.join("tx.fa"),
        alignments: dir.join("test.bam"),
        gene_map: None,
        output: dir.join(format!("{library}-{prior}")),
        lib_type: library.to_owned(),
        incompat_prior: prior,
        threads: 2,
    }
}

#[test]
fn every_format_and_auto_write_the_exact_expected_library_counts() {
    let dir = std::env::temp_dir().join(format!("nomics-quant-libraries-{}", std::process::id()));
    fixture(&dir, true);
    // (library, compatible fragments, concordant histogram, sense bias).
    for (format, compatible, concordant, bias) in [
        ("IU", 8, 2, 0.5),
        ("ISF", 3, 1, 0.5),
        ("ISR", 3, 1, 0.5),
        ("OU", 6, 0, 0.0),
        ("OSF", 2, 0, 0.0),
        ("OSR", 2, 0, 0.0),
        ("MU", 8, 2, 0.5),
        ("MSF", 4, 1, 0.5),
        ("MSR", 4, 1, 0.5),
        ("U", 6, 6, 0.5),
        ("SF", 1, 3, 0.5),
        ("SR", 1, 3, 0.5),
        ("A", 8, 2, 0.5),
    ] {
        let o = options(&dir, format, 0.0);
        let q = run(&o).unwrap();
        let expected = if format == "A" { "IU" } else { format };
        let mut json = serde_json::json!({
            "read_files": [], "expected_format": expected,
            "compatible_fragment_ratio": compatible as f64 / 10.0,
            "num_compatible_fragments": compatible, "num_incompatible_fragments": 10 - compatible,
            "num_assigned_fragments": compatible,
            "num_frags_with_concordant_consistent_mappings": concordant,
            "num_frags_with_inconsistent_or_orphan_mappings": 10 - concordant,
            "strand_mapping_bias": bias,
        });
        for (f, count) in [
            ("IU", 0),
            ("ISF", 1),
            ("ISR", 1),
            ("OU", 0),
            ("OSF", 0),
            ("OSR", 0),
            ("MU", 0),
            ("MSF", 1),
            ("MSR", 1),
            ("U", 0),
            ("SF", 3),
            ("SR", 3),
        ] {
            json[f] = count.into();
        }
        let actual: serde_json::Value = serde_json::from_slice(
            &std::fs::read(o.output.join("lib_format_counts.json")).unwrap(),
        )
        .unwrap();
        assert_eq!(actual, json, "{format}");
        assert_eq!(q.num_mapped, compatible);
        assert!((q.counts.iter().sum::<f64>() - compatible as f64).abs() < 1e-10);
    }
    // A positive prior keeps incompatible hits without changing strand judgment.
    let q = run(&options(&dir, "ISF", 0.1)).unwrap();
    assert_eq!(
        (q.num_mapped, q.num_compatible, q.num_incompatible),
        (10, 3, 7)
    );
    assert!((q.counts.iter().sum::<f64>() - 10.0).abs() < 1e-10);
    assert_eq!(run(&options(&dir, "ISF", -1.0)).unwrap().num_mapped, 3);
    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn auto_without_unique_proper_pairs_keeps_the_source_unresolved_behavior() {
    let dir = std::env::temp_dir().join(format!("nomics-quant-auto-empty-{}", std::process::id()));
    fixture(&dir, false);
    let o = options(&dir, "A", 0.0);
    let q = run(&o).unwrap();
    assert_eq!(q.expected_format, "A");
    assert_eq!(
        (q.num_mapped, q.num_compatible, q.num_incompatible),
        (6, 0, 0)
    );
    let actual: serde_json::Value =
        serde_json::from_slice(&std::fs::read(o.output.join("lib_format_counts.json")).unwrap())
            .unwrap();
    assert_eq!(actual["compatible_fragment_ratio"], 1.0);
    assert_eq!(actual["num_compatible_fragments"], 6);
    assert_eq!(actual["num_frags_with_concordant_consistent_mappings"], 0);
    assert_eq!(actual["num_frags_with_inconsistent_or_orphan_mappings"], 0);
    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn incompatible_prior_changes_multitranscript_posterior_weights() {
    let dir = std::env::temp_dir().join(format!("nomics-quant-soft-prior-{}", std::process::id()));
    fixture(&dir, true);
    let pair = |name: &str, tid: i32, forward: bool| {
        let mut a = record(name, tid, 1 | 64 | if forward { 0 } else { 16 });
        let mut b = record(name, tid, 1 | 128 | if forward { 16 } else { 0 });
        a.next_ref_id = tid;
        a.next_pos = 200;
        a.tlen = 175;
        b.pos = 200;
        b.next_ref_id = tid;
        b.next_pos = 100;
        b.tlen = -175;
        [a, b]
    };
    let mut records = Vec::new();
    for i in 0..100 {
        records.extend(pair(&format!("f{i}"), 0, true));
        records.extend(pair(&format!("r{i}"), 1, false));
    }
    records.extend(pair("ambiguous", 0, true));
    records.extend(pair("ambiguous", 1, false));
    write_records(&dir, records);
    let unpenalized = run(&options(&dir, "ISF", 1.0)).unwrap();
    assert!((unpenalized.counts[0] - 100.5).abs() < 1e-10);
    let penalized = run(&options(&dir, "ISF", 0.1)).unwrap();
    assert!(penalized.counts[0] > 100.9 && penalized.counts[0] < 101.0);
    assert_eq!(
        (
            penalized.num_mapped,
            penalized.num_compatible,
            penalized.num_incompatible
        ),
        (201, 101, 100)
    );
    assert!((penalized.counts.iter().sum::<f64>() - 201.0).abs() < 1e-10);
    let dropped = run(&options(&dir, "ISF", 0.0)).unwrap();
    assert_eq!(
        (dropped.num_mapped, dropped.counts[0], dropped.counts[1]),
        (101, 101.0, 0.0)
    );
    std::fs::remove_dir_all(dir).unwrap();
}
