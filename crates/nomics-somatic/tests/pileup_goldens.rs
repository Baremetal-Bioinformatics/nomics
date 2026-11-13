//! Literal byte gate against GATK 4.6.2.0 chr21 tumour and normal outputs.
use nomics_somatic::get_pileup_summaries::run;
use std::{fs, path::PathBuf};
#[test]
fn pileup_goldens() {
    let home = PathBuf::from(std::env::var_os("HOME").unwrap_or_default());
    let data = std::env::var_os("NOMICS_PILEUP_GOLDENS")
        .map(PathBuf::from)
        .unwrap_or_else(|| home.join("omics/runs/somatic/gps-chr21"));
    let reference = std::env::var_os("NOMICS_PILEUP_REFERENCE")
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            home.join(
                "omics/ref/GATK.GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta",
            )
        });
    let required = [
        "HCC1395_WES.chr21.bam",
        "HCC1395_WES.chr21.bam.bai",
        "HCC1395BL_WES.chr21.bam",
        "HCC1395BL_WES.chr21.bam.bai",
        "af-only-gnomad.chr21.vcf.gz",
        "af-only-gnomad.chr21.vcf.gz.tbi",
        "chr21.bed",
        "HCC1395_WES.chr21.pileups.table",
        "HCC1395BL_WES.chr21.pileups.table",
    ];
    if required.iter().any(|f| !data.join(f).is_file())
        || !reference.is_file()
        || !reference.with_extension("fasta.fai").is_file()
    {
        eprintln!(
            "SKIP pileup goldens: chr21 fixture or reference absent ({})",
            data.display()
        );
        return;
    }
    let output = std::env::temp_dir().join(format!(
        "nomics-pileup-goldens-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    fs::create_dir(&output).unwrap();
    for (sample, rows, processed, filters) in [
        (
            "HCC1395_WES",
            3084,
            1327905,
            [155946, 0, 0, 0, 704, 268210, 0, 0, 2589, 0, 0],
        ),
        (
            "HCC1395BL_WES",
            3115,
            1695291,
            [209062, 0, 0, 0, 2306, 434260, 0, 0, 1992, 0, 0],
        ),
    ] {
        let table = output.join(format!("{sample}.table"));
        let stats = run(
            &data.join(format!("{sample}.chr21.bam")),
            &data.join("af-only-gnomad.chr21.vcf.gz"),
            &reference,
            &data.join("chr21.bed"),
            &table,
        )
        .unwrap();
        assert_eq!(stats.rows, rows);
        assert_eq!(stats.processed, processed);
        assert_eq!(stats.filtered, filters);
        assert!(
            fs::read(&table).unwrap()
                == fs::read(data.join(format!("{sample}.chr21.pileups.table"))).unwrap(),
            "{sample}: table bytes differ"
        );
    }
    fs::remove_dir_all(output).unwrap();
}

#[test]
fn sparse_intervals_feature_priority_and_empty_loci_from_java() {
    use nomics_io::{
        bam::{self, BinRule, Header, Record, Reference, Value},
        bgzf::{BgzfWriter, Profile, Stored},
        htsjdk_index,
    };
    use std::{fs::File, io::Write};
    let dir = std::env::temp_dir().join(format!("nomics-pileup-small-{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    let header = Header {
        text: b"@HD\tSO:coordinate\n@RG\tID:z\tSM:Z\n@RG\tID:a\tSM:A\n".to_vec(),
        refs: vec![
            Reference {
                name: b"c1".to_vec(),
                len: 10,
            },
            Reference {
                name: b"c2".to_vec(),
                len: 10,
            },
        ],
    };
    let bam = dir.join("in.bam");
    let mut w = BgzfWriter::new(File::create(&bam).unwrap(), Stored, Profile::HTSJDK);
    header.write(&mut w).unwrap();
    let mut r = Record {
        ref_id: 0,
        pos: 0,
        mapq: 60,
        flag: 0,
        next_ref_id: -1,
        next_pos: -1,
        tlen: 0,
        read_name: b"r".to_vec(),
        cigar: vec![1 << 4, 2 << 4 | 2, 7 << 4],
        seq: b"ANGAAAAA".to_vec(),
        qual: vec![0; 8],
        tags: vec![(*b"RG", Value::String(b"a".to_vec()))],
    };
    let mut raw = Vec::new();
    r.encode(BinRule::Htsjdk, &mut raw).unwrap();
    bam::write_raw_record(&mut w, &raw).unwrap();
    // A supplementary read cannot turn an uncovered c2 site into a row.
    r.ref_id = 1;
    r.flag = 0x800;
    r.cigar = vec![8 << 4];
    raw.clear();
    r.encode(BinRule::Htsjdk, &mut raw).unwrap();
    bam::write_raw_record(&mut w, &raw).unwrap();
    w.finish().unwrap();
    htsjdk_index::write_bai_htsjdk(&bam, &dir.join("in.bam.bai")).unwrap();
    let vcf = dir.join("in.vcf.gz");
    let mut w = BgzfWriter::new(File::create(&vcf).unwrap(), Stored, Profile::HTSJDK);
    writeln!(w,"##fileformat=VCFv4.2\n##contig=<ID=c1,length=10>\n##contig=<ID=c2,length=10>\n##INFO=<ID=AF,Number=A,Type=Float,Description=\"AF\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();
    for line in [
        "c1\t1\t.\tAA\tA\t.\t.\tAF=0.1",
        "c1\t2\t.\tA\tC\t.\t.\tAF=0.1",
        "c1\t3\t.\tA\tC\t.\t.\tAF=0.1",
        "c1\t4\t.\tA\tC\t.\t.\tAF=0.1",
        "c1\t5\t.\tA\tC\t.\tFAIL\tAF=0.1",
        "c1\t6\t.\tA\tC\t.\t.\tAF=0.01",
        "c1\t7\t.\tA\tC,G\t.\t.\tAF=0.1,0.1",
        "c1\t8\t.\tA\tC\t.\t.\t.",
        "c1\t10\t.\tA\tC\t.\t.\tAF=0.1",
        "c2\t1\t.\tA\tC\t.\t.\tAF=.",
    ] {
        writeln!(w, "{line}").unwrap();
    }
    w.finish().unwrap();
    htsjdk_index::write_tbi_htsjdk_vcf(&vcf, &dir.join("in.vcf.gz.tbi")).unwrap();
    let reference = dir.join("ref.fasta");
    fs::write(&reference, b">c1\nAAAAAAAAAA\n>c2\nAAAAAAAAAA\n").unwrap();
    fs::write(
        dir.join("ref.fasta.fai"),
        "c1\t10\t4\t10\t11\nc2\t10\t19\t10\t11\n",
    )
    .unwrap();
    let bed = dir.join("in.bed");
    fs::write(&bed, "c2\t0\t10\nc1\t5\t10\nc1\t0\t3\nc1\t2\t5\n").unwrap();
    let out = dir.join("out.table");
    let stats = run(&bam, &vcf, &reference, &bed, &out).unwrap();
    // Derived from GetPileupSummaries.apply and LocusIteratorByState:
    // earlier indel hides SNP2; D and N-base rows have zero counts; G is
    // other-alt only; strict AF .01 excluded; absent AF skipped; empty c2 skipped.
    assert_eq!(
        fs::read_to_string(out).unwrap(),
        "#<METADATA>SAMPLE=A\ncontig\tposition\tref_count\talt_count\tother_alt_count\tallele_frequency\nc1\t3\t0\t0\t0\t0.1\nc1\t4\t0\t0\t0\t0.1\nc1\t5\t0\t0\t1\t0.1\nc1\t10\t1\t0\t0\t0.1\n"
    );
    assert_eq!(stats.filtered[4], 1);
    // No intervals/read coverage writes metadata plus header, without rows.
    fs::write(&bed, "c2\t9\t10\n").unwrap();
    assert_eq!(
        run(&bam, &vcf, &reference, &bed, &dir.join("empty.table"))
            .unwrap()
            .rows,
        0
    );
    fs::remove_dir_all(dir).unwrap();
}
