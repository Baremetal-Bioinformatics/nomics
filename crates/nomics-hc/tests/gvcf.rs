use nomics_hc::annotate::{
    AnnotationCall, AnnotationGenotype, AnnotationRead, annotate, annotate_mode,
};
use nomics_hc::call::add_non_ref_likelihoods;
use nomics_hc::jrandom::JavaRandom;

#[test]
fn no_call_annotations_and_format_follow_gatk_guards() {
    use nomics_hc::genotype::NO_CALL_ALLELE;
    use nomics_hc::vcf::{Genotype, GvcfGenotype, Record, encode};
    let alleles = vec![b"T".to_vec(), b"G".to_vec(), b"<NON_REF>".to_vec()];
    let mut call = AnnotationCall {
        alleles: alleles.clone(),
        matrix_alleles: alleles,
        log10_p_error: -0.001,
        genotype: AnnotationGenotype {
            alleles: [NO_CALL_ALLELE; 2],
            pl: vec![0; 6],
            gq: 0,
            mle_ac: vec![0; 2],
            mle_af: vec![f64::NAN; 2],
        },
        contig: "chr3".into(),
        pos: 23921573,
        end: 23921573,
        reads: vec![],
    };
    let ann = annotate_mode(&call, &mut JavaRandom::gatk(), true);
    assert_eq!(
        ann.info.iter().map(|(k, _)| k.as_str()).collect::<Vec<_>>(),
        ["MLEAC", "MLEAF"]
    );
    assert!(ann.ad.is_empty());
    assert_eq!(ann.dp, -1);
    assert_eq!(ann.sb, None);
    let record = Record {
        contig: 0,
        pos: call.pos,
        id: None,
        ref_allele: "T".into(),
        alts: vec!["G".into(), "<NON_REF>".into()],
        qual: 0.01,
        info: ann.info,
        genotype: Genotype {
            alleles: vec![NO_CALL_ALLELE; 2],
            ad: ann.ad,
            dp: ann.dp,
            gq: 0,
            pl: vec![0; 6],
            gvcf: Some(GvcfGenotype::default()),
        },
    };
    let mut encoded = Vec::new();
    encode(&record, &[("chr3".into(), 200_000_000)], &mut encoded);
    assert_eq!(
        String::from_utf8(encoded).unwrap(),
        "chr3\t23921573\t.\tT\tG,<NON_REF>\t0.01\t.\tMLEAC=0,0;MLEAF=NaN,NaN\tGT:GQ:PL\t./.:0:0,0,0,0,0,0\n"
    );
    call.reads.push(AnnotationRead {
        flag: 0,
        mapq: 60,
        reverse: false,
        start: call.pos,
        cigar: vec![1 << 4],
        bases: b"T".to_vec(),
        quals: vec![30],
        likelihoods: vec![0., -10., -10.],
    });
    let ann = annotate_mode(&call, &mut JavaRandom::gatk(), true);
    assert!(ann.info.iter().any(|(k, _)| k == "DP"));
    assert!(ann.info.iter().any(|(k, _)| k == "RAW_MQandDP"));
    assert!(!ann.info.iter().any(|(k, _)| k == "ExcessHet"));
    assert!(ann.ad.is_empty());
    assert_eq!((ann.dp, ann.sb), (-1, None));
}

#[test]
fn non_ref_uses_median_of_strictly_worse_likelihoods() {
    let mut matrix = vec![
        vec![0., 0., -1.],
        vec![-2., 0., -1.],
        vec![-4., 0., -7.],
        vec![-6., 0., -3.],
    ];
    add_non_ref_likelihoods(&mut matrix);
    assert_eq!(matrix[4], [-4., 0., -5.]);
}

#[test]
fn gvcf_annotations_do_not_advance_the_shared_random_stream() {
    let call = AnnotationCall {
        alleles: vec![b"A".to_vec(), b"C".to_vec()],
        matrix_alleles: vec![b"A".to_vec(), b"C".to_vec()],
        log10_p_error: -100.,
        genotype: AnnotationGenotype {
            alleles: [1, 1],
            pl: vec![1000, 3, 0],
            gq: 3,
            mle_ac: vec![2],
            mle_af: vec![1.],
        },
        contig: "chr1".into(),
        pos: 10,
        end: 10,
        reads: vec![AnnotationRead {
            flag: 0,
            mapq: 60,
            reverse: false,
            start: 10,
            cigar: vec![1 << 4],
            bases: b"C".to_vec(),
            quals: vec![30],
            likelihoods: vec![-100., 0.],
        }],
    };
    let mut baseline = JavaRandom::gatk();
    let mut gvcf_rng = baseline.clone();
    let mut normal_rng = baseline.clone();
    let gvcf = annotate_mode(&call, &mut gvcf_rng, true);
    let normal = annotate(&call, &mut normal_rng);
    assert!(normal.info.iter().any(|(k, _)| k == "QD"));
    assert!(!gvcf.info.iter().any(|(k, _)| k == "QD"));
    assert_eq!(gvcf.sb, Some([0, 0, 1, 0]));
    let expected: Vec<_> = (0..20).map(|_| baseline.next_int(1000)).collect();
    assert_eq!(
        (0..20).map(|_| gvcf_rng.next_int(1000)).collect::<Vec<_>>(),
        expected
    );
    assert_ne!(
        (0..20)
            .map(|_| normal_rng.next_int(1000))
            .collect::<Vec<_>>(),
        expected
    );
}
