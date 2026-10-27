//! Hand-made BAM -> CRAM -> decode with samtools: the records and header
//! must come back unchanged (bar the `M5:` we add), and our `.crai` must
//! equal `samtools index`'s. Skipped when samtools is not on PATH.
//!
//! The records cover what the writer has to get right: substitutions against
//! A/C/G/T/N and IUPAC reference bases, lower-case reference, every CIGAR op
//! a CRAM can carry, reads running off the end of a contig, attached and
//! detached mates (including a start tie that htslib and htsjdk resolve
//! differently), secondary/supplementary, placed and unplaced unmapped reads,
//! missing qualities, every BAM tag type (a `Z` with a tab, `B` arrays), and
//! a mapped read without MD/NM. Four records per slice, so pairs also span
//! containers.

use std::path::{Path, PathBuf};
use std::process::Command;

use nomics_io::bam::{ArrayValue, Header, IntType, Record, RecordSet, Reference, Value};
use nomics_io::cram::{CramOptions, crai_bytes, write_cram};
use nomics_io::fasta::Fasta;

fn samtools_ok() -> bool {
    Command::new("samtools")
        .arg("--version")
        .output()
        .is_ok_and(|o| o.status.success())
}

fn contig(seed: u64, n: usize) -> Vec<u8> {
    let mut x = seed;
    (0..n)
        .map(|_| {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            b"ACGT"[(x >> 60) as usize & 3]
        })
        .collect()
}

fn write_fasta(dir: &Path, contigs: &[(&str, Vec<u8>)]) -> PathBuf {
    let fa = dir.join("ref.fa");
    let (mut text, mut fai) = (Vec::new(), String::new());
    for (name, seq) in contigs {
        text.extend_from_slice(format!(">{name} test\n").as_bytes());
        let offset = text.len();
        for line in seq.chunks(60) {
            text.extend_from_slice(line);
            text.push(b'\n');
        }
        fai.push_str(&format!("{name}\t{}\t{offset}\t60\t61\n", seq.len()));
    }
    std::fs::write(&fa, text).unwrap();
    std::fs::write(dir.join("ref.fa.fai"), fai).unwrap();
    fa
}

fn parse_cigar(s: &str) -> Vec<u32> {
    let mut out = Vec::new();
    let mut n = 0u32;
    for c in s.bytes() {
        if c.is_ascii_digit() {
            n = n * 10 + (c - b'0') as u32;
        } else {
            let op = b"MIDNSHP=X".iter().position(|&o| o == c).unwrap() as u32;
            out.push(n << 4 | op);
            n = 0;
        }
    }
    out
}

/// A read whose aligned bases copy `reference` (upper-cased; past the end:
/// `A`), inserted/clipped bases are `T`, then `edits` (read offset, base).
fn read_seq(reference: &[u8], pos: i32, cigar: &[u32], edits: &[(usize, u8)]) -> Vec<u8> {
    let mut seq = Vec::new();
    let mut r = pos.max(0) as usize;
    for &c in cigar {
        let (op, len) = (c & 0xf, (c >> 4) as usize);
        match op {
            0 | 7 | 8 => {
                for k in 0..len {
                    seq.push(
                        reference
                            .get(r + k)
                            .map_or(b'A', |b| b.to_ascii_uppercase()),
                    );
                }
                r += len;
            }
            1 | 4 => seq.extend(std::iter::repeat_n(b'T', len)),
            2 | 3 => r += len,
            _ => {}
        }
    }
    for &(i, b) in edits {
        seq[i] = b;
    }
    seq
}

struct Rec {
    name: &'static str,
    flag: u16,
    ref_id: i32,
    pos: i32,
    mapq: u8,
    cigar: &'static str,
    edits: Vec<(usize, u8)>,
    next: (i32, i32, i32),
    tags: Vec<([u8; 2], Value)>,
    no_qual: bool,
    unmapped_len: usize,
}

fn rec(name: &'static str, flag: u16, ref_id: i32, pos: i32, cigar: &'static str) -> Rec {
    Rec {
        name,
        flag,
        ref_id,
        pos,
        mapq: if flag & 4 != 0 { 0 } else { 60 },
        cigar,
        edits: vec![],
        next: (-1, -1, 0),
        tags: vec![],
        no_qual: false,
        unmapped_len: 0,
    }
}

fn std_tags() -> Vec<([u8; 2], Value)> {
    vec![
        (*b"MD", Value::String(b"30".to_vec())),
        (*b"RG", Value::String(b"rg1".to_vec())),
        (*b"NM", Value::Int(IntType::U8, 0)),
        (*b"AS", Value::Int(IntType::U8, 30)),
    ]
}

fn records(contigs: &[(&str, Vec<u8>)]) -> Vec<Record> {
    let mut v: Vec<Rec> = Vec::new();
    let mate = |mut r: Rec, next_ref: i32, next_pos: i32, tlen: i32| {
        r.next = (next_ref, next_pos, tlen);
        r.tags = std_tags();
        r
    };
    // Contig c1 (300 bp: IUPAC R at 50, N run 100..106, lower case 200..220).
    v.push(mate(rec("pairA", 99, 0, 9, "30M"), 0, 59, 80));
    let mut t = mate(rec("pairA", 147, 0, 59, "30M"), 0, 9, -80);
    t.edits = vec![(3, b'G'), (4, b'N'), (29, b'C')];
    v.push(t);
    // A start tie: htslib gives the shorter-ending read +TLEN, htsjdk 5 the
    // first-of-pair read. Stored detached either way.
    v.push(mate(rec("tie", 163, 0, 20, "20M10S"), 0, 20, -30));
    v.push(mate(rec("tie", 83, 0, 20, "30M"), 0, 20, 30));
    let mut t = rec("iupac", 0, 0, 45, "3S12M");
    t.edits = vec![(8, b'A'), (9, b'C'), (14, b'R')]; // ref R at 50 vs A
    t.tags = std_tags();
    v.push(t);
    let mut t = rec("ops", 16, 0, 60, "3H2S10M2I5M3D10M4S5H");
    t.edits = vec![(12, b'G')];
    t.tags = std_tags();
    v.push(t);
    let mut t = rec("skip_pad", 0, 0, 70, "10M50N5M1P2I5M");
    t.tags = std_tags();
    v.push(t);
    let mut t = rec("nrun", 0, 0, 95, "20M");
    t.edits = vec![(6, b'N'), (7, b'A'), (8, b'n'.to_ascii_uppercase())];
    t.tags = std_tags();
    v.push(t);
    v.push(mate(rec("other_contig", 97, 0, 120, "25M"), 1, 9, 0));
    v.push(mate(rec("placed", 73, 0, 150, "30M"), 0, 150, 0));
    let mut t = mate(rec("placed", 133, 0, 150, ""), 0, 150, 0);
    t.unmapped_len = 30;
    t.tags = vec![(*b"RG", Value::String(b"rg1".to_vec()))];
    v.push(t);
    let mut t = rec("lower", 0, 0, 195, "30M");
    t.tags = std_tags();
    t.edits = vec![(10, b'C'), (11, b'G')];
    v.push(t);
    let mut t = rec("alltags", 0, 0, 210, "20M");
    t.tags = vec![
        (*b"XA", Value::Char(b'q')),
        (*b"Xc", Value::Int(IntType::I8, -5)),
        (*b"XC", Value::Int(IntType::U8, 250)),
        (*b"Xs", Value::Int(IntType::I16, -300)),
        (*b"XS", Value::Int(IntType::U16, 60000)),
        (*b"Xi", Value::Int(IntType::I32, -70000)),
        (*b"XI", Value::Int(IntType::U32, 4_000_000_000)),
        (*b"Xf", Value::Float(1.5f32.to_bits())),
        (*b"XZ", Value::String(b"tab\there".to_vec())),
        (*b"XH", Value::Hex(b"1AE301".to_vec())),
        (
            *b"XB",
            Value::Array(ArrayValue::Int(IntType::I8, vec![-1, 0, 1])),
        ),
        (
            *b"Xb",
            Value::Array(ArrayValue::Int(IntType::U16, vec![1, 65535])),
        ),
        (
            *b"Xg",
            Value::Array(ArrayValue::Float(vec![0.25f32.to_bits()])),
        ),
        (*b"Xe", Value::Array(ArrayValue::Int(IntType::I32, vec![]))),
    ];
    v.push(t);
    v.push(rec("notags", 0, 0, 215, "20M"));
    let mut t = rec("noqual", 16, 0, 220, "20M");
    t.no_qual = true;
    t.tags = std_tags();
    v.push(t);
    v.push(mate(rec("supp", 2145, 0, 230, "15M15S"), 1, 49, 0));
    v.push(mate(rec("supp", 337, 0, 240, "30M"), 1, 49, 0));
    let mut t = rec("offend", 0, 0, 285, "25M");
    t.tags = std_tags();
    v.push(t);
    // Contig c2 (120 bp).
    v.push(mate(rec("other_contig", 145, 1, 9, "25M"), 0, 120, 0));
    v.push(mate(rec("supp", 81, 1, 49, "30M"), 0, 230, 0));
    for i in 0..5 {
        v.push(mate(rec("dense", 99, 1, 60 + i, "20M"), 1, 70 + i, 30));
        v.push(mate(rec("dense", 147, 1, 70 + i, "20M"), 1, 60 + i, -30));
    }
    // Unplaced unmapped tail.
    for (name, flag, no_qual) in [("u1", 77, false), ("u1", 141, true), ("u2", 4, false)] {
        let mut t = rec(name, flag, -1, -1, "");
        t.unmapped_len = 12;
        t.no_qual = no_qual;
        t.tags = vec![(*b"RG", Value::String(b"rg1".to_vec()))];
        v.push(t);
    }

    v.into_iter()
        .enumerate()
        .map(|(i, r)| {
            let cigar = parse_cigar(r.cigar);
            let seq = if r.flag & 4 != 0 {
                (0..r.unmapped_len).map(|k| b"ACGTN"[(k + i) % 5]).collect()
            } else {
                read_seq(&contigs[r.ref_id as usize].1, r.pos, &cigar, &r.edits)
            };
            let qual = if r.no_qual {
                vec![0xff; seq.len()]
            } else {
                (0..seq.len())
                    .map(|k| 2 + ((k * 7 + i) % 39) as u8)
                    .collect()
            };
            // Distinct names per dense pair.
            let name = if r.name == "dense" {
                format!("dense{}", (r.pos - 60) % 10).into_bytes()
            } else {
                r.name.as_bytes().to_vec()
            };
            Record {
                ref_id: r.ref_id,
                pos: r.pos,
                mapq: r.mapq,
                flag: r.flag,
                next_ref_id: r.next.0,
                next_pos: r.next.1,
                tlen: r.next.2,
                read_name: name,
                cigar,
                seq,
                qual,
                tags: r.tags,
            }
        })
        .collect()
}

fn run(args: &[&str]) -> Vec<u8> {
    let out = Command::new("samtools").args(args).output().unwrap();
    assert!(
        out.status.success(),
        "samtools {args:?}: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    out.stdout
}

#[test]
fn bam_to_cram_round_trips_through_samtools() {
    if !samtools_ok() {
        eprintln!("samtools not found; skipping");
        return;
    }
    let dir = std::env::temp_dir().join(format!("nomics-cram-rt-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let mut c1 = contig(1, 300);
    c1[50] = b'R';
    for b in &mut c1[100..106] {
        *b = b'N';
    }
    for b in &mut c1[200..220] {
        *b = b.to_ascii_lowercase();
    }
    let contigs = vec![("c1", c1), ("c2", contig(2, 120))];
    let fa = write_fasta(&dir, &contigs);

    let text = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c1\tLN:300\n@SQ\tSN:c2\tLN:120\n\
                @RG\tID:rg1\tSM:s\n@PG\tID:test\tPN:test\n";
    let header = Header {
        text: text.as_bytes().to_vec(),
        refs: vec![
            Reference {
                name: b"c1".to_vec(),
                len: 300,
            },
            Reference {
                name: b"c2".to_vec(),
                len: 120,
            },
        ],
    };
    let fa_s = fa.to_str().unwrap();
    let fasta = Fasta::load(&fa, 2).unwrap();
    let encode = |recs: &[Record]| {
        let mut raw = Vec::new();
        let mut body = Vec::new();
        for r in recs {
            body.clear();
            r.encode(nomics_io::bam::BinRule::Sam, &mut body).unwrap();
            raw.extend_from_slice(&(body.len() as i32).to_le_bytes());
            raw.extend_from_slice(&body);
        }
        RecordSet::from_parts(header.clone(), &raw).unwrap()
    };

    // Unsorted input (AP stored absolute where a slice goes backwards):
    // records only.
    let unsorted = encode(&records(&contigs));
    let ubam = dir.join("unsorted.bam");
    unsorted.write(&ubam, 1).unwrap();
    let (ucram, _) = write_cram(
        Vec::new(),
        &header,
        &fasta,
        &CramOptions::default(),
        |sink| {
            for i in 0..unsorted.len() {
                sink.push(unsorted.record(i))?;
            }
            Ok(())
        },
    )
    .unwrap();
    let ucram_path = dir.join("unsorted.cram");
    std::fs::write(&ucram_path, ucram).unwrap();
    assert_eq!(
        run(&["view", "--no-PG", "-T", fa_s, ucram_path.to_str().unwrap()]),
        run(&["view", "--no-PG", ubam.to_str().unwrap()])
    );

    let mut sorted = records(&contigs);
    sorted.sort_by_key(|r| (r.ref_id < 0, r.ref_id, r.pos));
    let set = encode(&sorted);
    let bam = dir.join("in.bam");
    set.write(&bam, 1).unwrap();

    let opts = CramOptions {
        records_per_slice: 4,
        threads: 3,
        ..Default::default()
    };
    let cram = dir.join("out.cram");
    let (bytes, summary) = write_cram(Vec::new(), &header, &fasta, &opts, |sink| {
        for i in 0..set.len() {
            sink.push(set.record(i))?;
        }
        Ok(())
    })
    .unwrap();
    assert_eq!(summary.records, set.len() as u64);
    assert!(bytes.starts_with(b"CRAM\x03\x00"));
    std::fs::write(&cram, &bytes).unwrap();
    let crai = dir.join("out.cram.crai");
    std::fs::write(&crai, crai_bytes(&summary.index)).unwrap();

    let bam_s = bam.to_str().unwrap();
    let cram_s = cram.to_str().unwrap();
    let want = run(&["view", "--no-PG", bam_s]);
    let got = run(&["view", "--no-PG", "-T", fa_s, cram_s]);
    assert_eq!(
        String::from_utf8_lossy(&got),
        String::from_utf8_lossy(&want),
        "records differ"
    );

    // Header: ours adds M5 to each @SQ; nothing else changes.
    let hdr = String::from_utf8(run(&["view", "-H", "--no-PG", "-T", fa_s, cram_s])).unwrap();
    let m5 = |seq: &[u8]| {
        let mut m = nomics_io::cram::md5::Md5::new();
        m.update_upper(seq);
        nomics_io::cram::md5::hex(&m.finish())
    };
    let expected = text
        .replace("LN:300\n", &format!("LN:300\tM5:{}\n", m5(&contigs[0].1)))
        .replace("LN:120\n", &format!("LN:120\tM5:{}\n", m5(&contigs[1].1)));
    assert_eq!(hdr, expected);

    // The index samtools builds for our CRAM equals ours.
    let copy = dir.join("copy.cram");
    std::fs::copy(&cram, &copy).unwrap();
    run(&["index", copy.to_str().unwrap()]);
    let gunzip = |p: &Path| {
        let mut out = Vec::new();
        nomics_io::gzip::gunzip_to(&std::fs::read(p).unwrap(), &mut |b: &[u8]| {
            out.extend_from_slice(b);
            Ok(())
        })
        .unwrap();
        String::from_utf8(out).unwrap()
    };
    let ours = gunzip(&crai);
    assert_eq!(ours, gunzip(&dir.join("copy.cram.crai")));
    assert!(
        ours.lines().count() >= 8,
        "expected several slices:\n{ours}"
    );

    // Region queries through our index match the BAM's through its .bai.
    run(&["index", bam_s]);
    for region in ["c1:1-40", "c1:100-160", "c1:280-300", "c2:60-75", "c2", "*"] {
        let a = run(&[
            "view",
            "--no-PG",
            "-T",
            fa_s,
            "-X",
            cram_s,
            crai.to_str().unwrap(),
            region,
        ]);
        let b = run(&["view", "--no-PG", bam_s, region]);
        assert!(!a.is_empty(), "region {region} returned nothing");
        assert_eq!(a, b, "region {region}");
    }
    // No decoder warnings (e.g. a slice past the contig end).
    let out = Command::new("samtools")
        .args(["view", "-c", "-T", fa_s, cram_s])
        .output()
        .unwrap();
    assert!(
        out.stderr.is_empty(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );

    // sarek-md profile: decodes as htslib's own `samtools view -C` does.
    // `samtools calmd` gives most reads htslib-equal MD/NM (dropped and
    // regenerated); IUPAC/N-run/off-end reads keep theirs verbatim.
    let calmd_bytes = run(&["calmd", "-b", bam_s, fa_s]);
    let calmd = dir.join("calmd.bam");
    std::fs::write(&calmd, calmd_bytes).unwrap();
    let calmd_s = calmd.to_str().unwrap();
    let cset = RecordSet::load(&calmd, 1).unwrap();
    let (sh, sopts) = nomics_io::cram::Profile::SarekMd
        .prepare(
            &cset.header,
            &fa,
            Some(("test", "cram_roundtrip")),
            CramOptions {
                records_per_slice: 4,
                ..Default::default()
            },
        )
        .unwrap();
    let (sbytes, _) = write_cram(Vec::new(), &sh, &fasta, &sopts, |sink| {
        for i in 0..cset.len() {
            sink.push(cset.record(i))?;
        }
        Ok(())
    })
    .unwrap();
    let ours_md = dir.join("sarek_md.cram");
    std::fs::write(&ours_md, sbytes).unwrap();
    let st_md = dir.join("samtools_md.cram");
    run(&[
        "view",
        "-C",
        "-T",
        fa_s,
        "--output-fmt-option",
        "version=3.0",
        "-o",
        st_md.to_str().unwrap(),
        calmd_s,
    ]);
    let ours_txt = run(&["view", "--no-PG", "-T", fa_s, ours_md.to_str().unwrap()]);
    let st_txt = run(&["view", "--no-PG", "-T", fa_s, st_md.to_str().unwrap()]);
    // Two records where htslib's own CRAM is lossy and we are not (see
    // docs/GATES.md): `tie` (TLEN against the 5' convention, attached by
    // htslib's encoder and rebuilt with the other sign by its decoder) and
    // `alltags` (a Z value with a TAB, which htslib's TAB-stopped encoding
    // truncates). Everything else must decode identically.
    let keep = |t: &[u8]| -> String {
        String::from_utf8_lossy(t)
            .lines()
            .filter(|l| !l.starts_with("tie\t") && !l.starts_with("alltags\t"))
            .map(|l| format!("{l}\n"))
            .collect()
    };
    assert_eq!(
        keep(&ours_txt),
        keep(&st_txt),
        "sarek-md decode differs from samtools view -C"
    );
    assert_ne!(ours_txt, st_txt, "the two lossy htslib cases should show");
    let moved = String::from_utf8_lossy(&st_txt)
        .lines()
        .filter(|l| l.ends_with("RG:Z:rg1") && l.contains("\tMD:Z:"))
        .count();
    assert!(
        moved > 5,
        "expected MD/NM regenerated on decode for most reads"
    );

    if std::env::var_os("NOMICS_KEEP_TEST_DIR").is_none() {
        std::fs::remove_dir_all(&dir).unwrap();
    } else {
        eprintln!("kept {}", dir.display());
    }
}

/// One-contig scratch setup for the edge-case tests below.
struct Scratch {
    dir: PathBuf,
    fa: PathBuf,
    fasta: Fasta,
    header: Header,
    seq: Vec<u8>,
}

fn scratch(tag: &str, len: usize) -> Scratch {
    let dir = std::env::temp_dir().join(format!("nomics-cram-{tag}-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let seq = contig(9, len);
    let fa = write_fasta(&dir, &[("c1", seq.clone())]);
    let fasta = Fasta::load(&fa, 1).unwrap();
    let header = Header {
        text: format!("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c1\tLN:{len}\n@RG\tID:rg1\tSM:s\n")
            .into_bytes(),
        refs: vec![Reference {
            name: b"c1".to_vec(),
            len: len as u32,
        }],
    };
    Scratch {
        dir,
        fa,
        fasta,
        header,
        seq,
    }
}

fn bare(name: &str, flag: u16, ref_id: i32, pos: i32, cigar: Vec<u32>, seq: Vec<u8>) -> Record {
    let qual = vec![30; seq.len()];
    Record {
        ref_id,
        pos,
        mapq: if flag & 4 != 0 { 0 } else { 60 },
        flag,
        next_ref_id: -1,
        next_pos: -1,
        tlen: 0,
        read_name: name.as_bytes().to_vec(),
        cigar,
        seq,
        qual,
        tags: vec![(*b"RG", Value::String(b"rg1".to_vec()))],
    }
}

fn raw_set(header: &Header, recs: &[Record]) -> RecordSet {
    let mut raw = Vec::new();
    let mut body = Vec::new();
    for r in recs {
        body.clear();
        r.encode(nomics_io::bam::BinRule::Sam, &mut body).unwrap();
        raw.extend_from_slice(&(body.len() as i32).to_le_bytes());
        raw.extend_from_slice(&body);
    }
    RecordSet::from_parts(header.clone(), &raw).unwrap()
}

fn encode_set(s: &Scratch, set: &RecordSet, opts: &CramOptions) -> std::io::Result<Vec<u8>> {
    write_cram(Vec::new(), &set.header, &s.fasta, opts, |sink| {
        (0..set.len()).try_for_each(|i| sink.push(set.record(i)))
    })
    .map(|(b, _)| b)
}

fn sarek_md_opts(s: &Scratch) -> CramOptions {
    nomics_io::cram::Profile::SarekMd
        .prepare(
            &s.header,
            &s.fa,
            Some(("test", "test")),
            CramOptions::default(),
        )
        .unwrap()
        .1
}

/// Decode `cram` and the BAM `set` with samtools and assert equal records.
fn assert_decodes_like(s: &Scratch, set: &RecordSet, cram: &[u8], label: &str) {
    let bam = s.dir.join(format!("{label}.bam"));
    set.write(&bam, 1).unwrap();
    let c = s.dir.join(format!("{label}.cram"));
    std::fs::write(&c, cram).unwrap();
    let fa = s.fa.to_str().unwrap();
    let got = run(&["view", "--no-PG", "-T", fa, c.to_str().unwrap()]);
    let want = run(&["view", "--no-PG", bam.to_str().unwrap()]);
    assert!(!want.is_empty(), "{label}: empty BAM");
    assert_eq!(
        String::from_utf8_lossy(&got),
        String::from_utf8_lossy(&want),
        "{label}"
    );
}

/// More than 65,535 CIGAR ops: BAM stores `<l_seq>S<span>N` plus `CG:B:I`;
/// the CRAM must carry the real CIGAR (as samtools does) and no CG tag.
#[test]
fn long_cigar_from_cg_tag() {
    if !samtools_ok() {
        return;
    }
    let pairs = 35_000usize;
    let s = scratch("cg", 2 * pairs + 100);
    let mut real = Vec::new();
    let mut seq = Vec::new();
    for k in 0..pairs {
        real.push(1 << 4); // 1M
        real.push(1 << 4 | 1); // 1I
        seq.push(s.seq[k].to_ascii_uppercase());
        seq.push(b'T');
    }
    let l_seq = seq.len() as u32;
    let placeholder = vec![l_seq << 4 | 4, (pairs as u32) << 4 | 3];
    let mut r = bare("long", 0, 0, 0, placeholder, seq);
    r.tags.push((
        *b"CG",
        Value::Array(ArrayValue::Int(
            IntType::U32,
            real.iter().map(|&c| c as i64).collect(),
        )),
    ));
    let set = raw_set(&s.header, &[r]);
    for (label, opts) in [
        ("cg_lossless", CramOptions::default()),
        ("cg_md", sarek_md_opts(&s)),
    ] {
        let cram = encode_set(&s, &set, &opts).unwrap();
        if label == "cg_lossless" {
            assert_decodes_like(&s, &set, &cram, label);
        } else {
            // sarek-md decodes as samtools' own CRAM of the BAM does (MD/NM
            // regenerated): compare with that.
            std::fs::write(s.dir.join("cg_md.cram"), &cram).unwrap();
            let bam = s.dir.join("cg.bam");
            set.write(&bam, 1).unwrap();
            let st = s.dir.join("cg_st.cram");
            let fa = s.fa.to_str().unwrap();
            run(&[
                "view",
                "-C",
                "-T",
                fa,
                "--output-fmt-option",
                "version=3.0",
                "-o",
                st.to_str().unwrap(),
                bam.to_str().unwrap(),
            ]);
            let ours = run(&[
                "view",
                "--no-PG",
                "-T",
                fa,
                s.dir.join("cg_md.cram").to_str().unwrap(),
            ]);
            let theirs = run(&["view", "--no-PG", "-T", fa, st.to_str().unwrap()]);
            assert!(ours == theirs, "cg_md decodes unlike samtools view -C");
        }
        let c = s.dir.join(format!("{label}.cram"));
        let text = run(&["view", "-T", s.fa.to_str().unwrap(), c.to_str().unwrap()]);
        assert!(
            !String::from_utf8_lossy(&text).contains("CG:B"),
            "{label}: CG kept"
        );
    }
    std::fs::remove_dir_all(&s.dir).unwrap();
}

/// A slice of only SEQ=*/QUAL=* unmapped records still declares QS.
#[test]
fn slice_without_any_bases() {
    if !samtools_ok() {
        return;
    }
    let s = scratch("noseq", 200);
    let recs: Vec<Record> = (0..3)
        .map(|i| bare(&format!("u{i}"), 4, -1, -1, vec![], vec![]))
        .collect();
    let set = raw_set(&s.header, &recs);
    let cram = encode_set(&s, &set, &CramOptions::default()).unwrap();
    assert_decodes_like(&s, &set, &cram, "noseq");
    std::fs::remove_dir_all(&s.dir).unwrap();
}

/// Malformed records are errors (not panics) in every tag mode.
#[test]
fn malformed_records_are_errors() {
    let s = scratch("bad", 300);
    let seq1 = vec![s.seq[10].to_ascii_uppercase()];
    let mut mate = bare("mate", 1 | 8, 0, 10, vec![1 << 4], seq1.clone());
    mate.next_ref_id = 0;
    mate.next_pos = i32::MAX;
    let cases: Vec<(&str, Record)> = vec![
        (
            "CIGAR longer than SEQ",
            bare("short", 0, 0, 10, vec![100 << 4], seq1.clone()),
        ),
        (
            "negative position",
            bare("neg", 0, 0, -5, vec![1 << 4], seq1.clone()),
        ),
        ("mate position overflow", mate),
    ];
    for (what, r) in cases {
        let set = raw_set(&s.header, &[r]);
        for opts in [CramOptions::default(), sarek_md_opts(&s)] {
            let tags = opts.tags;
            let got = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                encode_set(&s, &set, &opts)
            }));
            match got {
                Ok(Err(_)) => {}
                Ok(Ok(_)) => panic!("{what} ({tags:?}): encoded"),
                Err(_) => panic!("{what} ({tags:?}): panicked"),
            }
        }
    }
    std::fs::remove_dir_all(&s.dir).unwrap();
}

/// Slices also end at `max_slice_bytes`.
#[test]
fn slices_capped_by_bytes() {
    if !samtools_ok() {
        return;
    }
    let s = scratch("bytes", 2000);
    let recs: Vec<Record> = (0..50)
        .map(|i| {
            let at = (i * 20) as usize;
            let seq = s.seq[at..at + 50]
                .iter()
                .map(|b| b.to_ascii_uppercase())
                .collect();
            bare(&format!("r{i}"), 0, 0, i * 20, vec![50 << 4], seq)
        })
        .collect();
    let set = raw_set(&s.header, &recs);
    let opts = CramOptions {
        max_slice_bytes: 1000,
        ..Default::default()
    };
    let (cram, summary) = write_cram(Vec::new(), &set.header, &s.fasta, &opts, |sink| {
        (0..set.len()).try_for_each(|i| sink.push(set.record(i)))
    })
    .unwrap();
    assert!(summary.containers >= 5, "{} containers", summary.containers);
    assert_decodes_like(&s, &set, &cram, "bytes");
    std::fs::remove_dir_all(&s.dir).unwrap();
}

/// Boundary cases of htslib's `bam_tag2cigar`: our CRAM must decode as
/// samtools reads the BAM (it applies the same rule when reading).
#[test]
fn cg_boundary_cases_match_htslib() {
    if !samtools_ok() {
        return;
    }
    let s = scratch("cgb", 400);
    let m = |n: u32| n << 4;
    let ins = |n: u32| n << 4 | 1;
    let soft = |n: u32| n << 4 | 4;
    let skip = |n: u32| n << 4 | 3;
    let cg = |ty: IntType, ops: &[u32]| {
        (
            *b"CG",
            Value::Array(ArrayValue::Int(ty, ops.iter().map(|&c| c as i64).collect())),
        )
    };
    let bases = |at: usize, n: usize| -> Vec<u8> {
        s.seq[at..at + n]
            .iter()
            .map(|b| b.to_ascii_uppercase())
            .collect()
    };
    let mut recs = Vec::new();
    // (name, cigar, CG tags, expanded?)
    let cases: Vec<(&str, Vec<u32>, Vec<([u8; 2], Value)>)> = vec![
        // CG shorter than the placeholder: untouched (CG kept).
        (
            "short_cg",
            vec![soft(10), skip(10)],
            vec![cg(IntType::U32, &[m(10)])],
        ),
        // Signed B:i is accepted.
        (
            "signed_cg",
            vec![soft(10), skip(8)],
            vec![cg(IntType::I32, &[m(5), ins(2), m(3)])],
        ),
        // Exactly as long as the placeholder.
        (
            "equal_len_cg",
            vec![soft(10), skip(6)],
            vec![cg(IntType::U32, &[m(6), ins(4)])],
        ),
        // Wrong subtype: untouched.
        (
            "ushort_cg",
            vec![soft(10), skip(10)],
            vec![cg(IntType::U16, &[m(10)])],
        ),
        // First op not a whole-read soft clip: untouched.
        (
            "partial_clip",
            vec![soft(5), m(5)],
            vec![cg(IntType::U32, &[m(5), ins(5)])],
        ),
        // The first CG tag is not B: untouched, both tags kept.
        (
            "first_cg_z",
            vec![soft(10), skip(10)],
            vec![
                (*b"CG", Value::String(b"x".to_vec())),
                cg(IntType::U32, &[m(10)]),
            ],
        ),
    ];
    for (i, (name, cigar, tags)) in cases.into_iter().enumerate() {
        let mut r = bare(
            name,
            0,
            0,
            20 + 30 * i as i32,
            cigar,
            bases(20 + 30 * i, 10),
        );
        r.tags.extend(tags);
        recs.push(r);
    }
    let set = raw_set(&s.header, &recs);
    let cram = encode_set(&s, &set, &CramOptions::default()).unwrap();
    assert_decodes_like(&s, &set, &cram, "cg_boundaries");
    let text = String::from_utf8(run(&[
        "view",
        "--no-PG",
        "-T",
        s.fa.to_str().unwrap(),
        s.dir.join("cg_boundaries.cram").to_str().unwrap(),
    ]))
    .unwrap();
    for (name, expanded) in [
        ("short_cg", false),
        ("signed_cg", true),
        ("equal_len_cg", true),
        ("ushort_cg", false),
        ("partial_clip", false),
        ("first_cg_z", false),
    ] {
        let line = text
            .lines()
            .find(|l| l.starts_with(&format!("{name}\t")))
            .unwrap();
        assert_eq!(!line.contains("\tCG:"), expanded, "{name}: {line}");
    }
    std::fs::remove_dir_all(&s.dir).unwrap();
}

/// A record bigger than `max_slice_bytes` is its own slice; one bigger than
/// `MAX_RECORD_BYTES` is an error.
#[test]
fn oversized_records() {
    let s = scratch("big", 3000);
    let up = |a: usize, n: usize| -> Vec<u8> {
        s.seq[a..a + n]
            .iter()
            .map(|b| b.to_ascii_uppercase())
            .collect()
    };
    let recs = vec![
        bare("small1", 0, 0, 10, vec![50 << 4], up(10, 50)),
        bare("big", 0, 0, 100, vec![2000 << 4], up(100, 2000)),
        bare("small2", 0, 0, 2200, vec![50 << 4], up(2200, 50)),
    ];
    let set = raw_set(&s.header, &recs);
    let opts = CramOptions {
        max_slice_bytes: 1000,
        ..Default::default()
    };
    let (cram, summary) = write_cram(Vec::new(), &set.header, &s.fasta, &opts, |sink| {
        (0..set.len()).try_for_each(|i| sink.push(set.record(i)))
    })
    .unwrap();
    assert_eq!(summary.containers, 3, "the big record forms its own slice");
    if samtools_ok() {
        assert_decodes_like(&s, &set, &cram, "big");
    }
    // Over the hard limit: rejected before any parsing.
    let mut huge = set.record(0).to_vec();
    huge.resize(nomics_io::cram::MAX_RECORD_BYTES + 1, 0);
    let got = write_cram(Vec::new(), &set.header, &s.fasta, &opts, |sink| {
        sink.push(&huge)
    });
    assert!(got.is_err());
    std::fs::remove_dir_all(&s.dir).unwrap();
}

/// The unplaced tail spread over several slices must all come back through
/// the index (`samtools view '*'`, and htsjdk's unmapped query when the GATK
/// image is present). Regression: with span 0 on unplaced slices, htslib's
/// `cram_index_query` skipped all but the last of them.
#[test]
fn unplaced_tail_over_several_slices_is_indexed() {
    if !samtools_ok() {
        return;
    }
    let s = scratch("unplaced", 500);
    let mut recs: Vec<Record> = (0..6)
        .map(|i| {
            let at = 10 + 40 * i as usize;
            let seq = s.seq[at..at + 30]
                .iter()
                .map(|b| b.to_ascii_uppercase())
                .collect();
            bare(&format!("m{i}"), 0, 0, at as i32, vec![30 << 4], seq)
        })
        .collect();
    let unplaced = 47;
    for i in 0..unplaced {
        recs.push(bare(
            &format!("u{i}"),
            4,
            -1,
            -1,
            vec![],
            b"ACGTNACGTA".to_vec(),
        ));
    }
    let set = raw_set(&s.header, &recs);
    let opts = CramOptions {
        records_per_slice: 10,
        ..Default::default()
    };
    let (cram, summary) = write_cram(Vec::new(), &set.header, &s.fasta, &opts, |sink| {
        (0..set.len()).try_for_each(|i| sink.push(set.record(i)))
    })
    .unwrap();
    let unplaced_slices = summary.index.iter().filter(|e| e.ref_id == -1).count();
    assert!(unplaced_slices >= 4, "{unplaced_slices} unplaced slices");
    assert_decodes_like(&s, &set, &cram, "unplaced");
    let c = s.dir.join("unplaced.cram");
    std::fs::write(c.with_extension("cram.crai"), crai_bytes(&summary.index)).unwrap();
    let bam = s.dir.join("unplaced.bam");
    run(&["index", bam.to_str().unwrap()]);
    let (fa, c_s, bam_s) = (
        s.fa.to_str().unwrap(),
        c.to_str().unwrap(),
        bam.to_str().unwrap(),
    );
    let ours = run(&["view", "--no-PG", "-T", fa, c_s, "*"]);
    let theirs = run(&["view", "--no-PG", bam_s, "*"]);
    assert_eq!(String::from_utf8_lossy(&ours).lines().count(), unplaced);
    assert_eq!(ours, theirs, "'*' query content");
    for e in summary.index.iter().filter(|e| e.ref_id == -1) {
        assert_eq!(
            (e.start, e.span),
            (0, 1),
            "unplaced slice extent as htslib writes it"
        );
    }

    // Our index must also equal what samtools builds.
    let copy = s.dir.join("reindex.cram");
    std::fs::copy(&c, &copy).unwrap();
    run(&["index", copy.to_str().unwrap()]);
    let gunzip = |p: PathBuf| {
        let mut out = Vec::new();
        nomics_io::gzip::gunzip_to(&std::fs::read(p).unwrap(), &mut |b: &[u8]| {
            out.extend_from_slice(b);
            Ok(())
        })
        .unwrap();
        out
    };
    assert_eq!(
        gunzip(c.with_extension("cram.crai")),
        gunzip(copy.with_extension("cram.crai"))
    );
    // htsjdk 4 (GATK 4.6.2.0) unmapped query, when its image is present.
    let image = "community.wave.seqera.io/library/gatk4_htslib_samtools:6fd46722be032bdd";
    let have = Command::new("docker")
        .args(["image", "inspect", image])
        .output()
        .is_ok_and(|o| o.status.success());
    if !have {
        eprintln!("GATK image not present: htsjdk unmapped query skipped");
    } else {
        let dict = run(&["dict", fa]);
        std::fs::write(s.fa.with_extension("dict"), dict).unwrap();
        let d = std::fs::canonicalize(&s.dir).unwrap();
        let dm = format!("{}:{}", d.display(), d.display());
        let out = Command::new("docker")
            .args([
                "run",
                "--rm",
                "-v",
                &dm,
                image,
                "gatk",
                "PrintReads",
                "-I",
                &d.join("unplaced.cram").display().to_string(),
                "-R",
                &d.join("ref.fa").display().to_string(),
                "-L",
                "unmapped",
                "-O",
                &d.join("gatk_unmapped.sam").display().to_string(),
                "--disable-read-filter",
                "WellformedReadFilter",
            ])
            .output()
            .unwrap();
        assert!(
            out.status.success(),
            "{}",
            String::from_utf8_lossy(&out.stderr)
        );
        let sam = std::fs::read_to_string(d.join("gatk_unmapped.sam")).unwrap();
        assert_eq!(
            sam.lines().filter(|l| !l.starts_with('@')).count(),
            unplaced,
            "htsjdk unmapped query"
        );
    }
    std::fs::remove_dir_all(&s.dir).unwrap();
}

/// Codec trials: enough full slices that most containers follow an earlier
/// trial container's codec choices. The bytes must not depend on the thread
/// count (the choice is a function of container indices and data only), and
/// the records must decode unchanged.
#[test]
fn codec_trials_are_thread_independent() {
    let s = scratch("trials", 200_000);
    let mut x = 7u64;
    let mut rnd = || {
        x = x
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        x >> 33
    };
    let recs: Vec<Record> = (0..90_000)
        .map(|i| {
            let at = (i * 2) as usize;
            let mut seq: Vec<u8> = s.seq[at..at + 30]
                .iter()
                .map(|b| b.to_ascii_uppercase())
                .collect();
            if rnd() % 7 == 0 {
                seq[rnd() as usize % 30] = b"ACGT"[rnd() as usize % 4];
            }
            let mut r = bare(&format!("r{i}"), 0, 0, at as i32, vec![30 << 4], seq);
            // Quality distributions that drift along the file.
            let spread = 2 + (i / 20_000) as u64 * 8;
            r.qual = (0..30).map(|_| (10 + rnd() % spread) as u8).collect();
            r
        })
        .collect();
    let set = raw_set(&s.header, &recs);
    let encode = |threads| {
        let opts = CramOptions {
            records_per_slice: 1000,
            threads,
            ..Default::default()
        };
        encode_set(&s, &set, &opts).unwrap()
    };
    let one = encode(1);
    for threads in [2, 5] {
        assert!(encode(threads) == one, "{threads} threads differ from 1");
    }
    if samtools_ok() {
        assert_decodes_like(&s, &set, &one, "trials");
    }
}

/// Finished containers waiting for an earlier one stay within the reorder
/// window, however the workers race.
#[test]
fn reorder_buffer_is_bounded() {
    let s = scratch("reorder", 50_000);
    let recs: Vec<Record> = (0..20_000)
        .map(|i| {
            let at = (i * 2) as usize;
            let seq = s.seq[at..at + 30].iter().map(|b| b.to_ascii_uppercase()).collect();
            bare(&format!("r{i}"), 0, 0, at as i32, vec![30 << 4], seq)
        })
        .collect();
    let set = raw_set(&s.header, &recs);
    for threads in [1, 3, 8] {
        let opts = CramOptions {
            records_per_slice: 50,
            threads,
            ..Default::default()
        };
        let (_, summary) = write_cram(Vec::new(), &set.header, &s.fasta, &opts, |sink| {
            (0..set.len()).try_for_each(|i| sink.push(set.record(i)))
        })
        .unwrap();
        assert_eq!(summary.containers, 400);
        assert!(
            summary.reorder_peak as u64 <= nomics_io::cram::reorder_window(threads),
            "{threads} threads: {} pending",
            summary.reorder_peak
        );
    }
}
