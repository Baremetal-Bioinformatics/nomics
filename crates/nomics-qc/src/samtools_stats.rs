//! samtools 1.24 `stats --threads 1 --reference REF INPUT` (BAM mode).
//! `--reference` is the global decoder option, NOT stats' `--ref-seq`.
//! Source audit and boundary quirks: docs/spec/qc-bam-stats.md.
mod report;
mod ring;
mod stats;

use nomics_io::{bam, bgzf::ParallelBgzfReader};
pub use stats::{Record, Stats};
use std::{
    fs::File,
    io::{self, BufReader},
    path::Path,
};

fn invalid(s: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, s)
}

/// Stream BAM through nomics-io. The reference is only a provenance name in
/// this BAM command; upstream never opens it (unlike `stats -r`).
pub fn bam_stats(
    bam_path: &Path,
    reference: &Path,
    out: &Path,
    input_name: Option<&str>,
    reference_name: Option<&str>,
) -> io::Result<()> {
    let basename = |p: &Path| {
        p.file_name()
            .and_then(|s| s.to_str())
            .map(str::to_owned)
            .ok_or_else(|| invalid("input basename is not UTF-8"))
    };
    let input_name = input_name
        .map(str::to_owned)
        .map(Ok)
        .unwrap_or_else(|| basename(bam_path))?;
    let reference_name = reference_name
        .map(str::to_owned)
        .map(Ok)
        .unwrap_or_else(|| basename(reference))?;
    let file = BufReader::with_capacity(1 << 20, File::open(bam_path)?);
    let mut reader = ParallelBgzfReader::new(file, 2);
    let _header = bam::Header::read(&mut reader)?;
    let mut stats = Stats::new();
    let mut raw = Vec::new();
    let mut cigar = Vec::new();
    while bam::read_raw_record(&mut reader, &mut raw)? {
        let v = bam::RecordView::new(&raw)?;
        let quals = v
            .quals()
            .ok_or_else(|| invalid("truncated BAM sequence/qualities"))?;
        cigar.clear();
        cigar.extend(v.cigar());
        let mut nm = None;
        let mut cg_seen = false;
        let mut barcode_seen = [false; 4];
        let mut barqual_seen = [false; 4];
        let mut barcode = [None; 4];
        let mut barqual = [None; 4];
        let mut aux = v.aux().ok_or_else(|| invalid("truncated BAM tags"))?;
        while aux.len() >= 3 {
            let tag = &aux[..2];
            let ty = aux[2];
            let b = &aux[3..];
            let len = (if ty == b'd' {
                Some(8)
            } else {
                bam::aux_value_len(ty, b)
            })
            .filter(|&n| n <= b.len())
            .ok_or_else(|| invalid("invalid BAM auxiliary field"))?;
            if tag == b"NM" && nm.is_none() {
                nm = Some(match ty {
                    b'c' => b[0] as i8 as i64,
                    b'C' => b[0] as i64,
                    b's' => i16::from_le_bytes(b[..2].try_into().unwrap()) as i64,
                    b'S' => u16::from_le_bytes(b[..2].try_into().unwrap()) as i64,
                    b'i' => i32::from_le_bytes(b[..4].try_into().unwrap()) as i64,
                    b'I' => u32::from_le_bytes(b[..4].try_into().unwrap()) as i64,
                    _ => 0,
                });
            }
            if tag == b"CG"
                && !cg_seen
                && v.ref_id() >= 0
                && v.pos() >= 0
                && ty == b'B'
                && matches!(b[0], b'I' | b'i')
                && cigar.first().copied() == Some((v.len() as u32) << 4 | 4)
            {
                let n = u32::from_le_bytes(b[1..5].try_into().unwrap()) as usize;
                if n >= cigar.len() && n < (1 << 29) {
                    cigar.clear();
                    cigar.extend(
                        b[5..len]
                            .as_chunks::<4>()
                            .0
                            .iter()
                            .map(|s| u32::from_le_bytes(*s)),
                    );
                }
            }
            if tag == b"CG" {
                cg_seen = true;
            }
            {
                for (i, (bt, qt)) in [
                    (b"BC", b"QT"),
                    (b"CR", b"CY"),
                    (b"OX", b"BZ"),
                    (b"RX", b"QX"),
                ]
                .iter()
                .enumerate()
                {
                    if tag == *bt && !barcode_seen[i] {
                        barcode_seen[i] = true;
                        if matches!(ty, b'Z' | b'H') {
                            barcode[i] = Some(&b[..len - 1]);
                        }
                    }
                    if tag == *qt && !barqual_seen[i] {
                        barqual_seen[i] = true;
                        if matches!(ty, b'Z' | b'H') {
                            barqual[i] = Some(&b[..len - 1]);
                        }
                    }
                }
            }
            aux = &b[len..];
        }
        if v.flag() & 4 == 0 && !v.is_empty() && !cigar.is_empty() {
            let consumed: u64 = cigar
                .iter()
                .filter(|c| matches!(**c & 15, 0 | 1 | 4 | 7 | 8))
                .map(|c| (c >> 4) as u64)
                .sum();
            if consumed != v.len() as u64 {
                return Err(invalid("CIGAR and query sequence lengths differ"));
            }
        }
        if raw[8] == 0 {
            return Err(invalid("empty BAM qname field"));
        }
        let name = &raw[32..32 + raw[8] as usize];
        let name = name.split(|&b| b == 0).next().unwrap_or_default();
        let end = v.aux_offset() - v.len();
        let packed = &raw[end - v.len().div_ceil(2)..end];
        stats.add(&Record {
            name,
            packed,
            quals,
            cigar: &cigar,
            tid: v.ref_id(),
            pos: v.pos() as i64,
            mate_tid: v.mate_ref_id(),
            mate_pos: v.mate_pos() as i64,
            tlen: v.tlen() as i64,
            flag: v.flag(),
            mapq: v.mapq(),
            nm,
            barcode,
            barqual,
        })?;
    }
    let text = stats.render(&input_name, &reference_name)?;
    std::fs::write(out, text)
}

#[cfg(test)]
mod tests {
    use super::*;
    use nomics_io::{
        bam::{ArrayValue, BinRule, IntType, Value},
        bgzf::{BgzfWriter, Profile, Stored},
    };
    use std::io::Write;
    #[test]
    fn bam_adapter_promotes_cg_and_keeps_packed_padding_without_loading_reference() {
        let dir = std::env::temp_dir().join(format!("nomics-bam-adapter-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let input = dir.join("physical.bam");
        let output = dir.join("out.stats");
        let mut stream = Vec::new();
        bam::Header {
            text: b"@HD\tVN:1.6\n@RG\tID:a\n@RG\tID:b\n".to_vec(),
            refs: vec![
                bam::Reference {
                    name: b"one".to_vec(),
                    len: 1000,
                },
                bam::Reference {
                    name: b"unused".to_vec(),
                    len: 1000,
                },
            ],
        }
        .write(&mut stream)
        .unwrap();
        let mut r = bam::Record {
            ref_id: 0,
            pos: 10,
            mapq: 60,
            flag: 0,
            next_ref_id: -1,
            next_pos: -1,
            tlen: 0,
            read_name: b"read".to_vec(),
            cigar: vec![52, 51],
            seq: b"ACG".to_vec(),
            qual: vec![255; 3],
            tags: vec![
                (*b"RG", Value::String(b"a".to_vec())),
                (
                    *b"CG",
                    Value::Array(ArrayValue::Int(IntType::U32, vec![16, 16, 16])),
                ),
                (*b"NM", Value::Int(IntType::U8, 2)),
            ],
        };
        let mut raw = Vec::new();
        r.encode(BinRule::Sam, &mut raw).unwrap();
        let v = bam::RecordView::new(&raw).unwrap();
        let end = v.aux_offset() - v.len();
        raw[end - 1] |= 0xf; // nonzero padding nibble participates in CHK.
        let packed_crc = nomics_io::crc32::crc32(&raw[end - 2..end]);
        bam::write_raw_record(&mut stream, &raw).unwrap();
        r.pos = 20;
        r.tags = vec![(*b"RG", Value::String(b"b".to_vec()))];
        r.cigar = vec![48];
        raw.clear();
        r.encode(BinRule::Sam, &mut raw).unwrap();
        bam::write_raw_record(&mut stream, &raw).unwrap();
        let mut writer = BgzfWriter::new(Vec::new(), Stored, Profile::HTSJDK);
        writer.write_all(&stream).unwrap();
        std::fs::write(&input, writer.finish().unwrap()).unwrap();
        bam_stats(
            &input,
            &dir.join("nonexistent.fasta"),
            &output,
            Some("logical.bam"),
            Some("logical.fa"),
        )
        .unwrap();
        let text = std::fs::read_to_string(output).unwrap();
        assert!(text.contains("--reference logical.fa logical.bam\n"));
        assert!(text.contains("SN\tsequences:\t2\n"));
        assert!(text.contains("SN\tmismatches:\t2\t"));
        assert!(text.contains("SN\tbases mapped (cigar):\t6\t"));
        let expected = packed_crc.wrapping_add(nomics_io::crc32::crc32(&[0x12, 0x40]));
        assert_eq!(
            text.lines()
                .find(|s| s.starts_with("CHK\t"))
                .unwrap()
                .split('\t')
                .nth(2)
                .unwrap(),
            format!("{expected:08x}")
        );
        assert!(!text.contains("MPC\t"));
        std::fs::remove_dir_all(dir).unwrap();
    }
}
