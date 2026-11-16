//! File-order BAM projection. Integer storage widths and tag storage order are
//! normalized to SAM values; float payloads remain raw bits. AR/CR values are fixed
//! here even though attached inside callMutations. Only called-allele XA is excluded.
use super::{Aligned, cigar};
use crate::mutect2_regions::Inputs;
use nomics_io::{
    bam::{self, ArrayValue, Header, IntType, Record, RecordView, Value},
    bgzf::BgzfReader,
};
use std::{
    collections::{HashMap, HashSet},
    fs::File,
    io::{self, BufReader, BufWriter, Write},
    path::Path,
};

fn value(v: &Value) -> String {
    match v {
        Value::Int(_, v) => format!("i:{v}"),
        Value::Char(v) => format!("A:{}", *v as char),
        Value::Float(v) => format!("f:{v:08x}"),
        Value::String(v) => format!("Z:{}", String::from_utf8_lossy(v)),
        Value::Hex(v) => format!("H:{}", String::from_utf8_lossy(v)),
        Value::Array(ArrayValue::Int(t, v)) => {
            let ty = match t {
                IntType::I8 => 'c',
                IntType::U8 => 'C',
                IntType::I16 => 's',
                IntType::U16 => 'S',
                IntType::I32 => 'i',
                IntType::U32 => 'I',
            };
            format!(
                "B:{ty},{}",
                v.iter().map(i64::to_string).collect::<Vec<_>>().join(",")
            )
        }
        Value::Array(ArrayValue::Float(v)) => format!(
            "B:f,{}",
            v.iter()
                .map(|v| format!("{v:08x}"))
                .collect::<Vec<_>>()
                .join(",")
        ),
    }
}
fn line(w: &mut impl Write, header: &Header, r: &Record) -> io::Result<()> {
    let reference = |id: i32| {
        if id < 0 {
            "*".into()
        } else {
            String::from_utf8_lossy(&header.refs[id as usize].name).into_owned()
        }
    };
    write!(
        w,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        String::from_utf8_lossy(&r.read_name),
        r.flag,
        reference(r.ref_id),
        r.pos + 1,
        r.mapq,
        cigar(&r.cigar),
        reference(r.next_ref_id),
        r.next_pos + 1,
        r.tlen,
        String::from_utf8_lossy(&r.seq),
        r.qual
            .iter()
            .map(|q| format!("{q:02x}"))
            .collect::<String>()
    )?;
    let mut tags: Vec<_> = r.tags.iter().collect();
    tags.sort_by_key(|(t, _)| *t);
    for (t, v) in tags {
        write!(w, "\t{}:{}", String::from_utf8_lossy(t), value(v))?;
    }
    writeln!(w)
}
/// Project all biological records, without sorting: BAM file order is the gate.
pub fn project(bam: &Path, output: &Path) -> io::Result<usize> {
    let mut stream = BgzfReader::new(BufReader::new(File::open(bam)?));
    let header = Header::read(&mut stream)?;
    let mut raw = Vec::new();
    let mut w = BufWriter::new(File::create(output)?);
    let mut count = 0;
    while bam::read_raw_record(&mut stream, &mut raw)? {
        let r = RecordView::new(&raw)?;
        if r.string_tag(b"RG") == Some(b"ArtificialHaplotypeRG") {
            continue;
        }
        line(&mut w, &header, &Record::decode(&raw)?)?;
        count += 1;
    }
    w.flush()?;
    Ok(count)
}
pub(super) fn write(
    inputs: &Inputs<'_>,
    mut aligned: Vec<Aligned>,
    output: &Path,
) -> io::Result<()> {
    // Recover only tags for evaluated records. Stage one retains no extra BAM
    // payload, so existing traversal and assembly avoid this second I/O pass.
    let wanted: HashSet<_> = aligned.iter().map(|a| a.source).collect();
    let mut originals = HashMap::new();
    let headers: Vec<_> = inputs
        .input
        .iter()
        .map(|p| nomics_io::bam_query::BamQuery::open(p, &[]).map(|q| q.header))
        .collect::<io::Result<_>>()?;
    let pg = super::header::translations(&headers, b"@PG")?;
    let rg = super::header::translations(&headers, b"@RG")?;
    for (file, path) in inputs.input.iter().enumerate() {
        let mut stream = BgzfReader::new(BufReader::new(File::open(path)?));
        let _ = Header::read(&mut stream)?;
        let mut raw = Vec::new();
        while bam::read_raw_record(&mut stream, &mut raw)? {
            let key = (file, stream.virtual_offset());
            if wanted.contains(&key) {
                let mut record = Record::decode(&raw)?;
                for (tag, value) in &mut record.tags {
                    if let Value::String(id) = value {
                        let translations = if tag == b"PG" {
                            Some(&pg[file])
                        } else if tag == b"RG" {
                            Some(&rg[file])
                        } else {
                            None
                        };
                        if let Some(map) = translations {
                            *id = map
                                .get(id)
                                .ok_or_else(|| io::Error::other("tag ID missing from header"))?
                                .clone();
                        }
                    }
                }
                originals.insert(key, record);
            }
        }
    }
    aligned.sort_by(|a, b| nomics_hc::finalize::coordinate_cmp(&a.read, &b.read));
    let mut w = BufWriter::new(File::create(output)?);
    for a in aligned {
        let r = &a.read;
        let original = originals
            .get(&a.source)
            .ok_or_else(|| io::Error::other("original BAM record missing"))?;
        let mut record = original.clone();
        let indels = a.indel_tags || record.tags.iter().any(|(t, _)| matches!(t, b"BI" | b"BD"));
        record.tags.retain(|(t, _)| !matches!(t, b"BI" | b"BD"));
        if indels {
            for (tag, q) in [(*b"BI", a.ins), (*b"BD", a.del)] {
                record
                    .tags
                    .push((tag, Value::String(q.iter().map(|b| b + 33).collect())));
            }
        }
        if let Some(hc) = a.hc {
            record.tags.retain(|(t, _)| t != b"HC");
            record
                .tags
                .push((*b"HC", Value::Int(IntType::I32, i64::from(hc))));
        }
        if let Some(xa) = &a.xa {
            let mut value = record
                .tags
                .iter()
                .find_map(|(tag, value)| {
                    if tag == b"XA"
                        && let Value::String(v) = value
                    {
                        let mut v = v.clone();
                        v.extend_from_slice(b", ");
                        return Some(v);
                    }
                    None
                })
                .unwrap_or_default();
            value.extend_from_slice(xa.as_bytes());
            record.tags.retain(|(t, _)| t != b"XA");
            record.tags.push((*b"XA", Value::String(value)));
        }
        record.tags.retain(|(t, _)| !matches!(t, b"AR" | b"CR"));
        record
            .tags
            .push((*b"AR", Value::String(a.alignment_region.into_bytes())));
        record
            .tags
            .push((*b"CR", Value::String(a.callable_region.into_bytes())));
        record.ref_id = r.contig;
        record.pos = r.start as i32 - 1;
        record.flag = r.flag;
        record.mapq = r.mapq;
        record.next_ref_id = r.mate_contig;
        record.next_pos = r.mate_start as i32 - 1;
        record.tlen = r.tlen;
        record.cigar = r.cigar.clone();
        record.seq = r.bases.clone();
        record.qual = r.quals.clone();
        line(&mut w, &headers[0], &record)?;
    }
    w.flush()
}
