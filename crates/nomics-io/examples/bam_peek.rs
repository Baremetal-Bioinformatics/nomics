//! Print the first N records of a BAM with raw tag order and integer type codes
//! (which `samtools view` hides).
//!
//! `cargo run --release -p nomics-io --example bam_peek -- <file.bam> [N] [skip]`

use std::fs::File;
use std::io::BufReader;

use nomics_io::bam::{ArrayValue, Header, IntType, Record, Value, read_raw_record};
use nomics_io::bgzf::BgzfReader;

fn ty(it: IntType) -> char {
    match it {
        IntType::I8 => 'c',
        IntType::U8 => 'C',
        IntType::I16 => 's',
        IntType::U16 => 'S',
        IntType::I32 => 'i',
        IntType::U32 => 'I',
    }
}

fn main() {
    let mut args = std::env::args().skip(1);
    let path = args.next().expect("usage: bam_peek <bam> [n] [skip]");
    let n: usize = args.next().map_or(5, |s| s.parse().unwrap());
    let skip: usize = args.next().map_or(0, |s| s.parse().unwrap());
    let mut r = BgzfReader::new(BufReader::new(File::open(&path).unwrap()));
    let header = Header::read(&mut r).unwrap();
    for line in String::from_utf8_lossy(&header.text)
        .lines()
        .filter(|l| !l.starts_with("@SQ"))
    {
        println!("{}", &line[..line.len().min(160)]);
    }
    let mut raw = Vec::new();
    let mut i = 0;
    while i < skip + n && read_raw_record(&mut r, &mut raw).unwrap() {
        i += 1;
        if i <= skip {
            continue;
        }
        let rec = Record::decode(&raw).unwrap();
        let bin = u16::from_le_bytes([raw[10], raw[11]]);
        let tags: Vec<String> = rec
            .tags
            .iter()
            .map(|(t, v)| {
                let t = String::from_utf8_lossy(t);
                match v {
                    Value::Int(it, x) => format!("{t}:{}:{x}", ty(*it)),
                    Value::String(s) => format!("{t}:Z:{}", String::from_utf8_lossy(s)),
                    Value::Char(c) => format!("{t}:A:{}", *c as char),
                    Value::Float(b) => format!("{t}:f:{}", f32::from_bits(*b)),
                    Value::Hex(h) => format!("{t}:H:{}", String::from_utf8_lossy(h)),
                    Value::Array(ArrayValue::Int(it, xs)) => format!("{t}:B:{}{:?}", ty(*it), xs),
                    Value::Array(ArrayValue::Float(xs)) => format!("{t}:B:f[{}]", xs.len()),
                }
            })
            .collect();
        println!(
            "{} flag={} ref={} pos={} bin={bin} {}",
            String::from_utf8_lossy(&rec.read_name),
            rec.flag,
            rec.ref_id,
            rec.pos,
            tags.join(" ")
        );
    }
}
