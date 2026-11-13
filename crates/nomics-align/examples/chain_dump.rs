//! Print `c` and `C` lines (oracle dump format) for a chunk's read pairs.
//!
//! usage: chain_dump <index-hint> <R1.fastq.gz> <R2.fastq.gz> <out-prefix>
//! writes <out-prefix>.c.txt and <out-prefix>.C.txt
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use nomics_align::chain::{mem_chain, mem_chain_flt};
use nomics_align::index::Index;
use nomics_align::seed::{SmemAux, nt4};
use nomics_io::gzip::GzReader;

fn seqs(p: &str) -> Vec<Vec<u8>> {
    let r = BufReader::with_capacity(1 << 20, GzReader::open(Path::new(p)).unwrap());
    r.lines()
        .skip(1)
        .step_by(4)
        .map(|l| l.unwrap().bytes().map(nt4).collect())
        .collect()
}

/// glibc `printf("%a", (double)x)` for normal and zero values.
fn hexfloat(x: f64) -> String {
    if x == 0.0 {
        return if x.is_sign_negative() {
            "-0x0p+0".into()
        } else {
            "0x0p+0".into()
        };
    }
    let bits = x.to_bits();
    let sign = if bits >> 63 != 0 { "-" } else { "" };
    let exp = ((bits >> 52) & 0x7ff) as i64 - 1023;
    let mut mant = format!("{:013x}", bits & ((1 << 52) - 1));
    while mant.ends_with('0') {
        mant.pop();
    }
    let dot = if mant.is_empty() {
        String::new()
    } else {
        format!(".{mant}")
    };
    format!("{sign}0x1{dot}p{exp:+}")
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let idx = Index::load(Path::new(&a[1])).unwrap();
    let (r1, r2) = (seqs(&a[2]), seqs(&a[3]));
    let mut oc = BufWriter::new(File::create(format!("{}.c.txt", a[4])).unwrap());
    let mut of = BufWriter::new(File::create(format!("{}.C.txt", a[4])).unwrap());
    let mut aux = SmemAux::default();
    for (pair, (s1, s2)) in r1.iter().zip(&r2).enumerate() {
        for (mate, s) in [s1, s2].into_iter().enumerate() {
            let id = pair << 1 | mate;
            let mut chains = mem_chain(&idx, s, &mut aux);
            for (i, c) in chains.iter().enumerate() {
                write!(oc, "c {id} {i} {} {} {} |", c.pos, c.rid, c.seeds.len()).unwrap();
                for s in &c.seeds {
                    write!(oc, " {},{},{},{}", s.rbeg, s.qbeg, s.len, s.score).unwrap();
                }
                writeln!(oc).unwrap();
            }
            mem_chain_flt(&mut chains);
            for (i, c) in chains.iter().enumerate() {
                write!(
                    of,
                    "C {id} {i} {} {} {} {} {} {} {} |",
                    c.pos,
                    c.rid,
                    c.kept,
                    c.first,
                    c.w,
                    hexfloat(c.frac_rep as f64),
                    c.seeds.len()
                )
                .unwrap();
                for s in &c.seeds {
                    write!(of, " {},{},{},{}", s.rbeg, s.qbeg, s.len, s.score).unwrap();
                }
                writeln!(of).unwrap();
            }
        }
    }
}
