//! Indexed FASTA (`.fai`) loaded fully into memory as contiguous per-contig
//! bases (line breaks removed, case preserved).

use std::fs;
use std::io;
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FaiEntry {
    pub name: String,
    pub len: u64,
    pub offset: u64,
    pub line_bases: u64,
    pub line_width: u64,
}

fn invalid(msg: String) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg)
}

pub fn read_fai(path: &Path) -> io::Result<Vec<FaiEntry>> {
    let text = fs::read_to_string(path)?;
    text.lines()
        .filter(|l| !l.is_empty())
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            let num = |i: usize| -> io::Result<u64> {
                f.get(i)
                    .and_then(|s| s.parse().ok())
                    .ok_or_else(|| invalid(format!("bad .fai line: {l}")))
            };
            let e = FaiEntry {
                name: f[0].to_owned(),
                len: num(1)?,
                offset: num(2)?,
                line_bases: num(3)?,
                line_width: num(4)?,
            };
            if e.line_bases == 0 || e.line_width < e.line_bases {
                return Err(invalid(format!("bad .fai geometry: {l}")));
            }
            Ok(e)
        })
        .collect()
}

pub struct Fasta {
    pub index: Vec<FaiEntry>,
    seqs: Vec<Vec<u8>>,
}

impl Fasta {
    /// Load `fasta` using `<fasta>.fai`, extracting contigs on `threads` threads.
    pub fn load(fasta: &Path, threads: usize) -> io::Result<Fasta> {
        let mut fai = fasta.as_os_str().to_owned();
        fai.push(".fai");
        let index = read_fai(Path::new(&fai))?;
        let data = fs::read(fasta)?;
        for e in &index {
            let full_lines = e.len / e.line_bases;
            let rem = e.len % e.line_bases;
            let span = full_lines * e.line_width + rem;
            if e.offset + span > data.len() as u64 {
                return Err(invalid(format!(
                    "{} extends past end of {}",
                    e.name,
                    fasta.display()
                )));
            }
        }
        let next = AtomicUsize::new(0);
        let mut seqs: Vec<Vec<u8>> = vec![Vec::new(); index.len()];
        let slots: Vec<std::sync::Mutex<&mut Vec<u8>>> =
            seqs.iter_mut().map(std::sync::Mutex::new).collect();
        std::thread::scope(|scope| {
            for _ in 0..threads.max(1) {
                scope.spawn(|| {
                    loop {
                        let i = next.fetch_add(1, Ordering::Relaxed);
                        let Some(e) = index.get(i) else { break };
                        let mut seq = Vec::with_capacity(e.len as usize);
                        let mut remaining = e.len;
                        let mut at = e.offset as usize;
                        while remaining > 0 {
                            let n = remaining.min(e.line_bases) as usize;
                            seq.extend_from_slice(&data[at..at + n]);
                            remaining -= n as u64;
                            at += e.line_width as usize;
                        }
                        **slots[i].lock().unwrap() = seq;
                    }
                });
            }
        });
        drop(slots);
        Ok(Fasta { index, seqs })
    }

    pub fn len(&self) -> usize {
        self.seqs.len()
    }

    /// Upper-case every base and turn anything that is not `ACGTN` into `N`,
    /// as htsjdk's reader does with `preserveCase=false, preserveIUPAC=false`.
    pub fn normalize_upper_acgtn(&mut self) {
        std::thread::scope(|scope| {
            for seq in self.seqs.iter_mut() {
                scope.spawn(move || {
                    for b in seq.iter_mut() {
                        *b = match b.to_ascii_uppercase() {
                            x @ (b'A' | b'C' | b'G' | b'T' | b'N') => x,
                            _ => b'N',
                        };
                    }
                });
            }
        });
    }

    pub fn is_empty(&self) -> bool {
        self.seqs.is_empty()
    }

    /// Bases of contig `i` (0-based positions).
    pub fn contig(&self, i: usize) -> &[u8] {
        &self.seqs[i]
    }

    pub fn index_of(&self, name: &str) -> Option<usize> {
        self.index.iter().position(|e| e.name == name)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn loads_wrapped_contigs() {
        let dir = std::env::temp_dir().join(format!("nomics-fa-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let fa = dir.join("r.fa");
        // chr1: 10 bases at 4 per line; chr2: 8 bases at 4 per line (exact multiple); chrM: 3 bases.
        let text = ">chr1 desc\nACGT\nacgt\nNN\n>chr2\nTTTT\nGGGG\n>chrM\nCAT\n";
        std::fs::write(&fa, text).unwrap();
        let fai = "chr1\t10\t11\t4\t5\nchr2\t8\t30\t4\t5\nchrM\t3\t46\t4\t5\n";
        std::fs::write(dir.join("r.fa.fai"), fai).unwrap();
        let f = Fasta::load(&fa, 3).unwrap();
        std::fs::remove_dir_all(&dir).unwrap();
        assert_eq!(f.len(), 3);
        assert_eq!(f.contig(0), b"ACGTacgtNN");
        assert_eq!(f.contig(1), b"TTTTGGGG");
        assert_eq!(f.contig(2), b"CAT");
        assert_eq!(f.index_of("chrM"), Some(2));
    }
}
