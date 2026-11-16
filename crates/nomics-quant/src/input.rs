//! Bounded, allocation-reusing text scans for large transcriptomes/annotations.
use std::collections::HashSet;
use std::io::{self, BufRead, BufReader, Read};
use std::path::Path;

pub(crate) fn reader(path: &Path) -> io::Result<Box<dyn BufRead + Send>> {
    let mut file = std::fs::File::open(path)?;
    let mut magic = [0; 2];
    let n = file.read(&mut magic)?;
    if n == 2 && magic == [31, 139] {
        Ok(Box::new(BufReader::with_capacity(
            1 << 20,
            nomics_io::gzip::GzReader::open(path)?,
        )))
    } else {
        Ok(Box::new(BufReader::with_capacity(
            1 << 20,
            std::fs::File::open(path)?,
        )))
    }
}

pub(crate) fn validate_fasta(path: &Path, header: &nomics_io::bam::Header) -> io::Result<()> {
    let mut missing: HashSet<&[u8]> = header.refs.iter().map(|r| r.name.as_slice()).collect();
    let mut r = reader(path)?;
    let mut line = Vec::new();
    loop {
        let buffer = r.fill_buf()?;
        if buffer.is_empty() {
            break;
        }
        if buffer[0] != b'>' {
            r.skip_until(b'\n')?;
            continue;
        }
        line.clear();
        r.read_until(b'\n', &mut line)?;
        let name = line[1..]
            .split(|b| b.is_ascii_whitespace())
            .next()
            .unwrap_or_default();
        missing.remove(name);
    }
    // Diagnose in @SQ order, never HashSet iteration order.
    for reference in &header.refs {
        if missing.contains(reference.name.as_slice()) {
            return Err(crate::invalid(format!(
                "BAM reference {} absent from FASTA",
                String::from_utf8_lossy(&reference.name)
            )));
        }
    }
    Ok(())
}
