//! Minimal checked ustar transport. No extraction paths are written to disk.
use nomics_io::{crc32::Crc32, deflate::ZlibDeflater, gzip::GzReader};
use nomics_mem::HostBuf;
use std::io::{self, Read};
use std::path::Path;

pub(super) fn invalid(message: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message.into())
}
fn text(bytes: &[u8]) -> io::Result<&str> {
    std::str::from_utf8(bytes.split(|&b| b == 0).next().unwrap())
        .map_err(|e| invalid(e.to_string()))
}
fn octal(bytes: &[u8]) -> io::Result<usize> {
    usize::from_str_radix(text(bytes)?.trim(), 8).map_err(|e| invalid(e.to_string()))
}
pub fn read(
    path: &Path,
    mut visit: impl FnMut(&str, &[u8]) -> io::Result<()>,
) -> io::Result<()> {
    let mut reader = GzReader::open(path)?;
    let mut names = std::collections::HashSet::new();
    loop {
        let mut h = [0u8; 512];
        reader.read_exact(&mut h)?;
        if h.iter().all(|&b| b == 0) {
            break;
        }
        let checksum: usize = h
            .iter()
            .enumerate()
            .map(|(i, &b)| {
                if (148..156).contains(&i) {
                    32
                } else {
                    b as usize
                }
            })
            .sum();
        if checksum != octal(&h[148..156])? {
            return Err(invalid("tar checksum mismatch"));
        }
        if &h[257..262] != b"ustar" {
            return Err(invalid("expected ustar header"));
        }
        let name = text(&h[..100])?;
        let prefix = text(&h[345..500])?;
        let name = if prefix.is_empty() {
            name.to_owned()
        } else {
            format!("{prefix}/{name}")
        };
        let name = name.strip_prefix("./").unwrap_or(&name);
        if name.contains('/') || name == ".." {
            return Err(invalid("only top-level tar leaves are supported"));
        }
        if !names.insert(name.to_owned()) {
            return Err(invalid("duplicate tar leaf"));
        }
        if !matches!(h[156], 0 | b'0') {
            return Err(invalid("only regular tar members are supported"));
        }
        let size = octal(&h[124..136])?;
        let mut data = HostBuf::new(size)?;
        reader.read_exact(data.as_mut_slice())?;
        visit(name, data.as_slice())?;
        let padding = (512 - size % 512) % 512;
        reader.read_exact(&mut h[..padding])?;
    }
    // Drain the gzip reader, validating CRC/ISIZE even beyond the tar terminator.
    io::copy(&mut reader, &mut io::sink())?;
    Ok(())
}
fn put_octal(field: &mut [u8], value: usize) -> io::Result<()> {
    let s = format!("{value:0width$o}", width = field.len() - 1);
    if s.len() >= field.len() {
        return Err(invalid("ustar field overflow"));
    }
    field[..s.len()].copy_from_slice(s.as_bytes());
    Ok(())
}
pub fn write(path: &Path, members: &[(String, Vec<u8>)]) -> io::Result<()> {
    let mut tar = Vec::new();
    for (name, data) in members {
        if name.len() > 100 || !name.is_ascii() {
            return Err(invalid("ustar name exceeds 100 ASCII bytes"));
        }
        let mut h = [0u8; 512];
        h[..name.len()].copy_from_slice(name.as_bytes());
        put_octal(&mut h[100..108], 0o644)?;
        put_octal(&mut h[108..116], 0)?;
        put_octal(&mut h[116..124], 0)?;
        put_octal(&mut h[124..136], data.len())?;
        put_octal(&mut h[136..148], 0)?;
        h[148..156].fill(b' ');
        h[156] = b'0';
        h[257..263].copy_from_slice(b"ustar\0");
        h[263..265].copy_from_slice(b"00");
        put_octal(&mut h[329..337], 0)?;
        put_octal(&mut h[337..345], 0)?;
        let sum: usize = h.iter().map(|&x| x as usize).sum();
        h[148..154].copy_from_slice(format!("{sum:06o}").as_bytes());
        h[154] = 0;
        h[155] = b' ';
        tar.extend_from_slice(&h);
        tar.extend_from_slice(data);
        tar.resize(tar.len().next_multiple_of(512), 0);
    }
    tar.resize(tar.len() + 1024, 0);
    let mut gzip = vec![31, 139, 8, 0, 0, 0, 0, 0, 0, 255];
    ZlibDeflater::new(2).compress(&tar, &mut gzip);
    gzip.extend_from_slice(&Crc32::new().update(&tar).finish().to_le_bytes());
    gzip.extend_from_slice(&(tar.len() as u32).to_le_bytes());
    std::fs::write(path, gzip)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn ustar_roundtrip_and_gzip_integrity() {
        let path = std::env::temp_dir().join(format!(
            "nomics-orientation-tar-{}-{}.tar.gz",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let members = vec![(
            "./sample.orientation_priors".to_owned(),
            b"prior table\n".to_vec(),
        )];
        write(&path, &members).unwrap();
        let mut got = Vec::new();
        read(&path, |n, b| {
            got.push((n.to_owned(), b.to_vec()));
            Ok(())
        })
        .unwrap();
        assert_eq!(
            got,
            vec![(
                "sample.orientation_priors".to_owned(),
                b"prior table\n".to_vec()
            )]
        );
        let mut bytes = std::fs::read(&path).unwrap();
        let last = bytes.len() - 8;
        bytes[last] ^= 1;
        std::fs::write(&path, bytes).unwrap();
        assert!(read(&path, |_, _| Ok(())).is_err());
        std::fs::remove_file(path).unwrap();
    }
}
