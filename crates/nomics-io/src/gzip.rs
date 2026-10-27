//! gzip (RFC 1952) decoding: single- or multi-member files, which includes
//! BGZF. Members are verified against their CRC32 and ISIZE.
//!
//! [`gunzip_to`] streams decoded bytes to a callback with bounded memory.
//! [`GzReader`] runs that on a background thread and exposes a `Read`, so two
//! FASTQ files decompress in parallel while their records are paired up.
//! [`gunzip_to`] takes the whole compressed input as a slice; [`GzReader`]
//! maps the file rather than reading it, so that slice is page cache.

use std::fs;
use std::io::{self, Read};
use std::path::Path;
use std::sync::mpsc::{Receiver, SyncSender, sync_channel};
use std::thread::JoinHandle;

use crate::crc32::Crc32;
use crate::inflate::{ByteSink, inflate_to};

fn invalid(msg: String) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg)
}

/// Length of the member header starting at `h`.
fn header_len(h: &[u8], offset: usize) -> io::Result<usize> {
    let bad = || invalid(format!("bad gzip header at byte {offset}"));
    if h.len() < 10 || h[0] != 0x1f || h[1] != 0x8b || h[2] != 8 || h[3] & 0xe0 != 0 {
        return Err(bad());
    }
    let flg = h[3];
    let mut at = 10;
    if flg & 4 != 0 {
        let xlen =
            u16::from_le_bytes(h.get(at..at + 2).ok_or_else(bad)?.try_into().unwrap()) as usize;
        at += 2 + xlen;
    }
    for bit in [8u8, 16] {
        if flg & bit != 0 {
            let nul = h
                .get(at..)
                .and_then(|r| r.iter().position(|&c| c == 0))
                .ok_or_else(bad)?;
            at += nul + 1;
        }
    }
    if flg & 2 != 0 {
        at += 2;
    }
    if at > h.len() { Err(bad()) } else { Ok(at) }
}

/// Decode every gzip member in `data`, passing decoded bytes to `sink`.
/// Returns the total decoded length.
pub fn gunzip_to(data: &[u8], sink: ByteSink) -> io::Result<u64> {
    let mut at = 0;
    let mut total = 0u64;
    while at < data.len() {
        let start = at + header_len(&data[at..], at)?;
        let mut crc = Crc32::new();
        let mut size = 0u64;
        let used = inflate_to(&data[start..], &mut |chunk| {
            crc = crc.update(chunk);
            size += chunk.len() as u64;
            sink(chunk)
        })?;
        let trailer = data
            .get(start + used..start + used + 8)
            .ok_or_else(|| invalid(format!("truncated gzip trailer at byte {}", start + used)))?;
        let want_crc = u32::from_le_bytes(trailer[..4].try_into().unwrap());
        let want_size = u32::from_le_bytes(trailer[4..].try_into().unwrap());
        if crc.finish() != want_crc || size as u32 != want_size {
            return Err(invalid(format!(
                "gzip member at byte {at}: CRC32/ISIZE mismatch"
            )));
        }
        total += size;
        at = start + used + 8;
    }
    Ok(total)
}

/// The first decoded bytes of a gzip file (at most `limit`, fewer if the
/// first member is shorter), reading only the start of the file.
pub fn head(path: impl AsRef<Path>, limit: usize) -> io::Result<Vec<u8>> {
    let mut data = Vec::new();
    fs::File::open(path.as_ref())?
        .take(4 << 20)
        .read_to_end(&mut data)?;
    if data.is_empty() {
        return Ok(Vec::new());
    }
    let start = header_len(&data, 0)?;
    let mut out = Vec::new();
    match crate::inflate::inflate(&data[start..], &mut out, limit) {
        Ok(_)
        | Err(crate::inflate::InflateError::OutputTooLarge)
        | Err(crate::inflate::InflateError::Truncated) => Ok(out),
        Err(e) => Err(invalid(format!("{}: {e}", path.as_ref().display()))),
    }
}

/// A gzip file decoded on a background thread, read through `Read`.
pub struct GzReader {
    rx: Receiver<Vec<u8>>,
    chunk: Vec<u8>,
    pos: usize,
    worker: Option<JoinHandle<io::Result<u64>>>,
}

impl GzReader {
    pub fn open(path: impl AsRef<Path>) -> io::Result<GzReader> {
        // Mapped, not read: the compressed file stays reclaimable page cache
        // (a 30x FASTQ pair is ~60 GB).
        let data = nomics_mem::MappedFile::open(path.as_ref())?;
        let (tx, rx): (SyncSender<Vec<u8>>, _) = sync_channel(8);
        let worker = std::thread::spawn(move || {
            gunzip_to(&data, &mut |chunk| {
                // A closed receiver means the reader was dropped early; stop quietly.
                tx.send(chunk.to_vec())
                    .map_err(|_| io::Error::new(io::ErrorKind::BrokenPipe, "reader dropped"))
            })
        });
        Ok(GzReader {
            rx,
            chunk: Vec::new(),
            pos: 0,
            worker: Some(worker),
        })
    }

    /// Surface a decode error from the worker once the channel has drained.
    fn finish(&mut self) -> io::Result<()> {
        match self.worker.take() {
            Some(w) => w
                .join()
                .map_err(|_| io::Error::other("gzip worker panicked"))?
                .map(|_| ()),
            None => Ok(()),
        }
    }
}

impl Read for GzReader {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        while self.pos == self.chunk.len() {
            match self.rx.recv() {
                Ok(c) => {
                    self.chunk = c;
                    self.pos = 0;
                }
                Err(_) => {
                    self.finish()?;
                    return Ok(0);
                }
            }
        }
        let n = buf.len().min(self.chunk.len() - self.pos);
        buf[..n].copy_from_slice(&self.chunk[self.pos..self.pos + n]);
        self.pos += n;
        Ok(n)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bgzf::{BgzfWriter, Profile};
    use crate::deflate::HtsjdkDeflater;
    use std::io::Write;

    // Python: gzip.compress(b"hello gzip\n", mtime=0) (one member, no FNAME).
    const HELLO: [u8; 31] = [
        0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x02, 0x13, 0xcb, 0x48, 0xcd, 0xc9, 0xc9,
        0x57, 0x48, 0xaf, 0xca, 0x2c, 0xe0, 0x02, 0x00, 0x39, 0x7c, 0x63, 0x56, 0x0b, 0x00, 0x00,
        0x00,
    ];

    fn collect(data: &[u8]) -> io::Result<Vec<u8>> {
        let mut out = Vec::new();
        gunzip_to(data, &mut |c| {
            out.extend_from_slice(c);
            Ok(())
        })?;
        Ok(out)
    }

    #[test]
    fn single_member() {
        assert_eq!(collect(&HELLO).unwrap(), b"hello gzip\n");
    }

    #[test]
    fn concatenated_members_and_bgzf() {
        let mut two = HELLO.to_vec();
        two.extend_from_slice(&HELLO);
        assert_eq!(collect(&two).unwrap(), b"hello gzip\nhello gzip\n");

        let payload: Vec<u8> = (0..200_000u32).map(|i| (i % 251) as u8).collect();
        let mut w = BgzfWriter::new(Vec::new(), HtsjdkDeflater::new(), Profile::HTSJDK);
        w.write_all(&payload).unwrap();
        assert_eq!(collect(&w.finish().unwrap()).unwrap(), payload);
    }

    #[test]
    fn bad_crc_is_rejected() {
        let mut bad = HELLO;
        bad[23] ^= 1;
        assert!(collect(&bad).is_err());
    }

    #[test]
    fn reader_thread_streams_a_file() {
        let dir = std::env::temp_dir().join(format!("nomics-gz-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("x.gz");
        let payload: Vec<u8> = (0..3_000_000u32)
            .map(|i| (i.wrapping_mul(7919) >> 5) as u8)
            .collect();
        let mut w = BgzfWriter::new(Vec::new(), HtsjdkDeflater::new(), Profile::HTSJDK);
        w.write_all(&payload).unwrap();
        std::fs::write(&path, w.finish().unwrap()).unwrap();
        let mut got = Vec::new();
        GzReader::open(&path)
            .unwrap()
            .read_to_end(&mut got)
            .unwrap();
        std::fs::remove_dir_all(&dir).unwrap();
        assert!(got == payload);
    }
}
