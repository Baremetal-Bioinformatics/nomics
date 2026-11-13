//! A denser sampled suffix array derived exactly from bwa's sampled SA.
//!
//! The cache can sample every 1st, 2nd, or 4th BWT row. Values are stored as
//! packed 40-bit little-endian integers; all ones represents bwa's implicit
//! `SA[0] = -1`. Lookups LF-walk to the nearest cached row, so changing the
//! interval only changes the number of walks, never the returned position.
//!
//! Building each output entry independently would repeat about 15.5 LF steps
//! per entry with bwa's interval-32 SA. Instead, the builder starts at every
//! bwa sample and follows LF until the next bwa-sampled row. These paths
//! partition the complete LF cycle, so every BWT row is visited exactly once.
//! Path discovery is parallel; a writer packs the resulting random entries
//! into the single `HostBuf` without unsafe shared mutation.

use std::fs::File;
use std::io::{self, Read, Write};
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::mpsc;

use nomics_mem::HostBuf;

use crate::index::{Bwt, NONE, Sa};

const LEGACY_MAGIC: &[u8; 8] = b"NOMSA4v1";
const MAGIC: &[u8; 8] = b"NOMSAIv2";
const LEGACY_HEADER_BYTES: u64 = 40;
const HEADER_BYTES: u64 = 48;
const BYTES: usize = 5;
const ALL_ONES: u64 = (1 << 40) - 1;
const SOURCE_TASK: usize = 4096;
const RECORD_BATCH: usize = 1 << 15;

pub struct DenseSa {
    data: HostBuf,
    n: u64,
    interval: u64,
}

fn invalid(msg: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg.into())
}

fn validate_interval(interval: u64) -> io::Result<()> {
    if matches!(interval, 1 | 2 | 4) {
        Ok(())
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("dense SA interval must be 1, 2, or 4, not {interval}"),
        ))
    }
}

fn entry_count(seq_len: u64, interval: u64) -> io::Result<u64> {
    seq_len
        .checked_div(interval)
        .and_then(|n| n.checked_add(1))
        .ok_or_else(|| invalid("dense SA entry count overflow"))
}

fn data_bytes(n: u64) -> io::Result<usize> {
    usize::try_from(n)
        .ok()
        .and_then(|n| n.checked_mul(BYTES))
        .ok_or_else(|| invalid("dense SA is too large for this platform"))
}

/// Identifies the source index: primary, length and a hash of bwa's samples.
fn fingerprint(bwt: &Bwt, sa: &Sa) -> [u64; 3] {
    let mut h = 0xcbf2_9ce4_8422_2325u64;
    for &v in sa.entries() {
        h = (h ^ v).wrapping_mul(0x0000_0100_0000_01b3);
    }
    [bwt.primary, bwt.seq_len, h]
}

#[inline]
fn preceding_sa(value: u64, seq_len: u64) -> u64 {
    match value {
        NONE => seq_len - 1,
        0 => NONE,
        _ => value - 1,
    }
}

#[inline]
fn put_entry(out: &mut [u8], i: u64, value: u64) {
    let at = i as usize * BYTES;
    out[at..at + BYTES].copy_from_slice(&(value & ALL_ONES).to_le_bytes()[..BYTES]);
}

#[inline]
fn get_entry(data: &[u8], i: u64) -> u64 {
    let at = i as usize * BYTES;
    let b = &data[at..at + BYTES];
    let value = u64::from_le_bytes([b[0], b[1], b[2], b[3], b[4], 0, 0, 0]);
    if value == ALL_ONES { NONE } else { value }
}

impl DenseSa {
    /// Invert the LF cycle from bwa's samples on `threads` worker threads.
    pub fn build(bwt: &Bwt, sa: &Sa, interval: u64, threads: usize) -> io::Result<Self> {
        validate_interval(interval)?;
        if bwt.seq_len == 0 {
            return Err(invalid("cannot build a dense SA for an empty BWT"));
        }
        if bwt.seq_len > ALL_ONES {
            return Err(invalid("BWT positions do not fit in 40 bits"));
        }
        let n = entry_count(bwt.seq_len, interval)?;
        let mut data = HostBuf::new(data_bytes(n)?)?;
        data.set_label("dense SA");
        let workers = threads.max(1);
        let next_source = AtomicUsize::new(0);
        let source_count = sa.entries().len();
        let (send, recv) = mpsc::sync_channel::<Vec<(u64, u64)>>(workers * 2);

        let (visited, written) = std::thread::scope(|scope| -> io::Result<(u64, u64)> {
            let out = data.as_mut_slice();
            let writer = scope.spawn(move || {
                let mut written = 0u64;
                while let Ok(records) = recv.recv() {
                    written += records.len() as u64;
                    for (i, value) in records {
                        put_entry(out, i, value);
                    }
                }
                written
            });

            let mut handles = Vec::with_capacity(workers);
            for _ in 0..workers {
                let send = send.clone();
                let next_source = &next_source;
                handles.push(scope.spawn(move || -> io::Result<u64> {
                    let mut records = Vec::with_capacity(RECORD_BATCH);
                    let mut visited = 0u64;
                    loop {
                        let first = next_source.fetch_add(SOURCE_TASK, Ordering::Relaxed);
                        if first >= source_count {
                            break;
                        }
                        let end = (first + SOURCE_TASK).min(source_count);
                        for source in first..end {
                            let mut row = source as u64 * sa.intv;
                            let mut value = sa.entries()[source];
                            let mut path_steps = 0u64;
                            loop {
                                visited += 1;
                                path_steps += 1;
                                if row % interval == 0 {
                                    records.push((row / interval, value));
                                    if records.len() == RECORD_BATCH {
                                        send.send(std::mem::take(&mut records)).map_err(|_| {
                                            io::Error::new(
                                                io::ErrorKind::BrokenPipe,
                                                "dense SA writer stopped",
                                            )
                                        })?;
                                        records = Vec::with_capacity(RECORD_BATCH);
                                    }
                                }
                                row = bwt.inv_psi(row);
                                value = preceding_sa(value, bwt.seq_len);
                                if row % sa.intv == 0 {
                                    if value != sa.entries()[(row / sa.intv) as usize] {
                                        return Err(invalid(format!(
                                            "inconsistent SA/LF cycle at sampled row {row}"
                                        )));
                                    }
                                    break;
                                }
                                if path_steps > bwt.seq_len + 1 {
                                    return Err(invalid("SA/LF walk did not reach a sampled row"));
                                }
                            }
                        }
                    }
                    if !records.is_empty() {
                        send.send(records).map_err(|_| {
                            io::Error::new(io::ErrorKind::BrokenPipe, "dense SA writer stopped")
                        })?;
                    }
                    Ok(visited)
                }));
            }
            drop(send);

            let mut visited = 0u64;
            for handle in handles {
                visited += handle
                    .join()
                    .map_err(|_| invalid("dense SA worker panicked"))??;
            }
            let written = writer
                .join()
                .map_err(|_| invalid("dense SA writer panicked"))?;
            Ok((visited, written))
        })?;

        if visited != bwt.seq_len + 1 || written != n {
            return Err(invalid(format!(
                "incomplete SA/LF inversion: visited {visited} of {} rows, wrote {written} of {n} entries",
                bwt.seq_len + 1
            )));
        }
        Ok(Self { data, n, interval })
    }

    /// Load the requested interval from `path`, or build and save it if absent.
    /// Version-1 interval-4 caches remain readable.
    pub fn load_or_build(
        path: &Path,
        bwt: &Bwt,
        sa: &Sa,
        interval: u64,
        threads: usize,
    ) -> io::Result<Self> {
        validate_interval(interval)?;
        match File::open(path) {
            Ok(mut file) => return Self::load(&mut file, path, bwt, sa, interval),
            Err(e) if e.kind() == io::ErrorKind::NotFound => {}
            Err(e) => return Err(e),
        }

        let dense = Self::build(bwt, sa, interval, threads)?;
        let tmp = path.with_extension(format!("tmp.{}", std::process::id()));
        {
            let mut file = io::BufWriter::with_capacity(1 << 24, File::create(&tmp)?);
            file.write_all(MAGIC)?;
            for word in [
                interval,
                bwt.primary,
                bwt.seq_len,
                fingerprint(bwt, sa)[2],
                dense.n,
            ] {
                file.write_all(&word.to_le_bytes())?;
            }
            file.write_all(dense.data.as_slice())?;
            file.flush()?;
        }
        std::fs::rename(&tmp, path)?;
        Ok(dense)
    }

    fn load(
        file: &mut File,
        path: &Path,
        bwt: &Bwt,
        sa: &Sa,
        requested_interval: u64,
    ) -> io::Result<Self> {
        let mut magic = [0u8; 8];
        file.read_exact(&mut magic)?;
        let (interval, header_bytes, words) = if &magic == LEGACY_MAGIC {
            let mut bytes = [0u8; 32];
            file.read_exact(&mut bytes)?;
            let words = bytes
                .chunks_exact(8)
                .map(|b| u64::from_le_bytes(b.try_into().unwrap()))
                .collect::<Vec<_>>();
            (4, LEGACY_HEADER_BYTES, words)
        } else if &magic == MAGIC {
            let mut bytes = [0u8; 40];
            file.read_exact(&mut bytes)?;
            let mut words = bytes
                .chunks_exact(8)
                .map(|b| u64::from_le_bytes(b.try_into().unwrap()));
            let interval = words.next().unwrap();
            (interval, HEADER_BYTES, words.collect::<Vec<_>>())
        } else {
            return Err(invalid(format!(
                "{}: unrecognized dense SA cache format",
                path.display()
            )));
        };

        validate_interval(interval).map_err(|e| invalid(format!("{}: {e}", path.display())))?;
        let fp = fingerprint(bwt, sa);
        let n = entry_count(bwt.seq_len, interval)?;
        let expected_size = header_bytes
            .checked_add(u64::try_from(data_bytes(n)?).unwrap())
            .ok_or_else(|| invalid("dense SA cache size overflow"))?;
        let mismatch = interval != requested_interval
            || words.as_slice() != [fp[0], fp[1], fp[2], n]
            || file.metadata()?.len() != expected_size;
        if mismatch {
            return Err(invalid(format!(
                "{}: dense SA cache does not match this index and interval {requested_interval}; use the matching cache path or delete it to rebuild",
                path.display()
            )));
        }
        let mut data = HostBuf::new(data_bytes(n)?)?;
        data.set_label("dense SA");
        file.read_exact(data.as_mut_slice())?;
        Ok(Self { data, n, interval })
    }

    /// `bwt_sa(k)` using this cache's sampling interval.
    #[inline]
    pub fn sa(&self, bwt: &Bwt, mut k: u64) -> u64 {
        let mut steps = 0u64;
        while k & (self.interval - 1) != 0 {
            steps += 1;
            k = bwt.inv_psi(k);
        }
        steps.wrapping_add(get_entry(self.data.as_slice(), k / self.interval))
    }

    pub fn interval(&self) -> u64 {
        self.interval
    }

    pub fn bytes(&self) -> usize {
        self.data.len()
    }

    pub fn huge_page_bytes(&self) -> Option<usize> {
        self.data.huge_page_bytes()
    }

    pub fn lock(&self) -> std::io::Result<()> {
        self.data.lock()
    }

    pub fn collapse_huge(&self) -> std::io::Result<()> {
        self.data.collapse_huge()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn packed_entries_round_trip() {
        let mut data = [0u8; 20];
        for (i, value) in [NONE, 0, 1, ALL_ONES - 1].into_iter().enumerate() {
            put_entry(&mut data, i as u64, value);
            assert_eq!(get_entry(&data, i as u64), value);
        }
    }

    #[test]
    fn predecessor_wraps_like_bwa_sa() {
        assert_eq!(preceding_sa(NONE, 12), 11);
        assert_eq!(preceding_sa(11, 12), 10);
        assert_eq!(preceding_sa(1, 12), 0);
        assert_eq!(preceding_sa(0, 12), NONE);
    }

    #[test]
    fn only_supported_intervals_are_accepted() {
        for interval in [1, 2, 4] {
            validate_interval(interval).unwrap();
        }
        for interval in [0, 3, 8] {
            assert!(validate_interval(interval).is_err());
        }
    }
}
