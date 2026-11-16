//! Exact k-mer acceleration table for BWA bi-interval extension.
//!
//! A full-k endpoint alone is not sufficient for `bwt_smem1`: its forward
//! pass saves the interval immediately before every occurrence-count change,
//! and those saved intervals are the starting set for the backward pass.  The
//! table therefore contains the endpoint for every word of every length from
//! 2 through `k`.  Looking up those prefixes reproduces the same change events
//! without touching the BWT occurrence blocks.

use std::fs::{File, OpenOptions};
use std::io::{self, Read, Seek, Write};
use std::path::{Path, PathBuf};

use nomics_mem::HostBuf;

use crate::index::Bwt;
use crate::seed::{Intv, extend};

const MAGIC: &[u8; 8] = b"NMKMER01";
const VERSION: u32 = 1;
pub const MIN_K: usize = 2;
pub const MAX_K: usize = 13;
const HEADER_LEN: u64 = 96;

/// Endpoints for all DNA words of lengths 2 through `k`, in length-major and
/// then base-4 lexical order. Each endpoint is three native little-endian
/// `u64`s and the whole payload is huge-page-backed for direct GPU reads.
pub struct KmerTable {
    k: usize,
    offsets: [u64; MAX_K + 1],
    intervals: HostBuf,
}

impl KmerTable {
    /// Load a cache file, or build it atomically when it does not exist.
    pub fn load_or_build(path: &Path, bwt: &Bwt, k: usize) -> io::Result<Self> {
        validate_k(k)?;
        match Self::load(path, bwt, k) {
            Ok(table) => Ok(table),
            Err(e) if e.kind() == io::ErrorKind::NotFound => {
                let table = Self::build(bwt, k)?;
                table.store(path, bwt)?;
                Ok(table)
            }
            Err(e) => Err(e),
        }
    }

    pub fn load(path: &Path, bwt: &Bwt, k: usize) -> io::Result<Self> {
        validate_k(k)?;
        let mut file = File::open(path)?;
        let mut header = [0u8; HEADER_LEN as usize];
        file.read_exact(&mut header)?;
        let mut at = 0usize;
        let magic = take::<8>(&header, &mut at);
        let version = u32::from_le_bytes(take::<4>(&header, &mut at));
        let file_k = u32::from_le_bytes(take::<4>(&header, &mut at)) as usize;
        let primary = u64::from_le_bytes(take::<8>(&header, &mut at));
        let seq_len = u64::from_le_bytes(take::<8>(&header, &mut at));
        let mut l2 = [0u64; 5];
        for value in &mut l2 {
            *value = u64::from_le_bytes(take::<8>(&header, &mut at));
        }
        let words_len = u64::from_le_bytes(take::<8>(&header, &mut at));
        let fingerprint = u64::from_le_bytes(take::<8>(&header, &mut at));
        let entries = u64::from_le_bytes(take::<8>(&header, &mut at));
        let expected_entries = entry_count(k)?;
        let expected_size = HEADER_LEN
            .checked_add(expected_entries.checked_mul(24).ok_or_else(too_large)?)
            .ok_or_else(too_large)?;
        let metadata_len = file.metadata()?.len();
        if &magic != MAGIC
            || version != VERSION
            || file_k != k
            || primary != bwt.primary
            || seq_len != bwt.seq_len
            || l2 != bwt.l2
            || words_len != bwt.words().len() as u64
            || fingerprint != bwt_fingerprint(bwt)
            || entries != expected_entries
            || metadata_len != expected_size
        {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("{} is not the requested BWT k-mer table", path.display()),
            ));
        }
        let bytes = usize::try_from(expected_entries.checked_mul(24).ok_or_else(too_large)?)
            .map_err(|_| too_large())?;
        let mut intervals = HostBuf::new(bytes)?;
        intervals.set_label("k-mer intervals");
        file.read_exact(intervals.as_mut_slice())?;
        Ok(Self {
            k,
            offsets: level_offsets(k)?,
            intervals,
        })
    }

    fn build(bwt: &Bwt, k: usize) -> io::Result<Self> {
        validate_k(k)?;
        let offsets = level_offsets(k)?;
        let entries = entry_count(k)?;
        let bytes = usize::try_from(entries.checked_mul(24).ok_or_else(too_large)?)
            .map_err(|_| too_large())?;
        let mut intervals = HostBuf::new(bytes)?;
        intervals.set_label("k-mer intervals");
        let words = intervals.as_u64s_mut();
        for len in MIN_K..=k {
            let count = 1u64 << (2 * len);
            let current_word = usize::try_from(offsets[len] * 3).map_err(|_| too_large())?;
            let (earlier, current_and_later) = words.split_at_mut(current_word);
            let current =
                &mut current_and_later[..usize::try_from(count * 3).map_err(|_| too_large())?];
            for code in 0..count {
                let prefix = code >> 2;
                let ik = if len == MIN_K {
                    set_single(bwt, prefix as u8)
                } else {
                    read_intv(earlier, offsets[len - 1] + prefix)
                };
                let next = extend(bwt, &ik, false)[(3 - (code as u8 & 3)) as usize];
                let out = &mut current[code as usize * 3..][..3];
                out.copy_from_slice(&next.x);
            }
        }
        Ok(Self {
            k,
            offsets,
            intervals,
        })
    }

    fn store(&self, path: &Path, bwt: &Bwt) -> io::Result<()> {
        if !cfg!(target_endian = "little") {
            return Err(io::Error::new(
                io::ErrorKind::Unsupported,
                "k-mer cache format requires a little-endian host",
            ));
        }
        if let Some(parent) = path.parent()
            && !parent.as_os_str().is_empty()
        {
            std::fs::create_dir_all(parent)?;
        }
        let tmp = temporary_path(path);
        let mut file = OpenOptions::new()
            .create(true)
            .truncate(true)
            .write(true)
            .open(&tmp)?;
        file.write_all(MAGIC)?;
        file.write_all(&VERSION.to_le_bytes())?;
        file.write_all(&(self.k as u32).to_le_bytes())?;
        file.write_all(&bwt.primary.to_le_bytes())?;
        file.write_all(&bwt.seq_len.to_le_bytes())?;
        for value in bwt.l2 {
            file.write_all(&value.to_le_bytes())?;
        }
        file.write_all(&(bwt.words().len() as u64).to_le_bytes())?;
        file.write_all(&bwt_fingerprint(bwt).to_le_bytes())?;
        file.write_all(&entry_count(self.k)?.to_le_bytes())?;
        debug_assert_eq!(file.stream_position()?, HEADER_LEN);
        file.write_all(self.intervals.as_slice())?;
        file.sync_all()?;
        std::fs::rename(&tmp, path)?;
        Ok(())
    }

    #[inline]
    pub fn k(&self) -> usize {
        self.k
    }

    pub fn bytes(&self) -> usize {
        self.intervals.len()
    }

    /// Endpoint for the base-4 encoded word `code` of the given length.
    #[inline]
    pub(crate) fn interval(&self, len: usize, code: u64) -> Intv {
        debug_assert!((MIN_K..=self.k).contains(&len));
        debug_assert!(code < 1u64 << (2 * len));
        read_intv(self.intervals.as_u64s(), self.offsets[len] + code)
    }

    /// Borrow the same huge-page-backed payload in the CUDA ABI form.
    pub fn gpu_table(&self) -> nomics_gpu::GpuKmerTable<'_> {
        nomics_gpu::GpuKmerTable {
            k: self.k as u32,
            level_offsets: self.offsets,
            intervals: self.intervals.as_u64s(),
        }
    }
}

fn validate_k(k: usize) -> io::Result<()> {
    if (MIN_K..=MAX_K).contains(&k) {
        Ok(())
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("k-mer table k must be in {MIN_K}..={MAX_K}, got {k}"),
        ))
    }
}

fn level_offsets(k: usize) -> io::Result<[u64; MAX_K + 1]> {
    let mut offsets = [0u64; MAX_K + 1];
    let mut entries = 0u64;
    for (len, offset) in offsets.iter_mut().enumerate().take(k + 1).skip(MIN_K) {
        *offset = entries;
        entries = entries
            .checked_add(1u64 << (2 * len))
            .ok_or_else(too_large)?;
    }
    Ok(offsets)
}

fn entry_count(k: usize) -> io::Result<u64> {
    let offsets = level_offsets(k)?;
    offsets[k]
        .checked_add(1u64 << (2 * k))
        .ok_or_else(too_large)
}

#[inline]
fn read_intv(words: &[u64], entry: u64) -> Intv {
    let at = entry as usize * 3;
    Intv {
        x: [words[at], words[at + 1], words[at + 2]],
        info: 0,
    }
}

#[inline]
fn set_single(bwt: &Bwt, c: u8) -> Intv {
    let c = c as usize;
    Intv {
        x: [bwt.l2[c] + 1, bwt.l2[3 - c] + 1, bwt.l2[c + 1] - bwt.l2[c]],
        info: 0,
    }
}

fn bwt_fingerprint(bwt: &Bwt) -> u64 {
    // FNV-1a is an identity guard against accidentally reusing another
    // reference's cache, not a security boundary.
    let mut hash = 0xcbf2_9ce4_8422_2325u64;
    for &word in bwt.words() {
        hash ^= word as u64;
        hash = hash.wrapping_mul(0x100_0000_01b3);
    }
    hash
}

fn temporary_path(path: &Path) -> PathBuf {
    let mut name = path.as_os_str().to_owned();
    name.push(format!(".tmp.{}", std::process::id()));
    PathBuf::from(name)
}

fn too_large() -> io::Error {
    io::Error::new(io::ErrorKind::InvalidInput, "k-mer table size overflow")
}

fn take<const N: usize>(bytes: &[u8], at: &mut usize) -> [u8; N] {
    let out = bytes[*at..*at + N].try_into().unwrap();
    *at += N;
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exact_table_sizes() {
        assert_eq!(entry_count(10).unwrap() * 24, 32 * 1024 * 1024 - 128);
        assert_eq!(entry_count(11).unwrap() * 24, 128 * 1024 * 1024 - 128);
        assert_eq!(entry_count(12).unwrap() * 24, 512 * 1024 * 1024 - 128);
        assert_eq!(entry_count(13).unwrap() * 24, 2 * 1024 * 1024 * 1024 - 128);
    }
}
