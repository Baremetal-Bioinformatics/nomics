//! Host memory shared in place by CPU and GPU.
//!
//! On GB10 the GPU reads ordinary host memory at near-device bandwidth, but only
//! when that memory is backed by 2 MiB pages. These rules were measured on the
//! Spark in an earlier project and are enforced here:
//!
//! - Allocate with anonymous `mmap`, 2 MiB aligned, and `madvise(MADV_HUGEPAGE)`
//!   *before* first touch (the system THP mode is `madvise`).
//! - Never `cudaHostRegister` (15x slower) and never `cudaMalloc` large buffers
//!   (it is system memory on GB10 and fragments into small pages).
//! - Allocate the largest long-lived buffer (the index) first, then verify its
//!   backing with [`HostBuf::huge_page_bytes`].
//!
//! This is one of two crates allowed to use `unsafe` (the other is `nomics-gpu`).
#![allow(unsafe_code)]

mod array;
pub mod diag;
pub use array::{HostArray, HostElement};
pub use diag::{CanaryDamage, check_all_canaries, diag};

use std::io;
use std::ptr::NonNull;

/// Huge-page size on GB10 / aarch64 Linux with 4 KiB base pages.
pub const HUGE_PAGE: usize = 2 << 20;

static SMALL_PAGE_BUFFERS: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);

/// GPU buffers that [`HostBuf::touch_huge`] found under 90% on 2 MiB pages so far.
pub fn small_page_buffers() -> usize {
    SMALL_PAGE_BUFFERS.load(std::sync::atomic::Ordering::Relaxed)
}

/// Bytes currently mapped by live [`HostBuf`]s (diagnostics).
static HOSTBUF_BYTES: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);

/// Bytes currently mapped by live [`HostBuf`]s, touched or not.
pub fn hostbuf_bytes() -> usize {
    HOSTBUF_BYTES.load(std::sync::atomic::Ordering::Relaxed)
}

/// This process's resident anonymous memory (`RssAnon`, Linux), in bytes:
/// the part comparable with system memory use, without mapped input files.
pub fn rss_anon_bytes() -> Option<u64> {
    let status = std::fs::read_to_string("/proc/self/status").ok()?;
    let line = status.lines().find(|l| l.starts_with("RssAnon:"))?;
    let kb: u64 = line[8..].trim().trim_end_matches("kB").trim().parse().ok()?;
    Some(kb * 1024)
}

/// Maximum resident set size of this process so far, in bytes.
pub fn peak_rss_bytes() -> Option<u64> {
    let mut usage = std::mem::MaybeUninit::<libc::rusage>::uninit();
    if unsafe { libc::getrusage(libc::RUSAGE_SELF, usage.as_mut_ptr()) } != 0 {
        return None;
    }
    let rss = unsafe { usage.assume_init() }.ru_maxrss as u64;
    #[cfg(target_os = "macos")]
    return Some(rss);
    #[cfg(not(target_os = "macos"))]
    return Some(rss * 1024);
}

/// A 2 MiB-aligned anonymous mapping, advised for transparent huge pages.
pub struct HostBuf {
    ptr: NonNull<u8>,
    len: usize,
    map_len: usize,
    map_base: *mut libc::c_void,
    /// Guarded layout ([`HostBuf::with_guards`]): the end of the accessible
    /// pages; 0 for the plain layout.
    guard_end: usize,
}

// The buffer is plain memory with a single owner; sharing follows &/&mut rules.
unsafe impl Send for HostBuf {}
unsafe impl Sync for HostBuf {}

impl HostBuf {
    /// Map `len` bytes (rounded up to a whole number of huge pages), zero-filled
    /// lazily by the kernel. Pages are not touched here, so the huge-page advice
    /// takes effect on first write. With `NOMICS_DIAG` set the buffer is
    /// mapped with guard pages and canaries ([`HostBuf::with_guards`]).
    pub fn new(len: usize) -> io::Result<Self> {
        if diag::guards() {
            return Self::with_guards(len);
        }
        let rounded = len.max(1).next_multiple_of(HUGE_PAGE);
        // Over-map by one huge page so a 2 MiB-aligned start always fits.
        let map_len = rounded + HUGE_PAGE;
        let base = unsafe {
            libc::mmap(
                std::ptr::null_mut(),
                map_len,
                libc::PROT_READ | libc::PROT_WRITE,
                libc::MAP_PRIVATE | libc::MAP_ANONYMOUS,
                -1,
                0,
            )
        };
        if base == libc::MAP_FAILED {
            return Err(io::Error::last_os_error());
        }
        let aligned = (base as usize).next_multiple_of(HUGE_PAGE);
        HOSTBUF_BYTES.fetch_add(map_len, std::sync::atomic::Ordering::Relaxed);
        #[cfg(target_os = "linux")]
        unsafe {
            if libc::madvise(aligned as *mut libc::c_void, rounded, libc::MADV_HUGEPAGE) != 0 {
                let err = io::Error::last_os_error();
                libc::munmap(base, map_len);
                HOSTBUF_BYTES.fetch_sub(map_len, std::sync::atomic::Ordering::Relaxed);
                return Err(err);
            }
        }
        Ok(Self {
            ptr: NonNull::new(aligned as *mut u8).expect("mmap returned null"),
            len,
            map_len,
            map_base: base,
            guard_end: 0,
        })
    }

    /// The diagnostic layout (`NOMICS_DIAG`, see [`diag`]): the same 2 MiB
    /// aligned start and `MADV_HUGEPAGE` data range, but only the pages
    /// holding `len` bytes are accessible. The slack from `len` to the end
    /// of its page and one page before the start are canaries; everything
    /// else in the mapping (at least one page on each side) is `PROT_NONE`.
    /// Only the last, partial 2 MiB of the data loses huge-page backing (its
    /// range is split by the tail guard). The buffer is registered for
    /// [`check_all_canaries`] and the fault handler.
    pub fn with_guards(len: usize) -> io::Result<Self> {
        let page = page_size();
        let rounded = len.max(1).next_multiple_of(HUGE_PAGE);
        // Two huge pages of slack: the aligned start always leaves at least
        // a canary and a guard page before it, and a guard page after.
        let map_len = rounded + 2 * HUGE_PAGE;
        let base = unsafe {
            libc::mmap(
                std::ptr::null_mut(),
                map_len,
                libc::PROT_READ | libc::PROT_WRITE,
                libc::MAP_PRIVATE | libc::MAP_ANONYMOUS,
                -1,
                0,
            )
        };
        if base == libc::MAP_FAILED {
            return Err(io::Error::last_os_error());
        }
        let b = base as usize;
        let aligned = (b + 2 * page).next_multiple_of(HUGE_PAGE);
        let end = aligned + len.next_multiple_of(page);
        debug_assert!(end + page <= b + map_len);
        let fail = |err: io::Error| {
            unsafe { libc::munmap(base, map_len) };
            Err(err)
        };
        #[cfg(target_os = "linux")]
        if end > aligned
            && unsafe {
                libc::madvise(
                    aligned as *mut libc::c_void,
                    end - aligned,
                    libc::MADV_HUGEPAGE,
                )
            } != 0
        {
            return fail(io::Error::last_os_error());
        }
        let none = |from: usize, to: usize| -> io::Result<()> {
            if to > from
                && unsafe { libc::mprotect(from as *mut libc::c_void, to - from, libc::PROT_NONE) }
                    != 0
            {
                return Err(io::Error::last_os_error());
            }
            Ok(())
        };
        if let Err(e) = none(b, aligned - page).and_then(|()| none(end, b + map_len)) {
            return fail(e);
        }
        // SAFETY: both ranges are read-write pages of this new mapping.
        unsafe {
            diag::fill(aligned - page, aligned);
            diag::fill(aligned + len, end);
        }
        HOSTBUF_BYTES.fetch_add(map_len, std::sync::atomic::Ordering::Relaxed);
        diag::register(&diag::Layout {
            map_base: b,
            map_len,
            data: aligned,
            len,
            end,
        });
        Ok(Self {
            ptr: NonNull::new(aligned as *mut u8).expect("mmap returned null"),
            len,
            map_len,
            map_base: base,
            guard_end: end,
        })
    }

    /// Whether this buffer has the diagnostic guard layout.
    pub fn is_guarded(&self) -> bool {
        self.guard_end != 0
    }

    /// Bytes from [`Self::as_ptr`] that may be accessed (e.g. wrapped for a
    /// GPU): whole 2 MiB pages normally, only the pages holding `len` with
    /// guards (the next page is `PROT_NONE`).
    pub fn accessible_len(&self) -> usize {
        if self.is_guarded() {
            self.guard_end - self.ptr.as_ptr() as usize
        } else {
            self.len.max(1).next_multiple_of(HUGE_PAGE)
        }
    }

    /// Name this buffer in diagnostic reports (guarded buffers only; a no-op
    /// otherwise). [`Self::touch_huge`] sets its label too.
    pub fn set_label(&self, label: &str) {
        if self.is_guarded() {
            diag::set_label(self.ptr.as_ptr() as usize, label);
        }
    }

    /// Check this buffer's canaries (guarded buffers; always `Ok` otherwise),
    /// naming it `name` in the damage report.
    pub fn check_canaries(&self, name: &str) -> Result<(), CanaryDamage> {
        if !self.is_guarded() {
            return Ok(());
        }
        diag::check_layout(name, self.ptr.as_ptr() as usize, self.len, self.guard_end)
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    pub fn as_ptr(&self) -> *const u8 {
        self.ptr.as_ptr()
    }

    pub fn as_mut_ptr(&mut self) -> *mut u8 {
        self.ptr.as_ptr()
    }

    pub fn as_slice(&self) -> &[u8] {
        unsafe { std::slice::from_raw_parts(self.ptr.as_ptr(), self.len) }
    }

    pub fn as_mut_slice(&mut self) -> &mut [u8] {
        unsafe { std::slice::from_raw_parts_mut(self.ptr.as_ptr(), self.len) }
    }

    /// View the buffer as `u64` words (the buffer is 2 MiB aligned).
    pub fn as_u64s(&self) -> &[u64] {
        unsafe { std::slice::from_raw_parts(self.ptr.as_ptr().cast(), self.len / 8) }
    }

    pub fn as_u64s_mut(&mut self) -> &mut [u64] {
        unsafe { std::slice::from_raw_parts_mut(self.ptr.as_ptr().cast(), self.len / 8) }
    }

    /// View the buffer as `u32` words (the buffer is 2 MiB aligned).
    pub fn as_u32s(&self) -> &[u32] {
        unsafe { std::slice::from_raw_parts(self.ptr.as_ptr().cast(), self.len / 4) }
    }

    pub fn as_u32s_mut(&mut self) -> &mut [u32] {
        unsafe { std::slice::from_raw_parts_mut(self.ptr.as_ptr().cast(), self.len / 4) }
    }

    /// Lock the buffer in RAM (`mlock`). The GPU reads some buffers in place
    /// without touching CPU page tables, so the kernel sees them as idle and
    /// may swap them out under page-cache pressure; they come back on 4 KiB
    /// pages, which makes GPU access several times slower. Locked pages are
    /// never swapped, so they keep their 2 MiB backing.
    pub fn lock(&self) -> io::Result<()> {
        // SAFETY: the range is this buffer's own live mapping.
        let rc = unsafe { libc::mlock(self.ptr.as_ptr() as *const libc::c_void, self.len) };
        if rc == 0 {
            Ok(())
        } else {
            Err(io::Error::last_os_error())
        }
    }

    /// Ask the kernel to back this buffer with 2 MiB pages now
    /// (`MADV_COLLAPSE`, Linux 6.1+): it compacts memory and promotes the
    /// range synchronously. `MADV_HUGEPAGE` before first touch is not enough
    /// when the page cache holds the free memory at load time; the faults then
    /// fall back to 4 KiB pages and khugepaged is too slow to repair a 20 GB
    /// index. A no-op error elsewhere; callers treat failure as advisory.
    pub fn collapse_huge(&self) -> io::Result<()> {
        #[cfg(target_os = "linux")]
        {
            const MADV_COLLAPSE: libc::c_int = 25;
            // SAFETY: the range is this buffer's own live mapping; collapse
            // preserves contents.
            let rc = unsafe {
                libc::madvise(
                    self.ptr.as_ptr() as *mut libc::c_void,
                    self.len,
                    MADV_COLLAPSE,
                )
            };
            if rc == 0 {
                Ok(())
            } else {
                Err(io::Error::last_os_error())
            }
        }
        #[cfg(not(target_os = "linux"))]
        Err(io::Error::from(io::ErrorKind::Unsupported))
    }

    /// Touch the whole buffer, ask for 2 MiB backing (`collapse_huge`), and
    /// check that it took: for buffers the GPU reads in place. Below 90% huge
    /// pages it warns (the first few times per process, naming `label`) and
    /// counts the buffer in [`small_page_buffers`]. Never fails: 4 KiB backing
    /// is slower, and was once suspected of wrong GPU reads under page
    /// migration (2026-10-20), so it must be visible, not fatal.
    pub fn touch_huge(&mut self, label: &str) {
        self.set_label(label);
        if std::env::var_os("NOMICS_TIMING").is_some() {
            eprintln!(
                "[mem] hostbuf {label}: {} MiB (live {:.2} GB)",
                self.len >> 20,
                hostbuf_bytes() as f64 / 1e9
            );
        }
        self.as_mut_slice().fill(0);
        let _ = self.collapse_huge();
        // Buffers under one huge page cannot be huge-page backed.
        if self.len < HUGE_PAGE {
            return;
        }
        let Some(huge) = self.huge_page_bytes() else {
            return;
        };
        // With guards the last, partial huge page is split by the tail guard.
        let expected = if self.is_guarded() {
            self.len - self.len % HUGE_PAGE
        } else {
            self.len
        };
        if huge.saturating_mul(10) < expected.saturating_mul(9) {
            let n = SMALL_PAGE_BUFFERS.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            if n < 5 {
                eprintln!(
                    "warning: {label}: only {} of {} MiB on 2 MiB pages (memory fragmented; evict the page cache before timing runs)",
                    huge >> 20,
                    self.len >> 20
                );
            }
        }
    }

    /// Bytes of this buffer currently backed by transparent huge pages, read
    /// from `/proc/self/smaps`. `None` where smaps is unavailable (macOS).
    pub fn huge_page_bytes(&self) -> Option<usize> {
        let smaps = std::fs::read_to_string("/proc/self/smaps").ok()?;
        let (start, end) = (
            self.ptr.as_ptr() as usize,
            self.ptr.as_ptr() as usize + self.len,
        );
        let mut total = 0usize;
        let mut in_range = false;
        for line in smaps.lines() {
            if let Some((range, _)) = line.split_once(' ')
                && let Some((lo, hi)) = range.split_once('-')
                && let (Ok(lo), Ok(hi)) =
                    (usize::from_str_radix(lo, 16), usize::from_str_radix(hi, 16))
            {
                in_range = lo < end && hi > start;
                continue;
            }
            if in_range && let Some(rest) = line.strip_prefix("AnonHugePages:") {
                let kb: usize = rest.trim().trim_end_matches("kB").trim().parse().ok()?;
                total += kb * 1024;
            }
        }
        Some(total)
    }
}

impl Drop for HostBuf {
    fn drop(&mut self) {
        if self.is_guarded() {
            diag::unregister(self.ptr.as_ptr() as usize);
        }
        unsafe {
            libc::munmap(self.map_base, self.map_len);
        }
        HOSTBUF_BYTES.fetch_sub(self.map_len, std::sync::atomic::Ordering::Relaxed);
    }
}

/// CPU-only memory of a fixed length on ordinary pages: an anonymous mapping
/// that is returned to the system the moment it is dropped. For large,
/// short-lived blocks freed on another thread than the one that made them
/// (a chunk's alignment records), which an allocator may keep resident
/// until that thread allocates again. Zero-filled lazily; never for the GPU.
pub struct PageBuf {
    ptr: NonNull<u8>,
    len: usize,
}

// Plain memory with a single owner; sharing follows &/&mut rules.
unsafe impl Send for PageBuf {}
unsafe impl Sync for PageBuf {}

impl PageBuf {
    pub fn new(len: usize) -> io::Result<Self> {
        if len == 0 {
            return Ok(Self {
                ptr: NonNull::dangling(),
                len: 0,
            });
        }
        let base = unsafe {
            libc::mmap(
                std::ptr::null_mut(),
                len,
                libc::PROT_READ | libc::PROT_WRITE,
                libc::MAP_PRIVATE | libc::MAP_ANONYMOUS,
                -1,
                0,
            )
        };
        if base == libc::MAP_FAILED {
            return Err(io::Error::last_os_error());
        }
        Ok(Self {
            ptr: NonNull::new(base.cast()).expect("mmap returned null"),
            len,
        })
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    pub fn as_slice(&self) -> &[u8] {
        unsafe { std::slice::from_raw_parts(self.ptr.as_ptr(), self.len) }
    }

    pub fn as_mut_slice(&mut self) -> &mut [u8] {
        unsafe { std::slice::from_raw_parts_mut(self.ptr.as_ptr(), self.len) }
    }

    /// Return the whole pages inside `offset..offset + len` to the system
    /// now (`MADV_DONTNEED`); they read as zeros afterwards. Partial pages at
    /// either end are kept, so neighbouring data is never touched.
    pub fn release(&mut self, offset: usize, len: usize) {
        let page = page_size();
        let end = offset.saturating_add(len).min(self.len);
        let Some(start) = offset.checked_next_multiple_of(page) else {
            return;
        };
        let end = end - end % page;
        if start < end {
            unsafe {
                libc::madvise(
                    self.ptr.as_ptr().add(start).cast(),
                    end - start,
                    libc::MADV_DONTNEED,
                );
            }
        }
    }
}

/// The base page size, for laying out data that [`PageBuf::release`] frees.
pub fn page_size() -> usize {
    static PAGE: std::sync::OnceLock<usize> = std::sync::OnceLock::new();
    *PAGE.get_or_init(|| {
        let n = unsafe { libc::sysconf(libc::_SC_PAGESIZE) };
        if n > 0 { n as usize } else { 4096 }
    })
}

impl Drop for PageBuf {
    fn drop(&mut self) {
        if self.len > 0 {
            unsafe {
                libc::munmap(self.ptr.as_ptr().cast(), self.len);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    /// An offset whose page rounding overflows releases nothing (review
    /// finding 2026-10-22: it used to wrap to 0 and discard the first page).
    #[test]
    fn page_buf_release_overflowing_offset_keeps_data() {
        let page = page_size();
        let mut buf = PageBuf::new(page).unwrap();
        buf.as_mut_slice().fill(7);
        buf.release(usize::MAX, 1);
        buf.release(usize::MAX - 3, 10);
        assert!(buf.as_mut_slice().iter().all(|&b| b == 7));
    }

    use super::*;

    #[test]
    fn page_buf_is_zeroed_and_writable() {
        let mut buf = PageBuf::new(3 * 4096 + 5).unwrap();
        assert_eq!(buf.len(), 3 * 4096 + 5);
        assert!(buf.as_slice().iter().all(|&b| b == 0));
        buf.as_mut_slice().fill(7);
        assert!(buf.as_slice().iter().all(|&b| b == 7));
        assert!(PageBuf::new(0).unwrap().is_empty());
    }

    #[test]
    fn page_buf_release_frees_only_whole_pages_inside() {
        let page = page_size();
        let mut buf = PageBuf::new(4 * page).unwrap();
        buf.as_mut_slice().fill(9);
        // Pages 1 and 2 lie inside; the partial pages 0 and 3 stay.
        buf.release(page - 1, 2 * page + 2);
        let s = buf.as_slice();
        assert!(s[..page].iter().all(|&b| b == 9));
        assert!(s[3 * page..].iter().all(|&b| b == 9));
        #[cfg(target_os = "linux")]
        assert!(s[page..3 * page].iter().all(|&b| b == 0));
        buf.release(0, 0);
        buf.release(3 * page, 10 * page);
    }

    #[test]
    fn aligned_and_writable() {
        let mut buf = HostBuf::new(3 * HUGE_PAGE + 17).unwrap();
        assert_eq!(buf.as_ptr() as usize % HUGE_PAGE, 0);
        assert_eq!(buf.len(), 3 * HUGE_PAGE + 17);
        buf.as_mut_slice().fill(0xA5);
        assert!(buf.as_slice().iter().all(|&b| b == 0xA5));
    }

    #[test]
    fn zero_length_maps_one_page() {
        let buf = HostBuf::new(0).unwrap();
        assert!(buf.is_empty());
        assert_eq!(buf.as_ptr() as usize % HUGE_PAGE, 0);
    }

    #[test]
    fn guarded_layout_keeps_alignment_and_intact_canaries() {
        let page = page_size();
        for len in [0, 1, 17, page, page + 3, 3 * HUGE_PAGE + 17, 2 * HUGE_PAGE] {
            let mut buf = HostBuf::with_guards(len).unwrap();
            assert!(buf.is_guarded());
            assert_eq!(buf.as_ptr() as usize % HUGE_PAGE, 0);
            assert_eq!(buf.len(), len);
            assert_eq!(buf.accessible_len(), len.next_multiple_of(page));
            buf.as_mut_slice().fill(0xA5);
            assert!(buf.as_slice().iter().all(|&b| b == 0xA5));
            buf.touch_huge("guarded layout test");
            assert!(buf.as_slice().iter().all(|&b| b == 0));
            assert_eq!(buf.check_canaries("guarded layout test"), Ok(()));
        }
        // The plain layout is unchanged.
        let plain = HostBuf::new(17).unwrap();
        if !diag() {
            assert!(!plain.is_guarded());
            assert_eq!(plain.accessible_len(), HUGE_PAGE);
            assert_eq!(plain.check_canaries("plain"), Ok(()));
        }
    }

    #[test]
    fn canary_damage_is_detected_and_reported() {
        let len = 3 * page_size() + 5;
        let mut buf = HostBuf::with_guards(len).unwrap();
        buf.set_label("canary damage test");
        let data = buf.as_mut_ptr();
        // Writing the slack past `len` does not fault: it damages the canary.
        unsafe { data.add(len + 2).write_volatile(0x00) };
        let err = buf.check_canaries("damage test").unwrap_err();
        assert_eq!(err.name, "damage test");
        assert_eq!(err.offset, (len + 2) as isize);
        assert_eq!(err.damaged, 1);
        assert_eq!(err.bytes[0].0, 0x00);
        let text = err.to_string();
        assert!(
            text.contains("tail canary damaged") && text.contains("offset +"),
            "{text}"
        );
        // The registry sees it under the buffer's label.
        assert!(
            diag::verify_all_canaries()
                .iter()
                .any(|d| d.name == "canary damage test" && d.offset == (len + 2) as isize)
        );
        assert!(diag::describe_all().contains("canary damage test"));
        // Repair, then damage the head canary instead.
        unsafe {
            let at = data.add(len + 2) as usize;
            data.add(len + 2).write_volatile(diag::canary_byte(at));
        }
        assert_eq!(buf.check_canaries("damage test"), Ok(()));
        unsafe { data.sub(8).write_volatile(0x11) };
        let err = buf.check_canaries("damage test").unwrap_err();
        assert_eq!(err.offset, -8);
        assert!(err.to_string().contains("head canary damaged"));
        drop(buf);
        // Dropped buffers leave the registry.
        assert!(!diag::describe_all().contains("canary damage test"));
    }

    /// Run in a child process by [`write_past_end_hits_guard_page`].
    #[test]
    fn guard_page_child() {
        let Ok(side) = std::env::var("NOMICS_MEM_GUARD_CHILD") else {
            return;
        };
        let len = page_size() + 100;
        let mut buf = HostBuf::with_guards(len).unwrap();
        buf.set_label("guard child");
        let data = buf.as_mut_ptr();
        unsafe {
            // The slack is writable.
            data.add(len).write_volatile(1);
            match side.as_str() {
                "tail" => data
                    .add(len.next_multiple_of(page_size()))
                    .write_volatile(1),
                _ => data.sub(page_size() + 1).write_volatile(1),
            }
        }
        unreachable!("the guard page did not fault");
    }

    #[test]
    fn write_past_end_hits_guard_page() {
        use std::os::unix::process::ExitStatusExt;
        for side in ["tail", "head"] {
            let out = std::process::Command::new(std::env::current_exe().unwrap())
                .args([
                    "tests::guard_page_child",
                    "--exact",
                    "--nocapture",
                    "--test-threads=1",
                ])
                .env("NOMICS_MEM_GUARD_CHILD", side)
                .output()
                .unwrap();
            let signal = out.status.signal();
            let stderr = String::from_utf8_lossy(&out.stderr);
            assert!(
                signal == Some(libc::SIGSEGV) || signal == Some(libc::SIGBUS),
                "{side}: child exited with {:?}, stderr: {stderr}",
                out.status
            );
            let part = if side == "tail" {
                "tail guard page"
            } else {
                "head guard page"
            };
            assert!(
                stderr.contains(&format!("in the {part} of HostBuf 'guard child'")),
                "{side}: {stderr}"
            );
        }
    }

    static CHAIN_HITS: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
    static CHAIN_PAGE: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);

    /// The "previous" handler of [`fault_handler_chains_child`]: counts, and
    /// makes its own PROT_NONE page writable so the fault is recoverable.
    extern "C" fn chain_previous(_: libc::c_int, _: *mut libc::siginfo_t, _: *mut libc::c_void) {
        CHAIN_HITS.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
        let page = CHAIN_PAGE.load(std::sync::atomic::Ordering::SeqCst);
        if page != 0 {
            unsafe {
                libc::mprotect(
                    page as *mut libc::c_void,
                    page_size(),
                    libc::PROT_READ | libc::PROT_WRITE,
                )
            };
        }
    }

    /// Run in a child process by [`fault_handler_chains_foreign_signals`].
    #[test]
    fn fault_handler_chains_child() {
        if std::env::var_os("NOMICS_MEM_CHAIN_CHILD").is_none() {
            return;
        }
        use std::sync::atomic::Ordering::SeqCst;
        unsafe {
            for sig in [libc::SIGSEGV, libc::SIGBUS] {
                let mut action: libc::sigaction = std::mem::zeroed();
                action.sa_sigaction = chain_previous as *const () as usize;
                action.sa_flags = libc::SA_SIGINFO;
                libc::sigemptyset(&mut action.sa_mask);
                assert_eq!(libc::sigaction(sig, &action, std::ptr::null_mut()), 0);
            }
        }
        // The first guarded buffer installs ours over chain_previous.
        let _buf = HostBuf::with_guards(100).unwrap();
        assert!(diag::fault_handler_installed());
        // A software-sent SIGSEGV goes to the previous handler; ours stays.
        unsafe { libc::raise(libc::SIGSEGV) };
        assert_eq!(CHAIN_HITS.load(SeqCst), 1);
        assert!(diag::fault_handler_installed());
        // A real fault outside every registered mapping: chained, recovered
        // by the previous handler, and ours stays installed.
        let page = unsafe {
            libc::mmap(
                std::ptr::null_mut(),
                page_size(),
                libc::PROT_NONE,
                libc::MAP_PRIVATE | libc::MAP_ANONYMOUS,
                -1,
                0,
            )
        };
        assert_ne!(page, libc::MAP_FAILED);
        CHAIN_PAGE.store(page as usize, SeqCst);
        unsafe { (page as *mut u8).write_volatile(7) };
        assert_eq!(unsafe { (page as *const u8).read_volatile() }, 7);
        assert_eq!(CHAIN_HITS.load(SeqCst), 2);
        assert!(diag::fault_handler_installed());
        println!("chain ok");
    }

    #[test]
    fn fault_handler_chains_foreign_signals() {
        let out = std::process::Command::new(std::env::current_exe().unwrap())
            .args([
                "tests::fault_handler_chains_child",
                "--exact",
                "--nocapture",
                "--test-threads=1",
            ])
            .env("NOMICS_MEM_CHAIN_CHILD", "1")
            .output()
            .unwrap();
        let stdout = String::from_utf8_lossy(&out.stdout);
        assert!(
            out.status.success() && stdout.contains("chain ok"),
            "child: {:?}\n{stdout}\n{}",
            out.status,
            String::from_utf8_lossy(&out.stderr)
        );
    }

    #[test]
    fn diag_modes_parse() {
        use diag::{Mode, parse_mode};
        assert_eq!(parse_mode(None), Mode::Off);
        assert_eq!(parse_mode(Some("")), Mode::Off);
        assert_eq!(parse_mode(Some("0")), Mode::Off);
        assert_eq!(parse_mode(Some("regions")), Mode::Regions);
        assert_eq!(parse_mode(Some("REGIONS")), Mode::Regions);
        assert_eq!(parse_mode(Some("1")), Mode::Full);
        assert_eq!(parse_mode(Some("full")), Mode::Full);
    }

    #[cfg(target_os = "linux")]
    #[test]
    fn touched_buffer_is_huge_page_backed() {
        let mut buf = HostBuf::new(64 * HUGE_PAGE).unwrap();
        buf.as_mut_slice().fill(1);
        let huge = buf.huge_page_bytes().expect("smaps readable on linux");
        // THP can be partially denied under memory pressure; require a majority.
        assert!(
            huge >= buf.len() / 2,
            "only {huge} of {} bytes on huge pages",
            buf.len()
        );
    }
}

/// A whole file mapped read-only, for large inputs read once front to back
/// (compressed FASTQs). Its pages are reclaimable page cache rather than
/// anonymous memory, so they do not add to the process's committed memory.
/// The file must not be truncated or rewritten while mapped.
pub struct MappedFile {
    ptr: *const u8,
    len: usize,
}

// Read-only memory with a single owner.
unsafe impl Send for MappedFile {}
unsafe impl Sync for MappedFile {}

impl MappedFile {
    pub fn open(path: impl AsRef<std::path::Path>) -> io::Result<Self> {
        use std::os::fd::AsRawFd;
        let file = std::fs::File::open(path)?;
        let len = usize::try_from(file.metadata()?.len())
            .map_err(|_| io::Error::other("file larger than the address space"))?;
        if len == 0 {
            return Ok(Self {
                ptr: NonNull::<u8>::dangling().as_ptr(),
                len: 0,
            });
        }
        let base = unsafe {
            libc::mmap(
                std::ptr::null_mut(),
                len,
                libc::PROT_READ,
                libc::MAP_PRIVATE,
                file.as_raw_fd(),
                0,
            )
        };
        if base == libc::MAP_FAILED {
            return Err(io::Error::last_os_error());
        }
        // Read-ahead hint only; failure is harmless.
        unsafe {
            libc::madvise(base, len, libc::MADV_SEQUENTIAL);
        }
        Ok(Self {
            ptr: base as *const u8,
            len,
        })
    }
}

impl std::ops::Deref for MappedFile {
    type Target = [u8];
    fn deref(&self) -> &[u8] {
        unsafe { std::slice::from_raw_parts(self.ptr, self.len) }
    }
}

impl Drop for MappedFile {
    fn drop(&mut self) {
        if self.len > 0 {
            unsafe {
                libc::munmap(self.ptr as *mut libc::c_void, self.len);
            }
        }
    }
}
