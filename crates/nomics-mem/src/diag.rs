//! `NOMICS_DIAG=1`: memory-corruption diagnostics for host buffers the GPU
//! reads and writes in place (docs/DIAG.md).
//!
//! With diagnostics on, every [`HostBuf`](crate::HostBuf) is mapped as
//!
//! ```text
//! [head guard ...][head canary page][data: 2 MiB aligned, len][slack canary to page end][tail guard ...]
//! ```
//!
//! The guards are `PROT_NONE`: a CPU access faults at once (reported here by
//! a SIGSEGV/SIGBUS handler that names the buffer), and on GB10 a GPU access
//! should fault too, because ATS walks the CPU page tables, and surface as a
//! CUDA error (`nomics-gpu` makes GPU errors fatal in this mode). The
//! canaries catch writes that stop short of a guard: they are checked after
//! every GPU ticket wait ([`check_all_canaries`]). Every guarded buffer is
//! listed in a registry, which the fault handler and the region checks in
//! `nomics-align` use to name the buffers around an address.
//!
//! With diagnostics off (the default) nothing here runs: buffers keep the
//! original layout and no registry exists.
//!
//! `NOMICS_DIAG=regions` keeps the original HostBuf layout (no guards,
//! canaries or registry, so buffer adjacency is as in a normal run) and
//! enables only the checks that do not move memory: the region invariants
//! (`nomics-align`), the checked `mem_pair` report and fatal GPU errors.
use std::sync::Mutex;
use std::sync::atomic::{AtomicPtr, AtomicUsize, Ordering};

/// What `NOMICS_DIAG` turns on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Mode {
    /// Unset, empty or `0`.
    Off,
    /// `regions`: region invariants and fatal GPU errors, original layout.
    Regions,
    /// Anything else (e.g. `1`): also guard pages, canaries and the registry.
    Full,
}

/// Parse a `NOMICS_DIAG` value.
pub fn parse_mode(value: Option<&str>) -> Mode {
    match value.map(str::trim) {
        None | Some("") | Some("0") => Mode::Off,
        Some(v) if v.eq_ignore_ascii_case("regions") => Mode::Regions,
        Some(_) => Mode::Full,
    }
}

/// The `NOMICS_DIAG` mode, read once; the first read prints which one is on.
pub fn mode() -> Mode {
    static MODE: std::sync::OnceLock<Mode> = std::sync::OnceLock::new();
    *MODE.get_or_init(|| {
        let mode = parse_mode(std::env::var("NOMICS_DIAG").ok().as_deref());
        match mode {
            Mode::Off => {}
            Mode::Regions => eprintln!(
                "NOMICS_DIAG=regions: region invariants, checked mem_pair and fatal GPU errors; original HostBuf layout (no guard pages, canaries or registry)"
            ),
            Mode::Full => eprintln!(
                "NOMICS_DIAG (full): HostBufs get guard pages and canaries (page {} bytes); region invariant failures and GPU errors abort",
                crate::page_size()
            ),
        }
        mode
    })
}

/// Whether any diagnostics are on (`NOMICS_DIAG` set, either mode): the
/// region checks and fatal GPU errors.
pub fn diag() -> bool {
    mode() != Mode::Off
}

/// Whether HostBufs get guard pages, canaries and the registry (full mode).
pub fn guards() -> bool {
    mode() == Mode::Full
}

/// The canary byte at address `addr`: an 8-byte pattern by address, so a
/// copy displaced by a multiple of 8 does not pass as intact.
#[inline]
pub(crate) fn canary_byte(addr: usize) -> u8 {
    const PATTERN: [u8; 8] = [0xDE, 0xAD, 0xC0, 0xDE, 0x5A, 0xFE, 0xB0, 0xF5];
    PATTERN[addr & 7]
}

/// One registered buffer. All fields are atomics so the fault handler can
/// read them without a lock (it may run while another thread holds one).
struct Slot {
    map_base: AtomicUsize,
    map_len: AtomicUsize,
    data: AtomicUsize,
    len: AtomicUsize,
    /// End of the accessible pages (data + len rounded up to a page).
    end: AtomicUsize,
    /// An interned, leaked, immutable `&'static str` (null: none), so the
    /// fault handler reads a label with one atomic load.
    label: AtomicPtr<&'static str>,
}

const SLOTS: usize = 8192;

#[allow(clippy::declare_interior_mutable_const)]
const EMPTY: Slot = Slot {
    map_base: AtomicUsize::new(0),
    map_len: AtomicUsize::new(0),
    data: AtomicUsize::new(0),
    len: AtomicUsize::new(0),
    end: AtomicUsize::new(0),
    label: AtomicPtr::new(std::ptr::null_mut()),
};

static REGISTRY: [Slot; SLOTS] = [EMPTY; SLOTS];

/// Serializes registration, removal and canary checks, so a buffer is never
/// unmapped while a check reads it. The fault handler takes no lock.
static LOCK: Mutex<()> = Mutex::new(());

fn lock() -> std::sync::MutexGuard<'static, ()> {
    LOCK.lock().unwrap_or_else(|e| e.into_inner())
}

/// Labels are few distinct strings, interned once and never freed: each is
/// a leaked `&'static str` behind a leaked, immutable cell, published to a
/// slot as one pointer.
fn intern(label: &str) -> *mut &'static str {
    struct Cell(*mut &'static str);
    // SAFETY: the cells are leaked and never written after creation.
    unsafe impl Send for Cell {}
    static LABELS: Mutex<Vec<Cell>> = Mutex::new(Vec::new());
    let mut labels = LABELS.lock().unwrap_or_else(|e| e.into_inner());
    // SAFETY: every cell points to a leaked, immutable &'static str.
    if let Some(c) = labels.iter().find(|c| unsafe { *c.0 } == label) {
        return c.0;
    }
    let s: &'static str = Box::leak(label.to_owned().into_boxed_str());
    let cell = Box::leak(Box::new(s)) as *mut &'static str;
    labels.push(Cell(cell));
    cell
}

fn slot_label(slot: &Slot) -> &'static str {
    let cell = slot.label.load(Ordering::Acquire);
    if cell.is_null() {
        return "unlabelled HostBuf";
    }
    // SAFETY: a non-null label is an interned cell: leaked, immutable, and
    // holding a leaked &'static str; one load gives a consistent (ptr, len).
    unsafe { *cell }
}

pub(crate) struct Layout {
    pub map_base: usize,
    pub map_len: usize,
    pub data: usize,
    pub len: usize,
    pub end: usize,
}

pub(crate) fn register(l: &Layout) {
    install_fault_handler();
    let _g = lock();
    for slot in &REGISTRY {
        if slot.map_base.load(Ordering::Relaxed) == 0 {
            slot.map_len.store(l.map_len, Ordering::Relaxed);
            slot.data.store(l.data, Ordering::Relaxed);
            slot.len.store(l.len, Ordering::Relaxed);
            slot.end.store(l.end, Ordering::Relaxed);
            slot.label.store(std::ptr::null_mut(), Ordering::Relaxed);
            slot.map_base.store(l.map_base, Ordering::Release);
            return;
        }
    }
    static WARNED: std::sync::Once = std::sync::Once::new();
    WARNED.call_once(|| {
        eprintln!("NOMICS_DIAG: more than {SLOTS} live HostBufs; later ones are guarded but not canary-checked")
    });
}

fn find(data: usize) -> Option<&'static Slot> {
    REGISTRY
        .iter()
        .find(|s| s.map_base.load(Ordering::Acquire) != 0 && s.data.load(Ordering::Relaxed) == data)
}

pub(crate) fn unregister(data: usize) {
    let _g = lock();
    if let Some(slot) = find(data) {
        slot.map_base.store(0, Ordering::Release);
    }
}

pub(crate) fn set_label(data: usize, label: &str) {
    let label = intern(label);
    let _g = lock();
    if let Some(slot) = find(data) {
        slot.label.store(label, Ordering::Release);
    }
}

/// Damage found in a buffer's canaries.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CanaryDamage {
    pub name: String,
    pub data: usize,
    pub len: usize,
    /// Offset of the first damaged byte from the buffer's data start
    /// (negative: the head canary).
    pub offset: isize,
    /// Damaged bytes in the canary area.
    pub damaged: usize,
    /// Up to 32 bytes from the first damaged one: (found, expected).
    pub bytes: Vec<(u8, u8)>,
}

impl std::fmt::Display for CanaryDamage {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let side = if self.offset < 0 { "head" } else { "tail" };
        write!(
            f,
            "HostBuf '{}' (data {:#x}, len {}) {side} canary damaged: {} bytes, first at offset {:+} ({:#x}); found/expected:",
            self.name,
            self.data,
            self.len,
            self.damaged,
            self.offset,
            self.data as isize + self.offset
        )?;
        for (found, expected) in &self.bytes {
            write!(f, " {found:02x}/{expected:02x}")?;
        }
        Ok(())
    }
}

/// Fill `[from, to)` with the canary pattern.
///
/// # Safety
/// The range must be writable memory owned by the caller.
pub(crate) unsafe fn fill(from: usize, to: usize) {
    for addr in from..to {
        unsafe { std::ptr::write_volatile(addr as *mut u8, canary_byte(addr)) };
    }
}

/// Check `[from, to)` against the canary pattern.
///
/// # Safety
/// The range must be readable, mapped for the duration of the call.
/// (first damaged address, damaged bytes, up to 32 (found, expected) pairs).
type Damage = (usize, usize, Vec<(u8, u8)>);

unsafe fn scan(from: usize, to: usize) -> Option<Damage> {
    // Fast path: aligned words (the pattern repeats every 8 bytes by
    // address, so every aligned word must equal WORD). Bytes on a mismatch.
    let word = u64::from_ne_bytes(std::array::from_fn(canary_byte));
    let head = from.next_multiple_of(8).min(to);
    let tail = head.max(to - to % 8);
    let mut intact = (from..head)
        .all(|a| unsafe { std::ptr::read_volatile(a as *const u8) } == canary_byte(a))
        && (tail..to).all(|a| unsafe { std::ptr::read_volatile(a as *const u8) } == canary_byte(a));
    let mut a = head;
    while intact && a < tail {
        intact = unsafe { std::ptr::read_volatile(a as *const u64) } == word;
        a += 8;
    }
    if intact {
        return None;
    }
    let mut first = None;
    let mut damaged = 0;
    let mut bytes = Vec::new();
    let mut addr = from;
    while addr < to {
        let found = unsafe { std::ptr::read_volatile(addr as *const u8) };
        let expected = canary_byte(addr);
        if found != expected {
            damaged += 1;
            first.get_or_insert(addr);
        }
        if first.is_some() && bytes.len() < 32 {
            bytes.push((found, expected));
        }
        addr += 1;
    }
    first.map(|f| (f, damaged, bytes))
}

pub(crate) fn check_layout(
    name: &str,
    data: usize,
    len: usize,
    end: usize,
) -> Result<(), CanaryDamage> {
    let page = crate::page_size();
    let damage = |at: usize, damaged, bytes| CanaryDamage {
        name: name.to_owned(),
        data,
        len,
        offset: at as isize - data as isize,
        damaged,
        bytes,
    };
    // SAFETY: the head canary page and the slack are this buffer's own
    // read-write pages, mapped while the caller holds the buffer (or the
    // registry lock, which unmapping takes).
    if let Some((at, n, bytes)) = unsafe { scan(data - page, data) } {
        return Err(damage(at, n, bytes));
    }
    if let Some((at, n, bytes)) = unsafe { scan(data + len, end) } {
        return Err(damage(at, n, bytes));
    }
    Ok(())
}

/// Check every registered buffer's canaries; returns the damage found.
pub fn verify_all_canaries() -> Vec<CanaryDamage> {
    let _g = lock();
    let mut out = Vec::new();
    for slot in &REGISTRY {
        if slot.map_base.load(Ordering::Acquire) == 0 {
            continue;
        }
        let (data, len, end) = (
            slot.data.load(Ordering::Relaxed),
            slot.len.load(Ordering::Relaxed),
            slot.end.load(Ordering::Relaxed),
        );
        if let Err(d) = check_layout(slot_label(slot), data, len, end) {
            out.push(d);
        }
    }
    out
}

/// With `NOMICS_DIAG`: check every live HostBuf's canaries and abort the
/// process, naming each damaged buffer, if any is damaged. `context` says
/// where the check ran (e.g. which GPU wait). A no-op otherwise.
pub fn check_all_canaries(context: &str) {
    if !guards() {
        return;
    }
    SWEEPS.fetch_add(1, Ordering::Relaxed);
    let damage = verify_all_canaries();
    if damage.is_empty() {
        return;
    }
    let mut report = format!("NOMICS_DIAG: canary damage detected after {context}:\n");
    for d in &damage {
        report.push_str(&format!("  {d}\n"));
    }
    report.push_str(&describe_all());
    fatal(&report);
}

static SWEEPS: AtomicUsize = AtomicUsize::new(0);

/// One line on what diagnostics checked so far (for the end of a run), or
/// `None` without `NOMICS_DIAG`.
pub fn summary() -> Option<String> {
    diag().then(|| {
        let live = REGISTRY
            .iter()
            .filter(|s| s.map_base.load(Ordering::Acquire) != 0)
            .count();
        format!(
            "NOMICS_DIAG: {} canary sweeps after GPU waits found no damage; {live} guarded HostBufs live",
            SWEEPS.load(Ordering::Relaxed)
        )
    })
    .map(|line| {
        if guards() {
            line
        } else {
            "NOMICS_DIAG=regions: no region invariant or GPU failure (no guard pages or canaries in this mode)".to_owned()
        }
    })
}

/// Print `report` to stderr in one write and abort the process (a core
/// dump keeps the state; a panic would only end one thread).
pub fn fatal(report: &str) -> ! {
    use std::io::Write as _;
    let thread = std::thread::current();
    let text = format!(
        "\n==== NOMICS_DIAG FATAL (thread {}) ====\n{report}\n==== aborting ====\n",
        thread.name().unwrap_or("?")
    );
    let _ = std::io::stderr().lock().write_all(text.as_bytes());
    std::process::abort();
}

/// One line per live registered HostBuf, ordered by address.
pub fn describe_all() -> String {
    describe_near(0, usize::MAX)
}

/// The `n` registered HostBufs whose mappings lie nearest `addr` (all of
/// them for `n == usize::MAX`), ordered by address, one line each.
pub fn describe_near(addr: usize, n: usize) -> String {
    let _g = lock();
    let mut rows: Vec<(usize, usize, String)> = Vec::new();
    for slot in &REGISTRY {
        let base = slot.map_base.load(Ordering::Acquire);
        if base == 0 {
            continue;
        }
        let map_len = slot.map_len.load(Ordering::Relaxed);
        let (data, len, end) = (
            slot.data.load(Ordering::Relaxed),
            slot.len.load(Ordering::Relaxed),
            slot.end.load(Ordering::Relaxed),
        );
        let distance = if addr < base {
            base - addr
        } else {
            addr.saturating_sub(base + map_len)
        };
        rows.push((
            distance,
            data,
            format!(
                "  HostBuf '{}': data {:#x}..{:#x} (len {}), accessible to {:#x}, mapping {:#x}..{:#x}",
                slot_label(slot),
                data,
                data + len,
                len,
                end,
                base,
                base + map_len
            ),
        ));
    }
    let total = rows.len();
    if n < rows.len() {
        rows.sort_by_key(|r| r.0);
        rows.truncate(n);
    }
    rows.sort_by_key(|r| r.1);
    let mut out = if n == usize::MAX {
        format!("live HostBufs ({total}):\n")
    } else {
        format!("HostBufs nearest {addr:#x} ({} of {total}):\n", rows.len())
    };
    for (_, _, line) in rows {
        out.push_str(&line);
        out.push('\n');
    }
    out
}

// ---- Fault handler ----

static OLD_SEGV: AtomicUsize = AtomicUsize::new(0);
static OLD_BUS: AtomicUsize = AtomicUsize::new(0);

/// Install a SIGSEGV/SIGBUS handler (once) that, for a fault inside a
/// registered buffer's mapping, prints which buffer and which part of it
/// (guard, canary, data) before the previous handler (Rust's, then the
/// default action: a core dump) runs.
fn install_fault_handler() {
    static ONCE: std::sync::Once = std::sync::Once::new();
    ONCE.call_once(|| unsafe {
        for (sig, old) in [(libc::SIGSEGV, &OLD_SEGV), (libc::SIGBUS, &OLD_BUS)] {
            let mut action: libc::sigaction = std::mem::zeroed();
            action.sa_sigaction = on_fault as *const () as usize;
            action.sa_flags = libc::SA_SIGINFO | libc::SA_ONSTACK;
            libc::sigemptyset(&mut action.sa_mask);
            let saved = Box::leak(Box::new(std::mem::zeroed::<libc::sigaction>()));
            if libc::sigaction(sig, &action, saved) == 0 {
                old.store(saved as *mut libc::sigaction as usize, Ordering::Release);
            }
        }
    });
}

/// A fixed buffer formatted without allocating (signal-safe).
struct Msg {
    buf: [u8; 512],
    at: usize,
}

impl Msg {
    fn str(&mut self, s: &str) -> &mut Self {
        for &b in s.as_bytes() {
            if self.at < self.buf.len() {
                self.buf[self.at] = b;
                self.at += 1;
            }
        }
        self
    }
    fn hex(&mut self, v: usize) -> &mut Self {
        self.str("0x");
        let mut started = false;
        for shift in (0..16).rev() {
            let d = (v >> (shift * 4)) & 15;
            if d != 0 || started || shift == 0 {
                started = true;
                let c = b"0123456789abcdef"[d];
                self.str(std::str::from_utf8(&[c]).unwrap_or("?"));
            }
        }
        self
    }
    fn dec(&mut self, mut v: usize) -> &mut Self {
        let mut digits = [0u8; 20];
        let mut n = 0;
        loop {
            digits[n] = b'0' + (v % 10) as u8;
            n += 1;
            v /= 10;
            if v == 0 {
                break;
            }
        }
        for i in (0..n).rev() {
            self.str(std::str::from_utf8(&digits[i..=i]).unwrap_or("?"));
        }
        self
    }
    fn emit(&self) {
        unsafe { libc::write(2, self.buf.as_ptr().cast(), self.at) };
    }
}

#[cfg(any(target_os = "linux", target_os = "android"))]
unsafe fn fault_addr(info: *const libc::siginfo_t) -> usize {
    unsafe { (*info).si_addr() as usize }
}

#[cfg(not(any(target_os = "linux", target_os = "android")))]
unsafe fn fault_addr(info: *const libc::siginfo_t) -> usize {
    unsafe { (*info).si_addr as usize }
}

/// Whether a signal came from the kernel for a faulting access (not sent by
/// `kill`/`raise`/`sigqueue`): positive kernel codes, below the user-sent
/// range macOS uses (`SI_USER` = 0x10001); Linux user codes are <= 0.
fn kernel_fault(info: *const libc::siginfo_t) -> bool {
    let code = unsafe { (*info).si_code };
    code > 0 && code < 0x10000
}

/// The registered part of a buffer's mapping containing `addr`: guards and
/// canaries only (`None` for data, or outside every registered mapping).
fn guarded_part(addr: usize) -> Option<(&'static Slot, &'static str)> {
    let page = crate::page_size();
    for slot in &REGISTRY {
        let base = slot.map_base.load(Ordering::Acquire);
        let map_len = slot.map_len.load(Ordering::Relaxed);
        if base == 0 || addr < base || addr >= base + map_len {
            continue;
        }
        let (data, len, end) = (
            slot.data.load(Ordering::Relaxed),
            slot.len.load(Ordering::Relaxed),
            slot.end.load(Ordering::Relaxed),
        );
        let part = if addr < data.wrapping_sub(page) {
            "head guard page"
        } else if addr < data {
            "head canary page"
        } else if addr < data + len {
            return None;
        } else if addr < end {
            "tail slack canary"
        } else {
            "tail guard page"
        };
        return Some((slot, part));
    }
    None
}

extern "C" fn on_fault(sig: libc::c_int, info: *mut libc::siginfo_t, ctx: *mut libc::c_void) {
    let addr = unsafe { fault_addr(info) };
    if kernel_fault(info)
        && let Some((slot, part)) = guarded_part(addr)
    {
        let data = slot.data.load(Ordering::Relaxed);
        let len = slot.len.load(Ordering::Relaxed);
        let mut m = Msg {
            buf: [0; 512],
            at: 0,
        };
        m.str("\n==== NOMICS_DIAG: ")
            .str(if sig == libc::SIGBUS {
                "SIGBUS"
            } else {
                "SIGSEGV"
            })
            .str(" at ")
            .hex(addr)
            .str(" in the ")
            .str(part)
            .str(" of HostBuf '")
            .str(slot_label(slot))
            .str("' (data ")
            .hex(data)
            .str(", len ")
            .dec(len)
            .str(", offset ");
        if addr >= data {
            m.str("+").dec(addr - data);
        } else {
            m.str("-").dec(data - addr);
        }
        m.str(") ====\n");
        m.emit();
        // Ours: the default action (core dump) when the access runs again.
        unsafe { libc::signal(sig, libc::SIG_DFL) };
        return;
    }
    // Not ours: hand over to the previous handler, staying installed.
    let old = if sig == libc::SIGBUS {
        &OLD_BUS
    } else {
        &OLD_SEGV
    };
    let old = old.load(Ordering::Acquire) as *const libc::sigaction;
    let (handler, flags) = if old.is_null() {
        (libc::SIG_DFL, 0)
    } else {
        unsafe { ((*old).sa_sigaction, (*old).sa_flags) }
    };
    if handler == libc::SIG_IGN {
        return;
    }
    if handler == libc::SIG_DFL {
        // Default action: a re-run fault, or a re-raise for a sent signal
        // (blocked in this handler, so it is delivered once we return).
        unsafe {
            libc::signal(sig, libc::SIG_DFL);
            if !kernel_fault(info) {
                libc::raise(sig);
            }
        }
        return;
    }
    unsafe {
        if flags & libc::SA_SIGINFO != 0 {
            let f: extern "C" fn(libc::c_int, *mut libc::siginfo_t, *mut libc::c_void) =
                std::mem::transmute(handler);
            f(sig, info, ctx);
        } else {
            let f: extern "C" fn(libc::c_int) = std::mem::transmute(handler);
            f(sig);
        }
    }
}

/// Whether SIGSEGV is currently handled by this module's fault handler.
pub fn fault_handler_installed() -> bool {
    let mut current: libc::sigaction = unsafe { std::mem::zeroed() };
    unsafe { libc::sigaction(libc::SIGSEGV, std::ptr::null(), &mut current) };
    current.sa_sigaction == on_fault as *const () as usize
}
