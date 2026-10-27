//! GATK 4.6.2.0 read-threading local assembly.
//!
//! Gate: `tests/assembly.rs` (all tier0 regions vs `GatkAssemble.java`).

pub(crate) mod cigar;
mod graph;
mod sw;

#[doc(hidden)]
pub use sw::bench as sw_bench;

use cigar::{calculate_cigar, cigar_reference_len, hard_clip_soft_clips};
use graph::assemble_graph;

/// A finalized read supplied to local assembly.
#[derive(Clone, Debug)]
pub struct AsmRead {
    pub name: Vec<u8>,
    pub flag: u16,
    pub start: i64,
    /// BAM encoded (`length << 4 | op`).
    pub cigar: Vec<u32>,
    pub bases: Vec<u8>,
    pub quals: Vec<u8>,
}

/// One entry of `AssemblyResultSet.getHaplotypeList()`.
#[derive(Clone, Debug)]
pub struct Haplotype {
    pub bases: Vec<u8>,
    pub is_ref: bool,
    pub alignment_start_hap_wrt_ref: i32,
    pub cigar: Vec<u32>,
    pub contig: String,
    pub start: i64,
    pub end: i64,
    pub score: f64,
    pub kmer_size: i32,
}

/// Assemble the finalized reads for an active region using HaplotypeCaller's
/// default read-threading configuration.
pub fn assemble(
    contig: &str,
    region_start: i64,
    region_end: i64,
    padded_ref_start: i64,
    padded_ref: &[u8],
    reads: &[AsmRead],
) -> Vec<Haplotype> {
    let padded_start = region_start.saturating_sub(100).max(1);
    let padded_end = region_end + 100;
    let from = usize::try_from(padded_start - padded_ref_start)
        .expect("padded reference begins after assembly span");
    let wanted = usize::try_from(padded_end - padded_start + 1).expect("invalid assembly span");
    let to = from.saturating_add(wanted).min(padded_ref.len());
    let ref_bases = padded_ref[from..to].to_vec();
    let actual_end = padded_start + i64::try_from(ref_bases.len()).unwrap_or(0) - 1;
    let alignment_start =
        i32::try_from(padded_start - padded_ref_start).expect("reference offset does not fit i32");

    let ref_haplotype = Haplotype {
        bases: ref_bases.clone(),
        is_ref: true,
        alignment_start_hap_wrt_ref: alignment_start,
        cigar: vec![(u32::try_from(ref_bases.len()).expect("reference too long") << 4)],
        contig: contig.to_owned(),
        start: padded_start,
        end: actual_end,
        score: f64::NAN,
        kmer_size: 0,
    };

    // runLocalAssembly hard-clips soft clips once more immediately before
    // adding reads to the read-threading graph.
    let clipped: Vec<(Vec<u8>, Vec<u8>)> = reads
        .iter()
        .map(|read| hard_clip_soft_clips(&read.bases, &read.quals, &read.cigar))
        .collect();

    let mut result: Vec<Haplotype> = Vec::new();
    let mut any_graph = false;
    let mut kmer_sizes = vec![10_usize, 25];
    let mut next_extra = 35_usize;
    let mut attempt = 0_usize;
    while attempt < kmer_sizes.len() {
        let kmer = kmer_sizes[attempt];
        attempt += 1;
        // Expansion attempts are 35,45,55,65,75,85; only the sixth
        // expansion (k=85), not the sixth graph overall, relaxes both the
        // low-complexity and non-unique-reference checks.
        let allow_low_complexity = kmer == 85;
        let graph = assemble_graph(&ref_bases, &clipped, kmer, allow_low_complexity);
        if debug_assembly() {
            eprintln!(
                "  [asm {contig}:{padded_start}] k={kmer} -> {}",
                graph.as_ref().map_or("no graph".to_owned(), |p| format!("{} paths", p.len()))
            );
        }
        if let Some(paths) = graph {
            any_graph = true;
            for path in paths {
                if debug_assembly() {
                    eprintln!(
                        "    [path {contig}:{padded_start}] {:?} {}",
                        path.score,
                        String::from_utf8_lossy(&path.bases)
                    );
                }
                if result.iter().any(|h| h.bases == path.bases) || path.bases == ref_bases {
                    continue;
                }
                let Some(cigar) = calculate_cigar(&ref_bases, &path.bases) else {
                    continue;
                };
                // pathIsTooDivergentFromReference: a CIGAR `N` (reference
                // skip) operator, not `N` bases. Paths through reference `N`s
                // (padding that reaches an assembly gap) are kept.
                if cigar.iter().any(|&c| c & 0xf == 3) {
                    continue;
                }
                let reference_len = cigar_reference_len(&cigar);
                if reference_len < 30 || reference_len != ref_bases.len() {
                    continue;
                }
                result.push(Haplotype {
                    bases: path.bases,
                    is_ref: false,
                    alignment_start_hap_wrt_ref: alignment_start,
                    cigar,
                    contig: contig.to_owned(),
                    start: padded_start,
                    end: actual_end,
                    score: path.score,
                    kmer_size: i32::try_from(kmer).expect("kmer too large"),
                });
            }
            // The true reference replaces the non-reference reference path at
            // the first successful graph boundary and retains that position.
            if !result.iter().any(|h| h.is_ref) {
                result.push(ref_haplotype.clone());
            }
        }

        if attempt == kmer_sizes.len() && !any_graph && next_extra <= 85 {
            kmer_sizes.push(next_extra);
            next_extra += 10;
        }
    }

    if result.iter().all(|h| !h.is_ref) {
        result.push(ref_haplotype);
    }
    result
}

/// `NOMICS_DEBUG_ASSEMBLY=1` prints each k-mer attempt per region.
pub(super) fn debug_assembly() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("NOMICS_DEBUG_ASSEMBLY").is_some())
}
