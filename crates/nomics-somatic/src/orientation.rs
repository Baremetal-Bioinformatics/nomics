//! GATK 4.6.2.0 LearnReadOrientationModel, with ordered binary64 arithmetic.
mod javafmt;
pub(crate) mod jmath;
#[allow(clippy::excessive_precision)]
mod log_table;
pub(crate) mod math;
pub mod tar;

use nomics_mem::HostArray;
use std::collections::{BTreeMap, BTreeSet};
use std::io;
use std::path::{Path, PathBuf};
use tar::invalid;
const DEPTH: usize = 200;
type Vector = [f64; 12];
#[derive(Clone, Copy)]
struct Record {
    n: u32,
    k: u32,
    x: u32,
    alt: usize,
}
#[derive(Default)]
struct Records {
    chunks: Vec<Vec<Record>>,
    len: usize,
}
impl Records {
    fn push(&mut self, r: Record) {
        // Small chunks avoid a large reallocating ordinary-heap buffer.
        if self.chunks.last().is_none_or(|v| v.len() == 256) {
            self.chunks.push(Vec::with_capacity(256));
        }
        self.chunks.last_mut().unwrap().push(r);
        self.len += 1;
    }
    fn iter(&self) -> impl Iterator<Item = Record> + '_ {
        self.chunks.iter().flat_map(|v| v.iter().copied())
    }
}
struct Sample {
    reference: BTreeMap<String, [f64; DEPTH]>,
    singleton: BTreeMap<String, [f64; DEPTH]>,
    records: [Records; 64],
    ref_files: usize,
    alt_files: usize,
    table: bool,
}
impl Default for Sample {
    fn default() -> Self {
        Self {
            reference: BTreeMap::new(),
            singleton: BTreeMap::new(),
            records: std::array::from_fn(|_| Records::default()),
            ref_files: 0,
            alt_files: 0,
            table: false,
        }
    }
}
fn base(b: u8) -> io::Result<usize> {
    match b {
        b'A' => Ok(0),
        b'C' => Ok(1),
        b'G' => Ok(2),
        b'T' => Ok(3),
        _ => Err(invalid("non-ACGT base")),
    }
}
fn context(s: &str) -> io::Result<usize> {
    if s.len() != 3 {
        return Err(invalid("expected 3-base context"));
    }
    Ok(base(s.as_bytes()[0])? * 16 + base(s.as_bytes()[1])? * 4 + base(s.as_bytes()[2])?)
}
pub(crate) fn kmer(i: usize) -> String {
    [b"ACGT"[i / 16], b"ACGT"[i / 4 % 4], b"ACGT"[i % 4]]
        .map(char::from)
        .iter()
        .collect()
}
fn rc(i: usize) -> usize {
    (3 - i % 4) * 16 + (3 - i / 4 % 4) * 4 + (3 - i / 16)
}
fn flat(reference: usize) -> Vector {
    let mut p = [0.1; 12];
    p[reference] = 0.0;
    p[reference + 4] = 0.0;
    p
}
fn integer(s: &str) -> io::Result<u32> {
    s.parse()
        .map_err(|e| invalid(format!("invalid count {s}: {e}")))
}
fn load(
    name: &str,
    data: &[u8],
    samples: &mut BTreeMap<String, Sample>,
    seen: &mut BTreeSet<(String, &'static str)>,
) -> io::Result<()> {
    let s = std::str::from_utf8(data).map_err(|e| invalid(e.to_string()))?;
    if name.ends_with(".alt_table") {
        let sample = s
            .lines()
            .find_map(|l| l.strip_prefix("#<METADATA>SAMPLE="))
            .ok_or_else(|| invalid("missing SAMPLE"))?;
        if !seen.insert((sample.to_owned(), "table")) {
            return Err(invalid(
                "multiple alt tables for one sample in one archive: filesystem order is ungated",
            ));
        }
        let sample = samples.entry(sample.to_owned()).or_default();
        sample.table = true;
        let mut lines = s.lines().filter(|l| !l.starts_with('#') && !l.is_empty());
        if lines.next() != Some("context\tref_count\talt_count\tref_f1r2\talt_f1r2\tdepth\talt") {
            return Err(invalid("unexpected alt table columns"));
        }
        for line in lines {
            let fields: Vec<_> = line.split('\t').collect();
            if fields.len() != 7 {
                return Err(invalid("malformed alt row"));
            }
            let c = context(fields[0])?;
            let n = integer(fields[1])?
                .checked_add(integer(fields[2])?)
                .ok_or_else(|| invalid("depth overflow"))?;
            let k = integer(fields[2])?;
            let x = integer(fields[4])?;
            let _ref_f1r2 = integer(fields[3])?;
            let _ignored_depth = integer(fields[5])?; // GATK recomputes ref+alt, ignoring this column.
            if x > k || fields[6].len() != 1 {
                return Err(invalid("invalid alt counts/base"));
            }
            sample.records[c].push(Record {
                n,
                k,
                x,
                alt: base(fields[6].as_bytes()[0])?,
            });
        }
    } else if name.ends_with(".ref_histogram") || name.ends_with(".alt_histogram") {
        let is_ref = name.ends_with(".ref_histogram");
        let mut lines = s.lines();
        if lines.next() != Some("## htsjdk.samtools.metrics.StringHeader") {
            return Err(invalid("expected histogram StringHeader"));
        }
        let sample = lines
            .next()
            .and_then(|l| l.strip_prefix("# "))
            .ok_or_else(|| invalid("missing histogram sample"))?;
        if !seen.insert((sample.to_owned(), if is_ref { "ref" } else { "alt" })) {
            return Err(invalid(
                "multiple histograms for one sample in one archive: filesystem order is ungated",
            ));
        }
        let sample = samples.entry(sample.to_owned()).or_default();
        if is_ref {
            sample.ref_files += 1;
        } else {
            sample.alt_files += 1;
        }
        let target = if is_ref {
            &mut sample.reference
        } else {
            &mut sample.singleton
        };
        let header = lines
            .find(|l| l.starts_with("depth\t"))
            .ok_or_else(|| invalid("missing histogram columns"))?;
        let labels: Vec<_> = header.split('\t').skip(1).collect();
        if labels.len() != if is_ref { 64 } else { 384 } {
            return Err(invalid("incomplete histograms"));
        }
        for &label in &labels {
            target.entry(label.to_owned()).or_insert([0.0; DEPTH]);
        }
        for line in lines.take_while(|l| !l.is_empty()) {
            let mut fields = line.split('\t');
            let d = integer(fields.next().unwrap())? as usize;
            if !(1..=DEPTH).contains(&d) {
                return Err(invalid("expected histogram depth 1..200"));
            }
            for &label in &labels {
                let count = fields
                    .next()
                    .ok_or_else(|| invalid("short histogram row"))?
                    .parse::<f64>()
                    .map_err(|e| invalid(e.to_string()))?;
                if count < 0.0 || !count.is_finite() {
                    return Err(invalid("invalid histogram count"));
                }
                target.get_mut(label).unwrap()[d - 1] += count;
            }
            if fields.next().is_some() {
                return Err(invalid("long histogram row"));
            }
        }
    }
    Ok(())
}
// [AF artifact, homref, germline, somatic, homvar, orientation F1, F2, balanced].
fn likelihood(r: Record) -> [f64; 8] {
    let bb = math::beta_binomial;
    [
        bb(r.n, r.k, 1.0, 9.0),
        bb(r.n, r.k, 3.0, 10000.0),
        bb(r.n, r.k, 5.0, 5.0),
        bb(r.n, r.k, 2.0, 5.0),
        bb(r.n, r.k, 10000.0, 3.0),
        bb(r.k, r.x, 100.0, 1.0),
        bb(r.k, r.x, 1.0, 100.0),
        bb(r.k, r.x, 10.0, 10.0),
    ]
}
fn responsibility(reference: usize, alt: usize, l: &[f64], log_prior: &Vector) -> Vector {
    let mut r = [f64::NEG_INFINITY; 12];
    for s in 0..12 {
        if s < 8 && (s % 4 == reference || s % 4 != alt) {
            continue;
        }
        let af = if s < 8 { 0 } else { s - 7 };
        let orientation = if s < 4 {
            5
        } else if s < 8 {
            6
        } else {
            7
        };
        r[s] = (log_prior[s] + l[af]) + l[orientation];
    }
    let mut max = 0;
    for s in 1..12 {
        if r[s] > r[max] {
            max = s;
        }
    }
    let mut sum = 1.0;
    for s in 0..12 {
        if s != max && r[s] != f64::NEG_INFINITY {
            sum += jmath::exp(r[s] - r[max]);
        }
    }
    let log_sum = r[max] + if sum != 1.0 { jmath::log(sum) } else { 0.0 };
    r.map(|v| jmath::exp(v - log_sum))
}
/// ReadOrientationFilter uses the same E step with HOM_REF excluded, then
/// takes the maximum of the two orientation states (not their sum).
pub(crate) fn artifact_probability(
    reference: usize,
    alt: usize,
    n: u32,
    k: u32,
    x: u32,
    prior: &[f64; 12],
) -> f64 {
    let l = likelihood(Record { n, k, x, alt });
    let mut log_prior = prior.map(jmath::log);
    log_prior[8] = f64::NEG_INFINITY;
    let r = responsibility(reference, alt, &l, &log_prior);
    r[alt].max(r[alt + 4])
}
fn add(target: &mut Vector, values: Vector, scale: f64) {
    for s in 0..12 {
        target[s] += values[s] * scale;
    }
}
struct Learned {
    p: Vector,
    examples: i32,
    alt_examples: i32,
    iterations: usize,
}
fn learn(
    reference: usize,
    refs: &[f64; DEPTH],
    singletons: &[(usize, usize, [f64; DEPTH])],
    records: &Records,
) -> io::Result<Learned> {
    learn_with_limit(reference, refs, singletons, records, 20)
}
fn learn_with_limit(
    reference: usize,
    refs: &[f64; DEPTH],
    singletons: &[(usize, usize, [f64; DEPTH])],
    records: &Records,
    limit: usize,
) -> io::Result<Learned> {
    let mut cache = HostArray::<f64>::zeroed(records.len * 8)?;
    for (i, r) in records.iter().enumerate() {
        cache[i * 8..i * 8 + 8].copy_from_slice(&likelihood(r));
    }
    let ref_l: Vec<_> = (1..=DEPTH)
        .map(|n| {
            likelihood(Record {
                n: n as u32,
                k: 0,
                x: 0,
                alt: reference,
            })
        })
        .collect();
    let single_l: Vec<Vec<_>> = singletons
        .iter()
        .map(|&(alt, x, _)| {
            (1..=DEPTH)
                .map(|n| {
                    likelihood(Record {
                        n: n as u32,
                        k: 1,
                        x: x as u32,
                        alt,
                    })
                })
                .collect()
        })
        .collect();
    let pseudo = flat(reference);
    let mut p = pseudo;
    let mut iterations = 0;
    loop {
        let logs = p.map(jmath::log);
        let mut design = [0.0; 12];
        for (i, r) in records.iter().enumerate() {
            add(
                &mut design,
                responsibility(reference, r.alt, &cache[i * 8..i * 8 + 8], &logs),
                1.0,
            );
        }
        let mut singles = [0.0; 12];
        for (j, (alt, _, hist)) in singletons.iter().enumerate() {
            let mut effective = [0.0; 12];
            for d in 0..DEPTH {
                add(
                    &mut effective,
                    responsibility(reference, *alt, &single_l[j][d], &logs),
                    hist[d],
                );
            }
            for s in 0..12 {
                singles[s] += effective[s];
            }
        }
        let mut reference_counts = [0.0; 12];
        for d in 0..DEPTH {
            add(
                &mut reference_counts,
                responsibility(reference, reference, &ref_l[d], &logs),
                refs[d],
            );
        }
        let mut new = [0.0; 12];
        for s in 0..12 {
            new[s] = ((design[s] + singles[s]) + reference_counts[s]) + pseudo[s];
        }
        let total: f64 = new.iter().copied().sum();
        for v in &mut new {
            *v /= total;
        }
        let mut dist = 0.0;
        for s in 0..12 {
            let diff = p[s] - new[s];
            dist += diff * diff;
        }
        p = new;
        iterations += 1;
        if dist.sqrt() <= 1e-4 || iterations >= limit {
            break;
        }
    }
    let sum = |h: &[f64; DEPTH]| h.iter().copied().sum::<f64>() as i32;
    let alt_examples = singletons
        .iter()
        .fold(records.len as i32, |n, (_, _, h)| n.wrapping_add(sum(h)));
    Ok(Learned {
        p,
        examples: alt_examples.wrapping_add(sum(refs)),
        alt_examples,
        iterations,
    })
}
pub(crate) fn hash(s: &str) -> u32 {
    let h = s
        .bytes()
        .fold(0u32, |h, b| h.wrapping_mul(31).wrapping_add(b as u32));
    h ^ (h >> 16)
}
fn table(sample_name: &str, sample: &Sample) -> io::Result<Vec<u8>> {
    if sample.ref_files == 0 || (sample.alt_files != 0 && sample.ref_files != sample.alt_files) {
        return Err(invalid("unmatched ref/alt histogram files"));
    }
    let mut priors: Vec<_> = (0..64)
        .map(|i| Learned {
            p: flat(i / 4 % 4),
            examples: 0,
            alt_examples: 0,
            iterations: 0,
        })
        .collect();
    // SequenceUtil's breadth-first generator moves AAA to the back before
    // returning. Canonical contexts keep their first encounter in this order.
    let mut canonical = Vec::new();
    for i in (1..64).chain(std::iter::once(0)) {
        let c = i.min(rc(i));
        if !canonical.contains(&c) {
            canonical.push(c);
        }
    }
    for c in canonical {
        let rev = rc(c);
        let ctx = kmer(c);
        let reverse = kmer(rev);
        let mut refs = [0.0; DEPTH];
        for (d, v) in refs.iter_mut().enumerate() {
            *v = sample.reference.get(&ctx).map_or(0.0, |h| h[d])
                + sample.reference.get(&reverse).map_or(0.0, |h| h[d]);
        }
        let mut singles = Vec::new();
        if sample.alt_files > 0 {
            for alt in 0..4 {
                if alt == c / 4 % 4 {
                    continue;
                }
                for orientation in 0..2 {
                    let label = format!(
                        "{ctx}_{}_{}",
                        b"ACGT"[alt] as char,
                        if orientation == 0 { "F1R2" } else { "F2R1" }
                    );
                    let rc_label = format!(
                        "{reverse}_{}_{}",
                        b"ACGT"[3 - alt] as char,
                        if orientation == 0 { "F2R1" } else { "F1R2" }
                    );
                    let h = std::array::from_fn(|d| {
                        sample.singleton.get(&label).map_or(0.0, |h| h[d])
                            + sample.singleton.get(&rc_label).map_or(0.0, |h| h[d])
                    });
                    singles.push((alt, usize::from(orientation == 0), h));
                }
            }
        }
        let mut records = Records::default();
        for r in sample.records[c].iter() {
            records.push(r);
        }
        for r in sample.records[rev].iter() {
            records.push(Record {
                x: r.k - r.x,
                alt: 3 - r.alt,
                ..r
            });
        }
        if refs.iter().all(|&v| v == 0.0) || records.len == 0 {
            continue;
        }
        let learned = learn(c / 4 % 4, &refs, &singles, &records)?;
        eprintln!(
            "Context {ctx}: {} ref and {} alt examples, {} EM iterations",
            learned.examples.wrapping_sub(learned.alt_examples),
            learned.alt_examples,
            learned.iterations
        );
        priors[rev] = Learned {
            p: std::array::from_fn(|s| {
                learned.p[if s < 8 {
                    (1 - s / 4) * 4 + 3 - s % 4
                } else {
                    s
                }]
            }),
            ..learned
        };
        priors[c] = learned;
    }
    let mut out = format!(
        "#<METADATA>SAMPLE={sample_name}\ncontext\trev_comp\tf1r2_a\tf1r2_c\tf1r2_g\tf1r2_t\tf2r1_a\tf2r1_c\tf2r1_g\tf2r1_t\thom_ref\tgermline_het\tsomatic_het\thom_var\tnum_examples\tnum_alt_examples\n"
    );
    // HashMap(64) resizes to 128 at 49 entries. Bucket chains retain insertion
    // order; reproduce SequenceUtil.generateAllKmers (AAA is moved to last).
    let mut order: Vec<_> = (1..64).chain(std::iter::once(0)).collect();
    order.sort_by_key(|&i| hash(&kmer(i)) & 127);
    for i in order {
        out.push_str(&format!("{}\t{}", kmer(i), kmer(rc(i))));
        for p in priors[i].p {
            out.push('\t');
            out.push_str(&javafmt::double_to_string(p));
        }
        out.push_str(&format!(
            "\t{}\t{}\n",
            priors[i].examples, priors[i].alt_examples
        ));
    }
    Ok(out.into_bytes())
}
pub(crate) fn url_encode(s: &str) -> String {
    let mut out = String::new();
    for b in s.bytes() {
        match b {
            b'A'..=b'Z' | b'a'..=b'z' | b'0'..=b'9' | b'.' | b'-' | b'*' | b'_' => {
                out.push(b as char)
            }
            b' ' => out.push('+'),
            _ => out.push_str(&format!("%{b:02X}")),
        }
    }
    out
}
/// Learn default orientation priors. Inputs are reduced in the supplied order.
/// The output gate is ordered member names and exact member bytes, not headers.
pub fn run(inputs: &[PathBuf], output: &Path) -> io::Result<()> {
    if inputs.is_empty() || !output.as_os_str().to_string_lossy().ends_with(".tar.gz") {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "need input archive(s) and an output ending .tar.gz",
        ));
    }
    let mut samples = BTreeMap::new();
    for input in inputs {
        let mut seen = BTreeSet::new();
        tar::read(input, |name, data| {
            load(name, data, &mut samples, &mut seen)
        })?;
    }
    let any_alt = samples.values().any(|s| s.alt_files > 0);
    if samples
        .values()
        .any(|s| any_alt && (s.ref_files == 0 || s.ref_files != s.alt_files))
    {
        return Err(invalid(
            "ref and alt histogram sample sets/counts must agree",
        ));
    }
    // GATK serializes several samples through temp-dir listFiles(), whose
    // ordering is not specified by Java and differs between filesystems.
    // Refuse that ungated regime rather than inventing an archive member order.
    if samples.values().filter(|s| s.table).count() > 1 {
        return Err(invalid(
            "multiple output samples are not gated: GATK member order depends on the filesystem",
        ));
    }
    let mut members = Vec::new();
    for (name, sample) in samples {
        if sample.table {
            members.push((
                format!("./{}.orientation_priors", url_encode(&name)),
                table(&name, &sample)?,
            ));
        }
    }
    tar::write(output, &members)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn java_arithmetic_oracles() {
        // Independent Java harness, compiled with javac 17 against the pinned
        // Commons Math source (no Rust calculations used to construct values).
        // Oracle recipe and arithmetic are recorded in the spec, section 1.5.
        for line in include_str!("orientation/oracle.txt").lines() {
            if line.starts_with("bb ") {
                let f: Vec<_> = line.split_whitespace().skip(1).collect();
                let got = math::beta_binomial(
                    f[0].parse().unwrap(),
                    f[1].parse().unwrap(),
                    f[2].parse().unwrap(),
                    f[3].parse().unwrap(),
                );
                assert_eq!(
                    got.to_bits(),
                    u64::from_str_radix(f[4], 16).unwrap(),
                    "{line}"
                );
            } else if line.starts_with("log ")
                || line.starts_with("log1p ")
                || line.starts_with("beta ")
            {
                let f: Vec<_> = line.split_whitespace().collect();
                let input = || f64::from_bits(u64::from_str_radix(f[1], 16).unwrap());
                let got = match f[0] {
                    "log" => math::log(input()),
                    "log1p" => math::log1p(input()),
                    _ => math::log_beta(f[1].parse().unwrap(), f[2].parse().unwrap()),
                };
                assert_eq!(
                    got.to_bits(),
                    u64::from_str_radix(f[f.len() - 1], 16).unwrap(),
                    "{line}"
                );
            } else {
                let end = line.find(']').unwrap();
                let f: Vec<u32> = line[6..end]
                    .split(", ")
                    .map(|s| s.parse().unwrap())
                    .collect();
                let r = Record {
                    n: f[0],
                    k: f[1],
                    x: f[2],
                    alt: if f[1] == 0 { 1 } else { 0 },
                };
                let got = responsibility(1, r.alt, &likelihood(r), &flat(1).map(jmath::log));
                for (v, bits) in got.iter().zip(line[end + 1..].split_whitespace()) {
                    assert_eq!(
                        v.to_bits(),
                        u64::from_str_radix(bits, 16).unwrap(),
                        "{line}"
                    );
                }
            }
        }
    }
    #[test]
    fn independent_java_em_trajectory() {
        let mut refs = [0.0; DEPTH];
        refs[29] = 1000.0;
        let mut singletons = Vec::new();
        for alt in [0, 2, 3] {
            for o in 0..2 {
                let mut h = [0.0; DEPTH];
                if alt == 0 && o == 0 {
                    h[19] = 7.0;
                }
                singletons.push((alt, usize::from(o == 0), h));
            }
        }
        let mut records = Records::default();
        for (n, k, x) in [(30, 2, 2), (30, 2, 0), (100, 50, 23)] {
            records.push(Record { n, k, x, alt: 0 });
        }
        let lines: Vec<_> = include_str!("orientation/em_oracle.txt").lines().collect();
        for (i, line) in lines.iter().enumerate() {
            let got = learn_with_limit(1, &refs, &singletons, &records, i + 1).unwrap();
            assert_eq!(got.iterations, i + 1);
            assert_eq!(got.examples, 1010);
            assert_eq!(got.alt_examples, 10);
            for (actual, bits) in got.p.iter().zip(line.split_whitespace().skip(2)) {
                assert_eq!(
                    actual.to_bits(),
                    u64::from_str_radix(bits, 16).unwrap(),
                    "{line}"
                );
            }
        }
        assert_eq!(
            learn(1, &refs, &singletons, &records).unwrap().iterations,
            lines.len()
        );
    }
    #[test]
    fn singleton_only_context_keeps_flat_default_and_zero_counts() {
        let mut sample = Sample {
            ref_files: 1,
            alt_files: 1,
            table: true,
            ..Sample::default()
        };
        let mut h = [0.0; DEPTH];
        h[29] = 1000.0;
        sample.reference.insert("ACA".to_owned(), h);
        sample.singleton.insert("ACA_A_F1R2".to_owned(), h);
        let bytes = table("test sample", &sample).unwrap();
        let text = std::str::from_utf8(&bytes).unwrap();
        let row: Vec<_> = text
            .lines()
            .find(|l| l.starts_with("ACA\t"))
            .unwrap()
            .split('\t')
            .collect();
        assert_eq!(
            &row[2..14],
            [
                "0.1", "0.0", "0.1", "0.1", "0.1", "0.0", "0.1", "0.1", "0.1", "0.1", "0.1", "0.1"
            ]
        );
        assert_eq!(&row[14..], ["0", "0"]);
        assert_eq!(url_encode("A B/é"), "A+B%2F%C3%A9");
    }
    #[test]
    fn hashmap_order_and_reverse_complement() {
        let mut order: Vec<_> = (1..64).chain(std::iter::once(0)).collect();
        order.sort_by_key(|&i| hash(&kmer(i)) & 127);
        assert_eq!(
            order[..7].iter().map(|&i| kmer(i)).collect::<Vec<_>>(),
            ["ATT", "CTT", "GTT", "TAT", "AAA", "CAA", "AAC"]
        );
        for i in 0..64 {
            assert_eq!(rc(rc(i)), i);
        }
    }
}
