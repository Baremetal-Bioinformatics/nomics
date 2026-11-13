//! `bwa` 0.7.19-r1273. Only `bwa mem` with sarek's options is handled:
//! `-K 100000000 -Y -R <@RG line> -t <n> <index prefix> <r1.fq.gz> <r2.fq.gz>`,
//! SAM on stdout. Options are parsed as `main_mem`'s `getopt` does (glibc:
//! options may follow positionals; values attached or separate).
//!
//! Output (docs/spec/bwa.md §0.3): `@HD VN:1.5 SO:unsorted GO:query`, one
//! `@SQ` per `.ann` contig, the unescaped `-R` line, then
//! `@PG ID:bwa PN:bwa VN:0.7.19-r1273 CL:<argv joined by spaces>` — the same
//! header bwa writes, so a downstream `samtools sort` keeps its usual
//! provenance; then the alignments, batch by batch.

use std::io::{BufReader, Read, Write};
use std::path::{Path, PathBuf};

use super::args::Unsupported;
use super::{Decision, threads, with_suffix};

/// Options taking a value in bwa 0.7.19 `main_mem`'s getopt string
/// `"51qpaMCSPVYjuk:c:v:s:r:t:R:A:B:O:E:U:w:L:d:T:Q:D:m:I:N:o:f:W:x:G:h:y:K:X:H:F:z:"`;
/// the other letters there are flags.
const WITH_VALUE: &str = "kcvsrtRABOEUwLdTQDmINofWxGhyKXHFz";

struct Plan {
    index: PathBuf,
    r1: PathBuf,
    r2: PathBuf,
    rg_line: String,
    rg_id: String,
    threads: Option<usize>,
    pg: String,
}

pub fn decide(args: &[String]) -> Decision {
    let label = match args.first() {
        Some(sub) if !sub.starts_with('-') => format!("bwa {sub}"),
        _ => "bwa".to_owned(),
    };
    match plan(args) {
        Ok(p) => Decision::Handle {
            label,
            run: Box::new(move || run(p)),
        },
        Err(why) => Decision::Delegate { label, why },
    }
}

/// `bwa_escape`: `\t`, `\n`, `\r`, `\\` -> the characters.
fn unescape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    let mut it = s.chars();
    while let Some(c) = it.next() {
        if c == '\\' {
            match it.next() {
                Some('t') => out.push('\t'),
                Some('n') => out.push('\n'),
                Some('r') => out.push('\r'),
                Some('\\') => out.push('\\'),
                Some(o) => {
                    out.push('\\');
                    out.push(o);
                }
                None => out.push('\\'),
            }
        } else {
            out.push(c);
        }
    }
    out
}

fn plan(args: &[String]) -> Result<Plan, Unsupported> {
    match args.first().map(String::as_str) {
        Some("mem") => {}
        Some(sub) => return Err(Unsupported::new(format!("bwa {sub}"))),
        None => return Err(Unsupported::new("usage/version query")),
    }
    let rest = &args[1..];
    let (mut k, mut y, mut rg, mut t) = (None, false, None, None);
    let mut positional = Vec::new();
    let mut i = 0;
    while i < rest.len() {
        let a = &rest[i];
        i += 1;
        if a == "--" {
            positional.extend(rest[i..].iter().cloned());
            break;
        }
        if !a.starts_with('-') || a == "-" {
            positional.push(a.clone());
            continue;
        }
        // An option cluster, e.g. `-Y`, `-t20`, `-YK100000000`.
        let chars: Vec<char> = a[1..].chars().collect();
        let mut j = 0;
        while j < chars.len() {
            let c = chars[j];
            j += 1;
            if WITH_VALUE.contains(c) {
                let v: String = if j < chars.len() {
                    let v = chars[j..].iter().collect();
                    j = chars.len();
                    v
                } else {
                    let Some(v) = rest.get(i) else {
                        return Err(Unsupported::new(format!("-{c} without a value")));
                    };
                    i += 1;
                    v.clone()
                };
                match c {
                    'K' => k = Some(v),
                    'R' => rg = Some(v),
                    't' => t = Some(v),
                    _ => return Err(Unsupported::new(format!("-{c}"))),
                }
            } else if c == 'Y' {
                y = true;
            } else {
                return Err(Unsupported::new(format!("-{c}")));
            }
        }
    }
    if k.as_deref() != Some("100000000") {
        return Err(Unsupported::new(match k {
            Some(v) => format!("-K {v}"),
            None => "no -K (batch size depends on -t)".into(),
        }));
    }
    if !y {
        return Err(Unsupported::new("no -Y (hard-clipped supplementaries)"));
    }
    let threads = match t {
        None => None,
        Some(v) => Some(
            v.parse::<usize>()
                .map_err(|_| Unsupported::new(format!("-t {v}")))?,
        ),
    };
    let Some(rg) = rg else {
        return Err(Unsupported::new("no -R"));
    };
    let rg_line = unescape(&rg);
    if !rg_line.starts_with("@RG") || rg.contains('\t') {
        return Err(Unsupported::new("-R not an escaped @RG line"));
    }
    let rg_id = rg_line
        .split('\t')
        .find_map(|f| f.strip_prefix("ID:"))
        .ok_or_else(|| Unsupported::new("-R without ID"))?
        .to_owned();
    let [index, r1, r2] = positional.as_slice() else {
        return Err(Unsupported::new(match positional.len() {
            2 => "single-end or interleaved input".to_owned(),
            n => format!("{n} positional arguments"),
        }));
    };
    // bwa_idx_infer_prefix: <hint>.64.bwt, else <hint>.bwt.
    let hint = Path::new(index);
    let prefix = if with_suffix(hint, ".64.bwt").exists() {
        with_suffix(hint, ".64")
    } else if with_suffix(hint, ".bwt").exists() {
        hint.to_path_buf()
    } else {
        return Err(Unsupported::new(format!("no BWA index at {index}")));
    };
    if with_suffix(&prefix, ".alt").exists() {
        return Err(Unsupported::new("ALT contigs (.alt file)"));
    }
    for r in [r1, r2] {
        let mut magic = [0u8; 2];
        let gz = std::fs::File::open(r)
            .and_then(|mut f| f.read_exact(&mut magic))
            .is_ok()
            && magic == [0x1f, 0x8b];
        if !gz {
            return Err(Unsupported::new(format!("{r} is not gzip FASTQ")));
        }
    }
    // main.c: "@PG\tID:bwa\tPN:bwa\tVN:" PACKAGE_VERSION "\tCL:" argv joined by spaces.
    let argv0 = std::env::args().next().unwrap_or_else(|| "bwa".into());
    let mut pg = format!("@PG\tID:bwa\tPN:bwa\tVN:0.7.19-r1273\tCL:{argv0}");
    for a in args {
        pg.push(' ');
        pg.push_str(a);
    }
    Ok(Plan {
        index: hint.to_path_buf(),
        r1: PathBuf::from(r1),
        r2: PathBuf::from(r2),
        rg_line,
        rg_id,
        threads,
        pg,
    })
}

fn run(p: Plan) -> Result<(), String> {
    let t = std::time::Instant::now();
    let idx = nomics_align::index::Index::load(&p.index)
        .map_err(|e| format!("index {}: {e}", p.index.display()))?;
    eprintln!(
        "[nomics-shim bwa] index loaded in {:.1}s",
        t.elapsed().as_secs_f64()
    );
    let stdout = std::io::stdout();
    let mut out = std::io::BufWriter::with_capacity(1 << 22, stdout.lock());
    let mut header = String::from("@HD\tVN:1.5\tSO:unsorted\tGO:query\n");
    for a in &idx.bns.anns {
        header.push_str(&format!("@SQ\tSN:{}\tLN:{}\n", a.name, a.len));
    }
    header.push_str(&p.rg_line);
    header.push('\n');
    header.push_str(&p.pg);
    header.push('\n');
    out.write_all(header.as_bytes()).map_err(|e| e.to_string())?;
    let open = |path: &Path| -> Result<_, String> {
        nomics_io::gzip::GzReader::open(path)
            .map(|r| BufReader::with_capacity(1 << 20, r))
            .map_err(|e| format!("{}: {e}", path.display()))
    };
    let (mut r1, mut r2) = (open(&p.r1)?, open(&p.r2)?);
    let pairs = nomics_align::mem::align_stream_sam_each(
        &idx,
        p.rg_id.as_bytes(),
        &mut r1,
        &mut r2,
        threads(p.threads),
        |sam, log| {
            for l in log {
                eprintln!("{l}");
            }
            out.write_all(&sam)
        },
    )
    .map_err(|e| e.to_string())?;
    out.flush().map_err(|e| e.to_string())?;
    eprintln!(
        "[nomics-shim bwa] {pairs} pairs in {:.1}s",
        t.elapsed().as_secs_f64()
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unescapes_read_group() {
        assert_eq!(unescape(r"@RG\tID:a\tSM:b"), "@RG\tID:a\tSM:b");
    }

    #[test]
    fn options_outside_the_set_delegate() {
        let v = |s: &str| s.split_whitespace().map(str::to_owned).collect::<Vec<_>>();
        assert_eq!(
            plan(&v("mem -K 100000000 -Y -M -R x idx a b")).err().unwrap().0,
            "-M"
        );
        assert_eq!(plan(&v("index ref.fa")).err().unwrap().0, "bwa index");
        assert!(plan(&v("mem -K 10 -Y idx a b")).is_err());
        assert!(plan(&[]).is_err());
    }
}
