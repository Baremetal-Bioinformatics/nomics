//! Sarek's interval scatter: `CREATE_INTERVALS_BED` splits the `--intervals`
//! BED into groups, and BaseRecalibrator, ApplyBQSR and HaplotypeCaller run
//! once per group. The gather and merge steps take the groups in Nextflow's
//! task-completion order, which is recorded from the reference run.

use std::collections::HashMap;
use std::io;
use std::path::{Path, PathBuf};

/// One interval group: the BED file name sarek gives it and its lines.
pub struct Group {
    pub name: String,
    pub lines: Vec<String>,
}

/// `CREATE_INTERVALS_BED` (nf-core/sarek 3.10.0), an awk script over the BED:
///
/// ```awk
/// t = $5; if (t == "") t = ($3 - $2) / nucleotides_per_second
/// if (name == "" || (chunk > 600 && (chunk + t) > longest * 1.05)) {
///     name = sprintf("%s_%d-%d.bed", $1, $2+1, $3); chunk = 0; longest = 0 }
/// if (t > longest) longest = t
/// chunk += t
/// print $0 > name
/// ```
///
/// awk does this arithmetic in doubles, as here. Lines sent to a name that
/// was already used are appended to that file.
pub fn split_bed(text: &str, nucleotides_per_second: f64) -> io::Result<Vec<Group>> {
    let mut groups: Vec<Group> = Vec::new();
    let mut by_name: HashMap<String, usize> = HashMap::new();
    let mut current: Option<usize> = None;
    let (mut chunk, mut longest) = (0.0f64, 0.0f64);
    for line in text.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        let num = |i: usize| -> io::Result<f64> {
            let s = f.get(i).copied().unwrap_or("");
            s.trim().parse::<f64>().map_err(|_| {
                io::Error::new(io::ErrorKind::InvalidData, format!("bad BED line: {line}"))
            })
        };
        let t = match f.get(4).copied().unwrap_or("") {
            "" => (num(2)? - num(1)?) / nucleotides_per_second,
            _ => num(4)?,
        };
        if current.is_none() || (chunk > 600.0 && chunk + t > longest * 1.05) {
            let name = format!("{}_{}-{}.bed", f[0], num(1)? as i64 + 1, num(2)? as i64);
            let next = groups.len();
            let at = *by_name.entry(name.clone()).or_insert(next);
            if at == next {
                groups.push(Group {
                    name,
                    lines: Vec::new(),
                });
            }
            current = Some(at);
            chunk = 0.0;
            longest = 0.0;
        }
        if t > longest {
            longest = t;
        }
        chunk += t;
        groups[current.expect("set above")]
            .lines
            .push(line.to_owned());
    }
    Ok(groups)
}

/// Write each group as `<dir>/<name>` and return the paths.
pub fn write_groups(groups: &[Group], dir: &Path) -> io::Result<Vec<PathBuf>> {
    std::fs::create_dir_all(dir)?;
    groups
        .iter()
        .map(|g| {
            let path = dir.join(&g.name);
            let mut text = g.lines.join("\n");
            text.push('\n');
            std::fs::write(&path, text)?;
            Ok(path)
        })
        .collect()
}

/// How the reference run gathered and merged the groups.
pub struct RunOrder {
    pub gather: Vec<usize>,
    pub merge: Vec<usize>,
    /// samtools merge's `time(NULL)` seed, when recorded.
    pub merge_seed: Option<i64>,
}

/// The reference run's gather and merge orders, from a file of
/// `gather <group name>` and `merge <group name>` lines in order, plus an
/// optional `merge-seed <unix seconds>`. Without a file both orders are the
/// groups' BED order.
pub fn orders(groups: &[Group], path: Option<&Path>) -> io::Result<RunOrder> {
    let identity: Vec<usize> = (0..groups.len()).collect();
    let Some(path) = path else {
        return Ok(RunOrder {
            gather: identity.clone(),
            merge: identity,
            merge_seed: None,
        });
    };
    let index: HashMap<&str, usize> = groups
        .iter()
        .enumerate()
        .map(|(i, g)| (g.name.as_str(), i))
        .collect();
    let (mut gather, mut merge, mut merge_seed) = (Vec::new(), Vec::new(), None);
    for line in std::fs::read_to_string(path)?.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let bad = || {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("bad scatter order line: {line}"),
            )
        };
        let (kind, name) = line.split_once(char::is_whitespace).ok_or_else(bad)?;
        if kind == "merge-seed" {
            merge_seed = Some(name.trim().parse::<i64>().map_err(|_| bad())?);
            continue;
        }
        let name = name.trim().trim_end_matches(".bed");
        let i = *index.get(format!("{name}.bed").as_str()).ok_or_else(bad)?;
        match kind {
            "gather" => gather.push(i),
            "merge" => merge.push(i),
            _ => return Err(bad()),
        }
    }
    for order in [&gather, &merge] {
        let mut sorted = order.clone();
        sorted.sort_unstable();
        if sorted != identity {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("{}: each order must list every group once", path.display()),
            ));
        }
    }
    Ok(RunOrder {
        gather,
        merge,
        merge_seed,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn splits_like_the_awk_script() {
        // 100 Mb per line at 200,000 nt/s is 500 s: the second line takes the
        // chunk past 600 s, and the third starts a new group because
        // 1000 + 500 > 500 * 1.05.
        let bed = "chr1\t0\t100000000\nchr1\t100000000\t200000000\nchr2\t9\t100000009\n";
        let g = split_bed(bed, 200000.0).unwrap();
        assert_eq!(g.len(), 2);
        assert_eq!(g[0].name, "chr1_1-100000000.bed");
        assert_eq!(g[0].lines.len(), 2);
        assert_eq!(g[1].name, "chr2_10-100000009.bed");
    }

    #[test]
    fn small_beds_are_one_group() {
        let bed = "chr1\t10\t20\nchr2\t10\t20\n";
        let g = split_bed(bed, 200000.0).unwrap();
        assert_eq!(g.len(), 1);
        assert_eq!(g[0].lines.len(), 2);
    }
}
