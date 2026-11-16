//! Header merge of the per-chunk inputs and the MarkDuplicates `@PG` line
//! (spec §1.2-1.4, §1.6).

use std::io;

use nomics_io::bam::{Header, Reference};

fn invalid(msg: String) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg)
}

/// One parsed header line: record type and `(tag, value)` fields in order.
#[derive(Clone, Debug, PartialEq, Eq)]
struct Line {
    kind: String,
    fields: Vec<(String, String)>,
}

impl Line {
    fn parse(s: &str) -> Line {
        let mut parts = s.split('\t');
        let kind = parts.next().unwrap_or("").to_owned();
        let fields = parts
            .map(|f| match f.split_once(':') {
                Some((k, v)) => (k.to_owned(), v.to_owned()),
                None => (f.to_owned(), String::new()),
            })
            .collect();
        Line { kind, fields }
    }

    fn get(&self, key: &str) -> Option<&str> {
        self.fields
            .iter()
            .find(|(k, _)| k == key)
            .map(|(_, v)| v.as_str())
    }

    fn set(&mut self, key: &str, value: &str) {
        match self.fields.iter_mut().find(|(k, _)| k == key) {
            Some(f) => f.1 = value.to_owned(),
            None => self.fields.push((key.to_owned(), value.to_owned())),
        }
    }

    /// htsjdk `SAMTextHeaderCodec` output: ID (or SN, LN) first, then the rest in stored order.
    fn render(&self, first: &[&str]) -> String {
        let mut out = self.kind.clone();
        for key in first {
            if let Some(v) = self.get(key) {
                out.push_str(&format!("\t{key}:{v}"));
            }
        }
        for (k, v) in &self.fields {
            if !first.contains(&k.as_str()) {
                out.push_str(&format!("\t{k}:{v}"));
            }
        }
        out
    }
}

fn lines(h: &Header) -> Vec<Line> {
    String::from_utf8_lossy(&h.text)
        .trim_end_matches('\0')
        .lines()
        .filter(|l| !l.is_empty())
        .map(Line::parse)
        .collect()
}

/// `SamFileHeaderMerger.positiveFourDigitBase36Str`: plain base 36, digits 0-9A-Z.
fn base36(mut n: u32) -> String {
    const DIGITS: &[u8; 36] = b"0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    if n == 0 {
        return "0".into();
    }
    let mut s = Vec::new();
    while n > 0 {
        s.push(DIGITS[(n % 36) as usize]);
        n /= 36;
    }
    s.reverse();
    String::from_utf8(s).unwrap()
}

/// Records with equal ID and equal attributes collapse; otherwise later
/// distinct records get `ID.<base36(counter)>` (`mergeHeaderRecords`).
/// `taken` holds IDs already assigned; returns, per input record, its new ID.
fn merge_level(
    records: &[(usize, Line)],
    taken: &mut Vec<String>,
    counter: &mut u32,
    out: &mut Vec<Line>,
) -> Vec<String> {
    let mut assigned: Vec<(Line, String)> = Vec::new();
    let mut ids = Vec::with_capacity(records.len());
    for (_, rec) in records {
        if let Some((_, id)) = assigned.iter().find(|(l, _)| l == rec) {
            ids.push(id.clone());
            continue;
        }
        let id = rec.get("ID").unwrap_or("").to_owned();
        let new_id = if !taken.contains(&id) {
            *counter += 1;
            id
        } else {
            loop {
                let candidate = format!("{id}.{}", base36(*counter));
                *counter += 1;
                if !taken.contains(&candidate) {
                    break candidate;
                }
            }
        };
        taken.push(new_id.clone());
        let mut line = rec.clone();
        line.set("ID", &new_id);
        assigned.push((rec.clone(), new_id.clone()));
        out.push(line);
        ids.push(new_id);
    }
    ids
}

/// The merged header text plus reference list for the MarkDuplicates output.
pub fn merged_header(
    inputs: &[Header],
    md_version: &str,
    md_command_line: &str,
) -> io::Result<Header> {
    let first = inputs.first().ok_or_else(|| invalid("no inputs".into()))?;
    for (i, h) in inputs.iter().enumerate() {
        if h.refs != first.refs {
            return Err(invalid(format!(
                "input {i} has a different sequence dictionary"
            )));
        }
    }
    let parsed: Vec<Vec<Line>> = inputs.iter().map(lines).collect();

    // @SQ lines come from the first input.
    let sq: Vec<&Line> = parsed[0].iter().filter(|l| l.kind == "@SQ").collect();

    // @RG: merge identical records; distinct records with a clashing ID would
    // need renaming and RG-tag rewriting, which this pipeline never produces.
    let mut rg: Vec<Line> = Vec::new();
    for l in parsed.iter().flatten().filter(|l| l.kind == "@RG") {
        if !rg.contains(l) {
            if rg.iter().any(|r| r.get("ID") == l.get("ID")) {
                return Err(invalid(format!(
                    "conflicting @RG for ID {:?} (renaming not implemented)",
                    l.get("ID")
                )));
            }
            rg.push(l.clone());
        }
    }
    rg.sort_by(|a, b| a.get("ID").cmp(&b.get("ID")));

    // @PG: roots first, then children whose PP points into the previous level of
    // the same input; one counter across all levels (mergeProgramGroups).
    let mut pg_out: Vec<Line> = Vec::new();
    let mut taken: Vec<String> = Vec::new();
    let mut counter = 0u32;
    let mut translation: Vec<Vec<(String, String)>> = vec![Vec::new(); inputs.len()];
    let mut current: Vec<(usize, Line)> = parsed
        .iter()
        .enumerate()
        .flat_map(|(i, ls)| {
            ls.iter()
                .filter(|l| l.kind == "@PG" && l.get("PP").is_none())
                .map(move |l| (i, l.clone()))
        })
        .collect();
    let mut done: Vec<Vec<String>> = vec![Vec::new(); inputs.len()];
    while !current.is_empty() {
        // PP values of this level are translated through the parent's new ID.
        let translated: Vec<(usize, Line)> = current
            .iter()
            .map(|(i, l)| {
                let mut l = l.clone();
                if let Some(pp) = l.get("PP").map(str::to_owned)
                    && let Some((_, new)) = translation[*i].iter().find(|(old, _)| *old == pp)
                {
                    l.set("PP", &new.clone());
                }
                (*i, l)
            })
            .collect();
        let new_ids = merge_level(&translated, &mut taken, &mut counter, &mut pg_out);
        for ((i, l), new_id) in current.iter().zip(new_ids) {
            let old = l.get("ID").unwrap_or("").to_owned();
            translation[*i].push((old.clone(), new_id));
            done[*i].push(old);
        }
        current = parsed
            .iter()
            .enumerate()
            .flat_map(|(i, ls)| {
                let parents: Vec<String> = current
                    .iter()
                    .filter(|(j, _)| *j == i)
                    .filter_map(|(_, l)| l.get("ID").map(str::to_owned))
                    .collect();
                let done_i = done[i].clone();
                ls.iter()
                    .filter(move |l| {
                        l.kind == "@PG"
                            && l.get("PP")
                                .is_some_and(|pp| parents.iter().any(|p| p == pp))
                            && !done_i.iter().any(|d| Some(d.as_str()) == l.get("ID"))
                    })
                    .map(move |l| (i, l.clone()))
            })
            .collect();
    }
    pg_out.sort_by(|a, b| a.get("ID").cmp(&b.get("ID")));
    let md_id = if taken.iter().any(|t| t == "MarkDuplicates") {
        let mut n = pg_out.len() as u32;
        loop {
            let c = format!("MarkDuplicates.{}", base36(n));
            n += 1;
            if !taken.contains(&c) {
                break c;
            }
        }
    } else {
        "MarkDuplicates".to_owned()
    };
    let md = Line {
        kind: "@PG".into(),
        fields: vec![
            ("ID".into(), md_id),
            ("VN".into(), md_version.into()),
            ("CL".into(), md_command_line.into()),
            ("PN".into(), "MarkDuplicates".into()),
        ],
    };

    let comments: Vec<&Line> = parsed
        .iter()
        .flatten()
        .filter(|l| l.kind == "@CO")
        .collect();
    let mut text = String::from("@HD\tVN:1.6\tGO:none\tSO:coordinate\n");
    for l in sq {
        text += &l.render(&["SN", "LN"]);
        text.push('\n');
    }
    for l in &rg {
        text += &l.render(&["ID"]);
        text.push('\n');
    }
    for l in pg_out.iter().chain(std::iter::once(&md)) {
        text += &l.render(&["ID"]);
        text.push('\n');
    }
    for l in comments {
        text += &l.kind;
        for (k, v) in &l.fields {
            text += &format!(
                "\t{k}{}{v}",
                if v.is_empty() && !k.contains(':') {
                    ""
                } else {
                    ":"
                }
            );
        }
        text.push('\n');
    }
    Ok(Header {
        text: text.into_bytes(),
        refs: first
            .refs
            .iter()
            .map(|r| Reference {
                name: r.name.clone(),
                len: r.len,
            })
            .collect(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn input(chunk: &str) -> Header {
        let text = format!(
            "@HD\tVN:1.5\tSO:coordinate\n@SQ\tSN:chr1\tLN:100\n@RG\tID:L1\tPU:L1\tLB:lib\n\
             @PG\tID:bwa\tPN:bwa\tVN:0.7\tCL:bwa mem {chunk}\n@PG\tID:samtools\tPN:samtools\tPP:bwa\tVN:1.22\tCL:sort {chunk}\n"
        );
        Header {
            text: text.into_bytes(),
            refs: vec![Reference {
                name: b"chr1".to_vec(),
                len: 100,
            }],
        }
    }

    #[test]
    fn base36_digits() {
        assert_eq!(
            [0, 1, 9, 10, 11, 13, 23, 35, 36].map(base36),
            ["0", "1", "9", "A", "B", "D", "N", "Z", "10"]
        );
    }

    #[test]
    fn twelve_chunks_get_htsjdk_ids() {
        let inputs: Vec<Header> = (1..=12).map(|i| input(&format!("c{i:02}"))).collect();
        let h = merged_header(&inputs, "Version:4.6.2.0", "MarkDuplicates --INPUT x").unwrap();
        let text = String::from_utf8(h.text).unwrap();
        let ids: Vec<&str> = text
            .lines()
            .filter(|l| l.starts_with("@PG"))
            .map(|l| l.split('\t').nth(1).unwrap())
            .collect();
        let want: Vec<String> = ["bwa"]
            .into_iter()
            .map(String::from)
            .chain((1..=11).map(|n| format!("bwa.{}", base36(n))))
            .chain(["samtools".to_owned()])
            .chain((13..=23).map(|n| format!("samtools.{}", base36(n))))
            .chain(["MarkDuplicates".to_owned()])
            .map(|s| format!("ID:{s}"))
            .collect();
        assert_eq!(ids, want);
        // samtools.D is the second input's sort; its PP must point at that input's bwa (bwa.1).
        assert!(
            text.contains("@PG\tID:samtools.D\tPN:samtools\tPP:bwa.1\tVN:1.22\tCL:sort c02\n"),
            "{text}"
        );
        assert!(text.starts_with(
            "@HD\tVN:1.6\tGO:none\tSO:coordinate\n@SQ\tSN:chr1\tLN:100\n@RG\tID:L1\tPU:L1\tLB:lib\n"
        ));
        assert!(text.ends_with("@PG\tID:MarkDuplicates\tVN:Version:4.6.2.0\tCL:MarkDuplicates --INPUT x\tPN:MarkDuplicates\n"));
    }
}
