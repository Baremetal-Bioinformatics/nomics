//! Option tables and parsers shared by the shims.
//!
//! Each tool's supported set is a table of [`Spec`]s. Anything outside the
//! table (an unknown option, an abbreviation, a positional where none is
//! accepted, a repeated scalar, a value outside the accepted list) yields an
//! [`Unsupported`] naming it, and the caller delegates.

use std::collections::BTreeMap;
use std::fmt;

/// Why an invocation is delegated: the first item outside the supported set.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Unsupported(pub String);

impl Unsupported {
    pub fn new(s: impl Into<String>) -> Self {
        Unsupported(s.into())
    }
}

impl fmt::Display for Unsupported {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Kind {
    /// Takes one value; may be given once.
    Value,
    /// Takes one value; may be repeated (a collection argument).
    Multi,
    /// Barclay boolean: `--flag`, `--flag true|false` or `--flag=true|false`.
    Bool,
    /// getopt flag without a value (fastp).
    Flag,
}

pub struct Spec {
    /// Accepted spellings including dashes, e.g. `["--INPUT", "-I"]`.
    pub names: &'static [&'static str],
    pub kind: Kind,
    /// Key in the parsed map (the canonical long name).
    pub key: &'static str,
}

pub const fn spec(names: &'static [&'static str], kind: Kind, key: &'static str) -> Spec {
    Spec { names, kind, key }
}

/// Parsed options: key -> values in command-line order (`"true"`/`"false"`
/// for booleans and flags).
#[derive(Default, Debug)]
pub struct Parsed {
    pub map: BTreeMap<&'static str, Vec<String>>,
}

impl Parsed {
    pub fn one(&self, key: &str) -> Option<&str> {
        self.map.get(key).and_then(|v| v.first()).map(String::as_str)
    }

    pub fn all(&self, key: &str) -> &[String] {
        self.map.get(key).map_or(&[], Vec::as_slice)
    }

    pub fn has(&self, key: &str) -> bool {
        self.map.contains_key(key)
    }

    /// A required value.
    pub fn req(&self, key: &str) -> Result<&str, Unsupported> {
        self.one(key)
            .ok_or_else(|| Unsupported::new(format!("missing --{key}")))
    }

    /// A boolean; `default` when absent.
    pub fn flag(&self, key: &str, default: bool) -> bool {
        self.one(key).map_or(default, |v| v == "true")
    }

    /// Require `key`, when given, to be one of `allowed` (exact spelling).
    pub fn only(&self, key: &str, allowed: &[&str]) -> Result<(), Unsupported> {
        match self.one(key) {
            Some(v) if !allowed.contains(&v) => Err(Unsupported::new(format!("--{key} {v}"))),
            _ => Ok(()),
        }
    }
}

fn find<'a>(table: &'a [Spec], name: &str) -> Option<&'a Spec> {
    table.iter().find(|s| s.names.contains(&name))
}

fn parse_bool(v: &str) -> Option<&'static str> {
    if v.eq_ignore_ascii_case("true") {
        Some("true")
    } else if v.eq_ignore_ascii_case("false") {
        Some("false")
    } else {
        None
    }
}

/// Parse `args` against `table`. With `getopt`, short options may carry
/// their value attached (`-w12`), as getopt allows; otherwise (Barclay)
/// a value is always the next token or follows `=` on a long name.
pub fn parse(args: &[String], table: &[Spec], getopt: bool) -> Result<Parsed, Unsupported> {
    let mut out = Parsed::default();
    let mut i = 0;
    while i < args.len() {
        let tok = &args[i];
        i += 1;
        if !tok.starts_with('-') || tok == "-" {
            return Err(Unsupported::new(format!("positional argument {tok:?}")));
        }
        let (name, inline): (&str, Option<String>) = if tok.starts_with("--") {
            match tok.split_once('=') {
                Some((n, v)) => (n, Some(v.to_owned())),
                None => (tok.as_str(), None),
            }
        } else if getopt && tok.len() > 2 && find(table, &tok[..2]).is_some() {
            (&tok[..2], Some(tok[2..].to_owned()))
        } else {
            (tok.as_str(), None)
        };
        let Some(s) = find(table, name) else {
            return Err(Unsupported::new(name.to_owned()));
        };
        let value = match s.kind {
            Kind::Value | Kind::Multi => match inline {
                Some(v) => v,
                None => {
                    let Some(v) = args.get(i) else {
                        return Err(Unsupported::new(format!("{name} without a value")));
                    };
                    i += 1;
                    v.clone()
                }
            },
            Kind::Bool => match inline {
                Some(v) => parse_bool(&v)
                    .ok_or_else(|| Unsupported::new(format!("{name}={v}")))?
                    .to_owned(),
                None => match args.get(i).and_then(|v| parse_bool(v)) {
                    Some(b) => {
                        i += 1;
                        b.to_owned()
                    }
                    None => "true".to_owned(),
                },
            },
            Kind::Flag => {
                if inline.is_some() {
                    return Err(Unsupported::new(format!("{name} with a value")));
                }
                "true".to_owned()
            }
        };
        let slot = out.map.entry(s.key).or_default();
        if !slot.is_empty() && s.kind != Kind::Multi {
            return Err(Unsupported::new(format!("{name} given more than once")));
        }
        slot.push(value);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    const T: &[Spec] = &[
        spec(&["--INPUT", "-I"], Kind::Multi, "INPUT"),
        spec(&["--OUTPUT", "-O"], Kind::Value, "OUTPUT"),
        spec(&["--CREATE_INDEX"], Kind::Bool, "CREATE_INDEX"),
        spec(&["--thread", "-w"], Kind::Value, "thread"),
        spec(&["--detect_adapter_for_pe"], Kind::Flag, "detect_adapter_for_pe"),
    ];

    fn v(s: &[&str]) -> Vec<String> {
        s.iter().map(|x| x.to_string()).collect()
    }

    #[test]
    fn barclay_spellings() {
        let p = parse(
            &v(&["--INPUT", "a.bam", "-I", "b.bam", "--OUTPUT=o.bam", "--CREATE_INDEX", "true"]),
            T,
            false,
        )
        .unwrap();
        assert_eq!(p.all("INPUT"), ["a.bam", "b.bam"]);
        assert_eq!(p.one("OUTPUT"), Some("o.bam"));
        assert!(p.flag("CREATE_INDEX", false));
        let p = parse(&v(&["--CREATE_INDEX", "--OUTPUT", "o"]), T, false).unwrap();
        assert!(p.flag("CREATE_INDEX", false));
        let p = parse(&v(&["--CREATE_INDEX", "FALSE"]), T, false).unwrap();
        assert!(!p.flag("CREATE_INDEX", true));
    }

    #[test]
    fn rejects_outside_the_table() {
        assert_eq!(
            parse(&v(&["--OUTPUT", "a", "--OUTPUT", "b"]), T, false).unwrap_err(),
            Unsupported::new("--OUTPUT given more than once")
        );
        assert_eq!(
            parse(&v(&["--ASSUME_SORT_ORDER", "x"]), T, false).unwrap_err(),
            Unsupported::new("--ASSUME_SORT_ORDER")
        );
        assert!(parse(&v(&["x.bam"]), T, false).is_err());
        assert!(parse(&v(&["--OUTPUT"]), T, false).is_err());
        // Barclay does not split short options; getopt does.
        assert!(parse(&v(&["-w12"]), T, false).is_err());
        assert_eq!(parse(&v(&["-w12"]), T, true).unwrap().one("thread"), Some("12"));
        assert!(parse(&v(&["--detect_adapter_for_pe=1"]), T, true).is_err());
    }
}
