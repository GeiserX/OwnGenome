//! The command lines of `og couple`, read as argparse read the Python's:
//! `--name VALUE` or `--name=VALUE`, repeated options kept in order. An
//! option the verb does not take is refused before anything is read.

use std::collections::HashMap;
use std::io::Write;

/// How an option takes its value.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Kind {
    /// One value; the last one given wins.
    One,
    /// One value each time, every one kept in order.
    Append,
}

/// What a verb takes.
pub struct Spec {
    pub prog: &'static str,
    pub usage: &'static str,
    pub options: &'static [(&'static str, Kind)],
    pub required: &'static [&'static str],
}

/// A parsed command line.
#[derive(Debug, Default)]
pub struct Args {
    vals: HashMap<&'static str, Vec<String>>,
}

impl Args {
    /// The value of a `One` option.
    pub fn one(&self, name: &str) -> Option<&str> {
        self.vals
            .get(name)
            .and_then(|v| v.last())
            .map(String::as_str)
    }

    /// Every value of an `Append` option (none when absent).
    pub fn many(&self, name: &str) -> &[String] {
        self.vals.get(name).map_or(&[], Vec::as_slice)
    }
}

/// A value, not an option: argparse takes `-` and negative numbers as values.
fn takes_value(a: &str) -> bool {
    a == "-" || !a.starts_with('-') || a[1..].parse::<f64>().is_ok()
}

/// Parses `args` by `spec`. `Err(code)` once a refusal (exit 2) or the
/// usage (`--help`, exit 0) is printed.
pub fn parse(
    spec: &Spec,
    args: &[String],
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> Result<Args, u8> {
    let mut refuse = |msg: String| -> u8 {
        let _ = write!(err, "{}", spec.usage);
        let _ = writeln!(err, "{}: error: {msg}", spec.prog);
        2
    };
    let mut got = Args::default();
    let mut extra: Vec<&str> = Vec::new();
    let mut i = 0;
    while i < args.len() {
        let a = args[i].as_str();
        i += 1;
        if a == "-h" || a == "--help" {
            let _ = write!(out, "{}", spec.usage);
            return Err(0);
        }
        let (name, inline) = match a.split_once('=') {
            Some((n, v)) if n.starts_with("--") => (n, Some(v)),
            _ => (a, None),
        };
        let Some(&(known, kind)) = spec.options.iter().find(|(o, _)| *o == name) else {
            extra.push(a);
            continue;
        };
        let value = match inline {
            Some(v) => v.to_owned(),
            None => match args.get(i) {
                Some(v) if takes_value(v) => {
                    i += 1;
                    v.clone()
                }
                _ => return Err(refuse(format!("argument {name}: expected one argument"))),
            },
        };
        let slot = got.vals.entry(known).or_default();
        if kind == Kind::One {
            slot.clear();
        }
        slot.push(value);
    }
    let missing: Vec<&str> = spec
        .required
        .iter()
        .copied()
        .filter(|r| !got.vals.contains_key(r))
        .collect();
    if !missing.is_empty() {
        return Err(refuse(format!(
            "the following arguments are required: {}",
            missing.join(", ")
        )));
    }
    if !extra.is_empty() {
        return Err(refuse(format!(
            "unrecognized arguments: {}",
            extra.join(" ")
        )));
    }
    Ok(got)
}

#[cfg(test)]
mod tests {
    use super::*;

    const SPEC: Spec = Spec {
        prog: "og couple t",
        usage: "usage: og couple t\n",
        options: &[("--out", Kind::One), ("--in", Kind::Append)],
        required: &["--out"],
    };

    fn p(args: &[&str]) -> (Result<Args, u8>, String) {
        let v: Vec<String> = args.iter().map(|s| (*s).to_owned()).collect();
        let (mut out, mut err) = (Vec::new(), Vec::new());
        let r = parse(&SPEC, &v, &mut out, &mut err);
        (r, String::from_utf8(err).unwrap())
    }

    #[test]
    fn values_repeats_and_refusals() {
        let (a, _) = p(&["--in", "a", "--out=o", "--in=b", "--out", "-"]);
        let a = a.unwrap();
        assert_eq!(a.one("--out"), Some("-"));
        assert_eq!(a.many("--in"), ["a", "b"]);
        let (r, err) = p(&["--in", "a"]);
        assert_eq!(r.unwrap_err(), 2);
        assert!(err.ends_with("error: the following arguments are required: --out\n"));
        // An option the verb does not take is refused, before any work.
        let (r, err) = p(&["--out", "o", "--og-scaffold-probe"]);
        assert_eq!(r.unwrap_err(), 2);
        assert!(err.ends_with("error: unrecognized arguments: --og-scaffold-probe\n"));
        let (r, err) = p(&["--out"]);
        assert_eq!(r.unwrap_err(), 2);
        assert!(err.ends_with("argument --out: expected one argument\n"));
    }
}
