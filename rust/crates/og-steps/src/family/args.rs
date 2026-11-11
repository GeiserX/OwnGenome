//! The command lines of `og family`, read as argparse read the Python's:
//! `--name VALUE` or `--name=VALUE`, repeated options (`append`), options
//! that take every value up to the next option (`nargs="*"`), flags, and
//! fixed choices. No abbreviated options: no caller abbreviates one.

use std::collections::HashMap;

use super::Stop;

/// How an option takes its values.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Kind {
    /// One value; the last one given wins.
    One,
    /// One value each time, every one kept (`action="append"`).
    Append,
    /// Every value up to the next option (`nargs="*"`); the last use wins.
    Star,
    /// No value (`action="store_true"`).
    Flag,
    /// One of these values; the default when absent.
    Choice(&'static [&'static str], &'static str),
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
    /// The value of a `One` or `Choice` option.
    pub fn one(&self, name: &str) -> Option<&str> {
        self.vals
            .get(name)
            .and_then(|v| v.last())
            .map(String::as_str)
    }

    /// Every value of an `Append` or `Star` option (none when absent).
    pub fn many(&self, name: &str) -> &[String] {
        self.vals.get(name).map_or(&[], Vec::as_slice)
    }

    /// Whether a `Flag` was given.
    pub fn flag(&self, name: &str) -> bool {
        self.vals.contains_key(name)
    }
}

fn takes_value(a: &str) -> bool {
    a == "-"
        || !a.starts_with('-')
        || a.strip_prefix('-')
            .is_some_and(|d| d.parse::<f64>().is_ok())
}

/// Parses `args` by `spec`; a refusal prints argparse's message (exit 2),
/// `--help` the usage (exit 0).
pub fn parse(spec: &Spec, args: &[String]) -> Result<Args, Stop> {
    let usage_error = |msg: String| -> Stop {
        eprint!("{}", spec.usage);
        eprintln!("{}: error: {msg}", spec.prog);
        Stop::Exit(2)
    };
    let mut out = Args::default();
    let mut extra: Vec<String> = Vec::new();
    let mut i = 0;
    while i < args.len() {
        let a = &args[i];
        i += 1;
        if a == "-h" || a == "--help" {
            print!("{}", spec.usage);
            return Err(Stop::Exit(0));
        }
        let (name, inline) = match a.split_once('=') {
            Some((n, v)) if n.starts_with("--") => (n, Some(v.to_owned())),
            _ => (a.as_str(), None),
        };
        let Some(&(known, kind)) = spec.options.iter().find(|(o, _)| *o == name) else {
            extra.push(a.clone());
            continue;
        };
        let one = |inline: Option<String>, i: &mut usize| -> Result<String, Stop> {
            match inline {
                Some(v) => Ok(v),
                None => match args.get(*i) {
                    Some(v) if takes_value(v) => {
                        *i += 1;
                        Ok(v.clone())
                    }
                    _ => Err(usage_error(format!(
                        "argument {name}: expected one argument"
                    ))),
                },
            }
        };
        match kind {
            Kind::Flag => {
                if inline.is_some() {
                    return Err(usage_error(format!(
                        "argument {name}: ignored explicit argument '{}'",
                        inline.unwrap_or_default()
                    )));
                }
                out.vals.insert(known, Vec::new());
            }
            Kind::One => {
                let v = one(inline, &mut i)?;
                out.vals.insert(known, vec![v]);
            }
            Kind::Append => {
                let v = one(inline, &mut i)?;
                out.vals.entry(known).or_default().push(v);
            }
            Kind::Choice(choices, _) => {
                let v = one(inline, &mut i)?;
                if !choices.contains(&v.as_str()) {
                    let list: Vec<String> = choices.iter().map(|c| format!("'{c}'")).collect();
                    return Err(usage_error(format!(
                        "argument {name}: invalid choice: '{v}' (choose from {})",
                        list.join(", ")
                    )));
                }
                out.vals.insert(known, vec![v]);
            }
            Kind::Star => {
                let mut vs = Vec::new();
                if let Some(v) = inline {
                    vs.push(v);
                } else {
                    while let Some(v) = args.get(i).filter(|v| takes_value(v)) {
                        vs.push(v.clone());
                        i += 1;
                    }
                }
                out.vals.insert(known, vs);
            }
        }
    }
    let missing: Vec<&str> = spec
        .required
        .iter()
        .copied()
        .filter(|r| !out.vals.contains_key(r))
        .collect();
    if !missing.is_empty() {
        return Err(usage_error(format!(
            "the following arguments are required: {}",
            missing.join(", ")
        )));
    }
    if !extra.is_empty() {
        return Err(usage_error(format!(
            "unrecognized arguments: {}",
            extra.join(" ")
        )));
    }
    for (name, kind) in spec.options {
        if let Kind::Choice(_, default) = kind {
            out.vals
                .entry(name)
                .or_insert_with(|| vec![(*default).to_owned()]);
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    const SPEC: Spec = Spec {
        prog: "og family t",
        usage: "usage: og family t\n",
        options: &[
            ("--out", Kind::One),
            ("--gvcf", Kind::Append),
            ("--have", Kind::Star),
            ("--which", Kind::Flag),
            ("--panel", Kind::Choice(&["yes", "no"], "no")),
        ],
        required: &["--out"],
    };

    fn p(args: &[&str]) -> Result<Args, Stop> {
        let v: Vec<String> = args.iter().map(|s| (*s).to_owned()).collect();
        parse(&SPEC, &v)
    }

    #[test]
    fn kinds() {
        let a = p(&[
            "--have",
            "a",
            "b",
            "--gvcf",
            "x=1",
            "--out",
            "o",
            "--gvcf=y=2",
            "--which",
        ])
        .unwrap();
        assert_eq!(a.many("--have"), ["a", "b"]);
        assert_eq!(a.many("--gvcf"), ["x=1", "y=2"]);
        assert_eq!(a.one("--out"), Some("o"));
        assert!(a.flag("--which"));
        assert_eq!(a.one("--panel"), Some("no"));
        let a = p(&["--out", "o", "--have", "--panel", "yes"]).unwrap();
        assert!(a.many("--have").is_empty() && a.vals.contains_key("--have"));
        assert_eq!(a.one("--panel"), Some("yes"));
        assert!(matches!(p(&["--have", "a"]), Err(Stop::Exit(2))));
        assert!(matches!(
            p(&["--out", "o", "--panel", "maybe"]),
            Err(Stop::Exit(2))
        ));
        assert!(matches!(p(&["--out", "o", "stray"]), Err(Stop::Exit(2))));
        assert!(matches!(p(&["--out"]), Err(Stop::Exit(2))));
    }
}
