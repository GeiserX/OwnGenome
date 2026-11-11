//! A small argparse: `--name VALUE` or `--name=VALUE` anywhere, positionals
//! in order, `--` ends the options, `-h`/`--help` prints the usage. No
//! abbreviated options: no caller abbreviates one.

use std::collections::HashMap;

/// What a tool takes.
pub struct Spec {
    /// `og sv inv-filter`, for messages.
    pub prog: &'static str,
    /// The usage text, printed for `--help` and above an error.
    pub usage: &'static str,
    /// Every option, each taking one value.
    pub options: &'static [&'static str],
    /// The positional arguments, by name.
    pub positionals: &'static [&'static str],
    /// The options that must be given.
    pub required: &'static [&'static str],
}

/// A parsed command line.
#[derive(Debug, Default)]
pub struct Parsed {
    opts: HashMap<&'static str, String>,
    /// The positionals, as many as the spec names.
    pub pos: Vec<String>,
}

/// Why parsing stopped.
#[derive(Debug, PartialEq, Eq)]
pub enum ArgError {
    /// `-h` or `--help`.
    Help,
    /// A command line argparse would refuse, and why.
    Usage(String),
}

impl Parsed {
    /// The value of `--name`, when given (the last one wins, as in argparse).
    pub fn get(&self, name: &str) -> Option<&str> {
        self.opts.get(name).map(String::as_str)
    }

    /// `--name` read by `parse`, or `default` when it is absent.
    pub fn typed<T>(
        &self,
        name: &str,
        default: T,
        kind: &str,
        parse: impl Fn(&str) -> Option<T>,
    ) -> Result<T, ArgError> {
        match self.get(name) {
            None => Ok(default),
            Some(v) => parse(v).ok_or_else(|| {
                ArgError::Usage(format!(
                    "argument {name}: invalid {kind} value: {}",
                    crate::py_repr(v)
                ))
            }),
        }
    }
}

/// Parses `args` (the words after the tool's own) by `spec`.
pub fn parse(spec: &Spec, args: &[String]) -> Result<Parsed, ArgError> {
    let mut out = Parsed::default();
    let mut extra = Vec::new();
    let mut i = 0;
    let mut options_done = false;
    while i < args.len() {
        let a = &args[i];
        i += 1;
        if options_done || a == "-" || !a.starts_with('-') {
            if out.pos.len() < spec.positionals.len() {
                out.pos.push(a.clone());
            } else {
                extra.push(a.clone());
            }
            continue;
        }
        if a == "--" {
            options_done = true;
            continue;
        }
        if a == "-h" || a == "--help" {
            return Err(ArgError::Help);
        }
        let (name, inline) = match a.split_once('=') {
            Some((n, v)) if n.starts_with("--") => (n, Some(v.to_owned())),
            _ => (a.as_str(), None),
        };
        let Some(&known) = spec.options.iter().find(|&&o| o == name) else {
            extra.push(a.clone());
            continue;
        };
        let value = match inline {
            Some(v) => v,
            None => match args.get(i) {
                Some(v) if v == "-" || !v.starts_with('-') || is_negative_number(v) => {
                    i += 1;
                    v.clone()
                }
                _ => {
                    return Err(ArgError::Usage(format!(
                        "argument {name}: expected one argument"
                    )));
                }
            },
        };
        out.opts.insert(known, value);
    }
    let missing: Vec<&str> = spec.positionals[out.pos.len()..]
        .iter()
        .copied()
        .chain(
            spec.required
                .iter()
                .copied()
                .filter(|r| !out.opts.contains_key(r)),
        )
        .collect();
    if !missing.is_empty() {
        return Err(ArgError::Usage(format!(
            "the following arguments are required: {}",
            missing.join(", ")
        )));
    }
    if !extra.is_empty() {
        return Err(ArgError::Usage(format!(
            "unrecognized arguments: {}",
            extra.join(" ")
        )));
    }
    Ok(out)
}

/// argparse takes `-5` as a value when no option looks like a number.
fn is_negative_number(s: &str) -> bool {
    s.strip_prefix('-')
        .is_some_and(|d| !d.is_empty() && d.parse::<f64>().is_ok())
}

/// Prints the outcome of a failed parse as argparse does and returns the
/// exit code: 0 for `--help`, 2 for a refused command line.
pub fn report(spec: &Spec, e: &ArgError) -> u8 {
    match e {
        ArgError::Help => {
            print!("{}", spec.usage);
            0
        }
        ArgError::Usage(msg) => {
            eprint!("{}", spec.usage);
            eprintln!("{}: error: {msg}", spec.prog);
            2
        }
    }
}

/// `int(s)` as argparse's `type=int` reads it.
pub fn int(s: &str) -> Option<i128> {
    pycompat::num::parse_int(s).ok()
}

/// `float(s)` as argparse's `type=float` reads it.
pub fn float(s: &str) -> Option<f64> {
    pycompat::num::parse_float(s).ok()
}

#[cfg(test)]
mod tests {
    use super::*;

    const SPEC: Spec = Spec {
        prog: "og t",
        usage: "usage: og t [--n N] A B\n",
        options: &["--n", "--out"],
        positionals: &["a", "b"],
        required: &[],
    };

    fn p(args: &[&str]) -> Result<Parsed, ArgError> {
        parse(
            &SPEC,
            &args.iter().map(|s| (*s).to_owned()).collect::<Vec<_>>(),
        )
    }

    #[test]
    fn options_anywhere_and_both_spellings() {
        let r = p(&["x", "--n", "3", "y", "--out=o", "--n=4"]).unwrap();
        assert_eq!(r.pos, ["x", "y"]);
        assert_eq!(r.get("--n"), Some("4"));
        assert_eq!(r.get("--out"), Some("o"));
        let r = p(&["--n", "-5", "-", "--", "-z"]).unwrap();
        assert_eq!(r.get("--n"), Some("-5"));
        assert_eq!(r.pos, ["-", "-z"]);
    }

    #[test]
    fn refusals() {
        assert_eq!(
            p(&["x"]).unwrap_err(),
            ArgError::Usage("the following arguments are required: b".into())
        );
        assert_eq!(
            p(&["x", "y", "z", "--nn", "1"]).unwrap_err(),
            ArgError::Usage("unrecognized arguments: z --nn 1".into())
        );
        assert_eq!(
            p(&["x", "y", "--n"]).unwrap_err(),
            ArgError::Usage("argument --n: expected one argument".into())
        );
        assert_eq!(p(&["--help"]).unwrap_err(), ArgError::Help);
        let r = p(&["x", "y", "--n", "q"]).unwrap();
        assert_eq!(
            r.typed("--n", 0, "int", int).unwrap_err(),
            ArgError::Usage("argument --n: invalid int value: 'q'".into())
        );
    }
}
