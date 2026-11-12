//! argparse's reading of a step tool's command line: `--opt value` or
//! `--opt=value` (the last one wins), flags without a value, choices and
//! required options. Abbreviated options are not read (no caller uses them).

use std::collections::HashMap;

/// What an option takes.
pub(crate) enum Kind {
    /// A value.
    Value,
    /// A value from a fixed list.
    Choice(&'static [&'static str]),
    /// No value: present or not.
    Flag,
}

/// One option: its name, what it takes, and whether it is required.
pub(crate) struct Opt(pub &'static str, pub Kind, pub bool);

/// The options read: name to value (`""` for a flag that was given).
#[derive(Debug)]
pub(crate) struct Parsed(HashMap<&'static str, String>);

impl Parsed {
    /// The value of `name`, when given.
    pub(crate) fn get(&self, name: &str) -> Option<&str> {
        self.0.get(name).map(String::as_str)
    }

    /// The value of `name`, or `""`.
    pub(crate) fn text(&self, name: &str) -> &str {
        self.get(name).unwrap_or("")
    }

    /// True when the flag `name` was given.
    pub(crate) fn flag(&self, name: &str) -> bool {
        self.0.contains_key(name)
    }
}

/// True for an argument argparse reads as an option rather than a value:
/// it starts with `-` and is not a negative number.
fn looks_like_option(a: &str) -> bool {
    a.starts_with('-') && a.len() > 1 && pycompat::num::parse_float(a).is_err()
}

/// Reads `argv` against `opts`; the error is argparse's message.
pub(crate) fn parse(argv: &[String], opts: &[Opt]) -> Result<Parsed, String> {
    let mut got = HashMap::new();
    let mut it = argv.iter();
    while let Some(arg) = it.next() {
        let (key, inline) = match arg.split_once('=') {
            Some((k, v)) if k.starts_with("--") => (k, Some(v.to_owned())),
            _ => (arg.as_str(), None),
        };
        let Some(Opt(name, kind, _)) = opts.iter().find(|o| o.0 == key) else {
            return Err(format!("unrecognized arguments: {arg}"));
        };
        let value = match kind {
            Kind::Flag => {
                if inline.is_some() {
                    return Err(format!(
                        "argument {name}: ignored explicit argument {:?}",
                        inline.unwrap_or_default()
                    ));
                }
                String::new()
            }
            Kind::Value | Kind::Choice(_) => match inline {
                Some(v) => v,
                None => match it.next() {
                    Some(v) if !looks_like_option(v) => v.clone(),
                    _ => return Err(format!("argument {name}: expected one argument")),
                },
            },
        };
        if let Kind::Choice(choices) = kind
            && !choices.contains(&value.as_str())
        {
            let list: Vec<String> = choices.iter().map(|c| crate::repr(c)).collect();
            return Err(format!(
                "argument {name}: invalid choice: {} (choose from {})",
                crate::repr(&value),
                list.join(", ")
            ));
        }
        got.insert(*name, value);
    }
    let missing: Vec<&str> = opts
        .iter()
        .filter(|o| o.2 && !got.contains_key(o.0))
        .map(|o| o.0)
        .collect();
    if !missing.is_empty() {
        return Err(format!(
            "the following arguments are required: {}",
            missing.join(", ")
        ));
    }
    Ok(Parsed(got))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn v(a: &[&str]) -> Vec<String> {
        a.iter().map(|s| (*s).to_owned()).collect()
    }

    const OPTS: &[Opt] = &[
        Opt("--out", Kind::Value, true),
        Opt("--sex", Kind::Choice(&["", "male", "female"]), false),
        Opt("--zero", Kind::Flag, false),
        Opt("--warn", Kind::Value, false),
    ];

    #[test]
    fn reads_values_flags_and_choices() {
        let p = parse(
            &v(&[
                "--out", "a", "--sex=", "--zero", "--warn", "-0.5", "--out=b",
            ]),
            OPTS,
        )
        .unwrap();
        assert_eq!(p.get("--out"), Some("b"));
        assert_eq!(p.get("--sex"), Some(""));
        assert!(p.flag("--zero"));
        assert_eq!(p.text("--warn"), "-0.5");
    }

    #[test]
    fn refuses_what_argparse_refuses() {
        assert!(
            parse(&v(&[]), OPTS)
                .unwrap_err()
                .contains("required: --out")
        );
        assert!(
            parse(&v(&["--out"]), OPTS)
                .unwrap_err()
                .contains("expected one argument")
        );
        assert!(
            parse(&v(&["--out", "--zero"]), OPTS)
                .unwrap_err()
                .contains("expected one argument")
        );
        assert!(
            parse(&v(&["--out", "a", "--sex", "x"]), OPTS)
                .unwrap_err()
                .contains("invalid choice")
        );
        assert!(
            parse(&v(&["--out", "a", "--ou", "b"]), OPTS)
                .unwrap_err()
                .contains("unrecognized")
        );
        assert!(parse(&v(&["--out", "a", "--zero=1"]), OPTS).is_err());
    }
}
