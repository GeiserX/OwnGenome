//! The command lines of the ported tools, read the way argparse reads them:
//! `--name VALUE` or `--name=VALUE` (the last one given wins) and bare
//! `--flag`s. A refused command line is exit 2, as argparse does.

use std::collections::HashMap;

/// The options of one command line.
#[derive(Debug, Default)]
pub(crate) struct Args {
    values: HashMap<&'static str, String>,
    flags: Vec<&'static str>,
}

impl Args {
    /// Reads `argv` against the options a verb takes: `valued` names take a
    /// value, `flags` take none. An unknown word is an error.
    pub(crate) fn parse(
        argv: &[String],
        valued: &[&'static str],
        flags: &[&'static str],
    ) -> Result<Self, String> {
        let mut out = Args::default();
        let mut it = argv.iter();
        while let Some(word) = it.next() {
            let (name, inline) = match word.split_once('=') {
                Some((n, v)) if n.starts_with("--") => (n, Some(v.to_owned())),
                _ => (word.as_str(), None),
            };
            if let Some(&known) = valued.iter().find(|&&v| v == name) {
                let value = match inline {
                    Some(v) => v,
                    None => it
                        .next()
                        .filter(|v| !v.starts_with("--"))
                        .cloned()
                        .ok_or_else(|| format!("argument {known}: expected one argument"))?,
                };
                out.values.insert(known, value);
            } else if let Some(&known) = flags.iter().find(|&&f| f == name && inline.is_none()) {
                out.flags.push(known);
            } else {
                return Err(format!("unrecognized arguments: {word}"));
            }
        }
        Ok(out)
    }

    pub(crate) fn get(&self, name: &str) -> Option<&str> {
        self.values.get(name).map(String::as_str)
    }

    pub(crate) fn flag(&self, name: &str) -> bool {
        self.flags.contains(&name)
    }

    /// The value of an option the verb cannot run without.
    pub(crate) fn required(&self, name: &str) -> Result<&str, String> {
        self.get(name)
            .ok_or_else(|| format!("the following arguments are required: {name}"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn argv(s: &str) -> Vec<String> {
        s.split_whitespace().map(str::to_owned).collect()
    }

    #[test]
    fn values_flags_and_the_last_one_wins() {
        let a = Args::parse(
            &argv("--out a --flagged --out=b"),
            &["--out"],
            &["--flagged"],
        )
        .unwrap();
        assert_eq!(a.get("--out"), Some("b"));
        assert!(a.flag("--flagged"));
    }

    #[test]
    fn a_missing_value_or_an_unknown_word_is_refused() {
        assert!(Args::parse(&argv("--out"), &["--out"], &[]).is_err());
        assert!(Args::parse(&argv("--out --x"), &["--out"], &[]).is_err());
        assert!(Args::parse(&argv("stray"), &["--out"], &[]).is_err());
        assert!(Args::parse(&argv("--flagged=1"), &[], &["--flagged"]).is_err());
    }
}
