//! The pins: `versions.env`, and the few that live in `scripts/setup.sh`
//! and `scripts/lib/common.sh`.

use crate::{Lookup, failed};
use indexmap::IndexMap;
use regex::Regex;
use std::path::Path;
use std::sync::LazyLock;

/// `{NAME: value}` and `{NAME: note}` of a `versions.env`.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Pins {
    pub values: IndexMap<String, String>,
    pub notes: IndexMap<String, String>,
}

static PIN_LINE: LazyLock<Regex> = LazyLock::new(|| {
    Regex::new(r#"^([A-Z][A-Z0-9_]*)="([^"]*)"\s*(?:#\s*(.*))?$"#).expect("pin pattern")
});

/// The first sentence of a note, without its closing `.` or `;`.
pub fn first_sentence(text: &str) -> String {
    static S: LazyLock<Regex> =
        LazyLock::new(|| Regex::new(r"^(.+?[.;])(?:\s|$)").expect("sentence"));
    let first = S
        .captures(text)
        .map_or(text, |c| c.get(1).expect("group").as_str());
    first.trim_end_matches(['.', ';']).to_owned()
}

impl Pins {
    /// Reads `versions.env`. A pin's note is the trailing comment on its
    /// line, or a `# hold:` comment block directly above it.
    pub fn parse(text: &str) -> Pins {
        let mut p = Pins::default();
        let mut block: Vec<String> = vec![];
        for line in text.lines() {
            let s = line.trim();
            if let Some(c) = s.strip_prefix('#') {
                block.push(c.trim_start_matches('#').trim().to_owned());
                continue;
            }
            if let Some(c) = PIN_LINE.captures(s) {
                let name = c[1].to_owned();
                p.values.insert(name.clone(), c[2].to_owned());
                let mut note = c.get(3).map_or("", |m| m.as_str()).to_owned();
                if note.is_empty() && block.first().is_some_and(|b| b.starts_with("hold:")) {
                    note = block.join(" ");
                }
                if !note.is_empty() {
                    p.notes.insert(name, first_sentence(&note));
                }
            }
            block.clear();
        }
        p
    }

    pub fn read(path: &Path) -> std::io::Result<Pins> {
        Ok(Pins::parse(&std::fs::read_to_string(path)?))
    }

    /// A pin that must be set.
    pub fn need(&self, name: &str) -> Lookup<&str> {
        match self.values.get(name) {
            Some(v) if !v.is_empty() => Ok(v),
            _ => failed(format!("{name} is not set in versions.env")),
        }
    }

    pub fn get(&self, name: &str) -> Option<&str> {
        self.values.get(name).map(String::as_str)
    }
}

/// Only the `NAME="value"` lines of a `versions.env`, as a release's file
/// is read: no notes.
pub fn plain_pins(text: &str) -> IndexMap<String, String> {
    static P: LazyLock<Regex> =
        LazyLock::new(|| Regex::new(r#"(?m)^([A-Z][A-Z0-9_]*)="([^"]*)""#).expect("plain pin"));
    P.captures_iter(text)
        .map(|c| (c[1].to_owned(), c[2].to_owned()))
        .collect()
}

/// `HLA_DB_RELEASE`'s default in a `scripts/lib/common.sh` text
/// (`HLA_DB_RELEASE=${HLA_DB_RELEASE:-x.y.z}`).
pub fn hla_release_in(text: &str) -> Option<String> {
    static P: LazyLock<Regex> = LazyLock::new(|| {
        Regex::new(r"(?m)^HLA_DB_RELEASE=\$\{HLA_DB_RELEASE:-(\d+(?:\.\d+)+)\}").expect("hla")
    });
    P.captures(text).map(|c| c[1].to_owned())
}

/// [`hla_release_in`] of a file; a missing line is a failed lookup.
pub fn hla_db_release(common: &Path) -> Lookup<String> {
    let text = std::fs::read_to_string(common).unwrap_or_default();
    hla_release_in(&text).ok_or_else(|| {
        crate::LookupError::Failed(format!(
            "no HLA_DB_RELEASE=${{HLA_DB_RELEASE:-x.y.z}} line in {}",
            common.display()
        ))
    })
}

/// The `vep_cache_url` printf template of `scripts/lib/common.sh`.
pub fn vep_cache_template(common: &Path) -> Lookup<String> {
    static P: LazyLock<Regex> =
        LazyLock::new(|| Regex::new(r"vep_cache_url\(\)\s*\{\s*printf '([^']+)'").expect("vep"));
    let text = std::fs::read_to_string(common).unwrap_or_default();
    P.captures(&text).map(|c| c[1].to_owned()).ok_or_else(|| {
        crate::LookupError::Failed("no vep_cache_url template in scripts/lib/common.sh".into())
    })
}

/// The value of the first line matching `pattern` (one group) in a file.
pub fn setting(path: &Path, pattern: &Regex, what: &str) -> Lookup<String> {
    let text = std::fs::read_to_string(path).unwrap_or_default();
    pattern
        .captures(&text)
        .map(|c| c[1].to_owned())
        .ok_or_else(|| crate::LookupError::Failed(format!("no {what} in {}", path.display())))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn notes_come_from_the_line_or_a_hold_block() {
        let p = Pins::parse(
            "# hold: pypgx 0.27 broke. Keep it.\n# more\nPYPGX_IMAGE=\"a:1\"\n\n# not a hold\nX_IMAGE=\"b:2\"  # trailing note; and more\nY=\"\"\n",
        );
        assert_eq!(p.get("PYPGX_IMAGE"), Some("a:1"));
        assert_eq!(p.notes["PYPGX_IMAGE"], "hold: pypgx 0.27 broke");
        assert_eq!(p.notes["X_IMAGE"], "trailing note");
        assert!(p.need("Y").is_err() && p.need("Z").is_err());
        assert_eq!(first_sentence("no stop"), "no stop");
    }

    #[test]
    fn common_sh_values() {
        assert_eq!(
            hla_release_in("x\nHLA_DB_RELEASE=${HLA_DB_RELEASE:-3.65.0}\n").as_deref(),
            Some("3.65.0")
        );
        assert_eq!(hla_release_in("HLA_DB_RELEASE=3.65.0"), None);
        let pins = plain_pins("A=\"1\" # x\n  B=\"2\"\nC_D=\"\"\n");
        assert_eq!(pins.get("A").map(String::as_str), Some("1"));
        assert!(!pins.contains_key("B"));
        assert_eq!(pins.get("C_D").map(String::as_str), Some(""));
    }
}
