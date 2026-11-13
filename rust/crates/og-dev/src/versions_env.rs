//! A strict reader of `versions.env`, the one list of images and coupled
//! data versions.
//!
//! The bash steps source the file, so a line could hold anything bash
//! accepts. The dev verbs never run it: they read it, and only in one shape,
//! `NAME="value"` with an optional trailing `# comment`, where the value
//! holds no `"`, `$` or backquote (nothing for bash to expand). A blank line
//! or a comment line is skipped. Any other line, and a name set twice, is an
//! error naming the line, so a check never reads a file it does not
//! understand as one with fewer pins.

use std::path::Path;
use std::sync::LazyLock;

use regex::Regex;

/// One `NAME="value"` line.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Pin {
    pub name: String,
    pub value: String,
    /// What follows the closing quote, as written (the comment, if any).
    pub rest: String,
    /// The 1-based line number.
    pub line: usize,
}

static LINE: LazyLock<Regex> = LazyLock::new(|| {
    Regex::new(r#"^([A-Z][A-Z0-9_]*)="([^"$`]*)"(\s*(?:#.*)?)$"#).expect("a valid pattern")
});

/// Reads `versions.env` text. Lines may end in `\n` or `\r\n`.
pub fn parse(text: &str) -> Result<Vec<Pin>, String> {
    let mut pins: Vec<Pin> = Vec::new();
    for (i, raw) in text.split('\n').enumerate() {
        let line = raw.strip_suffix('\r').unwrap_or(raw);
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let Some(c) = LINE.captures(line) else {
            return Err(format!(
                "versions.env:{}: not a NAME=\"value\" line (a value holds no \", $ or backquote): {line}",
                i + 1
            ));
        };
        let name = c[1].to_owned();
        if let Some(first) = pins.iter().find(|p| p.name == name) {
            return Err(format!(
                "versions.env:{}: {name} is set again (first on line {})",
                i + 1,
                first.line
            ));
        }
        pins.push(Pin {
            name,
            value: c[2].to_owned(),
            rest: c[3].to_owned(),
            line: i + 1,
        });
    }
    Ok(pins)
}

/// Reads the `versions.env` at `path`.
pub fn read(path: &Path) -> Result<Vec<Pin>, String> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| format!("cannot read {}: {e}", path.display()))?;
    parse(&text)
}

/// The value of `name`, if it is set.
pub fn get<'a>(pins: &'a [Pin], name: &str) -> Option<&'a Pin> {
    pins.iter().find(|p| p.name == name)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reads_values_and_comments() {
        let p = parse("# head\n\nA_IMAGE=\"x/y:1\"  # minimap2=2.31\r\nB=\"\"\n").unwrap();
        assert_eq!(p.len(), 2);
        assert_eq!(p[0].name, "A_IMAGE");
        assert_eq!(p[0].value, "x/y:1");
        assert_eq!(p[0].rest, "  # minimap2=2.31");
        assert_eq!(p[0].line, 3);
        assert_eq!(p[1].value, "");
    }

    #[test]
    fn refuses_what_bash_would_expand_or_a_second_value() {
        for bad in [
            "A=\"${B}\"",
            "A=\"$(date)\"",
            "A=\"`x`\"",
            "A=x",
            "export A=\"x\"",
            "A=\"x\" B",
            "a=\"x\"",
        ] {
            assert!(parse(bad).is_err(), "{bad} was read");
        }
        let e = parse("A=\"1\"\nA=\"2\"\n").unwrap_err();
        assert!(e.contains("A is set again (first on line 1)"), "{e}");
    }

    #[test]
    fn the_repository_file_reads() {
        let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../../versions.env");
        let pins = read(&path).unwrap();
        assert!(get(&pins, "SAMTOOLS_IMAGE").is_some());
    }
}
