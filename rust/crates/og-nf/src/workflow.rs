//! A workflow file (`workflows/*.nf`): its first workflow's name, its text
//! with comments blanked (strings respected), and the `emit:` lines that
//! name a channel (`name = ch_x`).

use crate::groovy::strip_comments;
use regex::Regex;
use std::sync::LazyLock;

static WORKFLOW: LazyLock<Regex> = LazyLock::new(|| re(r"(?m)^workflow\s+(\w+)"));
static EMIT: LazyLock<Regex> = LazyLock::new(|| re(r"(?m)^\s*(\w+)\s*=\s*(ch_\w+)\s*$"));

fn re(s: &str) -> Regex {
    Regex::new(s).expect("a valid built-in pattern")
}

/// A workflow file, read.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Workflow {
    pub name: String,
    /// The whole file with its comments blanked.
    pub text: String,
    /// `(emit, channel)` of each `emit = ch_x` line after the last `emit:`.
    pub emits: Vec<(String, String)>,
}

/// The workflow of `src`, or `None` when the file names no workflow or has
/// no `emit:` block.
pub fn read(src: &str) -> Option<Workflow> {
    let text = strip_comments(src);
    let name = WORKFLOW.captures(&text)?[1].to_owned();
    let at = text.rfind("emit:")?;
    let emits = EMIT
        .captures_iter(&text[at + "emit:".len()..])
        .map(|c| (c[1].to_owned(), c[2].to_owned()))
        .collect();
    Some(Workflow { name, text, emits })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reads_name_and_emits() {
        let src = "// workflow NOT\nworkflow PGX {\n  main:\n  ch_a = X.out.a // emit: no\n  emit:\n  a = ch_a\n  b = ch_b.mix(ch_a)\n  versions = ch_versions\n}\n";
        let w = read(src).unwrap();
        assert_eq!(w.name, "PGX");
        assert_eq!(
            w.emits,
            [
                ("a".into(), "ch_a".into()),
                ("versions".into(), "ch_versions".into())
            ]
        );
        assert!(read("workflow A {\n main:\n x = 1\n}\n").is_none());
        assert!(read("emit:\n a = ch_a\n").is_none());
    }
}
