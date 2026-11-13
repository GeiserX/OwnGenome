//! `og y`: the Y chromosome tools.
//!
//! - `str-bed`, `str-report` ([`crate::ystr`], step 45): Y-STR counts.
//! - `tree-positions`, `tree-place` ([`crate::ytree`], step 37): the
//!   placement on the current YFull YTree.
//! - `private-snps` ([`crate::yprivate`], step 37): the chrY SNVs no public
//!   SNP table names.
//! - `yleaf-positions` ([`crate::yleaf`], step 37): Yleaf's own marker
//!   positions, for the pileup made beside it.

use std::io::Write;

const VERBS: &str = "str-bed|str-report|tree-positions|tree-place|private-snps|yleaf-positions";

/// Routes `og y VERB ARGS`; the exit code.
pub fn cli(verb: Option<&str>, rest: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    match verb {
        Some("str-bed" | "str-report") => crate::ystr::cli(verb, rest, out, err),
        Some("tree-positions") => crate::ytree::positions_cli(rest, out, err),
        Some("tree-place") => crate::ytree::place_cli(rest, out, err),
        Some("private-snps") => crate::yprivate::cli(rest, out, err),
        Some("yleaf-positions") => crate::yleaf::cli(rest, out, err),
        other => crate::unknown_verb("y", VERBS, other, err),
    }
}

/// `str(pathlib.Path(dir) / name)`: pathlib drops empty and `.` parts of
/// `dir` and keeps a leading `/` (two for a leading `//`).
pub(crate) fn path_join(dir: &str, name: &str) -> String {
    let mut s = String::new();
    if dir.starts_with("//") && !dir.starts_with("///") {
        s.push_str("//");
    } else if dir.starts_with('/') {
        s.push('/');
    }
    for part in dir.split('/').filter(|p| !p.is_empty() && *p != ".") {
        s.push_str(part);
        s.push('/');
    }
    s.push_str(name);
    s
}

/// True when `path` is a regular file, symlinks followed (`Path.is_file()`).
pub(crate) fn is_file(path: &str) -> bool {
    std::fs::metadata(path).is_ok_and(|m| m.is_file())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn path_join_as_pathlib() {
        assert_eq!(path_join(".", "t.json"), "t.json");
        assert_eq!(path_join("", "t.json"), "t.json");
        assert_eq!(path_join("a//b/./", "t.json"), "a/b/t.json");
        assert_eq!(path_join("/x/", "t.json"), "/x/t.json");
        assert_eq!(path_join("//x", "t.json"), "//x/t.json");
        assert_eq!(path_join("///x", "t.json"), "/x/t.json");
        assert_eq!(path_join("../d", "t.json"), "../d/t.json");
    }
}
