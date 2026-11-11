//! A task name as Nextflow writes it: `WORKFLOW:PROCESS (tag)`.

use pycompat::text::is_space;

/// `'SV:MANTA (s1)'` → `("MANTA", Some("s1"))`; `'FASTP'` → `("FASTP",
/// None)`. The process loses its workflow prefix (everything up to the last
/// `:`). The tag is what the outermost parentheses at the end hold, so
/// `'X (a) (b)'` has the tag `a) (b`. `None` when the name is not one word
/// with an optional parenthesised tag after it (`'X y'`, `''`).
///
/// This is the regex `^\s*(\S+?)\s*(?:\((.*)\))?\s*$` of the Python it
/// replaces, with its backtracking worked out: the process is the shortest
/// run of non-space characters after which the rest is spaces, or spaces, a
/// tag in parentheses and spaces.
pub fn task_name(name: &str) -> Option<(&str, Option<&str>)> {
    let start = name.find(|c: char| !is_space(c))?;
    let word = &name[start..];
    let run = word.find(is_space).unwrap_or(word.len());
    for (i, c) in word[..run].char_indices() {
        let end = i + c.len_utf8();
        if let Some(tag) = tail(&word[end..]) {
            let full = &word[..end];
            let process = full.rsplit_once(':').map_or(full, |(_, p)| p);
            return Some((process, tag));
        }
    }
    None
}

/// What may follow the process: `Some(None)` for spaces only,
/// `Some(Some(tag))` for spaces, `(tag)` and spaces, `None` otherwise. The
/// tag holds no line break (the regex's `.`).
fn tail(rest: &str) -> Option<Option<&str>> {
    let r = rest.trim_start_matches(is_space);
    if r.is_empty() {
        return Some(None);
    }
    let inner = r.strip_prefix('(')?;
    // Only the last ')' can have nothing but spaces after it.
    let close = inner.rfind(')')?;
    let tag = &inner[..close];
    (inner[close + 1..].chars().all(is_space) && !tag.contains('\n')).then_some(Some(tag))
}

#[cfg(test)]
mod tests {
    use super::task_name;

    #[test]
    fn names() {
        assert_eq!(task_name("SV:MANTA (s1)"), Some(("MANTA", Some("s1"))));
        assert_eq!(task_name("FASTP"), Some(("FASTP", None)));
        assert_eq!(task_name("  A:B:C  "), Some(("C", None)));
        assert_eq!(task_name("X (a) (b)"), Some(("X", Some("a) (b"))));
        assert_eq!(task_name("X(t)"), Some(("X", Some("t"))));
        assert_eq!(task_name("X ()"), Some(("X", Some(""))));
        assert_eq!(task_name("A(b)c"), Some(("A(b)c", None)));
        assert_eq!(task_name("X (a"), None);
        assert_eq!(task_name("X y"), None);
        assert_eq!(task_name(""), None);
        assert_eq!(task_name("   "), None);
        assert_eq!(task_name("X (a\nb)"), None);
        assert_eq!(task_name("P:Q (s1:part) \n"), Some(("Q", Some("s1:part"))));
    }
}
