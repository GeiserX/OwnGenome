//! Escaping for the HTML report: Python's `html.escape`, with `quote=True`.
//!
//! Every text the report takes from summary.json goes through [`e`] before
//! it lands in the page, so a sample name or a note cannot open a tag.

use std::borrow::Cow;

/// The apostrophe as `html.escape` writes it. The `fault-escape-apos-39`
/// feature plants the other common spelling, which the e2e cases and the
/// tests would no longer find.
#[cfg(not(feature = "fault-escape-apos-39"))]
const APOS: &str = "&#x27;";
#[cfg(feature = "fault-escape-apos-39")]
const APOS: &str = "&#39;";

/// `html.escape(s)`: `&`, `<`, `>`, `"` and `'` as `&amp;`, `&lt;`, `&gt;`,
/// `&quot;` and `&#x27;`; everything else as it is.
pub fn e(s: &str) -> Cow<'_, str> {
    if !s.contains(['&', '<', '>', '"', '\'']) {
        return Cow::Borrowed(s);
    }
    let mut out = String::with_capacity(s.len() + 16);
    for c in s.chars() {
        match c {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            '\'' => out.push_str(APOS),
            c => out.push(c),
        }
    }
    Cow::Owned(out)
}

#[cfg(test)]
mod tests {
    use super::e;

    #[test]
    fn escapes_the_five_characters_and_nothing_else() {
        assert_eq!(e("plain ÄÖ 😀"), "plain ÄÖ 😀");
        #[cfg(not(feature = "fault-escape-apos-39"))]
        assert_eq!(e(r#"a&b<c>d"e'f"#), "a&amp;b&lt;c&gt;d&quot;e&#x27;f");
    }
}
