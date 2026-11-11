//! The little of Groovy that the readers of `.nf` and `.config` files need:
//! where its comments and strings are, so a `//` inside a URL is not a
//! comment and a `}` inside a string never closes a block.
//!
//! Strings are `'…'`, `"…"`, `'''…'''` and `"""…"""`, each with backslash
//! escapes. Slashy strings (`/…/`) are not strings here: a `/` is an
//! operator, which is what every `.nf` file of this repository needs. A
//! one-line string left open ends at its line, so one stray quote cannot
//! swallow the rest of a file.

/// What a byte of the source is part of.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Part {
    Code,
    Comment,
    /// A string, its quotes included.
    Str,
}

/// The [`Part`] of every byte of `src`.
pub fn parts(src: &str) -> Vec<Part> {
    #[derive(Clone, Copy)]
    enum State {
        Code,
        Line,
        Block,
        Str(u8, bool),
    }
    let b = src.as_bytes();
    let mut out = vec![Part::Code; b.len()];
    let mut state = State::Code;
    let mut i = 0;
    while i < b.len() {
        let c = b[i];
        let at = |k: usize| b.get(i + k).copied();
        match state {
            State::Code => match c {
                b'/' if at(1) == Some(b'/') => {
                    state = State::Line;
                    out[i] = Part::Comment;
                }
                b'/' if at(1) == Some(b'*') => {
                    state = State::Block;
                    out[i] = Part::Comment;
                    out[i + 1] = Part::Comment;
                    i += 1;
                }
                b'\'' | b'"' => {
                    let triple = at(1) == Some(c) && at(2) == Some(c);
                    let n = if triple { 3 } else { 1 };
                    out[i..i + n].fill(Part::Str);
                    state = State::Str(c, triple);
                    i += n - 1;
                }
                _ => {}
            },
            State::Line => {
                if c == b'\n' {
                    state = State::Code;
                } else {
                    out[i] = Part::Comment;
                }
            }
            State::Block => {
                out[i] = Part::Comment;
                if c == b'*' && at(1) == Some(b'/') {
                    out[i + 1] = Part::Comment;
                    i += 1;
                    state = State::Code;
                }
            }
            State::Str(q, triple) => {
                out[i] = Part::Str;
                if c == b'\\' && i + 1 < b.len() {
                    out[i + 1] = Part::Str;
                    i += 1;
                } else if c == q && !triple {
                    state = State::Code;
                } else if c == q && at(1) == Some(q) && at(2) == Some(q) {
                    out[i + 1] = Part::Str;
                    out[i + 2] = Part::Str;
                    i += 2;
                    state = State::Code;
                } else if c == b'\n' && !triple {
                    state = State::Code;
                }
            }
        }
        i += 1;
    }
    out
}

/// `src` with each `//` and `/* */` comment blanked: every byte of it a
/// space, its newlines kept, so byte offsets and line numbers stay those of
/// `src`. Strings are left as they are.
pub fn strip_comments(src: &str) -> String {
    let parts = parts(src);
    let out: Vec<u8> = src
        .bytes()
        .zip(parts)
        .map(|(c, p)| {
            if p == Part::Comment && c != b'\n' {
                b' '
            } else {
                c
            }
        })
        .collect();
    // A comment starts at an ASCII `/` and ends at a newline or an ASCII
    // `/`, so only whole characters were blanked.
    String::from_utf8(out).expect("blanking whole characters keeps UTF-8")
}

/// The index of the `}` that closes the `{` at `open`, counting only braces
/// outside strings and comments; `None` when the block never closes or
/// `open` is not a `{`.
pub fn closing_brace(src: &str, open: usize) -> Option<usize> {
    let b = src.as_bytes();
    if b.get(open) != Some(&b'{') {
        return None;
    }
    let parts = parts(&src[open..]);
    let mut depth = 0usize;
    for (i, (&c, p)) in b[open..].iter().zip(parts).enumerate() {
        if p != Part::Code {
            continue;
        }
        match c {
            b'{' => depth += 1,
            b'}' => {
                depth -= 1;
                if depth == 0 {
                    return Some(open + i);
                }
            }
            _ => {}
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn comments_go_strings_stay() {
        let src = "a = 'http://x' // gone\nb = \"//kept\" /* gone\nalso */ c\n";
        let out = strip_comments(src);
        assert_eq!(out.len(), src.len());
        assert!(out.contains("'http://x'"));
        assert!(out.contains("\"//kept\""));
        assert!(!out.contains("gone"));
        assert!(!out.contains("also"));
        assert!(out.ends_with(" c\n"));
        assert_eq!(out.lines().count(), src.lines().count());
    }

    #[test]
    fn triple_quoted_strings_hold_quotes_and_slashes() {
        let src = "s = \"\"\"it's // not \"a\" comment\"\"\" // but this is\n";
        let out = strip_comments(src);
        assert!(out.contains("it's // not \"a\" comment"));
        assert!(!out.contains("but this is"));
    }

    #[test]
    fn a_brace_inside_a_string_or_comment_does_not_close() {
        let src = "x { a = '}' // }\n b = { y -> \"}\" } }";
        assert_eq!(closing_brace(src, 2), Some(src.len() - 1));
        assert_eq!(closing_brace(src, 0), None);
        assert_eq!(closing_brace("x { never", 2), None);
    }

    #[test]
    fn an_escaped_quote_does_not_end_a_string() {
        let src = "{a = 'it\\'s } here' }";
        assert_eq!(closing_brace(src, 0), Some(src.len() - 1));
    }

    #[test]
    fn an_unclosed_quote_ends_at_its_line() {
        let src = "{ a = 'oops\n b = 1 }";
        assert_eq!(closing_brace(src, 0), Some(src.len() - 1));
    }
}
