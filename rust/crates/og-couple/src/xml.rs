//! A small streaming XML reader for Orphadata's products: start tags with
//! their attributes, end tags and text, one at a time, so a 20 MB product
//! is never held whole. Comments, processing instructions and the DOCTYPE
//! are skipped; CDATA reads as text; the five named entities and numeric
//! character references are decoded.

use og_pedigree::Die;
use std::io::BufRead;

/// One piece of the document.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Event<'a> {
    /// A start tag: the name and the raw attribute text. A self-closing
    /// tag gives a `Start` and then an `End`.
    Start(&'a str, &'a str),
    End(&'a str),
    /// Text between tags, decoded.
    Text(String),
}

/// Decodes `&amp;`, `&lt;`, `&gt;`, `&quot;`, `&apos;`, `&#N;` and `&#xN;`.
pub fn decode(s: &str) -> String {
    if !s.contains('&') {
        return s.to_owned();
    }
    let mut out = String::with_capacity(s.len());
    let mut rest = s;
    while let Some(i) = rest.find('&') {
        out.push_str(&rest[..i]);
        rest = &rest[i..];
        let Some(j) = rest.find(';') else { break };
        let ent = &rest[1..j];
        let ch = match ent {
            "amp" => Some('&'),
            "lt" => Some('<'),
            "gt" => Some('>'),
            "quot" => Some('"'),
            "apos" => Some('\''),
            _ => ent
                .strip_prefix("#x")
                .map(|h| u32::from_str_radix(h, 16))
                .or_else(|| ent.strip_prefix('#').map(str::parse::<u32>))
                .and_then(Result::ok)
                .and_then(char::from_u32),
        };
        match ch {
            Some(c) => {
                out.push(c);
                rest = &rest[j + 1..];
            }
            None => {
                out.push('&');
                rest = &rest[1..];
            }
        }
    }
    out.push_str(rest);
    out
}

/// The value of attribute `name` in a start tag's attribute text.
pub fn attr(attrs: &str, name: &str) -> Option<String> {
    let mut rest = attrs;
    while let Some(eq) = rest.find('=') {
        let key = rest[..eq].trim();
        let after = rest[eq + 1..].trim_start();
        let q = after.chars().next()?;
        if q != '"' && q != '\'' {
            return None;
        }
        let end = after[1..].find(q)?;
        if key == name {
            return Some(decode(&after[1..=end]));
        }
        rest = &after[end + 2..];
    }
    None
}

fn bad(path: &str, what: &str) -> Die {
    Die::one(format!("{path}: not well-formed XML ({what})"))
}

/// Calls `f` with every event of the XML file at `path` (plain or gzip). A
/// file with no root element, or one that ends before its root closes (a
/// truncated download), is refused.
pub fn each_event(path: &str, mut f: impl FnMut(Event<'_>) -> Result<(), Die>) -> Result<(), Die> {
    let mut r = og_io::open(path).map_err(|e| Die::io(path, &e))?;
    // Open elements, and whether a root has started.
    let (mut depth, mut root) = (0usize, false);
    let mut buf: Vec<u8> = Vec::new();
    let io = |e: std::io::Error| Die::io(path, &e);
    loop {
        buf.clear();
        r.read_until(b'<', &mut buf).map_err(io)?;
        let open = buf.last() == Some(&b'<');
        if open {
            buf.pop();
        }
        if !buf.is_empty() {
            let text = std::str::from_utf8(&buf).map_err(|_| bad(path, "text is not UTF-8"))?;
            f(Event::Text(decode(text)))?;
        }
        if !open {
            return match (root, depth) {
                (false, _) => Err(bad(path, "no root element")),
                (true, 0) => Ok(()),
                _ => Err(bad(path, "the file ends before its root element closes")),
            };
        }
        buf.clear();
        r.read_until(b'>', &mut buf).map_err(io)?;
        // A comment or CDATA section may hold '>': read on to its own end.
        let ends = |b: &[u8], tail: &[u8]| b.len() > tail.len() && b[..b.len() - 1].ends_with(tail);
        if buf.starts_with(b"!--") {
            while !ends(&buf, b"--") {
                if r.read_until(b'>', &mut buf).map_err(io)? == 0 {
                    return Err(bad(path, "a comment does not end"));
                }
            }
            continue;
        }
        if buf.starts_with(b"![CDATA[") {
            while !ends(&buf, b"]]") {
                if r.read_until(b'>', &mut buf).map_err(io)? == 0 {
                    return Err(bad(path, "a CDATA section does not end"));
                }
            }
            let body = &buf[8..buf.len() - 3];
            let text = std::str::from_utf8(body).map_err(|_| bad(path, "text is not UTF-8"))?;
            f(Event::Text(text.to_owned()))?;
            continue;
        }
        if buf.last() != Some(&b'>') {
            return Err(bad(path, "a tag does not end"));
        }
        buf.pop();
        if buf.starts_with(b"?") || buf.starts_with(b"!") {
            continue;
        }
        let tag = std::str::from_utf8(&buf).map_err(|_| bad(path, "a tag is not UTF-8"))?;
        if let Some(name) = tag.strip_prefix('/') {
            depth = depth
                .checked_sub(1)
                .ok_or_else(|| bad(path, "an end tag with no start tag"))?;
            f(Event::End(name.trim()))?;
            continue;
        }
        let (tag, closed) = match tag.strip_suffix('/') {
            Some(t) => (t, true),
            None => (tag, false),
        };
        let split = tag
            .find(|c: char| c.is_ascii_whitespace())
            .unwrap_or(tag.len());
        let (name, attrs) = tag.split_at(split);
        if root && depth == 0 {
            return Err(bad(path, "a second root element"));
        }
        root = true;
        f(Event::Start(name, attrs))?;
        if closed {
            f(Event::End(name))?;
        } else {
            depth += 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn events_entities_and_attributes() {
        let dir = std::env::temp_dir().join(format!("og-couple-xml-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("t.xml");
        std::fs::write(
            &p,
            "<?xml version=\"1.0\"?><!-- a > b --><R date=\"2026-01-02 03:04:05\" x='y'>\
             <N lang=\"en\">A &amp; B&#39;s &#x3b1;</N><E/><![CDATA[<raw>]]></R>",
        )
        .unwrap();
        let mut got = Vec::new();
        each_event(&p.to_string_lossy(), |e| {
            got.push(match e {
                Event::Start(n, a) => format!("S {n}|{}", attr(a, "date").unwrap_or_default()),
                Event::End(n) => format!("E {n}"),
                Event::Text(t) => format!("T {t}"),
            });
            Ok(())
        })
        .unwrap();
        std::fs::remove_dir_all(&dir).unwrap();
        assert_eq!(
            got,
            [
                "S R|2026-01-02 03:04:05",
                "S N|",
                "T A & B's \u{3b1}",
                "E N",
                "S E|",
                "E E",
                "T <raw>",
                "E R"
            ]
        );
    }

    #[test]
    fn a_truncated_or_rootless_file_is_refused() {
        let dir = std::env::temp_dir().join(format!("og-couple-xml-cut-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let run = |name: &str, body: &str| {
            let p = dir.join(name);
            std::fs::write(&p, body).unwrap();
            each_event(&p.to_string_lossy(), |_| Ok(()))
                .err()
                .map(|e| e.lines.join(" "))
        };
        let cut = run("cut.xml", "<?xml version=\"1.0\"?><R><N>a</N><D>");
        assert!(cut.unwrap().contains("ends before its root element closes"));
        let none = run("none.xml", "<?xml version=\"1.0\"?>\n");
        assert!(none.unwrap().contains("no root element"));
        let two = run("two.xml", "<R/><S/>");
        assert!(two.unwrap().contains("a second root element"));
        let stray = run("stray.xml", "<R></R></S>");
        assert!(stray.unwrap().contains("an end tag with no start tag"));
        // Control: a whole document passes, with a comment after its root.
        assert_eq!(
            run("whole.xml", "<R><N>a</N><E/></R>\n<!-- end -->\n"),
            None
        );
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
