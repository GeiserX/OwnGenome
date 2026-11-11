//! The four filters that sit between a third-party tool and the next one:
//!
//! - [`svlen`] (`og sv svlen-fix`): Manta's sequence-resolved records get
//!   their net length change before AnnotSV reads them (step 05).
//! - [`inv`] (`og sv inv-filter`): an inversion Manta saw one junction of
//!   is not a whole-span inversion (step 04).
//! - [`refblock`] (`og pgx refblock-fill`): PharmCAT's multi-base positions
//!   called 0/0 from the gVCF's reference blocks (step 07).
//! - [`mito`] (`og pgx mito-merge`): the MT-RNR1 outside call PharmCAT
//!   gets from the chrM calls (step 07).
//!
//! Each one writes what the Python it replaced wrote, byte for byte
//! (`tests/parity`), except where it fixes that Python on purpose; each
//! module says where. Every input error is a message and the tool's
//! documented exit code, never a crash.

pub mod args;
pub mod inv;
pub mod mito;
pub mod refblock;
pub mod svlen;
mod zip;

use std::io::{self, BufRead};

/// Text lines as Python reads them: a file with universal newlines (`\n`,
/// `\r\n` and a lone `\r` each end a line and read as `\n`), and
/// `sys.stdin`, which Python opens with `newline="\n"` on Unix, so only
/// `\n` ends a line and a `\r` stays in it. UTF-8, strict.
pub enum TextLines<R> {
    Universal(og_io::Lines<R>),
    Stdin(R),
}

impl<R: BufRead> TextLines<R> {
    /// A file's lines (universal newlines).
    pub fn file(r: R) -> Self {
        TextLines::Universal(og_io::Lines::new(r, og_io::Utf8::Strict))
    }

    /// `sys.stdin`'s lines (only `\n` ends one).
    pub fn stdin(r: R) -> Self {
        TextLines::Stdin(r)
    }

    /// The next line with its `\n`, if it has one; `false` at the end.
    pub fn read_line(&mut self, out: &mut String) -> io::Result<bool> {
        match self {
            TextLines::Universal(l) => l.read_line(out),
            TextLines::Stdin(r) => {
                out.clear();
                let mut raw = Vec::new();
                if r.read_until(b'\n', &mut raw)? == 0 {
                    return Ok(false);
                }
                out.push_str(
                    std::str::from_utf8(&raw)
                        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?,
                );
                Ok(true)
            }
        }
    }
}

/// Python's `repr()` of a str, for messages that quote a path or a value
/// as Python did: single quotes unless the text holds one and no double.
pub fn py_repr(s: &str) -> String {
    let quote = if s.contains('\'') && !s.contains('"') {
        '"'
    } else {
        '\''
    };
    let mut out = String::with_capacity(s.len() + 2);
    out.push(quote);
    for c in s.chars() {
        match c {
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if c == quote => {
                out.push('\\');
                out.push(c);
            }
            c if (c as u32) < 0x20 || c as u32 == 0x7f => {
                out.push_str(&format!("\\x{:02x}", c as u32));
            }
            c => out.push(c),
        }
    }
    out.push(quote);
    out
}

/// An `OSError` as Python prints it: `[Errno 2] No such file or directory:
/// 'x.vcf'`. Errors without an errno print as Rust has them.
pub fn os_error(e: &io::Error, path: &str) -> String {
    match e.raw_os_error() {
        Some(n) => {
            let text = io::Error::from_raw_os_error(n).to_string();
            let text = text
                .strip_suffix(&format!(" (os error {n})"))
                .unwrap_or(&text)
                .to_owned();
            format!("[Errno {n}] {text}: {}", py_repr(path))
        }
        None => e.to_string(),
    }
}

/// Python's `len()` of a str: code points, not bytes.
pub(crate) fn chars(s: &str) -> usize {
    s.chars().count()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn repr_quotes_as_python() {
        assert_eq!(py_repr("a.vcf"), "'a.vcf'");
        assert_eq!(py_repr("it's"), "\"it's\"");
        assert_eq!(py_repr("a'\"b"), "'a\\'\"b'");
        assert_eq!(py_repr("a\tb\\"), "'a\\tb\\\\'");
        assert_eq!(py_repr("\u{1}"), "'\\x01'");
    }

    #[test]
    fn os_errors_read_as_python_prints_them() {
        let e = std::fs::File::open("/no/such/file.vcf").unwrap_err();
        assert_eq!(
            os_error(&e, "/no/such/file.vcf"),
            "[Errno 2] No such file or directory: '/no/such/file.vcf'"
        );
    }
}
