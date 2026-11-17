//! Reading the inputs line by line, plain or gzip, the way the step's
//! Python read them: universal newlines, and a byte that is no UTF-8 read
//! as U+FFFD unless the file must be strict UTF-8 (the carrier list).

use og_io::{Lines, Utf8};
use og_pedigree::Die;

/// Calls `f` with every line of `path`, without its ending. One line buffer
/// is reused, so memory stays at the longest line.
pub fn each_line(
    path: &str,
    utf8: Utf8,
    mut f: impl FnMut(&str) -> Result<(), Die>,
) -> Result<(), Die> {
    let read = og_io::open(path).map_err(|e| Die::io(path, &e))?;
    let mut lines = Lines::new(read, utf8);
    let mut line = String::new();
    while lines.read_line(&mut line).map_err(|e| Die::io(path, &e))? {
        if line.ends_with('\n') {
            line.pop();
        }
        f(&line)?;
    }
    Ok(())
}

/// Calls `f` with the tab-separated fields of every data line (not starting
/// with `#`) of a VCF.
pub fn each_record(path: &str, mut f: impl FnMut(&[&str]) -> Result<(), Die>) -> Result<(), Die> {
    each_line(path, Utf8::Replace, |line| {
        if line.starts_with('#') {
            return Ok(());
        }
        let fields: Vec<&str> = line.split('\t').collect();
        f(&fields)
    })
}

/// A whole number as Python's `int()` reads a string of ASCII digits (the
/// callers check `isdigit()` first); a value past `i64` is refused.
pub fn int(s: &str, what: &str) -> Result<i64, Die> {
    pycompat::num::parse_int(s)
        .ok()
        .and_then(|v| i64::try_from(v).ok())
        .ok_or_else(|| Die::one(format!("{what}: '{s}' is not a whole number")))
}

/// `str.isdigit()` for ASCII text: not empty, digits only.
pub fn is_digits(s: &str) -> bool {
    !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit())
}

/// `info_map`: the `;`-separated `key=value` pairs of an INFO field, a
/// later key over an earlier one.
pub fn info_get<'a>(info: &'a str, key: &str) -> Option<&'a str> {
    let mut found = None;
    for kv in info.split(';') {
        let (k, v) = kv.split_once('=').unwrap_or((kv, ""));
        if k == key {
            found = Some(v);
        }
    }
    found
}
