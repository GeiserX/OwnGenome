//! Reading files as Python's `open(path)` does, with the exception text
//! Python would put in the report (`FileNotFoundError: [Errno 2] No such
//! file or directory: 'x.tsv'`), because the recommendations quote it.

use pycompat::json::{self, Value};
use std::fmt;
use std::io;

/// A Python exception as `f"{type(e).__name__}: {e}"` shows it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PyErr {
    /// The exception class: `ValueError`, `FileNotFoundError`, `Error`
    /// (`csv.Error`), ...
    pub class: String,
    /// `str(e)`.
    pub msg: String,
}

impl PyErr {
    pub fn new(class: &str, msg: impl Into<String>) -> Self {
        PyErr {
            class: class.to_owned(),
            msg: msg.into(),
        }
    }

    /// `ValueError(msg)`.
    pub fn value(msg: impl Into<String>) -> Self {
        PyErr::new("ValueError", msg)
    }
}

impl fmt::Display for PyErr {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}: {}", self.class, self.msg)
    }
}

/// `repr(s)` of a string.
pub fn repr(s: &str) -> String {
    pycompat::value::py_repr(&Value::Str(s.to_owned())).unwrap_or_else(|_| format!("'{s}'"))
}

/// An `OSError` from `open(path)`: the subclass Python raises for the errno
/// and `[Errno N] strerror: 'path'`.
pub fn os_error(e: &io::Error, path: &str) -> PyErr {
    let Some(errno) = e.raw_os_error() else {
        return PyErr::new("OSError", e.to_string());
    };
    let text = e.to_string();
    let strerror = text
        .strip_suffix(&format!(" (os error {errno})"))
        .unwrap_or(&text);
    let class = match errno {
        2 => "FileNotFoundError",
        1 | 13 => "PermissionError",
        20 => "NotADirectoryError",
        21 => "IsADirectoryError",
        _ => "OSError",
    };
    PyErr::new(class, format!("[Errno {errno}] {strerror}: {}", repr(path)))
}

/// The `UnicodeDecodeError` text of Python's UTF-8 codec for `bytes`.
fn decode_error(bytes: &[u8], e: &std::str::Utf8Error) -> PyErr {
    let at = e.valid_up_to();
    let msg = match e.error_len() {
        None if bytes.len() - at > 1 => format!(
            "'utf-8' codec can't decode bytes in position {at}-{}: unexpected end of data",
            bytes.len() - 1
        ),
        None => format!(
            "'utf-8' codec can't decode byte 0x{:02x} in position {at}: unexpected end of data",
            bytes[at]
        ),
        Some(_) => {
            let b = bytes[at];
            let why = if (0x80..0xc2).contains(&b) || b > 0xf4 {
                "invalid start byte"
            } else {
                "invalid continuation byte"
            };
            format!("'utf-8' codec can't decode byte 0x{b:02x} in position {at}: {why}")
        }
    };
    PyErr::new("UnicodeDecodeError", msg)
}

/// The whole of `path` as `open(path).read()` gives it: strict UTF-8, and
/// `\r\n` and a lone `\r` read as `\n` (universal newlines).
pub fn read_text(path: &str) -> Result<String, PyErr> {
    let bytes = std::fs::read(path).map_err(|e| os_error(&e, path))?;
    let text = std::str::from_utf8(&bytes).map_err(|e| decode_error(&bytes, &e))?;
    Ok(if text.contains('\r') {
        text.replace("\r\n", "\n").replace('\r', "\n")
    } else {
        text.to_owned()
    })
}

/// The lines of `text`, each with its `\n`, as iterating a Python file
/// gives them.
pub fn lines(text: &str) -> impl Iterator<Item = io::Result<&str>> {
    text.split_inclusive('\n').map(Ok)
}

/// `json.load(open(path))`, with Python's message on a failure.
pub fn load_json(path: &str) -> Result<Value, PyErr> {
    let text = read_text(path)?;
    if text.starts_with('\u{feff}') {
        return Err(PyErr::new(
            "JSONDecodeError",
            "Unexpected UTF-8 BOM (decode using utf-8-sig): line 1 column 1 (char 0)",
        ));
    }
    json::loads(&text).map_err(|e| match e {
        json::JsonError::Decode(what, at) => {
            // Python counts code points, and the column from the last newline.
            let before = &text[..at.min(text.len())];
            let pos = before.chars().count();
            let line = before.matches('\n').count() + 1;
            let col = match before.rfind('\n') {
                Some(i) => before[i + 1..].chars().count() + 1,
                None => pos + 1,
            };
            PyErr::new(
                "JSONDecodeError",
                format!("{what}: line {line} column {col} (char {pos})"),
            )
        }
        other => PyErr::new(other.exception(), other.to_string()),
    })
}

/// A `csv.Error` as Python words it (no line number).
pub fn csv_error(e: &pycompat::csv::Error) -> PyErr {
    match e {
        pycompat::csv::Error::NewlineInUnquotedField { .. } => PyErr::new(
            "Error",
            "new-line character seen in unquoted field - do you need to open the file with newline=''?",
        ),
        pycompat::csv::Error::FieldTooLarge { .. } => PyErr::new(
            "Error",
            format!(
                "field larger than field limit ({})",
                pycompat::csv::FIELD_SIZE_LIMIT
            ),
        ),
        pycompat::csv::Error::Io(e) => PyErr::new("OSError", e.to_string()),
    }
}

/// `os.path.basename(path)`.
pub fn basename(path: &str) -> &str {
    path.rsplit('/').next().unwrap_or(path)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn os_errors_read_as_python_words_them() {
        let e = io::Error::from_raw_os_error(2);
        let got = os_error(&e, "m.txt").to_string();
        assert_eq!(
            got,
            "FileNotFoundError: [Errno 2] No such file or directory: 'm.txt'"
        );
    }

    #[test]
    fn decode_errors_name_the_byte() {
        let b: Vec<u8> = [&b"Gene\n"[..], &[0xff, 0xfe]].concat();
        let e = std::str::from_utf8(&b).unwrap_err();
        assert_eq!(
            decode_error(&b, &e).to_string(),
            "UnicodeDecodeError: 'utf-8' codec can't decode byte 0xff in position 5: invalid start byte"
        );
    }

    #[test]
    fn json_errors_give_line_and_column() {
        let dir = std::env::temp_dir().join(format!("og-pgx-json-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("r.json");
        std::fs::write(&p, "not json\n").unwrap();
        let e = load_json(p.to_str().unwrap()).unwrap_err();
        assert_eq!(
            e.to_string(),
            "JSONDecodeError: Expecting value: line 1 column 1 (char 0)"
        );
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn basename_is_the_last_part() {
        assert_eq!(basename("a/b/report.json"), "report.json");
        assert_eq!(basename("report.json"), "report.json");
        assert_eq!(basename("dir/"), "");
    }
}
