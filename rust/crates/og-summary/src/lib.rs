//! The summary cluster of OwnGenome: the step tools that turn a third-party
//! tool's output into one of our stable tables.
//!
//! - [`qc`]: `og qc sample`, step 33's identity and contamination verdict
//!   from somalier and VerifyBamID2 (`qc/<id>_sample_qc.tsv`).
//! - [`prs`]: `og prs table`, step 25's PRS table and step 26's ancestry
//!   table from pgsc_calc's output (`prs/<id>_prs_summary.tsv`,
//!   `ancestry/<id>_ancestry.tsv`).
//! - [`mito`]: `og mito annotate|deletions|copies`, step 20's chrM tables.
//! - [`collect`]: `og summary collect`, summary.json, ported one group of
//!   sections at a time beside the Python collector (`--only a`).
//!
//! Each tool exits 0 when it wrote its table, 1 with `ERROR: <message>`
//! when an input cannot be read or the table cannot be written, and 2 with
//! the usage for a command line it cannot read.

pub mod collect;
pub mod mito;
pub mod prs;
pub mod qc;

mod args;

use pycompat::json::Value;
use pycompat::num::{self, ParseError};
use pycompat::value::py_repr;
use std::io::{self, Write};
use std::path::Path;

/// `repr(s)` of a string, as Python's error messages quote it.
pub(crate) fn repr(s: &str) -> String {
    py_repr(&Value::Str(s.to_owned())).unwrap_or_else(|_| format!("'{s}'"))
}

/// `int(s)`, with Python's message when it is not an integer.
pub(crate) fn int(s: &str) -> Result<i128, String> {
    num::parse_int(s).map_err(|e| match e {
        ParseError::Invalid => format!("invalid literal for int() with base 10: {}", repr(s)),
        ParseError::TooLarge => format!("{} does not fit 128 bits", repr(s.trim())),
    })
}

/// `float(s)`, with Python's message when it is not a number.
pub(crate) fn float(s: &str) -> Result<f64, String> {
    num::parse_float(s).map_err(|_| format!("could not convert string to float: {}", repr(s)))
}

/// `int(x)` of a float: towards zero. Python's int has no bound; a value
/// that is not finite, or past 128 bits, is an error here as it is there
/// (OverflowError, ValueError), with a message instead of a traceback.
pub(crate) fn trunc(x: f64) -> Result<i128, String> {
    if x.is_nan() {
        return Err("cannot convert float NaN to integer".into());
    }
    if x.is_infinite() {
        return Err("cannot convert float infinity to integer".into());
    }
    let t = x.trunc();
    // 2^127: the first float past i128.
    if t.abs() >= 170_141_183_460_469_231_731_687_303_715_884_105_728.0 {
        return Err(format!("{} does not fit 128 bits", pycompat::fmt::repr(x)));
    }
    Ok(t as i128)
}

/// `os.path.basename(p)`: the text after the last `/`.
pub(crate) fn basename(p: &str) -> &str {
    p.rsplit('/').next().unwrap_or(p)
}

/// `os.path.join(a, b)`: one `/` between them unless `a` is empty or
/// already ends with one.
pub(crate) fn join(a: &str, b: &str) -> String {
    if a.is_empty() || a.ends_with('/') {
        format!("{a}{b}")
    } else {
        format!("{a}/{b}")
    }
}

/// `os.path.isfile(p)`.
pub(crate) fn is_file(p: impl AsRef<Path>) -> bool {
    p.as_ref().is_file()
}

/// The message of an I/O error as Python prints an `OSError` naming a file.
pub(crate) fn os_error(e: &io::Error, path: &str) -> String {
    let what = match e.kind() {
        io::ErrorKind::NotFound => "No such file or directory".to_owned(),
        io::ErrorKind::PermissionDenied => "Permission denied".to_owned(),
        io::ErrorKind::IsADirectory => "Is a directory".to_owned(),
        _ => {
            let s = e.to_string();
            return match s.split_once(" (os error ") {
                Some((msg, code)) => {
                    format!(
                        "[Errno {}] {msg}: {}",
                        code.trim_end_matches(')'),
                        repr(path)
                    )
                }
                None => format!("{s}: {}", repr(path)),
            };
        }
    };
    format!(
        "[Errno {}] {what}: {}",
        e.raw_os_error().unwrap_or(0),
        repr(path)
    )
}

/// Opens `path` as Python's `open(path, errors="replace")` or
/// `gzip.open(path, "rt", errors="replace")` reads it: universal newlines,
/// invalid UTF-8 as U+FFFD, gzip told by its magic bytes.
pub(crate) fn text_lines(path: &str) -> Result<og_io::Lines<Box<dyn io::BufRead>>, String> {
    let r = og_io::open(path).map_err(|e| os_error(&e, path))?;
    Ok(og_io::Lines::new(r, og_io::Utf8::Replace))
}

/// The whole text of `path`, as `open(path, errors="replace").read()`.
pub(crate) fn read_text(path: &str) -> Result<String, String> {
    let mut out = String::new();
    let mut lines = text_lines(path)?;
    let mut line = String::new();
    while lines
        .read_line(&mut line)
        .map_err(|e| format!("{path}: {e}"))?
    {
        out.push_str(&line);
    }
    Ok(out)
}

/// Writes `body` to `path` through a temporary file renamed into place.
pub(crate) fn write_file(path: &str, body: &[u8]) -> Result<(), String> {
    let write = || -> io::Result<()> {
        let mut f = og_io::AtomicFile::create(path)?;
        f.write_all(body)?;
        f.commit()
    };
    write().map_err(|e| format!("cannot write {path}: {e}"))
}

/// The exit code of a step tool: 0, or 1 after `ERROR: <message>`.
pub(crate) fn finish(result: Result<(), String>, err: &mut impl Write) -> u8 {
    match result {
        Ok(()) => 0,
        Err(e) => {
            let _ = writeln!(err, "ERROR: {e}");
            1
        }
    }
}

#[cfg(test)]
pub(crate) mod testutil {
    use std::path::PathBuf;

    /// A scratch directory, removed when dropped.
    pub struct Dir(pub PathBuf);

    impl Dir {
        pub fn new(tag: &str) -> Dir {
            let p = std::env::temp_dir().join(format!(
                "og-summary-{tag}-{}-{:?}",
                std::process::id(),
                std::thread::current().id()
            ));
            let _ = std::fs::remove_dir_all(&p);
            std::fs::create_dir_all(&p).unwrap();
            Dir(p)
        }

        /// Writes `body` to `name` (parents made) and returns its path.
        pub fn put(&self, name: &str, body: impl AsRef<[u8]>) -> String {
            let p = self.0.join(name);
            std::fs::create_dir_all(p.parent().unwrap()).unwrap();
            std::fs::write(&p, body).unwrap();
            p.display().to_string()
        }

        /// Writes `text` gzipped to `name` and returns its path.
        pub fn put_gz(&self, name: &str, text: &str) -> String {
            use std::io::Write;
            let mut e = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
            e.write_all(text.as_bytes()).unwrap();
            self.put(name, e.finish().unwrap())
        }

        pub fn path(&self, name: &str) -> String {
            self.0.join(name).display().to_string()
        }

        pub fn read(&self, name: &str) -> String {
            std::fs::read_to_string(self.0.join(name)).unwrap()
        }
    }

    impl Drop for Dir {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }

    /// Runs a step tool's `cli` and returns (exit, stdout, stderr).
    pub fn run(
        cli: impl FnOnce(&[String], &mut Vec<u8>, &mut Vec<u8>) -> u8,
        args: &[&str],
    ) -> (u8, String, String) {
        let argv: Vec<String> = args.iter().map(|s| (*s).to_owned()).collect();
        let (mut out, mut err) = (Vec::new(), Vec::new());
        let code = cli(&argv, &mut out, &mut err);
        (
            code,
            String::from_utf8(out).unwrap(),
            String::from_utf8(err).unwrap(),
        )
    }
}
