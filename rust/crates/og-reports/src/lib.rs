//! The opt-in step reports: each one reads one third-party tool's output
//! and writes OwnGenome's table of it.
//!
//! - [`vntyper`] (`og vntyper report`, step 52): VNtyper's MUC1 VNTR call.
//! - [`locityper`] (`og locityper report`, step 53): Locityper's genotypes.
//! - [`mocha`] (`og mocha report`, step 51): MoChA's mosaic alterations.
//! - [`mei`] (`og mei report`, `og mei annotsv-input`, step 50): MEGAnE's
//!   mobile element insertions and the known disease insertions.
//! - [`viral`] (`og viral report`, step 47): viral DNA in the unmapped reads.
//! - [`ystr`] (`og y str-bed`, `og y str-report`, step 45): Y-STR counts.
//! - [`exomiser`] (`og exomiser prepare`, `og exomiser report`, step 48).
//!
//! Each one writes what the Python it replaced wrote, cell for cell
//! (`tests/parity`), except where its module names a fix. What every tool
//! changes on purpose: an input it cannot read is one `ERROR:` line on
//! stderr and the tool's exit code, never a traceback, and each output is
//! written whole, through a temporary file and a rename, or not at all.

pub mod exomiser;
pub mod locityper;
pub mod mei;
pub mod mocha;
pub mod viral;
pub mod vntyper;
pub mod ystr;

use og_filters::args::{self, Parsed, Spec};
use og_io::Utf8;
use std::fmt;
use std::io::{self, Write};

/// Why a tool stopped: the message after `ERROR: `.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Error(pub String);

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

impl From<og_tsv::Error> for Error {
    /// A file that cannot be opened is worded as Python's `OSError`; any
    /// other error as og-tsv words it.
    fn from(e: og_tsv::Error) -> Self {
        match &e {
            og_tsv::Error::Io {
                path,
                line: None,
                source,
            } => Error::io(source, path),
            _ => Error(e.to_string()),
        }
    }
}

impl Error {
    /// An I/O error on `path`, worded as Python's `OSError`
    /// (`[Errno 2] No such file or directory: 'x.tsv'`).
    pub fn io(e: &io::Error, path: &str) -> Error {
        if e.kind() == io::ErrorKind::InvalidData {
            Error(format!("{path}: {e}"))
        } else {
            Error(og_filters::os_error(e, path))
        }
    }
}

/// The lines of `path` as Python's `open()` iterates them, one at a time:
/// universal newlines, each line with its `\n` (none on a last line without
/// one). Plain or gzip. A read error ends the iteration with that error.
pub(crate) fn lines(path: &str, utf8: Utf8) -> Result<LineIter, Error> {
    let file = og_io::open(path).map_err(|e| Error::io(&e, path))?;
    Ok(LineIter {
        lines: Some(og_io::Lines::new(file, utf8)),
        path: path.to_owned(),
    })
}

/// [`lines`]: the lines of one file.
pub(crate) struct LineIter {
    lines: Option<og_io::Lines<Box<dyn io::BufRead>>>,
    path: String,
}

impl Iterator for LineIter {
    type Item = Result<String, Error>;

    fn next(&mut self) -> Option<Self::Item> {
        let lines = self.lines.as_mut()?;
        let mut line = String::new();
        match lines.read_line(&mut line) {
            Ok(true) => Some(Ok(line)),
            Ok(false) => {
                self.lines = None;
                None
            }
            Err(e) => {
                self.lines = None;
                Some(Err(Error::io(&e, &self.path)))
            }
        }
    }
}

/// Every line of a small file, as [`lines`] reads them.
pub(crate) fn read_lines(path: &str, utf8: Utf8) -> Result<Vec<String>, Error> {
    lines(path, utf8)?.collect()
}

/// A header-keyed TSV, read whole: its header (empty for an empty file)
/// and its rows.
pub(crate) fn read_table(
    path: &str,
    opts: og_tsv::Options,
) -> Result<(Vec<String>, Vec<og_tsv::Row>), Error> {
    let mut t = opts.open(path)?;
    let header = t.header()?.map(<[String]>::to_vec).unwrap_or_default();
    let rows = t.collect::<Result<Vec<_>, _>>()?;
    Ok((header, rows))
}

/// Writes `body` to `path` whole, through a temporary file and a rename.
pub(crate) fn write_atomic(path: &str, body: &[u8]) -> Result<(), Error> {
    let fail = |e: io::Error| Error::io(&e, path);
    let mut f = og_io::AtomicFile::create(path).map_err(fail)?;
    f.write_all(body).map_err(fail)?;
    f.commit().map_err(fail)
}

/// Python's `s.strip()`.
pub(crate) fn strip(s: &str) -> &str {
    pycompat::text::strip(s, None)
}

/// `f"{x:.{prec}f}"`.
pub(crate) fn fixed(x: f64, prec: usize) -> String {
    pycompat::fmt::fixed(x, prec)
}

/// `s.isdigit()` for ASCII digits only. Python's also takes other Unicode
/// digits (`²`), which its `int()` then refuses with a traceback; here they
/// are not digits.
pub(crate) fn digits(s: &str) -> bool {
    !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit())
}

/// `int(s)` for a string [`digits`] accepted.
pub(crate) fn int(s: &str) -> Option<i128> {
    pycompat::num::parse_int(s).ok()
}

/// `float(s)`, or `None` where Python raises `ValueError`.
pub(crate) fn float(s: &str) -> Option<f64> {
    pycompat::num::parse_float(s).ok()
}

/// What a tool's run gave: its stdout text, or the message and exit code it
/// stops with.
pub(crate) type Outcome = Result<String, (Error, u8)>;

/// Parses `rest` by `spec` and runs `body`, printing its stdout or its
/// `ERROR:` line; the exit code.
pub(crate) fn run(
    spec: &Spec,
    rest: &[String],
    out: &mut dyn Write,
    err: &mut dyn Write,
    body: impl FnOnce(&Parsed) -> Result<Outcome, args::ArgError>,
) -> u8 {
    let a = match args::parse(spec, rest) {
        Ok(a) => a,
        Err(e) => return args::report(spec, &e),
    };
    match body(&a) {
        Err(e) => args::report(spec, &e),
        Ok(Ok(text)) => match out.write_all(text.as_bytes()).and_then(|()| out.flush()) {
            Ok(()) => 0,
            Err(e) => {
                let _ = writeln!(err, "ERROR: cannot write stdout: {e}");
                1
            }
        },
        Ok(Err((e, code))) => {
            let _ = writeln!(err, "ERROR: {e}");
            code
        }
    }
}

/// `--name` as argparse's `type=int` reads it, or `default`.
pub(crate) fn int_opt(a: &Parsed, name: &str, default: i64) -> Result<i64, args::ArgError> {
    a.typed(name, default, "int", |s| {
        args::int(s).and_then(|v| i64::try_from(v).ok())
    })
}

/// The value of a required option (the parser has checked it is there).
pub(crate) fn opt<'a>(a: &'a Parsed, name: &str) -> &'a str {
    a.get(name).unwrap_or_default()
}

/// The usage line of an area with several verbs, and its refusal.
pub(crate) fn unknown_verb(area: &str, verbs: &str, verb: Option<&str>, err: &mut dyn Write) -> u8 {
    let _ = writeln!(err, "usage: og {area} {verbs} [args]");
    match verb {
        Some(v) => {
            let _ = writeln!(err, "og {area}: error: unknown command '{v}'");
        }
        None => {
            let _ = writeln!(
                err,
                "og {area}: error: the following arguments are required: command"
            );
        }
    }
    2
}

#[cfg(test)]
pub(crate) mod testutil {
    use std::path::{Path, PathBuf};

    /// A scratch directory removed when dropped.
    pub struct Dir(pub PathBuf);

    impl Dir {
        pub fn new(tag: &str) -> Dir {
            let p = std::env::temp_dir().join(format!(
                "og-reports-{tag}-{}-{:?}",
                std::process::id(),
                std::thread::current().id()
            ));
            let _ = std::fs::remove_dir_all(&p);
            std::fs::create_dir_all(&p).unwrap();
            Dir(p)
        }

        /// Writes `text` to `name` (parents made) and returns its path.
        pub fn put(&self, name: &str, text: &str) -> String {
            let p = self.0.join(name);
            std::fs::create_dir_all(p.parent().unwrap()).unwrap();
            std::fs::write(&p, text).unwrap();
            p.display().to_string()
        }

        pub fn path(&self, name: &str) -> String {
            self.0.join(name).display().to_string()
        }

        pub fn read(&self, name: &str) -> String {
            std::fs::read_to_string(self.0.join(name)).unwrap()
        }

        pub fn root(&self) -> &Path {
            &self.0
        }
    }

    impl Drop for Dir {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }
}
