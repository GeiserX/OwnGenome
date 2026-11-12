//! Copy number from short-read depth (replaces `bin/hba_copy_number.py`,
//! `bin/lpa_kiv2.py`, the commands of `bin/copy_number.py` and the
//! copy-aware markers of `bin/traits_report.py`).
//!
//! One method, one catalog row per locus (`bin/cn_loci.tsv`,
//! `docs/cn-loci.md`): the depth over the row's target segments against the
//! depth over single-copy control segments, as `samtools bedcov` sums it,
//!
//! ```text
//! estimate = scale x T / C - offset
//! ```
//!
//! guarded by a minimum control depth, flank agreement and the share of the
//! target depth from reads of mapping quality 0, then called by the row's
//! `call` rule.
//!
//! - [`catalog`]: the catalog, read and checked once per process.
//! - [`engine`]: the guards and the call ([`engine::assess`]), and the PSV
//!   sites of a `call psv` row ([`engine::psv_assess`]).
//! - [`input`]: `samtools bedcov` and `samtools mpileup` tables, BAM headers.
//! - [`layout`]: each step's `metric<TAB>value` table and the interval table.
//! - [`markers`]: the copy-aware markers (LPA KIV-2 4925G>A, SMN1
//!   g.27134T>G) read from the pileup at every GRCh38 copy of the site.
//! - [`smn`]: the one rule that reads the SMN1 psCN at exon 7 from step 35's
//!   SMN table.
//! - [`cli`]: `og cn bed|call|intervals|check-bam|marker-regions|marker-read`.
//!
//! Each command writes what the Python wrote, byte for byte (`tests/parity`),
//! except where it fixes it on purpose; the module that does says where.

pub mod catalog;
pub mod cli;
pub mod engine;
pub mod input;
pub mod layout;
pub mod markers;
pub mod smn;

use std::fmt;
use std::io;

/// A problem the command reports as `ERROR: <text>` with its exit code.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Error(pub String);

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for Error {}

impl Error {
    /// An I/O error on `path`, worded as Python's `OSError`.
    pub fn io(e: &io::Error, path: &str) -> Error {
        if e.kind() == io::ErrorKind::InvalidData {
            Error(format!("{path}: {e}"))
        } else {
            Error(og_filters::os_error(e, path))
        }
    }
}

/// `int(s)` with Python's message when it is not one.
pub(crate) fn py_int(s: &str) -> Result<i64, Error> {
    pycompat::num::parse_int(s)
        .ok()
        .and_then(|v| i64::try_from(v).ok())
        .ok_or_else(|| {
            Error(format!(
                "invalid literal for int() with base 10: {}",
                og_filters::py_repr(s)
            ))
        })
}

/// `float(s)` with Python's message when it is not one.
pub(crate) fn py_float(s: &str) -> Result<f64, Error> {
    pycompat::num::parse_float(s).map_err(|_| {
        Error(format!(
            "could not convert string to float: {}",
            og_filters::py_repr(s)
        ))
    })
}

/// `f"{x:.{prec}f}"`.
pub(crate) fn fixed(x: f64, prec: usize) -> String {
    pycompat::fmt::fixed(x, prec)
}

/// `"." if x is None else f"{x:.2f}"`.
pub(crate) fn f2(x: Option<f64>) -> String {
    x.map_or_else(|| ".".to_owned(), |v| fixed(v, 2))
}

/// `"." if x is None else str(x)`.
pub(crate) fn d<T: fmt::Display>(x: Option<T>) -> String {
    x.map_or_else(|| ".".to_owned(), |v| v.to_string())
}

/// `round(x)` of a finite float, as an integer (half to even).
pub(crate) fn round_i(x: f64) -> i64 {
    x.round_ties_even() as i64
}

/// The text lines of a file as Python's `open(path, errors="replace")`
/// reads them (universal newlines), each without its `\n`.
pub(crate) fn read_lines(path: &str) -> Result<Vec<String>, Error> {
    let file = std::fs::File::open(path).map_err(|e| Error::io(&e, path))?;
    lines_of(io::BufReader::new(file), path)
}

/// The lines of the file at `path`, as Python's `open(path)` reads them:
/// universal newlines, strict UTF-8.
pub(crate) fn read_lines_strict(path: &str) -> Result<Vec<String>, Error> {
    let file = std::fs::File::open(path).map_err(|e| Error::io(&e, path))?;
    og_io::Lines::new(io::BufReader::new(file), og_io::Utf8::Strict)
        .map(|l| l.map_err(|e| Error::io(&e, path)))
        .collect()
}

/// The lines of `r`, as [`read_lines`].
pub(crate) fn lines_of(r: impl io::BufRead, path: &str) -> Result<Vec<String>, Error> {
    let mut out = Vec::new();
    each_line_of(r, path, |l| {
        out.push(l.to_owned());
        Ok(())
    })?;
    Ok(out)
}

/// Streams the lines of `r` (as [`read_lines`] reads them) through `each`,
/// one at a time, so a long pileup or header is never held whole.
pub(crate) fn each_line_of(
    r: impl io::BufRead,
    path: &str,
    mut each: impl FnMut(&str) -> Result<(), Error>,
) -> Result<(), Error> {
    let mut lines = og_io::Lines::new(r, og_io::Utf8::Replace);
    let mut line = String::new();
    while lines
        .read_line(&mut line)
        .map_err(|e| Error::io(&e, path))?
    {
        each(line.strip_suffix('\n').unwrap_or(&line))?;
    }
    Ok(())
}

/// [`each_line_of`] over the file at `path`.
pub(crate) fn each_line(
    path: &str,
    each: impl FnMut(&str) -> Result<(), Error>,
) -> Result<(), Error> {
    let file = std::fs::File::open(path).map_err(|e| Error::io(&e, path))?;
    each_line_of(io::BufReader::new(file), path, each)
}

/// Python's `str.strip()`.
pub(crate) fn strip(s: &str) -> &str {
    pycompat::text::strip(s, None)
}
