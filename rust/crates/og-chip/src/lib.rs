//! Clonal haematopoiesis (CHIP) of step 49 (replaces the commands of
//! `bin/chip_report.py`):
//!
//! - [`prepare`] (`og chip prepare`): each table gene's MANE Select
//!   transcript from Ensembl's GFF3, as the GFF3 `bcftools csq` reads, and
//!   its coding exons, 10 bases wider on each side, as a BED.
//! - [`report`] (`og chip report`): each protein change of a table gene with
//!   its class and status, and per gene the depth and the detection limit.
//! - [`table`] (`og chip check-table`): the schema of `assets/chip_genes.tsv`.
//!
//! Each command writes what the Python wrote, byte for byte (`tests/parity`).
//! What changes on purpose: every input problem is an `ERROR:` line and exit
//! 1, never a traceback, and an output is written whole or not at all.

pub mod binom;
pub mod prepare;
pub mod report;
pub mod table;

use og_filters::args::{self, Spec};
use std::fmt;
use std::io::{self, Write};

/// Reads a CHIP gene needs with the change, for a call and for the
/// detection limit.
pub const MIN_ALT_READS: u32 = 3;
/// The detection limit is the VAF found this often at the mean depth.
pub const DETECT_POWER: f64 = 0.8;

/// A problem with an input or an output: printed as `ERROR: <text>`, exit 1.
#[derive(Debug, PartialEq, Eq)]
pub struct Error(pub String);

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

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

/// `int(s)`, or an error naming the file and the value.
pub(crate) fn int(s: &str, path: &str) -> Result<i128, Error> {
    pycompat::num::parse_int(s).map_err(|_| {
        Error(format!(
            "{path}: {} is not a whole number",
            og_filters::py_repr(s)
        ))
    })
}

/// The text lines of `path` as Python's `open()` reads them, each with its
/// `\n` (none on a last line without one), passed to `each` in order.
pub(crate) fn each_line(
    path: &str,
    mut each: impl FnMut(&str) -> Result<(), Error>,
) -> Result<(), Error> {
    let file = og_io::open(path).map_err(|e| Error::io(&e, path))?;
    let mut lines = og_io::Lines::new(file, og_io::Utf8::Strict);
    let mut line = String::new();
    while lines
        .read_line(&mut line)
        .map_err(|e| Error::io(&e, path))?
    {
        each(&line)?;
    }
    Ok(())
}

/// Writes `body` to `path` whole, through a temporary file and a rename.
pub(crate) fn write_atomic(path: &str, body: &[u8]) -> Result<(), Error> {
    let fail = |e: io::Error| Error(format!("cannot write {path}: {e}"));
    let mut f = og_io::AtomicFile::create(path).map_err(fail)?;
    f.write_all(body).map_err(fail)?;
    f.commit().map_err(fail)
}

const PREPARE: Spec = Spec {
    prog: "og chip prepare",
    usage: "usage: og chip prepare --gff GFF --genes GENES --fai FAI --out-gff OUT_GFF --out-bed OUT_BED\n",
    options: &["--gff", "--genes", "--fai", "--out-gff", "--out-bed"],
    positionals: &[],
    required: &["--gff", "--genes", "--fai", "--out-gff", "--out-bed"],
};

const REPORT: Spec = Spec {
    prog: "og chip report",
    usage: "usage: og chip report --genes GENES --gff GFF --bed BED --calls CALLS --depth DEPTH --sample SAMPLE --out OUT --depth-out DEPTH_OUT\n",
    options: &[
        "--genes",
        "--gff",
        "--bed",
        "--calls",
        "--depth",
        "--sample",
        "--out",
        "--depth-out",
    ],
    positionals: &[],
    required: &[
        "--genes",
        "--gff",
        "--bed",
        "--calls",
        "--depth",
        "--sample",
        "--out",
        "--depth-out",
    ],
};

const CHECK: Spec = Spec {
    prog: "og chip check-table",
    usage: "usage: og chip check-table TABLE\n",
    options: &[],
    positionals: &["table"],
    required: &[],
};

const USAGE: &str = "usage: og chip prepare|report|check-table [args]\n";

/// `og chip VERB ARGS`: the exit code (0, 1 for a problem, 2 for a refused
/// command line).
pub fn cli(verb: &str, rest: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let spec = match verb {
        "prepare" => &PREPARE,
        "report" => &REPORT,
        "check-table" => &CHECK,
        "-h" | "--help" => {
            let _ = out.write_all(USAGE.as_bytes());
            return 0;
        }
        other => {
            let _ = writeln!(err, "{USAGE}og chip: error: unknown command '{other}'");
            return 2;
        }
    };
    let a = match args::parse(spec, rest) {
        Ok(a) => a,
        Err(e) => return args::report(spec, &e),
    };
    let opt = |name: &str| a.get(name).unwrap_or_default().to_owned();
    let done = match verb {
        "check-table" => {
            let path = &a.pos[0];
            match table::check_table(path) {
                Ok(problems) => {
                    for p in &problems {
                        let _ = writeln!(err, "{path}: {p}");
                    }
                    return u8::from(!problems.is_empty());
                }
                Err(e) => Err(e),
            }
        }
        "prepare" => checked_table(&opt("--genes")).and_then(|()| {
            prepare::prepare(
                &opt("--gff"),
                &opt("--genes"),
                &opt("--fai"),
                &opt("--out-gff"),
                &opt("--out-bed"),
            )
        }),
        _ => checked_table(&opt("--genes")).and_then(|()| {
            let files = report::Files {
                genes: opt("--genes"),
                gff: opt("--gff"),
                bed: opt("--bed"),
                calls: opt("--calls"),
                depth: opt("--depth"),
                out: opt("--out"),
                depth_out: opt("--depth-out"),
            };
            let lines = report::report(&files, &opt("--sample"))?;
            out.write_all(lines.as_bytes())
                .map_err(|e| Error(format!("cannot write stdout: {e}")))
        }),
    };
    match done {
        Ok(()) => 0,
        Err(e) => {
            let _ = writeln!(err, "ERROR: {e}");
            1
        }
    }
}

/// prepare and report stop on a table that fails its schema check.
fn checked_table(path: &str) -> Result<(), Error> {
    let problems = table::check_table(path)?;
    if problems.is_empty() {
        Ok(())
    } else {
        Err(Error(format!("{path}: {}", problems.join("; "))))
    }
}
