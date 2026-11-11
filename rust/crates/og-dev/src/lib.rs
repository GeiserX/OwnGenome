//! `og dev`: verbs for development and CI, never run by a pipeline step.
//!
//! - `og dev parity-compare KIND ORACLE NEW` ([`parity`]): is the file a
//!   port wrote the same, by meaning, as the one the Python it replaces
//!   wrote? `scripts/ci/oracle-diff.sh` calls it for every output of every
//!   case in `tests/parity/cases.tsv`.
//! - `og dev contract` ([`contract`]): the output contract against the code
//!   (`--static`) and against a run's results (`--outputs`).
//! - `og dev schema SCHEMA DOC` ([`schema`]): a JSON document against its
//!   JSON Schema (`summary.json` against `tests/schema/summary.schema.json`).
//! - `og dev report-only` ([`report_only`]): every process that feeds only
//!   the report carries the report-only error strategy.

pub mod contract;
pub mod parity;
pub mod report_only;
pub mod schema;
mod zip;

use std::io::Write;

const USAGE: &str = "usage: og dev parity-compare KIND ORACLE NEW\n       og dev parity-compare --self-test\n       og dev contract --static [--repo DIR] [--contract FILE]\n       og dev contract --outputs DIR --leg bash|nextflow --sample S [--observe]\n       og dev schema SCHEMA.json DOC.json\n       og dev report-only [--repo DIR]\n";

/// Python's `str(OSError)`: `[Errno 2] No such file or directory: 'path'`.
pub(crate) fn py_os_error(e: &std::io::Error, path: &str) -> String {
    let quoted = pycompat::value::py_repr(&pycompat::json::Value::Str(path.to_owned()))
        .unwrap_or_else(|_| format!("'{path}'"));
    match e.raw_os_error() {
        Some(code) => {
            let text = e.to_string();
            let msg = text
                .strip_suffix(&format!(" (os error {code})"))
                .unwrap_or(&text);
            format!("[Errno {code}] {msg}: {quoted}")
        }
        None => format!("{e}: {quoted}"),
    }
}

/// Runs a verb that writes to stdout and stderr, flushing both.
fn run(f: fn(&[String], &mut dyn Write, &mut dyn Write) -> u8, args: &[String]) -> u8 {
    let mut out = std::io::stdout().lock();
    let mut err = std::io::stderr().lock();
    let code = f(args, &mut out, &mut err);
    let _ = out.flush();
    code
}

/// Runs `og dev ARGS`; returns the exit code.
pub fn main(args: &[String]) -> u8 {
    match args.first().map(String::as_str) {
        Some("parity-compare") => parity::cli(&args[1..]),
        Some("contract") => run(contract::cli, &args[1..]),
        Some("schema") => run(schema::cli, &args[1..]),
        Some("report-only") => run(report_only::cli, &args[1..]),
        Some("-h" | "--help") => {
            print!("{USAGE}");
            0
        }
        Some(other) => {
            eprint!("{USAGE}");
            eprintln!("og dev: error: unknown verb '{other}'");
            2
        }
        None => {
            eprint!("{USAGE}");
            eprintln!("og dev: error: the following arguments are required: verb");
            2
        }
    }
}
