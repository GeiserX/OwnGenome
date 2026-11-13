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
//! - `og dev helpers` ([`helpers`]): the helper binaries each script and
//!   module calls, and the image the call runs in
//!   (`scripts/ci/check-container-helpers.sh`).
//! - `og dev remnants` ([`remnants`]): does each Python file kept only for
//!   its importers (`tests/parity/remnants.tsv`) hold nothing else?
//! - `og dev scan personal-data` ([`personal_data`]): no tracked file holds
//!   a home path, a private address or a genotype.
//! - `og dev scan images`, `og dev versions coupled` ([`images`]): no file
//!   names an image outside `versions.env`; the alignment task's image holds
//!   the releases of step 02's two images.
//! - `og dev versions check-run`, `og dev versions expand` ([`versions`]): a
//!   run used the images `conf/containers.config` names; download URLs with
//!   their `versions.env` variables filled.
//! - `og dev pgs` ([`pgs`]): the PGS Catalog's answers against
//!   `assets/pgs_scores.tsv`.

pub mod contract;
pub mod helpers;
pub mod images;
mod lookaround;
pub mod parity;
pub mod personal_data;
pub mod pgs;
pub mod remnants;
pub mod report_only;
pub mod schema;
pub mod shlex;
pub mod versions;
pub mod versions_env;
mod zip;

use std::io::Write;

const USAGE: &str = "usage: og dev parity-compare KIND ORACLE NEW\n       og dev parity-compare --self-test\n       og dev contract --static [--repo DIR] [--contract FILE]\n       og dev contract --outputs DIR --leg bash|nextflow --sample S [--observe]\n       og dev schema SCHEMA.json DOC.json\n       og dev report-only [--repo DIR]\n       og dev remnants [--list FILE] [--root DIR]\n       og dev scan personal-data|images [args]\n       og dev versions check-run|expand|coupled [args]\n       og dev pgs trait|label-match|metric-match [args]\n       og dev helpers --list|--static|--images|--orad-calls [--root DIR]\n       og dev helpers --orad --help-file FILE [--orad-version V] [--root DIR]\n";

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

fn unknown(verb: &str) -> u8 {
    eprint!("{USAGE}");
    eprintln!("og dev: error: unknown verb '{verb}'");
    2
}

/// Runs `og dev ARGS`; returns the exit code.
pub fn main(args: &[String]) -> u8 {
    match args.first().map(String::as_str) {
        Some("parity-compare") => parity::cli(&args[1..]),
        Some("remnants") => remnants::cli(&args[1..]),
        Some("contract") => run(contract::cli, &args[1..]),
        Some("schema") => run(schema::cli, &args[1..]),
        Some("report-only") => run(report_only::cli, &args[1..]),
        Some("scan") => match args.get(1).map(String::as_str) {
            Some("personal-data") => run(personal_data::cli, &args[2..]),
            Some("images") => run(images::cli_scan, &args[2..]),
            _ => unknown(args.get(1).map_or("scan", String::as_str)),
        },
        Some("versions") => match args.get(1).map(String::as_str) {
            Some("check-run") => run(versions::cli_check_run, &args[2..]),
            Some("expand") => run(versions::cli_expand, &args[2..]),
            Some("coupled") => run(images::cli_coupled, &args[2..]),
            _ => unknown(args.get(1).map_or("versions", String::as_str)),
        },
        Some("pgs") => run(pgs::cli, &args[1..]),
        Some("helpers") => run(helpers::cli, &args[1..]),
        Some("-h" | "--help") => {
            print!("{USAGE}");
            0
        }
        Some(other) => unknown(other),
        None => {
            eprint!("{USAGE}");
            eprintln!("og dev: error: the following arguments are required: verb");
            2
        }
    }
}
