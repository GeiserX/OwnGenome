//! `og dev`: verbs for development and CI, never run by a pipeline step.
//!
//! - `og dev parity-compare KIND ORACLE NEW` ([`parity`]): is the file a
//!   port wrote the same, by meaning, as the one the Python it replaces
//!   wrote? `scripts/ci/oracle-diff.sh` calls it for every output of every
//!   case in `tests/parity/cases.tsv`.
//! - `og dev remnants` ([`remnants`]): does each Python file kept only for
//!   its importers (`tests/parity/remnants.tsv`) hold nothing else?

pub mod parity;
pub mod remnants;
mod zip;

const USAGE: &str = "usage: og dev parity-compare KIND ORACLE NEW\n       og dev parity-compare --self-test\n       og dev remnants [--list FILE] [--root DIR]\n";

/// Runs `og dev ARGS`; returns the exit code.
pub fn main(args: &[String]) -> u8 {
    match args.first().map(String::as_str) {
        Some("parity-compare") => parity::cli(&args[1..]),
        Some("remnants") => remnants::cli(&args[1..]),
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
