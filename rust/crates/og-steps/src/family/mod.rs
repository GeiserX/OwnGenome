//! `og family`: the family analysis of step 38 (scripts/38-family.sh) and
//! of the Nextflow FAMILY workflow, port of `bin/family_report.py` and
//! `bin/family_segments.py`.
//!
//! - `og family plan`: which family analyses run, are skipped or wait for a
//!   member's own run (`plan.rs`);
//! - `og family prepare`: members, somalier's and WhatsHap's pedigrees and
//!   the plan for the members with a gVCF (`prepare.rs`);
//! - `og family check`: the relatedness somalier measured against the one
//!   the pedigree implies (`check.rs`);
//! - `og family report`: the allele each parent passed on, de novo
//!   candidates, compound heterozygotes and the chromosome map, from one
//!   read of the joint VCF (`report.rs`);
//! - `og family segments`: the chromosome map alone ([`segments`]).
//!
//! The pedigree rules live in `og-pedigree`. Exit codes are the Python's:
//! 2 for a problem in the input, 1 when `check` finds a pair that
//! disagrees with the pedigree.

mod args;
mod check;
mod plan;
mod prepare;
mod report;
pub mod segments;
mod vcf;

use og_pedigree::Die;

const USAGE: &str = "usage: og family {plan,prepare,check,report,segments} ...\n";

/// Runs `og family ARGS`; returns the exit code.
pub fn cli(args: &[String]) -> u8 {
    let Some(verb) = args.first() else {
        eprint!("{USAGE}");
        eprintln!("og family: error: the following arguments are required: verb");
        return 2;
    };
    let rest = &args[1..];
    let out = match verb.as_str() {
        "plan" => args::parse(&plan::SPEC, rest).and_then(|a| plan::run(&a)),
        "prepare" => args::parse(&prepare::SPEC, rest).and_then(|a| prepare::run(&a)),
        "check" => args::parse(&check::SPEC, rest).and_then(|a| check::run(&a)),
        "report" => args::parse(&report::SPEC, rest).and_then(|a| report::run(&a)),
        "segments" => args::parse(&segments::SPEC, rest).and_then(|a| segments::run(&a)),
        "-h" | "--help" => {
            print!("{USAGE}");
            return 0;
        }
        other => {
            eprint!("{USAGE}");
            eprintln!("og family: error: unknown verb '{other}'");
            return 2;
        }
    };
    match out {
        Ok(()) => 0,
        Err(Stop::Die(d)) => d.report(),
        Err(Stop::Exit(code)) => code,
    }
}

/// Why a verb stopped early.
#[derive(Debug)]
pub enum Stop {
    /// `ERROR:` lines and an exit code.
    Die(Die),
    /// Already reported (usage, `--help`): just this exit code.
    Exit(u8),
}

impl From<Die> for Stop {
    fn from(d: Die) -> Stop {
        Stop::Die(d)
    }
}

type Result<T> = std::result::Result<T, Stop>;

/// Writes `text` to `path` through a temporary file renamed into place.
fn write_file(path: &std::path::Path, text: &[u8]) -> std::result::Result<(), Die> {
    use std::io::Write;
    let shown = path.display().to_string();
    let fail = |e: std::io::Error| Die::one(format!("cannot write {shown}: {e}"));
    let mut f = og_io::AtomicFile::create(path).map_err(fail)?;
    f.write_all(text).map_err(fail)?;
    f.commit().map_err(fail)
}

/// Prints `text` on stdout; a closed pipe is not an error.
fn print_out(text: &str) {
    use std::io::Write;
    let _ = std::io::stdout().lock().write_all(text.as_bytes());
}
