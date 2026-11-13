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
//! - `og dev fixture` ([`fixture`]): writers of test fixtures (the
//!   compound-het reads, a trimmed PharmCAT report).
//! - `og dev evidence` ([`evidence`]): the measurements of the public
//!   evidence workflow on GIAB's HG002.
//! - `og dev renovate-check LOG` ([`renovate`]): what a Renovate dry run
//!   found, against every pin in the repository.
//! - `og dev build cn-intervals|callability|sv-rule-regions` ([`build`]):
//!   rebuild a committed asset table from its pinned downloads.

pub mod build;
pub mod contract;
pub mod evidence;
pub mod fixture;
pub mod helpers;
pub mod parity;
pub mod remnants;
pub mod renovate;
pub mod report_only;
pub mod schema;
pub mod shlex;
mod zip;

use std::io::Write;

const USAGE: &str = "usage: og dev parity-compare KIND ORACLE NEW\n       og dev parity-compare --self-test\n       og dev contract --static [--repo DIR] [--contract FILE]\n       og dev contract --outputs DIR --leg bash|nextflow --sample S [--observe]\n       og dev schema SCHEMA.json DOC.json\n       og dev report-only [--repo DIR]\n       og dev remnants [--list FILE] [--root DIR]\n       og dev fixture comphet|trim-pharmcat [args]\n       og dev evidence clashes|hg002-gvcf [args]\n       og dev renovate-check LOG [--pins FILE] [--repo DIR]\n       og dev helpers --list|--static|--images|--orad-calls [--root DIR]\n       og dev helpers --orad --help-file FILE [--orad-version V] [--root DIR]\n       og dev build cn-intervals|callability|sv-rule-regions [args]\n";

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
    run_with_writers(f, args, &mut out, &mut err)
}

pub(crate) fn run_with_writers(
    f: fn(&[String], &mut dyn Write, &mut dyn Write) -> u8,
    args: &[String],
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> u8 {
    let mut out = CheckedWriter {
        inner: out,
        failure: None,
    };
    let mut err = CheckedWriter {
        inner: err,
        failure: None,
    };
    let code = f(args, &mut out, &mut err);
    let _ = out.flush();
    if let Some(reason) = &out.failure {
        let _ = writeln!(err, "og dev: stdout: {reason}");
    }
    let _ = err.flush();
    if code == 0 && (out.failure.is_some() || err.failure.is_some()) {
        1
    } else {
        code
    }
}

/// Remember errors even when an older verb ignores its Write result.
struct CheckedWriter<'a> {
    inner: &'a mut dyn Write,
    failure: Option<String>,
}
impl Write for CheckedWriter<'_> {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        let result = match self.inner.write(bytes) {
            Ok(0) if !bytes.is_empty() => Err(std::io::Error::from(std::io::ErrorKind::WriteZero)),
            result => result,
        };
        if let Err(e) = &result {
            self.failure.get_or_insert_with(|| e.to_string());
        }
        result
    }
    fn flush(&mut self) -> std::io::Result<()> {
        let result = self.inner.flush();
        if let Err(e) = &result {
            self.failure.get_or_insert_with(|| e.to_string());
        }
        result
    }
}

/// Runs `og dev ARGS`; returns the exit code.
pub fn main(args: &[String]) -> u8 {
    match args.first().map(String::as_str) {
        Some("parity-compare") => parity::cli(&args[1..]),
        Some("remnants") => remnants::cli(&args[1..]),
        Some("contract") => run(contract::cli, &args[1..]),
        Some("schema") => run(schema::cli, &args[1..]),
        Some("report-only") => run(report_only::cli, &args[1..]),
        Some("fixture") => run(fixture::cli, &args[1..]),
        Some("evidence") => run(evidence::cli, &args[1..]),
        Some("renovate-check") => run(renovate::cli, &args[1..]),
        Some("build") => run(build::cli, &args[1..]),
        Some("helpers") => run(helpers::cli, &args[1..]),
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

#[cfg(test)]
mod stream_tests {
    use super::*;
    struct Fails {
        flush_only: bool,
    }
    impl Write for Fails {
        fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
            if self.flush_only {
                Ok(bytes.len())
            } else {
                Err(std::io::Error::other("write control"))
            }
        }
        fn flush(&mut self) -> std::io::Result<()> {
            if self.flush_only {
                Err(std::io::Error::other("flush control"))
            } else {
                Ok(())
            }
        }
    }
    pub(crate) fn assert_status_errors(
        verb: fn(&[String], &mut dyn Write, &mut dyn Write) -> u8,
        args: &[String],
        writes: bool,
    ) {
        for flush_only in [false, true] {
            let mut bad = Fails { flush_only };
            let mut err = Vec::new();
            assert_eq!(
                run_with_writers(verb, args, &mut bad, &mut err),
                u8::from(writes || flush_only),
                "{}",
                String::from_utf8_lossy(&err)
            );
        }
    }

    #[test]
    fn added_verbs_propagate_write_and_flush_errors() {
        let verbs = [fixture::cli, evidence::cli, renovate::cli];
        for verb in verbs {
            for flush_only in [false, true] {
                let mut bad = Fails { flush_only };
                let mut err = Vec::new();
                assert_eq!(
                    run_with_writers(verb, &["--help".into()], &mut bad, &mut err),
                    1
                );
                assert!(String::from_utf8(err).unwrap().contains("stdout:"));
            }
        }
    }
    #[test]
    fn shared_runner_checks_stderr_and_an_empty_stdout() {
        fn clean(_: &[String], _: &mut dyn Write, err: &mut dyn Write) -> u8 {
            let _ = writeln!(err, "status");
            0
        }
        for flush_only in [false, true] {
            assert_eq!(
                run_with_writers(clean, &[], &mut Vec::new(), &mut Fails { flush_only }),
                1
            );
        }
        assert_eq!(
            run_with_writers(clean, &[], &mut Fails { flush_only: true }, &mut Vec::new()),
            1
        );
    }
}
