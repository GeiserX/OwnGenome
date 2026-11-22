//! `og dev fixture <kind>`: the inputs the clinical, consumer-export and
//! HLA e2e cases make for themselves, and the two consumer-export checks
//! that read the gVCF a second, independent way.
//!
//! Each kind replaces a Python heredoc of a case under `tests/e2e/`. It takes
//! the arguments the heredoc took, in the same order, and writes the same
//! bytes, so a case keeps every count and threshold it asserts.
//!
//! Exit 0 when written, 1 when an input cannot be read or has the wrong
//! shape, 2 on a usage error.

use std::fs::File;
use std::io::{self, BufRead, BufWriter, Write};

mod consumer;
mod hla;
mod vep;

/// Why a kind stopped.
#[derive(Debug)]
pub enum Fail {
    /// The arguments are wrong (exit 2).
    Usage(String),
    /// An input is missing or has the wrong shape (exit 1).
    Run(String),
}

impl From<io::Error> for Fail {
    fn from(e: io::Error) -> Self {
        Fail::Run(e.to_string())
    }
}

/// One kind: its name, its usage line and the function that runs it.
type Run = fn(&[String], &mut dyn Write, &mut dyn Write) -> Result<(), Fail>;

const TABLE: &[(&str, &str, Run)] = &[
    (
        "vep-af",
        "HITS.tsv OUT.vcf CHROM POS [--filler] [--want WANT.tsv] [--rsids]",
        vep::af,
    ),
    (
        "consumer-expected",
        "GVCF.gz SITES.tsv.gz PAR.bed",
        consumer::expected,
    ),
    (
        "consumer-roundtrip",
        "LIFTED.tsv CALLS.tsv SITES.tsv.gz CHIP.txt",
        consumer::roundtrip,
    ),
    ("hla-spike", "HLA_DAT OUT_DIR ALLELE...", hla::spike),
];

/// The usage lines of every kind, one per line.
pub fn usage() -> String {
    TABLE
        .iter()
        .map(|(name, args, _)| format!("       og dev fixture {name} {args}\n"))
        .collect()
}

/// Runs `og dev fixture KIND ARGS` when KIND is one of these kinds; `None`
/// for a kind this module does not have.
pub fn cli(kind: &str, args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> Option<u8> {
    let (name, line, run) = TABLE.iter().find(|(name, _, _)| *name == kind)?;
    if matches!(args, [a] if a == "-h" || a == "--help") {
        return Some(match writeln!(out, "usage: og dev fixture {name} {line}") {
            Ok(()) => 0,
            Err(e) => {
                let _ = writeln!(err, "og dev fixture {name}: {e}");
                1
            }
        });
    }
    Some(
        match run(args, out, err).and_then(|()| out.flush().map_err(Fail::from)) {
            Ok(()) => 0,
            Err(Fail::Usage(msg)) => {
                let _ = writeln!(
                    err,
                    "usage: og dev fixture {name} {line}\nog dev fixture {name}: error: {msg}"
                );
                2
            }
            Err(Fail::Run(msg)) => {
                let _ = writeln!(err, "og dev fixture {name}: {msg}");
                1
            }
        },
    )
}

/// A [`Fail::Run`] from any message.
fn bad(msg: impl Into<String>) -> Fail {
    Fail::Run(msg.into())
}

/// `int(s)` of a table cell (Python's `int` allows surrounding blanks), or
/// an error naming the cell.
fn cell_int(s: &str) -> Result<i64, Fail> {
    s.trim()
        .parse()
        .map_err(|_| bad(format!("'{s}' is not a whole number")))
}

/// Calls `f` on each line of `path` (plain or gzip) without its `\n`, as
/// `for l in open(path): l.rstrip("\n")` reads it.
fn each_line(path: &str, mut f: impl FnMut(&str) -> Result<(), Fail>) -> Result<(), Fail> {
    let mut input = og_io::open(path).map_err(|e| bad(format!("{path}: {e}")))?;
    let mut buf = Vec::new();
    loop {
        buf.clear();
        if input
            .read_until(b'\n', &mut buf)
            .map_err(|e| bad(format!("{path}: {e}")))?
            == 0
        {
            return Ok(());
        }
        if buf.last() == Some(&b'\n') {
            buf.pop();
        }
        let line = std::str::from_utf8(&buf).map_err(|_| bad(format!("{path}: not UTF-8")))?;
        f(line)?;
    }
}

/// A file to write, buffered, with the path in the error.
fn create(path: &str) -> Result<BufWriter<File>, Fail> {
    File::create(path)
        .map(BufWriter::new)
        .map_err(|e| bad(format!("{path}: {e}")))
}

/// The parts of a GT, split on `/` and `|` as `re.split(r"[/|]", gt)` splits
/// it (empty parts kept).
fn gt_parts(gt: &str) -> Vec<&str> {
    gt.split(['/', '|']).collect()
}

/// Python's `repr` of a list of plain strings: `['A', 'G']`.
fn repr_list(items: &[String]) -> String {
    let quoted: Vec<String> = items.iter().map(|x| format!("'{x}'")).collect();
    format!("[{}]", quoted.join(", "))
}

#[cfg(test)]
pub(crate) mod testutil {
    use std::path::PathBuf;

    /// A fresh scratch folder for one test.
    pub fn dir(name: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!("og-clinical-{name}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).expect("scratch");
        d
    }

    /// Runs a kind; returns its exit code, stdout and stderr.
    pub fn run(kind: &str, args: &[&str]) -> (u8, String, String) {
        let args: Vec<String> = args.iter().map(|s| s.to_string()).collect();
        let (mut out, mut err) = (Vec::new(), Vec::new());
        let code = super::cli(kind, &args, &mut out, &mut err).expect("a kind of this module");
        (
            code,
            String::from_utf8(out).expect("utf8"),
            String::from_utf8(err).expect("utf8"),
        )
    }

    /// A path in `d` as a string argument.
    pub fn arg(d: &std::path::Path, name: &str) -> String {
        d.join(name).to_str().expect("utf8").to_string()
    }

    /// Writes `text` gzipped to `path`.
    pub fn gz(path: &std::path::Path, text: &str) {
        use std::io::Write;
        let f = std::fs::File::create(path).expect("create");
        let mut w = flate2::write::GzEncoder::new(f, flate2::Compression::default());
        w.write_all(text.as_bytes()).expect("write");
        w.finish().expect("finish");
    }
}

#[cfg(test)]
mod tests {
    use super::testutil::run;

    #[test]
    fn every_kind_prints_its_usage_and_refuses_bad_arguments() {
        for (name, line, _) in super::TABLE {
            let (code, out, _) = run(name, &["--help"]);
            assert_eq!(code, 0, "{name}");
            assert_eq!(out, format!("usage: og dev fixture {name} {line}\n"));
            let (code, _, err) = run(name, &[]);
            assert_eq!(code, 2, "{name}: {err}");
            assert!(
                err.starts_with(&format!("usage: og dev fixture {name} ")),
                "{err}"
            );
        }
    }

    #[test]
    fn fixture_help_and_an_unknown_kind_list_every_clinical_kind() {
        for (args, want) in [(["--help"], 0), (["bogus"], 2)] {
            let args: Vec<String> = args.iter().map(|s| s.to_string()).collect();
            let (mut out, mut err) = (Vec::new(), Vec::new());
            assert_eq!(crate::fixture::cli(&args, &mut out, &mut err), want);
            let text = String::from_utf8([out, err].concat()).expect("utf8");
            for (name, line, _) in super::TABLE {
                let usage = format!("og dev fixture {name} {line}\n");
                assert!(text.contains(&usage), "{name}: {text}");
            }
        }
    }

    #[test]
    fn a_kind_of_another_module_is_not_taken() {
        let (mut out, mut err) = (Vec::new(), Vec::new());
        assert!(super::cli("rfc1", &[], &mut out, &mut err).is_none());
    }
}
