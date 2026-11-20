//! `og sv survivor-prep`: one caller's structural variants made ready for
//! SURVIVOR's merge in step 22, so that only a call the sample carries
//! counts as that caller's vote.
//!
//! SURVIVOR counts a caller as support for every record of its input,
//! whatever the genotype says. Before this verb, step 22 kept every record
//! whose FILTER was `PASS` or `.`, so a PASS Delly record genotyped `0/0`
//! (Delly tested the site and found the reference), a Delly call whose own
//! sample filter (FORMAT `FT`) says `LowQual`, or a CNVpytor segment no
//! quality rule had looked at, each made a Manta-only call read "2
//! callers". A record gives no vote when:
//!
//! - its FILTER is not `PASS` or `.`;
//! - its first sample's `FT` is there and is not `PASS` or `.`;
//! - its first sample's genotype calls an allele but no ALT: `0/0`, `0|0`,
//!   haploid `0`, `0/.`. A no-call (`./.`), a record without `GT` and a
//!   sites-only record keep their vote: nothing in them says the call is
//!   absent, and SURVIVOR writes `./.` itself;
//! - for `--caller cnvpytor`, its e-value (`pytorP1`) is not below
//!   [`MAX_EVAL`] or its share of MAPQ 0 reads (`pytorQ0`) is not below
//!   [`MAX_Q0`]: the rule step 05's `og sv add-callers` applies to the same
//!   calls. A missing or unreadable value is no evidence, so the call does
//!   not vote.
//!
//! The kept records are written as plain text with one sample column named
//! after the caller: SURVIVOR names its output columns after the input
//! samples, and three inputs that all name the sample would give one name
//! three times. A sites-only input gets `GT` `./.`. The header passes
//! through; a header without `GT` or without the CNVpytor tags is read the
//! same as one with them. A summary line goes to stderr.

use std::fmt;
use std::fs::File;
use std::io::{self, BufRead, BufWriter, Write};

use og_filters::args::{self, ArgError, Spec};

use super::add_callers::{MAX_EVAL, MAX_Q0};
use super::{first_gt, info_value};

const USAGE: &str = "usage: og sv survivor-prep --caller NAME IN.vcf[.gz]|- OUT.vcf|-
";

const SPEC: Spec = Spec {
    prog: "og sv survivor-prep",
    usage: USAGE,
    options: &["--caller"],
    positionals: &["IN", "OUT"],
    required: &["--caller"],
};

const GT_HEADER: &str = "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">";

/// Why a record gives no vote.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NoVote {
    /// FILTER is not `PASS` or `.`.
    Filter,
    /// The first sample's `FT` is not `PASS` or `.`.
    SampleFilter,
    /// The first sample's genotype calls an allele, none of them ALT.
    Reference,
    /// A CNVpytor call that fails the quality rule.
    Quality,
}

/// Whether the genotype calls an allele (`0`, `1`, ...) but no ALT: `0/0`,
/// `0|0`, `0`, `0/.`. `./.`, `.` and an empty value call none.
pub fn reference_only(gt: &str) -> bool {
    let mut called_ref = false;
    for a in gt.split(['/', '|']) {
        if a.is_empty() || !a.bytes().all(|b| b.is_ascii_digit()) {
            continue;
        }
        if a.bytes().all(|b| b == b'0') {
            called_ref = true;
        } else {
            return false;
        }
    }
    called_ref
}

/// The value of the first sample's FORMAT key `key`, when FORMAT has it.
fn sample_value<'a>(fields: &[&'a str], key: &str) -> Option<&'a str> {
    if fields.len() < 10 {
        return None;
    }
    let i = fields[8].split(':').position(|k| k == key)?;
    Some(fields[9].split(':').nth(i).unwrap_or("."))
}

/// A number of INFO, or NaN (missing, a flag or not a number), which no
/// threshold passes.
fn info_number(info: &str, key: &str) -> f64 {
    pycompat::num::parse_float(info_value(info, key)).unwrap_or(f64::NAN)
}

/// Why the record `fields` (at least 8) gives no vote, or `None` when it
/// does. `cnvpytor` applies CNVpytor's quality rule.
pub fn verdict(fields: &[&str], cnvpytor: bool) -> Option<NoVote> {
    if !matches!(fields[6], "PASS" | ".") {
        return Some(NoVote::Filter);
    }
    if let Some(ft) = sample_value(fields, "FT")
        && !matches!(ft, "PASS" | "." | "")
    {
        return Some(NoVote::SampleFilter);
    }
    if fields.len() >= 10 && reference_only(first_gt(fields)) {
        return Some(NoVote::Reference);
    }
    if cnvpytor {
        let info = fields[7];
        let ok = info_number(info, "pytorP1") < MAX_EVAL && info_number(info, "pytorQ0") < MAX_Q0;
        if !ok {
            return Some(NoVote::Quality);
        }
    }
    None
}

/// What one run read and kept.
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct Counts {
    pub records: usize,
    pub kept: usize,
    pub filter: usize,
    pub sample_filter: usize,
    pub reference: usize,
    pub quality: usize,
}

impl fmt::Display for Counts {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{} records, {} vote; no vote: {} FILTER not PASS, {} sample FT not PASS, {} genotyped reference",
            self.records, self.kept, self.filter, self.sample_filter, self.reference
        )?;
        if self.quality > 0 {
            write!(f, ", {} below CNVpytor's quality rule", self.quality)?;
        }
        Ok(())
    }
}

/// A failure that ends the run with exit 1.
#[derive(Debug)]
pub enum Error {
    Read(String, io::Error),
    Write(String, io::Error),
    Short(String, usize),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::Read(p, e) => write!(f, "cannot read {p}: {e}"),
            Error::Write(p, e) => write!(f, "cannot write {p}: {e}"),
            Error::Short(p, n) => write!(f, "{p} line {n}: a record needs at least 8 columns"),
        }
    }
}

/// Reads `input` (named `path`, for messages) and writes the records that
/// vote to `out`, with one sample column named `caller`.
pub fn prep(
    input: impl BufRead,
    path: &str,
    caller: &str,
    out: &mut impl Write,
) -> Result<Counts, Error> {
    let w = |e| Error::Write("the output".to_owned(), e);
    let cnvpytor = caller == "cnvpytor";
    let mut c = Counts::default();
    let mut sites = false;
    let mut has_gt = false;
    for (i, line) in input.lines().enumerate() {
        let line = line.map_err(|e| Error::Read(path.to_owned(), e))?;
        let line = line.trim_end_matches('\r');
        if line.is_empty() {
            continue;
        }
        if line.starts_with("##") {
            has_gt |= line.starts_with("##FORMAT=<ID=GT,");
            writeln!(out, "{line}").map_err(w)?;
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if line.starts_with("#CHROM") {
            sites = f.len() < 10;
            if sites && !has_gt {
                writeln!(out, "{GT_HEADER}").map_err(w)?;
            }
            let head = f.get(..8).unwrap_or(&f).join("\t");
            writeln!(out, "{head}\tFORMAT\t{caller}").map_err(w)?;
            continue;
        }
        if f.len() < 8 {
            return Err(Error::Short(path.to_owned(), i + 1));
        }
        c.records += 1;
        match verdict(&f, cnvpytor) {
            Some(NoVote::Filter) => c.filter += 1,
            Some(NoVote::SampleFilter) => c.sample_filter += 1,
            Some(NoVote::Reference) => c.reference += 1,
            Some(NoVote::Quality) => c.quality += 1,
            None => {
                c.kept += 1;
                let site = f[..8].join("\t");
                if sites || f.len() < 10 {
                    writeln!(out, "{site}\tGT\t./.").map_err(w)?;
                } else {
                    writeln!(out, "{site}\t{}\t{}", f[8], f[9]).map_err(w)?;
                }
            }
        }
    }
    out.flush().map_err(w)?;
    Ok(c)
}

/// Runs the verb on files: `out_path` `-` is stdout.
fn run(input: &str, out_path: &str, caller: &str, stdout: &mut dyn Write) -> Result<Counts, Error> {
    let reader = og_io::open(input).map_err(|e| Error::Read(input.to_owned(), e))?;
    if out_path == "-" {
        let mut w = BufWriter::new(stdout);
        return prep(reader, input, caller, &mut w);
    }
    // Written beside the target and renamed, so a failed run leaves no
    // half file for SURVIVOR to read.
    let part = format!("{out_path}.part");
    let file = File::create(&part).map_err(|e| Error::Write(part.clone(), e))?;
    let mut w = BufWriter::new(file);
    let c = prep(reader, input, caller, &mut w).inspect_err(|_| {
        let _ = std::fs::remove_file(&part);
    })?;
    drop(w);
    std::fs::rename(&part, out_path).map_err(|e| Error::Write(out_path.to_owned(), e))?;
    Ok(c)
}

/// `og sv survivor-prep --caller NAME IN OUT`: exit 0, 1 for an input or
/// output problem, 2 for a usage error.
pub fn cli(argv: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let p = match args::parse(&SPEC, argv) {
        Ok(p) => p,
        Err(ArgError::Help) => {
            let _ = write!(out, "{USAGE}");
            return 0;
        }
        Err(ArgError::Usage(m)) => {
            let _ = writeln!(err, "{USAGE}og sv survivor-prep: error: {m}");
            return 2;
        }
    };
    let caller = p.get("--caller").unwrap_or("");
    if caller.is_empty() || caller.bytes().any(|b| b.is_ascii_whitespace() || b == b':') {
        let _ = writeln!(
            err,
            "{USAGE}og sv survivor-prep: error: argument --caller: a caller name is one word, got {caller:?}"
        );
        return 2;
    }
    match run(&p.pos[0], &p.pos[1], caller, out) {
        Ok(c) => {
            let _ = writeln!(err, "og sv survivor-prep {caller}: {c}");
            0
        }
        Err(e) => {
            let _ = writeln!(err, "ERROR: {e}");
            1
        }
    }
}
