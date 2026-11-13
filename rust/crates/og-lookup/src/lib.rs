//! The genotype lookup: "what is my genotype at rs4988235, at
//! chr2:135851076, or in CFTR?", answered from the files the pipeline
//! already wrote (`docs/lookup.md`).
//!
//! - `og lookup query --gvcf G [--annotated A] [--summary S] [--sites T]
//!   [--table T]... [--genes GTF] Q...` answers every query of the command
//!   line in one process. `scripts/owngenome-lookup.sh` runs it.
//! - `og lookup resolve` and `og lookup report` are its two halves, with
//!   the command lines `bin/genome_lookup.py` had.
//!
//! A gVCF or annotated VCF with a `.tbi` or `.csi` beside it is read through
//! the index, region by region ([`og_io::tabix`]); one without is streamed
//! whole, keeping only the region's records. An rsID the tables do not
//! place is looked for in the annotated VCF in one pass for every such
//! rsID of the command line.
//!
//! The output is what the Python wrote, byte for byte (`tests/parity`),
//! with one deliberate change: a gene's "bases have a genotype call" count
//! is the union of the called records, so overlapping records (a deletion
//! and the reference block that follows it) no longer count a base twice.

mod query;
mod report;
mod resolve;
mod summary;
mod vcf;

use std::io::{self, Write};
use vcf::{Region, Span};

/// Why a query or a command stopped: the message after `ERROR: ` and the
/// exit code (1 an unreadable file, 2 a bad query or value, 3 not found).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Error {
    pub msg: String,
    pub code: u8,
}

impl Error {
    fn value(msg: impl Into<String>) -> Error {
        Error {
            msg: msg.into(),
            code: 2,
        }
    }
    fn not_found(msg: impl Into<String>) -> Error {
        Error {
            msg: msg.into(),
            code: 3,
        }
    }
    fn other(msg: impl Into<String>) -> Error {
        Error {
            msg: msg.into(),
            code: 1,
        }
    }
    fn io(e: &io::Error, path: &str) -> Error {
        Error::other(format!("{path}: {e}"))
    }
}

/// Python's `repr()` of a string.
fn repr(s: &str) -> String {
    pycompat::value::py_repr(&pycompat::json::Value::Str(s.to_owned())).unwrap_or_default()
}

/// The files a lookup reads.
#[derive(Debug, Default)]
pub struct Inputs {
    pub gvcf: Option<String>,
    pub annotated: Option<String>,
    pub summary: Option<String>,
    pub sites: Option<String>,
    pub tables: Vec<String>,
    pub genes: Option<String>,
}

const USAGE: &str = "usage: og lookup query --gvcf GVCF [--annotated VCF] [--summary JSON] [--sites TSV] [--table TSV]... [--genes GTF] [--] QUERY...
       og lookup resolve --query Q [--sites TSV] [--table TSV]... [--genes GTF] [--annotated VCF]
       og lookup report --query Q --kind site|gene --chrom C --start N --end N [--ref BASE] [--source TEXT] --gvcf GVCF [--annotated VCF] [--summary JSON]
";

/// A command line read as `--name VALUE` / `--name=VALUE` options (the
/// repeatable ones gathered) and positionals.
struct Args {
    opts: Vec<(String, String)>,
    pos: Vec<String>,
}

impl Args {
    fn parse(args: &[String], known: &[&str], positionals: bool) -> Result<Args, String> {
        let mut a = Args {
            opts: Vec::new(),
            pos: Vec::new(),
        };
        let mut it = args.iter();
        while let Some(x) = it.next() {
            if x == "--" && positionals {
                a.pos.extend(it.by_ref().cloned());
                break;
            }
            if x.starts_with("--") {
                let (name, val) = match x.split_once('=') {
                    Some((n, v)) => (n.to_owned(), Some(v.to_owned())),
                    None => (x.clone(), None),
                };
                if !known.contains(&name.as_str()) {
                    return Err(format!("unrecognized arguments: {x}"));
                }
                let val = match val {
                    Some(v) => v,
                    None => it
                        .next()
                        .cloned()
                        .ok_or_else(|| format!("argument {name}: expected one argument"))?,
                };
                a.opts.push((name, val));
            } else if positionals {
                a.pos.push(x.clone());
            } else {
                return Err(format!("unrecognized arguments: {x}"));
            }
        }
        Ok(a)
    }

    /// The last value of `name`.
    fn get(&self, name: &str) -> Option<&str> {
        self.opts
            .iter()
            .rev()
            .find(|(n, _)| n == name)
            .map(|(_, v)| v.as_str())
    }

    /// The last value of `name`, `None` when absent or empty.
    fn path(&self, name: &str) -> Option<String> {
        self.get(name).filter(|v| !v.is_empty()).map(str::to_owned)
    }

    fn all(&self, name: &str) -> Vec<String> {
        self.opts
            .iter()
            .filter(|(n, _)| n == name)
            .map(|(_, v)| v.clone())
            .collect()
    }

    fn required(&self, names: &[&str]) -> Result<(), String> {
        let missing: Vec<&str> = names
            .iter()
            .copied()
            .filter(|n| self.get(n).is_none())
            .collect();
        if missing.is_empty() {
            Ok(())
        } else {
            Err(format!(
                "the following arguments are required: {}",
                missing.join(", ")
            ))
        }
    }

    fn inputs(&self) -> Inputs {
        Inputs {
            gvcf: self.path("--gvcf"),
            annotated: self.path("--annotated"),
            summary: self.path("--summary"),
            sites: self.path("--sites"),
            tables: self.all("--table"),
            genes: self.path("--genes"),
        }
    }
}

/// `og lookup VERB ARGS`: the exit code.
pub fn cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let usage = |err: &mut dyn Write, msg: &str| {
        let _ = writeln!(err, "{USAGE}og lookup: error: {msg}");
        2
    };
    let (verb, rest) = match args.split_first() {
        Some((v, r)) => (v.as_str(), r),
        None => return usage(err, "the following arguments are required: verb"),
    };
    if rest.iter().any(|a| a == "-h" || a == "--help") || verb == "-h" || verb == "--help" {
        let _ = out.write_all(USAGE.as_bytes());
        return 0;
    }
    let parsed = match verb {
        "query" => Args::parse(
            rest,
            &[
                "--gvcf",
                "--annotated",
                "--summary",
                "--sites",
                "--table",
                "--genes",
            ],
            true,
        )
        .and_then(|a| {
            a.required(&["--gvcf"])?;
            if a.pos.is_empty() {
                return Err("the following arguments are required: QUERY".into());
            }
            Ok(a)
        }),
        "resolve" => Args::parse(
            rest,
            &["--query", "--sites", "--table", "--genes", "--annotated"],
            false,
        )
        .and_then(|a| a.required(&["--query"]).map(|()| a)),
        "report" => Args::parse(
            rest,
            &[
                "--query",
                "--kind",
                "--chrom",
                "--start",
                "--end",
                "--ref",
                "--source",
                "--gvcf",
                "--annotated",
                "--summary",
            ],
            false,
        )
        .and_then(|a| {
            a.required(&["--query", "--kind", "--chrom", "--start", "--end", "--gvcf"])?;
            match a.get("--kind") {
                Some("site" | "gene") => Ok(a),
                Some(k) => Err(format!(
                    "argument --kind: invalid choice: {} (choose from 'site', 'gene')",
                    repr(k)
                )),
                None => unreachable!("--kind is required"),
            }
        }),
        other => {
            return usage(
                err,
                &format!(
                    "unknown lookup verb {} (choose from query, resolve, report)",
                    repr(other)
                ),
            );
        }
    };
    let a = match parsed {
        Ok(a) => a,
        Err(msg) => return usage(err, &msg),
    };
    match verb {
        "query" => run_query(&a, out, err),
        "resolve" => run_resolve(&a, out, err),
        _ => run_report(&a, out, err),
    }
}

fn fail(err: &mut dyn Write, e: &Error) -> u8 {
    let _ = writeln!(err, "ERROR: {}", e.msg);
    e.code
}

fn run_resolve(a: &Args, out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let q = a.get("--query").unwrap_or_default().to_owned();
    match resolve::resolve_all(&[q], &a.inputs()).remove(0) {
        Ok(rows) => {
            let text: String = rows.iter().map(resolve::Row::line).collect();
            match out.write_all(text.as_bytes()) {
                Ok(()) => 0,
                Err(e) => fail(err, &Error::other(format!("cannot write stdout: {e}"))),
            }
        }
        Err(e) => fail(err, &e),
    }
}

/// `--start`/`--end` as argparse's `type=int` reads them.
fn position(a: &Args, name: &str) -> Result<u64, String> {
    let v = a.get(name).unwrap_or_default();
    pycompat::num::parse_int(v)
        .ok()
        .and_then(|n| u64::try_from(n).ok())
        .ok_or_else(|| format!("argument {name}: invalid int value: {}", repr(v)))
}

fn run_report(a: &Args, out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let (start, end) = match (position(a, "--start"), position(a, "--end")) {
        (Ok(s), Ok(e)) => (s, e),
        (Err(m), _) | (_, Err(m)) => {
            let _ = writeln!(err, "{USAGE}og lookup: error: {m}");
            return 2;
        }
    };
    let inp = a.inputs();
    let row = resolve::Row {
        kind: if a.get("--kind") == Some("site") {
            "site"
        } else {
            "gene"
        },
        chrom: a.get("--chrom").unwrap_or_default().to_owned(),
        start,
        end,
        ref_: a.get("--ref").unwrap_or(".").to_owned(),
        source: a.get("--source").unwrap_or_default().to_owned(),
    };
    let summary = inp.summary.as_deref().map(summary::load);
    let query = a.get("--query").unwrap_or_default();
    match answer(query, &row, &inp, summary.as_ref(), &mut io::sink()) {
        Ok(text) => write_out(out, err, &text, 0),
        Err((text, e)) => {
            let rc = write_out(out, err, &text, 0);
            if rc != 0 { rc } else { fail(err, &e) }
        }
    }
}

fn write_out(out: &mut dyn Write, err: &mut dyn Write, text: &str, rc: u8) -> u8 {
    match out.write_all(text.as_bytes()).and_then(|()| out.flush()) {
        Ok(()) => rc,
        Err(e) => fail(err, &Error::other(format!("cannot write stdout: {e}"))),
    }
}

/// The answer for one region, or what was written before the error and
/// the error.
fn answer(
    query: &str,
    row: &resolve::Row,
    inp: &Inputs,
    summary: Option<&summary::Summary>,
    note: &mut dyn Write,
) -> Result<String, (String, Error)> {
    let chrom = query::norm_chrom(&row.chrom);
    let span = Span {
        chrom: &chrom,
        start: row.start,
        end: row.end,
    };
    let gvcf_path = inp.gvcf.as_deref().unwrap_or_default();
    let read = |p: &str, note: &mut dyn Write| vcf::read(p, span, note);
    let gvcf = read(gvcf_path, note).map_err(|e| (String::new(), e))?;
    let ann: Option<Region> = match &inp.annotated {
        Some(p) => Some(read(p, note).map_err(|e| (String::new(), e))?),
        None => None,
    };
    let ask = report::Ask {
        query,
        row,
        gvcf: &gvcf,
        ann: ann.as_ref(),
        summary,
    };
    let mut text = String::new();
    match report::write(&ask, &mut text) {
        Ok(()) => Ok(text),
        Err(e) => Err((text, e)),
    }
}

/// Every query in turn: its answer and a blank line, or its error on
/// stderr and a blank line. The exit code is the last failing query's.
fn run_query(a: &Args, out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let inp = a.inputs();
    let summary = inp.summary.as_deref().map(summary::load);
    let mut rc = 0;
    for (q, res) in a.pos.iter().zip(resolve::resolve_all(&a.pos, &inp)) {
        match res {
            Err(e) => rc = fail(err, &e),
            Ok(rows) => {
                for row in rows {
                    let (text, e) = match answer(q, &row, &inp, summary.as_ref(), err) {
                        Ok(t) => (t, None),
                        Err((t, e)) => (t, Some(e)),
                    };
                    let w = write_out(out, err, &text, 0);
                    if w != 0 {
                        return w;
                    }
                    if let Some(e) = e {
                        rc = fail(err, &e);
                    }
                    let w = write_out(out, err, "\n", 0);
                    if w != 0 {
                        return w;
                    }
                }
                continue;
            }
        }
        let w = write_out(out, err, "\n", 0);
        if w != 0 {
            return w;
        }
    }
    rc
}
