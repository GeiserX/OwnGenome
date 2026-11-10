//! The command lines of `og clinvar`.
//!
//! ```text
//! og clinvar join-adjacent --mode prepare|finish|tag|build [--fasta REF] [--window N] [VCF]
//! og clinvar hits [HITS_VCF]
//! og clinvar deep-intronic [--no-gvcf] [--min-gq N] SITES VCF_CALLS GVCF_CALLS
//! ```
//!
//! A missing file, or `-`, is standard input; plain, gzip and BGZF input are
//! read alike. A refused command line is exit 2, a file that cannot be read
//! or written exit 1.

use std::fmt::Display;
use std::io::{self, BufRead, BufWriter, Write};

use og_intervals::fasta::Indexed;

use crate::join::{self, Mode, WINDOW};
use crate::{hits, sites};

const USAGE: &str = "usage: og clinvar join-adjacent --mode prepare|finish|tag|build [--fasta REF] [--window N] [VCF]
       og clinvar hits [HITS_VCF]
       og clinvar deep-intronic [--no-gvcf] [--min-gq N] SITES VCF_CALLS GVCF_CALLS";

fn refuse(msg: &str) -> u8 {
    eprintln!("{USAGE}");
    eprintln!("og clinvar: error: {msg}");
    2
}

fn fail(what: &str, e: &dyn Display) -> u8 {
    eprintln!("og clinvar {what}: {e}");
    1
}

/// Options with a value, flags, and the other words in order.
struct Line {
    values: Vec<(&'static str, String)>,
    flags: Vec<&'static str>,
    files: Vec<String>,
}

impl Line {
    fn parse(
        args: &[String],
        valued: &[&'static str],
        flags: &[&'static str],
    ) -> Result<Line, String> {
        let mut l = Line {
            values: Vec::new(),
            flags: Vec::new(),
            files: Vec::new(),
        };
        let mut it = args.iter();
        while let Some(w) = it.next() {
            let (name, inline) = match w.split_once('=') {
                Some((n, v)) if n.starts_with("--") => (n, Some(v.to_owned())),
                _ => (w.as_str(), None),
            };
            if let Some(&k) = valued.iter().find(|&&k| k == name) {
                let v = match inline {
                    Some(v) => v,
                    None => it
                        .next()
                        .cloned()
                        .ok_or_else(|| format!("argument {k}: expected one argument"))?,
                };
                l.values.retain(|(n, _)| *n != k);
                l.values.push((k, v));
            } else if let Some(&k) = flags.iter().find(|&&k| k == name && inline.is_none()) {
                l.flags.push(k);
            } else if w.starts_with("--") {
                return Err(format!("unrecognized arguments: {w}"));
            } else {
                l.files.push(w.clone());
            }
        }
        Ok(l)
    }

    fn get(&self, k: &str) -> Option<&str> {
        self.values
            .iter()
            .find(|(n, _)| *n == k)
            .map(|(_, v)| v.as_str())
    }

    /// A whole-number option, or its default.
    fn number(&self, k: &str, default: u64) -> Result<u64, String> {
        match self.get(k) {
            None => Ok(default),
            Some(v) => v
                .parse()
                .ok()
                .filter(|_| v.bytes().all(|c| c.is_ascii_digit()))
                .ok_or_else(|| format!("argument {k}: '{v}' is not a whole number")),
        }
    }
}

/// Standard input for `None` or `-`, else the file.
fn input(path: Option<&str>) -> io::Result<Box<dyn BufRead>> {
    match path {
        None | Some("-") => og_io::reader(io::stdin()),
        Some(p) => og_io::open(p),
    }
}

/// Runs `og clinvar <verb> ...` and returns the exit code.
pub fn main(args: &[String]) -> u8 {
    let Some(verb) = args.first() else {
        return refuse("the following arguments are required: verb");
    };
    let rest = &args[1..];
    let stdout = io::stdout();
    let mut out = BufWriter::with_capacity(1 << 16, stdout.lock());
    let code = match verb.as_str() {
        "join-adjacent" => join_adjacent(rest, &mut out),
        "hits" => match Line::parse(rest, &[], &[]) {
            Err(e) => return refuse(&e),
            Ok(l) if l.files.len() > 1 => return refuse("og clinvar hits reads one hits VCF"),
            Ok(l) => match input(l.files.first().map(String::as_str)) {
                Err(e) => fail("hits", &e),
                Ok(i) => hits::hits(i, &mut out).map_or_else(|e| fail("hits", &e), |()| 0),
            },
        },
        "deep-intronic" => deep_intronic(rest, &mut out),
        other => return refuse(&format!("unknown clinvar verb '{other}'")),
    };
    match out.flush() {
        Err(e) if code == 0 => fail(verb, &e),
        _ => code,
    }
}

fn join_adjacent(args: &[String], out: &mut impl Write) -> u8 {
    let l = match Line::parse(args, &["--mode", "--fasta", "--window"], &[]) {
        Ok(l) => l,
        Err(e) => return refuse(&e),
    };
    let window = match l.number("--window", WINDOW) {
        Ok(w) => w,
        Err(e) => return refuse(&e),
    };
    if l.files.len() > 1 {
        return refuse("join-adjacent reads one VCF");
    }
    let mode = l.get("--mode");
    match (mode, l.get("--fasta")) {
        (Some("build"), None) => return refuse("--mode build needs --fasta"),
        (Some("prepare" | "finish" | "tag"), Some(_)) => {
            return refuse("--fasta goes with --mode build only");
        }
        (None, _) => return refuse("the following arguments are required: --mode"),
        (Some(m), _) if !["prepare", "finish", "tag", "build"].contains(&m) => {
            return refuse(&format!(
                "argument --mode: invalid choice: '{m}' (choose from prepare, finish, tag, build)"
            ));
        }
        _ => {}
    }
    let what = "join-adjacent";
    let input = match input(l.files.first().map(String::as_str)) {
        Ok(i) => i,
        Err(e) => return fail(what, &e),
    };
    let done = match (mode, l.get("--fasta")) {
        (Some("prepare"), _) => {
            join::join(Mode::Prepare, window, input, out).map_err(|e| e.to_string())
        }
        (Some("finish"), _) => join::finish(input, out).map_err(|e| e.to_string()),
        (Some("tag"), _) => join::tag(input, window, out).map_err(|e| e.to_string()),
        (Some("build"), Some(fasta)) => match Indexed::open(fasta) {
            Ok(mut ix) => {
                join::join(Mode::Build(&mut ix), window, input, out).map_err(|e| e.to_string())
            }
            Err(e) => Err(format!("{fasta}: {e}")),
        },
        _ => unreachable!("the mode was checked above"),
    };
    done.map_or_else(|e| fail(what, &e), |()| 0)
}

fn deep_intronic(args: &[String], out: &mut impl Write) -> u8 {
    let l = match Line::parse(args, &["--min-gq"], &["--no-gvcf"]) {
        Ok(l) => l,
        Err(e) => return refuse(&e),
    };
    let min_gq = match l.number("--min-gq", 30) {
        Ok(v) => v,
        Err(e) => return refuse(&e),
    };
    let [s, v, g] = l.files.as_slice() else {
        return refuse("deep-intronic takes SITES VCF_CALLS GVCF_CALLS");
    };
    let mut texts = Vec::new();
    for p in [s, v, g] {
        let mut buf = Vec::new();
        match input(Some(p)).and_then(|mut i| i.read_to_end(&mut buf)) {
            Ok(_) => texts.push(buf),
            Err(e) => return fail("deep-intronic", &format!("{p}: {e}")),
        }
    }
    let o = sites::Options {
        min_gq: l.get("--min-gq").unwrap_or("30").to_owned(),
        min_gq_value: min_gq,
        gvcf: !l.flags.contains(&"--no-gvcf"),
    };
    sites::deep_intronic(&texts[0], &texts[1], &texts[2], &o, out)
        .map_or_else(|e| fail("deep-intronic", &e), |()| 0)
}
