//! `og consumer`: step 40, a consumer raw-data file from the gVCF.
//!
//! Genealogy and health services that match relatives or read health traits
//! take the text file a genotyping chip company gives its customers. These
//! verbs write that file from the pipeline's own calls, for the sites of the
//! 23andMe v5 chip, so the genome can be used there. They write files and
//! nothing else: uploading one is the user's decision
//! (docs/40-export-consumer.md).
//!
//! - `og consumer build-sites --template T.txt.gz --dbsnp VCF|- --out
//!   SITES.tsv.gz --dropped DROPPED.tsv [--source TEXT]` ([`sites`]): the
//!   site table, the 23andMe v5 list joined by rsID with dbSNP's GRCh38
//!   VCF, which it reads itself, BGZF on many cores.
//! - `og consumer export --gvcf G.g.vcf.gz --sites SITES.tsv.gz --sample ID
//!   [--sex male|female] [--par PAR.bed] --out-dir DIR [--date YYYY-MM-DD]`
//!   ([`export`]): every site genotyped from the gVCF, in the 23andMe layout
//!   and the AncestryDNA layout, each zipped, plus a summary table.
//! - `og consumer check --layout 23andme|ancestrydna FILE` ([`check`]): the
//!   strict reader of both layouts, on a `.txt` or a `.zip`.
//!
//! Each verb writes what `bin/export_consumer.py` wrote, checked by
//! `scripts/ci/oracle-diff.sh` (tests/parity). The deliberate differences
//! are named where they are made.

pub mod check;
pub mod export;
pub mod sites;
pub mod zip;

use std::fmt::Write as _;

/// 23andMe's chromosome names, in the order its files list them.
pub const CHROMS: [&str; 25] = [
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17",
    "18", "19", "20", "21", "22", "X", "Y", "MT",
];
/// The index of X, Y and MT in [`CHROMS`].
pub const X: u8 = 22;
pub const Y: u8 = 23;
pub const MT: u8 = 24;

/// The place of a 23andMe chromosome name in [`CHROMS`].
pub fn chrom_index(name: &str) -> Option<u8> {
    CHROMS
        .iter()
        .position(|&c| c == name)
        .and_then(|i| u8::try_from(i).ok())
}

/// The GRCh38 contig of the pipeline for a chromosome of [`CHROMS`].
pub fn gvcf_contig(chrom: u8) -> String {
    if chrom == MT {
        "chrM".to_owned()
    } else {
        format!("chr{}", CHROMS[usize::from(chrom)])
    }
}

/// A verb that stops: `ERROR: <msg>` on stderr and the exit code.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Fail {
    pub code: u8,
    pub msg: String,
}

impl Fail {
    /// Exit 2: a refused input.
    pub fn input(msg: impl Into<String>) -> Fail {
        Fail {
            code: 2,
            msg: msg.into(),
        }
    }

    /// Exit 1: the data says no.
    pub fn data(msg: impl Into<String>) -> Fail {
        Fail {
            code: 1,
            msg: msg.into(),
        }
    }

    /// Exit 1 on a read or write error, naming the path.
    pub fn io(path: &str, e: &std::io::Error) -> Fail {
        Fail::data(format!("{path}: {e}"))
    }
}

/// Python's `repr()` of a str: single quotes unless the text has a `'` and no
/// `"`; backslash, the quote, tab, newline and carriage return escaped;
/// other control characters and the unprintable Latin-1 ones as `\xNN`.
/// Characters above U+00FF are kept as they are (U+FFFD, which an undecodable
/// byte reads as, is printable).
pub fn py_repr(s: &str) -> String {
    let quote = if s.contains('\'') && !s.contains('"') {
        '"'
    } else {
        '\''
    };
    let mut out = String::with_capacity(s.len() + 2);
    out.push(quote);
    for c in s.chars() {
        match c {
            '\\' => out.push_str("\\\\"),
            '\t' => out.push_str("\\t"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            c if c == quote => {
                out.push('\\');
                out.push(c);
            }
            '\0'..='\x1f' | '\x7f'..='\u{a0}' | '\u{ad}' => {
                let _ = write!(out, "\\x{:02x}", u32::from(c));
            }
            c => out.push(c),
        }
    }
    out.push(quote);
    out
}

/// Python's `repr()` of a list of str: `['a', 'b']`.
pub fn py_repr_list<S: AsRef<str>>(items: &[S]) -> String {
    let parts: Vec<String> = items.iter().map(|s| py_repr(s.as_ref())).collect();
    format!("[{}]", parts.join(", "))
}

/// `s[:n]`, counted in characters.
pub fn head(s: &str, n: usize) -> &str {
    s.char_indices().nth(n).map_or(s, |(i, _)| &s[..i])
}

/// `os.path.join(dir, name)` for a relative or absolute `name`.
pub fn path_join(dir: &str, name: &str) -> String {
    if name.starts_with('/') || dir.is_empty() {
        name.to_owned()
    } else if dir.ends_with('/') {
        format!("{dir}{name}")
    } else {
        format!("{dir}/{name}")
    }
}

/// `os.path.basename(path)`.
pub fn basename(path: &str) -> &str {
    path.rsplit('/').next().unwrap_or(path)
}

/// The options of one verb, read as argparse reads them: `--name value` or
/// `--name=value` in any order (the last one wins), and positionals.
pub struct Args {
    opts: Vec<(&'static str, String)>,
    pub positional: Vec<String>,
}

impl Args {
    /// `names` are the options the verb takes, each taking a value.
    pub fn parse(args: &[String], names: &[&'static str]) -> Result<Args, String> {
        let mut opts = Vec::new();
        let mut positional = Vec::new();
        let mut it = args.iter();
        while let Some(a) = it.next() {
            if let Some(body) = a.strip_prefix("--") {
                let (name, inline) = match body.split_once('=') {
                    Some((n, v)) => (n, Some(v.to_owned())),
                    None => (body, None),
                };
                let Some(&known) = names.iter().find(|&&n| n == name) else {
                    return Err(format!("unrecognized arguments: {a}"));
                };
                let value = match inline {
                    Some(v) => v,
                    None => it
                        .next()
                        .cloned()
                        .ok_or_else(|| format!("argument --{known}: expected one argument"))?,
                };
                opts.push((known, value));
            } else {
                positional.push(a.clone());
            }
        }
        Ok(Args { opts, positional })
    }

    /// The last value given for `name`.
    pub fn get(&self, name: &str) -> Option<&str> {
        self.opts
            .iter()
            .rev()
            .find(|(n, _)| *n == name)
            .map(|(_, v)| v.as_str())
    }

    /// The value of a required option.
    pub fn need(&self, name: &str) -> Result<&str, String> {
        self.get(name)
            .ok_or_else(|| format!("the following arguments are required: --{name}"))
    }
}

const USAGE: &str = "usage: og consumer build-sites --template T --dbsnp VCF|- --out SITES.tsv.gz --dropped DROPPED.tsv [--source TEXT]
       og consumer export --gvcf G.g.vcf.gz --sites SITES.tsv.gz --sample ID [--sex male|female] [--par PAR.bed] --out-dir DIR [--date YYYY-MM-DD]
       og consumer check --layout 23andme|ancestrydna FILE
";

/// A wrong command line: the usage and the error on stderr, exit 2.
fn usage_error(verb: &str, msg: &str) -> u8 {
    eprint!("{USAGE}");
    eprintln!("og consumer{}: error: {msg}", verb);
    2
}

/// `og consumer <verb> [args]`: the exit code.
pub fn cli(args: &[String]) -> u8 {
    let Some(verb) = args.first() else {
        return usage_error("", "the following arguments are required: verb");
    };
    let rest = &args[1..];
    let result = match verb.as_str() {
        "build-sites" => sites::cli(rest),
        "export" => export::cli(rest),
        "check" => check::cli(rest),
        "-h" | "--help" => {
            print!("{USAGE}");
            return 0;
        }
        other => return usage_error("", &format!("unknown verb '{other}'")),
    };
    match result {
        Ok(code) => code,
        Err(Usage(msg)) => usage_error(&format!(" {verb}"), &msg),
    }
}

/// A command line the verb refuses.
pub struct Usage(pub String);

impl From<String> for Usage {
    fn from(s: String) -> Usage {
        Usage(s)
    }
}

/// Prints a [`Fail`] as Python's `die()` did and returns its exit code.
pub fn report(f: &Fail) -> u8 {
    eprintln!("ERROR: {}", f.msg);
    f.code
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn repr_matches_python() {
        assert_eq!(py_repr("a\tb"), "'a\\tb'");
        assert_eq!(py_repr("it's"), "\"it's\"");
        assert_eq!(py_repr("'\""), "'\\'\"'");
        assert_eq!(py_repr("x\r"), "'x\\r'");
        assert_eq!(py_repr("\u{1}\u{7f}\u{fffd}"), "'\\x01\\x7f\u{fffd}'");
        assert_eq!(py_repr_list(&["rsid", "chrom"]), "['rsid', 'chrom']");
    }

    #[test]
    fn head_counts_characters() {
        assert_eq!(head("ñandú", 3), "ñan");
        assert_eq!(head("ab", 80), "ab");
    }

    #[test]
    fn joins_as_os_path_does() {
        assert_eq!(path_join(".", "S"), "./S");
        assert_eq!(path_join("out/", "S"), "out/S");
        assert_eq!(path_join("out", "/abs"), "/abs");
        assert_eq!(path_join("", "S"), "S");
        assert_eq!(basename("./a/S_x.zip"), "S_x.zip");
    }

    #[test]
    fn args_as_argparse_reads_them() {
        let a: Vec<String> = ["--sex", "male", "--out-dir=o", "f", "--sex", ""]
            .iter()
            .map(|s| (*s).to_owned())
            .collect();
        let p = Args::parse(&a, &["sex", "out-dir"]).unwrap();
        assert_eq!(p.get("sex"), Some(""));
        assert_eq!(p.get("out-dir"), Some("o"));
        assert_eq!(p.positional, ["f"]);
        assert!(Args::parse(&a, &["sex"]).is_err());
        assert!(p.need("gvcf").is_err());
    }
}
