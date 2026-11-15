//! `og sv add-callers`: the VCF AnnotSV classifies, with the calls of the
//! other SV callers that no call in it covers (port of
//! `bin/cnvpytor_to_vcf.py`).
//!
//! AnnotSV classifies only the records of the one VCF it is given, and
//! step 05 gives it Manta's calls. A large deletion that CNVpytor (read
//! depth) or Delly calls, but Manta calls as something else or not at all,
//! would never get an ACMG class. This writes the base VCF again with,
//! after it, the deletions and duplications of the other callers that no
//! call already in the file covers:
//!
//! - Delly (step 19): FILTER PASS, SVTYPE DEL or DUP, a genotype the sample
//!   carries, and for a call over [`DEPTH_MIN_SIZE`] bp a read depth that
//!   shows the change ([`depth_contradicts`]);
//! - CNVpytor (step 18): DEL and DUP with e-val (`pytorP1`) below
//!   `--max-eval` and a share of MAPQ 0 reads (`pytorQ0`) below `--max-q0`.
//!
//! A call is covered when a record of the same type on the same contig
//! overlaps it by at least `--min-overlap` of the longer of the two. Only a
//! base record that is a call the sample carries ([`call_holds`]) covers.
//! Delly's calls go in before CNVpytor's, so a deletion both make is added
//! once. Each record names its caller in INFO SVCALLER. The output is sorted
//! by the base header's `##contig` order, then POS, keeping the input order
//! of ties.
//!
//! The base records pass through byte for byte but for INFO, which gets
//! `;SVCALLER=<name>` (or is replaced, when it is `.` or empty). Fixes
//! against the Python, each with a golden in `tests/parity`: inputs are
//! read as gzip by their magic bytes, not their suffix; a base record whose
//! POS is not a whole number, a record of fewer than 8 columns and a
//! missing file are an `ERROR:` line and exit 1, not a traceback.

use std::fmt;
use std::io::{self, Write};

use og_filters::TextLines;
use og_filters::args::{self, ArgError, Spec};
use og_pedigree::Par;

use super::{
    DEPTH_MIN_SIZE, OverlapIndex, Rdcn, Sex, SvKind, call_holds, depth_contradicts, first_gt,
};

const USAGE: &str = "usage: og sv add-callers --base BASE.vcf[.gz]|- [--base-name NAME]
       [--delly DELLY.vcf.gz] [--sex male|female] [--par-bed PAR.bed]
       [--cnvpytor CNVPYTOR.vcf.gz]
       [--max-eval E] [--max-q0 Q] [--min-overlap F] > OUT.vcf
";

const SPEC: Spec = Spec {
    prog: "og sv add-callers",
    usage: USAGE,
    options: &[
        "--base",
        "--base-name",
        "--delly",
        "--sex",
        "--par-bed",
        "--cnvpytor",
        "--max-eval",
        "--max-q0",
        "--min-overlap",
    ],
    positionals: &[],
    required: &["--base"],
};

/// The header lines a record of the output may need, in the order they go
/// in when the base lacks them.
const STANDARD_INFO: [(&str, &str); 7] = [
    (
        "END",
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">",
    ),
    (
        "SVTYPE",
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">",
    ),
    (
        "SVLEN",
        "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">",
    ),
    (
        "IMPRECISE",
        "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">",
    ),
    (
        "CIPOS",
        "##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS\">",
    ),
    (
        "CIEND",
        "##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END\">",
    ),
    (
        "SVCALLER",
        "##INFO=<ID=SVCALLER,Number=1,Type=String,Description=\"The SV caller of this record: the input of step 05 (manta), or delly or cnvpytor for a call no record of that input covers\">",
    ),
];

/// The INFO keys an added record keeps, in this order.
const KEPT: [&str; 4] = ["END", "SVTYPE", "CIPOS", "CIEND"];

/// What the command line asks for.
#[derive(Debug, Clone)]
pub struct Opts {
    pub base: String,
    pub base_name: String,
    pub delly: Option<String>,
    pub sex: Sex,
    pub par_bed: Option<String>,
    pub cnvpytor: Option<String>,
    pub max_eval: f64,
    pub max_q0: f64,
    pub min_overlap: f64,
}

/// Why the step stopped.
#[derive(Debug)]
pub enum Error {
    /// A file that cannot be opened or read.
    Read(String, io::Error),
    /// The PAR file could not be read.
    Par(String),
    /// A record of fewer than 8 columns: file and record number.
    Short(String, u64),
    /// A base record whose POS is not a whole number: file, record number, POS.
    BadPos(String, u64, String),
    /// The output could not be written.
    Write(io::Error),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::Read(p, e) => write!(f, "cannot read {p}: {}", og_filters::os_error(e, p)),
            Error::Par(m) => write!(f, "{m}"),
            Error::Short(p, n) => write!(
                f,
                "{p}: record {n} has fewer than the 8 columns of a VCF record"
            ),
            Error::BadPos(p, n, pos) => write!(
                f,
                "{p}: record {n}: POS {} is not a whole number",
                og_filters::py_repr(pos)
            ),
            Error::Write(e) => write!(f, "cannot write the output: {e}"),
        }
    }
}

/// INFO as a dict of its items: the last item of a key wins; a flag has no
/// value.
struct Info<'a>(Vec<(&'a str, Option<&'a str>)>);

impl<'a> Info<'a> {
    fn parse(info: &'a str) -> Info<'a> {
        if info.is_empty() || info == "." {
            return Info(Vec::new());
        }
        Info(
            info.split(';')
                .map(|it| match it.split_once('=') {
                    Some((k, v)) => (k, Some(v)),
                    None => (it, None),
                })
                .collect(),
        )
    }

    /// `None`: no such key; `Some(None)`: a flag.
    fn get(&self, key: &str) -> Option<Option<&'a str>> {
        self.0
            .iter()
            .rev()
            .find(|(k, _)| *k == key)
            .map(|&(_, v)| v)
    }

    /// The value as Python's `str()` of the dict entry reads it: a flag is
    /// `True`.
    fn text(&self, key: &str) -> &'a str {
        match self.get(key) {
            Some(Some(v)) => v,
            Some(None) => "True",
            None => "",
        }
    }

    /// `float()` of the entry: a flag is 1, a missing or unreadable value
    /// NaN, which no threshold passes.
    fn float(&self, key: &str) -> f64 {
        match self.get(key) {
            Some(Some(v)) => pycompat::num::parse_float(v).unwrap_or(f64::NAN),
            Some(None) => 1.0,
            None => f64::NAN,
        }
    }
}

fn int(s: &str) -> Option<i64> {
    pycompat::num::parse_int(s)
        .ok()
        .and_then(|v| i64::try_from(v).ok())
}

/// The type of a record: SVTYPE, or its symbolic ALT.
fn sv_kind(fields: &[&str], info: &Info) -> SvKind {
    let t = info.text("SVTYPE");
    if t.starts_with("DUP") || fields[4].starts_with("<DUP") {
        SvKind::Dup
    } else if t == "DEL" || fields[4] == "<DEL>" {
        SvKind::Del
    } else {
        SvKind::Other
    }
}

/// POS and END of a record; `None` when either is not a whole number.
fn span(fields: &[&str], info: &Info) -> Option<(i64, i64)> {
    let end = match info.get("END")? {
        Some(v) => int(v)?,
        // int(True) in the Python.
        None => 1,
    };
    Some((int(fields[1])?, end))
}

/// The records of a VCF, one call of `each` per record with its number
/// (blank lines skipped); the header lines go to `header` when given.
fn read_vcf(
    path: &str,
    mut header: Option<&mut Vec<String>>,
    mut each: impl FnMut(&str, u64) -> Result<(), Error>,
) -> Result<(), Error> {
    let read = |e| Error::Read(path.to_owned(), e);
    let input = og_io::open(path).map_err(read)?;
    let mut lines = if path == "-" {
        TextLines::stdin(input)
    } else {
        TextLines::file(input)
    };
    let mut line = String::new();
    let mut n = 0u64;
    while lines.read_line(&mut line).map_err(read)? {
        let body = line.trim_end_matches('\n');
        if body.is_empty() {
            continue;
        }
        if body.starts_with('#') {
            if let Some(h) = header.as_deref_mut() {
                h.push(body.to_owned());
            }
            continue;
        }
        n += 1;
        each(body, n)?;
    }
    Ok(())
}

/// The fields of a record, at least 8.
fn fields<'a>(line: &'a str, path: &str, n: u64) -> Result<Vec<&'a str>, Error> {
    let f: Vec<&str> = line.split('\t').collect();
    if f.len() < 8 {
        return Err(Error::Short(path.to_owned(), n));
    }
    Ok(f)
}

/// One record of the output: its line and where it sorts.
struct Out {
    line: String,
    /// The contig, `line[..chrom]`.
    chrom: usize,
    /// The base record's number, for a POS error; 0 for an added record.
    base_n: u64,
}

/// Whether a record of another caller is one to add: its fields, INFO,
/// type and span.
type Keep<'k> = dyn FnMut(&[&str], &Info, SvKind, (i64, i64)) -> bool + 'k;

/// Counts of one other caller, for the summary.
#[derive(Default)]
struct Taken {
    kept: usize,
    added: usize,
}

/// Runs the step: the VCF to `out`, and the summary line it returns.
pub fn run(o: &Opts, out: impl Write) -> Result<String, Error> {
    let mut header = Vec::new();
    let mut records: Vec<Out> = Vec::new();
    let mut known = OverlapIndex::default();
    let mut base_lines = Vec::new();
    read_vcf(&o.base, Some(&mut header), |line, n| {
        fields(line, &o.base, n)?;
        base_lines.push((line.to_owned(), n));
        Ok(())
    })?;
    let n_samples = header
        .iter()
        .rfind(|h| h.starts_with("#CHROM"))
        .map_or(0, |h| h.split('\t').count().saturating_sub(9));
    for (line, n) in base_lines {
        let f: Vec<&str> = line.split('\t').collect();
        let info = Info::parse(f[7]);
        let kind = sv_kind(&f, &info);
        if kind != SvKind::Other
            && let Some((s, e)) = span(&f, &info)
            && call_holds(&f, n_samples)
        {
            known.add(f[0], kind, s, e);
        }
        let new_info = if f[7].is_empty() || f[7] == "." {
            format!("SVCALLER={}", o.base_name)
        } else {
            format!("{};SVCALLER={}", f[7], o.base_name)
        };
        let mut g = f.clone();
        g[7] = &new_info;
        let joined = g.join("\t");
        records.push(Out {
            chrom: f[0].len(),
            line: joined,
            base_n: n,
        });
    }

    let mut summary = Vec::new();
    let mut added_any = false;
    let mut take = |caller: &str,
                    path: &str,
                    keep: &mut Keep,
                    records: &mut Vec<Out>|
     -> Result<Taken, Error> {
        let mut t = Taken::default();
        read_vcf(path, None, |line, n| {
            let f = fields(line, path, n)?;
            let info = Info::parse(f[7]);
            let kind = sv_kind(&f, &info);
            if kind == SvKind::Other {
                return Ok(());
            }
            let Some(sp) = span(&f, &info) else {
                return Ok(());
            };
            if !keep(&f, &info, kind, sp) {
                return Ok(());
            }
            t.kept += 1;
            if known.covers(f[0], kind, sp.0, sp.1, o.min_overlap) {
                return Ok(());
            }
            known.add(f[0], kind, sp.0, sp.1);
            let mut items: Vec<String> = KEPT
                .iter()
                .filter_map(|k| match info.get(k) {
                    Some(Some(v)) => Some(format!("{k}={v}")),
                    _ => None,
                })
                .collect();
            if info.get("IMPRECISE").is_some() {
                items.push("IMPRECISE".to_owned());
            }
            items.push(format!("SVCALLER={caller}"));
            let mut row = [
                f[0],
                f[1],
                f[2],
                "N",
                kind.alt(),
                f[5],
                "PASS",
                &items.join(";"),
            ]
            .join("\t");
            if n_samples > 0 {
                row.push_str("\tGT\t");
                row.push_str(first_gt(&f));
                for _ in 1..n_samples {
                    row.push_str("\t./.");
                }
            }
            records.push(Out {
                chrom: f[0].len(),
                line: row,
                base_n: 0,
            });
            t.added += 1;
            Ok(())
        })?;
        Ok(t)
    };

    if let Some(delly) = o.delly.as_deref() {
        let par = Par::read(o.par_bed.as_deref()).map_err(|d| Error::Par(d.lines.join("\n")))?;
        let mut no_depth = 0usize;
        let mut keep = |f: &[&str], _: &Info, kind: SvKind, (s, e): (i64, i64)| {
            if f[6] != "PASS" || !call_holds(f, usize::from(f.len() > 9)) {
                return false;
            }
            let (keys, vals): (Vec<&str>, Vec<&str>) = if f.len() > 9 {
                (f[8].split(':').collect(), f[9].split(':').collect())
            } else {
                (Vec::new(), Vec::new())
            };
            let rdcn = Rdcn::from_format(&keys, &vals);
            if depth_contradicts(kind, f[0], e - s, rdcn, o.sex, par.contains(f[0], s)) {
                no_depth += 1;
                return false;
            }
            true
        };
        let t = take("delly", delly, &mut keep, &mut records)?;
        added_any |= t.added > 0;
        summary.push(format!(
            "delly: {} PASS DEL/DUP calls, {} added, and {no_depth} PASS calls longer than {DEPTH_MIN_SIZE} bp set aside because read depth (RDCN) does not show the change",
            t.kept, t.added
        ));
    }
    if let Some(cnvpytor) = o.cnvpytor.as_deref() {
        let mut keep = |_: &[&str], info: &Info, _: SvKind, _: (i64, i64)| {
            info.float("pytorP1") < o.max_eval && info.float("pytorQ0") < o.max_q0
        };
        let t = take("cnvpytor", cnvpytor, &mut keep, &mut records)?;
        added_any |= t.added > 0;
        summary.push(format!(
            "cnvpytor: {} DEL/DUP calls with e-val < {} and q0 < {}, {} added",
            t.kept,
            pycompat::fmt::general(o.max_eval, 6),
            pycompat::fmt::general(o.max_q0, 6),
            t.added
        ));
    }

    // Contig order: the header's ##contig lines (a contig listed twice
    // takes its last place, as the Python's dict does), then the contigs of
    // records in record order.
    let mut order: std::collections::HashMap<String, usize> = std::collections::HashMap::new();
    let contigs: Vec<&str> = header
        .iter()
        .filter_map(|h| h.strip_prefix("##contig=<ID="))
        .map(|r| r.split(',').next().unwrap_or("").trim_end_matches('>'))
        .collect();
    for (i, c) in contigs.iter().enumerate() {
        order.insert((*c).to_owned(), i);
    }
    let mut keyed = Vec::with_capacity(records.len());
    for r in records {
        let chrom = &r.line[..r.chrom];
        let next = order.len();
        let rank = *order.entry(chrom.to_owned()).or_insert(next);
        let pos_text = r.line[r.chrom + 1..].split('\t').next().unwrap_or("");
        let pos = int(pos_text)
            .ok_or_else(|| Error::BadPos(o.base.clone(), r.base_n, pos_text.to_owned()))?;
        keyed.push(((rank, pos), r));
    }
    keyed.sort_by_key(|(k, _)| *k);

    let present: Vec<&str> = header
        .iter()
        .filter_map(|h| h.strip_prefix("##INFO=<ID="))
        .map(|r| r.split(',').next().unwrap_or(""))
        .collect();
    let note = format!(
        "Callers beside {}: {}",
        o.base_name,
        if summary.is_empty() {
            "none given".to_owned()
        } else {
            summary.join("; ")
        }
    );
    let has_gt = header.iter().any(|h| h.starts_with("##FORMAT=<ID=GT,"));
    let mut out = io::BufWriter::with_capacity(1 << 16, out);
    let w = |out: &mut io::BufWriter<_>| -> io::Result<()> {
        for h in &header {
            if h.starts_with("##svcallers_added=") {
                continue;
            }
            if h.starts_with("#CHROM") {
                writeln!(out, "##svcallers_added=\"{note}\"")?;
                for (key, line) in STANDARD_INFO {
                    if !present.contains(&key) {
                        writeln!(out, "{line}")?;
                    }
                }
                if added_any && n_samples > 0 && !has_gt {
                    writeln!(
                        out,
                        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
                    )?;
                }
            }
            writeln!(out, "{h}")?;
        }
        for (_, r) in &keyed {
            out.write_all(r.line.as_bytes())?;
            out.write_all(b"\n")?;
        }
        out.flush()
    };
    w(&mut out).map_err(Error::Write)?;
    Ok(note)
}

/// Reads the command line; `Err` is the exit code after the message.
fn opts(argv: &[String]) -> Result<Opts, u8> {
    let p = args::parse(&SPEC, argv).map_err(|e| args::report(&SPEC, &e))?;
    let sex = p.get("--sex").unwrap_or("");
    let Some(sex) = Sex::parse(sex) else {
        let e = ArgError::Usage(format!(
            "argument --sex: invalid choice: {} (choose from 'male', 'female', '')",
            og_filters::py_repr(sex)
        ));
        return Err(args::report(&SPEC, &e));
    };
    let float = |name, default| {
        p.typed(name, default, "float", args::float)
            .map_err(|e| args::report(&SPEC, &e))
    };
    // An empty path is "not given", as the Python's `if args.delly:` reads it.
    let given = |name| p.get(name).filter(|v| !v.is_empty()).map(str::to_owned);
    Ok(Opts {
        base: p.get("--base").unwrap_or("").to_owned(),
        base_name: p.get("--base-name").unwrap_or("manta").to_owned(),
        delly: given("--delly"),
        sex,
        par_bed: given("--par-bed"),
        cnvpytor: given("--cnvpytor"),
        max_eval: float("--max-eval", 0.01)?,
        max_q0: float("--max-q0", 0.5)?,
        min_overlap: float("--min-overlap", 0.5)?,
    })
}

/// `og sv add-callers ARGS`: exit 0, 1 for an input or output problem, 2
/// for a usage error.
pub fn cli(argv: &[String]) -> u8 {
    let o = match opts(argv) {
        Ok(o) => o,
        Err(code) => return code,
    };
    match run(&o, io::stdout().lock()) {
        Ok(note) => {
            eprintln!("{note}");
            0
        }
        Err(e) => {
            eprintln!("ERROR: {e}");
            1
        }
    }
}

/// [`run`] into a string, for tests.
#[cfg(test)]
pub(crate) fn run_text(o: &Opts) -> Result<(String, String), Error> {
    let mut out = Vec::new();
    let note = run(o, &mut out)?;
    Ok((String::from_utf8(out).expect("utf-8"), note))
}
