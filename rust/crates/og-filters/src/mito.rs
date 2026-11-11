//! `og pgx mito-merge --jar pharmcat.jar --mito chrM.vcf.gz [--outside F]
//! [--depth F] [--min-af 0.95] [--min-depth 20] --out F`: an MT-RNR1
//! outside call for PharmCAT from the chrM calls (replaces
//! `bin/pgx_mito_merge.py`).
//!
//! PharmCAT 3.4.0 reads MT-RNR1 only as an outside call, so the gene was No
//! Result on every sample although step 20 calls chrM. PharmCAT's own
//! MT-RNR1 alleles are read from the jar. A chrM record with a named
//! single-base allele's position, REF and ALT counts when its FILTER is
//! PASS (or `.`) and its allele fraction (FORMAT AF, else INFO AF) is
//! `--min-af` or more. One such allele is the call; of several, the one
//! PharmCAT ranks highest (increased risk, then uncertain, then normal).
//! With none, `Reference`, but only when `--depth` shows `--min-depth`
//! reads or more at every named position. Anything that makes the answer
//! unsure holds every call back, with a line saying why ([`Held`]).
//!
//! `--outside` (step 36's calls) is copied to `--out` without its MT-RNR1
//! line, and the call is added. The summary goes to stdout. Exit 2 when
//! the jar has no MT-RNR1 definitions or an input cannot be read.
//!
//! Fixes against the Python, each with a golden in `tests/parity`:
//! - a PASS change of several bases with REF and ALT of one length (an
//!   MNV) at a named position was held back as "an insertion or deletion";
//!   it is now called "a multi-base change". The outcome is the same:
//!   nothing is passed;
//! - with `--min-depth 0` and a named position missing from the depth
//!   file, the Python crashed (`KeyError`) on the lowest depth; that
//!   position now counts as depth 0, as the threshold check already did.

use std::collections::{BTreeSet, HashMap};
use std::fmt;
use std::io::{self, Write};

use og_io::{Lines, Utf8};
use pycompat::fmt::{fixed, general};
use pycompat::json::Value;

use crate::args::{self, ArgError, Spec};

pub const MIN_AF: f64 = 0.95;
pub const MIN_DEPTH: i128 = 20;
const GENE: &str = "MT-RNR1";
const RESOURCE: &str = "org/pharmgkb/pharmcat/phenotype/MT_RNR1.json";
const CHROMS: [&str; 4] = ["chrM", "MT", "M", "chrMT"];
const NEAR_DELETION: i128 = 10;
const PREFIX: &str = "pgx_mito_merge";

const SPEC: Spec = Spec {
    prog: "og pgx mito-merge",
    usage: "usage: og pgx mito-merge --jar JAR --mito MITO [--outside OUTSIDE] [--depth DEPTH] [--min-af MIN_AF] [--min-depth MIN_DEPTH] --out OUT\n",
    options: &[
        "--jar",
        "--mito",
        "--outside",
        "--depth",
        "--min-af",
        "--min-depth",
        "--out",
    ],
    positionals: &[],
    required: &["--jar", "--mito", "--out"],
};

/// A message for `pgx_mito_merge: ERROR:` and exit 2.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Die(pub String);

impl fmt::Display for Die {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

/// An `OSError` with its Python class name, as `type(e).__name__: e`.
fn os_class(e: &io::Error, path: &str) -> String {
    let class = match e.kind() {
        io::ErrorKind::NotFound => "FileNotFoundError",
        io::ErrorKind::PermissionDenied => "PermissionError",
        io::ErrorKind::IsADirectory => "IsADirectoryError",
        io::ErrorKind::NotADirectory => "NotADirectoryError",
        io::ErrorKind::InvalidData => "UnicodeDecodeError",
        _ => "OSError",
    };
    format!("{class}: {}", crate::os_error(e, path))
}

/// `int(s)`, or Python's ValueError text.
fn py_int(s: &str) -> Result<i128, String> {
    pycompat::num::parse_int(s).map_err(|_| {
        format!(
            "ValueError: invalid literal for int() with base 10: {}",
            crate::py_repr(s)
        )
    })
}

/// A named single-base allele: `(pos, ref, alt)` and its name and function.
type Snv = HashMap<(i128, char, char), (String, String)>;

/// PharmCAT's MT-RNR1 definitions.
#[derive(Debug, Default)]
pub struct Alleles {
    snv: Snv,
    /// Positions of the other named alleles (deletions and the like), sorted.
    other: Vec<i128>,
    has_reference: bool,
}

/// `^m\.(\d+)` and what follows the digits.
fn name_pos(name: &str) -> Option<(i128, &str)> {
    let rest = name.strip_prefix("m.")?;
    let digits = rest.bytes().take_while(u8::is_ascii_digit).count();
    let pos = pycompat::num::parse_int(&rest[..digits]).ok()?;
    Some((pos, &rest[digits..]))
}

/// `^m\.(\d+)([ACGT])>([ACGT])$`, where `$` also matches before a final `\n`.
fn snv_name(name: &str) -> Option<(i128, char, char)> {
    let (pos, rest) = name_pos(name)?;
    let rest = rest.strip_suffix('\n').unwrap_or(rest);
    let b = rest.as_bytes();
    let base = |c: u8| b"ACGT".contains(&c);
    (b.len() == 3 && base(b[0]) && b[1] == b'>' && base(b[2]))
        .then(|| (pos, char::from(b[0]), char::from(b[2])))
}

/// Reads the MT-RNR1 definitions from the JSON PharmCAT ships in its jar.
pub fn parse_alleles(json: &Value, jar: &str) -> Result<Alleles, Die> {
    let shape = || {
        Die(format!(
            "cannot read {RESOURCE} from {jar} (TypeError: not PharmCAT's gene definition shape)"
        ))
    };
    let Value::Object(top) = json else {
        return Err(shape());
    };
    let list = match top.get("namedAlleles") {
        Some(Value::Array(v)) => v.as_slice(),
        Some(v) if pycompat::value::truthy(v) => return Err(shape()),
        _ => &[],
    };
    let mut a = Alleles::default();
    let mut other = BTreeSet::new();
    for item in list {
        let text = |key: &str| -> Result<String, Die> {
            match item {
                Value::Object(o) => match o.get(key) {
                    Some(Value::Str(s)) => Ok(s.clone()),
                    Some(v) if pycompat::value::truthy(v) => Err(shape()),
                    _ => Ok(String::new()),
                },
                v if pycompat::value::truthy(v) => Err(shape()),
                _ => Ok(String::new()),
            }
        };
        let name = text("name")?;
        let func = text("functionValue")?;
        if name == "Reference" {
            a.has_reference = true;
        } else if let Some(key) = snv_name(&name) {
            a.snv.insert(key, (name, func));
        } else if let Some((pos, _)) = name_pos(&name) {
            other.insert(pos);
        }
    }
    if a.snv.is_empty() {
        return Err(Die(format!(
            "{RESOURCE} in {jar} names no MT-RNR1 single-base allele"
        )));
    }
    a.other = other.into_iter().collect();
    Ok(a)
}

/// Reads PharmCAT's MT-RNR1 definitions from the jar.
pub fn read_alleles(jar: &str) -> Result<Alleles, Die> {
    let why = |what: String| Die(format!("cannot read {RESOURCE} from {jar} ({what})"));
    let bytes = crate::zip::member(jar, RESOURCE).map_err(|e| match e {
        crate::zip::ZipError::Io(e) => why(os_class(&e, jar)),
        crate::zip::ZipError::Bad(m) => why(format!("BadZipFile: {m}")),
        crate::zip::ZipError::Missing => why(format!(
            "KeyError: \"There is no item named {} in the archive\"",
            crate::py_repr(RESOURCE)
        )),
    })?;
    let text = String::from_utf8(bytes).map_err(|e| why(format!("UnicodeDecodeError: {e}")))?;
    let json = pycompat::json::loads(&text).map_err(|e| why(format!("{}: {e}", e.exception())))?;
    parse_alleles(&json, jar)
}

/// One ALT of a chrM record.
#[derive(Debug, Clone, PartialEq)]
pub struct Rec {
    pub pos: i128,
    pub ref_: String,
    pub alt: String,
    /// `PASS` for PASS or `.`, else the FILTER.
    pub filter: String,
    pub af: Option<f64>,
}

/// The AF of each ALT: FORMAT AF of the first sample, else INFO AF.
fn allele_fractions(info: &str, format: &str, sample: &str, n_alt: usize) -> Vec<Option<f64>> {
    let keys: Vec<&str> = if format.is_empty() {
        Vec::new()
    } else {
        format.split(':').collect()
    };
    let vals: Vec<&str> = if sample.is_empty() {
        Vec::new()
    } else {
        sample.split(':').collect()
    };
    let raw = match keys.iter().position(|&k| k == "AF") {
        Some(i) if i < vals.len() => Some(vals[i]),
        _ => info
            .split(';')
            .filter_map(|kv| kv.strip_prefix("AF="))
            .next_back(),
    };
    (0..n_alt)
        .map(|i| {
            raw.and_then(|r| r.split(',').nth(i))
                .and_then(|v| pycompat::num::parse_float(v).ok())
        })
        .collect()
}

/// The chrM records, one per ALT, and the count of PASS chrM records.
pub fn read_chrm(path: &str) -> Result<(Vec<Rec>, usize), Die> {
    let die = |what: String| Die(format!("cannot read {path} ({what})"));
    let input = og_io::open(path).map_err(|e| die(os_class(&e, path)))?;
    let mut recs = Vec::new();
    let mut n_pass = 0;
    for line in Lines::new(input, Utf8::Strict) {
        let line = line.map_err(|e| die(os_class(&e, path)))?;
        if line.starts_with('#') {
            continue;
        }
        let c: Vec<&str> = line.split('\t').collect();
        if c.len() < 8 {
            return Err(Die(format!("{path}: a record with {} columns", c.len())));
        }
        if !CHROMS.contains(&c[0]) {
            continue;
        }
        let ok = c[6] == "PASS" || c[6] == ".";
        n_pass += usize::from(ok);
        let alts: Vec<&str> = c[4].split(',').collect();
        let afs = allele_fractions(
            c[7],
            c.get(8).copied().unwrap_or(""),
            c.get(9).copied().unwrap_or(""),
            alts.len(),
        );
        for (alt, af) in alts.into_iter().zip(afs) {
            if alt == "*" || alt == "." || alt.starts_with('<') {
                continue;
            }
            recs.push(Rec {
                pos: py_int(c[1]).map_err(die)?,
                ref_: c[3].to_owned(),
                alt: alt.to_owned(),
                filter: if ok {
                    "PASS".to_owned()
                } else {
                    c[6].to_owned()
                },
                af,
            });
        }
    }
    Ok((recs, n_pass))
}

/// `samtools depth` rows on chrM: position to depth.
pub fn read_depth(path: &str) -> Result<HashMap<i128, i128>, Die> {
    let die = |what: String| Die(format!("cannot read {path} ({what})"));
    let file = std::fs::File::open(path).map_err(|e| die(os_class(&e, path)))?;
    let mut out = HashMap::new();
    for line in Lines::new(io::BufReader::new(file), Utf8::Strict) {
        let line = line.map_err(|e| die(os_class(&e, path)))?;
        let c = pycompat::text::split(&line, None, None);
        if c.len() >= 3 && CHROMS.contains(&c[0]) {
            out.insert(py_int(c[1]).map_err(die)?, py_int(c[2]).map_err(die)?);
        }
    }
    Ok(out)
}

/// How PharmCAT ranks a function: increased risk first.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Rank {
    Increased,
    Uncertain,
    Normal,
    Other,
}

impl Rank {
    fn of(function: &str) -> Rank {
        let f = pycompat::text::lower(function);
        [
            ("increased", Rank::Increased),
            ("uncertain", Rank::Uncertain),
            ("normal", Rank::Normal),
        ]
        .into_iter()
        .find_map(|(w, r)| f.starts_with(w).then_some(r))
        .unwrap_or(Rank::Other)
    }
}

/// Why no call is passed although a record was read.
#[derive(Debug, Clone, PartialEq)]
pub enum Held {
    /// A homoplasmic PASS change at a named position PharmCAT does not name.
    UnnamedChange {
        shown: String,
    },
    Filtered {
        name: String,
        shown: String,
        filter: String,
    },
    NoAf {
        name: String,
        shown: String,
    },
    Heteroplasmic {
        name: String,
        shown: String,
        min_af: f64,
    },
    /// A PASS indel, or an MNV, at or near a named position.
    Near {
        shown: String,
        pos: i128,
        mnv: bool,
    },
}

impl fmt::Display for Held {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Held::UnnamedChange { shown } => {
                write!(
                    f,
                    "{shown}: a change at an MT-RNR1 position that PharmCAT does not name"
                )
            }
            Held::Filtered {
                name,
                shown,
                filter,
            } => write!(f, "{name} ({shown}) is filtered ({filter})"),
            Held::NoAf { name, shown } => write!(f, "{name} ({shown}) has no allele fraction"),
            Held::Heteroplasmic {
                name,
                shown,
                min_af,
            } => write!(
                f,
                "{name} ({shown}) is heteroplasmic (AF below {}); PharmCAT takes one homoplasmic allele",
                general(*min_af, 6)
            ),
            Held::Near { shown, pos, mnv } => {
                let what = if *mnv {
                    "a multi-base change"
                } else {
                    "an insertion or deletion"
                };
                write!(
                    f,
                    "{what} ({shown}) at or near m.{pos}, which this script does not map to PharmCAT's names"
                )
            }
        }
    }
}

/// A homoplasmic named allele.
#[derive(Debug, Clone, PartialEq, PartialOrd)]
struct Called {
    rank: Rank,
    pos: i128,
    name: String,
    function: String,
    shown: String,
}

/// The MT-RNR1 call, if any, and the lines that say why.
pub fn decide(
    recs: &[Rec],
    n_pass: usize,
    a: &Alleles,
    min_af: f64,
    depth: Option<&HashMap<i128, i128>>,
    min_depth: i128,
) -> (Option<String>, Vec<String>) {
    let named: BTreeSet<i128> = a.snv.keys().map(|k| k.0).collect();
    let mut calls: Vec<Called> = Vec::new();
    let mut held: Vec<Held> = Vec::new();
    for r in recs {
        let shown = format!(
            "chrM:{} {}>{}{}",
            r.pos,
            r.ref_,
            r.alt,
            r.af.map_or(" (no AF)".to_owned(), |af| format!(" AF {}", fixed(af, 3)))
        );
        let one = |s: &str| {
            let mut c = s.chars();
            c.next().filter(|_| c.next().is_none())
        };
        if let (Some(rb), Some(ab)) = (one(&r.ref_), one(&r.alt)) {
            let Some((name, function)) = a.snv.get(&(r.pos, rb, ab)) else {
                let at_named = named.contains(&r.pos) || a.other.contains(&r.pos);
                if at_named && r.filter == "PASS" && r.af.is_none_or(|af| af >= min_af) {
                    held.push(Held::UnnamedChange { shown });
                }
                continue;
            };
            let (name, shown) = (name.clone(), shown);
            match r.af {
                _ if r.filter != "PASS" => held.push(Held::Filtered {
                    name,
                    shown,
                    filter: r.filter.clone(),
                }),
                None => held.push(Held::NoAf { name, shown }),
                Some(af) if af < min_af => held.push(Held::Heteroplasmic {
                    name,
                    shown,
                    min_af,
                }),
                Some(_) => calls.push(Called {
                    rank: Rank::of(function),
                    pos: r.pos,
                    name,
                    function: function.clone(),
                    shown,
                }),
            }
            continue;
        }
        if r.filter != "PASS" {
            continue;
        }
        let ref_len = i128::try_from(crate::chars(&r.ref_)).unwrap_or(i128::MAX);
        let start = r.pos;
        // An insertion sits between pos and pos + 1.
        let end = r.pos + ref_len.max(2) - 1;
        let near = named
            .iter()
            .filter(|&&p| start <= p && p <= end)
            .chain(
                a.other
                    .iter()
                    .filter(|&&p| start - NEAR_DELETION <= p && p <= end + NEAR_DELETION),
            )
            .min();
        if let Some(&pos) = near {
            let mnv = crate::chars(&r.ref_) == crate::chars(&r.alt);
            held.push(Held::Near { shown, pos, mnv });
        }
    }
    let mut lines = Vec::new();
    if n_pass == 0 {
        lines.push(format!(
            "the chrM call set has no PASS record (chrM not sequenced or not called): no {GENE} call passed, it stays No Result"
        ));
        return (None, lines);
    }
    if !held.is_empty() {
        lines.extend(held.iter().map(|h| format!("held back: {h}")));
        lines.push(format!(
            "no {GENE} call passed to PharmCAT: it stays No Result; read the records above in the chrM call set"
        ));
        return (None, lines);
    }
    if !calls.is_empty() {
        calls.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
        let best = &calls[0];
        lines.push(format!(
            "{GENE} {} passed to PharmCAT as an outside call (homoplasmic, {}; {})",
            best.name, best.shown, best.function
        ));
        for c in &calls[1..] {
            lines.push(format!(
                "also homoplasmic, not passed (PharmCAT takes one {GENE} allele): {} ({}; {})",
                c.name, c.shown, c.function
            ));
        }
        return (Some(best.name.clone()), lines);
    }
    if !a.has_reference {
        lines.push(format!(
            "no named {GENE} allele, and PharmCAT's definitions have no Reference: no call passed"
        ));
        return (None, lines);
    }
    let Some(depth) = depth else {
        lines.push(format!(
            "no named {GENE} allele is called, but there is no read depth to show its positions were covered (no BAM): Reference not passed, {GENE} stays No Result"
        ));
        return (None, lines);
    };
    let positions: BTreeSet<i128> = named.iter().chain(&a.other).copied().collect();
    let at = |p: &i128| depth.get(p).copied().unwrap_or(0);
    let low: Vec<(i128, i128)> = positions
        .iter()
        .filter(|p| at(p) < min_depth)
        .map(|p| (*p, at(p)))
        .collect();
    if !low.is_empty() {
        let mut shown = low
            .iter()
            .take(5)
            .map(|(p, d)| format!("chrM:{p} depth {d}"))
            .collect::<Vec<_>>()
            .join(", ");
        if low.len() > 5 {
            shown.push_str(", ...");
        }
        lines.push(format!(
            "no named {GENE} allele is called, but {} of its {} positions have depth below {min_depth} ({shown}): Reference not passed, {GENE} stays No Result",
            low.len(),
            positions.len()
        ));
        return (None, lines);
    }
    let lowest = positions.iter().map(at).min().unwrap_or(0);
    lines.push(format!(
        "{GENE} Reference passed to PharmCAT as an outside call: none of PharmCAT's {} named single-base alleles is called, and all {} named positions have depth {min_depth} or more (lowest {lowest}); {n_pass} PASS chrM records read",
        a.snv.len(),
        positions.len()
    ));
    (Some("Reference".to_owned()), lines)
}

/// The options of one run.
struct Opts {
    jar: String,
    mito: String,
    outside: Option<String>,
    depth: Option<String>,
    min_af: f64,
    min_depth: i128,
    out: String,
}

fn say(line: &str) {
    println!("{PREFIX}: {line}");
}

fn merge(o: &Opts) -> Result<(), Die> {
    let alleles = read_alleles(&o.jar)?;
    let (recs, n_pass) = read_chrm(&o.mito)?;
    let depth = o.depth.as_deref().map(read_depth).transpose()?;
    let (call, lines) = decide(
        &recs,
        n_pass,
        &alleles,
        o.min_af,
        depth.as_ref(),
        o.min_depth,
    );
    for l in &lines {
        say(l);
    }
    let mut out = String::new();
    if let Some(path) = &o.outside {
        let die = |e: io::Error| {
            Die(format!(
                "cannot read {path} ({})",
                crate::os_error(&e, path)
            ))
        };
        let file = std::fs::File::open(path).map_err(die)?;
        let mut lines = Lines::new(io::BufReader::new(file), Utf8::Strict);
        let mut line = String::new();
        while lines.read_line(&mut line).map_err(die)? {
            if pycompat::num::py_strip(&line).is_empty() {
                continue;
            }
            let first = line.split('\t').next().unwrap_or_default();
            if pycompat::num::py_strip(first) == GENE {
                say(&format!(
                    "dropped the {GENE} line of {path}: the chrM calls decide it"
                ));
                continue;
            }
            out.push_str(&line);
            if !line.ends_with('\n') {
                out.push('\n');
            }
        }
    }
    if let Some(c) = call {
        out.push_str(&format!("{GENE}\t{c}\n"));
    }
    let write = || -> io::Result<()> {
        let mut f = og_io::AtomicFile::create(&o.out)?;
        f.write_all(out.as_bytes())?;
        f.commit()
    };
    write().map_err(|e| {
        Die(format!(
            "cannot write {} ({})",
            o.out,
            crate::os_error(&e, &o.out)
        ))
    })
}

/// `og pgx mito-merge`: exit 0, or 2 for a usage error or unreadable input.
pub fn cli(argv: &[String]) -> u8 {
    let parse = || -> Result<Opts, ArgError> {
        let a = args::parse(&SPEC, argv)?;
        let own = |k: &str| a.get(k).map(str::to_owned);
        Ok(Opts {
            jar: own("--jar").unwrap_or_default(),
            mito: own("--mito").unwrap_or_default(),
            outside: own("--outside").filter(|s| !s.is_empty()),
            depth: own("--depth").filter(|s| !s.is_empty()),
            min_af: a.typed("--min-af", MIN_AF, "float", args::float)?,
            min_depth: a.typed("--min-depth", MIN_DEPTH, "int", args::int)?,
            out: own("--out").unwrap_or_default(),
        })
    };
    let o = match parse() {
        Ok(o) => o,
        Err(e) => return args::report(&SPEC, &e),
    };
    match merge(&o) {
        Ok(()) => 0,
        Err(Die(msg)) => {
            let _ = io::stdout().flush();
            eprintln!("{PREFIX}: ERROR: {msg}");
            2
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn alleles() -> Alleles {
        let json = r#"{"namedAlleles": [
            {"name": "Reference", "functionValue": "Normal risk"},
            {"name": "m.827A>G", "functionValue": "Normal risk"},
            {"name": "m.961T>del", "functionValue": "Uncertain risk"},
            {"name": "m.1095T>C", "functionValue": "Increased risk"},
            {"name": "m.1555A>G", "functionValue": "Increased risk"}]}"#;
        parse_alleles(&pycompat::json::loads(json).unwrap(), "j").unwrap()
    }

    fn rec(pos: i128, r: &str, alt: &str, af: Option<f64>) -> Rec {
        Rec {
            pos,
            ref_: r.into(),
            alt: alt.into(),
            filter: "PASS".into(),
            af,
        }
    }

    #[test]
    fn red_first_an_af_at_the_floor_is_homoplasmic() {
        let a = alleles();
        let (call, lines) = decide(&[rec(1555, "A", "G", Some(0.95))], 1, &a, 0.95, None, 20);
        assert_eq!(call.as_deref(), Some("m.1555A>G"), "{lines:?}");
        let (call, lines) = decide(&[rec(1555, "A", "G", Some(0.949))], 1, &a, 0.95, None, 20);
        assert_eq!(call, None);
        assert_eq!(
            lines[0],
            "held back: m.1555A>G (chrM:1555 A>G AF 0.949) is heteroplasmic (AF below 0.95); PharmCAT takes one homoplasmic allele"
        );
    }

    #[test]
    fn two_calls_the_increased_risk_one_wins() {
        let a = alleles();
        let r = [
            rec(827, "A", "G", Some(1.0)),
            rec(1555, "A", "G", Some(0.99)),
        ];
        let (call, lines) = decide(&r, 2, &a, 0.95, None, 20);
        assert_eq!(call.as_deref(), Some("m.1555A>G"));
        assert!(lines[1].starts_with("also homoplasmic, not passed"));
    }

    #[test]
    fn an_mnv_is_named_as_one() {
        let a = alleles();
        let (call, lines) = decide(&[rec(1554, "TA", "CG", Some(1.0))], 1, &a, 0.95, None, 20);
        assert_eq!(call, None);
        assert_eq!(
            lines[0],
            "held back: a multi-base change (chrM:1554 TA>CG AF 1.000) at or near m.1555, which this script does not map to PharmCAT's names"
        );
        let (_, lines) = decide(&[rec(955, "AC", "A", Some(1.0))], 1, &a, 0.95, None, 20);
        assert!(lines[0].starts_with(
            "held back: an insertion or deletion (chrM:955 AC>A AF 1.000) at or near m.961"
        ));
    }

    #[test]
    fn reference_needs_depth_everywhere() {
        let a = alleles();
        let r = [rec(263, "A", "G", Some(1.0))];
        let mut d: HashMap<i128, i128> = [827, 961, 1095, 1555]
            .into_iter()
            .map(|p| (p, 30))
            .collect();
        let (call, _) = decide(&r, 1, &a, 0.95, Some(&d), 20);
        assert_eq!(call.as_deref(), Some("Reference"));
        d.insert(961, 19);
        let (call, lines) = decide(&r, 1, &a, 0.95, Some(&d), 20);
        assert_eq!(call, None);
        assert!(lines[0].contains("1 of its 4 positions have depth below 20 (chrM:961 depth 19)"));
        let (call, _) = decide(&r, 1, &a, 0.95, None, 20);
        assert_eq!(call, None);
    }

    #[test]
    fn af_from_format_then_info() {
        assert_eq!(
            allele_fractions("AF=0.5", "GT:AF", "0/1:0.9,0.1", 2),
            [Some(0.9), Some(0.1)]
        );
        assert_eq!(
            allele_fractions("X=1;AF=0.3,0.7", "GT", "0/1", 2),
            [Some(0.3), Some(0.7)]
        );
        assert_eq!(allele_fractions(".", "", "", 1), [None]);
    }

    #[test]
    fn names_as_python_reads_them() {
        assert_eq!(snv_name("m.1555A>G"), Some((1555, 'A', 'G')));
        assert_eq!(snv_name("m.1555A>G\n"), Some((1555, 'A', 'G')));
        assert_eq!(snv_name("m.961T>del"), None);
        assert_eq!(name_pos("m.961T>del").map(|p| p.0), Some(961));
        assert_eq!(name_pos("Reference"), None);
    }
}
