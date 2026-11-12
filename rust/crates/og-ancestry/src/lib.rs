//! Step 26b, fine-grained ancestry (`docs/26b-ancestry-fine.md`): the
//! commands around ADMIXTURE and RFMix, as `og ancestry-fine <verb>`.
//!
//! - [`panel`]: `panel-map` writes the panel's sample table once, at setup;
//!   `check-map` refuses a table given with `ANCESTRY_FINE_SAMPLE_MAP` that
//!   has a group too small to model.
//! - [`order`]: the group of each column of a supervised ADMIXTURE `.Q`.
//! - [`sample_bed`]: ADMIXTURE's input for one person, and the bootstrap
//!   script the ADMIXTURE image runs.
//! - [`chroms`]: the chromosomes long enough for RFMix to paint.
//! - [`report`]: the fractions table, the painting (BED and PNG) and the
//!   plain-text summary.
//!
//! Each command writes what `bin/ancestry_fine_report.py` wrote (checked by
//! `tests/parity`), with three deliberate changes: every output is written
//! under a temporary name and renamed into place only once every input has
//! been read and checked, so a failure leaves no partial ADMIXTURE input; a
//! missing or malformed input is an `ERROR:` line and exit 1, never a
//! traceback; and the PNG's pixels are the contract, not its compressed
//! bytes.
//!
//! A group name is ours: a set of reference populations we put together
//! and named, not a nationality.

pub mod chroms;
pub mod order;
pub mod panel;
pub mod report;
pub mod sample_bed;

use std::fmt;
use std::io;
use std::path::Path;

/// Population of the panel (gnomAD's HGDP + 1000 Genomes, v2 metadata
/// labels) to (continental group, regional group). Admixed populations,
/// between the continental groups, are in no group and so never a
/// reference: ACB, ASW, CLM, MXL, PEL, PUR, Hazara and Uygur.
const GROUPS: &[(&str, &str, &str)] = &[
    // Africa
    ("YRI", "AFR", "WAF"),
    ("ESN", "AFR", "WAF"),
    ("GWD", "AFR", "WAF"),
    ("MSL", "AFR", "WAF"),
    ("Yoruba", "AFR", "WAF"),
    ("Mandenka", "AFR", "WAF"),
    ("LWK", "AFR", "EAF"),
    ("BantuKenya", "AFR", "EAF"),
    ("Biaka", "AFR", "CAF"),
    ("Mbuti", "AFR", "CAF"),
    ("San", "AFR", "SAF"),
    ("BantuSouthAfrica", "AFR", "SAF"),
    // Europe, the Middle East and North Africa
    ("CEU", "WEA", "NWE"),
    ("GBR", "WEA", "NWE"),
    ("Orcadian", "WEA", "NWE"),
    ("French", "WEA", "NWE"),
    ("FIN", "WEA", "NEE"),
    ("Russian", "WEA", "NEE"),
    ("IBS", "WEA", "IBE"),
    ("Basque", "WEA", "IBE"),
    ("TSI", "WEA", "ITA"),
    ("Tuscan", "WEA", "ITA"),
    ("BergamoItalian", "WEA", "ITA"),
    ("Sardinian", "WEA", "SAR"),
    ("Adygei", "WEA", "CAU"),
    ("Druze", "WEA", "LEV"),
    ("Palestinian", "WEA", "LEV"),
    ("Bedouin", "WEA", "BDN"),
    ("Mozabite", "WEA", "NAF"),
    // Central and South Asia
    ("GIH", "CSA", "CSA"),
    ("PJL", "CSA", "CSA"),
    ("BEB", "CSA", "CSA"),
    ("STU", "CSA", "CSA"),
    ("ITU", "CSA", "CSA"),
    ("Balochi", "CSA", "CSA"),
    ("Brahui", "CSA", "CSA"),
    ("Makrani", "CSA", "CSA"),
    ("Pathan", "CSA", "CSA"),
    ("Sindhi", "CSA", "CSA"),
    ("Burusho", "CSA", "CSA"),
    ("Kalash", "CSA", "CSA"),
    // East Asia and Siberia
    ("CHB", "EAS", "EAS"),
    ("CHS", "EAS", "EAS"),
    ("CDX", "EAS", "EAS"),
    ("JPT", "EAS", "EAS"),
    ("KHV", "EAS", "EAS"),
    ("Han", "EAS", "EAS"),
    ("NorthernHan", "EAS", "EAS"),
    ("Japanese", "EAS", "EAS"),
    ("Dai", "EAS", "EAS"),
    ("Daur", "EAS", "EAS"),
    ("Hezhen", "EAS", "EAS"),
    ("Lahu", "EAS", "EAS"),
    ("Miao", "EAS", "EAS"),
    ("Mongolian", "EAS", "EAS"),
    ("Naxi", "EAS", "EAS"),
    ("Oroqen", "EAS", "EAS"),
    ("She", "EAS", "EAS"),
    ("Tu", "EAS", "EAS"),
    ("Tujia", "EAS", "EAS"),
    ("Xibo", "EAS", "EAS"),
    ("Yi", "EAS", "EAS"),
    ("Cambodian", "EAS", "EAS"),
    ("Yakut", "EAS", "SIB"),
    // The Americas (HGDP's Native American populations)
    ("Karitiana", "AMR", "AMR"),
    ("Maya", "AMR", "AMR"),
    ("Pima", "AMR", "AMR"),
    ("Surui", "AMR", "AMR"),
    ("Colombian", "AMR", "AMR"),
    // Oceania
    ("Bougainville", "OCE", "OCE"),
    ("PapuanHighlands", "OCE", "OCE"),
    ("PapuanSepik", "OCE", "OCE"),
];

/// The (continental, regional) groups of a population, if it has them.
pub fn groups_of(population: &str) -> Option<(&'static str, &'static str)> {
    GROUPS
        .iter()
        .find(|(p, _, _)| *p == population)
        .map(|&(_, c, r)| (c, r))
}

/// The two levels of the model. Their names are the first column of
/// `levels.tsv` and part of every output file name.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Level {
    Continent,
    Region,
}

impl Level {
    /// Both levels, in the order the report lists them.
    pub const ALL: [Level; 2] = [Level::Continent, Level::Region];

    pub fn parse(s: &str) -> Option<Level> {
        match s {
            "continent" => Some(Level::Continent),
            "region" => Some(Level::Region),
            _ => None,
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            Level::Continent => "continent",
            Level::Region => "region",
        }
    }

    /// The column of the level's group in `samples.tsv`.
    fn column(self) -> usize {
        match self {
            Level::Continent => 2,
            Level::Region => 3,
        }
    }

    fn labels(self) -> &'static [(&'static str, &'static str)] {
        match self {
            Level::Continent => &[
                ("AFR", "Africa"),
                ("WEA", "West Eurasia (Europe, Middle East, North Africa)"),
                ("CSA", "Central and South Asia"),
                ("EAS", "East Asia"),
                ("AMR", "Americas (Native American)"),
                ("OCE", "Oceania"),
            ],
            Level::Region => &[
                ("WAF", "West Africa"),
                ("EAF", "East Africa"),
                ("CAF", "Central African rainforest"),
                ("SAF", "Southern Africa"),
                ("NWE", "Northwestern Europe"),
                ("NEE", "Northeastern Europe"),
                ("IBE", "Iberia"),
                ("ITA", "Italy"),
                ("SAR", "Sardinia"),
                ("CAU", "Caucasus"),
                ("LEV", "Levant"),
                ("BDN", "Bedouin"),
                ("NAF", "North Africa"),
                ("CSA", "Central and South Asia"),
                ("EAS", "East Asia"),
                ("SIB", "Siberia"),
                ("AMR", "Americas (Native American)"),
                ("OCE", "Oceania"),
            ],
        }
    }
}

/// The readable name of a group code at a level; the code itself when it
/// has none.
pub fn label(level: Level, code: &str) -> &str {
    level
        .labels()
        .iter()
        .find(|(c, _)| *c == code)
        .map_or(code, |&(_, l)| l)
}

/// Painting colour of each continental group.
pub const COLOURS: &[(&str, [u8; 3])] = &[
    ("AFR", [230, 159, 0]),
    ("WEA", [0, 114, 178]),
    ("CSA", [0, 158, 115]),
    ("EAS", [204, 121, 167]),
    ("AMR", [213, 94, 0]),
    ("OCE", [86, 180, 233]),
];

/// The fewest unrelated reference people a group needs to be modelled.
/// RFMix fits its CRF weight on people it simulates from `floor(n / 10) +
/// 1` of each group's references and leaves those parents out of the
/// reference meanwhile, so a group of one has no reference haplotype at all
/// and RFMix can segfault; two is its hard floor. We ask for 10 because a
/// group is meant to stand for a population.
pub const MIN_REFERENCES: usize = 10;

/// Why a command stopped: printed as `ERROR: <message>`, exit 1.
#[derive(Debug, PartialEq, Eq)]
pub struct Fail(pub String);

impl fmt::Display for Fail {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

impl From<String> for Fail {
    fn from(s: String) -> Self {
        Fail(s)
    }
}

pub type Result<T> = std::result::Result<T, Fail>;

/// An I/O error's text without Rust's ` (os error N)`.
fn reason(e: &io::Error) -> String {
    let s = e.to_string();
    match s.rfind(" (os error ") {
        Some(i) if s.ends_with(')') => s[..i].to_owned(),
        _ => s,
    }
}

fn cannot_read(path: &Path, e: &io::Error) -> Fail {
    Fail(format!("cannot read {}: {}", path.display(), reason(e)))
}

fn cannot_write(path: &Path, e: &io::Error) -> Fail {
    Fail(format!("cannot write {}: {}", path.display(), reason(e)))
}

/// The lines of a text file as Python's `open(path)` iterates them
/// (universal newlines, UTF-8), each without its ending. gzip is read
/// through when the file starts with its magic.
fn read_lines(path: &Path) -> Result<Vec<String>> {
    let mut out = Vec::new();
    for_each_line(path, |l| {
        out.push(l.to_owned());
        Ok(())
    })?;
    Ok(out)
}

/// Calls `f` on each line of `path` (as [`read_lines`] reads them), one at
/// a time, so a file of any size streams.
fn for_each_line(path: &Path, mut f: impl FnMut(&str) -> Result<()>) -> Result<()> {
    let input = og_io::open(path).map_err(|e| cannot_read(path, &e))?;
    let mut lines = og_io::Lines::new(input, og_io::Utf8::Strict);
    let mut line = String::new();
    while lines
        .read_line(&mut line)
        .map_err(|e| cannot_read(path, &e))?
    {
        f(line.strip_suffix('\n').unwrap_or(&line))?;
    }
    Ok(())
}

/// The raw lines of a text file with their endings (`f.readlines()`),
/// read by `f`, one at a time.
fn for_each_raw_line(path: &Path, mut f: impl FnMut(&str) -> Result<()>) -> Result<()> {
    let input = og_io::open(path).map_err(|e| cannot_read(path, &e))?;
    let mut lines = og_io::Lines::new(input, og_io::Utf8::Strict);
    let mut line = String::new();
    while lines
        .read_line(&mut line)
        .map_err(|e| cannot_read(path, &e))?
    {
        f(&line)?;
    }
    Ok(())
}

/// Python's `s.split()`: runs of whitespace separate, none at the ends.
fn split_ws(s: &str) -> impl Iterator<Item = &str> {
    s.split(pycompat::num::py_isspace).filter(|w| !w.is_empty())
}

/// Python's `int(s)`, to a type that holds every value a caller needs.
fn int(s: &str) -> Option<i64> {
    pycompat::num::parse_int(s)
        .ok()
        .and_then(|n| i64::try_from(n).ok())
}

/// Python's `repr()` of a str.
fn repr(s: &str) -> String {
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
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if c == quote => {
                out.push('\\');
                out.push(c);
            }
            c if (c as u32) < 0x20 || c as u32 == 0x7f => {
                out.push_str(&format!("\\x{:02x}", c as u32));
            }
            c => out.push(c),
        }
    }
    out.push(quote);
    out
}

/// Python's `repr()` of a list of str.
fn repr_list<S: AsRef<str>>(items: &[S]) -> String {
    let parts: Vec<String> = items.iter().map(|s| repr(s.as_ref())).collect();
    format!("[{}]", parts.join(", "))
}

/// How chromosomes sort: numbers in numeric order, then the rest by name;
/// a `chr` prefix is ignored.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
enum ChromKey<'a> {
    Number(i128),
    Name(&'a str),
}

fn chrom_key(c: &str) -> ChromKey<'_> {
    let bare = c.strip_prefix("chr").unwrap_or(c);
    if pycompat::text::isdigit(bare)
        && let Ok(n) = pycompat::num::parse_int(bare)
    {
        return ChromKey::Number(n);
    }
    ChromKey::Name(bare)
}

/// A file written in full under a temporary name, renamed into place by
/// [`Pending::commit_all`] once every output of a command is ready.
struct Pending {
    files: Vec<(std::path::PathBuf, og_io::AtomicFile)>,
}

impl Pending {
    fn new() -> Self {
        Pending { files: Vec::new() }
    }

    /// A new output at `path`, written through the returned writer.
    fn create(&mut self, path: std::path::PathBuf) -> Result<&mut og_io::AtomicFile> {
        let f = og_io::AtomicFile::create(&path).map_err(|e| cannot_write(&path, &e))?;
        self.files.push((path, f));
        Ok(&mut self.files.last_mut().expect("just pushed").1)
    }

    /// A whole output at `path`.
    fn put(&mut self, path: std::path::PathBuf, data: &[u8]) -> Result<()> {
        use std::io::Write;
        let shown = path.clone();
        self.create(path)?
            .write_all(data)
            .map_err(|e| cannot_write(&shown, &e))
    }

    /// Renames each file into place, in the order they were created. Until
    /// this is called no output exists; dropped without it, the temporary
    /// files are removed.
    fn commit_all(self) -> Result<()> {
        for (path, f) in self.files {
            f.commit().map_err(|e| cannot_write(&path, &e))?;
        }
        Ok(())
    }
}

/// A command line, parsed as the Python's argparse took it: `--name VALUE`
/// or `--name=VALUE`, an option given twice keeps the last value except
/// where the command collects every one (`--level`).
struct Args {
    given: Vec<(&'static str, String)>,
}

/// What a verb takes.
struct Spec {
    verb: &'static str,
    usage: &'static str,
    options: &'static [&'static str],
    required: &'static [&'static str],
}

enum ArgError {
    Help,
    Usage(String),
}

impl Args {
    fn parse(spec: &Spec, args: &[String]) -> std::result::Result<Args, ArgError> {
        let mut given = Vec::new();
        let mut extra = Vec::new();
        let mut i = 0;
        while i < args.len() {
            let a = &args[i];
            i += 1;
            if a == "-h" || a == "--help" {
                return Err(ArgError::Help);
            }
            let (name, inline) = match a.split_once('=') {
                Some((n, v)) if n.starts_with("--") => (n, Some(v.to_owned())),
                _ => (a.as_str(), None),
            };
            let Some(&known) = spec.options.iter().find(|&&o| o == name) else {
                extra.push(a.clone());
                continue;
            };
            let value = match inline {
                Some(v) => v,
                None => match args.get(i) {
                    Some(v) if !v.starts_with('-') || v == "-" || is_negative_number(v) => {
                        i += 1;
                        v.clone()
                    }
                    _ => {
                        return Err(ArgError::Usage(format!(
                            "argument {name}: expected one argument"
                        )));
                    }
                },
            };
            given.push((known, value));
        }
        let missing: Vec<&str> = spec
            .required
            .iter()
            .copied()
            .filter(|r| !given.iter().any(|(n, _)| n == r))
            .collect();
        if !missing.is_empty() {
            return Err(ArgError::Usage(format!(
                "the following arguments are required: {}",
                missing.join(", ")
            )));
        }
        if !extra.is_empty() {
            return Err(ArgError::Usage(format!(
                "unrecognized arguments: {}",
                extra.join(" ")
            )));
        }
        Ok(Args { given })
    }

    /// The last value of `--name`.
    fn get(&self, name: &str) -> Option<&str> {
        self.given
            .iter()
            .rev()
            .find(|(n, _)| *n == name)
            .map(|(_, v)| v.as_str())
    }

    /// A required option's value.
    fn req(&self, name: &str) -> &str {
        self.get(name).expect("checked by parse")
    }

    /// Every value of `--name`, in order.
    fn all(&self, name: &str) -> Vec<&str> {
        self.given
            .iter()
            .filter(|(n, _)| *n == name)
            .map(|(_, v)| v.as_str())
            .collect()
    }

    /// `--name` as argparse's `type=int`, or `default`.
    fn int(&self, name: &str, default: i64) -> std::result::Result<i64, ArgError> {
        self.typed(name, default, "int", int)
    }

    /// `--name` as argparse's `type=float`, or `default`.
    fn float(&self, name: &str, default: f64) -> std::result::Result<f64, ArgError> {
        self.typed(name, default, "float", |s| {
            pycompat::num::parse_float(s).ok()
        })
    }

    fn typed<T>(
        &self,
        name: &str,
        default: T,
        kind: &str,
        parse: impl Fn(&str) -> Option<T>,
    ) -> std::result::Result<T, ArgError> {
        match self.get(name) {
            None => Ok(default),
            Some(v) => parse(v).ok_or_else(|| {
                ArgError::Usage(format!(
                    "argument {name}: invalid {kind} value: {}",
                    repr(v)
                ))
            }),
        }
    }
}

fn is_negative_number(s: &str) -> bool {
    s.strip_prefix('-')
        .is_some_and(|d| !d.is_empty() && d.parse::<f64>().is_ok())
}

/// A verb's body: the parsed options and stdout.
type Run = fn(&Args, &mut dyn io::Write) -> Result<()>;

const USAGE: &str =
    "usage: og ancestry-fine {panel-map,check-map,order,sample-bed,chroms,report} ...\n";

/// `og ancestry-fine <verb> [args]`: the exit code. Output goes to `out`,
/// messages to `err`.
pub fn cli(args: &[String], out: &mut dyn io::Write, err: &mut dyn io::Write) -> u8 {
    let Some(verb) = args.first() else {
        let _ = write!(err, "{USAGE}");
        let _ = writeln!(
            err,
            "og ancestry-fine: error: the following arguments are required: verb"
        );
        return 2;
    };
    let rest = &args[1..];
    let (spec, run): (&Spec, Run) = match verb.as_str() {
        "panel-map" => (&panel::PANEL_MAP, panel::panel_map_main),
        "check-map" => (&panel::CHECK_MAP, panel::check_map_main),
        "order" => (&order::SPEC, order::main),
        "sample-bed" => (&sample_bed::SPEC, sample_bed::main),
        "chroms" => (&chroms::SPEC, chroms::main),
        "report" => (&report::SPEC, report::main),
        "-h" | "--help" => {
            let _ = write!(out, "{USAGE}");
            return 0;
        }
        other => {
            let _ = write!(err, "{USAGE}");
            let _ = writeln!(
                err,
                "og ancestry-fine: error: argument verb: invalid choice: {}",
                repr(other)
            );
            return 2;
        }
    };
    let parsed = Args::parse(spec, rest).and_then(|a| {
        check_types(spec, &a)?;
        Ok(a)
    });
    let a = match parsed {
        Ok(a) => a,
        Err(ArgError::Help) => {
            let _ = write!(out, "{}", spec.usage);
            return 0;
        }
        Err(ArgError::Usage(msg)) => {
            let _ = write!(err, "{}", spec.usage);
            let _ = writeln!(err, "og ancestry-fine {}: error: {msg}", spec.verb);
            return 2;
        }
    };
    match run(&a, out) {
        Ok(()) => 0,
        Err(Fail(msg)) => {
            let _ = out.flush();
            let _ = writeln!(err, "ERROR: {msg}");
            1
        }
    }
}

/// The typed options argparse checked before the command ran.
fn check_types(spec: &Spec, a: &Args) -> std::result::Result<(), ArgError> {
    for &name in spec.options {
        match name {
            "--paint-per-group" | "--bootstrap" | "--min-sites" => {
                a.int(name, 0)?;
            }
            "--min-cm" => {
                a.float(name, 0.0)?;
            }
            _ => {}
        }
    }
    Ok(())
}

#[cfg(test)]
pub(crate) mod testutil {
    use std::path::PathBuf;

    /// A fresh directory under the target's temp dir, removed on drop.
    pub struct Dir(pub PathBuf);

    impl Dir {
        pub fn new(name: &str) -> Dir {
            let p = std::env::temp_dir().join(format!(
                "og-ancestry-{name}-{}-{:?}",
                std::process::id(),
                std::thread::current().id()
            ));
            let _ = std::fs::remove_dir_all(&p);
            std::fs::create_dir_all(&p).unwrap();
            Dir(p)
        }

        pub fn put(&self, name: &str, text: &str) -> String {
            let p = self.0.join(name);
            if let Some(d) = p.parent() {
                std::fs::create_dir_all(d).unwrap();
            }
            std::fs::write(&p, text).unwrap();
            p.to_string_lossy().into_owned()
        }

        pub fn path(&self, name: &str) -> String {
            self.0.join(name).to_string_lossy().into_owned()
        }

        pub fn read(&self, name: &str) -> String {
            std::fs::read_to_string(self.0.join(name)).unwrap()
        }
    }

    impl Drop for Dir {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }

    /// Runs `og ancestry-fine ARGS`: (exit code, stdout, stderr).
    pub fn run(args: &[&str]) -> (u8, String, String) {
        let args: Vec<String> = args.iter().map(|s| (*s).to_owned()).collect();
        let (mut out, mut err) = (Vec::new(), Vec::new());
        let rc = super::cli(&args, &mut out, &mut err);
        (
            rc,
            String::from_utf8(out).unwrap(),
            String::from_utf8(err).unwrap(),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn labels_and_groups() {
        assert_eq!(
            label(Level::Continent, "WEA"),
            "West Eurasia (Europe, Middle East, North Africa)"
        );
        assert_eq!(label(Level::Region, "IBE"), "Iberia");
        assert_eq!(label(Level::Region, "XYZ"), "XYZ");
        assert_eq!(groups_of("Basque"), Some(("WEA", "IBE")));
        assert_eq!(groups_of("ASW"), None);
        assert_eq!(GROUPS.len(), 72);
    }

    #[test]
    fn the_minimum_is_ten() {
        assert_eq!(MIN_REFERENCES, 10);
    }

    #[test]
    fn chromosomes_sort_numbers_first() {
        let mut v = vec!["chrX", "chr10", "chr2", "chrM", "1"];
        v.sort_by(|a, b| chrom_key(a).cmp(&chrom_key(b)));
        assert_eq!(v, ["1", "chr2", "chr10", "chrM", "chrX"]);
    }

    #[test]
    fn reprs_as_python() {
        assert_eq!(repr_list(&["y0", "it's"]), "['y0', \"it's\"]");
        assert_eq!(repr_list::<&str>(&[]), "[]");
    }

    #[test]
    fn unknown_verbs_and_missing_options_are_usage_errors() {
        let (rc, _, err) = testutil::run(&["paint"]);
        assert_eq!(rc, 2);
        assert!(err.contains("invalid choice: 'paint'"), "{err}");
        let (rc, _, err) = testutil::run(&["order", "--pop", "x"]);
        assert_eq!(rc, 2);
        assert!(
            err.contains("the following arguments are required: --q"),
            "{err}"
        );
        let (rc, _, err) =
            testutil::run(&["chroms", "--counts", "c", "--maps", "m", "--min-cm", "x"]);
        assert_eq!(rc, 2);
        assert!(
            err.contains("argument --min-cm: invalid float value: 'x'"),
            "{err}"
        );
        let (rc, out, _) = testutil::run(&["report", "--help"]);
        assert_eq!(rc, 0);
        assert!(out.starts_with("usage: og ancestry-fine report"));
    }
}
