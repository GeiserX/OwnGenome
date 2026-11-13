//! Step 41: the red cell, platelet and neutrophil antigens one genome
//! predicts (replaces `bin/bloodtype_report.py`).
//!
//! `og bloodtype report` reads the sample's genotype at each defining
//! variant of `assets/bloodtype_markers.tsv` and writes one row per blood
//! group system: ABO (with A1 and A2), RhD, RhCE, Kell, Duffy, Kidd, MNS,
//! Lutheran, Diego, Dombrock, Colton, Vel, Langereis, Junior, H (Bombay and
//! para-Bombay), secretor status, the platelet antigens HPA-1 to -35 and the
//! neutrophil antigens HNA-1 (a, b, c), HNA-3, HNA-4 and HNA-5. Columns:
//! system, phenotype, genotype, confidence, note. The systems where the
//! person lacks a high-prevalence antigen come first. It is a prediction
//! from the DNA, not a serological typing. The evidence file has one row per
//! site and per RhD depth region: what the VCF and the reads said.
//!
//! - [`markers`]: the marker table and the RhD depth regions (the RHD row
//!   of `bin/cn_loci.tsv`).
//! - [`genotypes`]: the VCF at the marker sites, through its index when it
//!   has one, and the `samtools bedcov` depth.
//! - [`rules`]: one row per system.
//!
//! `og bloodtype check-markers` checks that the marker table and the rules
//! agree; `og bloodtype systems` prints the systems in their base order,
//! before the report moves the ones that lack a high-prevalence antigen to
//! the top.
//!
//! What changes on purpose: the VCF is read through its `.tbi` or `.csi`
//! (318 positions, not the whole file), the catalog is read and the RHD row
//! assessed once, a rule site missing from the table is an `ERROR:` line,
//! not a traceback, and each output is written whole or not at all.

pub mod genotypes;
pub mod markers;
pub mod rules;

use og_cn::catalog::Catalog;
use og_filters::args::{self, Spec};
use std::collections::{BTreeSet, HashMap, HashSet};
use std::fmt;
use std::io::{self, Write};

/// The table's columns.
pub const TABLE_HEADER: [&str; 5] = ["system", "phenotype", "genotype", "confidence", "note"];
/// The evidence table's columns.
pub const EVIDENCE_HEADER: [&str; 12] = [
    "system",
    "name",
    "variant",
    "rsid",
    "chrom",
    "pos",
    "ref",
    "alt",
    "status",
    "alt_copies",
    "depth",
    "note",
];

/// A problem with an input or an output: printed as `ERROR: <text>`, exit 1.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Error(pub String);

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

impl Error {
    /// An I/O error on `path`, worded as Python's `OSError`.
    pub fn io(e: &io::Error, path: &str) -> Error {
        if e.kind() == io::ErrorKind::InvalidData {
            Error(format!("{path}: {e}"))
        } else {
            Error(og_filters::os_error(e, path))
        }
    }
}

impl From<og_cn::Error> for Error {
    fn from(e: og_cn::Error) -> Error {
        Error(e.0)
    }
}

/// The inputs of one report.
#[derive(Debug, Clone, Default)]
pub struct Files {
    pub vcf: String,
    pub markers: String,
    /// `samtools bedcov -Q 20 -j` over the sites and the RHD regions.
    pub depth: Option<String>,
    /// The copy-number catalog; `None` is the one beside og.
    pub cn_loci: Option<String>,
    pub out: String,
    pub evidence: String,
}

fn catalog(path: Option<&str>) -> Result<Catalog, Error> {
    match path {
        Some(p) => Ok(Catalog::load(p)?),
        None => {
            let p = og_cn::catalog::default_path()
                .ok_or_else(|| Error("no bin/cn_loci.tsv beside og; give --cn-loci".to_owned()))?;
            Ok(Catalog::load(&p.to_string_lossy())?)
        }
    }
}

/// A TSV body as Python's `csv.writer` (tab, QUOTE_MINIMAL, `\n`) writes it.
fn tsv<R: AsRef<str>>(header: &[&str], rows: impl IntoIterator<Item = Vec<R>>) -> Vec<u8> {
    let mut w = pycompat::csv::Writer::new(Vec::new(), '\t');
    // Writing to a Vec cannot fail.
    let _ = w.write_row(header);
    for r in rows {
        let _ = w.write_row(r);
    }
    w.into_inner()
}

/// Writes `body` to `path` whole, through a temporary file and a rename.
fn write_atomic(path: &str, body: &[u8]) -> Result<(), Error> {
    let fail = |e: io::Error| Error(format!("cannot write {path}: {e}"));
    let mut f = og_io::AtomicFile::create(path).map_err(fail)?;
    f.write_all(body).map_err(fail)?;
    f.commit().map_err(fail)
}

/// Runs the report: writes the table and the evidence, and returns the rows
/// and how the VCF was read.
pub fn report(files: &Files) -> Result<(Vec<rules::Row>, genotypes::Access), Error> {
    let cat = catalog(files.cn_loci.as_deref())?;
    let rhd_row = markers::rhd_row(&cat)?;
    let regions = markers::rhd_regions(rhd_row);
    let sites = markers::read_markers(&files.markers)?;
    let (mut sites, access) = genotypes::read_vcf(&files.vcf, sites)?;
    let depth: Option<HashMap<String, f64>> = match &files.depth {
        Some(p) => Some(genotypes::read_depth(p)?),
        None => None,
    };
    genotypes::apply_depth(&mut sites, depth.as_ref());
    let rhd = rules::rhd_depth(rhd_row, &regions, depth.as_ref());
    let rows = rules::predict(&sites, &rhd, rhd_row)?;
    let table = tsv(
        &TABLE_HEADER,
        rows.iter().map(|r| {
            vec![
                r.system.as_str(),
                r.phenotype.as_str(),
                r.genotype.as_str(),
                r.confidence,
                r.note.as_str(),
            ]
        }),
    );
    let evidence = tsv(
        &EVIDENCE_HEADER,
        rules::evidence_rows(&sites, &regions, depth.as_ref(), &rhd)
            .into_iter()
            .map(Vec::from),
    );
    write_atomic(&files.out, &table)?;
    write_atomic(&files.evidence, &evidence)?;
    Ok((rows, access))
}

/// The disagreements between the marker table at `path` (with the RhD
/// depth regions of the catalog) and the rules, one sentence each.
pub fn check_markers(path: &str, cn_loci: Option<&str>) -> Result<Vec<String>, Error> {
    let cat = catalog(cn_loci)?;
    let regions = markers::rhd_regions(markers::rhd_row(&cat)?);
    let sites = markers::read_markers(path)?;
    let mut problems = Vec::new();
    let mut fail = |what: &str, names: BTreeSet<String>| {
        if !names.is_empty() {
            let list: Vec<String> = names.into_iter().collect();
            problems.push(format!("{what}: {}", list.join(", ")));
        }
    };
    let nulls: BTreeSet<&str> = sites
        .values()
        .filter(|m| m.is_null())
        .map(|m| m.name.as_str())
        .collect();
    let ruled: BTreeSet<&str> = rules::PAIRS
        .iter()
        .map(|(n, _)| n.as_str())
        .chain(rules::SPECIAL_SITES)
        .chain(nulls.iter().copied())
        .collect();
    let owned = |it: &mut dyn Iterator<Item = &str>| it.map(str::to_owned).collect::<BTreeSet<_>>();
    fail(
        "sites no rule reads",
        owned(
            &mut sites
                .keys()
                .map(String::as_str)
                .filter(|n| !ruled.contains(n)),
        ),
    );
    fail(
        "rule sites the table lacks",
        owned(&mut ruled.iter().copied().filter(|n| !sites.contains_key(*n))),
    );
    fail(
        "phenotype templates that do not name both antigens of their site",
        owned(&mut rules::PAIRS.iter().filter_map(|(n, p)| {
            let m = sites.get(n)?;
            let both = p.template.contains(&format!("{{{}}}", m.ref_means))
                && p.template.contains(&format!("{{{}}}", m.alt_means));
            (!both).then_some(n.as_str())
        })),
    );
    let systems: HashSet<&str> = sites.values().map(|m| m.system.as_str()).collect();
    fail(
        "systems with no site",
        owned(
            &mut rules::SYSTEMS
                .iter()
                .map(|(s, _)| s.as_str())
                .filter(|s| !systems.contains(s)),
        ),
    );
    let mut seen: HashMap<(&str, i64), &str> = HashMap::new();
    let mut shared = BTreeSet::new();
    for m in sites.values() {
        if let Some(first) = seen.insert((m.chrom.as_str(), m.pos), m.name.as_str()) {
            shared.insert(format!("{first} and {}", m.name));
        }
    }
    fail("sites at one position", shared);
    let want: BTreeSet<String> = ["RHD_FLANK_UP", "RHD_FLANK_DOWN", "RHD_GENE"]
        .into_iter()
        .map(str::to_owned)
        .chain((1..=10).map(|e| format!("RHD_EX{e}")))
        .collect();
    fail(
        "RhD depth regions the RHD row of the catalog lacks",
        want.iter()
            .filter(|n| !regions.contains_key(*n))
            .cloned()
            .collect(),
    );
    fail(
        "systems read from null rows with none",
        owned(
            &mut rules::NULL_SYSTEMS
                .iter()
                .map(|s| s.system)
                .filter(|s| !sites.values().any(|m| m.is_null() && m.system == *s)),
        ),
    );
    fail(
        "Kell or Kidd null rows with no null_on",
        owned(&mut sites.values().filter_map(|m| {
            (m.is_null() && rules::NULL_HOSTS.contains(&m.system.as_str()) && m.null_on.is_empty())
                .then_some(m.name.as_str())
        })),
    );
    let mut labels: HashMap<&str, HashSet<&str>> = HashMap::new();
    for m in sites.values().filter(|m| !m.is_null()) {
        labels
            .entry(m.system.as_str())
            .or_default()
            .extend([m.ref_means.as_str(), m.alt_means.as_str()]);
    }
    fail(
        "null_on alleles that are no antigen of the row's system",
        owned(&mut sites.values().filter_map(|m| {
            let own = labels.get(m.system.as_str());
            let off = m
                .null_on
                .split_whitespace()
                .any(|a| !own.is_some_and(|l| l.contains(a)));
            (m.is_null() && off).then_some(m.name.as_str())
        })),
    );
    fail(
        "high_prevalence antigens the row does not give",
        owned(&mut sites.values().filter_map(|m| {
            let h = m.high_prevalence.as_str();
            (!h.is_empty() && h != m.ref_means && h != m.alt_means).then_some(m.name.as_str())
        })),
    );
    Ok(problems)
}

const REPORT: Spec = Spec {
    prog: "og bloodtype report",
    usage: "usage: og bloodtype report --sample SAMPLE --vcf VCF --markers MARKERS [--depth DEPTH] [--cn-loci CN_LOCI] --out OUT --evidence EVIDENCE\n",
    options: &[
        "--sample",
        "--vcf",
        "--markers",
        "--depth",
        "--cn-loci",
        "--out",
        "--evidence",
    ],
    positionals: &[],
    required: &["--sample", "--vcf", "--markers", "--out", "--evidence"],
};

const CHECK: Spec = Spec {
    prog: "og bloodtype check-markers",
    usage: "usage: og bloodtype check-markers MARKERS [--cn-loci CN_LOCI]\n",
    options: &["--cn-loci"],
    positionals: &["markers"],
    required: &[],
};

const SYSTEMS: Spec = Spec {
    prog: "og bloodtype systems",
    usage: "usage: og bloodtype systems\n",
    options: &[],
    positionals: &[],
    required: &[],
};

const USAGE: &str = "usage: og bloodtype report|check-markers|systems [args]\n";

/// `og bloodtype VERB ARGS`: the exit code (0, 1 for a problem, 2 for a
/// refused command line).
pub fn cli(verb: &str, rest: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let spec = match verb {
        "report" => &REPORT,
        "check-markers" => &CHECK,
        "systems" => &SYSTEMS,
        "-h" | "--help" => {
            let _ = out.write_all(USAGE.as_bytes());
            return 0;
        }
        other => {
            let _ = writeln!(err, "{USAGE}og bloodtype: error: unknown command '{other}'");
            return 2;
        }
    };
    let a = match args::parse(spec, rest) {
        Ok(a) => a,
        Err(e) => return args::report(spec, &e),
    };
    let opt = |name: &str| a.get(name).unwrap_or_default().to_owned();
    let done = match verb {
        "systems" => {
            let mut text = String::new();
            for (s, _) in rules::SYSTEMS.iter() {
                text.push_str(s);
                text.push('\n');
            }
            out.write_all(text.as_bytes())
                .map_err(|e| Error(format!("cannot write stdout: {e}")))
        }
        "check-markers" => match check_markers(&a.pos[0], a.get("--cn-loci")) {
            Ok(problems) => {
                for p in &problems {
                    let _ = writeln!(err, "{}: {p}", a.pos[0]);
                }
                return u8::from(!problems.is_empty());
            }
            Err(e) => Err(e),
        },
        _ => {
            let files = Files {
                vcf: opt("--vcf"),
                markers: opt("--markers"),
                depth: a.get("--depth").map(str::to_owned),
                cn_loci: a.get("--cn-loci").map(str::to_owned),
                out: opt("--out"),
                evidence: opt("--evidence"),
            };
            report(&files).and_then(|(rows, _)| {
                let mut text = format!(
                    "Blood groups of {}, predicted from the genome (not a serological typing):\n",
                    opt("--sample")
                );
                for r in &rows {
                    text.push_str(&format!(
                        "  {:<9} {:<34} {}\n",
                        r.system, r.phenotype, r.confidence
                    ));
                }
                out.write_all(text.as_bytes())
                    .map_err(|e| Error(format!("cannot write stdout: {e}")))
            })
        }
    };
    match done {
        Ok(()) => 0,
        Err(e) => {
            let _ = writeln!(err, "ERROR: {e}");
            1
        }
    }
}
