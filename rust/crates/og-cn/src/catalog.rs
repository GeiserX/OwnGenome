//! `bin/cn_loci.tsv`: one row per locus, read and checked once per process.
//!
//! The checks and their messages are those of `copy_number.py`
//! (`parse_row`, `load_catalog`), plus one it lacked: a `more=` reading
//! needs a numbered reading to count past (the Python failed on `max()` of
//! an empty sequence when a call reached it, naming no row).

use crate::{Error, py_float, py_int};
use indexmap::IndexMap;
use og_filters::py_repr;
use std::collections::HashSet;
use std::path::{Path, PathBuf};

/// The catalog's columns, in order.
pub const COLUMNS: [&str; 22] = [
    "locus",
    "kind",
    "chrom",
    "targets",
    "controls",
    "details",
    "target_depth",
    "control_depth",
    "scale",
    "offset",
    "bedcov",
    "mapq0_min",
    "mapq0_max",
    "min_control_depth",
    "max_control_ratio",
    "call",
    "interval",
    "readings",
    "input",
    "owner",
    "cross_check",
    "pmids",
];

/// The findings a `call psv` row must have a reading for.
pub const PSV_FINDINGS: [&str; 4] = ["none", "deletion", "conversion", "pseudogene_deletion"];

/// What a row counts.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kind {
    Gene,
    Interval,
    Organelle,
}

/// Where a row's depths come from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Input {
    Bedcov,
    CnvpytorLevel,
}

/// How the target segments combine into T.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TargetMode {
    /// Added up: each segment is one reference copy.
    Sum,
    /// Length-weighted mean.
    Pooled,
}

/// How the control segments combine into C.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ControlMode {
    /// Plain mean of the segments.
    Flanks,
    /// Length-weighted mean.
    Pooled,
}

/// The `call` column.
#[derive(Debug, Clone, PartialEq)]
pub enum CallRule {
    /// The nearest whole number, when the estimate lies within `max_off` of it.
    Round { max_off: f64 },
    /// The number of cut points the estimate reaches.
    Cut(Vec<f64>),
    /// low, typical or high by the 95% interval against two bounds.
    Class { low: f64, high: f64 },
    /// The module total as `Round`, then the PSV sites.
    Psv { max_off: f64 },
}

/// One `NAME:start-end[:label]` segment, 0-based half-open.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Segment {
    pub name: String,
    pub start: i64,
    pub end: i64,
    pub label: String,
}

/// One paralogous sequence variant site of a `call psv` row.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PsvSite {
    pub name: String,
    /// 0-based position on the gene.
    pub gene_pos: i64,
    /// 0-based position on the pseudogene.
    pub pseudo_pos: i64,
    /// The gene's bases there (`b"AC"` for In2G).
    pub gene: Vec<u8>,
    /// The pseudogene's bases there.
    pub pseudo: Vec<u8>,
    pub text: String,
}

/// The depth-ratio settings of a `bedcov` row.
#[derive(Debug, Clone, PartialEq)]
pub struct Depth {
    pub target_mode: TargetMode,
    pub control_mode: ControlMode,
    pub scale: f64,
    pub offset: f64,
    /// The `bedcov` passes, as written (`.` for the plain pass).
    pub passes: Vec<String>,
    pub mapq0_min: Option<f64>,
    pub mapq0_max: Option<f64>,
    pub min_control: f64,
    pub max_ratio: Option<f64>,
    /// A 95% interval from the per-segment spread (`interval cv`).
    pub ci_cv: bool,
}

/// One catalog row.
#[derive(Debug, Clone, PartialEq)]
pub struct Locus {
    pub id: String,
    pub kind: Kind,
    pub chrom: String,
    pub input: Input,
    pub targets: Vec<Segment>,
    pub controls: Vec<Segment>,
    pub details: Vec<Segment>,
    pub call: CallRule,
    /// `KEY=text` readings, in the row's order.
    pub readings: IndexMap<String, String>,
    pub psv: Vec<PsvSite>,
    /// `target_depth` as written (a `call psv` row needs `sum` whatever its input).
    pub target_depth: String,
    /// Present for a `bedcov` row.
    pub depth: Option<Depth>,
}

impl Locus {
    /// The bedcov settings; every row `og cn` counts has them.
    pub fn depth(&self) -> &Depth {
        self.depth
            .as_ref()
            .expect("a bedcov row has its depth settings")
    }

    /// `[(chrom, start, end, name)]` of the row's segments, sorted by
    /// position (controls, targets, details; stable).
    pub fn regions(&self) -> Vec<(&str, i64, i64, &str)> {
        let mut out: Vec<_> = self
            .controls
            .iter()
            .chain(&self.targets)
            .chain(&self.details)
            .map(|s| (self.chrom.as_str(), s.start, s.end, s.name.as_str()))
            .collect();
        out.sort_by_key(|r| (r.1, r.2));
        out
    }

    /// The row's reading for a count: its own line, else `more=` past the
    /// highest numbered line, else empty.
    pub fn reading_count(&self, n: i64) -> String {
        if let Some(t) = self.readings.get(&n.to_string()) {
            return t.clone();
        }
        match (self.readings.get("more"), self.max_numbered()) {
            (Some(more), Some(max)) if n > max => more.replace("{n}", &n.to_string()),
            _ => String::new(),
        }
    }

    /// The row's reading for a word key (a class or a finding), or empty.
    pub fn reading_word(&self, key: &str) -> String {
        self.readings.get(key).cloned().unwrap_or_default()
    }

    fn max_numbered(&self) -> Option<i64> {
        self.readings
            .keys()
            .filter(|k| !k.is_empty() && k.bytes().all(|c| c.is_ascii_digit()))
            .filter_map(|k| k.parse().ok())
            .max()
    }

    /// An interval row `og cn intervals` calls: kind interval, bedcov, one
    /// plain pass (LCE3BC, with a MAPQ >= 1 pass, stays a gene-dosage row).
    pub fn is_interval_row(&self) -> bool {
        self.kind == Kind::Interval && self.input == Input::Bedcov && self.depth().passes == ["."]
    }
}

/// The catalog, rows in file order.
#[derive(Debug, Clone)]
pub struct Catalog {
    pub path: String,
    pub loci: IndexMap<String, Locus>,
}

impl Catalog {
    /// The row named `id`, when it is a bedcov row.
    pub fn bedcov(&self, id: &str) -> Option<&Locus> {
        self.loci.get(id).filter(|l| l.input == Input::Bedcov)
    }

    /// Reads and checks the catalog at `path`.
    pub fn load(path: &str) -> Result<Catalog, Error> {
        let lines = crate::read_lines_strict(path)?;
        Catalog::parse(path, lines.iter().map(String::as_str))
    }

    /// Checks the catalog's lines (each without its `\n`); `path` names it
    /// in messages.
    pub fn parse<'a>(path: &str, lines: impl Iterator<Item = &'a str>) -> Result<Catalog, Error> {
        let mut rows = lines.filter(|l| !crate::strip(l).is_empty() && !l.starts_with('#'));
        let header_ok = rows
            .next()
            .is_some_and(|h| h.split('\t').eq(COLUMNS.iter().copied()));
        let mut loci = IndexMap::new();
        let mut seen: HashSet<String> = HashSet::new();
        let mut any = false;
        for line in rows {
            if !header_ok {
                break;
            }
            any = true;
            let cells: Vec<&str> = line.split('\t').collect();
            if cells.len() != COLUMNS.len() {
                return Err(Error(format!(
                    "{path}: row {} has the wrong number of columns",
                    py_repr(cells[0])
                )));
            }
            if loci.contains_key(cells[0]) {
                return Err(Error(format!("{path}: locus {} is listed twice", cells[0])));
            }
            let row = parse_row(&cells)?;
            for s in row.targets.iter().chain(&row.controls).chain(&row.details) {
                if !seen.insert(s.name.clone()) {
                    return Err(Error(format!(
                        "{path}: segment name {} is used by two loci",
                        s.name
                    )));
                }
            }
            loci.insert(row.id.clone(), row);
        }
        if !header_ok || !any {
            return Err(Error(format!(
                "{path}: the header must be {}",
                COLUMNS.join(" ")
            )));
        }
        Ok(Catalog {
            path: path.to_owned(),
            loci,
        })
    }
}

/// The catalog beside the running binary: the first folder above it that
/// holds `cn_loci.tsv` or `bin/cn_loci.tsv` (`bin/.og/og-linux-*` in a
/// checkout, `/pgp-bin/.og/...` in a container, `rust/target/*/og`).
pub fn default_path() -> Option<PathBuf> {
    let exe = std::env::current_exe().ok()?;
    exe.ancestors().skip(1).find_map(|d| {
        [d.join("cn_loci.tsv"), d.join("bin").join("cn_loci.tsv")]
            .into_iter()
            .find(|p| Path::is_file(p))
    })
}

fn cell<'a>(cells: &[&'a str], col: &str) -> &'a str {
    cells[COLUMNS.iter().position(|c| *c == col).expect("a column")]
}

/// A number column: `.` gives `default`; anything not a finite number is
/// refused.
fn num(cells: &[&str], col: &str, default: Option<f64>) -> Result<Option<f64>, Error> {
    let v = cell(cells, col);
    if v == "." {
        return Ok(default);
    }
    match pycompat::num::parse_float(v) {
        Ok(x) if x.is_finite() => Ok(Some(x)),
        _ => Err(Error(format!(
            "{}: {col} is {}, not a number or '.'",
            cells[0],
            py_repr(v)
        ))),
    }
}

fn segments(cells: &[&str], col: &str) -> Result<Vec<Segment>, Error> {
    let v = cell(cells, col);
    if v == "." {
        return Ok(Vec::new());
    }
    let loc = cells[0];
    let mut out = Vec::new();
    for part in v.split(';') {
        let mut bits = part.splitn(3, ':');
        let name = bits.next().unwrap_or("");
        let span = bits.next();
        let label = bits.next().unwrap_or("");
        let bad = || {
            Error(format!(
                "{loc}: {col} segment {} is not NAME:start-end[:label]",
                py_repr(part)
            ))
        };
        let span: Vec<&str> = span.ok_or_else(bad)?.split('-').collect();
        let [s, e] = span[..] else { return Err(bad()) };
        let (start, end) = match (py_int(s), py_int(e)) {
            (Ok(s), Ok(e)) => (s, e),
            _ => return Err(bad()),
        };
        if name.is_empty() || end <= start || start < 0 {
            return Err(Error(format!(
                "{loc}: {col} segment {} has no name or an empty span",
                py_repr(part)
            )));
        }
        out.push(Segment {
            name: name.to_owned(),
            start,
            end,
            label: label.to_owned(),
        });
    }
    Ok(out)
}

fn readings(cells: &[&str]) -> Result<IndexMap<String, String>, Error> {
    let v = cell(cells, "readings");
    let mut out = IndexMap::new();
    if v == "." {
        return Ok(out);
    }
    for part in v.split('|') {
        match part.split_once('=') {
            Some((k, t)) if !k.is_empty() && !t.is_empty() => {
                out.insert(k.to_owned(), t.to_owned());
            }
            _ => {
                return Err(Error(format!(
                    "{}: reading {} is not KEY=text",
                    cells[0],
                    py_repr(part)
                )));
            }
        }
    }
    Ok(out)
}

fn call(cells: &[&str]) -> Result<CallRule, Error> {
    let raw = cell(cells, "call");
    let bad = || {
        Error(format!(
            "{}: call {} is not round:D, cut:A,B..., class:LOW,HIGH or psv:D",
            cells[0],
            py_repr(raw)
        ))
    };
    let (mode, args) = raw.split_once(':').ok_or_else(bad)?;
    let vals: Vec<f64> = args
        .split(',')
        .map(|x| py_float(x).ok().filter(|v| v.is_finite()))
        .collect::<Option<_>>()
        .unwrap_or_default();
    if vals.is_empty() || vals.windows(2).any(|w| w[0] > w[1]) {
        return Err(bad());
    }
    match (mode, &vals[..]) {
        ("round", [d]) => Ok(CallRule::Round { max_off: *d }),
        ("psv", [d]) => Ok(CallRule::Psv { max_off: *d }),
        ("class", [lo, hi]) => Ok(CallRule::Class {
            low: *lo,
            high: *hi,
        }),
        ("cut", _) => Ok(CallRule::Cut(vals)),
        _ => Err(bad()),
    }
}

fn psv_sites(loc: &str, details: &[Segment]) -> Result<Vec<PsvSite>, Error> {
    let mut out = Vec::new();
    for seg in details {
        let label = &seg.label;
        let bits: Vec<&str> = label.splitn(3, ' ').collect();
        let parsed = py_int(bits[0])
            .ok()
            .and_then(|pos| bits.get(1).map(|b| (pos, b)));
        let Some((pos, alleles)) = parsed else {
            return Err(Error(format!(
                "{loc}: PSV {} label {} is not 'PSEUDO_POS GENE>PSEUDO text'",
                seg.name,
                py_repr(label)
            )));
        };
        let (gene, pseudo) = alleles.split_once('>').unwrap_or((alleles, ""));
        let acgt = |s: &str| s.bytes().all(|c| b"ACGT".contains(&c));
        let shared = gene.bytes().any(|c| pseudo.as_bytes().contains(&c));
        if seg.end - seg.start != 1
            || pos < 0
            || gene.is_empty()
            || pseudo.is_empty()
            || !acgt(gene)
            || !acgt(pseudo)
            || shared
            || bits.len() < 3
        {
            return Err(Error(format!(
                "{loc}: PSV {} must be one base with distinct gene and pseudogene bases and a text (got {}-{}, {})",
                seg.name,
                seg.start,
                seg.end,
                py_repr(label)
            )));
        }
        out.push(PsvSite {
            name: seg.name.clone(),
            gene_pos: seg.start,
            pseudo_pos: pos,
            gene: unique(gene),
            pseudo: unique(pseudo),
            text: bits[2].to_owned(),
        });
    }
    Ok(out)
}

/// The distinct bases of `s`, in order (the Python kept them as a set).
fn unique(s: &str) -> Vec<u8> {
    let mut out: Vec<u8> = Vec::new();
    for c in s.bytes() {
        if !out.contains(&c) {
            out.push(c);
        }
    }
    out
}

/// Checks one row (its 22 cells) as `copy_number.parse_row` does, in the
/// same order, so the first problem is the one the Python named.
pub fn parse_row(cells: &[&str]) -> Result<Locus, Error> {
    let loc = cells[0];
    let kind = match cell(cells, "kind") {
        "gene" => Kind::Gene,
        "interval" => Kind::Interval,
        "organelle" => Kind::Organelle,
        k => {
            return Err(Error(format!(
                "{loc}: kind {} (one of gene, interval, organelle)",
                py_repr(k)
            )));
        }
    };
    let input = match cell(cells, "input") {
        "bedcov" => Input::Bedcov,
        "cnvpytor-level" => Input::CnvpytorLevel,
        i => {
            return Err(Error(format!(
                "{loc}: input {} (one of bedcov, cnvpytor-level)",
                py_repr(i)
            )));
        }
    };
    let chrom = cell(cells, "chrom");
    if !is_chrom(chrom) {
        return Err(Error(format!(
            "{loc}: chrom {} (GRCh38 names, chr1 to chrM)",
            py_repr(chrom)
        )));
    }
    let targets = segments(cells, "targets")?;
    let controls = segments(cells, "controls")?;
    let details = segments(cells, "details")?;
    let call = call(cells)?;
    let readings = readings(cells)?;
    let is_psv = matches!(call, CallRule::Psv { .. });
    let psv = if is_psv {
        psv_sites(loc, &details)?
    } else {
        Vec::new()
    };
    let target_depth = cell(cells, "target_depth");
    if is_psv
        && (psv.is_empty()
            || target_depth != "sum"
            || !PSV_FINDINGS.iter().all(|f| readings.contains_key(*f)))
    {
        return Err(Error(format!(
            "{loc}: call psv needs PSV sites in details, target_depth sum (gene and pseudogene) and a reading for each of {}",
            PSV_FINDINGS.join(", ")
        )));
    }
    if targets.is_empty() {
        return Err(Error(format!("{loc}: no target segment")));
    }
    let depth = if input == Input::Bedcov {
        Some(depth(cells, &controls, &call)?)
    } else {
        None
    };
    let mut names = HashSet::new();
    if !targets
        .iter()
        .chain(&controls)
        .chain(&details)
        .all(|s| names.insert(s.name.as_str()))
    {
        return Err(Error(format!("{loc}: a segment name is used twice")));
    }
    let row = Locus {
        id: loc.to_owned(),
        kind,
        chrom: chrom.to_owned(),
        input,
        targets,
        controls,
        details,
        call,
        readings,
        psv,
        target_depth: target_depth.to_owned(),
        depth,
    };
    if row.readings.contains_key("more") && row.max_numbered().is_none() {
        return Err(Error(format!(
            "{loc}: reading more= needs a numbered reading (0=, 1=, ...) to count past"
        )));
    }
    Ok(row)
}

fn depth(cells: &[&str], controls: &[Segment], call: &CallRule) -> Result<Depth, Error> {
    let loc = cells[0];
    if controls.is_empty() {
        return Err(Error(format!("{loc}: a bedcov row needs control segments")));
    }
    let target_mode = match cell(cells, "target_depth") {
        "sum" => Some(TargetMode::Sum),
        "pooled" => Some(TargetMode::Pooled),
        _ => None,
    };
    let control_mode = match cell(cells, "control_depth") {
        "flanks" => Some(ControlMode::Flanks),
        "pooled" => Some(ControlMode::Pooled),
        _ => None,
    };
    let (Some(target_mode), Some(control_mode)) = (target_mode, control_mode) else {
        return Err(Error(format!(
            "{loc}: target_depth sum|pooled and control_depth flanks|pooled"
        )));
    };
    let scale = num(cells, "scale", None)?;
    let offset = num(cells, "offset", Some(0.0))?.unwrap_or(0.0);
    let scale = match scale {
        Some(s) if s > 0.0 => s,
        _ => return Err(Error(format!("{loc}: scale must be a positive number"))),
    };
    let min_control = num(cells, "min_control_depth", Some(0.0))?.unwrap_or(0.0);
    let max_ratio = num(cells, "max_control_ratio", None)?;
    let mapq0_min = num(cells, "mapq0_min", None)?;
    let mapq0_max = num(cells, "mapq0_max", None)?;
    let passes: Vec<String> = cell(cells, "bedcov")
        .split(',')
        .map(str::to_owned)
        .collect();
    if (mapq0_min.is_some() || mapq0_max.is_some()) && passes.len() != 2 {
        return Err(Error(format!(
            "{loc}: a MAPQ 0 guard needs two bedcov passes (all reads, MAPQ >= 1)"
        )));
    }
    let interval = cell(cells, "interval");
    if matches!(call, CallRule::Class { .. }) && interval != "cv" {
        return Err(Error(format!("{loc}: call class needs interval cv")));
    }
    if interval != "cv" && interval != "." {
        return Err(Error(format!(
            "{loc}: interval {} (cv or '.')",
            py_repr(interval)
        )));
    }
    Ok(Depth {
        target_mode,
        control_mode,
        scale,
        offset,
        passes,
        mapq0_min,
        mapq0_max,
        min_control,
        max_ratio,
        ci_cv: interval == "cv",
    })
}

/// A GRCh38 primary chromosome name: chr1 to chr22, chrX, chrY, chrM.
fn is_chrom(c: &str) -> bool {
    match c.strip_prefix("chr") {
        Some("X" | "Y" | "M") => true,
        Some(n) => n
            .parse::<u8>()
            .is_ok_and(|v| (1..=22).contains(&v) && n == v.to_string()),
        None => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn row(over: &[(&str, &str)]) -> Vec<String> {
        let mut base: Vec<(&str, String)> = COLUMNS
            .iter()
            .map(|c| {
                let v = match *c {
                    "locus" => "G",
                    "kind" => "gene",
                    "chrom" => "chr1",
                    "targets" => "G_A:1000-2000",
                    "controls" => "G_L:0-500;G_R:3000-3500",
                    "target_depth" => "sum",
                    "control_depth" => "flanks",
                    "scale" => "2",
                    "offset" => "0",
                    "bedcov" => ".",
                    "min_control_depth" => "8",
                    "call" => "round:0.3",
                    "input" => "bedcov",
                    _ => ".",
                };
                (*c, v.to_owned())
            })
            .collect();
        for (k, v) in over {
            base.iter_mut().find(|(c, _)| c == k).unwrap().1 = (*v).to_owned();
        }
        base.into_iter().map(|(_, v)| v).collect()
    }

    fn parse(over: &[(&str, &str)]) -> Result<Locus, Error> {
        let r = row(over);
        let cells: Vec<&str> = r.iter().map(String::as_str).collect();
        parse_row(&cells)
    }

    #[test]
    fn a_plain_row_loads() {
        let l = parse(&[]).unwrap();
        assert_eq!(l.depth().scale, 2.0);
        assert_eq!(l.regions()[0], ("chr1", 0, 500, "G_L"));
    }

    #[test]
    fn broken_rows_are_refused_with_the_python_text() {
        let cases: &[(&[(&str, &str)], &str)] = &[
            (
                &[("kind", "cnv")],
                "G: kind 'cnv' (one of gene, interval, organelle)",
            ),
            (
                &[("chrom", "chrx")],
                "G: chrom 'chrx' (GRCh38 names, chr1 to chrM)",
            ),
            (
                &[("targets", "G_A:5-5")],
                "G: targets segment 'G_A:5-5' has no name or an empty span",
            ),
            (
                &[("targets", "G_A:5")],
                "G: targets segment 'G_A:5' is not NAME:start-end[:label]",
            ),
            (
                &[("call", "class:30,45")],
                "G: call class needs interval cv",
            ),
            (
                &[("call", "cut:1.5,0.5")],
                "G: call 'cut:1.5,0.5' is not round:D, cut:A,B..., class:LOW,HIGH or psv:D",
            ),
            (
                &[("controls", "G_A:0-500")],
                "G: a segment name is used twice",
            ),
            (&[("scale", "0")], "G: scale must be a positive number"),
            (
                &[("scale", "nan")],
                "G: scale is 'nan', not a number or '.'",
            ),
            (
                &[("mapq0_min", "0.02")],
                "G: a MAPQ 0 guard needs two bedcov passes (all reads, MAPQ >= 1)",
            ),
            (&[("readings", "x")], "G: reading 'x' is not KEY=text"),
        ];
        for (over, want) in cases {
            assert_eq!(parse(over).unwrap_err().0, *want, "{over:?}");
        }
    }

    #[test]
    fn more_without_a_numbered_reading_is_a_catalog_error() {
        let e = parse(&[("readings", "more={n} copies")]).unwrap_err();
        assert!(e.0.contains("more= needs a numbered reading"), "{e}");
        let l = parse(&[("readings", "2=two|more={n} copies")]).unwrap();
        assert_eq!(l.reading_count(5), "5 copies");
        assert_eq!(l.reading_count(2), "two");
        assert_eq!(l.reading_count(1), "");
    }

    #[test]
    fn chrom_names() {
        for c in ["chr1", "chr22", "chrX", "chrY", "chrM"] {
            assert!(is_chrom(c), "{c}");
        }
        for c in ["chr0", "chr23", "chr01", "chrx", "1", "chrUn"] {
            assert!(!is_chrom(c), "{c}");
        }
    }

    #[test]
    fn the_repository_catalog_loads() {
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../../bin/cn_loci.tsv");
        let cat = Catalog::load(path).unwrap();
        assert!(cat.bedcov("HBA").is_some() && cat.bedcov("LPA_KIV2").is_some());
        assert_eq!(cat.loci["CYP21A2"].psv.len(), 6);
    }
}
