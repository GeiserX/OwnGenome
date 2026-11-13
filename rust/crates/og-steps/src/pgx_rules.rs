//! Step 36's rules (`docs/36-pgx-consensus.md`): what reaches PharmCAT as
//! an outside call, and step 08's T1K and OptiType concordance table.
//!
//! ```text
//! og pgx outside-calls --calls OUT_CALLS.tsv --consensus OUT_CONSENSUS.tsv
//!                      [--hla T1K_genotype.tsv] [--optitype RESULT.tsv] [--require-optitype]
//!                      [--pypgx PYPGX_summary.tsv] [--cyrius CYRIUS.tsv]
//!                      [--depth-check CYP2D6_depth_check.tsv]
//! og pgx hla-concordance --t1k T1K_genotype.tsv --optitype RESULT.tsv --out OUT.tsv
//! ```
//!
//! PharmCAT types neither HLA nor CYP2D6 from a VCF, so this decides what
//! its outside-call file (`gene<TAB>diplotype`, PharmCAT's `-po`) holds:
//!
//! - HLA-A and HLA-B: T1K's two alleles at two fields (`*57:01`) when every
//!   allele has a quality above 0; one allele reads as homozygous. With
//!   `--optitype`, only when OptiType gives the same pair; with
//!   `--require-optitype` and no OptiType result, neither gene.
//! - CYP2D6: only when pypgx (step 32) and Cyrius (step 21) give the same
//!   diplotype and the depth check passed. One depth-based call alone can be
//!   a false whole-gene deletion that changes a prescription.
//!
//! The consensus table (`Gene Result Outside_call Reason Evidence`) says
//! for each of the three genes what each caller said and what was passed.
//! The two files are written together or not at all: both go to temporary
//! files first, and when the second rename fails the first file is removed,
//! so an outside call never stands without the consensus row behind it.
//!
//! This is `bin/pgx_outside_calls.py`'s outside calls and concordance table
//! rewritten; `tests/parity` holds it to the Python's bytes.

use std::fmt::Write as _;
use std::fs::File;
use std::io::{self, BufReader, Write};

use og_io::{Lines, Utf8};
use pycompat::csv::{DictReader, Quoting, Reader};

/// The consensus table's columns.
pub const CONSENSUS_HEADER: [&str; 5] = ["Gene", "Result", "Outside_call", "Reason", "Evidence"];
/// The concordance table's columns.
pub const CONCORDANCE_HEADER: [&str; 5] = ["Gene", "T1K", "OptiType", "Agree", "Note"];
/// The HLA genes step 36 can pass to PharmCAT.
pub const HLA_GENES: [&str; 2] = ["HLA-A", "HLA-B"];
/// The class I genes OptiType types.
pub const OPTITYPE_GENES: [&str; 3] = ["HLA-A", "HLA-B", "HLA-C"];
/// Diplotype texts that are no call, compared lower-cased.
pub const NO_CALL: [&str; 9] = [
    "",
    "none",
    "n/a",
    "na",
    "failed",
    "indeterminate",
    "not called",
    "unknown",
    ".",
];
/// The reason a lone CYP2D6 caller is not passed on.
pub const ONE_CALLER: &str =
    "one caller only: a second caller must agree (Cyrius: setup.sh --cyrius)";
const NO_OPTITYPE: &str = "the second caller was asked for (HLA_QC), but there is no OptiType result: step 08 did not finish it";

/// Why a tool stopped: the message after `ERROR: `.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Error(pub String);

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl Error {
    /// An I/O error on `path`, worded as Python's `OSError`.
    pub fn io(e: &io::Error, path: &str) -> Error {
        if e.kind() == io::ErrorKind::InvalidData {
            Error(format!("'utf-8' codec can't decode {path}: {e}"))
        } else {
            Error(og_filters::os_error(e, path))
        }
    }
}

// --- reading ----------------------------------------------------------------

/// The lines of `path` as Python's `open()` iterates them: universal
/// newlines, each with its `\n` (none on a last line without one). Plain
/// text only: Python's `open` does not inflate gzip.
pub(crate) struct TextLines {
    lines: Option<Lines<BufReader<File>>>,
}

impl TextLines {
    pub(crate) fn open(path: &str, utf8: Utf8) -> Result<TextLines, Error> {
        let f = File::open(path).map_err(|e| Error::io(&e, path))?;
        if f.metadata().is_ok_and(|m| m.is_dir()) {
            return Err(Error::io(&io::Error::from_raw_os_error(21), path));
        }
        Ok(TextLines {
            lines: Some(Lines::new(BufReader::new(f), utf8)),
        })
    }
}

impl Iterator for TextLines {
    type Item = io::Result<String>;

    fn next(&mut self) -> Option<Self::Item> {
        let lines = self.lines.as_mut()?;
        let mut line = String::new();
        match lines.read_line(&mut line) {
            Ok(true) => Some(Ok(line)),
            Ok(false) => {
                self.lines = None;
                None
            }
            Err(e) => {
                self.lines = None;
                Some(Err(e))
            }
        }
    }
}

/// A csv reader error as `ERROR: {e}` words it.
fn csv_error(e: pycompat::csv::Error, path: &str) -> Error {
    match e {
        pycompat::csv::Error::Io(e) => Error::io(&e, path),
        other => Error(other.to_string()),
    }
}

/// `line.rstrip("\n")`.
fn chomp(line: &str) -> &str {
    line.strip_suffix('\n').unwrap_or(line)
}

/// Python's `s.strip()`.
pub(crate) fn strip(s: &str) -> &str {
    pycompat::text::strip(s, None)
}

/// `d.lower() in NO_CALL`.
fn no_call(d: &str) -> bool {
    NO_CALL.contains(&pycompat::text::lower(d).as_str())
}

// --- HLA --------------------------------------------------------------------

/// `'HLA-B*57:01:01:02'` cut to `'*57:01'`.
pub fn two_fields(allele: &str, gene: &str) -> String {
    let name = allele.strip_prefix(gene).unwrap_or(allele);
    let name = if name.starts_with('*') {
        name.to_owned()
    } else {
        format!("*{name}")
    };
    name.split(':').take(2).collect::<Vec<_>>().join(":")
}

/// One T1K allele: its name and quality (-1 when the quality is not a number).
#[derive(Debug, Clone, PartialEq)]
pub struct Allele {
    pub name: String,
    pub quality: f64,
}

impl Allele {
    /// T1K says to ignore an allele with quality 0 or below; a quality
    /// that is not a finite number is no better.
    fn unusable(&self) -> bool {
        !self.quality.is_finite() || self.quality <= 0.0
    }
}

/// T1K's alleles per gene, for the genes asked for. T1K's genotype TSV:
/// gene, number of distinct alleles, then allele, abundance, quality twice
/// (`.` when there is no second allele). The last row of a gene wins.
pub type T1k = Vec<(String, Vec<Allele>)>;

/// The alleles of `gene` in `t1k`, none when T1K wrote no row for it.
fn alleles_of<'a>(t1k: &'a T1k, gene: &str) -> &'a [Allele] {
    t1k.iter()
        .find(|(g, _)| g == gene)
        .map_or(&[], |(_, a)| a.as_slice())
}

/// Reads T1K's genotype TSV for `genes`.
pub fn read_t1k(path: &str, genes: &[&str]) -> Result<T1k, Error> {
    let mut out: T1k = Vec::new();
    for line in TextLines::open(path, Utf8::Strict)? {
        let line = line.map_err(|e| Error::io(&e, path))?;
        let fields: Vec<&str> = chomp(&line).split('\t').collect();
        if fields.len() < 5 || !genes.contains(&fields[0]) {
            continue;
        }
        let mut alleles = Vec::new();
        for i in [2, 5] {
            if i + 2 < fields.len() && !matches!(fields[i], "" | ".") {
                let quality = pycompat::num::parse_float(fields[i + 2]).unwrap_or(-1.0);
                alleles.push(Allele {
                    name: fields[i].to_owned(),
                    quality,
                });
            }
        }
        match out.iter_mut().find(|(g, _)| g == fields[0]) {
            Some(slot) => slot.1 = alleles,
            None => out.push((fields[0].to_owned(), alleles)),
        }
    }
    Ok(out)
}

/// OptiType's two-field calls per class I gene, in [`OPTITYPE_GENES`]
/// order. Empty for a file step 08 wrote instead (a `#` line saying why
/// OptiType did not run): no call for any gene.
pub type OptiType = Vec<(&'static str, Vec<String>)>;

/// Reads OptiType's result TSV: an index column, then `A1 A2 B1 B2 C1 C2
/// Reads Objective`, alleles as `A*01:01`.
pub fn read_optitype(path: &str) -> Result<OptiType, Error> {
    let mut lines = Vec::new();
    for line in TextLines::open(path, Utf8::Strict)? {
        let line = line.map_err(|e| Error::io(&e, path))?;
        if !strip(&line).is_empty() && !line.starts_with('#') {
            lines.push(Ok::<_, io::Error>(line));
        }
    }
    let mut rows = Vec::new();
    for row in Reader::new(lines.into_iter(), '\t', Quoting::Minimal) {
        rows.push(row.map_err(|e| csv_error(e, path))?);
        if rows.len() == 2 {
            break;
        }
    }
    if rows.len() < 2 {
        return Ok(Vec::new());
    }
    let col = |name: &str| rows[0].iter().rposition(|h| strip(h) == name);
    let mut out = Vec::new();
    for gene in OPTITYPE_GENES {
        let letter = &gene["HLA-".len()..];
        let mut calls = Vec::new();
        for k in [format!("{letter}1"), format!("{letter}2")] {
            let v = col(&k)
                .and_then(|i| rows[1].get(i))
                .map_or("", |v| strip(v));
            if !v.is_empty() && !no_call(v) {
                let v = if v.starts_with(&format!("{letter}*")) {
                    &v[letter.len()..]
                } else {
                    v
                };
                calls.push(two_fields(v, gene));
            }
        }
        if calls.len() == 1 {
            calls.push(calls[0].clone());
        }
        out.push((gene, calls));
    }
    Ok(out)
}

/// OptiType's calls of `gene`, none when it made none.
fn optitype_of<'a>(ot: &'a OptiType, gene: &str) -> &'a [String] {
    ot.iter()
        .find(|(g, _)| *g == gene)
        .map_or(&[], |(_, c)| c.as_slice())
}

/// T1K's alleles of one gene at two fields, two of them: one allele reads
/// as homozygous.
fn t1k_pair(gene: &str, alleles: &[Allele]) -> Vec<String> {
    let names: Vec<String> = alleles.iter().map(|a| two_fields(&a.name, gene)).collect();
    if names.len() == 1 {
        vec![names[0].clone(), names[0].clone()]
    } else {
        names
    }
}

/// Two allele pairs are the same pair, in either order.
fn same_pair(a: &[String], b: &[String]) -> bool {
    let (mut a, mut b) = (a.to_vec(), b.to_vec());
    a.sort();
    b.sort();
    a.len() == 2 && a == b
}

/// What step 36 decided for one gene.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Decision {
    /// Passed to PharmCAT as this diplotype.
    Passed(String),
    /// Not passed: the Result cell (`not typed` or `indeterminate`).
    Withheld(&'static str),
}

/// One row of the consensus table.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Row {
    pub gene: &'static str,
    pub decision: Decision,
    pub reason: String,
    pub evidence: String,
}

impl Row {
    fn withheld(gene: &'static str, result: &'static str, reason: &str, evidence: &str) -> Row {
        Row {
            gene,
            decision: Decision::Withheld(result),
            reason: reason.to_owned(),
            evidence: evidence.to_owned(),
        }
    }

    /// The table's cells.
    pub fn cells(&self) -> [&str; 5] {
        let (result, passed) = match &self.decision {
            Decision::Passed(d) => (d.as_str(), "yes"),
            Decision::Withheld(r) => (*r, "no"),
        };
        [self.gene, result, passed, &self.reason, &self.evidence]
    }
}

/// The decision for one HLA gene. `t1k`: T1K's calls, or `None` when step
/// 08 did not run; `optitype`: OptiType's, or `None` without a result;
/// `need_optitype`: the second caller was asked for.
pub fn hla_row(
    gene: &'static str,
    t1k: Option<&T1k>,
    optitype: Option<&OptiType>,
    need_optitype: bool,
) -> Row {
    let Some(t1k) = t1k else {
        return Row::withheld(gene, "not typed", "HLA typing (step 08) did not run", "-");
    };
    let alleles = alleles_of(t1k, gene);
    if alleles.is_empty() {
        return Row::withheld(gene, "not typed", "T1K reported no allele", "-");
    }
    let mut evidence = alleles
        .iter()
        .map(|a| {
            format!(
                "{} (quality {})",
                a.name,
                pycompat::fmt::general(a.quality, 6)
            )
        })
        .collect::<Vec<_>>()
        .join("; ");
    if let Some(ot) = optitype {
        let calls = optitype_of(ot, gene);
        let said = if calls.is_empty() {
            "no call".to_owned()
        } else {
            calls.join("/")
        };
        evidence = format!("T1K: {evidence}; OptiType: {said}");
    }
    if alleles.iter().any(Allele::unusable) {
        let why = "T1K gives an allele quality 0 or below (T1K says to ignore it)";
        return Row::withheld(gene, "indeterminate", why, &evidence);
    }
    let names = t1k_pair(gene, alleles);
    let dip = names.join("/");
    match optitype {
        None if need_optitype => Row::withheld(gene, "indeterminate", NO_OPTITYPE, &evidence),
        Some(ot) => {
            let calls = optitype_of(ot, gene);
            if calls.is_empty() {
                let why = "OptiType, the second caller, made no call";
                Row::withheld(gene, "indeterminate", why, &evidence)
            } else if !same_pair(&names, calls) {
                let why = "T1K and OptiType disagree at two fields";
                Row::withheld(gene, "indeterminate", why, &evidence)
            } else {
                Row {
                    gene,
                    decision: Decision::Passed(dip),
                    reason: "T1K and OptiType (step 08), which agree at two fields".to_owned(),
                    evidence,
                }
            }
        }
        None => {
            let mut reason = "T1K (step 08)".to_owned();
            if alleles.len() == 1 {
                reason.push_str(": one allele, read as homozygous");
            }
            Row {
                gene,
                decision: Decision::Passed(dip),
                reason,
                evidence,
            }
        }
    }
}

/// Step 08's concordance table: per class I gene the two-field calls of
/// both callers and whether they agree (`yes`, `no` or `not compared`).
pub fn concordance_rows(t1k: &T1k, optitype: &OptiType) -> Vec<[String; 5]> {
    let mut rows = Vec::new();
    for gene in OPTITYPE_GENES {
        let alleles = alleles_of(t1k, gene);
        let ot = optitype_of(optitype, gene);
        let t_txt = if alleles.is_empty() {
            "no call".to_owned()
        } else {
            t1k_pair(gene, alleles).join("/")
        };
        let o_txt = if ot.is_empty() {
            "no call".to_owned()
        } else {
            ot.join("/")
        };
        let (agree, mut note) = if alleles.is_empty() || ot.is_empty() {
            ("not compared", "a caller made no call".to_owned())
        } else if alleles.iter().any(Allele::unusable) {
            (
                "not compared",
                "T1K gives an allele quality 0 or below".to_owned(),
            )
        } else if same_pair(&t1k_pair(gene, alleles), ot) {
            ("yes", "-".to_owned())
        } else {
            ("no", "the callers disagree at two fields".to_owned())
        };
        if HLA_GENES.contains(&gene) && agree != "yes" {
            note.push_str("; step 36 does not pass it to PharmCAT");
        }
        rows.push([gene.to_owned(), t_txt, o_txt, agree.to_owned(), note]);
    }
    rows
}

// --- CYP2D6 -----------------------------------------------------------------

/// A diplotype that is a call: not a no-call text, two alleles, one call.
pub fn is_call(diplotype: &str) -> bool {
    let d = strip(diplotype);
    !no_call(d) && d.contains('/') && !d.contains(';')
}

/// Equal up to the order of the two alleles.
pub fn same_diplotype(a: &str, b: &str) -> bool {
    fn key(d: &str) -> Vec<&str> {
        let mut k: Vec<&str> = d.split('/').map(strip).collect();
        k.sort_unstable();
        k
    }
    key(a) == key(b)
}

/// The rows of a tab-separated table with a header, as `csv.DictReader`
/// reads them; `stop` ends the reading early with its value.
fn dict_rows<T>(
    path: &str,
    mut stop: impl FnMut(&pycompat::csv::DictRow) -> Option<T>,
) -> Result<Option<T>, Error> {
    let lines = TextLines::open(path, Utf8::Strict)?;
    for row in DictReader::new(lines, '\t', Quoting::Minimal) {
        let row = row.map_err(|e| csv_error(e, path))?;
        if let Some(v) = stop(&row) {
            return Ok(Some(v));
        }
    }
    Ok(None)
}

/// `(row.get(key) or "").strip()`.
fn cell(row: &pycompat::csv::DictRow, key: &str) -> String {
    strip(row.get(key).flatten().unwrap_or("")).to_owned()
}

/// pypgx's CYP2D6 diplotype (`""` when its summary has no CYP2D6 row).
fn read_pypgx(path: &str) -> Result<String, Error> {
    let got = dict_rows(path, |row| {
        (cell(row, "Gene") == "CYP2D6").then(|| cell(row, "Diplotype"))
    })?;
    Ok(got.unwrap_or_default())
}

/// `(genotype, filter)` of the first sample row of Cyrius's TSV.
fn read_cyrius(path: &str) -> Result<(String, String), Error> {
    let got = dict_rows(path, |row| {
        Some((cell(row, "Genotype"), cell(row, "Filter")))
    })?;
    Ok(got.unwrap_or_default())
}

/// `(status, message)` of the CYP2D6 depth check's `key<TAB>value` table.
fn read_depth_check(path: &str) -> Result<(String, String), Error> {
    let (mut status, mut message) = (String::new(), String::new());
    for line in TextLines::open(path, Utf8::Strict)? {
        let line = line.map_err(|e| Error::io(&e, path))?;
        let (k, v) = chomp(&line).split_once('\t').unwrap_or((chomp(&line), ""));
        match k {
            "status" => v.clone_into(&mut status),
            "message" => v.clone_into(&mut message),
            _ => {}
        }
    }
    Ok((status, message))
}

/// The decision for CYP2D6 from the files given (`None`: the step did not
/// run, or the option was empty).
pub fn cyp2d6_row(
    pypgx: Option<&str>,
    cyrius: Option<&str>,
    depth: Option<&str>,
) -> Result<Row, Error> {
    let pg = pypgx.map(read_pypgx).transpose()?;
    let cy = cyrius.map(read_cyrius).transpose()?;
    let dc = depth.map(read_depth_check).transpose()?;

    let pg_ok = pg.as_deref().is_some_and(is_call);
    let cy_ok = cy.as_ref().is_some_and(|(g, f)| is_call(g) && f == "PASS");
    let pg_txt = match &pg {
        None => "not run".to_owned(),
        Some(p) if pg_ok => p.clone(),
        Some(p) if p.is_empty() => "no call (no CYP2D6 row)".to_owned(),
        Some(p) => format!("no call ({p})"),
    };
    let cy_txt = match &cy {
        None => "not run".to_owned(),
        Some((g, _)) if cy_ok => g.clone(),
        Some((g, f)) => {
            let said: Vec<&str> = [g.as_str(), f.as_str()]
                .into_iter()
                .filter(|x| !x.is_empty())
                .collect();
            if said.is_empty() {
                "no call (no row)".to_owned()
            } else {
                format!("no call ({})", said.join(", "))
            }
        }
    };
    let dc_txt = match &dc {
        None => "not run",
        Some((s, _)) if s.is_empty() => "unreadable",
        Some((s, _)) => s.as_str(),
    };
    let evidence = format!("pypgx: {pg_txt}; Cyrius: {cy_txt}; depth check: {dc_txt}");
    let withheld = |reason: &str| Ok(Row::withheld("CYP2D6", "indeterminate", reason, &evidence));

    if let Some((status, message)) = &dc
        && status != "ok"
    {
        return withheld(if message.is_empty() {
            "the CYP2D6 depth check failed"
        } else {
            message
        });
    }
    if !pg_ok && !cy_ok {
        return withheld("no caller made a call");
    }
    if !(pg_ok && cy_ok) {
        return withheld(ONE_CALLER);
    }
    let (pg, cy) = (pg.unwrap_or_default(), cy.unwrap_or_default());
    if !same_diplotype(&pg, &cy.0) {
        return withheld("pypgx and Cyrius disagree");
    }
    if dc.is_none() {
        return withheld("the CYP2D6 depth check did not run");
    }
    Ok(Row {
        gene: "CYP2D6",
        decision: Decision::Passed(pg),
        reason: "pypgx and Cyrius agree; the depth check passed".to_owned(),
        evidence,
    })
}

// --- the outside-call file and the consensus table ---------------------------

/// What `og pgx outside-calls` reads.
#[derive(Debug, Default)]
pub struct Inputs {
    pub hla: Option<String>,
    pub optitype: Option<String>,
    pub require_optitype: bool,
    pub pypgx: Option<String>,
    pub cyrius: Option<String>,
    pub depth_check: Option<String>,
}

/// The consensus rows of the three genes, in table order.
pub fn decide(inputs: &Inputs) -> Result<Vec<Row>, Error> {
    fn given(o: &Option<String>) -> Option<&str> {
        o.as_deref().filter(|p| !p.is_empty())
    }
    let t1k = given(&inputs.hla)
        .map(|p| read_t1k(p, &HLA_GENES))
        .transpose()?;
    let optitype = given(&inputs.optitype).map(read_optitype).transpose()?;
    let mut rows: Vec<Row> = HLA_GENES
        .iter()
        .map(|g| hla_row(g, t1k.as_ref(), optitype.as_ref(), inputs.require_optitype))
        .collect();
    rows.push(cyp2d6_row(
        given(&inputs.pypgx),
        given(&inputs.cyrius),
        given(&inputs.depth_check),
    )?);
    Ok(rows)
}

/// A csv table: the header and the rows, `\t`-separated, `\n`-ended.
fn csv_table<'a>(header: &[&str], rows: impl IntoIterator<Item = [&'a str; 5]>) -> Vec<u8> {
    let mut w = pycompat::csv::Writer::new(Vec::new(), '\t');
    // Writing to a Vec cannot fail.
    let _ = w.write_row(header);
    for r in rows {
        let _ = w.write_row(r);
    }
    w.into_inner()
}

/// Writes both files whole, or neither: each goes to a temporary file
/// first, and when the second rename fails the first file is removed.
pub fn write_both(first: (&str, &[u8]), second: (&str, &[u8])) -> Result<(), Error> {
    let stage = |(path, body): (&str, &[u8])| {
        let fail = |e: io::Error| Error::io(&e, path);
        let mut f = og_io::AtomicFile::create(path).map_err(fail)?;
        f.write_all(body).map_err(fail)?;
        Ok::<_, Error>(f)
    };
    let a = stage(first)?;
    let b = stage(second)?;
    a.commit().map_err(|e| Error::io(&e, first.0))?;
    b.commit().map_err(|e| {
        let _ = std::fs::remove_file(first.0);
        Error::io(&e, second.0)
    })
}

/// `og pgx outside-calls`: writes both files and returns the log.
pub fn outside_calls(
    inputs: &Inputs,
    calls_path: &str,
    consensus_path: &str,
) -> Result<String, Error> {
    let rows = decide(inputs)?;
    let mut calls = String::new();
    for r in &rows {
        if let Decision::Passed(dip) = &r.decision {
            if dip.is_empty() || dip.contains(['\t', '\n']) {
                return Err(Error(format!(
                    "refusing to write a malformed call ({}, {})",
                    og_filters::py_repr(r.gene),
                    og_filters::py_repr(dip)
                )));
            }
            let _ = writeln!(calls, "{}\t{dip}", r.gene);
        }
    }
    let table = csv_table(&CONSENSUS_HEADER, rows.iter().map(Row::cells));
    write_both((calls_path, calls.as_bytes()), (consensus_path, &table))?;
    let mut log = String::new();
    for r in &rows {
        let [gene, result, passed, reason, _] = r.cells();
        let said = if passed == "yes" {
            "passed to PharmCAT"
        } else {
            "not passed"
        };
        let _ = writeln!(log, "{gene}: {result} ({said}; {reason})");
    }
    Ok(log)
}

/// `og pgx hla-concordance`: writes the table and returns the log.
pub fn hla_concordance(t1k: &str, optitype: &str, out: &str) -> Result<String, Error> {
    let rows = concordance_rows(&read_t1k(t1k, &OPTITYPE_GENES)?, &read_optitype(optitype)?);
    let table = csv_table(
        &CONCORDANCE_HEADER,
        rows.iter().map(|r| r.each_ref().map(String::as_str)),
    );
    let fail = |e: io::Error| Error::io(&e, out);
    let mut f = og_io::AtomicFile::create(out).map_err(fail)?;
    f.write_all(&table).map_err(fail)?;
    f.commit().map_err(fail)?;
    let mut log = String::new();
    for [gene, t, o, agree, note] in &rows {
        let note = if note == "-" {
            String::new()
        } else {
            format!(" ({note})")
        };
        let _ = writeln!(log, "{gene}: T1K {t}, OptiType {o}: agree {agree}{note}");
    }
    Ok(log)
}

// --- the command line -------------------------------------------------------

const OUTSIDE_USAGE: &str =
    "usage: og pgx outside-calls --calls CALLS --consensus CONSENSUS [--hla HLA]
                            [--optitype OPTITYPE] [--require-optitype] [--pypgx PYPGX]
                            [--cyrius CYRIUS] [--depth-check DEPTH_CHECK]
";
const CONCORDANCE_USAGE: &str =
    "usage: og pgx hla-concordance --t1k T1K --optitype OPTITYPE --out OUT\n";

/// A parsed command line: `--name VALUE` options (the last one wins) and
/// flags.
pub(crate) struct Opts {
    values: Vec<(&'static str, String)>,
    flags: Vec<&'static str>,
}

impl Opts {
    /// The last value of `--name`.
    pub(crate) fn get(&self, name: &str) -> Option<&str> {
        self.values
            .iter()
            .rev()
            .find(|(n, _)| *n == name)
            .map(|(_, v)| v.as_str())
    }

    /// Every value of `--name`, in order (argparse's `append`).
    pub(crate) fn all(&self, name: &str) -> Vec<String> {
        self.values
            .iter()
            .filter(|(n, _)| *n == name)
            .map(|(_, v)| v.clone())
            .collect()
    }

    pub(crate) fn flag(&self, name: &str) -> bool {
        self.flags.contains(&name)
    }
}

/// Parses `args` as argparse would for these options and flags; `Err` is
/// the exit code after the usage and the error are printed.
pub(crate) fn parse(
    args: &[String],
    usage: &str,
    prog: &str,
    options: &[&'static str],
    flags: &[&'static str],
    required: &[&str],
    err: &mut dyn Write,
) -> Result<Opts, u8> {
    let mut o = Opts {
        values: Vec::new(),
        flags: Vec::new(),
    };
    let refuse = |err: &mut dyn Write, msg: &str| {
        let _ = write!(err, "{usage}");
        let _ = writeln!(err, "{prog}: error: {msg}");
        2
    };
    let mut extra = Vec::new();
    let mut i = 0;
    while i < args.len() {
        let a = &args[i];
        i += 1;
        if a == "-h" || a == "--help" {
            let _ = io::stdout().write_all(usage.as_bytes());
            return Err(0);
        }
        let (name, inline) = match a.split_once('=') {
            Some((n, v)) if n.starts_with("--") => (n, Some(v.to_owned())),
            _ => (a.as_str(), None),
        };
        if let Some(&f) = flags.iter().find(|&&f| f == name) {
            if inline.is_some() {
                return Err(refuse(
                    err,
                    &format!("argument {f}: ignored explicit argument"),
                ));
            }
            o.flags.push(f);
        } else if let Some(&n) = options.iter().find(|&&n| n == name) {
            let v = match inline {
                Some(v) => v,
                None => match args.get(i) {
                    Some(v) if v == "-" || !v.starts_with('-') => {
                        i += 1;
                        v.clone()
                    }
                    _ => return Err(refuse(err, &format!("argument {n}: expected one argument"))),
                },
            };
            o.values.push((n, v));
        } else {
            extra.push(a.clone());
        }
    }
    let missing: Vec<&str> = required
        .iter()
        .copied()
        .filter(|r| o.get(r).is_none())
        .collect();
    if !missing.is_empty() {
        let msg = format!(
            "the following arguments are required: {}",
            missing.join(", ")
        );
        return Err(refuse(err, &msg));
    }
    if !extra.is_empty() {
        return Err(refuse(
            err,
            &format!("unrecognized arguments: {}", extra.join(" ")),
        ));
    }
    Ok(o)
}

/// Prints the log or the error; the exit code.
pub(crate) fn finish(r: Result<String, Error>, out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    match r {
        Ok(log) => match out.write_all(log.as_bytes()).and_then(|()| out.flush()) {
            Ok(()) => 0,
            Err(e) => {
                let _ = writeln!(err, "ERROR: cannot write stdout: {e}");
                1
            }
        },
        Err(e) => {
            let _ = writeln!(err, "ERROR: {e}");
            1
        }
    }
}

/// `og pgx outside-calls ARGS`.
pub fn outside_calls_cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let o = match parse(
        args,
        OUTSIDE_USAGE,
        "og pgx outside-calls",
        &[
            "--calls",
            "--consensus",
            "--hla",
            "--optitype",
            "--pypgx",
            "--cyrius",
            "--depth-check",
        ],
        &["--require-optitype"],
        &["--calls", "--consensus"],
        err,
    ) {
        Ok(o) => o,
        Err(code) => return code,
    };
    let s = |n: &str| o.get(n).map(str::to_owned);
    let inputs = Inputs {
        hla: s("--hla"),
        optitype: s("--optitype"),
        require_optitype: o.flag("--require-optitype"),
        pypgx: s("--pypgx"),
        cyrius: s("--cyrius"),
        depth_check: s("--depth-check"),
    };
    let r = outside_calls(
        &inputs,
        o.get("--calls").unwrap_or_default(),
        o.get("--consensus").unwrap_or_default(),
    );
    finish(r, out, err)
}

/// `og pgx hla-concordance ARGS`.
pub fn hla_concordance_cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let req = ["--t1k", "--optitype", "--out"];
    let o = match parse(
        args,
        CONCORDANCE_USAGE,
        "og pgx hla-concordance",
        &req,
        &[],
        &req,
        err,
    ) {
        Ok(o) => o,
        Err(code) => return code,
    };
    let get = |n| o.get(n).unwrap_or_default();
    finish(
        hla_concordance(get("--t1k"), get("--optitype"), get("--out")),
        out,
        err,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tmp(tag: &str) -> std::path::PathBuf {
        let p = std::env::temp_dir().join(format!("og-pgx-rules-{tag}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&p);
        std::fs::create_dir_all(&p).unwrap();
        p
    }

    fn put(dir: &std::path::Path, name: &str, text: &str) -> String {
        let p = dir.join(name);
        std::fs::write(&p, text).unwrap();
        p.display().to_string()
    }

    #[test]
    fn two_fields_cuts_and_prefixes() {
        assert_eq!(two_fields("HLA-B*57:01:01:02", "HLA-B"), "*57:01");
        assert_eq!(two_fields("57:01:01", "HLA-B"), "*57:01");
        assert_eq!(two_fields("*02", "HLA-A"), "*02");
    }

    #[test]
    fn calls_and_diplotypes() {
        assert!(is_call(" *1/*2 "));
        assert!(!is_call("Indeterminate"));
        assert!(!is_call("*1/*2;*1/*3"));
        assert!(!is_call("*1"));
        assert!(same_diplotype("*2/*1", "*1/ *2"));
        assert!(!same_diplotype("*1/*4", "*1/*2"));
    }

    const T1K: &str = "HLA-A\t2\tHLA-A*02:01:01\t10\t60\tHLA-A*24:02:01\t9\t60\n\
                       HLA-B\t2\tHLA-B*07:02:01\t10\t60\tHLA-B*44:02:01\t9\t48\n";

    fn optitype(a: &str, b: &str) -> String {
        format!(
            "\tA1\tA2\tB1\tB2\tC1\tC2\tReads\tObjective\n0\t{a}\t{b}\tB*07:02\tB*44:02\tC*05:01\tC*07:02\t100\t90.0\n"
        )
    }

    #[test]
    fn optitype_agreeing_passes_and_disagreeing_withholds() {
        let d = tmp("ot");
        let t1k = read_t1k(&put(&d, "t.tsv", T1K), &HLA_GENES).unwrap();
        let same = read_optitype(&put(&d, "o1.tsv", &optitype("A*24:02", "A*02:01"))).unwrap();
        let row = hla_row("HLA-A", Some(&t1k), Some(&same), true);
        assert_eq!(row.decision, Decision::Passed("*02:01/*24:02".to_owned()));
        // The red-first check of the plan: an OptiType disagreement must
        // never be passed through.
        let other = read_optitype(&put(&d, "o2.tsv", &optitype("A*01:01", "A*02:01"))).unwrap();
        let row = hla_row("HLA-A", Some(&t1k), Some(&other), true);
        assert_eq!(row.decision, Decision::Withheld("indeterminate"));
        assert_eq!(row.reason, "T1K and OptiType disagree at two fields");
        let row = hla_row("HLA-A", Some(&t1k), None, true);
        assert_eq!(row.reason, NO_OPTITYPE);
        let _ = std::fs::remove_dir_all(d);
    }

    #[test]
    fn a_quality_that_is_not_a_finite_number_withholds() {
        let d = tmp("nan");
        let ot = read_optitype(&put(&d, "o.tsv", &optitype("A*24:02", "A*02:01"))).unwrap();
        for q in ["nan", "inf"] {
            let text = format!("HLA-A\t2\tHLA-A*02:01:01\t10\t{q}\tHLA-A*24:02:01\t9\t60\n");
            let t1k = read_t1k(&put(&d, "t.tsv", &text), &HLA_GENES).unwrap();
            let row = hla_row("HLA-A", Some(&t1k), Some(&ot), true);
            assert_eq!(row.decision, Decision::Withheld("indeterminate"), "{q}");
            let rows = concordance_rows(&t1k, &ot);
            let a = rows.iter().find(|r| r[0] == "HLA-A").unwrap();
            assert_eq!(a[3], "not compared", "{q}");
        }
        let _ = std::fs::remove_dir_all(d);
    }

    #[test]
    fn a_failed_depth_check_withholds_cyp2d6() {
        let d = tmp("dc");
        let pg = put(&d, "p.tsv", "Gene\tDiplotype\nCYP2D6\t*1/*2\n");
        // Cyrius's header and its call on separate lines, so no line pairs a
        // star-allele call with a sample label (scripts/ci/check-personal-data.sh).
        let head = "Sample\tGenotype\tFilter\n";
        let cy = put(&d, "c.tsv", &format!("{head}S\t*2/*1\tPASS\n"));
        let ok = put(&d, "ok.tsv", "status\tok\nmessage\t\n");
        let bad = put(&d, "bad.tsv", "status\tfail\nmessage\tmulti-mapped\n");
        let row = cyp2d6_row(Some(&pg), Some(&cy), Some(&ok)).unwrap();
        assert_eq!(row.decision, Decision::Passed("*1/*2".to_owned()));
        let row = cyp2d6_row(Some(&pg), Some(&cy), Some(&bad)).unwrap();
        assert_eq!(row.decision, Decision::Withheld("indeterminate"));
        assert_eq!(row.reason, "multi-mapped");
        let _ = std::fs::remove_dir_all(d);
    }

    #[test]
    fn both_files_or_neither() {
        let d = tmp("both");
        let calls = d.join("calls.tsv").display().to_string();
        let missing = d
            .join("no-such-dir")
            .join("consensus.tsv")
            .display()
            .to_string();
        let r = outside_calls(&Inputs::default(), &calls, &missing);
        assert!(r.is_err());
        assert!(
            !std::path::Path::new(&calls).exists(),
            "the calls file stands alone"
        );
        let left: Vec<_> = std::fs::read_dir(&d).unwrap().collect();
        assert!(left.is_empty(), "a temporary file was left: {left:?}");
        let _ = std::fs::remove_dir_all(d);
    }
}
