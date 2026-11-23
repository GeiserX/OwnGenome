//! `og sv postfilter`: the checks a class 4 or 5 structural variant needs
//! before anyone reads it (port of `bin/annotsv_postfilter.py`).
//!
//! AnnotSV classes every SV of its input, and some calls reach class 4 or
//! 5 without being a real rare change: a call one caller made, a
//! "deletion" whose breakpoints sit on the exon edges of one transcript (a
//! processed transcript or cDNA aligned across an intron), an SV many
//! people carry, a call the read depth does not support, and calls in
//! regions where short-read SVs mislead. Step 05 and the Nextflow ANNOTSV
//! write `<sample>_sv_class4_5_checked.tsv` from AnnotSV's class 4 and 5
//! table with these columns after AnnotSV's own ([`KEEP`]):
//!
//! - `Caller`: INFO SVCALLER, which `og sv add-callers` writes;
//! - `Exon_edge`, `Exon_edge_detail`: for a deletion, `yes` when its start
//!   lies within 5 bp of an exon's end and its end within 5 bp of a later
//!   exon's start of one transcript of `--genes`;
//! - `Population_AF`, `Population_source`, `Common`: AnnotSV's benign-set
//!   frequency for the type and its source; common at 1% or more. AnnotSV
//!   writes a frequency only when a benign record of 1% or more holds the
//!   whole call, so a blank is no evidence that the call is rare:
//!   `Population_source` then says `no frequency data` ([`NO_FREQUENCY`]);
//! - `Benign_match`, after `Common`, only when AnnotSV's table has a
//!   `B_*_coord` column: whether a benign record AnnotSV lists for the type
//!   is about the call's size or much larger ([`benign_match`]);
//! - `AnnotSV_class`: AnnotSV's own class; `ACMG_class` is the class after
//!   these checks;
//! - `Depth_support`: duphold's DHFFC (< 0.7 for a deletion) or DHBFC
//!   (> 1.3 for a duplication), else for a Delly call over 1,000 bp its
//!   RDCN by the same rule `og sv add-callers` applies
//!   ([`depth_contradicts`]); `og sv add-callers` keeps RDCN on each Delly
//!   call it adds. A Delly call over 1,000 bp with no RDCN value (a VCF not
//!   written by Delly 2.x) reads `unknown`. `no` moves the call to class 3;
//! - `Rule`, `Rule_detail`: the region rules of `--rules`
//!   (`assets/sv_rule_regions.tsv`) and the depth rule that fired.
//!
//! With `--consensus` (step 22's VCF), `Support` and `Support_callers`
//! follow, both empty when the consensus cannot tell ([`Membership`]);
//! then AnnotSV's FORMAT and sample columns. A summary line goes to
//! stderr.
//!
//! Fixes against the Python: inputs are read as gzip by their magic bytes;
//! a missing `--class45` file is an `ERROR:` line and exit 1, not a
//! traceback; an RDCN that is not a whole number of copies is no evidence
//! (`not checked`), where the Python read it as `no`; a rule row whose start
//! or end is not a number is an `ERROR:` line naming it; a failed write is
//! a message and leaves no partial table.

use std::collections::{HashMap, HashSet};
use std::fmt;
use std::io::{self, BufRead, Write};
use std::rc::Rc;

use og_filters::args::{self, ArgError, Spec};
use og_io::{Lines, Utf8};
use og_pedigree::Par;
use pycompat::csv::{DictReader, DictRow, Quoting, Writer};
use pycompat::fmt::{general, percent};

use super::{DEPTH_MIN_SIZE, Rdcn, Sex, SvKind, chrom_key, depth_contradicts, info_value};

const USAGE: &str = "usage: og sv postfilter --class45 CLASS45.tsv [--genes GENES.bed]
       [--rules RULES.tsv] [--consensus CONSENSUS.vcf.gz]
       [--sex male|female] [--par-bed PAR.bed] --out CHECKED.tsv
";

const SPEC: Spec = Spec {
    prog: "og sv postfilter",
    usage: USAGE,
    options: &[
        "--class45",
        "--genes",
        "--rules",
        "--consensus",
        "--sex",
        "--par-bed",
        "--out",
    ],
    positionals: &[],
    required: &["--class45", "--out"],
};

/// AnnotSV's columns the table keeps.
pub const KEEP: [&str; 10] = [
    "AnnotSV_ID",
    "ID",
    "SV_chrom",
    "SV_start",
    "SV_end",
    "SV_length",
    "SV_type",
    "Gene_name",
    "ACMG_class",
    "AnnotSV_ranking_score",
];
/// The columns the checks add.
pub const ADDED: [&str; 10] = [
    "Caller",
    "Exon_edge",
    "Exon_edge_detail",
    "Population_AF",
    "Population_source",
    "Common",
    "AnnotSV_class",
    "Depth_support",
    "Rule",
    "Rule_detail",
];
/// The columns `--consensus` adds.
pub const SUPPORT: [&str; 2] = ["Support", "Support_callers"];
/// The column the benign records' coordinates add, after `Common`.
pub const BENIGN_MATCH: &str = "Benign_match";
/// `Population_source` of a call AnnotSV gives no benign-set frequency.
pub const NO_FREQUENCY: &str = "no frequency data";

const EDGE_BP: i64 = 5;
const COMMON_AF: f64 = 0.01;
/// A benign record matches the call's size when each covers this share of
/// the other, the rule of the couple screen's common deletions.
const SAME_OVERLAP_PERCENT: i64 = 50;
const BENIGN_TYPES: [&str; 4] = ["B_loss", "B_gain", "B_ins", "B_inv"];
const DHFFC_MAX: f64 = 0.7;
const DHBFC_MIN: f64 = 1.3;
const REPEAT_SLACK_BP: i64 = 10;
const LOCUS_FRACTION: f64 = 0.5;
const REQUIRED: [&str; 5] = ["SV_chrom", "SV_start", "SV_end", "SV_type", "ACMG_class"];

/// A region rule of the rule table.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum RuleKind {
    TandemRepeat,
    IgTrLocus,
    Mhc,
    GeneAlias,
}

impl RuleKind {
    fn parse(s: &str) -> Option<RuleKind> {
        Some(match s {
            "tandem_repeat" => RuleKind::TandemRepeat,
            "ig_tr_locus" => RuleKind::IgTrLocus,
            "mhc" => RuleKind::Mhc,
            "gene_alias" => RuleKind::GeneAlias,
            _ => return None,
        })
    }

    fn name(self) -> &'static str {
        match self {
            RuleKind::TandemRepeat => "tandem_repeat",
            RuleKind::IgTrLocus => "ig_tr_locus",
            RuleKind::Mhc => "mhc",
            RuleKind::GeneAlias => "gene_alias",
        }
    }
}

/// One rule region: 0-based start, end, rule, name.
#[derive(Debug, Clone)]
struct Region {
    start: i64,
    end: i64,
    rule: RuleKind,
    name: String,
}

/// One transcript of the genes file.
#[derive(Debug, Clone)]
struct Transcript {
    start: i64,
    end: i64,
    tx: String,
    gene: String,
    minus: bool,
    /// (start0, end1), sorted.
    exons: Vec<(i64, i64)>,
}

/// The transcripts of the contigs asked for, in file order, with an index
/// by start for the window a deletion can reach.
#[derive(Debug, Default)]
struct Transcripts {
    by_chrom: HashMap<String, Chrom>,
}

#[derive(Debug, Default)]
struct Chrom {
    /// In file order: the first transcript that fits is the one named.
    list: Vec<Transcript>,
    /// (start, index in `list`), sorted.
    starts: Vec<(i64, usize)>,
    max_len: i64,
}

impl Transcripts {
    fn chrom(&self, c: &str) -> &[Transcript] {
        self.by_chrom.get(c).map_or(&[], |c| &c.list)
    }

    /// The indexes of the transcripts of `c` with `te >= lo` and
    /// `ts <= hi`, in file order.
    fn window(&self, c: &str, lo: i64, hi: i64) -> Vec<usize> {
        let Some(ch) = self.by_chrom.get(c) else {
            return Vec::new();
        };
        let from = ch.starts.partition_point(|&(s, _)| s < lo - ch.max_len);
        let to = ch.starts.partition_point(|&(s, _)| s <= hi);
        let mut ids: Vec<usize> = ch.starts[from..to.max(from)]
            .iter()
            .map(|&(_, i)| i)
            .filter(|&i| ch.list[i].end >= lo)
            .collect();
        ids.sort_unstable();
        ids
    }
}

/// What the command line asks for.
#[derive(Debug, Clone)]
pub struct Opts {
    pub class45: String,
    pub genes: Option<String>,
    pub rules: Option<String>,
    pub consensus: Option<String>,
    pub sex: Sex,
    pub par_bed: Option<String>,
    pub out: String,
}

/// Why the step stopped: one message, exit 1.
#[derive(Debug)]
pub struct Error(String);

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

fn read_err(path: &str, e: &io::Error) -> Error {
    Error(format!(
        "cannot read {path}: {}",
        og_filters::os_error(e, path)
    ))
}

/// The lines of a text file as Python's `open(path, errors="replace")`
/// reads them (gzip by its magic), each with its `\n`.
fn lines(path: &str) -> Result<Lines<Box<dyn BufRead>>, Error> {
    let r = og_io::open(path).map_err(|e| read_err(path, &e))?;
    Ok(Lines::new(r, Utf8::Replace))
}

fn for_each_line(
    path: &str,
    mut each: impl FnMut(u64, &str) -> Result<(), Error>,
) -> Result<(), Error> {
    let mut l = lines(path)?;
    let mut line = String::new();
    let mut n = 0;
    while l.read_line(&mut line).map_err(|e| read_err(path, &e))? {
        n += 1;
        each(n, &line)?;
    }
    Ok(())
}

fn int(s: &str) -> Option<i64> {
    pycompat::num::parse_int(s)
        .ok()
        .and_then(|v| i64::try_from(v).ok())
}

/// `float(v)` that is not NaN.
fn number(v: Option<&str>) -> Option<f64> {
    let x = pycompat::num::parse_float(v?).ok()?;
    (!x.is_nan()).then_some(x)
}

fn event_kind(kind: &str) -> String {
    let kind = pycompat::text::upper(kind);
    if kind.starts_with("DUP") {
        "DUP".to_owned()
    } else {
        kind
    }
}

/// One class 4 or 5 row: AnnotSV's cells by name.
struct Row {
    cells: DictRow,
}

impl Row {
    /// `row.get(key) or ""`: a missing column or a short row read as empty.
    fn get(&self, key: &str) -> &str {
        self.cells.get(key).flatten().unwrap_or("")
    }

    fn kind(&self) -> String {
        event_kind(self.get("SV_type"))
    }

    fn sv_kind(&self) -> SvKind {
        match self.kind().as_str() {
            "DEL" => SvKind::Del,
            "DUP" => SvKind::Dup,
            _ => SvKind::Other,
        }
    }

    fn start_end(&self) -> Option<(i64, i64)> {
        Some((int(self.get("SV_start"))?, int(self.get("SV_end"))?))
    }
}

/// The FORMAT and sample columns of AnnotSV's header: the column after
/// FORMAT, unless there is none or it is `Annotation_mode`.
fn sample_columns(keys: &[&str]) -> Option<(String, String)> {
    let i = keys.iter().position(|&k| k == "FORMAT")? + 1;
    let s = keys.get(i)?;
    (*s != "Annotation_mode").then(|| ("FORMAT".to_owned(), (*s).to_owned()))
}

/// `{FORMAT key: value}` of the sample column (the last key of a name wins).
fn sample_format<'a>(row: &'a Row, sample: Option<&(String, String)>) -> HashMap<&'a str, &'a str> {
    let Some((_, s)) = sample else {
        return HashMap::new();
    };
    row.get("FORMAT")
        .split(':')
        .zip(row.get(s).split(':'))
        .collect()
}

fn read_transcripts(path: &str, chroms: &HashSet<String>) -> Result<Transcripts, Error> {
    let mut t = Transcripts::default();
    for_each_line(path, |_, line| {
        let c: Vec<&str> = line.trim_end_matches('\n').split('\t').collect();
        if c.len() < 10 || line.starts_with('#') || !chroms.contains(chrom_key(c[0])) {
            return Ok(());
        }
        let list = |s: &str| -> Option<Vec<i64>> {
            s.split(',').filter(|x| !x.is_empty()).map(int).collect()
        };
        let (Some(start), Some(end), Some(starts), Some(ends)) =
            (int(c[1]), int(c[2]), list(c[8]), list(c[9]))
        else {
            return Ok(());
        };
        let mut exons: Vec<(i64, i64)> = starts.into_iter().zip(ends).collect();
        exons.sort_unstable();
        let ch = t.by_chrom.entry(chrom_key(c[0]).to_owned()).or_default();
        ch.starts.push((start, ch.list.len()));
        ch.max_len = ch.max_len.max(end - start);
        ch.list.push(Transcript {
            start,
            end,
            tx: c[5].to_owned(),
            gene: c[4].to_owned(),
            minus: c[3] == "-",
            exons,
        });
        Ok(())
    })?;
    for ch in t.by_chrom.values_mut() {
        ch.starts.sort_unstable();
    }
    Ok(t)
}

/// `(yes|no, detail)` for a deletion from `start` to `end` (AnnotSV's
/// SV_start, the base before the deleted bases, and SV_end, the last
/// deleted base): a deletion that removes an intron exactly starts on the
/// last base of an exon and ends on the base before the next exon.
fn exon_edge(chrom: &str, start: i64, end: i64, tx: &Transcripts) -> (&'static str, String) {
    let c = chrom_key(chrom);
    let all = tx.chrom(c);
    for i in tx.window(c, start - EDGE_BP, end + EDGE_BP) {
        let t = &all[i];
        let n = t.exons.len();
        for (a, &(_, e_end)) in t.exons.iter().enumerate() {
            if (start - e_end).abs() > EDGE_BP {
                continue;
            }
            for b in a + 1..n {
                if (end - t.exons[b].0).abs() <= EDGE_BP {
                    let (x, y) = if t.minus {
                        (n - b, n - a)
                    } else {
                        (a + 1, b + 1)
                    };
                    return (
                        "yes",
                        format!("{} ({}): exon {x} and exon {y} edges", t.tx, t.gene),
                    );
                }
            }
        }
    }
    ("no", String::new())
}

/// The prefix of AnnotSV's benign-set columns for the call's type.
fn benign_columns(row: &Row) -> Option<&'static str> {
    Some(match row.kind().as_str() {
        "DEL" => "B_loss",
        "DUP" => "B_gain",
        "INS" => "B_ins",
        "INV" => "B_inv",
        _ => return None,
    })
}

/// AnnotSV's benign-set frequency for the type, and its source.
fn population(row: &Row) -> (Option<f64>, String) {
    let Some(col) = benign_columns(row) else {
        return (None, String::new());
    };
    let raw = pycompat::num::py_strip(row.get(&format!("{col}_AFmax")));
    match pycompat::num::parse_float(raw) {
        Ok(af) => (
            Some(af),
            pycompat::num::py_strip(row.get(&format!("{col}_source"))).to_owned(),
        ),
        Err(_) => (None, String::new()),
    }
}

/// How the benign records AnnotSV lists for the call's type
/// (`B_loss_coord` for a deletion, and so on) relate to the call. AnnotSV
/// lists only records that hold the whole call (bedtools `-f 1.0`), and
/// `B_*_AFmax` is the highest frequency among them, so a small call inside
/// a large common SV takes the large one's frequency.
///
/// - `same`: a listed record and the call each cover at least half of the
///   other ([`SAME_OVERLAP_PERCENT`]);
/// - `contained`: every listed record is more than twice as long as the
///   call, so its frequency is the larger SV's, not this call's;
/// - `none`: no record is listed;
/// - `.`: a type with no benign set;
/// - `not checked`: the table has no coordinates column for the type, or a
///   position is not a number.
fn benign_match(row: &Row) -> &'static str {
    let Some(col) = benign_columns(row) else {
        return ".";
    };
    let Some(cell) = row.cells.get(&format!("{col}_coord")) else {
        return "not checked";
    };
    let coords = pycompat::num::py_strip(cell.unwrap_or(""));
    if coords.is_empty() {
        return "none";
    }
    let Some((a, b)) = row.start_end() else {
        return "not checked";
    };
    // An insertion can span no base: give it one, so it can overlap.
    let (s, e) = (a.min(b), a.max(b).max(a.min(b) + 1));
    let chrom = chrom_key(row.get("SV_chrom"));
    let mut read = false;
    for c in coords.split(';') {
        let Some((bc, span)) = c.trim().rsplit_once(':') else {
            continue;
        };
        let Some((Some(bs), Some(be))) = span.split_once('-').map(|(x, y)| (int(x), int(y))) else {
            continue;
        };
        if chrom_key(bc) != chrom {
            continue;
        }
        read = true;
        let ov = e.min(be) - s.max(bs);
        if ov > 0
            && ov * 100 >= SAME_OVERLAP_PERCENT * (e - s)
            && ov * 100 >= SAME_OVERLAP_PERCENT * (be - bs)
        {
            return "same";
        }
    }
    if read { "contained" } else { "not checked" }
}

/// `(yes|no|unknown|not checked|., detail)` of the read depth of a deletion or
/// duplication.
fn depth_support(
    row: &Row,
    caller: &str,
    fmt: &HashMap<&str, &str>,
    sex: Sex,
    par: &Par,
) -> (&'static str, String) {
    let kind = row.sv_kind();
    if kind == SvKind::Other {
        return (".", String::new());
    }
    if caller == "cnvpytor" {
        return ("yes", "a read-depth call (CNVpytor)".to_owned());
    }
    let (key, ok): (&str, fn(f64) -> bool) = if kind == SvKind::Del {
        ("DHFFC", |v| v < DHFFC_MAX)
    } else {
        ("DHBFC", |v| v > DHBFC_MIN)
    };
    if let Some(v) = number(fmt.get(key).copied()) {
        return (
            if ok(v) { "yes" } else { "no" },
            format!("duphold {key} {}", general(v, 6)),
        );
    }
    let rdcn = number(fmt.get("RDCN").copied());
    let len = {
        let declared = match row.get("SV_length") {
            "" => Some(0),
            s => int(s).map(i64::abs),
        };
        match declared {
            Some(0) => row.start_end().map_or(0, |(s, e)| (e - s).abs()),
            Some(n) => n,
            None => 0,
        }
    };
    if let Some(v) = rdcn
        && len > DEPTH_MIN_SIZE
        && let Rdcn::Copies(n) = Rdcn::from_number(v)
    {
        // The rule of og sv add-callers, so a call reads the same in both
        // places.
        let chrom = format!(
            "chr{}",
            pycompat::text::upper(chrom_key(row.get("SV_chrom")))
        );
        let start = int(row.get("SV_start")).unwrap_or(0);
        let no = depth_contradicts(
            kind,
            &chrom,
            len,
            Rdcn::Copies(n),
            sex,
            par.contains(&chrom, start),
        );
        return (
            if no { "no" } else { "yes" },
            format!("Delly RDCN {}", general(v, 6)),
        );
    }
    if let Some(v) = rdcn
        && len > DEPTH_MIN_SIZE
    {
        // A fractional, infinite or negative RDCN is no count of copies.
        return (
            "not checked",
            format!(
                "Delly RDCN {} is not a whole number of copies",
                general(v, 6)
            ),
        );
    }
    if caller == "delly" && rdcn.is_none() && len > DEPTH_MIN_SIZE {
        return ("unknown", "no Delly RDCN in the record".to_owned());
    }
    (
        "not checked",
        "no duphold or Delly depth value in the record".to_owned(),
    )
}

fn read_rules(path: &str, chroms: &HashSet<String>) -> Result<HashMap<String, Vec<Region>>, Error> {
    let mut out: HashMap<String, Vec<Region>> = HashMap::new();
    for_each_line(path, |n, line| {
        if line.starts_with('#') || pycompat::num::py_strip(line).is_empty() {
            return Ok(());
        }
        let c: Vec<&str> = line.trim_end_matches('\n').split('\t').collect();
        let rule = c.get(3).and_then(|r| RuleKind::parse(r));
        let (Some(rule), true) = (rule, c.len() >= 6) else {
            return Err(Error(format!(
                "{path} line {n}: want chrom, start, end, rule (gene_alias, ig_tr_locus, mhc, tandem_repeat), action, name, source"
            )));
        };
        if chroms.contains(chrom_key(c[0])) {
            let (Some(start), Some(end)) = (int(c[1]), int(c[2])) else {
                return Err(Error(format!(
                    "{path} line {n}: start and end must be whole numbers"
                )));
            };
            out.entry(chrom_key(c[0]).to_owned())
                .or_default()
                .push(Region {
                    start,
                    end,
                    rule,
                    name: c[5].to_owned(),
                });
        }
        Ok(())
    })?;
    Ok(out)
}

/// The rule hits of one row, its new Gene_name, and whether the class
/// goes to 3.
fn region_rules(
    row: &Row,
    rules: &HashMap<String, Vec<Region>>,
    tx: Option<&Transcripts>,
) -> (Vec<(&'static str, String)>, Option<String>, bool) {
    let Some((start, end)) = row.start_end() else {
        return (Vec::new(), None, false);
    };
    let end = end.max(start + 1);
    let chrom = chrom_key(row.get("SV_chrom"));
    let genes: Vec<&str> = row
        .get("Gene_name")
        .split(';')
        .filter(|g| !g.is_empty())
        .collect();
    let mut hits: Vec<(&'static str, String)> = Vec::new();
    let mut seen: Vec<RuleKind> = Vec::new();
    let mut remap: Vec<(String, String)> = Vec::new();
    for r in rules.get(chrom).map_or(&[][..], Vec::as_slice) {
        if r.start >= end + REPEAT_SLACK_BP || r.end <= start - REPEAT_SLACK_BP {
            continue;
        }
        let at = format!("{}:{}-{}", row.get("SV_chrom"), r.start + 1, r.end);
        let fresh = !seen.contains(&r.rule);
        match r.rule {
            RuleKind::TandemRepeat
                if fresh
                    && start >= r.start - REPEAT_SLACK_BP
                    && end <= r.end + REPEAT_SLACK_BP =>
            {
                seen.push(r.rule);
                hits.push((r.rule.name(), format!(
                    "inside the coding tandem repeat {} at {at} (TRExplorer): a change in the repeat's length, which SV callers do not size; read it with a repeat genotyper",
                    r.name
                )));
            }
            RuleKind::IgTrLocus | RuleKind::Mhc if fresh => {
                let frac = (end.min(r.end) - start.max(r.start)) as f64 / (end - start) as f64;
                if frac > LOCUS_FRACTION {
                    seen.push(r.rule);
                    let pct = percent(frac, 0);
                    hits.push((r.rule.name(), if r.rule == RuleKind::Mhc {
                        format!("{pct} inside the MHC ({at}), whose haplotypes differ widely between people; gnomAD-SV marks such SVs low confidence (IGH_MHC_OVERLAP)")
                    } else {
                        format!("{pct} inside the {} locus ({at}), where blood cells cut DNA out when they rearrange their receptor genes (V(D)J recombination)", r.name)
                    }));
                }
            }
            RuleKind::GeneAlias if end.min(r.end) > start.max(r.start) => {
                let Some(tx) = tx else { continue };
                let mut parts = r.name.split('|').chain(std::iter::repeat(""));
                let (alias, pseudo, pseudo_id, gene_id) = (
                    parts.next().unwrap_or(""),
                    parts.next().unwrap_or(""),
                    parts.next().unwrap_or(""),
                    parts.next().unwrap_or(""),
                );
                if !genes.contains(&alias) || remap.iter().any(|(a, _)| a == alias) {
                    continue;
                }
                let real = tx.chrom(chrom).iter().any(|t| {
                    t.gene == alias
                        && t.end > start
                        && t.start < end
                        && !(t.start >= r.start && t.end <= r.end)
                });
                if !real {
                    remap.push((alias.to_owned(), pseudo.to_owned()));
                    hits.push((r.rule.name(), format!(
                        "{alias} here is the pseudogene {pseudo} ({pseudo_id}), which carries {alias} as an alias; the gene {alias} ({gene_id}) lies elsewhere"
                    )));
                }
            }
            _ => {}
        }
    }
    let mut downgrade = hits.iter().any(|(rule, _)| *rule != "gene_alias");
    let mut new_genes = None;
    if !remap.is_empty() {
        let mapped = |g: &str| -> String {
            remap
                .iter()
                .find(|(a, _)| a == g)
                .map_or_else(|| g.to_owned(), |(_, p)| p.clone())
        };
        let mut out: Vec<String> = Vec::new();
        for g in &genes {
            let m = mapped(g);
            if !out.contains(&m) {
                out.push(m);
            }
        }
        new_genes = Some(out.join(";"));
        if genes.iter().all(|g| remap.iter().any(|(a, _)| a == g)) {
            downgrade = true;
        } else {
            for (rule, d) in &mut hits {
                if *rule == "gene_alias" {
                    d.push_str("; the class may rest on the other genes, so it stays");
                }
            }
        }
    }
    (hits, new_genes, downgrade)
}

/// SURVIVOR's own membership, indexed by contig, type, caller and source
/// record ID. Several caller IDs can point to the same set of votes.
#[derive(Default)]
pub struct ConsensusSupport {
    names: Vec<String>,
    members: HashMap<MemberKey, Option<Rc<[String]>>>,
    /// Callers that voted in an event whose member ID cannot be read: a
    /// record of theirs that no key names may still be in an event.
    unnamed: HashSet<String>,
}

/// Whether a checked row's source record is a voting member of a
/// consensus event.
#[derive(Debug, PartialEq, Eq)]
pub enum Membership<'a> {
    /// A voting member: the event's voting callers, in column order.
    Member(&'a [String]),
    /// In no event, and the consensus names every member of this caller:
    /// the call is the row's caller's alone.
    NotMember,
    /// Not known: the row has no ID (empty or `.`), its caller is not a
    /// column of the consensus, its ID is a member of two events, or a vote
    /// of its caller carries no ID this can read.
    Unknown,
}

#[derive(Eq, Hash, PartialEq)]
struct MemberKey {
    chrom: String,
    kind: String,
    caller: String,
    id: String,
}

impl MemberKey {
    fn new(chrom: &str, kind: &str, caller: &str, id: &str) -> Self {
        Self {
            chrom: chrom_key(chrom).to_owned(),
            kind: event_kind(kind),
            caller: caller.to_owned(),
            // SURVIVOR 1.0.7 replaces colons in FORMAT/ID with underscores.
            id: id.replace(':', "_"),
        }
    }
}

impl ConsensusSupport {
    /// Read a header or record from the step 22 VCF. A record whose FILTER
    /// is not `PASS` or `.` is no event. A PASS record whose votes cannot
    /// be read (no SVTYPE, a SUPP_VEC or SUPP that does not fit the
    /// columns, no FORMAT/ID) marks each caller it may count as unnamed.
    pub fn add_line(&mut self, line: &str) {
        if line.starts_with("#CHROM\t") {
            self.names = line.split('\t').skip(9).map(str::to_owned).collect();
            if self.names.iter().collect::<HashSet<_>>().len() != self.names.len() {
                self.names.clear();
            }
            return;
        }
        if line.starts_with('#') {
            return;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        let [chrom, _, _, _, _, _, "PASS" | ".", info, rest @ ..] = fields.as_slice() else {
            return;
        };
        let bits = info_value(info, "SUPP_VEC").as_bytes();
        let bits_ok =
            bits.len() == self.names.len() && bits.iter().all(|&b| matches!(b, b'0' | b'1'));
        let voters = bits.iter().filter(|&&b| b == b'1').count();
        let kind = info_value(info, "SVTYPE");
        let id_index = rest
            .first()
            .and_then(|format| format.split(':').position(|key| key == "ID"));
        let samples = rest.get(1..).unwrap_or_default();
        let readable = bits_ok
            && !kind.is_empty()
            && samples.len() == self.names.len()
            && voters > 0
            && info_value(info, "SUPP").parse() == Ok(voters);
        let (true, Some(id_index)) = (readable, id_index) else {
            for (i, name) in self.names.iter().enumerate() {
                if !bits_ok || bits[i] == b'1' {
                    self.unnamed.insert(name.clone());
                }
            }
            return;
        };
        let callers: Rc<[String]> = self
            .names
            .iter()
            .zip(bits)
            .filter(|&(_, &bit)| bit == b'1')
            .map(|(name, _)| name.clone())
            .collect();
        for ((caller, sample), &bit) in self.names.iter().zip(samples).zip(bits) {
            if bit != b'1' {
                continue;
            }
            match sample.split(':').nth(id_index) {
                Some(id) if !matches!(id, "" | ".") => {
                    self.members
                        .entry(MemberKey::new(chrom, kind, caller, id))
                        // A source ID in two events has no unique one.
                        .and_modify(|old| *old = None)
                        .or_insert_with(|| Some(Rc::clone(&callers)));
                }
                _ => {
                    self.unnamed.insert(caller.clone());
                }
            }
        }
    }

    /// Whether the source record `id` of `caller` is a voting member of an
    /// event of this contig and type, and that event's voting callers.
    pub fn membership(&self, chrom: &str, kind: &str, caller: &str, id: &str) -> Membership<'_> {
        if matches!(id, "" | ".") || !self.names.iter().any(|name| name == caller) {
            return Membership::Unknown;
        }
        match self.members.get(&MemberKey::new(chrom, kind, caller, id)) {
            Some(Some(callers)) => Membership::Member(callers),
            Some(None) => Membership::Unknown,
            None if self.unnamed.contains(caller) => Membership::Unknown,
            None => Membership::NotMember,
        }
    }
}

/// The checked table: header and rows of cells.
struct Checked {
    cols: Vec<String>,
    rows: Vec<HashMap<String, String>>,
}

/// Runs the step; returns the summary line.
pub fn run(o: &Opts) -> Result<String, Error> {
    let par = match o.par_bed.as_deref() {
        Some(p) if !std::path::Path::new(p).is_file() => {
            return Err(Error(format!("no PAR file {p}")));
        }
        p => Par::read(p).map_err(|d| Error(d.lines.join("\n")))?,
    };
    let meta = std::fs::metadata(&o.class45).map_err(|e| read_err(&o.class45, &e))?;
    let mut keys: Vec<String> = Vec::new();
    let mut rows: Vec<Row> = Vec::new();
    if meta.len() > 0 {
        // The csv reader wants each line with its ending, as Python's file
        // gives it, or a quoted cell that spans lines loses its newline.
        let mut l = lines(&o.class45)?;
        let kept = std::iter::from_fn(move || {
            let mut s = String::new();
            l.read_line(&mut s)
                .map(|more| more.then_some(s))
                .transpose()
        });
        let mut reader = DictReader::new(kept, '\t', Quoting::Minimal);
        let bad = |e: pycompat::csv::Error| Error(format!("{}: {e}", o.class45));
        for r in reader.by_ref() {
            rows.push(Row {
                cells: r.map_err(bad)?,
            });
        }
        if let Some(first) = rows.first() {
            keys = first.cells.iter().map(|(k, _)| k.to_owned()).collect();
        }
    }
    if !rows.is_empty() && !REQUIRED.iter().all(|c| keys.iter().any(|k| k == c)) {
        return Err(Error(format!(
            "{} lacks one of SV_chrom, SV_start, SV_end, SV_type, ACMG_class",
            o.class45
        )));
    }
    let chroms: HashSet<String> = rows
        .iter()
        .map(|r| chrom_key(r.get("SV_chrom")).to_owned())
        .collect();
    let tx = match o.genes.as_deref() {
        Some(p) if !std::path::Path::new(p).is_file() => {
            return Err(Error(format!("no genes file {p}")));
        }
        Some(p) => Some(read_transcripts(p, &chroms)?),
        None => None,
    };
    let rules = match o.rules.as_deref() {
        Some(p) if !std::path::Path::new(p).is_file() => {
            return Err(Error(format!("no rule table {p}")));
        }
        Some(p) => Some(read_rules(p, &chroms)?),
        None => None,
    };
    let key_refs: Vec<&str> = keys.iter().map(String::as_str).collect();
    let sample = sample_columns(&key_refs);
    let no_rules = HashMap::new();
    // Benign_match needs the records' coordinates, which AnnotSV 3 writes;
    // a table without them keeps the columns it had.
    let coords = BENIGN_TYPES
        .iter()
        .any(|t| keys.iter().any(|k| *k == format!("{t}_coord")));

    let mut checked = Checked {
        cols: KEEP
            .iter()
            .chain(ADDED.iter())
            .map(|s| (*s).to_owned())
            .collect(),
        rows: Vec::with_capacity(rows.len()),
    };
    if coords && let Some(i) = checked.cols.iter().position(|c| c == "Common") {
        checked.cols.insert(i + 1, BENIGN_MATCH.to_owned());
    }
    for r in &rows {
        let mut o_row: HashMap<String, String> = KEEP
            .iter()
            .map(|k| ((*k).to_owned(), r.get(k).to_owned()))
            .collect();
        if let Some((f, s)) = &sample {
            o_row.insert(f.clone(), r.get(f).to_owned());
            o_row.insert(s.clone(), r.get(s).to_owned());
        }
        let caller = match info_value(r.get("INFO"), "SVCALLER") {
            "" => ".".to_owned(),
            c => c.to_owned(),
        };
        let (edge, edge_detail) = if r.kind() != "DEL" {
            (".", String::new())
        } else {
            match (&tx, r.start_end()) {
                (None, _) | (_, None) => ("not checked", String::new()),
                (Some(t), Some((s, e))) => exon_edge(r.get("SV_chrom"), s, e, t),
            }
        };
        let (af, src) = population(r);
        let fmt = sample_format(r, sample.as_ref());
        let (depth, depth_detail) = depth_support(r, &caller, &fmt, o.sex, &par);
        let (mut hits, new_genes, mut downgrade) =
            region_rules(r, rules.as_ref().unwrap_or(&no_rules), tx.as_ref());
        if depth == "no" {
            hits.insert(
                0,
                (
                    "no_depth_support",
                    format!("read depth does not show the change ({depth_detail})"),
                ),
            );
            downgrade = true;
        }
        if let Some(g) = new_genes {
            o_row.insert("Gene_name".to_owned(), g);
        }
        let class = r.get("ACMG_class");
        if downgrade && number(Some(class)).is_some_and(|c| c > 3.0) {
            o_row.insert("ACMG_class".to_owned(), "3".to_owned());
        }
        let set = |row: &mut HashMap<String, String>, k: &str, v: String| {
            row.insert(k.to_owned(), v);
        };
        set(&mut o_row, "Caller", caller);
        set(&mut o_row, "Exon_edge", edge.to_owned());
        set(&mut o_row, "Exon_edge_detail", edge_detail);
        set(
            &mut o_row,
            "Population_AF",
            af.map_or(String::new(), |a| general(a, 6)),
        );
        set(
            &mut o_row,
            "Population_source",
            if af.is_some() {
                src
            } else {
                NO_FREQUENCY.to_owned()
            },
        );
        if coords {
            set(&mut o_row, BENIGN_MATCH, benign_match(r).to_owned());
        }
        set(
            &mut o_row,
            "Common",
            if af.is_some_and(|a| a >= COMMON_AF) {
                "yes"
            } else {
                "no"
            }
            .to_owned(),
        );
        set(&mut o_row, "AnnotSV_class", class.to_owned());
        set(&mut o_row, "Depth_support", depth.to_owned());
        set(
            &mut o_row,
            "Rule",
            hits.iter().map(|h| h.0).collect::<Vec<_>>().join(";"),
        );
        set(
            &mut o_row,
            "Rule_detail",
            hits.iter()
                .map(|h| h.1.as_str())
                .collect::<Vec<_>>()
                .join(" | "),
        );
        checked.rows.push(o_row);
    }

    if let Some(path) = o.consensus.as_deref() {
        support(&mut checked.rows, path)?;
        checked.cols.extend(SUPPORT.iter().map(|s| (*s).to_owned()));
    }
    // Last: a sample named like a column above would lose its GT and PS.
    if let Some((f, s)) = &sample {
        if let Some(clash) = [f, s].into_iter().find(|k| checked.cols.contains(k)) {
            return Err(Error(format!(
                "the sample column {clash} has the name of a checked-table column, so its GT and PS would be lost; rename the sample"
            )));
        }
        checked.cols.push(f.clone());
        checked.cols.push(s.clone());
    }
    write(&o.out, &checked)?;

    let count = |col: &str, v: &str| checked.rows.iter().filter(|r| r[col] == v).count();
    let moved = checked
        .rows
        .iter()
        .filter(|r| r["ACMG_class"] != r["AnnotSV_class"])
        .count();
    Ok(format!(
        "Class 4-5 SVs checked: {}; deletions on exon edges: {}; common (benign-set AF >= 1%): {}; moved to class 3: {moved}{}; depth support no: {}, not checked: {}{}",
        checked.rows.len(),
        if tx.is_some() {
            count("Exon_edge", "yes").to_string()
        } else {
            "not checked".to_owned()
        },
        count("Common", "yes"),
        if rules.is_some() {
            ""
        } else {
            " (no rule table)"
        },
        count("Depth_support", "no"),
        count("Depth_support", "not checked"),
        match count("Depth_support", "unknown") {
            0 => String::new(),
            n => format!(", unknown: {n}"),
        },
    ))
}

/// Support and Support_callers of each row from its SURVIVOR membership
/// ([`ConsensusSupport::membership`]): a member has its event's voting
/// callers, a row in no event its own caller alone, and a row whose
/// membership is not known both cells empty, as without the file.
fn support(rows: &mut [HashMap<String, String>], path: &str) -> Result<(), Error> {
    if !std::path::Path::new(path).is_file() {
        for r in rows {
            r.insert("Support".to_owned(), String::new());
            r.insert("Support_callers".to_owned(), String::new());
        }
        return Ok(());
    }
    let mut consensus = ConsensusSupport::default();
    for_each_line(path, |_, line| {
        consensus.add_line(line.trim_end_matches('\n'));
        Ok(())
    })?;
    for r in rows {
        let (s, c) =
            match consensus.membership(&r["SV_chrom"], &r["SV_type"], &r["Caller"], &r["ID"]) {
                Membership::Member(callers) => (callers.len().to_string(), callers.join(",")),
                Membership::NotMember => ("1".to_owned(), r["Caller"].clone()),
                Membership::Unknown => (String::new(), String::new()),
            };
        r.insert("Support".to_owned(), s);
        r.insert("Support_callers".to_owned(), c);
    }
    Ok(())
}

/// Writes the table whole or not at all.
fn write(path: &str, t: &Checked) -> Result<(), Error> {
    let err = |e: io::Error| Error(format!("cannot write {path}: {e}"));
    let f = og_io::AtomicFile::create(path).map_err(err)?;
    let mut w = Writer::new(f, '\t');
    w.write_row(&t.cols).map_err(err)?;
    for r in &t.rows {
        w.write_row(t.cols.iter().map(|c| r.get(c).map_or("", String::as_str)))
            .map_err(err)?;
    }
    let mut f = w.into_inner();
    f.flush().map_err(err)?;
    f.commit().map_err(err)
}

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
    // An empty path is "not given", as the Python's `if a.genes:` reads it.
    let given = |name| p.get(name).filter(|v| !v.is_empty()).map(str::to_owned);
    Ok(Opts {
        class45: p.get("--class45").unwrap_or("").to_owned(),
        genes: given("--genes"),
        rules: given("--rules"),
        consensus: given("--consensus"),
        sex,
        par_bed: given("--par-bed"),
        out: p.get("--out").unwrap_or("").to_owned(),
    })
}

/// `og sv postfilter ARGS`: exit 0, 1 for an input or output problem, 2
/// for a usage error.
pub fn cli(argv: &[String]) -> u8 {
    let o = match opts(argv) {
        Ok(o) => o,
        Err(code) => return code,
    };
    match run(&o) {
        Ok(summary) => {
            eprintln!("{summary}");
            0
        }
        Err(e) => {
            eprintln!("ERROR: {e}");
            1
        }
    }
}

/// The written table's header, for tests.
#[cfg(test)]
pub(crate) fn header_of(path: &str) -> Vec<String> {
    std::fs::read_to_string(path)
        .unwrap_or_default()
        .lines()
        .next()
        .unwrap_or("")
        .split('\t')
        .map(str::to_owned)
        .collect()
}
