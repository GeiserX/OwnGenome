//! Two hits in one gene across callers: a called deletion plus a small
//! variant.
//!
//! Small variants alone miss a gene hit twice when one hit is a deletion:
//! an exon deletion of one copy of MUTYH plus a pathogenic variant on the
//! other, an HBA deletion plus an HBA point variant, one SMN1 copy plus an
//! SMN1 variant, a 1q21.1 deletion plus the RBM8A low-expression allele
//! (TAR syndrome: Albers et al., Nat Genet 2012, PMID 22366785). The
//! deletions come from the engines that call them: AnnotSV's class 4 and 5
//! deletions (step 05), step 35's interval table, its HBA count and
//! Parascopy's SMN1 count.
//!
//! The two are in trans (biallelic) only when the reads show it: the
//! variant sits inside the deletion and is called homozygous (one copy is
//! left, and it carries the variant), or both calls carry read-backed phase
//! (a `|` GT and the same PS) on different copies. On the same copy (cis)
//! they are no second hit; otherwise they are possibly in trans.
//!
//! The readers here take from step 05's and step 35's tables only what the
//! rule needs. The AnnotSV and Parascopy sections themselves are read by
//! their own modules.

use super::clinvar::obj;
use super::text::{self, ReadError, cell, columns, digits, rows};
use indexmap::IndexMap;
use pycompat::csv::DictRow;
use pycompat::json::Value;
use pycompat::text::isdigit;
use std::collections::{BTreeSet, HashMap};

/// Gene spans: chrom, start, end (1-based, inclusive).
pub type Spans = IndexMap<String, (String, i128, i128)>;

/// A called deletion that reaches a gene.
#[derive(Debug, Clone, PartialEq)]
pub struct Deletion {
    pub chrom: String,
    pub start: i128,
    pub end: i128,
    pub text: String,
    /// A variant inside the span is on the one copy left. False for the
    /// HBA count, which cannot say which gene of the pair is gone.
    pub inside: bool,
    /// The call's GT and phase set (AnnotSV only).
    pub gt: Option<String>,
    pub ps: Option<String>,
}

/// A small variant of a gene, as the rule reads it.
#[derive(Debug, Clone, PartialEq)]
pub struct Variant {
    pub gene: String,
    pub key: VKey,
    pub chrom: String,
    pub pos: i128,
    pub text: String,
    pub genotype: String,
    pub copies: usize,
    pub ploidy: usize,
    pub gt: String,
    pub ps: Option<String>,
    /// The severity name of a hypomorph allele.
    pub allele: Option<String>,
}

/// A variant key: a small variant `(chrom, pos, ref, alt)`.
pub type VKey = (String, i128, String, String);

/// How a deletion and a small variant in one gene pair up.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Verdict {
    /// Both copies are hit.
    Biallelic,
    /// Phase unknown, or the calls disagree: possibly both copies.
    PossiblyBiallelic,
    /// Both on the same copy: no second hit.
    Cis,
}

impl Verdict {
    pub fn word(self) -> &'static str {
        match self {
            Verdict::Biallelic => "biallelic",
            Verdict::PossiblyBiallelic => "possibly biallelic",
            Verdict::Cis => "cis",
        }
    }
}

/// The best pair of a gene.
#[derive(Debug, Clone, PartialEq)]
pub struct Entry {
    pub gene: String,
    pub verdict: Verdict,
    pub why: String,
    pub deletion: String,
    pub variant: String,
    pub genotype: String,
    pub allele: Option<String>,
}

impl Entry {
    pub fn to_json(&self) -> IndexMap<String, Value> {
        let Value::Object(m) = obj([
            ("gene", Value::Str(self.gene.clone())),
            ("verdict", Value::Str(self.verdict.word().into())),
            ("why", Value::Str(self.why.clone())),
            ("deletion", Value::Str(self.deletion.clone())),
            ("variant", Value::Str(self.variant.clone())),
            ("genotype", Value::Str(self.genotype.clone())),
            (
                "allele",
                self.allele.clone().map_or(Value::Null, Value::Str),
            ),
        ]) else {
            unreachable!()
        };
        m
    }
}

/// The copies a phased GT puts the alt allele on, or `None` when the GT
/// is not phased or a copy is not called (the missing copy could carry the
/// other allele, so no cis or trans can be read from it).
pub fn phase_haps(gt: Option<&str>) -> Option<BTreeSet<usize>> {
    let gt = gt.filter(|g| !g.is_empty() && g.contains('|'))?;
    let alleles: Vec<&str> = gt.split('|').collect();
    if alleles.contains(&".") {
        return None;
    }
    Some(
        alleles
            .iter()
            .enumerate()
            .filter(|(_, a)| **a != "0")
            .map(|(i, _)| i)
            .collect(),
    )
}

/// The verdict of one deletion and one small variant in one gene, and why.
pub fn pair_phase(dl: &Deletion, v: &Variant) -> (Verdict, String) {
    let inside = dl.inside && dl.chrom == v.chrom && dl.start <= v.pos && v.pos <= dl.end;
    // Read-backed phase first: on the same copy (cis) is no second hit,
    // wherever the variant sits.
    let hd = phase_haps(dl.gt.as_deref()).filter(|h| !h.is_empty());
    let hv = phase_haps(Some(&v.gt)).filter(|h| !h.is_empty());
    let mut phased: Option<(Verdict, String)> = None;
    if let (Some(hd), Some(hv)) = (&hd, &hv)
        && let Some(ps) = dl.ps.as_deref().filter(|p| !p.is_empty())
        && Some(ps) == v.ps.as_deref()
    {
        if hd == hv {
            return (
                Verdict::Cis,
                format!("in cis: read-backed phase puts both on the same copy (phase set {ps})"),
            );
        }
        phased = Some(if hd.is_disjoint(hv) {
            (
                Verdict::Biallelic,
                format!(
                    "in trans: read-backed phase puts them on different copies (phase set {ps})"
                ),
            )
        } else {
            (
                Verdict::PossiblyBiallelic,
                format!(
                    "read-backed phase puts the variant on both copies (phase set {ps}); check the reads"
                ),
            )
        });
    }
    if inside && v.copies >= 2 {
        return (
            Verdict::Biallelic,
            "the variant is called homozygous inside the deletion: it is hemizygous, on the one \
             copy the deletion leaves"
                .into(),
        );
    }
    if inside && v.ploidy == 2 {
        return (
            Verdict::PossiblyBiallelic,
            "a heterozygous call inside a called deletion: one of the two calls may be wrong; \
             check the reads"
                .into(),
        );
    }
    if let Some(p) = phased {
        return p;
    }
    (
        Verdict::PossiblyBiallelic,
        "phase unknown, possibly in trans (no read-backed phase links the two calls)".into(),
    )
}

/// What [`two_hits`] finds.
#[derive(Debug, Default)]
pub struct TwoHits {
    /// The best pair of each gene hit twice, in the order genes were found.
    pub entries: IndexMap<String, Entry>,
    /// Variants called homozygous inside a deletion (hemizygous): the
    /// deletion's text.
    pub hemi: HashMap<(String, VKey), String>,
    /// The best pair of each variant that has one: `<verdict> with
    /// <deletion>: <why>`.
    pub pairs: HashMap<(String, VKey), String>,
}

/// The genes hit twice.
pub fn two_hits(variants: &[Variant], dels: &HashMap<String, Vec<Deletion>>) -> TwoHits {
    let mut out = TwoHits::default();
    let mut best: HashMap<(String, VKey), Verdict> = HashMap::new();
    for v in variants {
        let k = (v.gene.clone(), v.key.clone());
        for dl in dels.get(&v.gene).map_or(&[][..], Vec::as_slice) {
            let (verdict, why) = pair_phase(dl, v);
            if verdict == Verdict::Biallelic && why.starts_with("the variant is called homozygous")
            {
                out.hemi.insert(k.clone(), dl.text.clone());
            }
            if verdict == Verdict::Cis {
                continue;
            }
            let e = Entry {
                gene: v.gene.clone(),
                verdict,
                why: why.clone(),
                deletion: dl.text.clone(),
                variant: v.text.clone(),
                genotype: if out.hemi.contains_key(&k) {
                    "hemizygous".into()
                } else {
                    v.genotype.clone()
                },
                allele: v.allele.clone(),
            };
            match out.entries.get(&v.gene) {
                Some(old) if old.verdict <= verdict => {}
                _ => {
                    out.entries.insert(v.gene.clone(), e);
                }
            }
            if best.get(&k).is_none_or(|b| verdict < *b) {
                best.insert(k.clone(), verdict);
                out.pairs.insert(
                    k.clone(),
                    format!("{} with {}: {why}", verdict.word(), dl.text),
                );
            }
        }
    }
    out
}

// --- the deletions --------------------------------------------------------------------

/// One AnnotSV row as the two-hit and ACMG rules read it.
#[derive(Debug, Clone, PartialEq)]
pub struct SvRow {
    /// `SV_chrom:SV_start-SV_end` as AnnotSV writes it.
    pub sv: String,
    pub sv_type: Option<String>,
    pub gt: Option<String>,
    pub ps: Option<String>,
    pub acmg_class: Option<String>,
    pub genes: Vec<String>,
    pub caller: Option<String>,
    /// Callers supporting it in step 22's consensus; `None` without one.
    pub support: Option<i128>,
    pub common: bool,
}

impl SvRow {
    /// `str(acmg_class)`.
    pub fn class(&self) -> &str {
        self.acmg_class.as_deref().unwrap_or("None")
    }
}

/// Step 05's checked class 4 and 5 table: its class 4 and 5 rows, `None`
/// when there is no table, or the error text the rules print.
pub fn read_annotsv(
    at: Option<(String, String)>,
    consensus: Option<String>,
) -> Result<Option<Vec<SvRow>>, String> {
    let Some((p, _)) = at else {
        return Ok(None);
    };
    let size = std::fs::metadata(&p)
        .map_err(|e| crate::os_error(&e, &p))?
        .len();
    let rows = if size > 0 {
        rows(&p).map_err(|e| e.text())?
    } else {
        Vec::new()
    };
    const NEED: [&str; 9] = [
        "ACMG_class",
        "Caller",
        "Common",
        "Exon_edge",
        "Population_AF",
        "SV_chrom",
        "SV_end",
        "SV_start",
        "SV_type",
    ];
    if let Some(r) = rows.first()
        && NEED.iter().any(|k| r.get(k).is_none())
    {
        let list: Vec<String> = NEED.iter().map(|k| crate::repr(k)).collect();
        return Err(format!(
            "{} lacks a column of [{}]",
            crate::basename(&p),
            list.join(", ")
        ));
    }
    let support = consensus_support(&rows, consensus.as_deref());
    let py = |v: Option<&str>| v.unwrap_or("None").to_owned();
    let out = rows
        .iter()
        .zip(support)
        .map(|(r, support)| {
            let fmt = sample_format(r);
            let some = |k: &str| {
                fmt.get(k)
                    .filter(|v| !v.is_empty())
                    .map(|v| (*v).to_owned())
            };
            SvRow {
                sv: format!(
                    "{}:{}-{}",
                    py(cell(r, "SV_chrom")),
                    py(cell(r, "SV_start")),
                    py(cell(r, "SV_end"))
                ),
                sv_type: cell(r, "SV_type").map(str::to_owned),
                gt: some("GT"),
                ps: some("PS"),
                acmg_class: cell(r, "ACMG_class").map(str::to_owned),
                genes: text::cell_or(r, "Gene_name")
                    .split(';')
                    .filter(|g| !g.is_empty())
                    .map(str::to_owned)
                    .collect(),
                caller: cell(r, "Caller").map(str::to_owned),
                support,
                common: cell(r, "Common") == Some("yes"),
            }
        })
        .filter(|x| matches!(x.class(), "4" | "5"))
        .collect();
    Ok(Some(out))
}

/// The first sample's FORMAT values of an AnnotSV row: the column after
/// FORMAT in the header, none when that is `Annotation_mode` or absent.
fn sample_format(r: &DictRow) -> HashMap<&str, &str> {
    let keys: Vec<&str> = r.iter().map(|(k, _)| k).collect();
    let Some(i) = keys.iter().position(|k| *k == "FORMAT") else {
        return HashMap::new();
    };
    match keys.get(i + 1) {
        Some(&k) if k != "Annotation_mode" => {
            text::format_map(text::cell_or(r, "FORMAT"), text::cell_or(r, k))
        }
        _ => HashMap::new(),
    }
}

/// How far a consensus record may sit from an AnnotSV call, each end.
const CONSENSUS_BP: i128 = 1000;

fn chrom_key(c: &str) -> &str {
    if c.as_bytes()
        .get(..3)
        .is_some_and(|prefix| prefix.eq_ignore_ascii_case(b"chr"))
    {
        // An ASCII prefix ends at a UTF-8 boundary.
        &c[3..]
    } else {
        c
    }
}

fn info_value<'a>(info: &'a str, key: &str) -> &'a str {
    info.split(';')
        .find_map(|kv| {
            let (k, v) = kv.split_once('=').unwrap_or((kv, ""));
            (k == key).then_some(v)
        })
        .unwrap_or("")
}

fn dup_or(t: &str) -> String {
    let t = t.to_uppercase();
    if t.starts_with("DUP") {
        "DUP".into()
    } else {
        t
    }
}

/// The support of each row from step 22's consensus: the nearest record of
/// the same contig and type within 1 kb at both ends; 1 without one; `None`
/// for every row when there is no consensus or it cannot be read.
fn consensus_support(rows: &[DictRow], path: Option<&str>) -> Vec<Option<i128>> {
    let none = || vec![None; rows.len()];
    let Some(path) = path.filter(|p| crate::is_file(p)) else {
        return none();
    };
    let Ok(lines) = text::lines(path) else {
        return none();
    };
    // (contig, type) to the records' (pos, end, supp).
    type Records = HashMap<(String, String), Vec<(i128, i128, i128)>>;
    let mut by: Records = HashMap::new();
    for l in lines.iter().filter(|l| !l.starts_with('#')) {
        let c: Vec<&str> = l.split('\t').collect();
        if c.len() < 8 {
            continue;
        }
        let num = |s: &str| crate::int(s).ok();
        let Some(pos) = num(c[1]) else { continue };
        let end = match info_value(c[7], "END") {
            "" => Some(pos),
            e => num(e),
        };
        let supp = match info_value(c[7], "SUPP") {
            "" => Some(0),
            e => num(e),
        };
        let (Some(end), Some(supp)) = (end, supp) else {
            continue;
        };
        by.entry((
            chrom_key(c[0]).to_owned(),
            dup_or(info_value(c[7], "SVTYPE")),
        ))
        .or_default()
        .push((pos, end, supp));
    }
    rows.iter()
        .map(|r| {
            let (Some(Ok(start)), Some(Ok(end))) = (
                cell(r, "SV_start").map(crate::int),
                cell(r, "SV_end").map(crate::int),
            ) else {
                return Some(1);
            };
            let key = (
                chrom_key(cell(r, "SV_chrom").unwrap_or("")).to_owned(),
                dup_or(cell(r, "SV_type").unwrap_or("")),
            );
            let mut best: Option<(i128, i128)> = None;
            for &(p, e, supp) in by.get(&key).map_or(&[][..], Vec::as_slice) {
                let (dp, de) = ((p - start).abs(), (e - end).abs());
                if dp <= CONSENSUS_BP && de <= CONSENSUS_BP && best.is_none_or(|(d, _)| dp + de < d)
                {
                    best = Some((dp + de, supp));
                }
            }
            Some(best.map_or(1, |(_, s)| s))
        })
        .collect()
}

/// What step 35 wrote that the rule reads.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Step35 {
    /// The interval table's deletions and duplications; `None` without it.
    pub intervals: Option<Vec<IndexMap<String, String>>>,
    /// The HBA table (metric to value); `None` without it.
    pub hba: Option<IndexMap<String, String>>,
    /// Parascopy's resolved SMN1 copy number at exon 7.
    pub smn1: Option<i128>,
}

/// SMN1 exon 7, c.840, GRCh38 1-based.
pub const SMN1_EXON7: (&str, i128) = ("chr5", 70951946);
/// The tables in paralogs/ that are not Parascopy catalog tables.
const NOT_CATALOG_TABLES: [&str; 2] = ["hba", "paraphase"];
const PARASCOPY_COLUMNS: [&str; 9] = [
    "chrom",
    "start",
    "end",
    "agCN",
    "agCN_filter",
    "agCN_qual",
    "psCN",
    "psCN_qual",
    "psCN_filter",
];
const INTERVAL_COLUMNS: [&str; 9] = [
    "locus", "name", "chrom", "start", "end", "status", "copies", "expected", "change",
];
const DOSAGE_COLUMNS: [&str; 6] = ["locus", "status", "copies", "estimate", "reason", "reading"];

fn has_all(path: &str, cols: &[&str]) -> Result<bool, String> {
    let have = columns(path).map_err(|e| e.text())?;
    Ok(cols.iter().all(|c| have.iter().any(|h| h == c)))
}

fn short_row(path: &str, rows: &[DictRow], cols: &[&str]) -> Result<(), String> {
    if let Some(i) = rows
        .iter()
        .position(|r| cols.iter().any(|k| cell(r, k).is_none()))
    {
        return Err(format!(
            "{} line {} lacks a field",
            crate::basename(path),
            i + 2
        ));
    }
    Ok(())
}

/// Step 35's HBA or CYP21A2 table as metric to value; a table without a
/// status row reads as not assessed.
fn read_hba(path: &str) -> Result<IndexMap<String, String>, String> {
    let mut out = IndexMap::new();
    for l in text::lines(path).map_err(|e| e.text())? {
        if let Some((k, v)) = l.split_once('\t')
            && k != "metric"
        {
            out.insert(k.to_owned(), v.to_owned());
        }
    }
    if !matches!(
        out.get("status").map(String::as_str),
        Some("called" | "not assessed")
    ) {
        out = IndexMap::from([
            ("status".to_owned(), "not assessed".to_owned()),
            (
                "message".to_owned(),
                format!("not assessed: {} has no status row", crate::basename(path)),
            ),
        ]);
    }
    Ok(out)
}

/// The psCN of SMN1 at exon 7 when Parascopy resolved it (filter PASS,
/// quality 20 or more, a whole number), as `parascopy_cn_at` reads it.
fn smn1_cn(rows: &[DictRow]) -> Option<i128> {
    let (chrom, pos) = SMN1_EXON7;
    for r in rows {
        let (Some(c), Some(Ok(s)), Some(Ok(e))) = (
            cell(r, "chrom"),
            cell(r, "start").map(crate::int),
            cell(r, "end").map(crate::int),
        ) else {
            continue;
        };
        let mut regions = vec![(c.to_owned(), s + 1, e)];
        for h in text::cell(r, "homologous_regions")
            .filter(|v| !v.is_empty())
            .unwrap_or("*")
            .split(',')
        {
            let h = pycompat::text::strip(h, None);
            if let Some((hc, rest)) = h.split_once(':') {
                let a: String = rest.chars().take_while(char::is_ascii_digit).collect();
                let rest = &rest[a.len()..];
                if let Some(rest) = rest.strip_prefix('-') {
                    let b: String = rest.chars().take_while(char::is_ascii_digit).collect();
                    if !hc.is_empty() && !a.is_empty() && !b.is_empty() {
                        regions.push((hc.to_owned(), crate::int(&a).ok()?, crate::int(&b).ok()?));
                    }
                }
            }
        }
        let Some(idx) = regions
            .iter()
            .position(|(c, b, e)| c == chrom && *b <= pos && pos <= *e)
        else {
            continue;
        };
        let cns: Vec<&str> = text::cell_or(r, "psCN").split(',').collect();
        let quals: Vec<&str> = text::cell_or(r, "psCN_qual").split(',').collect();
        let cn = cns.get(idx).copied().unwrap_or("?");
        let q = quals
            .get(idx)
            .and_then(|q| pycompat::num::parse_float(q).ok());
        let filt = cell(r, "psCN_filter")
            .filter(|v| !v.is_empty())
            .unwrap_or(".");
        let ok = filt == "PASS" && q.is_some_and(|q| q >= 20.0) && isdigit(cn);
        return if ok { crate::int(cn).ok() } else { None };
    }
    None
}

/// What the rule takes from step 35's tables, in the order the Parascopy
/// section reads them, so the same fault stops it. `Ok(None)`: step 35
/// wrote no table.
pub fn read_step35(dir: &str, s: &str, assets: &str) -> Result<Option<Step35>, String> {
    let paralogs = crate::join(dir, "paralogs");
    let mut tables: Vec<(String, String)> = Vec::new();
    if let Ok(entries) = std::fs::read_dir(&paralogs) {
        for e in entries.flatten() {
            let Ok(name) = e.file_name().into_string() else {
                continue;
            };
            let Some(t) = name
                .strip_prefix(&format!("{s}_"))
                .and_then(|n| n.strip_suffix("_copy_number.tsv"))
            else {
                continue;
            };
            let p = crate::join(&paralogs, &name);
            if !NOT_CATALOG_TABLES.contains(&t) && crate::is_file(&p) {
                tables.push((t.to_owned(), p));
            }
        }
    }
    tables.sort_by(|a, b| pycompat::mtime::py_cmp(&a.1, &b.1));
    let side = |n: &str| crate::join(&paralogs, &format!("{s}_{n}"));
    let (hp, ip, cp, dp) = (
        side("hba_copy_number.tsv"),
        side("cn_intervals.tsv"),
        side("cyp21a2_psv.tsv"),
        side("cn_dosage.tsv"),
    );
    let smn = tables
        .iter()
        .find(|(t, _)| t == "smn")
        .map(|(_, p)| p.clone());
    let any =
        smn.is_some() || !tables.is_empty() || [&hp, &ip, &cp, &dp].iter().any(crate::is_file);
    if !any {
        return Ok(None);
    }
    let mut smn1 = None;
    if let Some(p) = &smn {
        let r = rows(p).map_err(|e| e.text())?;
        if !has_all(p, &PARASCOPY_COLUMNS)? {
            return Err(format!(
                "{} lacks the chrom, start, end, agCN, agCN_filter, agCN_qual, psCN, psCN_qual or \
                 psCN_filter column",
                crate::basename(p)
            ));
        }
        smn1 = smn1_cn(&r);
    }
    let others: Vec<&(String, String)> = tables.iter().filter(|(t, _)| t != "smn").collect();
    if !others.is_empty() {
        let cat = super::assets::asset_rows(assets, "paralog_loci.tsv")?;
        for r in &cat {
            if r.get("table").is_none() {
                return Err("'table'".into());
            }
        }
    }
    let hba = if crate::is_file(&hp) {
        Some(read_hba(&hp)?)
    } else {
        None
    };
    let intervals = if crate::is_file(&ip) {
        let r = rows(&ip).map_err(|e| e.text())?;
        if !has_all(&ip, &INTERVAL_COLUMNS)? {
            return Err(format!(
                "{} lacks the locus, name, status, copies, expected or change column",
                crate::basename(&ip)
            ));
        }
        short_row(&ip, &r, &INTERVAL_COLUMNS)?;
        Some(
            r.iter()
                .filter(|x| matches!(cell(x, "change"), Some("deletion" | "duplication")))
                .map(|x| {
                    [
                        "locus", "name", "chrom", "start", "end", "copies", "expected", "change",
                    ]
                    .iter()
                    .map(|k| ((*k).to_owned(), text::cell_or(x, k).to_owned()))
                    .collect()
                })
                .collect(),
        )
    } else {
        None
    };
    if crate::is_file(&cp) {
        read_hba(&cp)?;
    }
    if crate::is_file(&dp) {
        if !has_all(&dp, &DOSAGE_COLUMNS)? {
            return Err(format!(
                "{} lacks the {} or reading column",
                crate::basename(&dp),
                DOSAGE_COLUMNS[..5].join(", ")
            ));
        }
        let r = rows(&dp).map_err(|e| e.text())?;
        short_row(&dp, &r, &DOSAGE_COLUMNS)?;
    }
    for (_, p) in others {
        rows(p).map_err(|e: ReadError| e.text())?;
        columns(p).map_err(|e| e.text())?;
    }
    Ok(Some(Step35 {
        intervals,
        hba,
        smn1,
    }))
}

/// The deletion inputs of the sample, read once for every rule that needs
/// them.
#[derive(Debug, Clone)]
pub struct DeletionInputs {
    pub annotsv: Result<Option<Vec<SvRow>>, String>,
    pub step35: Result<Option<Step35>, String>,
}

/// The genes a called deletion reaches: the deletions of each gene, the
/// engines read, and the engines not read with why.
pub fn gene_deletions(
    inp: &DeletionInputs,
    spans: &Spans,
) -> (HashMap<String, Vec<Deletion>>, Vec<String>, Vec<String>) {
    let mut out: HashMap<String, Vec<Deletion>> = HashMap::new();
    let (mut read, mut unread) = (Vec::new(), Vec::new());
    let over = |chrom: &str, a: i128, b: i128| -> Vec<String> {
        spans
            .iter()
            .filter(|(_, (c, gs, ge))| c == chrom && *gs <= b && *ge >= a)
            .map(|(g, _)| g.clone())
            .collect()
    };
    let mut add = |g: String, x: Deletion| out.entry(g).or_default().push(x);
    match &inp.annotsv {
        Err(e) => unread.push(format!("AnnotSV (step 05): {e}")),
        Ok(None) => unread.push("AnnotSV (step 05): no table".into()),
        Ok(Some(rows)) => {
            read.push("AnnotSV class 4 and 5 deletions (step 05)".into());
            for x in rows {
                let (chrom, span) = x.sv.split_once(':').unwrap_or((&x.sv, ""));
                let (a, b) = span.split_once('-').unwrap_or((span, ""));
                if x.sv_type.as_deref() != Some("DEL") || x.common || !(isdigit(a) && isdigit(b)) {
                    continue;
                }
                // A call genotyped 0/0 is no deletion in this sample.
                if let Some(gt) = x.gt.as_deref()
                    && gt.replace('|', "/").split('/').all(|a| a == "0")
                {
                    continue;
                }
                let chrom = if chrom.starts_with("chr") {
                    chrom.to_owned()
                } else {
                    format!("chr{chrom}")
                };
                let (Ok(a_n), Ok(b_n)) = (crate::int(a), crate::int(b)) else {
                    continue;
                };
                // SV_start is the base before the deleted bases.
                let first = a_n + 1;
                let mut genes: BTreeSet<String> = over(&chrom, first, b_n).into_iter().collect();
                genes.extend(x.genes.iter().filter(|g| !spans.contains_key(*g)).cloned());
                for g in genes {
                    add(
                        g,
                        Deletion {
                            chrom: chrom.clone(),
                            start: first,
                            end: b_n,
                            inside: true,
                            text: format!(
                                "DEL {chrom}:{a}-{b} (AnnotSV class {}, {})",
                                x.class(),
                                x.caller.as_deref().unwrap_or("None")
                            ),
                            gt: x.gt.clone(),
                            ps: x.ps.clone(),
                        },
                    );
                }
            }
        }
    }
    let st = match &inp.step35 {
        Ok(Some(st)) => st,
        Ok(None) => {
            unread.push("step 35's copy-number tables: none".into());
            return (out, read, unread);
        }
        Err(e) => {
            unread.push(format!("step 35's copy-number tables: {e}"));
            return (out, read, unread);
        }
    };
    match &st.intervals {
        Some(rows) => {
            read.push("step 35's interval table (read depth)".into());
            for r in rows {
                let f = |k: &str| r.get(k).map_or("", String::as_str);
                if f("change") != "deletion" || !(isdigit(f("start")) && isdigit(f("end"))) {
                    continue;
                }
                let (Ok(a), Ok(b)) = (crate::int(f("start")), crate::int(f("end"))) else {
                    continue;
                };
                for g in over(f("chrom"), a + 1, b) {
                    add(
                        g,
                        Deletion {
                            chrom: f("chrom").to_owned(),
                            start: a + 1,
                            end: b,
                            inside: true,
                            gt: None,
                            ps: None,
                            text: format!(
                                "{} deletion, {} of {} copies (read depth, step 35)",
                                f("locus"),
                                f("copies"),
                                f("expected")
                            ),
                        },
                    );
                }
            }
        }
        None => unread.push("step 35's interval table: none".into()),
    }
    let hba = st.hba.clone().unwrap_or_default();
    let n = hba.get("alpha_globin_genes").map_or("", String::as_str);
    if hba.get("status").map(String::as_str) == Some("called") && isdigit(n) {
        read.push("step 35's HBA1/HBA2 count".into());
        // 0 genes: no alpha-globin gene is left to carry a variant.
        if digits(Some(n)).is_some_and(|k| (1..4).contains(&k)) {
            for g in ["HBA1", "HBA2"] {
                if let Some((c, a, b)) = spans.get(g) {
                    add(
                        g.to_owned(),
                        Deletion {
                            chrom: c.clone(),
                            start: *a,
                            end: *b,
                            inside: false,
                            gt: None,
                            ps: None,
                            text: format!("{n} of 4 alpha-globin genes (read depth, step 35)"),
                        },
                    );
                }
            }
        }
    } else {
        unread.push(format!(
            "step 35's HBA1/HBA2 count: {}",
            hba.get("message")
                .filter(|m| !m.is_empty())
                .map_or("not called", String::as_str)
        ));
    }
    match st.smn1 {
        Some(cn) => {
            read.push("Parascopy's SMN1 count (step 35)".into());
            if cn == 1
                && let Some((c, a, b)) = spans.get("SMN1")
            {
                add(
                    "SMN1".into(),
                    Deletion {
                        chrom: c.clone(),
                        start: *a,
                        end: *b,
                        inside: true,
                        gt: None,
                        ps: None,
                        text: "one SMN1 copy (Parascopy, step 35)".into(),
                    },
                );
            }
        }
        None => unread.push("Parascopy's SMN1 count: not resolved".into()),
    }
    (out, read, unread)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn chromosome_prefix_is_safe_for_non_ascii_contigs() {
        for (input, expected) in [
            ("éé", "éé"),
            ("é", "é"),
            ("☃1", "☃1"),
            ("chr☃", "☃"),
            ("CHR1", "1"),
            ("cHrX", "X"),
            ("chr", ""),
            ("ch", "ch"),
            ("", ""),
        ] {
            assert_eq!(chrom_key(input), expected);
        }
    }

    fn del(gt: Option<&str>, ps: Option<&str>) -> Deletion {
        Deletion {
            chrom: "chr1".into(),
            start: 100,
            end: 200,
            text: "DEL".into(),
            inside: true,
            gt: gt.map(str::to_owned),
            ps: ps.map(str::to_owned),
        }
    }

    fn var(pos: i128, gt: &str, ps: Option<&str>) -> Variant {
        let (copies, ploidy) = super::super::clinvar::allele_copies(gt);
        Variant {
            gene: "G".into(),
            key: ("chr1".into(), pos, "A".into(), "T".into()),
            chrom: "chr1".into(),
            pos,
            text: format!("chr1:{pos} A>T"),
            genotype: "het".into(),
            copies,
            ploidy,
            gt: gt.into(),
            ps: ps.map(str::to_owned),
            allele: None,
        }
    }

    #[test]
    fn phase_reads_only_called_phased_copies() {
        assert_eq!(phase_haps(Some("0|1")), Some(BTreeSet::from([1])));
        assert_eq!(phase_haps(Some("1|0")), Some(BTreeSet::from([0])));
        assert_eq!(phase_haps(Some("0/1")), None);
        assert_eq!(phase_haps(Some(".|1")), None);
        assert_eq!(phase_haps(None), None);
    }

    #[test]
    fn read_backed_phase_decides_outside_the_deletion() {
        // Outside the deletion's span: phase alone decides.
        let trans = pair_phase(&del(Some("0|1"), Some("7")), &var(500, "1|0", Some("7")));
        assert_eq!(trans.0, Verdict::Biallelic);
        assert!(trans.1.starts_with("in trans"));
        let cis = pair_phase(&del(Some("0|1"), Some("7")), &var(500, "0|1", Some("7")));
        assert_eq!(cis.0, Verdict::Cis);
        // Another phase set links nothing.
        let other = pair_phase(&del(Some("0|1"), Some("7")), &var(500, "1|0", Some("8")));
        assert_eq!(other.0, Verdict::PossiblyBiallelic);
        assert!(other.1.starts_with("phase unknown"));
    }

    #[test]
    fn a_phased_cis_pair_is_no_second_hit_even_inside_the_deletion() {
        let v = var(150, "0|1", Some("7"));
        let cis = pair_phase(&del(Some("0|1"), Some("7")), &v);
        assert_eq!(cis.0, Verdict::Cis);
        let found = two_hits(
            &[v],
            &HashMap::from([("G".to_owned(), vec![del(Some("0|1"), Some("7"))])]),
        );
        assert!(found.entries.is_empty());
    }

    #[test]
    fn homozygous_inside_is_hemizygous() {
        let v = var(150, "1/1", None);
        let found = two_hits(
            std::slice::from_ref(&v),
            &HashMap::from([("G".to_owned(), vec![del(None, None)])]),
        );
        let e = &found.entries["G"];
        assert_eq!(e.verdict, Verdict::Biallelic);
        assert_eq!(e.genotype, "hemizygous");
        assert!(found.hemi.contains_key(&("G".to_owned(), v.key.clone())));
    }

    #[test]
    fn the_best_pair_of_a_gene_wins() {
        let v = var(500, "1|0", Some("7"));
        let found = two_hits(
            std::slice::from_ref(&v),
            &HashMap::from([(
                "G".to_owned(),
                vec![del(None, None), del(Some("0|1"), Some("7"))],
            )]),
        );
        assert_eq!(found.entries["G"].verdict, Verdict::Biallelic);
        assert!(found.pairs[&("G".to_owned(), v.key)].starts_with("biallelic with DEL: in trans"));
    }
}
