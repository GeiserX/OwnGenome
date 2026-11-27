//! The ACMG secondary-findings tier of the clinical section: the variants
//! that are secondary findings by the rule of their gene
//! (`assets/acmg_sf.tsv`, ACMG SF v3.3, Lee et al., Genet Med 2025).
//!
//! Three sources: step 06's ClinVar P/LP hits, the clinical filter's rare
//! HIGH-impact records (step 23, ranked MANE Select, then canonical, then
//! any transcript, as `og vcf csq-rank` ranks them) and step 05's class 4
//! and 5 SVs in the genes whose rule allows them. A variant of an SF gene
//! that its rule does not report (one MUTYH variant, a heterozygous HFE
//! p.C282Y) is listed under `not_reported` with the reason. A called
//! deletion counts as the second hit of a gene whose rule needs two
//! ([`two_hits`](super::two_hits)), unless read-backed phase puts both on
//! one copy.
//!
//! A gene whose rule needs two hits counts the copies of the gene hit, not
//! the alt copies of its records. A homozygote hits both copies. Phased
//! heterozygotes (a `|` GT with a phase set, PS) count the copies their
//! phase set covers: two in cis count once, two in trans twice. The two alt
//! alleles of one split call (`1/.` and `./1` at one position) sit on
//! different copies. Any other heterozygote counts once with its phase
//! unknown, so a pair of them is reported as possibly biallelic. A record
//! with no alt copy called counts nothing, and a combined allele of step 06
//! counts once with its parts, unless it is unphased and two of its parts
//! are candidates on their own (they may sit on different copies).
//!
//! A ClinVar hit is mapped to its gene's rule through GENEINFO, which names
//! the gene, not the disease. So each ClinVar row keeps the condition
//! ClinVar names for the allele (step 06 copies `CLNDN` and `CLNDISDB`) and
//! the disease ACMG lists the gene for (`phenotype`). When every condition
//! named is another disease (a RET variant asserted for Hirschsprung
//! disease only, not for MEN2), the row carries `condition_not_acmg: true`.
//! It is still reported: the mark is for the reader, never a filter.
//!
//! A HIGH-impact record is a predicted loss of function, so it counts only
//! in a gene where loss of function causes the disease ACMG lists it for:
//! the `lof` column, built by `og dev build acmg-lof` from ClinGen's
//! haploinsufficiency scores and the reviewed hand rows of
//! `assets/acmg_lof_overrides.tsv`. A PCSK9 stop-gain (PCSK9 causes familial
//! hypercholesterolemia by gain of function) is listed under `not_reported`
//! with the reason; TGFBR1, TGFBR2 and RBM20 count it with a
//! `mechanism_note`. A HIGH-impact record of unknown frequency, or with no
//! GQ, DP or AD to judge its call (`call_unverified`), counts but is not
//! counted as rare, and a pair that needs one carries a `pair_note`.
//!
//! The `limit` column applies two footnotes of Table 1. TTN counts only a
//! truncating variant in its exons with high PSI (SF v3.2 footnote f); the
//! regions come from cardiodb's exon table, through `og dev build
//! acmg-lof`; a record lies in one when any base its REF spans does. A
//! ClinVar hit is truncating when its record in the clinical VCF is; a
//! hit that file does not hold is reported with a `rule_condition` to
//! check it. An LMNA variant that may be a missense one
//! carries the case-evidence condition in `rule_condition`.

use super::assets::{
    ACMG_SF_VERSION, Key, Limit, Lof, Rule, acmg_sf_rules, asset_rows, carrier_genes,
};
use super::clinvar::{
    COMMON_AF, Condition, Floor, Hit, HitRecords, allele_copies, condition_fields, floor_values,
    obj, opt_int, opt_s,
};
use super::text::{self, digits, format_map, gt_field, info_map};
use super::two_hits::{
    DeletionInputs, Entry, Spans, Variant, Verdict, gene_deletions, phase_haps, two_hits_with,
};
use indexmap::IndexMap;
use pycompat::json::Value;
use pycompat::text::isdigit;
use std::collections::{BTreeSet, HashMap, HashSet};

/// The genotype check's floors: step 06's GQ and DP, the alt read share a
/// heterozygote must lie within and a homozygous alt call must reach.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ClinicalFloor {
    pub call: Floor,
    pub het_ab: (f64, f64),
    pub hom_ab: f64,
}

impl Default for ClinicalFloor {
    fn default() -> Self {
        ClinicalFloor {
            call: Floor::default(),
            het_ab: (0.2, 0.8),
            hom_ab: 0.8,
        }
    }
}

/// The length of the `[0-9.]+` run at the start of `s`.
fn run(s: &str) -> usize {
    s.bytes()
        .take_while(|c| c.is_ascii_digit() || *c == b'.')
        .count()
}

/// `re.search(KEY([0-9.]+)-([0-9.]+))` (with `pair`) or `KEY([0-9.]+)`:
/// the first place the whole pattern matches.
fn search<'a>(line: &'a str, key: &str, pair: bool) -> Option<(&'a str, Option<&'a str>)> {
    line.match_indices(key).find_map(|(i, _)| {
        let rest = &line[i + key.len()..];
        let n = run(rest);
        if n == 0 {
            return None;
        }
        if !pair {
            return Some((&rest[..n], None));
        }
        let hi = rest[n..].strip_prefix('-')?;
        let m = run(hi);
        (m > 0).then(|| (&rest[..n], Some(&hi[..m])))
    })
}

fn float(s: &str) -> Result<f64, String> {
    crate::float(s).map_err(|e| format!("ValueError: {e}"))
}

impl ClinicalFloor {
    /// Reads the floors `og vcf csq-rank` wrote into the header.
    pub fn read(lines: &[String]) -> Result<ClinicalFloor, String> {
        let mut f = ClinicalFloor::default();
        for l in lines.iter().take_while(|l| l.starts_with("##")) {
            if !l.starts_with("##PGP_CALL_FLOOR=") {
                continue;
            }
            floor_values(l, true, &mut f.call);
            if let Some((lo, Some(hi))) = search(l, "HET_AB=", true) {
                f.het_ab = (float(lo)?, float(hi)?);
            }
            if let Some((v, _)) = search(l, "HOM_AB=", false) {
                f.hom_ab = float(v)?;
            }
        }
        Ok(f)
    }
}

/// One ranked CSQ entry.
#[derive(Debug, Clone, PartialEq)]
pub struct CsqEntry {
    pub impact: String,
    pub gene: String,
    pub consequence: String,
}

/// `csq_rank`: the transcript class (MANE, CANONICAL or ANY), its entries,
/// and the frequency word (rare, common, unknown, or none when the CSQ has
/// no frequency field) from MAX_AF, else gnomADe_AF and gnomADg_AF, over
/// every entry.
pub fn csq_rank(fmt: &[String], csq: &str) -> (&'static str, Vec<CsqEntry>, &'static str) {
    let ix: HashMap<&str, usize> = fmt
        .iter()
        .enumerate()
        .map(|(i, k)| (k.as_str(), i))
        .collect();
    let entries: Vec<Vec<&str>> = csq
        .split(',')
        .filter(|t| !t.is_empty())
        .map(|t| t.split('|').collect::<Vec<_>>())
        .filter(|v| v.len() == fmt.len())
        .collect();
    let get = |v: &[&str], k: &str| ix.get(k).map_or("", |&i| v[i]).to_owned();
    let mane: Vec<&Vec<&str>> = entries
        .iter()
        .filter(|v| !get(v, "MANE_SELECT").is_empty())
        .collect();
    let canon: Vec<&Vec<&str>> = entries
        .iter()
        .filter(|v| get(v, "CANONICAL") == "YES")
        .collect();
    let (cls, chosen) = if !mane.is_empty() {
        ("MANE", mane)
    } else if !canon.is_empty() {
        ("CANONICAL", canon)
    } else {
        ("ANY", entries.iter().collect())
    };
    let fields: Vec<&str> = if ix.contains_key("MAX_AF") {
        vec!["MAX_AF"]
    } else {
        ["gnomADe_AF", "gnomADg_AF"]
            .into_iter()
            .filter(|f| ix.contains_key(f))
            .collect()
    };
    let afs: Vec<f64> = entries
        .iter()
        .flat_map(|v| {
            fields
                .iter()
                .flat_map(|f| {
                    get(v, f)
                        .split('&')
                        .filter_map(|a| pycompat::num::parse_float(a).ok())
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>()
        })
        .collect();
    let freq = if fields.is_empty() {
        "none"
    } else {
        match text::py_max(afs) {
            None => "unknown",
            Some(m) if m < COMMON_AF => "rare",
            Some(_) => "common",
        }
    };
    let out = chosen
        .iter()
        .map(|v| CsqEntry {
            impact: get(v, "IMPACT"),
            gene: get(v, "SYMBOL"),
            consequence: get(v, "Consequence").replace('&', ","),
        })
        .collect();
    (cls, out, freq)
}

/// The genotype check of one call: weak, the reasons, GQ, DP and AD.
pub fn call_check(
    fmt: &str,
    sample: &str,
    floor: &ClinicalFloor,
) -> (
    bool,
    Vec<String>,
    Option<i128>,
    Option<i128>,
    Option<String>,
) {
    let f = format_map(fmt, sample);
    let gq = digits(f.get("GQ").copied());
    let dp = digits(f.get("DP").copied());
    let ad = f.get("AD").copied().unwrap_or("");
    let mut why = Vec::new();
    if let Some(g) = gq
        && g < floor.call.gq
    {
        why.push(format!("GQ {g} < {}", floor.call.gq));
    }
    if let Some(d) = dp
        && d < floor.call.dp
    {
        why.push(format!("DP {d} < {}", floor.call.dp));
    }
    let gt = f.get("GT").copied().unwrap_or("").replace('|', "/");
    let al: Vec<&str> = gt.split('/').collect();
    let depths: Vec<&str> = ad.split(',').collect();
    if al.len() == 2
        && al.iter().all(|a| isdigit(a))
        && depths.len() > 1
        && depths.iter().all(|x| isdigit(x))
    {
        let n: Vec<i128> = depths.iter().filter_map(|x| crate::int(x).ok()).collect();
        let tot: i128 = n.iter().sum();
        if tot != 0 {
            let ab = n[1..].iter().sum::<i128>() as f64 / tot as f64;
            let (lo, hi) = floor.het_ab;
            if (al[0] == "0") != (al[1] == "0") && !(lo <= ab && ab <= hi) {
                why.push(format!(
                    "alt reads {} of the depth on a heterozygote",
                    pycompat::fmt::fixed(ab, 2)
                ));
            } else if al[0] == al[1] && al[1] != "0" && ab < floor.hom_ab {
                why.push(format!(
                    "alt reads {} of the depth on a homozygous call",
                    pycompat::fmt::fixed(ab, 2)
                ));
            }
        }
    }
    (
        !why.is_empty(),
        why,
        gq,
        dp,
        (!ad.is_empty()).then(|| ad.to_owned()),
    )
}

/// One record of the clinical VCF with a HIGH-impact ranked entry.
#[derive(Debug, Clone, PartialEq)]
pub struct HighRecord {
    pub genotype: String,
    pub ps: Option<String>,
    /// Gene to its first HIGH consequence.
    pub genes: IndexMap<String, String>,
    pub transcript: &'static str,
    pub frequency: &'static str,
    pub weak: bool,
    pub why: Vec<String>,
    pub gq: Option<i128>,
    pub dp: Option<i128>,
    pub ad: Option<String>,
}

impl HighRecord {
    /// The call has nothing to judge it by: no number in GQ or DP, and no
    /// read total in AD on a call of two numbered alleles. Step 23 writes
    /// such a call `CALL_CHECK=unknown` (`og vcf csq-rank`, the same rule);
    /// the card marks it `call_unverified` and does not count it as rare.
    /// A vendor VCF can carry such calls; DeepVariant writes GQ, DP and AD.
    pub fn unverified(&self) -> bool {
        if self.weak || self.gq.is_some() || self.dp.is_some() {
            return false;
        }
        let gt = self.genotype.replace('|', "/");
        let al: Vec<&str> = gt.split('/').collect();
        let depths: Vec<&str> = self.ad.as_deref().unwrap_or("").split(',').collect();
        let judged = al.len() == 2
            && al.iter().all(|a| isdigit(a))
            && depths.len() > 1
            && depths.iter().all(|x| isdigit(x))
            && depths
                .iter()
                .filter_map(|x| crate::int(x).ok())
                .sum::<i128>()
                != 0;
        !judged
    }
}

/// The consequences that truncate the protein, as ACMG counts a truncating
/// TTN variant: nonsense, frameshift and canonical splice.
const TRUNCATING: [&str; 4] = [
    "stop_gained",
    "frameshift_variant",
    "splice_donor_variant",
    "splice_acceptor_variant",
];

/// Whether a ranked consequence (terms joined by `,`) truncates.
fn truncating(consequence: &str) -> bool {
    consequence.split(',').any(|c| TRUNCATING.contains(&c))
}

/// `chr2` and `2` name one chromosome.
fn same_chrom(a: &str, b: &str) -> bool {
    a.strip_prefix("chr").unwrap_or(a) == b.strip_prefix("chr").unwrap_or(b)
}

/// Whether a candidate counts under its gene's mechanism (`lof`) and
/// footnote (`limit`): `Err` is the reason it is listed apart, `Ok` the
/// condition a reported row must carry, if any. `high` says the candidate
/// is a HIGH-impact record; `consequence` is its ranked HIGH consequence
/// (for a ClinVar hit, that of its record in the clinical VCF, `None` when
/// that VCF has no HIGH-impact record of it); `clinical` is `None` when no
/// clinical VCF was read, else whether it holds a record of the variant.
fn gate(
    gene: &str,
    rule: &Rule,
    key: &Key,
    high: bool,
    consequence: Option<&str>,
    clinical: Option<bool>,
) -> Result<Option<String>, String> {
    if high && rule.lof == Lof::No {
        return Err(format!(
            "a predicted loss-of-function call does not count in {gene}: loss of function is not \
             how {gene} causes the disease ACMG lists it for (lof no in assets/acmg_sf.tsv, from \
             ClinGen's haploinsufficiency score); only a ClinVar P/LP variant counts"
        ));
    }
    let trunc = consequence.is_some_and(truncating);
    match &rule.limit {
        Limit::None => Ok(None),
        Limit::TruncatingIn { regions, name } => {
            let rule = format!(
                "SF v3.2 footnote f: {gene} is reported only for a nonsense, frameshift or \
                 canonical splice variant in an exon with high PSI, and OwnGenome counts one \
                 only in {name}"
            );
            // The bases the record's REF spans, POS to POS + len(REF) - 1: a
            // deletion that starts before a region and deletes into it (its
            // canonical splice site) lies in it.
            let last = key.1 + key.2.len().max(1) as i128 - 1;
            let inside = regions
                .iter()
                .any(|(c, start, end)| same_chrom(&key.0, c) && *start <= last && key.1 <= *end);
            if !inside {
                Err(format!("{rule}: this variant lies outside them"))
            } else if trunc {
                Ok(None)
            } else if high {
                Err(format!("{rule}: this call does not truncate the protein"))
            } else {
                match clinical {
                    Some(true) => Err(format!(
                        "{rule}: its record in the clinical VCF does not truncate the protein"
                    )),
                    Some(false) => Ok(Some(format!(
                        "{rule}. The clinical VCF holds no record of this variant, so its \
                         consequence is not known here: check that it truncates the protein."
                    ))),
                    None => Ok(Some(format!(
                        "{rule}. With no clinical VCF its consequence is not known here: check \
                         that it truncates the protein."
                    ))),
                }
            }
        }
        Limit::MissenseNeedsCaseEvidence if trunc => Ok(None),
        Limit::MissenseNeedsCaseEvidence if high => Ok(None),
        Limit::MissenseNeedsCaseEvidence => Ok(Some(format!(
            "SF v3.3 footnote e: ACMG reports a missense variant of {gene} only with case \
             evidence of cardiac disease, which OwnGenome cannot see; check that before you \
             read this as a finding"
        ))),
    }
}

/// The clinical VCF read once: its records count and its HIGH-impact
/// records. An error is the note of a file that cannot be read.
#[derive(Debug, Clone)]
pub struct ClinicalVcf {
    pub records: usize,
    /// Keyed by `(chrom, pos, ref, alt)`; an error is the collector's note
    /// (a POS that is no integer on a HIGH-impact record, a bad floor).
    pub high: Result<IndexMap<Key, HighRecord>, String>,
    /// The `(chrom, pos, ref, alt)` of every record, HIGH-impact or not, so
    /// a ClinVar hit the file holds with a consequence that does not
    /// truncate is told apart from one it does not hold at all. Step 23's
    /// clinical VCF has a few hundred records.
    pub held: HashSet<Key>,
}

impl ClinicalVcf {
    /// One pass over the file, line by line: a whole-genome clinical VCF
    /// is never held in memory.
    pub fn read(path: &str) -> Result<ClinicalVcf, String> {
        let mut records = 0;
        let mut scan = HighScan::default();
        text::each_line(path, |l| {
            if !l.starts_with('#') {
                records += 1;
            }
            scan.line(l);
        })
        .map_err(|e| e.unreadable(path))?;
        let held = std::mem::take(&mut scan.held);
        Ok(ClinicalVcf {
            records,
            high: scan.finish(),
            held,
        })
    }
}

/// `csq_high_impact_genes`, fed one line at a time. The floors come from
/// the header (`clinical_call_floor`), which ends at the first line that
/// does not start with `##`.
#[derive(Default)]
struct HighScan {
    header: Vec<String>,
    floor: Option<Result<ClinicalFloor, String>>,
    fmt: Option<Vec<String>>,
    out: IndexMap<Key, HighRecord>,
    err: Option<String>,
    held: HashSet<Key>,
}

impl HighScan {
    fn line(&mut self, line: &str) {
        if self.floor.is_none() {
            if line.starts_with("##") {
                self.header.push(line.to_owned());
            } else {
                self.floor = Some(ClinicalFloor::read(&std::mem::take(&mut self.header)));
            }
        }
        if !line.starts_with('#') {
            let f: Vec<&str> = line.splitn(6, '\t').collect();
            if let (Some(c), Some(Ok(p)), Some(r), Some(a)) = (
                f.first(),
                f.get(1).map(|p| crate::int(p)),
                f.get(3),
                f.get(4),
            ) {
                self.held
                    .insert(((*c).to_owned(), p, (*r).to_owned(), (*a).to_owned()));
            }
        }
        if self.err.is_some() {
            return;
        }
        if line.starts_with("##INFO=<ID=CSQ") {
            self.fmt = line.split_once("Format: ").and_then(|(_, f)| {
                let f = f.split('"').next().unwrap_or("");
                (!f.is_empty()).then(|| f.split('|').map(str::to_owned).collect())
            });
            return;
        }
        let Some(Ok(floor)) = &self.floor else { return };
        let Some(fm) = self.fmt.as_ref() else { return };
        if line.starts_with('#')
            || !fm.iter().any(|c| c == "IMPACT")
            || !fm.iter().any(|c| c == "SYMBOL")
        {
            return;
        }
        if let Err(e) = high_record(fm, line, floor, &mut self.out) {
            self.err = Some(e);
        }
    }

    fn finish(mut self) -> Result<IndexMap<Key, HighRecord>, String> {
        let floor = self
            .floor
            .take()
            .unwrap_or_else(|| ClinicalFloor::read(&self.header));
        floor?;
        match self.err {
            Some(e) => Err(e),
            None => Ok(self.out),
        }
    }
}

/// One data line of `csq_high_impact_genes`.
fn high_record(
    fm: &[String],
    line: &str,
    floor: &ClinicalFloor,
    out: &mut IndexMap<Key, HighRecord>,
) -> Result<(), String> {
    let r: Vec<String> = line.split('\t').map(str::to_owned).collect();
    if r.len() < 8 {
        return Ok(());
    }
    let (cls, entries, freq) = csq_rank(fm, info_map(&r[7]).get("CSQ").copied().unwrap_or(""));
    let mut genes: IndexMap<String, String> = IndexMap::new();
    for e in entries.into_iter().filter(|e| e.impact == "HIGH") {
        genes.entry(e.gene).or_insert(e.consequence);
    }
    if genes.is_empty() {
        return Ok(());
    }
    let gt = text::gt_of(&r).to_owned();
    let (weak, why, gq, dp, ad) = if r.len() > 9 {
        call_check(&r[8], &r[9], floor)
    } else {
        (false, Vec::new(), None, None, None)
    };
    let pos = text::pos(&r)?;
    out.insert(
        (r[0].clone(), pos, r[3].clone(), r[4].clone()),
        HighRecord {
            genotype: gt,
            ps: gt_field(&r, "PS"),
            genes,
            transcript: cls,
            frequency: freq,
            weak,
            why,
            gq,
            dp,
            ad,
        },
    );
    Ok(())
}

/// A source of a candidate variant.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum Kind {
    ClinVar,
    High,
    Sv,
}

impl Kind {
    fn source(self) -> &'static str {
        match self {
            Kind::ClinVar => "ClinVar",
            Kind::High => "HIGH impact",
            Kind::Sv => "AnnotSV",
        }
    }
}

/// A variant key the verdicts are counted by: a small variant, an SV by
/// its span, or a gene's called deletion.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
enum AKey {
    Small(Key),
    Sv(String, String, String, String),
    Deletion(String),
}

/// One candidate: the fields it is printed with (in their order), its key,
/// copies and ploidy, and its call's GT and phase set.
#[derive(Debug, Clone)]
struct Cand {
    kind: Kind,
    gene: String,
    fields: IndexMap<String, Value>,
    key: AKey,
    copies: usize,
    ploidy: usize,
    gt: String,
    ps: Option<String>,
    /// The parts of a combined allele (step 06's `PGP_JOINED_FROM`).
    joined: Vec<Key>,
    /// The reason its gene's mechanism or footnote leaves it uncounted
    /// ([`gate`]); such a candidate is listed apart and counts nothing.
    apart: Option<String>,
}

impl Cand {
    /// A HIGH-impact call of unknown frequency, or with nothing to judge
    /// the call by: it counts, but a pair that needs it carries a note.
    fn unsure(&self) -> bool {
        self.kind == Kind::High
            && (matches!(
                self.fields.get("frequency_unknown"),
                Some(Value::Bool(true))
            ) || matches!(self.fields.get("call_unverified"), Some(Value::Bool(true))))
    }
}

/// One distinct variant of a gene as the verdict counts it: its alt
/// copies, ploidy, GT and phase set, whether it stands for a combined
/// allele and its parts, and the alt copies its confident calls give (`firm`:
/// 0 when every call of it is low confidence: below step 06's call floor, or
/// a weak HIGH-impact call), with the GT and phase set of the confident call
/// that gives them.
#[derive(Debug, Clone, PartialEq)]
struct Allele {
    copies: usize,
    ploidy: usize,
    gt: String,
    ps: Option<String>,
    combined: bool,
    firm: usize,
    firm_gt: String,
    firm_ps: Option<String>,
    /// Some call of it is more than a HIGH-impact call of unknown frequency
    /// or with no GQ, DP or AD ([`Cand::unsure`]). An SV or a deletion is.
    sure: bool,
}

impl Allele {
    fn of(copies: usize, ploidy: usize) -> Allele {
        Allele {
            copies,
            ploidy,
            gt: String::new(),
            ps: None,
            combined: false,
            firm: 0,
            firm_gt: String::new(),
            firm_ps: None,
            sure: true,
        }
    }

    /// Count a confident call of this allele: the most copies one gives,
    /// with that call's GT and phase set.
    fn confirm(&mut self, copies: usize, gt: &str, ps: &Option<String>) {
        if copies > self.firm {
            self.firm = copies;
            self.firm_gt = gt.to_owned();
            self.firm_ps = ps.clone();
        }
    }

    /// The copies read-backed phase puts the alt allele on, with its phase
    /// set: only a `|` GT with every copy called and a non-empty PS.
    fn phased(&self) -> Option<(&str, BTreeSet<usize>)> {
        let ps = self.ps.as_deref().filter(|p| !p.is_empty())?;
        let h = phase_haps(Some(&self.gt)).filter(|h| !h.is_empty())?;
        Some((ps, h))
    }

    /// The copy a half call (`1/.`, `./1`) puts its alt allele on.
    fn half(&self) -> Option<usize> {
        let a: Vec<&str> = self.gt.split(['/', '|']).collect();
        match a.as_slice() {
            [x, "."] if !matches!(*x, "0" | ".") => Some(0),
            [".", y] if !matches!(*y, "0" | ".") => Some(1),
            _ => None,
        }
    }
}

/// How a gene whose rule needs two hits met it: `biallelic` when a call
/// proves both copies hit, `possibly biallelic` when the count needs a
/// heterozygote of unknown phase; and why.
type Pair = (&'static str, String);

/// The copies of a gene hit, over its counted variants: the count, the
/// proof of two copies if a call gives one, and the phase sets in which
/// two or more variants sit on one copy (cis).
fn copies_hit(kept: &IndexMap<AKey, Allele>) -> (usize, Option<String>, Vec<String>) {
    let mut n = 0;
    let mut proof: Option<String> = None;
    let mut sets: IndexMap<&str, (BTreeSet<usize>, usize)> = IndexMap::new();
    let mut halves: IndexMap<(&str, i128), (BTreeSet<usize>, usize)> = IndexMap::new();
    for (k, a) in kept {
        let AKey::Small(key) = k else {
            n += 1;
            continue;
        };
        if a.ploidy == 1 {
            n += a.copies;
        } else if a.copies >= 2 {
            n += 2;
            proof.get_or_insert_with(|| {
                format!(
                    "homozygous: one call puts {}:{} {}>{} on both copies",
                    key.0, key.1, key.2, key.3
                )
            });
        } else if let Some((ps, h)) = a.phased() {
            let s = sets.entry(ps).or_default();
            s.0.extend(h);
            s.1 += 1;
        } else if let Some(i) = a.half() {
            let s = halves.entry((key.0.as_str(), key.1)).or_default();
            s.0.insert(i);
            s.1 += 1;
        } else {
            n += 1;
        }
    }
    let mut cis = Vec::new();
    for (ps, (h, count)) in &sets {
        n += h.len();
        if h.len() >= 2 {
            proof.get_or_insert_with(|| {
                format!(
                    "in trans: read-backed phase puts them on different copies (phase set {ps})"
                )
            });
        } else if *count >= 2 {
            cis.push((*ps).to_owned());
        }
    }
    for ((chrom, pos), (h, count)) in &halves {
        if h.len() >= 2 {
            n += 2 + count - h.len();
            proof.get_or_insert_with(|| {
                format!("in trans: the two alt alleles of one call at {chrom}:{pos}")
            });
        } else {
            n += count;
        }
    }
    (n, proof, cis)
}

/// What [`sf_verdicts`] gives a gene: for each distinct variant `None`
/// (reported) or the reason it is not, the partial-deficiency note, how a
/// gene whose rule needs two hits met it, and the note of a pair that needs
/// a HIGH-impact call of unknown frequency or with no GQ, DP or AD.
type Verdicts = (
    HashMap<AKey, Option<String>>,
    Option<String>,
    Option<Pair>,
    Option<String>,
);

/// `sf_verdicts`: for each distinct variant of a gene, `None` (reported) or
/// the reason it is not, the partial-deficiency note of a reported gene
/// with a hypomorph allele, and how a gene whose rule needs two hits met
/// it. `deletion` is the gene's pair with a called deletion, if any.
fn sf_verdicts(
    gene: &str,
    rule: &Rule,
    alleles: &IndexMap<AKey, Allele>,
    deletion: Option<&Entry>,
) -> Verdicts {
    let sev = |k: &AKey| match k {
        AKey::Small(k) => rule.severity.get(k),
        _ => None,
    };
    let mut out: IndexMap<AKey, Option<String>> = IndexMap::new();
    let mut kept: IndexMap<AKey, Allele> = IndexMap::new();
    for (k, cp) in alleles {
        match sev(k) {
            Some(a) if a.kind.as_deref() == Some("pseudodeficiency") => {
                out.insert(
                    k.clone(),
                    Some(format!(
                        "{}: a pseudodeficiency allele, not counted",
                        a.name
                    )),
                );
            }
            _ => {
                kept.insert(k.clone(), cp.clone());
            }
        }
    }
    let (n, proof, cis) = copies_hit(&kept);
    let ri = rule.report_if.as_str();
    let hemi = kept.values().any(|a| a.ploidy == 1);
    let counted = ri == "biallelic" || (ri == "hemizygous" && !hemi);
    let combined = kept.values().any(|a| a.combined);
    let one = if let Some(ps) = cis.first() {
        format!(
            "two heterozygous variants in cis: read-backed phase puts them on the same copy \
             (phase set {ps})"
        )
    } else if combined {
        "one heterozygous allele: a combined allele and its parts count once".to_owned()
    } else {
        "one heterozygous variant".to_owned()
    };
    for (k, a) in &kept {
        let good = match ri {
            "any" | "het_female_too" => true,
            "biallelic" => n >= 2,
            "hemizygous" => hemi || n >= 2,
            _ => {
                matches!(k, AKey::Small(s) if rule.variants.contains_key(s))
                    && a.copies >= 2
                    && a.ploidy == 2
            }
        };
        let why = if good {
            None
        } else if ri == "hom_variant" {
            Some(format!("{gene} is reported only for: {}", rule.table_text))
        } else if ri == "hemizygous" {
            Some(format!("{one}; {gene} is reported for {}", rule.table_text))
        } else {
            Some(format!(
                "{one}; {gene} is reported with {}",
                rule.table_text
            ))
        };
        out.insert(k.clone(), why);
    }
    // The count again with only the copies confident calls give: a pair
    // that needs a low-confidence call is only possible until the call is
    // confirmed. A small allele no confident call gives drops out; any other
    // takes the copies, GT and phase set of the confident call that gives
    // them, so a phase only a low call gives proves nothing. A deletion
    // counts as it is: its pair is judged by `two_hits_with`.
    let firm: IndexMap<AKey, Allele> = kept
        .iter()
        .filter(|(k, a)| !matches!(k, AKey::Small(_)) || a.firm > 0)
        .map(|(k, a)| {
            let mut f = a.clone();
            if matches!(k, AKey::Small(_)) {
                f.copies = a.firm;
                f.gt = a.firm_gt.clone();
                f.ps = a.firm_ps.clone();
            }
            (k.clone(), f)
        })
        .collect();
    let (n_firm, proof_firm, _) = copies_hit(&firm);
    let low = |why: String| {
        (
            Verdict::LowConfidence.word(),
            format!(
                "a call it rests on is low confidence (below the call floor, or a weak \
                 HIGH-impact call), so confirm it first; {why}"
            ),
        )
    };
    let unknown = || {
        "phase unknown: no read-backed phase links the variants, so they may sit on the same \
         copy (cis)"
            .to_owned()
    };
    // The best the evidence gives: a confident proof; then the deletion pair
    // when it is at least possibly biallelic (`two_hits_with` already lowers
    // a pair on a low call); then, with no deletion pair, two confident
    // copies of unknown phase; then the low word for a pair that needs a
    // low call.
    let pair = (counted && n >= 2).then(|| {
        let del = deletion.map(|e| (e.verdict, format!("with {}: {}", e.deletion, e.why)));
        match (proof_firm, proof, del) {
            (Some(f), _, _) => ("biallelic", f),
            (None, _, Some((v @ (Verdict::Biallelic | Verdict::PossiblyBiallelic), why))) => {
                (v.word(), why)
            }
            (None, Some(p), Some(_)) => low(p),
            (None, Some(_), None) if n_firm >= 2 => (Verdict::PossiblyBiallelic.word(), unknown()),
            (None, Some(p), None) => low(p),
            (None, None, Some((v, why))) => (v.word(), why),
            (None, None, None) if n_firm >= 2 => (Verdict::PossiblyBiallelic.word(), unknown()),
            (None, None, None) => low(unknown()),
        }
    });
    // A HIGH-impact call of unknown frequency (or with nothing to judge its
    // call by) counts toward the pair, with a note when the pair needs it:
    // the count again without such alleles.
    let sure: IndexMap<AKey, Allele> = kept
        .iter()
        .filter(|(_, a)| a.sure)
        .map(|(k, a)| (k.clone(), a.clone()))
        .collect();
    let pair_note = (pair.is_some() && copies_hit(&sure).0 < 2).then(|| {
        "the pair needs a HIGH-impact call of unknown frequency or with no GQ, DP or AD to \
         judge it: check that call (its frequency in gnomAD, its reads) before you count both \
         copies"
            .to_owned()
    });
    let hypo: Vec<&str> = out
        .iter()
        .filter(|(_, why)| why.is_none())
        .filter_map(|(k, _)| sev(k).filter(|a| a.kind.as_deref() == Some("hypomorph")))
        .map(|a| a.name.as_str())
        .collect();
    let partial = (!hypo.is_empty()).then(|| {
        format!(
            "partial, not profound: {} is a hypomorph allele, so the deficiency it causes is partial",
            hypo.join("; ")
        )
    });
    (out.into_iter().collect(), partial, pair, pair_note)
}

/// `CHROM:POS:REF:ALT` joined by `|`, as step 06 writes `PGP_JOINED_FROM`.
fn joined_keys(v: &str) -> Vec<Key> {
    v.split('|')
        .filter_map(|p| {
            let mut it = p.rsplitn(4, ':');
            let (alt, r, pos, chrom) = (it.next()?, it.next()?, it.next()?, it.next()?);
            let pos = crate::int(pos).ok()?;
            Some((chrom.to_owned(), pos, r.to_owned(), alt.to_owned()))
        })
        .collect()
}

/// The key each small candidate of a gene is counted under. A combined
/// allele and the parts of it that are candidates count once, under the
/// combined allele's key, when the combined call is phased or at most one
/// part is a candidate. An unphased combined allele with two or more parts
/// that are candidates counts through its parts, which may sit on
/// different copies.
fn counting_keys(cands: &[&Cand]) -> HashMap<(String, AKey), AKey> {
    let mut map = HashMap::new();
    let present: HashSet<(&str, &Key)> = cands
        .iter()
        .filter_map(|x| match &x.key {
            AKey::Small(k) => Some((x.gene.as_str(), k)),
            _ => None,
        })
        .collect();
    for x in cands.iter().filter(|x| !x.joined.is_empty()) {
        let parts: Vec<&Key> = x
            .joined
            .iter()
            .filter(|k| present.contains(&(x.gene.as_str(), *k)))
            .collect();
        let phased = Allele {
            gt: x.gt.clone(),
            ps: x.ps.clone(),
            ..Allele::of(x.copies, x.ploidy)
        }
        .phased()
        .is_some();
        if phased || parts.len() <= 1 {
            for p in parts {
                map.entry((x.gene.clone(), AKey::Small(p.clone())))
                    .or_insert_with(|| x.key.clone());
            }
        } else {
            map.entry((x.gene.clone(), x.key.clone()))
                .or_insert_with(|| AKey::Small(parts[0].clone()));
        }
    }
    map
}

/// `sf_structural`: step 05's class 4 and 5 SVs in genes whose rule allows
/// them, and the class 4 or 5 deletions over a gene's sv_regions. A common
/// SV is left out, but not one inside a much larger common SV, nor one with
/// no frequency data ([`SvRow::common`](super::two_hits::SvRow::common)).
fn sf_structural(
    inp: &DeletionInputs,
    rules: &IndexMap<String, Rule>,
) -> (Vec<Cand>, bool, Option<String>) {
    let rows = match &inp.annotsv {
        Err(e) => {
            return (
                Vec::new(),
                false,
                Some(format!("the AnnotSV table could not be read: {e}")),
            );
        }
        Ok(None) => return (Vec::new(), false, None),
        Ok(Some(rows)) => rows,
    };
    let mut out = Vec::new();
    for x in rows {
        if x.common || !matches!(x.class(), "4" | "5") {
            continue;
        }
        let (chrom, span) = x.sv.split_once(':').unwrap_or((&x.sv, ""));
        let (a, b) = span.split_once('-').unwrap_or((span, ""));
        let chrom = if chrom.starts_with("chr") {
            chrom.to_owned()
        } else {
            format!("chr{chrom}")
        };
        let mut hits: IndexMap<String, Option<String>> = IndexMap::new();
        for g in &x.genes {
            if rules.get(g).is_some_and(|r| r.sv_cn) {
                hits.insert(g.clone(), None);
            }
        }
        let sv_type = x.sv_type.as_deref().unwrap_or("None");
        let (an, bn) = (
            isdigit(a).then(|| crate::int(a).ok()).flatten(),
            isdigit(b).then(|| crate::int(b).ok()).flatten(),
        );
        for (g, rule) in rules {
            for (c, rs, re, name) in &rule.sv_regions {
                if let (Some(an), Some(bn)) = (an, bn)
                    && sv_type == "DEL"
                    && *c == chrom
                    && an <= *re
                    && bn >= *rs
                    && !hits.contains_key(g)
                {
                    hits.insert(g.clone(), Some(name.clone()));
                }
            }
        }
        hits.sort_keys();
        for (g, via) in hits {
            let mut f = IndexMap::new();
            f.insert("gene".into(), Value::Str(g.clone()));
            f.insert(
                "variant".into(),
                Value::Str(format!("{sv_type} {chrom}:{a}-{b}")),
            );
            f.insert("acmg_class".into(), Value::Str(x.class().to_owned()));
            f.insert("caller".into(), opt_s(x.caller.as_deref()));
            f.insert("support".into(), opt_int(x.support));
            f.insert("via".into(), opt_s(via.as_deref()));
            out.push(Cand {
                kind: Kind::Sv,
                gene: g,
                fields: f,
                key: AKey::Sv(
                    chrom.clone(),
                    a.to_owned(),
                    b.to_owned(),
                    sv_type.to_owned(),
                ),
                copies: 1,
                ploidy: 2,
                gt: String::new(),
                ps: None,
                joined: Vec::new(),
                apart: None,
            });
        }
    }
    (out, true, None)
}

/// The disease ACMG lists a gene for: its name (the `phenotype` column of
/// `assets/acmg_sf.tsv`) and the keys a ClinVar condition is matched by
/// (`phenotype_keys`): database ids such as `OMIM:171400`, Table 1's MIM
/// disorder numbers, and phrases looked for in ClinVar's condition names.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Phenotype {
    pub name: String,
    pub ids: Vec<String>,
    pub phrases: Vec<String>,
}

/// Database prefixes a key of `phenotype_keys` can be an id of.
const ID_PREFIXES: [&str; 4] = ["OMIM:", "MedGen:", "MONDO:", "Orphanet:"];

/// Condition names that ClinVar gives a whole panel of genes, not one
/// disease (normalized as [`words`] writes them). They say nothing for or
/// against the ACMG disease, so they never set the mark.
const PANEL_TERMS: [&str; 7] = [
    "not provided",
    "not specified",
    "see cases",
    "inborn genetic diseases",
    "hereditary cancer predisposing syndrome",
    "cardiovascular phenotype",
    "cardiac arrhythmia",
];

/// What follows the gene in ClinVar's name for every disease of one gene,
/// such as "RYR1-related disorder" (normalized as [`words`] writes them).
const UMBRELLA_TAILS: [&str; 5] = [
    "related disorder",
    "related disorders",
    "related condition",
    "related conditions",
    "related cancer predisposition",
];

/// Whether `name` only names `gene`, as "RYR1-related disorder" does on a
/// RYR1 hit. Like a panel term it says nothing for or against the ACMG
/// disease. Another gene's such name ("FBXO11-related disorders" on an
/// MSH6 hit) is not this.
fn umbrella(gene: &str, name: &str) -> bool {
    let (g, n) = (words(gene), words(name));
    !g.is_empty()
        && n.strip_prefix(g.as_str())
            .and_then(|tail| tail.strip_prefix(' '))
            .is_some_and(|tail| UMBRELLA_TAILS.contains(&tail))
}

/// A name as the matching reads it: lower case, and every run of
/// characters other than letters and digits one space.
pub fn words(s: &str) -> String {
    s.to_lowercase()
        .split(|c: char| !c.is_alphanumeric())
        .filter(|w| !w.is_empty())
        .collect::<Vec<_>>()
        .join(" ")
}

impl Phenotype {
    /// Reads one row's two cells; `None` when the row has no keys (an
    /// `acmg_sf.tsv` without the columns, or a `.` cell).
    fn of(name: Option<&str>, keys: Option<&str>) -> Option<Phenotype> {
        let mut p = Phenotype {
            name: name
                .filter(|n| !matches!(*n, "" | "."))
                .unwrap_or("")
                .to_owned(),
            ..Phenotype::default()
        };
        for k in keys.unwrap_or("").split(';').map(str::trim) {
            if k.is_empty() || k == "." {
                continue;
            }
            if ID_PREFIXES.iter().any(|p| k.starts_with(p)) {
                p.ids.push(k.to_owned());
            } else {
                p.phrases.push(words(k));
            }
        }
        (!p.ids.is_empty() || !p.phrases.is_empty()).then_some(p)
    }

    /// Whether ClinVar's condition is this disease: one of its ids is a key,
    /// or a key phrase starts a word of its name ("multiple endocrine
    /// neoplasia type 2" is in "multiple endocrine neoplasia type 2a").
    pub fn matches(&self, c: &Condition) -> bool {
        if c.ids.iter().any(|id| self.ids.contains(id)) {
            return true;
        }
        let name = words(&c.name);
        self.phrases.iter().any(|p| {
            name.match_indices(p.as_str())
                .any(|(i, _)| i == 0 || name.as_bytes()[i - 1] == b' ')
        })
    }

    /// The mark on a hit of `gene`: every condition ClinVar names for the
    /// allele is another disease. A placeholder, a panel term or a name
    /// that only names `gene` counts neither way, so a hit whose only
    /// condition is "not provided", "Hereditary cancer-predisposing
    /// syndrome" or "RYR1-related disorder" (on RYR1) is not marked. A
    /// condition with ids and no name never sets the mark, but its ids can
    /// still be the ACMG disease.
    pub fn other_disease(&self, gene: &str, conditions: &[Condition]) -> bool {
        let other = conditions.iter().any(|c| {
            !c.placeholder()
                && !PANEL_TERMS.contains(&words(&c.name).as_str())
                && !umbrella(gene, &c.name)
        });
        other && !conditions.iter().any(|c| self.matches(c))
    }
}

/// `{gene: phenotype}` of `acmg_sf.tsv`: the genes whose row has keys.
pub fn phenotypes(assets: &str) -> Result<HashMap<String, Phenotype>, String> {
    let mut out = HashMap::new();
    for r in asset_rows(assets, "acmg_sf.tsv")? {
        let cell = |k: &str| r.get(k).flatten();
        if let (Some(g), Some(p)) = (
            cell("gene"),
            Phenotype::of(cell("phenotype"), cell("phenotype_keys")),
        ) {
            out.insert(g.to_owned(), p);
        }
    }
    Ok(out)
}

fn zyg_map(gt: &str) -> String {
    match gt {
        "0/1" | "1/0" | "0|1" | "1|0" => "het".into(),
        "1/1" | "1|1" => "hom".into(),
        _ => gt.into(),
    }
}

/// What the ACMG tier reads.
pub struct Inputs<'a> {
    pub assets: &'a str,
    /// Step 06's floor and hits, or the note of a hits VCF that cannot be
    /// read; `None` without one.
    pub hits: Option<&'a Result<HitRecords, String>>,
    /// The clinical VCF, or its note; `None` without one.
    pub clinical: Option<&'a Result<ClinicalVcf, String>>,
    pub deletions: &'a DeletionInputs,
}

/// `acmg_sf_or_note`: the tier, or `Ok` of the unavailable form when a
/// list or an input cannot be read. An `Err` is the collector's note of the
/// whole clinical section.
pub fn acmg_sf(i: &Inputs) -> Result<Value, String> {
    match tier(i)? {
        Ok(v) => Ok(v),
        Err(why) => Ok(obj([
            ("version", Value::Str(ACMG_SF_VERSION.into())),
            ("unavailable", Value::Str(why)),
            ("clinvar_hits", Value::Array(Vec::new())),
            ("high_impact", Value::Array(Vec::new())),
            ("structural", Value::Array(Vec::new())),
            ("not_reported", Value::Array(Vec::new())),
        ])),
    }
}

/// `acmg_sf_tier`: the outer error is the section's note, the inner one
/// the reason the tier is unavailable.
fn tier(i: &Inputs) -> Result<Result<Value, String>, String> {
    let rules = match acmg_sf_rules(i.assets) {
        Ok(r) => r,
        Err(e) => return Ok(Err(e)),
    };
    // The list was read just above, so its phenotypes read too.
    let pheno = phenotypes(i.assets).unwrap_or_default();
    let mut clinvar: Vec<Cand> = Vec::new();
    if let Some(h) = i.hits {
        let (floor, recs) = match h {
            Ok(x) => x,
            Err(e) => return Ok(Err(e.clone())),
        };
        for r in recs.iter().filter(|r| r.len() >= 8) {
            let row = Hit::parse(r, *floor)?;
            let gt = text::gt_of(r).to_owned();
            let (copies, ploidy) = allele_copies(&gt);
            let joined = joined_keys(
                info_map(&r[7])
                    .get("PGP_JOINED_FROM")
                    .copied()
                    .unwrap_or(""),
            );
            for g in row.gene.split(',') {
                if !rules.contains_key(g) {
                    continue;
                }
                let mut f = IndexMap::new();
                f.insert("gene".into(), Value::Str(g.to_owned()));
                f.insert("variant".into(), Value::Str(row.variant()));
                f.insert("genotype".into(), Value::Str(row.genotype.clone()));
                f.insert("significance".into(), Value::Str(row.significance.clone()));
                f.insert("stars".into(), Value::Int(row.stars.into()));
                f.insert("low_confidence".into(), Value::Bool(row.low_confidence));
                if row.combined_allele != "." {
                    f.insert(
                        "combined_allele".into(),
                        Value::Str(row.combined_allele.clone()),
                    );
                }
                // The condition ClinVar names, and whether it is another
                // disease than the one ACMG lists the gene for: GENEINFO
                // alone puts a RET Hirschsprung hit in the MEN2 slot. The
                // hit stays a finding either way; the mark only says so.
                if let Some(c) = &row.conditions {
                    for (k, v) in condition_fields(c) {
                        f.insert(k.into(), v);
                    }
                    if let Some(p) = pheno.get(g) {
                        f.insert("acmg_phenotype".into(), Value::Str(p.name.clone()));
                        f.insert(
                            "condition_not_acmg".into(),
                            Value::Bool(p.other_disease(g, c)),
                        );
                    }
                }
                clinvar.push(Cand {
                    kind: Kind::ClinVar,
                    gene: g.to_owned(),
                    fields: f,
                    key: AKey::Small(row.key()),
                    copies,
                    ploidy,
                    gt: gt.clone(),
                    ps: gt_field(r, "PS"),
                    joined: joined.clone(),
                    apart: None,
                });
            }
        }
    }
    let mut high: Vec<Cand> = Vec::new();
    if let Some(c) = i.clinical {
        let c = match c {
            Ok(c) => c,
            Err(e) => return Ok(Err(e.clone())),
        };
        let mut recs: Vec<(&Key, &HighRecord)> =
            c.high.as_ref().map_err(Clone::clone)?.iter().collect();
        recs.sort_by(|a, b| a.0.cmp(b.0));
        for ((ch, pos, r, a), x) in recs {
            if x.frequency == "common" {
                continue;
            }
            let mut genes: Vec<(&String, &String)> = x.genes.iter().collect();
            genes.sort();
            for (g, cons) in genes {
                if !rules.contains_key(g) {
                    continue;
                }
                let (copies, ploidy) = allele_copies(&x.genotype);
                let mut f = IndexMap::new();
                f.insert("gene".into(), Value::Str(g.clone()));
                f.insert("variant".into(), Value::Str(format!("{ch}:{pos} {r}>{a}")));
                f.insert("genotype".into(), Value::Str(zyg_map(&x.genotype)));
                f.insert("consequence".into(), Value::Str(cons.clone()));
                f.insert("transcript".into(), Value::Str(x.transcript.into()));
                f.insert("frequency".into(), Value::Str(x.frequency.into()));
                f.insert(
                    "frequency_unknown".into(),
                    Value::Bool(x.frequency == "unknown"),
                );
                f.insert("weak_call".into(), Value::Bool(x.weak));
                if x.unverified() {
                    f.insert("call_unverified".into(), Value::Bool(true));
                }
                f.insert(
                    "call_note".into(),
                    if x.why.is_empty() {
                        Value::Null
                    } else {
                        Value::Str(x.why.join("; "))
                    },
                );
                f.insert("gq".into(), opt_int(x.gq));
                f.insert("dp".into(), opt_int(x.dp));
                f.insert("ad".into(), opt_s(x.ad.as_deref()));
                high.push(Cand {
                    kind: Kind::High,
                    gene: g.clone(),
                    fields: f,
                    key: AKey::Small((ch.clone(), *pos, r.clone(), a.clone())),
                    copies,
                    ploidy,
                    gt: x.genotype.clone(),
                    ps: x.ps.clone(),
                    joined: Vec::new(),
                    apart: None,
                });
            }
        }
    }
    // Each gene's mechanism and footnote: a HIGH-impact call in a gene where
    // loss of function is not the disease mechanism, or a TTN variant that
    // is no truncating one in the exons with high PSI, is listed apart and
    // counts nothing. A ClinVar hit reads its consequence from its record
    // in the clinical VCF.
    let clinical_vcf = match i.clinical {
        Some(Ok(c)) if c.high.is_ok() => Some(c),
        _ => None,
    };
    let clinical_high = clinical_vcf.and_then(|c| c.high.as_ref().ok());
    let mut conditions: HashMap<(Kind, String, Key), String> = HashMap::new();
    for x in clinvar.iter_mut().chain(high.iter_mut()) {
        let AKey::Small(k) = &x.key else { continue };
        let is_high = x.kind == Kind::High;
        let consequence = if is_high {
            match x.fields.get("consequence") {
                Some(Value::Str(s)) => Some(s.clone()),
                _ => None,
            }
        } else {
            clinical_high
                .and_then(|h| h.get(k))
                .and_then(|r| r.genes.get(&x.gene))
                .cloned()
        };
        let held = if is_high {
            Some(true)
        } else {
            clinical_vcf.map(|c| c.held.contains(k))
        };
        match gate(
            &x.gene,
            &rules[&x.gene],
            k,
            is_high,
            consequence.as_deref(),
            held,
        ) {
            Err(why) => x.apart = Some(why),
            Ok(Some(c)) => {
                conditions.insert((x.kind, x.gene.clone(), k.clone()), c);
            }
            Ok(None) => {}
        }
    }
    let (svs, sv_checked, sv_note) = sf_structural(i.deletions, &rules);

    // One verdict per gene, over its distinct variants. A record with no
    // alt copy called is no evidence of a hit and is not counted.
    let small: Vec<&Cand> = clinvar
        .iter()
        .chain(&high)
        .filter(|x| x.copies > 0 && x.apart.is_none())
        .collect();
    let units = counting_keys(&small);
    let unit = |x: &Cand| {
        units
            .get(&(x.gene.clone(), x.key.clone()))
            .cloned()
            .unwrap_or_else(|| x.key.clone())
    };
    let mut alleles: IndexMap<String, IndexMap<AKey, Allele>> = IndexMap::new();
    // The candidates counted under their own key first, so a combined
    // allele's own call stands for its parts.
    let (own, folded): (Vec<&Cand>, Vec<&Cand>) =
        small.iter().copied().partition(|x| unit(x) == x.key);
    // A low-confidence call: a step 06 hit below its call floor (the
    // ClinVar row's low_confidence), or a weak HIGH-impact call (its
    // weak_call: the clinical call floor or allele balance). Neither makes
    // its allele firm.
    let is_low = |x: &Cand| {
        matches!(x.fields.get("low_confidence"), Some(Value::Bool(true)))
            || matches!(x.fields.get("weak_call"), Some(Value::Bool(true)))
    };
    for x in &own {
        let p = if x.ploidy == 0 { 2 } else { x.ploidy };
        let a = alleles
            .entry(x.gene.clone())
            .or_default()
            .entry(x.key.clone())
            .or_insert(Allele {
                gt: x.gt.clone(),
                ps: x.ps.clone(),
                sure: false,
                ..Allele::of(x.copies, p)
            });
        if !is_low(x) {
            a.confirm(x.copies, &x.gt, &x.ps);
        }
        a.sure |= !x.unsure();
    }
    for x in &folded {
        let p = if x.ploidy == 0 { 2 } else { x.ploidy };
        let a = alleles
            .entry(x.gene.clone())
            .or_default()
            .entry(unit(x))
            .or_insert(Allele {
                gt: x.gt.clone(),
                ps: x.ps.clone(),
                sure: false,
                ..Allele::of(x.copies, p)
            });
        a.copies = a.copies.max(x.copies);
        a.combined = true;
        if !is_low(x) {
            a.confirm(x.copies, &x.gt, &x.ps);
        }
        a.sure |= !x.unsure();
    }
    for x in &svs {
        alleles
            .entry(x.gene.clone())
            .or_default()
            .entry(x.key.clone())
            .or_insert(Allele::of(1, 2));
    }
    let (spans, two_hit_note): (Spans, Option<String>) = match carrier_genes(i.assets) {
        Ok(g) => (g.into_iter().map(|(k, v)| (k, v.span)).collect(), None),
        Err(e) => (
            IndexMap::new(),
            Some(format!(
                "A deletion as the second hit: only AnnotSV's gene names were read; step 35's \
                 deletions need the gene spans of {e}"
            )),
        ),
    };
    let (dels, _, _) = gene_deletions(i.deletions, &spans);
    // The confident call of a small allele whose own call is low, when it
    // gives as many copies: a combined allele counts as confident when any
    // call of it is, and its pair with a deletion then rests on that call's
    // GT and phase set, not on the low call's.
    let confirmed = |x: &Cand| {
        alleles
            .get(&x.gene)
            .and_then(|a| a.get(&x.key))
            .filter(|a| is_low(x) && a.firm >= x.copies)
            .map(|a| (a.firm_gt.clone(), a.firm_ps.clone()))
    };
    let mut seen = HashSet::new();
    let mut variants = Vec::new();
    for x in &own {
        let AKey::Small(k) = &x.key else { continue };
        let pseudo = rules[&x.gene].severity_kind(k) == Some("pseudodeficiency");
        if !pseudo && seen.insert((x.gene.clone(), k.clone())) {
            let (gt, ps) = confirmed(x).unwrap_or_else(|| (x.gt.clone(), x.ps.clone()));
            variants.push(Variant {
                gene: x.gene.clone(),
                key: k.clone(),
                chrom: k.0.clone(),
                pos: k.1,
                text: match &x.fields["variant"] {
                    Value::Str(s) => s.clone(),
                    _ => String::new(),
                },
                genotype: match &x.fields["genotype"] {
                    Value::Str(s) => s.clone(),
                    _ => String::new(),
                },
                copies: x.copies,
                ploidy: if x.ploidy == 0 { 2 } else { x.ploidy },
                gt,
                ps,
                allele: None,
            });
        }
    }
    let low: HashSet<(String, Key)> = own
        .iter()
        .filter(|x| is_low(x) && confirmed(x).is_none())
        .filter_map(|x| match &x.key {
            AKey::Small(k) => Some((x.gene.clone(), k.clone())),
            _ => None,
        })
        .collect();
    let mut found = two_hits_with(&variants, &dels, |v| {
        low.contains(&(v.gene.clone(), v.key.clone()))
    });
    // The card's low calls include a weak HIGH-impact call, so its pairs
    // name both kinds, as its own low-confidence word does.
    let (floor, either) = (
        "the variant's call is below the call floor (low confidence)",
        "the variant's call is low confidence (below the call floor, or a weak HIGH-impact call)",
    );
    for e in found.entries.values_mut() {
        e.why = e.why.replacen(floor, either, 1);
    }
    for t in found.pairs.values_mut() {
        *t = t.replacen(floor, either, 1);
    }
    let mut two: HashMap<(String, Key), String> = HashMap::new();
    for g in found.entries.keys() {
        if matches!(rules[g].report_if.as_str(), "biallelic" | "hemizygous") {
            if let Some(a) = alleles.get_mut(g) {
                a.insert(AKey::Deletion(g.clone()), Allele::of(1, 2));
            }
            two.extend(
                found
                    .pairs
                    .iter()
                    .filter(|(k, _)| &k.0 == g)
                    .map(|(k, t)| (k.clone(), t.clone())),
            );
        }
    }
    let verdict: HashMap<&String, _> = alleles
        .iter()
        .map(|(g, a)| {
            let del = found
                .entries
                .get(g)
                .filter(|_| a.contains_key(&AKey::Deletion(g.clone())));
            (g, sf_verdicts(g, &rules[g], a, del))
        })
        .collect();
    let mut reported: [Vec<Value>; 3] = Default::default();
    let mut not_reported = Vec::new();
    for x in clinvar.into_iter().chain(high).chain(svs) {
        if x.copies == 0 {
            not_reported.push(obj([
                ("gene", Value::Str(x.gene.clone())),
                ("variant", x.fields["variant"].clone()),
                (
                    "genotype",
                    x.fields.get("genotype").cloned().unwrap_or(Value::Null),
                ),
                ("source", Value::Str(x.kind.source().into())),
                (
                    "reason",
                    Value::Str(format!("no alt copy called (GT {}): not counted", x.gt)),
                ),
            ]));
            continue;
        }
        if let Some(reason) = &x.apart {
            not_reported.push(obj([
                ("gene", Value::Str(x.gene.clone())),
                ("variant", x.fields["variant"].clone()),
                (
                    "genotype",
                    x.fields.get("genotype").cloned().unwrap_or(Value::Null),
                ),
                ("source", Value::Str(x.kind.source().into())),
                ("reason", Value::Str(reason.clone())),
            ]));
            continue;
        }
        let key = unit(&x);
        let (whys, partial, pair, pair_note) = &verdict[&x.gene];
        let why = whys.get(&key).cloned().flatten();
        let rule = &rules[&x.gene];
        let mut f = x.fields;
        if let AKey::Small(k) = &key
            && found.hemi.contains_key(&(x.gene.clone(), k.clone()))
        {
            f.insert("genotype".into(), Value::Str("hemizygous".into()));
        }
        let slot = match x.kind {
            Kind::ClinVar => 0,
            Kind::High => 1,
            Kind::Sv => 2,
        };
        match why {
            None => {
                f.insert("rule".into(), Value::Str(rule.table_text.clone()));
                f.insert("rule_note".into(), opt_s(rule.note.as_deref()));
                if let AKey::Small(k) = &x.key
                    && let Some(c) = conditions.get(&(x.kind, x.gene.clone(), k.clone()))
                {
                    f.insert("rule_condition".into(), Value::Str(c.clone()));
                }
                if x.kind == Kind::High && rule.lof == Lof::Mark {
                    f.insert(
                        "mechanism_note".into(),
                        Value::Str(format!(
                            "{} counts a predicted loss-of-function call with a mark: its \
                             disease is mostly caused by missense variants, so a truncating \
                             one is a less certain cause; check it with a clinician",
                            x.gene
                        )),
                    );
                }
                f.insert("severity".into(), opt_s(partial.as_deref()));
                if x.kind != Kind::Sv
                    && let AKey::Small(k) = &key
                    && let Some(t) = two.get(&(x.gene.clone(), k.clone()))
                {
                    f.insert("two_hit".into(), Value::Str(t.clone()));
                }
                if let Some((word, phase)) = pair {
                    f.insert("pair".into(), Value::Str((*word).into()));
                    f.insert("phase".into(), Value::Str(phase.clone()));
                    if let Some(n) = pair_note {
                        f.insert("pair_note".into(), Value::Str(n.clone()));
                    }
                }
                reported[slot].push(Value::Object(f));
            }
            Some(reason) => not_reported.push(obj([
                ("gene", Value::Str(x.gene.clone())),
                ("variant", f["variant"].clone()),
                (
                    "genotype",
                    f.get("genotype").cloned().unwrap_or(Value::Null),
                ),
                ("source", Value::Str(x.kind.source().into())),
                ("reason", Value::Str(reason)),
            ])),
        }
    }
    let [clinvar_rows, high_rows, sv_rows] = reported;
    let flag =
        |x: &Value, k: &str| matches!(x, Value::Object(m) if m.get(k) == Some(&Value::Bool(true)));
    let sure = high_rows
        .iter()
        .filter(|x| {
            !(flag(x, "weak_call") || flag(x, "frequency_unknown") || flag(x, "call_unverified"))
        })
        .count();
    let unverified = high_rows
        .iter()
        .filter(|x| flag(x, "call_unverified"))
        .count();
    let weak = high_rows.iter().filter(|x| flag(x, "weak_call")).count();
    let unknown = high_rows
        .iter()
        .filter(|x| flag(x, "frequency_unknown"))
        .count();
    // A HIGH-impact call with nothing to judge it by is counted apart, and
    // only a run that has one gets the count.
    let mut tail = vec![("high_impact_frequency_unknown", Value::Int(unknown as i128))];
    if unverified > 0 {
        tail.push(("high_impact_unverified", Value::Int(unverified as i128)));
    }
    Ok(Ok(obj([
        ("version", Value::Str(ACMG_SF_VERSION.into())),
        ("genes_on_list", Value::Int(rules.len() as i128)),
        ("rules", Value::Str("assets/acmg_sf.tsv".into())),
        ("clinvar_hits", Value::Array(clinvar_rows)),
        ("clinvar_checked", Value::Bool(i.hits.is_some())),
        ("high_impact", Value::Array(high_rows)),
        ("high_impact_checked", Value::Bool(i.clinical.is_some())),
        ("high_impact_rare", Value::Int(sure as i128)),
        ("high_impact_weak", Value::Int(weak as i128)),
    ]
    .into_iter()
    .chain(tail)
    .chain([
        ("structural", Value::Array(sv_rows)),
        ("structural_checked", Value::Bool(sv_checked)),
        ("structural_note", opt_s(sv_note.as_deref())),
        ("not_reported", Value::Array(not_reported)),
        ("two_hit_note", opt_s(two_hit_note.as_deref())),
    ]))))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fmt(s: &str) -> Vec<String> {
        s.split('|').map(str::to_owned).collect()
    }

    #[test]
    fn mane_then_canonical_then_any() {
        let f = fmt("Allele|IMPACT|SYMBOL|Consequence|CANONICAL|MANE_SELECT|gnomADe_AF");
        let (c, e, q) = csq_rank(
            &f,
            "T|HIGH|A|stop_gained|YES||0.5,T|LOW|A|synonymous_variant||NM_1|0.001",
        );
        assert_eq!((c, q), ("MANE", "common"));
        assert_eq!(e.len(), 1);
        assert_eq!(e[0].impact, "LOW");
        let (c, _, q) = csq_rank(&f, "T|HIGH|A|stop_gained&splice|YES||");
        assert_eq!((c, q), ("CANONICAL", "unknown"));
        let f2 = fmt("Allele|IMPACT|SYMBOL|Consequence");
        let (c, e, q) = csq_rank(&f2, "T|HIGH|A|a&b");
        assert_eq!((c, q), ("ANY", "none"));
        assert_eq!(e[0].consequence, "a,b");
    }

    #[test]
    fn the_genotype_check() {
        let f = ClinicalFloor::default();
        let (weak, why, gq, dp, ad) = call_check("GT:GQ:DP:AD", "0/1:20:30:25,5", &f);
        assert!(weak);
        assert_eq!(
            why,
            [
                "GQ 20 < 30",
                "alt reads 0.17 of the depth on a heterozygote"
            ]
        );
        assert_eq!((gq, dp, ad.as_deref()), (Some(20), Some(30), Some("25,5")));
        let (weak, why, ..) = call_check("GT:AD", "1|1:5,15", &f);
        assert!(weak);
        assert_eq!(why, ["alt reads 0.75 of the depth on a homozygous call"]);
        let (weak, ..) = call_check("GT:GQ:DP:AD", "0/1:40:30:15,15", &f);
        assert!(!weak);
    }

    #[test]
    fn a_condition_matches_by_id_or_by_a_phrase_at_a_word_start() {
        let p = Phenotype::of(
            Some("MEN2"),
            Some("OMIM:171400; multiple endocrine neoplasia type 2; men2"),
        )
        .unwrap();
        let c = |name: &str, ids: &[&str]| Condition {
            name: name.into(),
            ids: ids.iter().map(|s| (*s).to_owned()).collect(),
        };
        assert!(p.matches(&c("Multiple endocrine neoplasia, type 2a", &[])));
        assert!(p.matches(&c("MEN2B", &[])));
        assert!(p.matches(&c("Something else", &["OMIM:171400"])));
        assert!(!p.matches(&c("Hirschsprung disease 1", &["OMIM:142623"])));
        // A phrase inside a word is no match.
        assert!(!p.matches(&c("xmen2", &[])));
        let hd = [c(
            "Hirschsprung disease, susceptibility to, 1",
            &["OMIM:142623"],
        )];
        assert!(p.other_disease("RET", &hd));
        // A placeholder or a panel term counts neither way.
        assert!(!p.other_disease("RET", &[c("not provided", &[])]));
        assert!(!p.other_disease("RET", &[c("Hereditary cancer-predisposing syndrome", &[])]));
        assert!(!p.other_disease("RET", &[c("Cardiac arrhythmia", &[])]));
        assert!(p.other_disease("RET", &[hd[0].clone(), c("not provided", &[])]));
        // One condition of the ACMG disease is enough.
        assert!(!p.other_disease("RET", &[hd[0].clone(), c("MEN2A", &[])]));
        assert!(!p.other_disease("RET", &[]));
        // A name that only names the hit's own gene counts neither way;
        // another gene's such name, or a disease named after it, does.
        for own in [
            "RET-related disorder",
            "RET-related disorders",
            "RET related condition",
            "RET-related conditions",
            "RET-related cancer predisposition",
        ] {
            assert!(!p.other_disease("RET", &[c(own, &[])]), "{own}");
        }
        assert!(p.other_disease("RET", &[c("GLA-related disorder", &[])]));
        assert!(p.other_disease("RET", &[c("RET-related myopathy", &[])]));
        assert!(p.other_disease("RET", &[c("RET-related disorder", &[]), hd[0].clone()]));
        assert!(!umbrella("", "related disorder"));
        // A condition with ids and no name never sets the mark, and its ids
        // can still be the ACMG disease.
        assert!(!p.other_disease("RET", &[c("", &["OMIM:142623"])]));
        assert!(!p.other_disease("RET", &[hd[0].clone(), c("", &["OMIM:171400"])]));
        assert!(p.other_disease("RET", &[hd[0].clone(), c("", &["OMIM:1"])]));
        assert_eq!(Phenotype::of(Some("x"), Some(".")), None);
        assert_eq!(
            words("Ehlers-Danlos syndrome, type IV"),
            "ehlers danlos syndrome type iv"
        );
    }

    fn rec(gt: &str, gq: Option<i128>, dp: Option<i128>, ad: Option<&str>) -> HighRecord {
        HighRecord {
            genotype: gt.into(),
            ps: None,
            genes: IndexMap::new(),
            transcript: "MANE",
            frequency: "rare",
            weak: false,
            why: Vec::new(),
            gq,
            dp,
            ad: ad.map(str::to_owned),
        }
    }

    /// A call with no GQ, DP or AD read total is unverified, as step 23's
    /// `CALL_CHECK=unknown`; any one of them makes it judged.
    #[test]
    fn a_call_with_nothing_to_judge_it_by_is_unverified() {
        assert!(rec("0/1", None, None, None).unverified());
        assert!(rec("0/1", None, None, Some(".")).unverified());
        assert!(rec("0/1", None, None, Some("0,0")).unverified());
        assert!(rec("1", None, None, Some("3,4")).unverified());
        assert!(!rec("0/1", Some(50), None, None).unverified());
        assert!(!rec("0/1", None, Some(30), None).unverified());
        assert!(!rec("0|1", None, None, Some("15,15")).unverified());
        let mut w = rec("0/1", None, None, None);
        w.weak = true;
        assert!(!w.unverified());
    }

    #[test]
    fn the_gate_of_a_gene_mechanism_and_footnote() {
        let rule = |lof, limit| Rule {
            report_if: "any".into(),
            variants: IndexMap::new(),
            sv_cn: false,
            sv_regions: Vec::new(),
            severity: IndexMap::new(),
            table_text: "All P and LP".into(),
            note: None,
            lof,
            limit,
        };
        let k = |chrom: &str, pos| (chrom.to_owned(), pos, "C".to_owned(), "T".to_owned());
        let no = rule(Lof::No, Limit::None);
        assert!(
            gate(
                "PCSK9",
                &no,
                &k("chr1", 5),
                true,
                Some("stop_gained"),
                Some(true)
            )
            .is_err()
        );
        assert_eq!(
            gate("PCSK9", &no, &k("chr1", 5), false, None, Some(false)),
            Ok(None)
        );
        let ttn = rule(
            Lof::Yes,
            Limit::TruncatingIn {
                regions: vec![("chr2".into(), 100, 200), ("chr2".into(), 300, 400)],
                name: "the TTN exons with PSI over 90%".into(),
            },
        );
        let stop = |pos| {
            gate(
                "TTN",
                &ttn,
                &k("chr2", pos),
                true,
                Some("stop_gained"),
                Some(true),
            )
        };
        // Each end of each region counts; the base beyond it does not.
        for pos in [100, 150, 200, 300, 400] {
            assert_eq!(stop(pos), Ok(None), "{pos}");
        }
        for pos in [99, 201, 250, 299, 401] {
            let e = stop(pos).unwrap_err();
            assert!(
                e.contains("outside") && e.contains("footnote f"),
                "{pos}: {e}"
            );
        }
        // A record counts by the bases its REF spans, not by its POS alone: a
        // frameshift deletion that starts 3 bp before a region and deletes
        // into it (its canonical splice site) counts; one whose REF ends a
        // base before the region does not.
        let del = |pos, r: &str| {
            gate(
                "TTN",
                &ttn,
                &("chr2".to_owned(), pos, r.to_owned(), "C".to_owned()),
                true,
                Some("frameshift_variant"),
                Some(true),
            )
        };
        assert_eq!(del(97, "CAAAAA"), Ok(None), "97 CAAAAA reaches 101");
        assert_eq!(del(97, "CAAA"), Ok(None), "97 CAAA reaches 100");
        assert!(del(97, "CAA").is_err(), "97 CAA ends at 99");
        assert_eq!(del(399, "CAAAAA"), Ok(None), "399 CAAAAA starts inside");
        assert!(del(401, "CAAAAA").is_err(), "401 starts past the end");
        assert_eq!(
            gate(
                "TTN",
                &ttn,
                &k("2", 150),
                true,
                Some("frameshift_variant"),
                Some(true)
            ),
            Ok(None)
        );
        assert!(
            gate(
                "TTN",
                &ttn,
                &k("chr3", 150),
                true,
                Some("stop_gained"),
                Some(true)
            )
            .is_err()
        );
        assert!(
            gate(
                "TTN",
                &ttn,
                &k("chr2", 150),
                true,
                Some("start_lost"),
                Some(true)
            )
            .is_err()
        );
        // A ClinVar hit: its record in the clinical VCF decides; a record
        // that does not truncate sets it apart, no record at all (or no
        // clinical VCF) reports it with the condition to check.
        let e = gate("TTN", &ttn, &k("chr2", 150), false, None, Some(true)).unwrap_err();
        assert!(e.contains("does not truncate"), "{e}");
        assert_eq!(
            gate(
                "TTN",
                &ttn,
                &k("chr2", 150),
                false,
                Some("stop_gained,splice_region_variant"),
                Some(true)
            ),
            Ok(None)
        );
        assert!(
            matches!(gate("TTN", &ttn, &k("chr2", 150), false, None, Some(false)), Ok(Some(c)) if c.contains("holds no record") && c.contains("not known"))
        );
        assert!(
            matches!(gate("TTN", &ttn, &k("chr2", 150), false, None, None), Ok(Some(c)) if c.contains("no clinical VCF") && c.contains("not known"))
        );
        // Outside the regions even a ClinVar hit with no record is apart.
        assert!(gate("TTN", &ttn, &k("chr2", 250), false, None, Some(false)).is_err());
        let lmna = rule(Lof::Yes, Limit::MissenseNeedsCaseEvidence);
        assert!(
            matches!(gate("LMNA", &lmna, &k("chr1", 1), false, None, Some(false)), Ok(Some(c)) if c.contains("footnote e"))
        );
        assert_eq!(
            gate(
                "LMNA",
                &lmna,
                &k("chr1", 1),
                false,
                Some("stop_gained"),
                Some(true)
            ),
            Ok(None)
        );
        assert_eq!(
            gate(
                "LMNA",
                &lmna,
                &k("chr1", 1),
                true,
                Some("start_lost"),
                Some(true)
            ),
            Ok(None)
        );
    }

    #[test]
    fn the_floors_of_the_header() {
        let l = |s: &str| vec![s.to_owned(), "#CHROM".to_owned()];
        let f = ClinicalFloor::read(&l(
            "##PGP_CALL_FLOOR=GQ=20,DP=5,HET_AB=0.25-0.75,HOM_AB=0.9",
        ))
        .unwrap();
        assert_eq!(f.call, Floor { gq: 20, dp: 5 });
        assert_eq!((f.het_ab, f.hom_ab), ((0.25, 0.75), 0.9));
        assert!(ClinicalFloor::read(&l("##PGP_CALL_FLOOR=HOM_AB=0..9")).is_err());
    }
}
