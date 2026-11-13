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

use super::assets::{ACMG_SF_VERSION, Key, Rule, acmg_sf_rules, carrier_genes};
use super::clinvar::{
    COMMON_AF, Floor, Hit, HitRecords, allele_copies, floor_values, obj, opt_int, opt_s,
};
use super::text::{self, digits, format_map, gt_field, info_map};
use super::two_hits::{DeletionInputs, Spans, Variant, gene_deletions, two_hits};
use indexmap::IndexMap;
use pycompat::json::Value;
use pycompat::text::isdigit;
use std::collections::{HashMap, HashSet};

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

/// The clinical VCF read once: its records count and its HIGH-impact
/// records. An error is the note of a file that cannot be read.
#[derive(Debug, Clone)]
pub struct ClinicalVcf {
    pub records: usize,
    /// Keyed by `(chrom, pos, ref, alt)`; an error is the collector's note
    /// (a POS that is no integer on a HIGH-impact record, a bad floor).
    pub high: Result<IndexMap<Key, HighRecord>, String>,
}

impl ClinicalVcf {
    pub fn read(path: &str) -> Result<ClinicalVcf, String> {
        let lines = text::lines(path).map_err(|e| e.unreadable(path))?;
        let records = lines.iter().filter(|l| !l.starts_with('#')).count();
        Ok(ClinicalVcf {
            records,
            high: high_impact(&lines),
        })
    }
}

/// `csq_high_impact_genes`.
fn high_impact(lines: &[String]) -> Result<IndexMap<Key, HighRecord>, String> {
    let floor = ClinicalFloor::read(lines)?;
    let mut fmt: Option<Vec<String>> = None;
    let mut out = IndexMap::new();
    for line in lines {
        if line.starts_with("##INFO=<ID=CSQ") {
            fmt = line.split_once("Format: ").and_then(|(_, f)| {
                let f = f.split('"').next().unwrap_or("");
                (!f.is_empty()).then(|| f.split('|').map(str::to_owned).collect())
            });
            continue;
        }
        let Some(fm) = fmt.as_ref() else { continue };
        if line.starts_with('#')
            || !fm.iter().any(|c| c == "IMPACT")
            || !fm.iter().any(|c| c == "SYMBOL")
        {
            continue;
        }
        let r: Vec<String> = line.split('\t').map(str::to_owned).collect();
        if r.len() < 8 {
            continue;
        }
        let (cls, entries, freq) = csq_rank(fm, info_map(&r[7]).get("CSQ").copied().unwrap_or(""));
        let mut genes: IndexMap<String, String> = IndexMap::new();
        for e in entries.into_iter().filter(|e| e.impact == "HIGH") {
            genes.entry(e.gene).or_insert(e.consequence);
        }
        if genes.is_empty() {
            continue;
        }
        let gt = text::gt_of(&r).to_owned();
        let (weak, why, gq, dp, ad) = if r.len() > 9 {
            call_check(&r[8], &r[9], &floor)
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
    }
    Ok(out)
}

/// A source of a candidate variant.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
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
}

/// `sf_verdicts`: for each distinct variant of a gene, `None` (reported) or
/// the reason it is not, and the partial-deficiency note of a reported gene
/// with a hypomorph allele.
fn sf_verdicts(
    gene: &str,
    rule: &Rule,
    alleles: &IndexMap<AKey, (usize, usize)>,
) -> (HashMap<AKey, Option<String>>, Option<String>) {
    let sev = |k: &AKey| match k {
        AKey::Small(k) => rule.severity.get(k),
        _ => None,
    };
    let mut out: IndexMap<AKey, Option<String>> = IndexMap::new();
    let mut kept: IndexMap<AKey, (usize, usize)> = IndexMap::new();
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
                kept.insert(k.clone(), *cp);
            }
        }
    }
    let n: usize = kept.values().map(|(c, _)| c).sum();
    let ri = rule.report_if.as_str();
    let hemi = kept.values().any(|(_, p)| *p == 1);
    for (k, (c, p)) in &kept {
        let good = match ri {
            "any" | "het_female_too" => true,
            "biallelic" => n >= 2,
            "hemizygous" => hemi || n >= 2,
            _ => matches!(k, AKey::Small(s) if rule.variants.contains_key(s)) && *c >= 2 && *p == 2,
        };
        let why = if good {
            None
        } else if ri == "hom_variant" {
            Some(format!("{gene} is reported only for: {}", rule.table_text))
        } else if ri == "hemizygous" {
            Some(format!(
                "one heterozygous variant; {gene} is reported for {}",
                rule.table_text
            ))
        } else {
            Some(format!(
                "one heterozygous variant; {gene} is reported with {}",
                rule.table_text
            ))
        };
        out.insert(k.clone(), why);
    }
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
    (out.into_iter().collect(), partial)
}

/// `sf_structural`: step 05's class 4 and 5 SVs in genes whose rule allows
/// them, and the class 4 or 5 deletions over a gene's sv_regions.
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
            });
        }
    }
    (out, true, None)
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
                clinvar.push(Cand {
                    kind: Kind::ClinVar,
                    gene: g.to_owned(),
                    fields: f,
                    key: AKey::Small(row.key()),
                    copies,
                    ploidy,
                    gt: gt.clone(),
                    ps: gt_field(r, "PS"),
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
                });
            }
        }
    }
    let (svs, sv_checked, sv_note) = sf_structural(i.deletions, &rules);

    // One verdict per gene, over its distinct variants.
    let mut alleles: IndexMap<String, IndexMap<AKey, (usize, usize)>> = IndexMap::new();
    for x in clinvar.iter().chain(&high) {
        let p = if x.ploidy == 0 { 2 } else { x.ploidy };
        alleles
            .entry(x.gene.clone())
            .or_default()
            .entry(x.key.clone())
            .or_insert((x.copies.max(1), p));
    }
    for x in &svs {
        alleles
            .entry(x.gene.clone())
            .or_default()
            .entry(x.key.clone())
            .or_insert((1, 2));
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
    let mut seen = HashSet::new();
    let mut variants = Vec::new();
    for x in clinvar.iter().chain(&high) {
        let AKey::Small(k) = &x.key else { continue };
        let pseudo = rules[&x.gene].severity_kind(k) == Some("pseudodeficiency");
        if !pseudo && seen.insert((x.gene.clone(), k.clone())) {
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
                gt: x.gt.clone(),
                ps: x.ps.clone(),
                allele: None,
            });
        }
    }
    let found = two_hits(&variants, &dels);
    let mut two: HashMap<(String, Key), String> = HashMap::new();
    for g in found.entries.keys() {
        if matches!(rules[g].report_if.as_str(), "biallelic" | "hemizygous") {
            if let Some(a) = alleles.get_mut(g) {
                a.insert(AKey::Deletion(g.clone()), (1, 2));
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
        .map(|(g, a)| (g, sf_verdicts(g, &rules[g], a)))
        .collect();
    let mut reported: [Vec<Value>; 3] = Default::default();
    let mut not_reported = Vec::new();
    for x in clinvar.into_iter().chain(high).chain(svs) {
        let (whys, partial) = &verdict[&x.gene];
        let why = whys.get(&x.key).cloned().flatten();
        let rule = &rules[&x.gene];
        let mut f = x.fields;
        if let AKey::Small(k) = &x.key
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
                f.insert("severity".into(), opt_s(partial.as_deref()));
                if x.kind != Kind::Sv
                    && let AKey::Small(k) = &x.key
                    && let Some(t) = two.get(&(x.gene.clone(), k.clone()))
                {
                    f.insert("two_hit".into(), Value::Str(t.clone()));
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
        .filter(|x| !(flag(x, "weak_call") || flag(x, "frequency_unknown")))
        .count();
    let weak = high_rows.iter().filter(|x| flag(x, "weak_call")).count();
    let unknown = high_rows
        .iter()
        .filter(|x| flag(x, "frequency_unknown"))
        .count();
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
        ("high_impact_frequency_unknown", Value::Int(unknown as i128)),
        ("structural", Value::Array(sv_rows)),
        ("structural_checked", Value::Bool(sv_checked)),
        ("structural_note", opt_s(sv_note.as_deref())),
        ("not_reported", Value::Array(not_reported)),
        ("two_hit_note", opt_s(two_hit_note.as_deref())),
    ])))
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
