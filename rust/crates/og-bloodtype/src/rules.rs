//! The rules: from the sites to one row per blood group system.
//!
//! Phase is not read: two variants of one system are placed on the
//! chromosomes the common alleles carry them on, and the note says so.
//! Confidence: high (every site called or covered by reads), medium (a site
//! read as reference with no depth to confirm it, or a null row of the
//! system not determined), low (a site in a region where short reads are
//! unreliable, named in the note, or markers that disagree), not determined.

use crate::Error;
use crate::genotypes::{Genotyped, Site, Status};
use crate::markers::Regions;
use og_cn::catalog::{CallRule, Locus};
use og_cn::engine::{self, AssessError, Means, Reason, Sex};
use pycompat::fmt::fixed;
use std::collections::HashMap;
use std::sync::LazyLock;

/// The exons an RHD-CE(3-7)-D hybrid takes from RHCE.
pub const RHD_HYBRID_EXONS: [i64; 5] = [3, 4, 5, 6, 7];
pub const NOT_DETERMINED: &str = "not determined";

/// A two-allele site: its label and its phenotype template, where
/// `{antigen}` is replaced by `+` or `-`. The antigen names are the marker
/// table's `ref_means` and `alt_means`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Pair {
    pub label: String,
    pub template: String,
}

/// HPA-7 to -35 (HPA-15 has its own line): one site each, the a antigen on
/// the reference allele and the low-prevalence b antigen on the alternative.
pub const HPA_RARE: [&str; 28] = [
    "HPA7_1297",
    "HPA8_1984",
    "HPA9_2602",
    "HPA10_263",
    "HPA11_1976",
    "HPA12_119",
    "HPA13_2483",
    "HPA14_1909",
    "HPA16_497",
    "HPA17_662",
    "HPA18_2235",
    "HPA19_487",
    "HPA20_1949",
    "HPA21_1960",
    "HPA22_584",
    "HPA23_1942",
    "HPA24_1508",
    "HPA25_3347",
    "HPA26_1818",
    "HPA27_2614",
    "HPA28_2311",
    "HPA29_98",
    "HPA30_2511",
    "HPA31_368",
    "HPA32_521",
    "HPA33_1373",
    "HPA34_349",
    "HPA35_1514",
];

/// The number of an `HPA<k>_<pos>` name.
fn hpa_number(name: &str) -> &str {
    name[3..].split('_').next().unwrap_or_default()
}

/// Every two-allele site, by marker name, in the order the rules list them.
pub static PAIRS: LazyLock<Vec<(String, Pair)>> = LazyLock::new(|| {
    let fixed: [(&str, &str, &str); 22] = [
        ("RHCE_307", "C/c", "C{C}c{c}"),
        ("RHCE_676", "E/e", "E{E}e{e}"),
        ("KEL_578", "K/k", "K{K}k{k}"),
        ("KEL_841", "Kp", "Kp(a{Kpa}b{Kpb})"),
        ("KEL_1790", "Js", "Js(a{Jsa}b{Jsb})"),
        ("JK_838", "Jk", "Jk(a{Jka}b{Jkb})"),
        ("GYPA_59", "M/N", "M{M}N{N}"),
        ("GYPB_143", "S/s", "S{S}s{s}"),
        ("LU_230", "Lu", "Lu(a{Lua}b{Lub})"),
        ("DI_2561", "Di", "Di(a{Dia}b{Dib})"),
        ("DO_793", "Do", "Do(a{Doa}b{Dob})"),
        ("CO_134", "Co", "Co(a{Coa}b{Cob})"),
        ("HPA1_176", "HPA-1", "HPA-1(a{1a}b{1b})"),
        ("HPA2_482", "HPA-2", "HPA-2(a{2a}b{2b})"),
        ("HPA3_2621", "HPA-3", "HPA-3(a{3a}b{3b})"),
        ("HPA4_506", "HPA-4", "HPA-4(a{4a}b{4b})"),
        ("HPA5_1600", "HPA-5", "HPA-5(a{5a}b{5b})"),
        ("HPA6_1544", "HPA-6", "HPA-6(a{6a}b{6b})"),
        ("HPA15_2108", "HPA-15", "HPA-15(a{15a}b{15b})"),
        ("HNA3_461", "HNA-3", "HNA-3(a{3a}b{3b})"),
        ("HNA4_230", "HNA-4", "HNA-4(a{4a}b{4b})"),
        ("HNA5_2372", "HNA-5", "HNA-5(a{5a}b{5b})"),
    ];
    let mut out: Vec<(String, Pair)> = fixed
        .iter()
        .map(|(n, l, t)| {
            (
                (*n).to_owned(),
                Pair {
                    label: (*l).to_owned(),
                    template: (*t).to_owned(),
                },
            )
        })
        .collect();
    for n in HPA_RARE {
        let k = hpa_number(n);
        out.push((
            n.to_owned(),
            Pair {
                label: format!("HPA-{k}"),
                template: format!("HPA-{k}(a{{{k}a}}b{{{k}b}})"),
            },
        ));
    }
    out
});

/// The two-allele site called `name`.
pub fn pair(name: &str) -> Option<&'static Pair> {
    PAIRS.iter().find(|(n, _)| n == name).map(|(_, p)| p)
}

/// A system read from its null rows alone.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct NullSystem {
    pub system: &'static str,
    /// The phenotype with the antigen.
    pub plus: &'static str,
    /// The phenotype without it.
    pub minus: &'static str,
    /// The allele label of a working copy.
    pub ok: &'static str,
    /// The allele label of a null copy.
    pub null: &'static str,
}

pub const NULL_SYSTEMS: [NullSystem; 5] = [
    NullSystem {
        system: "Vel",
        plus: "Vel+",
        minus: "Vel-",
        ok: "Vel",
        null: "null",
    },
    NullSystem {
        system: "Langereis",
        plus: "Lan+",
        minus: "Lan-",
        ok: "Lan",
        null: "null",
    },
    NullSystem {
        system: "Junior",
        plus: "Jr(a+)",
        minus: "Jr(a-)",
        ok: "Jra",
        null: "null",
    },
    NullSystem {
        system: "H",
        plus: "H+",
        minus: "H-",
        ok: "H",
        null: "h",
    },
    NullSystem {
        system: "Secretor",
        plus: "secretor",
        minus: "non-secretor",
        ok: "Se",
        null: "se",
    },
];

fn null_system(system: &str) -> Option<&'static NullSystem> {
    NULL_SYSTEMS.iter().find(|s| s.system == system)
}

/// Systems with two-allele sites whose null rows silence a chromosome.
pub const NULL_HOSTS: [&str; 2] = ["Kell", "Kidd"];

/// HNA-1 (FCGR3B): FCGR3B*01 (HNA-1a) and *02 (HNA-1b) differ at c.108,
/// c.114, c.194, c.244 and c.316; *03 (HNA-1c) is *02 with c.233C>A and
/// carries HNA-1b as well. GRCh38 carries the *02 base at four of the five
/// and the *01 base at c.194, so the call counts HNA-1a on the four and only
/// checks c.194.
pub const HNA1_SITES: [&str; 4] = ["FCGR3B_108", "FCGR3B_114", "FCGR3B_244", "FCGR3B_316"];

/// The sites the special rules (ABO, RhD, Duffy, HNA-1) read by name.
pub const SPECIAL_SITES: [&str; 16] = [
    "ABO_261",
    "ABO_802",
    "ABO_803",
    "ABO_796",
    "ABO_1061",
    "RHD_PSI_DUP",
    "RHD_PSI_609",
    "FY_125",
    "FY_GATA",
    "FY_265",
    "FCGR3B_108",
    "FCGR3B_114",
    "FCGR3B_244",
    "FCGR3B_316",
    "FCGR3B_194",
    "FCGR3B_233",
];

/// A system and its two-allele sites (`None` for ABO, RhD, Duffy, HNA-1
/// and the systems read from null rows alone).
pub type System = (String, Option<Vec<&'static str>>);

/// The systems in their base order (the report then lists first the ones
/// that lack a high-prevalence antigen).
pub static SYSTEMS: LazyLock<Vec<System>> = LazyLock::new(|| {
    let s = |name: &str, sites: &[&'static str]| (name.to_owned(), Some(sites.to_vec()));
    let none = |name: &str| (name.to_owned(), None);
    let mut out = vec![
        none("ABO"),
        none("RhD"),
        s("RhCE", &["RHCE_307", "RHCE_676"]),
        s("Kell", &["KEL_578", "KEL_841", "KEL_1790"]),
        none("Duffy"),
        s("Kidd", &["JK_838"]),
        s("MNS", &["GYPA_59", "GYPB_143"]),
        s("Lutheran", &["LU_230"]),
        s("Diego", &["DI_2561"]),
        s("Dombrock", &["DO_793"]),
        s("Colton", &["CO_134"]),
        none("Vel"),
        none("Langereis"),
        none("Junior"),
        none("H"),
        none("Secretor"),
        s("HPA-1", &["HPA1_176"]),
        s("HPA-2", &["HPA2_482"]),
        s("HPA-3", &["HPA3_2621"]),
        s("HPA-4", &["HPA4_506"]),
        s("HPA-5", &["HPA5_1600"]),
        s("HPA-6", &["HPA6_1544"]),
    ];
    let rare = |below: bool| {
        HPA_RARE
            .iter()
            .filter(move |n| {
                let k: u32 = hpa_number(n).parse().unwrap_or(0);
                if below { k < 15 } else { k > 15 }
            })
            .map(|n| (format!("HPA-{}", hpa_number(n)), Some(vec![*n])))
    };
    out.extend(rare(true));
    out.push(s("HPA-15", &["HPA15_2108"]));
    out.extend(rare(false));
    out.extend([
        none("HNA-1"),
        s("HNA-3", &["HNA3_461"]),
        s("HNA-4", &["HNA4_230"]),
        s("HNA-5", &["HNA5_2372"]),
    ]);
    out
});

/// One row of the table.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Row {
    pub system: String,
    pub phenotype: String,
    pub genotype: String,
    pub confidence: &'static str,
    pub note: String,
    /// The high-prevalence antigens the person lacks: such a row is listed first.
    pub lacks_high: Vec<String>,
}

impl Row {
    fn new(
        system: &str,
        phenotype: impl Into<String>,
        genotype: impl Into<String>,
        confidence: &'static str,
        note: impl Into<String>,
    ) -> Row {
        Row {
            system: system.to_owned(),
            phenotype: phenotype.into(),
            genotype: genotype.into(),
            confidence,
            note: note.into(),
            lacks_high: Vec::new(),
        }
    }

    fn not_determined(system: &str, note: impl Into<String>) -> Row {
        Row::new(system, NOT_DETERMINED, ".", NOT_DETERMINED, note)
    }

    fn add_note(&mut self, part: &str) {
        self.note = join_parts([self.note.as_str(), part]);
    }
}

/// `"; ".join(x for x in parts if x)`.
fn join_parts<'a>(parts: impl IntoIterator<Item = &'a str>) -> String {
    parts
        .into_iter()
        .filter(|p| !p.is_empty())
        .collect::<Vec<_>>()
        .join("; ")
}

/// `[x] * n`, empty for n <= 0.
fn rep(x: &str, n: i64) -> impl Iterator<Item = &str> {
    std::iter::repeat_n(x, usize::try_from(n).unwrap_or(0))
}

fn sign(present: bool) -> &'static str {
    if present { "+" } else { "-" }
}

fn status_part(s: &Site) -> String {
    let base = format!("{}: {}", s.row.hgvs, s.status.as_str());
    if s.note.is_empty() {
        base
    } else {
        format!("{base} ({})", s.note)
    }
}

fn is_stray(s: &Site) -> bool {
    s.note.starts_with("another allele")
}

/// Confidence from the sites a call used: low for a caveat or a stray
/// allele, medium for a site read as reference without depth, else high.
/// Another ALT at a null row's position is not that null allele and leaves
/// the antigens read, so it is only named in the note.
fn level_of(used: &[&Site]) -> &'static str {
    let mut lv = "high";
    for s in used {
        if !s.row.caveat.is_empty() || (is_stray(s) && !s.row.is_null()) {
            return "low";
        }
        if s.status == Status::Assumed {
            lv = "medium";
        }
    }
    lv
}

fn notes_of(used: &[&Site], extra: Vec<String>) -> String {
    let mut out = extra;
    for s in used {
        if !s.row.caveat.is_empty() && !out.contains(&s.row.caveat) {
            out.push(s.row.caveat.clone());
        }
        if is_stray(s) {
            out.push(format!("{}: {}", s.row.hgvs, s.note));
        }
    }
    if used.iter().any(|s| s.status == Status::Assumed) {
        out.push(
            "sites with no VCF record read as reference (no BAM depth to confirm them)".to_owned(),
        );
    }
    out.join("; ")
}

fn unknown(system: &str, missing: &[&Site], extra: &str) -> Row {
    let why = missing
        .iter()
        .map(|s| status_part(s))
        .collect::<Vec<_>>()
        .join("; ");
    Row::not_determined(system, join_parts([extra, why.as_str()]))
}

/// The copies left of each antigen of a two-allele site, alt first.
type Left = Vec<(String, i64)>;

fn count(left: &Left, antigen: &str) -> i64 {
    left.iter()
        .rev()
        .find(|(a, _)| a == antigen)
        .map_or(0, |(_, c)| *c)
}

fn site_left(s: &Site, n: i64) -> Left {
    vec![
        (s.row.alt_means.clone(), n),
        (s.row.ref_means.clone(), 2 - n),
    ]
}

/// (phenotype part, genotype part) of one two-allele site, or `None`.
/// `left` holds the antigens' copies after null copies took theirs (default:
/// the site's own).
fn pair_part(name: &str, s: &Site, left: Option<&Left>) -> Option<(String, String)> {
    let own;
    let left = match (left, s.n) {
        (Some(l), _) => l,
        (None, Some(n)) => {
            own = site_left(s, n);
            &own
        }
        (None, None) => return None,
    };
    let (alt, rf) = (&s.row.alt_means, &s.row.ref_means);
    let (na, nr) = (count(left, alt), count(left, rf));
    let tmpl = &pair(name)
        .expect("a two-allele site has a template")
        .template;
    let pheno = tmpl
        .replace(&format!("{{{alt}}}"), sign(na > 0))
        .replace(&format!("{{{rf}}}"), sign(nr > 0));
    let mut alleles: Vec<&str> = rep(alt, na).chain(rep(rf, nr)).collect();
    alleles.sort_unstable();
    let pad = 2 - alleles.len() as i64;
    alleles.extend(rep("null", pad));
    Some((pheno, alleles.join("/")))
}

fn null_rows<'a>(system: &str, sites: &'a Genotyped) -> Vec<&'a Site> {
    sites
        .values()
        .filter(|s| s.row.system == system && s.row.is_null())
        .collect()
}

/// A system's null rows read together.
struct NullCopies<'a> {
    /// Null copies, at most two.
    n: i64,
    /// The null sites with copies.
    found: Vec<&'a Site>,
    /// The null sites not determined.
    missing: Vec<&'a Site>,
    extra: Vec<String>,
}

fn null_copies<'a>(nulls: &[&'a Site]) -> NullCopies<'a> {
    let missing: Vec<&Site> = nulls.iter().copied().filter(|s| s.n.is_none()).collect();
    let found: Vec<&Site> = nulls
        .iter()
        .copied()
        .filter(|s| s.n.is_some_and(|n| n != 0))
        .collect();
    let n: i64 = found.iter().filter_map(|s| s.n).sum();
    let mut extra = Vec::new();
    if !found.is_empty() {
        let list: Vec<String> = found
            .iter()
            .map(|s| format!("{} x{}", s.row.hgvs, s.n.unwrap_or(0)))
            .collect();
        extra.push(format!(
            "null allele{}: {}",
            if found.len() > 1 { "s" } else { "" },
            list.join(", ")
        ));
    }
    if n >= 2 && !found.iter().any(|s| s.n == Some(2)) {
        extra.push(
            "two different null alleles, assumed to sit on different chromosomes (phase is not read)"
                .to_owned(),
        );
    }
    if n > 2 {
        extra.push(format!(
            "{n} null copies, more than two chromosomes carry: a rare allele or a misplaced read"
        ));
    }
    if !missing.is_empty() {
        let first: Vec<String> = missing
            .iter()
            .take(3)
            .map(|s| format!("{}: {}", s.row.hgvs, s.status.as_str()))
            .collect();
        extra.push(format!(
            "{} null row{} not determined ({}{})",
            missing.len(),
            if missing.len() > 1 { "s" } else { "" },
            first.join("; "),
            if missing.len() > 3 { "; ..." } else { "" }
        ));
    }
    NullCopies {
        n: n.min(2),
        found,
        missing,
        extra,
    }
}

/// `level_of` the sites used, then medium when a null row was not
/// determined, and low for two different null alleles (no phase) or more
/// than two copies.
fn null_confidence(used: &[&Site], found: &[&Site], missing: &[&Site]) -> &'static str {
    let conf = level_of(used);
    let copies: i64 = found.iter().filter_map(|s| s.n).sum();
    if copies > 2 || (copies == 2 && found.len() > 1) {
        return "low";
    }
    if conf == "high" && !missing.is_empty() {
        return "medium";
    }
    conf
}

fn hp_prefix(system: &str) -> String {
    if system.starts_with("HPA-") || system.starts_with("HNA-") {
        format!("{}-", system.split('-').next().unwrap_or_default())
    } else {
        String::new()
    }
}

/// Marks the row with the high-prevalence antigens it lacks (the table's
/// `high_prevalence` column), which [`predict`] lists first.
fn with_lacks(mut row: Row, lacks: &[String], used: &[&Site]) -> Row {
    let hp: Vec<&str> = used
        .iter()
        .map(|s| s.row.high_prevalence.as_str())
        .filter(|h| !h.is_empty())
        .collect();
    let mut gone: Vec<String> = Vec::new();
    for a in lacks {
        if !gone.contains(a) && hp.contains(&a.as_str()) {
            gone.push(a.clone());
        }
    }
    if !gone.is_empty() {
        let prefix = hp_prefix(&row.system);
        let names: Vec<String> = gone.iter().map(|a| format!("{prefix}{a}")).collect();
        row.add_note(&format!(
            "lacks {}, carried by most people (high prevalence)",
            names.join(", ")
        ));
        row.lacks_high = gone;
    }
    row
}

/// The table's meaning column for the antigens the person lacks, as a note
/// part: `lacks HPA-1a: <meaning>`. Only rows with a meaning (HPA, HNA) give one.
fn meaning(system: &str, lacks: &[String], used: &[&Site]) -> Vec<String> {
    let text = used
        .iter()
        .map(|s| s.row.meaning.as_str())
        .find(|m| !m.is_empty())
        .unwrap_or_default();
    if text.is_empty() || lacks.is_empty() {
        return Vec::new();
    }
    let prefix = format!("{}-", system.split('-').next().unwrap_or_default());
    let names: Vec<String> = lacks.iter().map(|a| format!("{prefix}{a}")).collect();
    vec![format!("lacks {}: {text}", names.join(", "))]
}

fn system_pairs(system: &str, names: &[&str], sites: &Genotyped) -> Row {
    let (mut parts, mut genos) = (Vec::new(), Vec::new());
    let (mut used, mut missing): (Vec<&Site>, Vec<&Site>) = (Vec::new(), Vec::new());
    let nulls = null_rows(system, sites);
    let nc = null_copies(&nulls);
    // What each null copy silences: one allele at every site, the one of its
    // null_on column first (the chromosome the ISBT table finds it on), else
    // the other.
    let mut left: Vec<(&str, Left)> = Vec::new();
    for &name in names {
        if let Some(n) = sites[name].n {
            left.push((name, site_left(&sites[name], n)));
        }
    }
    let queue: Vec<&Site> = nc
        .found
        .iter()
        .flat_map(|x| rep("", x.n.unwrap_or(0)).map(move |_| *x))
        .take(usize::try_from(nc.n).unwrap_or(0))
        .collect();
    for x in queue {
        let prefer: Vec<&str> = x.row.null_on.split_whitespace().collect();
        for (name, cnt) in &mut left {
            let r = &sites[*name].row;
            let pick = prefer
                .iter()
                .find(|a| count(cnt, a) > 0)
                .map(|a| (*a).to_owned())
                .or_else(|| {
                    [&r.ref_means, &r.alt_means]
                        .into_iter()
                        .find(|a| count(cnt, a) > 0)
                        .cloned()
                })
                .expect("a site holds two copies, more than the null copies take");
            if let Some(c) = cnt.iter_mut().rev().find(|(a, _)| *a == pick) {
                c.1 -= 1;
            }
        }
    }
    // Two null copies silence both chromosomes: a pair site with no call
    // still lacks both antigens.
    if nc.n >= 2 {
        for &name in names {
            if !left.iter().any(|(n, _)| *n == name) {
                let r = &sites[name].row;
                left.push((
                    name,
                    vec![(r.alt_means.clone(), 0), (r.ref_means.clone(), 0)],
                ));
            }
        }
    }
    let left_of = |name: &str| left.iter().find(|(n, _)| *n == name).map(|(_, l)| l);
    for &name in names {
        let s = &sites[name];
        match pair_part(name, s, left_of(name)) {
            None => {
                let label = &pair(name).expect("a two-allele site").label;
                parts.push(format!("{label} {NOT_DETERMINED}"));
                missing.push(s);
            }
            Some((p, g)) => {
                parts.push(p);
                genos.push(g);
                used.push(s);
            }
        }
    }
    if used.is_empty() {
        return unknown(system, &missing, "");
    }
    let mut used_all = used.clone();
    used_all.extend(nulls.iter().copied().filter(|s| s.n.is_some()));
    let conf = if missing.is_empty() {
        null_confidence(&used_all, &nc.found, &nc.missing)
    } else {
        "low"
    };
    let mut extra: Vec<String> = missing.iter().map(|s| status_part(s)).collect();
    extra.extend(nc.extra.iter().cloned());
    let placed: Vec<&Site> = nc
        .found
        .iter()
        .copied()
        .filter(|s| !s.row.null_on.is_empty())
        .collect();
    let on = |s: &Site| {
        s.row
            .null_on
            .split_whitespace()
            .collect::<Vec<_>>()
            .join(", ")
    };
    if placed.len() == 1 {
        extra.push(format!(
            "each null copy is placed with {} first, the alleles the ISBT table finds it with (phase is not read)",
            on(placed[0])
        ));
    } else if !placed.is_empty() {
        let each: Vec<String> = placed
            .iter()
            .map(|s| format!("{}: {}", s.row.hgvs, on(s)))
            .collect();
        extra.push(format!(
            "each null copy is placed with the alleles the ISBT table finds it with first ({}; phase is not read)",
            each.join("; ")
        ));
    }
    let mut lacks: Vec<String> = Vec::new();
    for &name in names {
        if let Some(l) = left_of(name) {
            lacks.extend(l.iter().filter(|(_, c)| *c == 0).map(|(a, _)| a.clone()));
        }
    }
    if nc.n >= 2 {
        let antigen = &nc.found[0].row.ref_means;
        lacks.push(antigen.clone());
        extra.insert(
            0,
            format!("two null copies: no {system} antigen on the red cells ({antigen}-)"),
        );
    }
    let geno = if nc.n >= 2 && missing.is_empty() {
        "null/null".to_owned()
    } else {
        genos.join(" ")
    };
    extra.extend(meaning(system, &lacks, &used));
    let row = Row::new(
        system,
        parts.join(" "),
        geno,
        conf,
        notes_of(&used_all, extra),
    );
    with_lacks(row, &lacks, &used_all)
}

/// A system read from its null rows alone (Vel, Langereis, Junior, H, secretor).
fn system_null(ns: &NullSystem, sites: &Genotyped) -> Row {
    let nulls = null_rows(ns.system, sites);
    let nc = null_copies(&nulls);
    let used: Vec<&Site> = nulls.iter().copied().filter(|s| s.n.is_some()).collect();
    if used.is_empty() {
        let first: Vec<&Site> = nc.missing.iter().copied().take(3).collect();
        return unknown(
            ns.system,
            &first,
            &format!("no null row of {} was determined", ns.system),
        );
    }
    let conf = null_confidence(&used, &nc.found, &nc.missing);
    let lacks = if nc.n >= 2 {
        vec![nulls[0].row.ref_means.clone()]
    } else {
        Vec::new()
    };
    let geno: Vec<&str> = rep(ns.ok, 2 - nc.n).chain(rep(ns.null, nc.n)).collect();
    let row = Row::new(
        ns.system,
        if nc.n >= 2 { ns.minus } else { ns.plus },
        geno.join("/"),
        conf,
        notes_of(&used, nc.extra.clone()),
    );
    with_lacks(row, &lacks, &used)
}

fn site<'a>(sites: &'a Genotyped, name: &str) -> &'a Site {
    &sites[name]
}

fn system_hna1(sites: &Genotyped) -> Row {
    let mut used: Vec<&Site> = HNA1_SITES.iter().map(|n| site(sites, n)).collect();
    let missing: Vec<&Site> = used.iter().copied().filter(|s| s.n.is_none()).collect();
    if !missing.is_empty() {
        return unknown("HNA-1", &missing, "");
    }
    let counts: Vec<i64> = used.iter().filter_map(|s| s.n).collect();
    if counts.iter().any(|c| *c != counts[0]) {
        let list: Vec<String> = used
            .iter()
            .map(|s| format!("{} x{}", s.row.hgvs, s.n.unwrap_or(0)))
            .collect();
        return unknown(
            "HNA-1",
            &[],
            &format!(
                "the FCGR3B markers do not agree on HNA-1a or HNA-1b ({}): an FCGR3A read placed on FCGR3B, an extra FCGR3B copy, or a rare allele; {}",
                list.join(", "),
                used[0].row.caveat
            ),
        );
    }
    let n_a = counts[0];
    let n_b = 2 - n_a;
    let (s233, s194) = (site(sites, "FCGR3B_233"), site(sites, "FCGR3B_194"));
    let mut extra = Vec::new();
    let n_c = s233.n;
    match n_c {
        None => extra.push(format!(
            "HNA-1c not determined ({}: {})",
            s233.row.hgvs,
            s233.status.as_str()
        )),
        Some(_) => used.push(s233),
    }
    if let Some(n194) = s194.n
        && 2 - n194 != n_a
    {
        extra.push(format!(
            "{} (x{n194}) does not match the other FCGR3B markers",
            s194.row.hgvs
        ));
    }
    let has_c = n_c.is_some_and(|c| c != 0);
    if has_c && n_c.unwrap_or(0) > n_b {
        extra.push(
            "HNA-1c with no HNA-1b allele left for it: often an extra FCGR3B copy carrying FCGR3B*03"
                .to_owned(),
        );
    }
    let n_c2 = n_c.unwrap_or(0).min(n_b);
    // Without c.233 an HNA-1b allele may be *02 or *03 (*03 is *02 with c.233C>A).
    let b_allele = if n_c.is_some() {
        "FCGR3B*02"
    } else {
        "FCGR3B*02 or *03"
    };
    let alleles: Vec<&str> = rep("FCGR3B*01", n_a)
        .chain(rep(b_allele, n_b - n_c2))
        .chain(rep("FCGR3B*03", n_c2))
        .collect();
    let c_txt = n_c.map_or(String::new(), |c| format!("c{}", sign(c > 0)));
    let has_b = n_b > 0 || has_c;
    let lacks: Vec<String> = [("1a", n_a > 0), ("1b", has_b), ("1c", has_c)]
        .into_iter()
        .filter(|(x, has)| !(*has || (*x == "1c" && n_c.is_none())))
        .map(|(x, _)| x.to_owned())
        .collect();
    extra.extend(meaning("HNA-1", &lacks, &used));
    Row::new(
        "HNA-1",
        format!("HNA-1(a{}b{}{c_txt})", sign(n_a > 0), sign(has_b)),
        alleles.join("/"),
        level_of(&used),
        notes_of(&used, extra),
    )
}

fn system_abo(sites: &Genotyped) -> Row {
    let need: Vec<&Site> = ["ABO_261", "ABO_802", "ABO_803"]
        .iter()
        .map(|n| site(sites, n))
        .collect();
    let missing: Vec<&Site> = need.iter().copied().filter(|s| s.n.is_none()).collect();
    if !missing.is_empty() {
        return unknown("ABO", &missing, "");
    }
    let [n_func, n_o2, n_b] = [0, 1, 2].map(|i| need[i].n.unwrap_or(0));
    let non_o = n_func - n_o2;
    if non_o < 0 || n_b > non_o {
        return unknown(
            "ABO",
            &[],
            &format!(
                "the markers do not add up to two ABO alleles (c.261 insertion x{n_func}, c.802G>A x{n_o2}, c.803G>C x{n_b})"
            ),
        );
    }
    let n_a = non_o - n_b;
    let (mut used, mut extra) = (need.clone(), Vec::new());
    let mut alleles: Vec<&str> = Vec::new();
    let mut a_label = "";
    if n_a != 0 {
        let s1061 = site(sites, "ABO_1061");
        match s1061.n {
            None => {
                alleles.extend(rep("A", n_a));
                a_label = "A";
                extra.push(format!(
                    "A1 or A2 not determined ({}: {})",
                    s1061.row.hgvs,
                    s1061.status.as_str()
                ));
            }
            Some(n1061) => {
                // A1 or A2 only: c.1061delC's caveat goes to the note, not to
                // the confidence of the group.
                if !s1061.row.caveat.is_empty() {
                    extra.push(format!(
                        "A1 or A2 rests on {}, {}",
                        s1061.row.hgvs, s1061.row.caveat
                    ));
                }
                if s1061.status == Status::Assumed {
                    used.push(s1061);
                }
                let n_a2 = n1061.min(n_a);
                alleles.extend(rep("A1", n_a - n_a2).chain(rep("A2", n_a2)));
                a_label = if n_a > n_a2 { "A1" } else { "A2" };
            }
        }
    }
    alleles.extend(
        rep("B", n_b)
            .chain(rep("O.01", 2 - n_func))
            .chain(rep("O.02", n_o2)),
    );
    let s796 = site(sites, "ABO_796");
    let mut low = false;
    if let Some(n796) = s796.n
        && n796 != n_b
    {
        low = true;
        extra.push(format!(
            "c.796C>A (x{n796}) and c.803G>C (x{n_b}) disagree: a rare A or B allele"
        ));
    }
    let mut pheno = format!("{a_label}{}", if n_b != 0 { "B" } else { "" });
    if pheno.is_empty() {
        pheno = "O".to_owned();
    }
    let conf = if low { "low" } else { level_of(&used) };
    if n_b != 0 && n_a != 0 {
        extra.push("A and B on different chromosomes assumed (no phase)".to_owned());
    }
    Row::new(
        "ABO",
        pheno,
        alleles.join("/"),
        conf,
        notes_of(&used, extra),
    )
}

fn system_duffy(sites: &Genotyped) -> Row {
    let (s125, sg, s265) = (
        site(sites, "FY_125"),
        site(sites, "FY_GATA"),
        site(sites, "FY_265"),
    );
    let missing: Vec<&Site> = [s125, sg].into_iter().filter(|s| s.n.is_none()).collect();
    if !missing.is_empty() {
        return unknown("Duffy", &missing, "");
    }
    let (nb, ng) = (s125.n.unwrap_or(0), sg.n.unwrap_or(0));
    let sil_b = ng.min(nb);
    let sil_a = ng - sil_b;
    let (a_expr, b_expr) = ((2 - nb) - sil_a, nb - sil_b);
    let (mut used, mut extra) = (vec![s125, sg], Vec::new());
    let mut weak = 0;
    if let Some(n265) = s265.n {
        used.push(s265);
        weak = n265.min(b_expr);
    }
    let b_sign = if b_expr - weak > 0 {
        "+"
    } else if weak != 0 {
        "+w"
    } else {
        "-"
    };
    let alleles: Vec<&str> = rep("FY*A", 2 - nb - sil_a)
        .chain(rep("FY*A GATA", sil_a))
        .chain(rep("FY*B", b_expr - weak))
        .chain(rep("FY*X", weak))
        .chain(rep("FY*B GATA", sil_b))
        .collect();
    if ng != 0 {
        extra.push(
            "c.-67T>C (GATA box) stops Duffy expression on red cells only; with no phase it is placed on the Fyb chromosome first, where it usually lies (FY*02N.01)"
                .to_owned(),
        );
    }
    Row::new(
        "Duffy",
        format!("Fy(a{}b{b_sign})", sign(a_expr > 0)),
        alleles.join("/"),
        level_of(&used),
        notes_of(&used, extra),
    )
}

/// RhD's depth reading: the RHD row of the catalog, assessed once.
#[derive(Debug, Clone, PartialEq)]
pub enum RhdDepth {
    /// No BAM was given.
    NoBam,
    /// A control segment has no read while the others are covered.
    Broken(String),
    /// No depth over the row's segments (the BAM has no such contig).
    NoData,
    /// Too little depth beside RHD: the flank depth.
    LowFlank(f64),
    /// The gene's copies, the copies of each exon, the flank depth.
    Copies {
        copies: f64,
        exons: Vec<(i64, f64)>,
        flank: f64,
    },
}

impl RhdDepth {
    /// The flank depth, where one was measured.
    fn flank(&self) -> Option<f64> {
        match self {
            RhdDepth::LowFlank(f) | RhdDepth::Copies { flank: f, .. } => Some(*f),
            _ => None,
        }
    }
}

/// Assesses the RHD row on the bedcov means (`og_cn::engine::assess`).
pub fn rhd_depth(row: &Locus, regions: &Regions, depth: Option<&HashMap<String, f64>>) -> RhdDepth {
    let Some(depth) = depth else {
        return RhdDepth::NoBam;
    };
    let m: Means = regions
        .keys()
        .filter_map(|n| depth.get(n).map(|d| (n.clone(), *d)))
        .collect();
    let b = match engine::assess(row, &m, None, Sex::Unknown) {
        Ok(b) => b,
        Err(AssessError::BrokenControl(e) | AssessError::NeedsMapq1(e)) => {
            return RhdDepth::Broken(e);
        }
    };
    let flank = b.control_depth.unwrap_or(0.0);
    match (b.reason, b.estimate) {
        (Some(Reason::NoData), _) => RhdDepth::NoData,
        (Some(Reason::LowControlDepth), _) | (_, None) => RhdDepth::LowFlank(flank),
        (_, Some(copies)) => {
            let scale = row.depth().scale;
            let exons = row
                .details
                .iter()
                .filter_map(|s| {
                    let d = m.get(&s.name)?;
                    let e = s.name.strip_prefix("RHD_EX")?.parse().ok()?;
                    Some((e, scale * 2.0 / 2.0 * d / flank))
                })
                .collect();
            RhdDepth::Copies {
                copies,
                exons,
                flank,
            }
        }
    }
}

/// The RHD row's minimum flank depth, written as the old constant was (`10`).
fn min_flank(row: &Locus) -> String {
    let v = row.depth().min_control;
    if v == v.trunc() {
        format!("{}", v as i64)
    } else {
        pycompat::fmt::repr(v)
    }
}

fn system_rhd(sites: &Genotyped, rhd: &RhdDepth, row: &Locus) -> Result<Row, Error> {
    let (copies, exons) = match rhd {
        RhdDepth::NoBam => {
            return Ok(Row::not_determined(
                "RhD",
                "needs the BAM: the common RhD-negative allele (RHD*01N.01) deletes the whole RHD gene, which a VCF does not show",
            ));
        }
        RhdDepth::Broken(e) => return Ok(Row::not_determined("RhD", e.clone())),
        RhdDepth::NoData => {
            return Ok(Row::not_determined(
                "RhD",
                format!(
                    "no read depth over RHD and its flanks (the BAM has no {} contig)",
                    row.chrom
                ),
            ));
        }
        RhdDepth::LowFlank(f) => {
            return Ok(Row::not_determined(
                "RhD",
                format!(
                    "too few reads beside RHD ({}x) for a copy number (at least {}x)",
                    fixed(*f, 1),
                    min_flank(row)
                ),
            ));
        }
        RhdDepth::Copies { copies, exons, .. } => (*copies, exons),
    };
    let CallRule::Cut(cuts) = &row.call else {
        return Err(Error(format!(
            "the RHD row of the catalog calls by {:?}, not by cut points",
            row.call
        )));
    };
    let call = engine::cut_copies(copies, cuts);
    let (sd, s609) = (site(sites, "RHD_PSI_DUP"), site(sites, "RHD_PSI_609"));
    let (mut used, mut extra): (Vec<&Site>, Vec<String>) = (Vec::new(), Vec::new());
    let mut n_psi = 0;
    if call == 0 {
        // No RHD left: the pseudogene markers have nothing to say.
    } else if let (Some(a), Some(b)) = (sd.n, s609.n) {
        used.extend([sd, s609]);
        n_psi = a.min(b).min(call);
    } else {
        let why: Vec<String> = [sd, s609]
            .into_iter()
            .filter(|s| s.n.is_none())
            .map(|s| format!("{}: {}", s.row.hgvs, s.status.as_str()))
            .collect();
        extra.push(format!(
            "RHD pseudogene (RHD*08N.01) not checked: {}",
            why.join("; ")
        ));
    }
    let exon = |e: i64| exons.iter().find(|(x, _)| *x == e).map(|(_, c)| *c);
    let hybrid: Vec<i64> = RHD_HYBRID_EXONS
        .into_iter()
        .filter(|e| call >= 1 && exon(*e).unwrap_or(1.0) < 0.5)
        .collect();
    let mut sorted: Vec<(i64, f64)> = Vec::new();
    for &(e, c) in exons {
        // A repeated exon keeps its last value, at its first place, as a dict does.
        match sorted.iter_mut().find(|(x, _)| *x == e) {
            Some(slot) => slot.1 = c,
            None => sorted.push((e, c)),
        }
    }
    sorted.sort_by_key(|(e, _)| *e);
    let exon_txt: Vec<String> = sorted
        .iter()
        .map(|(e, c)| format!("exon {e} {}", fixed(*c, 2)))
        .collect();
    let geno: Vec<&str> = rep("RHD", call - n_psi)
        .chain(rep("RHD*08N.01", n_psi))
        .chain(rep("RHD deleted", 2 - call))
        .collect();
    let c2 = fixed(copies, 2);
    let basis = format!("{} (RHD copies by read depth {c2})", geno.join("/"));
    if !hybrid.is_empty() {
        let list: Vec<String> = hybrid.iter().map(i64::to_string).collect();
        return Ok(Row::new(
            "RhD",
            NOT_DETERMINED,
            basis,
            NOT_DETERMINED,
            format!(
                "RHD exon{} {} without RHD reads while the gene has {c2} copies: the pattern of an RHD-RHCE hybrid gene (such as RHD*03N.01) or a partial deletion. Copies per exon: {}",
                if hybrid.len() > 1 { "s" } else { "" },
                list.join(", "),
                exon_txt.join(", ")
            ),
        ));
    }
    let margin = (copies - 0.5).abs().min((copies - 1.5).abs());
    let mut conf = if !used.is_empty() {
        level_of(&used)
    } else if call == 0 {
        "high"
    } else {
        "medium"
    };
    if margin < 0.2 {
        conf = "low";
        extra.push(format!(
            "read depth {c2} copies is near a copy-number boundary"
        ));
    }
    if copies > 2.5 {
        extra.push(format!("read depth gives {c2} RHD copies, more than two"));
    }
    extra.push("weak and partial D (RHD missense alleles) are not assessed".to_owned());
    Ok(Row::new(
        "RhD",
        if call - n_psi > 0 { "D+" } else { "D-" },
        basis,
        conf,
        notes_of(&used, extra),
    ))
}

/// H on red cells needs FUT1: without it (two FUT1 null copies) the ABO
/// alleles make no A or B there. Names Bombay or para-Bombay from the
/// secretor row, and says so in the ABO row.
fn bombay(rows: &mut [Row]) {
    let at = |rows: &[Row], s: &str| rows.iter().position(|r| r.system == s);
    let (Some(h), Some(se), Some(abo)) = (at(rows, "H"), at(rows, "Secretor"), at(rows, "ABO"))
    else {
        return;
    };
    let (hs, ss) = (
        null_system("H").expect("H"),
        null_system("Secretor").expect("Secretor"),
    );
    if rows[h].phenotype != hs.minus {
        return;
    }
    if rows[se].phenotype == ss.minus {
        rows[h].phenotype = "H- (Bombay, Oh)".to_owned();
    } else if rows[se].phenotype == ss.plus {
        rows[h].phenotype = "H- (para-Bombay)".to_owned();
        rows[h].add_note(
            "para-Bombay: a working FUT2 makes H in secretions, and red cells carry little or none",
        );
    }
    rows[abo].add_note(
        "no FUT1 (H-): the red cells carry no H (or very little, para-Bombay), so these ABO alleles make little or no A or B on them and the cells usually group as O",
    );
}

/// Every site a rule reads by name, missing from the marker table.
fn missing_sites(sites: &Genotyped) -> Vec<&'static str> {
    let named = SYSTEMS
        .iter()
        .flat_map(|(_, n)| n.iter().flatten().copied())
        .chain(SPECIAL_SITES);
    let mut out: Vec<&'static str> = named.filter(|n| !sites.contains_key(*n)).collect();
    out.dedup();
    out
}

/// One row per system: the systems that lack a high-prevalence antigen
/// first, each list in the base order of [`SYSTEMS`].
pub fn predict(sites: &Genotyped, rhd: &RhdDepth, rhd_row: &Locus) -> Result<Vec<Row>, Error> {
    let gone = missing_sites(sites);
    if !gone.is_empty() {
        return Err(Error(format!(
            "the marker table has no row {}: the rules read it",
            gone.join(", ")
        )));
    }
    let mut rows = Vec::with_capacity(SYSTEMS.len());
    for (system, names) in SYSTEMS.iter() {
        rows.push(match (system.as_str(), names) {
            ("ABO", _) => system_abo(sites),
            ("RhD", _) => system_rhd(sites, rhd, rhd_row)?,
            ("Duffy", _) => system_duffy(sites),
            ("HNA-1", _) => system_hna1(sites),
            (s, Some(names)) => system_pairs(s, names, sites),
            (s, None) => system_null(null_system(s).expect("a system read from null rows"), sites),
        });
    }
    if let Some(se) = rows.iter_mut().find(|r| r.system == "Secretor")
        && se.phenotype != NOT_DETERMINED
    {
        se.add_note(
            "weak secretor alleles (such as FUT2 c.418A>T, also written c.385A>T) and a deletion of the whole FUT2 gene are not assessed",
        );
    }
    bombay(&mut rows);
    let (mut first, rest): (Vec<Row>, Vec<Row>) =
        rows.into_iter().partition(|r| !r.lacks_high.is_empty());
    first.extend(rest);
    Ok(first)
}

/// The evidence table's rows: one per site, then one per RhD depth region,
/// in the column order of [`crate::EVIDENCE_HEADER`].
pub fn evidence_rows(
    sites: &Genotyped,
    regions: &Regions,
    depth: Option<&HashMap<String, f64>>,
    rhd: &RhdDepth,
) -> Vec<[String; 12]> {
    let dot = |x: Option<String>| x.unwrap_or_else(|| ".".to_owned());
    let mut out = Vec::with_capacity(sites.len() + regions.len());
    for s in sites.values() {
        let r = &s.row;
        out.push([
            r.system.clone(),
            r.name.clone(),
            r.hgvs.clone(),
            r.rsid.clone(),
            r.chrom.clone(),
            r.pos.to_string(),
            r.ref_.clone(),
            r.alt.clone(),
            s.status.as_str().to_owned(),
            dot(s.n.map(|n| n.to_string())),
            dot(s.depth.map(|d| fixed(d, 1))),
            s.note.clone(),
        ]);
    }
    let flank = rhd.flank().filter(|f| *f != 0.0);
    for (name, r) in regions {
        let d = depth.and_then(|m| m.get(name).copied());
        let cpy = match (d, flank) {
            (Some(d), Some(f)) if !name.starts_with("RHD_FLANK") => Some(fixed(2.0 * d / f, 2)),
            _ => None,
        };
        let note = if cpy.is_some() {
            "copies = 2 x depth / flank depth"
        } else {
            ""
        };
        out.push([
            "RhD".to_owned(),
            name.clone(),
            r.label.clone(),
            String::new(),
            r.chrom.clone(),
            format!("{}-{}", r.start + 1, r.end),
            String::new(),
            String::new(),
            if d.is_some() { "depth" } else { "no depth" }.to_owned(),
            dot(cpy),
            dot(d.map(|d| fixed(d, 1))),
            note.to_owned(),
        ]);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn systems_in_report_order() {
        let names: Vec<&str> = SYSTEMS.iter().map(|(s, _)| s.as_str()).collect();
        assert_eq!(names.len(), 55);
        assert_eq!(&names[..3], ["ABO", "RhD", "RhCE"]);
        let at = |s: &str| names.iter().position(|n| *n == s).unwrap();
        assert!(at("HPA-14") < at("HPA-15") && at("HPA-15") < at("HPA-16"));
        assert_eq!(
            &names[names.len() - 4..],
            ["HNA-1", "HNA-3", "HNA-4", "HNA-5"]
        );
    }

    #[test]
    fn rare_hpa_templates() {
        assert_eq!(pair("HPA21_1960").unwrap().template, "HPA-21(a{21a}b{21b})");
        assert_eq!(pair("HPA21_1960").unwrap().label, "HPA-21");
    }

    #[test]
    fn join_skips_empty_parts() {
        assert_eq!(join_parts(["", "a", "", "b"]), "a; b");
    }
}
