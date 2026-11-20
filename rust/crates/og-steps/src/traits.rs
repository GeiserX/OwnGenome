//! Step 42: bounded site genotypes, haplotype rules and read-backed phase.
use indexmap::IndexMap;
use og_io::{AtomicFile, Lines, Utf8};
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::io::Write;

const COLUMNS: &[&str] = &[
    "chrom", "pos", "rsid", "gene", "ref", "alt", "trait", "hom_ref", "het", "hom_alt", "ancestry",
    "note", "pmids",
];
const OUT: &[&str] = &[
    "rsid", "gene", "chrom", "pos", "trait", "genotype", "status", "reading", "ancestry", "note",
    "pmids",
];
const RULE_OUT: &[&str] = &[
    "rule",
    "gene",
    "trait",
    "rsids",
    "chrom",
    "pos",
    "genotypes",
    "gq",
    "phase",
    "result",
    "status",
    "reading",
    "ancestry",
    "note",
    "pmids",
];
type Row = BTreeMap<String, String>;
type Key = (String, String);
type Genotypes = HashMap<Key, SiteGt>;
type Phased = HashMap<Key, (String, String, String)>;
type Result<T> = std::result::Result<T, String>;

#[derive(Clone, Debug)]
struct Site {
    row: Row,
    line: usize,
}
impl Site {
    fn key(&self) -> Key {
        (self.row["chrom"].clone(), self.row["pos"].clone())
    }
}
#[derive(Clone, Debug)]
struct Hap {
    name: String,
    alleles: Vec<String>,
    line: usize,
}
#[derive(Clone, Debug)]
struct Combo {
    key: String,
    label: String,
    reading: String,
    line: usize,
}
#[derive(Clone, Debug)]
struct Rule {
    id: String,
    meta: Row,
    line: usize,
    sites: Vec<Site>,
    haps: Vec<Hap>,
    combos: Vec<Combo>,
}
/// A call whose GQ is below this is not read (decision D10, the floor of
/// `og pgx refblock-fill` and the Y private SNPs). A missing GQ stays readable.
const MIN_GQ: u32 = 20;
/// The FILTER values that keep a call readable: a pass, no filter (a
/// reference block), and DeepVariant's RefCall (a candidate it called 0/0).
const READABLE_FILTERS: &[&str] = &["", ".", "PASS", "RefCall"];
#[derive(Clone, Debug, PartialEq, Eq)]
enum Status {
    Called,
    NotCovered,
    OtherAllele,
    LowGq,
    Unphased,
    OtherHaplotype,
}
impl Status {
    fn text(&self) -> &'static str {
        match self {
            Self::Called => "called",
            Self::NotCovered => "not_covered",
            Self::OtherAllele => "other_allele",
            Self::LowGq => "low_gq",
            Self::Unphased => "unphased",
            Self::OtherHaplotype => "other_haplotype",
        }
    }
}
/// Why a call is too weak to read.
#[derive(Clone, Debug, PartialEq, Eq)]
enum Weak {
    /// Its GQ, below [`MIN_GQ`].
    Gq(String),
    /// The caller's FILTER value.
    Filtered(String),
}
impl Weak {
    /// The record's GQ and FILTER columns: `None` when the call is readable.
    fn of(gq: &str, filter: &str) -> Option<Self> {
        if !READABLE_FILTERS.contains(&filter) {
            return Some(Self::Filtered(filter.to_owned()));
        }
        gq.parse::<f64>()
            .is_ok_and(|q| q.is_finite() && q < f64::from(MIN_GQ))
            .then(|| Self::Gq(gq.to_owned()))
    }
    fn why(&self) -> String {
        match self {
            Self::Gq(q) => format!("genotype quality {q}, below {MIN_GQ}"),
            Self::Filtered(f) => format!("the caller filtered the call: {f}"),
        }
    }
}
#[derive(Clone, Debug)]
struct SiteGt {
    alleles: Vec<String>,
    status: Status,
    gq: String,
    /// Set when the call is too weak to read: the status is then LowGq.
    weak: Option<Weak>,
    snv: bool,
}
impl Default for SiteGt {
    fn default() -> Self {
        Self {
            alleles: vec![],
            status: Status::NotCovered,
            gq: ".".into(),
            weak: None,
            snv: false,
        }
    }
}
fn row(keys: &[&str], values: &[String]) -> Row {
    keys.iter()
        .zip(values)
        .map(|(k, v)| (k.to_string(), v.clone()))
        .collect()
}
fn records(path: &str) -> Result<Vec<(usize, Vec<String>)>> {
    // Catalogs are plain TSV: quotes are literal and never hide a tab.
    let mut lines = Lines::new(
        og_io::open(path).map_err(|e| og_filters::os_error(&e, path))?,
        Utf8::Strict,
    );
    let mut text = String::new();
    let mut number = 0;
    let mut out = Vec::new();
    while lines.read_line(&mut text).map_err(|e| e.to_string())? {
        number += 1;
        if text.trim_end_matches('\n').is_empty() || text.starts_with('#') {
            continue;
        }
        let cells = text
            .trim_end_matches('\n')
            .split('\t')
            .map(str::to_owned)
            .collect();
        out.push((number, cells));
    }
    Ok(out)
}
fn read_table(path: &str) -> Result<Vec<Site>> {
    let mut input = records(path)?.into_iter();
    let (n, header) = input
        .next()
        .ok_or_else(|| format!("{path}: no header line"))?;
    if header != COLUMNS {
        return Err(format!("{path}:{n}: unexpected traits header"));
    }
    input
        .map(|(line, c)| {
            if c.len() != COLUMNS.len() {
                return Err(format!(
                    "{path}:{line}: {} columns, expected {}",
                    c.len(),
                    COLUMNS.len()
                ));
            }
            Ok(Site {
                row: row(COLUMNS, &c),
                line,
            })
        })
        .collect()
}
fn read_rules(path: &str) -> Result<Vec<Rule>> {
    let mut rules: IndexMap<String, Rule> = IndexMap::new();
    for (line, c) in records(path)? {
        let fields: &[&str] = match c[0].as_str() {
            "rule" => &["gene", "trait", "needs_phase", "ancestry", "note", "pmids"],
            "site" => &["rsid", "chrom", "pos", "ref", "alt"],
            "haplotype" => &["name", "alleles"],
            "combination" => &["diplotypes", "label", "reading"],
            kind => {
                return Err(format!(
                    "{path}:{line}: kind '{kind}' is not one of rule, site, haplotype, combination"
                ));
            }
        };
        if c.len() != fields.len() + 2 {
            return Err(format!(
                "{path}:{line}: a {} line has {} columns, expected {}",
                c[0],
                c.len(),
                fields.len() + 2
            ));
        }
        let id = &c[1];
        if c[0] == "rule" {
            if rules.contains_key(id) {
                return Err(format!("{path}:{line}: rule {id} has two rule lines"));
            }
            rules.insert(
                id.clone(),
                Rule {
                    id: id.clone(),
                    meta: row(fields, &c[2..]),
                    line,
                    sites: vec![],
                    haps: vec![],
                    combos: vec![],
                },
            );
        } else {
            let r = rules.get_mut(id).ok_or_else(|| {
                format!(
                    "{path}:{line}: a {} line of rule '{id}' before its rule line",
                    c[0]
                )
            })?;
            match c[0].as_str() {
                "site" => r.sites.push(Site {
                    row: row(fields, &c[2..]),
                    line,
                }),
                "haplotype" => r.haps.push(Hap {
                    name: c[2].clone(),
                    alleles: c[3].split(',').map(str::to_owned).collect(),
                    line,
                }),
                _ => r.combos.push(Combo {
                    key: c[2].clone(),
                    label: c[3].clone(),
                    reading: c[4].clone(),
                    line,
                }),
            }
        }
    }
    if rules.is_empty() {
        return Err(format!("{path}: no rule line"));
    }
    Ok(rules.into_values().collect())
}
fn positive(s: &str) -> bool {
    s.as_bytes()
        .first()
        .is_some_and(|c| matches!(c, b'1'..=b'9'))
        && s.bytes().all(|c| c.is_ascii_digit())
}
fn pmids(s: &str) -> bool {
    s.split(',').all(positive)
}
fn text_errors(w: &str, text: &str, e: &mut Vec<String>) {
    if text.to_lowercase().contains("snpedia") {
        e.push(format!(
            "{w}: mentions SNPedia; the table takes no text from it"
        ));
    }
    if text.contains('"') || !text.is_ascii() {
        e.push(format!(
            "{w}: a double quote or a non-ASCII character (the reports read the table as plain TSV)"
        ));
    }
}
fn site_errors(w: &str, s: &Site, e: &mut Vec<String>) {
    let r = &s.row;
    if !(1..=22).any(|n| r["chrom"] == format!("chr{n}")) {
        e.push(format!("{w}: chrom '{}' is not chr1 to chr22", r["chrom"]));
    }
    if !positive(&r["pos"]) {
        e.push(format!("{w}: pos '{}' is not a positive integer", r["pos"]));
    }
    if !r["rsid"].strip_prefix("rs").is_some_and(positive) {
        e.push(format!("{w}: rsid is not rs<number>"));
    }
    for k in ["ref", "alt"] {
        if !["A", "C", "G", "T"].contains(&r[k].as_str()) {
            e.push(format!("{w}: {k} '{}' is not one of A, C, G, T", r[k]));
        }
    }
    if r["ref"] == r["alt"] {
        e.push(format!("{w}: ref and alt are the same base"));
    }
}
fn check_rows(rows: &[Site]) -> Vec<String> {
    let mut e = vec![];
    let mut ids = HashSet::new();
    let mut sites = HashSet::new();
    for s in rows {
        let r = &s.row;
        let w = format!("line {} ({})", s.line, r["rsid"]);
        site_errors(&w, s, &mut e);
        if !ids.insert(r["rsid"].as_str()) {
            e.push(format!("{w}: rsid appears twice"));
        }
        if !sites.insert(s.key()) {
            e.push(format!(
                "{w}: position {}:{} appears twice",
                r["chrom"], r["pos"]
            ));
        }
        for k in ["gene", "trait", "hom_ref", "het", "hom_alt", "note"] {
            if r[k].trim().is_empty() {
                e.push(format!("{w}: {k} is empty"));
            }
        }
        if !["yes", "no"].contains(&r["ancestry"].as_str()) {
            e.push(format!(
                "{w}: ancestry '{}' is not yes or no",
                r["ancestry"]
            ));
        }
        if !pmids(&r["pmids"]) {
            e.push(format!(
                "{w}: pmids '{}' is not a comma-separated list of PubMed ids",
                r["pmids"]
            ));
        }
        text_errors(
            &w,
            &COLUMNS
                .iter()
                .map(|k| r[*k].as_str())
                .collect::<Vec<_>>()
                .join(" "),
            &mut e,
        );
    }
    for (id, a, b) in [("rs429358", "T", "C"), ("rs7412", "C", "T")] {
        if let Some(s) = rows.iter().rev().find(|s| s.row["rsid"] == id) {
            let r = &s.row;
            if r["chrom"] != "chr19" || r["ref"] != a || r["alt"] != b {
                e.push(format!(
                    "APOE row {id} must be chr19 {a}>{b}, is {} {}>{}",
                    r["chrom"], r["ref"], r["alt"]
                ));
            }
        } else {
            e.push(format!("APOE row {id} is missing"));
        }
    }
    e
}
fn allowed(rule: &Rule, gt: &[Vec<String>]) -> Vec<String> {
    let mut out = vec![];
    for (i, a) in rule.haps.iter().enumerate() {
        for b in &rule.haps[i..] {
            if gt.iter().enumerate().all(|(k, g)| {
                g.len() == 2
                    && ((g[0] == a.alleles[k] && g[1] == b.alleles[k])
                        || (g[1] == a.alleles[k] && g[0] == b.alleles[k]))
            }) {
                out.push(format!("{}/{}", a.name, b.name));
            }
        }
    }
    out
}
fn check_rules(rules: &[Rule], rows: Option<&[Site]>) -> Vec<String> {
    let mut e = vec![];
    for r in rules {
        let n0 = e.len();
        let w = format!("rule {} (line {})", r.id, r.line);
        if r.id.is_empty() || !r.id.bytes().all(|b| b.is_ascii_alphanumeric() || b == b'_') {
            e.push(format!("{w}: the rule id is not letters, digits and _"));
        }
        for k in ["gene", "trait", "note"] {
            if r.meta[k].trim().is_empty() {
                e.push(format!("{w}: {k} is empty"));
            }
        }
        for k in ["needs_phase", "ancestry"] {
            if !["yes", "no"].contains(&r.meta[k].as_str()) {
                e.push(format!("{w}: {k} '{}' is not yes or no", r.meta[k]));
            }
        }
        if !pmids(&r.meta["pmids"]) {
            e.push(format!(
                "{w}: pmids '{}' is not a comma-separated list of PubMed ids",
                r.meta["pmids"]
            ));
        }
        if !(2..=8).contains(&r.sites.len()) {
            e.push(format!("{w}: {} sites, expected 2 to 8", r.sites.len()));
        }
        let mut ids = HashSet::new();
        let mut positions = HashSet::new();
        for s in &r.sites {
            let sw = format!("{w}, site line {} ({})", s.line, s.row["rsid"]);
            site_errors(&sw, s, &mut e);
            if !ids.insert(&s.row["rsid"]) || !positions.insert(s.key()) {
                e.push(format!("{sw}: the site appears twice in the rule"));
            }
            if let Some(t) =
                rows.and_then(|rows| rows.iter().rev().find(|t| t.row["rsid"] == s.row["rsid"]))
                && ["chrom", "pos", "ref", "alt"]
                    .iter()
                    .any(|k| t.row[*k] != s.row[*k])
            {
                e.push(format!(
                    "{sw}: the traits table has {}:{} {}>{}, the rule {}:{} {}>{}",
                    t.row["chrom"],
                    t.row["pos"],
                    t.row["ref"],
                    t.row["alt"],
                    s.row["chrom"],
                    s.row["pos"],
                    s.row["ref"],
                    s.row["alt"]
                ));
            }
        }
        if r.sites
            .iter()
            .map(|s| &s.row["chrom"])
            .collect::<HashSet<_>>()
            .len()
            > 1
        {
            e.push(format!("{w}: the sites are on more than one chromosome"));
        }
        let mut names = Vec::new();
        let mut tuples = HashSet::new();
        for h in &r.haps {
            let hw = format!("{w}, haplotype line {} ({})", h.line, h.name);
            if h.name.is_empty()
                || !h
                    .name
                    .bytes()
                    .all(|b| b.is_ascii_alphanumeric() || b"*.".contains(&b))
            {
                e.push(format!("{hw}: the name is not letters, digits, * and ."));
            }
            if names.contains(&h.name.as_str()) {
                e.push(format!("{hw}: the name appears twice"));
            }
            names.push(h.name.as_str());
            if h.alleles.len() != r.sites.len() {
                e.push(format!(
                    "{hw}: {} alleles for {} sites",
                    h.alleles.len(),
                    r.sites.len()
                ));
                continue;
            }
            for (a, s) in h.alleles.iter().zip(&r.sites) {
                if a != &s.row["ref"] && a != &s.row["alt"] {
                    e.push(format!(
                        "{hw}: allele {a} at {} is neither its ref {} nor its alt {}",
                        s.row["rsid"], s.row["ref"], s.row["alt"]
                    ));
                }
            }
            if !tuples.insert(&h.alleles) {
                e.push(format!("{hw}: the same alleles as another haplotype"));
            }
        }
        if names.len() < 2 {
            e.push(format!("{w}: fewer than two haplotypes"));
        }
        let mut keys = BTreeSet::new();
        for c in &r.combos {
            let cw = format!("{w}, combination line {} ({})", c.line, c.key);
            for (k, v) in [("label", &c.label), ("reading", &c.reading)] {
                if v.trim().is_empty() {
                    e.push(format!("{cw}: {k} is empty"));
                }
            }
            text_errors(&cw, &format!("{} {} {}", c.key, c.label, c.reading), &mut e);
            if !keys.insert(c.key.clone()) {
                e.push(format!("{cw}: the diplotypes appear twice"));
            }
            for d in c.key.split('|') {
                let pair: Vec<_> = d.split('/').collect();
                if pair.len() != 2 || pair.iter().any(|p| !names.contains(p)) {
                    e.push(format!(
                        "{cw}: '{d}' is not two of the rule's haplotypes joined by /"
                    ));
                } else if names.iter().position(|n| n == &pair[0])
                    > names.iter().position(|n| n == &pair[1])
                {
                    e.push(format!(
                        "{cw}: '{d}' is not in the order of the haplotype lines (write {}/{})",
                        pair[1], pair[0]
                    ));
                }
            }
        }
        text_errors(
            &w,
            &format!("{} {} {}", r.meta["gene"], r.meta["trait"], r.meta["note"]),
            &mut e,
        );
        if e.len() != n0 {
            continue;
        }
        let mut want = BTreeSet::new();
        if r.meta["needs_phase"] == "yes" {
            for (i, a) in names.iter().enumerate() {
                for b in &names[i..] {
                    want.insert(format!("{a}/{b}"));
                }
            }
        } else {
            for mut n in 0..3usize.pow(r.sites.len() as u32) {
                let mut gt = vec![];
                for s in &r.sites {
                    let v = n % 3;
                    n /= 3;
                    gt.push(vec![
                        s.row[if v == 2 { "alt" } else { "ref" }].clone(),
                        s.row[if v == 0 { "ref" } else { "alt" }].clone(),
                    ]);
                }
                let key = allowed(r, &gt).join("|");
                if !key.is_empty() {
                    want.insert(key);
                }
            }
        }
        for k in want.difference(&keys) {
            e.push(format!("{w}: no combination line for {k}"));
        }
        for k in keys.difference(&want) {
            e.push(format!(
                "{w}: the combination line {k} can never be read {}",
                if k.contains('|') {
                    "(a rule that needs phase reads one diplotype)"
                } else {
                    "from these haplotypes"
                }
            ));
        }
    }
    if let Some(rows) = rows {
        if let Some(r) = rules.iter().find(|r| r.id == "APOE") {
            for s in &r.sites {
                if !rows.iter().any(|t| t.row["rsid"] == s.row["rsid"]) {
                    e.push(format!(
                        "rule APOE: site {} has no traits table row",
                        s.row["rsid"]
                    ));
                }
            }
        } else {
            e.push("rule APOE is missing: its line goes into the traits table, where the reports read it".into());
        }
    }
    e
}

// A locus may be reused by rules with different alternate alleles. The
// selected record is shared, but each site's allele vocabulary is its own.
fn genotype(site: &Site, records: &Genotypes) -> SiteGt {
    let mut gt = records.get(&site.key()).cloned().unwrap_or_default();
    if gt.status == Status::NotCovered {
        return gt;
    }
    gt.alleles.sort_by_key(|a| {
        (
            if a == &site.row["ref"] {
                0
            } else if a == &site.row["alt"] {
                1
            } else {
                2
            },
            a.clone(),
        )
    });
    gt.status = if gt.weak.is_some() {
        Status::LowGq
    } else if gt.alleles.len() == 2
        && gt
            .alleles
            .iter()
            .all(|a| a == &site.row["ref"] || a == &site.row["alt"])
    {
        Status::Called
    } else {
        Status::OtherAllele
    };
    gt
}

fn read_genotypes(path: &str, sites: &[&Site]) -> Result<Genotypes> {
    let mut wanted: HashMap<Key, Vec<&Site>> = HashMap::new();
    for s in sites {
        wanted.entry(s.key()).or_default().push(s);
    }
    let mut result = HashMap::new();
    let mut lines = Lines::new(
        og_io::open(path).map_err(|e| og_filters::os_error(&e, path))?,
        Utf8::Strict,
    );
    let mut line = String::new();
    let mut n = 0;
    while lines.read_line(&mut line).map_err(|e| e.to_string())? {
        n += 1;
        let line = line.trim_end_matches('\n');
        if line.is_empty() {
            continue;
        }
        let c: Vec<_> = line.split('\t').collect();
        if c.len() < 5 {
            return Err(format!(
                "{path}:{n}: {} columns, expected CHROM POS REF ALT GENOTYPE [GQ [FILTER]]",
                c.len()
            ));
        }
        let key = (c[0].to_owned(), c[1].to_owned());
        let Some(sites) = wanted.get(&key) else {
            continue;
        };
        let site = sites[0];
        for s in sites {
            if c[2].chars().next() != s.row["ref"].chars().next() {
                return Err(format!(
                    "{} at {}:{}: the gVCF has REF {}, the table {} (the GRCh38 base). Is the gVCF on GRCh38?",
                    s.row["rsid"], c[0], c[1], c[2], s.row["ref"]
                ));
            }
        }
        let mut al: Vec<_> = c[4].split(['/', '|']).map(str::to_owned).collect();
        let snv = al.iter().all(|a| a.chars().count() == 1);
        if result.get(&key).is_some_and(|old: &SiteGt| old.snv || !snv) {
            continue;
        }
        al.sort_by_key(|a| {
            (
                if a == &site.row["ref"] {
                    0
                } else if a == &site.row["alt"] {
                    1
                } else {
                    2
                },
                a.clone(),
            )
        });
        let gq = c
            .get(5)
            .filter(|g| !g.is_empty())
            .unwrap_or(&".")
            .to_string();
        let weak = Weak::of(&gq, c.get(6).copied().unwrap_or("."));
        let status = if weak.is_some() {
            Status::LowGq
        } else if al.len() == 2
            && al
                .iter()
                .all(|a| a == &site.row["ref"] || a == &site.row["alt"])
        {
            Status::Called
        } else {
            Status::OtherAllele
        };
        result.insert(
            key,
            SiteGt {
                alleles: al,
                status,
                gq,
                weak,
                snv,
            },
        );
    }
    Ok(result)
}
fn read_phased(path: &str, wanted: &HashSet<Key>) -> Result<Phased> {
    let mut lines = Lines::new(
        og_io::open(path).map_err(|e| og_filters::os_error(&e, path))?,
        Utf8::Strict,
    );
    let mut text = String::new();
    let mut out = HashMap::new();
    while lines.read_line(&mut text).map_err(|e| e.to_string())? {
        if text.starts_with('#') {
            continue;
        }
        let c: Vec<_> = text.trim_end_matches('\n').split('\t').collect();
        if c.len() < 10 {
            continue;
        }
        let key = (c[0].to_owned(), c[1].to_owned());
        if !wanted.contains(&key) {
            continue;
        }
        let fmt: HashMap<_, _> = c[8].split(':').zip(c[9].split(':')).collect();
        let (Some(gt), Some(ps)) = (fmt.get("GT"), fmt.get("PS")) else {
            continue;
        };
        if ps.is_empty() || *ps == "." || !gt.contains('|') {
            continue;
        }
        let bases: Vec<_> = std::iter::once(c[3]).chain(c[4].split(',')).collect();
        let idx: Option<Vec<_>> = gt
            .split('|')
            .map(|i| {
                pycompat::num::parse_int(i)
                    .ok()
                    .and_then(|i| i.to_string().parse::<isize>().ok())
                    .and_then(|i| {
                        let n = if i < 0 { bases.len() as isize + i } else { i };
                        usize::try_from(n).ok().and_then(|n| bases.get(n).copied())
                    })
            })
            .collect();
        if let Some(al) = idx
            && al.len() == 2
        {
            out.insert(key, (al[0].into(), al[1].into(), ps.to_string()));
        }
    }
    Ok(out)
}
/// The two haplotypes of a key that is one diplotype (`a/b`). An empty key
/// or several diplotypes (`a/b|c/d`, whose reading already says what phase
/// would decide) give `None`.
fn assumed_copies<'a>(r: &'a Rule, key: &str) -> Option<(&'a Hap, &'a Hap)> {
    if key.is_empty() || key.contains('|') {
        return None;
    }
    let (a, b) = key.split_once('/')?;
    let hap = |n: &str| r.haps.iter().find(|h| h.name == n);
    Some((hap(a)?, hap(b)?))
}
/// "x and y", or "x, y and z".
fn join_and(items: &[String]) -> String {
    match items {
        [] => String::new(),
        [one] => one.clone(),
        [rest @ .., last] => format!("{} and {last}", rest.join(", ")),
    }
}
/// Where the changed alleles (the alt of each heterozygous site) sit, given
/// the alleles of one copy: the other copy carries the rest.
fn placement(r: &Rule, het: &[usize], copy: &[String]) -> String {
    let alt = |i: &usize| format!("{} {}", r.sites[*i].row["rsid"], r.sites[*i].row["alt"]);
    let (on, off): (Vec<usize>, Vec<usize>) = het
        .iter()
        .partition(|i| copy[**i] == r.sites[**i].row["alt"]);
    let all = join_and(&het.iter().map(alt).collect::<Vec<_>>());
    if on.is_empty() || off.is_empty() {
        format!("{all} sit on one copy")
    } else if het.len() == 2 {
        format!("{all} sit on different copies")
    } else {
        format!(
            "{} sit on one copy and {} on the other",
            join_and(&on.iter().map(alt).collect::<Vec<_>>()),
            join_and(&off.iter().map(alt).collect::<Vec<_>>())
        )
    }
}
/// The reading when read-backed phase puts the alleles on copies the rule
/// does not describe, where the genotypes alone read `label`.
fn phased_apart(
    r: &Rule,
    het: &[usize],
    copies: [(&Vec<String>, Option<usize>); 2],
    label: &str,
) -> String {
    let rs: Vec<_> = r.sites.iter().map(|s| s.row["rsid"].as_str()).collect();
    let name = |(h, n): (&Vec<String>, Option<usize>)| {
        n.map_or_else(
            || format!("{} at {}", h.join(","), rs.join(", ")),
            |n| r.haps[n].name.clone(),
        )
    };
    let changed = copies
        .iter()
        .filter(|(h, _)| h.iter().zip(&r.sites).any(|(a, s)| a == &s.row["alt"]))
        .count();
    let unknown = copies.iter().filter(|(_, n)| n.is_none()).count();
    format!(
        "Read-backed phase shows that {}: the copies are {} and {}{}. This rule does not describe {}, so the pair is not read; without phase these genotypes would read {label}.",
        placement(r, het, copies[0].0),
        name(copies[0]),
        name(copies[1]),
        if changed == 1 {
            ", so one copy is changed and the other is not"
        } else {
            ""
        },
        if unknown == 1 {
            "that copy"
        } else {
            "those copies"
        }
    )
}
fn rule_line(r: &Rule, gt: &Genotypes, phased: Option<&Phased>) -> Row {
    let g: Vec<_> = r.sites.iter().map(|s| genotype(s, gt)).collect();
    let rs: Vec<_> = r.sites.iter().map(|s| s.row["rsid"].as_str()).collect();
    let mut out = r.meta.clone();
    for (k, v) in [
        ("rule", r.id.clone()),
        ("rsids", rs.join(",")),
        ("chrom", r.sites[0].row["chrom"].clone()),
        (
            "pos",
            r.sites
                .iter()
                .map(|s| s.row["pos"].as_str())
                .collect::<Vec<_>>()
                .join(","),
        ),
        (
            "genotypes",
            g.iter()
                .map(|g| {
                    if g.alleles.is_empty() {
                        ".".into()
                    } else {
                        g.alleles.join("/")
                    }
                })
                .collect::<Vec<_>>()
                .join(","),
        ),
        (
            "gq",
            g.iter()
                .map(|g| g.gq.as_str())
                .collect::<Vec<_>>()
                .join(","),
        ),
    ] {
        out.insert(k.into(), v);
    }
    let finish = |mut out: Row, phase: &str, result: &str, status: Status, reading: String| {
        for (k, v) in [
            ("phase", phase),
            ("result", result),
            ("status", status.text()),
            ("reading", &reading),
        ] {
            out.insert(k.into(), v.into());
        }
        out
    };
    if g.iter().any(|g| g.status != Status::Called) {
        let need = if rs.len() == 2 {
            format!("both {} and {}", rs[0], rs[1])
        } else {
            format!(
                "every one of {} and {}",
                rs[..rs.len() - 1].join(", "),
                rs[rs.len() - 1]
            )
        };
        let weak: Vec<_> = r
            .sites
            .iter()
            .zip(&g)
            .filter_map(|(s, g)| g.weak.as_ref().map(|w| (s, w)))
            .collect();
        let (status, reading) = if g.iter().any(|g| g.status == Status::NotCovered) {
            (
                Status::NotCovered,
                format!("Not read: needs a genotype at {need}."),
            )
        } else if !weak.is_empty() {
            (
                Status::LowGq,
                format!(
                    "A low-quality call at {}: not read; check the reads.",
                    weak.iter()
                        .map(|(s, w)| format!("{} ({})", s.row["rsid"], w.why()))
                        .collect::<Vec<_>>()
                        .join(" and ")
                ),
            )
        } else {
            (
                Status::OtherAllele,
                format!("Not read: needs a genotype at {need}."),
            )
        };
        return finish(out, ".", ".", status, reading);
    }
    let mut phase = if r.meta["needs_phase"] == "no" {
        "not_used"
    } else {
        "not_needed"
    };
    let het: Vec<_> = g
        .iter()
        .enumerate()
        .filter(|(_, g)| g.alleles[0] != g.alleles[1])
        .map(|(i, _)| i)
        .collect();
    // The two copies read-backed phase gives, when one phase block holds
    // every heterozygous site; the homozygous sites are on both copies.
    let from_reads = || -> Option<(Vec<String>, Vec<String>)> {
        let calls: Vec<_> = het
            .iter()
            .map(|i| phased.and_then(|p| p.get(&r.sites[*i].key())))
            .collect::<Option<_>>()?;
        if !calls.iter().all(|v| v.2 == calls[0].2) {
            return None;
        }
        let mut h1: Vec<_> = g.iter().map(|g| g.alleles[0].clone()).collect();
        let mut h2: Vec<_> = g.iter().map(|g| g.alleles[1].clone()).collect();
        for (i, c) in het.iter().zip(calls) {
            h1[*i] = c.0.clone();
            h2[*i] = c.1.clone();
        }
        Some((h1, h2))
    };
    let why = if phased.is_some() {
        "no read or read pair covered them together"
    } else {
        "no phasing ran (no BAM, or phasing failed)"
    };
    let mut caveat = String::new();
    let key = if r.meta["needs_phase"] == "no" {
        let key = allowed(r, &g.iter().map(|g| g.alleles.clone()).collect::<Vec<_>>()).join("|");
        // One diplotype read without phase from two or more heterozygous
        // sites assumes which alleles share a copy: read-backed phase checks
        // it, and without it the reading says what it assumes.
        let assumed = assumed_copies(r, &key);
        match (het.len() >= 2, assumed) {
            (true, Some((a, b))) => match from_reads() {
                Some((h1, h2)) => {
                    phase = "read_backed";
                    let n1 = r.haps.iter().position(|h| h.alleles == h1);
                    let n2 = r.haps.iter().position(|h| h.alleles == h2);
                    if let (Some(x), Some(y)) = (n1, n2) {
                        format!("{}/{}", r.haps[x.min(y)].name, r.haps[x.max(y)].name)
                    } else {
                        let label = r
                            .combos
                            .iter()
                            .rev()
                            .find(|c| c.key == key)
                            .map_or(key.as_str(), |c| c.label.as_str());
                        return finish(
                            out,
                            phase,
                            ".",
                            Status::OtherHaplotype,
                            phased_apart(r, &het, [(&h1, n1), (&h2, n2)], label),
                        );
                    }
                }
                None => {
                    caveat = format!(
                        " Phase unknown: this reading assumes {} (the copies {} and {}); {why}, so that is not checked.",
                        placement(r, &het, &a.alleles),
                        a.name,
                        b.name
                    );
                    key
                }
            },
            _ => key,
        }
    } else {
        let mut h1: Vec<_> = g.iter().map(|g| g.alleles[0].clone()).collect();
        let mut h2: Vec<_> = g.iter().map(|g| g.alleles[1].clone()).collect();
        if het.len() >= 2 {
            if let Some((p1, p2)) = from_reads() {
                (h1, h2) = (p1, p2);
                phase = "read_backed";
            } else {
                return finish(
                    out,
                    "unphased",
                    ".",
                    Status::Unphased,
                    format!(
                        "Heterozygous at {} of the {} sites, and {why}, so which alleles share a copy is not known: not read.",
                        het.len(),
                        r.sites.len()
                    ),
                );
            }
        }
        let n1 = r.haps.iter().position(|h| h.alleles == h1);
        let n2 = r.haps.iter().position(|h| h.alleles == h2);
        if let (Some(a), Some(b)) = (n1, n2) {
            format!("{}/{}", r.haps[a.min(b)].name, r.haps[a.max(b)].name)
        } else {
            let unknown = [(&h1, n1), (&h2, n2)]
                .into_iter()
                .filter(|(_, n)| n.is_none())
                .map(|(h, _)| h.join(","))
                .collect::<Vec<_>>()
                .join(" and ");
            return finish(
                out,
                phase,
                ".",
                Status::OtherHaplotype,
                format!(
                    "A haplotype this rule does not describe ({unknown} at {}): not read.",
                    rs.join(", ")
                ),
            );
        }
    };
    match r
        .combos
        .iter()
        .rev()
        .find(|c| !key.is_empty() && c.key == key)
    {
        Some(c) => finish(
            out,
            phase,
            &c.label,
            Status::Called,
            format!("{}{caveat}", c.reading),
        ),
        None => finish(
            out,
            phase,
            ".",
            Status::OtherHaplotype,
            "A combination of haplotypes this rule does not describe: not read.".into(),
        ),
    }
}
fn stage_tsv(path: &str, columns: &[&str], rows: &[Row]) -> Result<AtomicFile> {
    let mut file = AtomicFile::create(path).map_err(|e| format!("cannot write {path}: {e}"))?;
    writeln!(file, "{}", columns.join("\t")).map_err(|e| format!("cannot write {path}: {e}"))?;
    for row in rows {
        let fields: Vec<_> = columns.iter().map(|k| row[*k].as_str()).collect();
        writeln!(file, "{}", fields.join("\t")).map_err(|e| format!("cannot write {path}: {e}"))?;
    }
    file.flush()
        .map_err(|e| format!("cannot write {path}: {e}"))?;
    Ok(file)
}
fn phase_vcf(path: &str, sample: &str, rules: &[Rule], gt: &Genotypes) -> Result<usize> {
    let mut sites = BTreeMap::new();
    for r in rules {
        let calls: Option<Vec<_>> = r
            .sites
            .iter()
            .map(|s| {
                let g = genotype(s, gt);
                (g.status == Status::Called).then_some(g)
            })
            .collect();
        // A rule that reads without phase sends its sites only when the
        // genotypes read one diplotype that assumes which alleles share a
        // copy: read-backed phase then checks that assumption.
        let calls = calls.filter(|c| {
            r.meta["needs_phase"] == "yes"
                || assumed_copies(
                    r,
                    &allowed(r, &c.iter().map(|g| g.alleles.clone()).collect::<Vec<_>>()).join("|"),
                )
                .is_some()
        });
        if let Some(calls) = calls {
            let het: Vec<_> = r
                .sites
                .iter()
                .zip(calls)
                .filter(|(_, g)| g.alleles[0] != g.alleles[1])
                .map(|(s, _)| s)
                .collect();
            if het.len() >= 2 {
                for s in het {
                    let chrom = s.row["chrom"]
                        .strip_prefix("chr")
                        .and_then(|s| s.parse::<u32>().ok())
                        .ok_or("invalid chromosome")?;
                    let pos = s.row["pos"].parse::<u64>().map_err(|e| e.to_string())?;
                    sites.entry((chrom, pos)).or_insert(s);
                }
            }
        }
    }
    let mut out = AtomicFile::create(path).map_err(|e| format!("cannot write {path}: {e}"))?;
    let write = |out: &mut AtomicFile| -> std::io::Result<()> {
        writeln!(out, "##fileformat=VCFv4.2")?;
        for chrom in sites.keys().map(|(c, _)| c).collect::<BTreeSet<_>>() {
            writeln!(out, "##contig=<ID=chr{chrom}>")?;
        }
        writeln!(
            out,
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
        )?;
        writeln!(
            out,
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}"
        )?;
        for s in sites.values() {
            writeln!(
                out,
                "{}\t{}\t{}\t{}\t{}\t.\tPASS\t.\tGT\t0/1",
                s.row["chrom"], s.row["pos"], s.row["rsid"], s.row["ref"], s.row["alt"]
            )?;
        }
        Ok(())
    };
    write(&mut out).map_err(|e| format!("cannot write {path}: {e}"))?;
    out.commit()
        .map_err(|e| format!("cannot write {path}: {e}"))?;
    Ok(sites.len())
}
fn run(verb: &str, opts: &HashMap<String, String>) -> Result<u8> {
    let required = |k: &str| {
        opts.get(k)
            .map(String::as_str)
            .ok_or_else(|| format!("missing --{k}"))
    };
    if verb == "check" {
        let table = required("table")?;
        let rows = read_table(table)?;
        let mut errs = check_rows(&rows);
        if let Some(path) = opts.get("rules") {
            errs.extend(check_rules(&read_rules(path)?, Some(&rows)));
        }
        for e in &errs {
            eprintln!("ERROR: {e}");
        }
        if errs.is_empty() {
            println!(
                "{table}{}: schema OK",
                opts.get("rules")
                    .map(|r| format!(" and {r}"))
                    .unwrap_or_default()
            );
        }
        return Ok(u8::from(!errs.is_empty()));
    }
    let rules_path = required("rules")?;
    let rules = read_rules(rules_path)?;
    let rows = if verb == "table" {
        read_table(required("table")?)?
    } else {
        vec![]
    };
    let mut errs = if verb == "table" {
        check_rows(&rows)
            .into_iter()
            .map(|e| format!("{}: {e}", opts["table"]))
            .collect::<Vec<_>>()
    } else {
        vec![]
    };
    errs.extend(
        check_rules(&rules, if verb == "table" { Some(&rows) } else { None })
            .into_iter()
            .map(|e| format!("{rules_path}: {e}")),
    );
    if !errs.is_empty() {
        for e in errs {
            eprintln!("ERROR: {e}");
        }
        return Ok(1);
    }
    let sites: Vec<_> = rows
        .iter()
        .chain(rules.iter().flat_map(|r| &r.sites))
        .collect();
    let gt = read_genotypes(required("genotypes")?, &sites)?;
    let sample = required("sample")?;
    if verb == "phase-sites" {
        let path = required("out")?;
        let n = phase_vcf(path, sample, &rules, &gt)?;
        println!("{sample}: {n} heterozygous rule sites to phase. Written: {path}");
        return Ok(0);
    }
    let wanted = sites.iter().map(|s| s.key()).collect();
    let phased = opts
        .get("phased")
        .map(|p| read_phased(p, &wanted))
        .transpose()?;
    let lines: Vec<_> = rules
        .iter()
        .map(|r| rule_line(r, &gt, phased.as_ref()))
        .collect();
    let last = rules.iter().find(|r| r.id == "APOE").and_then(|r| {
        rows.iter()
            .rposition(|s| r.sites.iter().any(|t| s.row["rsid"] == t.row["rsid"]))
    });
    let mut out = vec![];
    let mut counts = BTreeMap::<String, usize>::new();
    for (i, s) in rows.iter().enumerate() {
        let g = genotype(s, &gt);
        let mut r = s.row.clone();
        r.insert(
            "genotype".into(),
            if g.alleles.is_empty() {
                ".".into()
            } else {
                g.alleles.join("/")
            },
        );
        r.insert("status".into(), g.status.text().into());
        let reading = match g.status {
            Status::NotCovered => {
                "No genotype at this position in the gVCF (no reads, or a no-call): not read."
                    .into()
            }
            Status::LowGq => format!(
                "A low-quality call at this position ({}): not read; check the reads.",
                g.weak.as_ref().map(Weak::why).unwrap_or_default()
            ),
            Status::Called => s.row[["hom_ref", "het", "hom_alt"]
                [g.alleles.iter().filter(|a| *a == &s.row["alt"]).count()]]
            .clone(),
            _ => "An allele this table does not describe at this position: not read.".into(),
        };
        r.insert("reading".into(), reading);
        out.push(r);
        if last == Some(i) {
            let mut r = lines
                .iter()
                .find(|r| r["rule"] == "APOE")
                .expect("validated APOE")
                .clone();
            r.insert("rsid".into(), r["rsids"].replace(',', "+"));
            r.insert("genotype".into(), r["result"].clone());
            out.push(r);
        }
    }
    // Compute both tables before opening either output. Each rename is atomic.
    let out_path = required("out")?;
    let rules_path = required("rules-out")?;
    // Refuse known rename failures before either previous output is replaced.
    // Publication still uses two renames, not one transaction.
    for path in [out_path, rules_path] {
        if std::path::Path::new(path).is_dir() {
            return Err(format!("cannot write {path}: destination is a directory"));
        }
    }
    let traits_file = stage_tsv(out_path, OUT, &out)?;
    let rules_file = stage_tsv(rules_path, RULE_OUT, &lines)?;
    traits_file
        .commit()
        .map_err(|e| format!("cannot write {out_path}: {e}"))?;
    rules_file
        .commit()
        .map_err(|e| format!("cannot write {rules_path}: {e}"))?;
    for r in &out {
        *counts.entry(r["status"].clone()).or_default() += 1;
    }
    println!(
        "{sample}: {} table rows; {} (APOE line included). Written: {}",
        rows.len(),
        counts
            .iter()
            .map(|(k, v)| format!("{k} {v}"))
            .collect::<Vec<_>>()
            .join(", "),
        opts["out"]
    );
    println!(
        "{sample}: {} genotype rules; {}. Written: {}",
        lines.len(),
        lines
            .iter()
            .map(|r| format!(
                "{} {}",
                r["rule"],
                r[if r["status"] == "called" {
                    "result"
                } else {
                    "status"
                }]
            ))
            .collect::<Vec<_>>()
            .join(", "),
        opts["rules-out"]
    );
    Ok(0)
}
/// Run the three traits verbs. Usage errors exit 2; data and I/O errors exit 1.
pub fn cli(args: &[String]) -> u8 {
    let usage = "usage: og traits table|phase-sites|check [--table TABLE] --rules RULES [--genotypes G --sample S --out OUT --rules-out OUT]";
    if args.iter().any(|a| a == "--help" || a == "-h") {
        println!("{usage}");
        return 0;
    }
    let Some(verb) = args
        .first()
        .map(String::as_str)
        .filter(|v| ["table", "phase-sites", "check"].contains(v))
    else {
        eprintln!("{usage}");
        return 2;
    };
    let mut opts = HashMap::new();
    let mut rest = args[1..].iter();
    while let Some(arg) = rest.next() {
        let (key, value) = if let Some((k, v)) = arg.split_once('=') {
            (k, Some(v.to_string()))
        } else {
            (arg.as_str(), rest.next().cloned())
        };
        let Some(key) = key.strip_prefix("--").filter(|k| {
            [
                "table",
                "rules",
                "genotypes",
                "phased",
                "sample",
                "out",
                "rules-out",
            ]
            .contains(k)
        }) else {
            eprintln!("{usage}\nog traits: error: unknown argument {arg}");
            return 2;
        };
        let Some(value) = value.filter(|v| !v.starts_with("--")) else {
            eprintln!("{usage}\nog traits: error: --{key} needs a value");
            return 2;
        };
        opts.insert(key.into(), value);
    }
    let needs: &[&str] = match verb {
        "check" => &["table"],
        "phase-sites" => &["rules", "genotypes", "sample", "out"],
        _ => &["table", "rules", "genotypes", "sample", "out", "rules-out"],
    };
    if let Some(k) = needs.iter().find(|k| !opts.contains_key(**k)) {
        eprintln!("{usage}\nog traits: error: missing --{k}");
        return 2;
    }
    match run(verb, &opts) {
        Ok(code) => code,
        Err(e) => {
            eprintln!("ERROR: {e}");
            1
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;
    use std::sync::atomic::{AtomicU64, Ordering};
    static SEQ: AtomicU64 = AtomicU64::new(0);
    struct Scratch(PathBuf);
    impl Scratch {
        fn new() -> Self {
            let p = std::env::temp_dir().join(format!(
                "og-traits-{}-{}",
                std::process::id(),
                SEQ.fetch_add(1, Ordering::Relaxed)
            ));
            std::fs::create_dir_all(&p).unwrap();
            Self(p)
        }
        fn put(&self, name: &str, text: &str) -> String {
            let p = self.0.join(name);
            std::fs::write(&p, text).unwrap();
            p.to_string_lossy().into_owned()
        }
    }
    impl Drop for Scratch {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }
    fn assets() -> (Vec<Site>, Vec<Rule>) {
        let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../../assets");
        (
            read_table(root.join("traits.tsv").to_str().unwrap()).unwrap(),
            read_rules(root.join("genotype_rules.tsv").to_str().unwrap()).unwrap(),
        )
    }
    fn genotypes(r: &Rule, choices: &[u8]) -> Genotypes {
        r.sites
            .iter()
            .zip(choices)
            .map(|(s, v)| {
                let al = vec![
                    s.row[if *v == 2 { "alt" } else { "ref" }].clone(),
                    s.row[if *v == 0 { "ref" } else { "alt" }].clone(),
                ];
                (
                    s.key(),
                    SiteGt {
                        alleles: al,
                        status: Status::Called,
                        gq: "40".into(),
                        weak: None,
                        snv: true,
                    },
                )
            })
            .collect()
    }
    #[test]
    fn committed_catalogs_and_every_rule_genotype() {
        let (rows, rules) = assets();
        assert!(check_rows(&rows).is_empty());
        assert_eq!(check_rules(&rules, Some(&rows)), Vec::<String>::new());
        for r in &rules {
            for mut n in 0..3usize.pow(r.sites.len() as u32) {
                let mut choices = vec![];
                for _ in &r.sites {
                    choices.push((n % 3) as u8);
                    n /= 3;
                }
                let gt = genotypes(r, &choices);
                let line = rule_line(r, &gt, None);
                assert!(RULE_OUT.iter().all(|k| line.contains_key(*k)));
            }
        }
    }
    #[test]
    fn apoe_grid_keeps_haplotype_order_and_missing_calls() {
        let (_, rules) = assets();
        let r = rules.iter().find(|r| r.id == "APOE").unwrap();
        let expected = [
            ["e3/e3", "e2/e3", "e2/e2"],
            ["e3/e4", "e2/e4", "e1/e2"],
            ["e4/e4", "e1/e4", "e1/e1"],
        ];
        for (a, labels) in expected.iter().enumerate() {
            for (b, label) in labels.iter().enumerate() {
                let line = rule_line(r, &genotypes(r, &[a as u8, b as u8]), None);
                assert_eq!(line["status"], "called");
                assert_eq!(line["result"], *label);
            }
        }
        let mut gt = genotypes(r, &[1, 0]);
        gt.remove(&r.sites[1].key());
        let line = rule_line(r, &gt, None);
        assert_eq!(line["status"], "not_covered");
        assert_eq!(line["phase"], ".");
    }
    #[test]
    fn schema_faults_are_refused_before_enumeration() {
        let (rows, rules) = assets();
        for (k, v, needle) in [
            ("ref", "N", "not one of"),
            ("alt", rows[0].row["ref"].as_str(), "same base"),
            ("pmids", "PMC123", "PubMed"),
            ("chrom", "chrX", "chr1 to chr22"),
            ("het", " ", "empty"),
            ("note", "From SNPedia.", "SNPedia"),
        ] {
            let mut bad = rows.clone();
            bad[0].row.insert(k.into(), v.into());
            assert!(check_rows(&bad).iter().any(|e| e.contains(needle)), "{k}");
        }
        let mut bad = rows.clone();
        bad.remove(bad.iter().position(|s| s.row["rsid"] == "rs7412").unwrap());
        assert!(check_rows(&bad).iter().any(|e| e.contains("missing")));
        let mut bad = rules.clone();
        bad[0].combos.pop();
        assert!(
            check_rules(&bad, Some(&rows))
                .iter()
                .any(|e| e.contains("no combination"))
        );
        let mut bad = rules.clone();
        bad[0].haps[0].alleles[0] = "N".into();
        assert!(
            check_rules(&bad, None)
                .iter()
                .any(|e| e.contains("neither"))
        );
        let mut bad = rules.clone();
        bad[0].sites.clear();
        assert!(!check_rules(&bad, None).is_empty());
        let mut bad = rules.clone();
        bad[0].sites[0].row.insert("pos".into(), "9".into());
        assert!(
            check_rules(&bad, Some(&rows))
                .iter()
                .any(|e| e.contains("traits table has"))
        );
        let mut bad = rules.clone();
        bad[0].combos[0].key = "absent/absent".into();
        assert!(
            check_rules(&bad, None)
                .iter()
                .any(|e| e.contains("not two"))
        );
        let mut bad = rules.clone();
        bad[0].combos[1].key = "e3/e2".into();
        assert!(check_rules(&bad, None).iter().any(|e| e.contains("order")));
        assert!(
            check_rules(&rules[1..], Some(&rows))
                .iter()
                .any(|e| e.contains("rule APOE is missing"))
        );
    }
    #[test]
    fn stream_prefers_first_snv_but_checks_every_ref_and_bounds_memory() {
        let (rows, _) = assets();
        let s = &rows[0];
        let scratch = Scratch::new();
        let prefix = format!("{}\t{}\t", s.row["chrom"], s.row["pos"]);
        let a = &s.row["ref"];
        let b = &s.row["alt"];
        for (gt, status) in [
            (format!("{a}|{b}"), Status::Called),
            (format!("./{b}"), Status::OtherAllele),
            (a.clone(), Status::OtherAllele),
            ("é/é".into(), Status::OtherAllele),
        ] {
            let path = scratch.put("g", &format!("{prefix}{a}\t{b}\t{gt}\t99\n"));
            let result = read_genotypes(&path, &[s]).unwrap();
            assert_eq!(result[&s.key()].status, status);
        }
        let text = format!(
            "{prefix}{a}T\t{a}\t{a}T/{a}\n{prefix}{a}\t{b}\t{b}/{a}\t42\n{prefix}{a}\t{b}\t{b}/{b}\t1\n"
        );
        let path = scratch.put("g", &text);
        let result = read_genotypes(&path, &[s]).unwrap();
        assert_eq!(result[&s.key()].alleles, vec![a.clone(), b.clone()]);
        assert_eq!(result[&s.key()].gq, "42");
        let path = scratch.put("bad", &(text + &format!("{prefix}N\t{b}\tN/N\n")));
        assert!(read_genotypes(&path, &[s]).unwrap_err().contains("GRCh38"));
        let path = scratch.put("many", &"chr99\t1\tN\t.\tN/N\n".repeat(10000));
        assert!(read_genotypes(&path, &[s]).unwrap().is_empty());
    }
    #[test]
    fn shared_position_keeps_each_rules_allele_vocabulary() {
        let (rows, _) = assets();
        let first = &rows[0];
        let mut second = first.clone();
        let alt = ["A", "C", "G", "T"]
            .into_iter()
            .find(|a| *a != first.row["ref"] && *a != first.row["alt"])
            .unwrap();
        second.row.insert("alt".into(), alt.into());
        let scratch = Scratch::new();
        let path = scratch.put(
            "g",
            &format!(
                "{}\t{}\t{}\t{}\t{}/{}\n",
                first.row["chrom"], first.row["pos"], first.row["ref"], alt, first.row["ref"], alt
            ),
        );
        let gt = read_genotypes(&path, &[first, &second]).unwrap();
        assert_eq!(genotype(first, &gt).status, Status::OtherAllele);
        assert_eq!(genotype(&second, &gt).status, Status::Called);
    }

    #[test]
    fn phase_requires_one_block_and_only_multiple_hets_are_exported() {
        let (_, rules) = assets();
        let r = rules.iter().find(|r| r.id == "TAS2R38").unwrap();
        let gt = genotypes(r, &[1, 1, 1]);
        assert_eq!(rule_line(r, &gt, None)["phase"], "unphased");
        assert!(
            rule_line(r, &gt, Some(&HashMap::new()))["reading"].contains("no read or read pair")
        );
        let scratch = Scratch::new();
        let mut vcf = String::new();
        let hap1 = &r.haps[0];
        let hap2 = &r.haps[1];
        for (i, s) in r.sites.iter().enumerate() {
            let a = usize::from(hap1.alleles[i] == s.row["alt"]);
            let b = usize::from(hap2.alleles[i] == s.row["alt"]);
            vcf += &format!(
                "{}\t{}\t.\t{}\t{}\t.\tPASS\t.\tGT:PS\t{a}|{b}:7\n",
                s.row["chrom"], s.row["pos"], s.row["ref"], s.row["alt"]
            );
        }
        let p = scratch.put("phased", &vcf);
        let mut phased = read_phased(&p, &r.sites.iter().map(Site::key).collect()).unwrap();
        assert_eq!(rule_line(r, &gt, Some(&phased))["phase"], "read_backed");
        phased.get_mut(&r.sites[0].key()).unwrap().2 = "8".into();
        assert_eq!(rule_line(r, &gt, Some(&phased))["status"], "unphased");
        let path = scratch.put("phase", "");
        assert_eq!(phase_vcf(&path, "S", &rules, &gt).unwrap(), 3);
        let text = std::fs::read_to_string(&path).unwrap();
        assert!(text.contains("141972804\t"));
        assert!(!text.contains("rs429358"));
        assert_eq!(
            phase_vcf(&path, "S", &rules, &genotypes(r, &[1, 0, 0])).unwrap(),
            0
        );
        let gt = genotypes(r, &[0, 0, 0]);
        assert_eq!(rule_line(r, &gt, None)["phase"], "not_needed");
    }
    /// A rule read without phase (HFE, SERPINA1, GC, CAD9p21, MAPT) that
    /// gets one diplotype from two heterozygous sites assumes which alleles
    /// share a copy. Unphased, the reading says so; read-backed phase in one
    /// block confirms it, or puts the alleles on copies the rule does not
    /// describe and the pair is not read.
    #[test]
    fn a_double_heterozygote_read_without_phase_says_what_it_assumes() {
        let (_, rules) = assets();
        let rule = |id: &str| rules.iter().find(|r| r.id == id).unwrap();
        // Phase of each site: 0|1 puts the alt on the second copy.
        let phased = |r: &Rule, gts: &[&str], ps: [&str; 2]| -> Phased {
            r.sites
                .iter()
                .zip(gts)
                .zip(ps)
                .map(|((s, gt), ps)| {
                    let (a, b) = gt.split_once('|').unwrap();
                    let base = |i: &str| s.row[if i == "1" { "alt" } else { "ref" }].clone();
                    (s.key(), (base(a), base(b), ps.to_owned()))
                })
                .collect()
        };
        for (id, label, place) in [
            (
                "HFE",
                "C282Y/H63D",
                "sit on different copies (the copies C282Y and H63D)",
            ),
            (
                "SERPINA1",
                "PiSZ",
                "sit on different copies (the copies S and Z)",
            ),
            (
                "GC",
                "Gc1S/Gc2",
                "sit on different copies (the copies Gc1S and Gc2)",
            ),
            (
                "CAD9p21",
                "one risk copy",
                "sit on one copy (the copies nonrisk and risk)",
            ),
            ("MAPT", "H1/H2", "sit on one copy (the copies H1 and H2)"),
        ] {
            let r = rule(id);
            let gt = genotypes(r, &[1, 1]);
            let plain = r.combos.iter().find(|c| c.label == label).unwrap();
            let line = rule_line(r, &gt, None);
            assert_eq!(
                (&line["status"][..], &line["result"][..], &line["phase"][..]),
                ("called", label, "not_used"),
                "{id}"
            );
            assert!(line["reading"].starts_with(&plain.reading), "{id}");
            assert!(
                line["reading"].contains("Phase unknown"),
                "{id}: {}",
                line["reading"]
            );
            assert!(line["reading"].contains(place), "{id}: {}", line["reading"]);
            assert!(line["reading"].contains("no phasing ran"), "{id}");
            let rs = &r.sites;
            assert!(
                line["reading"].contains(&format!(
                    "{} {} and {} {}",
                    rs[0].row["rsid"], rs[0].row["alt"], rs[1].row["rsid"], rs[1].row["alt"]
                )),
                "{id}: {}",
                line["reading"]
            );
            // Phased as the reading assumes: the plain reading, no caveat.
            let same = if place.contains("different") {
                ["1|0", "0|1"]
            } else {
                ["1|0", "1|0"]
            };
            let other = if place.contains("different") {
                ["1|0", "1|0"]
            } else {
                ["1|0", "0|1"]
            };
            let line = rule_line(r, &gt, Some(&phased(r, &same, ["5", "5"])));
            assert_eq!(
                (&line["status"][..], &line["result"][..], &line["phase"][..]),
                ("called", label, "read_backed"),
                "{id}"
            );
            assert_eq!(line["reading"], plain.reading, "{id}");
            // Two phase blocks are no phase: the caveat stays.
            let line = rule_line(r, &gt, Some(&phased(r, &same, ["5", "6"])));
            assert_eq!(line["phase"], "not_used", "{id}");
            assert!(line["reading"].contains("no read or read pair"), "{id}");
            // Phased the other way: not the assumed pair, not read.
            let line = rule_line(r, &gt, Some(&phased(r, &other, ["5", "5"])));
            assert_eq!(
                (&line["status"][..], &line["result"][..], &line["phase"][..]),
                ("other_haplotype", ".", "read_backed"),
                "{id}"
            );
            assert!(line["reading"].contains("not read"), "{id}");
            assert!(
                line["reading"].contains(&format!("would read {label}.")),
                "{id}: {}",
                line["reading"]
            );
            // One heterozygous site assumes nothing.
            for choices in [[1, 0], [0, 1], [2, 1], [1, 2]] {
                let line = rule_line(r, &genotypes(r, &choices), None);
                assert!(
                    !line["reading"].contains("Phase unknown"),
                    "{id} {choices:?}"
                );
            }
        }
        // HFE with C282Y and H63D on one copy: one copy changed, the other wt.
        let r = rule("HFE");
        let line = rule_line(
            r,
            &genotypes(r, &[1, 1]),
            Some(&phased(r, &["0|1", "0|1"], ["5", "5"])),
        );
        assert_eq!(
            line["reading"],
            "Read-backed phase shows that rs1800562 A and rs1799945 G sit on one copy: the copies \
             are wt and A,G at rs1800562, rs1799945, so one copy is changed and the other is not. \
             This rule does not describe that copy, so the pair is not read; without phase these \
             genotypes would read C282Y/H63D."
        );
        // APOE's double heterozygote line already names the other pair.
        let apoe = rule("APOE");
        let line = rule_line(apoe, &genotypes(apoe, &[1, 1]), None);
        assert_eq!(line["result"], "e2/e4");
        assert!(!line["reading"].contains("Phase unknown"));
        // WhatsHap gets the two HFE sites, not APOE's.
        let scratch = Scratch::new();
        let path = scratch.put("phase", "");
        let mut gt = genotypes(r, &[1, 1]);
        gt.extend(genotypes(apoe, &[1, 1]));
        assert_eq!(phase_vcf(&path, "S", &rules, &gt).unwrap(), 2);
        let text = std::fs::read_to_string(&path).unwrap();
        assert!(text.contains("rs1800562") && text.contains("rs1799945"));
        assert!(!text.contains("rs429358"));
        assert_eq!(
            phase_vcf(&path, "S", &rules, &genotypes(r, &[1, 2])).unwrap(),
            0
        );
    }
    /// The call floor (decision D10): a site whose GQ is below 20, or whose
    /// record the caller filtered, is not read. GQ 20, a missing GQ and
    /// DeepVariant's RefCall stay readable.
    #[test]
    fn low_quality_calls_at_rule_sites_are_not_read() {
        let (rows, rules) = assets();
        let rule = |id: &str| rules.iter().find(|r| r.id == id).unwrap();
        let sites: Vec<&Site> = rows
            .iter()
            .chain(rules.iter().flat_map(|r| &r.sites))
            .collect();
        let scratch = Scratch::new();
        let read = |text: &str, id: &str| {
            let path = scratch.put("g", text);
            rule_line(rule(id), &read_genotypes(&path, &sites).unwrap(), None)
        };
        // APOE e2/e4: heterozygous at both sites.
        let apoe = |tail: &str| {
            format!("chr19\t44908684\tT\tC\tT/C{tail}\nchr19\t44908822\tC\tT\tC/T{tail}\n")
        };
        for (tail, status, result) in [
            ("\t0", "low_gq", "."),
            ("\t1", "low_gq", "."),
            ("\t19", "low_gq", "."),
            ("\t20", "called", "e2/e4"),
            ("\t99", "called", "e2/e4"),
            ("\t.", "called", "e2/e4"),
            ("", "called", "e2/e4"),
            ("\t50\tPASS", "called", "e2/e4"),
            ("\t50\t.", "called", "e2/e4"),
            ("\t50\tLowQual", "low_gq", "."),
        ] {
            let line = read(&apoe(tail), "APOE");
            assert_eq!(line["status"], status, "APOE with '{tail}'");
            assert_eq!(line["result"], result, "APOE with '{tail}'");
            assert_eq!(line["genotypes"], "T/C,C/T", "APOE with '{tail}'");
            if status == "low_gq" {
                assert!(line["reading"].contains("not read"), "{}", line["reading"]);
                assert!(line["reading"].contains("rs429358"), "{}", line["reading"]);
            }
        }
        // One weak site is enough to stop the rule, and the reading names it.
        let line = read(
            "chr19\t44908684\tT\tC\tC/C\t1\nchr19\t44908822\tC\tT\tC/C\t50\n",
            "APOE",
        );
        assert_eq!((&line["status"][..], &line["result"][..]), ("low_gq", "."));
        assert!(line["reading"].contains("rs429358"), "{}", line["reading"]);
        assert!(!line["reading"].contains("rs7412"), "{}", line["reading"]);
        // A site with no call at all still reads not_covered.
        let line = read("chr19\t44908684\tT\tC\tT/C\t1\n", "APOE");
        assert_eq!(line["status"], "not_covered");
        // HFE C282Y/H63D: a LowQual call at C282Y is not read; RefCall and
        // a pass are.
        let hfe = |filter: &str| {
            format!(
                "chr6\t26092913\tG\tA\tG/A\t50\t{filter}\nchr6\t26090951\tC\tG\tC/G\t50\tRefCall\n"
            )
        };
        let line = read(&hfe("LowQual"), "HFE");
        assert_eq!((&line["status"][..], &line["result"][..]), ("low_gq", "."));
        assert!(line["reading"].contains("LowQual"), "{}", line["reading"]);
        assert_eq!(read(&hfe("PASS"), "HFE")["status"], "called");
        // SERPINA1 PiMZ at GQ 0 is not read; at GQ 30 it is.
        let serpina1 = |gq: u32| {
            format!("chr14\t94380925\tT\tA\tT/T\t{gq}\nchr14\t94378610\tC\tT\tC/T\t{gq}\n")
        };
        assert_eq!(read(&serpina1(0), "SERPINA1")["status"], "low_gq");
        assert_eq!(read(&serpina1(30), "SERPINA1")["status"], "called");
        // A low-quality heterozygous site is not sent to WhatsHap.
        let tas2r38 = |gq: u32| {
            format!(
                "chr7\t141973545\tC\tG\tC/G\t{gq}\nchr7\t141972905\tG\tA\tG/A\t{gq}\nchr7\t141972804\tT\tC\tT/C\t{gq}\n"
            )
        };
        for (gq, n) in [(5, 0), (40, 3)] {
            let path = scratch.put("g", &tas2r38(gq));
            let gt = read_genotypes(&path, &sites).unwrap();
            let out = scratch.put("phase.vcf", "");
            assert_eq!(phase_vcf(&out, "S", &rules, &gt).unwrap(), n, "GQ {gq}");
        }
    }

    /// The traits table: a weak table row and the APOE line built from weak
    /// sites read low_gq with a "not read" reading.
    #[test]
    fn low_quality_calls_in_the_traits_table_are_not_read() {
        let scratch = Scratch::new();
        let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../../assets");
        let gt = scratch.put(
            "g",
            "chr2\t135851076\tG\tA\tG/G\t3\n\
             chr19\t40848628\tA\tT\tA/T\t40\tLowQual\n\
             chr2\t135859184\tC\tT\tC/T\t20\n\
             chr19\t44908684\tT\tC\tC/C\t0\nchr19\t44908822\tC\tT\tC/C\t0\n",
        );
        let out = scratch.0.join("traits.tsv").to_string_lossy().into_owned();
        let rules_out = scratch.0.join("rules.tsv").to_string_lossy().into_owned();
        let opts = [
            ("table", root.join("traits.tsv")),
            ("rules", root.join("genotype_rules.tsv")),
        ]
        .into_iter()
        .map(|(k, v)| (k.to_string(), v.to_string_lossy().into_owned()))
        .chain([
            ("genotypes".into(), gt),
            ("sample".into(), "S".into()),
            ("out".into(), out.clone()),
            ("rules-out".into(), rules_out),
        ])
        .collect();
        assert_eq!(run("table", &opts).unwrap(), 0);
        let text = std::fs::read_to_string(&out).unwrap();
        let row = |rsid: &str| -> Vec<String> {
            text.lines()
                .map(|l| l.split('\t').map(str::to_owned).collect::<Vec<_>>())
                .find(|c| c[0] == rsid)
                .unwrap_or_else(|| panic!("no row {rsid}"))
        };
        // rsid gene chrom pos trait genotype status reading ...
        for rsid in ["rs4988235", "rs1801272", "rs429358", "rs7412"] {
            let r = row(rsid);
            assert_eq!(r[6], "low_gq", "{rsid}: {r:?}");
            assert!(r[7].contains("not read"), "{rsid}: {r:?}");
        }
        assert_eq!(row("rs4988235")[5], "G/G");
        assert!(row("rs1801272")[7].contains("LowQual"));
        assert_eq!(row("rs182549")[6], "called");
        let apoe = row("rs429358+rs7412");
        assert_eq!((&apoe[5][..], &apoe[6][..]), (".", "low_gq"));
    }

    #[test]
    fn unwritable_rules_output_keeps_existing_traits_table() {
        let scratch = Scratch::new();
        let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../../assets");
        let target = scratch.put("out", "previous complete table\n");
        let gt = scratch.put("g", "");
        let opts = [
            (
                "table",
                root.join("traits.tsv").to_string_lossy().into_owned(),
            ),
            (
                "rules",
                root.join("genotype_rules.tsv")
                    .to_string_lossy()
                    .into_owned(),
            ),
            ("genotypes", gt),
            ("sample", "S".into()),
            ("out", target.clone()),
            (
                "rules-out",
                scratch
                    .0
                    .join("missing/rules.tsv")
                    .to_string_lossy()
                    .into_owned(),
            ),
        ]
        .into_iter()
        .map(|(k, v)| (k.into(), v))
        .collect();
        assert!(run("table", &opts).unwrap_err().contains("cannot write"));
        assert_eq!(
            std::fs::read_to_string(target).unwrap(),
            "previous complete table\n"
        );
    }

    #[test]
    fn catalog_quotes_are_literal_and_embedded_tabs_are_columns() {
        let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../../assets");
        let text = std::fs::read_to_string(root.join("traits.tsv")).unwrap();
        let rules = std::fs::read_to_string(root.join("genotype_rules.tsv")).unwrap();
        let scratch = Scratch::new();
        let path = scratch.put(
            "quoted-field",
            &text.replacen("\tMCM6\t", "\t\"MCM6\"\t", 1),
        );
        assert!(
            check_rows(&read_table(&path).unwrap())
                .iter()
                .any(|e| e.contains("quote"))
        );
        let path = scratch.put(
            "quoted-rule",
            &rules.replacen("\tAPOE\tAPOE\t", "\tAPOE\t\"APOE\"\t", 1),
        );
        assert!(!check_rules(&read_rules(&path).unwrap(), None).is_empty());
        let path = scratch.put(
            "embedded-tab",
            &text.replacen("\tMCM6\t", "\t\"MCM6\textra\"\t", 1),
        );
        assert!(
            read_table(&path)
                .unwrap_err()
                .contains("14 columns, expected 13")
        );
    }

    #[test]
    fn directory_destinations_keep_both_previous_tables() {
        let scratch = Scratch::new();
        let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../../assets");
        let gt = scratch.put("g", "");
        for directory_key in ["out", "rules-out"] {
            let previous = scratch.put("previous.tsv", "previous complete table\n");
            let mut opts: HashMap<String, String> = [
                (
                    "table",
                    root.join("traits.tsv").to_string_lossy().into_owned(),
                ),
                (
                    "rules",
                    root.join("genotype_rules.tsv")
                        .to_string_lossy()
                        .into_owned(),
                ),
                ("genotypes", gt.clone()),
                ("sample", "S".into()),
                ("out", previous.clone()),
                ("rules-out", previous.clone()),
            ]
            .into_iter()
            .map(|(k, v)| (k.into(), v))
            .collect();
            opts.insert(
                directory_key.into(),
                scratch.0.to_string_lossy().into_owned(),
            );
            assert!(run("table", &opts).unwrap_err().contains("cannot write"));
            assert!(
                std::fs::read_to_string(&previous).unwrap() == "previous complete table\n",
                "{directory_key}: previous output changed"
            );
            assert_eq!(std::fs::read_dir(&scratch.0).unwrap().count(), 2);
        }
    }

    #[test]
    fn plain_header_unicode_rules_and_atomic_write_errors() {
        let scratch = Scratch::new();
        let mut text = COLUMNS
            .iter()
            .map(|k| format!("\"{k}\""))
            .collect::<Vec<_>>()
            .join("\t");
        text.push_str("\r\n");
        let path = scratch.put("quoted", &text);
        assert!(
            read_table(&path)
                .unwrap_err()
                .contains("unexpected traits header")
        );
        let path = scratch.put("r", "rule\té\tG\tT\tno\tyes\tN\t1\n");
        assert!(!check_rules(&read_rules(&path).unwrap(), None).is_empty());
        assert!(stage_tsv(scratch.0.join("absent/x").to_str().unwrap(), OUT, &[]).is_err());
        assert_eq!(cli(&["table".into()]), 2);
        assert_eq!(
            cli(&[
                "check".into(),
                "--table".into(),
                scratch.0.join("absent").to_string_lossy().into_owned()
            ]),
            1
        );
    }
}
