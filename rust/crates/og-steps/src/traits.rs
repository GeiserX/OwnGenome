//! Step 42: bounded site genotypes, haplotype rules and read-backed phase.
use indexmap::IndexMap;
use og_io::{AtomicFile, Lines, Utf8};
use pycompat::csv::{Quoting, Reader, Writer};
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
#[derive(Clone, Debug, PartialEq, Eq)]
enum Status {
    Called,
    NotCovered,
    OtherAllele,
    Unphased,
    OtherHaplotype,
}
impl Status {
    fn text(&self) -> &'static str {
        match self {
            Self::Called => "called",
            Self::NotCovered => "not_covered",
            Self::OtherAllele => "other_allele",
            Self::Unphased => "unphased",
            Self::OtherHaplotype => "other_haplotype",
        }
    }
}
#[derive(Clone, Debug)]
struct SiteGt {
    alleles: Vec<String>,
    status: Status,
    gq: String,
    snv: bool,
}
impl Default for SiteGt {
    fn default() -> Self {
        Self {
            alleles: vec![],
            status: Status::NotCovered,
            gq: ".".into(),
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
    // Parse even the header as a CSV record. Comments do not enter the parser.
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
        let cells = Reader::new(
            std::iter::once(Ok::<_, std::io::Error>(text.clone())),
            '\t',
            Quoting::Minimal,
        )
        .next()
        .ok_or("empty record")?
        .map_err(|e| e.to_string())?;
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
    gt.status = if gt.alleles.len() == 2
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
                "{path}:{n}: {} columns, expected CHROM POS REF ALT GENOTYPE [GQ]",
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
        let status = if al.len() == 2
            && al
                .iter()
                .all(|a| a == &site.row["ref"] || a == &site.row["alt"])
        {
            Status::Called
        } else {
            Status::OtherAllele
        };
        let gq = c
            .get(5)
            .filter(|g| !g.is_empty())
            .unwrap_or(&".")
            .to_string();
        result.insert(
            key,
            SiteGt {
                alleles: al,
                status,
                gq,
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
        let status = if g.iter().any(|g| g.status == Status::NotCovered) {
            Status::NotCovered
        } else {
            Status::OtherAllele
        };
        return finish(
            out,
            ".",
            ".",
            status,
            format!("Not read: needs a genotype at {need}."),
        );
    }
    let mut phase = if r.meta["needs_phase"] == "no" {
        "not_used"
    } else {
        "not_needed"
    };
    let key = if r.meta["needs_phase"] == "no" {
        allowed(r, &g.iter().map(|g| g.alleles.clone()).collect::<Vec<_>>()).join("|")
    } else {
        let het: Vec<_> = g
            .iter()
            .enumerate()
            .filter(|(_, g)| g.alleles[0] != g.alleles[1])
            .map(|(i, _)| i)
            .collect();
        let mut h1: Vec<_> = g.iter().map(|g| g.alleles[0].clone()).collect();
        let mut h2: Vec<_> = g.iter().map(|g| g.alleles[1].clone()).collect();
        if het.len() >= 2 {
            let calls: Option<Vec<_>> = het
                .iter()
                .map(|i| phased.and_then(|p| p.get(&r.sites[*i].key())))
                .collect();
            if let Some(calls) = calls.filter(|c| c.iter().all(|v| v.2 == c[0].2)) {
                for (i, c) in het.iter().zip(calls) {
                    h1[*i] = c.0.clone();
                    h2[*i] = c.1.clone();
                }
                phase = "read_backed";
            } else {
                let why = if phased.is_some() {
                    "no read or read pair covered them together"
                } else {
                    "no phasing ran (no BAM, or phasing failed)"
                };
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
        Some(c) => finish(out, phase, &c.label, Status::Called, c.reading.clone()),
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
    {
        let mut writer = Writer::new(&mut file, '\t');
        writer
            .write_row(columns.iter().copied())
            .map_err(|e| format!("cannot write {path}: {e}"))?;
        for row in rows {
            writer
                .write_row(columns.iter().map(|k| row[*k].as_str()))
                .map_err(|e| format!("cannot write {path}: {e}"))?;
        }
    }
    file.flush()
        .map_err(|e| format!("cannot write {path}: {e}"))?;
    Ok(file)
}
fn phase_vcf(path: &str, sample: &str, rules: &[Rule], gt: &Genotypes) -> Result<usize> {
    let mut sites = BTreeMap::new();
    for r in rules.iter().filter(|r| r.meta["needs_phase"] == "yes") {
        let calls: Option<Vec<_>> = r
            .sites
            .iter()
            .map(|s| {
                let g = genotype(s, gt);
                (g.status == Status::Called).then_some(g)
            })
            .collect();
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
    fn csv_header_unicode_rules_and_atomic_write_errors() {
        let scratch = Scratch::new();
        let mut text = COLUMNS
            .iter()
            .map(|k| format!("\"{k}\""))
            .collect::<Vec<_>>()
            .join("\t");
        text.push_str("\r\n");
        let path = scratch.put("quoted", &text);
        assert!(read_table(&path).unwrap().is_empty());
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
