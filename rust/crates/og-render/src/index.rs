//! `og report index`: the search index of one run, every gene symbol and
//! rsID its results name, and the query rule the live site's search box
//! (`assets/live/plugins/search.js`) applies to it in the browser.
//!
//! ```text
//! og report index build --summary summary.json [--sample-dir DIR] -o INDEX.json
//! og report index query --index INDEX.json QUERY
//! ```
//!
//! What it reads:
//! - every section of summary.json in state ok or stale: each dict in its
//!   data that names a gene or an rsID, by `og lookup`'s names rule
//!   ([`og_lookup::names::names_of`]), so the search and a lookup cannot
//!   disagree. A dict inside one that already named a name adds only the
//!   names its parents did not;
//! - the gene column of step 23 (clinical filter) and step 31 (slivar),
//!   whose sections keep counts only: from the summary table on the bash
//!   path (`GENE`; `SYMBOL` and `Existing_variation`), else from the VCF the
//!   Nextflow report stages (INFO/PGP_GENE, else the genes of the CSQ
//!   entries `csq_rank` picks; the rsIDs of every CSQ entry).
//!
//! A section that is missing, failed or unreadable gives no entry, and a
//! table that cannot be read gives none either: one warning line on stderr,
//! and the index is still written.
//!
//! The index:
//! `{schema_version: 1, sample, generated_utc, terms: {TERM: [{section,
//! anchor, title, step, state, source, count, mentions: [{text, locus}],
//! veil}]}}`. TERM is a gene symbol in upper case or an rsID in lower case.
//! One item per section that names the term; count is how many variants (a
//! mention with a locus counts once per locus) or other mentions it has
//! there, mentions the first [`MENTIONS_KEPT`] of them. anchor is the
//! section's id in the HTML report, or null for a key no section shows.
//! veil, only on an item that names a result a child's page keeps closed
//! (its gene, or a locus inside it), is that result's id ([`Closable`]).
//!
//! A query ([`query`], and the same rule in search.js, which a cargo test
//! runs under node) is an rsID, a gene symbol (a prefix also lists the terms
//! it starts: BRCA gives BRCA1 and BRCA2) or a position, which the index
//! cannot answer: a position is a genotype question, for
//! `scripts/owngenome-lookup.sh`.
//!
//! The veil needs the findings a child's page keeps closed, which the live
//! results build. Until they are ported, the command line passes none, so
//! `og report index build` writes no veil field; [`build`] takes them.

use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::io::Write;

use indexmap::IndexMap;
use og_lookup::names::{flat, is_rsid, names_of};
use og_lookup::{Query, parse_query};
use og_summary::collect::acmg::csq_rank;
use og_summary::collect::text::{ReadError, cell, cell_or, each_line, info_map, rows};
use pycompat::json::{Dump, Value, dumps, loads};
use pycompat::text;
use pycompat::value::{py_str, truthy};

use crate::chapters::section_of;

pub const SCHEMA_VERSION: i128 = 1;
/// Mentions kept per term and section; count says how many there were.
pub const MENTIONS_KEPT: usize = 20;
/// Terms listed for a prefix query.
pub const SUGGESTIONS: usize = 20;
/// What a position query gets instead of items.
const POSITION_HINT: &str = "The search finds gene symbols and rsIDs in these results. A position is a \
question about your genotype: run scripts/owngenome-lookup.sh <sample> {q}, which reads it from the gVCF on \
this computer.";

const USAGE: &str = "usage: og report index build --summary JSON [--sample-dir DIR] -o INDEX.json
       og report index query --index INDEX.json QUERY";

/// A result a child's page keeps closed: its id, its gene (upper case) and
/// its locus (`chrN`, start, end) when it has one.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Closable {
    pub id: String,
    pub gene: Option<String>,
    pub locus: Option<(String, i128, i128)>,
}

/// One section's entry for one term while the index is built.
struct Item {
    section: String,
    anchor: Value,
    title: Value,
    step: Value,
    state: Value,
    source: Value,
    count: i128,
    mentions: Vec<(String, Option<String>)>,
    /// One mention per variant (locus) or per text: the ClinVar list, the
    /// secondary findings and a step's table can name the same variant.
    seen: HashSet<(bool, String)>,
}

/// A mention of a step's table: its text, its locus and the names in it.
type Record = (String, String, BTreeSet<String>, BTreeSet<String>);

/// Every dict of `obj` that names a gene or rsID no dict above it named,
/// with those new names.
fn walk<'a>(
    obj: &'a Value,
    blocked: &BTreeSet<String>,
    out: &mut Vec<(&'a Value, BTreeSet<String>)>,
) {
    match obj {
        Value::Object(o) => {
            let (genes, rsids) = names_of(obj);
            let names: BTreeSet<String> = genes.into_iter().chain(rsids).collect();
            let new: BTreeSet<String> = names.difference(blocked).cloned().collect();
            if !new.is_empty() {
                out.push((obj, new));
            }
            let inner: BTreeSet<String> = blocked.union(&names).cloned().collect();
            for v in o.values() {
                walk(v, &inner, out);
            }
        }
        Value::Array(a) => {
            for v in a {
                walk(v, blocked, out);
            }
        }
        _ => {}
    }
}

/// `true` for a character of Python's `\w`.
fn is_word(c: char) -> bool {
    c.is_alphanumeric() || c == '_'
}

/// `chrom:pos` of a mention that carries them (`chrom` and `pos` keys, else
/// a `variant` that starts with `chrom:pos`), else `None`.
fn locus_of(d: &Value) -> Option<String> {
    let Value::Object(o) = d else {
        return None;
    };
    if let Some(Value::Str(c)) = o.get("chrom")
        && !c.is_empty()
    {
        match o.get("pos") {
            Some(Value::Int(p)) => return Some(format!("{c}:{p}")),
            Some(Value::Str(p)) if !p.is_empty() && p.chars().all(text::is_decimal) => {
                return Some(format!("{c}:{p}"));
            }
            _ => {}
        }
    }
    let variant = match o.get("variant") {
        Some(v) if truthy(v) => py_str(v).unwrap_or_default(),
        _ => return None,
    };
    // ^((?:chr)?[0-9XYMT]{1,2}):([0-9]+)\b
    let rest = variant.strip_prefix("chr").unwrap_or(&variant);
    let chrom_len = rest
        .bytes()
        .take(2)
        .take_while(|b| b.is_ascii_digit() || b"XYMT".contains(b))
        .count();
    if chrom_len == 0 || rest.as_bytes().get(chrom_len) != Some(&b':') {
        return None;
    }
    let after = &rest[chrom_len + 1..];
    let digits = after.bytes().take_while(u8::is_ascii_digit).count();
    if digits == 0 || after[digits..].chars().next().is_some_and(is_word) {
        return None;
    }
    let head = &variant[..variant.len() - rest.len() + chrom_len];
    Some(format!("{head}:{}", &after[..digits]))
}

/// The rsIDs of an Existing_variation value (`rs123&COSV1` gives `rs123`).
fn rs_ids(value: &str) -> BTreeSet<String> {
    value
        .split(['&', ','])
        .filter(|x| is_rsid(x))
        .map(text::lower)
        .collect()
}

fn variant_text(chrom: &str, pos: &str, r: &str, a: &str) -> String {
    format!("{chrom}:{pos} {r}>{a}")
}

/// `{x.upper()} - {"", "."}`.
fn one_gene(g: &str) -> BTreeSet<String> {
    let g = text::upper(g);
    if g.is_empty() || g == "." {
        BTreeSet::new()
    } else {
        BTreeSet::from([g])
    }
}

/// The variants of step 23 or 31's summary table: the gene column (`GENE`
/// or `SYMBOL`), and for step 31 the rsIDs of `Existing_variation`.
fn table_records(path: &str, gene_col: &str, with_rsids: bool) -> Result<Vec<Record>, ReadError> {
    Ok(rows(path)?
        .iter()
        .map(|r| {
            let c = |k: &str| cell(r, k).unwrap_or("None").to_owned();
            let (chrom, pos) = (c("CHROM"), c("POS"));
            let rsids = if with_rsids {
                rs_ids(cell_or(r, "Existing_variation"))
            } else {
                BTreeSet::new()
            };
            (
                variant_text(&chrom, &pos, &c("REF"), &c("ALT")),
                format!("{chrom}:{pos}"),
                one_gene(cell_or(r, gene_col)),
                rsids,
            )
        })
        .collect())
}

/// The variants of step 23 or 31's VCF: the gene of `gene_tag` when the
/// record has one, else the genes of the CSQ entries `csq_rank` picks; with
/// `with_rsids`, the rsIDs of every CSQ entry's Existing_variation.
fn vcf_records(
    path: &str,
    gene_tag: Option<&str>,
    with_rsids: bool,
) -> Result<Vec<Record>, ReadError> {
    let mut fmt: Option<Vec<String>> = None;
    let mut out = Vec::new();
    each_line(path, |line| {
        if line.starts_with("##INFO=<ID=CSQ") {
            fmt = line.split_once("Format: ").and_then(|(_, f)| {
                let f: String = f.chars().take_while(|&c| c != '"').collect();
                (!f.is_empty()).then(|| f.split('|').map(str::to_owned).collect())
            });
            return;
        }
        if line.starts_with('#') {
            return;
        }
        let r: Vec<&str> = line.split('\t').collect();
        if r.len() < 8 {
            return;
        }
        let info = info_map(r[7]);
        let csq = info.get("CSQ").copied().unwrap_or("");
        let tagged = gene_tag
            .and_then(|t| info.get(t))
            .filter(|g| !g.is_empty() && **g != ".");
        let genes = match (tagged, &fmt) {
            (Some(g), _) => BTreeSet::from([text::upper(g)]),
            (None, Some(f)) => csq_rank(f, csq)
                .1
                .iter()
                .filter(|e| !e.gene.is_empty())
                .map(|e| text::upper(&e.gene))
                .collect(),
            (None, None) => BTreeSet::new(),
        };
        let mut rsids = BTreeSet::new();
        if with_rsids
            && let Some(f) = &fmt
            && let Some(i) = f.iter().position(|k| k == "Existing_variation")
        {
            for t in csq.split(',') {
                let v: Vec<&str> = t.split('|').collect();
                if v.len() == f.len() {
                    rsids.extend(rs_ids(v[i]));
                }
            }
        }
        out.push((
            variant_text(r[0], r[1], r[3], r[4]),
            format!("{}:{}", r[0], r[1]),
            genes,
            rsids,
        ));
    })?;
    Ok(out)
}

/// The variants of a section whose genes are in its step's table, not in
/// summary.json: step 23 (`clinical`) and step 31 (`slivar`). `None` for
/// any other section.
fn table_of(key: &str, path: &str) -> Option<Result<Vec<Record>, ReadError>> {
    let tsv = path.ends_with(".tsv");
    Some(match (key, tsv) {
        ("clinical", true) => table_records(path, "GENE", false),
        ("clinical", false) => vcf_records(path, Some("PGP_GENE"), false),
        ("slivar", true) => table_records(path, "SYMBOL", true),
        ("slivar", false) => vcf_records(path, None, true),
        _ => return None,
    })
}

/// `os.path.isabs(src)` and `os.path.join(sample_dir, src)`.
fn table_path(sample_dir: &str, src: &str) -> String {
    if src.starts_with('/') || sample_dir.is_empty() {
        src.to_owned()
    } else if sample_dir.ends_with('/') {
        format!("{sample_dir}{src}")
    } else {
        format!("{sample_dir}/{src}")
    }
}

/// Every `chrom:pos` or `chrom:start-end` in `text`, as (`chrN`, low,
/// high): a chromosome `1`-`99`, X, Y, M or MT, with an optional `chr`,
/// not right after a letter or digit; commas in the numbers are dropped.
fn loci(s: &str) -> Vec<(String, i128, i128)> {
    let b = s.as_bytes();
    let mut out = Vec::new();
    let mut i = 0;
    while i < b.len() {
        let free = i == 0 || !b[i - 1].is_ascii_alphanumeric();
        match free.then(|| locus_at(b, i)).flatten() {
            Some((end, found)) => {
                out.push(found);
                i = end;
            }
            None => i += 1,
        }
    }
    out
}

/// The number of a `[0-9][0-9,]*` run at `b[i..]`, and where it ends.
fn number_at(b: &[u8], i: usize) -> Option<(usize, i128)> {
    if !b.get(i).is_some_and(u8::is_ascii_digit) {
        return None;
    }
    let end = i + b[i..]
        .iter()
        .take_while(|c| c.is_ascii_digit() || **c == b',')
        .count();
    let n = b[i..end]
        .iter()
        .filter(|c| c.is_ascii_digit())
        .fold(0i128, |n, c| {
            n.saturating_mul(10).saturating_add(i128::from(c - b'0'))
        });
    Some((end, n))
}

/// The locus that starts at `b[i]`, if one does, and where it ends.
fn locus_at(b: &[u8], i: usize) -> Option<(usize, (String, i128, i128))> {
    let starts = if b[i..].starts_with(b"chr") {
        vec![i + 3, i]
    } else {
        vec![i]
    };
    for j in starts {
        let two = b.get(j..j + 2).filter(|d| d.iter().all(u8::is_ascii_digit));
        let one = b.get(j..j + 1).filter(|d| d[0].is_ascii_digit());
        let names: [Option<&[u8]>; 6] = [two, one, Some(b"X"), Some(b"Y"), Some(b"M"), Some(b"MT")];
        for name in names.into_iter().flatten() {
            let k = j + name.len();
            if !b[j..].starts_with(name) || b.get(k) != Some(&b':') {
                continue;
            }
            let Some((end, lo)) = number_at(b, k + 1) else {
                continue;
            };
            let (end, hi) = match b.get(end) {
                Some(b'-') => number_at(b, end + 1).unwrap_or((end, lo)),
                _ => (end, lo),
            };
            let chrom = format!("chr{}", String::from_utf8_lossy(name));
            return Some((end, (chrom, lo.min(hi), lo.max(hi))));
        }
    }
    None
}

/// The words of `text` in upper case (`[A-Z0-9][A-Z0-9._-]*`), each also
/// without its trailing `.`, `_` or `-`: "in KCNH2." names KCNH2, and HLA-A
/// stays whole.
fn words(s: &str) -> HashSet<String> {
    let up = text::upper(s);
    let mut out = HashSet::new();
    for w in up
        .split(|c: char| {
            !(c.is_ascii_uppercase() || c.is_ascii_digit() || matches!(c, '.' | '_' | '-'))
        })
        .filter(|w| !w.is_empty())
    {
        // A word starts at its first letter or digit.
        let start = w.find(|c: char| c.is_ascii_uppercase() || c.is_ascii_digit());
        if let Some(st) = start {
            let w = &w[st..];
            out.insert(w.to_owned());
            out.insert(w.trim_end_matches(['.', '_', '-']).to_owned());
        }
    }
    out
}

/// The id of the closed result a row belongs to, or `None`: the result
/// whose locus the row's text names exactly comes first, else the first
/// result whose gene the text names or whose span holds a locus of it.
pub fn veil_of(spans: &[Closable], text: &str) -> Option<String> {
    if spans.is_empty() {
        return None;
    }
    let locs = loci(text);
    if let Some(c) = spans
        .iter()
        .find(|c| c.locus.as_ref().is_some_and(|l| locs.contains(l)))
    {
        return Some(c.id.clone());
    }
    let ws = words(text);
    spans
        .iter()
        .find(|c| {
            c.gene
                .as_ref()
                .is_some_and(|g| !g.is_empty() && ws.contains(g))
                || c.locus.as_ref().is_some_and(|(lc, a, b)| {
                    locs.iter().any(|(c2, s, e)| c2 == lc && s <= b && e >= a)
                })
        })
        .map(|c| c.id.clone())
}

/// The index of `summary` (a summary.json value). `sample_dir`, when given,
/// is where the step 23 and 31 tables are read from (the paths summary.json
/// records); `spans` are the results a child's page keeps closed; a table
/// that cannot be read is said on `warn`.
pub fn build(
    summary: &Value,
    sample_dir: Option<&str>,
    spans: &[Closable],
    warn: &mut dyn Write,
) -> Result<Value, String> {
    let Value::Object(top) = summary else {
        return Err("summary.json is not an object".to_owned());
    };
    let sections = match top.get("sections") {
        Some(Value::Object(m)) => Some(m),
        Some(v) if truthy(v) => return Err("summary.sections is not an object".to_owned()),
        _ => None,
    };
    let mut terms: BTreeMap<String, Vec<Item>> = BTreeMap::new();
    for (key, sec) in sections.into_iter().flatten() {
        let Value::Object(s) = sec else {
            return Err(format!("summary.sections.{key} is not an object"));
        };
        if !matches!(s.get("state"), Some(Value::Str(st)) if st == "ok" || st == "stale") {
            continue;
        }
        // The ACMG card of its own section is the card the clinical section
        // also holds (`clinical.data.acmg_sf`), which the index already reads.
        // Skip it so a gene is not listed twice, until the report readers
        // move from the nested card to this section (i411-display-sf).
        if key == "acmg_sf" {
            continue;
        }
        let field = |k: &str| s.get(k).cloned().unwrap_or(Value::Null);
        let mut add = |term: &str, txt: &str, locus: Option<&str>| {
            let items = terms.entry(term.to_owned()).or_default();
            let at = match items.iter().position(|x| x.section == *key) {
                Some(at) => at,
                None => {
                    items.push(Item {
                        section: key.clone(),
                        anchor: section_of(key)
                            .map_or(Value::Null, |x| Value::Str(x.id.to_owned())),
                        title: s
                            .get("title")
                            .filter(|t| truthy(t))
                            .cloned()
                            .unwrap_or_else(|| Value::Str(key.clone())),
                        step: field("step"),
                        state: field("state"),
                        source: field("source"),
                        count: 0,
                        mentions: Vec::new(),
                        seen: HashSet::new(),
                    });
                    items.len() - 1
                }
            };
            let item = &mut items[at];
            let mark = match locus {
                Some(l) => (true, l.to_owned()),
                None => (false, txt.to_owned()),
            };
            if !item.seen.insert(mark) {
                return;
            }
            item.count += 1;
            if item.mentions.len() < MENTIONS_KEPT {
                item.mentions
                    .push((txt.to_owned(), locus.map(str::to_owned)));
            }
        };
        let mut found = Vec::new();
        if let Some(d) = s.get("data").filter(|d| truthy(d)) {
            walk(d, &BTreeSet::new(), &mut found);
        }
        for (d, names) in found {
            let (txt, locus) = (flat(d, 160), locus_of(d));
            for n in &names {
                add(n, &txt, locus.as_deref());
            }
        }
        let src = match s.get("source") {
            Some(Value::Str(x)) if !x.is_empty() => x,
            _ => continue,
        };
        let Some(dir) = sample_dir else { continue };
        let path = table_path(dir, src);
        if !std::path::Path::new(&path).is_file() {
            continue;
        }
        let records = match table_of(key, &path) {
            None => continue,
            Some(Ok(r)) => r,
            Some(Err(e)) => {
                let _ = writeln!(
                    warn,
                    "[WARN] {key}: cannot read {src} for the search index: {}",
                    e.text()
                );
                Vec::new()
            }
        };
        for (txt, locus, genes, rsids) in records {
            for n in genes.union(&rsids) {
                add(n, &txt, Some(&locus));
            }
        }
    }
    let mut out = IndexMap::new();
    for (term, items) in terms {
        let list = items
            .into_iter()
            .map(|it| {
                // The veil of each item, from the term and that item's own
                // mentions: a locus another section names does not veil it.
                let shown: Vec<String> = std::iter::once(term.clone())
                    .chain(
                        it.mentions
                            .iter()
                            .map(|(t, l)| format!("{t} {}", l.as_deref().unwrap_or(""))),
                    )
                    .collect();
                let veil = veil_of(spans, &shown.join(" "));
                let mut o = IndexMap::new();
                o.insert("section".to_owned(), Value::Str(it.section));
                o.insert("anchor".to_owned(), it.anchor);
                o.insert("title".to_owned(), it.title);
                o.insert("step".to_owned(), it.step);
                o.insert("state".to_owned(), it.state);
                o.insert("source".to_owned(), it.source);
                o.insert("count".to_owned(), Value::Int(it.count));
                let mentions = it
                    .mentions
                    .into_iter()
                    .map(|(t, l)| {
                        let mut m = IndexMap::new();
                        m.insert("text".to_owned(), Value::Str(t));
                        m.insert("locus".to_owned(), l.map_or(Value::Null, Value::Str));
                        Value::Object(m)
                    })
                    .collect();
                o.insert("mentions".to_owned(), Value::Array(mentions));
                if let Some(id) = veil {
                    o.insert("veil".to_owned(), Value::Str(id));
                }
                Value::Object(o)
            })
            .collect();
        out.insert(term, Value::Array(list));
    }
    let mut idx = IndexMap::new();
    idx.insert("schema_version".to_owned(), Value::Int(SCHEMA_VERSION));
    idx.insert(
        "sample".to_owned(),
        top.get("sample").cloned().unwrap_or(Value::Null),
    );
    idx.insert(
        "generated_utc".to_owned(),
        top.get("generated_utc").cloned().unwrap_or(Value::Null),
    );
    idx.insert("terms".to_owned(), Value::Object(out));
    Ok(Value::Object(idx))
}

/// What a query asks.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kind {
    Rsid,
    Gene,
    Site,
    Error,
}

impl Kind {
    /// The word search.js uses for it.
    pub fn word(self) -> &'static str {
        match self {
            Kind::Rsid => "rsid",
            Kind::Gene => "gene",
            Kind::Site => "site",
            Kind::Error => "error",
        }
    }
}

/// The answer to one query.
#[derive(Debug, Clone, PartialEq)]
pub struct Answer {
    pub kind: Kind,
    /// The term looked up: an rsID in lower case, a gene in upper case.
    pub term: Option<String>,
    /// The term's items in the index.
    pub items: Vec<Value>,
    /// For a gene: up to [`SUGGESTIONS`] other terms that start with it, in
    /// order.
    pub suggestions: Vec<String>,
    /// For a position or a query that is none of these: what to do instead.
    pub hint: Option<String>,
}

/// The answer to query `q` on `index`: an rsID, a gene symbol (whose
/// suggestions are the other terms that start with it) or a position,
/// which the index cannot answer. search.js's `match()` is the same rule.
pub fn query(index: &Value, q: &str) -> Answer {
    let none = |kind, hint| Answer {
        kind,
        term: None,
        items: Vec::new(),
        suggestions: Vec::new(),
        hint: Some(hint),
    };
    let parsed = match parse_query(q) {
        Ok(p) => p,
        Err(e) => return none(Kind::Error, e.msg),
    };
    let terms = match index {
        Value::Object(o) => match o.get("terms") {
            Some(Value::Object(t)) => Some(t),
            _ => None,
        },
        _ => None,
    };
    let (kind, term) = match parsed {
        Query::Site(..) => {
            return none(
                Kind::Site,
                POSITION_HINT.replace("{q}", text::strip(q, None)),
            );
        }
        Query::Rsid(r) => (Kind::Rsid, text::lower(&r)),
        Query::Gene(g) => (Kind::Gene, text::upper(&g)),
    };
    let items = match terms.and_then(|t| t.get(&term)) {
        Some(Value::Array(a)) => a.clone(),
        _ => Vec::new(),
    };
    let mut suggestions = Vec::new();
    if kind == Kind::Gene {
        let mut keys: Vec<&String> = terms.into_iter().flat_map(|t| t.keys()).collect();
        keys.sort();
        suggestions = keys
            .into_iter()
            .filter(|k| **k != term && k.starts_with(term.as_str()))
            .take(SUGGESTIONS)
            .cloned()
            .collect();
    }
    Answer {
        kind,
        term: Some(term),
        items,
        suggestions,
        hint: None,
    }
}

/// `str(v)` of a field of an item, `None` when it is missing.
fn field(item: &Value, key: &str) -> String {
    match item {
        Value::Object(o) => o
            .get(key)
            .map_or_else(|| "None".to_owned(), |v| py_str(v).unwrap_or_default()),
        _ => "None".to_owned(),
    }
}

/// The lines `og report index query` prints for an answer.
pub fn query_lines(a: &Answer) -> Vec<String> {
    if matches!(a.kind, Kind::Error | Kind::Site) {
        return vec![a.hint.clone().unwrap_or_default()];
    }
    let term = a.term.as_deref().unwrap_or_default();
    let mut out = Vec::new();
    if a.items.is_empty() {
        out.push(format!("none of these results mentions {term}"));
    } else {
        out.push(format!("{term} is in {} section(s):", a.items.len()));
        for it in &a.items {
            let stale = if field(it, "state") == "stale" {
                " (stale: from an earlier run)"
            } else {
                ""
            };
            out.push(format!(
                "  {} (step {}){stale}: {} mention(s)",
                field(it, "title"),
                field(it, "step"),
                field(it, "count")
            ));
            let mentions = match it {
                Value::Object(o) => match o.get("mentions") {
                    Some(Value::Array(m)) => m.as_slice(),
                    _ => &[],
                },
                _ => &[],
            };
            for m in mentions {
                out.push(format!("    - {}", field(m, "text")));
            }
            let count = match it {
                Value::Object(o) => match o.get("count") {
                    Some(Value::Int(n)) => *n,
                    _ => 0,
                },
                _ => 0,
            };
            let shown = i128::try_from(mentions.len()).unwrap_or(i128::MAX);
            if count > shown {
                out.push(format!(
                    "    - ... and {} more in {}",
                    count - shown,
                    field(it, "source")
                ));
            }
        }
    }
    if !a.suggestions.is_empty() {
        out.push(format!(
            "terms that start with it: {}",
            a.suggestions.join(", ")
        ));
    }
    out
}

/// The options of one command: `--name value` pairs and positionals.
fn options(
    args: &[String],
    names: &[&str],
) -> Result<(BTreeMap<String, String>, Vec<String>), String> {
    let mut opts = BTreeMap::new();
    let mut pos = Vec::new();
    let mut it = args.iter();
    while let Some(arg) = it.next() {
        if arg == "--" {
            pos.extend(it.by_ref().cloned());
            break;
        }
        if !arg.starts_with('-') || arg == "-" {
            pos.push(arg.clone());
            continue;
        }
        let (name, inline) = match arg.split_once('=') {
            Some((n, v)) if n.starts_with("--") => (n.to_owned(), Some(v.to_owned())),
            _ => (arg.clone(), None),
        };
        let name = if name == "--out" {
            "-o".to_owned()
        } else {
            name
        };
        if !names.contains(&name.as_str()) {
            return Err(format!("unrecognized arguments: {arg}"));
        }
        let value = match inline {
            Some(v) => v,
            None => it
                .next()
                .cloned()
                .ok_or_else(|| format!("argument {name}: expected one argument"))?,
        };
        opts.insert(name, value);
    }
    Ok((opts, pos))
}

fn read_json(path: &str) -> Result<Value, String> {
    let t = std::fs::read_to_string(path).map_err(|e| format!("{path}: {e}"))?;
    loads(&t).map_err(|e| format!("{path}: {e}"))
}

fn write(path: &str, body: &str) -> std::io::Result<()> {
    let mut f = og_io::AtomicFile::create(path)?;
    f.write_all(body.as_bytes())?;
    f.commit()
}

/// `og report index build|query ARGS`: 0 when the index is written or a
/// query names an rsID or a gene, 2 for a position, a query that is none of
/// them or a command line it cannot read, 1 when a file cannot be read or
/// written.
pub fn cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    if args.is_empty() {
        let _ = writeln!(err, "{USAGE}\nog report index: error: give build or query");
        return 2;
    }
    if args[0] == "-h" || args[0] == "--help" {
        let _ = writeln!(out, "{USAGE}");
        return 0;
    }
    let usage = |err: &mut dyn Write, why: &str| {
        let _ = writeln!(err, "{USAGE}\nog report index: error: {why}");
        2
    };
    let fail = |err: &mut dyn Write, why: String| {
        let _ = writeln!(err, "og report index: {why}");
        1
    };
    let rest = &args[1..];
    if rest.iter().any(|a| a == "-h" || a == "--help") {
        let _ = writeln!(out, "{USAGE}");
        return 0;
    }
    match args[0].as_str() {
        "build" => {
            let (o, pos) = match options(rest, &["--summary", "--sample-dir", "-o"]) {
                Ok(x) => x,
                Err(why) => return usage(err, &why),
            };
            if let Some(p) = pos.first() {
                return usage(err, &format!("unrecognized arguments: {p}"));
            }
            let (Some(summary), Some(dest)) = (o.get("--summary"), o.get("-o")) else {
                return usage(err, "give --summary and -o");
            };
            let s = match read_json(summary) {
                Ok(s) => s,
                Err(why) => return fail(err, why),
            };
            let idx = match build(&s, o.get("--sample-dir").map(String::as_str), &[], err) {
                Ok(i) => i,
                Err(why) => return fail(err, format!("{summary}: {why}")),
            };
            let opts = Dump {
                indent: Some(1),
                ..Dump::default()
            };
            let body = match dumps(&idx, &opts) {
                Ok(b) => b + "\n",
                Err(e) => return fail(err, format!("{dest}: {e}")),
            };
            if let Err(e) = write(dest, &body) {
                return fail(err, format!("{dest}: {e}"));
            }
            let n = match &idx {
                Value::Object(o) => match o.get("terms") {
                    Some(Value::Object(t)) => t.len(),
                    _ => 0,
                },
                _ => 0,
            };
            let _ = writeln!(out, "{n} terms");
            0
        }
        "query" => {
            let (o, pos) = match options(rest, &["--index"]) {
                Ok(x) => x,
                Err(why) => return usage(err, &why),
            };
            let (Some(path), [q]) = (o.get("--index"), pos.as_slice()) else {
                return usage(err, "give --index and one QUERY");
            };
            let idx = match read_json(path) {
                Ok(i) => i,
                Err(why) => return fail(err, why),
            };
            let a = query(&idx, q);
            for line in query_lines(&a) {
                let _ = writeln!(out, "{line}");
            }
            if matches!(a.kind, Kind::Rsid | Kind::Gene) {
                0
            } else {
                2
            }
        }
        other => usage(
            err,
            &format!("argument cmd: invalid choice: '{other}' (choose from 'build', 'query')"),
        ),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn v(t: &str) -> Value {
        loads(t).unwrap()
    }

    #[test]
    fn a_locus_comes_from_chrom_and_pos_or_the_variant() {
        assert_eq!(
            locus_of(&v(r#"{"chrom": "chr7", "pos": 117559590}"#)).as_deref(),
            Some("chr7:117559590")
        );
        assert_eq!(
            locus_of(&v(r#"{"chrom": "chr7", "pos": "12"}"#)).as_deref(),
            Some("chr7:12")
        );
        assert_eq!(
            locus_of(&v(r#"{"variant": "chrM:11778 G>A"}"#)).as_deref(),
            Some("chrM:11778")
        );
        assert_eq!(
            locus_of(&v(r#"{"variant": "13:2000_A"}"#)),
            None,
            "a word right after the position"
        );
        assert_eq!(locus_of(&v(r#"{"variant": "chr123:5"}"#)), None);
        assert_eq!(locus_of(&v(r#"{"variant": "NM_1:c.5"}"#)), None);
        assert_eq!(
            locus_of(&v(r#"{"chrom": "", "pos": 5, "variant": "X:5-9"}"#)).as_deref(),
            Some("X:5")
        );
    }

    #[test]
    fn existing_variation_keeps_the_rsids() {
        assert_eq!(
            rs_ids("rs123&COSV1,RS9").into_iter().collect::<Vec<_>>(),
            ["rs123", "rs9"]
        );
        assert!(rs_ids("").is_empty());
    }

    #[test]
    fn loci_and_words_read_a_row_like_the_live_veil() {
        assert_eq!(
            loci("CFTR chr7:1,000-900 x13:5 (MT:3) chrX:7"),
            [
                ("chr7".to_owned(), 900, 1000),
                ("chrMT".to_owned(), 3, 3),
                ("chrX".to_owned(), 7, 7)
            ]
        );
        let w = words("in kcnh2. HLA-A");
        assert!(
            w.contains("KCNH2") && w.contains("KCNH2.") && w.contains("HLA-A"),
            "{w:?}"
        );
    }

    #[test]
    fn the_veil_takes_the_exact_locus_first_then_gene_or_overlap() {
        let spans = [
            Closable {
                id: "a".into(),
                gene: Some("BRCA2".into()),
                locus: None,
            },
            Closable {
                id: "b".into(),
                gene: None,
                locus: Some(("chr13".into(), 2000, 2000)),
            },
            Closable {
                id: "c".into(),
                gene: None,
                locus: Some(("chr7".into(), 100, 500)),
            },
        ];
        assert_eq!(
            veil_of(&spans, "BRCA2 chr13:2000 G>A").as_deref(),
            Some("b")
        );
        assert_eq!(veil_of(&spans, "BRCA2 gene=BRCA2").as_deref(), Some("a"));
        assert_eq!(veil_of(&spans, "TP53 chr7:300").as_deref(), Some("c"));
        assert_eq!(veil_of(&spans, "TP53 chr7:501"), None);
        assert_eq!(veil_of(&[], "BRCA2"), None);
    }
}
