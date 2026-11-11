//! `og str catalog`: one repeat catalog for ExpansionHunter and Stranger
//! from STRchive (Hiatt et al. 2025, PMID 40140942).
//!
//! ```text
//! og str catalog --eh STRchive-disease-loci-<v>.hg38.expansionhunter.json \
//!                --loci STRchive-loci-<v>.json --release <v> --out catalog.json \
//!                [--ranges str_ranges.tsv] [--extra assets/str_extra_loci.json]
//! ```
//!
//! STRchive publishes the disease repeats as an ExpansionHunter catalog and,
//! apart, every locus with its benign, intermediate and pathogenic sizes and
//! its motifs. Stranger needs both in one file, so this writes one entry per
//! ExpansionHunter locus: STRchive's fields, the locus named as
//! ExpansionHunter's bundled catalog names it (the gene in capitals,
//! STRchive's id where one gene has two loci), the pathogenic repeat of a
//! locus of several as `VariantId` `LocusId`, Stranger's `NormalMax` and
//! `PathologicMin` (or `NoThresholds` with the reason: a reference repeat
//! longer than [`SHORT_READ_LIMIT`] bases, no pathogenic size, or a
//! contraction), the disease fields Stranger copies, and for a locus where
//! the motif decides, the motifs `og str motif` counts. At RFC1, AAAGG is
//! kept benign at the sizes short reads see (Cortese et al. 2019, PMID
//! 30926972), though STRchive lists it pathogenic. `--extra` appends
//! repeats STRchive leaves out, each checked first. `--ranges` writes the
//! size ranges as a table. `docs/09-str-expansions.md` has the rest.
//!
//! The output is what `bin/strchive_convert.py` wrote, byte for byte
//! (`tests/parity`, tool `str-catalog`). Deliberate differences, all on
//! input the Python could not take: every input error exits 2 with
//! `ERROR: ...` and writes nothing, where some (a `null` coordinate, a
//! locus structure without a repeat, a value of the wrong type) were a
//! Python traceback and exit 1; a JSON input that is not a list, or a
//! field that is not text where text belongs, is named in the message.

use indexmap::IndexMap;
use og_filters::{os_error, py_repr};
use pycompat::json::{self, Dump, JsonError, Value};
use std::cmp::Ordering;
use std::io::Write;

/// Bases: twice a 150-base read.
pub const SHORT_READ_LIMIT: i128 = 300;
/// A motif read in the reads is a short tandem repeat unit, not a VNTR's.
const MIN_PERIOD: usize = 2;
const MAX_PERIOD: usize = 6;
/// Listed pathogenic, benign at the sizes short reads see.
const SIZE_DEPENDENT: &[(&str, &[&str])] = &[]; // RED CONTROL: the RFC1 correction dropped

const EXTRA_FIELDS: [&str; 9] = [
    "LocusId",
    "LocusStructure",
    "ReferenceRegion",
    "VariantType",
    "NoThresholds",
    "Disease",
    "Source",
    "SourceId",
    "SourceDisplay",
];

const RANGE_COLUMNS: &str = "locus\tsource_id\tnormal_max\tintermediate_min\tintermediate_max\tpathogenic_min\tno_thresholds";

/// An input the catalog cannot be built from; printed as `ERROR: <text>`.
#[derive(Debug, Clone, PartialEq)]
pub struct Bad(pub String);

type R<T> = Result<T, Bad>;

/// Catalog entries, each a JSON object in its key order.
type Catalog = Vec<IndexMap<String, Value>>;

fn bad<T>(msg: impl Into<String>) -> R<T> {
    Err(Bad(msg.into()))
}

/// `d[key]` on a JSON object: a missing key reads as Python's `KeyError`.
fn get<'a>(d: &'a IndexMap<String, Value>, key: &str) -> R<&'a Value> {
    d.get(key).ok_or_else(|| Bad(py_repr(key)))
}

/// `d.get(key)`, with JSON null and a missing key alike as `None`.
fn opt<'a>(d: &'a IndexMap<String, Value>, key: &str) -> Option<&'a Value> {
    d.get(key).filter(|v| !matches!(v, Value::Null))
}

fn obj<'a>(v: &'a Value, what: &str) -> R<&'a IndexMap<String, Value>> {
    match v {
        Value::Object(o) => Ok(o),
        other => bad(format!("{what} is not an object: {}", repr(other))),
    }
}

fn text<'a>(v: &'a Value, what: &str) -> R<&'a str> {
    match v {
        Value::Str(s) => Ok(s),
        other => bad(format!("{what} is not text: {}", repr(other))),
    }
}

/// Python truthiness of a JSON value.
fn truthy(v: Option<&Value>) -> bool {
    match v {
        None | Some(Value::Null) | Some(Value::Bool(false)) => false,
        Some(Value::Int(0)) => false,
        Some(Value::Float(f)) => *f != 0.0,
        Some(Value::Str(s)) => !s.is_empty(),
        Some(Value::Array(a)) => !a.is_empty(),
        Some(Value::Object(o)) => !o.is_empty(),
        Some(_) => true,
    }
}

/// Python's `repr()` of a JSON value as `json.load` returns it.
pub fn repr(v: &Value) -> String {
    match v {
        Value::Null => "None".into(),
        Value::Bool(true) => "True".into(),
        Value::Bool(false) => "False".into(),
        Value::Int(i) => i.to_string(),
        Value::Float(f) => pycompat::fmt::repr(*f),
        Value::Str(s) => py_repr(s),
        Value::Array(a) => format!("[{}]", a.iter().map(repr).collect::<Vec<_>>().join(", ")),
        Value::Object(o) => format!(
            "{{{}}}",
            o.iter()
                .map(|(k, v)| format!("{}: {}", py_repr(k), repr(v)))
                .collect::<Vec<_>>()
                .join(", ")
        ),
    }
}

/// Python's `str()` of a JSON value: text as it is, the rest as `repr`.
fn py_str(v: &Value) -> String {
    match v {
        Value::Str(s) => s.clone(),
        other => repr(other),
    }
}

/// A Python number: an int stays an int, a float makes the result a float.
#[derive(Debug, Clone, Copy)]
enum Num {
    I(i128),
    F(f64),
}

impl Num {
    fn of(v: &Value, what: &str) -> R<Num> {
        match v {
            Value::Int(i) => Ok(Num::I(*i)),
            Value::Float(f) => Ok(Num::F(*f)),
            Value::Bool(b) => Ok(Num::I(i128::from(*b))),
            other => bad(format!("{what} is not a number: {}", repr(other))),
        }
    }

    fn f(self) -> f64 {
        match self {
            // Exact for every coordinate and size below 2^53.
            Num::I(i) => i as f64,
            Num::F(f) => f,
        }
    }

    fn sub(self, o: Num) -> Num {
        match (self, o) {
            (Num::I(a), Num::I(b)) => Num::I(a - b),
            _ => Num::F(self.f() - o.f()),
        }
    }

    fn cmp(self, o: Num) -> Ordering {
        match (self, o) {
            (Num::I(a), Num::I(b)) => a.cmp(&b),
            _ => self.f().partial_cmp(&o.f()).unwrap_or(Ordering::Equal),
        }
    }

    fn min(self, o: Num) -> Num {
        if o.cmp(self) == Ordering::Less {
            o
        } else {
            self
        }
    }

    fn max(self, o: Num) -> Num {
        if o.cmp(self) == Ordering::Greater {
            o
        } else {
            self
        }
    }
}

impl std::fmt::Display for Num {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Num::I(i) => write!(f, "{i}"),
            Num::F(x) => f.write_str(&pycompat::fmt::repr(*x)),
        }
    }
}

/// `as_list(x)`: a list as it is, anything else as a list of one.
fn as_list(v: &Value) -> Vec<&Value> {
    match v {
        Value::Array(a) => a.iter().collect(),
        other => vec![other],
    }
}

/// `chrom, start, end` of `chr..:start-end` (the pattern
/// `^(chr[^:]+):(\d+)-(\d+)$`, whose `$` also takes one final newline).
fn region(v: &Value) -> R<(String, i128, i128)> {
    let t = text(v, "a reference region")?;
    let fail = || Bad(format!("not a region: {}", py_repr(t)));
    let body = t.strip_suffix('\n').unwrap_or(t);
    let (chrom, rest) = body.split_once(':').ok_or_else(fail)?;
    if !chrom.starts_with("chr") || chrom.len() <= 3 {
        return Err(fail());
    }
    let (a, b) = rest.split_once('-').ok_or_else(fail)?;
    let num = |s: &str| -> R<i128> {
        if s.is_empty() || !s.chars().all(pycompat::text::is_decimal) {
            return Err(fail());
        }
        pycompat::num::parse_int(s).map_err(|_| fail())
    };
    Ok((chrom.to_owned(), num(a)?, num(b)?))
}

/// The motifs of a locus structure: each `(...)` followed by `*`.
pub fn motifs_of(structure: &str) -> Vec<&str> {
    let mut out = Vec::new();
    let mut i = 0;
    while let Some(open) = structure[i..].find('(').map(|k| k + i) {
        match structure[open + 1..].find(')').map(|k| k + open + 1) {
            Some(close) if structure[close + 1..].starts_with('*') => {
                out.push(&structure[open + 1..close]);
                i = close + 2;
            }
            _ => i = open + 1,
        }
    }
    out
}

fn structure_of(e: &IndexMap<String, Value>) -> R<&str> {
    text(get(e, "LocusStructure")?, "LocusStructure")
}

/// `a.upper() in forms(b.upper())`: the same motif, rotated or on the other
/// strand, by the motif counter's own rule (`og str motif`).
fn same_motif(a: &str, b: &str) -> bool {
    let (a, b) = (pycompat::text::upper(a), pycompat::text::upper(b));
    !b.is_empty()
        && og_bam::motif::forms(b.as_bytes())
            .iter()
            .any(|f| f == a.as_bytes())
}

/// The strings of a JSON value iterated as Python would: a list's items, a
/// text's characters, nothing for a missing or empty one.
fn strings(v: Option<&Value>, what: &str) -> R<Vec<String>> {
    match v {
        None => Ok(Vec::new()),
        Some(x) if !truthy(Some(x)) => Ok(Vec::new()),
        Some(Value::Array(a)) => a.iter().map(|x| text(x, what).map(str::to_owned)).collect(),
        Some(Value::Str(s)) => Ok(s.chars().map(String::from).collect()),
        Some(other) => bad(format!("{what} is not a list: {}", repr(other))),
    }
}

fn chars(s: &str) -> usize {
    s.chars().count()
}

/// STRchive id -> LocusId: the gene in capitals, or the STRchive id where
/// one gene has two loci.
fn names(
    eh: &[&IndexMap<String, Value>],
    by: &IndexMap<String, &IndexMap<String, Value>>,
) -> R<IndexMap<String, String>> {
    let mut genes: IndexMap<String, Vec<String>> = IndexMap::new();
    for e in eh {
        let sid = text(get(e, "LocusId")?, "LocusId")?;
        let gene = text(get(by[sid], "gene")?, "gene")?;
        genes
            .entry(pycompat::text::upper(gene))
            .or_default()
            .push(sid.to_owned());
    }
    let mut out = IndexMap::new();
    for (gene, ids) in &genes {
        for i in ids {
            out.insert(
                i.clone(),
                if ids.len() == 1 {
                    gene.clone()
                } else {
                    i.clone()
                },
            );
        }
    }
    Ok(out)
}

/// Which of a locus's repeats is the pathogenic one.
fn pathogenic_index(e: &IndexMap<String, Value>, rec: &IndexMap<String, Value>) -> R<usize> {
    let regions = as_list(get(e, "ReferenceRegion")?);
    if regions.len() == 1 {
        return Ok(0);
    }
    let sid = py_str(get(e, "LocusId")?);
    let ls = structure_of(e)?;
    let mot = motifs_of(ls);
    if mot.len() != regions.len() {
        return bad(format!(
            "{sid}: {} motifs in {ls} for {} regions",
            mot.len(),
            regions.len()
        ));
    }
    let mut want = Vec::new();
    if truthy(rec.get("locus_structure")) {
        for s in as_list(get(rec, "locus_structure")?) {
            let s = obj(s, "a locus_structure item")?;
            if s.get("type") == Some(&Value::Str("pathogenic_repeat".into())) {
                want.push(text(get(s, "motif")?, "motif")?.to_owned());
            }
        }
    }
    let hits: Vec<usize> = (0..mot.len())
        .filter(|&i| want.iter().any(|w| same_motif(mot[i], w)))
        .collect();
    let overlap = |i: usize| -> R<Num> {
        let (_, a, b) = region(regions[i])?;
        let stop = Num::of(get(rec, "stop_hg38")?, "stop_hg38")?;
        let start = Num::of(get(rec, "start_hg38")?, "start_hg38")?;
        Ok(Num::I(0).max(Num::I(b).min(stop).sub(Num::I(a).max(start))))
    };
    let pool: Vec<usize> = if hits.is_empty() {
        (0..regions.len()).collect()
    } else {
        hits
    };
    let ov = pool.iter().map(|&i| overlap(i)).collect::<R<Vec<_>>>()?;
    let Some(best) = ov.iter().copied().reduce(|x, y| x.max(y)) else {
        return bad(format!("{sid}: no reference region"));
    };
    let top: Vec<usize> = pool
        .iter()
        .zip(&ov)
        .filter(|(_, o)| o.cmp(best) == Ordering::Equal)
        .map(|(i, _)| *i)
        .collect();
    if top.len() != 1 {
        let list = top
            .iter()
            .map(usize::to_string)
            .collect::<Vec<_>>()
            .join(", ");
        return bad(format!(
            "{sid}: cannot tell its pathogenic repeat (candidates [{list}])"
        ));
    }
    Ok(top[0])
}

/// (benign, pathogenic) motif lists for a locus where the motif decides.
fn motif_classes(
    sid: &str,
    rec: &IndexMap<String, Value>,
) -> R<Option<(Vec<String>, Vec<String>)>> {
    let up = |key: &str| -> R<Vec<String>> {
        Ok(strings(rec.get(key), key)?
            .iter()
            .map(|m| pycompat::text::upper(m))
            .collect())
    };
    let path = up("pathogenic_motif_reference_orientation")?;
    let refm = up("reference_motif_reference_orientation")?;
    let ben = up("benign_motif_reference_orientation")?;
    if path.is_empty() || path.iter().chain(&refm).any(|m| m.contains('N')) {
        return Ok(None);
    }
    let size = chars(&path[0]);
    if !(MIN_PERIOD..=MAX_PERIOD).contains(&size) || path.iter().any(|m| chars(m) != size) {
        return Ok(None);
    }
    let mut moved: Vec<&str> = SIZE_DEPENDENT
        .iter()
        .find(|(k, _)| *k == sid)
        .map(|(_, v)| v.to_vec())
        .unwrap_or_default();
    moved.sort_unstable();
    let pathogenic: Vec<String> = path
        .iter()
        .filter(|m| !moved.contains(&m.as_str()))
        .cloned()
        .collect();
    let mut benign: Vec<String> = Vec::new();
    for m in refm
        .iter()
        .chain(&ben)
        .cloned()
        .chain(moved.iter().map(|s| (*s).to_owned()))
    {
        if chars(&m) == size
            && !m.contains('N')
            && !pathogenic.iter().any(|p| same_motif(&m, p))
            && !benign.iter().any(|b| same_motif(&m, b))
        {
            benign.push(m);
        }
    }
    if pathogenic
        .iter()
        .all(|p| refm.iter().any(|r| same_motif(p, r)))
    {
        // The pathogenic motif is the reference one: the size decides.
        return Ok(None);
    }
    Ok(Some((benign, pathogenic)))
}

fn strs(v: &[String]) -> Value {
    Value::Array(v.iter().map(|s| Value::Str(s.clone())).collect())
}

/// The converted catalog: one entry per ExpansionHunter locus.
pub fn convert(eh: &Value, loci: &Value, release: &str) -> R<Vec<IndexMap<String, Value>>> {
    let eh: Vec<&IndexMap<String, Value>> = match eh {
        Value::Array(a) => a
            .iter()
            .map(|e| obj(e, "an ExpansionHunter locus"))
            .collect::<R<_>>()?,
        other => {
            return bad(format!(
                "the ExpansionHunter catalog is not a list: {}",
                repr(other)
            ));
        }
    };
    let by = by_id(loci)?;
    let mut missing = Vec::new();
    for e in &eh {
        let sid = text(get(e, "LocusId")?, "LocusId")?;
        if !by.contains_key(sid) {
            missing.push(sid);
        }
    }
    if !missing.is_empty() {
        return bad(format!(
            "ExpansionHunter loci without a STRchive record: {}",
            missing.join(", ")
        ));
    }
    let ids = names(&eh, &by)?;
    let mut out = Vec::new();
    let mut seen: Vec<String> = Vec::new();
    for e in &eh {
        let sid = text(get(e, "LocusId")?, "LocusId")?;
        let rec = by[sid];
        let lid = &ids[sid];
        let rr = get(e, "ReferenceRegion")?;
        let regions = as_list(rr);
        let mut entry = IndexMap::new();
        entry.insert("LocusId".to_owned(), Value::Str(lid.clone()));
        entry.insert(
            "LocusStructure".to_owned(),
            get(e, "LocusStructure")?.clone(),
        );
        entry.insert("ReferenceRegion".to_owned(), rr.clone());
        entry.insert("VariantType".to_owned(), get(e, "VariantType")?.clone());
        let pi = pathogenic_index(e, rec)?;
        let ls = structure_of(e)?;
        let mot = motifs_of(ls);
        if regions.len() > 1 {
            let mut vids: Vec<String> = Vec::new();
            for (i, m) in mot.iter().enumerate() {
                let v = if i == pi {
                    lid.clone()
                } else {
                    format!("{lid}_{m}")
                };
                let v = if vids.contains(&v) {
                    format!("{v}_{}", i + 1)
                } else {
                    v
                };
                vids.push(v);
            }
            entry.insert("VariantId".to_owned(), strs(&vids));
        }
        let names_here: Vec<String> = match entry.get("VariantId") {
            Some(Value::Array(a)) => a.iter().map(py_str).collect(),
            _ => vec![lid.clone()],
        };
        for v in names_here {
            if seen.contains(&v) {
                return bad(format!("two variants are named {v}"));
            }
            seen.push(v);
        }
        let Some(&pregion) = regions.get(pi) else {
            return bad(format!("{sid}: no reference region"));
        };
        let (_, start, end) = region(pregion)?;
        let Some(pm) = mot.get(pi) else {
            return bad(format!("{sid}: no repeat in {ls}"));
        };
        let size = i128::try_from(chars(pm)).unwrap_or(1).max(1);
        let span = Num::of(get(rec, "stop_hg38")?, "stop_hg38")?
            .sub(Num::of(get(rec, "start_hg38")?, "start_hg38")?);
        let bmax = opt(rec, "benign_max");
        let pmin = opt(rec, "pathogenic_min");
        if span.cmp(Num::I(SHORT_READ_LIMIT)) == Ordering::Greater {
            entry.insert(
                "NoThresholds".to_owned(),
                Value::Str(format!(
                    "the reference repeat is {span} bases, longer than {SHORT_READ_LIMIT}: short reads cannot size it"
                )),
            );
        } else if let Some(pmin) = pmin {
            let pn = Num::of(pmin, "pathogenic_min")?;
            let contraction = match bmax {
                Some(b) => pn.cmp(Num::of(b, "benign_max")?) != Ordering::Greater,
                None => false,
            };
            if let (true, Some(b)) = (contraction, bmax) {
                entry.insert(
                    "NoThresholds".to_owned(),
                    Value::Str(format!(
                        "pathogenic from {} units, below the benign maximum {}: a contraction, which Stranger does not grade",
                        py_str(pmin),
                        py_str(b)
                    )),
                );
            } else {
                let normal = match bmax {
                    Some(b) => b.clone(),
                    None => Value::Int((end - start).div_euclid(size)),
                };
                if Num::of(&normal, "benign_max")?.cmp(pn) != Ordering::Less {
                    return bad(format!(
                        "{sid}: normal maximum {} is not below pathogenic {}",
                        py_str(&normal),
                        py_str(pmin)
                    ));
                }
                entry.insert("NormalMax".to_owned(), normal);
                entry.insert("PathologicMin".to_owned(), pmin.clone());
            }
        } else {
            entry.insert(
                "NoThresholds".to_owned(),
                Value::Str("STRchive gives no pathogenic size".to_owned()),
            );
        }
        if regions.len() > 1 {
            entry.insert("PathologicRegion".to_owned(), pregion.clone());
        }
        let disease = match rec.get("disease_id") {
            d @ Some(v) if truthy(d) => v.clone(),
            _ => Value::Str(sid.split('_').next().unwrap_or(sid).to_owned()),
        };
        entry.insert("Disease".to_owned(), disease);
        if truthy(rec.get("inheritance")) {
            let inh = strings(rec.get("inheritance"), "inheritance")?;
            entry.insert("InheritanceMode".to_owned(), Value::Str(inh.join("/")));
        }
        entry.insert("Source".to_owned(), Value::Str("STRchive".to_owned()));
        entry.insert("SourceId".to_owned(), Value::Str(sid.to_owned()));
        entry.insert(
            "SourceDisplay".to_owned(),
            Value::Str(format!("STRchive-{release}")),
        );
        if let Some((b, p)) = motif_classes(sid, rec)? {
            entry.insert("BenignMotifs".to_owned(), strs(&b));
            entry.insert("PathogenicMotifs".to_owned(), strs(&p));
        }
        out.push(entry);
    }
    Ok(out)
}

/// `{r["id"]: r for r in loci}`: the last record of a repeated id wins.
fn by_id(loci: &Value) -> R<IndexMap<String, &IndexMap<String, Value>>> {
    let Value::Array(recs) = loci else {
        return bad(format!(
            "the STRchive loci file is not a list: {}",
            repr(loci)
        ));
    };
    let mut by = IndexMap::new();
    for r in recs {
        let r = obj(r, "a STRchive record")?;
        by.insert(text(get(r, "id")?, "id")?.to_owned(), r);
    }
    Ok(by)
}

/// The `--extra` entries, checked, to append to the converted catalog.
pub fn extra_entries(
    extra: &Value,
    cat: &[IndexMap<String, Value>],
) -> R<Vec<IndexMap<String, Value>>> {
    let list = match extra {
        Value::Array(a) if !a.is_empty() => a,
        _ => return bad("the extra catalog is not a non-empty list of loci"),
    };
    let mut taken: Vec<String> = Vec::new();
    for e in cat {
        taken.push(py_str(&e["LocusId"]));
        if let Some(v) = e.get("VariantId") {
            taken.extend(as_list(v).into_iter().map(py_str));
        }
    }
    let mut out = Vec::new();
    for item in list {
        let d = match item {
            Value::Object(o) => Some(o),
            _ => None,
        };
        let lid = d.and_then(|o| o.get("LocusId"));
        let where_ = format!(
            "extra locus {}",
            if truthy(lid) {
                py_str(lid.unwrap_or(&Value::Null))
            } else {
                repr(item)
            }
        );
        let missing: Vec<&str> = EXTRA_FIELDS
            .iter()
            .copied()
            .filter(|k| !truthy(d.and_then(|o| o.get(*k))))
            .collect();
        let Some(e) = d.filter(|_| missing.is_empty()) else {
            return bad(format!("{where_}: no {}", missing.join(", ")));
        };
        let lid = get(e, "LocusId")?;
        let ls = structure_of(e)?;
        let mot = motifs_of(ls);
        let regions = as_list(get(e, "ReferenceRegion")?);
        let types = as_list(get(e, "VariantType")?);
        if mot.is_empty() || mot.len() != regions.len() || types.len() != regions.len() {
            return bad(format!(
                "{where_}: {} repeats in {ls}, {} regions, {} variant types",
                mot.len(),
                regions.len(),
                types.len()
            ));
        }
        let one = Value::Array(vec![lid.clone()]);
        let vids = as_list(e.get("VariantId").unwrap_or(&one));
        if regions.len() > 1 && vids.len() != regions.len() {
            return bad(format!(
                "{where_}: {n} repeats need {n} VariantId names",
                n = regions.len()
            ));
        }
        for (m, r) in mot.iter().zip(&regions) {
            let (_, a, b) = region(r)?;
            let whole = !m.is_empty()
                && m.bytes().all(|c| b"ACGT".contains(&c))
                && b > a
                && (b - a) % i128::try_from(m.len()).unwrap_or(1) == 0;
            if !whole {
                return bad(format!(
                    "{where_}: region {} is not a whole number of {m} motifs",
                    py_str(r)
                ));
            }
        }
        let mut names = vec![py_str(lid)];
        for v in &vids {
            let v = py_str(v);
            if !names.contains(&v) {
                names.push(v);
            }
        }
        for v in names {
            if taken.contains(&v) {
                return bad(format!("{where_}: the name {v} is already in the catalog"));
            }
            taken.push(v);
        }
        if e.contains_key("NormalMax") || e.contains_key("PathologicMin") {
            return bad(format!(
                "{where_}: an extra locus is read by its counts and takes no Stranger thresholds"
            ));
        }
        out.push(e.clone());
    }
    Ok(out)
}

/// The `--ranges` table: one row per catalog entry.
pub fn ranges(cat: &[IndexMap<String, Value>], loci: &Value, release: &str) -> R<String> {
    let by = by_id(loci)?;
    let empty = IndexMap::new();
    let cell = |v: Option<&Value>| {
        v.filter(|x| !matches!(x, Value::Null))
            .map(py_str)
            .unwrap_or_default()
    };
    let mut s = format!(
        "# STRchive {release} repeat ranges, one row per ExpansionHunter locus (bin/strchive_convert.py)\n{RANGE_COLUMNS}\n"
    );
    for e in cat {
        let source_id = get(e, "SourceId")?;
        let r = if e.get("Source") == Some(&Value::Str("STRchive".into())) {
            let k = text(source_id, "SourceId")?;
            *by.get(k).ok_or_else(|| Bad(py_repr(k)))?
        } else {
            &empty
        };
        let row = [
            text(get(e, "LocusId")?, "LocusId")?.to_owned(),
            text(source_id, "SourceId")?.to_owned(),
            cell(e.get("NormalMax")),
            cell(r.get("intermediate_min")),
            cell(r.get("intermediate_max")),
            cell(e.get("PathologicMin")),
            match e.get("NoThresholds") {
                Some(v) => text(v, "NoThresholds")?.to_owned(),
                None => String::new(),
            },
        ];
        s.push_str(&row.join("\t"));
        s.push('\n');
    }
    Ok(s)
}

/// Python's `JSONDecodeError` text: what, line, column and character.
fn decode_error(e: &JsonError, doc: &str) -> String {
    match e {
        JsonError::Decode(what, at) => {
            let at = (*at).min(doc.len());
            let head = doc.get(..at).unwrap_or(doc);
            let pos = head.chars().count();
            let line = head.matches('\n').count() + 1;
            let col = match head.rfind('\n') {
                Some(nl) => head[nl..].chars().count(),
                None => pos + 1,
            };
            format!("{what}: line {line} column {col} (char {pos})")
        }
        other => other.to_string(),
    }
}

fn load(path: &str) -> R<Value> {
    let bytes = std::fs::read(path).map_err(|e| Bad(os_error(&e, path)))?;
    let doc = String::from_utf8(bytes).map_err(|e| {
        let at = e.utf8_error().valid_up_to();
        let b = e.as_bytes().get(at).copied().unwrap_or(0);
        Bad(format!(
            "'utf-8' codec can't decode byte 0x{b:02x} in position {at}: invalid start byte"
        ))
    })?;
    json::loads(&doc).map_err(|e| Bad(decode_error(&e, &doc)))
}

fn dump(cat: &[IndexMap<String, Value>]) -> String {
    let v = Value::Array(cat.iter().cloned().map(Value::Object).collect());
    let mut s = json::dumps(
        &v,
        &Dump {
            indent: Some(1),
            ..Dump::default()
        },
    )
    .unwrap_or_default();
    s.push('\n');
    s
}

#[derive(Default)]
struct Args {
    eh: Option<String>,
    loci: Option<String>,
    release: Option<String>,
    out: Option<String>,
    ranges: Option<String>,
    extra: Option<String>,
}

const USAGE: &str = "usage: og str catalog --eh EH --loci LOCI --release RELEASE --out OUT [--ranges RANGES] [--extra EXTRA]";

fn parse(args: &[String]) -> Result<Args, String> {
    let mut a = Args::default();
    let mut it = args.iter();
    while let Some(arg) = it.next() {
        let (key, inline) = match arg.split_once('=') {
            Some((k, v)) if k.starts_with("--") => (k, Some(v.to_owned())),
            _ => (arg.as_str(), None),
        };
        let slot = match key {
            "--eh" => &mut a.eh,
            "--loci" => &mut a.loci,
            "--release" => &mut a.release,
            "--out" => &mut a.out,
            "--ranges" => &mut a.ranges,
            "--extra" => &mut a.extra,
            other => return Err(format!("unrecognized arguments: {other}")),
        };
        let v = match inline {
            Some(v) => v,
            None => it
                .next()
                .cloned()
                .ok_or_else(|| format!("argument {key}: expected one argument"))?,
        };
        *slot = Some(v);
    }
    let need: Vec<&str> = [
        ("--eh", &a.eh),
        ("--loci", &a.loci),
        ("--release", &a.release),
        ("--out", &a.out),
    ]
    .iter()
    .filter(|(_, v)| v.is_none())
    .map(|(k, _)| *k)
    .collect();
    if !need.is_empty() {
        return Err(format!(
            "the following arguments are required: {}",
            need.join(", ")
        ));
    }
    Ok(a)
}

/// The catalog and, with `--ranges`, the ranges table, both checked before
/// anything is written.
fn build(a: &Args, release: &str) -> R<(Catalog, Option<String>)> {
    let eh = load(a.eh.as_deref().unwrap_or_default())?;
    let loci = load(a.loci.as_deref().unwrap_or_default())?;
    let mut cat = convert(&eh, &loci, release)?;
    if let Some(x) = &a.extra {
        let more = extra_entries(&load(x)?, &cat)?;
        cat.extend(more);
    }
    let table = match &a.ranges {
        Some(_) => Some(ranges(&cat, &loci, release)?),
        None => None,
    };
    Ok((cat, table))
}

/// Writes both files under temporary names, then renames the catalog and
/// then the ranges into place.
fn write_all(out: &str, ranges: Option<(&str, &str)>, catalog: &str) -> std::io::Result<()> {
    let mut c = og_io::AtomicFile::create(out)?;
    c.write_all(catalog.as_bytes())?;
    let r = match ranges {
        Some((path, text)) => {
            let mut f = og_io::AtomicFile::create(path)?;
            f.write_all(text.as_bytes())?;
            Some(f)
        }
        None => None,
    };
    c.commit()?;
    if let Some(f) = r {
        f.commit()?;
    }
    Ok(())
}

/// The command line; returns the exit code: 0, 1 when a file cannot be
/// written, 2 on any input error (writing nothing) or a wrong command line.
pub fn cli(args: &[String]) -> u8 {
    let a = match parse(args) {
        Ok(a) => a,
        Err(m) => {
            eprintln!("{USAGE}\nog str catalog: error: {m}");
            return 2;
        }
    };
    let release = a.release.clone().unwrap_or_default();
    let out = a.out.clone().unwrap_or_default();
    let (cat, table) = match build(&a, &release) {
        Ok(x) => x,
        Err(Bad(m)) => {
            eprintln!("ERROR: {m}");
            return 2;
        }
    };
    let ranges = a.ranges.as_deref().zip(table.as_deref());
    if let Err(e) = write_all(&out, ranges, &dump(&cat)) {
        eprintln!("ERROR: {e}");
        return 1;
    }
    let graded = cat.iter().filter(|e| e.contains_key("NormalMax")).count();
    let motif = cat
        .iter()
        .filter(|e| e.contains_key("PathogenicMotifs"))
        .count();
    let extra = cat
        .iter()
        .filter(|e| e.get("Source") != Some(&Value::Str("STRchive".into())))
        .count();
    let from = match &a.extra {
        Some(x) => format!(", {extra} from {x}"),
        None => String::new(),
    };
    println!(
        "STRchive {release}: {} loci, {graded} with Stranger thresholds, {motif} read by motif{from}: {out}",
        cat.len()
    );
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    fn v(s: &str) -> Value {
        json::loads(s).unwrap()
    }

    const RFC1_EH: &str = r#"[{"LocusId": "CANVAS_RFC1", "LocusStructure": "(AAAAG)*",
        "ReferenceRegion": "chr4:39348424-39348483", "VariantType": "Repeat"}]"#;
    const RFC1_LOCI: &str = r#"[{"id": "CANVAS_RFC1", "gene": "RFC1", "start_hg38": 39348424,
        "stop_hg38": 39348483, "benign_max": null, "pathogenic_min": 400,
        "disease_id": "CANVAS", "inheritance": ["AR"],
        "pathogenic_motif_reference_orientation": ["AAGGG", "ACAGG", "AAAGG"],
        "reference_motif_reference_orientation": ["AAAAG"],
        "benign_motif_reference_orientation": ["AAAGGG"]}]"#;

    #[test]
    fn rfc1_keeps_aaagg_benign() {
        let cat = convert(&v(RFC1_EH), &v(RFC1_LOCI), "v1").unwrap();
        let e = &cat[0];
        assert_eq!(e["LocusId"], Value::Str("RFC1".into()));
        assert_eq!(e["BenignMotifs"], v(r#"["AAAAG", "AAAGG"]"#));
        assert_eq!(e["PathogenicMotifs"], v(r#"["AAGGG", "ACAGG"]"#));
        assert_eq!(e["NormalMax"], Value::Int(11));
        assert_eq!(e["InheritanceMode"], Value::Str("AR".into()));
    }

    #[test]
    fn a_null_coordinate_is_an_input_error() {
        let loci = RFC1_LOCI.replace("\"stop_hg38\": 39348483", "\"stop_hg38\": null");
        let err = convert(&v(RFC1_EH), &v(&loci), "v1").unwrap_err();
        assert_eq!(err, Bad("stop_hg38 is not a number: None".into()));
    }

    #[test]
    fn motifs_and_regions_parse_as_the_python_patterns() {
        assert_eq!(motifs_of("(CAG)*CAACAG(CCG)*"), ["CAG", "CCG"]);
        assert_eq!(motifs_of("((A)*(B)C(D)*"), ["(A", "D"]);
        assert_eq!(
            region(&v("\"chr1:10-20\"")).unwrap(),
            ("chr1".into(), 10, 20)
        );
        assert_eq!(region(&v("\"chr1:10-20\\n\"")).unwrap().2, 20);
        for bad in [
            "\"1:10-20\"",
            "\"chr:10-20\"",
            "\"chr1:10-\"",
            "\"chr1:a-2\"",
        ] {
            assert!(region(&v(bad)).is_err(), "{bad}");
        }
    }

    #[test]
    fn repr_is_python_s() {
        assert_eq!(
            repr(&v(r#"{"a": [1, 2.0, null, true, "it's"]}"#)),
            "{'a': [1, 2.0, None, True, \"it's\"]}"
        );
    }

    #[test]
    fn json_errors_read_as_python_s() {
        let e = json::loads("[\n  1,\n  x]").unwrap_err();
        assert_eq!(
            decode_error(&e, "[\n  1,\n  x]"),
            "Expecting value: line 3 column 3 (char 9)"
        );
    }
}
