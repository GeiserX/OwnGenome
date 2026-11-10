//! PharmCAT's `report.json`, read into typed gene calls.
//!
//! PharmCAT 3.x writes `genes` as a flat `{gene: GeneReport}` map; 2.x
//! nested it by guidance source (`{source: {gene: GeneReport}}`); older
//! versions used a list, and some a `geneResults` list. All four are read.
//! A field of an unexpected JSON type reads as missing, where the Python
//! stopped with a traceback.

use crate::status::{Status, is_non_normal, is_normal};
use pycompat::json::Value;
use pycompat::text::{lower, strip};
use pycompat::value::{py_str, truthy};
use std::collections::HashMap;

/// The allele names that mean "no call".
const NO_CALL_NAMES: [&str; 5] = ["", "unknown", "none", "?", "n/a"];

/// The genes step 36 (`bin/pgx_outside_calls.py`) can pass to PharmCAT. Any
/// other outside call, such as MT-RNR1 from step 20's chrM calls, never
/// appears in its table, so the table cannot confirm or refute it.
pub const STEP36_GENES: [&str; 3] = ["CYP2D6", "HLA-A", "HLA-B"];

/// `str(v)` for a value Python formats into text; a value whose `repr` this
/// port does not carry reads as its JSON-free best effort.
pub fn text(v: &Value) -> String {
    py_str(v).unwrap_or_default()
}

/// `d.get(key)` when `d` is an object.
pub fn field<'a>(d: &'a Value, key: &str) -> Option<&'a Value> {
    match d {
        Value::Object(o) => o.get(key),
        _ => None,
    }
}

/// `d.get(key)` when it is true, as Python's `d.get(key) or ...` reads it.
pub fn truthy_field<'a>(d: &'a Value, key: &str) -> Option<&'a Value> {
    field(d, key).filter(|v| truthy(v))
}

/// The elements of a list field; anything else has none.
pub fn list<'a>(d: &'a Value, key: &str) -> &'a [Value] {
    match field(d, key) {
        Some(Value::Array(a)) => a,
        _ => &[],
    }
}

/// One gene of the report: what PharmCAT called and how to read it.
///
/// `labels` are all the diplotypes PharmCAT lists for the gene. More than
/// one means the data cannot tell them apart (positions missing from the
/// VCF, for example); when their phenotypes differ the gene is ambiguous and
/// its first diplotype must not be read as the result.
#[derive(Debug, Clone, PartialEq)]
pub struct GeneCall {
    pub gene: String,
    /// PharmCAT's `callSource` is `OUTSIDE`: another tool gave the call.
    pub outside: bool,
    /// The first diplotype, with `(1 of N possible)` when there are more.
    pub diplotype: String,
    pub phenotype: String,
    pub called: bool,
    pub labels: Vec<String>,
    /// The distinct phenotypes of the possible diplotypes, in order.
    pub phenotypes: Vec<String>,
    /// Where the phenotype came from, when not the gene level.
    pub note: Option<&'static str>,
    /// PharmCAT's `activityScore` of each possible diplotype, when it is a
    /// finite number.
    pub scores: Vec<Option<f64>>,
}

impl GeneCall {
    /// A call with one possible diplotype.
    pub fn single(gene: String, diplotype: String, phenotype: String, called: bool) -> Self {
        GeneCall {
            gene,
            outside: false,
            labels: vec![diplotype.clone()],
            phenotypes: vec![phenotype.clone()],
            diplotype,
            phenotype,
            called,
            note: None,
            scores: Vec::new(),
        }
    }

    /// The phenotype as the table and the text show it, with its note.
    pub fn shown(&self) -> String {
        match self.note {
            Some(note) => format!("{} [{note}]", self.phenotype),
            None => self.phenotype.clone(),
        }
    }

    /// Called, with possible diplotypes of different phenotypes.
    pub fn ambiguous(&self) -> bool {
        self.called && self.phenotypes.len() > 2
    }

    pub fn status(&self) -> Status {
        if !self.called {
            Status::NotCalled
        } else if self.ambiguous() {
            Status::Ambiguous
        } else if is_normal(&self.phenotype) {
            Status::Normal
        } else if is_non_normal(&self.phenotype) {
            Status::NonNormal
        } else {
            Status::Unclassified
        }
    }

    /// An outside call of a gene step 36 passes on.
    pub fn step36_outside(&self) -> bool {
        self.outside && self.called && STEP36_GENES.contains(&self.gene.as_str())
    }
}

/// `float(value)` of an activity score, when it is a finite number (a NaN
/// would fall into the ultrarapid band).
pub fn score(v: &Value) -> Option<f64> {
    let x = match v {
        Value::Bool(b) => f64::from(u8::from(*b)),
        // An activity score is a small number; Python's float() of a huge
        // int is inf or an error, neither of which converts.
        #[allow(clippy::cast_precision_loss)]
        Value::Int(i) => *i as f64,
        Value::Float(x) => *x,
        Value::Str(s) => pycompat::num::parse_float(s).ok()?,
        _ => return None,
    };
    x.is_finite().then_some(x)
}

/// (label, phenotype, called) of one diplotype object.
fn one_diplotype(dip: &Value) -> (String, String, bool) {
    // A present allele: its name, or "" (a missing or null name, or an
    // allele that is not an object).
    let names: Vec<String> = ["allele1", "allele2"]
        .iter()
        .filter_map(|a| truthy_field(dip, a))
        .map(|a| {
            field(a, "name")
                .filter(|n| truthy(n))
                .map(text)
                .unwrap_or_default()
        })
        .collect();
    let label = match truthy_field(dip, "label") {
        Some(l) => text(l),
        None => {
            let joined = names
                .iter()
                .map(|n| if n.is_empty() { "?" } else { n.as_str() })
                .collect::<Vec<_>>()
                .join("/");
            if joined.is_empty() {
                "?".to_owned()
            } else {
                joined
            }
        }
    };
    let called = if names.is_empty() {
        !lower(&label).starts_with("unknown")
    } else {
        names
            .iter()
            .any(|n| !NO_CALL_NAMES.contains(&lower(strip(n, None)).as_str()))
    };
    let phenos: Vec<String> = match field(dip, "phenotypes") {
        Some(Value::Array(a)) => a.iter().filter(|p| truthy(p)).map(text).collect(),
        // Python iterates a string by its characters.
        Some(Value::Str(s)) => s.chars().map(String::from).collect(),
        Some(Value::Object(o)) => o.keys().filter(|k| !k.is_empty()).cloned().collect(),
        _ => Vec::new(),
    };
    let mut phenotype = if !phenos.is_empty() {
        phenos.join("; ")
    } else if let Some(p) = truthy_field(dip, "phenotype") {
        text(p)
    } else if called {
        "no phenotype assigned".to_owned()
    } else {
        "No Result".to_owned()
    };
    if !called && matches!(lower(strip(&phenotype, None)).as_str(), "no result" | "n/a") {
        phenotype = "No Result".to_owned();
    }
    (label, phenotype, called)
}

/// The call of a gene object with a diplotype array, else `None`.
fn diplotype_call(name: String, g: &Value) -> Option<GeneCall> {
    if !matches!(g, Value::Object(_)) {
        return None;
    }
    let source =
        truthy_field(g, "sourceDiplotypes").or_else(|| truthy_field(g, "recommendationDiplotypes"));
    let dips: Vec<&Value> = match source {
        Some(Value::Array(a)) => a.iter().filter(|d| matches!(d, Value::Object(_))).collect(),
        _ => Vec::new(),
    };
    let parsed: Vec<(String, String, bool)> = dips.iter().map(|d| one_diplotype(d)).collect();
    let (first_label, first_phenotype, called) = parsed.first()?.clone();
    let mut phenotypes: Vec<String> = Vec::new();
    for (_, p, _) in &parsed {
        if !phenotypes.contains(p) {
            phenotypes.push(p.clone());
        }
    }
    let mut label = first_label;
    let mut phenotype = first_phenotype;
    if parsed.len() > 1 {
        label = format!("{label} (1 of {} possible)", parsed.len());
        if phenotypes.len() > 1 {
            let more = if phenotypes.len() > 5 { ", ..." } else { "" };
            phenotype = format!(
                "ambiguous: one of {}{more}",
                phenotypes[..phenotypes.len().min(5)].join(", ")
            );
        }
    }
    Some(GeneCall {
        gene: name,
        outside: matches!(field(g, "callSource"), Some(Value::Str(s)) if s == "OUTSIDE"),
        diplotype: label,
        phenotype,
        called,
        labels: parsed.into_iter().map(|p| p.0).collect(),
        phenotypes,
        note: None,
        scores: dips
            .iter()
            .map(|d| field(d, "activityScore").and_then(score))
            .collect(),
    })
}

/// What [`parse_genes`] found.
#[derive(Debug, Clone, PartialEq)]
pub enum Parsed {
    /// The gene calls, one per gene, in the report's order.
    Ok(Vec<GeneCall>),
    /// A gene section in a known layout, but no gene diplotype in it.
    Empty,
    /// No gene section in any known layout.
    UnknownFormat,
}

/// The gene calls of a report.
pub fn parse_genes(data: &Value) -> Parsed {
    let Value::Object(top) = data else {
        return Parsed::UnknownFormat;
    };
    let mut calls = Vec::new();
    match top.get("genes") {
        Some(Value::Object(genes)) => {
            for (key, val) in genes {
                let Value::Object(o) = val else { continue };
                if o.contains_key("sourceDiplotypes") || o.contains_key("recommendationDiplotypes")
                {
                    // flat (3.x): the key is the gene
                    calls.extend(diplotype_call(key.clone(), val));
                } else {
                    // nested (2.x): the key is the guidance source
                    for (gene, g) in o {
                        calls.extend(diplotype_call(gene.clone(), g));
                    }
                }
            }
        }
        Some(Value::Array(entries)) => {
            for entry in entries {
                let Value::Object(o) = entry else { continue };
                let name = o
                    .get("geneSymbol")
                    .or_else(|| o.get("gene"))
                    .map_or_else(|| "Unknown".to_owned(), text);
                if let Some(c) = diplotype_call(name.clone(), entry) {
                    calls.push(c);
                    continue;
                }
                let dl = o.get("diplotype").map_or_else(|| "N/A".to_owned(), text);
                let ph = o.get("phenotype").map_or_else(|| "N/A".to_owned(), text);
                if dl != "N/A" || ph != "N/A" {
                    let called = dl != "N/A" && dl != "Unknown/Unknown";
                    calls.push(GeneCall::single(name, dl, ph, called));
                }
            }
        }
        _ => match top.get("geneResults") {
            Some(Value::Array(results)) => {
                for gr in results.iter().filter(|g| matches!(g, Value::Object(_))) {
                    let get = |k: &str, d: &str| field(gr, k).map_or_else(|| d.to_owned(), text);
                    let dl = get("diplotype", "N/A");
                    let called = dl != "N/A" && dl != "Unknown/Unknown";
                    calls.push(GeneCall::single(
                        get("gene", "Unknown"),
                        dl,
                        get("phenotype", "N/A"),
                        called,
                    ));
                }
            }
            _ => return Parsed::UnknownFormat,
        },
    }
    let mut out: Vec<GeneCall> = Vec::with_capacity(calls.len());
    for c in calls {
        if !out.iter().any(|o| o.gene == c.gene) {
            out.push(c);
        }
    }
    if out.is_empty() {
        Parsed::Empty
    } else {
        Parsed::Ok(out)
    }
}

/// `{gene: [drug]}` from each gene's `relatedDrugs` (PharmCAT's gene-drug
/// links).
pub fn related_drugs(data: &Value) -> HashMap<String, Vec<String>> {
    let Some(Value::Object(genes)) = field(data, "genes") else {
        return HashMap::new();
    };
    genes
        .iter()
        .filter_map(|(name, g)| {
            let Some(Value::Array(rel)) = field(g, "relatedDrugs") else {
                return None;
            };
            let names: Vec<String> = rel
                .iter()
                .filter(|d| matches!(d, Value::Object(_)))
                .filter_map(|d| truthy_field(d, "name"))
                .map(text)
                .collect();
            (!names.is_empty()).then(|| (name.clone(), names))
        })
        .collect()
}

/// The entry of PharmCAT's `variants` list for `rsid` in `gene` (3.x flat
/// report), or `None` when the report lists none.
pub fn position_read<'a>(data: &'a Value, gene: &str, rsid: &str) -> Option<&'a Value> {
    let g = field(field(data, "genes")?, gene)?;
    let variants = match truthy_field(g, "variants") {
        Some(Value::Array(a)) => a.as_slice(),
        _ => &[],
    };
    variants.iter().find(|v| {
        matches!(v, Value::Object(_))
            && truthy_field(v, "dbSnpId")
                .map(text)
                .unwrap_or_default()
                .split(|c: char| c == ';' || c == ',' || pycompat::text::is_space(c))
                .any(|id| id == rsid)
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use pycompat::json::loads;

    fn calls(json: &str) -> Vec<GeneCall> {
        match parse_genes(&loads(json).unwrap()) {
            Parsed::Ok(c) => c,
            other => panic!("{other:?}"),
        }
    }

    fn dip(a1: &str, a2: &str, phen: &str) -> String {
        format!(
            r#"{{"sourceDiplotypes": [{{"allele1": {{"name": "{a1}"}}, "allele2": {{"name": "{a2}"}}, "label": "{a1}/{a2}", "phenotypes": ["{phen}"]}}]}}"#
        )
    }

    #[test]
    fn every_layout_parses() {
        let g2 = format!(
            r#"{{"CYP2C19": {}, "CYP2D6": {}}}"#,
            dip("*1", "*2", "Intermediate Metabolizer"),
            dip("*1", "*1", "Normal Metabolizer")
        );
        assert_eq!(calls(&format!(r#"{{"genes": {g2}}}"#)).len(), 2);
        assert_eq!(
            calls(&format!(r#"{{"genes": {{"CPIC": {g2}, "DPWG": {g2}}}}}"#)).len(),
            2
        );
        let one = dip("*1", "*2", "Intermediate Metabolizer");
        let listed = format!(
            r#"{{"genes": [{}, "gene": "CYP2C19"}}]}}"#,
            &one[..one.len() - 1]
        );
        assert_eq!(calls(&listed)[0].gene, "CYP2C19");
        let rec = r#"{"genes": {"DPYD": {"recommendationDiplotypes": [{"label": "c.1905+1G>A/Reference", "phenotypes": ["Intermediate Metabolizer"]}]}}}"#;
        let c = &calls(rec)[0];
        assert_eq!(
            (c.diplotype.as_str(), c.status()),
            ("c.1905+1G>A/Reference", Status::NonNormal)
        );
        let old = r#"{"geneResults": [{"gene": "TPMT", "diplotype": "*1/*3A", "phenotype": "Intermediate Metabolizer"}, {"gene": "X", "diplotype": "Unknown/Unknown"}]}"#;
        let c = calls(old);
        assert!(c[0].called && !c[1].called);
        assert_eq!(c[1].phenotype, "N/A");
    }

    #[test]
    fn nothing_parsed_is_not_an_all_clear() {
        assert_eq!(
            parse_genes(&loads(r#"{"metadata": {}}"#).unwrap()),
            Parsed::UnknownFormat
        );
        assert_eq!(parse_genes(&loads("[]").unwrap()), Parsed::UnknownFormat);
        assert_eq!(
            parse_genes(&loads(r#"{"genes": {"CYP2D6": {"geneSymbol": "CYP2D6"}}}"#).unwrap()),
            Parsed::Empty
        );
    }

    #[test]
    fn possible_diplotypes_with_different_phenotypes_are_ambiguous() {
        let c = &calls(
            r#"{"genes": {"CYP2C19": {"sourceDiplotypes": [
              {"allele1": {"name": "*1"}, "allele2": {"name": "*1"}, "label": "*1/*1", "phenotypes": ["Normal Metabolizer"]},
              {"allele1": {"name": "*1"}, "allele2": {"name": "*2"}, "label": "*1/*2", "phenotypes": ["Intermediate Metabolizer"]}]}}}"#,
        )[0];
        assert_eq!(c.status(), Status::Ambiguous);
        assert_eq!(c.diplotype, "*1/*1 (1 of 2 possible)");
        assert_eq!(
            c.phenotype,
            "ambiguous: one of Normal Metabolizer, Intermediate Metabolizer"
        );
        assert_eq!(c.labels, ["*1/*1", "*1/*2"]);
    }

    #[test]
    fn one_diplotype_reads_no_call_names_and_missing_phenotypes() {
        let v = loads(
            r#"{"allele1": {"name": "Unknown"}, "allele2": {"name": " n/a "}, "phenotypes": []}"#,
        )
        .unwrap();
        assert_eq!(
            one_diplotype(&v),
            ("Unknown/ n/a ".to_owned(), "No Result".to_owned(), false)
        );
        let v = loads(r#"{"allele1": {"name": "*1"}, "allele2": null}"#).unwrap();
        assert_eq!(
            one_diplotype(&v),
            ("*1".to_owned(), "no phenotype assigned".to_owned(), true)
        );
        let v = loads(r#"{"label": "Unknown/Unknown", "phenotype": "n/a"}"#).unwrap();
        assert_eq!(
            one_diplotype(&v),
            ("Unknown/Unknown".to_owned(), "No Result".to_owned(), false)
        );
        let v = loads(r#"{"allele1": {"name": null}, "allele2": {}}"#).unwrap();
        // An empty allele object is not a present allele.
        assert_eq!(one_diplotype(&v).0, "?");
    }

    #[test]
    fn activity_scores_are_finite_numbers() {
        for (v, want) in [
            ("1", Some(1.0)),
            ("1.5", Some(1.5)),
            ("\"1.0\"", Some(1.0)),
            ("\"\u{2265}3\"", None),
            ("\"NaN\"", None),
            ("null", None),
            ("true", Some(1.0)),
            ("[1]", None),
        ] {
            assert_eq!(score(&loads(v).unwrap()), want, "{v}");
        }
    }

    #[test]
    fn a_position_is_found_by_its_rsid() {
        let d = loads(
            r#"{"genes": {"CYP4F2": {"variants": [{"dbSnpId": "rs1; rs2108622", "call": "C/T"}, {"dbSnpId": "rs2108622x"}]}}}"#,
        )
        .unwrap();
        assert_eq!(
            position_read(&d, "CYP4F2", "rs2108622").map(|v| text(field(v, "call").unwrap())),
            Some("C/T".to_owned())
        );
        assert!(position_read(&d, "CYP4F2", "rs2").is_none());
        assert!(position_read(&d, "VKORC1", "rs9923231").is_none());
    }
}
