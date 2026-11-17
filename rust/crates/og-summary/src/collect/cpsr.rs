//! Cancer predisposition (step 17, CPSR): the class counts, and CPSR's own
//! class beside the final class and ClinVar's, so a final P/LP class that
//! ClinVar set over CPSR's own benign call is not hidden.

use super::clinvar::{CdnaKey, Floor, cdna_checks, check_word, not_in_table, opt_int, opt_s};
use super::text::{self, cell, digits, stripped};
use indexmap::IndexMap;
use og_bam::cdna::{cpsr_variant, norm_class};
use pycompat::csv::DictRow;
use pycompat::json::Value;
use std::collections::HashMap;

const PLP: [&str; 2] = ["pathogenic", "likely pathogenic"];
const CLASS_COLUMNS: [&str; 3] = [
    "CLASSIFICATION",
    "FINAL_CLASSIFICATION",
    "CPSR_CLASSIFICATION",
];

#[derive(Debug, Clone, PartialEq)]
pub struct Cpsr {
    pub html_report: bool,
    /// The class column and the counts, when the table was read.
    pub table: Option<Table>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Table {
    /// `None`: the table has no class column.
    pub column: Option<String>,
    pub classification: IndexMap<String, usize>,
    pub compare: Option<Compare>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Compare {
    pub plp: Vec<Value>,
    pub disagreements: Vec<Value>,
    pub own_benign: usize,
    pub own: Option<String>,
    pub floor_dp: i128,
    pub cdna: bool,
}

/// Where the inputs beside the CPSR table are.
pub struct Inputs {
    pub html: Option<String>,
    /// The classification table (full path).
    pub tsv: Option<String>,
    /// Step 06's floor, the default without a hits VCF; an error is the
    /// note of a hits VCF that cannot be read.
    pub floor: Result<Floor, String>,
    /// Step 17's cDNA read check (full path).
    pub cdna: Option<String>,
}

impl Cpsr {
    pub fn read(i: &Inputs) -> Result<Cpsr, String> {
        let html_report = i.html.is_some();
        let Some(tsv) = &i.tsv else {
            return Ok(Cpsr {
                html_report,
                table: None,
            });
        };
        let rows = text::rows(tsv).map_err(|e| match e {
            text::ReadError::Os(_) => e.unreadable(tsv),
            text::ReadError::Other(n) => n,
        })?;
        // `rows[0].keys()`: the header, and the None key of a long first row.
        let cols: Vec<&str> = rows
            .first()
            .map(|r| r.iter().map(|(k, _)| k).collect())
            .unwrap_or_default();
        let col = CLASS_COLUMNS.iter().find(|c| cols.contains(c)).copied();
        let Some(col) = col else {
            return Ok(Cpsr {
                html_report,
                table: Some(Table {
                    column: None,
                    classification: IndexMap::new(),
                    compare: None,
                }),
            });
        };
        let mut counts: IndexMap<String, usize> = IndexMap::new();
        for r in &rows {
            let k = cell(r, col).filter(|v| !v.is_empty()).unwrap_or(".");
            *counts.entry(k.to_owned()).or_default() += 1;
        }
        counts.sort_by(|_, a, _, b| b.cmp(a));
        let compare = compare(i, &rows, &cols, col)?;
        Ok(Cpsr {
            html_report,
            table: Some(Table {
                column: Some(col.to_owned()),
                classification: counts,
                compare: Some(compare),
            }),
        })
    }

    pub fn to_json(&self) -> Value {
        let mut m: IndexMap<String, Value> = IndexMap::new();
        m.insert("html_report".into(), Value::Bool(self.html_report));
        let Some(t) = &self.table else {
            return Value::Object(m);
        };
        m.insert("classification_column".into(), opt_s(t.column.as_deref()));
        if t.column.is_some() {
            m.insert(
                "classification".into(),
                Value::Object(
                    t.classification
                        .iter()
                        .map(|(k, v)| (k.clone(), Value::Int(*v as i128)))
                        .collect(),
                ),
            );
        }
        if let Some(c) = &t.compare {
            let field = |x: &Value, k: &str| match x {
                Value::Object(m) => m.get(k).cloned(),
                _ => None,
            };
            let low = c
                .plp
                .iter()
                .filter(|x| field(x, "low_confidence") == Some(Value::Bool(true)))
                .count();
            m.insert("plp".into(), Value::Array(c.plp.clone()));
            m.insert(
                "disagreements".into(),
                Value::Array(c.disagreements.clone()),
            );
            m.insert(
                "final_plp_own_benign".into(),
                Value::Int(c.own_benign as i128),
            );
            m.insert("own_class_column".into(), opt_s(c.own.as_deref()));
            m.insert("low_confidence".into(), Value::Int(low as i128));
            m.insert("call_floor_dp".into(), Value::Int(c.floor_dp));
            if c.cdna {
                let flagged = c
                    .plp
                    .iter()
                    .filter(|x| field(x, "cdna").as_ref().and_then(check_word) == Some("flagged"))
                    .count();
                m.insert("cdna_flagged".into(), Value::Int(flagged as i128));
            }
        }
        Value::Object(m)
    }
}

/// `cpsr_compare`.
fn compare(i: &Inputs, rows: &[DictRow], cols: &[&str], col: &str) -> Result<Compare, String> {
    let own = (cols.contains(&"CPSR_CLASSIFICATION") && col != "CPSR_CLASSIFICATION")
        .then_some("CPSR_CLASSIFICATION");
    let floor = i.floor.clone()?;
    let cdna: Option<HashMap<CdnaKey, Value>> = match &i.cdna {
        Some(p) => Some(cdna_checks(p)?),
        None => None,
    };
    let has_clinvar_col = cols.contains(&"CLINVAR_CLASSIFICATION");
    let (mut plp, mut disagree, mut own_benign) = (Vec::new(), Vec::new(), 0);
    for r in rows {
        let final_ = norm_class(cell(r, col));
        let mine = own.map(|o| norm_class(cell(r, o))).unwrap_or_default();
        let authority = stripped(r, "ASSERTION_AUTHORITY");
        let clinvar_raw = if has_clinvar_col {
            stripped(r, "CLINVAR_CLASSIFICATION")
        } else if authority.to_lowercase() == "clinvar" {
            stripped(r, col)
        } else {
            ""
        };
        let clinvar = norm_class(Some(clinvar_raw));
        let is = |c: &str| PLP.contains(&c);
        if is(&final_) && matches!(mine.as_str(), "benign" | "likely benign") {
            own_benign += 1;
        }
        let is_plp = is(&final_) || is(&mine);
        let differ =
            own.is_some() && !clinvar.is_empty() && !mine.is_empty() && is(&clinvar) != is(&mine);
        if !(is_plp || differ) {
            continue;
        }
        let v = cpsr_variant(cell(r, "VAR_ID"), cell(r, "GENOMIC_CHANGE"))
            .map_err(|e| format!("ValueError: {e}"))?;
        let dp = digits(Some(stripped(r, "DP_CONTROL")));
        let or_dot =
            |x: Option<&str>| Value::Str(x.filter(|v| !v.is_empty()).unwrap_or(".").to_owned());
        let mut x: IndexMap<String, Value> = IndexMap::new();
        x.insert(
            "variant".into(),
            match &v {
                Some(h) => Value::Str(format!("{}:{} {}>{}", h.chrom, h.pos, h.reference, h.alt)),
                None => or_dot(cell(r, "GENOMIC_CHANGE")),
            },
        );
        x.insert("gene".into(), or_dot(cell(r, "SYMBOL")));
        x.insert("final".into(), or_dot(cell(r, col)));
        x.insert(
            "cpsr".into(),
            own.map_or(Value::Null, |o| or_dot(cell(r, o))),
        );
        x.insert(
            "clinvar".into(),
            if clinvar.is_empty() {
                Value::Null
            } else {
                Value::Str(clinvar_raw.to_owned())
            },
        );
        x.insert(
            "authority".into(),
            opt_s(Some(authority).filter(|a| !a.is_empty())),
        );
        x.insert(
            "acmg_codes".into(),
            opt_s(Some(stripped(r, "ACMG_CODE")).filter(|a| !a.is_empty())),
        );
        x.insert("dp".into(), opt_int(dp));
        x.insert(
            "low_confidence".into(),
            Value::Bool(dp.is_some_and(|d| d < floor.dp)),
        );
        if let Some(t) = &cdna
            && is_plp
        {
            let found = v.as_ref().and_then(|h| {
                t.get(&(
                    Some(h.chrom.clone()),
                    i128::from(h.pos),
                    Some(h.reference.clone()),
                    Some(h.alt.clone()),
                ))
            });
            x.insert(
                "cdna".into(),
                found
                    .filter(|c| pycompat::value::truthy(c))
                    .cloned()
                    .unwrap_or_else(not_in_table),
            );
        }
        let x = Value::Object(x);
        if differ {
            disagree.push(x.clone());
        }
        if is_plp {
            plp.push(x);
        }
    }
    Ok(Compare {
        plp,
        disagreements: disagree,
        own_benign,
        own: own.map(str::to_owned),
        floor_dp: floor.dp,
        cdna: cdna.is_some(),
    })
}
