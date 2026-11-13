//! The clinical variant filter (step 23) with its ACMG secondary-findings
//! tier, and the variant prioritization (step 31, slivar) with its
//! compound-het pairs phased from the reads.
//!
//! The bash pipeline writes summary tables; the Nextflow report gets the
//! VCFs instead, so each section reads the table when there is one and
//! counts the VCF's records otherwise.

use super::acmg::{self, ClinicalVcf};
use super::clinvar::obj;
use super::text::{self, cell, digits};
use indexmap::IndexMap;
use pycompat::csv::DictRow;
use pycompat::json::Value;
use std::collections::BTreeSet;

/// The trans pairs the report lists by name; the rest are counted.
pub const COMPHET_TRANS_SHOWN: usize = 20;

fn rows_or_note(path: &str) -> Result<Vec<DictRow>, String> {
    text::rows(path).map_err(|e| e.note())
}

fn int(n: usize) -> Value {
    Value::Int(n as i128)
}

/// The data of the clinical section.
#[derive(Debug, Clone, PartialEq)]
pub struct Clinical(pub Value);

impl Clinical {
    /// `sec_clinical`: `table` is step 23's summary table, `vcf` the
    /// clinical VCF read once (both full paths).
    pub fn read(
        table: Option<&str>,
        vcf: Option<&Result<ClinicalVcf, String>>,
        tier: &acmg::Inputs,
    ) -> Result<Clinical, String> {
        let Some(table) = table else {
            // The Nextflow report: the record count only (one record per
            // variant, as in the table).
            let v = match vcf {
                Some(Ok(v)) => v,
                Some(Err(e)) => return Err(e.clone()),
                None => return Ok(Clinical(Value::Null)),
            };
            return Ok(Clinical(obj([
                ("variants", int(v.records)),
                ("acmg_sf", acmg::acmg_sf(tier)?),
            ])));
        };
        let rows = rows_or_note(table)?;
        let mut by: IndexMap<String, usize> = IndexMap::new();
        for r in &rows {
            let k = cell(r, "IMPACT").filter(|v| !v.is_empty()).unwrap_or(".");
            *by.entry(k.to_owned()).or_default() += 1;
        }
        let genes: BTreeSet<&str> = rows
            .iter()
            .filter_map(|r| cell(r, "GENE"))
            .filter(|g| !matches!(*g, "" | "."))
            .collect();
        let mut m: IndexMap<String, Value> = IndexMap::new();
        m.insert("variants".into(), int(rows.len()));
        m.insert(
            "by_impact".into(),
            Value::Object(by.into_iter().map(|(k, v)| (k, int(v))).collect()),
        );
        m.insert("genes".into(), int(genes.len()));
        m.insert("acmg_sf".into(), acmg::acmg_sf(tier)?);
        // Step 23's tier and call check columns; an older table has neither.
        if let Some(r) = rows.first() {
            if r.get("FREQ_TIER").is_some() {
                let n = rows
                    .iter()
                    .filter(|r| cell(r, "FREQ_TIER") == Some("unknown"))
                    .count();
                m.insert("frequency_unknown".into(), int(n));
            }
            if r.get("CALL_CHECK").is_some() {
                let n = rows
                    .iter()
                    .filter(|r| cell(r, "CALL_CHECK") == Some("weak"))
                    .count();
                m.insert("weak_calls".into(), int(n));
            }
        }
        Ok(Clinical(Value::Object(m)))
    }

    pub fn to_json(&self) -> Value {
        self.0.clone()
    }
}

/// The data of the slivar section.
#[derive(Debug, Clone, PartialEq)]
pub struct Slivar(pub Value);

impl Slivar {
    /// `sec_slivar` with step 31's summary table (`table`) or, for the
    /// Nextflow report, its prioritized VCF (`vcf`); full paths.
    pub fn read(table: Option<&str>, vcf: Option<&str>, s: &str) -> Result<Slivar, String> {
        let mut m: IndexMap<String, Value> = IndexMap::new();
        let Some(p) = table else {
            let Some(v) = vcf else {
                return Ok(Slivar(Value::Null));
            };
            let n = text::vcf_records(v).map_err(|e| e.unreadable(v))?.len();
            m.insert("prioritized".into(), int(n));
            comphet_phase(
                &crate::join(text::dirname(v), &format!("{s}_compound_het_pairs.tsv")),
                &mut m,
            )?;
            return Ok(Slivar(Value::Object(m)));
        };
        let rows = rows_or_note(p)?;
        m.insert("prioritized".into(), int(rows.len()));
        // Step 31's frequency tier: prioritized for the compound-het search, not rare.
        if rows.first().is_some_and(|r| r.get("FREQ_TIER").is_some()) {
            let n = rows
                .iter()
                .filter(|r| cell(r, "FREQ_TIER") == Some("unknown"))
                .count();
            m.insert("frequency_unknown".into(), int(n));
        }
        let dir = text::dirname(p);
        let ch = crate::join(dir, &format!("{s}_compound_hets.tsv"));
        if crate::is_file(&ch) {
            let c = rows_or_note(&ch)?;
            let genes: std::collections::HashSet<&str> = c
                .iter()
                .filter_map(|r| cell(r, "GENE"))
                .filter(|g| !g.is_empty())
                .collect();
            m.insert("compound_het_variants".into(), int(c.len()));
            m.insert("compound_het_genes".into(), int(genes.len()));
        }
        comphet_phase(
            &crate::join(dir, &format!("{s}_compound_het_pairs.tsv")),
            &mut m,
        )?;
        Ok(Slivar(Value::Object(m)))
    }

    pub fn to_json(&self) -> Value {
        self.0.clone()
    }
}

/// Step 31's pairs table: the pairs by phase from the reads, and the trans
/// ones (both copies of the gene hit) by name. A run from before the table
/// adds none of these keys.
fn comphet_phase(path: &str, m: &mut IndexMap<String, Value>) -> Result<(), String> {
    if !crate::is_file(path) {
        return Ok(());
    }
    let rows = rows_or_note(path)?;
    let phase = obj(["trans", "cis", "unknown"].map(|k| {
        (
            k,
            int(rows.iter().filter(|r| cell(r, "PHASE") == Some(k)).count()),
        )
    }));
    let py = |r: &DictRow, k: &str| cell(r, k).unwrap_or("None").to_owned();
    let n = |r: &DictRow, k: &str| digits(cell(r, k)).unwrap_or(0);
    let trans: Vec<Value> = rows
        .iter()
        .filter(|r| cell(r, "PHASE") == Some("trans"))
        .take(COMPHET_TRANS_SHOWN)
        .map(|r| {
            let gene = match r.get("GENE") {
                None => Value::Str(String::new()),
                Some(v) => v.map_or(Value::Null, |g| Value::Str(g.to_owned())),
            };
            obj([
                ("gene", gene),
                (
                    "variants",
                    Value::Str(format!(
                        "{}:{} {}>{} + {}:{} {}>{}",
                        py(r, "CHROM"),
                        py(r, "POS_1"),
                        py(r, "REF_1"),
                        py(r, "ALT_1"),
                        py(r, "CHROM"),
                        py(r, "POS_2"),
                        py(r, "REF_2"),
                        py(r, "ALT_2")
                    )),
                ),
                (
                    "reads",
                    Value::Str(format!(
                        "{} of {}",
                        py(r, "TRANS_READS"),
                        n(r, "TRANS_READS") + n(r, "CIS_READS")
                    )),
                ),
            ])
        })
        .collect();
    m.insert("compound_het_pairs".into(), int(rows.len()));
    m.insert("compound_het_phase".into(), phase);
    m.insert("compound_het_trans".into(), Value::Array(trans));
    Ok(())
}
