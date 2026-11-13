//! The pharmacogenomics sections: PharmCAT (step 07), the CPIC table (step
//! 27), pypgx (step 32) and Cyrius (step 21), and the summary's `cyp2d6`
//! block, which sets the three CYP2D6 callers side by side by step 36's
//! rule (`og_steps::pgx_rules`).

use super::clinvar::{obj, opt_int, opt_s};
use super::model::{Section, State};
use super::sections::Data;
use super::text::{self, cell, stripped};
use indexmap::IndexMap;
use og_steps::pgx_rules::{NO_CALL, is_call, same_diplotype};
use og_steps::pypgx_summary::{LOW_GQ, NOT_CALLABLE};
use pycompat::csv::DictRow;
use pycompat::json::Value;
use pycompat::text::{lower, strip, upper};

/// The label of pypgx's CYP2D6 call when it is the only usable one.
pub const PROBABLE: &str = "probable (one caller: pypgx)";
/// The rows of the CYP2D6 depth check the reports show.
pub const DEPTH_CHECK_KEYS: [&str; 5] = [
    "status",
    "gene_mapq1_share",
    "gene_mapq20_share",
    "body_mapq1_share",
    "body_mapq20_share",
];

fn rows_or_note(path: &str) -> Result<Vec<DictRow>, String> {
    text::rows(path).map_err(|e| e.note())
}

/// `row.get(key, "")`: `""` without the column, `None` for a short row.
fn get_default(r: &DictRow, key: &str) -> Value {
    match r.get(key) {
        None => Value::Str(String::new()),
        Some(v) => opt_s(v),
    }
}

/// `row[key]` as JSON: a `KeyError` note without the column.
fn need(r: &DictRow, key: &str) -> Result<Value, String> {
    text::need(r, key).map(opt_s)
}

// --- PharmCAT ------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq)]
pub struct PharmCat {
    pub version: String,
    pub html_report: bool,
}

/// `json.load` of a file opened as text: universal newlines, and Python's
/// `JSONDecodeError` text (`MSG: line L column C (char N)`, in code points).
fn json_load(raw: &str) -> Result<Value, String> {
    let text = raw.replace("\r\n", "\n").replace('\r', "\n");
    if text.starts_with('\u{feff}') {
        return Err(
            "Unexpected UTF-8 BOM (decode using utf-8-sig): line 1 column 1 (char 0)".into(),
        );
    }
    pycompat::json::loads(&text).map_err(|e| match e {
        pycompat::json::JsonError::Decode(what, at) => {
            let before = &text[..at.min(text.len())];
            let pos = before.chars().count();
            let line = before.matches('\n').count() + 1;
            let col = match before.rfind('\n') {
                Some(i) => before[i + 1..].chars().count() + 1,
                None => pos + 1,
            };
            format!("{what}: line {line} column {col} (char {pos})")
        }
        other => other.to_string(),
    })
}

impl PharmCat {
    pub fn read(path: &str) -> Result<PharmCat, String> {
        let raw = std::fs::read(path)
            .map_err(|e| format!("{}: {}", crate::basename(path), crate::os_error(&e, path)))?;
        let body = String::from_utf8(raw)
            .map_err(|e| format!("{}: 'utf-8' codec can't decode: {e}", crate::basename(path)))?;
        let data = json_load(&body).map_err(|e| format!("{}: {e}", crate::basename(path)))?;
        let version = match &data {
            Value::Object(m) => match m.get("pharmcatVersion") {
                None => String::new(),
                Some(v) => pycompat::value::py_str(v).map_err(|e| e.to_string())?,
            },
            _ => String::new(),
        };
        let html = path.strip_suffix(".json").map(|p| format!("{p}.html"));
        Ok(PharmCat {
            version,
            html_report: html.is_some_and(crate::is_file),
        })
    }

    pub fn to_json(&self) -> Value {
        obj([
            ("version", Value::Str(self.version.clone())),
            ("html_report", Value::Bool(self.html_report)),
        ])
    }
}

// --- CPIC ----------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq)]
pub struct CpicGene {
    pub gene: Value,
    pub diplotype: Value,
    pub phenotype: Value,
    pub status: Value,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Cpic {
    pub genes: Vec<CpicGene>,
    pub parse_failed: bool,
    pub warnings: Vec<String>,
}

impl Cpic {
    /// `sec_cpic`: the phenotype table and, beside it, the recommendation
    /// text's parse failure and warnings.
    pub fn read(path: &str, s: &str) -> Result<Cpic, String> {
        let genes: Vec<CpicGene> = rows_or_note(path)?
            .iter()
            .map(|r| CpicGene {
                gene: get_default(r, "Gene"),
                diplotype: get_default(r, "Diplotype"),
                phenotype: get_default(r, "Phenotype"),
                status: get_default(r, "Status"),
            })
            .collect();
        let rec = crate::join(
            text::dirname(path),
            &format!("{s}_cpic_recommendations.txt"),
        );
        let (mut warnings, mut failed) = (Vec::new(), false);
        if crate::is_file(&rec) {
            for line in text::lines(&rec).map_err(|e| e.note())? {
                if line.contains("PARSING FAILED") {
                    failed = true;
                }
                let t = strip(&line, None);
                if let Some(w) = t.strip_prefix("WARNING:") {
                    warnings.push(strip(w, None).to_owned());
                }
            }
        }
        let parse_failed = failed || genes.is_empty();
        Ok(Cpic {
            genes,
            parse_failed,
            warnings,
        })
    }

    fn count(&self, status: &str) -> Value {
        Value::Int(
            self.genes
                .iter()
                .filter(|g| matches!(&g.status, Value::Str(s) if s == status))
                .count() as i128,
        )
    }

    /// The CYP2D6 row, if the table has one.
    pub fn cyp2d6(&self) -> Option<&CpicGene> {
        self.genes
            .iter()
            .find(|g| matches!(&g.gene, Value::Str(s) if s == "CYP2D6"))
    }

    pub fn to_json(&self) -> Value {
        obj([
            (
                "genes",
                Value::Array(
                    self.genes
                        .iter()
                        .map(|g| {
                            obj([
                                ("gene", g.gene.clone()),
                                ("diplotype", g.diplotype.clone()),
                                ("phenotype", g.phenotype.clone()),
                                ("status", g.status.clone()),
                            ])
                        })
                        .collect(),
                ),
            ),
            ("non_normal", self.count("non-normal")),
            ("ambiguous", self.count("ambiguous")),
            ("not_called", self.count("not called")),
            ("unclassified", self.count("unclassified")),
            ("parse_failed", Value::Bool(self.parse_failed)),
            (
                "warnings",
                Value::Array(self.warnings.iter().cloned().map(Value::Str).collect()),
            ),
        ])
    }
}

// --- pypgx ---------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq)]
pub struct Pypgx {
    pub genes_total: usize,
    pub genes_called: usize,
    pub cyp2d6: String,
    /// Gene (as JSON's key text) to the copy-number call or why none.
    pub cnv_call: IndexMap<String, String>,
    pub not_callable: Vec<Value>,
    pub low_gq: Vec<Value>,
    /// phenotype, defining_filter, quality of the CYP2D6 row.
    pub cyp2d6_evidence: Option<[String; 3]>,
    pub depth_check: Option<IndexMap<String, String>>,
    /// conflicts, one_tool_only, concordant.
    pub comparison: Option<[usize; 3]>,
}

/// A JSON object key from a cell that may be `None` (`json.dumps` writes `null`).
fn key_text(v: Option<&str>) -> String {
    v.unwrap_or("null").to_owned()
}

impl Pypgx {
    pub fn read(path: &str, s: &str) -> Result<Pypgx, String> {
        let rows = rows_or_note(path)?;
        let called = rows
            .iter()
            .filter(|r| !NO_CALL.contains(&lower(stripped(r, "Diplotype")).as_str()))
            .count();
        let cyp_row = rows.iter().find(|r| cell(r, "Gene") == Some("CYP2D6"));
        let cyp = cyp_row.and_then(|r| cell(r, "Diplotype")).unwrap_or("");
        let mut cnv_call = IndexMap::new();
        for r in &rows {
            let status = stripped(r, "Status");
            if upper(stripped(r, "Source")) != "BAM" || status == NOT_CALLABLE {
                continue;
            }
            let v = stripped(r, "CNV_call");
            let gene = key_text(text::need(r, "Gene")?);
            let call = if !matches!(v, "" | ".") && !v.starts_with("N/A") {
                v.to_owned()
            } else if lower(status) == "failed" || upper(stripped(r, "Diplotype")) == "FAILED" {
                "not called: pypgx failed".to_owned()
            } else if status.contains("depth check") || v.contains("depth check") {
                "not called: the CYP2D6 depth check failed".to_owned()
            } else {
                "not called: pypgx gave no copy-number call".to_owned()
            };
            cnv_call.insert(gene, call);
        }
        let mut not_callable = Vec::new();
        for r in rows
            .iter()
            .filter(|r| stripped(r, "Status") == NOT_CALLABLE)
        {
            not_callable.push(need(r, "Gene")?);
        }
        let mut low_gq = Vec::new();
        for r in rows.iter().filter(|r| stripped(r, "Quality") == LOW_GQ) {
            let gene = need(r, "Gene")?;
            let min = stripped(r, "Min_GQ");
            let min_gq = if pycompat::text::isdigit(min) {
                opt_int(crate::int(min).ok())
            } else {
                Value::Null
            };
            low_gq.push(obj([
                ("gene", gene),
                ("diplotype", Value::Str(stripped(r, "Diplotype").to_owned())),
                ("min_gq", min_gq),
                (
                    "weakest",
                    Value::Str(stripped(r, "Weakest_variant").to_owned()),
                ),
            ]));
        }
        let cyp2d6_evidence = cyp_row.map(|r| {
            [
                stripped(r, "Phenotype").to_owned(),
                stripped(r, "Defining_filter").to_owned(),
                stripped(r, "Quality").to_owned(),
            ]
        });
        let dir = text::dirname(path);
        let dc = crate::join(dir, &format!("{s}_cyp2d6_depth_check.tsv"));
        let depth_check = if crate::is_file(&dc) {
            let mut kv: IndexMap<String, String> = IndexMap::new();
            for line in text::lines(&dc).map_err(|e| e.note())? {
                if let Some((k, v)) = line.split_once('\t') {
                    kv.insert(k.to_owned(), v.to_owned());
                }
            }
            Some(
                DEPTH_CHECK_KEYS
                    .iter()
                    .filter_map(|k| kv.get(*k).map(|v| ((*k).to_owned(), v.clone())))
                    .collect(),
            )
        } else {
            None
        };
        let comp = crate::join(dir, &format!("{s}_pharmcat_comparison.tsv"));
        let comparison = if crate::is_file(&comp) {
            let c = rows_or_note(&comp)?;
            let n =
                |f: &dyn Fn(Option<&str>) -> bool| c.iter().filter(|r| f(cell(r, "Match"))).count();
            Some([
                n(&|m| m == Some("No")),
                n(&|m| matches!(m, Some("pypgx only" | "PharmCAT only"))),
                n(&|m| m == Some("Yes")),
            ])
        } else {
            None
        };
        Ok(Pypgx {
            genes_total: rows.len(),
            genes_called: called,
            cyp2d6: if cyp.is_empty() { "not called" } else { cyp }.to_owned(),
            cnv_call,
            not_callable,
            low_gq,
            cyp2d6_evidence,
            depth_check,
            comparison,
        })
    }

    pub fn to_json(&self) -> Value {
        let mut m: IndexMap<String, Value> = IndexMap::new();
        m.insert("genes_total".into(), Value::Int(self.genes_total as i128));
        m.insert("genes_called".into(), Value::Int(self.genes_called as i128));
        m.insert("cyp2d6".into(), Value::Str(self.cyp2d6.clone()));
        m.insert(
            "cnv_call".into(),
            Value::Object(
                self.cnv_call
                    .iter()
                    .map(|(k, v)| (k.clone(), Value::Str(v.clone())))
                    .collect(),
            ),
        );
        m.insert(
            "not_callable".into(),
            Value::Array(self.not_callable.clone()),
        );
        m.insert("low_gq".into(), Value::Array(self.low_gq.clone()));
        if let Some([p, f, q]) = &self.cyp2d6_evidence {
            m.insert(
                "cyp2d6_evidence".into(),
                obj([
                    ("phenotype", Value::Str(p.clone())),
                    ("defining_filter", Value::Str(f.clone())),
                    ("quality", Value::Str(q.clone())),
                ]),
            );
        }
        if let Some(d) = &self.depth_check {
            m.insert(
                "depth_check".into(),
                Value::Object(
                    d.iter()
                        .map(|(k, v)| (k.clone(), Value::Str(v.clone())))
                        .collect(),
                ),
            );
        }
        if let Some([c, o, y]) = self.comparison {
            m.insert(
                "comparison".into(),
                obj([
                    ("conflicts", Value::Int(c as i128)),
                    ("one_tool_only", Value::Int(o as i128)),
                    ("concordant", Value::Int(y as i128)),
                ]),
            );
        }
        Value::Object(m)
    }
}

// --- Cyrius --------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq)]
pub struct Cyrius {
    pub genotype: String,
    /// `None` for a row cut short of the column.
    pub filter: Option<String>,
}

impl Cyrius {
    pub fn read(path: &str) -> Result<Cyrius, String> {
        let rows = rows_or_note(path)?;
        let r = rows.first();
        let genotype = r.and_then(|r| cell(r, "Genotype")).unwrap_or("");
        let filter = match r.map(|r| r.get("Filter")) {
            None | Some(None) => Some(String::new()),
            Some(Some(v)) => v.map(str::to_owned),
        };
        Ok(Cyrius {
            genotype: if genotype.is_empty() {
                "none"
            } else {
                genotype
            }
            .to_owned(),
            filter,
        })
    }

    pub fn to_json(&self) -> Value {
        obj([
            ("genotype", Value::Str(self.genotype.clone())),
            ("filter", opt_s(self.filter.as_deref())),
        ])
    }
}

// --- the cyp2d6 block ----------------------------------------------------------------

/// Step 36's CYP2D6 row, as the summary keeps it.
#[derive(Debug, Clone, PartialEq)]
pub struct Consensus {
    pub result: String,
    pub passed_to_pharmcat: bool,
    pub reason: String,
    pub source: Option<String>,
}

/// pypgx's call alone, labelled.
#[derive(Debug, Clone, PartialEq)]
pub struct Probable {
    pub diplotype: String,
    pub phenotype: String,
}

/// The `cyp2d6` block of summary.json.
#[derive(Debug, Clone, PartialEq)]
pub struct Cyp2d6 {
    pub pharmcat: Value,
    pub pypgx: Value,
    pub cyrius: Value,
    pub agree: Option<bool>,
    pub consensus: Option<Consensus>,
    pub probable: Option<Probable>,
}

/// `probable_cyp2d6`: pypgx's call when it is the only usable one.
fn probable(
    diplotype: &str,
    evidence: Option<&[String; 3]>,
    cyrius: Option<&Cyrius>,
    depth_status: Option<&str>,
) -> Option<Probable> {
    if !is_call(diplotype) {
        return None;
    }
    if let Some(c) = cyrius
        && is_call(&c.genotype)
        && c.filter.as_deref() == Some("PASS")
    {
        return None;
    }
    if depth_status.is_some_and(|d| d != "ok") {
        return None;
    }
    let [phenotype, filter, quality] = evidence.cloned().unwrap_or_default();
    // A row without the evidence passes None, which neither list holds.
    let have = evidence.is_some();
    if !have
        || !matches!(filter.as_str(), "PASS" | "none")
        || !matches!(quality.as_str(), "ok" | "N/A")
    {
        return None;
    }
    Some(Probable {
        diplotype: strip(diplotype, None).to_owned(),
        phenotype: strip(&phenotype, None).to_owned(),
    })
}

/// `read_consensus`: step 36's CYP2D6 row, or `None` when the table is
/// missing, older than a caller's result it was made from, or has no
/// CYP2D6 row.
fn read_consensus(
    found: Option<(String, String)>,
    dir: &str,
    sources: [Option<&str>; 2],
) -> Result<Option<Consensus>, String> {
    let Some((p, rel)) = found else {
        return Ok(None);
    };
    let mtime = |x: &str| pycompat::mtime::mtime(x).map_err(|e| crate::os_error(&e, x));
    let own = mtime(&p)?;
    for src in sources.into_iter().flatten() {
        if mtime(&crate::join(dir, src))? > own {
            return Ok(None);
        }
    }
    let rows = text::rows(&p).map_err(|e| e.text())?;
    let Some(row) = rows.iter().find(|r| stripped(r, "Gene") == "CYP2D6") else {
        return Ok(None);
    };
    Ok(Some(Consensus {
        result: text::cell_or(row, "Result").to_owned(),
        passed_to_pharmcat: cell(row, "Outside_call") == Some("yes"),
        reason: text::cell_or(row, "Reason").to_owned(),
        source: Some(rel),
    }))
}

impl Cyp2d6 {
    /// `cyp2d6_block` over the three caller sections.
    pub fn new(
        cpic: &Section,
        pypgx: &Section,
        cyrius: &Section,
        consensus_at: Option<(String, String)>,
        dir: &str,
    ) -> Cyp2d6 {
        let pc = match &cpic.data {
            Some(Data::Cpic(c)) => c.cyp2d6(),
            _ => None,
        };
        let pd = match (&pypgx.data, pypgx.state) {
            (Some(Data::Pypgx(p)), State::Ok) => Some(p),
            _ => None,
        };
        let cy = match (&cyrius.data, cyrius.state) {
            (Some(Data::Cyrius(c)), State::Ok) => Some(c),
            _ => None,
        };
        let failed = |s: &Section| {
            (s.state == State::Failed).then(|| Value::Str(format!("failed (step {})", s.spec.step)))
        };
        let pg = pd.map(|p| p.cyp2d6.as_str());
        let pg_ok = pg.is_some_and(is_call);
        let cy_ok = cy.is_some_and(|c| is_call(&c.genotype) && c.filter.as_deref() == Some("PASS"));
        let cy_txt = cy.map(|c| {
            let mut t = c.genotype.clone();
            if is_call(&t) && !cy_ok {
                let f = c
                    .filter
                    .as_deref()
                    .filter(|f| !f.is_empty())
                    .unwrap_or("missing");
                t.push_str(&format!(" (not usable: Filter {f})"));
            }
            t
        });
        let pharmcat = match pc {
            Some(g) => {
                let st = &g.status;
                if matches!(st, Value::Str(s) if s == "not called" || s == "ambiguous") {
                    st.clone()
                } else {
                    g.diplotype.clone()
                }
            }
            None => failed(cpic).unwrap_or(Value::Null),
        };
        let agree = match (pg_ok && cy_ok, pg, cy) {
            (true, Some(a), Some(c)) => Some(same_diplotype(a, &c.genotype)),
            _ => None,
        };
        let consensus = match read_consensus(
            consensus_at,
            dir,
            [pypgx.source.as_deref(), cyrius.source.as_deref()],
        ) {
            Ok(c) => c,
            Err(e) => Some(Consensus {
                result: "unreadable".into(),
                passed_to_pharmcat: false,
                reason: e,
                source: None,
            }),
        };
        let probable = match (pg_ok, pd) {
            (true, Some(p)) => probable(
                &p.cyp2d6,
                p.cyp2d6_evidence.as_ref(),
                cy,
                p.depth_check
                    .as_ref()
                    .and_then(|d| d.get("status"))
                    .map(String::as_str),
            ),
            _ => None,
        };
        Cyp2d6 {
            pharmcat,
            pypgx: failed(pypgx).unwrap_or_else(|| opt_s(pg)),
            cyrius: failed(cyrius).unwrap_or_else(|| opt_s(cy_txt.as_deref())),
            agree,
            consensus,
            probable,
        }
    }

    pub fn to_json(&self) -> Value {
        obj([
            (
                "calls",
                obj([
                    ("PharmCAT", self.pharmcat.clone()),
                    ("pypgx", self.pypgx.clone()),
                    ("Cyrius", self.cyrius.clone()),
                ]),
            ),
            ("agree", self.agree.map_or(Value::Null, Value::Bool)),
            (
                "consensus",
                self.consensus.as_ref().map_or(Value::Null, |c| {
                    obj([
                        ("result", Value::Str(c.result.clone())),
                        ("passed_to_pharmcat", Value::Bool(c.passed_to_pharmcat)),
                        ("reason", Value::Str(c.reason.clone())),
                        ("source", opt_s(c.source.as_deref())),
                    ])
                }),
            ),
            (
                "probable",
                self.probable.as_ref().map_or(Value::Null, |p| {
                    obj([
                        ("diplotype", Value::Str(p.diplotype.clone())),
                        ("phenotype", Value::Str(p.phenotype.clone())),
                        ("label", Value::Str(PROBABLE.into())),
                    ])
                }),
            ),
        ])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn json_errors_read_as_python_words_them() {
        assert_eq!(
            json_load("{\"pharmcatVersion\": ").unwrap_err(),
            "Expecting value: line 1 column 21 (char 20)"
        );
        assert_eq!(
            json_load("{\r\n\"é\": x").unwrap_err(),
            "Expecting value: line 2 column 6 (char 7)"
        );
        assert!(
            json_load("\u{feff}{}")
                .unwrap_err()
                .starts_with("Unexpected UTF-8 BOM")
        );
        assert!(json_load("[1]").is_ok());
    }
}
