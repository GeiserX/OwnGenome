//! What the caller passes on stdin: `{"results": <the results body the page
//! gets>}`, and the person's sex.
//!
//! - Sex: the sex step 16 (indexcov) infers, from `{run}/{run}_sex_check.tsv`,
//!   when it is male or female; else the results body's `person.inferred_sex`,
//!   then `person.declared_sex`; else unknown.
//! - Findings: `results.findings`, one object per finding, in the page's
//!   order. Each keeps the fields the map places and labels it with.
//! - The gated APOE result: the traits card keeps APOE behind its own
//!   warning (`hidden` on the card `traits-panel`). It is not a row of
//!   `results.findings`, so the map adds it as one finding, `apoe`, with
//!   `gate` set, on the APOE gene. Only that the card holds it is read: the
//!   hidden text, which names the genotype, never enters the map.

use pycompat::json::{self, Value};
use std::io::BufRead;
use std::path::Path;

/// Where step 16 writes the sex it infers, under GENOME_DIR. A row of
/// `assets/output_contract.tsv` (the contract-path test checks it).
pub const SEX_PATH: &str = "{sample}/{sample}_sex_check.tsv";

/// The card that keeps APOE behind its warning, and the finding the map
/// makes of it.
pub const GATED_CARD: &str = "traits-panel";
pub const APOE_ID: &str = "apoe";

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Sex {
    Male,
    Female,
    /// Neither measured nor declared: no flag that depends on sex is set.
    Unknown,
}

impl Sex {
    pub fn word(self) -> &'static str {
        match self {
            Sex::Male => "male",
            Sex::Female => "female",
            Sex::Unknown => "unknown",
        }
    }

    fn from_word(s: &str) -> Sex {
        match s.trim().to_ascii_lowercase().as_str() {
            "male" => Sex::Male,
            "female" => Sex::Female,
            _ => Sex::Unknown,
        }
    }
}

/// One finding of the results body, as the map reads it.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Finding {
    pub id: String,
    /// `act`, `medicine`, `review`, `carrier`, `trust` or `info`.
    pub level: String,
    /// Set for a finding the page keeps behind its own warning (`apoe`).
    pub gate: Option<String>,
    /// `chrom:pos` or `chrom:start-end`, 1-based, as the page's reads link.
    pub viewer: Option<String>,
    pub gene: Option<String>,
    pub title: String,
    /// The evidence row `Variant`, when there is one.
    pub hgvs: Option<String>,
    /// The evidence row `Genotype`, when there is one.
    pub genotype: Option<String>,
    pub steps: Vec<String>,
}

/// The findings part of stdin.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Findings {
    /// No results body, or one without `findings`.
    Missing,
    /// `findings` is there but is not a list: the reason.
    Unreadable(String),
    Read(Vec<Finding>),
}

/// The facts the map needs from stdin.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Facts {
    /// The sex the results body gives (the TSV is read apart: [`sex_file`]).
    pub sex: Sex,
    pub findings: Findings,
}

fn text(v: Option<&Value>) -> Option<String> {
    match v {
        Some(Value::Str(s)) if !s.trim().is_empty() => Some(s.clone()),
        _ => None,
    }
}

/// The value of the evidence row named `name`, if any.
fn evidence(f: &Object, name: &str) -> Option<String> {
    let Some(Value::Array(rows)) = f.get("evidence") else {
        return None;
    };
    rows.iter().find_map(|r| match r {
        Value::Array(pair) if pair.len() == 2 && pair[0] == Value::Str(name.to_owned()) => {
            text(pair.get(1))
        }
        _ => None,
    })
}

fn finding(v: &Value) -> Option<Finding> {
    let Value::Object(f) = v else {
        return None;
    };
    let id = text(f.get("id"))?;
    let steps = match f.get("steps") {
        Some(Value::Array(a)) => a
            .iter()
            .filter_map(|s| match s {
                Value::Str(s) => Some(s.clone()),
                Value::Int(i) => Some(i.to_string()),
                _ => None,
            })
            .collect(),
        _ => Vec::new(),
    };
    Some(Finding {
        id,
        level: text(f.get("level")).unwrap_or_else(|| "info".to_owned()),
        gate: text(f.get("gate")),
        viewer: text(f.get("viewer")),
        gene: text(f.get("gene")),
        title: text(f.get("title")).unwrap_or_default(),
        hgvs: evidence(f, "Variant"),
        genotype: evidence(f, "Genotype"),
        steps,
    })
}

type Object = indexmap::IndexMap<String, Value>;

/// The card of the results body that keeps APOE behind its warning, when
/// it holds something there.
fn gated_card(results: &Object) -> Option<&Object> {
    let Some(Value::Array(chapters)) = results.get("chapters") else {
        return None;
    };
    chapters
        .iter()
        .filter_map(|ch| match ch {
            Value::Object(ch) => match ch.get("sections") {
                Some(Value::Array(cards)) => Some(cards),
                _ => None,
            },
            _ => None,
        })
        .flatten()
        .find_map(|c| match c {
            Value::Object(c)
                if c.get("id") == Some(&Value::Str(GATED_CARD.to_owned()))
                    && c.get("hidden").is_some_and(pycompat::value::truthy) =>
            {
                Some(c)
            }
            _ => None,
        })
}

/// The gated APOE finding: placed on the APOE gene, with the card's steps
/// and never a genotype.
fn apoe(card: &Object) -> Finding {
    let steps = match card.get("steps") {
        Some(Value::Array(s)) => s
            .iter()
            .filter_map(|x| match x {
                Value::Str(s) => Some(s.clone()),
                _ => None,
            })
            .collect(),
        _ => Vec::new(),
    };
    Finding {
        id: APOE_ID.to_owned(),
        level: "info".to_owned(),
        gate: Some(APOE_ID.to_owned()),
        gene: Some("APOE".to_owned()),
        title: "APOE".to_owned(),
        steps,
        ..Finding::default()
    }
}

/// Parses the stdin body. Empty input is no facts; anything else must be a
/// JSON object.
pub fn parse(bytes: &[u8]) -> Result<Facts, String> {
    let text = std::str::from_utf8(bytes).map_err(|_| "stdin is not UTF-8".to_owned())?;
    if text.trim().is_empty() {
        return Ok(Facts {
            sex: Sex::Unknown,
            findings: Findings::Missing,
        });
    }
    let v = json::loads(text).map_err(|e| format!("stdin is not JSON: {e}"))?;
    let Value::Object(top) = v else {
        return Err("stdin is not a JSON object".to_owned());
    };
    let results = match top.get("results") {
        Some(Value::Object(r)) => Some(r),
        _ => None,
    };
    let person = results.and_then(|r| r.get("person"));
    let word = |key: &str| match person {
        Some(Value::Object(p)) => match p.get(key) {
            Some(Value::Str(s)) => Sex::from_word(s),
            _ => Sex::Unknown,
        },
        _ => Sex::Unknown,
    };
    let sex = match word("inferred_sex") {
        Sex::Unknown => word("declared_sex"),
        known => known,
    };
    let findings = match results.map(|r| (r, r.get("findings"))) {
        None | Some((_, None)) => Findings::Missing,
        Some((r, Some(Value::Array(list)))) => {
            let mut out: Vec<Finding> = list.iter().filter_map(finding).collect();
            let has_apoe = out
                .iter()
                .any(|f| f.id == APOE_ID || f.gate.as_deref() == Some(APOE_ID));
            if !has_apoe && let Some(card) = gated_card(r) {
                out.push(apoe(card));
            }
            Findings::Read(out)
        }
        Some((_, Some(_))) => Findings::Unreadable("results.findings is not a list".to_owned()),
    };
    Ok(Facts { sex, findings })
}

/// The sex in `{run}/{run}_sex_check.tsv` under `root`, when the file is
/// there and says male or female. A missing or unreadable file, or another
/// word, gives `None`, and the results body decides.
pub fn sex_file(root: &Path, run: &str) -> Option<Sex> {
    let path = root.join(SEX_PATH.replace("{sample}", run));
    let mut lines = og_io::open(&path).ok()?.lines();
    let head = lines.next()?.ok()?;
    let col = head
        .split('\t')
        .position(|h| h.trim_start_matches('#') == "inferred_sex")?;
    let row = lines.next()?.ok()?;
    match Sex::from_word(row.split('\t').nth(col)?) {
        Sex::Unknown => None,
        known => Some(known),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sex_from_reads_then_declared() {
        let p = |s: &str| parse(s.as_bytes()).map(|f| f.sex);
        assert_eq!(p(""), Ok(Sex::Unknown));
        assert_eq!(p("{}"), Ok(Sex::Unknown));
        assert_eq!(
            p(r#"{"results":{"person":{"inferred_sex":"female","declared_sex":"male"}}}"#),
            Ok(Sex::Female)
        );
        assert_eq!(
            p(r#"{"results":{"person":{"inferred_sex":"unknown","declared_sex":"Male"}}}"#),
            Ok(Sex::Male)
        );
        assert_eq!(
            p(r#"{"results":{"person":{"inferred_sex":null,"declared_sex":null}}}"#),
            Ok(Sex::Unknown)
        );
        assert!(p("[1]").is_err());
        assert!(p("{").is_err());
    }

    #[test]
    fn findings_and_the_gated_card() {
        let p = |s: &str| parse(s.as_bytes()).unwrap().findings;
        assert_eq!(p(""), Findings::Missing);
        assert_eq!(p(r#"{"results":{}}"#), Findings::Missing);
        assert!(matches!(
            p(r#"{"results":{"findings":{}}}"#),
            Findings::Unreadable(_)
        ));
        let body = r#"{"results":{"findings":[
            {"id":"sf-x","level":"act","title":"T","gene":"BRCA2","viewer":"chr13:5",
             "steps":["24"],"evidence":[["Variant","c.1A>G"],["Genotype","het"]]},
            {"level":"act"},
            {"id":"q","title":"Q"}],
          "chapters":[{"sections":[{"id":"traits-panel","steps":["42"],
                                    "hidden":{"summary":"S","text":"secret"}}]}]}}"#;
        let Findings::Read(f) = p(body) else {
            panic!("not read")
        };
        // The finding with no id cannot be referenced and is dropped.
        assert_eq!(f.len(), 3);
        assert_eq!(f[0].hgvs.as_deref(), Some("c.1A>G"));
        assert_eq!(f[0].genotype.as_deref(), Some("het"));
        assert_eq!(f[0].steps, ["24"]);
        assert_eq!((f[1].id.as_str(), f[1].level.as_str()), ("q", "info"));
        assert_eq!(f[2].id, "apoe");
        assert_eq!(f[2].gate.as_deref(), Some("apoe"));
        assert_eq!(
            (f[2].gene.as_deref(), f[2].genotype.as_deref()),
            (Some("APOE"), None)
        );
        assert_eq!(f[2].steps, ["42"]);
        // No hidden part on the card: no APOE finding.
        let Findings::Read(f) = p(&body.replace(
            r#""hidden":{"summary":"S","text":"secret"}"#,
            r#""hidden":null"#,
        )) else {
            panic!("not read")
        };
        assert_eq!(f.len(), 2);
    }
}
