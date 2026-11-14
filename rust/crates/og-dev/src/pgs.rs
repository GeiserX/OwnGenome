//! `og dev pgs`: the readers of `scripts/ci/check-pgs-labels.sh` that were
//! Python. The script still asks the PGS Catalog with curl; these verbs
//! read its answer on stdin and never touch the network.
//!
//! ```text
//! og dev pgs trait --id PGS000018 < score.json
//! og dev pgs label-match LABEL TRAIT
//! og dev pgs metric-match --id PGS --ppm PPM --ancestry GROUP --effect "OR 1.61" < performance.json
//! ```
//!
//! - `trait` prints the score's `trait_reported`, or exits 1 when the
//!   answer is not JSON, is for another id (the API answers `{}` with HTTP
//!   200 for an unknown one) or has no trait.
//! - `label-match` exits 0 when every word of LABEL is a word of TRAIT,
//!   after lower-casing and dropping apostrophes and punctuation.
//! - `metric-match` prints what the catalog says of the evaluation
//!   (`score, ancestries, effects`) and exits 0 when it belongs to the
//!   score, was measured only in samples of that ancestry group and reports
//!   that odds or hazard ratio; 1 when it does not; 2 when the answer is not
//!   an evaluation with that id.
//!
//! Where it differs from the Python on purpose: a JSON `true` is never the
//! number 1, so `"estimate": true` no longer matches `OR 1`; and an answer
//! whose samples or effects are not the objects and strings the catalog
//! sends is exit 2 (not an evaluation) instead of a traceback read as a
//! mismatch.

use std::io::{Read, Write};
use std::sync::LazyLock;

use pycompat::json::{Value, loads};
use pycompat::value::{py_str, truthy};
use regex::Regex;

const USAGE: &str = "usage: og dev pgs trait --id PGS_ID < ANSWER\n       og dev pgs label-match LABEL TRAIT\n       og dev pgs metric-match --id PGS_ID --ppm PPM_ID --ancestry GROUP --effect 'OR|HR N' < ANSWER\n";

/// The PGS Catalog `ancestry_broad` terms of each 1000 Genomes group
/// pgsc_calc reports (its ancestry_distribution uses the same codes).
const GROUPS: [(&str, &[&str]); 5] = [
    ("EUR", &["European"]),
    (
        "AFR",
        &[
            "African unspecified",
            "African American or Afro-Caribbean",
            "Sub-Saharan African",
        ],
    ),
    ("AMR", &["Hispanic or Latin American"]),
    ("EAS", &["East Asian"]),
    ("SAS", &["South Asian"]),
];

fn group(term: &str) -> Option<&'static str> {
    GROUPS
        .iter()
        .find(|(_, terms)| terms.contains(&term))
        .map(|(g, _)| *g)
}

fn field<'a>(d: &'a Value, key: &str) -> Option<&'a Value> {
    match d {
        Value::Object(o) => o.get(key),
        _ => None,
    }
}

fn is_str(v: Option<&Value>, want: &str) -> bool {
    matches!(v, Some(Value::Str(s)) if s == want)
}

/// The trait of score `id` in `body`, or None.
pub fn trait_of(body: &str, id: &str) -> Option<String> {
    let d = loads(body).ok()?;
    if !matches!(d, Value::Object(_)) || !is_str(field(&d, "id"), id) {
        return None;
    }
    match field(&d, "trait_reported") {
        Some(Value::Str(t)) => {
            let t = pycompat::num::py_strip(t);
            (!t.is_empty()).then(|| t.to_owned())
        }
        _ => None,
    }
}

static NOT_WORD: LazyLock<Regex> =
    LazyLock::new(|| Regex::new("[^a-z0-9]+").expect("a valid pattern"));

fn words(s: &str) -> std::collections::BTreeSet<String> {
    let s = pycompat::text::lower(s).replace(['\u{2019}', '\''], "");
    NOT_WORD
        .replace_all(&s, " ")
        .split_whitespace()
        .map(str::to_owned)
        .collect()
}

/// True when every word of `label` appears in `trait_`.
pub fn label_matches(label: &str, trait_: &str) -> bool {
    let (l, t) = (words(label), words(trait_));
    !l.is_empty() && l.is_subset(&t)
}

static EFFECT: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"^(OR|HR) ([0-9]+(?:\.[0-9]+)?)$").expect("a valid pattern"));

/// What the catalog says of an evaluation, and whether it matches. Err is
/// the reason the answer is not an evaluation with that id.
pub fn metric(
    body: &str,
    id: &str,
    ppm: &str,
    anc: &str,
    eff: &str,
) -> Result<(String, bool), String> {
    let d = loads(body).map_err(|_| "(not JSON)".to_owned())?;
    if !matches!(d, Value::Object(_)) || !is_str(field(&d, "id"), ppm) {
        return Err("(no such evaluation)".into());
    }
    let odd = || "(not an evaluation the check can read)".to_owned();
    let empty = Value::Null;
    let or_empty =
        |v: Option<&'_ Value>| -> Value { v.filter(|v| truthy(v)).cloned().unwrap_or(Value::Null) };
    let samples = match or_empty(field(&or_empty(field(&d, "sampleset")), "samples")) {
        Value::Null => Vec::new(),
        Value::Array(a) => a,
        _ => return Err(odd()),
    };
    let mut terms = Vec::new();
    for s in &samples {
        if !matches!(s, Value::Object(_)) {
            return Err(odd());
        }
        match or_empty(field(s, "ancestry_broad")) {
            Value::Null => terms.push(String::new()),
            Value::Str(t) => terms.push(t),
            _ => return Err(odd()),
        }
    }
    let found: std::collections::BTreeSet<Option<&str>> = terms.iter().map(|t| group(t)).collect();
    let sizes = match or_empty(field(
        &or_empty(field(&d, "performance_metrics")),
        "effect_sizes",
    )) {
        Value::Null => Vec::new(),
        Value::Array(a) => a,
        _ => return Err(odd()),
    };
    let mut effects = Vec::new();
    for e in &sizes {
        if !matches!(e, Value::Object(_)) {
            return Err(odd());
        }
        let name = field(e, "name_short").unwrap_or(&empty).clone();
        let est = field(e, "estimate").unwrap_or(&empty).clone();
        effects.push((name, est));
    }
    let pgs = field(&d, "associated_pgs_id").unwrap_or(&empty);
    let show = |v: &Value| py_str(v).map_err(|_| odd());
    let effects_said = effects
        .iter()
        .map(|(n, v)| Ok(format!("{} {}", show(n)?, show(v)?)))
        .collect::<Result<Vec<_>, String>>()?
        .join("; ");
    let terms_said = terms.join("; ");
    let said = format!(
        "{}, {}, {}",
        show(pgs)?,
        if terms_said.is_empty() {
            "no samples"
        } else {
            &terms_said
        },
        if effects_said.is_empty() {
            "no effect size"
        } else {
            &effects_said
        },
    );
    let ok = is_str(Some(pgs), id)
        && group_named(anc)
        && found.len() == 1
        && found.contains(&Some(anc))
        && EFFECT.captures(eff).is_some_and(|m| {
            let want: f64 = m[2].parse().unwrap_or(f64::NAN);
            effects.iter().any(|(n, v)| {
                is_str(Some(n), &m[1])
                    && match v {
                        Value::Int(i) => (*i as f64 - want).abs() < 1e-9,
                        Value::Float(x) => (x - want).abs() < 1e-9,
                        _ => false, // a bool is never a number
                    }
            })
        });
    Ok((said, ok))
}

fn group_named(anc: &str) -> bool {
    GROUPS.iter().any(|(g, _)| *g == anc)
}

fn stdin() -> Result<String, String> {
    let mut s = String::new();
    std::io::stdin()
        .read_to_string(&mut s)
        .map_err(|e| e.to_string())?;
    Ok(s)
}

/// `--name value` pairs.
fn opts<'a>(args: &'a [String], names: &[&str]) -> Option<Vec<&'a str>> {
    let mut got = vec![None; names.len()];
    let mut it = args.iter();
    while let Some(a) = it.next() {
        let i = names.iter().position(|n| n == a)?;
        got[i] = Some(it.next()?.as_str());
    }
    got.into_iter().collect()
}

/// `og dev pgs ARGS`.
pub fn cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let usage = |err: &mut dyn Write| {
        let _ = write!(err, "{USAGE}");
        2
    };
    match args.first().map(String::as_str) {
        Some("trait") => {
            let Some(v) = opts(&args[1..], &["--id"]) else {
                return usage(err);
            };
            match stdin().ok().and_then(|b| trait_of(&b, v[0])) {
                Some(t) => {
                    let _ = writeln!(out, "{t}");
                    0
                }
                None => 1,
            }
        }
        Some("label-match") => match &args[1..] {
            [label, trait_] => u8::from(!label_matches(label, trait_)),
            _ => usage(err),
        },
        Some("metric-match") => {
            let Some(v) = opts(&args[1..], &["--id", "--ppm", "--ancestry", "--effect"]) else {
                return usage(err);
            };
            let body = match stdin() {
                Ok(b) => b,
                Err(e) => {
                    let _ = writeln!(err, "ERROR: {e}");
                    return 2;
                }
            };
            match metric(&body, v[0], v[1], v[2], v[3]) {
                Ok((said, ok)) => {
                    let _ = writeln!(out, "{said}");
                    u8::from(!ok)
                }
                Err(why) => {
                    let _ = writeln!(out, "{why}");
                    2
                }
            }
        }
        _ => usage(err),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn trait_needs_the_same_id_and_a_trait() {
        assert_eq!(
            trait_of(
                r#"{"id": "PGS000018", "trait_reported": " Coronary artery disease "}"#,
                "PGS000018"
            )
            .as_deref(),
            Some("Coronary artery disease")
        );
        assert_eq!(trait_of("{}", "PGS000018"), None);
        assert_eq!(
            trait_of(r#"{"id": "PGS000020", "trait_reported": "x"}"#, "PGS000018"),
            None
        );
        assert_eq!(
            trait_of(r#"{"id": "PGS000018", "trait_reported": ""}"#, "PGS000018"),
            None
        );
        assert_eq!(trait_of("not json", "PGS000018"), None);
    }

    #[test]
    fn labels_match_by_words() {
        assert!(label_matches(
            "Type 2 diabetes (T2D)",
            "Type 2 diabetes (T2D)"
        ));
        assert!(label_matches("type_2_diabetes", "Type 2 diabetes (T2D)"));
        assert!(label_matches(
            "Late-onset Alzheimer\u{2019}s disease",
            "Late-onset Alzheimer's disease"
        ));
        assert!(!label_matches(
            "inflammatory_bowel_disease",
            "Type 2 diabetes (T2D)"
        ));
        assert!(!label_matches("", "anything"));
    }

    const ANSWER: &str = r#"{"id": "PPM000591", "associated_pgs_id": "PGS000018",
        "sampleset": {"samples": [{"ancestry_broad": "European"}]},
        "performance_metrics": {"effect_sizes": [{"name_short": "HR", "estimate": 1.53}]}}"#;

    #[test]
    fn an_evaluation_matches_its_score_ancestry_and_effect() {
        let (said, ok) = metric(ANSWER, "PGS000018", "PPM000591", "EUR", "HR 1.53").unwrap();
        assert!(ok);
        assert_eq!(said, "PGS000018, European, HR 1.53");
        for (id, anc, eff) in [
            ("PGS000004", "EUR", "HR 1.53"),
            ("PGS000018", "AFR", "HR 1.53"),
            ("PGS000018", "EUR", "HR 1.54"),
            ("PGS000018", "EUR", "OR 1.53"),
            ("PGS000018", "XYZ", "HR 1.53"),
            ("PGS000018", "EUR", "HR"),
        ] {
            assert!(
                !metric(ANSWER, id, "PPM000591", anc, eff).unwrap().1,
                "{id} {anc} {eff}"
            );
        }
        assert_eq!(
            metric("{}", "PGS000018", "PPM000591", "EUR", "HR 1.53").unwrap_err(),
            "(no such evaluation)"
        );
        assert_eq!(
            metric("<html>", "PGS000018", "PPM000591", "EUR", "HR 1.53").unwrap_err(),
            "(not JSON)"
        );
    }

    #[test]
    fn a_bool_is_never_a_number() {
        let body = ANSWER.replace("1.53", "true");
        let (said, ok) = metric(&body, "PGS000018", "PPM000591", "EUR", "HR 1").unwrap();
        assert!(!ok, "{said}");
        assert_eq!(said, "PGS000018, European, HR True");
        let body = ANSWER.replace("1.53", "2");
        assert!(
            metric(&body, "PGS000018", "PPM000591", "EUR", "HR 2")
                .unwrap()
                .1
        );
    }

    #[test]
    fn no_samples_no_effects() {
        let body = r#"{"id": "PPM1", "associated_pgs_id": null, "sampleset": null}"#;
        let (said, ok) = metric(body, "PGS1", "PPM1", "EUR", "OR 1").unwrap();
        assert!(!ok);
        assert_eq!(said, "None, no samples, no effect size");
    }
}
