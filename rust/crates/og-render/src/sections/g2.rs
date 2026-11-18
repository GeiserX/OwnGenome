//! A genotype rule's result (several sites read together), and APOL1's G2
//! beside it: G2 (rs71785313) is read by its own helper, so the rule's
//! G0/G0 or G0/G1 says nothing of it, and the report prints the combined
//! answer of `summary.json` (`traits.rules[APOL1].g2`) under the rule.

use pycompat::json::Value;

use crate::model::{Inline, Malformed, get, is_str, need, py, truthy};

const PATH: &str = "sections.traits.data.rules[]";

/// How a rule's sites were phased, in words; `None` for a phase with no
/// words.
fn phase_text(phase: &Value) -> Option<&'static str> {
    let Value::Str(p) = phase else { return None };
    Some(match p.as_str() {
        "read_backed" => "phased from the reads",
        "not_needed" => "no phase needed",
        "not_used" => "phase not used",
        "unphased" => "not phased",
        _ => return None,
    })
}

/// A rule's result: its label when called, else its status; and how it was
/// phased.
pub fn rule_result(r: &Value) -> Result<String, Malformed> {
    let st = need(r, "status", PATH)?;
    let res = if is_str(Some(st), "called") {
        py(need(r, "result", PATH)?)?
    } else {
        py(st)?.replace('_', " ")
    };
    Ok(match phase_text(need(r, "phase", PATH)?) {
        Some(ph) => format!("{res} ({ph})"),
        None => res,
    })
}

/// The rule's G2 object, when it has one (APOL1 only).
fn g2(r: &Value) -> Result<Option<&Value>, Malformed> {
    Ok(get(r, "g2", PATH)?.filter(|v| truthy(v)))
}

/// The APOL1 rule's result with G2, or `None` for any other rule.
pub fn g2_result(r: &Value) -> Result<Option<String>, Malformed> {
    let Some(g) = g2(r)? else { return Ok(None) };
    let res = need(g, "result", "the rule's g2")?;
    Ok(Some(if is_str(Some(res), ".") {
        "not read".to_owned()
    } else {
        py(res)?
    }))
}

/// The text report's G2 line under the APOL1 rule.
pub fn g2_text(r: &Value) -> Result<Option<String>, Malformed> {
    let Some(g) = g2(r)? else { return Ok(None) };
    let f = |k: &str| need(g, k, "the rule's g2").and_then(py);
    let res = g2_result(r)?.unwrap_or_default();
    Ok(Some(format!(
        "      With G2 (rs71785313): {res}. {} {}",
        f("text")?,
        f("reading")?
    )))
}

/// The G2 result under the rule's result in the HTML table: nothing
/// without one.
pub fn g2_result_html(r: &Value) -> Result<Vec<Inline>, Malformed> {
    Ok(match g2_result(r)?.filter(|x| !x.is_empty()) {
        Some(res) => vec![Inline::Break, Inline::Text(format!("With G2: {res}"))],
        None => Vec::new(),
    })
}

/// The G2 reading and its explanation under the rule's reading in the HTML
/// table: nothing without G2.
pub fn g2_reading_html(r: &Value) -> Result<Vec<Inline>, Malformed> {
    let Some(g) = g2(r)? else {
        return Ok(Vec::new());
    };
    let f = |k: &str| need(g, k, "the rule's g2").and_then(py);
    Ok(vec![
        Inline::Break,
        Inline::Text(format!("With G2: {}", f("reading")?)),
        Inline::Small(f("text")?),
    ])
}
