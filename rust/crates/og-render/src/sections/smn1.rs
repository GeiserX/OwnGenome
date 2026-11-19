//! SMN1 read two ways: the copy-aware marker (g.27134T>G, the silent
//! carrier sign, and LPA's KIV-2 4925G>A uses the same reading), and
//! Paraphase's count beside Parascopy's.

use pycompat::json::Value;

use crate::model::{Malformed, get, is_str, need, py, truthy};

/// What the SMN1 marker can and cannot say.
pub const SMN1_MARKER_CAVEAT: &str = "The marker raises the chance that two SMN1 copies sit on one \
    chromosome (2+0, a silent carrier); it does not prove 2+0, and a carrier test in a clinical \
    laboratory confirms it.";

/// A copy-aware marker's state in a few words: present or uncertain with
/// the reads and the copies that carry it, not seen with the reads, too
/// little depth, or not assessed.
pub fn marker_state(m: &Value) -> Result<String, Malformed> {
    let at = |k: &str| get(m, k, "the marker");
    let st = match at("status")?.filter(|v| truthy(v)) {
        Some(v) => py(v)?,
        None => "not assessed".to_owned(),
    };
    if st == "not assessed" {
        return Ok(st);
    }
    let dot = |k: &str| at(k).and_then(|v| v.map_or(Ok(".".to_owned()), py));
    let mut out = format!(
        "{st}, {} of {} reads carry it",
        dot("alt_reads")?,
        dot("depth")?
    );
    if let Some(c) = at("copies")?
        && !matches!(c, Value::Null)
        && !is_str(Some(c), "")
        && !is_str(Some(c), ".")
    {
        out += &format!(", about {} copies", py(c)?);
    }
    Ok(out)
}

/// SMN1 from Paraphase beside Parascopy (step 35), with the verdict; `None`
/// when Parascopy has no SMN1 count.
pub fn smn1_text(m: &Value) -> Result<Option<String>, Malformed> {
    let at = |k: &str| get(m, k, "sections.paralogs.data.smn1");
    let Some(pc) = at("parascopy")?.filter(|v| !matches!(v, Value::Null)) else {
        return Ok(None);
    };
    let pp = match at("paraphase")?.filter(|v| !matches!(v, Value::Null)) {
        Some(v) => py(v)?,
        None => "no call".to_owned(),
    };
    let agree = at("agree")?;
    let verdict = if is_str(agree, "agree") {
        "the two agree"
    } else if is_str(agree, "disagree") {
        "they DISAGREE: check both before using either"
    } else {
        "no comparison: one of them has no call"
    };
    let pc = py(need(pc, "text", "sections.paralogs.data.smn1.parascopy")?)?;
    Ok(Some(format!("Paraphase {pp}, Parascopy {pc}: {verdict}")))
}
