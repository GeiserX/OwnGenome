//! Step 9b's motif count of a flagged repeat locus: RFC1's AARRG unit, and
//! any other locus whose motif decides more than its length.

use pycompat::json::Value;

use crate::model::{Malformed, get, py};

/// `<locus> motif count: <verdict>: <note>`, as both reports print it.
pub fn motif_text(locus: &str, m: &Value) -> Result<String, Malformed> {
    let at = |k: &str| get(m, k, "the motif count");
    let verdict = at("verdict")?.map_or(Ok(".".to_owned()), py)?;
    let note = at("note")?.map_or(Ok(String::new()), py)?;
    Ok(format!("{locus} motif count: {verdict}: {note}"))
}

/// RFC1's motif count.
pub fn rfc1_text(m: &Value) -> Result<String, Malformed> {
    motif_text("RFC1", m)
}
