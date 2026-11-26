//! Step 35's HBA1 and HBA2 count from read depth (alpha-thalassaemia), one
//! line in both reports.

use pycompat::json::Value;

use crate::model::{Malformed, get, is_str, py};

/// `HBA1 and HBA2 together (alpha-globin genes, from depth): ...`.
pub fn hba_line(h: &Value) -> Result<String, Malformed> {
    let at = |k: &str| get(h, k, "sections.parascopy.data.hba");
    let or = |k: &str, d: &str| at(k).and_then(|v| v.map_or(Ok(d.to_owned()), py));
    let value = if is_str(at("status")?, "called") {
        format!(
            "{} ({}; estimate {})",
            or("alpha_globin_genes", ".")?,
            or("message", ".")?,
            or("estimate", ".")?
        )
    } else {
        or("message", "not assessed")?
    };
    Ok(format!(
        "HBA1 and HBA2 together (alpha-globin genes, from depth): {value}"
    ))
}
