//! Step 35's CYP21A2 call (congenital adrenal hyperplasia carriers): the
//! copies from depth and the six sites where CYP21A2 differs from its
//! pseudogene CYP21A1P, one line in both reports.

use pycompat::json::Value;

use crate::model::{Malformed, get, is_str, py, truthy};

const LABEL: &str = "CYP21A2 (RCCX module, from depth and six PSV sites)";

/// The CYP21A2 line; not assessed only when the table says so or is not
/// there.
pub fn cyp21a2_line(c: Option<&Value>) -> Result<String, Malformed> {
    let Some(c) = c.filter(|c| truthy(c)) else {
        return Ok(format!(
            "{LABEL}: not assessed (step 35's <id>_cyp21a2_psv.tsv is not in this run's output)"
        ));
    };
    let at = |k: &str| get(c, k, "sections.parascopy.data.cyp21a2");
    let or = |k: &str| at(k).and_then(|v| v.map_or(Ok(".".to_owned()), py));
    if !is_str(at("status")?, "called") {
        let why = match at("message")?
            .filter(|v| truthy(v))
            .or(at("reason")?.filter(|v| truthy(v)))
        {
            Some(v) => py(v)?,
            None => ".".to_owned(),
        };
        let why = why.strip_prefix("not assessed: ").unwrap_or(&why);
        return Ok(format!("{LABEL}: not assessed ({why})"));
    }
    let mut value = format!(
        "{} CYP21A2 copies, module total {} of {}, finding {}",
        or("cyp21a2_copies")?,
        or("module_total")?,
        or("expected_total")?,
        or("finding")?
    );
    if let Some(v) = at("conversion_sites")?.filter(|v| !is_str(Some(v), ".")) {
        value += &format!(" at {}", py(v)?);
    }
    if let Some(v) = at("reading")?.filter(|v| !is_str(Some(v), ".")) {
        value += &format!(". Reading: {}", py(v)?);
    }
    Ok(format!("{LABEL}: {value}"))
}
