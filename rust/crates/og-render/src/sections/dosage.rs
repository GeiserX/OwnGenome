//! Step 35's gene-dosage loci from read depth (HP, LCE3B/LCE3C,
//! CCL3L1/CCL3L3), one line each in both reports.

use pycompat::json::Value;

use crate::model::{Malformed, is_str, need, py, truthy};

/// The display name of a dosage locus.
fn name(locus: &str) -> &str {
    match locus {
        "HP" => "HP (haptoglobin)",
        "LCE3BC" => "LCE3B/LCE3C",
        "CCL3L1" => "CCL3L1/CCL3L3",
        other => other,
    }
}

/// One line per dosage row, or one line saying why there is none.
pub fn dosage_lines(rows: Option<&Value>) -> Result<Vec<String>, Malformed> {
    let rows = match rows {
        None | Some(Value::Null) => {
            return Ok(vec![
                "Gene-dosage loci (HP, LCE3B/LCE3C, CCL3L1): not assessed (step 35's \
                 <id>_cn_dosage.tsv is not in this run's output)"
                    .to_owned(),
            ]);
        }
        Some(r) if !truthy(r) => {
            return Ok(vec![
                "Gene-dosage loci (HP, LCE3B/LCE3C, CCL3L1): not assessed (step 35's dosage \
                 table has no row)"
                    .to_owned(),
            ]);
        }
        Some(Value::Array(r)) => r,
        Some(_) => {
            return Err(Malformed(
                "sections.parascopy.data.dosage is not a list".into(),
            ));
        }
    };
    let path = "sections.parascopy.data.dosage[]";
    let mut out = Vec::new();
    for r in rows {
        let f = |k: &str| need(r, k, path);
        let label = format!("{} copies (from depth)", name(&py(f("locus")?)?));
        let value = if is_str(Some(f("status")?), "called") {
            let mut v = format!("{} (estimate {})", py(f("copies")?)?, py(f("estimate")?)?);
            if !is_str(Some(f("reading")?), ".") {
                v += &format!(". Reading: {}", py(f("reading")?)?);
            }
            v
        } else {
            format!("{} ({})", py(f("status")?)?, py(f("reason")?)?)
        };
        out.push(format!("{label}: {value}"));
    }
    Ok(out)
}
