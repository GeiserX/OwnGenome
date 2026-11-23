//! Step 16b's coverage: the mean depth over the genome, the first number of
//! the quality-control card in both reports.

use pycompat::fmt::fixed;
use pycompat::json::Value;

use crate::model::{Malformed, Summary, num};

/// The mean depth as both reports print it (`30.1x`), or `None` when the
/// summary has no number for it.
pub fn mean_depth(s: &Summary) -> Result<Option<String>, Malformed> {
    let md = s.section("coverage").get("mean_depth")?;
    Ok(num(md).map(|x| format!("{}x", fixed(x, 1))))
}

/// The scope marks of the results that carry one, in summary order.
pub(super) fn partial_lines(s: &Summary) -> Vec<String> {
    let Some(Value::Object(sections)) = s.get("sections") else {
        return Vec::new();
    };
    sections
        .iter()
        .filter_map(|(key, section)| {
            if !s.section(key).state.ran() {
                return None;
            }
            let Value::Object(fields) = section else {
                return None;
            };
            match (fields.get("title"), fields.get("partial")) {
                (Some(Value::Str(title)), Some(Value::Str(mark))) if !mark.is_empty() => {
                    Some(format!("{title}: {mark}"))
                }
                _ => None,
            }
        })
        .collect()
}
