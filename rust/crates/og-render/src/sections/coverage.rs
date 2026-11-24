//! Step 16b's coverage: the mean depth over the genome, the first number of
//! the quality-control card in both reports.

use pycompat::fmt::fixed;
use pycompat::json::Value;

use crate::model::{Malformed, Section, Summary, num};

/// The mean depth as both reports print it (`30.1x`), or `None` when the
/// summary has no number for it.
pub fn mean_depth(s: &Summary) -> Result<Option<String>, Malformed> {
    let md = s.section("coverage").get("mean_depth")?;
    Ok(num(md).map(|x| format!("{}x", fixed(x, 1))))
}

/// The recorded scope of a result that was read.
pub(crate) fn partial_mark<'a>(s: &'a Summary, sec: &Section<'_>) -> Option<&'a str> {
    if !sec.state.ran() {
        return None;
    }
    let Value::Object(fields) = s.need2("sections", sec.key).ok()? else {
        return None;
    };
    match fields.get("partial") {
        Some(Value::Str(mark)) if !mark.is_empty() => Some(mark),
        _ => None,
    }
}

/// The scope marks of the results that carry one, in summary order.
pub(super) fn partial_lines(s: &Summary) -> Vec<String> {
    let Some(Value::Object(sections)) = s.get("sections") else {
        return Vec::new();
    };
    sections
        .iter()
        .filter_map(|(key, section)| {
            let Value::Object(fields) = section else {
                return None;
            };
            let Some(Value::Str(title)) = fields.get("title") else {
                return None;
            };
            partial_mark(s, &s.section(key)).map(|mark| format!("{title}: {mark}"))
        })
        .collect()
}
