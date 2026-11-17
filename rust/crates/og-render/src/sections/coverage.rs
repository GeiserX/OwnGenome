//! Step 16b's coverage: the mean depth over the genome, the first number of
//! the quality-control card in both reports.

use pycompat::fmt::fixed;

use crate::model::{Malformed, Summary, num};

/// The mean depth as both reports print it (`30.1x`), or `None` when the
/// summary has no number for it.
pub fn mean_depth(s: &Summary) -> Result<Option<String>, Malformed> {
    let md = s.section("coverage").get("mean_depth")?;
    Ok(num(md).map(|x| format!("{}x", fixed(x, 1))))
}
