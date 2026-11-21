//! The live page's charts.
//!
//! Step 24 writes `<run>_charts.json` beside the report. The browser may
//! fetch it (`/runs/<run>/<run>_charts.json`,
//! [`crate::plugins::RUN_GLOBS`]), and `state.json` says under
//! `plugins.charts` whether it is there and when it was written;
//! `assets/live/plugins/charts.js` draws it.

use pycompat::json::Value;
use std::path::Path;

/// `plugins.charts` of `state.json`: `{file, mtime}` (mtime in ns) of the
/// run's `<run>_charts.json`, or `{}` while there is none.
pub fn state(run_dir: &Path) -> Value {
    crate::search::file_state(run_dir, "_charts.json")
}
