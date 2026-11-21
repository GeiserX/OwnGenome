//! The live page's search box.
//!
//! Step 24 writes `<run>_search_index.json` beside the report: every gene
//! symbol and rsID the results mention, with the sections that name them.
//! The browser may fetch that file (`/runs/<run>/<run>_search_index.json`,
//! [`crate::plugins::RUN_GLOBS`]), and `state.json` says under
//! `plugins.search` whether it is there and when it was written.
//! `assets/live/plugins/search.js` filters it in the browser, so a query is
//! in no request and no log.

use pycompat::json::Value;
use std::path::Path;

/// `{file, mtime}` (mtime in ns) of the file `<folder name><suffix>` in
/// `run_dir`, or `{}` while there is none.
pub fn file_state(run_dir: &Path, suffix: &str) -> Value {
    let base = run_dir
        .components()
        .next_back()
        .map(|c| c.as_os_str().to_string_lossy().into_owned())
        .unwrap_or_default();
    let name = format!("{base}{suffix}");
    let p = run_dir.join(&name);
    let mut o = indexmap::IndexMap::new();
    if p.is_file()
        && let Ok(ns) = pycompat::mtime::mtime_ns(&p)
    {
        o.insert("file".to_owned(), Value::Str(name));
        o.insert("mtime".to_owned(), Value::Int(ns));
    }
    Value::Object(o)
}

/// `plugins.search` of `state.json`: `{file, mtime}` of the run's
/// `<run>_search_index.json`, or `{}` while there is none.
pub fn state(run_dir: &Path) -> Value {
    file_state(run_dir, "_search_index.json")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_index_file_or_nothing() {
        let d = std::env::temp_dir().join(format!("og-live-search-{}", std::process::id()));
        let run = d.join("S1");
        std::fs::create_dir_all(&run).unwrap();
        assert_eq!(state(&run), Value::Object(indexmap::IndexMap::new()));
        std::fs::write(run.join("S1_search_index.json"), "{}").unwrap();
        let v = state(&run);
        let Value::Object(o) = &v else { panic!() };
        assert_eq!(o["file"], Value::Str("S1_search_index.json".into()));
        assert!(matches!(o["mtime"], Value::Int(n) if n > 0));
        std::fs::remove_dir_all(&d).unwrap();
    }
}
