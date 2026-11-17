//! The scan of a site's files for anything that would load or run code
//! from outside the files themselves.

use regex::Regex;
use std::collections::BTreeMap;
use std::path::Path;
use std::sync::LazyLock;

static URL: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r#"(?:https?:)?//[A-Za-z0-9.-]+\.[a-z]{2,}[^\s'"]*"#).unwrap());
static SCRIPT_TAG: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"(?i)<script\b[^>]*>").unwrap());
static ON_ATTR: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"(?i)\son[a-z]+\s*=").unwrap());
static EVAL: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"\beval\s*\(|new\s+Function\s*\(").unwrap());

/// The SVG namespace name: a URL in the text, never fetched.
const SVG_NS: &str = "http://www.w3.org/2000/svg";

/// (file, problem) for each outside URL in any file, each inline
/// `<script>` (one without `src=`) and `on*=` attribute in an `.html` file,
/// and each `eval(` or `new Function(` in a `.js` file. `texts` maps a
/// file's name to its text; the result follows the names' order.
pub fn asset_problems(texts: &BTreeMap<String, String>) -> Vec<(String, String)> {
    let mut out = Vec::new();
    for (name, t) in texts {
        let mut add = |p: &str| out.push((name.clone(), p.to_string()));
        for u in URL
            .find_iter(t)
            .map(|m| m.as_str())
            .filter(|u| *u != SVG_NS)
        {
            add(u);
        }
        if name.ends_with(".html") {
            for tag in SCRIPT_TAG.find_iter(t).map(|m| m.as_str()) {
                if !tag.contains("src=") {
                    add(tag);
                }
            }
            for a in ON_ATTR.find_iter(t) {
                add(a.as_str());
            }
        }
        if name.ends_with(".js") && EVAL.is_match(t) {
            add("eval");
        }
    }
    out
}

/// The text of every file in `dir` and its subfolders, by its path relative
/// to `dir` (`plugins/charts.js`); names that start with `.` are skipped.
pub fn read_assets(dir: &Path) -> std::io::Result<BTreeMap<String, String>> {
    let mut out = BTreeMap::new();
    let mut todo = vec![dir.to_path_buf()];
    while let Some(d) = todo.pop() {
        for e in std::fs::read_dir(&d)? {
            let p = e?.path();
            if p.file_name()
                .is_some_and(|n| n.to_string_lossy().starts_with('.'))
            {
                continue;
            }
            if p.is_dir() {
                todo.push(p);
            } else if p.is_file() {
                let rel = p
                    .strip_prefix(dir)
                    .expect("under dir")
                    .to_string_lossy()
                    .into_owned();
                out.insert(rel, std::fs::read_to_string(&p)?);
            }
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn one(name: &str, text: &str) -> Vec<(String, String)> {
        asset_problems(&BTreeMap::from([(name.to_string(), text.to_string())]))
    }

    #[test]
    fn each_kind_of_problem_is_named() {
        assert_eq!(
            one("a.js", "fetch('https://cdn.example.com/x.js');"),
            [("a.js".into(), "https://cdn.example.com/x.js".into())]
        );
        assert_eq!(one("a.css", "url(//fonts.example.org/f.woff)").len(), 1);
        assert_eq!(
            one(
                "a.html",
                "<SCRIPT>alert(1)</SCRIPT><script src=\"x.js\"></script>"
            ),
            [("a.html".into(), "<SCRIPT>".into())]
        );
        assert_eq!(
            one("a.html", "<a href=x ONCLICK = 'f()'>"),
            [("a.html".into(), " ONCLICK =".into())]
        );
        assert_eq!(
            one("a.js", "var f = new Function ('x');"),
            [("a.js".into(), "eval".into())]
        );
        assert_eq!(
            one("a.js", "window.eval(s)"),
            [("a.js".into(), "eval".into())]
        );
    }

    #[test]
    fn what_is_allowed_is_not_named() {
        assert!(
            one(
                "a.js",
                "document.createElementNS('http://www.w3.org/2000/svg', 'g'); medieval(x)"
            )
            .is_empty()
        );
        assert!(one("a.js", "el.onclick = f; '<script>'").is_empty()); // html rules apply to .html only
        assert!(
            one(
                "a.html",
                "<script src=\"live.js\" defer></script><a href=\"/runs/\">"
            )
            .is_empty()
        );
        assert!(one("a.js", "fetch('/state.json'); var u = 'a//b';").is_empty());
    }
}
