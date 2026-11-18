//! The scan of a site's files for anything that would load or run code
//! from outside the files themselves.

use regex::Regex;
use std::collections::BTreeMap;
use std::path::Path;
use std::sync::LazyLock;

static URL: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r#"(?i)(?:https?:)?//[a-z0-9.-]+\.[a-z]{2,}[^\s'"]*"#).unwrap());
static SCRIPT_TAG: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"(?i)<script\b[^>]*>").unwrap());
/// A `src` attribute of a tag (not `data-src`).
static SRC_ATTR: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"(?i)\ssrc\s*=").unwrap());
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
                if !SRC_ATTR.is_match(tag) {
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
/// A link to a file is read; a link to a folder, a broken link or an entry
/// that cannot be looked at is an error, so no file is left out unseen.
pub fn read_assets(dir: &Path) -> std::io::Result<BTreeMap<String, String>> {
    let mut out = BTreeMap::new();
    let mut todo = vec![dir.to_path_buf()];
    while let Some(d) = todo.pop() {
        for e in std::fs::read_dir(&d)? {
            let e = e?;
            let p = e.path();
            if p.file_name()
                .is_some_and(|n| n.to_string_lossy().starts_with('.'))
            {
                continue;
            }
            let kind = e.file_type()?;
            let real = if kind.is_symlink() {
                std::fs::metadata(&p)?.file_type()
            } else {
                kind
            };
            if real.is_dir() && kind.is_symlink() {
                return Err(std::io::Error::other(format!(
                    "{}: a link to a folder is not scanned",
                    p.display()
                )));
            }
            if real.is_dir() {
                todo.push(p);
            } else if real.is_file() {
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

    #[test]
    fn hosts_and_src_are_read_in_any_case() {
        assert_eq!(
            one("a.js", "fetch('HTTPS://CDN.Example.COM/x.js');"),
            [("a.js".into(), "HTTPS://CDN.Example.COM/x.js".into())]
        );
        assert_eq!(one("a.css", "url(//cdn.example.com/f.woff)").len(), 1);
        // data-src is not src: the script is inline.
        assert_eq!(
            one("a.html", "<script data-src=\"x.js\">f()</script>"),
            [("a.html".into(), "<script data-src=\"x.js\">".into())]
        );
        assert!(one("a.html", "<script SRC = \"x.js\"></script>").is_empty());
    }

    #[test]
    fn read_assets_reads_subfolders_and_refuses_a_linked_folder() {
        let d = crate::Scratch::new(&std::env::temp_dir(), "assets");
        let site = d.path().join("site");
        std::fs::create_dir_all(site.join("plugins")).unwrap();
        std::fs::write(site.join("a.js"), "1").unwrap();
        std::fs::write(site.join("plugins/b.js"), "2").unwrap();
        std::fs::write(site.join(".hidden"), "3").unwrap();
        std::fs::write(d.path().join("outside.js"), "4").unwrap();
        std::os::unix::fs::symlink(d.path().join("outside.js"), site.join("linked.js")).unwrap();
        let got = read_assets(&site).unwrap();
        assert_eq!(
            got.into_iter().collect::<Vec<_>>(),
            [
                ("a.js".to_string(), "1".to_string()),
                ("linked.js".into(), "4".into()),
                ("plugins/b.js".into(), "2".into())
            ]
        );
        std::os::unix::fs::symlink(d.path(), site.join("up")).unwrap();
        let e = read_assets(&site).unwrap_err().to_string();
        assert!(e.contains("a link to a folder is not scanned"), "{e}");
        std::fs::remove_file(site.join("up")).unwrap();
        std::os::unix::fs::symlink(d.path().join("gone.js"), site.join("broken.js")).unwrap();
        assert!(read_assets(&site).is_err());
    }
}
