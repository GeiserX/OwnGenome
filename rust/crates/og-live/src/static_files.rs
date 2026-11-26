//! The live site's static files: the assets under `assets/live` (the
//! scripts, the style and the plugin scripts) and the page templates there,
//! which `og live serve` fills in and never serves as they are.
//!
//! Also the small text rules the pages need: a file's content type, HTML
//! escaping (Python's `html.escape`), `html.unescape` of the report's
//! links, `os.path.normpath` and a query string's values.

use pycompat::text::html_escape;
use std::fs;
use std::path::Path;

/// The page templates of `assets/live`: served filled in, never as an asset.
pub const PAGES: [&str; 5] = [
    "index.html",
    "login.html",
    "home.html",
    "family.html",
    "finding.html",
];

/// What a browser's page shows while there is no run under the root.
pub const NO_RUN: &str = "<li>No run under this folder yet.</li>";

/// The content type of a file by its extension (with the dot, as given:
/// the assets keep its case, a run file is asked with it lowered); `None`
/// for one the site has no type for.
pub fn content_type(ext: &str) -> Option<&'static str> {
    Some(match ext {
        ".html" => "text/html; charset=utf-8",
        ".js" => "text/javascript; charset=utf-8",
        ".css" => "text/css; charset=utf-8",
        ".json" => "application/json",
        ".txt" | ".tsv" | ".md" => "text/plain; charset=utf-8",
        ".svg" => "image/svg+xml",
        ".png" => "image/png",
        ".pdf" => "application/pdf",
        _ => return None,
    })
}

/// The type of an asset or run file: [`content_type`], else bytes.
pub fn type_or_bytes(ext: &str) -> &'static str {
    content_type(ext).unwrap_or("application/octet-stream")
}

/// `os.path.splitext(p)[1]`: the last `.` of the file name and what
/// follows, `""` when the name has none or only starts with dots.
pub fn ext(path: &str) -> &str {
    let name = path.rsplit('/').next().unwrap_or(path);
    let lead = name.len() - name.trim_start_matches('.').len();
    match name[lead..].rfind('.') {
        Some(i) => &name[lead + i..],
        None => "",
    }
}

/// An asset's name as the site serves it: a file of `assets/live`, or of
/// its `plugins/` folder, whose name starts with a letter or digit and
/// holds only letters, digits, `.`, `_` and `-`.
pub fn valid_asset(name: &str) -> bool {
    let rest = name.strip_prefix("plugins/").unwrap_or(name);
    crate::runs::valid_run(rest)
}

/// The bytes and type of asset `name` under `dir`; `None` for a name that
/// is not an asset's, a page template, or a file that is not inside `dir`.
pub fn asset(dir: &Path, name: &str) -> Option<(&'static str, Vec<u8>)> {
    if !valid_asset(name) || PAGES.contains(&name) {
        return None;
    }
    let p = dir.join(name);
    if !p.is_file() || !crate::family::inside(&p, dir) {
        return None;
    }
    let body = fs::read(&p).ok()?;
    Some((type_or_bytes(ext(name)), body))
}

/// `html.escape(text)`: `&`, `<`, `>`, `"` and `'`.
pub fn esc(text: &str) -> String {
    html_escape(text).into_owned()
}

/// The script tags of the plugin scripts in `assets/plugins`, by name.
pub fn plugin_tags(assets: &Path) -> String {
    let mut js: Vec<String> = fs::read_dir(assets.join("plugins"))
        .map(|rd| {
            rd.filter_map(|e| e.ok()?.file_name().into_string().ok())
                .filter(|f| valid_asset(f) && f.ends_with(".js"))
                .collect()
        })
        .unwrap_or_default();
    js.sort();
    js.iter()
        .map(|f| format!("<script src=\"/assets/plugins/{}\" defer></script>", esc(f)))
        .collect::<Vec<_>>()
        .join("\n")
}

/// The page of run `name`: `family.html` for a family (when the assets have
/// it), else `index.html`, filled in for the live site.
pub fn run_page(assets: &Path, name: &str, family: bool) -> std::io::Result<String> {
    let fam = assets.join("family.html");
    let page = if family && fam.is_file() {
        fs::read_to_string(fam)?
    } else {
        fs::read_to_string(assets.join("index.html"))?
    };
    Ok(page
        .replace("{{assets}}", "/assets")
        .replace("{{run}}", &esc(name))
        .replace("{{mode}}", "live")
        .replace("{{plugins}}", &plugin_tags(assets)))
}

/// The home page: `home.html` with the plain list of runs a browser without
/// scripts shows; `runs` is (name, kind) by name.
pub fn home_page(assets: &Path, runs: &[(String, String)]) -> std::io::Result<String> {
    let items = runs
        .iter()
        .map(|(n, kind)| {
            format!(
                "<li><a href=\"/runs/{0}/\">{0}</a> <span class=\"kind\">{1}</span></li>",
                esc(n),
                esc(kind)
            )
        })
        .collect::<Vec<_>>()
        .join("\n");
    let page = fs::read_to_string(assets.join("home.html"))?;
    Ok(page
        .replace("{{assets}}", "/assets")
        .replace("{{data}}", "/runs.json")
        .replace("{{runs}}", if items.is_empty() { NO_RUN } else { &items }))
}

/// The login page with `message` and the address to go back to.
pub fn login_page(assets: &Path, message: &str, next: &str) -> std::io::Result<String> {
    let page = fs::read_to_string(assets.join("login.html"))?;
    Ok(page
        .replace("{{message}}", &esc(message))
        .replace("{{next}}", &esc(next)))
}

/// A finding's page (`finding.html`) for run `name` and finding `id`.
pub fn finding_page(assets: &Path, name: &str, id: &str) -> std::io::Result<String> {
    let page = fs::read_to_string(assets.join("finding.html"))?;
    Ok(page
        .replace("{{assets}}", "/assets")
        .replace("{{run}}", &esc(name))
        .replace("{{finding}}", &esc(id))
        .replace("{{mode}}", "live"))
}

/// `html.unescape` of the references a report's link can hold: the five
/// named ones HTML escaping writes and the numeric ones. Any other `&...;`
/// stays as it is.
pub fn unescape(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    let mut rest = text;
    while let Some(i) = rest.find('&') {
        out.push_str(&rest[..i]);
        rest = &rest[i..];
        let end = rest[1..].find(';').map(|j| j + 1);
        let decoded = end.and_then(|e| {
            let name = &rest[1..e];
            let c = match name {
                "amp" => Some('&'),
                "lt" => Some('<'),
                "gt" => Some('>'),
                "quot" => Some('"'),
                "apos" => Some('\''),
                _ => name.strip_prefix('#').and_then(|n| {
                    let code = match n.strip_prefix(['x', 'X']) {
                        Some(h) => u32::from_str_radix(h, 16).ok(),
                        None => n.parse::<u32>().ok(),
                    };
                    code.and_then(char::from_u32)
                }),
            };
            c.map(|c| (c, e))
        });
        match decoded {
            Some((c, e)) => {
                out.push(c);
                rest = &rest[e + 1..];
            }
            None => {
                out.push('&');
                rest = &rest[1..];
            }
        }
    }
    out.push_str(rest);
    out
}

/// `posixpath.normpath`: `.` and empty parts dropped, `..` taken against
/// the part before it (kept at the start of a relative path), `.` for
/// nothing left.
pub fn normpath(path: &str) -> String {
    if path.is_empty() {
        return ".".into();
    }
    let lead = if path.starts_with("//") && !path.starts_with("///") {
        "//"
    } else if path.starts_with('/') {
        "/"
    } else {
        ""
    };
    let mut parts: Vec<&str> = Vec::new();
    for p in path.split('/') {
        match p {
            "" | "." => {}
            ".." => {
                if !lead.is_empty() && parts.is_empty() {
                    continue;
                }
                if parts.last().is_some_and(|l| *l != "..") {
                    parts.pop();
                } else {
                    parts.push("..");
                }
            }
            _ => parts.push(p),
        }
    }
    let joined = format!("{lead}{}", parts.join("/"));
    if joined.is_empty() {
        ".".into()
    } else {
        joined
    }
}

/// The first value of `key` in a query string, as `urllib.parse.parse_qs`
/// reads it (`+` is a space, `%XX` a byte, an empty value left out).
pub fn query_value(query: &str, key: &str) -> Option<String> {
    query
        .split('&')
        .filter_map(|pair| {
            let (k, v) = pair.split_once('=')?;
            let k = crate::viewer::unquote(&k.replace('+', " "));
            let v = crate::viewer::unquote(&v.replace('+', " "));
            (k == key && !v.is_empty()).then_some(v)
        })
        .next()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn extensions_as_splitext_gives_them() {
        assert_eq!(ext("a/b.tsv"), ".tsv");
        assert_eq!(ext("x.tar.gz"), ".gz");
        assert_eq!(ext(".hidden"), "");
        assert_eq!(ext("dir.d/noext"), "");
        assert_eq!(type_or_bytes(".HTML"), "application/octet-stream");
    }

    #[test]
    fn asset_names() {
        assert!(valid_asset("live.css") && valid_asset("plugins/ops.js"));
        assert!(!valid_asset("../x") && !valid_asset("plugins/") && !valid_asset(".x"));
        assert!(!valid_asset("a/b.js") && !valid_asset("plugins/plugins/a.js"));
    }

    #[test]
    fn normpath_as_python() {
        for (a, b) in [
            ("a/./b", "a/b"),
            ("a/../b", "b"),
            ("../a", "../a"),
            ("a/../../b", "../b"),
            ("/../a", "/a"),
            ("", "."),
            ("a/..", "."),
            ("cpsr/", "cpsr"),
            ("//a", "//a"),
            ("///a", "/a"),
        ] {
            assert_eq!(normpath(a), b, "{a}");
        }
    }

    #[test]
    fn unescape_the_references_of_a_link() {
        assert_eq!(unescape("a&amp;b&#47;c&#x2F;d"), "a&b/c/d");
        assert_eq!(unescape("&nbsp; & &zz"), "&nbsp; & &zz");
    }

    #[test]
    fn query_values() {
        assert_eq!(
            query_value("next=%2Fruns%2FX%2F&a=1", "next").as_deref(),
            Some("/runs/X/")
        );
        assert_eq!(query_value("next=", "next"), None);
        assert_eq!(query_value("a=b+c", "a").as_deref(), Some("b c"));
    }
}
