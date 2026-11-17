//! `og dev site-check SITE [ZENSICAL_TOML]`: the post-build check of the
//! Zensical docs site, for the rules `zensical build --strict` does not
//! enforce. It replaces `.github/scripts/site-check.py` with the same checks
//! and the same ERROR lines.
//!
//! Run it from the repository root (the Docs workflow does): the config is
//! read from `./zensical.toml` unless its path is given. No config file is
//! an error, never a skipped pass.
//!
//! Exit 1 on any of these, each printed as one `ERROR` line:
//! - a link or image whose target is not in the built site: a nav entry
//!   whose page does not exist, a missing image, any other dangling href or
//!   src, in any quoting and with HTML entities decoded. A full URL into the
//!   site itself (`https://geiserx.github.io/<repo>/...`) is checked the same
//!   way. A redirect stub's own links (canonical and body) are checked too;
//! - a link whose `#fragment` is not an id on the target page, including
//!   links that come in through a snippet include;
//! - an href ending in `.md` in a page body: a correct page link is
//!   rewritten to a directory URL, so a `.md` left in the HTML is a raw-HTML
//!   or unprocessed link;
//! - a root-absolute link (`/getting-started/`) outside the site's own path
//!   prefix. The prefix comes from `site_url`, else the home page's canonical
//!   link; with neither, that is an error rather than a check that cannot
//!   fire;
//! - a page that no navigation entry links to (a stray page);
//! - a `.py`, `.pyc` or `__pycache__` entry in the site;
//! - an awesome list's count token left in a page (`<!-- n:proyectos -->`,
//!   `n:categorias`, `n:retirados`);
//! - a raw `<title>`, or an unclosed `<textarea>`, `<xmp>`, `<noembed>`,
//!   `<noframes>` (a self-closed `<textarea/>` counts as unclosed: HTML
//!   ignores that slash), or any `<plaintext>`, in a page's article outside
//!   inline SVG and MathML: a browser shows nothing after it;
//! - a "Copy as Markdown" button (`data-md-copy-url`) whose Markdown copy
//!   does not exist in the site;
//! - a link that cannot be resolved to a file (a `%00` in its path);
//! - an `index.html` that is not a content page (no `md-content` block, or
//!   a redirect stub);
//! - a `site_url`, or the canonical link standing in for it, that is not an
//!   absolute http(s) URL.
//!
//! And, from `zensical.toml` (Zensical drops what it does not know without
//! a word):
//! - a plugin name Zensical 0.1.0 does not know;
//! - a theme feature name the Zensical UI 0.0.41 does not know;
//! - `theme.variant` other than `modern` or `classic`, and a modern config
//!   whose `index.html` does not load the modern stylesheet;
//! - a `redirect_maps` entry whose stub was not written, a `redirect_maps`
//!   that is not a table, and `theme.features` that is not a list of
//!   strings;
//! - an enabled llmstxt plugin with no `site/llms.txt`, or
//!   `content.action.copy` in the features with a page that carries no copy
//!   button.
//!
//! Plugins are read in every form Zensical accepts (a table, or a list of
//! names and one-entry tables; `material/redirects` is `redirects`), and a
//! plugin with `enabled = false` is off. With directory URLs off (the offline
//! plugin turns them off) every `.html` file is a page.
//!
//! The HTML is read the way Python's `html.parser` (3.14) reads it: any
//! quoting, entities decoded in attribute values, nothing taken from
//! `<script>`/`<style>` text or comments. Numeric character references are
//! all decoded; of the named ones, only `&amp;`, `&lt;`, `&gt;`, `&quot;`,
//! `&apos;` and `&nbsp;` are, the ones a Markdown build writes. Any other
//! named reference stays as written, so a link or id that holds one is
//! reported, never passed in silence.

use indexmap::IndexMap;
use pycompat::json::Value;
use pycompat::text::{lower, split, strip};
use regex::Regex;
use std::collections::{HashMap, HashSet};
use std::io::Write;
use std::path::{Component, Path, PathBuf};
use std::rc::Rc;
use std::sync::LazyLock;

const USAGE: &str = "usage: og dev site-check SITE_DIR [ZENSICAL_TOML]";

/// Tags whose raw use in an article hides the rest of the page.
const SWALLOW: [&str; 6] = [
    "title",
    "textarea",
    "xmp",
    "plaintext",
    "noembed",
    "noframes",
];
/// HTML ignores the slash in `<x/>` except on these and in SVG/MathML:
/// `<textarea/>` opens a textarea.
const VOID: [&str; 14] = [
    "area", "base", "br", "col", "embed", "hr", "img", "input", "link", "meta", "param", "source",
    "track", "wbr",
];
/// The names of an awesome list's count tokens.
const TOKENS: [&str; 3] = ["n:proyectos", "n:categorias", "n:retirados"];
/// The `<title>` of a Zensical redirect stub.
const STUB: &str = "Redirecting...";

/// Keys of `_PLUGIN_UNSUPPORTED_OPTIONS` in zensical/config.py (0.1.0): every
/// other name is dropped by Zensical without a message.
const KNOWN_PLUGINS: [&str; 31] = [
    "api-autonav",
    "autorefs",
    "awesome-nav",
    "blog",
    "callouts",
    "exclude",
    "gh-admonitions",
    "glightbox",
    "include-markdown",
    "literate-nav",
    "llmstxt",
    "macros",
    "markdown-exec",
    "markdownextradata",
    "meta",
    "mike",
    "minify",
    "mkdocs-audio",
    "mkdocs-autoapi",
    "mkdocs-nav-weight",
    "mkdocs-video",
    "mkdocstrings",
    "offline",
    "optimize",
    "privacy",
    "redirects",
    "rss",
    "search",
    "social",
    "table-reader",
    "tags",
];

/// Feature names the shipped UI tests for (zensical_ui 0.0.41, read from its
/// templates and `bundle.*.min.js`). Regenerate after a bump with:
///   grep -rhoE '"[a-z]+(\.[a-z]+)+" in features' <site-packages>/zensical_ui/templates | sort -u
///   grep -oE '"(navigation|content|search|toc|header|announce)\.[a-z.]+"' .../bundle.*.min.js | sort -u
const KNOWN_FEATURES: [&str; 29] = [
    "announce.dismiss",
    "content.action.copy",
    "content.action.edit",
    "content.action.view",
    "content.code.annotate",
    "content.code.copy",
    "content.code.select",
    "content.footnote.tooltips",
    "content.lazy",
    "content.tabs.link",
    "content.tooltips",
    "header.autohide",
    "navigation.expand",
    "navigation.footer",
    "navigation.indexes",
    "navigation.instant",
    "navigation.instant.prefetch",
    "navigation.instant.preview",
    "navigation.instant.progress",
    "navigation.path",
    "navigation.prune",
    "navigation.sections",
    "navigation.tabs",
    "navigation.tabs.sticky",
    "navigation.top",
    "navigation.tracking",
    "search.highlight",
    "toc.follow",
    "toc.integrate",
];

/// `og dev site-check SITE [ZENSICAL_TOML]`: 0 when the site passes, 1 on
/// any ERROR line, 2 on a usage error.
pub fn cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    if args.first().is_some_and(|a| a == "-h" || a == "--help") {
        let _ = writeln!(out, "{USAGE}");
        return 0;
    }
    if args.is_empty() || args.len() > 2 {
        let _ = writeln!(err, "{USAGE}");
        return 2;
    }
    let config = args.get(1).map_or("zensical.toml", String::as_str);
    check(&args[0], config, out)
}

/// Runs every check over the built site and config, prints one ERROR line
/// each and the summary line; returns 1 if any check failed.
pub fn check(site_dir: &str, config_path: &str, out: &mut dyn Write) -> u8 {
    let site = resolve_path(Path::new(site_dir));
    let index = site.join("index.html");
    if !index.is_file() {
        let _ = writeln!(out, "ERROR {site_dir}: no index.html (did the build run?)");
        return 1;
    }
    let mut errors: HashSet<String> = HashSet::new();
    let cfg = read_config(Path::new(config_path), &mut errors);
    let plugins = plugin_table(cfg.get("plugins"), &mut errors);
    let theme = match cfg.get("theme") {
        Some(Value::Str(name)) => {
            let mut t = IndexMap::new();
            t.insert("name".to_owned(), Value::Str(name.clone()));
            t
        }
        Some(Value::Object(t)) => t.clone(),
        _ => IndexMap::new(),
    };
    // The offline plugin turns directory URLs off (zensical/config.py,
    // _apply_offline_plugin).
    let directory_urls = !matches!(cfg.get("use_directory_urls"), Some(Value::Bool(false)))
        && !plugin_on(&plugins, "offline");

    let mut pages_cache = Pages::default();
    let configured = cfg.get("site_url").filter(|v| truthy(v));
    let from_canonical = configured.is_none();
    let mut site_url = match configured {
        Some(v) => py_str(v),
        None => pages_cache
            .get(&index)
            .canonical
            .clone()
            .unwrap_or_default(),
    };
    let origin = match urlparse(&site_url).and_then(|u| u.port().map(|p| (u, p))) {
        Ok(found) => Some(found),
        Err(()) => {
            errors.insert(format!(
                "site_url {} is not a valid URL",
                repr_str(&site_url)
            ));
            site_url.clear();
            None
        }
    };
    let (origin, origin_port) = origin.unwrap_or_default();
    if !site_url.is_empty()
        && (!matches!(origin.scheme.as_str(), "http" | "https") || origin.hostname().is_none())
    {
        errors.insert(format!(
            "site_url {} is not an absolute http(s) URL",
            repr_str(&site_url)
        ));
    }
    let host = lower(&origin.hostname().unwrap_or_default());
    let mut prefix = unquote(&origin.path);
    if from_canonical && prefix.ends_with(".html") {
        // A canonical .../index.html: its directory.
        prefix = prefix
            .rsplit_once('/')
            .map_or(prefix.clone(), |(d, _)| d.to_owned());
    }
    if site_url.is_empty() {
        errors.insert(
            "no site_url in zensical.toml and no canonical link in index.html: the \
             root-absolute link check has no prefix to check against"
                .to_owned(),
        );
        prefix = "/".to_owned();
    }
    if !prefix.ends_with('/') {
        prefix.push('/');
    }
    let ctx = Ctx {
        site: site.clone(),
        host,
        origin_port,
        prefix: prefix.clone(),
    };

    for p in walk(&site) {
        let rel = p.strip_prefix(&site).unwrap_or(&p).to_path_buf();
        let suffix = py_suffix(&p);
        let cache_dir = rel.components().any(|c| c.as_os_str() == "__pycache__");
        if suffix == ".py" || suffix == ".pyc" || cache_dir {
            errors.insert(format!(
                "{}: source file published (add its directory to exclude_docs)",
                rel.display()
            ));
        }
    }

    let mut nav_targets: HashSet<PathBuf> = HashSet::new();
    let mut pages: Vec<PathBuf> = Vec::new();
    let features: Vec<String> = match theme.get("features").filter(|v| truthy(v)) {
        None => Vec::new(),
        Some(Value::Array(items)) if items.iter().all(|f| matches!(f, Value::Str(_))) => items
            .iter()
            .filter_map(|f| match f {
                Value::Str(s) => Some(s.clone()),
                _ => None,
            })
            .collect(),
        Some(other) => {
            errors.insert(format!(
                "zensical.toml: theme.features must be a list of strings, got {}",
                repr(other)
            ));
            Vec::new()
        }
    };
    let copy_on = features.iter().any(|f| f == "content.action.copy");
    let mut html_files: Vec<PathBuf> = walk(&site)
        .into_iter()
        .filter(|p| {
            !p.is_dir()
                && p.file_name()
                    .is_some_and(|n| n.to_string_lossy().ends_with(".html"))
        })
        .collect();
    html_files.sort();
    for page in html_files {
        let rel = page.strip_prefix(&site).unwrap_or(&page).to_path_buf();
        let rel_s = rel.display().to_string();
        let doc = pages_cache.get(&page);
        let named = if directory_urls {
            rel.file_name().is_some_and(|n| n == "index.html")
        } else {
            rel_s != "404.html"
        };
        let is_page = named && doc.content && !doc.stub;
        if is_page {
            pages.push(page.clone());
        }
        for tag in &doc.swallowed {
            errors.insert(format!(
                "{rel_s}: raw <{tag}> in the page body hides everything after it in a browser"
            ));
        }
        for tok in &doc.tokens {
            errors.insert(format!(
                "{rel_s}: build token {tok} was not filled in (extension not loaded?)"
            ));
        }
        for (tag, url, in_head) in &doc.urls {
            match ctx.resolve(&page, url) {
                Target::Outside => {}
                Target::Bad => {
                    errors.insert(format!("{rel_s}: {} is not a valid URL", repr_str(url)));
                }
                Target::Abs => {
                    errors.insert(format!("{rel_s}: absolute link {url} (outside {prefix})"));
                }
                Target::File(t) if !t.starts_with(&site) => {
                    // ../LICENSE from the top page: the file exists in the
                    // checkout, not on Pages.
                    errors.insert(format!(
                        "{rel_s}: {url} leaves the built site (resolves to {})",
                        t.display()
                    ));
                }
                Target::File(t) if !t.exists() => {
                    errors.insert(format!(
                        "{rel_s}: {url} -> missing {}",
                        t.strip_prefix(&site).unwrap_or(&t).display()
                    ));
                }
                Target::File(t) => {
                    // <head> holds <link rel="alternate" type="text/markdown"
                    // href="./index.md"> from the llmstxt plugin, a
                    // legitimate .md href; this rule reads the body only.
                    if !in_head && unquote(&url_path(url)).ends_with(".md") {
                        errors.insert(format!(
                            "{rel_s}: {url} links a Markdown source, not a page URL (raw HTML link?)"
                        ));
                    }
                    // A text directive (#id:~:text=...) is not part of the
                    // id; #top needs no id. A redirect stub's ids are not
                    // checked: the browser lands on the stub's target.
                    let fragment = urlparse(strip(url, None))
                        .map(|u| u.fragment)
                        .unwrap_or_default();
                    let fragment = unquote(&fragment);
                    let frag = fragment.split(":~:").next().unwrap_or_default();
                    if matches!(tag.as_str(), "a" | "area")
                        && !frag.is_empty()
                        && lower(frag) != "top"
                        && py_suffix(&t) == ".html"
                    {
                        let target = pages_cache.get(&t);
                        if !target.stub && !target.ids.contains(frag) {
                            errors.insert(format!(
                                "{rel_s}: {url} -> no id {} on {}",
                                repr_str(frag),
                                t.strip_prefix(&site).unwrap_or(&t).display()
                            ));
                        }
                    }
                }
            }
        }
        for url in &doc.copies {
            let fine = match ctx.resolve(&page, url) {
                Target::File(t) => t.starts_with(&site) && py_suffix(&t) == ".md" && t.is_file(),
                _ => false,
            };
            if !fine {
                errors.insert(format!(
                    "{rel_s}: Copy as Markdown points at {}, which is not a Markdown copy in the site (llmstxt sections?)",
                    repr_str(url)
                ));
            }
        }
        if is_page && copy_on && doc.copies.is_empty() {
            errors.insert(format!(
                "{rel_s}: content.action.copy is on but the page has no Copy as Markdown button"
            ));
        }
        for url in &doc.nav {
            if let Target::File(t) = ctx.resolve(&page, url) {
                nav_targets.insert(t);
            }
        }
    }

    if !pages.contains(&index) {
        errors.insert(
            "index.html is not a content page (no md-content block, or a redirect stub)".to_owned(),
        );
    }
    for page in &pages {
        if !nav_targets.contains(&resolve_path(page)) {
            errors.insert(format!(
                "{}: page is in no navigation entry (stray page)",
                page.strip_prefix(&site).unwrap_or(page).display()
            ));
        }
    }

    if !cfg.is_empty() {
        for name in plugins.keys() {
            if !KNOWN_PLUGINS.contains(&canonical_plugin(name)) {
                errors.insert(format!(
                    "zensical.toml: plugin {} is unknown to Zensical and silently ignored",
                    repr_str(name)
                ));
            }
        }
        for f in &features {
            if !KNOWN_FEATURES.contains(&f.as_str()) {
                errors.insert(format!(
                    "zensical.toml: feature {} is unknown to the Zensical UI and silently ignored",
                    repr_str(f)
                ));
            }
        }
        // Zensical's default (config.py, set_default).
        let modern = Value::Str("modern".to_owned());
        let variant = theme.get("variant").unwrap_or(&modern);
        let variant_name = match variant {
            Value::Str(s) if s == "modern" || s == "classic" => Some(s.as_str()),
            _ => None,
        };
        match variant_name {
            None => {
                errors.insert(format!(
                    "zensical.toml: theme.variant {} is not modern or classic (Zensical falls back to classic silently)",
                    repr(variant)
                ));
            }
            Some("modern")
                if !pages_cache
                    .get(&index)
                    .styles
                    .iter()
                    .any(|h| h.contains("assets/stylesheets/modern/")) =>
            {
                errors.insert(
                    "index.html does not load the modern stylesheet although theme.variant is modern"
                        .to_owned(),
                );
            }
            Some(_) => {}
        }
        let empty = Value::Object(IndexMap::new());
        let redirects = if plugin_on(&plugins, "redirects") {
            plugin_opts(&plugins, "redirects")
        } else {
            None
        };
        let maps = redirects
            .and_then(|r| r.get("redirect_maps"))
            .filter(|v| truthy(v))
            .unwrap_or(&empty);
        let maps = match maps {
            Value::Object(m) => m.clone(),
            other => {
                errors.insert(format!(
                    "zensical.toml: redirects redirect_maps must be a table, got {}",
                    repr(other)
                ));
                IndexMap::new()
            }
        };
        for src in maps.keys() {
            let trimmed = src.trim_start_matches('/');
            let mut s = trimmed.strip_suffix(".md").unwrap_or(trimmed).to_owned();
            // A README page is its folder's index: sub/README.md ->
            // sub/index.html in every mode.
            if s == "README" || s.ends_with("/README") {
                s.truncate(s.len() - "README".len());
                s.push_str("index");
            }
            let stub = if !directory_urls || s == "index" || s.ends_with("/index") {
                resolve_path(&site.join(format!("{s}.html")))
            } else {
                resolve_path(&site.join(&s).join("index.html"))
            };
            let shown = stub.strip_prefix(&site).unwrap_or(&stub).to_path_buf();
            if !stub.is_file() || !pages_cache.get(&stub).stub {
                errors.insert(format!(
                    "redirect_maps {}: no redirect stub at {}",
                    repr_str(src),
                    shown.display()
                ));
            }
        }
        if plugin_on(&plugins, "llmstxt") && !site.join("llms.txt").is_file() {
            errors.insert(
                "llmstxt plugin is configured but site/llms.txt was not written".to_owned(),
            );
        }
    }

    let mut sorted: Vec<String> = errors.into_iter().collect();
    sorted.sort();
    for e in &sorted {
        let _ = writeln!(out, "ERROR {e}");
    }
    let _ = writeln!(
        out,
        "site-check: {} pages, {} errors, prefix {prefix}",
        pages.len(),
        sorted.len()
    );
    u8::from(!sorted.is_empty())
}

/// Each file is parsed once: fragment checks read target pages again.
#[derive(Default)]
struct Pages(HashMap<PathBuf, Rc<Page>>);

impl Pages {
    fn get(&mut self, path: &Path) -> Rc<Page> {
        self.0
            .entry(path.to_path_buf())
            .or_insert_with(|| {
                // read_text(errors="replace") in text mode: newlines translated.
                let text = std::fs::read(path)
                    .map(|b| {
                        String::from_utf8_lossy(&b)
                            .replace("\r\n", "\n")
                            .replace('\r', "\n")
                    })
                    .unwrap_or_default();
                Rc::new(Page::parse(&text))
            })
            .clone()
    }
}

/// Where a URL lands.
enum Target {
    /// Outside this site (another host, mailto:, data:, another repo's site).
    Outside,
    /// A rooted path outside the site's prefix.
    Abs,
    /// A URL that cannot be parsed or resolved.
    Bad,
    /// The file the URL lands on, resolved.
    File(PathBuf),
}

struct Ctx {
    site: PathBuf,
    host: String,
    origin_port: Option<u16>,
    prefix: String,
}

impl Ctx {
    fn resolve(&self, page: &Path, url: &str) -> Target {
        let url = strip(url, None);
        let Ok(u) = urlparse(url) else {
            return Target::Bad;
        };
        let Ok(port) = u.port() else {
            return Target::Bad;
        };
        if !matches!(u.scheme.as_str(), "" | "http" | "https") {
            return Target::Outside; // mailto:, data:, javascript:
        }
        let path;
        if !u.netloc.is_empty() {
            // https://host/... or //host/...: in this site only on its own
            // host and prefix.
            let same_port =
                port.is_none() || port == self.origin_port || matches!(port, Some(80 | 443));
            if lower(&u.hostname().unwrap_or_default()) != self.host || !same_port {
                return Target::Outside;
            }
            let p = unquote(&u.path);
            path = if p.is_empty() { "/".to_owned() } else { p };
            if !path.starts_with(&self.prefix) && format!("{path}/") != self.prefix {
                // Another site on the same host (geiserx.github.io/<other repo>).
                return Target::Outside;
            }
        } else if !u.scheme.is_empty() {
            return Target::Outside;
        } else {
            path = unquote(&u.path);
            if path.is_empty() {
                return Target::File(page.to_path_buf()); // "", "#x", "?q#x": the page itself
            }
        }
        if path.contains('\0') {
            return Target::Bad;
        }
        let mut target = if path.starts_with('/') {
            if !path.starts_with(&self.prefix) && format!("{path}/") != self.prefix {
                return Target::Abs;
            }
            let rest = path.strip_prefix(&self.prefix).unwrap_or("");
            self.site.join(rest)
        } else {
            page.parent().unwrap_or(Path::new("/")).join(&path)
        };
        if path.ends_with('/') || target.is_dir() {
            target = target.join("index.html");
        }
        Target::File(resolve_path(&target))
    }
}

/// Python's `Path.resolve()` (strict=False): an absolute path with every
/// symlink followed and `.` and `..` removed; parts that do not exist are
/// kept as written.
fn resolve_path(p: &Path) -> PathBuf {
    let abs = if p.is_absolute() {
        p.to_path_buf()
    } else {
        std::env::current_dir()
            .unwrap_or_else(|_| PathBuf::from("/"))
            .join(p)
    };
    let mut queue: std::collections::VecDeque<std::ffi::OsString> =
        std::collections::VecDeque::new();
    let push_all = |queue: &mut std::collections::VecDeque<std::ffi::OsString>, p: &Path| {
        let parts: Vec<std::ffi::OsString> = p
            .components()
            .map(|c| match c {
                Component::RootDir => std::ffi::OsString::from("/"),
                other => other.as_os_str().to_owned(),
            })
            .collect();
        for part in parts.into_iter().rev() {
            queue.push_front(part);
        }
    };
    push_all(&mut queue, &abs);
    let mut out = PathBuf::from("/");
    let mut links = 0;
    while let Some(part) = queue.pop_front() {
        if part == "/" {
            out = PathBuf::from("/");
        } else if part == "." {
        } else if part == ".." {
            out.pop();
        } else {
            let candidate = out.join(&part);
            let is_link =
                std::fs::symlink_metadata(&candidate).is_ok_and(|m| m.file_type().is_symlink());
            if is_link && links < 40 {
                links += 1;
                match std::fs::read_link(&candidate) {
                    Ok(target) => push_all(&mut queue, &target),
                    Err(_) => out = candidate,
                }
            } else {
                out = candidate;
            }
        }
    }
    out
}

/// Every entry under `root`, recursively, without entering symlinked
/// directories (Python's `rglob("*")`).
fn walk(root: &Path) -> Vec<PathBuf> {
    let mut out = Vec::new();
    let mut stack = vec![root.to_path_buf()];
    while let Some(dir) = stack.pop() {
        let Ok(entries) = std::fs::read_dir(&dir) else {
            continue;
        };
        for entry in entries.flatten() {
            let path = entry.path();
            if entry.file_type().is_ok_and(|t| t.is_dir()) {
                stack.push(path.clone());
            }
            out.push(path);
        }
    }
    out
}

/// Python's `PurePath.suffix`: the name's last dotted part, `""` for none
/// or for a leading dot only.
fn py_suffix(p: &Path) -> String {
    let name = p
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default();
    match name.rfind('.') {
        Some(i) if i > 0 && i < name.len() - 1 => name[i..].to_owned(),
        _ => String::new(),
    }
}

// ---------------------------------------------------------------- config

/// The `[project]` table of zensical.toml; an unreadable or empty config is
/// an error.
fn read_config(path: &Path, errors: &mut HashSet<String>) -> IndexMap<String, Value> {
    let shown = path.display();
    if !path.is_file() {
        errors.insert(format!(
            "{shown}: config file not found (run from the repository root, or pass its path)"
        ));
        return IndexMap::new();
    }
    let parsed = std::fs::read(path)
        .map_err(|e| e.to_string())
        .and_then(|b| String::from_utf8(b).map_err(|e| e.to_string()))
        .and_then(|text| toml::parse(&text));
    let doc = match parsed {
        Ok(doc) => doc,
        Err(e) => {
            errors.insert(format!("{shown}: not valid TOML: {e}"));
            return IndexMap::new();
        }
    };
    // Every config check reads [project]: without it they would all pass by
    // skipping.
    match doc.get("project") {
        Some(Value::Object(p)) if !p.is_empty() => p.clone(),
        _ => {
            errors.insert(format!(
                "{shown}: no [project] table, so no config check can run"
            ));
            IndexMap::new()
        }
    }
}

/// `plugins` in any form Zensical accepts (a table, or a list of names and
/// one-entry tables) as one table keyed by the names as written.
fn plugin_table(raw: Option<&Value>, errors: &mut HashSet<String>) -> IndexMap<String, Value> {
    match raw {
        None => IndexMap::new(),
        Some(Value::Object(t)) => t.clone(),
        Some(Value::Array(items)) => {
            let mut table = IndexMap::new();
            for item in items {
                match item {
                    Value::Str(name) => {
                        table.insert(name.clone(), Value::Object(IndexMap::new()));
                    }
                    Value::Object(one) if one.len() == 1 => {
                        for (k, v) in one {
                            table.insert(k.clone(), v.clone());
                        }
                    }
                    other => {
                        errors.insert(format!(
                            "zensical.toml: plugins entry {} is neither a name nor a one-entry table",
                            repr(other)
                        ));
                    }
                }
            }
            table
        }
        Some(other) => {
            errors.insert(format!(
                "zensical.toml: plugins is a {}, not a table or a list",
                type_name(other)
            ));
            IndexMap::new()
        }
    }
}

/// The plugin name Zensical resolves a config key to (`material/search` is
/// `search`).
fn canonical_plugin(name: &str) -> &str {
    let name = name.strip_prefix("material/").unwrap_or(name);
    if ["social/", "privacy/", "optimize/"]
        .iter()
        .any(|p| name.starts_with(p))
    {
        name.split('/').next().unwrap_or(name)
    } else {
        name
    }
}

/// True when the plugin is configured under any of its key spellings and
/// not `enabled = false`.
fn plugin_on(plugins: &IndexMap<String, Value>, name: &str) -> bool {
    for (key, opts) in plugins {
        if canonical_plugin(key) == name {
            return !matches!(opts, Value::Object(o) if matches!(o.get("enabled"), Some(Value::Bool(false))));
        }
    }
    false
}

/// The plugin's option table, under whichever key spelling the config uses.
fn plugin_opts<'a>(
    plugins: &'a IndexMap<String, Value>,
    name: &str,
) -> Option<&'a IndexMap<String, Value>> {
    plugins.iter().find_map(|(key, opts)| match opts {
        Value::Object(o) if canonical_plugin(key) == name => Some(o),
        _ => None,
    })
}

/// Python's truth value of a TOML value.
fn truthy(v: &Value) -> bool {
    match v {
        Value::Null => false,
        Value::Bool(b) => *b,
        Value::Int(i) => *i != 0,
        Value::Float(f) => *f != 0.0,
        Value::Str(s) => !s.is_empty(),
        Value::Array(a) => !a.is_empty(),
        Value::Object(o) => !o.is_empty(),
    }
}

/// Python's `type(v).__name__`.
fn type_name(v: &Value) -> &'static str {
    match v {
        Value::Null => "NoneType",
        Value::Bool(_) => "bool",
        Value::Int(_) => "int",
        Value::Float(_) => "float",
        Value::Str(_) => "str",
        Value::Array(_) => "list",
        Value::Object(_) => "dict",
    }
}

/// `repr(v)`, recursing into lists and tables so every string inside goes
/// through [`write_str_repr`]. Scalars other than strings never fail in
/// `py_repr`.
fn repr(v: &Value) -> String {
    let mut out = String::new();
    write_repr(&mut out, v);
    out
}

fn write_repr(out: &mut String, v: &Value) {
    match v {
        Value::Str(s) => write_str_repr(out, s),
        Value::Array(a) => {
            out.push('[');
            for (i, x) in a.iter().enumerate() {
                if i > 0 {
                    out.push_str(", ");
                }
                write_repr(out, x);
            }
            out.push(']');
        }
        Value::Object(o) => {
            out.push('{');
            for (i, (k, x)) in o.iter().enumerate() {
                if i > 0 {
                    out.push_str(", ");
                }
                write_str_repr(out, k);
                out.push_str(": ");
                write_repr(out, x);
            }
            out.push('}');
        }
        scalar => out.push_str(&pycompat::value::py_repr(scalar).unwrap_or_default()),
    }
}

/// The characters Python's `str.isprintable()` rejects above U+00FF: the
/// general categories Cc, Cf, Co, Cn, Zl, Zp and Zs (Cs cannot occur in a
/// Rust `char`).
static UNPRINTABLE: LazyLock<Regex> = LazyLock::new(|| {
    Regex::new(r"^[\p{Cc}\p{Cf}\p{Co}\p{Cn}\p{Zl}\p{Zp}\p{Zs}]$").expect("static regex")
});

/// Python's `repr(s)` of a string: single quotes unless the text has a `'`
/// and no `"`; `\\`, the quote, `\t`, `\n`, `\r` escaped; other
/// unprintable characters as `\xNN` below U+0100, `\uNNNN` below U+10000
/// and `\UNNNNNNNN` above; printable characters as written.
fn write_str_repr(out: &mut String, s: &str) {
    let quote = if s.contains('\'') && !s.contains('"') {
        '"'
    } else {
        '\''
    };
    out.push(quote);
    let mut buf = [0u8; 4];
    for c in s.chars() {
        let code = u32::from(c);
        match c {
            '\\' => out.push_str("\\\\"),
            '\t' => out.push_str("\\t"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            c if c == quote => {
                out.push('\\');
                out.push(c);
            }
            '\0'..='\x1f' | '\x7f'..='\u{a0}' | '\u{ad}' => {
                out.push_str(&format!("\\x{code:02x}"));
            }
            '\u{100}'.. if UNPRINTABLE.is_match(c.encode_utf8(&mut buf)) => {
                if code < 0x1_0000 {
                    out.push_str(&format!("\\u{code:04x}"));
                } else {
                    out.push_str(&format!("\\U{code:08x}"));
                }
            }
            c => out.push(c),
        }
    }
    out.push(quote);
}

/// `repr(s)` of a string.
fn repr_str(s: &str) -> String {
    repr(&Value::Str(s.to_owned()))
}

/// `str(v)`.
fn py_str(v: &Value) -> String {
    match v {
        Value::Str(s) => s.clone(),
        other => repr(other),
    }
}

// ------------------------------------------------------------------- URLs

/// The parts of Python's `urlparse` the checks read.
#[derive(Default)]
struct Url {
    scheme: String,
    netloc: String,
    path: String,
    fragment: String,
}

impl Url {
    /// `(hostname, port text)` as Python's `_hostinfo` splits the netloc.
    fn hostinfo(&self) -> (&str, Option<&str>) {
        let hostinfo = self
            .netloc
            .rsplit_once('@')
            .map_or(self.netloc.as_str(), |(_, h)| h);
        let (hostname, port) = match hostinfo.split_once('[') {
            Some((_, bracketed)) => {
                let (h, after) = bracketed.split_once(']').unwrap_or((bracketed, ""));
                (h, after.split_once(':').map_or("", |(_, p)| p))
            }
            None => hostinfo.split_once(':').unwrap_or((hostinfo, "")),
        };
        (hostname, (!port.is_empty()).then_some(port))
    }

    /// `.hostname`: lowercased, `None` when empty (an IPv6 zone keeps its case).
    fn hostname(&self) -> Option<String> {
        let (h, _) = self.hostinfo();
        if h.is_empty() {
            return None;
        }
        Some(match h.split_once('%') {
            Some((addr, zone)) => format!("{}%{zone}", lower(addr)),
            None => lower(h),
        })
    }

    /// `.port`: an error where Python raises ValueError.
    fn port(&self) -> Result<Option<u16>, ()> {
        match self.hostinfo().1 {
            None => Ok(None),
            Some(p) if !p.is_empty() && p.bytes().all(|b| b.is_ascii_digit()) => {
                let digits = p.trim_start_matches('0');
                if digits.len() > 5 {
                    return Err(());
                }
                digits
                    .parse::<u32>()
                    .map_or(Ok(Some(0)), |n| u16::try_from(n).map(Some).map_err(|_| ()))
            }
            Some(_) => Err(()),
        }
    }
}

/// Python's `urllib.parse.urlparse(url)`: an error where it raises
/// ValueError (a bad bracketed host).
fn urlparse(url: &str) -> Result<Url, ()> {
    let mut url: String = url
        .trim_start_matches(|c: char| c <= ' ')
        .chars()
        .filter(|c| !matches!(c, '\t' | '\r' | '\n'))
        .collect();
    let mut scheme = String::new();
    if let Some(i) = url.find(':') {
        let head = &url[..i];
        if i > 0
            && head.as_bytes()[0].is_ascii_alphabetic()
            && head
                .bytes()
                .all(|b| b.is_ascii_alphanumeric() || matches!(b, b'+' | b'-' | b'.'))
        {
            scheme = head.to_ascii_lowercase();
            url = url[i + 1..].to_owned();
        }
    }
    let mut netloc = String::new();
    if url.starts_with("//") {
        let end = url[2..].find(['/', '?', '#']).map_or(url.len(), |e| e + 2);
        netloc = url[2..end].to_owned();
        url = url[end..].to_owned();
        let open = netloc.contains('[');
        let close = netloc.contains(']');
        if open != close {
            return Err(());
        }
        if open {
            check_bracketed_netloc(&netloc)?;
        }
    }
    let mut fragment = String::new();
    if let Some((u, f)) = url.split_once('#') {
        fragment = f.to_owned();
        url = u.to_owned();
    }
    if let Some((u, _query)) = url.split_once('?') {
        url = u.to_owned();
    }
    // urlparse also splits ;params off the last path segment.
    const USES_PARAMS: [&str; 16] = [
        "", "ftp", "hdl", "prospero", "http", "imap", "https", "shttp", "rtsp", "rtsps", "rtspu",
        "sip", "sips", "mms", "sftp", "tel",
    ];
    if USES_PARAMS.contains(&scheme.as_str()) && url.contains(';') {
        let cut = match url.rfind('/') {
            Some(slash) => url[slash..].find(';').map(|i| i + slash),
            None => url.find(';'),
        };
        if let Some(i) = cut {
            url.truncate(i);
        }
    }
    Ok(Url {
        scheme,
        netloc,
        path: url,
        fragment,
    })
}

/// Python's `_check_bracketed_netloc`: the host inside `[...]` must be an
/// IPv6 address or an IPvFuture literal.
fn check_bracketed_netloc(netloc: &str) -> Result<(), ()> {
    let host_port = netloc.rsplit_once('@').map_or(netloc, |(_, h)| h);
    let host = match host_port.split_once('[') {
        Some((before, bracketed)) => {
            if !before.is_empty() {
                return Err(());
            }
            let (h, port) = bracketed.split_once(']').unwrap_or((bracketed, ""));
            if !port.is_empty() && !port.starts_with(':') {
                return Err(());
            }
            h
        }
        None => host_port.split_once(':').map_or(host_port, |(h, _)| h),
    };
    if host.starts_with(['v', 'V']) {
        let rest = &host[1..];
        let (hex, tail) = rest.split_once('.').ok_or(())?;
        if hex.is_empty() || !hex.bytes().all(|b| b.is_ascii_hexdigit()) || tail.is_empty() {
            return Err(());
        }
        return Ok(());
    }
    let addr = match host.split_once('%') {
        Some((a, zone)) if !zone.is_empty() && !zone.contains('/') => a,
        Some(_) => return Err(()),
        None => host,
    };
    addr.parse::<std::net::Ipv6Addr>()
        .map(|_| ())
        .map_err(|_| ())
}

/// The path of a URL, or `""` when there is none or the URL cannot be
/// parsed.
fn url_path(url: &str) -> String {
    urlparse(strip(url, None))
        .map(|u| u.path)
        .unwrap_or_default()
}

/// Python's `urllib.parse.unquote`: `%XX` escapes decoded as UTF-8, invalid
/// bytes replaced.
fn unquote(s: &str) -> String {
    if !s.contains('%') {
        return s.to_owned();
    }
    let b = s.as_bytes();
    let mut bytes = Vec::with_capacity(b.len());
    let mut i = 0;
    while i < b.len() {
        if b[i] == b'%'
            && i + 2 < b.len()
            && b[i + 1].is_ascii_hexdigit()
            && b[i + 2].is_ascii_hexdigit()
        {
            let hex = std::str::from_utf8(&b[i + 1..i + 3]).unwrap_or("00");
            bytes.push(u8::from_str_radix(hex, 16).unwrap_or(0));
            i += 3;
        } else {
            bytes.push(b[i]);
            i += 1;
        }
    }
    String::from_utf8_lossy(&bytes).into_owned()
}

// ------------------------------------------------------------------- HTML

/// The parts of one built HTML file the checks read.
#[derive(Default)]
struct Page {
    /// `(tag, url, in <head>)` of every href and src.
    urls: Vec<(String, String, bool)>,
    ids: HashSet<String>,
    copies: Vec<String>,
    nav: Vec<String>,
    styles: Vec<String>,
    canonical: Option<String>,
    /// The `<head>` title text, once one opens.
    title: Option<String>,
    in_title: bool,
    /// Carries the theme's md-content block: a page, not a copied asset.
    content: bool,
    /// Raw tags that hide the rest of the page.
    swallowed: Vec<String>,
    tokens: Vec<String>,
    in_head: bool,
    article: usize,
    template: usize,
    /// `true` inside `<svg>`/`<math>`, `false` again in `<foreignObject>`.
    foreign: Vec<bool>,
    /// Swallowing tags waiting for their end tag.
    open: Vec<String>,
    /// A redirect stub: its links are checked, it is not a page.
    stub: bool,
}

impl Page {
    fn parse(text: &str) -> Page {
        let mut page = Page {
            in_head: true,
            foreign: vec![false],
            ..Page::default()
        };
        tokenize(text, &mut page);
        let open = std::mem::take(&mut page.open);
        page.swallowed.extend(open);
        page.stub = strip(page.title.as_deref().unwrap_or(""), None) == STUB;
        page
    }

    fn in_foreign(&self) -> bool {
        self.foreign.last().copied().unwrap_or(false)
    }

    fn start(&mut self, tag: &str, attrs: Vec<(String, Option<String>)>) {
        self.tag(tag, attrs, false);
    }

    /// A self-closed tag: empty only for a void element or in SVG/MathML;
    /// HTML ignores the slash elsewhere, so `<textarea/>` opens a textarea.
    fn start_end(&mut self, tag: &str, attrs: Vec<(String, Option<String>)>) {
        let foreign = self.in_foreign() || tag == "svg" || tag == "math";
        let empty = VOID.contains(&tag) || foreign;
        self.tag(tag, attrs, empty);
    }

    fn tag(&mut self, tag: &str, attrs: Vec<(String, Option<String>)>, empty: bool) {
        // A later duplicate attribute wins, as in a dict built from the list.
        let mut a: HashMap<String, String> = HashMap::new();
        for (k, v) in attrs {
            a.insert(k, strip(v.as_deref().unwrap_or(""), None).to_owned());
        }
        if tag == "body" {
            self.in_head = false;
        }
        if tag == "title" && self.in_head && self.title.is_none() && !empty {
            self.title = Some(String::new());
            self.in_title = true;
        }
        if !empty {
            if tag == "article" {
                self.article += 1;
            } else if tag == "template" {
                self.template += 1;
            } else if tag == "svg" || tag == "math" {
                self.foreign.push(true);
            } else if tag == "foreignobject" && self.in_foreign() {
                self.foreign.push(false);
            }
        }
        if self.template == 0 {
            // ids inside <template> are not in the document.
            let keys: &[&str] = if tag == "a" { &["id", "name"] } else { &["id"] };
            for key in keys {
                if let Some(v) = a.get(*key).filter(|v| !v.is_empty()) {
                    self.ids.insert(v.clone());
                }
            }
        }
        for key in ["href", "src"] {
            if let Some(v) = a.get(key) {
                self.urls.push((tag.to_owned(), v.clone(), self.in_head));
            }
        }
        if let Some(v) = a.get("data-md-copy-url") {
            self.copies.push(v.clone());
        }
        let rel_attr = a.get("rel").map_or("", String::as_str);
        let rel = split(rel_attr, None, None);
        if tag == "link" && rel.contains(&"canonical") {
            self.canonical = a.get("href").cloned();
        }
        if tag == "link" && rel.contains(&"stylesheet") {
            self.styles.push(a.get("href").cloned().unwrap_or_default());
        }
        let class_attr = a.get("class").map_or("", String::as_str);
        let classes = split(class_attr, None, None);
        if classes.contains(&"md-content") {
            self.content = true;
        }
        // The toc's #heading entries carry md-nav__link too: they point at
        // the page itself, so only links with a path count toward
        // reachability.
        if tag == "a"
            && (classes.contains(&"md-nav__link") || classes.contains(&"md-tabs__link"))
            && let Some(href) = a.get("href")
            && !url_path(href).is_empty()
        {
            self.nav.push(href.clone());
        }
        // In an article and in HTML context (an SVG <title> is the drawing's
        // name), a raw <title> or <plaintext> always hides text; the others
        // only when never closed.
        if self.article > 0 && !self.in_foreign() && SWALLOW.contains(&tag) {
            if tag == "title" || tag == "plaintext" {
                self.swallowed.push(tag.to_owned());
            } else if !empty {
                self.open.push(tag.to_owned());
            }
        }
    }

    fn end(&mut self, tag: &str) {
        if tag == "title" {
            self.in_title = false;
        }
        if tag == "article" && self.article > 0 {
            self.article -= 1;
        } else if tag == "template" && self.template > 0 {
            self.template -= 1;
        } else if self.foreign.len() > 1 {
            let closes = match tag {
                "svg" | "math" => Some(true),
                "foreignobject" => Some(false),
                _ => None,
            };
            if closes == Some(self.in_foreign()) {
                self.foreign.pop();
            }
        }
        if let Some(i) = self.open.iter().position(|t| t == tag) {
            self.open.remove(i);
        }
    }

    fn data(&mut self, data: &str) {
        if self.in_title
            && let Some(title) = self.title.as_mut()
        {
            title.push_str(data);
        }
    }

    /// An awesome list's count token is a comment: report one left in an
    /// article.
    fn comment(&mut self, data: &str) {
        if self.article > 0 && TOKENS.contains(&strip(data, None)) {
            self.tokens.push(format!("<!--{data}-->"));
        }
    }
}

/// Raw text elements: only their own end tag ends them.
const RAWTEXT: [&str; 6] = ["script", "style", "xmp", "iframe", "noembed", "noframes"];
/// Escapable raw text elements: as RAWTEXT, with character references.
const RCDATA: [&str; 2] = ["textarea", "title"];

fn is_ws(b: u8) -> bool {
    matches!(b, b'\t' | b'\n' | b'\r' | 0x0c | b' ')
}

/// Python 3.14 `html.parser.HTMLParser` (convert_charrefs=True) over a
/// whole document, calling the page's handlers in the same order.
fn tokenize(text: &str, page: &mut Page) {
    let b = text.as_bytes();
    let n = b.len();
    let mut i = 0;
    // (element, escapable) while inside a raw text element.
    let mut cdata: Option<(String, bool)> = None;
    let unescape_if = |s: &str, escapable: bool| {
        if escapable { unescape(s) } else { s.to_owned() }
    };
    while i < n {
        let escapable = cdata.as_ref().is_none_or(|c| c.1);
        let j = match &cdata {
            None => find_byte(b, b'<', i).unwrap_or(n),
            Some((elem, _)) => match find_cdata_end(b, i, elem) {
                Some(j) => j,
                None => break,
            },
        };
        if i < j {
            page.data(&unescape_if(&text[i..j], escapable));
        }
        i = j;
        if i == n {
            break;
        }
        let rest = &b[i..];
        let k = if rest.len() > 1 && rest[1].is_ascii_alphabetic() {
            parse_starttag(text, i, page, &mut cdata)
        } else if rest.starts_with(b"</") {
            let k = parse_endtag(text, i, page);
            if k.is_some() && cdata.is_some() {
                cdata = None;
            }
            k
        } else if rest.starts_with(b"<!--") {
            parse_comment(text, i, page)
        } else if rest.starts_with(b"<?") {
            find_byte(b, b'>', i + 2).map(|p| p + 1)
        } else if rest.starts_with(b"<!") {
            parse_declaration(text, i, page)
        } else {
            page.data("<");
            Some(i + 1)
        };
        i = match k {
            Some(k) => k,
            None => {
                // The end of the document inside an unfinished construct.
                if rest.len() > 1 && rest[1].is_ascii_alphabetic() {
                } else if rest.starts_with(b"</") {
                    if i + 2 == n {
                        page.data("</");
                    } else if !rest[2].is_ascii_alphabetic() {
                        page.comment(&text[i + 2..]);
                    }
                } else if rest.starts_with(b"<!--") {
                    let body = &text[i + 4..];
                    let body = ["--!", "--", "-"]
                        .iter()
                        .find_map(|s| body.strip_suffix(s))
                        .unwrap_or(body);
                    page.comment(body);
                } else if rest.starts_with(b"<![CDATA[")
                    || (rest.len() >= 9 && rest[..9].eq_ignore_ascii_case(b"<!doctype"))
                {
                } else if rest.starts_with(b"<!") {
                    page.comment(&text[i + 2..]);
                }
                n
            }
        };
    }
    if i < n {
        let escapable = cdata.as_ref().is_none_or(|c| c.1);
        page.data(&unescape_if(&text[i..n], escapable));
    }
}

fn find_byte(b: &[u8], needle: u8, from: usize) -> Option<usize> {
    b.get(from..)?
        .iter()
        .position(|&c| c == needle)
        .map(|p| p + from)
}

fn find_bytes(b: &[u8], needle: &[u8], from: usize) -> Option<usize> {
    b.get(from..)?
        .windows(needle.len())
        .position(|w| w == needle)
        .map(|p| p + from)
}

/// The next `</elem` followed by whitespace, `/` or `>` (ASCII case
/// insensitive); never for `plaintext`, which runs to the end.
fn find_cdata_end(b: &[u8], from: usize, elem: &str) -> Option<usize> {
    if elem == "plaintext" {
        return None;
    }
    let e = elem.as_bytes();
    let mut i = from;
    while let Some(p) = find_bytes(b, b"</", i) {
        let name = p + 2;
        if b.len() > name + e.len()
            && b[name..name + e.len()].eq_ignore_ascii_case(e)
            && (is_ws(b[name + e.len()]) || matches!(b[name + e.len()], b'/' | b'>'))
        {
            return Some(p);
        }
        i = p + 1;
    }
    None
}

/// `locatetagend` from `i` (the tag name's first letter): the end of the
/// tag, past its `>` when it has one.
fn locate_tag_end(b: &[u8], mut i: usize) -> usize {
    let n = b.len();
    i += 1;
    while i < n && !(is_ws(b[i]) || b[i] == b'/' || b[i] == b'>') {
        i += 1;
    }
    while i < n && (is_ws(b[i]) || b[i] == b'/') {
        i += 1;
    }
    while i < n && attr_may_start(b, i) {
        i = attr_name_end(b, i);
        i = attr_value_end(b, i).unwrap_or(i);
        while i < n && (is_ws(b[i]) || b[i] == b'/') {
            i += 1;
        }
    }
    if i < n && b[i] == b'>' {
        i += 1;
    }
    i
}

/// An attribute name may start at `i`: the byte before is a quote,
/// whitespace or `/`, and the byte itself is not whitespace, `/` or `>`.
fn attr_may_start(b: &[u8], i: usize) -> bool {
    i > 0
        && i < b.len()
        && (is_ws(b[i - 1]) || matches!(b[i - 1], b'\'' | b'"' | b'/'))
        && !(is_ws(b[i]) || b[i] == b'/' || b[i] == b'>')
}

fn attr_name_end(b: &[u8], mut i: usize) -> usize {
    i += 1;
    while i < b.len() && !(is_ws(b[i]) || matches!(b[i], b'/' | b'=' | b'>')) {
        i += 1;
    }
    i
}

/// The end of `[ws]*=[ws]*value` after an attribute name, as the regex
/// with its backtracking matches it; `None` when there is no value.
fn attr_value_end(b: &[u8], name_end: usize) -> Option<usize> {
    let n = b.len();
    let mut j = name_end;
    while j < n && is_ws(b[j]) {
        j += 1;
    }
    if j >= n || b[j] != b'=' {
        return None;
    }
    let after_eq = j + 1;
    let mut m = after_eq;
    while m < n && is_ws(b[m]) {
        m += 1;
    }
    match b.get(m) {
        Some(&q @ (b'\'' | b'"')) => match find_byte(b, q, m + 1) {
            Some(close) => Some(close + 1),
            // An unclosed quote: the regex backs off to an empty bare
            // value right after `=`, which needs a whitespace there.
            None if m > after_eq => Some(after_eq),
            None => None,
        },
        _ => {
            let mut k = m;
            while k < n && !(is_ws(b[k]) || b[k] == b'>') {
                k += 1;
            }
            Some(k)
        }
    }
}

/// `(?:[\t\n\r\f ]|/(?!>))*` from `i`.
fn skip_ws_slash(b: &[u8], mut i: usize) -> usize {
    while i < b.len() && (is_ws(b[i]) || (b[i] == b'/' && b.get(i + 1) != Some(&b'>'))) {
        i += 1;
    }
    i
}

fn parse_starttag(
    text: &str,
    i: usize,
    page: &mut Page,
    cdata: &mut Option<(String, bool)>,
) -> Option<usize> {
    let b = text.as_bytes();
    let endpos = locate_tag_end(b, i + 1);
    if b[endpos - 1] != b'>' {
        return None;
    }
    let mut name_end = i + 2;
    while name_end < b.len() && !(is_ws(b[name_end]) || b[name_end] == b'/' || b[name_end] == b'>')
    {
        name_end += 1;
    }
    let tag = lower(&text[i + 1..name_end]);
    let mut k = skip_ws_slash(b, name_end);
    let mut attrs = Vec::new();
    while k < endpos && attr_may_start(b, k) {
        let name_stop = attr_name_end(b, k);
        let name = lower(&text[k..name_stop]);
        let value = match attr_value_end(b, name_stop) {
            None => {
                k = name_stop;
                None
            }
            Some(end) => {
                let eq = name_stop + text[name_stop..end].find('=').unwrap_or(0);
                let raw = strip_ws(&text[eq + 1..end]);
                let raw = match raw.as_bytes().first() {
                    Some(&q @ (b'\'' | b'"'))
                        if raw.len() >= 2 && raw.as_bytes()[raw.len() - 1] == q =>
                    {
                        &raw[1..raw.len() - 1]
                    }
                    _ => raw,
                };
                k = end;
                Some(if raw.is_empty() {
                    String::new()
                } else {
                    unescape_attr(raw)
                })
            }
        };
        attrs.push((name, value));
        k = skip_ws_slash(b, k);
    }
    let end = strip(&text[k..endpos], None);
    if end != ">" && end != "/>" {
        page.data(&text[i..endpos]);
        return Some(endpos);
    }
    if end == "/>" {
        page.start_end(&tag, attrs);
    } else {
        page.start(&tag, attrs);
        if RAWTEXT.contains(&tag.as_str()) || tag == "plaintext" {
            *cdata = Some((tag, false));
        } else if RCDATA.contains(&tag.as_str()) {
            *cdata = Some((tag, true));
        }
    }
    Some(endpos)
}

/// The value part of `[ws]*=[ws]*value` without the whitespace around `=`
/// (the leading whitespace before the value).
fn strip_ws(s: &str) -> &str {
    s.trim_start_matches(['\t', '\n', '\r', '\x0c', ' '])
}

fn parse_endtag(text: &str, i: usize, page: &mut Page) -> Option<usize> {
    let b = text.as_bytes();
    find_byte(b, b'>', i + 2)?;
    if !b.get(i + 2).is_some_and(u8::is_ascii_alphabetic) {
        if b.get(i + 2) == Some(&b'>') {
            return Some(i + 3); // </> is ignored
        }
        let pos = find_byte(b, b'>', i + 2)?;
        page.comment(&text[i + 2..pos]);
        return Some(pos + 1);
    }
    let j = locate_tag_end(b, i + 2);
    if b[j - 1] != b'>' {
        return None;
    }
    let mut name_end = i + 3;
    while name_end < b.len() && !(is_ws(b[name_end]) || b[name_end] == b'/' || b[name_end] == b'>')
    {
        name_end += 1;
    }
    page.end(&lower(&text[i + 2..name_end]));
    Some(j)
}

fn parse_comment(text: &str, i: usize, page: &mut Page) -> Option<usize> {
    let b = text.as_bytes();
    let body = i + 4;
    // An empty comment is abruptly closed by the first ">" or "->".
    if b.get(body) == Some(&b'>') {
        page.comment("");
        return Some(body + 1);
    }
    if b.get(body) == Some(&b'-') && b.get(body + 1) == Some(&b'>') {
        page.comment("");
        return Some(body + 2);
    }
    let mut p = body;
    while let Some(q) = find_bytes(b, b"--", p) {
        if b.get(q + 2) == Some(&b'>') {
            page.comment(&text[body..q]);
            return Some(q + 3);
        }
        if b.get(q + 2) == Some(&b'!') && b.get(q + 3) == Some(&b'>') {
            page.comment(&text[body..q]);
            return Some(q + 4);
        }
        p = q + 1;
    }
    None
}

fn parse_declaration(text: &str, i: usize, page: &mut Page) -> Option<usize> {
    let b = text.as_bytes();
    if b[i..].starts_with(b"<![CDATA[") {
        return find_bytes(b, b"]]>", i + 9).map(|j| j + 3);
    }
    if b.len() >= i + 9 && b[i..i + 9].eq_ignore_ascii_case(b"<!doctype") {
        return find_byte(b, b'>', i + 9).map(|p| p + 1);
    }
    let pos = find_byte(b, b'>', i + 2)?;
    page.comment(&text[i + 2..pos]);
    Some(pos + 1)
}

/// The named character references decoded: the ones a Markdown build
/// writes. `legacy` ones are also decoded without their `;`.
const ENTITIES: [(&str, &str, bool); 6] = [
    ("amp", "&", true),
    ("lt", "<", true),
    ("gt", ">", true),
    ("quot", "\"", true),
    ("apos", "'", false),
    ("nbsp", "\u{a0}", true),
];

fn entity(name: &str) -> Option<&'static str> {
    let (bare, semi) = match name.strip_suffix(';') {
        Some(b) => (b, true),
        None => (name, false),
    };
    ENTITIES
        .iter()
        .find(|(n, _, legacy)| *n == bare && (semi || *legacy))
        .map(|(_, v, _)| *v)
}

/// A numeric reference's text (`#38`, `#x26`, with or without `;`) as
/// `html.unescape` decodes it.
fn numeric_ref(s: &str) -> String {
    let body = s.trim_end_matches(';');
    let (digits, radix) = match body.strip_prefix(['x', 'X']) {
        Some(h) => (h, 16),
        None => (body, 10),
    };
    let num = u64::from_str_radix(digits, radix).unwrap_or(u64::MAX);
    const CP1252: [u32; 32] = [
        0x20ac, 0x81, 0x201a, 0x0192, 0x201e, 0x2026, 0x2020, 0x2021, 0x02c6, 0x2030, 0x0160,
        0x2039, 0x0152, 0x8d, 0x017d, 0x8f, 0x90, 0x2018, 0x2019, 0x201c, 0x201d, 0x2022, 0x2013,
        0x2014, 0x02dc, 0x2122, 0x0161, 0x203a, 0x0153, 0x9d, 0x017e, 0x0178,
    ];
    let pick = |c: u32| char::from_u32(c).map(String::from).unwrap_or_default();
    match num {
        0 => "\u{fffd}".to_owned(),
        0x0d => "\r".to_owned(),
        0x80..=0x9f => pick(CP1252[usize::try_from(num - 0x80).unwrap_or(0)]),
        0xd800..=0xdfff => "\u{fffd}".to_owned(),
        n if n > 0x10ffff => "\u{fffd}".to_owned(),
        0x1..=0x8 | 0xb | 0xe..=0x1f | 0x7f | 0xfdd0..=0xfdef => String::new(),
        n if n & 0xfffe == 0xfffe => String::new(),
        n => pick(u32::try_from(n).unwrap_or(0xfffd)),
    }
}

/// `html.unescape(s)` for data: `&` followed by a reference.
fn unescape(s: &str) -> String {
    if !s.contains('&') {
        return s.to_owned();
    }
    let mut out = String::with_capacity(s.len());
    let mut rest = s;
    while let Some(amp) = rest.find('&') {
        out.push_str(&rest[..amp]);
        let after = &rest[amp + 1..];
        let ab = after.as_bytes();
        if ab.first() == Some(&b'#') {
            let hex = matches!(ab.get(1), Some(b'x' | b'X'));
            let start = if hex { 2 } else { 1 };
            let mut end = start;
            while end < ab.len()
                && (if hex {
                    ab[end].is_ascii_hexdigit()
                } else {
                    ab[end].is_ascii_digit()
                })
            {
                end += 1;
            }
            if end > start {
                let end = if ab.get(end) == Some(&b';') {
                    end + 1
                } else {
                    end
                };
                out.push_str(&numeric_ref(&after[1..end]));
                rest = &after[end..];
                continue;
            }
            out.push('&');
            rest = after;
            continue;
        }
        // [^\t\n\f <&#;]{1,32};?
        let mut end = 0;
        for (count, (pos, c)) in after.char_indices().enumerate() {
            if count == 32 || matches!(c, '\t' | '\n' | '\x0c' | ' ' | '<' | '&' | '#' | ';') {
                break;
            }
            end = pos + c.len_utf8();
        }
        if end == 0 {
            out.push('&');
            rest = after;
            continue;
        }
        let with_semi = if after[end..].starts_with(';') {
            end + 1
        } else {
            end
        };
        let name = &after[..with_semi];
        if let Some(v) = entity(name) {
            out.push_str(v);
        } else if let Some((x, v)) = (2..name.len())
            .rev()
            .filter(|&x| name.is_char_boundary(x))
            .find_map(|x| entity(&name[..x]).map(|v| (x, v)))
        {
            out.push_str(v);
            out.push_str(&name[x..]);
        } else {
            out.push('&');
            out.push_str(name);
        }
        rest = &after[with_semi..];
    }
    out.push_str(rest);
    out
}

/// Character references in an attribute value: numeric ones always; a
/// named one only when it is exact and not followed by `=`.
fn unescape_attr(s: &str) -> String {
    if !s.contains('&') {
        return s.to_owned();
    }
    let mut out = String::with_capacity(s.len());
    let mut rest = s;
    while let Some(amp) = rest.find('&') {
        out.push_str(&rest[..amp]);
        let after = &rest[amp + 1..];
        let ab = after.as_bytes();
        let body_end = if ab.first() == Some(&b'#') {
            let hex = matches!(ab.get(1), Some(b'x' | b'X'));
            let start = if hex { 2 } else { 1 };
            let mut end = start;
            while end < ab.len()
                && (if hex {
                    ab[end].is_ascii_hexdigit()
                } else {
                    ab[end].is_ascii_digit()
                })
            {
                end += 1;
            }
            (end > start).then_some(end)
        } else if ab.first().is_some_and(u8::is_ascii_alphabetic) {
            let mut end = 1;
            while end < ab.len() && ab[end].is_ascii_alphanumeric() {
                end += 1;
            }
            Some(end)
        } else {
            None
        };
        let Some(end) = body_end else {
            out.push('&');
            rest = after;
            continue;
        };
        let term = matches!(ab.get(end), Some(b';' | b'='));
        let full = if term { end + 1 } else { end };
        let reference = &after[..full];
        if reference.starts_with('#') || (!reference.ends_with('=') && entity(reference).is_some())
        {
            out.push_str(&unescape(&format!("&{reference}")));
        } else {
            out.push('&');
            out.push_str(reference);
        }
        rest = &after[full..];
    }
    out.push_str(rest);
    out
}

// ------------------------------------------------------------------- TOML

/// A TOML 1.0 reader for the config: tables, arrays of tables, dotted
/// keys, inline tables, the four string forms, integers, floats and
/// booleans. Dates and times are kept as their text.
mod toml {
    use indexmap::IndexMap;
    use pycompat::json::Value;
    use std::collections::HashSet;

    type Table = IndexMap<String, Value>;

    pub fn parse(src: &str) -> Result<Table, String> {
        let mut p = Parser {
            s: src.as_bytes(),
            src,
            i: 0,
        };
        let mut root = Table::new();
        let mut current: Vec<String> = Vec::new();
        let mut defined: HashSet<Vec<String>> = HashSet::new();
        // Keys given a value by `key = value`: an inline table or array is
        // complete, so no header or dotted key may extend it.
        let mut fixed: HashSet<Vec<String>> = HashSet::new();
        loop {
            p.skip_blank_lines();
            if p.i >= p.s.len() {
                return Ok(root);
            }
            if p.peek() == Some(b'[') {
                let array = p.s.get(p.i + 1) == Some(&b'[');
                p.i += if array { 2 } else { 1 };
                p.skip_ws();
                let path = p.key()?;
                p.skip_ws();
                let close: &[u8] = if array { b"]]" } else { b"]" };
                if !p.s[p.i..].starts_with(close) {
                    return Err(p.err("expected the end of a table header"));
                }
                p.i += close.len();
                p.end_of_line()?;
                if (1..=path.len()).any(|n| fixed.contains(&path[..n])) {
                    return Err(p.err(&format!("cannot extend {}, defined inline", path.join("."))));
                }
                if array {
                    let (last, parent) = path.split_last().ok_or_else(|| p.err("empty header"))?;
                    let table = table_at(&mut root, parent).map_err(|e| p.err(&e))?;
                    match table
                        .entry(last.clone())
                        .or_insert_with(|| Value::Array(Vec::new()))
                    {
                        Value::Array(items) => items.push(Value::Object(Table::new())),
                        _ => return Err(p.err(&format!("{last} is not an array of tables"))),
                    }
                    // A new table of the array: its keys and sub-tables start over.
                    fixed.retain(|k| !k.starts_with(&path));
                    defined.retain(|k| !(k.len() > path.len() && k.starts_with(&path)));
                } else {
                    if !defined.insert(path.clone()) {
                        return Err(p.err(&format!("table [{}] defined twice", path.join("."))));
                    }
                    table_at(&mut root, &path).map_err(|e| p.err(&e))?;
                }
                current = path;
                continue;
            }
            let keys = p.key()?;
            p.skip_ws();
            if p.peek() != Some(b'=') {
                return Err(p.err("expected '=' after a key"));
            }
            p.i += 1;
            p.skip_ws();
            let value = p.value()?;
            p.end_of_line()?;
            let full: Vec<String> = current.iter().chain(&keys).cloned().collect();
            if (1..full.len()).any(|n| fixed.contains(&full[..n])) {
                return Err(p.err(&format!("cannot extend {}, defined inline", full.join("."))));
            }
            let table = table_at(&mut root, &current).map_err(|e| p.err(&e))?;
            insert_dotted(table, &keys, value).map_err(|e| p.err(&e))?;
            fixed.insert(full);
        }
    }

    /// The table at `path`, made where missing; through an array of
    /// tables, its last table.
    fn table_at<'r>(root: &'r mut Table, path: &[String]) -> Result<&'r mut Table, String> {
        let mut t = root;
        for key in path {
            let v = t
                .entry(key.clone())
                .or_insert_with(|| Value::Object(Table::new()));
            t = match v {
                Value::Object(o) => o,
                Value::Array(items) => match items.last_mut() {
                    Some(Value::Object(o)) => o,
                    _ => return Err(format!("{key} is not a table")),
                },
                _ => return Err(format!("{key} is not a table")),
            };
        }
        Ok(t)
    }

    fn insert_dotted(table: &mut Table, keys: &[String], value: Value) -> Result<(), String> {
        let (last, parents) = keys.split_last().ok_or("empty key")?;
        let mut t = table;
        for key in parents {
            let v = t
                .entry(key.clone())
                .or_insert_with(|| Value::Object(Table::new()));
            t = match v {
                Value::Object(o) => o,
                _ => return Err(format!("{key} is not a table")),
            };
        }
        if t.contains_key(last) {
            return Err(format!("key {last} defined twice"));
        }
        t.insert(last.clone(), value);
        Ok(())
    }

    struct Parser<'a> {
        s: &'a [u8],
        src: &'a str,
        i: usize,
    }

    impl Parser<'_> {
        fn err(&self, msg: &str) -> String {
            let line = self.s[..self.i.min(self.s.len())]
                .iter()
                .filter(|&&b| b == b'\n')
                .count()
                + 1;
            format!("{msg} (at line {line})")
        }

        fn peek(&self) -> Option<u8> {
            self.s.get(self.i).copied()
        }

        fn skip_ws(&mut self) {
            while matches!(self.peek(), Some(b' ' | b'\t')) {
                self.i += 1;
            }
        }

        fn skip_comment(&mut self) {
            if self.peek() == Some(b'#') {
                while !matches!(self.peek(), None | Some(b'\n')) {
                    self.i += 1;
                }
            }
        }

        /// Whitespace, comments and newlines.
        fn skip_blank_lines(&mut self) {
            loop {
                self.skip_ws();
                self.skip_comment();
                match self.peek() {
                    Some(b'\n') => self.i += 1,
                    Some(b'\r') if self.s.get(self.i + 1) == Some(&b'\n') => self.i += 2,
                    _ => return,
                }
            }
        }

        fn end_of_line(&mut self) -> Result<(), String> {
            self.skip_ws();
            self.skip_comment();
            match self.peek() {
                None => Ok(()),
                Some(b'\n') => {
                    self.i += 1;
                    Ok(())
                }
                Some(b'\r') if self.s.get(self.i + 1) == Some(&b'\n') => {
                    self.i += 2;
                    Ok(())
                }
                _ => Err(self.err("expected the end of the line")),
            }
        }

        /// A dotted key: bare, basic or literal parts joined by `.`.
        fn key(&mut self) -> Result<Vec<String>, String> {
            let mut parts = Vec::new();
            loop {
                self.skip_ws();
                let part = match self.peek() {
                    Some(b'"') => {
                        if self.s[self.i..].starts_with(b"\"\"\"") {
                            return Err(self.err("a key cannot be a multi-line string"));
                        }
                        self.basic_string()?
                    }
                    Some(b'\'') => {
                        if self.s[self.i..].starts_with(b"'''") {
                            return Err(self.err("a key cannot be a multi-line string"));
                        }
                        self.literal_string()?
                    }
                    _ => {
                        let start = self.i;
                        while self
                            .peek()
                            .is_some_and(|b| b.is_ascii_alphanumeric() || b == b'_' || b == b'-')
                        {
                            self.i += 1;
                        }
                        if start == self.i {
                            return Err(self.err("expected a key"));
                        }
                        self.src[start..self.i].to_owned()
                    }
                };
                parts.push(part);
                self.skip_ws();
                if self.peek() == Some(b'.') {
                    self.i += 1;
                } else {
                    return Ok(parts);
                }
            }
        }

        fn value(&mut self) -> Result<Value, String> {
            let rest = &self.s[self.i..];
            match self.peek() {
                Some(b'"') if rest.starts_with(b"\"\"\"") => self.ml_basic_string().map(Value::Str),
                Some(b'"') => self.basic_string().map(Value::Str),
                Some(b'\'') if rest.starts_with(b"'''") => self.ml_literal_string().map(Value::Str),
                Some(b'\'') => self.literal_string().map(Value::Str),
                Some(b'[') => self.array(),
                Some(b'{') => self.inline_table(),
                Some(_) if rest.starts_with(b"true") && !word_continues(rest, 4) => {
                    self.i += 4;
                    Ok(Value::Bool(true))
                }
                Some(_) if rest.starts_with(b"false") && !word_continues(rest, 5) => {
                    self.i += 5;
                    Ok(Value::Bool(false))
                }
                Some(_) => self.scalar(),
                None => Err(self.err("expected a value")),
            }
        }

        fn array(&mut self) -> Result<Value, String> {
            self.i += 1;
            let mut items = Vec::new();
            loop {
                self.skip_blank_lines();
                if self.peek() == Some(b']') {
                    self.i += 1;
                    return Ok(Value::Array(items));
                }
                items.push(self.value()?);
                self.skip_blank_lines();
                match self.peek() {
                    Some(b',') => self.i += 1,
                    Some(b']') => {
                        self.i += 1;
                        return Ok(Value::Array(items));
                    }
                    _ => return Err(self.err("expected ',' or ']' in an array")),
                }
            }
        }

        fn inline_table(&mut self) -> Result<Value, String> {
            self.i += 1;
            let mut table = Table::new();
            self.skip_ws();
            if self.peek() == Some(b'}') {
                self.i += 1;
                return Ok(Value::Object(table));
            }
            loop {
                self.skip_blank_lines();
                let keys = self.key()?;
                self.skip_ws();
                if self.peek() != Some(b'=') {
                    return Err(self.err("expected '=' in an inline table"));
                }
                self.i += 1;
                self.skip_ws();
                let v = self.value()?;
                insert_dotted(&mut table, &keys, v).map_err(|e| self.err(&e))?;
                self.skip_blank_lines();
                match self.peek() {
                    Some(b',') => self.i += 1,
                    Some(b'}') => {
                        self.i += 1;
                        return Ok(Value::Object(table));
                    }
                    _ => return Err(self.err("expected ',' or '}' in an inline table")),
                }
            }
        }

        fn literal_string(&mut self) -> Result<String, String> {
            self.i += 1;
            let start = self.i;
            while let Some(b) = self.peek() {
                if b == b'\'' {
                    let s = self.src[start..self.i].to_owned();
                    self.i += 1;
                    return Ok(s);
                }
                if b == b'\n' {
                    break;
                }
                self.i += 1;
            }
            Err(self.err("unterminated string"))
        }

        fn ml_literal_string(&mut self) -> Result<String, String> {
            self.i += 3;
            self.skip_first_newline();
            let start = self.i;
            let close = super::find_bytes(self.s, b"'''", start)
                .ok_or_else(|| self.err("unterminated string"))?;
            let mut end = close;
            // Up to two quotes right before the delimiter belong to the string.
            while end + 3 < self.s.len() && self.s[end + 3] == b'\'' && end + 3 - close < 2 {
                end += 1;
            }
            let s = self.src[start..end].to_owned();
            self.i = end + 3;
            Ok(s)
        }

        fn skip_first_newline(&mut self) {
            if self.peek() == Some(b'\n') {
                self.i += 1;
            } else if self.s[self.i..].starts_with(b"\r\n") {
                self.i += 2;
            }
        }

        fn basic_string(&mut self) -> Result<String, String> {
            self.i += 1;
            let mut out = String::new();
            loop {
                match self.peek() {
                    None | Some(b'\n') => return Err(self.err("unterminated string")),
                    Some(b'"') => {
                        self.i += 1;
                        return Ok(out);
                    }
                    Some(b'\\') => self.escape(&mut out, false)?,
                    Some(_) => self.push_char(&mut out),
                }
            }
        }

        fn ml_basic_string(&mut self) -> Result<String, String> {
            self.i += 3;
            self.skip_first_newline();
            let mut out = String::new();
            loop {
                if self.s[self.i..].starts_with(b"\"\"\"") {
                    let mut extra = 0;
                    while extra < 2 && self.s.get(self.i + 3 + extra) == Some(&b'"') {
                        extra += 1;
                    }
                    out.push_str(&"\"".repeat(extra));
                    self.i += 3 + extra;
                    return Ok(out);
                }
                match self.peek() {
                    None => return Err(self.err("unterminated string")),
                    Some(b'\\') => self.escape(&mut out, true)?,
                    Some(_) => self.push_char(&mut out),
                }
            }
        }

        fn push_char(&mut self, out: &mut String) {
            let c = self.src[self.i..].chars().next().unwrap_or('\u{fffd}');
            out.push(c);
            self.i += c.len_utf8();
        }

        fn escape(&mut self, out: &mut String, multiline: bool) -> Result<(), String> {
            self.i += 1;
            let c = self.peek().ok_or_else(|| self.err("unterminated string"))?;
            self.i += 1;
            match c {
                b'b' => out.push('\u{8}'),
                b't' => out.push('\t'),
                b'n' => out.push('\n'),
                b'f' => out.push('\u{c}'),
                b'r' => out.push('\r'),
                b'e' => out.push('\u{1b}'),
                b'"' => out.push('"'),
                b'\\' => out.push('\\'),
                b'u' | b'U' => {
                    let len = if c == b'u' { 4 } else { 8 };
                    let hex = self
                        .src
                        .get(self.i..self.i + len)
                        .ok_or_else(|| self.err("short unicode escape"))?;
                    let code = u32::from_str_radix(hex, 16)
                        .ok()
                        .and_then(char::from_u32)
                        .ok_or_else(|| self.err("invalid unicode escape"))?;
                    out.push(code);
                    self.i += len;
                }
                b' ' | b'\t' | b'\r' | b'\n' if multiline => {
                    // A line-ending backslash trims the whitespace and
                    // newlines that follow.
                    self.i -= 1;
                    let mut j = self.i;
                    while matches!(self.s.get(j), Some(b' ' | b'\t')) {
                        j += 1;
                    }
                    if !matches!(self.s.get(j), Some(b'\n' | b'\r')) {
                        return Err(self.err("invalid escape"));
                    }
                    while matches!(self.s.get(j), Some(b' ' | b'\t' | b'\n' | b'\r')) {
                        j += 1;
                    }
                    self.i = j;
                }
                _ => return Err(self.err("invalid escape")),
            }
            Ok(())
        }

        /// A number, a date or a time.
        fn scalar(&mut self) -> Result<Value, String> {
            let start = self.i;
            while let Some(b) = self.peek() {
                let date_space = b == b' '
                    && self.i - start == 10
                    && self.s.get(self.i + 1).is_some_and(u8::is_ascii_digit)
                    && self.s[start..self.i].iter().filter(|&&c| c == b'-').count() == 2;
                if b.is_ascii_alphanumeric()
                    || matches!(b, b'_' | b'+' | b'-' | b'.' | b':')
                    || date_space
                {
                    self.i += 1;
                } else {
                    break;
                }
            }
            let tok = &self.src[start..self.i];
            if tok.is_empty() {
                return Err(self.err("expected a value"));
            }
            let tb = tok.as_bytes();
            let is_date = tb.len() >= 10 && tb[..4].iter().all(u8::is_ascii_digit) && tb[4] == b'-';
            let is_time = tb.len() >= 8 && tb[..2].iter().all(u8::is_ascii_digit) && tb[2] == b':';
            if is_date || is_time {
                return Ok(Value::Str(tok.to_owned()));
            }
            number(tok).ok_or_else(|| self.err(&format!("invalid value {tok}")))
        }
    }

    fn word_continues(rest: &[u8], at: usize) -> bool {
        rest.get(at)
            .is_some_and(|b| b.is_ascii_alphanumeric() || matches!(b, b'_' | b'-'))
    }

    fn number(tok: &str) -> Option<Value> {
        let (sign, body) = match tok.as_bytes().first() {
            Some(b'+') => (1.0, &tok[1..]),
            Some(b'-') => (-1.0, &tok[1..]),
            _ => (1.0, tok),
        };
        match body {
            "inf" => return Some(Value::Float(sign * f64::INFINITY)),
            "nan" => return Some(Value::Float(f64::NAN)),
            _ => {}
        }
        let underscores_ok =
            |s: &str| !s.starts_with('_') && !s.ends_with('_') && !s.contains("__");
        for (prefix, radix) in [("0x", 16), ("0o", 8), ("0b", 2)] {
            if let Some(digits) = tok.strip_prefix(prefix) {
                if !underscores_ok(digits) {
                    return None;
                }
                return i128::from_str_radix(&digits.replace('_', ""), radix)
                    .ok()
                    .filter(|_| !digits.is_empty())
                    .map(Value::Int);
            }
        }
        if !body.bytes().next().is_some_and(|b| b.is_ascii_digit()) {
            return None;
        }
        let parts_ok = body.split(['.', 'e', 'E']).all(|p| {
            let p = p.trim_start_matches(['+', '-']);
            !p.is_empty() && underscores_ok(p) && p.bytes().all(|b| b.is_ascii_digit() || b == b'_')
        });
        if !parts_ok {
            return None;
        }
        let clean = body.replace('_', "");
        let int_part = clean.split(['.', 'e', 'E']).next().unwrap_or("");
        if int_part.len() > 1 && int_part.starts_with('0') {
            return None;
        }
        if clean.contains(['.', 'e', 'E']) {
            clean.parse::<f64>().ok().map(|f| Value::Float(sign * f))
        } else {
            clean
                .parse::<i128>()
                .ok()
                .map(|n| Value::Int(if sign < 0.0 { -n } else { n }))
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn reads_every_form_the_config_uses() {
            let doc = parse(
                "# c\n[project]\nsite_url = \"https://x.io/R/\" # tail\nstrict = true\n\
                 exclude = \"\"\"\nA*.md\nb/\n\"\"\"\nnav = [\n  { Home = \"index.md\" },\n  { \"Get\" = [ { a = 'a.md' }, ] },\n]\n\
                 [project.theme]\nfeatures = [\"x.y\", 'z']\nicon.repo = \"gh\"\n\
                 [[project.theme.palette]]\nscheme = \"slate\"\ntoggle.icon = \"sun\"\n\
                 [[project.theme.palette]]\nscheme = \"default\"\n\
                 [project.plugins.search]\n[project.plugins.llmstxt.sections]\nDocs = [\"*.md\"]\n\
                 [project.markdown_extensions]\n\"pymdownx.snippets\".base_path = [\".\", \"docs\"]\nn = 1_000\nf = -0.5e2\nd = 1979-05-27 07:32:00Z\n",
            )
            .unwrap();
            let Value::Object(p) = &doc["project"] else {
                panic!()
            };
            assert_eq!(p["site_url"], Value::Str("https://x.io/R/".into()));
            assert_eq!(p["strict"], Value::Bool(true));
            assert_eq!(p["exclude"], Value::Str("A*.md\nb/\n".into()));
            let Value::Object(theme) = &p["theme"] else {
                panic!()
            };
            let Value::Array(palette) = &theme["palette"] else {
                panic!()
            };
            assert_eq!(palette.len(), 2);
            let Value::Object(plugins) = &p["plugins"] else {
                panic!()
            };
            assert_eq!(plugins.keys().collect::<Vec<_>>(), ["search", "llmstxt"]);
            let Value::Object(ext) = &p["markdown_extensions"] else {
                panic!()
            };
            assert_eq!(ext["n"], Value::Int(1000));
            assert_eq!(ext["f"], Value::Float(-50.0));
            assert_eq!(ext["d"], Value::Str("1979-05-27 07:32:00Z".into()));
        }

        #[test]
        fn refuses_broken_toml() {
            for bad in [
                "[project\nx = 1\n",
                "x = \n",
                "x = 1\nx = 2\n",
                "[a]\n[a]\n",
                "x = \"open\n",
                "x = [1, 2\n",
                "x = 1 y\n",
                "x = 01\n",
                "x = 1__0\n",
                // An inline value is complete: no header or dotted key extends it.
                "x = [{a = 1}]\n[x.b]\n",
                "x = [1]\n[[x]]\n",
                "x = {a = 1}\n[x.b]\n",
                "x = {a = 1}\nx.b = 2\n",
            ] {
                assert!(parse(bad).is_err(), "{bad:?} parsed");
            }
        }

        #[test]
        fn each_table_of_an_array_starts_over() {
            let doc =
                parse("[[a]]\nx = {p = 1}\n[a.b]\n[[a]]\nx = {p = 2}\n[a.b]\ny = 1\n").unwrap();
            let Value::Array(items) = &doc["a"] else {
                panic!()
            };
            assert_eq!(items.len(), 2);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn attrs_of(html: &str) -> Page {
        Page::parse(html)
    }

    #[test]
    fn attribute_values_in_every_quoting_with_entities() {
        let p = attrs_of(
            "<body><a href=a.html>x</a><a href='b.html'>x</a><a href=\"c&amp;d.html\">x</a>\
             <a HREF = \"e&#47;f.html\" >x</a><a href=\"g&ampx=1\">x</a><img src=i.png/>",
        );
        let urls: Vec<&str> = p.urls.iter().map(|u| u.1.as_str()).collect();
        assert_eq!(
            urls,
            [
                "a.html", "b.html", "c&d.html", "e/f.html", "g&ampx=1", "i.png/"
            ]
        );
    }

    #[test]
    fn script_style_and_comments_hold_no_links() {
        let p = attrs_of(
            "<body><script>var s = '<a href=\"x.html\">';</script><style>a{}</style>\
             <!-- <a href=\"y.html\"> --><a href=\"z.html\">z</a>",
        );
        let urls: Vec<&str> = p.urls.iter().map(|u| u.1.as_str()).collect();
        assert_eq!(urls, ["z.html"]);
    }

    #[test]
    fn swallowing_tags_only_in_html_context_of_an_article() {
        let svg = attrs_of("<body><article><svg><title>Drawing</title></svg><p>x</p></article>");
        assert!(svg.swallowed.is_empty());
        let raw = attrs_of("<body><article><p>a <title> b</p></article>");
        assert_eq!(raw.swallowed, ["title"]);
        let closed = attrs_of("<body><article><textarea>t</textarea></article>");
        assert!(closed.swallowed.is_empty());
        let self_closed = attrs_of("<body><article><textarea/></article>");
        assert_eq!(self_closed.swallowed, ["textarea"]);
        let outside = attrs_of("<body><textarea>never closed");
        assert!(outside.swallowed.is_empty());
    }

    #[test]
    fn head_title_marks_a_redirect_stub() {
        assert!(attrs_of("<head><title>Redirecting...</title></head>").stub);
        assert!(!attrs_of("<head><title>Home</title></head>").stub);
    }

    #[test]
    fn template_ids_and_count_tokens() {
        let p = attrs_of(
            "<body><template><p id=\"hidden\"></p></template><p id=\"shown\"></p><a name=\"anchor\"></a>\
             <article><!-- n:proyectos --><!-- n:other --></article>",
        );
        assert!(p.ids.contains("shown") && p.ids.contains("anchor") && !p.ids.contains("hidden"));
        assert_eq!(p.tokens, ["<!-- n:proyectos -->"]);
    }

    #[test]
    fn urlparse_as_python() {
        let u = urlparse("https://Geiserx.github.io:443/OwnGenome/a;p?q#f").unwrap();
        assert_eq!(u.scheme, "https");
        assert_eq!(u.hostname().as_deref(), Some("geiserx.github.io"));
        assert_eq!(u.port(), Ok(Some(443)));
        assert_eq!(u.path, "/OwnGenome/a");
        assert_eq!(u.fragment, "f");
        assert!(urlparse("http://[::1/").is_err());
        assert!(urlparse("http://[zz]/").is_err());
        assert!(urlparse("http://h:99999/").unwrap().port().is_err());
        assert!(urlparse("http://h:x/").unwrap().port().is_err());
        assert_eq!(unquote("a%20b%2Fc%zz"), "a b/c%zz");
        assert_eq!(unquote("%00"), "\0");
    }

    #[test]
    fn unescape_as_python() {
        assert_eq!(unescape("a &amp; b &lt;&#62;&#x3C; &amp"), "a & b <>< &");
        assert_eq!(unescape("&#128; &#0; &bogus;"), "\u{20ac} \u{fffd} &bogus;");
        assert_eq!(unescape_attr("x&amp=1&amp;y&#38;"), "x&amp=1&y&");
    }
}
