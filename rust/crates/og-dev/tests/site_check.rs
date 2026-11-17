//! `og dev site-check` on a small built site: green as written, and red with
//! the expected ERROR line for each fault planted in a copy of it.
use og_dev::site_check;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};

static NEXT: AtomicUsize = AtomicUsize::new(0);

/// A scratch folder holding `site/` and `zensical.toml`, removed on drop.
struct Site(PathBuf);

impl Drop for Site {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

const CONFIG: &str = r#"
[project]
site_name = "Demo"
site_url = "https://geiserx.github.io/OwnGenome/"
nav = [{ Home = "index.md" }, { Guide = "guide.md" }]

[project.theme]
variant = "modern"
features = ["navigation.tabs", "content.action.copy"]

[project.plugins.search]

[project.plugins.llmstxt.sections]
Docs = ["*.md"]

[project.plugins.redirects.redirect_maps]
"old.md" = "guide.md"
"#;

const NAV: &str = r#"<nav><a class="md-nav__link" href="/OwnGenome/">Home</a>
<a class="md-tabs__link" href="/OwnGenome/guide/">Guide</a></nav>"#;

fn index_html() -> String {
    format!(
        r##"<!doctype html>
<html><head><meta charset="utf-8"><title>Demo</title>
<link rel="canonical" href="https://geiserx.github.io/OwnGenome/">
<link rel="stylesheet" href="assets/stylesheets/modern/main.css">
<link rel="alternate" type="text/markdown" href="./index.md">
</head><body>{NAV}<a class="md-nav__link" href="#intro">Intro</a>
<div class="md-content"><article>
<button data-md-copy-url="index.md"></button>
<h1 id="intro">Intro</h1>
<p><a href='guide/#install'>install</a> <a href=guide/#top>top</a>
<a href="guide/#install:~:text=Install">directive</a>
<img src="images/logo.svg"> <a href="https://example.com/x.md">out</a>
<a href="mailto:a@b.c">mail</a> <a href="https://geiserx.github.io/Other/x/">other repo</a>
<a href="https://geiserx.github.io/OwnGenome/guide/">full URL</a> <a href="#intro">self</a>
<a href="guide/?q=1#install">query</a> <a href="caf&amp;e/">entity</a></p>
<svg><title>Drawing</title></svg>
<template><p id="hidden"></p></template>
<script>var s = '<a href="nowhere.html">';</script>
<!-- <a href="nowhere.html"> -->
<!-- n:other -->
</article></div></body></html>
"##
    )
}

fn guide_html() -> String {
    format!(
        r#"<html><head><title>Guide</title></head><body>{NAV}
<div class="md-content"><article>
<button data-md-copy-url="index.md"></button>
<h2 id="install">Install</h2><textarea>closed</textarea>
<a href="../">home</a></article></div></body></html>
"#
    )
}

const STUB: &str = r#"<html><head><title>Redirecting...</title>
<link rel="canonical" href="../guide/"></head>
<body><a href="../guide/">here</a></body></html>"#;

impl Site {
    /// The clean site: two pages, a redirect stub, a 404 page, the
    /// Markdown copies, an image, the modern stylesheet and llms.txt.
    fn new() -> Self {
        let root = std::env::temp_dir().join(format!(
            "og-site-check-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        let s = Site(root);
        s.write("zensical.toml", CONFIG);
        s.write("site/index.html", &index_html());
        s.write("site/index.md", "# Intro\n");
        s.write("site/guide/index.html", &guide_html());
        s.write("site/guide/index.md", "# Guide\n");
        s.write(
            "site/caf&e/index.html",
            "<html><body>asset page</body></html>",
        );
        s.write("site/old/index.html", STUB);
        s.write("site/404.html", "<html><body>not found</body></html>");
        s.write("site/images/logo.svg", "<svg/>");
        s.write("site/assets/stylesheets/modern/main.css", "");
        s.write("site/llms.txt", "# Demo\n");
        s
    }

    fn path(&self, rel: &str) -> PathBuf {
        self.0.join(rel)
    }

    fn write(&self, rel: &str, text: &str) {
        let p = self.path(rel);
        std::fs::create_dir_all(p.parent().unwrap()).unwrap();
        std::fs::write(p, text).unwrap();
    }

    fn edit(&self, rel: &str, from: &str, to: &str) {
        let text = std::fs::read_to_string(self.path(rel)).unwrap();
        assert!(text.contains(from), "{rel} has no {from:?}");
        self.write(rel, &text.replacen(from, to, 1));
    }

    fn remove(&self, rel: &str) {
        let p = self.path(rel);
        if p.is_dir() {
            std::fs::remove_dir_all(p).unwrap();
        } else {
            std::fs::remove_file(p).unwrap();
        }
    }

    fn run(&self) -> (u8, String) {
        run_on(&self.path("site"), &self.path("zensical.toml"))
    }
}

fn run_on(site: &Path, config: &Path) -> (u8, String) {
    let mut out = Vec::new();
    let mut err = Vec::new();
    let args = [
        site.to_str().unwrap().to_owned(),
        config.to_str().unwrap().to_owned(),
    ];
    let code = site_check::cli(&args, &mut out, &mut err);
    assert!(err.is_empty(), "{}", String::from_utf8_lossy(&err));
    (code, String::from_utf8(out).unwrap())
}

/// The planted fault turns the check red with this ERROR line.
fn red(site: &Site, line: &str) {
    let (code, out) = site.run();
    assert_eq!(code, 1, "no failure for {line:?}:\n{out}");
    assert!(
        out.lines()
            .any(|l| l.starts_with("ERROR ") && l.contains(line)),
        "no ERROR line with {line:?}:\n{out}"
    );
}

#[test]
fn clean_site_passes() {
    let site = Site::new();
    let (code, out) = site.run();
    assert_eq!(
        (code, out.as_str()),
        (0, "site-check: 2 pages, 0 errors, prefix /OwnGenome/\n")
    );
}

#[test]
fn dangling_links_and_images() {
    let cases: [(&str, &str, &str); 5] = [
        (
            "/OwnGenome/\">Home</a>",
            "/OwnGenome/\">Home</a><a class=\"md-nav__link\" href=\"missing/\">x</a>",
            "missing/ -> missing missing/index.html",
        ),
        (
            "images/logo.svg",
            "images/gone.png",
            "images/gone.png -> missing images/gone.png",
        ),
        (
            "/OwnGenome/guide/\">full",
            "/OwnGenome/gone/\">full",
            "https://geiserx.github.io/OwnGenome/gone/ -> missing gone/index.html",
        ),
        (
            "href=\"#intro\">self",
            "href=\"../../LICENSE\">self",
            "../../LICENSE leaves the built site",
        ),
        (
            "href=\"#intro\">self",
            "href=\"a%00b\">self",
            "index.html: 'a%00b' is not a valid URL",
        ),
    ];
    for (from, to, line) in cases {
        let site = Site::new();
        site.edit("site/index.html", from, to);
        red(&site, line);
    }
}

#[test]
fn missing_fragment() {
    let site = Site::new();
    site.edit("site/index.html", "guide/#install'", "guide/#nowhere'");
    red(
        &site,
        "index.html: guide/#nowhere -> no id 'nowhere' on guide/index.html",
    );
}

#[test]
fn markdown_source_link_in_the_body() {
    let site = Site::new();
    site.write("site/notes.md", "x");
    site.edit(
        "site/index.html",
        "href=\"#intro\">self",
        "href=\"notes.md\">self",
    );
    red(&site, "index.html: notes.md links a Markdown source");
}

#[test]
fn rooted_link_outside_the_prefix() {
    let site = Site::new();
    site.edit(
        "site/index.html",
        "href=\"#intro\">self",
        "href=\"/getting-started/\">self",
    );
    red(
        &site,
        "index.html: absolute link /getting-started/ (outside /OwnGenome/)",
    );
}

#[test]
fn stray_page() {
    let site = Site::new();
    site.write("site/stray/index.html", &guide_html());
    site.write("site/stray/index.md", "x");
    red(
        &site,
        "stray/index.html: page is in no navigation entry (stray page)",
    );
}

#[test]
fn published_python() {
    let site = Site::new();
    site.write("site/hooks/count.py", "");
    red(&site, "hooks/count.py: source file published");
    let site = Site::new();
    site.write("site/hooks/__pycache__/count.txt", "");
    red(&site, "hooks/__pycache__: source file published");
}

#[test]
fn count_token_left_in_a_page() {
    let site = Site::new();
    site.edit(
        "site/index.html",
        "<!-- n:other -->",
        "<!-- n:proyectos -->",
    );
    red(
        &site,
        "index.html: build token <!-- n:proyectos --> was not filled in",
    );
}

#[test]
fn raw_tags_that_hide_the_page() {
    let cases = [
        (
            "<p><a href='guide",
            "<p>a <title> b<a href='guide",
            "raw <title>",
        ),
        (
            "<p><a href='guide",
            "<p><textarea/><a href='guide",
            "raw <textarea>",
        ),
        (
            "<p><a href='guide",
            "<p><xmp>open<a href='guide",
            "raw <xmp>",
        ),
        (
            "<p><a href='guide",
            "<p><plaintext><a href='guide",
            "raw <plaintext>",
        ),
    ];
    for (from, to, line) in cases {
        let site = Site::new();
        site.edit("site/index.html", from, to);
        red(
            &site,
            &format!("index.html: {line} in the page body hides everything"),
        );
    }
}

#[test]
fn copy_buttons() {
    let site = Site::new();
    site.edit(
        "site/index.html",
        "data-md-copy-url=\"index.md\"",
        "data-md-copy-url=\"gone.md\"",
    );
    red(&site, "index.html: Copy as Markdown points at 'gone.md'");
    let site = Site::new();
    site.edit(
        "site/guide/index.html",
        "data-md-copy-url=\"index.md\"",
        "class=\"x\"",
    );
    red(
        &site,
        "guide/index.html: content.action.copy is on but the page has no Copy as Markdown button",
    );
}

#[test]
fn home_page_must_be_a_content_page() {
    let site = Site::new();
    site.edit("site/index.html", "class=\"md-content\"", "class=\"other\"");
    red(&site, "index.html is not a content page");
    let site = Site::new();
    site.remove("site/index.html");
    let (code, out) = site.run();
    assert_eq!(code, 1);
    assert!(
        out.ends_with("site: no index.html (did the build run?)\n"),
        "{out}"
    );
}

#[test]
fn site_url_and_its_prefix() {
    let site = Site::new();
    site.edit(
        "zensical.toml",
        "https://geiserx.github.io/OwnGenome/",
        "ftp://x/",
    );
    red(&site, "site_url 'ftp://x/' is not an absolute http(s) URL");
    let site = Site::new();
    site.edit(
        "zensical.toml",
        "https://geiserx.github.io/OwnGenome/",
        "http://[::1/",
    );
    red(&site, "site_url 'http://[::1/' is not a valid URL");
    // No site_url: the canonical link stands in, with the same prefix.
    let site = Site::new();
    site.edit(
        "zensical.toml",
        "site_url = \"https://geiserx.github.io/OwnGenome/\"",
        "",
    );
    let (code, out) = site.run();
    assert_eq!(code, 0, "{out}");
    assert!(out.ends_with("prefix /OwnGenome/\n"), "{out}");
    // Neither: an error, not a check that cannot fire.
    site.edit(
        "site/index.html",
        "<link rel=\"canonical\" href=\"https://geiserx.github.io/OwnGenome/\">",
        "",
    );
    red(
        &site,
        "no site_url in zensical.toml and no canonical link in index.html",
    );
}

#[test]
fn redirect_stubs() {
    let site = Site::new();
    site.remove("site/old");
    red(
        &site,
        "redirect_maps 'old.md': no redirect stub at old/index.html",
    );
    // A stub's own links are checked.
    let site = Site::new();
    site.edit(
        "site/old/index.html",
        "href=\"../guide/\"></head>",
        "href=\"../nowhere/\"></head>",
    );
    red(
        &site,
        "old/index.html: ../nowhere/ -> missing nowhere/index.html",
    );
    let site = Site::new();
    site.edit(
        "zensical.toml",
        "[project.plugins.redirects.redirect_maps]\n\"old.md\" = \"guide.md\"",
        "[project.plugins.redirects]\nredirect_maps = 3",
    );
    red(
        &site,
        "zensical.toml: redirects redirect_maps must be a table, got 3",
    );
}

#[test]
fn config_names_zensical_ignores() {
    let cases = [
        (
            "[project.plugins.search]",
            "[project.plugins.redirect]",
            "plugin 'redirect' is unknown to Zensical",
        ),
        (
            "\"navigation.tabs\"",
            "\"search.suggest\"",
            "feature 'search.suggest' is unknown to the Zensical UI",
        ),
        (
            "variant = \"modern\"",
            "variant = \"fancy\"",
            "theme.variant 'fancy' is not modern or classic",
        ),
        (
            "features = [\"navigation.tabs\", \"content.action.copy\"]",
            "features = \"navigation.tabs\"",
            "theme.features must be a list of strings, got 'navigation.tabs'",
        ),
    ];
    for (from, to, line) in cases {
        let site = Site::new();
        site.edit("zensical.toml", from, to);
        red(&site, line);
    }
    let site = Site::new();
    site.edit("site/index.html", "stylesheets/modern/", "stylesheets/");
    red(
        &site,
        "index.html does not load the modern stylesheet although theme.variant is modern",
    );
}

#[test]
fn plugins_in_list_form() {
    let site = Site::new();
    site.edit(
        "zensical.toml",
        "site_name = \"Demo\"",
        "site_name = \"Demo\"\nplugins = [\"search\", { \"material/llmstxt\" = {} }, 3]",
    );
    site.edit(
        "zensical.toml",
        "[project.plugins.search]\n\n[project.plugins.llmstxt.sections]\nDocs = [\"*.md\"]\n\n[project.plugins.redirects.redirect_maps]\n\"old.md\" = \"guide.md\"",
        "",
    );
    red(
        &site,
        "zensical.toml: plugins entry 3 is neither a name nor a one-entry table",
    );
    site.remove("site/llms.txt");
    red(
        &site,
        "llmstxt plugin is configured but site/llms.txt was not written",
    );
}

#[test]
fn llms_txt_only_when_the_plugin_is_on() {
    let site = Site::new();
    site.remove("site/llms.txt");
    red(
        &site,
        "llmstxt plugin is configured but site/llms.txt was not written",
    );
    site.edit(
        "zensical.toml",
        "[project.plugins.llmstxt.sections]",
        "[project.plugins.llmstxt]\nenabled = false\n[project.plugins.llmstxt.sections]",
    );
    let (code, out) = site.run();
    assert_eq!(code, 0, "{out}");
}

#[test]
fn config_that_cannot_be_read() {
    let site = Site::new();
    site.remove("zensical.toml");
    red(&site, "zensical.toml: config file not found");
    let site = Site::new();
    site.write("zensical.toml", "[project\n");
    red(&site, "zensical.toml: not valid TOML");
    let site = Site::new();
    site.write("zensical.toml", "[other]\nx = 1\n");
    red(
        &site,
        "zensical.toml: no [project] table, so no config check can run",
    );
}

#[test]
fn usage() {
    let mut out = Vec::new();
    let mut err = Vec::new();
    assert_eq!(site_check::cli(&[], &mut out, &mut err), 2);
    assert!(
        String::from_utf8(err)
            .unwrap()
            .starts_with("usage: og dev site-check")
    );
}
