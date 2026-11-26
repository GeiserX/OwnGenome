//! The docs theme on this checkout: the current page in the side menu stays
//! readable, and code blocks carry no empty line anchors.
//!
//! 1. Zensical paints the active navigation item in `--md-accent-fg-color`
//!    on a background of `--md-accent-fg-color--transparent` laid over the
//!    page. For each colour scheme, the accent that scheme ends up with (its
//!    own block in `docs/stylesheets/theme.css`, else `:root`) must reach
//!    WCAG AA for normal text, 4.5:1, against that background. The reverse
//!    holds too: where Zensical fills a control with the accent (the "Back to
//!    top" button and the side menu button on hover or focus, a linked tag,
//!    the "Copied" note), the text on it, `--md-accent-bg-color` or the page
//!    colour, must reach 4.5:1 on the accent.
//! 2. `zensical.toml` must not turn on `pymdownx.highlight`'s
//!    `anchor_linenums`. The site shows no line numbers, so each code line
//!    would get an empty link with no name that still takes a Tab stop.

use std::fs;
use std::path::{Path, PathBuf};

fn repo() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../..")
}

fn read(rel: &str) -> String {
    fs::read_to_string(repo().join(rel)).unwrap_or_else(|e| panic!("{rel}: {e}"))
}

/// The CSS with every `/* ... */` comment removed.
fn strip_comments(css: &str) -> String {
    let mut out = String::with_capacity(css.len());
    let mut rest = css;
    while let Some(start) = rest.find("/*") {
        out.push_str(&rest[..start]);
        rest = match rest[start + 2..].find("*/") {
            Some(end) => &rest[start + 2 + end + 2..],
            None => "",
        };
    }
    out.push_str(rest);
    out
}

/// The body of the first rule whose selector is exactly `selector`.
fn block<'a>(css: &'a str, selector: &str) -> Option<&'a str> {
    let mut rest = css;
    while let Some(open) = rest.find('{') {
        let head = rest[..open].rsplit('}').next().unwrap_or("").trim();
        let close = open + rest[open..].find('}')?;
        if head == selector {
            return Some(&rest[open + 1..close]);
        }
        rest = &rest[close + 1..];
    }
    None
}

/// The value of custom property `name` in a rule body, the last one if set twice.
fn prop(body: &str, name: &str) -> Option<String> {
    body.split(';')
        .filter_map(|d| d.split_once(':'))
        .filter(|(k, _)| k.trim() == name)
        .map(|(_, v)| v.trim().to_owned())
        .next_back()
}

/// `#rrggbb` or `#rrggbbaa` as red, green, blue (0-255) and alpha (0-1).
fn hex(v: &str) -> [f64; 4] {
    let h = v
        .strip_prefix('#')
        .unwrap_or_else(|| panic!("not a hex colour: {v}"));
    assert!(h.len() == 6 || h.len() == 8, "not #rrggbb(aa): {v}");
    let byte = |i: usize| {
        f64::from(u8::from_str_radix(&h[i..i + 2], 16).unwrap_or_else(|_| panic!("bad hex: {v}")))
    };
    let a = if h.len() == 8 { byte(6) / 255.0 } else { 1.0 };
    [byte(0), byte(2), byte(4), a]
}

/// `top` laid over the opaque `under`.
fn over(top: [f64; 4], under: [f64; 4]) -> [f64; 4] {
    let a = top[3];
    [
        top[0] * a + under[0] * (1.0 - a),
        top[1] * a + under[1] * (1.0 - a),
        top[2] * a + under[2] * (1.0 - a),
        1.0,
    ]
}

/// WCAG 2 relative luminance of an opaque colour.
fn luminance(c: [f64; 4]) -> f64 {
    let lin = |v: f64| {
        let s = v / 255.0;
        if s <= 0.040_45 {
            s / 12.92
        } else {
            ((s + 0.055) / 1.055).powf(2.4)
        }
    };
    0.2126 * lin(c[0]) + 0.7152 * lin(c[1]) + 0.0722 * lin(c[2])
}

fn contrast(a: [f64; 4], b: [f64; 4]) -> f64 {
    let (la, lb) = (luminance(a), luminance(b));
    (la.max(lb) + 0.05) / (la.min(lb) + 0.05)
}

/// The scheme's accent and the active item's background, as the theme sets them.
/// `page_default` is the page colour Zensical uses when the block sets none;
/// `None` means the block must set its own.
fn active_item(css: &str, scheme: &str, page_default: Option<&str>) -> (String, [f64; 4], f64) {
    let root = block(css, ":root").expect("theme.css has a :root block");
    let sel = format!("[data-md-color-scheme=\"{scheme}\"]");
    let own = block(css, &sel).unwrap_or_else(|| panic!("theme.css has no {sel} block"));
    let pick = |name: &str| prop(own, name).or_else(|| prop(root, name));
    let accent = pick("--md-accent-fg-color").expect("an accent colour");
    let tint = pick("--md-accent-fg-color--transparent").expect("a transparent accent");
    let page = pick("--md-default-bg-color")
        .or_else(|| page_default.map(str::to_owned))
        .unwrap_or_else(|| panic!("{sel} sets no --md-default-bg-color"));
    let bg = over(hex(&tint), hex(&page));
    let ratio = contrast(hex(&accent), bg);
    (accent, bg, ratio)
}

#[test]
fn the_active_menu_item_reads_in_both_schemes() {
    let css = strip_comments(&read("docs/stylesheets/theme.css"));
    // The default scheme keeps Zensical's white page; the slate block sets
    // its own page colour.
    let mut low = Vec::new();
    for (scheme, page) in [("slate", None), ("default", Some("#ffffff"))] {
        let (accent, bg, ratio) = active_item(&css, scheme, page);
        if ratio < 4.5 {
            low.push(format!(
                "{scheme}: {accent} on rgb({:.0},{:.0},{:.0}) is {ratio:.2}:1",
                bg[0], bg[1], bg[2]
            ));
        }
    }
    assert!(
        low.is_empty(),
        "the active menu item is below WCAG AA (4.5:1): {}",
        low.join("; ")
    );
}

/// The colour a scheme ends up with for `name`: its own block, then `:root`,
/// then `fallback` (Zensical's own value when the theme sets none).
fn effective(css: &str, scheme: &str, name: &str, fallback: Option<&str>) -> String {
    let root = block(css, ":root").expect("theme.css has a :root block");
    let sel = format!("[data-md-color-scheme=\"{scheme}\"]");
    let own = block(css, &sel).unwrap_or_else(|| panic!("theme.css has no {sel} block"));
    prop(own, name)
        .or_else(|| prop(root, name))
        .or_else(|| fallback.map(str::to_owned))
        .unwrap_or_else(|| panic!("{sel} sets no {name}"))
}

#[test]
fn text_on_the_accent_reads_in_both_schemes() {
    let css = strip_comments(&read("docs/stylesheets/theme.css"));
    let mut low = Vec::new();
    for (scheme, page) in [("slate", None), ("default", Some("#ffffff"))] {
        let accent = effective(&css, scheme, "--md-accent-fg-color", None);
        // Zensical sets --md-accent-bg-color to #fff for a custom accent.
        // The "Back to top" button, the side menu button and linked tags
        // paint it on the accent on hover and focus.
        let on_accent = effective(&css, scheme, "--md-accent-bg-color", Some("#ffffff"));
        // The "Copied" note writes the page colour on the accent.
        let page = effective(&css, scheme, "--md-default-bg-color", page);
        for (what, text) in [("--md-accent-bg-color", on_accent), ("page colour", page)] {
            let ratio = contrast(hex(&text), hex(&accent));
            if ratio < 4.5 {
                low.push(format!(
                    "{scheme}: {what} {text} on {accent} is {ratio:.2}:1"
                ));
            }
        }
    }
    assert!(
        low.is_empty(),
        "text on the accent is below WCAG AA (4.5:1): {}",
        low.join("; ")
    );
}

#[test]
fn the_contrast_maths_matches_known_pairs() {
    // Black on white is 21:1, and a colour on itself is 1:1.
    assert!((contrast(hex("#000000"), hex("#ffffff")) - 21.0).abs() < 1e-9);
    assert!((contrast(hex("#2563eb"), hex("#2563eb")) - 1.0).abs() < 1e-9);
    // The two failing pairs the old theme had: #2563eb on its own 10% tint
    // over each page.
    let slate = over(hex("#2563eb1a"), hex("#1c1917"));
    let light = over(hex("#2563eb1a"), hex("#ffffff"));
    let (s, l) = (
        contrast(hex("#2563eb"), slate),
        contrast(hex("#2563eb"), light),
    );
    assert!((3.10..3.15).contains(&s), "slate {s}");
    assert!((4.45..4.5).contains(&l), "default {l}");
}

#[test]
fn code_lines_carry_no_empty_anchors() {
    let toml = read("zensical.toml");
    let on: Vec<&str> = toml
        .lines()
        .map(str::trim)
        .filter(|l| !l.starts_with('#'))
        .filter(|l| l.contains("anchor_linenums"))
        .filter(|l| {
            let bare: String = l
                .split('#')
                .next()
                .unwrap_or("")
                .split_whitespace()
                .collect();
            bare.contains("anchor_linenums=true")
        })
        .collect();
    assert!(
        on.is_empty(),
        "zensical.toml turns on anchor_linenums, which adds an empty, unnamed link to every code line: {on:?}"
    );
}
