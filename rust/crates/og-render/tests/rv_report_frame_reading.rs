//! The report's frame for a reader on a phone or a keyboard (og report
//! render), on a synthetic summary and on the demo summaries DEMO-001,
//! DEMO-002 and DEMO-003:
//!
//! - both reports say "Research software, not a clinical report" under the
//!   title, before the first chapter, and keep the full disclaimer at the end;
//! - the first element of the HTML report a keyboard reaches is the link
//!   "Skip to results" to `<main id="main" tabindex="-1">`, before the contents;
//! - the contents sit in a closed details element, and every table sits in a
//!   labelled scroll region (`.tscroll`), no two named alike, with the CSS
//!   that wraps long words in `main` and keeps a card inside a phone's width.
//!
//! The demo summaries are written by `og dev demo-sample` and `og summary
//! collect`, then rendered here by this crate. That needs the og binary: the
//! `OG` variable, else the `og` beside this test's own target folder (cargo
//! test --workspace builds it). Without one the demo test says SKIP, and
//! fails in CI. The browser check of the same reports at 320, 390 and 720 px
//! is rust/crates/og-browser/tests/rv_report_frame_reading.rs.
//!
//! Every test here fails on the report before this change (no research-only
//! line, no skip link, bare tables, the contents open); the PR shows that run.

use og_render::model::Summary;
use og_render::{html, text};
use std::path::{Path, PathBuf};
use std::process::Command;

const REPO: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../..");
const RESEARCH: &str = "Research software, not a clinical report";
/// The line under the HTML report's title.
const LINE: &str = "Research software, not a clinical report: confirm any finding in a certified laboratory before acting on it.";
/// The same line under the text report's title, on two lines.
const TEXT_LINES: [&str; 2] = [
    "  Research software, not a clinical report: confirm any finding in a",
    "  certified laboratory before acting on it.",
];
const SKIP: &str = "<a class=\"skip\" href=\"#main\">Skip to results</a>";
const REGION: &str = "<div class=\"tscroll\" role=\"region\" aria-label=\"";

/// A summary with the quality sections, a manifest (whose table the footer
/// shows) and one step not run.
fn synthetic() -> String {
    r#"{
 "schema_version": 2, "sample": "S1", "generated_utc": "2026-11-06 12:00 UTC",
 "run": {"started_utc": "2026-11-06 10:00 UTC", "declared_sex": "male"},
 "databases": {"clinvar_release": "2026-10-01", "cpsr_clinvar_release": null, "hla_database": "IPD-IMGT/HLA 3.60"},
 "manifest": [["run", "written_utc", "2026-11-06 12:00:00 UTC"], ["run", "written_by", "run-all.sh"],
              ["run", "git_commit", "0123456789abcdef0123456789abcdef01234567"], ["image", "PYTHON_IMAGE", "python:3.15"]],
 "sections": {
  "coverage": {"title": "Coverage (mosdepth)", "state": "ok", "step": "16b", "source": "qc/S1.mosdepth.summary.txt",
               "file_date": "2026-11-06 11:00:00", "note": null, "data": {"mean_depth": 31.25}},
  "variants": {"title": "Variant calling", "state": "ok", "step": "03", "source": "x", "file_date": null, "note": null,
               "data": {"total": 20000, "pass": 14545, "snps": 11111, "indels": 6667}}
 },
 "cyp2d6": null, "not_run": ["ClinVar screen"], "not_assessed": ["Repeat expansions longer than a read"]
}"#
    .to_owned()
}

fn render(json: &str) -> (String, String) {
    let s = Summary::parse(json).expect("a summary");
    (
        text::report(&s).expect("the text report"),
        html::report(&s).expect("the HTML report"),
    )
}

fn at(hay: &str, needle: &str, what: &str) -> usize {
    hay.find(needle)
        .unwrap_or_else(|| panic!("{what}: no {needle:?}"))
}

/// The text report: the research-only line under the title, before the
/// first chapter banner and the first section, and the disclaimer at the end.
fn text_problems(txt: &str, what: &str) -> Vec<String> {
    let mut out = Vec::new();
    let Some(line) = txt.find(RESEARCH) else {
        return vec![format!("{what}: the text report has no research-only line")];
    };
    let title = at(txt, "  GENOMICS ANALYSIS REPORT\n", what);
    if line < title {
        out.push(format!(
            "{what}: the research-only line comes before the title"
        ));
    }
    for first in ["\n####", "\n## "] {
        if txt.find(first).is_some_and(|i| i < line) {
            out.push(format!(
                "{what}: {first:?} comes before the research-only line"
            ));
        }
    }
    let joined: Vec<&str> = TEXT_LINES.iter().map(|l| l.trim()).collect();
    assert_eq!(joined.join(" "), LINE, "the text and HTML lines differ");
    let under = format!("  GENOMICS ANALYSIS REPORT\n{}\n", TEXT_LINES.join("\n"));
    if !txt.contains(&under) {
        out.push(format!(
            "{what}: the two research-only lines are not right under the title"
        ));
    }
    if !txt.contains("DISCLAIMER: This is NOT a clinical report.")
        || txt.find("DISCLAIMER:").is_some_and(|d| d < line)
    {
        out.push(format!("{what}: the full disclaimer is not at the end"));
    }
    out
}

/// The tags in `page` from its body on, as (offset, tag) pairs.
fn tags(page: &str) -> Vec<(usize, &str)> {
    let body = page.find("<body>").map_or(0, |i| i + "<body>".len());
    let mut out = Vec::new();
    let mut rest = &page[body..];
    let mut off = body;
    while let Some(i) = rest.find('<') {
        let Some(j) = rest[i..].find('>') else { break };
        out.push((off + i, &rest[i..i + j + 1]));
        off += i + j + 1;
        rest = &rest[i + j + 1..];
    }
    out
}

/// Whether a keyboard can reach `tag` (an opening tag): a link with an
/// address, a form control, a summary, or anything with a tabindex of 0 or more.
fn focusable(tag: &str) -> bool {
    let name: String = tag[1..]
        .chars()
        .take_while(|c| c.is_ascii_alphanumeric())
        .collect();
    match name.as_str() {
        "a" => tag.contains(" href="),
        "button" | "input" | "select" | "textarea" | "summary" => true,
        _ => tag
            .split_once(" tabindex=\"")
            .is_some_and(|(_, v)| !v.starts_with('-')),
    }
}

/// The HTML report: the research-only line in the header under the title and
/// before the first chapter; the skip link first; the contents folded; every
/// table in a labelled scroll region; the reading CSS in the page.
fn html_problems(page: &str, what: &str) -> Vec<String> {
    let mut out = Vec::new();
    let line = format!("  <p class=\"research-only\">{LINE}</p>");
    match (
        page.find(&line),
        page.find("  <h1>"),
        page.find("</header>"),
    ) {
        (Some(l), Some(h1), Some(end)) if h1 < l && l < end => {}
        _ => out.push(format!(
            "{what}: no research-only line under the title in the header"
        )),
    }
    if let (Some(l), Some(ch)) = (page.find(&line), page.find("<section class=\"chapter\""))
        && ch < l
    {
        out.push(format!(
            "{what}: a chapter comes before the research-only line"
        ));
    }
    if !page.contains("<strong>Disclaimer:</strong>") {
        out.push(format!("{what}: the full disclaimer is gone"));
    }

    let first = tags(page).into_iter().find(|(_, t)| focusable(t));
    let skip_end = page
        .find(SKIP)
        .map(|i| i + SKIP.len() - "Skip to results</a>".len());
    match first {
        Some((i, t)) if Some(i + t.len()) == skip_end => {}
        other => out.push(format!(
            "{what}: the first element a keyboard reaches is {:?}, not the skip link",
            other.map(|(_, t)| t)
        )),
    }
    if page.matches("<main id=\"main\" tabindex=\"-1\">").count() != 1
        || page.matches("<main").count() != 1
    {
        out.push(format!(
            "{what}: no single <main id=\"main\" tabindex=\"-1\">"
        ));
    }
    if !page
        .contains("<details class=\"toc-fold\">\n<summary>Contents</summary>\n<nav class=\"toc\"")
        || !page.contains("</nav>\n</details>")
    {
        out.push(format!(
            "{what}: the contents are not in a closed details element"
        ));
    }

    let tables = page.matches("<table").count();
    let mut names: Vec<&str> = Vec::new();
    for (i, _) in page.match_indices("<table") {
        let before = &page[..i];
        let label = before
            .rfind(REGION)
            .map(|r| &before[r + REGION.len()..])
            .and_then(|l| l.strip_suffix("\" tabindex=\"0\">"));
        if let Some(l) = label {
            // A screen reader lists the regions by name: two alike are one
            // to a listener (a card with several tables names each).
            if names.contains(&l) {
                out.push(format!("{what}: two scroll regions are named {l:?}"));
            }
            names.push(l);
        }
        if !label.is_some_and(|l| !l.is_empty() && !l.contains('"')) {
            let from = before.len().saturating_sub(120);
            out.push(format!(
                "{what}: a table outside a labelled scroll region: {:?}",
                &page[from..(i + 40).min(page.len())]
            ));
        }
    }
    if page.matches("</table></div>").count() != tables {
        out.push(format!(
            "{what}: not every table's scroll region closes after it"
        ));
    }
    for rule in [
        ".tscroll { overflow-x: auto; max-width: 100%; }",
        "main { overflow-wrap: anywhere; }",
        ".sec, .card { min-width: 0; }",
        ".skip:focus { top: 8px; }",
        "@media print { .toc-fold, .skip { display: none; } .tscroll { overflow: visible; } }",
    ] {
        if !page.contains(rule) {
            out.push(format!("{what}: the page's CSS lacks {rule:?}"));
        }
    }
    out
}

#[test]
fn the_synthetic_report_says_research_only_first_and_skips_to_the_results() {
    let (txt, page) = render(&synthetic());
    let mut problems = text_problems(&txt, "synthetic");
    problems.extend(html_problems(&page, "synthetic"));
    assert!(problems.is_empty(), "{}", problems.join("\n"));
    // The footer's manifest table is the one table here.
    assert_eq!(page.matches("<table").count(), 1, "{page}");
    assert!(page.contains(&format!("{REGION}Run manifest\" tabindex=\"0\"><table>")));
}

/// Each detector above fails on a page without the fix (main's frame), so a
/// green run means the frame holds, not that the detectors are blind.
#[test]
fn the_detectors_fail_on_the_old_frame() {
    let (txt, page) = render(&synthetic());
    let old_txt = TEXT_LINES
        .iter()
        .fold(txt.clone(), |t, l| t.replace(&format!("{l}\n"), ""));
    assert!(!text_problems(&old_txt, "old").is_empty());
    let late = old_txt.replace(
        "  DISCLAIMER:",
        &format!("{}\n  DISCLAIMER:", TEXT_LINES.join("\n")),
    );
    assert!(!text_problems(&late, "late").is_empty());

    let mangle = |from: &str, to: &str| {
        assert!(page.contains(from), "{from}");
        page.replacen(from, to, 1)
    };
    let plants = [
        mangle(&format!("  <p class=\"research-only\">{LINE}</p>\n"), ""),
        mangle(&format!("{SKIP}\n"), ""),
        mangle("<main id=\"main\" tabindex=\"-1\">", "<main>"),
        mangle(
            "<details class=\"toc-fold\">\n<summary>Contents</summary>\n",
            "",
        ),
        mangle(REGION, "<div class=\"x\" aria-label=\""),
        mangle("</table></div>", "</table>"),
        mangle("main { overflow-wrap: anywhere; }", ""),
        // A second region named like the first.
        mangle(
            "</table></div>",
            "</table></div>\n<div class=\"tscroll\" role=\"region\" aria-label=\"Run manifest\" tabindex=\"0\"><table></table></div>",
        ),
        // A link before the skip link is reached first.
        mangle("<body>\n", "<body>\n<a href=\"#x\">x</a>\n"),
    ];
    for (n, p) in plants.iter().enumerate() {
        assert!(!html_problems(p, "plant").is_empty(), "plant {n} passed");
    }
}

/// The og binary that writes the demo summaries, or `None` after a SKIP line
/// (a panic in CI).
fn og_or_skip() -> Option<PathBuf> {
    let from_env = std::env::var_os("OG")
        .filter(|v| !v.is_empty())
        .map(PathBuf::from);
    let beside = std::env::current_exe()
        .ok()
        .and_then(|e| Some(e.parent()?.parent()?.join("og")));
    let found = from_env.or(beside).filter(|p| p.is_file());
    if found.is_none() {
        assert!(
            std::env::var("CI").as_deref() != Ok("true"),
            "no og binary to write the demo summaries: CI must run this test (set OG)"
        );
        println!("[SKIP] no og binary to write the demo summaries (set OG)");
    }
    found
}

fn og(bin: &Path, args: &[&str]) {
    let out = Command::new(bin)
        .args(args)
        .current_dir(REPO)
        .output()
        .expect("og runs");
    assert!(
        out.status.success(),
        "og {args:?}: {}",
        String::from_utf8_lossy(&out.stderr)
    );
}

/// The summary.json of each demo person, written in `root`.
fn demo_summaries(bin: &Path, root: &Path) -> Vec<(&'static str, String)> {
    let r = root.to_str().unwrap();
    let mut out = Vec::new();
    for (sample, extra) in [
        ("DEMO-001", &[][..]),
        ("DEMO-002", &["--sex", "female", "--seed", "2"][..]),
        ("DEMO-003", &["--seed", "3"][..]),
    ] {
        let mut args = vec![
            "dev",
            "demo-sample",
            "--genome-dir",
            r,
            "--records",
            "3000",
            "--repo",
            REPO,
            "--sample",
            sample,
        ];
        args.extend(extra);
        og(bin, &args);
        let d = root.join(sample);
        let json = d.join("summary.json");
        og(
            bin,
            &[
                "summary",
                "collect",
                "--sample",
                sample,
                "--sample-dir",
                d.to_str().unwrap(),
                "--assets",
                &format!("{REPO}/assets"),
                "--out",
                json.to_str().unwrap(),
            ],
        );
        out.push((sample, std::fs::read_to_string(&json).unwrap()));
    }
    out
}

#[test]
fn the_demo_reports_say_research_only_first_and_put_every_table_in_a_region() {
    let Some(bin) = og_or_skip() else { return };
    let root = Path::new(env!("CARGO_TARGET_TMPDIR"))
        .join(format!("rv-report-frame-reading-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&root);
    std::fs::create_dir_all(&root).unwrap();
    let mut problems = Vec::new();
    for (sample, json) in demo_summaries(&bin, &root) {
        let (txt, page) = render(&json);
        problems.extend(text_problems(&txt, sample));
        problems.extend(html_problems(&page, sample));
        // The demo fills every card: many tables, so the region check is not empty.
        let tables = page.matches("<table").count();
        if tables < 10 {
            problems.push(format!("{sample}: only {tables} tables"));
        }
        for sid in ["cpic", "consumer-export"] {
            if !page.contains(&format!(" id=\"{sid}\">")) {
                problems.push(format!("{sample}: no #{sid} section"));
            }
        }
    }
    std::fs::remove_dir_all(&root).unwrap();
    assert!(problems.is_empty(), "{}", problems.join("\n"));
}
