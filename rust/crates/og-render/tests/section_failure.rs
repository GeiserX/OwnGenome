//! One section that cannot be written no longer stops the report: its
//! writer panics on purpose here, and both reports still come out, with
//! that section's error line in place and every other section written; the
//! verb says on stderr how many parts failed and still exits 0.

use og_render::chapters::{self, Failed, HtmlFn, TextFn};
use og_render::model::{HtmlPart, Malformed, Summary};
use og_render::{html, render, text};

const SUMMARY: &str = r#"{
 "schema_version": 2, "sample": "S1", "generated_utc": "2026-11-06 12:00 UTC",
 "run": {"started_utc": "2026-11-06 10:00 UTC", "declared_sex": "male"},
 "databases": {"clinvar_release": "2026-10-01", "cpsr_clinvar_release": null, "hla_database": null},
 "manifest": [["run", "git_commit", "0123456789abcdef0123456789abcdef01234567"]],
 "sections": {
  "coverage": {"title": "Coverage (mosdepth)", "state": "ok", "step": "16b", "source": "x", "file_date": null,
               "note": null, "data": {"mean_depth": 31.25}},
  "variants": {"title": "Variant calling", "state": "ok", "step": "03", "source": "x", "file_date": null, "note": null,
               "data": {"total": 20000, "pass": 14545, "snps": 11111, "indels": 6667}}
 },
 "cyp2d6": null, "not_run": [], "not_assessed": []
}"#;

const WHY: &str = "internal error: planted fault in the variants writer";

fn boom_text(_: &Summary) -> Result<Vec<String>, Malformed> {
    panic!("planted fault in the variants writer")
}

fn boom_html(_: &Summary) -> Result<HtmlPart, Malformed> {
    panic!("planted fault in the variants writer")
}

/// The real writers, with the variants section's replaced by one that panics.
fn planted(sid: &str) -> Option<(TextFn, HtmlFn)> {
    if sid == "variants" {
        Some((boom_text, boom_html))
    } else {
        chapters::writer(sid)
    }
}

fn frame_of<'a>(page: &'a str, sid: &str) -> &'a str {
    page.split_once(&format!("id=\"{sid}\">"))
        .map_or("", |x| x.1.split("</section>").next().unwrap_or(""))
}

fn want_failed() -> Vec<Failed> {
    vec![Failed {
        part: "variants".to_owned(),
        why: WHY.to_owned(),
    }]
}

#[test]
fn a_section_that_panics_writes_its_error_line_and_the_rest_of_the_text_report() {
    let s = Summary::parse(SUMMARY).expect("a summary");
    let (txt, failed) = text::report_using(&s, planted).expect("the text report");
    assert_eq!(failed, want_failed());
    assert!(
        txt.contains(&format!(
            "## Variant calling\n---\n  [could not be read: {WHY}]\n"
        )),
        "{txt}"
    );
    assert!(
        txt.contains("## Quality control"),
        "the other sections are written:\n{txt}"
    );
    assert!(
        txt.contains("## Not Assessed by This Pipeline"),
        "the report ends as always:\n{txt}"
    );
    let (clean, none) = text::report_using(&s, chapters::writer).expect("the text report");
    assert!(
        none.is_empty() && !clean.contains("could not be read"),
        "control: the real writers fail nothing"
    );
}

#[test]
fn a_section_that_panics_writes_its_error_line_and_the_rest_of_the_html_report() {
    let s = Summary::parse(SUMMARY).expect("a summary");
    let (page, failed) = html::report_using(&s, None, planted).expect("the HTML report");
    assert_eq!(failed, want_failed());
    let v = frame_of(&page, "variants");
    assert!(v.contains(&format!("could not be read: {WHY}")), "{v}");
    assert!(
        v.contains("<span class=\"badge badge-red\">Failed</span>"),
        "{v}"
    );
    assert!(
        page.contains("<a href=\"#variants\"><span class=\"dot dot-failed\""),
        "the contents say failed"
    );
    assert!(
        frame_of(&page, "qc").contains("Mean depth"),
        "the other sections are written"
    );
    assert!(page.trim_end().ends_with("</html>"), "the page is whole");
}

#[test]
fn the_verb_counts_the_failed_parts_on_stderr_and_still_writes_both_reports() {
    let d = std::env::temp_dir().join(format!("og-render-section-failure-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).expect("a temp folder");
    let sum = d.join("summary.json");
    std::fs::write(&sum, SUMMARY).expect("written");
    let (txt, page) = (d.join("r.txt"), d.join("r.html"));
    let args: Vec<String> = [
        "--summary",
        sum.to_str().unwrap_or(""),
        "-o",
        txt.to_str().unwrap_or(""),
        "-o",
    ]
    .iter()
    .map(|a| (*a).to_owned())
    .chain([page.to_string_lossy().into_owned()])
    .collect();
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = render::cli_using(&args, &mut out, &mut err, planted);
    let err = String::from_utf8(err).expect("utf-8");
    assert_eq!(code, 0, "{err}");
    assert_eq!(
        err,
        format!(
            "og report render: WARNING: 1 part of the report could not be written and says so in place: \
             variants ({WHY})\n"
        )
    );
    let t = std::fs::read_to_string(&txt).expect("the text report");
    let h = std::fs::read_to_string(&page).expect("the HTML report");
    assert!(t.contains(&format!("[could not be read: {WHY}]")) && t.contains("## Quality control"));
    assert!(frame_of(&h, "variants").contains(WHY) && frame_of(&h, "qc").contains("Mean depth"));
    let _ = std::fs::remove_dir_all(&d);
}
