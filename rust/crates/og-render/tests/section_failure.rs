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

/// A fresh folder for one test, removed when the test ends, even on a
/// failed assertion.
struct Scratch(std::path::PathBuf);

impl Scratch {
    fn new(name: &str) -> Scratch {
        let d = std::env::temp_dir().join(format!("og-render-{name}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).expect("a temp folder");
        Scratch(d)
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

fn args_of(parts: &[&std::path::Path]) -> Vec<String> {
    parts
        .iter()
        .map(|p| p.to_string_lossy().into_owned())
        .collect()
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
    let scratch = Scratch::new("section-failure");
    let d = &scratch.0;
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
}

#[test]
fn charts_that_cannot_be_made_leave_an_empty_charts_file_and_the_rest_of_the_run() {
    // A ClinVar section whose hits are not a table: its chart cannot be made.
    let body = SUMMARY.replacen(
        "\"sections\": {",
        "\"sections\": {\n  \"clinvar\": {\"title\": \"ClinVar\", \"state\": \"ok\", \"step\": \"07\", \"source\": \"x\", \
         \"file_date\": null, \"note\": null, \"data\": []},",
        1,
    );
    let scratch = Scratch::new("charts-failure");
    let d = &scratch.0;
    let sum = d.join("summary.json");
    std::fs::write(&sum, &body).expect("written");
    let (page, charts) = (d.join("r.html"), d.join("charts.json"));
    let mut args = args_of(&[
        std::path::Path::new("--summary"),
        &sum,
        std::path::Path::new("-o"),
        &page,
    ]);
    args.extend(args_of(&[std::path::Path::new("--charts-json"), &charts]));
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = render::cli_using(&args, &mut out, &mut err, chapters::writer);
    let err = String::from_utf8(err).expect("utf-8");
    assert_eq!(code, 0, "{err}");
    assert!(
        err.contains("og report render: WARNING: ") && err.contains("charts ("),
        "{err}"
    );
    let c = std::fs::read_to_string(&charts).expect("the charts file is written");
    assert!(
        c.contains("\"charts\": []") && c.contains("\"sample\": \"S1\""),
        "{c}"
    );
    let h = std::fs::read_to_string(&page).expect("the HTML report");
    assert!(
        frame_of(&h, "qc").contains("Mean depth"),
        "the other sections are written"
    );
}

#[test]
fn a_search_index_that_cannot_be_built_leaves_an_empty_index_and_the_rest_of_the_run() {
    // A section that is not an object, and sections that are not an
    // object: the index cannot be built from either.
    let cases = [
        (
            SUMMARY.replacen(
                "\"sections\": {",
                "\"sections\": {\n  \"clinvar\": \"x\",",
                1,
            ),
            "unreadable 1",
        ),
        (
            SUMMARY
                .split_once("\"sections\": {")
                .map(|(head, rest)| {
                    let tail = rest.split_once("\n },").map_or("", |x| x.1);
                    format!("{head}\"sections\": [1],{tail}")
                })
                .unwrap_or_default(),
            "none readable",
        ),
    ];
    for (n, (body, counted)) in cases.iter().enumerate() {
        let scratch = Scratch::new(&format!("index-failure-{n}"));
        let d = &scratch.0;
        let sum = d.join("summary.json");
        std::fs::write(&sum, body).expect("written");
        let (txt, page) = (d.join("r.txt"), d.join("r.html"));
        let (index, charts) = (d.join("index.json"), d.join("charts.json"));
        let flag = |f: &str| std::path::PathBuf::from(f);
        let args = args_of(&[
            &flag("--summary"),
            &sum,
            &flag("-o"),
            &txt,
            &flag("-o"),
            &page,
            &flag("--search-index"),
            &index,
            &flag("--charts-json"),
            &charts,
        ]);
        let (mut out, mut err) = (Vec::new(), Vec::new());
        let code = render::cli_using(&args, &mut out, &mut err, chapters::writer);
        let (out, err) = (
            String::from_utf8(out).expect("utf-8"),
            String::from_utf8(err).expect("utf-8"),
        );
        assert_eq!(code, 0, "case {n}: {err}");
        assert!(
            err.contains("og report render: WARNING: ") && err.contains("terms ("),
            "case {n}: {err}"
        );
        assert!(
            out.starts_with("Sections: ") && out.contains(counted),
            "case {n}: {out}"
        );
        let i = std::fs::read_to_string(&index).expect("the search index is written");
        assert!(
            i.contains("\"schema_version\": 1")
                && i.contains("\"sample\": \"S1\"")
                && i.contains("\"terms\": {}"),
            "case {n}: {i}"
        );
        assert!(charts.is_file(), "case {n}: the charts file is written");
        assert!(txt.is_file(), "case {n}: the text report is written");
        let h = std::fs::read_to_string(&page).expect("the HTML report");
        assert!(
            h.trim_end().ends_with("</html>"),
            "case {n}: the page is whole"
        );
    }
}

/// In the child process only (see the test below): installs og's panic
/// hook, panics inside `guarded`, then on a worker thread outside it.
const HOOK_CHILD: &str = "OG_TEST_SECTION_FAILURE_HOOK_CHILD";

#[test]
fn hook_child_entry() {
    if std::env::var_os(HOOK_CHILD).is_none() {
        return;
    }
    chapters::abort_outside_guarded();
    let caught = chapters::guarded(|| -> Result<(), Malformed> { panic!("planted in a section") });
    assert!(caught.is_err());
    println!("the section's panic was caught");
    // A worker whose panic only ends its thread: with unwinding the caller
    // reads that as the end of the work and goes on.
    let _ = std::thread::spawn(|| panic!("planted in a worker")).join();
    println!("the worker's panic was missed");
    std::process::exit(0);
}

#[test]
fn og_s_panic_hook_spares_a_report_section_and_aborts_on_a_worker_thread() {
    let out = std::process::Command::new(std::env::current_exe().expect("the test binary"))
        .args([
            "hook_child_entry",
            "--exact",
            "--nocapture",
            "--test-threads=1",
        ])
        .env(HOOK_CHILD, "1")
        .output()
        .expect("the child runs");
    let so = String::from_utf8_lossy(&out.stdout);
    let se = String::from_utf8_lossy(&out.stderr);
    assert!(so.contains("the section's panic was caught"), "{so}\n{se}");
    assert!(se.contains("planted in a worker"), "{so}\n{se}");
    assert!(
        !so.contains("the worker's panic was missed") && !out.status.success(),
        "a panic outside a section must abort: {:?}\n{so}\n{se}",
        out.status
    );
    #[cfg(unix)]
    {
        use std::os::unix::process::ExitStatusExt;
        assert_eq!(out.status.signal(), Some(6), "{so}\n{se}");
    }
}
