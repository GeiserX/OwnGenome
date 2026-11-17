//! The report's frame and its first sections (og report render), on
//! synthetic summaries: the chapters, the page's safety rules (no fetch, no
//! script, every contents link resolves), the run facts with the OwnGenome
//! release and the database re-check, the new-since badges, escaping, the
//! caveats under their card titles, and a malformed section that no longer
//! stops the report. Each detector is shown red on a planted bad input, and
//! three rules are guarded by planted faults of the crate (`fault-*`
//! features; the controls in .github/workflows/guard.yml run them).

use og_render::chapters::{CHAPTERS, Chapter, Sec, key_problems, summary_keys};
use og_render::model::Summary;
use og_render::{html, render, text};

const COMMIT: &str = "0123456789abcdef0123456789abcdef01234567";
const BADGE: &str = "New since 2026-01-01T00:00:00Z";

/// A summary of a run whose quality steps and variant calling ran, with
/// `extra` top-level keys and `sections` merged in (JSON object bodies).
fn summary(extra: &str, manifest: &str, sections: &str) -> String {
    format!(
        r#"{{
 "schema_version": 2, "sample": "S1", "generated_utc": "2026-11-06 12:00 UTC",
 "run": {{"started_utc": "2026-11-06 10:00 UTC", "declared_sex": "male"}},
 "databases": {{"clinvar_release": "2026-10-01", "cpsr_clinvar_release": null, "hla_database": "IPD-IMGT/HLA 3.60"}},
 "manifest": [{manifest}],
 "sections": {{
  "coverage": {{"title": "Coverage (mosdepth)", "state": "ok", "step": "16b", "source": "qc/S1.mosdepth.summary.txt",
               "file_date": "2026-11-06 11:00:00", "note": null, "data": {{"mean_depth": 31.25}}}},
  "sex_check": {{"title": "Sex check (indexcov)", "state": "ok", "step": "16", "source": "x", "file_date": "2026-11-06 11:30:00",
                "note": null, "data": {{"inferred_sex": "male"}}}},
  "sample_qc": {{"title": "Sample identity", "state": "ok", "step": "33", "source": "x", "file_date": null, "note": null,
                "data": {{"inferred_sex": "male", "sex_check": "ok", "contamination": "ok", "freemix": 0.001,
                         "freemix_warn_above": 0.02}}}},
  "read_count": {{"title": "Reads kept by alignment", "state": "missing", "step": "02", "source": null, "file_date": null,
                 "note": null, "data": {{}}}},
  "variants": {{"title": "Variant calling", "state": "ok", "step": "03", "source": "x", "file_date": null, "note": null,
               "data": {{"total": 20000, "pass": 14545, "snps": 11111, "indels": 6667}}}}{sections}
 }},
 "cyp2d6": null, "not_run": ["ClinVar screen"], "not_assessed": ["Repeat expansions longer than a read"]{extra}
}}"#
    )
}

fn base_manifest() -> String {
    format!(
        r#"["run", "written_utc", "2026-11-06 12:00:00 UTC"], ["run", "written_by", "run-all.sh"],
           ["run", "git_commit", "{COMMIT}"], ["image", "PYTHON_IMAGE", "python:3.15"]"#
    )
}

fn render_both(json: &str) -> (String, String) {
    let s = Summary::parse(json).expect("a summary");
    (
        text::report(&s).expect("the text report"),
        html::report(&s).expect("the HTML report"),
    )
}

fn plain() -> (String, String) {
    render_both(&summary("", &base_manifest(), ""))
}

/// The frame of section `sid` in the HTML report, empty when it has none.
fn frame_of<'a>(page: &'a str, sid: &str) -> &'a str {
    let open = format!("id=\"{sid}\">");
    page.split_once(&open).map_or("", |(_, rest)| {
        rest.split("</section>").next().unwrap_or("")
    })
}

/// Every `name="value"` value in `page`.
fn attrs<'a>(page: &'a str, name: &str) -> Vec<&'a str> {
    let open = format!("{name}=\"");
    page.match_indices(&open)
        .filter_map(|(i, _)| {
            let rest = &page[i + open.len()..];
            rest.find('"').map(|j| &rest[..j])
        })
        .collect()
}

/// Ids given twice, contents links that reach no id, or no contents links.
fn anchor_problems(page: &str) -> Vec<String> {
    let ids: Vec<&str> = page
        .match_indices(" id=\"")
        .map(|(i, _)| {
            let rest = &page[i + 5..];
            &rest[..rest.find('"').unwrap_or(0)]
        })
        .collect();
    let toc = page.split("</nav>").next().unwrap_or("");
    let links: Vec<&str> = attrs(toc, "href")
        .into_iter()
        .filter_map(|h| h.strip_prefix('#'))
        .collect();
    let mut out: Vec<String> = Vec::new();
    let mut sorted = ids.clone();
    sorted.sort_unstable();
    sorted.dedup();
    for i in sorted {
        if ids.iter().filter(|x| **x == i).count() > 1 {
            out.push(format!("id twice: {i}"));
        }
    }
    for l in &links {
        if !ids.contains(l) {
            out.push(format!("contents link to nothing: {l}"));
        }
    }
    if links.is_empty() {
        out.push("no contents links".to_owned());
    }
    out
}

/// Every http(s) src or href outside a link marked as opening a website.
fn external_problems(page: &str) -> Vec<String> {
    let mut rest = page.to_owned();
    while let Some(i) = rest.find("<span class=\"offsite\"><a href=\"") {
        let Some(j) = rest[i..].find("(opens a website") else {
            break;
        };
        let Some(k) = rest[i + j..].find("</span></span>") else {
            break;
        };
        rest.replace_range(i..i + j + k + "</span></span>".len(), "");
    }
    ["href", "src"]
        .iter()
        .flat_map(|n| attrs(&rest, n))
        .filter(|v| v.starts_with("http:") || v.starts_with("https:"))
        .map(str::to_owned)
        .collect()
}

/// The section anchors of the HTML report, in page order.
fn shown_sections(page: &str) -> Vec<&str> {
    page.lines()
        .filter(|l| l.starts_with("  <section class=\"sec"))
        .flat_map(|l| attrs(l, "id"))
        .collect()
}

/// The chapter titles of the text report, in order.
fn text_chapters(txt: &str) -> Vec<String> {
    txt.lines()
        .filter_map(|l| l.strip_prefix("#  "))
        .filter_map(|l| l.split_once(". ").map(|(_, t)| t))
        .filter_map(|t| t.split_once(": ").map(|(t, _)| t.to_owned()))
        .collect()
}

#[test]
fn the_release_stands_next_to_the_commit_in_both_reports() {
    let with = format!(
        r#"{}, ["run", "pipeline_version", "1.0.0"]"#,
        base_manifest()
    );
    let (txt, page) = render_both(&summary("", &with, ""));
    assert!(
        txt.contains(&format!("(OwnGenome 1.0.0, pipeline commit {COMMIT},")),
        "{txt}"
    );
    assert!(page.contains("OwnGenome release: 1.0.0 &nbsp;|&nbsp;"));
    assert!(page.contains("Pipeline commit: 0123456789abcdef"));
    // Red control: without the row neither report names a release.
    let (txt, page) = plain();
    assert!(!txt.contains("OwnGenome 1.") && txt.contains("(pipeline commit 0123456789abcdef"));
    assert!(!page.contains("OwnGenome release:"));
    // No manifest: the footer says so and the text report has no manifest line.
    let (txt, page) = render_both(&summary("", "", ""));
    assert!(!txt.contains("Run manifest:"));
    assert!(page.contains("<p>No run_manifest.tsv in the sample folder.</p>"));
    assert!(page.contains("Pipeline commit: unknown"));
}

#[test]
fn every_summary_key_is_in_exactly_one_section() {
    let keys: Vec<&str> = summary_keys().map(|s| s.key).collect();
    assert_eq!(keys.len(), 48);
    assert_eq!(key_problems(&CHAPTERS, &keys), Vec::<String>::new());
    // Red controls: a key taken out, and a key in two sections, are caught.
    let cut: Vec<&str> = keys.iter().copied().filter(|k| *k != "telomere").collect();
    assert_eq!(
        key_problems(&CHAPTERS, &cut),
        ["unknown telomere".to_owned()]
    );
    let mut twice = CHAPTERS.to_vec();
    twice.push(Chapter {
        id: "x",
        title: "X",
        question: "?",
        sections: &[Sec {
            id: "x",
            title: "X",
            keys: &["clinvar"],
        }],
    });
    assert_eq!(key_problems(&twice, &keys), ["twice clinvar".to_owned()]);
    let mut short = CHAPTERS.to_vec();
    short.pop();
    let mut missing = key_problems(&short, &keys);
    missing.sort();
    assert_eq!(missing, ["missing consumer_export", "missing viral"]);
}

#[test]
fn the_page_fetches_nothing_and_runs_no_script() {
    let (_, page) = plain();
    assert!(page.contains(
        "<meta http-equiv=\"Content-Security-Policy\" content=\"default-src 'none'; img-src data:; \
         style-src 'unsafe-inline'\">"
    ));
    assert!(!page.to_lowercase().contains("<script"));
    assert_eq!(external_problems(&page), Vec::<String>::new());
    assert!(
        page.contains("<span class=\"offsite\"><a href=\"https://github.com/GeiserX/OwnGenome\">")
    );
    // Red controls: an unmarked link, an image from off the machine, and a
    // link marked the wrong way are each caught.
    for bad in [
        "<a href=\"https://example.org/x\">x</a>",
        "<img src=\"https://example.org/i.png\">",
        "<a class=\"offsite\" href=\"https://example.org/y\">y</a>",
    ] {
        assert!(
            !external_problems(&format!("{page}{bad}")).is_empty(),
            "{bad}"
        );
    }
}

#[test]
fn ids_are_unique_and_every_contents_link_reaches_one() {
    let (_, page) = plain();
    assert_eq!(anchor_problems(&page), Vec::<String>::new());
    // Red control: a second id and a dead contents link are caught.
    let broken = page
        .replace("id=\"variants\"", "id=\"qc\"")
        .replace("href=\"#not-run\"", "href=\"#gone\"");
    assert_eq!(
        anchor_problems(&broken).len(),
        3,
        "{:?}",
        anchor_problems(&broken)
    );
}

#[test]
fn both_reports_follow_the_chapter_order() {
    let (txt, page) = plain();
    let sids: Vec<&str> = CHAPTERS
        .iter()
        .flat_map(|c| c.sections.iter().map(|s| s.id))
        .collect();
    let shown = shown_sections(&page);
    // The sections ported so far, in chapter order, then the limits chapter.
    let n = shown.len();
    assert!(
        n >= 4 && shown[..2] == ["qc", "variants"] && shown[n - 2..] == ["not-run", "not-assessed"],
        "{shown:?}"
    );
    let in_order: Vec<&str> = sids.iter().copied().filter(|s| shown.contains(s)).collect();
    assert_eq!(&shown[..n - 2], in_order.as_slice());
    assert_eq!(text_chapters(&txt), ["QUALITY"]);
    let q = txt.split("#  1. QUALITY").nth(1).unwrap_or("");
    assert!(q.find("## Quality control") < q.find("## Variant Calling (DeepVariant)"));
    assert!(q.contains("## Variant Calling (DeepVariant)"));
    // Red control: the order check sees a page whose sections were swapped.
    let swapped = page
        .replace("id=\"qc\"", "id=\"tmp\"")
        .replace("id=\"variants\"", "id=\"qc\"")
        .replace("id=\"tmp\"", "id=\"variants\"");
    assert_ne!(shown_sections(&swapped)[..n - 2], in_order[..]);
}

#[test]
fn the_text_quality_block_carries_the_badge_when_coverage_changed() {
    let since =
        r#", "new_since": {"previous_utc": "2026-01-01T00:00:00Z", "sections": ["coverage"]}"#;
    let (txt, page) = render_both(&summary(since, &base_manifest(), ""));
    let qc = txt.split("## Quality control").nth(1).unwrap_or("");
    assert!(
        qc.starts_with(&format!("\n---\n  [{BADGE}]\n")),
        "{}",
        &qc[..qc.len().min(200)]
    );
    assert_eq!(txt.matches(BADGE).count(), 1);
    assert!(frame_of(&page, "qc").contains(BADGE));
}

#[test]
fn the_badge_sits_only_on_the_section_whose_result_changed() {
    let since =
        r#", "new_since": {"previous_utc": "2026-01-01T00:00:00Z", "sections": ["variants"]}"#;
    let (txt, page) = render_both(&summary(since, &base_manifest(), ""));
    let marked = |page: &str| -> Vec<&'static str> {
        ["qc", "variants"]
            .into_iter()
            .filter(|s| frame_of(page, s).contains(BADGE))
            .collect()
    };
    assert_eq!(marked(&page), ["variants"]);
    assert_eq!(page.matches(BADGE).count(), 1);
    assert!(txt.contains(&format!(
        "## Variant Calling (DeepVariant)\n---\n  [{BADGE}]\n"
    )));
    // Red control: a badge on every section is caught.
    let every = r#", "new_since": {"previous_utc": "2026-01-01T00:00:00Z", "sections": ["coverage", "variants"]}"#;
    let (_, all) = render_both(&summary(every, &base_manifest(), ""));
    assert_ne!(marked(&all), ["variants"]);
    // A section gone missing keeps its badge on its HTML card, and the text
    // report leaves the section out.
    let gone = summary(since, &base_manifest(), "").replace(
        "\"state\": \"ok\", \"step\": \"03\"",
        "\"state\": \"missing\", \"step\": \"03\"",
    );
    let (txt, page) = render_both(&gone);
    let v = frame_of(&page, "variants");
    assert!(v.contains(BADGE) && v.contains("Not run"), "{v}");
    assert!(!txt.contains("## Variant Calling") && !txt.contains(BADGE));
}

#[test]
fn no_badge_without_a_previous_run() {
    let (txt, page) = plain();
    assert!(!txt.contains("New since") && !page.contains("New since"));
    let none = r#", "new_since": {"previous_utc": null, "sections": ["variants"]}"#;
    let (txt, _) = render_both(&summary(none, &base_manifest(), ""));
    assert!(txt.contains("  [New since the previous run]"));
}

#[test]
fn the_recheck_is_one_line_of_run_facts_in_both_reports_and_no_caveat() {
    let rc = r#", "recheck": {"checked_utc": "2026-02-02T00:00:00Z", "affected_steps": ["06", "24"], "rerun": null,
      "items": [{"item": "clinvar_file_date", "label": "ClinVar", "class": "refresh-available"},
                {"item": "VEP_CACHE_RELEASE", "label": "VEP", "class": "no-newer-release"}]}"#;
    let (txt, page) = render_both(&summary(rc, &base_manifest(), ""));
    let line = "Database re-check: 2026-02-02T00:00:00Z: ClinVar refresh-available";
    let (_, plain_page) = plain();
    assert!(page.split("</header>").next().unwrap_or("").contains(line));
    assert!(
        txt.split(&"#".repeat(80))
            .next()
            .unwrap_or("")
            .contains(line)
    );
    assert_eq!(
        page.matches("<p class=\"caveat\">").count(),
        plain_page.matches("<p class=\"caveat\">").count()
    );
    let nothing = r#", "recheck": {"checked_utc": null, "items": []}"#;
    let (txt, _) = render_both(&summary(nothing, &base_manifest(), ""));
    assert!(txt.contains("  Database re-check: unknown date: nothing newer than this run\n"));
}

#[test]
fn summary_text_is_escaped_in_the_page() {
    let json = summary("", &base_manifest(), "")
        .replace("\"sample\": \"S1\"", r#""sample": "A&B <x> \"q\" 's""#);
    let (txt, page) = render_both(&json);
    assert!(
        page.contains("Sample: <strong>A&amp;B &lt;x&gt; &quot;q&quot; &#x27;s</strong>"),
        "{}",
        page.split("<main>").nth(1).unwrap_or("")
    );
    assert!(txt.contains("  Sample: A&B <x> \"q\" 's\n"));
    assert!(!page.contains("<x>"));
}

#[test]
fn each_caveat_sits_under_its_card_title_after_the_step_line() {
    let (_, page) = plain();
    let qc = frame_of(&page, "qc");
    let want = "\n  <div class=\"card\">\n    <h2>Quality Control</h2>\n    <p class=\"sec-meta\">Steps 16b, 16, 33, 02 \
                &middot; <span class=\"badge badge-yellow\">Part not run</span> &middot; file date 2026-11-06 11:30:00 \
                &middot; <a class=\"anchor\" href=\"#qc\">link to this section</a></p>\n    <p class=\"caveat\">Mean depth is";
    assert!(qc.starts_with(want), "{qc}");
    for note in [
        og_render::notes::COVERAGE_NOTE,
        og_render::notes::VARIANTS_NOTE,
    ] {
        let note = og_render::escape::e(note);
        assert_eq!(page.matches(note.as_ref()).count(), 1, "{note}");
        assert_eq!(
            page.matches(&format!("<p class=\"caveat\">{note}</p>"))
                .count(),
            1
        );
    }
    // Red control: a caveat left in a card's body is caught.
    let moved = page.replace(
        &format!("<p class=\"caveat\">{}", og_render::notes::VARIANTS_NOTE),
        &format!("<p class=\"note\">{}", og_render::notes::VARIANTS_NOTE),
    );
    assert_eq!(
        moved
            .matches(&format!(
                "<p class=\"caveat\">{}",
                og_render::notes::VARIANTS_NOTE
            ))
            .count(),
        0
    );
    // The karyotype caveat comes before the coverage one, as in the card.
    let k = summary("", &base_manifest(), "").replace(
        r#""data": {"inferred_sex": "male"}"#,
        r#""data": {"inferred_sex": "male", "karyotype": "47,XXY", "cn_chrX": 2, "cn_chrY": 1}"#,
    );
    let (txt, page) = render_both(&k);
    let qc = frame_of(&page, "qc");
    assert!(qc.find("Named from the copy number") < qc.find("Mean depth is"));
    assert!(qc.contains("<span class=\"label\">Sex chromosomes</span><span class=\"value\"><span class=\"badge badge-red\">47,XXY (indexcov"));
    assert!(
        txt.contains(
            "  Sex chromosomes: 47,XXY (indexcov copy number chrX 2, chrY 1)\n  Named from"
        )
    );
}

#[test]
fn print_rules_keep_the_limits_and_the_header_carries_the_provenance() {
    let (_, page) = plain();
    assert!(page.contains("@media print") && page.contains(".sec:not(:target):not(.always)"));
    assert!(page.contains("class=\"sec always\" id=\"not-assessed\""));
    let header = page.split("</header>").next().unwrap_or("");
    for x in [
        "Sample: <strong>S1</strong>",
        "Pipeline commit:",
        "ClinVar file: 2026-10-01",
        "ClinVar in CPSR's bundle: unknown",
        "HLA database:",
        "Latest run: 2026-11-06 10:00 UTC",
    ] {
        assert!(header.contains(x), "{x}");
    }
}

#[test]
fn a_malformed_section_reads_as_unreadable_and_the_rest_is_written() {
    // Python stopped with a KeyError and wrote no report at all.
    let bad = summary("", &base_manifest(), "").replace(r#""total": 20000, "#, "");
    let (txt, page) = render_both(&bad);
    assert!(
        txt.contains(
            "## Variant calling\n---\n  [could not be read: sections.variants.data has no total]\n"
        ),
        "{txt}"
    );
    assert!(txt.contains("## Quality control"));
    // The HTML card prints None for the missing count, as the Python did.
    assert!(
        frame_of(&page, "variants").contains(
            "<span class=\"label\">Total variants</span><span class=\"value\">None</span>"
        )
    );
    let read = summary("", &base_manifest(), "").replace(
        r#""state": "missing", "step": "02", "source": null, "file_date": null,
                 "note": null, "data": {}"#,
        r#""state": "ok", "step": "02", "source": null, "file_date": null,
                 "note": null, "data": {"input_reads": 5}"#,
    );
    let (txt, page) = render_both(&read);
    assert!(txt.contains("[could not be read: sections.read_count.data has no primary_records]"));
    assert!(
        frame_of(&page, "qc")
            .contains("could not be read: sections.read_count.data has no primary_records")
    );
    // The section that could not be read is failed in the contents and its
    // badge, though each of its keys says ok.
    assert!(
        page.contains("<a href=\"#qc\"><span class=\"dot dot-failed\" title=\"Failed\">"),
        "{page}"
    );
    assert!(frame_of(&page, "qc").contains("<span class=\"badge badge-red\">Failed</span>"));
    assert!(page.contains("<a href=\"#variants\"><span class=\"dot dot-ok\""));
    assert!(frame_of(&page, "variants").contains("Total variants"));
    // The four quality keys and the variants gone (Python stopped with a
    // KeyError): the text leaves both blocks and the chapter out, the
    // quality card still says what could not be measured.
    let mut v = pycompat::json::loads(&summary("", &base_manifest(), "")).expect("json");
    if let pycompat::json::Value::Object(top) = &mut v
        && let Some(pycompat::json::Value::Object(secs)) = top.get_mut("sections")
    {
        for k in [
            "coverage",
            "sex_check",
            "sample_qc",
            "read_count",
            "variants",
        ] {
            secs.shift_remove(k);
        }
    }
    let gone = pycompat::json::dumps(&v, &pycompat::json::Dump::default()).expect("dumped");
    let (txt, page) = render_both(&gone);
    assert!(!txt.contains("## ") || !txt.contains("QUALITY"), "{txt}");
    assert!(frame_of(&page, "qc").contains("Mean depth</span><span class=\"value\">not available"));
    assert!(frame_of(&page, "variants").contains("Not run"));
}

fn tmp(tag: &str) -> std::path::PathBuf {
    let d = std::env::temp_dir().join(format!("og-render-{tag}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).expect("a temp folder");
    d
}

fn cli(args: &[&str]) -> (u8, String, String) {
    let args: Vec<String> = args.iter().map(|a| (*a).to_owned()).collect();
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = render::cli(&args, &mut out, &mut err);
    (
        code,
        String::from_utf8(out).expect("utf-8"),
        String::from_utf8(err).expect("utf-8"),
    )
}

#[test]
fn the_verb_writes_both_reports_and_refuses_what_it_does_not_take() {
    let d = tmp("cli");
    let sum = d.join("summary.json");
    std::fs::write(&sum, summary("", &base_manifest(), "")).expect("written");
    let (txt, page) = (d.join("r.txt"), d.join("r.html"));
    let (code, out, err) = cli(&[
        "--summary",
        sum.to_str().unwrap_or_default(),
        "-o",
        txt.to_str().unwrap_or_default(),
        &format!("--out={}", page.display()),
    ]);
    assert_eq!(
        (code, out.as_str(), err.as_str()),
        (0, "Sections: missing 1, ok 4\n", "")
    );
    let (want_txt, want_page) = plain();
    assert_eq!(std::fs::read_to_string(&txt).ok(), Some(want_txt));
    assert_eq!(std::fs::read_to_string(&page).ok(), Some(want_page));
    let names: Vec<String> = std::fs::read_dir(&d)
        .expect("listed")
        .filter_map(|e| e.ok().map(|e| e.file_name().to_string_lossy().into_owned()))
        .collect();
    assert_eq!(names.len(), 3, "no temporary file is left: {names:?}");
    // Refused before anything is read: an unknown option, an option of the
    // parts still in Python, a report name of another kind, no summary.
    for (args, why) in [
        (
            vec!["--og-scaffold-probe"],
            "unrecognized arguments: --og-scaffold-probe",
        ),
        (
            vec!["--summary", "x.json", "--search-index", "i.json"],
            "--search-index is not ported yet",
        ),
        (
            vec!["--summary", "x.json", "-o", "r.pdf"],
            "-o r.pdf: the name must end in .txt or .html",
        ),
        (vec!["-o", "r.txt"], "give --summary"),
    ] {
        let (code, out, err) = cli(&args);
        assert_eq!((code, out.as_str()), (2, ""), "{args:?}");
        assert!(
            err.contains(why) && err.starts_with("usage: og report render"),
            "{err}"
        );
    }
    // A summary that cannot be read exits 1 and writes nothing.
    let gone = d.join("gone.txt");
    let (code, _, err) = cli(&[
        "--summary",
        "no-such.json",
        "-o",
        gone.to_str().unwrap_or_default(),
    ]);
    assert_eq!(code, 1, "{err}");
    assert!(!gone.exists());
    let _ = std::fs::remove_dir_all(&d);
}
