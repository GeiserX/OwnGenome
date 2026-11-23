use og_render::model::Summary;
use og_render::{html, text};

fn reports(fields: &str, sections: &str) -> [String; 2] {
    let summary = Summary::parse(&format!(
        r#"{{"sample": "S", "generated_utc": "2026-01-02 03:04 UTC",
            "run": {{"started_utc": null, "declared_sex": "male"}},
            "databases": {{"clinvar_release": null, "cpsr_clinvar_release": null, "hla_database": null}},
            "cyp2d6": null, "manifest": [], "not_run": [], "not_assessed": [],
            "sections": {{{sections}}}{fields}}}"#
    ))
    .unwrap();
    [
        text::report(&summary).unwrap(),
        html::report(&summary).unwrap(),
    ]
}

#[test]
fn identity_uncertainty_qualifies_a_possible_match() {
    let [text, html] = reports(
        "",
        r#""sample_qc": {"state": "ok", "data": {
            "identity": "uncertain", "identity_reason": "only 300 shared sites",
            "same_person_as": "S2 (0.99)"}}"#,
    );
    assert!(text.contains("Sample identity: uncertain (only 300 shared sites)"));
    assert!(html.contains(r#"class="badge badge-yellow">uncertain (only 300 shared sites)<"#));
    assert!(html.contains(r#"class="badge badge-yellow">S2 (0.99)<"#));
    for report in [text, html] {
        assert!(report.contains("Possible match to"));
        assert!(report.contains("S2 (0.99)"));
        assert!(!report.contains("The same person as"));
    }
    for report in reports(
        "",
        r#""sample_qc": {"state": "ok", "data": {
            "identity": "checked", "same_person_as": "S2 (0.99)"}}"#,
    ) {
        assert!(!report.contains("Sample identity"));
        assert!(report.contains("The same person as"));
    }
}

#[test]
fn contamination_uncertainty_qualifies_the_estimate() {
    let [text, html] = reports(
        "",
        r#""sample_qc": {"state": "ok", "data": {
          "contamination": "ok", "freemix": "0.0020", "freemix_warn_above": "0.03",
          "contamination_confidence": "uncertain",
          "contamination_confidence_reason": "FREEMIX rests on 300 markers & needs more"}}"#,
    );
    assert!(text.contains("uncertain (FREEMIX 0.0020): FREEMIX rests on 300 markers & needs more"));
    assert!(html.contains(
        r#"class="badge badge-yellow">uncertain (FREEMIX 0.0020): FREEMIX rests on 300 markers &amp; needs more<"#
    ));
    for report in [text, html] {
        assert!(!report.contains("warning above 0.03"));
    }
    let [text, html] = reports(
        "",
        r#""sample_qc": {"state": "ok", "data": {
          "contamination": "ok", "freemix": "0.0020", "freemix_warn_above": "0.03"}}"#,
    );
    for report in [text, html] {
        assert!(report.contains("0.0020 (warning above 0.03)"));
        assert!(!report.contains("uncertain (FREEMIX"));
    }
}

#[test]
fn contamination_warning_precedes_uncertainty() {
    let [text, html] = reports(
        "",
        r#""sample_qc": {"state": "ok", "data": {
          "contamination": "warn", "freemix": "0.0930", "freemix_warn_above": "0.03",
          "contamination_confidence": "uncertain",
          "contamination_confidence_reason": "FREEMIX rests on 300 markers"}}"#,
    );
    let warning = "0.0930 (above 0.03: possible contamination, see docs/33-sample-qc.md). The estimate is uncertain: FREEMIX rests on 300 markers";
    assert!(text.contains(warning));
    assert!(html.contains(&format!(r#"class="badge badge-yellow">{warning}<"#)));
}

#[test]
fn selected_steps_show_the_reason_the_summary_records() {
    let sections = r#"
        "carrier": {"title": "Carrier screen <A & B>", "state": "missing",
          "note": "selected, no result: --tools named carrier & step 41a wrote nothing"},
        "roh": {"title": "Runs of homozygosity", "state": "failed",
          "note": "step 11 failed"}"#;
    let [text, html] = reports(r#", "selected_no_result": ["carrier"]"#, sections);
    let reason = "Carrier screen <A & B>: selected, no result: --tools named carrier & step 41a wrote nothing";
    let escaped = "Carrier screen &lt;A &amp; B&gt;: selected, no result: --tools named carrier &amp; step 41a wrote nothing";
    assert!(text.contains(reason));
    assert!(html.contains(escaped));
    let [text, html] = reports("", sections);
    assert!(!text.contains(reason));
    assert!(!html.contains(escaped));
}

#[test]
fn coverage_marks_use_the_scope_of_each_section() {
    let sections = r#"
        "callability": {"title": "Callability", "state": "ok",
          "partial": "partial: covers 2 of 22 autosomes"},
        "prs": {"title": "Polygenic risk scores", "state": "stale",
          "partial": "partial: not known, from an earlier run"},
        "ancestry_fine": {"title": "Fine-grained ancestry", "state": "ok",
          "partial": "partial: limited to the run's intervals (chr20:1-9)"},
        "telomere": {"title": "Telomere content", "state": "ok",
          "partial": "partial: covers 2 of 22 autosomes"},
        "roh": {"title": "Runs of homozygosity", "state": "ok",
          "partial": "partial: covers 2 of 22 autosomes"},
        "archaic": {"title": "Archaic ancestry", "state": "ok",
          "partial": "partial: covers 2 of 22 autosomes"},
        "variants": {"title": "Variant calling", "state": "ok"}"#;
    let [text, html] = reports(
        r#", "scope": {"partial": true, "mark": "partial: covers 2 of 22 autosomes"}"#,
        sections,
    );
    for report in [&text, &html] {
        for title in [
            "Callability",
            "Telomere content",
            "Runs of homozygosity",
            "Archaic ancestry",
        ] {
            assert!(report.contains(&format!("{title}: partial: covers 2 of 22 autosomes")));
        }
        assert!(report.contains("Polygenic risk scores: partial: not known, from an earlier run"));
        assert!(!report.contains("Variant calling: partial:"));
    }
    assert!(
        text.contains("Fine-grained ancestry: partial: limited to the run's intervals (chr20:1-9)")
    );
    assert!(html.contains(
        "Fine-grained ancestry: partial: limited to the run&#x27;s intervals (chr20:1-9)"
    ));
}

#[test]
fn partial_scope_is_visible_inside_the_prs_section() {
    let [text, html] = reports(
        "",
        r#""prs": {"title": "Polygenic risk scores", "state": "stale",
          "partial": "partial: not known, from an earlier run",
          "data": {"scores": [], "adjusted": true}}"#,
    );
    let prs_head = "## Polygenic Risk Scores\n---\n  [partial: not known, from an earlier run]\n";
    let prs_card = html.split_once(r#"id="prs">"#).unwrap().1;
    let prs_card = prs_card.split_once("</section>").unwrap().0;
    let marked =
        prs_card.contains(r#"<div class="stale">partial: not known, from an earlier run</div>"#);
    assert_eq!((text.contains(prs_head), marked), (true, true));
}
