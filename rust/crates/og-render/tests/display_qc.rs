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
    assert!(html.contains(">uncertain (only 300 shared sites)<"));
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
fn selected_steps_show_the_reason_the_summary_records() {
    let sections = r#"
        "carrier": {"title": "Carrier screen", "state": "missing",
          "note": "selected, no result: --tools named carrier, but step 41a wrote nothing"},
        "roh": {"title": "Runs of homozygosity", "state": "failed",
          "note": "step 11 failed"}"#;
    let reason =
        "Carrier screen: selected, no result: --tools named carrier, but step 41a wrote nothing";
    for report in reports(r#", "selected_no_result": ["carrier"]"#, sections) {
        assert!(report.contains(reason));
    }
    for report in reports("", sections) {
        assert!(!report.contains(reason));
    }
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
    for report in [text, html.replace("&#x27;", "'")] {
        for title in [
            "Callability",
            "Telomere content",
            "Runs of homozygosity",
            "Archaic ancestry",
        ] {
            assert!(report.contains(&format!("{title}: partial: covers 2 of 22 autosomes")));
        }
        assert!(report.contains("Polygenic risk scores: partial: not known, from an earlier run"));
        assert!(report.contains(
            "Fine-grained ancestry: partial: limited to the run's intervals (chr20:1-9)"
        ));
        assert!(!report.contains("Variant calling: partial:"));
    }
}
