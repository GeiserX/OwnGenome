//! What the callability and ClinVar collectors put in summary.json, shown in
//! both reports: the share at one decimal (0.8988 is 89.9%, never the 90% a
//! gene needs), the coding exon that was not read, the owned genes whose own
//! step gave no result, the introns caveat, and the conflicting ClinVar
//! records to review, apart from every count.

use og_render::model::Summary;
use og_render::{html, text};

fn reports(sections: &str) -> [String; 2] {
    let summary = Summary::parse(&format!(
        r#"{{"sample": "S", "generated_utc": "2026-01-02 03:04 UTC",
            "run": {{"started_utc": null, "declared_sex": "male"}},
            "databases": {{"clinvar_release": "2026-09-28", "cpsr_clinvar_release": null, "hla_database": null}},
            "cyp2d6": null, "manifest": [], "not_run": [], "not_assessed": [],
            "sections": {{{sections}}}}}"#
    ))
    .unwrap();
    [
        text::report(&summary).unwrap(),
        html::report(&summary).unwrap(),
    ]
}

/// A callability section: `not_read` and `not_read_owned` as given, and
/// `extra` fields after them.
fn callability(not_read: &str, owned: &str, extra: &str) -> String {
    format!(
        r#""callability": {{"title": "Gene callability (mosdepth)", "state": "ok", "data": {{
            "min_depth": 10, "min_mapq": 20, "min_share": 0.9, "catalog_genes": 120, "read": 110,
            "not_read": [{not_read}], "not_read_owned": [{owned}], "not_on_reference": []{extra}}}}}"#
    )
}

const NEAR: &str = r#"{"gene": "GENEA", "hgnc_id": "HGNC:1", "reasons": ["paralog"],
    "share": 0.8988, "status": "not read"}"#;

#[test]
fn a_share_under_the_bar_never_reads_as_reaching_it() {
    let edge = r#"{"gene": "GENEB", "hgnc_id": "HGNC:2", "reasons": ["paralog"],
        "share": 0.8999, "status": "not read"}"#;
    let low = r#"{"gene": "GENEC", "hgnc_id": "HGNC:3", "reasons": ["paralog"],
        "share": 0.1667, "status": "not read"}"#;
    for report in reports(&callability(&format!("{NEAR}, {edge}, {low}"), "", "")) {
        assert!(report.contains("GENEA (89.9%)"), "{report}");
        assert!(report.contains("GENEB (89.9%)"), "{report}");
        assert!(report.contains("GENEC (16.7%)"), "{report}");
        assert!(!report.contains("(90%)") && !report.contains("(90.0%)"));
    }
}

#[test]
fn a_gene_whose_body_passed_names_the_exon_that_did_not() {
    let exon = r#"{"gene": "GENED", "hgnc_id": "HGNC:4", "reasons": ["paralog"],
        "share": 0.95, "status": "exon not read",
        "worst_exon": {"exon": "GENED:exon3", "share": 0.41}}"#;
    let [text, html] = reports(&callability(exon, "", r#", "exons_checked": 812"#));
    for report in [&text, &html] {
        assert!(
            report.contains("GENED (95.0%; exon GENED:exon3 not read, 41.0%)"),
            "{report}"
        );
        assert!(report.contains("the 812 coding exons checked must pass too"));
    }
    assert!(text.contains("introns included"));
}

#[test]
fn without_an_exon_table_the_caveat_says_no_exon_was_checked() {
    for report in reports(&callability(NEAR, "", "")) {
        assert!(report.contains("introns included"), "{report}");
        assert!(report.contains("no exon was checked"), "{report}");
    }
}

#[test]
fn an_owned_gene_whose_step_did_not_run_says_so() {
    let owned = r#"
        {"gene": "SMN1", "owner": "parascopy", "owner_step": "35", "owner_status": "not run",
         "share": 0.3712, "status": "not read"},
        {"gene": "LPA", "owner": "cn_loci:LPA_KIV2", "owner_step": "46",
         "owner_status": "from an earlier run", "share": 0.5, "status": "not read"},
        {"gene": "CYP2D6", "owner": "pypgx", "owner_step": "32", "owner_status": "ran",
         "share": 0.8, "status": "not read"}"#;
    let [text, html] = reports(&callability(NEAR, owned, ""));
    let gave_none = "Not read here, and the step that calls them gave no result in this run (2): \
        SMN1 (37.1%) by parascopy, owner step 35 not run, \
        LPA (50.0%) by cn_loci:LPA_KIV2, owner step 46 from an earlier run";
    let called = "Not read here, called by a step of their own (1): CYP2D6 (80.0%) by pypgx";
    for report in [&text, &html] {
        assert!(report.contains(gave_none), "{report}");
        assert!(report.contains(called), "{report}");
    }
    assert!(html.contains("Not read, their own step gave no result"));
    assert!(html.contains(r#"class="badge badge-red">2<"#));
    // A summary from before owner_status lists every owned gene as called.
    let old = r#"{"gene": "SMN1", "owner": "parascopy", "share": 0.3712, "status": "not read"}"#;
    let [text, html] = reports(&callability(NEAR, old, ""));
    assert!(text.contains("called by a step of their own (1): SMN1 (37.1%) by parascopy"));
    assert!(!text.contains("gave no result"));
    assert!(!html.contains("their own step gave no result"));
}

/// A ClinVar section with one P/LP hit, and `conflicting` when given.
fn clinvar(conflicting: Option<&str>) -> String {
    let extra = conflicting.map_or(String::new(), |c| format!(r#", "conflicting": {c}"#));
    format!(
        r#""clinvar": {{"title": "ClinVar screen", "state": "ok", "data": {{
            "count": 1, "by_stars": {{"4": 0, "3": 0, "2": 1, "1": 0, "0": 0}},
            "hits": [{{"chrom": "chr2", "pos": 47403191, "ref": "C", "alt": "T", "genotype": "het",
              "gene": "DEMOA", "significance": "Pathogenic",
              "review_status": "criteria provided, multiple submitters, no conflicts", "stars": 2}}]{extra}}}}}"#
    )
}

const CONFLICTING: &str = r#"{"checked": true, "note": null, "count": 1, "hits": [
    {"chrom": "chr5", "pos": 31400123, "ref": "G", "alt": "A", "genotype": "het", "gene": "DEMOD",
     "clinvar_id": "90004", "stars": 1, "review_status": "criteria provided, conflicting classifications",
     "submissions": "Pathogenic(2), Uncertain_significance(1)", "plp_submissions": 2,
     "clinvar_condition": "Demo condition D", "clinvar_condition_ids": ["MONDO:MONDO:0000001"],
     "dp": 9, "gq": 48, "low_confidence": true}]}"#;

const TITLE: &str = "Conflicting classifications with a P/LP submission: review, not counted";

#[test]
fn the_conflicting_records_get_a_review_list_of_their_own() {
    let [text, html] = reports(&clinvar(Some(CONFLICTING)));
    assert!(text.contains(&format!("  {TITLE}: 1")), "{text}");
    assert!(text.contains(
        "    DEMOD chr5:31400123 G>A het: Pathogenic(2), Uncertain_significance(1); 1 star; \
         Demo condition D; low-confidence call (GQ 48, DP 9): check the reads"
    ));
    assert!(html.contains(&format!("{TITLE} (1 of 1)")), "{html}");
    assert!(html.contains("<td>Pathogenic(2), Uncertain_significance(1)</td>"));
    for report in [&text, &html] {
        assert!(report.contains("They are not in the P/LP count"));
    }
}

/// `report` without the lines that hold `marks`: what is left must be the
/// report of the same summary without the list, so the list touches no count.
fn without(report: &str, marks: &[&str]) -> String {
    report
        .lines()
        .filter(|l| !marks.iter().any(|m| l.contains(m)))
        .collect::<Vec<_>>()
        .join("\n")
}

#[test]
fn the_review_list_never_touches_the_counts() {
    let [text, _] = reports(&clinvar(Some(CONFLICTING)));
    let [plain, plain_html] = reports(&clinvar(None));
    assert!(!plain.contains(TITLE) && !plain_html.contains(TITLE));
    assert_eq!(
        without(&text, &[TITLE, "DEMOD", "They are not in the P/LP count"]),
        without(&plain, &[]),
    );
    assert!(text.contains("  Pathogenic/Likely Pathogenic hits: 1"));
    assert!(
        text.contains("By review status: 4 stars 0, 3 stars 0, 2 stars 1, 1 star 0, 0 stars 0")
    );
    let [_, html] = reports(&clinvar(Some(CONFLICTING)));
    assert!(html.contains(
        r#"<span class="label">ClinVar matches</span><span class="value"><span class="badge badge-yellow">1<"#
    ));
    assert!(html.contains("ClinVar Hits (best-reviewed first, 1 of 1)"));
}

#[test]
fn a_list_that_was_not_checked_never_reads_as_zero() {
    let skipped = r#"{"checked": false, "note": "not checked: run setup.sh --refresh clinvar",
        "count": null, "hits": []}"#;
    let [text, html] = reports(&clinvar(Some(skipped)));
    let line = format!("{TITLE}: not checked: run setup.sh --refresh clinvar");
    assert!(text.contains(&line), "{text}");
    assert!(
        html.contains(r#"class="badge badge-yellow">not checked: run setup.sh --refresh clinvar<"#)
    );
    assert!(!text.contains(&format!("{TITLE}: 0")));
    let none = r#"{"checked": true, "note": null, "count": 0, "hits": []}"#;
    let [text, html] = reports(&clinvar(Some(none)));
    assert!(text.contains(&format!("{TITLE}: 0")));
    assert!(html.contains(r#"class="badge badge-gray">0<"#));
}
