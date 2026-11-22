//! The PGx warnings and evidence shown in the text and HTML reports.

use og_render::model::Summary;
use og_render::sections::{hla, pgx};
use og_render::{html, text};
use pycompat::json::loads;

fn summary(sections: &str, top: &str) -> Summary {
    Summary::parse(&format!(
        r#"{{"schema_version": 2, "sample": "TEST", "generated_utc": "2026-11-22 12:00 UTC",
 "run": {{"started_utc": null, "declared_sex": null}},
 "databases": {{"clinvar_release": null, "cpsr_clinvar_release": null, "hla_database": null}},
 "manifest": [], "sections": {{{sections}}}, "cyp2d6": {top}, "not_run": [], "not_assessed": []}}"#
    ))
    .expect("summary")
}

#[test]
fn html_shows_a_cpic_warning_even_when_every_call_is_normal() {
    let s = summary(
        r#""cpic": {"state": "ok", "source": "cpic/TEST_cpic_summary.tsv", "data": {
        "parse_failed": false, "non_normal": 0, "ambiguous": 0, "not_called": 0,
        "genes": [{"gene": "CYP2C19", "diplotype": "*1/*1", "phenotype": "Normal Metabolizer", "status": "normal"}],
        "warnings": ["Synthetic drug contraindication warning"]}}"#,
        "null",
    );
    let txt = text::report(&s).expect("text report");
    let page = html::report(&s).expect("HTML report");
    let warning = "WARNING: Synthetic drug contraindication warning";
    assert!(txt.contains(warning));
    assert!(page.contains(warning));
    assert!(!page.contains("<th>Diplotype</th>"));
}

#[test]
fn reports_count_and_name_a_cpic_call_with_unread_positions() {
    let s = summary(
        r#""cpic": {"state": "ok", "source": "cpic/TEST_cpic_summary.tsv", "data": {
        "parse_failed": false, "non_normal": 0, "ambiguous": 0, "not_called": 0, "incomplete": 1,
        "genes": [{"gene": "CYP2C19", "diplotype": "*1/*1",
            "phenotype": "Normal Metabolizer [incomplete: positions not read: rs4244285]",
            "status": "incomplete"}], "warnings": []}}"#,
        "null",
    );
    let txt = text::report(&s).expect("text report");
    let page = html::report(&s).expect("HTML report");
    let label = "Genes called with positions PharmCAT never read (no drug guidance)";
    assert!(txt.contains(&format!("{label}: 1")));
    assert!(page.contains(&format!("{label}</span><span class=\"value\">1</span>")));
    let reason = "incomplete: positions not read: rs4244285";
    assert!(
        txt.lines()
            .any(|line| line.contains("CYP2C19") && line.contains(reason))
    );
    assert!(page.contains(&format!(
        "<td>CYP2C19</td><td>*1/*1</td><td>Normal Metabolizer [{reason}]</td>"
    )));
}

fn depth_line(shares: &str) -> String {
    let s = summary(
        &format!(
            r#""pypgx": {{"state": "ok", "data": {{"depth_check": {{
            "status": "ok", "gene_depth_all": "18.50", "flank_depth_all": "29.50",
            "gene_depth_mapq1": "16.50", "flank_depth_mapq1": "29.40"{shares}
            }}}}}}"#
        ),
        "null",
    );
    let lines = pgx::pypgx_lines(&s.section("pypgx")).expect("pypgx lines");
    assert_eq!(lines.len(), 1);
    lines.into_iter().next().expect("depth check line")
}

#[test]
fn depth_check_shows_mean_depths_and_reports_unmeasured_shares_once() {
    let line = depth_line("");
    assert!(
        line.contains("mean depth, all reads: CYP2D6 18.50, flanks 29.50"),
        "{line}"
    );
    assert!(
        line.contains("MAPQ >= 1: CYP2D6 16.50, flanks 29.40"),
        "{line}"
    );
    assert!(
        line.contains("MAPQ shares: not measured in this run"),
        "{line}"
    );
    assert_eq!(line.matches("not measured").count(), 1, "{line}");
}

#[test]
fn measured_depth_check_shows_mean_depths_and_all_four_shares() {
    let line = depth_line(
        r#", "gene_mapq1_share": "0.89", "gene_mapq20_share": "0.75",
            "body_mapq1_share": "0.98", "body_mapq20_share": "0.92""#,
    );
    assert!(line.contains("CYP2D6 18.50, flanks 29.50"), "{line}");
    assert!(
        line.contains("Cyrius's region with REP6 0.89 / 0.75, gene body without REP6 0.98 / 0.92"),
        "{line}"
    );
}

#[test]
fn cyp2d6_consensus_shows_weak_pypgx_evidence_even_when_callers_agree() {
    let s = summary(
        r#""pypgx": {"state": "ok", "data": {
        "genes_called": 1, "genes_total": 1, "cyp2d6": "*1/*4",
        "cyp2d6_evidence": {"defining_filter": "Synthetic filter", "quality": "low-confidence: check the reads"}}}"#,
        r#"{"calls": {"pypgx": "*1/*4", "Cyrius": "*4/*1"}, "agree": true,
        "consensus": {"result": "*1/*4", "passed_to_pharmcat": true,
        "reason": "pypgx and Cyrius agree; the depth check passed"}}"#,
    );
    let txt = text::report(&s).expect("text report");
    let page = html::report(&s).expect("HTML report");
    let flag = "Weak pypgx CYP2D6 evidence: Synthetic filter, low-confidence: check the reads";
    assert!(txt.contains(flag), "{txt}");
    assert!(page.contains(flag));
}

#[test]
fn a_tied_hla_call_with_good_quality_can_still_be_ambiguous_for_pharmcat() {
    let locus = loads(
        r#"{"gene": "HLA-B", "alleles": ["B*57:01,B*57:02", "B*07:02"],
        "quality": ["1", "1"], "candidates": [2, 1],
        "low_confidence": false, "withheld_from_pharmcat": true}"#,
    )
    .expect("locus");
    let line = hla::hla_text(&locus).expect("HLA report line");
    assert!(line.contains("ambiguous for PharmCAT"), "{line}");
    assert!(line.contains("not passed to PharmCAT"), "{line}");
}
