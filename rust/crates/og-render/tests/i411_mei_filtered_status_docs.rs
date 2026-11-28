//! The MEI section explains every status step 50 writes (#411 sv-34).
//!
//! `og mei report` reads a catalog site as seen, not seen, not callable, or
//! `call filtered out (FLAGS)` when the only calls MEGAnE made there are ones
//! it flagged. The note under the section explained the first three only, so
//! a reader saw "call filtered out (LC)" with nothing to say what it means.
//! The note now says it: MEGAnE called an insertion there but flagged it, so
//! it is not a finding, and not proof of absence either.
//!
//! The report-render parity goldens (tests/parity/cases.tsv) hold the same
//! note, byte for byte.

use og_render::model::Summary;
use og_render::notes::MEI_NOTE;
use og_render::sections::mei;
use og_render::{html, text};

/// The sentence that explains the filtered-out status.
const FILTERED: &str = "Call filtered out: MEGAnE called an insertion there but flagged it (LC is low \
    confidence), so it is not a finding, and not proof of absence either.";

/// A summary whose MEI section ran, with one site of each status.
fn summary() -> String {
    r#"{
 "schema_version": 2, "sample": "S1", "generated_utc": "2026-11-20 12:00 UTC",
 "run": {"started_utc": "2026-11-20 10:00 UTC", "declared_sex": "male"},
 "databases": {"clinvar_release": "2026-10-01", "cpsr_clinvar_release": null, "hla_database": null},
 "manifest": [],
 "sections": {
  "coverage": {"title": "Coverage (mosdepth)", "state": "ok", "step": "16b", "source": "qc/S1.mosdepth.summary.txt",
               "file_date": "2026-11-20 11:00:00", "note": null, "data": {"mean_depth": 31.25}},
  "mei": {"title": "Mobile element insertions (MEGAnE)", "state": "ok", "step": "50",
          "source": "mei/S1_mei_sites.tsv", "file_date": "2026-11-20 11:00:00", "note": null,
          "data": {"calls": 2, "by_element": {"Alu": 1, "SVA": 1}, "sites": [
    {"gene": "FKTN", "element": "SVA", "condition": "Fukuyama congenital muscular dystrophy",
     "inheritance": "autosomal recessive", "status": "call filtered out (LC)", "mean_depth": "30.0",
     "call": "chr9:105639660 SVA_F 1 copy LC"},
    {"gene": "TAF1", "element": "SVA", "condition": "X-linked dystonia-parkinsonism", "inheritance": "X-linked",
     "status": "not callable", "mean_depth": "0.0", "call": "."},
    {"gene": "BBS1", "element": "SVA", "condition": "Bardet-Biedl syndrome", "inheritance": "autosomal recessive",
     "status": "not seen", "mean_depth": "30.0", "call": "."}]}}
 },
 "cyp2d6": null, "not_run": [], "not_assessed": []
}"#
    .to_owned()
}

/// The note names each status a catalog row can have, and says what a
/// filtered-out call means.
#[test]
fn the_note_explains_every_status() {
    assert!(MEI_NOTE.contains(FILTERED), "{MEI_NOTE}");
    for status in [
        "Not callable:",
        "Not seen is not proof of absence",
        "Call filtered out:",
    ] {
        assert!(
            MEI_NOTE.contains(status),
            "{status} is not explained: {MEI_NOTE}"
        );
    }
}

/// Both reports show the filtered-out row and the note that explains it.
#[test]
fn both_reports_explain_a_filtered_out_call() {
    let s = Summary::parse(&summary()).expect("a summary");
    let lines = mei::lines(&s.section("mei")).unwrap();
    assert_eq!(
        lines[1],
        "FKTN SVA (Fukuyama congenital muscular dystrophy; autosomal recessive): call filtered out (LC) \
         (depth 30.0, chr9:105639660 SVA_F 1 copy LC)"
    );
    assert_eq!(lines.last().unwrap(), MEI_NOTE);

    let txt = text::report(&s).expect("the text report");
    assert!(txt.contains("call filtered out (LC)"), "{txt}");
    assert!(
        txt.contains(FILTERED),
        "the text report lacks the sentence:\n{txt}"
    );

    let page = html::report(&s).expect("the HTML report");
    let card = page.split_once("id=\"mei\">").map_or("", |(_, rest)| {
        rest.split("</section>").next().unwrap_or("")
    });
    assert!(card.contains("<td>call filtered out (LC)</td>"), "{card}");
    assert!(
        card.contains(FILTERED),
        "the HTML card lacks the sentence:\n{card}"
    );
}
