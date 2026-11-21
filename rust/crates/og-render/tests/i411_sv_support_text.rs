//! The reports give the right reason when an SV's caller support is not
//! known (#411 sv-01).
//!
//! `og summary collect` writes an SV's `support` as `null` in two cases:
//! step 22's consensus did not run, or it ran and cannot tell (an `SV_VCF`
//! run's `sv_vcf` rows, a row with no record ID, an ID in two events, a
//! vote with no readable ID). The reports gave both the reason "no step 22
//! consensus", so a user whose step 22 did run read that it did not. The
//! AnnotSV header counted only single-caller rows, so with the consensus
//! set and every row unknown it read "single caller: 0".
//!
//! Now the AnnotSV section's `consensus` field picks the reason, in the
//! AnnotSV block and in the secondary-findings rows, and the header counts
//! the unknown rows when the consensus ran. Without the consensus the old
//! text stays.

use og_render::escape::e;
use og_render::model::Summary;
use og_render::sections::{acmg_sf, annotsv};
use og_render::{html, text};

const CANNOT_TELL: &str = "not known (step 22's consensus cannot tell)";
const NO_CONSENSUS: &str = "not known (no step 22 consensus)";

/// A summary with an AnnotSV section of two class 4-5 SVs (one seen by one
/// caller, one whose support is `second`) and an SF tier with one
/// structural row whose support is null. `consensus` is the AnnotSV
/// section's `consensus` field, as JSON.
fn summary(consensus: &str, second: &str) -> Summary {
    Summary::parse(&summary_json(consensus, second, true)).expect("a summary")
}

/// The summary's JSON; without `with_annotsv` it has no AnnotSV section.
fn summary_json(consensus: &str, second: &str, with_annotsv: bool) -> String {
    let annotsv = format!(
        r#""annotsv": {{"title": "Structural Variants, Class 4 and 5 (AnnotSV)", "state": "ok", "step": "05",
   "source": "annotsv/S1_sv_class4_5_checked.tsv", "file_date": "2026-11-21 11:00:00", "note": null,
   "data": {{"count": 2, "consensus": {consensus}, "consensus_note": null, "single_caller": 1,
    "exon_edge": 0, "exon_edge_checked": true, "common": 0, "downgraded": [], "rows": [
    {{"type": "DEL", "sv": "chr1:1000-5000", "acmg_class": "5", "genes": ["G1"], "caller": "manta",
     "support": 1, "support_callers": ["manta"], "exon_edge": "no", "population_af": null}},
    {{"type": "DUP", "sv": "chr2:2000-9000", "acmg_class": "4", "genes": ["G2"], "caller": "sv_vcf",
     "support": {second}, "support_callers": ["manta", "delly"], "exon_edge": "no", "population_af": null}}]}}}},"#
    );
    let annotsv = if with_annotsv { annotsv.as_str() } else { "" };
    format!(
        r#"{{
 "schema_version": 2, "sample": "S1", "generated_utc": "2026-11-21 12:00 UTC",
 "run": {{"started_utc": "2026-11-21 10:00 UTC", "declared_sex": "female"}},
 "databases": {{"clinvar_release": "2026-10-01", "cpsr_clinvar_release": null, "hla_database": null}},
 "manifest": [],
 "sections": {{
  {annotsv}
  "clinical": {{"title": "Clinical Variant Filter", "state": "ok", "step": "23",
   "source": "clinical/S1_clinical.tsv", "file_date": "2026-11-21 11:00:00", "note": null,
   "data": {{"variants": 3, "acmg_sf": {{"version": "SF v3.3", "genes_on_list": 84, "clinvar_hits": [],
    "high_impact": [], "structural_checked": true, "structural": [
    {{"gene": "G1", "variant": "DEL chr1:1000-5000", "acmg_class": "5", "caller": "sv_vcf", "support": null}}]}}}}}}
 }},
 "cyp2d6": null, "not_run": [], "not_assessed": []
}}"#
    )
}

/// The block of section `sid` in the HTML report.
fn frame<'a>(page: &'a str, sid: &str) -> &'a str {
    page.split_once(&format!("id=\"{sid}\">"))
        .map_or("", |(_, rest)| {
            rest.split("</section>").next().unwrap_or("")
        })
}

/// Both reports, text and HTML.
fn reports(s: &Summary) -> (String, String) {
    (
        text::report(s).expect("the text report"),
        html::report(s).expect("the HTML report"),
    )
}

/// The consensus ran: a null support is one it cannot tell, in the AnnotSV
/// rows, the header and the SF row; "no step 22 consensus" is gone.
#[test]
fn a_consensus_that_cannot_tell_says_so() {
    let s = summary(r#""sv_merged/S1_sv_consensus.vcf.gz""#, "null");
    let (txt, page) = reports(&s);

    let head = "  Class 4-5 SVs: 2 (single caller: 1, not known: 1; deletions on exon edges: 0; common: 0)";
    assert!(txt.lines().any(|l| l == head), "{txt}");
    assert!(
        txt.lines().any(
            |l| l.starts_with("    DUP  chr2:2000-9000  class 4  G2  support: ")
                && l.ends_with(CANNOT_TELL)
        ),
        "{txt}"
    );
    assert!(
        txt.contains("G1       DEL chr1:1000-5000  .  AnnotSV class 5, caller support not known (step 22's consensus cannot tell)"),
        "{txt}"
    );
    // The consensus ran, so its "did not run" note is not there either.
    assert!(!txt.contains("did not run"), "{txt}");

    let sv = frame(&page, "structural-variants");
    assert!(sv.contains(&format!("<td>{}</td>", e(CANNOT_TELL))), "{sv}");
    assert!(
        sv.contains(&*e(
            "single caller: 1, not known: 1; deletions on exon edges: 0"
        )),
        "{sv}"
    );
    let clinical = frame(&page, "clinical-filter");
    assert!(
        clinical.contains(&*e(
            "AnnotSV class 5, caller support not known (step 22's consensus cannot tell)"
        )),
        "{clinical}"
    );
    for (name, report) in [("text", txt.as_str()), ("HTML", page.as_str())] {
        assert!(
            !report.contains("no step 22 consensus"),
            "the {name} report still says the consensus did not run:\n{report}"
        );
    }
}

/// The consensus ran and saw every row: the header has no unknown count.
#[test]
fn a_consensus_that_saw_every_row_adds_no_unknown_count() {
    let s = summary(r#""sv_merged/S1_sv_consensus.vcf.gz""#, "2");
    let lines =
        annotsv::lines(&og_render::sections::clinvar::D::of(&s.section("annotsv"))).unwrap();
    assert_eq!(
        lines[0],
        "Class 4-5 SVs: 2 (single caller: 1; deletions on exon edges: 0; common: 0)"
    );
    let (txt, _) = reports(&s);
    assert!(txt.contains("support: 2 callers (manta, delly)"), "{txt}");
}

/// No consensus: the old reason stays everywhere, with the note that says
/// step 22 did not run, and the header reads "not known" as before.
#[test]
fn without_the_consensus_the_old_reason_stays() {
    let s = summary("null", "null");
    let (txt, page) = reports(&s);

    let head =
        "  Class 4-5 SVs: 2 (single caller: not known; deletions on exon edges: 0; common: 0)";
    assert!(txt.lines().any(|l| l == head), "{txt}");
    assert!(
        txt.lines().any(
            |l| l.starts_with("    DUP  chr2:2000-9000  class 4  G2  support: ")
                && l.ends_with(NO_CONSENSUS)
        ),
        "{txt}"
    );
    assert!(
        txt.contains("G1       DEL chr1:1000-5000  .  AnnotSV class 5, caller support not known (no step 22 consensus)"),
        "{txt}"
    );
    assert!(txt.contains("step 22's consensus did not run"), "{txt}");
    let sv = frame(&page, "structural-variants");
    assert!(
        sv.contains(&format!("<td>{}</td>", e(NO_CONSENSUS))),
        "{sv}"
    );
    for (name, report) in [("text", txt.as_str()), ("HTML", page.as_str())] {
        assert!(
            !report.contains("cannot tell"),
            "the {name} report gives the consensus's reason without a consensus:\n{report}"
        );
    }
}

/// The SF rows read the AnnotSV section's flag: with no AnnotSV section at
/// all, the old reason stays.
#[test]
fn sf_rows_without_an_annotsv_section_keep_the_old_reason() {
    let s = Summary::parse(&summary_json(
        r#""sv_merged/S1_sv_consensus.vcf.gz""#,
        "null",
        false,
    ))
    .expect("a summary");
    let txt = acmg_sf::text(&s).unwrap().join("\n");
    assert!(
        txt.contains("AnnotSV class 5, caller support not known (no step 22 consensus)"),
        "{txt}"
    );
    assert!(!txt.contains("cannot tell"), "{txt}");
}
