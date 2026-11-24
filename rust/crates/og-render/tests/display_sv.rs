//! The reports show what the structural-variant, repeat, insertion,
//! Locityper and LPA steps now say (#411: sv-23, sv-24, sv-29, sv-33,
//! sv-34, sv-36, sv-37, sv-39, D14).
//!
//! Each of those steps learned to say something new in its own table, and
//! summary.json carries it, but the reports still printed the old line:
//!
//! - AnnotSV (step 05) writes "no frequency data" when it gave no
//!   benign-set AF, and per-record evidence (`benign_match`,
//!   `clinical_common`). The reports printed nothing for a missing AF, so it
//!   read as rare, and showed the aggregate AF as if it were the SV's own.
//! - Locityper (step 53) reads "uncalibrated" outside HLA and KIR, with
//!   what the rule alone would say; the reports hid that reason.
//! - The HTT range of reduced penetrance (36-39) reached the summary but no
//!   line: a reader saw Stranger's plain "full_mutation".
//! - LPA (step 46) says its 95% interval covers read sampling only; the
//!   headline still said a bare "95% interval".
//! - MEI (step 50) "call filtered out (LC)" and "not callable" reach both
//!   reports with the call and the note that explains them.
//!
//! A summary without the new keys keeps its old lines (the control of each
//! test), so the report-render parity goldens do not move.

use og_render::escape::e;
use og_render::model::Summary;
use og_render::{html, text};

/// The head of a summary, then `sections`, then the tail.
fn summary(sections: &str) -> Summary {
    Summary::parse(&format!(
        r#"{{
 "schema_version": 2, "sample": "S1", "generated_utc": "2026-11-24 12:00 UTC",
 "run": {{"started_utc": "2026-11-24 10:00 UTC", "declared_sex": "female"}},
 "databases": {{"clinvar_release": "2026-10-01", "cpsr_clinvar_release": null, "hla_database": null}},
 "manifest": [],
 "sections": {{ {sections} }},
 "cyp2d6": null, "not_run": [], "not_assessed": []
}}"#
    ))
    .expect("a summary")
}

/// A section of state ok from `step`, with `data`.
fn section(key: &str, title: &str, step: &str, data: &str) -> String {
    format!(
        r#""{key}": {{"title": "{title}", "state": "ok", "step": "{step}", "source": "{key}/S1.tsv",
           "file_date": "2026-11-24 11:00:00", "note": null, "data": {data}}}"#
    )
}

/// What the reports print for an SV whose benign-set AF AnnotSV left blank.
const NO_FREQUENCY: &str = "benign-set AF: no frequency data (unknown, not the same as rare)";

/// The start of the caveat that tells the aggregate AF from the per-record check.
const EVIDENCE_NOTE: &str = "The benign-set AF is AnnotSV's aggregate: the highest AF of every \
    benign record of the call's type it lists, larger SVs around the call included";

/// Both reports, text and HTML.
fn reports(s: &Summary) -> (String, String) {
    (
        text::report(s).expect("the text report"),
        html::report(s).expect("the HTML report"),
    )
}

/// The block of section `sid` in the HTML report.
fn frame<'a>(page: &'a str, sid: &str) -> &'a str {
    page.split_once(&format!("id=\"{sid}\">"))
        .map_or("", |(_, rest)| {
            rest.split("</section>").next().unwrap_or("")
        })
}

// --- AnnotSV: unknown frequency, aggregate AF and per-record evidence ---------------

/// An AnnotSV section with `rows` (JSON objects) as its class 4-5 rows.
fn annotsv(rows: &[&str]) -> Summary {
    summary(&section(
        "annotsv",
        "Structural Variants, Class 4 and 5 (AnnotSV)",
        "05",
        &format!(
            r#"{{"count": {}, "consensus": null, "consensus_note": null, "single_caller": {},
               "exon_edge": 0, "exon_edge_checked": true, "common": 1, "downgraded": [], "rows": [{}]}}"#,
            rows.len(),
            rows.len(),
            rows.join(", ")
        ),
    ))
}

/// A class 4 deletion over G1 from Manta, with `extra` fields.
fn sv_row(sv: &str, extra: &str) -> String {
    format!(
        r#"{{"type": "DEL", "sv": "{sv}", "acmg_class": "4", "genes": ["G1"], "caller": "manta",
            "support": 1, "support_callers": ["manta"], "exon_edge": "no", {extra}}}"#
    )
}

/// The rows of step 05 since #601: one with no frequency data, one whose
/// only benign record is a larger common SV around it, one with a same-size
/// common record.
fn evidence_rows() -> [String; 3] {
    [
        sv_row(
            "chr1:1000-5000",
            r#""population_af": null, "population_source": "no frequency data", "common": false,
               "benign_match": "none", "clinical_common": "no""#,
        ),
        sv_row(
            "chr2:2000-2500",
            r#""population_af": 0.05, "population_source": "gnomAD-SV_v4.1", "common": true,
               "benign_match": "contained", "clinical_common": "no""#,
        ),
        sv_row(
            "chr3:1000-9000",
            r#""population_af": 0.05, "population_source": "gnomAD-SV_v4.1", "common": true,
               "benign_match": "same", "clinical_common": "yes""#,
        ),
    ]
}

#[test]
fn an_sv_without_frequency_data_says_so_instead_of_reading_as_rare() {
    let [blank, ..] = evidence_rows();
    let s = annotsv(&[&blank]);
    let (txt, page) = reports(&s);
    assert!(txt.contains(&format!("         {NO_FREQUENCY}\n")), "{txt}");
    let card = frame(&page, "structural-variants");
    assert!(card.contains(NO_FREQUENCY), "{card}");
    assert!(card.contains(&*e(EVIDENCE_NOTE)), "{card}");

    // Control: a row from before step 05 wrote the marker (a blank source)
    // keeps its old lines, with no frequency line and no new note.
    let old = sv_row(
        "chr1:1000-5000",
        r#""population_af": null, "population_source": """#,
    );
    let (txt, page) = reports(&annotsv(&[&old]));
    assert!(!txt.contains("no frequency data"), "{txt}");
    assert!(!txt.contains(EVIDENCE_NOTE), "{txt}");
    assert!(!page.contains("no frequency data"), "{page}");
}

#[test]
fn the_aggregate_af_is_told_apart_from_a_same_size_record() {
    let rows = evidence_rows();
    let s = annotsv(&rows.iter().map(String::as_str).collect::<Vec<_>>());
    let (txt, page) = reports(&s);
    for want in [
        "    DEL  chr2:2000-2500  class 4  G1  support: single caller (manta)\n         \
         benign-set AF 0.05 (gnomAD-SV_v4.1; the highest AF of the benign records listed, larger \
         SVs around it included), common\n         \
         benign records: only larger SVs around it, same-size record at 1% or more: no\n",
        "    DEL  chr3:1000-9000  class 4  G1  support: single caller (manta)\n         \
         benign-set AF 0.05 (gnomAD-SV_v4.1; the highest AF of the benign records listed, larger \
         SVs around it included), common\n         \
         benign records: one of about the same size, same-size record at 1% or more: yes\n",
        "  With a same-size benign record at 1% or more: 1 (the common count above uses the \
         benign-set AF, AnnotSV's aggregate)\n",
        EVIDENCE_NOTE,
    ] {
        assert!(txt.contains(want), "the text report lacks {want:?}:\n{txt}");
    }
    // The head keeps its words, which e2e cases read.
    assert!(
        txt.contains(
            "Class 4-5 SVs: 3 (single caller: not known; deletions on exon edges: 0; common: 1)"
        ),
        "{txt}"
    );
    let card = frame(&page, "structural-variants");
    assert!(
        card.contains(
            "<td>benign-set AF 0.05 (gnomAD-SV_v4.1; the highest AF of the benign records listed, \
             larger SVs around it included), common; benign records: one of about the same size, \
             same-size record at 1% or more: yes</td>"
        ),
        "{card}"
    );
    assert!(
        card.contains(&*e(
            "With a same-size benign record at 1% or more: 1 (the common count above uses the \
             benign-set AF, AnnotSV's aggregate)"
        )),
        "{card}"
    );

    // Control: the same AF on a row without the evidence keys reads as before.
    let old = sv_row(
        "chr2:2000-2500",
        r#""population_af": 0.05, "population_source": "gnomAD-SV_v4.1", "common": true"#,
    );
    let (txt, _) = reports(&annotsv(&[&old]));
    assert!(
        txt.contains("         benign-set AF 0.05 (gnomAD-SV_v4.1), common\n"),
        "{txt}"
    );
    assert!(!txt.contains("same-size"), "{txt}");
}

// --- MEI: a filtered call and a site with too few reads ---------------------------

#[test]
fn both_reports_show_a_filtered_call_and_a_site_not_assessed() {
    let s = summary(&section(
        "mei",
        "Mobile element insertions (MEGAnE)",
        "50",
        r#"{"calls": 1, "by_element": {"SVA": 1}, "sites": [
            {"gene": "FKTN", "element": "SVA", "condition": "Fukuyama congenital muscular dystrophy",
             "inheritance": "autosomal recessive", "status": "call filtered out (LC)", "mean_depth": "30.0",
             "call": "chr9:105639660 SVA_F 1 copy LC"},
            {"gene": "TAF1", "element": "SVA", "condition": "X-linked dystonia-parkinsonism",
             "inheritance": "X-linked", "status": "not callable", "mean_depth": "4.0", "call": "."}]}"#,
    ));
    let (txt, page) = reports(&s);
    assert!(
        txt.contains(
            "  FKTN SVA (Fukuyama congenital muscular dystrophy; autosomal recessive): call filtered \
             out (LC) (depth 30.0, chr9:105639660 SVA_F 1 copy LC)\n"
        ),
        "{txt}"
    );
    assert!(
        txt.contains(
            "  TAF1 SVA (X-linked dystonia-parkinsonism; X-linked): not callable (depth 4.0)\n"
        ),
        "{txt}"
    );
    let card = frame(&page, "mei");
    assert!(card.contains("<td>call filtered out (LC)</td>"), "{card}");
    assert!(card.contains("<td>not callable</td>"), "{card}");
    assert!(card.contains("Not callable: too few reads"), "{card}");
}

// --- HTT: reduced penetrance ---------------------------------------------------------

/// An expansions section whose HTT locus has `count` repeats, flagged by
/// Stranger as a full mutation, with `extra` keys on the locus.
fn htt(count: &str, extra: &str) -> Summary {
    let locus = format!(
        r#"{{"locus": "HTT", "repeat_count": "{count}", "status": "full_mutation", "normal_max": "26",
            "pathologic_min": "36", "support": "SPANNING/SPANNING", "adsp": "12/9", "ci": "17-17/{count}-{count}",
            "sizing": "sized", "allele_sizing": "sized/sized", "zone": "", "zone_label": "",
            "zone_citation": ""{extra}}}"#
    );
    summary(&section(
        "expansions",
        "Repeat expansions (ExpansionHunter)",
        "09",
        &format!(
            r#"{{"records": 1, "key_loci": [], "loci": [{locus}], "stranger": true,
                "stranger_outdated": false, "unsized": [], "intermediate": [], "readings": [],
                "flagged": [{locus}], "no_status": 0}}"#
        ),
    ))
}

const HTT_PENETRANCE: &str = r#", "penetrance": "36-39",
    "penetrance_label": "reduced penetrance (may or may not cause Huntington disease; full penetrance from 40)",
    "penetrance_citation": "GeneReviews: Huntington Disease, https://www.ncbi.nlm.nih.gov/books/NBK1305/""#;

const HTT_LINE: &str =
    "36-39: reduced penetrance (may or may not cause Huntington disease; full penetrance from 40)";

#[test]
fn an_htt_count_of_reduced_penetrance_says_so() {
    let (txt, page) = reports(&htt("17/38", HTT_PENETRANCE));
    assert!(
        txt.contains(&format!(
            "    HTT      17/38  full_mutation\n      {HTT_LINE}\n"
        )),
        "{txt}"
    );
    assert!(
        txt.contains(&format!(
            "normal max 26; {HTT_LINE}; GeneReviews: Huntington Disease, \
             https://www.ncbi.nlm.nih.gov/books/NBK1305/"
        )),
        "{txt}"
    );
    let card = frame(&page, "expansions");
    assert!(
        card.contains(&format!("<p class=\"note\">HTT {}</p>", e(HTT_LINE))),
        "{card}"
    );

    // Control: a count of 40 carries no penetrance keys and reads as before.
    let (txt, page) = reports(&htt("17/40", ""));
    assert!(!txt.contains("reduced penetrance"), "{txt}");
    assert!(!page.contains("reduced penetrance"), "{page}");
}

// --- Locityper: uncalibrated outside HLA and KIR -------------------------------------

/// A Locityper section with one HLA locus that passes and one CYP2D6
/// locus with `filter` and `reason`.
fn locityper(filter: &str, reason: &str) -> Summary {
    let locus = |name: &str, filter: &str, reason: &str| {
        format!(
            r#"{{"locus": "{name}", "region": "chr22:1-2", "status": "called", "genotype": "a.1,b.2",
                "haplotype_1": "a.1", "haplotype_2": "b.2", "quality": "40", "total_reads": "900",
                "unexplained_reads": "10", "weighted_distance": "4", "warnings": ".", "filter": "{filter}",
                "reason": "{reason}", "primary": "*1/*2"}}"#
        )
    };
    let passed = usize::from(filter == "pass") + 1;
    summary(&section(
        "locityper",
        "Locityper second opinion at hard loci",
        "53",
        &format!(
            r#"{{"called": 2, "passed": {passed}, "loci": [{}, {}]}}"#,
            locus("HLA-DRB1", "pass", "."),
            locus("CYP2D6", filter, reason)
        ),
    ))
}

const UNCALIBRATED_WHY: &str = "the authors calibrated this filter on HLA and KIR loci only; by its \
    rule alone it passes";

#[test]
fn an_uncalibrated_locus_shows_why_and_never_a_green_pass() {
    let (txt, page) = reports(&locityper("uncalibrated", UNCALIBRATED_WHY));
    assert!(
        txt.contains(
            "  Loci called: 2 of 2; passing the authors' filter: 1; uncalibrated (outside HLA and \
             KIR, so neither pass nor fail): 1\n"
        ),
        "{txt}"
    );
    assert!(
        txt.contains(&format!(
            "weighted distance 4; filter uncalibrated: {UNCALIBRATED_WHY})"
        )),
        "{txt}"
    );
    let card = frame(&page, "locityper");
    assert!(
        card.contains("<td><span class=\"badge badge-gray\">uncalibrated</span></td>"),
        "{card}"
    );
    assert_eq!(
        card.matches("<span class=\"badge badge-green\">pass</span>")
            .count(),
        1,
        "only HLA-DRB1 passes:\n{card}"
    );
    assert!(card.contains(&*e(UNCALIBRATED_WHY)), "{card}");
    assert!(
        card.contains("2 of 2 (1 pass the filter, 1 uncalibrated outside HLA and KIR)"),
        "{card}"
    );

    // Control: a table from before reads pass at CYP2D6, with no new words.
    let (txt, _) = reports(&locityper("pass", "."));
    assert!(
        txt.contains("  Loci called: 2 of 2; passing the authors' filter: 2\n"),
        "{txt}"
    );
    assert!(!txt.contains("uncalibrated"), "{txt}");
}

// --- LPA: the interval covers read sampling only ---------------------------------------

/// An LPA section estimated at 29 repeats, with step 46's `message`.
fn lpa_summary(message: &str) -> Summary {
    summary(&section(
        "lpa",
        "Lp(a): LPA KIV-2 repeat count and risk variants",
        "46",
        &format!(
            r#"{{"status": "estimated", "kiv2_copies": "29.4", "ci95_low": "23.8", "ci95_high": "35.8",
                "class": "typical", "kiv2_depth": "65.36", "flank_depth": "29.28", "mapq0_share": "0.900",
                "message": "{message}", "advice": "a blood test measures Lp(a)", "null_allele": null,
                "variants": [], "flag": "", "flag_text": "no genetic sign of a high Lp(a)"}}"#
        ),
    ))
}

const SAMPLING: &str = "about 29 KIV-2 repeats in total over both copies (95% interval for read \
    sampling, 23.8 to 35.8; aligner or library bias not included): within or overlapping the range \
    most people have";

const CALIBRATED: &str = "about 29 KIV-2 repeats in total over both copies (95% interval 23.8 to \
    35.8, read sampling plus the error measured against assembled genomes): within or overlapping \
    the range most people have";

#[test]
fn the_kiv2_headline_says_the_interval_is_read_sampling() {
    let head = "KIV-2 repeats (both chromosomes): 29.4 (95% interval for read sampling only, 23.8 \
                to 35.8), typical";
    let (txt, page) = reports(&lpa_summary(SAMPLING));
    assert!(txt.contains(&format!("  {head}\n")), "{txt}");
    assert!(
        frame(&page, "lpa").contains(&format!("<p>{head}</p>")),
        "{page}"
    );

    let (txt, _) = reports(&lpa_summary(CALIBRATED));
    assert!(
        txt.contains(
            "KIV-2 repeats (both chromosomes): 29.4 (95% interval for read sampling and the error \
             measured against assembled genomes, 23.8 to 35.8), typical"
        ),
        "{txt}"
    );

    // Control: a message from before step 46 said what it covers.
    let (txt, _) = reports(&lpa_summary("ok"));
    assert!(
        txt.contains("KIV-2 repeats (both chromosomes): 29.4 (95% interval 23.8 to 35.8), typical"),
        "{txt}"
    );
}
