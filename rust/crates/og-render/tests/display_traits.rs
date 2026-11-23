use og_render::model::Summary;
use og_render::{html, text};

fn reports(section: &str, data: &str) -> [String; 2] {
    let summary = Summary::parse(&format!(
        r#"{{
            "schema_version": 2, "sample": "S1", "generated_utc": "2026-01-02 03:04 UTC",
            "run": {{"started_utc": "2026-01-02 02:00 UTC", "declared_sex": "male"}},
            "databases": {{"clinvar_release": null, "cpsr_clinvar_release": null, "hla_database": null}},
            "cyp2d6": null, "manifest": [], "not_run": [], "not_assessed": [],
            "sections": {{
                "{section}": {{"title": "{section}", "step": "42", "state": "ok",
                    "source": "{section}/S1.tsv", "file_date": null, "note": null, "data": {data}}}
            }}
        }}"#
    ))
    .unwrap();
    [
        text::report(&summary).unwrap(),
        html::report(&summary).unwrap(),
    ]
}

fn assert_unread_trait_reading(reading: &str) {
    let data = format!(
        r#"{{"called": 0, "not_covered": 0, "other_allele": 0, "low_gq": 1, "rows": [
        {{"rsid": "rs1", "trait": "Example trait", "gene": "GENE", "genotype": "A/G",
         "status": "low_gq", "reading": "{reading}", "note": "Study note",
         "ancestry": "yes", "pmids": "1234"}}]}}"#
    );
    let [text, html] = reports("traits", &data);
    for report in [&text, &html] {
        assert!(report.contains("Low-quality genotype calls"));
        assert!(!report.contains("A/G"));
    }
    assert_eq!(
        [
            text.contains(&format!("  Example trait (GENE rs1): {reading}\n")),
            html.contains(&format!("</span></td><td>.</td><td>{reading}<br>")),
        ],
        [true, true],
        "text and HTML should both show the producer's reading: {reading}"
    );
}

#[test]
fn traits_show_the_genotype_quality_reason() {
    assert_unread_trait_reading(
        "A low-quality call at this position (genotype quality 12, below 20): not read; check the reads.",
    );
}

#[test]
fn traits_show_the_caller_filter_reason() {
    assert_unread_trait_reading(
        "A low-quality call at this position (the caller filtered the call: LowQual): not read; check the reads.",
    );
}

#[test]
fn blood_groups_display_the_supplied_notes() {
    let data = r#"{"abo": "A1", "rhd": "not determined", "systems": [
        {"system": "ABO", "phenotype": "A1", "genotype": "A1/O.01", "confidence": "high",
         "note": "Rare ABO alleles are outside this panel"},
        {"system": "RhD", "phenotype": "not determined", "genotype": ".", "confidence": ".",
         "note": "RhD needs read depth"},
        {"system": "Duffy", "phenotype": "Fy(a+b-)", "genotype": "FY*A/FY*B", "confidence": "medium",
         "note": "Duffy phase was unknown"},
        {"system": "Kell", "phenotype": "K-k+", "genotype": "k/k", "confidence": "high",
         "note": "."}]}"#;
    for report in reports("bloodtype", data) {
        for note in [
            "Rare ABO alleles",
            "RhD needs read depth",
            "Duffy phase was unknown",
        ] {
            assert!(report.contains(note), "missing {note}");
        }
    }
}

#[test]
fn chip_shows_the_assessable_share_and_preserves_the_call_note() {
    let data = r#"{"calls": [
        {"gene": "GENE", "protein_change": "p.X1Y", "vaf": "0.300", "alt_reads": "9",
         "depth": "30", "note": "possible clonal haematopoiesis (not a diagnosis): 9 of 30 reads also fit an inherited heterozygous variant", "status": "chip"}
        ], "median_depth": 30, "median_limit": 0.18, "genes": [
        {"gene": "GENE", "assessable_fraction": "0.800"},
        {"gene": "ZERO", "assessable_fraction": "0.000"},
        {"gene": "UNQUANTIFIED", "assessable_fraction": "."}
        ]}"#;
    for report in reports("chip", data) {
        assert!(report.contains("GENE 80%") && report.contains("ZERO 0%"));
        assert!(report.contains(
            "GENE p.X1Y: VAF 0.300 (9 of 30 reads), possible clonal haematopoiesis (not a diagnosis): 9 of 30 reads also fit an inherited heterozygous variant"
        ));
        assert!(!report.contains("UNQUANTIFIED 0%"));
    }
}

#[test]
fn viral_reads_below_threshold_do_not_read_as_a_positive_or_a_negative() {
    let data = r#"{"hhv6": "no integrated HHV-6", "human_depth": "30.00", "genomes": 24, "found": [
        {"virus": "HHV-6B", "call": "below threshold", "reads": "7", "breadth": "0.0041",
         "relative_depth": "0.000", "copies_per_cell": ".",
         "note": "7 counted reads, under the 10 a call needs: too few to tell the virus from contamination at this depth."}]}"#;
    for report in reports("viral", data) {
        assert!(report.contains("below threshold"));
        assert!(report.contains("under the 10 a call needs"));
        assert!(report.contains("no integration detected"));
        assert!(!report.contains("No read on any"));
    }
}

#[test]
fn chip_explains_why_the_whole_gene_has_no_read_support_limit() {
    let data = r#"{"median_depth": 60, "median_limit": 0.09,
    "genes": [{"gene": "DNMT3A", "assessable_fraction": "0.985"}],
    "not_covered": ["JAK2"]}"#;
    for report in reports("chip", data) {
        assert!(report.contains(
            "No read-support limit for the whole gene (too shallow or not on the reference): JAK2."
        ));
    }
}

#[test]
fn ebv_reads_below_threshold_preserve_the_hhv6_headline() {
    let data = r#"{"hhv6": "no integrated HHV-6", "human_depth": "30.00", "genomes": 24, "found": [
        {"virus": "EBV", "call": "below threshold", "reads": "4", "breadth": "0.0041",
         "relative_depth": "0.000", "copies_per_cell": ".",
         "note": "4 counted reads, under the 10 a call needs: too few to tell the virus from contamination at this depth."}]}"#;
    let [text, html] = reports("viral", data);
    assert_eq!(
        [
            text.contains("HHV-6: no integrated HHV-6"),
            html.contains(r#"HHV-6</span><span class="value">no integrated HHV-6</span>"#),
        ],
        [true, true]
    );
    for report in [text, html] {
        assert!(report.contains("EBV: 4 counted reads"));
    }
}
