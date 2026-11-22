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

#[test]
fn traits_show_low_quality_rows_as_unread() {
    let data = r#"{"called": 0, "not_covered": 0, "other_allele": 0, "low_gq": 1, "rows": [
        {"rsid": "rs1", "trait": "Example trait", "gene": "GENE", "genotype": "A/G",
         "status": "low_gq", "reading": "Do not show this interpretation", "note": "weak genotype call",
         "ancestry": "yes", "pmids": "1234"}]}"#;
    for report in reports("traits", data) {
        assert!(report.contains("Low-quality genotype calls"));
        assert!(report.contains("Example trait") && report.contains("not read"));
        assert!(!report.contains("Do not show this interpretation"));
        assert!(!report.contains("A/G"));
    }
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
        assert!(report.contains("also fit an inherited heterozygous variant"));
        assert!(report.contains("not a diagnosis"));
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
