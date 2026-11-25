//! The report's links to the genome viewer: every findings table of the HTML
//! report ends each row with a Viewer cell that opens the row's position,
//! with the viewer note under the table. Moved from check 5 of
//! tests/test_viewer.py, which went with the Python report and live server.
//! A control strips the links from the same page: the check must then name
//! every table, so it can fail.

use og_render::html::{self, VIEWER_URL};
use og_render::model::{Cell, Summary};
use og_render::notes::VIEWER_NOTE;
use og_render::sections::clinvar::viewer_of;

/// Each findings table, by the start of its header row, and the viewer
/// address its one row (made by [`findings_summary`]) must link to.
const FINDINGS_TABLES: [(&str, &str, &str); 10] = [
    ("ClinVar hits", "<th>Chr</th><th>Position</th>", "chr1:451"),
    (
        "secondary findings",
        "<th>Gene</th><th>Variant</th><th>Genotype</th><th>Evidence</th>",
        "chr2:120",
    ),
    (
        "mitochondrial variants",
        "<th>Variant</th><th>AF</th><th>ClinVar</th>",
        "chrM:3243",
    ),
    (
        "deep intronic sites",
        "<th>Gene</th><th>Variant</th><th>Position</th><th>Call</th>",
        "chr1:900",
    ),
    (
        "cDNA scan",
        "<th>Locus</th><th>Junctions</th><th>Reads</th>",
        "chr21:100-900",
    ),
    (
        "AnnotSV class 4-5",
        "<th>SV</th><th>Class</th><th>Genes</th>",
        "chr1:1000-1500",
    ),
    (
        "AnnotSV moved to class 3",
        "<th>SV</th><th>AnnotSV class</th>",
        "chr20:2000-2400",
    ),
    (
        "slivar pairs in trans",
        "<th>Gene</th><th>Pair in trans</th>",
        "chr1:50-90",
    ),
    (
        "mosaic alterations",
        "<th>Type</th><th>Arm</th><th>Region</th>",
        "chr7:1-5000000",
    ),
    (
        "Exomiser ranking",
        "<th>Rank</th><th>Gene</th><th>Inheritance</th>",
        "chr3:700",
    ),
];

/// A summary with one row in every findings table.
fn findings_summary() -> String {
    let sv = r#""type": "DEL", "genes": ["GENEF"], "caller": "manta", "support": 2,
        "support_callers": ["manta", "delly"], "exon_edge": "no", "exon_edge_detail": "",
        "population_af": null, "population_source": "", "common": false"#;
    let sections = [
        (
            "clinvar",
            "ClinVar",
            "06",
            r#"{"hits": [{"chrom": "chr1", "pos": 451, "ref": "T", "alt": "G", "genotype": "het", "gene": "GENEB",
              "significance": "Pathogenic", "review_status": "reviewed by expert panel", "stars": 3}],
             "deep_intronic": [{"gene": "GENEE", "hgvs": "NM_000001.1:c.100-500A>G", "chrom": "chr1", "pos": 900,
              "ref": "A", "alt": "G", "call": "het", "detail": "GQ 40"}],
             "cdna_scan": {"note": "", "loci": [{"chrom": "chr21", "start": 100, "end": 900, "junctions": 3, "reads": 24}]}}"#
                .to_owned(),
        ),
        (
            "clinical",
            "Clinical",
            "23",
            r#"{"variants": 1, "acmg_sf": {"version": "v3.3", "high_impact": [],
              "clinvar_hits": [{"gene": "GENEC", "variant": "chr2:120 A>G", "genotype": "0/1", "significance": "Pathogenic"}]}}"#
                .to_owned(),
        ),
        (
            "mito",
            "Mitochondria",
            "20",
            r#"{"pass": 1, "heteroplasmic": 0, "annotation": {"source": "synthetic",
              "variants": [{"pos": 3243, "ref": "A", "alt": "G"}]}}"#
                .to_owned(),
        ),
        (
            "annotsv",
            "AnnotSV",
            "05",
            format!(
                r#"{{"count": 1, "single_caller": 0, "consensus": true, "exon_edge": 0, "exon_edge_checked": true, "common": 0,
              "rows": [{{{sv}, "rule": "", "rule_detail": "", "sv": "chr1:1000-1500", "acmg_class": "5", "annotsv_class": "5"}}],
              "downgraded": [{{{sv}, "sv": "chr20:2000-2400", "acmg_class": "3", "annotsv_class": "4",
                "rule": "ig_tr_locus", "rule_detail": "inside a locus"}}]}}"#
            ),
        ),
        (
            "slivar",
            "slivar",
            "31",
            r#"{"prioritized": 1, "compound_het_pairs": 1, "compound_het_phase": {"trans": 1, "cis": 0, "unknown": 0},
              "compound_het_trans": [{"gene": "G2", "variants": "chr1:50 A>G + chr1:90 C>T", "reads": "19 of 20"}]}"#
                .to_owned(),
        ),
        (
            "mosaic",
            "Mosaic alterations",
            "51",
            r#"{"status": "assessed", "calls_reported": "1", "floor": "0.05", "depth": "30",
              "calls": [{"type": "Loss", "arms": "7q", "chrom": "chr7", "start": "1", "end": "5000000",
                "length_mb": "5.0", "cell_fraction": "0.2"}]}"#
                .to_owned(),
        ),
        (
            "exomiser",
            "Exomiser",
            "48",
            r#"{"terms": [{"id": "HP:0001250", "label": "Seizure"}], "genes": [
              {"rank": "1", "gene": "GENEG", "moi": "AD", "exomiser_score": "0.9", "phenotype_score": "0.8",
               "variants": [{"variant": "chr3:700 C>T", "genotype": "0/1", "consequence": "missense_variant",
                 "hgvs": "p.(Arg1Cys)", "acmg": "LIKELY_PATHOGENIC", "disease": "a disease"}]}]}"#
                .to_owned(),
        ),
    ];
    let mut s = String::new();
    for (key, title, step, data) in sections {
        s += &format!(
            r#", "{key}": {{"title": "{title}", "state": "ok", "step": "{step}", "source": "{key}/S_{key}.tsv",
             "file_date": "2026-11-06 11:00:00", "note": null, "data": {data}}}"#
        );
    }
    format!(
        r#"{{
 "schema_version": 2, "sample": "S", "generated_utc": "2026-11-06 12:00 UTC",
 "run": {{"started_utc": "2026-11-06 10:00 UTC", "declared_sex": "male"}},
 "databases": {{"clinvar_release": "2026-10-01", "cpsr_clinvar_release": null, "hla_database": null}},
 "manifest": [],
 "sections": {{
  "coverage": {{"title": "Coverage (mosdepth)", "state": "ok", "step": "16b", "source": "qc/S.mosdepth.summary.txt",
               "file_date": "2026-11-06 11:00:00", "note": null, "data": {{"mean_depth": 31.25}}}}{s}
 }},
 "cyp2d6": null, "not_run": [], "not_assessed": []
}}"#
    )
}

fn page() -> String {
    html::report(&Summary::parse(&findings_summary()).expect("a summary")).expect("the HTML report")
}

/// Whether `row` ends with a cell that is one viewer link.
fn ends_with_link(row: &str) -> bool {
    let Some(at) = row.rfind("<td>") else {
        return false;
    };
    let Some(inner) = row[at..]
        .strip_prefix("<td><a class=\"locus\" href=\"")
        .and_then(|c| c.strip_suffix("</a></td>"))
    else {
        return false;
    };
    let Some((attrs, text)) = inner.split_once('>') else {
        return false;
    };
    !attrs.contains('<') && !text.is_empty() && !text.contains('<')
}

/// The name of each findings table missing from `html`, without the Viewer
/// column last, with a row whose last cell is not a viewer link to the
/// expected address, or without the viewer note in its card.
fn tables_without_links(html: &str) -> Vec<String> {
    let mut bad = Vec::new();
    let chunks: Vec<&str> = html.split("<table").skip(1).collect();
    for (name, head, at) in FINDINGS_TABLES {
        let Some(chunk) = chunks
            .iter()
            .find(|c| c.split("</tr>").next().unwrap_or("").contains(head))
        else {
            bad.push(format!("{name}: no table"));
            continue;
        };
        let table = chunk.split("</table>").next().unwrap_or("");
        let after = &html[html.find(table).unwrap() + table.len()..];
        if !table
            .split("</tr>")
            .next()
            .unwrap_or("")
            .ends_with("<th>Viewer</th>")
        {
            bad.push(format!(
                "{name}: the header does not end with the Viewer column"
            ));
        }
        let rows: Vec<&str> = table
            .split("<tr><td")
            .skip(1)
            .map(|r| r.rsplit_once("</tr>").map_or(r, |(a, _)| a))
            .collect();
        let link = format!(
            "<td><a class=\"locus\" href=\"{VIEWER_URL}#{at}\" target=\"_blank\" rel=\"noopener noreferrer\">{at}</a></td>"
        );
        if rows.is_empty()
            || rows.iter().any(|r| !ends_with_link(r))
            || !rows.iter().any(|r| r.ends_with(&link))
        {
            bad.push(format!(
                "{name}: a row whose last cell is not a viewer link to {at}"
            ));
        }
        if !after
            .split("<div class=\"card")
            .next()
            .unwrap_or("")
            .contains(VIEWER_NOTE.trim())
        {
            bad.push(format!(
                "{name}: no viewer note in the card under the table"
            ));
        }
    }
    bad
}

#[test]
fn every_findings_table_ends_each_row_with_a_viewer_link() {
    let html = page();
    assert_eq!(tables_without_links(&html), Vec::<String>::new());
    assert!(
        VIEWER_NOTE.contains("scripts/owngenome-view.sh"),
        "{VIEWER_NOTE}"
    );
}

#[test]
fn control_with_the_links_removed_the_check_names_every_table() {
    // Each viewer link becomes an empty cell, as a report without links would have.
    let mut plain = String::new();
    let mut rest = page();
    while let Some(at) = rest.find("<td><a class=\"locus\"") {
        plain.push_str(&rest[..at]);
        let end = rest[at..].find("</a></td>").unwrap() + at + "</a></td>".len();
        plain.push_str("<td></td>");
        rest = rest[end..].to_owned();
    }
    plain.push_str(&rest);
    let missed = tables_without_links(&plain);
    assert_eq!(missed.len(), FINDINGS_TABLES.len(), "{missed:?}");
    assert!(
        missed.iter().all(|m| m.contains("not a viewer link")),
        "{missed:?}"
    );
}

#[test]
fn a_rows_viewer_cell_spans_its_first_chromosome_and_is_empty_without_one() {
    let cell = |t: &str| match viewer_of(t).expect("a cell") {
        Cell::Viewer(at) => at,
        _ => panic!("not a Viewer cell"),
    };
    assert_eq!(
        cell("chr1:50 A>G + chr1:90 C>T").as_deref(),
        Some("chr1:50-90")
    );
    assert_eq!(cell("chrM:3243 A>G").as_deref(), Some("chrM:3243"));
    assert_eq!(cell("no locus <x>"), None);
}
