//! The PharmCAT reports of `tests/fixtures/pharmcat`, read by the library:
//! the checks of `tests/test_cpic_parser.py` that need more than the
//! command's output.

use og_pgx::guidance::drug_guidance;
use og_pgx::report::parse_genes;
use og_pgx::{Parsed, Status};
use pycompat::json::{Value, loads};

fn fixture(name: &str) -> Value {
    let path = format!(
        "{}/../../../tests/fixtures/pharmcat/{name}",
        env!("CARGO_MANIFEST_DIR")
    );
    loads(&std::fs::read_to_string(&path).unwrap()).unwrap()
}

fn calls(data: &Value) -> Vec<og_pgx::GeneCall> {
    match parse_genes(data) {
        Parsed::Ok(c) => c,
        other => panic!("{other:?}"),
    }
}

#[test]
fn the_hg002_reports_parse_and_list_their_possible_diplotypes() {
    for version in ["3.2.0", "3.4.0"] {
        let data = fixture(&format!("report-{version}.json"));
        let c = calls(&data);
        assert!(c.len() >= 10, "{version}: {} genes", c.len());
        let several: Vec<(&str, usize)> = c
            .iter()
            .filter(|c| c.labels.len() > 1)
            .map(|c| (c.gene.as_str(), c.labels.len()))
            .collect();
        assert_eq!(several, [("CYP2B6", 4), ("CYP2C19", 528)], "{version}");
        for (gene, _) in several {
            let g = c.iter().find(|c| c.gene == gene).unwrap();
            assert_eq!(g.status(), Status::Ambiguous, "{version} {gene}");
        }
        let guidance = drug_guidance(&data);
        assert!(
            guidance.values().map(Vec::len).sum::<usize>() > 0,
            "{version}: the drugs section is read"
        );
    }
}

#[test]
fn the_docs_example_has_three_non_normal_and_seven_unclassified_genes() {
    let c = calls(&fixture("pharmcat-docs-example.json"));
    let with = |s: Status| -> Vec<&str> {
        let mut g: Vec<&str> = c
            .iter()
            .filter(|c| c.status() == s)
            .map(|c| c.gene.as_str())
            .collect();
        g.sort_unstable();
        g
    };
    assert_eq!(with(Status::NonNormal), ["CYP2D6", "HLA-B", "NAT2"]);
    assert_eq!(
        with(Status::Unclassified),
        [
            "CACNA1S", "CFTR", "CYP4F2", "IFNL3", "MT-RNR1", "RYR1", "VKORC1"
        ]
    );
}
