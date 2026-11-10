//! Every path og-map reads under GENOME_DIR is a `file` row of
//! `assets/output_contract.tsv`, so an output renamed there breaks this test
//! instead of silently turning a map layer to `missing`.

use std::path::Path;

#[test]
fn every_input_path_is_a_contract_row() {
    let contract =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../../assets/output_contract.tsv");
    let text = std::fs::read_to_string(&contract).expect("assets/output_contract.tsv reads");
    let rows: Vec<&str> = text
        .lines()
        .filter(|l| !l.starts_with('#'))
        .filter_map(|l| {
            let mut f = l.split('\t');
            (f.next() == Some("file")).then(|| f.next()).flatten()
        })
        .collect();
    assert!(rows.len() > 100, "the contract lists {} files", rows.len());
    let missing: Vec<&&str> = og_map::INPUT_PATHS
        .iter()
        .filter(|p| !rows.contains(p))
        .collect();
    assert!(
        missing.is_empty(),
        "og-map reads paths that are not in assets/output_contract.tsv: {missing:?}"
    );
}
