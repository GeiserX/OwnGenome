//! `assets/traits.tsv`, the curated traits panel of step 42: it holds at
//! least 30 rows and passes `og traits check`, and a copy with a planted
//! wrong row fails it. The row-by-row schema rules are the unit tests in
//! `og_steps::traits`; the reference bases, PubMed ids and dbSNP records are
//! `og dev traits-sources`, which the e2e case `traits-panel-1` runs on the
//! fixture reference.
//!
//! This was the schema part of `tests/test_traits_table.py`.

use std::path::PathBuf;

const REPO: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../..");

fn table() -> String {
    format!("{REPO}/assets/traits.tsv")
}

/// The table's data rows: no `#` line, no blank line, no header.
fn rows(text: &str) -> Vec<&str> {
    text.lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .skip(1)
        .collect()
}

fn check(path: &str) -> u8 {
    og_steps::traits::cli(&["check".to_owned(), "--table".to_owned(), path.to_owned()])
}

fn scratch(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("og-traits-table-{}-{name}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

#[test]
fn the_table_has_at_least_30_rows() {
    let text = std::fs::read_to_string(table()).unwrap();
    let n = rows(&text).len();
    assert!(n >= 30, "assets/traits.tsv has {n} rows");
}

#[test]
fn the_table_passes_og_traits_check() {
    assert_eq!(check(&table()), 0);
}

#[test]
fn control_a_planted_bad_row_fails_og_traits_check() {
    let text = std::fs::read_to_string(table()).unwrap();
    let lactase = rows(&text)
        .into_iter()
        .find(|l| l.split('\t').nth(2) == Some("rs4988235"))
        .expect("the lactase row rs4988235");
    // The same rsID twice and a reference base that is no base.
    let mut cols: Vec<&str> = lactase.split('\t').collect();
    cols[4] = "N";
    let planted = format!("{text}{}\n", cols.join("\t"));
    let d = scratch("planted");
    let p = d.join("traits.tsv");
    std::fs::write(&p, planted).unwrap();
    assert_eq!(check(&p.to_string_lossy()), 1);
    let _ = std::fs::remove_dir_all(&d);
}
