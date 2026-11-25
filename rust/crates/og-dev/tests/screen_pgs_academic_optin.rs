//! `og dev pgs licence` against the PGS Catalog's saved answers for the 15
//! scores of `assets/pgs_scores.tsv` (`tests/fixtures/screen-pgs-academic-optin/`,
//! the REST `id`, `trait_reported` and `license` fields as fetched on
//! 2026-11-25). Every row's `licence` column must be what the catalog says,
//! and a copy with one row's licence flipped must not pass.

use og_dev::pgs::licence_of;
use std::collections::HashMap;
use std::path::{Path, PathBuf};

fn repo(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../..").join(rel)
}

/// pgs_id -> licence column of a score list's text.
fn column(list: &str) -> Vec<(String, String)> {
    let mut idx = None;
    let mut out = Vec::new();
    for line in list.lines().filter(|l| !l.starts_with('#')) {
        let c: Vec<&str> = line.split('\t').collect();
        match idx {
            None => idx = c.iter().position(|h| *h == "licence"),
            Some(i) => out.push((c[0].to_owned(), c.get(i).copied().unwrap_or("").to_owned())),
        }
    }
    assert!(idx.is_some(), "the list has no licence column");
    out
}

/// The saved answer of each id.
fn answers() -> HashMap<String, String> {
    let text = std::fs::read_to_string(repo(
        "tests/fixtures/screen-pgs-academic-optin/score-licences.jsonl",
    ))
    .expect("saved answers");
    text.lines()
        .map(|l| {
            let id = l.split("\"id\":\"").nth(1).and_then(|r| r.split('"').next());
            (id.expect("an id").to_owned(), l.to_owned())
        })
        .collect()
}

/// The rows whose licence column is not what the saved answer says.
fn wrong(list: &str) -> Vec<String> {
    let a = answers();
    column(list)
        .into_iter()
        .filter(|(id, lic)| a.get(id).and_then(|b| licence_of(b, id)) != Some(lic.as_str()))
        .map(|(id, lic)| format!("{id} {lic}"))
        .collect()
}

#[test]
fn every_licence_of_the_list_is_the_catalogs() {
    let list = std::fs::read_to_string(repo("assets/pgs_scores.tsv")).expect("list");
    assert_eq!(column(&list).len(), 15);
    assert_eq!(wrong(&list), Vec::<String>::new());
    let academic: Vec<String> = column(&list)
        .into_iter()
        .filter(|(_, l)| l == "academic-research")
        .map(|(id, _)| id)
        .collect();
    assert_eq!(academic, ["PGS000014", "PGS000016", "PGS000027", "PGS000017"]);
}

#[test]
fn a_flipped_row_and_a_new_licence_text_fail() {
    let list = std::fs::read_to_string(repo("assets/pgs_scores.tsv")).expect("list");
    // One row flipped each way.
    let flipped = list
        .replacen("\tacademic-research\n", "\tauthor\n", 1)
        .replacen("Breast cancer\tfemale\tPPM000004\tEUR\tOR 1.61\tauthor", "Breast cancer\tfemale\tPPM000004\tEUR\tOR 1.61\tacademic-research", 1);
    assert_eq!(wrong(&flipped), ["PGS000014 author", "PGS000004 academic-research"]);
    // An empty cell is no licence.
    let empty = list.replacen("\tauthor\n", "\t\n", 1);
    assert_eq!(wrong(&empty), ["PGS000018 "]);
}

#[test]
fn the_reader_knows_two_texts_and_nothing_else() {
    let author = r#"{"id": "PGS1", "license": "PGS obtained from the Catalog should be cited appropriately, and used in accordance with any licensing restrictions set by the authors. See EBI Terms of Use (https://www.ebi.ac.uk/about/terms-of-use/) for additional details."}"#;
    assert_eq!(licence_of(author, "PGS1"), Some("author"));
    // Line breaks and doubled spaces (JSON escapes in the answer) do not change the text.
    assert_eq!(
        licence_of(&author.replace("cited ", "cited\\n  "), "PGS1"),
        Some("author")
    );
    let academic = r#"{"id": "PGS1", "license": "Freely available to the academic community for research use. Parties interested in using the scores for commercial purposes should contact the Broad Office."}"#;
    assert_eq!(licence_of(academic, "PGS1"), Some("academic-research"));
    for body in [
        // a restriction the catalog adds to the default text
        &author.replace("additional details.", "additional details. Non-commercial use only."),
        r#"{"id": "PGS1", "license": "CC BY-NC-SA 4.0"}"#,
        r#"{"id": "PGS1", "license": ""}"#,
        r#"{"id": "PGS1"}"#,
        r#"{"id": "PGS2", "license": "Freely available to the academic community for research use."}"#,
        "{}",
        "<html>",
    ] {
        assert_eq!(licence_of(body, "PGS1"), None, "{body}");
    }
}
