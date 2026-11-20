//! A run with the ancestry panel whose input lists variants only: og prs
//! table withholds every percentile and says why in `percentile_note`, so
//! the summary has no percentile and `adjusted` is false. Both reports must
//! say the percentile was withheld, not that no panel was installed. A run
//! without the panel (no note on any score) keeps the no-panel line.

use std::path::Path;

use og_render::model::Summary;
use og_render::{html, text};

const WITHHELD: &str = "Percentile withheld: the input lists variants only";
const NO_PANEL: &str = "no ancestry reference panel was installed";

/// The golden summary with every percentile gone and `adjusted` false; its
/// T2D score keeps its note when `note` is set.
fn raw(note: bool) -> Summary {
    let p = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../../tests/parity/golden/summary-collect/all/all.json");
    let mut j = std::fs::read_to_string(p).unwrap();
    for (from, to) in [
        ("\"adjusted\": true", "\"adjusted\": false"),
        ("\"percentile\": \"75\"", "\"percentile\": \"\""),
    ] {
        assert_eq!(j.matches(from).count(), 1, "{from}");
        j = j.replace(from, to);
    }
    if !note {
        let from = "\"percentile_note\": \"percentile withheld: the input lists variants only\"";
        assert_eq!(j.matches(from).count(), 1);
        j = j.replace(from, "\"percentile_note\": \"\"");
    }
    Summary::parse(&j).unwrap()
}

#[test]
fn a_withheld_percentile_is_not_called_a_missing_panel() {
    let s = raw(true);
    let txt = text::report(&s).unwrap();
    let page = html::report(&s).unwrap();
    for r in [&txt, &page] {
        assert!(r.contains(WITHHELD), "{r}");
        assert!(!r.contains(NO_PANEL), "{r}");
        assert!(
            !r.contains("scripts/setup.sh --ancestry-panel adds it"),
            "{r}"
        );
    }
}

#[test]
fn no_note_keeps_the_no_panel_line() {
    let s = raw(false);
    let txt = text::report(&s).unwrap();
    let page = html::report(&s).unwrap();
    for r in [&txt, &page] {
        assert!(r.contains(&format!("Raw score only: {NO_PANEL}")), "{r}");
        assert!(!r.contains(WITHHELD), "{r}");
    }
}
