//! A sample whose genetic ancestry sits between two groups (`mixed`, with
//! `second_population`, from og summary collect) is ranked among the first
//! group only. Both reports, and the live site (which reads the same note),
//! name the two groups, as step 25's own note does. A low-confidence match
//! that is not between two groups keeps the short caveat.

use std::path::Path;

use og_render::model::Summary;
use og_render::sections::prs;
use og_render::{html, text};
use pycompat::json::{Value, loads};

const BETWEEN: &str = "the genetic ancestry sits between EUR and AFR, and the percentile ranks it among EUR samples only";
const SHORT: &str = "The ancestry match is low-confidence: read the percentile with care.";

/// The golden summary (population EUR, mixed True, second_population AFR)
/// with its `mixed` flag set to `mixed`.
fn summary(mixed: &str) -> String {
    let p = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../../tests/parity/golden/summary-collect/all/all.json");
    let j = std::fs::read_to_string(p).unwrap();
    let from = "\"mixed\": \"True\"";
    assert_eq!(j.matches(from).count(), 1);
    j.replace(from, &format!("\"mixed\": \"{mixed}\""))
}

#[test]
fn a_sample_between_two_groups_is_named_between_them() {
    let j = summary("True");
    let s = Summary::parse(&j).unwrap();
    let at = |v: &Value, k: &str| match v {
        Value::Object(m) => m.get(k).cloned().unwrap(),
        other => panic!("{other:?}"),
    };
    let root = loads(&j).unwrap();
    let data = at(&at(&at(&root, "sections"), "prs"), "data");
    let note = prs::note(&data).unwrap();
    for r in [text::report(&s).unwrap(), html::report(&s).unwrap(), note] {
        assert!(r.contains(BETWEEN), "{r}");
        assert!(!r.contains(SHORT), "{r}");
    }
}

#[test]
fn a_low_confidence_match_not_between_two_groups_keeps_the_short_caveat() {
    let j = summary("False");
    let s = Summary::parse(&j).unwrap();
    for r in [text::report(&s).unwrap(), html::report(&s).unwrap()] {
        assert!(r.contains(SHORT), "{r}");
        assert!(!r.contains("sits between"), "{r}");
    }
}
