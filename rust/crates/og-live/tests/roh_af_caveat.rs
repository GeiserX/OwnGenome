//! The ROH caveat of the live site follows the allele frequencies step 11
//! used (`af_used`), as the reports do: the results page's ROH card and the
//! home page's ROH line both say which panel group, and only a run with no
//! frequencies (no `af_used`, or `0.4`) says none was used. Synthetic
//! summaries.

use og_live::state;
use og_render::model::Summary;
use og_render::notes;
use pycompat::json::{self, Value};

/// Step 11's section with `af_used`, or none.
fn roh(af_used: Option<&str>) -> String {
    let af = af_used.map_or(String::new(), |a| format!(r#", "af_used": "{a}""#));
    format!(
        r#"{{"state": "ok", "title": "Runs of Homozygosity", "data": {{"segments": 12, "total_mb": 40.5,
        "largest_mb": 6.2, "long_mb": 6.2, "f_roh": 0.002, "autosomal_over_5mb": [{{"region": "chr3:1-6200000", "mb": 6.2}}]{af}}}}}"#
    )
}

fn summary(sec: &str) -> String {
    format!(r#"{{"databases": {{}}, "sections": {{"roh": {sec}}}}}"#)
}

/// The caveats of the results page's ROH card.
fn card_caveats(sec: &str) -> Vec<String> {
    let v = json::loads(&summary(sec)).unwrap();
    let out = og_live::results::build("S", &v, None).unwrap();
    let Value::Object(root) = &out else { panic!("{out:?}") };
    let Some(Value::Array(chapters)) = root.get("chapters") else { panic!("no chapters") };
    for ch in chapters {
        let Value::Object(ch) = ch else { continue };
        let Some(Value::Array(cards)) = ch.get("sections") else { continue };
        for card in cards {
            let Value::Object(card) = card else { continue };
            if !matches!(card.get("keys"), Some(Value::Array(k)) if k.iter().any(|k| matches!(k, Value::Str(s) if s == "roh"))) {
                continue;
            }
            let Some(Value::Array(c)) = card.get("caveat") else { panic!("{card:?}") };
            return c
                .iter()
                .map(|c| match c {
                    Value::Str(s) => s.clone(),
                    other => panic!("{other:?}"),
                })
                .collect();
        }
    }
    panic!("no ROH card in {out:?}")
}

/// The home page's ROH line.
fn home_line(sec: &str) -> String {
    let s = Summary::parse(&summary(sec)).unwrap();
    let raw = json::loads(sec).unwrap();
    state::headline("roh", Some(&raw), &s).unwrap().0
}

const NONE_USED: &str = "No population allele frequencies were used";

#[test]
fn a_run_with_a_panel_group_names_it_on_both_pages() {
    let sec = roh(Some("nfe from the gnomAD v4.1 panel"));
    let want = notes::text(notes::roh_note_key(Some(&Value::Str("nfe".into())))).unwrap();
    assert!(!want.contains(NONE_USED), "{want}");
    let caveats = card_caveats(&sec);
    assert!(caveats.iter().any(|c| c == want), "{caveats:?}");
    assert!(!caveats.iter().any(|c| c.contains(NONE_USED)), "{caveats:?}");
    let line = home_line(&sec);
    assert!(line.contains(want), "{line}");
    assert!(!line.contains(NONE_USED), "{line}");
}

#[test]
fn the_global_panel_and_an_unknown_mode_each_get_their_own_caveat() {
    for (af, key) in [("global", "ROH_GLOBAL_NOTE"), ("xyz", "ROH_UNKNOWN_AF_NOTE")] {
        let want = notes::text(key).unwrap();
        let sec = roh(Some(af));
        assert!(card_caveats(&sec).iter().any(|c| c == want), "{af}");
        assert!(home_line(&sec).contains(want), "{af}");
    }
}

#[test]
fn a_run_with_no_frequencies_says_none_was_used() {
    for af in [None, Some("0.4")] {
        let sec = roh(af);
        assert!(
            card_caveats(&sec).iter().any(|c| c == notes::ROH_NOTE),
            "{af:?}"
        );
        assert!(home_line(&sec).contains(NONE_USED), "{af:?}");
    }
}
