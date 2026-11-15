//! A heteroplasmic MT-RNR1 allele reaches PharmCAT (issue #411, pgx-04 and
//! pgx-32; floor D10).
//!
//! CPIC defines no heteroplasmy level at which aminoglycosides become safe,
//! and follows the homoplasmic guidance for any detected MT-RNR1 variant. So
//! a PASS named allele at AF 0.05 or more is the call, whatever its
//! heteroplasmy, and the log carries the level. Below 0.05 it is not called,
//! the log says it is below the detection floor, and no `Reference` is
//! passed either: the gene is not called absent.

use std::collections::HashMap;

use og_filters::mito::{Alleles, MIN_AF, Rec, decide, parse_alleles};

const INC: &str = "Increased risk of aminoglycoside-induced hearing loss";

fn alleles() -> Alleles {
    let json = r#"{"namedAlleles": [
        {"name": "Reference", "functionValue": "Normal risk of aminoglycoside-induced hearing loss"},
        {"name": "m.827A>G", "functionValue": "Normal risk of aminoglycoside-induced hearing loss"},
        {"name": "m.1494C>T", "functionValue": "Increased risk of aminoglycoside-induced hearing loss"},
        {"name": "m.1555A>G", "functionValue": "Increased risk of aminoglycoside-induced hearing loss"}]}"#;
    parse_alleles(&pycompat::json::loads(json).unwrap(), "pharmcat.jar").unwrap()
}

fn rec(pos: i128, r: &str, alt: &str, af: f64) -> Rec {
    Rec {
        pos,
        ref_: r.into(),
        alt: alt.into(),
        filter: "PASS".into(),
        af: Some(af),
    }
}

/// Two homoplasmic background records and `extra`, with the default floor
/// and depth 500 at every named position.
fn run(extra: &[Rec]) -> (Option<String>, Vec<String>) {
    let mut recs = vec![rec(263, "A", "G", 0.998), rec(750, "A", "G", 0.999)];
    recs.extend_from_slice(extra);
    let depth: HashMap<i128, i128> = (648..=1601).map(|p| (p, 500)).collect();
    let n = recs.len();
    decide(&recs, n, &alleles(), MIN_AF, Some(&depth), 20)
}

#[test]
fn the_default_floor_is_af_0_05() {
    assert!((MIN_AF - 0.05).abs() < f64::EPSILON, "MIN_AF is {MIN_AF}");
}

#[test]
fn red_first_m1555_at_af_0_40_is_the_call() {
    let (call, lines) = run(&[rec(1555, "A", "G", 0.40)]);
    assert_eq!(call.as_deref(), Some("m.1555A>G"), "{lines:?}");
    assert_eq!(
        lines[0],
        format!(
            "MT-RNR1 m.1555A>G passed to PharmCAT as an outside call (heteroplasmic, chrM:1555 A>G AF 0.400; {INC}); CPIC sets no heteroplasmy level at which aminoglycosides are safe, so the homoplasmic guidance applies"
        )
    );
}

#[test]
fn red_first_m1494_from_0_06_to_0_94_is_the_call() {
    for af in [0.05, 0.06, 0.10, 0.50, 0.94, 0.949] {
        let (call, lines) = run(&[rec(1494, "C", "T", af)]);
        assert_eq!(call.as_deref(), Some("m.1494C>T"), "AF {af}: {lines:?}");
        assert!(
            lines[0].contains("(heteroplasmic, chrM:1494 C>T AF "),
            "{lines:?}"
        );
    }
}

#[test]
fn af_0_95_and_above_still_reads_homoplasmic() {
    let (call, lines) = run(&[rec(1555, "A", "G", 0.95)]);
    assert_eq!(call.as_deref(), Some("m.1555A>G"));
    assert_eq!(
        lines[0],
        format!(
            "MT-RNR1 m.1555A>G passed to PharmCAT as an outside call (homoplasmic, chrM:1555 A>G AF 0.950; {INC})"
        )
    );
}

#[test]
fn red_first_af_0_03_is_below_the_floor_not_absent() {
    let (call, lines) = run(&[rec(1555, "A", "G", 0.03)]);
    assert_eq!(call, None, "neither the allele nor Reference: {lines:?}");
    assert_eq!(
        lines,
        [
            "held back: m.1555A>G (chrM:1555 A>G AF 0.030) is below the detection floor (AF below 0.05): not called, and not absent either",
            "no MT-RNR1 call passed to PharmCAT: it stays No Result; read the records above in the chrM call set",
        ]
    );
}

#[test]
fn red_first_a_heteroplasmic_increased_risk_allele_beats_a_homoplasmic_normal_one() {
    let (call, lines) = run(&[rec(827, "A", "G", 1.0), rec(1555, "A", "G", 0.20)]);
    assert_eq!(call.as_deref(), Some("m.1555A>G"), "{lines:?}");
    assert!(lines[0].contains("(heteroplasmic, chrM:1555 A>G AF 0.200; "));
    assert_eq!(
        lines[1],
        "also homoplasmic, not passed (PharmCAT takes one MT-RNR1 allele): m.827A>G (chrM:827 A>G AF 1.000; Normal risk of aminoglycoside-induced hearing loss)"
    );
}

#[test]
fn red_first_a_heteroplasmic_unnamed_change_at_a_named_position_holds_reference_back() {
    let (call, lines) = run(&[rec(1555, "A", "C", 0.30)]);
    assert_eq!(call, None, "{lines:?}");
    assert_eq!(
        lines[0],
        "held back: chrM:1555 A>C AF 0.300: a change at an MT-RNR1 position that PharmCAT does not name"
    );
    // Below the floor it is noise, as any other record below it.
    let (call, _) = run(&[rec(1555, "A", "C", 0.03)]);
    assert_eq!(call.as_deref(), Some("Reference"));
}

#[test]
fn red_first_a_nan_af_is_no_af_not_a_call() {
    let (call, lines) = run(&[rec(1555, "A", "G", f64::NAN)]);
    assert_eq!(call, None, "{lines:?}");
    assert_eq!(
        lines[0],
        "held back: m.1555A>G (chrM:1555 A>G AF nan) has no allele fraction"
    );
}

#[test]
fn a_named_allele_below_the_floor_does_not_hide_a_call() {
    let (call, lines) = run(&[rec(827, "A", "G", 0.03), rec(1555, "A", "G", 0.99)]);
    assert_eq!(call.as_deref(), Some("m.1555A>G"), "{lines:?}");
    assert_eq!(
        lines,
        [
            format!(
                "MT-RNR1 m.1555A>G passed to PharmCAT as an outside call (homoplasmic, chrM:1555 A>G AF 0.990; {INC})"
            ),
            "not passed: m.827A>G (chrM:827 A>G AF 0.030) is below the detection floor (AF below 0.05): not called, and not absent either".to_owned(),
        ]
    );
    // Any other held record still withholds the call.
    let mut filtered = rec(1494, "C", "T", 0.50);
    filtered.filter = "strand_bias".into();
    let (call, _) = run(&[
        rec(827, "A", "G", 0.03),
        filtered,
        rec(1555, "A", "G", 0.99),
    ]);
    assert_eq!(call, None);
}
