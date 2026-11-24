//! A same-person relatedness from a pair that shares few sites is only a
//! possible match, even when another pair of the run makes identity
//! `checked`. Step 33's table lists every pair at 0.9 or more in
//! `same_person_as`, as before, and names the ones on fewer than
//! `MIN_IDENTITY_SITES` shared sites (or an unknown number) again in
//! `possible_same_person_as`, a row written only when it has a name, so a
//! run without such a pair keeps its 25 rows. Synthetic somalier tables.

use og_summary::qc::{FREEMIX_WARN, Inputs, MIN_IDENTITY_SITES, table};
use std::fs;
use std::path::PathBuf;

/// somalier relate's samples.tsv with one sample, S, genotyped at 15,000
/// sites (well over the floor) and a clear male chrX.
const SAMPLES: &str = "#family_id\tsample_id\tpaternal_id\tmaternal_id\tsex\tphenotype\t\
original_pedigree_sex\tgt_depth_mean\tn_hom_ref\tn_het\tn_hom_alt\tp_middling_ab\t\
X_depth_mean\tX_n\tX_hom_ref\tX_het\tX_hom_alt\tY_depth_mean\tY_n\n\
S\tS\t-9\t-9\t1\t-9\t-9\t30.0\t5000\t6000\t4000\t0.010\t15.0\t39\t0\t0\t39\t14.0\t5\n";

const PAIRS: &str = "#sample_a\tsample_b\trelatedness\tibs0\tibs2\tn\n";

struct Dir(PathBuf);

impl Dir {
    fn new(tag: &str) -> Dir {
        let d = std::env::temp_dir().join(format!("og-thin-pair-{}-{tag}", std::process::id()));
        let _ = fs::remove_dir_all(&d);
        fs::create_dir_all(&d).unwrap();
        Dir(d)
    }
    fn put(&self, name: &str, body: &str) -> String {
        let p = self.0.join(name);
        fs::write(&p, body).unwrap();
        p.to_str().unwrap().to_owned()
    }
}

impl Drop for Dir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

/// Step 33's rows for S with these pairs.
fn rows(tag: &str, pairs: &str) -> Vec<(&'static str, String)> {
    let d = Dir::new(tag);
    let samples = d.put("s.samples.tsv", SAMPLES);
    let pairs = d.put("s.pairs.tsv", pairs);
    table(&Inputs {
        sample: "S",
        samples_tsv: &samples,
        pairs_tsv: Some(&pairs),
        freemix_warn: FREEMIX_WARN,
        ..Inputs::default()
    })
    .unwrap()
}

fn value<'a>(rows: &'a [(&'static str, String)], key: &str) -> Option<&'a str> {
    rows.iter()
        .find(|(k, _)| *k == key)
        .map(|(_, v)| v.as_str())
}

#[test]
fn the_floor_is_the_identity_floor() {
    assert_eq!(MIN_IDENTITY_SITES, 1000);
}

#[test]
fn a_thin_pair_is_named_a_possible_match_beside_a_well_compared_one() {
    // T: the same person on 12,000 shared sites. U: 0.97 on 300 shared
    // sites. V: unrelated on 15,000, which also makes identity `checked`.
    let t = rows(
        "mixed",
        &format!(
            "{PAIRS}S\tT\t0.99\t0\t300\t12000\nU\tS\t0.97\t0\t300\t300\nS\tV\t0.1\t9\t300\t15000\n"
        ),
    );
    assert_eq!(value(&t, "identity"), Some("checked"));
    assert_eq!(value(&t, "same_person_as"), Some("T (0.99), U (0.97)"));
    assert_eq!(value(&t, "possible_same_person_as"), Some("U (0.97)"));
    // Right after same_person_as, before the evidence rows.
    let at = |k: &str| t.iter().position(|(key, _)| *key == k).unwrap();
    assert_eq!(at("possible_same_person_as"), at("same_person_as") + 1);
    assert_eq!(t.len(), 26);
}

#[test]
fn a_well_compared_pair_stays_firm_and_the_table_keeps_25_rows() {
    // 1,000 shared sites is at the floor, so firm.
    let t = rows(
        "firm",
        &format!("{PAIRS}S\tT\t0.99\t0\t300\t12000\nU\tS\t0.95\t0\t300\t1000\n"),
    );
    assert_eq!(value(&t, "identity"), Some("checked"));
    assert_eq!(value(&t, "same_person_as"), Some("T (0.99), U (0.95)"));
    assert_eq!(value(&t, "possible_same_person_as"), None);
    assert_eq!(t.len(), 25);
    // No same-person pair at all: no row either.
    let t = rows("none", &format!("{PAIRS}S\tV\t0.1\t9\t300\t15000\n"));
    assert_eq!(value(&t, "same_person_as"), Some(""));
    assert_eq!(value(&t, "possible_same_person_as"), None);
    assert_eq!(t.len(), 25);
}

#[test]
fn one_site_below_the_floor_or_an_unknown_count_is_thin() {
    // 999 shared sites, while V on 1,000 makes identity checked.
    let t = rows(
        "999",
        &format!("{PAIRS}U\tS\t0.97\t0\t300\t999\nS\tV\t0.1\t9\t300\t1000\n"),
    );
    assert_eq!(value(&t, "identity"), Some("checked"));
    assert_eq!(value(&t, "possible_same_person_as"), Some("U (0.97)"));
    // An older pairs table without n: no pair is known to share enough.
    let t = rows(
        "no-n",
        "#sample_a\tsample_b\trelatedness\nS\tT\t0.95\nU\tS\t0.991\n",
    );
    assert_eq!(value(&t, "identity"), Some("uncertain"));
    assert_eq!(value(&t, "same_person_as"), Some("T (0.95), U (0.99)"));
    assert_eq!(
        value(&t, "possible_same_person_as"),
        Some("T (0.95), U (0.99)")
    );
}
