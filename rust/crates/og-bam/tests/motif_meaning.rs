//! What `og str motif`'s goldens must mean, by meaning: the checks of the
//! old `tests/test_rfc1_motif.py`, in cargo test.
//!
//! `tests/parity.rs` holds og to its goldens byte for byte
//! (`tests/data/parity/want/<case>`, or `og/<case>` for a deliberate fix).
//! This test says what those goldens must say about the made-up reads at
//! the GRCh38 RFC1 repeat, so a golden rewritten with another verdict fails
//! here even when og agrees with it:
//!
//! - benign: AAAAG reads across the repeat and soft-clipped, benign_motif;
//! - pathogenic: AAGGG reads soft-clipped, unmapped and placed at their
//!   mate, rescued from another contig and on the reverse strand, plus three
//!   expanded AAAAG reads; the eight reads across an unexpanded (AAAAG)11
//!   allele are not counted: pathogenic_motif with both shares;
//! - mixed (AAAAG and AAAGG): benign_motif with both shares; acagg:
//!   pathogenic_motif; unknown (AAGAG): uncertain_motif; few (two reads):
//!   not_resolved;
//! - controls, each failing when the rescue it names is off: `--no-softclip`
//!   counts the 6 clipped reads fewer, `--no-mates` the 3 rescued mates;
//! - `--skip` and a VCF without an RFC1 record: not_checked.

use std::collections::HashMap;
use std::path::{Path, PathBuf};

fn parity() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/parity")
}

/// og's golden of a parity case when it has one (a deliberate fix), else
/// the oracle's.
fn golden(case: &str) -> PathBuf {
    let og = parity().join("og").join(case);
    if og.is_dir() {
        og
    } else {
        parity().join("want").join(case)
    }
}

/// The one row of the case's motif.tsv, by column.
fn table(case: &str) -> HashMap<String, String> {
    let path = golden(case).join("files/motif.tsv");
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    let mut lines = text.lines();
    let head = lines.next().unwrap_or_default().split('\t');
    let row = lines.next().unwrap_or_default().split('\t');
    head.zip(row)
        .map(|(k, v)| (k.to_owned(), v.to_owned()))
        .collect()
}

fn exit_code(case: &str) -> i32 {
    let path = golden(case).join("exit");
    std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("{}: {e}", path.display()))
        .trim()
        .parse()
        .unwrap()
}

fn get<'a>(t: &'a HashMap<String, String>, k: &str) -> &'a str {
    t.get(k).map_or("<missing>", String::as_str)
}

#[test]
fn each_read_set_gets_its_verdict() {
    for (kind, verdict, reads) in [
        ("benign", "benign_motif", "10"),
        ("pathogenic", "pathogenic_motif", "16"),
        ("mixed", "benign_motif", "10"),
        ("acagg", "pathogenic_motif", "5"),
        ("unknown", "uncertain_motif", "5"),
        ("few", "not_resolved", "2"),
    ] {
        let case = format!("motif-{kind}");
        let t = table(&case);
        assert_eq!(
            (
                exit_code(&case),
                get(&t, "verdict"),
                get(&t, "repeat_reads")
            ),
            (0, verdict, reads),
            "{kind}: exit 0, {verdict}, {reads} repeat reads"
        );
    }
}

#[test]
fn the_pathogenic_reads_and_their_shares() {
    let t = table("motif-pathogenic");
    assert_eq!(
        ["chrom", "start", "end", "dominant_motif", "motif_shares"].map(|k| get(&t, k)),
        [
            "chr4",
            "39348425",
            "39348479",
            "AAGGG",
            "AAGGG=0.81,AAAAG=0.19"
        ],
        "pathogenic: the locus from the VCF, AAGGG dominant, both shares; the 8 reference-allele reads not counted"
    );
    assert!(
        get(&t, "note").contains("repeat-primed PCR"),
        "pathogenic: the note asks for repeat-primed PCR: {}",
        get(&t, "note")
    );
    let benign = table("motif-benign");
    assert!(
        get(&benign, "note").contains("the flag is not supported"),
        "benign: the note says the flag is not supported: {}",
        get(&benign, "note")
    );
    assert_eq!(
        get(&table("motif-mixed"), "motif_shares"),
        "AAAAG=0.60,AAAGG=0.40",
        "mixed: AAAAG and AAAGG counted apart"
    );
    let unknown = table("motif-unknown");
    assert_eq!(
        (
            get(&unknown, "dominant_motif"),
            get(&unknown, "dominant_class")
        ),
        ("AAGAG", "unknown"),
        "unknown: AAGAG named by its smallest form, class unknown"
    );
}

#[test]
fn the_rescue_controls_and_the_unchecked_cases() {
    assert_eq!(
        get(&table("motif-pathogenic-no-softclip"), "repeat_reads"),
        "10",
        "control: --no-softclip loses the 6 clipped reads"
    );
    assert_eq!(
        get(&table("motif-pathogenic-no-mates"), "repeat_reads"),
        "13",
        "control: --no-mates loses the 3 rescued mates"
    );
    let skip = table("motif-skip");
    assert_eq!(
        (get(&skip, "verdict"), get(&skip, "note")),
        ("not_checked", "no BAM"),
        "--skip: not_checked with the reason"
    );
    assert_eq!(
        get(&table("motif-no-rfc1-record"), "verdict"),
        "not_checked",
        "a VCF without an RFC1 record: not_checked"
    );
}
