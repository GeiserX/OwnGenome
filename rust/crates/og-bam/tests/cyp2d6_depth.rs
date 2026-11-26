//! The CYP2D6 depth check (`og pgx cyp2d6-bed`, `-measure`, `-check`), by
//! meaning.
//!
//! `tests/parity.rs` holds og to the files under `tests/data/parity/`, byte
//! for byte: what the Python tool og replaced wrote on each case
//! (`want/<case>/`), or what og writes where it fixes that tool on purpose
//! (`og/<case>/`). This test says what those files must mean, so a golden
//! rewritten with another verdict fails here even when og agrees with it.
//!
//! Five depth tables (`tests/fixtures/pgx/depth/<case>.q0.bed` and `.q1.bed`,
//! the mosdepth region files of all reads and of MAPQ >= 1):
//!
//! - `clean_no_alt`: the fixture's CYP2D slice mapped to the no-ALT
//!   reference: ok;
//! - `multimapped_with_alt`: the same reads mapped to the Broad hg38 FASTA
//!   with ALT contigs: unreliable, multi-mapped;
//! - `deleted_both`: CYP2D6 deleted on both copies (*5/*5), almost no reads
//!   at the gene, normal flanks: ok, the deletion is real;
//! - `deleted_one`: one copy deleted, half the flank depth: ok;
//! - `no_flank_reads`: no reads at all: unreliable.
//!
//! Also: the check file has the status, the message, the four depths and the
//! MAPQ shares; `cyp2d6-bed` prints CYP2D6 between its two flanks; an input
//! without a CYP2D6 row, with regions other than those (a zero-length one
//! included), a region listed twice with two depths, or an output that
//! cannot be written exits 2 and writes nothing. The MAPQ shares beside the
//! verdict (`cyp2d6-measure`, then `cyp2d6-check --measure`): on synthetic
//! reads at CYP2D6, reads in REP6 at MAPQ 0 lower the share over Cyrius's
//! region but not over the gene body; without the measure file every share
//! it gives reads 'not measured'; the shares never change the status; a file
//! that is not measure's exits 2.

use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};

const OK: &str = "CYP2D6 depth is fit for a copy-number call";
const UNRELIABLE: &str = "CYP2D6 copy number unreliable: reads are multi-mapped \
                          (was this BAM aligned to a reference with ALT contigs?)";
const NO_DEPTH: &str = "CYP2D6 copy number unreliable: no reads in the flanks of CYP2D6 \
                        (chr22:42.05-42.25 Mb), so its depth cannot be compared";
const ROWS: [&str; 12] = [
    "metric",
    "status",
    "message",
    "gene_depth_all",
    "gene_depth_mapq1",
    "flank_depth_all",
    "flank_depth_mapq1",
    "gene_mapq1_share",
    "gene_mapq20_share",
    "body_depth_all",
    "body_mapq1_share",
    "body_mapq20_share",
];

const CASES: [(&str, &str, &str); 5] = [
    ("clean_no_alt", "ok", OK),
    ("multimapped_with_alt", "unreliable", UNRELIABLE),
    ("deleted_both", "ok", OK),
    ("deleted_one", "ok", OK),
    ("no_flank_reads", "unreliable", NO_DEPTH),
];

fn parity() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/parity")
}

/// og's golden of a case when it has one (a deliberate fix), else the
/// oracle's.
fn golden(case: &str) -> PathBuf {
    let og = parity().join("og").join(case);
    if og.is_dir() {
        og
    } else {
        parity().join("want").join(case)
    }
}

fn exit_code(case: &str) -> i32 {
    let text = fs::read_to_string(golden(case).join("exit")).unwrap();
    text.trim().parse().unwrap()
}

/// The files the case writes, by name.
fn written(case: &str) -> Vec<String> {
    let Ok(dir) = fs::read_dir(golden(case).join("files")) else {
        return Vec::new();
    };
    let mut names: Vec<String> = dir
        .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    names.sort();
    names
}

/// A two-column table: its rows by key, and the keys in file order.
fn table(case: &str, name: &str) -> (BTreeMap<String, String>, Vec<String>) {
    let text = fs::read_to_string(golden(case).join("files").join(name)).unwrap();
    let mut rows = BTreeMap::new();
    let mut order = Vec::new();
    for line in text.lines() {
        let (k, v) = line.split_once('\t').unwrap_or((line, ""));
        rows.insert(k.to_owned(), v.to_owned());
        order.push(k.to_owned());
    }
    (rows, order)
}

fn get<'a>(rows: &'a BTreeMap<String, String>, key: &str) -> Option<&'a str> {
    rows.get(key).map(String::as_str)
}

#[test]
fn each_depth_table_gives_its_verdict_and_the_shares() {
    for (case, status, message) in CASES {
        for (suffix, measured) in [("", false), ("-measured", true)] {
            let c = format!("check-{case}{suffix}");
            assert_eq!(exit_code(&c), 0, "{c}: exit code");
            let (rows, order) = table(&c, "check.tsv");
            assert_eq!(get(&rows, "status"), Some(status), "{c}: status");
            assert_eq!(get(&rows, "message"), Some(message), "{c}: message");
            assert_eq!(order, ROWS, "{c}: header and the eleven rows");
            let shares: Vec<Option<&str>> = [
                "gene_mapq20_share",
                "body_depth_all",
                "body_mapq1_share",
                "body_mapq20_share",
            ]
            .iter()
            .map(|k| get(&rows, k))
            .collect();
            if measured {
                assert!(
                    !shares.contains(&Some("not measured")),
                    "{c}: the measure file gives every share: {shares:?}"
                );
            } else {
                assert_eq!(
                    shares,
                    [Some("not measured"); 4],
                    "{c}: without the measure file the other shares are not measured"
                );
            }
        }
    }
}

#[test]
fn clean_no_alt_depths_and_share() {
    let (rows, _) = table("check-clean_no_alt", "check.tsv");
    assert_eq!(
        get(&rows, "gene_mapq1_share"),
        Some("0.91"),
        "clean_no_alt: the verdict's MAPQ >= 1 share over Cyrius's region (16.79 / 18.53)"
    );
    let depths: Vec<Option<&str>> = [
        "gene_depth_all",
        "gene_depth_mapq1",
        "flank_depth_all",
        "flank_depth_mapq1",
    ]
    .iter()
    .map(|k| get(&rows, k))
    .collect();
    assert_eq!(
        depths,
        [Some("18.53"), Some("16.79"), Some("29.68"), Some("29.67")],
        "clean_no_alt: depths, length-weighted"
    );
}

#[test]
fn cyp2d6_bed_prints_the_gene_between_two_flanks() {
    let bed = fs::read_to_string(golden("cyp2d6-bed").join("stdout")).unwrap();
    let names: Vec<&str> = bed
        .lines()
        .map(|l| l.split('\t').nth(3).unwrap_or_default())
        .collect();
    assert_eq!(
        names,
        ["flank", "CYP2D6", "flank"],
        "cyp2d6-bed: CYP2D6 between two flanks"
    );
}

#[test]
fn bad_inputs_exit_2_and_write_nothing() {
    for (c, why) in [
        ("check-no-gene", "an input without a CYP2D6 row"),
        (
            "check-other-region",
            "regions other than the ones cyp2d6-bed prints",
        ),
        ("check-zero-length", "a zero-length region"),
        (
            "check-duplicate-row",
            "a region listed twice with two depths",
        ),
        ("check-bad-measure", "a measure file that is not measure's"),
        ("check-out-unwritable", "an output that cannot be written"),
    ] {
        assert_eq!(
            (exit_code(c), written(c)),
            (2, Vec::<String>::new()),
            "{why}: exit 2 and nothing written"
        );
    }
}

/// The synthetic reads of `parity/in/cyp2d6_synthetic.bam`: Cyrius's region tiled, MAPQ 0
/// in REP6 and around the body, the body at MAPQ 60 except some at MAPQ 10.
#[test]
fn measure_shares_on_synthetic_reads() {
    let text = fs::read_to_string(
        golden("measure-cyp2d6_synthetic")
            .join("files")
            .join("shares.tsv"),
    )
    .unwrap();
    let m: BTreeMap<(String, u32), f64> = text
        .lines()
        .skip(1)
        .map(|l| {
            let c: Vec<&str> = l.split('\t').collect();
            (
                (c[0].to_owned(), c[1].parse().unwrap()),
                c[2].parse().unwrap(),
            )
        })
        .collect();
    let at = |r: &str, q: u32| m[&(r.to_owned(), q)];
    assert_eq!(
        at("body", 1),
        at("body", 0),
        "measure: the MAPQ 0 reads outside the body leave the body's MAPQ >= 1 depth whole"
    );
    assert_eq!(
        ["gene", "body"].map(|k| at(k, 20) < at(k, 1)),
        [true, true],
        "measure: the MAPQ 10 reads count at MAPQ >= 1, not at >= 20"
    );
    assert!(
        at("gene", 1) < at("gene", 0),
        "measure: over Cyrius's region the MAPQ 0 reads lower the MAPQ >= 1 depth"
    );
}
