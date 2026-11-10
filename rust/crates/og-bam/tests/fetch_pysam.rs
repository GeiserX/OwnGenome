//! The fetch smoke: og-bam on a slice of the public test fixture writes,
//! byte for byte, what pysam writes (`tests/data/fetch.pysam.tsv`, made by
//! `tests/pysam_fetch.py`; `tests/data/README.md` says how).

use std::fs::File;
use std::io::BufReader;
use std::path::Path;

#[test]
fn fetch_equals_pysam_on_the_fixture_slice() {
    let data = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let regions = BufReader::new(File::open(data.join("regions.tsv")).unwrap());
    let mut got = Vec::new();
    og_bam::fetch_probe(&data.join("cyp2d6.bam"), regions, &mut got).unwrap();
    let want = std::fs::read_to_string(data.join("fetch.pysam.tsv")).unwrap();
    let got = String::from_utf8(got).unwrap();
    assert!(!want.is_empty(), "the pysam file is empty");
    for (i, (g, w)) in got.lines().zip(want.lines()).enumerate() {
        assert_eq!(g, w, "line {} differs from pysam", i + 1);
    }
    assert_eq!(
        got.lines().count(),
        want.lines().count(),
        "line count differs from pysam"
    );
}

#[test]
fn an_unknown_contig_is_an_error() {
    let data = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let r = og_bam::fetch_probe(
        &data.join("cyp2d6.bam"),
        "chrNope\t0\t10\n".as_bytes(),
        Vec::new(),
    );
    assert!(r.is_err());
}
