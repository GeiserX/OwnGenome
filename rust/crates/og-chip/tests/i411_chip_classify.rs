//! Issue #411, the CHIP classes and statuses of `og chip report`.
//!
//! Before the fix:
//! - a call under the 2% VAF that defines CHIP still read `chip`, though
//!   the step and its doc said CHIP starts at 2%;
//! - a passing call whose reads fit an inherited heterozygous variant as
//!   well as a clone (10 of 30 reads, under the 35% cut) read `chip` with
//!   nothing to say so, and the summary line counted it like any other;
//! - a `chip` row and the summary line named the result CHIP, a clinical
//!   category that needs a blood count and a doctor, with no word that it
//!   is a possible finding only;
//! - ASXL1 counted a truncating change in any exon, where the CHIP studies
//!   count the ones in its last two exons (12 and 13 of the MANE Select
//!   transcript).

use std::collections::HashMap;
use std::fs;
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};

static N: AtomicUsize = AtomicUsize::new(0);

const TABLE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../../assets/chip_genes.tsv"
);
const DNMT3A: &str = "ENST00000321117";
const ASXL1: &str = "ENST00000375687";
/// ASXL1's MANE Select exons in Ensembl 116 (chr20, + strand).
const ASXL1_EXONS: [(u64, u64); 13] = [
    (32_358_331, 32_358_832),
    (32_366_384, 32_366_466),
    (32_367_727, 32_367_729),
    (32_369_015, 32_369_123),
    (32_428_128, 32_428_248),
    (32_428_325, 32_428_422),
    (32_429_338, 32_429_431),
    (32_429_901, 32_430_053),
    (32_431_321, 32_431_484),
    (32_431_583, 32_431_679),
    (32_432_880, 32_432_985),
    (32_433_284, 32_433_917),
    (32_434_432, 32_439_319),
];

fn gff() -> String {
    let mut s = "##gff-version 3\n".to_owned();
    s += "chr2\tensembl\tgene\t25234000\t25240300\t.\t-\t.\tID=gene:DNMT3A;Name=DNMT3A\n";
    s += &format!(
        "chr2\tensembl\tmRNA\t25234000\t25240300\t.\t-\t.\tID=transcript:{DNMT3A};Parent=gene:DNMT3A;transcript_id={DNMT3A}\n"
    );
    s += &format!(
        "chr2\tensembl\texon\t25234000\t25234500\t.\t-\t.\tParent=transcript:{DNMT3A};rank=23\n###\n"
    );
    s += "chr20\tensembl\tgene\t32358331\t32439319\t.\t+\t.\tID=gene:ASXL1;Name=ASXL1\n";
    s += &format!(
        "chr20\tensembl\tmRNA\t32358331\t32439319\t.\t+\t.\tID=transcript:{ASXL1};Parent=gene:ASXL1;transcript_id={ASXL1}\n"
    );
    for (i, (a, b)) in ASXL1_EXONS.iter().enumerate() {
        s += &format!(
            "chr20\tensembl\texon\t{a}\t{b}\t.\t+\t.\tParent=transcript:{ASXL1};rank={}\n",
            i + 1
        );
    }
    s + "###\n"
}

/// One bcftools query line: a call with its reads.
fn line(
    chrom: &str,
    pos: u64,
    csq: &str,
    gene: &str,
    tx: &str,
    aa: &str,
    (ad, af, dp): (&str, &str, &str),
) -> String {
    format!(
        "{chrom}\t{pos}\tC\tT\tPASS\t{csq}|{gene}|{tx}|protein_coding|+|{aa}|{pos}C>T\t{ad}\t{af}\t{dp}\n"
    )
}

/// DNMT3A R882H with these reads.
fn r882(ad: &str, af: &str, dp: &str) -> String {
    line(
        "chr2",
        25_234_373,
        "missense",
        "DNMT3A",
        DNMT3A,
        "882R>882H",
        (ad, af, dp),
    )
}

/// `og chip report` on the calls: its stdout and its rows by column.
fn report(calls: &str) -> (String, Vec<HashMap<String, String>>) {
    let n = N.fetch_add(1, Ordering::Relaxed);
    let dir: PathBuf =
        std::env::temp_dir().join(format!("og-chip-classify-{}-{n}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    let put = |f: &str, body: &str| {
        let p = dir.join(f);
        fs::write(&p, body).unwrap();
        p.to_str().unwrap().to_owned()
    };
    let gff = put("chip.gff3", &gff());
    let bed = put(
        "chip.bed",
        "chr2\t25233990\t25234510\tDNMT3A\nchr20\t32434421\t32439319\tASXL1\n",
    );
    let calls = put("calls.tsv", calls);
    let depth = put("depth.tsv", "");
    let out = dir.join("S_chip.tsv").to_str().unwrap().to_owned();
    let dout = dir.join("S_chip_depth.tsv").to_str().unwrap().to_owned();
    let args: Vec<String> = [
        "--genes",
        TABLE,
        "--gff",
        &gff,
        "--bed",
        &bed,
        "--calls",
        &calls,
        "--depth",
        &depth,
        "--sample",
        "S",
        "--out",
        &out,
        "--depth-out",
        &dout,
    ]
    .iter()
    .map(|s| (*s).to_owned())
    .collect();
    let (mut so, mut se) = (Vec::new(), Vec::new());
    let code = og_chip::cli("report", &args, &mut so, &mut se);
    assert_eq!(code, 0, "{}", String::from_utf8_lossy(&se));
    let text = fs::read_to_string(&out).unwrap();
    let mut lines = text.lines();
    let head: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let rows = lines
        .map(|l| {
            head.iter()
                .zip(l.split('\t'))
                .map(|(k, v)| ((*k).to_owned(), v.to_owned()))
                .collect()
        })
        .collect();
    let _ = fs::remove_dir_all(&dir);
    (String::from_utf8(so).unwrap(), rows)
}

fn one(calls: &str) -> HashMap<String, String> {
    let (_, mut rows) = report(calls);
    assert_eq!(rows.len(), 1, "{rows:?}");
    rows.remove(0)
}

#[test]
fn red_first_a_call_under_2_percent_is_not_counted() {
    // 1% of 600 reads: 6 reads, enough support, but under the 2% of CHIP.
    let r = one(&r882("594,6", "0.010", "600"));
    assert_eq!(r["status"], "low_vaf", "{r:?}");
    assert_eq!(r["class"], "missense", "the class still says what it is");
    assert!(r["note"].contains("under the 2% VAF"), "{r:?}");
    // Just under and at the floor.
    assert_eq!(one(&r882("589,11", "0.019", "600"))["status"], "low_vaf");
    assert_eq!(one(&r882("588,12", "0.020", "600"))["status"], "chip");
    // The VAF from AD when AF is missing: 10 of 1000 reads is 1%.
    assert_eq!(one(&r882("990,10", ".", "1000"))["status"], "low_vaf");
    let (out, _) = report(&r882("594,6", "0.010", "600"));
    assert!(out.starts_with("S: 0 CHIP calls "), "{out}");
}

#[test]
fn red_first_a_call_that_fits_an_inherited_variant_keeps_chip_with_a_note() {
    // 10 of 30 reads: P(X <= 10) is 0.049 for a heterozygous variant.
    let r = one(&r882("20,10", "0.333", "30"));
    assert_eq!(r["status"], "chip", "no new status: {r:?}");
    assert!(
        r["note"].contains("10 of 30 reads also fit an inherited heterozygous variant"),
        "{r:?}"
    );
    // 9 of 30 reads: 0.021, still a fit.
    let r = one(&r882("21,9", "0.300", "30"));
    assert_eq!(r["status"], "chip");
    assert!(r["note"].contains("9 of 30 reads also fit"), "{r:?}");
    // 6 of 46 reads: 0.0000009, a clone; and 8 of 40: 0.00011.
    for (ad, af, dp) in [("40,6", "0.130", "46"), ("32,8", "0.200", "40")] {
        let r = one(&r882(ad, af, dp));
        assert_eq!(r["status"], "chip");
        assert!(!r["note"].contains("inherited"), "{r:?}");
    }
    // The summary line counts the calls that fit apart.
    let calls = r882("20,10", "0.333", "30")
        + &line(
            "chr2",
            25_234_374,
            "missense",
            "DNMT3A",
            DNMT3A,
            "882R>882C",
            ("40,6", "0.130", "46"),
        );
    let (out, rows) = report(&calls);
    assert_eq!(rows.len(), 2);
    assert!(
        out.starts_with("S: 2 CHIP calls in 13 genes (possible clonal haematopoiesis, not a diagnosis; 1 of them also fits an inherited variant; 0 other protein changes listed, not counted)\n"),
        "{out}"
    );
    let (out, _) = report(&r882("40,6", "0.130", "46"));
    assert!(
        !out.contains("fits an inherited"),
        "only when one fits: {out}"
    );
}

#[test]
fn red_first_a_chip_row_says_it_is_a_possible_finding_not_a_diagnosis() {
    let r = one(&r882("40,6", "0.130", "46"));
    assert_eq!(r["status"], "chip");
    assert_eq!(
        r["note"],
        "possible clonal haematopoiesis (not a diagnosis): hotspot R882"
    );
    let (out, _) = report(&r882("40,6", "0.130", "46"));
    assert!(
        out.starts_with(
            "S: 1 CHIP call in 13 genes (possible clonal haematopoiesis, not a diagnosis; "
        ),
        "{out}"
    );
}

#[test]
fn a_call_at_35_percent_or_more_says_why_it_is_left_out() {
    let r = one(&r882("20,19", "0.48", "39"));
    assert_eq!(r["status"], "likely_germline");
    assert!(
        r["note"].starts_with("VAF 0.48: at or above the 35% cut")
            && r["note"].contains("not counted"),
        "{r:?}"
    );
}

#[test]
fn red_first_asxl1_counts_truncating_changes_in_its_last_two_exons() {
    let table = fs::read_to_string(TABLE).unwrap();
    let row: Vec<&str> = table
        .lines()
        .find(|l| l.starts_with("ASXL1\t"))
        .unwrap()
        .split('\t')
        .collect();
    assert_eq!(&row[1..3], ["exons:12-13", "none"], "{row:?}");
    let asxl1 = |pos: u64, csq: &str, aa: &str| {
        one(&line(
            "chr20",
            pos,
            csq,
            "ASXL1",
            ASXL1,
            aa,
            ("40,6", "0.130", "46"),
        ))
    };
    // G646fs (c.1934dupG), the ASXL1 hotspot, in exon 13.
    let r = asxl1(32_434_638, "frameshift", "645G>645GX");
    assert_eq!(
        (r["status"].as_str(), r["exon"].as_str()),
        ("chip", "13"),
        "{r:?}"
    );
    let r = asxl1(32_433_500, "stop_gained", "500R>500*");
    assert_eq!(
        (r["status"].as_str(), r["exon"].as_str()),
        ("chip", "12"),
        "{r:?}"
    );
    // A stop in exon 5 is listed, not counted.
    let r = asxl1(32_428_200, "stop_gained", "120R>120*");
    assert_eq!(
        (r["status"].as_str(), r["exon"].as_str()),
        ("other_change", "5"),
        "{r:?}"
    );
    assert!(r["note"].contains("the table counts exons 12-13"), "{r:?}");
}

#[test]
fn the_germline_fit_is_a_binomial_tail_at_half() {
    use og_chip::report::{fits_inherited, het_lower_tail};
    let near = |a: f64, b: f64| (a - b).abs() < 1e-4;
    assert!(
        near(het_lower_tail(10, 30), 0.04937),
        "{}",
        het_lower_tail(10, 30)
    );
    assert!(
        near(het_lower_tail(9, 30), 0.02139),
        "{}",
        het_lower_tail(9, 30)
    );
    assert!(
        near(het_lower_tail(15, 30), 0.57223),
        "{}",
        het_lower_tail(15, 30)
    );
    assert_eq!(het_lower_tail(30, 30), 1.0);
    // Deep runs stay finite: 2^-5000 underflows a plain product.
    let deep = het_lower_tail(2400, 5000);
    assert!(deep.is_finite() && deep > 0.0 && deep < 0.01, "{deep}");
    assert!(
        near(het_lower_tail(2500, 5000), 0.50564),
        "{}",
        het_lower_tail(2500, 5000)
    );
    assert!(fits_inherited(10, 30) && fits_inherited(5, 20));
    assert!(!fits_inherited(6, 46) && !fits_inherited(12, 60) && !fits_inherited(0, 0));
}
