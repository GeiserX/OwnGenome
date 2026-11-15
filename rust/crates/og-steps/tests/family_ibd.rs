//! `og family ibd`: the DNA segments two sequenced relatives share (step
//! 38). It replaces `tests/test_ibd_report.py`.
//!
//! Each case writes synthetic Beagle VCFs, hap-ibd segment files, PLINK
//! maps and a pedigree, runs `og family ibd`, and checks what it writes:
//!
//! 1. the output schema: the exact header of ibd_pairs.tsv and
//!    ibd_segments.tsv, one pairs row for every pair of members (a pair
//!    with no segment too), one segments row per shared stretch;
//! 2. a parent and child: hap-ibd's two segments (the parent's haplotype
//!    changes at a crossover) are one stretch over the whole analysed map,
//!    shared fraction 1, nothing in both copies, first degree;
//! 3. full siblings: both copies shared over a quarter of the map is IBD2,
//!    relatedness from (IBD1 / 4 + IBD2 / 2), first degree;
//! 4. a short segment only, in hap-ibd's other sample order, with opposite
//!    homozygotes on both sides: distant; a pair with no segment: no shared
//!    segment; a member without a pedigree row: not in the pedigree;
//! 5. two chromosomes: the analysed map is the sum of first to last marker
//!    of each, read on the genetic map;
//! 6. the report says it in words, and that one chromosome is a rough
//!    estimate;
//! 7. damaged input stops it (exit 2, an `ERROR:` line): a hap-ibd row with
//!    7 columns, a row on another chromosome, a sample the phased VCF does
//!    not hold, a map row without its bp column;
//! 8. no chromosome phased: headers only, and the report says so;
//! 9. stretches broken by phasing: two joined across a gap with no opposite
//!    homozygote, not across one with two; a stretch reaches the
//!    chromosome's first and last marker across one opposite homozygote,
//!    not across three;
//! 10. a missing input file (a map, a phased VCF, a hap-ibd file, the
//!     chromosome list) stops with one ERROR line naming it, exit 2, and
//!     writes nothing.
//!
//! The cases that read stderr run `og family` in a child process: this
//! test binary itself, with [`CHILD`] set (see [`child_entry`]).

use flate2::Compression;
use flate2::write::GzEncoder;
use std::collections::HashMap;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

const MB: i64 = 1_000_000;
const CHILD: &str = "OG_TEST_FAMILY_IBD_CHILD";

const PAIRS_HEAD: &str = "sample_a\tsample_b\tdeclared\texpected_relatedness\tsegments\tshared_cm\t\
ibd2_cm\tlongest_cm\tanalysed_cm\tshared_fraction\trelatedness\tdegree";
const SEGS_HEAD: &str = "sample_a\tsample_b\tchrom\tstart\tend\tcm\tibd2_cm";
const PEDIGREE: &str =
    "sample\tsex\tfather\tmother\tpartner\nA\tmale\tB\tX\t.\nC\tfemale\tB\tX\t.\nB\tmale\t.\t.\tX\n";

/// In the child process only: runs `og family` with the arguments in
/// [`CHILD`] (separated by U+001F) and exits with its code, so the parent
/// reads its stderr. In a normal test run it does nothing.
#[test]
fn child_entry() {
    let Ok(joined) = std::env::var(CHILD) else {
        return;
    };
    let args: Vec<String> = joined.split('\u{1f}').map(str::to_owned).collect();
    let code = og_steps::family::cli(&args);
    std::process::exit(i32::from(code));
}

/// `og family ARGS` in a child process: (exit code, stderr).
fn family_child(args: &[String]) -> (i32, String) {
    let out = Command::new(std::env::current_exe().unwrap())
        .args(["child_entry", "--exact", "--nocapture", "--test-threads=1"])
        .env(CHILD, args.join("\u{1f}"))
        .output()
        .unwrap();
    (
        out.status.code().unwrap_or(-1),
        String::from_utf8_lossy(&out.stderr).into_owned(),
    )
}

/// A hap-ibd row: sample 1, haplotype, sample 2, haplotype, chrom, start,
/// end and cM (None leaves the cM column out).
type Seg = (&'static str, u8, &'static str, u8, &'static str, i64, i64, Option<f64>);

fn scratch(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("og-family-ibd-{}-{name}", std::process::id()));
    let _ = fs::remove_dir_all(&d);
    fs::create_dir_all(&d).unwrap();
    d
}

fn gz(path: &Path, text: &str) {
    let mut e = GzEncoder::new(fs::File::create(path).unwrap(), Compression::fast());
    e.write_all(text.as_bytes()).unwrap();
    e.finish().unwrap();
}

/// 1 cM per Mb from 0 to `map_end` Mb; markers every Mb from `first` to
/// `last` Mb; `dos(sample, Mb)` is the sample's ALT count there.
fn write_chrom(
    d: &Path,
    chrom: &str,
    samples: &[&str],
    (first, last): (i64, i64),
    rows: &[Seg],
    map_end: i64,
    dos: &dyn Fn(&str, i64) -> u8,
) {
    for sub in ["maps", "phased", "ibd"] {
        fs::create_dir_all(d.join(sub)).unwrap();
    }
    let map: String = (0..=map_end)
        .map(|i| format!("{chrom}\t.\t{:.4}\t{}\n", i as f64, (i * MB).max(1)))
        .collect();
    fs::write(d.join(format!("maps/{chrom}.map")), map).unwrap();
    let mut vcf = format!(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n",
        samples.join("\t")
    );
    for i in first..=last {
        let gts: Vec<&str> = samples
            .iter()
            .map(|s| ["0|0", "0|1", "1|1"][usize::from(dos(s, i))])
            .collect();
        vcf += &format!("{chrom}\t{}\t.\tA\tG\t.\tPASS\t.\tGT\t{}\n", i * MB, gts.join("\t"));
    }
    gz(&d.join(format!("phased/{chrom}.vcf.gz")), &vcf);
    let ibd: String = rows
        .iter()
        .map(|(a, ha, b, hb, c, s, e, cm)| {
            let cm = cm.map(|x| format!("\t{x:?}")).unwrap_or_default();
            format!("{a}\t{ha}\t{b}\t{hb}\t{c}\t{s}\t{e}{cm}\n")
        })
        .collect();
    gz(&d.join(format!("ibd/{chrom}.ibd.gz")), &ibd);
}

/// The command line for `d`, after writing the chromosome list and the
/// pedigree, and removing `gone` (an input file) when given.
fn args(d: &Path, chroms: &[&str], pedigree: &str, gone: Option<&str>) -> Vec<String> {
    let list: String = chroms.iter().map(|c| format!("{c}\n")).collect();
    fs::write(d.join("chroms.txt"), list).unwrap();
    if let Some(g) = gone {
        fs::remove_file(d.join(g)).unwrap();
    }
    fs::write(d.join("pedigree.tsv"), pedigree).unwrap();
    let p = |x: &str| d.join(x).display().to_string();
    [
        "ibd",
        "--pedigree",
        &p("pedigree.tsv"),
        "--chroms",
        &p("chroms.txt"),
        "--maps",
        &p("maps"),
        "--phased",
        &p("phased"),
        "--ibd",
        &p("ibd"),
        "--out-dir",
        &p("out"),
    ]
    .map(str::to_owned)
    .to_vec()
}

fn run(d: &Path, chroms: &[&str], pedigree: &str) -> u8 {
    og_steps::family::cli(&args(d, chroms, pedigree, None))
}

type Row = HashMap<String, String>;

fn table(path: &Path) -> (String, Vec<Row>) {
    let text = fs::read_to_string(path).unwrap();
    let mut lines = text.lines();
    let head = lines.next().unwrap_or_default().to_owned();
    let keys: Vec<&str> = head.split('\t').collect();
    let rows = lines
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            keys.iter()
                .zip(l.split('\t'))
                .map(|(k, v)| ((*k).to_owned(), v.to_owned()))
                .collect()
        })
        .collect();
    (head, rows)
}

fn by_pair(rows: &[Row]) -> HashMap<(String, String), Row> {
    rows.iter()
        .map(|r| ((r["sample_a"].clone(), r["sample_b"].clone()), r.clone()))
        .collect()
}

fn get<'a>(p: &'a HashMap<(String, String), Row>, a: &str, b: &str) -> &'a Row {
    &p[&(a.to_owned(), b.to_owned())]
}

const SAMPLES: [&str; 4] = ["A", "B", "C", "D"];

/// D is 1/1 where A, B and C are 0/0, except over its short segment with B.
fn dos1(s: &str, mb: i64) -> u8 {
    match s {
        "D" if (10..=13).contains(&mb) => 1,
        "D" => 2,
        _ => 0,
    }
}

fn one(_: &str, _: i64) -> u8 {
    1
}

const CHR1: [Seg; 5] = [
    ("A", 1, "B", 1, "chr1", MB, 50 * MB, Some(49.0)), // parent-child, then the
    ("A", 1, "B", 2, "chr1", 50 * MB, 99 * MB, Some(49.0)), // parent's other haplotype
    ("A", 1, "C", 1, "chr1", MB, 99 * MB, Some(98.0)), // siblings: one copy everywhere,
    ("A", 2, "C", 2, "chr1", 20 * MB, 44_500_000, Some(24.5)), // both over a quarter
    ("D", 2, "B", 1, "chr1", 10 * MB, 13 * MB, Some(3.0)), // short, in the other order
];

fn pairs_of(rows: &[Row]) -> Vec<(String, String)> {
    rows.iter()
        .map(|r| (r["sample_a"].clone(), r["sample_b"].clone()))
        .collect()
}

fn owned(v: &[(&str, &str)]) -> Vec<(String, String)> {
    v.iter()
        .map(|(a, b)| ((*a).to_owned(), (*b).to_owned()))
        .collect()
}

#[test]
fn one_and_two_chromosomes() {
    let d = scratch("main");
    write_chrom(&d, "chr1", &SAMPLES, (1, 99), &CHR1, 100, &dos1);
    assert_eq!(run(&d, &["chr1"], PEDIGREE), 0, "1 chromosome: exit 0");
    let out = d.join("out");
    let (head, pairs) = table(&out.join("ibd_pairs.tsv"));
    assert_eq!(head, PAIRS_HEAD, "1 ibd_pairs.tsv header");
    assert_eq!(
        pairs_of(&pairs),
        owned(&[("A", "B"), ("A", "C"), ("A", "D"), ("B", "C"), ("B", "D"), ("C", "D")]),
        "1 one row for every pair of the four members"
    );
    let (shead, segs) = table(&out.join("ibd_segments.tsv"));
    assert_eq!(shead, SEGS_HEAD, "1 ibd_segments.tsv header");
    assert_eq!(
        pairs_of(&segs),
        owned(&[("A", "B"), ("A", "C"), ("B", "D")]),
        "1 one segments row per shared stretch"
    );
    let p = by_pair(&pairs);
    let ab = get(&p, "A", "B");
    let f = |r: &Row, ks: &[&str]| ks.iter().map(|k| r[*k].clone()).collect::<Vec<_>>();
    assert_eq!(ab["declared"], "parent-child");
    assert_eq!(ab["expected_relatedness"], "0.500");
    assert_eq!(ab["segments"], "1", "2 the two hap-ibd segments are one stretch");
    assert_eq!(f(ab, &["shared_cm", "analysed_cm"]), ["98.00", "98.00"]);
    assert_eq!(ab["longest_cm"], "98.00");
    assert_eq!(f(ab, &["shared_fraction", "ibd2_cm"]), ["1.0000", "0.00"]);
    assert_eq!(
        f(ab, &["relatedness", "degree"]),
        ["0.500", "first degree (parent-child)"]
    );
    let keys: Vec<&str> = SEGS_HEAD.split('\t').collect();
    assert_eq!(
        f(&segs[0], &keys),
        ["A", "B", "chr1", "1000000", "99000000", "98.00", "0.00"],
        "2 the stretch row"
    );

    let ac = get(&p, "A", "C");
    assert_eq!(ac["declared"], "full siblings");
    assert_eq!(ac["ibd2_cm"], "24.50", "3 both copies over a quarter");
    assert_eq!(ac["relatedness"], "0.625", "3 (IBD1 / 4 + IBD2 / 2) * 2");
    assert_eq!(ac["degree"], "first degree (full siblings)");

    let bd = get(&p, "B", "D");
    assert_eq!(f(bd, &["segments", "shared_cm"]), ["1", "3.00"]);
    assert_eq!(bd["degree"], "distant or unrelated (short segments only)");
    assert_eq!(bd["declared"], "not in the pedigree");
    assert_eq!(
        f(get(&p, "C", "D"), &["segments", "degree"]),
        ["0", "no shared segment"]
    );
    let md = fs::read_to_string(out.join("ibd_report.md")).unwrap();
    assert!(
        md.contains(
            "| A and B | parent-child | 98.0 cM | 100.0% | 98.0 cM | 1 | 0.0 cM | first degree (parent-child) |"
        ),
        "6 the report has the pair's line"
    );
    assert!(md.contains("rough estimate"), "6 one chromosome is a rough estimate");
    assert!(
        md.contains("Analysed: 1 chromosome (chr1), 98.0 cM"),
        "6 the report names what was analysed"
    );

    // 5. two chromosomes
    write_chrom(
        &d,
        "chr2",
        &SAMPLES,
        (2, 49),
        &[("B", 2, "A", 2, "chr2", 2 * MB, 49 * MB, Some(47.0))],
        50,
        &dos1,
    );
    assert_eq!(run(&d, &["chr1", "chr2"], PEDIGREE), 0, "5 two chromosomes: exit 0");
    let (_, pairs) = table(&out.join("ibd_pairs.tsv"));
    let p = by_pair(&pairs);
    let ab = get(&p, "A", "B");
    assert_eq!(ab["analysed_cm"], "145.00", "5 the analysed map is the sum (98 + 47 cM)");
    assert_eq!(
        f(ab, &["segments", "shared_cm", "longest_cm", "shared_fraction"]),
        ["2", "145.00", "98.00", "1.0000"]
    );
    fs::remove_dir_all(d).unwrap();
}

/// 9. P and Q: 1-30 and 35-60 Mb joined (nothing between), 70-99 Mb not
/// (Q is 1/1 at 62 and 66 Mb); P and R: 5-50 Mb, R 1/1 at 1-3 and 51-99
/// Mb, so no reach to the ends; P and S: 3-97 Mb reaches 1 and 99 Mb, S 1/1
/// at 99 Mb only.
fn dos9(s: &str, mb: i64) -> u8 {
    let alt = match s {
        "Q" => mb == 62 || mb == 66,
        "R" => mb <= 3 || mb >= 51,
        "S" => mb == 99,
        _ => false,
    };
    if alt { 2 } else { 0 }
}

#[test]
fn stretches_broken_by_phasing() {
    let d = scratch("stitch");
    write_chrom(
        &d,
        "chr3",
        &["P", "Q", "R", "S"],
        (1, 99),
        &[
            ("P", 1, "Q", 1, "chr3", MB, 30 * MB, Some(29.0)),
            ("P", 1, "Q", 2, "chr3", 35 * MB, 60 * MB, Some(25.0)),
            ("P", 1, "Q", 1, "chr3", 70 * MB, 99 * MB, Some(29.0)),
            ("P", 2, "R", 1, "chr3", 5 * MB, 50 * MB, Some(45.0)),
            ("P", 2, "S", 2, "chr3", 3 * MB, 97 * MB, Some(94.0)),
        ],
        100,
        &dos9,
    );
    let ped = "sample\tsex\tfather\tmother\tpartner\nP\tmale\t.\t.\t.\n";
    assert_eq!(run(&d, &["chr3"], ped), 0, "9 stitching: exit 0");
    let (_, segs) = table(&d.join("out/ibd_segments.tsv"));
    let got: Vec<[String; 4]> = segs
        .iter()
        .map(|x| ["sample_a", "sample_b", "start", "end"].map(|k| x[k].clone()))
        .collect();
    let want: Vec<[String; 4]> = [
        ["P", "Q", "1000000", "60000000"],
        ["P", "Q", "70000000", "99000000"],
        ["P", "R", "5000000", "50000000"],
        ["P", "S", "1000000", "99000000"],
    ]
    .iter()
    .map(|r| r.map(str::to_owned))
    .collect();
    assert_eq!(got, want, "9 stitching: the stretches");
    let (_, pairs) = table(&d.join("out/ibd_pairs.tsv"));
    let p = by_pair(&pairs);
    let pq = get(&p, "P", "Q");
    assert_eq!(
        ["segments", "shared_cm", "longest_cm"].map(|k| pq[k].as_str()),
        ["2", "88.00", "59.00"],
        "9 P and Q, two stretches, 88 cM, longest 59"
    );
    fs::remove_dir_all(d).unwrap();
}

#[test]
fn damaged_input_exits_2_and_says_why() {
    let cases: [(&str, Vec<Seg>, Option<&str>); 4] = [
        (
            "a hap-ibd row with 7 columns",
            vec![("A", 1, "B", 1, "chr1", MB, 5 * MB, None)],
            None,
        ),
        (
            "a row on another chromosome",
            vec![("A", 1, "B", 1, "chr9", MB, 5 * MB, Some(4.0))],
            None,
        ),
        (
            "a sample the phased VCF does not hold",
            vec![("A", 1, "Z", 1, "chr1", MB, 5 * MB, Some(4.0))],
            None,
        ),
        (
            "a map row without its bp column",
            CHR1.to_vec(),
            Some("chr1\t.\t3.0\n"),
        ),
    ];
    for (i, (name, rows, mapline)) in cases.into_iter().enumerate() {
        let d = scratch(&format!("bad{i}"));
        write_chrom(&d, "chr1", &SAMPLES, (1, 99), &rows, 100, &one);
        if let Some(line) = mapline {
            let p = d.join("maps/chr1.map");
            let mut text = fs::read_to_string(&p).unwrap();
            text += line;
            fs::write(&p, text).unwrap();
        }
        let (code, err) = family_child(&args(&d, &["chr1"], PEDIGREE, None));
        assert_eq!(code, 2, "7 {name}: exit 2 ({err})");
        assert!(err.starts_with("ERROR: "), "7 {name}: says why: {err:?}");
        fs::remove_dir_all(d).unwrap();
    }
}

#[test]
fn no_chromosome_phased() {
    let d = scratch("none");
    assert_eq!(run(&d, &[], PEDIGREE), 0, "8 no chromosome: exit 0");
    let out = d.join("out");
    assert_eq!(
        fs::read_to_string(out.join("ibd_pairs.tsv")).unwrap(),
        format!("{PAIRS_HEAD}\n")
    );
    assert_eq!(
        fs::read_to_string(out.join("ibd_segments.tsv")).unwrap(),
        format!("{SEGS_HEAD}\n")
    );
    assert!(
        fs::read_to_string(out.join("ibd_report.md"))
            .unwrap()
            .contains("No chromosome was phased")
    );
    fs::remove_dir_all(d).unwrap();
}

#[test]
fn a_missing_input_file_stops_with_one_error_line_naming_it() {
    for (i, gone) in ["maps/chr1.map", "phased/chr1.vcf.gz", "ibd/chr1.ibd.gz", "chroms.txt"]
        .into_iter()
        .enumerate()
    {
        let d = scratch(&format!("gone{i}"));
        write_chrom(&d, "chr1", &SAMPLES, (1, 99), &CHR1, 100, &one);
        let (code, err) = family_child(&args(&d, &["chr1"], PEDIGREE, Some(gone)));
        assert_eq!(code, 2, "10 missing {gone}: exit 2 ({err})");
        let want = format!("ERROR: cannot read {}: ", d.join(gone).display());
        assert!(
            err.starts_with(&want) && err.matches('\n').count() == 1,
            "10 missing {gone}: one ERROR line naming the file: {err:?}"
        );
        assert!(!d.join("out").exists(), "10 missing {gone}: nothing written");
        fs::remove_dir_all(d).unwrap();
    }
}
