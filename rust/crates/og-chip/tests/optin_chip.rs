//! Step 49 (clonal haematopoiesis): `og chip` on synthetic bcftools output
//! and on `assets/chip_genes.tsv`. Moved from tests/test_chip_report.py
//! (its parts A and 1 to 6); the summary's reading of the two tables is in
//! og-summary's `tests/optin_summary.rs`, the report lines in og-render's
//! `tests/optin_report.rs`.
//!
//! 1. check-table passes the real table and fails a malformed token or
//!    spec, a gene listed twice, a gene that counts no class and a bad
//!    PubMed id;
//! 2. the CHIP classes of DNMT3A, TET2, JAK2, PPM1D and TP53;
//! 3. the germline exclusion, the filter and the read support;
//! 4. what is left out;
//! 5. the depth table and the summary line;
//! 6. prepare on a synthetic Ensembl GFF3, and the report command on its
//!    output.

use og_chip::binom::{detection_limit, p_at_least};
use std::collections::HashMap;
use std::fs;
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};

static N: AtomicUsize = AtomicUsize::new(0);

const TABLE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../../assets/chip_genes.tsv"
);
const OUT_COLUMNS: &str = "gene\tchrom\tpos\tref\talt\tconsequence\tprotein_change\texon\tvaf\talt_reads\tdepth\tfilter\tstatus\tclass\tnote";
const DEPTH_COLUMNS: [&str; 9] = [
    "gene",
    "transcript",
    "chrom",
    "coding_bases",
    "mean_depth",
    "detection_limit_vaf",
    "status",
    "callable_fraction",
    "assessable_fraction",
];

/// A fresh folder of the test's own.
struct Dir(PathBuf);

impl Dir {
    fn new(tag: &str) -> Dir {
        let n = N.fetch_add(1, Ordering::Relaxed);
        let d =
            std::env::temp_dir().join(format!("og-optin-chip-{}-{tag}-{n}", std::process::id()));
        let _ = fs::remove_dir_all(&d);
        fs::create_dir_all(&d).unwrap();
        Dir(d)
    }
    fn path(&self, name: &str) -> String {
        self.0.join(name).to_str().unwrap().to_owned()
    }
    fn put(&self, name: &str, text: &str) -> String {
        let p = self.path(name);
        fs::write(&p, text).unwrap();
        p
    }
}

impl Drop for Dir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

/// `og chip VERB ARGS`: exit code, stdout, stderr.
fn og(verb: &str, args: &[&str]) -> (u8, String, String) {
    let args: Vec<String> = args.iter().map(|s| (*s).to_owned()).collect();
    let (mut o, mut e) = (Vec::new(), Vec::new());
    let code = og_chip::cli(verb, &args, &mut o, &mut e);
    (
        code,
        String::from_utf8(o).unwrap(),
        String::from_utf8(e).unwrap(),
    )
}

/// The genes of the real table.
fn table_genes() -> Vec<String> {
    let text = fs::read_to_string(TABLE).unwrap();
    text.lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .skip(1)
        .map(|l| l.split('\t').next().unwrap().to_owned())
        .collect()
}

/// The rows of a TSV as maps by column.
fn rows(text: &str) -> Vec<HashMap<String, String>> {
    let mut lines = text.lines();
    let head: Vec<&str> = lines.next().unwrap_or("").split('\t').collect();
    lines
        .map(|l| {
            head.iter()
                .zip(l.split('\t'))
                .map(|(k, v)| ((*k).to_owned(), v.to_owned()))
                .collect()
        })
        .collect()
}

/// Synthetic MANE transcripts (exons as start, end, rank).
/// gene, transcript, contig, strand, exons.
type Gene = (
    &'static str,
    &'static str,
    &'static str,
    &'static str,
    &'static [(u64, u64, u32)],
);

const GENES: [Gene; 5] = [
    (
        "DNMT3A",
        "ENST00000321117",
        "chr2",
        "-",
        &[(25_234_000, 25_234_500, 23), (25_240_000, 25_240_300, 22)],
    ),
    (
        "TET2",
        "ENST00000380013",
        "chr4",
        "+",
        &[(105_200_000, 105_210_000, 3)],
    ),
    (
        "JAK2",
        "ENST00000381652",
        "chr9",
        "+",
        &[(5_073_700, 5_073_800, 14)],
    ),
    (
        "PPM1D",
        "ENST00000305921",
        "chr17",
        "+",
        &[(60_600_000, 60_600_200, 5), (60_663_000, 60_666_280, 6)],
    ),
    (
        "TP53",
        "ENST00000269305",
        "chr17",
        "-",
        &[(7_673_700, 7_674_000, 7)],
    ),
];

fn gff() -> String {
    let mut s = "##gff-version 3\n".to_owned();
    for (g, tx, c, st, ex) in GENES {
        let (first, last) = (ex[0].0, ex[ex.len() - 1].1);
        s += &format!("{c}\tensembl\tgene\t{first}\t{last}\t.\t{st}\t.\tID=gene:{g};Name={g}\n");
        s += &format!(
            "{c}\tensembl\tmRNA\t{first}\t{last}\t.\t{st}\t.\tID=transcript:{tx};Parent=gene:{g};transcript_id={tx}\n"
        );
        for (a, b, r) in ex {
            s += &format!(
                "{c}\tensembl\texon\t{a}\t{b}\t.\t{st}\t.\tParent=transcript:{tx};rank={r}\n"
            );
        }
        s += "###\n";
    }
    s
}

/// One bcftools query line of a call.
struct Call<'a> {
    chrom: &'a str,
    pos: u64,
    csq: &'a str,
    gene: &'a str,
    tx: &'a str,
    aa: &'a str,
    flt: &'a str,
    ad: &'a str,
    af: &'a str,
    dp: &'a str,
}

fn call<'a>(
    chrom: &'a str,
    pos: u64,
    csq: &'a str,
    gene: &'a str,
    tx: &'a str,
    aa: &'a str,
) -> Call<'a> {
    Call {
        chrom,
        pos,
        csq,
        gene,
        tx,
        aa,
        flt: "PASS",
        ad: "40,6",
        af: "0.13",
        dp: "46",
    }
}

impl Call<'_> {
    fn line(&self) -> String {
        format!(
            "{}\t{}\tC\tT\t{}\t{}|{}|{}|protein_coding|-|{}|{}C>T\t{}\t{}\t{}\n",
            self.chrom,
            self.pos,
            self.flt,
            self.csq,
            self.gene,
            self.tx,
            self.aa,
            self.pos,
            self.ad,
            self.af,
            self.dp
        )
    }
}

const D: &str = "missense";
const T: &str = "ENST00000321117";

fn r882<'a>() -> Call<'a> {
    call("chr2", 25_234_373, D, "DNMT3A", T, "882R>882H")
}

/// The rows `og chip report` writes for the calls, or its stderr.
fn rows_of(calls: &[Call<'_>]) -> Result<Vec<HashMap<String, String>>, String> {
    let d = Dir::new("rows");
    let gff = d.put("chip.gff3", &gff());
    let bed = d.put("chip.bed", "chr2\t25233990\t25234510\tDNMT3A\n");
    let empty = d.put("empty.tsv", "");
    let text: String = calls.iter().map(Call::line).collect();
    let calls = d.put("calls.tsv", &text);
    let (out, dout) = (d.path("S_chip.tsv"), d.path("S_chip_depth.tsv"));
    let (code, _, err) = og(
        "report",
        &[
            "--genes",
            TABLE,
            "--gff",
            &gff,
            "--bed",
            &bed,
            "--calls",
            &calls,
            "--depth",
            &empty,
            "--sample",
            "S",
            "--out",
            &out,
            "--depth-out",
            &dout,
        ],
    );
    if code != 0 {
        return Err(err);
    }
    Ok(rows(&fs::read_to_string(&out).unwrap()))
}

/// (status, class, protein change) of the one call, `None` when it is left
/// out.
fn status_of(c: &Call<'_>) -> Option<(String, String, String)> {
    let r = rows_of(std::slice::from_ref(c)).unwrap();
    r.first().map(|r| {
        (
            r["status"].clone(),
            r["class"].clone(),
            r["protein_change"].clone(),
        )
    })
}

fn status(c: &Call<'_>) -> String {
    status_of(c).map_or_else(|| "<left out>".to_owned(), |s| s.0)
}

/// The test's own oracle for the limit: the smallest VAF found 8 times in
/// 10 at `depth` with 3 reads needed, by bisection; none under 3 reads.
fn limit(depth: u64) -> Option<f64> {
    if depth < 3 {
        return None;
    }
    let (mut lo, mut hi) = (0.0f64, 1.0f64);
    for _ in 0..60 {
        let mid = (lo + hi) / 2.0;
        if p_at_least(3, depth, mid) >= 0.8 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    Some(hi)
}

#[test]
fn a_the_detection_limit_is_an_independent_binomial() {
    let (l30, l60) = (limit(30).unwrap(), limit(60).unwrap());
    assert!(0.13 < l30 && l30 < 0.16, "{l30}");
    assert!((p_at_least(3, 30, l30) - 0.8).abs() < 1e-6);
    assert!(l60 < l30);
    assert_eq!(limit(2), None);
    assert_eq!(detection_limit(2.9).unwrap(), None);
    let og30 = detection_limit(30.0).unwrap().unwrap();
    assert!((og30 - l30).abs() < 1e-6, "og {og30}, the test's own {l30}");
    assert!(table_genes().iter().any(|g| g == "DNMT3A"));
}

#[test]
fn check_table_passes_the_real_table_and_fails_each_kind_of_damage() {
    let (code, _, err) = og("check-table", &[TABLE]);
    assert_eq!((code, err.as_str()), (0, ""));
    let d = Dir::new("table");
    let head = "gene\ttruncating\tmissense\tnote\tpmids\n";
    for (what, body) in [
        (
            "a malformed missense token",
            "DNMT3A\tall\tR88X2\tn\t25426837\n",
        ),
        (
            "a malformed truncating spec",
            "DNMT3A\texons:9-2\tnone\tn\t25426837\n",
        ),
        (
            "a gene listed twice",
            "TET2\tall\tnone\tn\t25426837\nTET2\tall\tnone\tn\t25426837\n",
        ),
        (
            "a gene that counts no class",
            "TET2\tnone\tnone\tn\t25426837\n",
        ),
        ("a bad PubMed id", "TET2\tall\tnone\tn\tPMID1\n"),
    ] {
        let t = d.put("t.tsv", &format!("{head}{body}"));
        let (code, _, err) = og("check-table", &[&t]);
        assert!(code == 1 && err.contains("t.tsv: "), "{what}: {code} {err}");
    }
}

#[test]
fn the_chip_classes() {
    let s = |c: Call<'_>| status_of(&c).unwrap();
    assert_eq!(
        s(r882()),
        ("chip".into(), "missense".into(), "R882H".into())
    );
    assert_eq!(
        status(&call("chr2", 25_234_380, D, "DNMT3A", T, "880M>880L")),
        "other_change"
    );
    assert_eq!(
        s(call(
            "chr4",
            105_205_000,
            "frameshift",
            "TET2",
            "ENST00000380013",
            "1200QRS>1200P*"
        )),
        ("chip".into(), "truncating".into(), "Q1200fs".into())
    );
    assert_eq!(
        status(&call(
            "chr2",
            25_234_502,
            "splice_acceptor",
            "DNMT3A",
            T,
            ""
        )),
        "chip"
    );
    let jak2 = |aa| call("chr9", 5_073_770, D, "JAK2", "ENST00000381652", aa);
    assert_eq!(status(&jak2("617V>617F")), "chip");
    assert_eq!(status(&jak2("617V>617L")), "other_change");
    let ppm1d = |pos, aa| call("chr17", pos, "stop_gained", "PPM1D", "ENST00000305921", aa);
    assert_eq!(status(&ppm1d(60_664_000, "500R>500*")), "chip", "last exon");
    assert_eq!(
        status(&ppm1d(60_600_100, "100R>100*")),
        "other_change",
        "exon 5"
    );
    let tp53 = |aa| call("chr17", 7_673_800, D, "TP53", "ENST00000269305", aa);
    assert_eq!(status(&tp53("248R>248Q")), "chip", "DNA-binding domain");
    assert_eq!(status(&tp53("340E>340K")), "other_change");
}

#[test]
fn germline_exclusion_filter_and_read_support() {
    let with = |ad, af, dp| Call {
        ad,
        af,
        dp,
        ..r882()
    };
    assert_eq!(status(&with("20,19", "0.48", "39")), "likely_germline");
    assert_eq!(status(&with("1,35", "0.97", "36")), "likely_germline");
    assert_eq!(status(&with("25,13", "0.34", "38")), "chip");
    assert_eq!(
        status(&Call {
            flt: "weak_evidence",
            ..r882()
        }),
        "filtered"
    );
    assert_eq!(status(&with("30,2", "0.06", "32")), "low_support");
    let r = rows_of(&[with("30,10", ".", "40")]).unwrap();
    assert_eq!(r[0]["vaf"], "0.250", "the VAF falls back on AD");
}

#[test]
fn what_is_left_out() {
    assert_eq!(
        status_of(&call("chr2", 25_234_373, "synonymous", "DNMT3A", T, "882R")),
        None
    );
    assert_eq!(
        status_of(&call("chr2", 25_236_000, "intron", "DNMT3A", T, "")),
        None
    );
    assert_eq!(
        status_of(&call(
            "chr2",
            25_234_373,
            "*missense",
            "DNMT3A",
            T,
            "882R>882H"
        )),
        None,
        "downstream of a stop"
    );
    assert_eq!(
        status_of(&call(
            "chr2",
            25_234_373,
            D,
            "DNMT3A",
            "ENST00000264709",
            "838R>838H"
        )),
        None,
        "another transcript"
    );
}

#[test]
fn the_depth_table_and_the_summary_line() {
    let d = Dir::new("depth");
    let gff = d.put("chip.gff3", &gff());
    let bed = d.put("b.bed", "chr2\t100\t110\tDNMT3A\nchr4\t200\t210\tTET2\n");
    let depth: String = (101..106).map(|p| format!("chr2\t{p}\t30\n")).collect();
    let dep = d.put("d.tsv", &depth);
    let empty = d.put("empty.tsv", "");
    let (out, dout) = (d.path("S_chip.tsv"), d.path("S_chip_depth.tsv"));
    let run = |depth: &str| {
        og(
            "report",
            &[
                "--genes",
                TABLE,
                "--gff",
                &gff,
                "--bed",
                &bed,
                "--calls",
                &empty,
                "--depth",
                depth,
                "--sample",
                "S",
                "--out",
                &out,
                "--depth-out",
                &dout,
            ],
        )
    };
    let (code, stdout, err) = run(&dep);
    assert_eq!(code, 0, "{err}");
    let text = fs::read_to_string(&dout).unwrap();
    assert_eq!(text.lines().next().unwrap(), DEPTH_COLUMNS.join("\t"));
    let by: HashMap<String, HashMap<String, String>> = rows(&text)
        .into_iter()
        .map(|r| (r["gene"].clone(), r))
        .collect();
    assert_eq!(
        by["DNMT3A"]["mean_depth"], "15.0",
        "an unread base is depth 0"
    );
    assert_eq!(
        (
            by["TET2"]["status"].as_str(),
            by["TET2"]["mean_depth"].as_str()
        ),
        ("too_shallow", "0.0")
    );
    assert_eq!(by["JAK2"]["status"], "not_on_reference");
    let mut genes: Vec<&String> = by.keys().collect();
    genes.sort();
    let mut want = table_genes();
    want.sort();
    assert_eq!(
        genes,
        want.iter().collect::<Vec<_>>(),
        "a depth row per table gene"
    );
    assert_eq!(
        [
            by["DNMT3A"]["status"].as_str(),
            by["DNMT3A"]["detection_limit_vaf"].as_str(),
            by["DNMT3A"]["callable_fraction"].as_str()
        ],
        ["partly_covered", "", "0.500"]
    );
    assert!(
        stdout
            .lines()
            .last()
            .unwrap_or("")
            .contains("no result either way"),
        "{stdout}"
    );
    let (_, stdout, _) = run(&empty);
    assert!(stdout.contains("no result either way"), "{stdout}");
}

const ENSEMBL: &str = "##gff-version 3
2\tensembl\tgene\t1000\t5000\t.\t-\t.\tID=gene:G1;Name=DNMT3A;biotype=protein_coding
2\tensembl\tmRNA\t1000\t5000\t.\t-\t.\tID=transcript:T1;Parent=gene:G1;biotype=protein_coding;tag=basic;transcript_id=T1
2\tensembl\tCDS\t1100\t1200\t.\t-\t0\tID=CDS:P1;Parent=transcript:T1
2\tensembl\tmRNA\t1000\t5000\t.\t-\t.\tID=transcript:T2;Parent=gene:G1;biotype=protein_coding;tag=basic,MANE_Select;transcript_id=T2
2\tensembl\texon\t1000\t1300\t.\t-\t.\tParent=transcript:T2;rank=2
2\tensembl\tCDS\t1050\t1300\t.\t-\t0\tID=CDS:P2;Parent=transcript:T2
2\tensembl\texon\t4000\t5000\t.\t-\t.\tParent=transcript:T2;rank=1
2\tensembl\tCDS\t4000\t4100\t.\t-\t0\tID=CDS:P2;Parent=transcript:T2
9\tensembl\tgene\t100\t900\t.\t+\t.\tID=gene:G2;Name=JAK2;biotype=protein_coding
9\tensembl\tmRNA\t100\t900\t.\t+\t.\tID=transcript:T3;Parent=gene:G2;biotype=protein_coding;tag=MANE_Select;transcript_id=T3
9\tensembl\texon\t100\t900\t.\t+\t.\tParent=transcript:T3;rank=1
9\tensembl\tCDS\t200\t800\t.\t+\t0\tID=CDS:P3;Parent=transcript:T3
";

#[test]
fn prepare_then_report() {
    let d = Dir::new("prepare");
    let genes = d.put(
        "genes.tsv",
        "gene\ttruncating\tmissense\tnote\tpmids\nDNMT3A\tall\tR882\tn\t25426837\nJAK2\tnone\tV617F\tn\t15858187\n",
    );
    let gff = d.put("in.gff3", ENSEMBL);
    let fai = d.put("ref.fai", "chr2\t10000\t6\t60\t61\n");
    let (og_gff, og_bed) = (d.path("o.gff3"), d.path("o.bed"));
    let prep = |gff: &str| {
        og(
            "prepare",
            &[
                "--gff",
                gff,
                "--genes",
                &genes,
                "--fai",
                &fai,
                "--out-gff",
                &og_gff,
                "--out-bed",
                &og_bed,
            ],
        )
    };
    let (code, _, err) = prep(&gff);
    assert_eq!(code, 0, "{err}");
    let g = fs::read_to_string(&og_gff).unwrap();
    assert!(g.contains("transcript:T2") && g.contains("rank=2"), "{g}");
    assert!(!g.contains("T1"), "the other transcript is dropped: {g}");
    assert!(
        g.contains("\nchr2\t") && g.contains("\nchr9\t") && !g.contains("\n2\t"),
        "{g}"
    );
    assert_eq!(
        fs::read_to_string(&og_bed).unwrap(),
        "chr2\t1039\t1310\tDNMT3A\nchr2\t3989\t4110\tDNMT3A\n",
        "coding exons padded by 10, on the reference's contigs only"
    );
    let bad = d.put(
        "bad.gff3",
        &ENSEMBL.replace("tag=MANE_Select;", "tag=basic;"),
    );
    let (code, _, err) = prep(&bad);
    assert!(
        code == 1 && err.contains("JAK2"),
        "no MANE Select: {code} {err}"
    );
    prep(&gff);

    let calls = d.put(
        "calls.tsv",
        "chr2\t1150\tC\tT\tPASS\tmissense|DNMT3A|T2|protein_coding|-|882R>882H|1150C>T\t40,6\t0.130\t46\n",
    );
    let depth: String = (1040..1311)
        .chain(3990..4111)
        .map(|p| format!("chr2\t{p}\t30\n"))
        .collect();
    let depth = d.put("depth.tsv", &depth);
    let (out, dout) = (d.path("S_chip.tsv"), d.path("S_chip_depth.tsv"));
    let report = |calls: &str| {
        og(
            "report",
            &[
                "--genes",
                &genes,
                "--gff",
                &og_gff,
                "--bed",
                &og_bed,
                "--calls",
                calls,
                "--depth",
                &depth,
                "--sample",
                "S",
                "--out",
                &out,
                "--depth-out",
                &dout,
            ],
        )
    };
    let (code, stdout, err) = report(&calls);
    assert!(
        code == 0 && stdout.contains("1 CHIP call ") && stdout.contains("about 14%"),
        "{stdout}{err}"
    );
    let text = fs::read_to_string(&out).unwrap();
    let lines: Vec<&str> = text.lines().collect();
    assert_eq!(lines[0], OUT_COLUMNS);
    assert!(
        lines.len() == 2 && lines[1].contains("\tchip\t"),
        "{lines:?}"
    );
    let dl = fs::read_to_string(&dout).unwrap();
    let dl: Vec<&str> = dl.lines().collect();
    assert!(
        dl.len() == 3
            && dl[1].contains("\t30.0\t")
            && dl[2].split('\t').any(|c| c == "not_on_reference"),
        "{dl:?}"
    );
    let (code, _, err) = report(&d.path("none.tsv"));
    assert!(
        code == 1 && err.starts_with("ERROR: ") && !err.contains("Traceback"),
        "a missing calls file: {code} {err}"
    );
}
