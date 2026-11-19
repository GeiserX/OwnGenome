//! `og y str-bed` and `og y str-report`: Y-STR repeat counts from HipSTR's
//! calls (step 45), then the summary and both reports. It replaces
//! `tests/test_ystr_report.py`.
//!
//! Each case writes a synthetic HipSTR VCF (one sample, only the records it
//! plants), runs `og y str-report` on assets/ystr_panel.tsv and checks the
//! table:
//!
//! 1. the panel: 111 markers in FamilyTreeDNA's order (the Y-12, Y-25 and
//!    Y-37 heads checked by name), every multi-copy marker kind multicopy,
//!    the BED holds every region of the single and parts markers once, and
//!    one region over the whole locus for each span marker;
//! 2. counts: reference copies plus GB / motif length; a partial repeat as
//!    N.r; a negative GB; the sum of a marker's parts (DYS449); a span
//!    marker's parts' copies plus its one GB (DYS510, DYS389II at 250 bp);
//!    an offset from a panel that gives one;
//! 3. confidence: high; low for a posterior under 0.9, for over 15 % of
//!    reads with stutter or with a flank indel, and for under 10 reads; no
//!    call for a region without a record or without a genotype, for two
//!    lengths on the haploid Y, and for a marker whose other part is
//!    missing; a region without a record carries the reason HipSTR's log
//!    gives for it;
//! 4. not callable: multi-copy markers, markers the marker set has no
//!    region for (DYS439), and a region longer than the reads can span
//!    (DYS452 at 150 bp reads, callable at 250 bp);
//! 5. the CSV: one row per marker in the panel's order;
//! 6. a damaged input: a VCF with two samples, a panel out of order;
//! 7. the reports: og summary collect reads the table (counts per
//!    confidence) and marks one with an unknown confidence unreadable; the
//!    HTML card lists the called markers only, with the note, and shows an
//!    unreadable table as Unreadable.

use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

use og_render::model::Summary;
use og_render::{html, text};
use og_reports::y::cli;
use pycompat::json::{Dump, Value, dumps};

fn repo() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../..")
}

fn panel_path() -> String {
    repo().join("assets/ystr_panel.tsv").display().to_string()
}

/// A scratch folder that is removed when the test ends, also when it fails.
struct Tmp(PathBuf);

impl std::ops::Deref for Tmp {
    type Target = Path;
    fn deref(&self) -> &Path {
        &self.0
    }
}

impl AsRef<Path> for Tmp {
    fn as_ref(&self) -> &Path {
        &self.0
    }
}

impl AsRef<std::ffi::OsStr> for Tmp {
    fn as_ref(&self) -> &std::ffi::OsStr {
        self.0.as_os_str()
    }
}

impl Tmp {
    /// A new folder that only this test owns. `create_dir` fails when the
    /// path exists, so a run never shares or deletes another run's folder.
    fn new(kind: &str, tag: &str) -> Tmp {
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let d = std::env::temp_dir().join(format!(
            "og-ancestry-{kind}-{}-{tag}-{nanos}",
            std::process::id()
        ));
        std::fs::create_dir(&d).unwrap();
        Tmp(d)
    }
}

impl Drop for Tmp {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// A fresh folder for one test.
fn tmp(tag: &str) -> Tmp {
    Tmp::new("ystr", tag)
}

/// One region of a panel row: chrom, start, end, period, copies, name.
type Region = (String, i64, i64, i64, f64, String);

/// One row of the panel, as the test reads it.
struct Row {
    panel: String,
    marker: String,
    kind: String,
    regions: Vec<Region>,
}

/// The panel's rows; a span marker has one region over its whole locus, its
/// parts' copies rounded half up and summed.
fn read_panel(path: &str) -> Vec<Row> {
    let text = std::fs::read_to_string(path).unwrap();
    let lines: Vec<Vec<&str>> = text
        .lines()
        .filter(|l| !l.trim().is_empty() && !l.starts_with('#'))
        .map(|l| l.split('\t').collect())
        .collect();
    let head = &lines[0];
    let col = |name: &str| head.iter().position(|h| *h == name).unwrap();
    let (cp, cm, ck, cr) = (col("panel"), col("marker"), col("kind"), col("regions"));
    lines[1..]
        .iter()
        .map(|c| {
            let mut regions: Vec<Region> = if c[cr] == "." {
                Vec::new()
            } else {
                c[cr]
                    .split(';')
                    .map(|x| {
                        let p: Vec<&str> = x.splitn(4, '/').collect();
                        let (chrom, span) = p[0].split_once(':').unwrap();
                        let (s, e) = span.split_once('-').unwrap();
                        (
                            chrom.to_owned(),
                            s.parse().unwrap(),
                            e.parse().unwrap(),
                            p[1].parse().unwrap(),
                            p[2].parse().unwrap(),
                            p[3].to_owned(),
                        )
                    })
                    .collect()
            };
            if c[ck] == "span" {
                let first = regions[0].clone();
                let last = regions[regions.len() - 1].clone();
                let copies: f64 = regions.iter().map(|x| (x.4 + 0.5).trunc()).sum();
                regions = vec![(first.0, first.1, last.2, first.3, copies, c[cm].to_owned())];
            }
            Row {
                panel: c[cp].to_owned(),
                marker: c[cm].to_owned(),
                kind: c[ck].to_owned(),
                regions,
            }
        })
        .collect()
}

const FORMAT: &str = "GT:GB:Q:PQ:DP:DSNP:PSNP:PDP:GLDIFF:DSTUTTER:DFLANKINDEL";

/// One HipSTR call: GT, GB, Q, DP, DSTUTTER, DFLANKINDEL.
type Call = [&'static str; 6];

const GOOD: Call = ["0", "0", "0.99", "20", "0", "0"];

fn good(gb: &'static str) -> Call {
    [
        if gb == "0" { "0" } else { "1" },
        gb,
        "0.99",
        "20",
        "0",
        "0",
    ]
}

/// A HipSTR VCF with these records (by HipSTR region name), each written at
/// its region's start.
fn vcf(d: &Path, records: &[(&str, Call)], samples: &[&str]) -> PathBuf {
    let panel = read_panel(&panel_path());
    let mut lines = vec![
        "##fileformat=VCFv4.1".to_owned(),
        "##source=HipSTR".to_owned(),
        format!(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}",
            samples.join("\t")
        ),
    ];
    for (name, [gt, gb, q, dp, st, fl]) in records {
        let (c, s) = panel
            .iter()
            .flat_map(|r| r.regions.iter())
            .find(|x| x.5 == *name)
            .map(|x| (x.0.clone(), x.1))
            .unwrap_or_else(|| panic!("no panel region is named {name}"));
        let cell = format!("{gt}:{gb}:{q}:{q}:{dp}:0:0|0:{dp}|0:.:{st}:{fl}");
        let mut row = vec![
            c,
            (s - 1).to_string(),
            (*name).to_owned(),
            "A".into(),
            ".".into(),
            ".".into(),
            "PASS".into(),
            "PERIOD=4".into(),
            FORMAT.into(),
        ];
        row.extend(samples.iter().map(|_| cell.clone()));
        lines.push(row.join("\t"));
    }
    let n = std::fs::read_dir(d).unwrap().count();
    let p = d.join(format!("calls{n}.vcf"));
    std::fs::write(&p, lines.join("\n") + "\n").unwrap();
    p
}

/// What one `og y str-report` run wrote.
struct Out {
    by: BTreeMap<String, BTreeMap<String, String>>,
    rows: Vec<BTreeMap<String, String>>,
    csv: Vec<Vec<String>>,
    stdout: String,
}

fn og(verb: &str, args: &[String]) -> (u8, String, String) {
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = cli(Some(verb), args, &mut out, &mut err);
    (
        code,
        String::from_utf8(out).unwrap(),
        String::from_utf8(err).unwrap(),
    )
}

/// Runs `og y str-report` on these records; the exit code and stderr, and
/// the tables when it exits 0.
fn run_with(
    d: &Path,
    records: &[(&str, Call)],
    panel: &str,
    read_length: Option<u32>,
    samples: &[&str],
    log: Option<&str>,
) -> (u8, String, Option<Out>) {
    let v = vcf(d, records, samples);
    let tsv = format!("{}.tsv", v.display());
    let csv = format!("{}.csv", v.display());
    let mut args: Vec<String> = [
        "--panel",
        panel,
        "--vcf",
        v.to_str().unwrap(),
        "--tsv",
        &tsv,
        "--csv",
        &csv,
    ]
    .map(str::to_owned)
    .to_vec();
    if let Some(n) = read_length {
        args.extend(["--read-length".to_owned(), n.to_string()]);
    }
    if let Some(text) = log {
        let lp = format!("{}.log", v.display());
        std::fs::write(&lp, text).unwrap();
        args.extend(["--log".to_owned(), lp]);
    }
    let (code, stdout, stderr) = og("str-report", &args);
    if code != 0 {
        return (code, stderr, None);
    }
    let t = std::fs::read_to_string(&tsv).unwrap();
    let mut lines = t.lines();
    let head: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let rows: Vec<BTreeMap<String, String>> = lines
        .map(|l| {
            head.iter()
                .zip(l.split('\t'))
                .map(|(h, v)| ((*h).to_owned(), v.to_owned()))
                .collect()
        })
        .collect();
    let by = rows
        .iter()
        .map(|r| (r["marker"].clone(), r.clone()))
        .collect();
    let csv = std::fs::read_to_string(&csv)
        .unwrap()
        .lines()
        .map(|l| l.split(',').map(str::to_owned).collect())
        .collect();
    (
        code,
        stderr,
        Some(Out {
            by,
            rows,
            csv,
            stdout,
        }),
    )
}

fn run(d: &Path, records: &[(&str, Call)]) -> Out {
    run_at(d, records, None)
}

fn run_at(d: &Path, records: &[(&str, Call)], read_length: Option<u32>) -> Out {
    let (code, err, out) = run_with(d, records, &panel_path(), read_length, &["S1"], None);
    assert_eq!(code, 0, "{err}");
    out.unwrap()
}

fn cell<'a>(o: &'a Out, marker: &str, col: &str) -> &'a str {
    &o.by[marker][col]
}

#[test]
fn the_panel_is_familytreedna_s_111_in_order_and_the_bed_holds_each_region_once() {
    let panel = read_panel(&panel_path());
    let names: Vec<&str> = panel.iter().map(|r| r.marker.as_str()).collect();
    assert_eq!(names.len(), 111);
    assert_eq!(
        names[..12],
        [
            "DYS393", "DYS390", "DYS19", "DYS391", "DYS385a", "DYS385b", "DYS426", "DYS388",
            "DYS439", "DYS389I", "DYS392", "DYS389II"
        ]
    );
    assert_eq!(
        names[12..25],
        [
            "DYS458", "DYS459a", "DYS459b", "DYS455", "DYS454", "DYS447", "DYS437", "DYS448",
            "DYS449", "DYS464a", "DYS464b", "DYS464c", "DYS464d"
        ]
    );
    assert_eq!(
        names[25..37],
        [
            "DYS460",
            "Y-GATA-H4",
            "YCAIIa",
            "YCAIIb",
            "DYS456",
            "DYS607",
            "DYS576",
            "DYS570",
            "CDYa",
            "CDYb",
            "DYS442",
            "DYS438"
        ]
    );
    assert_eq!(names[110], "DYS435");
    let sizes: Vec<usize> = ["Y-12", "Y-25", "Y-37", "Y-67", "Y-111"]
        .iter()
        .map(|p| panel.iter().filter(|r| r.panel == *p).count())
        .collect();
    assert_eq!(sizes, [12, 13, 12, 30, 44]);
    let multi: BTreeSet<&str> = panel
        .iter()
        .filter(|r| r.kind == "multicopy")
        .map(|r| r.marker.as_str())
        .collect();
    let want: BTreeSet<&str> = [
        "DYS385a",
        "DYS385b",
        "DYS459a",
        "DYS459b",
        "DYS464a",
        "DYS464b",
        "DYS464c",
        "DYS464d",
        "YCAIIa",
        "YCAIIb",
        "CDYa",
        "CDYb",
        "DYF395S1a",
        "DYF395S1b",
        "DYS413a",
        "DYS413b",
    ]
    .into();
    assert_eq!(multi, want);

    let d = tmp("bed");
    let bed = d.join("regions.bed");
    let (code, _, err) = og(
        "str-bed",
        &[
            "--panel".into(),
            panel_path(),
            "--out".into(),
            bed.display().to_string(),
        ],
    );
    assert_eq!(code, 0, "{err}");
    let text = std::fs::read_to_string(&bed).unwrap();
    let bl: Vec<Vec<&str>> = text.lines().map(|l| l.split('\t').collect()).collect();
    // Each row's name and coordinates, with repeats kept, against the
    // panel's regions, each once.
    let mut got: Vec<(String, i64, i64, String)> = bl
        .iter()
        .map(|x| {
            (
                x[0].to_owned(),
                x[1].parse().unwrap(),
                x[2].parse().unwrap(),
                x[5].to_owned(),
            )
        })
        .collect();
    got.sort_unstable();
    let regions: Vec<(String, i64, i64, String)> = panel
        .iter()
        .flat_map(|r| {
            r.regions
                .iter()
                .map(|x| (x.0.clone(), x.1, x.2, x.5.clone()))
        })
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    assert_eq!(got, regions, "the BED has every region once, at its place");
    assert!(bl.iter().all(|x| x.len() == 6
        && 0 < x[1].parse::<i64>().unwrap()
        && x[1].parse::<i64>().unwrap() < x[2].parse::<i64>().unwrap()));
    let starts: Vec<i64> = bl.iter().map(|x| x[1].parse().unwrap()).collect();
    let mut sorted = starts.clone();
    sorted.sort_unstable();
    assert_eq!(starts, sorted, "the BED is sorted");
    let row = |name: &str| -> Vec<Vec<&str>> {
        bl.iter()
            .filter(|x| x[5].starts_with(name))
            .cloned()
            .collect()
    };
    assert_eq!(
        row("DYS393"),
        [["chrY", "3263111", "3263158", "4", "12", "DYS393"]]
    );
    let spans: BTreeSet<&str> = panel
        .iter()
        .filter(|r| r.kind == "span")
        .map(|r| r.marker.as_str())
        .collect();
    assert_eq!(
        spans,
        ["DYS389II", "DYS448", "DYS450", "DYS504", "DYS510", "DYS532"].into()
    );
    assert_eq!(
        bl.len(),
        86,
        "76 single, 6 spans, the 4 parts of DYS449 and DYS552"
    );
    assert_eq!(
        row("DYS450"),
        [["chrY", "8258259", "8258341", "5", "14", "DYS450"]],
        "one region from its first part's start to its last part's end, 9 + 5 copies"
    );
    assert_eq!(
        row("DYS510"),
        [["chrY", "15188021", "15188138", "4", "25", "DYS510"]],
        "three parts are one region, 17 + 5 + 3 copies (4.5 halves up)"
    );
    assert!(
        bl.iter()
            .all(|x| !(x[5].contains('.') && spans.contains(x[5].split('.').next().unwrap()))),
        "no part of a span marker is in the BED"
    );
    let parts: Vec<&str> = row("DYS449").iter().map(|x| x[5]).collect();
    assert_eq!(parts, ["DYS449.1", "DYS449.2"]);
}

#[test]
fn counts_are_reference_copies_plus_the_length_difference() {
    let d = tmp("counts");
    let o = run(
        &d,
        &[
            ("DYS393", good("4")),
            ("DYS390", good("-8")),
            ("DYS391", good("2")),
            ("DYS392", good("-1")),
            ("DYS449.1", good("4")),
            ("DYS449.2", good("-4")),
            ("DYS510", good("4")),
            ("DYS389II", good("4")),
            ("DYS448", good("6")),
            ("DYS19/DYS394", GOOD),
        ],
    );
    assert_eq!(cell(&o, "DYS19", "count"), "16", "a reference call");
    assert_eq!(cell(&o, "DYS393", "count"), "13", "GB +4 at a 4-base motif");
    assert_eq!(cell(&o, "DYS390", "count"), "22", "GB -8");
    assert_eq!(cell(&o, "DYS391", "count"), "11.2", "a partial repeat");
    assert_eq!(
        cell(&o, "DYS392", "count"),
        "12.2",
        "GB -1 at a 3-base motif"
    );
    assert_eq!(cell(&o, "DYS449", "count"), "36", "the sum of its parts");
    assert_eq!(
        cell(&o, "DYS510", "count"),
        "26",
        "a span's parts plus its GB"
    );
    assert_eq!(
        [
            cell(&o, "DYS510", "confidence"),
            cell(&o, "DYS510", "reads"),
            cell(&o, "DYS510", "note")
        ],
        ["high", "20", "."]
    );
    for (m, bp) in [("DYS389II", 164), ("DYS448", 156)] {
        assert_eq!(cell(&o, m, "confidence"), "not callable", "{m}");
        assert_eq!(cell(&o, m, "count"), "", "{m}");
        assert!(
            cell(&o, m, "note").contains(&format!("a {bp} bp repeat")),
            "{m}"
        );
    }
    assert!(
        o.stdout.trim().starts_with("Y-STR markers: "),
        "{}",
        o.stdout
    );
    let o = run_at(&d, &[("DYS389II", good("4"))], Some(250));
    assert_eq!(cell(&o, "DYS389II", "count"), "30");
    assert!(cell(&o, "DYS389II", "note").contains("covers DYS389I"));

    let off = d.join("offset_panel.tsv");
    let text: String = std::fs::read_to_string(panel_path())
        .unwrap()
        .lines()
        .map(|l| {
            let mut c: Vec<&str> = l.split('\t').collect();
            if c.len() == 8 && c[2] == "DYS393" {
                c[6] = "-1";
            }
            c.join("\t") + "\n"
        })
        .collect();
    std::fs::write(&off, text).unwrap();
    let (code, err, o) = run_with(
        &d,
        &[("DYS393", good("4"))],
        off.to_str().unwrap(),
        None,
        &["S1"],
        None,
    );
    assert_eq!(code, 0, "{err}");
    assert_eq!(
        cell(&o.unwrap(), "DYS393", "count"),
        "12",
        "a panel offset is added"
    );
}

#[test]
fn confidence_follows_the_posterior_the_reads_and_the_log() {
    let d = tmp("confidence");
    let o = run(
        &d,
        &[
            ("DYS393", GOOD),
            ("DYS390", ["0", "0", "0.6", "20", "0", "0"]),
            ("DYS391", ["0", "0", "0.99", "20", "4", "0"]),
            ("DYS426", ["0", "0", "0.99", "20", "0", "4"]),
            ("DYS388", ["0", "0", "0.99", "6", "0", "0"]),
            ("DYS392", [".", ".", ".", "0", "0", "0"]),
            ("DYS458", ["0/1", "0/4", "0.99", "20", "0", "0"]),
            ("DYS449.1", GOOD),
            ("DYS456", ["0", "0", "0.99", "20", "3", "3"]),
        ],
    );
    assert_eq!(
        [
            cell(&o, "DYS393", "confidence"),
            cell(&o, "DYS393", "reads"),
            cell(&o, "DYS393", "q")
        ],
        ["high", "20", "0.99"]
    );
    assert_eq!(cell(&o, "DYS390", "confidence"), "low");
    assert!(cell(&o, "DYS390", "note").contains("posterior 0.6 below 0.9"));
    assert_eq!(
        cell(&o, "DYS390", "count"),
        "24",
        "a low call keeps its count"
    );
    assert_eq!(
        [
            cell(&o, "DYS391", "confidence"),
            cell(&o, "DYS391", "stutter_reads")
        ],
        ["low", "4"]
    );
    assert_eq!(cell(&o, "DYS426", "confidence"), "low");
    assert_eq!(cell(&o, "DYS388", "confidence"), "low");
    assert!(cell(&o, "DYS388", "note").contains("6 reads, fewer than 10"));
    assert_eq!(
        cell(&o, "DYS456", "confidence"),
        "high",
        "15 % is not above"
    );
    assert_eq!(
        [
            cell(&o, "DYS392", "confidence"),
            cell(&o, "DYS392", "count")
        ],
        ["no call", ""]
    );
    assert_eq!(cell(&o, "DYS458", "confidence"), "no call");
    assert!(cell(&o, "DYS458", "note").contains("two lengths"));
    assert_eq!(cell(&o, "DYS449", "confidence"), "no call");
    assert!(cell(&o, "DYS449", "note").contains("DYS449.2"));
    assert_eq!(cell(&o, "DYS19", "confidence"), "no call");
    assert!(cell(&o, "DYS19", "note").contains("see HipSTR's log"));

    // HipSTR's log: a "Processing region" line per region (0-based start),
    // then why it stopped.
    let log = "Processing region chrY 9684379 9684443
Aligning reads to each candidate haplotype
Aborting genotyping of the locus as the sequence downstream of the repeat is too repetitive for accurate genotyping
\tFlanking sequence = AAACACTATATATATATAACACTATATATATAATAC
Processing region chrY 18041621 18041651
Skipping locus with too few reads: TOTAL=3, MIN=5
Processing region chrY 3263110 3263158
Locus timing:
Processing region chrY 8258258 8258341
Skipping locus with too few informative reads for stutter training: TOTAL=2, MIN=5
Processing region chrY 12500447 12500611
Skipping region as the reference allele length exceeds the threshold (164 vs 150)
";
    let (code, err, o) = run_with(
        &d,
        &[("DYS390", GOOD)],
        &panel_path(),
        Some(250),
        &["S1"],
        Some(log),
    );
    assert_eq!(code, 0, "{err}");
    let o = o.unwrap();
    assert_eq!(
        cell(&o, "DYS19", "note"),
        "HipSTR wrote no record for DYS19/DYS394; its log: Aborting genotyping of the locus as the \
         sequence downstream of the repeat is too repetitive for accurate genotyping"
    );
    assert_eq!(
        cell(&o, "DYS425", "note"),
        "HipSTR wrote no record for DYS425; its log: Skipping locus with too few reads: TOTAL=3, MIN=5"
    );
    assert!(cell(&o, "DYS450", "note").contains("too few informative reads for stutter training"));
    assert!(cell(&o, "DYS389II", "note").contains("exceeds the threshold (164 vs 150)"));
    assert_eq!(
        cell(&o, "DYS393", "note"),
        "HipSTR wrote no record for DYS393, and its log gives no reason"
    );
    assert_eq!(
        cell(&o, "DYS391", "note"),
        "HipSTR wrote no record for DYS391, and its log gives no reason"
    );
    assert_eq!(cell(&o, "DYS390", "note"), ".");
}

#[test]
fn multi_copy_absent_and_too_long_markers_are_not_callable() {
    let d = tmp("notcallable");
    let o = run(&d, &[("DYS452", GOOD)]);
    for m in ["DYS464a", "DYS464d", "DYF395S1a", "DYS385b", "CDYa"] {
        assert_eq!(cell(&o, m, "confidence"), "not callable", "{m}");
        assert!(cell(&o, m, "note").contains("multi-copy"), "{m}");
        assert_eq!(cell(&o, m, "count"), "", "{m}");
    }
    assert_eq!(cell(&o, "DYS439", "confidence"), "not callable");
    assert!(cell(&o, "DYS439", "note").contains("no region"));
    assert_eq!(
        [
            cell(&o, "DYS452", "confidence"),
            cell(&o, "DYS452", "count")
        ],
        ["not callable", ""],
        "155 bp at 150 bp reads"
    );
    let o = run_at(&d, &[("DYS452", GOOD)], Some(250));
    assert_eq!(
        [
            cell(&o, "DYS452", "confidence"),
            cell(&o, "DYS452", "count")
        ],
        ["high", "31"]
    );
}

#[test]
fn the_csv_and_the_tsv_follow_the_panel_and_damaged_input_stops() {
    let d = tmp("csv");
    let names: Vec<String> = read_panel(&panel_path())
        .into_iter()
        .map(|r| r.marker)
        .collect();
    let o = run(&d, &[("DYS393", good("4"))]);
    assert_eq!(o.csv[0], ["marker", "value", "confidence"]);
    let csv_names: Vec<&str> = o.csv[1..].iter().map(|r| r[0].as_str()).collect();
    assert_eq!(csv_names, names);
    assert_eq!(o.csv[1], ["DYS393", "13", "high"]);
    assert_eq!(o.csv[5], ["DYS385a", "", "not callable"]);
    let tsv_names: Vec<&str> = o.rows.iter().map(|r| r["marker"].as_str()).collect();
    assert_eq!(tsv_names, names);

    let (code, err, _) = run_with(
        &d,
        &[("DYS393", GOOD)],
        &panel_path(),
        None,
        &["S1", "S2"],
        None,
    );
    assert_eq!(code, 1);
    assert!(err.contains("2 samples"), "{err}");
    let bad = d.join("bad_panel.tsv");
    let text = std::fs::read_to_string(panel_path()).unwrap();
    let lines: Vec<&str> = text.split_inclusive('\n').collect();
    let (head, body): (Vec<&str>, Vec<&str>) = lines
        .iter()
        .partition(|x| x.starts_with('#') || x.starts_with("order\t"));
    let mut out: String = head.concat();
    out += body[1];
    out += body[0];
    out += &body[2..].concat();
    std::fs::write(&bad, out).unwrap();
    let (code, err, _) = run_with(&d, &[], bad.to_str().unwrap(), None, &["S1"], None);
    assert_eq!(code, 1);
    assert!(err.contains("order column"), "{err}");
}

/// The summary og summary collect writes for sample `s` in `dir`, as text.
fn collect(s: &str, dir: &Path) -> String {
    let input = og_summary::collect::Input {
        sample: s.to_owned(),
        sample_dir: dir.display().to_string(),
        declared_sex: None,
        assets: repo().join("assets").display().to_string(),
        group: og_summary::collect::Group::All,
        tables: repo().join("bin").display().to_string(),
    };
    let v = og_summary::collect::post::collect(&input).unwrap();
    dumps(&v, &Dump::default()).unwrap()
}

fn section(summary: &str, key: &str) -> Value {
    let v = pycompat::json::loads(summary).unwrap();
    pycompat::value::get(pycompat::value::get(&v, "sections").unwrap().unwrap(), key)
        .unwrap()
        .unwrap()
        .clone()
}

fn get<'a>(v: &'a Value, key: &str) -> &'a Value {
    pycompat::value::get(v, key).unwrap().unwrap()
}

#[test]
fn the_summary_counts_the_confidences_and_the_card_lists_the_called_markers() {
    let d = tmp("reports");
    let sd = d.join("sample");
    std::fs::create_dir_all(sd.join("y_str")).unwrap();
    let panel = read_panel(&panel_path());
    let v = vcf(
        &d,
        &[
            ("DYS393", good("4")),
            ("DYS390", ["0", "0", "0.6", "20", "0", "0"]),
        ],
        &["S1"],
    );
    let tsv = sd.join("y_str/S1_y_str.tsv");
    let (code, _, err) = og(
        "str-report",
        &[
            "--panel".into(),
            panel_path(),
            "--vcf".into(),
            v.display().to_string(),
            "--tsv".into(),
            tsv.display().to_string(),
            "--csv".into(),
            sd.join("y_str/S1_y_str.csv").display().to_string(),
        ],
    );
    assert_eq!(code, 0, "{err}");
    let summary = collect("S1", &sd);
    let sec = section(&summary, "y_str");
    assert_eq!(
        get(&sec, "source"),
        &Value::Str("y_str/S1_y_str.tsv".into())
    );
    let data = get(&sec, "data");
    let unplaced = panel
        .iter()
        .filter(|r| r.kind == "multicopy" || r.kind == "absent")
        .count() as i64;
    // Not callable: the multi-copy and absent markers, and DYS452, DYS389II
    // and DYS448, too long for 150 bp reads.
    assert_eq!(
        [
            get(data, "high"),
            get(data, "low"),
            get(data, "not_callable"),
            get(data, "no_call")
        ],
        [
            &Value::Int(1),
            &Value::Int(1),
            &Value::Int(i128::from(unplaced + 3)),
            &Value::Int(i128::from(111 - 2 - unplaced - 3))
        ]
    );
    let Value::Array(markers) = get(data, "markers") else {
        panic!("markers is not a list")
    };
    let order: Vec<&Value> = markers.iter().map(|m| get(m, "marker")).collect();
    let names: Vec<Value> = panel.iter().map(|r| Value::Str(r.marker.clone())).collect();
    assert_eq!(order, names.iter().collect::<Vec<_>>());
    let nobody = section(&collect("nobody", &d), "y_str");
    assert_eq!(get(&nobody, "state"), &Value::Str("missing".into()));

    let s = Summary::parse(&summary).unwrap();
    let page = html::report(&s).unwrap();
    assert!(page.contains("<td>DYS393</td><td>13</td><td>high</td>"));
    assert!(page.contains("<td>DYS390</td><td>24</td><td>low</td>"));
    assert!(
        !page.contains("DYS464a"),
        "markers without a count are left out"
    );
    assert!(
        page.contains("fixed number of repeats"),
        "the convention note"
    );
    let txt = text::report(&s).unwrap();
    assert!(
        txt.contains("  Y-12   DYS393        13     (high)"),
        "{txt}"
    );

    let t = std::fs::read_to_string(&tsv)
        .unwrap()
        .replacen("\thigh\t", "\tgreat\t", 1);
    std::fs::write(&tsv, t).unwrap();
    let summary = collect("S1", &sd);
    assert_eq!(
        get(&section(&summary, "y_str"), "state"),
        &Value::Str("unreadable".into()),
        "an unknown confidence is unreadable"
    );
    let page = html::report(&Summary::parse(&summary).unwrap()).unwrap();
    let card = page
        .split("<h2>Y-STR Markers (HipSTR)</h2>")
        .nth(1)
        .expect("the Y-STR card");
    let card = &card[..card.find("  </div>").unwrap()];
    assert!(card.contains("Unreadable"), "{card}");
}
