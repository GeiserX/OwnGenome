//! Step 45's counts are on FamilyTreeDNA's scale where a published counting
//! rule exists.
//!
//! A HipSTR VCF that calls every region of `assets/ystr_panel.tsv` at the
//! reference length (GB 0, every read at that length) is the GRCh38
//! reference itself, so `og y str-report` must give the reference's
//! FamilyTreeDNA counts. `REFERENCE` holds them: the GRCh38 chrY sequence of
//! each marker counted by its NIST STRBase repeat structure (the panel's
//! motif column), which FamilyTreeDNA follows for these markers, and for
//! Y-GATA-H4 ISOGG's FamilyTreeDNA-NIST conversion (one less than NIST).
//! For example DYS19 is `[TCTA]12 CCTA [TCTA]3` there: FamilyTreeDNA counts
//! the TCTA blocks, 15, while HipSTR's region holds 16 repeat units. DYS449
//! is `[TTTC]15 N50 [TTTC]14`, 29, while HipSTR's two parts hold 36 units.
//! DYS389II is `[TAGA]9 [CAGA]3 N48 [TAGA]13 [CAGA]5`, 30, while its span's
//! two parts round to 17 + 12 copies.
//!
//! The reads are 250 bp long here, so DYS389II and DYS448 (over 150 bp) are
//! counted too.

use og_reports::y::cli;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

/// Marker, FamilyTreeDNA count of the GRCh38 reference.
const REFERENCE: [(&str, &str); 24] = [
    ("DYS393", "12"),
    ("DYS390", "24"),
    ("DYS19", "15"),
    ("DYS391", "11"),
    ("DYS389I", "12"),
    ("DYS392", "13"),
    ("DYS389II", "30"),
    ("DYS458", "16"),
    ("DYS448", "19"),
    ("DYS449", "29"),
    ("DYS460", "10"),
    ("Y-GATA-H4", "11"),
    ("DYS456", "15"),
    ("DYS576", "17"),
    ("DYS570", "17"),
    ("DYS438", "12"),
    ("DYS481", "22"),
    ("DYS505", "12"),
    ("DYS549", "13"),
    ("DYS522", "10"),
    ("DYS533", "12"),
    ("DYS635", "23"),
    ("DYS643", "11"),
    ("DYS461", "12"),
];

fn panel() -> String {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../../assets/ystr_panel.tsv")
        .display()
        .to_string()
}

/// A scratch folder, removed when the test ends.
struct Tmp(PathBuf);

impl Drop for Tmp {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

fn tmp(tag: &str) -> Tmp {
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let d = std::env::temp_dir().join(format!("og-i411-ystr-{tag}-{}-{nanos}", std::process::id()));
    std::fs::create_dir(&d).unwrap();
    Tmp(d)
}

fn og(verb: &str, args: &[String]) -> (u8, String) {
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = cli(Some(verb), args, &mut out, &mut err);
    (code, String::from_utf8(err).unwrap())
}

/// Every region of the panel called with this GB, from 20 reads that all
/// show it, then the table by marker: the count, the confidence and the note.
fn reference_table(d: &Path, gb: &str) -> BTreeMap<String, [String; 3]> {
    let s = |p: &Path| p.display().to_string();
    let bed = d.join("regions.bed");
    let (code, err) = og(
        "str-bed",
        &["--panel".into(), panel(), "--out".into(), s(&bed)],
    );
    assert_eq!(code, 0, "{err}");
    let mut vcf = String::from(
        "##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n",
    );
    for l in std::fs::read_to_string(&bed).unwrap().lines() {
        let c: Vec<&str> = l.split('\t').collect();
        let pos: i64 = c[1].parse().unwrap();
        vcf += &format!(
            "{}\t{}\t{}\tA\t.\t.\tPASS\tPERIOD={}\tGT:GB:Q:DP:DSTUTTER:DFLANKINDEL:MALLREADS\t0:{gb}:0.99:20:0:0:{gb}|20\n",
            c[0],
            pos - 1,
            c[5],
            c[3]
        );
    }
    let v = d.join("calls.vcf");
    std::fs::write(&v, vcf).unwrap();
    let (tsv, csv) = (d.join("o.tsv"), d.join("o.csv"));
    let (code, err) = og(
        "str-report",
        &[
            "--panel".into(),
            panel(),
            "--vcf".into(),
            s(&v),
            "--tsv".into(),
            s(&tsv),
            "--csv".into(),
            s(&csv),
            "--read-length".into(),
            "250".into(),
        ],
    );
    assert_eq!(code, 0, "{err}");
    let text = std::fs::read_to_string(&tsv).unwrap();
    text.lines()
        .skip(1)
        .map(|l| {
            let c: Vec<&str> = l.split('\t').collect();
            (c[2].to_owned(), [c[3].into(), c[4].into(), c[9].into()])
        })
        .collect()
}

/// The panel rows: marker and offset cell.
fn offsets() -> Vec<(String, String)> {
    std::fs::read_to_string(panel())
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('#') && !l.starts_with("order\t"))
        .map(|l| {
            let c: Vec<&str> = l.split('\t').collect();
            (c[2].to_owned(), c[6].to_owned())
        })
        .collect()
}

/// Red first: before the offsets, DYS19 read 16 and DYS449 36.
#[test]
fn the_reference_reads_familytreedna_s_counts() {
    let d = tmp("ref");
    let t = reference_table(&d.0, "0");
    for (m, want) in REFERENCE {
        assert_eq!(
            [t[m][0].as_str(), t[m][1].as_str()],
            [want, "high"],
            "{m}: the GRCh38 reference on FamilyTreeDNA's scale"
        );
    }
    for (m, _) in REFERENCE {
        assert!(
            !t[m][2].contains(og_reports::ystr::REFERENCE_SCALE),
            "{m}: a count on FamilyTreeDNA's scale does not say otherwise"
        );
    }
}

/// An offset no published count backs is a guess: every marker with one is
/// in `REFERENCE`, and every called marker without one says that its count
/// is on the reference's scale.
#[test]
fn every_offset_is_checked_and_every_other_count_says_its_scale() {
    let checked: Vec<&str> = REFERENCE.iter().map(|(m, _)| *m).collect();
    let rows = offsets();
    let with: Vec<&str> = rows
        .iter()
        .filter(|(_, o)| o != ".")
        .map(|(m, _)| m.as_str())
        .collect();
    let mut a = with.clone();
    let mut b = checked.clone();
    a.sort_unstable();
    b.sort_unstable();
    assert_eq!(a, b, "the markers with an offset are the checked ones");

    let d = tmp("scale");
    let t = reference_table(&d.0, "0");
    for (m, o) in &rows {
        let [count, conf, note] = &t[m];
        if conf == "not callable" {
            continue;
        }
        assert!(!count.is_empty(), "{m}: called");
        assert_eq!(
            o == ".",
            note.contains(og_reports::ystr::REFERENCE_SCALE),
            "{m}: offset {o}, note {note}"
        );
    }
}

/// A change in the reads moves the count from the reference's FamilyTreeDNA
/// count: one repeat more at every 4-base marker.
#[test]
fn a_change_counts_from_the_familytreedna_reference() {
    let d = tmp("plus");
    let t = reference_table(&d.0, "4");
    assert_eq!(t["DYS19"][0], "16");
    assert_eq!(t["DYS393"][0], "13");
    assert_eq!(t["Y-GATA-H4"][0], "12");
    // DYS449's two parts each gain a repeat: 29 + 2.
    assert_eq!(t["DYS449"][0], "31");
}
