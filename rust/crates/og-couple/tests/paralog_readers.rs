//! The readers of step 35's SMN table and catalog tables: the summary
//! (`og summary collect`), the SMN1 rule beside Paraphase (`og_cn::smn`,
//! which the summary's Paraphase side-by-side and the couple screen apply),
//! and the couple screen. These are the reader checks of two Python tests
//! the port replaced:
//!
//! - tests/test_paralog_loci.py section 4: the not-assessed row step 35
//!   writes (its own lines, cut from scripts/35-paralogs.sh and run) reads
//!   as not assessed with the reason in the summary, in the SMN1 rule and in
//!   the couple screen; the old row with "." for its place loses the reason;
//!   a called table still reads one SMN1 copy as a carrier; a catalog table
//!   without agCN_filter or agCN_qual is unreadable, not "not resolved";
//! - tests/test_paraphase_report.py section 4: SMN1 from the row holding
//!   exon 7, or from its homologous region on a row in SMN2; quality under
//!   20 or a filter other than PASS is no call; short rows are skipped or
//!   not resolved, never a crash; Paraphase and Parascopy agree, disagree,
//!   or are not compared when one has no call.

mod common;

use common::{Dir, HEAD, PARASCOPY_HEAD, og, read, repo, rows};
use og_summary::collect::{Group, Input, collect};
use pycompat::json::Value;

fn at<'a>(v: &'a Value, path: &[&str]) -> &'a Value {
    path.iter().fold(v, |v, k| match v {
        Value::Object(m) => m.get(*k).unwrap_or_else(|| panic!("no {k} in {path:?}")),
        _ => panic!("{path:?}: {v:?}"),
    })
}

fn text(v: &Value) -> &str {
    match v {
        Value::Str(s) => s,
        _ => panic!("{v:?}"),
    }
}

/// The group c sections og summary collect writes for sample S in `d`.
fn summary(d: &str) -> Value {
    collect(&Input {
        sample: "S".into(),
        sample_dir: d.into(),
        declared_sex: None,
        assets: repo("assets"),
        group: Group::C,
        tables: repo("bin"),
    })
    .unwrap()
    .to_json()
}

fn smn1(path: &str) -> og_cn::smn::SmnRead {
    og_cn::smn::read(Some(path), "chr5", &[]).unwrap()
}

#[test]
fn step_35s_not_assessed_smn_table_names_its_reason_everywhere() {
    let t = Dir::new("readers");
    let d = t.path("S");
    let smn_t = t.put("S/paralogs/S_smn_copy_number.tsv", "");
    // the lines of step 35 that write a not-assessed table, run as they are
    let src = read(&repo("scripts/35-paralogs.sh"));
    let start = src
        .find("\n  chrom=${region%%:*}")
        .expect("step 35 holds the lines");
    let end = start
        + src[start..].find("\n  } > \"$out\"\n").expect("their end")
        + "\n  } > \"$out\"\n".len();
    let snippet = format!(
        "header() {{ printf 'chrom\\tstart\\tend\\tlocus\\tagCN_filter\\tagCN\\tagCN_qual\\tpsCN_filter\\tpsCN\\tpsCN_qual\\thomologous_regions\\n'; }}\n{}",
        &src[start + 1..end]
    );
    let reason = "no read over chr5:70951946-70951946";
    let st = std::process::Command::new("bash")
        .args(["-c", &snippet])
        .env("region", "chr5:70951946-70951946")
        .env("model", "SMN1.gz")
        .env("reason", reason)
        .env("out", &smn_t)
        .status()
        .unwrap();
    assert!(st.success());
    let s = summary(&d);
    let data = at(&s, &["sections", "parascopy", "data"]);
    assert_eq!(
        text(at(data, &["reading"])),
        format!("SMN1 and SMN2 not assessed: {reason}: no carrier reading")
    );
    assert_eq!(at(data, &["smn1", "cn"]), &Value::Null);
    let p = smn1(&smn_t);
    assert!(p.pscn.is_none() && p.text.contains(reason), "{p:?}");
    // the couple screen: a gap that names the reason, not a carrier
    let ped = t.put(
        "pedigree.tsv",
        "sample\tsex\tfather\tmother\tpartner\nA\tfemale\t.\t.\tB\nB\tmale\t.\t.\tA\n",
    );
    let none = t.put("none.vcf", HEAD);
    let genes = repo("assets/carrier_genes.tsv");
    let out = t.path("couple");
    let (ca, cb, sa) = (
        format!("A={none}"),
        format!("B={none}"),
        format!("A={smn_t}"),
    );
    let (rc, log) = og(
        "screen",
        &[
            "--pedigree",
            &ped,
            "--genes",
            &genes,
            "--clinvar",
            &ca,
            "--clinvar",
            &cb,
            "--smn",
            &sa,
            "--out-dir",
            &out,
        ],
    );
    assert_eq!(rc, 0, "{log}");
    assert!(rows(&format!("{out}/couple_A_B.tsv")).is_empty());
    let gap = read(&format!("{out}/couple_A_B.md"))
        .lines()
        .find(|l| l.contains("SMN1 copy number of A: "))
        .unwrap_or_default()
        .to_owned();
    assert!(
        gap.contains(reason) && gap.contains("not resolved"),
        "{gap}"
    );
    // control: the old not-assessed row, with . for its place, loses the reason there
    t.put(
        "S/paralogs/S_smn_copy_number.tsv",
        &format!("{PARASCOPY_HEAD}.\t.\t.\tSMN1\tnot assessed: {reason}\t.\t.\t.\t.\t.\t.\n"),
    );
    assert!(!smn1(&smn_t).text.contains(reason));
    t.put(
        "S/paralogs/S_smn_copy_number.tsv",
        &format!(
            "{PARASCOPY_HEAD}chr5\t70925000\t70953000\tSMN1\tPASS\t4\t40\tPASS\t1,3\t35,30\tchr5:70049000-70077000:-\n"
        ),
    );
    let s = summary(&d);
    let data = at(&s, &["sections", "parascopy", "data"]);
    let reading = text(at(data, &["reading"]));
    assert_eq!(at(data, &["smn1", "cn"]), &Value::Int(1));
    assert!(
        reading.contains("carrier") && !reading.contains("not assessed"),
        "{reading}"
    );
    // a catalog table without agCN_filter or agCN_qual: unreadable, never "not resolved"
    let cols: Vec<&str> = PARASCOPY_HEAD.trim_end().split('\t').collect();
    let vals = |c: &str| match c {
        "chrom" => "chr1",
        "start" => "155231478",
        "end" => "155244699",
        "locus" => "GBA",
        "agCN_filter" | "psCN_filter" => "PASS",
        "agCN" => "4",
        "agCN_qual" => "40",
        "psCN" => "2,2",
        "psCN_qual" => "30,30",
        _ => "*",
    };
    let gba1 = |s: &Value| -> Value {
        match at(s, &["sections", "parascopy", "data", "loci"]) {
            Value::Array(a) => a
                .iter()
                .find(|x| text(at(x, &["table"])) == "gba1")
                .cloned()
                .expect("a gba1 locus"),
            v => panic!("{v:?}"),
        }
    };
    for drop in ["agCN_filter", "agCN_qual"] {
        let kept: Vec<&str> = cols.iter().copied().filter(|c| *c != drop).collect();
        let line: Vec<&str> = kept.iter().map(|c| vals(c)).collect();
        t.put(
            "S/paralogs/S_gba1_copy_number.tsv",
            &format!("{}\n{}\n", kept.join("\t"), line.join("\t")),
        );
        let loc = gba1(&summary(&d));
        assert_eq!(
            text(at(&loc, &["status"])),
            "unreadable",
            "without {drop}: {loc:?}"
        );
        assert_eq!(at(&loc, &["cn"]), &Value::Null);
    }
    t.put(
        "S/paralogs/S_gba1_copy_number.tsv",
        &format!(
            "{PARASCOPY_HEAD}chr1\t155231478\t155244699\tGBA\tPASS\t4\t40\tPASS\t2,2\t30,30\t*\n"
        ),
    );
    let loc = gba1(&summary(&d));
    assert_eq!(text(at(&loc, &["status"])), "called", "control: {loc:?}");
    assert_eq!(at(&loc, &["cn"]), &Value::Int(4));
}

#[test]
fn smn1_from_parascopy_and_beside_paraphase() {
    let t = Dir::new("smn1");
    let p = t.put(
        "s1.tsv",
        &format!("{PARASCOPY_HEAD}chr5\t70925000\t70953000\tSMN1\tPASS\t4\t40\tPASS\t2,2\t35,30\tchr5:70049000-70077000:-\n"),
    );
    let s1 = smn1(&p);
    assert!(
        s1.pscn == Some(2) && s1.text.contains("quality 35"),
        "the row in SMN1, first psCN: {s1:?}"
    );
    let p2 = t.put(
        "s2.tsv",
        &format!("{PARASCOPY_HEAD}chr5\t70049000\t70077000\tSMN1\tPASS\t4\t40\tPASS\t3,1\t30,25\tchr5:70925000-70953000:-\n"),
    );
    assert_eq!(
        smn1(&p2).pscn,
        Some(1),
        "the homologous region of a row in SMN2, second psCN"
    );
    let p3 = t.put(
        "s3.tsv",
        &format!(
            "{PARASCOPY_HEAD}chr5\t70925000\t70953000\tSMN1\tPASS\t4\t40\tPASS\t2,2\t12,30\t*\n"
        ),
    );
    let s3 = smn1(&p3);
    assert!(
        s3.pscn.is_none() && s3.text.contains("not resolved"),
        "quality 12 is no call: {s3:?}"
    );
    let p4 = t.put(
        "s4.tsv",
        &format!("{PARASCOPY_HEAD}chr5\t70925000\t70953000\tSMN1\tPASS\t4\t40\tFewReliable\t2,2\t35,30\t*\n"),
    );
    assert_eq!(smn1(&p4).pscn, None, "a filter other than PASS is no call");
    // A short row (fewer fields than the header, as a cut-off file leaves) is
    // skipped, and the next row still read; one that reaches exon 7 is not resolved.
    let p5 = t.put(
        "s5.tsv",
        &format!(
            "{PARASCOPY_HEAD}chr5\nchr5\t70925000\nchr5\t70925000\t70953000\tSMN1\tPASS\t4\t40\tPASS\t2,2\t35,30\t*\n"
        ),
    );
    assert_eq!(smn1(&p5).pscn, Some(2));
    let p6 = t.put(
        "s6.tsv",
        &format!("{PARASCOPY_HEAD}chr5\t70925000\t70953000\n"),
    );
    let s6 = smn1(&p6);
    assert!(
        s6.pscn.is_none() && s6.text.contains("not resolved"),
        "{s6:?}"
    );
    // agree, disagree, or no comparison: the summary's Paraphase section
    let d = t.path("S");
    let compare = |paraphase: &str, smn_table: &str| -> Value {
        t.put(
            "S/paralogs/S_paraphase_copy_number.tsv",
            &format!(
                "region\tgene\tcopy_number\tstatus\tdetail\ttotal_cn\tregion_depth\tgenome_depth\n\
                 smn1\tSMN1\t{paraphase}\tcalled\t.\t4\t120\t31.5\n"
            ),
        );
        t.put("S/paralogs/S_smn_copy_number.tsv", &read(smn_table));
        at(
            &summary(&d),
            &["sections", "paralogs", "data", "smn1", "agree"],
        )
        .clone()
    };
    assert_eq!(compare("2", &p), Value::Str("agree".into()));
    assert_eq!(compare("2", &p2), Value::Str("disagree".into()));
    assert_eq!(compare(".", &p), Value::Null);
    assert_eq!(compare("2", &p3), Value::Null);
}
