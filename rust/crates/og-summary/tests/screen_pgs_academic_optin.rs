//! `og summary collect` lists, under the PRS section's `skipped`, the scores
//! `assets/pgs_scores.tsv` marks academic-research that the run's PRS table
//! does not hold, with the reason, so a reader sees why four traits are
//! missing. Nothing changes for a run that scored them, or without the list.

use pycompat::json::{Value, loads};
use std::fs;
use std::path::{Path, PathBuf};

const LIST: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../../assets/pgs_scores.tsv");
const REASON: &str = "academic-only licence (PRS_ACADEMIC=true to include)";

/// The ids of the shipped list and their licence column, in its order.
fn listed() -> Vec<(String, String, String)> {
    let text = fs::read_to_string(LIST).expect("list");
    let mut head: Vec<String> = Vec::new();
    let mut out = Vec::new();
    for line in text.lines().filter(|l| !l.starts_with('#')) {
        let c: Vec<&str> = line.split('\t').collect();
        if head.is_empty() {
            head = c.iter().map(|s| (*s).to_owned()).collect();
            continue;
        }
        let at = |k: &str| c[head.iter().position(|h| h == k).expect(k)].to_owned();
        out.push((at("pgs_id"), at("trait_reported"), at("licence")));
    }
    out
}

struct Case(PathBuf);

impl Case {
    /// A sample folder whose PRS table holds `ids`; `with_list` stages the
    /// shipped list in the assets folder.
    fn new(tag: &str, ids: &[&str], with_list: bool) -> Case {
        let d = std::env::temp_dir().join(format!("og-sum-academic-{tag}-{}", std::process::id()));
        let _ = fs::remove_dir_all(&d);
        fs::create_dir_all(d.join("S/prs")).expect("dir");
        fs::create_dir_all(d.join("assets")).expect("assets");
        if with_list {
            fs::copy(LIST, d.join("assets/pgs_scores.tsv")).expect("list");
        }
        let mut t = String::from(
            "Condition\tPGS_ID\tScore_SUM\tVariants_Matched\tVariants_Total\tPercentile\tInput\n",
        );
        for id in ids {
            t += &format!("x\t{id}\t0.5\t9\t9\tNA\tgvcf\n");
        }
        fs::write(d.join("S/prs/S_prs_summary.tsv"), t).expect("table");
        Case(d)
    }

    /// The PRS section's data.
    fn prs(&self) -> Value {
        let p = |x: &Path| x.display().to_string();
        let argv: Vec<String> = [
            "--sample",
            "S",
            "--sample-dir",
            &p(&self.0.join("S")),
            "--assets",
            &p(&self.0.join("assets")),
            "--only",
            "c",
            "--out",
            &p(&self.0.join("c.json")),
        ]
        .iter()
        .map(|s| (*s).to_owned())
        .collect();
        let (mut out, mut err) = (Vec::new(), Vec::new());
        let rc = og_summary::collect::cli(&argv, &mut out, &mut err);
        assert_eq!(rc, 0, "{}", String::from_utf8_lossy(&err));
        let j = loads(&fs::read_to_string(self.0.join("c.json")).expect("json")).expect("parse");
        let mut v = &j;
        for k in ["sections", "prs", "data"] {
            v = match v {
                Value::Object(o) => o.get(k).unwrap_or_else(|| panic!("no {k} in {j:?}")),
                _ => panic!("{k}: not an object"),
            };
        }
        v.clone()
    }
}

impl Drop for Case {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn get<'a>(v: &'a Value, k: &str) -> Option<&'a Value> {
    match v {
        Value::Object(o) => o.get(k),
        _ => None,
    }
}

fn text(v: Option<&Value>) -> &str {
    match v {
        Some(Value::Str(s)) => s,
        _ => "",
    }
}

#[test]
fn a_default_run_lists_the_four_academic_traits_with_the_reason() {
    let all = listed();
    let author: Vec<&str> = all
        .iter()
        .filter(|r| r.2 == "author")
        .map(|r| r.0.as_str())
        .collect();
    assert_eq!(author.len(), 11);
    let c = Case::new("default", &author, true);
    let prs = c.prs();
    let Some(Value::Array(skipped)) = get(&prs, "skipped") else {
        panic!("no skipped list: {prs:?}");
    };
    let got: Vec<(String, String, String)> = skipped
        .iter()
        .map(|s| {
            (
                text(get(s, "pgs_id")).to_owned(),
                text(get(s, "condition")).to_owned(),
                text(get(s, "reason")).to_owned(),
            )
        })
        .collect();
    let want: Vec<(String, String, String)> = all
        .iter()
        .filter(|r| r.2 == "academic-research")
        .map(|r| (r.0.clone(), r.1.clone(), REASON.to_owned()))
        .collect();
    assert_eq!(got, want);
    assert_eq!(got[0].1, "Type 2 diabetes (T2D)");
    let Some(Value::Array(scores)) = get(&prs, "scores") else {
        panic!("no scores")
    };
    assert_eq!(scores.len(), 11);
}

#[test]
fn a_run_that_scored_them_or_has_no_list_says_nothing() {
    let all: Vec<String> = listed().into_iter().map(|r| r.0).collect();
    let ids: Vec<&str> = all.iter().map(String::as_str).collect();
    assert!(get(&Case::new("all", &ids, true).prs(), "skipped").is_none());
    assert!(get(&Case::new("nolist", &["PGS000018"], false).prs(), "skipped").is_none());
}

#[test]
fn only_the_academic_scores_the_run_lacks_are_listed() {
    let c = Case::new("two", &["PGS000018", "PGS000014", "PGS000017"], true);
    let Some(Value::Array(skipped)) = get(&c.prs(), "skipped").cloned() else {
        panic!("no skipped list")
    };
    let ids: Vec<&str> = skipped.iter().map(|s| text(get(s, "pgs_id"))).collect();
    assert_eq!(ids, ["PGS000016", "PGS000027"]);
}
