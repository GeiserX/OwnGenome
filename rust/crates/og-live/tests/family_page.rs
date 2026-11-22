//! `og live family` (og-live's family module) on tests/fixtures/families_live
//! and on the demo's invented trio (tests/fixtures/live/demo/DEMO-FAMILY).
//! Moved here from tests/test_live_family_depth.py and from the checks of
//! tests/test_owngenome_live_family.py that read the family plugin's data
//! (that file keeps the checks of the Python site that serves it):
//!
//! - the family folders are listed as runs of kind family, a sample's folder
//!   is not; family_a's state: its members, each child with its parents and
//!   map, the relationship check, the shared segments, the couple screen
//!   counted from step 39's table, the plan and exactly the result files
//!   that exist (a file no step writes and a link out of the folder are not
//!   among them); a members file with a header and no rows is an error;
//! - each child's stretches are its `<child>_segments.tsv` row for row, and
//!   the chromosome lengths the painting scales by are GRCh38's chr1 to
//!   chr22 and chrX; the compound-het pairs come in trans first, and every
//!   de novo row is there; a family without the trio files says so and
//!   lists nothing;
//! - the members side by side: one row per gene or system any member's
//!   summary.json has, '.' where a member's summary lacks the row, "no
//!   summary yet" for a member without one, a stale value flagged in its
//!   own cell, every blood group system step 41 typed in step 41's order,
//!   and a summary of a shape the report never writes costing its member's
//!   column only; a summary.json linking out of its run is never read;
//! - each child's trio status, and the compound-het line told by the child's
//!   step 23 candidates;
//! - a file of the family that cannot be read: no family.json, the family
//!   listed unreadable.
//!
//! Each check runs again on a planted mistake and must fail there. The
//! wording check for family.json runs in tests/results_journal.rs.

use og_live::family::{
    Given, candidates_given, cli, families, family_page, family_state, member_summary, read_tsv,
    side_by_side, side_by_side_in, site_page, step41_order, with_candidates,
};
use pycompat::json::{Dump, Value, dumps, loads};
use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::{Duration, SystemTime};

fn repo() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../..")
}

fn demo() -> PathBuf {
    repo().join("tests/fixtures/live/demo/DEMO-FAMILY")
}

fn tmp(name: &str) -> PathBuf {
    static N: AtomicUsize = AtomicUsize::new(0);
    let d = Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!(
        "family-page-{name}-{}-{}",
        std::process::id(),
        N.fetch_add(1, Ordering::Relaxed)
    ));
    let _ = fs::remove_dir_all(&d);
    fs::create_dir_all(&d).unwrap();
    fs::canonicalize(&d).unwrap()
}

fn copy_tree(from: &Path, to: &Path) {
    fs::create_dir_all(to).unwrap();
    for e in fs::read_dir(from).unwrap() {
        let e = e.unwrap();
        let (src, dst) = (e.path(), to.join(e.file_name()));
        let ft = e.file_type().unwrap();
        if ft.is_symlink() {
            std::os::unix::fs::symlink(fs::read_link(&src).unwrap(), &dst).unwrap();
        } else if ft.is_dir() {
            copy_tree(&src, &dst);
        } else {
            fs::copy(&src, &dst).unwrap();
        }
    }
}

/// tests/fixtures/families_live as a GENOME_DIR, with a sample run in no family.
fn genome() -> PathBuf {
    let root = tmp("genome");
    copy_tree(&repo().join("tests/fixtures/families_live"), &root);
    fs::create_dir_all(root.join("loner/logs")).unwrap();
    root
}

fn j(text: &str) -> Value {
    loads(text).unwrap_or_else(|e| panic!("{e}: {text}"))
}

fn show(v: &Value) -> String {
    dumps(v, &Dump::default()).unwrap()
}

fn at<'a>(v: &'a Value, path: &[&str]) -> &'a Value {
    static NULL: Value = Value::Null;
    let mut v = v;
    for k in path {
        v = match v {
            Value::Object(o) => o.get(*k).unwrap_or(&NULL),
            Value::Array(a) => k
                .parse::<usize>()
                .ok()
                .and_then(|i| a.get(i))
                .unwrap_or(&NULL),
            _ => &NULL,
        };
    }
    v
}

fn list(v: &Value) -> &[Value] {
    match v {
        Value::Array(a) => a,
        _ => &[],
    }
}

fn st(v: &Value) -> &str {
    match v {
        Value::Str(s) => s,
        _ => "",
    }
}

fn state(root: &Path, name: &str) -> Value {
    family_state(name, &root.join(name)).unwrap()
}

// --- the family's state ---------------------------------------------------------------------

const CHILDREN: &str = r#"[{"child": "child_a", "father": "father_a", "mother": "mother_a",
  "segments": "/family/family_a/child_a_segments.svg"},
 {"child": "child_b", "father": "father_a", "mother": "mother_a", "segments": null},
 {"child": "father_a", "father": "gf_a", "mother": "gm_a", "segments": null},
 {"child": "mother_a", "father": null, "mother": "gm_m", "segments": null}]"#;
const FILES: [&str; 9] = [
    "child_a_segments.svg",
    "couple_father_a_mother_a.md",
    "couple_father_a_mother_a.tsv",
    "family_plan.tsv",
    "family_report.md",
    "ibd_pairs.tsv",
    "ibd_report.md",
    "members.tsv",
    "relationships.tsv",
];

fn children_problems(s: &Value) -> Vec<String> {
    let mut out = Vec::new();
    if *at(s, &["children"]) != j(CHILDREN) {
        out.push(format!("children {}", show(at(s, &["children"]))));
    }
    let seq: Vec<(&str, bool)> = list(at(s, &["members"]))
        .iter()
        .map(|m| {
            (
                st(at(m, &["sample"])),
                *at(m, &["sequenced"]) == Value::Bool(true),
            )
        })
        .collect();
    let want = [
        ("child_a", true),
        ("child_b", true),
        ("father_a", true),
        ("mother_a", true),
        ("gf_a", true),
        ("gm_a", true),
        ("gm_m", false),
    ];
    if seq != want {
        out.push(format!("members {seq:?}"));
    }
    out
}

#[test]
fn the_family_folders_are_runs_of_kind_family() {
    let root = genome();
    let fams = families(&root);
    let names: Vec<&str> = fams.iter().map(|f| st(at(f, &["name"]))).collect();
    assert_eq!(
        names,
        ["family_a", "family_b"],
        "child_a and loner are not families"
    );
    assert!(
        fams.iter()
            .all(|f| st(at(f, &["kind"])) == "family"
                && *at(f, &["unreadable"]) == Value::Bool(false))
    );
    assert_eq!(
        st(at(&fams[0], &["dir"])),
        root.join("family_a").to_str().unwrap()
    );
    assert_eq!(list(at(&fams[0], &["members"])).len(), 7);
    // the command line lists the same
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = cli(
        &["families".into(), root.to_string_lossy().into_owned()],
        &mut out,
        &mut err,
    );
    assert_eq!(code, 0);
    assert_eq!(j(&String::from_utf8(out).unwrap()), Value::Array(fams));
}

#[test]
fn family_a_state_is_what_the_steps_wrote() {
    let root = genome();
    let s = state(&root, "family_a");
    assert_eq!(children_problems(&s), Vec::<String>::new());
    let rel: Vec<(&str, &str, &str)> = list(at(&s, &["relationships"]))
        .iter()
        .map(|r| {
            (
                st(at(r, &["a"])),
                st(at(r, &["b"])),
                st(at(r, &["verdict"])),
            )
        })
        .collect();
    assert_eq!(
        rel,
        [("child_a", "father_a", "ok"), ("child_a", "child_b", "ok")],
        "the relationship check as step 38 wrote it"
    );
    assert_eq!(
        *at(&s, &["ibd_pairs"]),
        j(
            r#"[{"a": "child_a", "b": "child_b", "declared": "full siblings", "shared_cm": "2550.0", "degree": "first (full siblings)"}]"#
        ),
        "the shared segments as step 38 wrote them"
    );
    assert_eq!(
        *at(&s, &["couples"]),
        j(
            r#"[{"a": "father_a", "b": "mother_a", "report": "/family/family_a/couple_father_a_mother_a.md",
               "table": "/family/family_a/couple_father_a_mother_a.tsv", "screened": true, "shared": 1, "x_linked": 1},
              {"a": "gf_a", "b": "gm_a", "report": null, "table": null, "screened": false, "shared": 0, "x_linked": 0}]"#
        ),
        "the couple screen counted from step 39's table, the grandparents not screened"
    );
    let plan: Vec<&str> = list(at(&s, &["plan"]))
        .iter()
        .map(|r| st(at(r, &["analysis"])))
        .collect();
    assert_eq!(plan, ["joint_call", "couple_carrier"]);
    let files: Vec<&str> = list(at(&s, &["files"])).iter().map(st).collect();
    assert_eq!(
        files, FILES,
        "exactly the result files that exist (notes.txt is written by no step)"
    );
    assert_eq!(st(at(&s, &["pedigree"])), "/family/family_a/pedigree.svg");
    let b = state(&root, "family_b");
    assert_eq!(
        *at(&b, &["couples"]),
        j(
            r#"[{"a": "partner_x", "b": "partner_y", "report": null, "table": null, "screened": false, "shared": 0, "x_linked": 0}]"#
        ),
        "family_b (pedigree.tsv only): a couple not screened yet"
    );
}

#[test]
fn control_a_swapped_parent_in_the_state() {
    let root = genome();
    let p = root.join("family_a/members.tsv");
    let text = fs::read_to_string(&p).unwrap();
    fs::write(
        &p,
        text.replace(
            "child_a\tfemale\tfather_a\tmother_a",
            "child_a\tfemale\tgf_a\tmother_a",
        ),
    )
    .unwrap();
    assert!(
        !children_problems(&state(&root, "family_a")).is_empty(),
        "control: a swapped parent was NOT caught"
    );
}

#[test]
fn a_link_out_of_the_folder_is_not_listed_or_read() {
    let root = genome();
    let rep = root.join("family_a/ibd_report.md");
    fs::remove_file(&rep).unwrap();
    std::os::unix::fs::symlink(repo().join("README.md"), &rep).unwrap();
    let outside = root.join("outside-folder-name.tsv");
    fs::write(&outside, "sample_a\tsample_b\tverdict\nx\ty\tok\n").unwrap();
    let rel = root.join("family_a/relationships.tsv");
    fs::remove_file(&rel).unwrap();
    std::os::unix::fs::symlink(&outside, &rel).unwrap();
    let s = state(&root, "family_a");
    let files: Vec<&str> = list(at(&s, &["files"])).iter().map(st).collect();
    assert!(
        !files.contains(&"ibd_report.md") && !files.contains(&"relationships.tsv"),
        "{files:?}"
    );
    assert!(
        list(at(&s, &["relationships"])).is_empty(),
        "a table linking out of the folder is not read"
    );
    assert!(!show(&s).contains("outside-folder-name"));
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = cli(
        &[
            "files".into(),
            root.join("family_a").to_string_lossy().into_owned(),
        ],
        &mut out,
        &mut err,
    );
    assert_eq!(code, 0);
    assert_eq!(
        j(&String::from_utf8(out).unwrap()),
        *at(&s, &["files"]),
        "og live family files is the state's list"
    );
}

#[test]
fn a_members_file_with_no_rows_is_an_error() {
    let root = genome();
    let empty = root.join("family_empty");
    fs::create_dir_all(&empty).unwrap();
    let head = fs::read_to_string(root.join("family_a/members.tsv")).unwrap();
    fs::write(
        empty.join("members.tsv"),
        head.lines().next().unwrap().to_owned() + "\n",
    )
    .unwrap();
    let s = state(&root, "family_empty");
    assert!(
        st(at(&s, &["error"])).starts_with("the pedigree could not be read: "),
        "{}",
        show(&s)
    );
    let fam = families(&root)
        .into_iter()
        .find(|f| st(at(f, &["name"])) == "family_empty")
        .unwrap();
    assert!(
        list(at(&fam, &["members"])).is_empty() && *at(&fam, &["unreadable"]) == Value::Bool(true)
    );
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = cli(
        &["svg".into(), empty.to_string_lossy().into_owned()],
        &mut out,
        &mut err,
    );
    assert_eq!(code, 2, "no pedigree to draw");
    let page = family_page("family_empty", &empty).unwrap();
    assert_eq!(*at(&page, &["layout"]), Value::Null);
}

/// chmod 0 blocks the read unless the test runs as root.
fn blocked(p: &Path) -> bool {
    use std::os::unix::fs::PermissionsExt;
    fs::set_permissions(p, fs::Permissions::from_mode(0o000)).unwrap();
    fs::read(p).is_err()
}

fn unblock(p: &Path) {
    use std::os::unix::fs::PermissionsExt;
    fs::set_permissions(p, fs::Permissions::from_mode(0o644)).unwrap();
}

#[test]
fn a_file_that_cannot_be_read_gives_no_page_and_an_unreadable_family() {
    let root = genome();
    let d = root.join("family_a");
    let rel = d.join("relationships.tsv");
    if !blocked(&rel) {
        unblock(&rel);
        eprintln!("[SKIP] chmod 0 does not stop this user reading");
        return;
    }
    let page = family_page("family_a", &d);
    let fams = families(&root);
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = cli(
        &["page".into(), d.to_string_lossy().into_owned()],
        &mut out,
        &mut err,
    );
    unblock(&rel);
    assert!(page.is_err(), "family.json is not built");
    let err = String::from_utf8(err).unwrap();
    assert!(
        code == 1 && out.is_empty() && err.contains("PermissionError") && !err.contains("/"),
        "{code} {err}"
    );
    let a = &fams[0];
    assert!(
        *at(a, &["unreadable"]) == Value::Bool(true) && list(at(a, &["members"])).len() == 7,
        "{}",
        show(a)
    );
    let reopened = families(&root);
    assert_eq!(
        *at(&reopened[0], &["unreadable"]),
        Value::Bool(false),
        "readable again at once"
    );
    let mem = d.join("members.tsv");
    if blocked(&mem) {
        let fams = families(&root);
        unblock(&mem);
        assert!(
            list(at(&fams[0], &["members"])).is_empty()
                && *at(&fams[0], &["unreadable"]) == Value::Bool(true)
        );
    } else {
        unblock(&mem);
    }
}

// --- the depth of the page: stretches, trio tables ----------------------------------------

fn segment_problems(page: &Value, d: &Path) -> Vec<String> {
    let rows = read_tsv(&d.join("DEMO-003_segments.tsv")).unwrap();
    let got = list(at(page, &["segments", "DEMO-003"]));
    let mut out = Vec::new();
    let key = |r: &Value| ["side", "chrom", "start", "end"].map(|k| st(at(r, &[k])).to_owned());
    let want: Vec<[String; 4]> = rows
        .iter()
        .map(|r| ["side", "chrom", "start", "end"].map(|k| r.get(k).cloned().unwrap_or_default()))
        .collect();
    if rows.is_empty() || got.iter().map(key).collect::<Vec<_>>() != want {
        out.push(format!(
            "{} stretches against {} rows",
            got.len(),
            rows.len()
        ));
    }
    let chroms: Vec<&str> = list(at(page, &["chrom_lengths"]))
        .iter()
        .map(|c| st(at(c, &["0"])))
        .collect();
    let want: Vec<String> = (1..=22)
        .map(|i| format!("chr{i}"))
        .chain(["chrX".into()])
        .collect();
    if chroms != want {
        out.push(format!("lengths {chroms:?}"));
    }
    if *at(page, &["chrom_lengths", "0", "1"]) != Value::Int(248956422) {
        out.push("chr1's length".into());
    }
    out
}

fn pair_problems(page: &Value) -> Vec<String> {
    let rank = |p: &str| match p {
        "in trans" => 0,
        "in cis" => 2,
        _ => 1,
    };
    let phases: Vec<&str> = list(at(page, &["compound_het"]))
        .iter()
        .map(|r| st(at(r, &["phase"])))
        .collect();
    let mut sorted = phases.clone();
    sorted.sort_by_key(|p| rank(p));
    let mut out = Vec::new();
    if phases.first() != Some(&"in trans") || phases != sorted {
        out.push(format!("order {phases:?}"));
    }
    if list(at(page, &["de_novo"])).len() != 2 {
        out.push("de novo rows".into());
    }
    out
}

#[test]
fn stretches_and_trio_tables_as_the_steps_wrote_them() {
    let page = family_page("DEMO-FAMILY", &demo()).unwrap();
    assert_eq!(segment_problems(&page, &demo()), Vec::<String>::new());
    let mut wrong = page.clone();
    if let Value::Object(o) = &mut wrong
        && let Some(Value::Object(s)) = o.get_mut("segments")
        && let Some(Value::Array(a)) = s.get_mut("DEMO-003")
    {
        a.reverse();
    }
    assert!(
        !segment_problems(&wrong, &demo()).is_empty(),
        "control: the stretches in another order"
    );
    assert_eq!(
        pair_problems(&page),
        Vec::<String>::new(),
        "the compound-het pairs in trans first, both de novo rows"
    );
    let mut planted = page.clone();
    let file_order: Vec<Value> = read_tsv(&demo().join("compound_het.tsv"))
        .unwrap()
        .iter()
        .map(|r| {
            Value::Object(
                r.iter()
                    .map(|(k, v)| (k.clone(), Value::Str(v.clone())))
                    .collect(),
            )
        })
        .collect();
    if let Value::Object(o) = &mut planted {
        o.insert("compound_het".into(), Value::Array(file_order));
    }
    assert!(
        !pair_problems(&planted).is_empty(),
        "control: the pairs left in file order"
    );
    // the de novo row's own columns only: the parents' depth and quality are not on the page
    let row = at(&page, &["de_novo", "0"]);
    assert!(
        *at(row, &["father_dp"]) == Value::Null && st(at(row, &["child_ab"])) == "0.47",
        "{}",
        show(row)
    );
}

fn bare_ok(p: &Value) -> bool {
    *at(p, &["trio_files"]) == j(r#"{"de_novo.tsv": false, "compound_het.tsv": false}"#)
        && list(at(p, &["de_novo"])).is_empty()
        && list(at(p, &["compound_het"])).is_empty()
        && *at(p, &["segments"]) == j("{}")
}

#[test]
fn a_family_without_the_trio_files_says_so() {
    let bare = tmp("bare");
    fs::copy(demo().join("members.tsv"), bare.join("members.tsv")).unwrap();
    assert!(bare_ok(&family_page("bare", &bare).unwrap()));
    fs::copy(demo().join("de_novo.tsv"), bare.join("de_novo.tsv")).unwrap();
    assert!(
        !bare_ok(&family_page("bare", &bare).unwrap()),
        "control: a planted de_novo.tsv"
    );
}

// --- the members side by side -------------------------------------------------------------

fn summary(cyp: &str, cpic: &[(&str, &str, &str)], hap: &str, y: Option<&str>) -> Value {
    let genes: Vec<String> = cpic
        .iter()
        .map(|(g, d, p)| format!(r#"{{"gene": "{g}", "diplotype": "{d}", "phenotype": "{p}"}}"#))
        .collect();
    let y = y
        .map(|y| format!(r#", "y_tree": {{"state": "ok", "data": {{"haplogroup": "{y}"}}}}"#))
        .unwrap_or_default();
    j(&format!(
        r#"{{"cyp2d6": {{"calls": {{"pypgx": "{cyp}", "Cyrius": "{cyp}"}}, "agree": true}}, "sections": {{
           "cpic": {{"state": "ok", "data": {{"genes": [{}]}}}},
           "haplogroup": {{"state": "ok", "data": {{"haplogroup": "{hap}"}}}},
           "bloodtype": {{"state": "ok", "data": {{"abo": "A1", "rhd": "D+"}}}}{y}}}}}"#,
        genes.join(", ")
    ))
}

fn set(v: &mut Value, path: &[&str], x: Value) {
    let mut v = v;
    for k in &path[..path.len() - 1] {
        v = match v {
            Value::Object(o) => o
                .entry((*k).to_owned())
                .or_insert_with(|| Value::Object(Default::default())),
            _ => panic!("not an object at {k}"),
        };
    }
    if let Value::Object(o) = v {
        o.insert(path[path.len() - 1].to_owned(), x);
    }
}

fn stale_summary(state: &str) -> Value {
    let mut s = summary(
        "*1/*1",
        &[("DPYD", "Reference/Reference", "Normal Metabolizer")],
        "H1",
        None,
    );
    set(
        &mut s,
        &["sections", "bloodtype"],
        j(&format!(
            r#"{{"state": "{state}", "data": {{"abo": "A1", "rhd": "D+", "systems": [{{"system": "ABO", "phenotype": "A1"}},
               {{"system": "RhD", "phenotype": "D+"}}, {{"system": "Kell", "phenotype": "K-k+"}},
               {{"system": "Duffy", "phenotype": "Fy(a+b-)"}}]}}}}"#
        )),
    );
    s
}

/// Every side-by-side source present and complete, with PharmCAT's CYP2D6
/// call, the Y haplogroup from y_haplogroup alone (the fallback), and one
/// section, `key`, in the given state.
fn source_summary(key: &str, state: &str) -> Value {
    let mut s = summary(
        "*1/*1",
        &[("DPYD", "Reference/Reference", "Normal Metabolizer")],
        "H1",
        None,
    );
    set(
        &mut s,
        &["cyp2d6", "calls", "PharmCAT"],
        Value::Str("*1/*1".into()),
    );
    set(
        &mut s,
        &["sections", "hla"],
        j(
            r#"{"state": "ok", "data": {"loci": [{"gene": "HLA-A", "alleles": ["A*01:01", "A*02:01"]}]}}"#,
        ),
    );
    set(
        &mut s,
        &["sections", "y_haplogroup"],
        j(r#"{"state": "ok", "data": {"haplogroup": "R-M269"}}"#),
    );
    set(
        &mut s,
        &["sections", key, "state"],
        Value::Str(state.into()),
    );
    s
}

fn names(v: &[&str]) -> Vec<String> {
    v.iter().map(|s| (*s).to_owned()).collect()
}

fn sums(pairs: &[(&str, Option<Value>)]) -> HashMap<String, Option<Value>> {
    pairs
        .iter()
        .map(|(k, v)| ((*k).to_owned(), v.clone()))
        .collect()
}

fn abc() -> HashMap<String, Option<Value>> {
    sums(&[
        (
            "A",
            Some(summary(
                "*1/*4",
                &[
                    ("CYP2D6", "*1/*4", "Intermediate Metabolizer"),
                    ("CYP2C19", "*1/*2", "Intermediate Metabolizer"),
                    ("DPYD", "Reference/Reference", "Normal Metabolizer"),
                ],
                "H1",
                Some("I-M253"),
            )),
        ),
        (
            "B",
            Some(summary(
                "*1/*1",
                &[("DPYD", "Reference/Reference", "Normal Metabolizer")],
                "U5",
                None,
            )),
        ),
        ("C", None),
    ])
}

fn row<'a>(sb: &'a Value, group: &str, label: &str) -> Option<&'a Value> {
    list(at(sb, &["rows"]))
        .iter()
        .find(|r| st(at(r, &["group"])) == group && st(at(r, &["label"])) == label)
}

fn values(sb: &Value, group: &str, label: &str) -> Vec<String> {
    row(sb, group, label)
        .map(|r| {
            list(at(r, &["values"]))
                .iter()
                .map(|v| st(v).to_owned())
                .collect()
        })
        .unwrap_or_default()
}

fn side_problems(sb: &Value) -> Vec<String> {
    let want: [(&str, &str, [&str; 3]); 6] = [
        (
            "Medicines",
            "CYP2D6 (callers)",
            [
                "pypgx *1/*4, Cyrius *1/*4 (callers agree)",
                "pypgx *1/*1, Cyrius *1/*1 (callers agree)",
                "no summary yet",
            ],
        ),
        (
            "Medicines",
            "CYP2C19",
            ["*1/*2 \"Intermediate Metabolizer\"", ".", "no summary yet"],
        ),
        (
            "Medicines",
            "CYP2D6",
            ["*1/*4 \"Intermediate Metabolizer\"", ".", "no summary yet"],
        ),
        (
            "Medicines",
            "DPYD",
            [
                "Reference/Reference \"Normal Metabolizer\"",
                "Reference/Reference \"Normal Metabolizer\"",
                "no summary yet",
            ],
        ),
        (
            "Haplogroups",
            "Y chromosome",
            ["I-M253", ".", "no summary yet"],
        ),
        ("Blood groups", "ABO", ["A1", "A1", "no summary yet"]),
    ];
    let mut out: Vec<String> = want
        .iter()
        .filter(|(g, l, v)| values(sb, g, l) != *v)
        .map(|(g, l, v)| format!("{g} {l}: {:?} != {v:?}", values(sb, g, l)))
        .collect();
    if *at(sb, &["members"]) != j(r#"["A", "B", "C"]"#) {
        out.push("members".into());
    }
    let order = ["Medicines", "HLA", "Haplogroups", "Blood groups"];
    let groups: Vec<usize> = list(at(sb, &["rows"]))
        .iter()
        .map(|r| {
            order
                .iter()
                .position(|g| *g == st(at(r, &["group"])))
                .unwrap()
        })
        .collect();
    if !groups.is_sorted() {
        out.push(format!("groups out of order {groups:?}"));
    }
    out
}

#[test]
fn the_members_side_by_side() {
    let sb = side_by_side(&names(&["A", "B", "C"]), &abc());
    assert_eq!(side_problems(&sb), Vec::<String>::new());
    let mut dropped = abc();
    dropped.insert("A".into(), None);
    assert!(
        !side_problems(&side_by_side(&names(&["A", "B", "C"]), &dropped)).is_empty(),
        "control: member A's summary dropped"
    );
}

#[test]
fn a_summary_linking_out_of_its_run_is_never_read() {
    let root = tmp("root");
    let run = root.join("A");
    fs::create_dir_all(&run).unwrap();
    let outside = root
        .parent()
        .unwrap()
        .join(format!("outside-{}.json", std::process::id()));
    fs::write(&outside, show(abc()["A"].as_ref().unwrap())).unwrap();
    std::os::unix::fs::symlink(&outside, run.join("summary.json")).unwrap();
    assert_eq!(
        member_summary(&root, Some(&run)),
        None,
        "a summary.json linking out of its run is not read"
    );
    // control: the same file inside the run is read
    fs::remove_file(run.join("summary.json")).unwrap();
    fs::copy(&outside, run.join("summary.json")).unwrap();
    fs::remove_file(&outside).unwrap();
    assert!(
        member_summary(&root, Some(&run)).is_some(),
        "control: the inside check reads a plain file"
    );
    // a run folder that links out of the root: never read either
    let elsewhere = tmp("elsewhere");
    fs::copy(run.join("summary.json"), elsewhere.join("summary.json")).unwrap();
    std::os::unix::fs::symlink(&elsewhere, root.join("B")).unwrap();
    assert_eq!(member_summary(&root, Some(&root.join("B"))), None);
    fs::write(run.join("summary.json"), "[1, 2]").unwrap();
    assert_eq!(
        member_summary(&root, Some(&run)),
        None,
        "a summary that is not an object"
    );
    fs::write(run.join("summary.json"), "{not json").unwrap();
    assert_eq!(member_summary(&root, Some(&run)), None);
}

fn stale_problems(sb: &Value) -> Vec<String> {
    let mut out: Vec<String> = list(at(sb, &["rows"]))
        .iter()
        .filter_map(|r| {
            let want = j(&format!(
                "[{}, false, false]",
                st(at(r, &["group"])) == "Blood groups"
            ));
            (*at(r, &["stale"]) != want)
                .then(|| format!("{}: {}", st(at(r, &["label"])), show(at(r, &["stale"]))))
        })
        .collect();
    if list(at(sb, &["rows"])).is_empty() {
        out.push("no rows".into());
    }
    out
}

fn blood_problems(sb: &Value) -> Vec<String> {
    let rows: Vec<&Value> = list(at(sb, &["rows"]))
        .iter()
        .filter(|r| st(at(r, &["group"])) == "Blood groups")
        .collect();
    let labels: Vec<&str> = rows.iter().map(|r| st(at(r, &["label"]))).collect();
    let mut out = Vec::new();
    if labels != ["ABO", "RhD", "Kell", "Duffy"] {
        out.push(format!("blood rows {labels:?}"));
    }
    if values(sb, "Blood groups", "Kell")
        .first()
        .map(String::as_str)
        != Some("K-k+")
        || values(sb, "Blood groups", "Duffy")
            .first()
            .map(String::as_str)
            != Some("Fy(a+b-)")
    {
        out.push("values".into());
    }
    out
}

fn abc7(a: Value) -> HashMap<String, Option<Value>> {
    sums(&[
        ("A", Some(a)),
        ("B", Some(summary("*1/*1", &[], "U5", None))),
        ("C", None),
    ])
}

#[test]
fn a_stale_value_is_flagged_in_its_own_cell() {
    let sb = side_by_side(&names(&["A", "B", "C"]), &abc7(stale_summary("stale")));
    assert_eq!(stale_problems(&sb), Vec::<String>::new());
    let complete = side_by_side(&names(&["A", "B", "C"]), &abc7(stale_summary("ok")));
    assert!(
        !stale_problems(&complete).is_empty(),
        "control: the stale section made complete"
    );
    let label = |r: &Value| st(at(r, &["label"])).to_owned();
    for (key, want) in [
        ("cpic", vec!["CYP2D6 (callers)", "DPYD"]),
        ("hla", vec!["HLA-A"]),
        ("haplogroup", vec!["Mitochondrial"]),
        ("y_haplogroup", vec!["Y chromosome"]),
    ] {
        let sb = side_by_side(
            &names(&["A"]),
            &sums(&[("A", Some(source_summary(key, "stale")))]),
        );
        let mut got: Vec<String> = list(at(&sb, &["rows"]))
            .iter()
            .filter(|r| *at(r, &["stale", "0"]) == Value::Bool(true))
            .map(label)
            .collect();
        got.sort();
        assert_eq!(
            got, want,
            "{key}: each stale source flags its own cells only"
        );
    }
    let sb = side_by_side(
        &names(&["A"]),
        &sums(&[("A", Some(source_summary("cpic", "ok")))]),
    );
    assert!(
        list(at(&sb, &["rows"]))
            .iter()
            .all(|r| *at(r, &["stale", "0"]) == Value::Bool(false)),
        "control: the CPIC section made complete flags nothing"
    );
}

#[test]
fn every_blood_group_system_in_step_41s_order() {
    let sb = side_by_side(&names(&["A", "B", "C"]), &abc7(stale_summary("stale")));
    assert_eq!(blood_problems(&sb), Vec::<String>::new());
    let mut emptied = stale_summary("stale");
    set(
        &mut emptied,
        &["sections", "bloodtype", "data", "systems"],
        j("[]"),
    );
    assert!(
        !blood_problems(&side_by_side(&names(&["A", "B", "C"]), &abc7(emptied))).is_empty(),
        "control: the systems list emptied"
    );
    let with_systems = |drop: &str, kidd: bool| {
        let mut s = stale_summary("ok");
        let mut systems: Vec<Value> =
            list(at(&s, &["sections", "bloodtype", "data", "systems"])).to_vec();
        if kidd {
            systems.push(j(r#"{"system": "Kidd", "phenotype": "Jk(a+b+)"}"#));
        }
        systems.retain(|x| st(at(x, &["system"])) != drop);
        set(
            &mut s,
            &["sections", "bloodtype", "data", "systems"],
            Value::Array(systems),
        );
        s
    };
    let merged = |a: Value, b: Value, order: &[String]| -> Vec<String> {
        let sb = side_by_side_in(
            &names(&["A", "B"]),
            &sums(&[("A", Some(a)), ("B", Some(b))]),
            order,
        );
        list(at(&sb, &["rows"]))
            .iter()
            .filter(|r| st(at(r, &["group"])) == "Blood groups")
            .map(|r| st(at(r, &["label"])).to_owned())
            .collect()
    };
    let order = step41_order();
    assert_eq!(
        merged(with_systems("Kell", false), stale_summary("ok"), &order),
        ["ABO", "RhD", "Kell", "Duffy"]
    );
    let got = merged(
        with_systems("Duffy", true),
        with_systems("Kell", true),
        &order,
    );
    assert_eq!(
        got,
        ["ABO", "RhD", "Kell", "Duffy", "Kidd"],
        "a system only a later member has keeps step 41's order"
    );
    assert_ne!(
        merged(with_systems("Duffy", true), with_systems("Kell", true), &[]),
        ["ABO", "RhD", "Kell", "Duffy", "Kidd"],
        "control: step 41's order taken away"
    );
    assert_eq!(
        &order[..5],
        ["ABO", "RhD", "RhCE", "Kell", "Duffy"],
        "the page's order is og bloodtype systems'"
    );
}

type Summaries = HashMap<String, Option<Value>>;

fn malformed_problems(read: &dyn Fn(&[String], &Summaries) -> Value) -> Vec<String> {
    let mut out = Vec::new();
    let three = names(&["A", "B", "C"]);
    for bad in [
        r#"{"sections": ["not", "a", "mapping"]}"#,
        r#"{"cyp2d6": "CYP2D6 *1/*4"}"#,
        r#"{"sections": {"cpic": {"state": "ok", "data": {"genes": [{"gene": ["CYP2D6"]}]}}}}"#,
    ] {
        let mut s = abc();
        s.insert("B".into(), Some(j(bad)));
        let got = read(&three, &s);
        let mut none = abc();
        none.insert("B".into(), None);
        let want = side_by_side(&three, &none);
        for r in list(at(&want, &["rows"])) {
            let (g, l) = (st(at(r, &["group"])), st(at(r, &["label"])));
            let cols = values(&got, g, l);
            if cols.first().map(String::as_str) != Some(st(at(r, &["values", "0"])))
                || cols.get(1).map(String::as_str) != Some("summary not readable")
            {
                out.push(format!("{bad}: {l} {cols:?}"));
            }
        }
    }
    out
}

#[test]
fn a_summary_of_a_shape_never_written_costs_its_column_only() {
    assert_eq!(
        malformed_problems(&|n, s| side_by_side(n, s)),
        Vec::<String>::new()
    );
    // the planted mistake: a bad summary read as one with no rows, never said
    let quiet = |n: &[String], s: &HashMap<String, Option<Value>>| {
        let mut s = s.clone();
        s.insert("B".into(), Some(j("{}")));
        side_by_side(n, &s)
    };
    assert!(
        !malformed_problems(&quiet).is_empty(),
        "control: a bad summary shown as '.' in every row"
    );
}

// --- each child's trio status ---------------------------------------------------------------

fn trio_problems(status: &Value, want: &[(&str, &str, Option<bool>, &str)]) -> Vec<String> {
    want.iter()
        .filter_map(|(child, key, assessed, words)| {
            let got = at(status, &[child, key]);
            let a = match at(got, &["assessed"]) {
                Value::Bool(b) => Some(*b),
                _ => None,
            };
            (a != *assessed || !st(at(got, &["why"])).contains(words))
                .then(|| format!("{child} {key}: {}", show(got)))
        })
        .collect()
}

struct Trio {
    dir: PathBuf,
    plan: PathBuf,
    plan_text: String,
    secret: PathBuf,
}

/// The demo's trio with a second child, DEMO-004, whose mother DEMO-005 is
/// not sequenced, and compound_het.tsv a link out of the folder.
fn trio() -> Trio {
    let work = tmp("trio");
    let dir = work.join("trio");
    copy_tree(&demo(), &dir);
    let mut m = fs::read_to_string(dir.join("members.tsv")).unwrap();
    m.push_str(
        "DEMO-004\tfemale\tDEMO-001\tDEMO-005\t\tyes\tyes\nDEMO-005\tfemale\t\t\t\tyes\tno\n",
    );
    fs::write(dir.join("members.tsv"), m).unwrap();
    let plan = dir.join("family_plan.tsv");
    let plan_text = fs::read_to_string(&plan).unwrap();
    fs::write(
        &plan,
        plan_text.replace(
            "DEMO-001,DEMO-002,DEMO-003\t",
            "DEMO-001,DEMO-002,DEMO-003,DEMO-004\t",
        ),
    )
    .unwrap();
    let secret = work.join("elsewhere-private-name.tsv");
    fs::copy(dir.join("compound_het.tsv"), &secret).unwrap();
    fs::remove_file(dir.join("compound_het.tsv")).unwrap();
    std::os::unix::fs::symlink(&secret, dir.join("compound_het.tsv")).unwrap();
    Trio {
        dir,
        plan,
        plan_text,
        secret,
    }
}

fn status(t: &Trio) -> (Value, Value) {
    let p = family_page("trio", &t.dir).unwrap();
    (at(&p, &["trio_status"]).clone(), p)
}

const ASSESSED: [(&str, &str, Option<bool>, &str); 4] = [
    ("DEMO-003", "de_novo", Some(true), "both parents"),
    (
        "DEMO-004",
        "de_novo",
        Some(false),
        "DEMO-005 (the mother) is not",
    ),
    (
        "DEMO-003",
        "compound_het",
        Some(false),
        "link out of the family folder",
    ),
    (
        "DEMO-004",
        "compound_het",
        Some(false),
        "link out of the family folder",
    ),
];

#[test]
fn each_childs_trio_status() {
    let t = trio();
    let (s, page) = status(&t);
    assert_eq!(trio_problems(&s, &ASSESSED), Vec::<String>::new());
    assert_eq!(
        *at(&page, &["trio_file_state"]),
        j(r#"{"de_novo.tsv": "here", "compound_het.tsv": "rejected"}"#)
    );
    assert!(
        !show(&page).contains("elsewhere-private-name"),
        "the link's target is not in family.json"
    );
    fs::write(&t.plan, &t.plan_text).unwrap();
    assert_eq!(
        trio_problems(
            &status(&t).0,
            &[(
                "DEMO-004",
                "de_novo",
                Some(false),
                "DEMO-004 is not in the joint call"
            )]
        ),
        Vec::<String>::new()
    );
    fs::remove_file(t.dir.join("de_novo.tsv")).unwrap();
    assert_eq!(
        trio_problems(
            &status(&t).0,
            &[(
                "DEMO-003",
                "de_novo",
                Some(false),
                "de_novo.tsv is not here"
            )]
        ),
        Vec::<String>::new(),
        "a missing de_novo.tsv, said apart from a refused one"
    );
    fs::remove_file(t.dir.join("compound_het.tsv")).unwrap();
    fs::copy(&t.secret, t.dir.join("compound_het.tsv")).unwrap();
    assert!(
        !trio_problems(&status(&t).0, &ASSESSED[2..]).is_empty(),
        "control: the link replaced by a plain file"
    );
    fs::copy(demo().join("de_novo.tsv"), t.dir.join("de_novo.tsv")).unwrap();
    fs::write(
        &t.plan,
        t.plan_text.replace(
            "DEMO-001,DEMO-002,DEMO-003\t",
            "DEMO-002,DEMO-003,DEMO-004\t",
        ),
    )
    .unwrap();
    let want = [
        (
            "DEMO-003",
            "de_novo",
            Some(false),
            "DEMO-001 (the father) is not",
        ),
        (
            "DEMO-003",
            "compound_het",
            Some(true),
            "a parent is sequenced",
        ),
        (
            "DEMO-004",
            "de_novo",
            Some(false),
            "DEMO-001 (the father) is not",
        ),
        (
            "DEMO-004",
            "compound_het",
            Some(false),
            "need a parent sequenced and in the joint call",
        ),
    ];
    assert_eq!(
        trio_problems(&status(&t).0, &want),
        Vec::<String>::new(),
        "a sequenced parent left out of the joint call"
    );
    fs::write(&t.plan, &t.plan_text).unwrap();
    assert!(
        !trio_problems(&status(&t).0, &want[..1]).is_empty(),
        "control: the parent back in the joint call"
    );
}

fn left_out_of_the_joint_call() -> Value {
    let t = trio();
    fs::remove_file(t.dir.join("compound_het.tsv")).unwrap();
    fs::copy(&t.secret, t.dir.join("compound_het.tsv")).unwrap();
    fs::write(
        &t.plan,
        t.plan_text.replace(
            "DEMO-001,DEMO-002,DEMO-003\t",
            "DEMO-002,DEMO-003,DEMO-004\t",
        ),
    )
    .unwrap();
    status(&t).0
}

fn given(pairs: &[(&str, Option<Given>)]) -> HashMap<String, Option<Given>> {
    pairs.iter().map(|(k, v)| ((*k).to_owned(), *v)).collect()
}

#[test]
fn compound_het_pairs_by_the_childs_step_23_candidates() {
    let st9 = left_out_of_the_joint_call();
    let told = |g: &[(&str, Option<Given>)]| {
        let mut s = st9.clone();
        with_candidates(&mut s, &given(g));
        s
    };
    let mut bad = trio_problems(
        &told(&[
            ("DEMO-003", Some(Given::Absent)),
            ("DEMO-004", Some(Given::Here)),
        ]),
        &[
            (
                "DEMO-003",
                "compound_het",
                Some(false),
                "wrote no candidates for DEMO-003",
            ),
            ("DEMO-004", "compound_het", Some(false), "need a parent"),
        ],
    );
    bad.extend(trio_problems(
        &told(&[("DEMO-003", Some(Given::Here))]),
        &[(
            "DEMO-003",
            "compound_het",
            Some(true),
            "candidates (step 23) are here",
        )],
    ));
    for (g, words) in [
        (Some(Given::Newer), "newer than compound_het.tsv"),
        (None, "run is not on this site"),
    ] {
        bad.extend(trio_problems(
            &told(&[("DEMO-003", g)]),
            &[("DEMO-003", "compound_het", None, words)],
        ));
    }
    assert_eq!(bad, Vec::<String>::new());
    assert!(
        !trio_problems(
            &told(&[("DEMO-003", None)]),
            &[(
                "DEMO-003",
                "compound_het",
                Some(false),
                "wrote no candidates"
            )]
        )
        .is_empty(),
        "control: the candidates not known read as absent"
    );

    let root = tmp("candidates");
    let run3 = root.join("DEMO-003");
    fs::create_dir_all(run3.join("clinical")).unwrap();
    let vcf = run3.join("clinical/DEMO-003_clinical.vcf.gz");
    let none_yet = candidates_given(&root, Some(&run3), "DEMO-003", None);
    fs::write(&vcf, "").unwrap();
    let t = |s: u64| SystemTime::UNIX_EPOCH + Duration::from_secs(s);
    fs::File::options()
        .write(true)
        .open(&vcf)
        .unwrap()
        .set_modified(t(1000))
        .unwrap();
    let here = candidates_given(&root, Some(&run3), "DEMO-003", Some(t(2000)));
    let newer = candidates_given(&root, Some(&run3), "DEMO-003", Some(t(500)));
    let outside = root
        .parent()
        .unwrap()
        .join(format!("outside-{}.vcf.gz", std::process::id()));
    fs::write(&outside, "").unwrap();
    fs::remove_file(&vcf).unwrap();
    std::os::unix::fs::symlink(&outside, &vcf).unwrap();
    let linked = candidates_given(&root, Some(&run3), "DEMO-003", None);
    let got = [
        none_yet,
        here,
        newer,
        linked,
        candidates_given(&root, None, "DEMO-003", None),
    ];
    fs::remove_file(&outside).unwrap();
    assert_eq!(
        got,
        [
            Some(Given::Absent),
            Some(Given::Here),
            Some(Given::Newer),
            Some(Given::Absent),
            None
        ]
    );
}

#[test]
fn the_site_page_adds_runs_candidates_and_summaries() {
    let t = trio();
    let root = t.dir.parent().unwrap().to_path_buf();
    for n in ["DEMO-001", "DEMO-002", "DEMO-003"] {
        fs::create_dir_all(root.join(n).join("clinical")).unwrap();
    }
    fs::write(
        root.join("DEMO-001/summary.json"),
        show(abc()["A"].as_ref().unwrap()),
    )
    .unwrap();
    fs::write(
        root.join("DEMO-002/summary.json"),
        r#"{"sections": ["not", "a", "mapping"]}"#,
    )
    .unwrap();
    fs::remove_file(t.dir.join("compound_het.tsv")).unwrap();
    fs::copy(&t.secret, t.dir.join("compound_het.tsv")).unwrap();
    let dir = |n: &str| Value::Str(root.join(n).to_string_lossy().into_owned());
    let mut input = j(r#"{"runs": {"DEMO-001": {"url": "/runs/DEMO-001/", "child": false}}}"#);
    set(&mut input, &["candidates", "DEMO-003"], dir("DEMO-003"));
    for n in ["DEMO-001", "DEMO-002", "DEMO-003"] {
        set(&mut input, &["summaries", n], dir(n));
    }
    let page = site_page("trio", &t.dir, &root, &input).unwrap();
    assert_eq!(
        *at(&page, &["runs"]),
        *at(&input, &["runs"]),
        "the member runs as the site gave them"
    );
    let keys: Vec<&str> = match &page {
        Value::Object(o) => o.keys().map(String::as_str).collect(),
        _ => Vec::new(),
    };
    assert_eq!(keys.last(), Some(&"runs"));
    assert_eq!(
        trio_problems(
            at(&page, &["trio_status"]),
            &[
                (
                    "DEMO-003",
                    "compound_het",
                    Some(false),
                    "wrote no candidates for DEMO-003"
                ),
                (
                    "DEMO-004",
                    "compound_het",
                    None,
                    "DEMO-004's run is not on this site"
                ),
            ]
        ),
        Vec::<String>::new(),
        "candidates told only for the children the site named"
    );
    let sb = at(&page, &["side_by_side"]);
    assert_eq!(
        *at(sb, &["members"]),
        j(r#"["DEMO-001", "DEMO-002", "DEMO-003", "DEMO-004"]"#)
    );
    assert_eq!(
        values(sb, "Blood groups", "ABO"),
        [
            "A1",
            "summary not readable",
            "no summary yet",
            "no summary yet"
        ]
    );
    // the command line reads the same input on stdin's stand-in, a file
    let f = root.join("input.json");
    fs::write(&f, show(&input)).unwrap();
    let args: Vec<String> = [
        "page",
        t.dir.to_str().unwrap(),
        "--name",
        "trio",
        "--root",
        root.to_str().unwrap(),
        "--input",
        f.to_str().unwrap(),
    ]
    .iter()
    .map(|a| (*a).to_owned())
    .collect();
    let (mut out, mut err) = (Vec::new(), Vec::new());
    assert_eq!(
        cli(&args, &mut out, &mut err),
        0,
        "{}",
        String::from_utf8_lossy(&err)
    );
    assert_eq!(j(&String::from_utf8(out).unwrap()), page);
    // without the input: the page as family_page builds it, no runs, nobody side by side
    let plain = family_page("trio", &t.dir).unwrap();
    assert!(*at(&plain, &["runs"]) == Value::Null && *at(&plain, &["side_by_side"]) == Value::Null);
}

#[test]
fn ok_reads_in_words_on_the_page_and_as_written_in_the_state() {
    let page = family_page("DEMO-FAMILY", &demo()).unwrap();
    for rows in [
        at(&page, &["relationships"]),
        at(&page, &["state", "relationships"]),
    ] {
        assert!(
            list(rows)
                .iter()
                .all(|r| st(at(r, &["verdict"])) == "matches the declared relation"),
            "{}",
            show(rows)
        );
    }
    let s = family_state("DEMO-FAMILY", &demo()).unwrap();
    assert!(
        list(at(&s, &["relationships"]))
            .iter()
            .all(|r| st(at(r, &["verdict"])) == "ok")
    );
}
