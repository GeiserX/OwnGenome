use super::super::group_c::{self, GROUP_C};
use super::super::model::{GROUP_A, Spec};
use super::super::{Group, Input, group_b};
use super::*;
use crate::testutil::{Dir, run};
use pycompat::json::{Value, loads};

fn input(d: &Dir) -> Input {
    Input {
        sample: "S".into(),
        sample_dir: d.path("S"),
        declared_sex: None,
        assets: d.path("assets"),
        group: Group::All,
        tables: d.path("tables"),
    }
}

fn summary(d: &Dir) -> Value {
    std::fs::create_dir_all(d.path("S")).unwrap();
    collect(&input(d)).unwrap()
}

fn at<'a>(v: &'a Value, path: &[&str]) -> &'a Value {
    path.iter().fold(v, |v, k| match v {
        Value::Object(m) => m.get(*k).unwrap_or_else(|| panic!("no {k} in {path:?}")),
        Value::Array(a) => &a[k.parse::<usize>().unwrap()],
        _ => panic!("{path:?}"),
    })
}

fn has(v: &Value, key: &str) -> bool {
    matches!(v, Value::Object(m) if m.contains_key(key))
}

fn texts(v: &Value) -> Vec<String> {
    match v {
        Value::Array(a) => a
            .iter()
            .map(|x| match x {
                Value::Str(s) => s.clone(),
                _ => panic!("{x:?}"),
            })
            .collect(),
        _ => panic!("{v:?}"),
    }
}

fn word(v: &Value) -> &str {
    match v {
        Value::Str(s) => s,
        _ => panic!("{v:?}"),
    }
}

#[test]
fn the_order_holds_every_section_once() {
    let mut keys: Vec<&str> = GROUP_A
        .iter()
        .chain(group_b::GROUP_B.iter())
        .chain(GROUP_C.iter())
        .map(|s: &Spec| s.key)
        .collect();
    assert_eq!(keys.len(), ORDER.len());
    keys.sort_unstable();
    let mut order = ORDER.to_vec();
    order.sort_unstable();
    assert_eq!(keys, order);
}

#[test]
fn an_empty_folder_lists_every_step_and_every_line() {
    let d = Dir::new("post-empty");
    let v = summary(&d);
    let Value::Object(m) = &v else { panic!() };
    let keys: Vec<&str> = m.keys().map(String::as_str).collect();
    assert_eq!(
        keys,
        [
            "schema_version",
            "sample",
            "generated_utc",
            "run",
            "databases",
            "manifest",
            "sections",
            "cyp2d6",
            "not_run",
            "not_assessed"
        ]
    );
    let Value::Object(sections) = at(&v, &["sections"]) else {
        panic!()
    };
    assert_eq!(
        sections.keys().map(String::as_str).collect::<Vec<_>>(),
        ORDER
    );
    let not_run = texts(at(&v, &["not_run"]));
    assert_eq!(not_run.len(), 45);
    assert_eq!(not_run[0], "Coverage (mosdepth)");
    assert!(!not_run.contains(&"Reads kept by alignment".to_owned()));
    assert!(!not_run.contains(&"Carrier list".to_owned()));
    assert!(!not_run.contains(&"Couple carrier screen".to_owned()));
    let lines = texts(at(&v, &["not_assessed"]));
    let mut expected: Vec<&str> = NOT_ASSESSED.to_vec();
    expected.extend([
        PRS_NOT_ADJUSTED,
        HBA_NOT_ASSESSED,
        MEI_NOT_ASSESSED,
        VNTYPER_NOT_ASSESSED,
        MCA_NOT_ASSESSED,
    ]);
    assert_eq!(lines, expected);
}

#[test]
fn a_failed_step_says_so_under_not_run() {
    let d = Dir::new("post-failed");
    d.put(
        "S/logs/run_status.tsv",
        "meta\tstarted_epoch\t1577836800\nstep\t06\tfailed (exit 1)\nstep\t02\tfailed (exit 1)\n",
    );
    let v = summary(&d);
    let not_run = texts(at(&v, &["not_run"]));
    assert!(not_run.contains(&"ClinVar screen (step 06 failed)".to_owned()));
    // The read count is no step of its own, failed or not; the carrier
    // list shares step 06 with ClinVar.
    assert!(!not_run.iter().any(|x| x.starts_with("Reads kept")));
    assert!(!not_run.iter().any(|x| x.starts_with("Carrier list")));
    assert_eq!(word(at(&v, &["sections", "carrier", "state"])), "failed");
}

#[test]
fn step_35s_hba_count_assesses_hba() {
    let d = Dir::new("post-assessed");
    d.put(
        "S/paralogs/S_hba_copy_number.tsv",
        "metric\tvalue\nstatus\tcalled\ncopies\t4\nreading\tfour\n",
    );
    let v = summary(&d);
    assert_eq!(word(at(&v, &["sections", "parascopy", "state"])), "ok");
    let lines = texts(at(&v, &["not_assessed"]));
    assert!(!lines.contains(&HBA_NOT_ASSESSED.to_owned()));
    assert!(lines.contains(&VNTYPER_NOT_ASSESSED.to_owned()));
}

#[test]
fn a_stale_vntyper_reading_is_not_assessed() {
    // Only an ok section counts for VNtyper; any other state leaves the line.
    let mk = |state: State, status: &str| {
        let mut s = Section::new(
            GROUP_C
                .iter()
                .copied()
                .find(|s| s.key == "vntyper")
                .unwrap(),
            Ok(None),
            None,
        );
        s.state = state;
        s.data = Some(Data::Json(
            loads(&format!(r#"{{"status": "{status}"}}"#)).unwrap(),
        ));
        s
    };
    let lines = not_assessed(&[mk(State::Ok, "positive")]);
    assert!(!lines.contains(&VNTYPER_NOT_ASSESSED));
    let lines = not_assessed(&[mk(State::Stale, "positive")]);
    assert!(lines.contains(&VNTYPER_NOT_ASSESSED));
    let lines = not_assessed(&[mk(State::Ok, "failed")]);
    assert!(lines.contains(&VNTYPER_NOT_ASSESSED));
}

#[test]
fn paraphase_needs_both_hba_genes_called() {
    let paralogs = |genes: &str| {
        let mut s = Section::new(
            GROUP_C
                .iter()
                .copied()
                .find(|s| s.key == "paralogs")
                .unwrap(),
            Ok(None),
            None,
        );
        s.state = State::Ok;
        s.data = Some(Data::Json(
            loads(&format!(r#"{{"genes": {genes}}}"#)).unwrap(),
        ));
        s
    };
    let both = r#"[{"region": "hba", "gene": "HBA1", "status": "called"},
                   {"region": "hba", "gene": "HBA2", "status": "called"}]"#;
    let one = r#"[{"region": "hba", "gene": "HBA1", "status": "called"},
                  {"region": "hba", "gene": "HBA2", "status": "no call"}]"#;
    // A row without a region, the shape that stopped the Python collector.
    let odd = r#"[{"gene": "HBA1", "status": "called"}, 7]"#;
    assert!(!not_assessed(&[paralogs(both)]).contains(&HBA_NOT_ASSESSED));
    assert!(not_assessed(&[paralogs(one)]).contains(&HBA_NOT_ASSESSED));
    assert!(not_assessed(&[paralogs(odd)]).contains(&HBA_NOT_ASSESSED));
}

#[test]
fn a_part_of_the_post_pass_that_fails_stays_in_its_section() {
    // An error and a panic alike come back as the section's key and why.
    assert_eq!(group_c::contain("cnvpytor", || Ok(())), None);
    assert_eq!(
        group_c::contain("expansions", || Err("no table".into())),
        Some(("expansions", "no table".to_owned()))
    );
    let panicked = group_c::contain("locityper", || panic!("planted"));
    assert_eq!(
        panicked,
        Some((
            "locityper",
            "internal error: the post pass of locityper stopped".to_owned()
        ))
    );
    // The section is then unreadable, with no data, and names its file.
    let mut s = Section::new(
        GROUP_C
            .iter()
            .copied()
            .find(|s| s.key == "locityper")
            .unwrap(),
        Ok(Some((
            "locityper/S_locityper.tsv".into(),
            1_600_000_000.0,
            Data::Json(Value::Object(indexmap::IndexMap::new())),
        ))),
        None,
    );
    s.make_unreadable("internal error: the post pass of locityper stopped");
    let j = s.to_json();
    assert_eq!(word(at(&j, &["state"])), "unreadable");
    assert_eq!(
        word(at(&j, &["note"])),
        "S_locityper.tsv: internal error: the post pass of locityper stopped"
    );
    assert_eq!(at(&j, &["source"]), &Value::Null);
    assert_eq!(at(&j, &["data"]), &Value::Object(indexmap::IndexMap::new()));
}

#[test]
fn an_infinite_telomere_read_length_keeps_the_summary() {
    // The Python collector stopped on this table with an OverflowError and
    // wrote no summary.json; here only the telomere section is unreadable.
    let d = Dir::new("post-inf");
    d.put(
        "S/telomere/S/S/S_summary.tsv",
        "PID\tsample\ttotal_reads\tread_lengths\tintratel_reads\ttel_content\nS\ttumor\t1000\tinf\t40\t1.5\n",
    );
    d.put(
        "S/mosdepth/S.mosdepth.summary.txt",
        "chrom\tlength\tbases\tmean\tmin\tmax\ntotal\t100\t3000\t30\t0\t60\n",
    );
    let v = summary(&d);
    assert_eq!(
        word(at(&v, &["sections", "telomere", "state"])),
        "unreadable"
    );
    assert_eq!(word(at(&v, &["sections", "coverage", "state"])), "ok");
}

#[test]
fn a_declared_female_drops_cnvpytor_chry_calls_across_groups() {
    // The declared sex comes from step 33's table (group a) and changes
    // CNVpytor's count (group c).
    let d = Dir::new("post-sex");
    d.put(
        "S/qc/S_sample_qc.tsv",
        "key\tvalue\nsample\tS\ndeclared_sex\tfemale\ninferred_sex\tfemale\nsex_check\tok\n",
    );
    d.put(
        "S/cnvpytor/S_cnvs.txt",
        "deletion\tchrY:1-200000\t200000\t0.4\t0.00001\t1\t1\t1\t0.00001\t0.00001\t0.1\n\
deletion\tchr2:1-200000\t200000\t0.4\t0.00001\t1\t1\t1\t0.00001\t0.00001\t0.1\n",
    );
    let v = summary(&d);
    assert_eq!(word(at(&v, &["run", "declared_sex"])), "female");
    let cnv = at(&v, &["sections", "cnvpytor"]);
    assert_eq!(word(at(cnv, &["state"])), "ok");
    assert_eq!(word(at(cnv, &["data", "declared_sex"])), "female");
    assert!(has(at(cnv, &["data", "dropped"]), "chry_female"));
}

const PREV: &str = r#"{"schema_version": 2, "generated_utc": "2020-01-01 00:00 UTC",
 "sections": {"coverage": {"state": "ok", "data": {}}, "clinvar": {"state": "missing", "data": {}}}}"#;

#[test]
fn new_since_compares_with_the_newest_kept_run() {
    let d = Dir::new("post-history");
    d.put("S/history/20200101T000000/summary.json", PREV);
    // A newer folder without a summary.json is a half-copied one: passed over.
    d.put("S/history/20990101T000000/run_manifest.tsv", "run\tx\ty\n");
    let v = summary(&d);
    assert_eq!(
        at(&v, &["new_since", "previous_utc"]),
        &Value::Str("2020-01-01 00:00 UTC".into())
    );
    assert!(texts(at(&v, &["new_since", "sections"])).contains(&"coverage".to_owned()));
}

#[test]
fn a_previous_run_that_cannot_be_read_leaves_new_since_out() {
    for body in ["{not json", "[]", r#"{"sections": 7}"#] {
        let d = Dir::new("post-history-bad");
        d.put("S/history/20200101T000000/summary.json", body);
        let v = summary(&d);
        assert!(!has(&v, "new_since"), "{body}");
    }
}

#[test]
fn the_recheck_block_comes_from_its_state_file() {
    let d = Dir::new("post-recheck");
    d.put(
        "S/recheck/state.json",
        r#"{"checked_utc": "2026-01-02 03:04:05 UTC", "items": [{"item": "clinvar", "class": "newer", "x": 1}, 5],
            "affected_steps": ["06"]}"#,
    );
    let v = summary(&d);
    assert_eq!(
        word(at(&v, &["recheck", "checked_utc"])),
        "2026-01-02 03:04:05 UTC"
    );
    assert_eq!(word(at(&v, &["recheck", "items", "0", "item"])), "clinvar");
    // A check older than the run's manifest is of an earlier run.
    d.put(
        "S/run_manifest.tsv",
        "run\twritten_utc\t2027-01-01 00:00:00 UTC\n",
    );
    assert!(!has(&summary(&d), "recheck"));
}

#[test]
fn without_only_the_command_writes_summary_json_in_the_folder() {
    let d = Dir::new("post-cli");
    std::fs::create_dir_all(d.path("S")).unwrap();
    let dir = d.path("S");
    let assets = d.path("assets");
    let (code, out, err) = run(
        crate::collect::cli,
        &["--sample", "S", "--sample-dir", &dir, "--assets", &assets],
    );
    assert_eq!((code, err.as_str()), (0, ""));
    let target = format!("{dir}/summary.json");
    assert_eq!(out, format!("Summary written: {target}\n"));
    let v = loads(&d.read("S/summary.json")).unwrap();
    assert!(has(&v, "not_assessed"));
    // A part of the sections never takes summary.json's place.
    let (code, _, err) = run(
        crate::collect::cli,
        &["--sample", "S", "--sample-dir", &dir, "--only", "a"],
    );
    assert_eq!(code, 2);
    assert!(err.contains("--only needs --out"), "{err}");
}
