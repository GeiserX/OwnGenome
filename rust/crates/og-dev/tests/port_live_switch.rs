//! `og dev steps-check` on this checkout, and on planted mistakes. The
//! checks moved here from the Python test of the run state reader, which
//! went with the Python live server; each must fail on its planted mistake,
//! because a check that cannot fail is no check.

use og_dev::steps_check::{Code, Outcome, Row, check, checks, cli, read_code, read_rows};
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::OnceLock;

fn repo() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../..")
        .canonicalize()
        .unwrap()
}

fn real() -> &'static (Vec<Row>, Code) {
    static R: OnceLock<(Vec<Row>, Code)> = OnceLock::new();
    R.get_or_init(|| {
        (
            read_rows(&repo().join("assets/steps.tsv")).unwrap(),
            read_code(&repo()).unwrap(),
        )
    })
}

/// The names of the checks that fail.
fn failing(out: &[Outcome]) -> Vec<&'static str> {
    out.iter()
        .filter(|(_, r)| r.is_err())
        .map(|(n, _)| *n)
        .collect()
}

const IDLE: &str = "every row is a step run-all.sh runs or a process of main.nf";

#[test]
fn this_checkout_passes_every_check() {
    let out = check(&repo()).unwrap();
    assert_eq!(out.len(), 10, "{out:?}");
    assert_eq!(failing(&out), Vec::<&str>::new(), "{out:?}");
}

/// The row a plan adds by mistake: a step number that nothing runs.
fn planted_row() -> Row {
    Row {
        step: "99".into(),
        label: "Planted".into(),
        tool: String::new(),
        recorded: false,
        processes: vec![],
        sections: vec![],
        upstream: vec![],
        hours: None,
        hours_from: String::new(),
        hours_low: String::new(),
    }
}

#[test]
fn a_row_for_a_step_nothing_runs_fails() {
    let (rows, code) = real();
    let mut planted = rows.clone();
    planted.push(planted_row());
    assert_eq!(failing(&checks(&planted, code)), [IDLE]);
    // The same row with a process of a module is a step that runs.
    let mut with_process = rows.clone();
    with_process.push(Row {
        processes: vec!["NO_SUCH_PROCESS".into()],
        ..planted_row()
    });
    assert!(!failing(&checks(&with_process, code)).contains(&IDLE));
}

#[test]
fn each_check_fails_on_its_planted_mistake() {
    let (rows, code) = real();
    let edit = |f: &dyn Fn(&mut Row)| -> Vec<Row> {
        rows.iter()
            .cloned()
            .map(|mut r| {
                f(&mut r);
                r
            })
            .collect()
    };
    let cases: Vec<(&str, Vec<Row>, Code)> = vec![
        (
            "every run-all.sh plan step has its row, label and tool",
            rows.iter().filter(|r| r.step != "49").cloned().collect(),
            code.clone(),
        ),
        (
            "every run-all.sh plan step has its row, label and tool",
            edit(&|r| {
                if r.step == "07" {
                    r.label = "Renamed".into();
                }
            }),
            code.clone(),
        ),
        (
            "the recorded steps are what run-all.sh records",
            edit(&|r| {
                if r.step == "04b" {
                    r.recorded = false;
                }
            }),
            code.clone(),
        ),
        (
            "the summary's recorded steps are the recorded rows",
            rows.clone(),
            Code {
                run_all_steps: code
                    .run_all_steps
                    .iter()
                    .filter(|s| *s != "33")
                    .cloned()
                    .collect(),
                ..code.clone()
            },
        ),
        (
            "every summary section belongs to the row of its step",
            edit(&|r| r.sections.retain(|k| k != "mosaic_sex_loss")),
            code.clone(),
        ),
        (
            "every summary section belongs to the row of its step",
            edit(&|r| {
                if r.step == "11" {
                    r.sections.push("no_such_section".into());
                }
            }),
            code.clone(),
        ),
        (
            "every module process belongs to one step",
            edit(&|r| r.processes.retain(|p| p != "Y_STR")),
            code.clone(),
        ),
        (
            "every module process belongs to one step",
            edit(&|r| {
                if r.step == "45" {
                    r.processes.push("ROH".into());
                }
            }),
            code.clone(),
        ),
        (
            "the tools are main.nf's known tools",
            edit(&|r| {
                if r.step == "34" {
                    r.tool.clear();
                }
            }),
            code.clone(),
        ),
        (
            "CPSR's read check counts under step 17, ClinVar's under step 06",
            edit(&|r| r.processes.retain(|p| p != "CPSR_CDNA_CHECK")),
            code.clone(),
        ),
        (
            "CPSR's read check counts under step 17, ClinVar's under step 06",
            rows.clone(),
            Code {
                clinical: code.clinical.replace(
                    "include { CLINVAR_CDNA_CHECK as CPSR_CDNA_CHECK }",
                    "include { CLINVAR_CDNA_CHECK }",
                ),
                ..code.clone()
            },
        ),
        (
            "one row per step, every upstream step has a row",
            rows.iter().cloned().chain([rows[0].clone()]).collect(),
            code.clone(),
        ),
        (
            "one row per step, every upstream step has a row",
            edit(&|r| {
                if r.step == "11" {
                    r.upstream.push("98".into());
                }
            }),
            code.clone(),
        ),
        (
            "every measured row has hours_low below hours, no other row has one",
            edit(&|r| {
                if r.step == "03" {
                    r.hours_low.clear();
                }
            }),
            code.clone(),
        ),
        (
            "every measured row has hours_low below hours, no other row has one",
            edit(&|r| {
                if r.step == "04" {
                    r.hours_low = format!("{}", r.hours.unwrap());
                }
            }),
            code.clone(),
        ),
        (
            "every measured row has hours_low below hours, no other row has one",
            edit(&|r| {
                if r.hours_from != "measured" && r.step == "14" {
                    r.hours_low = "0.1".into();
                }
            }),
            code.clone(),
        ),
    ];
    for (name, rows, code) in cases {
        let got = failing(&checks(&rows, &code));
        assert!(
            got.contains(&name),
            "the planted mistake for '{name}' passed: {got:?}"
        );
    }
    // Control: the unchanged catalogue passes, so the cases test the mistake.
    assert_eq!(failing(&checks(rows, code)), Vec::<&str>::new());
}

/// A checkout that is this one but for a planted row in assets/steps.tsv:
/// links to every other file.
fn planted_checkout(dir: &Path) {
    let r = repo();
    fs::create_dir_all(dir.join("assets")).unwrap();
    for e in fs::read_dir(&r).unwrap() {
        let e = e.unwrap();
        if e.file_name() != "assets" {
            std::os::unix::fs::symlink(e.path(), dir.join(e.file_name())).unwrap();
        }
    }
    for e in fs::read_dir(r.join("assets")).unwrap() {
        let e = e.unwrap();
        if e.file_name() != "steps.tsv" {
            std::os::unix::fs::symlink(e.path(), dir.join("assets").join(e.file_name())).unwrap();
        }
    }
    let mut steps = fs::read_to_string(r.join("assets/steps.tsv")).unwrap();
    steps.push_str("99\tPlanted\t\tno\t\t\t\t\t\t\t\n");
    fs::write(dir.join("assets/steps.tsv"), steps).unwrap();
}

#[test]
fn the_command_passes_this_checkout_and_fails_a_planted_row() {
    let run = |root: &Path| {
        let (mut out, mut err) = (Vec::new(), Vec::new());
        let args = ["--root".to_owned(), root.to_string_lossy().into_owned()];
        let code = cli(&args, &mut out, &mut err);
        (
            code,
            String::from_utf8(out).unwrap(),
            String::from_utf8(err).unwrap(),
        )
    };
    let (code, out, err) = run(&repo());
    assert_eq!(code, 0, "{out}{err}");
    assert_eq!(
        out.lines().filter(|l| l.starts_with("ok ")).count(),
        10,
        "{out}"
    );

    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"))
        .join(format!("steps-check-planted-{}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    planted_checkout(&dir);
    let (code, out, err) = run(&dir);
    let _ = fs::remove_dir_all(&dir);
    assert_eq!(code, 1, "{out}{err}");
    let fails: Vec<&str> = out.lines().filter(|l| l.starts_with("FAIL")).collect();
    assert_eq!(fails.len(), 1, "{out}");
    assert!(
        fails[0].contains(IDLE) && fails[0].contains("\"99\""),
        "{out}"
    );

    let (code, _, err) = run(Path::new("/nonexistent-og-steps-check"));
    assert_eq!(code, 2, "{err}");
    let (mut o, mut e) = (Vec::new(), Vec::new());
    assert_eq!(cli(&["--nope".into()], &mut o, &mut e), 2);
}
