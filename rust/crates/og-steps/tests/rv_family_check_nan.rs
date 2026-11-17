//! `og family check` on a pair somalier gave no measurement for.
//!
//! somalier relate writes `nan` (or `inf`) as the relatedness of a pair it
//! had nothing to divide by. Every comparison with NaN is false, so the
//! degree bins read it as unrelated, and a pair the pedigree declares
//! unrelated or distant came out `ok`. Such a pair, and one with counts
//! that cannot be counts (negative, or more opposite homozygotes than
//! sites), is `not_checked` with the reason and a `WARNING:` line, and the
//! step still exits 0. A finite relatedness keeps its verdict.
//!
//! The case that reads stderr runs `og family` in a child process: this
//! test binary itself, with [`CHILD`] set (see [`child_entry`]).

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

const CHILD: &str = "OG_TEST_RV_FAMILY_CHECK_NAN_CHILD";

const MEMBERS_HEAD: &str = "sample\tsex\tfather\tmother\tpartner\tin_pedigree\tsequenced\n";
const PAIRS_HEAD: &str = "#sample_a\tsample_b\trelatedness\tibs0\tibs2\tn\n";

/// Two sequenced partners: declared unrelated, expected relatedness 0.
const PARTNERS: &str = "a\tmale\t\t\tb\tyes\tyes\nb\tfemale\t\t\ta\tyes\tyes\n";

/// Two sequenced first cousins (expected 0.125, third degree): the
/// grandparents, their two children and the children's partners are in
/// the pedigree, not sequenced.
const COUSINS: &str = "gf\tmale\t\t\tgm\tyes\tno\n\
gm\tfemale\t\t\tgf\tyes\tno\n\
p1\tmale\tgf\tgm\ts1\tyes\tno\n\
s1\tfemale\t\t\tp1\tyes\tno\n\
p2\tfemale\tgf\tgm\ts2\tyes\tno\n\
s2\tmale\t\t\tp2\tyes\tno\n\
c1\tmale\tp1\ts1\t\tyes\tyes\n\
c2\tfemale\ts2\tp2\t\tyes\tyes\n";

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

fn scratch(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!(
        "og-rv-family-check-nan-{}-{name}",
        std::process::id()
    ));
    let _ = fs::remove_dir_all(&d);
    fs::create_dir_all(&d).unwrap();
    d
}

/// Writes the members and somalier's pairs, runs `og family check` in a
/// child process, and returns (exit code, stderr, the one row of
/// relationships.tsv as (column, value) pairs).
fn check(name: &str, members: &str, pair: &str) -> (i32, String, Vec<(String, String)>) {
    let d = scratch(name);
    let (m, p, out) = (
        d.join("members.tsv"),
        d.join("pairs.tsv"),
        d.join("relationships.tsv"),
    );
    fs::write(&m, format!("{MEMBERS_HEAD}{members}")).unwrap();
    fs::write(&p, format!("{PAIRS_HEAD}{pair}\n")).unwrap();
    let args = [
        "check",
        "--members",
        path(&m),
        "--pairs",
        path(&p),
        "--out",
        path(&out),
    ];
    let res = Command::new(std::env::current_exe().unwrap())
        .args(["child_entry", "--exact", "--nocapture", "--test-threads=1"])
        .env(CHILD, args.join("\u{1f}"))
        .output()
        .unwrap();
    let code = res.status.code().unwrap_or(-1);
    let err = String::from_utf8_lossy(&res.stderr).into_owned();
    let text = fs::read_to_string(&out).unwrap_or_default();
    let mut lines = text.lines();
    let head: Vec<&str> = lines.next().unwrap_or_default().split('\t').collect();
    let rows: Vec<&str> = lines.collect();
    assert_eq!(
        rows.len(),
        1,
        "{name}: one pair, one row: {text}\nstderr: {err}"
    );
    let row = head
        .iter()
        .zip(rows[0].split('\t'))
        .map(|(k, v)| ((*k).to_owned(), v.to_owned()))
        .collect();
    let _ = fs::remove_dir_all(&d);
    (code, err, row)
}

fn path(p: &Path) -> &str {
    p.to_str().unwrap()
}

fn get<'a>(row: &'a [(String, String)], k: &str) -> &'a str {
    row.iter()
        .find(|(c, _)| c == k)
        .map_or_else(|| panic!("no column {k} in {row:?}"), |(_, v)| v.as_str())
}

/// One pair that must be not_checked, with exit 0 and a WARNING line that
/// carries `reason`.
fn assert_not_checked(name: &str, members: &str, pair: &str, reason: &str) {
    let (code, err, row) = check(name, members, pair);
    assert_eq!(code, 0, "{name}: exit code; stderr: {err}");
    assert_eq!(get(&row, "verdict"), "not_checked", "{name}: {row:?}");
    assert_eq!(get(&row, "reason"), reason, "{name}: {row:?}");
    assert!(
        err.lines()
            .any(|l| l.starts_with("WARNING: ") && l.ends_with(reason)),
        "{name}: no WARNING line with the reason; stderr: {err}"
    );
}

#[test]
fn nan_relatedness_of_declared_unrelated_pair_is_not_checked() {
    assert_not_checked(
        "nan-partners",
        PARTNERS,
        "a\tb\tnan\t10\t500\t1000",
        "somalier gave relatedness 'nan' for this pair: not measured",
    );
}

#[test]
fn inf_relatedness_of_declared_unrelated_pair_is_not_checked() {
    assert_not_checked(
        "inf-partners",
        PARTNERS,
        "a\tb\tinf\t10\t500\t1000",
        "somalier gave relatedness 'inf' for this pair: not measured",
    );
    assert_not_checked(
        "minus-inf-partners",
        PARTNERS,
        "a\tb\t-inf\t10\t500\t1000",
        "somalier gave relatedness '-inf' for this pair: not measured",
    );
}

#[test]
fn nan_relatedness_of_declared_third_degree_pair_is_not_checked() {
    let (code, _, row) = check("cousins-ok", COUSINS, "c1\tc2\t0.12\t5\t500\t5000");
    assert_eq!(code, 0);
    assert_eq!(get(&row, "expected_relatedness"), "0.125", "{row:?}");
    assert_eq!(get(&row, "verdict"), "ok", "{row:?}");
    assert_not_checked(
        "nan-cousins",
        COUSINS,
        "c1\tc2\tNaN\t5\t500\t5000",
        "somalier gave relatedness 'NaN' for this pair: not measured",
    );
}

#[test]
fn counts_that_cannot_be_counts_are_not_checked() {
    let reason = |ibs0: i32, ibs2: i32, n: i32| {
        format!(
            "somalier gave ibs0 {ibs0}, ibs2 {ibs2} and {n} sites for this pair: not counts of \
             sites, not measured"
        )
    };
    assert_not_checked(
        "ibs0-over-n",
        PARTNERS,
        "a\tb\t0.01\t2000\t500\t1000",
        &reason(2000, 500, 1000),
    );
    assert_not_checked(
        "negative-ibs2",
        PARTNERS,
        "a\tb\t0.01\t10\t-1\t1000",
        &reason(10, -1, 1000),
    );
    assert_not_checked(
        "negative-n",
        PARTNERS,
        "a\tb\t0.01\t-5\t500\t-10",
        &reason(-5, 500, -10),
    );
}

#[test]
fn finite_relatedness_keeps_its_verdict() {
    let (code, err, row) = check("finite-partners", PARTNERS, "a\tb\t0.01\t600\t100\t5000");
    assert_eq!(code, 0, "stderr: {err}");
    assert_eq!(get(&row, "verdict"), "ok", "{row:?}");
    assert_eq!(get(&row, "reason"), "unrelated, as declared", "{row:?}");
    assert!(!err.contains("WARNING"), "stderr: {err}");
    // Too few sites keeps its own reason, also beside a nan relatedness or
    // more opposite homozygotes than sites (the parity case family/check-few).
    for (name, pair) in [
        ("few-partners", "a\tb\t0.01\t6\t10\t50"),
        ("few-nan-partners", "a\tb\tnan\t6\t10\t50"),
        ("few-ibs0-over-n-partners", "a\tb\t0.01\t600\t10\t50"),
    ] {
        let (code, _, row) = check(name, PARTNERS, pair);
        assert_eq!(code, 0, "{name}");
        assert_eq!(get(&row, "verdict"), "not_checked", "{name}: {row:?}");
        assert!(
            get(&row, "reason").starts_with("somalier compared the two at 50 sites"),
            "{name}: {row:?}"
        );
    }
    // A measured first degree for declared partners still stops the step.
    let (code, _, row) = check("close-partners", PARTNERS, "a\tb\t0.5\t1\t900\t5000");
    assert_eq!(code, 1);
    assert_eq!(get(&row, "verdict"), "mismatch", "{row:?}");
}
