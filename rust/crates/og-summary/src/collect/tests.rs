use super::model::{RUN_ALL_STEPS, RunStatus, State, is_stale};
use super::sections::{Variants, karyotype};
use super::*;
use crate::testutil::{Dir, run};
use pycompat::json::{Value, loads};

const CATALOG: &str = "# a comment line\ngene\thgnc_id\tchrom\tstart\tend\treasons\towner\n\
A1\tHGNC:1\tchr1\t0\t100\tsegdup\t.\n\
B2\tHGNC:2\tchr1\t200\t300\tCMRG,segdup\tparascopy\n\
C3\tHGNC:3\tchrY\t0\t10\tsegdup\t.\n\
D4\tHGNC:4\tchr2_alt\t.\t.\tALT-only\t.\n";

fn assets(d: &Dir) -> String {
    d.put("assets/callability_genes.tsv", CATALOG);
    d.path("assets")
}

fn summary(d: &Dir, sample_dir: &str) -> Value {
    let input = Input {
        sample: "S".into(),
        sample_dir: d.path(sample_dir),
        declared_sex: None,
        assets: d.path("assets"),
    };
    collect(&input).unwrap().to_json()
}

fn section<'a>(v: &'a Value, key: &str) -> &'a Value {
    let Value::Object(m) = v else { panic!() };
    let Value::Object(s) = &m["sections"] else {
        panic!()
    };
    &s[key]
}

fn get<'a>(v: &'a Value, path: &str) -> &'a Value {
    path.split('.').fold(v, |v, k| match v {
        Value::Object(m) => m.get(k).unwrap_or_else(|| panic!("no {k} in {path}")),
        _ => panic!("{path}"),
    })
}

fn s(x: &str) -> Value {
    Value::Str(x.into())
}

fn set_mtime(path: &str, secs: i64) {
    let f = std::fs::File::options().write(true).open(path).unwrap();
    f.set_modified(std::time::UNIX_EPOCH + std::time::Duration::from_secs(secs as u64))
        .unwrap();
}

#[test]
fn every_section_is_missing_in_an_empty_folder() {
    let d = Dir::new("collect-empty");
    std::fs::create_dir_all(d.path("S")).unwrap();
    let v = summary(&d, "S");
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
            "sections"
        ]
    );
    assert_eq!(get(&v, "run.status_file"), &Value::Bool(false));
    assert_eq!(get(&v, "manifest"), &Value::Array(vec![]));
    for spec in GROUP_A {
        let sec = section(&v, spec.key);
        assert_eq!(get(sec, "state"), &s("missing"), "{}", spec.key);
        assert_eq!(get(sec, "data"), &Value::Object(Default::default()));
    }
}

#[test]
fn a_file_older_than_the_run_is_stale_unless_its_step_was_ok() {
    let d = Dir::new("collect-stale");
    let started = 1_577_836_800; // 2020-01-01 00:00 UTC
    d.put(
        "S/logs/run_status.tsv",
        format!("meta\tstarted_epoch\t{started}\nstep\t16b\tskipped\nstep\t33\tok\n"),
    );
    let old = d.put("S/qc/S_sample_qc.tsv", "key\tvalue\ninferred_sex\tmale\n");
    set_mtime(&old, started - 3600);
    let cov = d.put(
        "S/mosdepth/S.mosdepth.summary.txt",
        "chrom\tlength\tbases\tmean\ntotal\t1\t1\t31.25\n",
    );
    set_mtime(&cov, started - 60);
    let vcf = d.put_gz("S/vcf/S.vcf.gz", "#h\n");
    set_mtime(&vcf, started - 60);
    let v = summary(&d, "S");
    let c = section(&v, "coverage");
    assert_eq!(get(c, "state"), &s("stale"));
    assert_eq!(
        get(c, "note"),
        &s(
            "from an earlier run (2019-12-31 23:59 UTC); step 16b was skipped in the run of 2020-01-01 00:00 UTC"
        )
    );
    // The data of a stale section is kept.
    assert_eq!(get(c, "data.mean_depth"), &Value::Float(31.25));
    // Step 33 was ok in the run: its older table is current.
    assert_eq!(get(section(&v, "sample_qc"), "state"), &s("ok"));
    // Step 03 is not a run-all step: never stale.
    assert_eq!(get(section(&v, "variants"), "state"), &s("ok"));
}

#[test]
fn a_file_written_in_the_second_the_run_started_is_current() {
    // The planted fault `summary-stale-off-by-one` (`<=` for `<`) turns
    // this red.
    let started = 1_577_836_800.0;
    assert!(!is_stale(started, started));
    assert!(is_stale(started - 0.5, started));
    let d = Dir::new("collect-same-second");
    d.put(
        "S/logs/run_status.tsv",
        "meta\tstarted_epoch\t1577836800\nstep\t16b\tskipped\n",
    );
    let cov = d.put(
        "S/coverage/S.mosdepth.summary.txt",
        "chrom\tmean\ntotal\t30\n",
    );
    set_mtime(&cov, 1_577_836_800);
    let v = summary(&d, "S");
    assert_eq!(get(section(&v, "coverage"), "state"), &s("ok"));
    assert_eq!(
        get(section(&v, "coverage"), "source"),
        &s("coverage/S.mosdepth.summary.txt")
    );
    assert_eq!(
        get(section(&v, "coverage"), "file_date"),
        &s("2020-01-01 00:00 UTC")
    );
}

#[test]
fn a_missing_file_of_a_failed_step_is_failed() {
    let d = Dir::new("collect-failed");
    d.put(
        "S/logs/run_status.tsv",
        "step\t33\tfailed (exit 1)\nmeta\tstarted_utc\tsome day\n",
    );
    let v = summary(&d, "S");
    let q = section(&v, "sample_qc");
    assert_eq!(get(q, "state"), &s("failed"));
    assert_eq!(
        get(q, "note"),
        &s("step 33 failed in the run of some day; its log is in logs/")
    );
    assert_eq!(get(&v, "run.started_utc"), &s("some day"));
    assert_eq!(get(section(&v, "coverage"), "state"), &s("missing"));
}

#[test]
fn an_empty_start_time_is_an_unknown_date() {
    let d = Dir::new("collect-failed-empty-date");
    d.put(
        "S/logs/run_status.tsv",
        "meta\tstarted_utc\t\nstep\t33\tfailed (exit 1)\n",
    );
    let v = summary(&d, "S");
    assert_eq!(
        get(section(&v, "sample_qc"), "note"),
        &s("step 33 failed in the run of unknown date; its log is in logs/")
    );
}

#[test]
fn the_run_status_rows() {
    let st = RunStatus::parse(
        "# comment\nmeta\tstarted_epoch\t100\nmeta\tstarted_epoch\tsoon\nmeta\tdeclared_sex\tfemale\nstep\t16\tok\nstep\t33\n",
    );
    assert_eq!(st.started_epoch, Some(100.0));
    assert_eq!(st.started_utc.as_deref(), Some("1970-01-01 00:01 UTC"));
    assert_eq!(st.declared_sex.as_deref(), Some("female"));
    assert_eq!(st.steps.get("16").map(String::as_str), Some("ok"));
    assert!(!st.steps.contains_key("33"));
    // A start time that is no date is no start time.
    assert_eq!(
        RunStatus::parse("meta\tstarted_epoch\tnan\n").started_epoch,
        None
    );
}

#[test]
fn the_declared_sex_comes_from_the_first_place_that_names_one() {
    let d = Dir::new("collect-declared");
    d.put(
        "S/run_manifest.tsv",
        "run\tdeclared_sex\tmale\ndata\tclinvar_file_date\t2026-01-01\tx\n",
    );
    d.put(
        "S/qc/S_sample_qc.tsv",
        "key\tvalue\ninferred_sex\tmale\ndeclared_sex\tfemale\n",
    );
    let v = summary(&d, "S");
    assert_eq!(get(&v, "run.declared_sex"), &s("male"));
    assert_eq!(get(&v, "databases.clinvar_release"), &s("2026-01-01\tx"));
    std::fs::remove_file(d.path("S/run_manifest.tsv")).unwrap();
    let v = summary(&d, "S");
    assert_eq!(get(&v, "run.declared_sex"), &s("female"));
}

#[test]
fn callability_against_the_catalog() {
    let d = Dir::new("collect-callability");
    assets(&d);
    d.put_gz(
        "S/mosdepth/S_callability.bed.gz",
        "#chrom\tstart\tend\tregion\t10X\nchr1\t0\t100\tA1\t95\nchr1\t200\t300\tB2\t50\n\n",
    );
    let v = summary(&d, "S");
    let c = section(&v, "callability");
    assert_eq!(get(c, "state"), &s("ok"), "{c:?}");
    assert_eq!(get(c, "data.min_depth"), &Value::Int(10));
    assert_eq!(get(c, "data.catalog_genes"), &Value::Int(4));
    assert_eq!(get(c, "data.read"), &Value::Int(1));
    let want = loads(
        r#"{"not_read": [{"gene": "C3", "hgnc_id": "HGNC:3", "reasons": ["segdup"], "share": null, "status": "not on this reference"},
                         {"gene": "D4", "hgnc_id": "HGNC:4", "reasons": ["ALT-only"], "share": null, "status": "not on this reference"}],
            "not_read_owned": [{"gene": "B2", "owner": "parascopy", "share": 0.5, "status": "not read"}],
            "not_on_reference": ["C3", "D4"]}"#,
    )
    .unwrap();
    for k in ["not_read", "not_read_owned", "not_on_reference"] {
        assert_eq!(get(c, &format!("data.{k}")), get(&want, k), "{k}");
    }
}

#[test]
fn callability_errors_name_the_file() {
    let d = Dir::new("collect-callability-bad");
    assets(&d);
    let cases = [
        (
            "chrom\tstart\tend\tregion\t10X\n",
            "S_callability.bed.gz is not mosdepth's thresholds table of one depth",
        ),
        (
            "#chrom\tstart\tend\tregion\t10X\nchr1\t0\t100\tA1\n",
            "S_callability.bed.gz: a row has 4 columns, not 5",
        ),
        (
            "#chrom\tstart\tend\tregion\t10X\n",
            "S_callability.bed.gz has no gene rows",
        ),
        (
            "#chrom\tstart\tend\tregion\t10X\nchr1\t0\t100\tA1\t95\n",
            "S_callability.bed.gz has no row for B2: it was written from another assets/callability_genes.tsv; run step 16b again",
        ),
        (
            "#chrom\tstart\tend\tregion\t10X\nchr1\tx\t100\tA1\t95\n",
            "ValueError: invalid literal for int() with base 10: 'x'",
        ),
    ];
    for (table, note) in cases {
        d.put_gz("S/mosdepth/S_callability.bed.gz", table);
        let v = summary(&d, "S");
        let c = section(&v, "callability");
        assert_eq!(get(c, "state"), &s("unreadable"));
        assert_eq!(get(c, "note"), &s(note));
    }
    std::fs::remove_dir_all(d.path("assets")).unwrap();
    let v = summary(&d, "S");
    let note = get(section(&v, "callability"), "note");
    assert_eq!(
        note,
        &s(&format!(
            "callability_genes.tsv is not in {}: the report needs the repository's assets/ folder (render_report.py --assets)",
            d.path("assets")
        ))
    );
}

#[test]
fn karyotypes_from_whole_copy_numbers_only() {
    assert_eq!(karyotype(Some(1.0), Some(0.0)), Some("45,X"));
    assert_eq!(karyotype(Some(2.2), Some(0.8)), Some("47,XXY"));
    assert_eq!(karyotype(Some(2.0), Some(0.0)), None);
    assert_eq!(karyotype(Some(1.0), Some(1.0)), None);
    // Half to even, then 0.5 from the whole number: none.
    assert_eq!(karyotype(Some(1.5), Some(0.0)), None);
    assert_eq!(karyotype(Some(6.0), Some(0.0)), None);
    assert_eq!(karyotype(Some(f64::INFINITY), Some(0.0)), None);
    assert_eq!(karyotype(None, Some(0.0)), None);
}

#[test]
fn sex_check_reads_the_last_sample_row() {
    let d = Dir::new("collect-sex");
    d.put(
        "S/indexcov/S-indexcov.ped",
        "#family_id\tsample_id\tsex\tCNchrX\tCNchrY\r\nf\tA\t2\t2.0\t0.0\r\n\r\nf\tS\t1\t1.02\t0.04\r\n",
    );
    let v = summary(&d, "S");
    let c = section(&v, "sex_check");
    assert_eq!(get(c, "source"), &s("indexcov/S-indexcov.ped"));
    assert_eq!(get(c, "data.inferred_sex"), &s("male"));
    assert_eq!(get(c, "data.cn_chrX"), &s("1.02"));
    assert_eq!(get(c, "data.karyotype"), &s("45,X"));
    d.put(
        "S/indexcov/indexcov-indexcov.ped",
        "#family_id\tsample_id\n",
    );
    let v = summary(&d, "S");
    assert_eq!(
        get(section(&v, "sex_check"), "note"),
        &s("indexcov-indexcov.ped has no sex column or no sample row")
    );
}

#[test]
fn sample_qc_and_mosaic_from_one_table() {
    let d = Dir::new("collect-qc");
    d.put(
        "S/qc/S_sample_qc.tsv",
        "key\tvalue\nsample\tS\ninferred_sex\tmale\nx_depth_ratio\t0.51\ny_depth_ratio\t0.01\nx_sites\t400\nshort\n\
mosaic_call\tdetected\nmosaic_depth_ratio\t0.8\nmosaic_floor\tinf\nmosaic_x_het_sites\t12\n",
    );
    let v = summary(&d, "S");
    let q = section(&v, "sample_qc");
    assert_eq!(get(q, "data.short"), &Value::Null);
    assert_eq!(get(q, "data.karyotype"), &s("45,X"));
    let m = section(&v, "mosaic_sex_loss");
    assert_eq!(get(m, "state"), &s("ok"));
    let Value::Object(data) = get(m, "data") else {
        panic!()
    };
    let keys: Vec<&str> = data.keys().map(String::as_str).collect();
    assert_eq!(&keys[..4], ["call", "depth_ratio", "floor", "x_het_sites"]);
    assert_eq!(data["depth_ratio"], Value::Float(0.8));
    assert_eq!(data["floor"], Value::Null);
    assert_eq!(data["x_het_sites"], Value::Int(12));
    assert_eq!(data["loss_af"], Value::Null);
    // Without mosaic rows the section is missing; with a call it does not
    // know, unreadable.
    d.put("S/qc/S_sample_qc.tsv", "key\tvalue\ninferred_sex\tmale\n");
    assert_eq!(
        get(section(&summary(&d, "S"), "mosaic_sex_loss"), "state"),
        &s("missing")
    );
    d.put(
        "S/qc/S_sample_qc.tsv",
        "key\tvalue\ninferred_sex\tmale\nmosaic_call\tmaybe\n",
    );
    assert_eq!(
        get(section(&summary(&d, "S"), "mosaic_sex_loss"), "note"),
        &s("S_sample_qc.tsv has no mosaic_call row (detected, not_detected, gain, not_estimated)")
    );
    d.put("S/qc/S_sample_qc.tsv", "key\tvalue\nsex_check\tok\n");
    assert_eq!(
        get(section(&summary(&d, "S"), "sample_qc"), "note"),
        &s("S_sample_qc.tsv has no inferred_sex row")
    );
}

#[test]
fn read_count_needs_whole_numbers_and_a_check() {
    let d = Dir::new("collect-reads");
    d.put(
        "S/qc/S_read_count.tsv",
        "key\tvalue\ninput_reads\t100\nprimary_records\t99\nmissing_reads\t-1\nread_count_check\tmismatch\n",
    );
    let v = summary(&d, "S");
    assert_eq!(
        get(section(&v, "read_count"), "data.missing_reads"),
        &Value::Int(-1)
    );
    d.put("S/qc/S_read_count.tsv", "key\tvalue\ninput_reads\t1.5\n");
    assert_eq!(
        get(section(&summary(&d, "S"), "read_count"), "note"),
        &s("S_read_count.tsv has no whole number for input_reads")
    );
    d.put(
        "S/qc/S_read_count.tsv",
        "key\tvalue\ninput_reads\t1\nprimary_records\t1\nmissing_reads\t0\n",
    );
    assert_eq!(
        get(section(&summary(&d, "S"), "read_count"), "note"),
        &s("S_read_count.tsv has no read_count_check row (ok or mismatch)")
    );
}

#[test]
fn variant_counts() {
    let d = Dir::new("collect-variants");
    let vcf = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\n\
chr1\t1\t.\tA\tG\t50\tPASS\n\
chr1\t2\t.\tA\tAT,<NON_REF>\t50\tLowQual\tinfo\n\
chr1\t3\t.\tAT\tCG\t50\tPASS\r\n\
chr1\t4\t.\tA\t<DEL>\t50\tPASS\n\
chr1\t5\t.\tA\t*,C\t50\tPASS\n\
chr1\t6\t.\tA\tG\n\
chr1\t7\t.\tA\t\u{e9}\t50\tPASS\rchr1\t8\t.\tA\tGT,C\t50\t.\n";
    let p = d.put_gz("v.vcf.gz", vcf);
    assert_eq!(
        Variants::read(&p).unwrap(),
        Variants {
            total: 7,
            pass: 5,
            snps: 4,
            indels: 2
        }
    );
    // Two gzip members, as BGZF writes them.
    let mut two = std::fs::read(&p).unwrap();
    two.extend(std::fs::read(&p).unwrap());
    let p2 = d.put("two.vcf.gz", two);
    assert_eq!(Variants::read(&p2).unwrap().total, 14);
    let full = std::fs::read(&p).unwrap();
    let cut = d.put("cut.vcf.gz", &full[..full.len() - 12]);
    assert_eq!(
        Variants::read(&cut).unwrap_err(),
        "cut.vcf.gz: Compressed file ended before the end-of-stream marker was reached"
    );
}

#[test]
fn the_command_line() {
    let d = Dir::new("collect-cli");
    std::fs::create_dir_all(d.path("S")).unwrap();
    let out = d.path("out.json");
    let (code, stdout, _) = run(
        cli,
        &[
            "--sample",
            "S",
            "--sample-dir",
            &d.path("S"),
            "--out",
            &out,
            "--only",
            "a",
        ],
    );
    assert_eq!(code, 0);
    assert_eq!(stdout, format!("Summary written: {out}\n"));
    let text = d.read("out.json");
    assert!(text.starts_with("{\n \"schema_version\": 2,\n \"sample\": \"S\",\n"));
    assert!(text.ends_with("}\n"));
    // --only is required until every section is ported.
    let (code, _, err) = run(
        cli,
        &["--sample", "S", "--sample-dir", &d.path("S"), "--out", &out],
    );
    assert_eq!(code, 2);
    assert!(err.contains("required: --only"), "{err}");
    // So is --out: a part of the sections must not land on the sample's
    // complete summary.json.
    std::fs::write(d.path("S/summary.json"), "complete\n").unwrap();
    let (code, _, err) = run(
        cli,
        &["--sample", "S", "--sample-dir", &d.path("S"), "--only", "a"],
    );
    assert_eq!(code, 2);
    assert!(err.contains("required: --out"), "{err}");
    assert_eq!(d.read("S/summary.json"), "complete\n");
    let (code, _, err) = run(
        cli,
        &[
            "--sample",
            "S",
            "--sample-dir",
            "x",
            "--out",
            &d.path("no/such/dir/x.json"),
            "--only",
            "a",
        ],
    );
    assert_eq!(code, 1);
    assert!(err.starts_with("ERROR: cannot write"), "{err}");
}

/// A run status that is there and cannot be read is no reason to call an
/// older result current: the collector stops and writes nothing.
#[cfg(unix)]
#[test]
fn an_unreadable_run_status_stops_the_collector() {
    use std::os::unix::fs::PermissionsExt;
    let d = Dir::new("collect-unreadable-status");
    std::fs::create_dir_all(d.path("S/logs")).unwrap();
    let status = d.path("S/logs/run_status.tsv");
    std::fs::write(&status, "meta\tstarted_epoch\t1700000000\n").unwrap();
    std::fs::set_permissions(&status, std::fs::Permissions::from_mode(0o000)).unwrap();
    if std::fs::read(&status).is_ok() {
        return; // root reads every file: nothing to prove here
    }
    let out = d.path("out.json");
    let (code, stdout, err) = run(
        cli,
        &[
            "--sample",
            "S",
            "--sample-dir",
            &d.path("S"),
            "--out",
            &out,
            "--only",
            "a",
        ],
    );
    std::fs::set_permissions(&status, std::fs::Permissions::from_mode(0o644)).unwrap();
    assert_eq!(code, 1);
    assert_eq!(stdout, "");
    assert!(err.contains("run_status.tsv"), "{err}");
    assert!(!std::path::Path::new(&out).exists());
}

#[test]
fn states_have_their_words() {
    let w: Vec<&str> = [
        State::Ok,
        State::Missing,
        State::Failed,
        State::Stale,
        State::Unreadable,
    ]
    .iter()
    .map(|s| s.word())
    .collect();
    assert_eq!(w, ["ok", "missing", "failed", "stale", "unreadable"]);
}

#[test]
fn run_all_steps_are_the_recorded_rows_of_steps_tsv() {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../../assets/steps.tsv");
    let text = std::fs::read_to_string(p).unwrap();
    let mut lines = text
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty());
    let head: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let (step, rec) = (
        head.iter().position(|h| *h == "step").unwrap(),
        head.iter().position(|h| *h == "recorded").unwrap(),
    );
    let mut want: Vec<&str> = lines
        .map(|l| l.split('\t').collect::<Vec<_>>())
        .filter(|c| c.get(rec) == Some(&"yes"))
        .map(|c| c[step])
        .collect();
    want.sort_unstable();
    let mut got = RUN_ALL_STEPS.to_vec();
    got.sort_unstable();
    assert_eq!(got, want);
}
