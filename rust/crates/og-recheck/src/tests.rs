//! The check on synthetic sample folders (`tests/fixtures/recheck`), with
//! no network: the answers come from the offline file, or from a server on
//! 127.0.0.1.

use super::*;
use crate::answers::{Offline, Online, PGS_API};
use crate::rerun::{compare, launch_dir_command};
use og_fresh::http::GitHub;
use og_fresh::pins::Pins;
use std::fs;
use std::io::{BufRead, BufReader, Read, Write};
use std::net::TcpListener;
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};

const CLINVAR_STEPS: &[&str] = &["06", "13", "20", "23", "24", "27", "30", "31", "38", "39"];

fn repo() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../..")
        .canonicalize()
        .unwrap()
}

fn fix(name: &str) -> PathBuf {
    repo().join("tests/fixtures/recheck").join(name)
}

fn steps() -> Vec<Step> {
    og_runlog::steps::read_steps(&repo().join("assets/steps.tsv")).unwrap()
}

/// The pins "here" are the run's own: the version rows of the fixture
/// manifest.
fn this_pins() -> IndexMap<String, String> {
    fs::read_to_string(fix("run_manifest.tsv"))
        .unwrap()
        .lines()
        .filter_map(|l| {
            let c: Vec<&str> = l.split('\t').collect();
            (c[0] == "version").then(|| (c[1].to_owned(), c[2].to_owned()))
        })
        .collect()
}

fn tmp(name: &str) -> PathBuf {
    static N: AtomicUsize = AtomicUsize::new(0);
    let d = std::env::temp_dir().join(format!(
        "og-recheck-{}-{name}-{}",
        std::process::id(),
        N.fetch_add(1, Ordering::Relaxed)
    ));
    let _ = fs::remove_dir_all(&d);
    fs::create_dir_all(&d).unwrap();
    d
}

fn gz(path: &Path, text: &str) {
    fs::create_dir_all(path.parent().unwrap()).unwrap();
    let mut e = flate2_free_gzip(text);
    fs::write(path, &mut e).unwrap();
}

/// A gzip file without a gzip writer: one stored deflate block.
fn flate2_free_gzip(text: &str) -> Vec<u8> {
    let data = text.as_bytes();
    assert!(data.len() < 65_535);
    let mut out = vec![0x1f, 0x8b, 8, 0, 0, 0, 0, 0, 0, 255, 1];
    let len = u16::try_from(data.len()).unwrap();
    out.extend(len.to_le_bytes());
    out.extend((!len).to_le_bytes());
    out.extend(data);
    out.extend(crc32(data).to_le_bytes());
    out.extend(u32::try_from(data.len()).unwrap().to_le_bytes());
    out
}

fn crc32(data: &[u8]) -> u32 {
    let mut c = !0u32;
    for &b in data {
        c ^= u32::from(b);
        for _ in 0..8 {
            c = if c & 1 == 1 {
                (c >> 1) ^ 0xEDB8_8320
            } else {
                c >> 1
            };
        }
    }
    !c
}

struct Genome {
    root: PathBuf,
}

impl Genome {
    fn sample(&self) -> PathBuf {
        self.root.join("S1")
    }
    fn manifest(&self) -> PathBuf {
        self.sample().join("run_manifest.tsv")
    }
    fn status(&self) -> PathBuf {
        self.sample().join("logs/run_status.tsv")
    }
    fn add_manifest(&self, rows: &[(&str, &str)]) {
        let mut f = fs::OpenOptions::new()
            .append(true)
            .open(self.manifest())
            .unwrap();
        for (k, v) in rows {
            writeln!(f, "version\t{k}\t{v}").unwrap();
        }
    }
    fn add_steps(&self, rows: &[(&str, &str)]) {
        let mut f = fs::OpenOptions::new()
            .append(true)
            .open(self.status())
            .unwrap();
        for (k, v) in rows {
            writeln!(f, "step\t{k}\t{v}").unwrap();
        }
    }
}

#[derive(Default)]
struct Opt<'a> {
    installed_clinvar: Option<&'a str>,
    running: bool,
    installed_pgs: Option<&'a str>,
    /// The run's ClinVar row as given (`Some("")` drops it).
    run_clinvar: Option<&'a str>,
    command: Option<String>,
}

fn genome(name: &str, o: Opt<'_>) -> Genome {
    let root = tmp(name).join("genome");
    let s = root.join("S1");
    fs::create_dir_all(s.join("logs")).unwrap();
    let mut man = fs::read_to_string(fix("run_manifest.tsv")).unwrap();
    if let Some(rc) = o.run_clinvar {
        man = man
            .lines()
            .filter(|l| !l.starts_with("data\tclinvar_file_date\t"))
            .map(|l| format!("{l}\n"))
            .collect();
        if !rc.is_empty() {
            man +=
                &format!("data\tclinvar_file_date\t{rc} (clinvar/clinvar_pathogenic_chr.vcf.gz)\n");
        }
    }
    fs::write(s.join("run_manifest.tsv"), man).unwrap();
    gz(
        &root.join("clinvar/clinvar_pathogenic_chr.vcf.gz"),
        &format!(
            "##fileformat=VCFv4.1\n##fileDate={}\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
            o.installed_clinvar.unwrap_or("20260901")
        ),
    );
    if let Some(p) = o.installed_pgs {
        gz(
            &root.join("prs_scores/PGS000001.txt.gz"),
            &format!("#HmPOS_date={p}\n"),
        );
    }
    let mut st = format!(
        "meta\tstarted_epoch\t1000\nmeta\tstarted_utc\t2026-09-02 10:00:00 UTC\nmeta\tcommand\t{}\nmeta\tswitches\tGENOME_DIR TOOLS\n",
        o.command.unwrap_or_else(|| "echo replayed".into())
    );
    if !o.running {
        st += "meta\tfinished_utc\t2026-09-02 12:00:00 UTC\n";
    }
    fs::write(s.join("logs/run_status.tsv"), st).unwrap();
    Genome { root }
}

/// The offline answers, with some entries replaced.
fn answers(over: &[(&str, &str)]) -> Offline {
    let mut text = fs::read_to_string(fix("offline.json")).unwrap();
    let mut v = match pycompat::json::loads(&text).unwrap() {
        Value::Object(m) => m,
        _ => unreachable!(),
    };
    for (k, j) in over {
        v.insert((*k).to_owned(), pycompat::json::loads(j).unwrap());
    }
    text = pycompat::json::dumps(&Value::Object(v), &Default::default()).unwrap();
    let p = tmp("answers").join("answers.json");
    fs::write(&p, text).unwrap();
    Offline::read(&p).unwrap()
}

fn run_check(
    g: &Genome,
    a: &dyn Answers,
    pins: &IndexMap<String, String>,
    hla: Option<&str>,
) -> Check {
    check("S1", &g.root, a, &steps(), &Here { pins, hla })
}

fn by_item(c: &Check) -> IndexMap<&str, &Item> {
    c.items.iter().map(|i| (i.item.as_str(), i)).collect()
}

fn cli(args: &[&str]) -> (u8, String, String) {
    let args: Vec<String> = args.iter().map(|a| (*a).to_owned()).collect();
    let (mut o, mut e) = (Vec::new(), Vec::new());
    let rc = cli::run(&args, &mut o, &mut e);
    (
        rc,
        String::from_utf8(o).unwrap(),
        String::from_utf8(e).unwrap(),
    )
}

fn s(p: &Path) -> &str {
    p.to_str().unwrap()
}

#[test]
fn clinvar_newer_upstream_is_refresh_available_with_06_and_its_downstream() {
    let g = genome("cv", Opt::default());
    let c = run_check(&g, &answers(&[]), &this_pins(), Some("3.66.0"));
    let cv = by_item(&c)["clinvar_file_date"];
    assert_eq!(cv.class, Class::RefreshAvailable);
    assert_eq!(cv.steps, CLINVAR_STEPS);
    assert_eq!(c.affected_steps, CLINVAR_STEPS);
}

#[test]
fn a_newer_clinvar_installed_here_is_rerun_ready() {
    let g = genome(
        "cv-here",
        Opt {
            installed_clinvar: Some("20261005"),
            ..Opt::default()
        },
    );
    let c = run_check(&g, &answers(&[]), &this_pins(), Some("3.66.0"));
    assert_eq!(by_item(&c)["clinvar_file_date"].class, Class::RerunReady);
    assert_eq!(c.rerun_ready_steps, CLINVAR_STEPS);
}

#[test]
fn vep_117_in_the_newest_release_is_in_pipeline_release_else_upstream_only() {
    let g = genome("vep", Opt::default());
    let rel = r#"{"tag": "v9.9.0", "pins": {"VEP_CACHE_RELEASE": "117"}}"#;
    let c = run_check(
        &g,
        &answers(&[("release", rel)]),
        &this_pins(),
        Some("3.66.0"),
    );
    let v = by_item(&c)["VEP_CACHE_RELEASE"];
    assert_eq!(v.class, Class::InPipelineRelease);
    assert_eq!(v.in_release, Some(("v9.9.0".into(), "117".into())));
    assert!(v.action.contains("v9.9.0"));
    let rel = r#"{"tag": "v9.9.0", "pins": {"VEP_CACHE_RELEASE": "116"}}"#;
    let c = run_check(
        &g,
        &answers(&[("release", rel)]),
        &this_pins(),
        Some("3.66.0"),
    );
    assert_eq!(by_item(&c)["VEP_CACHE_RELEASE"].class, Class::UpstreamOnly);
}

#[test]
fn pgs_files_newer_upstream_or_here() {
    let g = genome("pgs", Opt::default());
    let c = run_check(&g, &answers(&[]), &this_pins(), Some("3.66.0"));
    let p = by_item(&c)["pgs:PGS000001"];
    assert_eq!(p.class, Class::UpstreamOnly);
    assert_eq!(p.steps, ["24", "25", "26"]);
    let g = genome(
        "pgs-here",
        Opt {
            installed_pgs: Some("2026-09-01"),
            ..Opt::default()
        },
    );
    let c = run_check(&g, &answers(&[]), &this_pins(), Some("3.66.0"));
    assert_eq!(by_item(&c)["pgs:PGS000001"].class, Class::RerunReady);
    assert_eq!(c.rerun_ready_steps, ["24", "25", "26"]);
}

#[test]
fn what_cannot_be_compared_is_an_error() {
    let g = genome(
        "impossible",
        Opt {
            run_clinvar: Some("20260230"),
            ..Opt::default()
        },
    );
    let c = run_check(&g, &answers(&[]), &this_pins(), Some("3.66.0"));
    let cv = by_item(&c)["clinvar_file_date"];
    assert_eq!(cv.class, Class::Error);
    assert!(cv.error.as_deref().unwrap().contains("which ClinVar file"));
    let g = genome(
        "nodate",
        Opt {
            run_clinvar: Some(""),
            ..Opt::default()
        },
    );
    let c = run_check(&g, &answers(&[]), &this_pins(), Some("3.66.0"));
    assert_eq!(by_item(&c)["clinvar_file_date"].class, Class::Error);
    let g = genome("latest", Opt::default());
    let c = run_check(
        &g,
        &answers(&[("PHARMCAT_IMAGE", r#"{"upstream": "latest"}"#)]),
        &this_pins(),
        Some("3.66.0"),
    );
    let ph = by_item(&c)["PHARMCAT_IMAGE"];
    assert_eq!(ph.class, Class::Error);
    assert!(ph.error.as_deref().unwrap().contains("names no release"));
    let c = run_check(&g, &answers(&[]), &this_pins(), None);
    let hl = by_item(&c)["hla_database"];
    assert_eq!(hl.class, Class::Error);
    assert!(hl.error.as_deref().unwrap().contains("HLA_DB_RELEASE"));
    let c = run_check(
        &g,
        &answers(&[("release", r#"{"pins": {}}"#)]),
        &this_pins(),
        Some("3.66.0"),
    );
    assert!(c.release.as_ref().is_err_and(|e| e.contains("tag")));
}

/// A server on 127.0.0.1 that answers every request with `status`, or the
/// PGS score answer and a HEAD with this `Last-Modified` when one is given.
fn serve(last_modified: Option<&'static str>) -> u16 {
    let l = TcpListener::bind("127.0.0.1:0").unwrap();
    let port = l.local_addr().unwrap().port();
    std::thread::spawn(move || {
        for mut st in l.incoming().flatten() {
            let mut r = BufReader::new(st.try_clone().unwrap());
            let mut line = String::new();
            r.read_line(&mut line).unwrap();
            loop {
                let mut h = String::new();
                if r.read_line(&mut h).unwrap_or(0) == 0 || h.trim().is_empty() {
                    break;
                }
            }
            let reply = match last_modified {
                None => "HTTP/1.1 500 Error\r\nContent-Length: 0\r\nConnection: close\r\n\r\n"
                    .to_owned(),
                Some(lm) if line.starts_with("HEAD") => {
                    format!(
                        "HTTP/1.1 200 OK\r\nLast-Modified: {lm}\r\nContent-Length: 0\r\nConnection: close\r\n\r\n"
                    )
                }
                Some(_) => {
                    let body = format!(
                        r#"{{"id": "PGS000001", "ftp_harmonized_scoring_files": {{"GRCh38": {{"positions": "http://127.0.0.1:{port}/f/PGS000001_hmPOS_GRCh38.txt.gz"}}}}}}"#
                    );
                    format!(
                        "HTTP/1.1 200 OK\r\nContent-Type: application/json\r\nContent-Length: {}\r\nConnection: close\r\n\r\n{body}",
                        body.len()
                    )
                }
            };
            let _ = st.write_all(reply.as_bytes());
            let mut rest = [0u8; 1];
            let _ = st.read(&mut rest);
        }
    });
    port
}

fn online(port: u16) -> Online {
    let mut o = answers::online(
        Pins::default(),
        "o/r",
        &format!("http://127.0.0.1:{port}/rest/score"),
        &format!("http://127.0.0.1:{port}"),
        &repo(),
    );
    o.client = og_fresh::http::Client::new(0);
    o.client.github = GitHub {
        api: format!("http://127.0.0.1:{port}"),
        token: None,
    };
    o
}

/// The online release and PGS lookups, the offline answers for the rest.
struct Mixed(Online, Offline);

impl Answers for Mixed {
    fn clinvar(&self) -> og_fresh::Lookup<Option<String>> {
        self.1.clinvar()
    }
    fn pinned(&self, key: &str) -> og_fresh::Lookup<Option<String>> {
        self.1.pinned(key)
    }
    fn pgs(&self, id: &str) -> og_fresh::Lookup<Option<String>> {
        self.0.pgs(id)
    }
    fn release(&self) -> og_fresh::Lookup<answers::Release> {
        self.0.release()
    }
    fn hosts(&self) -> Vec<String> {
        self.0.hosts()
    }
    fn offline(&self) -> bool {
        false
    }
}

#[test]
fn a_failed_lookup_is_an_error_and_the_command_exits_1() {
    let g = genome("fail", Opt::default());
    let c = run_check(
        &g,
        &answers(&[("PHARMCAT_IMAGE", r#"{"error": "HTTP 500"}"#)]),
        &this_pins(),
        Some("3.66.0"),
    );
    let bad = by_item(&c)["PHARMCAT_IMAGE"];
    assert_eq!(bad.class, Class::Error);
    assert!(bad.error.as_deref().unwrap().contains("500"));
    let port = serve(None);
    let m = Mixed(online(port), answers(&[]));
    let c = run_check(&g, &m, &this_pins(), Some("3.66.0"));
    let it = by_item(&c);
    assert_eq!(it["pgs:PGS000001"].class, Class::Error);
    assert!(
        it["pgs:PGS000001"]
            .error
            .as_deref()
            .unwrap()
            .contains("500")
    );
    assert_eq!(it["VEP_CACHE_RELEASE"].class, Class::Error);
    assert_eq!(c.hosts, [format!("127.0.0.1:{port}")]);
    let p = tmp("a500").join("a500.json");
    let text = fs::read_to_string(fix("offline.json")).unwrap().replace(
        r#""PHARMCAT_IMAGE": {"upstream": "3.4.0"}"#,
        r#""PHARMCAT_IMAGE": {"error": "HTTP 500"}"#,
    );
    fs::write(&p, text).unwrap();
    let (rc, _, _) = cli(&[
        "--sample",
        "S1",
        "--genome-dir",
        s(&g.root),
        "--offline",
        s(&p),
    ]);
    assert_eq!(rc, 1);
}

#[test]
fn a_malformed_last_modified_of_a_pgs_file_is_an_error_not_a_crash() {
    let g = genome("lastmod", Opt::default());
    let m = Mixed(online(serve(Some("not a date"))), answers(&[]));
    let c = run_check(&g, &m, &this_pins(), Some("3.66.0"));
    let p = by_item(&c)["pgs:PGS000001"];
    assert_eq!(p.class, Class::Error);
    assert!(
        p.error.as_deref().unwrap().contains("Last-Modified"),
        "{p:?}"
    );
    let ok = Mixed(
        online(serve(Some("Mon, 05 Oct 2026 10:00:00 GMT"))),
        answers(&[]),
    );
    assert_eq!(ok.pgs("PGS000001").unwrap().as_deref(), Some("20261005"));
}

#[test]
fn state_json_is_replaced_whole_and_log_tsv_only_grows() {
    let g = genome("files", Opt::default());
    let off = fix("offline.json");
    let args = [
        "--sample",
        "S1",
        "--genome-dir",
        s(&g.root),
        "--offline",
        s(&off),
    ];
    cli(&args);
    let d = g.sample().join("recheck");
    let first = fs::read_to_string(d.join("log.tsv")).unwrap();
    cli(&args);
    let log = fs::read_to_string(d.join("log.tsv")).unwrap();
    let mut names: Vec<String> = fs::read_dir(&d)
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    names.sort();
    assert!(log.starts_with(&first) && log.lines().count() > first.lines().count());
    assert!(log.starts_with("utc\tevent\t"));
    assert_eq!(names, ["log.tsv", "state.json"]);
    let state = fs::read_to_string(d.join("state.json")).unwrap();
    assert!(state.contains("\"sample\": \"S1\"") && state.ends_with("}\n"));
}

#[test]
fn rerun_refuses_while_a_run_is_going_on() {
    let g = genome(
        "running",
        Opt {
            running: true,
            ..Opt::default()
        },
    );
    let off = fix("offline.json");
    let (rc, _, err) = cli(&[
        "--sample",
        "S1",
        "--genome-dir",
        s(&g.root),
        "--offline",
        s(&off),
        "--rerun",
    ]);
    assert_eq!(rc, 3);
    assert!(err.contains("going on"), "{err}");
}

#[test]
fn rerun_refuses_while_another_rerun_holds_the_lock() {
    let g = genome("lock", Opt::default());
    fs::create_dir_all(g.sample().join("recheck")).unwrap();
    let held = fs::OpenOptions::new()
        .append(true)
        .create(true)
        .open(g.sample().join("recheck/rerun.lock"))
        .unwrap();
    held.try_lock().unwrap();
    let off = fix("offline.json");
    let (rc, _, err) = cli(&[
        "--sample",
        "S1",
        "--genome-dir",
        s(&g.root),
        "--offline",
        s(&off),
        "--rerun",
    ]);
    assert_eq!(rc, 3);
    assert!(err.contains("another --rerun"), "{err}");
}

#[test]
fn after_rerun_the_check_is_taken_again_on_the_manifest_the_rerun_wrote() {
    let base = tmp("after");
    let newman = base.join("new_manifest.tsv");
    let text: String = fs::read_to_string(fix("run_manifest.tsv"))
        .unwrap()
        .lines()
        .map(|l| {
            if l.starts_with("data\tclinvar_file_date\t") {
                "data\tclinvar_file_date\t20261005 (clinvar/clinvar_pathogenic_chr.vcf.gz)\n"
                    .to_owned()
            } else {
                format!("{l}\n")
            }
        })
        .collect();
    fs::write(&newman, text).unwrap();
    let g = genome("after-g", Opt::default());
    let cmd = format!("cp '{}' '{}'", s(&newman), s(&g.manifest()));
    let g = Genome {
        root: {
            let st = fs::read_to_string(g.status())
                .unwrap()
                .replace("echo replayed", &cmd);
            fs::write(g.status(), st).unwrap();
            gz(
                &g.root.join("clinvar/clinvar_pathogenic_chr.vcf.gz"),
                "##fileformat=VCFv4.1\n##fileDate=20261005\n#CHROM\n",
            );
            g.root
        },
    };
    let off = fix("offline.json");
    let (rc, out, _) = cli(&[
        "--sample",
        "S1",
        "--genome-dir",
        s(&g.root),
        "--offline",
        s(&off),
        "--rerun",
    ]);
    assert_eq!(rc, 0, "{out}");
    let st = fs::read_to_string(g.sample().join("recheck/state.json")).unwrap();
    let v = pycompat::json::loads(&st).unwrap();
    let items = match pycompat::value::get(&v, "items").unwrap().unwrap() {
        Value::Array(a) => a.clone(),
        _ => unreachable!(),
    };
    let cv = &items[0];
    assert_eq!(
        pycompat::value::get(cv, "class").unwrap(),
        Some(&Value::Str("no-newer-release".into()))
    );
    assert_eq!(
        pycompat::value::get(cv, "ran_with").unwrap(),
        Some(&Value::Str("20261005".into()))
    );
    let rerun = pycompat::value::get(&v, "rerun").unwrap().unwrap();
    assert_eq!(
        pycompat::value::get(rerun, "rc").unwrap(),
        Some(&Value::Int(0))
    );
    let log = fs::read_to_string(g.sample().join("recheck/log.tsv")).unwrap();
    assert!(log.contains("\trerun-start\t") && log.contains("\trerun-end\t"));
}

#[test]
fn refresh_clinvar_refuses_while_another_sample_runs_or_a_run_holds_the_lock() {
    let base = tmp("refresh");
    let ran = base.join("ran");
    let setup = base.join("setup.sh");
    fs::write(&setup, format!("#!/bin/sh\necho \"$@\" >> '{}'\n", s(&ran))).unwrap();
    std::process::Command::new("chmod")
        .arg("755")
        .arg(&setup)
        .status()
        .unwrap();
    let g = genome("refresh-g", Opt::default());
    fs::create_dir_all(g.root.join("S2/logs")).unwrap();
    fs::write(
        g.root.join("S2/logs/run_status.tsv"),
        "meta\tstarted_epoch\t1000\nmeta\tstarted_utc\t2026-09-02 10:00:00 UTC\n",
    )
    .unwrap();
    let off = fix("offline.json");
    let args = [
        "--sample",
        "S1",
        "--genome-dir",
        s(&g.root),
        "--offline",
        s(&off),
        "--rerun",
        "--refresh-clinvar",
        "--setup-script",
        s(&setup),
    ];
    let (rc, _, err) = cli(&args);
    assert_eq!(rc, 3);
    assert!(err.contains("S2"), "{err}");
    assert!(!ran.exists());
    fs::remove_dir_all(g.root.join("S2")).unwrap();
    let held = fs::OpenOptions::new()
        .append(true)
        .create(true)
        .open(g.root.join(CLINVAR_LOCK))
        .unwrap();
    // Blocking: a child of another test thread may hold a copy of this
    // file's lock from the call above until it execs.
    held.lock_shared().unwrap();
    let (rc, _, err) = cli(&args);
    assert_eq!(rc, 3);
    assert!(err.contains(".refresh.lock"), "{err}");
    assert!(!ran.exists());
}

/// A shell line that appends `held` or `free` to `out`: whether another
/// process can take a shared `flock(2)` on `lock`, as `flock(1)` in
/// `run-all.sh` takes it, asking `tries` times: a lock copied into a child
/// of another test thread between its fork and its exec must not read as
/// held.
fn probe(lock: &str, out: &Path, tries: u32) -> String {
    format!(
        "perl -e 'use Fcntl qw(:flock); for (1..{tries}) {{ open(F, \">>\", $ARGV[0]) or die; if (flock(F, LOCK_SH | LOCK_NB)) {{ print \"free\\n\"; exit }} close F; select(undef, undef, undef, 0.02) }} print \"held\\n\"' \"{lock}\" >> '{}'",
        s(out)
    )
}

#[test]
fn refresh_clinvar_holds_the_lock_while_setup_refreshes_and_lets_it_go_before_the_replay() {
    let base = tmp("refresh-held");
    let seen = base.join("seen");
    let setup = base.join("setup.sh");
    fs::write(
        &setup,
        format!(
            "#!/bin/sh\n{}\n",
            probe("$3/clinvar/.refresh.lock", &seen, 1)
        ),
    )
    .unwrap();
    std::process::Command::new("chmod")
        .arg("755")
        .arg(&setup)
        .status()
        .unwrap();
    let g = genome("refresh-held-g", Opt::default());
    let replay = probe(s(&g.root.join(CLINVAR_LOCK)), &seen, 20);
    let st = fs::read_to_string(g.status())
        .unwrap()
        .replace("echo replayed", &replay);
    fs::write(g.status(), st).unwrap();
    let off = fix("offline.json");
    let (rc, out, err) = cli(&[
        "--sample",
        "S1",
        "--genome-dir",
        s(&g.root),
        "--offline",
        s(&off),
        "--rerun",
        "--refresh-clinvar",
        "--setup-script",
        s(&setup),
    ]);
    assert_eq!(rc, 0, "{out}{err}");
    assert_eq!(fs::read_to_string(&seen).unwrap(), "held\nfree\n");
}

#[test]
fn the_rerun_comparison_names_each_difference() {
    let st = steps();
    let cmp = compare(
        &["06".into(), "13".into(), "14".into()],
        vec!["06".into(), "07".into()],
        &st,
    );
    assert_eq!(cmp.redone, ["06", "07"]);
    assert_eq!(cmp.predicted_not_redone, ["13"]);
    assert_eq!(cmp.redone_not_predicted, ["07"]);
    assert_eq!(cmp.not_redone_script_steps, ["14"]);
}

#[test]
fn launch_dir_gives_the_last_command_with_resume() {
    let d = tmp("launch");
    fs::create_dir_all(d.join(".nextflow")).unwrap();
    fs::write(
        d.join(".nextflow/history"),
        "2026-09-01 10:00:00\t1h\tfirst_run\tOK\tabc\tsess1\tnextflow run main.nf --input a.csv\n\
         2026-09-02 10:00:00\t2h\tlast_run\tOK\tabc\tsess2\tnextflow run main.nf --input b.csv\n",
    )
    .unwrap();
    assert_eq!(
        launch_dir_command(&d).as_deref(),
        Some("nextflow run main.nf --input b.csv -resume")
    );
    fs::write(
        d.join(".nextflow/history"),
        "2026-09-02 10:00:00\t2h\tlast_run\tOK\tabc\tsess2\tnextflow run main.nf --input 'data -resume old.csv'\n",
    )
    .unwrap();
    assert_eq!(
        launch_dir_command(&d).as_deref(),
        Some("nextflow run main.nf --input 'data -resume old.csv' -resume")
    );
    let (rc, _, err) = cli(&["--launch-dir", s(&tmp("empty"))]);
    assert_eq!(rc, 2);
    assert!(err.contains("no .nextflow/history"));
}

fn with(k: &str, v: &str) -> IndexMap<String, String> {
    let mut p = this_pins();
    p.insert(k.into(), v.into());
    p
}

#[test]
fn a_new_digest_at_the_same_tag_is_rerun_ready() {
    let reference = this_pins()["PHARMCAT_IMAGE"].clone();
    let (d1, d2) = (
        format!("{reference}@sha256:{}", "a".repeat(64)),
        format!("{reference}@sha256:{}", "b".repeat(64)),
    );
    let g = genome("digest", Opt::default());
    g.add_manifest(&[("MEGANE_IMAGE", &d1)]);
    let new = run_check(
        &g,
        &answers(&[]),
        &with("MEGANE_IMAGE", &d2),
        Some("3.66.0"),
    );
    let same = run_check(
        &g,
        &answers(&[]),
        &with("MEGANE_IMAGE", &d1),
        Some("3.66.0"),
    );
    let mg = by_item(&new)["MEGANE_IMAGE"];
    assert_eq!(mg.class, Class::RerunReady);
    assert!(mg.steps.iter().any(|x| x == "50"));
    assert_eq!(mg.refs.as_ref().unwrap().1.as_deref(), Some(d2.as_str()));
    assert_eq!(by_item(&same)["MEGANE_IMAGE"].class, Class::NoNewerRelease);

    let gp = genome("digest-pharmcat", Opt::default());
    let man: String = fs::read_to_string(gp.manifest())
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with("version\tPHARMCAT_IMAGE\t"))
        .map(|l| format!("{l}\n"))
        .collect();
    fs::write(
        gp.manifest(),
        man + &format!("version\tPHARMCAT_IMAGE\t{d1}\n"),
    )
    .unwrap();
    let ph = run_check(
        &gp,
        &answers(&[]),
        &with("PHARMCAT_IMAGE", &d2),
        Some("3.66.0"),
    );
    let ph1 = run_check(
        &gp,
        &answers(&[]),
        &with("PHARMCAT_IMAGE", &d1),
        Some("3.66.0"),
    );
    let item = by_item(&ph)["PHARMCAT_IMAGE"].clone();
    assert_eq!(item.class, Class::RerunReady);
    assert_eq!(by_item(&ph1)["PHARMCAT_IMAGE"].class, Class::NoNewerRelease);
    let mut one = ph.clone();
    one.items = vec![item];
    one.affected_steps = vec![];
    let mut text = Vec::new();
    print_check(&one, &mut text).unwrap();
    assert!(
        String::from_utf8(text)
            .unwrap()
            .contains("same tag, new image")
    );
}

#[test]
fn a_pin_with_no_lookup_the_manifest_lacks_is_unknown_and_exits_0() {
    let d1 = format!(
        "{}@sha256:{}",
        this_pins()["PHARMCAT_IMAGE"],
        "a".repeat(64)
    );
    let g = genome("unknown", Opt::default());
    let c = run_check(
        &g,
        &answers(&[]),
        &with("MEGANE_IMAGE", &d1),
        Some("3.66.0"),
    );
    let mu = by_item(&c)["MEGANE_IMAGE"];
    assert_eq!(mu.class, Class::Unknown);
    assert!(!c.affected_steps.iter().any(|x| x == "50"));
    assert!(!c.rerun_ready_steps.iter().any(|x| x == "50"));
    assert!(mu.action.contains("older version"));
    assert_eq!(c.errors(), 0);
}

#[test]
fn a_pin_read_only_by_a_step_the_run_did_not_run_is_left_out() {
    let d1 = format!(
        "{}@sha256:{}",
        this_pins()["PHARMCAT_IMAGE"],
        "a".repeat(64)
    );
    let d2 = format!(
        "{}@sha256:{}",
        this_pins()["PHARMCAT_IMAGE"],
        "b".repeat(64)
    );
    // A default run (no TOOLS): run-all.sh writes rows for the steps it ran
    // and none for step 50.
    let gd = genome("not-in-tools", Opt::default());
    gd.add_manifest(&[("MEGANE_IMAGE", &d1)]);
    gd.add_steps(&[("06", "planned"), ("06", "ok")]);
    let nd = run_check(
        &gd,
        &answers(&[]),
        &with("MEGANE_IMAGE", &d2),
        Some("3.66.0"),
    );
    assert!(!by_item(&nd).contains_key("MEGANE_IMAGE"));
    assert!(!nd.rerun_ready_steps.iter().any(|x| x == "50"));
    let gs = genome("skipped", Opt::default());
    gs.add_manifest(&[("MEGANE_IMAGE", &d1)]);
    gs.add_steps(&[("06", "ok"), ("50", "skipped (no BAM)")]);
    let sk = run_check(
        &gs,
        &answers(&[]),
        &with("MEGANE_IMAGE", &d2),
        Some("3.66.0"),
    );
    assert!(!by_item(&sk).contains_key("MEGANE_IMAGE"));
    let gr = genome("ran", Opt::default());
    gr.add_manifest(&[("MEGANE_IMAGE", &d1)]);
    gr.add_steps(&[("06", "ok"), ("50", "planned"), ("50", "ok")]);
    let rn = run_check(
        &gr,
        &answers(&[]),
        &with("MEGANE_IMAGE", &d2),
        Some("3.66.0"),
    );
    assert_eq!(by_item(&rn)["MEGANE_IMAGE"].class, Class::RerunReady);
}

/// Steps 50 to 53 read pins of versions.env that change what they find: a
/// bump of each must mark its step and the report, and nothing else.
const OPT_IN_PINS: &[(&str, &str)] = &[
    ("MEGANE_IMAGE", "50"),
    ("MOCHA_PANEL_VERSION", "51"),
    ("MOCHA_IMAGE", "51"),
    ("VNTYPER_IMAGE", "52"),
    ("LOCITYPER_PANEL_VERSION", "53"),
    ("LOCITYPER_IMAGE", "53"),
];

fn pins_ok(steps: &[Step]) -> bool {
    OPT_IN_PINS.iter().all(|(k, n)| {
        let mut got = affected_steps(steps, &[k]);
        got.sort();
        let mut want = vec![(*n).to_owned(), REPORT_STEP.to_owned()];
        want.sort();
        got == want
    })
}

#[test]
fn a_bumped_pin_of_steps_50_to_53_marks_its_step() {
    let st = steps();
    assert!(pins_ok(&st));
    let blank: Vec<Step> = st
        .iter()
        .cloned()
        .map(|mut x| {
            if x.step == "53" {
                x.reads.clear();
            }
            x
        })
        .collect();
    assert!(!pins_ok(&blank), "passed with row 53 without reads");
    let no_report: Vec<Step> = st
        .iter()
        .filter(|x| x.step != REPORT_STEP)
        .cloned()
        .collect();
    assert!(!pins_ok(&no_report), "passed without the report step");
}

#[test]
fn a_run_on_an_older_locityper_panel_names_step_53_rerun_ready() {
    let g = genome("53", Opt::default());
    g.add_manifest(&[("LOCITYPER_PANEL_VERSION", "panel-1.0")]);
    let c = run_check(
        &g,
        &answers(&[]),
        &with("LOCITYPER_PANEL_VERSION", "panel-2.0"),
        Some("3.66.0"),
    );
    let same = run_check(
        &g,
        &answers(&[]),
        &with("LOCITYPER_PANEL_VERSION", "panel-1.0"),
        Some("3.66.0"),
    );
    assert!(c.rerun_ready_steps.iter().any(|x| x == "53"));
    assert!(!same.rerun_ready_steps.iter().any(|x| x == "53"));
    assert!(
        by_item(&c)["LOCITYPER_PANEL_VERSION"]
            .action
            .contains("./scripts/setup.sh --locityper-data $GENOME_DIR")
    );
}

#[test]
fn usage_errors_exit_2() {
    assert_eq!(
        cli(&["--sample", "S1", "--genome-dir", "/nonexistent"]).0,
        2
    );
    assert_eq!(
        cli(&["--sample", "S1", "--genome-dir", "/x", "--refresh-clinvar"]).0,
        2
    );
    assert_eq!(cli(&["--bogus"]).0, 2);
    assert_eq!(
        cli(&["--affected", "clinvar_file_date"]).1,
        format!("{}\n", CLINVAR_STEPS.join(" "))
    );
    assert_eq!(PGS_API, "https://www.pgscatalog.org/rest/score");
}

#[test]
fn small_readers() {
    assert_eq!(date8("2026-10-05").as_deref(), Some("20261005"));
    assert_eq!(date8("20260230"), None);
    assert_eq!(date8("x"), None);
    assert_eq!(
        pgs_field("pgs_id=P;HmPOS_date=2022-07-29", "HmPOS_date").as_deref(),
        Some("2022-07-29")
    );
    assert!(key_matches("pgs:PGS1", "pgs:*@self") && !key_matches("pgsx", "pgs:*"));
    assert_eq!(
        pin_value("A_IMAGE", Some("r/x:1.2@sha256:a")).as_deref(),
        Some("1.2")
    );
    assert_eq!(
        pin_value("hla_database", Some("IPD-IMGT/HLA 3.66.0")).as_deref(),
        Some("3.66.0")
    );
}

#[test]
fn state_json_written_item_by_item_is_byte_for_byte_the_whole_document() {
    use crate::files::{write_check, write_json};
    let g = genome("stream", Opt::default());
    let mut c = run_check(&g, &answers(&[]), &this_pins(), Some("3.66.0"));
    let d = tmp("stream-out");
    let same = |c: &Check, name: &str| {
        let (a, b) = (
            d.join(format!("{name}-whole.json")),
            d.join(format!("{name}-items.json")),
        );
        write_json(&a, &c.to_json()).unwrap();
        write_check(&b, c).unwrap();
        let (a, b) = (
            fs::read_to_string(a).unwrap(),
            fs::read_to_string(b).unwrap(),
        );
        assert_eq!(a, b, "{name}");
        a
    };
    let text = same(&c, "items");
    assert!(
        text.contains("\n  {\n   \"item\": \"clinvar_file_date\""),
        "{text}"
    );
    c.rerun = Some(pycompat::json::loads(r#"{"rc": 0, "command": "x \"y\"\n"}"#).unwrap());
    same(&c, "rerun");
    c.items.clear();
    assert!(same(&c, "empty").contains("\"items\": [],"));
}

#[test]
fn a_lookup_that_fails_after_the_replay_makes_the_command_exit_1() {
    let base = tmp("after-fail");
    let newman = base.join("new_manifest.tsv");
    let mut text = fs::read_to_string(fix("run_manifest.tsv")).unwrap();
    text += "data\tpgs:PGS000099\tpgs_id=PGS000099;HmPOS_date=2022-07-29\n";
    fs::write(&newman, text).unwrap();
    let g = genome("after-fail-g", Opt::default());
    let cmd = format!("cp '{}' '{}'", s(&newman), s(&g.manifest()));
    let st = fs::read_to_string(g.status())
        .unwrap()
        .replace("echo replayed", &cmd);
    fs::write(g.status(), st).unwrap();
    let off = fix("offline.json");
    let (rc, out, _) = cli(&[
        "--sample",
        "S1",
        "--genome-dir",
        s(&g.root),
        "--offline",
        s(&off),
        "--rerun",
    ]);
    let state = fs::read_to_string(g.sample().join("recheck/state.json")).unwrap();
    assert!(state.contains("pgs:PGS000099") && state.contains("\"class\": \"error\""));
    assert_eq!(rc, 1, "{out}");
}
