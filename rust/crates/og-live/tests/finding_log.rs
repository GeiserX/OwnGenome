//! The access log names no finding. Opening a finding's page on the live
//! site writes one access-log row; its section column is `-`, never the
//! finding's id, because the id is built from the gene and the position
//! (docs/live-site.md, "Privacy and what is stored": the log never holds a
//! result).
//!
//! The route is `bin/owngenome_live_results.py` while the Python site
//! serves; the check runs that site on a made-up run with one ClinVar hit in
//! an invented gene. A control runs the same site with the route's old log
//! line (the finding id in the section column): the check must fail there.
//! A check that cannot fail is no check.

use pycompat::json::{Value, loads};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Duration;

const RUN: &str = "FIND-001";
/// The invented gene and position of the run's one ClinVar hit.
const GENE: &str = "GENEQ";
const POS: &str = "7654321";
/// The route's log line as this fix writes it, and as it was before.
const FIXED: &str = "app.log(\"finding\", run[\"name\"], \"-\", 200)";
const OLD: &str = "app.log(\"finding\", run[\"name\"], fid, 200)";

fn repo() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../..")
}

/// A new empty folder of this test run.
fn tmp() -> PathBuf {
    static N: AtomicUsize = AtomicUsize::new(0);
    let d = Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!(
        "finding-log-{}-{}",
        std::process::id(),
        N.fetch_add(1, Ordering::Relaxed)
    ));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn free_port() -> u16 {
    std::net::TcpListener::bind("127.0.0.1:0")
        .unwrap()
        .local_addr()
        .unwrap()
        .port()
}

/// One GET on the local site: status and body; `None` while it does not
/// answer.
fn get(port: u16, path: &str) -> Option<(u16, String)> {
    use std::io::{Read, Write};
    let mut c = std::net::TcpStream::connect(("127.0.0.1", port)).ok()?;
    c.set_read_timeout(Some(Duration::from_secs(10))).ok()?;
    write!(c, "GET {path} HTTP/1.0\r\nHost: 127.0.0.1\r\n\r\n").ok()?;
    let mut raw = Vec::new();
    c.read_to_end(&mut raw).ok()?;
    let raw = String::from_utf8(raw).ok()?;
    let (head, body) = raw.split_once("\r\n\r\n")?;
    let status = head.split(' ').nth(1)?.parse().ok()?;
    Some((status, body.to_owned()))
}

/// The site's process, stopped when the test ends, pass or fail.
struct Site(std::process::Child);

impl Drop for Site {
    fn drop(&mut self) {
        let _ = self.0.kill();
        let _ = self.0.wait();
    }
}

fn key<'a>(v: &'a Value, k: &str) -> Option<&'a Value> {
    match v {
        Value::Object(m) => m.get(k),
        _ => None,
    }
}

/// A made-up run folder under `root`: step 06 done (its task ended on
/// 2026-01-02, before `summary.json` was written) and a summary with one
/// ClinVar hit in [`GENE`] at [`POS`].
fn make_run(root: &Path) {
    let run = root.join(RUN);
    std::fs::create_dir_all(run.join("logs")).unwrap();
    std::fs::create_dir_all(root.join("pipeline_info")).unwrap();
    std::fs::write(
        run.join("logs/run_status.tsv"),
        "meta\tstarted_epoch\t1767348000\nstep\t06\tok\n",
    )
    .unwrap();
    std::fs::write(
        root.join("pipeline_info/trace_1.txt"),
        format!(
            "task_id\tname\tstatus\tsubmit\tduration\n\
             1\tCLINVAR_SCREEN ({RUN})\tCOMPLETED\t2026-01-02 10:00:00.000\t1m\n"
        ),
    )
    .unwrap();
    std::fs::write(
        run.join("summary.json"),
        format!(
            "{{\"sections\": {{\"clinvar\": {{\"state\": \"ok\", \"data\": {{\"hits\": [\
             {{\"gene\": \"{GENE}\", \"chrom\": \"chr2\", \"pos\": {POS}, \"ref\": \"C\", \"alt\": \"T\", \
             \"significance\": \"Pathogenic\", \"stars\": 2, \"genotype\": \"het\", \"af\": 0.0001}}]}}}}}}}}\n"
        ),
    )
    .unwrap();
}

/// Starts the Python site on a new made-up run (`plugin_dir` in place of
/// `bin/` for its plugins when given), opens the run's one finding page,
/// and returns the finding's id and the access log.
fn open_finding(plugin_dir: Option<&Path>) -> Result<(String, String), String> {
    let base = tmp();
    let root = base.join("root");
    let state_dir = base.join("state");
    make_run(&root);
    let (port, iport) = (free_port(), free_port());
    let exe = std::env::var("OG_PYTHON").unwrap_or_else(|_| "python3".into());
    let err_log = base.join("site.stderr");
    let mut cmd = Command::new(&exe);
    cmd.arg(repo().join("bin/owngenome_live.py"))
        .arg("--root")
        .arg(&root)
        .arg("--state-dir")
        .arg(&state_dir)
        .args([
            "--port",
            &port.to_string(),
            "--ingest-port",
            &iport.to_string(),
            "--quiet",
        ])
        .env_remove("OWNGENOME_LIVE_PASSWORD_HASH")
        .stdout(Stdio::null())
        .stderr(std::fs::File::create(&err_log).unwrap());
    if let Some(d) = plugin_dir {
        cmd.arg("--plugin-dir").arg(d);
    }
    let mut site = Site(
        cmd.spawn()
            .map_err(|e| format!("cannot run {exe}: {e}; set OG_PYTHON"))?,
    );
    let why = || std::fs::read_to_string(&err_log).unwrap_or_default();
    let results = format!("/runs/{RUN}/results.json");
    let mut body = None;
    for _ in 0..200 {
        if let Some(code) = site.0.try_wait().unwrap() {
            return Err(format!("the site exited ({code}): {}", why()));
        }
        if let Some((200, b)) = get(port, &results) {
            body = Some(b);
            break;
        }
        std::thread::sleep(Duration::from_millis(100));
    }
    let body = body.ok_or_else(|| format!("results.json: no answer; {}", why()))?;
    let res = loads(&body).map_err(|e| format!("results.json: {e:?}"))?;
    let fid = match key(&res, "findings") {
        Some(Value::Array(found)) => match found.first().and_then(|f| key(f, "id")) {
            Some(Value::Str(id)) => id.clone(),
            other => return Err(format!("the first finding has no id: {other:?}")),
        },
        other => return Err(format!("results.json has no findings: {other:?}")),
    };
    match get(port, &format!("/runs/{RUN}/finding/{fid}/")) {
        Some((200, _)) => {}
        other => return Err(format!("the finding page answered {other:?}")),
    }
    let log = std::fs::read_to_string(state_dir.join("access.log"))
        .map_err(|e| format!("access.log: {e}"))?;
    Ok((fid, log))
}

/// What is wrong with the access log after the finding page opened: the
/// page must write one `finding` row whose section is `-`, and no row may
/// hold the finding's id, its gene or its position.
fn log_problems(fid: &str, log: &str) -> Vec<String> {
    let mut out = Vec::new();
    if !fid.contains(&GENE.to_lowercase()) || !fid.contains(POS) {
        out.push(format!(
            "the finding id {fid} names no gene and position: the check would see nothing"
        ));
    }
    let rows: Vec<Vec<&str>> = log
        .lines()
        .skip(1)
        .map(|l| l.split('\t').collect())
        .collect();
    let opened: Vec<&Vec<&str>> = rows
        .iter()
        .filter(|r| r.len() == 7 && r[3] == "finding" && r[6] == "200")
        .collect();
    if opened.len() != 1 || opened[0][4] != RUN || opened[0][5] != "-" {
        out.push(format!(
            "want one finding row of {RUN} with section '-': {opened:?}"
        ));
    }
    let lower = log.to_lowercase();
    for (what, s) in [
        ("id", fid.to_owned()),
        ("gene", GENE.to_lowercase()),
        ("position", POS.into()),
    ] {
        if lower.contains(&s) {
            out.push(format!("the access log holds the finding's {what}"));
        }
    }
    out
}

/// While `bin/owngenome_live.py` lives, opening a finding writes no result
/// into the access log. The control runs the site with the route's old log
/// line, which wrote the id there: the check must fail on it.
#[test]
fn opening_a_finding_page_writes_no_finding_id_into_the_access_log() {
    let bin = repo().join("bin");
    let route = bin.join("owngenome_live_results.py");
    if !bin.join("owngenome_live.py").is_file() || !route.is_file() {
        eprintln!("bin/owngenome_live.py is gone: the Python site check retired with it");
        return;
    }
    let (fid, log) = open_finding(None).unwrap_or_else(|e| panic!("{e}"));
    let bad = log_problems(&fid, &log);
    assert!(bad.is_empty(), "{bad:?}\n{log}");

    let src = std::fs::read_to_string(&route).unwrap();
    assert!(
        src.matches(FIXED).count() == 1,
        "the control found no route log line to put back"
    );
    let dir = tmp().join("plugins");
    std::fs::create_dir_all(&dir).unwrap();
    std::fs::write(
        dir.join("owngenome_live_results.py"),
        src.replacen(FIXED, OLD, 1),
    )
    .unwrap();
    let (fid, log) = open_finding(Some(&dir)).unwrap_or_else(|e| panic!("control: {e}"));
    assert!(
        !log_problems(&fid, &log).is_empty(),
        "the check passed the route that logs the finding id:\n{log}"
    );
}
