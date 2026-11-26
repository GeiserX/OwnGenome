//! The genome map on `og live serve`, with a fake og that records its
//! arguments and stdin and answers as told (the mode file beside it).
//!
//! 1. Without the results the map is never built: 503, og not run (the
//!    results hold the veil a child's page keeps closed);
//! 2. with them, `/runs/<run>/map.json` is og's bytes unchanged (equal to a
//!    direct og call with the same arguments and stdin), as
//!    application/gzip with the site's headers;
//! 3. og gets `map build --root ROOT --run RUN --assets <repo>/assets
//!    --cache-dir STATE/map`, and its stdin is `{"results": the body
//!    /runs/<run>/results.json sends}`;
//! 4. the gene files are og's bytes of `og map genes --chrom` and
//!    `--names`;
//! 5. og exit 3 is 404 with og's reason, 127 or no og 503 with the setup
//!    line, any other exit 500 with none of og's bytes;
//! 6. a family's run and an unknown run are 404 and og is not run;
//! 7. one access-log row `map` per request, naming only the run or
//!    `genes`; there is no POST route.

#[path = "port_live_server_util.rs"]
mod util;

use og_live::plugins::ResultsFns;
use pycompat::json::{Value, loads};
use std::collections::BTreeSet;
use std::fs;
use std::io::Write;
use std::os::unix::fs::PermissionsExt;
use std::path::PathBuf;
use std::process::{Command, Stdio};
use util::*;

const FAKE_OG: &str = r#"#!/bin/sh
d=$(dirname "$0")
n=$(ls "$d/calls" | wc -l | tr -d ' ')
mkdir -p "$d/calls/$n"
for a in "$@"; do printf '%s\n' "$a"; done > "$d/calls/$n/args"
cat > "$d/calls/$n/stdin"
mode=$(cat "$d/mode" 2>/dev/null || echo ok)
case "$mode" in
nomap) echo "og map build: no map for this run: the planted reason" >&2; exit 3;;
notbuilt) echo "ERROR: og is not built for this host; run scripts/setup.sh --og" >&2; exit 127;;
crash) printf 'half a document'; printf 'og map build: the planted first line\nand a second line\n' >&2; exit 1;;
esac
cat "$d/payload.bin"
cat "$d/calls/$n/args" "$d/calls/$n/stdin" | cksum
printf '\r\n\000end'
"#;

const STEP06_OK: &str =
    r#"{"event": "pgp_step", "step": "06", "state": "ok", "utcTime": "2026-01-02T10:06:40Z"}"#;

const FINDING: &str = "clinvar-plantgene1-1";

/// The stub results: their answer carries every argument, so a site that
/// built them unlike results.json sends og another body.
fn build(name: &str, summary: &Value, files: Option<&BTreeSet<String>>) -> Result<Value, String> {
    let files: Vec<String> = files
        .map(|f| f.iter().cloned().collect())
        .unwrap_or_default();
    let text = format!(
        r#"{{"schema_version": 1, "sample": {}, "chapters": [], "note": null,
            "findings": [{{"id": "{FINDING}", "gene": "PLANTGENE1", "title": "planted finding", "level": "info"}}],
            "summary": {}, "files": {}}}"#,
        og_live::state::dumps(&Value::Str(name.into())),
        og_live::state::dumps(summary),
        og_live::state::dumps(&Value::Array(files.into_iter().map(Value::Str).collect())),
    );
    loads(&text).map_err(|e| e.to_string())
}

fn overview(_: &Value) -> Result<Value, String> {
    loads(r#"{"attention": 0, "levels": {}}"#).map_err(|e| e.to_string())
}

struct Map {
    w: Work,
    og: PathBuf,
}

impl Map {
    fn new(name: &str) -> Map {
        let w = Work::new(name);
        make_root(&w);
        let d = w.root().join(RUN);
        write(
            &d.join(format!("{RUN}_report.html")),
            "<a href=\"tables/linked.tsv\">a linked table</a>",
        );
        let ev = fixture_events();
        store(&w, RUN, &ev.iter().map(String::as_str).collect::<Vec<_>>());
        store(&w, RUN, &[STEP06_OK]);
        write(
            &d.join("summary.json"),
            r#"{"sections": {}, "planted": true}"#,
        );
        let fam = w.root().join("FAM-A");
        fs::create_dir_all(&fam).unwrap();
        fs::copy(
            repo().join("tests/fixtures/live/demo/DEMO-FAMILY/members.tsv"),
            fam.join("members.tsv"),
        )
        .unwrap();
        let bin = w.dir.join("bin");
        fs::create_dir_all(bin.join("calls")).unwrap();
        let og = bin.join("og");
        fs::write(&og, FAKE_OG).unwrap();
        fs::set_permissions(&og, fs::Permissions::from_mode(0o755)).unwrap();
        let mut payload = b"\x1f\x8b".to_vec();
        for _ in 0..1200 {
            payload.extend(0u8..=255);
        }
        fs::write(bin.join("payload.bin"), payload).unwrap();
        Map { w, og }
    }

    fn config(&self, results: bool) -> og_live::server::Config {
        let mut c = config(&self.w);
        c.og = self.og.clone();
        c.results = results.then_some(ResultsFns { build, overview });
        c
    }

    fn mode(&self, m: &str) {
        fs::write(self.og.with_file_name("mode"), m).unwrap();
    }

    /// The fake og's calls since the last read: (arguments, stdin).
    fn calls(&self) -> Vec<(Vec<String>, String)> {
        let d = self.og.with_file_name("calls");
        let mut names: Vec<usize> = fs::read_dir(&d)
            .unwrap()
            .filter_map(|e| e.ok()?.file_name().to_str()?.parse().ok())
            .collect();
        names.sort();
        let out = names
            .iter()
            .map(|n| {
                let c = d.join(n.to_string());
                let args = fs::read_to_string(c.join("args"))
                    .unwrap()
                    .lines()
                    .map(str::to_owned)
                    .collect();
                (args, fs::read_to_string(c.join("stdin")).unwrap())
            })
            .collect();
        fs::remove_dir_all(&d).unwrap();
        fs::create_dir_all(&d).unwrap();
        out
    }

    /// The fake og run directly with `args` and `stdin`: its stdout.
    fn direct(&self, args: &[String], stdin: &str) -> Vec<u8> {
        let mut p = Command::new(&self.og)
            .args(args)
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .spawn()
            .unwrap();
        p.stdin.take().unwrap().write_all(stdin.as_bytes()).unwrap();
        let out = p.wait_with_output().unwrap().stdout;
        self.calls();
        out
    }
}

fn sha(b: &[u8]) -> String {
    // A content digest without a hashing crate: equal bytes, equal length.
    format!(
        "{}:{:x}",
        b.len(),
        b.iter()
            .fold(0xcbf29ce484222325u64, |h, x| (h ^ u64::from(*x))
                .wrapping_mul(0x100000001b3))
    )
}

#[test]
fn without_the_results_the_map_is_never_built() {
    let m = Map::new("map-veil");
    let srv = start(m.config(false));
    let r = get(srv.port(), &format!("/runs/{RUN}/map.json"));
    assert_eq!(r.status, 503, "{}", r.text());
    assert_eq!(r.text(), og_live::map::NO_RESULTS);
    assert!(m.calls().is_empty(), "og ran without the results' veil");
    // The gene files hold no result: they are og's as ever.
    assert_eq!(get(srv.port(), "/map/genes/chr7.json").status, 200);
    assert_eq!(m.calls().len(), 1);
}

#[test]
fn the_map_is_ogs_bytes_with_the_sites_headers() {
    let m = Map::new("map");
    let srv = start(m.config(true));
    let p = srv.port();
    let root = m.w.root();
    let state = m.w.state();
    let n_rows = log_rows(&m.w).len();

    // 2, 3. the person's document
    let r = get(p, &format!("/runs/{RUN}/map.json"));
    let calls = m.calls();
    let results = get(p, &format!("/runs/{RUN}/results.json"));
    assert_eq!(r.status, 200, "{}", r.text());
    assert_eq!(calls.len(), 1);
    let direct = m.direct(&calls[0].0, &calls[0].1);
    assert!(r.body.len() > 300_000);
    assert_eq!(
        sha(&r.body),
        sha(&direct),
        "the map is not og's bytes unchanged"
    );
    assert_eq!(r.header("content-type"), Some("application/gzip"));
    assert_eq!(r.header("content-security-policy"), Some(og_live::csp::CSP));
    assert!(policy_ok(&r));
    let want: Vec<String> = [
        "map",
        "build",
        "--root",
        &root.to_string_lossy(),
        "--run",
        RUN,
        "--assets",
        &repo().join("assets").to_string_lossy(),
        "--cache-dir",
        &state.join("map").to_string_lossy(),
    ]
    .iter()
    .map(|x| (*x).to_owned())
    .collect();
    assert_eq!(calls[0].0, want);
    let body = results.json();
    let mut o = indexmap::IndexMap::new();
    o.insert("results".to_owned(), body.clone());
    assert_eq!(
        loads(&calls[0].1).unwrap(),
        Value::Object(o),
        "og's stdin is not results.json's body"
    );
    let text = og_live::state::dumps(&body);
    assert!(
        text.contains(FINDING)
            && text.contains("\"planted\": true")
            && text.contains("tables/linked.tsv"),
        "{text}"
    );

    // 4. the gene files
    for (path, extra) in [
        ("/map/genes/chr7.json", vec!["--chrom", "chr7"]),
        ("/map/genes/names.json", vec!["--names"]),
    ] {
        let r = get(p, path);
        let calls = m.calls();
        assert_eq!(
            (r.status, r.header("content-type")),
            (200, Some("application/gzip")),
            "{path}"
        );
        assert!(policy_ok(&r));
        let mut want: Vec<String> = [
            "map",
            "genes",
            "--root",
            &root.to_string_lossy(),
            "--cache-dir",
            &state.join("map").to_string_lossy(),
        ]
        .iter()
        .map(|x| (*x).to_owned())
        .collect();
        want.extend(extra.iter().map(|x| (*x).to_owned()));
        assert_eq!(calls.len(), 1);
        assert_eq!(calls[0].0, want, "{path}");
        assert_eq!(sha(&r.body), sha(&m.direct(&calls[0].0, "")), "{path}");
    }
    for bad in [
        "/map/genes/chr123.json",
        "/map/genes/x.json",
        "/map/genes/chr1.txt",
    ] {
        assert_eq!(get(p, bad).status, 404, "{bad}");
    }

    // 5. og's exit codes
    m.mode("nomap");
    let r = get(p, &format!("/runs/{RUN}/map.json"));
    assert_eq!(
        (r.status, r.text().as_str()),
        (
            404,
            "og map build: no map for this run: the planted reason\n"
        )
    );
    m.mode("notbuilt");
    let r = get(p, "/map/genes/chr7.json");
    assert_eq!((r.status, r.text().as_str()), (503, og_live::map::NO_OG));
    fs::rename(&m.og, m.og.with_extension("away")).unwrap();
    let r = get(p, "/map/genes/names.json");
    assert_eq!((r.status, r.text().as_str()), (503, og_live::map::NO_OG));
    fs::rename(m.og.with_extension("away"), &m.og).unwrap();
    m.mode("crash");
    let r = get(p, &format!("/runs/{RUN}/map.json"));
    assert_eq!((r.status, r.text().as_str()), (500, og_live::map::FAILED));
    assert!(!r.text().contains("half a document"));
    m.mode("ok");
    m.calls();

    // 6. runs with no person's map
    for run in ["FAM-A", "NO-SUCH-RUN"] {
        assert_eq!(
            get(p, &format!("/runs/{run}/map.json")).status,
            404,
            "{run}"
        );
    }
    assert!(
        m.calls().is_empty(),
        "og ran for a run with no person's map"
    );

    // 7. the access log, and no POST route
    let post = request(
        p,
        "POST",
        &format!("/runs/{RUN}/map.json"),
        &[("Content-Length", "0")],
    );
    assert_eq!(post.status, 404);
    let rows: Vec<Vec<String>> = log_rows(&m.w)[n_rows..]
        .iter()
        .filter(|r| r[3] == "map")
        .cloned()
        .collect();
    assert_eq!(rows.len(), 9, "{rows:?}");
    for r in &rows {
        let pair = (r[4].as_str(), r[5].as_str());
        assert!(
            [
                (RUN, "-"),
                ("-", "genes"),
                ("FAM-A", "-"),
                ("NO-SUCH-RUN", "-")
            ]
            .contains(&pair),
            "{r:?}"
        );
    }
    let log = fs::read_to_string(state.join("access.log")).unwrap();
    for w in ["PLANTGENE1", FINDING, "chr7", "planted"] {
        assert!(!log.contains(w), "{w} in the access log");
    }
}

#[test]
fn a_finding_page_is_served_only_for_a_finding_of_the_run() {
    let m = Map::new("finding");
    let srv = start(m.config(true));
    let p = srv.port();
    let r = get(p, &format!("/runs/{RUN}/finding/{FINDING}/"));
    assert_eq!(r.status, 200);
    assert!(r.text().contains(FINDING) || r.text().contains(RUN));
    for bad in [
        format!("/runs/{RUN}/finding/clinvar-other-1/"),
        format!("/runs/{RUN}/finding/NOT-WELL-FORMED/"),
        format!("/runs/FAM-A/finding/{FINDING}/"),
        format!("/runs/NO-SUCH/finding/{FINDING}/"),
    ] {
        assert_eq!(get(p, &bad).status, 404, "{bad}");
    }
    let rows = log_rows(&m.w);
    assert!(
        rows.iter()
            .any(|r| r[3] == "finding" && r[4] == RUN && r[6] == "200")
    );
    assert!(
        !fs::read_to_string(m.w.state().join("access.log"))
            .unwrap()
            .contains(FINDING)
    );
    // One row for the page served, and its section column says nothing
    // about the result: not the id, not the gene.
    let opened: Vec<_> = rows
        .iter()
        .filter(|r| r[3] == "finding" && r[6] == "200")
        .collect();
    assert_eq!(opened.len(), 1, "{rows:?}");
    assert_eq!(
        opened[0][5], "-",
        "the finding row's section names a result"
    );
    let log = fs::read_to_string(m.w.state().join("access.log"))
        .unwrap()
        .to_lowercase();
    assert!(!log.contains("plantgene1"), "the gene in the access log");
}
