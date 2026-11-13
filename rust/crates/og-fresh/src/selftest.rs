//! `og dev freshness --self-test`: each rule against a planted case, on
//! stubbed answers, a fake GitHub on 127.0.0.1 and, last, the real
//! registries with a bogus image, a missing tag and a bogus ClinVar URL.

use crate::http::{Body, Client, Fetch, GitHub, Headers, Request, Scratch};
use crate::issue::{self, Action, LABEL, MARKER};
use crate::lookups::{self, Ctx, image};
use crate::pins::Pins;
use crate::registry::{Online, Registry, split_image};
use crate::report::Report;
use crate::version::{change_kind, newer_tags};
use crate::{Lookup, LookupError, failed};
use indexmap::IndexMap;
use pycompat::json::{Dump, Value, dumps, loads};
use std::io::{BufRead, BufReader, Read, Write};
use std::net::{TcpListener, TcpStream};
use std::path::Path;
use std::sync::{Arc, Mutex};

/// A stubbed answer: the URL part it answers, the body (`None`: HTTP
/// 404) and its headers.
pub type StubAnswer = (&'static str, Option<Body>, Vec<(String, String)>);

/// A fetcher that answers from a table: the first key found in the URL
/// gives the body (`None`: HTTP 404), and every GitHub path answers
/// `github`.
pub struct Stub {
    pub answers: Vec<StubAnswer>,
    pub github: GitHub,
    pub sent: Mutex<Vec<String>>,
}

impl Stub {
    pub fn new(answers: Vec<(&'static str, Option<Body>)>) -> Stub {
        Stub {
            answers: answers.into_iter().map(|(k, b)| (k, b, vec![])).collect(),
            github: GitHub {
                api: "https://api.github.invalid".into(),
                token: None,
            },
            sent: Mutex::new(vec![]),
        }
    }
}

impl Fetch for Stub {
    fn github_api(&self) -> &GitHub {
        &self.github
    }

    fn request(&self, req: &Request<'_>) -> Lookup<(Body, Headers)> {
        self.sent.lock().expect("sent").push(req.url.to_owned());
        for (key, body, headers) in &self.answers {
            if req.url.contains(key) {
                return match body {
                    None => failed(format!("{}: HTTP 404", req.url)),
                    Some(b) => Ok((b.clone(), Headers::from_pairs(headers.clone()))),
                };
            }
        }
        failed(format!("{}: not stubbed", req.url))
    }
}

/// A registry that answers fixed tags and digests.
pub struct FakeRegistry {
    pub tags: fn(&str) -> Vec<String>,
    pub digests: fn(&str) -> String,
    pub asked: Mutex<Vec<String>>,
}

impl Registry for FakeRegistry {
    fn list_tags(&self, _host: &str, repo: &str) -> Lookup<Vec<String>> {
        Ok((self.tags)(repo))
    }
    fn manifest_digest(&self, _host: &str, _repo: &str, reference: &str) -> Lookup<String> {
        self.asked.lock().expect("asked").push(reference.into());
        Ok((self.digests)(reference))
    }
}

fn strs(xs: &[&str]) -> Vec<String> {
    xs.iter().map(|x| (*x).to_owned()).collect()
}

fn obj(pairs: Vec<(&str, Value)>) -> Value {
    Value::Object(pairs.into_iter().map(|(k, v)| (k.to_owned(), v)).collect())
}

fn s(x: &str) -> Value {
    Value::Str(x.into())
}

/// The issues API in memory: just enough for `issue::update`.
#[derive(Default)]
pub struct FakeGitHub {
    pub issues: Vec<IndexMap<String, Value>>,
    pub labels: Vec<String>,
    /// Creates the issue, then answers 502.
    pub lose_create_answer: bool,
}

fn reply(stream: &mut TcpStream, code: u16, body: &[u8], extra: &str) {
    let reason = match code {
        200 => "OK",
        201 => "Created",
        404 => "Not Found",
        _ => "Error",
    };
    let _ = write!(
        stream,
        "HTTP/1.1 {code} {reason}\r\nContent-Type: application/json\r\n{extra}Connection: close\r\n\r\n"
    );
    let _ = stream.write_all(body);
}

fn json_reply(stream: &mut TcpStream, code: u16, v: &Value) {
    let raw = dumps(v, &Dump::default()).expect("no floats");
    reply(
        stream,
        code,
        raw.as_bytes(),
        &format!("Content-Length: {}\r\n", raw.len()),
    );
}

fn handle(mut stream: TcpStream, gh: &Mutex<FakeGitHub>) {
    let mut reader = BufReader::new(stream.try_clone().expect("clone"));
    let mut line = String::new();
    if reader.read_line(&mut line).is_err() {
        return;
    }
    let mut parts = line.split_whitespace();
    let (method, target) = (
        parts.next().unwrap_or("").to_owned(),
        parts.next().unwrap_or("").to_owned(),
    );
    let mut len = 0usize;
    loop {
        let mut h = String::new();
        if reader.read_line(&mut h).is_err() || h.trim().is_empty() {
            break;
        }
        if let Some((n, v)) = h.split_once(':')
            && n.eq_ignore_ascii_case("content-length")
        {
            len = v.trim().parse().unwrap_or(0);
        }
    }
    let mut raw = vec![0; len];
    let _ = reader.read_exact(&mut raw);
    let body = loads(&String::from_utf8_lossy(&raw)).unwrap_or(Value::Object(IndexMap::new()));
    let get = |k: &str| match &body {
        Value::Object(m) => m.get(k).cloned().unwrap_or(Value::Null),
        _ => Value::Null,
    };
    let path = target.split('?').next().unwrap_or("").to_owned();
    let mut g = gh.lock().expect("fake github");
    match (method.as_str(), path.as_str()) {
        ("GET", "/empty") => reply(&mut stream, 200, b"", "Content-Length: 0\r\n"),
        // Promises 100 bytes and sends 10.
        ("GET", "/truncated") => reply(
            &mut stream,
            200,
            b"{\"a\": 1, \"",
            "Content-Length: 100\r\n",
        ),
        ("GET", p) if p.starts_with("/repos/o/r/labels/") => {
            let name = p.rsplit('/').next().unwrap_or("");
            if g.labels.iter().any(|l| l == name) {
                json_reply(&mut stream, 200, &obj(vec![("name", s(name))]));
            } else {
                json_reply(&mut stream, 404, &obj(vec![]));
            }
        }
        ("GET", "/repos/o/r/issues") => {
            let list = g
                .issues
                .iter()
                .filter(
                    |i| matches!(i.get("labels"), Some(Value::Array(l)) if l.contains(&s(LABEL))),
                )
                .map(|i| Value::Object(i.clone()))
                .collect();
            json_reply(&mut stream, 200, &Value::Array(list));
        }
        ("POST", "/repos/o/r/labels") => {
            if let Value::Str(n) = get("name") {
                g.labels.push(n);
            }
            json_reply(&mut stream, 201, &body);
        }
        ("POST", "/repos/o/r/issues") => {
            let n = g
                .issues
                .iter()
                .filter_map(|i| match i.get("number") {
                    Some(Value::Int(n)) => Some(*n),
                    _ => None,
                })
                .max()
                .unwrap_or(0)
                + 1;
            let labels = match get("labels") {
                Value::Null => Value::Array(vec![]),
                l => l,
            };
            let issue: IndexMap<String, Value> = [
                ("number", Value::Int(n)),
                ("state", s("open")),
                ("title", get("title")),
                ("body", get("body")),
                ("labels", labels),
                ("html_url", s(&format!("https://x/{n}"))),
            ]
            .into_iter()
            .map(|(k, v)| (k.to_owned(), v))
            .collect();
            g.issues.push(issue.clone());
            if g.lose_create_answer {
                json_reply(&mut stream, 502, &obj(vec![]));
            } else {
                json_reply(&mut stream, 201, &Value::Object(issue));
            }
        }
        ("PATCH", p) => {
            let n: i128 = p
                .rsplit('/')
                .next()
                .and_then(|x| x.parse().ok())
                .unwrap_or(-1);
            let found = g
                .issues
                .iter_mut()
                .find(|i| i.get("number") == Some(&Value::Int(n)));
            match (found, &body) {
                (Some(i), Value::Object(m)) => {
                    for (k, v) in m {
                        i.insert(k.clone(), v.clone());
                    }
                    let v = Value::Object(i.clone());
                    json_reply(&mut stream, 200, &v);
                }
                _ => json_reply(&mut stream, 404, &obj(vec![])),
            }
        }
        _ => json_reply(&mut stream, 404, &obj(vec![])),
    }
}

/// Serves `gh` on 127.0.0.1 until the process ends; returns its base URL.
pub fn serve(gh: Arc<Mutex<FakeGitHub>>) -> std::io::Result<String> {
    let l = TcpListener::bind("127.0.0.1:0")?;
    let port = l.local_addr()?.port();
    std::thread::spawn(move || {
        for s in l.incoming().flatten() {
            let gh = gh.clone();
            std::thread::spawn(move || handle(s, &gh));
        }
    });
    Ok(format!("http://127.0.0.1:{port}"))
}

/// The data-release rows (IPD-IMGT/HLA, IPD-KIR, pgsc_calc, dbSNP and
/// Parascopy) on stubbed answers: their states and their errors.
pub fn data_rows(
    answers: Vec<(&'static str, Option<Body>)>,
    release: Option<&str>,
) -> (Vec<String>, usize) {
    let mut answers = answers;
    let gh = match release {
        Some(r) => Body::Json(obj(vec![("tag_name", s(r))])),
        None => Body::Json(obj(vec![])),
    };
    answers.push(("api.github.invalid", Some(gh)));
    let stub = Stub::new(answers);
    let scratch = Scratch::new().expect("a private scratch folder");
    let d = scratch.path();
    let (common, setup) = (d.join("common.sh"), d.join("setup.sh"));
    let _ = std::fs::write(&common, "HLA_DB_RELEASE=${HLA_DB_RELEASE:-3.65.0}\n");
    let _ = std::fs::write(&setup, "PARASCOPY_DATA_RECORD=15019940\n");
    let pins = Pins::parse(
        "KIR_DB_RELEASE=\"2.15.0\"\nPGSC_CALC_VERSION=\"v2.3.0\"\nCONSUMER_DBSNP_BUILD=\"b157\"\nPARASCOPY_DATA_VERSION=\"1.7\"\n",
    );
    let reg = FakeRegistry {
        tags: |_| vec![],
        digests: |_| String::new(),
        asked: Mutex::new(vec![]),
    };
    let ctx = Ctx {
        fetch: &stub,
        registry: &reg,
        pins: &pins,
        versions: Path::new("/nonexistent"),
        setup: &setup,
        common: &common,
        now: jiff::Timestamp::now(),
    };
    let mut r = Report::default();
    for (what, f) in [
        ("IPD-IMGT/HLA", lookups::imgt as lookups::LookupFn),
        ("IPD-KIR", lookups::kir),
        ("pgsc_calc", lookups::pgsc_calc),
        ("dbSNP", lookups::dbsnp),
        ("Parascopy data", lookups::parascopy_data),
    ] {
        match f(&ctx) {
            Ok(found) => r.db_rows.push(found.row),
            Err(e) => r.db_error(what, &e),
        }
    }
    drop(scratch);
    (
        r.db_rows.into_iter().map(|x| x.state).collect(),
        r.errors.len(),
    )
}

fn text(t: &str) -> Option<Body> {
    Some(Body::Text(t.into()))
}

fn zenodo(name: &str) -> Option<Body> {
    // Another assembly's archive beside it must be ignored.
    let files = Value::Array(vec![
        obj(vec![("key", s(name))]),
        obj(vec![("key", s("GRCh37_v9.9.tar.gz"))]),
    ]);
    Some(Body::Json(obj(vec![(
        "hits",
        obj(vec![(
            "hits",
            Value::Array(vec![obj(vec![("id", Value::Int(1)), ("files", files)])]),
        )]),
    )])))
}

/// The checks; each prints `ok` or `FAIL`. Returns the exit code.
pub fn run(out: &mut dyn Write) -> u8 {
    let mut fails = 0;
    let mut check = |ok: bool, what: &str| {
        let _ = writeln!(out, "{}{what}", if ok { "ok   " } else { "FAIL " });
        if !ok {
            fails += 1;
        }
        fails
    };

    // 1. Version order: the sort -V traps.
    let bio = strs(&[
        "1.3.6--h43da1c4_0",
        "1.3.10--h5ca1c30_0",
        "1.3.9--hdeadbe_2",
        "1.3.6--habc123_1",
        "latest",
    ]);
    check(
        newer_tags("1.3.6--h43da1c4_0", &bio).ok()
            == Some((
                Some("1.3.10--h5ca1c30_0".into()),
                Some("1.3.10--h5ca1c30_0".into()),
            )),
        "biocontainers: 1.3.10 ranks above 1.3.9, whatever the build hash",
    );
    check(
        change_kind("1.3.6--h43da1c4_0", "1.3.6--habc123_1") == "build",
        "biocontainers: a new build is a build",
    );
    let reg = FakeRegistry {
        tags: |repo| {
            if repo == "library/python" {
                strs(&[
                    "3.11",
                    "3.9",
                    "3.14",
                    "3.14-slim",
                    "3.15.0rc1",
                    "3.11.9",
                    "3.15-rc",
                    "windowsservercore",
                    "3.14.0",
                ])
            } else {
                strs(&["1.20", "1.24"])
            }
        },
        digests: |_| String::new(),
        asked: Mutex::new(vec![]),
    };
    let mut rl = Report::default();
    for (var, r) in [
        ("PYTHON_IMAGE", "python:3.11"),
        ("SAMTOOLS_IMAGE", "staphb/samtools:1.20"),
    ] {
        if let Some(a) = image(&reg, r) {
            rl.add_image(var, "", a);
        }
    }
    check(
        rl.line_rows == vec![("PYTHON_IMAGE".into(), "3.11".into(), "3.14".into())]
            && rl
                .image_rows
                .iter()
                .map(|r| r.0.as_str())
                .collect::<Vec<_>>()
                == ["SAMTOOLS_IMAGE"],
        "python:3.11 is a release-line pin, not behind; samtools 1.20 -> 1.24 is behind",
    );
    let vep = strs(&[
        "release_116.0",
        "release_116.2",
        "release_117.0",
        "latest",
        "116.2",
    ]);
    check(
        newer_tags("release_116.0", &vep).ok()
            == Some((Some("release_117.0".into()), Some("release_116.2".into()))),
        "VEP: major and same-major found",
    );
    check(
        change_kind("release_116.0", "release_117.0") == "major",
        "VEP 116 -> 117 is a major",
    );
    check(
        newer_tags("1.21", &strs(&["1.21", "1.20", "1.21-2"])).ok() == Some((None, None)),
        "current pin has nothing newer",
    );
    check(
        newer_tags("latest", &strs(&["latest"])).is_err(),
        "a pin with no version is an error",
    );
    check(
        split_image("a/b@sha256:00").tag.is_none()
            && split_image("a/b@sha256:00").digest.as_deref() == Some("sha256:00"),
        "digest-only pin has no tag",
    );
    let reg = FakeRegistry {
        tags: |_| strs(&["3.11.17", "3.11.18", "3.14.0"]),
        digests: |_| String::new(),
        asked: Mutex::new(vec![]),
    };
    let mut rt = Report::default();
    if let Some(a) = image(&reg, "python:3.11.17@sha256:ab") {
        rt.add_image("PYTHON_IMAGE", "", a);
    }
    check(
        rt.image_rows
            .iter()
            .map(|r| (r.0.as_str(), r.1.as_str(), r.2.as_str()))
            .collect::<Vec<_>>()
            == [("PYTHON_IMAGE", "3.11.17", "3.14.0")],
        "a tag pinned to a digest is still checked for newer tags",
    );
    let reg = FakeRegistry {
        tags: |_| strs(&["latest"]),
        digests: |r| match r {
            "latest" => format!("sha256:{}", "1".repeat(64)),
            "3.11.17" => format!("sha256:{}", "2".repeat(64)),
            _ => format!("sha256:{}", "3".repeat(64)),
        },
        asked: Mutex::new(vec![]),
    };
    let tagged = lookups::digest(&reg, &format!("python:3.11.17@sha256:{}", "2".repeat(64)));
    let only = lookups::digest(&reg, &format!("a/b@sha256:{}", "1".repeat(64)));
    let rebuilt = lookups::digest(&reg, &format!("python:3.11.17@sha256:{}", "4".repeat(64)));
    let state = |x: &Option<Lookup<(String, String, String)>>| match x {
        Some(Ok((_, _, s))) => s.clone(),
        _ => String::new(),
    };
    check(
        reg.asked
            .lock()
            .expect("asked")
            .iter()
            .any(|r| r == "3.11.17")
            && state(&tagged).contains("still points at the pinned digest"),
        "a tag pinned to a digest is compared with its own tag, not with latest",
    );
    check(
        state(&only) == "same as publisher's latest",
        "a digest-only pin is compared with latest",
    );
    check(
        state(&rebuilt).contains("was rebuilt"),
        "a rebuilt tag is reported",
    );
    let mut cross = Stub::new(vec![]);
    cross.answers.push((
        "quay.io",
        Some(Body::Json(obj(vec![(
            "tags",
            Value::Array(vec![s("1.0")]),
        )]))),
        vec![(
            "Link".into(),
            "<https://elsewhere.example/v2/x/tags/list?last=1.0>; rel=\"next\"".into(),
        )],
    ));
    let refused = Online::list_tags_with(&cross, "quay.io", "biocontainers/x", None).is_err();
    check(
        refused
            && !cross
                .sent
                .lock()
                .expect("sent")
                .iter()
                .any(|u| u.contains("elsewhere.example")),
        "a tag-list next page on another host is refused, and nothing is sent there",
    );
    let shape = Stub::new(vec![(
        "rest.ensembl.org",
        Some(Body::Json(Value::Array(vec![]))),
    )]);
    let pins = Pins::parse(
        "VEP_CACHE_RELEASE=\"116\"\nVEP_IMAGE=\"ensemblorg/ensembl-vep:release_116.0\"\n",
    );
    let ctx = Ctx {
        fetch: &shape,
        registry: &reg,
        pins: &pins,
        versions: Path::new("/nonexistent"),
        setup: Path::new("/nonexistent"),
        common: Path::new("/nonexistent"),
        now: jiff::Timestamp::now(),
    };
    let mut rs = Report::default();
    if let Err(e) = lookups::ensembl(&ctx) {
        rs.db_error("odd shape", &e);
    }
    check(
        rs.errors
            .first()
            .is_some_and(|e| e.0 == "odd shape" && e.1.contains("unexpected answer"))
            && rs.db_rows.first().is_some_and(|r| r.state == "ERROR"),
        "an answer in an unexpected shape is an error row, not a crash",
    );

    // The rows for data releases, on stubbed answers: a newer upstream is
    // behind, the same one is current, and an answer without the version is
    // an error row, never "current".
    let newer = vec![
        (
            "IMGTHLA/Latest",
            text("# date: 2026-10-14\n# version: IPD-IMGT/HLA 3.66.0\n"),
        ),
        ("IPDKIR/Latest", text("# version: IPD-KIR 2.16.0\n")),
        (
            "snp/archive",
            text("<a href=\"b156/\">b156/</a> <a href=\"b158/\">b158/</a>"),
        ),
        ("zenodo.org", zenodo("GRCh38_v1.8.tar.gz")),
    ];
    let same = vec![
        ("IMGTHLA/Latest", text("# version: IPD-IMGT/HLA 3.65.0\n")),
        ("IPDKIR/Latest", text("# version: IPD-KIR 2.15.0\n")),
        ("snp/archive", text("<a href=\"b157/\">b157/</a>")),
        ("zenodo.org", zenodo("GRCh38_v1.7.tar.gz")),
    ];
    let broken = vec![
        ("IMGTHLA/Latest", text("# version: something else\n")),
        ("IPDKIR/Latest", None),
        ("snp/archive", text("<html>no builds</html>")),
        ("zenodo.org", zenodo("CHM13_v1.16.tar.gz")),
    ];
    let (st, e) = data_rows(newer, Some("v2.4.0"));
    check(
        st == ["behind"; 5] && e == 0,
        &format!("HLA, KIR, pgsc_calc, dbSNP and Parascopy: a newer upstream is behind ({st:?})"),
    );
    let (st, e) = data_rows(same, Some("v2.3.0"));
    check(
        st == ["current"; 5] && e == 0,
        &format!("the same release upstream is current ({st:?})"),
    );
    let (st, e) = data_rows(broken, None);
    check(
        st == ["ERROR"; 5] && e == 5,
        &format!("an answer without the release is an error row for each of the five ({st:?})"),
    );

    // 2. A fake server: an empty or truncated 200 must fail, and the issue
    //    logic.
    let gh = Arc::new(Mutex::new(FakeGitHub::default()));
    let base = match serve(gh.clone()) {
        Ok(b) => b,
        Err(e) => {
            check(false, &format!("the fake GitHub starts ({e})"));
            let _ = writeln!(out, "SELF-TEST FAILED");
            return 1;
        }
    };
    let client = Client::new(0);
    let (empty_url, trunc_url) = (format!("{base}/empty"), format!("{base}/truncated"));
    let mut empty = Request::get(&empty_url, crate::http::Want::Json);
    empty.retries = Some(0);
    check(
        client.request(&empty).is_err(),
        "an empty 200 answer fails closed",
    );
    let mut trunc = Request::get(&trunc_url, crate::http::Want::Json);
    trunc.retries = Some(0);
    check(
        client.request(&trunc).is_err(),
        "a truncated answer fails closed",
    );
    let mut late = Request::get("https://192.0.2.1/", crate::http::Want::Json);
    late.retries = Some(0);
    late.timeout = Some(5);
    check(
        client.request(&late).is_err(),
        "no answer (timeout) fails closed",
    );

    let mut c = Client::new(crate::http::RETRIES);
    c.backoff = 0;
    c.github = GitHub {
        api: base.clone(),
        token: Some("fake".into()),
    };
    {
        let mut g = gh.lock().expect("fake github");
        g.issues = vec![
            [
                ("number", Value::Int(7)),
                ("state", s("open")),
                ("title", s("other")),
                ("body", s("image test failed")),
                ("labels", Value::Array(vec![s(LABEL)])),
                ("html_url", s("https://x/7")),
            ]
            .into_iter()
            .map(|(k, v)| (k.to_owned(), v))
            .collect(),
        ];
    }
    let a1 = issue::update(&c, "o/r", "report body\n");
    let a2 = issue::update(&c, "o/r", "report body\n");
    let mine = |g: &FakeGitHub| -> Vec<i128> {
        g.issues
            .iter()
            .filter(|i| matches!(i.get("body"), Some(Value::Str(b)) if b.contains(MARKER)))
            .filter_map(|i| match i.get("number") {
                Some(Value::Int(n)) => Some(*n),
                _ => None,
            })
            .collect()
    };
    let n1 = a1.as_ref().map(|x| x.1).unwrap_or(-1);
    {
        let g = gh.lock().expect("fake github");
        check(
            a1.as_ref().map(|x| x.0).ok() == Some(Action::Created)
                && a2.ok() == Some((Action::Edited, n1))
                && mine(&g).len() == 1,
            "first run creates the issue, second run edits the same one",
        );
        check(
            g.issues.first().and_then(|i| i.get("body")) == Some(&s("image test failed")),
            "a freshness-labelled issue without the marker is untouched",
        );
        check(
            g.labels.iter().any(|l| l == LABEL),
            "a missing label is created",
        );
    }
    {
        let mut g = gh.lock().expect("fake github");
        if let Some(i) = g
            .issues
            .iter_mut()
            .find(|i| i.get("number") == Some(&Value::Int(n1)))
        {
            i.insert("state".into(), s("closed"));
        }
    }
    let a3 = issue::update(&c, "o/r", "report body\n");
    {
        let g = gh.lock().expect("fake github");
        let reopened = g
            .issues
            .iter()
            .find(|i| i.get("number") == Some(&Value::Int(n1)))
            .and_then(|i| i.get("state"))
            == Some(&s("open"));
        check(
            a3.ok() == Some((Action::ReopenedAndEdited, n1)) && reopened,
            "a closed report is reopened, not duplicated",
        );
    }
    {
        let mut g = gh.lock().expect("fake github");
        g.issues.clear();
        g.lose_create_answer = true;
    }
    let lost = issue::update(&c, "o/r", "report body\n");
    {
        let g = gh.lock().expect("fake github");
        let n = mine(&g).len();
        check(
            matches!(&lost, Err(LookupError::Failed(m)) if m.contains("HTTP 502")) && n == 1,
            &format!(
                "a create whose answer is lost fails the run and is not sent twice ({n} issue(s))"
            ),
        );
    }
    let mut no_token = Client::new(0);
    no_token.github = GitHub {
        api: base,
        token: None,
    };
    check(
        issue::update(&no_token, "o/r", "x").is_err(),
        "--update-issue without GITHUB_TOKEN fails",
    );

    // 3. Controls on the real registries: a bogus image name, a pinned tag
    //    that does not exist and a bogus ClinVar URL must make the run exit 1,
    //    and the real image beside them must not be an error.
    let scratch = Scratch::new().expect("a private scratch folder");
    let d = scratch.path();
    let (ve, se, o) = (d.join("versions.env"), d.join("setup.sh"), d.join("out.md"));
    let _ = std::fs::write(
        &ve,
        "SAMTOOLS_IMAGE=\"staphb/samtools:1.20\"\n\
         BOGUS_IMAGE=\"quay.io/biocontainers/pgp-freshness-no-such-image:1.0--h0_0\"\n\
         BADTAG_IMAGE=\"staphb/samtools:0.0.0-pgp-no-such-tag\"\n",
    );
    let _ = std::fs::write(
        &se,
        "CLINVAR_URL=\"https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/pgp-no-such-file.vcf.gz\"\n",
    );
    let args: Vec<String> = [
        "--versions",
        &ve.to_string_lossy(),
        "--setup",
        &se.to_string_lossy(),
        "--sections",
        "clinvar,images",
        "--out",
        &o.to_string_lossy(),
    ]
    .iter()
    .map(|x| (*x).to_owned())
    .collect();
    let (mut so, mut se_) = (Vec::new(), Vec::new());
    let rc = crate::cli::run(&args, &mut so, &mut se_);
    let report = std::fs::read_to_string(&o).unwrap_or_default();
    drop(scratch);
    let shown = String::from_utf8_lossy(&so).into_owned() + &String::from_utf8_lossy(&se_);
    let mut lines = String::new();
    for l in shown.lines() {
        lines += &format!("     | {l}\n");
    }
    check(
        rc == 1,
        &format!("the run with a bogus image and a bogus URL exits 1\n{lines}"),
    );
    let errors = report.split_once("### Errors").map_or("", |x| x.1);
    check(
        errors.contains("BOGUS_IMAGE"),
        "the bogus image is an error",
    );
    check(
        errors.contains("BADTAG_IMAGE") && errors.contains("pinned tag is not in"),
        "a missing pinned tag is an error",
    );
    check(
        errors.contains("ClinVar") && errors.contains("HTTP 404"),
        "the bogus ClinVar URL is an error",
    );
    check(
        !errors.contains("SAMTOOLS_IMAGE") && report.contains("SAMTOOLS_IMAGE"),
        "the real image is looked up without error",
    );
    let fails = check(true, "every check ran");
    let _ = writeln!(
        out,
        "SELF-TEST {}",
        if fails == 0 {
            "OK".to_owned()
        } else {
            format!("FAILED: {fails} check(s)")
        }
    );
    u8::from(fails != 0)
}
