//! The live site's plugins, compiled into og: each adds a part of
//! `state.json` (`plugins.<name>`), routes, files of a run folder the
//! browser may fetch, or a change to the served report. Their page scripts
//! are `assets/live/plugins/<name>.js`, which call
//! `window.PGP.register({name, onState, onEvent, panel})`.
//!
//! | plugin  | state.json | routes | run files | report |
//! |---------|------------|--------|-----------|--------|
//! | charts  | [`crate::charts`] | | `*_charts.json` | |
//! | family  | [`crate::family`] | `/family/<family>/<file>`, `/runs/<family>/family.json` | | |
//! | map     | | [`crate::map`] | | |
//! | ops     | [`crate::ops`] | `/metrics` | | |
//! | recheck | [`crate::recheck`] | | `recheck/state.json` | |
//! | results | | `/runs/<run>/results.json`, `/runs/<run>/finding/<id>/` | | |
//! | search  | [`crate::search`] | | `*_search_index.json` | |
//! | viewer  | | `/runs/<run>/viewer/...` ([`crate::viewer`]) | | its links point here |
//!
//! The results are a seam ([`ResultsFns`]): `og live results` fills it when
//! it lands. Until then `results.json` carries a note and no card, a
//! finding's page is 404 and the map is not built.

use crate::family::{self, inside, realpath};
use crate::ops::{self, Asker, Ctx};
use crate::runs::{self, Run};
use crate::server::{Req, quote};
use crate::state::{Live, get, now, obj, s};
use crate::static_files as sf;
use crate::viewer::{self, LIVE_CSP, unquote};
use indexmap::IndexMap;
use og_runlog::host::{Host, System};
use og_runlog::state::State;
use pycompat::json::{Dump, Value, dumps};
use std::collections::BTreeSet;
use std::fs::{self, File};
use std::io::{self, Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};

/// The files of a run folder the plugins let the browser fetch, as
/// `fnmatch` patterns.
pub const RUN_GLOBS: [&str; 3] = ["*_charts.json", "recheck/state.json", "*_search_index.json"];

/// `results.json` of run `name` from its summary, given the files the
/// run's report links (`None`: any).
pub type BuildFn = fn(&str, &Value, Option<&BTreeSet<String>>) -> Result<Value, String>;

/// Each person's results: what `og live results` gives the site.
#[derive(Clone, Copy)]
pub struct ResultsFns {
    /// `results.json` of run `name` from its summary: the chapters, the
    /// cards and the findings (`[{id, level, ...}]`). The files the run's
    /// report links (`None`: any).
    pub build: BuildFn,
    /// The home page's overview of a summary: `attention`, `levels`, `top`,
    /// `lead`, `lead_closed`, `screened`.
    pub overview: fn(&Value) -> Result<Value, String>,
}

/// The note of `results.json` while the site has no results.
pub const NO_RESULTS_NOTE: &str =
    "The results are not on og live serve yet: the report and the Progress tab still work.";
/// The note of `results.json` while no step has finished.
pub const NO_STEP_NOTE: &str =
    "No step has finished yet: the results appear here as the steps end.";

const TEXT: &str = "text/plain; charset=utf-8";

/// A finding's id: a lower-case letter or digit, then up to 79 of them or
/// `-`.
pub fn valid_finding_id(id: &str) -> bool {
    let b = id.as_bytes();
    let ok = |c: &u8| c.is_ascii_lowercase() || c.is_ascii_digit();
    !b.is_empty() && b.len() <= 80 && ok(&b[0]) && b[1..].iter().all(|c| ok(c) || *c == b'-')
}

// --- state.json ---------------------------------------------------------------------------

/// `plugins` of `state.json`: each plugin's part, in the order of the
/// plugins' names.
pub fn state_hooks(
    live: &Live,
    run: &Run,
    d: &Path,
    st: &State<'_>,
    _v: &Value,
    events: &[Value],
) -> Value {
    let mut o = IndexMap::new();
    o.insert("charts".into(), crate::charts::state(d));
    o.insert("family".into(), family_hook(live, d));
    o.insert("ops".into(), ops_hook(live, d, st, events));
    o.insert("recheck".into(), crate::recheck::state(d));
    o.insert("search".into(), crate::search::state(d));
    let _ = run;
    Value::Object(o)
}

/// `plugins.ops`: the run health card, or the line it shows instead.
pub fn ops_hook(live: &Live, d: &Path, st: &State<'_>, events: &[Value]) -> Value {
    let host = System;
    let ctx = Ctx {
        now: now(),
        steps: &live.steps,
        host: &host,
        tz: &live.tz,
    };
    ops::hook(Some(d), st, events, &live.cfg.ops, &ctx).to_json()
}

fn folder_name(d: &Path) -> String {
    d.components()
        .next_back()
        .map(|c| c.as_os_str().to_string_lossy().into_owned())
        .unwrap_or_default()
}

/// `plugins.family`: the family's state on a family's page; on a person's
/// page the families that list it and whether it is a child
/// ([`member_state`]).
pub fn family_hook(live: &Live, d: &Path) -> Value {
    let name = folder_name(d);
    if family::is_family(d) && inside(d, &live.root) {
        return match family::family_state(&name, d) {
            Ok(v) => v,
            Err(e) => {
                eprintln!(
                    "og live serve: the family state of {name} could not be read: {}",
                    e.0
                );
                obj([("error", s("the family's files cannot be read"))])
            }
        };
    }
    member_state(&live.root, &name)
}

/// On the page of sample `sample`: the families that list it, and whether
/// its pages treat it as a child. With no readable family listing it while
/// a family cannot be read, the answer is unknown, never "in no family".
pub fn member_state(root: &Path, sample: &str) -> Value {
    let fams = family::families(root);
    let child = runs::child_of(&fams, sample);
    let mine: Vec<Value> = fams
        .iter()
        .filter(|f| match get(f, "members") {
            Some(Value::Array(ms)) => ms
                .iter()
                .any(|m| matches!(get(m, "sample"), Some(Value::Str(x)) if x == sample)),
            _ => false,
        })
        .map(|f| {
            let n = get(f, "name").cloned().unwrap_or(Value::Null);
            let url = Value::Str(format!(
                "/runs/{}/",
                match &n {
                    Value::Str(x) => x.clone(),
                    _ => String::new(),
                }
            ));
            obj([("name", n), ("url", url)])
        })
        .collect();
    if !mine.is_empty() {
        return obj([
            ("kind", s("member")),
            ("families", Value::Array(mine)),
            ("child", Value::Bool(child)),
        ]);
    }
    if child {
        obj([
            ("kind", s("member")),
            ("families", Value::Array(Vec::new())),
            ("unknown", Value::Bool(true)),
        ])
    } else {
        Value::Null
    }
}

// --- the report ---------------------------------------------------------------------------

/// The served report: its viewer links and note pointed at this site's
/// genome viewer of the run, when the viewer is set up.
pub fn report_filter(live: &Live, run: &str, html: String) -> String {
    if viewer_ready(live).is_some() {
        return html;
    }
    let here = format!("/runs/{}/viewer/", quote(run));
    html.replace(
        &format!("href=\"{}#", og_render::html::VIEWER_URL),
        &format!("href=\"{here}#"),
    )
    .replace(og_render::notes::VIEWER_NOTE, LIVE_NOTE)
}

/// The note the served report carries for the site's viewer.
pub const LIVE_NOTE: &str = "    <p class=\"viewer\">The Viewer column opens the reads at each row in this site's genome viewer \
     (docs/live-site.md).</p>";

// --- the routes ---------------------------------------------------------------------------

/// `/runs/<name>/<tail>` with a valid run name: `(name, tail)`.
fn run_route<'p>(path: &'p str, tail: &str) -> Option<&'p str> {
    let rest = path.strip_prefix("/runs/")?;
    let name = rest.strip_suffix(tail)?;
    runs::valid_run(name).then_some(name)
}

/// Answers `path` when a plugin has a route for it: `Ok(true)` once
/// answered, `Ok(false)` for a path no plugin takes (or a method its route
/// does not).
pub fn route(req: &mut Req<'_>, path: &str) -> io::Result<bool> {
    let get_ = req.method == "GET";
    // family
    if let Some(rest) = path.strip_prefix("/family/")
        && let Some((name, f)) = rest.split_once('/')
        && runs::valid_run(name)
        && runs::valid_run(f)
    {
        if !get_ {
            return Ok(false);
        }
        family_file(req, name, f)?;
        return Ok(true);
    }
    if let Some(name) = run_route(path, "/family.json") {
        if !get_ {
            return Ok(false);
        }
        family_json(req, name)?;
        return Ok(true);
    }
    // map
    if let Some(name) = run_route(path, "/map.json") {
        if !get_ {
            return Ok(false);
        }
        crate::map::build(req, name)?;
        return Ok(true);
    }
    if let Some(which) = path
        .strip_prefix("/map/genes/")
        .and_then(|w| w.strip_suffix(".json"))
        .filter(|w| crate::map::valid_genes(w))
    {
        if !get_ {
            return Ok(false);
        }
        crate::map::genes(req, which)?;
        return Ok(true);
    }
    // ops
    if path == "/metrics" {
        if !get_ {
            return Ok(false);
        }
        metrics(req)?;
        return Ok(true);
    }
    // results
    if let Some(name) = run_route(path, "/results.json") {
        if !get_ {
            return Ok(false);
        }
        results(req, name)?;
        return Ok(true);
    }
    if let Some(rest) = path.strip_prefix("/runs/")
        && let Some((name, tail)) = rest.split_once("/finding/")
        && runs::valid_run(name)
        && let Some(id) = tail.strip_suffix('/')
        && !id.is_empty()
        && id.chars().count() <= 200
        && !id.contains('/')
    {
        if !get_ {
            return Ok(false);
        }
        finding(req, name, id)?;
        return Ok(true);
    }
    // viewer
    if let Some(rest) = path.strip_prefix("/runs/") {
        if let Some(run) = rest
            .strip_suffix("/viewer")
            .filter(|r| !r.is_empty() && !r.contains('/'))
        {
            if !get_ {
                return Ok(false);
            }
            viewer_slash(req, run)?;
            return Ok(true);
        }
        if let Some((run, tail)) = rest.split_once('/')
            && !run.is_empty()
            && let Some(file) = tail.strip_prefix("viewer/")
        {
            if !get_ {
                return Ok(false);
            }
            viewer_route(req, run, file)?;
            return Ok(true);
        }
    }
    Ok(false)
}

// --- family -------------------------------------------------------------------------------

/// The family folder of `name`, when it is one inside the root.
fn family_dir(req: &Req<'_>, name: &str) -> Option<PathBuf> {
    req.run_dir(name)
        .filter(|d| family::is_family(d) && inside(d, &req.live.root))
}

/// `GET /family/<family>/<file>`: the pedigree drawn now, or one of the
/// family's files listed in `state.json`; any other name is 404.
fn family_file(req: &mut Req<'_>, name: &str, f: &str) -> io::Result<()> {
    let Some(d) = family_dir(req, name) else {
        req.log("family", name, f, 404);
        return req.send_body(404, TEXT, "not found\n");
    };
    if f == "pedigree.svg" {
        let svg = family::members_of(&d).and_then(|p| og_pedigree::layout::draw(&p, name));
        return match svg {
            Ok(svg) => {
                req.log("family", name, f, 200);
                req.send_body(200, "image/svg+xml", svg)
            }
            Err(_) => {
                req.log("family", name, f, 404);
                req.send_body(404, TEXT, "not found\n")
            }
        };
    }
    let files = family::members_of(&d)
        .map(|p| family::served_files(&d, &p))
        .unwrap_or_default();
    let p = d.join(f);
    let body = (files.iter().any(|x| x == f) && inside(&p, &d))
        .then(|| fs::read(&p).ok())
        .flatten();
    let Some(body) = body else {
        req.log("family", name, f, 404);
        return req.send_body(404, TEXT, "not found\n");
    };
    let ctype = match sf::ext(f) {
        ".svg" => "image/svg+xml",
        _ => TEXT,
    };
    req.log("family", name, f, 200);
    req.send_body(200, ctype, body)
}

/// `GET /runs/<family>/family.json`: the family page's data, with what
/// only the site knows: each member run's progress, attention count and
/// whether its pages treat it as a child.
fn family_json(req: &mut Req<'_>, name: &str) -> io::Result<()> {
    let Some(d) = family_dir(req, name) else {
        req.log("family", name, "family.json", 404);
        return req.send_body(404, TEXT, "not found\n");
    };
    match family_page(req.live, name, &d) {
        Ok(v) => {
            req.log("family", name, "family.json", 200);
            req.send_json(&v)
        }
        Err(kind) => {
            eprintln!("og live serve: family.json of {name} not built: {kind}");
            req.log("family", name, "family.json", 500);
            req.send_body(500, TEXT, "family not readable\n")
        }
    }
}

fn family_page(live: &Live, name: &str, d: &Path) -> Result<Value, String> {
    let st = family::family_state(name, d).map_err(|e| e.0)?;
    let listed = live.runs();
    let (mut runs_o, mut candidates, mut summaries) =
        (IndexMap::new(), IndexMap::new(), IndexMap::new());
    let members = match get(&st, "members") {
        Some(Value::Array(m)) => m.clone(),
        _ => Vec::new(),
    };
    for p in &members {
        let Some(Value::Str(sample)) = get(p, "sample") else {
            continue;
        };
        let Some(r) = listed.get(sample) else {
            continue;
        };
        let rd = r.folder(&live.root);
        let rd_s = Value::Str(rd.to_string_lossy().into_owned());
        summaries.insert(sample.clone(), rd_s.clone());
        if r.is_family() {
            continue;
        }
        candidates.insert(sample.clone(), rd_s);
        let (bst, events) = live.base_state(r);
        let v = live.base_json(&bst, events.len());
        let planned: Vec<_> = bst
            .steps
            .iter()
            .filter(|x| x.state.as_str() != "not planned")
            .collect();
        let done = planned
            .iter()
            .filter(|x| x.state.as_str() == "done")
            .count();
        let (mut n, mtime) = live.attention(&rd);
        if live.gated(sample, &rd, &v, Some(&listed)) {
            n = None;
        }
        runs_o.insert(
            sample.clone(),
            obj([
                ("url", Value::Str(format!("/runs/{sample}/"))),
                ("state", s(bst.run.as_str())),
                ("planned", Value::Int(planned.len() as i128)),
                ("done", Value::Int(done as i128)),
                ("attention", n.unwrap_or(Value::Null)),
                ("results", Value::Bool(mtime.is_some())),
                ("child", Value::Bool(live.child(sample, Some(&listed)))),
            ]),
        );
    }
    let input = obj([
        ("runs", Value::Object(runs_o)),
        ("candidates", Value::Object(candidates)),
        ("summaries", Value::Object(summaries)),
    ]);
    family::site_page(name, d, &live.root, &input).map_err(|e| e.0)
}

// --- ops ----------------------------------------------------------------------------------

/// `GET /metrics`: every run's health in Prometheus' text format. Without a
/// password the site is loopback only and `/metrics` is open, like every
/// page; with one it needs a login or `--metrics-token-file`'s bearer token.
fn metrics(req: &mut Req<'_>) -> io::Result<()> {
    let who = Asker {
        login_required: req.live.cfg.auth.hash.is_some(),
        logged_in: req.sid.is_some(),
        authorization: req.header("authorization"),
    };
    let token = req.live.cfg.metrics_token_file.clone();
    match ops::metrics_allowed(&who, || ops::read_token(token.as_deref())) {
        Err(e) => {
            eprintln!(
                "og live serve: --metrics-token-file cannot be read: {}",
                e.kind()
            );
            req.log("metrics", "-", "-", 503);
            return req.send_body(503, TEXT, "the metrics token file cannot be read\n");
        }
        Ok(false) => {
            req.log("metrics", "-", "-", 401);
            return req.send_body(401, TEXT, "login or bearer token required\n");
        }
        Ok(true) => {}
    }
    let live = req.live;
    let mut cards = Vec::new();
    let mut errors = 0;
    for (name, r) in live.runs() {
        let (st, events) = live.base_state(&r);
        match ops::hook(
            Some(&r.folder(&live.root)),
            &st,
            &events,
            &live.cfg.ops,
            &Ctx {
                now: now(),
                steps: &live.steps,
                host: &System,
                tz: &live.tz,
            },
        ) {
            ops::OpsState::Ok(o) => cards.push((name, o)),
            ops::OpsState::Error(_) => {
                errors += 1;
                eprintln!("og live serve: {name}: run health not read");
            }
        }
    }
    let refs: Vec<(&str, &ops::RunOps)> = cards
        .iter()
        .map(|(n, o)| (n.as_str(), o.as_ref()))
        .collect();
    let text = ops::metrics_text(
        &refs,
        ops::disk_warn_ratio(&live.cfg.ops),
        System.physical_memory(),
        errors,
    );
    req.log("metrics", "-", "-", 200);
    req.send_body(200, ops::CONTENT_TYPE, text)
}

// --- results ------------------------------------------------------------------------------

/// What `GET /runs/<run>/results.json` sends, built the same way for the
/// map: the results of the run's summary, or a note while no step is done.
pub fn results_body(live: &Live, run: &Run, fns: &ResultsFns) -> Result<Value, String> {
    let d = run.folder(&live.root);
    let (st, events) = live.base_state(run);
    let v = live.base_json(&st, events.len());
    match live.summary(&run.name, &d, &v) {
        Some(sm) => {
            let files = if run.dir.is_some() {
                runs::run_files(&run.name, run.dir.as_deref())
            } else {
                BTreeSet::new()
            };
            (fns.build)(&run.name, &sm, Some(&files))
        }
        None => Ok(note(&run.name, NO_STEP_NOTE)),
    }
}

fn note(name: &str, text: &str) -> Value {
    obj([
        ("schema_version", Value::Int(1)),
        ("sample", s(name)),
        ("chapters", Value::Array(Vec::new())),
        ("note", s(text)),
    ])
}

/// `GET /runs/<run>/results.json`: the person's results page data. A reader
/// error gives no cards and a note, never a 500.
fn results(req: &mut Req<'_>, name: &str) -> io::Result<()> {
    let live = req.live;
    let Some(run) = live.run(name).filter(|r| !r.is_family()) else {
        req.log("results", name, "-", 404);
        return req.send_body(404, TEXT, "not found\n");
    };
    let out = match live.cfg.results {
        None => note(&run.name, NO_RESULTS_NOTE),
        Some(fns) => {
            let d = run.folder(&live.root);
            let (st, events) = live.base_state(&run);
            let v = live.base_json(&st, events.len());
            match results_body(live, &run, &fns) {
                Ok(mut out) => {
                    // A child whose report the site keeps closed: the
                    // flagged-card count of summary.json, which the family
                    // page shows once every result to act on is opened.
                    if live.summary(&run.name, &d, &v).is_some()
                        && live.gated(&run.name, &d, &v, None)
                        && let Value::Object(o) = &mut out
                    {
                        o.insert(
                            "saved_attention".into(),
                            live.attention(&d).0.unwrap_or(Value::Null),
                        );
                    }
                    out
                }
                Err(_) => {
                    eprintln!("og live serve: results of {} not built", run.name);
                    note(
                        &run.name,
                        "The results could not be read; the report and the Progress tab still work.",
                    )
                }
            }
        }
    };
    req.log("results", &run.name, "-", 200);
    req.send_json(&out)
}

/// `GET /runs/<run>/finding/<id>/`: one finding's page, 404 unless `id` is
/// a finding of the run's results now; for a child with a result to act on,
/// the page for any well-formed id (a 404 would tell a guessed id apart
/// from a result the child has).
fn finding(req: &mut Req<'_>, name: &str, id: &str) -> io::Result<()> {
    let live = req.live;
    let run = live.run(name).filter(|r| !r.is_family());
    let mut ok = false;
    if let (Some(run), Some(fns)) = (&run, live.cfg.results)
        && valid_finding_id(id)
    {
        let d = run.folder(&live.root);
        let (st, events) = live.base_state(run);
        let v = live.base_json(&st, events.len());
        let found = match live.summary(&run.name, &d, &v) {
            Some(sm) => (fns.build)(&run.name, &sm, None).map(|b| match get(&b, "findings") {
                Some(Value::Array(f)) => f.clone(),
                _ => Vec::new(),
            }),
            None => Ok(Vec::new()),
        };
        match found {
            Ok(found) => {
                let is = |f: &Value, k: &str, want: &str| matches!(get(f, k), Some(Value::Str(x)) if x == want);
                ok = found.iter().any(|f| is(f, "id", id))
                    || (found.iter().any(|f| is(f, "level", "act")) && live.child(&run.name, None));
            }
            Err(e) => eprintln!("og live serve: finding of {} not checked: {e}", run.name),
        }
    }
    let Some(run) = run.filter(|_| ok) else {
        req.log("finding", name, "-", 404);
        return req.send_body(404, TEXT, "not found\n");
    };
    match sf::finding_page(&live.cfg.assets, &run.name, id) {
        Ok(page) => {
            req.log("finding", &run.name, "-", 200);
            req.send_body(200, "text/html; charset=utf-8", page)
        }
        Err(_) => req.deny(500, None),
    }
}

// --- viewer -------------------------------------------------------------------------------

const LOCAL_TEXT: &str = "served from this computer only.";
const LIVE_TEXT: &str = "served by the run site, behind its login.";

/// Why the viewer cannot run, or `None` when it can.
pub fn viewer_ready(live: &Live) -> Option<&'static str> {
    let (Some(fasta), Some(igv)) = (&live.cfg.viewer_fasta, &live.cfg.viewer_igv_js) else {
        return Some(
            "the genome viewer is not set up: start the site with --viewer-fasta and --viewer-igv-js \
             (scripts/owngenome-live.sh passes both)",
        );
    };
    let mut fai = fasta.clone().into_os_string();
    fai.push(".fai");
    if !fasta.is_file() || !Path::new(&fai).is_file() {
        return Some("the genome viewer's reference or its .fai is missing");
    }
    if !igv.is_file() {
        return Some("igv.min.js is missing");
    }
    None
}

/// The viewer's headers: same-origin resources, its own policy, framing by
/// this site only.
fn viewer_head(
    req: &mut Req<'_>,
    code: u16,
    ctype: &str,
    len: u64,
    extra: &[(&str, String)],
) -> io::Result<()> {
    let mut all = vec![("Cross-Origin-Resource-Policy", "same-origin".to_owned())];
    all.extend(extra.iter().cloned());
    req.head(
        code,
        ctype,
        Some(len),
        &all,
        LIVE_CSP,
        crate::csp::SAME_ORIGIN,
    )
}

fn viewer_send(req: &mut Req<'_>, code: u16, ctype: &str, body: &[u8]) -> io::Result<()> {
    viewer_head(req, code, ctype, body.len() as u64, &[])?;
    req.body(body)
}

/// `src`, or the one byte range the request asks for (206, or 416 past the
/// end), under the viewer's headers.
fn viewer_file(req: &mut Req<'_>, ctype: &str, src: &Path) -> io::Result<()> {
    let mut f = File::open(src)?;
    let size = f.metadata()?.len();
    let rng = viewer::parse_range(req.header("range"), size);
    if rng == viewer::Range::Bad {
        return viewer_head(
            req,
            416,
            TEXT,
            0,
            &[("Content-Range", format!("bytes */{size}"))],
        );
    }
    let (start, end) = match rng {
        viewer::Range::Part(a, b) => (a, b),
        _ => (0, size.saturating_sub(1)),
    };
    let length = if size == 0 { 0 } else { end - start + 1 };
    let mut extra = vec![("Accept-Ranges", "bytes".to_owned())];
    if rng != viewer::Range::Whole {
        extra.push(("Content-Range", format!("bytes {start}-{end}/{size}")));
    }
    let code = if rng == viewer::Range::Whole {
        200
    } else {
        206
    };
    viewer_head(req, code, ctype, length, &extra)?;
    if req.head_only || length == 0 {
        return Ok(());
    }
    f.seek(SeekFrom::Start(start))?;
    let mut left = length;
    let mut buf = vec![0u8; (1usize << 20).min(length as usize)];
    while left > 0 {
        let want = buf.len().min(left as usize);
        let n = f.read(&mut buf[..want])?;
        if n == 0 {
            break;
        }
        req.out().write_all(&buf[..n])?;
        left -= n as u64;
    }
    req.out().flush()
}

/// `(config, files)` of one run's viewer: the run folder's reads, VCF and SV
/// VCF that lie inside the run and the root, each track with its index
/// there too. Files 0 and 1 are the reference and its `.fai`, from any
/// folder.
fn run_config(live: &Live, run: &str, d: &Path) -> Result<(Value, Vec<Option<PathBuf>>), String> {
    let top = realpath(d);
    let root = live.root.clone();
    let ok = |p: &Path| {
        let real = realpath(p);
        [&top, &root].iter().all(|t| {
            let mut pre = t.to_string_lossy().trim_end_matches('/').to_owned();
            pre.push('/');
            real.to_string_lossy().starts_with(&pre)
        })
    };
    let sfiles = viewer::sample_files(d, run);
    let keep = |p: Option<PathBuf>| p.filter(|p| ok(p));
    let fasta = live.cfg.viewer_fasta.clone().unwrap_or_default();
    let cfg = viewer::build_config(
        run,
        &fasta,
        keep(sfiles.reads).as_deref(),
        keep(sfiles.vcf).as_deref(),
        keep(sfiles.sv_vcf).as_deref(),
        None,
    )?;
    let mut files: Vec<Option<PathBuf>> = cfg
        .files
        .iter()
        .enumerate()
        .map(|(n, f)| (n < 2 || ok(f)).then(|| f.clone()))
        .collect();
    let mut json = cfg.json;
    let nums = |t: &Value| -> Vec<usize> {
        ["url", "indexURL"]
            .iter()
            .filter_map(|k| match get(t, k) {
                Some(Value::Str(u)) => u.split('/').nth(1).and_then(|n| n.parse().ok()),
                _ => None,
            })
            .collect()
    };
    if let Value::Object(o) = &mut json
        && let Some(Value::Array(tracks)) = o.get_mut("tracks")
    {
        for t in tracks.iter() {
            let ns = nums(t);
            if !ns
                .iter()
                .all(|&n| files.get(n).is_some_and(Option::is_some))
            {
                for n in ns {
                    if let Some(f) = files.get_mut(n) {
                        *f = None;
                    }
                }
            }
        }
        tracks.retain(|t| {
            nums(t)
                .iter()
                .all(|&n| files.get(n).is_some_and(Option::is_some))
        });
    }
    Ok((json, files))
}

/// `GET /runs/<run>/viewer`: to the page.
fn viewer_slash(req: &mut Req<'_>, raw: &str) -> io::Result<()> {
    let run = unquote(raw);
    if req.run_dir(&run).is_none() {
        return req.deny(404, None);
    }
    let to = format!("/runs/{}/viewer/", quote(&run));
    req.redirect(301, &to)
}

/// `GET /runs/<run>/viewer/<file>`: the page, its script, igv.js, the
/// config and the config's files by byte range.
fn viewer_route(req: &mut Req<'_>, raw: &str, rest: &str) -> io::Result<()> {
    let run = unquote(raw);
    let live = req.live;
    let page = rest.is_empty() || rest == "index.html";
    let Some(d) = req.run_dir(&run) else {
        req.log("viewer", "-", "-", 404);
        return req.deny(404, None);
    };
    if let Some(why) = viewer_ready(live) {
        req.log("viewer", &run, "-", 404);
        return viewer_send(req, 404, TEXT, why.as_bytes());
    }
    let dir = live.cfg.repo_assets.join("viewer");
    if page {
        let Ok(html) = fs::read_to_string(dir.join("index.html")) else {
            return req.deny(404, None);
        };
        req.log("viewer", &run, "-", 200);
        return viewer_send(
            req,
            200,
            "text/html; charset=utf-8",
            html.replace(LOCAL_TEXT, LIVE_TEXT).as_bytes(),
        );
    }
    if rest == "viewer.js" {
        let Ok(js) = fs::read(dir.join("viewer.js")) else {
            return req.deny(404, None);
        };
        return viewer_send(req, 200, "text/javascript; charset=utf-8", &js);
    }
    if rest == "igv.min.js" {
        let igv = live.cfg.viewer_igv_js.clone().unwrap_or_default();
        return viewer_file(req, "text/javascript; charset=utf-8", &igv);
    }
    let (cfg, files) = match run_config(live, &run, &d) {
        Ok(x) => x,
        Err(e) => {
            req.log("viewer", &run, "-", 409);
            return req.deny(
                409,
                Some(&format!("the genome viewer cannot open this run: {e}")),
            );
        }
    };
    if rest == "config.json" {
        let text = dumps(
            &cfg,
            &Dump {
                indent: Some(1),
                ..Dump::default()
            },
        )
        .unwrap_or_default();
        return viewer_send(req, 200, "application/json", text.as_bytes());
    }
    if let Some(dm) = rest.strip_prefix("data/")
        && let Some((n, fname)) = dm.split_once('/')
        && !n.is_empty()
        && n.bytes().all(|b| b.is_ascii_digit())
        && !fname.is_empty()
        && let Ok(n) = n.parse::<usize>()
        && let Some(Some(f)) = files.get(n)
        && f.file_name().map(|x| x.to_string_lossy().into_owned()) == Some(unquote(fname))
    {
        let f = f.clone();
        return viewer_file(req, "application/octet-stream", &f);
    }
    req.deny(404, None)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finding_ids() {
        assert!(valid_finding_id("clinvar-brca1-1") && valid_finding_id("0"));
        assert!(!valid_finding_id("") && !valid_finding_id("-a") && !valid_finding_id("A"));
        assert!(!valid_finding_id(&"a".repeat(81)) && valid_finding_id(&"a".repeat(80)));
    }

    #[test]
    fn run_routes() {
        assert_eq!(run_route("/runs/S1/map.json", "/map.json"), Some("S1"));
        assert_eq!(run_route("/runs/a/b/map.json", "/map.json"), None);
        assert_eq!(run_route("/runs/.x/map.json", "/map.json"), None);
    }
}
