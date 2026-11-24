//! The live site's run state: what `state.json`, `/runs.json` and the
//! event stream of `og live serve` say about each run.
//!
//! - [`Store`]: the stored messages of each run, `STATE_DIR/runs/<run>/
//!   events.jsonl`, read again whenever the file changes, so a message
//!   another process appends reaches every open stream.
//! - [`messages`], [`current_run`], [`run_finished`]: what a browser is sent
//!   of each message, and which messages are the newest run's.
//! - [`headline`]: one line per finished section, from og-render's own
//!   helpers, so it reads as the report does.
//! - [`Live`]: the state every request shares: the runs, each run's state
//!   ([`Live::base_state`], [`Live::state`]), the summary the headlines read
//!   and the home page's list ([`Live::runs_index`]).
//!
//! No path, task tag or parameter of a message reaches a page.

use crate::family::{inside, realpath};
use crate::runs::{self, Run};
use crate::server::{AccessLog, Config};
use indexmap::IndexMap;
use jiff::tz::TimeZone;
use og_render::model::{Malformed, Summary};
use og_render::sections::clinvar::D;
use og_runlog::files;
use og_runlog::state::{Inputs, RunState, State, trace_files};
use og_runlog::steps::Step;
use og_runlog::tasks;
use pycompat::json::{self, Dump, Value};
use pycompat::value::{py_str, truthy};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::{Arc, Mutex, MutexGuard, PoisonError};
use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};

pub type Obj = IndexMap<String, Value>;

/// The lock's data, also after another request panicked holding it.
pub(crate) fn lock<T>(m: &Mutex<T>) -> MutexGuard<'_, T> {
    m.lock().unwrap_or_else(PoisonError::into_inner)
}

pub(crate) fn s(text: &str) -> Value {
    Value::Str(text.to_owned())
}

pub(crate) fn obj<const N: usize>(pairs: [(&str, Value); N]) -> Value {
    Value::Object(pairs.into_iter().map(|(k, v)| (k.to_owned(), v)).collect())
}

pub(crate) fn get<'a>(v: &'a Value, key: &str) -> Option<&'a Value> {
    match v {
        Value::Object(o) => o.get(key),
        _ => None,
    }
}

fn get_str<'a>(v: &'a Value, key: &str) -> Option<&'a str> {
    match get(v, key) {
        Some(Value::Str(x)) => Some(x),
        _ => None,
    }
}

/// `json.dumps(value)` with Python's default separators.
pub fn dumps(value: &Value) -> String {
    json::dumps(value, &Dump::default()).unwrap_or_else(|_| "null".into())
}

/// `json.dumps(value, separators=(",", ":"))`.
pub fn dumps_compact(value: &Value) -> String {
    json::dumps(
        value,
        &Dump {
            separators: Some((",", ":")),
            ..Dump::default()
        },
    )
    .unwrap_or_else(|_| "null".into())
}

/// Seconds since the epoch.
pub fn now() -> f64 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_or(0.0, |d| d.as_secs_f64())
}

/// A file's (mtime in ns, size); `None` when it cannot be read.
pub fn file_stamp(p: &Path) -> Option<(i128, u64)> {
    let m = fs::metadata(p).ok()?;
    Some((pycompat::mtime::mtime_ns(p).ok()?, m.len()))
}

// --- the stored messages ------------------------------------------------------------------

struct Cached {
    stamp: Option<(i128, u64)>,
    events: Arc<Vec<Value>>,
}

/// `STATE_DIR/runs/<run>/events.jsonl` of every run, kept in memory and read
/// again when the file's time or size changes.
pub struct Store {
    dir: PathBuf,
    runs: Mutex<HashMap<String, Cached>>,
}

/// How often an open stream looks at its run's file.
const POLL: Duration = Duration::from_millis(200);

impl Store {
    pub fn new(state_dir: &Path) -> Store {
        Store {
            dir: state_dir.join("runs"),
            runs: Mutex::new(HashMap::new()),
        }
    }

    /// The events file of `run`.
    pub fn path(&self, run: &str) -> PathBuf {
        self.dir.join(run).join("events.jsonl")
    }

    /// The runs with an events file, by name.
    pub fn known(&self) -> Vec<String> {
        let mut out: Vec<String> = fs::read_dir(&self.dir)
            .map(|rd| {
                rd.filter_map(|e| e.ok()?.file_name().into_string().ok())
                    .filter(|r| runs::valid_run(r) && self.path(r).is_file())
                    .collect()
            })
            .unwrap_or_default();
        out.sort();
        out
    }

    /// The stored messages of `run`: each line that is a JSON object.
    pub fn events(&self, run: &str) -> Arc<Vec<Value>> {
        let p = self.path(run);
        let stamp = file_stamp(&p);
        {
            let runs = lock(&self.runs);
            if let Some(c) = runs.get(run)
                && c.stamp == stamp
            {
                return Arc::clone(&c.events);
            }
        }
        let events = Arc::new(files::read_events(Some(&p)).unwrap_or_default());
        lock(&self.runs).insert(
            run.to_owned(),
            Cached {
                stamp,
                events: Arc::clone(&events),
            },
        );
        events
    }

    /// The messages of `run` after the first `after`, once there are any or
    /// `timeout` has passed (then none), or `stop` is set.
    pub fn wait(
        &self,
        run: &str,
        after: usize,
        timeout: Duration,
        stop: &AtomicBool,
    ) -> Vec<Value> {
        let end = Instant::now() + timeout;
        loop {
            let evs = self.events(run);
            if evs.len() > after {
                return evs[after..].to_vec();
            }
            let left = end.saturating_duration_since(Instant::now());
            if left.is_zero() || stop.load(Ordering::SeqCst) {
                return Vec::new();
            }
            std::thread::sleep(left.min(POLL));
        }
    }
}

// --- what a message says ------------------------------------------------------------------

fn event_is(ev: &Value, name: &str) -> bool {
    get_str(ev, "event") == Some(name)
}

/// Whether the run is over: run-all.sh posts `pgp_run_finished` as its last
/// message, so a run is over only while that is the newest message (a
/// message after it is a later run of the same sample).
pub fn run_finished(events: &[Value]) -> bool {
    events
        .last()
        .is_some_and(|e| event_is(e, "pgp_run_finished"))
}

/// The messages of the newest run of a sample: those after the last
/// `pgp_run_finished` that is not the newest message, and from the last
/// weblog `started` when that comes later (a bare nextflow run posts no
/// `pgp_run_finished`).
pub fn current_run(events: &[Value]) -> &[Value] {
    let mut start = 0;
    if events.len() > 1 {
        for (i, ev) in events[..events.len() - 1].iter().enumerate() {
            if event_is(ev, "pgp_run_finished") {
                start = i + 1;
            }
        }
    }
    if let Some(i) = (start..events.len())
        .rev()
        .find(|&i| event_is(&events[i], "started"))
    {
        start = i;
    }
    &events[start..]
}

/// `'SV:MANTA (s1)'` as `(MANTA, s1)`: the process without its workflow
/// prefix and the tag (og-nf's reading); `None` for a name of two words.
pub fn task_name(name: &str) -> (Option<String>, Option<String>) {
    match og_nf::name::task_name(name) {
        Some((p, tag)) => (Some(p.to_owned()), tag.map(str::to_owned)),
        None => (None, None),
    }
}

/// Which step runs each process.
pub fn step_of_process(steps: &[Step]) -> HashMap<&str, &str> {
    let mut out = HashMap::new();
    for s in steps {
        for p in &s.processes {
            out.insert(p.as_str(), s.step.as_str());
        }
    }
    out
}

fn text_of(v: &Value) -> String {
    match v {
        Value::Str(x) => x.clone(),
        other => py_str(other).unwrap_or_default(),
    }
}

/// What the browser is sent for each stored message, in order: `(id,
/// payload)`. The id is the message's place in the run's messages (the
/// first of `events` is `start + 1`). No path, parameter or tag leaves here.
pub fn messages(events: &[Value], steps: &[Step], start: usize) -> Vec<(usize, Value)> {
    let owner = step_of_process(steps);
    events
        .iter()
        .enumerate()
        .map(|(i, ev)| {
            let mut p = Obj::new();
            p.insert(
                "event".into(),
                get(ev, "event").cloned().unwrap_or(Value::Null),
            );
            p.insert(
                "utc".into(),
                get(ev, "utcTime").cloned().unwrap_or(Value::Null),
            );
            if let Some(tr @ Value::Object(_)) = get(ev, "trace") {
                let name = [get(tr, "name"), get(tr, "process")]
                    .into_iter()
                    .flatten()
                    .find(|v| truthy(v))
                    .map(text_of)
                    .unwrap_or_default();
                let (proc, _) = task_name(&name);
                let step = proc
                    .as_deref()
                    .and_then(|p| owner.get(p))
                    .map_or(Value::Null, |st| s(st));
                p.insert("process".into(), proc.map_or(Value::Null, Value::Str));
                p.insert("step".into(), step);
                p.insert(
                    "status".into(),
                    get(tr, "status").cloned().unwrap_or(Value::Null),
                );
                p.insert(
                    "attempt".into(),
                    get(tr, "attempt").cloned().unwrap_or(Value::Null),
                );
            }
            if event_is(ev, "pgp_step") {
                p.insert(
                    "step".into(),
                    get(ev, "step").cloned().unwrap_or(Value::Null),
                );
                p.insert(
                    "state".into(),
                    get(ev, "state").cloned().unwrap_or(Value::Null),
                );
            }
            if event_is(ev, "pgp_run_finished") {
                p.insert(
                    "exit".into(),
                    get(ev, "exit").cloned().unwrap_or(Value::Null),
                );
            }
            p.retain(|_, v| *v != Value::Null);
            (start + i + 1, Value::Object(p))
        })
        .collect()
}

// --- headlines ----------------------------------------------------------------------------

/// The first line, and the section's caveat when the helper ends with it
/// (never re-worded), as two sentences the page can cut after the first.
fn with_note(lines: &[String], note: &str) -> Option<String> {
    let first = lines.first()?;
    if lines.last().map(String::as_str) != Some(note) {
        return Some(first.clone());
    }
    let dot = if first.trim_end().ends_with(['.', '!', '?']) {
        ""
    } else {
        "."
    };
    Some(format!("{first}{dot} {note}"))
}

fn first(lines: Vec<String>) -> Result<String, Malformed> {
    lines
        .into_iter()
        .next()
        .ok_or_else(|| Malformed("no line".into()))
}

/// A number Python compares with `> 0`; an error for any other value.
fn positive(v: &Value) -> Result<bool, Malformed> {
    match v {
        Value::Int(n) => Ok(*n > 0),
        Value::Float(f) => Ok(*f > 0.0),
        Value::Bool(b) => Ok(*b),
        _ => Err(Malformed("not a number".into())),
    }
}

fn or_zero(v: Option<&Value>) -> Value {
    v.filter(|v| truthy(v)).cloned().unwrap_or(Value::Int(0))
}

/// The keys [`headline`] has a line of its own for.
pub const HEADLINE_KEYS: [&str; 11] = [
    "clinvar", "carrier", "couple", "roh", "pypgx", "annotsv", "viral", "archaic", "lpa", "cpsr",
    "mito",
];

/// `(text, attention)` of section `key` from its helper, `None` for a key
/// without one. An error is a summary the helper cannot read.
fn line_of(key: &str, sm: &Summary) -> Result<Option<(Option<String>, bool)>, Malformed> {
    use og_render::notes::{COUPLE_NOTE, ROH_NOTE, VIRAL_NOTE};
    use og_render::sections as sec;
    let x = sm.section(key);
    let d = D::of(&x);
    let data = x.data;
    let flag = |k: &str| get(data, k).is_some_and(truthy);
    Ok(Some(match key {
        "clinvar" => {
            let n = or_zero(get(data, "count"));
            (
                Some(format!(
                    "Pathogenic/Likely Pathogenic hits: {}. {}",
                    py_str(&n).map_err(|e| Malformed(e.to_string()))?,
                    og_render::notes::CLINVAR_FREQ_NOTE
                )),
                positive(&n)?,
            )
        }
        "carrier" => (
            Some(first(sec::carrier::lines(&d)?)?),
            flag("carriers") || flag("hom_or_hemizygous"),
        ),
        "couple" => (
            with_note(&sec::couple::lines(&d)?, COUPLE_NOTE),
            flag("shared") || flag("x_linked"),
        ),
        "roh" => {
            let lines = sec::roh::lines(&x)?;
            let f = lines
                .iter()
                .find(|l| l.starts_with("**"))
                .cloned()
                .map_or_else(|| first(lines.clone()), Ok)?;
            let text = if f == ROH_NOTE {
                Some(f)
            } else {
                with_note(&[f, ROH_NOTE.to_owned()], ROH_NOTE)
            };
            (text, flag("upd_flag"))
        }
        "pypgx" => (
            sec::pgx::pypgx_lines(&x)?.into_iter().next(),
            flag("low_gq") || flag("not_callable"),
        ),
        "annotsv" => (Some(first(sec::annotsv::lines(&d)?)?), flag("count")),
        "viral" => (
            with_note(&sec::viral::lines(&x)?, VIRAL_NOTE),
            flag("found"),
        ),
        "archaic" => (Some(first(sec::archaic::lines(data)?)?), false),
        "lpa" => (Some(first(sec::lpa::lines(&x)?)?), false),
        "cpsr" => {
            let empty = Value::Object(IndexMap::new());
            let dbs = sm.get("databases").filter(|v| truthy(v)).unwrap_or(&empty);
            let lines = sec::cpsr::lines(&d, &D::new(dbs, "databases"))?;
            let Some(l0) = lines.into_iter().next() else {
                return Ok(Some((None, false)));
            };
            let plp = get(data, "plp").filter(|v| truthy(v));
            let n = match plp {
                None => 0,
                Some(Value::Array(a)) => a.len(),
                Some(Value::Object(o)) => o.len(),
                Some(Value::Str(t)) => t.chars().count(),
                Some(_) => return Err(Malformed("plp has no length".into())),
            };
            (
                Some(format!("Pathogenic/Likely Pathogenic rows: {n}. {l0}")),
                n > 0,
            )
        }
        "mito" => {
            let lines = sec::mito::mito_lines(&x)?;
            let f = lines
                .iter()
                .find(|l| l.starts_with("ClinVar pathogenic"))
                .or(lines.first())
                .cloned();
            let a = get(data, "annotation").filter(|v| truthy(v));
            let plp = or_zero(a.and_then(|a| get(a, "clinvar_plp")));
            (f, positive(&plp)?)
        }
        _ => return Ok(None),
    }))
}

/// `(text, tone)` of one summary section, tone `info` or `attention`;
/// `None` while it has no file. The text comes from og-render's helpers, so
/// it reads as the report does; no tone says a section is fine (a ClinVar
/// screen with no hit is `info`, never an all-clear).
pub fn headline(
    key: &str,
    sec: Option<&Value>,
    summary: &Summary,
) -> Option<(String, &'static str)> {
    let sec = sec?;
    let state = get(sec, "state");
    let state_s = match state {
        None | Some(Value::Null) => return None,
        Some(Value::Str(x)) if x == "missing" => return None,
        Some(v) => text_of(v),
    };
    let title = get(sec, "title")
        .filter(|v| truthy(v))
        .map_or_else(|| key.to_owned(), text_of);
    if state_s == "failed" || state_s == "unreadable" {
        let note = get(sec, "note")
            .filter(|v| truthy(v))
            .map_or_else(|| state_s.clone(), text_of);
        return Some((format!("{title}: {note}"), "attention"));
    }
    if state_s == "stale" {
        let note = og_render::notes::stale_note(&summary.section(key)).unwrap_or_default();
        return Some((format!("{title}: {note}"), "attention"));
    }
    let (text, attention) = line_of(key, summary)
        .ok()
        .flatten()
        .unwrap_or((None, false));
    let text = text
        .unwrap_or_else(|| format!("{title}: finished, see the report"))
        .replace("**", "");
    Some((text, if attention { "attention" } else { "info" }))
}

/// The report's anchor of a summary key: the id of the section that lists
/// it among its keys.
pub fn anchor_of(key: &str) -> Option<&'static str> {
    og_render::chapters::CHAPTERS
        .iter()
        .flat_map(|c| c.sections.iter())
        .find(|s| s.keys.contains(&key))
        .map(|s| s.id)
}

// --- the summary the headlines read --------------------------------------------------------

/// `YYYY-MM-DDTHH:MM:SS` (the first 19 characters) as UTC epoch seconds.
fn utc_seconds(text: &str) -> Option<f64> {
    let t: String = text.chars().take(19).collect();
    let (date, time) = t.split_once('T')?;
    let mut d = date.split('-');
    let mut h = time.split(':');
    let num = |x: Option<&str>, max: usize| -> Option<i64> {
        let x = x?;
        (!x.is_empty() && x.len() <= max && x.bytes().all(|b| b.is_ascii_digit()))
            .then(|| x.parse().ok())
            .flatten()
    };
    let (y, mo, dd) = (num(d.next(), 4)?, num(d.next(), 2)?, num(d.next(), 2)?);
    let (hh, mi, ss) = (num(h.next(), 2)?, num(h.next(), 2)?, num(h.next(), 2)?);
    if d.next().is_some() || h.next().is_some() || y < 1 {
        return None;
    }
    let dt = jiff::civil::DateTime::new(
        i16::try_from(y).ok()?,
        i8::try_from(mo).ok()?,
        i8::try_from(dd).ok()?,
        i8::try_from(hh).ok()?,
        i8::try_from(mi).ok()?,
        i8::try_from(ss).ok()?,
        0,
    )
    .ok()?;
    let z = dt.to_zoned(TimeZone::UTC).ok()?;
    Some(z.timestamp().as_second() as f64)
}

/// The done steps of a state's JSON, by number.
fn done_steps(st: &Value) -> Vec<&Value> {
    match get(st, "steps") {
        Some(Value::Array(a)) => a
            .iter()
            .filter(|s| get_str(s, "state") == Some("done"))
            .collect(),
        _ => Vec::new(),
    }
}

/// `<run>/summary.json` when it is newer than every done step's end, else
/// `None`. A summary.json that links out of the run folder is never read.
pub fn summary_file(d: &Path, st: &Value) -> Option<Value> {
    let p = d.join("summary.json");
    if !inside(&p, &realpath(d)) {
        return None;
    }
    let mtime = pycompat::mtime::mtime(&p).ok()?;
    for s in done_steps(st) {
        let e = get_str(s, "ended_utc").filter(|e| !e.is_empty())?;
        if utc_seconds(e)? > mtime {
            return None;
        }
    }
    let v = json::loads(&String::from_utf8_lossy(&fs::read(&p).ok()?)).ok()?;
    matches!(get(&v, "sections"), Some(Value::Object(_))).then_some(v)
}

/// The cache key of a run's summary: the done steps, the message count and
/// summary.json's stamp.
#[derive(Debug, Clone, PartialEq)]
pub struct SummaryKey {
    pub done: Vec<String>,
    pub events: Option<Value>,
    pub stamp: Option<(i128, u64)>,
}

struct SummaryEntry {
    key: SummaryKey,
    headlines: Option<Value>,
    summary: Option<Arc<Value>>,
}

/// A summary.json's stamp and the overview read from it.
type Stamped = ((i128, u64), Value);

/// The overview of a summary for the home page (`levels`, `attention`,
/// `top`, `lead`, ...), and the results built from one: what
/// [`crate::plugins::ResultsFns`] gives when the site has them.
pub type Overview = Option<Value>;

// --- the shared state ---------------------------------------------------------------------

/// The state both the requests and the event streams share.
pub struct Live {
    pub cfg: Config,
    /// The root, links followed.
    pub root: PathBuf,
    pub store: Store,
    pub steps: Vec<Step>,
    pub tz: TimeZone,
    pub access: AccessLog,
    /// The Host names a request may carry.
    pub allowed: BTreeSet<String>,
    /// Event streams open now, and how many may be.
    pub streams: AtomicUsize,
    pub max_streams: AtomicUsize,
    pub stopping: AtomicBool,
    summaries: Mutex<HashMap<String, SummaryEntry>>,
    /// summary.json path -> (its stamp, its overview).
    attention_cache: Mutex<HashMap<PathBuf, Stamped>>,
    /// run -> (the summary results.json reads now, its overview).
    summary_overviews: Mutex<HashMap<String, (Arc<Value>, Overview)>>,
    /// run -> (the summary key, the overview of that summary).
    live_overview: Mutex<HashMap<String, (SummaryKey, Overview)>>,
}

impl Live {
    pub fn new(cfg: Config, steps: Vec<Step>, access: AccessLog) -> Live {
        let root = realpath(&cfg.root);
        let mut allowed: BTreeSet<String> = ["127.0.0.1", "localhost", "::1"]
            .iter()
            .map(|h| (*h).to_owned())
            .collect();
        if cfg.bind != "0.0.0.0" && cfg.bind != "::" {
            allowed.insert(
                cfg.bind
                    .trim_matches(|c| c == '[' || c == ']')
                    .to_lowercase(),
            );
        }
        Live {
            store: Store::new(&cfg.state_dir),
            max_streams: AtomicUsize::new(cfg.max_streams),
            root,
            steps,
            tz: TimeZone::system(),
            access,
            allowed,
            streams: AtomicUsize::new(0),
            stopping: AtomicBool::new(false),
            summaries: Mutex::new(HashMap::new()),
            attention_cache: Mutex::new(HashMap::new()),
            summary_overviews: Mutex::new(HashMap::new()),
            live_overview: Mutex::new(HashMap::new()),
            cfg,
        }
    }

    /// Every run the site knows now.
    pub fn runs(&self) -> BTreeMap<String, Run> {
        runs::list(&self.root, &self.store.known())
    }

    /// Run `name`, when it is a valid name of a run the site knows.
    pub fn run(&self, name: &str) -> Option<Run> {
        if !runs::valid_run(name) {
            return None;
        }
        self.runs().remove(name)
    }

    /// The state of `run`'s newest run, without headlines or plugins (the
    /// home page's view), and its messages.
    pub fn base_state(&self, run: &Run) -> (State<'_>, Arc<Vec<Value>>) {
        let d = run.folder(&self.root);
        let events = self.store.events(&run.name);
        let st = self.base_state_of(&run.name, &d, &events);
        (st, events)
    }

    /// [`Live::base_state`] of the messages `events`.
    pub fn base_state_of(&self, name: &str, d: &Path, events: &[Value]) -> State<'_> {
        self.base_state_of_at(name, d, events, None)
    }

    fn base_state_of_at(
        &self,
        name: &str,
        d: &Path,
        events: &[Value],
        replay: Option<f64>,
    ) -> State<'_> {
        let cur = current_run(events);
        let status = files::read_run_status(&d.join("logs").join("run_status.tsv"))
            .ok()
            .flatten();
        let started = status.as_ref().and_then(|s| s.started());
        let pt_path = d.join("nextflow").join("process_status.tsv");
        let mut pt = files::read_process_status(&pt_path).ok().flatten();
        if pt.is_some()
            && let Some(st) = started
            && pycompat::mtime::mtime(&pt_path).is_ok_and(|m| m < st)
        {
            pt = None; // an earlier run's table
        }
        let events_path = self.store.path(name);
        let inp = Inputs {
            status: status.as_ref(),
            process_status: pt.as_ref(),
            trace_tasks: tasks::trace_tasks(&trace_files(d, started), name, &self.tz),
            event_tasks: tasks::event_tasks(cur, name),
            events: cur,
            sources: files::summary_sources(d),
            heard: match replay {
                Some(at) => Some(
                    cur.last()
                        .and_then(|ev| get_str(ev, "utcTime"))
                        .and_then(utc_seconds)
                        .unwrap_or(at),
                ),
                None => og_runlog::state::last_heard(d, &[events_path.as_path()]),
            },
        };
        let mut st = og_runlog::state::build(
            name,
            &self.steps,
            &inp,
            replay.unwrap_or_else(now),
            &self.tz,
        );
        if run_finished(events) && matches!(st.run, RunState::Running | RunState::Stopped) {
            st.run = RunState::Finished;
        }
        st
    }

    /// [`Live::base_state`] as JSON, with `events`: the message count the
    /// page paints a message ahead of `state.json` past.
    pub fn base_json(&self, st: &State<'_>, events: usize) -> Value {
        let mut v = st.to_json();
        if let Value::Object(o) = &mut v {
            o.insert("events".into(), Value::Int(events as i128));
        }
        v
    }

    /// `state.json` of `run`: the base state, each section's report anchor,
    /// the headlines, the summary's stamp and what the plugins add.
    pub fn state(&self, run: &Run) -> Value {
        let d = run.folder(&self.root);
        let (st, events) = self.base_state(run);
        let plugins = crate::plugins::state_hooks(self, run, &d, &st, &Value::Null, &events);
        self.state_json(run, &d, &st, events.len(), plugins)
    }

    /// The same page state at a replayed time, with the demo's invented host.
    pub(crate) fn replay_state(
        &self,
        run: &Run,
        at: f64,
        host: &dyn og_runlog::host::Host,
    ) -> Value {
        let d = run.folder(&self.root);
        let events = self.store.events(&run.name);
        let st = self.base_state_of_at(&run.name, &d, &events, Some(at));
        let ctx = crate::ops::Ctx {
            now: at,
            steps: &self.steps,
            host,
            tz: &self.tz,
        };
        let plugins = crate::plugins::state_hooks_with(self, &d, || {
            crate::ops::hook(Some(&d), &st, &events, &self.cfg.ops, &ctx).to_json()
        });
        self.state_json(run, &d, &st, events.len(), plugins)
    }

    fn state_json(
        &self,
        run: &Run,
        d: &Path,
        st: &State<'_>,
        events: usize,
        plugins: Value,
    ) -> Value {
        let mut v = self.base_json(st, events);
        let mut anchors = Obj::new();
        for k in st.steps.iter().flat_map(|s| s.step.sections.iter()) {
            if let Some(a) = anchor_of(k) {
                anchors.insert(k.clone(), s(a));
            }
        }
        let headlines = self.headlines(&run.name, d, &v);
        let stamp = file_stamp(&d.join("summary.json"))
            .map_or(Value::Null, |(m, n)| Value::Str(format!("{m}-{n}")));
        if let Value::Object(o) = &mut v {
            o.insert("anchors".into(), Value::Object(anchors));
            o.insert("headlines".into(), headlines);
            o.insert("results_stamp".into(), stamp);
            o.insert("plugins".into(), plugins);
        }
        v
    }

    fn summary_key(&self, d: Option<&Path>, st: &Value) -> SummaryKey {
        SummaryKey {
            done: done_steps(st)
                .iter()
                .filter_map(|s| get_str(s, "step").map(str::to_owned))
                .collect(),
            events: get(st, "events").cloned(),
            stamp: d.and_then(|d| file_stamp(&d.join("summary.json"))),
        }
    }

    /// The run's summary for the steps that are done: its summary.json when
    /// that is newer than the end of every done step, else collected now.
    /// Read again only when the done steps, the message count or
    /// summary.json change; `None` while nothing is done.
    pub fn summary(&self, sample: &str, d: &Path, st: &Value) -> Option<Arc<Value>> {
        let key = self.summary_key(Some(d), st);
        if let Some(e) = lock(&self.summaries).get(sample)
            && e.key == key
        {
            return e.summary.clone();
        }
        let summary = if !key.done.is_empty() && d.is_dir() {
            summary_file(d, st)
                .or_else(|| self.collect(sample, d))
                .map(Arc::new)
        } else {
            None
        };
        lock(&self.summaries).insert(
            sample.to_owned(),
            SummaryEntry {
                key,
                headlines: None,
                summary: summary.clone(),
            },
        );
        summary
    }

    /// `og summary collect` of the folder, in this process; `None` (and a
    /// line on stderr) when it cannot be read.
    fn collect(&self, sample: &str, d: &Path) -> Option<Value> {
        use og_summary::collect::{Group, Input, group_c, post};
        let input = Input {
            sample: sample.to_owned(),
            sample_dir: d.to_string_lossy().into_owned(),
            declared_sex: None,
            assets: self.cfg.repo_assets.to_string_lossy().into_owned(),
            group: Group::All,
            tables: bin_tables(group_c::util::default_tables(), &self.cfg.repo_assets),
        };
        let got = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| post::collect(&input)));
        match got {
            Ok(Ok(v)) => Some(v),
            Ok(Err(_)) | Err(_) => {
                eprintln!("og live serve: summary of {sample} not read");
                None
            }
        }
    }

    /// `{section: {text, tone, step}}` for the sections of the done steps.
    pub fn headlines(&self, sample: &str, d: &Path, st: &Value) -> Value {
        let summary = self.summary(sample, d, st);
        let key = self.summary_key(Some(d), st);
        if let Some(e) = lock(&self.summaries).get(sample)
            && e.key == key
            && let Some(h) = &e.headlines
        {
            return h.clone();
        }
        let mut out = Obj::new();
        let parsed = summary
            .as_ref()
            .and_then(|sm| Summary::parse(&dumps(sm)).ok());
        if let (Some(sm), Some(parsed)) = (&summary, &parsed) {
            let sections = get(sm, "sections");
            for step in done_steps(st) {
                let Some(Value::Array(keys)) = get(step, "sections") else {
                    continue;
                };
                for k in keys.iter().filter_map(|k| match k {
                    Value::Str(x) => Some(x.as_str()),
                    _ => None,
                }) {
                    let sec = sections.and_then(|x| get(x, k));
                    if let Some((text, tone)) = headline(k, sec, parsed) {
                        out.insert(
                            k.to_owned(),
                            obj([
                                ("text", Value::Str(text)),
                                ("tone", s(tone)),
                                ("step", get(step, "step").cloned().unwrap_or(Value::Null)),
                            ]),
                        );
                    }
                }
            }
        }
        let out = Value::Object(out);
        if let Some(e) = lock(&self.summaries).get_mut(sample)
            && e.key == key
            && e.summary.as_ref().map(Arc::as_ptr) == summary.as_ref().map(Arc::as_ptr)
        {
            e.headlines = Some(out.clone());
        }
        out
    }

    /// The overview of `<run>/summary.json` and its time, cached by its
    /// stamp; `(None, None)` without one. A summary.json that links out of
    /// its run is never read. Without the results the overview is `None`.
    pub fn overview(&self, d: &Path) -> (Overview, Option<f64>) {
        let p = d.join("summary.json");
        if !inside(&p, &realpath(d)) || !inside(&p, &self.root) {
            return (None, None);
        }
        let Ok(mtime) = pycompat::mtime::mtime(&p) else {
            return (None, None);
        };
        let Some(fns) = &self.cfg.results else {
            return (None, Some(mtime));
        };
        let stamp = file_stamp(&p);
        if let Some(st) = stamp
            && let Some((h, v)) = lock(&self.attention_cache).get(&p)
            && *h == st
        {
            return (Some(v.clone()), Some(mtime));
        }
        let n = fs::read(&p)
            .ok()
            .and_then(|b| json::loads(&String::from_utf8_lossy(&b)).ok())
            .ok_or_else(|| "not JSON".to_owned())
            .and_then(|v| (fns.overview)(&v));
        match n {
            Ok(v) => {
                if let Some(st) = stamp {
                    lock(&self.attention_cache).insert(p, (st, v.clone()));
                }
                (Some(v), Some(mtime))
            }
            Err(_) => {
                let run = d.file_name().map(|n| n.to_string_lossy().into_owned());
                eprintln!(
                    "og live serve: {}/summary.json not read",
                    run.unwrap_or_default()
                );
                (None, Some(mtime))
            }
        }
    }

    /// The count of cards worth a look, from `<run>/summary.json` alone.
    pub fn attention(&self, d: &Path) -> (Option<Value>, Option<f64>) {
        let (o, mtime) = self.overview(d);
        (o.and_then(|o| get(&o, "attention").cloned()), mtime)
    }

    /// Whether `name`'s pages treat it as a child, over the families the
    /// site lists now (or those of `runs`).
    pub fn child(&self, name: &str, runs: Option<&BTreeMap<String, Run>>) -> bool {
        let fams = match runs {
            Some(r) => runs::families_of(r),
            None => runs::families_of(&self.runs()),
        };
        runs::child_of(&fams, name)
    }

    /// The overview of the summary results.json reads now, cached with that
    /// summary; `None` while there is none.
    pub fn summary_overview(&self, name: &str, d: &Path, st: &Value) -> Overview {
        let sm = self.summary(name, d, st)?;
        let fns = self.cfg.results.as_ref()?;
        if let Some((s0, o)) = lock(&self.summary_overviews).get(name)
            && Arc::ptr_eq(s0, &sm)
        {
            return o.clone();
        }
        let o = match (fns.overview)(&sm) {
            Ok(v) => Some(v),
            Err(_) => {
                eprintln!("og live serve: overview of {name} not built");
                None
            }
        };
        lock(&self.summary_overviews).insert(name.to_owned(), (sm, o.clone()));
        o
    }

    /// Whether `name`'s own page keeps its report closed: a child whose
    /// results hold a result to act on, counted from the summary its
    /// results.json reads (summary.json while that has none).
    pub fn gated(
        &self,
        name: &str,
        d: &Path,
        st: &Value,
        runs: Option<&BTreeMap<String, Run>>,
    ) -> bool {
        if !self.child(name, runs) {
            return false;
        }
        let o = self
            .summary_overview(name, d, st)
            .or_else(|| self.overview(d).0);
        has_act(o.as_ref())
    }

    /// `/runs.json`: each run's progress, and the families with their
    /// members, for the home page. No result beyond the attention count,
    /// the counts by level, the titles of the first findings and the lead;
    /// no path.
    pub fn runs_index(&self) -> Value {
        let runs = self.runs();
        lock(&self.summaries).retain(|k, _| runs.contains_key(k));
        lock(&self.live_overview).retain(|k, _| runs.contains_key(k));
        lock(&self.summary_overviews).retain(|k, _| runs.contains_key(k));
        let keep: BTreeSet<PathBuf> = runs
            .values()
            .map(|r| r.folder(&self.root).join("summary.json"))
            .collect();
        lock(&self.attention_cache).retain(|k, _| keep.contains(k));
        let mut people = Vec::new();
        let mut fams = Vec::new();
        for (name, r) in &runs {
            let url = format!("/runs/{name}/");
            if r.is_family() {
                fams.push(obj([
                    ("name", s(name)),
                    ("url", Value::Str(url)),
                    ("members", Value::Array(r.members.clone())),
                    ("unreadable", Value::Bool(r.unreadable)),
                ]));
                continue;
            }
            let d = r.folder(&self.root);
            let (st, events) = self.base_state(r);
            let v = self.base_json(&st, events.len());
            let planned: Vec<_> = st
                .steps
                .iter()
                .filter(|x| x.state.as_str() != "not planned")
                .collect();
            let count = |w: &str| planned.iter().filter(|x| x.state.as_str() == w).count() as i128;
            let run_v = get(&v, "run").cloned().unwrap_or(Value::Null);
            let field = |k: &str| get(&run_v, k).cloned().unwrap_or(Value::Null);
            let mut item = Obj::new();
            item.insert("name".into(), s(name));
            item.insert("kind".into(), s(&r.kind));
            item.insert("url".into(), Value::Str(url));
            item.insert("state".into(), field("state"));
            item.insert("started_utc".into(), field("started_utc"));
            item.insert("last_heard_utc".into(), field("last_heard_utc"));
            item.insert("finished_utc".into(), field("finished_utc"));
            item.insert("declared_sex".into(), field("declared_sex"));
            item.insert("planned".into(), Value::Int(planned.len() as i128));
            item.insert("done".into(), Value::Int(count("done")));
            item.insert("running".into(), Value::Int(count("running")));
            item.insert("failed".into(), Value::Int(count("failed")));
            let (mut o, mtime) = self.overview(&d);
            // While a run goes, the summary its page last read (results.json's
            // own) is newer than summary.json: count from it while it is
            // current, so the home never says less than the person page.
            if st.run == RunState::Running && self.cfg.results.is_some() {
                let key = self.summary_key(Some(&d), &v);
                let hit = lock(&self.summaries)
                    .get(name)
                    .filter(|e| e.key.done == key.done && e.key.stamp == key.stamp)
                    .map(|e| (e.key.clone(), e.summary.clone()));
                if let Some((hk, Some(sm))) = hit {
                    let cached = lock(&self.live_overview)
                        .get(name)
                        .filter(|(k, _)| *k == hk)
                        .map(|(_, o)| o.clone());
                    let lo = cached.unwrap_or_else(|| {
                        let fns = self.cfg.results.as_ref().map(|f| f.overview);
                        let lo = fns.and_then(|f| match f(&sm) {
                            Ok(v) => Some(v),
                            Err(_) => {
                                eprintln!("og live serve: overview of {name} not built");
                                None
                            }
                        });
                        lock(&self.live_overview).insert(name.clone(), (hk, lo.clone()));
                        lo
                    });
                    if lo.is_some() {
                        o = lo;
                    }
                }
            }
            // A child's results to act on decide whether its own page is
            // closed: when the summary its results.json reads, or the
            // overview above, has one, that summary's overview replaces it,
            // so the home page closes the child exactly when its page does.
            if self.child(name, Some(&runs)) {
                let so = self.summary_overview(name, &d, &v);
                if so.is_some() && (has_act(so.as_ref()) || has_act(o.as_ref())) {
                    o = so;
                }
            }
            item.insert(
                "attention".into(),
                o.as_ref()
                    .and_then(|x| get(x, "attention").cloned())
                    .unwrap_or(Value::Null),
            );
            item.insert("results".into(), Value::Bool(mtime.is_some()));
            if let Some(o) = &o {
                for k in ["levels", "top", "lead", "lead_closed"] {
                    item.insert(k.into(), get(o, k).cloned().unwrap_or(Value::Null));
                }
                item.insert(
                    "screened".into(),
                    get(o, "screened").cloned().unwrap_or(Value::Bool(false)),
                );
            }
            people.push(Value::Object(item));
        }
        obj([
            ("runs", Value::Array(people)),
            ("families", Value::Array(fams)),
        ])
    }
}

/// Whether an overview counts a result to act on.
pub fn has_act(o: Option<&Value>) -> bool {
    o.and_then(|o| get(o, "levels"))
        .and_then(|l| get(l, "act"))
        .is_some_and(truthy)
}

/// The folder of the `bin/` tables the summary reads: og-summary's own
/// choice (`OG_BIN_TABLES`, else the checkout the binary sits in), else the
/// `bin/` beside the site's `assets/`. An og built outside a checkout would
/// otherwise read `bin/` under the current folder and fail the repeat
/// expansion headline.
fn bin_tables(default: String, repo_assets: &Path) -> String {
    if Path::new(&default).join("hla_alleles.tsv").is_file() {
        return default;
    }
    match repo_assets.parent().map(|r| r.join("bin")) {
        Some(b) if b.join("hla_alleles.tsv").is_file() => b.to_string_lossy().into_owned(),
        _ => default,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_bin_tables_fall_back_to_the_checkout_of_the_assets() {
        let repo = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../..");
        let got = bin_tables("bin-nowhere".into(), &repo.join("assets"));
        assert!(
            Path::new(&got)
                .join("str_intermediate_ranges.tsv")
                .is_file(),
            "{got}"
        );
        let mine = repo.join("bin").to_string_lossy().into_owned();
        assert_eq!(bin_tables(mine.clone(), Path::new("/nowhere/assets")), mine);
        assert_eq!(
            bin_tables("bin-nowhere".into(), Path::new("/nowhere/assets")),
            "bin-nowhere"
        );
    }

    fn ev(name: &str) -> Value {
        obj([("event", s(name))])
    }

    #[test]
    fn the_newest_run() {
        let done = [ev("started"), ev("pgp_run_finished")];
        assert!(run_finished(&done));
        let rerun = [ev("started"), ev("pgp_run_finished"), ev("pgp_step")];
        assert!(!run_finished(&rerun) && !run_finished(&[]));
        assert_eq!(current_run(&rerun), &rerun[2..]);
        let bare = [
            ev("started"),
            ev("completed"),
            ev("started"),
            ev("process_started"),
        ];
        assert_eq!(current_run(&bare), &bare[2..]);
        assert_eq!(current_run(&done), &done[..]);
    }

    #[test]
    fn task_names() {
        let t = |n: &str| task_name(n);
        assert_eq!(
            t("SV:MANTA (s1)"),
            (Some("MANTA".into()), Some("s1".into()))
        );
        assert_eq!(
            t("  A:B:CLINVAR_SCREEN  "),
            (Some("CLINVAR_SCREEN".into()), None)
        );
        assert_eq!(t("ROH(x y)"), (Some("ROH".into()), Some("x y".into())));
        assert_eq!(t("a b"), (None, None));
        assert_eq!(t("A:"), (Some(String::new()), None));
        assert_eq!(t(""), (None, None));
        assert_eq!(t("P (a) (b)"), (Some("P".into()), Some("a) (b".into())));
    }

    #[test]
    fn utc_times() {
        assert_eq!(utc_seconds("1970-01-01T00:01:00Z"), Some(60.0));
        assert_eq!(
            utc_seconds("2026-01-02T10:07:00.123"),
            utc_seconds("2026-01-02T10:07:00")
        );
        assert_eq!(utc_seconds("2026-01-02 10:07:00"), None);
        assert_eq!(utc_seconds(""), None);
    }

    /// The carrier headline of a summary whose carrier section is `data`.
    fn carrier_headline(data: &str) -> (String, &'static str) {
        let sec = format!(r#"{{"state": "ok", "title": "Carrier list", "data": {data}}}"#);
        let summary = og_render::model::Summary::parse(&format!(
            r#"{{"databases": {{}}, "sections": {{"carrier": {sec}}}}}"#
        ))
        .unwrap();
        let raw = pycompat::json::loads(&sec).unwrap();
        headline("carrier", Some(&raw), &summary).unwrap()
    }

    #[test]
    fn a_carrier_list_of_only_hom_or_hemizygous_rows_asks_for_attention() {
        // A man with only a hemizygous GLA row: no carrier, but no unchanged
        // copy either, so the headline asks for attention as the Python
        // server's does (bin/owngenome_live.py, _h_carrier).
        let gla = r#"{"gene": "GLA", "condition": "Fabry disease", "inheritance": "XL",
            "variant": "chrX:101398000 C>T", "genotype": "hemizygous", "zygosity": "hemizygous",
            "significance": "Pathogenic", "gq": 50, "dp": 30, "low_confidence": false}"#;
        let only = carrier_headline(&format!(
            r#"{{"genes_on_list": 3, "carriers": [], "hom_or_hemizygous": [{gla}], "set_aside": []}}"#
        ));
        assert_eq!(only.1, "attention", "{}", only.0);
        assert!(only.0.contains("1 with no unchanged copy"), "{}", only.0);
        let none = carrier_headline(
            r#"{"genes_on_list": 3, "carriers": [], "hom_or_hemizygous": [], "set_aside": []}"#,
        );
        assert_eq!(none.1, "info", "{}", none.0);
        // A summary written before the list existed reads as before.
        let old = carrier_headline(r#"{"genes_on_list": 3, "carriers": [], "set_aside": []}"#);
        assert_eq!(old.1, "info", "{}", old.0);
    }

    #[test]
    fn a_note_joined_as_two_sentences() {
        let l = |a: &[&str]| a.iter().map(|x| (*x).to_owned()).collect::<Vec<_>>();
        assert_eq!(with_note(&l(&["A", "N."]), "N.").unwrap(), "A. N.");
        assert_eq!(with_note(&l(&["A.", "N."]), "N.").unwrap(), "A. N.");
        assert_eq!(with_note(&l(&["A", "B"]), "N.").unwrap(), "A");
        assert_eq!(with_note(&[], "N."), None);
    }
}
