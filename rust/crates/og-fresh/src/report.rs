//! The report: every lookup run (in parallel), then the Markdown.

use crate::lookups::{Ctx, DATABASES, ImageState, clinvar, digest, image};
use crate::{DbRow, LookupError};
use std::sync::atomic::{AtomicUsize, Ordering};

pub const SECTIONS: &[&str] = &["databases", "clinvar", "images", "digests"];

/// At most this many lookups at once; requests to one host still go one at
/// a time (`http::Client`).
pub const PARALLEL: usize = 8;

/// What the lookups found, in report order.
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct Report {
    pub db_rows: Vec<DbRow>,
    /// (variable, pinned, newest, kind, note)
    pub image_rows: Vec<(String, String, String, &'static str, String)>,
    pub current: Vec<String>,
    /// Only a new build of the pinned version.
    pub rebuilds: Vec<String>,
    /// Tags that float on a release line: (variable, pinned, newest line).
    pub line_rows: Vec<(String, String, String)>,
    /// (variable, pinned, what the tag points at now, state)
    pub digest_rows: Vec<(String, String, String, String)>,
    /// (lookup, message, the upstream did not answer): every failed
    /// lookup, in report order.
    pub errors: Vec<(String, String, bool)>,
}

impl Report {
    /// A failed lookup.
    pub fn error(&mut self, what: &str, e: &LookupError) {
        self.errors
            .push((what.into(), e.to_string(), e.is_unreachable()));
    }

    /// The failed lookups that count against the run: all of them, or with
    /// `unreachable_warns` those whose upstream did answer.
    pub fn hard_errors(&self, unreachable_warns: bool) -> usize {
        self.errors
            .iter()
            .filter(|e| !(unreachable_warns && e.2))
            .count()
    }

    /// A failed database lookup: an error, and an ERROR row in its place.
    pub fn db_error(&mut self, what: &str, e: &LookupError) {
        self.error(what, e);
        self.db_rows.push(DbRow {
            item: what.into(),
            pinned: String::new(),
            upstream: String::new(),
            state: "ERROR".into(),
            todo: "lookup failed, see Errors".into(),
        });
    }

    /// Adds one image's answer.
    pub fn add_image(
        &mut self,
        var: &str,
        note: &str,
        answer: Result<(String, ImageState), LookupError>,
    ) {
        match answer {
            Err(e) => self.error(var, &e),
            Ok((_, ImageState::Current)) => self.current.push(var.into()),
            Ok((_, ImageState::Rebuild)) => self.rebuilds.push(var.into()),
            Ok((tag, ImageState::Line(n))) => self.line_rows.push((var.into(), tag, n)),
            Ok((tag, ImageState::Behind(rows))) => {
                for (n, kind) in rows {
                    self.image_rows
                        .push((var.into(), tag.clone(), n, kind, note.into()));
                }
            }
        }
    }
}

/// `f(i)` for every `i < n` on up to [`PARALLEL`] threads, results in order.
pub fn parallel<T: Send>(n: usize, f: impl Fn(usize) -> T + Sync) -> Vec<T> {
    let next = AtomicUsize::new(0);
    let mut out: Vec<Option<T>> = (0..n).map(|_| None).collect();
    let slots: Vec<std::sync::Mutex<&mut Option<T>>> =
        out.iter_mut().map(std::sync::Mutex::new).collect();
    std::thread::scope(|s| {
        for _ in 0..PARALLEL.min(n) {
            s.spawn(|| {
                loop {
                    let i = next.fetch_add(1, Ordering::Relaxed);
                    if i >= n {
                        break;
                    }
                    let v = f(i);
                    **slots[i].lock().expect("slot") = Some(v);
                }
            });
        }
    });
    drop(slots);
    out.into_iter()
        .map(|v| v.expect("every index ran"))
        .collect()
}

/// Runs the lookups of `sections`.
pub fn build(c: &Ctx<'_>, sections: &[&str]) -> Report {
    let mut r = Report::default();
    let has = |s: &str| sections.contains(&s);
    let mut dbs: Vec<(&str, crate::lookups::LookupFn)> = vec![];
    if has("databases") || has("clinvar") {
        dbs.push(("ClinVar", clinvar));
    }
    if has("databases") {
        dbs.extend_from_slice(DATABASES);
    }
    let images: Vec<(&String, &String)> = c
        .pins
        .values
        .iter()
        .filter(|(k, _)| k.ends_with("_IMAGE"))
        .collect::<std::collections::BTreeMap<_, _>>()
        .into_iter()
        .collect();
    enum Job<'a> {
        Db(&'a str, crate::lookups::LookupFn),
        Image(&'a str, &'a str),
        Digest(&'a str, &'a str),
    }
    let mut jobs: Vec<Job<'_>> = dbs.iter().map(|(w, f)| Job::Db(w, *f)).collect();
    if has("images") {
        jobs.extend(images.iter().map(|(k, v)| Job::Image(k, v)));
    }
    if has("digests") {
        jobs.extend(images.iter().map(|(k, v)| Job::Digest(k, v)));
    }
    enum Done {
        Db(Result<crate::Found, LookupError>),
        Image(Option<Result<(String, ImageState), LookupError>>),
        Digest(Option<Result<(String, String, String), LookupError>>),
    }
    let done = parallel(jobs.len(), |i| match &jobs[i] {
        Job::Db(_, f) => Done::Db(f(c)),
        Job::Image(_, v) => Done::Image(image(c.registry, v)),
        Job::Digest(_, v) => Done::Digest(digest(c.registry, v)),
    });
    for (job, d) in jobs.iter().zip(done) {
        match (job, d) {
            (Job::Db(..), Done::Db(Ok(f))) => r.db_rows.push(f.row),
            (Job::Db(w, _), Done::Db(Err(e))) => r.db_error(w, &e),
            (Job::Image(var, _), Done::Image(Some(a))) => {
                let note = c.pins.notes.get(*var).map_or("", String::as_str);
                r.add_image(var, note, a);
            }
            (Job::Digest(var, _), Done::Digest(Some(a))) => match a {
                Ok((p, d, s)) => r.digest_rows.push(((*var).into(), p, d, s)),
                Err(e) => r.error(var, &e),
            },
            _ => {}
        }
    }
    r
}

/// A table cell: `|` escaped, line breaks as spaces.
pub fn esc(text: &str) -> String {
    text.replace('|', "\\|").replace('\n', " ")
}

/// The Markdown report of `sections`: every failed lookup under Errors.
pub fn render(r: &Report, sections: &[&str]) -> String {
    render_with(r, sections, false)
}

/// What a database row says when its upstream did not answer and that is
/// a warning.
const NO_ANSWER_STATE: &str = "NO ANSWER";

/// The Markdown report of `sections`. With `unreachable_warns`, a lookup
/// whose upstream did not answer is listed under its own heading as a
/// warning, and only the others are under Errors.
pub fn render_with(r: &Report, sections: &[&str], unreachable_warns: bool) -> String {
    let has = |s: &str| sections.contains(&s);
    let warned = |what: &str| unreachable_warns && r.errors.iter().any(|e| e.2 && e.0 == what);
    let mut out: Vec<String> = vec![];
    if has("databases") || has("clinvar") {
        out.extend(
            [
                "### Databases, bundles and tools Renovate cannot see",
                "",
                "| Item | Pinned | Upstream | State | What to do |",
                "|---|---|---|---|---|",
            ]
            .map(String::from),
        );
        for d in &r.db_rows {
            let mut d = d.clone();
            if d.state == "ERROR" && warned(&d.item) {
                d.state = NO_ANSWER_STATE.into();
                d.todo = "upstream did not answer, see Upstreams that did not answer".into();
            }
            let cells = [&d.item, &d.pinned, &d.upstream, &d.state, &d.todo].map(|c| esc(c));
            out.push(format!("| {} |", cells.join(" | ")));
        }
        out.push(String::new());
    }
    if has("images") {
        let (majors, rest): (Vec<_>, Vec<_>) = r.image_rows.iter().partition(|x| x.3 == "major");
        out.extend(
            [
                "### Images behind on a major version",
                "",
                "Renovate never opens a major update on its own: it waits for approval on its \
                 dependency dashboard, or renovate.json disables the image. Each needs its \
                 revalidation run.",
                "",
            ]
            .map(String::from),
        );
        if majors.is_empty() {
            out.push("None.".into());
        } else {
            out.push("| Variable | Pinned | Newest | Note |".into());
            out.push("|---|---|---|---|".into());
            for (v, p, n, _, note) in majors {
                out.push(format!("| {v} | `{p}` | `{n}` | {} |", esc(note)));
            }
        }
        out.extend(
            [
                "",
                "### Images behind on a minor, patch or build",
                "",
                "Renovate opens most of these as pull requests. The images renovate.json disables, \
                 holds or sends to the dependency dashboard (PCGR, Python) wait there. They are \
                 listed so a stalled or silent Renovate shows.",
                "",
            ]
            .map(String::from),
        );
        if rest.is_empty() {
            out.push("None.".into());
        } else {
            out.push("| Variable | Pinned | Newest | Change | Note |".into());
            out.push("|---|---|---|---|---|".into());
            for (v, p, n, k, note) in rest {
                out.push(format!("| {v} | `{p}` | `{n}` | {k} | {} |", esc(note)));
            }
        }
        out.push(String::new());
        if !r.line_rows.is_empty() {
            let lines: Vec<String> = r
                .line_rows
                .iter()
                .map(|(v, p, n)| format!("{v} `{p}` (newest line `{n}`)"))
                .collect();
            out.push(format!(
                "Floating on a release line (the tag follows its own patch releases; a newer line is \
                 a deliberate move, so it is not counted as behind): {}.",
                lines.join(", ")
            ));
            out.push(String::new());
        }
        if !r.rebuilds.is_empty() {
            out.push(format!(
                "Only a new build of the pinned version (not listed above): {}.",
                r.rebuilds.join(", ")
            ));
            out.push(String::new());
        }
        let cur = if r.current.is_empty() {
            "none".to_owned()
        } else {
            r.current.join(", ")
        };
        out.push(format!("Current: {} images ({cur}).", r.current.len()));
        out.push(String::new());
    }
    if has("digests") {
        out.extend(
            [
                "### Images pinned by digest",
                "",
                "| Variable | Pinned | What the tag points at now (`latest` for a digest-only pin) | State |",
                "|---|---|---|---|",
            ]
            .map(String::from),
        );
        for (v, p, d, s) in &r.digest_rows {
            out.push(format!("| {v} | `{p}` | `{d}` | {} |", esc(s)));
        }
        out.push(String::new());
    }
    let (soft, hard): (Vec<_>, Vec<_>) = r.errors.iter().partition(|e| unreachable_warns && e.2);
    if !hard.is_empty() {
        out.extend(
            [
                "### Errors",
                "",
                "A lookup that fails is an error, never \"current\". This run exits 1.",
                "",
                "| Lookup | Error |",
                "|---|---|",
            ]
            .map(String::from),
        );
        for (w, m, _) in hard {
            out.push(format!("| {} | {} |", esc(w), esc(m)));
        }
        out.push(String::new());
    }
    if !soft.is_empty() {
        out.extend(
            [
                "### Upstreams that did not answer",
                "",
                "Warnings: these lookups got no connection, a timeout, HTTP \
                 5xx or 429 after the retries, or HTTP 403, so nothing was compared. They do \
                 not fail this \
                 run (`--unreachable-warns`); the monthly report counts them as errors.",
                "",
                "| Lookup | Error |",
                "|---|---|",
            ]
            .map(String::from),
        );
        for (w, m, _) in soft {
            out.push(format!("| {} | {} |", esc(w), esc(m)));
        }
        out.push(String::new());
    }
    out.join("\n")
}
