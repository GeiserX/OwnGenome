//! The one freshness issue.

use crate::http::Fetch;
use crate::{Lookup, LookupError, array, failed, field, need_str, opt_field};
use indexmap::IndexMap;
use pycompat::json::{Dump, Value, dumps};

pub const LABEL: &str = "freshness";
pub const MARKER: &str = "<!-- freshness-report: edited by .github/workflows/freshness.yml -->";
pub const TITLE: &str = "Freshness report: databases, bundles, images and links";

/// What [`update`] did to the issue.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Action {
    Created,
    Edited,
    ReopenedAndEdited,
}

impl Action {
    pub fn as_str(self) -> &'static str {
        match self {
            Action::Created => "created",
            Action::Edited => "edited",
            Action::ReopenedAndEdited => "reopened and edited",
        }
    }
}

fn s(x: &str) -> Value {
    Value::Str(x.into())
}

fn json(v: &Value) -> Vec<u8> {
    dumps(v, &Dump::default()).expect("no floats").into_bytes()
}

fn call(
    f: &dyn Fetch,
    path: &str,
    method: &str,
    data: Option<&[u8]>,
    retries: Option<u32>,
) -> Lookup<Value> {
    let gh = f.github_api();
    let url = gh.url(path);
    let mut req = gh.request(&url, method, data);
    req.retries = retries;
    f.json_with(req)
}

/// A body cut to fit an issue (65,536 characters), with the marker first.
pub fn issue_body(body: &str) -> String {
    let mut b = if body.contains(MARKER) {
        body.to_owned()
    } else {
        format!("{MARKER}\n{body}")
    };
    if b.chars().count() > 65_000 {
        b = b.chars().take(64_000).collect::<String>()
            + "\n\n(Cut at 64,000 characters; the full report is in the run log.)\n";
    }
    b
}

fn number(i: &Value) -> Lookup<i128> {
    match field(i, "number")? {
        Value::Int(n) => Ok(*n),
        _ => Err(LookupError::Shape("an issue number is not a number".into())),
    }
}

/// Edits the one open issue that carries [`MARKER`]; reopens the newest
/// closed one when none is open; creates it when none exists. An issue
/// with the label but without the marker (another workflow's) is never
/// touched. Returns what it did and the issue's number.
pub fn update(f: &dyn Fetch, repo: &str, body: &str) -> Lookup<(Action, i128)> {
    if f.github_api().token.is_none() {
        return failed("--update-issue needs GITHUB_TOKEN");
    }
    let body = issue_body(body);
    if let Err(e) = call(
        f,
        &format!("/repos/{repo}/labels/{LABEL}"),
        "GET",
        None,
        Some(1),
    ) {
        if !e.to_string().contains("HTTP 404") {
            return Err(e);
        }
        let label = Value::Object(IndexMap::from([
            ("name".into(), s(LABEL)),
            ("color".into(), s("c5def5")),
            (
                "description".into(),
                s("Monthly report of pins and databases behind upstream"),
            ),
        ]));
        call(
            f,
            &format!("/repos/{repo}/labels"),
            "POST",
            Some(&json(&label)),
            None,
        )?;
    }
    let mut found: Vec<(i128, bool)> = vec![];
    for page in 1.. {
        let batch = call(
            f,
            &format!("/repos/{repo}/issues?labels={LABEL}&state=all&per_page=100&page={page}"),
            "GET",
            None,
            None,
        )?;
        let batch = array(&batch, "issue list")
            .map_err(|_| LookupError::Failed("issue list is not a list".into()))?;
        for i in batch {
            let is_pr = opt_field(i, "pull_request")?.is_some();
            let mine = match opt_field(i, "body")? {
                Some(Value::Str(b)) => b.contains(MARKER),
                _ => false,
            };
            if !is_pr && mine {
                let open = matches!(field(i, "state")?, Value::Str(st) if st == "open");
                found.push((number(i)?, open));
            }
        }
        if batch.len() < 100 {
            break;
        }
    }
    found.sort();
    let open: Vec<i128> = found.iter().filter(|x| x.1).map(|x| x.0).collect();
    let closed: Vec<i128> = found.iter().filter(|x| !x.1).map(|x| x.0).collect();
    let mut payload = IndexMap::from([
        ("title".to_owned(), s(TITLE)),
        ("body".to_owned(), s(&body)),
    ]);
    let (target, action) = if let Some(&n) = open.first() {
        if open.len() > 1 {
            println!(
                "WARNING: {} open freshness issues; editing #{n} only.",
                open.len()
            );
        }
        (Some(n), Action::Edited)
    } else if let Some(&n) = closed.last() {
        payload.insert("state".into(), s("open"));
        (Some(n), Action::ReopenedAndEdited)
    } else {
        (None, Action::Created)
    };
    let res = match target {
        None => {
            payload.insert("labels".into(), Value::Array(vec![s(LABEL)]));
            // No retry: a create whose answer was lost may have landed, and
            // a second POST would open a second issue. A failed create fails
            // the run.
            let data = json(&Value::Object(payload));
            call(
                f,
                &format!("/repos/{repo}/issues"),
                "POST",
                Some(&data),
                Some(0),
            )?
        }
        Some(n) => {
            let data = json(&Value::Object(payload));
            call(
                f,
                &format!("/repos/{repo}/issues/{n}"),
                "PATCH",
                Some(&data),
                None,
            )?
        }
    };
    println!(
        "Issue {}: {}",
        action.as_str(),
        need_str(&res, "html_url", "issue URL")?
    );
    Ok((action, number(&res)?))
}
