//! `og recheck`: which of a finished sample's results a newer ClinVar or
//! catalog would change, and the replay that redoes only those steps.
//!
//! It reads the sample's `run_manifest.tsv`: the databases and pinned data
//! the run used. For each it asks upstream for the newest release through
//! `og-fresh`'s typed lookups (never the freshness report's text), plus one
//! PGS Catalog lookup per score file, and the newest release of this
//! pipeline. Each item gets one [`Class`]. The result goes to
//! `<sample>/recheck/state.json` (rewritten whole) and one row per item to
//! `<sample>/recheck/log.tsv` (only ever appended). `docs/recheck.md` says
//! how to schedule it and how to read it.

pub mod answers;
pub mod cli;
pub mod files;
pub mod rerun;

use answers::{Answers, Release};
use indexmap::IndexMap;
use og_fresh::registry::image_tag;
use og_fresh::version::ver;
use og_runlog::state::{RunFiles, Sources, StepState};
use og_runlog::steps::Step;
use pycompat::json::Value;
use regex::Regex;
use std::collections::BTreeSet;
use std::path::Path;
use std::sync::LazyLock;

pub const SCHEMA_VERSION: i128 = 1;
pub const DEFAULT_REPO: &str = "GeiserX/OwnGenome";
pub const RAW: &str = "https://raw.githubusercontent.com";
/// The HTML report; the text report is written after it.
pub const REPORT_STEP: &str = "24";
/// ClinVar's `##fileDate` is the release; NCBI's Last-Modified comes a day
/// or two later.
pub const CLINVAR_SLACK_DAYS: i64 = 2;
/// Held exclusive by `--refresh-clinvar`, shared by `scripts/run-all.sh`
/// and `scripts/06-clinvar-screen.sh` (with `flock(1)`).
pub const CLINVAR_LOCK: &str = "clinvar/.refresh.lock";
pub const CLINVAR_FILES: &[&str] = &[
    "clinvar/clinvar_pathogenic_chr.vcf.gz",
    "clinvar/clinvar_chr.vcf.gz",
    "clinvar/clinvar.vcf.gz",
];
/// The looked-up pins: (manifest key, label).
pub const PINNED: &[(&str, &str)] = &[
    ("VEP_CACHE_RELEASE", "Ensembl / VEP cache"),
    ("PCGR_DATA_BUNDLE", "PCGR / CPSR data bundle"),
    ("PYPGX_BUNDLE_VERSION", "pypgx and its bundle"),
    ("PHARMCAT_IMAGE", "PharmCAT"),
    ("hla_database", "IPD-IMGT/HLA"),
];
/// The opt-in steps' pins: the setup.sh option that installs their data and
/// image (the main setup leaves them out). VNtyper's reference data is
/// inside its image, which the run pulls, so it needs no option.
pub const OPT_IN_SETUP: &[(&str, Option<&str>)] = &[
    ("MEGANE_IMAGE", Some("--megane-kmers")),
    ("MOCHA_PANEL_VERSION", Some("--mocha-data")),
    ("MOCHA_IMAGE", Some("--mocha-data")),
    ("LOCITYPER_PANEL_VERSION", Some("--locityper-data")),
    ("LOCITYPER_IMAGE", Some("--locityper-data")),
    ("VNTYPER_IMAGE", None),
];

/// What one item's comparison says.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Class {
    /// What is installed here is newer than what the run used: a rerun
    /// picks it up.
    RerunReady,
    /// ClinVar only: NCBI has a newer file than the one installed.
    RefreshAvailable,
    /// A newer pin is in a pipeline release: update the pipeline.
    InPipelineRelease,
    /// Upstream has a newer release that no pipeline release pins.
    UpstreamOnly,
    /// The lookup answered and nothing is newer.
    NoNewerRelease,
    /// A pin with no upstream lookup that the run's manifest does not name.
    Unknown,
    /// The lookup failed (never read as up to date).
    Error,
}

impl Class {
    pub fn as_str(self) -> &'static str {
        match self {
            Class::RerunReady => "rerun-ready",
            Class::RefreshAvailable => "refresh-available",
            Class::InPipelineRelease => "in-pipeline-release",
            Class::UpstreamOnly => "upstream-only",
            Class::NoNewerRelease => "no-newer-release",
            Class::Unknown => "unknown",
            Class::Error => "error",
        }
    }
}

/// A database file the run read, or a pin.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kind {
    File,
    Pin,
}

/// One item of a check.
#[derive(Debug, Clone, PartialEq)]
pub struct Item {
    pub item: String,
    pub label: String,
    pub kind: Kind,
    pub ran_with: Option<String>,
    pub installed: Option<String>,
    /// Images: the full references (a new digest at the same tag counts).
    pub refs: Option<(Option<String>, Option<String>)>,
    /// The answer, when the lookup gave one (it may name nothing).
    pub upstream: Option<Option<String>>,
    pub error: Option<String>,
    pub class: Class,
    /// The newest pipeline release's tag and its pin, when it pins a newer one.
    pub in_release: Option<(String, String)>,
    pub steps: Vec<String>,
    pub action: String,
}

impl Item {
    fn new(
        item: &str,
        label: &str,
        kind: Kind,
        ran_with: Option<String>,
        installed: Option<String>,
    ) -> Item {
        Item {
            item: item.into(),
            label: label.into(),
            kind,
            ran_with,
            installed,
            refs: None,
            upstream: None,
            error: None,
            class: Class::Error,
            in_release: None,
            steps: vec![],
            action: String::new(),
        }
    }

    pub fn to_json(&self) -> Value {
        let opt = |v: &Option<String>| v.as_ref().map_or(Value::Null, |s| Value::Str(s.clone()));
        let mut m = IndexMap::new();
        m.insert("item".into(), Value::Str(self.item.clone()));
        m.insert("label".into(), Value::Str(self.label.clone()));
        m.insert(
            "kind".into(),
            Value::Str(
                match self.kind {
                    Kind::File => "file",
                    Kind::Pin => "pin",
                }
                .into(),
            ),
        );
        m.insert("ran_with".into(), opt(&self.ran_with));
        m.insert("installed".into(), opt(&self.installed));
        if let Some((r, i)) = &self.refs {
            m.insert("ran_with_ref".into(), opt(r));
            m.insert("installed_ref".into(), opt(i));
        }
        if let Some(u) = &self.upstream {
            m.insert("upstream".into(), opt(u));
        }
        if let Some(e) = &self.error {
            m.insert("error".into(), Value::Str(e.clone()));
        }
        m.insert("class".into(), Value::Str(self.class.as_str().into()));
        if let Some((tag, pin)) = &self.in_release {
            m.insert(
                "in_owngenome_release".into(),
                Value::Object(IndexMap::from([
                    ("tag".to_owned(), Value::Str(tag.clone())),
                    ("pin".to_owned(), Value::Str(pin.clone())),
                ])),
            );
        }
        m.insert(
            "steps".into(),
            Value::Array(self.steps.iter().map(|s| Value::Str(s.clone())).collect()),
        );
        m.insert("action".into(), Value::Str(self.action.clone()));
        Value::Object(m)
    }
}

/// A JSON scalar as text (`117` and `"117"` alike); `None` for null and
/// for what has no text.
pub fn value_text(v: &Value) -> Option<String> {
    match v {
        Value::Str(s) => Some(s.clone()),
        Value::Int(i) => Some(i.to_string()),
        Value::Float(f) => Some(pycompat::fmt::repr(*f)),
        Value::Bool(b) => Some(if *b { "True" } else { "False" }.into()),
        _ => None,
    }
}

/// The first date in `text` as `YYYYMMDD` (`20261005`, `2026-10-05`), or
/// `None`; an impossible date (`20260230`) is no date.
pub fn date8(text: &str) -> Option<String> {
    static D: LazyLock<Regex> =
        LazyLock::new(|| Regex::new(r"(\d{4})-?(\d{2})-?(\d{2})").expect("date"));
    let c = D.captures(text)?;
    let n = |i: usize| c[i].parse::<i16>().ok();
    let (y, m, d) = (n(1)?, n(2)?, n(3)?);
    if y < 1 {
        return None;
    }
    jiff::civil::Date::new(y, i8::try_from(m).ok()?, i8::try_from(d).ok()?).ok()?;
    Some(format!("{}{}{}", &c[1], &c[2], &c[3]))
}

fn days_between(later: &str, earlier: &str) -> Option<i64> {
    let p = |s: &str| -> Option<jiff::civil::Date> {
        jiff::civil::Date::new(
            s.get(..4)?.parse().ok()?,
            s.get(4..6)?.parse().ok()?,
            s.get(6..8)?.parse().ok()?,
        )
        .ok()
    };
    let (a, b) = (p(later)?, p(earlier)?);
    Some(i64::from(b.until(a).ok()?.get_days()))
}

/// The value of `#key=` or `##key=` in the leading comment lines of a
/// (gzipped) text file.
pub fn header_value(path: &Path, key: &str) -> Option<String> {
    let mut r = og_io::open(path).ok()?;
    let prefix_ok = |line: &str| -> Option<String> {
        let rest = line.trim_start_matches('#');
        (rest.len() < line.len()).then_some(())?;
        rest.strip_prefix(key)?
            .strip_prefix('=')
            .map(|v| v.trim().to_owned())
    };
    let mut buf = Vec::new();
    for n in 0.. {
        buf.clear();
        match std::io::BufRead::read_until(&mut r, b'\n', &mut buf) {
            Ok(0) | Err(_) => return None,
            Ok(_) => {}
        }
        let line = String::from_utf8_lossy(&buf);
        let line = line.strip_suffix('\n').unwrap_or(&line);
        if !line.starts_with('#') || n > 500 {
            return None;
        }
        if let Some(v) = prefix_ok(line) {
            return Some(v);
        }
    }
    None
}

/// The date of the ClinVar file installed in `genome_dir`.
pub fn installed_clinvar(genome_dir: &Path) -> Option<String> {
    for f in CLINVAR_FILES {
        let p = genome_dir.join(f);
        if p.is_file()
            && let Some(d) = header_value(&p, "fileDate")
        {
            return date8(&d);
        }
    }
    None
}

/// A field of a manifest PGS value (`pgs_id=..;HmPOS_date=..`).
pub fn pgs_field(value: &str, field: &str) -> Option<String> {
    value
        .split(';')
        .find_map(|kv| kv.strip_prefix(field)?.strip_prefix('=').map(str::to_owned))
}

/// A `reads` entry matches a key: `KEY`, or a prefix `PREFIX*`, either with
/// `@self`.
pub fn key_matches(key: &str, read: &str) -> bool {
    let read = read.strip_suffix("@self").unwrap_or(read);
    key == read || read.strip_suffix('*').is_some_and(|p| key.starts_with(p))
}

/// The steps that read one of `keys`, every step downstream of them, then
/// the HTML report, in catalogue order. A step that reads a key only for a
/// table of its own (`KEY@self`) is marked, and nothing below it.
pub fn affected_steps(steps: &[Step], keys: &[&str]) -> Vec<String> {
    let mut direct: IndexMap<&str, bool> = IndexMap::new();
    for s in steps {
        let hits: Vec<&String> = s
            .reads
            .iter()
            .filter(|r| keys.iter().any(|k| key_matches(k, r)))
            .collect();
        if !hits.is_empty() {
            direct.insert(&s.step, hits.iter().all(|r| r.ends_with("@self")));
        }
    }
    if direct.is_empty() {
        return vec![];
    }
    let mut hit: BTreeSet<&str> = direct
        .iter()
        .filter(|(_, own)| !**own)
        .map(|(n, _)| *n)
        .collect();
    let mut grew = true;
    while grew {
        grew = false;
        for s in steps {
            if !hit.contains(s.step.as_str()) && s.upstream.iter().any(|u| hit.contains(u.as_str()))
            {
                hit.insert(&s.step);
                grew = true;
            }
        }
    }
    hit.extend(direct.keys());
    hit.insert(REPORT_STEP);
    steps
        .iter()
        .filter(|s| hit.contains(s.step.as_str()))
        .map(|s| s.step.clone())
        .collect()
}

/// A manifest or versions.env value as the release it names.
pub fn pin_value(key: &str, value: Option<&str>) -> Option<String> {
    let v = value?;
    if key.ends_with("_IMAGE") {
        return Some(image_tag(v).to_owned());
    }
    if key == "hla_database" {
        static H: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"(\d+\.\d+\.\d+)").expect("hla"));
        return H.captures(v).map(|c| c[1].to_owned());
    }
    Some(v.to_owned())
}

fn truthy(v: &Option<String>) -> Option<&str> {
    v.as_deref().filter(|s| !s.is_empty())
}

/// The class of a looked-up item; may set its error or the release that
/// pins a newer one.
pub fn classify(item: &mut Item, release: &Result<Release, String>, clinvar: bool) -> Class {
    if item.error.is_some() {
        return Class::Error;
    }
    let (ran, inst) = (truthy(&item.ran_with), truthy(&item.installed));
    if let (Some(i), Some(r)) = (inst, ran) {
        if ver(i) > ver(r) {
            return Class::RerunReady;
        }
        // An image re-pinned at the same tag (a new digest): the tags
        // compare equal, the references do not.
        if let Some((Some(rr), Some(ir))) = &item.refs
            && ver(i) == ver(r)
            && !rr.is_empty()
            && !ir.is_empty()
            && rr != ir
        {
            return Class::RerunReady;
        }
    }
    // ClinVar's date is in the manifest whenever step 06 had a file, so a run
    // without one cannot be compared with the file here. A pin with no run
    // value (HLA data not installed then) compares this checkout's.
    if clinvar && ran.is_none() {
        item.error = Some("the manifest does not say which ClinVar file the run used".into());
        return Class::Error;
    }
    let Some(base) = inst.or(ran).map(str::to_owned) else {
        item.error = Some("the manifest does not say which release the run used".into());
        return Class::Error;
    };
    let up = item.upstream.clone().flatten().unwrap_or_default();
    if ver(&up).is_empty() {
        item.error = Some("the answer names no release".into());
        return Class::Error;
    }
    if clinvar {
        return match days_between(&up, &base) {
            Some(d) if d > CLINVAR_SLACK_DAYS => Class::RefreshAvailable,
            Some(_) => Class::NoNewerRelease,
            None => {
                item.error = Some("the answer names no release".into());
                Class::Error
            }
        };
    }
    if ver(&up) <= ver(&base) {
        return Class::NoNewerRelease;
    }
    if item.kind == Kind::File {
        return Class::UpstreamOnly;
    }
    let (tag, pins) = match release {
        Err(e) => {
            item.error = Some(format!(
                "the pipeline's newest release could not be read: {e}"
            ));
            return Class::Error;
        }
        Ok(r) => r,
    };
    if let Some(rel) =
        pin_value(&item.item, pins.get(&item.item).map(String::as_str)).filter(|r| !r.is_empty())
        && ver(&rel) > ver(&base)
    {
        item.in_release = Some((tag.clone(), rel));
        return Class::InPipelineRelease;
    }
    Class::UpstreamOnly
}

/// What to do about an item.
pub fn action(item: &Item, sample: &str) -> String {
    let rerun = format!("bin/og recheck --sample {sample} --rerun");
    let inst = item.installed.as_deref().unwrap_or("None");
    match item.class {
        Class::RerunReady if item.kind == Kind::Pin => {
            match OPT_IN_SETUP.iter().find(|(k, _)| *k == item.item) {
                Some((_, opt)) => {
                    let how = match opt {
                        Some(o) => format!(
                            "run ./scripts/setup.sh {o} $GENOME_DIR (the main setup leaves this step's data out)"
                        ),
                        None => "the run pulls the image".into(),
                    };
                    format!("this checkout pins {inst}: {how}, then rerun: {rerun}")
                }
                None => format!(
                    "this checkout pins {inst}: run ./scripts/setup.sh (it installs the data), then rerun: {rerun}"
                ),
            }
        }
        Class::RerunReady => format!("rerun: {rerun}"),
        Class::RefreshAvailable => format!(
            "refresh it with ./scripts/setup.sh --refresh clinvar $GENOME_DIR, then rerun \
             ({rerun} --refresh-clinvar does both)"
        ),
        Class::InPipelineRelease => {
            let (tag, pin) = item.in_release.clone().unwrap_or_default();
            format!("update the pipeline to {tag} (it pins {pin}), run setup.sh, then rerun")
        }
        Class::UpstreamOnly if item.item.starts_with("pgs:") => format!(
            "download {} again from the PGS Catalog (docs/25-prs.md), then rerun",
            &item.item[4..]
        ),
        Class::UpstreamOnly => {
            "no pipeline release pins it yet; nothing to do until one does".into()
        }
        Class::Error => "the lookup failed; check again later".into(),
        Class::Unknown => {
            "the run's manifest does not name this pin (a run of an older version): nothing to \
                           compare; the next run records it"
                .into()
        }
        Class::NoNewerRelease => String::new(),
    }
}

/// The steps this run did not run: the ones `run-all.sh` recorded as
/// skipped (no BAM, missing data), and the ones it wrote no row for (an
/// opt-in step not in TOOLS, "not planned"). Empty when
/// `logs/run_status.tsv` records no step at all (a run of an older
/// version), so nothing is left out on a guess.
pub fn steps_not_run(sample: &str, sample_dir: &Path, steps: &[Step]) -> BTreeSet<String> {
    let Ok(files) = RunFiles::read(
        sample_dir,
        &Sources {
            traces: &[],
            events: None,
            stale: false,
        },
    ) else {
        return BTreeSet::new();
    };
    if files.status.as_ref().is_none_or(|s| s.steps.is_empty()) {
        return BTreeSet::new();
    }
    let st = files.state(sample, steps, now_epoch(), &jiff::tz::TimeZone::system());
    st.steps
        .iter()
        .filter(|s| matches!(s.state, StepState::Skipped | StepState::NotPlanned))
        .map(|s| s.step.step.clone())
        .collect()
}

pub fn now_epoch() -> f64 {
    let t = jiff::Timestamp::now();
    t.as_second() as f64 + f64::from(t.subsec_microsecond()) / 1e6
}

/// Now, as `2026-10-06T00:00:00Z`.
pub fn now_utc() -> String {
    jiff::Timestamp::now()
        .strftime("%Y-%m-%dT%H:%M:%SZ")
        .to_string()
}

/// What the check reads from this checkout.
pub struct Here<'a> {
    /// This checkout's `versions.env` pins.
    pub pins: &'a IndexMap<String, String>,
    /// This checkout's `HLA_DB_RELEASE`, when it could be read.
    pub hla: Option<&'a str>,
}

/// One check.
#[derive(Debug, Clone, PartialEq)]
pub struct Check {
    pub sample: String,
    pub checked_utc: String,
    pub manifest_written_utc: Option<String>,
    pub offline: bool,
    pub hosts: Vec<String>,
    /// The newest pipeline release's tag, or why it could not be read.
    pub release: Result<String, String>,
    pub items: Vec<Item>,
    pub affected_steps: Vec<String>,
    pub rerun_ready_steps: Vec<String>,
    /// The last rerun, kept until another one replaces it.
    pub rerun: Option<Value>,
}

impl Check {
    pub fn to_json(&self) -> Value {
        self.json_with(Value::Array(self.items.iter().map(Item::to_json).collect()))
    }

    /// The document with `items` as given.
    pub fn json_with(&self, items: Value) -> Value {
        let strs = |v: &[String]| Value::Array(v.iter().map(|s| Value::Str(s.clone())).collect());
        let mut m = IndexMap::new();
        m.insert("schema_version".into(), Value::Int(SCHEMA_VERSION));
        m.insert("sample".into(), Value::Str(self.sample.clone()));
        m.insert("checked_utc".into(), Value::Str(self.checked_utc.clone()));
        m.insert(
            "manifest_written_utc".into(),
            self.manifest_written_utc
                .clone()
                .map_or(Value::Null, Value::Str),
        );
        m.insert("offline".into(), Value::Bool(self.offline));
        m.insert("hosts".into(), strs(&self.hosts));
        let (k, v) = match &self.release {
            Ok(t) => ("tag", t),
            Err(e) => ("error", e),
        };
        m.insert(
            "pipeline_release".into(),
            Value::Object(IndexMap::from([(k.to_owned(), Value::Str(v.clone()))])),
        );
        m.insert("items".into(), items);
        m.insert("affected_steps".into(), strs(&self.affected_steps));
        m.insert("rerun_ready_steps".into(), strs(&self.rerun_ready_steps));
        if let Some(r) = &self.rerun {
            m.insert("rerun".into(), r.clone());
        }
        Value::Object(m)
    }

    pub fn errors(&self) -> usize {
        self.items
            .iter()
            .filter(|i| i.class == Class::Error)
            .count()
    }
}

/// Checks `sample` of `genome_dir` against `answers`.
pub fn check(
    sample: &str,
    genome_dir: &Path,
    answers: &dyn Answers,
    steps: &[Step],
    here: &Here<'_>,
) -> Check {
    let sample_dir = genome_dir.join(sample);
    let man = og_runlog::files::read_manifest(Some(&sample_dir.join("run_manifest.tsv")))
        .unwrap_or_default();
    let get = |sec: &str, key: &str| {
        man.get(&(sec.to_owned(), key.to_owned()))
            .map(String::as_str)
    };
    let release = answers.release().map_err(|e| e.to_string());
    let mut items: Vec<Item> = vec![];
    let finish =
        |mut item: Item, answer: Option<og_fresh::Lookup<Option<String>>>, clinvar: bool| -> Item {
            if let Some(a) = answer {
                match a {
                    Ok(v) => item.upstream = Some(v),
                    Err(e) => item.error = Some(e.to_string()),
                }
                item.class = classify(&mut item, &release, clinvar);
            }
            item.steps = affected_steps(steps, &[&item.item]);
            item.action = action(&item, sample);
            item
        };

    let cv = Item::new(
        "clinvar_file_date",
        "ClinVar (step 06's file)",
        Kind::File,
        get("data", "clinvar_file_date").and_then(date8),
        installed_clinvar(genome_dir),
    );
    items.push(finish(cv, Some(answers.clinvar()), true));
    for (key, label) in PINNED {
        let ran_raw = if *key == "hla_database" {
            get("data", key)
        } else {
            get("version", key)
        };
        let inst = if *key == "hla_database" {
            here.hla.map(str::to_owned)
        } else {
            pin_value(key, here.pins.get(*key).map(String::as_str))
        };
        let mut item = Item::new(key, label, Kind::Pin, pin_value(key, ran_raw), inst);
        if key.ends_with("_IMAGE") {
            let nonempty = |v: Option<&str>| v.filter(|s| !s.is_empty()).map(str::to_owned);
            item.refs = Some((
                nonempty(get("version", key)),
                nonempty(here.pins.get(*key).map(String::as_str)),
            ));
        }
        if *key == "hla_database" && here.hla.is_none() {
            item.error = Some(
                "could not read this checkout's HLA_DB_RELEASE (scripts/lib/common.sh)".into(),
            );
        }
        let answer = answers.pinned(key);
        items.push(finish(item, Some(answer), false));
    }
    // The other versions.env pins a step reads (steps 50 to 53): no upstream
    // lookup exists for them, so nothing is asked; a pin here that differs
    // from the run's (an image by its full reference, so a new digest at the
    // same tag counts) makes its step rerun-ready. A pin the run's manifest
    // lacks is unknown, never no-newer-release. A pin read only by steps this
    // run did not run is left out: a rerun would not run that step either.
    let skipped = steps_not_run(sample, &sample_dir, steps);
    let mut seen: Vec<&str> = vec![];
    for r in steps.iter().flat_map(|s| s.reads.iter()) {
        let key = r.strip_suffix("@self").unwrap_or(r);
        if seen.contains(&key) {
            continue;
        }
        seen.push(key);
        if PINNED.iter().any(|(k, _)| *k == key) || !here.pins.contains_key(key) {
            continue;
        }
        let readers: Vec<&str> = steps
            .iter()
            .filter(|s| s.reads.iter().any(|r| key_matches(key, r)))
            .map(|s| s.step.as_str())
            .collect();
        if !readers.is_empty() && readers.iter().all(|n| skipped.contains(*n)) {
            continue;
        }
        let ran_ref = get("version", key)
            .filter(|s| !s.is_empty())
            .map(str::to_owned);
        let inst_ref = here.pins.get(key).filter(|s| !s.is_empty()).cloned();
        let mut item = Item::new(
            key,
            &format!("{key} (no upstream lookup)"),
            Kind::Pin,
            pin_value(key, ran_ref.as_deref()),
            pin_value(key, inst_ref.as_deref()),
        );
        item.class = match (&ran_ref, &inst_ref) {
            (Some(r), Some(i)) if r != i => Class::RerunReady,
            (Some(_), Some(_)) => Class::NoNewerRelease,
            _ => Class::Unknown,
        };
        if key.ends_with("_IMAGE") {
            item.refs = Some((ran_ref, inst_ref));
        }
        items.push(finish(item, None, false));
    }
    let mut pgs: Vec<(&String, &String)> = man
        .iter()
        .filter(|((sec, key), _)| sec == "data" && key.starts_with("pgs:"))
        .map(|((_, k), v)| (k, v))
        .collect();
    pgs.sort();
    for (key, value) in pgs {
        let id = &key[4..];
        let now = header_value(
            &genome_dir.join("prs_scores").join(format!("{id}.txt.gz")),
            "HmPOS_date",
        );
        let item = Item::new(
            key,
            &format!("PGS Catalog score {id}"),
            Kind::File,
            pgs_field(value, "HmPOS_date").as_deref().and_then(date8),
            now.as_deref().and_then(date8),
        );
        items.push(finish(item, Some(answers.pgs(id)), false));
    }
    let named = |classes: &[Class]| -> Vec<String> {
        let keys: Vec<&str> = items
            .iter()
            .filter(|i| classes.contains(&i.class))
            .map(|i| i.item.as_str())
            .collect();
        affected_steps(steps, &keys)
    };
    let affected = named(&[Class::RerunReady, Class::RefreshAvailable]);
    let ready = named(&[Class::RerunReady]);
    Check {
        sample: sample.into(),
        checked_utc: now_utc(),
        manifest_written_utc: get("run", "written_utc").map(str::to_owned),
        offline: answers.offline(),
        hosts: answers.hosts(),
        release: release.map(|(t, _)| t),
        items,
        affected_steps: affected,
        rerun_ready_steps: ready,
        rerun: None,
    }
}

/// The check as `og recheck` prints it.
pub fn print_check(c: &Check, out: &mut dyn std::io::Write) -> std::io::Result<()> {
    writeln!(
        out,
        "Database re-check of {} ({}):",
        c.sample, c.checked_utc
    )?;
    let or_unknown = |v: &Option<String>| truthy(v).unwrap_or("unknown").to_owned();
    for i in &c.items {
        let mut detail = match &i.error {
            Some(e) if !e.is_empty() => e.clone(),
            _ => format!(
                "run used {}, here {}, upstream {}",
                or_unknown(&i.ran_with),
                or_unknown(&i.installed),
                or_unknown(&i.upstream.clone().flatten())
            ),
        };
        if i.error.is_none()
            && let Some((Some(r), Some(h))) = &i.refs
            && r != h
            && i.ran_with == i.installed
        {
            detail += &format!("; same tag, new image: run used {r}, here {h}");
        }
        writeln!(out, "  {}: {} ({detail})", i.label, i.class.as_str())?;
        if i.class != Class::NoNewerRelease && !i.action.is_empty() {
            writeln!(out, "    {}", i.action)?;
        }
        if matches!(i.class, Class::RerunReady | Class::RefreshAvailable) && !i.steps.is_empty() {
            writeln!(
                out,
                "    steps it changes: {} and the text report",
                i.steps.join(", ")
            )?;
        }
    }
    if c.affected_steps.is_empty() {
        writeln!(out, "Nothing here is newer than what the run used.")
    } else {
        writeln!(
            out,
            "A rerun would redo steps {} and the text report.",
            c.affected_steps.join(", ")
        )
    }
}

#[cfg(test)]
mod tests;
