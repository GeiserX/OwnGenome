//! The typed summary: section states, the run status, the manifest, and
//! how each is written into summary.json.
//!
//! Every object is built in the key order the Python collector inserted
//! its keys, because summary.json is compared with key order.

use super::sections::Data;
use super::time::utc_minute;
use indexmap::IndexMap;
use pycompat::json::Value;
use std::collections::HashMap;

/// summary.json's `schema_version`.
pub const SCHEMA_VERSION: i128 = 2;

/// Steps run-all.sh records in `logs/run_status.tsv`: the rows of
/// `assets/steps.tsv` with `recorded=yes` (a test keeps the two equal). A
/// section of another step (variant calling, step 03, which a run reuses on
/// purpose) is never stale.
pub const RUN_ALL_STEPS: &[&str] = &[
    "03a", "03b", "03c", "03d", "04", "04b", "05", "06", "07", "08", "09", "09b", "10", "11", "12",
    "13", "14", "15", "16", "16b", "17", "18", "19", "20", "21", "22", "23", "25", "26", "26b",
    "27", "28", "29", "30", "31", "32", "33", "35", "35b", "37", "38", "39", "40", "41", "42",
    "43", "44", "45", "46", "47", "48", "49", "50", "51", "52", "53",
];

/// What a section's answer is worth.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum State {
    /// Read from its file.
    Ok,
    /// No file: the step has not run.
    Missing,
    /// No file, and the step failed in the last run.
    Failed,
    /// The file is older than the last run, and its step was not ok in it.
    Stale,
    /// The file is there but could not be read; the note says why.
    Unreadable,
}

impl State {
    /// The word summary.json uses.
    pub fn word(self) -> &'static str {
        match self {
            State::Ok => "ok",
            State::Missing => "missing",
            State::Failed => "failed",
            State::Stale => "stale",
            State::Unreadable => "unreadable",
        }
    }
}

/// One section's place in summary.json.
#[derive(Debug, Clone, Copy)]
pub struct Spec {
    pub key: &'static str,
    pub title: &'static str,
    pub step: &'static str,
}

/// The sections `og summary collect --only a` reads, in summary.json's
/// order (the start of the Python `SECTIONS`).
pub const GROUP_A: [Spec; 7] = [
    Spec {
        key: "coverage",
        title: "Coverage (mosdepth)",
        step: "16b",
    },
    Spec {
        key: "callability",
        title: "Gene callability (mosdepth)",
        step: "16b",
    },
    Spec {
        key: "sex_check",
        title: "Sex check (indexcov)",
        step: "16",
    },
    Spec {
        key: "sample_qc",
        title: "Sample identity and contamination (somalier, VerifyBamID2)",
        step: "33",
    },
    Spec {
        key: "mosaic_sex_loss",
        title: "Mosaic loss of chromosome Y or X",
        step: "33",
    },
    Spec {
        key: "read_count",
        title: "Reads kept by alignment",
        step: "02",
    },
    Spec {
        key: "variants",
        title: "Variant calling",
        step: "03",
    },
];

/// One section of the summary.
#[derive(Debug)]
pub struct Section {
    pub spec: Spec,
    pub state: State,
    /// The file it was read from, relative to the sample folder.
    pub source: Option<String>,
    pub file_date: Option<String>,
    pub note: Option<String>,
    /// Present for `ok` and `stale` only.
    pub data: Option<Data>,
}

impl Section {
    /// Sets the state from what the reader returned, the file's time and
    /// the run status, as the Python collector did.
    pub fn new(
        spec: Spec,
        read: Result<Option<(String, f64, Data)>, String>,
        run: Option<&RunStatus>,
    ) -> Section {
        let mut s = Section {
            spec,
            state: State::Missing,
            source: None,
            file_date: None,
            note: None,
            data: None,
        };
        match read {
            Err(note) => {
                s.state = State::Unreadable;
                s.note = Some(note);
            }
            Ok(None) => {
                if let Some(run) = run
                    && run
                        .steps
                        .get(spec.step)
                        .is_some_and(|r| r.starts_with("failed"))
                {
                    s.state = State::Failed;
                    s.note = Some(format!(
                        "step {} failed in the run of {}; its log is in logs/",
                        spec.step,
                        run.started_utc
                            .as_deref()
                            .filter(|t| !t.is_empty())
                            .unwrap_or("unknown date")
                    ));
                }
            }
            Ok(Some((source, mtime, data))) => {
                s.state = State::Ok;
                s.source = Some(source);
                s.file_date = utc_minute(mtime);
                s.data = Some(data);
                if let Some(run) = run
                    && let Some(started) = run.started_epoch
                    && RUN_ALL_STEPS.contains(&spec.step)
                {
                    let result = run.steps.get(spec.step).map(String::as_str);
                    if result != Some("ok") && is_stale(mtime, started) {
                        s.state = State::Stale;
                        s.note = Some(format!(
                            "from an earlier run ({}); step {} was {} in the run of {}",
                            utc_minute(mtime).unwrap_or_default(),
                            spec.step,
                            result.filter(|r| !r.is_empty()).unwrap_or("not run"),
                            run.started_utc.as_deref().unwrap_or("None"),
                        ));
                    }
                }
            }
        }
        s
    }

    pub fn to_json(&self) -> Value {
        let mut m = IndexMap::new();
        m.insert("title".into(), Value::Str(self.spec.title.into()));
        m.insert("step".into(), Value::Str(self.spec.step.into()));
        m.insert("state".into(), Value::Str(self.state.word().into()));
        m.insert("source".into(), opt_str(self.source.as_deref()));
        m.insert("file_date".into(), opt_str(self.file_date.as_deref()));
        m.insert("note".into(), opt_str(self.note.as_deref()));
        m.insert(
            "data".into(),
            self.data
                .as_ref()
                .map_or_else(|| Value::Object(IndexMap::new()), Data::to_json),
        );
        Value::Object(m)
    }
}

/// A file written before the run started is from an earlier run. A file
/// written in the very second the run started is the run's own.
pub fn is_stale(mtime: f64, started_epoch: f64) -> bool {
    mtime < started_epoch
}

/// `logs/run_status.tsv`, which run-all.sh writes.
#[derive(Debug, Default)]
pub struct RunStatus {
    pub started_epoch: Option<f64>,
    pub started_utc: Option<String>,
    pub declared_sex: Option<String>,
    /// Each step's result word (`ok`, `skipped ...`, `failed (exit 1)`).
    pub steps: HashMap<String, String>,
}

impl RunStatus {
    /// Reads the table: `meta` and `step` rows of at least three cells;
    /// lines starting with `#` are comments. A `started_epoch` that is not
    /// a number leaves the one before it. The start time is written from
    /// the epoch; without one, from a `meta started_utc` row.
    pub fn parse(text: &str) -> RunStatus {
        let mut st = RunStatus::default();
        for line in text.lines() {
            let c: Vec<&str> = line.split('\t').collect();
            if c.len() < 3 || c[0].starts_with('#') {
                continue;
            }
            match (c[0], c[1]) {
                ("meta", "started_epoch") => {
                    if let Ok(v) = pycompat::num::parse_float(c[2]) {
                        st.started_epoch = Some(v);
                    }
                }
                ("meta", "started_utc") => st.started_utc = Some(c[2].to_owned()),
                ("meta", "declared_sex") => st.declared_sex = Some(c[2].to_owned()),
                ("step", step) => {
                    st.steps.insert(step.to_owned(), c[2].to_owned());
                }
                _ => {}
            }
        }
        // A start time that names no date (NaN, infinity) is no start time:
        // Python stopped with a traceback here and wrote no summary.
        if let Some(e) = st.started_epoch {
            match utc_minute(e) {
                Some(t) => st.started_utc = Some(t),
                None => st.started_epoch = None,
            }
        }
        st
    }
}

/// `run_manifest.tsv` (bin/write_manifest.sh): section, key, value; a
/// value with tabs keeps them.
#[derive(Debug, Default)]
pub struct Manifest(pub Vec<[String; 3]>);

impl Manifest {
    pub fn parse(text: &str) -> Manifest {
        Manifest(
            text.lines()
                .filter_map(|line| {
                    let c: Vec<&str> = line.split('\t').collect();
                    (c.len() >= 3 && !c[0].starts_with('#'))
                        .then(|| [c[0].to_owned(), c[1].to_owned(), c[2..].join("\t")])
                })
                .collect(),
        )
    }

    /// The value of the first row of `section` and `key`.
    pub fn value(&self, section: &str, key: &str) -> Option<&str> {
        self.0
            .iter()
            .find(|r| r[0] == section && r[1] == key)
            .map(|r| r[2].as_str())
    }

    pub fn to_json(&self) -> Value {
        Value::Array(
            self.0
                .iter()
                .map(|r| Value::Array(r.iter().map(|c| Value::Str(c.clone())).collect()))
                .collect(),
        )
    }
}

/// The summary of one sample, as far as the ported sections go.
#[derive(Debug)]
pub struct Summary {
    pub sample: String,
    pub generated_utc: String,
    pub run: Option<RunStatus>,
    pub manifest: Option<Manifest>,
    pub declared_sex: Option<String>,
    pub sections: Vec<Section>,
}

impl Summary {
    pub fn to_json(&self) -> Value {
        let manifest = self.manifest.as_ref();
        let db = |key: &str| opt_str(manifest.and_then(|m| m.value("data", key)));
        let mut run = IndexMap::new();
        run.insert("status_file".into(), Value::Bool(self.run.is_some()));
        run.insert(
            "started_utc".into(),
            opt_str(self.run.as_ref().and_then(|r| r.started_utc.as_deref())),
        );
        run.insert("declared_sex".into(), opt_str(self.declared_sex.as_deref()));
        let mut databases = IndexMap::new();
        databases.insert("clinvar_release".into(), db("clinvar_file_date"));
        databases.insert("cpsr_clinvar_release".into(), db("cpsr_clinvar_release"));
        databases.insert("hla_database".into(), db("hla_database"));
        let mut m = IndexMap::new();
        m.insert("schema_version".into(), Value::Int(SCHEMA_VERSION));
        m.insert("sample".into(), Value::Str(self.sample.clone()));
        m.insert(
            "generated_utc".into(),
            Value::Str(self.generated_utc.clone()),
        );
        m.insert("run".into(), Value::Object(run));
        m.insert("databases".into(), Value::Object(databases));
        m.insert(
            "manifest".into(),
            manifest.map_or(Value::Array(Vec::new()), Manifest::to_json),
        );
        m.insert(
            "sections".into(),
            Value::Object(
                self.sections
                    .iter()
                    .map(|s| (s.spec.key.to_owned(), s.to_json()))
                    .collect(),
            ),
        );
        Value::Object(m)
    }
}

pub fn opt_str(s: Option<&str>) -> Value {
    s.map_or(Value::Null, |s| Value::Str(s.to_owned()))
}
