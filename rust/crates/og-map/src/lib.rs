//! The genome map: one person's map document, built on request from the
//! run folder and the repository's reference assets (`README.md` holds the
//! format).
//!
//! Phase 1: the reference bands, the chromosome lengths, centromeres and
//! masked spans, the findings layer with the veil tags (`veil`), and the
//! runs-of-homozygosity layer. The gene files (`genes`) are built once per
//! GENCODE release, apart from any run.

pub mod doc;
pub mod facts;
pub mod findings;
pub mod genes;
pub mod layers;
pub mod predicate;
pub mod reference;
pub mod stamp;
pub mod text;
pub mod veil;

use facts::{Facts, Sex};
use findings::Placed;
use layers::roh::{self, Roh};
use reference::Reference;
use std::fmt;
use std::path::{Path, PathBuf};

/// The document format version.
pub const VERSION: u32 = 1;

/// Every path under GENOME_DIR the map reads, as written in
/// `assets/output_contract.tsv` (`{sample}` is the run).
pub const INPUT_PATHS: &[&str] = &[roh::PATHS[0], roh::PATHS[1], facts::SEX_PATH];

/// A layer's state: three answers, never two. A missing input is
/// `missing`; an input that cannot be read is `unreadable`; a read input
/// with no rows is `ok` with a count of 0.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum State {
    Ok,
    Missing,
    Unreadable,
}

impl State {
    pub fn word(self) -> &'static str {
        match self {
            State::Ok => "ok",
            State::Missing => "missing",
            State::Unreadable => "unreadable",
        }
    }
}

/// What a layer says about itself, beside its rows.
#[derive(Debug, Clone, PartialEq)]
pub struct Layer {
    pub state: State,
    pub step: &'static str,
    pub step_name: &'static str,
    /// The file read or looked for, relative to GENOME_DIR.
    pub source: String,
    /// Why it is missing, or `file:line: reason` when unreadable.
    pub note: Option<String>,
    pub count: usize,
    /// (rows read, rows on the map), when the file was read.
    pub totals: Option<(usize, usize)>,
    pub caveat: String,
    /// The sentence for a read file with no rows.
    pub empty: Option<String>,
}

impl Layer {
    pub fn new(step: &'static str, step_name: &'static str) -> Layer {
        Layer {
            state: State::Missing,
            step,
            step_name,
            source: String::new(),
            note: None,
            count: 0,
            totals: None,
            caveat: String::new(),
            empty: None,
        }
    }
}

/// One person's map, in memory: what `og map build` writes.
#[derive(Debug, Clone)]
pub struct Model {
    pub run: String,
    pub stamp: String,
    pub made_by: String,
    pub sex: Sex,
    pub reference: Reference,
    /// The GENCODE release installed under the root (`v50`), if any.
    pub gencode: Option<String>,
    pub findings: Placed,
    pub roh: Roh,
}

/// Why no document was built.
#[derive(Debug)]
pub enum BuildError {
    /// The command line or stdin is wrong (exit 2).
    Usage(String),
    /// There is no map for this run (exit 3).
    NoMap(String),
    /// Something that should not fail did (exit 1).
    Internal(String),
}

impl fmt::Display for BuildError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            BuildError::Usage(s) | BuildError::NoMap(s) | BuildError::Internal(s) => f.write_str(s),
        }
    }
}

impl std::error::Error for BuildError {}

/// One `og map build` request.
#[derive(Debug, Clone)]
pub struct Request<'a> {
    /// GENOME_DIR: one folder per run.
    pub root: &'a Path,
    pub run: &'a str,
    /// The repository's `assets/`.
    pub assets: &'a Path,
    /// The og tree hash: the code that builds the document.
    pub made_by: &'a str,
    /// The facts body read from stdin.
    pub stdin: &'a [u8],
}

/// A run name as the live site accepts it (`validate_sample` in
/// scripts/lib/common.sh): a letter or digit, then letters, digits, `.`, `_`
/// or `-`.
pub fn valid_run(name: &str) -> bool {
    let mut b = name.bytes();
    b.next().is_some_and(|c| c.is_ascii_alphanumeric())
        && b.all(|c| c.is_ascii_alphanumeric() || matches!(c, b'.' | b'_' | b'-'))
}

impl Request<'_> {
    fn check(&self) -> Result<(), BuildError> {
        if !valid_run(self.run) {
            return Err(BuildError::Usage(format!("not a run name: {:?}", self.run)));
        }
        if !self.root.join(self.run).is_dir() {
            return Err(BuildError::NoMap(format!(
                "no map for this run: no folder {} under the root",
                self.run
            )));
        }
        Ok(())
    }

    /// Every file the document depends on, present or not.
    pub fn inputs(&self) -> Vec<PathBuf> {
        let mut v = vec![
            self.assets.join(reference::BAND_FILE),
            self.root
                .join(facts::SEX_PATH.replace("{sample}", self.run)),
            // The genes file that places gene-only findings; while none is
            // installed, the folder setup installs it in, so its arrival
            // changes the stamp.
            genes::locate(self.root).map_or_else(|_| self.root.join("reference"), |i| i.gtf),
        ];
        v.extend(
            roh::candidates(self.root, self.run)
                .into_iter()
                .map(|(_, p)| p),
        );
        v
    }

    /// The stamp of this request's inputs.
    pub fn stamp(&self) -> String {
        stamp::stamp(VERSION, self.made_by, self.stdin, &self.inputs())
    }
}

/// Builds the person document's model.
pub fn build_person(req: &Request<'_>) -> Result<Model, BuildError> {
    req.check()?;
    build_with_stamp(req, req.stamp())
}

/// [`build_person`] with the stamp already taken, so a cached document
/// carries the stamp it is stored under. The caller has run `req.check()`.
fn build_with_stamp(req: &Request<'_>, stamp: String) -> Result<Model, BuildError> {
    let Facts { sex, findings } = facts::parse(req.stdin).map_err(BuildError::Usage)?;
    let sex = facts::sex_file(req.root, req.run).unwrap_or(sex);
    let reference = Reference::load(req.assets).map_err(|e| BuildError::Internal(e.to_string()))?;
    let findings = findings::read(&findings, req.root, req.run, &reference);
    let gencode = genes::locate(req.root).ok().map(|i| i.release);
    let roh = roh::read(req.root, req.run, &reference, sex);
    Ok(Model {
        run: req.run.to_owned(),
        stamp,
        made_by: req.made_by.to_owned(),
        sex,
        reference,
        gencode,
        findings,
        roh,
    })
}

/// The person document's bytes: from the cache under `cache_dir` when its
/// stamp is there, else built and stored. A cache that cannot be written
/// is reported on stderr and the document is still returned.
pub fn build_cached(req: &Request<'_>, cache_dir: &Path) -> Result<Vec<u8>, BuildError> {
    req.check()?;
    let stamp = req.stamp();
    let path = stamp::cache_path(cache_dir, req.run, &stamp);
    if let Some(doc) = stamp::lookup(&path) {
        return Ok(doc);
    }
    let model = build_with_stamp(req, stamp)?;
    let doc = doc::write(&model).map_err(|e| BuildError::Internal(e.to_string()))?;
    if let Err(e) = stamp::store(&path, &doc) {
        eprintln!("og map build: the cache was not written: {e}");
    }
    Ok(doc)
}
