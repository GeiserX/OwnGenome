//! A run, read from the files that say how it goes, never from a result.
//!
//! - [`steps`]: `assets/steps.tsv`, the catalogue of steps, their processes,
//!   sections, upstream steps and usual hours.
//! - [`files`]: `logs/run_status.tsv` (run-all.sh), `nextflow/process_status.tsv`
//!   (main.nf), weblog messages (one JSON object per line), the `source` of
//!   each `summary.json` section, `run_manifest.tsv`.
//! - [`tasks`]: the sample's tasks from trace rows and weblog messages, one
//!   per task, the last attempt winning.
//! - [`state`]: where the run stands, step by step (`og run state`).
//! - [`eta`]: this box's step history and the hours a run has left
//!   (`og run eta`).
//! - [`headroom`]: each task's peak memory against its request, and how full
//!   the work and results disks are (`og run headroom`).
//! - [`delta`]: what changed since the previous run of the sample
//!   (`og run delta`), and [`delta::newest_history`], the one rule for which
//!   history folder is "the previous run".
//! - [`cli`]: the `og run` verbs.
//!
//! Every reader takes what a half-written or missing file gives it: a run
//! file that is not there is a value (`None`), not an error.

pub mod cli;
pub mod delta;
pub mod eta;
pub mod files;
pub mod headroom;
pub mod host;
pub mod state;
pub mod steps;
pub mod tasks;
mod util;
