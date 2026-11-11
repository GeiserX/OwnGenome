//! What Nextflow writes about a run, read as OwnGenome reads it.
//!
//! - [`trace`]: a trace file (`trace_*.txt`), one header-keyed row per task,
//!   streamed; old traces (the 14 default columns) and new ones (with the
//!   columns `nextflow.config` adds) alike.
//! - [`name`]: a task name (`'SV:MANTA (s1)'`) split into its process and
//!   its tag.
//! - [`units`]: durations, memory sizes and times as a trace or a weblog
//!   message writes them, human (`1h 2m 3s`, `18.3 GB`) or raw (ms, bytes).
//! - [`time`]: the clock formats the run files use, parsed and written as
//!   Python's `datetime` does: naive trace times in the local zone, ISO
//!   8601 UTC out.

pub mod name;
pub mod time;
pub mod trace;
pub mod units;

pub use name::task_name;
pub use trace::{Row, for_each_row};
