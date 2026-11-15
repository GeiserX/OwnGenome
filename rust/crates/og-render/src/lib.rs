//! The reports of one run, from its `summary.json`: the HTML and text
//! reports (`og report render`, [`render`]), the run's search index
//! (`og report index`, [`index`]), the chart specs (`og report charts`,
//! [`charts`]) and the report's sections, one module each ([`sections`]).
//!
//! Not ported yet: each verb says so and exits 2 until its unit of the
//! 1.0.0 plan lands. The crate exists now so those units add code without
//! editing the same files.

pub mod charts;
pub mod index;
pub mod render;
pub mod sections;
