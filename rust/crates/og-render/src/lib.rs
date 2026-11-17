//! The reports of one run, from its `summary.json`: the HTML and text
//! reports (`og report render`, [`render`]), the run's search index
//! (`og report index`, [`index`]), the chart specs (`og report charts`,
//! [`charts`]) and the report's sections, one module each ([`sections`]).
//!
//! The report is a pure function of the summary. [`model`] reads it with
//! Python's rules and holds the document the sections build: cards of typed
//! lines. [`chapters`] lays the sections out, [`text`] and [`html`] write the
//! two reports, [`notes`] holds the caveats, [`escape`] the HTML escaping and
//! [`css`] the page's style.
//!
//! Ported so far: the frame (chapters, contents, section frames, the run
//! facts, what the report leaves out) and the quality-control and variant
//! sections, and the search index. A section not ported yet is left out of
//! both reports; the charts verb says it is not ported and exits 2.

pub mod chapters;
pub mod charts;
pub mod css;
pub mod escape;
pub mod html;
pub mod index;
pub mod model;
pub mod notes;
pub mod render;
pub mod sections;
pub mod text;
