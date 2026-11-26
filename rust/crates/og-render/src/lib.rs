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
//! Every part of the report is ported: the frame (chapters, contents,
//! section frames, the run facts, what the report leaves out), every section
//! [`chapters::writer`] names, the search index and the charts. A section
//! that cannot be written says so in its place in both reports, and the rest
//! of the report is still written ([`render`]). `og report charts` writes the
//! chart specs on its own, and `og report render --charts-json` writes the
//! same file for the steps.

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
