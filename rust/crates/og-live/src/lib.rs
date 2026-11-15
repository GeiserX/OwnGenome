//! The live site: its server (`og live serve`, [`server`]), each person's
//! results page (`og live results`, [`results`]), the family page
//! (`og live family`, [`family`]), the public demo as static pages
//! (`og live demo`, [`demo`]) and the reads viewer (`og view serve`,
//! [`viewer`]).
//!
//! Not ported yet: each verb says so and exits 2 until its unit of the
//! 1.0.0 plan lands. Every module those units name is declared here now,
//! so they add code without editing the same files.

pub mod auth;
pub mod charts;
pub mod csp;
pub mod demo;
pub mod family;
pub mod headlines;
pub mod ingest;
pub mod map;
pub mod ops;
pub mod plugins;
pub mod recheck;
pub mod results;
pub mod runs;
pub mod search;
pub mod server;
pub mod state;
pub mod static_files;
pub mod veil;
pub mod viewer;
pub mod widgets;
