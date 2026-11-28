//! What `scripts/setup.sh` runs once per install, in the `og` binary:
//!
//! - [`fetch`] (`og setup fetch`, `og setup source`, `og setup path`): every
//!   download of setup, checked and kept, from the table
//!   `assets/data-sources.tsv`; [`plan`] (`og setup plan`) lists what a
//!   setup will download before it starts, [`prune`] (`og setup prune`)
//!   removes the kept raw downloads.
//! - [`strchive`] (`og str catalog`): the STRchive repeat catalog of steps
//!   09 and 9b, one file for ExpansionHunter and Stranger.
//! - [`unzip`] (`og setup unzip`): Exomiser's data zips unpacked, each
//!   member streamed to disk.
//! - [`genetic_map`] (`og setup genetic-map`): the PLINK GRCh38 genetic map
//!   per chromosome for steps 26b and 51, written atomically.
//! - [`revel`] (`og setup revel`): the REVEL release projected onto GRCh38,
//!   the score table of step 30.
//! - [`megane_kmers`] (`og setup megane-kmers`): MEGAnE's k-mer index of the
//!   reference for step 50.
//!
//! [`zip`] reads the archives `unzip` and `genetic_map` take. Step 11's
//! allele-frequency panel needs no tool of its own: `og setup fetch`
//! downloads and checks it, and setup.sh indexes it with tabix.

pub mod fetch;
pub mod genetic_map;
pub mod megane_kmers;
pub mod plan;
pub mod prune;
pub mod revel;
pub mod strchive;
pub mod unzip;
pub mod zip;
