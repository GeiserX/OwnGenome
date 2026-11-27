//! What `scripts/setup.sh` runs once per install, in the `og` binary:
//!
//! - [`install`] (`og setup install`, `og setup ported`, `og setup flags`):
//!   setup.sh's installers run by og, with the same arguments, files,
//!   messages and exit codes; a flag og has not ported yet is run by
//!   setup.sh's bash. [`components`] is the one table of setup's flags;
//!   [`engine`] builds the container runs, [`images`] pulls and tags the
//!   images, [`lock`] is the installers' lock and [`disk`] the disk checks
//!   before a large install. [`check`] is `og setup check`, to come.
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
//!
//! [`zip`] reads the archives the last two take. [`revel`] (`og setup
//! revel`) and [`roh_af`] (`og setup roh-af`) are not ported yet.

pub mod check;
pub mod components;
pub mod disk;
pub mod engine;
pub mod fetch;
pub mod genetic_map;
pub mod images;
pub mod install;
pub mod lock;
pub mod megane_kmers;
pub mod plan;
pub mod prune;
pub mod revel;
pub mod roh_af;
pub mod strchive;
pub mod unzip;
pub mod zip;
