//! What `scripts/setup.sh` runs once per install, in the `og` binary:
//!
//! - [`strchive`] (`og str catalog`): the STRchive repeat catalog of steps
//!   09 and 9b, one file for ExpansionHunter and Stranger.
//! - [`unzip`] (`og setup unzip`): Exomiser's data zips unpacked, each
//!   member streamed to disk.
//! - [`genetic_map`] (`og setup genetic-map`): the PLINK GRCh38 genetic map
//!   per chromosome for steps 26b and 51, written atomically.
//!
//! [`zip`] reads the archives the last two take.

pub mod genetic_map;
pub mod strchive;
pub mod unzip;
pub mod zip;
