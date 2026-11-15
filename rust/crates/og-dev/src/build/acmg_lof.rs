//! `og dev build acmg-lof`: the loss-of-function column of
//! `assets/acmg_sf.tsv`, built from pinned ClinGen scores.
//!
//! Not ported yet: until its unit of the 1.0.0 plan lands, the verb says so
//! on stderr and exits 2, before reading any argument.

use std::io::Write;

/// `og dev build acmg-lof ARGS`: not ported yet, so one line on stderr and exit 2.
pub fn cli(_args: &[String], _out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let _ = writeln!(err, "og dev build acmg-lof: not ported yet");
    2
}
