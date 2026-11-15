//! `og dev demo-sample`: the made-up demo sample the docs, the pictures and
//! the live demo use, which `tests/demo/make_demo_sample.py` writes today.
//!
//! Not ported yet: until its unit of the 1.0.0 plan lands, the verb says so
//! on stderr and exits 2, before reading any argument.

pub mod ancestry_paint;
pub mod clinical;
pub mod clinical_clinvar;
pub mod family;
pub mod lineage;
pub mod pgx;
pub mod qc;
pub mod sv;
pub mod traits;

use std::io::Write;

/// `og dev demo-sample ARGS`: not ported yet, so one line on stderr and exit 2.
pub fn cli(_args: &[String], _out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let _ = writeln!(err, "og dev demo-sample: not ported yet");
    2
}
