//! `og setup chrm-shift` and `og vcf shift-back`: the chrM reference
//! rotated by 8000 bp for the shifted-reference pass, and the calls on it
//! mapped back to chrM's own positions.
//!
//! Not ported yet: until its unit of the 1.0.0 plan lands, each verb says
//! so on stderr and exits 2, before reading any argument.

use std::io::Write;

/// `og setup chrm-shift ARGS`: not ported yet, so one line on stderr and exit 2.
pub fn chrm_shift_cli(_args: &[String], _out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let _ = writeln!(err, "og setup chrm-shift: not ported yet");
    2
}

/// `og vcf shift-back ARGS`: not ported yet, so one line on stderr and exit 2.
pub fn shift_back_cli(_args: &[String], _out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let _ = writeln!(err, "og vcf shift-back: not ported yet");
    2
}
