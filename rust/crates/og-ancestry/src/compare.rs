//! `og archaic compare` and `og archaic genomes`: step 43's archaic totals
//! beside published calls, which `bin/archaic_compare.py` writes today.
//!
//! Not ported yet: until its unit of the 1.0.0 plan lands, each verb says
//! so on stderr and exits 2, before reading any argument.

use std::io::Write;

/// `og archaic compare ARGS`: not ported yet, so one line on stderr and exit 2.
pub fn compare_cli(_args: &[String], _out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let _ = writeln!(err, "og archaic compare: not ported yet");
    2
}

/// `og archaic genomes ARGS`: not ported yet, so one line on stderr and exit 2.
pub fn genomes_cli(_args: &[String], _out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let _ = writeln!(err, "og archaic genomes: not ported yet");
    2
}
