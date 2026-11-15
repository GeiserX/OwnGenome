//! `og couple screen`: step 39's couple carrier screen, which
//! `bin/couple_report.py screen` runs today.
//!
//! Not ported yet: until its unit of the 1.0.0 plan lands, the verb says so
//! on stderr and exits 2, before reading any argument.

use std::io::Write;

/// `og couple screen ARGS`: not ported yet, so one line on stderr and exit 2.
pub fn cli(_args: &[String], _out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let _ = writeln!(err, "og couple screen: not ported yet");
    2
}
