//! `og longread clair3-model`: the Clair3 model that matches the basecaller
//! of a long-read sample.
//!
//! Not ported yet: until its unit of the 1.0.0 plan lands, the verb says so
//! on stderr and exits 2, before reading any argument.

use std::io::Write;

/// `og longread clair3-model ARGS`: not ported yet, so one line on stderr and exit 2.
pub fn clair3_model_cli(_args: &[String], _out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let _ = writeln!(err, "og longread clair3-model: not ported yet");
    2
}
