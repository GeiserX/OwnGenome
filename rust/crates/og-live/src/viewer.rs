//! `og view serve`: igv.js and one sample's reads and calls on 127.0.0.1,
//! which `assets/viewer/serve.py` serves today.
//!
//! Not ported yet: until its unit of the 1.0.0 plan lands, the verb says so
//! on stderr and exits 2, before reading any argument.

use std::io::Write;

/// `og view serve ARGS`: not ported yet, so one line on stderr and exit 2.
pub fn cli(_args: &[String], _out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let _ = writeln!(err, "og view serve: not ported yet");
    2
}
