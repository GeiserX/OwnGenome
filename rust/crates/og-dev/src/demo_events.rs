//! `og dev demo-events`: the replayed event log of the live demo, which
//! `tests/demo/make_demo_events.py` writes today.
//!
//! Not ported yet: until its unit of the 1.0.0 plan lands, the verb says so
//! on stderr and exits 2, before reading any argument.

use std::io::Write;

/// `og dev demo-events ARGS`: not ported yet, so one line on stderr and exit 2.
pub fn cli(_args: &[String], _out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let _ = writeln!(err, "og dev demo-events: not ported yet");
    2
}
