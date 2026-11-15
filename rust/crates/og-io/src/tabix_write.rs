//! `og vcf index`: a `.tbi` index for a BGZF VCF, which step 38 writes
//! today with `pysam.tabix_index`.
//!
//! Not ported yet: until its unit of the 1.0.0 plan lands, the verb says so
//! on stderr and exits 2, before reading any argument.

use std::io::Write;

/// `og vcf index ARGS`: not ported yet, so one line on stderr and exit 2.
pub fn cli(_args: &[String], _out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let _ = writeln!(err, "og vcf index: not ported yet");
    2
}
