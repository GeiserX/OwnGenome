//! Not ported yet: scripts/setup.sh's bash runs these flags. The slice
//! setup-port-base fills this module:
//! the default run (setup.sh <dir>), --no-clinical-data, --everything,
//! --refresh clinvar and --strchive-data (the STRchive catalog the default
//! run builds, with STR_EXTRA_LOCI). --og stays in setup.sh: it builds og.
//! It lists each flag it runs in [`PORTED`] and answers it in [`run`];
//! nothing else in the crate changes.

use super::Run;
use std::io::Write;

/// The keys this module runs.
pub const PORTED: &[&str] = &[];

/// Runs KEY when it is one of [`PORTED`]; its exit status.
pub fn run(_key: &str, _r: &Run, _out: &mut dyn Write, _err: &mut dyn Write) -> Option<u8> {
    None
}
