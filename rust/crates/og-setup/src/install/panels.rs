//! Not ported yet: scripts/setup.sh's bash runs these flags. The slice
//! setup-port-panels fills this module:
//! --ancestry-panel, --ancestry-fine-panel, --mocha-data, --ancient-panel,
//! --archaic and --consumer-sites.
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
