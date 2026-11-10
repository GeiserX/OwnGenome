//! Half-open intervals on named sequences.
//!
//! - [`bed`]: BED rows (0-based, end exclusive) read as OwnGenome's Python
//!   read them: whitespace-split, `#`, `track` and `browser` lines skipped,
//!   rows of fewer than three fields skipped; each sequence's rows sorted and
//!   merged into disjoint intervals.
//! - [`mask`]: a FASTA copied byte for byte with every base inside a row set
//!   to `N`, in one pass over the file.

pub mod bed;
pub mod mask;

pub use bed::{BedError, Rows};
pub use mask::{MaskError, Summary, mask_fasta};
