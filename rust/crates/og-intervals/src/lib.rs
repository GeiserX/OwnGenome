//! Half-open intervals on named sequences.
//!
//! - [`bed`]: BED rows (0-based, end exclusive) read as OwnGenome's Python
//!   read them: whitespace-split, `#`, `track` and `browser` lines skipped,
//!   rows of fewer than three fields skipped; each sequence's rows sorted and
//!   merged into disjoint intervals.
//! - [`mask`]: a FASTA copied byte for byte with every base inside a row set
//!   to `N`, in one pass over the file.
//! - [`sets`]: sorted interval sets (union, intersection, lookups by binary
//!   search, BED files such as the PAR regions read into them) and an
//!   overlap index with reciprocal overlap.
//! - [`fasta`]: random access into a FASTA through its `.fai`.

pub mod bed;
pub mod fasta;
pub mod mask;
pub mod sets;

pub use bed::{BedError, Rows};
pub use mask::{MaskError, Summary, mask_fasta};
