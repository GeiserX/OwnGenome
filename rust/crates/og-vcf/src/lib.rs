//! VCF read in place.
//!
//! A record is a borrowed line: [`Record::parse`] finds the tabs once and
//! every column after that is a slice of the line, so reading a record
//! copies nothing. INFO is read as an ordered list of items and written back
//! byte for byte ([`write_info`]); a tool that drops or adds items keeps the
//! others in their order. [`Gt`] reads a genotype, [`CsqLayout`] finds the
//! fields of VEP's CSQ annotation, and [`Header`] keeps the header lines as
//! they were, so a filter passes them through untouched.
//!
//! ```
//! use og_vcf::{Record, write_info};
//!
//! let line = b"chr1\t100\t.\tA\tG\t50\tPASS\tDP=9;DB;AF=0.5\tGT:DP\t0|1:9";
//! let rec = Record::parse(line).unwrap();
//! assert_eq!(rec.pos(), Ok(100));
//! // Drop DB; the other items keep their order.
//! let mut info = Vec::new();
//! write_info(rec.info_items().filter(|(k, _)| *k != b"DB"), &mut info);
//! let mut out = Vec::new();
//! rec.write_with_info(&info, &mut out);
//! assert_eq!(out, b"chr1\t100\t.\tA\tG\t50\tPASS\tDP=9;AF=0.5\tGT:DP\t0|1:9");
//! ```

mod csq;
mod gt;
mod header;
mod record;

pub use csq::{CsqLayout, csq_field};
pub use gt::Gt;
pub use header::{Header, Lines};
pub use record::{Error, InfoItems, Record, write_info};
