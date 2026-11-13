//! The helpers the step scripts and the Nextflow modules call, one module
//! per tool.
//!
//! - [`constraint`]: `og table constraint-join`, gnomAD v4.1 gene constraint
//!   columns appended to a table.
//! - [`prs`]: `og prs format`, step 25's PGS Catalog scoring files written
//!   as the custom GRCh38 files pgsc_calc reads.
//! - [`consumer`]: `og consumer build-sites|export|check`, step 40's
//!   23andMe and AncestryDNA raw-data files.
//! - [`family`]: `og family plan|prepare|check|report|segments`, the family
//!   analysis of step 38 and the Nextflow FAMILY workflow.
//! - [`pgx_rules`]: `og pgx outside-calls|hla-concordance`, step 36's
//!   PharmCAT outside calls and step 08's T1K and OptiType table.
//! - [`pypgx_summary`]: `og pgx pypgx-summary`, step 32's summary of
//!   pypgx's results.
//! - [`paraphase`]: `og paraphase check-bam|table`, step 35b's read-group
//!   check and its two tables.
//! - [`sv`]: `og sv add-callers|postfilter`, step 05's SV calls into
//!   AnnotSV and the checks of its class 4 and 5 calls.

pub mod constraint;
pub mod consumer;
pub mod family;
pub mod paraphase;
pub mod pgx_rules;
pub mod prs;
pub mod pypgx_summary;
pub mod sv;
