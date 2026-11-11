//! `og pgx cpic-report` and `og pgx step36-outside`: the pipeline's one
//! reader of PharmCAT's `report.json` (`docs/27-cpic-lookup.md`).
//!
//! PharmCAT turns a VCF into star-allele diplotypes and phenotypes per
//! pharmacogene and matches them to CPIC, DPWG and FDA guidance. This crate
//! writes the per-gene phenotype table and the drug recommendations a person
//! reads, and carries the pipeline's own safety rules:
//!
//! - a report that parses to no gene is a failure, never an all-clear;
//! - several possible diplotypes with different phenotypes make a gene
//!   [`Status::Ambiguous`], with no drug guidance;
//! - an outside call of a step 36 gene (CYP2D6, HLA-A, HLA-B) counts only
//!   when step 36's consensus passed that same call;
//! - a gene left Indeterminate at gene level takes the phenotype its drug
//!   annotations agree on, with a note;
//! - CYP4F2 *3 and VKORC1 -1639 get a labelled single-position fallback when
//!   PharmCAT gives no single result;
//! - an opt-in medication list phenoconverts CYP2D6 by CPIC's activity-score
//!   rule and notes inhibitors of CYP2C19, CYP2C9, CYP3A4 and CYP3A5;
//! - a pathogenic ClinVar hit in a disease gene PharmCAT does not call adds a
//!   WARNING naming the drug.
//!
//! The text and the tables are what `bin/pgx_parse.py` wrote, byte for byte
//! on every input it read without a traceback (`tests/parity`).

mod cli;
pub mod guidance;
pub mod pyio;
pub mod report;
pub mod sections;
pub mod status;
pub mod tables;

pub use cli::{CpicArgs, cpic_report, main, step36_outside};
pub use report::{GeneCall, Parsed};
pub use status::{Status, is_non_normal, is_normal};
