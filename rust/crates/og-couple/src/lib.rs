//! Step 39, the couple carrier screen: `og couple screen` ([`screen`]) and
//! `og couple genes` ([`genes`]).
//!
//! For every two samples a pedigree marks as partners, the screen collects
//! each partner's carrier findings in the genes of the carrier list
//! (`assets/carrier_genes.tsv`): step 06's ClinVar hits, Delly deletions over
//! an exon, SMN1 copy number, CYP21A2 from depth, and the FMR1 and FXN
//! repeats. It then says, gene by gene, what the two findings mean for a
//! child, and writes one table and one report per couple.
//!
//! Each input is read once, line by line. The rules other steps own are
//! theirs, never copied: the read-depth rule for a long Delly deletion is
//! `og_steps::sv::depth_contradicts` (step 05), SMN1 at exon 7 is
//! `og_cn::smn::read` (step 35), the pedigree is `og_pedigree`'s (step 38),
//! and a ClinVar hit reads as `og_summary`'s `clinvar_row` reads it.

mod args;
mod cn;
mod files;
mod findings;
pub mod genes;
mod list;
mod prose;
pub mod screen;
mod xml;
