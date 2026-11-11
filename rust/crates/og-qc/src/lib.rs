//! Sample quality checks of step 33.
//!
//! - [`sex_loss`]: `og qc sex-loss`, the share of cells that lost the Y (a
//!   man) or one X (a woman), from mosdepth's windows and the VCF's allele
//!   fractions, added as `mosaic_*` rows to the step's key/value table.

pub mod sex_loss;
