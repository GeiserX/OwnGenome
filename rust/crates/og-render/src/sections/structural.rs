//! The structural-variants section: the callers' counts (Manta, Delly,
//! CNVpytor and step 22's consensus, [`sv_callers`]), then AnnotSV's class 4
//! and 5 calls ([`annotsv`]).

use crate::model::{HtmlPart, Malformed, Summary};
use crate::sections::{annotsv, sv_callers};

/// The text report's blocks: each caller that ran, then AnnotSV's.
pub fn text(s: &Summary) -> Result<Vec<String>, Malformed> {
    let mut w = sv_callers::text_blocks(s)?;
    w.extend(annotsv::text(s)?);
    Ok(w)
}

/// The HTML report's cards: the callers' card, always shown, then AnnotSV's
/// when the step ran or was meant to.
pub fn html(s: &Summary) -> Result<HtmlPart, Malformed> {
    let mut cards = vec![sv_callers::card(s)?];
    cards.extend(annotsv::html(s)?.cards);
    Ok(HtmlPart {
        cards,
        files: Vec::new(),
    })
}
