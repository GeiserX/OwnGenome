//! The sample's own carrier list (no partner): the ClinVar P/LP hits in
//! the genes of `assets/carrier_genes.tsv`, the rows with no unchanged copy
//! (`hom_or_hemizygous`, never counted as carriers), the hypomorph and
//! pseudodeficiency alleles set aside, and the two-hit genes
//! ([`super::two_hits`]).

use crate::model::{Card, HtmlPart, Line, Malformed, State, Summary};
use crate::notes::CARRIER_NOTE;
use crate::sections::clinvar::D;
use crate::sections::two_hits;
use crate::{html, text};

const TITLE: &str = "Carrier List (this sample)";

fn hemi(x: &D<'_>) -> Result<String, Malformed> {
    Ok(match x.opt("hemizygous")? {
        Some(h) => format!(" ({})", h.text()?),
        None => String::new(),
    })
}

/// The list's lines, the same in both reports.
pub fn lines(c: &D<'_>) -> Result<Vec<String>, Malformed> {
    if c.flag("unavailable")? {
        return Ok(vec![format!(
            "Carrier list: not made, {}",
            c.s("unavailable")?
        )]);
    }
    let carriers = c.at("carriers")?.items()?;
    // Homozygous or hemizygous: no unchanged copy, so listed apart and never
    // counted as a carrier. A summary written before the list has none.
    let both = match c.opt("hom_or_hemizygous")? {
        Some(b) => b.items()?,
        None => Vec::new(),
    };
    let set_aside = c.at("set_aside")?;
    let mut head = format!(
        "Carrier list ({} genes): {} finding(s)",
        c.s("genes_on_list")?,
        carriers.len()
    );
    if !both.is_empty() {
        head += &format!(", {} with no unchanged copy", both.len());
    }
    if set_aside.truthy() {
        head += &format!(", {} set aside", set_aside.len()?);
    }
    let mut out = vec![head];
    for x in &carriers {
        out.push(format!(
            "{} ({}; {}): {} {}{}, ClinVar {}",
            x.s("gene")?,
            x.s("condition")?,
            x.s("inheritance")?,
            x.s("variant")?,
            x.s("genotype")?,
            hemi(x)?,
            x.s("significance")?
        ));
    }
    for x in &both {
        let copies = if x.at("genotype")?.is("hom") {
            "two changed copies"
        } else {
            "the only copy changed"
        };
        out.push(format!(
            "{} ({}; {}): {} {}{}, ClinVar {}: {copies}, not the pattern of a carrier",
            x.s("gene")?,
            x.s("condition")?,
            x.s("inheritance")?,
            x.s("variant")?,
            x.s("genotype")?,
            hemi(x)?,
            x.s("significance")?
        ));
    }
    for x in set_aside.items()? {
        out.push(format!(
            "Set aside: {} {} {}{}: {} is a {} allele, not a carrier finding on its own",
            x.s("gene")?,
            x.s("variant")?,
            x.s("genotype")?,
            hemi(&x)?,
            x.s("name")?,
            x.s("kind")?
        ));
    }
    out.extend(two_hits::lines(c.get("two_hits")?.as_ref())?);
    Ok(out)
}

/// The text report's block.
pub fn text(s: &Summary) -> Result<Vec<String>, Malformed> {
    let Some(head) = text::head(s, "carrier", Some(TITLE))? else {
        return Ok(Vec::new());
    };
    let mut w = head.lines;
    if head.ran {
        let c = D::of(&s.section("carrier"));
        w.extend(lines(&c)?.iter().map(|l| format!("  {l}")));
        w.push(format!("  {CARRIER_NOTE}"));
        w.push(String::new());
    }
    Ok(w)
}

/// The HTML report's card, when the list was made or its step failed.
pub fn html(s: &Summary) -> Result<HtmlPart, Malformed> {
    let sec = s.section("carrier");
    if sec.state == State::Missing {
        return Ok(HtmlPart::default());
    }
    let body = if sec.state.ran() {
        lines(&D::of(&sec))?.into_iter().map(Line::Para).collect()
    } else {
        Vec::new()
    };
    let card = html::card(
        &sec,
        Card {
            title: TITLE.to_owned(),
            full: true,
            caveats: vec![CARRIER_NOTE],
            body,
            ..Card::default()
        },
    )?;
    Ok(HtmlPart {
        cards: vec![card],
        files: Vec::new(),
    })
}
