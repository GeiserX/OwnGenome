//! Step 16b's gene callability: how many of the hard-to-read genes of
//! `assets/callability_genes.tsv` the reads covered, and which were not
//! read, apart from the ones a step of their own calls. A share prints at one
//! decimal, so 0.8988 reads 89.9%, never the 90% a gene needs; a gene whose
//! body passed names the coding exon that did not; an owned gene says when
//! its own step did not run.

use pycompat::fmt::fixed;

use crate::model::{Card, Colour, HtmlPart, Inline, Line, Malformed, State, Summary};
use crate::notes::CALLABILITY_NOTE;
use crate::sections::clinvar::D;
use crate::{html, text};

const TITLE: &str = "Genes This Genome Could Not Read Reliably";

/// `share` as a percent at one decimal. A share under `min` never prints at
/// or above it: 0.8999 is 89.9%, not 90.0%. The bar is compared as printed,
/// since `100.0 * 0.55` is `55.00000000000001` and `55.0` would pass under it.
pub fn percent(share: f64, min: f64) -> String {
    let p = 100.0 * share;
    let t = fixed(p, 1);
    let bar = fixed(100.0 * min, 1).parse::<f64>().unwrap_or(100.0 * min);
    if share < min && t.parse::<f64>().is_ok_and(|x| x >= bar) {
        let down = (p * 10.0).floor() / 10.0;
        return fixed(if down >= bar { bar - 0.1 } else { down }, 1);
    }
    t
}

/// `GENE (37.0%)`, `GENE (95.0%; exon GENE:exon3 not read, 41.0%)`, or
/// `GENE (not on this reference)`.
fn gene(g: &D<'_>, min: f64) -> Result<String, Malformed> {
    if g.at("status")?.is("not on this reference") {
        return Ok(format!("{} (not on this reference)", g.s("gene")?));
    }
    let mut out = format!("{} ({}%", g.s("gene")?, percent(g.at("share")?.num()?, min));
    if let Some(w) = g.opt("worst_exon")? {
        out += &format!(
            "; exon {} not read, {}%",
            w.s("exon")?,
            percent(w.at("share")?.num()?, min)
        );
    }
    out.push(')');
    Ok(out)
}

/// An owned gene and its owner, with the owner's step when it did not run
/// (`owner_status`, from the summaries that carry it).
fn owned_gene(g: &D<'_>, min: f64) -> Result<String, Malformed> {
    let mut out = format!("{} by {}", gene(g, min)?, g.s("owner")?);
    if let (Some(step), Some(status)) = (g.opt("owner_step")?, g.opt("owner_status")?)
        && !status.is("ran")
    {
        out += &format!(", owner step {} {}", step.text()?, status.text()?);
    }
    Ok(out)
}

/// Whether an owned gene's step ran: true when the summary does not say.
fn owner_ran(g: &D<'_>) -> Result<bool, Malformed> {
    Ok(g.opt("owner_status")?.is_none_or(|s| s.is("ran")))
}

/// The share's caveat: it counts introns, so the exon check says what it
/// adds, or that this run had none.
fn introns_line(d: &D<'_>) -> Result<String, Malformed> {
    let head = "The share counts every base of a gene, introns included, so a gene can pass while a \
                coding exon was not read";
    Ok(match d.get("exons_checked")? {
        Some(n) => format!(
            "{head}: the {} coding exons checked must pass too, and a gene with one that did not \
             is listed as not read with that exon named.",
            n.text()?
        ),
        None => {
            format!("{head}: this run has no exon table from step 16b, so no exon was checked.")
        }
    })
}

/// The lines of the text block and the card.
pub fn lines(d: &D<'_>) -> Result<Vec<String>, Malformed> {
    let min = d.at("min_share")?.num()?;
    let mut out = vec![format!(
        "{} of the {} hard-to-read genes of assets/callability_genes.tsv were read: {}% or more of their \
         bases had {} or more reads of mapping quality {} or more.",
        d.s("read")?,
        d.s("catalog_genes")?,
        fixed(100.0 * min, 0),
        d.s("min_depth")?,
        d.s("min_mapq")?
    )];
    out.push(introns_line(d)?);
    let nr = d.at("not_read")?.items()?;
    let listed = if nr.is_empty() {
        "none".to_owned()
    } else {
        nr.iter()
            .map(|g| gene(g, min))
            .collect::<Result<Vec<_>, _>>()?
            .join(", ")
    };
    out.push(format!(
        "Not read, and no other step calls them ({}): {listed}",
        nr.len()
    ));
    let (mut ran, mut not_ran) = (Vec::new(), Vec::new());
    for g in d.at("not_read_owned")?.items()? {
        let line = owned_gene(&g, min)?;
        if owner_ran(&g)? {
            ran.push(line);
        } else {
            not_ran.push(line);
        }
    }
    if !not_ran.is_empty() {
        out.push(format!(
            "Not read here, and the step that calls them gave no result in this run ({}): {}",
            not_ran.len(),
            not_ran.join(", ")
        ));
    }
    if !ran.is_empty() {
        out.push(format!(
            "Not read here, called by a step of their own ({}): {}",
            ran.len(),
            ran.join(", ")
        ));
    }
    Ok(out)
}

/// The owned genes whose step gave no result in this run, when the summary
/// says which ran.
fn owner_not_run(d: &D<'_>) -> Result<Option<usize>, Malformed> {
    let owned = d.at("not_read_owned")?.items()?;
    if !owned.iter().any(|g| g.has("owner_status")) {
        return Ok(None);
    }
    let mut n = 0;
    for g in &owned {
        if !owner_ran(g)? {
            n += 1;
        }
    }
    Ok(Some(n))
}

/// The text report's block.
pub fn text(s: &Summary) -> Result<Vec<String>, Malformed> {
    let Some(head) = text::head(s, "callability", Some(TITLE))? else {
        return Ok(Vec::new());
    };
    let mut w = head.lines;
    if head.ran {
        let d = D::of(&s.section("callability"));
        w.extend(lines(&d)?.iter().map(|l| format!("  {l}")));
        w.push(format!("  {CALLABILITY_NOTE}"));
        w.push(String::new());
    }
    Ok(w)
}

fn count_badge(n: usize) -> Inline {
    Inline::Badge(
        n.to_string(),
        if n > 0 { Colour::Red } else { Colour::Green },
    )
}

/// The HTML report's card, when the step ran or failed.
pub fn html(s: &Summary) -> Result<HtmlPart, Malformed> {
    let sec = s.section("callability");
    if sec.state == State::Missing {
        return Ok(HtmlPart::default());
    }
    let d = D::of(&sec);
    let mut body = Vec::new();
    if sec.state.ran() {
        let not_read = d.at("not_read")?.len()?;
        body.push(Line::Stat(
            "Genes read".to_owned(),
            Inline::Text(format!("{} of {}", d.s("read")?, d.s("catalog_genes")?)),
        ));
        body.push(Line::Stat(
            "Not read, no other caller".to_owned(),
            count_badge(not_read),
        ));
        if let Some(n) = owner_not_run(&d)? {
            body.push(Line::Stat(
                "Not read, their own step gave no result".to_owned(),
                count_badge(n),
            ));
        }
        body.extend(lines(&d)?.into_iter().map(Line::Note));
    }
    let card = html::card(
        &sec,
        Card {
            title: TITLE.to_owned(),
            caveats: vec![CALLABILITY_NOTE],
            body,
            ..Card::default()
        },
    )?;
    Ok(HtmlPart {
        cards: vec![card],
        files: Vec::new(),
    })
}
