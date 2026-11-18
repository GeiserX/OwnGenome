//! Step 06's cDNA marks, shown in the ClinVar block of both reports: the
//! read check of each hit (alt reads clipped at an exon edge) and the hits at a
//! known lab-construct site. The scan of every exon for reads that jump from
//! one exon to the next is `sections::cdna_scan`.

use crate::model::Malformed;
use crate::sections::clinvar::D;

/// Step 06's read check of one ClinVar (or CPSR) hit. A hit at a known
/// construct site always shows its read counts when they were measured.
pub fn cdna_text(h: &D<'_>) -> Result<String, Malformed> {
    let c = h.at("cdna")?;
    let counts = format!(
        "{} of {} alt reads and {} of {} ref reads clipped at an exon edge",
        c.s("alt_clipped")?,
        c.s("alt_reads")?,
        c.s("ref_clipped")?,
        c.s("ref_reads")?
    );
    let site = match c.opt("construct")? {
        Some(x) => format!("known lab-construct site ({}): ", x.text()?),
        None => String::new(),
    };
    let check = c.at("check")?;
    if check.is("flagged") {
        return Ok(format!("{site}{} ({counts})", c.s("note")?));
    }
    if check.is("clear") {
        return Ok(if site.is_empty() {
            "clear".to_owned()
        } else {
            format!("{site}clear ({counts})")
        });
    }
    let note = c.at("note")?;
    let mut out = format!("{site}not checked");
    if note.truthy() {
        out += &format!(": {}", note.text()?);
    }
    let measured = c.get("alt_reads")?.is_some_and(|x| !x.is_none());
    if !site.is_empty() && measured {
        out += &format!(" ({counts})");
    }
    Ok(out)
}

/// The hits at a known lab-construct site below the first `shown`: both
/// reports print them too, so every hit the overview counts has a row.
pub fn construct_past_cap<'a>(hits: &[D<'a>], shown: usize) -> Result<Vec<D<'a>>, Malformed> {
    let mut out = Vec::new();
    for h in hits.iter().skip(shown) {
        let site = match h.opt("cdna")? {
            Some(c) => c.flag("construct")?,
            None => false,
        };
        if site {
            out.push(h.clone());
        }
    }
    Ok(out)
}

/// One such hit as a line.
pub fn construct_past_cap_text(h: &D<'_>) -> Result<String, Malformed> {
    Ok(format!(
        "{} {}:{} {}>{} {}: {}",
        h.s("gene")?,
        h.s("chrom")?,
        h.s("pos")?,
        h.s("ref")?,
        h.s("alt")?,
        h.s("genotype")?,
        cdna_text(h)?
    ))
}
