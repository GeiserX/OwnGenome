//! Step 06's cDNA marks, shown in the ClinVar block of both reports: the
//! read check of each hit (alt reads clipped at an exon edge), the hits at a
//! known lab-construct site, and the scan of every exon for reads that jump
//! from one exon to the next.

use crate::escape::e;
use crate::model::{Line, Malformed};
use crate::sections::clinvar::{D, VIEWER_TH, viewer_td};

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

/// The scan's lines, the same in both reports; none when it did not run.
pub fn scan_lines(d: &D<'_>) -> Result<Vec<String>, Malformed> {
    let Some(scan) = d.get("cdna_scan")?.filter(|x| !x.is_none()) else {
        return Ok(Vec::new());
    };
    let loci = scan.at("loci")?;
    let note = scan.at("note")?;
    let note_text = note.string()?;
    if note_text.starts_with("not run") {
        return Ok(vec![format!("cDNA scan of every exon: {note_text}")]);
    }
    // The note says what was read (exon regions, contigs): a BAM short of
    // contigs is no whole-exome negative.
    let loci = loci.items()?;
    let mut first = format!(
        "cDNA scan of every exon: {} loci with reads across two or more exon junctions",
        loci.len()
    );
    if note.truthy() {
        first += &format!(" ({note_text})");
    }
    let mut out = vec![first];
    for x in loci.iter().take(20) {
        out.push(format!(
            "  {}:{}-{}: {} junctions, {} reads",
            x.s("chrom")?,
            x.s("start")?,
            x.s("end")?,
            x.s("junctions")?,
            x.s("reads")?
        ));
    }
    Ok(out)
}

/// The HTML table of the scan's first 20 loci, with the Viewer note.
pub fn scan_table(scan: &D<'_>) -> Result<Vec<Line>, Malformed> {
    let mut out = vec![
        Line::Raw("    <table>".to_owned()),
        Line::Raw(format!(
            "      <tr><th>Locus</th><th>Junctions</th><th>Reads</th>{VIEWER_TH}</tr>"
        )),
    ];
    for x in scan.at("loci")?.items()?.iter().take(20) {
        let (chrom, start, end) = (x.at("chrom")?, x.at("start")?, x.at("end")?);
        out.push(Line::Raw(format!(
            "      <tr><td>{}:{}-{}</td><td>{}</td><td>{}</td>{}</tr>",
            e(chrom.string()?),
            e(&start.text()?),
            e(&end.text()?),
            e(&x.s("junctions")?),
            e(&x.s("reads")?),
            viewer_td(chrom.v, start.v, Some(end.v))?
        )));
    }
    out.push(Line::Raw("    </table>".to_owned()));
    out.push(Line::Raw(crate::notes::VIEWER_NOTE.to_owned()));
    Ok(out)
}
