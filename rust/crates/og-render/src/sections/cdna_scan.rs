//! Step 06's genome-wide cDNA scan: the loci whose reads jump from one
//! exon to the next, a sign of cDNA or plasmid contamination or of a
//! processed copy of a gene.
//!
//! The scan is part of the ClinVar section's data (`clinvar.cdna_scan`),
//! and `sections::clinvar` shows it: [`text_lines`] in its text block,
//! [`stat`] on its screening card, and [`html_lines`] on its card of deep
//! intronic sites.

use pycompat::json::Value;

use crate::model::{
    Cell, Colour, Inline, Layout, Line, Malformed, Section, Table, need, py, truthy,
};
use crate::notes::CDNA_SCAN_NOTE;

const PATH: &str = "sections.clinvar.data.cdna_scan";

/// The most loci either report lists.
pub const SHOWN: usize = 20;

/// The scan, when step 06 wrote one: its note and loci.
fn scan<'a>(d: &Section<'a>) -> Result<Option<(String, &'a [Value])>, Malformed> {
    let Some(scan) = d.get("cdna_scan")?.filter(|v| !matches!(v, Value::Null)) else {
        return Ok(None);
    };
    let note = py(need(scan, "note", PATH)?)?;
    let Value::Array(loci) = need(scan, "loci", PATH)? else {
        return Err(Malformed(format!("{PATH}.loci is not a list")));
    };
    Ok(Some((note, loci)))
}

/// One locus as `chrom:start-end`.
fn locus(x: &Value) -> Result<String, Malformed> {
    let p = format!("{PATH}.loci[]");
    Ok(format!(
        "{}:{}-{}",
        py(need(x, "chrom", &p)?)?,
        py(need(x, "start", &p)?)?,
        py(need(x, "end", &p)?)?
    ))
}

/// The scan's lines as both reports print them: the count of loci, then
/// the first [`SHOWN`] loci; one line when the scan did not run, none
/// when step 06 wrote no scan.
pub fn lines(d: &Section<'_>) -> Result<Vec<String>, Malformed> {
    let Some((note, loci)) = scan(d)? else {
        return Ok(Vec::new());
    };
    if note.starts_with("not run") {
        return Ok(vec![format!("cDNA scan of every exon: {note}")]);
    }
    // The note says what was read: a BAM short of contigs is no
    // whole-exome negative.
    let mut out = vec![format!(
        "cDNA scan of every exon: {} loci with reads across two or more exon junctions{}",
        loci.len(),
        if note.is_empty() {
            String::new()
        } else {
            format!(" ({note})")
        }
    )];
    let p = format!("{PATH}.loci[]");
    for x in loci.iter().take(SHOWN) {
        out.push(format!(
            "  {}: {} junctions, {} reads",
            locus(x)?,
            py(need(x, "junctions", &p)?)?,
            py(need(x, "reads", &p)?)?
        ));
    }
    Ok(out)
}

/// The text report's lines of the scan in the ClinVar block, with the
/// caveat when a locus was found.
pub fn text_lines(d: &Section<'_>) -> Result<Vec<String>, Malformed> {
    let mut out: Vec<String> = lines(d)?.into_iter().map(|l| format!("  {l}")).collect();
    if scan(d)?.is_some_and(|(_, loci)| !loci.is_empty()) && d.get("cdna_scan")?.is_some_and(truthy)
    {
        out.push(format!("  {CDNA_SCAN_NOTE}"));
    }
    Ok(out)
}

/// The screening card's count of loci with spliced reads, red when there is
/// one; `None` when the scan did not run.
pub fn stat(d: &Section<'_>) -> Result<Option<Line>, Malformed> {
    if !d.get("cdna_scan")?.is_some_and(truthy) {
        return Ok(None);
    }
    let Some((note, loci)) = scan(d)? else {
        return Ok(None);
    };
    if note.starts_with("not run") {
        return Ok(None);
    }
    let n = loci.len();
    Ok(Some(Line::Stat(
        "Loci with spliced reads (cDNA scan)".to_owned(),
        Inline::Badge(
            n.to_string(),
            if n > 0 { Colour::Red } else { Colour::Gray },
        ),
    )))
}

/// The scan's part of the deep intronic card: the count line, then a table
/// of the loci with their Viewer links, and the caveat that goes under the
/// card's title when there is a table.
pub fn html_lines(d: &Section<'_>) -> Result<(Vec<Line>, Option<&'static str>), Malformed> {
    let lines = lines(d)?;
    let Some(first) = lines.first() else {
        return Ok((Vec::new(), None));
    };
    let mut out = vec![Line::Para(first.clone())];
    if lines.len() == 1 {
        return Ok((out, None));
    }
    let Some((_, loci)) = scan(d)? else {
        return Ok((out, None));
    };
    let p = format!("{PATH}.loci[]");
    let mut rows = Vec::new();
    for x in loci.iter().take(SHOWN) {
        rows.push(vec![
            Cell::Text(locus(x)?),
            Cell::Text(py(need(x, "junctions", &p)?)?),
            Cell::Text(py(need(x, "reads", &p)?)?),
            Cell::viewer(
                need(x, "chrom", &p)?,
                need(x, "start", &p)?,
                Some(need(x, "end", &p)?),
            )?,
        ]);
    }
    out.push(Line::Table(Table {
        head: ["Locus", "Junctions", "Reads"].map(str::to_owned).to_vec(),
        rows,
        viewer: true,
        layout: Layout::Indented,
    }));
    Ok((out, Some(CDNA_SCAN_NOTE)))
}
