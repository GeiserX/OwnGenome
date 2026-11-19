//! Step 06's known deep intronic sites, reported on every run: they lie far
//! from any exon, so most filters miss them.

use crate::model::{Cell, Layout, Line, Malformed, Table};
use crate::sections::clinvar::D;

/// One site's call: `het, GQ 40, DP 30`, `ref, GQ 50 (gVCF block)`,
/// `no call: why`.
pub fn call(x: &D<'_>) -> Result<String, Malformed> {
    let c = x.at("call")?;
    let detail = x.at("detail")?;
    if c.is("no call") {
        let mut out = "no call".to_owned();
        if detail.truthy() {
            out += &format!(": {}", detail.text()?);
        }
        return Ok(out);
    }
    let mut out = c.string()?.to_owned();
    if detail.truthy() {
        out += &format!(", {}", detail.text()?);
    }
    Ok(out)
}

/// The table of the sites, with their Viewer links, and the Viewer note.
pub fn table(sites: &D<'_>) -> Result<Vec<Line>, Malformed> {
    let st =
        |x: &D<'_>, k: &str| -> Result<String, Malformed> { Ok(x.at(k)?.string()?.to_owned()) };
    let mut rows = Vec::new();
    for x in sites.items()? {
        let (chrom, pos) = (x.at("chrom")?, x.at("pos")?);
        rows.push(vec![
            Cell::Text(st(&x, "gene")?),
            Cell::Text(st(&x, "hgvs")?),
            Cell::Text(format!(
                "{}:{} {}>{}",
                chrom.string()?,
                pos.text()?,
                st(&x, "ref")?,
                st(&x, "alt")?
            )),
            Cell::Text(call(&x)?),
            Cell::viewer(chrom.v, pos.v, None)?,
        ]);
    }
    Ok(vec![
        Line::Table(Table {
            head: ["Gene", "Variant", "Position", "Call"]
                .map(str::to_owned)
                .to_vec(),
            rows,
            viewer: true,
            layout: Layout::Indented,
        }),
        Line::ViewerNote,
    ])
}
