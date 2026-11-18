//! Step 06's known deep intronic sites, reported on every run: they lie far
//! from any exon, so most filters miss them.

use crate::escape::e;
use crate::model::{Line, Malformed};
use crate::notes::VIEWER_NOTE;
use crate::sections::clinvar::{D, VIEWER_TH, ev, viewer_td};

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

/// The HTML table of the sites, with the Viewer note.
pub fn table(sites: &D<'_>) -> Result<Vec<Line>, Malformed> {
    let mut out = vec![
        Line::Raw("    <table>".to_owned()),
        Line::Raw(format!(
            "      <tr><th>Gene</th><th>Variant</th><th>Position</th><th>Call</th>{VIEWER_TH}</tr>"
        )),
    ];
    for x in sites.items()? {
        let (chrom, pos) = (x.at("chrom")?, x.at("pos")?);
        out.push(Line::Raw(format!(
            "      <tr><td>{}</td><td>{}</td><td>{}:{} {}&gt;{}</td><td>{}</td>{}</tr>",
            ev(&x.at("gene")?)?,
            ev(&x.at("hgvs")?)?,
            ev(&chrom)?,
            e(&pos.text()?),
            ev(&x.at("ref")?)?,
            ev(&x.at("alt")?)?,
            e(&call(&x)?),
            viewer_td(chrom.v, pos.v, None)?
        )));
    }
    out.push(Line::Raw("    </table>".to_owned()));
    out.push(Line::Raw(VIEWER_NOTE.to_owned()));
    Ok(out)
}
