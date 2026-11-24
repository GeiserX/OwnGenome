//! Step 05's AnnotSV class 4 and 5 SVs, with how many callers saw each
//! (step 22's consensus), the exon-edge and benign-set checks, and the
//! calls step 05's rules moved to class 3 with the reason.
//!
//! It is one block of the structural-variants section, beside the SV
//! callers' counts (`sv_callers`): [`text`] and [`html`] write the AnnotSV
//! block and card, and the section's writer puts them after the callers'.

use pycompat::fmt::{general, ljust};

use crate::model::{
    Card, Cell, Colour, HtmlPart, Inline, Layout, Line, Malformed, State, Summary, Table,
};
use crate::notes::ANNOTSV_NOTE;
use crate::sections::clinvar::{D, locus_text, viewer_of};
use crate::{html, text};

const TITLE: &str = "Structural Variants, Class 4 and 5 (AnnotSV)";

/// Why an SV's caller support is not known: step 22's consensus ran
/// (`consensus`, the section's flag) and cannot tell, or it did not run.
pub fn unknown_reason(consensus: bool) -> &'static str {
    if consensus {
        "step 22's consensus cannot tell"
    } else {
        "no step 22 consensus"
    }
}

/// How many callers saw the SV; `consensus` is the section's flag.
pub fn support_text(x: &D<'_>, consensus: bool) -> Result<String, Malformed> {
    let s = x.at("support")?;
    if s.is_none() {
        return Ok(format!("not known ({})", unknown_reason(consensus)));
    }
    if s.is_num(1) {
        return Ok(format!("single caller ({})", x.s("caller")?));
    }
    let callers = x.at("support_callers")?.items()?;
    let names: Vec<&str> = callers.iter().map(D::string).collect::<Result<_, _>>()?;
    Ok(format!("{} callers ({})", s.text()?, names.join(", ")))
}

/// The first three genes, and how many more.
pub fn genes_text(x: &D<'_>) -> Result<String, Malformed> {
    let g = x.at("genes")?;
    if !g.truthy() {
        return Ok("no gene".to_owned());
    }
    let items = g.items()?;
    let first: Vec<&str> = items
        .iter()
        .take(3)
        .map(D::string)
        .collect::<Result<_, _>>()?;
    let mut out = first.join(", ");
    if items.len() > 3 {
        out += &format!(" and {} more", items.len() - 3);
    }
    Ok(out)
}

/// What step 05 writes in `Population_source` when AnnotSV gave no
/// benign-set AF: the frequency is unknown, which is not the same as rare.
pub const NO_FREQUENCY: &str = "no frequency data";

/// What the reports say of an SV whose benign-set AF AnnotSV left blank.
pub const NO_FREQUENCY_TEXT: &str =
    "benign-set AF: no frequency data (unknown, not the same as rare)";

/// Where the benign-set AF comes from, said beside it on a row that also
/// carries step 05's per-record evidence.
const AGGREGATE: &str =
    "the highest AF of the benign records listed, larger SVs around it included";

/// The caveat under the AnnotSV block when a row carries step 05's
/// per-record benign evidence or has no frequency data: the benign-set AF
/// (and "common", read from it) is AnnotSV's aggregate, the same-size
/// record is the one-by-one check the two-hit and secondary-findings rules
/// use. The live card shows it too, without the report's common count, so
/// it names no count.
pub const BENIGN_EVIDENCE_NOTE: &str = "The benign-set AF is AnnotSV's aggregate: the highest AF of \
    every benign record of the call's type it lists, larger SVs around the call included, so it \
    is not this SV's own frequency, and \"common\" is read from it. Step 05 also checks the records \
    one by one: a same-size record at 1% or more is a record of the same type that covers at \
    least half of the call, and the call at least half of it, with an AF of 1% or more \
    (unknown when the evidence cannot tell); only those calls are set aside by the two-hit and \
    secondary-findings rules. No frequency data means AnnotSV gave no AF: the frequency is \
    unknown, which is not the same as rare (docs/05-annotsv.md).";

/// Step 05's `Benign_match` in words.
pub fn benign_match_text(v: &str) -> &str {
    match v {
        "same" => "one of about the same size",
        "contained" => "only larger SVs around it",
        "none" => "none listed",
        "." => "no benign set for this SV type",
        other => other,
    }
}

/// Whether the row carries step 05's per-record benign evidence
/// (`benign_match` or `clinical_common`, from tables that have the columns).
pub fn has_evidence(x: &D<'_>) -> Result<bool, Malformed> {
    Ok(x.opt("benign_match")?.is_some() || x.opt("clinical_common")?.is_some())
}

/// Whether AnnotSV gave the row no benign-set AF, as step 05 marks it.
pub fn no_frequency(x: &D<'_>) -> Result<bool, Malformed> {
    Ok(x.get("population_af")?.is_none_or(|v| v.is_none())
        && x.get("population_source")?
            .is_some_and(|v| v.is(NO_FREQUENCY)))
}

/// The exon-edge, benign-set and rule flags of one class 4-5 SV.
pub fn flags(x: &D<'_>) -> Result<Vec<String>, Malformed> {
    let mut out = Vec::new();
    let edge = x.at("exon_edge")?;
    if edge.is("yes") {
        out.push(format!(
            "on exon edges: {}",
            x.at("exon_edge_detail")?.string()?
        ));
    } else if edge.is("not checked") {
        out.push("exon edges not checked".to_owned());
    }
    let evidence = has_evidence(x)?;
    let af = x.at("population_af")?;
    if !af.is_none() {
        let mut f = format!("benign-set AF {}", general(af.num()?, 6));
        let src = x.at("population_source")?;
        match (src.truthy(), evidence) {
            (true, false) => f += &format!(" ({})", src.text()?),
            (true, true) => f += &format!(" ({}; {AGGREGATE})", src.text()?),
            (false, true) => f += &format!(" ({AGGREGATE})"),
            (false, false) => {}
        }
        if x.at("common")?.truthy() {
            f += ", common";
        }
        out.push(f);
    } else if no_frequency(x)? {
        out.push(NO_FREQUENCY_TEXT.to_owned());
    }
    let mut parts = Vec::new();
    if let Some(m) = x.opt("benign_match")? {
        parts.push(format!(
            "benign records: {}",
            benign_match_text(m.string()?)
        ));
    }
    if let Some(c) = x.opt("clinical_common")? {
        parts.push(format!("same-size record at 1% or more: {}", c.text()?));
    }
    if !parts.is_empty() {
        out.push(parts.join(", "));
    }
    if let Some(detail) = x.opt("rule_detail")? {
        let class = x.s("acmg_class")?;
        if class == "4" || class == "5" {
            out.push(format!("rule {}: {}", x.s("rule")?, detail.text()?));
        }
    }
    Ok(out)
}

/// `v or default`, as text.
fn or_text(x: &D<'_>, key: &str, default: &str) -> Result<String, Malformed> {
    let v = x.at(key)?;
    if v.truthy() {
        v.text()
    } else {
        Ok(default.to_owned())
    }
}

/// One SV step 05's rules moved to class 3: where, from which class, and
/// why.
pub fn moved_text(x: &D<'_>) -> Result<String, Malformed> {
    Ok(format!(
        "{} {}  {}  AnnotSV class {} -> {}  ({}: {})",
        x.s("type")?,
        x.s("sv")?,
        genes_text(x)?,
        x.s("annotsv_class")?,
        x.s("acmg_class")?,
        or_text(x, "rule", "rule not named")?,
        or_text(x, "rule_detail", "no detail")?
    ))
}

/// What the count of step 05's `clinical_common` is called, in both
/// reports and on the live card.
pub const SAME_SIZE_LABEL: &str = "Calls with a confirmed same-size benign record at 1% or more";

/// The class 4-5 rows that carry step 05's `clinical_common`: how many read
/// yes (confirmed), then `, not assessed: N` for the rows the evidence
/// could not tell, so a 0 never reads as "none" when it is "not known".
/// `None` when no row carries the key (a table from before the column), so
/// older reports keep their lines.
pub fn same_size_common(d: &D<'_>) -> Result<Option<String>, Malformed> {
    let Some(rows) = d.opt("rows")? else {
        return Ok(None);
    };
    let mut any = false;
    let (mut yes, mut unknown) = (0, 0);
    for x in rows.items()? {
        if let Some(c) = x.opt("clinical_common")? {
            any = true;
            if c.is("yes") {
                yes += 1;
            } else if !c.is("no") {
                unknown += 1;
            }
        }
    }
    Ok(any.then(|| {
        if unknown > 0 {
            format!("{yes}, not assessed: {unknown}")
        } else {
            yes.to_string()
        }
    }))
}

/// Whether the block shows [`BENIGN_EVIDENCE_NOTE`]: a class 4-5 row has
/// per-record evidence or no frequency data.
pub fn shows_evidence_note(d: &D<'_>) -> Result<bool, Malformed> {
    let Some(rows) = d.opt("rows")? else {
        return Ok(false);
    };
    for x in rows.items()? {
        if has_evidence(&x)? || no_frequency(&x)? {
            return Ok(true);
        }
    }
    Ok(false)
}

/// The class 4-5 lines, the same in both reports.
pub fn lines(d: &D<'_>) -> Result<Vec<String>, Malformed> {
    let count = d.at("count")?;
    let mut head = format!("Class 4-5 SVs: {}", count.text()?);
    if count.truthy() {
        let consensus = d.at("consensus")?.truthy();
        let single = if consensus {
            // The consensus ran, but a row it cannot tell has no support:
            // count those apart, so "single caller: 0" never reads as
            // "every SV had two callers". A row without the key counts too.
            let mut unknown = 0;
            if let Some(rows) = d.opt("rows")? {
                for x in rows.items()? {
                    if x.get("support")?.is_none_or(|v| v.is_none()) {
                        unknown += 1;
                    }
                }
            }
            let single = d.s("single_caller")?;
            if unknown > 0 {
                format!("{single}, not known: {unknown}")
            } else {
                single
            }
        } else {
            "not known".to_owned()
        };
        let edge = if d.at("exon_edge_checked")?.truthy() {
            d.s("exon_edge")?
        } else {
            "not checked".to_owned()
        };
        head += &format!(
            " (single caller: {single}; deletions on exon edges: {edge}; common: {})",
            d.s("common")?
        );
    }
    let mut out = vec![head];
    if let Some(n) = same_size_common(d)? {
        out.push(format!(
            "{SAME_SIZE_LABEL}: {n} (the common count above uses the benign-set AF, AnnotSV's \
             aggregate)"
        ));
    }
    if let Some(down) = d.opt("downgraded")? {
        out.push(format!(
            "Moved to class 3 by step 05's rules, not counted above: {}",
            down.len()?
        ));
    }
    if count.truthy() && !d.at("consensus")?.truthy() {
        out.push(match d.opt("consensus_note")? {
            Some(n) => n.text()?,
            None => "Caller support: step 22's consensus did not run, so how many callers saw each SV is \
                     not known."
                .to_owned(),
        });
    }
    Ok(out)
}

/// The AnnotSV block of the text report.
pub fn text(s: &Summary) -> Result<Vec<String>, Malformed> {
    let Some(head) = text::head(s, "annotsv", Some(TITLE))? else {
        return Ok(Vec::new());
    };
    let mut w = head.lines;
    if !head.ran {
        return Ok(w);
    }
    let d = D::of(&s.section("annotsv"));
    w.extend(lines(&d)?.iter().map(|l| format!("  {l}")));
    let consensus = d.at("consensus")?.truthy();
    for x in d.at("rows")?.items()? {
        w.push(format!(
            "    {} {}  class {}  {}  support: {}",
            ljust(&x.s("type")?, 4),
            x.s("sv")?,
            x.s("acmg_class")?,
            genes_text(&x)?,
            support_text(&x, consensus)?
        ));
        for f in flags(&x)? {
            w.push(format!("         {f}"));
        }
    }
    let down = d.opt("downgraded")?;
    if let Some(down) = &down {
        w.push(
            "  Moved to class 3 by step 05's rules (from AnnotSV's class, with the reason):"
                .to_owned(),
        );
        for x in down.items()? {
            w.push(format!("    {}", moved_text(&x)?));
        }
    }
    if d.at("count")?.truthy() || down.is_some() {
        w.push(format!("  {ANNOTSV_NOTE}"));
    }
    if shows_evidence_note(&d)? {
        w.push(format!("  {BENIGN_EVIDENCE_NOTE}"));
    }
    w.push(String::new());
    Ok(w)
}

/// `type sv` of a row, as the tables show it.
fn sv_cell(x: &D<'_>) -> Result<String, Malformed> {
    Ok(format!(
        "{} {}",
        x.at("type")?.string()?,
        x.at("sv")?.string()?
    ))
}

/// The AnnotSV card of the HTML report, when the step ran or failed.
pub fn html(s: &Summary) -> Result<HtmlPart, Malformed> {
    let sec = s.section("annotsv");
    if sec.state == State::Missing {
        return Ok(HtmlPart::default());
    }
    let d = D::of(&sec);
    let done = if sec.state == State::Stale {
        Inline::Badge("Stale".to_owned(), Colour::Yellow)
    } else {
        Inline::Badge("Complete".to_owned(), Colour::Green)
    };
    let mut body = vec![Line::Stat("Status".to_owned(), done)];
    if d.truthy() {
        body.extend(lines(&d)?.into_iter().map(Line::Note));
    }
    let rows = d.opt("rows")?;
    if let Some(rows) = &rows {
        let consensus = d.at("consensus")?.truthy();
        let mut cells = Vec::new();
        for x in rows.items()? {
            let checks = flags(&x)?.join("; ");
            cells.push(vec![
                Cell::Text(sv_cell(&x)?),
                Cell::Text(x.s("acmg_class")?),
                Cell::Text(genes_text(&x)?),
                Cell::Text(support_text(&x, consensus)?),
                Cell::Text(if checks.is_empty() {
                    "none".to_owned()
                } else {
                    checks
                }),
                viewer_of(&locus_text(&x.at("sv")?)?)?,
            ]);
        }
        body.push(Line::Table(Table {
            head: ["SV", "Class", "Genes", "Callers", "Checks"]
                .map(str::to_owned)
                .to_vec(),
            rows: cells,
            viewer: true,
            layout: Layout::OneLine,
        }));
    }
    let down = d.opt("downgraded")?;
    if let Some(down) = &down {
        body.push(Line::Raw(
            "    <h3>Moved to class 3 by step 05's rules</h3>".to_owned(),
        ));
        let mut cells = Vec::new();
        for x in down.items()? {
            cells.push(vec![
                Cell::Text(sv_cell(&x)?),
                Cell::Text(x.s("annotsv_class")?),
                Cell::Text(x.s("acmg_class")?),
                Cell::Text(genes_text(&x)?),
                Cell::Text(or_text(&x, "rule", "not named")?),
                Cell::Text(or_text(&x, "rule_detail", "no detail")?),
                viewer_of(&locus_text(&x.at("sv")?)?)?,
            ]);
        }
        body.push(Line::Table(Table {
            head: ["SV", "AnnotSV class", "Class", "Genes", "Rule", "Why"]
                .map(str::to_owned)
                .to_vec(),
            rows: cells,
            viewer: true,
            layout: Layout::OneLine,
        }));
    }
    let mut caveats = Vec::new();
    if rows.is_some() || down.is_some() {
        body.push(Line::ViewerNote);
        caveats.push(ANNOTSV_NOTE);
    }
    if d.truthy() && shows_evidence_note(&d)? {
        caveats.push(BENIGN_EVIDENCE_NOTE);
    }
    let card = html::card(
        &sec,
        Card {
            title: TITLE.to_owned(),
            full: true,
            caveats,
            body,
            ..Card::default()
        },
    )?;
    Ok(HtmlPart {
        cards: vec![card],
        files: Vec::new(),
    })
}

#[cfg(test)]
mod tests {
    use super::html;
    use crate::escape::e;
    use crate::model::Summary;
    use crate::notes::{ANNOTSV_NOTE, VIEWER_NOTE};

    /// The structural-variants writer is not wired yet, so no report case
    /// reaches this card: this renders it through the section frame.
    fn page(data: &str) -> String {
        let s = Summary::parse(&format!(
            r#"{{"sections": {{"annotsv": {{"state": "ok", "data": {data}}}}}}}"#
        ))
        .unwrap();
        let part = html(&s).unwrap();
        crate::html::frame("structural-variants", &part, None, false).join("\n")
    }

    const ROWS: &str = r#""rows": [
        {"type": "DEL", "sv": "chr1:1000-5000", "acmg_class": "5", "genes": ["A<B", "C&D", "G1", "G2"],
         "support": 2, "support_callers": ["manta", "delly"], "exon_edge": "yes",
         "exon_edge_detail": "exon 2", "population_af": 0.001, "population_source": "gnomAD",
         "common": false, "rule": "R1", "rule_detail": "loss <of> function & more"},
        {"type": "INS", "sv": "chr3:7", "acmg_class": "4", "genes": [], "support": null,
         "exon_edge": "no", "population_af": null}]"#;

    const DOWN: &str = r#""downgraded": [
        {"type": "DUP", "sv": "chr2:10-90", "annotsv_class": "4", "acmg_class": "3", "genes": [],
         "rule": "", "rule_detail": ""},
        {"type": "DEL", "sv": "chr2:500-900", "annotsv_class": "5", "acmg_class": "3", "genes": ["X&Y"],
         "rule": "B2", "rule_detail": "common <AF> & benign"}]"#;

    const HEAD: &str = r#""count": 2, "consensus": true, "single_caller": 0, "exon_edge_checked": true,
        "exon_edge": 1, "common": 0"#;

    fn viewer(at: &str) -> String {
        format!(
            "<td><a class=\"locus\" href=\"http://127.0.0.1:8765/#{at}\" target=\"_blank\" \
             rel=\"noopener noreferrer\">{at}</a></td>"
        )
    }

    #[test]
    fn both_tables_escape_their_cells_and_the_viewer_line_follows_them_once() {
        let p = page(&format!("{{{HEAD}, {ROWS}, {DOWN}}}"));
        let rows = format!(
            "    <table><tr><th>SV</th><th>Class</th><th>Genes</th><th>Callers</th><th>Checks</th>\
             <th>Viewer</th></tr>\
             <tr><td>DEL chr1:1000-5000</td><td>5</td><td>A&lt;B, C&amp;D, G1 and 1 more</td>\
             <td>2 callers (manta, delly)</td><td>on exon edges: exon 2; benign-set AF 0.001 (gnomAD); \
             rule R1: loss &lt;of&gt; function &amp; more</td>{}</tr>\
             <tr><td>INS chr3:7</td><td>4</td><td>no gene</td><td>not known (step 22&#x27;s consensus cannot tell)</td>\
             <td>none</td>{}</tr></table>",
            viewer("chr1:1000-5000"),
            viewer("chr3:7")
        );
        let down = format!(
            "    <table><tr><th>SV</th><th>AnnotSV class</th><th>Class</th><th>Genes</th><th>Rule</th>\
             <th>Why</th><th>Viewer</th></tr>\
             <tr><td>DUP chr2:10-90</td><td>4</td><td>3</td><td>no gene</td><td>not named</td>\
             <td>no detail</td>{}</tr>\
             <tr><td>DEL chr2:500-900</td><td>5</td><td>3</td><td>X&amp;Y</td><td>B2</td>\
             <td>common &lt;AF&gt; &amp; benign</td>{}</tr></table>",
            viewer("chr2:10-90"),
            viewer("chr2:500-900")
        );
        let lines: Vec<&str> = p.lines().collect();
        let at = |l: &str| lines.iter().position(|x| *x == l);
        let (r, h, d, v) = (
            at(&rows).expect(&p),
            at("    <h3>Moved to class 3 by step 05's rules</h3>").expect(&p),
            at(&down).expect(&p),
            at(VIEWER_NOTE).expect(&p),
        );
        assert!(r < h && h + 1 == d && d + 1 == v, "{p}");
        assert_eq!(p.matches("<p class=\"viewer\">").count(), 1, "{p}");
        assert!(
            p.contains(&format!("    <p class=\"caveat\">{}</p>", e(ANNOTSV_NOTE))),
            "{p}"
        );
        // Red control: an unescaped gene would not match the row.
        assert!(!p.contains("A<B"));
    }

    #[test]
    fn every_benign_match_reads_in_words() {
        for (value, words) in [
            ("same", "one of about the same size"),
            ("contained", "only larger SVs around it"),
            ("none", "none listed"),
            (".", "no benign set for this SV type"),
            ("not checked", "not checked"),
        ] {
            assert_eq!(super::benign_match_text(value), words);
        }
    }

    #[test]
    fn one_table_alone_still_gets_the_viewer_line_and_none_gets_neither() {
        let p = page(&format!("{{{HEAD}, {DOWN}}}"));
        assert!(!p.contains("<th>Checks</th>"), "{p}");
        assert_eq!(p.matches("<p class=\"viewer\">").count(), 1, "{p}");
        let p = page(r#"{"count": 0}"#);
        assert!(!p.contains("<table"), "{p}");
        assert!(!p.contains("<p class=\"viewer\">"), "{p}");
        assert!(!p.contains("class=\"caveat\""), "{p}");
    }
}
