//! The HTML report: one page with no script and no fetch, a link that skips
//! to the results, a table of contents with a state dot per section (folded
//! on a phone), the header with the research-only line under the title, the
//! chapters, what the report leaves out, the disclaimer and the run manifest.
//! Every table sits in its own scroll region, so the page fits a phone.
//!
//! Each section is a frame: an anchor, a line with the step, its state, the
//! newest file date and "New since" under the first card title, each card's
//! caveats under its title, the cards, then links to the step's own files.

use pycompat::json::Value;

use crate::chapters::{
    CHAPTERS, Failed, FrameState, LIMITS, Sec, TERMS, WriterFn, guarded, section_of, step_of,
    writer,
};
use crate::charts;
use crate::css::{CSS, READING};
use crate::escape::e;
use crate::model::{
    Card, Cell, Colour, HtmlPart, Inline, Layout, Line, Malformed, Section, State, Summary, Table,
    py,
};
use crate::notes::{NOTES, VIEWER_NOTE, stale_note};
use crate::sections::run;

/// The Content-Security-Policy: nothing is fetched, images only inline,
/// styles only from the page.
pub const CSP: &str = "default-src 'none'; img-src data:; style-src 'unsafe-inline'";

/// A link that leaves this computer, marked so a reader knows before
/// opening it.
pub const OFFSITE: &str = "opens a website";

/// The line under the report's title that says what the report is before
/// any result: the full disclaimer stays at the end.
pub const RESEARCH_ONLY: &str = "Research software, not a clinical report: confirm any finding in a \
certified laboratory before acting on it.";

/// The first thing a keyboard reaches: a link past the contents to the
/// results (`<main id="main">`), visible only while it has the focus.
pub const SKIP_LINK: &str = "<a class=\"skip\" href=\"#main\">Skip to results</a>";

const PRINT_HINT: &str = "<p class=\"print-hint\">To print one section, open it from the contents (or its \
\"link to this section\"), then print: the printout holds that section, this header, the list of what is not \
assessed and the disclaimer. A chapter prints the same way. In a browser without CSS :has() the whole report \
prints.</p>";

/// A coloured badge.
pub fn badge(text: &str, colour: Colour) -> String {
    format!(
        "<span class=\"badge badge-{}\">{}</span>",
        colour.class(),
        e(text)
    )
}

/// A link to an address off this computer, with its mark.
pub fn offsite_link(url: &str, text: &str, why: &str) -> String {
    format!(
        "<span class=\"offsite\"><a href=\"{}\">{}</a> <span class=\"offsite-note\">({})</span></span>",
        e(url),
        e(text),
        e(why)
    )
}

fn inline(v: &Inline) -> String {
    match v {
        Inline::Text(t) => e(t).into_owned(),
        Inline::Fixed(t) => (*t).to_owned(),
        Inline::Badge(t, c) => badge(t, *c),
        Inline::Quiet(t) => format!(
            "<span style=\"font-weight:normal;font-size:13px\">{}</span>",
            e(t)
        ),
        Inline::Break => "<br>".to_owned(),
        Inline::Small(t) => format!(
            "<br><span style=\"font-size:12px;color:#666\">{}</span>",
            e(t)
        ),
    }
}

fn line(out: &mut Vec<String>, l: &Line) {
    match l {
        Line::Stat(label, value) => out.push(format!(
            "    <div class=\"stat\"><span class=\"label\">{}</span><span class=\"value\">{}</span></div>",
            e(label),
            inline(value)
        )),
        Line::Table(t) => table(out, t),
        Line::Note(t) => out.push(format!("    <p class=\"note\">{}</p>", e(t))),
        Line::Para(t) => out.push(format!("    <p>{}</p>", e(t))),
        Line::List(items) => {
            out.push("    <ul class=\"plain\">".to_owned());
            out.extend(items.iter().map(|t| format!("      <li>{}</li>", e(t))));
            out.push("    </ul>".to_owned());
        }
        Line::Warn(t) => out.push(format!("    <div class=\"stale\">{}</div>", e(t))),
        Line::ViewerNote => out.push(VIEWER_NOTE.to_owned()),
        Line::Folded(summary, t) => {
            out.push(format!("    <details><summary>{}</summary>", e(summary)));
            table(out, t);
            if let Some(last) = out.last_mut() {
                last.push_str("</details>");
            }
        }
        Line::Raw(t) => out.push(t.clone()),
    }
}

/// The genome viewer on this computer, which a Viewer cell opens.
pub const VIEWER_URL: &str = "http://127.0.0.1:8765/";

fn cell(c: &Cell) -> String {
    match c {
        Cell::Text(t) => format!("<td>{}</td>", e(t)),
        Cell::Badge(t, colour) => format!("<td>{}</td>", badge(t, *colour)),
        Cell::Sub(t, sub) => format!("<td>{}{}</td>", e(t), inline(&Inline::Small(sub.clone()))),
        Cell::Parts(parts) => format!("<td>{}</td>", parts.iter().map(inline).collect::<String>()),
        Cell::Viewer(None) => "<td></td>".to_owned(),
        Cell::Viewer(Some(at)) => format!(
            "<td><a class=\"locus\" href=\"{VIEWER_URL}#{0}\" target=\"_blank\" rel=\"noopener noreferrer\">{0}</a></td>",
            e(at)
        ),
    }
}

fn table(out: &mut Vec<String>, t: &Table) {
    let mut head: String = t
        .head
        .iter()
        .map(|h| format!("<th>{}</th>", e(h)))
        .collect();
    if t.viewer {
        head.push_str("<th>Viewer</th>");
    }
    let rows = t
        .rows
        .iter()
        .map(|r| format!("<tr>{}</tr>", r.iter().map(cell).collect::<String>()));
    match t.layout {
        Layout::OneLine => out.push(format!(
            "    <table><tr>{head}</tr>{}</table>",
            rows.collect::<String>()
        )),
        Layout::Indented => {
            out.push("    <table>".to_owned());
            out.push(format!("      <tr>{head}</tr>"));
            out.extend(rows.map(|r| format!("      {r}")));
            out.push("    </table>".to_owned());
        }
        Layout::RowLines => {
            out.push("    <table>".to_owned());
            out.push(format!("      <tr>{head}</tr>"));
            out.push(rows.collect::<Vec<_>>().join("\n"));
            out.push("    </table>".to_owned());
        }
    }
}

/// Puts each table of `lines` in a scroll region, so a table wider than its
/// card scrolls inside it instead of widening the page on a phone. The
/// region takes the focus, so a keyboard can scroll it too. It is named by
/// the title of the card it sits in (the last card title line above it),
/// else by `label`.
fn scroll_tables(lines: &mut [String], label: &str) {
    let mut name = e(label).into_owned();
    for l in lines.iter_mut() {
        if let Some(title) = l
            .strip_prefix("    <h2>")
            .and_then(|t| t.strip_suffix("</h2>"))
        {
            name = title.to_owned();
        } else if l.contains("table>") {
            let open = format!(
                "<div class=\"tscroll\" role=\"region\" aria-label=\"{name}\" tabindex=\"0\"><table>"
            );
            *l = l
                .replace("<table>", &open)
                .replace("</table>", "</table></div>");
        }
    }
}

/// The card of one summary section: a section that did not run shows
/// "Not run", a failed one "Failed" and an unreadable one "Unreadable",
/// in place of its lines; a stale, failed or unreadable one says so in a box.
pub fn card(sec: &Section<'_>, mut c: Card) -> Result<Card, Malformed> {
    let note = stale_note(sec)?;
    if !note.is_empty() {
        c.stale.push(note);
    }
    let status = match sec.state {
        State::Missing => Some(("Not run", Colour::Gray)),
        State::Failed => Some(("Failed", Colour::Red)),
        State::Unreadable => Some(("Unreadable", Colour::Gray)),
        _ => None,
    };
    if let Some((text, colour)) = status {
        c.caveats.clear();
        c.body = vec![Line::Stat(
            "Status".to_owned(),
            Inline::Badge(text.to_owned(), colour),
        )];
    }
    Ok(c)
}

/// The step line under a section's first card title, with the section's
/// state badge.
fn meta(
    s: &Summary,
    sid: &str,
    keys: &[&str],
    state: FrameState,
    since: Option<&str>,
) -> Result<String, Malformed> {
    let secs: Vec<Section<'_>> = keys.iter().map(|k| s.section(k)).collect();
    let mut steps: Vec<&str> = Vec::new();
    for k in keys {
        let step = step_of(k);
        if !steps.contains(&step) {
            steps.push(step);
        }
    }
    let (label, colour) = state.badge();
    let mut dates = Vec::new();
    for x in &secs {
        if let Some(d) = x.file_date.filter(|d| crate::model::truthy(d)) {
            dates.push(py(d)?);
        }
    }
    let mut out = format!(
        "    <p class=\"sec-meta\">Step{} {} &middot; {}",
        if steps.len() > 1 { "s" } else { "" },
        e(&steps.join(", ")),
        badge(label, colour)
    );
    if let Some(newest) = dates.iter().max() {
        out += &format!(" &middot; file date {}", e(newest));
    }
    // The fault-new-since-html-badge control drops the badge.
    if let Some(since) = since.filter(|_| !cfg!(feature = "fault-new-since-html-badge")) {
        out += &format!(
            " &middot; {}",
            badge(&format!("New since {since}"), Colour::Yellow)
        );
    }
    out += &format!(" &middot; <a class=\"anchor\" href=\"#{sid}\">link to this section</a></p>");
    Ok(out)
}

/// A card that shows a Viewer column with rows also says how the viewer
/// opens: each section places [`Line::ViewerNote`] itself, so this catches
/// one that forgot.
fn viewer_note_placed(c: &Card) -> bool {
    let shows_viewer = c.body.iter().any(|l| match l {
        Line::Table(t) | Line::Folded(_, t) => t.viewer && !t.rows.is_empty(),
        _ => false,
    });
    !shows_viewer || c.body.contains(&Line::ViewerNote)
}

/// One section's frame.
pub(crate) fn frame(sid: &str, part: &HtmlPart, meta: Option<String>, always: bool) -> Vec<String> {
    debug_assert!(
        part.cards.iter().all(viewer_note_placed),
        "{sid}: a card shows a Viewer column without the line saying how it opens"
    );
    let wide = part.cards.len() > 1 || part.cards.iter().any(|c| c.full);
    let mut out = vec![format!(
        "  <section class=\"sec{}{}\" id=\"{sid}\">",
        if wide { " full-width" } else { "" },
        if always { " always" } else { "" }
    )];
    let mut meta = meta;
    for c in &part.cards {
        out.push(format!(
            "  <div class=\"card{}\">",
            if c.full { " full-width" } else { "" }
        ));
        out.push(format!("    <h2>{}</h2>", e(&c.title)));
        out.extend(meta.take());
        out.extend(
            c.caveats
                .iter()
                .map(|n| format!("    <p class=\"caveat\">{}</p>", e(n))),
        );
        out.extend(
            c.stale
                .iter()
                .map(|n| format!("    <div class=\"stale\">{}</div>", e(n))),
        );
        for l in &c.body {
            line(&mut out, l);
        }
        out.push("  </div>".to_owned());
    }
    if !part.files.is_empty() {
        let links: Vec<String> = part
            .files
            .iter()
            .map(|(h, t)| format!("<a href=\"{}\">{}</a>", e(h), e(t)))
            .collect();
        out.push(format!(
            "    <p class=\"files\">Files of this step, beside this report: {}</p>",
            links.join(" &middot; ")
        ));
    }
    out.push("  </section>".to_owned());
    out
}

/// The part of a section whose summary data the report could not read: one
/// card that says why, so the rest of the report is still written.
fn unreadable(sec: &Sec, why: &Malformed) -> HtmlPart {
    HtmlPart {
        cards: vec![Card {
            title: sec.title.to_owned(),
            stale: vec![format!("could not be read: {why}")],
            body: vec![Line::Stat(
                "Status".to_owned(),
                Inline::Badge("Unreadable".to_owned(), Colour::Gray),
            )],
            ..Card::default()
        }],
        files: Vec::new(),
    }
}

/// The two cards of the limits chapter: the steps that did not run, and
/// what no step assesses.
fn limits(s: &Summary) -> Result<(HtmlPart, HtmlPart), Malformed> {
    let not_run = s.texts("not_run")?;
    let body = if not_run.is_empty() {
        Line::Para("All major steps completed.".to_owned())
    } else {
        Line::List(not_run)
    };
    let one = |title: &str, body: Line| HtmlPart {
        cards: vec![Card {
            title: title.to_owned(),
            body: vec![body],
            ..Card::default()
        }],
        files: Vec::new(),
    };
    Ok((
        one("Steps Not Run", body),
        one(
            "Not Assessed by This Pipeline",
            Line::List(s.texts("not_assessed")?),
        ),
    ))
}

fn chapter_open(body: &mut Vec<String>, cid: &str, title: &str, question: &str) {
    body.push(format!("<section class=\"chapter\" id=\"ch-{cid}\">"));
    body.push(format!("  <h2 class=\"chapter-title\">{}</h2>", e(title)));
    body.push(format!("  <p class=\"question\">{}</p>", e(question)));
    body.push("  <div class=\"grid\">".to_owned());
}

/// The caveats a note line is made of: one, or several joined by spaces
/// (step 11's F_ROH and band notes); empty when the line holds anything
/// else.
fn notes_in(line: &str) -> Vec<&'static str> {
    let Some(mut rest) = line
        .strip_prefix("<p class=\"note\">")
        .and_then(|r| r.strip_suffix("</p>"))
    else {
        return Vec::new();
    };
    let mut found = Vec::new();
    while !rest.is_empty() {
        let next = NOTES.iter().map(|(_, n)| *n).find_map(|n| {
            let en = e(n);
            if rest == en {
                Some((n, ""))
            } else {
                rest.strip_prefix(&*en)
                    .and_then(|r| r.strip_prefix(' '))
                    .map(|r| (n, r))
            }
        });
        let Some((n, r)) = next else {
            return Vec::new();
        };
        found.push(n);
        rest = r;
    }
    found
}

/// The sections whose card prints its caveat as its last paragraph, before
/// the caveat moves up under the title: the chart goes where that paragraph
/// was, at the end of the card (render_report.py's ancestry card ends with
/// `<p>{ANCESTRY_FINE_NOTE}</p>`).
const FIGURE_LAST: &[&str] = &["ancestry-fine"];

/// Puts a chart's figure line into its section's lines, under its caveat:
/// right after the line whose sentence the caveat is part of (the PRS
/// caption), at the end of the first card for the sections of
/// [`FIGURE_LAST`], else under the first card's title, after that card's
/// step line and caveats.
fn place_figure(lines: &mut Vec<String>, sid: &str, caveat: &str, figure: String) {
    if FIGURE_LAST.contains(&sid) {
        let h2 = lines.iter().position(|x| x.trim().starts_with("<h2>"));
        let end = h2.and_then(|h| {
            lines[h..]
                .iter()
                .position(|x| x == "  </div>")
                .map(|i| h + i)
        });
        if let Some(i) = end {
            lines.insert(i, figure);
            return;
        }
    }
    let cav = e(caveat);
    let inside = lines.iter().position(|x| {
        let t = x.trim();
        t.contains(&*cav) && !t.starts_with("<p class=\"caveat\">") && notes_in(t).is_empty()
    });
    let at = match inside {
        Some(i) => i + 1,
        None => {
            let Some(h2) = lines.iter().position(|x| x.trim().starts_with("<h2>")) else {
                return;
            };
            let mut i = h2 + 1;
            if lines
                .get(i)
                .is_some_and(|x| x.starts_with("    <p class=\"sec-meta\">"))
            {
                i += 1;
            }
            while lines
                .get(i)
                .is_some_and(|x| x.starts_with("    <p class=\"caveat\">"))
            {
                i += 1;
            }
            i
        }
    };
    lines.insert(at, figure);
}

/// The closing appendix: one line per gene symbol or rsID of the search
/// index, with a link to each section that names it when the report shows
/// that section. No script: the browser's find searches it; the live site's
/// search box filters the same index.
fn terms_appendix(index: &Value, shown: &[&str]) -> Result<Vec<String>, Malformed> {
    let mut body = vec![
        "  <div class=\"card\">".to_owned(),
        "    <p class=\"note\">Every gene symbol and rsID a section of this report names, with the sections \
         that name it and how many times. A gene listed here is a gene the results mention, not a finding: open \
         the section to see what it says. The genes of the clinical filter (step 23) and variant prioritization \
         (step 31) come from their tables, which the report counts but does not list.</p>"
            .to_owned(),
    ];
    let path = "index";
    let terms = crate::model::get(index, "terms", path)?.filter(|v| crate::model::truthy(v));
    let terms = match terms {
        None => None,
        Some(Value::Object(m)) => Some(m),
        Some(_) => return Err(Malformed("index.terms is not an object".to_owned())),
    };
    match terms {
        None => {
            body.push("    <p>No section of this report names a gene or an rsID.</p>".to_owned())
        }
        Some(terms) => {
            body.push("    <ul class=\"plain terms\">".to_owned());
            for (term, items) in terms {
                let Value::Array(items) = items else {
                    return Err(Malformed(format!("index.terms.{term} is not a list")));
                };
                let mut links = Vec::new();
                for it in items {
                    let ipath = format!("index.terms.{term}[]");
                    let need = |k: &str| crate::model::need(it, k, &ipath);
                    let mut label = e(&py(need("title")?)?).into_owned();
                    if crate::model::is_str(Some(need("state")?), "stale") {
                        label.push_str(" (stale)");
                    }
                    let anchor = crate::model::get(it, "anchor", &ipath)?
                        .filter(|a| crate::model::truthy(a))
                        .map(py)
                        .transpose()?;
                    let mut link = match anchor {
                        Some(a) if shown.contains(&a.as_str()) => {
                            format!("<a href=\"#{}\">{label}</a>", e(&a))
                        }
                        _ => label,
                    };
                    let count = need("count")?;
                    if crate::model::num(Some(count)).is_some_and(|n| n > 1.0) {
                        link.push_str(&format!(" ({})", py(count)?));
                    }
                    links.push(link);
                }
                body.push(format!(
                    "      <li><strong>{}</strong>: {}</li>",
                    e(term),
                    links.join(", ")
                ));
            }
            body.push("    </ul>".to_owned());
        }
    }
    body.push("  </div>".to_owned());
    Ok(body)
}

/// The search index of a summary read alone: the step 23 and 31 tables,
/// which need the sample folder, are left out.
pub fn summary_index(s: &Summary) -> Result<Value, Malformed> {
    crate::index::build(s.root(), None, &[], &mut std::io::sink()).map_err(Malformed)
}

/// The HTML report of one summary, with its closing appendix from the
/// summary's own search index.
pub fn report(s: &Summary) -> Result<String, Malformed> {
    let idx = guarded(|| summary_index(s)).ok();
    report_using(s, idx.as_ref(), writer).map(|(page, _)| page)
}

/// The HTML report of one summary with the section writers of `writers`,
/// the closing appendix listing `index` (none without one), and the parts
/// that could not be written. A section whose data cannot be read,
/// or whose writer stops on an internal error, shows as failed with why;
/// charts or an appendix that cannot be made are left out; the rest of the
/// report is still written.
pub fn report_using(
    s: &Summary,
    index: Option<&Value>,
    writers: WriterFn,
) -> Result<(String, Vec<Failed>), Malformed> {
    let mut failed = Vec::new();
    let mut fail = |part: &str, why: &Malformed| {
        failed.push(Failed {
            part: part.to_owned(),
            why: why.to_string(),
        });
    };
    let specs = guarded(|| charts::build(s)).unwrap_or_else(|why| {
        fail("charts", &why);
        Vec::new()
    });
    let mut shown_ids: Vec<&str> = Vec::new();
    let ns = s.new_since();
    // The contents sit in a details element: closed under 900 px, so a phone
    // opens on the results; above that the CSS shows them without the
    // summary line, in a browser that can style ::details-content.
    let mut toc = vec![
        "<details class=\"toc-fold\">".to_owned(),
        "<summary>Contents</summary>".to_owned(),
        "<nav class=\"toc\" aria-label=\"Contents\">".to_owned(),
        "  <h2 class=\"toc-title\">Contents</h2>".to_owned(),
        "  <ol>".to_owned(),
    ];
    let mut body = Vec::new();
    for chapter in &CHAPTERS {
        let mut shown = Vec::new();
        for sec in chapter.sections {
            if let Some((_, html)) = writers(sec.id) {
                // A section whose data could not be read shows as failed in
                // the contents and its badge, whatever its keys' states say.
                let (part, failed) = match guarded(|| html(s)) {
                    Ok(part) => (part, false),
                    Err(why) => {
                        fail(sec.id, &why);
                        (unreadable(sec, &why), true)
                    }
                };
                if !part.cards.is_empty() {
                    shown.push((sec, part, failed));
                }
            }
        }
        if shown.is_empty() {
            continue;
        }
        toc.push(format!(
            "    <li><a href=\"#ch-{}\">{}</a>",
            chapter.id,
            e(chapter.title)
        ));
        toc.push("      <ol>".to_owned());
        chapter_open(&mut body, chapter.id, chapter.title, chapter.question);
        for (sec, part, failed) in shown {
            let state = if failed {
                FrameState::Failed
            } else {
                let states: Vec<State> = sec.keys.iter().map(|k| s.section(k).state).collect();
                FrameState::of(&states)
            };
            toc.push(format!(
                "        <li><a href=\"#{}\"><span class=\"dot dot-{}\" title=\"{}\"></span>{}</a></li>",
                sec.id,
                state.class(),
                state.badge().0,
                e(sec.title)
            ));
            let since = if ns.any_of(sec.keys) {
                Some(ns.previous()?)
            } else {
                None
            };
            let m = meta(s, sec.id, sec.keys, state, since.as_deref())?;
            let mut lines = frame(sec.id, &part, Some(m), false);
            for sp in &specs {
                if section_of(sp.section).is_some_and(|x| x.id == sec.id)
                    && let Some(fig) = charts::figure(sp)
                {
                    place_figure(&mut lines, sec.id, sp.caveat, fig);
                }
            }
            body.extend(lines);
            shown_ids.push(sec.id);
        }
        toc.push("      </ol>".to_owned());
        toc.push("    </li>".to_owned());
        body.push("  </div>".to_owned());
        body.push("</section>".to_owned());
    }
    let (cid, title, question) = LIMITS;
    let (not_run, not_assessed) = limits(s)?;
    toc.extend([
        format!("    <li><a href=\"#ch-{cid}\">{}</a>", e(title)),
        "      <ol>".to_owned(),
        "        <li><a href=\"#not-run\">Steps not run</a></li>".to_owned(),
        "        <li><a href=\"#not-assessed\">Not assessed</a></li>".to_owned(),
        "      </ol>".to_owned(),
        "    </li>".to_owned(),
        "  </ol>".to_owned(),
        "</nav>".to_owned(),
        "</details>".to_owned(),
    ]);
    chapter_open(&mut body, cid, title, question);
    body.extend(frame("not-run", &not_run, None, false));
    body.extend(frame("not-assessed", &not_assessed, None, true));
    body.push("  </div>".to_owned());
    body.push("</section>".to_owned());
    if let Some(index) = index {
        match guarded(|| terms_appendix(index, &shown_ids)) {
            Ok(cards) => {
                let (cid, title, question) = TERMS;
                let at = toc.len() - 3;
                toc.insert(
                    at,
                    format!("    <li><a href=\"#ch-{cid}\">{}</a></li>", e(title)),
                );
                chapter_open(&mut body, cid, title, question);
                body.push("  <section class=\"sec\" id=\"terms-list\">".to_owned());
                body.extend(cards);
                body.push("  </section>".to_owned());
                body.push("  </div>".to_owned());
                body.push("</section>".to_owned());
            }
            Err(why) => fail("terms", &why),
        }
    }

    let mut out: Vec<String> = [
        "<!DOCTYPE html>",
        "<html lang=\"en\">",
        "<head>",
        "<meta charset=\"UTF-8\">",
    ]
    .map(str::to_owned)
    .to_vec();
    out.push(format!(
        "<meta http-equiv=\"Content-Security-Policy\" content=\"{CSP}\">"
    ));
    out.push(
        "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\">".to_owned(),
    );
    out.push("<title>OwnGenome Report</title>".to_owned());
    out.push(format!("<style>{CSS}{READING}</style>"));
    out.extend(["</head>", "<body>", SKIP_LINK, "<div class=\"layout\">"].map(str::to_owned));
    out.extend(toc);
    out.push("<main id=\"main\" tabindex=\"-1\">".to_owned());
    let mut header = run::html_header(s)?;
    let at = header
        .iter()
        .position(|l| l.trim_start().starts_with("<h1>"))
        .map_or(1, |i| i + 1)
        .min(header.len());
    header.insert(
        at,
        format!("  <p class=\"research-only\">{}</p>", e(RESEARCH_ONLY)),
    );
    out.extend(header);
    out.push(PRINT_HINT.to_owned());
    scroll_tables(&mut body, "Table");
    out.extend(body);
    out.extend(
        [
            "<div class=\"disclaimer\">",
            "  <strong>Disclaimer:</strong> This report is for educational and research purposes only.",
            "  It is not a clinical diagnosis. Always discuss genomic findings with a qualified healthcare",
            "  professional before making any medical decisions. Variants of Uncertain Significance (VUS)",
            "  are not clinically actionable.",
            "</div>",
            "",
            "<div class=\"footer\">",
        ]
        .map(str::to_owned),
    );
    out.push(format!(
        "  Generated by {}",
        offsite_link("https://github.com/GeiserX/OwnGenome", "OwnGenome", OFFSITE)
    ));
    out.push("  &mdash; 100% local analysis, no data uploaded".to_owned());
    let mut manifest = run::html_manifest(s)?;
    scroll_tables(&mut manifest, "Run manifest");
    out.extend(manifest);
    out.extend(["</div>", "</main>", "</div>", "</body>", "</html>"].map(str::to_owned));
    let mut page = out.join("\n");
    page.push('\n');
    Ok((page, failed))
}
