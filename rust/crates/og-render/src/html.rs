//! The HTML report: one page with no script and no fetch, a table of
//! contents with a state dot per section, the chapters, what the report
//! leaves out, the disclaimer and the run manifest.
//!
//! Each section is a frame: an anchor, a line with the step, its state, the
//! newest file date and "New since" under the first card title, each card's
//! caveats under its title, the cards, then links to the step's own files.

use crate::chapters::{CHAPTERS, FrameState, LIMITS, Sec, step_of, writer};
use crate::css::CSS;
use crate::escape::e;
use crate::model::{Card, Colour, HtmlPart, Inline, Line, Malformed, Section, State, Summary, py};
use crate::notes::stale_note;
use crate::sections::run;

/// The Content-Security-Policy: nothing is fetched, images only inline,
/// styles only from the page.
pub const CSP: &str = "default-src 'none'; img-src data:; style-src 'unsafe-inline'";

/// A link that leaves this computer, marked so a reader knows before
/// opening it.
pub const OFFSITE: &str = "opens a website";

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
        Inline::Badge(t, c) => badge(t, *c),
        Inline::Quiet(t) => format!(
            "<span style=\"font-weight:normal;font-size:13px\">{}</span>",
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
        Line::Note(t) => out.push(format!("    <p class=\"note\">{}</p>", e(t))),
        Line::Para(t) => out.push(format!("    <p>{}</p>", e(t))),
        Line::List(items) => {
            out.push("    <ul class=\"plain\">".to_owned());
            out.extend(items.iter().map(|t| format!("      <li>{}</li>", e(t))));
            out.push("    </ul>".to_owned());
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

/// The step line under a section's first card title.
fn meta(s: &Summary, sid: &str, keys: &[&str], since: Option<&str>) -> Result<String, Malformed> {
    let secs: Vec<Section<'_>> = keys.iter().map(|k| s.section(k)).collect();
    let mut steps: Vec<&str> = Vec::new();
    for k in keys {
        let step = step_of(k);
        if !steps.contains(&step) {
            steps.push(step);
        }
    }
    let (label, colour) = FrameState::of(&secs.iter().map(|x| x.state).collect::<Vec<_>>()).badge();
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

/// One section's frame.
fn frame(sid: &str, part: &HtmlPart, meta: Option<String>, always: bool) -> Vec<String> {
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

/// The HTML report of one summary.
pub fn report(s: &Summary) -> Result<String, Malformed> {
    let ns = s.new_since();
    let mut toc = vec![
        "<nav class=\"toc\" aria-label=\"Contents\">".to_owned(),
        "  <h2 class=\"toc-title\">Contents</h2>".to_owned(),
        "  <ol>".to_owned(),
    ];
    let mut body = Vec::new();
    for chapter in &CHAPTERS {
        let mut shown = Vec::new();
        for sec in chapter.sections {
            if let Some((_, html)) = writer(sec.id) {
                let part = html(s).unwrap_or_else(|why| unreadable(sec, &why));
                if !part.cards.is_empty() {
                    shown.push((sec, part));
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
        for (sec, part) in shown {
            let states: Vec<State> = sec.keys.iter().map(|k| s.section(k).state).collect();
            let state = FrameState::of(&states);
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
            let m = meta(s, sec.id, sec.keys, since.as_deref())?;
            body.extend(frame(sec.id, &part, Some(m), false));
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
    ]);
    chapter_open(&mut body, cid, title, question);
    body.extend(frame("not-run", &not_run, None, false));
    body.extend(frame("not-assessed", &not_assessed, None, true));
    body.push("  </div>".to_owned());
    body.push("</section>".to_owned());

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
    out.push(format!("<style>{CSS}</style>"));
    out.extend(["</head>", "<body>", "<div class=\"layout\">"].map(str::to_owned));
    out.extend(toc);
    out.push("<main>".to_owned());
    out.extend(run::html_header(s)?);
    out.push(PRINT_HINT.to_owned());
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
    out.extend(run::html_manifest(s)?);
    out.extend(["</div>", "</main>", "</div>", "</body>", "</html>"].map(str::to_owned));
    let mut page = out.join("\n");
    page.push('\n');
    Ok(page)
}
