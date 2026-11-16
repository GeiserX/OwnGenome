//! The text report: the run facts, then each chapter under a banner with
//! its question, then the steps that did not run, what no step assesses,
//! the manifest line and the disclaimer.

use crate::chapters::{CHAPTERS, Sec, writer};
use crate::model::{Malformed, State, Summary, py};
use crate::notes::stale_note;
use crate::sections::run;

/// A section block's opening lines, and whether its step ran (its result
/// was read, from this run or an earlier one).
#[derive(Debug)]
pub struct Head {
    pub lines: Vec<String>,
    pub ran: bool,
}

/// `  [New since <date>]` when one of `keys` changed since the previous run.
pub fn new_since(s: &Summary, keys: &[&str]) -> Result<Option<String>, Malformed> {
    let ns = s.new_since();
    if keys.iter().any(|k| ns.lists(k)) {
        Ok(Some(format!("  [New since {}]", ns.previous()?)))
    } else {
        Ok(None)
    }
}

/// The opening of summary key `key`'s block: its title (or `title`), a rule,
/// the new-since mark and the state note. `None` when the step did not run,
/// so the block is left out.
pub fn head(s: &Summary, key: &str, title: Option<&str>) -> Result<Option<Head>, Malformed> {
    let sec = s.section(key);
    if sec.state == State::Missing {
        return Ok(None);
    }
    let title = match title {
        Some(t) => t.to_owned(),
        None => py(sec.title)?,
    };
    let mut lines = vec![format!("## {title}"), "---".to_owned()];
    if let Some(line) = new_since(s, &[key])? {
        lines.push(line);
    }
    let note = stale_note(&sec)?;
    if !note.is_empty() {
        lines.push(format!("  [{note}]"));
    }
    if sec.state == State::Failed {
        lines.push(String::new());
    }
    Ok(Some(Head {
        lines,
        ran: sec.state.ran(),
    }))
}

/// The block of a section whose summary data the report could not read: its
/// title and why, so the rest of the report is still written.
fn unreadable(sec: &Sec, why: &Malformed) -> Vec<String> {
    vec![
        format!("## {}", sec.title),
        "---".to_owned(),
        format!("  [could not be read: {why}]"),
        String::new(),
    ]
}

/// The text report of one summary.
pub fn report(s: &Summary) -> Result<String, Malformed> {
    let mut lines = run::text_head(s)?;
    for (n, chapter) in CHAPTERS.iter().enumerate() {
        let mut body = Vec::new();
        for sec in chapter.sections {
            if let Some((text, _)) = writer(sec.id) {
                body.extend(text(s).unwrap_or_else(|why| unreadable(sec, &why)));
            }
        }
        if !body.is_empty() {
            let banner = "#".repeat(80);
            lines.push(banner.clone());
            lines.push(format!(
                "#  {}. {}: {}",
                n + 1,
                pycompat::text::upper(chapter.title),
                chapter.question
            ));
            lines.push(banner);
            lines.push(String::new());
            lines.extend(body);
        }
    }
    lines.push("## Steps Not Run".to_owned());
    lines.push("---".to_owned());
    let not_run = s.texts("not_run")?;
    if not_run.is_empty() {
        lines.push("  All major steps completed.".to_owned());
    }
    lines.extend(not_run.iter().map(|t| format!("  - {t}")));
    lines.push(String::new());
    lines.push("## Not Assessed by This Pipeline".to_owned());
    lines.push("---".to_owned());
    lines.extend(s.texts("not_assessed")?.iter().map(|t| format!("  - {t}")));
    lines.push(String::new());
    lines.extend(run::text_manifest(s)?);
    let rule = "=".repeat(80);
    lines.push(rule.clone());
    lines.push("  DISCLAIMER: This is NOT a clinical report. Discuss findings with a".to_owned());
    lines.push("  qualified healthcare professional before making medical decisions.".to_owned());
    lines.push(rule);
    let mut out = lines.join("\n");
    out.push('\n');
    Ok(out)
}
