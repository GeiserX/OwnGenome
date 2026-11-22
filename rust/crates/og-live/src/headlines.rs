//! The journal lead and the home page's open and closed versions.

use crate::results::{Finding, Level, data, list, str_val, string, text, val};
use crate::state::{get, obj, s};
use og_render::model::{Summary, truthy};
use pycompat::{fmt::fixed, json::Value};
use std::collections::BTreeSet;

pub(crate) struct LeadPart {
    pub text: String,
    bold: bool,
    pub names_act: bool,
}

impl LeadPart {
    fn new(text: impl Into<String>, bold: bool) -> Self {
        Self {
            text: text.into(),
            bold,
            names_act: false,
        }
    }

    pub(crate) fn json(self) -> Value {
        let mut part = vec![s(&self.text), Value::Bool(self.bold)];
        if self.names_act {
            part.push(s("act"));
        }
        Value::Array(part)
    }
}

fn screened(summary: &Summary) -> bool {
    ["clinvar", "clinical", "cpsr"]
        .iter()
        .any(|key| summary.section(key).state.ran())
}

fn names(found: &[&Finding], limit: usize) -> String {
    let mut titles = found
        .iter()
        .take(limit)
        .map(|f| f.title.as_str())
        .collect::<Vec<_>>();
    let Some(last) = titles.pop() else {
        return String::new();
    };
    let mut names = if titles.is_empty() {
        last.to_owned()
    } else {
        format!("{} and {last}", titles.join(", "))
    };
    if found.len() > limit {
        names.push_str(&format!(", and {} more", found.len() - limit));
    }
    names
}

fn count(n: usize, one: &str, many: &str, capitalize: bool) -> String {
    let word = [
        "zero", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine",
    ]
    .get(n)
    .map_or_else(|| n.to_string(), |word| (*word).to_owned());
    let mut text = format!("{word} {}", if n == 1 { one } else { many });
    if capitalize {
        text[..1].make_ascii_uppercase();
    }
    text
}

pub(crate) fn lead(summary: &Summary, found: &[Finding]) -> Vec<LeadPart> {
    let by = |lv| found.iter().filter(|f| f.level == lv).collect::<Vec<_>>();
    let mut lead = Vec::new();
    let trust = by(Level::Trust);
    if !trust.is_empty() {
        lead.push(LeadPart::new("Read the data-quality section first: ", true));
        lead.push(LeadPart::new(format!("{}. ", names(&trust, 3)), false));
    }
    let act = by(Level::Act);
    if !act.is_empty() {
        lead.push(LeadPart::new(
            "Of everything this genome was screened for, ",
            false,
        ));
        lead.push(LeadPart::new(
            format!(
                "{} go to a doctor",
                count(act.len(), "finding should", "findings should", false)
            ),
            true,
        ));
        lead.push(LeadPart {
            text: format!(": {}. ", names(&act, 3)),
            bold: false,
            names_act: true,
        });
    } else if screened(summary) {
        lead.push(LeadPart::new(
            "No finding in what these steps assessed needs a doctor; ",
            true,
        ));
        lead.push(LeadPart::new(
            "sections 8 and 9 say what they do not see. ",
            false,
        ));
    } else {
        lead.push(LeadPart::new(
            "The clinical screens have no result yet; ",
            true,
        ));
        lead.push(LeadPart::new(
            "section 9 lists the steps that did not run. ",
            false,
        ));
    }
    let medicine = by(Level::Medicine);
    if !medicine.is_empty() {
        let genes = medicine
            .iter()
            .filter(|f| truthy(&f.gene))
            .map(|f| text(&f.gene))
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        let verb = if medicine.len() == 1 { "bears" } else { "bear" };
        lead.push(LeadPart::new(
            count(medicine.len(), "medicine result", "medicine results", true),
            true,
        ));
        lead.push(LeadPart::new(
            format!(
                " {verb} on how drugs work ({}): tell whoever prescribes. ",
                genes.join(", ")
            ),
            false,
        ));
    }
    let carrier = by(Level::Carrier);
    if !carrier.is_empty() {
        lead.push(LeadPart::new(
            format!("A carrier of {}", names(&carrier, 2)),
            true,
        ));
        lead.push(LeadPart::new(
            ": carriers are not affected; it matters when a partner carries the same gene. ",
            false,
        ));
    }
    let review = by(Level::Review);
    if !review.is_empty() {
        lead.push(LeadPart::new(
            format!(
                "{} worth a look",
                count(review.len(), "result is", "results are", true)
            ),
            true,
        ));
        lead.push(LeadPart::new(format!(": {}. ", names(&review, 3)), false));
    }
    let coverage = data(summary, "coverage");
    let qc = data(summary, "sample_qc");
    if trust.is_empty()
        && let Some(depth) = match val(coverage, "mean_depth") {
            Value::Int(n) => Some(*n as f64),
            Value::Float(n) => Some(*n),
            Value::Bool(n) => Some(f64::from(u8::from(*n))),
            _ => None,
        }
    {
        let checked = string(qc, "contamination") == "ok" && string(qc, "sex_check") == "ok";
        let note = if checked {
            "the identity and contamination checks raised no flag."
        } else {
            "the identity and contamination checks have no result here: section 5 says why."
        };
        lead.push(LeadPart::new(
            format!("Mean depth {}×; {note}", fixed(depth, 1)),
            false,
        ));
    }
    if let Some(last) = lead.last_mut() {
        last.text.truncate(last.text.trim_end().len());
    }
    lead
}

fn closed_segment(t: &str) -> String {
    const CLOSED: &str = "closed until opened on purpose";
    if t.trim_start().starts_with(':') {
        return format!(": {CLOSED}. ");
    }
    if let Some((start, end)) = t.find('(').zip(t.rfind(')')).filter(|(a, b)| a < b) {
        return format!("{}({CLOSED}){}", &t[..start], &t[end + 1..]);
    }
    if let Some((start, _)) = t.split_once(" of ") {
        return format!("{start} result, {CLOSED}");
    }
    format!("{CLOSED}. ")
}

pub(crate) fn overview(summary: &Summary, doc: &Value) -> Value {
    let veil = val(doc, "veil");
    let named = val(veil, "findings");
    let mut top: Vec<_> = list(val(doc, "findings"))
        .iter()
        .filter(|f| string(f, "level") != "info" && get(named, string(f, "id")).is_none())
        .collect();
    top.sort_by_key(|f| match string(f, "level") {
        "trust" => 0,
        "act" => 1,
        "medicine" => 2,
        "review" => 3,
        _ => 4,
    });
    let tags = list(val(veil, "lead"));
    let mut lead = String::new();
    let mut closed = String::new();
    for (i, part) in list(val(doc, "lead")).iter().map(list).enumerate() {
        let text = part.first().map_or("", str_val);
        lead.push_str(text);
        if part.get(2).is_some_and(|x| str_val(x) == "act") {
            closed.push_str(": closed until opened on purpose. ");
        } else if tags.get(i).is_some_and(truthy) {
            closed.push_str(&closed_segment(text));
        } else {
            closed.push_str(text);
        }
    }
    obj([
        ("attention", val(doc, "attention").clone()),
        ("levels", val(doc, "levels").clone()),
        ("screened", Value::Bool(screened(summary))),
        (
            "top",
            Value::Array(
                top.into_iter()
                    .take(4)
                    .map(|f| obj(["id", "level", "title"].map(|key| (key, val(f, key).clone()))))
                    .collect(),
            ),
        ),
        ("lead", s(&lead)),
        ("lead_closed", s(&closed)),
    ])
}
