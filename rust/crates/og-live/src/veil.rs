//! Match result references by exact locus, then by gene or overlapping locus.
//! The browser's report gate uses findings. These tags also keep related chart
//! marks and the home page's closed lead consistent with that gate.

use crate::headlines::LeadPart;
use crate::results::{Finding, Level, data, list, rows, str_val, string, text, val};
use crate::state::{Obj, obj, s};
use og_render::{
    charts,
    index::{Closable, veil_of},
    model::Summary,
};
use pycompat::json::Value;

struct Locus {
    chrom: String,
    start: i128,
    end: i128,
}

impl Locus {
    fn new(chrom: &str, a: i128, b: i128) -> Self {
        Self {
            chrom: if chrom.starts_with("chr") {
                chrom.to_owned()
            } else {
                format!("chr{chrom}")
            },
            start: a.min(b),
            end: a.max(b),
        }
    }

    fn parts(self) -> (String, i128, i128) {
        (self.chrom, self.start, self.end)
    }
}

fn loci(text: &str) -> Vec<Locus> {
    let bytes = text.as_bytes();
    let position = |start: usize| -> Option<(usize, i128)> {
        if !bytes.get(start).is_some_and(u8::is_ascii_digit) {
            return None;
        }
        let len = bytes[start..]
            .iter()
            .take_while(|b| b.is_ascii_digit() || **b == b',')
            .count();
        let end = start + len;
        Some((end, text[start..end].replace(',', "").parse().ok()?))
    };
    let mut out = Vec::new();
    let mut from = 0;
    for (colon, _) in text.match_indices(':') {
        let mut start = colon;
        while start > from && bytes[start - 1].is_ascii_alphanumeric() {
            start -= 1;
        }
        let chrom = &text[start..colon];
        let name = chrom.strip_prefix("chr").unwrap_or(chrom);
        if !(matches!(name, "X" | "Y" | "M" | "MT")
            || ((1..=2).contains(&name.len()) && name.bytes().all(|b| b.is_ascii_digit())))
        {
            continue;
        }
        let Some((end, a)) = position(colon + 1) else {
            continue;
        };
        let range = (bytes.get(end) == Some(&b'-'))
            .then(|| position(end + 1))
            .flatten();
        let (next, b) = range.unwrap_or((end, a));
        out.push(Locus::new(chrom, a, b));
        from = next;
    }
    out
}

pub(crate) fn spans(found: &[Finding]) -> Vec<Closable> {
    found
        .iter()
        .filter(|f| f.level == Level::Act)
        .map(|f| Closable {
            id: f.id.clone(),
            gene: match &f.gene {
                Value::Null => None,
                v => {
                    let gene = text(v).to_uppercase();
                    (!gene.is_empty()).then_some(gene)
                }
            },
            locus: f
                .viewer
                .iter()
                .chain(&f.detail)
                .find_map(|t| loci(t).into_iter().next())
                .map(Locus::parts),
        })
        .collect()
}

fn span_values(spans: &[Closable]) -> Value {
    Value::Array(
        spans
            .iter()
            .map(|span| {
                let (chrom, start, end) = span.locus.as_ref().map_or(
                    (Value::Null, Value::Null, Value::Null),
                    |(chrom, start, end)| (s(chrom), Value::Int(*start), Value::Int(*end)),
                );
                Value::Array(vec![
                    s(&span.id),
                    span.gene.as_deref().map_or(Value::Null, s),
                    chrom,
                    start,
                    end,
                ])
            })
            .collect(),
    )
}

pub(crate) fn spans_json(found: &[Finding]) -> Value {
    span_values(&spans(found))
}

fn id(spans: &[Closable], text: &str, own: Option<Locus>) -> Option<String> {
    if let Some(locus) = own {
        let text = format!("{text} {}:{}-{}", locus.chrom, locus.start, locus.end);
        return veil_of(spans, &text);
    }
    veil_of(spans, text)
}

/// Tag the legacy search index without rebuilding its terms in Python.
pub(crate) fn tag_search(terms: &mut Value, spans: &[Closable]) -> Result<(), String> {
    let Value::Object(terms) = terms else {
        return Err("search terms must be an object".into());
    };
    for (term, items) in terms {
        let Value::Array(items) = items else {
            return Err("search items must be an array".into());
        };
        for item in items {
            let mut shown = term.clone();
            for mention in list(val(item, "mentions")) {
                shown.push(' ');
                shown.push_str(string(mention, "text"));
                shown.push(' ');
                shown.push_str(string(mention, "locus"));
            }
            let Value::Object(item) = item else {
                return Err("search item must be an object".into());
            };
            if let Some(id) = veil_of(spans, &shown) {
                item.insert("veil".into(), Value::Str(id));
            } else {
                item.shift_remove("veil");
            }
        }
    }
    Ok(())
}

fn joined(v: &Value) -> String {
    list(v).iter().map(text).collect::<Vec<_>>().join(" ")
}

fn explicit(row: &Value) -> Option<Locus> {
    let chrom = string(row, "chrom");
    if chrom.is_empty() {
        return None;
    }
    let pos = match val(row, "pos") {
        Value::Int(n) => *n,
        Value::Float(n) if n.is_finite() => *n as i128,
        v => str_val(v).trim().parse().ok()?,
    };
    Some(Locus::new(chrom, pos, pos))
}

fn tags(name: &str, map: &mut Obj, ids: impl IntoIterator<Item = Option<String>>) {
    let ids: Vec<_> = ids.into_iter().collect();
    if ids.iter().any(Option::is_some) {
        map.insert(name.into(), tags_value(ids));
    }
}

fn tags_value(ids: impl IntoIterator<Item = Option<String>>) -> Value {
    Value::Array(
        ids.into_iter()
            .map(|id| id.map_or(Value::Null, Value::Str))
            .collect(),
    )
}

pub(crate) fn build(
    summary: &Summary,
    found: &[Finding],
    cards: &[Obj],
    lead: &[LeadPart],
) -> Value {
    let spans = spans(found);
    if spans.is_empty() {
        return obj([]);
    }
    let mut row_tags = Obj::new();
    let mut viewer = Obj::new();
    let mut lines = Obj::new();
    let mut facts = Obj::new();
    let mut widgets = Obj::new();
    let mut marks = Obj::new();
    let mut findings = Obj::new();
    let hit_tags: Vec<_> = rows(val(data(summary, "clinvar"), "hits"))
        .iter()
        .take(og_render::sections::clinvar::HITS_SHOWN)
        .map(|h| id(&spans, string(h, "gene"), explicit(h)))
        .collect();
    tags("clinvar", &mut row_tags, hit_tags.iter().cloned());
    for card in cards {
        let name = str_val(&card["id"]);
        if name != "clinvar" {
            tags(
                name,
                &mut row_tags,
                list(val(&card["table"], "rows"))
                    .iter()
                    .map(|r| id(&spans, &joined(r), None)),
            );
        }
        if let Some(items) = card.get("viewer") {
            tags(
                name,
                &mut viewer,
                list(items).iter().enumerate().map(|(i, locus)| {
                    (name == "clinvar")
                        .then(|| hit_tags.get(i).cloned().flatten())
                        .flatten()
                        .or_else(|| id(&spans, str_val(locus), None))
                }),
            );
        }
        tags(
            name,
            &mut lines,
            list(&card["lines"])
                .iter()
                .map(|line| id(&spans, str_val(line), None)),
        );
        tags(
            name,
            &mut facts,
            list(&card["facts"])
                .iter()
                .map(|fact| id(&spans, &joined(fact), None)),
        );
        if let Some(widget) = card.get("widget") {
            let items = match string(widget, "kind") {
                "vaf" => list(val(widget, "calls")),
                "copies" => list(val(widget, "genes")),
                _ => &[],
            };
            tags(
                name,
                &mut widgets,
                items.iter().map(|item| {
                    id(
                        &spans,
                        &format!("{} {}", string(item, "gene"), string(item, "change")),
                        None,
                    )
                }),
            );
        }
    }
    if let Ok(specs) = og_render::chapters::guarded(|| charts::build(summary)) {
        for spec in specs {
            let tagged: Obj = spec
                .marks
                .iter()
                .filter(|m| m.kind != Some(charts::MarkKind::Roh))
                .filter_map(|m| {
                    let own = m
                        .chrom
                        .as_deref()
                        .zip(m.pos)
                        .map(|(c, p)| Locus::new(c, p, m.end.unwrap_or(p)));
                    id(&spans, &format!("{} {}", m.label, m.detail), own)
                        .map(|id| (m.id.clone(), Value::Str(id)))
                })
                .collect();
            if !tagged.is_empty() {
                marks.insert(spec.id.to_string(), Value::Object(tagged));
            }
        }
    }
    for f in found.iter().filter(|f| f.level != Level::Act) {
        let t = format!(
            "{} {} {} {}",
            str_val(&f.gene),
            f.detail.as_deref().unwrap_or(""),
            f.title,
            f.viewer.as_deref().unwrap_or("")
        );
        if let Some(link) = id(&spans, &t, None) {
            findings.insert(f.id.clone(), Value::Str(link));
        }
    }
    let lead: Vec<_> = lead
        .iter()
        .map(|part| {
            if part.names_act {
                None
            } else {
                id(&spans, &part.text, None)
            }
        })
        .collect();
    let mut value: Obj = [
        ("rows", row_tags),
        ("viewer", viewer),
        ("lines", lines),
        ("facts", facts),
        ("widget", widgets),
        ("marks", marks),
        ("findings", findings),
    ]
    .into_iter()
    .filter(|(_, v)| !v.is_empty())
    .map(|(key, v)| (key.into(), Value::Object(v)))
    .collect();
    if lead.iter().any(Option::is_some) {
        value.insert("lead".into(), tags_value(lead));
    }
    value.insert("spans".into(), span_values(&spans));
    Value::Object(value)
}
