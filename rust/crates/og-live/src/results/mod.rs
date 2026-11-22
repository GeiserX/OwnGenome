//! The live results and home-page overview of a run's summary.

mod clinical;
mod clinvar;
mod lineage;
mod pgx;
mod qc;
mod sv;
mod traits;

use crate::plugins::ResultsFns;
use crate::state::{Obj, get, obj, s};
use crate::{headlines, veil};
use og_render::chapters::{CHAPTERS, FrameState, Sec, step_of};
use og_render::model::{Malformed, State, Summary, num, py, truthy};
use og_render::notes;
use pycompat::json::{Dump, Value, dumps};
use std::collections::{BTreeSet, HashSet};
use std::io::Write;

const EMPTY_CLINVAR: &str = "No P/LP match in this ClinVar file at the call floor. ClinVar lists only part of what can cause disease, so this is information, and the caveats below still apply.";
const TABLE_LIMIT: usize = 200;
const USAGE: &str = "usage: og live results SUMMARY [--sample NAME] [--file PATH]... [--out PATH]\n       og live results SUMMARY --overview | --spans | --tag-search";

/// The result builders used by the live server.
pub const FNS: ResultsFns = ResultsFns { build, overview };

type Read<T> = Result<T, Malformed>;

pub(crate) fn val<'a>(v: &'a Value, key: &str) -> &'a Value {
    get(v, key).unwrap_or(&Value::Null)
}

pub(crate) fn str_val(v: &Value) -> &str {
    match v {
        Value::Str(s) => s,
        _ => "",
    }
}

pub(crate) fn string<'a>(v: &'a Value, key: &str) -> &'a str {
    str_val(val(v, key))
}

pub(crate) fn text(v: &Value) -> String {
    py(v).unwrap_or_default()
}

pub(crate) fn text_or(v: &Value, default: &str) -> String {
    if truthy(v) {
        text(v)
    } else {
        default.to_owned()
    }
}

fn default(v: &Value, key: &str, fallback: &str) -> String {
    get(v, key).map_or_else(|| fallback.to_owned(), text)
}

pub(crate) fn data<'a>(summary: &'a Summary, section: &'a str) -> &'a Value {
    let part = summary.section(section);
    if part.state.ran() && matches!(part.data, Value::Object(_)) {
        part.data
    } else {
        &Value::Null
    }
}

pub(crate) fn list(v: &Value) -> &[Value] {
    match v {
        Value::Array(items) => items,
        _ => &[],
    }
}

pub(crate) fn rows(v: &Value) -> &[Value] {
    if let Value::Array(items) = v {
        items
    } else if truthy(v) {
        std::slice::from_ref(v)
    } else {
        &[]
    }
}

fn array<'a>(v: &'a Value, key: &str) -> Read<&'a [Value]> {
    match val(v, key) {
        Value::Array(items) => Ok(items),
        other if !truthy(other) => Ok(&[]),
        _ => Err(Malformed(format!("{key} is not a list"))),
    }
}

fn fields(v: &Value) -> impl Iterator<Item = (&str, &Value)> {
    let fields = match v {
        Value::Object(fields) => Some(fields),
        _ => None,
    };
    fields
        .into_iter()
        .flatten()
        .map(|(key, value)| (key.as_str(), value))
}

fn row(v: &Value, keys: &[&str]) -> Vec<String> {
    keys.iter().map(|key| cell(val(v, key))).collect()
}

fn plain(text: &str) -> String {
    let folded = text.replace("**", "");
    let mut out = String::new();
    for word in folded.split_whitespace() {
        if !out.is_empty() {
            let passed_check = out.ends_with(':')
                && word.strip_prefix("ok").is_some_and(|tail| {
                    !tail.starts_with(|c: char| c.is_alphanumeric() || c == '_')
                });
            out.push(' ');
            if passed_check {
                out.push_str("passed");
                out.push_str(&word[2..]);
                continue;
            }
        }
        out.push_str(word);
    }
    out
}

fn cell(v: &Value) -> String {
    match v {
        Value::Null => ".".to_owned(),
        Value::Str(t) if t.is_empty() => ".".to_owned(),
        Value::Float(n) => pycompat::fmt::general(*n, 6),
        Value::Array(a) if !a.is_empty() => a.iter().map(cell).collect::<Vec<_>>().join(", "),
        Value::Array(_) => ".".to_owned(),
        _ => text(v),
    }
}

fn int(n: usize) -> Value {
    Value::Int(n as i128)
}

fn strings(items: impl IntoIterator<Item = impl AsRef<str>>) -> Value {
    Value::Array(items.into_iter().map(|item| s(item.as_ref())).collect())
}

fn table(columns: &[&str], rows: Vec<Vec<String>>) -> Value {
    if rows.is_empty() {
        return Value::Null;
    }
    let more = rows.len().saturating_sub(TABLE_LIMIT);
    let mut table = Obj::from([
        ("columns".into(), strings(columns)),
        (
            "rows".into(),
            Value::Array(rows.into_iter().take(TABLE_LIMIT).map(strings).collect()),
        ),
    ]);
    if more > 0 {
        table.insert("more".into(), int(more));
    }
    Value::Object(table)
}

#[derive(Default)]
struct Adapter {
    facts: Vec<Value>,
    table: Option<Value>,
    lines: Vec<String>,
    attention: bool,
    links: Vec<&'static str>,
    extra: Obj,
}

impl Adapter {
    fn fact(&mut self, label: impl AsRef<str>, value: impl AsRef<str>) {
        self.facts
            .push(strings([plain(label.as_ref()), plain(value.as_ref())]));
    }

    fn facts(&mut self, d: &Value, pairs: &[(&str, &str)]) {
        for (label, key) in pairs {
            self.fact(label, cell(val(d, key)));
        }
    }

    fn table(&mut self, columns: &[&str], rows: Vec<Vec<String>>) {
        self.table = Some(table(columns, rows));
    }

    fn fields_table(&mut self, columns: &[&str], items: &[Value], keys: &[&str]) {
        self.table(columns, items.iter().map(|x| row(x, keys)).collect());
    }

    fn extra(&mut self, key: &str, value: Value) {
        if truthy(&value) {
            self.extra.insert(key.into(), value);
        }
    }
}

fn adapter(id: &str, summary: &Summary) -> Read<Adapter> {
    match id {
        "qc" | "callability" | "variants" => qc::card(id, summary),
        "clinvar" => clinvar::card(summary),
        "clinical-filter" | "cpsr" | "slivar" | "exomiser" | "carrier-list" | "couple"
        | "paralogs" => clinical::card(id, summary),
        "structural-variants"
        | "mei"
        | "expansions"
        | "vntyper"
        | "mito"
        | "lpa"
        | "chip"
        | "mosaic"
        | "mosaic-mca" => sv::card(id, summary),
        "pharmcat" | "cpic" | "cyp2d6" | "hla" | "locityper" => pgx::card(id, summary),
        "ancestry-fine" | "mt-haplogroup" | "y-chromosome" | "archaic" | "ancient" | "roh" => {
            lineage::card(id, summary)
        }
        "traits-panel" | "prs" | "bloodtype" | "telomere" | "viral" | "consumer-export" => {
            traits::card(id, summary)
        }
        _ => Ok(Adapter::default()),
    }
}

fn report_links(kinds: &[&str], paths: Option<&BTreeSet<String>>) -> Value {
    let mut links = Vec::new();
    if let Some(paths) = paths {
        for kind in kinds {
            for path in paths {
                let leaf = path.rsplit('/').next().unwrap_or(path);
                let matched = match *kind {
                    "multiqc" | "fastp" => leaf.starts_with(kind) && leaf.ends_with(".html"),
                    "_fastp" => leaf.ends_with("_fastp.html"),
                    "vep" => path.split('/').any(|p| p == "vep") && leaf.ends_with("summary.html"),
                    "cpsr" => path.starts_with("cpsr/"),
                    "pharmcat" => path.rsplit_once('/').is_some_and(|(parent, _)| {
                        parent.rsplit('/').next() == Some("pharmcat") && leaf.ends_with(".html")
                    }),
                    _ => false,
                };
                if matched && !links.contains(&path) {
                    links.push(path);
                    if links.len() == 4 {
                        return strings(links);
                    }
                }
            }
        }
    }
    strings(links)
}

fn card(section: &Sec, chapter: &str, summary: &Summary, files: Option<&BTreeSet<String>>) -> Obj {
    let raw_sections = summary.get("sections").unwrap_or(&Value::Null);
    let sections: Vec<_> = section
        .keys
        .iter()
        .map(|key| summary.section(key))
        .collect();
    let states: Vec<_> = sections.iter().map(|sec| sec.state).collect();
    let state = match FrameState::of(&states) {
        FrameState::Ok => "ran",
        other => other.class(),
    };
    let mut steps = Vec::new();
    let mut state_notes = Vec::new();
    let mut headline = None;
    let mut flagged = false;
    let caveats: Vec<_> = notes::section_notes(section.id)
        .iter()
        .filter_map(|name| notes::text(name))
        .collect();
    for sec in sections {
        let raw = get(raw_sections, sec.key);
        let step = raw
            .and_then(|v| get(v, "step"))
            .filter(|v| truthy(v))
            .cloned()
            .unwrap_or_else(|| s(step_of(sec.key)));
        if step != s("?") && !steps.contains(&step) {
            steps.push(step);
        }
        if matches!(sec.state, State::Failed | State::Unreadable | State::Stale) {
            state_notes.push(format!(
                "{}: {}",
                text_or(sec.title, sec.key),
                plain(&notes::stale_note(&sec).unwrap_or_default())
            ));
        }
        if sec.state.ran()
            && let Some((line, tone)) = crate::state::headline(sec.key, raw, summary)
        {
            flagged |= tone == "attention";
            let line = plain(&line);
            if headline.is_none() && !line.is_empty() && !line.ends_with("finished, see the report")
            {
                headline = Some(line);
            }
        }
    }
    if let Some(line) = &mut headline
        && let Some(stripped) = caveats
            .iter()
            .map(|c| plain(c))
            .filter(|c| {
                line.strip_suffix(c)
                    .is_some_and(|start| start.ends_with(' '))
            })
            .max_by_key(String::len)
    {
        line.truncate(line.len() - stripped.len());
        line.truncate(line.trim_end().len());
    }
    let adapter = if matches!(state, "missing" | "failed") {
        Adapter::default()
    } else {
        og_render::chapters::guarded(|| adapter(section.id, summary)).unwrap_or_else(|err| {
            Adapter {
                lines: vec![format!(
                    "This card could not be drawn ({}); the report shows the section.",
                    error_kind(&err)
                )],
                ..Adapter::default()
            }
        })
    };
    flagged |= adapter.attention;
    let attention = flagged || !state_notes.is_empty();
    let lines = adapter.lines.iter().map(|line| plain(line)).filter(|line| {
        !line.is_empty()
            && headline.as_ref() != Some(line)
            && !caveats.iter().any(|c| plain(c) == *line)
            && !notes::NOTES.iter().any(|(_, c)| plain(c) == *line)
    });
    let mut card = Obj::from([
        ("id".into(), s(section.id)),
        ("title".into(), s(section.title)),
        ("keys".into(), strings(section.keys)),
        ("steps".into(), Value::Array(steps)),
        ("state".into(), s(state)),
        ("anchor".into(), s(section.id)),
        ("attention".into(), Value::Bool(attention)),
        (
            "headline".into(),
            headline.as_deref().map_or(Value::Null, s),
        ),
        ("facts".into(), Value::Array(adapter.facts)),
        ("table".into(), adapter.table.unwrap_or(Value::Null)),
        ("lines".into(), strings(lines)),
        ("caveat".into(), strings(caveats)),
        (
            "empty_text".into(),
            if section.id == "clinvar" {
                s(EMPTY_CLINVAR)
            } else {
                Value::Null
            },
        ),
        ("links".into(), report_links(&adapter.links, files)),
        ("notes".into(), strings(state_notes)),
        ("flagged".into(), Value::Bool(flagged)),
        ("chapter".into(), s(chapter)),
    ]);
    card.extend(adapter.extra);
    card
}

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum Level {
    Act,
    Medicine,
    Review,
    Carrier,
    Trust,
    Info,
}

impl Level {
    const ALL: [Self; 6] = [
        Self::Act,
        Self::Medicine,
        Self::Review,
        Self::Carrier,
        Self::Trust,
        Self::Info,
    ];

    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Act => "act",
            Self::Medicine => "medicine",
            Self::Review => "review",
            Self::Carrier => "carrier",
            Self::Trust => "trust",
            Self::Info => "info",
        }
    }

    fn words(self) -> &'static str {
        match self {
            Self::Act => "Take to a doctor",
            Self::Medicine => "Tell your prescriber",
            Self::Review => "Worth a look",
            Self::Carrier => "For family planning",
            Self::Trust => "Data quality",
            Self::Info => "For information",
        }
    }

    fn section(self) -> usize {
        match self {
            Self::Act => 1,
            Self::Medicine => 2,
            Self::Carrier => 3,
            Self::Review => 4,
            Self::Trust => 5,
            Self::Info => 8,
        }
    }

    fn plan(self) -> &'static [&'static str] {
        match self {
            Self::Act => &[
                "Bring this page to a doctor or a genetic counsellor.",
                "Ask for a clinical test in a certified laboratory to confirm it: this pipeline is research software.",
                "Ask what it means for parents, brothers, sisters and children: close relatives may carry it too.",
            ],
            Self::Medicine => &[
                "Show this to whoever prescribes for you, and to your pharmacist.",
                "Do not change a medicine on your own: the guidance is written for the prescriber.",
                "Ask for a clinical pharmacogenetic test before a prescription changes.",
            ],
            Self::Review => &[
                "Read the evidence and its caveat below.",
                "Talk it over with a doctor or a genetic counsellor if it bears on your health or a pregnancy.",
                "A measurement or a clinical test settles it: this result alone is not a diagnosis.",
            ],
            Self::Carrier => &[
                "A carrier is not affected by the condition.",
                "It matters when a partner carries a change in the same gene: compare with the couple screen.",
                "Relatives planning a pregnancy may carry it too.",
            ],
            Self::Trust => &[
                "Read nothing else in this report until this is understood.",
                "Check the sample's label and the run's inputs, then run it again.",
            ],
            Self::Info => &[],
        }
    }
}

pub(crate) struct Finding {
    pub(crate) id: String,
    pub(crate) level: Level,
    card: &'static str,
    pub(crate) title: String,
    text: String,
    pub(crate) gene: Value,
    pub(crate) detail: Option<String>,
    evidence: Vec<Value>,
    pub(crate) viewer: Option<String>,
    af: Value,
}

impl Finding {
    fn new(
        level: Level,
        card: &'static str,
        key: &[&str],
        title: impl AsRef<str>,
        text: impl AsRef<str>,
    ) -> Self {
        Self {
            id: slug(key),
            level,
            card,
            title: plain(title.as_ref()),
            text: plain(text.as_ref()),
            gene: Value::Null,
            detail: None,
            evidence: Vec::new(),
            viewer: None,
            af: Value::Null,
        }
    }

    fn gene(mut self, gene: &Value) -> Self {
        self.gene = gene.clone();
        self
    }

    fn detail(mut self, detail: impl AsRef<str>) -> Self {
        let detail = plain(detail.as_ref());
        self.detail = (!detail.is_empty()).then_some(detail);
        self
    }

    fn viewer(mut self, viewer: String) -> Self {
        self.viewer = (!viewer.is_empty()).then_some(viewer);
        self
    }

    fn evidence(mut self, label: &str, value: &Value) -> Self {
        if !matches!(value, Value::Null) && !matches!(str_val(value), ".") && value != &s("") {
            self.evidence
                .push(strings([plain(label), plain(&cell(value))]));
        }
        self
    }

    fn evidence_fields(mut self, value: &Value, pairs: &[(&str, &str)]) -> Self {
        for (label, key) in pairs {
            self = self.evidence(label, val(value, key));
        }
        self
    }

    fn json(&self, cards: &[Obj], found: &[Finding]) -> Value {
        let card = cards
            .iter()
            .find(|c| c.get("id").is_some_and(|v| str_val(v) == self.card));
        let card_list = |key| {
            card.and_then(|c| c.get(key))
                .cloned()
                .unwrap_or_else(|| Value::Array(Vec::new()))
        };
        obj([
            ("id", s(&self.id)),
            ("level", s(self.level.name())),
            ("card", s(self.card)),
            ("title", s(&self.title)),
            ("text", s(&self.text)),
            ("gene", self.gene.clone()),
            ("detail", self.detail.as_deref().map_or(Value::Null, s)),
            ("evidence", Value::Array(self.evidence.clone())),
            ("viewer", self.viewer.as_deref().map_or(Value::Null, s)),
            ("af", self.af.clone()),
            ("section", int(self.level.section())),
            ("steps", card_list("steps")),
            ("caveat", card_list("caveat")),
            ("anchor", s(self.card)),
            ("plan", strings(self.level.plan())),
            ("level_words", s(self.level.words())),
            (
                "related",
                strings(
                    found
                        .iter()
                        .filter(|other| {
                            other.id != self.id
                                && (other.card == self.card
                                    || (truthy(&self.gene) && other.gene == self.gene))
                        })
                        .take(6)
                        .map(|other| &other.id),
                ),
            ),
        ])
    }
}

fn slug(parts: &[&str]) -> String {
    let raw = parts
        .iter()
        .copied()
        .filter(|p| !p.is_empty())
        .collect::<Vec<_>>()
        .join("-");
    let mut id = String::new();
    let mut gap = false;
    for ch in raw.to_lowercase().chars() {
        if ch.is_ascii_alphanumeric() {
            if gap && !id.is_empty() {
                id.push('-');
            }
            id.push(ch);
            gap = false;
        } else {
            gap = true;
        }
    }
    if id.len() > 80 {
        // IDs are bookmarked URLs. Keep the original SHA-1 suffix for long alleles.
        id = format!(
            "{}-{:08x}",
            id[..71].trim_end_matches('-'),
            sha1_prefix(raw.as_bytes())
        );
    }
    if id.is_empty() { "finding".into() } else { id }
}

// Only the first digest word is part of the existing finding ID format.
fn sha1_prefix(bytes: &[u8]) -> u32 {
    let mut message = bytes.to_vec();
    message.push(0x80);
    message.resize((message.len() + 8).div_ceil(64) * 64 - 8, 0);
    message.extend_from_slice(&(bytes.len() as u64 * 8).to_be_bytes());
    let mut digest = [
        0x67452301_u32,
        0xefcdab89,
        0x98badcfe,
        0x10325476,
        0xc3d2e1f0,
    ];
    for block in message.chunks_exact(64) {
        let mut words = [0_u32; 80];
        for (word, bytes) in words.iter_mut().zip(block.chunks_exact(4)) {
            *word = u32::from_be_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]);
        }
        for i in 16..80 {
            words[i] = (words[i - 3] ^ words[i - 8] ^ words[i - 14] ^ words[i - 16]).rotate_left(1);
        }
        let [mut a, mut b, mut c, mut d, mut e] = digest;
        for (i, word) in words.into_iter().enumerate() {
            let (f, k) = match i {
                0..20 => ((b & c) | (!b & d), 0x5a827999),
                20..40 => (b ^ c ^ d, 0x6ed9eba1),
                40..60 => ((b & c) | (b & d) | (c & d), 0x8f1bbcdc),
                _ => (b ^ c ^ d, 0xca62c1d6),
            };
            let next = a
                .rotate_left(5)
                .wrapping_add(f)
                .wrapping_add(e)
                .wrapping_add(k)
                .wrapping_add(word);
            [a, b, c, d, e] = [next, a, b.rotate_left(30), c, d];
        }
        for (word, x) in digest.iter_mut().zip([a, b, c, d, e]) {
            *word = word.wrapping_add(x);
        }
    }
    digest[0]
}

fn error_kind(error: &Malformed) -> &'static str {
    if error.0.contains("missing") {
        "KeyError"
    } else if error.0.contains("not a number") || error.0.contains("not an integer") {
        "ValueError"
    } else {
        "TypeError"
    }
}

fn add_row(
    found: &mut Vec<Finding>,
    card: &'static str,
    level: Level,
    read: impl FnOnce() -> Read<Option<Finding>>,
) {
    match read() {
        Ok(Some(f)) => found.push(f),
        Ok(None) => {}
        Err(err) => found.push(unreadable(card, level, found.len(), error_kind(&err))),
    }
}

fn unreadable(card: &'static str, level: Level, index: usize, kind: &str) -> Finding {
    Finding::new(
        level,
        card,
        &["unreadable", card, &index.to_string()],
        "A row of this section could not be read",
        format!(
            "One row of this section holds a value of a kind the page cannot read ({kind}). The other rows show as usual; the report and its files have the row as the step wrote it."
        ),
    )
}

fn object_row(v: &Value) -> Read<()> {
    if matches!(v, Value::Object(_)) {
        Ok(())
    } else {
        Err(Malformed("row is not an object".to_owned()))
    }
}

fn finding_rows(summary: &Summary) -> Vec<Finding> {
    let mut found = Vec::new();
    let excluded = clinical::clinical_findings(summary, &mut found);
    clinvar::findings(summary, &excluded, &mut found);
    sv::findings(summary, &mut found);
    clinical::family_findings(summary, &mut found);
    pgx::findings(summary, &mut found);
    qc::findings(summary, &mut found);
    let mut seen = HashSet::new();
    found.retain(|f| seen.insert(f.id.clone()));
    found
}

const JOURNAL: [(&str, &str, &str); 10] = [
    ("act", "To act on", "Is there anything to take to a doctor?"),
    (
        "medicines",
        "Medicines",
        "Which medicines should a prescriber know about?",
    ),
    (
        "family",
        "For your family",
        "What could I pass on to a child?",
    ),
    (
        "look",
        "Worth a look",
        "What deserves a conversation, but is not a diagnosis?",
    ),
    (
        "trust",
        "Can I trust this?",
        "This is about the data, not about the person.",
    ),
    ("origins", "Origins", "Where do my ancestors come from?"),
    (
        "traits",
        "Traits and scores",
        "What do studies say about traits and risk scores?",
    ),
    (
        "assessed",
        "Assessed, no finding",
        "What each step looked at, and what it does not see.",
    ),
    (
        "cannot",
        "What it cannot see",
        "What no step of this run looks at.",
    ),
    (
        "sources",
        "Sources and files",
        "The databases, the report and the files behind this page.",
    ),
];

fn journal(cards: &mut [Obj], found: &mut Vec<Finding>) -> (Value, Value) {
    let have: HashSet<_> = found.iter().map(|f| f.card).collect();
    let defs = CHAPTERS.iter().flat_map(|chapter| chapter.sections);
    for (card, sec) in cards.iter().zip(defs) {
        if truthy(&card["flagged"])
            && !have.contains(sec.id)
            && !matches!(str_val(&card["state"]), "missing" | "failed")
        {
            let level = if str_val(&card["chapter"]) == "quality" {
                Level::Trust
            } else {
                Level::Review
            };
            let mut finding = Finding::new(
                level,
                sec.id,
                &["card", sec.id],
                sec.title,
                text_or(
                    &card["headline"],
                    &format!("The report flags this section: {}.", sec.title),
                ),
            );
            finding.evidence = list(&card["facts"])
                .iter()
                .take(4)
                .filter(|v| list(v).get(1) != Some(&s(".")))
                .cloned()
                .collect();
            found.push(finding);
        }
    }
    let mut info = 0;
    for card in cards.iter_mut() {
        let id = str_val(&card["id"]);
        let related: Vec<_> = found.iter().filter(|f| f.card == id).collect();
        let top = related.iter().map(|f| f.level).min();
        let chapter = str_val(&card["chapter"]);
        let home = match chapter {
            "quality" => 5,
            "clinical" => top.map_or(8, Level::section),
            "carrier" => 3,
            "medicines" => 2,
            "ancestry" => 6,
            "traits" => 7,
            _ => 10,
        };
        let level = if matches!(str_val(&card["state"]), "missing" | "failed") {
            None
        } else {
            Some(top.unwrap_or(if chapter == "quality" {
                Level::Trust
            } else if truthy(&card["flagged"]) {
                Level::Review
            } else {
                Level::Info
            }))
        };
        info += usize::from(level == Some(Level::Info));
        card.insert("level".into(), level.map_or(Value::Null, |v| s(v.name())));
        card.insert("section".into(), int(home));
        card.insert("findings".into(), strings(related.iter().map(|f| &f.id)));
    }
    let levels = Level::ALL
        .into_iter()
        .map(|lv| {
            (
                lv.name().to_owned(),
                int(if lv == Level::Info {
                    info
                } else {
                    found.iter().filter(|f| f.level == lv).count()
                }),
            )
        })
        .collect();
    let sections = JOURNAL
        .iter()
        .enumerate()
        .map(|(index, (key, title, question))| {
            let n = index + 1;
            obj([
                ("n", int(n)),
                ("key", s(key)),
                ("title", s(title)),
                ("question", s(question)),
                (
                    "findings",
                    strings(
                        found
                            .iter()
                            .filter(|f| f.level.section() == n)
                            .map(|f| &f.id),
                    ),
                ),
                (
                    "cards",
                    strings(
                        cards
                            .iter()
                            .filter(|c| c["section"] == int(n))
                            .map(|c| str_val(&c["id"])),
                    ),
                ),
            ])
        })
        .collect();
    (Value::Object(levels), Value::Array(sections))
}

fn parsed(summary: &Value) -> Result<Summary, String> {
    if get(summary, "sections")
        .is_some_and(|value| truthy(value) && !matches!(value, Value::Object(_)))
    {
        return Err("expected a summary sections object".to_owned());
    }
    let normalized;
    let input = if get(summary, "sections").is_some() {
        summary
    } else if !truthy(summary) || matches!(summary, Value::Object(_)) {
        let mut root: Obj = fields(summary)
            .map(|(k, v)| (k.to_owned(), v.clone()))
            .collect();
        root.insert("sections".into(), Value::Object(Obj::new()));
        normalized = Value::Object(root);
        &normalized
    } else {
        return Err("expected a summary object".to_owned());
    };
    let json = dumps(input, &Dump::default()).map_err(|e| format!("{e:?}"))?;
    Summary::parse(&json).map_err(|e| e.to_string())
}

/// Build the live site's results JSON from an already read summary.
pub fn build(
    sample: &str,
    summary: &Value,
    files: Option<&BTreeSet<String>>,
) -> Result<Value, String> {
    Ok(build_summary(sample, &parsed(summary)?, files))
}

fn build_summary(sample: &str, summary: &Summary, files: Option<&BTreeSet<String>>) -> Value {
    let mut cards: Vec<_> = CHAPTERS
        .iter()
        .flat_map(|chapter| {
            chapter
                .sections
                .iter()
                .map(|sec| card(sec, chapter.id, summary, files))
        })
        .collect();
    let mut found = finding_rows(summary);
    let (levels, journal) = journal(&mut cards, &mut found);
    let lead = headlines::lead(summary, &found);
    let veil = veil::build(summary, &found, &cards, &lead);
    let findings = found.iter().map(|f| f.json(&cards, &found)).collect();
    let attention = cards.iter().filter(|c| truthy(&c["attention"])).count();
    let mut cards = cards.into_iter();
    let chapters = CHAPTERS
        .iter()
        .map(|chapter| {
            let cards: Vec<_> = cards.by_ref().take(chapter.sections.len()).collect();
            let attention = cards.iter().filter(|c| truthy(&c["attention"])).count();
            obj([
                ("key", s(chapter.id)),
                ("title", s(chapter.title)),
                ("question", s(chapter.question)),
                (
                    "sections",
                    Value::Array(cards.into_iter().map(Value::Object).collect()),
                ),
                ("attention", int(attention)),
            ])
        })
        .collect();
    let root = summary.root();
    let sex = summary.section("sex_check").data;
    let databases = fields(val(root, "databases"))
        .filter(|(_, v)| {
            matches!(
                v,
                Value::Str(_) | Value::Int(_) | Value::Float(_) | Value::Bool(_)
            )
        })
        .map(|(k, v)| (k.to_owned(), s(&plain(&cell(v)))))
        .collect();
    obj([
        ("schema_version", Value::Int(1)),
        ("sample", s(sample)),
        ("generated_utc", val(root, "generated_utc").clone()),
        (
            "person",
            obj([
                ("sample", s(sample)),
                (
                    "declared_sex",
                    val(val(root, "run"), "declared_sex").clone(),
                ),
                ("inferred_sex", val(sex, "inferred_sex").clone()),
                ("karyotype", val(sex, "karyotype").clone()),
            ]),
        ),
        ("chapters", Value::Array(chapters)),
        ("findings", Value::Array(findings)),
        ("levels", levels),
        ("journal", journal),
        (
            "lead",
            Value::Array(lead.into_iter().map(|p| p.json()).collect()),
        ),
        (
            "level_words",
            Value::Object(
                Level::ALL
                    .into_iter()
                    .map(|lv| (lv.name().into(), s(lv.words())))
                    .collect(),
            ),
        ),
        ("databases", Value::Object(databases)),
        ("attention", int(attention)),
        (
            "not_assessed",
            strings(
                list(val(root, "not_assessed"))
                    .iter()
                    .map(|v| plain(&text(v))),
            ),
        ),
        (
            "not_run",
            strings(list(val(root, "not_run")).iter().map(|v| plain(&text(v)))),
        ),
        (
            "new_since",
            get(root, "new_since")
                .filter(|v| truthy(v))
                .cloned()
                .unwrap_or(Value::Null),
        ),
        (
            "charts",
            obj([
                ("clinvar", s("clinvar")),
                ("genome", s("clinvar")),
                ("mito", s("mito")),
                ("ancestry", s("ancestry-fine")),
                ("prs", s("prs")),
            ]),
        ),
        ("veil", veil),
    ])
}

/// Build the live home page's counts and its open and closed lead.
pub fn overview(summary: &Value) -> Result<Value, String> {
    let summary = parsed(summary)?;
    Ok(headlines::overview(
        &summary,
        &build_summary("", &summary, None),
    ))
}

/// The results a report's search index needs to tag.
pub fn closable_of(root: &Value) -> Vec<og_render::index::Closable> {
    parsed(root).map_or_else(
        |_| Vec::new(),
        |summary| veil::spans(&finding_rows(&summary)),
    )
}

#[derive(Default)]
struct Args {
    summary: String,
    sample: Option<String>,
    files: BTreeSet<String>,
    out: Option<String>,
    overview: bool,
    spans: bool,
    tag_search: bool,
}

fn args(argv: &[String]) -> Result<Args, String> {
    let mut args = Args::default();
    let mut argv = argv.iter();
    while let Some(arg) = argv.next() {
        let next = |argv: &mut std::slice::Iter<'_, String>| {
            argv.next()
                .cloned()
                .ok_or_else(|| format!("{arg} needs a value"))
        };
        match arg.as_str() {
            "--overview" => args.overview = true,
            "--spans" => args.spans = true,
            "--tag-search" => args.tag_search = true,
            "--sample" => args.sample = Some(next(&mut argv)?),
            "--file" => {
                args.files.insert(next(&mut argv)?);
            }
            "--out" => args.out = Some(next(&mut argv)?),
            "--summary" if args.summary.is_empty() => args.summary = next(&mut argv)?,
            path if (!path.starts_with('-') || path == "-") && args.summary.is_empty() => {
                args.summary.clone_from(arg)
            }
            _ => return Err(format!("unexpected argument: {arg}")),
        }
    }
    if args.summary.is_empty() {
        return Err("a summary is required".into());
    }
    if [args.overview, args.spans, args.tag_search]
        .into_iter()
        .filter(|flag| *flag)
        .count()
        > 1
    {
        return Err("--overview, --spans and --tag-search are mutually exclusive".into());
    }
    Ok(args)
}

/// Dispatch the live results verb. The summary may be read from stdin with `-`.
pub fn cli(argv: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    if argv
        .iter()
        .any(|arg| matches!(arg.as_str(), "-h" | "--help"))
    {
        let _ = writeln!(out, "{USAGE}");
        return 0;
    }
    let args = match args(argv) {
        Ok(args) => args,
        Err(why) => {
            let _ = writeln!(err, "{USAGE}\nog live results: {why}");
            return 2;
        }
    };
    let run = || -> Result<(), String> {
        let text = if args.summary == "-" {
            std::io::read_to_string(std::io::stdin())
        } else {
            std::fs::read_to_string(&args.summary)
        }
        .map_err(|e| e.to_string())?;
        let root = pycompat::json::loads(&text).map_err(|e| e.to_string())?;
        let result = if args.tag_search {
            let summary = parsed(val(&root, "summary"))?;
            let mut terms = val(&root, "terms").clone();
            veil::tag_search(&mut terms, &veil::spans(&finding_rows(&summary)))?;
            terms
        } else {
            let summary = parsed(&root)?;
            if args.spans {
                veil::spans_json(&finding_rows(&summary))
            } else if args.overview {
                headlines::overview(&summary, &build_summary("", &summary, None))
            } else {
                let sample = args
                    .sample
                    .as_deref()
                    .unwrap_or_else(|| string(summary.root(), "sample"));
                build_summary(sample, &summary, Some(&args.files))
            }
        };
        let body = dumps(
            &result,
            &Dump {
                indent: Some(2),
                ..Dump::default()
            },
        )
        .map_err(|e| format!("{e:?}"))?;
        if let Some(path) = args.out {
            let mut file = og_io::AtomicFile::create(path).map_err(|e| e.to_string())?;
            writeln!(file, "{body}")
                .and_then(|_| file.commit())
                .map_err(|e| e.to_string())
        } else {
            writeln!(out, "{body}").map_err(|e| e.to_string())
        }
    };
    match run() {
        Ok(()) => 0,
        Err(why) => {
            let _ = writeln!(err, "og live results: {why}");
            1
        }
    }
}
