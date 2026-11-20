//! `og report charts`: the report's five charts, as data and as inline SVG,
//! built once from summary.json for both places a chart is shown: the HTML
//! report puts each SVG in its section, under that section's caveat
//! ([`figure`]), and the live page's charts plugin reads every spec with its
//! SVG from `<sample>_charts.json` ([`json`]).
//!
//! ```text
//! og report charts --summary IN.json [--charts-json OUT.json]
//! ```
//!
//! It writes what `bin/render_report.py --summary IN.json --charts-json
//! OUT.json` writes, and prints the same last line (the summary's sections
//! counted by state), so the two run with the same arguments.
//!
//! The five charts, each drawn only when its section ran:
//!
//! - `clinvar-stars` (bars): ClinVar hits per review star, 4 to 0, and a
//!   grey bar for the low-confidence calls, which the star counts leave out.
//!   Yellow and grey only, never green; with no hit a grey frame says that
//!   no hit is not an all-clear.
//! - `genome` (ideogram), drawn in the ROH section so only when ROH ran:
//!   chromosomes 1 to 22, X and Y of GRCh38 with the runs of homozygosity of
//!   5 MB or more, a tick for each ClinVar hit, an outlined span for each
//!   AnnotSV class 4 or 5 SV and a bar under the chromosome for each MoChA
//!   call. What has no place on them is counted in `not_drawn`.
//! - `mito` (scatter): the heteroplasmic chrM sites by position and allele
//!   fraction, the D-loop poly-C length variants marked apart.
//! - `ancestry-fine` (ranges): each fine-grained ancestry group the report
//!   lists, its share as a dot and its bootstrap range as a line.
//! - `prs` (strip): PRS percentiles on a 0 to 100 strip, only when pgsc_calc
//!   had an ancestry panel. Without one a raw score has no scale, so the
//!   chart is the report's own note (`text`); a score of the other sex, or
//!   one with no percentile, says "no percentile" and is never placed.
//!
//! Every caveat comes from [`crate::notes`] by its name, so the chart and the
//! report say the same thing. The SVG has no script, no style attribute and
//! no outside reference: colours are presentation attributes, links are
//! `<a href="#anchor">` inside the page and each mark carries a `<title>` and
//! `data-mark="<id>"`.

mod build;
mod svg;

use std::io::Write;

use indexmap::IndexMap;
use pycompat::json::{Dump, Value, dumps};

use crate::model::{Malformed, Summary, need, py};

pub use build::build;
pub use svg::draw;

/// The `schema_version` of `<sample>_charts.json`.
pub const SCHEMA_VERSION: i128 = 1;

/// GRCh38 primary assembly lengths (UCSC hg38 chrom.sizes), in the order the
/// genome chart draws them. The live site's family page uses them too.
pub const GRCH38: [(&str, i128); 24] = [
    ("chr1", 248_956_422),
    ("chr2", 242_193_529),
    ("chr3", 198_295_559),
    ("chr4", 190_214_555),
    ("chr5", 181_538_259),
    ("chr6", 170_805_979),
    ("chr7", 159_345_973),
    ("chr8", 145_138_636),
    ("chr9", 138_394_717),
    ("chr10", 133_797_422),
    ("chr11", 135_086_622),
    ("chr12", 133_275_309),
    ("chr13", 114_364_328),
    ("chr14", 107_043_718),
    ("chr15", 101_991_189),
    ("chr16", 90_338_345),
    ("chr17", 83_257_441),
    ("chr18", 80_373_285),
    ("chr19", 58_617_616),
    ("chr20", 64_444_167),
    ("chr21", 46_709_983),
    ("chr22", 50_818_468),
    ("chrX", 156_040_895),
    ("chrY", 57_227_415),
];

/// The length of a GRCh38 primary chromosome, `None` for any other name.
pub fn chrom_len(chrom: &str) -> Option<i128> {
    GRCH38.iter().find(|(c, _)| *c == chrom).map(|(_, n)| *n)
}

/// The length of chrM.
pub const MITO_LEN: i128 = 16569;

/// The shortest run of homozygosity the genome chart draws, in MB: the runs
/// F_ROH counts (the roh section's own floor).
pub const ROH_MIN_MB: f64 = og_summary::collect::group_c::roh::ROH_MIN_BP / 1e6;

/// What a mark's colour means: yellow for something to look at, grey for
/// what is not counted, blue for plain information. No green anywhere: a
/// chart is not an all-clear.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Tone {
    Attention,
    Muted,
    Info,
}

impl Tone {
    pub fn name(self) -> &'static str {
        match self {
            Tone::Attention => "attention",
            Tone::Muted => "muted",
            Tone::Info => "info",
        }
    }

    pub fn colour(self) -> &'static str {
        match self {
            Tone::Attention => "#e0b94a",
            Tone::Muted => "#8a909b",
            Tone::Info => "#5b9ee6",
        }
    }
}

/// How a chart is drawn.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kind {
    Bars,
    Ranges,
    Strip,
    Text,
    Ideogram,
    Scatter,
}

impl Kind {
    pub fn name(self) -> &'static str {
        match self {
            Kind::Bars => "bars",
            Kind::Ranges => "ranges",
            Kind::Strip => "strip",
            Kind::Text => "text",
            Kind::Ideogram => "ideogram",
            Kind::Scatter => "scatter",
        }
    }
}

/// What a mark of the genome or mito chart stands for.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MarkKind {
    Roh,
    Hit,
    Sv,
    Mca,
    Site,
    Polyc,
}

impl MarkKind {
    pub fn name(self) -> &'static str {
        match self {
            MarkKind::Roh => "roh",
            MarkKind::Hit => "hit",
            MarkKind::Sv => "sv",
            MarkKind::Mca => "mca",
            MarkKind::Site => "site",
            MarkKind::Polyc => "polyc",
        }
    }
}

/// One mark of a chart. `value`, `low` and `high` keep the summary's own
/// type (an int stays an int); `null` where they do not apply.
#[derive(Debug, Clone, PartialEq)]
pub struct Mark {
    pub id: String,
    pub label: String,
    pub value: Value,
    pub low: Value,
    pub high: Value,
    pub chrom: Option<String>,
    pub pos: Option<i128>,
    pub end: Option<i128>,
    pub tone: Tone,
    /// The id of the report element the mark links to.
    pub anchor: Option<String>,
    /// `chrom:pos` or `chrom:start-end`, for the genome viewer.
    pub locus: Option<String>,
    pub detail: String,
    /// The heading of the ancestry level the mark belongs to.
    pub group: Option<String>,
    pub kind: Option<MarkKind>,
    /// A genome chart mark of another section carries that section's caveat.
    pub caveat: Option<String>,
}

impl Mark {
    /// A mark with nothing but its id, label and value; tone `info`.
    pub fn new(id: String, label: String, value: Value) -> Mark {
        Mark {
            id,
            label,
            value,
            low: Value::Null,
            high: Value::Null,
            chrom: None,
            pos: None,
            end: None,
            tone: Tone::Info,
            anchor: None,
            locus: None,
            detail: String::new(),
            group: None,
            kind: None,
            caveat: None,
        }
    }

    fn to_value(&self) -> Value {
        let s = |x: &Option<String>| x.clone().map_or(Value::Null, Value::Str);
        let i = |x: Option<i128>| x.map_or(Value::Null, Value::Int);
        let mut m = IndexMap::new();
        m.insert("id".into(), Value::Str(self.id.clone()));
        m.insert("label".into(), Value::Str(self.label.clone()));
        m.insert("value".into(), self.value.clone());
        m.insert("low".into(), self.low.clone());
        m.insert("high".into(), self.high.clone());
        m.insert("chrom".into(), s(&self.chrom));
        m.insert("pos".into(), i(self.pos));
        m.insert("end".into(), i(self.end));
        m.insert("tone".into(), Value::Str(self.tone.name().into()));
        m.insert("anchor".into(), s(&self.anchor));
        m.insert("locus".into(), s(&self.locus));
        m.insert("detail".into(), Value::Str(self.detail.clone()));
        if let Some(g) = &self.group {
            m.insert("group".into(), Value::Str(g.clone()));
        }
        if let Some(k) = self.kind {
            m.insert("kind".into(), Value::Str(k.name().into()));
        }
        if let Some(c) = &self.caveat {
            m.insert("caveat".into(), Value::Str(c.clone()));
        }
        Value::Object(m)
    }
}

/// Something a chart has no place for, and how many.
#[derive(Debug, Clone, PartialEq)]
pub struct NotDrawn {
    pub what: String,
    pub count: i128,
}

/// One chart: its data and its SVG.
#[derive(Debug, Clone, PartialEq)]
pub struct Spec {
    pub id: &'static str,
    /// The summary section whose data it shows.
    pub section: &'static str,
    pub kind: Kind,
    pub title: &'static str,
    /// The name of its caveat in [`crate::notes`], and that caveat's text.
    pub caveat_key: &'static str,
    pub caveat: &'static str,
    pub axis: IndexMap<String, Value>,
    pub marks: Vec<Mark>,
    pub not_drawn: Vec<NotDrawn>,
    /// The text a `text` chart shows.
    pub text: Option<String>,
    /// The chart drawn ([`draw`]).
    pub svg: String,
}

impl Spec {
    /// The spec as `<sample>_charts.json` holds it.
    pub fn to_value(&self) -> Value {
        let mut m = IndexMap::new();
        m.insert("id".into(), Value::Str(self.id.into()));
        m.insert("section".into(), Value::Str(self.section.into()));
        m.insert("kind".into(), Value::Str(self.kind.name().into()));
        m.insert("title".into(), Value::Str(self.title.into()));
        m.insert("caveat_key".into(), Value::Str(self.caveat_key.into()));
        m.insert("caveat".into(), Value::Str(self.caveat.into()));
        m.insert("axis".into(), Value::Object(self.axis.clone()));
        m.insert(
            "marks".into(),
            Value::Array(self.marks.iter().map(Mark::to_value).collect()),
        );
        let nd = self.not_drawn.iter().map(|x| {
            let mut o = IndexMap::new();
            o.insert("what".into(), Value::Str(x.what.clone()));
            o.insert("count".into(), Value::Int(x.count));
            Value::Object(o)
        });
        m.insert("not_drawn".into(), Value::Array(nd.collect()));
        if let Some(t) = &self.text {
            m.insert("text".into(), Value::Str(t.clone()));
        }
        m.insert("svg".into(), Value::Str(self.svg.clone()));
        Value::Object(m)
    }
}

/// The HTML report's line for a chart, which goes into its section right
/// under its caveat; `None` for a `text` chart, which is the note the card
/// already prints (the live page shows it).
pub fn figure(spec: &Spec) -> Option<String> {
    (spec.kind != Kind::Text).then(|| format!("    <figure class=\"chart\">{}</figure>", spec.svg))
}

/// `<sample>_charts.json`: the schema version, the sample, when the summary
/// was made and the charts.
pub fn json(s: &Summary) -> Result<Value, Malformed> {
    let mut m = IndexMap::new();
    m.insert("schema_version".into(), Value::Int(SCHEMA_VERSION));
    m.insert("sample".into(), s.need("sample")?.clone());
    m.insert("generated_utc".into(), s.need("generated_utc")?.clone());
    let charts = build(s)?.iter().map(Spec::to_value).collect();
    m.insert("charts".into(), Value::Array(charts));
    Ok(Value::Object(m))
}

const USAGE: &str = "usage: og report charts --summary IN.json [--charts-json OUT.json]";

/// The arguments: the summary, and where to write the charts.
fn parse(args: &[String]) -> Result<(String, Option<String>), String> {
    let (mut summary, mut out) = (None, None);
    let mut it = args.iter();
    while let Some(arg) = it.next() {
        let (name, inline) = match arg.split_once('=') {
            Some((n, v)) if n.starts_with("--") => (n, Some(v.to_owned())),
            _ => (arg.as_str(), None),
        };
        let mut value = || {
            inline
                .clone()
                .or_else(|| it.next().cloned())
                .ok_or_else(|| format!("argument {name}: expected one argument"))
        };
        match name {
            "--summary" => summary = Some(value()?),
            "--charts-json" => out = Some(value()?),
            n => return Err(format!("unrecognized arguments: {n}")),
        }
    }
    Ok((summary.ok_or("give --summary")?, out))
}

/// `Sections: missing 28, ok 20`: the sections counted by state, as `og
/// report render` and `bin/render_report.py` end.
fn state_counts(s: &Summary) -> Result<String, Malformed> {
    let Some(Value::Object(secs)) = s.get("sections") else {
        return Err(Malformed("summary.sections is not an object".into()));
    };
    let mut counts: Vec<(String, usize)> = Vec::new();
    for (key, sec) in secs {
        let state = py(need(sec, "state", &format!("sections.{key}"))?)?;
        match counts.iter_mut().find(|(k, _)| *k == state) {
            Some((_, n)) => *n += 1,
            None => counts.push((state, 1)),
        }
    }
    counts.sort();
    let parts: Vec<String> = counts.iter().map(|(k, n)| format!("{k} {n}")).collect();
    Ok(format!("Sections: {}", parts.join(", ")))
}

fn write(path: &str, body: &str) -> std::io::Result<()> {
    let mut f = og_io::AtomicFile::create(path)?;
    f.write_all(body.as_bytes())?;
    f.commit()
}

/// `og report charts ARGS`.
pub fn cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    if args.iter().any(|a| a == "-h" || a == "--help") {
        let _ = writeln!(out, "{USAGE}");
        return 0;
    }
    let (path, dest) = match parse(args) {
        Ok(a) => a,
        Err(why) => {
            let _ = writeln!(err, "{USAGE}\nog report charts: error: {why}");
            return 2;
        }
    };
    let run = || -> Result<String, String> {
        let text = std::fs::read_to_string(&path).map_err(|e| format!("{path}: {e}"))?;
        let s = Summary::parse(&text).map_err(|e| format!("{path}: {e}"))?;
        if let Some(dest) = &dest {
            let v = json(&s).map_err(|e| format!("{path}: {e}"))?;
            let opts = Dump {
                indent: Some(1),
                ..Dump::default()
            };
            let body = dumps(&v, &opts).map_err(|e| format!("{path}: {e:?}"))? + "\n";
            write(dest, &body).map_err(|e| format!("{dest}: {e}"))?;
        }
        state_counts(&s).map_err(|e| format!("{path}: {e}"))
    };
    match run() {
        Ok(line) => {
            let _ = writeln!(out, "{line}");
            0
        }
        Err(why) => {
            let _ = writeln!(err, "og report charts: {why}");
            1
        }
    }
}
