//! The map document, version 1: columnar JSON in one gzip member.
//!
//! - The JSON is `json.dumps(doc, separators=(",", ":"))`, keys in the
//!   order of the format (`README.md`).
//! - The gzip member has mtime 0 and no file name, so the same input always
//!   gives the same bytes.
//! - Per chromosome, rows sorted by start: starts as gaps from the previous
//!   row's start (`ds`, the first one absolute) and ends as lengths (`l`).

use crate::layers::roh::{self, ART_NAMES};
use crate::reference::CHROMS;
use crate::text::{self, Piece};
use crate::veil;
use crate::{Layer, Model};
use flate2::{Compression, GzBuilder};
use indexmap::IndexMap;
use pycompat::json::{self, Dump, Value};
use std::io::{self, Write};

/// The bin size of the density columns, in bases.
pub const BIN: u64 = 1_000_000;

/// An object with keys in the order given.
pub fn obj<const N: usize>(pairs: [(&str, Value); N]) -> Value {
    Value::Object(pairs.into_iter().map(|(k, v)| (k.to_owned(), v)).collect())
}

fn int(x: impl Into<i128>) -> Value {
    Value::Int(x.into())
}

fn str(s: &str) -> Value {
    Value::Str(s.to_owned())
}

fn ints<T: Into<i128> + Copy>(xs: impl IntoIterator<Item = T>) -> Value {
    Value::Array(xs.into_iter().map(|x| Value::Int(x.into())).collect())
}

/// The gap column of sorted starts: the first start, then each start minus
/// the one before it.
pub fn gaps(starts: impl IntoIterator<Item = u64>) -> Vec<u64> {
    let mut prev = 0;
    starts
        .into_iter()
        .map(|s| {
            let d = s - prev;
            prev = s;
            d
        })
        .collect()
}

fn pieces(p: Vec<Piece>) -> Value {
    Value::Array(
        p.into_iter()
            .map(|x| obj([("t", Value::Str(x.t)), ("b", Value::Bool(x.b))]))
            .collect(),
    )
}

fn layer(l: &Layer) -> Value {
    let mut m: IndexMap<String, Value> = IndexMap::new();
    m.insert("state".into(), str(l.state.word()));
    m.insert("step".into(), str(l.step));
    m.insert("step_name".into(), str(l.step_name));
    m.insert("source".into(), str(&l.source));
    if let Some(n) = &l.note {
        m.insert("note".into(), str(n));
    }
    m.insert("count".into(), int(l.count as u64));
    if let Some((all, kept)) = l.totals {
        m.insert(
            "totals".into(),
            obj([("all", int(all as u64)), ("kept", int(kept as u64))]),
        );
    }
    m.insert("caveat".into(), str(&l.caveat));
    if let Some(e) = &l.empty {
        m.insert("empty".into(), str(e));
    }
    Value::Object(m)
}

fn reference(m: &Model) -> Value {
    let chroms = m
        .reference
        .chroms
        .iter()
        .map(|c| {
            obj([
                ("n", str(c.name)),
                ("len", int(c.len)),
                ("cen", ints([c.cen.0, c.cen.1])),
                (
                    "masked",
                    Value::Array(c.masked.iter().map(|&(s, e)| ints([s, e])).collect()),
                ),
            ])
        })
        .collect();
    let bands = m
        .reference
        .chroms
        .iter()
        .map(|c| {
            let b = &c.bands;
            let col = obj([
                ("ds", ints(gaps(b.iter().map(|x| x.s)))),
                ("l", ints(b.iter().map(|x| x.e - x.s))),
                (
                    "name",
                    Value::Array(b.iter().map(|x| str(&x.name)).collect()),
                ),
                (
                    "stain",
                    Value::Array(b.iter().map(|x| str(&x.stain)).collect()),
                ),
            ]);
            (c.name.to_owned(), col)
        })
        .collect();
    obj([
        ("chroms", Value::Array(chroms)),
        ("bands", Value::Object(bands)),
        ("gencode", m.gencode.as_deref().map_or(Value::Null, str)),
    ])
}

/// The person's ROH columns, its bins and the sentences of its long runs.
fn roh_columns(m: &Model) -> (Value, Value, Vec<Piece>, IndexMap<String, Value>) {
    let rows = &m.roh.rows;
    let mut by_chrom: IndexMap<String, Value> = IndexMap::new();
    let mut chrom_text: IndexMap<String, Value> = IndexMap::new();
    let (mut n_long, mut bp_long) = (0usize, 0u64);
    for part in rows.chunk_by(|a, b| a.chrom == b.chrom) {
        let name = m.reference.chroms[part[0].chrom].name;
        by_chrom.insert(
            name.to_owned(),
            obj([
                ("ds", ints(gaps(part.iter().map(|r| r.s)))),
                ("l", ints(part.iter().map(|r| r.e - r.s))),
                ("markers", ints(part.iter().map(|r| r.markers))),
                (
                    "q",
                    Value::Array(part.iter().map(|r| Value::Float(r.q)).collect()),
                ),
                ("long", ints(part.iter().map(|r| u8::from(r.long)))),
                ("art", ints(part.iter().map(|r| r.art as u8))),
            ]),
        );
        let long: Vec<u64> = part.iter().filter(|r| r.long).map(|r| r.e - r.s).collect();
        if !long.is_empty() {
            let bp = long.iter().sum();
            n_long += long.len();
            bp_long += bp;
            chrom_text.insert(
                name.to_owned(),
                pieces(text::roh_long_chrom(long.len(), bp)),
            );
        }
    }
    let mut bins: IndexMap<String, Value> = IndexMap::new();
    bins.insert("size".into(), int(BIN));
    for (c, b) in roh::bins(rows, &m.reference, BIN) {
        bins.insert(m.reference.chroms[c].name.to_owned(), ints(b));
    }
    let genome = if m.roh.layer.state == crate::State::Ok {
        text::roh_long(n_long, bp_long)
    } else {
        Vec::new()
    };
    (
        Value::Object(by_chrom),
        Value::Object(bins),
        genome,
        chrom_text,
    )
}

fn opt(s: Option<&str>) -> Value {
    s.map_or(Value::Null, str)
}

fn strs<'a>(xs: impl IntoIterator<Item = &'a str>) -> Value {
    Value::Array(xs.into_iter().map(str).collect())
}

/// The person's findings on the map: one flat column per field, rows in
/// map order. `veil` is the index in `dict.veil` of the first closable
/// finding whose span the row overlaps, else -1; `veils` every such index.
fn findings_columns(m: &Model) -> Value {
    let rows = &m.findings.rows;
    let tags = m.findings.veil.tags(rows);
    obj([
        ("id", strs(rows.iter().map(|r| r.finding.id.as_str()))),
        ("level", strs(rows.iter().map(|r| r.finding.level.as_str()))),
        (
            "gate",
            Value::Array(
                rows.iter()
                    .map(|r| opt(r.finding.gate.as_deref()))
                    .collect(),
            ),
        ),
        ("chrom", strs(rows.iter().map(|r| CHROMS[r.chrom]))),
        (
            "pos",
            Value::Array(
                rows.iter()
                    .map(|r| r.pos.map_or(Value::Null, int))
                    .collect(),
            ),
        ),
        ("s", ints(rows.iter().map(|r| r.s))),
        ("e", ints(rows.iter().map(|r| r.e))),
        (
            "gene",
            Value::Array(
                rows.iter()
                    .map(|r| opt(r.finding.gene.as_deref()))
                    .collect(),
            ),
        ),
        ("title", strs(rows.iter().map(|r| r.finding.title.as_str()))),
        (
            "hgvs",
            Value::Array(
                rows.iter()
                    .map(|r| opt(r.finding.hgvs.as_deref()))
                    .collect(),
            ),
        ),
        ("zyg", strs(rows.iter().map(|r| r.zyg))),
        // A copy (which parent's) needs the family: phase 2.
        ("copy", Value::Array(vec![Value::Null; rows.len()])),
        (
            "steps",
            Value::Array(
                rows.iter()
                    .map(|r| strs(r.finding.steps.iter().map(String::as_str)))
                    .collect(),
            ),
        ),
        (
            "href",
            Value::Array(rows.iter().map(|r| Value::Str(r.href())).collect()),
        ),
        (
            "reads",
            Value::Array(rows.iter().map(|r| Value::Str(r.reads())).collect()),
        ),
        (
            "phrase",
            Value::Array(rows.iter().map(|r| Value::Str(r.phrase())).collect()),
        ),
        ("veil", ints(veil::column(&tags))),
        (
            "veils",
            Value::Array(
                m.findings
                    .veil
                    .all_tags(rows)
                    .into_iter()
                    .map(ints)
                    .collect(),
            ),
        ),
    ])
}

/// The document as a JSON value.
pub fn to_value(m: &Model) -> Value {
    let (roh_cols, roh_bins, genome, chrom_text) = roh_columns(m);
    let has_y = match m.sex {
        crate::Sex::Male => Value::Bool(true),
        crate::Sex::Female => Value::Bool(false),
        crate::Sex::Unknown => Value::Null,
    };
    let person = obj([
        ("id", str(&m.run)),
        ("label", str(&m.run)),
        ("sex", str(m.sex.word())),
        ("has_y", has_y),
    ]);
    let mut per: IndexMap<String, Value> = IndexMap::new();
    per.insert(
        m.run.clone(),
        obj([
            ("findings", findings_columns(m)),
            ("roh", roh_cols),
            ("roh_bins", roh_bins),
        ]),
    );
    let findings_text = Value::Object(
        text::findings_count()
            .into_iter()
            .map(|(k, p)| (k.to_owned(), pieces(p)))
            .collect(),
    );
    obj([
        ("v", int(crate::VERSION)),
        ("build", str("GRCh38")),
        ("kind", str("person")),
        ("run", str(&m.run)),
        ("stamp", str(&m.stamp)),
        ("made_by", str(&m.made_by)),
        ("people", Value::Array(vec![person])),
        ("ref", reference(m)),
        (
            "layers",
            obj([
                ("findings", layer(&m.findings.layer)),
                ("roh", layer(&m.roh.layer)),
            ]),
        ),
        ("per", Value::Object(per)),
        (
            "dict",
            obj([
                ("art", strs(ART_NAMES)),
                ("veil", strs(m.findings.veil.ids.iter().map(String::as_str))),
            ]),
        ),
        (
            "text",
            obj([
                ("genome", pieces(genome)),
                ("chrom", Value::Object(chrom_text)),
                ("findings", findings_text),
            ]),
        ),
    ])
}

/// The document's bytes: the JSON in one gzip member with mtime 0.
pub fn write(m: &Model) -> io::Result<Vec<u8>> {
    let opts = Dump {
        separators: Some((",", ":")),
        ..Dump::default()
    };
    let text = json::dumps(&to_value(m), &opts)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    let mut gz = GzBuilder::new()
        .mtime(0)
        .write(Vec::with_capacity(text.len() / 4), Compression::best());
    gz.write_all(text.as_bytes())?;
    gz.finish()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gap_column() {
        assert_eq!(gaps([5, 5, 9, 20]), [5, 0, 4, 11]);
        assert_eq!(gaps([]), Vec::<u64>::new());
    }
}
