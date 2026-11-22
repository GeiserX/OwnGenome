//! Drawing data for the live cards. Numeric values are finite or null.

use crate::results::{list, str_val, string, text_or, val};
use crate::state::{obj, s};
use og_render::model::{num, truthy};
use pycompat::json::Value;

fn number(v: &Value) -> Option<f64> {
    num(Some(v))
        .or_else(|| str_val(v).trim().parse().ok())
        .filter(|n| n.is_finite())
}

fn nullable(v: &Value) -> Value {
    number(v).map_or(Value::Null, Value::Float)
}

pub(crate) fn str_gauge(loci: impl IntoIterator<Item = (Value, String)>) -> Value {
    let loci = loci
        .into_iter()
        .map(|(locus, status)| {
            obj([
                ("locus", val(&locus, "locus").clone()),
                ("repeat_count", val(&locus, "repeat_count").clone()),
                ("normal_max", val(&locus, "normal_max").clone()),
                ("pathologic_min", val(&locus, "pathologic_min").clone()),
                ("ci", val(&locus, "ci").clone()),
                ("allele_sizing", val(&locus, "allele_sizing").clone()),
                ("status", s(&status)),
            ])
        })
        .collect();
    obj([("kind", s("str")), ("loci", Value::Array(loci))])
}

pub(crate) fn lpa(d: &Value) -> Value {
    if string(d, "status") != "estimated" {
        return Value::Null;
    }
    let (Some(value), Some(low), Some(high)) = (
        number(val(d, "kiv2_copies")),
        number(val(d, "ci95_low")),
        number(val(d, "ci95_high")),
    ) else {
        return Value::Null;
    };
    obj([
        ("kind", s("ci")),
        ("label", s("KIV-2 repeats, both chromosomes")),
        ("value", Value::Float(value)),
        ("low", Value::Float(low)),
        ("high", Value::Float(high)),
        ("unit", s("repeats")),
        (
            "max",
            Value::Int(60.max(((high * 1.25 / 10.0 + 1.0) as i128).saturating_mul(10))),
        ),
        ("class", val(d, "class").clone()),
    ])
}

pub(crate) fn chip(genes: &[Value], calls: &[Value]) -> Value {
    if genes.is_empty() && calls.is_empty() {
        return Value::Null;
    }
    let rendered_genes = genes
        .iter()
        .map(|g| {
            obj([
                ("gene", val(g, "gene").clone()),
                ("limit", nullable(val(g, "detection_limit_vaf"))),
                ("depth", nullable(val(g, "mean_depth"))),
                ("status", val(g, "status").clone()),
            ])
        })
        .collect();
    let calls = calls
        .iter()
        .filter_map(|c| {
            let vaf = number(val(c, "vaf"))?;
            let gene = val(c, "gene");
            let limit = genes
                .iter()
                .rev()
                .find(|g| val(g, "gene") == gene)
                .map_or(Value::Null, |g| nullable(val(g, "detection_limit_vaf")));
            Some(obj([
                ("gene", gene.clone()),
                ("change", val(c, "protein_change").clone()),
                ("vaf", Value::Float(vaf)),
                ("limit", limit),
            ]))
        })
        .collect();
    obj([
        ("kind", s("vaf")),
        ("max", Value::Float(0.5)),
        ("genes", Value::Array(rendered_genes)),
        ("calls", Value::Array(calls)),
    ])
}

fn copies_value(v: &Value) -> Value {
    number(v)
        .filter(|n| n.fract() == 0.0 && (0.0..=20.0).contains(n))
        .map_or(Value::Null, |n| Value::Int(n as i128))
}

pub(crate) fn copies(genes: &[Value], parascopy: &Value) -> Value {
    const GENES: [&str; 5] = ["SMN1", "SMN2", "HBA1", "HBA2", "CYP21A2"];
    let mut dots: Vec<_> = genes
        .iter()
        .filter(|g| {
            let name = string(g, "gene").split(' ').next().unwrap_or_default();
            GENES.contains(&name)
        })
        .map(|g| {
            obj([
                ("gene", val(g, "gene").clone()),
                ("copies", copies_value(val(g, "copy_number"))),
                ("status", s(&text_or(val(g, "status"), "."))),
            ])
        })
        .collect();
    if truthy(val(parascopy, "rows")) {
        for gene in ["smn1", "smn2"] {
            let x = val(parascopy, gene);
            if !truthy(x) {
                continue;
            }
            let copies = copies_value(val(x, "cn"));
            let status = if matches!(copies, Value::Null) {
                text_or(val(x, "text"), "not called")
            } else {
                "called".to_owned()
            };
            dots.push(obj([
                (
                    "gene",
                    s(&format!("{} (Parascopy)", gene.to_ascii_uppercase())),
                ),
                ("copies", copies),
                ("status", s(&status)),
            ]));
        }
    }
    if dots.is_empty() {
        Value::Null
    } else {
        obj([("kind", s("copies")), ("genes", Value::Array(dots))])
    }
}

pub(crate) fn y_path(tree: &Value) -> Value {
    let path = list(val(tree, "path"));
    if path.is_empty() {
        return Value::Null;
    }
    obj([
        ("kind", s("ypath")),
        ("tree", val(tree, "tree").clone()),
        ("tip", val(tree, "haplogroup").clone()),
        (
            "path",
            Value::Array(
                path.iter()
                    .map(|p| {
                        obj([
                            ("level", val(p, "Level").clone()),
                            ("branch", val(p, "Branch").clone()),
                            ("derived", nullable(val(p, "Derived"))),
                            ("ancestral", nullable(val(p, "Ancestral"))),
                            ("no_call", nullable(val(p, "No_call"))),
                        ])
                    })
                    .collect(),
            ),
        ),
    ])
}
