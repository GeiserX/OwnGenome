use super::*;
use crate::widgets;
use og_render::sections::{ancestry_fine, archaic, haplogroups, roh};
use pycompat::fmt::fixed;

fn share(v: &Value) -> String {
    num(Some(v)).map_or_else(|| ".".to_owned(), |n| format!("{}%", fixed(n * 100.0, 0)))
}

fn explore(d: &Value) -> Read<Value> {
    Ok(Value::Array(
        array(d, "explore")?
            .iter()
            .map(|x| {
                obj([
                    ("name", val(x, "name").clone()),
                    ("url", val(x, "url").clone()),
                ])
            })
            .collect(),
    ))
}

// Historical summaries may hold raw SVG instead of the base64 used by the collector.
fn base64(bytes: &[u8]) -> String {
    let alphabet = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let mut result = String::with_capacity(bytes.len().div_ceil(3) * 4);
    for chunk in bytes.chunks(3) {
        let first = chunk[0];
        let second = chunk.get(1).copied().unwrap_or(0);
        let third = chunk.get(2).copied().unwrap_or(0);
        result.push(char::from(alphabet[usize::from(first >> 2)]));
        result.push(char::from(
            alphabet[usize::from((first & 3) << 4 | second >> 4)],
        ));
        result.push(if chunk.len() >= 2 {
            char::from(alphabet[usize::from((second & 15) << 2 | third >> 6)])
        } else {
            '='
        });
        result.push(if chunk.len() == 3 {
            char::from(alphabet[usize::from(third & 63)])
        } else {
            '='
        });
    }
    result
}

pub(super) fn card(id: &str, summary: &Summary) -> Read<Adapter> {
    let mut out = Adapter::default();
    match id {
        "ancestry-fine" => {
            let d = summary.section("ancestry_fine").data;
            let mut rows = Vec::new();
            for (level, label) in ancestry_fine::LEVELS {
                rows.extend(ancestry_fine::rows(d, level)?.into_iter().map(|r| {
                    vec![
                        label.into(),
                        cell(val(&r, "group")),
                        share(val(&r, "fraction")),
                        format!("{} to {}", share(val(&r, "low")), share(val(&r, "high"))),
                    ]
                }));
            }
            out.table(&["Level", "Group", "Share", "Bootstrap range"], rows);
            out.fact("Chromosome painting", ancestry_fine::painted(d)?);
        }
        "mt-haplogroup" => {
            let sec = summary.section("haplogroup");
            out.fact("Haplogroup", cell(val(sec.data, "haplogroup")));
            out.fact(
                "Two mtDNA haplogroups in the reads",
                haplogroups::contamination_text(&sec)?,
            );
            out.extra("explore", explore(sec.data)?);
        }
        "y-chromosome" => {
            let y = summary.section("y_haplogroup").data;
            let t = summary.section("y_tree").data;
            if truthy(y) {
                out.facts(
                    y,
                    &[
                        ("Haplogroup (Yleaf)", "haplogroup"),
                        ("QC score", "qc_score"),
                    ],
                );
                out.extra("explore", explore(y)?);
            }
            if truthy(t) {
                out.fact(
                    format!("Placement on {}", text(val(t, "tree"))),
                    cell(val(t, "haplogroup")),
                );
                out.fields_table(
                    &["Level", "Branch", "Derived", "Ancestral", "No call"],
                    array(t, "path")?,
                    &["Level", "Branch", "Derived", "Ancestral", "No_call"],
                );
                out.extra("widget", widgets::y_path(t));
            }
            let markers = array(summary.section("y_str").data, "markers")?;
            if !markers.is_empty() {
                out.fact(
                    "Y-STR markers read",
                    markers
                        .iter()
                        .filter(|m| truthy(val(m, "count")))
                        .count()
                        .to_string(),
                );
            }
        }
        "archaic" => {
            let d = summary.section("archaic").data;
            out.fact(
                "Neanderthal",
                format!("{}%", text(val(d, "neanderthal_percent"))),
            );
            out.fact(
                "Denisovan",
                format!("{}%", text(val(d, "denisovan_percent"))),
            );
            out.lines = archaic::lines(d)?;
        }
        "ancient" => {
            let d = summary.section("ancient").data;
            let nearest = array(d, "nearest")?;
            out.fields_table(
                &["Rank", "Group", "Dates", "Country", "Distance"],
                &nearest[..nearest.len().min(10)],
                &["rank", "group", "dates", "country", "distance"],
            );
            out.fact(
                "Sites used",
                format!(
                    "{} of {}",
                    text(val(d, "sites_used")),
                    text(val(d, "panel_sites"))
                ),
            );
            let svg = string(d, "plot_svg");
            if !svg.trim().is_empty() {
                let body = if svg.trim_start().starts_with('<') {
                    base64(svg.as_bytes())
                } else {
                    svg.to_owned()
                };
                out.extra(
                    "image",
                    obj([
                        (
                            "alt",
                            s("PCA of ancient genomes with this genome projected"),
                        ),
                        ("src", s(&format!("data:image/svg+xml;base64,{body}"))),
                    ]),
                );
            }
        }
        "roh" => {
            let sec = summary.section("roh");
            let d = sec.data;
            out.fact("Segments", cell(val(d, "segments")));
            out.fact("Total", format!("{} Mb", text(val(d, "total_mb"))));
            out.fact("Largest", format!("{} Mb", text(val(d, "largest_mb"))));
            out.facts(d, &[("F_ROH", "f_roh"), ("Nearest band", "relationship")]);
            out.lines = roh::lines(&sec)?;
            out.attention = truthy(val(d, "upd_flag"));
        }
        _ => unreachable!(),
    }
    Ok(out)
}
