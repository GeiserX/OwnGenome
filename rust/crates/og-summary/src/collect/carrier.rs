//! The sample's own carrier list: step 06's ClinVar hits in the genes of
//! `assets/carrier_genes.tsv`, one person, no partner needed. An allele the
//! list marks hypomorph or pseudodeficiency is set aside, never a carrier
//! finding on its own. The genes of autosomal recessive conditions that a
//! called deletion and a hit both reach are listed as two hits, and a hit
//! called homozygous inside a deletion is shown as hemizygous.

use super::assets::carrier_genes;
use super::clinvar::{Floor, Hit, HitRecords, allele_copies, carrier_genotype, obj};
use super::text::{self, gt_field};
use super::two_hits::{DeletionInputs, Spans, Variant, gene_deletions, two_hits};
use indexmap::IndexMap;
use pycompat::json::Value;

/// The section's data: the list, or why it could not be made.
#[derive(Debug, Clone, PartialEq)]
pub struct Carrier(pub Value);

impl Carrier {
    /// `carrier_list_or_note`. An `Err` is the collector's note (a hit
    /// whose POS is no integer).
    pub fn read(
        assets: &str,
        hits: &Result<HitRecords, String>,
        deletions: &DeletionInputs,
    ) -> Result<Carrier, String> {
        let unavailable = |e: String| Carrier(obj([("unavailable", Value::Str(e))]));
        let genes = match carrier_genes(assets) {
            Ok(g) => g,
            Err(e) => return Ok(unavailable(e)),
        };
        let recs = match hits {
            Ok((_, r)) => r,
            Err(e) => return Ok(unavailable(e.clone())),
        };
        let mut items: Vec<IndexMap<String, Value>> = Vec::new();
        // Which list each item is on: carriers (false) or set aside (true).
        let mut aside: Vec<bool> = Vec::new();
        let mut variants: Vec<(Variant, usize)> = Vec::new();
        for r in recs.iter().filter(|r| r.len() >= 8) {
            let row = Hit::parse(r, Floor::default())?;
            let gt = text::gt_of(r).to_owned();
            let (copies, ploidy) = allele_copies(&gt);
            let key = row.key();
            for g in row.gene.split(',') {
                let Some(x) = genes.get(g) else { continue };
                let mut item = IndexMap::new();
                item.insert("gene".into(), Value::Str(g.to_owned()));
                item.insert("variant".into(), Value::Str(row.variant()));
                item.insert("genotype".into(), Value::Str(carrier_genotype(&gt, &row)));
                item.insert("significance".into(), Value::Str(row.significance.clone()));
                item.insert("stars".into(), Value::Int(row.stars.into()));
                item.insert("inheritance".into(), Value::Str(x.inheritance.clone()));
                item.insert("condition".into(), Value::Str(x.condition.clone()));
                let sev = x.severity.get(&key);
                if let Some(s) = sev {
                    item.insert(
                        "kind".into(),
                        s.kind.clone().map_or(Value::Null, Value::Str),
                    );
                    item.insert("name".into(), Value::Str(s.name.clone()));
                }
                let genotype = match &item["genotype"] {
                    Value::Str(s) => s.clone(),
                    _ => String::new(),
                };
                items.push(item);
                aside.push(sev.is_some());
                if sev.is_none_or(|s| s.kind.as_deref() == Some("hypomorph")) {
                    variants.push((
                        Variant {
                            gene: g.to_owned(),
                            key: key.clone(),
                            chrom: row.chrom.clone(),
                            pos: row.pos,
                            text: row.variant(),
                            genotype,
                            copies,
                            ploidy,
                            gt: gt.clone(),
                            ps: gt_field(r, "PS"),
                            allele: sev.map(|s| s.name.clone()),
                        },
                        items.len() - 1,
                    ));
                }
            }
        }
        let spans: Spans = genes
            .iter()
            .map(|(g, x)| (g.clone(), x.span.clone()))
            .collect();
        let (dels, read, unread) = gene_deletions(deletions, &spans);
        let vs: Vec<Variant> = variants.iter().map(|(v, _)| v.clone()).collect();
        let found = two_hits(&vs, &dels);
        for (v, at) in &variants {
            if let Some(t) = found.hemi.get(&(v.gene.clone(), v.key.clone())) {
                let item = &mut items[*at];
                item.insert("genotype".into(), Value::Str("hemizygous".into()));
                item.insert(
                    "hemizygous".into(),
                    Value::Str(format!("called homozygous inside {t}")),
                );
            }
        }
        let mut rows: Vec<(&String, _)> = found.entries.iter().collect();
        rows.sort_by(|a, b| a.0.cmp(b.0));
        let two: Vec<Value> = rows
            .into_iter()
            .filter(|(g, _)| genes[*g].inheritance == "AR")
            .map(|(g, e)| {
                let mut m = e.to_json();
                m.insert("condition".into(), Value::Str(genes[g].condition.clone()));
                m.insert("inheritance".into(), Value::Str("AR".into()));
                Value::Object(m)
            })
            .collect();
        let (mut carriers, mut set_aside) = (Vec::new(), Vec::new());
        for (item, a) in items.into_iter().zip(aside) {
            if a {
                set_aside.push(Value::Object(item));
            } else {
                carriers.push(Value::Object(item));
            }
        }
        let strs = |v: Vec<String>| Value::Array(v.into_iter().map(Value::Str).collect());
        Ok(Carrier(obj([
            ("genes_on_list", Value::Int(genes.len() as i128)),
            ("list", Value::Str("assets/carrier_genes.tsv".into())),
            ("carriers", Value::Array(carriers)),
            ("set_aside", Value::Array(set_aside)),
            (
                "two_hits",
                obj([
                    ("genes", Value::Array(two)),
                    ("deletions_read", strs(read)),
                    ("deletions_not_read", strs(unread)),
                ]),
            ),
        ])))
    }

    pub fn to_json(&self) -> Value {
        self.0.clone()
    }
}
