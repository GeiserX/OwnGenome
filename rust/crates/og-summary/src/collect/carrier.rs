//! The sample's own carrier list: step 06's ClinVar hits in the genes of
//! `assets/carrier_genes.tsv`, one person, no partner needed. An allele the
//! list marks hypomorph or pseudodeficiency is set aside, never a carrier
//! finding on its own. The genes of autosomal recessive conditions that a
//! called deletion and a hit both reach are listed as two hits, and a hit
//! called homozygous inside a deletion is shown as hemizygous. A hit that
//! is part of step 06's combined allele, or the combined allele itself,
//! carries the note the ClinVar table prints for it.
//!
//! Each row keeps its zygosity and step 06's call quality: GQ, DP and the
//! low-confidence mark against the hits file's own call floor. A carrier
//! has one changed copy, so a row called homozygous or hemizygous is listed
//! apart (`hom_or_hemizygous`), never as a carrier. Two different
//! heterozygous hits in one recessive gene are counted as copies, the way
//! the secondary-findings card counts them: read-backed phase on different
//! copies is two hits, phase unknown is possible two hits, both on one copy
//! is one hit. A pair that needs a low-confidence call is possible two hits
//! (low-confidence call); a pair the confident calls prove alone keeps its
//! word, whatever weak call sits beside it.

use super::assets::carrier_genes;
use super::clinvar::{
    Floor, Hit, HitRecords, VarKey, allele_copies, carrier_genotype, obj, opt_int, zygosity_word,
};
use super::text::{self, gt_field};
use super::two_hits::{DeletionInputs, Spans, Variant, gene_deletions, phase_haps, two_hits_with};
use indexmap::IndexMap;
use pycompat::json::Value;
use std::collections::{BTreeSet, HashMap, HashSet};

/// The section's data: the list, or why it could not be made.
#[derive(Debug, Clone, PartialEq)]
pub struct Carrier(pub Value);

/// One heterozygous carrier hit of a recessive gene, as the pair rule
/// counts it.
#[derive(Clone)]
struct Het {
    key: VarKey,
    gt: String,
    ps: Option<String>,
    /// The parts of a combined allele (step 06's `PGP_JOINED_FROM`).
    joined: Vec<VarKey>,
    low: bool,
}

impl Het {
    /// The copies read-backed phase puts the alt allele on, with its phase
    /// set: only a `|` GT with every copy called and a non-empty PS.
    fn phased(&self) -> Option<(&str, BTreeSet<usize>)> {
        let ps = self.ps.as_deref().filter(|p| !p.is_empty())?;
        let h = phase_haps(Some(&self.gt)).filter(|h| !h.is_empty())?;
        Some((ps, h))
    }

    /// The copy a half call (`1/.`, `./1`) puts its alt allele on.
    fn half(&self) -> Option<usize> {
        let a: Vec<&str> = self.gt.split(['/', '|']).collect();
        match a.as_slice() {
            [x, "."] if !matches!(*x, "0" | ".") => Some(0),
            [".", y] if !matches!(*y, "0" | ".") => Some(1),
            _ => None,
        }
    }
}

/// `CHROM:POS:REF:ALT` joined by `|`, as step 06 writes `PGP_JOINED_FROM`.
fn joined_keys(r: &[String]) -> Vec<VarKey> {
    let info = r.get(7).map(String::as_str).unwrap_or("");
    let v = text::info_map(info)
        .get("PGP_JOINED_FROM")
        .copied()
        .unwrap_or("");
    v.split('|')
        .filter_map(|p| {
            let mut it = p.rsplitn(4, ':');
            let (alt, rf, pos, chrom) = (it.next()?, it.next()?, it.next()?, it.next()?);
            let pos = crate::int(pos).ok()?;
            Some((chrom.to_owned(), pos, rf.to_owned(), alt.to_owned()))
        })
        .collect()
}

/// How two or more different heterozygous hits of one gene pair up: the
/// pair word (`two hits`, `possible two hits`, none when they sit on one
/// copy) and why. `None` when fewer than two distinct alleles are left
/// once a combined allele and its parts count once. The hets are counted
/// twice, all of them and the confident ones alone: two hits when the
/// confident ones prove it; else possible two hits when the confident ones
/// still count two copies of unknown phase, whatever a low-confidence call
/// adds; else the low-confidence word when a pair needs a low-confidence
/// call.
fn het_pair(hets: &[Het]) -> Option<(Option<String>, String)> {
    // A combined allele and its parts count once, under the combined
    // allele, when its call is phased or at most one part is here; an
    // unphased one with two or more parts here counts through its parts,
    // which may sit on different copies.
    let present: HashSet<&VarKey> = hets.iter().map(|h| &h.key).collect();
    let mut gone: HashSet<&VarKey> = HashSet::new();
    for h in hets.iter().filter(|h| !h.joined.is_empty()) {
        let parts: Vec<&VarKey> = h.joined.iter().filter(|k| present.contains(k)).collect();
        if h.phased().is_some() || parts.len() <= 1 {
            gone.extend(parts);
        } else {
            gone.insert(&h.key);
        }
    }
    let kept: Vec<&Het> = hets.iter().filter(|h| !gone.contains(&h.key)).collect();
    if kept.len() < 2 {
        return None;
    }
    // A kept allele is firm when a call of it is confident: its own, or a
    // call folded into it (a part under its combined allele, or an unphased
    // combined allele under its first part), as the secondary-findings card
    // counts it. The confident count takes that call's GT and phase set: a
    // phase only a low-confidence call gives proves nothing.
    let mut firm: HashMap<&VarKey, &Het> = kept
        .iter()
        .filter(|h| !h.low)
        .map(|h| (&h.key, *h))
        .collect();
    for h in hets.iter().filter(|h| !h.joined.is_empty()) {
        if gone.contains(&h.key) {
            if !h.low
                && let Some(k) = h.joined.iter().find(|k| present.contains(k))
            {
                firm.entry(k).or_insert(h);
            }
        } else if let Some(p) = hets.iter().find(|p| !p.low && h.joined.contains(&p.key)) {
            firm.entry(&h.key).or_insert(p);
        }
    }
    // The copies a set of hets hits: the count, the proof of two copies if a
    // call gives one, and the first phase set holding two or more of them on
    // one copy (cis).
    let count = |hs: &[&Het]| -> (usize, Option<String>, Option<String>) {
        let mut n = 0;
        let mut sets: IndexMap<&str, (BTreeSet<usize>, usize)> = IndexMap::new();
        let mut halves: IndexMap<(&str, i128), (BTreeSet<usize>, usize)> = IndexMap::new();
        for h in hs {
            if let Some((ps, s)) = h.phased() {
                let e = sets.entry(ps).or_default();
                e.0.extend(s);
                e.1 += 1;
            } else if let Some(i) = h.half() {
                let e = halves.entry((h.key.0.as_str(), h.key.1)).or_default();
                e.0.insert(i);
                e.1 += 1;
            } else {
                n += 1;
            }
        }
        let mut proof: Option<String> = None;
        let mut cis: Option<String> = None;
        for (ps, (h, count)) in &sets {
            n += h.len();
            if h.len() >= 2 {
                proof.get_or_insert_with(|| {
                    format!(
                        "in trans: read-backed phase puts them on different copies (phase set {ps})"
                    )
                });
            } else if *count >= 2 {
                cis.get_or_insert_with(|| (*ps).to_owned());
            }
        }
        for ((chrom, pos), (h, count)) in &halves {
            if h.len() >= 2 {
                n += 2 + count - h.len();
                proof.get_or_insert_with(|| {
                    format!("in trans: the two alt alleles of one call at {chrom}:{pos}")
                });
            } else {
                n += count;
            }
        }
        (n, proof, cis)
    };
    let (n, proof, cis) = count(&kept);
    let sure: Vec<Het> = kept
        .iter()
        .filter_map(|h| {
            firm.get(&h.key).map(|c| Het {
                gt: c.gt.clone(),
                ps: c.ps.clone(),
                ..(*h).clone()
            })
        })
        .collect();
    let (n_firm, proof_firm, _) = count(&sure.iter().collect::<Vec<_>>());
    // A pair that needs a low-confidence call is only possible until the
    // call is confirmed.
    let low = |why: String| {
        Some((
            Some("possible two hits (low-confidence call)".to_owned()),
            format!(
                "a call below the call floor (low confidence) is one of them, so confirm it first; {why}"
            ),
        ))
    };
    let unknown = "phase unknown: no read-backed phase links the variants, so they may sit on the \
                   same copy (cis)"
        .to_owned();
    match (proof_firm, proof, cis) {
        (Some(f), _, _) => Some((Some("two hits".to_owned()), f)),
        (None, Some(_), _) if n_firm >= 2 => Some((Some("possible two hits".to_owned()), unknown)),
        (None, Some(p), _) => low(p),
        (None, None, _) if n >= 2 && n_firm >= 2 => {
            Some((Some("possible two hits".to_owned()), unknown))
        }
        (None, None, _) if n >= 2 => low(unknown),
        (None, None, Some(ps)) => Some((
            None,
            format!(
                "two heterozygous variants in cis: read-backed phase puts them on the same copy \
                 (phase set {ps}), so one copy is hit"
            ),
        )),
        (None, None, None) => Some((
            None,
            "the variants sit on one copy, so one copy is hit".to_owned(),
        )),
    }
}

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
        // The hits file's own call floor (`##PGP_CALL_FLOOR`), as the
        // ClinVar section judges the same rows.
        let (floor, recs): (Floor, _) = match hits {
            Ok((f, r)) => (*f, r),
            Err(e) => return Ok(unavailable(e.clone())),
        };
        let mut items: Vec<IndexMap<String, Value>> = Vec::new();
        // Which list each item is on: carriers (false) or set aside (true).
        let mut aside: Vec<bool> = Vec::new();
        let mut variants: Vec<(Variant, usize)> = Vec::new();
        // The parts of each variant's combined allele, by the same index.
        let mut vjoined: Vec<Vec<VarKey>> = Vec::new();
        let mut low: HashSet<(String, VarKey)> = HashSet::new();
        // The heterozygous carrier hits of each recessive gene, with the
        // items they are.
        let mut hets: IndexMap<String, Vec<(Het, usize)>> = IndexMap::new();
        for r in recs.iter().filter(|r| r.len() >= 8) {
            let row = Hit::parse(r, floor)?;
            let gt = text::gt_of(r).to_owned();
            let (copies, ploidy) = allele_copies(&gt);
            let key = row.key();
            for g in row.gene.split(',') {
                let Some(x) = genes.get(g) else { continue };
                let mut item = IndexMap::new();
                item.insert("gene".into(), Value::Str(g.to_owned()));
                item.insert("variant".into(), Value::Str(row.variant()));
                item.insert("genotype".into(), Value::Str(carrier_genotype(&gt, &row)));
                item.insert(
                    "zygosity".into(),
                    Value::Str(zygosity_word(copies, ploidy).into()),
                );
                item.insert("significance".into(), Value::Str(row.significance.clone()));
                item.insert("stars".into(), Value::Int(row.stars.into()));
                item.insert("gq".into(), opt_int(row.gq));
                item.insert("dp".into(), opt_int(row.dp));
                item.insert("low_confidence".into(), Value::Bool(row.low_confidence));
                item.insert("inheritance".into(), Value::Str(x.inheritance.clone()));
                item.insert("condition".into(), Value::Str(x.condition.clone()));
                if row.combined_allele != "." {
                    item.insert(
                        "combined_allele".into(),
                        Value::Str(row.combined_allele.clone()),
                    );
                }
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
                if row.low_confidence {
                    low.insert((g.to_owned(), key.clone()));
                }
                if sev.is_none() && x.inheritance == "AR" && copies == 1 && ploidy == 2 {
                    hets.entry(g.to_owned()).or_default().push((
                        Het {
                            key: key.clone(),
                            gt: gt.clone(),
                            ps: gt_field(r, "PS"),
                            joined: joined_keys(r),
                            low: row.low_confidence,
                        },
                        items.len() - 1,
                    ));
                }
                // A record with no alt copy called is no hit to pair.
                if copies > 0 && sev.is_none_or(|s| s.kind.as_deref() == Some("hypomorph")) {
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
                    vjoined.push(joined_keys(r));
                }
            }
        }
        for list in hets.values() {
            // The same record twice (a duplicate line) is one allele.
            let mut seen = HashSet::new();
            let distinct: Vec<Het> = list
                .iter()
                .filter(|(h, _)| seen.insert(h.key.clone()))
                .map(|(h, _)| h.clone())
                .collect();
            let Some((word, why)) = het_pair(&distinct) else {
                continue;
            };
            for (_, at) in list {
                let item = &mut items[*at];
                if let Some(w) = &word {
                    item.insert("pair".into(), Value::Str(w.clone()));
                }
                item.insert("phase".into(), Value::Str(why.clone()));
            }
        }
        let spans: Spans = genes
            .iter()
            .map(|(g, x)| (g.clone(), x.span.clone()))
            .collect();
        let (dels, read, unread) = gene_deletions(deletions, &spans);
        // A low-confidence row pairs with a deletion by a confident call of
        // its combined allele, as the secondary-findings card counts it: the
        // first confident row of the gene, with as many copies, that is a
        // part of it, or an unphased combined allele whose first part present
        // it is (the card counts such an allele under that part). A phased
        // combined allele stands for its parts, so it never lends them its
        // phase. That row's GT and phase set stand for it, and the pair no
        // longer needs a low-confidence call.
        let mut vs: Vec<Variant> = variants.iter().map(|(v, _)| v.clone()).collect();
        let mut low_pair = low.clone();
        let is_low = |v: &Variant| low.contains(&(v.gene.clone(), v.key.clone()));
        let phased = |v: &Variant| {
            v.ps.as_deref().is_some_and(|p| !p.is_empty())
                && phase_haps(Some(&v.gt)).is_some_and(|h| !h.is_empty())
        };
        // The first part of combined row `j` present among the gene's
        // variants, in `PGP_JOINED_FROM` order.
        let first_part = |j: usize| {
            let g = &variants[j].0.gene;
            vjoined[j]
                .iter()
                .find(|k| variants.iter().any(|(x, _)| &x.gene == g && &x.key == *k))
        };
        for (i, (v, _)) in variants.iter().enumerate() {
            if !is_low(v) {
                continue;
            }
            let sure = variants.iter().enumerate().find(|(j, (c, _))| {
                c.gene == v.gene
                    && !is_low(c)
                    && c.copies >= v.copies
                    && (vjoined[i].contains(&c.key)
                        || (!phased(c) && first_part(*j) == Some(&v.key)))
            });
            if let Some((_, (c, _))) = sure {
                vs[i].gt = c.gt.clone();
                vs[i].ps = c.ps.clone();
                low_pair.remove(&(v.gene.clone(), v.key.clone()));
            }
        }
        let found = two_hits_with(&vs, &dels, |v| {
            low_pair.contains(&(v.gene.clone(), v.key.clone()))
        });
        for (v, at) in &variants {
            if let Some(t) = found.hemi.get(&(v.gene.clone(), v.key.clone())) {
                let item = &mut items[*at];
                item.insert("genotype".into(), Value::Str("hemizygous".into()));
                item.insert("zygosity".into(), Value::Str("hemizygous".into()));
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
        let (mut carriers, mut both, mut set_aside) = (Vec::new(), Vec::new(), Vec::new());
        for (item, a) in items.into_iter().zip(aside) {
            let zyg = match item.get("zygosity") {
                Some(Value::Str(z)) => z.clone(),
                _ => String::new(),
            };
            if a {
                set_aside.push(Value::Object(item));
            } else if matches!(zyg.as_str(), "hom" | "hemizygous") {
                // No unchanged copy left: not the pattern of a carrier.
                both.push(Value::Object(item));
            } else {
                carriers.push(Value::Object(item));
            }
        }
        let strs = |v: Vec<String>| Value::Array(v.into_iter().map(Value::Str).collect());
        Ok(Carrier(obj([
            ("genes_on_list", Value::Int(genes.len() as i128)),
            ("list", Value::Str("assets/carrier_genes.tsv".into())),
            ("carriers", Value::Array(carriers)),
            ("hom_or_hemizygous", Value::Array(both)),
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
