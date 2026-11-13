//! Step 09 (ExpansionHunter) and 9b (Stranger, the RFC1 and other motif
//! counts): every locus with its counts, how well each allele was sized,
//! Stranger's grade, the intermediate zone, and what the counts of CFTR's
//! TG and poly-T tracts and UGT1A1's TA repeat mean.

use super::super::model::State;
use super::util::{Obj, bin_rows, cells, count, get_default, s, strs, table, vcf_lines};
use super::{Ctx, Out};
use indexmap::IndexMap;
use pycompat::json::Value;

/// The five loci the report names first.
pub const EH_LOCI: [&str; 5] = ["HTT", "FMR1", "C9ORF72", "ATXN1", "DMPK"];
/// Fewer spanning reads than this: a low-confidence size.
const MIN_SPANNING: i128 = 2;
const MOTIF_KEYS: [&str; 7] = [
    "verdict",
    "note",
    "dominant_motif",
    "dominant_class",
    "dominant_share",
    "motif_shares",
    "repeat_reads",
];
/// The repeats STRchive leaves out, read when step 09 sized them.
const EXTRA_SETUP: &str = "STR_EXTRA_LOCI=yes ./scripts/setup.sh --strchive-data <genome_dir>";
const CFTR_R117H: &str = "rs78655421";
const CFTR_SOURCE: &str = "Deignan et al. 2020, PMID 32404922";
const UGT1A1_SOURCE: &str = "Gammal et al. 2016, PMID 26417955";
const HTTLPR_READING: &str = "not tested by any step. CPIC's 2023 SSRI guideline (Bousman et al. 2023, PMID 37032427) found that the data on SLC6A4 do not support its use in antidepressant prescribing, so a result would carry no recommendation.";

/// How well ExpansionHunter sized an allele, worst first; a locus takes
/// the worst of its alleles.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Sizing {
    NoCall,
    FlankingOnly,
    InRepeat,
    LowConfidence,
    Sized,
    Unknown,
}

impl Sizing {
    pub fn word(self) -> &'static str {
        match self {
            Sizing::NoCall => "no_call",
            Sizing::FlankingOnly => "flanking_only",
            Sizing::InRepeat => "in_repeat",
            Sizing::LowConfidence => "low_confidence",
            Sizing::Sized => "sized",
            Sizing::Unknown => "unknown",
        }
    }

    fn of(w: &str) -> Option<Sizing> {
        [
            Sizing::NoCall,
            Sizing::FlankingOnly,
            Sizing::InRepeat,
            Sizing::LowConfidence,
            Sizing::Sized,
            Sizing::Unknown,
        ]
        .into_iter()
        .find(|z| z.word() == w)
    }

    /// Sized well enough to read a count.
    fn read(self) -> bool {
        matches!(self, Sizing::Sized | Sizing::Unknown)
    }
}

fn str_int(x: &str) -> Option<i128> {
    crate::int(x).ok()
}

/// One allele's sizing from ExpansionHunter's SO, ADSP and REPCI fields.
fn allele_sizing(so: &str, adsp: &str, ci: &str, cn: &str, path_min: Option<i128>) -> Sizing {
    let Some(n) = str_int(cn) else {
        return Sizing::NoCall;
    };
    match so {
        "" | "." => return Sizing::Unknown,
        "FLANKING" => return Sizing::FlankingOnly,
        "INREPEAT" => return Sizing::InRepeat,
        _ => {}
    }
    let hi = if ci.is_empty() {
        None
    } else {
        str_int(ci.rsplit('-').next().unwrap_or(ci))
    };
    let few = str_int(adsp).is_some_and(|a| a < MIN_SPANNING);
    let reaches = matches!((path_min, hi), (Some(m), Some(h)) if n < m && m <= h);
    if few || reaches {
        Sizing::LowConfidence
    } else {
        Sizing::Sized
    }
}

/// A per-allele FORMAT value as `n` parts; one value serves every allele.
fn split_n(x: &str, n: usize) -> Vec<String> {
    let parts: Vec<&str> = if x.is_empty() {
        Vec::new()
    } else {
        x.split('/').collect()
    };
    if parts.len() == n {
        parts.into_iter().map(str::to_owned).collect()
    } else if parts.len() == 1 {
        vec![parts[0].to_owned(); n]
    } else {
        vec![String::new(); n]
    }
}

/// `bin/str_intermediate_ranges.tsv` by locus.
fn ranges(c: &Ctx<'_>) -> Result<IndexMap<String, IndexMap<String, String>>, String> {
    let mut out = IndexMap::new();
    for r in bin_rows(&c.table("str_intermediate_ranges.tsv"))? {
        let k = r
            .get("locus")
            .cloned()
            .ok_or_else(|| super::util::key_error("locus"))?;
        out.insert(k, r);
    }
    Ok(out)
}

fn mtime(p: &str) -> Result<f64, String> {
    pycompat::mtime::mtime(p).map_err(|e| super::util::os_note(&crate::os_error(&e, p)))
}

/// The motif table `name` unless it is older than Stranger's file.
fn motif_table(c: &Ctx<'_>, name: &str, after: &str) -> Result<Option<super::util::Table>, String> {
    let Some(p) = c.find1(format!("expansion_hunter/{}_{name}", c.s())) else {
        return Ok(None);
    };
    if mtime(&p.path)? < mtime(after)? {
        return Ok(None);
    }
    table(&p.path, false).map(Some)
}

fn motif_row(r: &pycompat::csv::DictRow) -> Value {
    MOTIF_KEYS
        .iter()
        .fold(Obj::new(), |o, k| o.set(k, get_default(r, k, ".")))
        .val()
}

/// `{locus: reason}` of the loci the catalog step 09 ran on gives no
/// thresholds (its NoThresholds field), from the copy beside its VCF.
fn not_graded(c: &Ctx<'_>) -> IndexMap<String, String> {
    let mut out = IndexMap::new();
    let Some(p) = c.find1(format!("expansion_hunter/{}_eh_catalog.json", c.s())) else {
        return out;
    };
    let Ok(text) = std::fs::read(&p.path).map(String::from_utf8) else {
        return out;
    };
    let Ok(Ok(Value::Array(cat))) = text.map(|t| pycompat::json::loads(&t)) else {
        return out;
    };
    for e in cat {
        let Value::Object(m) = e else { continue };
        let (Some(id), Some(why)) = (m.get("LocusId"), m.get("NoThresholds")) else {
            continue;
        };
        if pycompat::value::truthy(id)
            && pycompat::value::truthy(why)
            && let (Value::Str(id), Ok(why)) = (id, pycompat::value::py_str(why))
        {
            out.insert(id.clone(), why);
        }
    }
    out
}

/// One locus as the reports read it.
struct Locus {
    json: Obj,
    locus: String,
}

/// Step 09's ExpansionHunter VCF, or Stranger's copy of it with a grade per
/// locus when it is not older than the VCF.
pub fn expansions(c: &Ctx<'_>) -> Out {
    let s_ = c.s();
    let eh = c.find(&[
        format!("expansion_hunter/{s_}_eh.vcf"),
        format!("expansion_hunter/{s_}.vcf"),
        "expansion_hunter/*_eh.vcf".to_owned(),
    ]);
    let mut st = c.find(&[
        format!("expansion_hunter/{s_}_eh_stranger.vcf"),
        "expansion_hunter/*_eh_stranger.vcf".to_owned(),
    ]);
    let outdated = match (&st, &eh) {
        (Some(a), Some(b)) => mtime(&a.path)? < mtime(&b.path)?,
        _ => false,
    };
    if outdated {
        st = None;
    }
    let Some(p) = st.clone().or(eh) else {
        return Ok(None);
    };
    let ranges = ranges(c)?;
    let stp = st.as_ref().map(|f| f.path.as_str());
    let rfc1 = match stp {
        Some(after) => {
            motif_table(c, "rfc1_motif.tsv", after)?.and_then(|t| t.rows.first().map(motif_row))
        }
        None => None,
    };
    let mut motifs: IndexMap<String, Value> = IndexMap::new();
    if let Some(after) = stp
        && let Some(t) = motif_table(c, "str_motifs.tsv", after)?
    {
        for r in &t.rows {
            if let Some(l) = r.get("locus").flatten().filter(|l| !l.is_empty()) {
                motifs.insert(l.to_owned(), motif_row(r));
            }
        }
    }
    let ungraded = if st.is_some() {
        not_graded(c)
    } else {
        IndexMap::new()
    };
    let mut records = Vec::new();
    vcf_lines(&p.path, |l| records.push(cells(l)))?;
    let (mut loci, mut flagged, mut no_status, mut graded_none) =
        (Vec::<Locus>::new(), Vec::new(), 0, Vec::<String>::new());
    let mut tested = 0;
    for r in &records {
        if r.len() < 10 {
            continue;
        }
        tested += 1;
        let info: IndexMap<&str, &str> = r[7]
            .split(';')
            .map(|kv| kv.split_once('=').unwrap_or((kv, "")))
            .collect();
        let get = |k: &str| info.get(k).copied().filter(|v| !v.is_empty());
        let lid = get("REPID")
            .or(get("VARID"))
            .map_or_else(|| format!("{}:{}", r[0], r[1]), str::to_owned);
        let rep = get("VARID").map_or_else(|| lid.clone(), str::to_owned);
        let fmt: Vec<&str> = r[8].split(':').collect();
        let val: Vec<&str> = r[9].split(':').collect();
        let field = |k: &str| {
            fmt.iter()
                .position(|f| *f == k)
                .and_then(|i| val.get(i))
                .copied()
                .unwrap_or("")
        };
        let repcn = field("REPCN");
        let counts: Vec<&str> = if repcn.is_empty() || repcn == "." {
            Vec::new()
        } else {
            repcn.split('/').collect()
        };
        let n = counts.len().max(1);
        let path_min = info.get("STR_PATHOLOGIC_MIN").and_then(|v| str_int(v));
        let cn_list: Vec<&str> = if counts.is_empty() {
            vec!["."]
        } else {
            counts.clone()
        };
        let sizing: Vec<Sizing> = split_n(field("SO"), n)
            .iter()
            .zip(split_n(field("ADSP"), n))
            .zip(split_n(field("REPCI"), n))
            .zip(cn_list)
            .map(|(((so, ad), ci), cn)| allele_sizing(so, &ad, &ci, cn, path_min))
            .collect();
        let status: Vec<&str> = if st.is_some() {
            info.get("STR_STATUS")
                .copied()
                .unwrap_or("")
                .split(',')
                .filter(|x| !x.is_empty() && *x != ".")
                .collect()
        } else {
            Vec::new()
        };
        let worst = sizing.iter().copied().min().unwrap_or(Sizing::NoCall);
        let dots = |k: &str| {
            info.get(k)
                .copied()
                .unwrap_or("")
                .trim_matches('.')
                .to_owned()
        };
        let mut item = Obj::new()
            .set("locus", s(&rep))
            .set(
                "repeat_count",
                s(if repcn.is_empty() { "." } else { repcn }),
            )
            .set("status", s(&status.join(",")))
            .set("normal_max", s(&dots("STR_NORMAL_MAX")))
            .set("pathologic_min", s(&dots("STR_PATHOLOGIC_MIN")))
            .set("support", s(field("SO")))
            .set("adsp", s(field("ADSP")))
            .set("ci", s(field("REPCI")))
            .set("sizing", s(worst.word()))
            .set(
                "allele_sizing",
                s(&sizing
                    .iter()
                    .map(|z| z.word())
                    .collect::<Vec<_>>()
                    .join("/")),
            )
            .set("zone", s(""))
            .set("zone_label", s(""))
            .set("zone_citation", s(""));
        let sized: Vec<i128> = counts
            .iter()
            .zip(&sizing)
            .filter(|(_, z)| z.read())
            .filter_map(|(c, _)| str_int(c))
            .collect();
        if let Some(rng) = ranges.get(&rep) {
            let bound = |k: &str| -> Result<i128, String> {
                let v = rng.get(k).ok_or_else(|| super::util::key_error(k))?;
                crate::int(v).map_err(|e| super::super::sections::py_error("ValueError", &e))
            };
            if !sized.is_empty() {
                let (lo, hi) = (bound("intermediate_min")?, bound("intermediate_max")?);
                if sized.iter().any(|c| lo <= *c && *c <= hi) {
                    let cell =
                        |k: &str| rng.get(k).cloned().ok_or_else(|| super::util::key_error(k));
                    item.put(
                        "zone",
                        s(&format!(
                            "{}-{}",
                            cell("intermediate_min")?,
                            cell("intermediate_max")?
                        )),
                    );
                    item.put("zone_label", s(&cell("label")?));
                    item.put("zone_citation", s(&cell("citation")?));
                }
            }
        }
        let abnormal = status.iter().any(|x| *x != "normal");
        if rep == "RFC1" && st.is_some() && abnormal {
            item.put(
                "rfc1_motif",
                rfc1.clone().unwrap_or_else(|| {
                    Obj::new()
                        .set("verdict", s("not_checked"))
                        .set("note", s("no motif count beside Stranger's output (step 9b runs it with the sample's BAM): confirm the flag with repeat-primed PCR"))
                        .val()
                }),
            );
        }
        if rep != "RFC1"
            && abnormal
            && let Some(m) = motifs.get(&rep)
        {
            item.put("motif", m.clone());
        }
        if st.is_some()
            && status.is_empty()
            && let Some(why) = ungraded.get(&lid)
        {
            item.put("not_graded", s(why));
        }
        if st.is_some() {
            if status.is_empty() && item.0.contains_key("not_graded") {
                if !graded_none.contains(&lid) {
                    graded_none.push(lid.clone());
                }
            } else if status.is_empty() {
                no_status += 1;
            } else if abnormal {
                let mut f = Obj::new();
                for k in ["locus", "repeat_count", "status", "sizing"] {
                    f.put(k, item.0[k].clone());
                }
                for k in ["rfc1_motif", "motif"] {
                    if let Some(v) = item.0.get(k) {
                        f.put(k, v.clone());
                    }
                }
                flagged.push(f.val());
            }
        }
        loci.push(Locus {
            json: item,
            locus: rep,
        });
    }
    let by_unsorted = by_locus(loci.iter().map(|l| (&l.locus, &l.json)));
    loci.sort_by(|a, b| {
        let rank = |l: &str| {
            EH_LOCI
                .iter()
                .position(|k| *k == l)
                .unwrap_or(EH_LOCI.len())
        };
        rank(&a.locus)
            .cmp(&rank(&b.locus))
            .then_with(|| a.locus.cmp(&b.locus))
    });
    let key_loci = EH_LOCI
        .iter()
        .map(|k| {
            Obj::new()
                .set("locus", s(k))
                .set(
                    "repeat_count",
                    by_unsorted
                        .get(*k)
                        .map_or_else(|| s("not in output"), |l| l.repeat_count.clone()),
                )
                .val()
        })
        .collect();
    let word = |l: &Locus, k: &str| match l.json.0.get(k) {
        Some(Value::Str(x)) => x.clone(),
        _ => String::new(),
    };
    let unsized_: Vec<String> = loci
        .iter()
        .filter(|l| !Sizing::of(&word(l, "sizing")).is_some_and(Sizing::read))
        .map(|l| l.locus.clone())
        .collect();
    let intermediate: Vec<String> = loci
        .iter()
        .filter(|l| !word(l, "zone").is_empty())
        .map(|l| l.locus.clone())
        .collect();
    let readings = readings(c, &by_unsorted, State::Ok, State::Ok)?;
    let mut out = Obj::new()
        .set("records", count(tested))
        .set("key_loci", Value::Array(key_loci))
        .set(
            "loci",
            Value::Array(loci.into_iter().map(|l| l.json.val()).collect()),
        )
        .set("stranger", Value::Bool(st.is_some()))
        .set("stranger_outdated", Value::Bool(outdated))
        .set("unsized", strs(&unsized_))
        .set("intermediate", strs(&intermediate))
        .set("readings", readings);
    if st.is_some() {
        out.put("flagged", Value::Array(flagged));
        out.put("no_status", count(no_status));
        out.put(
            "not_graded",
            Value::Object(
                graded_none
                    .iter()
                    .map(|k| (k.clone(), s(&ungraded[k])))
                    .collect(),
            ),
        );
    }

    Ok(Some((p, out.val())))
}

/// What a reading needs of a locus: its count and how well it was sized.
#[derive(Debug, Clone)]
struct Sized {
    repeat_count: Value,
    sizing: Option<Sizing>,
}

/// `{locus: l}` over loci in this order, the last of a name kept.
fn by_locus<'a>(loci: impl Iterator<Item = (&'a String, &'a Obj)>) -> IndexMap<String, Sized> {
    loci.map(|(k, o)| {
        let get = |f: &str| o.0.get(f).cloned().unwrap_or(Value::Null);
        let sizing = match get("sizing") {
            Value::Str(w) => Sizing::of(&w),
            _ => None,
        };
        (
            k.clone(),
            Sized {
                repeat_count: get("repeat_count"),
                sizing,
            },
        )
    })
    .collect()
}

/// The two counts of a sized record, sorted.
fn counts(item: Option<&Sized>) -> Option<[i128; 2]> {
    let item = item?;
    if !item.sizing.is_some_and(Sizing::read) {
        return None;
    }
    let Value::Str(rc) = &item.repeat_count else {
        return None;
    };
    let mut c: Vec<Option<i128>> = rc.split('/').map(str_int).collect();
    if c.len() == 1 && c[0].is_some() {
        c.push(c[0]);
    }
    match c[..] {
        [Some(a), Some(b)] => Some(if a <= b { [a, b] } else { [b, a] }),
        _ => None,
    }
}

fn text_of(v: &Value) -> String {
    match v {
        Value::Str(x) => x.clone(),
        other => pycompat::value::py_str(other).unwrap_or_default(),
    }
}

fn sizing_text(item: &Sized) -> String {
    item.sizing.map_or("None", Sizing::word).to_owned()
}

/// What a section state other than ok says about another step's table.
fn other_state(st: State) -> &'static str {
    match st {
        State::Stale => "from an earlier run",
        State::Failed => "its step failed in the last run",
        _ => "unreadable",
    }
}

/// R117H's genotype in step 42's traits table: (genotype, why not read,
/// a note on a table of an earlier run).
fn r117h(c: &Ctx<'_>, state: State) -> (Option<String>, String, String) {
    if matches!(state, State::Failed | State::Unreadable) {
        return (
            None,
            format!("step 42's traits table is {}", other_state(state)),
            String::new(),
        );
    }
    let row = c
        .find1(format!("traits/{}_traits.tsv", c.s()))
        .and_then(|p| table(&p.path, false).ok())
        .and_then(|t| {
            t.rows
                .into_iter()
                .find(|r| r.get("rsid").flatten() == Some(CFTR_R117H))
        });
    let called = row.as_ref().map(|r| {
        (
            super::util::or(r, "status", ""),
            super::util::or(r, "genotype", ""),
        )
    });
    match called {
        Some((st, gt)) if st == "called" => (
            Some(gt),
            String::new(),
            if state == State::Stale {
                " (R117H from step 42's traits table of an earlier run)".into()
            } else {
                String::new()
            },
        ),
        _ => (
            None,
            "no step 42 traits row for rs78655421 was called".into(),
            String::new(),
        ),
    }
}

fn a_copies(gt: &str) -> usize {
    gt.split('/').filter(|x| *x == "A").count()
}

fn cftr(c: &Ctx<'_>, by: &IndexMap<String, Sized>, traits: State) -> Value {
    let out = Obj::new()
        .set("locus", s("CFTR"))
        .set("repeat", s("intron 9 (legacy IVS8) TG and poly-T tract"))
        .set("source", s(CFTR_SOURCE));
    let not_sized = |o: Obj, reading: String| {
        o.set("status", s("not_sized"))
            .set("genotype", s("."))
            .set("phase", s("."))
            .set("reading", s(&reading))
            .val()
    };
    let (Some(tg), Some(t)) = (by.get("CFTR_TG"), by.get("CFTR_T")) else {
        return not_sized(
            out,
            format!(
                "not in step 09's output: {EXTRA_SETUP} adds the tract to the catalog (docs/09-str-expansions.md)."
            ),
        );
    };
    let (Some(a), Some(b)) = (counts(Some(tg)), counts(Some(t))) else {
        return not_sized(
            out,
            format!(
                "ExpansionHunter did not size both tracts (TG {}, {}; poly-T {}, {}): not read.",
                text_of(&tg.repeat_count),
                sizing_text(tg),
                text_of(&t.repeat_count),
                sizing_text(t)
            ),
        );
    };
    let pairs: Vec<String> = if a[0] == a[1] || b[0] == b[1] {
        vec![format!("TG{}-{}T/TG{}-{}T", a[0], b[0], a[1], b[1])]
    } else {
        vec![
            format!("TG{}-{}T/TG{}-{}T", a[0], b[0], a[1], b[1]),
            format!("TG{}-{}T/TG{}-{}T", a[1], b[0], a[0], b[1]),
        ]
    };
    let genotype = if pairs.len() == 1 {
        pairs[0].clone()
    } else {
        format!("TG{}/TG{} {}T/{}T", a[0], a[1], b[0], b[1])
    };
    let ints = |x: [i128; 2]| Value::Array(x.iter().map(|v| Value::Int(*v)).collect());
    let out = out
        .set("status", s("read"))
        .set("genotype", s(&genotype))
        .set("pairings", strs(&pairs))
        .set("tg", ints(a))
        .set("poly_t", ints(b));
    let five = b.iter().filter(|n| **n == 5).count();
    let mut other: Vec<i128> = b
        .iter()
        .copied()
        .filter(|n| !matches!(n, 5 | 7 | 9))
        .collect();
    other.sort_unstable();
    other.dedup();
    let odd = if other.is_empty() {
        String::new()
    } else {
        format!(
            " A poly-T length of {} is none of the usual 5T, 7T and 9T.",
            other
                .iter()
                .map(|n| format!("{n}T"))
                .collect::<Vec<_>>()
                .join(", ")
        )
    };
    let (r117, why, note) = r117h(c, traits);
    let r117_text = s(r117.as_deref().filter(|x| !x.is_empty()).unwrap_or("."));
    if five == 0 {
        let mut line = format!(
            "No 5T copy ({}T/{}T): the tract adds no 5T allele of its own. The poly-T length still matters beside R117H, which is milder and of variable penetrance on a 7T or 9T copy and CF-causing with 5T in cis.",
            b[0], b[1]
        );
        match &r117 {
            None => line += &format!(" R117H was not read ({why})."),
            Some(g) if a_copies(g) == 2 => {
                let copies = if b[0] == b[1] {
                    format!("both {}T copies", b[0])
                } else {
                    format!("both copies ({}T and {}T)", b[0], b[1])
                };
                line += &format!(" R117H ({g}) is read: it sits on {copies}.{note}");
            }
            Some(g) if a_copies(g) > 0 => {
                let copies = if b[0] == b[1] {
                    format!("{}T", b[0])
                } else {
                    format!("{}T or {}T", b[0], b[1])
                };
                line += &format!(
                    " R117H ({g}) is read: with no 5T copy it sits on a {copies} copy.{note}"
                );
            }
            Some(_) => line += &format!(" No R117H copy.{note}"),
        }
        return out
            .set("five_t", Value::Bool(false))
            .set("phase", s("not_needed"))
            .set("r117h", r117_text)
            .set("reading", s(&(line + &odd)))
            .val();
    }
    let tg_text = if pairs.len() == 1 {
        let tgs: Vec<i128> = if a[0] == a[1] { vec![a[0]] } else { a.to_vec() };
        format!(
            "{}{}",
            tgs.iter()
                .map(|x| format!("TG{x}"))
                .collect::<Vec<_>>()
                .join(" and "),
            if five == 1 {
                " on the 5T copy"
            } else {
                " on the 5T copies"
            }
        )
    } else {
        format!(
            "TG{} or TG{} on the 5T copy, not told by the counts ({})",
            a[0],
            a[1],
            pairs.join(" or ")
        )
    };
    let line = format!(
        "5T: {five} copy{} with the 5T poly-T tract, {tg_text}. 5T lowers the splicing of exon 10 and is a variable-penetrance allele: with TG12 or TG13 on the same copy it more often causes CBAVD or another CFTR-related disorder when the other copy carries a CF-causing variant.",
        if five > 1 { "s" } else { "" }
    );
    let (phase, text) = match &r117 {
        None => (
            "not_read",
            format!(" R117H was not read ({why}), so whether it sits on the 5T copy is not known."),
        ),
        Some(g) => {
            let n = a_copies(g);
            let (p, t) = if n == 0 {
                (
                    "not_needed",
                    " No R117H copy: the 5T is read on its own.".to_owned(),
                )
            } else if n == 2 || five == 2 {
                (
                    "in_cis",
                    format!(
                        " R117H ({g}) and 5T share a copy: R117H with 5T in cis is a CF-causing allele of variable penetrance."
                    ),
                )
            } else {
                (
                    "unknown",
                    " One R117H copy and one 5T copy: phase unknown. R117H lies 17.6 kb from the tract, farther than a read pair spans, so whether they share a copy (in cis, a CF-causing allele of variable penetrance) or not (in trans) is not read; testing a parent tells them apart.".to_owned(),
                )
            };
            (p, t + &note)
        }
    };
    out.set("five_t", Value::Bool(true))
        .set("phase", s(phase))
        .set("r117h", r117_text)
        .set("reading", s(&(line + &text + &odd)))
        .val()
}

fn ugt1a1_name(n: i128) -> Option<&'static str> {
    match n {
        6 => Some("*36"),
        7 => Some("*1"),
        8 => Some("*28"),
        9 => Some("*37"),
        _ => None,
    }
}

fn ugt1a1(c: &Ctx<'_>, by: &IndexMap<String, Sized>, cpic: State) -> Value {
    let out = Obj::new()
        .set("locus", s("UGT1A1"))
        .set("repeat", s("promoter TA repeat"))
        .set("source", s(UGT1A1_SOURCE));
    let not_sized = |o: Obj, reading: String| {
        o.set("status", s("not_sized"))
            .set("genotype", s("."))
            .set("pharmcat", s("."))
            .set("agreement", s("not_compared"))
            .set("reading", s(&reading))
            .val()
    };
    let Some(item) = by.get("UGT1A1") else {
        return not_sized(
            out,
            format!(
                "not in step 09's output: {EXTRA_SETUP} adds the repeat to the catalog (docs/09-str-expansions.md). PharmCAT (step 7) still reads UGT1A1 *28."
            ),
        );
    };
    let Some(cn) = counts(Some(item)) else {
        return not_sized(
            out,
            format!(
                "ExpansionHunter did not size the repeat ({}, {}): not read.",
                text_of(&item.repeat_count),
                sizing_text(item)
            ),
        );
    };
    let names: Vec<String> = cn
        .iter()
        .map(|n| ugt1a1_name(*n).map_or_else(|| format!("TA{}", n - 1), str::to_owned))
        .collect();
    let genotype = names.join("/");
    let dec = names.iter().filter(|n| *n == "*28" || *n == "*37").count();
    let ta = cn
        .iter()
        .map(|n| format!("TA{}", n - 1))
        .collect::<Vec<_>>()
        .join("/");
    let out = out
        .set("status", s("read"))
        .set("genotype", s(&genotype))
        .set(
            "ta_units",
            Value::Array(cn.iter().map(|v| Value::Int(*v)).collect()),
        );
    let mut reading = format!(
        "{genotype} ({ta}, {}/{} TA units as ExpansionHunter counts them): ",
        cn[0], cn[1]
    );
    if !cn.iter().all(|n| ugt1a1_name(*n).is_some()) {
        reading += "a TA length outside TA5 to TA8 is not a named allele here, so its function is not read and it is not compared with PharmCAT (step 7).";
        return out
            .set("pharmcat", s("."))
            .set("agreement", s("not_compared"))
            .set("reading", s(&reading))
            .val();
    }
    reading += &if dec == 0 {
        "no decreased-function TA allele.".to_owned()
    } else {
        format!(
            "{dec} decreased-function allele{} (*28 and *37 are decreased function).",
            if dec > 1 { "s" } else { "" }
        )
    };
    let (mut pc, mut agreement) = (".".to_owned(), "not_compared");
    let readable = !matches!(cpic, State::Failed | State::Unreadable);
    let row = c
        .find1(format!("cpic/{}_phenotypes.tsv", c.s()))
        .filter(|_| readable)
        .and_then(|p| table(&p.path, false).ok())
        .and_then(|t| {
            t.rows
                .into_iter()
                .find(|r| r.get("Gene").flatten() == Some("UGT1A1"))
        });
    if let Some(row) = row {
        let d = super::util::or(&row, "Diplotype", "");
        let d = pycompat::num::py_strip(&d);
        pc = if d.is_empty() {
            ".".into()
        } else {
            d.to_owned()
        };
        let alleles: Vec<&str> = pc.split('/').collect();
        if alleles.len() == 2
            && !pc.contains("possible")
            && super::util::or(&row, "Status", "") != "not called"
        {
            let n = alleles
                .iter()
                .filter(|a| {
                    a.split(['+', ' ', '[', ']'])
                        .any(|x| x == "*28" || x == "*37")
                })
                .count();
            agreement = if n == dec { "agrees" } else { "differs" };
        }
    }
    let old = if cpic == State::Stale && agreement != "not_compared" {
        " (step 27's table is from an earlier run)"
    } else {
        ""
    };
    reading += &match agreement {
        "agrees" => {
            format!(" PharmCAT (step 7) called {pc}: the same number of *28 and *37 copies{old}.")
        }
        "differs" => format!(
            " PharmCAT (step 7) called {pc}, a different number of *28 and *37 copies{old}: read the reads at chr2:233760234 before using either."
        ),
        _ if !readable => format!(
            " Not compared with PharmCAT (step 27's table is {}).",
            other_state(cpic)
        ),
        _ => " Not compared with PharmCAT (no unambiguous UGT1A1 call from step 27).".to_owned(),
    };
    out.set("pharmcat", s(&pc))
        .set("agreement", s(agreement))
        .set("reading", s(&reading))
        .val()
}

/// What the counts of the extra repeats mean, and the 5-HTTLPR line no
/// step can test.
fn readings(
    c: &Ctx<'_>,
    by: &IndexMap<String, Sized>,
    traits: State,
    cpic: State,
) -> Result<Value, String> {
    Ok(Value::Array(vec![
        cftr(c, by, traits),
        ugt1a1(c, by, cpic),
        Obj::new()
            .set("locus", s("SLC6A4"))
            .set("repeat", s("5-HTTLPR promoter length polymorphism"))
            .set("status", s("not_tested"))
            .set("genotype", s("."))
            .set("reading", s(HTTLPR_READING))
            .set("source", s("Bousman et al. 2023, PMID 37032427"))
            .val(),
    ]))
}

/// The post pass: the readings again, with the traits and CPIC states, over
/// the loci as the section lists them.
pub fn post(c: &Ctx<'_>, v: &mut Value, traits: State, cpic: State) {
    let Value::Object(m) = v else { return };
    let loci: Vec<(String, Obj)> = match m.get("loci") {
        Some(Value::Array(xs)) => xs
            .iter()
            .filter_map(|x| match x {
                Value::Object(o) => match o.get("locus") {
                    Some(Value::Str(l)) => Some((l.clone(), Obj(o.clone()))),
                    _ => None,
                },
                _ => None,
            })
            .collect(),
        _ => return,
    };
    let by = by_locus(loci.iter().map(|(k, o)| (k, o)));
    if let Ok(r) = readings(c, &by, traits, cpic) {
        m.insert("readings".into(), r);
    }
}
