//! The whole summary.json: every section of groups a, b and c, read over
//! one walk of the folder, then the post pass.
//!
//! The post pass is what joins sections after each is read:
//!
//! - the repeat readings with the traits and CPIC states, Locityper's
//!   primary calls with HLA, pypgx and Cyrius, and CNVpytor's chrY filter
//!   with the declared sex ([`super::post_c`]);
//! - the `cyp2d6` block over CPIC, pypgx and Cyrius;
//! - `not_run`: the steps whose section has no file;
//! - `not_assessed`: what no step looked at, with the lines that depend on
//!   the polygenic scores, HBA1/HBA2, mobile elements, VNtyper and MoChA;
//! - `scope`: what part of the genome the run called, and the `partial`
//!   mark on each genome-wide section of a run that called less, and on
//!   each stale one, whose earlier run's scope is not known;
//! - `recheck`: the last database re-check of the sample;
//! - `new_since`: what changed since the previous run kept in `history/`.
//!
//! Each part that changes a section runs on its own: when it returns an
//! error, only its section becomes `unreadable`, with a note naming the
//! file, and the summary is still written. The parts read with defaults, so
//! a row of an unexpected shape is no error; a panic is not caught, and it
//! ends og with an error. The Python collector ran this pass outside
//! its per-section error handling, so one failure there lost summary.json
//! and both reports. A previous run that cannot be read leaves `new_since`
//! out, as before; the comparison itself has no error to hide.

use super::group_c::Others;
use super::model::{Section, State, Summary};
use super::sections::Data;
use super::{Folder, Input, declared, group_b, post_c, qc_declared, read_a, read_c, scan, time};
use og_runlog::{delta, files};
use pycompat::json::Value;
use pycompat::value::truthy;
use std::path::Path;

/// summary.json's sections, in its order.
pub const ORDER: [&str; 49] = [
    "coverage",
    "callability",
    "sex_check",
    "sample_qc",
    "mosaic_sex_loss",
    "read_count",
    "variants",
    "clinvar",
    "couple",
    "carrier",
    "acmg_sf",
    "pharmcat",
    "cpic",
    "pypgx",
    "cyrius",
    "parascopy",
    "paralogs",
    "hla",
    "prs",
    "ancestry_fine",
    "manta",
    "delly",
    "cnvpytor",
    "sv_consensus",
    "annotsv",
    "expansions",
    "telomere",
    "roh",
    "haplogroup",
    "y_haplogroup",
    "y_tree",
    "y_str",
    "consumer_export",
    "bloodtype",
    "lpa",
    "vntyper",
    "locityper",
    "mito",
    "traits",
    "archaic",
    "ancient",
    "viral",
    "exomiser",
    "chip",
    "mei",
    "mosaic",
    "cpsr",
    "clinical",
    "slivar",
];

/// Sections that are no step of their own, so never in `not_run`: the read
/// count is written by alignment, which a run from a BAM or a VCF never
/// does; the couple screen, the carrier list and the ACMG secondary-findings
/// card come with other steps (39 and 06), whose own sections say not run.
const NOT_A_STEP: [&str; 4] = ["read_count", "couple", "carrier", "acmg_sf"];

/// What no step assesses.
const NOT_ASSESSED: [&str; 5] = [
    "Mosaic and low-fraction small variants (only chrM heteroplasmy is reported, and with step 49 the CHIP genes)",
    "Methylation and imprinting",
    "Repeat expansions outside the ExpansionHunter catalog",
    "Phase: two variants in one gene are compound-het candidates, not confirmed",
    "Variants of uncertain significance (not interpreted)",
];

/// Added unless step 25 ran with an ancestry reference panel.
const PRS_NOT_ADJUSTED: &str =
    "Polygenic score percentiles: no ancestry reference panel was used, so the scores are raw sums";

/// Added in its place when the panel ran but every percentile was withheld,
/// because the input lists variants only (a score's `percentile_note`).
const PRS_WITHHELD: &str = "Polygenic score percentiles: the input lists variants only, so the \
percentiles were withheld and the scores are raw sums";

/// Added unless step 35's depth count (short reads) or step 35b
/// (Paraphase, PacBio HiFi, both genes) made a call.
const HBA_NOT_ASSESSED: &str = "HBA1/HBA2 copy number (alpha-thalassaemia deletions): neither step 35's read-depth count \
(short reads) nor step 35b (Paraphase, PacBio HiFi) made a call";

/// Added when step 50 (opt-in, MEGAnE) did not run. A failed step 50 is
/// listed under not run, and an unreadable table under its section.
const MEI_NOT_ASSESSED: &str = "Mobile element insertions (Alu, L1 and SVA, the FKTN, TAF1 and BBS1 SVA insertions among them): \
not looked for (step 50, MEGAnE, is opt-in)";

/// Added unless step 52 (VNtyper) gave a positive, possible or negative
/// reading: no other step sees a frameshift inside the MUC1 VNTR.
const VNTYPER_NOT_ASSESSED: &str = "MUC1 VNTR frameshifts (ADTKD-MUC1): step 52 (VNtyper, opt-in) gave no reading; \
no standard caller sees them";

/// Added unless step 51 (MoChA) assessed the sample.
const MCA_NOT_ASSESSED: &str = "Autosomal mosaic chromosomal alterations (large gains, losses and copy-neutral LOH in a share \
of the blood cells): step 51 (MoChA) did not assess this sample";

/// Reads the whole summary of a sample folder, over one walk of the
/// folder. An error is a run status or a manifest that is there and cannot
/// be read.
pub fn collect(input: &Input) -> Result<Value, String> {
    let folder = Folder::open(input)?;
    let vcf = folder
        .index
        .first(&[format!(
            "vcf/{}.vcf.gz",
            pycompat::mtime::escape(&input.sample)
        )])
        .map(|rel| crate::join(&input.sample_dir, rel));
    // Group a reads the sample VCF, the input that takes time, on its own
    // thread, and so does the read of its contigs; groups b and c read
    // their sections meanwhile, one at a time, so the peak memory stays
    // near the biggest section's.
    let (a, b, mut c, contigs) = std::thread::scope(|scope| {
        let contigs = scope.spawn(|| vcf.as_deref().and_then(|p| scan::vcf_contigs(p).ok()));
        let a = scope.spawn(|| read_a(input, &folder));
        let b = group_b::read(input, &folder);
        let c = read_c(input, &folder);
        let stopped = |specs: &[super::model::Spec]| -> Vec<Section> {
            specs
                .iter()
                .map(|s| {
                    Section::new(
                        *s,
                        Err("internal error: the section reader stopped".into()),
                        folder.run.as_ref(),
                    )
                })
                .collect()
        };
        let a = a.join().unwrap_or_else(|_| stopped(&super::model::GROUP_A));
        (a, b, c, contigs.join().ok().flatten())
    });
    let intervals = folder
        .manifest
        .as_ref()
        .and_then(|m| m.value("run", "intervals"))
        .map(str::trim)
        .filter(|v| !v.is_empty())
        .map(str::to_owned);
    let declared_sex = declared(input, &folder, qc_declared(&a));
    post_c(
        input,
        &folder.index,
        &mut c,
        &others(&b),
        declared_sex.as_deref(),
    );
    let cyp2d6 = group_b::cyp2d6(input, &folder.index, &b);
    let sections = in_order(a.into_iter().chain(b).chain(c).collect());
    let not_run = not_run(&sections);
    let not_assessed = not_assessed(&sections);
    let now = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_or(0.0, |d| d.as_secs_f64());
    let summary = Summary {
        sample: input.sample.clone(),
        generated_utc: time::utc_minute(now).unwrap_or_default(),
        run: folder.run,
        manifest: folder.manifest,
        declared_sex,
        sections,
        cyp2d6,
    };
    let mut v = summary.to_json();
    drop(summary);
    let Value::Object(m) = &mut v else {
        return Err("internal error: the summary is not an object".into());
    };
    m.insert(
        "not_run".into(),
        Value::Array(not_run.into_iter().map(Value::Str).collect()),
    );
    m.insert(
        "not_assessed".into(),
        Value::Array(
            not_assessed
                .into_iter()
                .map(|x| Value::Str(x.to_owned()))
                .collect(),
        ),
    );
    if let Some(scope) = scope(contigs.as_deref(), intervals.as_deref()) {
        mark_partial(m, &scope);
        m.insert("scope".into(), scope);
    }
    let dir = Path::new(&input.sample_dir);
    if let Some(r) = delta::read_recheck(&dir.join("recheck").join("state.json")) {
        m.insert("recheck".into(), Value::Object(r));
    }
    if let Some(n) = new_since(dir, &v)
        && let Value::Object(m) = &mut v
    {
        m.insert("new_since".into(), n);
    }
    Ok(v)
}

/// The sections whose result is a share, a sum or a placement over the
/// whole genome: a run over a part of it gives a different result that
/// reads the same.
pub const GENOME_WIDE: [&str; 6] = [
    "callability",
    "prs",
    "ancestry_fine",
    "telomere",
    "roh",
    "archaic",
];

/// `scope`: what part of the genome the run called, from the contigs of
/// the sample VCF's records (`chr20` and `20` alike) and the intervals the
/// run manifest records (`run intervals`). `None` (no `scope` key) without a
/// sample VCF read to its end: the scope is then not known, and nothing
/// is marked.
///
/// The run is partial when the calls cover fewer than the 22 autosomes, or
/// when the run was limited to intervals. chrX is recorded but makes no run
/// partial: the genome-wide sections read the autosomes.
fn scope(contigs: Option<&[String]>, intervals: Option<&str>) -> Option<Value> {
    let names: Vec<&str> = contigs?
        .iter()
        .map(|c| c.strip_prefix("chr").unwrap_or(c))
        .collect();
    let autosomes: Vec<String> = (1..=22)
        .map(|i| i.to_string())
        .filter(|n| names.contains(&n.as_str()))
        .map(|n| format!("chr{n}"))
        .collect();
    let n = autosomes.len();
    let mark = if n < 22 {
        Some(format!("partial: covers {n} of 22 autosomes"))
    } else {
        intervals.map(|i| format!("partial: limited to the run's intervals ({i})"))
    };
    let mut out = indexmap::IndexMap::new();
    out.insert("autosomes_covered".to_owned(), Value::Int(n as i128));
    out.insert(
        "autosomes".to_owned(),
        Value::Array(autosomes.into_iter().map(Value::Str).collect()),
    );
    out.insert("chrx".to_owned(), Value::Bool(names.contains(&"X")));
    out.insert(
        "intervals".to_owned(),
        intervals.map_or(Value::Null, |i| Value::Str(i.to_owned())),
    );
    out.insert("partial".to_owned(), Value::Bool(mark.is_some()));
    out.insert("mark".to_owned(), mark.map_or(Value::Null, Value::Str));
    Some(Value::Object(out))
}

/// The mark of a stale genome-wide section: its result is from an earlier
/// run, and that run's scope was not recorded.
const STALE_MARK: &str = "partial: not known, from an earlier run";

/// Adds a `partial` mark to the genome-wide sections. A result of this run
/// (ok) takes the scope's mark when the run was partial. A stale result is
/// from an earlier run, whose scope this run's VCF does not tell: it takes
/// `STALE_MARK` whatever this run covered, so it never reads as complete.
fn mark_partial(summary: &mut indexmap::IndexMap<String, Value>, scope: &Value) {
    let mark = match at(Some(scope), "mark") {
        Some(Value::Str(m)) => Some(m.clone()),
        _ => None,
    };
    let Some(Value::Object(sections)) = summary.get_mut("sections") else {
        return;
    };
    for key in GENOME_WIDE {
        let Some(Value::Object(sec)) = sections.get_mut(key) else {
            continue;
        };
        let text = if is(sec.get("state"), "stale") {
            Some(STALE_MARK.to_owned())
        } else if is(sec.get("state"), "ok") {
            mark.clone()
        } else {
            None
        };
        if let Some(text) = text {
            sec.insert("partial".into(), Value::Str(text));
        }
    }
}

/// What the post pass of group c needs from group b: the CPIC state,
/// pypgx's CYP2D6 call and Cyrius's genotype and Filter.
fn others(b: &[Section]) -> Others {
    let sec = |k: &str| b.iter().find(|s| s.spec.key == k);
    let live = |s: &Section| matches!(s.state, State::Ok | State::Stale);
    Others {
        cpic: sec("cpic").map(|s| s.state),
        traits: None,
        pypgx: sec("pypgx")
            .filter(|s| live(s))
            .and_then(|s| match &s.data {
                Some(Data::Pypgx(p)) => Some((s.state, p.cyp2d6.clone())),
                _ => None,
            }),
        cyrius: sec("cyrius")
            .filter(|s| live(s))
            .and_then(|s| match &s.data {
                Some(Data::Cyrius(c)) => {
                    Some((s.state, Some(c.genotype.clone()), c.filter.clone()))
                }
                _ => None,
            }),
    }
}

/// The sections in summary.json's order ([`ORDER`]).
fn in_order(mut sections: Vec<Section>) -> Vec<Section> {
    sections.sort_by_key(|s| ORDER.iter().position(|k| *k == s.spec.key));
    sections
}

/// The steps whose section has no file, with " (step N failed)" for a step
/// that failed in the last run.
fn not_run(sections: &[Section]) -> Vec<String> {
    sections
        .iter()
        .filter(|s| matches!(s.state, State::Missing | State::Failed))
        .filter(|s| !NOT_A_STEP.contains(&s.spec.key))
        .map(|s| {
            if s.state == State::Failed {
                format!("{} (step {} failed)", s.spec.title, s.spec.step)
            } else {
                s.spec.title.to_owned()
            }
        })
        .collect()
}

/// What no step assessed for this sample.
fn not_assessed(sections: &[Section]) -> Vec<&'static str> {
    let sec = |k: &str| sections.iter().find(|s| s.spec.key == k);
    let state = |k: &str| sec(k).map_or(State::Missing, |s| s.state);
    // A section's data is there when it is ok or stale; else it reads {}.
    let field = |k: &str, key: &str| -> Option<&Value> {
        match sec(k).and_then(|s| s.data.as_ref()) {
            Some(Data::Json(Value::Object(m))) => m.get(key),
            _ => None,
        }
    };
    let mut out = NOT_ASSESSED.to_vec();
    if !field("prs", "adjusted").is_some_and(truthy) {
        let withheld = match field("prs", "scores") {
            Some(Value::Array(rows)) => rows
                .iter()
                .any(|r| at(Some(r), "percentile_note").is_some_and(truthy)),
            _ => false,
        };
        out.push(if withheld {
            PRS_WITHHELD
        } else {
            PRS_NOT_ADJUSTED
        });
    }
    // Step 35b on a region subset (--paraphase_regions) may leave hba out.
    let hba_genes: Vec<&str> = match field("paralogs", "genes") {
        Some(Value::Array(genes)) => genes
            .iter()
            .filter(|g| is(at(Some(g), "region"), "hba") && is(at(Some(g), "status"), "called"))
            .filter_map(|g| match at(Some(g), "gene") {
                Some(Value::Str(x)) => Some(x.as_str()),
                _ => None,
            })
            .collect(),
        _ => Vec::new(),
    };
    let paraphase_hba = ["HBA1", "HBA2"].iter().all(|g| hba_genes.contains(g));
    if !(paraphase_hba || is(at(field("parascopy", "hba"), "status"), "called")) {
        out.push(HBA_NOT_ASSESSED);
    }
    if state("mei") == State::Missing {
        out.push(MEI_NOT_ASSESSED);
    }
    let vntyper_read = state("vntyper") == State::Ok
        && ["positive", "possible", "negative"]
            .iter()
            .any(|w| is(field("vntyper", "status"), w));
    if !vntyper_read {
        out.push(VNTYPER_NOT_ASSESSED);
    }
    if !is(field("mosaic", "status"), "assessed") {
        out.push(MCA_NOT_ASSESSED);
    }
    out
}

/// True when `v` is the string `word`.
fn is(v: Option<&Value>, word: &str) -> bool {
    matches!(v, Some(Value::Str(x)) if x == word)
}

/// `v[key]` when `v` is an object that has it.
fn at<'a>(v: Option<&'a Value>, key: &str) -> Option<&'a Value> {
    match v {
        Some(Value::Object(m)) => m.get(key),
        _ => None,
    }
}

/// `{previous_utc, sections}`: what changed since the newest run kept in
/// `history/` that holds a summary.json (generate-report.sh keeps one per
/// run), or `None` without one. A previous run whose summary or manifest
/// cannot be read, or whose summary is not shaped as a summary, counts as
/// none.
fn new_since(dir: &Path, cur: &Value) -> Option<Value> {
    let prev_dir = delta::newest_history(dir)?;
    let bytes = std::fs::read(prev_dir.join("summary.json")).ok()?;
    let prev = pycompat::json::loads(std::str::from_utf8(&bytes).ok()?).ok()?;
    let prev_manifest = files::read_manifest(Some(&prev_dir.join("run_manifest.tsv"))).ok()?;
    let manifest = files::read_manifest(Some(&dir.join("run_manifest.tsv"))).ok()?;
    let changes = delta::compare(&prev, &prev_manifest, cur, &manifest).ok()?;
    let previous = match &prev {
        Value::Object(m) => m.get("generated_utc").cloned().unwrap_or(Value::Null),
        _ => Value::Null,
    };
    let mut out = indexmap::IndexMap::new();
    out.insert("previous_utc".to_owned(), previous);
    out.insert(
        "sections".to_owned(),
        Value::Array(
            delta::changed_sections(&changes)
                .into_iter()
                .map(Value::Str)
                .collect(),
        ),
    );
    Some(Value::Object(out))
}

#[cfg(test)]
mod tests;
