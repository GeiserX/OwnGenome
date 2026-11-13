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
//! - `recheck`: the last database re-check of the sample;
//! - `new_since`: what changed since the previous run kept in `history/`.
//!
//! Each part that changes a section runs on its own: when it fails, only
//! its section becomes `unreadable`, with a note naming the file, and the
//! summary is still written. The Python collector ran this pass outside
//! its per-section error handling, so one failure there lost summary.json
//! and both reports. A previous run that cannot be read leaves `new_since`
//! out, as before; the comparison itself has no error to hide.

use super::group_c::Others;
use super::model::{Section, State, Summary};
use super::sections::Data;
use super::{Folder, Input, declared, group_b, post_c, qc_declared, read_a, read_c, time};
use og_runlog::{delta, files};
use pycompat::json::Value;
use pycompat::value::truthy;
use std::path::Path;

/// summary.json's sections, in its order.
pub const ORDER: [&str; 48] = [
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
/// does; the couple screen and the carrier list come with other steps.
const NOT_A_STEP: [&str; 3] = ["read_count", "couple", "carrier"];

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

/// Reads the whole summary of a sample folder. The three groups are read
/// at once over one walk of the folder. An error is a run status or a
/// manifest that is there and cannot be read.
pub fn collect(input: &Input) -> Result<Value, String> {
    let folder = Folder::open(input)?;
    let (a, b, mut c) = std::thread::scope(|scope| {
        let a = scope.spawn(|| read_a(input, &folder));
        let b = scope.spawn(|| group_b::read(input, &folder));
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
        let a = a
            .join()
            .unwrap_or_else(|_| stopped(&super::model::GROUP_A));
        let b = b.join().unwrap_or_else(|_| stopped(&group_b::GROUP_B));
        (a, b, c)
    });
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
    let dir = Path::new(&input.sample_dir);
    if let Some(r) = delta::read_recheck(&dir.join("recheck").join("state.json")) {
        m.insert("recheck".into(), Value::Object(r));
    }
    if let Some(n) = new_since(dir, &v) {
        if let Value::Object(m) = &mut v {
            m.insert("new_since".into(), n);
        }
    }
    Ok(v)
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
                Some(Data::Cyrius(c)) => Some((s.state, Some(c.genotype.clone()), c.filter.clone())),
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
        out.push(PRS_NOT_ADJUSTED);
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
