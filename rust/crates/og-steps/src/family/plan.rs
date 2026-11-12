//! `og family plan`: checks a families file and picks one family, then
//! says, for each analysis of steps 38 and 39, whether it runs, is skipped
//! or waits for a member's own run, and why.

use std::collections::{BTreeSet, HashMap, HashSet};
use std::path::Path;

use og_pedigree::{Die, Members, build_members, family_rows, pedigree_text};
use pycompat::fmt::ljust;

use super::args::{Args, Kind, Spec};
use super::{Result, print_out, write_file};

pub const SPEC: Spec = Spec {
    prog: "og family plan",
    usage: "usage: og family plan [-h] --families FAMILIES [--family FAMILY] [--sample SAMPLE] [--sex SEX]\n                      [--have-gvcf [ID ...]] [--have-clinvar [ID ...]] [--ibd-panel {yes,no}]\n                      [--write-pedigree DIR] [--which] [--out OUT]\n",
    options: &[
        ("--families", Kind::One),
        ("--family", Kind::One),
        ("--sample", Kind::One),
        ("--sex", Kind::One),
        ("--have-gvcf", Kind::Star),
        ("--have-clinvar", Kind::Star),
        ("--ibd-panel", Kind::Choice(&["yes", "no"], "no")),
        ("--write-pedigree", Kind::One),
        ("--which", Kind::Flag),
        ("--out", Kind::One),
    ],
    required: &["--families"],
};

/// The columns of `family_plan.tsv`.
pub const PLAN_COLUMNS: [&str; 5] = ["family", "analysis", "members", "decision", "reason"];

/// One row of `family_plan.tsv`.
pub struct PlanRow {
    pub family: String,
    pub analysis: &'static str,
    pub members: String,
    pub decision: &'static str,
    pub reason: String,
}

/// The rows of `family_plan.tsv`: for each analysis of steps 38 and 39,
/// whether it runs, is skipped, or waits for a member's own run, and why.
/// `have_gvcf` and `have_clinvar` are the samples with a gVCF and with
/// step 06's ClinVar hits.
pub fn plan_family(
    name: &str,
    members: &Members,
    have_gvcf: &HashSet<String>,
    have_clinvar: &HashSet<String>,
    ibd_panel: bool,
) -> Vec<PlanRow> {
    let rows_of: Vec<&str> = members
        .iter()
        .filter(|(_, m)| m.row)
        .map(|(s, _)| s.as_str())
        .collect::<BTreeSet<&str>>()
        .into_iter()
        .collect();
    let seq: Vec<&str> = rows_of
        .iter()
        .copied()
        .filter(|s| have_gvcf.contains(*s))
        .collect();
    let lacking: Vec<&str> = rows_of
        .iter()
        .copied()
        .filter(|s| !have_gvcf.contains(*s))
        .collect();
    let no_gvcf = |ss: &[&str]| -> String {
        ss.iter()
            .map(|s| format!("{s} has no gVCF yet: run run-all.sh {s} first"))
            .collect::<Vec<_>>()
            .join("; ")
    };
    let mut out = Vec::new();
    let mut row = |analysis: &'static str, who: &[&str], decision: &'static str, reason: String| {
        out.push(PlanRow {
            family: name.to_owned(),
            analysis,
            members: if who.is_empty() {
                ".".to_owned()
            } else {
                who.join(",")
            },
            decision,
            reason,
        });
    };

    let joint: (&'static str, String) = if rows_of.len() < 2 {
        (
            "skipped",
            format!(
                "the family has {} sample with a row ({}); the family analysis needs two or more",
                rows_of.len(),
                if rows_of.is_empty() {
                    "none".to_owned()
                } else {
                    rows_of.join(", ")
                }
            ),
        )
    } else if seq.len() < 2 {
        ("waiting", no_gvcf(&lacking))
    } else {
        let mut why = format!(
            "GLnexus joint call, somalier relationship check and the family report of {}",
            seq.join(", ")
        );
        if !lacking.is_empty() {
            why.push_str(&format!(
                "; left out: {} (run step 38 again after its run)",
                no_gvcf(&lacking)
            ));
        }
        ("runs", why)
    };
    let runs = joint.0 == "runs";
    row(
        "joint_call",
        if runs { &seq } else { &rows_of },
        joint.0,
        joint.1.clone(),
    );

    let in_seq = |p: &Option<String>| p.as_deref().is_some_and(|p| seq.contains(&p));
    let trios: Vec<&str> = seq
        .iter()
        .copied()
        .filter(|s| in_seq(&members[*s].father) && in_seq(&members[*s].mother))
        .collect();
    let duos: Vec<&str> = seq
        .iter()
        .copied()
        .filter(|s| {
            !trios.contains(s) && (in_seq(&members[*s].father) || in_seq(&members[*s].mother))
        })
        .collect();
    if !runs {
        for a in ["trio_phasing", "inheritance_map", "shared_segments"] {
            row(a, &rows_of, joint.0, joint.1.clone());
        }
    } else {
        if !trios.is_empty() {
            row(
                "trio_phasing",
                &trios,
                "runs",
                format!(
                    "WhatsHap phases {} by inheritance from both parents",
                    trios.join(", ")
                ),
            );
        } else {
            let has_one = if duos.is_empty() {
                String::new()
            } else {
                format!(" ({} has one)", duos.join(", "))
            };
            row(
                "trio_phasing",
                &duos,
                "skipped",
                format!(
                    "no child has both parents sequenced{has_one}; phasing by inheritance needs a trio"
                ),
            );
        }
        let mut kids: Vec<&str> = trios.iter().chain(duos.iter()).copied().collect();
        kids.sort_unstable();
        if !kids.is_empty() {
            row(
                "inheritance_map",
                &kids,
                "runs",
                format!(
                    "the allele each sequenced parent passed on, and the chromosome map, of {}",
                    kids.join(", ")
                ),
            );
        } else {
            row(
                "inheritance_map",
                &[],
                "skipped",
                "no child has a sequenced parent".to_owned(),
            );
        }
        if ibd_panel {
            row(
                "shared_segments",
                &seq,
                "runs",
                "Beagle phases the joint call against step 26b's panel and hap-ibd lists the DNA each pair shares"
                    .to_owned(),
            );
        } else {
            row(
                "shared_segments",
                &seq,
                "skipped",
                "no phasing panel: install it with setup.sh --ancestry-fine-panel".to_owned(),
            );
        }
    }

    let couples: BTreeSet<(String, String)> = rows_of
        .iter()
        .filter_map(|s| {
            let p = members[*s].partner.as_deref()?;
            members.get(p).filter(|pm| pm.row)?;
            let (a, b) = if *s < p { (*s, p) } else { (p, *s) };
            Some((a.to_owned(), b.to_owned()))
        })
        .collect();
    if couples.is_empty() {
        row(
            "couple_carrier",
            &[],
            "skipped",
            "no two samples of the family are marked as partners".to_owned(),
        );
    }
    // Step 39 reads each partner's ClinVar hits, never a gVCF: a partner run
    // from a VCF alone still gets the screen.
    for (a, b) in &couples {
        let pair = [a.as_str(), b.as_str()];
        let lack: Vec<&str> = pair
            .iter()
            .copied()
            .filter(|s| !have_clinvar.contains(*s))
            .collect();
        if !lack.is_empty() {
            let why: Vec<String> = lack
                .iter()
                .map(|s| {
                    format!(
                        "{s} has no ClinVar hits yet (step 06 writes clinvar/{s}_clinvar_hits.vcf): run run-all.sh {s} first"
                    )
                })
                .collect();
            row("couple_carrier", &pair, "waiting", why.join("; "));
        } else {
            row(
                "couple_carrier",
                &pair,
                "runs",
                format!("carrier screen of {a} and {b}"),
            );
        }
    }
    out
}

/// `family_plan.tsv` as text.
pub fn plan_text(rows: &[PlanRow]) -> String {
    let mut s = PLAN_COLUMNS.join("\t");
    s.push('\n');
    for r in rows {
        s.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\n",
            r.family, r.analysis, r.members, r.decision, r.reason
        ));
    }
    s
}

pub fn run(a: &Args) -> Result<()> {
    let families = a.one("--families").unwrap_or_default();
    let sample = a.one("--sample");
    let sex = a
        .one("--sex")
        .map(|s| pycompat::text::lower(pycompat::text::strip(s, None)))
        .unwrap_or_default();
    if !sex.is_empty() && sex != "male" && sex != "female" {
        return Err(Die::one(format!(
            "--sex {}: the sex is male or female",
            a.one("--sex").unwrap_or_default()
        ))
        .into());
    }
    if !sex.is_empty() && sample.is_none_or(str::is_empty) {
        return Err(Die::one("--sex needs --sample").into());
    }
    let sample = sample.filter(|s| !s.is_empty());
    let family = a.one("--family").filter(|s| !s.is_empty());
    let (name, rows) = family_rows(families, family, sample)?;
    if let Some(s) = sample
        && rows.is_empty()
    {
        if a.flag("--which") {
            return Ok(());
        }
        return Err(Die::one(format!("{s} is in no family of {families}")).into());
    }
    let name = name.unwrap_or_default();
    let mut declared = HashMap::new();
    if !sex.is_empty()
        && let Some(s) = sample
    {
        declared.insert(s.to_owned(), sex);
    }
    let (members, errors) =
        build_members(&rows, &declared, families, "the sex given to run-all.sh");
    if !errors.is_empty() {
        return Err(Die::new(errors).into());
    }
    if let Some(s) = sample
        && !rows.is_empty()
        && !members.contains_key(s)
    {
        return Err(Die::one(format!("{families}: '{s}' is not in family '{name}'")).into());
    }
    if a.flag("--which") {
        print_out(&format!("{name}\n"));
        return Ok(());
    }
    if let Some(dir) = a.one("--write-pedigree") {
        let d = Path::new(dir).join(&name);
        std::fs::create_dir_all(&d)
            .map_err(|e| Die::one(format!("cannot create {}: {e}", d.display())))?;
        write_file(&d.join("pedigree.tsv"), pedigree_text(&rows).as_bytes())?;
    }
    let set = |k: &str| -> HashSet<String> { a.many(k).iter().cloned().collect() };
    let plan = plan_family(
        &name,
        &members,
        &set("--have-gvcf"),
        &set("--have-clinvar"),
        a.one("--ibd-panel") == Some("yes"),
    );
    if let Some(out) = a.one("--out") {
        write_file(Path::new(out), plan_text(&plan).as_bytes())?;
    }
    let mut text = String::new();
    for r in &plan {
        text.push_str(&format!(
            "{} {} {} {}\n",
            r.family,
            ljust(r.analysis, 16),
            ljust(r.decision, 8),
            r.reason
        ));
    }
    print_out(&text);
    Ok(())
}
