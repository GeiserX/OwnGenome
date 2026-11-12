//! `og family prepare`: reads and checks the pedigree, or one family of a
//! families file, and checks that each gVCF names its own sample. Writes
//! `members.tsv` (every person of the pedigree), `family.ped` (somalier
//! relate's pedigree: the sequenced members), `trios.ped` (WhatsHap's: the
//! children whose two parents are sequenced) and `family_plan.tsv`.

use std::collections::{HashMap, HashSet};
use std::path::Path;

use indexmap::IndexMap;
use og_pedigree::{Die, Kinship, build_members, declared_relationship, family_rows, pedigree_text};
use pycompat::fmt::fixed;

use super::args::{Args, Kind, Spec};
use super::plan::{plan_family, plan_text};
use super::{Result, print_out, write_file};

pub const SPEC: Spec = Spec {
    prog: "og family prepare",
    usage: "usage: og family prepare [-h] --pedigree PEDIGREE [--gvcf ID=PATH] [--declared-sex ID=SEX] [--family FAMILY]\n                         [--have-clinvar [ID ...]] [--ibd-panel {yes,no}] [--pedigree-out PEDIGREE_OUT]\n                         --out-dir OUT_DIR\n",
    options: &[
        ("--pedigree", Kind::One),
        ("--gvcf", Kind::Append),
        ("--declared-sex", Kind::Append),
        ("--family", Kind::One),
        ("--have-clinvar", Kind::Star),
        ("--ibd-panel", Kind::Choice(&["yes", "no"], "no")),
        ("--pedigree-out", Kind::One),
        ("--out-dir", Kind::One),
    ],
    required: &["--pedigree", "--out-dir"],
};

/// The sample columns of a VCF's `#CHROM` line (reads the header only);
/// `None` when it has none.
fn gvcf_samples(path: &str) -> std::result::Result<Option<Vec<String>>, Die> {
    let r = og_io::open(path).map_err(|e| Die::io(path, &e))?;
    let mut lines = og_io::Lines::new(r, og_io::Utf8::Replace);
    let mut line = String::new();
    while lines.read_line(&mut line).map_err(|e| Die::io(path, &e))? {
        if line.starts_with("#CHROM") {
            let l = line.trim_end_matches('\n');
            return Ok(Some(l.split('\t').skip(9).map(str::to_owned).collect()));
        }
        if !line.starts_with('#') {
            break;
        }
    }
    Ok(None)
}

pub fn run(a: &Args) -> Result<()> {
    let pedigree = a.one("--pedigree").unwrap_or_default();
    let mut declared = HashMap::new();
    for kv in a.many("--declared-sex") {
        let (k, v) = kv.split_once('=').unwrap_or((kv, ""));
        let v = pycompat::text::lower(pycompat::text::strip(v, None));
        if !v.is_empty() {
            if v != "male" && v != "female" {
                return Err(
                    Die::one(format!("--declared-sex {kv}: the sex is male or female")).into(),
                );
            }
            declared.insert(k.to_owned(), v);
        }
    }
    let family = a.one("--family").filter(|s| !s.is_empty());
    let (name, rows) = family_rows(pedigree, family, None)?;
    let name = name.unwrap_or_default();
    let (mut members, mut errors) = build_members(&rows, &declared, pedigree, "the samplesheet");
    let mut gvcfs: IndexMap<String, String> = IndexMap::new();
    for kv in a.many("--gvcf") {
        let (k, v) = kv.split_once('=').unwrap_or((kv, ""));
        gvcfs.insert(k.to_owned(), v.to_owned());
    }
    let mut seq: Vec<String> = gvcfs
        .keys()
        .filter(|s| members.get(*s).is_some_and(|m| m.row))
        .cloned()
        .collect();
    seq.sort();
    if seq.len() < 2 {
        errors.push(format!(
            "the pedigree names {} of the processed samples ({}); the family analysis needs two \
             or more, each with a row in {pedigree} and a gVCF",
            seq.len(),
            if seq.is_empty() {
                "none".to_owned()
            } else {
                seq.join(", ")
            }
        ));
    }
    for s in &seq {
        let path = &gvcfs[s];
        let names = gvcf_samples(path)?;
        if names.as_deref() != Some(std::slice::from_ref(s)) {
            let base = Path::new(path)
                .file_name()
                .map(|n| n.to_string_lossy().into_owned())
                .unwrap_or_default();
            let (plural, listed) = match names.as_deref() {
                Some(n) if !n.is_empty() => (if n.len() > 1 { "s" } else { "" }, n.join(", ")),
                _ => ("", "none (no #CHROM line)".to_owned()),
            };
            errors.push(format!(
                "the gVCF of '{s}' ({base}) holds the sample{plural} {listed}, not '{s}'. The joint call names each \
                 column after the gVCF's own sample. Rename it: echo {s} > name.txt; bcftools reheader \
                 -s name.txt -o {s}.g.vcf.gz {base}; bcftools index -t {s}.g.vcf.gz"
            ));
        }
    }
    if !errors.is_empty() {
        return Err(Die::new(errors).into());
    }
    for s in &seq {
        if let Some(m) = members.get_mut(s) {
            m.sequenced = true;
        }
    }
    let out_dir = Path::new(a.one("--out-dir").unwrap_or_default());
    std::fs::create_dir_all(out_dir)
        .map_err(|e| Die::one(format!("cannot create {}: {e}", out_dir.display())))?;

    let mut names: Vec<&String> = members.keys().collect();
    names.sort();
    let mut text = String::from("sample\tsex\tfather\tmother\tpartner\tin_pedigree\tsequenced\n");
    for s in names {
        let m = &members[s];
        let o = |v: &Option<String>| v.clone().unwrap_or_default();
        let yn = |b: bool| if b { "yes" } else { "no" };
        text.push_str(&format!(
            "{s}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            o(&m.sex),
            o(&m.father),
            o(&m.mother),
            o(&m.partner),
            yn(m.row),
            yn(m.sequenced)
        ));
    }
    write_file(&out_dir.join("members.tsv"), text.as_bytes())?;

    let in_seq = |p: &Option<String>| p.as_deref().is_some_and(|p| seq.iter().any(|s| s == p));
    let code = |sex: &Option<String>, none: &'static str| match sex.as_deref() {
        Some("male") => "1",
        Some("female") => "2",
        _ => none,
    };
    let mut ped = String::new();
    for s in &seq {
        let m = &members[s];
        let parent = |p: &Option<String>| {
            if in_seq(p) {
                p.clone().unwrap_or_default()
            } else {
                "-9".to_owned()
            }
        };
        ped.push_str(&format!(
            "family\t{s}\t{}\t{}\t{}\t-9\n",
            parent(&m.father),
            parent(&m.mother),
            code(&m.sex, "-9")
        ));
    }
    write_file(&out_dir.join("family.ped"), ped.as_bytes())?;
    let trios: Vec<&String> = seq
        .iter()
        .filter(|s| in_seq(&members[*s].father) && in_seq(&members[*s].mother))
        .collect();
    let mut tp = String::new();
    for s in &trios {
        let m = &members[*s];
        tp.push_str(&format!(
            "family\t{s}\t{}\t{}\t{}\t0\n",
            m.father.as_deref().unwrap_or_default(),
            m.mother.as_deref().unwrap_or_default(),
            code(&m.sex, "0")
        ));
    }
    write_file(&out_dir.join("trios.ped"), tp.as_bytes())?;
    if let Some(p) = a.one("--pedigree-out") {
        write_file(Path::new(p), pedigree_text(&rows).as_bytes())?;
    }
    // The same decision plan writes before the step runs (run-all.sh).
    let have: HashSet<String> = seq.iter().cloned().collect();
    let clinvar: HashSet<String> = a.many("--have-clinvar").iter().cloned().collect();
    let plan = plan_family(
        &name,
        &members,
        &have,
        &clinvar,
        a.one("--ibd-panel") == Some("yes"),
    );
    write_file(
        &out_dir.join("family_plan.tsv"),
        plan_text(&plan).as_bytes(),
    )?;

    let mut say = format!(
        "Family: {} people in the pedigree, {} sequenced ({})\n",
        members.len(),
        seq.len(),
        seq.join(", ")
    );
    let mut kin = Kinship::new(&members);
    for (i, x) in seq.iter().enumerate() {
        for y in &seq[i + 1..] {
            let e = 2.0 * kin.phi(x, y);
            say.push_str(&format!(
                "  {x} and {y}: declared {} (expected relatedness {})\n",
                declared_relationship(x, y, &members, e),
                fixed(e, 3)
            ));
        }
    }
    for s in &seq {
        let ps: Vec<&str> = members[s]
            .parents()
            .filter(|p| seq.iter().any(|q| q == p))
            .collect();
        if !ps.is_empty() {
            say.push_str(&format!("  child {s}: {} sequenced\n", ps.join(" and ")));
        }
    }
    let tr: Vec<&str> = trios.iter().map(|s| s.as_str()).collect();
    say.push_str(&format!(
        "Trios WhatsHap phases: {}\n",
        if tr.is_empty() {
            "none (no child has both parents sequenced)".to_owned()
        } else {
            tr.join(", ")
        }
    ));
    print_out(&say);
    Ok(())
}
