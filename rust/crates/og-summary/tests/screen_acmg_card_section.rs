//! The ACMG secondary-findings card is its own section, `acmg_sf`, and
//! follows step 06: a run without the VEP cache (no step 23, no
//! `clinical/` folder) still gets it from the ClinVar hits.
//!
//! 1. Red first: a folder with step 06's hits (a BRCA2 P/LP heterozygote,
//!    two MUTYH heterozygotes) and no `clinical/` folder. Before the
//!    section existed it got no card at all.
//! 2. Agreement: every summary-collect parity folder of group b
//!    (`tests/parity/gen/summary-clinical.sh`) that has both the hits and
//!    the clinical VCF, collected as it is and again without `clinical/`.
//!    Without step 23 the card may lose only what the HIGH-impact source
//!    brings: (a) every ClinVar row of the card without VEP is in the card
//!    with VEP, unchanged; (b) every gene reported without VEP is reported
//!    with VEP; (c) a row not reported without VEP stays as it is, or is
//!    reported with VEP because a HIGH-impact row of its gene joined it.

use og_summary::collect::{Group, Input, group_b, post};
use pycompat::json::Value;
use std::path::{Path, PathBuf};
use std::process::Command;

fn repo() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../..")
}

fn scratch(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!(
        "og-screen-acmg-{name}-{}",
        std::process::id()
    ));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn put(p: &Path, text: &str) {
    std::fs::create_dir_all(p.parent().unwrap()).unwrap();
    std::fs::write(p, text).unwrap();
}

fn input(dir: &Path, assets: &Path, group: Group) -> Input {
    Input {
        sample: "S".into(),
        sample_dir: dir.to_str().unwrap().to_owned(),
        declared_sex: None,
        assets: assets.to_str().unwrap().to_owned(),
        group,
        tables: og_summary::collect::group_c::util::default_tables(),
    }
}

/// The sections of group b, as summary.json writes them.
fn sections_b(dir: &Path, assets: &Path) -> Value {
    let v = group_b::collect(&input(dir, assets, Group::B))
        .unwrap()
        .to_json();
    at(&v, &["sections"]).clone()
}

fn at<'a>(v: &'a Value, path: &[&str]) -> &'a Value {
    path.iter().fold(v, |v, k| match v {
        Value::Object(m) => m
            .get(*k)
            .unwrap_or_else(|| panic!("no {k} in {path:?}: {v:?}")),
        _ => panic!("{path:?} is no object at {k}"),
    })
}

fn get<'a>(v: &'a Value, k: &str) -> Option<&'a Value> {
    match v {
        Value::Object(m) => m.get(k),
        _ => None,
    }
}

fn s(v: &Value) -> &str {
    match v {
        Value::Str(x) => x,
        _ => panic!("no string: {v:?}"),
    }
}

fn rows(card: &Value, list: &str) -> Vec<Value> {
    match get(card, list) {
        Some(Value::Array(a)) => a.clone(),
        None => Vec::new(),
        Some(x) => panic!("{list} is no list: {x:?}"),
    }
}

fn gene_variant(r: &Value) -> (String, String) {
    (
        s(at(r, &["gene"])).to_owned(),
        s(at(r, &["variant"])).to_owned(),
    )
}

const HITS_HEAD: &str = "##fileformat=VCFv4.2\n##PGP_CALL_FLOOR=GQ=20,DP=8\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n";
const PL: &str = "CLNREVSTAT=criteria_provided,_multiple_submitters,_no_conflicts";

/// Step 06's hits: a BRCA2 P/LP heterozygote and two MUTYH heterozygotes.
fn hits() -> String {
    format!(
        "{HITS_HEAD}\
chr13\t32340300\t1\tC\tT\t.\tPASS\tGENEINFO=BRCA2:675;CLNSIG=Pathogenic;{PL}\tGT:GQ:DP\t0/1:60:35\n\
chr1\t45331556\t2\tC\tT\t.\tPASS\tGENEINFO=MUTYH:4595;CLNSIG=Pathogenic;{PL}\tGT:GQ:DP\t0/1:50:30\n\
chr1\t45332400\t3\tG\tA\t.\tPASS\tGENEINFO=MUTYH:4595;CLNSIG=Likely_pathogenic;{PL}\tGT:GQ:DP\t0/1:55:32\n"
    )
}

#[test]
fn a_run_without_step_23_gets_the_card_from_the_clinvar_hits() {
    let d = scratch("red");
    let sd = d.join("S");
    put(&sd.join("clinvar/S_clinvar_hits.vcf"), &hits());
    let assets = repo().join("assets");
    let secs = sections_b(&sd, &assets);
    assert!(
        get(&secs, "clinical").is_some_and(|c| s(at(c, &["state"])) == "missing"),
        "no clinical section: step 23 did not run"
    );
    let sec = get(&secs, "acmg_sf").expect("summary.json has no acmg_sf section");
    assert_eq!(s(at(sec, &["state"])), "ok", "{sec:?}");
    assert_eq!(s(at(sec, &["title"])), "ACMG secondary findings");
    assert_eq!(
        s(at(sec, &["source"])),
        "clinvar/S_clinvar_hits.vcf",
        "the section's file is step 06's"
    );
    let card = at(sec, &["data"]);
    assert_eq!(at(card, &["genes_on_list"]), &Value::Int(84));
    let reported: Vec<(String, String)> = rows(card, "clinvar_hits")
        .iter()
        .map(|r| {
            let rule = s(at(r, &["rule"])).to_owned();
            (s(at(r, &["gene"])).to_owned(), rule)
        })
        .collect();
    assert_eq!(
        reported,
        [
            ("BRCA2".to_owned(), "All P and LP".to_owned()),
            ("MUTYH".to_owned(), "P and LP (2 variants)".to_owned()),
            ("MUTYH".to_owned(), "P and LP (2 variants)".to_owned()),
        ],
        "BRCA2 and the MUTYH pair are reported under their rules"
    );
    assert!(rows(card, "not_reported").is_empty(), "{card:?}");
    assert_eq!(at(card, &["clinvar_checked"]), &Value::Bool(true));
    assert_eq!(at(card, &["high_impact_checked"]), &Value::Bool(false));
    let note = s(at(card, &["high_impact_note"]));
    assert!(
        note.starts_with("Rare HIGH-impact variants were not screened: step 23 (clinical filter) did not run")
            && note.contains("not in ClinVar"),
        "{note}"
    );

    // The run status gives step 23's reason; a failed step 23 is named.
    let skip = "skipped (needs VEP, data not installed: vep_cache/homo_sapiens/116_GRCh38/info.txt)";
    for (word, want) in [
        (skip, format!("step 23 (clinical filter) was {skip}")),
        ("failed (exit 1)", "step 23 (clinical filter) failed (exit 1)".to_owned()),
    ] {
        put(
            &sd.join("logs/run_status.tsv"),
            &format!("meta\tstarted_utc\t2020-01-01 00:00 UTC\nstep\t06\tok\nstep\t23\t{word}\n"),
        );
        let secs = sections_b(&sd, &assets);
        let note = s(at(&secs, &["acmg_sf", "data", "high_impact_note"])).to_owned();
        assert!(note.contains(&want), "{note}");
    }
}

#[test]
fn the_section_follows_step_06_and_is_never_a_not_run_line_of_its_own() {
    let d = scratch("state");
    let sd = d.join("S");
    std::fs::create_dir_all(&sd).unwrap();
    let assets = repo().join("assets");
    let secs = sections_b(&sd, &assets);
    assert_eq!(s(at(&secs, &["acmg_sf", "state"])), "missing");
    put(
        &sd.join("logs/run_status.tsv"),
        "meta\tstarted_utc\t2020-01-01 00:00 UTC\nstep\t06\tfailed (exit 2)\n",
    );
    let secs = sections_b(&sd, &assets);
    assert_eq!(s(at(&secs, &["acmg_sf", "state"])), "failed");
    assert!(s(at(&secs, &["acmg_sf", "note"])).starts_with("step 06 failed"));
    // The whole summary: ClinVar screen names the failed step once.
    let all = post::collect(&input(&sd, &assets, Group::All)).unwrap();
    let not_run: Vec<&str> = match at(&all, &["not_run"]) {
        Value::Array(a) => a.iter().map(s).collect(),
        x => panic!("{x:?}"),
    };
    assert!(not_run.contains(&"ClinVar screen (step 06 failed)"), "{not_run:?}");
    assert!(
        !not_run.iter().any(|l| l.starts_with("ACMG secondary findings")),
        "{not_run:?}"
    );
}

/// The parity folders of group b, written by their generator.
fn parity_folders(d: &Path) -> PathBuf {
    let script = repo().join("tests/parity/gen/summary-clinical.sh");
    let st = Command::new("bash")
        .arg(&script)
        .current_dir(d)
        .env("NO_ORACLE", "1")
        .status()
        .expect("bash runs");
    assert!(st.success(), "{} failed", script.display());
    for f in ["cb_vep.vcf.gz", "cb_vep.vcf.gz.tbi"] {
        std::fs::copy(repo().join("tests/parity/inputs/summary-collect").join(f), d.join(f)).unwrap();
    }
    d.join("assets")
}

/// Copies `from` to `to`, all but its `clinical/` folder (symlinks kept).
fn copy_without_clinical(from: &Path, to: &Path) {
    let st = Command::new("cp")
        .arg("-PR")
        .arg(from)
        .arg(to)
        .status()
        .unwrap();
    assert!(st.success());
    std::fs::remove_dir_all(to.join("clinical")).unwrap();
}

/// The three rules of the module doc, for one folder: the reasons it
/// breaks them, none when it keeps them.
fn agreement(with: &Value, without: &Value) -> Vec<String> {
    let mut bad = Vec::new();
    let high_genes: Vec<String> = rows(with, "high_impact")
        .iter()
        .map(|r| gene_variant(r).0)
        .collect();
    let with_cv = rows(with, "clinvar_hits");
    let without_cv = rows(without, "clinvar_hits");
    // (a) every ClinVar row without VEP, unchanged with VEP.
    for r in &without_cv {
        if !with_cv.contains(r) {
            bad.push(format!("(a) ClinVar row without VEP not in the card with VEP: {r:?}"));
        }
    }
    // (a') a ClinVar row only with VEP is one (c) explains.
    let reported_with: Vec<(String, String)> = ["clinvar_hits", "high_impact", "structural"]
        .iter()
        .flat_map(|l| rows(with, l))
        .map(|r| gene_variant(&r))
        .collect();
    let not_without: Vec<(String, String)> =
        rows(without, "not_reported").iter().map(gene_variant).collect();
    for r in &with_cv {
        if !without_cv.contains(r) {
            let gv = gene_variant(r);
            let moved = not_without.contains(&gv) && high_genes.contains(&gv.0);
            // A row whose fields grew (a two-hit line) is the same row.
            let same = without_cv.iter().any(|w| gene_variant(w) == gv);
            if !(moved || same) {
                bad.push(format!("(a) ClinVar row only with VEP, not explained: {r:?}"));
            }
        }
    }
    // (b) every gene reported without VEP is reported with VEP.
    for (g, _) in ["clinvar_hits", "structural"]
        .iter()
        .flat_map(|l| rows(without, l))
        .map(|r| gene_variant(&r))
    {
        if !reported_with.iter().any(|(h, _)| *h == g) {
            bad.push(format!("(b) {g} reported without VEP, not with VEP"));
        }
    }
    // (c) a not-reported row stays, or is reported because a HIGH-impact
    // row of its gene joined it.
    let not_with = rows(with, "not_reported");
    for r in rows(without, "not_reported") {
        let gv = gene_variant(&r);
        let stays = not_with.contains(&r);
        let joined = reported_with.contains(&gv) && high_genes.contains(&gv.0);
        if !(stays || joined) {
            bad.push(format!("(c) not reported without VEP, neither kept nor joined: {r:?}"));
        }
    }
    bad
}

#[test]
fn without_the_clinical_vcf_the_card_loses_only_what_the_high_impact_source_brings() {
    let d = scratch("agree");
    let assets = parity_folders(&d);
    let mut compared = Vec::new();
    let mut moved = 0;
    let mut names: Vec<_> = std::fs::read_dir(&d)
        .unwrap()
        .map(|e| e.unwrap().path())
        .filter(|p| p.join("clinvar/S_clinvar_hits.vcf").is_file() && p.join("clinical").is_dir())
        .collect();
    names.sort();
    for dir in names {
        let name = dir.file_name().unwrap().to_str().unwrap().to_owned();
        let bare = d.join(format!("{name}.noclinical"));
        copy_without_clinical(&dir, &bare);
        let with = sections_b(&dir, &assets);
        let without = sections_b(&bare, &assets);
        let wsec = at(&with, &["acmg_sf"]);
        let nsec = at(&without, &["acmg_sf"]);
        if s(at(wsec, &["state"])) != "ok" {
            // A hits VCF that cannot be read is the same note either way.
            assert_eq!(wsec, nsec, "{name}");
            continue;
        }
        let (w, n) = (at(wsec, &["data"]), at(nsec, &["data"]));
        assert_eq!(at(n, &["high_impact_checked"]), &Value::Bool(false), "{name}");
        assert!(
            s(at(n, &["high_impact_note"])).contains("did not run")
                || s(at(n, &["high_impact_note"])).contains("wrote no clinical VCF"),
            "{name}: {n:?}"
        );
        if at(w, &["high_impact_checked"]) == &Value::Bool(false) {
            // The clinical VCF is there and cannot be read: the card keeps
            // its ClinVar rows and says why.
            assert!(
                s(at(w, &["high_impact_note"])).contains("clinical VCF cannot be read"),
                "{name}: {w:?}"
            );
            assert_eq!(rows(w, "clinvar_hits"), rows(n, "clinvar_hits"), "{name}");
            continue;
        }
        assert!(get(w, "high_impact_note").is_none(), "{name}: {w:?}");
        // The clinical section's own tier, with the same inputs, agrees.
        assert_eq!(
            Some(w),
            get(at(&with, &["clinical", "data"]), "acmg_sf"),
            "{name}: the card and the clinical section's tier differ"
        );
        let bad = agreement(w, n);
        assert!(bad.is_empty(), "{name}:\n{}", bad.join("\n"));
        let not_with = rows(w, "not_reported");
        let n_moved = rows(n, "not_reported")
            .iter()
            .filter(|r| !not_with.contains(r))
            .count();
        eprintln!(
            "{name}: ClinVar rows {} without VEP, {} with; not reported {} without, {} with; {n_moved} joined",
            rows(n, "clinvar_hits").len(),
            rows(w, "clinvar_hits").len(),
            rows(n, "not_reported").len(),
            not_with.len()
        );
        moved += n_moved;
        compared.push(name);
    }
    let _ = std::fs::remove_dir_all(&d);
    assert!(
        compared.len() >= 5,
        "only {compared:?} compared: the generator lost its folders"
    );
    assert!(
        moved > 0,
        "no row moved from not reported to reported with VEP: the folders no longer test rule (c)"
    );
}
