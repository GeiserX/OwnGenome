//! `og prs format` leaves out the scores `assets/pgs_scores.tsv` marks
//! `academic-research` unless it gets `--academic`, whichever way the scores
//! are chosen (the folder, or `--ids`).

use og_steps::prs::{ACADEMIC_REASON, Args, Labels, run};
use std::fs;
use std::path::{Path, PathBuf};

/// The four scores the PGS Catalog publishes for academic research use only
/// (its REST `license` field; `scripts/ci/check-pgs-labels.sh` checks it),
/// in the list's order.
const ACADEMIC: [&str; 4] = ["PGS000014", "PGS000016", "PGS000027", "PGS000017"];

fn list() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../../assets/pgs_scores.tsv")
}

fn ids() -> Vec<String> {
    fs::read_to_string(list())
        .expect("assets/pgs_scores.tsv")
        .lines()
        .filter_map(|l| l.split('\t').next())
        .filter(|id| id.starts_with("PGS") && id[3..].chars().all(|c| c.is_ascii_digit()))
        .map(str::to_owned)
        .collect()
}

/// A folder under the system temp folder with one small GRCh38 score per
/// id, removed on drop.
struct Case(PathBuf);

impl Case {
    fn new(tag: &str, ids: &[String]) -> Case {
        let d = std::env::temp_dir().join(format!("og-prs-academic-{tag}-{}", std::process::id()));
        let _ = fs::remove_dir_all(&d);
        fs::create_dir_all(d.join("in")).expect("dir");
        for id in ids {
            fs::write(
                d.join(format!("in/{id}_hmPOS_GRCh38.txt")),
                format!("#pgs_id={id}\n#HmPOS_build=GRCh38\nhm_chr\thm_pos\teffect_allele\teffect_weight\n1\t5\tA\t1\n"),
            )
            .expect("score");
        }
        Case(d)
    }

    fn run(&self, extra: &[&str]) -> Result<(String, Vec<String>), String> {
        self.run_with_labels(Some(&list()), extra)
    }

    fn run_with_labels(
        &self,
        labels: Option<&Path>,
        extra: &[&str],
    ) -> Result<(String, Vec<String>), String> {
        let mut v: Vec<String> = vec![
            "--scores".into(),
            self.0.join("in").display().to_string(),
            "--out".into(),
            self.0.join("pgs").display().to_string(),
            "--alleles".into(),
            self.0.join("alleles.tsv").display().to_string(),
        ];
        if let Some(labels) = labels {
            v.extend(["--labels".to_owned(), labels.display().to_string()]);
        }
        v.extend(extra.iter().map(|s| (*s).to_owned()));
        let _ = fs::remove_dir_all(self.0.join("pgs"));
        let mut out = Vec::new();
        run(&Args::parse(&v)?, &mut out)?;
        let mut written: Vec<String> = fs::read_dir(self.0.join("pgs"))
            .expect("pgs")
            .flatten()
            .map(|e| e.file_name().to_string_lossy().into_owned())
            .filter_map(|n| n.strip_suffix(".txt.gz").map(str::to_owned))
            .collect();
        written.sort();
        Ok((String::from_utf8(out).expect("utf8"), written))
    }
}

impl Drop for Case {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

#[test]
fn the_list_marks_exactly_the_four_academic_scores() {
    let l = Labels::read(&list()).expect("read the list");
    assert_eq!(l.academic_ids(), ACADEMIC);
    assert_eq!(ids().len(), 15, "the list has 15 scores");
}

#[test]
fn a_folder_with_all_15_scores_11_by_default_and_15_with_academic() {
    let all = ids();
    let c = Case::new("folder", &all);
    let (said, written) = c.run(&[]).expect("default");
    assert_eq!(written.len(), 11, "{written:?}");
    for id in ACADEMIC {
        assert!(
            !written.contains(&id.to_owned()),
            "{id} was scored by default"
        );
        assert!(
            said.contains(&format!("  {id}: not scored, {ACADEMIC_REASON} (")),
            "{id} not named on stdout:\n{said}"
        );
    }
    assert!(said.contains(
        "  PGS000014: not scored, academic-only licence (setup.sh --prs-academic, then PRS_ACADEMIC=true, to include) (Type 2 diabetes (T2D))\n"
    ));
    // The scores' own files never name the four.
    let notes = fs::read_to_string(c.0.join("pgs/score_notes.tsv")).expect("notes");
    assert!(ACADEMIC.iter().all(|id| !notes.contains(id)), "{notes}");

    let (said, written) = c.run(&["--academic"]).expect("--academic");
    assert_eq!(written.len(), 15, "{written:?}");
    assert!(!said.contains("not scored"), "{said}");
}

#[test]
fn ids_naming_the_four_skip_them_even_when_their_files_are_gone() {
    // Step 25 passes every id of the list; a new setup no longer downloads
    // the four, so their files are missing, and that is not an error.
    let all = ids();
    let kept: Vec<String> = all
        .iter()
        .filter(|id| !ACADEMIC.contains(&id.as_str()))
        .cloned()
        .collect();
    let c = Case::new("ids", &kept);
    let joined = all.join(",");
    let (said, written) = c.run(&["--ids", &joined]).expect("default");
    assert_eq!(written, {
        let mut k = kept.clone();
        k.sort();
        k
    });
    assert_eq!(said.matches(": not scored, ").count(), 4, "{said}");
    // With --academic the missing files are an error again.
    let err = c
        .run(&["--ids", &joined, "--academic"])
        .expect_err("missing");
    assert!(
        err.starts_with("no scoring file for PGS000014, PGS000016, PGS000027, PGS000017 in "),
        "{err}"
    );
}

#[test]
fn only_academic_scores_asked_for_is_an_error_that_names_them() {
    let mut four: Vec<String> = ACADEMIC.iter().map(|s| (*s).to_owned()).collect();
    four.sort();
    let c = Case::new("only", &four);
    let err = c.run(&[]).expect_err("nothing to score");
    assert_eq!(
        err,
        "every score asked for is free for academic research use only (PGS000014, PGS000016, PGS000017, PGS000027); --academic scores them"
    );
    assert_eq!(c.run(&["--academic"]).expect("--academic").1, four);
}

#[test]
fn shipped_academic_scores_stay_opt_in_without_licence_labels() {
    let all = ids();
    let c = Case::new("no-licence-labels", &all);
    let legacy_list = c.0.join("legacy.tsv");
    fs::write(
        &legacy_list,
        "pgs_id\ttrait_reported\nPGS000014\tDiabetes\n",
    )
    .expect("legacy labels");
    let mut expected: Vec<_> = all
        .iter()
        .filter(|id| !ACADEMIC.contains(&id.as_str()))
        .cloned()
        .collect();
    expected.sort();

    for labels in [None, Some(legacy_list.as_path())] {
        let (said, written) = c.run_with_labels(labels, &[]).expect("default");
        assert_eq!(written, expected, "without licence labels");
        assert_eq!(said.matches(": not scored, ").count(), 4, "{said}");
        let (_, written) = c.run_with_labels(labels, &["--academic"]).expect("opt-in");
        assert_eq!(written.len(), 15);
    }
}

#[test]
fn a_file_whose_header_names_an_academic_score_is_refused_by_default() {
    // The written file takes its id from #pgs_id, not from the file name, so
    // the gate must hold on that id too, before anything is written.
    let c = Case::new("renamed", &["PGS000018".to_owned()]);
    fs::write(
        c.0.join("in/PGS900014_hmPOS_GRCh38.txt"),
        "#pgs_id=PGS000014\n#HmPOS_build=GRCh38\nhm_chr\thm_pos\teffect_allele\teffect_weight\n1\t5\tA\t1\n",
    )
    .expect("score");
    for extra in [&[][..], &["--ids", "PGS000018,PGS900014"][..]] {
        let err = c.run(extra).expect_err("refused");
        assert_eq!(
            err,
            "PGS900014_hmPOS_GRCh38.txt is PGS000014 by its #pgs_id line, a score free for academic research use only; --academic scores it"
        );
        assert!(!c.0.join("pgs/PGS000014.txt.gz").exists());
    }
    let (_, written) = c.run(&["--academic"]).expect("--academic");
    assert_eq!(written, ["PGS000014", "PGS000018"]);
}
