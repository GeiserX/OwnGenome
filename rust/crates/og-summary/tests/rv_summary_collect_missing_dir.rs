//! `og summary collect` refuses a sample folder it cannot list: a missing
//! one, a file, or one without read permission exits 1 and leaves `--out`
//! as it was. An empty folder that exists is still collected, with every
//! section missing.

use og_summary::collect::cli;
use std::fs;
use std::path::PathBuf;

const GOOD: &str = "{\"a good summary\": true}\n";

struct Dir(PathBuf);

impl Dir {
    fn new(tag: &str) -> Dir {
        let p = std::env::temp_dir().join(format!(
            "og-summary-missing-dir-{tag}-{}",
            std::process::id()
        ));
        let _ = fs::remove_dir_all(&p);
        fs::create_dir_all(&p).unwrap();
        Dir(p)
    }

    fn path(&self, name: &str) -> String {
        self.0.join(name).display().to_string()
    }
}

impl Drop for Dir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn assets() -> String {
    concat!(env!("CARGO_MANIFEST_DIR"), "/../../../assets").into()
}

/// Runs the verb; gives the exit code and stderr.
fn run(sample_dir: &str, extra: &[&str]) -> (u8, String) {
    let mut argv: Vec<String> = [
        "--sample",
        "S",
        "--sample-dir",
        sample_dir,
        "--assets",
        &assets(),
    ]
    .iter()
    .map(|s| (*s).to_owned())
    .collect();
    argv.extend(extra.iter().map(|s| (*s).to_owned()));
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = cli(&argv, &mut out, &mut err);
    (code, String::from_utf8(err).unwrap())
}

/// The folder is refused: exit 1, the error names it, and `out` keeps the
/// good summary, byte for byte.
fn assert_refused(sample_dir: &str, out: &str, extra: &[&str]) {
    let mut args = vec!["--out", out];
    args.extend_from_slice(extra);
    let (code, err) = run(sample_dir, &args);
    assert_eq!(code, 1, "stderr: {err}");
    assert!(
        err.contains(&format!("ERROR: {sample_dir}: not a readable folder")),
        "stderr: {err}"
    );
    assert_eq!(fs::read_to_string(out).unwrap(), GOOD);
}

#[test]
fn a_missing_sample_folder_leaves_the_summary_as_it_was() {
    let d = Dir::new("missing");
    let out = d.path("summary.json");
    fs::write(&out, GOOD).unwrap();
    let typo = d.path("SAMPEL");
    assert_refused(&typo, &out, &[]);
    // Each group reads the same folder, so each refuses it.
    for g in ["a", "b", "c"] {
        assert_refused(&typo, &out, &["--only", g]);
    }
    assert!(!PathBuf::from(&typo).exists());
}

#[test]
fn without_out_a_missing_folder_writes_nothing() {
    let d = Dir::new("no-out");
    let typo = d.path("SAMPEL");
    let (code, err) = run(&typo, &[]);
    assert_eq!(code, 1, "stderr: {err}");
    // It stops on the folder, not on a write into a folder that is not there.
    assert!(err.contains("not a readable folder"), "stderr: {err}");
    assert!(!PathBuf::from(&typo).exists());
}

#[test]
fn a_file_is_not_a_sample_folder() {
    let d = Dir::new("file");
    let out = d.path("summary.json");
    fs::write(&out, GOOD).unwrap();
    let file = d.path("S.vcf.gz");
    fs::write(&file, "x").unwrap();
    assert_refused(&file, &out, &[]);
}

#[cfg(unix)]
#[test]
fn an_unreadable_sample_folder_leaves_the_summary_as_it_was() {
    use std::os::unix::fs::PermissionsExt;
    let d = Dir::new("unreadable");
    let out = d.path("summary.json");
    fs::write(&out, GOOD).unwrap();
    let locked = d.path("S");
    fs::create_dir_all(format!("{locked}/qc")).unwrap();
    fs::set_permissions(&locked, fs::Permissions::from_mode(0o000)).unwrap();
    // Root lists any folder: there is nothing to refuse.
    if fs::read_dir(&locked).is_err() {
        assert_refused(&locked, &out, &[]);
    }
    fs::set_permissions(&locked, fs::Permissions::from_mode(0o755)).unwrap();
}

#[test]
fn an_empty_sample_folder_is_still_collected() {
    let d = Dir::new("empty");
    let sample = d.path("S");
    fs::create_dir_all(&sample).unwrap();
    let out = d.path("summary.json");
    fs::write(&out, GOOD).unwrap();
    let (code, err) = run(&sample, &["--out", &out, "--only", "a"]);
    assert_eq!(code, 0, "stderr: {err}");
    let body = fs::read_to_string(&out).unwrap();
    assert_ne!(body, GOOD);
    assert!(body.contains("\"coverage\""), "{body}");
    assert!(body.contains("\"state\": \"missing\""), "{body}");
    assert!(!body.contains("\"state\": \"ok\""), "{body}");
}
