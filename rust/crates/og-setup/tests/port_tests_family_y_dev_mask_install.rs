//! `setup.sh --yleaf-data` keeps step 37's chrY callable mask between
//! runs. It replaces `tests/test_y_private_mask_install.py`.
//!
//! The mask `og y private-snps` keeps to lives in `<yleaf data>/ytree`
//! beside the current tree, checked against `YCALLABLE_MASK_SHA256`. This
//! test runs the real `scripts/setup.sh` and `scripts/lib/common.sh` from a
//! copy of the repository whose `versions.env` pins synthetic tree, table
//! and mask files, with a fake curl that serves them and logs every URL it
//! is asked for (no network). Yleaf's own tables are already in place, so
//! only the tree folder is in play:
//!
//! 1. a first run downloads the three tree files and the mask;
//! 2. a second run asks for nothing;
//! 3. a mask that fails its sha256 (one byte added) is downloaded again,
//!    and nothing else is;
//! 4. a new tree pin (an older version.txt) downloads the three tree files
//!    and keeps the mask, with no request for it;
//! 5. a tree download that fails keeps the checked mask, and no
//!    version.txt, so no step reads the tree; the next run asks for the
//!    tree files only;
//! 6. a missing mask beside current tables is downloaded on its own.

use flate2::Compression;
use flate2::write::GzEncoder;
use std::fs;
use std::io::Write;
use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};
use std::process::Command;

const FAKE_CURL: &str = r#"#!/usr/bin/env bash
# Serves $SERVE/<URL's base name> to -o FILE (or stdout) and logs the URL to
# $CURL_LOG. A URL matching $CURL_FAIL, or a file not in $SERVE, exits 22.
out="" url=""
while [ $# -gt 0 ]; do
  case "$1" in
    -o) out=$2; shift 2 ;;
    -C|--connect-timeout|--speed-limit|--speed-time|--retry|--max-time) shift 2 ;;
    -*) shift ;;
    *) url=$1; shift ;;
  esac
done
echo "$url" >> "$CURL_LOG"
if [ -n "${CURL_FAIL:-}" ] && [[ "$url" =~ $CURL_FAIL ]]; then exit 22; fi
src="$SERVE/$(basename "$url")"
[ -f "$src" ] || exit 22
if [ -z "$out" ] || [ "$out" = "-" ]; then cat "$src"; else cat "$src" > "$out"; fi
"#;

fn repo() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../..")
}

/// The file's sha256 (through `sha256sum`, else `shasum -a 256`, as
/// `scripts/lib/common.sh` computes it), or None when it is not there.
fn sha256(path: &Path) -> Option<String> {
    if !path.is_file() {
        return None;
    }
    let out = Command::new("sha256sum")
        .arg(path)
        .output()
        .or_else(|_| {
            Command::new("shasum")
                .args(["-a", "256"])
                .arg(path)
                .output()
        })
        .unwrap();
    assert!(out.status.success(), "sha256 of {}", path.display());
    let text = String::from_utf8(out.stdout).unwrap();
    Some(text.split_whitespace().next().unwrap().to_owned())
}

/// The value of `NAME="..."` in versions.env.
fn pin(env: &str, name: &str) -> String {
    let prefix = format!("{name}=\"");
    let line = env
        .lines()
        .find(|l| l.starts_with(&prefix))
        .unwrap_or_else(|| panic!("versions.env has {name}"));
    line[prefix.len()..].split('"').next().unwrap().to_owned()
}

/// versions.env with `NAME="..."` set to `value`; the rest of the line
/// (a comment) stays. Panics unless exactly one line matches.
fn set_pin(env: &str, name: &str, value: &str) -> String {
    let prefix = format!("{name}=\"");
    let mut hits = 0;
    let out: String = env
        .split_inclusive('\n')
        .map(
            |l| match l.strip_prefix(&prefix).and_then(|r| r.split_once('"')) {
                Some((_, rest)) => {
                    hits += 1;
                    format!("{prefix}{value}\"{rest}")
                }
                None => l.to_owned(),
            },
        )
        .collect();
    assert_eq!(hits, 1, "versions.env has one {name} line");
    out
}

struct World {
    work: PathBuf,
    root: PathBuf,
    serve: PathBuf,
    genome: PathBuf,
    bindir: PathBuf,
    log: PathBuf,
}

impl World {
    /// Runs setup.sh --yleaf-data; returns (exit code, the base names of
    /// the URLs asked for, the output).
    fn setup(&self, fail: &str) -> (i32, Vec<String>, String) {
        fs::write(&self.log, "").unwrap();
        let path = format!(
            "{}:{}",
            self.bindir.display(),
            std::env::var("PATH").unwrap_or_default()
        );
        let p = Command::new("bash")
            .arg(self.root.join("scripts/setup.sh"))
            .arg("--yleaf-data")
            .arg(&self.genome)
            .env("PATH", path)
            .env("SERVE", &self.serve)
            .env("CURL_LOG", &self.log)
            .env("CURL_FAIL", fail)
            .env("CONTAINER_ENGINE", "true")
            .env("FETCH_TRIES", "1")
            .env("FETCH_WAIT", "0")
            .output()
            .unwrap();
        let asked = fs::read_to_string(&self.log)
            .unwrap()
            .lines()
            .filter(|l| !l.trim().is_empty())
            .map(|l| l.trim().rsplit('/').next().unwrap().to_owned())
            .collect();
        let mut out = String::from_utf8_lossy(&p.stdout).into_owned();
        out += &String::from_utf8_lossy(&p.stderr);
        let tail = out[out.len().saturating_sub(800)..].to_owned();
        (p.status.code().unwrap_or(-1), asked, tail)
    }
}

#[test]
fn setup_keeps_the_checked_chry_mask() {
    let work = std::env::temp_dir().join(format!("og-ymask-{}", std::process::id()));
    let _ = fs::remove_dir_all(&work);
    let w = World {
        root: work.join("repo"),
        serve: work.join("serve"),
        genome: work.join("genome"),
        bindir: work.join("bin"),
        log: work.join("curl.log"),
        work,
    };
    fs::create_dir_all(w.root.join("scripts/lib")).unwrap();
    fs::create_dir_all(&w.serve).unwrap();
    fs::create_dir_all(&w.bindir).unwrap();
    fs::copy(
        repo().join("scripts/setup.sh"),
        w.root.join("scripts/setup.sh"),
    )
    .unwrap();
    fs::copy(
        repo().join("scripts/lib/common.sh"),
        w.root.join("scripts/lib/common.sh"),
    )
    .unwrap();
    let mut env = fs::read_to_string(repo().join("versions.env")).unwrap();
    let release = pin(&env, "YCALLABLE_MASK_RELEASE");
    let yleaf = pin(&env, "YLEAF_DATA_VERSION");
    let mask_name = format!("{release}.chrY.mask.fasta.gz");
    let mut gz = GzEncoder::new(Vec::new(), Compression::default());
    gz.write_all(b">chrY\nNNPPPPNN\n").unwrap();
    // The synthetic files, under the base names of the URLs setup.sh asks for.
    let files: [(&str, Vec<u8>); 4] = [
        ("current_tree.json", b"{\"id\": \"synthetic\"}\n".to_vec()),
        (
            "ftdna_positions_hg38.txt",
            b"A1\tchrY\t100\tA\tG\n".to_vec(),
        ),
        ("new_positions.txt", b"B1\tchrY\t200\tC\tT\n".to_vec()),
        (&mask_name, gz.finish().unwrap()),
    ];
    for (name, data) in &files {
        fs::write(w.serve.join(name), data).unwrap();
    }
    for (var, name) in [
        ("YTREE_SHA256", "current_tree.json"),
        ("YTREE_FTDNA_SHA256", "ftdna_positions_hg38.txt"),
        ("YTREE_YFULL_SHA256", "new_positions.txt"),
        ("YCALLABLE_MASK_SHA256", mask_name.as_str()),
    ] {
        env = set_pin(&env, var, &sha256(&w.serve.join(name)).unwrap());
    }
    fs::write(w.root.join("versions.env"), env).unwrap();
    let data = w.genome.join(format!("reference/yleaf-{yleaf}/data"));
    for rel in ["hg38/new_positions.txt", "hg_prediction_tables/tree.json"] {
        let p = data.join(rel);
        fs::create_dir_all(p.parent().unwrap()).unwrap();
        fs::write(p, "seeded\n").unwrap();
    }
    let curl = w.bindir.join("curl");
    fs::write(&curl, FAKE_CURL).unwrap();
    fs::set_permissions(&curl, fs::Permissions::from_mode(0o755)).unwrap();
    let yt = data.join("ytree");
    let mask = yt.join("callable_mask_chrY.fasta.gz");
    let version = yt.join("version.txt");
    let want_mask = sha256(&w.serve.join(&mask_name));
    let tree3 = [
        "current_tree.json",
        "ftdna_positions_hg38.txt",
        "new_positions.txt",
    ];
    let with_mask: Vec<&str> = tree3.iter().copied().chain([mask_name.as_str()]).collect();

    let (rc, asked, out) = w.setup("");
    assert!(
        rc == 0 && asked == with_mask && version.is_file() && sha256(&mask) == want_mask,
        "1. first run: the three tree files and the mask are downloaded: {rc} {asked:?} {out}"
    );
    let (rc, asked, out) = w.setup("");
    assert!(
        rc == 0 && asked.is_empty(),
        "2. second run: no request at all: {rc} {asked:?} {out}"
    );
    let mut f = fs::OpenOptions::new().append(true).open(&mask).unwrap();
    f.write_all(b"x").unwrap();
    drop(f);
    let version_before = fs::read_to_string(&version).unwrap();
    let (rc, asked, out) = w.setup("");
    assert!(
        rc == 0
            && asked == [mask_name.as_str()]
            && sha256(&mask) == want_mask
            && fs::read_to_string(&version).unwrap() == version_before,
        "3. a mask that fails its sha256 is downloaded again, and nothing else: {rc} {asked:?} {out}"
    );
    fs::write(&version, "tree\tan older pin\n").unwrap();
    let (rc, asked, out) = w.setup("");
    assert!(
        rc == 0
            && asked == tree3
            && sha256(&mask) == want_mask
            && fs::read_to_string(&version).unwrap() == version_before,
        "4. a new tree pin downloads the tree files and keeps the mask: {rc} {asked:?} {out}"
    );
    fs::write(&version, "tree\tan older pin\n").unwrap();
    let (rc, asked, out) = w.setup("current_tree");
    let left: Vec<String> = fs::read_dir(&yt)
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    assert!(
        rc == 0
            && !asked.contains(&mask_name)
            && sha256(&mask) == want_mask
            && !version.exists()
            && left == ["callable_mask_chrY.fasta.gz"],
        "5. a failed tree download keeps the checked mask and leaves no version.txt: \
         {rc} {asked:?} {left:?} {out}"
    );
    let (rc, asked, out) = w.setup("");
    assert!(
        rc == 0 && asked == tree3 && sha256(&mask) == want_mask,
        "5. the next run asks for the tree files only: {rc} {asked:?} {out}"
    );
    fs::remove_file(&mask).unwrap();
    let (rc, asked, out) = w.setup("");
    assert!(
        rc == 0 && asked == [mask_name.as_str()] && sha256(&mask) == want_mask,
        "6. a missing mask beside current tables is downloaded on its own: {rc} {asked:?} {out}"
    );
    fs::remove_dir_all(&w.work).unwrap();
}
