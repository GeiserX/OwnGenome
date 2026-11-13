//! og against the Python it replaced: each row of
//! `tests/data/parity/cases.tsv` runs `og pgx cyp2d6-*` or `og str motif` on
//! the inputs `make_parity.py` made, and each row of
//! `tests/data/reads/cases.tsv` runs `og cdna check|scan`, `og comphet
//! phase` or `og bam contigs` on the inputs `make_reads_parity.py` made.
//! Each must give the Python's exit code,
//! stdout and output files (`want/<case>/`), byte for byte. stderr is not
//! compared: it carries Python's exception texts.
//!
//! Where og differs from the Python on purpose, `og/<case>/` holds what og
//! must write instead, and the case fails when the Python's answer equals
//! it, so a fix that no longer differs cannot stay listed.

use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};

/// A scratch directory of this process: emptied when made, removed when
/// dropped, so a stale one or a panic leaves nothing behind.
struct Scratch(PathBuf);

impl Scratch {
    fn new(name: &str) -> Self {
        let dir = std::env::temp_dir().join(format!("{name}-{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).expect("scratch dir");
        Self(dir)
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn data() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data")
}

fn repo() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../..")
}

/// exit, stdout and every file of an answer directory, by name.
fn answer(dir: &Path) -> BTreeMap<String, Vec<u8>> {
    let mut out = BTreeMap::new();
    for item in ["exit", "stdout"] {
        out.insert(
            item.to_owned(),
            fs::read(dir.join(item)).unwrap_or_default(),
        );
    }
    walk(&dir.join("files"), &dir.join("files"), &mut out);
    out
}

fn walk(root: &Path, dir: &Path, out: &mut BTreeMap<String, Vec<u8>>) {
    let Ok(entries) = fs::read_dir(dir) else {
        return;
    };
    for e in entries {
        let p = e.unwrap().path();
        if p.is_dir() {
            walk(root, &p, out);
        } else {
            let rel = p.strip_prefix(root).unwrap().to_string_lossy().into_owned();
            out.insert(format!("files/{rel}"), fs::read(&p).unwrap());
        }
    }
}

fn run_og(tool: &str, args: &[String]) -> (u8, Vec<u8>) {
    let (mut stdout, mut stderr) = (Vec::new(), Vec::new());
    let rc = match tool {
        "str-motif" => og_bam::motif::cli(args, &mut stdout, &mut stderr),
        "cdna-check" => og_bam::cdna::check_cli(args, &mut stdout, &mut stderr),
        "cdna-scan" => og_bam::cdna::scan_cli(args, &mut stdout, &mut stderr),
        "comphet-phase" => og_bam::comphet::cli(args, &mut stdout, &mut stderr),
        "bam-contigs" => og_bam::contigs_cli(args, &mut stdout, &mut stderr),
        t => og_bam::cyp2d6::cli(
            t.strip_prefix("cyp2d6-").expect("a cyp2d6 verb"),
            args,
            &mut stdout,
            &mut stderr,
        ),
    };
    (rc, stdout)
}

fn show(b: &[u8]) -> String {
    String::from_utf8_lossy(b).into_owned()
}

#[test]
fn og_writes_what_the_python_wrote() {
    compare("parity", 60);
}

#[test]
fn og_reads_at_the_variant_as_the_python_did() {
    compare("reads", 40);
}

/// Runs every case of `tests/data/<set>/cases.tsv`; at least `min` must run
/// and at least one deliberate fix must be checked.
fn compare(set: &str, min: usize) {
    let parity = data().join(set);
    let cases = fs::read_to_string(parity.join("cases.tsv")).expect("cases.tsv");
    let scratch = Scratch::new(&format!("og-bam-{set}"));
    let (mut n, mut fixes, mut bad) = (0, 0, Vec::new());
    for line in cases
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
    {
        let mut f = line.split('\t');
        let (name, tool) = (f.next().unwrap(), f.next().unwrap());
        let out = scratch.0.join(name);
        fs::create_dir_all(&out).unwrap();
        let args: Vec<String> = f
            .map(|a| {
                a.replace("{d}", &data().to_string_lossy())
                    .replace("{r}", &repo().to_string_lossy())
                    .replace("{o}", &out.to_string_lossy())
            })
            .collect();
        let (rc, stdout) = run_og(tool, &args);
        let got_dir = scratch.0.join(format!("{name}.got"));
        fs::create_dir_all(&got_dir).unwrap();
        fs::write(got_dir.join("exit"), format!("{rc}\n")).unwrap();
        fs::write(got_dir.join("stdout"), &stdout).unwrap();
        fs::rename(&out, got_dir.join("files")).unwrap();
        let got = answer(&got_dir);
        let python = answer(&parity.join("want").join(name));
        let fix = parity.join("og").join(name);
        let want = if fix.is_dir() {
            fixes += 1;
            let w = answer(&fix);
            if w == python {
                bad.push(format!(
                    "{name}: og/ equals the Python's answer, so it fixes nothing"
                ));
            }
            w
        } else {
            python
        };
        n += 1;
        for key in want
            .keys()
            .chain(got.keys())
            .collect::<std::collections::BTreeSet<_>>()
        {
            match (want.get(key), got.get(key)) {
                (Some(w), Some(g)) if w == g => {}
                (w, g) => bad.push(format!(
                    "{name}: {key}: want {:?}, og wrote {:?}",
                    w.map(|b| show(b)),
                    g.map(|b| show(b))
                )),
            }
        }
    }
    assert!(n >= min, "only {n} cases ran");
    assert!(fixes >= 1, "no deliberate fix was checked");
    assert!(
        bad.is_empty(),
        "{} of {n} cases differ:\n{}",
        bad.len(),
        bad.join("\n")
    );
}
