//! An index older than its VCF is not trusted: a contig the VCF gained after
//! the index was built is not in the index, so an indexed read would answer
//! "no call" for a variant the file holds. The lookup reads the whole file
//! instead and says so; a fresh index is still used.

use std::fs::{self, File};
use std::io::Write;
use std::num::NonZero;
use std::path::{Path, PathBuf};
use std::time::{Duration, SystemTime};

fn input(name: &str) -> String {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../../tests/parity/inputs/lookup")
        .join(name)
        .to_string_lossy()
        .into_owned()
}

/// Runs `og lookup query chrM:73` over `gvcf`: stdout, stderr.
fn chrm_73(gvcf: &str) -> (String, String) {
    let args: Vec<String> = [
        "query",
        "--sites",
        &input("sites.tsv.gz"),
        "--genes",
        &input("gencode.v50.basic.genes.gtf"),
        "--table",
        &input("traits.tsv"),
        "--gvcf",
        gvcf,
        "--",
        "chrM:73",
    ]
    .map(String::from)
    .to_vec();
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let rc = og_lookup::cli(&args, &mut out, &mut err);
    let (out, err) = (
        String::from_utf8(out).unwrap(),
        String::from_utf8(err).unwrap(),
    );
    assert_eq!(rc, 0, "{out}{err}");
    (out, err)
}

/// Writes `text` to `path` as BGZF.
fn bgzip(path: &Path, text: &str) {
    let mut w =
        og_io::bgzf::Writer::new(File::create(path).unwrap(), NonZero::<usize>::MIN).unwrap();
    w.write_all(text.as_bytes()).unwrap();
    w.finish().unwrap();
}

/// Builds the `.tbi` of the BGZF VCF `vcf`.
fn index(vcf: &Path, tbi: &Path) {
    let built = og_io::tabix_write::build(File::open(vcf).unwrap()).unwrap();
    built.write_tbi(tbi, true).unwrap();
}

/// Sets the modification time of `path` a minute after (or before) that of `of`.
fn touch(path: &Path, of: &Path, later: bool) {
    let base = fs::metadata(of).unwrap().modified().unwrap();
    let minute = Duration::from_secs(60);
    let t: SystemTime = if later { base + minute } else { base - minute };
    File::options()
        .write(true)
        .open(path)
        .unwrap()
        .set_modified(t)
        .unwrap();
}

#[test]
fn an_index_older_than_its_vcf_is_not_trusted() {
    let dir: PathBuf = Path::new(env!("CARGO_TARGET_TMPDIR")).join("rv_lookup_stale_index");
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    let (gz, tbi) = (dir.join("g.vcf.gz"), dir.join("g.vcf.gz.tbi"));
    let gz_s = gz.to_string_lossy().into_owned();
    let full = fs::read_to_string(input("g.vcf")).unwrap();
    assert!(full.contains("\nchrM\t73\t"), "the fixture calls chrM:73");

    // The answer from the whole plain file.
    let (want, _) = chrm_73(&input("g.vcf"));
    assert!(want.contains("genotype: G/G"), "{want}");

    // The index of an earlier VCF without chrM, then the VCF rewritten with it.
    let old: String = full
        .split_inclusive('\n')
        .filter(|l| !l.starts_with("chrM\t"))
        .collect();
    bgzip(&gz, &old);
    index(&gz, &tbi);
    let stale = fs::read(&tbi).unwrap();
    bgzip(&gz, &full);
    touch(&tbi, &gz, false);

    let (out, err) = chrm_73(&gz_s);
    assert_eq!(out, want, "the stale index hid the chrM call");
    assert!(
        err.contains(&format!(
            "NOTE: the index {} is older than {gz_s}; reading the whole file",
            tbi.display()
        )),
        "{err}"
    );

    // The same index with a newer time is read, which shows the indexed path
    // is the one taken: it lacks chrM, so the answer is not the call.
    touch(&tbi, &gz, true);
    let (out, err) = chrm_73(&gz_s);
    assert_ne!(out, want, "a fresh index is read, not the whole file");
    assert!(err.is_empty(), "{err}");

    // A fresh index of the current VCF answers the call, with no note.
    index(&gz, &tbi);
    assert_ne!(fs::read(&tbi).unwrap(), stale);
    touch(&tbi, &gz, true);
    let (out, err) = chrm_73(&gz_s);
    assert_eq!(out, want);
    assert!(err.is_empty(), "{err}");

    fs::remove_dir_all(&dir).unwrap();
}
