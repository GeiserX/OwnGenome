//! og's clinical installers (install/clinical.rs) and the raw-download
//! helpers (install/keep.rs) against the bash they replace: clinical_disk_ok,
//! adopt and done_with run from scripts/lib/common.sh's own text with the
//! same inputs and compare byte for byte; install_vep_cache and
//! install_pcgr_bundle run on small archives served from file:// URLs.
//!
//! The bash side runs from this checkout; each test that reads common.sh's
//! text goes once setup's bash is gone.

use og_setup::disk::{self, Probe};
use og_setup::fetch::{Kind, digest};
use og_setup::install::{Run, Vars, clinical, keep};
use std::collections::BTreeMap;
use std::fs;
use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};
use std::process::Command;

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../..")
        .canonicalize()
        .unwrap()
}

fn scratch(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("og-setup-clinical-{name}-{}", std::process::id()));
    let _ = fs::remove_dir_all(&d);
    fs::create_dir_all(&d).unwrap();
    d.canonicalize().unwrap()
}

/// The text of shell function NAME in scripts/lib/common.sh.
fn common_fn(name: &str) -> String {
    let t = fs::read_to_string(root().join("scripts/lib/common.sh")).unwrap();
    let start = t
        .find(&format!("\n{name}() {{\n"))
        .unwrap_or_else(|| panic!("common.sh has no {name}()"));
    let body = &t[start + 1..];
    let end = body.find("\n}\n").unwrap();
    body[..end + 3].to_owned()
}

/// bash -c SCRIPT in DIR under OG_SETUP=bash, with PATH (fake tools first)
/// and ENV: exit code, stdout, stderr.
fn bash(
    script: &str,
    dir: &Path,
    path_first: Option<&Path>,
    env: &[(&str, &str)],
) -> (i32, String, String) {
    let path = std::env::var("PATH").unwrap_or_default();
    let path = match path_first {
        Some(p) => format!("{}:{path}", p.display()),
        None => path,
    };
    let mut c = Command::new("bash");
    c.arg("-c").arg(script).current_dir(dir).env_clear();
    c.env("PATH", path).env("OG_SETUP", "bash");
    for (k, v) in env {
        c.env(k, v);
    }
    let o = c.output().unwrap();
    (
        o.status.code().unwrap_or(255),
        String::from_utf8_lossy(&o.stdout).into_owned(),
        String::from_utf8_lossy(&o.stderr).into_owned(),
    )
}

fn run_for(genome: &Path, downloads: &Path, keep_downloads: bool) -> Run {
    let mut m = BTreeMap::new();
    m.insert("PATH".to_owned(), std::env::var("PATH").unwrap_or_default());
    m.insert(
        "OG_DOWNLOADS".to_owned(),
        downloads.to_string_lossy().into_owned(),
    );
    m.insert("FETCH_TRIES".to_owned(), "1".to_owned());
    m.insert("FETCH_WAIT".to_owned(), "0".to_owned());
    let mut err = Vec::new();
    Run::new(
        &root(),
        Vars::from_map(m),
        "./scripts/setup.sh",
        &genome.to_string_lossy(),
        keep_downloads,
        &mut err,
    )
    .unwrap_or_else(|c| panic!("Run::new exited {c}: {}", String::from_utf8_lossy(&err)))
}

// --- clinical_disk_ok ---------------------------------------------------------

/// What the fake df and du of the bash side print.
struct Fixed(Option<String>, Option<String>);

impl Probe for Fixed {
    fn free_kib(&self, _: &Path) -> Option<String> {
        self.0.clone()
    }
    fn used_kib(&self, _: &Path) -> Option<String> {
        self.1.clone()
    }
}

/// og's and bash's clinical_disk_ok on one disk: status and stdout, with
/// free KiB FREE (None: df prints nothing usable), an archive.part whose
/// allocated KiB is USED (None: no partial download).
fn disk_both(name: &str, free: Option<&str>, used: Option<&str>, peak: u64) -> [(bool, String); 2] {
    let d = scratch(name);
    let bin = d.join("bin");
    fs::create_dir_all(&bin).unwrap();
    let df = match free {
        Some(f) => format!(
            "#!/bin/sh\necho 'Filesystem 1024-blocks Used Available Capacity Mounted on'\necho 'fakefs 1 0 {f} 0% /'\n"
        ),
        None => "#!/bin/sh\nexit 1\n".to_owned(),
    };
    let du = format!(
        "#!/bin/sh\nprintf '%s\\t%s\\n' '{}' \"$2\"\n",
        used.unwrap_or("0")
    );
    for (n, t) in [("df", df), ("du", du)] {
        fs::write(bin.join(n), t).unwrap();
        fs::set_permissions(bin.join(n), fs::Permissions::from_mode(0o755)).unwrap();
    }
    let dir = d.join("vep_cache");
    let archive = d.join("dl/cache.tar.gz");
    fs::create_dir_all(d.join("dl")).unwrap();
    if used.is_some() {
        fs::write(d.join("dl/cache.tar.gz.part"), b"x").unwrap();
    }
    let script = format!(
        "{}\nclinical_disk_ok {} {peak} 'VEP 116 cache' {}",
        common_fn("clinical_disk_ok"),
        dir.display(),
        archive.display()
    );
    let (code, said, _) = bash(&script, &d, Some(&bin), &[("GENOME_DIR", "/g")]);
    let mut out = Vec::new();
    let ok = disk::clinical_disk_ok(
        &Fixed(free.map(str::to_owned), used.map(str::to_owned)),
        &dir,
        peak,
        "VEP 116 cache",
        Some(&archive),
        "/g",
        &mut out,
    );
    let _ = fs::remove_dir_all(&d);
    [(code == 0, said), (ok, String::from_utf8(out).unwrap())]
}

/// The skip lines are bash's to the byte, at and around the 5 GB margin,
/// with a partial download counted, at the rounding ties of awk's %.1f,
/// and when df cannot be read.
#[test]
fn the_disk_check_is_clinical_disk_ok_s() {
    let need = (56 + 5) * 1_048_576u64;
    let cases: Vec<(String, Option<String>, Option<String>)> = vec![
        ("exact".into(), Some(need.to_string()), None),
        ("one-short".into(), Some((need - 1).to_string()), None),
        ("zero".into(), Some("0".into()), None),
        // 0.25 and 0.75 GB: ties, where a rounding rule shows.
        ("tie-low".into(), Some("262144".into()), None),
        ("tie-high".into(), Some("786432".into()), None),
        (
            "tie-big".into(),
            Some((40 * 1_048_576 + 262_144).to_string()),
            None,
        ),
        (
            "part-fits".into(),
            Some((need - 100).to_string()),
            Some("100".into()),
        ),
        (
            "part-short".into(),
            Some((need - 100).to_string()),
            Some("99".into()),
        ),
        ("no-df".into(), None, None),
        ("bad-df".into(), Some("12a".into()), None),
        ("bad-du".into(), Some(need.to_string()), Some("x".into())),
    ];
    for (name, free, used) in &cases {
        let [b, o] = disk_both(name, free.as_deref(), used.as_deref(), 56);
        assert_eq!(o, b, "{name}: og (left) and bash (right) differ");
    }
}

// --- adopt and done_with ------------------------------------------------------

/// One state of a file: absent, a checked download ("good"), a file that
/// fails the check ("bad"), a folder, or a partial download of N bytes.
#[derive(Clone, Copy, Debug)]
enum F {
    No,
    Good,
    Bad,
    Dir,
    Part(usize),
}

fn plant(p: &Path, f: F) {
    match f {
        F::No => {}
        F::Good => fs::write(p, "good").unwrap(),
        F::Bad => fs::write(p, "bad").unwrap(),
        F::Dir => fs::create_dir_all(p).unwrap(),
        F::Part(n) => fs::write(p, vec![b'p'; n]).unwrap(),
    }
}

/// Every file under DIR with its content, sorted.
fn tree(dir: &Path) -> Vec<String> {
    let mut v = Vec::new();
    let mut stack = vec![dir.to_path_buf()];
    while let Some(d) = stack.pop() {
        for e in fs::read_dir(&d).unwrap().flatten() {
            let p = e.path();
            let rel = p.strip_prefix(dir).unwrap().display().to_string();
            if p.is_dir() {
                v.push(format!("{rel}/"));
                stack.push(p);
            } else {
                v.push(format!(
                    "{rel} {}",
                    fs::read_to_string(&p).unwrap_or_default()
                ));
            }
        }
    }
    v.sort();
    v
}

/// adopt against common.sh's, over every mix of what OLD, OLD.part, NEW and
/// NEW.part can be, with and without a checksum: the same messages and the
/// same files afterwards. fetch's check is "the file holds good".
#[test]
fn adopt_is_common_sh_s() {
    let d = scratch("adopt");
    let fake_fetch = "fetch() { [ -f \"$2\" ] && [ \"$(cat \"$2\")\" = good ]; }\n";
    let script_base = format!("{fake_fetch}{}\n", common_fn("adopt"));
    let olds = [F::No, F::Good, F::Bad, F::Dir];
    let parts = [F::No, F::Part(3), F::Part(9)];
    let news = [F::No, F::Good, F::Bad];
    let mut n = 0;
    for sum in [true, false] {
        for &o in &olds {
            for &op in &parts {
                for &nw in &news {
                    for &np in &parts {
                        let mut seen = Vec::new();
                        for side in ["bash", "og"] {
                            let _ = fs::remove_dir_all(&d);
                            fs::create_dir_all(d.join("old")).unwrap();
                            let old = d.join("old/x.tgz");
                            let new = d.join("dl/x.tgz");
                            plant(&old, o);
                            plant(&d.join("old/x.tgz.part"), op);
                            if !matches!(nw, F::No) || !matches!(np, F::No) {
                                fs::create_dir_all(d.join("dl")).unwrap();
                            }
                            plant(&new, nw);
                            plant(&d.join("dl/x.tgz.part"), np);
                            let (o_s, n_s) = (old.to_string_lossy(), new.to_string_lossy());
                            let said = if side == "bash" {
                                let args = if sum {
                                    " https://x/x.tgz sha256 v"
                                } else {
                                    " https://x/x.tgz"
                                };
                                let s = format!("{script_base}adopt {o_s} {n_s}{args}");
                                let (code, said, e) = bash(&s, &d, None, &[]);
                                assert_eq!(code, 0, "bash adopt failed: {e}");
                                said
                            } else {
                                let mut out = Vec::new();
                                let verify =
                                    |f: &str| fs::read_to_string(f).is_ok_and(|t| t == "good");
                                keep::adopt_with(&o_s, &n_s, true, sum, &verify, &mut out);
                                String::from_utf8(out).unwrap()
                            };
                            seen.push((said, tree(&d)));
                        }
                        assert_eq!(
                            seen[1], seen[0],
                            "adopt with old {o:?} old.part {op:?} new {nw:?} new.part {np:?} sum {sum}: og (left), bash (right)"
                        );
                        n += 1;
                    }
                }
            }
        }
    }
    assert_eq!(n, 216);
    let _ = fs::remove_dir_all(&d);
}

/// done_with against common.sh's: kept by default, deleted with its .part
/// and .checksum under SETUP_KEEP_DOWNLOADS=false, and the same status.
#[test]
fn done_with_is_common_sh_s() {
    let d = scratch("done");
    for keep_all in [true, false] {
        for present in [true, false] {
            let mut seen = Vec::new();
            for side in ["bash", "og"] {
                let _ = fs::remove_dir_all(&d);
                fs::create_dir_all(d.join("dl")).unwrap();
                let f = d.join("dl/a.tgz");
                if present {
                    for x in ["", ".part", ".checksum"] {
                        fs::write(format!("{}{x}", f.display()), "x").unwrap();
                    }
                }
                let keep_v = if keep_all { "true" } else { "false" };
                let fs_ = f.to_string_lossy().into_owned();
                let (ok, said) = if side == "bash" {
                    let s = format!("{}\ndone_with {fs_}", common_fn("done_with"));
                    let (code, said, _) = bash(&s, &d, None, &[("SETUP_KEEP_DOWNLOADS", keep_v)]);
                    (code == 0, said)
                } else {
                    let r = run_for(&d.join("g"), &d.join("dl"), keep_all);
                    let mut out = Vec::new();
                    let ok = keep::done_with(&r, &[&fs_], &mut out);
                    (ok, String::from_utf8(out).unwrap())
                };
                seen.push((ok, said, tree(&d.join("dl"))));
            }
            assert_eq!(
                seen[1], seen[0],
                "done_with keep {keep_all} present {present}"
            );
        }
    }
    let _ = fs::remove_dir_all(&d);
}

// --- the installers, on small archives ------------------------------------------

/// A gzipped tar of SRC's entries NAMES, at OUT.
fn tgz(src: &Path, names: &[&str], out: &Path) {
    let s = Command::new("tar")
        .arg("-czf")
        .arg(out)
        .arg("-C")
        .arg(src)
        .args(names)
        .status()
        .unwrap();
    assert!(s.success());
}

/// A VEP cache tarball holding homo_sapiens/TOP (with info.txt unless
/// NO_INFO), served at file://SERVED/NAME with its CHECKSUMS beside it:
/// (url, sums url).
fn served_cache(served: &Path, name: &str, top: &str, no_info: bool) -> (String, String) {
    let src = served.join(format!("src-{name}"));
    let c = src.join("homo_sapiens").join(top);
    fs::create_dir_all(c.join("1")).unwrap();
    fs::write(c.join("1/all_vars.gz"), b"fake cache\n").unwrap();
    if !no_info {
        fs::write(c.join("info.txt"), b"species homo_sapiens\n").unwrap();
    }
    let tar = served.join(name);
    tgz(&src, &["homo_sapiens"], &tar);
    let sum = digest(Kind::Sum, &tar).unwrap();
    fs::write(served.join("CHECKSUMS"), format!("{sum} 1 {name}\n")).unwrap();
    (
        format!("file://{}", tar.display()),
        format!("file://{}", served.join("CHECKSUMS").display()),
    )
}

fn staging_left(dir: &Path) -> Vec<String> {
    fs::read_dir(dir)
        .map(|rd| {
            rd.flatten()
                .map(|e| e.file_name().to_string_lossy().into_owned())
                .filter(|n| n.starts_with(".extract."))
                .collect()
        })
        .unwrap_or_default()
}

#[test]
fn a_cache_is_unpacked_checked_and_moved_into_place() {
    let d = scratch("cache");
    let served = d.join("served");
    fs::create_dir_all(&served).unwrap();
    let name = "homo_sapiens_vep_116_GRCh38.tar.gz";
    let (url, sums) = served_cache(&served, name, "116_GRCh38", false);
    let g = d.join("g");
    let dir = format!("{}/vep_cache", g.display());
    fs::create_dir_all(format!("{dir}/.extract.old/x")).unwrap();
    let dl = d.join("dl");
    let tarball = format!("{}/{name}", dl.display());
    let r = run_for(&g, &dl, true);
    let (mut out, mut err) = (Vec::new(), Vec::new());
    assert!(
        clinical::install_cache_from(&r, &dir, "116", &url, &tarball, &sums, &mut out, &mut err),
        "{}",
        String::from_utf8_lossy(&err)
    );
    let fin = format!("{dir}/homo_sapiens/116_GRCh38");
    assert!(Path::new(&format!("{fin}/info.txt")).is_file());
    assert!(Path::new(&format!("{fin}/1/all_vars.gz")).is_file());
    assert!(
        staging_left(Path::new(&dir)).is_empty(),
        "the staging folders are gone"
    );
    assert!(
        Path::new(&tarball).is_file(),
        "the checked tarball stays in the downloads folder"
    );
    assert_eq!(
        String::from_utf8(out).unwrap(),
        format!("  Extracting {name}...\n")
    );
    // Again: nothing to do, and the same line as bash.
    let mut out = Vec::new();
    assert!(clinical::install_cache_from(
        &r, &dir, "116", &url, &tarball, &sums, &mut out, &mut err
    ));
    assert_eq!(
        String::from_utf8(out).unwrap(),
        "[OK] VEP 116 cache already present.\n"
    );
    // SETUP_KEEP_DOWNLOADS=false deletes the tarball once installed.
    fs::remove_dir_all(format!("{dir}/homo_sapiens")).unwrap();
    let r = run_for(&g, &dl, false);
    let mut out = Vec::new();
    assert!(clinical::install_cache_from(
        &r, &dir, "116", &url, &tarball, &sums, &mut out, &mut err
    ));
    assert!(!Path::new(&tarball).exists());
    assert!(String::from_utf8(out).unwrap().ends_with(&format!(
        "  Deleted the raw download {tarball} (SETUP_KEEP_DOWNLOADS=false).\n"
    )));
    let _ = fs::remove_dir_all(&d);
}

/// A tarball whose top folder is another release, or that has no
/// info.txt, is refused: nothing is moved into place and the staging
/// folder goes.
#[test]
fn a_tarball_of_another_layout_is_refused() {
    for (case, top, no_info) in [
        ("other-release", "115_GRCh38", false),
        ("no-info", "116_GRCh38", true),
    ] {
        let d = scratch(case);
        let served = d.join("served");
        fs::create_dir_all(&served).unwrap();
        let name = "homo_sapiens_vep_116_GRCh38.tar.gz";
        let (url, sums) = served_cache(&served, name, top, no_info);
        let g = d.join("g");
        let dir = format!("{}/vep_cache", g.display());
        let tarball = format!("{}/dl/{name}", d.display());
        let r = run_for(&g, &d.join("dl"), true);
        let (mut out, mut err) = (Vec::new(), Vec::new());
        assert!(!clinical::install_cache_from(
            &r, &dir, "116", &url, &tarball, &sums, &mut out, &mut err
        ));
        assert!(
            !Path::new(&format!("{dir}/homo_sapiens/116_GRCh38")).exists(),
            "{case}: moved into place"
        );
        assert!(
            !Path::new(&format!("{dir}/homo_sapiens/115_GRCh38")).exists(),
            "{case}: moved into place"
        );
        assert!(
            staging_left(Path::new(&dir)).is_empty(),
            "{case}: staging left"
        );
        assert!(
            String::from_utf8_lossy(&err).contains(&format!(
                "ERROR: could not unpack {tarball} into a 116_GRCh38 cache.\n"
            )),
            "{case}: {}",
            String::from_utf8_lossy(&err)
        );
        let _ = fs::remove_dir_all(&d);
    }
}

/// A cache folder there without info.txt stops with bash's message, before
/// any download.
#[test]
fn a_cache_folder_without_info_txt_stops() {
    let d = scratch("noinfo");
    let g = d.join("g");
    let dir = format!("{}/vep_cache", g.display());
    let fin = format!("{dir}/homo_sapiens/116_GRCh38");
    fs::create_dir_all(&fin).unwrap();
    let dl = d.join("dl");
    let r = run_for(&g, &dl, true);
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let tarball = format!("{}/x.tar.gz", dl.display());
    assert!(!clinical::install_cache_from(
        &r,
        &dir,
        "116",
        "file:///nonexistent/x.tar.gz",
        &tarball,
        "file:///nonexistent/CHECKSUMS",
        &mut out,
        &mut err
    ));
    assert_eq!(
        String::from_utf8(err).unwrap(),
        format!(
            "ERROR: {fin} exists but has no info.txt (an interrupted extraction?).\n  Remove it and run again.\n"
        )
    );
    assert!(out.is_empty());
    assert!(!dl.exists(), "nothing was downloaded");
    let _ = fs::remove_dir_all(&d);
}

/// The PCGR bundle: checked by sha256, unpacked, the staging folder becomes
/// the bundle; a tarball without data/grch38, or a bundle folder already
/// there without data, is refused.
#[test]
fn the_pcgr_bundle_is_checked_and_moved_into_place() {
    let d = scratch("pcgr");
    let src = d.join("src");
    fs::create_dir_all(src.join("data/grch38")).unwrap();
    fs::write(src.join("data/grch38/README"), b"bundle\n").unwrap();
    let tar = d.join("pcgr.tgz");
    tgz(&src, &["data"], &tar);
    let sha = digest(Kind::Sha256, &tar).unwrap();
    let url = format!("file://{}", tar.display());
    let g = d.join("g");
    let dir = format!("{}/pcgr_data", g.display());
    let tarball = format!("{}/dl/pcgr.tgz", d.display());
    let r = run_for(&g, &d.join("dl"), true);
    let (mut out, mut err) = (Vec::new(), Vec::new());
    assert!(
        clinical::install_bundle_from(
            &r, &dir, "20260620", &url, &tarball, &sha, &mut out, &mut err
        ),
        "{}",
        String::from_utf8_lossy(&err)
    );
    assert!(Path::new(&format!("{dir}/20260620/data/grch38/README")).is_file());
    assert_eq!(
        String::from_utf8(out).unwrap(),
        "[OK] PCGR/CPSR 20260620 bundle installed.\n"
    );
    assert!(staging_left(Path::new(&dir)).is_empty());
    // A bundle folder without data is never taken for a bundle.
    fs::remove_dir_all(format!("{dir}/20260620/data")).unwrap();
    let mut err = Vec::new();
    assert!(!clinical::install_bundle_from(
        &r,
        &dir,
        "20260620",
        &url,
        &tarball,
        &sha,
        &mut Vec::new(),
        &mut err
    ));
    // fetch first says it uses the kept, checked archive.
    assert!(String::from_utf8(err).unwrap().ends_with(&format!(
        "ERROR: {dir}/20260620 exists without data; remove it and retry.\n"
    )));
    assert!(staging_left(Path::new(&dir)).is_empty());
    // A tarball without data/grch38 is refused.
    let src2 = d.join("src2");
    fs::create_dir_all(src2.join("data/grch37")).unwrap();
    fs::write(src2.join("data/grch37/README"), b"x\n").unwrap();
    let tar2 = d.join("pcgr2.tgz");
    tgz(&src2, &["data"], &tar2);
    let sha2 = digest(Kind::Sha256, &tar2).unwrap();
    fs::remove_dir_all(format!("{dir}/20260620")).unwrap();
    let tarball2 = format!("{}/dl/pcgr2.tgz", d.display());
    let mut err = Vec::new();
    assert!(!clinical::install_bundle_from(
        &r,
        &dir,
        "20260620",
        &format!("file://{}", tar2.display()),
        &tarball2,
        &sha2,
        &mut Vec::new(),
        &mut err
    ));
    assert!(String::from_utf8_lossy(&err).contains(&format!(
        "ERROR: could not unpack {tarball2} into PCGR data/grch38."
    )));
    assert!(!Path::new(&format!("{dir}/20260620")).exists());
    assert!(staging_left(Path::new(&dir)).is_empty());
    let _ = fs::remove_dir_all(&d);
}

/// `og setup install --part` takes the three helpers' arguments only.
#[test]
fn the_part_forms_are_checked() {
    for bad in [
        vec!["--part"],
        vec!["--part", "vep-cache", "dir"],
        vec!["--part", "disk-ok", "dir", "56"],
        vec!["--part", "pcgr-bundle"],
        vec!["--part", "nope", "x"],
    ] {
        let args: Vec<String> = bad.iter().map(|s| (*s).to_owned()).collect();
        let (mut out, mut err) = (Vec::new(), Vec::new());
        assert_eq!(
            clinical::part_cli(&args, &mut out, &mut err),
            Some(2),
            "{bad:?}"
        );
        assert!(String::from_utf8_lossy(&err).starts_with("usage: og setup install --part"));
    }
    let args: Vec<String> = ["--kir-data", "/g"]
        .iter()
        .map(|s| (*s).to_owned())
        .collect();
    assert_eq!(
        clinical::part_cli(&args, &mut Vec::new(), &mut Vec::new()),
        None
    );
}

/// common.sh's three helpers hand the work to og unless OG_SETUP=bash, as
/// their first line.
#[test]
fn the_common_sh_helpers_are_thin_calls() {
    for (f, part) in [
        ("install_vep_cache", "vep-cache"),
        ("clinical_disk_ok", "disk-ok"),
        ("install_pcgr_bundle", "pcgr-bundle"),
    ] {
        let body = common_fn(f);
        let first = body.lines().nth(1).unwrap_or_default();
        assert_eq!(
            first,
            format!(
                "  [ \"${{OG_SETUP:-og}}\" = bash ] || {{ og_setup install --part {part} \"$@\"; return; }}"
            ),
            "{f}"
        );
    }
}
