//! The installers of install/tools.rs against the bash they replace: the
//! pins and container scripts og holds against setup.sh's text, the awk
//! programs of `--locityper-data` against awk itself, the engine memory
//! check of `--megane-kmers` against docker's and podman's `info`, and
//! whole runs of setup.sh's bash next to og's on the paths a test machine
//! cannot reach (an engine too small, a PyPGx tag at another commit).
//!
//! The bash side runs from this checkout; each test that reads setup.sh's
//! text says so and goes once setup.sh is a wrapper.

use og_setup::install::{self, Run, Vars, tools};
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
    let d = std::env::temp_dir().join(format!("og-setup-tools-{name}-{}", std::process::id()));
    let _ = fs::remove_dir_all(&d);
    fs::create_dir_all(&d).unwrap();
    d.canonicalize().unwrap()
}

fn setup_sh_text() -> String {
    fs::read_to_string(root().join("scripts/setup.sh")).unwrap()
}

/// The body of setup.sh's bash function NAME (`NAME() {` or `NAME() (`
/// to the closing line), or `None` once setup.sh no longer has it.
fn bash_fn(name: &str) -> Option<String> {
    let text = setup_sh_text();
    let start = text
        .lines()
        .position(|l| l == format!("{name}() {{") || l == format!("{name}() ("))?;
    let mut body = Vec::new();
    for l in text.lines().skip(start) {
        body.push(l);
        if body.len() > 1 && (l == "}" || l == ")") {
            break;
        }
    }
    Some(body.join("\n"))
}

fn write_exe(path: &Path, body: &str) {
    fs::write(path, body).unwrap();
    fs::set_permissions(path, fs::Permissions::from_mode(0o755)).unwrap();
}

/// `awk PROGRAM FILES...` with AWK_ARGS before the program: its stdout.
fn awk(args: &[&str], program: &str, files: &[&Path]) -> String {
    let o = Command::new("awk")
        .args(args)
        .arg(program)
        .args(files)
        .output()
        .unwrap();
    String::from_utf8_lossy(&o.stdout).into_owned()
}

// --- what og holds against setup.sh's text ------------------------------------------

#[test]
fn the_parascopy_pins_are_setup_sh_s() {
    let text = setup_sh_text();
    let Some(i) = text.find("PARASCOPY_FILES=\"") else {
        return; // setup.sh is a wrapper now: og's list is the only one.
    };
    let rest = &text[i + "PARASCOPY_FILES=\"".len()..];
    let list = &rest[..rest.find('"').unwrap()];
    let got: Vec<(String, String)> = list
        .lines()
        .map(|l| {
            let mut w = l.split_whitespace();
            (w.next().unwrap().to_owned(), w.next().unwrap().to_owned())
        })
        .collect();
    let want: Vec<(String, String)> = tools::PARASCOPY_FILES
        .iter()
        .map(|(a, b)| ((*a).to_owned(), (*b).to_owned()))
        .collect();
    assert_eq!(got, want);
}

/// The text between `-c '` and the next `'` in the N-th `bash -euo
/// pipefail -c '` of install_locityper_data.
fn locityper_script(n: usize) -> Option<String> {
    let body = bash_fn("install_locityper_data")?;
    let mut rest = body.as_str();
    let mut found = None;
    for _ in 0..=n {
        let i = rest.find("bash -euo pipefail -c '")? + "bash -euo pipefail -c '".len();
        let tail = &rest[i..];
        let j = tail.find('\'')?;
        found = Some(tail[..j].to_owned());
        rest = &tail[j + 1..];
    }
    found
}

#[test]
fn the_locityper_container_scripts_are_setup_sh_s() {
    let Some(panel) = locityper_script(0) else {
        return; // setup.sh is a wrapper now.
    };
    assert_eq!(panel, tools::LOCITYPER_PANEL_SCRIPT);
    assert_eq!(
        locityper_script(1).unwrap(),
        tools::LOCITYPER_BUILD_SCRIPT
    );
}

// --- the awk programs of --locityper-data, against awk ---------------------------------

/// setup.sh's programs, as install_locityper_data spells them.
const AWK_BAD: &str = r#"NR == FNR {len[$1] = $2; next}
    !/^#/ && NF && (NF < 4 || $2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ || $4 == "" || $2 + 0 >= $3 + 0 || !($1 in len) || $3 + 0 > len[$1] + 0) {
      print FILENAME ":" FNR; exit}"#;
const AWK_SIZE: &str = r#"{n += $2} END {n = n < 10000000 ? 10000000 : n; printf "%.0f", n}"#;
const AWK_KEEP: &str = r#"!/^#/ && NF {print $1, $2, $3, $4}"#;
const AWK_PAD: &str = r#"{s = $2 - 250000; if (s < 0) s = 0; print $1, s, $3 + 250000}"#;
const AWK_MERGE: &str = r#"$1 == c && $2 <= e {if ($3 > e) e = $3; next} {if (c != "") print c ":" s + 1 "-" e; c = $1; s = $2; e = $3}
                  END {if (c != "") print c ":" s + 1 "-" e}"#;

#[test]
fn the_awk_programs_are_setup_sh_s() {
    let Some(body) = bash_fn("install_locityper_data") else {
        return; // setup.sh is a wrapper now: these copies are the record.
    };
    for p in [AWK_BAD, AWK_SIZE, AWK_KEEP, AWK_PAD, AWK_MERGE] {
        assert!(body.contains(p), "install_locityper_data no longer holds:\n{p}");
    }
}

/// GRCh38's primary contigs (the no-ALT set's chromosomes), a 235 Mb
/// reference, and one under the 10M floor.
fn fais() -> Vec<(&'static str, String)> {
    let grch38: [(&str, u64); 25] = [
        ("chr1", 248_956_422),
        ("chr2", 242_193_529),
        ("chr3", 198_295_559),
        ("chr4", 190_214_555),
        ("chr5", 181_538_259),
        ("chr6", 170_805_979),
        ("chr7", 159_345_973),
        ("chr8", 145_138_636),
        ("chr9", 138_394_717),
        ("chr10", 133_797_422),
        ("chr11", 135_086_622),
        ("chr12", 133_275_309),
        ("chr13", 114_364_328),
        ("chr14", 107_043_718),
        ("chr15", 101_991_189),
        ("chr16", 90_338_345),
        ("chr17", 83_257_441),
        ("chr18", 80_373_285),
        ("chr19", 58_617_616),
        ("chr20", 64_444_167),
        ("chr21", 46_709_983),
        ("chr22", 50_818_468),
        ("chrX", 156_040_895),
        ("chrY", 57_227_415),
        ("chrM", 16_569),
    ];
    let fai = |rows: &[(&str, u64)]| -> String {
        rows.iter()
            .map(|(c, n)| format!("{c}\t{n}\t6\t60\t61\n"))
            .collect()
    };
    vec![
        ("GRCh38", fai(&grch38)),
        (
            "235Mb",
            fai(&[("chr20", 64_444_167), ("chr21", 46_709_983), ("chr1", 123_845_850)]),
        ),
        ("floor", fai(&[("chr1", 1_000_000)])),
        ("one line, no newline", "chr1\t2500000000\t6\t60\t61".to_owned()),
    ]
}

#[test]
fn jellyfish_s_hash_size_is_awk_s() {
    let d = scratch("hsize");
    for (name, text) in fais() {
        let f = d.join("ref.fai");
        fs::write(&f, &text).unwrap();
        let want = awk(&[], AWK_SIZE, &[&f]);
        assert_eq!(tools::jellyfish_size(&text), want, "{name}");
    }
    let total = tools::jellyfish_size(&fais()[0].1);
    assert_eq!(total, "3088286401", "GRCh38's primary contigs");
    assert_eq!(tools::jellyfish_size(&fais()[1].1), "235000000");
    let _ = fs::remove_dir_all(&d);
}

#[test]
fn the_loci_check_refuses_the_rows_awk_refuses() {
    let d = scratch("loci");
    let fai = d.join("ref.fai");
    fs::write(&fai, "chr1\t1000\t6\t60\t61\nchr2\t500\t6\t60\t61\n").unwrap();
    let head = "# a comment\n\nchr1\t10\t20\tGOOD\n";
    let bad_rows = [
        ("three columns", "chr1\t10\t20\n"),
        ("start not a number", "chr1\tx\t20\tA\n"),
        ("end not a number", "chr1\t10\t2e1\tA\n"),
        ("no name", "chr1\t10\t20\t\n"),
        ("start not below end", "chr1\t20\t20\tA\n"),
        ("contig not in the reference", "chr3\t10\t20\tA\n"),
        ("end beyond the contig", "chr2\t10\t501\tA\n"),
        ("a negative start", "chr1\t-1\t20\tA\n"),
        ("spaces, not tabs", "chr1 10 20 A\n"),
    ];
    let fai_text = fs::read_to_string(&fai).unwrap();
    let mut refused = 0;
    for (what, row) in bad_rows {
        let loci = d.join("loci.bed");
        let text = format!("{head}{row}chr1\t30\t40\tLATER\n");
        fs::write(&loci, &text).unwrap();
        let want = awk(&["-F\t"], AWK_BAD, &[&fai, &loci]);
        let name = loci.to_string_lossy().into_owned();
        let got = tools::bad_locus(&fai_text, &text, &name).unwrap_or_default();
        assert_eq!(got, want.trim_end_matches('\n'), "{what}");
        assert_eq!(got, format!("{name}:4"), "{what}: awk and og both name row 4");
        refused += 1;
    }
    assert!(refused >= 6);
    // Good rows, an end equal to the contig's length, comments and blank
    // lines: nothing refused by either.
    let loci = d.join("good.bed");
    let text = format!("{head}chr2\t0\t500\tEDGE\n# x\tq\n\nchr1\t1\t2\tB\textra\n");
    fs::write(&loci, &text).unwrap();
    assert_eq!(awk(&["-F\t"], AWK_BAD, &[&fai, &loci]), "");
    assert_eq!(tools::bad_locus(&fai_text, &text, "good.bed"), None);
    // The rows kept in loci.bed are awk's.
    let kept = awk(&["-F\t", "-v", "OFS=\t"], AWK_KEEP, &[&loci]);
    let og: String = tools::kept_loci(&text)
        .iter()
        .map(|l| format!("{l}\n"))
        .collect();
    assert_eq!(og, kept);
    let _ = fs::remove_dir_all(&d);
}

#[test]
fn the_panel_regions_are_the_shell_pipelines() {
    let d = scratch("regions");
    let loci = d.join("loci.bed.new");
    // Out of order, overlapping once padded, on several contigs, one start
    // under 250 kb, two rows with the same start.
    let text = "chr2\t1000000\t1000100\tA\nchr10\t5000000\t5000100\tB\nchr1\t100\t200\tC\nchr2\t1400000\t1400100\tD\nchr2\t3000000\t3000100\tE\nchrX\t200000\t200500\tF\nchr2\t3000000\t3200000\tG\nchr10\t1000\t2000\tH\n";
    fs::write(&loci, text).unwrap();
    let pipeline = format!(
        "awk -F'\\t' -v OFS='\\t' '{AWK_PAD}' \"$1\" | sort -k1,1 -k2,2n | awk -F'\\t' '{AWK_MERGE}' | tr '\\n' ' '"
    );
    let o = Command::new("bash")
        .args(["-c", &pipeline, "_"])
        .arg(&loci)
        .output()
        .unwrap();
    let want = String::from_utf8_lossy(&o.stdout).into_owned();
    let kept = tools::kept_loci(text);
    assert_eq!(tools::panel_regions(&kept), want);
    assert!(want.contains("chr2:750001-1650100 "), "{want}");
    let _ = fs::remove_dir_all(&d);
}

// --- the memory check of --megane-kmers -------------------------------------------------

/// An engine's answers to `info --format FIELD`: docker's (a number at
/// .MemTotal) or podman's (a template error there, the number at
/// .Host.MemTotal).
fn docker_info(bytes: &'static str) -> impl Fn(&str) -> Option<String> {
    move |f: &str| (f == "{{.MemTotal}}").then(|| bytes.to_owned())
}

fn podman_info(bytes: &'static str) -> impl Fn(&str) -> Option<String> {
    move |f: &str| (f == "{{.Host.MemTotal}}").then(|| bytes.to_owned())
}

#[test]
fn the_engine_memory_is_read_from_docker_and_from_podman() {
    assert_eq!(tools::engine_mem_gb(&docker_info("17179869184")), Some(17));
    assert_eq!(tools::engine_mem_gb(&podman_info("2147483648")), Some(2));
    // Docker Desktop and colima answer in bytes of the VM; 999 MB is 0 GB,
    // which the check reads as "unknown" and lets through.
    assert_eq!(tools::engine_mem_gb(&docker_info("999000000")), Some(0));
    // A podman that answers .MemTotal with no number (exit 0) is asked .Host.
    let odd = |f: &str| {
        Some(if f == "{{.MemTotal}}" {
            "<no value>".to_owned()
        } else {
            "8000000000".to_owned()
        })
    };
    assert_eq!(tools::engine_mem_gb(&odd), Some(8));
    assert_eq!(tools::engine_mem_gb(&|_: &str| None), None);
}

/// A container engine for the bash-and-og runs: `info --format` answers as
/// docker (FIELD .MemTotal) or podman (.Host.MemTotal) with BYTES; every
/// other call succeeds and prints nothing; each call is logged to LOG.
fn stub_engine(dir: &Path, name: &str, field: &str, bytes: &str) -> PathBuf {
    let log = dir.join(format!("{name}.log"));
    let other = if field == "{{.MemTotal}}" {
        "{{.Host.MemTotal}}"
    } else {
        "{{.MemTotal}}"
    };
    let path = dir.join(name);
    write_exe(
        &path,
        &format!(
            "#!/usr/bin/env bash\nprintf '%s\\n' \"$*\" >> '{}'\nif [ \"$1\" = info ]; then\n  case \"$*\" in\n    *'{field}'*) echo {bytes}; exit 0 ;;\n    *'{other}'*) echo 'Error: template: info:1:2: cannot evaluate field' >&2; exit 125 ;;\n  esac\nfi\nexit 0\n",
            log.display()
        ),
    );
    path
}

/// The environment both sides start with.
fn base_env(g: &Path, extra: &[(&str, String)]) -> BTreeMap<String, String> {
    let mut m = BTreeMap::new();
    m.insert("PATH".to_owned(), std::env::var("PATH").unwrap_or_default());
    m.insert("HOME".to_owned(), std::env::var("HOME").unwrap_or_default());
    m.insert("GENOME_DIR".to_owned(), g.to_string_lossy().into_owned());
    m.insert("SETUP_DISK_CHECK".to_owned(), "false".to_owned());
    m.insert("OG".to_owned(), "/usr/bin/true".to_owned());
    m.insert("FETCH_TRIES".to_owned(), "1".to_owned());
    for (k, v) in extra {
        m.insert((*k).to_owned(), v.clone());
    }
    m
}

/// setup.sh's bash on FLAG G: exit code, stdout, stderr.
fn bash_side(flag: &str, g: &Path, env: &BTreeMap<String, String>) -> (u8, String, String) {
    let mut c = Command::new("bash");
    c.arg(root().join("scripts/setup.sh"))
        .args([flag, &g.to_string_lossy()])
        .env_clear()
        .envs(env)
        .env("OG_SETUP", "bash");
    let o = c.output().unwrap();
    (
        o.status.code().unwrap_or(255) as u8,
        String::from_utf8_lossy(&o.stdout).into_owned(),
        String::from_utf8_lossy(&o.stderr).into_owned(),
    )
}

/// og's installer KEY on G with ENV: exit code, stdout, stderr.
fn og_side(key: &str, g: &Path, env: &BTreeMap<String, String>) -> (u8, String, String) {
    let script = root().join("scripts/setup.sh").to_string_lossy().into_owned();
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let run = Run::new(
        &root(),
        Vars::from_map(env.clone()),
        &script,
        &g.to_string_lossy(),
        true,
        &mut err,
    )
    .unwrap_or_else(|c| panic!("Run::new exited {c}: {}", String::from_utf8_lossy(&err)));
    let code = install::tools::run(key, &run, &mut out, &mut err).unwrap();
    (
        code,
        String::from_utf8_lossy(&out).into_owned(),
        String::from_utf8_lossy(&err).into_owned(),
    )
}

/// A reference of one tiny contig, with its .fai.
fn tiny_reference(g: &Path) {
    let r = g.join("reference");
    fs::create_dir_all(&r).unwrap();
    fs::write(r.join("GRCh38_no_alt_analysis_set.fasta"), ">chr1\nACGT\n").unwrap();
    fs::write(
        r.join("GRCh38_no_alt_analysis_set.fasta.fai"),
        "chr1\t4\t6\t4\t5\n",
    )
    .unwrap();
}

#[test]
fn an_engine_too_small_for_megane_is_refused_as_bash_refuses_it() {
    for (engine, field) in [("docker", "{{.MemTotal}}"), ("podman", "{{.Host.MemTotal}}")] {
        let d = scratch(&format!("megane-{engine}"));
        let stub = stub_engine(&d, engine, field, "2147483648");
        let mut sides = Vec::new();
        for side in ["bash", "og"] {
            let g = d.join(side);
            tiny_reference(&g);
            let env = base_env(&g, &[("CONTAINER_ENGINE", stub.to_string_lossy().into_owned())]);
            let got = if side == "bash" {
                bash_side("--megane-kmers", &g, &env)
            } else {
                og_side("megane-kmers", &g, &env)
            };
            let g_s = g.to_string_lossy().into_owned();
            sides.push((
                got.0,
                got.1.replace(&g_s, "G"),
                got.2.replace(&g_s, "G"),
            ));
            assert!(
                !g.join("reference/megane").exists(),
                "{engine} {side}: the refusal made reference/megane"
            );
        }
        assert_eq!(sides[0], sides[1], "{engine}: bash first, og second");
        assert_eq!(sides[1].0, 1, "{engine}");
        let stub_s = stub.to_string_lossy().into_owned();
        assert!(
            sides[1].1.contains(&format!(
                "[WARN] MEGAnE's k-mer set needs about 6 GB of RAM and {stub_s} has 2 GB."
            )),
            "{engine}: {}",
            sides[1].1
        );
        let _ = fs::remove_dir_all(&d);
    }
}

// --- --pypgx-data and the host's git ---------------------------------------------------

/// A git that clones by writing COMMIT into the clone's `commit` file and
/// answers `rev-parse HEAD` from it (setup-everything's fake git): `-C`
/// alone only for the fresh clone in setup's stage, and an installed bundle
/// only with `-c safe.directory=<it>`, as git asks of a folder another user
/// may own.
fn fake_git(dir: &Path, commit: &str) -> String {
    let bin = dir.join("bin");
    fs::create_dir_all(&bin).unwrap();
    write_exe(
        &bin.join("git"),
        &format!(
            r#"#!/usr/bin/env bash
set -eu
case "$1" in
  clone) dest=${{@: -1}}; mkdir -p "$dest"; printf '%s\n' '{commit}' > "$dest/commit" ;;
  -C) [[ "$2" = */.pypgx.*/bundle ]] && [ "$3 $4" = 'rev-parse HEAD' ] || exit 1; cat "$2/commit" ;;
  -c) [ "$2" = "safe.directory=$4" ] && [ "$3 $5 $6" = '-C rev-parse HEAD' ] || exit 1; cat "$4/commit" ;;
  *) exit 1 ;;
esac
"#
        ),
    );
    format!(
        "{}:{}",
        bin.display(),
        std::env::var("PATH").unwrap_or_default()
    )
}

fn pin(k: &str) -> String {
    let t = fs::read_to_string(root().join("versions.env")).unwrap();
    og_setup::fetch::versions_pins(&t).get(k).cloned().unwrap()
}

/// Both sides of `--pypgx-data` on their own G with a fake git cloning
/// COMMIT; each side's (code, stdout, stderr) with G masked, then G.
fn pypgx_both(d: &Path, commit: &str, seed: &dyn Fn(&Path)) -> Vec<((u8, String, String), PathBuf)> {
    let path = fake_git(d, commit);
    ["bash", "og"]
        .iter()
        .map(|side| {
            let g = d.join(side);
            fs::create_dir_all(&g).unwrap();
            seed(&g);
            let env = base_env(&g, &[("PATH", path.clone())]);
            let got = if *side == "bash" {
                bash_side("--pypgx-data", &g, &env)
            } else {
                og_side("pypgx-data", &g, &env)
            };
            let g_s = g.to_string_lossy().into_owned();
            (
                (got.0, got.1.replace(&g_s, "G"), got.2.replace(&g_s, "G")),
                g,
            )
        })
        .collect()
}

fn leftovers(g: &Path) -> Vec<String> {
    fs::read_dir(g.join("reference"))
        .map(|rd| {
            rd.filter_map(Result::ok)
                .map(|e| e.file_name().to_string_lossy().into_owned())
                .filter(|n| n.starts_with(".pypgx.") || n.ends_with(".setup.lock"))
                .collect()
        })
        .unwrap_or_default()
}

#[test]
fn a_pypgx_tag_at_another_commit_is_refused_as_bash_refuses_it() {
    let d = scratch("pypgx-wrong");
    let sides = pypgx_both(&d, "0000000000000000000000000000000000000000", &|_| {});
    assert_eq!(sides[0].0, sides[1].0, "bash first, og second");
    let (code, _, err) = &sides[1].0;
    assert_eq!(*code, 1);
    assert!(
        err.contains(&format!(
            "ERROR: PyPGx bundle tag {} differs from its pinned commit.",
            pin("PYPGX_BUNDLE_VERSION")
        )),
        "{err}"
    );
    for (_, g) in &sides {
        assert!(!g.join("reference/pypgx-bundle").exists());
        assert!(leftovers(g).is_empty(), "{:?}", leftovers(g));
    }
    let _ = fs::remove_dir_all(&d);
}

#[test]
fn the_pinned_pypgx_tag_installs_once_and_another_bundle_is_left_alone() {
    let commit = pin("PYPGX_BUNDLE_COMMIT");
    let d = scratch("pypgx-right");
    let sides = pypgx_both(&d, &commit, &|_| {});
    assert_eq!(sides[0].0, sides[1].0, "bash first, og second");
    assert_eq!(sides[1].0.0, 0);
    for (_, g) in &sides {
        let c = fs::read_to_string(g.join("reference/pypgx-bundle/commit")).unwrap();
        assert_eq!(c.trim_end(), commit);
        assert!(leftovers(g).is_empty(), "{:?}", leftovers(g));
    }
    let _ = fs::remove_dir_all(&d);

    // A bundle there at the pin: kept, nothing said. One at another
    // commit: refused, kept.
    for (there, want_code) in [(commit.as_str(), 0u8), ("1111111111111111111111111111111111111111", 1)] {
        let d = scratch("pypgx-there");
        let seed = |g: &Path| {
            let b = g.join("reference/pypgx-bundle");
            fs::create_dir_all(&b).unwrap();
            fs::write(b.join("commit"), format!("{there}\n")).unwrap();
        };
        let sides = pypgx_both(&d, &commit, &seed);
        assert_eq!(sides[0].0, sides[1].0, "bash first, og second ({there})");
        assert_eq!(sides[1].0.0, want_code);
        for (_, g) in &sides {
            let c = fs::read_to_string(g.join("reference/pypgx-bundle/commit")).unwrap();
            assert_eq!(c.trim_end(), there);
        }
        let _ = fs::remove_dir_all(&d);
    }
}

// --- paths ---------------------------------------------------------------------------

#[test]
fn the_megane_path_is_setup_sh_s() {
    let Some(body) = bash_fn("megane_kmers_path") else {
        return; // setup.sh is a wrapper now.
    };
    let forms: [(&str, &str); 5] = [
        ("/g/reference/GRCh38_no_alt_analysis_set.fasta", ""),
        ("/g/reference/my.ref.fa.gz", ""),
        ("/g/reference/noext", ""),
        ("/g/reference/x.fasta", "kmers/set.mk"),
        ("/g/reference/x.fasta", "/elsewhere/set.mk"),
    ];
    for (reference, mk) in forms {
        let o = Command::new("bash")
            .args(["-c", &format!("{body}\nmegane_kmers_path")])
            .env("REF_FASTA", reference)
            .env("MEGANE_KMERS", mk)
            .env("GENOME_DIR", "/g")
            .output()
            .unwrap();
        let want = String::from_utf8_lossy(&o.stdout).trim_end().to_owned();
        let m = (!mk.is_empty()).then_some(mk);
        assert_eq!(tools::megane_kmers_path(reference, m, "/g"), want, "{reference} {mk}");
    }
}

#[test]
fn the_masked_reference_lies_beside_its_source() {
    let url = "https://example.invalid/x/GCA_exclusions.bed";
    assert_eq!(
        tools::masked_paths("/g/reference/GRCh38_no_alt_analysis_set.fasta", url),
        (
            "/g/reference/GRCh38_no_alt_analysis_set_maskedGRC.fasta".to_owned(),
            "/g/reference/GRCh38_no_alt_analysis_set_maskedGRC.dict".to_owned(),
            "/g/reference/GCA_exclusions.bed".to_owned(),
        )
    );
    // `${src%.*}` cuts at the last '.', even one in a folder name.
    assert_eq!(
        tools::masked_paths("/g/ref.d/chr21", url).0,
        "/g/ref_maskedGRC.fasta"
    );
}
