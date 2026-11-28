//! og setup install's frame against the bash it replaces: the argv parse
//! against scripts/setup.sh's own (every form), the container argv against
//! run_in, the pinned tag against pinned_tag, the missing list against
//! print_missing, the pins og holds against setup.sh's, and the component
//! table against the source table and versions.env.
//!
//! The bash side runs from this checkout; each test that reads setup.sh's
//! text says so and goes once setup.sh is a wrapper.

use og_setup::components::{self, Component, Kind as CKind};
use og_setup::engine::{Place, RunIn};
use og_setup::fetch::{Ctx, Sources};
use og_setup::install::{self, Missing, Mode};
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
    let d = std::env::temp_dir().join(format!("og-setup-core-{name}-{}", std::process::id()));
    let _ = fs::remove_dir_all(&d);
    fs::create_dir_all(&d).unwrap();
    d.canonicalize().unwrap()
}

fn setup_sh() -> PathBuf {
    root().join("scripts/setup.sh")
}

/// The versions.env value of K, as the usage text reads it.
fn versions(k: &str) -> String {
    let t = fs::read_to_string(root().join("versions.env")).unwrap();
    og_setup::fetch::versions_pins(&t)
        .get(k)
        .cloned()
        .unwrap_or_default()
}

/// bash setup.sh ARGS with OG_SETUP=bash and ENV only: code, stdout, stderr.
fn bash_setup(args: &[&str], env: &[(&str, &str)]) -> (u8, String, String) {
    let mut c = Command::new("bash");
    c.arg(setup_sh()).args(args).env_clear();
    c.env("PATH", std::env::var("PATH").unwrap_or_default());
    c.env("OG_SETUP", "bash");
    // A run that gets past the parse by mistake starts no container.
    c.env("CONTAINER_ENGINE", "false");
    for (k, v) in env {
        c.env(k, v);
    }
    let o = c.output().unwrap();
    (
        o.status.code().unwrap_or(255) as u8,
        String::from_utf8_lossy(&o.stdout).into_owned(),
        String::from_utf8_lossy(&o.stderr).into_owned(),
    )
}

/// One run of setup.sh: its arguments, and the variables it starts with.
type Form<A, V> = (Vec<A>, Vec<(&'static str, V)>);

fn og_parse(args: &[&str], env: &[(&str, &str)]) -> Result<install::Parsed, install::Stop> {
    let args: Vec<String> = args.iter().map(|s| (*s).to_owned()).collect();
    let m: BTreeMap<String, String> = env
        .iter()
        .map(|(k, v)| ((*k).to_owned(), (*v).to_owned()))
        .collect();
    let envf = move |k: &str| m.get(k).cloned();
    install::parse(&args, &envf, &versions, &setup_sh().to_string_lossy())
}

/// The forms that stop in the parse: og prints and exits as setup.sh does.
/// The first form (no arguments) compares the whole usage text, so a line
/// setup.sh's usage gains that USAGE lacks turns this red.
#[test]
fn the_forms_that_stop_match_setup_sh() {
    let cases: Vec<Form<&str, &str>> = vec![
        (vec![], vec![]),
        (vec!["--kir-data"], vec![]),
        (vec!["--no-keep-downloads"], vec![]),
        (vec!["--no-clinical-data"], vec![]),
        (vec!["--everything"], vec![]),
        (vec!["--refresh"], vec![]),
        (vec!["--refresh", "dbsnp", "/g"], vec![]),
        (vec!["--refresh", "--kir-data"], vec![]),
        (vec!["--bogus", "/g"], vec![]),
        (vec!["--help"], vec![]),
        (
            vec!["--kir-data", "/g"],
            vec![("SETUP_KEEP_DOWNLOADS", "yes")],
        ),
    ];
    for (args, env) in cases {
        let (code, out, err) = bash_setup(&args, &env);
        let Err(stop) = og_parse(&args, &env) else {
            panic!("og parsed {args:?} where setup.sh stops with {code}: {err}");
        };
        assert_eq!(
            (stop.code, stop.stdout.as_str(), stop.stderr.as_str()),
            (code, out.as_str(), err.as_str()),
            "{args:?} {env:?}"
        );
    }
}

/// A second flag after one og runs is refused; setup.sh would read it as
/// the install folder.
#[test]
fn a_second_flag_is_refused() {
    for args in [
        vec!["--kir-data", "--viral-data"],
        vec!["--kir-data", "/g", "--viral-data"],
        vec!["--pull-only", "--kir-data"],
    ] {
        let stop = og_parse(&args, &[]).unwrap_err();
        assert_eq!(stop.code, 2, "{args:?}");
        assert!(
            stop.stderr.contains("one flag at a time"),
            "{}",
            stop.stderr
        );
    }
}

/// After a flag og hands back to setup.sh, a second flag is setup.sh's to
/// read, as it always was: og parses the form and hands it on (it refused
/// `--vep-cache /data --kir-data` with exit 2).
#[test]
fn a_second_flag_after_a_flag_bash_runs_goes_to_bash() {
    for args in [
        vec!["--vep-cache", "/g", "--kir-data"],
        vec!["/g", "--kir-data"],
        vec!["--refresh", "clinvar", "/g", "--everything"],
    ] {
        let p = og_parse(&args, &[])
            .unwrap_or_else(|s| panic!("{args:?} stopped with {}: {}", s.code, s.stderr));
        assert!(
            !og_setup::install::ported(p.mode.key()),
            "{args:?}: og runs {}",
            p.mode.key()
        );
    }
}

/// Every form setup.sh hands on: a fake og stands in for og, so setup.sh's
/// own parse runs and its dispatch block shows what it hands on (the flag,
/// the folder, SETUP_KEEP_DOWNLOADS). og's parse must give the same.
#[test]
fn every_handed_form_matches_setup_sh() {
    let d = scratch("forms");
    let fake = d.join("og");
    let rec = d.join("record");
    fs::write(
        &fake,
        "#!/usr/bin/env bash\ncase \"$1 $2\" in\n  'setup ported'|'setup plan') exit 0 ;;\n  'setup install') shift 2; { printf '%s\\n' \"keep=$SETUP_KEEP_DOWNLOADS\"; printf '%s\\n' \"$@\"; } > \"$(dirname \"$0\")/record\"; exit 0 ;;\nesac\nexit 3\n",
    )
    .unwrap();
    fs::set_permissions(&fake, fs::Permissions::from_mode(0o755)).unwrap();
    let g = d.join("genome").to_string_lossy().into_owned();
    let mut forms: Vec<Form<String, String>> = Vec::new();
    let mut flags: Vec<String> = install::OPT_IN.iter().map(|k| format!("--{k}")).collect();
    flags.extend(["--clinical-data", "--sample-qc-data", "--pull-only"].map(String::from));
    for f in &flags {
        forms.push((vec![f.clone(), g.clone()], vec![]));
    }
    forms.push((
        vec!["--no-keep-downloads".into(), "--kir-data".into(), g.clone()],
        vec![],
    ));
    forms.push((
        vec![
            "--viral-data".into(),
            g.clone(),
            "--no-keep-downloads".into(),
        ],
        vec![],
    ));
    forms.push((vec!["--chip-data".into()], vec![("GENOME_DIR", g.clone())]));
    forms.push((vec!["--pull-only".into()], vec![]));
    forms.push((vec!["--somatic-data".into(), format!("{g}/")], vec![]));
    forms.push((
        vec!["--cnvpytor-data".into(), g.clone()],
        vec![("SETUP_KEEP_DOWNLOADS", "false".into())],
    ));
    assert!(forms.len() >= 33, "{} forms", forms.len());
    for (args, env) in forms {
        let _ = fs::remove_file(&rec);
        let a: Vec<&str> = args.iter().map(String::as_str).collect();
        let mut e: Vec<(&str, &str)> = env.iter().map(|(k, v)| (*k, v.as_str())).collect();
        e.push(("OG", fake.to_str().unwrap()));
        e.push(("OG_SETUP", "og"));
        let (code, _, err) = bash_setup(&a, &e);
        assert_eq!(code, 0, "setup.sh {args:?}: {err}");
        let got = fs::read_to_string(&rec)
            .unwrap_or_else(|_| panic!("setup.sh {args:?} handed nothing on: {err}"));
        let p = og_parse(&a, &e).unwrap_or_else(|s| panic!("og refused {args:?}: {}", s.stderr));
        let mut want = format!("keep={}\n", p.keep_downloads);
        // setup.sh hands on its own words after the flag; og reads the same.
        for w in &p.argv {
            want.push_str(w);
            want.push('\n');
        }
        // GENOME_DIR without a trailing slash is setup.sh's, which og strips too.
        assert_eq!(
            got.trim_end_matches('\n'),
            want.trim_end_matches('\n'),
            "{args:?}"
        );
        let folder = if p.genome_dir.is_empty() {
            None
        } else {
            Some(p.genome_dir.trim_end_matches('/'))
        };
        if p.mode != Mode::PullOnly {
            assert_eq!(folder, Some(g.as_str()), "{args:?}");
        }
    }
    let _ = fs::remove_dir_all(&d);
}

/// The modes setup.sh keeps running itself parse to their own mode.
#[test]
fn the_other_modes_parse() {
    let m = |a: &[&str]| og_parse(a, &[]).unwrap().mode;
    assert_eq!(m(&["/g"]), Mode::Full);
    assert_eq!(m(&["--everything", "/g"]), Mode::Everything);
    assert_eq!(m(&["--no-clinical-data", "/g"]), Mode::NoClinical);
    assert_eq!(m(&["--refresh", "clinvar", "/g"]), Mode::Refresh);
    assert_eq!(m(&["--og"]), Mode::Og);
    assert_eq!(
        m(&["--clinical-data", "/g"]),
        Mode::OptIn("clinical-data".into())
    );
    assert_eq!(
        og_parse(&["--og", "--kir-data"], &[]).unwrap().mode,
        Mode::Og
    );
}

/// run_in's argv for a fake engine that writes its arguments, one per line.
/// Paths go to bash as arguments, never inside the command text.
fn bash_run_in(genome: &str, downloads: &str, args: &[String]) -> Vec<String> {
    let d = scratch("runin-bash");
    let fake = d.join("engine");
    let rec = d.join("argv");
    fs::write(
        &fake,
        "#!/usr/bin/env bash\nprintf '%s\\n' \"$@\" > \"$(dirname \"$0\")/argv\"\n",
    )
    .unwrap();
    fs::set_permissions(&fake, fs::Permissions::from_mode(0o755)).unwrap();
    let o = Command::new("bash")
        .args([
            "-c",
            ". \"$1/scripts/lib/common.sh\" && shift && run_in \"$@\"",
            "_",
        ])
        .arg(root())
        .args(args)
        .env_clear()
        .env("PATH", std::env::var("PATH").unwrap_or_default())
        .env("GENOME_DIR", genome)
        .env("OG_DOWNLOADS", downloads)
        .env("CONTAINER_ENGINE", &fake)
        .output()
        .unwrap();
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let v = fs::read_to_string(&rec)
        .unwrap()
        .lines()
        .map(str::to_owned)
        .collect();
    let _ = fs::remove_dir_all(&d);
    v
}

#[test]
fn the_engine_argv_is_run_ins_for_every_mix_setup_uses() {
    let d = scratch("runin");
    let g = d.join("genome").to_string_lossy().into_owned();
    let inside = format!("{g}/downloads");
    let outside = d.join("shared-dl").to_string_lossy().into_owned();
    fs::create_dir_all(&outside).unwrap();
    let w = format!("{g}/reference/x");
    let img = "example/tool:1.0";
    // (run_in's words, og's RunIn, the words after the opt-outs, downloads)
    let ws = |w: &[&str]| -> Vec<String> { w.iter().map(|s| (*s).to_owned()).collect() };
    let mixes: Vec<(Vec<String>, RunIn, Vec<String>, &str)> = vec![
        (
            ws(&[img, "true"]),
            RunIn::default(),
            vec![img.into(), "true".into()],
            &inside,
        ),
        (
            ws(&["--rw", &w, img, "tool", "-o", "/genome/reference/x/out"]),
            RunIn {
                rw: vec![w.clone()],
                ..Default::default()
            },
            vec![
                img.into(),
                "tool".into(),
                "-o".into(),
                "/genome/reference/x/out".into(),
            ],
            &inside,
        ),
        (
            ws(&["--net", "--rw", &w, img, "fetch", "all"]),
            RunIn {
                net: true,
                rw: vec![w.clone()],
                ..Default::default()
            },
            vec![img.into(), "fetch".into(), "all".into()],
            &inside,
        ),
        (
            ws(&["--rw", &w, "-v", "/x/bin:/pgp-bin:ro", img, "og"]),
            RunIn {
                rw: vec![w.clone()],
                ..Default::default()
            },
            vec![
                "-v".into(),
                "/x/bin:/pgp-bin:ro".into(),
                img.into(),
                "og".into(),
            ],
            &inside,
        ),
        (
            ws(&["--rw", &w, "--cpus", "2", img, "tool"]),
            RunIn {
                rw: vec![w.clone()],
                ..Default::default()
            },
            vec!["--cpus".into(), "2".into(), img.into(), "tool".into()],
            &inside,
        ),
        (
            ws(&["--net", "--rw", &w, "-v", "/l:/lock.txt:ro", img, "fetch"]),
            RunIn {
                net: true,
                rw: vec![w.clone()],
                ..Default::default()
            },
            vec![
                "-v".into(),
                "/l:/lock.txt:ro".into(),
                img.into(),
                "fetch".into(),
            ],
            &inside,
        ),
        (
            ws(&["--rw", &format!("{outside}/a"), img, "tool"]),
            RunIn {
                rw: vec![format!("{outside}/a")],
                ..Default::default()
            },
            vec![img.into(), "tool".into()],
            &outside,
        ),
        (
            ws(&["--root", img, "id"]),
            RunIn {
                root: true,
                ..Default::default()
            },
            vec![img.into(), "id".into()],
            &outside,
        ),
    ];
    for (words, r, rest, dl) in mixes {
        let want = bash_run_in(&g, dl, &words);
        let place = Place {
            genome: g.clone(),
            downloads: Some(dl.to_owned()),
            sample: None,
            uid: rustix::process::getuid().as_raw(),
            gid: rustix::process::getgid().as_raw(),
        };
        assert_eq!(
            r.argv(&place, &rest).unwrap(),
            want,
            "run_in {}",
            words.join(" ")
        );
    }
    let place = Place {
        genome: g.clone(),
        downloads: Some(inside),
        sample: None,
        uid: 1,
        gid: 1,
    };
    let bad = RunIn {
        rw: vec!["/elsewhere".into()],
        ..Default::default()
    };
    assert!(
        bad.argv(&place, &[])
            .unwrap_err()
            .contains("is outside GENOME_DIR")
    );
    let _ = fs::remove_dir_all(&d);
}

#[test]
fn pinned_tags_match_the_shell() {
    let h = "0123456789abcdef".repeat(4);
    let shapes = [
        format!("ubuntu:24.04@sha256:{h}"),
        format!("ubuntu@sha256:{h}"),
        format!("quay.io/biocontainers/eklipse:1.8--hdfd78af_1@sha256:{h}"),
        format!("registry.example:5000/team/tool@sha256:{h}"),
        format!("registry.example:5000/team/tool:2.0@sha256:{h}"),
        "broadinstitute/gatk:4.7.0.0".to_owned(),
        format!("Upper/case@sha256:{h}"),
        format!("short@sha256:{}", &h[..63]),
    ];
    for s in shapes {
        let o = Command::new("bash")
            .args([
                "-c",
                ". \"$1/scripts/lib/common.sh\" && pinned_tag \"$2\"",
                "_",
            ])
            .arg(root())
            .arg(&s)
            .output()
            .unwrap();
        let want = String::from_utf8_lossy(&o.stdout).trim_end().to_owned();
        let got = og_setup::images::pinned_tag(&s).unwrap_or_default();
        assert_eq!(got, want, "{s}");
    }
}

/// print_missing of setup.sh (read from its text while it has it).
#[test]
fn the_missing_list_is_print_missings() {
    let text = fs::read_to_string(setup_sh()).unwrap();
    if !text.contains("print_missing() {") {
        return;
    }
    let mut m = Missing::default();
    m.note("PGS scores (step 25)", "./scripts/setup.sh /g");
    m.note("a|b with a bar", "x");
    let script = "eval \"$(sed -n '/^note_missing() /p; /^print_missing() {/,/^}/p' \"$1\")\"; MISSING=(); note_missing 'PGS scores (step 25)' './scripts/setup.sh /g'; note_missing 'a|b with a bar' x; print_missing";
    let o = Command::new("bash")
        .args(["-c", script, "_"])
        .arg(setup_sh())
        .output()
        .unwrap();
    // The shell splits each entry at its first bar, so "a|b with a bar"
    // prints as "a" with the retry "b with a bar|x": the list is fed only
    // by setup's own entries, which hold no bar.
    let bash = String::from_utf8_lossy(&o.stdout).into_owned();
    let mut first = Missing::default();
    first.note("PGS scores (step 25)", "./scripts/setup.sh /g");
    assert!(
        bash.starts_with(&first.text()),
        "{bash}\n---\n{}",
        first.text()
    );
    assert_eq!(install::SETUP_PARTIAL_EXIT, 3);
    assert_eq!(m.len(), 2);
}

/// The pins of `--sample-qc-data` live in setup.sh too while its bash
/// installer is there: both must agree.
#[test]
fn the_sample_qc_pins_match_setup_sh() {
    let text = fs::read_to_string(setup_sh()).unwrap();
    let Some(url) = text
        .lines()
        .find_map(|l| l.strip_prefix("SOMALIER_SITES_URL="))
    else {
        return;
    };
    use og_setup::install::small;
    assert_eq!(url, small::SOMALIER_SITES_URL);
    let sha = text
        .lines()
        .find_map(|l| l.strip_prefix("SOMALIER_SITES_SHA256="))
        .unwrap();
    assert_eq!(sha, small::SOMALIER_SITES_SHA256);
    for (u, s) in small::VB2_PANEL_FILES {
        assert!(
            text.lines().any(|l| l
                .trim_start_matches("VB2_PANEL_FILES=\"")
                .trim_end_matches('"')
                == format!("{u} {s}")),
            "setup.sh's VB2_PANEL_FILES lacks {u} {s}"
        );
    }
}

fn ctx_with(genome: &str) -> Ctx {
    let g = genome.to_owned();
    Ctx::new(
        &root(),
        Box::new(move |k| (k == "GENOME_DIR").then(|| g.clone())),
    )
    .unwrap()
}

/// The paths of a component's sources (each part of a family), or a folder.
fn source_paths(c: &Component, ctx: &Ctx, s: &Sources) -> Vec<PathBuf> {
    let mut v = Vec::new();
    for r in c.sources(s) {
        if r.parts.is_empty() {
            v.push(ctx.path(r, None).unwrap());
        } else {
            for p in &r.parts {
                v.push(ctx.path(r, Some(p)).unwrap());
            }
        }
    }
    v
}

/// The installed paths of C that none of its sources writes: a path of a
/// source, or a file inside a source that is a folder.
fn unwritten(c: &Component, ctx: &Ctx, s: &Sources, genome: &Path) -> Vec<PathBuf> {
    let written = source_paths(c, ctx, s);
    let reference = genome.join("reference/GRCh38_no_alt_analysis_set.fasta");
    c.installed_paths(ctx, genome, &reference)
        .unwrap()
        .into_iter()
        .filter(|p| {
            !written
                .iter()
                .any(|w| w == p || (w.to_string_lossy().ends_with('/') && p.starts_with(w)))
        })
        .collect()
}

#[test]
fn every_installed_path_of_a_ported_installer_is_one_it_writes() {
    let genome = Path::new("/g");
    let ctx = ctx_with("/g");
    let s = Sources::read(&root()).unwrap();
    let mut seen = 0;
    for c in components::TABLE {
        if !c.ported() || c.kind != CKind::Installer {
            continue;
        }
        seen += 1;
        assert!(!c.installed.is_empty(), "{} has no installed paths", c.flag);
        assert!(!c.sources(&s).is_empty(), "{} has no source row", c.flag);
        let bad = unwritten(c, &ctx, &s, genome);
        assert!(
            bad.is_empty(),
            "{}: nothing it downloads writes {bad:?}",
            c.flag
        );
    }
    assert_eq!(seen, 7, "the seven small data installers");
    // The check can fail: a row naming a path nothing writes.
    let planted = Component {
        installed: &["reference/viral-${VIRAL_REF_VERSION}/viral_genome.fa"],
        ..*components::by_flag("--viral-data").unwrap()
    };
    assert_eq!(unwritten(&planted, &ctx, &s, genome).len(), 1);
}

#[test]
fn every_row_expands_and_names_real_images() {
    let ctx = ctx_with("/g");
    let env = fs::read_to_string(root().join("versions.env")).unwrap();
    let images = og_setup::fetch::versions_pins(&env);
    let reference = Path::new("/g/reference/GRCh38_no_alt_analysis_set.fasta");
    for c in components::TABLE {
        let paths = c
            .installed_paths(&ctx, Path::new("/g"), reference)
            .unwrap_or_else(|e| panic!("{}: {e}", c.flag));
        for p in paths {
            let t = p.to_string_lossy();
            assert!(!t.contains('{') && !t.contains('$'), "{}: {t}", c.flag);
        }
        for i in c.images {
            assert!(
                images.contains_key(*i),
                "{}: {i} is not in versions.env",
                c.flag
            );
        }
        assert!(c.flag.strip_prefix("--") == Some(c.key), "{}", c.flag);
    }
    let masked = components::by_flag("--masked-reference").unwrap();
    assert_eq!(
        masked
            .installed_paths(&ctx, Path::new("/g"), reference)
            .unwrap()[0],
        PathBuf::from("/g/reference/GRCh38_no_alt_analysis_set_maskedGRC.fasta")
    );
    let mk = components::by_flag("--megane-kmers").unwrap();
    assert_eq!(
        mk.installed_paths(&ctx, Path::new("/g"), reference)
            .unwrap()[0],
        PathBuf::from("/g/reference/megane/GRCh38_no_alt_analysis_set.mk")
    );
    // MEGANE_KMERS moves it, as setup.sh's megane_kmers_path: a relative
    // path is under the install folder, and the .mi sits beside the .mk.
    for (set, want) in [
        ("custom/panel.mk", "/g/custom/panel"),
        ("/elsewhere/panel.mk", "/elsewhere/panel"),
    ] {
        let ctx = Ctx::new(
            &root(),
            Box::new(move |k| match k {
                "GENOME_DIR" => Some("/g".to_owned()),
                "MEGANE_KMERS" => Some(set.to_owned()),
                _ => None,
            }),
        )
        .unwrap();
        assert_eq!(
            mk.installed_paths(&ctx, Path::new("/g"), reference)
                .unwrap(),
            vec![
                PathBuf::from(format!("{want}.mk")),
                PathBuf::from(format!("{want}.mi"))
            ],
            "MEGANE_KMERS={set}"
        );
    }
    // --annotation-scores lists the files validate-setup.sh looks for in
    // annotations/ (it listed `annotation_scores`, which nothing writes).
    let validate = fs::read_to_string(root().join("scripts/validate-setup.sh")).unwrap();
    let scores = components::by_flag("--annotation-scores").unwrap();
    let paths = scores
        .installed_paths(&ctx, Path::new("/g"), reference)
        .unwrap();
    assert_eq!(paths.len(), 12, "{paths:?}");
    for p in paths {
        let name = p
            .strip_prefix("/g/annotations")
            .unwrap_or_else(|_| panic!("{p:?}"));
        let name = name.to_string_lossy();
        let base = name.strip_suffix(".tbi").unwrap_or(&name);
        // It checks the two SpliceAI kinds in one loop over ${kind}.
        let base = if base.starts_with("spliceai_scores.") {
            base.replace(".snv.", ".${kind}.")
                .replace(".indel.", ".${kind}.")
        } else {
            base.to_owned()
        };
        assert!(
            validate.contains(&format!("${{ANNOT_DIR}}/{base}")),
            "validate-setup.sh does not look for {name}"
        );
    }
}

/// The retry command of every row is one setup takes, and runs that row:
/// its words parse, and to the row's own mode.
#[test]
fn every_retry_command_parses_to_its_row() {
    for c in components::TABLE {
        let cmd = c.retry_command("./scripts/setup.sh", "/g");
        let words: Vec<&str> = cmd.split(' ').skip(1).collect();
        let p = og_parse(&words, &[])
            .unwrap_or_else(|s| panic!("{}: `{cmd}` stops with {}: {}", c.flag, s.code, s.stderr));
        assert_eq!(p.mode.key(), c.key, "{}: `{cmd}`", c.flag);
        if c.kind == CKind::Installer {
            assert!(
                matches!(p.mode, Mode::OptIn(_) | Mode::SampleQc),
                "{}: {:?}",
                c.flag,
                p.mode
            );
        }
    }
}

/// OPT_IN is the pattern of setup.sh's opt-in case, flag for flag and in
/// its order, while setup.sh has the case; each one is an installer row.
#[test]
fn opt_in_is_setup_sh_s_case() {
    let text = fs::read_to_string(setup_sh()).unwrap();
    let lines: Vec<&str> = text.lines().collect();
    if let Some(i) = lines.iter().position(|l| l.trim() == "OPT_IN=${1#--}") {
        let pattern = lines[i - 1].trim().strip_suffix(')').unwrap();
        let case: Vec<&str> = pattern
            .split('|')
            .map(|f| f.strip_prefix("--").unwrap())
            .collect();
        assert_eq!(case, install::OPT_IN, "setup.sh's opt-in case");
    }
    for k in install::OPT_IN {
        let c = components::by_key(k).unwrap_or_else(|| panic!("--{k} has no row"));
        assert_eq!(c.kind, CKind::Installer, "--{k}");
    }
}

#[test]
fn every_ported_key_is_a_row() {
    for k in install::small::PORTED {
        assert!(components::by_key(k).is_some(), "{k}");
        assert!(install::ported(k));
    }
    assert!(!install::ported("archaic"));
    assert!(!install::ported("setup"));
}

#[test]
fn quoting_is_printf_q() {
    for w in [
        "/data/genome dir",
        "plain/path",
        "it's",
        "a$b",
        "~x",
        "x~",
        "#y",
        "x#",
        "",
        "tab\there",
    ] {
        let o = Command::new("bash")
            .args(["-c", "printf '%q' \"$1\"", "_", w])
            .output()
            .unwrap();
        assert_eq!(
            install::shell_quote(w),
            String::from_utf8_lossy(&o.stdout),
            "{w:?}"
        );
    }
    // Control characters, as setup's bash (4.4 or later) writes them: its
    // named escapes, and a C1 control as the bytes of its UTF-8 form (og
    // wrote U+0085 as \205, which $'...' reads as one byte). macOS's
    // /bin/bash 3.2 garbles the C1 case, and setup does not run on it.
    let new_bash = Command::new("bash")
        .args([
            "-c",
            "(( BASH_VERSINFO[0] * 100 + BASH_VERSINFO[1] >= 404 ))",
        ])
        .status()
        .is_ok_and(|s| s.success());
    if !new_bash {
        eprintln!("skipped the control characters: bash is older than 4.4");
        return;
    }
    for w in [
        "a\u{85}b",
        "a\u{9f}b",
        "a\x07\x08\x0b\x0c\x1bb",
        "x\x01",
        "a\x7fb",
        "it's\x01",
        "back\\slash\x01",
        "new\nline",
    ] {
        let o = Command::new("bash")
            .args(["-c", "printf '%q' \"$1\"", "_", w])
            .env("LC_ALL", "C.UTF-8")
            .output()
            .unwrap();
        assert_eq!(
            install::shell_quote(w),
            String::from_utf8_lossy(&o.stdout),
            "{w:?}"
        );
    }
}

/// Only the usage template's own `$0` and `${NAME}` are expanded: a value
/// holding a placeholder is written as it is (one that held its own name
/// was expanded forever).
#[test]
fn usage_expands_only_its_template() {
    let (tx, rx) = std::sync::mpsc::channel();
    std::thread::spawn(move || {
        let _ = tx.send(install::usage("/s/${BWA_VERSION}/setup.sh", &|n| {
            format!("${{{n}}}")
        }));
    });
    let text = rx
        .recv_timeout(std::time::Duration::from_secs(20))
        .expect("usage did not return");
    assert!(text.contains("/s/${BWA_VERSION}/setup.sh"), "{text}");
    let plain = install::usage("./scripts/setup.sh", &|_| "v1".to_owned());
    assert!(!plain.contains("${") && !plain.contains("$0"), "{plain}");
    assert!(
        plain.contains("./scripts/setup.sh /data/genomics"),
        "{plain}"
    );
}

/// A variable that is not UTF-8 (a name or a value) does not stop setup;
/// og's environment leaves it out (it panicked). This test runs itself
/// as og, with such a variable.
#[test]
fn a_variable_that_is_not_utf8_is_left_out() {
    use std::os::unix::ffi::OsStrExt;
    const CHILD: &str = "OG_SETUP_TEST_VARS";
    if std::env::var_os(CHILD).is_some() {
        let v = install::Vars::from_process();
        assert_eq!(v.get(CHILD), Some("1"));
        assert_eq!(v.get("OG_SETUP_TEST_BAD"), None);
        return;
    }
    let o = Command::new(std::env::current_exe().unwrap())
        .args([
            "--exact",
            "a_variable_that_is_not_utf8_is_left_out",
            "--test-threads=1",
        ])
        .env(CHILD, "1")
        .env(
            "OG_SETUP_TEST_BAD",
            std::ffi::OsStr::from_bytes(b"\xff\xfe"),
        )
        .env(std::ffi::OsStr::from_bytes(b"OG_\xff"), "x")
        .output()
        .unwrap();
    assert!(
        o.status.success(),
        "{}{}",
        String::from_utf8_lossy(&o.stdout),
        String::from_utf8_lossy(&o.stderr)
    );
    assert!(String::from_utf8_lossy(&o.stdout).contains("1 passed"));
}

/// A flag og hands back runs setup.sh in og's place (exec): setup.sh's
/// bash has og's pid, so a SIGTERM sent to og ends it, and nothing keeps
/// installing after og is gone. This test runs itself as og, with a
/// stand-in checkout whose setup.sh prints its pid and waits.
#[test]
fn a_flag_handed_back_runs_in_og_s_place() {
    use std::io::BufRead;
    const CHILD: &str = "OG_SETUP_TEST_HANDBACK";
    if std::env::var_os(CHILD).is_some() {
        let args = ["--archaic".to_owned(), "/g".to_owned()];
        let code = install::cli(&args, &mut std::io::stdout(), &mut std::io::stderr());
        println!("og went on after the hand-back (exit {code})");
        std::process::exit(i32::from(code));
    }
    let d = scratch("handback-exec");
    let sh = d.join("scripts/setup.sh");
    fs::create_dir_all(sh.parent().unwrap()).unwrap();
    fs::write(
        &sh,
        "#!/usr/bin/env bash\necho \"setup.sh pid $$\"\nexec sleep 60\n",
    )
    .unwrap();
    fs::set_permissions(&sh, fs::Permissions::from_mode(0o755)).unwrap();
    let mut og = Command::new(std::env::current_exe().unwrap())
        .args([
            "--exact",
            "a_flag_handed_back_runs_in_og_s_place",
            "--test-threads=1",
            "--nocapture",
        ])
        .env(CHILD, "1")
        .env("PGP_ROOT", &d)
        .env_remove("OG_SETUP")
        .env_remove("OG_SETUP_SCRIPT")
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::null())
        .spawn()
        .unwrap();
    let mut said = String::new();
    let mut lines = std::io::BufReader::new(og.stdout.take().unwrap()).lines();
    let pid = loop {
        let Some(Ok(l)) = lines.next() else {
            break None;
        };
        said.push_str(&l);
        said.push('\n');
        if let Some(p) = l.split("setup.sh pid ").nth(1) {
            break p.trim().parse::<i32>().ok();
        }
    };
    let og_pid = i32::try_from(og.id()).unwrap();
    let _ = rustix::process::kill_process(
        rustix::process::Pid::from_raw(og_pid).unwrap(),
        rustix::process::Signal::TERM,
    );
    let status = og.wait().unwrap();
    // A setup.sh that ran as og's child is still there: stop it.
    let left = pid
        .filter(|&p| p != og_pid)
        .and_then(rustix::process::Pid::from_raw);
    if let Some(p) = left {
        let _ = rustix::process::kill_process(p, rustix::process::Signal::KILL);
    }
    let _ = fs::remove_dir_all(&d);
    assert!(pid.is_some(), "setup.sh never ran: {said}");
    assert_eq!(
        pid,
        Some(og_pid),
        "setup.sh ran as og's child, not in its place: {said}"
    );
    assert_eq!(
        std::os::unix::process::ExitStatusExt::signal(&status),
        Some(15),
        "{status}"
    );
}

/// Every process og's setup starts goes through og_setup::child, so a
/// signal stops it before a lock goes (lock.rs).
#[test]
fn every_child_is_started_through_child_rs() {
    let found = |text: &str| -> usize {
        let code = text.split("#[cfg(test)]").next().unwrap_or(text);
        code.lines()
            .filter(|l| {
                !l.trim_start().starts_with("//")
                    && [".spawn()", ".status()", ".output()"]
                        .iter()
                        .any(|c| l.contains(c))
            })
            .count()
    };
    // The check can fail.
    assert_eq!(found("let o = Command::new(\"du\").output();\n"), 1);
    let mut bad = vec![];
    let mut dirs = vec![root().join("rust/crates/og-setup/src")];
    while let Some(d) = dirs.pop() {
        for e in fs::read_dir(&d).unwrap() {
            let p = e.unwrap().path();
            if p.is_dir() {
                dirs.push(p);
            } else if p.extension().is_some_and(|x| x == "rs")
                && p.file_name().is_some_and(|n| n != "child.rs")
                && found(&fs::read_to_string(&p).unwrap()) > 0
            {
                bad.push(p);
            }
        }
    }
    assert!(
        bad.is_empty(),
        "start these through og_setup::child: {bad:?}"
    );
}
