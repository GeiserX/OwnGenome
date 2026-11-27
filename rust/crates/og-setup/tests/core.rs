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

/// A second flag is refused; setup.sh would read it as the install folder.
#[test]
fn a_second_flag_is_refused() {
    for args in [
        vec!["--kir-data", "--viral-data"],
        vec!["--kir-data", "/g", "--viral-data"],
        vec!["/g", "--kir-data"],
        vec!["--refresh", "clinvar", "--everything"],
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
        format!(
            "#!/usr/bin/env bash\ncase \"$1 $2\" in\n  'setup ported'|'setup plan') exit 0 ;;\n  'setup install') shift 2; {{ printf '%s\\n' \"keep=$SETUP_KEEP_DOWNLOADS\"; printf '%s\\n' \"$@\"; }} > '{}'; exit 0 ;;\nesac\nexit 3\n",
            rec.display()
        ),
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
fn bash_run_in(genome: &str, downloads: &str, args: &str) -> Vec<String> {
    let d = scratch("runin-bash");
    let fake = d.join("engine");
    let rec = d.join("argv");
    fs::write(
        &fake,
        format!(
            "#!/usr/bin/env bash\nprintf '%s\\n' \"$@\" > '{}'\n",
            rec.display()
        ),
    )
    .unwrap();
    fs::set_permissions(&fake, fs::Permissions::from_mode(0o755)).unwrap();
    let script = format!(
        ". '{}/scripts/lib/common.sh' && run_in {args}",
        root().display()
    );
    let o = Command::new("bash")
        .args(["-c", &script])
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
    let mixes: Vec<(String, RunIn, Vec<String>, &str)> = vec![
        (
            format!("{img} true"),
            RunIn::default(),
            vec![img.into(), "true".into()],
            &inside,
        ),
        (
            format!("--rw {w} {img} tool -o /genome/reference/x/out"),
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
            format!("--net --rw {w} {img} fetch all"),
            RunIn {
                net: true,
                rw: vec![w.clone()],
                ..Default::default()
            },
            vec![img.into(), "fetch".into(), "all".into()],
            &inside,
        ),
        (
            format!("--rw {w} -v /x/bin:/pgp-bin:ro {img} og"),
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
            format!("--rw {w} --cpus 2 {img} tool"),
            RunIn {
                rw: vec![w.clone()],
                ..Default::default()
            },
            vec!["--cpus".into(), "2".into(), img.into(), "tool".into()],
            &inside,
        ),
        (
            format!("--net --rw {w} -v /l:/lock.txt:ro {img} fetch"),
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
            format!("--rw {outside}/a {img} tool"),
            RunIn {
                rw: vec![format!("{outside}/a")],
                ..Default::default()
            },
            vec![img.into(), "tool".into()],
            &outside,
        ),
        (
            format!("--root {img} id"),
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
        assert_eq!(r.argv(&place, &rest).unwrap(), want, "run_in {words}");
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
                &format!(
                    ". '{}/scripts/lib/common.sh' && pinned_tag \"$1\"",
                    root().display()
                ),
                "_",
                &s,
            ])
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
    let script = format!(
        "eval \"$(sed -n '/^note_missing() /p; /^print_missing() {{/,/^}}/p' '{}')\"; MISSING=(); note_missing 'PGS scores (step 25)' './scripts/setup.sh /g'; note_missing 'a|b with a bar' x; print_missing",
        setup_sh().display()
    );
    let o = Command::new("bash").args(["-c", &script]).output().unwrap();
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
}
