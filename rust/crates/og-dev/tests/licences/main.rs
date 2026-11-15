//! `og dev licences notices` on a planted workspace, through real
//! `cargo metadata`. The workspace links crates the binary must credit
//! (one per kind of licence text) and crates it must leave out (build, dev,
//! procedural macro, another target, a feature that is off). Its Windows
//! build links a crate outside the allow list and one with no licence text,
//! and must fail naming both.
use og_dev::licences;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};

static NEXT: AtomicUsize = AtomicUsize::new(0);

struct Temp(PathBuf);

impl Temp {
    fn new() -> Self {
        let p = std::env::temp_dir().join(format!(
            "og-licences-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir_all(&p).unwrap();
        Self(p)
    }
    fn write(&self, rel: &str, text: &str) {
        let p = self.0.join(rel);
        std::fs::create_dir_all(p.parent().unwrap()).unwrap();
        std::fs::write(p, text).unwrap();
    }
    fn path(&self, rel: &str) -> String {
        self.0.join(rel).to_str().unwrap().into()
    }
}

impl Drop for Temp {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

const GOOD_TEXT: &str = "Copyright (c) the good authors\nshared permissive text\n";

/// A crate at vendor/NAME: its manifest lines after [package], its lib, and
/// any files (path relative to the crate, text).
fn vendor(t: &Temp, name: &str, manifest: &str, files: &[(&str, &str)]) {
    t.write(
        &format!("vendor/{name}/Cargo.toml"),
        &format!("[package]\nname = \"{name}\"\nversion = \"0.1.0\"\nedition = \"2021\"\n{manifest}\n"),
    );
    t.write(&format!("vendor/{name}/src/lib.rs"), "");
    for (rel, text) in files {
        t.write(&format!("vendor/{name}/{rel}"), text);
    }
}

/// The planted workspace: ws/app (the binary) and its dependencies in
/// vendor/, outside the workspace folder so they are not members.
fn workspace() -> Temp {
    let t = Temp::new();
    t.write("ws/Cargo.toml", "[workspace]\nmembers = [\"app\"]\nresolver = \"2\"\n");
    t.write(
        "ws/app/Cargo.toml",
        r#"[package]
name = "app"
version = "0.1.0"
edition = "2021"
license = "GPL-3.0-or-later"

[features]
default = ["extra"]
extra = ["dep:extra"]

[dependencies]
good = { path = "../../vendor/good" }
dual = { path = "../../vendor/dual" }
csys = { path = "../../vendor/csys" }
pm = { path = "../../vendor/pm" }
extra = { path = "../../vendor/extra", optional = true }

[build-dependencies]
buildonly = { path = "../../vendor/buildonly" }

[dev-dependencies]
devonly = { path = "../../vendor/devonly" }

[target.'cfg(windows)'.dependencies]
bad = { path = "../../vendor/bad" }
notext = { path = "../../vendor/notext" }
"#,
    );
    t.write("ws/app/src/lib.rs", "");
    vendor(&t, "good", "license = \"MIT\"", &[("LICENSE", GOOD_TEXT)]);
    // The same text as good: it must come once, for both.
    vendor(&t, "extra", "license = \"MIT OR Apache-2.0\"", &[("LICENSE-MIT", GOOD_TEXT)]);
    vendor(
        &t,
        "dual",
        "license = \"MIT/Apache-2.0\"\nauthors = [\"Dual Author <dual@example.org>\"]",
        &[],
    );
    vendor(
        &t,
        "csys",
        "license = \"MIT\"\nlinks = \"csys\"",
        &[
            ("LICENSE-MIT", "csys crate text\n"),
            // cargo wants a build script beside `links`.
            ("build.rs", "fn main() {}\n"),
            ("vendored/clib/COPYING", "the C library's own notice\n"),
            ("vendored/clib/src/clib.c", "int clib(void) { return 0; }\n"),
        ],
    );
    let sspl = [("LICENSE", "server side text\n")];
    vendor(&t, "pm", "license = \"SSPL-1.0\"\n[lib]\nproc-macro = true\n[dependencies]\npmdep = { path = \"../pmdep\" }", &sspl);
    vendor(&t, "pmdep", "license = \"SSPL-1.0\"", &sspl);
    vendor(&t, "buildonly", "license = \"SSPL-1.0\"", &sspl);
    vendor(&t, "devonly", "license = \"SSPL-1.0\"", &sspl);
    vendor(&t, "bad", "license = \"SSPL-1.0\"", &sspl);
    vendor(&t, "notext", "license = \"Apache-2.0\"", &[]);
    let cargo = std::env::var_os("CARGO").unwrap_or_else(|| "cargo".into());
    let lock = std::process::Command::new(cargo)
        .args(["generate-lockfile", "--offline", "--manifest-path"])
        .arg(Path::new(&t.path("ws/Cargo.toml")))
        .output()
        .unwrap();
    assert!(lock.status.success(), "{}", String::from_utf8_lossy(&lock.stderr));
    t
}

fn notices(t: &Temp, targets: &[&str]) -> (u8, String, String) {
    let out_file = t.path("THIRD-PARTY-LICENSES.txt");
    let _ = std::fs::remove_file(&out_file);
    let mut args: Vec<String> = vec!["notices".into(), "--manifest-path".into(), t.path("ws/Cargo.toml")];
    args.extend(targets.iter().map(|s| s.to_string()));
    args.extend(["--package".into(), "app".into(), "--out".into(), out_file.clone()]);
    let mut out = Vec::new();
    let mut err = Vec::new();
    let code = licences::cli(&args, &mut out, &mut err);
    let text = std::fs::read_to_string(&out_file).unwrap_or_default();
    (code, text, String::from_utf8(err).unwrap())
}

const LINUX: &str = "x86_64-unknown-linux-gnu";

#[test]
fn every_linked_crate_is_credited_and_nothing_else() {
    let t = workspace();
    let (code, text, err) = notices(&t, &["--target", LINUX]);
    assert_eq!(code, 0, "{err}");
    let listed: Vec<&str> = text
        .lines()
        .skip_while(|l| !l.starts_with("Crates ("))
        .skip(2)
        .take_while(|l| !l.is_empty())
        .collect();
    assert_eq!(
        listed,
        [
            "csys 0.1.0  MIT",
            "dual 0.1.0  MIT/Apache-2.0",
            "extra 0.1.0  MIT OR Apache-2.0",
            "good 0.1.0  MIT",
        ],
        "{text}"
    );
    // The shared text once, for both crates that ship it.
    assert_eq!(text.matches("shared permissive text").count(), 1, "{text}");
    assert!(text.contains("For extra 0.1.0: LICENSE-MIT\nFor good 0.1.0: LICENSE\n"), "{text}");
    // A crate that builds C carries its library's notice from below its root.
    assert!(text.contains("For csys 0.1.0: vendored/clib/COPYING\n\nthe C library's own notice"), "{text}");
    // No file, MIT allowed: the MIT text with the crate's authors.
    assert!(text.contains("For dual 0.1.0: MIT, no licence file in the crate\n\nCopyright (c) Dual Author <dual@example.org>\n\nPermission is hereby granted"), "{text}");
    for absent in ["app 0.1.0", "pm 0.1.0", "pmdep", "buildonly", "devonly", "bad 0.1.0", "notext", "server side text"] {
        assert!(!text.contains(absent), "{absent} must not be credited:\n{text}");
    }
}

#[test]
fn a_build_without_default_features_leaves_the_optional_crate_out() {
    let t = workspace();
    let (code, text, err) = notices(&t, &["--target", LINUX, "--target-no-default", "aarch64-apple-darwin"]);
    assert_eq!(code, 0, "{err}");
    assert!(text.contains("Targets: aarch64-apple-darwin, x86_64-unknown-linux-gnu\n"), "{text}");
    assert!(text.contains("extra 0.1.0  MIT OR Apache-2.0  (only x86_64-unknown-linux-gnu)\n"), "{text}");
    assert!(text.contains("good 0.1.0  MIT\n"), "{text}");
}

#[test]
fn a_crate_outside_the_allow_list_or_without_a_text_fails() {
    let t = workspace();
    let (code, text, err) = notices(&t, &["--target", LINUX, "--target", "x86_64-pc-windows-gnu"]);
    assert_eq!(code, 1, "{err}");
    assert!(text.is_empty(), "no file on an error, got:\n{text}");
    assert!(err.contains("bad 0.1.0 is licensed \"SSPL-1.0\", which the allow list"), "{err}");
    assert!(err.contains("notext 0.1.0 ships no licence file, and \"Apache-2.0\" cannot be met with the MIT text"), "{err}");
    assert_eq!(err.lines().count(), 2, "{err}");
}

#[test]
fn usage_errors_exit_2() {
    let mut out = Vec::new();
    for args in [&[][..], &["notices"][..], &["notices", "--target"][..], &["assets"][..]] {
        let mut err = Vec::new();
        let args: Vec<String> = args.iter().map(|s| s.to_string()).collect();
        assert_eq!(licences::cli(&args, &mut out, &mut err), 2, "{args:?}");
        assert!(String::from_utf8(err).unwrap().contains("usage: og dev licences notices"));
    }
}
