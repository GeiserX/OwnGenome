//! Test support for OwnGenome's browser checks: a headless Chrome driven
//! over the DevTools protocol, a static file server for the pages and the
//! scan of the live site's assets. No product code depends on it.
//!
//! - [`find_chrome`] and [`chrome_or_skip`]: which Chrome to run, and what a
//!   test does without one (skip on a laptop, fail in CI).
//! - [`Browser`]: one headless Chrome and one page: [`Browser::go`],
//!   [`Browser::eval`], [`Browser::wait`] and [`Browser::page_errors`].
//! - [`StaticServer`]: a directory on 127.0.0.1 with the headers a test
//!   asks for (a Content-Security-Policy).
//! - [`asset_problems`]: outside URLs, inline scripts, `on*` attributes and
//!   `eval` in a site's files.
//!
//! How the checks run: `tests/browser/README.md`.

mod assets;
mod cdp;
mod serve;

pub use assets::{asset_problems, read_assets};
pub use cdp::{Browser, js_str};
pub use pycompat::json::Value;
pub use serve::StaticServer;

use std::ffi::OsStr;
use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};

/// The names tried on `PATH` when `CHROME` is not set, then the macOS app.
const NAMES: [&str; 4] = [
    "google-chrome",
    "google-chrome-stable",
    "chromium",
    "chromium-browser",
];
const MAC_APP: &str = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome";

/// The Chrome to run: `CHROME` when it is set, else the first of
/// `google-chrome`, `google-chrome-stable`, `chromium`, `chromium-browser`
/// on `PATH`, else Google Chrome's macOS app. `None` when there is none.
pub fn find_chrome() -> Option<PathBuf> {
    find_chrome_in(
        std::env::var_os("CHROME").as_deref(),
        std::env::var_os("PATH").as_deref(),
    )
}

/// [`find_chrome`] with `CHROME` and `PATH` given. A `CHROME` that is set
/// but names no executable gives `None`: a wrong path fails loudly instead
/// of falling back to another Chrome.
pub fn find_chrome_in(chrome: Option<&OsStr>, path: Option<&OsStr>) -> Option<PathBuf> {
    if let Some(c) = chrome.filter(|c| !c.is_empty()) {
        return resolve(Path::new(c), path);
    }
    NAMES
        .iter()
        .find_map(|n| resolve(Path::new(n), path))
        .or_else(|| resolve(Path::new(MAC_APP), None))
}

/// A path with a slash as it is, a bare name looked up on `path`; only an
/// executable file counts.
fn resolve(name: &Path, path: Option<&OsStr>) -> Option<PathBuf> {
    let exe = |p: &Path| {
        p.metadata()
            .is_ok_and(|m| m.is_file() && m.permissions().mode() & 0o111 != 0)
    };
    if name.components().count() > 1 {
        return exe(name).then(|| name.to_path_buf());
    }
    std::env::split_paths(path?)
        .map(|d| d.join(name))
        .find(|p| exe(p))
}

/// The Chrome for a test of `what`, or `None` after a `[SKIP]` line when
/// there is none. In CI (`CI=true`) a missing Chrome panics, so a runner
/// without one can never pass the browser checks by skipping them.
pub fn chrome_or_skip(what: &str) -> Option<PathBuf> {
    let ci = std::env::var("CI").is_ok_and(|v| v == "true");
    match without_chrome(find_chrome(), ci, what) {
        Ok(c) => c,
        Err(e) => panic!("{e}"),
    }
}

/// The decision of [`chrome_or_skip`]: the Chrome found, `Ok(None)` after a
/// `[SKIP]` line, or the failure in CI.
fn without_chrome(found: Option<PathBuf>, ci: bool, what: &str) -> Result<Option<PathBuf>, String> {
    match found {
        Some(c) => Ok(Some(c)),
        None if ci => Err(format!(
            "[FAIL] no Chrome or Chromium to {what}: CI must run the browser checks (set CHROME)"
        )),
        None => {
            println!("[SKIP] no Chrome or Chromium to {what}");
            Ok(None)
        }
    }
}

/// The repository's root: three levels above this crate.
pub fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .ancestors()
        .nth(3)
        .expect("rust/crates/og-browser sits three levels below the root")
        .to_path_buf()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    /// A fresh scratch folder for one test.
    fn scratch(name: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!("og-browser-{}-{name}", std::process::id()));
        let _ = fs::remove_dir_all(&d);
        fs::create_dir_all(&d).unwrap();
        d
    }

    fn executable(p: &Path, mode: u32) {
        fs::write(p, "#!/bin/sh\n").unwrap();
        fs::set_permissions(p, fs::Permissions::from_mode(mode)).unwrap();
    }

    #[test]
    fn chrome_names_an_executable_or_nothing() {
        let d = scratch("chrome-env");
        let c = d.join("my-chrome");
        executable(&c, 0o755);
        assert_eq!(find_chrome_in(Some(c.as_os_str()), None), Some(c.clone()));
        // Set but wrong: no fallback to a Chrome on PATH.
        let bin = d.join("bin");
        fs::create_dir(&bin).unwrap();
        executable(&bin.join("chromium"), 0o755);
        let missing = d.join("nothing-here");
        assert_eq!(
            find_chrome_in(Some(missing.as_os_str()), Some(bin.as_os_str())),
            None
        );
        // Not executable.
        let plain = d.join("plain");
        executable(&plain, 0o644);
        assert_eq!(find_chrome_in(Some(plain.as_os_str()), None), None);
        fs::remove_dir_all(&d).unwrap();
    }

    #[test]
    fn without_chrome_set_a_name_on_path_is_found() {
        let d = scratch("path");
        let (a, b) = (d.join("a"), d.join("b"));
        fs::create_dir(&a).unwrap();
        fs::create_dir(&b).unwrap();
        executable(&b.join("chromium"), 0o755);
        executable(&a.join("google-chrome-stable"), 0o644); // not executable: skipped
        let path = std::env::join_paths([&a, &b]).unwrap();
        assert_eq!(find_chrome_in(None, Some(&path)), Some(b.join("chromium")));
        assert_eq!(
            find_chrome_in(Some(OsStr::new("")), Some(&path)),
            Some(b.join("chromium"))
        );
        // Nothing on PATH: only the macOS app can still be found.
        let found = find_chrome_in(None, Some(d.join("empty").as_os_str()));
        assert!(
            found.is_none() || found == Some(PathBuf::from(MAC_APP)),
            "{found:?}"
        );
        fs::remove_dir_all(&d).unwrap();
    }

    #[test]
    fn a_missing_chrome_skips_locally_and_fails_in_ci() {
        assert_eq!(without_chrome(None, false, "load the page"), Ok(None));
        let e = without_chrome(None, true, "load the page").unwrap_err();
        assert!(e.contains("no Chrome or Chromium to load the page"), "{e}");
        let c = PathBuf::from("/x/chrome");
        assert_eq!(without_chrome(Some(c.clone()), true, "x"), Ok(Some(c)));
    }

    #[test]
    fn the_root_holds_the_live_assets() {
        assert!(repo_root().join("assets/live/live.js").is_file());
    }
}
