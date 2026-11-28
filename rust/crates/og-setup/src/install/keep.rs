//! The raw downloads an installer keeps, as scripts/lib/common.sh handles
//! them: [`adopt`] moves a download an older setup left beside its install
//! into the downloads folder (OG_DOWNLOADS), and [`done_with`] keeps an
//! installed download there, or deletes it under SETUP_KEEP_DOWNLOADS=false.
//!
//! Each moves and deletes files as `mv` and `rm -rf` do. When the native
//! call fails (another filesystem, no permission), the host's own `mv`,
//! `rm` or `mkdir` runs it, so the message on stderr is the one bash prints.

use super::Run;
use std::io::Write;
use std::path::Path;
use std::process::{Command, Stdio};

/// `[ -e PATH ]`: there, following a symlink.
pub fn exists(p: &str) -> bool {
    Path::new(p).exists()
}

/// `[ -f PATH ]`: a regular file, following a symlink.
pub fn is_file(p: &str) -> bool {
    Path::new(p).is_file()
}

fn readable(p: &str) -> bool {
    rustix::fs::access(p, rustix::fs::Access::READ_OK).is_ok()
}

/// Runs a host tool with stdout and stderr shared; true when it exits 0.
fn host(out: &mut dyn Write, tool: &str, args: &[&str]) -> bool {
    let _ = out.flush();
    crate::child::status(Command::new(tool).args(args).stdin(Stdio::null()))
        .is_ok_and(|s| s.success())
}

/// `mv -- FROM TO`.
pub fn mv(from: &str, to: &str, out: &mut dyn Write) -> bool {
    std::fs::rename(from, to).is_ok() || host(out, "mv", &["--", from, to])
}

/// `rm -rf -- PATH`: true when PATH is gone.
pub fn rm_rf(p: &str, out: &mut dyn Write) -> bool {
    let path = Path::new(p);
    let gone = match std::fs::symlink_metadata(path) {
        Err(e) => e.kind() == std::io::ErrorKind::NotFound,
        Ok(m) if m.is_dir() => std::fs::remove_dir_all(path).is_ok(),
        Ok(_) => std::fs::remove_file(path).is_ok(),
    };
    gone || host(out, "rm", &["-rf", "--", p])
}

/// `mkdir -p DIR`.
pub fn mkdir_p(dir: &str, out: &mut dyn Write) -> bool {
    std::fs::create_dir_all(dir).is_ok() || host(out, "mkdir", &["-p", "--", dir])
}

/// What adopt checks a file with: `fetch URL FILE [KIND VALUE]`, run as
/// `og setup fetch` with its output discarded (`> /dev/null 2>&1`).
pub struct Check<'a> {
    pub url: &'a str,
    /// The checksum kind and value; `None` for adopt's URL-only form.
    pub sum: Option<(&'a str, &'a str)>,
}

/// adopt OLD NEW [URL [KIND VALUE]] of common.sh: a raw download (and its
/// `.part`) an older setup left at OLD moves to NEW in the downloads
/// folder. When both are there, one goes only once fetch's check of NEW
/// says which is the download. Messages go to OUT.
pub fn adopt(r: &Run, old: &str, new: &str, check: Option<&Check>, out: &mut dyn Write) {
    let verify = |file: &str| -> bool {
        let Some(c) = check else { return false };
        let mut cmd = Command::new(&r.og);
        cmd.args(["setup", "fetch", c.url, file]);
        if let Some((kind, value)) = c.sum {
            cmd.args([kind, value]);
        }
        for k in [
            "PGP_ROOT",
            "GENOME_DIR",
            "OG_DOWNLOADS",
            "SETUP_KEEP_DOWNLOADS",
        ] {
            if let Some(v) = r.vars.get(k) {
                cmd.env(k, v);
            }
        }
        cmd.stdin(Stdio::null())
            .stdout(Stdio::null())
            .stderr(Stdio::null());
        crate::child::status(&mut cmd).is_ok_and(|s| s.success())
    };
    let has_url = check.is_some();
    let has_sum = check.is_some_and(|c| c.sum.is_some());
    adopt_with(old, new, has_url, has_sum, &verify, out);
}

/// adopt with the check given: VERIFY(NEW) is `fetch URL NEW ...`, HAS_URL
/// whether a URL was given at all (`$# -gt 0`), HAS_SUM whether a checksum
/// came with it (`$# -ge 3`).
pub fn adopt_with(
    old: &str,
    new: &str,
    has_url: bool,
    has_sum: bool,
    verify: &dyn Fn(&str) -> bool,
    out: &mut dyn Write,
) {
    let old_part = format!("{old}.part");
    if exists(old) && exists(new) {
        // Only files are checked: fetch downloads into a NEW that is no file.
        let has_url = has_url && is_file(old) && is_file(new);
        if has_url && verify(new) {
            // Not kept: a second copy of a checked file in the downloads folder.
            rm_rf(old, out);
            rm_rf(&old_part, out);
            let _ = writeln!(
                out,
                "  Deleted {old}, left by an older setup: {new} is the same download, checked."
            );
            return;
        }
        let aside = format!("{new}.adopt");
        if has_url && mv(new, &aside, out) {
            if mv(old, new, out) {
                if verify(new) {
                    // Not kept: the copy that failed its check, replaced by the one that passed.
                    rm_rf(&aside, out);
                    let _ = writeln!(
                        out,
                        "  {new} failed its check; replaced it with {old}, left by an older setup, which passes."
                    );
                    return;
                }
                if mv(new, old, out) {
                    mv(&aside, new, out);
                }
            } else {
                mv(&aside, new, out);
            }
        }
        let _ = writeln!(
            out,
            "  [WARN] {old}, left by an older setup, and {new} are both here; neither is deleted, as no check says which is the download."
        );
        return;
    }
    for x in ["", ".part"] {
        let ox = format!("{old}{x}");
        let nx = format!("{new}{x}");
        if !exists(&ox) {
            continue;
        }
        if x == ".part" && exists(new) {
            // A whole NEW is there, which fetch checks. A partial copy beside
            // it would be an orphan: it goes once NEW passes the check, and
            // with no check, or one that fails, it stays where it is.
            if has_sum && is_file(new) && verify(new) {
                // Not kept: a partial copy of a checked file in the downloads folder.
                rm_rf(&old_part, out);
                let _ = writeln!(
                    out,
                    "  Deleted {old}.part, a partial download left by an older setup: {new} is the whole download, checked."
                );
            } else {
                let _ = writeln!(
                    out,
                    "  [WARN] {old}.part, a partial download left by an older setup, stays: {new} is already there and fetch checks it. Once setup has finished, {old}.part is safe to delete."
                );
            }
        } else if exists(&nx) {
            // Two partial downloads: neither is deleted unchecked. The larger
            // one moves to NEW.part, so fetch resumes the most progress.
            let size = |p: &str| std::fs::metadata(p).map(|m| m.len()).unwrap_or(0);
            let aside = format!("{nx}.adopt");
            let swapped = x == ".part"
                && is_file(&ox)
                && is_file(&nx)
                && readable(&ox)
                && readable(&nx)
                && size(&ox) > size(&nx)
                && mv(&nx, &aside, out)
                && (mv(&ox, &nx, out) || {
                    mv(&aside, &nx, out);
                    false
                })
                && mv(&aside, &ox, out);
            if swapped {
                let _ = writeln!(
                    out,
                    "  [WARN] {old}.part, left by an older setup, was the larger partial download: it resumes at {new}.part, and the smaller one is kept at {old}.part."
                );
            } else {
                let _ = writeln!(
                    out,
                    "  [WARN] {ox}, left by an older setup, and {nx} are both here; neither is deleted, and {nx} is the one fetch uses."
                );
            }
        } else {
            let parent = Path::new(new)
                .parent()
                .map(|p| p.to_string_lossy().into_owned())
                .filter(|p| !p.is_empty())
                .unwrap_or_else(|| ".".to_owned());
            if mkdir_p(&parent, out) && mv(&ox, &nx, out) {
                let _ = writeln!(out, "  Moved {ox}, left by an older setup, to {nx}.");
            }
        }
    }
}

/// done_with FILE... of common.sh: each installed raw download stays in the
/// downloads folder, unless SETUP_KEEP_DOWNLOADS is false; then it goes,
/// with its `.part` and the checksum saved beside it. The status is the
/// last file's, as the shell loop's is.
pub fn done_with(r: &Run, files: &[&str], out: &mut dyn Write) -> bool {
    if r.vars.nonempty("SETUP_KEEP_DOWNLOADS").unwrap_or("true") != "false" {
        return true;
    }
    let mut ok = true;
    for f in files {
        if !exists(f) {
            ok = true;
            continue;
        }
        ok = rm_rf(f, out)
            && rm_rf(&format!("{f}.part"), out)
            && rm_rf(&format!("{f}.checksum"), out);
        if ok {
            let _ = writeln!(
                out,
                "  Deleted the raw download {f} (SETUP_KEEP_DOWNLOADS=false)."
            );
        }
    }
    ok
}
