//! setup's images, as scripts/lib/common.sh and scripts/setup.sh handle
//! them: the list setup pulls ([`pipeline_images`]), a pull tried
//! FETCH_TRIES times ([`pull_image`]), the local tag a digest-pinned image
//! gets so a dangling-image cleanup keeps it ([`pinned_tag`],
//! [`keep_image`]), and the pull of every image ([`pull_images`],
//! `setup.sh --pull-only`) or of one opt-in step's ([`ensure_image`]).
//!
//! Every engine call is the one the shell made, with the same arguments, so
//! a fake engine logs the same lines for both, and Docker's own messages go
//! to the terminal as they did.

use crate::install::Vars;
use std::io::Write;
use std::path::Path;
use std::process::{Command, Stdio};

/// The container engine: CONTAINER_ENGINE, docker by default.
#[derive(Debug, Clone)]
pub struct Engine {
    pub cmd: String,
}

impl Engine {
    pub fn from_vars(v: &Vars) -> Engine {
        Engine {
            cmd: v
                .nonempty("CONTAINER_ENGINE")
                .unwrap_or("docker")
                .to_owned(),
        }
    }

    /// The call's exit status with its output discarded (`&>/dev/null`).
    pub fn quiet(&self, args: &[&str]) -> bool {
        let _ = std::io::stdout().flush();
        Command::new(&self.cmd)
            .args(args)
            .stdin(Stdio::null())
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .status()
            .is_ok_and(|s| s.success())
    }

    /// The call's exit status, its output on the terminal.
    pub fn loud(&self, args: &[&str]) -> bool {
        let _ = std::io::stdout().flush();
        Command::new(&self.cmd)
            .args(args)
            .stdin(Stdio::null())
            .status()
            .is_ok_and(|s| s.success())
    }

    /// `$(engine ARGS 2>/dev/null)`: stdout without its trailing newlines,
    /// and whether the call succeeded.
    pub fn capture(&self, args: &[&str]) -> (bool, String) {
        let _ = std::io::stdout().flush();
        match Command::new(&self.cmd)
            .args(args)
            .stdin(Stdio::null())
            .stderr(Stdio::null())
            .output()
        {
            Ok(o) => {
                let s = String::from_utf8_lossy(&o.stdout);
                (o.status.success(), s.trim_end_matches('\n').to_owned())
            }
            Err(_) => (false, String::new()),
        }
    }

    /// `engine image inspect IMAGE &>/dev/null`.
    pub fn has(&self, img: &str) -> bool {
        self.quiet(&["image", "inspect", img])
    }

    fn id(&self, img: &str) -> (bool, String) {
        self.capture(&["image", "inspect", "--format", "{{.Id}}", img])
    }
}

/// pipeline_images: the value of every `NAME_IMAGE=` line of versions.env,
/// one per line, except the lines marked `# optional`.
pub fn pipeline_images(root: &Path) -> Result<Vec<String>, String> {
    let p = root.join("versions.env");
    let text = std::fs::read_to_string(&p)
        .map_err(|e| format!("ERROR: cannot read {}: {e}", p.display()))?;
    let pins = crate::fetch::versions_pins(&text);
    Ok(text
        .lines()
        .filter(|l| !l.contains("# optional"))
        .filter_map(|l| {
            let (name, _) = l.split_once('=')?;
            let ok = name.len() > "_IMAGE".len()
                && name.ends_with("_IMAGE")
                && name
                    .bytes()
                    .all(|b| b.is_ascii_uppercase() || b.is_ascii_digit() || b == b'_');
            ok.then(|| pins.get(name).cloned().unwrap_or_default())
        })
        .collect())
}

/// pinned_tag: the local tag of an image pinned by digest
/// (NAME[:TAG]@sha256:HEX): localhost/pgp-pinned/NAME:sha256-HEX, the tag
/// dropped and a registry port's colon made a dash; `None` for an image
/// pinned by tag alone.
pub fn pinned_tag(img: &str) -> Option<String> {
    let (name, hex) = img.split_once("@sha256:")?;
    let first = name.bytes().next()?;
    if !(first.is_ascii_lowercase() || first.is_ascii_digit())
        || name.contains('@')
        || hex.len() != 64
        || !hex
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
    {
        return None;
    }
    let last = name.rsplit('/').next().unwrap_or(name);
    let name = if last.contains(':') {
        &name[..name.rfind(':').unwrap_or(name.len())]
    } else {
        name
    };
    Some(format!(
        "localhost/pgp-pinned/{}:sha256-{hex}",
        name.replace(':', "-")
    ))
}

/// tag_pinned_image: give IMAGE its pinned tag; true for an image pinned by
/// tag alone.
pub fn tag_pinned_image(e: &Engine, img: &str) -> bool {
    match pinned_tag(img) {
        None => true,
        Some(t) => e.loud(&["tag", img, &t]),
    }
}

/// pinned_tag_present: IMAGE needs no pinned tag, or its pinned tag names
/// the same image id (two failed reads count as equal, as `[ "" = "" ]`).
pub fn pinned_tag_present(e: &Engine, img: &str) -> bool {
    match pinned_tag(img) {
        None => true,
        Some(t) => e.id(&t).1 == e.id(img).1,
    }
}

/// image_is_dangling: 0 when the engine lists IMAGE's id as dangling, 1 when
/// it does not, 2 when it gives no id or no list.
pub fn image_is_dangling(e: &Engine, img: &str) -> u8 {
    let (ok, id) = e.id(img);
    if !ok || id.is_empty() {
        return 2;
    }
    let (ok, dangling) = e.capture(&[
        "images",
        "--no-trunc",
        "--quiet",
        "--filter",
        "dangling=true",
    ]);
    if !ok {
        return 2;
    }
    let found = id.lines().any(|want| dangling.lines().any(|l| l == want));
    u8::from(!found)
}

/// pull_image: `engine pull IMAGE`, tried FETCH_TRIES times (default 3),
/// FETCH_WAIT seconds apart (default 5, then 10, ...), with `sleep` as the
/// shell ran it.
pub fn pull_image(e: &Engine, img: &str, v: &Vars, err: &mut dyn Write) -> bool {
    let raw = v.get("FETCH_TRIES").unwrap_or("").to_owned();
    let t = if raw.is_empty() { "3" } else { raw.as_str() };
    let tries = if t.bytes().all(|b| b.is_ascii_digit()) {
        t.trim_start_matches('0').parse::<u64>().unwrap_or(0)
    } else {
        0
    };
    if tries == 0 {
        let _ = writeln!(
            err,
            "ERROR: FETCH_TRIES must be a whole number above 0, got '{raw}'"
        );
        return false;
    }
    let wait = v.get("FETCH_WAIT").unwrap_or("").to_owned();
    if !wait.bytes().all(|b| b.is_ascii_digit()) {
        let _ = writeln!(
            err,
            "ERROR: FETCH_WAIT must be a whole number of seconds, got '{wait}'"
        );
        return false;
    }
    for i in 1..=tries {
        if e.loud(&["pull", img]) {
            return true;
        }
        let _ = writeln!(
            err,
            "  Pull attempt {i}/{tries} failed: {img} (Docker's message is above)"
        );
        if i < tries {
            let secs = if wait.is_empty() {
                (i * 5).to_string()
            } else {
                wait.clone()
            };
            let ran = Command::new("sleep")
                .arg(&secs)
                .stdin(Stdio::null())
                .status()
                .is_ok_and(|s| s.success());
            if !ran && let Ok(n) = secs.parse::<u64>() {
                std::thread::sleep(std::time::Duration::from_secs(n));
            }
        }
    }
    false
}

/// keep_image: tag a digest-pinned image so a dangling-image cleanup keeps
/// it; when the tag fails, warn with what the engine lists.
pub fn keep_image(e: &Engine, img: &str, err: &mut dyn Write) -> bool {
    if tag_pinned_image(e, img) {
        return true;
    }
    let msg = match image_is_dangling(e, img) {
        0 => "Docker lists it as dangling, so a dangling-image cleanup would delete it",
        1 => {
            "it has no pgp-pinned tag, but Docker does not list it as dangling, so a dangling-image cleanup keeps it"
        }
        _ => "Docker did not list its dangling images, so a dangling-image cleanup may delete it",
    };
    let _ = writeln!(
        err,
        "  WARNING: Could not tag {img} as {}; {msg}.",
        pinned_tag(img).unwrap_or_default()
    );
    false
}

/// ensure_image: the image of an opt-in step (IMAGE, from the versions.env
/// variable VAR), pulled with its data so a run does not pull it midway.
pub fn ensure_image(
    e: &Engine,
    img: &str,
    step: &str,
    v: &Vars,
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> bool {
    if e.has(img) {
        let _ = writeln!(out, "[OK] {img} (step {step}) already pulled.");
    } else if pull_image(e, img, v, err) {
        let _ = writeln!(out, "[OK] {img} (step {step}) pulled.");
    } else {
        let _ = writeln!(
            out,
            "[WARN] Could not pull {img}. Step {step} needs it; run this command again."
        );
        return false;
    }
    keep_image(e, img, err)
}

/// pull_images: pull every image of [`pipeline_images`] that is not there,
/// tag the digest-pinned ones (those already there included), and print the
/// counts; false when one failed to pull or tag.
pub fn pull_images(
    e: &Engine,
    images: &[String],
    v: &Vars,
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> bool {
    let (mut pulled, mut skipped, mut failed, mut tagged) = (0, 0, 0, 0);
    for img in images {
        if e.has(img) {
            skipped += 1;
        } else {
            let _ = writeln!(out, "  Pulling: {img}...");
            if pull_image(e, img, v, err) {
                pulled += 1;
            } else {
                let _ = writeln!(
                    out,
                    "  WARNING: Failed to pull {img}; Docker's message is above."
                );
                failed += 1;
                continue;
            }
        }
        if pinned_tag_present(e, img) {
            continue;
        }
        if keep_image(e, img, err) {
            tagged += 1;
        } else {
            failed += 1;
        }
    }
    let _ = writeln!(
        out,
        "[OK] Docker images: {pulled} pulled, {skipped} already present, {failed} failed, {tagged} pgp-pinned tag(s) added."
    );
    if failed > 0 {
        let _ = writeln!(out);
        let _ = writeln!(
            out,
            "ERROR: {failed} Docker image(s) failed to pull or tag. Fix the issues above and re-run setup."
        );
        return false;
    }
    true
}

#[cfg(test)]
mod tests {
    use super::pinned_tag;

    #[test]
    fn pinned_tags_of_the_reference_shapes() {
        let h = "a".repeat(64);
        assert_eq!(pinned_tag("ubuntu:24.04"), None);
        assert_eq!(
            pinned_tag(&format!("ubuntu:24.04@sha256:{h}")).unwrap(),
            format!("localhost/pgp-pinned/ubuntu:sha256-{h}")
        );
        assert_eq!(
            pinned_tag(&format!("quay.io/biocontainers/eklipse@sha256:{h}")).unwrap(),
            format!("localhost/pgp-pinned/quay.io/biocontainers/eklipse:sha256-{h}")
        );
        assert_eq!(
            pinned_tag(&format!("localhost:5000/team/tool:1.2@sha256:{h}")).unwrap(),
            format!("localhost/pgp-pinned/localhost-5000/team/tool:sha256-{h}")
        );
        assert_eq!(
            pinned_tag(&format!("localhost:5000/tool@sha256:{h}")).unwrap(),
            format!("localhost/pgp-pinned/localhost-5000/tool:sha256-{h}")
        );
        assert_eq!(pinned_tag(&format!("Ubuntu@sha256:{h}")), None);
        assert_eq!(pinned_tag("ubuntu@sha256:abc"), None);
    }
}
