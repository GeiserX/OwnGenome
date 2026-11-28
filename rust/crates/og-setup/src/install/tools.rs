//! The reference-derived data and the tool installs og runs:
//! `--masked-reference`, `--megane-kmers`, `--locityper-data`,
//! `--bwa-index`, `--cyrius`, `--parascopy-data`, `--yleaf-data` and
//! `--pypgx-data`. Each does what its bash function in scripts/setup.sh
//! does: the same downloads (through [`Run::fetch`]), the same containers
//! with the same arguments ([`RunIn`]), the same files and stamps, the same
//! messages and the same exit status. tests/setup/ab.sh compares the trees
//! they leave with bash's; tests/tools.rs holds the checks against the
//! shell's own code.
//!
//! Three things are not containers: the masked reference is written by og
//! in process (og-intervals' mask_fasta, which `og setup mask-fasta` ran in
//! OG_IMAGE), so it no longer needs og's Linux build; the archives of
//! Parascopy and Yleaf are unpacked by the host's `tar`, and the PyPGx bundle
//! is cloned by the host's `git`, as setup.sh did. MEGAnE's k-mer set is
//! still built by og in OG_IMAGE under `--memory 6g`, after the check that
//! the engine has that memory (docker's `{{.MemTotal}}`, podman's
//! `{{.Host.MemTotal}}`).
//!
//! A folder setup.sh made with `mktemp -d` and removed in an EXIT trap (the
//! BWA index's and the PyPGx bundle's stage) is held here as a
//! [`lock::Lock`]: a unique new folder, removed when the installer returns
//! and when SIGINT, SIGTERM or SIGHUP ends og, as the trap did.

use super::{Run, nonempty, sha256};
use crate::engine::{self, Place, RunIn};
use crate::images::{self, Engine};
use crate::lock;
use sha2::{Digest, Sha256};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

/// The keys this module runs.
pub const PORTED: &[&str] = &[
    "masked-reference",
    "megane-kmers",
    "locityper-data",
    "bwa-index",
    "cyrius",
    "parascopy-data",
    "yleaf-data",
    "pypgx-data",
];

/// Parascopy's two archives of Zenodo record PARASCOPY_DATA_RECORD (the
/// `{part}` of row parascopy_data) and each one's md5; setup.sh holds the
/// same list in PARASCOPY_FILES while it has its bash installer
/// (tests/tools.rs compares them).
pub const PARASCOPY_FILES: [(&str, &str); 2] = [
    ("GRCh38", "a95bf674f43317d3a4c1b8ddbb140945"),
    ("models_GRCh38_1KGP", "244110d8fa883cf11334527ec2383498"),
];

/// The memory MEGAnE's k-mer build runs under, in GB (`--memory 6g`), and
/// below which the engine is refused.
pub const MEGANE_NEED_GB: u64 = 6;

/// The script `bash -c` runs in LOCITYPER_IMAGE to read the panel over the
/// loci (its arguments: the panel URL, the local .tbi, the output, the
/// regions), as setup.sh spelled it.
pub const LOCITYPER_PANEL_SCRIPT: &str = "\n        url=$1 idx=$2 out=$3; shift 3\n        tabix -h \"${url}##idx##${idx}\" \"$@\" | bgzip > \"${out}.part\"\n        tabix -f -p vcf \"${out}.part\"\n        mv \"${out}.part.tbi\" \"${out}.tbi\"; mv \"${out}.part\" \"$out\"";

/// The script `bash -c` runs in LOCITYPER_IMAGE to build the database (its
/// arguments: the folder, the reference, the loci, the threads, the log).
pub const LOCITYPER_BUILD_SCRIPT: &str = "\n        locityper target -d \"$1/db.part\" -r \"$2\" -j \"$1/counts.jf\" -v \"$1/panel.vcf.gz\" -g GRCh38 -L \"$3\" -@ \"$4\" > \"$5\" 2>&1\n        locityper augment -d \"$1/db.part\" -@ \"$4\" >> \"$5\" 2>&1";

/// Runs KEY when it is one of [`PORTED`]; its exit status.
pub fn run(key: &str, r: &Run, out: &mut dyn Write, err: &mut dyn Write) -> Option<u8> {
    Some(match key {
        "masked-reference" => masked_reference(r, out, err),
        "megane-kmers" => {
            if ensure(r, "MEGANE_IMAGE", "50", out, err) {
                megane_kmers(r, out, err)
            } else {
                1
            }
        }
        "locityper-data" => locityper_data(r, out, err),
        "bwa-index" => bwa_index(r, out, err),
        "cyrius" => {
            if r.vars.get("SETUP_EVERYTHING_SCOPE") != Some("row") {
                let _ = writeln!(out, "  {}", crate::components::CYRIUS_TERMS);
            }
            if ensure(r, "PYTHON_IMAGE", "21", out, err) {
                cyrius(r, out, err)
            } else {
                1
            }
        }
        "parascopy-data" => match parascopy_data(r, out, err) {
            0 => u8::from(!ensure(r, "PARASCOPY_IMAGE", "35", out, err)),
            rc => rc,
        },
        "yleaf-data" => match yleaf_data(r, out, err) {
            0 => u8::from(!ensure(r, "YLEAF_IMAGE", "37", out, err)),
            rc => rc,
        },
        "pypgx-data" => pypgx_data(r, out, err),
        _ => return None,
    })
}

// --- shell helpers, each as setup.sh's line did it ---------------------------------

fn p(s: &str) -> PathBuf {
    PathBuf::from(s)
}

fn show(p: &Path) -> String {
    p.to_string_lossy().into_owned()
}

/// An OS error as coreutils prints it: without Rust's " (os error N)".
fn os(e: &std::io::Error) -> String {
    let s = e.to_string();
    match s.find(" (os error ") {
        Some(i) => s[..i].to_owned(),
        None => s,
    }
}

/// `dirname PATH`.
pub fn dirname(path: &str) -> String {
    let t = path.trim_end_matches('/');
    if t.is_empty() {
        return if path.starts_with('/') { "/" } else { "." }.to_owned();
    }
    match t.rfind('/') {
        Some(i) => {
            let d = t[..i].trim_end_matches('/');
            if d.is_empty() {
                "/".to_owned()
            } else {
                d.to_owned()
            }
        }
        None => ".".to_owned(),
    }
}

/// `basename PATH`.
pub fn basename(path: &str) -> String {
    let t = path.trim_end_matches('/');
    if t.is_empty() {
        return if path.is_empty() {
            String::new()
        } else {
            "/".to_owned()
        };
    }
    t.rsplit('/').next().unwrap_or(t).to_owned()
}

/// `${S%.*}`: S without its last '.' and what follows it (S itself when it
/// has no '.').
pub fn strip_ext(s: &str) -> &str {
    s.rfind('.').map_or(s, |i| &s[..i])
}

/// `${S#PREFIX/}`.
fn strip_dir<'a>(s: &'a str, dir: &str) -> &'a str {
    s.strip_prefix(dir)
        .and_then(|r| r.strip_prefix('/'))
        .unwrap_or(s)
}

/// `$(cat FILE 2>/dev/null)`: its text without trailing newlines, "" when
/// it cannot be read.
fn cat(f: &Path) -> String {
    std::fs::read(f)
        .map(|b| {
            String::from_utf8_lossy(&b)
                .trim_end_matches('\n')
                .to_owned()
        })
        .unwrap_or_default()
}

/// `rm -rf PATH` (`rm -f` for a file): true when nothing is left.
fn rm_rf(path: &Path) -> bool {
    match std::fs::symlink_metadata(path) {
        Err(_) => true,
        Ok(m) if m.is_dir() => std::fs::remove_dir_all(path).is_ok(),
        Ok(_) => std::fs::remove_file(path).is_ok(),
    }
}

/// `rm -rf PATH` with rm's message on ERR when something stays.
fn rm_rf_said(path: &Path, err: &mut dyn Write) -> bool {
    let left = match std::fs::symlink_metadata(path) {
        Err(_) => return true,
        Ok(m) if m.is_dir() => std::fs::remove_dir_all(path),
        Ok(_) => std::fs::remove_file(path),
    };
    match left {
        Ok(()) => true,
        Err(e) => {
            let _ = writeln!(err, "rm: cannot remove '{}': {}", show(path), os(&e));
            false
        }
    }
}

/// `mkdir -p DIR`, with mkdir's message on ERR when it fails.
fn mkdir(dir: &Path, err: &mut dyn Write) -> bool {
    match std::fs::create_dir_all(dir) {
        Ok(()) => true,
        Err(e) => {
            let _ = writeln!(
                err,
                "mkdir: cannot create directory '{}': {}",
                show(dir),
                os(&e)
            );
            false
        }
    }
}

/// `mv -f FROM TO`: a rename, or the host's `mv` when the two are on
/// different filesystems; mv's message on ERR when it fails.
fn mv(from: &Path, to: &Path, err: &mut dyn Write) -> bool {
    match std::fs::rename(from, to) {
        Ok(()) => true,
        Err(e) if e.raw_os_error() == Some(rustix::io::Errno::XDEV.raw_os_error()) => {
            let _ = std::io::stdout().flush();
            Command::new("mv")
                .arg("-f")
                .arg("--")
                .arg(from)
                .arg(to)
                .stdin(Stdio::null())
                .status()
                .is_ok_and(|s| s.success())
        }
        Err(e) => {
            let _ = writeln!(
                err,
                "mv: cannot move '{}' to '{}': {}",
                show(from),
                show(to),
                os(&e)
            );
            false
        }
    }
}

/// `cp FROM TO`: the file's bytes, and its mode less the umask (077), as
/// cp gives a new file.
fn cp(from: &Path, to: &Path, err: &mut dyn Write) -> bool {
    use std::os::unix::fs::PermissionsExt;
    let r = std::fs::copy(from, to).and_then(|_| {
        let mode = std::fs::metadata(from)?.permissions().mode() & 0o700;
        std::fs::set_permissions(to, std::fs::Permissions::from_mode(mode))
    });
    match r {
        Ok(()) => true,
        Err(e) => {
            let _ = writeln!(
                err,
                "cp: cannot copy '{}' to '{}': {}",
                show(from),
                show(to),
                os(&e)
            );
            false
        }
    }
}

/// `printf '%s\n' TEXT > FILE`, with the shell's message on ERR.
fn write_line(f: &Path, text: &str, err: &mut dyn Write) -> bool {
    match std::fs::write(f, format!("{text}\n")) {
        Ok(()) => true,
        Err(e) => {
            let _ = writeln!(err, "{}: {}", show(f), os(&e));
            false
        }
    }
}

/// `[ A -nt B ]`: A is there and was modified after B (or B is missing).
fn newer(a: &Path, b: &Path) -> bool {
    let Ok(ma) = std::fs::metadata(a).and_then(|m| m.modified()) else {
        return false;
    };
    match std::fs::metadata(b).and_then(|m| m.modified()) {
        Ok(mb) => ma > mb,
        Err(_) => true,
    }
}

/// The files in DIR whose name starts with PREFIX (`DIR/PREFIX*`).
fn starting_with(dir: &Path, prefix: &str) -> Vec<PathBuf> {
    let Ok(rd) = std::fs::read_dir(dir) else {
        return Vec::new();
    };
    let mut v: Vec<PathBuf> = rd
        .filter_map(Result::ok)
        .filter(|e| e.file_name().to_string_lossy().starts_with(prefix))
        .map(|e| e.path())
        .collect();
    v.sort();
    v
}

/// The sha256 of BYTES (`... | _digest sha256 -`).
pub fn sha256_of(bytes: &[u8]) -> String {
    Sha256::digest(bytes)
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect()
}

/// A program of the host, its output on the terminal; true when it exits 0.
fn host(prog: &str, args: &[&str], out: &mut dyn Write) -> bool {
    let _ = out.flush();
    let _ = std::io::stdout().flush();
    Command::new(prog)
        .args(args)
        .stdin(Stdio::null())
        .status()
        .is_ok_and(|s| s.success())
}

/// `run_in [--net] --rw DIR... REST`: the engine's exit code.
fn run_in(
    r: &Run,
    net: bool,
    rw: &[&str],
    rest: &[String],
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> i32 {
    let _ = out.flush();
    let ri = RunIn {
        net,
        root: false,
        rw: rw.iter().map(|d| (*d).to_owned()).collect(),
    };
    let e = Engine::from_vars(&r.vars);
    match ri.run(&e.cmd, &Place::from_vars(&r.vars), rest) {
        Ok(c) => c,
        Err(m) => {
            let _ = writeln!(err, "{m}");
            1
        }
    }
}

/// `$(run_in ... REST)`: the exit code and stdout without its trailing
/// newlines; stderr goes to the terminal.
fn run_in_capture(
    r: &Run,
    rw: &[&str],
    rest: &[String],
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> (i32, String) {
    let _ = out.flush();
    let _ = std::io::stdout().flush();
    let ri = RunIn {
        net: false,
        root: false,
        rw: rw.iter().map(|d| (*d).to_owned()).collect(),
    };
    let e = Engine::from_vars(&r.vars);
    let argv = match ri.argv(&Place::from_vars(&r.vars), rest) {
        Ok(a) => a,
        Err(m) => {
            let _ = writeln!(err, "{m}");
            return (1, String::new());
        }
    };
    match Command::new(&e.cmd)
        .args(&argv)
        .stdin(Stdio::inherit())
        .stderr(Stdio::inherit())
        .output()
    {
        Ok(o) => (
            o.status.code().unwrap_or(1),
            String::from_utf8_lossy(&o.stdout)
                .trim_end_matches('\n')
                .to_owned(),
        ),
        Err(_) => (1, String::new()),
    }
}

/// `cpath HOST`, or its message on ERR.
fn cpath(r: &Run, host: &str, err: &mut dyn Write) -> Option<String> {
    match engine::cpath(&r.genome, r.vars.nonempty("OG_DOWNLOADS"), host) {
        Ok(c) => Some(c),
        Err(m) => {
            let _ = writeln!(err, "{m}");
            None
        }
    }
}

/// `cpath HOST` where HOST is inside the install folder by construction.
fn c(r: &Run, host: &str) -> String {
    engine::cpath(&r.genome, r.vars.nonempty("OG_DOWNLOADS"), host).unwrap_or_default()
}

/// `ensure_image VAR STEP`.
fn ensure(r: &Run, var: &str, step: &str, out: &mut dyn Write, err: &mut dyn Write) -> bool {
    let _ = out.flush();
    images::ensure_image(
        &Engine::from_vars(&r.vars),
        &r.version(var),
        step,
        &r.vars,
        out,
        err,
    )
}

/// `source_of ID --part PART FIELD` (path or url) of common.sh; its error
/// goes to ERR.
fn source_part(
    r: &Run,
    id: &str,
    part: Option<&str>,
    field: &str,
    err: &mut dyn Write,
) -> Option<String> {
    let got = (|| -> Result<String, String> {
        let ctx = r.ctx()?;
        let s = crate::fetch::Sources::read(&ctx.root)?;
        let row = s.get(id)?;
        match field {
            "path" => Ok(ctx.path(row, part)?.to_string_lossy().into_owned()),
            "url" => ctx.url(row, part),
            f => Err(format!("ERROR: unknown field {f}")),
        }
    })();
    match got {
        Ok(v) => Some(v),
        Err(m) => {
            let _ = writeln!(err, "{m}");
            None
        }
    }
}

/// `fetch URL DEST KIND VALUE > /dev/null 2>&1`.
fn fetch_quiet(r: &Run, url: &str, dest: &Path, kind: &str, want: &str) -> bool {
    r.fetch(url, dest, kind, want, &mut std::io::sink())
}

/// The lock of an installer, or `None` after `ERROR: ...` on ERR.
fn hold(path: &Path, out: &mut dyn Write, err: &mut dyn Write) -> Option<lock::Lock> {
    match lock::acquire(path, out) {
        Ok(l) => Some(l),
        Err(e) => {
            let _ = writeln!(err, "ERROR: {e}");
            None
        }
    }
}

/// `mktemp -d DIR/PREFIX.XXXXXX`, held: a new folder removed when the
/// guard drops and when a signal ends og (an EXIT trap's `rm -rf`).
/// `None` after mktemp's message on ERR.
fn stage(dir: &Path, prefix: &str, err: &mut dyn Write) -> Option<(lock::Lock, PathBuf)> {
    const CHARS: &[u8] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";
    let mut seed = u64::from(std::process::id())
        ^ std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map_or(0, |d| d.as_nanos() as u64);
    for _ in 0..100 {
        let mut name = format!("{prefix}.");
        for _ in 0..6 {
            seed = seed
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            name.push(CHARS[(seed >> 33) as usize % CHARS.len()] as char);
        }
        let path = dir.join(&name);
        if std::fs::symlink_metadata(&path).is_ok() {
            continue;
        }
        match lock::acquire(&path, &mut std::io::sink()) {
            Ok(l) => return Some((l, path)),
            Err(e) => {
                let _ = writeln!(
                    err,
                    "mktemp: failed to create directory via template '{}/{prefix}.XXXXXX': {}",
                    show(dir),
                    os(&e)
                );
                return None;
            }
        }
    }
    let _ = writeln!(
        err,
        "mktemp: failed to create directory via template '{}/{prefix}.XXXXXX': File exists",
        show(dir)
    );
    None
}

/// adopt OLD NEW [URL KIND VALUE] of common.sh: a raw download an older
/// setup left at OLD moves to NEW in the downloads folder (with its .part);
/// when NEW is there already, one of the two goes only once the check
/// (the fetch that follows, run on the file that is there) says which.
pub fn adopt(
    r: &Run,
    old: &Path,
    new: &Path,
    check: Option<(&str, &str, &str)>,
    out: &mut dyn Write,
) {
    let exists = |p: &Path| std::fs::symlink_metadata(p).is_ok();
    let is_file = |p: &Path| std::fs::metadata(p).is_ok_and(|m| m.is_file());
    let with = |p: &Path, x: &str| p.with_file_name(format!("{}{x}", basename(&show(p))));
    let (o, n) = (show(old), show(new));
    if exists(old) && exists(new) {
        // Only files are checked: fetch downloads into a NEW that is no file.
        let check = check.filter(|_| is_file(old) && is_file(new));
        if let Some((url, kind, want)) = check {
            if fetch_quiet(r, url, new, kind, want) {
                rm_rf(old);
                rm_rf(&with(old, ".part"));
                let _ = writeln!(
                    out,
                    "  Deleted {o}, left by an older setup: {n} is the same download, checked."
                );
                return;
            }
            let aside = with(new, ".adopt");
            if std::fs::rename(new, &aside).is_ok() {
                if std::fs::rename(old, new).is_ok() {
                    if fetch_quiet(r, url, new, kind, want) {
                        rm_rf(&aside);
                        let _ = writeln!(
                            out,
                            "  {n} failed its check; replaced it with {o}, left by an older setup, which passes."
                        );
                        return;
                    }
                    if std::fs::rename(new, old).is_ok() {
                        let _ = std::fs::rename(&aside, new);
                    }
                } else {
                    let _ = std::fs::rename(&aside, new);
                }
            }
        }
        let _ = writeln!(
            out,
            "  [WARN] {o}, left by an older setup, and {n} are both here; neither is deleted, as no check says which is the download."
        );
        return;
    }
    for x in ["", ".part"] {
        let ox = with(old, x);
        let nx = with(new, x);
        if !exists(&ox) {
            continue;
        }
        if x == ".part" && exists(new) {
            let passes = check.is_some_and(|(url, kind, want)| {
                is_file(new) && fetch_quiet(r, url, new, kind, want)
            });
            if passes {
                rm_rf(&ox);
                let _ = writeln!(
                    out,
                    "  Deleted {o}.part, a partial download left by an older setup: {n} is the whole download, checked."
                );
            } else {
                let _ = writeln!(
                    out,
                    "  [WARN] {o}.part, a partial download left by an older setup, stays: {n} is already there and fetch checks it. Once setup has finished, {o}.part is safe to delete."
                );
            }
        } else if exists(&nx) {
            let size = |p: &Path| std::fs::metadata(p).map_or(0, |m| m.len());
            let readable = |p: &Path| std::fs::File::open(p).is_ok();
            let aside = with(new, ".part.adopt");
            let moved = x == ".part"
                && is_file(&ox)
                && is_file(&nx)
                && readable(&ox)
                && readable(&nx)
                && size(&ox) > size(&nx)
                && std::fs::rename(&nx, &aside).is_ok()
                && (std::fs::rename(&ox, &nx).is_ok() || {
                    let _ = std::fs::rename(&aside, &nx);
                    false
                })
                && std::fs::rename(&aside, &ox).is_ok();
            if moved {
                let _ = writeln!(
                    out,
                    "  [WARN] {o}.part, left by an older setup, was the larger partial download: it resumes at {n}.part, and the smaller one is kept at {o}.part."
                );
            } else {
                let _ = writeln!(
                    out,
                    "  [WARN] {o}{x}, left by an older setup, and {n}{x} are both here; neither is deleted, and {n}{x} is the one fetch uses."
                );
            }
        } else if std::fs::create_dir_all(dirname(&n)).is_ok() && mv(&ox, &nx, &mut std::io::sink())
        {
            let _ = writeln!(out, "  Moved {o}{x}, left by an older setup, to {n}{x}.");
        }
    }
}

/// done_with FILE of common.sh: a raw download, installed. It stays in the
/// downloads folder, unless SETUP_KEEP_DOWNLOADS=false, which deletes it
/// with its .part and the checksum saved beside it.
pub fn done_with(r: &Run, f: &Path, out: &mut dyn Write) {
    if r.vars.nonempty("SETUP_KEEP_DOWNLOADS").unwrap_or("true") != "false" {
        return;
    }
    if std::fs::symlink_metadata(f).is_err() {
        return;
    }
    let part = f.with_file_name(format!("{}.part", basename(&show(f))));
    let sum = f.with_file_name(format!("{}.checksum", basename(&show(f))));
    let gone = rm_rf(f) & rm_rf(&part) & rm_rf(&sum);
    if gone {
        let _ = writeln!(
            out,
            "  Deleted the raw download {} (SETUP_KEEP_DOWNLOADS=false).",
            show(f)
        );
    }
}

/// need_og WHAT: og's Linux build in bin/.og, made when it is missing or
/// stale (`rust/og.sh check --required || rust/og.sh install`); when that
/// fails, the warning that og WHAT.
fn need_og(r: &Run, what: &str, out: &mut dyn Write) -> bool {
    let og_sh = r.root.join("rust/og.sh");
    let _ = out.flush();
    let quiet = Command::new(&og_sh)
        .args(["check", "--required"])
        .stdin(Stdio::null())
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .is_ok_and(|s| s.success());
    if quiet || host(&show(&og_sh), &["install"], out) {
        return true;
    }
    let _ = writeln!(
        out,
        "[WARN] Could not build og, which {what}: run {} --og and look at its error.",
        r.script
    );
    false
}

// --- --cyrius -----------------------------------------------------------------------

/// install_cyrius: Cyrius CYRIUS_VERSION in tools/cyrius-<v>, installed by
/// pip in PYTHON_IMAGE from scripts/cyrius-constraints.txt, hash-locked;
/// INSTALLED names the image and the lock file's sha256.
pub fn cyrius(r: &Run, out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let v = r.version("CYRIUS_VERSION");
    let img = r.version("PYTHON_IMAGE");
    let dest = format!("{}/tools/cyrius-{v}", r.genome);
    let part = format!("{dest}.part");
    let lock = show(&r.root.join("scripts/cyrius-constraints.txt"));
    let stamp = format!("python={img} lock={}", sha256(&p(&lock)));
    if cat(&p(&format!("{dest}/INSTALLED"))) == stamp {
        let _ = writeln!(out, "[OK] Cyrius {v} (step 21) already installed: {dest}");
        return 0;
    }
    let _ = writeln!(
        out,
        "Installing Cyrius {v} from {lock} (hash-locked) into {dest}"
    );
    if !rm_rf_said(&p(&part), err) || !mkdir(&p(&part), err) {
        return 1;
    }
    // net: pip installs the locked wheels from PyPI, each checked by its
    // hash (a package install, not a data source).
    let rest: Vec<String> = [
        "-v",
        &format!("{lock}:/lock.txt:ro"),
        &img,
        "pip",
        "install",
        "--no-cache-dir",
        "--disable-pip-version-check",
        "-q",
        "--require-hashes",
        "--no-deps",
        "--only-binary",
        ":all:",
        "--target",
        &c(r, &part),
        "-r",
        "/lock.txt",
    ]
    .iter()
    .map(|s| (*s).to_owned())
    .collect();
    let tools = format!("{}/tools", r.genome);
    if run_in(r, true, &[&tools], &rest, out, err) != 0 {
        rm_rf(&p(&part));
        let _ = writeln!(
            out,
            "[WARN] Could not install Cyrius. Step 21 needs it; run: {} --cyrius {}",
            r.script, r.genome
        );
        return 1;
    }
    if !write_line(&p(&format!("{part}/INSTALLED")), &stamp, err)
        || !rm_rf_said(&p(&dest), err)
        || !mv(&p(&part), &p(&dest), err)
    {
        return 1;
    }
    let _ = writeln!(out, "[OK] Cyrius {v} (step 21): {dest}");
    0
}

// --- --parascopy-data ---------------------------------------------------------------

/// install_parascopy_data: Parascopy's GRCh38 homology table and the model
/// parameters of 1000 Genomes populations (step 35), two Zenodo archives
/// each checked against its md5 ([`PARASCOPY_FILES`]), unpacked by the
/// host's tar into reference/parascopy-<v>, moved into place once whole.
pub fn parascopy_data(r: &Run, out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let v = r.version("PARASCOPY_DATA_VERSION");
    let dest = format!("{}/reference/parascopy-{v}", r.genome);
    let part = format!("{dest}.part");
    if nonempty(&p(&format!("{dest}/homology_table/GRCh38.bed.gz")))
        && nonempty(&p(&format!("{dest}/models_GRCh38_1KGP/EUR/SMN1.gz")))
    {
        let _ = writeln!(
            out,
            "[OK] Parascopy {v} homology table and models (step 35) already present."
        );
        return 0;
    }
    let urls = |err: &mut dyn Write, part_name: &str| -> Option<(String, String)> {
        let url = source_part(r, "parascopy_data", Some(part_name), "url", err)?;
        let raw = source_part(r, "parascopy_data", Some(part_name), "path", err)?;
        Some((url, raw))
    };
    for (name, md5) in PARASCOPY_FILES {
        let Some((url, raw)) = urls(err, name) else {
            return 1;
        };
        let old = p(&format!("{part}/{}", basename(&raw)));
        adopt(r, &old, &p(&raw), Some((&url, "md5", md5)), out);
    }
    rm_rf(&p(&part));
    let _ = std::fs::create_dir_all(&part);
    for (name, md5) in PARASCOPY_FILES {
        let Some((url, raw)) = urls(err, name) else {
            return 1;
        };
        if !r.fetch(&url, &p(&raw), "md5", md5, err)
            || !host("tar", &["-xzf", &raw, "-C", &part], out)
        {
            rm_rf(&p(&part));
            let _ = writeln!(
                out,
                "[WARN] Could not install Parascopy's data. Step 35 needs it; run: {} --parascopy-data {}",
                r.script, r.genome
            );
            return 1;
        }
        done_with(r, &p(&raw), out);
    }
    // Not kept: the older install, replaced.
    rm_rf(&p(&dest));
    let _ = mv(&p(&part), &p(&dest), err);
    let _ = writeln!(
        out,
        "[OK] Parascopy {v} homology table and models (step 35): {dest}"
    );
    0
}

// --- --yleaf-data -------------------------------------------------------------------

/// install_yleaf_data: the hg38 and hg_prediction_tables folders of
/// yleaf/data in Yleaf YLEAF_DATA_VERSION's tag archive (step 37), then the
/// current Y tree beside them in data/ytree: YFull's tree, Yleaf's two SNP
/// tables and the chrY strict mask, each at its pin. A failed download of
/// the current tree is a warning, not a failure; version.txt names every
/// pin, so a re-pin installs the tables again; a mask that matches its pin
/// is kept with no download.
pub fn yleaf_data(r: &Run, out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let v = r.version("YLEAF_DATA_VERSION");
    let dest = format!("{}/reference/yleaf-{v}", r.genome);
    let part = format!("{dest}.part");
    let tables = |d: &str| {
        nonempty(&p(&format!("{d}/data/hg38/new_positions.txt")))
            && nonempty(&p(&format!("{d}/data/hg_prediction_tables/tree.json")))
    };
    if tables(&dest) {
        let _ = writeln!(
            out,
            "[OK] Yleaf {v} marker tables (step 37) already present."
        );
    } else {
        let Some(url) = source_part(r, "yleaf_data", None, "url", err) else {
            return 1;
        };
        let Some(tgz) = source_part(r, "yleaf_data", None, "path", err) else {
            return 1;
        };
        let sha = r.version("YLEAF_DATA_SHA256");
        adopt(
            r,
            &p(&format!("{dest}.tar.gz")),
            &p(&tgz),
            Some((&url, "sha256", &sha)),
            out,
        );
        if !r.fetch(&url, &p(&tgz), "sha256", &sha, err) {
            let _ = writeln!(out, "[WARN] Could not download Yleaf {v}.");
            return 1;
        }
        rm_rf(&p(&part));
        let _ = std::fs::create_dir_all(&part);
        let top = format!("Yleaf-{v}/yleaf/data");
        // tar's exit status is not read: the two files below are.
        host(
            "tar",
            &[
                "-xzf",
                &tgz,
                "-C",
                &part,
                "--strip-components",
                "2",
                &format!("{top}/hg38"),
                &format!("{top}/hg_prediction_tables"),
            ],
            out,
        );
        if !tables(&part) {
            let _ = writeln!(
                out,
                "[WARN] The Yleaf {v} archive has no yleaf/data/hg38/new_positions.txt or hg_prediction_tables/tree.json."
            );
            rm_rf(&p(&part));
            return 1;
        }
        done_with(r, &p(&tgz), out);
        // Not kept: the older install, replaced.
        rm_rf(&p(&dest));
        let _ = mv(&p(&part), &p(&dest), err);
        let _ = writeln!(
            out,
            "[OK] Yleaf {v} marker tables and trees (step 37): {dest}/data"
        );
    }

    let yt = format!("{dest}/data/ytree");
    let ytp = format!("{yt}.part");
    let label = format!("YFull YTree {} (CC BY 4.0)", r.version("YTREE_VERSION"));
    let mut urls = Vec::new();
    for id in [
        "y_callable_mask",
        "ytree",
        "ytree_ftdna_positions",
        "ytree_yfull_positions",
    ] {
        let Some(u) = source_part(r, id, None, "url", err) else {
            return 1;
        };
        urls.push(u);
    }
    let [mask_url, tree_url, ftdna_url, yfull_url] = &urls[..] else {
        return 1;
    };
    let mask_sha = r.version("YCALLABLE_MASK_SHA256");
    let mask = p(&format!("{yt}/callable_mask_chrY.fasta.gz"));
    let mask_ok = nonempty(&mask) && sha256(&mask) == mask_sha;
    let want = format!(
        "tree\t{label}\npositions\tYleaf {} GRCh38 SNP tables (GPL-3.0)\npins\t{} {} {} {} {} {} {}",
        r.version("YTREE_POSITIONS_VERSION"),
        r.version("YTREE_COMMIT"),
        r.version("YTREE_SHA256"),
        r.version("YTREE_POSITIONS_COMMIT"),
        r.version("YTREE_FTDNA_SHA256"),
        r.version("YTREE_YFULL_SHA256"),
        r.version("YCALLABLE_MASK_RELEASE"),
        mask_sha
    );
    let version_txt = p(&format!("{yt}/version.txt"));
    if nonempty(&version_txt) && cat(&version_txt) == want {
        if mask_ok {
            let _ = writeln!(
                out,
                "[OK] Current Y tree ({label}, step 37) already present."
            );
            return 0;
        }
        // The tables are current but the mask is missing or fails its
        // sha256: fetch it alone. Without it step 37 skips the private SNP list.
        let _ = writeln!(
            out,
            "[WARN] {} is missing or fails its sha256; downloading it again.",
            show(&mask)
        );
        // Not kept: an installed mask that fails its sha256.
        rm_rf(&mask);
        if r.fetch(mask_url, &mask, "sha256", &mask_sha, err) {
            let _ = writeln!(
                out,
                "[OK] Current Y tree ({label}, step 37) present, chrY callable mask downloaded again."
            );
        } else {
            let _ = writeln!(
                out,
                "[WARN] Could not download the chrY callable mask; step 37 will skip the private SNP list. Run this again to retry."
            );
        }
        return 0;
    }
    rm_rf(&p(&ytp));
    let _ = std::fs::create_dir_all(&ytp);
    // A mask that matches its pin is kept: a new tree pin does not download
    // it again. The two position tables are fetched again with the tree.
    if mask_ok {
        let _ = cp(
            &mask,
            &p(&format!("{ytp}/callable_mask_chrY.fasta.gz")),
            err,
        );
    }
    let got = r.fetch(
        tree_url,
        &p(&format!("{ytp}/current_tree.json")),
        "sha256",
        &r.version("YTREE_SHA256"),
        err,
    ) && r.fetch(
        ftdna_url,
        &p(&format!("{ytp}/positions_ftdna_hg38.txt")),
        "sha256",
        &r.version("YTREE_FTDNA_SHA256"),
        err,
    ) && r.fetch(
        yfull_url,
        &p(&format!("{ytp}/positions_yfull_hg38.txt")),
        "sha256",
        &r.version("YTREE_YFULL_SHA256"),
        err,
    ) && (mask_ok
        || r.fetch(
            mask_url,
            &p(&format!("{ytp}/callable_mask_chrY.fasta.gz")),
            "sha256",
            &mask_sha,
            err,
        ));
    if !got {
        // Yleaf's own data is in place, so step 37 still runs on Yleaf's
        // tree. An older tree goes too; a mask that matches its pin stays.
        let _ = writeln!(
            out,
            "[WARN] Could not download the current Y tree ({label}); step 37 will place the sample on Yleaf's tree only. Run this again to retry."
        );
        rm_rf(&p(&ytp));
        if mask_ok {
            if let Ok(rd) = std::fs::read_dir(&yt) {
                for e in rd.filter_map(Result::ok) {
                    if e.file_name() != "callable_mask_chrY.fasta.gz" {
                        rm_rf(&e.path());
                    }
                }
            }
        } else {
            rm_rf(&p(&yt));
        }
        return 0;
    }
    let _ = write_line(&p(&format!("{ytp}/version.txt")), &want, err);
    // Not kept: the installed tree of an older pin, replaced.
    rm_rf(&p(&yt));
    let _ = mv(&p(&ytp), &p(&yt), err);
    let _ = writeln!(
        out,
        "[OK] Current Y tree ({label}) and Yleaf {}'s SNP positions (step 37): {yt}",
        r.version("YTREE_POSITIONS_VERSION")
    );
    0
}

// --- --masked-reference -------------------------------------------------------------

/// The masked reference of REF and its two companions: the FASTA, its
/// dictionary (`${dest%.*}.dict`), and the GRC's list beside them.
pub fn masked_paths(reference: &str, exclusions_url: &str) -> (String, String, String) {
    let dest = format!("{}_maskedGRC.fasta", strip_ext(reference));
    let dir = dirname(&dest);
    let dict = format!("{}.dict", strip_ext(&dest));
    let bed = format!(
        "{dir}/{}",
        exclusions_url.rsplit('/').next().unwrap_or(exclusions_url)
    );
    (dest, dict, bed)
}

/// What `og setup mask-fasta SRC BED DEST` did, in process: the summary
/// line, or the message it printed (paths named as the command named them).
fn mask_in_process(src: &str, bed: &str, part: &str) -> Result<String, String> {
    use og_intervals::{BedError, MaskError, Rows, mask_fasta};
    let (src_n, bed_n, part_n) = (basename(src), basename(bed), basename(part));
    let fail = |what: &str, e: &dyn std::fmt::Display| format!("og setup mask-fasta: {what}: {e}");
    let rows = og_io::open(bed)
        .map_err(BedError::Io)
        .and_then(Rows::read)
        .map_err(|e| fail(&bed_n, &e))?;
    let input = std::fs::File::open(src)
        .map(|f| std::io::BufReader::with_capacity(1 << 20, f))
        .map_err(|e| fail(&src_n, &e))?;
    let mut o = og_io::AtomicFile::create(part).map_err(|e| fail(&part_n, &e))?;
    let summary = match mask_fasta(input, &rows, &mut o) {
        Ok(s) => s,
        Err(MaskError::Io(e)) => return Err(fail(&format!("{src_n} -> {part_n}"), &e)),
        Err(e) => return Err(fail(&src_n, &e)),
    };
    o.commit().map_err(|e| fail(&part_n, &e))?;
    Ok(summary.to_string())
}

/// install_masked_reference: REF_FASTA with the false duplications of
/// GRCh38 (GRC_EXCLUSIONS_URL, sha256 checked) set to N, beside it as
/// `<REF without extension>_maskedGRC.fasta`, which must have
/// MASKED_REF_SHA256 or GRC_MASKED_REF_SHA256 or is deleted; then its .fai
/// and .dict by samtools in SAMTOOLS_IMAGE.
pub fn masked_reference(r: &Run, out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let src = r.vars.get("REF_FASTA").unwrap_or("").to_owned();
    if src.ends_with(".gz") {
        let _ = writeln!(
            out,
            "[WARN] --masked-reference needs an unpacked FASTA; REF_FASTA is {src}."
        );
        return 1;
    }
    if src.contains("_maskedGRC.") {
        let _ = writeln!(
            out,
            "[WARN] REF_FASTA ({src}) is already a masked reference. Unset it, so the mask is built from the default reference."
        );
        return 1;
    }
    let url = r.version("GRC_EXCLUSIONS_URL");
    let (dest, dict, bed) = masked_paths(&src, &url);
    let dir = dirname(&dest);
    let fai = format!("{dest}.fai");
    let part = format!("{dest}.part");
    let want = r
        .vars
        .nonempty("MASKED_REF_SHA256")
        .map(str::to_owned)
        .unwrap_or_else(|| r.version("GRC_MASKED_REF_SHA256"));
    if nonempty(&p(&dest)) && nonempty(&p(&fai)) && nonempty(&p(&dict)) && sha256(&p(&dest)) == want
    {
        let _ = writeln!(
            out,
            "[OK] Masked reference already present (sha256 matches): {dest}"
        );
    } else {
        for f in [src.clone(), format!("{src}.fai")] {
            if !nonempty(&p(&f)) {
                let _ = writeln!(
                    out,
                    "[WARN] {f} is missing: install the reference first ({} {}).",
                    r.script, r.genome
                );
                return 1;
            }
        }
        let bed_sha = r.version("GRC_EXCLUSIONS_SHA256");
        let have_bed = nonempty(&p(&bed)) && sha256(&p(&bed)) == bed_sha;
        if !have_bed && !r.fetch(&url, &p(&bed), "sha256", &bed_sha, err) {
            let _ = writeln!(
                out,
                "[WARN] Could not download the GRC's list of false duplications ({url})."
            );
            return 1;
        }
        // Only the partial output goes now: a masked reference already
        // there stays until its replacement has passed the sha256 check.
        rm_rf(&p(&part));
        let _ = writeln!(
            out,
            "Masking the false duplications of {} (under a minute for GRCh38)...",
            basename(&src)
        );
        let log = match mask_in_process(&src, &bed, &part) {
            Ok(l) => l,
            Err(m) => {
                let _ = writeln!(err, "{m}");
                rm_rf(&p(&part));
                let _ = writeln!(out, "[WARN] Could not write the masked reference.");
                return 1;
            }
        };
        let _ = writeln!(out, "  {log}");
        let got = sha256(&p(&part));
        if got != want {
            rm_rf(&p(&part));
            let _ = writeln!(
                out,
                "[WARN] The masked reference has sha256 {got}, not {want}: {} is not the file it is pinned for (NCBI's no-ALT analysis set). Removed it.",
                basename(&src)
            );
            return 1;
        }
        if !mv(&p(&part), &p(&dest), err) {
            return 1;
        }
        rm_rf(&p(&fai));
        rm_rf(&p(&dict));
        // The index and the dictionary go next to the FASTA, so its folder
        // is writable here.
        let img = r.version("SAMTOOLS_IMAGE");
        let s = |a: &[&str]| -> Vec<String> { a.iter().map(|x| (*x).to_owned()).collect() };
        let indexed = run_in(
            r,
            false,
            &[&dir],
            &s(&[&img, "samtools", "faidx", &c(r, &dest)]),
            out,
            err,
        ) == 0
            && run_in(
                r,
                false,
                &[&dir],
                &s(&[&img, "samtools", "dict", "-o", &c(r, &dict), &c(r, &dest)]),
                out,
                err,
            ) == 0;
        if !indexed {
            rm_rf(&p(&fai));
            rm_rf(&p(&dict));
            let _ = writeln!(
                out,
                "[WARN] Could not index the masked reference. Kept the FASTA (its sha256 matches) and removed its partial index; run setup again."
            );
            return 1;
        }
        let _ = writeln!(out, "[OK] Masked reference (sha256 matches): {dest}");
    }
    let _ = writeln!(
        out,
        "  Use it for a sample by realigning it (docs/realignment.md#the-masked-reference):"
    );
    let _ = writeln!(
        out,
        "    export REF_FASTA={}   # bash steps; the Nextflow run takes --reference {}",
        super::shell_quote(strip_dir(&dest, &r.genome)),
        super::shell_quote(&dest)
    );
    0
}

// --- --megane-kmers -----------------------------------------------------------------

/// megane_kmers_path: where step 50 reads MEGAnE's k-mer set of REF,
/// MEGANE_KERS overriding it (relative to the install folder G when it is
/// not absolute): reference/megane/<REF's name without extension>.mk.
pub fn megane_kmers_path(reference: &str, megane_kmers: Option<&str>, g: &str) -> String {
    let n = basename(reference);
    let n = n.strip_suffix(".gz").unwrap_or(&n);
    let n = strip_ext(n);
    let path = match megane_kmers.filter(|m| !m.is_empty()) {
        Some(m) => m.to_owned(),
        None => format!("{g}/reference/megane/{n}.mk"),
    };
    if path.starts_with('/') {
        path
    } else {
        format!("{g}/{path}")
    }
}

/// The engine's memory in whole GB (bytes / 1e9, cut down), as setup.sh
/// reads it: docker's `info --format {{.MemTotal}}`, and podman's
/// `{{.Host.MemTotal}}` when that is not a number; `None` when neither
/// answers with one. ASK runs `ENGINE info --format FIELD` and gives its
/// stdout when it exits 0.
pub fn engine_mem_gb(ask: &dyn Fn(&str) -> Option<String>) -> Option<u64> {
    let digits = |s: &str| !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit());
    let mut have = ask("{{.MemTotal}}").unwrap_or_default();
    if !digits(&have) {
        have = ask("{{.Host.MemTotal}}").unwrap_or_default();
    }
    // awk '$1 ~ /^[0-9]+$/ { printf "%d", $1 / 1e9 }', line by line.
    let mut s = String::new();
    for line in have.split('\n') {
        if let Some(f) = line.split_whitespace().next()
            && digits(f)
        {
            let n: f64 = f.parse().unwrap_or(0.0);
            s.push_str(&format!("{}", (n / 1e9).trunc() as i64));
        }
    }
    s.parse().ok()
}

/// install_megane_kmers: MEGAnE's k-mer set of the reference (step 50),
/// built by og (`og setup megane-kmers`) in OG_IMAGE with bin/ mounted,
/// under `--memory 6g`; one newer than the reference is kept. An engine
/// with less memory (and more than 0) is refused before og runs.
pub fn megane_kmers(r: &Run, out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let reference = r.vars.get("REF_FASTA").unwrap_or("").to_owned();
    let mk = megane_kmers_path(&reference, r.vars.nonempty("MEGANE_KMERS"), &r.genome);
    let dir = dirname(&mk);
    let stem = mk.strip_suffix(".mk").unwrap_or(&mk).to_owned();
    let mi = format!("{stem}.mi");
    if !nonempty(&p(&reference)) || !nonempty(&p(&format!("{reference}.fai"))) {
        let _ = writeln!(
            out,
            "[WARN] Step 50's k-mer set is built from the reference, and {reference} (with its .fai) is not there yet: run ./scripts/setup.sh {} first.",
            r.genome
        );
        return 1;
    }
    if nonempty(&p(&mk)) && nonempty(&p(&mi)) && newer(&p(&mk), &p(&reference)) {
        let _ = writeln!(out, "[OK] MEGAnE k-mer set (step 50) already present: {mk}");
        return 0;
    }
    let e = Engine::from_vars(&r.vars);
    let have = engine_mem_gb(&|field| {
        let (ok, s) = e.capture(&["info", "--format", field]);
        ok.then_some(s)
    });
    if let Some(have) = have
        && have > 0
        && have < MEGANE_NEED_GB
    {
        let _ = writeln!(
            out,
            "[WARN] MEGAnE's k-mer set needs about {MEGANE_NEED_GB} GB of RAM and {} has {have} GB.",
            e.cmd
        );
        let _ = writeln!(
            out,
            "       Build it on a machine with more memory (./scripts/setup.sh --megane-kmers <genome_dir> there) and copy {}/ here.",
            strip_dir(&dir, &r.genome)
        );
        return 1;
    }
    if !need_og(r, "builds MEGAnE's k-mer set", out) {
        return 1;
    }
    if !mkdir(&p(&dir), err) {
        return 1;
    }
    let _ = writeln!(
        out,
        "Building MEGAnE's k-mer set of {} (step 50; a few minutes, about 5 GB of RAM)...",
        basename(&reference)
    );
    // og writes each file under a hidden temporary name and renames it at
    // the end; a build killed on the way leaves that file behind.
    rm_rf(&p(&mk));
    rm_rf(&p(&mi));
    for prefix in [
        format!(".{}.tmp-", basename(&mk)),
        format!(".{}.mi.tmp-", basename(&stem)),
    ] {
        for f in starting_with(&p(&dir), &prefix) {
            rm_rf(&f);
        }
    }
    let bin = format!("{}:/pgp-bin:ro", show(&r.root.join("bin")));
    let mem = format!("{MEGANE_NEED_GB}g");
    let rest: Vec<String> = [
        "-v",
        &bin,
        "--memory",
        &mem,
        &r.version("OG_IMAGE"),
        "/pgp-bin/og",
        "setup",
        "megane-kmers",
        "--fa",
        &c(r, &reference),
        "--out",
        &c(r, &dir),
        "--prefix",
        &basename(&stem),
        "--mem",
        "5",
    ]
    .iter()
    .map(|s| (*s).to_owned())
    .collect();
    let (rc, log) = run_in_capture(r, &[&dir], &rest, out, err);
    if rc != 0 || !nonempty(&p(&mk)) || !nonempty(&p(&mi)) {
        rm_rf(&p(&mk));
        rm_rf(&p(&mi));
        let _ = writeln!(
            out,
            "[WARN] Could not build MEGAnE's k-mer set (step 50); og's error is above."
        );
        return 1;
    }
    let _ = writeln!(out, "  {log}");
    let _ = writeln!(
        out,
        "[OK] MEGAnE k-mer set (step 50, {} repeated 32-mers): {mk}",
        cat(&p(&mi))
    );
    0
}

// --- --locityper-data ---------------------------------------------------------------

/// awk's number of a string: its leading numeric part ("12ab" is 12), 0
/// when it has none.
pub fn awk_num(s: &str) -> f64 {
    let t = s.trim_start_matches([' ', '\t', '\n']);
    let b = t.as_bytes();
    let mut i = 0;
    if i < b.len() && (b[i] == b'+' || b[i] == b'-') {
        i += 1;
    }
    let digits_from = |mut j: usize| {
        while j < b.len() && b[j].is_ascii_digit() {
            j += 1;
        }
        j
    };
    let int_end = digits_from(i);
    let mut end = int_end;
    if end < b.len() && b[end] == b'.' {
        end = digits_from(end + 1);
    }
    if end == i || (end == i + 1 && b.get(i) == Some(&b'.')) {
        return 0.0;
    }
    if end < b.len() && (b[end] == b'e' || b[end] == b'E') {
        let mut k = end + 1;
        if k < b.len() && (b[k] == b'+' || b[k] == b'-') {
            k += 1;
        }
        let k2 = digits_from(k);
        if k2 > k {
            end = k2;
        }
    }
    t[..end].parse().unwrap_or(0.0)
}

/// The records of TEXT as awk reads them: split at '\n', a last line with
/// no newline still a record.
fn records(text: &str) -> Vec<&str> {
    let mut v: Vec<&str> = text.split('\n').collect();
    if text.ends_with('\n') || text.is_empty() {
        v.pop();
    }
    v
}

/// The first row of LOCI (named NAME) that is not chrom, start, end, name
/// (tab separated) with a whole-number start below its end inside a contig
/// of FAI, as `NAME:LINE`; `None` when every row is one. Comment lines
/// (`#`) and empty lines are skipped. setup.sh's awk, in Rust.
pub fn bad_locus(fai: &str, loci: &str, name: &str) -> Option<String> {
    let mut len: std::collections::HashMap<&str, &str> = std::collections::HashMap::new();
    for rec in records(fai) {
        let mut f = rec.split('\t');
        let chrom = f.next().unwrap_or("");
        len.insert(chrom, f.next().unwrap_or(""));
    }
    let whole = |s: &str| !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit());
    for (i, rec) in records(loci).iter().enumerate() {
        if rec.starts_with('#') || rec.is_empty() {
            continue;
        }
        let f: Vec<&str> = rec.split('\t').collect();
        let bad = f.len() < 4
            || !whole(f[1])
            || !whole(f[2])
            || f[3].is_empty()
            || awk_num(f[1]) >= awk_num(f[2])
            || !len.contains_key(f[0])
            || awk_num(f[2]) > awk_num(len[f[0]]);
        if bad {
            return Some(format!("{name}:{}", i + 1));
        }
    }
    None
}

/// The rows of LOCI as setup.sh keeps them in loci.bed: no comments, no
/// empty lines, the first four columns.
pub fn kept_loci(loci: &str) -> Vec<String> {
    records(loci)
        .into_iter()
        .filter(|r| !r.starts_with('#') && !r.is_empty())
        .map(|r| {
            let mut f: Vec<&str> = r.split('\t').collect();
            f.resize(4, "");
            f[..4].join("\t")
        })
        .collect()
}

/// jellyfish's hash size for a reference: the sum of its .fai's lengths
/// (the second column), never below 10M (`printf "%.0f"` of awk's sum).
pub fn jellyfish_size(fai: &str) -> String {
    let n: f64 = records(fai)
        .iter()
        .map(|r| r.split_whitespace().nth(1).map_or(0.0, awk_num))
        .sum();
    let n = if n < 10_000_000.0 { 10_000_000.0 } else { n };
    format!("{n:.0}")
}

/// The loci of LOCI (kept rows) 250 kb either side, merged where they meet
/// (sorted by contig name, then start), as `chrom:start-end` words, each
/// followed by a space: the regions setup.sh's awk | sort | awk | tr gave.
pub fn panel_regions(loci: &[String]) -> String {
    let mut iv: Vec<(String, f64, f64)> = loci
        .iter()
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            let s = (awk_num(f.get(1).unwrap_or(&"")) - 250_000.0).max(0.0);
            let e = awk_num(f.get(2).unwrap_or(&"")) + 250_000.0;
            ((*f.first().unwrap_or(&"")).to_owned(), s, e)
        })
        .collect();
    iv.sort_by(|a, b| a.0.as_bytes().cmp(b.0.as_bytes()).then(a.1.total_cmp(&b.1)));
    let mut s = String::new();
    let mut cur: Option<(String, f64, f64)> = None;
    for (chrom, a, b) in iv {
        match &mut cur {
            Some((c, _, e)) if *c == chrom && a <= *e => {
                if b > *e {
                    *e = b;
                }
            }
            _ => {
                if let Some((c, cs, ce)) = cur.take() {
                    s.push_str(&format!("{c}:{}-{} ", cs + 1.0, ce));
                }
                cur = Some((chrom, a, b));
            }
        }
    }
    if let Some((c, cs, ce)) = cur {
        s.push_str(&format!("{c}:{}-{} ", cs + 1.0, ce));
    }
    s
}

/// install_locityper_data: Locityper's database for step 53 in
/// reference/locityper-<panel> (LOCITYPER_DATA): jellyfish's k-mer counts
/// of the reference (counts.jf, counts.stamp), the HPRC panel over the loci
/// read by tabix over https through its pinned .tbi (kept in the downloads
/// folder by the regions), `locityper target` and `augment` into db/, the
/// loci in loci.bed; DONE names the inputs, so a changed one builds again.
pub fn locityper_data(r: &Run, out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let pv = r.version("LOCITYPER_PANEL_VERSION");
    let mut dir = r
        .vars
        .nonempty("LOCITYPER_DATA")
        .map(str::to_owned)
        .unwrap_or_else(|| format!("{}/reference/locityper-{pv}", r.genome));
    if !dir.starts_with('/') {
        dir = format!("{}/{dir}", r.genome);
    }
    let loci = r
        .vars
        .nonempty("LOCITYPER_LOCI")
        .map(str::to_owned)
        .unwrap_or_else(|| show(&r.root.join("assets/locityper_loci.bed")));
    if cpath(r, &dir, err).is_none() {
        return 1;
    }
    let reference = r.vars.get("REF_FASTA").unwrap_or("").to_owned();
    let fai = format!("{reference}.fai");
    if !nonempty(&p(&reference)) || !nonempty(&p(&fai)) {
        let _ = writeln!(
            out,
            "[WARN] Step 53's database is built from the reference, and {reference} (with its .fai) is not there yet: run ./scripts/setup.sh {} first.",
            r.genome
        );
        return 1;
    }
    if !nonempty(&p(&loci)) {
        let _ = writeln!(out, "[WARN] {loci} is not a file");
        return 1;
    }
    let _ = writeln!(
        out,
        "Building Locityper's database (step 53, panel {pv}) in {dir}"
    );
    if !ensure(r, "LOCITYPER_IMAGE", "53", out, err) {
        return 1;
    }
    if !mkdir(&p(&dir), err) {
        return 1;
    }
    let read =
        |f: &str| String::from_utf8_lossy(&std::fs::read(f).unwrap_or_default()).into_owned();
    let (fai_text, loci_text) = (read(&fai), read(&loci));
    if let Some(bad) = bad_locus(&fai_text, &loci_text, &loci) {
        let _ = writeln!(
            out,
            "[WARN] {bad}: want chrom, start, end, name (tab separated), start below end, inside a contig of {reference}"
        );
        return 1;
    }
    let kept = kept_loci(&loci_text);
    let new = format!("{dir}/loci.bed.new");
    let body: String = kept.iter().map(|l| format!("{l}\n")).collect();
    if let Err(e) = std::fs::write(&new, &body) {
        let _ = writeln!(err, "{new}: {}", os(&e));
        return 1;
    }
    let n = kept.len();
    if n == 0 {
        rm_rf(&p(&new));
        let _ = writeln!(out, "[WARN] {loci} lists no locus");
        return 1;
    }
    let img = r.version("LOCITYPER_IMAGE");
    let tbi_sha = r.version("LOCITYPER_PANEL_TBI_SHA256");
    // The reference's .fai and its sequence (a new FASTA with the same
    // contigs is a new reference too).
    let jstamp = format!(
        "fai={} fasta={} image={img}",
        sha256(&p(&fai)),
        sha256(&p(&reference))
    );
    let stamp = format!("{jstamp} panel={tbi_sha} loci={}", sha256(&p(&new)));
    let at = |name: &str| p(&format!("{dir}/{name}"));
    if cat(&at("DONE")) == stamp
        && at("db").is_dir()
        && nonempty(&at("counts.jf"))
        && nonempty(&at("loci.bed"))
    {
        rm_rf(&p(&new));
        let _ = writeln!(
            out,
            "[OK] Locityper's database (step 53, {n} loci) already built: {dir}"
        );
        return 0;
    }
    if !rm_rf_said(&at("DONE"), err) {
        return 1;
    }
    let threads = r.vars.nonempty("THREADS").unwrap_or("8").to_owned();
    let s = |a: &[&str]| -> Vec<String> { a.iter().map(|x| (*x).to_owned()).collect() };
    if !nonempty(&at("counts.jf")) || cat(&at("counts.stamp")) != jstamp {
        let _ = writeln!(
            out,
            "  k-mer counts of {} (jellyfish, several GB of RAM for GRCh38)...",
            basename(&reference)
        );
        rm_rf(&at("counts.jf"));
        rm_rf(&at("counts.stamp"));
        // The hash holds one entry per distinct k-mer, at most one per base,
        // so it is sized to the reference's length, never below 10M.
        let hsize = jellyfish_size(&fai_text);
        let counted = run_in(
            r,
            false,
            &[&dir],
            &s(&[
                &img,
                "jellyfish",
                "count",
                "--canonical",
                "--lower-count",
                "2",
                "--out-counter-len",
                "2",
                "--mer-len",
                "25",
                "--threads",
                &threads,
                "--size",
                &hsize,
                "--output",
                &c(r, &format!("{dir}/counts.jf.part")),
                &c(r, &reference),
            ]),
            out,
            err,
        ) == 0;
        if !counted {
            rm_rf(&at("counts.jf.part"));
            let _ = writeln!(
                out,
                "[WARN] jellyfish could not count the reference's k-mers"
            );
            return 1;
        }
        if !mv(&at("counts.jf.part"), &at("counts.jf"), err)
            || !write_line(&at("counts.stamp"), &jstamp, err)
        {
            return 1;
        }
    }

    let url = r.version("LOCITYPER_PANEL_URL");
    if !r.fetch(
        &format!("{url}.tbi"),
        &at("panel.remote.tbi"),
        "sha256",
        &tbi_sha,
        err,
    ) {
        let _ = writeln!(
            out,
            "[WARN] Could not download the panel's index; run: {} --locityper-data {}",
            r.script, r.genome
        );
        return 1;
    }
    let regions = panel_regions(&kept);
    let trimmed = regions.strip_suffix(' ').unwrap_or(&regions);
    // What the regions read: kept in the downloads folder by the panel URL,
    // its pin and the regions (with its sha256 beside it), so a rerun over
    // the same loci reads nothing again.
    let key = sha256_of(format!("{url} {tbi_sha} {trimmed}\n").as_bytes());
    let dl = r.vars.nonempty("OG_DOWNLOADS").unwrap_or("").to_owned();
    let saved = format!("{dl}/locityper-{pv}/panel-{}.vcf.gz", &key[..16]);
    let (saved_tbi, saved_sha) = (format!("{saved}.tbi"), format!("{saved}.sha256"));
    if nonempty(&p(&saved)) && nonempty(&p(&saved_tbi)) && sha256(&p(&saved)) == cat(&p(&saved_sha))
    {
        let _ = writeln!(
            out,
            "  The panel's haplotypes over these loci are already in {saved} (sha256 matches); not read again."
        );
        if !(cp(&p(&saved), &at("panel.vcf.gz.part"), err)
            && cp(&p(&saved_tbi), &at("panel.vcf.gz.tbi"), err)
            && mv(&at("panel.vcf.gz.part"), &at("panel.vcf.gz"), err))
        {
            return 1;
        }
    } else {
        let _ = writeln!(out, "  The panel's haplotypes over {trimmed}...");
        // net: the regions of the loci only, saved below.
        let mut rest = s(&[
            &img,
            "bash",
            "-euo",
            "pipefail",
            "-c",
            LOCITYPER_PANEL_SCRIPT,
            "_",
            &url,
            &c(r, &format!("{dir}/panel.remote.tbi")),
            &c(r, &format!("{dir}/panel.vcf.gz")),
        ]);
        rest.extend(regions.split_whitespace().map(str::to_owned));
        if run_in(r, true, &[&dir], &rest, out, err) != 0 {
            rm_rf(&at("panel.vcf.gz.part"));
            rm_rf(&at("panel.vcf.gz.part.tbi"));
            let _ = writeln!(
                out,
                "[WARN] Could not read the panel over https; run: {} --locityper-data {}",
                r.script, r.genome
            );
            return 1;
        }
        if r.vars.nonempty("SETUP_KEEP_DOWNLOADS").unwrap_or("true") != "false" {
            let sink = &mut std::io::sink();
            let kept = std::fs::create_dir_all(dirname(&saved)).is_ok()
                && cp(&at("panel.vcf.gz"), &p(&format!("{saved}.part")), sink)
                && cp(&at("panel.vcf.gz.tbi"), &p(&saved_tbi), sink)
                && write_line(&p(&saved_sha), &sha256(&p(&format!("{saved}.part"))), sink)
                && mv(&p(&format!("{saved}.part")), &p(&saved), sink);
            if !kept {
                let _ = writeln!(
                    out,
                    "[WARN] Could not keep the panel's haplotypes in {saved}; the next run reads them again."
                );
            }
        }
    }

    let _ = writeln!(
        out,
        "  locityper target and augment on {n} loci (the pairwise alignments take the longest)..."
    );
    rm_rf(&at("db"));
    rm_rf(&at("db.part"));
    let built = run_in(
        r,
        false,
        &[&dir],
        &s(&[
            &img,
            "bash",
            "-euo",
            "pipefail",
            "-c",
            LOCITYPER_BUILD_SCRIPT,
            "_",
            &c(r, &dir),
            &c(r, &reference),
            &c(r, &new),
            &threads,
            &c(r, &format!("{dir}/build.log")),
        ]),
        out,
        err,
    ) == 0;
    if !built {
        let _ = writeln!(
            out,
            "[WARN] Locityper could not build the database; the last lines of {dir}/build.log:"
        );
        host("tail", &["-n", "20", &format!("{dir}/build.log")], out);
        return 1;
    }
    let missing: Vec<String> = kept
        .iter()
        .map(|l| l.split('\t').nth(3).unwrap_or("").to_owned())
        .filter(|name| !p(&format!("{dir}/db.part/loci/{name}/success")).exists())
        .collect();
    if !missing.is_empty() {
        let _ = writeln!(
            out,
            "[WARN] Locityper built no database for: {} (see {dir}/build.log)",
            missing.join(" ")
        );
        return 1;
    }
    if !mv(&at("db.part"), &at("db"), err)
        || !mv(&p(&new), &at("loci.bed"), err)
        || !write_line(&at("DONE"), &stamp, err)
    {
        return 1;
    }
    let _ = out.flush();
    let du = Command::new("du")
        .args(["-sh", &dir])
        .stdin(Stdio::null())
        .output()
        .map(|o| {
            String::from_utf8_lossy(&o.stdout)
                .split('\t')
                .next()
                .unwrap_or("")
                .trim_end_matches('\n')
                .to_owned()
        })
        .unwrap_or_default();
    let _ = writeln!(
        out,
        "[OK] Locityper's database (step 53, {n} loci, panel {pv}): {dir} ({du})"
    );
    0
}

// --- --bwa-index --------------------------------------------------------------------

/// has_bwa_index PREFIX: the five files of a classic BWA index, not empty.
fn has_bwa_index(prefix: &str) -> bool {
    ["amb", "ann", "bwt", "pac", "sa"]
        .iter()
        .all(|ext| nonempty(&p(&format!("{prefix}.{ext}"))))
}

/// install_bwa_index: the classic BWA index of the reference for GRIDSS
/// (step 04b), under its lock, built by `bwa index` in BWA_IMAGE into a
/// hidden stage beside the reference, then moved beside it with `.sa` last,
/// so a half-moved index never looks complete.
pub fn bwa_index(r: &Run, out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let reference = r.vars.get("REF_FASTA").unwrap_or("").to_owned();
    let dir = dirname(&reference);
    if !nonempty(&p(&reference)) {
        let _ = writeln!(err, "ERROR: install the reference first: {reference}");
        return 1;
    }
    let Some(_held) = hold(&p(&format!("{reference}.bwa-index.lock")), out, err) else {
        return 1;
    };
    if has_bwa_index(&reference) {
        let _ = writeln!(out, "[OK] Classic BWA index already present: {reference}");
        return 0;
    }
    if !ensure(r, "BWA_IMAGE", "04b", out, err) {
        return 1;
    }
    let Some((_stage, st)) = stage(&p(&dir), ".bwa-index", err) else {
        return 1;
    };
    let prefix = format!("{}/index", show(&st));
    let _ = writeln!(
        out,
        "Building the classic BWA index for GRIDSS (about 5 GB on disk and 5.4 GB RAM)..."
    );
    let rest: Vec<String> = [
        r.version("BWA_IMAGE"),
        "bwa".into(),
        "index".into(),
        "-p".into(),
        c(r, &prefix),
        c(r, &reference),
    ]
    .to_vec();
    if run_in(r, false, &[&show(&st)], &rest, out, err) != 0 {
        return 1;
    }
    if !has_bwa_index(&prefix) {
        let _ = writeln!(err, "ERROR: bwa index did not write all five files.");
        return 1;
    }
    // Publish .sa last so a partial move cannot look complete to GRIDSS.
    if !rm_rf_said(&p(&format!("{reference}.sa")), err) {
        return 1;
    }
    for ext in ["amb", "ann", "bwt", "pac", "sa"] {
        if !mv(
            &p(&format!("{prefix}.{ext}")),
            &p(&format!("{reference}.{ext}")),
            err,
        ) {
            return 1;
        }
    }
    let _ = writeln!(out, "[OK] Classic BWA index: {reference}");
    0
}

// --- --pypgx-data -------------------------------------------------------------------

/// `git ARGS` of the host (PATH as the install's environment has it): its
/// stdout without trailing newlines when it exits 0.
fn git(r: &Run, args: &[&str], quiet: bool, out: &mut dyn Write) -> Option<String> {
    let _ = out.flush();
    let _ = std::io::stdout().flush();
    let mut cmd = Command::new("git");
    cmd.args(args).stdin(Stdio::null());
    if let Some(path) = r.vars.get("PATH") {
        cmd.env("PATH", path);
    }
    if quiet {
        cmd.stderr(Stdio::null());
    }
    let o = cmd.output().ok()?;
    o.status.success().then(|| {
        String::from_utf8_lossy(&o.stdout)
            .trim_end_matches('\n')
            .to_owned()
    })
}

/// install_pypgx_data: the upstream PyPGx bundle (step 32) in
/// reference/pypgx-bundle, under its lock: the tag PYPGX_BUNDLE_VERSION
/// cloned by the host's git into a hidden stage, refused unless its commit
/// is PYPGX_BUNDLE_COMMIT. A bundle there at the pinned commit is kept; one
/// at another commit is refused (moved aside by hand).
pub fn pypgx_data(r: &Run, out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let dir = format!("{}/reference", r.genome);
    let dest = format!("{dir}/pypgx-bundle");
    let v = r.version("PYPGX_BUNDLE_VERSION");
    let commit = r.version("PYPGX_BUNDLE_COMMIT");
    if !mkdir(&p(&dir), err) {
        return 1;
    }
    let Some(_held) = hold(&p(&format!("{dest}.setup.lock")), out, err) else {
        return 1;
    };
    if p(&dest).is_dir() {
        let safe = format!("safe.directory={dest}");
        let head = git(
            r,
            &["-c", &safe, "-C", &dest, "rev-parse", "HEAD"],
            true,
            out,
        );
        if head.as_deref() == Some(commit.as_str()) {
            return 0;
        }
        let _ = writeln!(
            err,
            "ERROR: move the existing PyPGx bundle aside before installing {v}: {dest}"
        );
        return 1;
    }
    let Some((_stage, st)) = stage(&p(&dir), ".pypgx", err) else {
        return 1;
    };
    let bundle = format!("{}/bundle", show(&st));
    let url = r.version("PYPGX_BUNDLE_URL");
    let _ = out.flush();
    let mut clone = Command::new("git");
    clone
        .args(["clone", "--depth", "1", "--branch", &v, &url, &bundle])
        .stdin(Stdio::null());
    if let Some(path) = r.vars.get("PATH") {
        clone.env("PATH", path);
    }
    let _ = std::io::stdout().flush();
    if !clone.status().is_ok_and(|s| s.success()) {
        return 1;
    }
    let head = git(r, &["-C", &bundle, "rev-parse", "HEAD"], false, out).unwrap_or_default();
    if head != commit {
        let _ = writeln!(
            err,
            "ERROR: PyPGx bundle tag {v} differs from its pinned commit."
        );
        return 1;
    }
    if !mv(&p(&bundle), &p(&dest), err) {
        return 1;
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn paths_are_cut_as_the_shell_cuts_them() {
        assert_eq!(dirname("/a/b/c.fa"), "/a/b");
        assert_eq!(dirname("c.fa"), ".");
        assert_eq!(dirname("/c.fa"), "/");
        assert_eq!(basename("/a/b/c.fa"), "c.fa");
        assert_eq!(strip_ext("/a.b/ref"), "/a");
        assert_eq!(strip_ext("/g/ref.fasta"), "/g/ref");
        assert_eq!(strip_dir("/g/reference/x", "/g"), "reference/x");
        assert_eq!(strip_dir("/other/x", "/g"), "/other/x");
    }

    #[test]
    fn awk_reads_the_leading_number() {
        assert_eq!(awk_num("248956422"), 248_956_422.0);
        assert_eq!(awk_num("12ab"), 12.0);
        assert_eq!(awk_num("ab"), 0.0);
        assert_eq!(awk_num(""), 0.0);
        assert_eq!(awk_num("1e3x"), 1000.0);
        assert_eq!(awk_num("."), 0.0);
    }

    #[test]
    fn a_held_stage_goes_when_dropped() {
        let d = std::env::temp_dir().join(format!("og-tools-stage-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        let (guard, st) = stage(&d, ".pypgx", &mut std::io::sink()).unwrap();
        let name = st.file_name().unwrap().to_string_lossy().into_owned();
        assert!(
            name.starts_with(".pypgx.") && name.len() == ".pypgx.".len() + 6,
            "{name}"
        );
        assert!(st.is_dir());
        drop(guard);
        assert!(!st.exists());
        let _ = std::fs::remove_dir_all(&d);
    }
}
