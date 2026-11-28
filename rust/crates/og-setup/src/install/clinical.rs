//! The clinical data og installs: `--vep-cache` (the VEP cache of
//! VEP_CACHE_RELEASE, for step 13 and the Nextflow VEP) and
//! `--clinical-data` (that cache, the PCGR/CPSR bundle and the VEP cache of
//! PCGR_VEP_CACHE_RELEASE for step 17). Each does what its bash does
//! (setup.sh's install_vep_cache_opt_in and install_clinical_data,
//! common.sh's install_vep_cache, clinical_disk_ok and
//! install_pcgr_bundle): the same disk check before each item (56, 15 and
//! 52 GB peak plus a 5 GB margin), the same staging folder (`.extract.*`,
//! moved into place only when complete), the same downloads (kept in the
//! downloads folder, through [`Run::fetch`]), the same messages, retry
//! commands and exit status. Unpacking is the host's `tar xzf`, as in bash,
//! so the tree it leaves is the one bash leaves on the same host.
//!
//! [`part_cli`] answers `og setup install --part`, the three common.sh
//! helpers above as og runs them: step 13 installs its cache at run time
//! through install_vep_cache, and a full setup's bash calls all three.

use super::keep::{self, Check, exists, is_file};
use super::{Missing, Run, Vars, shell_quote};
use crate::disk;
use crate::fetch::{self, Sources};
use std::io::Write;
use std::path::Path;
use std::process::{Command, Stdio};

/// The keys this module runs.
pub const PORTED: &[&str] = &["vep-cache", "clinical-data"];

/// Runs KEY when it is one of [`PORTED`]; its exit status.
pub fn run(key: &str, r: &Run, out: &mut dyn Write, err: &mut dyn Write) -> Option<u8> {
    Some(match key {
        "vep-cache" => vep_cache_opt_in(r, out, err),
        // Only a full setup lists what is missing; one flag exits with 1.
        "clinical-data" => clinical_data(r, &mut Missing::default(), out, err),
        _ => return None,
    })
}

fn is_dir(p: &str) -> bool {
    Path::new(p).is_dir()
}

/// One field of a row of assets/data-sources.tsv (`source_of ID FIELD`:
/// path, url or value); its error goes to ERR.
pub fn source_field(r: &Run, id: &str, field: &str, err: &mut dyn Write) -> Option<String> {
    if field != "value" {
        return r.source_of(id, field, err);
    }
    let v = (|| -> Result<String, String> {
        let ctx = r.ctx()?;
        let s = Sources::read(&ctx.root)?;
        ctx.value(&s, s.get(id)?, None)
    })();
    v.map_err(|m| {
        let _ = writeln!(err, "{m}");
    })
    .ok()
}

/// vep_cache_row RELEASE: the table row of the VEP cache of RELEASE.
pub fn vep_cache_row(r: &Run, rel: &str, err: &mut dyn Write) -> Option<&'static str> {
    if rel == r.version("VEP_CACHE_RELEASE") {
        Some("vep_cache")
    } else if rel == r.version("PCGR_VEP_CACHE_RELEASE") {
        Some("pcgr_vep_cache")
    } else {
        let _ = writeln!(
            err,
            "ERROR: VEP cache {rel} is neither VEP_CACHE_RELEASE nor PCGR_VEP_CACHE_RELEASE"
        );
        None
    }
}

/// The checksum of the PCGR bundle: PCGR_DATA_SHA256 from the environment
/// when set (setup.sh's override, for fixtures and verified mirrors), else
/// versions.env's pin.
pub fn pcgr_sha256(r: &Run) -> String {
    r.vars
        .nonempty("PCGR_DATA_SHA256")
        .map(str::to_owned)
        .unwrap_or_else(|| r.version("PCGR_DATA_SHA256"))
}

/// clinical_disk_ok DIR PEAK LABEL [ARCHIVE] with the shell's df and du.
pub fn disk_ok(
    r: &Run,
    dir: &str,
    peak: u64,
    label: &str,
    archive: &str,
    out: &mut dyn Write,
) -> bool {
    let a = (!archive.is_empty()).then(|| Path::new(archive));
    disk::clinical_disk_ok(&disk::Shell, Path::new(dir), peak, label, a, &r.genome, out)
}

/// `rm -rf DIR/.extract.*`: the staging trees of stopped extractions.
fn clear_staging(dir: &str, out: &mut dyn Write) {
    let Ok(rd) = std::fs::read_dir(dir) else {
        return;
    };
    for e in rd.flatten() {
        if e.file_name().to_string_lossy().starts_with(".extract.") {
            keep::rm_rf(&format!("{dir}/{}", e.file_name().to_string_lossy()), out);
        }
    }
}

/// `mktemp -d DIR/.extract.XXXXXX`: a new folder only its owner can enter.
pub fn staging(dir: &str, err: &mut dyn Write) -> Option<String> {
    use std::hash::{BuildHasher, Hasher};
    use std::os::unix::fs::DirBuilderExt;
    const ABC: &[u8] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";
    let mut last = None;
    for _ in 0..100 {
        let mut h = std::collections::hash_map::RandomState::new().build_hasher();
        h.write_u32(std::process::id());
        let mut n = h.finish();
        let suffix: String = (0..6)
            .map(|_| {
                let c = ABC[(n % ABC.len() as u64) as usize] as char;
                n /= ABC.len() as u64;
                c
            })
            .collect();
        let p = format!("{dir}/.extract.{suffix}");
        match std::fs::DirBuilder::new().mode(0o700).create(&p) {
            Ok(()) => return Some(p),
            Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => {}
            Err(e) => {
                last = Some(e);
                break;
            }
        }
    }
    let _ = writeln!(
        err,
        "mktemp: failed to create directory via template '{dir}/.extract.XXXXXX': {}",
        last.map_or_else(|| "File exists".to_owned(), |e| e.to_string())
    );
    None
}

/// `tar xzf ARCHIVE -C DIR`, the host's tar with its own messages.
pub fn untar(archive: &str, dir: &str, out: &mut dyn Write, err: &mut dyn Write) -> bool {
    let _ = out.flush();
    match crate::child::status(
        Command::new("tar")
            .args(["xzf", archive, "-C", dir])
            .stdin(Stdio::null()),
    ) {
        Ok(s) => s.success(),
        Err(e) => {
            let _ = writeln!(err, "ERROR: cannot run tar: {e}");
            false
        }
    }
}

/// install_vep_cache DIR RELEASE: the cache of RELEASE, from its table row.
pub fn install_vep_cache(
    r: &Run,
    dir: &str,
    rel: &str,
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> bool {
    let Some(row) = vep_cache_row(r, rel, err) else {
        return false;
    };
    let (Some(url), Some(tarball), Some(sums)) = (
        source_field(r, row, "url", err),
        source_field(r, row, "path", err),
        source_field(r, row, "value", err),
    ) else {
        return false;
    };
    install_cache_from(r, dir, rel, &url, &tarball, &sums, out, err)
}

fn base(p: &str) -> &str {
    p.rsplit('/').next().unwrap_or(p)
}

/// install_vep_cache with its download given: URL to TARBALL, checked
/// against the checksum list SUMS, unpacked into a staging folder and moved
/// to DIR/homo_sapiens/RELEASE_GRCh38 only when it holds that cache.
#[allow(clippy::too_many_arguments)]
pub fn install_cache_from(
    r: &Run,
    dir: &str,
    rel: &str,
    url: &str,
    tarball: &str,
    sums: &str,
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> bool {
    let fin = format!("{dir}/homo_sapiens/{rel}_GRCh38");
    if is_file(&format!("{fin}/info.txt")) {
        let _ = writeln!(out, "[OK] VEP {rel} cache already present.");
        return true;
    }
    if exists(&fin) {
        let _ = writeln!(
            err,
            "ERROR: {fin} exists but has no info.txt (an interrupted extraction?)."
        );
        let _ = writeln!(err, "  Remove it and run again.");
        return false;
    }
    if !keep::mkdir_p(&format!("{dir}/homo_sapiens"), out) {
        return false;
    }
    let check = Check {
        url,
        sum: Some(("sum", sums)),
    };
    keep::adopt(
        r,
        &format!("{dir}/{}", base(tarball)),
        tarball,
        Some(&check),
        out,
    );
    // fetch uses a kept tarball that matches its CHECKSUMS line with no
    // download (it reads that small list first), so an install stopped
    // during extraction unpacks it again.
    if !r.fetch(url, Path::new(tarball), "sum", sums, err) {
        return false;
    }
    let _ = writeln!(out, "  Extracting {}...", base(tarball));
    clear_staging(dir, out);
    let Some(tmp) = staging(dir, err) else {
        return false;
    };
    let unpacked = format!("{tmp}/homo_sapiens/{rel}_GRCh38");
    if !untar(tarball, &tmp, out, err) || !is_file(&format!("{unpacked}/info.txt")) {
        keep::rm_rf(&tmp, out);
        let _ = writeln!(
            err,
            "ERROR: could not unpack {tarball} into a {rel}_GRCh38 cache."
        );
        return false;
    }
    if !keep::mv(&unpacked, &fin, out) {
        keep::rm_rf(&tmp, out);
        return false;
    }
    keep::rm_rf(&tmp, out);
    keep::done_with(r, &[tarball], out)
}

/// install_pcgr_bundle DIR: the PCGR/CPSR bundle PCGR_DATA_BUNDLE, checked
/// against its sha256 before it is unpacked into a staging folder, which
/// becomes DIR/BUNDLE (upstream ships data/ at the root; CPSR reads
/// BUNDLE/data/).
pub fn install_pcgr_bundle(r: &Run, dir: &str, out: &mut dyn Write, err: &mut dyn Write) -> bool {
    let bundle = r.version("PCGR_DATA_BUNDLE");
    if is_dir(&format!("{dir}/{bundle}/data")) {
        let _ = writeln!(out, "[OK] PCGR/CPSR {bundle} bundle already present.");
        return true;
    }
    let _ = writeln!(
        out,
        "PCGR/CPSR reference data includes third-party databases with their own terms:"
    );
    let _ = writeln!(
        out,
        "  https://sigven.github.io/pcgr/articles/annotation_resources.html"
    );
    let _ = writeln!(
        out,
        "  Review those terms before use; --no-clinical-data skips this download."
    );
    let (Some(url), Some(tarball)) = (
        source_field(r, "pcgr_bundle", "url", err),
        source_field(r, "pcgr_bundle", "path", err),
    ) else {
        return false;
    };
    install_bundle_from(r, dir, &bundle, &url, &tarball, &pcgr_sha256(r), out, err)
}

/// install_pcgr_bundle with its download given (URL to TARBALL, SHA256).
#[allow(clippy::too_many_arguments)]
pub fn install_bundle_from(
    r: &Run,
    dir: &str,
    bundle: &str,
    url: &str,
    tarball: &str,
    sha256: &str,
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> bool {
    if !keep::mkdir_p(dir, out) {
        return false;
    }
    let check = Check {
        url,
        sum: Some(("sha256", sha256)),
    };
    keep::adopt(
        r,
        &format!("{dir}/{}", base(tarball)),
        tarball,
        Some(&check),
        out,
    );
    if !r.fetch(url, Path::new(tarball), "sha256", sha256, err) {
        return false;
    }
    clear_staging(dir, out);
    let Some(tmp) = staging(dir, err) else {
        return false;
    };
    if !untar(tarball, &tmp, out, err) || !is_dir(&format!("{tmp}/data/grch38")) {
        keep::rm_rf(&tmp, out);
        let _ = writeln!(
            err,
            "ERROR: could not unpack {tarball} into PCGR data/grch38."
        );
        return false;
    }
    // An incomplete destination is never mistaken for a completed bundle.
    let dest = format!("{dir}/{bundle}");
    if exists(&dest) {
        keep::rm_rf(&tmp, out);
        let _ = writeln!(
            err,
            "ERROR: {dest} exists without data; remove it and retry."
        );
        return false;
    }
    if !keep::mv(&tmp, &dest, out) {
        keep::rm_rf(&tmp, out);
        return false;
    }
    keep::done_with(r, &[tarball], out);
    let _ = writeln!(out, "[OK] PCGR/CPSR {bundle} bundle installed.");
    true
}

/// install_vep_cache_opt_in: `setup.sh --vep-cache`, the cache of
/// VEP_CACHE_RELEASE after its disk check (56 GB peak).
pub fn vep_cache_opt_in(r: &Run, out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let rel = r.version("VEP_CACHE_RELEASE");
    let dir = format!("{}/vep_cache", r.genome);
    let fin = format!("{dir}/homo_sapiens/{rel}_GRCh38");
    if is_file(&format!("{fin}/info.txt")) {
        let _ = writeln!(out, "[OK] VEP {rel} cache already present: {fin}");
        return 0;
    }
    let archive = source_field(r, "vep_cache", "path", err).unwrap_or_default();
    if !disk_ok(r, &dir, 56, &format!("VEP {rel} cache"), &archive, out) {
        return 1;
    }
    let _ = writeln!(
        out,
        "Installing the VEP {rel} cache: ~26 GB download, ~30 GB unpacked, both on disk until it ends."
    );
    if !install_vep_cache(r, &dir, &rel, out, err) {
        let _ = writeln!(
            out,
            "[WARN] Could not install the VEP {rel} cache. Run this command again: the download resumes."
        );
        return 1;
    }
    let _ = writeln!(out, "[OK] VEP {rel} cache: {fin}");
    0
}

/// install_clinical_data: `setup.sh --clinical-data`, and the clinical part
/// of a full setup. Each item gets its own disk check; a skipped or failed
/// one does not stop the others, makes the status 1 and is noted in MISSING
/// with the command that retries it.
pub fn clinical_data(
    r: &Run,
    missing: &mut Missing,
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> u8 {
    let g = &r.genome;
    let dir = format!("{g}/vep_cache");
    let vep = r.version("VEP_CACHE_RELEASE");
    let pcgr_vep = r.version("PCGR_VEP_CACHE_RELEASE");
    let bundle = r.version("PCGR_DATA_BUNDLE");
    let setup_cmd = shell_quote(&r.script);
    let mut rc = 0;
    let _ = writeln!(
        out,
        "Clinical data: about 57 GB to download, 65 GB unpacked; allow 1 to 3 hours."
    );
    for rel in [vep.as_str(), "bundle", pcgr_vep.as_str()] {
        if r.vars.get("SETUP_EVERYTHING_SCOPE") == Some("row") && rel == vep {
            continue;
        }
        if rel == "bundle" {
            let data = format!("{g}/pcgr_data/{bundle}/data");
            if is_dir(&data) {
                let _ = writeln!(out, "[OK] PCGR/CPSR bundle already present.");
                continue;
            }
            let archive = source_field(r, "pcgr_bundle", "path", err).unwrap_or_default();
            let pdir = format!("{g}/pcgr_data");
            if disk_ok(r, &pdir, 15, "PCGR/CPSR bundle", &archive, out) {
                if !install_pcgr_bundle(r, &pdir, out, err) {
                    rc = 1;
                    let _ = writeln!(
                        out,
                        "[WARN] Could not install PCGR/CPSR. Step 17 needs it; run: ./scripts/setup.sh --clinical-data \"{g}\""
                    );
                }
            } else {
                rc = 1;
            }
            if !is_dir(&data) {
                missing.note(
                    &format!("PCGR/CPSR bundle {bundle} (step 17)"),
                    &format!("{setup_cmd} --clinical-data \"{g}\""),
                );
            }
            continue;
        }
        let info = format!("{dir}/homo_sapiens/{rel}_GRCh38/info.txt");
        if is_file(&info) {
            let _ = writeln!(out, "[OK] VEP {rel} cache already present.");
            continue;
        }
        let peak = if rel == pcgr_vep { 52 } else { 56 };
        let row = vep_cache_row(r, rel, err).unwrap_or("");
        let archive = source_field(r, row, "path", err).unwrap_or_default();
        if disk_ok(r, &dir, peak, &format!("VEP {rel} cache"), &archive, out) {
            if !install_vep_cache(r, &dir, rel, out, err) {
                rc = 1;
                let _ = writeln!(
                    out,
                    "[WARN] Could not install VEP {rel}. Run: ./scripts/setup.sh --clinical-data \"{g}\" (interrupted downloads resume)."
                );
                if rel == vep {
                    let _ = writeln!(
                        out,
                        "  To retry only this cache: ./scripts/setup.sh --vep-cache \"{g}\""
                    );
                }
            }
        } else {
            rc = 1;
        }
        if !is_file(&info) {
            if rel == vep {
                missing.note(
                    &format!("VEP {rel} cache (steps 13, 30, 23 and 31)"),
                    &format!("{setup_cmd} --vep-cache \"{g}\""),
                );
            } else {
                missing.note(
                    &format!("VEP {rel} cache (step 17, CPSR)"),
                    &format!("{setup_cmd} --clinical-data \"{g}\""),
                );
            }
        }
    }
    rc
}

const PART_USAGE: &str = "usage: og setup install --part vep-cache DIR RELEASE
       og setup install --part disk-ok DIR PEAK_GB LABEL [ARCHIVE]
       og setup install --part pcgr-bundle DIR";

/// `og setup install --part NAME ARGS`: common.sh's install_vep_cache
/// (`vep-cache`), clinical_disk_ok (`disk-ok`) and install_pcgr_bundle
/// (`pcgr-bundle`), run by og for the bash that calls them (step 13, and a
/// full setup's clinical part). GENOME_DIR, OG_DOWNLOADS and
/// SETUP_KEEP_DOWNLOADS come from the environment, as common.sh's do.
/// `None` when ARGS is not a `--part` call.
pub fn part_cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> Option<u8> {
    if args.first().map(String::as_str) != Some("--part") {
        return None;
    }
    let a: Vec<&str> = args[1..].iter().map(String::as_str).collect();
    let shape_ok = matches!(
        a.as_slice(),
        ["vep-cache", _, _] | ["disk-ok", _, _, _] | ["disk-ok", _, _, _, _] | ["pcgr-bundle", _]
    );
    if !shape_ok {
        let _ = writeln!(err, "{PART_USAGE}");
        return Some(2);
    }
    let root = fetch::root_of(None);
    let vars = Vars::from_process();
    let script = vars
        .nonempty("OG_SETUP_SCRIPT")
        .map(str::to_owned)
        .unwrap_or_else(|| root.join("scripts/setup.sh").to_string_lossy().into_owned());
    let genome = vars.nonempty("GENOME_DIR").unwrap_or("").to_owned();
    let keep = vars.get("SETUP_KEEP_DOWNLOADS") != Some("false");
    let r = match Run::new(&root, vars, &script, &genome, keep, err) {
        Ok(r) => r,
        Err(code) => return Some(code),
    };
    let ok = match a.as_slice() {
        ["vep-cache", dir, rel] => install_vep_cache(&r, dir, rel, out, err),
        ["disk-ok", dir, peak, label, rest @ ..] => {
            let Some(peak) = peak
                .bytes()
                .all(|b| b.is_ascii_digit())
                .then(|| peak.parse::<u64>().ok())
                .flatten()
            else {
                let _ = writeln!(err, "ERROR: PEAK_GB must be a whole number, got '{peak}'");
                return Some(2);
            };
            disk_ok(
                &r,
                dir,
                peak,
                label,
                rest.first().copied().unwrap_or(""),
                out,
            )
        }
        ["pcgr-bundle", dir] => install_pcgr_bundle(&r, dir, out, err),
        _ => unreachable!("the shape was checked above"),
    };
    Some(u8::from(!ok))
}
