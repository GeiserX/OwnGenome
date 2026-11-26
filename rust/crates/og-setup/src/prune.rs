//! `og setup prune`: delete the raw downloads setup keeps.
//!
//! Setup keeps every raw download in the downloads folder (`OG_DOWNLOADS`,
//! or `$GENOME_DIR/downloads`), so a rerun, a second install or another use
//! of the same data does not download it again. Prune removes them on
//! request: the files of the table's rows under `downloads/` (and their
//! `.part`), all of them or the rows named. It never removes a file a step
//! reads (a row stored elsewhere in `GENOME_DIR`), and leaves any other file
//! in the downloads folder alone.

use crate::fetch::{Ctx, Source, Sources, flag, root_of, take};
use crate::plan::human;
use std::fs;
use std::io::ErrorKind;
use std::io::Write;
use std::path::{Path, PathBuf};

/// True when the folder that holds FILE resolves inside DL (the resolved
/// downloads folder). The leaf itself is not followed: a symlink there is
/// unlinked, never its target.
fn inside(file: &Path, dl: &Path) -> bool {
    file.parent()
        .and_then(|d| d.canonicalize().ok())
        .is_some_and(|d| d.starts_with(dl))
}

/// The rows prune works on: every download of the table, or the IDS named,
/// each of which must be a download (a row stored where a step reads it is
/// refused).
pub fn rows<'a>(s: &'a Sources, ids: &[String]) -> Result<Vec<&'a Source>, String> {
    if ids.is_empty() {
        return Ok(s.rows.values().filter(|r| r.is_download()).collect());
    }
    let mut v = vec![];
    for id in ids {
        let r = s.get(id)?;
        if !r.is_download() {
            return Err(format!(
                "ERROR: {id} is stored at {}, where a step reads it; prune removes only the downloads folder",
                r.path
            ));
        }
        v.push(r);
    }
    Ok(v)
}

/// The files prune would remove for ROWS: each existing download and
/// `.part`, with its size.
pub fn targets(ctx: &Ctx, rows: &[&Source]) -> Result<Vec<(PathBuf, u64)>, String> {
    let mut out = vec![];
    let dl = ctx.downloads()?;
    let dl = match dl.canonicalize() {
        Ok(d) => d,
        Err(e) if e.kind() == ErrorKind::NotFound => return Ok(out),
        Err(e) => return Err(format!("ERROR: cannot read {}: {e}", dl.display())),
    };
    let outside = |c: &Path| {
        format!(
            "ERROR: {} resolves outside the downloads folder {}; nothing was removed",
            c.display(),
            dl.display()
        )
    };
    for r in rows {
        let parts: Vec<Option<&str>> = if r.parts.is_empty() {
            vec![None]
        } else {
            r.parts.iter().map(|p| Some(p.as_str())).collect()
        };
        for p in parts {
            let f = ctx.path(r, p)?;
            if r.path.ends_with('/') {
                // A folder row: the files in it.
                let entries = match fs::read_dir(&f) {
                    Ok(x) => x,
                    Err(e) if e.kind() == ErrorKind::NotFound => continue,
                    Err(e) => return Err(format!("ERROR: cannot read {}: {e}", f.display())),
                };
                for e in entries {
                    let e = e.map_err(|e| format!("ERROR: cannot read {}: {e}", f.display()))?;
                    let m = e
                        .metadata()
                        .map_err(|x| format!("ERROR: cannot read {}: {x}", e.path().display()))?;
                    if m.is_file() {
                        if !inside(&e.path(), &dl) {
                            return Err(outside(&e.path()));
                        }
                        out.push((e.path(), m.len()));
                    }
                }
                continue;
            }
            for c in [f.clone(), PathBuf::from(format!("{}.part", f.display()))] {
                let m = match fs::symlink_metadata(&c) {
                    Ok(m) => m,
                    Err(e) if e.kind() == ErrorKind::NotFound => continue,
                    Err(e) => return Err(format!("ERROR: cannot read {}: {e}", c.display())),
                };
                if m.is_file() || m.file_type().is_symlink() {
                    if !inside(&c, &dl) {
                        return Err(outside(&c));
                    }
                    out.push((c, m.len()));
                }
            }
        }
    }
    let mut seen = std::collections::HashSet::new();
    out.retain(|(p, _)| seen.insert(p.clone()));
    Ok(out)
}

const USAGE: &str = "usage: og setup prune [ID...] [--dry-run] [--root DIR]
  removes the raw downloads of the table (all, or the sources named) from
  the downloads folder; --dry-run only lists them
";

/// `og setup prune`.
pub fn cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let mut a = args.to_vec();
    let dry = flag(&mut a, "--dry-run");
    let root = match take(&mut a, "--root") {
        Ok(r) => r,
        Err(m) => {
            let _ = writeln!(err, "og setup prune: {m}");
            return 2;
        }
    };
    if a.iter().any(|x| x.starts_with('-')) {
        let _ = write!(err, "{USAGE}");
        return 2;
    }
    let res = (|| -> Result<(), String> {
        let ctx = Ctx::from_env(&root_of(root.as_deref()))?;
        let s = Sources::read(&ctx.root)?;
        let rows = rows(&s, &a)?;
        let t = targets(&ctx, &rows)?;
        let mut freed = 0u64;
        for (p, n) in &t {
            if dry {
                writeln!(out, "would remove {} ({})", p.display(), human(*n))
                    .map_err(|e| e.to_string())?;
            } else {
                fs::remove_file(p)
                    .map_err(|e| format!("ERROR: cannot remove {}: {e}", p.display()))?;
                writeln!(out, "removed {} ({})", p.display(), human(*n))
                    .map_err(|e| e.to_string())?;
            }
            freed += n;
        }
        let verb = if dry { "would free" } else { "freed" };
        writeln!(
            out,
            "{} file{}, {verb} {}.",
            t.len(),
            if t.len() == 1 { "" } else { "s" },
            human(freed)
        )
        .map_err(|e| e.to_string())
    })();
    match res {
        Ok(()) => 0,
        Err(m) => {
            let _ = writeln!(err, "{m}");
            1
        }
    }
}
