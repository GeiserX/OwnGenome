//! `og setup plan`: before anything downloads, every source the chosen
//! installers need (`assets/data-sources.tsv`), whether it is already on
//! disk, its size, what is left to download and the free disk.
//!
//! Presence is a file at the row's path; `og setup fetch` checks its
//! checksum when it is reached, so the plan itself reads no file and makes
//! no request.

use crate::fetch::{Ctx, Source, Sources, flag, root_of, take};
use std::io::Write;
use std::path::{Path, PathBuf};

/// A size as people read it: kB, MB, GB, TB (powers of 1000).
pub fn human(b: u64) -> String {
    let units = ["bytes", "kB", "MB", "GB", "TB"];
    let mut v = b as f64;
    let mut u = 0;
    while v >= 1000.0 && u + 1 < units.len() {
        v /= 1000.0;
        u += 1;
    }
    if u == 0 {
        format!("{b} bytes")
    } else {
        format!("{v:.1} {}", units[u])
    }
}

/// The rows of INSTALLERS (any row naming one), in table order.
pub fn chosen<'a>(s: &'a Sources, installers: &[String]) -> Vec<&'a Source> {
    s.rows
        .values()
        .filter(|r| r.installers.iter().any(|i| installers.contains(i)))
        .collect()
}

/// The free bytes of the filesystem that holds PATH (its nearest existing
/// folder).
pub fn free_bytes(path: &Path) -> Option<(PathBuf, u64)> {
    let mut p = path.to_path_buf();
    while !p.exists() {
        p = p.parent()?.to_path_buf();
    }
    let st = rustix::fs::statvfs(&p).ok()?;
    Some((p, st.f_bavail.saturating_mul(st.f_frsize)))
}

/// One row of the plan.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Line {
    pub id: String,
    pub path: String,
    pub files: usize,
    pub present: usize,
    pub bytes: Option<u64>,
    /// Kept in the downloads folder (a raw download), not installed.
    pub download: bool,
}

/// The plan of ROWS: one line each, a family counted file by file.
pub fn lines(ctx: &Ctx, rows: &[&Source]) -> Result<Vec<Line>, String> {
    let mut out = vec![];
    for r in rows {
        let parts: Vec<Option<&str>> = if r.parts.is_empty() {
            vec![None]
        } else {
            r.parts.iter().map(|p| Some(p.as_str())).collect()
        };
        let mut present = 0;
        for p in &parts {
            let f = ctx.path(r, *p)?;
            // A folder row (a path ending in /) is there when it holds a file.
            let there = if r.path.ends_with('/') {
                has_file(&f)
            } else {
                f.is_file()
            };
            if there {
                present += 1;
            }
        }
        out.push(Line {
            id: r.id.clone(),
            path: ctx
                .expand(&r.path, None)
                .map_err(|e| format!("ERROR: {}: {e}", r.id))?,
            files: parts.len(),
            present,
            bytes: r.bytes,
            download: r.is_download(),
        });
    }
    Ok(out)
}

/// True when folder D holds a file.
pub fn has_file(d: &Path) -> bool {
    std::fs::read_dir(d).is_ok_and(|mut e| e.any(|x| x.is_ok_and(|x| x.path().is_file())))
}

/// How much NEED is over FREE, if it is.
pub fn shortfall(need: u64, free: u64) -> Option<u64> {
    (need > free).then(|| need - free)
}

/// True unless SETUP_KEEP_DOWNLOADS=false: raw downloads stay after the install.
pub fn keeping(ctx: &Ctx) -> bool {
    ctx.setting("SETUP_KEEP_DOWNLOADS").as_deref() != Some("false")
}

/// Writes the plan; returns the bytes still to download into the downloads
/// folder and that folder's free bytes (when readable).
pub fn write_plan(
    ctx: &Ctx,
    lines: &[Line],
    out: &mut dyn Write,
) -> Result<(u64, Option<u64>), String> {
    let downloads = ctx.downloads()?;
    let genome = ctx.genome_dir().ok();
    let w = |out: &mut dyn Write, s: String| out.write_all(s.as_bytes()).map_err(|e| e.to_string());
    w(out, format!("Downloads folder: {}\n", downloads.display()))?;
    w(out, "source\tfile\tsize\tstate\n".to_owned())?;
    let (mut need, mut have, mut unknown, mut n_need, mut kept, mut need_dl) =
        (0u64, 0u64, 0usize, 0usize, 0u64, 0u64);
    for l in lines {
        let size = l.bytes.map_or_else(|| "?".to_owned(), human);
        let state = if l.files == 1 {
            if l.present == 1 {
                "present".to_owned()
            } else {
                "missing".to_owned()
            }
        } else {
            format!("{} of {} files present", l.present, l.files)
        };
        w(out, format!("{}\t{}\t{size}\t{state}\n", l.id, l.path))?;
        let missing = l.files - l.present;
        if missing > 0 {
            n_need += 1;
            match l.bytes {
                // A family's size is shared out over its files.
                Some(b) => {
                    let n = b / l.files as u64 * missing as u64;
                    need += n;
                    if l.download {
                        need_dl += n;
                    }
                }
                None => unknown += 1,
            }
        }
        if l.download {
            kept += l.bytes.unwrap_or(0);
        }
        if l.present > 0
            && let Some(b) = l.bytes
        {
            have += b / l.files as u64 * l.present as u64;
        }
    }
    let unk = if unknown > 0 {
        format!(" ({unknown} of unknown size)")
    } else {
        String::new()
    };
    w(
        out,
        format!(
            "To download: {n_need} source{}, {}{unk}; already on disk: {}.\n",
            if n_need == 1 { "" } else { "s" },
            human(need),
            human(have)
        ),
    )?;
    if keeping(ctx) {
        w(
            out,
            format!(
                "Raw downloads kept after the install: {} (--no-keep-downloads or SETUP_KEEP_DOWNLOADS=false deletes each once installed; og setup prune deletes them later).\n",
                human(kept)
            ),
        )?;
    } else {
        w(
            out,
            "Raw downloads are deleted once installed (SETUP_KEEP_DOWNLOADS=false).\n".to_owned(),
        )?;
    }
    let mut dl_free = None;
    let mut seen: Vec<PathBuf> = vec![];
    for (i, d) in std::iter::once(downloads.clone()).chain(genome).enumerate() {
        let Some((at, free)) = free_bytes(&d) else {
            w(
                out,
                format!("Free disk: cannot read it for {}\n", d.display()),
            )?;
            continue;
        };
        if seen.iter().any(|s| same_fs(s, &at)) {
            continue;
        }
        w(
            out,
            format!("Free disk: {} in {}\n", human(free), at.display()),
        )?;
        // Only the downloads folder's own filesystem counts for --check.
        if i == 0 {
            dl_free = Some(free);
            if let Some(s) = shortfall(need_dl, free) {
                w(
                    out,
                    format!(
                        "[WARN] the raw downloads kept there need {} more than is free.\n",
                        human(s)
                    ),
                )?;
            }
        }
        seen.push(at);
    }
    Ok((need_dl, dl_free))
}

fn same_fs(a: &Path, b: &Path) -> bool {
    use std::os::unix::fs::MetadataExt;
    match (std::fs::metadata(a), std::fs::metadata(b)) {
        (Ok(x), Ok(y)) => x.dev() == y.dev(),
        _ => false,
    }
}

const USAGE: &str =
    "usage: og setup plan [--for INSTALLER,...|--all] [--skip INSTALLER,...] [--check] [--root DIR]
  --check: exit 1 when the raw downloads still to come do not fit in the
  free disk of the downloads folder (unless SETUP_KEEP_DOWNLOADS=false);
  an installer with no source then counts as nothing to download
  --skip: leave out the sources of these installers (setup.sh
  --no-clinical-data is --skip clinical-data)
  INSTALLER: setup (a plain setup), or a setup.sh flag without its dashes
  (clinical-data, vep-cache, archaic, ...); default setup
";

/// `og setup plan`.
pub fn cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let mut a = args.to_vec();
    let all = flag(&mut a, "--all");
    let check = flag(&mut a, "--check");
    let (list, root, skip) = match (
        take(&mut a, "--for"),
        take(&mut a, "--root"),
        take(&mut a, "--skip"),
    ) {
        (Ok(l), Ok(r), Ok(s)) => (l, r, s),
        (Err(m), _, _) | (_, Err(m), _) | (_, _, Err(m)) => {
            let _ = writeln!(err, "og setup plan: {m}");
            return 2;
        }
    };
    if !a.is_empty() {
        let _ = write!(err, "{USAGE}");
        return 2;
    }
    let res = (|| -> Result<(), String> {
        let ctx = Ctx::from_env(&root_of(root.as_deref()))?;
        let s = Sources::read(&ctx.root)?;
        let skip: Vec<&str> = skip
            .as_deref()
            .unwrap_or("")
            .split(',')
            .filter(|x| !x.is_empty())
            .collect();
        let rows: Vec<&Source> = if all {
            s.rows.values().collect()
        } else {
            let want: Vec<String> = list
                .as_deref()
                .unwrap_or("setup")
                .split(',')
                .filter(|x| !x.is_empty())
                .map(str::to_owned)
                .collect();
            for i in &want {
                if !check && !s.rows.values().any(|r| r.installers.contains(i)) {
                    return Err(format!(
                        "ERROR: no source of the table is installed by '{i}'"
                    ));
                }
            }
            chosen(&s, &want)
        };
        let rows: Vec<&Source> = rows
            .into_iter()
            .filter(|r| !r.installers.iter().any(|i| skip.contains(&i.as_str())))
            .collect();
        let l = lines(&ctx, &rows)?;
        let (need, free) = write_plan(&ctx, &l, out)?;
        if check
            && keeping(&ctx)
            && let Some(s) = free.and_then(|f| shortfall(need, f))
        {
            return Err(format!(
                "ERROR: the raw downloads of this setup need {} more than is free in {}: free some disk, set OG_DOWNLOADS to a larger disk, or run with --no-keep-downloads. The check counts every raw download that is not there, even one whose install is done (after og setup prune): SETUP_DISK_CHECK=false skips it",
                human(s),
                ctx.downloads()?.display()
            ));
        }
        Ok(())
    })();
    match res {
        Ok(()) => 0,
        Err(m) => {
            let _ = writeln!(err, "{m}");
            1
        }
    }
}

#[cfg(test)]
mod tests {
    use super::human;

    #[test]
    fn sizes_read_as_people_read_them() {
        assert_eq!(human(999), "999 bytes");
        assert_eq!(human(1_500), "1.5 kB");
        assert_eq!(human(26_000_000_000), "26.0 GB");
    }
}
