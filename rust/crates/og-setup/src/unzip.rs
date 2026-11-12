//! `og setup unzip ZIP DEST`: every member of ZIP unpacked under DEST, as
//! `zipfile.ZipFile(ZIP).extractall(DEST)` does it: directories made as
//! needed, files written in the order of the central directory (so the last
//! member of a repeated name wins), mode bits and times not restored.
//!
//! Each member is streamed to disk through one buffer, so memory stays
//! flat however large the member. One deliberate difference: a member
//! whose name climbs out of DEST (an absolute path or a `..` part) is
//! refused, exit 2, where `extractall` drops those parts silently and
//! writes the file somewhere else than its name says.

use crate::zip::{Archive, Entry};
use std::fs::{self, File};
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};

/// Where `e` goes under `dest`, or why it may not go anywhere. Empty and
/// `.` parts are dropped, as `extractall` drops them.
pub fn target(dest: &Path, e: &Entry) -> Result<PathBuf, String> {
    if e.name.starts_with('/') {
        return Err(format!("member {:?} has an absolute path", e.name));
    }
    let mut out = dest.to_path_buf();
    for part in e.name.split('/') {
        match part {
            "" | "." => {}
            ".." => return Err(format!("member {:?} climbs out of the folder", e.name)),
            p => out.push(p),
        }
    }
    Ok(out)
}

fn unpack(zip: &str, dest: &str) -> Result<(), (u8, String)> {
    let mut a = Archive::open(zip).map_err(|e| (2, format!("{zip}: {e}")))?;
    let dest = Path::new(dest);
    // Every name is checked before anything is written.
    let plan: Vec<(Entry, PathBuf)> = a
        .entries()
        .iter()
        .map(|e| target(dest, e).map(|t| (e.clone(), t)))
        .collect::<Result<_, _>>()
        .map_err(|m| (2, format!("{zip}: {m}")))?;
    let io = |p: &Path, e: io::Error| (1, format!("{}: {e}", p.display()));
    fs::create_dir_all(dest).map_err(|e| io(dest, e))?;
    for (e, t) in plan {
        if e.is_dir() {
            fs::create_dir_all(&t).map_err(|err| io(&t, err))?;
            continue;
        }
        if let Some(parent) = t.parent() {
            fs::create_dir_all(parent).map_err(|err| io(parent, err))?;
        }
        let mut out =
            BufWriter::with_capacity(1 << 16, File::create(&t).map_err(|err| io(&t, err))?);
        a.copy_to(&e, &mut out).map_err(|err| match err {
            crate::zip::ZipError::Io(x) => io(&t, x),
            crate::zip::ZipError::Bad(m) => (2, format!("{zip}: {m}")),
        })?;
        out.flush().map_err(|err| io(&t, err))?;
    }
    Ok(())
}

/// The command line; returns the exit code: 0, 1 on an I/O error, 2 on a
/// bad zip, a refused member name or a wrong command line.
pub fn cli(args: &[String]) -> u8 {
    let [zip, dest] = args else {
        eprintln!("usage: og setup unzip ZIP DEST");
        return 2;
    };
    match unpack(zip, dest) {
        Ok(()) => 0,
        Err((code, msg)) => {
            eprintln!("og setup unzip: {msg}");
            code
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::zip::tests::build;

    fn scratch(name: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!("og-unzip-{}-{name}", std::process::id()));
        let _ = fs::remove_dir_all(&d);
        fs::create_dir_all(&d).unwrap();
        d
    }

    fn run(dir: &Path, members: &[(&str, &[u8], bool)]) -> u8 {
        fs::write(dir.join("a.zip"), build(members, false)).unwrap();
        cli(&[
            dir.join("a.zip").display().to_string(),
            dir.join("out").display().to_string(),
        ])
    }

    #[test]
    fn unpacks_dirs_files_and_the_last_of_a_name() {
        let d = scratch("ok");
        let rc = run(
            &d,
            &[
                ("top/", b"", false),
                ("top/./x//y.txt", b"deep", true),
                ("top/a.txt", b"one", false),
                ("top/a.txt", b"two", true),
                ("empty/", b"", false),
            ],
        );
        assert_eq!(rc, 0);
        assert_eq!(fs::read(d.join("out/top/x/y.txt")).unwrap(), b"deep");
        assert_eq!(fs::read(d.join("out/top/a.txt")).unwrap(), b"two");
        assert!(d.join("out/empty").is_dir());
    }

    #[test]
    fn a_name_that_climbs_out_is_refused_and_nothing_is_written() {
        for bad in ["../evil.txt", "a/../../evil.txt", "/etc/evil.txt"] {
            let d = scratch("bad");
            let rc = run(&d, &[("fine.txt", b"x", false), (bad, b"y", false)]);
            assert_eq!(rc, 2, "{bad}");
            assert!(!d.join("out").exists(), "{bad}");
            assert!(!d.join("evil.txt").exists(), "{bad}");
        }
    }

    #[test]
    fn a_damaged_member_exits_2() {
        let d = scratch("crc");
        let mut z = build(&[("a.txt", b"some content", false)], false);
        z[30 + 5 + 2] ^= 1;
        fs::write(d.join("a.zip"), z).unwrap();
        let rc = cli(&[
            d.join("a.zip").display().to_string(),
            d.join("out").display().to_string(),
        ]);
        assert_eq!(rc, 2);
    }
}
