//! Run as root, [`og_io::AtomicFile`] hands what it writes to the owner
//! `NXF_OWNER` names, as Nextflow sets it in a docker task run as root with
//! `docker.fixOwnership` on. The output stays owner-only, and that owner can
//! read it. Without `NXF_OWNER`, or with a value that is not `uid:gid`, the
//! file stays root's, even in a folder another user owns.
//!
//! The test needs root: CI runs it with `sudo`. Run as any other user it
//! prints that it was skipped, unless `OG_IO_ROOT_TEST=1` is set: then it
//! fails, so a run meant to check this cannot pass without checking it.
//!
//! Each write runs in a child process (this test binary again, running only
//! `child_writes`), because the environment of a running test cannot be
//! changed without `unsafe`.

use og_io::AtomicFile;
use std::fs;
use std::io::Write;
use std::os::unix::fs::{MetadataExt, PermissionsExt, chown};
use std::path::Path;
use std::process::Command;

const CHILD: &str = "OG_IO_ROOT_CHILD";

/// Writes the file `OG_IO_ROOT_CHILD` names; does nothing without it.
#[test]
fn child_writes() {
    let Some(target) = std::env::var_os(CHILD) else {
        return;
    };
    let mut f = AtomicFile::create(&target).unwrap();
    f.write_all(b"new\n").unwrap();
    f.commit().unwrap();
}

/// Writes `target` in a child process with `NXF_OWNER` set to `owner`, or
/// unset, and returns the new file's owner and group.
fn write_as_child(target: &Path, owner: Option<&str>) -> (u32, u32) {
    let mut child = Command::new(std::env::current_exe().unwrap());
    child
        .args(["--exact", "child_writes", "--test-threads=1"])
        .env(CHILD, target);
    match owner {
        Some(o) => child.env("NXF_OWNER", o),
        None => child.env_remove("NXF_OWNER"),
    };
    let out = child.output().unwrap();
    assert!(
        out.status.success(),
        "the child failed: {}",
        String::from_utf8_lossy(&out.stdout)
    );
    let m = fs::metadata(target).unwrap();
    assert_eq!(m.permissions().mode() & 0o7777, 0o600, "{target:?}");
    assert_eq!(fs::read_to_string(target).unwrap(), "new\n");
    (m.uid(), m.gid())
}

#[test]
fn root_hands_the_file_to_nxf_owner() {
    let dir = std::env::temp_dir().join(format!("og-io-root-{}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    let probe = dir.join("probe");
    fs::write(&probe, "").unwrap();
    if fs::metadata(&probe).unwrap().uid() != 0 {
        fs::remove_dir_all(&dir).unwrap();
        assert!(
            std::env::var_os("OG_IO_ROOT_TEST").is_none(),
            "OG_IO_ROOT_TEST is set but the test does not run as root"
        );
        eprintln!("skipped: needs root");
        return;
    }

    // NXF_OWNER names the owner: the file is that owner's, still 0600.
    assert_eq!(
        write_as_child(&dir.join("named.tsv"), Some("4321:4322")),
        (4321, 4322),
        "the file NXF_OWNER names an owner for"
    );

    // No NXF_OWNER: root's, even in a folder another user owns.
    let user = dir.join("user");
    fs::create_dir(&user).unwrap();
    chown(&user, Some(4321), Some(4322)).unwrap();
    assert_eq!(
        write_as_child(&user.join("unnamed.tsv"), None),
        (0, 0),
        "a file in another user's folder without NXF_OWNER"
    );

    // An NXF_OWNER that is not uid:gid names no one.
    assert_eq!(
        write_as_child(&user.join("bad.tsv"), Some("alice:users")),
        (0, 0),
        "a file with an NXF_OWNER that is not uid:gid"
    );

    fs::remove_dir_all(&dir).unwrap();
}
