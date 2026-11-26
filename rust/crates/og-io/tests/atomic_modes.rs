//! The modes and temporary names of [`og_io::AtomicFile`]: what it writes is
//! owner-only whatever the umask, and two writers in one process never share
//! a temporary file.

use og_io::AtomicFile;
use std::fs;
use std::io::Write;
use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};
use std::process::Command;

const CHILD: &str = "OG_IO_UMASK_022_CHILD";

fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("og-io-{name}-{}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    dir
}

/// The permission bits in octal, so a failure reads `644`, not `420`.
fn mode(path: &Path) -> String {
    format!(
        "{:o}",
        fs::metadata(path).unwrap().permissions().mode() & 0o7777
    )
}

fn write(target: &Path, text: &str) {
    let mut f = AtomicFile::create(target).unwrap();
    f.write_all(text.as_bytes()).unwrap();
    f.commit().unwrap();
}

/// Runs itself again in a child under `umask 022`, the umask most shells
/// start with, and checks that the child ran this test and passed it.
#[test]
fn owner_only_under_umask_022() {
    if std::env::var_os(CHILD).is_none() {
        let exe = std::env::current_exe().unwrap();
        let out = Command::new("/bin/sh")
            .arg("-c")
            .arg("umask 022 && exec \"$0\" --exact owner_only_under_umask_022 --nocapture --test-threads=1")
            .arg(&exe)
            .env(CHILD, "1")
            .output()
            .unwrap();
        let stdout = String::from_utf8_lossy(&out.stdout);
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert!(
            out.status.success() && stdout.contains("1 passed"),
            "the umask 022 child failed or ran no test:\n{stdout}\n{stderr}"
        );
        return;
    }
    let dir = scratch("umask");
    // The control: under this umask a plain new file is 0644.
    let control = dir.join("control");
    fs::File::create(&control).unwrap();
    assert_eq!(
        mode(&control),
        "644",
        "the child does not run under umask 022"
    );

    // A new target.
    let new = dir.join("new.tsv");
    write(&new, "new\n");
    assert_eq!(mode(&new), "600", "a new target");

    // A 0600 target that is replaced stays 0600.
    let private = dir.join("private.tsv");
    fs::write(&private, "old\n").unwrap();
    fs::set_permissions(&private, fs::Permissions::from_mode(0o600)).unwrap();
    write(&private, "new\n");
    assert_eq!(mode(&private), "600", "a replaced 0600 target");
    assert_eq!(fs::read_to_string(&private).unwrap(), "new\n");

    // A 0644 target becomes owner-only.
    let open = dir.join("open.tsv");
    fs::write(&open, "old\n").unwrap();
    fs::set_permissions(&open, fs::Permissions::from_mode(0o644)).unwrap();
    write(&open, "new\n");
    assert_eq!(mode(&open), "600", "a replaced 0644 target");

    // A stricter target keeps its mode.
    let read_only = dir.join("read-only.tsv");
    fs::write(&read_only, "old\n").unwrap();
    fs::set_permissions(&read_only, fs::Permissions::from_mode(0o400)).unwrap();
    write(&read_only, "new\n");
    assert_eq!(mode(&read_only), "400", "a replaced 0400 target");
    assert_eq!(fs::read_to_string(&read_only).unwrap(), "new\n");

    fs::remove_dir_all(&dir).unwrap();
}

/// Two files open on one target in one process each write their own
/// temporary file: both commits succeed and the last one wins whole.
#[test]
fn two_writers_on_one_target() {
    let dir = scratch("two");
    let target = dir.join("out.tsv");
    let mut a = AtomicFile::create(&target).unwrap();
    let mut b = AtomicFile::create(&target).unwrap();
    a.write_all(b"from a\n").unwrap();
    b.write_all(b"from b\n").unwrap();
    a.commit().unwrap();
    assert_eq!(fs::read_to_string(&target).unwrap(), "from a\n");
    b.commit().unwrap();
    assert_eq!(fs::read_to_string(&target).unwrap(), "from b\n");
    // Nothing but the target is left.
    let names: Vec<_> = fs::read_dir(&dir)
        .unwrap()
        .map(|e| e.unwrap().file_name())
        .collect();
    assert_eq!(names, ["out.tsv"]);
    fs::remove_dir_all(&dir).unwrap();
}
