//! Run as root, [`og_io::AtomicFile`] hands what it writes in another user's
//! folder, or in a subfolder root made inside it, to that user, as in a
//! container run as root over a work folder bind-mounted from the host. The output stays owner-only, and the user who
//! owns the folder can read it.
//!
//! The test needs root, so it runs in a container (`docker run --user root`).
//! Run as any other user it prints that it was skipped, unless
//! `OG_IO_ROOT_TEST=1` is set: then it fails, so a run meant to check this
//! cannot pass without checking it.

use og_io::AtomicFile;
use std::fs;
use std::io::Write;
use std::os::unix::fs::{MetadataExt, PermissionsExt, chown};
use std::path::Path;

fn write(target: &Path) {
    let mut f = AtomicFile::create(target).unwrap();
    f.write_all(b"new\n").unwrap();
    f.commit().unwrap();
}

#[test]
fn root_hands_the_file_to_the_folder_owner() {
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

    // A folder another user owns: the file goes to that user, still 0600.
    let user = dir.join("user");
    fs::create_dir(&user).unwrap();
    chown(&user, Some(4321), Some(4322)).unwrap();
    let out = user.join("out.tsv");
    write(&out);
    let m = fs::metadata(&out).unwrap();
    assert_eq!((m.uid(), m.gid()), (4321, 4322), "the folder owner's file");
    assert_eq!(m.permissions().mode() & 0o7777, 0o600);
    assert_eq!(fs::read_to_string(&out).unwrap(), "new\n");

    // A subfolder made as root inside the user's folder: still the user's.
    let sub = user.join("sub");
    fs::create_dir(&sub).unwrap();
    let deep = sub.join("deep.tsv");
    write(&deep);
    let m = fs::metadata(&deep).unwrap();
    assert_eq!(
        (m.uid(), m.gid()),
        (4321, 4322),
        "a file in a root subfolder"
    );

    // Root's own folder (the temp folder, all root's up to `/` in a
    // container): the file stays root's.
    let own = dir.join("own.tsv");
    write(&own);
    assert_eq!(fs::metadata(&own).unwrap().uid(), 0, "root's file");

    fs::remove_dir_all(&dir).unwrap();
}
