//! [`og_io::AtomicFile`] creates its temporary file and never opens one that
//! is already there. This is the only test in its binary, so the first
//! temporary name this process tries is known: `.<target>.tmp-<pid>-0`.

use og_io::AtomicFile;
use std::fs;
use std::io::Write;
use std::os::unix::fs::symlink;

#[test]
fn a_link_at_the_temporary_name_is_left_alone() {
    let pid = std::process::id();
    let dir = std::env::temp_dir().join(format!("og-io-planted-{pid}"));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    let target = dir.join("out.tsv");
    let other = dir.join("other.tsv");
    fs::write(&other, "keep\n").unwrap();
    let missing = dir.join("missing.tsv");
    // Links at the first two names this process tries, and at the name
    // releases before 1.0.0 used (`.<target>.tmp-<pid>`).
    let planted = [
        (format!(".out.tsv.tmp-{pid}"), &other),
        (format!(".out.tsv.tmp-{pid}-0"), &other),
        (format!(".out.tsv.tmp-{pid}-1"), &missing),
    ];
    for (name, to) in &planted {
        symlink(to, dir.join(name)).unwrap();
    }

    let mut f = AtomicFile::create(&target).unwrap();
    f.write_all(b"new\n").unwrap();
    f.commit().unwrap();

    assert_eq!(fs::read_to_string(&other).unwrap(), "keep\n");
    assert!(
        !missing.exists(),
        "a file was created through a dangling link"
    );
    assert!(
        !fs::symlink_metadata(&target)
            .unwrap()
            .file_type()
            .is_symlink()
    );
    assert_eq!(fs::read_to_string(&target).unwrap(), "new\n");
    for (name, to) in &planted {
        assert_eq!(
            &fs::read_link(dir.join(name)).unwrap(),
            *to,
            "{name} changed"
        );
    }
    fs::remove_dir_all(&dir).unwrap();
}
