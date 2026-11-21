//! og is built with `panic = "unwind"` so that `og report render` can catch
//! a panic inside one section of a report. Everywhere else a panic must
//! still end og at once, as `panic = "abort"` did: a worker thread that
//! panics and is taken for the end of its input would otherwise let og
//! write a cut file and exit 0. og's main installs a panic hook that aborts
//! outside a report section; this test checks that the binary has it.

use std::process::Command;

#[cfg(unix)]
#[test]
fn a_panic_outside_a_report_section_aborts_og() {
    use std::os::unix::process::ExitStatusExt;
    // A stdout whose reader is gone: `og version` panics when it prints
    // ("failed printing to stdout"), on the main thread, outside any
    // report section. Unwinding would exit 101; the hook aborts (SIGABRT).
    let (reader, writer) = std::io::pipe().expect("a pipe");
    drop(reader);
    let out = Command::new(env!("CARGO_BIN_EXE_og"))
        .arg("version")
        .stdout(writer)
        .output()
        .expect("og runs");
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(err.contains("failed printing to stdout"), "{err}");
    assert_eq!(
        out.status.signal(),
        Some(6),
        "og must abort, not unwind: {:?}\n{err}",
        out.status
    );
}
