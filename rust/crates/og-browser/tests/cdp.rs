//! The driver itself, on two tiny pages: a clean one must load with no page
//! error, and one with a planted exception and `console.error` must report
//! both. With the `planted-page-error` feature the clean page throws too,
//! and its test must fail (the red control in guard.yml). Without Chrome
//! both say SKIP, and fail in CI.

use og_browser::{Browser, StaticServer, Value, chrome_or_skip, js_str};
use std::fs;
use std::path::PathBuf;
use std::time::Duration;

const CSP: &str = "default-src 'none'; script-src 'self'";

const PAGE: &str = "<!doctype html><html><head><script src=\"app.js\" defer></script></head><body></body></html>\n";

#[cfg(not(feature = "planted-page-error"))]
const CLEAN_JS: &str = "document.body.dataset.ok = 'yes'; window.answer = 6 * 7;\n";
#[cfg(feature = "planted-page-error")]
const CLEAN_JS: &str = "document.body.dataset.ok = 'yes'; window.answer = 6 * 7; null.planted;\n";

const ERRORS_JS: &str = "document.body.dataset.ok = 'yes'; console.error('planted console error', 7); \
     setTimeout(function () { throw new Error('planted page error'); }, 0);\n";

/// A fresh folder for one test, under cargo's scratch folder.
fn scratch(name: &str) -> PathBuf {
    let d = PathBuf::from(env!("CARGO_TARGET_TMPDIR"))
        .join(format!("og-browser-{name}-{}", std::process::id()));
    let _ = fs::remove_dir_all(&d);
    fs::create_dir_all(d.join("site")).unwrap();
    d
}

/// Serve PAGE with `js` as its script, open it and wait until it ran.
fn open(name: &str, js: &str) -> Option<(Browser, StaticServer, PathBuf)> {
    let chrome = chrome_or_skip("load the test page")?;
    let d = scratch(name);
    fs::write(d.join("site/index.html"), PAGE).unwrap();
    fs::write(d.join("site/app.js"), js).unwrap();
    let srv = StaticServer::start(&d.join("site"), &[("Content-Security-Policy", CSP)]).unwrap();
    let mut b = Browser::launch(&chrome, &d.join("chrome")).unwrap();
    b.go(&srv.url("/")).unwrap();
    let ok = b.wait(
        "document.body && document.body.dataset.ok",
        Duration::from_secs(20),
    );
    assert_eq!(
        ok,
        Some(Value::Str("yes".into())),
        "the page did not run its script"
    );
    Some((b, srv, d))
}

#[test]
fn a_clean_page_loads_with_no_page_error() {
    let Some((mut b, srv, d)) = open("clean", CLEAN_JS) else {
        return;
    };
    assert_eq!(b.eval("window.answer").unwrap(), Value::Int(42));
    assert_eq!(
        b.eval(&format!("{} + '!'", js_str("it's \"quoted\"")))
            .unwrap(),
        Value::Str("it's \"quoted\"!".into())
    );
    assert_eq!(b.eval("undefined").unwrap(), Value::Null);
    // An exception in our own expression is an Err, not a page error.
    assert!(
        b.eval("nothing.here")
            .unwrap_err()
            .contains("nothing is not defined")
    );
    // wait gives up on a value that never turns truthy.
    assert_eq!(
        b.wait("document.body.dataset.never", Duration::from_millis(500)),
        None
    );
    std::thread::sleep(Duration::from_millis(300));
    assert_eq!(
        b.page_errors(),
        Vec::<String>::new(),
        "the clean page reported errors"
    );
    drop((b, srv));
    let _ = fs::remove_dir_all(d);
}

#[test]
fn a_planted_exception_and_console_error_are_reported() {
    let Some((mut b, srv, d)) = open("errors", ERRORS_JS) else {
        return;
    };
    let mut errors = Vec::new();
    for _ in 0..50 {
        errors = b.page_errors();
        if errors.len() >= 2 {
            break;
        }
        std::thread::sleep(Duration::from_millis(100));
    }
    assert!(
        errors
            .iter()
            .any(|e| e == "console.error: planted console error 7"),
        "{errors:?}"
    );
    assert!(
        errors
            .iter()
            .any(|e| e.contains("Error: planted page error")),
        "{errors:?}"
    );
    drop((b, srv));
    let _ = fs::remove_dir_all(d);
}
