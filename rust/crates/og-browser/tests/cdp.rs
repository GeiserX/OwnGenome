//! The driver itself, on two tiny pages: a clean one must load with no page
//! error, and one with a planted exception and `console.error` must report
//! both. With the `planted-page-error` feature the clean page throws too,
//! and its test must fail (the red control in guard.yml). Without Chrome
//! both say SKIP, and fail in CI.

use og_browser::{Browser, Scratch, StaticServer, Value, chrome_or_skip, js_str};
use std::fs;
use std::path::Path;
use std::time::Duration;

const CSP: &str = "default-src 'none'; script-src 'self'";

const PAGE: &str = "<!doctype html><html><head><script src=\"app.js\" defer></script></head><body></body></html>\n";

#[cfg(not(feature = "planted-page-error"))]
const CLEAN_JS: &str = "document.body.dataset.ok = 'yes'; window.answer = 6 * 7;\n";
#[cfg(feature = "planted-page-error")]
const CLEAN_JS: &str = "document.body.dataset.ok = 'yes'; window.answer = 6 * 7; null.planted;\n";

const ERRORS_JS: &str = "document.body.dataset.ok = 'yes'; console.error('planted console error', 7); \
     setTimeout(function () { throw new Error('planted page error'); }, 0);\n";

/// Serve PAGE with `js` as its script, open it and wait until it ran. The
/// scratch folder comes first, so the caller drops it after the browser
/// and the server, also when the test panics.
fn open(name: &str, js: &str) -> Option<(Scratch, StaticServer, Browser)> {
    let chrome = chrome_or_skip("load the test page")?;
    let d = Scratch::new(Path::new(env!("CARGO_TARGET_TMPDIR")), name);
    let site = d.path().join("site");
    fs::create_dir(&site).unwrap();
    fs::write(site.join("index.html"), PAGE).unwrap();
    fs::write(site.join("app.js"), js).unwrap();
    let srv = StaticServer::start(&site, &[("Content-Security-Policy", CSP)]).unwrap();
    let mut b = Browser::launch(&chrome, &d.path().join("chrome")).unwrap();
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
    Some((d, srv, b))
}

#[test]
fn a_clean_page_loads_with_no_page_error() {
    let Some((_d, _srv, mut b)) = open("clean", CLEAN_JS) else {
        return;
    };
    assert_eq!(b.eval("window.answer").unwrap(), Value::Int(42));
    // Numbers JSON has no form for come back as floats; a BigInt is an error.
    assert_eq!(b.eval("1 / 0").unwrap(), Value::Float(f64::INFINITY));
    assert_eq!(b.eval("-1 / 0").unwrap(), Value::Float(f64::NEG_INFINITY));
    assert!(matches!(b.eval("0 / 0").unwrap(), Value::Float(x) if x.is_nan()));
    assert!(matches!(b.eval("-0").unwrap(), Value::Float(x) if x == 0.0 && x.is_sign_negative()));
    assert!(b.eval("10n").unwrap_err().contains("no JSON form"));
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
    // wait never runs past its timeout, even with a slow expression.
    let t = std::time::Instant::now();
    assert_eq!(
        b.wait(
            "(function () { var e = Date.now() + 3000; while (Date.now() < e) {} return 0; })()",
            Duration::from_millis(500)
        ),
        None
    );
    assert!(t.elapsed() < Duration::from_secs(2), "{:?}", t.elapsed());
    std::thread::sleep(Duration::from_millis(300));
    assert_eq!(
        b.page_errors(),
        Vec::<String>::new(),
        "the clean page reported errors"
    );
    // A navigation that cannot start is an error; the page stays usable.
    let e = b.go("http://127.0.0.1:1/").unwrap_err();
    assert!(e.contains("net::ERR_"), "{e}");
    assert_eq!(b.eval("6 * 7").unwrap(), Value::Int(42));
}

#[test]
fn a_planted_exception_and_console_error_are_reported() {
    let Some((_d, _srv, mut b)) = open("errors", ERRORS_JS) else {
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
}
