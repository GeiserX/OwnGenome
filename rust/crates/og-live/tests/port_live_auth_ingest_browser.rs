//! The login form of `og live serve` in a headless Chrome (og-browser;
//! skipped without Chrome, except in CI): a page opened without a session
//! goes to the login form, and the form leads back to that page, which
//! loads. Ported from check 7 of `tests/test_owngenome_live.py`.

#[path = "port_live_server_util.rs"]
mod util;

use og_browser::{Browser, Value, chrome_or_skip, js_str};
use og_live::auth::{self, Auth};
use std::sync::Arc;
use std::time::Duration;
use util::*;

const PASSWORD: &str = "synthetic test password 7";

fn href(b: &mut Browser) -> String {
    match b.eval("location.href + ' ' + document.documentElement.outerHTML.slice(0, 600)") {
        Ok(Value::Str(s)) => s,
        other => format!("{other:?}"),
    }
}

#[test]
fn a_page_without_a_session_goes_to_the_login_form_and_back() {
    let Some(chrome) = chrome_or_skip("load the live site's login form") else {
        return;
    };
    let w = Work::new("auth-browser");
    make_root(&w);
    let ev = fixture_events();
    store(&w, RUN, &ev.iter().map(String::as_str).collect::<Vec<_>>());
    let mut c = config(&w);
    c.auth = Arc::new(Auth::new(Some(
        auth::hash_password(PASSWORD, 1000).unwrap(),
    )));
    let srv = start(c);
    let base = format!("http://127.0.0.1:{}", srv.port());
    let mut b = Browser::launch(&chrome, &w.dir.join("chrome")).unwrap();
    b.go(&format!("{base}/runs/{RUN}/")).unwrap();
    let form = b.wait(
        "location.pathname === '/login' && !!document.getElementById('password')",
        Duration::from_secs(20),
    );
    assert!(form.is_some(), "no login form: {}", href(&mut b));
    let fill = format!(
        "document.getElementById('password').value = {}; document.querySelector('form').submit(); true",
        js_str(PASSWORD)
    );
    b.eval(&fill).unwrap();
    let back = b.wait(
        &format!(
            "location.pathname === '/runs/{RUN}/' && document.documentElement.dataset.pgpStatus === 'ready'"
        ),
        Duration::from_secs(30),
    );
    assert!(
        back.is_some(),
        "the form did not lead back to the run page: {}",
        href(&mut b)
    );
    assert!(b.page_errors().is_empty(), "{:?}", b.page_errors());
}
