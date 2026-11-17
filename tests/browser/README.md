# Browser checks

The browser checks load OwnGenome's pages in a real headless Chrome and test what the page's JavaScript does. They are Rust tests. The shared driver is the test-support crate `og-browser` in [rust/crates/og-browser](../../rust/crates/og-browser/). It needs no Python, no Node and no browser package: it starts the Chrome already on the machine and talks to it over the DevTools protocol.

## Run them

```bash
cd rust
cargo test -p og-browser                          # the driver's own tests and the ported checks
CHROME=/path/to/chrome cargo test -p og-browser   # a Chrome that is not on PATH
```

Which Chrome runs: `CHROME` when it is set, else the first of `google-chrome`, `google-chrome-stable`, `chromium` and `chromium-browser` on `PATH`, else `/Applications/Google Chrome.app`. A `CHROME` that names no executable means no Chrome: the driver does not fall back to another one.

Without a Chrome a check prints `[SKIP] no Chrome or Chromium to ...` and passes. In CI (`CI=true`) it fails instead, so a runner without Chrome cannot pass by skipping. The checks run in the Guard job `rust`, inside `cargo test --workspace`, on the runner's Chrome. Two red controls follow there: a page error planted in the clean test page (`--features planted-page-error`) must fail its test, and `CHROME=/nonexistent` must fail in CI and skip with `CI=` unset.

## What is checked so far

| Test | What it checks |
|---|---|
| `tests/cdp.rs` | the driver: a clean page loads with no page error and answers `eval`; a planted exception and a `console.error` are both reported |
| `tests/live_replay.rs` | the live page in replay mode on `tests/fixtures/live/frames.json`: it opens at the last frame, Play replays from the first, the map's roles and labels, a headline opens its step with Enter, the policy blocked nothing, a running replay opens at its first frame. Red control: without `frames.json` the first check fails |
| `tests/live_assets.rs` | `assets/live`, plugins included, has no outside URL, inline script, `on*` attribute or `eval`, with two planted mistakes the scan must catch |

The other browser checks still run in the Python tests that start the live site's Python server (`tests/test_owngenome_live.py`, `test_live_demo.py`, `test_live_viewer.py`, `test_owngenome_live_family.py`, `test_report_charts.py`, `test_viewer.py`). Each moves here when its server is ported to Rust.

## The API

Everything a check needs is in `og_browser`:

| Item | What it does |
|---|---|
| `chrome_or_skip(what) -> Option<PathBuf>` | the Chrome to run, or `None` after a `[SKIP]` line; panics in CI without one |
| `find_chrome()` | the Chrome to run, or `None` |
| `Browser::launch(chrome, profile)` | start a headless Chrome with its profile folder and attach to a blank page; tries three times, since a Chrome on a busy runner sometimes never answers its first call |
| `b.go(url)` | open a URL; it returns once the navigation starts |
| `b.eval(expr) -> Result<Value, String>` | the value of a JavaScript expression (`Value::Null` for `undefined`), or the exception it threw |
| `b.wait(expr, timeout) -> Option<Value>` | evaluate every 0.2 s until the value is truthy; `None` at the timeout |
| `b.page_errors() -> Vec<String>` | the page's uncaught exceptions and `console.error` calls so far |
| `b.call(method, params)` | any DevTools method on the page, for what the helpers above do not cover |
| `js_str(text)` | a string as a JavaScript literal, to build an expression |
| `StaticServer::start(dir, headers)` | serve a folder on a free port of 127.0.0.1 with extra headers (a Content-Security-Policy); `srv.url(path)` |
| `asset_problems(texts)`, `read_assets(dir)` | the asset scan above |
| `repo_root()` | the repository's root, to find `assets/` and `tests/fixtures/` |

Dropping a `Browser` closes Chrome; dropping a `StaticServer` stops it. Put each test's Chrome profile and site under `env!("CARGO_TARGET_TMPDIR")`.

How it talks to Chrome: Chrome starts with `--remote-debugging-pipe`, which reads commands on file descriptor 3 and writes answers on 4, one JSON message per NUL byte. `/bin/sh` moves the test's two pipes to those numbers before it starts Chrome. No port is opened and no unsafe code is needed.

## Write a new check

1. Add `rust/crates/og-browser/tests/<what>.rs`, one file per check, so two lanes never edit the same file.
2. Start with `let Some(chrome) = chrome_or_skip("...") else { return };`.
3. Serve the page with the server that will serve it for real (an `og` verb) or, for static pages, with `StaticServer`.
4. Give the check a red control: the same page with one thing broken must fail it. Run it red once and say where.
