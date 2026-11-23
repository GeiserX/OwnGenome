# Browser checks

The browser checks load OwnGenome's pages in a real headless Chrome and test what the page's JavaScript does. They are Rust tests. The shared driver is the test-support crate `og-browser` in [rust/crates/og-browser](../../rust/crates/og-browser/). It needs no Python, no Node and no browser package: it starts the Chrome already on the machine and talks to it over the DevTools protocol.

## Run them

```bash
cd rust
cargo test -p og-browser                          # the driver's own tests and the ported checks
CHROME=/path/to/chrome cargo test -p og-browser   # a Chrome that is not on PATH
```

Which Chrome runs: `CHROME` when it is set and not empty, else the first of `google-chrome`, `google-chrome-stable`, `chromium` and `chromium-browser` on `PATH`, else `/Applications/Google Chrome.app`. An empty `CHROME` counts as unset. A `CHROME` that is not empty but names no executable means no Chrome: the driver does not fall back to another one.

Without a Chrome a check prints `[SKIP] no Chrome or Chromium to ...` and passes. The test harness hides the output of a passing test, so the line shows only with `cargo test -p og-browser -- --nocapture`; without it a skipped check reads as a plain `ok`. In CI (`CI=true`) it fails instead, so a runner without Chrome cannot pass by skipping. The checks run in the Guard job `rust`, inside `cargo test --workspace`, on the runner's Chrome. Two red controls follow there: a page error planted in the clean test page (`--features planted-page-error`) must fail its test, and `CHROME=/nonexistent` must fail in CI and skip with `CI=` unset.

## What is checked so far

| Test | What it checks |
|---|---|
| `tests/cdp.rs` | the driver: a clean page loads with no page error and answers `eval`; a planted exception and a `console.error` are both reported |
| `tests/live_replay.rs` | the live page in replay mode on `tests/fixtures/live/frames.json`: it opens at the last frame, Play replays from the first, the map's roles and labels, a headline opens its step with Enter, the policy blocked nothing, a running replay opens at its first frame. Red control: without `frames.json` the page gets a 404 for it and the first check fails. The policy it is served with must equal the live site's `CSP` |
| `tests/live_assets.rs` | `assets/live`, plugins included, has no outside URL, inline script, `on*` attribute or `eval`, with two planted mistakes the scan must catch |
| `tests/rv_live_reads_sketch_label.rs` | the reads offers of the live demo, built by `og live demo` with the og of the test run: a finding's reads drawing says it is an illustration before the picture, and no reads offer in the demo (the finding page, a card's links, a chart mark's panel, a search hit with a place, the map's card) links to the viewer page; each says "Reads: not in the demo". Outside the demo the finding page keeps its link. After a search for CYP2D6 and the card's reads button, the top of the reads panel is in view, on a wide screen and on a phone. Red control: with a map.js that does not scroll to the panel, the panel stays out of view |
| `tests/rv_live_load_errors.rs` | a live page that cannot load its data says so, on a small demo site served by a stub that can fail one path, with the live site's policy: a person page whose `results.json` answers 500, a home page whose `runs.json` fails and a family page whose `family.json` answers `{}` each say what could not be loaded and offer Try again, which draws the page once the stub recovers; a closed page whose state cannot be read says "Results not loaded" and never "child", while DEMO-003, a real child, still reads "Results about a child"; while `frames.json` is held, the wire, run health and re-check panels say they wait for the run's data, and a state with no re-check then reads "No re-check yet". Each case was red on main before the fix |
| `tests/rv_viewer_shell.rs` | the genome viewer's page (`assets/viewer`) on a synthetic sample, served by a stub that answers the paths and the policy of `og view serve`: the page is dark at 1440x900 and 390x844, the igv.js search box and track label are readable, the variant track fits the live map's short frame, on a phone the page is no wider than the screen, it has one h1 naming the sample and one main landmark, and a 500 from `config.json` shows an alert with Try again, which starts the viewer once the server answers. `the_policy_is_og_view_serves` keeps the stub's policy equal to `og view serve`'s |
| `rust/crates/og-live/tests/port_viewer.rs` (`cargo test -p og-live`) | the local genome viewer, served by `og view serve` itself on a synthetic sample: igv.js starts on the config (status ready), the policy blocks nothing, and igv.js reads the reference, the BAM, its index and the VCF. Red control: with an empty igv.min.js the page reports an error |
| `rust/crates/og-live/tests/port_live_server_browser.rs` (`cargo test -p og-live`) | the live run page, served by `og live serve` on the synthetic DEMO-001 run: the step map, a box that repaints after a stored message, a stream refused with 503 and opened again, a box's detail with its tasks and its report anchors, and no policy violation or page error |
| `rust/crates/og-live/tests/port_live_auth_ingest_browser.rs` (`cargo test -p og-live`) | the live site's login form, served by `og live serve` with a password: a run page opened without a session goes to the form, and the form leads back to that page, which loads |
| `tests/browser/live_demo.rs` | the public site at a subpath, replay, search, report links, family, phone layout, and map views |
| `rust/crates/og-live/tests/demo.rs` | the live charts after sign-in, the ClinVar tick and its viewer link, and the missing-plugin control |

The other browser checks still run in the Python tests that start the live site's Python server (`test_live_viewer.py`, `test_owngenome_live_family.py`). They drive Chrome with the helpers of `tests/browser/cdp.py`. Each moves here when its server is ported to Rust.

## The API

Everything a check needs is in `og_browser`:

| Item | What it does |
|---|---|
| `chrome_or_skip(what) -> Option<PathBuf>` | the Chrome to run, or `None` after a `[SKIP]` line; panics in CI without one |
| `find_chrome()` | the Chrome to run, or `None` |
| `Browser::launch(chrome, profile)` | start a headless Chrome with its profile folder and attach to a blank page; tries three times, since a Chrome on a busy runner sometimes never answers its first call |
| `b.go(url)` | open a URL; it returns once the navigation starts. A navigation that cannot start (no server) is an error |
| `b.eval(expr) -> Result<Value, String>` | the value of a JavaScript expression (`Value::Null` for `undefined`, `Value::Float` for `Infinity`, `NaN` and `-0`), or the exception it threw; a `BigInt` is an error |
| `b.wait(expr, timeout) -> Option<Value>` | evaluate every 0.2 s until the value is truthy by Python's rule (`null`, `false`, `0`, `""`, `[]` and `{}` mean not yet), as in the Python checks; `None` at the timeout, never later |
| `b.page_errors() -> Vec<String>` | the page's uncaught exceptions and `console.error` calls so far |
| `b.call(method, params)` | any DevTools method on the page, for what the helpers above do not cover |
| `js_str(text)` | a string as a JavaScript literal, to build an expression |
| `StaticServer::start(dir, headers)` | serve a folder on a free port of 127.0.0.1 with extra headers (a Content-Security-Policy); `srv.url(path)` |
| `asset_problems(texts)`, `read_assets(dir)` | the asset scan above; `read_assets` refuses a link to a folder |
| `Scratch::new(parent, name)` | a fresh folder for one test, removed when dropped, also on a panic; `scratch.path()` |
| `repo_root()` | the repository's root, to find `assets/` and `tests/fixtures/` |

Dropping a `Browser` closes Chrome; dropping a `StaticServer` closes its connections and stops it. Put each test's Chrome profile and site in a `Scratch` under `env!("CARGO_TARGET_TMPDIR")`, declared before the server and the browser so it goes after them. `StaticServer` serves only files inside its folder, after links are followed.

How it talks to Chrome: Chrome starts with `--remote-debugging-pipe`, which reads commands on file descriptor 3 and writes answers on 4, one JSON message per NUL byte. `/bin/sh` moves the test's two pipes to those numbers before it starts Chrome. No port is opened and no unsafe code is needed.

## Write a new check

1. Add `rust/crates/og-browser/tests/<what>.rs`, one file per check, so two lanes never edit the same file.
2. Start with `let Some(chrome) = chrome_or_skip("...") else { return };`.
3. Serve the page with the server that will serve it for real (an `og` verb) or, for static pages, with `StaticServer`.
4. Give the check a red control: the same page with one thing broken must fail it. Run it red once and say where.
