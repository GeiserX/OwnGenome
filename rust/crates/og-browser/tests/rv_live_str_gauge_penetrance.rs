//! The repeat gauge of the live person page, in a headless Chrome, on the
//! demo the docs publish (og live demo, run with the og this test run
//! built), with its STR widget changed so a locus has no Stranger
//! thresholds.
//!
//! The range of reduced penetrance comes from og's own table, not from
//! Stranger, so `og summary collect` attaches it whether or not the locus
//! has Stranger thresholds. The gauge's legend names it in both cases. A
//! locus without thresholds and without a range keeps the old sentence.
//!
//! Red control: with the person.js of the head before the fix, the legend of
//! a locus without thresholds drops the range, so the first check fails.
//!
//! Without Chrome or a built og the test says SKIP, and fails in CI.

use og_browser::{Browser, Scratch, StaticServer, Value, chrome_or_skip, js_str, repo_root};
use pycompat::json::{dumps, loads};
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::Duration;

const RUN: &str = "DEMO-001";
const WAIT: Duration = Duration::from_secs(15);

/// The og to build the demo with: OG when set and not empty, else the og
/// beside this test's own folder (target/<profile>/og, which
/// `cargo test --workspace` builds).
fn find_og() -> Option<PathBuf> {
    if let Some(o) = std::env::var_os("OG").filter(|o| !o.is_empty()) {
        let p = PathBuf::from(o);
        return p.is_file().then_some(p);
    }
    let me = std::env::current_exe().ok()?;
    let p = me.parent()?.parent()?.join("og");
    p.is_file().then_some(p)
}

fn object(v: &mut Value) -> &mut indexmap::IndexMap<String, Value> {
    let Value::Object(o) = v else {
        panic!("expected an object, got {v:?}")
    };
    o
}

/// The `loci` of the demo's STR gauge, wherever results.json holds it.
fn str_loci(v: &mut Value) -> Option<&mut Vec<Value>> {
    match v {
        Value::Object(o) => {
            if matches!(o.get("kind"), Some(Value::Str(k)) if k == "str")
                && matches!(o.get("loci"), Some(Value::Array(_)))
            {
                let Some(Value::Array(loci)) = o.get_mut("loci") else {
                    unreachable!()
                };
                return Some(loci);
            }
            o.values_mut().find_map(str_loci)
        }
        Value::Array(a) => a.iter_mut().find_map(str_loci),
        _ => None,
    }
}

#[test]
fn the_str_gauge_names_reduced_penetrance_without_stranger_thresholds() {
    let Some(chrome) = chrome_or_skip("open the live demo's repeat gauge") else {
        return;
    };
    let Some(og) = find_og() else {
        assert_ne!(std::env::var("CI").as_deref(), Ok("true"), "build og in CI");
        println!("[SKIP] no og to build the live demo");
        return;
    };
    // Declared before the server and the browser: dropped after them, also on a panic.
    let scratch = Scratch::new(Path::new(env!("CARGO_TARGET_TMPDIR")), "str-gauge");
    let site = scratch.path().join("site");
    let built = Command::new(&og)
        .args(["live", "demo", "--out"])
        .current_dir(repo_root())
        .arg(&site)
        .arg("--work")
        .arg(scratch.path().join("work"))
        .output()
        .unwrap();
    assert!(
        built.status.success(),
        "og live demo failed:\n{}{}",
        String::from_utf8_lossy(&built.stdout),
        String::from_utf8_lossy(&built.stderr)
    );

    // The first locus loses both thresholds and gains a range of reduced
    // penetrance; the second loses both thresholds only (the control).
    let path = site.join("runs").join(RUN).join("results.json");
    let mut results = loads(&fs::read_to_string(&path).unwrap()).unwrap();
    let loci = str_loci(&mut results).expect("the demo has an STR gauge");
    assert!(loci.len() >= 3, "{loci:?}");
    for (i, l) in loci.iter_mut().take(2).enumerate() {
        let l = object(l);
        l.shift_remove("normal_max");
        l.shift_remove("pathologic_min");
        if i == 0 {
            l.insert("penetrance".into(), Value::Str("36-39".into()));
        }
    }
    fs::write(&path, dumps(&results, &Default::default()).unwrap()).unwrap();

    let srv = StaticServer::start(&site, &[]).unwrap();
    let mut b = Browser::launch(&chrome, &scratch.path().join("chrome")).unwrap();
    b.go(&srv.url(&format!("/runs/{RUN}/"))).unwrap();
    let rows = "#card-expansions [data-widget=str] .w-row .w-legend";
    let legends = b
        .wait(
            &format!(
                "(() => {{ const l = [...document.querySelectorAll({})].map(p => p.textContent); \
                 return l.length >= 3 ? JSON.stringify(l) : ''; }})()",
                js_str(rows)
            ),
            WAIT,
        )
        .unwrap_or_else(|| panic!("no STR gauge: {:?}", b.page_errors()));
    let Value::Str(legends) = legends else {
        panic!("{legends:?}")
    };
    let Value::Array(legends) = loads(&legends).unwrap() else {
        panic!("{legends}")
    };
    let legend = |i: usize| match &legends[i] {
        Value::Str(s) => s.clone(),
        v => panic!("{v:?}"),
    };

    let first = legend(0);
    assert!(
        first.starts_with("reduced penetrance 36-39 · No Stranger thresholds for this locus."),
        "a locus without thresholds dropped its range of reduced penetrance: {first:?}"
    );
    let second = legend(1);
    assert!(
        second.starts_with("No Stranger thresholds for this locus. Scale 0 to "),
        "{second:?}"
    );
    assert!(!second.contains("penetrance"), "{second:?}");
    let third = legend(2);
    assert!(third.starts_with("\"normal\" up to "), "{third:?}");
    assert!(b.page_errors().is_empty(), "{:?}", b.page_errors());
}
