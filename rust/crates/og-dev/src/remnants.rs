//! `og dev remnants`: a Python file kept after its tool moved to `og`,
//! only because other Python still imports names from it, may hold those
//! names and nothing else. `tests/parity/remnants.tsv` lists each remnant
//! with the names it may define, the files that import it, and the lane
//! that deletes it. The check fails when:
//!
//! - the remnant has no `# remnant:` header naming that lane;
//! - it defines a top-level name the row does not allow, or the row allows
//!   a name it no longer defines;
//! - an allowed name is used by no importer and nowhere else in the remnant;
//! - a listed importer is gone or no longer imports the module;
//! - a Python file under `bin/`, `scripts/` or `tests/` imports the module
//!   and is not listed.
//!
//! So a remnant only shrinks, and its last importer leaving shows it can go.

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::{Path, PathBuf};

use regex::Regex;

const USAGE: &str = "usage: og dev remnants [--list FILE] [--root DIR]\n";

/// One row of the list.
struct Remnant {
    file: String,
    names: BTreeSet<String>,
    importers: Vec<String>,
    lane: String,
}

fn read_list(path: &Path) -> Result<Vec<Remnant>, String> {
    let text =
        fs::read_to_string(path).map_err(|e| format!("cannot read {}: {e}", path.display()))?;
    let mut out = Vec::new();
    for (i, line) in text.lines().enumerate() {
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let c: Vec<&str> = line.split('\t').collect();
        let [file, names, importers, lane] = c.as_slice() else {
            return Err(format!(
                "{} line {}: {} fields, want 4: file, names, importers, lane",
                path.display(),
                i + 1,
                c.len()
            ));
        };
        let list = |s: &str| -> Vec<String> {
            s.split(',')
                .map(str::trim)
                .filter(|v| !v.is_empty())
                .map(str::to_owned)
                .collect()
        };
        out.push(Remnant {
            file: (*file).to_owned(),
            names: list(names).into_iter().collect(),
            importers: list(importers),
            lane: (*lane).to_owned(),
        });
    }
    Ok(out)
}

/// The names a Python file defines at top level: `def`, `class`, and
/// assignments, at column 0. Dunder names are skipped.
fn top_level(text: &str) -> BTreeMap<String, usize> {
    let def = Regex::new(r"^(?:async\s+def|def|class)\s+([A-Za-z_]\w*)").expect("regex");
    let assign = Regex::new(r"^([A-Za-z_]\w*)\s*(?::[^=]*)?=[^=]").expect("regex");
    let mut out = BTreeMap::new();
    for (i, line) in text.lines().enumerate() {
        let name = def
            .captures(line)
            .or_else(|| assign.captures(line))
            .map(|c| c[1].to_owned());
        if let Some(n) = name.filter(|n| !n.starts_with("__")) {
            out.entry(n).or_insert(i);
        }
    }
    out
}

/// How a file imports `module`: the names it uses from it, or `None` when
/// it does not import it. `import m` and `import m as x` count `x.name`;
/// `from m import a, b as c` counts `a` and `b`.
fn uses(text: &str, module: &str) -> Option<BTreeSet<String>> {
    let m = regex::escape(module);
    let plain = Regex::new(&format!(
        r"(?m)^[ \t]*import[ \t]+{m}(?:[ \t]+as[ \t]+([A-Za-z_]\w*))?[ \t]*(?:#.*)?$"
    ))
    .expect("regex");
    let from = Regex::new(&format!(
        r"(?ms)^[ \t]*from[ \t]+{m}[ \t]+import[ \t]+(\([^)]*\)|[^\n]*)"
    ))
    .expect("regex");
    let mut found = false;
    let mut used = BTreeSet::new();
    for c in plain.captures_iter(text) {
        found = true;
        let alias = c.get(1).map_or(module, |a| a.as_str());
        let attr =
            Regex::new(&format!(r"\b{}\.([A-Za-z_]\w*)", regex::escape(alias))).expect("regex");
        used.extend(attr.captures_iter(text).map(|a| a[1].to_owned()));
    }
    for c in from.captures_iter(text) {
        found = true;
        let list = c[1]
            .split('#')
            .next()
            .unwrap_or_default()
            .trim_matches(|ch| ch == '(' || ch == ')');
        for item in list.split(',') {
            if let Some(name) = item.split_whitespace().next() {
                used.insert(name.to_owned());
            }
        }
    }
    found.then_some(used)
}

/// The Python files under `dir`, skipping hidden folders.
fn python_files(dir: &Path, out: &mut Vec<PathBuf>) {
    let Ok(rd) = fs::read_dir(dir) else { return };
    for e in rd.flatten() {
        let p = e.path();
        let hidden = e.file_name().to_string_lossy().starts_with('.');
        if p.is_dir() && !hidden {
            python_files(&p, out);
        } else if p.extension().is_some_and(|x| x == "py") {
            out.push(p);
        }
    }
}

/// The problems of one remnant, relative to `root`.
fn check(root: &Path, r: &Remnant) -> Result<Vec<String>, String> {
    let path = root.join(&r.file);
    let text = fs::read_to_string(&path).map_err(|e| format!("cannot read {}: {e}", r.file))?;
    let module = Path::new(&r.file)
        .file_stem()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_default();
    let mut bad = Vec::new();
    let head = text.lines().take(5).find(|l| l.starts_with("# remnant:"));
    match head {
        None => bad.push(format!(
            "{}: no '# remnant:' header in its first 5 lines",
            r.file
        )),
        Some(h) if !h.contains(&r.lane) => {
            bad.push(format!(
                "{}: the remnant header does not name lane {}",
                r.file, r.lane
            ));
        }
        Some(_) => {}
    }
    let defined = top_level(&text);
    for n in defined.keys().filter(|n| !r.names.contains(*n)) {
        bad.push(format!(
            "{}: defines '{n}', which no importer is allowed to need",
            r.file
        ));
    }
    for n in r.names.iter().filter(|n| !defined.contains_key(*n)) {
        bad.push(format!(
            "{}: '{n}' is allowed but no longer defined; drop it from the list",
            r.file
        ));
    }

    let mut used: BTreeSet<String> = BTreeSet::new();
    for imp in &r.importers {
        let Ok(t) = fs::read_to_string(root.join(imp)) else {
            bad.push(format!(
                "{}: importer {imp} is gone; drop it from the list",
                r.file
            ));
            continue;
        };
        match uses(&t, &module) {
            Some(u) => used.extend(u),
            None => bad.push(format!(
                "{}: {imp} no longer imports {module}; drop it from the list",
                r.file
            )),
        }
    }
    // A name the remnant itself uses outside its own definition line.
    let lines: Vec<&str> = text.lines().collect();
    for (n, &at) in &defined {
        let word = Regex::new(&format!(r"\b{}\b", regex::escape(n))).expect("regex");
        if lines
            .iter()
            .enumerate()
            .any(|(i, l)| i != at && word.is_match(l))
        {
            used.insert(n.clone());
        }
    }
    for n in r
        .names
        .iter()
        .filter(|n| defined.contains_key(*n) && !used.contains(*n))
    {
        bad.push(format!(
            "{}: '{n}' is used by no importer; delete it",
            r.file
        ));
    }

    let mut files = Vec::new();
    for d in ["bin", "scripts", "tests"] {
        python_files(&root.join(d), &mut files);
    }
    files.sort();
    for f in files {
        let rel = f
            .strip_prefix(root)
            .unwrap_or(&f)
            .to_string_lossy()
            .into_owned();
        if rel == r.file || r.importers.contains(&rel) {
            continue;
        }
        if fs::read_to_string(&f).is_ok_and(|t| uses(&t, &module).is_some()) {
            bad.push(format!(
                "{}: {rel} imports {module} but is not listed; port it, or list it",
                r.file
            ));
        }
    }
    Ok(bad)
}

/// Runs `og dev remnants ARGS`: 0 when every remnant holds only what its
/// importers use, 1 when one does not, 2 on a usage or read error.
pub fn cli(args: &[String]) -> u8 {
    let (mut list, mut root) = (
        PathBuf::from("tests/parity/remnants.tsv"),
        PathBuf::from("."),
    );
    let mut it = args.iter();
    while let Some(a) = it.next() {
        match (a.as_str(), it.clone().next()) {
            ("--list", Some(v)) => {
                list = PathBuf::from(v);
                it.next();
            }
            ("--root", Some(v)) => {
                root = PathBuf::from(v);
                it.next();
            }
            ("-h" | "--help", _) => {
                print!("{USAGE}");
                return 0;
            }
            _ => {
                eprint!("{USAGE}");
                eprintln!("og dev remnants: error: unexpected argument '{a}'");
                return 2;
            }
        }
    }
    let list = if list.is_absolute() {
        list
    } else {
        root.join(list)
    };
    let rows = match read_list(&list) {
        Ok(r) => r,
        Err(e) => {
            eprintln!("og dev remnants: {e}");
            return 2;
        }
    };
    let mut failed = false;
    for r in &rows {
        match check(&root, r) {
            Ok(bad) if bad.is_empty() => println!(
                "ok      {}: {} names, {} importers, deleted in lane {}",
                r.file,
                r.names.len(),
                r.importers.len(),
                r.lane
            ),
            Ok(bad) => {
                failed = true;
                for b in bad {
                    println!("FAIL    {b}");
                }
            }
            Err(e) => {
                eprintln!("og dev remnants: {e}");
                return 2;
            }
        }
    }
    u8::from(failed)
}

#[cfg(test)]
mod tests {
    use super::*;

    struct Tree(PathBuf);

    impl Tree {
        fn new(tag: &str) -> Tree {
            let d = std::env::temp_dir().join(format!("og-remnants-{tag}-{}", std::process::id()));
            let _ = fs::remove_dir_all(&d);
            fs::create_dir_all(d.join("bin")).unwrap();
            fs::create_dir_all(d.join("tests")).unwrap();
            Tree(d)
        }
        fn put(&self, rel: &str, text: &str) {
            fs::write(self.0.join(rel), text).unwrap();
        }
        fn run(&self) -> (u8, Vec<String>) {
            let r = &read_list(&self.0.join("list.tsv")).unwrap()[0];
            let bad = check(&self.0, r).unwrap();
            (u8::from(!bad.is_empty()), bad)
        }
    }

    impl Drop for Tree {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.0);
        }
    }

    const LIB: &str = "#!/usr/bin/env python3\n# remnant: kept for a, b; deleted in lane R9\nimport os\n\
                       MISSING = {'.', ''}\nSEX: dict = {}\n\n\ndef missing(v):\n    return v in MISSING\n\n\
                       def die(m):\n    raise SystemExit(m)\n\n\nclass Thing:\n    pass\n";

    fn good(tag: &str) -> Tree {
        let t = Tree::new(tag);
        t.put("bin/lib.py", LIB);
        t.put(
            "bin/a.py",
            "import lib  # the rules\nlib.die(lib.missing(1))\n",
        );
        t.put(
            "bin/b.py",
            "def f():\n    from lib import (SEX,\n        Thing as T)\n",
        );
        t.put("list.tsv", "# file\tnames\timporters\tlane\nbin/lib.py\tMISSING,SEX,missing,die,Thing\tbin/a.py,bin/b.py\tR9\n");
        t
    }

    #[test]
    fn a_clean_remnant_passes() {
        let t = good("clean");
        assert_eq!(t.run(), (0, vec![]));
    }

    #[test]
    fn an_extra_function_fails() {
        let t = good("extra");
        t.put("bin/lib.py", &format!("{LIB}\n\ndef helper():\n    pass\n"));
        let (code, bad) = t.run();
        assert_eq!(code, 1);
        assert_eq!(
            bad,
            ["bin/lib.py: defines 'helper', which no importer is allowed to need"]
        );
    }

    #[test]
    fn a_dead_name_fails() {
        let t = good("dead");
        t.put("bin/a.py", "import lib as L\nL.die(1)\n");
        let (code, bad) = t.run();
        assert_eq!(code, 1);
        // MISSING stays used inside the remnant; missing has no user left.
        assert_eq!(
            bad,
            ["bin/lib.py: 'missing' is used by no importer; delete it"]
        );
    }

    #[test]
    fn a_removed_importer_fails() {
        let t = good("gone");
        t.put("bin/b.py", "print('ported')\n");
        let (code, bad) = t.run();
        assert_eq!(code, 1);
        assert!(
            bad.contains(
                &"bin/lib.py: bin/b.py no longer imports lib; drop it from the list".to_owned()
            ),
            "{bad:?}"
        );
        fs::remove_file(t.0.join("bin/b.py")).unwrap();
        let (_, bad) = t.run();
        assert!(
            bad.contains(
                &"bin/lib.py: importer bin/b.py is gone; drop it from the list".to_owned()
            ),
            "{bad:?}"
        );
    }

    #[test]
    fn a_new_importer_or_a_lost_header_fails() {
        let t = good("new");
        t.put("tests/test_x.py", "import lib\n");
        let (_, bad) = t.run();
        assert_eq!(
            bad,
            ["bin/lib.py: tests/test_x.py imports lib but is not listed; port it, or list it"]
        );
        let t = good("head");
        t.put(
            "bin/lib.py",
            &LIB.replace("# remnant: kept for a, b; deleted in lane R9\n", ""),
        );
        let (_, bad) = t.run();
        assert_eq!(
            bad,
            ["bin/lib.py: no '# remnant:' header in its first 5 lines"]
        );
    }

    #[test]
    fn names_and_imports_are_read_like_python() {
        let d = top_level(
            "A = 1\nB: int = 2\nC == 3\n  D = 4\ndef e():\n    F = 5\nclass G(object):\n__all__ = []\nasync def h():\n",
        );
        assert_eq!(
            d.keys().cloned().collect::<Vec<_>>(),
            ["A", "B", "G", "e", "h"]
        );
        assert_eq!(uses("import libx\n", "lib"), None);
        assert_eq!(uses("import lib\n", "lib"), Some(BTreeSet::new()));
        let u = uses("from lib import a, b as c  # x\n", "lib").unwrap();
        assert_eq!(u.into_iter().collect::<Vec<_>>(), ["a", "b"]);
    }
}
