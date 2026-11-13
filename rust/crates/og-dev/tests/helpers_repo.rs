//! `og dev helpers` on this checkout: rule 1 holds, and the table still
//! sees the calls of both the scripts and the modules, so a reader that
//! stopped understanding the tree cannot pass by finding nothing.

use std::path::Path;

fn run(args: &[&str]) -> (u8, String) {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../..");
    let mut a: Vec<String> = args.iter().map(|s| (*s).to_owned()).collect();
    a.push("--root".into());
    a.push(root.to_str().expect("a UTF-8 path").to_owned());
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = og_dev::helpers::cli(&a, &mut out, &mut err);
    let mut s = String::from_utf8(out).expect("UTF-8");
    s.push_str(&String::from_utf8(err).expect("UTF-8"));
    (code, s)
}

#[test]
fn rule_one_holds_on_this_checkout() {
    let (code, out) = run(&["--static"]);
    assert_eq!(code, 0, "{out}");
    assert!(
        out.starts_with("OK: no bgzip/tabix call in BCFTOOLS_IMAGE or SAMTOOLS_IMAGE."),
        "{out}"
    );
}

#[test]
fn the_table_sees_scripts_and_modules() {
    let (code, out) = run(&["--list"]);
    assert_eq!(code, 0, "{out}");
    let rows: Vec<&str> = out.lines().collect();
    let scripts = rows.iter().filter(|r| r.contains(" scripts/")).count();
    let modules = rows
        .iter()
        .filter(|r| r.contains(" modules/local/"))
        .count();
    assert!(scripts >= 50, "{scripts} script rows:\n{out}");
    assert!(modules >= 20, "{modules} module rows:\n{out}");
}
