//! `og y tree-place` and `og y tree-positions` on a synthetic tree,
//! synthetic position tables and synthetic pileups (no real sample). It
//! replaces `tests/test_y_tree_place.py`.
//!
//! The tree:
//!
//! ```text
//! (root) -> A (a1, a2/a2alt) -> B (b1, b2, b3) -> C (c1, c2)
//!                            -> X (x1)        -> Y (y1) -> Z (z1) -> W (w1)
//!                            -> P (p1)        -> Q (q1, q1b; one site)
//! ```
//!
//! 1. every SNP derived along A-B-C: C is the terminal branch at level 3,
//!    with the counts of each level, and nothing ancestral in its block;
//! 2. one ancestral SNP inside B (b3) and no read at C: B is the terminal
//!    branch, and b3 is listed as ancestral inside the terminal block;
//! 3. no-calls: C's SNPs read at depth 9 (below 10), or by a 50/50 mix of
//!    bases, are no-calls; a SNP no table places counts as a SNP but not as
//!    on the tables;
//! 4. a negative branch blocks the branches below it unless they outvote
//!    it (x1 ancestral with one derived SNP below blocks, with two it is
//!    outvoted and reported); a branch with two ancestral SNPs blocks
//!    however much derived support sits below it; two single-ancestral
//!    branches on one path block too; two SNP names at one position are
//!    one site;
//! 5. nothing derived anywhere: Hg NA, level 0;
//! 6. the first table wins for a name both tables place, and a SNP in the
//!    pseudoautosomal region is never asked for;
//! 7. the pileup reader skips read starts, ends and indels, and reads `.`
//!    and `,` as the reference base;
//! 8. a missing table stops the tool and names `setup.sh --yleaf-data`.

use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};

const TREE: &str = r#"{"id": "", "snps": "", "children": [
 {"id": "A", "snps": "a1, a2/a2alt, par", "children": [
  {"id": "B", "snps": "b1, b2, b3, bnot", "children": [
   {"id": "C", "snps": "c1, c2", "children": []}]},
  {"id": "X", "snps": "x1", "children": [
   {"id": "Y", "snps": "y1", "children": [
    {"id": "Z", "snps": "z1", "children": [
     {"id": "W", "snps": "w1", "children": []}]}]}]},
  {"id": "P", "snps": "p1", "children": [
   {"id": "Q", "snps": "q1, q1b", "children": []}]}]}]}"#;

/// name, position, ancestral, derived; bnot is on no table.
const FTDNA: [(&str, u32, char, char); 14] = [
    ("a1", 3000100, 'C', 'T'),
    ("a2alt", 3000200, 'G', 'A'),
    ("b1", 3000300, 'A', 'G'),
    ("b2", 3000400, 'T', 'C'),
    ("b3", 3000500, 'C', 'G'),
    ("c1", 3000600, 'G', 'T'),
    ("x1", 3000800, 'A', 'C'),
    ("y1", 3000900, 'T', 'A'),
    ("z1", 3001000, 'G', 'C'),
    ("w1", 3001100, 'C', 'A'),
    ("p1", 3001200, 'T', 'G'),
    ("q1", 3001300, 'A', 'T'),
    ("q1b", 3001300, 'A', 'T'),
    ("par", 2000000, 'A', 'G'),
];
const YFULL: [(&str, u32, char, char); 2] = [("a1", 3999999, 'C', 'T'), ("c2", 3000700, 'A', 'T')];

fn scratch(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("og-y-tree-{}-{name}", std::process::id()));
    let _ = fs::remove_dir_all(&d);
    fs::create_dir_all(&d).unwrap();
    d
}

fn write_tables(d: &Path) {
    let table = |rows: &[(&str, u32, char, char)]| -> String {
        rows.iter()
            .map(|(n, p, a, de)| format!("chry\t{n}\tR-{n}\t{p}\t{a}->{de}\t{a}\t{de}\n"))
            .collect()
    };
    fs::write(d.join("current_tree.json"), TREE).unwrap();
    fs::write(d.join("positions_ftdna_hg38.txt"), table(&FTDNA)).unwrap();
    fs::write(d.join("positions_yfull_hg38.txt"), table(&YFULL)).unwrap();
    fs::write(
        d.join("version.txt"),
        "tree\tYFull YTree test\npositions\ttest tables\n",
    )
    .unwrap();
}

fn pos_of(name: &str) -> (u32, char, char) {
    FTDNA
        .iter()
        .chain(YFULL.iter())
        .find(|r| r.0 == name)
        .map(|r| (r.1, r.2, r.3))
        .unwrap()
}

/// How a SNP reads in the pileup: derived, ancestral, a 50/50 mix, or
/// derived at another depth.
#[derive(Clone, Copy)]
enum St {
    D,
    A,
    Mix,
    DAt(usize),
}

fn pileup(path: &Path, states: &[(&str, St)]) {
    let mut text = String::new();
    for &(name, st) in states {
        let (pos, anc, der) = pos_of(name);
        let dp = if let St::DAt(n) = st { n } else { 12 };
        let bases: String = match st {
            St::Mix => format!(
                "{}{}",
                anc.to_string().repeat(dp / 2),
                der.to_string().repeat(dp - dp / 2)
            ),
            St::A => anc.to_string().repeat(dp),
            St::D | St::DAt(_) => der.to_string().repeat(dp),
        };
        text += &format!("chrY\t{pos}\tN\t{dp}\t{bases}\t{}\n", "I".repeat(dp));
    }
    fs::write(path, text).unwrap();
}

type Row = HashMap<String, String>;

struct Placed {
    code: u8,
    stdout: String,
    stderr: String,
    row: Row,
    path: Vec<Row>,
}

impl Placed {
    fn get(&self, k: &str) -> &str {
        self.row.get(k).map_or("", String::as_str)
    }
    fn branches(&self) -> Vec<&str> {
        self.path.iter().map(|r| r["Branch"].as_str()).collect()
    }
}

fn y(verb: &str, args: &[String]) -> (u8, String, String) {
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = og_reports::y::cli(Some(verb), args, &mut out, &mut err);
    (
        code,
        String::from_utf8(out).unwrap(),
        String::from_utf8(err).unwrap(),
    )
}

fn place(d: &Path, pileup_file: &Path, tag: &str, extra: &[&str]) -> Placed {
    let prefix = d.join(tag);
    let mut args: Vec<String> = [
        "--tree-dir",
        d.to_str().unwrap(),
        "--pileup",
        pileup_file.to_str().unwrap(),
        "--sample",
        "S",
        "--out-prefix",
        prefix.to_str().unwrap(),
    ]
    .map(str::to_owned)
    .to_vec();
    args.extend(extra.iter().map(|s| (*s).to_owned()));
    let (code, stdout, stderr) = y("tree-place", &args);
    let mut p = Placed {
        code,
        stdout,
        stderr,
        row: Row::new(),
        path: Vec::new(),
    };
    if code != 0 {
        return p;
    }
    let table = |suffix: &str| -> Vec<Row> {
        let text = fs::read_to_string(format!("{}{suffix}", prefix.display())).unwrap();
        let mut lines = text.lines();
        let head: Vec<&str> = lines.next().unwrap().split('\t').collect();
        lines
            .map(|l| {
                head.iter()
                    .zip(l.split('\t'))
                    .map(|(k, v)| ((*k).to_owned(), v.to_owned()))
                    .collect()
            })
            .collect()
    };
    p.row = table("_y_tree.txt").into_iter().next().unwrap_or_default();
    p.path = table("_y_tree_path.tsv");
    p
}

fn run(d: &Path, states: &[(&str, St)], tag: &str, extra: &[&str]) -> Placed {
    let p = d.join(format!("{tag}.pileup"));
    pileup(&p, states);
    place(d, &p, tag, extra)
}

fn with(base: &[(&'static str, St)], more: &[(&'static str, St)]) -> Vec<(&'static str, St)> {
    base.iter().chain(more).copied().collect()
}

const ALL_A: [(&str, St); 2] = [("a1", St::D), ("a2alt", St::D)];
const ALL_B: [(&str, St); 3] = [("b1", St::D), ("b2", St::D), ("b3", St::D)];

#[test]
fn tree_place_and_tree_positions() {
    let d = scratch("main");
    write_tables(&d);
    let ab: Vec<(&str, St)> = with(&ALL_A, &ALL_B);

    // 1. deepest branch
    let r = run(
        &d,
        &with(&ab, &[("c1", St::D), ("c2", St::D), ("x1", St::A)]),
        "deep",
        &[],
    );
    assert_eq!(r.code, 0, "1. the tool runs: {}", r.stderr);
    assert_eq!((r.get("Hg"), r.get("Level")), ("C", "3"), "1. C at level 3");
    assert_eq!(r.branches(), ["A", "B", "C"], "1. the path is A, B, C");
    let b = &r.path[1];
    assert_eq!(
        ["SNPs", "On_tables", "Derived", "Ancestral"].map(|k| b[k].as_str()),
        ["4", "3", "3", "0"],
        "1. B's counts"
    );
    assert_eq!(
        r.get("Ancestral_in_terminal"),
        ".",
        "1. nothing ancestral in the terminal block"
    );
    assert_eq!((r.get("Path_derived"), r.get("Path_ancestral")), ("7", "0"));
    assert_eq!(
        (r.get("Tree"), r.get("Positions")),
        ("YFull YTree test", "test tables"),
        "1. the result names its tree and tables"
    );
    assert!(
        r.stdout
            .contains("Current tree (YFull YTree test): C (level 3"),
        "1. the step line names the current tree: {}",
        r.stdout
    );

    // 2. one ancestral SNP in the terminal block
    let r = run(
        &d,
        &with(
            &ALL_A,
            &[("b1", St::D), ("b2", St::D), ("b3", St::A), ("x1", St::A)],
        ),
        "split",
        &[],
    );
    assert_eq!((r.get("Hg"), r.get("Ancestral_in_terminal")), ("B", "b3"));
    let last = r.path.last().unwrap();
    assert_eq!(
        ["Ancestral_SNPs", "Derived", "Ancestral"].map(|k| last[k].as_str()),
        ["b3", "2", "1"],
        "2. B's row lists b3 as ancestral"
    );

    // 3. no-calls
    let nocall = with(&ab, &[("c1", St::DAt(9)), ("c2", St::Mix), ("x1", St::A)]);
    let r = run(&d, &nocall, "nocall", &[]);
    assert_eq!(r.get("Hg"), "B", "3. depth 9 and a 50/50 mix are no-calls");
    assert_eq!(
        (r.path[1]["SNPs"].as_str(), r.path[1]["On_tables"].as_str()),
        ("4", "3"),
        "3. a SNP no table places is counted but not on the tables (bnot)"
    );
    let r2 = run(&d, &nocall, "nocall2", &["--min-depth", "5"]);
    assert!(
        r2.get("Hg") == "C" && r2.path.last().unwrap()["No_call"] == "1",
        "3. with --min-depth 5 the depth-9 read is derived: C, with one no-call"
    );

    // 4. a negative branch blocks the walk below it, unless outvoted
    let r = run(
        &d,
        &with(&ab, &[("x1", St::A), ("y1", St::D)]),
        "blocked",
        &[],
    );
    assert_eq!(
        (r.get("Hg"), r.get("Outvoted_ancestral")),
        ("B", "."),
        "4. one derived below: blocked"
    );
    let r = run(
        &d,
        // B has no read here: a clean B on another lineage would have more
        // support than Z, and win (i411_ytree_ranking.rs).
        &with(&ALL_A, &[("x1", St::A), ("y1", St::D), ("z1", St::D)]),
        "outvoted",
        &[],
    );
    assert_eq!(r.get("Hg"), "Z", "4. two derived below: Z");
    assert_eq!(r.branches(), ["A", "X", "Y", "Z"]);
    assert_eq!(
        r.get("Outvoted_ancestral"),
        "x1",
        "4. x1 is the outvoted site"
    );
    assert_eq!(r.path[1]["Ancestral"], "1");
    assert!(
        r.stdout.contains("outvoted") && r.stdout.contains("x1"),
        "4. the step line names the outvoted site: {}",
        r.stdout
    );
    let r = run(
        &d,
        &with(&ALL_A, &[("y1", St::D), ("z1", St::D)]),
        "unblocked",
        &[],
    );
    assert!(
        r.get("Hg") == "Z"
            && r.branches() == ["A", "X", "Y", "Z"]
            && r.get("Outvoted_ancestral") == ".",
        "4. X with no call: Z wins, nothing outvoted"
    );
    let r = run(
        &d,
        &with(
            &[("a1", St::A), ("a2alt", St::A)],
            &with(&ALL_B, &[("c1", St::D), ("c2", St::D)]),
        ),
        "twoanc",
        &[],
    );
    assert!(
        r.get("Hg") == "NA" && r.path.is_empty(),
        "4. two ancestral SNPs on A block B and C: Hg NA"
    );
    let r = run(
        &d,
        &with(
            &ab,
            &[("x1", St::A), ("y1", St::A), ("z1", St::D), ("w1", St::D)],
        ),
        "twosingles",
        &[],
    );
    assert_eq!(
        (r.get("Hg"), r.get("Outvoted_ancestral")),
        ("B", "."),
        "4. two single-ancestral branches on one path block Z and W"
    );
    let r = run(
        &d,
        &with(&ab, &[("p1", St::A), ("q1", St::D)]),
        "onesite",
        &[],
    );
    assert_eq!(
        (r.get("Hg"), r.get("Outvoted_ancestral")),
        ("B", "."),
        "4. two names at one derived site below are one site: Q is blocked"
    );

    // 5. nothing derived
    let r = run(&d, &[("a1", St::A), ("a2alt", St::A)], "none", &[]);
    assert!(
        r.get("Hg") == "NA" && r.get("Level") == "0" && r.path.is_empty(),
        "5. no derived call: Hg NA at level 0, an empty path"
    );
    assert!(
        r.stdout.contains("no branch has a derived call"),
        "5. the step line says so: {}",
        r.stdout
    );

    // 6. tables
    let out = d.join("pos.txt");
    let pargs = [
        "--tree-dir".to_owned(),
        d.display().to_string(),
        out.display().to_string(),
    ];
    let (code, _, err) = y("tree-positions", &pargs);
    assert_eq!(code, 0, "{err}");
    let text = fs::read_to_string(&out).unwrap();
    let pos: Vec<&str> = text.split_whitespace().collect();
    assert!(
        pos.contains(&"3000100") && !pos.contains(&"3999999"),
        "6. a1 at the first table's position: {pos:?}"
    );
    assert!(
        pos.contains(&"3000700"),
        "6. c2, on the second table only, is asked for"
    );
    assert!(
        !pos.contains(&"2000000"),
        "6. the pseudoautosomal SNP is never asked for"
    );
    assert!(
        text.lines().all(|l| l.split('\t').next() == Some("chrY")),
        "6. every line is chrY<TAB>pos"
    );

    // 7. the pileup reader: c1 (G>T) read as "^I..,$+2ACt-1Ga*" on
    // reference t is four T and one A, 80% T: a no-call at 90%; c2 is five
    // T of five: derived. Read as bases, the marks would change both.
    let p7 = d.join("marks.pileup");
    let mut text = String::new();
    for name in ["a1", "a2alt", "b1", "b2", "b3"] {
        let (pos, _, der) = pos_of(name);
        text += &format!(
            "chrY\t{pos}\tN\t12\t{}\t{}\n",
            der.to_string().repeat(12),
            "I".repeat(12)
        );
    }
    text += &format!("chrY\t{}\tt\t5\t^I..,$+2ACt-1Ga*\tIIIII\n", pos_of("c1").0);
    text += &format!("chrY\t{}\tt\t5\t^I.$,+2GG,t-1A.\tIIIII\n", pos_of("c2").0);
    fs::write(&p7, text).unwrap();
    let r = place(&d, &p7, "marks", &["--min-depth", "5"]);
    assert_eq!(r.code, 0, "{}", r.stderr);
    let rows = fs::read_to_string(d.join("marks_y_tree_path.tsv")).unwrap();
    let c_row: Vec<&str> = rows
        .lines()
        .map(|l| l.split('\t').collect::<Vec<_>>())
        .find(|x| x.get(1) == Some(&"C"))
        .map(|x| x[4..7].to_vec())
        .unwrap_or_default();
    assert_eq!(
        c_row,
        ["1", "0", "1"],
        "7. read marks and indels skipped, '.'/',' read as the reference"
    );

    // 8. a missing table
    fs::remove_file(d.join("positions_yfull_hg38.txt")).unwrap();
    let (code, _, err) = y("tree-positions", &pargs);
    assert!(
        code != 0 && err.contains("setup.sh --yleaf-data"),
        "a missing table stops the tool and names setup.sh --yleaf-data: {err}"
    );
    fs::remove_dir_all(d).unwrap();
}
