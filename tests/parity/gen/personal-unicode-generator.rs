//! Generate the scanner's conservative Unicode-version guard from official UCD
//! files. Run personal-unicode-generate.sh, which checks each source's SHA-256.
//! No network or Python is used by this generator or by the scanner.
use std::{env, fs, path::Path};

#[derive(Clone, Default, PartialEq, Eq)]
struct Properties {
    category: String,
    astral_case: bool,
    space: bool,
    case: [String; 4], // simple upper/lower/title mappings and full case fold
}
impl Properties {
    fn assigned(&self) -> bool {
        !self.category.is_empty()
    }
    fn word(&self) -> bool {
        self.category.starts_with(['L', 'N'])
    }
    fn digit(&self) -> bool {
        self.category == "Nd"
    }
    fn sensitive(&self) -> bool {
        self.category.starts_with(['L', 'M', 'N'])
            || self.category == "Pc"
            || self.space
            || self.case.iter().any(|s| !s.is_empty())
    }
    fn signature(&self) -> (bool, bool, bool, &[String; 4]) {
        (self.word(), self.digit(), self.space, &self.case)
    }
}
fn read(dir: &Path, version: &str) -> Vec<Properties> {
    let mut table = vec![Properties::default(); 0x110000];
    let data = fs::read_to_string(dir.join(version).join("UnicodeData.txt")).unwrap();
    let mut start = None;
    for line in data.lines() {
        let f: Vec<_> = line.split(';').collect();
        assert_eq!(f.len(), 15);
        let c = usize::from_str_radix(f[0], 16).unwrap();
        if f[1].ends_with(", First>") {
            assert!(start.replace(c).is_none());
            continue;
        }
        let a = if f[1].ends_with(", Last>") {
            start.take().unwrap()
        } else {
            c
        };
        for p in &mut table[a..=c] {
            assert!(!p.assigned());
            *p = Properties {
                category: f[2].into(),
                astral_case: false,
                space: f[2] == "Zs" || matches!(f[4], "WS" | "B" | "S"),
                case: [f[12].into(), f[13].into(), f[14].into(), String::new()],
            };
        }
    }
    assert!(start.is_none());
    // The astral rule uses every mapping, including full mappings and their
    // targets. It does not depend on any interpreter or regex pattern shape.
    for line in data.lines() {
        let f: Vec<_> = line.split(';').collect();
        let c = usize::from_str_radix(f[0], 16).unwrap();
        for mapping in &f[12..15] {
            mark_astral_case(&mut table, c, mapping);
        }
    }
    let special = fs::read_to_string(dir.join(version).join("SpecialCasing.txt")).unwrap();
    assert!(special.contains(&format!("SpecialCasing-{version}.txt")));
    for line in special.lines() {
        let line = line.split('#').next().unwrap().trim();
        if line.is_empty() {
            continue;
        }
        let f: Vec<_> = line.split(';').map(str::trim).collect();
        let c = usize::from_str_radix(f[0], 16).unwrap();
        for mapping in &f[1..4] {
            mark_astral_case(&mut table, c, mapping);
        }
    }
    let folding = fs::read_to_string(dir.join(version).join("CaseFolding.txt")).unwrap();
    assert!(folding.contains(&format!("CaseFolding-{version}.txt")));
    for line in folding.lines() {
        let line = line.split('#').next().unwrap().trim();
        if line.is_empty() {
            continue;
        }
        let f: Vec<_> = line.split(';').map(str::trim).collect();
        let source = usize::from_str_radix(f[0], 16).unwrap();
        mark_astral_case(&mut table, source, f[2]);
        if matches!(f[1], "C" | "F") {
            let c = usize::from_str_radix(f[0], 16).unwrap();
            assert!(table[c].assigned());
            table[c].case[3] = f[2].into();
        }
    }
    table
}
fn mark_astral_case(table: &mut [Properties], source: usize, mapping: &str) {
    for target in mapping.split_whitespace() {
        let target = usize::from_str_radix(target, 16).unwrap();
        if target != source {
            if source > 0xffff {
                table[source].astral_case = true;
            }
            if target > 0xffff {
                table[target].astral_case = true;
            }
        }
    }
}
fn main() {
    let args: Vec<_> = env::args().collect();
    assert_eq!(args.len(), 3, "usage: generator UCD_DIR OUTPUT");
    let versions = ["13.0.0", "15.0.0", "16.0.0", "17.0.0"];
    let tables: Vec<_> = versions
        .iter()
        .map(|v| read(Path::new(&args[1]), v))
        .collect();
    let mut refused = vec![false; 0x110000];
    for (c, bad) in refused.iter_mut().enumerate() {
        if char::from_u32(c as u32).is_none() {
            continue;
        }
        let old = &tables[0][c];
        *bad = tables.iter().any(|t| t[c].astral_case)
            || !tables[3][c].assigned()
            || (!old.assigned() && tables.iter().any(|t| t[c].sensitive()))
            || (old.assigned() && tables.iter().any(|t| t[c].signature() != old.signature()));
    }
    // Include both ends of a changed case relation. This also guards a stable
    // letter which acquires a newly assigned case partner in a later version.
    loop {
        let before = refused.clone();
        for table in &tables {
            for (c, p) in table.iter().enumerate() {
                for mapping in &p.case {
                    if mapping.is_empty() {
                        continue;
                    }
                    // Full folds can expand to several scalars. They affect
                    // the signature above, but are not regex case relations.
                    if mapping.contains(' ') {
                        continue;
                    }
                    let target = usize::from_str_radix(mapping, 16).unwrap();
                    if refused[c] || refused[target] {
                        refused[c] = true;
                        refused[target] = true;
                    }
                }
            }
        }
        if refused == before {
            break;
        }
    }
    let mut out = String::from(
        "// Generated by personal-unicode-generator.rs; do not edit.\n// UCD 13.0.0, 15.0.0, 16.0.0, 17.0.0; sources pinned in personal-unicode-sources.tsv.\nconst VERSION_DEPENDENT: &[(u32, u32)] = &[\n",
    );
    let mut c = 0;
    while c < refused.len() {
        if !refused[c] {
            c += 1;
            continue;
        }
        let start = c;
        while c + 1 < refused.len() && refused[c + 1] {
            c += 1;
        }
        out.push_str(&format!("    (0x{start:x}, 0x{c:x}),\n"));
        c += 1;
    }
    out.push_str("];\n");
    fs::write(&args[2], out).unwrap();
    println!(
        "generated guard for {} Unicode scalars",
        refused.iter().filter(|&&x| x).count()
    );
}
