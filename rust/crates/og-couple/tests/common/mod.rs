//! What the og-couple tests share: a scratch folder, the repository's
//! files, and `og couple screen|genes` run in the process.

#![allow(dead_code)]

use std::collections::HashMap;
use std::path::PathBuf;

/// The repository root (the tests read its assets).
pub fn repo(rel: &str) -> String {
    format!("{}/../../../{rel}", env!("CARGO_MANIFEST_DIR"))
}

/// A scratch folder, removed when dropped.
pub struct Dir(PathBuf);

impl Dir {
    pub fn new(name: &str) -> Dir {
        let p = std::env::temp_dir().join(format!("og-couple-{name}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&p);
        std::fs::create_dir_all(&p).unwrap();
        Dir(p)
    }

    /// The path of `rel` inside the folder.
    pub fn path(&self, rel: &str) -> String {
        self.0.join(rel).to_string_lossy().into_owned()
    }

    /// Writes `text` to `rel` (parents made) and gives its path.
    pub fn put(&self, rel: &str, text: &str) -> String {
        let p = self.0.join(rel);
        std::fs::create_dir_all(p.parent().unwrap()).unwrap();
        std::fs::write(&p, text).unwrap();
        p.to_string_lossy().into_owned()
    }
}

impl Drop for Dir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// Runs `og couple TOOL ARGS` in the process: the exit code, and stdout
/// followed by stderr.
pub fn og(tool: &str, args: &[&str]) -> (u8, String) {
    let argv: Vec<String> = args.iter().map(|s| (*s).to_owned()).collect();
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let rc = match tool {
        "screen" => og_couple::screen::cli(&argv, &mut out, &mut err),
        "genes" => og_couple::genes::cli(&argv, &mut out, &mut err),
        _ => panic!("{tool}"),
    };
    out.extend(err);
    (rc, String::from_utf8(out).unwrap())
}

/// The text of a file, or "" when it is not there.
pub fn read(path: &str) -> String {
    std::fs::read_to_string(path).unwrap_or_default()
}

/// The rows of a TSV with a header, by column name.
pub fn rows(path: &str) -> Vec<HashMap<String, String>> {
    let text = read(path);
    let mut lines = text.lines();
    let Some(head) = lines.next() else {
        return Vec::new();
    };
    let head: Vec<&str> = head.split('\t').collect();
    lines
        .map(|l| {
            head.iter()
                .map(|h| (*h).to_owned())
                .zip(l.split('\t').map(str::to_owned))
                .collect()
        })
        .collect()
}

/// The rows of a couple table by gene.
pub fn by_gene(path: &str) -> HashMap<String, HashMap<String, String>> {
    rows(path)
        .into_iter()
        .map(|r| (r["gene"].clone(), r))
        .collect()
}

/// A cell of a row of `by_gene`, or "" when the gene or the cell is missing.
pub fn cell<'a>(t: &'a HashMap<String, HashMap<String, String>>, gene: &str, k: &str) -> &'a str {
    t.get(gene)
        .and_then(|r| r.get(k))
        .map_or("", String::as_str)
}

/// One row of assets/carrier_genes.tsv.
#[derive(Debug, Clone)]
pub struct Gene {
    pub chrom: String,
    pub start: i64,
    pub end: i64,
    pub inheritance: String,
    pub lists: String,
    pub condition: String,
}

/// The rows of a carrier list, by gene.
pub fn genes(path: &str) -> HashMap<String, Gene> {
    let text = read(path);
    let mut lines = text
        .lines()
        .filter(|l| !l.trim().is_empty() && !l.starts_with('#'));
    let head: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let at = |k: &str| head.iter().position(|h| *h == k).unwrap();
    let (g, c, s, e, i, l, n) = (
        at("gene"),
        at("chrom"),
        at("start"),
        at("end"),
        at("inheritance"),
        at("lists"),
        at("condition"),
    );
    lines
        .map(|line| {
            let f: Vec<&str> = line.split('\t').collect();
            (
                f[g].to_owned(),
                Gene {
                    chrom: f[c].to_owned(),
                    start: f[s].parse().unwrap(),
                    end: f[e].parse().unwrap(),
                    inheritance: f[i].to_owned(),
                    lists: f[l].to_owned(),
                    condition: f[n].to_owned(),
                },
            )
        })
        .collect()
}

/// The header line of a single-sample VCF.
pub const HEAD: &str =
    "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n";

/// The header of step 35's SMN table.
pub const PARASCOPY_HEAD: &str = "chrom\tstart\tend\tlocus\tagCN_filter\tagCN\tagCN_qual\tpsCN_filter\tpsCN\tpsCN_qual\thomologous_regions\n";
