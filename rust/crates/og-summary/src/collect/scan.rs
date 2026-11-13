//! [`FolderIndex`]: one walk of a sample folder, then every section's file
//! found in it by pattern.
//!
//! The Python collector called `glob` once per candidate path of every
//! section (`first_existing`), about 100 directory listings per summary.
//! Here the folder is listed once, and each candidate is matched against
//! that list in the order Python tried them: candidate order, then sorted
//! path within one candidate, the first regular file wins.

use pycompat::mtime::{fnmatch, py_cmp};
use std::collections::HashSet;
use std::fs;
use std::path::{Path, PathBuf};

/// How deep the walk goes below the sample folder. The deepest file a
/// section reads today is three folders down; Nextflow's work folders are
/// not inside a sample folder.
pub const MAX_DEPTH: usize = 6;

/// Every regular file of a sample folder, by its path relative to the
/// folder (`/` between parts). Symlinks are followed, so the links the
/// HTML_REPORT task stages count as the files they point to; the depth cap
/// ends a symlink loop.
#[derive(Debug, Default)]
pub struct FolderIndex {
    files: HashSet<String>,
    /// The same paths, sorted as Python sorts them, for wildcard patterns.
    sorted: Vec<String>,
}

impl FolderIndex {
    /// Walks `root`. A folder that cannot be read counts as empty, as it
    /// does for `glob`; a missing root gives an empty index. A name that is
    /// not UTF-8 is left out: no pattern of a section can name it.
    pub fn walk(root: &Path) -> FolderIndex {
        let mut files = HashSet::new();
        let mut stack: Vec<(PathBuf, String, usize)> = vec![(root.to_path_buf(), String::new(), 0)];
        while let Some((dir, rel, depth)) = stack.pop() {
            let Ok(entries) = fs::read_dir(&dir) else {
                continue;
            };
            for e in entries.flatten() {
                let Ok(name) = e.file_name().into_string() else {
                    continue;
                };
                let path = e.path();
                let Ok(meta) = fs::metadata(&path) else {
                    continue; // a broken symlink: no file, as os.path.isfile says
                };
                let rel = if rel.is_empty() {
                    name
                } else {
                    format!("{rel}/{name}")
                };
                if meta.is_file() {
                    files.insert(rel);
                } else if meta.is_dir() && depth + 1 < MAX_DEPTH {
                    stack.push((path, rel, depth + 1));
                }
            }
        }
        let mut sorted: Vec<String> = files.iter().cloned().collect();
        sorted.sort_by(|a, b| py_cmp(a, b));
        FolderIndex { files, sorted }
    }

    /// `first_existing(base, patterns)`: the first file the patterns name,
    /// relative to the folder. A pattern is a relative path whose parts may
    /// hold `*`, `?` and `[...]`; a part that is exactly `**` stands for any
    /// number of folders. As in `glob`, a wildcard part never matches a
    /// name starting with `.`.
    pub fn first(&self, patterns: &[String]) -> Option<&str> {
        for p in patterns {
            if !has_magic(p) {
                if let Some(f) = self.files.get(p.as_str()) {
                    return Some(f);
                }
                continue;
            }
            let pat: Vec<&str> = p.split('/').collect();
            if let Some(f) = self.sorted.iter().find(|f| {
                let parts: Vec<&str> = f.split('/').collect();
                matches(&parts, &pat)
            }) {
                return Some(f);
            }
        }
        None
    }
}

fn has_magic(s: &str) -> bool {
    s.contains(['*', '?', '['])
}

/// One path, in parts, against one pattern, in parts.
fn matches(parts: &[&str], pat: &[&str]) -> bool {
    match (pat.first(), parts.first()) {
        (None, None) => true,
        (None, Some(_)) => false,
        (Some(&"**"), _) => {
            // Zero folders, or one more (not a hidden one) and try again.
            matches(parts, &pat[1..])
                || (parts.len() > 1 && !parts[0].starts_with('.') && matches(&parts[1..], pat))
        }
        (Some(_), None) => false,
        (Some(&p), Some(&name)) => {
            let ok = if has_magic(p) {
                (!name.starts_with('.') || p.starts_with('.')) && fnmatch(name, p)
            } else {
                name == p
            };
            ok && matches(&parts[1..], &pat[1..])
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testutil::Dir;

    fn pats(p: &[&str]) -> Vec<String> {
        p.iter().map(|s| (*s).to_owned()).collect()
    }

    #[test]
    fn candidates_in_order_then_sorted_paths() {
        let d = Dir::new("scan-order");
        d.put("indexcov/B-indexcov.ped", "x");
        d.put("indexcov/A-indexcov.ped", "x");
        d.put("indexcov/.C-indexcov.ped", "x");
        d.put("coverage/S.mosdepth.summary.txt", "x");
        std::fs::create_dir_all(d.path("mosdepth/S.mosdepth.summary.txt")).unwrap();
        let ix = FolderIndex::walk(&d.0);
        // A folder of the first candidate's name is no file.
        assert_eq!(
            ix.first(&pats(&[
                "mosdepth/S.mosdepth.summary.txt",
                "coverage/S.mosdepth.summary.txt"
            ])),
            Some("coverage/S.mosdepth.summary.txt")
        );
        assert_eq!(
            ix.first(&pats(&[
                "indexcov/indexcov-indexcov.ped",
                "indexcov/*-indexcov.ped"
            ])),
            Some("indexcov/A-indexcov.ped")
        );
        assert_eq!(
            ix.first(&pats(&["indexcov/.C-indexcov.ped"])),
            Some("indexcov/.C-indexcov.ped")
        );
        assert_eq!(ix.first(&pats(&["nothing/*"])), None);
    }

    #[test]
    fn wildcards_skip_hidden_names_and_recurse() {
        let d = Dir::new("scan-hidden");
        d.put("indexcov/.x-indexcov.ped", "x");
        d.put("manta/a/b/diploidSV.vcf.gz", "x");
        d.put("manta/.h/diploidSV.vcf.gz", "x");
        let ix = FolderIndex::walk(&d.0);
        assert_eq!(ix.first(&pats(&["indexcov/*-indexcov.ped"])), None);
        assert_eq!(
            ix.first(&pats(&["manta/**/diploidSV.vcf.gz"])),
            Some("manta/a/b/diploidSV.vcf.gz")
        );
    }

    #[test]
    fn symlinks_count_as_their_targets() {
        let d = Dir::new("scan-links");
        d.put("store/qc/S_read_count.tsv", "x");
        d.put("store/vcf.gz", "x");
        std::os::unix::fs::symlink("store/qc", d.path("qc")).unwrap();
        std::fs::create_dir_all(d.path("vcf")).unwrap();
        std::os::unix::fs::symlink("../store/vcf.gz", d.path("vcf/S.vcf.gz")).unwrap();
        std::os::unix::fs::symlink("nowhere", d.path("vcf/T.vcf.gz")).unwrap();
        // A loop ends at the depth cap.
        std::os::unix::fs::symlink(".", d.path("store/loop")).unwrap();
        let ix = FolderIndex::walk(&d.0);
        assert_eq!(
            ix.first(&pats(&["qc/S_read_count.tsv"])),
            Some("qc/S_read_count.tsv")
        );
        assert_eq!(ix.first(&pats(&["vcf/S.vcf.gz"])), Some("vcf/S.vcf.gz"));
        assert_eq!(ix.first(&pats(&["vcf/T.vcf.gz"])), None);
    }

    #[test]
    fn a_missing_folder_is_empty() {
        let ix = FolderIndex::walk(Path::new("/nonexistent/og-summary-scan"));
        assert_eq!(ix.first(&pats(&["a"])), None);
    }
}
