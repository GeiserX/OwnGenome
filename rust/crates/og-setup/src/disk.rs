//! The disk checks setup makes before a large install, as
//! scripts/lib/common.sh's clinical_disk_ok and setup.sh's
//! annotation_score_disk_ok make them: the free space of the destination's
//! filesystem (`df -Pk`), and what an interrupted download already holds,
//! against the peak the install needs plus a 5 GB (5 GiB) margin.
//!
//! A free space that cannot be read is a failed check, never permission to
//! fill the disk. [`Probe`] reads the disk; the shell's `df` and `du` are
//! the default, so a fake `df` of the tests answers both the same way.

use crate::install::keep;
use std::io::Write;
use std::path::Path;
use std::process::{Command, Stdio};

/// How a check reads the disk.
pub trait Probe {
    /// `df -Pk DIR | awk 'NR == 2 {print $4}'`: free KiB, as text.
    fn free_kib(&self, dir: &Path) -> Option<String>;
    /// `du -k FILE | awk 'NR == 1 {print $1}'`: allocated KiB, as text.
    fn used_kib(&self, file: &Path) -> Option<String>;
}

/// The shell's df and du, from PATH.
pub struct Shell;

fn field(out: &[u8], line: usize, col: usize) -> Option<String> {
    String::from_utf8_lossy(out)
        .lines()
        .nth(line)
        .and_then(|l| l.split_whitespace().nth(col).map(str::to_owned))
}

impl Probe for Shell {
    fn free_kib(&self, dir: &Path) -> Option<String> {
        let o = Command::new("df")
            .arg("-Pk")
            .arg(dir)
            .stdin(Stdio::null())
            .stderr(Stdio::null())
            .output()
            .ok()?;
        field(&o.stdout, 1, 3)
    }

    fn used_kib(&self, file: &Path) -> Option<String> {
        let o = Command::new("du")
            .arg("-k")
            .arg(file)
            .stdin(Stdio::null())
            .stderr(Stdio::null())
            .output()
            .ok()?;
        field(&o.stdout, 0, 0)
    }
}

fn number(s: Option<String>) -> Option<u64> {
    let s = s?;
    (!s.is_empty() && s.bytes().all(|b| b.is_ascii_digit()))
        .then(|| s.parse().ok())
        .flatten()
}

/// clinical_disk_ok DIR PEAK LABEL [ARCHIVE]: whether DIR's filesystem has
/// PEAK GB plus a 5 GB margin free, counting what ARCHIVE (or ARCHIVE.part)
/// already holds. DIR is made, and a stopped extraction in it
/// (`.extract.*`) removed first. GENOME names the install folder in the
/// retry line. Messages go to OUT.
pub fn clinical_disk_ok(
    probe: &dyn Probe,
    dir: &Path,
    peak: u64,
    label: &str,
    archive: Option<&Path>,
    genome: &str,
    out: &mut dyn Write,
) -> bool {
    // mkdir -p and rm -rf as the shell runs them: a failure prints their message.
    if !keep::mkdir_p(&dir.to_string_lossy(), out) {
        return false;
    }
    if let Ok(rd) = std::fs::read_dir(dir) {
        for e in rd.flatten() {
            if e.file_name().to_string_lossy().starts_with(".extract.")
                && !keep::rm_rf(&e.path().to_string_lossy(), out)
            {
                return false;
            }
        }
    }
    let Some(free) = number(probe.free_kib(dir)) else {
        let _ = writeln!(
            out,
            "[WARN] Cannot read free space for {label}; run: ./scripts/setup.sh --clinical-data \"{genome}\""
        );
        return false;
    };
    let mut have = 0;
    if let Some(a) = archive {
        let part = std::path::PathBuf::from(format!("{}.part", a.display()));
        for f in [a, part.as_path()] {
            if f.is_file() {
                let Some(h) = number(probe.used_kib(f)) else {
                    let _ = writeln!(
                        out,
                        "[WARN] Cannot read allocated space for {label}; retry the clinical install."
                    );
                    return false;
                };
                have = h;
                break;
            }
        }
    }
    if free + have < (peak + 5) * 1_048_576 {
        let _ = writeln!(
            out,
            "[WARN] Skipping {label}: needs {peak} GB peak plus a 5 GB margin; {:.1} GB free.",
            free as f64 / 1_048_576.0
        );
        let _ = writeln!(
            out,
            "  Install later: ./scripts/setup.sh --clinical-data \"{genome}\""
        );
        return false;
    }
    true
}

/// One data file of the annotation scores: the bytes it still needs.
#[derive(Debug, Clone)]
pub struct Want {
    /// The file's full size in bytes.
    pub bytes: u64,
    /// The bytes on disk already: the whole file when its size matches, or
    /// an interrupted `.part`.
    pub have_whole: bool,
    pub part: u64,
}

/// The bytes still to download for these files.
pub fn remaining(files: &[Want]) -> u64 {
    files
        .iter()
        .filter(|w| !w.have_whole)
        .map(|w| w.bytes.saturating_sub(w.part))
        .sum()
}

/// annotation_score_disk_ok's test: REMAINING bytes plus a 5 GiB margin
/// fit in FREE KiB. RETRY is the command that retries; messages go to ERR.
pub fn annotation_scores_fit(
    remaining: u64,
    free: Option<String>,
    retry: &str,
    err: &mut dyn Write,
) -> bool {
    if remaining == 0 {
        return true;
    }
    let Some(free) = number(free) else {
        let _ = writeln!(
            err,
            "[WARN] Cannot read free space for annotation scores. Retry: {retry}"
        );
        return false;
    };
    if free < remaining.div_ceil(1024) + 5 * 1_048_576 {
        let _ = writeln!(
            err,
            "[WARN] Skipping annotation scores: needs {remaining} more bytes for CADD and SpliceAI plus a 5 GiB margin."
        );
        let _ = writeln!(err, "  Install later: {retry}");
        return false;
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    struct Fixed(Option<&'static str>, Option<&'static str>);
    impl Probe for Fixed {
        fn free_kib(&self, _: &Path) -> Option<String> {
            self.0.map(str::to_owned)
        }
        fn used_kib(&self, _: &Path) -> Option<String> {
            self.1.map(str::to_owned)
        }
    }

    fn dir(name: &str) -> std::path::PathBuf {
        std::env::temp_dir().join(format!("og-disk-{name}-{}", std::process::id()))
    }

    #[test]
    fn the_clinical_margin_is_five_gb_exactly() {
        let d = dir("margin");
        let need = (56 + 5) * 1_048_576;
        let mut out = Vec::new();
        let enough = need.to_string().leak();
        assert!(clinical_disk_ok(
            &Fixed(Some(enough), None),
            &d,
            56,
            "VEP 116",
            None,
            "/g",
            &mut out
        ));
        let short = (need - 1).to_string().leak();
        assert!(!clinical_disk_ok(
            &Fixed(Some(short), None),
            &d,
            56,
            "VEP 116",
            None,
            "/g",
            &mut out
        ));
        assert_eq!(
            String::from_utf8(out).unwrap(),
            "[WARN] Skipping VEP 116: needs 56 GB peak plus a 5 GB margin; 61.0 GB free.\n  Install later: ./scripts/setup.sh --clinical-data \"/g\"\n"
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn a_partial_archive_counts_towards_the_peak() {
        let d = dir("part");
        std::fs::create_dir_all(&d).unwrap();
        let a = d.join("cache.tar.gz");
        std::fs::write(d.join("cache.tar.gz.part"), b"x").unwrap();
        let mut out = Vec::new();
        let free = ((15 + 5) * 1_048_576 - 100).to_string().leak();
        assert!(clinical_disk_ok(
            &Fixed(Some(free), Some("100")),
            &d,
            15,
            "PCGR",
            Some(&a),
            "/g",
            &mut out
        ));
        assert!(!clinical_disk_ok(
            &Fixed(Some(free), Some("99")),
            &d,
            15,
            "PCGR",
            Some(&a),
            "/g",
            &mut out
        ));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn an_unreadable_df_fails_the_check() {
        let d = dir("nodf");
        let mut out = Vec::new();
        assert!(!clinical_disk_ok(
            &Fixed(None, None),
            &d,
            1,
            "VEP 115",
            None,
            "/g",
            &mut out
        ));
        assert!(!clinical_disk_ok(
            &Fixed(Some("12a"), None),
            &d,
            1,
            "VEP 115",
            None,
            "/g",
            &mut out
        ));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn the_annotation_scores_need_every_byte_plus_five_gib() {
        let rem = remaining(&[
            Want {
                bytes: 10_240,
                have_whole: false,
                part: 1_024,
            },
            Want {
                bytes: 99,
                have_whole: true,
                part: 0,
            },
        ]);
        assert_eq!(rem, 9_216);
        let fit = (9 + 5 * 1_048_576).to_string();
        let mut err = Vec::new();
        assert!(annotation_scores_fit(rem, Some(fit), "r", &mut err));
        let short = (8 + 5 * 1_048_576).to_string();
        assert!(!annotation_scores_fit(rem, Some(short), "r", &mut err));
        // One byte more than whole KiB needs one KiB more.
        let fit = (10 + 5 * 1_048_576).to_string();
        assert!(annotation_scores_fit(rem + 1, Some(fit), "r", &mut err));
        let short = (9 + 5 * 1_048_576).to_string();
        assert!(!annotation_scores_fit(rem + 1, Some(short), "r", &mut err));
        assert!(
            annotation_scores_fit(0, None, "r", &mut err),
            "nothing left to download needs no disk"
        );
    }
}
