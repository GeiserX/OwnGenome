//! The CYP2D6 depth check needs real depth in the flanks (#411, pgx-09).
//!
//! `og pgx cyp2d6-check` looked only at the share of reads with MAPQ 0, so
//! flanks at 6x with every read mapped uniquely passed as "fit for a
//! copy-number call", and step 36 then counted pypgx and Cyrius agreeing on
//! those few reads as a confirmed CYP2D6 call. The flanks must now have
//! 10x or more at MAPQ >= 1, the reads the callers use, before the check
//! passes; under that the status is `unreliable` and the message, which step
//! 36 writes as its reason, says `agreement not counted (depth)`.
//!
//! These cases pin:
//!
//! 1. flanks at 6x, cleanly mapped: `unreliable`, `agreement not counted
//!    (depth)`;
//! 2. flanks at 12x: `ok`, the agreement counts;
//! 3. the edge: 9.99x at MAPQ >= 1 fails, 10x passes, and the floor reads
//!    the MAPQ >= 1 depth, not the depth of all reads;
//! 4. multi-mapped reads on shallow flanks keep the multi-mapped message,
//!    which names the cause the user can fix.

use std::fs;
use std::path::{Path, PathBuf};

use og_bam::cyp2d6::{MIN_FLANK_DEPTH, cli};

const OK: &str = "CYP2D6 depth is fit for a copy-number call";
const NOT_COUNTED: &str = "agreement not counted (depth)";
const MULTI: &str = "reads are multi-mapped";

struct Scratch(PathBuf);

impl Scratch {
    fn new(name: &str) -> Self {
        let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!("i411-cyp2d6-gate-{name}"));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        Scratch(dir)
    }

    fn path(&self, name: &str) -> String {
        self.0.join(name).display().to_string()
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

/// A mosdepth regions file over the regions `cyp2d6-bed` prints.
fn regions(gene: f64, flank: f64) -> String {
    format!(
        "chr22\t42050000\t42100000\tflank\t{flank:.2}\n\
         chr22\t42123192\t42132032\tCYP2D6\t{gene:.2}\n\
         chr22\t42200000\t42250000\tflank\t{flank:.2}\n"
    )
}

/// `(status, message)` of `og pgx cyp2d6-check` on the two depth tables,
/// each `(gene, flank)` mean depth, all reads first, then MAPQ >= 1.
fn check(name: &str, all: (f64, f64), q1: (f64, f64)) -> (String, String) {
    let s = Scratch::new(name);
    let (q0, q1p, out) = (s.path("q0.bed"), s.path("q1.bed"), s.path("check.tsv"));
    fs::write(&q0, regions(all.0, all.1)).unwrap();
    fs::write(&q1p, regions(q1.0, q1.1)).unwrap();
    let argv: Vec<String> = ["--all", &q0, "--mapq1", &q1p, "--out", &out]
        .iter()
        .map(|a| (*a).to_owned())
        .collect();
    let (mut o, mut e) = (Vec::new(), Vec::new());
    let code = cli("check", &argv, &mut o, &mut e);
    assert_eq!(code, 0, "{name}: {}", String::from_utf8_lossy(&e));
    let text = fs::read_to_string(&out).unwrap();
    let get = |k: &str| {
        text.lines()
            .find_map(|l| l.strip_prefix(&format!("{k}\t")))
            .unwrap_or_else(|| panic!("{name}: no {k} row in {text}"))
            .to_owned()
    };
    (get("status"), get("message"))
}

#[test]
fn flanks_at_6x_do_not_count_the_agreement() {
    let (status, message) = check("6x", (4.0, 6.0), (3.8, 6.0));
    assert_eq!(status, "unreliable");
    assert!(message.contains(NOT_COUNTED), "{message}");
    assert!(message.contains("under 10x at MAPQ >= 1"), "{message}");
}

#[test]
fn flanks_at_12x_count() {
    let (status, message) = check("12x", (7.5, 12.0), (7.0, 12.0));
    assert_eq!(status, "ok");
    assert_eq!(message, OK);
}

#[test]
fn the_floor_is_10x_at_mapq1() {
    assert!((MIN_FLANK_DEPTH - 10.0).abs() < f64::EPSILON);
    let (status, message) = check("9.99x", (6.0, 9.99), (6.0, 9.99));
    assert_eq!(status, "unreliable", "9.99x: {message}");
    assert!(message.contains(NOT_COUNTED), "{message}");
    let (status, _) = check("10x", (6.0, 10.0), (6.0, 10.0));
    assert_eq!(status, "ok", "10x at MAPQ >= 1 counts");
    // 11x of all reads, 9x of them at MAPQ >= 1 (a share of 0.82, not
    // multi-mapped): the callers see 9x, so the agreement does not count.
    let (status, message) = check("mapq1", (7.0, 11.0), (6.5, 9.0));
    assert_eq!(status, "unreliable", "{message}");
    assert!(message.contains(NOT_COUNTED), "{message}");
}

#[test]
fn multi_mapped_reads_are_named_before_the_depth() {
    // Half the flank reads at MAPQ 0, and 4x left at MAPQ >= 1.
    let (status, message) = check("multi", (5.0, 8.0), (2.0, 4.0));
    assert_eq!(status, "unreliable");
    assert!(message.contains(MULTI), "{message}");
    assert!(!message.contains(NOT_COUNTED), "{message}");
}
