//! The ALT contig refusal of `og cn check-bam` and `og cn call
//! --bam-header`: a header with an ALT or an HLA contig is refused (exit 4,
//! a REFUSED line), a header without one passes (the positive control).
//!
//! This was part 5 of `tests/test_callability.py`; its other parts are in
//! `rust/crates/og-render/tests/clinical.rs`. The Python fed the contig
//! names on stdin (`--header -`); in the process they come from a file,
//! which `check_bam` reads with the same line rule.

const REPO: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../..");

/// `og cn ARGS` in the process: (exit code, stderr).
fn og_cn(args: &[&str]) -> (u8, String) {
    let args: Vec<String> = args.iter().map(|a| (*a).to_owned()).collect();
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = og_cn::cli::main(&args, &mut out, &mut err);
    (code, String::from_utf8(err).unwrap())
}

/// A scratch folder, removed when dropped.
struct Dir(std::path::PathBuf);

impl Dir {
    fn new(name: &str) -> Dir {
        let p = std::env::temp_dir().join(format!("og-cn-{name}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&p);
        std::fs::create_dir_all(&p).unwrap();
        Dir(p)
    }

    fn put(&self, rel: &str, text: &str) -> String {
        let p = self.0.join(rel);
        std::fs::write(&p, text).unwrap();
        p.to_string_lossy().into_owned()
    }
}

impl Drop for Dir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

#[test]
fn a_header_with_an_alt_or_hla_contig_is_refused() {
    let t = Dir::new("alt-refusal");
    let noalt = t.put(
        "noalt.txt",
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr16\tLN:90338345\n@SQ\tSN:chrM\tLN:16569\n",
    );
    let alt = t.put(
        "alt.txt",
        "@SQ\tSN:chr16\tLN:90338345\n@SQ\tSN:chr22_KI270879v1_alt\tLN:304135\n@RG\tID:x\tSM:y\n",
    );
    let (rc, err) = og_cn(&["check-bam", "--header", &noalt]);
    assert_eq!(
        (rc, err.as_str()),
        (0, ""),
        "positive control: no ALT contig passes"
    );
    let (rc, err) = og_cn(&["check-bam", "--header", &alt]);
    assert_eq!(rc, 4, "{err}");
    assert!(
        err.starts_with("REFUSED:") && err.contains("chr22_KI270879v1_alt"),
        "{err}"
    );
    let names = t.put("names.txt", "chr1\nHLA-A*01:01:01:01\n");
    let (rc, err) = og_cn(&["check-bam", "--header", &names]);
    assert!(
        rc == 4 && err.contains("HLA-A"),
        "contig names with an HLA contig: {rc} {err}"
    );

    let bed = t.put("bedcov.tsv", "");
    let catalog = format!("{REPO}/bin/cn_loci.tsv");
    let call = |header: &str| {
        og_cn(&[
            "call",
            "--catalog",
            &catalog,
            "--locus",
            "HBA",
            "--bedcov",
            &bed,
            "--bedcov-mapq1",
            &bed,
            "--bam-header",
            header,
        ])
    };
    let (rc, err) = call(&alt);
    assert!(
        rc == 4 && err.starts_with("REFUSED:"),
        "call refuses before any depth: {rc} {err}"
    );
    let (rc, err) = call(&noalt);
    assert_eq!(
        rc, 0,
        "two empty tables: not assessed, and no refusal: {err}"
    );
}
