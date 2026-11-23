//! `og lookup` runs inside `docker run`, which copies the container's stdout
//! and stderr on separate threads. An error line written to stderr in pieces
//! (`ERROR: `, then the message, then the newline) let the blank line that
//! ends each answer on stdout land between them, so `ERROR: rs… was not
//! found` sometimes reached the log split over two lines. Each stderr line
//! must leave in one write.

use std::io::{self, Write};
use std::path::Path;

fn input(name: &str) -> String {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../../tests/parity/inputs/lookup")
        .join(name)
        .to_string_lossy()
        .into_owned()
}

/// Keeps every `write` call apart.
#[derive(Default)]
struct Writes(Vec<Vec<u8>>);

impl Write for Writes {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.0.push(buf.to_vec());
        Ok(buf.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

fn run(args: &[&str]) -> (u8, Vec<String>) {
    let args: Vec<String> = args.iter().map(|s| (*s).to_owned()).collect();
    let (mut out, mut err) = (Vec::new(), Writes::default());
    let rc = og_lookup::cli(&args, &mut out, &mut err);
    let writes = err
        .0
        .into_iter()
        .map(|w| String::from_utf8(w).unwrap())
        .collect();
    (rc, writes)
}

#[test]
fn an_unknown_rsid_error_is_one_write() {
    let (sites, genes, traits, gvcf) = (
        input("sites.tsv.gz"),
        input("gencode.v50.basic.genes.gtf"),
        input("traits.tsv"),
        input("g.vcf.gz"),
    );
    let (rc, writes) = run(&[
        "query", "--sites", &sites, "--genes", &genes, "--table", &traits, "--gvcf", &gvcf, "--",
        "rs999",
    ]);
    assert_eq!(rc, 3, "{writes:?}");
    let line = writes
        .iter()
        .find(|w| w.contains("ERROR"))
        .unwrap_or_else(|| panic!("no ERROR write: {writes:?}"));
    assert!(
        line.starts_with("ERROR: rs999 was not found in ") && line.ends_with('\n'),
        "the whole error line leaves in one write: {writes:?}"
    );
}

#[test]
fn a_usage_error_is_one_write() {
    let (rc, writes) = run(&["nosuchverb"]);
    assert_eq!(rc, 2, "{writes:?}");
    assert_eq!(
        writes.len(),
        1,
        "usage and error leave in one write: {writes:?}"
    );
    assert!(
        writes[0].contains("og lookup: error: unknown lookup verb"),
        "{writes:?}"
    );
}
