//! `og dev install-manifest`: a manifest of an install tree, and the
//! compare of two manifests that the Rust port of setup must pass
//! (`tests/setup/ab.sh`).
//!
//! 1. A manifest lists every file, folder and symlink with its type, size,
//!    sha256, mode and the digests its compare modes read, sorted, with no
//!    trace of where the tree sits; the `downloads` folder at the top,
//!    `*.log` files and `--exclude` paths are left out.
//! 2. Two copies of one tree compare green. One flipped byte, a removed
//!    file, an added file, a changed mode, a removed folder, a changed
//!    symlink target and a changed type each turn it red, with the path.
//! 3. Each compare mode lets through only the bytes it is for, and a real
//!    change of content inside a file of its kind is red: gzip (header time
//!    vs content), vcf (bcftools Date and fileDate vs a record or the
//!    bcftools version), index (offsets vs sequence names, and an index that
//!    does not fit its data), jellyfish (header vs counts), numeric (float
//!    noise vs a fourth digit), skip and perm.
//! 4. A stamp is compared byte for byte: a row that can match one is
//!    refused, and no row applies to one.
//! 5. A bad allowlist row, a file that is not a manifest and a missing tree
//!    are refused with exit 2; the committed allowlist is read cleanly.
//! 6. The control: the smallest file compared byte for byte is picked, one
//!    flipped bit is red and a second flip restores it. The log compare
//!    masks times and progress and is red on another missing list.

use flate2::{Compression, GzBuilder};
use og_dev::install_manifest::cli;
use std::fs;
use std::io::Write;
use std::os::unix::fs::{PermissionsExt, symlink};
use std::path::{Path, PathBuf};
use std::process::Command;

fn repo() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../..")
}

fn s(p: &Path) -> String {
    p.to_str().expect("a UTF-8 path").to_owned()
}

/// A fresh folder under the system temp folder, removed on drop.
struct Tmp(PathBuf);

impl Tmp {
    fn new(tag: &str) -> Tmp {
        use std::sync::atomic::{AtomicUsize, Ordering};
        static N: AtomicUsize = AtomicUsize::new(0);
        let p = std::env::temp_dir().join(format!(
            "og-install-manifest-{tag}-{}-{}",
            std::process::id(),
            N.fetch_add(1, Ordering::Relaxed)
        ));
        let _ = fs::remove_dir_all(&p);
        fs::create_dir_all(&p).unwrap();
        Tmp(p)
    }
}

impl Drop for Tmp {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn write(p: &Path, body: &[u8], mode: u32) {
    fs::create_dir_all(p.parent().unwrap()).unwrap();
    fs::write(p, body).unwrap();
    fs::set_permissions(p, fs::Permissions::from_mode(mode)).unwrap();
}

/// A small install tree: three files, an empty folder, a symlink, and what
/// a manifest leaves out (the downloads folder and a log).
fn tree(root: &Path) {
    write(&root.join("reference/x.txt"), b"hello", 0o644);
    write(&root.join("reference/sub/y.bin"), &[0, 1, 2, 3], 0o600);
    write(&root.join("tools/run.sh"), b"#!/bin/sh\n", 0o755);
    fs::create_dir_all(root.join("reference/empty")).unwrap();
    for d in ["reference", "reference/sub", "reference/empty", "tools"] {
        fs::set_permissions(root.join(d), fs::Permissions::from_mode(0o755)).unwrap();
    }
    symlink("sub/y.bin", root.join("reference/y.link")).unwrap();
    write(&root.join("downloads/big.gz"), b"raw", 0o644);
    write(
        &root.join("reference/sub/build.log"),
        b"built at 12:00",
        0o644,
    );
}

fn copy(from: &Path, to: &Path) {
    let st = Command::new("cp")
        .arg("-pR")
        .arg(from)
        .arg(to)
        .status()
        .unwrap();
    assert!(st.success(), "cp -pR failed");
}

fn run(args: &[&str]) -> (u8, String, String) {
    let args: Vec<String> = args.iter().map(|a| (*a).to_owned()).collect();
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = cli(&args, &mut out, &mut err);
    (
        code,
        String::from_utf8(out).unwrap(),
        String::from_utf8(err).unwrap(),
    )
}

/// The manifest of DIR, written to DIR's sibling DIR.tsv; its path.
fn manifest(dir: &Path, extra: &[&str]) -> PathBuf {
    let mut args = vec![s(dir)];
    args.extend(extra.iter().map(|a| (*a).to_owned()));
    let a: Vec<&str> = args.iter().map(String::as_str).collect();
    let (code, out, err) = run(&a);
    assert_eq!(code, 0, "install-manifest {dir:?} failed: {err}");
    let p = dir.with_extension("tsv");
    fs::write(&p, out).unwrap();
    p
}

fn compare(a: &Path, b: &Path, allow: Option<&Path>) -> (u8, String, String) {
    let (a, b) = (s(a), s(b));
    let mut args = vec!["--compare", &a, &b];
    let allow = allow.map(s);
    if let Some(al) = &allow {
        args.extend(["--allow", al]);
    }
    run(&args)
}

/// Two copies of `tree`, A and B, under one temp folder.
fn pair(tag: &str) -> (Tmp, PathBuf, PathBuf) {
    let t = Tmp::new(tag);
    let a = t.0.join("a");
    fs::create_dir_all(&a).unwrap();
    tree(&a);
    let b = t.0.join("b");
    copy(&a, &b);
    (t, a, b)
}

fn diff_lines(out: &str) -> Vec<&str> {
    out.lines().filter(|l| l.starts_with("DIFF\t")).collect()
}

const HEAD: &str = "path\tmode\tcomponent\treason\tcannot_compare\n";

fn allow_file(t: &Tmp, rows: &str) -> PathBuf {
    let p = t.0.join("allow.tsv");
    fs::write(&p, format!("# rows of this test\n{HEAD}{rows}")).unwrap();
    p
}

#[test]
fn manifest_lists_entries_and_leaves_out_downloads_and_logs() {
    let t = Tmp::new("list");
    let a = t.0.join("a");
    fs::create_dir_all(&a).unwrap();
    tree(&a);
    let text = fs::read_to_string(manifest(&a, &[])).unwrap();
    let rows: Vec<&str> = text.lines().filter(|l| !l.starts_with('#')).collect();
    assert_eq!(rows[0], "path\ttype\tsize\tsha256\tmode\tdigests", "{text}");
    let want = [
        "reference\td\t-\t-\t0755\t-",
        "reference/empty\td\t-\t-\t0755\t-",
        "reference/sub\td\t-\t-\t0755\t-",
        "reference/sub/y.bin\tf\t4\t054edec1d0211f624fed0cbca9d4f9400b0e491c43742af2c5b0abebf0c990d8\t0600\t-",
        "reference/x.txt\tf\t5\t2cf24dba5fb0a30e26e83b2ac5b9e29e1b161e5c1fa7425e73043362938b9824\t0644\t-",
        "reference/y.link\tl\t-\t->sub/y.bin\t-\t-",
        "tools\td\t-\t-\t0755\t-",
    ];
    assert_eq!(rows.len(), want.len() + 2, "{text}");
    assert_eq!(&rows[1..8], &want, "{text}");
    assert!(
        rows[8].starts_with("tools/run.sh\tf\t10\t") && rows[8].ends_with("\t0755\t-"),
        "{text}"
    );
    assert!(
        !text.contains("downloads"),
        "the downloads folder is listed:\n{text}"
    );
    assert!(!text.contains("build.log"), "a log is listed:\n{text}");
    assert!(
        !text.contains(&s(&a)),
        "the manifest names where the tree sits:\n{text}"
    );
    let ex = fs::read_to_string(manifest(&a, &["--exclude", "reference/sub"])).unwrap();
    assert!(!ex.contains("reference/sub"), "{ex}");
    assert!(ex.contains("reference/x.txt"), "{ex}");
}

#[test]
fn two_copies_compare_green_and_the_manifests_are_equal() {
    let (_t, a, b) = pair("copies");
    let (ma, mb) = (manifest(&a, &[]), manifest(&b, &[]));
    assert_eq!(
        fs::read_to_string(&ma).unwrap(),
        fs::read_to_string(&mb).unwrap()
    );
    let (code, out, err) = compare(&ma, &mb, None);
    assert_eq!(code, 0, "{out}{err}");
    assert!(diff_lines(&out).is_empty(), "{out}");
}

/// Applies `change` to B, then expects a red compare whose one difference
/// line names `kind` and `path`.
fn expect_red(tag: &str, change: impl Fn(&Path), kind: &str, path: &str) {
    let (_t, a, b) = pair(tag);
    change(&b);
    let (code, out, err) = compare(&manifest(&a, &[]), &manifest(&b, &[]), None);
    assert_eq!(code, 1, "{tag}: not red\n{out}{err}");
    let d = diff_lines(&out);
    assert_eq!(d.len(), 1, "{tag}: {out}");
    assert!(
        d[0].starts_with(&format!("DIFF\t{kind}\t{path}\t")),
        "{tag}: {out}"
    );
}

#[test]
fn a_flipped_byte_is_red_with_the_path() {
    expect_red(
        "byte",
        |b| write(&b.join("reference/sub/y.bin"), &[0, 1, 2, 4], 0o600),
        "content",
        "reference/sub/y.bin",
    );
}

#[test]
fn a_removed_file_is_red_with_the_path() {
    expect_red(
        "removed",
        |b| fs::remove_file(b.join("reference/x.txt")).unwrap(),
        "removed",
        "reference/x.txt",
    );
}

#[test]
fn an_added_file_is_red_with_the_path() {
    expect_red(
        "added",
        |b| write(&b.join("tools/new.txt"), b"new", 0o644),
        "added",
        "tools/new.txt",
    );
}

#[test]
fn a_changed_mode_is_red_with_the_path() {
    expect_red(
        "mode",
        |b| fs::set_permissions(b.join("tools/run.sh"), fs::Permissions::from_mode(0o644)).unwrap(),
        "mode",
        "tools/run.sh",
    );
}

#[test]
fn a_removed_folder_a_new_link_target_and_a_new_type_are_red() {
    expect_red(
        "dir",
        |b| fs::remove_dir(b.join("reference/empty")).unwrap(),
        "removed",
        "reference/empty",
    );
    expect_red(
        "link",
        |b| {
            fs::remove_file(b.join("reference/y.link")).unwrap();
            symlink("x.txt", b.join("reference/y.link")).unwrap();
        },
        "content",
        "reference/y.link",
    );
    expect_red(
        "type",
        |b| {
            fs::remove_file(b.join("reference/y.link")).unwrap();
            write(&b.join("reference/y.link"), &[0, 1, 2, 3], 0o600);
        },
        "type",
        "reference/y.link",
    );
}

/// Writes FILE in A and in B with the two bodies, compares under ROWS and
/// returns (exit code, output).
fn modes_case(tag: &str, file: &str, a_body: &[u8], b_body: &[u8], rows: &str) -> (u8, String) {
    let (t, a, b) = pair(tag);
    write(&a.join(file), a_body, 0o644);
    write(&b.join(file), b_body, 0o644);
    let allow = allow_file(&t, rows);
    let (code, out, err) = compare(&manifest(&a, &[]), &manifest(&b, &[]), Some(&allow));
    assert_ne!(code, 2, "{tag}: {out}{err}");
    (code, out)
}

fn gz(body: &[u8], mtime: u32) -> Vec<u8> {
    let mut e = GzBuilder::new()
        .mtime(mtime)
        .write(Vec::new(), Compression::default());
    e.write_all(body).unwrap();
    e.finish().unwrap()
}

fn bgzf(body: &[u8]) -> Vec<u8> {
    let mut w = og_io::bgzf::Writer::new(Vec::new(), std::num::NonZero::new(1).unwrap()).unwrap();
    w.write_all(body).unwrap();
    w.finish().unwrap()
}

#[test]
fn gzip_mode_ignores_the_header_time_and_sees_the_content() {
    let row = "reference/*.gz\tgzip\ttest\tgzip -c writes the time into the header\t-\n";
    let (code, out) = modes_case(
        "gz-green",
        "reference/t.gz",
        &gz(b"a\tb\n", 1),
        &gz(b"a\tb\n", 2),
        row,
    );
    assert_eq!(code, 0, "{out}");
    assert!(
        out.contains("ALLOWED\tcontent\treference/t.gz\tgzip equal: "),
        "{out}"
    );
    let (code, out) = modes_case(
        "gz-red",
        "reference/t.gz",
        &gz(b"a\tb\n", 1),
        &gz(b"a\tc\n", 2),
        row,
    );
    assert_eq!(code, 1, "a content change inside the gzip is green:\n{out}");
    assert!(
        diff_lines(&out)[0].starts_with("DIFF\tcontent\treference/t.gz\tgzip "),
        "{out}"
    );
    let (code, _) = modes_case(
        "gz-norow",
        "reference/t.gz",
        &gz(b"a\tb\n", 1),
        &gz(b"a\tb\n", 2),
        "",
    );
    assert_eq!(code, 1, "without a row the header time is a difference");
    let (code, out) = modes_case(
        "gz-bad",
        "reference/t.gz",
        &gz(b"a\tb\n", 1),
        b"\x1f\x8bnot gzip",
        row,
    );
    assert_eq!(code, 1, "{out}");
}

fn vcf(date: &str, version: &str, alt: &str) -> Vec<u8> {
    format!(
        "##fileformat=VCFv4.2\n##fileDate={date}\n##contig=<ID=chr1,length=1000>\n##bcftools_viewVersion={version}\n##bcftools_viewCommand=view -Oz -o x.vcf.gz in.vcf.gz; Date=Wed Nov 26 1{date}\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t10\t.\tA\t{alt}\t.\tPASS\t.\nchr1\t20\t.\tC\tT\t.\tPASS\t.\n"
    )
    .into_bytes()
}

#[test]
fn vcf_mode_masks_only_the_dates() {
    let row = "clinvar/*.vcf.gz\tvcf\ttest\tbcftools writes the date of the run\t-\n";
    let a = bgzf(&vcf("20261126", "1.22+htslib-1.22", "G"));
    let b = bgzf(&vcf("20261127", "1.22+htslib-1.22", "G"));
    let (code, out) = modes_case("vcf-green", "clinvar/c.vcf.gz", &a, &b, row);
    assert_eq!(code, 0, "{out}");
    assert!(
        out.contains("ALLOWED\tcontent\tclinvar/c.vcf.gz\tvcf equal: "),
        "{out}"
    );
    let b = bgzf(&vcf("20261127", "1.22+htslib-1.22", "T"));
    let (code, out) = modes_case("vcf-record", "clinvar/c.vcf.gz", &a, &b, row);
    assert_eq!(code, 1, "a changed record is green:\n{out}");
    let b = bgzf(&vcf("20261126", "1.23+htslib-1.23", "G"));
    let (code, out) = modes_case("vcf-version", "clinvar/c.vcf.gz", &a, &b, row);
    assert_eq!(code, 1, "another bcftools version is green:\n{out}");
    // Plain VCF too.
    let (code, out) = modes_case(
        "vcf-plain",
        "clinvar/c.vcf",
        &vcf("20261126", "1.22", "G"),
        &vcf("20261127", "1.22", "G"),
        "clinvar/*.vcf\tvcf\ttest\tbcftools writes the date of the run\t-\n",
    );
    assert_eq!(code, 0, "{out}");
}

/// A BGZF VCF and its .tbi built by og-io, written to DIR/NAME(.tbi).
fn vcf_with_index(dir: &Path, name: &str, body: &[u8]) {
    let data = dir.join(name);
    write(&data, &bgzf(body), 0o644);
    let idx = og_io::tabix_write::build(fs::File::open(&data).unwrap()).unwrap();
    let tbi = dir.join(format!("{name}.tbi"));
    let _ = fs::remove_file(&tbi);
    idx.write_tbi(&tbi, true).unwrap();
    fs::set_permissions(&tbi, fs::Permissions::from_mode(0o644)).unwrap();
}

#[test]
fn index_mode_compares_names_and_checks_the_index_fits_its_data() {
    let rows = "clinvar/*.vcf.gz\tvcf\ttest\tbcftools writes the date of the run\t-\nclinvar/*.tbi\tindex\ttest\toffsets move with the header's length\t-\n";
    let index_case = |tag: &str, b_body: &[u8], break_b: bool| {
        let (t, a, b) = pair(tag);
        vcf_with_index(
            &a.join("clinvar"),
            "c.vcf.gz",
            &vcf("20261126", "1.22", "G"),
        );
        vcf_with_index(&b.join("clinvar"), "c.vcf.gz", b_body);
        if break_b {
            // B's index stays, its data file becomes one short block.
            write(
                &b.join("clinvar/c.vcf.gz"),
                &bgzf(b"##fileformat=VCFv4.2\n"),
                0o644,
            );
        }
        let allow = allow_file(&t, rows);
        let (code, out, err) = compare(&manifest(&a, &[]), &manifest(&b, &[]), Some(&allow));
        assert_ne!(code, 2, "{out}{err}");
        (code, out)
    };
    // A longer date line moves every offset: the .tbi bytes differ.
    let (code, out) = index_case("tbi-green", &vcf("2026112612345", "1.22", "G"), false);
    assert_eq!(code, 0, "{out}");
    assert!(
        out.contains("ALLOWED\tcontent\tclinvar/c.vcf.gz.tbi\tindex equal: "),
        "{out}"
    );
    let renamed = String::from_utf8(vcf("20261126", "1.22", "G"))
        .unwrap()
        .replace("chr1\t", "chr2\t");
    let (code, out) = index_case("tbi-names", renamed.as_bytes(), false);
    assert_eq!(code, 1, "{out}");
    assert!(
        out.contains("DIFF\tcontent\tclinvar/c.vcf.gz.tbi\tindex "),
        "an index of another sequence is green:\n{out}"
    );
    let (code, out) = index_case("tbi-unfit", &vcf("20261126", "1.22", "G"), true);
    assert_eq!(code, 1, "{out}");
    assert!(
        out.contains(
            "DIFF\tcontent\tclinvar/c.vcf.gz.tbi\tindex: the index does not fit its data file"
        ),
        "{out}"
    );
}

/// A CSI index by hand: one sequence, one bin, one chunk ending at `end`.
fn csi(end: u64, aux: &[u8]) -> Vec<u8> {
    let mut v = b"CSI\x01".to_vec();
    v.extend(14i32.to_le_bytes());
    v.extend(5i32.to_le_bytes());
    v.extend((aux.len() as i32).to_le_bytes());
    v.extend(aux);
    v.extend(1i32.to_le_bytes()); // n_ref
    v.extend(1i32.to_le_bytes()); // n_bin
    v.extend(4681u32.to_le_bytes());
    v.extend(0u64.to_le_bytes()); // loffset
    v.extend(1i32.to_le_bytes()); // n_chunk
    v.extend(0u64.to_le_bytes());
    v.extend(end.to_le_bytes());
    bgzf(&v)
}

#[test]
fn index_mode_reads_csi_too() {
    let rows = "x/*.csi\tindex\ttest\toffsets move with the data's blocks\t-\n";
    let case = |tag: &str, a_csi: Vec<u8>, b_csi: Vec<u8>| {
        let (t, a, b) = pair(tag);
        for (d, c) in [(&a, &a_csi), (&b, &b_csi)] {
            write(&d.join("x/d.bcf"), &bgzf(&[7u8; 300]), 0o644);
            write(&d.join("x/d.bcf.csi"), c, 0o644);
        }
        let allow = allow_file(&t, rows);
        let (code, out, err) = compare(&manifest(&a, &[]), &manifest(&b, &[]), Some(&allow));
        assert_ne!(code, 2, "{out}{err}");
        (code, out)
    };
    let (code, out) = case("csi-green", csi(10, b""), csi(20, b""));
    assert_eq!(code, 0, "{out}");
    let (code, out) = case("csi-aux", csi(10, b""), csi(20, b"chrX"));
    assert_eq!(code, 1, "another header is green:\n{out}");
    let (code, out) = case("csi-unfit", csi(10, b""), csi(1 << 40, b""));
    assert_eq!(code, 1, "an index past its data is green:\n{out}");
}

fn jf(host: &str, counts: &[u8]) -> Vec<u8> {
    let json = format!("{{\"alignment\":8,\"hostname\":\"{host}\",\"format\":\"binary/sorted\"}}");
    let mut v = format!("{:09}{json}", json.len()).into_bytes();
    while v.len() % 8 != 0 {
        v.push(0);
    }
    v.extend(counts);
    v
}

#[test]
fn jellyfish_mode_skips_the_header_and_compares_the_counts() {
    let row = "locityper/*.jf\tjellyfish\ttest\tthe header holds the host, folder and time of the count\t-\n";
    let (code, out) = modes_case(
        "jf-green",
        "locityper/c.jf",
        &jf("a1b2", &[1, 2, 3, 4]),
        &jf("host-two", &[1, 2, 3, 4]),
        row,
    );
    assert_eq!(code, 0, "{out}");
    let (code, out) = modes_case(
        "jf-red",
        "locityper/c.jf",
        &jf("a1b2", &[1, 2, 3, 4]),
        &jf("host-two", &[1, 2, 3, 5]),
        row,
    );
    assert_eq!(code, 1, "a changed count is green:\n{out}");
}

#[test]
fn numeric_mode_tolerates_float_noise_only() {
    let row = "ancient/*.evec\tnumeric\ttest\tsmartpca sums floats on several threads\t-\n";
    let a = b"#eigvals: 3.123456 2.5\nI1 0.0123456 -0.0456789 Pop\n";
    let noisy = b"#eigvals: 3.1234561 2.5000001\nI1 0.01234561 -0.04567891 Pop\n";
    let (code, out) = modes_case("num-green", "ancient/p.evec", a, noisy, row);
    assert_eq!(code, 0, "{out}");
    let moved = b"#eigvals: 3.123456 2.5\nI1 0.0124456 -0.0456789 Pop\n";
    let (code, out) = modes_case("num-red", "ancient/p.evec", a, moved, row);
    assert_eq!(code, 1, "a change in the fourth digit is green:\n{out}");
    let relabelled = b"#eigvals: 3.123456 2.5\nI1 0.0123456 -0.0456789 Other\n";
    let (code, out) = modes_case("num-word", "ancient/p.evec", a, relabelled, row);
    assert_eq!(code, 1, "a changed word is green:\n{out}");
}

#[test]
fn skip_and_perm_allow_only_what_they_name() {
    let (t, a, b) = pair("skip");
    let allow = allow_file(
        &t,
        "reference/sub/*.bin\tskip\ttest\twritten with the time of the install\tthe bytes are a packed time, no part of them is stable\ntools/*.sh\tperm\ttest\tmode tightened on purpose\t-\n",
    );
    write(&b.join("reference/sub/y.bin"), b"other bytes", 0o600);
    fs::set_permissions(b.join("tools/run.sh"), fs::Permissions::from_mode(0o700)).unwrap();
    let (ma, mb) = (manifest(&a, &[]), manifest(&b, &[]));
    let (code, out, err) = compare(&ma, &mb, Some(&allow));
    assert_eq!(code, 0, "{out}{err}");
    assert!(
        out.contains("ALLOWED\tcontent\treference/sub/y.bin\tskip: "),
        "{out}"
    );
    assert!(
        out.contains("ALLOWED\tmode\ttools/run.sh\t0755 -> 0700\tmode tightened on purpose"),
        "{out}"
    );
    // The same changes outside the rows are red; skip does not allow a mode.
    write(&b.join("reference/x.txt"), b"jello", 0o644);
    fs::set_permissions(
        b.join("reference/sub/y.bin"),
        fs::Permissions::from_mode(0o644),
    )
    .unwrap();
    let (code, out, _) = compare(&ma, &manifest(&b, &[]), Some(&allow));
    assert_eq!(code, 1, "{out}");
    assert_eq!(
        diff_lines(&out),
        [
            "DIFF\tmode\treference/sub/y.bin\t0600 -> 0644",
            "DIFF\tcontent\treference/x.txt\tsize 5 -> 5, sha256 2cf24dba5fb0 -> 187c9bceeb91",
        ],
        "{out}"
    );
    // skip still needs the file on both sides.
    fs::remove_file(b.join("reference/sub/y.bin")).unwrap();
    let (code, out, _) = compare(&ma, &manifest(&b, &[]), Some(&allow));
    assert_eq!(code, 1, "{out}");
    assert!(
        out.contains("DIFF\tremoved\treference/sub/y.bin\t"),
        "{out}"
    );
}

#[test]
fn stamps_are_compared_byte_for_byte() {
    let t = Tmp::new("stamp-rows");
    let a = t.0.join("a");
    fs::create_dir_all(&a).unwrap();
    tree(&a);
    let m = manifest(&a, &[]);
    for row in [
        "reference/DONE\tskip\ttest\ta date\tno other part\n",
        "reference/*\tgzip\ttest\tall gzip\t-\n",
        "**\tperm\ttest\tall\t-\n",
        "x/chr*.stamp\tskip\ttest\tchromosome stamps\tno other part\n",
        "x/INSTALLED\tperm\ttest\tmode\t-\n",
    ] {
        let allow = allow_file(&t, row);
        let (code, out, err) = compare(&m, &m, Some(&allow));
        assert_eq!(
            code, 2,
            "a row over a stamp was accepted: {row:?}\n{out}{err}"
        );
        assert!(err.contains("stamp"), "{err}");
    }
    // A row whose name is no stamp shape but names one: it does not apply.
    let (t, a, b) = pair("stamp-skip");
    write(&a.join("reference/foo.stamp"), b"v1", 0o644);
    write(&b.join("reference/foo.stamp"), b"v2", 0o644);
    let allow = allow_file(
        &t,
        "reference/foo.stamp\tskip\ttest\ta stamp\tnothing else\n",
    );
    let (code, out, _) = compare(&manifest(&a, &[]), &manifest(&b, &[]), Some(&allow));
    assert_eq!(code, 1, "{out}");
    assert!(
        out.contains("DIFF\tcontent\treference/foo.stamp\t"),
        "{out}"
    );
}

#[test]
fn bad_input_is_refused() {
    let t = Tmp::new("bad");
    let a = t.0.join("a");
    fs::create_dir_all(&a).unwrap();
    tree(&a);
    let m = manifest(&a, &[]);
    for (name, body) in [
        ("no-reason", format!("{HEAD}x\tgzip\ttest\t\t-\n")),
        ("no-component", format!("{HEAD}x\tgzip\t\tvolatile\t-\n")),
        ("bad-mode", format!("{HEAD}x\tsize\ttest\tvolatile\t-\n")),
        (
            "two-modes",
            format!("{HEAD}x\tgzip,vcf\ttest\tvolatile\t-\n"),
        ),
        (
            "skip-without-why",
            format!("{HEAD}x\tskip\ttest\tvolatile\t-\n"),
        ),
        (
            "why-without-skip",
            format!("{HEAD}x\tgzip\ttest\tvolatile\tnothing else\n"),
        ),
        (
            "four-columns",
            "path\tmode\tcomponent\treason\nx\tgzip\ttest\tvolatile\n".to_owned(),
        ),
        ("no-header", "x\tgzip\ttest\tvolatile\t-\n".to_owned()),
    ] {
        let allow = t.0.join(format!("{name}.tsv"));
        fs::write(&allow, body).unwrap();
        let (code, out, err) = compare(&m, &m, Some(&allow));
        assert_eq!(code, 2, "{name}: {out}{err}");
    }
    let junk = t.0.join("junk.tsv");
    fs::write(&junk, "a\tb\n").unwrap();
    assert_eq!(compare(&m, &junk, None).0, 2);
    assert_eq!(run(&[&s(&t.0.join("missing"))]).0, 2);
    assert_eq!(run(&[]).0, 2);
}

#[test]
fn odd_names_are_escaped_and_read_back() {
    let (_t, a, b) = pair("odd");
    write(&a.join("tools/tab\there"), b"x", 0o644);
    write(&b.join("tools/tab\there"), b"x", 0o644);
    write(&a.join("tools/line\nbreak"), b"y", 0o644);
    write(&b.join("tools/line\nbreak"), b"y", 0o644);
    let (ma, mb) = (manifest(&a, &[]), manifest(&b, &[]));
    let text = fs::read_to_string(&ma).unwrap();
    assert!(text.contains("tools/tab\\there\tf\t"), "{text}");
    assert!(text.contains("tools/line\\nbreak\tf\t"), "{text}");
    assert_eq!(compare(&ma, &mb, None).0, 0);
    // The control flips the file the escaped path names.
    let (code, out, err) = run(&["--flip", &s(&b), "tools/tab\\there"]);
    assert_eq!(code, 0, "{out}{err}");
    let (code, out, _) = compare(&ma, &manifest(&b, &[]), None);
    assert_eq!(code, 1);
    assert!(out.contains("DIFF\tcontent\ttools/tab\\there\t"), "{out}");
}

#[test]
fn the_control_picks_a_byte_compared_file_and_one_bit_is_red() {
    let (t, a, b) = pair("control");
    write(&a.join("reference/t.gz"), &gz(b"z", 1), 0o644);
    write(&b.join("reference/t.gz"), &gz(b"z", 1), 0o644);
    let allow = allow_file(
        &t,
        "reference/sub/*.bin\tskip\ttest\tvolatile\tnothing stable in it\n",
    );
    let mb = manifest(&b, &[]);
    let (code, picked, _) = run(&["--control", &s(&mb), "--allow", &s(&allow)]);
    assert_eq!(code, 0);
    // y.bin (4 bytes) is smaller but skipped; x.txt (5 bytes) is next.
    assert_eq!(picked, "reference/x.txt\n");
    assert_eq!(run(&["--flip", &s(&b), "reference/x.txt"]).0, 0);
    let (code, out, _) = compare(&manifest(&a, &[]), &manifest(&b, &[]), Some(&allow));
    assert_eq!(code, 1, "{out}");
    assert_eq!(diff_lines(&out).len(), 1);
    assert!(
        diff_lines(&out)[0].starts_with("DIFF\tcontent\treference/x.txt\t"),
        "{out}"
    );
    assert_eq!(run(&["--flip", &s(&b), "reference/x.txt"]).0, 0);
    assert_eq!(fs::read(b.join("reference/x.txt")).unwrap(), b"hello");
    let empty = t.0.join("empty");
    fs::create_dir_all(&empty).unwrap();
    assert_eq!(
        run(&["--control", &s(&manifest(&empty, &[]))]).0,
        1,
        "no file to plant in"
    );
}

#[test]
fn the_log_compare_masks_times_and_sees_the_missing_list() {
    let t = Tmp::new("logs");
    let log = |name: &str, body: &str| {
        let p = t.0.join(name);
        fs::write(&p, body).unwrap();
        s(&p)
    };
    let base = "[OK] Docker is running.\nDownloading x 2026-11-26T10:00:00Z\r100  10  100  10    0     0   50      0 --:--:-- --:--:-- --:--:--  50\ndone in 3 s\n  Setup finished: core ready, 1 item(s) missing\nMissing (1):\n  - the KIR data\n      retry: ./scripts/setup.sh --kir-data /g\n";
    let a = log("a.log", base);
    let b = log(
        "b.log",
        &base
            .replace("10:00:00Z", "11:12:13Z")
            .replace("in 3 s", "in 41 s"),
    );
    let (code, out, err) = run(&["--compare-logs", &a, &b]);
    assert_eq!(code, 0, "{out}{err}");
    let c = log("c.log", &base.replace("--kir-data", "--viral-data"));
    let (code, out, _) = run(&["--compare-logs", &a, &c]);
    assert_eq!(code, 1, "{out}");
    assert!(out.contains("DIFF\tmissing\t"), "{out}");
    let d = log("d.log", &format!("{base}[WARN] one more line\n"));
    let (code, out, _) = run(&["--compare-logs", &a, &d]);
    assert_eq!(code, 1, "{out}");
    assert!(
        out.contains("DIFF\tlog\tB only\t[WARN] one more line"),
        "{out}"
    );
}

#[test]
fn the_committed_allowlist_reads_cleanly() {
    let t = Tmp::new("committed");
    let a = t.0.join("a");
    fs::create_dir_all(&a).unwrap();
    tree(&a);
    let m = manifest(&a, &[]);
    let allow = repo().join("tests/setup/install-manifest-allow.tsv");
    let (code, out, err) = compare(&m, &m, Some(&allow));
    assert_eq!(code, 0, "{out}{err}");
}
