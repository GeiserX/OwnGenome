//! `og dev install-manifest`: a manifest of an install tree, and the
//! compare of two manifests that the Rust port of setup must pass
//! (`tests/setup/ab.sh`).
//!
//! 1. A manifest lists every file, folder and symlink with its type, size,
//!    sha256, mode and the digests its compare modes read, sorted, with no
//!    trace of where the tree sits; the `downloads` folder at the top and
//!    `--exclude` paths are left out, and a `*.log` file is listed for its
//!    presence only.
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
//! 7. tests/setup/ab.sh, with a stub setup.sh: a side A that fails, a run
//!    that installs nothing, a control that would land in the seed and a
//!    file another user owns each make it red or stop it.

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
        "reference/sub/build.log\tf\t-\t-\t-\t-",
        "reference/sub/y.bin\tf\t4\t054edec1d0211f624fed0cbca9d4f9400b0e491c43742af2c5b0abebf0c990d8\t0600\t-",
        "reference/x.txt\tf\t5\t2cf24dba5fb0a30e26e83b2ac5b9e29e1b161e5c1fa7425e73043362938b9824\t0644\t-",
        "reference/y.link\tl\t-\t->sub/y.bin\t-\t-",
        "tools\td\t-\t-\t0755\t-",
    ];
    assert_eq!(rows.len(), want.len() + 2, "{text}");
    assert_eq!(&rows[1..9], &want, "{text}");
    assert!(
        rows[9].starts_with("tools/run.sh\tf\t10\t") && rows[9].ends_with("\t0755\t-"),
        "{text}"
    );
    assert!(
        !text.contains("downloads"),
        "the downloads folder is listed:\n{text}"
    );
    assert!(
        !text.contains(&s(&a)),
        "the manifest names where the tree sits:\n{text}"
    );
    let ex = fs::read_to_string(manifest(&a, &["--exclude", "reference/sub"])).unwrap();
    assert!(!ex.contains("reference/sub"), "{ex}");
    assert!(ex.contains("reference/x.txt"), "{ex}");
}

#[test]
fn a_log_counts_for_its_presence_only() {
    // Other bytes in a log are green; a log in one tree only is red.
    let (_t, a, b) = pair("log");
    write(
        &b.join("reference/sub/build.log"),
        b"built at 13:41 on host-two",
        0o600,
    );
    let (ma, mb) = (manifest(&a, &[]), manifest(&b, &[]));
    let (code, out, err) = compare(&ma, &mb, None);
    assert_eq!(code, 0, "{out}{err}");
    write(&b.join("reference/sub/results.log"), b"x", 0o644);
    let (code, out, _) = compare(&ma, &manifest(&b, &[]), None);
    assert_eq!(code, 1, "a log in B only is green:\n{out}");
    assert_eq!(
        diff_lines(&out),
        ["DIFF\tadded\treference/sub/results.log\tf in B only"],
        "{out}"
    );
    fs::remove_file(b.join("reference/sub/build.log")).unwrap();
    let (code, out, _) = compare(&ma, &manifest(&b, &[]), None);
    assert_eq!(code, 1, "{out}");
    assert!(
        out.contains("DIFF\tremoved\treference/sub/build.log\t"),
        "a log in A only is green:\n{out}"
    );
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
    // Plain gzip where BGZF was: the same content, but tabix cannot read it.
    let (code, out) = modes_case(
        "gz-kind",
        "reference/t.gz",
        &bgzf(b"a\tb\n"),
        &gz(b"a\tb\n", 0),
        row,
    );
    assert_eq!(code, 1, "plain gzip in place of BGZF is green:\n{out}");
    assert!(
        diff_lines(&out)[0].starts_with("DIFF\tcontent\treference/t.gz\tgzip bgzf:"),
        "{out}"
    );
    let (code, out) = modes_case(
        "gz-both-bgzf",
        "reference/t.gz",
        &bgzf(b"a\tb\n"),
        &bgzf(b"a\tb\n"),
        row,
    );
    assert_eq!(code, 0, "{out}");
    // The data block's BSIZE spans the end block too: the headers chain to
    // the end of the file, the members do not, and tabix cannot read it.
    let mut forged = bgzf(b"a\tb\n");
    let first = u16::from_le_bytes([forged[16], forged[17]]) as usize + 1;
    let total = (forged.len() - 1) as u16;
    assert!(first < forged.len());
    forged[16..18].copy_from_slice(&total.to_le_bytes());
    let (code, out) = modes_case(
        "gz-forged",
        "reference/t.gz",
        &bgzf(b"a\tb\n"),
        &forged,
        row,
    );
    assert_eq!(code, 1, "a BSIZE past its member is BGZF:\n{out}");
    assert!(
        diff_lines(&out)[0].starts_with("DIFF\tcontent\treference/t.gz\tgzip bgzf:"),
        "{out}"
    );
    // A block whose header, CRC and ISIZE agree but that decodes to more
    // than 64 KiB is not BGZF.
    let body = vec![b'a'; 70_000];
    let (code, out) = modes_case(
        "gz-big-member",
        "reference/t.gz",
        &bgzf(&body),
        &big_member(&body),
        row,
    );
    assert_eq!(code, 1, "a 70000-byte block is BGZF:\n{out}");
    // Two plain gzip members in one buffer: plain gzip, read to the end.
    let two = [gz(b"a\t", 0), gz(b"b\n", 0)].concat();
    let (code, out) = modes_case(
        "gz-two-members",
        "reference/t.gz",
        &two,
        &gz(b"a\tb\n", 0),
        row,
    );
    assert_eq!(code, 0, "{out}");
}

/// One gzip member with BGZF's header (its BSIZE right) holding all of
/// `body`, then BGZF's end block.
fn big_member(body: &[u8]) -> Vec<u8> {
    let mut e = GzBuilder::new()
        .extra(b"BC\x02\x00\x00\x00".to_vec())
        .write(Vec::new(), Compression::default());
    e.write_all(body).unwrap();
    let mut v = e.finish().unwrap();
    let bsize = (v.len() - 1) as u16;
    v[16..18].copy_from_slice(&bsize.to_le_bytes());
    let end = bgzf(b"");
    v.extend_from_slice(&end[end.len() - 28..]);
    v
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
    // The same records as plain gzip, or not compressed at all, where a
    // BGZF VCF was.
    for (tag, other) in [
        ("vcf-gzip", gz(&vcf("20261127", "1.22+htslib-1.22", "G"), 0)),
        ("vcf-uncompressed", vcf("20261127", "1.22+htslib-1.22", "G")),
    ] {
        let (code, out) = modes_case(tag, "clinvar/c.vcf.gz", &a, &other, row);
        assert_eq!(code, 1, "{tag}: another container is green:\n{out}");
    }
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

#[test]
fn an_index_of_fewer_records_is_red_beside_the_full_data() {
    let rows = "clinvar/*.tbi\tindex\ttest\toffsets move with the header's length\t-\n";
    // The short data fills exactly one BGZF block (65,280 bytes), so the
    // full data starts with the same block and every offset of the short
    // index still fits it: only what the index covers tells them apart.
    let base = String::from_utf8(vcf("20261126", "1.22", "G")).unwrap();
    let short = base.replacen(
        "##fileDate=20261126\n",
        &format!("##fileDate=20261126{}\n", "9".repeat(0xff00 - base.len())),
        1,
    );
    assert_eq!(short.len(), 0xff00);
    let full =
        short.clone() + "chr1\t5000000\t.\tG\tA\t.\tPASS\t.\nchr1\t5000100\t.\tT\tC\t.\tPASS\t.\n";
    let (t, a, b) = pair("tbi-stale");
    vcf_with_index(&a.join("clinvar"), "c.vcf.gz", full.as_bytes());
    // B: the index of the data less its last two records, beside the full
    // data.
    vcf_with_index(&b.join("clinvar"), "c.vcf.gz", short.as_bytes());
    write(&b.join("clinvar/c.vcf.gz"), &bgzf(full.as_bytes()), 0o644);
    let allow = allow_file(&t, rows);
    let (mb, ma) = (manifest(&b, &[]), manifest(&a, &[]));
    let text = fs::read_to_string(&mb).unwrap();
    assert!(
        !text.contains("index=invalid"),
        "the short index must fit the full data:\n{text}"
    );
    let (code, out, err) = compare(&ma, &mb, Some(&allow));
    assert_eq!(code, 1, "an index of fewer records is green:\n{out}{err}");
    assert!(
        out.contains("DIFF\tcontent\tclinvar/c.vcf.gz.tbi\tindex "),
        "{out}"
    );
}

/// Random letters, the same for one seed.
fn noise(seed: &mut u64, n: usize) -> String {
    const A: &[u8] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    (0..n)
        .map(|_| {
            *seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            A[(*seed >> 58) as usize] as char
        })
        .collect()
}

/// The bin numbers of each sequence of a .tbi, the pseudo-bin left out.
fn tbi_bins(path: &Path) -> Vec<Vec<u32>> {
    use std::io::Read;
    let mut p = Vec::new();
    flate2::read::MultiGzDecoder::new(fs::File::open(path).unwrap())
        .read_to_end(&mut p)
        .unwrap();
    let mut at = 4;
    let n_ref = int_at(&p, &mut at);
    at += 24;
    let l_nm = int_at(&p, &mut at);
    at += l_nm;
    let mut out = Vec::new();
    for _ in 0..n_ref {
        let mut bins = Vec::new();
        let n_bin = int_at(&p, &mut at);
        for _ in 0..n_bin {
            let bin = int_at(&p, &mut at) as u32;
            let n_chunk = int_at(&p, &mut at);
            at += 16 * n_chunk;
            if bin != 37450 {
                bins.push(bin);
            }
        }
        let n_intv = int_at(&p, &mut at);
        at += 8 * n_intv;
        bins.sort_unstable();
        out.push(bins);
    }
    out
}

fn int_at(p: &[u8], at: &mut usize) -> usize {
    let v = i32::from_le_bytes(p[*at..*at + 4].try_into().unwrap());
    *at += 4;
    v as usize
}

#[test]
fn bin_numbers_that_move_with_the_compressed_sizes_stay_out_of_the_index_digest() {
    // htslib (and og-io, which writes what tabix writes) merges a bin into
    // its parent when its chunks span less than 64 KiB of the compressed
    // file. 126 records of random text in one 16 kb window, about 130 kB:
    // after a short header they end in the second block (a span of one
    // block, merged); after a header 1,900 bytes longer they end in the
    // third (two blocks, kept). The record spanning the window's end makes
    // the parent bin exist.
    let body = |pad: usize| {
        let mut s = format!(
            "##fileformat=VCFv4.2\n##fileDate={}\n##contig=<ID=chr1,length=1000000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
            "9".repeat(pad)
        );
        let mut seed = 7;
        for i in 0..126 {
            s += &format!(
                "chr1\t{}\t.\tA\tG\t.\tPASS\tN={}\n",
                163_900 + i * 100,
                noise(&mut seed, 1000)
            );
        }
        s += &format!("chr1\t180200\t.\t{}\tA\t.\tPASS\t.\n", "A".repeat(50));
        s.into_bytes()
    };
    let rows = "clinvar/*.vcf.gz\tvcf\ttest\tthe date of the run\t-\nclinvar/*.tbi\tindex\ttest\toffsets move with the header's length\t-\n";
    let (t, a, b) = pair("tbi-merge");
    vcf_with_index(&a.join("clinvar"), "c.vcf.gz", &body(10));
    vcf_with_index(&b.join("clinvar"), "c.vcf.gz", &body(1910));
    let (ba, bb) = (
        tbi_bins(&a.join("clinvar/c.vcf.gz.tbi")),
        tbi_bins(&b.join("clinvar/c.vcf.gz.tbi")),
    );
    assert_ne!(ba, bb, "the case no longer moves a bin: {ba:?}");
    let allow = allow_file(&t, rows);
    let (code, out, err) = compare(&manifest(&a, &[]), &manifest(&b, &[]), Some(&allow));
    assert_eq!(
        code, 0,
        "the same records indexed after another header are red:\n{out}{err}"
    );
    assert!(
        out.contains("ALLOWED\tcontent\tclinvar/c.vcf.gz.tbi\tindex equal: "),
        "{out}"
    );
}

#[test]
fn the_two_lookup_fixture_indexes_get_their_own_digests() {
    let t = Tmp::new("fixtures");
    let d = t.0.join("d");
    let src = repo().join("tests/parity/inputs/lookup");
    for f in ["ann.vcf.gz", "ann.vcf.gz.tbi", "g.vcf.gz", "g.vcf.gz.tbi"] {
        write(&d.join(f), &fs::read(src.join(f)).unwrap(), 0o644);
    }
    let text = fs::read_to_string(manifest(&d, &[])).unwrap();
    let index = |name: &str| {
        text.lines()
            .find(|l| l.starts_with(&format!("{name}\t")))
            .and_then(|l| l.split("index=").nth(1))
            .unwrap_or_else(|| panic!("no index digest for {name}:\n{text}"))
            .to_owned()
    };
    let (ann, g) = (index("ann.vcf.gz.tbi"), index("g.vcf.gz.tbi"));
    assert!(ann != "invalid" && g != "invalid", "{text}");
    assert_ne!(ann, g, "two indexes of other files get one digest");
}

/// A CSI index by hand: one sequence, one bin, one chunk ending at `end`.
fn csi(end: u64, aux: &[u8]) -> Vec<u8> {
    bgzf(&csi_plain(end, aux))
}

/// The decompressed bytes of [`csi`].
fn csi_plain(end: u64, aux: &[u8]) -> Vec<u8> {
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
    v
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
    // The data is one block of 300 decoded bytes, then the BGZF end block.
    let data = bgzf(&[7u8; 300]);
    let eof = (data.len() as u64) << 16;
    let (code, out) = case("csi-block-end", csi(10, b""), csi(300, b""));
    assert_eq!(code, 0, "an end at its block's last byte is red:\n{out}");
    let (code, out) = case("csi-past-block", csi(10, b""), csi(301, b""));
    assert_eq!(code, 1, "an end past its block's bytes is green:\n{out}");
    let (code, out) = case("csi-eof", csi(10, b""), csi(eof, b""));
    assert_eq!(code, 0, "an end at the end of the file is red:\n{out}");
    let (code, out) = case("csi-past-eof", csi(10, b""), csi(eof | 1, b""));
    assert_eq!(code, 1, "an end inside no block is green:\n{out}");
    // An index over 64 MiB is not read: it is invalid, so it is red.
    let mut big = csi_plain(10, b"");
    big.resize(big.len() + (64 << 20), 0);
    let (code, out) = case("csi-big", csi(10, b""), bgzf(&big));
    assert_eq!(code, 1, "an index over 64 MiB is green:\n{out}");
    // The pseudo-bin's record counts and the count of records with no
    // position are compared; its offsets are not.
    let (code, out) = case("csi-meta", csi_meta(10, 5, None), csi_meta(20, 5, None));
    assert_eq!(code, 0, "{out}");
    let (code, out) = case("csi-mapped", csi_meta(10, 5, None), csi_meta(20, 4, None));
    assert_eq!(code, 1, "another record count is green:\n{out}");
    let (code, out) = case(
        "csi-no-coor",
        csi_meta(10, 5, Some(0)),
        csi_meta(10, 5, Some(3)),
    );
    assert_eq!(
        code, 1,
        "another count of unplaced records is green:\n{out}"
    );
}

/// A CSI like [`csi`] plus htslib's pseudo-bin (37450 at depth 5: its
/// offsets, then `mapped` records and none unmapped) and, when given, the
/// count of records with no position.
fn csi_meta(end: u64, mapped: u64, no_coor: Option<u64>) -> Vec<u8> {
    let mut v = b"CSI\x01".to_vec();
    v.extend(14i32.to_le_bytes());
    v.extend(5i32.to_le_bytes());
    v.extend(0i32.to_le_bytes());
    v.extend(1i32.to_le_bytes()); // n_ref
    v.extend(2i32.to_le_bytes()); // n_bin
    v.extend(4681u32.to_le_bytes());
    v.extend(0u64.to_le_bytes());
    v.extend(1i32.to_le_bytes());
    v.extend(0u64.to_le_bytes());
    v.extend(end.to_le_bytes());
    v.extend(37450u32.to_le_bytes());
    v.extend(0u64.to_le_bytes());
    v.extend(2i32.to_le_bytes());
    v.extend(0u64.to_le_bytes());
    v.extend(end.to_le_bytes());
    v.extend(mapped.to_le_bytes());
    v.extend(0u64.to_le_bytes());
    if let Some(n) = no_coor {
        v.extend(n.to_le_bytes());
    }
    bgzf(&v)
}

#[test]
fn an_index_beside_a_fifo_or_a_block_with_another_extra_field_is_invalid() {
    let digest = |d: &Path| {
        let text = fs::read_to_string(manifest(d, &[])).unwrap();
        text.lines()
            .find(|l| l.starts_with("x/d.bcf.csi\t"))
            .and_then(|l| l.split("index=").nth(1))
            .unwrap_or_else(|| panic!("no index digest:\n{text}"))
            .to_owned()
    };
    let t = Tmp::new("csi-odd");
    let good = t.0.join("good");
    write(&good.join("x/d.bcf"), &bgzf(&[7u8; 300]), 0o644);
    write(&good.join("x/d.bcf.csi"), &csi(10, b""), 0o644);
    assert_ne!(digest(&good), "invalid");
    // XLEN 12: an extra field htslib does not read as BGZF, so the trailer
    // is not where BSIZE says.
    let mut data = bgzf(&[7u8; 300]);
    data[10] = 12;
    let odd = t.0.join("xlen");
    write(&odd.join("x/d.bcf"), &data, 0o644);
    write(&odd.join("x/d.bcf.csi"), &csi(10, b""), 0o644);
    assert_eq!(digest(&odd), "invalid", "a block with XLEN 12 fits");
    // A FIFO where the data file is: invalid at once, never a blocked read.
    let fifo = t.0.join("fifo");
    write(&fifo.join("x/d.bcf.csi"), &csi(10, b""), 0o644);
    let st = Command::new("mkfifo")
        .arg(fifo.join("x/d.bcf"))
        .status()
        .unwrap();
    assert!(st.success());
    let (tx, rx) = std::sync::mpsc::channel();
    let f2 = fifo.clone();
    std::thread::spawn(move || {
        let _ = tx.send(digest(&f2));
    });
    let got = rx
        .recv_timeout(std::time::Duration::from_secs(20))
        .expect("the manifest blocks on a FIFO beside an index");
    assert_eq!(got, "invalid");
}

#[test]
fn an_index_beside_a_damaged_block_or_a_link_or_with_stray_bytes_is_invalid() {
    let digest = |d: &Path| {
        let text = fs::read_to_string(manifest(d, &[])).unwrap();
        text.lines()
            .find(|l| l.starts_with("x/d.bcf.csi\t"))
            .and_then(|l| l.split("index=").nth(1))
            .unwrap_or_else(|| panic!("no index digest:\n{text}"))
            .to_owned()
    };
    let t = Tmp::new("csi-damaged");
    let case = |name: &str, data: &[u8], index: &[u8]| {
        let d = t.0.join(name);
        write(&d.join("x/d.bcf"), data, 0o644);
        write(&d.join("x/d.bcf.csi"), index, 0o644);
        digest(&d)
    };
    let good = bgzf(&[7u8; 300]);
    let end = u16::from_le_bytes([good[16], good[17]]) as usize + 1;
    assert_ne!(case("good", &good, &csi(300, b"")), "invalid");
    // The trailer claims 301 bytes where the block holds 300: offset 301
    // reads past the data.
    let mut isize = good.clone();
    isize[end - 4..end].copy_from_slice(&301u32.to_le_bytes());
    assert_eq!(
        case("isize", &isize, &csi(301, b"")),
        "invalid",
        "a block whose ISIZE lies fits"
    );
    // A CRC that does not match the block's bytes.
    let mut crc = good.clone();
    crc[end - 8] ^= 1;
    assert_eq!(
        case("crc", &crc, &csi(10, b"")),
        "invalid",
        "a block with a bad CRC fits"
    );
    // A block that decodes to more than 64 KiB, its CRC and ISIZE right.
    assert_eq!(
        case("big", &big_member(&[7u8; 70_000]), &csi(10, b"")),
        "invalid",
        "a 70000-byte block fits"
    );
    // Four stray bytes after the last sequence: neither nothing nor the
    // eight-byte count of records with no position.
    let mut stray = csi_plain(10, b"");
    stray.extend([0u8; 4]);
    assert_eq!(
        case("stray", &good, &bgzf(&stray)),
        "invalid",
        "an index with stray bytes at its end is read"
    );
    // A data file that is a link: the manifest compares the link's text,
    // not the target's bytes, so the index cannot vouch for them.
    let link = t.0.join("link");
    write(&link.join("x/d.bcf.csi"), &csi(10, b""), 0o644);
    std::os::unix::fs::symlink(t.0.join("good/x/d.bcf"), link.join("x/d.bcf")).unwrap();
    assert_eq!(digest(&link), "invalid", "an index beside a link fits");
}

fn jf(host: &str, counts: &[u8]) -> Vec<u8> {
    jf_with(host, 2, counts)
}

/// A jellyfish file laid out as jellyfish 2.3.1 writes it: its host, folder
/// and time come from `host`, its `counter_len` from `counter`. The header
/// is padded with NULs to its alignment, and the 9-digit length counts
/// them (a real count in the pinned Locityper image wrote `000000407`
/// with 2 NULs inside, and with a host 36 bytes longer `000000447` with 6).
fn jf_with(host: &str, counter: u32, counts: &[u8]) -> Vec<u8> {
    let mut span = format!(
        "{{\"alignment\":8,\"counter_len\":{counter},\"format\":\"binary/sorted\",\"hostname\":\"{host}\",\"pwd\":\"/{host}\",\"time\":\"Thu Nov 26 03:56:51 {host}\"}}"
    )
    .into_bytes();
    while (9 + span.len()) % 8 != 0 {
        span.push(0);
    }
    let mut v = format!("{:09}", span.len()).into_bytes();
    v.extend(span);
    v.extend(counts);
    v
}

/// The NULs inside the header span of `jf_with(host, ..)`.
fn jf_pad(host: &str) -> usize {
    let v = jf_with(host, 2, &[]);
    v.iter().rev().take_while(|&&b| b == 0).count()
}

#[test]
fn jellyfish_mode_skips_the_header_and_compares_the_counts() {
    let row = "locityper/*.jf\tjellyfish\ttest\tthe header holds the host, folder and time of the count\t-\n";
    // The two hosts need another number of padding NULs.
    assert_ne!(jf_pad("a1b2"), jf_pad("host-two"));
    let (code, out) = modes_case(
        "jf-green",
        "locityper/c.jf",
        &jf("a1b2", &[1, 2, 3, 4]),
        &jf("host-two", &[1, 2, 3, 4]),
        row,
    );
    assert_eq!(
        code, 0,
        "two hosts that pad the header differently are red:\n{out}"
    );
    let (code, out) = modes_case(
        "jf-red",
        "locityper/c.jf",
        &jf("a1b2", &[1, 2, 3, 4]),
        &jf("host-two", &[1, 2, 3, 5]),
        row,
    );
    assert_eq!(code, 1, "a changed count is green:\n{out}");
    // The header's layout counts: the same bytes read with another counter
    // length are other counts.
    let (code, out) = modes_case(
        "jf-layout",
        "locityper/c.jf",
        &jf_with("a1b2", 2, &[1, 2, 3, 4]),
        &jf_with("a1b2", 4, &[1, 2, 3, 4]),
        row,
    );
    assert_eq!(code, 1, "a changed counter_len is green:\n{out}");
    let (code, out) = modes_case(
        "jf-empty",
        "locityper/c.jf",
        &jf_with("a1b2", 2, &[]),
        &jf_with("a1b2", 4, &[]),
        row,
    );
    assert_eq!(code, 1, "a changed layout with no counts is green:\n{out}");
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
        "reference/foo.stamp\tskip\ttest\ta stamp\tnothing else\n",
        "reference/foo*\tgzip\ttest\tall foo\t-\n",
        "reference/foo.*\tgzip\ttest\tall foo\t-\n",
        "x/*.st?mp\tgzip\ttest\tstamps\t-\n",
        "x/gnomad-sv.inst*\tgzip\ttest\tsv\t-\n",
        "x/?ONE\tperm\ttest\tmode\t-\n",
        "x/version.tx?\tperm\ttest\tmode\t-\n",
        "fine/pins.tsv\tperm\ttest\tmode\t-\n",
        "fine/*.tsv\tnumeric\ttest\tfloats\t-\n",
        "locityper/*.sha256\tskip\ttest\ta checksum\tnothing else\n",
        "**/*.sha2*\tperm\ttest\tmode\t-\n",
    ] {
        let allow = allow_file(&t, row);
        let (code, out, err) = compare(&m, &m, Some(&allow));
        assert_eq!(
            code, 2,
            "a row over a stamp was accepted: {row:?}\n{out}{err}"
        );
        assert!(err.contains("stamp"), "{err}");
    }
    // Names close to a stamp's that no stamp can have are rows like any.
    for row in [
        "x/*.stamps\tgzip\ttest\tnot a stamp\t-\n",
        "x/DONE.txt\tperm\ttest\tnot a stamp\t-\n",
        "x/*.release.gz\tgzip\ttest\tnot a stamp\t-\n",
        "x/?????\tperm\ttest\tfive bytes, no stamp has five\t-\n",
    ] {
        let allow = allow_file(&t, row);
        let (code, out, err) = compare(&m, &m, Some(&allow));
        assert_eq!(
            code, 0,
            "a row over no stamp was refused: {row:?}\n{out}{err}"
        );
    }
    // A stamp's bytes are compared even so.
    let (_t, a, b) = pair("stamp-bytes");
    write(&a.join("reference/foo.stamp"), b"v1", 0o644);
    write(&b.join("reference/foo.stamp"), b"v2", 0o644);
    let (code, out, _) = compare(&manifest(&a, &[]), &manifest(&b, &[]), None);
    assert_eq!(code, 1, "{out}");
    assert!(
        out.contains("DIFF\tcontent\treference/foo.stamp\t"),
        "{out}"
    );
}

#[test]
fn a_file_named_downloads_and_a_name_starting_with_hash_are_listed() {
    let (_t, a, b) = pair("names");
    fs::remove_dir_all(b.join("downloads")).unwrap();
    write(&b.join("downloads"), b"not a folder", 0o644);
    write(&b.join("#build-id"), b"1", 0o644);
    write(&b.join("reference/#x"), b"2", 0o644);
    let (ma, mb) = (manifest(&a, &[]), manifest(&b, &[]));
    let (code, out, err) = compare(&ma, &mb, None);
    assert_eq!(code, 1, "{out}{err}");
    for p in ["downloads", "#build-id", "reference/#x"] {
        assert!(out.contains(&format!("DIFF\tadded\t{p}\t")), "{p}: {out}");
    }
    assert_eq!(compare(&mb, &mb, None).0, 0);
    let (code, out, _) = compare(&mb, &ma, None);
    assert_eq!(code, 1, "{out}");
    assert!(out.contains("DIFF\tremoved\t#build-id\t"), "{out}");
}

#[test]
fn flip_stays_inside_the_tree_and_writes_only_a_file() {
    let t = Tmp::new("flip-out");
    let d = t.0.join("d");
    let outside = t.0.join("o.txt");
    write(&d.join("in/x.txt"), b"inside", 0o644);
    write(&outside, b"outside", 0o644);
    symlink(&outside, d.join("link")).unwrap();
    symlink(&t.0, d.join("up")).unwrap();
    for path in [
        "../o.txt".to_owned(),
        s(&outside),
        "link".to_owned(),
        "up/o.txt".to_owned(),
        "in/../../o.txt".to_owned(),
        "./in/x.txt".to_owned(),
        "in".to_owned(),
    ] {
        let (code, out, err) = run(&["--flip", &s(&d), &path]);
        assert_eq!(code, 2, "--flip {path} was not refused: {out}{err}");
    }
    assert_eq!(fs::read(&outside).unwrap(), b"outside");
    assert_eq!(fs::read(d.join("in/x.txt")).unwrap(), b"inside");
    assert_eq!(run(&["--flip", &s(&d), "in/x.txt"]).0, 0);
    assert_ne!(fs::read(d.join("in/x.txt")).unwrap(), b"inside");
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
            "empty-why-without-skip",
            format!("{HEAD}x\tgzip\ttest\tvolatile\t\n"),
        ),
        (
            "blank-why-without-skip",
            format!("{HEAD}x\tgzip\ttest\tvolatile\t  \n"),
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
    // A row with one digest kind twice: the second must not stand in for
    // the first.
    let (aa, bb) = ("a".repeat(64), "b".repeat(64));
    let twice = t.0.join("twice.tsv");
    fs::write(
        &twice,
        format!("#og-install-manifest 1\npath\ttype\tsize\tsha256\tmode\tdigests\nx.gz\tf\t3\tabc\t0644\tgzip=gzip:{aa},gzip=gzip:{bb}\n"),
    )
    .unwrap();
    let (code, out, err) = compare(&twice, &twice, None);
    assert_eq!(code, 2, "a row with two gzip digests is read:\n{out}{err}");
    // A digest the manifest never writes: two empty values must not compare
    // equal under a gzip row.
    for (name, digests) in [
        ("empty", "gzip=".to_owned()),
        ("no-container", format!("gzip={aa}")),
        ("short", "index=abc".to_owned()),
        ("vcf-bad", "vcf=bad".to_owned()),
        ("unknown", format!("size={aa}")),
    ] {
        let p = t.0.join(format!("digest-{name}.tsv"));
        fs::write(
            &p,
            format!("#og-install-manifest 1\npath\ttype\tsize\tsha256\tmode\tdigests\nx.gz\tf\t3\tabc\t0644\t{digests}\n"),
        )
        .unwrap();
        let (code, out, err) = compare(&p, &p, None);
        assert_eq!(code, 2, "{name}: {digests} is read:\n{out}{err}");
    }
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
fn the_control_leaves_out_what_the_seed_put_there() {
    let t = Tmp::new("control-seed");
    let seed = t.0.join("seed");
    write(&seed.join("ref.fa"), b"ACGT\n", 0o644);
    let b = t.0.join("b");
    copy(&seed, &b);
    write(&b.join("kir/kir.dat"), b"the component's own file\n", 0o644);
    write(&b.join("kir/run.log"), b"x", 0o644);
    let (ms, mb) = (manifest(&seed, &[]), manifest(&b, &[]));
    let (code, picked, _) = run(&["--control", &s(&mb)]);
    assert_eq!((code, picked.as_str()), (0, "ref.fa\n"));
    let (code, picked, err) = run(&["--control", &s(&mb), "--not-in", &s(&ms)]);
    assert_eq!(
        (code, picked.as_str()),
        (0, "kir/kir.dat\n"),
        "the control picks a seed file: {err}"
    );
    // A seed file the component wrote again with other bytes is its own.
    write(&b.join("ref.fa"), b"ACGTT\n", 0o644);
    fs::remove_file(b.join("kir/kir.dat")).unwrap();
    let (code, picked, _) = run(&["--control", &s(&manifest(&b, &[])), "--not-in", &s(&ms)]);
    assert_eq!((code, picked.as_str()), (0, "ref.fa\n"));
    // Nothing but the seed, and a log (presence only): no file to plant in.
    fs::remove_file(b.join("ref.fa")).unwrap();
    copy(&seed.join("ref.fa"), &b.join("ref.fa"));
    let (code, picked, _) = run(&["--control", &s(&manifest(&b, &[])), "--not-in", &s(&ms)]);
    assert_eq!(code, 1, "picked {picked}");
    assert_eq!(run(&["--control", &s(&mb), "--not-in"]).0, 2);
}

#[test]
fn a_file_that_grows_while_it_is_read_stops_the_manifest() {
    let t = Tmp::new("grows");
    let d = t.0.join("tree");
    let f = d.join("big.dat");
    // Over 16 MiB, so it is read as a stream.
    write(&f, &vec![b'x'; 32 << 20], 0o644);
    let stop = std::sync::Arc::new(std::sync::atomic::AtomicBool::new(false));
    let (s2, f2) = (stop.clone(), f.clone());
    // 4 KiB a millisecond, for at most 10 s (40 MiB): a reader that follows
    // the growth still reaches the end, and the disk never fills.
    // The manifest starts only once the writer is appending, so it cannot
    // finish before the writer runs: it reads 32 MiB, and an append comes
    // every millisecond of that read.
    let (started, wait) = std::sync::mpsc::channel();
    let writer = std::thread::spawn(move || {
        let mut o = fs::OpenOptions::new().append(true).open(&f2).unwrap();
        for i in 0..10_000 {
            if s2.load(std::sync::atomic::Ordering::Relaxed) {
                break;
            }
            o.write_all(&[b'y'; 4096]).unwrap();
            if i == 0 {
                started.send(()).unwrap();
            }
            std::thread::sleep(std::time::Duration::from_millis(1));
        }
    });
    wait.recv().unwrap();
    let (code, out, err) = run(&[&s(&d)]);
    stop.store(true, std::sync::atomic::Ordering::Relaxed);
    writer.join().unwrap();
    assert_eq!(code, 2, "a file that grew while read is listed:\n{out}");
    assert!(err.contains("changed while it was read"), "{err}");
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
    let e = log(
        "e.log",
        &base.replace("--kir-data /g", "--kir-data /tmp/2026-11-26"),
    );
    let f = log(
        "f.log",
        &base.replace("--kir-data /g", "--kir-data /tmp/2026-11-27"),
    );
    let (code, out, _) = run(&["--compare-logs", &e, &f]);
    assert_eq!(
        code, 1,
        "retry folders that differ by a date are green:\n{out}"
    );
    assert!(out.contains("DIFF\tmissing\t"), "{out}");
    // A date inside a file name is part of the name, not a time.
    let g = log(
        "g.log",
        &base.replace("Downloading x", "Downloading release-2026-11-26.tar.gz"),
    );
    let h = log(
        "h.log",
        &base.replace("Downloading x", "Downloading release-2026-11-27.tar.gz"),
    );
    let (code, out, _) = run(&["--compare-logs", &g, &h]);
    assert_eq!(
        code, 1,
        "file names that differ by a date are green:\n{out}"
    );
    // So is a date that is a whole folder of a path.
    let k = log(
        "k.log",
        &base.replace(
            "Downloading x",
            "Downloading https://example.org/releases/2026-11-26/m.tar.gz",
        ),
    );
    let l = log(
        "l.log",
        &base.replace(
            "Downloading x",
            "Downloading https://example.org/releases/2026-11-27/m.tar.gz",
        ),
    );
    let (code, out, _) = run(&["--compare-logs", &k, &l]);
    assert_eq!(
        code, 1,
        "paths that differ by a date folder are green:\n{out}"
    );
    // A timestamp that ends a sentence is still a time.
    let i = log(
        "i.log",
        &format!("{base}Finished at 2026-11-26 10:00:00.\n"),
    );
    let j = log(
        "j.log",
        &format!("{base}Finished at 2026-11-27 11:12:13.\n"),
    );
    let (code, out, _) = run(&["--compare-logs", &i, &j]);
    assert_eq!(code, 0, "{out}");
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

/// The og binary cargo built beside this test binary: ab.sh runs `og dev
/// install-manifest`. `cargo test --workspace` builds it; on its own, run
/// `cargo build -p og` first.
fn og_bin() -> PathBuf {
    let exe = std::env::current_exe().unwrap();
    let og = exe.parent().and_then(Path::parent).unwrap().join("og");
    assert!(
        og.is_file(),
        "{} is missing: build og first (cargo build -p og)",
        og.display()
    );
    og
}

/// A stub `scripts/setup.sh` for ab.sh's own tests; `AB_STUB` says what it
/// does to the install folder (its last argument).
const STUB_SETUP: &str = r#"#!/usr/bin/env bash
set -eu
dir=${!#}
case "$AB_STUB" in
  ok)
    mkdir -p "$dir/kir"
    printf 'the KIR data, longer than the seed file\n' >"$dir/kir/kir.dat"
    echo "[OK] KIR data" ;;
  fail)
    mkdir -p "$dir/kir"
    printf 'half of it\n' >"$dir/kir/kir.part"
    echo "ERROR: could not download the KIR data"
    exit 1 ;;
  readonly)
    mkdir -p "$dir/kir"
    printf 'the KIR data, longer than the seed file\n' >"$dir/kir/kir.dat"
    chmod 0444 "$dir/kir/kir.dat"
    echo "[OK] KIR data" ;;
  longname)
    mkdir -p "$dir/kir"
    printf 'the KIR data, longer than the seed file\n' >"$dir/kir/$(printf '%0255d' 0 | tr 0 n)"
    echo "[OK] KIR data" ;;
  newline)
    mkdir -p "$dir/kir"
    printf 'the KIR data, longer than the seed file\n' >"$dir/kir/k"$'\n'
    echo "[OK] KIR data" ;;
  images)
    mkdir -p "$dir/kir"
    printf 'the KIR data, longer than the seed file\n' >"$dir/kir/kir.dat"
    printf 'registry.example/example/kir:1.0\tid1\n' >>"$AB_ENGINE_STATE"
    echo "Status: Downloaded newer image for registry.example/example/kir:1.0"
    # A digest-pinned image and its pinned tag.
    printf 'localhost/pgp-pinned/example/pin:sha256-%064d\tid4\n' 7 >>"$AB_ENGINE_STATE"
    printf 'Status: Downloaded newer image for example/pin@sha256:%064d\n' 7
    # A second one, at another digest.
    printf 'localhost/pgp-pinned/example/pin2:sha256-%064d\tid6\n' 8 >>"$AB_ENGINE_STATE"
    printf 'Status: Downloaded newer image for example/pin2@sha256:%064d\n' 8
    # Another job on the engine adds, once, another tag of the same name
    # and the same name and tag from another registry.
    if [ ! -e "$AB_ENGINE_STATE.other" ]; then
      : >"$AB_ENGINE_STATE.other"
      printf 'registry.example/example/kir:2.0\tid2\nexample/kir:1.0\tid3\n' >>"$AB_ENGINE_STATE"
      # and another image's pinned tag at the same digest.
      printf 'localhost/pgp-pinned/example/other:sha256-%064d\tid5\n' 7 >>"$AB_ENGINE_STATE"
      # and a pinned tag with the name of one image side A logs and the
      # digest of the other.
      printf 'localhost/pgp-pinned/example/pin:sha256-%064d\tid7\n' 8 >>"$AB_ENGINE_STATE"
    fi ;;
  nothing)
    echo "nothing to do" ;;
  foreign)
    mkdir -p "$dir/kir"
    printf 'the KIR data, longer than the seed file\n' >"$dir/kir/kir.dat"
    printf 'written as another user, longer than the KIR data\n' >"$dir/kir/other.dat"
    if [ "$(id -u)" = 0 ]; then chown 65534 "$dir/kir/other.dat"; else sudo -n chown 0 "$dir/kir/other.dat"; fi ;;
esac
"#;

/// Runs ab.sh --bash-vs-bash kir-data in a scratch copy of the repo with
/// the stub setup.sh, a one-file seed and no container engine; its exit
/// code and output. A variable given empty in `env` is unset.
fn ab(tag: &str, stub: &str, env: &[(&str, &str)]) -> (i32, String) {
    ab_with(tag, stub, env, &|_| {})
}

/// [`ab`], with `prep` called on the scratch folder (the repo copy is in
/// its `repo` folder) before ab.sh runs.
fn ab_with(tag: &str, stub: &str, env: &[(&str, &str)], prep: &dyn Fn(&Path)) -> (i32, String) {
    let t = Tmp::new(tag);
    let root = t.0.join("repo");
    for f in [
        "tests/setup/ab.sh",
        "tests/setup/install-manifest-allow.tsv",
    ] {
        write(&root.join(f), &fs::read(repo().join(f)).unwrap(), 0o755);
    }
    write(&root.join("scripts/setup.sh"), STUB_SETUP.as_bytes(), 0o755);
    write(&root.join("rust/og.sh"), b"#!/bin/sh\nexit 0\n", 0o755);
    let seed = t.0.join("seed");
    write(&seed.join("ref.fa"), b"ACGT\n", 0o644);
    prep(&t.0);
    let mut cmd = Command::new("bash");
    cmd.arg(root.join("tests/setup/ab.sh"))
        .args(["--bash-vs-bash", "kir-data"])
        .arg(t.0.join("work"))
        .env("OG", og_bin())
        .env("CONTAINER_ENGINE", "true")
        .env("AB_SEED", &seed)
        .env("AB_STUB", stub)
        .env("TMPDIR", &t.0)
        .env_remove("OG_DOWNLOADS")
        .env_remove("AB_EXPECT_EXIT")
        .env_remove("AB_ALLOW");
    for &(k, v) in env {
        if v.is_empty() {
            cmd.env_remove(k);
        } else {
            cmd.env(k, v);
        }
    }
    let out = cmd.output().unwrap();
    let mut text = format!(
        "{}{}",
        String::from_utf8_lossy(&out.stdout),
        String::from_utf8_lossy(&out.stderr)
    );
    // The scratch folder goes when this returns: say here if the lock stayed.
    if root.join("bin/.og.ab-lock").exists() {
        text.push_str("\n[lock left]\n");
    }
    (out.status.code().unwrap_or(-1), text)
}

#[test]
fn ab_sh_is_green_on_a_real_install_and_plants_its_control_in_it() {
    let (code, out) = ab("ab-ok", "ok", &[]);
    assert_eq!(code, 0, "{out}");
    // The seed's ref.fa is smaller, but the seed did not come from the
    // component: a bit flipped there proves nothing about it.
    assert!(
        out.contains("control: one flipped bit in kir/kir.dat is red"),
        "the control is not in the component's file:\n{out}"
    );
}

#[test]
fn ab_sh_is_red_when_side_a_fails_unless_that_exit_is_expected() {
    let (code, out) = ab("ab-fail", "fail", &[]);
    assert_eq!(code, 1, "both sides failing is green:\n{out}");
    assert!(out.contains("[RED]  side-a-exit"), "{out}");
    let (code, out) = ab("ab-fail-expected", "fail", &[("AB_EXPECT_EXIT", "1")]);
    assert_eq!(code, 0, "an expected exit 1 is red:\n{out}");
    let (code, out) = ab("ab-ok-expected", "ok", &[("AB_EXPECT_EXIT", "1")]);
    assert_eq!(code, 1, "exit 0 where 1 was expected is green:\n{out}");
}

#[test]
fn ab_sh_runs_without_a_seed() {
    let (code, out) = ab("ab-noseed", "ok", &[("AB_SEED", "")]);
    assert_eq!(code, 0, "a run with no seed stops:\n{out}");
    assert!(
        out.contains("control: one flipped bit in kir/kir.dat is red"),
        "{out}"
    );
}

#[test]
fn ab_sh_plants_its_control_in_a_read_only_file() {
    let (code, out) = ab("ab-readonly", "readonly", &[]);
    assert_eq!(code, 0, "a read-only file stops the control:\n{out}");
    assert!(
        out.contains("control: one flipped bit in kir/kir.dat is red"),
        "{out}"
    );
}

#[test]
fn ab_sh_plants_its_control_in_a_name_ending_in_a_newline() {
    let (code, out) = ab("ab-newline", "newline", &[]);
    assert_eq!(code, 0, "{out}");
    assert!(
        out.contains("control: one flipped bit in kir/k\\n is red"),
        "{out}"
    );
}

#[test]
fn ab_sh_plants_its_control_in_a_file_with_a_255_byte_name() {
    let (code, out) = ab("ab-longname", "longname", &[]);
    assert_eq!(code, 0, "{out}");
    let want = format!("control: one flipped bit in kir/{} is red", "n".repeat(255));
    assert!(out.contains(&want), "{out}");
}

#[test]
fn ab_sh_refuses_a_checkout_another_run_uses() {
    let (code, out) = ab_with("ab-locked", "ok", &[], &|t| {
        fs::create_dir_all(t.join("repo/bin/.og.ab-lock")).unwrap();
    });
    assert_eq!(code, 2, "two runs share one bin/.og:\n{out}");
    assert!(
        out.contains("another ab.sh run uses this checkout"),
        "{out}"
    );
    // A refused run leaves the other run's lock; a run removes its own.
    assert!(out.contains("[lock left]"), "{out}");
    let (code, out) = ab("ab-unlocked", "ok", &[]);
    assert_eq!(code, 0, "{out}");
    assert!(
        !out.contains("[lock left]"),
        "the lock outlives its run:\n{out}"
    );
}

#[test]
fn ab_sh_leaves_another_tag_of_the_same_image_name_alone() {
    let t = Tmp::new("ab-engine");
    let engine = t.0.join("engine");
    write(
        &engine,
        b"#!/usr/bin/env bash\nst=$AB_ENGINE_STATE\ncase \"$1\" in\n  image) cat \"$st\" ;;\n  rmi) awk -F'\\t' -v t=\"$2\" '$1 != t' \"$st\" >\"$st.new\" && mv \"$st.new\" \"$st\" ;;\nesac\n",
        0o755,
    );
    let state = t.0.join("images.tsv");
    fs::write(&state, "").unwrap();
    let (_, out) = ab(
        "ab-images",
        "images",
        &[
            ("CONTAINER_ENGINE", &s(&engine)),
            ("AB_ENGINE_STATE", &s(&state)),
        ],
    );
    assert!(out.contains("removed 3 image tag(s)"), "{out}");
    let left = fs::read_to_string(&state).unwrap();
    let other = format!("\nlocalhost/pgp-pinned/example/other:sha256-{:064}\t", 7);
    let crossed = format!("\nlocalhost/pgp-pinned/example/pin:sha256-{:064}\t", 8);
    for tag in [
        "registry.example/example/kir:2.0\t",
        "\nexample/kir:1.0\t",
        other.as_str(),
        crossed.as_str(),
    ] {
        assert!(
            format!("\n{left}").contains(tag),
            "another job's tag {tag:?} was removed:\n{left}\n{out}"
        );
    }
}

#[test]
fn ab_sh_stops_when_the_component_installs_nothing() {
    let (code, out) = ab("ab-nothing", "nothing", &[]);
    assert_eq!(code, 2, "a run that installs nothing passes:\n{out}");
    assert!(out.contains("holds nothing the seed did not"), "{out}");
}

#[test]
fn ab_sh_is_red_on_a_file_another_user_owns() {
    let root = Command::new("id").arg("-u").output().unwrap().stdout == b"0\n";
    let sudo = !root
        && Command::new("sudo")
            .args(["-n", "true"])
            .stderr(std::process::Stdio::null())
            .status()
            .is_ok_and(|s| s.success());
    if !root && !sudo {
        eprintln!("skipped: making a file another user owns needs root or sudo -n");
        return;
    }
    // A volume mounted with owners ignored takes the chown and keeps the
    // owner: nothing to test there.
    let probe = Tmp::new("ab-owner-probe");
    let f = probe.0.join("f");
    write(&f, b"x", 0o644);
    let before = {
        use std::os::unix::fs::MetadataExt;
        fs::metadata(&f).unwrap().uid()
    };
    let chown = if root {
        Command::new("chown").arg("65534").arg(&f).status()
    } else {
        Command::new("sudo")
            .args(["-n", "chown", "0"])
            .arg(&f)
            .status()
    };
    let after = {
        use std::os::unix::fs::MetadataExt;
        fs::metadata(&f).unwrap().uid()
    };
    if !chown.is_ok_and(|s| s.success()) || after == before {
        eprintln!("skipped: the temp folder's filesystem does not keep file owners");
        return;
    }
    let (code, out) = ab("ab-foreign", "foreign", &[]);
    assert_eq!(code, 1, "a file another user owns is green:\n{out}");
    assert!(out.contains("[RED]  foreign-files"), "{out}");
    assert!(out.contains("kir/other.dat"), "{out}");
}

#[test]
fn the_linear_index_length_is_in_the_index_digest() {
    // One sequence, two records, one BGZF block each: the same name and
    // record count, and each index fits both data files. Only the window
    // of the last record differs, 100 against 5,000,000.
    let rows = "x/*.tbi\tindex\ttest\toffsets move with the header's length\t-\n";
    let body = |pos: u32| {
        format!(
            "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=10000000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t10\t.\tA\tG\t.\tPASS\t.\nchr1\t{pos}\t.\tC\tT\t.\tPASS\t.\n"
        )
        .into_bytes()
    };
    let (t, a, b) = pair("tbi-linear");
    vcf_with_index(&a.join("x"), "c.vcf.gz", &body(100));
    vcf_with_index(&b.join("x"), "c.vcf.gz", &body(5_000_000));
    let allow = allow_file(&t, rows);
    let (code, out, err) = compare(&manifest(&a, &[]), &manifest(&b, &[]), Some(&allow));
    assert_eq!(code, 1, "{out}{err}");
    assert!(
        out.contains("DIFF\tcontent\tx/c.vcf.gz.tbi\tindex "),
        "an index whose last record moved by 5 Mb is green:\n{out}"
    );
}

/// One gzip member laid out like a BGZF block, except that its extra field
/// holds a second subfield after BC (XLEN 12). BSIZE, the CRC and ISIZE
/// are right, so any gzip reader decodes it; htslib reads only XLEN 6 as
/// BGZF.
fn xlen12_block(body: &[u8]) -> Vec<u8> {
    let mut e = flate2::write::DeflateEncoder::new(Vec::new(), Compression::default());
    e.write_all(body).unwrap();
    let deflated = e.finish().unwrap();
    let mut crc = flate2::Crc::new();
    crc.update(body);
    let total = 12 + 12 + deflated.len() + 8;
    let mut v = vec![0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 12, 0];
    v.extend(b"BC");
    v.extend(2u16.to_le_bytes());
    v.extend(u16::try_from(total - 1).unwrap().to_le_bytes());
    v.extend(b"XY");
    v.extend(2u16.to_le_bytes());
    v.extend([0, 0]);
    v.extend(&deflated);
    v.extend(crc.sum().to_le_bytes());
    v.extend(u32::try_from(body.len()).unwrap().to_le_bytes());
    v
}

#[test]
fn an_index_into_a_well_formed_block_with_xlen_12_is_invalid() {
    let mut data = xlen12_block(&[7u8; 300]);
    let eof = bgzf(b"");
    data.extend_from_slice(&eof[eof.len() - 28..]);
    // A gzip reader takes the block whole: only XLEN tells it from BGZF.
    let mut plain = Vec::new();
    let mut dec = flate2::bufread::GzDecoder::new(&data[..]);
    std::io::Read::read_to_end(&mut dec, &mut plain).unwrap();
    assert_eq!(plain, [7u8; 300]);
    assert_eq!(
        dec.into_inner().len(),
        28,
        "the block is not one whole member"
    );
    let t = Tmp::new("csi-xlen12-whole");
    let d = t.0.join("d");
    write(&d.join("x/d.bcf"), &data, 0o644);
    write(&d.join("x/d.bcf.csi"), &csi(10, b""), 0o644);
    let text = fs::read_to_string(manifest(&d, &[])).unwrap();
    let row = text
        .lines()
        .find(|l| l.starts_with("x/d.bcf.csi\t"))
        .unwrap_or_else(|| panic!("no index row:\n{text}"));
    assert!(
        row.ends_with(",index=invalid"),
        "an index into a block with XLEN 12 fits:\n{row}"
    );
}

#[test]
fn two_broken_files_under_one_mode_row_differ_when_their_bytes_do() {
    let rows = "x/*.tbi\tindex\ttest\toffsets move with the header's length\t-\nx/*.gz\tgzip\ttest\tthe header time\t-\n";
    // Two stale indexes, neither of which fits a data file.
    let (code, out) = modes_case(
        "invalid-twice",
        "x/d.vcf.gz.tbi",
        b"not an index, side A",
        b"not an index, side B",
        rows,
    );
    assert!(
        out.contains("\tx/d.vcf.gz.tbi\tindex invalid -> invalid: "),
        "{out}"
    );
    assert_eq!(code, 1, "two different invalid indexes are green:\n{out}");
    // Two gzip files that do not decode.
    let (code, out) = modes_case(
        "bad-twice",
        "x/d.gz",
        b"\x1f\x8b\x08\x00 not deflate, side A",
        b"\x1f\x8b\x08\x00 not deflate, side B",
        rows,
    );
    assert!(out.contains("\tx/d.gz\tgzip bad -> bad: "), "{out}");
    assert_eq!(code, 1, "two different bad gzip files are green:\n{out}");
}

#[test]
fn a_manifest_that_lists_one_path_twice_is_refused() {
    let (t, a, _b) = pair("dup-path");
    let m = manifest(&a, &[]);
    let text = fs::read_to_string(&m).unwrap();
    let lines: Vec<&str> = text.lines().collect();
    // After the magic line and the column names: the first entry, twice.
    let mut dup: Vec<&str> = lines[..3].to_vec();
    dup.extend(&lines[2..]);
    let p = t.0.join("dup.tsv");
    fs::write(&p, dup.join("\n") + "\n").unwrap();
    for (x, y) in [(&m, &p), (&p, &m)] {
        let (code, out, err) = compare(x, y, None);
        assert_eq!(
            code, 2,
            "a manifest listing a path twice is read:\n{out}{err}"
        );
        assert!(err.contains("is listed twice"), "{err}");
    }
}

/// Whether DIR's filesystem keeps file owners, by what the system says of
/// its mount, not by ab.sh's own probe: a macOS volume mounted `noowners`
/// does not; any Linux one this suite runs on does.
fn mount_keeps_owners(dir: &Path) -> bool {
    if !cfg!(target_os = "macos") {
        return true;
    }
    let df = Command::new("df").arg("-P").arg(dir).output().unwrap();
    let df = String::from_utf8_lossy(&df.stdout);
    let mp = df
        .lines()
        .nth(1)
        .and_then(|l| l.split_whitespace().nth(5))
        .unwrap()
        .to_owned();
    let mount = Command::new("mount").output().unwrap();
    let line = String::from_utf8_lossy(&mount.stdout)
        .lines()
        .find(|l| l.contains(&format!(" on {mp} (")))
        .unwrap()
        .to_owned();
    !line.contains("noowners")
}

#[test]
fn ab_sh_checks_owners_only_where_the_filesystem_keeps_them() {
    // The real probe, in the temp folder: run this suite once with TMPDIR on
    // a folder that keeps owners and once on one that does not (a macOS
    // volume mounted noowners) to see both branches.
    let (code, out) = ab("ab-owners", "ok", &[]);
    assert_eq!(code, 0, "{out}");
    let root = Command::new("id").arg("-u").output().unwrap().stdout == b"0\n";
    let groups = Command::new("id").arg("-G").output().unwrap().stdout;
    let one_group = String::from_utf8_lossy(&groups).split_whitespace().count() < 2;
    if !mount_keeps_owners(&std::env::temp_dir()) {
        assert!(
            out.contains("[SKIP] foreign-files: this filesystem does not keep file owners\n"),
            "{out}"
        );
        assert!(!out.contains("[OK]   foreign-files"), "{out}");
        assert!(
            out.contains("every check is green, foreign-files skipped."),
            "{out}"
        );
    } else if !root && one_group {
        assert!(
            out.contains(
                "[SKIP] foreign-files: cannot tell whether this filesystem keeps file owners"
            ),
            "{out}"
        );
    } else {
        assert!(out.contains("[OK]   foreign-files\n"), "{out}");
        assert!(!out.contains("[SKIP]"), "{out}");
    }
}

#[test]
fn ab_sh_skips_the_owner_check_where_a_changed_owner_does_not_stay() {
    // A chgrp and a chown that do what they do on a volume that ignores
    // owners: succeed, and change nothing.
    let fake = Tmp::new("ab-noowners-bin");
    for c in ["chgrp", "chown"] {
        write(&fake.0.join(c), b"#!/bin/sh\nexit 0\n", 0o755);
    }
    let path = format!("{}:{}", s(&fake.0), std::env::var("PATH").unwrap());
    let (code, out) = ab("ab-noowners", "ok", &[("PATH", &path)]);
    assert_eq!(code, 0, "{out}");
    assert!(
        out.contains("[SKIP] foreign-files: this filesystem does not keep file owners\n"),
        "{out}"
    );
    assert!(!out.contains("[OK]   foreign-files"), "{out}");
    assert!(
        out.contains("every check is green, foreign-files skipped."),
        "{out}"
    );
}
