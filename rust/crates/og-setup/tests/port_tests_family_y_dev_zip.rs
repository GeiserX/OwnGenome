//! `og setup unzip` and `og setup genetic-map` on zips written here. It
//! replaces `tests/test_og_setup_zip.py`, which built the zips with
//! Python's `zipfile` and compared og with `extractall` and with the inline
//! programs `scripts/setup.sh` ran before og. The expected trees below are
//! what those wrote, by content:
//!
//! 1. unzip: stored and deflated members, ZIP64 records, nested and empty
//!    directories, a member named twice (the last one wins), names with `.`
//!    and empty parts, a non-ASCII name: the same tree, file for file and
//!    byte for byte, as `extractall` writes;
//! 2. unzip refuses, with exit 2 and nothing written, a member whose name
//!    climbs out of the folder (`../x`, `/x`), where `extractall` silently
//!    writes it elsewhere; and a member with a flipped byte fails its CRC
//!    check (exit 2);
//! 3. genetic-map: CHROM.map is the member byte for byte and
//!    CHROM.rfmix.map fields 1, 4 and 3 of each non-blank line (whitespace
//!    runs, tabs, CRLF and a last line without a newline), for two
//!    chromosomes in one call; without `--rfmix` only the .map is written;
//! 4. genetic-map stops with exit 2 on a row with fewer than four fields
//!    and on a chromosome the zip lacks, and leaves the files of that
//!    chromosome as they were;
//! 5. the control: a tree with one byte changed is told apart from the
//!    expected one, so the comparison can fail.
//!
//! The cases that read stderr run og in a child process: this test binary
//! itself, with [`CHILD`] set (see [`child_entry`]).

use flate2::write::DeflateEncoder;
use flate2::{Compression, Crc};
use std::collections::BTreeMap;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

const CHILD: &str = "OG_TEST_SETUP_ZIP_CHILD";

/// In the child process only: runs `og setup <tool>` with the arguments in
/// [`CHILD`] (separated by U+001F) and exits with its code, so the parent
/// reads its stderr. In a normal test run it does nothing.
#[test]
fn child_entry() {
    let Ok(joined) = std::env::var(CHILD) else {
        return;
    };
    let args: Vec<String> = joined.split('\u{1f}').map(str::to_owned).collect();
    let code = match args[0].as_str() {
        "unzip" => og_setup::unzip::cli(&args[1..]),
        _ => og_setup::genetic_map::cli(&args[1..]),
    };
    std::process::exit(i32::from(code));
}

/// `og setup ARGS` in a child process: (exit code, stderr).
fn og_setup(args: &[&str]) -> (i32, String) {
    let out = Command::new(std::env::current_exe().unwrap())
        .args(["child_entry", "--exact", "--nocapture", "--test-threads=1"])
        .env(CHILD, args.join("\u{1f}"))
        .output()
        .unwrap();
    (
        out.status.code().unwrap_or(-1),
        String::from_utf8_lossy(&out.stderr).into_owned(),
    )
}

#[derive(Clone, Copy, PartialEq)]
enum Method {
    /// A directory entry (`ZipInfo("name/")`, no content).
    Dir,
    Stored,
    Deflated,
}

fn le16(out: &mut Vec<u8>, v: u16) {
    out.extend_from_slice(&v.to_le_bytes());
}
fn le32(out: &mut Vec<u8>, v: u32) {
    out.extend_from_slice(&v.to_le_bytes());
}
fn le64(out: &mut Vec<u8>, v: u64) {
    out.extend_from_slice(&v.to_le_bytes());
}

/// A zip of `members`, laid out as Python's `zipfile` writes one: a local
/// header before each member, then the central directory, names flagged
/// UTF-8 when they are not ASCII. `zip64` puts every size and offset in
/// ZIP64 extra fields (all ones in the plain fields) and adds the ZIP64 end
/// record and its locator.
fn build(members: &[(&str, &[u8], Method)], zip64: bool) -> Vec<u8> {
    let (mut out, mut cd) = (Vec::new(), Vec::new());
    let version: u16 = if zip64 { 45 } else { 20 };
    for &(name, data, method) in members {
        let mut crc = Crc::new();
        crc.update(data);
        let body = if method == Method::Deflated {
            let mut e = DeflateEncoder::new(Vec::new(), Compression::default());
            e.write_all(data).unwrap();
            e.finish().unwrap()
        } else {
            data.to_vec()
        };
        let code: u16 = if method == Method::Deflated { 8 } else { 0 };
        let flags: u16 = if name.is_ascii() { 0 } else { 0x800 };
        let offset = out.len() as u64;
        let (size, csize) = (data.len() as u64, body.len() as u64);
        let small = |v: u64| if zip64 { u32::MAX } else { v as u32 };
        let mut local_extra = Vec::new();
        let mut central_extra = Vec::new();
        if zip64 {
            le16(&mut local_extra, 1);
            le16(&mut local_extra, 16);
            le64(&mut local_extra, size);
            le64(&mut local_extra, csize);
            le16(&mut central_extra, 1);
            le16(&mut central_extra, 24);
            le64(&mut central_extra, size);
            le64(&mut central_extra, csize);
            le64(&mut central_extra, offset);
        }
        le32(&mut out, 0x0403_4b50);
        for v in [version, flags, code, 0, 0x21] {
            le16(&mut out, v);
        }
        le32(&mut out, crc.sum());
        le32(&mut out, small(csize));
        le32(&mut out, small(size));
        le16(&mut out, name.len() as u16);
        le16(&mut out, local_extra.len() as u16);
        out.extend_from_slice(name.as_bytes());
        out.extend_from_slice(&local_extra);
        out.extend_from_slice(&body);
        le32(&mut cd, 0x0201_4b50);
        for v in [version, version, flags, code, 0, 0x21] {
            le16(&mut cd, v);
        }
        le32(&mut cd, crc.sum());
        le32(&mut cd, small(csize));
        le32(&mut cd, small(size));
        le16(&mut cd, name.len() as u16);
        le16(&mut cd, central_extra.len() as u16);
        for v in [0, 0, 0] {
            le16(&mut cd, v);
        }
        le32(&mut cd, if method == Method::Dir { 0x10 } else { 0 });
        le32(&mut cd, small(offset));
        cd.extend_from_slice(name.as_bytes());
        cd.extend_from_slice(&central_extra);
    }
    let cd_offset = out.len() as u64;
    out.extend_from_slice(&cd);
    let n = members.len() as u64;
    if zip64 {
        let eocd64 = out.len() as u64;
        le32(&mut out, 0x0606_4b50);
        le64(&mut out, 44);
        for v in [45, 45] {
            le16(&mut out, v);
        }
        le32(&mut out, 0);
        le32(&mut out, 0);
        for v in [n, n, cd.len() as u64, cd_offset] {
            le64(&mut out, v);
        }
        le32(&mut out, 0x0706_4b50);
        le32(&mut out, 0);
        le64(&mut out, eocd64);
        le32(&mut out, 1);
    }
    le32(&mut out, 0x0605_4b50);
    le32(&mut out, 0);
    let count = if zip64 { u16::MAX } else { n as u16 };
    le16(&mut out, count);
    le16(&mut out, count);
    le32(&mut out, cd.len() as u32);
    le32(&mut out, if zip64 { u32::MAX } else { cd_offset as u32 });
    le16(&mut out, 0);
    out
}

/// {relative path: content, or None for a directory} under `root`.
fn tree(root: &Path) -> BTreeMap<String, Option<Vec<u8>>> {
    let mut out = BTreeMap::new();
    let mut todo = vec![root.to_path_buf()];
    while let Some(d) = todo.pop() {
        for e in fs::read_dir(&d).unwrap() {
            let p = e.unwrap().path();
            let rel = p.strip_prefix(root).unwrap().to_string_lossy().into_owned();
            if p.is_dir() {
                out.insert(rel, None);
                todo.push(p);
            } else {
                out.insert(rel, Some(fs::read(&p).unwrap()));
            }
        }
    }
    out
}

fn scratch(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("og-setup-zip-{}-{name}", std::process::id()));
    let _ = fs::remove_dir_all(&d);
    fs::create_dir_all(&d).unwrap();
    d
}

/// 300,000 bytes that do not compress, then 2,000,000 that do.
fn big() -> Vec<u8> {
    let mut s = 0x9e37_79b9_7f4a_7c15u64;
    let mut v: Vec<u8> = (0..300_000)
        .map(|_| {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            (s >> 24) as u8
        })
        .collect();
    v.resize(2_300_000, b'A');
    v
}

fn file(v: &[u8]) -> Option<Vec<u8>> {
    Some(v.to_vec())
}

fn diff(
    a: &BTreeMap<String, Option<Vec<u8>>>,
    b: &BTreeMap<String, Option<Vec<u8>>>,
) -> Vec<String> {
    a.keys()
        .chain(b.keys())
        .filter(|k| a.get(*k) != b.get(*k))
        .cloned()
        .collect()
}

#[test]
fn unzip_writes_the_tree_extractall_writes() {
    let tmp = scratch("unzip");
    let big = big();
    let members: [(&str, &[u8], Method); 8] = [
        ("2602_phenotype/", b"", Method::Dir),
        ("2602_phenotype/phenotype.mv.db", &big, Method::Deflated),
        (
            "2602_phenotype/rw/stored.txt",
            b"stored as is\n",
            Method::Stored,
        ),
        ("2602_phenotype/empty/", b"", Method::Dir),
        ("2602_phenotype/twice.txt", b"first\n", Method::Deflated),
        ("2602_phenotype/twice.txt", b"second\n", Method::Stored),
        (
            "2602_phenotype/./odd//parts.txt",
            b"parts\n",
            Method::Deflated,
        ),
        ("2602_phenotype/café.txt", b"utf-8 name\n", Method::Deflated),
    ];
    // What zipfile.ZipFile(z).extractall(dest) wrote for these members.
    let want: BTreeMap<String, Option<Vec<u8>>> = [
        ("2602_phenotype", None),
        ("2602_phenotype/phenotype.mv.db", file(&big)),
        ("2602_phenotype/rw", None),
        ("2602_phenotype/rw/stored.txt", file(b"stored as is\n")),
        ("2602_phenotype/empty", None),
        ("2602_phenotype/twice.txt", file(b"first\n")),
        ("2602_phenotype/odd", None),
        ("2602_phenotype/odd/parts.txt", file(b"parts\n")),
        ("2602_phenotype/café.txt", file(b"utf-8 name\n")),
    ]
    .into_iter()
    .map(|(k, v)| (k.to_owned(), v))
    .collect();
    for zip64 in [false, true] {
        let z = tmp.join(format!("data{}.zip", u8::from(zip64)));
        fs::write(&z, build(&members, zip64)).unwrap();
        let dest = tmp.join(format!("og{}", u8::from(zip64)));
        let (code, err) = og_setup(&["unzip", z.to_str().unwrap(), dest.to_str().unwrap()]);
        let kind = if zip64 { "ZIP64" } else { "plain" };
        assert_eq!(code, 0, "unzip ({kind}): exit 0: {err}");
        let got = tree(&dest);
        assert!(
            got == want,
            "unzip ({kind}): the same tree as extractall, byte for byte: {:?}",
            diff(&got, &want)
        );
        if !zip64 {
            let mut changed = got.clone();
            changed.insert("2602_phenotype/twice.txt".into(), file(b"secone\n"));
            assert!(
                changed != want,
                "control: a tree with one byte changed differs"
            );
        }
    }

    for bad in ["../evil.txt", "top/../../evil.txt", "/abs/evil.txt"] {
        let z = tmp.join("bad.zip");
        fs::write(
            &z,
            build(
                &[
                    ("fine.txt", b"x", Method::Stored),
                    (bad, b"y", Method::Stored),
                ],
                false,
            ),
        )
        .unwrap();
        let dest = tmp.join("bad/out");
        let (code, err) = og_setup(&["unzip", z.to_str().unwrap(), dest.to_str().unwrap()]);
        assert!(
            code == 2 && !dest.exists() && !tmp.join("bad/evil.txt").exists(),
            "unzip refuses the member {bad:?} (exit 2, nothing written): {code} {err}"
        );
    }

    let z = tmp.join("crc.zip");
    let mut data = build(&[("a.txt", b"some content here", Method::Stored)], false);
    let at = data
        .windows(17)
        .position(|w| w == b"some content here")
        .unwrap();
    data[at] ^= 1;
    fs::write(&z, data).unwrap();
    let (code, err) = og_setup(&[
        "unzip",
        z.to_str().unwrap(),
        tmp.join("crc").to_str().unwrap(),
    ]);
    assert!(
        code == 2 && err.contains("CRC"),
        "unzip: a flipped byte fails the CRC check (exit 2): {code} {err}"
    );
    fs::remove_dir_all(tmp).unwrap();
}

const MAP20: &str = "chr20 . 0.0 60000\nchr20\t.\t0.12\t70001\n\n   \nchr20  .  1.5   80000\n";
const MAP21: &str = "chr21 . 0 1\r\nchr21 . 0.5 2\r\nchr21 . 1.25 3";

#[test]
fn genetic_map_cuts_the_maps_setup_cut_before_og() {
    let tmp = scratch("maps");
    let z = tmp.join("plink.GRCh38.map.zip");
    fs::write(
        &z,
        build(
            &[
                ("chr_in_chrom_field/", b"", Method::Dir),
                (
                    "chr_in_chrom_field/plink.chrchr20.GRCh38.map",
                    MAP20.as_bytes(),
                    Method::Deflated,
                ),
                (
                    "chr_in_chrom_field/plink.chrchr21.GRCh38.map",
                    MAP21.as_bytes(),
                    Method::Deflated,
                ),
                (
                    "chr_in_chrom_field/plink.chrchr22.GRCh38.map",
                    b"chr22 . 0.0 1\nchr22 . 0.4\n",
                    Method::Deflated,
                ),
            ],
            false,
        ),
    )
    .unwrap();
    let zs = z.to_str().unwrap();
    let ours = tmp.join("og-maps");
    fs::create_dir_all(&ours).unwrap();
    let (code, err) = og_setup(&[
        "genetic-map",
        "--zip",
        zs,
        "--out",
        ours.to_str().unwrap(),
        "--rfmix",
        "chr20",
        "chr21",
    ]);
    assert_eq!(code, 0, "genetic-map --rfmix chr20 chr21: exit 0: {err}");
    // What setup.sh's inline program wrote before og: the member as it is,
    // and fields 1, 4 and 3 of each non-blank line.
    let want: BTreeMap<String, Option<Vec<u8>>> = [
        ("chr20.map", MAP20),
        (
            "chr20.rfmix.map",
            "chr20\t60000\t0.0\nchr20\t70001\t0.12\nchr20\t80000\t1.5\n",
        ),
        ("chr21.map", MAP21),
        (
            "chr21.rfmix.map",
            "chr21\t1\t0\nchr21\t2\t0.5\nchr21\t3\t1.25\n",
        ),
    ]
    .into_iter()
    .map(|(k, v)| (k.to_owned(), file(v.as_bytes())))
    .collect();
    let got = tree(&ours);
    assert!(
        got == want,
        "genetic-map: the same .map and .rfmix.map files as the Python: {:?}",
        diff(&got, &want)
    );

    let only = tmp.join("only-map");
    fs::create_dir_all(&only).unwrap();
    let (code, err) = og_setup(&[
        "genetic-map",
        "--zip",
        zs,
        "--out",
        only.to_str().unwrap(),
        "chr20",
    ]);
    assert!(
        code == 0 && tree(&only).keys().collect::<Vec<_>>() == ["chr20.map"],
        "genetic-map without --rfmix writes only chr20.map: {code} {err}"
    );

    let keep = tmp.join("keep");
    fs::create_dir_all(&keep).unwrap();
    fs::write(keep.join("chr22.rfmix.map"), "old\n").unwrap();
    let before = tree(&keep);
    let (code, err) = og_setup(&[
        "genetic-map",
        "--zip",
        zs,
        "--out",
        keep.to_str().unwrap(),
        "--rfmix",
        "chr22",
    ]);
    assert!(
        code == 2 && err.contains("line 2"),
        "genetic-map: a row with three fields stops it (exit 2) and names the line: {code} {err}"
    );
    assert!(
        tree(&keep) == before,
        "genetic-map: and leaves that chromosome's files as they were"
    );
    let (code, err) = og_setup(&[
        "genetic-map",
        "--zip",
        zs,
        "--out",
        keep.to_str().unwrap(),
        "chrX",
    ]);
    assert_eq!(
        code, 2,
        "genetic-map: a chromosome the zip lacks stops it: {err}"
    );
    fs::remove_dir_all(tmp).unwrap();
}
