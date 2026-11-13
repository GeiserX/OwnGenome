//! `og dev parity-compare --self-test`: for every kind, a pair that must
//! compare the same (it differs only where the kind allows) and planted
//! differences that must each be reported. A comparator that cannot fail is
//! not a check.

use super::{Kind, compare};
use std::io::Write;
use std::path::PathBuf;

/// One pair of files and what the comparison must say.
struct Case {
    label: &'static str,
    kind: &'static str,
    oracle: Vec<u8>,
    new: Vec<u8>,
    same: bool,
}

fn gz(data: &[u8], level: u32) -> Vec<u8> {
    let mut e = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::new(level));
    e.write_all(data).expect("in memory");
    e.finish().expect("in memory")
}

fn png(rgba: &[u8], w: u32, h: u32, fast: bool) -> Vec<u8> {
    let mut out = Vec::new();
    {
        let mut e = png::Encoder::new(&mut out, w, h);
        e.set_color(png::ColorType::Rgba);
        e.set_depth(png::BitDepth::Eight);
        e.set_compression(if fast {
            png::Compression::Fastest
        } else {
            png::Compression::High
        });
        e.set_filter(if fast {
            png::Filter::NoFilter
        } else {
            png::Filter::Paeth
        });
        let mut w = e.write_header().expect("in memory");
        w.write_image_data(rgba).expect("in memory");
    }
    out
}

fn case(label: &'static str, kind: &'static str, oracle: &[u8], new: &[u8], same: bool) -> Case {
    Case {
        label,
        kind,
        oracle: oracle.to_vec(),
        new: new.to_vec(),
        same,
    }
}

fn cases() -> Vec<Case> {
    let tsv = b"# made by x\nid\tscore\tnote\nr1\t0.5\tok\nr2\t1.25\t\"a, b\"\n";
    let pixels: Vec<u8> = (0..64u8).collect();
    let mut planted_pixel = pixels.clone();
    planted_pixel[17] ^= 1;
    let zfiles: &[(&str, &[u8])] = &[("a.tsv", b"x\t1\n"), ("b.txt", b"hello\n")];
    let zchanged: &[(&str, &[u8])] = &[("a.tsv", b"x\t2\n"), ("b.txt", b"hello\n")];
    let zmissing: &[(&str, &[u8])] = &[("a.tsv", b"x\t1\n")];
    let vcf = b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t10\t.\tA\tG\t.\tPASS\tDP=3;AF=0.5\n";
    vec![
        case("bytes: equal files", "bytes", b"a\nb\n", b"a\nb\n", true),
        case("bytes: one changed byte", "bytes", b"a\nb\n", b"a\nc\n", false),
        case("bytes: a missing last byte", "bytes", b"a\nb\n", b"a\nb", false),
        case("gz-content: same content, other compression", "gz-content", &gz(b"x\n", 1), &gz(b"x\n", 9), true),
        case("gz-content: one changed byte inside", "gz-content", &gz(b"x\n", 6), &gz(b"y\n", 6), false),
        case("gz-content: new is not gzip", "gz-content", &gz(b"x\n", 6), b"x\n", false),
        case("content: gzip, same content, other compression", "content", &gz(b"x\n", 1), &gz(b"x\n", 9), true),
        case("content: gzip, one changed byte inside", "content", &gz(b"x\n", 6), &gz(b"y\n", 6), false),
        case("content: gzip expected, plain written", "content", &gz(b"x\n", 6), b"x\n", false),
        case("content: plain, one changed byte", "content", b"a\tb\n", b"a\tc\n", false),
        case("content: plain, gzip of the same bytes written", "content", b"x\n", &gz(b"x\n", 6), false),
        case("zip-content: same members, other timestamps and deflate", "zip-content", &crate::zip::write(zfiles, false, 1), &crate::zip::write(zfiles, true, 2), true),
        case("zip-content: one changed member", "zip-content", &crate::zip::write(zfiles, true, 1), &crate::zip::write(zchanged, true, 1), false),
        case("zip-content: a missing member", "zip-content", &crate::zip::write(zfiles, true, 1), &crate::zip::write(zmissing, true, 1), false),
        case("png-pixels: same pixels, other filter and compression", "png-pixels", &png(&pixels, 4, 4, true), &png(&pixels, 4, 4, false), true),
        case("png-pixels: one changed pixel", "png-pixels", &png(&pixels, 4, 4, true), &png(&planted_pixel, 4, 4, true), false),
        case("png-pixels: another size", "png-pixels", &png(&pixels, 4, 4, true), &png(&pixels, 8, 2, true), false),
        case("tsv: a quoted cell left unquoted", "tsv", tsv, b"# made by x\nid\tscore\tnote\nr1\t0.5\tok\nr2\t1.25\ta, b\n", true),
        case("tsv: one changed cell", "tsv", tsv, b"# made by x\nid\tscore\tnote\nr1\t0.50\tok\nr2\t1.25\t\"a, b\"\n", false),
        case("tsv: a missing row", "tsv", tsv, b"# made by x\nid\tscore\tnote\nr1\t0.5\tok\n", false),
        case("tsv: rows in another order", "tsv", tsv, b"# made by x\nid\tscore\tnote\nr2\t1.25\t\"a, b\"\nr1\t0.5\tok\n", false),
        case("tsv:unordered: rows in another order", "tsv:unordered", tsv, b"# made by x\nid\tscore\tnote\nr2\t1.25\t\"a, b\"\nr1\t0.5\tok\n", true),
        case("tsv:keys=id: one changed cell", "tsv:keys=id", tsv, b"# made by x\nid\tscore\tnote\nr2\t1.25\t\"a, b\"\nr1\t0.6\tok\n", false),
        case("tsv: a changed # line", "tsv", tsv, b"# made by y\nid\tscore\tnote\nr1\t0.5\tok\nr2\t1.25\t\"a, b\"\n", false),
        case("kv: equal, CRLF line endings", "kv", b"reads\t10\nbases\t20\n", b"reads\t10\r\nbases\t20\r\n", true),
        case("kv: reordered keys", "kv", b"reads\t10\nbases\t20\n", b"bases\t20\nreads\t10\n", false),
        case("kv: one changed value", "kv", b"reads\t10\nbases\t20\n", b"reads\t10\nbases\t21\n", false),
        case("vcf: equal, one gzip", "vcf", vcf, &gz(vcf, 6), true),
        case("vcf: INFO keys reordered", "vcf", vcf, b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t10\t.\tA\tG\t.\tPASS\tAF=0.5;DP=3\n", false),
        case("vcf: a changed header line", "vcf", vcf, b"##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t10\t.\tA\tG\t.\tPASS\tDP=3;AF=0.5\n", false),
        case("json: equal, other spacing", "json", br#"{"a": 1, "b": [0.1, "x"]}"#, br#"{"a":1,"b":[0.1,"x"]}"#, true),
        case("json: a float 1 ulp off", "json", br#"{"a": 0.1}"#, br#"{"a": 0.10000000000000002}"#, false),
        case("json: reordered keys", "json", br#"{"a": 1, "b": 2}"#, br#"{"b": 2, "a": 1}"#, false),
        case("json: int against float", "json", br#"{"a": 1}"#, br#"{"a": 1.0}"#, false),
        case("json: a missing array item", "json", br#"{"a": [1, 2]}"#, br#"{"a": [1]}"#, false),
        case("json:mask=generated_utc: only the masked key differs", "json:mask=generated_utc", br#"{"generated_utc": "x", "a": 1}"#, br#"{"generated_utc": "y", "a": 1}"#, true),
        case("json:paths=a: a difference outside the path", "json:paths=a", br#"{"a": {"b": 1}, "c": 1}"#, br#"{"a": {"b": 1}, "c": 2}"#, true),
        case("json:paths=a: a difference inside the path", "json:paths=a", br#"{"a": {"b": 1}}"#, br#"{"a": {"b": 2}}"#, false),
        case("json:tol=1e-12: a float inside the tolerance", "json:tol=1e-12", br#"{"eta": 0.1}"#, br#"{"eta": 0.10000000000000002}"#, true),
        case("json:tol=1e-12: a float outside the tolerance", "json:tol=1e-12", br#"{"eta": 0.1}"#, br#"{"eta": 0.1001}"#, false),
        case("mask: only a masked line differs", "mask:^Date: ", b"Date: 2026-01-01\nx\n", b"Date: 2026-11-10\nx\n", true),
        case("mask: an unmasked line differs", "mask:^Date: ", b"Date: 2026-01-01\nx\n", b"Date: 2026-01-01\ny\n", false),
        case("stderr-first-line: only a later line differs", "stderr-first-line", b"ERROR: no input\ntrace a\n", b"ERROR: no input\ntrace b\n", true),
        case("stderr-first-line: the first line differs", "stderr-first-line", b"ERROR: no input\n", b"ERROR: no output\n", false),
    ]
}

/// Runs every case; prints one line each; 0 when every one came out as it
/// must, 1 otherwise.
pub fn run() -> u8 {
    let dir: PathBuf =
        std::env::temp_dir().join(format!("og-parity-self-test-{}", std::process::id()));
    if let Err(e) = std::fs::create_dir_all(&dir) {
        eprintln!("og dev parity-compare --self-test: {}: {e}", dir.display());
        return 2;
    }
    let mut bad = 0;
    let all = cases();
    for (i, c) in all.iter().enumerate() {
        let (a, b) = (
            dir.join(format!("{i}.oracle")),
            dir.join(format!("{i}.new")),
        );
        let written = std::fs::write(&a, &c.oracle).and_then(|()| std::fs::write(&b, &c.new));
        let got = match (written, Kind::parse(c.kind)) {
            (Err(e), _) => Err(e.to_string()),
            (_, Err(e)) => Err(e),
            (Ok(()), Ok(k)) => compare(&k, &a, &b),
        };
        let line = match (&got, c.same) {
            (Ok(None), true) => format!("ok      same as it must be: {}", c.label),
            (Ok(Some(d)), false) => format!("ok      caught: {} ({d})", c.label),
            (Ok(None), false) => format!("MISSED  not reported: {}", c.label),
            (Ok(Some(d)), true) => format!(
                "WRONG   reported a difference it must allow: {} ({d})",
                c.label
            ),
            (Err(e), _) => format!("ERROR   {}: {e}", c.label),
        };
        let fine = line.starts_with("ok");
        bad += usize::from(!fine);
        println!("{line}");
    }
    let _ = std::fs::remove_dir_all(&dir);
    if bad > 0 {
        println!(
            "parity-compare self-test: {bad} of {} cases failed",
            all.len()
        );
        1
    } else {
        println!(
            "parity-compare self-test: all {} cases came out as they must",
            all.len()
        );
        0
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn self_test_passes() {
        assert_eq!(super::run(), 0);
    }
}
