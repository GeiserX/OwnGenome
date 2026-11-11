//! `og setup genetic-map --zip ZIP --out DIR [--rfmix] CHROM...`: the PLINK
//! GRCh38 genetic map of each chromosome, out of Beagle's
//! `plink.GRCh38.map.zip`, for step 26b (Beagle and RFMix) and step 51
//! (MoChA).
//!
//! For each CHROM it writes `DIR/CHROM.map`, the member
//! `chr_in_chrom_field/plink.chrCHROM.GRCh38.map` (CHROM is `chr20`, so the member is `plink.chrchr20...`) byte for byte, and with
//! `--rfmix` also `DIR/CHROM.rfmix.map`, RFMix's column order: chromosome,
//! position, cM (fields 1, 4 and 3 of each non-blank line, split on
//! whitespace as Python's `str.split()` splits). The zip is opened once.
//!
//! Each file is written under a temporary name and renamed into place, the
//! `.rfmix.map` last, and only after every row of the chromosome has been
//! read: a kill or an error leaves the old file or none, never half of one.
//! The Python this replaces wrote both with a plain `open().write()`, so a
//! kill could leave a truncated `.rfmix.map` that the next setup took as
//! done. A row with fewer than four fields now stops with its line number
//! (it was a traceback after `.map` was written and `.rfmix.map` emptied).

use crate::zip::{Archive, ZipError};
use og_io::AtomicFile;
use std::io::Write;
use std::path::Path;

/// Python's `str.splitlines()`: the lines without their ends.
fn splitlines(s: &str) -> Vec<&str> {
    let mut out = Vec::new();
    let mut start = 0;
    let mut it = s.char_indices().peekable();
    while let Some((i, c)) = it.next() {
        let end = match c {
            '\r' => {
                if let Some(&(_, '\n')) = it.peek() {
                    it.next();
                    i + 2
                } else {
                    i + 1
                }
            }
            '\n' | '\x0b' | '\x0c' | '\x1c' | '\x1d' | '\x1e' | '\u{85}' | '\u{2028}'
            | '\u{2029}' => i + c.len_utf8(),
            _ => continue,
        };
        out.push(&s[start..i]);
        start = end;
    }
    if start < s.len() {
        out.push(&s[start..]);
    }
    out
}

/// RFMix's columns of a PLINK map: `chrom\tpos\tcM\n` per non-blank line.
pub fn rfmix_rows(text: &str) -> Result<String, String> {
    let mut out = String::with_capacity(text.len());
    for (n, line) in splitlines(text).into_iter().enumerate() {
        let r = pycompat::text::split(line, None, None);
        if r.is_empty() {
            continue;
        }
        if r.len() < 4 {
            return Err(format!(
                "line {}: {} fields, a PLINK map row has 4",
                n + 1,
                r.len()
            ));
        }
        out.push_str(r[0]);
        out.push('\t');
        out.push_str(r[3]);
        out.push('\t');
        out.push_str(r[2]);
        out.push('\n');
    }
    Ok(out)
}

fn write(path: &Path, bytes: &[u8]) -> Result<AtomicFile, String> {
    let fail = |e: std::io::Error| format!("{}: {e}", path.display());
    let mut f = AtomicFile::create(path).map_err(fail)?;
    f.write_all(bytes).map_err(fail)?;
    Ok(f)
}

fn run(zip: &str, out: &str, rfmix: bool, chroms: &[String]) -> Result<(), (u8, String)> {
    let mut a = Archive::open(zip).map_err(|e| (2, format!("{zip}: {e}")))?;
    for c in chroms {
        let name = format!("chr_in_chrom_field/plink.chr{c}.GRCh38.map");
        let e = a.find(&name).cloned().ok_or_else(|| {
            (
                2,
                format!("{zip}: there is no item named '{name}' in the archive"),
            )
        })?;
        // The declared size bounds what copy_to writes, but it is not
        // checked yet: reserve at most 64 MiB from it (a map is a few MB).
        let mut bytes = Vec::with_capacity(usize::try_from(e.size.min(64 << 20)).unwrap_or(0));
        a.copy_to(&e, &mut bytes).map_err(|err| match err {
            ZipError::Io(x) => (1, format!("{zip}: {x}")),
            ZipError::Bad(m) => (2, format!("{zip}: {m}")),
        })?;
        let text = std::str::from_utf8(&bytes)
            .map_err(|err| (2, format!("{zip}: {name} is not UTF-8: {err}")))?;
        let rows = if rfmix {
            Some(rfmix_rows(text).map_err(|m| (2, format!("{zip}: {name}: {m}")))?)
        } else {
            None
        };
        let dir = Path::new(out);
        let map = write(&dir.join(format!("{c}.map")), &bytes).map_err(|m| (1, m))?;
        let rf = match &rows {
            Some(r) => {
                Some(write(&dir.join(format!("{c}.rfmix.map")), r.as_bytes()).map_err(|m| (1, m))?)
            }
            None => None,
        };
        map.commit()
            .map_err(|err| (1, format!("{}/{c}.map: {err}", dir.display())))?;
        if let Some(f) = rf {
            f.commit()
                .map_err(|err| (1, format!("{}/{c}.rfmix.map: {err}", dir.display())))?;
        }
    }
    Ok(())
}

const USAGE: &str = "usage: og setup genetic-map --zip ZIP --out DIR [--rfmix] CHROM...";

/// The command line; returns the exit code: 0, 1 on an I/O error, 2 on a
/// bad zip, a missing chromosome, a bad row or a wrong command line.
pub fn cli(args: &[String]) -> u8 {
    let (mut zip, mut out, mut rfmix, mut chroms) = (None, None, false, Vec::new());
    let mut it = args.iter();
    while let Some(a) = it.next() {
        match a.as_str() {
            "--zip" => zip = it.next().cloned(),
            "--out" => out = it.next().cloned(),
            "--rfmix" => rfmix = true,
            s if s.starts_with('-') => {
                eprintln!("{USAGE}\nog setup genetic-map: unknown option {s}");
                return 2;
            }
            s => chroms.push(s.to_owned()),
        }
    }
    let (Some(zip), Some(out)) = (zip, out) else {
        eprintln!("{USAGE}");
        return 2;
    };
    if chroms.is_empty() {
        eprintln!("{USAGE}");
        return 2;
    }
    match run(&zip, &out, rfmix, &chroms) {
        Ok(()) => 0,
        Err((code, msg)) => {
            eprintln!("og setup genetic-map: {msg}");
            code
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::zip::tests::build;
    use std::fs;
    use std::path::PathBuf;

    const MAP: &str = "chr20 . 0.0 60000\nchr20\t.\t0.12\t70001\n\nchr20  .  1.5   80000\n";

    fn scratch(name: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!("og-gmap-{}-{name}", std::process::id()));
        let _ = fs::remove_dir_all(&d);
        fs::create_dir_all(&d).unwrap();
        d
    }

    fn zip_with(dir: &Path, members: &[(&str, &[u8], bool)]) -> String {
        let p = dir.join("m.zip");
        fs::write(&p, build(members, false)).unwrap();
        p.display().to_string()
    }

    fn args(z: &str, d: &Path, rest: &[&str]) -> Vec<String> {
        let mut v = vec![
            "--zip".into(),
            z.into(),
            "--out".into(),
            d.display().to_string(),
        ];
        v.extend(rest.iter().map(|s| (*s).to_owned()));
        v
    }

    #[test]
    fn rfmix_columns_are_chrom_pos_cm_and_blank_lines_go() {
        assert_eq!(
            rfmix_rows(MAP).unwrap(),
            "chr20\t60000\t0.0\nchr20\t70001\t0.12\nchr20\t80000\t1.5\n"
        );
        assert_eq!(
            rfmix_rows("a b c d\r\ne f g h").unwrap(),
            "a\td\tc\ne\th\tg\n"
        );
        assert_eq!(splitlines("a\r\nb\rc\n\nd"), ["a", "b", "c", "", "d"]);
    }

    #[test]
    fn writes_the_member_byte_for_byte_and_the_rfmix_file() {
        let d = scratch("ok");
        let z = zip_with(
            &d,
            &[(
                "chr_in_chrom_field/plink.chrchr20.GRCh38.map",
                MAP.as_bytes(),
                true,
            )],
        );
        assert_eq!(cli(&args(&z, &d, &["--rfmix", "chr20"])), 0);
        assert_eq!(fs::read_to_string(d.join("chr20.map")).unwrap(), MAP);
        assert_eq!(
            fs::read_to_string(d.join("chr20.rfmix.map")).unwrap(),
            "chr20\t60000\t0.0\nchr20\t70001\t0.12\nchr20\t80000\t1.5\n"
        );
        assert_eq!(
            fs::read_dir(&d).unwrap().count(),
            3,
            "no temporary file left"
        );
    }

    #[test]
    fn a_truncated_row_writes_nothing_for_that_chromosome() {
        let d = scratch("short");
        let z = zip_with(
            &d,
            &[(
                "chr_in_chrom_field/plink.chrchr21.GRCh38.map",
                b"chr21 . 0.0 100\nchr21 . 0.5",
                false,
            )],
        );
        fs::write(d.join("chr21.rfmix.map"), "old\n").unwrap();
        assert_eq!(cli(&args(&z, &d, &["--rfmix", "chr21"])), 2);
        assert!(!d.join("chr21.map").exists());
        assert_eq!(
            fs::read_to_string(d.join("chr21.rfmix.map")).unwrap(),
            "old\n"
        );
    }

    #[test]
    fn a_missing_chromosome_exits_2() {
        let d = scratch("missing");
        let z = zip_with(
            &d,
            &[("chr_in_chrom_field/plink.chr1.GRCh38.map", b"x", false)],
        );
        assert_eq!(cli(&args(&z, &d, &["chrX"])), 2);
    }
}
