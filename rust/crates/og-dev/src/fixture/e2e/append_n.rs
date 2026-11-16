//! `append-n SRC OUT CONTIG LENGTH`: a copy of the reference SRC with one
//! more contig of LENGTH `N`s, 60 a line. Step 35's Parascopy cases need a
//! chr11 the fixture reference lacks, because a model names it.

use super::{Fail, create, exactly, n_lines, number};
use std::fs::{self, File};
use std::io::{Read, Write};
use std::os::unix::fs::MetadataExt;
use std::path::Path;

pub fn run(args: &[String]) -> Result<(), Fail> {
    let a = exactly(args, 4, "SRC OUT CONTIG LENGTH")?;
    let length = number(&a[3], "LENGTH")?;
    let src = Path::new(&a[0]);
    let mut input = File::open(src).map_err(|e| Fail::Run(format!("{}: {e}", src.display())))?;
    // OUT is truncated first: when it is SRC (or a link to it), the copy
    // would read nothing.
    let same = |o: fs::Metadata| {
        input
            .metadata()
            .is_ok_and(|i| (i.dev(), i.ino()) == (o.dev(), o.ino()))
    };
    if fs::metadata(&a[1]).is_ok_and(same) {
        return Err(Fail::Usage("append-n SRC and OUT are the same file".into()));
    }
    let mut out = create(Path::new(&a[1]))?;
    let mut buf = vec![0; 1 << 16];
    let mut last = b'\n';
    loop {
        let n = input.read(&mut buf)?;
        if n == 0 {
            break;
        }
        out.write_all(&buf[..n])?;
        last = buf[n - 1];
    }
    // A SRC whose last line has no newline would glue the header onto it.
    if last != b'\n' {
        out.write_all(b"\n")?;
    }
    writeln!(out, ">{}", a[2])?;
    if n_lines(&mut out, length, 0)? != 0 {
        out.write_all(b"\n")?;
    }
    out.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use crate::fixture::e2e::tests::{run, scratch};
    use std::fs;

    #[test]
    fn copies_the_reference_then_the_n_contig() {
        let d = scratch("append-n");
        fs::write(d.join("ref.fa"), ">chr1\nACGT\n").unwrap();
        let out = d.join("out.fa");
        let (code, _, err) = run(
            "append-n",
            &[
                d.join("ref.fa").to_str().unwrap(),
                out.to_str().unwrap(),
                "chr11",
                "125",
            ],
        );
        assert_eq!(code, Some(0), "{err}");
        let want = format!(">chr1\nACGT\n>chr11\n{n}\n{n}\nNNNNN\n", n = "N".repeat(60));
        assert_eq!(fs::read_to_string(&out).unwrap(), want);
        // 120 is two whole lines: no empty line after them.
        run(
            "append-n",
            &[
                d.join("ref.fa").to_str().unwrap(),
                out.to_str().unwrap(),
                "chr11",
                "120",
            ],
        );
        assert!(fs::read_to_string(&out).unwrap().ends_with("NNN\n"));
        assert!(!fs::read_to_string(&out).unwrap().ends_with("\n\n"));
        fs::remove_dir_all(d).unwrap();
    }

    #[test]
    fn a_last_line_without_newline_gets_one_before_the_header() {
        let d = scratch("append-n-nl");
        let out = d.join("out.fa");
        for (src, want) in [
            ("", ">chrN\nNN\n"),
            (">chr1\nACGT", ">chr1\nACGT\n>chrN\nNN\n"),
            (">chr1\nACGT\n", ">chr1\nACGT\n>chrN\nNN\n"),
        ] {
            fs::write(d.join("ref.fa"), src).unwrap();
            let (code, _, err) = run(
                "append-n",
                &[
                    d.join("ref.fa").to_str().unwrap(),
                    out.to_str().unwrap(),
                    "chrN",
                    "2",
                ],
            );
            assert_eq!(code, Some(0), "{err}");
            assert_eq!(fs::read_to_string(&out).unwrap(), want, "{src:?}");
        }
        fs::remove_dir_all(d).unwrap();
    }

    #[test]
    fn src_as_out_is_refused_and_left_whole() {
        let d = scratch("append-n-same");
        let src = d.join("ref.fa");
        fs::write(&src, ">chr1\nACGT\n").unwrap();
        std::os::unix::fs::symlink(&src, d.join("link.fa")).unwrap();
        for out in [src.clone(), d.join("link.fa"), d.join(".").join("ref.fa")] {
            let (code, _, err) = run(
                "append-n",
                &[src.to_str().unwrap(), out.to_str().unwrap(), "chrN", "2"],
            );
            assert_eq!(code, Some(2), "{err}");
            assert!(err.contains("same file"), "{err}");
            assert_eq!(fs::read_to_string(&src).unwrap(), ">chr1\nACGT\n");
        }
        fs::remove_dir_all(d).unwrap();
    }
}
