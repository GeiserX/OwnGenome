//! `append-n SRC OUT CONTIG LENGTH`: a copy of the reference SRC with one
//! more contig of LENGTH `N`s, 60 a line. Step 35's Parascopy cases need a
//! chr11 the fixture reference lacks, because a model names it.

use super::{Fail, create, exactly, n_lines, number};
use std::fs::File;
use std::io::{self, Write};
use std::path::Path;

pub fn run(args: &[String]) -> Result<(), Fail> {
    let a = exactly(args, 4, "SRC OUT CONTIG LENGTH")?;
    let length = number(&a[3], "LENGTH")?;
    let src = Path::new(&a[0]);
    let mut input = File::open(src).map_err(|e| Fail::Run(format!("{}: {e}", src.display())))?;
    let mut out = create(Path::new(&a[1]))?;
    io::copy(&mut input, &mut out)?;
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
}
