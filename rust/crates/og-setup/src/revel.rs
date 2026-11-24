//! `og setup revel ZIP TSV` projects the REVEL release CSV onto GRCh38.
//! Stream the ZIP member through the projection and use sort's external
//! merge so the full table does not have to fit in memory.

use crate::zip::{Archive, ZipError};
use og_io::AtomicFile;
use std::io::{self, BufWriter, Write};
use std::path::Path;
use std::process::{Command, Stdio};

const MEMBER: &str = "revel_with_transcript_ids";
const HEADER: &[u8] = b"#chr\tpos\tref\talt\tREVEL\n";

struct Projection<W> {
    out: W,
    pending: Vec<u8>,
    header: bool,
    records: usize,
}

impl<W: Write> Projection<W> {
    fn line(&mut self, line: &[u8]) -> io::Result<()> {
        let line = line.strip_suffix(b"\r").unwrap_or(line);
        if self.header {
            self.header = false;
            if line.starts_with(b"chr,hg19_pos,grch38_pos,ref,alt,aaref,aaalt,REVEL,") {
                return Ok(());
            }
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "unexpected REVEL CSV header",
            ));
        }
        let mut fields = line.split(|&b| b == b',');
        let [
            Some(chrom),
            Some(_),
            Some(pos),
            Some(reference),
            Some(alt),
            Some(_),
            Some(_),
            Some(score),
            Some(_),
        ] = std::array::from_fn(|_| fields.next())
        else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "short REVEL CSV row",
            ));
        };
        if pos == b"." {
            return Ok(());
        }
        self.out.write_all(b"chr")?;
        for field in [chrom, pos, reference, alt] {
            self.out.write_all(field)?;
            self.out.write_all(b"\t")?;
        }
        self.out.write_all(score)?;
        self.out.write_all(b"\n")?;
        self.records += 1;
        Ok(())
    }
}

impl<W: Write> Write for Projection<W> {
    fn write(&mut self, mut bytes: &[u8]) -> io::Result<usize> {
        let len = bytes.len();
        while let Some(end) = bytes.iter().position(|&b| b == b'\n') {
            let line = &bytes[..end];
            if self.pending.is_empty() {
                self.line(line)?;
            } else {
                let mut pending = std::mem::take(&mut self.pending);
                pending.extend_from_slice(line);
                self.line(&pending)?;
                pending.clear();
                self.pending = pending;
            }
            bytes = &bytes[end + 1..];
        }
        self.pending.extend_from_slice(bytes);
        Ok(len)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.out.flush()
    }
}

fn project(
    mut archive: Archive,
    entry: &crate::zip::Entry,
    out: impl Write,
) -> Result<(), ZipError> {
    let mut projection = Projection {
        out: BufWriter::new(out),
        pending: Vec::new(),
        header: true,
        records: 0,
    };
    archive.copy_to(entry, &mut projection)?;
    let pending = std::mem::take(&mut projection.pending);
    if !pending.is_empty() {
        projection.line(&pending)?;
    }
    if projection.records == 0 {
        return Err(ZipError::Bad("no GRCh38 REVEL records".into()));
    }
    projection.flush()?;
    Ok(())
}

fn convert(zip: &str, tsv: &Path) -> Result<(), ZipError> {
    let archive = Archive::open(zip)?;
    let entry = archive
        .find(MEMBER)
        .cloned()
        .ok_or_else(|| ZipError::Bad(format!("ZIP does not contain {MEMBER}")))?;
    let mut output = AtomicFile::create(tsv)?;
    output.write_all(HEADER)?;
    let dir = tsv
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    let mut sort = Command::new("sort")
        .env("LC_ALL", "C")
        .arg("-T")
        .arg(dir)
        .args(["-S", "256M", "-k1,1V", "-k2,2n"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()?;
    let input = sort.stdin.take().expect("piped stdin");
    let mut sorted = sort.stdout.take().expect("piped stdout");
    let (projected, copied) = std::thread::scope(|scope| {
        let writer = scope.spawn(|| project(archive, &entry, input));
        let copied = io::copy(&mut sorted, &mut output);
        // A write failure must close the read end before joining the producer.
        drop(sorted);
        (writer.join().expect("REVEL projection thread"), copied)
    });
    let status = sort.wait()?;
    copied?;
    projected?;
    if !status.success() {
        return Err(ZipError::Bad(format!("sort exited with {status}")));
    }
    output.commit()?;
    Ok(())
}

/// Write the sorted GRCh38 table to TSV, or report why the input could not be converted.
pub fn cli(args: &[String], _out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let [zip, tsv] = args else {
        let _ = writeln!(err, "usage: og setup revel ZIP TSV");
        return 2;
    };
    match convert(zip, Path::new(tsv)) {
        Ok(()) => 0,
        Err(error) => {
            let _ = writeln!(err, "og setup revel: {zip}: {error}");
            1
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::zip::tests::build;
    use std::fs;

    // First 20 records of REVEL v1.3 (Zenodo 7072866, ODbL 1.0).
    // EXPECTED was generated once by the former setup guide's awk and sort recipe.
    const SLICE: &[u8] = br#"chr,hg19_pos,grch38_pos,ref,alt,aaref,aaalt,REVEL,Ensembl_transcriptid
1,35142,35142,G,A,T,M,0.027,ENST00000417324
1,35142,35142,G,C,T,R,0.035,ENST00000417324
1,35142,35142,G,T,T,K,0.043,ENST00000417324
1,35143,35143,T,A,T,S,0.018,ENST00000417324
1,35143,35143,T,C,T,A,0.034,ENST00000417324
1,35143,35143,T,G,T,P,0.039,ENST00000417324
1,35144,35144,A,C,C,W,0.012,ENST00000417324
1,35145,35145,C,A,C,F,0.023,ENST00000417324
1,35145,35145,C,G,C,S,0.029,ENST00000417324
1,35145,35145,C,T,C,Y,0.016,ENST00000417324
1,35146,35146,A,C,C,G,0.031,ENST00000417324
1,35146,35146,A,G,C,R,0.016,ENST00000417324
1,35146,35146,A,T,C,S,0.025,ENST00000417324
1,35147,35147,T,A,L,F,0.004,ENST00000417324
1,35147,35147,T,G,L,F,0.004,ENST00000417324
1,35148,35148,A,G,L,S,0.010,ENST00000417324
1,35149,35149,A,C,L,V,0.029,ENST00000417324
1,35149,35149,A,T,L,I,0.022,ENST00000417324
1,35150,35150,T,A,R,S,0.038,ENST00000417324
1,35150,35150,T,G,R,S,0.055,ENST00000417324
"#;
    const EXPECTED: &[u8] = br#"#chr	pos	ref	alt	REVEL
chr1	35142	G	A	0.027
chr1	35142	G	C	0.035
chr1	35142	G	T	0.043
chr1	35143	T	A	0.018
chr1	35143	T	C	0.034
chr1	35143	T	G	0.039
chr1	35144	A	C	0.012
chr1	35145	C	A	0.023
chr1	35145	C	G	0.029
chr1	35145	C	T	0.016
chr1	35146	A	C	0.031
chr1	35146	A	G	0.016
chr1	35146	A	T	0.025
chr1	35147	T	A	0.004
chr1	35147	T	G	0.004
chr1	35148	A	G	0.010
chr1	35149	A	C	0.029
chr1	35149	A	T	0.022
chr1	35150	T	A	0.038
chr1	35150	T	G	0.055
"#;

    #[test]
    fn matches_the_setup_guides_revel_slice() {
        let dir = std::env::temp_dir().join(format!("og-revel-slice-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let zip = dir.join("scores.zip");
        let tsv = dir.join("scores.tsv");
        fs::write(&zip, build(&[(MEMBER, SLICE, true)], false)).unwrap();
        let args = [zip.display().to_string(), tsv.display().to_string()];
        let mut error = Vec::new();
        assert_eq!(
            cli(&args, &mut io::sink(), &mut error),
            0,
            "{}",
            String::from_utf8_lossy(&error)
        );
        assert_eq!(fs::read(&tsv).unwrap(), EXPECTED);

        // A malformed replacement fails and leaves the last complete table.
        fs::write(&zip, build(&[(MEMBER, b"wrong header\n", true)], false)).unwrap();
        assert_ne!(cli(&args, &mut io::sink(), &mut error), 0);
        assert_eq!(fs::read(&tsv).unwrap(), EXPECTED);
        fs::remove_dir_all(dir).unwrap();
    }
}
