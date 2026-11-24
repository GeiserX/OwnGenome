//! `og dev build roh-af`: the allele-frequency panel of step 11 (runs of
//! homozygosity), built from gnomAD's genome sites VCFs.
//!
//! `bcftools roh` needs the frequency of the alternate allele at each site
//! to tell a homozygous call that is expected anyway (a common allele) from
//! one that points to a run. Without a panel step 11 can only give it one
//! flat value (`--AF-dflt 0.4`) for every site. The panel gives it gnomAD's.
//!
//! ```text
//! og dev build roh-af rows --vcf FILE [--out FILE] [--threads N]
//! og dev build roh-af panel --rows FILE... --out FILE.tsv.gz [--check]
//! ```
//!
//! `rows` reads one gnomAD sites VCF (plain, gzip or BGZF; `-` is stdin, so
//! a download can stream in without being kept) and writes the panel rows
//! it holds, without a header. A site is kept when it is on chr1 to chr22,
//! its FILTER is `PASS`, REF and ALT are one base each (A, C, G or T), the
//! locus has one alternate allele in the whole callset (`n_alt_alleles=1`:
//! biallelic) and its global `AF` is 0.01 to 0.99. Each row is
//! `CHROM POS REF,ALT AF AF_afr AF_amr AF_eas AF_nfe AF_sas`, tab separated:
//! the global frequency, then the five continental groups (gnomAD's labels;
//! `nfe` is the European group of step 26's 1000 Genomes `EUR`). Each
//! frequency is copied from the VCF as written, or `.` when the VCF has
//! none (no allele number in that group). Columns 1 to 4 are the
//! `--AF-file` format of `bcftools roh`; it reads column 4 and nothing
//! after it.
//!
//! `panel` joins the rows of the chromosomes, in order, under the panel's
//! header and writes them as BGZF (`og_io::bgzf::Writer`, laid out as
//! `bgzip` lays it out, the bytes a function of the input only). It refuses
//! rows out of chromosome order, positions that do not go up, and a row
//! with a wrong number of columns. With `--check` it writes nothing and
//! exits 1 when `--out` differs from what it would write. The index is
//! `tabix -s1 -b2 -e2` of `BCFTOOLS_IMAGE`, made by `scripts/setup.sh`
//! after the download.
//!
//! How the published panel was made, and how to make it again, byte for
//! byte: `docs/00-reference-setup.md` ("ROH allele-frequency panel").

use super::{BuildError, Failure, chrom_index, input_error};
use std::fs;
use std::io::{BufRead, Write};
use std::num::NonZero;
use std::path::Path;

/// The gnomAD release the panel is built from (`ROH_AF_PANEL_URL` in
/// `versions.env` names the same one).
pub const GNOMAD_RELEASE: &str = "4.1";

/// The INFO keys of the five continental groups, in column order.
pub const GROUPS: [&str; 5] = ["AF_afr", "AF_amr", "AF_eas", "AF_nfe", "AF_sas"];

/// The global frequency range a site must fall in, both ends included.
const AF_MIN: f64 = 0.01;
const AF_MAX: f64 = 0.99;

/// Columns of a panel row.
const COLUMNS: usize = 4 + GROUPS.len();

/// The panel's header: what it is, where it comes from and its columns.
pub fn header() -> String {
    format!(
        "# OwnGenome ROH allele-frequency panel (step 11, bcftools roh --AF-file): biallelic PASS SNVs of chr1 to chr22 with a global AF of 0.01 to 0.99\n\
         # Source: gnomAD v{GNOMAD_RELEASE} genomes, sites VCFs gnomad.genomes.v{GNOMAD_RELEASE}.sites.chr1-22.vcf.bgz, released under CC0 1.0; built by og dev build roh-af\n\
         # AF is the global alternate allele frequency, AF_<group> that of a gnomAD genetic ancestry group (. where the group has no allele number)\n\
         #CHROM\tPOS\tREF,ALT\tAF\t{}\n",
        GROUPS.join("\t")
    )
}

/// One base, as the panel keeps them.
fn base(b: &[u8]) -> bool {
    matches!(b, [b'A' | b'C' | b'G' | b'T'])
}

/// The value of `key` in an INFO column, or `None` when it is absent.
fn info<'a>(info: &'a [u8], key: &str) -> Option<&'a [u8]> {
    info.split(|&c| c == b';').find_map(|kv| {
        kv.strip_prefix(key.as_bytes())
            .and_then(|rest| rest.strip_prefix(b"="))
    })
}

/// The panel row of one VCF record, or `None` when the site is not kept.
/// Reads the cheap columns first, so most records cost a few comparisons.
pub fn row(line: &[u8]) -> Result<Option<String>, String> {
    if line.first() == Some(&b'#') || line.is_empty() {
        return Ok(None);
    }
    let mut cols = line.splitn(9, |&c| c == b'\t');
    let mut next = || cols.next();
    let (Some(chrom), Some(pos), Some(_id), Some(r), Some(a), Some(_qual), Some(filter), Some(inf)) = (
        next(),
        next(),
        next(),
        next(),
        next(),
        next(),
        next(),
        next(),
    ) else {
        return Err(format!(
            "a VCF record with fewer than 8 columns: {}",
            String::from_utf8_lossy(&line[..line.len().min(80)])
        ));
    };
    if chrom_index(chrom).is_none_or(|i| i >= 22) || filter != b"PASS" || !base(r) || !base(a) {
        return Ok(None);
    }
    let Some(af) = info(inf, "AF") else {
        return Ok(None);
    };
    let text = std::str::from_utf8(af).map_err(|_| "an AF that is not text".to_owned())?;
    let value: f64 = text
        .parse()
        .map_err(|_| format!("an AF that is not a number: {text:?}"))?;
    if !(AF_MIN..=AF_MAX).contains(&value) {
        return Ok(None);
    }
    if info(inf, "n_alt_alleles") != Some(b"1") {
        return Ok(None);
    }
    let text_of = |b: &[u8]| String::from_utf8_lossy(b).into_owned();
    let mut out = format!(
        "{}\t{}\t{},{}\t{text}",
        text_of(chrom),
        text_of(pos),
        text_of(r),
        text_of(a)
    );
    for g in GROUPS {
        out.push('\t');
        out.push_str(&info(inf, g).map_or_else(|| ".".to_owned(), text_of));
    }
    out.push('\n');
    Ok(Some(out))
}

/// `rows`: the panel rows of one VCF, written to `out`. Returns the
/// number of records read and of rows written. An input with no `#CHROM`
/// line is no VCF (an empty download, an error page) and fails, so it
/// never becomes an empty part that `panel` would join without a word.
fn rows(
    input: Box<dyn BufRead + Send>,
    out: &mut dyn Write,
    what: &Path,
) -> Result<(u64, u64), BuildError> {
    let mut input = input;
    let (mut read, mut kept) = (0u64, 0u64);
    let mut header = false;
    let mut buf = Vec::with_capacity(1 << 14);
    loop {
        buf.clear();
        if input
            .read_until(b'\n', &mut buf)
            .map_err(|e| input_error(what, &e))?
            == 0
        {
            if !header {
                return Err(BuildError::Input(format!(
                    "{}: no #CHROM header line, so not a VCF (an empty or cut download?)",
                    what.display()
                )));
            }
            return Ok((read, kept));
        }
        let mut line = buf.as_slice();
        if let Some(l) = line.strip_suffix(b"\n") {
            line = l.strip_suffix(b"\r").unwrap_or(l);
        }
        header |= line.starts_with(b"#CHROM\t");
        if line.first() != Some(&b'#') {
            read += 1;
        }
        if let Some(r) =
            row(line).map_err(|m| BuildError::Input(format!("{}: {m}", what.display())))?
        {
            kept += 1;
            out.write_all(r.as_bytes())
                .map_err(|e| input_error(Path::new("output"), &e))?;
        }
    }
}

/// `panel`: the header and every row of `parts`, checked for order.
pub fn panel(parts: &[&Path]) -> Result<Vec<u8>, BuildError> {
    let mut text = header().into_bytes();
    let mut last: Option<(usize, u64)> = None;
    for part in parts {
        let bytes = fs::read(part).map_err(|e| input_error(part, &e))?;
        for (n, line) in bytes.split(|&c| c == b'\n').enumerate() {
            if line.is_empty() {
                continue;
            }
            let cols: Vec<&[u8]> = line.split(|&c| c == b'\t').collect();
            let place = || format!("{}:{}", part.display(), n + 1);
            if cols.len() != COLUMNS {
                return Err(BuildError::Data(format!(
                    "{}: {} columns, not {COLUMNS}",
                    place(),
                    cols.len()
                )));
            }
            let chrom = chrom_index(cols[0])
                .filter(|&i| i < 22)
                .ok_or_else(|| BuildError::Data(format!("{}: not chr1 to chr22", place())))?;
            let pos: u64 = std::str::from_utf8(cols[1])
                .ok()
                .and_then(|p| p.parse().ok())
                .ok_or_else(|| BuildError::Data(format!("{}: POS is not a number", place())))?;
            if last.is_some_and(|l| (chrom, pos) <= l) {
                return Err(BuildError::Data(format!(
                    "{}: out of order (chromosomes chr1 to chr22, positions going up, one row per position)",
                    place()
                )));
            }
            last = Some((chrom, pos));
            text.extend_from_slice(line);
            text.push(b'\n');
        }
    }
    let mut w = og_io::bgzf::Writer::new(Vec::new(), og_io::bgzf::workers())
        .map_err(|e| input_error(Path::new("BGZF"), &e))?;
    w.write_all(&text)
        .map_err(|e| input_error(Path::new("BGZF"), &e))?;
    w.finish().map_err(|e| input_error(Path::new("BGZF"), &e))
}

const USAGE: &str = "usage: og dev build roh-af rows --vcf FILE [--out FILE] [--threads N]\n       og dev build roh-af panel --rows FILE... --out FILE.tsv.gz [--check]\n";

/// `og dev build roh-af ARGS`; returns the exit code.
pub fn cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let result = match args.first().map(String::as_str) {
        Some("rows") => cli_rows(&args[1..], out, err),
        Some("panel") => cli_panel(&args[1..], err),
        Some("-h" | "--help") => {
            let _ = out.write_all(USAGE.as_bytes());
            return 0;
        }
        Some(other) => Err(Failure::Usage(format!("unknown mode '{other}'"))),
        None => Err(Failure::Usage(
            "the following arguments are required: rows or panel".to_owned(),
        )),
    };
    match result {
        Ok(code) => code,
        Err(Failure::Usage(m)) => {
            let _ = writeln!(err, "{USAGE}og dev build roh-af: error: {m}");
            2
        }
        Err(Failure::Build(e)) => {
            let _ = writeln!(err, "ERROR: {e}");
            e.code()
        }
    }
}

fn cli_rows(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> Result<u8, Failure> {
    let a = super::Args::parse(args, &["--vcf", "--out", "--threads"], &[])?;
    let vcf = Path::new(a.need("--vcf")?);
    let threads = match a.get("--threads") {
        None => NonZero::new(2).expect("2 is not 0"),
        Some(t) => t
            .parse::<usize>()
            .ok()
            .and_then(NonZero::new)
            .ok_or_else(|| format!("argument --threads: not a whole number above 0: {t}"))?,
    };
    let input = og_io::bgzf::open(vcf, threads).map_err(|e| input_error(vcf, &e))?;
    let (read, kept) = match a.get("--out") {
        None => rows(input, out, vcf)?,
        Some(path) => {
            let path = Path::new(path);
            let mut f = og_io::AtomicFile::create(path).map_err(|e| input_error(path, &e))?;
            let counts = rows(input, &mut f, vcf)?;
            f.commit().map_err(|e| input_error(path, &e))?;
            counts
        }
    };
    let _ = writeln!(
        err,
        "{}: {read} records read, {kept} panel rows",
        vcf.display()
    );
    Ok(0)
}

fn cli_panel(args: &[String], err: &mut dyn Write) -> Result<u8, Failure> {
    // --rows takes every following word that is not a flag.
    let mut parts: Vec<&Path> = Vec::new();
    let mut rest: Vec<String> = Vec::new();
    let mut it = args.iter().peekable();
    while let Some(arg) = it.next() {
        if arg == "--rows" {
            while let Some(p) = it.next_if(|p| !p.starts_with("--")) {
                parts.push(Path::new(p));
            }
        } else {
            rest.push(arg.clone());
        }
    }
    let a = super::Args::parse(&rest, &["--out"], &["--check"])?;
    let out = Path::new(a.need("--out")?);
    if parts.is_empty() {
        return Err(Failure::Usage(
            "the following arguments are required: --rows".to_owned(),
        ));
    }
    let bytes = panel(&parts)?;
    if a.has("--check") {
        let have = fs::read(out).map_err(|e| input_error(out, &e))?;
        if have == bytes {
            let _ = writeln!(err, "{}: as built", out.display());
            return Ok(0);
        }
        let _ = writeln!(
            err,
            "{}: differs from what og dev build roh-af panel writes",
            out.display()
        );
        return Ok(1);
    }
    super::write_atomic(out, &bytes)?;
    let _ = writeln!(err, "{}: written", out.display());
    Ok(0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;

    const INFO_ROW: &str = "AC=10;AN=100;AF=0.100000;AF_afr=0.200000;AF_amr=0.0500000;AF_eas=0.00000;AF_nfe=0.120000;n_alt_alleles=1;vep=A|b";

    fn rec(chrom: &str, pos: u32, r: &str, a: &str, filter: &str, info: &str) -> String {
        format!("{chrom}\t{pos}\t.\t{r}\t{a}\t.\t{filter}\t{info}")
    }

    fn s(v: &[&str]) -> Vec<String> {
        v.iter().map(|x| (*x).to_owned()).collect()
    }

    #[test]
    fn a_common_biallelic_snv_is_kept_with_its_groups() {
        let r = row(rec("chr22", 100, "A", "G", "PASS", INFO_ROW).as_bytes()).unwrap();
        // AF_sas is absent: the group had no allele number.
        assert_eq!(
            r.as_deref(),
            Some("chr22\t100\tA,G\t0.100000\t0.200000\t0.0500000\t0.00000\t0.120000\t.\n")
        );
    }

    #[test]
    fn every_other_site_is_left_out() {
        let multi = INFO_ROW.replace("n_alt_alleles=1", "n_alt_alleles=2");
        let rare = INFO_ROW.replace("AF=0.100000", "AF=0.00999");
        let fixed = INFO_ROW.replace("AF=0.100000", "AF=0.991");
        let no_count = INFO_ROW.replace(";n_alt_alleles=1", "");
        for (what, line) in [
            ("filtered", rec("chr1", 5, "A", "G", "AS_VQSR", INFO_ROW)),
            ("chrX", rec("chrX", 5, "A", "G", "PASS", INFO_ROW)),
            (
                "alt contig",
                rec("chr1_KI270706v1_random", 5, "A", "G", "PASS", INFO_ROW),
            ),
            ("indel", rec("chr1", 5, "AT", "A", "PASS", INFO_ROW)),
            ("N", rec("chr1", 5, "N", "A", "PASS", INFO_ROW)),
            (
                "multiallelic locus",
                rec("chr1", 5, "A", "G", "PASS", &multi),
            ),
            (
                "no n_alt_alleles",
                rec("chr1", 5, "A", "G", "PASS", &no_count),
            ),
            ("rare", rec("chr1", 5, "A", "G", "PASS", &rare)),
            ("near fixed", rec("chr1", 5, "A", "G", "PASS", &fixed)),
            ("header", "#CHROM\tPOS".to_owned()),
        ] {
            assert_eq!(row(line.as_bytes()).unwrap(), None, "{what}");
        }
        // The ends of the range are kept.
        for af in ["0.01", "0.99", "1.00000e-02"] {
            let i = INFO_ROW.replace("AF=0.100000", &format!("AF={af}"));
            assert!(
                row(rec("chr1", 5, "A", "G", "PASS", &i).as_bytes())
                    .unwrap()
                    .is_some(),
                "{af}"
            );
        }
    }

    #[test]
    fn a_key_is_matched_whole() {
        // AF_afr must not be read as AF, nor AF_afr_XX as AF_afr.
        let i = "AF_afr=0.5;AF_afr_XX=0.7;AF=0.3;n_alt_alleles=1";
        let r = row(rec("chr2", 7, "C", "T", "PASS", i).as_bytes())
            .unwrap()
            .unwrap();
        assert_eq!(r, "chr2\t7\tC,T\t0.3\t0.5\t.\t.\t.\t.\n");
        assert!(row(rec("chr2", 7, "C", "T", "PASS", "AF=x;n_alt_alleles=1").as_bytes()).is_err());
        assert!(row(b"chr2\t7\t.\tC").is_err());
    }

    fn tmp(name: &str) -> std::path::PathBuf {
        let d = std::env::temp_dir().join(format!("og-roh-af-{}-{name}", std::process::id()));
        fs::create_dir_all(&d).unwrap();
        d
    }

    #[test]
    fn rows_then_panel_round_trip_and_check() {
        let d = tmp("trip");
        let vcf = d.join("in.vcf");
        let body = [
            "##fileformat=VCFv4.2".to_owned(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_owned(),
            rec("chr1", 10, "A", "G", "PASS", INFO_ROW),
            rec("chr1", 20, "A", "G", "AC0", INFO_ROW),
            rec("chr2", 5, "G", "C", "PASS", INFO_ROW),
        ]
        .join("\n");
        fs::write(&vcf, body + "\n").unwrap();
        let part = d.join("part.tsv");
        let (mut o, mut e) = (Vec::new(), Vec::new());
        let argv = s(&[
            "rows",
            "--vcf",
            vcf.to_str().unwrap(),
            "--out",
            part.to_str().unwrap(),
        ]);
        assert_eq!(
            cli(&argv, &mut o, &mut e),
            0,
            "{}",
            String::from_utf8_lossy(&e)
        );
        assert!(String::from_utf8_lossy(&e).contains("3 records read, 2 panel rows"));
        let out = d.join("panel.tsv.gz");
        let argv = s(&[
            "panel",
            "--rows",
            part.to_str().unwrap(),
            "--out",
            out.to_str().unwrap(),
        ]);
        assert_eq!(cli(&argv, &mut o, &mut Vec::new()), 0);
        let mut text = String::new();
        og_io::open(&out)
            .unwrap()
            .read_to_string(&mut text)
            .unwrap();
        assert!(text.starts_with(&header()));
        assert_eq!(text.lines().filter(|l| !l.starts_with('#')).count(), 2);
        // The BGZF ends with the end-of-file block, so tabix takes it.
        assert!(fs::read(&out).unwrap().ends_with(&og_io::bgzf::EOF_BLOCK));
        let check = s(&[
            "panel",
            "--rows",
            part.to_str().unwrap(),
            "--out",
            out.to_str().unwrap(),
            "--check",
        ]);
        assert_eq!(cli(&check, &mut o, &mut Vec::new()), 0);
        // One changed byte in a row: --check exits 1 and writes nothing.
        let before = fs::read(&out).unwrap();
        fs::write(
            &part,
            fs::read_to_string(&part)
                .unwrap()
                .replace("chr2\t5", "chr2\t6"),
        )
        .unwrap();
        assert_eq!(cli(&check, &mut o, &mut Vec::new()), 1);
        assert_eq!(fs::read(&out).unwrap(), before);
        fs::remove_dir_all(&d).unwrap();
    }

    #[test]
    fn rows_refuses_input_with_no_vcf_header() {
        let d = tmp("empty");
        let part = d.join("part.tsv");
        for (name, body) in [
            ("empty", String::new()),
            (
                "no-chrom-line",
                rec("chr1", 10, "A", "G", "PASS", INFO_ROW) + "\n",
            ),
        ] {
            let vcf = d.join(format!("{name}.vcf"));
            fs::write(&vcf, body).unwrap();
            let (mut o, mut e) = (Vec::new(), Vec::new());
            let argv = s(&[
                "rows",
                "--vcf",
                vcf.to_str().unwrap(),
                "--out",
                part.to_str().unwrap(),
            ]);
            assert_ne!(cli(&argv, &mut o, &mut e), 0, "{name}");
            assert!(
                String::from_utf8_lossy(&e).contains("no #CHROM header line"),
                "{name}: {}",
                String::from_utf8_lossy(&e)
            );
            assert!(!part.exists(), "{name}: no part is written");
        }
        fs::remove_dir_all(&d).unwrap();
    }

    #[test]
    fn panel_refuses_rows_out_of_order_or_misshapen() {
        let d = tmp("order");
        let row1 = "chr2\t5\tA,G\t0.1\t.\t.\t.\t.\t.\n";
        let row0 = "chr1\t9\tA,G\t0.1\t.\t.\t.\t.\t.\n";
        for (what, body) in [
            ("chromosomes back", format!("{row1}{row0}")),
            ("same position", format!("{row0}{row0}")),
            (
                "chr10 before chr9",
                "chr10\t1\tA,G\t0.1\t.\t.\t.\t.\t.\nchr9\t1\tA,G\t0.1\t.\t.\t.\t.\t.\n".to_owned(),
            ),
            ("columns", "chr1\t9\tA,G\t0.1\n".to_owned()),
            ("chrX", "chrX\t9\tA,G\t0.1\t.\t.\t.\t.\t.\n".to_owned()),
        ] {
            let p = d.join("p.tsv");
            fs::write(&p, body).unwrap();
            assert!(matches!(panel(&[&p]), Err(BuildError::Data(_))), "{what}");
        }
        // chr9 then chr10 across two parts is in order.
        let (a, b) = (d.join("a.tsv"), d.join("b.tsv"));
        fs::write(&a, "chr9\t1\tA,G\t0.1\t.\t.\t.\t.\t.\n").unwrap();
        fs::write(&b, "chr10\t1\tA,G\t0.1\t.\t.\t.\t.\t.\n").unwrap();
        assert!(panel(&[&a, &b]).is_ok());
        fs::remove_dir_all(&d).unwrap();
    }

    #[test]
    fn usage_errors_exit_2() {
        for argv in [
            s(&[]),
            s(&["nope"]),
            s(&["rows"]),
            s(&["panel", "--out", "x"]),
            s(&["rows", "--vcf", "x", "--threads", "0"]),
        ] {
            assert_eq!(cli(&argv, &mut Vec::new(), &mut Vec::new()), 2, "{argv:?}");
        }
    }
}
