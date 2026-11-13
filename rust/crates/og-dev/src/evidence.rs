//! `og dev evidence`: the measurements of `.github/workflows/public-evidence.yml`
//! on public data (GIAB's HG002), which were Python inside the workflow.
//!
//! - `og dev evidence clashes RAW.vcf.gz SPLIT.vcf.gz DATA_DIR OUT.tsv`: the
//!   positions with two records or more, one of them an SNV, in the VCF as
//!   GIAB writes it and split with `bcftools norm -m-`, on the chromosomes
//!   of the archaic tables in DATA_DIR (`chr*.tsv.gz`), and how many rows of
//!   those tables sit at such a position. OUT holds `raw_positions`,
//!   `raw_on_table`, `split_positions` and `split_on_table`.
//! - `og dev evidence hg002-gvcf TRUTH.vcf.gz TRUTH.bed REF.fasta SITES.tsv OUT.vcf`:
//!   a gVCF of GIAB's genotypes plus a one-base reference block at every
//!   panel site inside the benchmark regions that no GIAB record covers,
//!   with the reference base from the FASTA (through its `.fai`). The
//!   output is unsorted: the workflow sorts it with bcftools.
//!
//! Both stream their VCFs. hg002-gvcf keeps, per chromosome, the merged
//! spans the GIAB records cover instead of every covered position, and
//! reads one reference base per site with a positioned read.
//!
//! Exit 0 when written, 1 when an input cannot be read or parsed, 2 on a
//! usage error.

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fs::{self, File};
use std::io::{BufRead, BufWriter, Write};
use std::os::unix::fs::FileExt;
use std::path::Path;

const USAGE: &str = "usage: og dev evidence clashes RAW.vcf.gz SPLIT.vcf.gz DATA_DIR OUT.tsv\n       og dev evidence hg002-gvcf TRUTH.vcf.gz TRUTH.bed REF.fasta SITES.tsv OUT.vcf\n";

/// Runs `og dev evidence ARGS`; returns the exit code.
pub fn cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let usage = |err: &mut dyn Write, msg: &str| {
        let _ = writeln!(err, "{USAGE}og dev evidence: error: {msg}");
        2
    };
    let rest = args.get(1..).unwrap_or(&[]);
    let result = match args.first().map(String::as_str) {
        Some("clashes") if rest.len() == 4 => clashes(&rest[0], &rest[1], &rest[2], &rest[3]),
        Some("hg002-gvcf") if rest.len() == 5 => {
            hg002_gvcf(&rest[0], &rest[1], &rest[2], &rest[3], &rest[4]).map(|(records, blocks)| {
                let _ = writeln!(
                    out,
                    "gVCF: {records} GIAB records, {blocks} reference blocks at panel sites"
                );
            })
        }
        Some("clashes" | "hg002-gvcf") => return usage(err, "wrong number of arguments"),
        Some("-h" | "--help") => {
            let _ = out.write_all(USAGE.as_bytes());
            return 0;
        }
        Some(other) => return usage(err, &format!("unknown evidence '{other}'")),
        None => return usage(err, "the following arguments are required: evidence"),
    };
    match result {
        Ok(()) => 0,
        Err(e) => {
            let _ = writeln!(err, "og dev evidence: {e}");
            1
        }
    }
}

fn open(path: &str) -> Result<Box<dyn BufRead>, String> {
    og_io::open(path).map_err(|e| format!("{path}: {e}"))
}

/// Each line of `path` without its line end, to `f`.
fn lines(path: &str, mut f: impl FnMut(&[u8]) -> Result<(), String>) -> Result<(), String> {
    let mut r = open(path)?;
    let mut buf = Vec::new();
    loop {
        buf.clear();
        let n = r
            .read_until(b'\n', &mut buf)
            .map_err(|e| format!("{path}: {e}"))?;
        if n == 0 {
            return Ok(());
        }
        let line = buf.strip_suffix(b"\n").unwrap_or(&buf);
        f(line.strip_suffix(b"\r").unwrap_or(line))?;
    }
}

fn pos(field: &[u8], path: &str) -> Result<u64, String> {
    std::str::from_utf8(field)
        .ok()
        .and_then(|s| s.trim().parse().ok())
        .ok_or_else(|| {
            format!(
                "{path}: position '{}' is not a number",
                String::from_utf8_lossy(field)
            )
        })
}

/// The positions of `path` on `chroms` with two records or more, one of
/// them an SNV. The VCF is sorted, so a position's records are adjacent.
fn clash_positions(
    path: &str,
    chroms: &BTreeSet<Vec<u8>>,
) -> Result<BTreeSet<(Vec<u8>, u64)>, String> {
    let mut hit = BTreeSet::new();
    let (mut key, mut n, mut snvs): (Option<(Vec<u8>, u64)>, u32, u32) = (None, 0, 0);
    lines(path, |line| {
        if line.first() == Some(&b'#') {
            return Ok(());
        }
        let c: Vec<&[u8]> = line.splitn(6, |&b| b == b'\t').collect();
        if c.len() < 5 || !chroms.contains(c[0]) {
            return Ok(());
        }
        let k = (c[0].to_vec(), pos(c[1], path)?);
        if key.as_ref() != Some(&k) {
            key = Some(k.clone());
            (n, snvs) = (0, 0);
        }
        n += 1;
        let one = |a: &[u8]| std::str::from_utf8(a).map_or(a.len(), |s| s.chars().count()) == 1;
        snvs += u32::from(one(c[3]) && one(c[4]));
        if n >= 2 && snvs >= 1 {
            hit.insert(k);
        }
        Ok(())
    })?;
    Ok(hit)
}

fn clashes(raw: &str, split: &str, data: &str, out: &str) -> Result<(), String> {
    let mut tables = BTreeMap::new();
    for e in fs::read_dir(data).map_err(|e| format!("{data}: {e}"))? {
        let name = e
            .map_err(|e| format!("{data}: {e}"))?
            .file_name()
            .into_string()
            .unwrap_or_default();
        if let Some(chrom) = name
            .strip_suffix(".tsv.gz")
            .filter(|c| c.starts_with("chr"))
        {
            tables.insert(chrom.as_bytes().to_vec(), Path::new(data).join(&name));
        }
    }
    let chroms: BTreeSet<Vec<u8>> = tables.keys().cloned().collect();
    let sets = [
        ("raw", clash_positions(raw, &chroms)?),
        ("split", clash_positions(split, &chroms)?),
    ];
    let mut on_table = [0u64; 2];
    for (chrom, table) in &tables {
        let want: Vec<BTreeSet<u64>> = sets
            .iter()
            .map(|(_, s)| {
                s.iter()
                    .filter(|(c, _)| c == chrom)
                    .map(|(_, p)| *p)
                    .collect()
            })
            .collect();
        if want.iter().all(BTreeSet::is_empty) {
            continue;
        }
        let path = table.to_string_lossy().into_owned();
        lines(&path, |line| {
            let c: Vec<&[u8]> = line.splitn(3, |&b| b == b'\t').collect();
            if c[0] != chrom.as_slice() || c.len() < 2 {
                return Ok(());
            }
            let p = pos(c[1], &path)?;
            for (k, w) in want.iter().enumerate() {
                on_table[k] += u64::from(w.contains(&p));
            }
            Ok(())
        })?;
    }
    let mut text = String::new();
    for (k, (name, s)) in sets.iter().enumerate() {
        text += &format!(
            "{name}_positions\t{}\n{name}_on_table\t{}\n",
            s.len(),
            on_table[k]
        );
    }
    fs::write(out, text).map_err(|e| format!("{out}: {e}"))
}

/// One `.fai` row: length, byte offset, bases per line, bytes per line.
struct Fai {
    len: u64,
    offset: u64,
    line_bases: u64,
    line_bytes: u64,
}

fn read_fai(path: &str) -> Result<HashMap<Vec<u8>, Fai>, String> {
    let mut fai = HashMap::new();
    lines(path, |line| {
        let f: Vec<&[u8]> = line.split(|&b| b == b'\t').collect();
        let n = |i: usize| {
            f.get(i)
                .ok_or_else(|| format!("{path}: a row has fewer than 5 columns"))
                .and_then(|x| pos(x, path))
        };
        let row = Fai {
            len: n(1)?,
            offset: n(2)?,
            line_bases: n(3)?,
            line_bytes: n(4)?,
        };
        if row.line_bases == 0 {
            return Err(format!("{path}: a row has 0 bases per line"));
        }
        fai.insert(f[0].to_vec(), row);
        Ok(())
    })?;
    Ok(fai)
}

/// Sorted, merged half-open spans `[start, end)` of 1-based positions.
fn merge(mut spans: Vec<(u64, u64)>) -> Vec<(u64, u64)> {
    spans.sort_unstable();
    let mut out: Vec<(u64, u64)> = Vec::with_capacity(spans.len());
    for (s, e) in spans {
        match out.last_mut() {
            Some(last) if s <= last.1 => last.1 = last.1.max(e),
            _ => out.push((s, e)),
        }
    }
    out
}

/// Is `p` inside one of the sorted, disjoint spans?
fn inside(spans: &[(u64, u64)], p: u64) -> bool {
    let i = spans.partition_point(|&(s, _)| s <= p);
    i > 0 && p < spans[i - 1].1
}

fn hg002_gvcf(
    truth: &str,
    bed: &str,
    reference: &str,
    sites: &str,
    out: &str,
) -> Result<(u64, u64), String> {
    let fai = read_fai(&format!("{reference}.fai"))?;
    let fasta = File::open(reference).map_err(|e| format!("{reference}: {e}"))?;
    // The benchmark regions, BED 0-based: (start, end) per chromosome,
    // sorted. A site p (1-based) is confident when the region with the
    // largest start below p ends at p or later, as the Python's bisect did.
    let mut conf: HashMap<Vec<u8>, Vec<(u64, u64)>> = HashMap::new();
    lines(bed, |line| {
        let f: Vec<&[u8]> = line
            .split(|b| b.is_ascii_whitespace())
            .filter(|x| !x.is_empty())
            .collect();
        if f.is_empty() {
            return Ok(());
        }
        if f.len() < 3 {
            return Err(format!("{bed}: a row has fewer than 3 columns"));
        }
        conf.entry(f[0].to_vec())
            .or_default()
            .push((pos(f[1], bed)?, pos(f[2], bed)?));
        Ok(())
    })?;
    for v in conf.values_mut() {
        v.sort_unstable();
    }
    let confident = |c: &[u8], p: u64| {
        conf.get(c).is_some_and(|v| {
            let i = v.partition_point(|&(s, _)| s < p);
            i > 0 && v[i - 1].1 >= p
        })
    };

    let file = File::create(out).map_err(|e| format!("{out}: {e}"))?;
    let mut w = BufWriter::with_capacity(1 << 20, file);
    let werr = |e: std::io::Error| format!("{out}: {e}");
    let mut covered: HashMap<Vec<u8>, Vec<(u64, u64)>> = HashMap::new();
    let mut records = 0u64;
    lines(truth, |line| {
        if line.starts_with(b"##") {
        } else if line.first() == Some(&b'#') {
            w.write_all(b"##INFO=<ID=END,Number=1,Type=Integer,Description=\"End of the reference block\">\n##ALT=<ID=*,Description=\"Any other allele\">\n")
                .map_err(werr)?;
        } else {
            let f: Vec<&[u8]> = line.splitn(6, |&b| b == b'\t').collect();
            if f.len() < 5 {
                return Err(format!("{truth}: a record has fewer than 5 columns"));
            }
            let p = pos(f[1], truth)?;
            covered
                .entry(f[0].to_vec())
                .or_default()
                .push((p, p + f[3].len() as u64));
            records += 1;
        }
        w.write_all(line)
            .and_then(|()| w.write_all(b"\n"))
            .map_err(werr)
    })?;
    let covered: HashMap<Vec<u8>, Vec<(u64, u64)>> =
        covered.into_iter().map(|(c, v)| (c, merge(v))).collect();

    let mut blocks = 0u64;
    let mut first = true;
    let mut base = [0u8; 1];
    lines(sites, |line| {
        if std::mem::take(&mut first) {
            return Ok(());
        }
        let f: Vec<&[u8]> = line.splitn(4, |&b| b == b'\t').collect();
        if f.len() < 3 {
            return Err(format!("{sites}: a row has fewer than 3 columns"));
        }
        let (c, p) = (f[1], pos(f[2], sites)?);
        let Some(row) = fai.get(c) else { return Ok(()) };
        if p == 0 || p > row.len || covered.get(c).is_some_and(|v| inside(v, p)) || !confident(c, p)
        {
            return Ok(());
        }
        let at =
            row.offset + ((p - 1) / row.line_bases) * row.line_bytes + (p - 1) % row.line_bases;
        fasta
            .read_exact_at(&mut base, at)
            .map_err(|e| format!("{reference}: {e}"))?;
        let b = base[0].to_ascii_uppercase();
        if !b"ACGT".contains(&b) {
            return Ok(());
        }
        blocks += 1;
        w.write_all(c)
            .and_then(|()| writeln!(w, "\t{p}\t.\t{}\t<*>\t.\t.\tEND={p}\tGT\t0/0", b as char))
            .map_err(werr)
    })?;
    w.flush().map_err(werr)?;
    Ok((records, blocks))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hg002_gvcf_keeps_truth_and_adds_only_valid_reference_blocks() {
        let repo = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../..");
        let inputs = repo.join("tests/parity/inputs/dev-evidence-hg002-gvcf");
        let temp = std::env::temp_dir().join(format!(
            "og-evidence-gvcf-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir(&temp).unwrap();
        let output = temp.join("out.vcf");
        let args = [
            "hg002-gvcf".to_owned(),
            inputs.join("truth.vcf.gz").to_str().unwrap().to_owned(),
            inputs
                .join("truth-past-end.bed")
                .to_str()
                .unwrap()
                .to_owned(),
            inputs.join("ref.fasta").to_str().unwrap().to_owned(),
            inputs
                .join("sites-past-end.tsv")
                .to_str()
                .unwrap()
                .to_owned(),
            output.to_str().unwrap().to_owned(),
        ];
        let (mut out, mut err) = (Vec::new(), Vec::new());
        assert_eq!(cli(&args, &mut out, &mut err), 0, "{err:?}");
        assert_eq!(
            out,
            b"gVCF: 4 GIAB records, 1 reference blocks at panel sites\n"
        );
        assert!(err.is_empty());
        let actual = fs::read(&output).unwrap();
        fs::remove_dir_all(&temp).unwrap();
        assert_eq!(
            actual,
            include_bytes!(
                "../../../../tests/parity/golden/dev-evidence-hg002-gvcf/past-the-end/out.vcf"
            )
        );
    }

    #[test]
    fn spans_merge_and_answer_inside() {
        let m = merge(vec![(10, 12), (5, 6), (11, 15), (15, 16)]);
        assert_eq!(m, [(5, 6), (10, 16)]);
        assert!(
            inside(&m, 5) && !inside(&m, 6) && inside(&m, 15) && !inside(&m, 16) && !inside(&m, 4)
        );
    }
}
