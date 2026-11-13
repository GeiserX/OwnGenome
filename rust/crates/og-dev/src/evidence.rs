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
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::os::unix::fs::{FileExt, PermissionsExt};
use std::path::Path;
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};

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

fn open(path: &str, gzip: bool) -> Result<Box<dyn BufRead>, String> {
    let mut file = BufReader::new(File::open(path).map_err(|e| format!("{path}: {e}"))?);
    let head = file.fill_buf().map_err(|e| format!("{path}: {e}"))?;
    if head.starts_with(&[0x1f, 0x8b]) != gzip {
        return Err(format!(
            "{path}: wrong compression (expected {})",
            if gzip { "gzip" } else { "plain" }
        ));
    }
    if gzip {
        og_io::reader(file).map_err(|e| format!("{path}: {e}"))
    } else {
        Ok(Box::new(file))
    }
}

/// Plain repository text with Python's universal-newline normalization.
pub(crate) fn read_text(path: &Path) -> Result<String, String> {
    use std::io::Read;
    let mut text = String::new();
    open(&path.to_string_lossy(), false)?
        .read_to_string(&mut text)
        .map_err(|e| format!("{}: invalid text: {e}", path.display()))?;
    Ok(text.replace("\r\n", "\n").replace('\r', "\n"))
}

/// Strict UTF-8 and universal newlines, with a 16 MiB record limit. This
/// local reader can use og-io's bounded reader once that API is available.
pub(crate) fn lines(path: &str, f: impl FnMut(&[u8]) -> Result<(), String>) -> Result<(), String> {
    input_lines(path, false, f)
}

fn input_lines(
    path: &str,
    gzip: bool,
    mut f: impl FnMut(&[u8]) -> Result<(), String>,
) -> Result<(), String> {
    const MAX_LINE: usize = 16 << 20;
    let mut r = open(path, gzip)?;
    let mut line = Vec::new();
    let mut after_cr = false;
    loop {
        let chunk = r.fill_buf().map_err(|e| format!("{path}: {e}"))?;
        if chunk.is_empty() {
            if !line.is_empty() {
                std::str::from_utf8(&line).map_err(|e| format!("{path}: invalid UTF-8: {e}"))?;
                f(&line)?;
            }
            return Ok(());
        }
        if after_cr {
            after_cr = false;
            if chunk[0] == b'\n' {
                r.consume(1);
                continue;
            }
        }
        let end = chunk.iter().position(|b| matches!(b, b'\r' | b'\n'));
        let n = end.unwrap_or(chunk.len());
        if line.len() + n > MAX_LINE {
            return Err(format!("{path}: record exceeds {MAX_LINE} bytes"));
        }
        line.extend_from_slice(&chunk[..n]);
        if end.is_some() {
            after_cr = chunk[n] == b'\r';
            r.consume(n + 1);
            std::str::from_utf8(&line).map_err(|e| format!("{path}: invalid UTF-8: {e}"))?;
            f(&line)?;
            line.clear();
        } else {
            r.consume(n);
        }
    }
}

/// Follow the final symlink, including relative and dangling link targets.
fn resolve_destination(destination: &Path) -> Result<PathBuf, String> {
    let mut path = destination.to_owned();
    for _ in 0..40 {
        match fs::read_link(&path) {
            Ok(target) => {
                path = if target.is_absolute() {
                    target
                } else {
                    path.parent().unwrap_or(Path::new(".")).join(target)
                }
            }
            Err(e)
                if matches!(
                    e.kind(),
                    std::io::ErrorKind::InvalidInput | std::io::ErrorKind::NotFound
                ) =>
            {
                return Ok(path);
            }
            Err(e) => return Err(format!("{}: {e}", path.display())),
        }
    }
    Err(format!(
        "{}: too many destination symlinks",
        destination.display()
    ))
}

fn validate_destination(path: &Path) -> Result<Option<fs::Permissions>, String> {
    match fs::metadata(path) {
        Ok(meta) => {
            if !meta.is_file() {
                return Err(format!(
                    "{}: destination is not a regular file",
                    path.display()
                ));
            }
            let permissions = meta.permissions();
            if permissions.mode() & 0o222 == 0 {
                return Err(format!("{}: read-only destination", path.display()));
            }
            // Check ACLs and effective-user access without truncating the old file.
            File::options()
                .write(true)
                .open(path)
                .map_err(|e| format!("{}: {e}", path.display()))?;
            Ok(Some(permissions))
        }
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => Ok(None),
        Err(e) => Err(format!("{}: {e}", path.display())),
    }
}

/// A unique, exclusively created sibling. Dropping it on any read, write
/// or rename error removes this temporary file. A multi-file set is staged
/// together, but its renames are separate transactions.
pub(crate) struct PendingOutput {
    path: PathBuf,
    destination: PathBuf,
}

impl PendingOutput {
    pub(crate) fn new(destination: &Path) -> Result<(Self, File), String> {
        static NEXT: AtomicU64 = AtomicU64::new(0);
        let destination = resolve_destination(destination)?;
        let permissions = validate_destination(&destination)?;
        let name = destination.file_name().ok_or("output has no file name")?;
        loop {
            let mut tmp = std::ffi::OsString::from(".");
            tmp.push(name);
            tmp.push(format!(
                ".tmp-{}-{}",
                std::process::id(),
                NEXT.fetch_add(1, Ordering::Relaxed)
            ));
            let path = destination.with_file_name(tmp);
            match File::options().write(true).create_new(true).open(&path) {
                Ok(file) => {
                    let pending = Self {
                        path,
                        destination: destination.to_owned(),
                    };
                    if let Some(permissions) = permissions {
                        file.set_permissions(permissions)
                            .map_err(|e| format!("{}: {e}", destination.display()))?;
                    }
                    return Ok((pending, file));
                }
                Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => continue,
                Err(e) => return Err(format!("{}: {e}", destination.display())),
            }
        }
    }

    pub(crate) fn commit(self) -> Result<(), String> {
        if let Some(permissions) = validate_destination(&self.destination)? {
            // A write can clear set-user-ID or set-group-ID bits. Apply the
            // final mode after every byte is staged, immediately before rename.
            fs::set_permissions(&self.path, permissions)
                .map_err(|e| format!("{}: {e}", self.destination.display()))?;
        }
        fs::rename(&self.path, &self.destination)
            .map_err(|e| format!("{}: {e}", self.destination.display()))
    }
}

impl Drop for PendingOutput {
    fn drop(&mut self) {
        let _ = fs::remove_file(&self.path);
    }
}

pub(crate) fn publish(path: &Path, bytes: &[u8]) -> Result<(), String> {
    let (pending, mut file) = PendingOutput::new(path)?;
    file.write_all(bytes)
        .and_then(|()| file.flush())
        .map_err(|e| format!("{}: {e}", path.display()))?;
    drop(file);
    pending.commit()
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
    input_lines(path, true, |line| {
        if line.first() == Some(&b'#') {
            return Ok(());
        }
        let c: Vec<&[u8]> = line.splitn(6, |&b| b == b'\t').collect();
        if !chroms.contains(c[0]) {
            return Ok(());
        }
        if c.len() < 6 {
            return Err(format!(
                "{path}: a selected record has fewer than 5 columns or missing trailing VCF columns"
            ));
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
        input_lines(&path, true, |line| {
            let c: Vec<&[u8]> = line.splitn(3, |&b| b == b'\t').collect();
            if c[0] != chrom.as_slice() {
                return Ok(());
            }
            if c.len() < 2 {
                return Err(format!("{path}: a selected row has fewer than 2 columns"));
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
    publish(Path::new(out), text.as_bytes())
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
    let mut magic = [0; 2];
    if fasta
        .read_at(&mut magic, 0)
        .map_err(|e| format!("{reference}: {e}"))?
        == 2
        && magic == [0x1f, 0x8b]
    {
        return Err(format!("{reference}: wrong compression (expected plain)"));
    }
    // The benchmark regions, BED 0-based: (start, end) per chromosome,
    // sorted. A site p (1-based) is confident when the region with the
    // largest start below p ends at p or later, as the Python's bisect did.
    let mut conf: HashMap<Vec<u8>, Vec<(u64, u64)>> = HashMap::new();
    lines(bed, |line| {
        let f: Vec<&[u8]> = line
            .split(|b| matches!(b, b' ' | b'\t'))
            .filter(|x| !x.is_empty())
            .collect();
        if f.iter().any(|field| {
            std::str::from_utf8(field)
                .expect("strict reader")
                .chars()
                .any(|ch| ch.is_whitespace() || matches!(ch, '\u{1c}'..='\u{1f}'))
        }) {
            return Err(format!("{bed}: unsupported whitespace in a field"));
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

    let (pending, file) = PendingOutput::new(Path::new(out))?;
    let mut w = BufWriter::with_capacity(1 << 20, file);
    let werr = |e: std::io::Error| format!("{out}: {e}");
    let mut covered: HashMap<Vec<u8>, Vec<(u64, u64)>> = HashMap::new();
    let mut records = 0u64;
    input_lines(truth, true, |line| {
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
            covered.entry(f[0].to_vec()).or_default().push((
                p,
                p.checked_add(
                    std::str::from_utf8(f[3])
                        .expect("strict reader")
                        .chars()
                        .count() as u64,
                )
                .ok_or_else(|| format!("{truth}: position span overflows"))?,
            ));
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
        std::str::from_utf8(&base).map_err(|e| format!("{reference}: invalid UTF-8: {e}"))?;
        let b = base[0].to_ascii_uppercase();
        if !b"ACGT".contains(&b) {
            return Ok(());
        }
        blocks += 1;
        w.write_all(c)
            .and_then(|()| writeln!(w, "\t{p}\t.\t{}\t<*>\t.\t.\tEND={p}\tGT\t0/0", b as char))
            .map_err(werr)
    })?;
    if first {
        return Err(format!("{sites}: missing header"));
    }
    let file = w.into_inner().map_err(|e| werr(e.into_error()))?;
    drop(file);
    pending.commit()?;
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

    struct Case {
        dir: std::path::PathBuf,
        kind: &'static str,
    }

    impl Case {
        fn new(kind: &'static str, label: &str) -> Self {
            let dir = std::env::temp_dir().join(format!(
                "og-evidence-{kind}-{label}-{}-{:?}",
                std::process::id(),
                std::thread::current().id()
            ));
            let _ = fs::remove_dir_all(&dir);
            fs::create_dir_all(&dir).unwrap();
            let src = Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("../../../tests/parity/inputs")
                .join(format!("dev-evidence-{kind}"));
            for e in fs::read_dir(src).unwrap() {
                let e = e.unwrap();
                fs::copy(e.path(), dir.join(e.file_name())).unwrap();
            }
            Self { dir, kind }
        }
        fn text(&self, name: &str) -> Vec<u8> {
            use std::io::Read;
            let mut bytes = Vec::new();
            open(self.dir.join(name).to_str().unwrap(), name.ends_with(".gz"))
                .unwrap()
                .read_to_end(&mut bytes)
                .unwrap();
            bytes
        }
        fn put(&self, name: &str, bytes: &[u8]) {
            if name.ends_with(".gz") {
                let mut gz =
                    flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
                gz.write_all(bytes).unwrap();
                fs::write(self.dir.join(name), gz.finish().unwrap()).unwrap();
            } else {
                fs::write(self.dir.join(name), bytes).unwrap();
            }
        }
        fn run(&self) -> Result<(), String> {
            let p = |n: &str| self.dir.join(n).to_str().unwrap().to_owned();
            let args = if self.kind == "clashes" {
                vec![
                    "clashes".to_owned(),
                    p("raw.vcf.gz"),
                    p("split.vcf.gz"),
                    self.dir.to_str().unwrap().to_owned(),
                    p("out"),
                ]
            } else {
                vec![
                    "hg002-gvcf".to_owned(),
                    p("truth.vcf.gz"),
                    p("truth.bed"),
                    p("ref.fasta"),
                    p("sites.tsv"),
                    p("out"),
                ]
            };
            let (mut out, mut err) = (Vec::new(), Vec::new());
            match cli(&args, &mut out, &mut err) {
                0 => Ok(()),
                1 => Err(String::from_utf8(err).unwrap()),
                code => panic!("unexpected CLI exit {code}: {err:?}"),
            }
        }
        fn fails_without_publication(&self) {
            fs::write(self.dir.join("out"), b"PREVIOUS OUTPUT\n").unwrap();
            assert!(self.run().is_err(), "{} must fail", self.dir.display());
            assert_eq!(
                fs::read(self.dir.join("out")).unwrap(),
                b"PREVIOUS OUTPUT\n"
            );
            let leftovers: Vec<_> = fs::read_dir(&self.dir)
                .unwrap()
                .filter_map(|e| {
                    let name = e.unwrap().file_name();
                    name.to_str()
                        .filter(|n| n.starts_with(".out."))
                        .map(str::to_owned)
                })
                .collect();
            assert!(leftovers.is_empty(), "{leftovers:?}");
        }
    }
    impl Drop for Case {
        fn drop(&mut self) {
            fs::remove_dir_all(&self.dir).unwrap();
        }
    }

    #[test]
    fn fai_offsets_reject_overflow_and_out_of_file_before_publication() {
        for (offset, bases, bytes, position) in [
            (18446744073709551598u64, 1, 2, 13),
            (1, 1, u64::MAX, 3),
            (u64::MAX, 1, 1, 1),
            (1000, 50, 51, 1),
        ] {
            let c = Case::new("hg002-gvcf", &format!("offset-{offset}-{position}"));
            c.put("ref.fasta.fai", format!("chrA\t50\t{offset}\t{bases}\t{bytes}\n").as_bytes());
            c.put("sites.tsv", format!("id\tchr\tpos\nx\tchrA\t{position}\n").as_bytes());
            c.put("truth.bed", b"chrA\t0\t50\n");
            c.put("truth.vcf.gz", b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\n");
            c.fails_without_publication();
            assert!(c.run().unwrap_err().contains(": invalid FASTA offset:"));
        }
    }

    #[test]
    fn late_vcf_headers_fail_before_publication() {
        for header in [b"##late=yes\n".as_slice(), b"#CHROM\tPOS\n"] {
            let c = Case::new("hg002-gvcf", &format!("late-{}", header.len()));
            c.put("truth.vcf.gz", &[c.text("truth.vcf.gz"), header.to_vec()].concat());
            c.fails_without_publication();
            assert!(c.run().unwrap_err().contains(": VCF header after first data record"));
        }
    }

    #[test]
    fn empty_gzip_inputs_and_unused_zero_width_fai_rows_fail() {
        for (kind, name) in [
            ("clashes", "raw.vcf.gz"), ("clashes", "split.vcf.gz"),
            ("clashes", "chrA.tsv.gz"), ("hg002-gvcf", "truth.vcf.gz"),
        ] {
            let c = Case::new(kind, name);
            fs::write(c.dir.join(name), b"").unwrap();
            c.fails_without_publication();
            assert!(c.run().unwrap_err().contains(": wrong compression (expected gzip)"));
        }
        let c = Case::new("hg002-gvcf", "unused-row");
        c.put("ref.fasta.fai", &[c.text("ref.fasta.fai"), b"chrUnused\t0\t0\t0\t0\n".to_vec()].concat());
        c.fails_without_publication();
        assert!(c.run().unwrap_err().contains(": a row has 0 bases per line"));
    }

    #[test]
    fn gzip_fasta_with_ascii_fname_is_refused_by_format_guard() {
        let c = Case::new("hg002-gvcf", "fname");
        let mut gz = flate2::GzBuilder::new().filename(vec![b'A'; 200])
            .write(Vec::new(), flate2::Compression::default());
        gz.write_all(&c.text("ref.fasta")).unwrap();
        fs::write(c.dir.join("ref.fasta"), gz.finish().unwrap()).unwrap();
        c.fails_without_publication();
        assert!(c.run().unwrap_err().contains(": wrong compression (expected plain)"));
    }

    #[test]
    fn non_writable_parents_preserve_evidence_destinations() {
        for kind in ["clashes", "hg002-gvcf"] {
            let c = Case::new(kind, "parent-mode");
            fs::write(c.dir.join("out"), b"PREVIOUS OUTPUT\n").unwrap();
            fs::set_permissions(&c.dir, fs::Permissions::from_mode(0o555)).unwrap();
            let result = c.run();
            fs::set_permissions(&c.dir, fs::Permissions::from_mode(0o755)).unwrap();
            assert!(result.unwrap_err().contains(": parent directory is not writable"));
            assert_eq!(fs::read(c.dir.join("out")).unwrap(), b"PREVIOUS OUTPUT\n");
        }
    }

    #[test]
    fn selected_truncated_clash_records_fail() {
        for name in ["raw.vcf.gz", "split.vcf.gz", "chrA.tsv.gz"] {
            let c = Case::new("clashes", name);
            c.put(
                name,
                if name == "chrA.tsv.gz" {
                    b"chrA"
                } else {
                    b"chrA\t3\n"
                },
            );
            c.fails_without_publication();
        }
    }

    #[test]
    fn missing_panel_header_and_blank_bed_fail() {
        for name in ["sites.tsv", "truth.bed"] {
            let c = Case::new("hg002-gvcf", name);
            let bytes = if name == "sites.tsv" {
                vec![]
            } else {
                [c.text(name), b"\n".to_vec()].concat()
            };
            c.put(name, &bytes);
            c.fails_without_publication();
        }
    }

    #[test]
    fn every_evidence_text_input_rejects_invalid_utf8() {
        for (kind, names) in [
            ("clashes", vec!["raw.vcf.gz", "split.vcf.gz", "chrA.tsv.gz"]),
            (
                "hg002-gvcf",
                vec!["truth.vcf.gz", "truth.bed", "ref.fasta.fai", "sites.tsv"],
            ),
        ] {
            for name in names {
                let c = Case::new(kind, name);
                // Invalid bytes in a header or filtered-out record must still fail.
                let mut bytes = c.text(name);
                bytes.extend_from_slice(b"#ignored\xff\n");
                c.put(name, &bytes);
                c.fails_without_publication();
            }
        }
    }

    #[test]
    fn every_evidence_text_input_has_universal_newlines() {
        for (kind, names) in [
            (
                "clashes",
                vec!["raw.vcf.gz", "split.vcf.gz", "chrA.tsv.gz", "chrB.tsv.gz"],
            ),
            (
                "hg002-gvcf",
                vec!["truth.vcf.gz", "truth.bed", "ref.fasta.fai", "sites.tsv"],
            ),
        ] {
            for name in names {
                for ending in [b"\r".as_slice(), b"\r\n".as_slice()] {
                    let c = Case::new(kind, &format!("{name}-{}", ending.len()));
                    c.run().unwrap();
                    let expected = fs::read(c.dir.join("out")).unwrap();
                    let bytes = c
                        .text(name)
                        .split(|b| *b == b'\n')
                        .collect::<Vec<_>>()
                        .join(ending);
                    c.put(name, &bytes);
                    c.run().unwrap();
                    assert_eq!(
                        fs::read(c.dir.join("out")).unwrap(),
                        expected,
                        "{kind}/{name}/{ending:?}"
                    );
                }
            }
        }
    }

    #[test]
    fn failed_gvcf_inputs_preserve_existing_output() {
        for name in ["truth.vcf.gz", "sites.tsv"] {
            for fault in ["missing", "directory", "malformed", "truncated"] {
                let c = Case::new("hg002-gvcf", &format!("{name}-{fault}"));
                let path = c.dir.join(name);
                match fault {
                    "missing" => fs::remove_file(path).unwrap(),
                    "directory" => {
                        fs::remove_file(&path).unwrap();
                        fs::create_dir(path).unwrap();
                    }
                    "malformed" => {
                        let mut bytes = c.text(name);
                        bytes.extend_from_slice(b"bad\n");
                        c.put(name, &bytes);
                    }
                    _ if name.ends_with(".gz") => {
                        let mut bytes = fs::read(&path).unwrap();
                        bytes.truncate(bytes.len() - 5);
                        fs::write(path, bytes).unwrap();
                    }
                    _ => c.put(name, b"header\nbad\n"),
                }
                c.fails_without_publication();
            }
        }
    }

    #[test]
    fn explicit_compression_is_required_at_every_input() {
        for (kind, names) in [
            ("clashes", vec!["raw.vcf.gz", "split.vcf.gz", "chrA.tsv.gz"]),
            (
                "hg002-gvcf",
                vec![
                    "truth.vcf.gz",
                    "truth.bed",
                    "ref.fasta.fai",
                    "ref.fasta",
                    "sites.tsv",
                ],
            ),
        ] {
            for name in names {
                let c = Case::new(kind, &format!("compression-{name}"));
                let bytes = c.text(name);
                let wrong = if name.ends_with(".gz") {
                    bytes
                } else {
                    let mut gz =
                        flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
                    gz.write_all(&bytes).unwrap();
                    gz.finish().unwrap()
                };
                fs::write(c.dir.join(name), wrong).unwrap();
                c.fails_without_publication();
                assert!(c.run().unwrap_err().contains(": wrong compression (expected "), "{name}");
            }
        }
    }

    #[test]
    fn bed_rejects_python_whitespace_inside_fields_and_as_separators() {
        for ch in [
            '\u{2028}', '\u{85}', '\u{b}', '\u{a0}', '\u{1c}', '\u{1d}', '\u{1e}', '\u{1f}',
        ] {
            for embedded in [true, false] {
                let c = Case::new("hg002-gvcf", &format!("space-{}-{embedded}", ch as u32));
                let text = String::from_utf8(c.text("truth.bed")).unwrap();
                let text = if embedded {
                    text.replacen("chrA", &format!("chr{ch}A"), 1)
                } else {
                    text.replace('\t', &ch.to_string())
                };
                c.put("truth.bed", text.as_bytes());
                c.fails_without_publication();
                assert!(c.run().unwrap_err().contains("unsupported whitespace"));
            }
        }
    }

    #[test]
    fn publication_preserves_permissions_and_follows_symlinks() {
        use std::os::unix::fs::{PermissionsExt, symlink};
        let c = Case::new("clashes", "publication");
        let target = c.dir.join("target");
        let link = c.dir.join("link");
        fs::write(&target, b"old").unwrap();
        fs::set_permissions(&target, fs::Permissions::from_mode(0o600)).unwrap();
        symlink("target", &link).unwrap();
        publish(&link, b"new").unwrap();
        assert!(
            fs::symlink_metadata(&link)
                .unwrap()
                .file_type()
                .is_symlink()
        );
        assert_eq!(fs::read(&target).unwrap(), b"new");
        assert_eq!(
            fs::metadata(&target).unwrap().permissions().mode() & 0o777,
            0o600
        );
        fs::set_permissions(&target, fs::Permissions::from_mode(0o4600)).unwrap();
        publish(&link, b"new").unwrap();
        assert_eq!(
            fs::metadata(&target).unwrap().permissions().mode() & 0o7777,
            0o4600
        );
        fs::set_permissions(&target, fs::Permissions::from_mode(0o444)).unwrap();
        assert!(publish(&link, b"bad").unwrap_err().contains("read-only"));
        assert_eq!(fs::read(&target).unwrap(), b"new");
        let dir_link = c.dir.join("dir-link");
        symlink(&c.dir, &dir_link).unwrap();
        for path in [&c.dir, &dir_link] {
            assert!(
                publish(path, b"bad")
                    .unwrap_err()
                    .contains("not a regular file")
            );
        }
    }

    #[test]
    fn approved_strict_evidence_shapes_fail_with_stable_errors() {
        for value in ["3_0", "٣", "-3", "18446744073709551616"] {
            let c = Case::new("hg002-gvcf", value);
            c.put(
                "truth.vcf.gz",
                format!("chrA\t{value}\t.\tA\tC\n").as_bytes(),
            );
            c.fails_without_publication();
            assert!(c.run().unwrap_err().contains("not a number"));
        }
        let c = Case::new("hg002-gvcf", "short-vcf");
        c.put("truth.vcf.gz", b"chrA\t3\t.\tA\n");
        c.fails_without_publication();
        assert!(c.run().unwrap_err().contains("fewer than 5 columns"));
        let c = Case::new("hg002-gvcf", "blank-fai");
        c.put("ref.fasta.fai", b"\n");
        c.fails_without_publication();
        assert!(c.run().unwrap_err().contains("fewer than 5 columns"));
        let c = Case::new("clashes", "newline-alt");
        c.put("raw.vcf.gz", b"chrA\t3\t.\tA\tC\nchrA\t3\t.\tAA\tA\n");
        c.fails_without_publication();
        assert!(c.run().unwrap_err().contains("VCF columns"));
    }

    #[test]
    fn approved_record_limit_and_final_lf_are_pinned() {
        let c = Case::new("clashes", "record-limit");
        c.put("raw.vcf.gz", &vec![b'#'; (16 << 20) + 1]);
        c.fails_without_publication();
        assert!(
            c.run()
                .unwrap_err()
                .contains("record exceeds 16777216 bytes")
        );
        let c = Case::new("hg002-gvcf", "final-lf");
        c.run().unwrap();
        let expected = fs::read(c.dir.join("out")).unwrap();
        let mut truth = c.text("truth.vcf.gz");
        assert_eq!(truth.pop(), Some(b'\n'));
        c.put("truth.vcf.gz", &truth);
        c.run().unwrap();
        let result = fs::read(c.dir.join("out")).unwrap();
        assert_eq!(result, expected);
        assert_eq!(result.last(), Some(&b'\n'));
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
