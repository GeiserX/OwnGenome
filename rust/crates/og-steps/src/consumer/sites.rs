//! `og consumer build-sites`: the site table of step 40.
//!
//! The 23andMe v5 site list (GRCh37 positions and reference bases, columns
//! chrom, pos, id, ref) joined with dbSNP's GRCh38 VCF by rsID. The GRCh38
//! position of a site is dbSNP's, not a liftover. A site is dropped, with
//! its reason, when it has no rsID (23andMe's own i-numbers), dbSNP has no
//! SNV for it, dbSNP places it twice or on another chromosome, or neither
//! strand of its GRCh37 base is one of dbSNP's alleles. MT sites take
//! dbSNP's rCRS position in both builds (the list numbers them on hg19's
//! chrM).
//!
//! dbSNP's file (about 30 GB, BGZF) is read once, inflated on many cores;
//! a line is parsed past its ID column only when the ID is on the list.
//!
//! Deliberately unlike `bin/export_consumer.py`: when dbSNP has two SNV
//! records of one rsID at one position, their ALT alleles are joined before
//! the strand test. The Python kept the first record's only, so a site whose
//! GRCh37 base was an ALT of the second record could be dropped or get the
//! wrong strand.

use std::collections::BTreeSet;
use std::io::{BufRead, Write};

use flate2::Compression;
use flate2::write::GzEncoder;
use og_io::{AtomicFile, Lines, Utf8};
use pycompat::num::parse_int;
use pycompat::text::isdigit;
use rustc_hash::FxHashMap;

use super::{Args, CHROMS, Fail, MT, Usage, chrom_index, head, py_repr, report};

/// One row of the template: rsID, chromosome, GRCh37 position and base.
struct Row {
    rsid: String,
    chrom: u8,
    pos37: i128,
    ref37: String,
}

/// One dbSNP SNV record of a wanted rsID: chromosome name, GRCh38 position,
/// REF and ALTs, upper case.
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord)]
struct Hit {
    chrom: &'static str,
    pos: i128,
    ref_: Vec<u8>,
    alts: Vec<Vec<u8>>,
}

/// `tsv_chrom`: chr1 / 1 / chrM / MT -> 23andMe's name, or None.
fn tsv_chrom(c: &str) -> Option<u8> {
    let c = c.strip_prefix("chr").unwrap_or(c);
    chrom_index(if c == "M" { "MT" } else { c })
}

/// dbSNP's GRCh38 RefSeq accession (without its version) -> 23andMe's name.
fn refseq38(acc: &[u8]) -> Option<&'static str> {
    if acc == b"NC_012920" {
        return Some("MT");
    }
    let n = acc.strip_prefix(b"NC_0000")?;
    if n.len() != 2 || !n.iter().all(u8::is_ascii_digit) {
        return None;
    }
    let i = usize::from(n[0] - b'0') * 10 + usize::from(n[1] - b'0');
    match i {
        1..=22 => Some(CHROMS[i - 1]),
        23 => Some("X"),
        24 => Some("Y"),
        _ => None,
    }
}

fn complement(b: &str) -> Option<&'static str> {
    match b {
        "A" => Some("T"),
        "C" => Some("G"),
        "G" => Some("C"),
        "T" => Some("A"),
        _ => None,
    }
}

/// Every line of `path` as Python's text mode reads it (universal
/// newlines, undecodable bytes as U+FFFD), with `f(number, line)`.
fn each_line(path: &str, mut f: impl FnMut(usize, &str) -> Result<(), Fail>) -> Result<(), Fail> {
    let input = og_io::open(path).map_err(|e| Fail::io(path, &e))?;
    let mut lines = Lines::new(input, Utf8::Replace);
    let mut line = String::new();
    let mut n = 0;
    while lines.read_line(&mut line).map_err(|e| Fail::io(path, &e))? {
        n += 1;
        f(n, &line)?;
    }
    Ok(())
}

fn read_template(path: &str) -> Result<Vec<Row>, Fail> {
    let mut rows = Vec::new();
    each_line(path, |n, line| {
        let line = line.trim_end_matches(['\r', '\n']);
        if line.is_empty() || line.starts_with('#') {
            return Ok(());
        }
        let c: Vec<&str> = line.split('\t').collect();
        if c.len() < 4 {
            return Err(Fail::input(format!(
                "{path} line {n}: want chrom, pos, id, ref; got {}",
                py_repr(head(line, 80))
            )));
        }
        let chrom = tsv_chrom(c[0]);
        let (Some(chrom), true) = (chrom, isdigit(c[1])) else {
            return Err(Fail::input(format!(
                "{path} line {n}: unknown chromosome or position: {}",
                py_repr(head(line, 80))
            )));
        };
        let pos37 = parse_int(c[1])
            .map_err(|_| Fail::data(format!("{path} line {n}: position {}", py_repr(c[1]))))?;
        rows.push(Row {
            rsid: c[2].to_owned(),
            chrom,
            pos37,
            ref37: c[3].to_uppercase(),
        });
        Ok(())
    })?;
    if rows.is_empty() {
        return Err(Fail::input(format!("{path} holds no site")));
    }
    Ok(rows)
}

/// `rsid.startswith("rs") and rsid[2:].isdigit()`.
fn is_rsid(s: &str) -> bool {
    s.strip_prefix("rs").is_some_and(isdigit)
}

/// What dbSNP holds for the wanted rsIDs.
#[derive(Default)]
struct Found {
    /// SNV records, per wanted rsID (the slot of [`Want`]).
    snv: Vec<Vec<Hit>>,
    /// rsIDs dbSNP has as a non-SNV record (indel, MNV).
    other: Vec<bool>,
}

/// `s.rstrip()` on bytes with Python's ASCII whitespace.
fn rstrip_ws(s: &[u8]) -> &[u8] {
    let n = s
        .iter()
        .rposition(|c| !matches!(c, b' ' | b'\t' | b'\n' | b'\r' | 0x0b | 0x0c | 0x1c..=0x1f))
        .map_or(0, |i| i + 1);
    &s[..n]
}

fn has_snv_class(info: &[u8]) -> bool {
    info.split(|&c| c == b';').any(|item| item == b"VC=SNV")
}

/// Reads dbSNP once: every SNV record of a wanted rsID on a primary
/// chromosome, and which wanted rsIDs it has only as something else.
fn scan_dbsnp(
    path: &str,
    want: &FxHashMap<Vec<u8>, usize>,
    workers: std::num::NonZero<usize>,
) -> Result<Found, Fail> {
    let mut found = Found {
        snv: vec![Vec::new(); want.len()],
        other: vec![false; want.len()],
    };
    let mut input = og_io::bgzf::open(path, workers).map_err(|e| Fail::io(path, &e))?;
    let mut line = Vec::with_capacity(1 << 12);
    loop {
        line.clear();
        if input
            .read_until(b'\n', &mut line)
            .map_err(|e| Fail::io(path, &e))?
            == 0
        {
            break;
        }
        if line.first() == Some(&b'#') {
            continue;
        }
        let mut body: &[u8] = &line;
        if let Some(b) = body.strip_suffix(b"\n") {
            body = b.strip_suffix(b"\r").unwrap_or(b);
        }
        // The ID column first: most lines stop here.
        let mut tabs = memchr::memchr_iter(b'\t', body);
        let (Some(_), Some(t2)) = (tabs.next(), tabs.next()) else {
            continue;
        };
        let id_end = tabs.next().unwrap_or(body.len());
        let Some(&slot) = want.get(&body[t2 + 1..id_end]) else {
            continue;
        };
        let cols: Vec<&[u8]> = body.splitn(9, |&c| c == b'\t').collect();
        if cols.len() < 8 {
            continue;
        }
        let acc = cols[0].split(|&c| c == b'.').next().unwrap_or_default();
        let Some(chrom) = refseq38(acc) else {
            // A patch or ALT scaffold: the primary assembly only.
            continue;
        };
        let ref_ = cols[3].to_ascii_uppercase();
        let alts: Vec<Vec<u8>> = cols[4]
            .split(|&c| c == b',')
            .map(<[u8]>::to_ascii_uppercase)
            .collect();
        if !has_snv_class(rstrip_ws(cols[7]))
            || ref_.len() != 1
            || alts.iter().any(|a| a.len() != 1)
        {
            found.other[slot] = true;
            continue;
        }
        let pos = std::str::from_utf8(cols[1])
            .ok()
            .and_then(|s| parse_int(s).ok())
            .ok_or_else(|| {
                Fail::data(format!(
                    "{path}: position {} is not a number",
                    String::from_utf8_lossy(cols[1])
                ))
            })?;
        let hit = Hit {
            chrom,
            pos,
            ref_,
            alts,
        };
        let hits = &mut found.snv[slot];
        if !hits.contains(&hit) {
            hits.push(hit);
        }
    }
    Ok(found)
}

/// A kept site: rsid, chrom, pos37, ref37, pos38, ref38, strand.
struct Kept {
    rsid: String,
    chrom: u8,
    pos37: i128,
    ref37: String,
    pos38: i128,
    ref38: String,
    strand: char,
}

/// A dropped site: rsid, chrom, the position written, reason.
struct Dropped {
    rsid: String,
    chrom: u8,
    pos37: i128,
    why: String,
}

/// The decision for one template row.
fn place(row: &Row, want: &FxHashMap<Vec<u8>, usize>, found: &Found) -> Result<Kept, Dropped> {
    let drop = |pos37: i128, why: String| Dropped {
        rsid: row.rsid.clone(),
        chrom: row.chrom,
        pos37,
        why,
    };
    let Some(&slot) = want.get(row.rsid.as_bytes()) else {
        return Err(drop(row.pos37, "no rsID (a 23andMe internal id)".into()));
    };
    let hits = &found.snv[slot];
    if hits.is_empty() {
        let why = if found.other[slot] {
            "not an SNV in dbSNP"
        } else {
            "not in dbSNP (merged or withdrawn)"
        };
        return Err(drop(row.pos37, why.into()));
    }
    // A pseudoautosomal site has a record on chrX and one on chrY: the
    // template's chromosome decides.
    let name = CHROMS[usize::from(row.chrom)];
    let same: Vec<&Hit> = hits.iter().filter(|h| h.chrom == name).collect();
    if same.is_empty() {
        let first = hits.iter().map(|h| h.chrom).min().unwrap_or_default();
        return Err(drop(row.pos37, format!("on chromosome {first} in GRCh38")));
    }
    let positions: BTreeSet<i128> = same.iter().map(|h| h.pos).collect();
    if positions.len() > 1 {
        return Err(drop(
            row.pos37,
            "more than one GRCh38 position in dbSNP".into(),
        ));
    }
    let first = same.iter().min().expect("not empty");
    let (pos38, ref38) = (first.pos, String::from_utf8_lossy(&first.ref_).into_owned());
    let (mut pos37, mut ref37) = (row.pos37, row.ref37.clone());
    if row.chrom == MT {
        // GRCh37's MT and GRCh38's chrM are the same sequence (the rCRS), but
        // the list numbers its MT sites the hg19 chrM way (the Yoruba
        // sequence, one or two bases off): the rCRS position is dbSNP's.
        pos37 = pos38;
        ref37.clone_from(&ref38);
    }
    // Every allele of the records at that position: REF, and the ALTs of
    // each record joined.
    let mut alleles: BTreeSet<String> = BTreeSet::from([ref38.clone()]);
    for h in &same {
        alleles.extend(
            h.alts
                .iter()
                .map(|a| String::from_utf8_lossy(a).into_owned()),
        );
    }
    let strand = if alleles.contains(&ref37) {
        '+'
    } else if complement(&ref37).is_some_and(|c| alleles.contains(c)) {
        '-'
    } else {
        // The MT position is dbSNP's here, as the Python wrote it.
        return Err(drop(
            pos37,
            "its GRCh37 base is not a dbSNP allele on either strand".into(),
        ));
    };
    Ok(Kept {
        rsid: row.rsid.clone(),
        chrom: row.chrom,
        pos37,
        ref37,
        pos38,
        ref38,
        strand,
    })
}

/// The options of `build-sites`.
pub struct BuildArgs<'a> {
    pub template: &'a str,
    pub dbsnp: &'a str,
    pub out: &'a str,
    pub dropped: &'a str,
    pub source: &'a str,
}

/// Builds the table and the dropped list, and returns the summary lines.
/// Exit 1 (after writing both) when no site is kept.
pub fn build_sites(a: &BuildArgs, workers: std::num::NonZero<usize>) -> Result<String, Fail> {
    let template = read_template(a.template)?;
    let mut want: FxHashMap<Vec<u8>, usize> = FxHashMap::default();
    for row in &template {
        if is_rsid(&row.rsid) {
            let next = want.len();
            want.entry(row.rsid.as_bytes().to_vec()).or_insert(next);
        }
    }
    let found = scan_dbsnp(a.dbsnp, &want, workers)?;
    let mut kept = Vec::new();
    let mut dropped = Vec::new();
    for row in &template {
        match place(row, &want, &found) {
            Ok(k) => kept.push(k),
            Err(d) => dropped.push(d),
        }
    }
    kept.sort_by(|x, y| (x.chrom, x.pos37, &x.rsid).cmp(&(y.chrom, y.pos37, &y.rsid)));

    let mut table = String::new();
    table.push_str("## consumer sites: 23andMe v5 chip, GRCh37 positions from the site list, GRCh38 positions from dbSNP\n");
    let source = if a.source.is_empty() { "." } else { a.source };
    table.push_str(&format!("## source: {source}\n"));
    table.push_str(&format!("## template_sites: {}\n", template.len()));
    table.push_str(&format!("## kept: {}\n", kept.len()));
    table.push_str(&format!("## dropped: {}\n", dropped.len()));
    table.push_str("#rsid\tchrom\tpos37\tref37\tpos38\tref38\tstrand\n");
    for k in &kept {
        table.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            k.rsid,
            CHROMS[usize::from(k.chrom)],
            k.pos37,
            k.ref37,
            k.pos38,
            k.ref38,
            k.strand
        ));
    }
    let mut drop_text = String::from("rsid\tchrom\tpos37\treason\n");
    for d in &dropped {
        drop_text.push_str(&format!(
            "{}\t{}\t{}\t{}\n",
            d.rsid,
            CHROMS[usize::from(d.chrom)],
            d.pos37,
            d.why
        ));
    }
    write_gzip(a.out, table.as_bytes())?;
    let mut f = AtomicFile::create(a.dropped).map_err(|e| Fail::io(a.dropped, &e))?;
    f.write_all(drop_text.as_bytes())
        .and_then(|()| f.commit())
        .map_err(|e| Fail::io(a.dropped, &e))?;

    // Reasons by count, most first; equal counts in the order first seen.
    let mut reasons: Vec<(&str, usize)> = Vec::new();
    for d in &dropped {
        match reasons.iter_mut().find(|(w, _)| *w == d.why) {
            Some(r) => r.1 += 1,
            None => reasons.push((&d.why, 1)),
        }
    }
    reasons.sort_by_key(|r| std::cmp::Reverse(r.1));
    let mut out = format!(
        "Site list: {} sites; kept {}, dropped {}\n",
        template.len(),
        kept.len(),
        dropped.len()
    );
    for (why, n) in reasons {
        out.push_str(&format!("  dropped, {why}: {n}\n"));
    }
    if kept.is_empty() {
        print!("{out}");
        return Err(Fail::data(
            "no site kept: is the dbSNP file the GRCh38 VCF?",
        ));
    }
    Ok(out)
}

/// Gzip with no name and mtime 0, so the same table gives the same bytes.
fn write_gzip(path: &str, data: &[u8]) -> Result<(), Fail> {
    let f = AtomicFile::create(path).map_err(|e| Fail::io(path, &e))?;
    let mut gz = GzEncoder::new(f, Compression::default());
    gz.write_all(data)
        .and_then(|()| gz.finish())
        .and_then(AtomicFile::commit)
        .map_err(|e| Fail::io(path, &e))
}

/// `og consumer build-sites ...`.
pub fn cli(args: &[String]) -> Result<u8, Usage> {
    let p = Args::parse(args, &["template", "dbsnp", "out", "dropped", "source"])?;
    if let Some(extra) = p.positional.first() {
        return Err(Usage(format!("unrecognized arguments: {extra}")));
    }
    let a = BuildArgs {
        template: p.need("template")?,
        dbsnp: p.need("dbsnp")?,
        out: p.need("out")?,
        dropped: p.need("dropped")?,
        source: p.get("source").unwrap_or(""),
    };
    Ok(match build_sites(&a, og_io::bgzf::workers()) {
        Ok(summary) => {
            print!("{summary}");
            0
        }
        Err(f) => report(&f),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn accessions() {
        assert_eq!(refseq38(b"NC_000001"), Some("1"));
        assert_eq!(refseq38(b"NC_000022"), Some("22"));
        assert_eq!(refseq38(b"NC_000023"), Some("X"));
        assert_eq!(refseq38(b"NC_000024"), Some("Y"));
        assert_eq!(refseq38(b"NC_012920"), Some("MT"));
        assert_eq!(refseq38(b"NC_000025"), None);
        assert_eq!(refseq38(b"NT_187361"), None);
    }

    #[test]
    fn snv_class_is_a_whole_item() {
        assert!(has_snv_class(b"RS=1;VC=SNV"));
        assert!(has_snv_class(b"VC=SNV"));
        assert!(!has_snv_class(b"RS=1;VC=SNVX"));
        assert!(!has_snv_class(b"RS=1;XVC=SNV"));
    }

    #[test]
    fn template_chromosomes() {
        assert_eq!(tsv_chrom("chr1"), Some(0));
        assert_eq!(tsv_chrom("chrM"), Some(MT));
        assert_eq!(tsv_chrom("MT"), Some(MT));
        assert_eq!(tsv_chrom("chrUn"), None);
    }
}
