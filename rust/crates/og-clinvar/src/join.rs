//! The combined allele of two neighbouring records.
//!
//! A caller can write one complex allele (a delins) as two records a few
//! bases apart on the same haplotype. Matched one by one, the parts can hit
//! two single-submitter ClinVar records and miss the better reviewed record
//! ClinVar holds for the whole allele, and `bcftools norm -m +any` joins only
//! records at the same position. So step 06 also matches, for each pair of
//! neighbours that can be one allele, the allele from the first part's start
//! to the second part's end: REF is the reference over that span, ALT is
//! ALT1, the reference bases between the parts, then ALT2.
//!
//! Two records can be one allele when the second starts at most `window`
//! bases after the first ends (a record that overlaps the other is never
//! joined) and their [`Hap`] signatures are equal.
//!
//! - [`Mode::Build`] reads the sample's split, left-aligned calls and writes
//!   the combined records, reading the reference bases from the FASTA.
//! - [`Mode::Prepare`] and [`finish`] are the two halves of the same
//!   thing around `bcftools +fill-from-fasta -- -c INFO/PGP_JOIN_SEQ`.
//! - [`tag`] marks the hits: each part of a matched combined allele gets
//!   `PGP_PART_OF`, and two hits in one gene that can be one allele but
//!   matched no combined record get `PGP_ONE_ALLELE`.

use std::collections::{BTreeMap, HashMap};
use std::fmt;
use std::io::{self, BufRead, Write};

use og_intervals::fasta::{FastaError, Indexed};

use crate::hap::{Hap, gt_of};
use crate::text::{info_add, info_get, nth, same, set_column, split, upper, whole};

/// The default distance: the second part starts at most this many bases
/// after the first ends.
pub const WINDOW: u64 = 10;

/// The INFO header lines step 06 adds, in the order it writes them.
const JOIN_SEQ_LINE: &[u8] = b"##INFO=<ID=PGP_JOIN_SEQ,Number=1,Type=String,Description=\"Step 06, temporary: the reference bases under REF\">";
const HEADER_LINES: [&[u8]; 3] = [
    b"##INFO=<ID=PGP_JOINED_FROM,Number=1,Type=String,Description=\"Step 06: the combined allele of these two neighbouring sample records, CHROM:POS:REF:ALT joined by |\">",
    b"##INFO=<ID=PGP_PART_OF,Number=1,Type=String,Description=\"Step 06: this hit is part of a combined allele that also matched ClinVar, CHROM:POS:REF:ALT (several joined by |)\">",
    b"##INFO=<ID=PGP_ONE_ALLELE,Number=1,Type=String,Description=\"Step 06: this hit and the one named, CHROM:POS:REF:ALT, are in one gene and can be one allele, but no combined ClinVar record matched: check ClinVar for one\">",
];
const PASS_LINE: &[u8] = b"##FILTER=<ID=PASS,Description=\"All filters passed\">";

/// Why a join could not finish.
#[derive(Debug)]
pub enum JoinError {
    Io(io::Error),
    Fasta(FastaError),
    /// The reference has no sequence of this name.
    NoSequence(String),
    /// A combined record of a VCF with this many samples: the record names
    /// one sample's genotype only.
    Samples(usize),
}

impl fmt::Display for JoinError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            JoinError::Io(e) => write!(f, "{e}"),
            JoinError::Fasta(e) => write!(f, "{e}"),
            JoinError::NoSequence(s) => write!(f, "the reference has no sequence {s}"),
            JoinError::Samples(n) => write!(
                f,
                "a combined allele takes one sample's genotype, and this VCF has {n} samples"
            ),
        }
    }
}

impl std::error::Error for JoinError {}

impl From<io::Error> for JoinError {
    fn from(e: io::Error) -> Self {
        JoinError::Io(e)
    }
}

impl From<FastaError> for JoinError {
    fn from(e: FastaError) -> Self {
        JoinError::Fasta(e)
    }
}

/// What the join writes.
pub enum Mode<'a> {
    /// Candidates for `bcftools +fill-from-fasta`: REF and ALT with `N` in
    /// the gap and `INFO/PGP_JOIN_SEQ=.` to fill.
    Prepare,
    /// The combined records, with the reference read from this FASTA.
    Build(&'a mut Indexed),
}

/// One record that can still pair with a later one.
struct Part {
    pos: i64,
    end: i64,
    ref_: Vec<u8>,
    alt: Vec<u8>,
    hap: Hap,
    gt: Vec<u8>,
}

/// A pair found and waiting for its turn in POS order.
struct Pair {
    chrom: Vec<u8>,
    first: (i64, Vec<u8>, Vec<u8>),
    second: (i64, Vec<u8>, Vec<u8>),
    /// Bases from the first part's end to the second's start.
    d: i64,
    fmt: &'static [u8],
    smp: Vec<u8>,
}

impl Pair {
    fn joined_from(&self) -> Vec<u8> {
        let (p1, r1, a1) = &self.first;
        let (p2, r2, a2) = &self.second;
        let c: &[u8] = &self.chrom;
        let (r1, a1, r2, a2): (&[u8], &[u8], &[u8], &[u8]) = (r1, a1, r2, a2);
        [
            c,
            b":",
            p1.to_string().as_bytes(),
            b":",
            r1,
            b":",
            a1,
            b"|",
            c,
            b":",
            p2.to_string().as_bytes(),
            b":",
            r2,
            b":",
            a2,
        ]
        .concat()
    }

    fn line(&self, ref_: &[u8], alt: &[u8], info: &[u8]) -> Vec<u8> {
        [
            &self.chrom[..],
            b"\t",
            self.first.0.to_string().as_bytes(),
            b"\t.\t",
            ref_,
            b"\t",
            alt,
            b"\t.\tPASS\t",
            info,
            b"\t",
            self.fmt,
            b"\t",
            &self.smp,
        ]
        .concat()
    }

    /// The candidate line `prepare` writes.
    fn prepare(&self) -> Vec<u8> {
        let n = vec![b'N'; (self.d - 1) as usize];
        let ref_ = [&self.first.1[..], &n, &self.second.1].concat();
        let alt = [&self.first.2[..], &n, &self.second.2].concat();
        let info = [&b"PGP_JOIN_SEQ=.;PGP_JOINED_FROM="[..], &self.joined_from()].concat();
        self.line(&ref_, &alt, &info)
    }

    /// The combined record over the reference bases `seq` (from the first
    /// part's start, as long as the candidate's REF), or none when a part
    /// disagrees with the reference, the reference is shorter, or ALT would
    /// equal REF.
    fn finish(&self, seq: &[u8]) -> Option<Vec<u8>> {
        let from = self.joined_from();
        let (p1, r1, a1) = &self.first;
        let (p2, r2, a2) = &self.second;
        combine(*p1, r1, a1, *p2, r2, a2, seq)
            .map(|(seq, alt)| self.line(&seq, &alt, &[&b"PGP_JOINED_FROM="[..], &from].concat()))
    }
}

/// REF and ALT of the combined allele of parts at `p1` and `p2`, over the
/// reference `seq` from `p1`; none when they cannot be one.
fn combine(
    p1: i64,
    r1: &[u8],
    a1: &[u8],
    p2: i64,
    r2: &[u8],
    a2: &[u8],
    seq: &[u8],
) -> Option<(Vec<u8>, Vec<u8>)> {
    let seq = upper(seq);
    let (r1, a1, r2, a2) = (upper(r1), upper(a1), upper(r2), upper(a2));
    let gap = p2 - (p1 + r1.len() as i64);
    if gap < 0 || seq.len() as i64 != p2 + r2.len() as i64 - p1 {
        return None;
    }
    let second = (p2 - p1) as usize;
    if seq[..r1.len()] != r1[..] || seq[second..] != r2[..] {
        return None;
    }
    let alt = [&a1[..], &seq[r1.len()..r1.len() + gap as usize], &a2].concat();
    (alt != seq).then_some((seq, alt))
}

/// Only these bases can be joined: no symbolic allele, no `*`, no `.`.
fn plain(allele: &[u8]) -> bool {
    !allele.is_empty() && allele.iter().all(|c| b"ACGTNacgtn".contains(c))
}

/// One line without its line ending; false at the end of the input.
fn read_line(input: &mut impl BufRead, buf: &mut Vec<u8>) -> io::Result<bool> {
    buf.clear();
    if input.read_until(b'\n', buf)? == 0 {
        return Ok(false);
    }
    if buf.last() == Some(&b'\n') {
        buf.pop();
    }
    Ok(true)
}

/// Writes a line and its line ending.
fn put(out: &mut impl Write, line: &[u8]) -> io::Result<()> {
    out.write_all(line)?;
    out.write_all(b"\n")
}

fn is_pgp_info(line: &[u8]) -> bool {
    [
        &b"JOIN_SEQ,"[..],
        b"JOINED_FROM,",
        b"PART_OF,",
        b"ONE_ALLELE,",
    ]
    .iter()
    .any(|t| {
        line.strip_prefix(b"##INFO=<ID=PGP_")
            .is_some_and(|r| r.starts_with(t))
    })
}

/// Writes the pending pairs that start before `limit` (all of them with
/// none), lowest POS first.
fn flush(
    pending: &mut BTreeMap<(i64, u64), Pair>,
    limit: Option<i64>,
    mode: &mut Mode<'_>,
    samples: usize,
    out: &mut impl Write,
) -> Result<(), JoinError> {
    while let Some(entry) = pending.first_entry() {
        if limit.is_some_and(|l| entry.key().0 >= l) {
            break;
        }
        let pair = entry.remove();
        match mode {
            Mode::Prepare => put(out, &pair.prepare())?,
            Mode::Build(fasta) => {
                if samples != 1 {
                    return Err(JoinError::Samples(samples));
                }
                let name = String::from_utf8_lossy(&pair.chrom).into_owned();
                let beg = (pair.first.0 - 1).max(0) as u64;
                let span =
                    (pair.second.0 + pair.second.1.len() as i64 - pair.first.0).max(0) as u64;
                let seq = fasta
                    .fetch(&name, beg, beg + span)?
                    .ok_or(JoinError::NoSequence(name))?;
                if let Some(rec) = pair.finish(&seq) {
                    put(out, &rec)?;
                }
            }
        }
    }
    Ok(())
}

/// Reads the sample's split, left-aligned calls (sorted by POS) and writes
/// the header and one record per pair of neighbours that can be one allele,
/// lowest POS first and, at one POS, in the order the pairs were found.
pub fn join(
    mut mode: Mode<'_>,
    window: u64,
    mut input: impl BufRead,
    out: &mut impl Write,
) -> Result<(), JoinError> {
    let window = window as i64;
    let build = matches!(mode, Mode::Build(_));
    let mut header: Vec<Vec<u8>> = Vec::new();
    let mut in_header = true;
    let mut samples = 0;
    let mut chrom: Option<Vec<u8>> = None;
    let mut parts: Vec<Part> = Vec::new();
    let mut pending: BTreeMap<(i64, u64), Pair> = BTreeMap::new();
    let mut found = 0u64;
    let mut line = Vec::new();
    while read_line(&mut input, &mut line)? {
        if line.starts_with(b"##") {
            if !is_pgp_info(&line) {
                if in_header {
                    header.push(line.clone())
                } else {
                    put(out, &line)?
                }
            }
            continue;
        }
        if line.starts_with(b"#CHROM") {
            if build && in_header && !header.iter().any(|h| h.starts_with(b"##FILTER=<ID=PASS,")) {
                // bcftools +fill-from-fasta, which build replaces, wrote it.
                header.insert(header.len().min(1), PASS_LINE.to_vec());
            }
            for h in header.drain(..) {
                put(out, &h)?;
            }
            if !build {
                put(out, JOIN_SEQ_LINE)?;
            }
            for h in HEADER_LINES {
                put(out, h)?;
            }
            put(out, &line)?;
            samples = split(&line, b'\t').len().saturating_sub(9);
            in_header = false;
            continue;
        }
        if in_header {
            for h in header.drain(..) {
                put(out, &h)?;
            }
            in_header = false;
        }
        let f = split(&line, b'\t');
        let c = nth(&f, 1);
        if chrom.as_deref().is_none_or(|x| !same(x, c)) {
            flush(&mut pending, None, &mut mode, samples, out)?;
            parts.clear();
            chrom = Some(c.to_vec());
        }
        let (ref_, alt) = (nth(&f, 4), nth(&f, 5));
        if !plain(ref_) || !plain(alt) || f.len() < 10 {
            continue;
        }
        let Some(hap) = Hap::of(f[8], f[9]) else {
            continue;
        };
        let Some(pos) = whole(f[1]).map(|p| p as i64) else {
            continue;
        };
        let gt = gt_of(f[8], f[9]).to_vec();
        let end = pos + ref_.len() as i64 - 1;
        let chrom_now = chrom.clone().unwrap_or_default();
        for b in &parts {
            let d = pos - b.end;
            if d < 1 || d > window || b.hap != hap {
                continue;
            }
            let (fmt, smp) = hap.format_columns(&hap.joined_gt(&b.gt, &gt));
            pending.insert(
                (b.pos, found),
                Pair {
                    chrom: chrom_now.clone(),
                    first: (b.pos, b.ref_.clone(), b.alt.clone()),
                    second: (pos, ref_.to_vec(), alt.to_vec()),
                    d,
                    fmt,
                    smp,
                },
            );
            found += 1;
        }
        // Keep only parts a later record can still pair with (the input is
        // sorted by POS), then write every pair that starts before them.
        parts.retain(|b| b.end >= pos - window);
        parts.push(Part {
            pos,
            end,
            ref_: ref_.to_vec(),
            alt: alt.to_vec(),
            hap,
            gt,
        });
        flush(&mut pending, Some(parts[0].pos), &mut mode, samples, out)?;
    }
    for h in header.drain(..) {
        put(out, &h)?;
    }
    flush(&mut pending, None, &mut mode, samples, out)?;
    Ok(())
}

/// `bcftools +fill-from-fasta`'s filled candidates in, the combined records
/// out: REF from INFO/PGP_JOIN_SEQ, ALT built, INFO only PGP_JOINED_FROM.
pub fn finish(mut input: impl BufRead, out: &mut impl Write) -> io::Result<()> {
    let mut line = Vec::new();
    while read_line(&mut input, &mut line)? {
        if line.starts_with(b"##INFO=<ID=PGP_JOIN_SEQ,") {
            continue;
        }
        if line.starts_with(b"#") {
            put(out, &line)?;
            continue;
        }
        let f = split(&line, b'\t');
        let info = nth(&f, 8);
        let seq = info_get(info, b"PGP_JOIN_SEQ");
        let from = info_get(info, b"PGP_JOINED_FROM");
        let ends = split(from, b'|');
        if ends.len() != 2 {
            continue;
        }
        let (a, b) = (split(ends[0], b':'), split(ends[1], b':'));
        let (Some(p1), Some(p2)) = (whole(nth(&a, 2)), whole(nth(&b, 2))) else {
            continue;
        };
        let Some((seq, alt)) = combine(
            p1 as i64,
            nth(&a, 3),
            nth(&a, 4),
            p2 as i64,
            nth(&b, 3),
            nth(&b, 4),
            seq,
        ) else {
            continue;
        };
        let mut rec = set_column(&line, 2, p1.to_string().as_bytes());
        rec = set_column(&rec, 4, &seq);
        rec = set_column(&rec, 5, &alt);
        rec = set_column(&rec, 8, &[&b"PGP_JOINED_FROM="[..], from].concat());
        put(out, &rec)?;
    }
    Ok(())
}

/// One hit of the step 06 hits VCF, as the tag pass sees it.
struct Hit {
    key: Vec<u8>,
    chrom: Vec<u8>,
    span: Option<(i64, i64)>,
    genes: Vec<Vec<u8>>,
    hap: Option<Hap>,
    combined: bool,
}

/// The gene symbols of GENEINFO (`SYMBOL:ID|SYMBOL:ID`), empty ones left out.
fn genes(info: &[u8]) -> Vec<Vec<u8>> {
    split(info_get(info, b"GENEINFO"), b'|')
        .into_iter()
        .map(|g| nth(&split(g, b':'), 1).to_vec())
        .filter(|g| !g.is_empty())
        .collect()
}

/// Marks the hits of step 06's hits VCF: `PGP_PART_OF` on each part of a
/// matched combined allele (naming every combined record it is part of, in
/// file order), and `PGP_ONE_ALLELE` on a hit with another hit in the same
/// gene, on the same haplotype, at most `window` bases away, when neither is
/// a combined record or a part of one (naming those hits in file order).
pub fn tag(input: impl BufRead, window: u64, out: &mut impl Write) -> io::Result<()> {
    let window = window as i64;
    let mut lines = Vec::new();
    for l in input.split(b'\n') {
        lines.push(l?);
    }
    let mut hits = Vec::new();
    let mut part_of: HashMap<Vec<u8>, Vec<u8>> = HashMap::new();
    for line in lines.iter().filter(|l| !l.starts_with(b"#")) {
        let f = split(line, b'\t');
        let (c, p, r, a) = (nth(&f, 1), nth(&f, 2), nth(&f, 4), nth(&f, 5));
        let key = [c, b":", p, b":", r, b":", a].concat();
        let from = info_get(nth(&f, 8), b"PGP_JOINED_FROM");
        for part in split(from, b'|') {
            part_of
                .entry(part.to_vec())
                .and_modify(|v| v.extend([&b"|"[..], &key].concat()))
                .or_insert_with(|| key.clone());
        }
        hits.push(Hit {
            key,
            chrom: c.to_vec(),
            span: whole(p).map(|p| (p as i64, p as i64 + r.len() as i64 - 1)),
            genes: genes(nth(&f, 8)),
            hap: if f.len() >= 10 {
                Hap::of(f[8], f[9])
            } else {
                None
            },
            combined: !from.is_empty(),
        });
    }
    let mut r = 0;
    for line in &lines {
        if line.starts_with(b"#") {
            put(out, line)?;
            continue;
        }
        let h = &hits[r];
        let f = split(line, b'\t');
        let info = nth(&f, 8);
        if let Some(of) = part_of.get(&h.key) {
            put(
                out,
                &set_column(
                    line,
                    8,
                    &info_add(info, &[&b"PGP_PART_OF="[..], of].concat()),
                ),
            )?;
        } else if let (false, Some(hap), Some((pos, end))) = (h.combined, &h.hap, h.span) {
            let mut other = Vec::new();
            for (j, o) in hits.iter().enumerate() {
                if j == r
                    || o.combined
                    || part_of.contains_key(&o.key)
                    || !same(&o.chrom, &h.chrom)
                    || o.hap.as_ref() != Some(hap)
                {
                    continue;
                }
                let Some((opos, oend)) = o.span else { continue };
                let d = if opos > end {
                    opos - end
                } else if pos > oend {
                    pos - oend
                } else {
                    continue;
                };
                if d > window || !h.genes.iter().any(|g| o.genes.contains(g)) {
                    continue;
                }
                if !other.is_empty() {
                    other.push(b'|');
                }
                other.extend_from_slice(&o.key);
            }
            if other.is_empty() {
                put(out, line)?;
            } else {
                put(
                    out,
                    &set_column(
                        line,
                        8,
                        &info_add(info, &[&b"PGP_ONE_ALLELE="[..], &other].concat()),
                    ),
                )?;
            }
        } else {
            put(out, line)?;
        }
        r += 1;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// GRCh38 chr22 from 30000200, as contig t from position 1.
    const SEQ: &str = "GCCTGTAATCCCAGCACTTTGGGAGGCTGAG";

    struct Ref(std::path::PathBuf);

    impl Ref {
        fn new(tag: &str) -> Ref {
            let dir = std::env::temp_dir().join(format!("og-clinvar-{}-{tag}", std::process::id()));
            std::fs::create_dir_all(&dir).unwrap();
            let fa = dir.join("ref.fa");
            std::fs::write(&fa, format!(">t\n{SEQ}\n")).unwrap();
            std::fs::write(
                dir.join("ref.fa.fai"),
                format!("t\t{}\t3\t{}\t{}\n", SEQ.len(), SEQ.len(), SEQ.len() + 1),
            )
            .unwrap();
            Ref(dir)
        }
    }

    impl Drop for Ref {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }

    fn vcf(samples: &str, records: &[&str]) -> String {
        format!(
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{samples}\n{}\n",
            records.join("\n")
        )
    }

    fn build(tag: &str, text: &str) -> Result<String, JoinError> {
        let r = Ref::new(tag);
        let mut fa = Indexed::open(r.0.join("ref.fa")).unwrap();
        let mut out = Vec::new();
        join(Mode::Build(&mut fa), WINDOW, text.as_bytes(), &mut out)?;
        Ok(String::from_utf8(out).unwrap())
    }

    fn records(s: &str) -> Vec<&str> {
        s.lines().filter(|l| !l.starts_with('#')).collect()
    }

    #[test]
    fn two_het_records_become_one_allele_with_the_bases_between() {
        let out = build(
            "pair",
            &vcf(
                "S",
                &[
                    "t\t4\t.\tT\tA\t50\tPASS\t.\tGT\t1/0",
                    "t\t9\t.\tT\tC\t50\tPASS\t.\tGT\t0/1",
                ],
            ),
        )
        .unwrap();
        assert_eq!(
            records(&out),
            ["t\t4\t.\tTGTAAT\tAGTAAC\t.\tPASS\tPGP_JOINED_FROM=t:4:T:A|t:9:T:C\tGT\t0/1"]
        );
        assert!(out.contains("##FILTER=<ID=PASS,"), "{out}");
        assert!(!out.contains("PGP_JOIN_SEQ"), "{out}");
    }

    #[test]
    fn a_part_that_is_not_the_reference_or_runs_past_its_end_is_dropped() {
        let out = build(
            "drop",
            &vcf(
                "S",
                &[
                    "t\t4\t.\tC\tA\t50\tPASS\t.\tGT\t0/1",
                    "t\t9\t.\tT\tC\t50\tPASS\t.\tGT\t0/1",
                ],
            ),
        )
        .unwrap();
        assert!(records(&out).is_empty(), "{out}");
        let out = build(
            "end",
            &vcf(
                "S",
                &[
                    "t\t26\t.\tG\tA\t50\tPASS\t.\tGT\t0/1",
                    "t\t30\t.\tAGTT\tA\t50\tPASS\t.\tGT\t0/1",
                ],
            ),
        )
        .unwrap();
        assert!(records(&out).is_empty(), "{out}");
    }

    #[test]
    fn what_the_reference_cannot_answer_stops_the_join() {
        let pair = [
            "u\t4\t.\tT\tA\t50\tPASS\t.\tGT\t0/1",
            "u\t9\t.\tT\tC\t50\tPASS\t.\tGT\t0/1",
        ];
        assert!(
            matches!(build("nochr", &vcf("S", &pair)), Err(JoinError::NoSequence(n)) if n == "u")
        );
        let two = [
            "t\t4\t.\tT\tA\t50\tPASS\t.\tGT\t0/1\t0/1",
            "t\t9\t.\tT\tC\t50\tPASS\t.\tGT\t0/1\t0/1",
        ];
        assert!(matches!(
            build("two", &vcf("S\tT", &two)),
            Err(JoinError::Samples(2))
        ));
        // Two samples and no pair: nothing to write, nothing refused.
        assert!(build("two-none", &vcf("S\tT", &two[..1])).is_ok());
    }

    #[test]
    fn pairs_come_out_in_pos_order() {
        let recs = [
            "t\t2\t.\tCCTGTAATCCCA\tC\t50\tPASS\t.\tGT\t0/1",
            "t\t6\t.\tT\tG\t50\tPASS\t.\tGT\t0/1",
            "t\t9\t.\tT\tC\t50\tPASS\t.\tGT\t0/1",
            "t\t15\t.\tC\tA\t50\tPASS\t.\tGT\t0/1",
        ];
        let mut out = Vec::new();
        join(Mode::Prepare, WINDOW, vcf("S", &recs).as_bytes(), &mut out).unwrap();
        let out = String::from_utf8(out).unwrap();
        let pos: Vec<&str> = records(&out)
            .iter()
            .map(|l| l.split('\t').nth(1).unwrap())
            .collect();
        assert_eq!(pos, ["2", "6", "6", "9"]);
    }

    #[test]
    fn tag_marks_parts_and_possible_single_alleles() {
        let hits = "#CHROM\n\
            c\t37\t.\tAGC\tA\t.\tPASS\tGENEINFO=G:1\tGT\t0/1\n\
            c\t38\t.\tGCCAT\tCATGGG\t.\tPASS\tPGP_JOINED_FROM=c:37:AGC:A|c:42:T:TGGG\tGT\t0/1\n\
            c\t42\t.\tT\tTGGG\t.\tPASS\t.\tGT\t0/1\n\
            c\t300\t.\tA\tG\t.\tPASS\tGENEINFO=H:1\tGT\t0/1\n\
            c\t305\t.\tC\tT\t.\tPASS\tGENEINFO=H:1\tGT\t1/0\n";
        let mut out = Vec::new();
        tag(hits.as_bytes(), WINDOW, &mut out).unwrap();
        let out = String::from_utf8(out).unwrap();
        let info: Vec<&str> = records(&out)
            .iter()
            .map(|l| l.split('\t').nth(7).unwrap())
            .collect();
        assert_eq!(
            info,
            [
                "GENEINFO=G:1;PGP_PART_OF=c:38:GCCAT:CATGGG",
                "PGP_JOINED_FROM=c:37:AGC:A|c:42:T:TGGG",
                "PGP_PART_OF=c:38:GCCAT:CATGGG",
                "GENEINFO=H:1;PGP_ONE_ALLELE=c:305:C:T",
                "GENEINFO=H:1;PGP_ONE_ALLELE=c:300:A:G",
            ]
        );
    }
}
