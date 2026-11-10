//! `og vcf csq-rank`: one transcript, one frequency tier and one call verdict
//! per VEP record, written as `PGP_*` INFO tags.
//!
//! A VCF goes in, the same VCF comes out with seven tags on each record:
//!
//! - **Transcript.** The CSQ entry of the MANE Select transcript
//!   (`MANE_SELECT` set), else of the canonical one (`CANONICAL=YES`), else
//!   of any transcript; within that set the most severe `IMPACT`, the first
//!   entry winning a tie. `PGP_IMPACT`, `PGP_CSQ` and `PGP_GENE` are that
//!   entry's `IMPACT`, `Consequence` and `SYMBOL`; `PGP_TX` names the set
//!   (`MANE`, `CANONICAL`, `ANY`, or `NONE` without a CSQ).
//! - **Frequency.** `MAX_AF` when the CSQ has that field, else `gnomADe_AF`
//!   and `gnomADg_AF`; the highest value over every entry (`&`-joined values
//!   too). `PGP_FREQ` is `rare` (below `--rare-af`), `common`, `unknown` (no
//!   value) or `none` (no frequency field at all).
//! - **Call.** From the first sample: GQ below `--min-gq`, DP below
//!   `--min-dp`, an alt read fraction outside `--het-lo`..`--het-hi` on a
//!   heterozygote or below `--hom-lo` on a homozygous alt call make the call
//!   `weak` (`PGP_CALL_WHY` says why, e.g. `GQ11,AB0.29`); `ok` when one of
//!   them could be judged, else `unknown`.
//!
//! The floors go into the header as `##PGP_CALL_FLOOR=`, which the summary
//! reads. Tags and header lines of an earlier pass are replaced.
//!
//! This is the rule `bin/csq_rank.awk` applied, and its output is the same
//! byte for byte (the parity cases in `tests/parity` hold it to mawk), with
//! three differences on purpose: the seven `##INFO` descriptions name
//! `og vcf csq-rank`, not the deleted awk file; a floor or `--keep` that is
//! not a number or a known word is a usage error (exit 2), where awk read it
//! as 0 or as `all`; and the input may be gzip.
//!
//! The work per record is a scan for tabs, `;`, `,` and `|`: no field is
//! copied except into the output, and only the CSQ fields the rule reads are
//! split out of each entry.

use std::collections::HashMap;
use std::fmt;
use std::io::{self, BufRead, Write};

use memchr::{memchr, memchr_iter, memmem};
use smallvec::SmallVec;

/// Which records are written.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Keep {
    /// Every record.
    All,
    /// Every record whose frequency tier is not `common`.
    NonCommon,
}

/// The thresholds, each with the text it was given as, which the
/// `##PGP_CALL_FLOOR` line repeats.
#[derive(Debug, Clone, PartialEq)]
pub struct Floors {
    pub rare_af: f64,
    pub min_gq: Floor,
    pub min_dp: Floor,
    pub het_lo: Floor,
    pub het_hi: Floor,
    pub hom_lo: Floor,
}

/// One threshold: its value and the text it was given as.
#[derive(Debug, Clone, PartialEq)]
pub struct Floor {
    pub value: f64,
    pub text: String,
}

impl Floor {
    /// Reads a decimal number (`30`, `0.2`, `.5`, `1e-3`, a sign allowed).
    pub fn parse(text: &str) -> Option<Floor> {
        let b = text.as_bytes();
        let body = b.strip_prefix(b"+").or(b.strip_prefix(b"-")).unwrap_or(b);
        let (mantissa, exp) = match body.iter().position(|&c| c == b'e' || c == b'E') {
            Some(i) => (&body[..i], Some(&body[i + 1..])),
            None => (body, None),
        };
        let digits = mantissa.iter().filter(|c| c.is_ascii_digit()).count();
        let dots = mantissa.iter().filter(|&&c| c == b'.').count();
        let mantissa_ok = digits > 0 && dots <= 1 && digits + dots == mantissa.len();
        let exp_ok = exp.is_none_or(|e| {
            let e = e.strip_prefix(b"+").or(e.strip_prefix(b"-")).unwrap_or(e);
            !e.is_empty() && e.iter().all(u8::is_ascii_digit)
        });
        if !(mantissa_ok && exp_ok) {
            return None;
        }
        let value: f64 = text.parse().ok()?;
        value.is_finite().then(|| Floor {
            value,
            text: text.to_owned(),
        })
    }

    fn default(value: f64, text: &str) -> Floor {
        Floor {
            value,
            text: text.to_owned(),
        }
    }
}

impl Default for Floors {
    /// GQ 30, DP 10, a heterozygote's alt fraction 0.2-0.8, a homozygote's
    /// at least 0.8, rare below 1%.
    fn default() -> Floors {
        Floors {
            rare_af: 0.01,
            min_gq: Floor::default(30.0, "30"),
            min_dp: Floor::default(10.0, "10"),
            het_lo: Floor::default(0.2, "0.2"),
            het_hi: Floor::default(0.8, "0.8"),
            hom_lo: Floor::default(0.8, "0.8"),
        }
    }
}

/// The transcript set the ranked entry came from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum TxSet {
    Mane,
    Canonical,
    Any,
    None,
}

impl TxSet {
    fn as_str(self) -> &'static str {
        match self {
            TxSet::Mane => "MANE",
            TxSet::Canonical => "CANONICAL",
            TxSet::Any => "ANY",
            TxSet::None => "NONE",
        }
    }
}

/// The frequency tier.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum FreqTier {
    Rare,
    Common,
    /// The CSQ has a frequency field, but no entry gives a value.
    Unknown,
    /// The CSQ has no frequency field at all.
    NoField,
}

impl FreqTier {
    fn as_str(self) -> &'static str {
        match self {
            FreqTier::Rare => "rare",
            FreqTier::Common => "common",
            FreqTier::Unknown => "unknown",
            FreqTier::NoField => "none",
        }
    }
}

/// The severity rank of an `IMPACT`: anything else ranks below MODIFIER.
fn severity(impact: &[u8]) -> u8 {
    match impact {
        b"HIGH" => 5,
        b"MODERATE" => 4,
        b"LOW" => 3,
        b"MODIFIER" => 2,
        _ => 1,
    }
}

/// The tags an earlier pass wrote, which this one replaces.
const TAGS: [&[u8]; 7] = [
    b"IMPACT",
    b"CSQ",
    b"GENE",
    b"TX",
    b"FREQ",
    b"CALL",
    b"CALL_WHY",
];

/// `PGP_<tag>=` for one of [`TAGS`].
fn is_old_tag(item: &[u8]) -> bool {
    item.strip_prefix(b"PGP_").is_some_and(|rest| {
        TAGS.iter()
            .any(|t| rest.starts_with(t) && rest.get(t.len()) == Some(&b'='))
    })
}

/// `##INFO=<ID=PGP_<tag>,` for one of [`TAGS`].
fn is_old_header(line: &[u8]) -> bool {
    line.starts_with(b"##PGP_CALL_FLOOR=")
        || line.strip_prefix(b"##INFO=<ID=PGP_").is_some_and(|rest| {
            TAGS.iter()
                .any(|t| rest.starts_with(t) && rest.get(t.len()) == Some(&b','))
        })
}

/// A whole number: one or more ASCII digits.
fn is_digits(s: &[u8]) -> bool {
    !s.is_empty() && s.iter().all(u8::is_ascii_digit)
}

/// A frequency value: `[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?`, so no sign, no
/// `nan` or `inf`, and at least one digit after a dot.
fn is_frequency(s: &[u8]) -> bool {
    let (mantissa, exp) = match s.iter().position(|&c| c == b'e' || c == b'E') {
        Some(i) => (&s[..i], Some(&s[i + 1..])),
        None => (s, None),
    };
    let mantissa_ok = match memchr(b'.', mantissa) {
        Some(d) => mantissa[..d].iter().all(u8::is_ascii_digit) && is_digits(&mantissa[d + 1..]),
        None => is_digits(mantissa),
    };
    mantissa_ok
        && exp.is_none_or(|e| is_digits(e.strip_prefix(b"+").or(e.strip_prefix(b"-")).unwrap_or(e)))
}

/// A number made of ASCII digits and at most one `.` and exponent, already
/// checked; as `strtod` reads it.
fn number(s: &[u8]) -> f64 {
    std::str::from_utf8(s)
        .ok()
        .and_then(|t| t.parse().ok())
        .unwrap_or(0.0)
}

/// Parts of `s` between `sep`, as awk's `split` gives them: none for an
/// empty string.
fn split(s: &[u8], sep: u8) -> impl Iterator<Item = &[u8]> {
    let parts = (!s.is_empty()).then(|| s.split(move |&c| c == sep));
    parts.into_iter().flatten()
}

/// Where the fields the rule reads sit in a CSQ entry.
#[derive(Debug, Default)]
struct Layout {
    /// Every name of every CSQ header line seen, with its 0-based position
    /// (a later line or a name given twice wins, as awk's `col[]`).
    col: HashMap<Vec<u8>, usize>,
    impact: Option<usize>,
    mane: Option<usize>,
    canonical: Option<usize>,
    consequence: Option<usize>,
    symbol: Option<usize>,
    /// The frequency fields, fixed at `#CHROM`.
    freq: SmallVec<[usize; 2]>,
    /// One past the highest position above.
    need: usize,
}

impl Layout {
    /// Takes the names of a `##INFO=<ID=CSQ,` line: after the last
    /// `Format: `, up to the first `"`, split on `|`.
    fn add_csq_line(&mut self, line: &[u8]) {
        let mut s = line;
        if let Some(at) = memmem::rfind(s, b"Format: ") {
            s = &s[at + b"Format: ".len()..];
        }
        if let Some(q) = memchr(b'"', s) {
            s = &s[..q];
        }
        for (i, name) in split(s, b'|').enumerate() {
            self.col.insert(name.to_vec(), i);
        }
        self.index();
    }

    /// Fixes the frequency fields: `MAX_AF`, else the two gnomAD ones.
    fn set_freq(&mut self) {
        self.freq.clear();
        if let Some(&i) = self.col.get(b"MAX_AF".as_slice()) {
            self.freq.push(i);
        } else {
            self.freq
                .extend(self.col.get(b"gnomADe_AF".as_slice()).copied());
            self.freq
                .extend(self.col.get(b"gnomADg_AF".as_slice()).copied());
        }
        self.index();
    }

    fn index(&mut self) {
        let get = |n: &[u8]| self.col.get(n).copied();
        self.impact = get(b"IMPACT");
        self.mane = get(b"MANE_SELECT");
        self.canonical = get(b"CANONICAL");
        self.consequence = get(b"Consequence");
        self.symbol = get(b"SYMBOL");
        self.need = [
            self.impact,
            self.mane,
            self.canonical,
            self.consequence,
            self.symbol,
        ]
        .into_iter()
        .flatten()
        .chain(self.freq.iter().copied())
        .map(|i| i + 1)
        .max()
        .unwrap_or(0);
    }
}

/// The fields of one CSQ entry up to `need`; a missing one reads as empty.
struct Entry<'a>(SmallVec<[&'a [u8]; 64]>);

impl<'a> Entry<'a> {
    fn new(entry: &'a [u8], need: usize) -> Entry<'a> {
        Entry(entry.splitn(need + 1, |&c| c == b'|').take(need).collect())
    }

    fn get(&self, i: Option<usize>) -> &'a [u8] {
        i.and_then(|i| self.0.get(i).copied()).unwrap_or(b"")
    }
}

/// Ranks a VCF stream.
pub struct Ranker {
    keep: Keep,
    floors: Floors,
    layout: Layout,
    /// The new INFO of the record being written, reused.
    info: Vec<u8>,
    /// What made the call weak, reused.
    why: Vec<u8>,
}

/// Why ranking stopped: reading or writing failed.
#[derive(Debug)]
pub struct Error(pub io::Error);

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}

impl std::error::Error for Error {}

impl Ranker {
    pub fn new(keep: Keep, floors: Floors) -> Ranker {
        Ranker {
            keep,
            floors,
            layout: Layout::default(),
            info: Vec::with_capacity(4096),
            why: Vec::new(),
        }
    }

    /// Reads `input` to its end and writes the ranked VCF to `out`. Each
    /// line is written with a `\n`, also a last line that had none.
    pub fn run(&mut self, mut input: impl BufRead, mut out: impl Write) -> Result<(), Error> {
        let mut line = Vec::with_capacity(1 << 16);
        loop {
            line.clear();
            if input.read_until(b'\n', &mut line).map_err(Error)? == 0 {
                break;
            }
            if line.last() == Some(&b'\n') {
                line.pop();
            }
            self.line(&line, &mut out).map_err(Error)?;
        }
        out.flush().map_err(Error)
    }

    /// One line: a header line, passed through, dropped or followed by the
    /// new header lines, or a record.
    pub fn line(&mut self, line: &[u8], out: &mut impl Write) -> io::Result<()> {
        if line.starts_with(b"##INFO=<ID=CSQ,") {
            self.layout.add_csq_line(line);
        } else if is_old_header(line) {
            return Ok(());
        } else if line.starts_with(b"#CHROM") {
            self.layout.set_freq();
            self.write_header(out)?;
        } else if !line.starts_with(b"#") {
            return self.record(line, out);
        }
        out.write_all(line)?;
        out.write_all(b"\n")
    }

    fn write_header(&self, out: &mut impl Write) -> io::Result<()> {
        const INFO: [(&str, &str, &str); 7] = [
            (
                "PGP_IMPACT",
                "1",
                "IMPACT of the ranked CSQ entry: MANE Select, else canonical, else any transcript",
            ),
            ("PGP_CSQ", "1", "Consequence of the ranked CSQ entry"),
            ("PGP_GENE", "1", "SYMBOL of the ranked CSQ entry"),
            (
                "PGP_TX",
                "1",
                "Transcript set of the ranked CSQ entry: MANE, CANONICAL or ANY",
            ),
            (
                "PGP_FREQ",
                "1",
                "Frequency tier: rare, common, unknown (no frequency for the allele) or none (no frequency field)",
            ),
            (
                "PGP_CALL",
                "1",
                "Genotype call check: ok, weak or unknown, against ##PGP_CALL_FLOOR",
            ),
            (
                "PGP_CALL_WHY",
                ".",
                "What made the call weak: GQ, DP or alt read fraction AB",
            ),
        ];
        for (id, number, what) in INFO {
            writeln!(
                out,
                "##INFO=<ID={id},Number={number},Type=String,Description=\"{what} (og vcf csq-rank)\">"
            )?;
        }
        let f = &self.floors;
        writeln!(
            out,
            "##PGP_CALL_FLOOR=GQ={},DP={},HET_AB={}-{},HOM_AB={}",
            f.min_gq.text, f.min_dp.text, f.het_lo.text, f.het_hi.text, f.hom_lo.text
        )
    }

    /// One record, or nothing when `--keep noncommon` drops it.
    fn record(&mut self, line: &[u8], out: &mut impl Write) -> io::Result<()> {
        // The first ten tabs end CHROM..the first sample; more are not needed.
        let mut tabs = [0usize; 10];
        let mut n = 0;
        for at in memchr_iter(b'\t', line) {
            if n == tabs.len() {
                break;
            }
            tabs[n] = at;
            n += 1;
        }
        let field = |k: usize| -> &[u8] {
            let start = if k == 0 { 0 } else { tabs[k - 1] + 1 };
            let end = if k < n { tabs[k] } else { line.len() };
            &line[start..end]
        };
        // awk's NF: an empty line has no field.
        let nf_at_least = |k: usize| !line.is_empty() && n + 1 >= k;

        // INFO without the old tags and `.`; the last CSQ= value.
        let info_in: &[u8] = if nf_at_least(8) { field(7) } else { b"" };
        self.info.clear();
        let mut csq: &[u8] = b"";
        for item in split(info_in, b';') {
            if let Some(v) = item.strip_prefix(b"CSQ=") {
                csq = v;
            }
            if item == b"." || is_old_tag(item) {
                continue;
            }
            if !self.info.is_empty() {
                self.info.push(b';');
            }
            self.info.extend_from_slice(item);
        }

        // The ranked entry: the strictly more severe one wins, so the first
        // of equals stays.
        let l = &self.layout;
        let (mut best_mane, mut best_canon, mut best_any) = (0u8, 0u8, 0u8);
        let (mut j_mane, mut j_canon, mut j_any) = (None, None, None);
        let mut known = false;
        let mut max_af = -1.0f64;
        for entry in split(csq, b',') {
            let e = Entry::new(entry, l.need);
            let r = severity(e.get(l.impact));
            if !e.get(l.mane).is_empty() && r >= best_mane {
                (best_mane, j_mane) = (r, Some(entry));
            }
            if e.get(l.canonical) == b"YES" && r >= best_canon {
                (best_canon, j_canon) = (r, Some(entry));
            }
            if r >= best_any {
                (best_any, j_any) = (r, Some(entry));
            }
            for &k in &l.freq {
                for v in split(e.get(Some(k)), b'&') {
                    if is_frequency(v) {
                        known = true;
                        max_af = max_af.max(number(v));
                    }
                }
            }
        }
        let (tx, ranked) = match (j_mane, j_canon, j_any) {
            (Some(e), _, _) => (TxSet::Mane, Some(e)),
            (None, Some(e), _) => (TxSet::Canonical, Some(e)),
            (None, None, Some(e)) => (TxSet::Any, Some(e)),
            (None, None, None) => (TxSet::None, None),
        };
        let (mut imp, mut cons, mut gene): (&[u8], &[u8], &[u8]) = (b"", b"", b"");
        if let Some(entry) = ranked {
            let e = Entry::new(entry, l.need);
            (imp, cons, gene) = (e.get(l.impact), e.get(l.consequence), e.get(l.symbol));
        }
        if imp.is_empty() {
            imp = b"NONE";
        }
        if cons.is_empty() {
            cons = b".";
        }
        if gene.is_empty() {
            gene = b".";
        }
        let freq = if l.freq.is_empty() {
            FreqTier::NoField
        } else if !known {
            FreqTier::Unknown
        } else if max_af < self.floors.rare_af {
            FreqTier::Rare
        } else {
            FreqTier::Common
        };
        if self.keep == Keep::NonCommon && freq == FreqTier::Common {
            return Ok(());
        }

        self.why.clear();
        let call = if nf_at_least(10) {
            self.judge(field(8), field(9))
        } else {
            "unknown"
        };

        // CHROM..FILTER, the new INFO, the rest.
        if nf_at_least(8) {
            out.write_all(&line[..tabs[6] + 1])?;
        } else if line.is_empty() {
            out.write_all(b"\t\t\t\t\t\t\t")?;
        } else {
            out.write_all(line)?;
            for _ in (n + 1)..8 {
                out.write_all(b"\t")?;
            }
        }
        if !self.info.is_empty() {
            out.write_all(&self.info)?;
            out.write_all(b";")?;
        }
        for (k, v) in [
            (&b"PGP_IMPACT="[..], imp),
            (b";PGP_CSQ=", cons),
            (b";PGP_GENE=", gene),
            (b";PGP_TX=", tx.as_str().as_bytes()),
            (b";PGP_FREQ=", freq.as_str().as_bytes()),
            (b";PGP_CALL=", call.as_bytes()),
        ] {
            out.write_all(k)?;
            out.write_all(v)?;
        }
        if !self.why.is_empty() {
            out.write_all(b";PGP_CALL_WHY=")?;
            out.write_all(&self.why)?;
        }
        if nf_at_least(9) {
            out.write_all(&line[tabs[7]..])?;
        }
        out.write_all(b"\n")
    }

    /// The call check of the first sample; fills `self.why`.
    fn judge(&mut self, format: &[u8], sample: &[u8]) -> &'static str {
        self.why.clear();
        let (mut gt, mut gq, mut dp, mut ad): (&[u8], &[u8], &[u8], &[u8]) = (b"", b"", b"", b"");
        let mut values = split(sample, b':');
        for key in split(format, b':') {
            let v = values.next().unwrap_or(b"");
            match key {
                b"GT" => gt = v,
                b"GQ" => gq = v,
                b"DP" => dp = v,
                b"AD" => ad = v,
                _ => {}
            }
        }
        let f = &self.floors;
        let mut judged = false;
        for (name, v, floor) in [(&b"GQ"[..], gq, &f.min_gq), (b"DP", dp, &f.min_dp)] {
            if is_digits(v) {
                judged = true;
                if number(v) < floor.value {
                    push_why(&mut self.why, name, v);
                }
            }
        }
        // The alt read fraction, on a diploid call of two numbered alleles.
        let mut alleles = split(gt, b'/').flat_map(|a| a.split(|&c| c == b'|'));
        let pair = match (alleles.next(), alleles.next(), alleles.next()) {
            (Some(a), Some(b), None) if is_digits(a) && is_digits(b) => {
                Some((number(a), number(b)))
            }
            _ => None,
        };
        let ad_ok = memchr(b',', ad).is_some() && split(ad, b',').all(is_digits);
        if let (Some((a, b)), true) = (pair, ad_ok) {
            let (mut total, mut alt) = (0.0f64, 0.0f64);
            for (i, v) in split(ad, b',').enumerate() {
                let x = number(v);
                total += x;
                if i > 0 {
                    alt += x;
                }
            }
            if total > 0.0 {
                judged = true;
                let ab = alt / total;
                let weak = if (a == 0.0) != (b == 0.0) {
                    ab < f.het_lo.value || ab > f.het_hi.value
                } else {
                    a == b && a != 0.0 && ab < f.hom_lo.value
                };
                if weak {
                    let text = pycompat::fmt::fixed(ab, 2);
                    push_why(&mut self.why, b"AB", text.as_bytes());
                }
            }
        }
        if !self.why.is_empty() {
            "weak"
        } else if judged {
            "ok"
        } else {
            "unknown"
        }
    }
}

fn push_why(why: &mut Vec<u8>, name: &[u8], value: &[u8]) {
    if !why.is_empty() {
        why.push(b',');
    }
    why.extend_from_slice(name);
    why.extend_from_slice(value);
}

const USAGE: &str =
    "usage: og vcf csq-rank [--keep all|noncommon] [--rare-af 0.01] [--min-gq 30] [--min-dp 10]
                       [--het-lo 0.2] [--het-hi 0.8] [--hom-lo 0.8] [--in FILE|-] [--out FILE|-]
";

/// `og vcf csq-rank ARGS`: the exit code. 0 done, 1 a read or write error,
/// 2 usage.
pub fn cli(args: &[String]) -> u8 {
    let usage = |msg: &str| {
        eprint!("{USAGE}");
        eprintln!("og vcf csq-rank: error: {msg}");
        2
    };
    let mut keep = Keep::All;
    let mut floors = Floors::default();
    let (mut input, mut output) = ("-".to_owned(), "-".to_owned());
    let mut it = args.iter();
    while let Some(a) = it.next() {
        if a == "-h" || a == "--help" {
            print!("{USAGE}");
            return 0;
        }
        let (key, inline) = match a.split_once('=') {
            Some((k, v)) if k.starts_with("--") => (k, Some(v.to_owned())),
            _ => (a.as_str(), None),
        };
        let known = [
            "--keep",
            "--rare-af",
            "--min-gq",
            "--min-dp",
            "--het-lo",
            "--het-hi",
            "--hom-lo",
            "--in",
            "--out",
        ];
        if !known.contains(&key) {
            return usage(&format!("unrecognized arguments: {a}"));
        }
        let Some(v) = inline.or_else(|| it.next().cloned()) else {
            return usage(&format!("argument {key}: expected one argument"));
        };
        match key {
            "--keep" => {
                keep = match v.as_str() {
                    "all" => Keep::All,
                    "noncommon" => Keep::NonCommon,
                    _ => {
                        return usage(&format!(
                            "argument --keep: invalid choice: '{v}' (choose from 'all', 'noncommon')"
                        ));
                    }
                }
            }
            "--in" => input = v,
            "--out" => output = v,
            _ => {
                let Some(f) = Floor::parse(&v) else {
                    return usage(&format!("argument {key}: not a number: '{v}'"));
                };
                match key {
                    "--rare-af" => floors.rare_af = f.value,
                    "--min-gq" => floors.min_gq = f,
                    "--min-dp" => floors.min_dp = f,
                    "--het-lo" => floors.het_lo = f,
                    "--het-hi" => floors.het_hi = f,
                    _ => floors.hom_lo = f,
                }
            }
        }
    }
    let reader = match og_io::open(&input) {
        Ok(r) => r,
        Err(e) => {
            eprintln!("og vcf csq-rank: {input}: {e}");
            return 1;
        }
    };
    let mut ranker = Ranker::new(keep, floors);
    let done = if output == "-" {
        let out = io::BufWriter::with_capacity(1 << 20, io::stdout().lock());
        ranker.run(reader, out)
    } else {
        match og_io::AtomicFile::create(&output) {
            Ok(mut f) => ranker
                .run(reader, &mut f)
                .and_then(|()| f.commit().map_err(Error)),
            Err(e) => Err(Error(e)),
        }
    };
    match done {
        Ok(()) => 0,
        Err(e) => {
            eprintln!("og vcf csq-rank: {e}");
            1
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const HEAD: &str = "##fileformat=VCFv4.2\n##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Feature|CANONICAL|MANE_SELECT|MAX_AF\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n";
    /// The seven records of the summary's test (tests/test_collect_summary.py).
    const BODY: &str = "chr2\t100\t.\tC\tT\t50\tPASS\tCSQ=T|stop_gained|HIGH|TTN|ENST2|||0.0001,T|intron_variant|MODIFIER|TTN|ENST1|YES|NM_1|0.0001\tGT:GQ:DP:AD\t0/1:50:30:15,15
chr13\t200\t.\tA\tAT\t50\tPASS\tCSQ=AT|frameshift_variant|HIGH|BRCA2|ENST3|YES|NM_2|0.0001\tGT:GQ:DP:AD\t0/1:12:.:3,2
chr19\t300\t.\tG\tA\t50\tPASS\tCSQ=A|stop_gained|HIGH|LDLR|ENST4|YES|NM_3|0.0001\tGT:GQ:DP:AD\t1/1:11:17:12,5
chr11\t400\t.\tG\tA\t50\tPASS\tCSQ=A|stop_gained|HIGH|MYBPC3|ENST5|YES|NM_4|0.0001\tGT:GQ:DP:AD\t1/1:50:30:0,30
chr14\t500\t.\tCTTT\tC\t50\tPASS\tCSQ=-|frameshift_variant|HIGH|MYH7|ENST6|YES|NM_5|\tGT:GQ:DP:AD\t0/1:50:30:15,15
chr11\t600\t.\tC\tT\t50\tPASS\tCSQ=T|stop_gained|HIGH|KCNQ1|ENST7|YES|NM_6|0.3\tGT:GQ:DP:AD\t0/1:50:30:15,15
chr3\t700\t.\tC\tT\t50\tPASS\tCSQ=T|missense_variant|MODERATE|SCN5A|ENST8|||0.0001,T|stop_gained|HIGH|SCN5A|ENST9|YES||0.0001\tGT:GQ:DP:AD\t0/1:50:30:15,15
";

    fn rank(keep: Keep, floors: Floors, vcf: &str) -> String {
        let mut out = Vec::new();
        Ranker::new(keep, floors)
            .run(vcf.as_bytes(), &mut out)
            .unwrap();
        String::from_utf8(out).unwrap()
    }

    /// The PGP_ tags of each record, by POS.
    fn tags(out: &str) -> Vec<(String, String)> {
        out.lines()
            .filter(|l| !l.starts_with('#'))
            .map(|l| {
                let f: Vec<&str> = l.split('\t').collect();
                let t: Vec<&str> = f[7].split(';').filter(|x| x.starts_with("PGP_")).collect();
                (f[1].to_owned(), t.join(";"))
            })
            .collect()
    }

    #[test]
    fn mane_then_canonical_frequency_and_call() {
        let out = rank(Keep::All, Floors::default(), &format!("{HEAD}{BODY}"));
        let want = [
            (
                "100",
                "PGP_IMPACT=MODIFIER;PGP_CSQ=intron_variant;PGP_GENE=TTN;PGP_TX=MANE;PGP_FREQ=rare;PGP_CALL=ok",
            ),
            (
                "200",
                "PGP_IMPACT=HIGH;PGP_CSQ=frameshift_variant;PGP_GENE=BRCA2;PGP_TX=MANE;PGP_FREQ=rare;PGP_CALL=weak;PGP_CALL_WHY=GQ12",
            ),
            (
                "300",
                "PGP_IMPACT=HIGH;PGP_CSQ=stop_gained;PGP_GENE=LDLR;PGP_TX=MANE;PGP_FREQ=rare;PGP_CALL=weak;PGP_CALL_WHY=GQ11,AB0.29",
            ),
            (
                "400",
                "PGP_IMPACT=HIGH;PGP_CSQ=stop_gained;PGP_GENE=MYBPC3;PGP_TX=MANE;PGP_FREQ=rare;PGP_CALL=ok",
            ),
            (
                "500",
                "PGP_IMPACT=HIGH;PGP_CSQ=frameshift_variant;PGP_GENE=MYH7;PGP_TX=MANE;PGP_FREQ=unknown;PGP_CALL=ok",
            ),
            (
                "600",
                "PGP_IMPACT=HIGH;PGP_CSQ=stop_gained;PGP_GENE=KCNQ1;PGP_TX=MANE;PGP_FREQ=common;PGP_CALL=ok",
            ),
            (
                "700",
                "PGP_IMPACT=HIGH;PGP_CSQ=stop_gained;PGP_GENE=SCN5A;PGP_TX=CANONICAL;PGP_FREQ=rare;PGP_CALL=ok",
            ),
        ];
        let got = tags(&out);
        let got: Vec<(&str, &str)> = got.iter().map(|(a, b)| (a.as_str(), b.as_str())).collect();
        assert_eq!(got, want);
        assert!(out.contains("##PGP_CALL_FLOOR=GQ=30,DP=10,HET_AB=0.2-0.8,HOM_AB=0.8\n#CHROM\t"));
        assert_eq!(out.matches("##INFO=<ID=PGP_").count(), 7);
        // The rest of a record is untouched.
        assert!(out.contains("\tGT:GQ:DP:AD\t0/1:12:.:3,2\n"));
    }

    #[test]
    fn noncommon_drops_only_the_common_record() {
        let out = rank(Keep::NonCommon, Floors::default(), &format!("{HEAD}{BODY}"));
        let pos: Vec<String> = tags(&out).into_iter().map(|(p, _)| p).collect();
        assert_eq!(pos, ["100", "200", "300", "400", "500", "700"]);
    }

    #[test]
    fn a_second_pass_changes_nothing() {
        let once = rank(Keep::All, Floors::default(), &format!("{HEAD}{BODY}"));
        assert_eq!(rank(Keep::All, Floors::default(), &once), once);
    }

    #[test]
    fn floors_are_echoed_as_given_and_applied() {
        let f = |t: &str| Floor::parse(t).unwrap();
        let floors = Floors {
            rare_af: 0.5,
            min_gq: f("0"),
            min_dp: f("0"),
            het_lo: f("0"),
            het_hi: f("1.0"),
            hom_lo: f("0"),
        };
        let out = rank(Keep::NonCommon, floors, &format!("{HEAD}{BODY}"));
        assert!(out.contains("##PGP_CALL_FLOOR=GQ=0,DP=0,HET_AB=0-1.0,HOM_AB=0\n"));
        assert!(!out.contains("PGP_CALL=weak"));
        assert!(out.contains("PGP_GENE=KCNQ1;PGP_TX=MANE;PGP_FREQ=rare"));
    }

    #[test]
    fn the_first_of_equally_severe_entries_wins() {
        let rec = "chr1\t1\t.\tA\tG\t50\tPASS\tCSQ=G|missense_variant|HIGH|FIRST|E1|YES||,G|stop_gained|HIGH|SECOND|E2|YES||\n";
        let out = rank(Keep::All, Floors::default(), &format!("{HEAD}{rec}"));
        assert!(out.contains("PGP_CSQ=missense_variant;PGP_GENE=FIRST;PGP_TX=CANONICAL;PGP_FREQ=unknown;PGP_CALL=unknown\n"), "{out}");
    }

    #[test]
    fn short_lines_and_odd_info() {
        // No CSQ, `.` INFO, fewer than 8 columns, an empty line.
        let body = "chr1\t1\t.\tA\tG\t50\tPASS\t.\nchr1\t2\nx\n\n";
        let out = rank(Keep::All, Floors::default(), &format!("{HEAD}{body}"));
        let tail =
            "PGP_IMPACT=NONE;PGP_CSQ=.;PGP_GENE=.;PGP_TX=NONE;PGP_FREQ=unknown;PGP_CALL=unknown";
        let records: Vec<&str> = out.lines().filter(|l| !l.starts_with('#')).collect();
        assert_eq!(
            records,
            [
                format!("chr1\t1\t.\tA\tG\t50\tPASS\t{tail}"),
                format!("chr1\t2\t\t\t\t\t\t{tail}"),
                format!("x\t\t\t\t\t\t\t{tail}"),
                format!("\t\t\t\t\t\t\t{tail}"),
            ]
        );
    }

    #[test]
    fn frequencies_and_floors_parse_strictly() {
        for (s, ok) in [
            ("0.5", true),
            (".5", true),
            ("1e-3", true),
            ("1E+2", true),
            ("5", true),
            ("5.", false),
            ("-0.1", false),
            ("nan", false),
            ("inf", false),
            ("", false),
            ("1e", false),
            (" 1", false),
        ] {
            assert_eq!(is_frequency(s.as_bytes()), ok, "{s:?}");
        }
        for (s, ok) in [
            ("30", true),
            ("0.2", true),
            ("5.", true),
            (".5", true),
            ("+1", true),
            ("1e2", true),
            ("abc", false),
            ("", false),
            ("nan", false),
            ("inf", false),
            ("1e999", false),
            ("1.2.3", false),
            ("0x10", false),
        ] {
            assert_eq!(Floor::parse(s).is_some(), ok, "{s:?}");
        }
    }
}
