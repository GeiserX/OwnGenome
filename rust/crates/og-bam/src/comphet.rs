//! `og comphet phase` (was `bin/comphet_phase.py`): are the two variants of
//! a compound-het candidate pair on the same copy of the gene (cis) or on
//! different copies (trans)?
//!
//! slivar compound-hets lists every pair of heterozygous candidates in a
//! gene. From one unphased genome a pair in cis looks the same as a pair in
//! trans, unless one read, or one read pair, covers both positions: in cis a
//! fragment shows both alt alleles or neither, in trans one of them. For each
//! pair this counts the fragments that show an allele at both positions and
//! calls trans or cis when at least [`MIN_READS`] do and at least
//! [`MIN_SHARE`] of them agree; otherwise unknown, with the reason.
//!
//! Reads that count: mapped, primary, not duplicates, not QC-failed, mapping
//! quality at least 20 ([`Filter::PHASE`]). The two mates of a pair are one
//! fragment; mates that disagree at a position are left out there. Each
//! variant's reads are fetched once.

use std::collections::{BTreeSet, HashMap};
use std::io::Write;
use std::path::Path;

use indexmap::IndexMap;
use og_io::{AtomicFile, Lines, Utf8};
use pycompat::num::parse_int;
use rustc_hash::FxHashMap;

use crate::Bam;
use crate::args::Args;
use crate::reads::{Filter, Seen, Site, allele_at, complex_allele};

/// The fewest fragments a verdict rests on.
pub const MIN_READS: u64 = 2;
/// The share of fragments that must agree.
pub const MIN_SHARE: f64 = 0.8;
/// The columns of the pairs table.
pub const COLUMNS: [&str; 13] = [
    "GENE",
    "CHROM",
    "POS_1",
    "REF_1",
    "ALT_1",
    "POS_2",
    "REF_2",
    "ALT_2",
    "DISTANCE",
    "PHASE",
    "CIS_READS",
    "TRANS_READS",
    "PHASE_NOTE",
];

/// One variant: contig, 1-based position, REF and ALT. Ordered as Python
/// orders the tuple: contig by code point, then position.
pub type Var = (String, i64, String, String);

/// One pair: its gene and its two variants, the lower first.
pub type Pair = (String, Var, Var);

fn int(s: &str) -> Result<i64, String> {
    parse_int(s)
        .ok()
        .and_then(|v| i64::try_from(v).ok())
        .ok_or_else(|| format!("invalid literal for int() with base 10: '{s}'"))
}

/// The pairs of slivar's VCF (plain or gzip), sorted. Each record's
/// INFO/slivar_comphet entries (`sample/gene/pair_id/chrom/pos/ref/alt`,
/// the gene may hold `/`) name its partner; the two records of a pair share
/// its id. A pair id that names other than two variants is skipped with a
/// warning on `warn`.
pub fn read_pairs(path: &str, warn: &mut dyn Write) -> Result<Vec<Pair>, String> {
    let fh = og_io::open(path).map_err(|e| format!("{path}: {e}"))?;
    let mut groups: IndexMap<(String, String), BTreeSet<Var>> = IndexMap::new();
    for line in Lines::new(fh, Utf8::Strict) {
        let line = line.map_err(|e| format!("{path}: {e}"))?;
        if line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 8 {
            continue;
        }
        let own: Var = (
            f[0].to_owned(),
            int(f[1])?,
            f[3].to_owned(),
            f[4].to_owned(),
        );
        // dict(): the last slivar_comphet key wins.
        let entries = f[7]
            .split(';')
            .rev()
            .map(|kv| kv.split_once('=').unwrap_or((kv, "")))
            .find(|(k, _)| *k == "slivar_comphet")
            .map_or("", |(_, v)| v);
        for entry in entries.split(',').filter(|e| !e.is_empty()) {
            let p: Vec<&str> = entry.split('/').collect();
            let n = p.len();
            if n < 7 {
                return Err(format!("unreadable slivar_comphet entry '{entry}'"));
            }
            let gene = p[1..n - 5].join("/");
            let other: Var = (
                p[n - 4].to_owned(),
                int(p[n - 3])?,
                p[n - 2].to_owned(),
                p[n - 1].to_owned(),
            );
            let set = groups.entry((gene, p[n - 5].to_owned())).or_default();
            set.insert(own.clone());
            set.insert(other);
        }
    }
    let mut pairs = BTreeSet::new();
    for ((gene, pid), vs) in groups {
        if vs.len() != 2 {
            let _ = writeln!(
                warn,
                "WARNING: pair {pid} in {gene} names {} variants, not 2; skipped",
                vs.len()
            );
            continue;
        }
        let mut it = vs.into_iter();
        let (a, b) = (it.next().expect("two"), it.next().expect("two"));
        pairs.insert((gene, a, b));
    }
    Ok(pairs.into_iter().collect())
}

/// The allele each fragment shows at a variant, by fragment name: a read
/// pair is one fragment (`join_mates`), or each read alone. A fragment whose
/// reads disagree is left out.
pub fn site_alleles<R: crate::reads::ReadView>(
    reads: impl IntoIterator<Item = Result<R, String>>,
    site: &Site,
    join_mates: bool,
) -> Result<FxHashMap<(Vec<u8>, bool), Seen>, String> {
    let mut out: FxHashMap<(Vec<u8>, bool), Seen> = FxHashMap::default();
    let mut bad = Vec::new();
    for r in reads {
        let r = r?;
        if !Filter::PHASE.keeps(&r) {
            continue;
        }
        let Some(a) = allele_at(&r, site) else {
            continue;
        };
        let key = (
            r.qname().to_vec(),
            !join_mates && r.has(crate::reads::flag::READ1),
        );
        if out.insert(key.clone(), a).is_some_and(|before| before != a) {
            bad.push(key);
        }
    }
    for k in bad {
        out.remove(&k);
    }
    Ok(out)
}

/// (cis, trans): the fragments that show both positions.
pub fn phase_counts<K: std::hash::Hash + Eq>(
    a1: &FxHashMap<K, Seen>,
    a2: &FxHashMap<K, Seen>,
) -> (u64, u64) {
    let (mut cis, mut trans) = (0, 0);
    for (k, x) in a1 {
        match a2.get(k) {
            Some(y) if y == x => cis += 1,
            Some(_) => trans += 1,
            None => {}
        }
    }
    (cis, trans)
}

/// (PHASE, PHASE_NOTE) from the two counts.
pub fn verdict(cis: u64, trans: u64) -> (&'static str, String) {
    let n = cis + trans;
    #[allow(clippy::cast_precision_loss)] // read counts are far below 2^52
    let share = |k: u64| k as f64 >= MIN_SHARE * n as f64;
    if n == 0 {
        (
            "unknown",
            "no read or read pair covers both positions".to_owned(),
        )
    } else if n < MIN_READS {
        (
            "unknown",
            format!("only {n} read or read pair covers both positions"),
        )
    } else if share(trans) {
        (
            "trans",
            format!("{trans} of {n} reads or read pairs show the two variants on different copies"),
        )
    } else if share(cis) {
        (
            "cis",
            format!("{cis} of {n} reads or read pairs show the two variants on the same copy"),
        )
    } else {
        (
            "unknown",
            format!("the reads disagree: {cis} cis, {trans} trans"),
        )
    }
}

type Fragments = FxHashMap<(Vec<u8>, bool), Seen>;

/// The fragments at each variant, each variant read once from the BAM.
pub struct Reads {
    bam: Bam,
    join_mates: bool,
    cache: HashMap<Var, Fragments>,
}

impl Reads {
    pub fn open(path: &str, join_mates: bool) -> Result<Reads, String> {
        Ok(Reads {
            bam: Bam::open(Path::new(path)).map_err(|e| format!("{path}: {e}"))?,
            join_mates,
            cache: HashMap::new(),
        })
    }

    /// The BAM's name for `chrom`: itself, else with `chr` removed or added.
    pub fn contig(&self, chrom: &str) -> Option<String> {
        let other = match chrom.strip_prefix("chr") {
            Some(bare) => bare.to_owned(),
            None => format!("chr{chrom}"),
        };
        [chrom.to_owned(), other]
            .into_iter()
            .find(|c| self.bam.tid(c).is_some())
    }

    fn at(&mut self, var: &Var) -> Result<&Fragments, String> {
        if !self.cache.contains_key(var) {
            let (chrom, pos, reference, alt) = var;
            let contig = self
                .contig(chrom)
                .ok_or_else(|| format!("the BAM has no contig {chrom}"))?;
            let site = Site::new(pos - 1, reference, alt);
            let (lo, hi) = site.window();
            if lo < 0 {
                return Err(format!("start out of range ({lo})"));
            }
            let join = self.join_mates;
            let reads = self
                .bam
                .fetch(&contig, lo, hi)
                .map_err(|e| format!("{contig}:{lo}-{hi}: {e}"))?
                .map(|r| r.map_err(|e| e.to_string()));
            let got = site_alleles(reads, &site, join)?;
            self.cache.insert(var.clone(), got);
        }
        Ok(&self.cache[var])
    }
}

/// One row of the pairs table.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Row {
    pub gene: String,
    pub v1: Var,
    pub v2: Var,
    pub phase: &'static str,
    /// (cis, trans), `None` when the pair was not read.
    pub counts: Option<(u64, u64)>,
    pub note: String,
}

impl Row {
    fn cells(&self) -> [String; 13] {
        let (v1, v2) = (&self.v1, &self.v2);
        let same = v1.0 == v2.0;
        let n = |f: fn((u64, u64)) -> u64| self.counts.map_or(".".to_owned(), |c| f(c).to_string());
        [
            self.gene.clone(),
            v1.0.clone(),
            v1.1.to_string(),
            v1.2.clone(),
            v1.3.clone(),
            v2.1.to_string(),
            v2.2.clone(),
            v2.3.clone(),
            if same {
                (v2.1 - v1.1).abs().to_string()
            } else {
                ".".to_owned()
            },
            self.phase.to_owned(),
            n(|c| c.0),
            n(|c| c.1),
            self.note.clone(),
        ]
    }
}

/// The phase of each pair: from `reads`, or every pair unknown with `skip`.
pub fn phase(
    pairs: Vec<Pair>,
    mut reads: Option<&mut Reads>,
    skip: Option<&str>,
) -> Result<Vec<Row>, String> {
    let mut rows = Vec::with_capacity(pairs.len());
    for (gene, v1, v2) in pairs {
        let unread = |note: String| ("unknown", None, note);
        let (phase, counts, note) = match (skip, reads.as_deref_mut()) {
            (Some(why), _) => unread(why.to_owned()),
            (None, None) => unread(String::new()),
            (None, Some(r)) => {
                if r.contig(&v1.0).is_none() {
                    unread(format!("the BAM has no contig {}", v1.0))
                } else if v1.0 != v2.0 {
                    unread("the two variants are on different contigs".to_owned())
                } else if complex_allele(&v1.2, &v1.3) || complex_allele(&v2.2, &v2.3) {
                    unread("a complex allele: not phased from the reads".to_owned())
                } else {
                    let a1 = r.at(&v1)?.clone();
                    let (cis, trans) = phase_counts(&a1, r.at(&v2)?);
                    let (p, note) = verdict(cis, trans);
                    (p, Some((cis, trans)), note)
                }
            }
        };
        rows.push(Row {
            gene,
            v1,
            v2,
            phase,
            counts,
            note,
        });
    }
    Ok(rows)
}

/// Each variant's PHASE over its pairs: trans when it is in trans with one
/// partner, else unknown when one pair is unknown, else cis. Keys are the
/// cells as the table writes them.
pub fn variant_phase(rows: &[Row]) -> HashMap<[String; 4], &'static str> {
    let mut seen: HashMap<[String; 4], Vec<&'static str>> = HashMap::new();
    for r in rows {
        for v in [&r.v1, &r.v2] {
            let key = [r.v1.0.clone(), v.1.to_string(), v.2.clone(), v.3.clone()];
            seen.entry(key).or_default().push(r.phase);
        }
    }
    seen.into_iter()
        .map(|(k, s)| {
            let best = ["trans", "unknown", "cis"]
                .into_iter()
                .find(|x| s.contains(x))
                .unwrap_or("unknown");
            (k, best)
        })
        .collect()
}

/// Step 31's per-variant table with a PHASE column after GT (last when it
/// has no GT), written atomically to `dst`.
pub fn add_phase_column(src: &str, dst: &str, rows: &[Row]) -> Result<(), String> {
    let best = variant_phase(rows);
    let fh = std::fs::File::open(src).map_err(|e| format!("{src}: {e}"))?;
    let mut lines = Lines::new(std::io::BufReader::new(fh), Utf8::Strict);
    let head_line = lines
        .next()
        .transpose()
        .map_err(|e| format!("{src}: {e}"))?;
    let head: Vec<String> = head_line
        .unwrap_or_default()
        .split('\t')
        .map(str::to_owned)
        .collect();
    let mut col = [0usize; 4];
    for (i, k) in ["CHROM", "POS", "REF", "ALT"].iter().enumerate() {
        col[i] = head
            .iter()
            .position(|h| h == k)
            .ok_or_else(|| format!("{src} has no {k} column"))?;
    }
    let at = head
        .iter()
        .position(|h| h == "GT")
        .map_or(head.len(), |i| i + 1);
    let fail = |e: std::io::Error| format!("cannot write {dst}: {e}");
    let mut out = AtomicFile::create(dst).map_err(fail)?;
    let mut h = head.clone();
    h.insert(at, "PHASE".to_owned());
    writeln!(out, "{}", h.join("\t")).map_err(fail)?;
    for (n, line) in lines.enumerate() {
        let line = line.map_err(|e| format!("{src}: {e}"))?;
        let mut f: Vec<&str> = line.split('\t').collect();
        let key = col
            .iter()
            .map(|&i| f.get(i).map(|s| (*s).to_owned()))
            .collect::<Option<Vec<String>>>()
            .ok_or_else(|| format!("{src}: line {} has fewer columns than its header", n + 2))?;
        let key: [String; 4] = key.try_into().expect("four columns");
        let phase = best.get(&key).copied().unwrap_or(".");
        f.insert(at.min(f.len()), phase);
        writeln!(out, "{}", f.join("\t")).map_err(fail)?;
    }
    out.commit().map_err(fail)
}

fn write_rows(path: &str, rows: &[Row]) -> Result<(), String> {
    let fail = |e: std::io::Error| format!("cannot write {path}: {e}");
    let mut f = AtomicFile::create(path).map_err(fail)?;
    writeln!(f, "{}", COLUMNS.join("\t")).map_err(fail)?;
    for r in rows {
        writeln!(f, "{}", r.cells().join("\t")).map_err(fail)?;
    }
    f.commit().map_err(fail)
}

/// `og comphet phase ARGS`: the exit code; a refused command line is exit 2.
pub fn cli(argv: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    const USAGE: &str = "usage: og comphet phase --vcf VCF --out OUT [--bam BAM] [--skip SKIP] \
                         [--variants VARIANTS] [--variants-out VARIANTS_OUT]";
    let refuse = |err: &mut dyn Write, e: &str| {
        let _ = writeln!(err, "{USAGE}\nog: error: {e}");
        2
    };
    let a = match Args::parse(
        argv,
        &[
            "--vcf",
            "--out",
            "--bam",
            "--skip",
            "--variants",
            "--variants-out",
        ],
        &[],
    ) {
        Ok(a) => a,
        Err(e) => return refuse(err, &e),
    };
    let missing: Vec<&str> = ["--vcf", "--out"]
        .into_iter()
        .filter(|k| a.get(k).is_none())
        .collect();
    if !missing.is_empty() {
        return refuse(
            err,
            &format!(
                "the following arguments are required: {}",
                missing.join(", ")
            ),
        );
    }
    let given = |k: &str| a.get(k).filter(|s| !s.is_empty());
    let (vcf, path) = (
        a.get("--vcf").unwrap_or_default(),
        a.get("--out").unwrap_or_default(),
    );
    let (bam, skip) = (given("--bam"), given("--skip"));
    if bam.is_some() == skip.is_some() {
        return refuse(err, "give --bam or --skip REASON");
    }
    let variants = given("--variants").zip(given("--variants-out"));
    if given("--variants").is_some() != given("--variants-out").is_some() {
        return refuse(err, "--variants and --variants-out go together");
    }
    let mut run = || -> Result<Vec<Row>, String> {
        let pairs = read_pairs(vcf, err)?;
        let mut reads = bam.map(|b| Reads::open(b, true)).transpose()?;
        let rows = phase(pairs, reads.as_mut(), skip)?;
        if let Some((src, dst)) = variants {
            add_phase_column(src, dst, &rows)?;
        }
        write_rows(path, &rows)?;
        Ok(rows)
    };
    let rows = match run() {
        Ok(r) => r,
        Err(e) => {
            let _ = writeln!(err, "ERROR: {e}");
            return 2;
        }
    };
    let n = |k: &str| rows.iter().filter(|r| r.phase == k).count();
    let tail = skip.map_or(String::new(), |s| format!(" ({s})"));
    match writeln!(
        out,
        "Phase of {} compound-het pairs from the reads: {} trans, {} cis, {} unknown{tail}",
        rows.len(),
        n("trans"),
        n("cis"),
        n("unknown")
    ) {
        Ok(()) => 0,
        Err(e) => {
            let _ = writeln!(err, "ERROR: {e}");
            2
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reads::SamRead;
    use crate::reads::flag;
    use crate::reads::tests::read;

    fn mate(name: &str, read1: bool, seq: &[u8]) -> SamRead {
        let mut r = read(name, 99, "3M", seq);
        r.flags = if read1 { flag::READ1 } else { flag::READ2 };
        r
    }

    fn alleles(reads: Vec<SamRead>, join: bool) -> Vec<(String, Seen)> {
        let site = Site::new(100, "C", "G");
        let mut v: Vec<_> = site_alleles(reads.into_iter().map(Ok), &site, join)
            .unwrap()
            .into_iter()
            .map(|((n, _), s)| (String::from_utf8(n).unwrap(), s))
            .collect();
        v.sort_by(|a, b| a.0.cmp(&b.0));
        v
    }

    #[test]
    fn mates_that_disagree_are_left_out() {
        let got = alleles(
            vec![
                mate("a", true, b"AGA"),
                mate("a", false, b"ACA"),
                mate("b", true, b"AGA"),
            ],
            true,
        );
        assert_eq!(got, vec![("b".to_owned(), Seen::Alt)]);
        // Each read alone: a's two reads are two fragments.
        let got = alleles(
            vec![mate("a", true, b"AGA"), mate("a", false, b"ACA")],
            false,
        );
        assert_eq!(got.len(), 2);
        // No stored sequence: skipped, the rest still counted.
        let got = alleles(vec![mate("a", true, b""), mate("b", true, b"AGA")], true);
        assert_eq!(got, vec![("b".to_owned(), Seen::Alt)]);
    }

    #[test]
    fn duplicates_supplementary_and_low_mapq_reads_do_not_count() {
        let mut d = mate("d", true, b"AGA");
        d.flags |= flag::DUPLICATE;
        let mut s = mate("s", true, b"AGA");
        s.flags |= flag::SUPPLEMENTARY;
        let mut q = mate("q", true, b"AGA");
        q.mapq = 19;
        let mut ok = mate("ok", true, b"AGA");
        ok.mapq = 20;
        assert_eq!(
            alleles(vec![d, s, q, ok], true),
            vec![("ok".to_owned(), Seen::Alt)]
        );
    }

    #[test]
    fn the_verdicts() {
        assert_eq!(verdict(0, 0).0, "unknown");
        assert_eq!(
            verdict(0, 1),
            (
                "unknown",
                "only 1 read or read pair covers both positions".to_owned()
            )
        );
        assert_eq!(verdict(2, 8).0, "trans");
        assert_eq!(
            verdict(3, 7),
            ("unknown", "the reads disagree: 3 cis, 7 trans".to_owned())
        );
        assert_eq!(verdict(9, 1).0, "cis");
        assert_eq!(
            verdict(22, 0),
            (
                "cis",
                "22 of 22 reads or read pairs show the two variants on the same copy".to_owned()
            )
        );
    }

    #[test]
    fn counts_cis_and_trans_over_shared_fragments() {
        let m = |v: &[(&str, Seen)]| -> FxHashMap<(Vec<u8>, bool), Seen> {
            v.iter()
                .map(|(n, s)| ((n.as_bytes().to_vec(), false), *s))
                .collect()
        };
        let a = m(&[("x", Seen::Alt), ("y", Seen::Ref), ("z", Seen::Alt)]);
        let b = m(&[("x", Seen::Alt), ("y", Seen::Alt), ("w", Seen::Ref)]);
        assert_eq!(phase_counts(&a, &b), (1, 1));
    }

    fn row(v1: Var, v2: Var, phase: &'static str) -> Row {
        Row {
            gene: "G".into(),
            v1,
            v2,
            phase,
            counts: None,
            note: String::new(),
        }
    }

    fn v(pos: i64, r: &str, a: &str) -> Var {
        ("c".into(), pos, r.into(), a.into())
    }

    #[test]
    fn per_variant_phase_trans_over_unknown_over_cis() {
        let rows = [
            row(v(1, "A", "C"), v(5, "G", "T"), "cis"),
            row(v(1, "A", "C"), v(9, "T", "A"), "trans"),
            row(v(5, "G", "T"), v(9, "T", "A"), "unknown"),
        ];
        let got = variant_phase(&rows);
        let k =
            |p: &str, r: &str, a: &str| ["c".to_owned(), p.to_owned(), r.to_owned(), a.to_owned()];
        assert_eq!(got.len(), 3);
        assert_eq!(got[&k("1", "A", "C")], "trans");
        assert_eq!(got[&k("5", "G", "T")], "unknown");
        assert_eq!(got[&k("9", "T", "A")], "trans");
    }

    #[test]
    fn slivar_pairs_are_read_once_sorted_and_genes_keep_their_slash() {
        let dir = std::env::temp_dir().join(format!("og-comphet-pairs-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("c.vcf");
        let info = |e: &str| format!("X=1;slivar_comphet={e}");
        let rec = |pos: u32, r: &str, a: &str, e: &str| {
            format!("chr1\t{pos}\t.\t{r}\t{a}\t.\tPASS\t{}\n", info(e))
        };
        let text = [
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n".to_owned(),
            rec(500, "G", "T", "s/HLA/A/1/chr1/100/A/C"),
            rec(100, "A", "C", "s/HLA/A/1/chr1/500/G/T,s/B/2/chr1/900/T/A"),
            rec(900, "T", "A", "s/B/2/chr1/100/A/C"),
            "chr1\t7\t.\tA\tC\n".to_owned(),
            rec(950, "T", "A", "s/B/3/chr1/951/T/A,s/B/3/chr1/952/T/A"),
        ]
        .concat();
        std::fs::write(&p, text).unwrap();
        let mut warn = Vec::new();
        let got = read_pairs(p.to_str().unwrap(), &mut warn).unwrap();
        std::fs::remove_dir_all(&dir).unwrap();
        assert_eq!(
            got,
            vec![
                (
                    "B".into(),
                    ("chr1".into(), 100, "A".into(), "C".into()),
                    ("chr1".into(), 900, "T".into(), "A".into())
                ),
                (
                    "HLA/A".into(),
                    ("chr1".into(), 100, "A".into(), "C".into()),
                    ("chr1".into(), 500, "G".into(), "T".into())
                ),
            ]
        );
        assert_eq!(
            String::from_utf8(warn).unwrap(),
            "WARNING: pair 3 in B names 3 variants, not 2; skipped\n"
        );
    }
}
