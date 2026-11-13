//! `og family report`: per child with a sequenced parent, the allele each
//! parent passed on at every site of the joint VCF where the genotypes
//! tell (`<child>_inheritance.tsv.gz`), de novo candidates when both
//! parents are sequenced (`de_novo.tsv`), pairs of candidate variants in
//! one gene in trans or in cis (`compound_het.tsv`), the chromosome map
//! (`<child>_segments.tsv`, `<child>_dnapainter.csv`, `<child>_segments.svg`)
//! and `family_report.md`, which says it in words.
//!
//! The joint VCF is read once: the inheritance and the map take their
//! genotypes from the same pass. With `--phased` the inheritance reads
//! WhatsHap's file (it needs the phase) while the map reads `--vcf` on a
//! second thread. Each child's inheritance file is compressed on a thread
//! of its own.

use std::collections::HashMap;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::mpsc::{SyncSender, sync_channel};
use std::thread::{self, JoinHandle};

use flate2::{Compression, GzBuilder};
use indexmap::IndexMap;
use og_io::AtomicFile;
use og_pedigree::{Die, Members, Par, cell, read_members, read_table};
use pycompat::fmt::fixed;

use super::args::{Args, Kind, Spec};
use super::segments::{self, Collector, Map};
use super::vcf::{self, Calls, resolve, set};
use super::{Result, print_out, write_file};

pub const SPEC: Spec = Spec {
    prog: "og family report",
    usage: "usage: og family report [-h] --members MEMBERS --relationships RELATIONSHIPS --vcf VCF [--phased PHASED]\n                        [--candidates ID=VCF] [--par-bed PAR_BED] --out-dir OUT_DIR\n",
    options: &[
        ("--members", Kind::One),
        ("--relationships", Kind::One),
        ("--vcf", Kind::One),
        ("--phased", Kind::One),
        ("--candidates", Kind::Append),
        ("--par-bed", Kind::One),
        ("--out-dir", Kind::One),
    ],
    required: &["--members", "--relationships", "--vcf", "--out-dir"],
};

/// De novo candidates: every one of the three must have these.
const MIN_GQ: i128 = 20;
const MIN_DP: i128 = 10;
const HET_AB: (f64, f64) = (0.25, 0.75);
/// Bytes of inheritance rows handed to a child's compressing thread at a time.
const CHUNK: usize = 1 << 18;
/// The level of the inheritance files (the content is what is compared).
const GZ_LEVEL: u32 = 6;

fn is_primary(c: &[u8]) -> bool {
    c == b"chrY" || segments::chrom_index(c).is_some()
}

type Genes = std::rc::Rc<Vec<String>>;
type CandKey = (Vec<u8>, i64, Vec<u8>, Vec<u8>);

/// `{(chrom, pos, ref, alt): genes}` of a candidate VCF: VEP's CSQ SYMBOL,
/// else snpEff's ANN gene, else INFO/GENE.
fn candidate_genes(path: &str) -> std::result::Result<HashMap<CandKey, Genes>, Die> {
    let r = og_io::open(path).map_err(|e| Die::io(path, &e))?;
    let mut lines = og_vcf::Lines::new(r);
    let (mut symbol, mut ann) = (None::<usize>, false);
    let mut out = HashMap::new();
    let mut number = 0;
    while let Some(line) = lines.next_line().map_err(|e| Die::io(path, &e))? {
        number += 1;
        if line.starts_with(b"##INFO=<ID=CSQ") {
            if let Some(k) = memchr::memmem::find(line, b"Format: ") {
                let rest = &line[k + 8..];
                let end = rest.iter().position(|&c| c == b'"').unwrap_or(rest.len());
                if end > 0 {
                    symbol = rest[..end]
                        .split(|&c| c == b'|')
                        .position(|f| f == b"SYMBOL");
                }
            }
            continue;
        }
        if line.starts_with(b"##INFO=<ID=ANN") {
            ann = true;
            continue;
        }
        if line.starts_with(b"#") {
            continue;
        }
        let c = vcf::columns(line);
        if c.len() < 8 {
            continue;
        }
        let mut info: HashMap<&[u8], &[u8]> = HashMap::new();
        for kv in c[7].split(|&b| b == b';') {
            match kv.iter().position(|&b| b == b'=') {
                Some(k) => info.insert(&kv[..k], &kv[k + 1..]),
                None => info.insert(kv, b""),
            };
        }
        let mut genes: Vec<String> = Vec::new();
        let mut add = |g: &[u8]| {
            let g = String::from_utf8_lossy(g).into_owned();
            if !genes.contains(&g) {
                genes.push(g);
            }
        };
        if let (Some(si), Some(csq)) = (symbol, info.get(&b"CSQ"[..])) {
            for entry in csq.split(|&b| b == b',') {
                if let Some(g) = entry
                    .split(|&b| b == b'|')
                    .nth(si)
                    .filter(|g| !g.is_empty())
                {
                    add(g);
                }
            }
        } else if let (true, Some(a)) = (ann, info.get(&b"ANN"[..])) {
            for entry in a.split(|&b| b == b',') {
                if let Some(g) = entry.split(|&b| b == b'|').nth(3).filter(|g| !g.is_empty()) {
                    add(g);
                }
            }
        } else if let Some(g) = info.get(&b"GENE"[..]).filter(|g| !g.is_empty()) {
            for x in g.split(|&b| b == b',').filter(|x| !x.is_empty()) {
                add(x);
            }
        }
        if !genes.is_empty() {
            let pos = vcf::pos(c[1], path, number)?;
            let genes = Genes::new(genes);
            for alt in c[4].split(|&b| b == b',') {
                out.insert(
                    (c[0].to_vec(), pos, c[3].to_vec(), alt.to_vec()),
                    genes.clone(),
                );
            }
        }
    }
    Ok(out)
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum How {
    Genotype,
    Ambiguous,
    Inconsistent,
}

impl How {
    fn as_bytes(self) -> &'static [u8] {
        match self {
            How::Genotype => b"genotype",
            How::Ambiguous => b"ambiguous",
            How::Inconsistent => b"inconsistent",
        }
    }
}

#[derive(Default)]
struct Counts {
    sites: u64,
    father: u64,
    mother: u64,
    both: u64,
    ambiguous: u64,
    inconsistent: u64,
    skipped_sex: u64,
}

/// A child's inheritance file, compressed on its own thread.
struct GzThread {
    tx: Option<SyncSender<Vec<u8>>>,
    handle: Option<JoinHandle<std::io::Result<AtomicFile>>>,
    buf: Vec<u8>,
    path: PathBuf,
}

impl GzThread {
    fn start(path: PathBuf) -> std::result::Result<GzThread, Die> {
        let file = AtomicFile::create(&path)
            .map_err(|e| Die::one(format!("cannot write {}: {e}", path.display())))?;
        let (tx, rx) = sync_channel::<Vec<u8>>(4);
        let handle = thread::spawn(move || {
            let mut gz = GzBuilder::new()
                .mtime(0)
                .operating_system(255)
                .write(file, Compression::new(GZ_LEVEL));
            for chunk in rx {
                gz.write_all(&chunk)?;
            }
            gz.finish()
        });
        Ok(GzThread {
            tx: Some(tx),
            handle: Some(handle),
            buf: Vec::with_capacity(CHUNK + 4096),
            path,
        })
    }

    fn flush(&mut self) -> std::result::Result<(), Die> {
        let chunk = std::mem::replace(&mut self.buf, Vec::with_capacity(CHUNK + 4096));
        let sent = self.tx.as_ref().is_some_and(|tx| tx.send(chunk).is_ok());
        if sent {
            Ok(())
        } else {
            Err(self
                .close()
                .err()
                .unwrap_or_else(|| self.fail("the writer stopped")))
        }
    }

    fn fail(&self, why: &str) -> Die {
        Die::one(format!("cannot write {}: {why}", self.path.display()))
    }

    /// Sends what is left, waits for the thread and returns the finished
    /// file, not yet renamed into place.
    fn close(&mut self) -> std::result::Result<AtomicFile, Die> {
        if !self.buf.is_empty()
            && let Some(tx) = &self.tx
        {
            let _ = tx.send(std::mem::take(&mut self.buf));
        }
        self.tx = None;
        match self.handle.take().map(JoinHandle::join) {
            Some(Ok(Ok(f))) => Ok(f),
            Some(Ok(Err(e))) => Err(self.fail(&e.to_string())),
            _ => Err(self.fail("the writer stopped")),
        }
    }
}

impl Drop for GzThread {
    /// On an early stop: end the thread, which drops its file unrenamed.
    fn drop(&mut self) {
        self.tx = None;
        if let Some(h) = self.handle.take() {
            let _ = h.join();
        }
    }
}

struct Child {
    name: String,
    sex: Option<String>,
    father: Option<String>,
    mother: Option<String>,
    f_seq: bool,
    m_seq: bool,
    col: usize,
    fcol: usize,
    mcol: usize,
    n: Counts,
    de_novo: Vec<Vec<String>>,
    cands: Vec<(Genes, String, Option<&'static str>)>,
    /// WhatsHap: per (chrom, PS) block, sites where the first haplotype
    /// holds the father's allele, and sites where it holds the mother's.
    blocks: HashMap<(Vec<u8>, Vec<u8>), [u64; 2]>,
    phased_hets: u64,
    out: GzThread,
}

/// The allele with index `i`: REF for 0, else the `i`-th ALT.
fn allele<'a>(
    cols: &[&'a [u8]],
    i: u32,
    path: &str,
    line: u64,
) -> std::result::Result<&'a [u8], Die> {
    let found = if i == 0 {
        Some(cols[3])
    } else {
        cols[4].split(|&c| c == b',').nth(i as usize - 1)
    };
    found.ok_or_else(|| {
        Die::one(format!(
            "{path} line {line}: a genotype names allele {i}, the record has fewer"
        ))
    })
}

/// `int(v)` of a FORMAT value, `None` when it is missing or not a number.
fn to_int(v: Option<&[u8]>) -> Option<i128> {
    let v = v.filter(|v| !v.is_empty() && *v != b".")?;
    pycompat::num::parse_int(std::str::from_utf8(v).ok()?).ok()
}

fn s_int(x: Option<i128>) -> String {
    x.map_or_else(|| ".".to_owned(), |v| v.to_string())
}

fn text(b: &[u8]) -> String {
    String::from_utf8_lossy(b).into_owned()
}

/// What the joint VCF pass needs.
struct Pass<'p> {
    path: &'p str,
    phased: bool,
    par: &'p Par,
    cands: &'p HashMap<String, HashMap<CandKey, Genes>>,
}

impl Pass<'_> {
    fn child(
        &self,
        ch: &mut Child,
        cols: &[&[u8]],
        pos: i64,
        gi: usize,
        format: &[u8],
        calls: &mut Calls<'_>,
    ) -> std::result::Result<(), Die> {
        let (path, line) = (self.path, calls.line);
        let chrom = cols[0];
        let (cg, cphased) = calls.gt_phased(ch.col)?;
        let Some(cg) = cg else { return Ok(()) };
        let fset = if ch.f_seq {
            calls.gt(ch.fcol)?.map(|g| set(&g))
        } else {
            None
        };
        let mset = if ch.m_seq {
            calls.gt(ch.mcol)?.map(|g| set(&g))
        } else {
            None
        };
        // a sequenced parent with no call here tells nothing
        let (f_known, m_known) = (fset.is_some(), mset.is_some());
        if !f_known && !m_known {
            return Ok(());
        }
        let is_x = chrom == b"chrX";
        let sexchr = (is_x || chrom == b"chrY") && !self.par.contains(&text(chrom), pos);
        let mut how = None;
        let (mut frm_f, mut frm_m) = (None::<u32>, None::<u32>);
        if sexchr {
            let sex = ch.sex.as_deref();
            if !matches!(sex, Some("male" | "female")) || (!is_x && sex == Some("female")) {
                ch.n.skipped_sex += 1;
                return Ok(());
            }
            if sex == Some("male") {
                if set(&cg).len() != 1 {
                    ch.n.skipped_sex += 1; // a heterozygous call on a single X or Y
                    return Ok(());
                }
                let a = cg[0];
                if (is_x && !m_known) || (!is_x && !f_known) {
                    return Ok(());
                }
                let src = if is_x { &mset } else { &fset };
                if src.as_ref().is_some_and(|s| s.contains(&a)) {
                    how = Some(How::Genotype);
                    if is_x {
                        frm_m = Some(a);
                    } else {
                        frm_f = Some(a);
                    }
                } else {
                    how = Some(How::Inconsistent);
                }
            }
        }
        let how = match how {
            Some(h) => h,
            None => {
                if cg.len() != 2 {
                    return Ok(());
                }
                let opts = resolve(&cg, fset.as_deref(), mset.as_deref());
                match opts.as_slice() {
                    [] => How::Inconsistent,
                    [(f, m)] => {
                        frm_f = Some(*f);
                        frm_m = Some(*m);
                        How::Genotype
                    }
                    _ => How::Ambiguous,
                }
            }
        };
        ch.n.sites += 1;
        match how {
            How::Genotype => {
                ch.n.both += u64::from(frm_f.is_some() && frm_m.is_some());
                ch.n.father += u64::from(frm_f.is_some());
                ch.n.mother += u64::from(frm_m.is_some());
            }
            How::Ambiguous => ch.n.ambiguous += 1,
            How::Inconsistent => ch.n.inconsistent += 1,
        }
        let cgt = calls.val(ch.col, gi)?.unwrap_or(b".");
        let fgt = if ch.f_seq {
            calls.val(ch.fcol, gi)?.unwrap_or(b".")
        } else {
            b"."
        };
        let mgt = if ch.m_seq {
            calls.val(ch.mcol, gi)?.unwrap_or(b".")
        } else {
            b"."
        };
        let from = |x: Option<u32>| -> std::result::Result<&[u8], Die> {
            match x {
                Some(i) => allele(cols, i, path, line),
                None => Ok(b"."),
            }
        };
        let (af, am) = (from(frm_f)?, from(frm_m)?);
        let buf = &mut ch.out.buf;
        for (k, f) in [
            chrom,
            cols[1],
            cols[3],
            cols[4],
            cgt,
            fgt,
            mgt,
            af,
            am,
            how.as_bytes(),
        ]
        .into_iter()
        .enumerate()
        {
            if k > 0 {
                buf.push(b'\t');
            }
            buf.extend_from_slice(f);
        }
        buf.push(b'\n');
        if buf.len() >= CHUNK {
            ch.out.flush()?;
        }
        let cset = set(&cg);
        let het = cset.len() == 2 && cg.contains(&0);
        // WhatsHap's phase of the child against the genotype origin
        if self.phased && cphased && het {
            ch.phased_hets += 1;
            let ps = vcf::format_index(format, b"PS")
                .and_then(|i| calls.val(ch.col, i).ok().flatten())
                .filter(|v| !v.is_empty() && *v != b".")
                .unwrap_or(b".");
            if how == How::Genotype
                && let Some(f) = frm_f
            {
                let b = ch.blocks.entry((chrom.to_vec(), ps.to_vec())).or_default();
                b[usize::from(cg[0] != f)] += 1;
            }
        }
        // de novo: an allele the child carries that neither parent has
        if how == How::Inconsistent && f_known && m_known {
            let (f, m) = (
                fset.as_deref().unwrap_or(&[]),
                mset.as_deref().unwrap_or(&[]),
            );
            // on a male's single X only the mother's alleles count, on his Y the father's;
            // a female's two X come one from each parent, so both parents' alleles count
            let seen: Vec<u32> = if sexchr && is_x {
                if ch.sex.as_deref() == Some("male") {
                    m.to_vec()
                } else {
                    [f, m].concat()
                }
            } else if sexchr {
                f.to_vec()
            } else {
                [f, m].concat()
            };
            let novel: Vec<u32> = cset
                .iter()
                .copied()
                .filter(|a| *a != 0 && !seen.contains(a))
                .collect();
            if !novel.is_empty() {
                let row = de_novo_row(ch, cols, format, gi, &novel, &cset, calls, path)?;
                ch.de_novo.push(row);
            }
        }
        // compound heterozygote candidates
        if het && let Some(table) = self.cands.get(&ch.name) {
            let ai = if cg[0] == 0 { cg[1] } else { cg[0] };
            let alt = allele(cols, ai, path, line)?;
            if let Some(genes) = table.get(&(chrom.to_vec(), pos, cols[3].to_vec(), alt.to_vec())) {
                let first = if cols[3] == alt {
                    0
                } else {
                    cols[4]
                        .split(|&c| c == b',')
                        .position(|x| x == alt)
                        .map_or(ai, |k| k as u32 + 1)
                };
                let origin = if frm_f == Some(first) {
                    Some("father")
                } else if frm_m == Some(first) {
                    Some("mother")
                } else {
                    None
                };
                ch.cands.push((
                    genes.clone(),
                    format!("{}:{pos}:{}>{}", text(chrom), text(cols[3]), text(alt)),
                    origin,
                ));
            }
        }
        Ok(())
    }
}

#[allow(clippy::too_many_arguments)]
fn de_novo_row(
    ch: &Child,
    cols: &[&[u8]],
    format: &[u8],
    gi: usize,
    novel: &[u32],
    cset: &[u32],
    calls: &mut Calls<'_>,
    path: &str,
) -> std::result::Result<Vec<String>, Die> {
    let line = calls.line;
    fn field<'a>(
        calls: &mut Calls<'a>,
        format: &[u8],
        col: usize,
        key: &[u8],
    ) -> std::result::Result<Option<&'a [u8]>, Die> {
        match vcf::format_index(format, key) {
            Some(i) => Ok(calls.val(col, i)?.filter(|v| !v.is_empty() && *v != b".")),
            None => Ok(None),
        }
    }
    let three = [ch.col, ch.fcol, ch.mcol];
    let mut dp = [None; 3];
    let mut gq = [None; 3];
    for (k, &c) in three.iter().enumerate() {
        dp[k] = to_int(field(calls, format, c, b"DP")?);
        gq[k] = to_int(field(calls, format, c, b"GQ")?);
    }
    let mut ab = None;
    if let Some(ad) = field(calls, format, ch.col, b"AD")? {
        let n: Vec<i128> = ad
            .split(|&c| c == b',')
            .map(|x| to_int(Some(x)).unwrap_or(0))
            .collect();
        let tot: i128 = n.iter().sum();
        if tot != 0 {
            ab = n.get(novel[0] as usize).map(|&v| v as f64 / tot as f64);
        }
    }
    let het = cset.len() == 2;
    let mut good = gq.iter().all(|x| x.is_some_and(|v| v >= MIN_GQ))
        && dp.iter().all(|x| x.is_some_and(|v| v >= MIN_DP));
    if het && ab.is_none_or(|ab| !(HET_AB.0 <= ab && ab <= HET_AB.1)) {
        good = false;
    }
    let mut alts = Vec::new();
    for &i in novel {
        alts.push(text(allele(cols, i, path, line)?));
    }
    let gt_of = |calls: &mut Calls<'_>, col: usize| -> std::result::Result<String, Die> {
        Ok(calls.val(col, gi)?.map_or_else(|| ".".to_owned(), text))
    };
    Ok(vec![
        ch.name.clone(),
        text(cols[0]),
        text(cols[1]),
        text(cols[3]),
        alts.join(","),
        gt_of(calls, ch.col)?,
        gt_of(calls, ch.fcol)?,
        gt_of(calls, ch.mcol)?,
        s_int(dp[0]),
        s_int(gq[0]),
        ab.map_or_else(|| ".".to_owned(), |v| fixed(v, 2)),
        s_int(dp[1]),
        s_int(gq[1]),
        s_int(dp[2]),
        s_int(gq[2]),
        if good { "pass" } else { "low" }.to_owned(),
    ])
}

fn pct(k: f64, n: f64) -> String {
    if n != 0.0 {
        format!("{}%", fixed(100.0 * k / n, 1))
    } else {
        "-".to_owned()
    }
}

pub fn run(a: &Args) -> Result<()> {
    let members_path = a.one("--members").unwrap_or_default();
    let members = read_members(members_path)?;
    let rel_path = a.one("--relationships").unwrap_or_default();
    let rels = read_table(rel_path, false)?;
    let par = Par::read(a.one("--par-bed"))?;
    let out_dir = PathBuf::from(a.one("--out-dir").unwrap_or_default());
    std::fs::create_dir_all(&out_dir)
        .map_err(|e| Die::one(format!("cannot create {}: {e}", out_dir.display())))?;
    let seq = segments::sequenced(&members);
    let is_seq = |p: &Option<String>| {
        p.as_deref()
            .and_then(|p| members.get(p))
            .is_some_and(|m| m.sequenced)
    };
    let mut cands: HashMap<String, HashMap<CandKey, Genes>> = HashMap::new();
    for kv in a.many("--candidates") {
        let (k, v) = kv.split_once('=').unwrap_or((kv, ""));
        if members.contains_key(k) {
            cands.insert(k.to_owned(), candidate_genes(v)?);
        }
    }
    let vcf_path = a.one("--vcf").unwrap_or_default();
    let phased_path = a.one("--phased");
    let path = phased_path.unwrap_or(vcf_path);

    // With --phased the map reads --vcf on its own thread.
    let seg_thread = phased_path.map(|_| {
        let (m, v, p) = (members.clone(), vcf_path.to_owned(), par.clone());
        thread::spawn(move || segments::collect(&m, &v, p))
    });

    let mut children = Vec::new();
    for s in &seq {
        let m = &members[s];
        if !(is_seq(&m.father) || is_seq(&m.mother)) {
            continue;
        }
        children.push(Child {
            name: s.clone(),
            sex: m.sex.clone(),
            father: m.father.clone(),
            mother: m.mother.clone(),
            f_seq: is_seq(&m.father),
            m_seq: is_seq(&m.mother),
            col: 0,
            fcol: 0,
            mcol: 0,
            n: Counts::default(),
            de_novo: Vec::new(),
            cands: Vec::new(),
            blocks: HashMap::new(),
            phased_hets: 0,
            out: GzThread::start(out_dir.join(format!("{s}_inheritance.tsv.gz")))?,
        });
    }
    for ch in &mut children {
        ch.out.buf.extend_from_slice(
            b"chrom\tpos\tref\talt\tchild_gt\tfather_gt\tmother_gt\tfrom_father\tfrom_mother\thow\n",
        );
    }

    let r = og_io::bgzf::open(path, og_io::bgzf::workers()).map_err(|e| Die::io(path, &e))?;
    let mut lines = og_vcf::Lines::new(r);
    let mut meta: Vec<Vec<u8>> = Vec::new();
    let names = segments::read_header(&mut lines, path, |l| {
        if l.starts_with(b"##contig") {
            meta.push(l.to_vec());
        }
    })?;
    let Some(names) = names else {
        return Err(Die::one(format!("{path} has no #CHROM line")).into());
    };
    let idx = vcf::sample_index(&names);
    let lost: Vec<&str> = seq
        .iter()
        .filter(|s| !idx.contains_key(*s))
        .map(String::as_str)
        .collect();
    if !lost.is_empty() {
        return Err(Die::one(format!(
            "{path} has no column for {} (its columns: {})",
            lost.join(", "),
            names.join(", ")
        ))
        .into());
    }
    for ch in &mut children {
        ch.col = idx[&ch.name];
        ch.fcol = ch
            .father
            .as_ref()
            .and_then(|f| idx.get(f))
            .copied()
            .unwrap_or(0);
        ch.mcol = ch
            .mother
            .as_ref()
            .and_then(|m| idx.get(m))
            .copied()
            .unwrap_or(0);
    }
    // Without --phased the map is built in this same pass.
    let mut inline: Option<Collector> = phased_path.is_none().then(|| {
        let mut c = Collector::new(&members, &idx, par.clone());
        for l in &meta {
            c.header(l);
        }
        c
    });
    let pass = Pass {
        path,
        phased: phased_path.is_some(),
        par: &par,
        cands: &cands,
    };
    let mut n = lines.number();
    while let Some(line) = lines.next_line().map_err(|e| Die::io(path, &e))? {
        n += 1;
        let cols = vcf::columns(line);
        if cols.len() < 2 {
            return Err(Die::one(format!("{path} line {n}: not a VCF record")).into());
        }
        let pos = vcf::pos(cols[1], path, n)?;
        let primary = is_primary(cols[0]);
        if primary && cols.len() < 9 {
            return Err(Die::one(format!("{path} line {n}: no FORMAT column")).into());
        }
        let format = cols.get(8).copied().unwrap_or(b"");
        let gi = if primary {
            vcf::format_index(format, b"GT")
        } else {
            None
        };
        let mut calls = Calls::new(&cols, gi.unwrap_or(0), path, n);
        if let Some(c) = &mut inline {
            c.record(path, pos, &mut calls, &cols, gi)?;
        }
        let Some(gi) = gi else { continue };
        for ch in &mut children {
            pass.child(ch, &cols, pos, gi, format, &mut calls)?;
        }
    }

    let map: Map = match (inline, seg_thread) {
        (Some(c), _) => c.finish(&members),
        (None, Some(h)) => {
            let (c, _) = h
                .join()
                .map_err(|_| Die::one("the chromosome map stopped"))??;
            c.finish(&members)
        }
        (None, None) => unreachable!("one of the two builds the map"),
    };

    // Every output is complete: write them.
    let mut dn = String::from(
        "child\tchrom\tpos\tref\talt\tchild_gt\tfather_gt\tmother_gt\tchild_dp\tchild_gq\tchild_ab\tfather_dp\tfather_gq\tmother_dp\tmother_gq\tquality\n",
    );
    for ch in &children {
        for row in &ch.de_novo {
            dn.push_str(&row.join("\t"));
            dn.push('\n');
        }
    }
    write_file(&out_dir.join("de_novo.tsv"), dn.as_bytes())?;
    let mut comp: Vec<[String; 7]> = Vec::new();
    for ch in &children {
        let mut by_gene: IndexMap<&str, Vec<(&str, Option<&str>)>> = IndexMap::new();
        for (genes, var, origin) in &ch.cands {
            for g in genes.iter() {
                by_gene.entry(g).or_default().push((var, *origin));
            }
        }
        by_gene.sort_keys();
        for (g, list) in &by_gene {
            for (i, (v1, o1)) in list.iter().enumerate() {
                for (v2, o2) in &list[i + 1..] {
                    let phase = match (o1, o2) {
                        (Some(a), Some(b)) if a == b => "in cis",
                        (Some(_), Some(_)) => "in trans",
                        _ => "unknown",
                    };
                    comp.push([
                        ch.name.clone(),
                        (*g).to_owned(),
                        (*v1).to_owned(),
                        (*v2).to_owned(),
                        o1.unwrap_or(".").to_owned(),
                        o2.unwrap_or(".").to_owned(),
                        phase.to_owned(),
                    ]);
                }
            }
        }
    }
    let mut ct = String::from("child\tgene\tvariant_a\tvariant_b\ta_from\tb_from\tphase\n");
    for row in &comp {
        ct.push_str(&row.join("\t"));
        ct.push('\n');
    }
    write_file(&out_dir.join("compound_het.tsv"), ct.as_bytes())?;
    for ch in &mut children {
        let f = ch.out.close()?;
        f.commit().map_err(|e| ch.out.fail(&e.to_string()))?;
    }
    map.write(&out_dir)?;
    let seg = map.markdown(&members);

    let md = markdown(
        &members,
        &rels,
        rel_path,
        &children,
        &comp,
        !cands.is_empty(),
        phased_path.is_some(),
        seg,
    )?;
    write_file(&out_dir.join("family_report.md"), md.as_bytes())?;
    print_out(&md);
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn markdown(
    members: &Members,
    rels: &[IndexMap<String, String>],
    rel_path: &str,
    children: &[Child],
    comp: &[[String; 7]],
    have_cands: bool,
    phased: bool,
    seg: Vec<String>,
) -> std::result::Result<String, Die> {
    let mut l: Vec<String> = vec!["# Family analysis".into(), String::new()];
    let mut seq: Vec<&str> = members
        .iter()
        .filter(|(_, m)| m.sequenced)
        .map(|(s, _)| s.as_str())
        .collect();
    seq.sort_unstable();
    let mut others: Vec<&str> = members
        .iter()
        .filter(|(_, m)| !m.sequenced)
        .map(|(s, _)| s.as_str())
        .collect();
    others.sort_unstable();
    l.push(format!(
        "Sequenced: {}. In the pedigree, not sequenced: {}.",
        seq.join(", "),
        if others.is_empty() {
            "none".to_owned()
        } else {
            others.join(", ")
        }
    ));
    l.extend([
        String::new(),
        "## Relationships".into(),
        String::new(),
        "Relatedness is twice the kinship coefficient: about the share of the genome the two have \
         inherited from common ancestors (1 the same person, 0.5 parent and child or full siblings, \
         0.25 second degree, 0 unrelated). somalier measures it at its sites; the pedigree gives the \
         expected value."
            .into(),
        String::new(),
        "| Pair | Declared | Expected | Measured | Opposite homozygotes | Sites | Verdict |".into(),
        "|---|---|---|---|---|---|---|".into(),
    ]);
    for r in rels {
        let c = |k: &str| cell(r, k, rel_path);
        let rate = c("ibs0_rate")?;
        let rate = pycompat::num::parse_float(rate)
            .map_err(|_| Die::one(format!("{rel_path}: ibs0_rate '{rate}' is not a number")))?;
        l.push(format!(
            "| {} and {} | {} | {} | {} | {} | {} | {}: {} |",
            c("sample_a")?,
            c("sample_b")?,
            c("declared")?,
            c("expected_relatedness")?,
            c("relatedness")?,
            pct(rate, 1.0),
            c("sites")?,
            c("verdict")?,
            c("reason")?
        ));
    }
    for ch in children {
        let n = &ch.n;
        let sites = n.sites as f64;
        l.extend([String::new(), format!("## {}", ch.name), String::new()]);
        let who = |p: &Option<String>, s: bool| match p {
            Some(p) => format!("{p} ({})", if s { "sequenced" } else { "not sequenced" }),
            None => "not named".to_owned(),
        };
        l.push(format!(
            "Father: {}. Mother: {}.",
            who(&ch.father, ch.f_seq),
            who(&ch.mother, ch.m_seq)
        ));
        let file = Path::new(&ch.out.path)
            .file_name()
            .map(|f| f.to_string_lossy().into_owned())
            .unwrap_or_default();
        l.extend([
            String::new(),
            format!(
                "Sites of the joint VCF where {} has a genotype and a sequenced parent does too: {}. `{file}` lists them.",
                ch.name, n.sites
            ),
            String::new(),
        ]);
        l.push(format!(
            "- the allele from the father is known at {} ({})",
            n.father,
            pct(n.father as f64, sites)
        ));
        l.push(format!(
            "- the allele from the mother is known at {} ({})",
            n.mother,
            pct(n.mother as f64, sites)
        ));
        l.push(format!(
            "- both parents heterozygous like the child, so no origin from the genotypes: {} ({})",
            n.ambiguous,
            pct(n.ambiguous as f64, sites)
        ));
        l.push(format!(
            "- inconsistent with the declared parents (a genotype error, a deletion, or a new mutation): {} ({})",
            n.inconsistent,
            pct(n.inconsistent as f64, sites)
        ));
        if n.skipped_sex > 0 {
            l.push(format!(
                "- left out on chrX or chrY outside the pseudoautosomal regions (sex unknown, or a \
                 heterozygous call on a single X or Y): {}",
                n.skipped_sex
            ));
        }
        for (who, s, k) in [
            ("father", ch.f_seq, n.father),
            ("mother", ch.m_seq, n.mother),
        ] {
            if !s {
                l.extend([
                    String::new(),
                    format!(
                        "The {who} is not sequenced. At {k} of these sites ({}) the genotypes give the allele \
                         the {who} passed on to {}: that is the half of the {who}'s genome {} carries. The other \
                         half was not passed on and is not seen here.",
                        pct(k as f64, sites),
                        ch.name,
                        ch.name
                    ),
                ]);
            }
        }
        if ch.f_seq && ch.m_seq {
            let good = ch
                .de_novo
                .iter()
                .filter(|r| r.last().is_some_and(|q| q == "pass"))
                .count();
            l.extend([
                String::new(),
                format!(
                    "De novo candidates (an allele neither parent has): {}, {good} of them with GQ >= {MIN_GQ} and \
                     depth >= {MIN_DP} in all three and, for a heterozygous call, an allele balance of {} to {} \
                     (`de_novo.tsv`). Most candidates of a genome are genotype errors; a candidate that \
                     matters needs confirming by another method.",
                    ch.de_novo.len(),
                    fixed(HET_AB.0, 2),
                    fixed(HET_AB.1, 2)
                ),
            ]);
        } else {
            l.extend([
                String::new(),
                "No de novo candidates are called: that needs both parents sequenced.".into(),
            ]);
        }
        if phased {
            let votes: u64 = ch.blocks.values().map(|b| b[0].max(b[1])).sum();
            let tot: u64 = ch.blocks.values().map(|b| b[0] + b[1]).sum();
            l.extend([
                String::new(),
                format!(
                    "WhatsHap phased {} heterozygous sites of {} by inheritance. At {tot} of them the \
                     genotypes also give the origin, and the phase agrees with it at {votes} ({}).",
                    ch.phased_hets,
                    ch.name,
                    pct(votes as f64, tot as f64)
                ),
            ]);
        }
    }
    l.extend([
        String::new(),
        "## Compound heterozygotes".into(),
        String::new(),
    ]);
    if !have_cands {
        l.push(
            "No candidate variants were given (the clinical filter, step 23, did not run for a child), \
             so no pair is reported."
                .into(),
        );
    } else {
        let t = comp.iter().filter(|r| r[6] == "in trans").count();
        let c = comp.iter().filter(|r| r[6] == "in cis").count();
        l.push(format!(
            "Pairs of the clinical filter's candidate variants a child carries heterozygous in one gene: \
             {}; in trans (one from each parent, both copies of the gene hit) {t}, in cis (both from one \
             parent) {c}, origin unknown {} (`compound_het.tsv`).",
            comp.len(),
            comp.len() - t - c
        ));
    }
    l.extend(seg);
    let mut s = l.join("\n");
    s.push('\n');
    Ok(s)
}
