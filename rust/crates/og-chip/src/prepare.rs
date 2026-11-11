//! `og chip prepare`: for each gene of the table, its MANE Select transcript
//! in Ensembl's GFF3 (the transcript tagged `MANE_Select`), written with the
//! gene line and the transcript's exons, CDS and UTRs under UCSC contig
//! names (`1` -> `chr1`, `MT` -> `chrM`): the file `bcftools csq` reads.
//! Its coding exons, 10 bases wider on each side so the splice sites are
//! inside, go to a BED of the contigs the reference has. A table gene with
//! no MANE Select transcript stops the run.
//!
//! One pass over the GFF3 with one line in memory. A line is split into its
//! nine columns, and column 9 is searched only for the one or two
//! attributes the line's place in the chain needs; only lines of the table's
//! genes are kept. The Python built a dictionary of column 9 for each of the
//! ~3 million lines.
//!
//! The output is the Python's, byte for byte, line order included.

use crate::table::read_table;
use crate::{Error, each_line, int, write_atomic};
use indexmap::IndexMap;
use memchr::memmem::FinderRev;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, Read};
use std::sync::mpsc::{Receiver, sync_channel};
use std::thread;

/// Bases each side of a coding exon: the splice donor and acceptor sites.
pub const PAD: i128 = 10;

/// The value of attribute `key` in GFF3 column 9 (the last one, when it is
/// there twice, as a dictionary built in order keeps it).
pub fn attr<'a>(col9: &'a str, key: &str) -> Option<&'a str> {
    col9.split(';')
        .filter_map(|kv| kv.split_once('='))
        .filter(|(k, _)| *k == key)
        .map(|(_, v)| v)
        .next_back()
}

/// The UCSC name of an Ensembl contig.
pub fn ucsc(seqid: &str) -> String {
    match seqid {
        "MT" | "M" | "chrM" => "chrM".to_owned(),
        s if s.starts_with("chr") => s.to_owned(),
        s => format!("chr{s}"),
    }
}

/// The nine columns of a GFF3 line (without its `\n`), or `None` when it has
/// another count.
pub fn columns(line: &str) -> Option<[&str; 9]> {
    let mut f = [""; 9];
    let mut n = 0;
    for part in line.split('\t') {
        if n == 9 {
            return None;
        }
        f[n] = part;
        n += 1;
    }
    (n == 9).then_some(f)
}

type Line = Vec<String>;

fn owned(f: &[&str; 9]) -> Line {
    f.iter().map(|s| (*s).to_owned()).collect()
}

/// What the pass over the GFF3 keeps.
#[derive(Default)]
struct Kept {
    /// Gene symbol -> its gene line.
    gene_line: IndexMap<String, Line>,
    /// Gene id -> symbol.
    gene_of: HashMap<String, String>,
    /// Gene symbol -> its MANE transcript line.
    tx_line: IndexMap<String, Line>,
    /// Transcript id -> its child lines, in file order.
    children: HashMap<String, Vec<Line>>,
}

/// The value of an attribute in column 9, as bytes: the last segment
/// between `;` that starts with `key=` (the finder's needle), as [`attr`]
/// reads it.
fn attr_bytes<'a>(col9: &'a [u8], key_eq: &FinderRev<'_>) -> Option<&'a [u8]> {
    let n = key_eq.needle().len();
    let mut end = col9.len();
    while let Some(at) = key_eq.rfind(&col9[..end]) {
        if at == 0 || col9[at - 1] == b';' {
            let v = &col9[at + n..];
            return Some(&v[..memchr::memchr(b';', v).unwrap_or(v.len())]);
        }
        end = at + n - 1;
    }
    None
}

/// The two attribute searches of the scan, built once.
struct Finders {
    name: FinderRev<'static>,
    parent: FinderRev<'static>,
}

impl Finders {
    fn new() -> Finders {
        Finders {
            name: FinderRev::new(b"Name="),
            parent: FinderRev::new(b"Parent="),
        }
    }
}

impl Kept {
    /// One line without its ending. Most lines are neither a table gene nor
    /// a child of one: they are told apart on bytes, and only a line that
    /// may be kept is decoded and split.
    ///
    /// The type is read between the second and third tabs and column 9 after
    /// the last one. That is right for every line of nine columns; a line
    /// of another count that looks like a candidate is dropped by
    /// [`columns`] after decoding, as the Python dropped it.
    fn scan(
        &mut self,
        line: &[u8],
        want: &HashSet<String>,
        find: &Finders,
        gff: &str,
    ) -> Result<(), Error> {
        if line.first() == Some(&b'#') {
            return Ok(());
        }
        let mut tabs = memchr::memchr_iter(b'\t', line);
        let (Some(t2), Some(t3)) = (tabs.nth(1), tabs.next()) else {
            return Ok(());
        };
        let Some(last) = memchr::memrchr(b'\t', line) else {
            return Ok(());
        };
        let col9 = &line[last + 1..];
        let str_in = |v: &[u8], keys: &dyn Fn(&str) -> bool| std::str::from_utf8(v).is_ok_and(keys);
        let candidate = if &line[t2 + 1..t3] == b"gene"
            && attr_bytes(col9, &find.name).is_some_and(|n| str_in(n, &|n| want.contains(n)))
        {
            true
        } else {
            attr_bytes(col9, &find.parent).is_some_and(|p| {
                str_in(p, &|p| {
                    self.gene_of.contains_key(p) || self.children.contains_key(p)
                })
            })
        };
        if candidate {
            let text = std::str::from_utf8(line).map_err(|e| Error(format!("{gff}: {e}")))?;
            if let Some(f) = columns(text) {
                self.line(&f, want);
            }
        }
        Ok(())
    }

    fn line(&mut self, f: &[&str; 9], want: &HashSet<String>) {
        let c9 = f[8];
        if f[2] == "gene"
            && let Some(name) = attr(c9, "Name")
            && want.contains(name)
        {
            self.gene_line.insert(name.to_owned(), owned(f));
            self.gene_of
                .insert(attr(c9, "ID").unwrap_or("").to_owned(), name.to_owned());
            return;
        }
        let Some(parent) = attr(c9, "Parent") else {
            return;
        };
        if let Some(g) = self.gene_of.get(parent)
            && attr(c9, "tag").is_some_and(|t| t.split(',').any(|x| x == "MANE_Select"))
        {
            self.tx_line.insert(g.clone(), owned(f));
            self.children
                .entry(attr(c9, "ID").unwrap_or("").to_owned())
                .or_default();
        } else if let Some(kids) = self.children.get_mut(parent) {
            kids.push(owned(f));
        }
    }
}

/// The input's bytes in chunks, decompressed on a thread of its own so the
/// scan does not wait for it (an Ensembl GFF3 is one gzip stream, about
/// 1.7 GB inflated).
struct Inflated {
    rx: Receiver<io::Result<Vec<u8>>>,
}

/// Bytes per chunk, and chunks in flight: about 1 MB of memory.
const CHUNK: usize = 1 << 18;
const IN_FLIGHT: usize = 3;

impl Inflated {
    fn open(path: &str) -> io::Result<Inflated> {
        let file = File::open(path)?;
        let (tx, rx) = sync_channel(IN_FLIGHT);
        thread::spawn(move || {
            let mut r = match og_io::reader(file) {
                Ok(r) => r,
                Err(e) => {
                    let _ = tx.send(Err(e));
                    return;
                }
            };
            loop {
                let mut buf = vec![0u8; CHUNK];
                let sent = match r.read(&mut buf) {
                    Ok(0) => return,
                    Ok(n) => {
                        buf.truncate(n);
                        tx.send(Ok(buf))
                    }
                    Err(e) if e.kind() == io::ErrorKind::Interrupted => continue,
                    Err(e) => tx.send(Err(e)),
                };
                if sent.is_err() {
                    return;
                }
            }
        });
        Ok(Inflated { rx })
    }
}

impl Iterator for Inflated {
    type Item = io::Result<Vec<u8>>;
    fn next(&mut self) -> Option<Self::Item> {
        self.rx.recv().ok()
    }
}

/// One line without its `\n`, split further where Python's universal
/// newlines would (`\r\n` and a lone `\r` end a line too).
fn emit(line: &[u8], each: &mut impl FnMut(&[u8]) -> Result<(), Error>) -> Result<(), Error> {
    let body = line.strip_suffix(b"\r").unwrap_or(line);
    if memchr::memchr(b'\r', body).is_none() {
        each(body)
    } else {
        body.split(|&c| c == b'\r').try_for_each(each)
    }
}

/// Each line of the chunks without its ending, in order. A line inside a
/// chunk is handed over in place; only a line across two chunks is copied.
fn each_raw_line(
    chunks: impl Iterator<Item = io::Result<Vec<u8>>>,
    mut each: impl FnMut(&[u8]) -> Result<(), Error>,
    path: &str,
) -> Result<(), Error> {
    let mut carry = Vec::new();
    for chunk in chunks {
        let chunk = chunk.map_err(|e| Error::io(&e, path))?;
        let mut start = 0;
        for nl in memchr::memchr_iter(b'\n', &chunk) {
            if carry.is_empty() {
                emit(&chunk[start..nl], &mut each)?;
            } else {
                carry.extend_from_slice(&chunk[start..nl]);
                emit(&carry, &mut each)?;
                carry.clear();
            }
            start = nl + 1;
        }
        carry.extend_from_slice(&chunk[start..]);
    }
    if carry.is_empty() {
        Ok(())
    } else {
        emit(&carry, &mut each)
    }
}

/// The contig names of a `.fai`: each line up to its first tab.
fn contigs(fai: &str) -> Result<HashSet<String>, Error> {
    let mut out = HashSet::new();
    each_line(fai, |ln| {
        out.insert(ln.split('\t').next().unwrap_or("").to_owned());
        Ok(())
    })?;
    Ok(out)
}

/// The GFF3 and BED texts for the table's genes.
pub fn build(gff: &str, genes: &str, fai: &str) -> Result<(String, String), Error> {
    let want: HashSet<String> = read_table(genes)?.into_iter().map(|r| r.gene).collect();
    let mut kept = Kept::default();
    let input = Inflated::open(gff).map_err(|e| Error::io(&e, gff))?;
    let find = Finders::new();
    each_raw_line(input, |ln| kept.scan(ln, &want, &find, gff), gff)?;
    let mut missing: Vec<&String> = want
        .iter()
        .filter(|g| !kept.tx_line.contains_key(*g))
        .collect();
    if missing.len() > want.len() {
        missing.sort();
        let names: Vec<&str> = missing.iter().map(|s| s.as_str()).collect();
        return Err(Error(format!(
            "{gff}: no MANE Select transcript for {} (an Ensembl GFF3 of release 108 or later tags them MANE_Select)",
            names.join(", ")
        )));
    }
    let contigs = contigs(fai)?;

    let mut order = Vec::with_capacity(kept.tx_line.len());
    for g in kept.tx_line.keys() {
        let line = &kept.gene_line[g];
        order.push((ucsc(&line[0]), int(&line[3], gff)?, g));
    }
    order.sort_by(|a, b| (&a.0, a.1).cmp(&(&b.0, b.1)));

    let mut out = String::from("##gff-version 3\n");
    let mut bed: Vec<(String, i128, i128, &str)> = Vec::new();
    for (_, _, g) in &order {
        let tx = &kept.tx_line[*g];
        let id = attr(&tx[8], "ID").ok_or_else(|| {
            Error(format!(
                "{gff}: the MANE Select transcript of {g} has no ID"
            ))
        })?;
        let mut kids = Vec::new();
        for f in &kept.children[id] {
            kids.push((int(&f[3], gff)?, f));
        }
        kids.sort_by_key(|k| k.0);
        for f in [&kept.gene_line[*g], tx]
            .into_iter()
            .chain(kids.iter().map(|k| k.1))
        {
            out.push_str(&ucsc(&f[0]));
            for c in &f[1..] {
                out.push('\t');
                out.push_str(c);
            }
            out.push('\n');
        }
        out.push_str("###\n");
        let chrom = ucsc(&tx[0]);
        if contigs.contains(&chrom) {
            for (start, f) in &kids {
                if f[2] == "CDS" {
                    let end = int(&f[4], gff)?;
                    bed.push((
                        chrom.clone(),
                        (start - 1 - PAD).max(0),
                        end + PAD,
                        g.as_str(),
                    ));
                }
            }
        }
    }
    bed.sort_by(|a, b| (&a.0, a.1).cmp(&(&b.0, b.1)));
    let bed: String = bed
        .iter()
        .map(|(c, s, e, g)| format!("{c}\t{s}\t{e}\t{g}\n"))
        .collect();
    Ok((out, bed))
}

/// `og chip prepare`: writes both files, each whole or not at all.
pub fn prepare(
    gff: &str,
    genes: &str,
    fai: &str,
    out_gff: &str,
    out_bed: &str,
) -> Result<(), Error> {
    let (g, b) = build(gff, genes, fai)?;
    write_atomic(out_gff, g.as_bytes())?;
    write_atomic(out_bed, b.as_bytes())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn attributes_as_a_dictionary_reads_them() {
        let c9 = "ID=transcript:T2;Parent=gene:G1;tag=basic,MANE_Select;x;Parent=gene:G9";
        assert_eq!(attr(c9, "ID"), Some("transcript:T2"));
        assert_eq!(attr(c9, "Parent"), Some("gene:G9"));
        assert_eq!(attr(c9, "x"), None);
        assert_eq!(attr("a=b=c", "a"), Some("b=c"));
    }

    #[test]
    fn contig_names() {
        assert_eq!(ucsc("2"), "chr2");
        assert_eq!(ucsc("MT"), "chrM");
        assert_eq!(ucsc("chrX"), "chrX");
        assert_eq!(ucsc("KI270728.1"), "chrKI270728.1");
    }

    /// A scratch directory of its own, removed when dropped.
    struct Scratch(std::path::PathBuf);
    impl Scratch {
        fn new(tag: &str) -> Scratch {
            let d = std::env::temp_dir().join(format!("og-chip-{tag}-{}", std::process::id()));
            std::fs::create_dir_all(&d).unwrap();
            Scratch(d)
        }
        fn file(&self, name: &str, body: &str) -> String {
            let p = self.0.join(name);
            std::fs::write(&p, body).unwrap();
            p.to_string_lossy().into_owned()
        }
    }
    impl Drop for Scratch {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }

    const GFF: &str = "##gff-version 3
2\te\tgene\t1000\t5000\t.\t-\t.\tID=gene:G1;Name=DNMT3A
2\te\tmRNA\t1000\t5000\t.\t-\t.\tID=transcript:T1;Parent=gene:G1;tag=basic;transcript_id=T1
2\te\tCDS\t1100\t1200\t.\t-\t0\tParent=transcript:T1
2\te\tmRNA\t1000\t5000\t.\t-\t.\tID=transcript:T2;Parent=gene:G1;tag=basic,MANE_Select;transcript_id=T2
2\te\texon\t4000\t5000\t.\t-\t.\tParent=transcript:T2;rank=1
2\te\tCDS\t4000\t4100\t.\t-\t0\tParent=transcript:T2
2\te\texon\t1000\t1300\t.\t-\t.\tParent=transcript:T2;rank=2
2\te\tCDS\t1050\t1300\t.\t-\t0\tParent=transcript:T2
";

    #[test]
    fn the_mane_transcript_and_its_padded_coding_exons() {
        let d = Scratch::new("mane");
        let genes = d.file(
            "g.tsv",
            "gene\ttruncating\tmissense\tnote\tpmids\nDNMT3A\tall\tR882\tn\t25426837\n",
        );
        let fai = d.file("r.fai", "chr2\t10000\t6\t60\t61\n");
        let (gff, bed) = build(&d.file("in.gff3", GFF), &genes, &fai).unwrap();
        assert!(gff.starts_with("##gff-version 3\nchr2\te\tgene\t1000"));
        assert!(gff.contains("transcript:T2") && !gff.contains("T1"));
        assert!(
            gff.find("\t1050\t").unwrap() < gff.find("\t4000\t").unwrap(),
            "children by start"
        );
        assert_eq!(bed, "chr2\t1039\t1310\tDNMT3A\nchr2\t3989\t4110\tDNMT3A\n");
        let off = d.file("o.fai", "chr9\t10000\t6\t60\t61\n");
        assert_eq!(build(&d.file("in.gff3", GFF), &genes, &off).unwrap().1, "");
    }

    #[test]
    fn red_first_a_gene_without_mane_stops_prepare() {
        let d = Scratch::new("nomane");
        let genes = d.file(
            "g.tsv",
            "gene\ttruncating\tmissense\tnote\tpmids\nDNMT3A\tall\tR882\tn\t25426837\nJAK2\tnone\tV617F\tn\t15858187\n",
        );
        let fai = d.file("r.fai", "chr2\t10000\t6\t60\t61\n");
        let gff = d.file("in.gff3", &GFF.replace("basic,MANE_Select", "basic"));
        let e = build(&gff, &genes, &fai).unwrap_err();
        assert!(e.0.ends_with("no MANE Select transcript for DNMT3A, JAK2 (an Ensembl GFF3 of release 108 or later tags them MANE_Select)"), "{e}");
    }

    #[test]
    fn lines_split_as_python_reads_them() {
        let mut got = Vec::new();
        let chunks = vec![
            Ok(b"a\r\nb\rc\n\nd".to_vec()),
            Ok(b"e\r".to_vec()),
            Ok(b"\nf\r".to_vec()),
        ];
        each_raw_line(
            chunks.into_iter(),
            |l| {
                got.push(String::from_utf8(l.to_vec()).unwrap());
                Ok(())
            },
            "x",
        )
        .unwrap();
        assert_eq!(got, ["a", "b", "c", "", "de", "f"]);
    }

    #[test]
    fn nine_columns_only() {
        assert!(columns("1\t2\t3\t4\t5\t6\t7\t8\t9").is_some());
        assert!(columns("1\t2\t3\t4\t5\t6\t7\t8").is_none());
        assert!(columns("1\t2\t3\t4\t5\t6\t7\t8\t9\t10").is_none());
    }
}
