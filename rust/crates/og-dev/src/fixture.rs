//! `og dev fixture`: writers of committed or generated test fixtures.
//!
//! - `og dev fixture comphet OUT_DIR [--swap]`: the compound-het phasing
//!   fixture of `tests/e2e/comphet-phasing-1-reads.sh` (made-up reads on a
//!   made-up contig). It replaces `tests/fixtures/comphet/make_fixture.py`
//!   and writes the same bytes: the contig is Python's `random.Random(31)`
//!   choosing among `ACGT`, so [`Mt19937`] is CPython's generator.
//! - `og dev fixture trim-pharmcat REPORT OUT [MAX_ANNOTATIONS]`: a PharmCAT
//!   `report.json` cut down to the fields `og pgx` reads. It replaces
//!   `tests/fixtures/pharmcat/trim_report.py` and writes the same bytes.
//!
//! Exit 0 when written, 1 when an input cannot be read or has the wrong
//! shape, 2 on a usage error.

use crate::py_os_error;
use indexmap::IndexMap;
use pycompat::json::{self, Dump, Value};
use std::fmt::Write as _;
use std::fs;
use std::io::Write;
use std::path::Path;

const USAGE: &str = "usage: og dev fixture comphet OUT_DIR [--swap]\n       og dev fixture trim-pharmcat REPORT.json OUT.json [MAX_ANNOTATIONS]\n";

/// Runs `og dev fixture ARGS`; returns the exit code.
pub fn cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let usage = |err: &mut dyn Write, msg: &str| {
        let _ = writeln!(err, "{USAGE}og dev fixture: error: {msg}");
        2
    };
    let result = match args.first().map(String::as_str) {
        Some("comphet") => {
            let rest: Vec<&String> = args[1..].iter().filter(|a| *a != "--swap").collect();
            if rest.len() != 1 {
                return usage(err, "comphet takes one OUT_DIR and an optional --swap");
            }
            comphet::write(Path::new(rest[0]), args[1..].iter().any(|a| a == "--swap"))
        }
        Some("trim-pharmcat") => {
            let rest = &args[1..];
            if !(2..=3).contains(&rest.len()) {
                return usage(
                    err,
                    "trim-pharmcat takes REPORT.json OUT.json [MAX_ANNOTATIONS]",
                );
            }
            let max = match rest.get(2).map(|m| m.trim().parse::<i64>()) {
                None => 4,
                Some(Ok(n)) => n,
                Some(Err(_)) => {
                    return usage(
                        err,
                        &format!("MAX_ANNOTATIONS '{}' is not a whole number", rest[2]),
                    );
                }
            };
            trim_pharmcat(Path::new(&rest[0]), Path::new(&rest[1]), max)
        }
        Some("-h" | "--help") => {
            let _ = out.write_all(USAGE.as_bytes());
            return 0;
        }
        Some(other) => return usage(err, &format!("unknown fixture '{other}'")),
        None => return usage(err, "the following arguments are required: fixture"),
    };
    match result {
        Ok(()) => 0,
        Err(e) => {
            let _ = writeln!(err, "og dev fixture: {e}");
            1
        }
    }
}

/// CPython's Mersenne Twister (`random.Random`), only as far as
/// `Random(int).choice(seq)` needs it.
pub struct Mt19937 {
    state: [u32; 624],
    index: usize,
}

impl Mt19937 {
    /// `random.Random(seed)` for a seed below 2^32: `init_by_array([seed])`.
    pub fn new(seed: u32) -> Self {
        let mut mt = Mt19937 {
            state: [0; 624],
            index: 624,
        };
        mt.init_genrand(19_650_218);
        let key = [seed];
        let (mut i, mut j) = (1usize, 0usize);
        for _ in 0..624.max(key.len()) {
            let prev = mt.state[i - 1] ^ (mt.state[i - 1] >> 30);
            mt.state[i] = (mt.state[i] ^ prev.wrapping_mul(1_664_525))
                .wrapping_add(key[j])
                .wrapping_add(j as u32);
            i += 1;
            j += 1;
            if i >= 624 {
                mt.state[0] = mt.state[623];
                i = 1;
            }
            if j >= key.len() {
                j = 0;
            }
        }
        for _ in 0..623 {
            let prev = mt.state[i - 1] ^ (mt.state[i - 1] >> 30);
            mt.state[i] = (mt.state[i] ^ prev.wrapping_mul(1_566_083_941)).wrapping_sub(i as u32);
            i += 1;
            if i >= 624 {
                mt.state[0] = mt.state[623];
                i = 1;
            }
        }
        mt.state[0] = 0x8000_0000;
        mt
    }

    fn init_genrand(&mut self, s: u32) {
        self.state[0] = s;
        for i in 1..624 {
            let prev = self.state[i - 1] ^ (self.state[i - 1] >> 30);
            self.state[i] = 1_812_433_253u32.wrapping_mul(prev).wrapping_add(i as u32);
        }
        self.index = 624;
    }

    /// `genrand_uint32`.
    pub fn next_u32(&mut self) -> u32 {
        if self.index >= 624 {
            for k in 0..624 {
                let y = (self.state[k] & 0x8000_0000) | (self.state[(k + 1) % 624] & 0x7fff_ffff);
                let mag = if y & 1 == 1 { 0x9908_b0df } else { 0 };
                self.state[k] = self.state[(k + 397) % 624] ^ (y >> 1) ^ mag;
            }
            self.index = 0;
        }
        let mut y = self.state[self.index];
        self.index += 1;
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c_5680;
        y ^= (y << 15) & 0xefc6_0000;
        y ^ (y >> 18)
    }

    /// `_randbelow(n)` for `0 < n < 2^32`: rejection sampling on
    /// `getrandbits(n.bit_length())`.
    pub fn below(&mut self, n: u32) -> u32 {
        let k = 32 - n.leading_zeros();
        loop {
            let r = self.next_u32() >> (32 - k);
            if r < n {
                return r;
            }
        }
    }
}

pub mod comphet {
    //! The compound-het phasing fixture: sample CHPH on contig chrP (80 kb
    //! of seeded random sequence), read pairs of 100 bp from 500 bp
    //! fragments, and three pairs of heterozygous variants:
    //!
    //! - GENEC: SNV at 1000 and a 2-base deletion at 1003, both on copy 1
    //!   (cis); single reads cover both. `--swap` moves the deletion to
    //!   copy 2 (trans).
    //! - GENET: SNVs at 5000 and 5400, one on each copy (trans); only a read
    //!   and its mate cover both.
    //! - GENEF: SNVs at 20000 and 70000, one on each copy, 50 kb apart: no
    //!   read pair covers both (unknown).
    //!
    //! Files: `reads.sam`, `vep.vcf` (with a VEP CSQ, for step 31) and
    //! `comphet.vcf` (the pairs as slivar names them, for `og comphet
    //! phase`).

    use super::Mt19937;
    use std::fmt::Write as _;
    use std::fs;
    use std::path::Path;

    pub const SAMPLE: &str = "CHPH";
    pub const CONTIG: &str = "chrP";
    pub const LENGTH: usize = 80_000;
    const READ: usize = 100;
    const FRAG: usize = 500;
    const FIELDS: &str = "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|CANONICAL|MANE_SELECT|Existing_variation|MAX_AF|gnomADe_AF|gnomADg_AF|CLIN_SIG";

    /// One variant: 1-based position, alleles, consequence and the copy
    /// (0 or 1) that carries it.
    #[derive(Debug, Clone, PartialEq, Eq)]
    pub struct Variant {
        pub gene: &'static str,
        pub pos: usize,
        pub alt_ref: String,
        pub alt: String,
        pub consequence: &'static str,
        pub copy: u8,
    }

    /// The contig: `random.Random(31)` choosing among `ACGT` 80,000 times.
    pub fn reference() -> Vec<u8> {
        let mut rng = Mt19937::new(31);
        (0..LENGTH)
            .map(|_| b"ACGT"[rng.below(4) as usize])
            .collect()
    }

    /// The base before `b` in `ACGT`, wrapping (`A` gives `T`).
    fn other(b: u8) -> u8 {
        let i = b"ACGT"
            .iter()
            .position(|&x| x == b)
            .expect("a base of ACGT");
        b"ACGT"[(i + 3) % 4]
    }

    pub fn variants(seq: &[u8], swap: bool) -> Vec<Variant> {
        let s = |a: usize, b: usize| String::from_utf8_lossy(&seq[a..b]).into_owned();
        let snv = |gene, pos: usize, consequence, copy| Variant {
            gene,
            pos,
            alt_ref: s(pos - 1, pos),
            alt: (other(seq[pos - 1]) as char).to_string(),
            consequence,
            copy,
        };
        vec![
            snv("GENEC", 1000, "stop_gained", 0),
            Variant {
                gene: "GENEC",
                pos: 1003,
                alt_ref: s(1002, 1005),
                alt: s(1002, 1003),
                consequence: "frameshift_variant",
                copy: u8::from(swap),
            },
            snv("GENET", 5000, "stop_gained", 0),
            snv("GENET", 5400, "stop_gained", 1),
            snv("GENEF", 20000, "stop_gained", 0),
            snv("GENEF", 70000, "stop_gained", 1),
        ]
    }

    /// The sequence and CIGAR of a read of `READ` bases from 0-based
    /// `start` on one copy.
    fn read_on(seq: &[u8], vs: &[Variant], copy: u8, start: usize) -> (String, String) {
        let mine = || vs.iter().filter(|v| v.copy == copy);
        let snv_at = |p: usize| {
            mine()
                .find(|v| v.alt_ref.len() == v.alt.len() && v.pos - 1 == p)
                .map(|v| v.alt.as_bytes()[0])
        };
        // A deletion is keyed on the first deleted base.
        let del_at = |p: usize| {
            mine()
                .find(|v| v.alt_ref.len() > v.alt.len() && v.pos + v.alt.len() - 1 == p)
                .map(|v| v.alt_ref.len() - v.alt.len())
        };
        let (mut out, mut ops, mut p) =
            (Vec::with_capacity(READ), Vec::<(u8, usize)>::new(), start);
        while out.len() < READ {
            if let Some(n) = del_at(p) {
                ops.push((b'D', n));
                p += n;
                continue;
            }
            out.push(snv_at(p).unwrap_or(seq[p]));
            ops.push((b'M', 1));
            p += 1;
        }
        let mut cigar = String::new();
        let mut i = 0;
        while i < ops.len() {
            let (op, mut n) = ops[i];
            let mut j = i + 1;
            while j < ops.len() && ops[j].0 == op {
                n += ops[j].1;
                j += 1;
            }
            let _ = write!(cigar, "{n}{}", op as char);
            i = j;
        }
        (String::from_utf8_lossy(&out).into_owned(), cigar)
    }

    pub fn sam(seq: &[u8], vs: &[Variant]) -> String {
        let mut s = format!(
            "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:{CONTIG}\tLN:{LENGTH}\n@RG\tID:{SAMPLE}\tSM:{SAMPLE}\n"
        );
        // Fragment starts (0-based) around each pair: read 1 covers the
        // first position, its mate the second one when 400 bp apart.
        let starts: [(&str, std::iter::StepBy<std::ops::Range<usize>>); 4] = [
            ("c", (920..991).step_by(7)),
            ("t", (4905..4996).step_by(10)),
            ("f", (19905..19996).step_by(10)),
            ("g", (69905..69996).step_by(10)),
        ];
        let q = "I".repeat(READ);
        for (region, rs) in starts {
            for copy in 0..2u8 {
                for st in rs.clone() {
                    let m = st + FRAG - READ;
                    let (s1, c1) = read_on(seq, vs, copy, st);
                    let (s2, c2) = read_on(seq, vs, copy, m);
                    let name = format!("{region}{copy}_{st}");
                    let _ = writeln!(
                        s,
                        "{name}\t99\t{CONTIG}\t{}\t60\t{c1}\t=\t{}\t{FRAG}\t{s1}\t{q}\tRG:Z:{SAMPLE}",
                        st + 1,
                        m + 1
                    );
                    let _ = writeln!(
                        s,
                        "{name}\t147\t{CONTIG}\t{}\t60\t{c2}\t=\t{}\t-{FRAG}\t{s2}\t{q}\tRG:Z:{SAMPLE}",
                        m + 1,
                        st + 1
                    );
                }
            }
        }
        s
    }

    fn vcf_head(info: &str) -> String {
        format!(
            "##fileformat=VCFv4.2\n##contig=<ID={CONTIG},length={LENGTH}>\n##FILTER=<ID=PASS,Description=\"All filters passed\">\n{info}\
             ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
             ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">\n\
             ##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">\n\
             ##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allele depths\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{SAMPLE}\n"
        )
    }

    pub fn vep_vcf(vs: &[Variant]) -> String {
        let mut out = vcf_head(&format!(
            "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: {FIELDS}\">\n"
        ));
        for v in vs {
            let (g, p) = (v.gene, v.pos);
            let _ = writeln!(
                out,
                "{CONTIG}\t{p}\t.\t{}\t{}\t50\tPASS\tCSQ={}|{}|HIGH|{g}|ENSG{g}|Transcript|ENST{p}|protein_coding|YES|NM_{p}.1|rs{p}|0.0001|0.0001|0.0001|\tGT:GQ:DP:AD\t0/1:50:30:15,15",
                v.alt_ref, v.alt, v.alt, v.consequence
            );
        }
        out
    }

    pub fn comphet_vcf(vs: &[Variant]) -> String {
        let mut out = vcf_head(
            "##INFO=<ID=slivar_comphet,Number=.,Type=String,Description=\"compound hets called by slivar\">\n",
        );
        for (i, v) in vs.iter().enumerate() {
            let mate = if i % 2 == 0 { &vs[i + 1] } else { &vs[i - 1] };
            let pid = match v.gene {
                "GENEC" => 1,
                "GENET" => 2,
                _ => 3,
            };
            let _ = writeln!(
                out,
                "{CONTIG}\t{}\t.\t{}\t{}\t50\tPASS\tslivar_comphet={SAMPLE}/{}/{pid}/{CONTIG}/{}/{}/{}\tGT:GQ:DP:AD\t0/1:50:30:15,15",
                v.pos, v.alt_ref, v.alt, v.gene, mate.pos, mate.alt_ref, mate.alt
            );
        }
        out
    }

    /// Writes the three files into `dir`, made if missing.
    pub fn write(dir: &Path, swap: bool) -> Result<(), String> {
        let seq = reference();
        let vs = variants(&seq, swap);
        fs::create_dir_all(dir).map_err(|e| format!("{}: {e}", dir.display()))?;
        for (name, text) in [
            ("reads.sam", sam(&seq, &vs)),
            ("vep.vcf", vep_vcf(&vs)),
            ("comphet.vcf", comphet_vcf(&vs)),
        ] {
            let path = dir.join(name);
            crate::evidence::publish(&path, text.as_bytes())?;
        }
        Ok(())
    }
}

/// `x.get(k)`: the value, or JSON null.
fn get(x: &IndexMap<String, Value>, k: &str) -> Value {
    x.get(k).cloned().unwrap_or(Value::Null)
}

/// Python's truth of a JSON value.
fn truthy(v: &Value) -> bool {
    pycompat::value::truthy(v)
}

/// The old trimmer used `v or []`: null and every false value are empty.
/// `v or []` as a list; a truthy value that is not a list is an error.
fn list<'a>(v: Option<&'a Value>, what: &str) -> Result<&'a [Value], String> {
    match v {
        Some(Value::Array(a)) => Ok(a),
        Some(x) if truthy(x) => Err(format!("{what} is not a list")),
        _ => Ok(&[]),
    }
}

fn object<'a>(v: &'a Value, what: &str) -> Result<&'a IndexMap<String, Value>, String> {
    match v {
        Value::Object(o) => Ok(o),
        _ => Err(format!("{what} is not an object")),
    }
}

/// The keys of `x` among `keys`, in that order.
fn pick(x: &IndexMap<String, Value>, keys: &[&str]) -> IndexMap<String, Value> {
    keys.iter()
        .filter_map(|k| x.get(*k).map(|v| ((*k).to_owned(), v.clone())))
        .collect()
}

/// One listed diplotype with the fields `og pgx` reads.
fn diplotype(x: &Value) -> Result<Value, String> {
    let Value::Object(x) = x else {
        // The Python dip() explicitly passes non-object listed diplotypes through.
        return Ok(x.clone());
    };
    let mut out = pick(
        x,
        &[
            "gene",
            "label",
            "phenotypes",
            "activityScore",
            "outsidePhenotype",
        ],
    );
    for a in ["allele1", "allele2"] {
        let Some(v) = x.get(a) else { continue };
        let kept = match v {
            Value::Object(o) if !o.is_empty() => {
                Value::Object(pick(o, &["gene", "name", "function"]))
            }
            v if !truthy(v) => v.clone(),
            _ => return Err(format!("diplotype {a} is not an object")),
        };
        out.insert(a.to_owned(), kept);
    }
    Ok(Value::Object(out))
}

/// `items[:n]` with Python's slice rule for a negative `n`.
fn head(items: &[Value], n: i64) -> &[Value] {
    let len = items.len() as i64;
    let end = if n < 0 { (len + n).max(0) } else { n.min(len) };
    &items[..end as usize]
}

/// The trimmed report: every gene with its listed diplotypes and related
/// drugs, and for each drug of each guidance source up to `max` of the
/// annotations PharmCAT matched, with only the fields `og pgx` reads.
pub fn trim(d: Value, max: i64) -> Result<Value, String> {
    let Value::Object(mut d) = d else {
        return Err("the report is not an object".to_owned());
    };
    let mut genes = IndexMap::new();
    let Some(Value::Object(input_genes)) = d.swap_remove("genes") else {
        return Err("genes is not an object".to_owned());
    };
    // Consume one source gene at a time: do not retain a second full report.
    for (g, v) in input_genes {
        let v = object(&v, &format!("gene {g}"))?;
        let mut gene = pick(v, &["geneSymbol", "chr", "phased", "callSource"]);
        let mut related = Vec::new();
        for r in list(v.get("relatedDrugs"), "relatedDrugs")? {
            let r = object(r, "a related drug")?;
            let name = r.get("name").ok_or("a related drug has no name")?.clone();
            related.push(Value::Object(IndexMap::from([
                ("name".to_owned(), name),
                ("id".to_owned(), get(r, "id")),
            ])));
        }
        gene.insert("relatedDrugs".to_owned(), Value::Array(related));
        // sourceDiplotypes is what og pgx reads; recommendationDiplotypes
        // only when sourceDiplotypes is absent, so it is kept only then.
        let k = if v.get("sourceDiplotypes").is_some_and(truthy) {
            "sourceDiplotypes"
        } else {
            "recommendationDiplotypes"
        };
        if v.contains_key(k) {
            let dips = list(v.get(k), k)?
                .iter()
                .map(diplotype)
                .collect::<Result<_, _>>()?;
            gene.insert(k.to_owned(), Value::Array(dips));
        }
        genes.insert(g, Value::Object(gene));
    }
    let mut drugs: IndexMap<String, Value> = IndexMap::new();
    let Some(Value::Object(input_drugs)) = d.swap_remove("drugs") else {
        return Err("drugs is not an object".to_owned());
    };
    for (source, by) in input_drugs {
        for (name, rep) in object(&by, &format!("drugs of {source}"))? {
            let rep = object(rep, &format!("drug {name}"))?;
            let mut gls = Vec::new();
            for gl in list(rep.get("guidelines"), "guidelines")? {
                let gl = object(gl, "a guideline")?;
                let mut anns = Vec::new();
                for a in head(list(gl.get("annotations"), "annotations")?, max) {
                    let a = object(a, "an annotation")?;
                    let mut gts = Vec::new();
                    for gt in list(a.get("genotypes"), "genotypes")? {
                        let gt = object(gt, "a genotype")?;
                        let dips = list(gt.get("diplotypes"), "diplotypes")?
                            .iter()
                            .filter_map(|x| match x {
                                Value::Object(x) => Some(Value::Object(IndexMap::from([
                                    ("gene".to_owned(), get(x, "gene")),
                                    ("label".to_owned(), get(x, "label")),
                                ]))),
                                // The Python comprehension has `if isinstance(x, dict)`.
                                _ => None,
                            })
                            .collect();
                        gts.push(Value::Object(IndexMap::from([(
                            "diplotypes".to_owned(),
                            Value::Array(dips),
                        )])));
                    }
                    anns.push(Value::Object(IndexMap::from([
                        ("phenotypes".to_owned(), get(a, "phenotypes")),
                        ("classification".to_owned(), get(a, "classification")),
                        (
                            "drugRecommendation".to_owned(),
                            get(a, "drugRecommendation"),
                        ),
                        ("genotypes".to_owned(), Value::Array(gts)),
                    ])));
                }
                if !anns.is_empty() {
                    gls.push(Value::Object(IndexMap::from([
                        ("name".to_owned(), get(gl, "name")),
                        ("source".to_owned(), get(gl, "source")),
                        ("annotations".to_owned(), Value::Array(anns)),
                    ])));
                }
            }
            if !gls.is_empty() {
                let entry = drugs
                    .entry(source.clone())
                    .or_insert_with(|| Value::Object(IndexMap::new()));
                if let Value::Object(e) = entry {
                    e.insert(
                        name.clone(),
                        Value::Object(IndexMap::from([
                            ("name".to_owned(), get(rep, "name")),
                            ("id".to_owned(), get(rep, "id")),
                            ("guidelines".to_owned(), Value::Array(gls)),
                        ])),
                    );
                }
            }
        }
    }
    let mut out = pick(
        &d,
        &["title", "timestamp", "pharmcatVersion", "dataVersion"],
    );
    out.insert("genes".to_owned(), Value::Object(genes));
    out.insert("drugs".to_owned(), Value::Object(drugs));
    let mut note = String::new();
    let _ = write!(
        note,
        "Trimmed by og dev fixture trim-pharmcat: every gene, and up to {max} matched annotations per drug guideline, with only the fields og pgx reads."
    );
    out.insert("_fixture_note".to_owned(), Value::Str(note));
    Ok(Value::Object(out))
}

fn trim_pharmcat(src: &Path, dst: &Path, max: i64) -> Result<(), String> {
    let text = fs::read(src).map_err(|e| py_os_error(&e, &src.to_string_lossy()))?;
    let text = String::from_utf8(text).map_err(|_| format!("{}: not UTF-8", src.display()))?;
    let d = json::loads(&text).map_err(|e| format!("{}: {e}", src.display()))?;
    drop(text);
    let out = trim(d, max)?;
    let opts = Dump {
        separators: Some((",", ":")),
        ..Dump::default()
    };
    let text = json::dumps(&out, &opts).map_err(|e| e.to_string())?;
    crate::evidence::publish(dst, text.as_bytes())
}

#[cfg(test)]
mod tests {
    use super::*;


    #[test]
    fn comphet_validates_later_destinations_before_publishing() {
        let d = std::env::temp_dir().join(format!("og-comphet-preflight-{}", std::process::id()));
        let _ = fs::remove_dir_all(&d);
        fs::create_dir_all(d.join("comphet.vcf")).unwrap();
        for name in ["reads.sam", "vep.vcf"] { fs::write(d.join(name), b"PREVIOUS OUTPUT").unwrap(); }
        assert!(comphet::write(&d, false).is_err());
        for name in ["reads.sam", "vep.vcf"] { assert_eq!(fs::read(d.join(name)).unwrap(), b"PREVIOUS OUTPUT"); }
        assert_eq!(fs::read_dir(&d).unwrap().count(), 3);
        fs::remove_dir_all(d).unwrap();
    }


    #[test]
    fn downloader_preflights_destinations_and_cleans_handled_failures() {
        use std::os::unix::fs::{PermissionsExt, symlink};
        use std::process::Command;
        let d = std::env::temp_dir().join(format!("og-download-safety-{}", std::process::id()));
        let _ = fs::remove_dir_all(&d);
        fs::create_dir_all(d.join("bin")).unwrap();
        let stub = d.join("bin/curl");
        fs::write(&stub, "#!/bin/sh\nprintf called >> \"$CALLED\"\nwhile [ \"$1\" != -o ]; do shift; done\nprintf partial > \"$2\"\nexit 22\n").unwrap();
        fs::set_permissions(&stub, fs::Permissions::from_mode(0o755)).unwrap();
        let script = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../../tests/smoke/vntyper_testdata.sh");
        for fault in ["directory", "directory-link", "curl", "crc", "sha", "signal"] {
            let out = d.join(fault);
            fs::create_dir(&out).unwrap();
            let dst = out.join("negative.bam");
            match fault {
                "directory" => fs::create_dir(&dst).unwrap(),
                "directory-link" => symlink(&d, &dst).unwrap(),
                _ => {
                    fs::write(&dst, b"PREVIOUS OUTPUT").unwrap();
                    symlink("negative.bam", out.join("negative.bam.part")).unwrap();
                }
            }
            if fault == "crc" {
                fs::write(&stub, "#!/bin/sh\nprintf called >> \"$CALLED\"\nwhile [ \"$1\" != -o ]; do shift; done\ndd if=/dev/zero of=\"$2\" bs=3376349 count=1 2>/dev/null\n").unwrap();
            } else if fault == "sha" {
                // Supply an inflater with correct size but the wrong SHA.
                fs::write(d.join("bin/gzip"), "#!/bin/sh\nprintf wrong-sha\n").unwrap();
                fs::set_permissions(d.join("bin/gzip"), fs::Permissions::from_mode(0o755)).unwrap();
            } else if fault == "signal" {
                fs::remove_file(d.join("bin/gzip")).unwrap();
                fs::write(&stub, "#!/bin/sh\nkill -TERM \"$PPID\"\nexit 22\n").unwrap();
            }
            let called = d.join("called");
            let _ = fs::remove_file(&called);
            let result = Command::new("sh").arg(&script).arg(&out).arg("negative.bam")
                .env("PATH", format!("{}:{}", d.join("bin").display(), std::env::var("PATH").unwrap()))
                .env("CALLED", &called).output().unwrap();
            assert!(!result.status.success(), "{fault}");
            if fault.starts_with("directory") { assert!(!called.exists(), "must reject before curl: {fault}"); }
            else {
                assert_eq!(fs::read(&dst).unwrap(), b"PREVIOUS OUTPUT", "{fault}");
                let entries: Vec<_> = fs::read_dir(&out).unwrap().map(|e| e.unwrap().file_name()).collect();
                assert_eq!(entries.len(), 2, "{fault}: {entries:?}");
            }
        }
        fs::remove_dir_all(d).unwrap();
    }

    #[test]
    fn approved_strict_pharmcat_containers_fail_at_the_cli() {
        let d = std::env::temp_dir().join(format!("og-pharmcat-shapes-{}", std::process::id()));
        fs::create_dir_all(&d).unwrap();
        for (input, class) in [
            (r#"{"genes":{"G":{"sourceDiplotypes":"x"}},"drugs":{}}"#, "sourceDiplotypes is not a list"),
            (r#"{"genes":{"G":{"recommendationDiplotypes":"x"}},"drugs":{}}"#, "recommendationDiplotypes is not a list"),
            (r#"{"genes":{"G":{"sourceDiplotypes":[{"allele1":"x"}]}},"drugs":{}}"#, "allele1 is not an object"),
            (r#"{"genes":{"G":{"sourceDiplotypes":[{"allele2":["x"]}]}},"drugs":{}}"#, "allele2 is not an object"),
            (r#"{"genes":{},"drugs":{"S":{"D":{"guidelines":[{"annotations":[{"genotypes":[{"diplotypes":{"x":1}}]}]}]}}}}"#, "diplotypes is not a list"),
        ] {
            fs::write(d.join("in"), input).unwrap();
            fs::write(d.join("out"), b"old").unwrap();
            let args = vec!["trim-pharmcat".into(), d.join("in").to_str().unwrap().into(), d.join("out").to_str().unwrap().into()];
            let (mut out, mut err) = (Vec::new(), Vec::new());
            assert_eq!(cli(&args, &mut out, &mut err), 1);
            assert!(String::from_utf8_lossy(&err).contains(class), "{err:?}");
            assert_eq!(fs::read(d.join("out")).unwrap(), b"old");
        }
        fs::remove_dir_all(d).unwrap();
    }

    #[test]
    fn mt19937_matches_cpython() {
        // random.Random(31): [r.getrandbits(32) for _ in range(3)] and the
        // first ten of "".join(r.choice("ACGT") ...) from a fresh seed.
        let mut r = Mt19937::new(31);
        let first: Vec<u32> = (0..3).map(|_| r.next_u32()).collect();
        assert_eq!(first, [52_734_659, 2_016_822_975, 482_749_413]);
        assert_eq!(&comphet::reference()[..10], b"ATATCACACC");
    }

    #[test]
    fn swap_moves_only_the_deletion() {
        let seq = comphet::reference();
        let (a, b) = (
            comphet::variants(&seq, false),
            comphet::variants(&seq, true),
        );
        assert_eq!(a[1].copy, 0);
        assert_eq!(b[1].copy, 1);
        assert_eq!(a[0], b[0]);
        assert_eq!(a[2..], b[2..]);
    }

    #[test]
    fn a_read_over_the_deletion_has_a_d_in_its_cigar() {
        let seq = comphet::reference();
        let vs = comphet::variants(&seq, false);
        let sam = comphet::sam(&seq, &vs);
        let line = sam
            .lines()
            .find(|l| l.starts_with("c0_920\t99"))
            .expect("read c0_920");
        assert!(line.contains("\t83M2D17M\t"), "{line}");
        let line = sam
            .lines()
            .find(|l| l.starts_with("c1_920\t99"))
            .expect("read c1_920");
        assert!(line.contains("\t100M\t"), "{line}");
    }

    #[test]
    fn trim_keeps_the_fields_og_pgx_reads() {
        let d = json::loads(
            r#"{"title":"t","extra":1,"genes":{"G":{"geneSymbol":"G","x":1,"relatedDrugs":[{"name":"d"}],
            "sourceDiplotypes":[],"recommendationDiplotypes":[{"label":"*1/*2","allele1":{"name":"*1","y":2},"allele2":null,"z":3},"str"]}},
            "drugs":{"CPIC":{"d":{"name":"d","guidelines":[{"name":"g","annotations":[{"genotypes":[{"diplotypes":[{"gene":"G"},5]}]},{},{}]},{"annotations":[]}]}}}}"#,
        )
        .expect("json");
        let out = trim(d.clone(), 2).expect("trim");
        let text = json::dumps(
            &out,
            &Dump {
                separators: Some((",", ":")),
                ..Dump::default()
            },
        )
        .expect("dumps");
        assert!(text.starts_with(r#"{"title":"t","genes":{"G":{"geneSymbol":"G","relatedDrugs":[{"name":"d","id":null}],"recommendationDiplotypes":[{"label":"*1/*2","allele1":{"name":"*1"},"allele2":null},"str"]}}"#), "{text}");
        assert!(text.contains(r#""annotations":[{"phenotypes":null,"classification":null,"drugRecommendation":null,"genotypes":[{"diplotypes":[{"gene":"G","label":null}]}]},{"phenotypes":null"#), "{text}");
        assert_eq!(text.matches("drugRecommendation").count(), 2, "{text}");
        // A negative count is Python's slice: all but the last.
        let out = trim(d, -2).expect("trim");
        let text = json::dumps(&out, &Dump::default()).expect("dumps");
        assert_eq!(text.matches("drugRecommendation").count(), 1, "{text}");
    }
}
