//! `og dev traits-sources`: the sources of the curated traits panel
//! (`assets/traits.tsv`, step 42) against the reference and the public
//! records it was curated from.
//!
//! ```text
//! og dev traits-sources [--table TABLE] [--fasta FASTA [--min-fasta-rows N]] [--online]
//! ```
//!
//! - `--fasta` reads each row's base from FASTA through its `.fai` and
//!   names every row whose `ref` differs. Rows on a contig the FASTA lacks
//!   are not checked; fewer than `--min-fasta-rows` (default 1) checked rows
//!   is a failure, so a FASTA with none of the table's contigs cannot pass.
//! - `--online` looks every PubMed id up in PubMed (NCBI E-utilities) and
//!   every rsID up in dbSNP (NCBI Variation Services), with curl: an id
//!   that does not resolve, an rsID that is merged or withdrawn, or one
//!   whose GRCh38 position, reference base or alleles differ from the row
//!   is a failure.
//!
//! Exit 0 when every check passed, 1 when one failed, 2 for a refused
//! command line. The table's schema is `og traits check`'s to judge; run it
//! first. This was the `--fasta` and `--online` part of
//! `tests/test_traits_table.py`.

use std::collections::BTreeSet;
use std::io::{Read, Seek, SeekFrom, Write};
use std::path::Path;
use std::time::Duration;

use pycompat::json::{Value, loads};
use pycompat::value::truthy;

const USAGE: &str = "usage: og dev traits-sources [--table TABLE] [--fasta FASTA [--min-fasta-rows N]] [--online]\n";
const TABLE: &str = "assets/traits.tsv";
const PUBMED: &str =
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&id=";
const DBSNP: &str = "https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/";

/// One row of the table: the columns these checks read.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Row {
    pub chrom: String,
    pub pos: u64,
    pub rsid: String,
    pub ref_: String,
    pub alt: String,
    pub pmids: String,
}

/// The rows of a traits table: `#` lines and blank lines skipped, the
/// first other line the header.
pub fn read_table(text: &str) -> Result<Vec<Row>, String> {
    let mut lines = text
        .lines()
        .enumerate()
        .filter(|(_, l)| !l.starts_with('#') && !l.trim().is_empty());
    let Some((_, header)) = lines.next() else {
        return Err("no header line".into());
    };
    let header: Vec<&str> = header.split('\t').collect();
    let col = |name: &str| {
        header
            .iter()
            .position(|h| *h == name)
            .ok_or_else(|| format!("no {name} column"))
    };
    let [chrom, pos, rsid, ref_, alt, pmids] =
        ["chrom", "pos", "rsid", "ref", "alt", "pmids"].map(col);
    let (chrom, pos, rsid, ref_, alt, pmids) = (chrom?, pos?, rsid?, ref_?, alt?, pmids?);
    lines
        .map(|(i, l)| {
            let c: Vec<&str> = l.split('\t').collect();
            let get = |k: usize| c.get(k).copied().unwrap_or("").to_owned();
            Ok(Row {
                chrom: get(chrom),
                pos: get(pos)
                    .parse()
                    .map_err(|_| format!("line {}: pos '{}' is not a number", i + 1, get(pos)))?,
                rsid: get(rsid),
                ref_: get(ref_),
                alt: get(alt),
                pmids: get(pmids),
            })
        })
        .collect()
}

/// (problems, rows checked) of the reference bases against FASTA, read
/// through its `.fai`.
pub fn fasta_errors(fasta: &Path, rows: &[Row]) -> Result<(Vec<String>, usize), String> {
    let fai_path = format!("{}.fai", fasta.display());
    let fai = std::fs::read_to_string(&fai_path).map_err(|e| format!("{fai_path}: {e}"))?;
    let mut index = std::collections::HashMap::new();
    for l in fai.lines().filter(|l| !l.is_empty()) {
        let c: Vec<&str> = l.split('\t').collect();
        let num = |k: usize| -> Result<u64, String> {
            c.get(k)
                .and_then(|v| v.parse().ok())
                .ok_or_else(|| format!("{fai_path}: not a .fai line: {l}"))
        };
        index.insert(c[0], (num(1)?, num(2)?, num(3)?, num(4)?));
    }
    let mut fh = std::fs::File::open(fasta).map_err(|e| format!("{}: {e}", fasta.display()))?;
    let (mut errs, mut checked) = (Vec::new(), 0);
    for r in rows {
        let Some(&(length, offset, bases, width)) = index.get(r.chrom.as_str()) else {
            continue;
        };
        checked += 1;
        let i = r.pos.checked_sub(1).ok_or("pos 0")?;
        let got = if i >= length || bases == 0 {
            "past the end of the contig".to_owned()
        } else {
            let mut b = [0u8];
            fh.seek(SeekFrom::Start(offset + (i / bases) * width + i % bases))
                .and_then(|_| fh.read_exact(&mut b))
                .map_err(|e| format!("{}: {e}", fasta.display()))?;
            (b[0].to_ascii_uppercase() as char).to_string()
        };
        if got != r.ref_ {
            errs.push(format!(
                "{} {}:{}: table ref {}, FASTA {got}",
                r.rsid, r.chrom, r.pos, r.ref_
            ));
        }
    }
    Ok((errs, checked))
}

fn field<'a>(v: &'a Value, key: &str) -> Option<&'a Value> {
    match v {
        Value::Object(o) => o.get(key),
        _ => None,
    }
}

fn text(v: Option<&Value>) -> String {
    match v {
        Some(Value::Str(s)) => s.clone(),
        Some(Value::Int(i)) => i.to_string(),
        Some(Value::Null) | None => String::new(),
        Some(other) => format!("{other:?}"),
    }
}

/// The sorted, distinct PubMed ids of the table.
pub fn pmids(rows: &[Row]) -> Vec<String> {
    let ids: BTreeSet<(u64, String)> = rows
        .iter()
        .flat_map(|r| r.pmids.split(','))
        .map(str::trim)
        .filter(|p| !p.is_empty())
        .map(|p| (p.parse().unwrap_or(u64::MAX), p.to_owned()))
        .collect();
    ids.into_iter().map(|(_, p)| p).collect()
}

/// The PubMed ids that do not resolve in an esummary answer.
pub fn pubmed_errors(answer: &Value, ids: &[String]) -> Vec<String> {
    let result = field(answer, "result");
    ids.iter()
        .filter_map(|p| {
            let rec = result.and_then(|r| field(r, p));
            let error = rec.and_then(|r| field(r, "error"));
            let titled = rec.and_then(|r| field(r, "title")).is_some_and(truthy);
            (error.is_some() || !titled).then(|| {
                let why = if error.is_some() {
                    text(error)
                } else {
                    "no record".into()
                };
                format!("PMID {p} does not resolve in PubMed ({why})")
            })
        })
        .collect()
}

/// `NC_000001` to `NC_000022` as `chr1` to `chr22`.
fn contig(seq_id: &str) -> String {
    let acc = seq_id.split('.').next().unwrap_or(seq_id);
    acc.strip_prefix("NC_0000")
        .and_then(|n| n.parse::<u32>().ok())
        .filter(|n| (1..=22).contains(n) && acc.len() == 9)
        .map_or_else(|| seq_id.to_owned(), |n| format!("chr{n}"))
}

/// What is wrong with a row against its dbSNP record, if anything.
pub fn dbsnp_error(r: &Row, snp: &Value) -> Option<String> {
    let Some(psd) = field(snp, "primary_snapshot_data").filter(|v| truthy(v)) else {
        return Some(format!(
            "{}: not a current dbSNP record (merged or withdrawn)",
            r.rsid
        ));
    };
    let placements = match field(psd, "placements_with_allele") {
        Some(Value::Array(a)) => a.as_slice(),
        _ => &[],
    };
    let grch38 = |pl: &&Value| {
        let first = match field(pl, "placement_annot")
            .and_then(|a| field(a, "seq_id_traits_by_assembly"))
        {
            Some(Value::Array(t)) => t.first(),
            _ => None,
        };
        field(pl, "is_ptlp").is_some_and(truthy)
            && first
                .and_then(|t| field(t, "assembly_name"))
                .is_some_and(|n| matches!(n, Value::Str(s) if s.starts_with("GRCh38")))
    };
    let Some(place) = placements.iter().find(grch38) else {
        return Some(format!("{}: no GRCh38 placement in dbSNP", r.rsid));
    };
    let spdi: Vec<&Value> = match field(place, "alleles") {
        Some(Value::Array(a)) => a
            .iter()
            .filter_map(|a| field(a, "allele").and_then(|x| field(x, "spdi")))
            .collect(),
        _ => vec![],
    };
    let Some(first) = spdi.first() else {
        return Some(format!("{}: no allele in dbSNP's GRCh38 placement", r.rsid));
    };
    let chrom = contig(&text(field(place, "seq_id")));
    let pos = match field(first, "position") {
        Some(Value::Int(p)) => p + 1,
        _ => return Some(format!("{}: dbSNP's allele has no position", r.rsid)),
    };
    let ref_ = text(field(first, "deleted_sequence"));
    let alts: BTreeSet<String> = spdi
        .iter()
        .map(|s| {
            (
                text(field(s, "deleted_sequence")),
                text(field(s, "inserted_sequence")),
            )
        })
        .filter(|(d, i)| d != i)
        .map(|(_, i)| i)
        .collect();
    if (chrom.as_str(), pos, ref_.as_str())
        != (r.chrom.as_str(), i128::from(r.pos), r.ref_.as_str())
    {
        Some(format!(
            "{}: table {}:{} {}, dbSNP GRCh38 {chrom}:{pos} {ref_}",
            r.rsid, r.chrom, r.pos, r.ref_
        ))
    } else if !alts.contains(&r.alt) {
        Some(format!(
            "{}: alt {} is not one of dbSNP's alleles {alts:?}",
            r.rsid, r.alt
        ))
    } else {
        None
    }
}

/// A JSON answer by HTTP GET.
pub type Fetch<'a> = dyn FnMut(&str) -> Result<Value, String> + 'a;

/// Problems with the sources: PubMed ids that do not resolve, and rsIDs
/// whose GRCh38 placement, reference base or alleles differ from dbSNP's.
/// `pause` is the wait between two dbSNP requests (NCBI allows three a
/// second without a key).
pub fn online_errors(
    rows: &[Row],
    get: &mut Fetch<'_>,
    pause: Duration,
    out: &mut dyn Write,
) -> Vec<String> {
    let ids = pmids(rows);
    let mut errs = match get(&format!("{PUBMED}{}", ids.join(","))) {
        Ok(answer) => pubmed_errors(&answer, &ids),
        Err(e) => vec![format!("lookup failed: {e}")],
    };
    let _ = writeln!(out, "  {} PubMed ids looked up", ids.len());
    for r in rows {
        std::thread::sleep(pause);
        let id = r.rsid.strip_prefix("rs").unwrap_or(&r.rsid);
        match get(&format!("{DBSNP}{id}")) {
            Ok(snp) => errs.extend(dbsnp_error(r, &snp)),
            Err(e) => errs.push(format!("{}: dbSNP lookup failed: {e}", r.rsid)),
        }
    }
    let _ = writeln!(out, "  {} rsIDs looked up in dbSNP", rows.len());
    errs
}

/// GET `url` with curl, four tries, and parse the answer as JSON.
fn curl_json(url: &str) -> Result<Value, String> {
    let mut last = String::new();
    for t in 0..4u64 {
        if t > 0 {
            std::thread::sleep(Duration::from_secs(5 * t));
        }
        match std::process::Command::new("curl")
            .args(["-sS", "--fail", "--max-time", "60", url])
            .output()
        {
            Ok(o) if o.status.success() => match loads(&String::from_utf8_lossy(&o.stdout)) {
                Ok(v) => return Ok(v),
                Err(e) => last = format!("not JSON: {e:?}"),
            },
            Ok(o) => last = String::from_utf8_lossy(&o.stderr).trim().to_owned(),
            Err(e) => last = format!("curl: {e}"),
        }
    }
    Err(format!("{url}: {last}"))
}

struct Opts {
    table: String,
    fasta: Option<String>,
    min_fasta_rows: usize,
    online: bool,
}

fn parse(args: &[String]) -> Result<Opts, String> {
    let mut o = Opts {
        table: TABLE.into(),
        fasta: None,
        min_fasta_rows: 1,
        online: false,
    };
    let mut it = args.iter();
    while let Some(a) = it.next() {
        let mut value = |name: &str| {
            it.next()
                .cloned()
                .ok_or_else(|| format!("argument {name}: expected one argument"))
        };
        match a.as_str() {
            "--table" => o.table = value("--table")?,
            "--fasta" => o.fasta = Some(value("--fasta")?),
            "--min-fasta-rows" => {
                let v = value("--min-fasta-rows")?;
                o.min_fasta_rows = v
                    .parse()
                    .map_err(|_| format!("argument --min-fasta-rows: invalid int value: '{v}'"))?;
            }
            "--online" => o.online = true,
            other => return Err(format!("unrecognized arguments: {other}")),
        }
    }
    if o.fasta.is_none() && !o.online {
        return Err("nothing to check: give --fasta, --online or both".into());
    }
    Ok(o)
}

/// `og dev traits-sources ARGS`.
pub fn cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    if args.iter().any(|a| a == "-h" || a == "--help") {
        let _ = out.write_all(USAGE.as_bytes());
        return 0;
    }
    let o = match parse(args) {
        Ok(o) => o,
        Err(e) => {
            let _ = writeln!(err, "{USAGE}og dev traits-sources: error: {e}");
            return 2;
        }
    };
    let rows = match std::fs::read_to_string(&o.table)
        .map_err(|e| e.to_string())
        .and_then(|t| read_table(&t))
    {
        Ok(r) => r,
        Err(e) => {
            let _ = writeln!(err, "ERROR: {}: {e}", o.table);
            return 1;
        }
    };
    let mut fails = 0usize;
    let mut check = |out: &mut dyn Write, desc: &str, problems: &[String]| {
        if problems.is_empty() {
            let _ = writeln!(out, "[PASS] {desc}");
        } else {
            fails += 1;
            let _ = writeln!(out, "[FAIL] {desc}");
            for p in problems {
                let _ = writeln!(out, "  {p}");
            }
        }
    };
    if let Some(fasta) = &o.fasta {
        match fasta_errors(Path::new(fasta), &rows) {
            Ok((errs, n)) => {
                let _ = writeln!(
                    out,
                    "  ref checked against {fasta}: {n} rows (the others are on contigs it lacks)"
                );
                check(
                    out,
                    &format!("every reference base matches the FASTA ({n} rows checked)"),
                    &errs,
                );
                let few = (n < o.min_fasta_rows)
                    .then(|| format!("only {n} rows were checked"))
                    .into_iter()
                    .collect::<Vec<_>>();
                check(
                    out,
                    &format!(
                        "at least {} rows were checked against the FASTA",
                        o.min_fasta_rows
                    ),
                    &few,
                );
            }
            Err(e) => check(out, "the FASTA could be read", &[e]),
        }
    }
    if o.online {
        let errs = online_errors(&rows, &mut curl_json, Duration::from_millis(400), out);
        check(
            out,
            "every PMID resolves and every rsID matches dbSNP's GRCh38 record",
            &errs,
        );
    }
    let _ = writeln!(
        out,
        "\nRESULT: {}",
        if fails == 0 {
            "ALL PASS".to_owned()
        } else {
            format!("{fails} FAILED")
        }
    );
    u8::from(fails > 0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    const REPO: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../..");

    fn shipped() -> Vec<Row> {
        read_table(&std::fs::read_to_string(format!("{REPO}/{TABLE}")).unwrap()).unwrap()
    }

    fn scratch(name: &str) -> PathBuf {
        let d =
            std::env::temp_dir().join(format!("og-traits-sources-{}-{name}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    /// The test FASTA's chr2 length, and where the lactase row's base sits
    /// in it: line 3, column 5, so the check must use the .fai's line width.
    const LEN: u64 = 150;
    const POS: u64 = 125;

    /// A 150-base chr2 (60 bases a line, as its .fai says), N everywhere but
    /// the lactase row's base at POS, and the lactase row moved to POS. No
    /// test writes a chromosome; the e2e case traits-panel-1 checks the real
    /// coordinates against the real reference.
    fn lactase_fasta(d: &Path) -> (PathBuf, Row) {
        let lct = shipped()
            .into_iter()
            .find(|r| r.rsid == "rs4988235")
            .unwrap();
        let lct = Row { pos: POS, ..lct };
        let mut seq = vec![b'N'; usize::try_from(LEN).unwrap()];
        seq[usize::try_from(POS - 1).unwrap()] = lct.ref_.as_bytes()[0];
        let body: Vec<u8> = seq
            .chunks(60)
            .flat_map(|c| c.iter().copied().chain([b'\n']))
            .collect();
        let fa = d.join("ref.fasta");
        std::fs::write(&fa, [b">chr2\n".as_slice(), &body].concat()).unwrap();
        std::fs::write(d.join("ref.fasta.fai"), format!("chr2\t{LEN}\t6\t60\t61\n")).unwrap();
        (fa, lct)
    }

    #[test]
    fn the_test_fasta_is_small() {
        let d = scratch("small");
        let (fa, _) = lactase_fasta(&d);
        let size = std::fs::metadata(&fa).unwrap().len();
        assert!(size < 4096, "{size} bytes");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn fasta_check_reads_exactly_the_row_position_and_stops_at_the_contig_end() {
        let d = scratch("exact");
        let (fa, lct) = lactase_fasta(&d);
        let at = |pos: u64| Row { pos, ..lct.clone() };
        let want = |pos: u64, got: &str| {
            vec![format!(
                "rs4988235 chr2:{pos}: table ref {}, FASTA {got}",
                lct.ref_
            )]
        };
        // One base either side holds N: an offset off by one fails.
        for pos in [POS - 1, POS + 1] {
            assert_eq!(fasta_errors(&fa, &[at(pos)]).unwrap(), (want(pos, "N"), 1));
        }
        // The last base reads; one past the declared contig does not.
        assert_eq!(fasta_errors(&fa, &[at(LEN)]).unwrap(), (want(LEN, "N"), 1));
        assert_eq!(
            fasta_errors(&fa, &[at(LEN + 1)]).unwrap(),
            (want(LEN + 1, "past the end of the contig"), 1)
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn the_shipped_table_reads() {
        let rows = shipped();
        assert!(rows.len() >= 30, "{}", rows.len());
        assert!(
            rows.iter()
                .all(|r| r.rsid.starts_with("rs") && !r.pmids.is_empty())
        );
    }

    #[test]
    fn fasta_check_the_right_base_passes() {
        let d = scratch("right");
        let (fa, lct) = lactase_fasta(&d);
        assert_eq!(fasta_errors(&fa, &[lct]).unwrap(), (vec![], 1));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn fasta_check_control_a_planted_wrong_base_fails_naming_the_row() {
        let d = scratch("wrong");
        let (fa, lct) = lactase_fasta(&d);
        let other = ["A", "C", "G", "T"]
            .into_iter()
            .find(|b| *b != lct.ref_)
            .unwrap();
        let wrong = Row {
            ref_: other.into(),
            ..lct.clone()
        };
        let (errs, n) = fasta_errors(&fa, &[wrong]).unwrap();
        assert_eq!(n, 1);
        assert_eq!(
            errs,
            [format!(
                "rs4988235 chr2:{}: table ref {other}, FASTA {}",
                lct.pos, lct.ref_
            )]
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn fasta_check_a_row_on_a_contig_the_fasta_lacks_is_not_checked() {
        let d = scratch("contig");
        let (fa, _) = lactase_fasta(&d);
        let other = shipped().into_iter().find(|r| r.chrom != "chr2").unwrap();
        assert_eq!(fasta_errors(&fa, &[other]).unwrap(), (vec![], 0));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn the_command_line_fails_on_a_wrong_base_and_too_few_rows() {
        let d = scratch("cli");
        let (fa, lct) = lactase_fasta(&d);
        let fa = fa.to_string_lossy().into_owned();
        // The shipped table without its other chr2 rows: one row on the
        // FASTA's contig, moved to the test FASTA's POS, every other row on a
        // contig it lacks. With `wrong`, the lactase row's ref is that base.
        let pos = POS.to_string();
        let text = std::fs::read_to_string(format!("{REPO}/{TABLE}")).unwrap();
        let write = |name: &str, wrong: Option<&str>| {
            let kept: String = text
                .lines()
                .filter_map(|l| {
                    let mut c: Vec<&str> = l.split('\t').collect();
                    let lactase = c.get(2) == Some(&"rs4988235");
                    if c.first() == Some(&"chr2") && !lactase {
                        return None;
                    }
                    if lactase {
                        c[1] = &pos;
                    }
                    if let (true, Some(w)) = (lactase, wrong) {
                        c[4] = w;
                    }
                    Some(c.join("\t") + "\n")
                })
                .collect();
            let p = d.join(name);
            std::fs::write(&p, kept).unwrap();
            p.to_string_lossy().into_owned()
        };
        let table = write("one.tsv", None);
        let run = |args: &[&str]| {
            let args: Vec<String> = args.iter().map(|a| (*a).to_owned()).collect();
            let (mut out, mut err) = (Vec::new(), Vec::new());
            let code = cli(&args, &mut out, &mut err);
            (
                code,
                String::from_utf8(out).unwrap(),
                String::from_utf8(err).unwrap(),
            )
        };
        let (code, out, _) = run(&["--table", &table, "--fasta", &fa]);
        assert_eq!(code, 0, "{out}");
        assert!(
            out.contains(&format!("ref checked against {fa}: 1 rows")),
            "{out}"
        );
        let (code, out, _) = run(&["--table", &table, "--fasta", &fa, "--min-fasta-rows", "2"]);
        assert_eq!(code, 1, "{out}");
        assert!(out.contains("only 1 rows were checked"), "{out}");
        let wrong = if lct.ref_ == "C" { "T" } else { "C" };
        let planted = write("planted.tsv", Some(wrong));
        let (code, out, _) = run(&["--table", &planted, "--fasta", &fa]);
        assert_eq!(code, 1, "{out}");
        assert!(
            out.contains(&format!(
                "rs4988235 chr2:{}: table ref {wrong}, FASTA {}",
                lct.pos, lct.ref_
            )),
            "{out}"
        );
        assert!(out.contains("RESULT: 1 FAILED"), "{out}");
        let (code, _, err) = run(&["--table", &table]);
        assert_eq!(code, 2);
        assert!(err.contains("nothing to check"), "{err}");
        let (code, _, err) = run(&["--fasta", &fa, "--bogus"]);
        assert_eq!(code, 2);
        assert!(err.contains("unrecognized arguments: --bogus"), "{err}");
        let _ = std::fs::remove_dir_all(&d);
    }

    fn row(chrom: &str, pos: u64, ref_: &str, alt: &str) -> Row {
        Row {
            chrom: chrom.into(),
            pos,
            rsid: "rs1801133".into(),
            ref_: ref_.into(),
            alt: alt.into(),
            pmids: "7647779".into(),
        }
    }

    /// A dbSNP refsnp answer with one GRCh38 placement (0-based position)
    /// and one GRCh37 placement before it.
    fn snp(seq_id: &str, pos0: i64, del: &str, ins: &[&str]) -> Value {
        snp_on("GRCh38.p14", seq_id, pos0, del, ins)
    }

    fn snp_on(assembly: &str, seq_id: &str, pos0: i64, del: &str, ins: &[&str]) -> Value {
        let alleles: Vec<String> = [del]
            .iter()
            .chain(ins)
            .map(|i| {
                format!(
                    r#"{{"allele": {{"spdi": {{"seq_id": "{seq_id}", "position": {pos0}, "deleted_sequence": "{del}", "inserted_sequence": "{i}"}}}}}}"#
                )
            })
            .collect();
        loads(&format!(
            r#"{{"primary_snapshot_data": {{"placements_with_allele": [
              {{"seq_id": "NC_000001.10", "is_ptlp": false, "placement_annot": {{"seq_id_traits_by_assembly": [{{"assembly_name": "GRCh37.p13"}}]}}, "alleles": []}},
              {{"seq_id": "{seq_id}", "is_ptlp": true, "placement_annot": {{"seq_id_traits_by_assembly": [{{"assembly_name": "{assembly}"}}]}}, "alleles": [{}]}}]}}}}"#,
            alleles.join(", ")
        ))
        .unwrap()
    }

    #[test]
    fn dbsnp_a_matching_record_passes() {
        let s = snp("NC_000001.11", 11_796_320, "G", &["A"]);
        assert_eq!(dbsnp_error(&row("chr1", 11_796_321, "G", "A"), &s), None);
    }

    #[test]
    fn dbsnp_controls_another_position_base_or_allele_fails() {
        let s = snp("NC_000001.11", 11_796_320, "G", &["A"]);
        assert_eq!(
            dbsnp_error(&row("chr1", 11_796_322, "G", "A"), &s).as_deref(),
            Some("rs1801133: table chr1:11796322 G, dbSNP GRCh38 chr1:11796321 G")
        );
        assert!(dbsnp_error(&row("chr1", 11_796_321, "C", "A"), &s).is_some());
        assert_eq!(
            dbsnp_error(&row("chr1", 11_796_321, "G", "T"), &s).as_deref(),
            Some("rs1801133: alt T is not one of dbSNP's alleles {\"A\"}")
        );
        let merged = loads(r#"{"merged_snapshot_data": {}}"#).unwrap();
        assert_eq!(
            dbsnp_error(&row("chr1", 1, "G", "A"), &merged).as_deref(),
            Some("rs1801133: not a current dbSNP record (merged or withdrawn)")
        );
        let only37 = snp_on("GRCh37.p13", "NC_000001.10", 11_856_377, "G", &["A"]);
        assert_eq!(
            dbsnp_error(&row("chr1", 1, "G", "A"), &only37).as_deref(),
            Some("rs1801133: no GRCh38 placement in dbSNP")
        );
    }

    #[test]
    fn contigs_are_named_as_the_table_names_them() {
        assert_eq!(contig("NC_000001.11"), "chr1");
        assert_eq!(contig("NC_000022.11"), "chr22");
        assert_eq!(contig("NC_000023.11"), "NC_000023.11");
        assert_eq!(contig("NT_187361.1"), "NT_187361.1");
    }

    #[test]
    fn pubmed_an_id_without_a_record_or_with_an_error_fails() {
        let answer = loads(
            r#"{"result": {"uids": ["7647779", "99999999"], "7647779": {"title": "A candidate genetic risk factor"},
                "99999999": {"uid": "99999999", "error": "cannot get document summary"}}}"#,
        )
        .unwrap();
        let ids = vec!["7647779".to_owned(), "99999999".into(), "123".into()];
        assert_eq!(
            pubmed_errors(&answer, &ids),
            [
                "PMID 99999999 does not resolve in PubMed (cannot get document summary)",
                "PMID 123 does not resolve in PubMed (no record)"
            ]
        );
    }

    #[test]
    fn online_asks_pubmed_once_and_dbsnp_per_row() {
        let rows = vec![
            row("chr1", 11_796_321, "G", "A"),
            Row {
                rsid: "rs671".into(),
                pmids: "99999999,7647779".into(),
                ..row("chr12", 111_803_962, "G", "A")
            },
        ];
        let mut asked = Vec::new();
        let mut get = |url: &str| -> Result<Value, String> {
            asked.push(url.to_owned());
            if url.starts_with(PUBMED) {
                return loads(r#"{"result": {"7647779": {"title": "t"}}}"#)
                    .map_err(|e| format!("{e:?}"));
            }
            if url.ends_with("/671") {
                return Err("timeout".into());
            }
            Ok(snp("NC_000001.11", 11_796_320, "G", &["A"]))
        };
        let mut out = Vec::new();
        let errs = online_errors(&rows, &mut get, Duration::ZERO, &mut out);
        assert_eq!(
            errs,
            [
                "PMID 99999999 does not resolve in PubMed (no record)",
                "rs671: dbSNP lookup failed: timeout"
            ]
        );
        assert_eq!(asked[0], format!("{PUBMED}7647779,99999999"));
        assert_eq!(
            &asked[1..],
            [format!("{DBSNP}1801133"), format!("{DBSNP}671")]
        );
        let out = String::from_utf8(out).unwrap();
        assert!(
            out.contains("2 PubMed ids looked up") && out.contains("2 rsIDs looked up in dbSNP")
        );
    }
}
