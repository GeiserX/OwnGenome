//! `og table constraint-join`: gnomAD v4.1 gene constraint columns appended
//! to a table.
//!
//! ```text
//! og table constraint-join --gene-col NAME [--constrained] CONSTRAINT.tsv TABLE.tsv
//! ```
//!
//! `CONSTRAINT.tsv` is `gnomad.v4.1.constraint_metrics.tsv` (plain or gzip):
//! its columns are found by header name (`gene`, `canonical`, `transcript`,
//! `lof.oe_ci.upper`, `lof.pLI`, `mis.z_score`). Only canonical rows count;
//! v4.1 lists an Ensembl and a RefSeq canonical row per gene, and the
//! Ensembl (`ENST`) one wins wherever it sits in the file. `TABLE.tsv` (or
//! `-`, stdin) has a header row, and `--gene-col` names its gene-symbol
//! column.
//!
//! Out: the table with `LOEUF`, `pLI` and `mis_z` appended, and with
//! `--constrained` `CONSTRAINED`: `YES` when LOEUF < 0.35 or pLI > 0.9 (the
//! loss-of-function intolerant genes), `NO` when either value is known and
//! neither says so, `.` when the gene has neither. A missing value is `.`.
//! A summary line goes to stderr.
//!
//! Exit 3: the constraint table is empty, or a needed column is missing.
//! Exit 4: rows carry a gene symbol but not one matched, which means a wrong
//! or broken constraint file. Exit 2: usage, or a file cannot be read.
//!
//! This is `bin/constraint_join.awk` rewritten; its output, messages and
//! exit codes are the same (the parity cases in `tests/parity` hold it to
//! mawk), except that an empty constraint file is reported as empty: awk
//! read the table as the constraint file then and blamed a missing column.

use std::collections::{HashMap, HashSet};
use std::io::{self, BufRead, Write};

use memchr::memchr_iter;

/// The constraint table's columns this reads, in the order they are checked.
const NEED: [&str; 6] = [
    "gene",
    "canonical",
    "transcript",
    "lof.oe_ci.upper",
    "lof.pLI",
    "mis.z_score",
];

/// One gene's values, each as written, or `.` for `NA` and empty.
#[derive(Debug, Clone, PartialEq, Eq)]
struct Gene {
    loeuf: Vec<u8>,
    pli: Vec<u8>,
    mis_z: Vec<u8>,
    /// From an Ensembl (`ENST`) transcript row.
    ensembl: bool,
}

/// Why the join stopped, with its exit code.
#[derive(Debug)]
pub enum Error {
    /// A message for stderr; exit 3 or 4.
    Data(u8, String),
    /// Reading or writing failed; exit 2 (the file is named in the text).
    Io(String),
}

impl Error {
    pub fn code(&self) -> u8 {
        match self {
            Error::Data(code, _) => *code,
            Error::Io(_) => 2,
        }
    }
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::Data(_, m) | Error::Io(m) => f.write_str(m),
        }
    }
}

/// The tab-separated fields of a line; an empty line has none.
fn fields(line: &[u8]) -> Vec<&[u8]> {
    if line.is_empty() {
        return Vec::new();
    }
    let mut out = Vec::with_capacity(16);
    let mut start = 0;
    for t in memchr_iter(b'\t', line) {
        out.push(&line[start..t]);
        start = t + 1;
    }
    out.push(&line[start..]);
    out
}

/// Field `i` (0-based), empty past the end.
fn get<'a>(f: &[&'a [u8]], i: usize) -> &'a [u8] {
    f.get(i).copied().unwrap_or(b"")
}

/// Reads one line without its `\n` into `buf`; false at the end.
fn next_line(input: &mut impl BufRead, buf: &mut Vec<u8>) -> io::Result<bool> {
    buf.clear();
    if input.read_until(b'\n', buf)? == 0 {
        return Ok(false);
    }
    if buf.last() == Some(&b'\n') {
        buf.pop();
    }
    Ok(true)
}

/// A number as awk reads a string (C's `strtod` on its longest numeric
/// prefix): leading white space, a sign, decimal or `0x` hexadecimal
/// digits, `inf`, `infinity` and `nan`; 0 when there is no number.
pub fn awk_number(s: &[u8]) -> f64 {
    let mut i = 0;
    while i < s.len() && matches!(s[i], b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r') {
        i += 1;
    }
    let mut sign = 1.0;
    if let Some(&c @ (b'+' | b'-')) = s.get(i) {
        if c == b'-' {
            sign = -1.0;
        }
        i += 1;
    }
    let rest = &s[i..];
    let lower = |n: usize| rest.get(..n).map(<[u8]>::to_ascii_lowercase);
    if lower(3).as_deref() == Some(b"inf") {
        return sign * f64::INFINITY;
    }
    if lower(3).as_deref() == Some(b"nan") {
        return f64::NAN;
    }
    let hex_digit = |c: u8| c.is_ascii_hexdigit();
    if let Some(h) = rest.strip_prefix(b"0x").or(rest.strip_prefix(b"0X"))
        && (h.first().copied().is_some_and(hex_digit)
            || (h.first() == Some(&b'.') && h.get(1).copied().is_some_and(hex_digit)))
    {
        return sign * hex_number(h);
    }
    // Decimal: digits, a dot, digits (one digit at least), an exponent.
    let digits = |from: usize| {
        rest[from..]
            .iter()
            .take_while(|c| c.is_ascii_digit())
            .count()
    };
    let int = digits(0);
    let mut end = int;
    let mut frac = 0;
    if rest.get(end) == Some(&b'.') {
        frac = digits(end + 1);
        if int + frac > 0 {
            end += 1 + frac;
        }
    }
    if int + frac == 0 {
        return 0.0;
    }
    if let Some(b'e' | b'E') = rest.get(end) {
        let mut k = end + 1;
        if let Some(b'+' | b'-') = rest.get(k) {
            k += 1;
        }
        let n = digits(k);
        if n > 0 {
            end = k + n;
        }
    }
    let text = std::str::from_utf8(&rest[..end]).unwrap_or("0");
    // Rust reads `5.` and `.5` as C does.
    sign * text.parse::<f64>().unwrap_or(0.0)
}

/// The value of hexadecimal digits with an optional fraction and `p`
/// exponent, after `0x`.
fn hex_number(h: &[u8]) -> f64 {
    let mut mantissa = 0.0f64;
    let mut exp: i32 = 0;
    let mut i = 0;
    let mut seen_dot = false;
    while let Some(&c) = h.get(i) {
        if c == b'.' && !seen_dot {
            seen_dot = true;
        } else if let Some(d) = (c as char).to_digit(16) {
            mantissa = mantissa * 16.0 + f64::from(d);
            if seen_dot {
                exp -= 4;
            }
        } else {
            break;
        }
        i += 1;
    }
    if let Some(b'p' | b'P') = h.get(i) {
        let mut k = i + 1;
        let mut neg = false;
        if let Some(&c @ (b'+' | b'-')) = h.get(k) {
            neg = c == b'-';
            k += 1;
        }
        let digits: Vec<u8> = h[k..]
            .iter()
            .copied()
            .take_while(u8::is_ascii_digit)
            .collect();
        if !digits.is_empty() {
            let p: i32 = std::str::from_utf8(&digits)
                .ok()
                .and_then(|t| t.parse().ok())
                .unwrap_or(i32::MAX / 2);
            exp = exp.saturating_add(if neg { -p } else { p });
        }
    }
    mantissa * 2f64.powi(exp)
}

/// Reads the constraint table: each gene's canonical values.
fn load(mut input: impl BufRead, name: &str) -> Result<HashMap<Vec<u8>, Gene>, Error> {
    let io_err = |e: io::Error| Error::Io(format!("og table constraint-join: {name}: {e}"));
    let mut line = Vec::new();
    if !next_line(&mut input, &mut line).map_err(io_err)? {
        return Err(Error::Data(
            3,
            format!("ERROR: the constraint table {name} is empty"),
        ));
    }
    let mut col: HashMap<&[u8], usize> = HashMap::new();
    let head = line.clone();
    for (i, f) in fields(&head).into_iter().enumerate() {
        col.insert(f, i);
    }
    let mut at = [0usize; 6];
    for (k, need) in NEED.iter().enumerate() {
        at[k] = *col.get(need.as_bytes()).ok_or_else(|| {
            Error::Data(
                3,
                format!("ERROR: column {need} missing from the constraint table"),
            )
        })?;
    }
    let [gene, canonical, transcript, loeuf, pli, mis_z] = at;
    let mut genes: HashMap<Vec<u8>, Gene> = HashMap::with_capacity(32_768);
    let value = |v: &[u8]| -> Vec<u8> {
        if v.is_empty() || v == b"NA" {
            b".".to_vec()
        } else {
            v.to_vec()
        }
    };
    while next_line(&mut input, &mut line).map_err(io_err)? {
        let f = fields(&line);
        if get(&f, canonical) != b"true" {
            continue;
        }
        let g = get(&f, gene);
        if g.is_empty() {
            continue;
        }
        let ensembl = get(&f, transcript).starts_with(b"ENST");
        if genes.get(g).is_some_and(|old| old.ensembl || !ensembl) {
            continue;
        }
        genes.insert(
            g.to_vec(),
            Gene {
                loeuf: value(get(&f, loeuf)),
                pli: value(get(&f, pli)),
                mis_z: value(get(&f, mis_z)),
                ensembl,
            },
        );
    }
    Ok(genes)
}

/// The counts the summary line reports.
#[derive(Debug, Default, PartialEq, Eq)]
pub struct Summary {
    pub rows: u64,
    pub with_gene: u64,
    pub matched: usize,
}

/// Joins: reads `constraint` (named `constraint_name` in messages) and
/// `table`, writes the joined table to `out`. The summary is for stderr.
pub fn join(
    constraint: impl BufRead,
    constraint_name: &str,
    mut table: impl BufRead,
    table_name: &str,
    gene_col: &str,
    constrained: bool,
    mut out: impl Write,
) -> Result<Summary, Error> {
    let genes = load(constraint, constraint_name)?;
    let read_err = |e: io::Error| Error::Io(format!("og table constraint-join: {table_name}: {e}"));
    let write_err = |e: io::Error| Error::Io(format!("og table constraint-join: writing: {e}"));
    let mut line = Vec::new();
    let mut sum = Summary::default();
    if !next_line(&mut table, &mut line).map_err(read_err)? {
        return Ok(sum);
    }
    let Some(gc) = fields(&line)
        .iter()
        .rposition(|f| *f == gene_col.as_bytes())
    else {
        return Err(Error::Data(
            3,
            format!("ERROR: the table has no column named '{gene_col}'"),
        ));
    };
    let tail: &[u8] = if constrained {
        b"\tLOEUF\tpLI\tmis_z\tCONSTRAINED\n"
    } else {
        b"\tLOEUF\tpLI\tmis_z\n"
    };
    out.write_all(&line).map_err(write_err)?;
    out.write_all(tail).map_err(write_err)?;
    let mut matched: HashSet<&[u8]> = HashSet::new();
    while next_line(&mut table, &mut line).map_err(read_err)? {
        sum.rows += 1;
        let g = get(&fields(&line), gc).to_vec();
        if g != b"." && !g.is_empty() {
            sum.with_gene += 1;
        }
        out.write_all(&line).map_err(write_err)?;
        match genes.get_key_value(g.as_slice()) {
            Some((key, v)) => {
                matched.insert(key.as_slice());
                for x in [&v.loeuf, &v.pli, &v.mis_z] {
                    out.write_all(b"\t").map_err(write_err)?;
                    out.write_all(x).map_err(write_err)?;
                }
                if constrained {
                    let known = |x: &[u8]| x != b".";
                    let call: &[u8] = if (known(&v.loeuf) && awk_number(&v.loeuf) < 0.35)
                        || (known(&v.pli) && awk_number(&v.pli) > 0.9)
                    {
                        b"YES"
                    } else if known(&v.loeuf) || known(&v.pli) {
                        b"NO"
                    } else {
                        b"."
                    };
                    out.write_all(b"\t").map_err(write_err)?;
                    out.write_all(call).map_err(write_err)?;
                }
            }
            None => {
                let none: &[u8] = if constrained {
                    b"\t.\t.\t.\t."
                } else {
                    b"\t.\t.\t."
                };
                out.write_all(none).map_err(write_err)?;
            }
        }
        out.write_all(b"\n").map_err(write_err)?;
    }
    out.flush().map_err(write_err)?;
    sum.matched = matched.len();
    Ok(sum)
}

const USAGE: &str =
    "usage: og table constraint-join --gene-col NAME [--constrained] CONSTRAINT.tsv TABLE.tsv\n";

/// `og table constraint-join ARGS`: the exit code.
pub fn cli(args: &[String]) -> u8 {
    let usage = |msg: &str| {
        eprint!("{USAGE}");
        eprintln!("og table constraint-join: error: {msg}");
        2
    };
    let mut gene_col = None;
    let mut constrained = false;
    let mut files = Vec::new();
    let mut it = args.iter();
    while let Some(a) = it.next() {
        match a.as_str() {
            "-h" | "--help" => {
                print!("{USAGE}");
                return 0;
            }
            "--constrained" => constrained = true,
            "--gene-col" => match it.next() {
                Some(v) => gene_col = Some(v.clone()),
                None => return usage("argument --gene-col: expected one argument"),
            },
            s if s.starts_with("--gene-col=") => {
                gene_col = Some(s["--gene-col=".len()..].to_owned())
            }
            s if s.starts_with("--") => return usage(&format!("unrecognized arguments: {s}")),
            _ => files.push(a.clone()),
        }
    }
    let Some(gene_col) = gene_col.filter(|g| !g.is_empty()) else {
        return usage("the following arguments are required: --gene-col");
    };
    let [constraint, table] = files.as_slice() else {
        return usage("give CONSTRAINT.tsv and TABLE.tsv");
    };
    let open = |p: &str| {
        og_io::open(p).map_err(|e| {
            eprintln!("og table constraint-join: {p}: {e}");
            2u8
        })
    };
    let c = match open(constraint) {
        Ok(r) => r,
        Err(code) => return code,
    };
    let t = match open(table) {
        Ok(r) => r,
        Err(code) => return code,
    };
    let out = io::BufWriter::with_capacity(1 << 16, io::stdout().lock());
    match join(c, constraint, t, table, &gene_col, constrained, out) {
        Ok(s) => {
            eprintln!(
                "gnomAD constraint: {} of {} rows carry a gene symbol; {} distinct genes matched",
                s.with_gene, s.rows, s.matched
            );
            if s.with_gene > 0 && s.matched == 0 {
                eprintln!(
                    "ERROR: no gene in the table matched the constraint file; check that it is the gnomAD v4.1 constraint table"
                );
                return 4;
            }
            0
        }
        Err(e) => {
            eprintln!("{e}");
            e.code()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const CONSTRAINT: &str = "gene\tgene_id\ttranscript\tcanonical\tmane_select\tlof.oe_ci.upper\tlof.pLI\tmis.z_score\tsyn.z_score
GENEA\tENSG1\tENST00000000002\tfalse\tfalse\t0.95\t0.01\t0.40\t0.2
GENEA\tENSG1\tENST00000000001\ttrue\ttrue\t0.21\t0.99\t3.52\t0.1
GENEB\t1001\tNM_000003.4\ttrue\tfalse\t0.80\t0.00\t1.10\t0.3
GENEB\tENSG2\tENST00000000004\ttrue\ttrue\t0.52\t0.20\t2.25\t0.4
GENEC\tENSG3\tENST00000000005\ttrue\ttrue\tNA\tNA\tNA\t0.5
GENEC\t1002\tNM_000006.2\ttrue\tfalse\t0.70\t0.10\t1.50\t0.6
";
    const TABLE: &str = "CHROM\tPOS\tGENE\nchr1\t100\tGENEA\nchr2\t200\tGENEB\nchr3\t300\tGENEC\nchr4\t400\tGENED\nchr5\t500\t.\n";

    fn run(constraint: &str, table: &str, constrained: bool) -> Result<(String, Summary), Error> {
        let mut out = Vec::new();
        let s = join(
            constraint.as_bytes(),
            "c.tsv",
            table.as_bytes(),
            "t.tsv",
            "GENE",
            constrained,
            &mut out,
        )?;
        Ok((String::from_utf8(out).unwrap(), s))
    }

    #[test]
    fn canonical_ensembl_rows_win_and_constrained_is_called() {
        let (out, s) = run(CONSTRAINT, TABLE, true).unwrap();
        assert_eq!(
            out,
            "CHROM\tPOS\tGENE\tLOEUF\tpLI\tmis_z\tCONSTRAINED
chr1\t100\tGENEA\t0.21\t0.99\t3.52\tYES
chr2\t200\tGENEB\t0.52\t0.20\t2.25\tNO
chr3\t300\tGENEC\t.\t.\t.\t.
chr4\t400\tGENED\t.\t.\t.\t.
chr5\t500\t.\t.\t.\t.\t.
"
        );
        assert_eq!(
            s,
            Summary {
                rows: 5,
                with_gene: 4,
                matched: 3
            }
        );
    }

    #[test]
    fn without_constrained_three_columns() {
        let (out, _) = run(CONSTRAINT, TABLE, false).unwrap();
        assert!(out.starts_with(
            "CHROM\tPOS\tGENE\tLOEUF\tpLI\tmis_z\nchr1\t100\tGENEA\t0.21\t0.99\t3.52\n"
        ));
    }

    #[test]
    fn errors_name_the_problem() {
        let old = CONSTRAINT.replace("mis.z_score", "missense.z_score");
        let e = run(&old, TABLE, false).unwrap_err();
        assert_eq!(e.code(), 3);
        assert_eq!(
            e.to_string(),
            "ERROR: column mis.z_score missing from the constraint table"
        );
        let e = run("", TABLE, false).unwrap_err();
        assert_eq!(
            (e.code(), e.to_string().as_str()),
            (3, "ERROR: the constraint table c.tsv is empty")
        );
        let e = run(CONSTRAINT, "CHROM\tSYMBOL\nchr1\tA\n", false).unwrap_err();
        assert_eq!(
            (e.code(), e.to_string().as_str()),
            (3, "ERROR: the table has no column named 'GENE'")
        );
        // No match is for the caller to judge (exit 4 in the cli).
        let (_, s) = run(CONSTRAINT, "GENE\nOTHER1\nOTHER2\n", false).unwrap();
        assert_eq!((s.with_gene, s.matched), (2, 0));
        // An empty table: no output, nothing counted.
        assert_eq!(run(CONSTRAINT, "", true).unwrap().0, "");
    }

    #[test]
    fn numbers_read_as_awk_reads_them() {
        for (s, want) in [
            ("0.2", 0.2),
            (" 0.2", 0.2),
            ("1e", 1.0),
            (".5e1x", 5.0),
            ("2.", 2.0),
            ("abc", 0.0),
            ("", 0.0),
            ("-0.5", -0.5),
            ("+3", 3.0),
            ("0x10", 16.0),
            ("0X1p3", 8.0),
            ("0x.8", 0.5),
            ("0x", 0.0),
            ("1e400", f64::INFINITY),
            ("-inf", f64::NEG_INFINITY),
            ("Infinity", f64::INFINITY),
        ] {
            assert_eq!(awk_number(s.as_bytes()), want, "{s:?}");
        }
        assert!(awk_number(b"nan").is_nan());
    }
}
