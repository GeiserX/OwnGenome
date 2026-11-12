//! `og consumer check`: the strict reader of the 23andMe and AncestryDNA
//! layouts, on a `.txt` file or a `.zip` that holds one. Every problem is
//! listed, one per line, with exit 1; a file that keeps its layout gets one
//! line and exit 0. The export runs the same reader on what it writes.

use rustc_hash::FxHashSet;

use super::{Args, Usage, chrom_index, head, py_repr, py_repr_list, zip};

pub const HEADER_23ANDME: &str = "# rsid\tchromosome\tposition\tgenotype";
pub const HEADER_ANCESTRY: &str = "rsid\tchromosome\tposition\tallele1\tallele2";
pub const ANCESTRY_FIRST: &str = "#AncestryDNA raw data download";

/// The two layouts.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Layout {
    TwentyThree,
    Ancestry,
}

impl Layout {
    fn name(self) -> &'static str {
        match self {
            Layout::TwentyThree => "23andme",
            Layout::Ancestry => "ancestrydna",
        }
    }
}

/// `^(rs|i)[1-9][0-9]*$`.
fn is_rsid(s: &str) -> bool {
    s.strip_prefix("rs")
        .or_else(|| s.strip_prefix('i'))
        .is_some_and(is_pos)
}

/// `^[1-9][0-9]*$`.
fn is_pos(s: &str) -> bool {
    let b = s.as_bytes();
    matches!(b.first(), Some(b'1'..=b'9')) && b.iter().all(u8::is_ascii_digit)
}

fn is_base(s: &str) -> bool {
    matches!(s, "A" | "C" | "G" | "T")
}

/// What is wrong with a 23andMe row of four columns, if anything.
fn row_23andme(c: &[&str]) -> Option<&'static str> {
    if !is_rsid(c[0]) {
        return Some("rsid is not rs<number> or i<number>");
    }
    if chrom_index(c[1]).is_none() {
        return Some("chromosome is not 1-22, X, Y or MT");
    }
    if !is_pos(c[2]) {
        return Some("position is not a positive whole number");
    }
    let g = c[3];
    if g == "--" {
        return None;
    }
    if !(1..=2).contains(&g.chars().count()) || !g.bytes().all(|b| b"ACGT".contains(&b)) {
        return Some("genotype is not one or two of A, C, G, T, or --");
    }
    if g.len() == 1 && !matches!(c[1], "X" | "Y" | "MT") {
        return Some("one letter on an autosome");
    }
    if g.len() == 2 && c[1] == "MT" {
        return Some("two letters on MT");
    }
    None
}

/// What is wrong with an AncestryDNA row of five columns, if anything.
fn row_ancestry(c: &[&str]) -> Option<&'static str> {
    if !is_rsid(c[0]) {
        return Some("rsid is not rs<number> or i<number>");
    }
    let code_ok = c[1]
        .parse::<u8>()
        .is_ok_and(|n| (1..=26).contains(&n) && c[1] == n.to_string());
    if !code_ok {
        return Some("chromosome is not a code from 1 to 26");
    }
    if !is_pos(c[2]) {
        return Some("position is not a positive whole number");
    }
    if (c[3], c[4]) == ("0", "0") {
        return None;
    }
    if !is_base(c[3]) || !is_base(c[4]) {
        return Some("an allele is not A, C, G or T (and not a 0 0 no-call)");
    }
    None
}

/// The rsIDs seen so far. A row that passed its checks has an rsID of
/// `rs` or `i` and digits, so up to 18 digits it is kept as one number
/// (bit 62 marks `i`); a longer one is kept as text.
#[derive(Default)]
struct Ids {
    short: FxHashSet<u64>,
    long: FxHashSet<Box<str>>,
}

impl Ids {
    /// False when `id` was seen before.
    fn insert(&mut self, id: &str) -> bool {
        let (digits, tag) = match id.strip_prefix("rs") {
            Some(d) => (d, 0),
            None => (id.strip_prefix('i').unwrap_or(id), 1u64 << 62),
        };
        match digits.parse::<u64>() {
            Ok(n) if digits.len() <= 18 => self.short.insert(n | tag),
            _ => self.long.insert(id.into()),
        }
    }
}

/// Where the reader is in the file.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Part {
    /// The comment lines before the column header.
    Comments,
    Rows,
    /// The column header is missing: the rest is not read.
    NoHeader,
}

/// The strict reader, fed one line at a time, so a file is checked as it is
/// written and never held whole. [`check_text`] feeds it a whole text.
pub struct Checker<'p> {
    path: &'p str,
    layout: Layout,
    /// Lines fed so far.
    n: usize,
    part: Part,
    ascii: bool,
    carriage_return: bool,
    /// Whether the last line fed ended with a newline.
    newline: bool,
    first_ok: bool,
    rows: usize,
    problems: Vec<String>,
    /// The last good row: its chromosome, the length and digits of its
    /// position (a whole number without leading zeros, compared by value),
    /// and its rsID.
    last: Option<(usize, usize)>,
    last_pos: String,
    last_id: String,
    out_of_order: Option<String>,
    seen: Ids,
    twice: Option<String>,
}

impl<'p> Checker<'p> {
    /// A reader of `layout` that names `path` in each problem; `rows` is
    /// how many rows to make room for.
    pub fn new(path: &'p str, layout: Layout, rows: usize) -> Checker<'p> {
        let mut seen = Ids::default();
        seen.short.reserve(rows);
        Checker {
            path,
            layout,
            n: 0,
            part: Part::Comments,
            ascii: true,
            carriage_return: false,
            newline: false,
            first_ok: false,
            rows: 0,
            problems: Vec::new(),
            last: None,
            last_pos: String::new(),
            last_id: String::new(),
            out_of_order: None,
            seen,
            twice: None,
        }
    }

    /// Reads one line, with its newline if it has one; only the last line
    /// of a file may lack it.
    pub fn line(&mut self, raw: &[u8]) {
        if !self.ascii {
            return;
        }
        let Some(raw) = raw
            .is_ascii()
            .then(|| std::str::from_utf8(raw).ok())
            .flatten()
        else {
            self.ascii = false;
            return;
        };
        self.n += 1;
        self.newline = raw.ends_with('\n');
        let line = raw.strip_suffix('\n').unwrap_or(raw);
        if line.contains('\r') {
            self.carriage_return = true;
        }
        let path = self.path;
        let (header, ncols) = match self.layout {
            Layout::TwentyThree => (HEADER_23ANDME, 4),
            Layout::Ancestry => (HEADER_ANCESTRY, 5),
        };
        if self.n == 1 {
            self.first_ok = match self.layout {
                Layout::TwentyThree => line.starts_with("# This data file"),
                Layout::Ancestry => self.newline && line == ANCESTRY_FIRST,
            };
        }
        match self.part {
            Part::NoHeader => return,
            Part::Comments => {
                let comment = match self.layout {
                    Layout::TwentyThree => line.starts_with('#') && line != HEADER_23ANDME,
                    Layout::Ancestry => line.starts_with('#'),
                };
                if !comment {
                    self.part = if line == header {
                        Part::Rows
                    } else {
                        self.problems.push(format!(
                            "{path}: the column header {} is missing after the comment lines",
                            py_repr(header)
                        ));
                        Part::NoHeader
                    };
                }
                return;
            }
            Part::Rows => {}
        }
        let n = self.n;
        let mut cols = [""; 5];
        let mut k = 0;
        for f in line.split('\t') {
            if let Some(slot) = cols.get_mut(k) {
                *slot = f;
            }
            k += 1;
        }
        if k != ncols {
            self.problems.push(format!(
                "{path} line {n}: {k} columns, want {ncols}: {}",
                py_repr(head(line, 60))
            ));
            return;
        }
        let c = &cols[..k];
        let err = match self.layout {
            Layout::TwentyThree => row_23andme(c),
            Layout::Ancestry => row_ancestry(c),
        };
        if let Some(e) = err {
            self.problems
                .push(format!("{path} line {n}: {e}: {}", py_repr(head(line, 60))));
            return;
        }
        self.rows += 1;
        let chrom = match self.layout {
            Layout::TwentyThree => usize::from(chrom_index(c[1]).unwrap_or(0)),
            Layout::Ancestry => c[1].parse().unwrap_or(0),
        };
        let key = (chrom, c[2].len());
        if self.out_of_order.is_none()
            && let Some(last) = self.last
            && (key, c[2]) < (last, self.last_pos.as_str())
        {
            self.out_of_order = Some(format!(
                "{path}: rows out of order at {} (after {})",
                c[0], self.last_id
            ));
        }
        self.last = Some(key);
        self.last_pos.clear();
        self.last_pos.push_str(c[2]);
        self.last_id.clear();
        self.last_id.push_str(c[0]);
        if self.twice.is_none() && !self.seen.insert(c[0]) {
            self.twice = Some(format!("{path}: {} appears twice", c[0]));
        }
    }

    /// Every problem of the lines read.
    pub fn finish(self) -> Vec<String> {
        let path = self.path;
        if !self.ascii {
            return vec![format!("{path}: not plain ASCII")];
        }
        let mut problems = Vec::new();
        if self.carriage_return {
            problems.push(format!("{path}: carriage returns (CRLF line ends)"));
        }
        if !self.newline {
            problems.push(format!("{path}: the last line has no newline (a cut file)"));
        }
        problems.extend(self.problems);
        match self.part {
            Part::Comments => problems.push(format!(
                "{path}: the column header {} is missing after the comment lines",
                py_repr(match self.layout {
                    Layout::TwentyThree => HEADER_23ANDME,
                    Layout::Ancestry => HEADER_ANCESTRY,
                })
            )),
            Part::Rows if self.rows == 0 => problems.push(format!("{path}: no data row")),
            _ => {}
        }
        if !self.first_ok {
            problems.push(match self.layout {
                Layout::TwentyThree => {
                    format!("{path}: the first line does not start with '# This data file'")
                }
                Layout::Ancestry => {
                    format!("{path}: the first line is not {}", py_repr(ANCESTRY_FIRST))
                }
            });
        }
        problems.extend(self.out_of_order);
        problems.extend(self.twice);
        problems
    }
}

/// Every problem of `data` read as `layout`; `path` names it in each one.
pub fn check_text(path: &str, data: &[u8], layout: Layout) -> Vec<String> {
    let mut c = Checker::new(path, layout, memchr::memchr_iter(b'\n', data).count() + 1);
    data.split_inclusive(|&b| b == b'\n')
        .for_each(|l| c.line(l));
    c.finish()
}

/// The text of a `.txt` file, or of the one `.txt` member of a `.zip`;
/// `Err` holds the problems that stop the reading.
fn read_layout(path: &str) -> Result<Vec<u8>, Vec<String>> {
    let bytes = std::fs::read(path).map_err(|e| vec![format!("{path}: {e}")])?;
    if !path.ends_with(".zip") {
        return Ok(bytes);
    }
    let members =
        zip::read_all(&bytes).map_err(|e| vec![format!("{path}: not a readable zip file: {e}")])?;
    let names: Vec<&str> = members.iter().map(|(n, _)| n.as_str()).collect();
    match members.as_slice() {
        [(name, data)] if name.ends_with(".txt") => Ok(data.clone()),
        _ => Err(vec![format!(
            "{path}: want one .txt member, found {}",
            py_repr_list(&names)
        )]),
    }
}

/// `og consumer check --layout 23andme|ancestrydna FILE`.
pub fn cli(args: &[String]) -> Result<u8, Usage> {
    let p = Args::parse(args, &["layout"])?;
    let layout = match p.need("layout")? {
        "23andme" => Layout::TwentyThree,
        "ancestrydna" => Layout::Ancestry,
        other => {
            return Err(Usage(format!(
                "argument --layout: invalid choice: {} (choose from '23andme', 'ancestrydna')",
                py_repr(other)
            )));
        }
    };
    let file = match p.positional.as_slice() {
        [f] => f,
        [] => return Err(Usage("the following arguments are required: file".into())),
        [_, extra, ..] => return Err(Usage(format!("unrecognized arguments: {extra}"))),
    };
    let problems = match read_layout(file) {
        Ok(text) => check_text(file, &text, layout),
        Err(p) => p,
    };
    for p in &problems {
        println!("{p}");
    }
    if !problems.is_empty() {
        return Ok(1);
    }
    println!("{file}: keeps the {} layout", layout.name());
    Ok(0)
}

#[cfg(test)]
mod tests {
    use super::*;

    const GOOD: &str = "# This data file was generated by a test\n#\n# rsid\tchromosome\tposition\tgenotype\n\
                        rs1\t1\t100\tAG\nrs2\t1\t200\t--\nrs3\tX\t50\tA\nrs4\tMT\t10\tG\n";

    fn bad(text: &str) -> Vec<String> {
        check_text("r.txt", text.as_bytes(), Layout::TwentyThree)
    }

    #[test]
    fn a_good_file_passes() {
        assert_eq!(bad(GOOD), Vec::<String>::new());
    }

    #[test]
    fn each_planted_fault_fails() {
        let cases = [
            GOOD.replace('\n', "\r\n"),
            GOOD.replace("# rsid\tchromosome\tposition\tgenotype\n", ""),
            GOOD.replace("AG", "AN"),
            GOOD.replace("rs2\t1\t200\t--", "rs2\t1\t200\tA"),
            GOOD.replace("rs1\t1\t100", "rs1\t1\t300"),
            GOOD.replace("rs2\t", "rs1\t"),
            GOOD[..GOOD.len() - 1].to_owned(),
            GOOD.replace("rs4\tMT\t10\tG", "rs4\tMT\t10\tGG"),
        ];
        for text in &cases {
            assert!(!bad(text).is_empty(), "{text:?} passed");
        }
    }

    #[test]
    fn messages_as_the_python_wrote_them() {
        assert_eq!(
            bad(&GOOD.replace("rs1\t1\t100", "rs1\t1\t300")),
            ["r.txt: rows out of order at rs2 (after rs1)"]
        );
        assert_eq!(
            bad(&GOOD.replace("AG", "AN")),
            ["r.txt line 4: genotype is not one or two of A, C, G, T, or --: 'rs1\\t1\\t100\\tAN'"]
        );
        assert_eq!(
            bad("x\n"),
            [
                "r.txt: the column header '# rsid\\tchromosome\\tposition\\tgenotype' is missing after the comment lines",
                "r.txt: the first line does not start with '# This data file'"
            ]
        );
    }

    #[test]
    fn positions_compare_by_value() {
        assert!(bad(&GOOD.replace("rs2\t1\t200", "rs2\t1\t1000")).is_empty());
        assert!(!bad(&GOOD.replace("rs2\t1\t200", "rs2\t1\t99")).is_empty());
    }

    #[test]
    fn rsids_are_told_apart_by_kind_and_value() {
        let t = |a: &str, b: &str| {
            bad(&GOOD
                .replace("rs1\t1\t100", &format!("{a}\t1\t100"))
                .replace("rs2\t1\t200", &format!("{b}\t1\t200")))
        };
        assert!(t("rs5", "i5").is_empty());
        assert!(t("rs123456789012345678", "rs123456789012345679").is_empty());
        assert!(t("rs1234567890123456789", "rs1234567890123456780").is_empty());
        assert_eq!(t("i7", "i7"), ["r.txt: i7 appears twice"]);
        let long = "rs12345678901234567890";
        assert_eq!(t(long, long), [format!("r.txt: {long} appears twice")]);
    }

    #[test]
    fn line_by_line_is_the_whole_text() {
        // Every way a text can end, read whole and fed in pieces.
        for text in [
            GOOD,
            &GOOD[..GOOD.len() - 1],
            "",
            "\n",
            "x",
            "# This data file\n",
        ] {
            let mut c = Checker::new("r.txt", Layout::TwentyThree, 0);
            for l in text.split_inclusive('\n') {
                c.line(l.as_bytes());
            }
            assert_eq!(c.finish(), bad(text), "{text:?}");
        }
        assert_eq!(
            check_text("a.txt", ANCESTRY_FIRST.as_bytes(), Layout::Ancestry),
            [
                "a.txt: the last line has no newline (a cut file)",
                "a.txt: the column header 'rsid\\tchromosome\\tposition\\tallele1\\tallele2' is missing after the comment lines",
                "a.txt: the first line is not '#AncestryDNA raw data download'"
            ]
        );
        assert_eq!(bad("caf\u{e9}\n"), ["r.txt: not plain ASCII"]);
    }

    #[test]
    fn ancestry_codes() {
        let t = format!(
            "{ANCESTRY_FIRST}\n{HEADER_ANCESTRY}\nrs1\t25\t70000\tA\tG\nrs2\t26\t5\t0\t0\n"
        );
        assert!(check_text("a.txt", t.as_bytes(), Layout::Ancestry).is_empty());
        let t = t.replace("\t26\t", "\t27\t");
        assert!(!check_text("a.txt", t.as_bytes(), Layout::Ancestry).is_empty());
    }
}
