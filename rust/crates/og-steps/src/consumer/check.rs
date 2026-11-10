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

/// A whole number written without leading zeros, compared by value.
fn num_key(s: &str) -> (usize, &str) {
    (s.len(), s)
}

/// The sort key of a row: chromosome order, then position.
type Key<'a> = (usize, (usize, &'a str));

/// The order and rsID checks, row by row: the first row out of order, and
/// the first rsID seen twice.
#[derive(Default)]
struct Order<'a> {
    last: Option<(Key<'a>, &'a str)>,
    out_of_order: Option<(&'a str, &'a str)>,
    seen: FxHashSet<&'a str>,
    twice: Option<&'a str>,
}

impl<'a> Order<'a> {
    fn row(&mut self, key: Key<'a>, id: &'a str) {
        if self.out_of_order.is_none()
            && let Some((last, last_id)) = self.last
            && key < last
        {
            self.out_of_order = Some((id, last_id));
        }
        self.last = Some((key, id));
        if self.twice.is_none() && !self.seen.insert(id) {
            self.twice = Some(id);
        }
    }

    fn problems(self, path: &str, problems: &mut Vec<String>) {
        if let Some((b, a)) = self.out_of_order {
            problems.push(format!("{path}: rows out of order at {b} (after {a})"));
        }
        if let Some(id) = self.twice {
            problems.push(format!("{path}: {id} appears twice"));
        }
    }
}

/// Every problem of `data` read as `layout`; `path` names it in each one.
pub fn check_text(path: &str, data: &[u8], layout: Layout) -> Vec<String> {
    let Ok(text) = std::str::from_utf8(data) else {
        return vec![format!("{path}: not plain ASCII")];
    };
    if !text.is_ascii() {
        return vec![format!("{path}: not plain ASCII")];
    }
    let (header, ncols) = match layout {
        Layout::TwentyThree => (HEADER_23ANDME, 4),
        Layout::Ancestry => (HEADER_ANCESTRY, 5),
    };
    let mut problems = Vec::new();
    if text.contains('\r') {
        problems.push(format!("{path}: carriage returns (CRLF line ends)"));
    }
    if !text.ends_with('\n') {
        problems.push(format!("{path}: the last line has no newline (a cut file)"));
    }
    let body = text.strip_suffix('\n').unwrap_or(text);
    // `text.split("\n")[:-1]` of an empty text is no line at all.
    let mut lines = (!text.is_empty())
        .then(|| body.split('\n'))
        .into_iter()
        .flatten()
        .peekable();
    let comment = |l: &str| match layout {
        Layout::TwentyThree => l.starts_with('#') && l != HEADER_23ANDME,
        Layout::Ancestry => l.starts_with('#'),
    };
    let mut n = 0usize;
    while lines.next_if(|l| comment(l)).is_some() {
        n += 1;
    }
    let mut order = Order::default();
    // One rsID per line at most: room for all of them up front, so the set
    // never grows (and never holds two tables while it does).
    order
        .seen
        .reserve(memchr::memchr_iter(b'\n', data).count() + 1);
    if lines.next() != Some(header) {
        problems.push(format!(
            "{path}: the column header {} is missing after the comment lines",
            py_repr(header)
        ));
    } else {
        n += 1;
        let mut rows = 0usize;
        let mut c: Vec<&str> = Vec::with_capacity(ncols + 1);
        for line in lines {
            n += 1;
            c.clear();
            c.extend(line.split('\t'));
            if c.len() != ncols {
                problems.push(format!(
                    "{path} line {n}: {} columns, want {ncols}: {}",
                    c.len(),
                    py_repr(head(line, 60))
                ));
                continue;
            }
            let err = match layout {
                Layout::TwentyThree => row_23andme(&c),
                Layout::Ancestry => row_ancestry(&c),
            };
            if let Some(e) = err {
                problems.push(format!("{path} line {n}: {e}: {}", py_repr(head(line, 60))));
                continue;
            }
            rows += 1;
            let chrom = match layout {
                Layout::TwentyThree => usize::from(chrom_index(c[1]).unwrap_or(0)),
                Layout::Ancestry => c[1].parse().unwrap_or(0),
            };
            order.row((chrom, num_key(c[2])), c[0]);
        }
        if rows == 0 {
            problems.push(format!("{path}: no data row"));
        }
    }
    match layout {
        Layout::TwentyThree => {
            if !text.starts_with("# This data file") {
                problems.push(format!(
                    "{path}: the first line does not start with '# This data file'"
                ));
            }
        }
        Layout::Ancestry => {
            if !text.starts_with(&format!("{ANCESTRY_FIRST}\n")) {
                problems.push(format!(
                    "{path}: the first line is not {}",
                    py_repr(ANCESTRY_FIRST)
                ));
            }
        }
    }
    order.problems(path, &mut problems);
    problems
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
    fn ancestry_codes() {
        let t = format!(
            "{ANCESTRY_FIRST}\n{HEADER_ANCESTRY}\nrs1\t25\t70000\tA\tG\nrs2\t26\t5\t0\t0\n"
        );
        assert!(check_text("a.txt", t.as_bytes(), Layout::Ancestry).is_empty());
        let t = t.replace("\t26\t", "\t27\t");
        assert!(!check_text("a.txt", t.as_bytes(), Layout::Ancestry).is_empty());
    }
}
