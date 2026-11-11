//! `og mei report` and `og mei annotsv-input` (step 50): the Alu, L1 and
//! SVA insertions MEGAnE calls, and the known disease insertions.
//!
//! ```text
//! og mei report --sample S --vcf S_megane.vcf --sites assets/mei_sites.tsv
//!               --depth S_mei_depth.tsv [--annotsv S_mei_annotsv.tsv]
//!               --out-calls S_mei.tsv --out-sites S_mei_sites.tsv
//! og mei annotsv-input --sample S --vcf S_megane.vcf --out S_mei_annotsv_input.vcf
//! ```
//!
//! MEGAnE (Kojima et al. 2023, Nat Genet, PMID 37169872) writes one record
//! per insertion: ALT `<INS:ME>`, the repeat class in INFO SVTYPE,
//! MEPRED=PASS when the inserted sequence was identified, the subfamily in
//! INFO MEI (or MEI_lbp and MEI_rbp) and the copy number in FORMAT CN.
//!
//! `--out-calls`: one row per Alu, L1 or SVA insertion with MEPRED=PASS
//! and at least one copy, sorted as a karyotype is read; with `--annotsv`
//! (AnnotSV's table of the `annotsv-input` VCF) its genes, location, benign
//! allele frequency and ACMG class. `--out-sites`: every row of the catalog
//! with its status: `seen` (an insertion of the same element inside the
//! window), `not seen` (none, with a mean depth of at least 10 over the
//! window), or `not callable`.
//!
//! `annotsv-input` writes the calls the table lists as the VCF AnnotSV
//! reads (`<INS:ME>`, `SVTYPE=INS`, `END=POS`, and the ID the table is
//! joined on).
//!
//! What changes on purpose: the depth of a catalog row is the bedcov line
//! of its own window (contig, start - 1, end), not the last line of the
//! same gene and contig, so two rows of one gene on one contig no longer
//! share a depth; the window lookup is a binary search over the sorted
//! calls; a POS, a summed depth or a catalog coordinate that is not a
//! number, and a catalog or AnnotSV row that ends early, stop with an
//! `ERROR:` line (exit 1) instead of a traceback; problems start `ERROR:`,
//! not `mei_report.py:`.

use crate::Error;
use og_filters::args::Spec;
use og_io::Utf8;
use og_tsv::{Options, Quoting};
use std::collections::HashMap;
use std::io::Write;

/// Mean depth over a catalog window before "not seen" is said.
const MIN_DEPTH: f64 = 10.0;
pub const CALL_HEADER: [&str; 12] = [
    "chrom",
    "pos",
    "end",
    "element",
    "subfamily",
    "copies",
    "filter",
    "confidence",
    "genes",
    "location",
    "benign_af",
    "acmg_class",
];
const ANNOTSV_COLUMNS: [&str; 4] = ["ID", "Annotation_mode", "Gene_name", "ACMG_class"];
pub const SITE_COLUMNS: [&str; 10] = [
    "gene",
    "chrom",
    "start",
    "end",
    "site",
    "element",
    "inheritance",
    "condition",
    "clinvar",
    "pmid",
];

/// The active element classes this report keeps.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Element {
    Alu,
    L1,
    Sva,
}

impl Element {
    fn from_svtype(s: &str) -> Option<Element> {
        match s {
            "SINE/Alu" => Some(Element::Alu),
            "LINE/L1" => Some(Element::L1),
            "Retroposon/SVA" => Some(Element::Sva),
            _ => None,
        }
    }

    fn name(self) -> &'static str {
        match self {
            Element::Alu => "Alu",
            Element::L1 => "L1",
            Element::Sva => "SVA",
        }
    }
}

/// AnnotSV's annotation of one call.
#[derive(Debug, Clone, PartialEq, Eq)]
struct Annotation {
    genes: String,
    location: Vec<String>,
    benign_af: String,
    acmg_class: String,
}

/// One kept insertion.
#[derive(Debug, Clone)]
pub struct Call {
    chrom: String,
    pos: i128,
    end: String,
    element: Element,
    subfamily: String,
    copies: String,
    filter: String,
    confidence: String,
    id: String,
    ann: Option<Annotation>,
}

fn die(msg: String) -> Error {
    Error(msg)
}

/// INFO as a dict: the last of a repeated key wins.
fn info_dict(text: &str) -> HashMap<&str, &str> {
    text.split(';')
        .filter(|x| !x.is_empty())
        .map(|x| x.split_once('=').unwrap_or((x, "")))
        .collect()
}

/// The first subfamily MEGAnE names: from MEI, else from either breakpoint.
fn subfamily(info: &HashMap<&str, &str>) -> String {
    for key in ["MEI", "MEI_lbp", "MEI_rbp"] {
        let v = info.get(key).copied().unwrap_or("");
        let name = v
            .split('|')
            .next()
            .unwrap_or("")
            .split(',')
            .next()
            .unwrap_or("");
        if !name.is_empty() && name != "NA" && name != "pA" {
            return name.into();
        }
    }
    ".".into()
}

/// The calls the table lists, and the `##contig` lines of the VCF.
pub fn read_calls(path: &str) -> Result<(Vec<Call>, Vec<String>), Error> {
    let mut calls = Vec::new();
    let mut contigs = Vec::new();
    let mut header = false;
    for (i, line) in crate::lines(path, Utf8::Strict)?.enumerate() {
        let line = &line?;
        let n = i + 1;
        if line.starts_with("##") {
            if line.starts_with("##contig=") {
                contigs.push(line.strip_suffix('\n').unwrap_or(line).to_owned());
            }
            continue;
        }
        if line.starts_with("#CHROM") {
            header = true;
            continue;
        }
        if crate::strip(line).is_empty() {
            continue;
        }
        if !header {
            return Err(die(format!(
                "{path}: a record before the #CHROM line (line {n}): not a VCF"
            )));
        }
        let f: Vec<&str> = line
            .strip_suffix('\n')
            .unwrap_or(line)
            .split('\t')
            .collect();
        if f.len() < 10 {
            return Err(die(format!(
                "{path} line {n}: {} columns, a VCF with one sample has 10",
                f.len()
            )));
        }
        let info = info_dict(f[7]);
        let (Some(svtype), Some(mepred)) = (info.get("SVTYPE"), info.get("MEPRED")) else {
            return Err(die(format!(
                "{path} line {n}: no SVTYPE or MEPRED in INFO: not MEGAnE's VCF"
            )));
        };
        let fmt: HashMap<&str, &str> = f[8].split(':').zip(f[9].split(':')).collect();
        let cn = fmt.get("CN").copied().unwrap_or(".");
        let Some(element) = Element::from_svtype(svtype) else {
            continue;
        };
        if *mepred != "PASS" || !crate::digits(cn) || crate::int(cn).is_none_or(|c| c < 1) {
            continue;
        }
        let Some(pos) = crate::int(f[1]) else {
            return Err(die(format!(
                "{path} line {n}: POS {} is not a whole number",
                og_filters::py_repr(f[1])
            )));
        };
        calls.push(Call {
            chrom: f[0].into(),
            pos,
            end: info.get("MEI_rpos").copied().unwrap_or(f[1]).into(),
            element,
            subfamily: subfamily(&info),
            copies: cn.into(),
            filter: f[6].into(),
            confidence: info.get("CN_conf").copied().unwrap_or(".").into(),
            id: String::new(),
            ann: None,
        });
    }
    if !header {
        return Err(die(format!("{path}: no #CHROM line: not a VCF")));
    }
    let mut seen: HashMap<String, usize> = HashMap::new();
    for c in &mut calls {
        let key = format!("{}_{}_{}", c.chrom, c.pos, c.element.name());
        let k = seen.entry(key.clone()).or_default();
        *k += 1;
        c.id = if *k == 1 { key } else { format!("{key}_{k}") };
    }
    Ok((calls, contigs))
}

/// The calls as the VCF AnnotSV reads: one `<INS:ME>` record each, END=POS.
fn annotsv_input(sample: &str, calls: &[Call], contigs: &[String]) -> String {
    let mut out = String::from("##fileformat=VCFv4.2\n");
    for c in contigs {
        out.push_str(c);
        out.push('\n');
    }
    out.push_str("##ALT=<ID=INS:ME,Description=\"Insertion of mobile element\">\n##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant\">\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
    out.push_str(&format!(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n"
    ));
    for c in calls {
        let gt = if c.copies == "1" { "0/1" } else { "1/1" };
        out.push_str(&format!(
            "{}\t{}\t{}\tN\t<INS:ME>\t.\t{}\tSVTYPE=INS;END={}\tGT\t{gt}\n",
            c.chrom, c.pos, c.id, c.filter, c.pos
        ));
    }
    out
}

/// AnnotSV's table of the calls, by ID. An empty file (AnnotSV had nothing
/// to annotate) is an empty map.
fn read_annotsv(path: &str) -> Result<HashMap<String, Annotation>, Error> {
    let (header, rows) = crate::read_table(path, Options::default().quoting(Quoting::None))?;
    if rows.is_empty() {
        return Ok(HashMap::new());
    }
    if ANNOTSV_COLUMNS
        .iter()
        .any(|c| !header.iter().any(|h| h == c))
    {
        return Err(die(format!(
            "{path}: not AnnotSV's table (a column of {} is missing)",
            ANNOTSV_COLUMNS.join(", ")
        )));
    }
    let mut out: HashMap<String, Annotation> = HashMap::new();
    for r in &rows {
        let id = r.text("ID")?;
        let a = out.entry(id.to_owned()).or_insert_with(|| Annotation {
            genes: ".".into(),
            location: Vec::new(),
            benign_af: ".".into(),
            acmg_class: ".".into(),
        });
        let mode = r.text("Annotation_mode")?;
        if mode == "full" {
            let genes: Vec<&str> = r
                .text("Gene_name")?
                .split(';')
                .filter(|g| !g.is_empty())
                .collect();
            a.genes = if genes.is_empty() {
                ".".into()
            } else {
                genes.join(",")
            };
            let b = crate::strip(r.get("B_ins_AFmax").unwrap_or(""));
            a.benign_af = if b.is_empty() { ".".into() } else { b.into() };
            let c = crate::strip(r.text("ACMG_class")?);
            a.acmg_class = if c.is_empty() { ".".into() } else { c.into() };
        } else if mode == "split" {
            let gene = r.text("Gene_name")?;
            let loc = crate::strip(r.get("Location").unwrap_or(""));
            if !gene.is_empty() && !loc.is_empty() {
                a.location.push(format!("{gene}:{loc}"));
            }
        }
    }
    Ok(out)
}

/// The catalog rows: each one's cells by [`SITE_COLUMNS`].
fn read_sites(path: &str) -> Result<Vec<Vec<String>>, Error> {
    let opts = Options::default().skip(|l| crate::strip(l).is_empty() || l.starts_with('#'));
    let (header, rows) = crate::read_table(path, opts)?;
    if rows.is_empty() || SITE_COLUMNS.iter().any(|c| !header.iter().any(|h| h == c)) {
        return Err(die(format!(
            "{path}: no rows, or a column of {} is missing",
            SITE_COLUMNS.join(", ")
        )));
    }
    rows.iter()
        .map(|r| {
            SITE_COLUMNS
                .iter()
                .map(|k| r.text(k).map(str::to_owned).map_err(Error::from))
                .collect()
        })
        .collect()
}

/// `samtools bedcov -c` over the catalog windows: the mean depth of each
/// window, keyed by (contig, 0-based start, end).
fn read_depth(path: &str) -> Result<HashMap<(String, i128, i128), f64>, Error> {
    let mut depth = HashMap::new();
    for (i, line) in crate::lines(path, Utf8::Strict)?.enumerate() {
        let line = &line?;
        let n = i + 1;
        if crate::strip(line).is_empty() {
            continue;
        }
        let f: Vec<&str> = line
            .strip_suffix('\n')
            .unwrap_or(line)
            .split('\t')
            .collect();
        let bad = || {
            die(format!(
                "{path} line {n}: not samtools bedcov -c output of a 4-column BED"
            ))
        };
        if f.len() < 6 || !crate::digits(f[1]) || !crate::digits(f[2]) {
            return Err(bad());
        }
        let (start, end) = (
            crate::int(f[1]).ok_or_else(bad)?,
            crate::int(f[2]).ok_or_else(bad)?,
        );
        let sum = crate::float(f[4]).ok_or_else(bad)?;
        let length = end - start;
        let mean = if length > 0 { sum / length as f64 } else { 0.0 };
        depth.insert((f[0].to_owned(), start, end), mean);
    }
    Ok(depth)
}

/// chr1 to chr22, then chrX, chrY and the rest.
fn chrom_rank(chrom: &str) -> i128 {
    let c = chrom.strip_prefix("chr").unwrap_or(chrom);
    if crate::digits(c) {
        return crate::int(c).unwrap_or(i128::MAX);
    }
    match c {
        "X" => 23,
        "Y" => 24,
        _ => 25,
    }
}

/// The status of a catalog row: (status, mean depth, call).
fn site_status(
    site: &[String],
    by_chrom: &HashMap<&str, Vec<&Call>>,
    depth: &HashMap<(String, i128, i128), f64>,
    sites_path: &str,
) -> Result<[String; 3], Error> {
    let coord = |k: usize| {
        crate::int(&site[k]).ok_or_else(|| {
            die(format!(
                "{sites_path}: {} {} is not a whole number ({})",
                SITE_COLUMNS[k],
                og_filters::py_repr(&site[k]),
                site[0]
            ))
        })
    };
    let (chrom, lo, hi) = (&site[1], coord(2)?, coord(3)?);
    let Some(&mean) = depth.get(&(chrom.clone(), lo - 1, hi)) else {
        return Ok(["not callable".into(), ".".into(), ".".into()]);
    };
    let shown = crate::fixed(mean, 1);
    let calls = by_chrom.get(chrom.as_str()).map_or(&[][..], Vec::as_slice);
    let first = calls.partition_point(|c| c.pos < lo);
    let hit = calls[first..]
        .iter()
        .take_while(|c| c.pos <= hi)
        .find(|c| c.element.name() == site[5]);
    if let Some(c) = hit {
        let cop = if c.copies == "1" { "y" } else { "ies" };
        return Ok([
            "seen".into(),
            shown,
            format!(
                "{}:{} {} {} cop{cop} {}",
                c.chrom, c.pos, c.subfamily, c.copies, c.filter
            ),
        ]);
    }
    if mean < MIN_DEPTH {
        return Ok(["not callable".into(), shown, ".".into()]);
    }
    Ok(["not seen".into(), shown, ".".into()])
}

/// The files of `og mei report`.
pub struct Files<'a> {
    pub vcf: &'a str,
    pub sites: &'a str,
    pub depth: &'a str,
    pub annotsv: Option<&'a str>,
    pub out_calls: &'a str,
    pub out_sites: &'a str,
}

/// `og mei report`: writes both tables; the stdout lines.
pub fn report(sample: &str, f: &Files<'_>) -> Result<String, Error> {
    let (mut calls, _) = read_calls(f.vcf)?;
    let ann = match f.annotsv {
        Some(p) => Some(read_annotsv(p)?),
        None => None,
    };
    for c in &mut calls {
        c.ann = ann.as_ref().and_then(|a| a.get(&c.id).cloned());
    }
    let sites = read_sites(f.sites)?;
    let depth = read_depth(f.depth)?;
    calls.sort_by(|a, b| {
        (chrom_rank(&a.chrom), &a.chrom, a.pos).cmp(&(chrom_rank(&b.chrom), &b.chrom, b.pos))
    });
    let mut body = CALL_HEADER.join("\t") + "\n";
    for c in &calls {
        let (genes, location, benign, acmg) = match &c.ann {
            Some(a) => (
                a.genes.clone(),
                if a.location.is_empty() {
                    ".".into()
                } else {
                    a.location.join(",")
                },
                a.benign_af.clone(),
                a.acmg_class.clone(),
            ),
            None => (".".into(), ".".into(), ".".into(), ".".into()),
        };
        let cells = [
            c.chrom.clone(),
            c.pos.to_string(),
            c.end.clone(),
            c.element.name().into(),
            c.subfamily.clone(),
            c.copies.clone(),
            c.filter.clone(),
            c.confidence.clone(),
            genes,
            location,
            benign,
            acmg,
        ];
        body.push_str(&cells.join("\t"));
        body.push('\n');
    }
    let mut by_chrom: HashMap<&str, Vec<&Call>> = HashMap::new();
    for c in &calls {
        by_chrom.entry(c.chrom.as_str()).or_default().push(c);
    }
    let mut sbody = SITE_COLUMNS.join("\t") + "\tstatus\tmean_depth\tcall\n";
    for s in &sites {
        let st = site_status(s, &by_chrom, &depth, f.sites)?;
        sbody.push_str(&s.join("\t"));
        for x in &st {
            sbody.push('\t');
            sbody.push_str(x);
        }
        sbody.push('\n');
    }
    crate::write_atomic(f.out_calls, body.as_bytes())?;
    crate::write_atomic(f.out_sites, sbody.as_bytes())?;
    let mut by: Vec<(Element, usize)> = Vec::new();
    for c in &calls {
        match by.iter_mut().find(|(e, _)| *e == c.element) {
            Some((_, n)) => *n += 1,
            None => by.push((c.element, 1)),
        }
    }
    by.sort_by_key(|(e, _)| e.name());
    let tally = if by.is_empty() {
        "none".to_owned()
    } else {
        by.iter()
            .map(|(e, n)| format!("{} {n}", e.name()))
            .collect::<Vec<_>>()
            .join(", ")
    };
    let mut text = format!("{sample}: {} insertions ({tally})\n", calls.len());
    if let Some(a) = &ann {
        let annotated = calls.iter().filter(|c| a.contains_key(&c.id)).count();
        let in_gene = calls
            .iter()
            .filter(|c| c.ann.as_ref().is_some_and(|x| x.genes != "."))
            .count();
        text.push_str(&format!(
            "{sample}: AnnotSV annotated {annotated} of them, {in_gene} inside a gene\n"
        ));
    }
    Ok(text)
}

const REPORT: Spec = Spec {
    prog: "og mei report",
    usage: "usage: og mei report [-h] --sample SAMPLE --vcf VCF --sites SITES --depth DEPTH [--annotsv ANNOTSV] --out-calls OUT_CALLS --out-sites OUT_SITES\n",
    options: &[
        "--sample",
        "--vcf",
        "--sites",
        "--depth",
        "--annotsv",
        "--out-calls",
        "--out-sites",
    ],
    positionals: &[],
    required: &[
        "--sample",
        "--vcf",
        "--sites",
        "--depth",
        "--out-calls",
        "--out-sites",
    ],
};

const INPUT: Spec = Spec {
    prog: "og mei annotsv-input",
    usage: "usage: og mei annotsv-input [-h] --sample SAMPLE --vcf VCF --out OUT\n",
    options: &["--sample", "--vcf", "--out"],
    positionals: &[],
    required: &["--sample", "--vcf", "--out"],
};

/// `og mei VERB ARGS`: the exit code (1 for a problem, 2 for a refused
/// command line).
pub fn cli(verb: Option<&str>, rest: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    match verb {
        Some("report") => crate::run(&REPORT, rest, out, err, |a| {
            let files = Files {
                vcf: crate::opt(a, "--vcf"),
                sites: crate::opt(a, "--sites"),
                depth: crate::opt(a, "--depth"),
                annotsv: a.get("--annotsv"),
                out_calls: crate::opt(a, "--out-calls"),
                out_sites: crate::opt(a, "--out-sites"),
            };
            Ok(report(crate::opt(a, "--sample"), &files).map_err(|e| (e, 1)))
        }),
        Some("annotsv-input") => crate::run(&INPUT, rest, out, err, |a| {
            let sample = crate::opt(a, "--sample");
            Ok((|| {
                let (calls, contigs) = read_calls(crate::opt(a, "--vcf"))?;
                let body = annotsv_input(sample, &calls, &contigs);
                crate::write_atomic(crate::opt(a, "--out"), body.as_bytes())?;
                Ok(format!(
                    "{sample}: {} insertions written for AnnotSV\n",
                    calls.len()
                ))
            })()
            .map_err(|e| (e, 1)))
        }),
        other => crate::unknown_verb("mei", "report|annotsv-input", other, err),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testutil::Dir;

    const HEAD: &str = "##fileformat=VCFv4.1\n##contig=<ID=chr4,length=190214555>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n";

    fn rec(chrom: &str, pos: u64, svtype: &str, cn: &str) -> String {
        format!(
            "{chrom}\t{pos}\t.\tN\t<INS:ME>\t.\tPASS\tSVTYPE={svtype};MEPRED=PASS;MEI=SVA_E,x|y;CN_conf=H\tGT:CN\t0/1:{cn}\n"
        )
    }

    const SITES: &str = "# catalog\ngene\tchrom\tstart\tend\tsite\telement\tinheritance\tcondition\tclinvar\tpmid\nFKTN\tchr4\t105620000\t105620600\t3' UTR\tSVA\tAR\tFukuyama\t.\t1\nFKTN\tchr4\t105630000\t105630600\tother\tSVA\tAR\tmade up\t.\t1\n";
    const DEPTH: &str = "chr4\t105619999\t105620600\tFKTN\t18030\t500\nchr4\t105629999\t105630600\tFKTN\t3005\t20\n";

    fn run(vcf: &str) -> (String, String) {
        let d = Dir::new("mei");
        let f = Files {
            vcf: &d.put("m.vcf", vcf),
            sites: &d.put("sites.tsv", SITES),
            depth: &d.put("depth.tsv", DEPTH),
            annotsv: None,
            out_calls: &d.path("calls.tsv"),
            out_sites: &d.path("sites_out.tsv"),
        };
        report("S", &f).unwrap();
        (d.read("calls.tsv"), d.read("sites_out.tsv"))
    }

    #[test]
    fn a_planted_sva_is_seen_and_each_window_has_its_own_depth() {
        let (calls, sites) = run(&format!(
            "{HEAD}{}{}",
            rec("chr4", 105620300, "Retroposon/SVA", "1"),
            rec("chr4", 105620300, "Retroposon/SVA", "2")
        ));
        assert_eq!(calls.lines().count(), 3);
        let rows: Vec<&str> = sites.lines().collect();
        assert!(
            rows[1].ends_with("\tseen\t30.0\tchr4:105620300 SVA_E 1 copy PASS"),
            "{}",
            rows[1]
        );
        // The second FKTN window has its own bedcov line: 5x, too few reads.
        assert!(rows[2].ends_with("\tnot callable\t5.0\t."), "{}", rows[2]);
    }

    /// Red first: a mean depth of exactly 10 is enough to say "not seen".
    #[test]
    fn red_first_a_depth_of_ten_is_callable() {
        let d = Dir::new("mei-ten");
        let depth: HashMap<_, _> = [(("chr4".to_owned(), 99, 200), 10.0)].into();
        let site: Vec<String> = ["G", "chr4", "100", "200", "s", "Alu", "AD", "c", ".", "1"]
            .map(String::from)
            .to_vec();
        let st = site_status(&site, &HashMap::new(), &depth, &d.path("x")).unwrap();
        assert_eq!(st[0], "not seen");
    }

    #[test]
    fn ids_and_annotsv_input() {
        let d = Dir::new("mei-ids");
        let p = d.put(
            "m.vcf",
            &format!(
                "{HEAD}{}{}{}",
                rec("chr4", 5, "SINE/Alu", "1"),
                rec("chr4", 5, "SINE/Alu", "2"),
                rec("chr4", 6, "DNA/hAT", "1")
            ),
        );
        let (calls, contigs) = read_calls(&p).unwrap();
        assert_eq!(
            calls.iter().map(|c| c.id.as_str()).collect::<Vec<_>>(),
            ["chr4_5_Alu", "chr4_5_Alu_2"]
        );
        let v = annotsv_input("S", &calls, &contigs);
        assert!(v.contains("##contig=<ID=chr4,length=190214555>\n"));
        assert!(
            v.ends_with("chr4\t5\tchr4_5_Alu_2\tN\t<INS:ME>\t.\tPASS\tSVTYPE=INS;END=5\tGT\t1/1\n")
        );
        assert_eq!(chrom_rank("chrX"), 23);
        assert_eq!(chrom_rank("chr10"), 10);
        assert_eq!(chrom_rank("chrUn_x"), 25);
    }
}
