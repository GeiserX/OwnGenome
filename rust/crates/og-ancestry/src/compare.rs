//! `og archaic compare` and `og archaic genomes`: step 43's Neanderthal
//! totals beside the calls Chen et al. 2020 released for the 1000 Genomes
//! people. Only the archaic jobs of `public-evidence.yml` run them; they are
//! not part of the pipeline.
//!
//! The released calls are Chen et al.'s file "Neanderthal sequence in 1000
//! genome.50kb.txt" (chr, start, end, LOD, maxLOD, size, pop, anc, ID,
//! whitespace separated), merged per person and chromosome.
//!
//! - `compare`: the same person run through step 43 on one or more inputs
//!   (a segments TSV each, `--mine LABEL=FILE`), beside the released figure
//!   for that person and for the people of one group, and how much of each
//!   lies under the released calls.
//! - `genomes`: the raw IBDmix output of one step 43 run, per archaic genome
//!   (Altai, Vindija, Denisova) at LOD 4 and 50 kb, and the Altai and
//!   Vindija union, each beside the released calls.
//!
//! Both print a markdown table and a few lines of text on stdout. Lengths
//! are base pairs of the union of each person's intervals on a chromosome
//! ([`og_intervals::sets::Set`]); a median of an even count is the mean of
//! the middle two, printed to the nearest whole number with ties to even.
//!
//! They print what `bin/archaic_compare.py` printed, byte for byte (checked
//! by `tests/parity`). Deliberate changes: a malformed number, a missing
//! file or an empty spans file is an `ERROR:` line and exit 1, never a
//! traceback; an interval with a negative start or an end before its start
//! stops the command, where the Python summed it as a negative length; a
//! `--mine` value without `=` is a usage error (exit 2), where the Python
//! tried to open a file with no name; and the chromosome of a spans row
//! with no tab is read without its line ending.

use super::{ArgError, Args, Fail, Result, Spec, for_each_line, int, split_ws};
use og_intervals::sets::Set;
use pycompat::fmt::fixed;
use std::collections::HashMap;
use std::io::Write;
use std::path::Path;

/// The archaic genomes IBDmix compares with, in the order `genomes` lists them.
const ARCHAICS: [&str; 3] = ["Altai", "Vindija", "Denisova"];
/// Chen et al. 2020's cutoffs, which step 43 keeps too: LOD 4 and 50 kb.
const MIN_LOD: f64 = 4.0;
const MIN_LENGTH: i64 = 50_000;
const SEGMENT_HEADER: [&str; 5] = ["chrom", "start", "end", "length_bp", "ancestry"];
const IBDMIX_HEADER: [&str; 5] = ["ID", "chrom", "start", "end", "slod"];

/// [`for_each_line`], with a file that cannot be read an `ERROR:` line
/// like the other messages the Python never printed.
fn lines(path: &Path, f: impl FnMut(&str) -> Result<()>) -> Result<()> {
    for_each_line(path, f).map_err(|Fail(m)| {
        Fail(if m.starts_with("cannot read ") {
            format!("ERROR: {m}")
        } else {
            m
        })
    })
}

/// Intervals per chromosome, as read (not merged).
type Rows = HashMap<String, Vec<(u64, u64)>>;

/// The chromosome name without a leading `chr`.
fn norm(c: &str) -> &str {
    c.strip_prefix("chr").unwrap_or(c)
}

/// How `--chroms` sorts: numbered chromosomes first, in numeric order, then
/// the rest by name.
fn chrom_order(c: &str) -> (u8, i128, &str) {
    if pycompat::text::isdigit(c)
        && let Ok(n) = pycompat::num::parse_int(c)
    {
        return (0, n, c);
    }
    (1, 0, c)
}

/// One interval of `path`: start and end as whole numbers, the start not
/// negative and the end not before it.
fn interval(path: &str, s: &str, e: &str) -> Result<(u64, u64)> {
    let num =
        |v: &str| int(v).ok_or_else(|| Fail(format!("ERROR: {path}: {v:?} is not a whole number")));
    let (s, e) = (num(s)?, num(e)?);
    if s < 0 || e < s {
        return Err(Fail(format!(
            "ERROR: {path}: the interval {s}-{e} has a negative start or ends before it starts"
        )));
    }
    Ok((s as u64, e as u64))
}

/// Chen et al.'s released calls: the intervals of each person per
/// chromosome, and each person's group (the `anc` column, the last row's
/// when they differ).
struct Chen {
    calls: HashMap<String, Rows>,
    group: HashMap<String, String>,
}

impl Chen {
    fn read(path: &str) -> Result<Chen> {
        let mut chen = Chen {
            calls: HashMap::new(),
            group: HashMap::new(),
        };
        lines(Path::new(path), |line| {
            let c: Vec<&str> = split_ws(line).collect();
            if c.len() < 9 || c[0] == "chr" {
                return Ok(());
            }
            let iv = interval(path, c[1], c[2])?;
            chen.calls
                .entry(c[8].to_owned())
                .or_default()
                .entry(norm(c[0]).to_owned())
                .or_default()
                .push(iv);
            chen.group.insert(c[8].to_owned(), c[7].to_owned());
            Ok(())
        })?;
        Ok(chen)
    }

    /// One person's intervals on one chromosome, merged.
    fn set(&self, person: &str, chrom: &str) -> Set {
        rows_set(self.calls.get(person), chrom)
    }

    /// One person's merged base pairs on one chromosome.
    fn bp(&self, person: &str, chrom: &str) -> u64 {
        self.set(person, chrom).bases()
    }
}

/// The merged intervals of `chrom` in `rows`; empty when there are none.
fn rows_set(rows: Option<&Rows>, chrom: &str) -> Set {
    rows.and_then(|r| r.get(chrom))
        .map(|v| Set::new(v.iter().copied()))
        .unwrap_or_default()
}

/// The Neanderthal rows of a step 43 segments TSV, per chromosome.
fn read_mine(path: &str, into: &mut Rows) -> Result<()> {
    let mut first = true;
    lines(Path::new(path), |line| {
        let c: Vec<&str> = line.split('\t').collect();
        if std::mem::take(&mut first) {
            if c.len() < 5 || c[..5] != SEGMENT_HEADER {
                return Err(Fail(format!("{path}: not a step 43 segments file")));
            }
            return Ok(());
        }
        if c.len() < 5 {
            return Err(Fail(format!(
                "ERROR: {path}: a row has {} columns, want at least 5",
                c.len()
            )));
        }
        if c[4] == "Neanderthal" {
            let iv = interval(path, c[1], c[2])?;
            into.entry(norm(c[0]).to_owned()).or_default().push(iv);
        }
        Ok(())
    })?;
    if first {
        return Err(Fail(format!("{path}: not a step 43 segments file")));
    }
    Ok(())
}

/// IBDmix's segments for SAMPLE in one of its output files, kept as step 43
/// keeps them (`og archaic report`: LOD 4 and 50 kb or longer).
fn ibdmix_segments(path: &str) -> Result<Vec<(u64, u64)>> {
    let mut out = Vec::new();
    let mut first = true;
    lines(Path::new(path), |line| {
        let c: Vec<&str> = line.split('\t').collect();
        if std::mem::take(&mut first) {
            let head = &c[..c.len().min(5)];
            if head != IBDMIX_HEADER {
                return Err(Fail(format!(
                    "ERROR: {path} is not IBDmix output (header {})",
                    head.join(" ")
                )));
            }
            return Ok(());
        }
        if c.len() < 5 || c[0] != "SAMPLE" {
            return Ok(());
        }
        let (s, e) = (whole(path, c[2])?, whole(path, c[3])?);
        let lod = pycompat::num::parse_float(c[4])
            .map_err(|_| Fail(format!("ERROR: {path}: {:?} is not a number", c[4])))?;
        if lod >= MIN_LOD && e - s >= MIN_LENGTH {
            out.push(interval(path, c[2], c[3])?);
        }
        Ok(())
    })?;
    if first {
        return Err(Fail(format!(
            "ERROR: {path} is not IBDmix output (header )"
        )));
    }
    Ok(out)
}

/// A whole number of `path`.
fn whole(path: &str, v: &str) -> Result<i64> {
    int(v).ok_or_else(|| Fail(format!("ERROR: {path}: {v:?} is not a whole number")))
}

/// `statistics.median` of sorted counts, as `f"{m:.0f}"` prints it.
fn median(sorted: &[u64]) -> String {
    let n = sorted.len();
    if n % 2 == 1 {
        sorted[n / 2].to_string()
    } else {
        fixed((sorted[n / 2 - 1] + sorted[n / 2]) as f64 / 2.0, 0)
    }
}

/// The table and the lines of `og archaic compare`.
fn cmd_compare(a: &Args, out: &mut String) -> Result<()> {
    let chen_path = a.req("--chen-file");
    let sample = a.req("--sample");
    let group = a.get("--group").unwrap_or("EUR");
    let chen = Chen::read(chen_path)?;
    if !chen.calls.contains_key(sample) {
        return Err(Fail(format!("{sample} is not in {chen_path}")));
    }
    let peers: Vec<&str> = chen
        .group
        .iter()
        .filter(|(_, g)| *g == group)
        .map(|(p, _)| p.as_str())
        .collect();
    if peers.is_empty() {
        return Err(Fail(format!("no person of group {group} in {chen_path}")));
    }
    let mut labels: Vec<&str> = Vec::new();
    let mut ivs: HashMap<&str, Rows> = HashMap::new();
    for item in a.all("--mine") {
        let (label, path) = item.split_once('=').expect("checked by compare_cli");
        if !labels.contains(&label) {
            labels.push(label);
        }
        read_mine(path, ivs.entry(label).or_default())?;
    }
    let chroms_arg = a.req("--chroms").replace(',', " ");
    let mut chroms: Vec<&str> = Vec::new();
    for c in split_ws(&chroms_arg).map(norm) {
        if !chroms.contains(&c) {
            chroms.push(c);
        }
    }
    chroms.sort_by(|x, y| chrom_order(x).cmp(&chrom_order(y)));
    if chroms.is_empty() {
        return Err(Fail("--chroms names no chromosome".to_owned()));
    }

    out.push_str(&format!(
        "| chromosome | Chen {sample} (bp) | Chen {group} median | Chen {group} max | {} |\n",
        labels
            .iter()
            .map(|l| format!("{l} (bp)"))
            .collect::<Vec<_>>()
            .join(" | ")
    ));
    out.push_str(&format!(
        "|---|---|---|---|{}\n",
        "---|".repeat(labels.len())
    ));
    // The "all" median and maximum come from person-level totals, not from
    // sums of per-chromosome medians or maxima, which no person need have.
    let mut sums: Vec<u64> = peers
        .iter()
        .map(|p| chroms.iter().map(|c| chen.bp(p, c)).sum())
        .collect();
    sums.sort_unstable();
    let mut tot_me = 0;
    let mut tot_mine = vec![0u64; labels.len()];
    for c in &chroms {
        let mut vals: Vec<u64> = peers.iter().map(|p| chen.bp(p, c)).collect();
        vals.sort_unstable();
        let me = chen.bp(sample, c);
        tot_me += me;
        let cells: Vec<String> = labels
            .iter()
            .zip(tot_mine.iter_mut())
            .map(|(l, tot)| {
                let v = rows_set(ivs.get(l), c).bases();
                *tot += v;
                v.to_string()
            })
            .collect();
        out.push_str(&format!(
            "| {c} | {me} | {} | {} | {} |\n",
            median(&vals),
            vals[vals.len() - 1],
            cells.join(" | ")
        ));
    }
    out.push_str(&format!(
        "| all | {tot_me} | {} | {} | {} |\n\n",
        median(&sums),
        sums[sums.len() - 1],
        tot_mine
            .iter()
            .map(u64::to_string)
            .collect::<Vec<_>>()
            .join(" | ")
    ));
    // The person-level spread over the same chromosomes.
    out.push_str(&format!(
        "Chen et al., {} {group} people, these {} chromosomes: lowest {}, median {}, highest {} bp; {sample} {tot_me} bp.\n",
        peers.len(),
        chroms.len(),
        sums[0],
        median(&sums),
        sums[sums.len() - 1],
    ));
    for (l, &tot) in labels.iter().zip(&tot_mine) {
        let ratio = if tot_me == 0 {
            f64::NAN
        } else {
            tot as f64 / tot_me as f64
        };
        let both: u64 = chroms
            .iter()
            .map(|c| {
                rows_set(ivs.get(l), c)
                    .intersect(&chen.set(sample, c))
                    .bases()
            })
            .sum();
        out.push_str(&format!(
            "{l}: {tot} bp, {} times the released {sample} figure; {both} bp under the released calls, {} bp only here, {} bp only released.\n",
            fixed(ratio, 2),
            tot - both,
            tot_me - both,
        ));
    }
    Ok(())
}

/// The table of `og archaic genomes`.
fn cmd_genomes(a: &Args, out: &mut String) -> Result<()> {
    let chen_path = a.req("--chen-file");
    let sample = a.req("--sample");
    let chen = Chen::read(chen_path)?;
    if !chen.calls.contains_key(sample) {
        return Err(Fail(format!("{sample} is not in {chen_path}")));
    }
    let spans = a.req("--spans");
    let mut chroms: Vec<String> = Vec::new();
    let mut header = true;
    lines(Path::new(spans), |line| {
        if std::mem::take(&mut header) || split_ws(line).next().is_none() {
            return Ok(());
        }
        let c = norm(line.split('\t').next().unwrap_or(""));
        if !chroms.iter().any(|x| x == c) {
            chroms.push(c.to_owned());
        }
        Ok(())
    })?;
    if header {
        return Err(Fail(format!("ERROR: {spans} is empty")));
    }
    let dir = Path::new(a.req("--ibdmix-dir"));
    // Per archaic genome, then the Altai and Vindija union: per chromosome.
    let mut per: Vec<Vec<Set>> = vec![Vec::new(); 4];
    for c in &chroms {
        let mut segs: Vec<Vec<(u64, u64)>> = Vec::new();
        for n in ARCHAICS {
            let p = dir.join(format!("chr{c}.{n}.txt"));
            segs.push(ibdmix_segments(&p.to_string_lossy())?);
        }
        let union: Vec<(u64, u64)> = segs[0].iter().chain(&segs[1]).copied().collect();
        segs.push(union);
        for (k, s) in segs.into_iter().enumerate() {
            per[k].push(Set::new(s));
        }
    }
    let released_sets: Vec<Set> = chroms.iter().map(|c| chen.set(sample, c)).collect();
    let released: u64 = released_sets.iter().map(Set::bases).sum();
    out.push_str(&format!(
        "| IBDmix output | bp | under the released {sample} calls | only here | only released |\n|---|---|---|---|---|\n"
    ));
    for (name, sets) in ["Altai", "Vindija", "Denisova", "Altai+Vindija"]
        .iter()
        .zip(&per)
    {
        let mine: u64 = sets.iter().map(Set::bases).sum();
        let both: u64 = sets
            .iter()
            .zip(&released_sets)
            .map(|(m, r)| m.intersect(r).bases())
            .sum();
        out.push_str(&format!(
            "| {name} | {mine} | {both} | {} | {} |\n",
            mine - both,
            released - both
        ));
    }
    out.push_str(&format!(
        "\nReleased {sample}: {released} bp on these {} chromosomes.\n",
        chroms.len()
    ));
    Ok(())
}

const COMPARE: Spec = Spec {
    verb: "compare",
    usage: "usage: og archaic compare --chen-file CHEN_FILE --sample SAMPLE [--group GROUP] --chroms CHROMS --mine LABEL=SEGMENTS.tsv [--mine LABEL=SEGMENTS.tsv ...]\n",
    options: &["--chen-file", "--sample", "--group", "--chroms", "--mine"],
    required: &["--chen-file", "--sample", "--chroms", "--mine"],
};

const GENOMES: Spec = Spec {
    verb: "genomes",
    usage: "usage: og archaic genomes --chen-file CHEN_FILE --sample SAMPLE --ibdmix-dir IBDMIX_DIR --spans SPANS\n",
    options: &["--chen-file", "--sample", "--ibdmix-dir", "--spans"],
    required: &["--chen-file", "--sample", "--ibdmix-dir", "--spans"],
};

/// Parses `args` for `spec`, runs `body` and prints what it wrote: the exit
/// code (0, 1 when it stopped, 2 for a wrong command line).
fn run(
    spec: &Spec,
    args: &[String],
    body: fn(&Args, &mut String) -> Result<()>,
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> u8 {
    let a = match Args::parse(spec, args).and_then(|a| {
        if spec.verb == "compare"
            && let Some(m) = a.all("--mine").into_iter().find(|m| !m.contains('='))
        {
            return Err(ArgError::Usage(format!(
                "argument --mine: {m:?} is not LABEL=SEGMENTS.tsv"
            )));
        }
        Ok(a)
    }) {
        Ok(a) => a,
        Err(ArgError::Help) => {
            let _ = write!(out, "{}", spec.usage);
            return 0;
        }
        Err(ArgError::Usage(msg)) => {
            let _ = write!(err, "{}", spec.usage);
            let _ = writeln!(err, "og archaic {}: error: {msg}", spec.verb);
            return 2;
        }
    };
    let mut text = String::new();
    match body(&a, &mut text) {
        Ok(()) => match out.write_all(text.as_bytes()).and_then(|()| out.flush()) {
            Ok(()) => 0,
            Err(e) => {
                let _ = writeln!(err, "ERROR: cannot write the output: {e}");
                1
            }
        },
        Err(Fail(msg)) => {
            let _ = writeln!(err, "{msg}");
            1
        }
    }
}

/// `og archaic compare ARGS`: the exit code. The table goes to `out`,
/// messages to `err`.
pub fn compare_cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    run(&COMPARE, args, cmd_compare, out, err)
}

/// `og archaic genomes ARGS`: the exit code. The table goes to `out`,
/// messages to `err`.
pub fn genomes_cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    run(&GENOMES, args, cmd_genomes, out, err)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::path::PathBuf;

    const CHEN: &str = "chr start   end LOD maxLOD  size    pop anc ID   \n\
        22\t100\t200\t9\t9\t100\tCEU\tEUR\tA\n\
        22\t150\t300\t9\t9\t150\tCEU\tEUR\tA\n\
        22\t1000\t1100\t9\t9\t100\tCEU\tEUR\tA\n\
        21\t10\t60\t9\t9\t50\tCEU\tEUR\tA\n\
        22\t100\t150\t9\t9\t50\tGBR\tEUR\tB\n\
        22\t100\t900\t9\t9\t800\tYRI\tAFR\tC\n";
    const MINE: &str = "chrom\tstart\tend\tlength_bp\tancestry\tmatched\tlod\tgenes\n\
        chr22\t100\t400\t300\tNeanderthal\tAltai\t5\t\n\
        chr22\t5000\t6000\t1000\tDenisovan\tDenisova\t5\t\n";

    /// A scratch folder of its own for each test.
    fn dir(name: &str) -> PathBuf {
        let d =
            std::env::temp_dir().join(format!("og-archaic-compare-{}-{name}", std::process::id()));
        let _ = fs::remove_dir_all(&d);
        fs::create_dir_all(&d).unwrap();
        d
    }

    fn put(d: &Path, name: &str, text: &str) -> String {
        let p = d.join(name);
        fs::write(&p, text).unwrap();
        p.to_string_lossy().into_owned()
    }

    fn strings(args: &[&str]) -> Vec<String> {
        args.iter().map(|s| (*s).to_owned()).collect()
    }

    /// Runs a verb: the exit code, stdout and stderr.
    fn call(
        f: fn(&[String], &mut dyn Write, &mut dyn Write) -> u8,
        args: &[&str],
    ) -> (u8, String, String) {
        let (mut out, mut err) = (Vec::new(), Vec::new());
        let code = f(&strings(args), &mut out, &mut err);
        (
            code,
            String::from_utf8(out).unwrap(),
            String::from_utf8(err).unwrap(),
        )
    }

    /// `compare` with sample A, group EUR, chromosomes 21 and 22 and the one
    /// segments file as `x`; `extra` comes last, so it can override.
    fn run_cmp(d: &Path, extra: &[&str]) -> (u8, String, String) {
        let chen = d.join("chen.txt").to_string_lossy().into_owned();
        let mine = format!("x={}", d.join("mine.tsv").display());
        let mut args = vec![
            "--chen-file",
            &chen,
            "--sample",
            "A",
            "--group",
            "EUR",
            "--chroms",
            "21 22",
            "--mine",
            &mine,
        ];
        args.extend_from_slice(extra);
        call(compare_cli, &args)
    }

    fn setup(name: &str) -> PathBuf {
        let d = dir(name);
        put(&d, "chen.txt", CHEN);
        put(&d, "mine.tsv", MINE);
        d
    }

    #[test]
    fn merged_bp_is_the_union() {
        let s = Set::new([(100, 200), (150, 300), (1000, 1100)]);
        assert_eq!(s.bases(), 300);
        assert_eq!(Set::new([]).bases(), 0);
    }

    #[test]
    fn intersect_bp_counts_shared_bases() {
        let a = Set::new([(0, 10), (20, 30)]);
        assert_eq!(a.intersect(&Set::new([(5, 25)])).bases(), 10);
        assert_eq!(Set::new([(0, 10)]).intersect(&Set::new([])).bases(), 0);
    }

    #[test]
    fn person_total_is_the_union() {
        let d = setup("union");
        let (code, out, err) = run_cmp(&d, &[]);
        assert_eq!(code, 0, "{err}");
        // chr22: 100-300 and 1000-1100 = 300; chr21: 50
        assert!(out.contains("| 22 | 300 |"), "{out}");
        assert!(out.contains("| 21 | 50 |"), "{out}");
        assert!(out.contains("| all | 350 |"), "{out}");
    }

    #[test]
    fn only_neanderthal_rows_of_mine_count() {
        let d = setup("neanderthal");
        let (_, out, _) = run_cmp(&d, &[]);
        assert!(
            out.contains("x: 300 bp, 0.86 times the released A figure;"),
            "{out}"
        );
    }

    #[test]
    fn overlap_with_the_released_calls() {
        let d = setup("overlap");
        let (_, out, _) = run_cmp(&d, &[]);
        // mine chr22 100-400; released chr22 100-300 and 1000-1100: 200 bp shared
        assert!(
            out.contains("200 bp under the released calls, 100 bp only here, 150 bp only released"),
            "{out}"
        );
    }

    #[test]
    fn whole_output() {
        let d = setup("whole");
        let (code, out, err) = run_cmp(&d, &[]);
        assert_eq!((code, err.as_str()), (0, ""));
        assert_eq!(
            out,
            "| chromosome | Chen A (bp) | Chen EUR median | Chen EUR max | x (bp) |\n\
             |---|---|---|---|---|\n\
             | 21 | 50 | 25 | 50 | 0 |\n\
             | 22 | 300 | 175 | 300 | 300 |\n\
             | all | 350 | 200 | 350 | 300 |\n\
             \n\
             Chen et al., 2 EUR people, these 2 chromosomes: lowest 50, median 200, highest 350 bp; A 350 bp.\n\
             x: 300 bp, 0.86 times the released A figure; 200 bp under the released calls, 100 bp only here, 150 bp only released.\n"
        );
    }

    #[test]
    fn all_row_maximum_is_a_person_total() {
        // A is largest on chr21 (100), B on chr22 (100): summed maxima would
        // be 200, no person has it
        let d = setup("max");
        put(
            &d,
            "chen.txt",
            "chr start end LOD maxLOD size pop anc ID\n\
             21\t0\t100\t9\t9\t100\tCEU\tEUR\tA\n22\t0\t10\t9\t9\t10\tCEU\tEUR\tA\n\
             21\t0\t10\t9\t9\t10\tGBR\tEUR\tB\n22\t0\t100\t9\t9\t100\tGBR\tEUR\tB\n",
        );
        let (_, out, _) = run_cmp(&d, &[]);
        assert!(out.contains("| all | 110 | 110 | 110 |"), "{out}");
    }

    #[test]
    fn all_row_median_is_a_person_total() {
        // each of A, B, D has 100 bp on two of three chromosomes: every person
        // totals 200, while summed per-chromosome medians would give 300
        let d = setup("median");
        put(
            &d,
            "chen.txt",
            "chr start end LOD maxLOD size pop anc ID\n\
             20\t0\t100\t9\t9\t100\tCEU\tEUR\tA\n21\t0\t100\t9\t9\t100\tCEU\tEUR\tA\n\
             21\t0\t100\t9\t9\t100\tGBR\tEUR\tB\n22\t0\t100\t9\t9\t100\tGBR\tEUR\tB\n\
             20\t0\t100\t9\t9\t100\tTSI\tEUR\tD\n22\t0\t100\t9\t9\t100\tTSI\tEUR\tD\n",
        );
        let (_, out, _) = run_cmp(&d, &["--chroms", "20 21 22"]);
        assert!(out.contains("| all | 200 | 200 | 200 |"), "{out}");
    }

    #[test]
    fn group_spread_uses_the_group_only() {
        let d = setup("group");
        let (_, out, _) = run_cmp(&d, &[]);
        assert!(out.contains("2 EUR people"), "{out}"); // A and B, not the AFR person C
        assert!(out.contains("highest 350"), "{out}");
    }

    #[test]
    fn median_of_an_even_count_rounds_half_to_even() {
        assert_eq!(median(&[1, 2]), "2"); // 1.5
        assert_eq!(median(&[2, 3]), "2"); // 2.5
        assert_eq!(median(&[3, 4]), "4"); // 3.5
        assert_eq!(median(&[7]), "7");
        assert_eq!(median(&[1, 5, 9]), "5");
    }

    #[test]
    fn chromosomes_sort_by_number_then_name_once_each() {
        let d = setup("sort");
        let (_, out, _) = run_cmp(&d, &["--chroms", "chr22,X 21 2 22 10"]);
        let rows: Vec<&str> = out
            .lines()
            .skip(2)
            .take_while(|l| l.starts_with('|'))
            .map(|l| l.split(" | ").next().unwrap())
            .collect();
        assert_eq!(rows, ["| 2", "| 10", "| 21", "| 22", "| X", "| all"]);
    }

    #[test]
    fn a_label_given_twice_unites_its_files() {
        let d = setup("labels");
        put(
            &d,
            "more.tsv",
            "chrom\tstart\tend\tlength_bp\tancestry\n22\t350\t500\t150\tNeanderthal\n",
        );
        let more = format!("x={}", d.join("more.tsv").display());
        let (_, out, _) = run_cmp(&d, &["--mine", &more]);
        // 100-400 and 350-500 merge to 100-500: 400 bp, 200 of them released
        assert!(out.contains("| 22 | 300 | 175 | 300 | 400 |"), "{out}");
        assert!(out.contains("x: 400 bp, 1.14 times"), "{out}");
    }

    #[test]
    fn no_released_figure_gives_nan() {
        let d = setup("nan");
        let (_, out, _) = run_cmp(&d, &["--chroms", "1"]);
        assert!(
            out.contains("x: 0 bp, nan times the released A figure;"),
            "{out}"
        );
    }

    #[test]
    fn unknown_sample_stops() {
        let d = setup("unknown");
        let chen = d.join("chen.txt").to_string_lossy().into_owned();
        let mine = format!("x={}", d.join("mine.tsv").display());
        let (code, out, err) = call(
            compare_cli,
            &[
                "--chen-file",
                &chen,
                "--sample",
                "Z",
                "--chroms",
                "22",
                "--mine",
                &mine,
            ],
        );
        assert_eq!((code, out.as_str()), (1, ""));
        assert_eq!(err, format!("Z is not in {chen}\n"));
    }

    #[test]
    fn empty_group_stops() {
        let d = setup("nogroup");
        let (code, out, err) = run_cmp(&d, &["--group", "NOPE"]);
        assert_eq!((code, out.as_str()), (1, ""));
        assert!(err.contains("NOPE"), "{err}");
    }

    #[test]
    fn wrong_segments_file_stops() {
        let d = setup("wrong");
        put(&d, "mine.tsv", "a\tb\n");
        let (code, _, err) = run_cmp(&d, &[]);
        assert_eq!(code, 1);
        assert!(err.ends_with(": not a step 43 segments file\n"), "{err}");
        put(&d, "mine.tsv", "");
        assert_eq!(run_cmp(&d, &[]).0, 1);
    }

    #[test]
    fn no_chromosome_stops() {
        let d = setup("nochrom");
        let (code, _, err) = run_cmp(&d, &["--chroms", " , "]);
        assert_eq!((code, err.as_str()), (1, "--chroms names no chromosome\n"));
    }

    #[test]
    fn a_reversed_or_negative_interval_stops() {
        let d = setup("reversed");
        put(
            &d,
            "mine.tsv",
            &format!(
                "{}chr22\t400\t100\t-300\tNeanderthal\n",
                &MINE[..MINE.find('\n').unwrap() + 1]
            ),
        );
        let (code, out, err) = run_cmp(&d, &[]);
        assert_eq!((code, out.as_str()), (1, ""));
        assert!(
            err.starts_with("ERROR: ") && err.contains("400-100"),
            "{err}"
        );
        put(&d, "chen.txt", "22\t-5\t10\t9\t9\t15\tCEU\tEUR\tA\n");
        assert_eq!(run_cmp(&d, &[]).0, 1);
        put(&d, "chen.txt", "22\tx\t10\t9\t9\t15\tCEU\tEUR\tA\n");
        let (code, _, err) = run_cmp(&d, &[]);
        assert_eq!(code, 1);
        assert!(err.contains("\"x\" is not a whole number"), "{err}");
    }

    #[test]
    fn command_line_errors_exit_2() {
        let d = setup("usage");
        let chen = d.join("chen.txt").to_string_lossy().into_owned();
        let (code, _, err) = call(
            compare_cli,
            &["--chen-file", &chen, "--sample", "A", "--chroms", "22"],
        );
        assert_eq!(code, 2);
        assert!(err.contains("required: --mine"), "{err}");
        let (code, _, err) = call(
            compare_cli,
            &[
                "--chen-file",
                &chen,
                "--sample",
                "A",
                "--chroms",
                "22",
                "--mine",
                "nolabel.tsv",
            ],
        );
        assert_eq!(code, 2);
        assert!(err.contains("is not LABEL=SEGMENTS.tsv"), "{err}");
        let (code, out, _) = call(genomes_cli, &["--help"]);
        assert_eq!(code, 0);
        assert!(out.starts_with("usage: og archaic genomes"), "{out}");
        assert_eq!(call(genomes_cli, &["--og-scaffold-probe"]).0, 2);
    }

    fn genomes_dir(name: &str) -> (PathBuf, String, String) {
        let d = dir(name);
        let hdr = "ID\tchrom\tstart\tend\tslod\n";
        for (n, body) in [
            ("Altai", "SAMPLE\t22\t100000\t300000\t5\n"),
            ("Vindija", "SAMPLE\t22\t1000000\t1100000\t5\n"),
            ("Denisova", "SAMPLE\t22\t5000000\t5100000\t5\n"),
        ] {
            put(&d, &format!("chr22.{n}.txt"), &format!("{hdr}{body}"));
        }
        let spans = put(&d, "spans.tsv", "chrom\tstart\tend\n22\t1\t9999\n");
        let chen = put(
            &d,
            "chen_big.txt",
            "chr start end LOD maxLOD size pop anc ID\n\
             22\t100000\t300000\t9\t9\t200000\tCEU\tEUR\tA\n\
             22\t1000000\t1100000\t9\t9\t100000\tCEU\tEUR\tA\n",
        );
        (d, spans, chen)
    }

    #[test]
    fn genomes_union_row() {
        let (d, spans, chen) = genomes_dir("genomes");
        let dir = d.to_string_lossy().into_owned();
        let (code, out, err) = call(
            genomes_cli,
            &[
                "--chen-file",
                &chen,
                "--sample",
                "A",
                "--ibdmix-dir",
                &dir,
                "--spans",
                &spans,
            ],
        );
        assert_eq!((code, err.as_str()), (0, ""));
        // released A on chr22: 300,000 bp; Altai and Vindija together are exactly that
        assert_eq!(
            out,
            "| IBDmix output | bp | under the released A calls | only here | only released |\n\
             |---|---|---|---|---|\n\
             | Altai | 200000 | 200000 | 0 | 100000 |\n\
             | Vindija | 100000 | 100000 | 0 | 200000 |\n\
             | Denisova | 100000 | 0 | 100000 | 300000 |\n\
             | Altai+Vindija | 300000 | 300000 | 0 | 0 |\n\
             \n\
             Released A: 300000 bp on these 1 chromosomes.\n"
        );
    }

    #[test]
    fn genomes_keeps_lod_4_and_50_kb_only() {
        let (d, spans, chen) = genomes_dir("cutoffs");
        put(
            &d,
            "chr22.Altai.txt",
            "ID\tchrom\tstart\tend\tslod\n\
             SAMPLE\t22\t100000\t150000\t4\n\
             SAMPLE\t22\t200000\t249999\t9\n\
             SAMPLE\t22\t300000\t400000\t3.99\n\
             OTHER\t22\t500000\t600000\t9\n\
             SAMPLE\t22\t700000\n",
        );
        let dir = d.to_string_lossy().into_owned();
        let (code, out, err) = call(
            genomes_cli,
            &[
                "--chen-file",
                &chen,
                "--sample",
                "A",
                "--ibdmix-dir",
                &dir,
                "--spans",
                &spans,
            ],
        );
        assert_eq!(code, 0, "{err}");
        assert!(
            out.contains("| Altai | 50000 | 50000 | 0 | 250000 |"),
            "{out}"
        );
    }

    #[test]
    fn genomes_refuses_what_is_not_ibdmix_output() {
        let (d, spans, chen) = genomes_dir("notibdmix");
        let bad = put(&d, "chr22.Vindija.txt", "chrom\tstart\n");
        let dir = d.to_string_lossy().into_owned();
        let args = [
            "--chen-file",
            &chen,
            "--sample",
            "A",
            "--ibdmix-dir",
            &dir,
            "--spans",
            &spans,
        ];
        let (code, out, err) = call(genomes_cli, &args);
        assert_eq!((code, out.as_str()), (1, ""));
        assert_eq!(
            err,
            format!("ERROR: {bad} is not IBDmix output (header chrom start)\n")
        );
        fs::remove_file(&bad).unwrap();
        let (code, _, err) = call(genomes_cli, &args);
        assert_eq!(code, 1);
        assert!(err.starts_with("ERROR: cannot read "), "{err}");
        put(&d, "spans.tsv", "");
        let (code, _, err) = call(genomes_cli, &args);
        assert_eq!(code, 1);
        assert!(err.ends_with(" is empty\n"), "{err}");
    }

    #[test]
    fn genomes_unknown_sample_stops() {
        let (d, spans, chen) = genomes_dir("genomes-unknown");
        let dir = d.to_string_lossy().into_owned();
        let (code, _, err) = call(
            genomes_cli,
            &[
                "--chen-file",
                &chen,
                "--sample",
                "Z",
                "--ibdmix-dir",
                &dir,
                "--spans",
                &spans,
            ],
        );
        assert_eq!((code, err), (1, format!("Z is not in {chen}\n")));
    }
}
