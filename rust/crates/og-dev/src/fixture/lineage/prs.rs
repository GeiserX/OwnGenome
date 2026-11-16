//! The PRS cases' score sites and their pgsc_calc match-log count.

use super::{Fail, Fasta, bad, cell_int, each_line, int, want};
use pycompat::csv::{DictReader, Quoting};
use std::io::Write;

/// The chr20 stretch the score sites come from (inside the fixture's slice).
const FROM: i64 = 10_050_000;
const TO: i64 = 10_450_000;
/// The step between two hom-ref candidates.
const STEP: usize = 7_919;

/// Whether no position of sorted `var` is within `gap` of `p`.
fn clear(var: &[i64], p: i64, gap: i64) -> bool {
    let i = var.partition_point(|&v| v < p - gap);
    var.get(i).is_none_or(|&v| v > p + gap)
}

/// `og dev fixture prs-sites BED TRUTH FASTA HOMREF [SNVS [--no-ambiguous]]`:
/// score sites on chr20 between 10,050,000 and 10,450,000. First HOMREF
/// hom-ref sites (`homref POS BASE`): every 7,919 bases, 200 bases inside a
/// confident interval of BED and more than 200 bases from any truth record,
/// on an `ACGT` reference base. Then SNVS of the truth set's heterozygous
/// or hom-alt SNVs, evenly spread (`snv POS ALT REF`); `--no-ambiguous`
/// leaves out A/T and C/G pairs, which pgsc_calc drops.
pub fn sites(args: &[String], out: &mut dyn Write, _err: &mut dyn Write) -> Result<(), Fail> {
    let no_amb = args.iter().any(|a| a == "--no-ambiguous");
    let pos: Vec<&String> = args.iter().filter(|a| *a != "--no-ambiguous").collect();
    if !(4..=5).contains(&pos.len()) || (no_amb && pos.len() != 5) {
        return Err(Fail::Usage(
            "takes BED TRUTH FASTA HOMREF [SNVS [--no-ambiguous]]".into(),
        ));
    }
    let (bed, truth, fasta) = (pos[0].as_str(), pos[1].as_str(), pos[2].as_str());
    let homref: usize = int(pos[3], "HOMREF")?;
    let n_snv: usize = pos.get(4).map(|s| int(s, "SNVS")).transpose()?.unwrap_or(0);
    let mut iv = Vec::new();
    each_line(bed, |l| {
        if l.starts_with("chr20\t") {
            let c: Vec<&str> = l.split_whitespace().collect();
            if c.len() < 3 {
                return Err(bad(format!("{bed}: a line with fewer than 3 columns")));
            }
            iv.push((cell_int(c[1])?, cell_int(c[2])?));
        }
        Ok(())
    })?;
    let (mut var, mut snv) = (Vec::new(), Vec::new());
    each_line(truth, |l| {
        if !l.starts_with("chr20\t") {
            return Ok(());
        }
        let f: Vec<&str> = l.split('\t').collect();
        if f.len() < 10 {
            return Err(bad(format!("{truth}: a record with fewer than 10 columns")));
        }
        let p = cell_int(f[1])?;
        var.push(p);
        let gt = f[9].split(':').next().unwrap_or("").replace('|', "/");
        let amb = matches!(
            (f[3], f[4]),
            ("A", "T") | ("T", "A") | ("C", "G") | ("G", "C")
        );
        if f[3].len() == 1
            && f[4].len() == 1
            && !(no_amb && amb)
            && FROM < p
            && p < TO
            && (gt == "0/1" || gt == "1/1")
        {
            snv.push((p, f[4].to_string(), f[3].to_string()));
        }
        Ok(())
    })?;
    var.sort_unstable();
    let mut fa = Fasta::open(fasta)?;
    if fa.len("chr20").is_none() {
        return Err(bad(format!("{fasta}.fai: no chr20")));
    }
    let mut picked = 0;
    for p in (FROM..TO).step_by(STEP) {
        if picked == homref {
            break;
        }
        let inside = iv.iter().any(|&(s, e)| s + 200 < p && p <= e - 200);
        if !inside || !clear(&var, p, 200) {
            continue;
        }
        if let Some(b @ (b'A' | b'C' | b'G' | b'T')) = fa.base("chr20", p as u64)? {
            writeln!(out, "homref\t{p}\t{}", b as char)?;
            picked += 1;
        }
    }
    if n_snv > 0 {
        let step = (snv.len() / n_snv).max(1);
        for (p, alt, r) in snv.iter().step_by(step).take(n_snv) {
            writeln!(out, "snv\t{p}\t{alt}\t{r}")?;
        }
    }
    Ok(())
}

/// `og dev fixture x-ref-sites BLOCKS VARIANTS FASTA CHROM N`: up to N
/// positions in the middle of the gVCF's reference blocks (`POS END GT`
/// rows of `bcftools query`, GT 0/0) of 50 bases or more, each more than
/// 2,000 bases after the one before, more than 100 bases from any call with
/// an ALT, on an `ACGT` base: `POS BASE` per line. The counts go to stderr.
pub fn x_ref_sites(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> Result<(), Fail> {
    want(args, 5)?;
    let (blocks, variants, fasta, chrom) = (&args[0], &args[1], &args[2], &args[3]);
    let n: usize = int(&args[4], "N")?;
    let mut var = Vec::new();
    each_line(variants, |l| {
        if !l.trim().is_empty() {
            var.push(cell_int(l)?);
        }
        Ok(())
    })?;
    let (mut rows, mut refs) = (0usize, Vec::new());
    each_line(blocks, |l| {
        if l.trim().is_empty() {
            return Ok(());
        }
        rows += 1;
        let r: Vec<&str> = l.split('\t').collect();
        if r.len() == 3 && !r[1].is_empty() && r[1] != "." && (r[2] == "0/0" || r[2] == "0|0") {
            refs.push((cell_int(r[0])?, cell_int(r[1])?));
        }
        Ok(())
    })?;
    writeln!(
        err,
        "{rows} gVCF records on the {chrom} slice, {} 0/0 blocks, {} calls with an ALT",
        refs.len(),
        var.len()
    )?;
    var.sort_unstable();
    let mut fa = Fasta::open(fasta)?;
    let (mut last, mut picked) = (0i64, 0usize);
    for (start, end) in refs {
        if picked == n {
            break;
        }
        let p = (start + end).div_euclid(2);
        if end - start >= 50
            && p - last > 2000
            && clear(&var, p, 100)
            && let Some(b @ (b'A' | b'C' | b'G' | b'T')) = fa.base(chrom, p as u64)?
        {
            writeln!(out, "{p}\t{}", b as char)?;
            last = p;
            picked += 1;
        }
    }
    Ok(())
}

/// One filter of `csv-count`: a column and the values it must (or must
/// not) hold.
struct Filter {
    key: String,
    values: Vec<String>,
    negate: bool,
}

/// `og dev fixture csv-count FILE KEY=V1,V2 KEY!=V1,V2 ...`: the rows of a
/// CSV (plain or gzip, Python's `csv.DictReader`) whose KEY is one of the
/// values (`=`) or none of them (`!=`, a missing column counts as none).
/// pgsc_calc's match log is the file the PRS cases count.
pub fn csv_count(args: &[String], out: &mut dyn Write, _err: &mut dyn Write) -> Result<(), Fail> {
    let Some((file, rest)) = args.split_first() else {
        return Err(Fail::Usage("takes FILE and filters".into()));
    };
    let filters = rest
        .iter()
        .map(|a| {
            let (key, values, negate) = if let Some((k, v)) = a.split_once("!=") {
                (k, v, true)
            } else if let Some((k, v)) = a.split_once('=') {
                (k, v, false)
            } else {
                return Err(Fail::Usage(format!(
                    "'{a}' is not KEY=VALUES or KEY!=VALUES"
                )));
            };
            Ok(Filter {
                key: key.to_string(),
                values: values.split(',').map(str::to_string).collect(),
                negate,
            })
        })
        .collect::<Result<Vec<_>, _>>()?;
    let mut n = 0u64;
    for row in DictReader::new(
        crate::smoke::csv_lines(file).map_err(bad)?,
        ',',
        Quoting::Minimal,
    ) {
        let row = row.map_err(|e| bad(format!("{file}: {e}")))?;
        if filters.iter().all(|f| {
            let v = row.get(&f.key).flatten();
            let hit = v.is_some_and(|v| f.values.iter().any(|x| x == v));
            hit != f.negate
        }) {
            n += 1;
        }
    }
    writeln!(out, "{n}")?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::super::testutil::{dir, fasta, run};

    #[test]
    fn csv_count_counts_matches_and_quoted_cells() {
        let d = dir("csv-count");
        let f = d.join("log.csv");
        std::fs::write(
            &f,
            "accession,chr_name,match_status,note\nPGS1,X,matched,\"a, b\"\nPGS1,23,matched,x\nPGS1,X,unmatched,x\nPGS1,20,matched,x\nPGS2,X,matched,x\n",
        )
        .expect("csv");
        let f = f.to_str().expect("utf8");
        let (code, out, _) = run(
            "csv-count",
            &[f, "accession=PGS1", "chr_name=X,23", "match_status=matched"],
        );
        assert_eq!((code, out.as_str()), (0, "2\n"));
        let (_, out, _) = run(
            "csv-count",
            &[
                f,
                "accession=PGS1",
                "chr_name!=X,23",
                "match_status=matched",
            ],
        );
        assert_eq!(out, "1\n");
        let (_, out, _) = run("csv-count", &["/dev/null", "accession=PGS1"]);
        assert_eq!(out, "0\n");
        let (code, _, _) = run("csv-count", &[f, "accession"]);
        assert_eq!(code, 2);
    }

    #[test]
    fn prs_sites_picks_hom_ref_then_spread_snvs() {
        let d = dir("prs-sites");
        let fa = d.join("ref.fa");
        fasta(&fa, &[("chr1", 100), ("chr20", 10_460_000)], 5);
        std::fs::write(
            d.join("conf.bed"),
            "chr20\t10049000\t10200000\nchr1\t0\t100\n",
        )
        .expect("bed");
        // A record near the 2nd candidate (10,057,919) knocks it out.
        let mut truth = String::from("##x\n");
        for (p, r, a, gt) in [
            (10_057_950, "A", "T", "0/1"),
            (10_100_000, "C", "G", "0|1"),
            (10_120_000, "C", "T", "1/1"),
            (10_130_000, "C", "TT", "0/1"),
            (10_140_000, "G", "A", "1|1"),
            (10_150_000, "G", "A", "0/0"),
        ] {
            truth.push_str(&format!("chr20\t{p}\t.\t{r}\t{a}\t50\tPASS\t.\tGT\t{gt}\n"));
        }
        std::fs::write(d.join("truth.vcf"), truth).expect("truth");
        let p = |n: &str| d.join(n).to_str().expect("utf8").to_string();
        let (code, out, err) = run(
            "prs-sites",
            &[&p("conf.bed"), &p("truth.vcf"), &p("ref.fa"), "3", "3"],
        );
        assert_eq!(code, 0, "{err}");
        let lines: Vec<&str> = out.lines().collect();
        assert_eq!(lines.len(), 6, "{out}");
        assert!(lines[0].starts_with("homref\t10050000\t"), "{out}");
        assert!(lines[1].starts_with("homref\t10065838\t"), "{out}");
        assert_eq!(
            &lines[3..],
            [
                "snv\t10057950\tT\tA",
                "snv\t10100000\tG\tC",
                "snv\t10120000\tT\tC"
            ]
        );
        let (_, out, _) = run(
            "prs-sites",
            &[
                &p("conf.bed"),
                &p("truth.vcf"),
                &p("ref.fa"),
                "0",
                "2",
                "--no-ambiguous",
            ],
        );
        assert_eq!(out, "snv\t10120000\tT\tC\nsnv\t10140000\tA\tG\n");
    }

    #[test]
    fn x_ref_sites_takes_block_middles_away_from_calls() {
        let d = dir("x-ref");
        let fa = d.join("ref.fa");
        fasta(&fa, &[("chrX", 20_000)], 9);
        std::fs::write(
            d.join("blocks.tsv"),
            "100\t200\t0/0\n150\t160\t0/0\n3000\t3100\t0|0\n4600\t4700\t0/0\n7000\t7040\t0/0\n9000\t9200\t0/1\n12000\t12100\t0/0\n13000\t.\t0/0\n",
        )
        .expect("blocks");
        std::fs::write(d.join("vars.txt"), "12040\n\n").expect("vars");
        let p = |n: &str| d.join(n).to_str().expect("utf8").to_string();
        let (code, out, err) = run(
            "x-ref-sites",
            &[&p("blocks.tsv"), &p("vars.txt"), &p("ref.fa"), "chrX", "5"],
        );
        assert_eq!(code, 0, "{err}");
        let pos: Vec<&str> = out
            .lines()
            .map(|l| l.split('\t').next().expect("pos"))
            .collect();
        assert_eq!(pos, ["3050"], "{out}");
        assert_eq!(
            err,
            "8 gVCF records on the chrX slice, 6 0/0 blocks, 1 calls with an ALT\n"
        );
    }
}
