//! The family case's somalier sites and its check of the chromosome map.

use super::{Fail, bad, cell_int, create, each_line, fai_contigs, want};
use pycompat::csv::{Quoting, Reader};
use std::collections::{BTreeSet, HashMap};
use std::io::Write;

/// `og dev fixture sites-vcf FAI SNVS OUT`: a sites VCF for somalier from
/// `CHROM POS REF ALT` rows (duplicates once, contigs FAI lacks left out),
/// in FAI's contig order, each `AF=0.5`. Prints the site count.
pub fn sites_vcf(args: &[String], out: &mut dyn Write, _err: &mut dyn Write) -> Result<(), Fail> {
    want(args, 3)?;
    let contigs = fai_contigs(&args[0])?;
    let order: HashMap<&str, usize> = contigs
        .iter()
        .enumerate()
        .map(|(i, c)| (c.0.as_str(), i))
        .collect();
    let mut keep = BTreeSet::new();
    each_line(&args[1], |l| {
        if l.trim().is_empty() {
            return Ok(());
        }
        let c: Vec<&str> = l.split('\t').collect();
        let [chrom, pos, r, a] = c[..] else {
            return Err(bad(format!("{}: a row without 4 columns", args[1])));
        };
        if let Some(&i) = order.get(chrom) {
            keep.insert((
                i,
                cell_int(pos)?,
                pos.to_string(),
                r.to_string(),
                a.to_string(),
            ));
        }
        Ok(())
    })?;
    let mut o = create(&args[2])?;
    writeln!(o, "##fileformat=VCFv4.2")?;
    for (c, len) in &contigs {
        writeln!(o, "##contig=<ID={c},length={len}>")?;
    }
    writeln!(o, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")?;
    for (i, _, pos, r, a) in &keep {
        writeln!(o, "{}\t{pos}\t.\t{r}\t{a}\t.\tPASS\tAF=0.5", contigs[*i].0)?;
    }
    o.flush()?;
    writeln!(out, "{} sites", keep.len())?;
    Ok(())
}

/// The DNA Painter CSV header step 38 writes.
const DNAPAINTER: [&str; 12] = [
    "Chr",
    "Start",
    "End",
    "cM",
    "SNPs",
    "Match",
    "Confidence",
    "Group",
    "Side",
    "Notes",
    "Colour",
    "Ancestor",
];

/// `og dev fixture segments-check SEGMENTS.tsv DNAPAINTER.csv`: checks
/// step 38's segments table and its DNA Painter CSV (the columns, each row
/// a paternal or maternal stretch of 20 sites or more on a `chr` contig,
/// and the CSV row for row the same stretches), then prints the paternal
/// and the maternal count.
pub fn segments_check(
    args: &[String],
    out: &mut dyn Write,
    _err: &mut dyn Write,
) -> Result<(), Fail> {
    want(args, 2)?;
    let (tsv, dp) = (&args[0], &args[1]);
    let mut head: Option<Vec<String>> = None;
    let mut rows: Vec<HashMap<String, String>> = Vec::new();
    each_line(tsv, |l| {
        let cells = l.split('\t').map(str::to_string);
        match &head {
            None => head = Some(cells.collect()),
            Some(h) => rows.push(h.iter().cloned().zip(cells).collect()),
        }
        Ok(())
    })?;
    let head = head.unwrap_or_default();
    let want_head = ["child", "side", "parent", "chrom", "start", "end"];
    if head.len() < 6 || head[..6] != want_head {
        return Err(bad(format!("{tsv}: header {head:?}")));
    }
    let get = |r: &HashMap<String, String>, k: &str| -> Result<String, Fail> {
        r.get(k)
            .cloned()
            .ok_or_else(|| bad(format!("{tsv}: a row without {k}")))
    };
    for r in &rows {
        let side = get(r, "side")?;
        let chrom = get(r, "chrom")?;
        let (s, e) = (cell_int(&get(r, "start")?)?, cell_int(&get(r, "end")?)?);
        if !(side == "paternal" || side == "maternal") || !chrom.starts_with("chr") {
            return Err(bad(format!("{tsv}: row {r:?}")));
        }
        if !(0 < s && s <= e && cell_int(&get(r, "sites")?)? >= 20) {
            return Err(bad(format!("{tsv}: row {r:?}")));
        }
    }
    let lines = crate::smoke::csv_lines(dp).map_err(bad)?;
    let csv = Reader::new(lines, ',', Quoting::Minimal)
        .collect::<Result<Vec<Vec<String>>, _>>()
        .map_err(|e| bad(format!("{dp}: {e}")))?;
    if csv.first().is_none_or(|h| h[..] != DNAPAINTER) {
        return Err(bad(format!("{dp}: header {:?}", csv.first())));
    }
    if csv.len() - 1 != rows.len() {
        return Err(bad(format!(
            "{dp}: {} rows, the table {}",
            csv.len() - 1,
            rows.len()
        )));
    }
    for (x, r) in csv[1..].iter().zip(&rows) {
        let chrom = get(r, "chrom")?;
        let want = [
            chrom.get(3..).unwrap_or("").to_string(),
            get(r, "start")?,
            get(r, "end")?,
            get(r, "side")?,
        ];
        let have = [0, 1, 2, 8].map(|i| x.get(i).cloned().unwrap_or_default());
        if x.len() < 9 || have != want {
            return Err(bad(format!("{dp}: row {x:?} against {r:?}")));
        }
    }
    let side = |s: &str| {
        rows.iter()
            .filter(|r| r.get("side").is_some_and(|x| x == s))
            .count()
    };
    writeln!(out, "{} {}", side("paternal"), side("maternal"))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::super::testutil::{dir, run};

    #[test]
    fn sites_vcf_sorts_by_the_fai_and_drops_doubles() {
        let d = dir("sites-vcf");
        std::fs::write(
            d.join("r.fai"),
            "chr2\t500\t1\t60\t61\nchr1\t900\t600\t60\t61\n",
        )
        .expect("fai");
        std::fs::write(
            d.join("snvs.tsv"),
            "chr1\t30\tA\tG\nchr2\t100\tC\tT\nchr1\t30\tA\tG\nchr9\t5\tA\tC\nchr2\t20\tG\tA\n\n",
        )
        .expect("snvs");
        let p = |n: &str| d.join(n).to_str().expect("utf8").to_string();
        let (code, out, err) = run("sites-vcf", &[&p("r.fai"), &p("snvs.tsv"), &p("o.vcf")]);
        assert_eq!((code, out.as_str()), (0, "3 sites\n"), "{err}");
        assert_eq!(
            std::fs::read_to_string(d.join("o.vcf")).expect("vcf"),
            "##fileformat=VCFv4.2\n##contig=<ID=chr2,length=500>\n##contig=<ID=chr1,length=900>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
             chr2\t20\t.\tG\tA\t.\tPASS\tAF=0.5\nchr2\t100\t.\tC\tT\t.\tPASS\tAF=0.5\nchr1\t30\t.\tA\tG\t.\tPASS\tAF=0.5\n"
        );
    }

    #[test]
    fn segments_check_counts_and_refuses_a_mismatch() {
        let d = dir("segments");
        std::fs::write(
            d.join("s.tsv"),
            "child\tside\tparent\tchrom\tstart\tend\tsites\nC\tpaternal\tF\tchr1\t10\t500\t25\nC\tmaternal\tM\tchr2\t5\t90\t20\nC\tpaternal\tF\tchr3\t1\t2\t40\n",
        )
        .expect("tsv");
        let head = "Chr,Start,End,cM,SNPs,Match,Confidence,Group,Side,Notes,Colour,Ancestor\n";
        std::fs::write(
            d.join("d.csv"),
            format!("{head}1,10,500,1.0,25,C,,,paternal,\"a, b\",,F\n2,5,90,1,20,C,,,maternal,,,M\n3,1,2,1,40,C,,,paternal,,,F\n"),
        )
        .expect("csv");
        let p = |n: &str| d.join(n).to_str().expect("utf8").to_string();
        let (code, out, err) = run("segments-check", &[&p("s.tsv"), &p("d.csv")]);
        assert_eq!((code, out.as_str()), (0, "2 1\n"), "{err}");
        std::fs::write(d.join("bad.csv"), format!("{head}1,10,500,1,25,C,,,maternal,,,F\n2,5,90,1,20,C,,,maternal,,,M\n3,1,2,1,40,C,,,paternal,,,F\n")).expect("csv");
        let (code, out, _) = run("segments-check", &[&p("s.tsv"), &p("bad.csv")]);
        assert_eq!((code, out.as_str()), (1, ""));
        std::fs::write(
            d.join("few.tsv"),
            "child\tside\tparent\tchrom\tstart\tend\tsites\nC\tpaternal\tF\tchr1\t10\t500\t19\n",
        )
        .expect("tsv");
        let (code, _, _) = run("segments-check", &[&p("few.tsv"), &p("d.csv")]);
        assert_eq!(code, 1);
    }
}
