//! GIAB's genotypes at an ancestry panel's SNVs, as a VCF a step reads.

use super::{Confident, Fail, bad, cell_int, each_line, fai_contigs, read_bed};
use std::collections::{BTreeSet, HashMap};
use std::io::Write;

/// One member: its name, its calls by (contig, position) and its
/// confident regions.
struct Member {
    name: String,
    calls: HashMap<(String, i64), (String, String, String)>,
    bed: HashMap<String, Confident>,
}

/// `og dev fixture truth-genotypes SITES FAI [--chrom C] NAME=CALLS,BED...`:
/// a VCF of each member's genotype at the panel's SNVs (`CHROM POS REF ALT`
/// rows of SITES) that lie inside every member's confident regions (BED).
/// A site a member's calls (a VCF, header optional) do not list is 0/0; a
/// site where a call has other alleles or a genotype other than 0/0, 0/1
/// or 1/1 is left out. `--chrom` keeps one contig. The site count goes to
/// stderr.
pub fn genotypes(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> Result<(), Fail> {
    let mut chrom: Option<&str> = None;
    let mut pos = Vec::new();
    let mut it = args.iter();
    while let Some(a) = it.next() {
        if a == "--chrom" {
            chrom = Some(
                it.next()
                    .ok_or_else(|| Fail::Usage("--chrom needs a contig".into()))?,
            );
        } else {
            pos.push(a.as_str());
        }
    }
    if pos.len() < 3 {
        return Err(Fail::Usage(
            "takes SITES FAI [--chrom C] NAME=CALLS,BED...".into(),
        ));
    }
    let (sites, fai) = (pos[0], pos[1]);
    let mut members = Vec::new();
    for m in &pos[2..] {
        let (name, files) = m
            .split_once('=')
            .ok_or_else(|| Fail::Usage(format!("'{m}' is not NAME=CALLS,BED")))?;
        let (calls_path, bed_path) = files
            .split_once(',')
            .ok_or_else(|| Fail::Usage(format!("'{m}' is not NAME=CALLS,BED")))?;
        let mut calls = HashMap::new();
        each_line(calls_path, |l| {
            if l.starts_with('#') {
                return Ok(());
            }
            let f: Vec<&str> = l.split('\t').collect();
            if f.len() < 10 {
                return Err(bad(format!(
                    "{calls_path}: a record with fewer than 10 columns"
                )));
            }
            let i = f[8]
                .split(':')
                .position(|k| k == "GT")
                .ok_or_else(|| bad(format!("{calls_path}: a record without GT")))?;
            let gt = f[9].split(':').nth(i).unwrap_or("").replace('|', "/");
            let gt = if gt == "1/0" { "0/1".to_string() } else { gt };
            calls.insert(
                (f[0].to_string(), cell_int(f[1])?),
                (f[3].to_string(), f[4].to_string(), gt),
            );
            Ok(())
        })?;
        members.push(Member {
            name: name.to_string(),
            calls,
            bed: read_bed(bed_path)?,
        });
    }
    let mut rows = Vec::new();
    each_line(sites, |l| {
        let c: Vec<&str> = l.split('\t').collect();
        if c.len() != 4 {
            return Err(bad(format!("{sites}: a row without 4 columns")));
        }
        if chrom.is_none_or(|k| k == c[0]) {
            rows.push([
                c[0].to_string(),
                c[1].to_string(),
                c[2].to_string(),
                c[3].to_string(),
            ]);
        }
        Ok(())
    })?;
    let chroms: BTreeSet<&str> = rows.iter().map(|r| r[0].as_str()).collect();
    writeln!(out, "##fileformat=VCFv4.2")?;
    writeln!(out, "##FILTER=<ID=PASS,Description=\"All filters passed\">")?;
    for (c, len) in fai_contigs(fai)? {
        if chroms.contains(c.as_str()) {
            writeln!(out, "##contig=<ID={c},length={len}>")?;
        }
    }
    writeln!(
        out,
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    )?;
    let names: Vec<&str> = members.iter().map(|m| m.name.as_str()).collect();
    writeln!(
        out,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}",
        names.join("\t")
    )?;
    let mut n = 0u64;
    'site: for [c, p, r, a] in &rows {
        let p = cell_int(p)?;
        if !members
            .iter()
            .all(|m| m.bed.get(c).is_some_and(|b| b.has(p)))
        {
            continue;
        }
        let mut gts = Vec::with_capacity(members.len());
        for m in &members {
            match m.calls.get(&(c.clone(), p)) {
                None => gts.push("0/0"),
                Some((tr, ta, g))
                    if tr == r && ta == a && matches!(g.as_str(), "0/0" | "0/1" | "1/1") =>
                {
                    gts.push(g.as_str());
                }
                Some(_) => continue 'site,
            }
        }
        writeln!(
            out,
            "{c}\t{p}\t.\t{r}\t{a}\t.\tPASS\t.\tGT\t{}",
            gts.join("\t")
        )?;
        n += 1;
    }
    writeln!(err, "{n} sites")?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::super::testutil::{dir, run};

    #[test]
    fn genotypes_follow_the_truth_inside_every_members_regions() {
        let d = dir("truth");
        let w = |n: &str, s: &str| std::fs::write(d.join(n), s).expect("write");
        w(
            "sites.tsv",
            "chr1\t100\tA\tG\nchr1\t200\tC\tT\nchr1\t300\tG\tA\nchr1\t400\tT\tC\nchr2\t50\tA\tC\n",
        );
        w(
            "ref.fa.fai",
            "chr1\t1000\t6\t60\t61\nchr2\t500\t1030\t60\t61\nchr3\t9\t1600\t60\t61\n",
        );
        w(
            "a.vcf",
            "##h\n#CHROM\nchr1\t100\t.\tA\tG\t.\t.\t.\tGT\t0|0\nchr1\t200\t.\tC\tT\t.\t.\t.\tGT:DP\t1|0:4\nchr1\t300\t.\tG\tC\t.\t.\t.\tGT\t0/1\nchr2\t50\t.\tA\tC\t.\t.\t.\tDP:GT\t3:1/1\n",
        );
        w("a.bed", "chr1\t0\t350\nchr2\t0\t100\n");
        w("b.vcf", "chr1\t100\t.\tA\tG\t.\t.\t.\tGT\t1/2\n");
        w("b.bed", "chr1\t50\t1000\n");
        let p = |n: &str| d.join(n).to_str().expect("utf8").to_string();
        let a = format!("A={},{}", p("a.vcf"), p("a.bed"));
        let (code, out, err) = run("truth-genotypes", &[&p("sites.tsv"), &p("ref.fa.fai"), &a]);
        assert_eq!(code, 0, "{err}");
        let body: Vec<&str> = out.lines().filter(|l| !l.starts_with('#')).collect();
        assert_eq!(
            body,
            [
                "chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/0",
                "chr1\t200\t.\tC\tT\t.\tPASS\t.\tGT\t0/1",
                "chr2\t50\t.\tA\tC\t.\tPASS\t.\tGT\t1/1"
            ]
        );
        assert!(
            out.contains("##contig=<ID=chr1,length=1000>\n##contig=<ID=chr2,length=500>\n##FORMAT"),
            "{out}"
        );
        assert!(out.contains("FORMAT\tA\n"), "{out}");
        assert_eq!(err, "3 sites\n");
        let b = format!("B={},{}", p("b.vcf"), p("b.bed"));
        let (_, out, err) = run(
            "truth-genotypes",
            &[&p("sites.tsv"), &p("ref.fa.fai"), "--chrom", "chr1", &a, &b],
        );
        let body: Vec<&str> = out.lines().filter(|l| !l.starts_with('#')).collect();
        assert_eq!(
            body,
            ["chr1\t200\t.\tC\tT\t.\tPASS\t.\tGT\t0/1\t0/0"],
            "{err}"
        );
        assert!(!out.contains("chr2"), "{out}");
    }
}
