//! `vep-af HITS.tsv OUT.vcf CHROM POS [--filler] [--want WANT.tsv] [--rsids]`:
//! a made-up VEP file with one record per row of step 06's hits table
//! (`<sample>_clinvar_hits.tsv`), each with a gnomAD exome AF in its CSQ
//! (3% on the planted hit at CHROM:POS, 0.01% times the row number on the
//! others) and half of it as the genome AF.
//!
//! - `--filler`: records every 97 bases from 2 Mb before to 2 Mb after the
//!   planted hit and 50 kb around each other hit, so a bgzipped copy spans
//!   many BGZF blocks and many windows of the linear index.
//! - `--want WANT.tsv`: `chrom:pos<TAB>af` of each hit, sorted.
//! - `--rsids`: each hit gets `rs<1000 + row>` as its ID and as VEP's
//!   Existing_variation, a seventh CSQ field.
//!
//! Stdout: `N VEP records`, and `, M at hits` with `--want`. The numbers
//! are Python's `repr` of the floats the heredoc computed, so the summary's
//! AFs read back as the same text.

use super::{Fail, bad, cell_int, create, each_line};
use indexmap::IndexMap;
use pycompat::fmt::{repr, round};
use std::io::Write;

/// The arguments.
struct Args<'a> {
    hits: &'a str,
    out: &'a str,
    chrom: &'a str,
    pos: i64,
    filler: bool,
    want: Option<&'a str>,
    rsids: bool,
}

fn parse(args: &[String]) -> Result<Args<'_>, Fail> {
    let mut pos = Vec::new();
    let (mut filler, mut want, mut rsids) = (false, None, false);
    let mut it = args.iter();
    while let Some(a) = it.next() {
        match a.as_str() {
            "--filler" => filler = true,
            "--rsids" => rsids = true,
            "--want" => {
                want = Some(
                    it.next()
                        .ok_or_else(|| Fail::Usage("--want takes a file".into()))?
                        .as_str(),
                );
            }
            _ => pos.push(a.as_str()),
        }
    }
    let [hits, out, chrom, p] = pos[..] else {
        return Err(Fail::Usage(format!(
            "takes HITS.tsv OUT.vcf CHROM POS, got {} positional arguments",
            pos.len()
        )));
    };
    let pos = p
        .trim()
        .parse()
        .map_err(|_| Fail::Usage(format!("POS '{p}' is not a whole number")))?;
    Ok(Args {
        hits,
        out,
        chrom,
        pos,
        filler,
        want,
        rsids,
    })
}

/// A record's key, as the heredoc's dict keyed it: `(chrom, pos, ref, alt)`.
type Key = (String, i64, String, String);

/// The allele VEP writes in CSQ: the ALT of an SNV or of an indel whose
/// first base is not the REF's, else the ALT without its anchor base (`-`
/// when nothing is left).
fn csq_allele(r: &str, alt: &str) -> String {
    let snv = r.chars().count() == 1 && alt.chars().count() == 1;
    let first = |s: &str| s.chars().next();
    if snv || first(alt) != first(r) {
        alt.to_string()
    } else {
        let rest: String = alt.chars().skip(1).collect();
        if rest.is_empty() { "-".into() } else { rest }
    }
}

pub fn af(args: &[String], out: &mut dyn Write, _err: &mut dyn Write) -> Result<(), Fail> {
    let a = parse(args)?;
    let mut rows: Vec<Vec<String>> = Vec::new();
    each_line(a.hits, |l| {
        rows.push(l.split('\t').map(str::to_string).collect());
        Ok(())
    })?;
    let mut recs: IndexMap<Key, String> = IndexMap::new();
    let mut expect = Vec::new();
    for (i, r) in rows.iter().skip(1).enumerate() {
        let [chrom, p, rf, alt, ..] = &r[..] else {
            return Err(bad(format!("{}: a row with fewer than 4 columns", a.hits)));
        };
        let p = cell_int(p)?;
        let n = i64::try_from(i).map_err(|_| bad("too many rows"))? + 1;
        let af = if *chrom == a.chrom && p == a.pos {
            0.03
        } else {
            #[allow(clippy::cast_precision_loss)]
            round(0.0001 * n as f64, 6).map_err(|e| bad(e.exception()))?
        };
        let allele = csq_allele(rf, alt);
        let (id, existing) = if a.rsids {
            (
                format!("rs{}", 1000 + n - 1),
                format!("|rs{}", 1000 + n - 1),
            )
        } else {
            (".".to_string(), String::new())
        };
        let line = format!(
            "{chrom}\t{p}\t{id}\t{rf}\t{alt}\t50\tPASS\tCSQ={allele}|x|LOW|G|{}|{}{existing}\tGT\t0/1",
            repr(af),
            repr(af / 2.0)
        );
        recs.insert((chrom.clone(), p, rf.clone(), alt.clone()), line);
        expect.push(format!("{chrom}:{p}\t{}", repr(af)));
    }
    if a.filler {
        let mut at: Vec<(String, i64)> = recs.keys().map(|k| (k.0.clone(), k.1)).collect();
        at.sort();
        at.dedup();
        let hit: std::collections::HashSet<(String, i64)> = at.iter().cloned().collect();
        for (chrom, pos) in &at {
            let span = if *chrom == a.chrom && *pos == a.pos {
                2_000_000
            } else {
                50_000
            };
            let mut q = (pos - span).max(1);
            while q < pos + span {
                if !hit.contains(&(chrom.clone(), q)) {
                    recs.entry((chrom.clone(), q, "A".into(), "C".into()))
                        .or_insert_with(|| {
                            format!(
                                "{chrom}\t{q}\t.\tA\tC\t50\tPASS\tCSQ=C|x|LOW|G|0.4|0.4\tGT\t0/1"
                            )
                        });
                }
                q += 97;
            }
        }
    }
    let mut contigs: Vec<&str> = recs.keys().map(|k| k.0.as_str()).collect();
    contigs.sort_by(|x, y| (x.chars().count(), *x).cmp(&(y.chars().count(), *y)));
    contigs.dedup();
    let order: std::collections::HashMap<&str, usize> =
        contigs.iter().enumerate().map(|(i, c)| (*c, i)).collect();
    let mut f = create(a.out)?;
    f.write_all(b"##fileformat=VCFv4.2\n")?;
    for c in &contigs {
        writeln!(f, "##contig=<ID={c}>")?;
    }
    let fields = if a.rsids {
        "Allele|Consequence|IMPACT|SYMBOL|gnomADe_AF|gnomADg_AF|Existing_variation"
    } else {
        "Allele|Consequence|IMPACT|SYMBOL|gnomADe_AF|gnomADg_AF"
    };
    writeln!(
        f,
        "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: {fields}\">"
    )?;
    f.write_all(b"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")?;
    f.write_all(b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n")?;
    let mut sorted: Vec<(&Key, &String)> = recs.iter().collect();
    sorted.sort_by_key(|(k, _)| (order[k.0.as_str()], k.1));
    for (_, line) in sorted {
        writeln!(f, "{line}")?;
    }
    f.flush()?;
    if let Some(w) = a.want {
        expect.sort();
        let mut wf = create(w)?;
        writeln!(wf, "{}", expect.join("\n"))?;
        wf.flush()?;
        writeln!(out, "{} VEP records, {} at hits", recs.len(), expect.len())?;
    } else {
        writeln!(out, "{} VEP records", recs.len())?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::super::testutil::{arg, dir, run};

    const HITS: &str = "chrom\tpos\tref\talt\tgene\n\
        chr2\t500\tA\tG\tX\n\
        chr10\t70\tAC\tA\tY\n\
        chr2\t100\tC\tCTT\tZ\n\
        chr2\t300\tG\tTA\tW\n";

    #[test]
    fn rsids_mode_writes_the_records_the_heredoc_wrote() {
        let d = dir("vep-rsids");
        std::fs::write(d.join("hits.tsv"), HITS).expect("hits");
        let (code, out, err) = run(
            "vep-af",
            &[
                &arg(&d, "hits.tsv"),
                &arg(&d, "out.vcf"),
                "chr2",
                "100",
                "--rsids",
            ],
        );
        assert_eq!(code, 0, "{err}");
        assert_eq!(out, "4 VEP records\n");
        let vcf = std::fs::read_to_string(d.join("out.vcf")).expect("vcf");
        assert_eq!(
            vcf,
            "##fileformat=VCFv4.2\n##contig=<ID=chr2>\n##contig=<ID=chr10>\n\
##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. \
Format: Allele|Consequence|IMPACT|SYMBOL|gnomADe_AF|gnomADg_AF|Existing_variation\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n\
chr2\t100\trs1002\tC\tCTT\t50\tPASS\tCSQ=TT|x|LOW|G|0.03|0.015|rs1002\tGT\t0/1\n\
chr2\t300\trs1003\tG\tTA\t50\tPASS\tCSQ=TA|x|LOW|G|0.0004|0.0002|rs1003\tGT\t0/1\n\
chr2\t500\trs1000\tA\tG\t50\tPASS\tCSQ=G|x|LOW|G|0.0001|5e-05|rs1000\tGT\t0/1\n\
chr10\t70\trs1001\tAC\tA\t50\tPASS\tCSQ=-|x|LOW|G|0.0002|0.0001|rs1001\tGT\t0/1\n"
        );
        assert!(!d.join("want.tsv").exists());
    }

    #[test]
    fn filler_mode_spreads_records_and_writes_the_wanted_afs() {
        let d = dir("vep-filler");
        std::fs::write(
            d.join("hits.tsv"),
            "c\tp\tr\ta\nchr1\t200\tA\tG\nchr1\t2000100\tT\tC\n",
        )
        .expect("hits");
        let (code, out, err) = run(
            "vep-af",
            &[
                &arg(&d, "hits.tsv"),
                &arg(&d, "out.vcf"),
                "chr1",
                "2000100",
                "--filler",
                "--want",
                &arg(&d, "want.tsv"),
            ],
        );
        assert_eq!(code, 0, "{err}");
        let vcf = std::fs::read_to_string(d.join("out.vcf")).expect("vcf");
        let recs: Vec<&str> = vcf.lines().filter(|l| !l.starts_with('#')).collect();
        // Around 200: from 1 below 50200 every 97 (518 positions). Around
        // 2000100: from 100 below 4000100 every 97 (41238 positions).
        let mut pos: Vec<i64> = recs
            .iter()
            .map(|l| l.split('\t').nth(1).expect("pos").parse().expect("int"))
            .collect();
        let n = pos.len();
        pos.dedup();
        assert_eq!(pos.len(), n, "one record per position");
        assert!(pos.windows(2).all(|w| w[0] < w[1]), "sorted");
        assert_eq!(n, 518 + 41_238 + 2);
        assert_eq!(out, format!("{n} VEP records, 2 at hits\n"));
        assert!(
            recs.contains(&"chr1\t200\t.\tA\tG\t50\tPASS\tCSQ=G|x|LOW|G|0.0001|5e-05\tGT\t0/1")
        );
        assert!(
            recs.contains(&"chr1\t2000100\t.\tT\tC\t50\tPASS\tCSQ=C|x|LOW|G|0.03|0.015\tGT\t0/1")
        );
        assert!(recs.contains(&"chr1\t1\t.\tA\tC\t50\tPASS\tCSQ=C|x|LOW|G|0.4|0.4\tGT\t0/1"));
        assert_eq!((pos[0], pos[n - 1]), (1, 100 + 97 * 41_237));
        assert_eq!(
            std::fs::read_to_string(d.join("want.tsv")).expect("want"),
            "chr1:200\t0.0001\nchr1:2000100\t0.03\n"
        );
    }

    #[test]
    fn a_row_without_four_columns_stops_the_kind() {
        let d = dir("vep-short");
        std::fs::write(d.join("hits.tsv"), "h\nchr1\t5\n").expect("hits");
        let (code, _, err) = run(
            "vep-af",
            &[&arg(&d, "hits.tsv"), &arg(&d, "o.vcf"), "chr1", "5"],
        );
        assert_eq!(code, 1);
        assert!(err.contains("fewer than 4 columns"), "{err}");
    }
}
