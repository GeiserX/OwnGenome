//! The two consumer-export checks of `tests/e2e/consumer-export-1.sh` that
//! read a gVCF a second way, independent of `og consumer export`.
//!
//! - `consumer-expected GVCF SITES PAR`: `rsid<TAB>letters` for each site of
//!   the 23andMe site table, from the gVCF records over it, looked up with a
//!   bisect among all records of the contig (not the export's walk).
//! - `consumer-roundtrip LIFTED CALLS SITES CHIP`: `called lifted moved
//!   compared nonref mismatches` for a chip file read back by
//!   `scripts/chip-to-vcf.sh` and lifted to GRCh38, against the case's own
//!   calls. chrM sites are left out.

use super::{Fail, bad, cell_int, each_line, gt_parts, repr_list};
use std::collections::{BTreeSet, HashMap};
use std::io::Write;

/// One gVCF record: start, end (`END=` or the REF's last base), REF and
/// ALTs upper-cased, FILTER and GT.
struct Rec {
    pos: i64,
    end: i64,
    r: String,
    alts: Vec<String>,
    filter: String,
    gt: String,
}

/// The number after the first `END=` that starts INFO or follows a `;`, as
/// `re.search(r"(?:^|;)END=(\d+)", info)` finds it.
fn info_end(info: &str) -> Option<i64> {
    let mut from = 0;
    while let Some(i) = info[from..].find("END=") {
        let at = from + i;
        if at == 0 || info.as_bytes()[at - 1] == b';' {
            let digits: String = info[at + 4..]
                .chars()
                .take_while(char::is_ascii_digit)
                .collect();
            if !digits.is_empty() {
                return digits.parse().ok();
            }
        }
        from = at + 1;
    }
    None
}

/// Python's `bisect.bisect_right(a, x)`, step for step.
fn bisect_right(a: &[i64], x: i64) -> usize {
    let (mut lo, mut hi) = (0, a.len());
    while lo < hi {
        let mid = usize::midpoint(lo, hi);
        if x < a[mid] {
            hi = mid;
        } else {
            lo = mid + 1;
        }
    }
    lo
}

/// An ALT that stands for no allele: `<*>`, `<NON_REF>` or `.`.
fn symbolic(a: &str) -> bool {
    matches!(a, "<*>" | "<NON_REF>" | ".")
}

/// A GT with an allele that is not `0`, `.` or empty.
fn nonref(gt: &str) -> bool {
    gt_parts(gt).iter().any(|i| !matches!(*i, "0" | "." | ""))
}

/// The first letter of allele `i` of `alleles`, as `alleles[int(i)][0]`.
fn letter(alleles: &[&str], i: &str) -> Result<String, Fail> {
    let k: i64 = cell_int(i)?;
    let n = i64::try_from(alleles.len()).map_err(|_| bad("too many alleles"))?;
    let at = if k < 0 { n + k } else { k };
    let a = usize::try_from(at)
        .ok()
        .and_then(|k| alleles.get(k))
        .ok_or_else(|| bad(format!("allele index '{i}' out of range")))?;
    a.chars()
        .next()
        .map(String::from)
        .ok_or_else(|| bad("an empty allele"))
}

pub fn expected(args: &[String], out: &mut dyn Write, _err: &mut dyn Write) -> Result<(), Fail> {
    let [gvcf, sites, par_bed] = args else {
        return Err(Fail::Usage(format!(
            "takes 3 arguments, got {}",
            args.len()
        )));
    };
    let mut par: Vec<(String, i64, i64)> = Vec::new();
    each_line(par_bed, |l| {
        if !l.trim().is_empty() {
            let c: Vec<&str> = l.split_whitespace().collect();
            if c.len() < 3 {
                return Err(bad(format!("{par_bed}: a row with fewer than 3 columns")));
            }
            par.push((c[0].to_string(), cell_int(c[1])? + 1, cell_int(c[2])?));
        }
        Ok(())
    })?;
    let mut recs: HashMap<String, Vec<Rec>> = HashMap::new();
    each_line(gvcf, |l| {
        if l.starts_with('#') {
            return Ok(());
        }
        let c: Vec<&str> = l.split('\t').collect();
        if c.len() < 10 {
            return Err(bad(format!("{gvcf}: a record with fewer than 10 columns")));
        }
        let pos = cell_int(c[1])?;
        let r = c[3].to_uppercase();
        let end = match info_end(c[7]) {
            Some(e) => e,
            None => pos + i64::try_from(r.chars().count()).map_err(|_| bad("a long REF"))? - 1,
        };
        let keys: Vec<&str> = c[8].split(':').collect();
        let vals: Vec<&str> = c[9].split(':').collect();
        let gt = keys
            .iter()
            .position(|k| *k == "GT")
            .and_then(|i| vals.get(i))
            .ok_or_else(|| {
                bad(format!(
                    "{gvcf}: a record without a GT at {}:{}",
                    c[0], c[1]
                ))
            })?;
        recs.entry(c[0].to_string()).or_default().push(Rec {
            pos,
            end,
            alts: c[4].to_uppercase().split(',').map(str::to_string).collect(),
            r,
            filter: c[6].to_string(),
            gt: (*gt).to_string(),
        });
        Ok(())
    })?;
    let starts: HashMap<&str, Vec<i64>> = recs
        .iter()
        .map(|(k, v)| (k.as_str(), v.iter().map(|r| r.pos).collect()))
        .collect();
    let empty: Vec<Rec> = Vec::new();
    each_line(sites, |l| {
        if l.starts_with('#') {
            return Ok(());
        }
        let c: Vec<&str> = l.split('\t').collect();
        let [rsid, chrom, _, _, pos38, ref38, strand] = c[..] else {
            return Err(bad(format!("{sites}: a row without 7 columns")));
        };
        let p = cell_int(pos38)?;
        let contig = if chrom == "MT" {
            "chrM".to_string()
        } else {
            format!("chr{chrom}")
        };
        let lst = recs.get(&contig).unwrap_or(&empty);
        let j = bisect_right(
            starts.get(contig.as_str()).map_or(&[][..], Vec::as_slice),
            p,
        );
        let cover: Vec<&Rec> = lst[j.saturating_sub(50)..j]
            .iter()
            .filter(|r| r.pos <= p && p <= r.end)
            .collect();
        let all_symbolic = |r: &Rec| r.alts.iter().all(|a| symbolic(a));
        let inside = cover
            .iter()
            .any(|r| r.pos < p && !all_symbolic(r) && r.r.chars().count() > 1 && nonref(&r.gt));
        let at: Vec<&&Rec> = cover
            .iter()
            .filter(|r| r.pos == p && !all_symbolic(r))
            .collect();
        let blk: Vec<&&Rec> = cover.iter().filter(|r| all_symbolic(r)).collect();
        let mut letters: Option<Vec<String>> = None;
        if inside {
        } else if let Some(r) = at.first() {
            let idx = gt_parts(&r.gt);
            let filtered =
                !matches!(r.filter.as_str(), "PASS" | ".") && idx.iter().any(|i| *i != "0");
            // The REF's first base is read only when the call is kept, as
            // the heredoc's `and` chain read it.
            let first = || {
                r.r.chars()
                    .next()
                    .map(String::from)
                    .ok_or_else(|| bad(format!("an empty REF at {contig}:{}", r.pos)))
            };
            if !idx.contains(&".") && !filtered && first()? == ref38 {
                let alleles: Vec<&str> = std::iter::once(r.r.as_str())
                    .chain(r.alts.iter().map(String::as_str))
                    .collect();
                let bases = idx
                    .iter()
                    .map(|i| letter(&alleles, i))
                    .collect::<Result<Vec<String>, Fail>>()?;
                if bases.iter().all(|b| "ACGT".contains(b.as_str())) {
                    letters = Some(bases);
                }
            }
        } else if let Some(b) = blk.last() {
            let idx = gt_parts(&b.gt);
            if idx.iter().all(|i| *i == "0") {
                letters = Some(vec![ref38.to_string(); idx.len()]);
            }
        }
        if let Some(l) = &letters {
            let in_par = par
                .iter()
                .any(|(c, s, e)| *c == contig && *s <= p && p <= *e);
            let haploid = chrom == "MT" || (matches!(chrom, "X" | "Y") && !in_par);
            if haploid && l.len() == 2 {
                letters = (l[0] == l[1]).then(|| vec![l[0].clone()]);
            }
        }
        if let Some(l) = &letters
            && strand == "-"
        {
            letters = Some(
                l.iter()
                    .map(|b| match b.as_str() {
                        "A" => Ok("T".to_string()),
                        "C" => Ok("G".to_string()),
                        "G" => Ok("C".to_string()),
                        "T" => Ok("A".to_string()),
                        other => Err(bad(format!("no complement of '{other}' at {rsid}"))),
                    })
                    .collect::<Result<Vec<String>, Fail>>()?,
            );
        }
        let call = match letters {
            Some(mut l) if !l.is_empty() => {
                l.sort();
                l.concat()
            }
            _ => "--".to_string(),
        };
        writeln!(out, "{rsid}\t{call}")?;
        Ok(())
    })?;
    Ok(())
}

/// The first letters of the alleles a GT names, without missing ones.
fn called_letters(alleles: &[&str], gt: &str) -> Result<BTreeSet<String>, Fail> {
    gt_parts(gt)
        .into_iter()
        .filter(|i| !matches!(*i, "." | ""))
        .map(|i| letter(alleles, i))
        .collect()
}

pub fn roundtrip(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> Result<(), Fail> {
    let [lifted, calls, sites, txt] = args else {
        return Err(Fail::Usage(format!(
            "takes 4 arguments, got {}",
            args.len()
        )));
    };
    let mut table: HashMap<String, Vec<String>> = HashMap::new();
    each_line(sites, |l| {
        if !l.starts_with('#') {
            let c: Vec<String> = l.split('\t').map(str::to_string).collect();
            if c.len() < 6 {
                return Err(bad(format!("{sites}: a row with fewer than 6 columns")));
            }
            table.insert(c[0].clone(), c);
        }
        Ok(())
    })?;
    let mut dv: HashMap<(String, i64), BTreeSet<String>> = HashMap::new();
    each_line(calls, |l| {
        let c: Vec<&str> = l.split('\t').collect();
        let [chrom, pos, r, alt, flt, gt] = c[..] else {
            return Err(bad(format!("{calls}: a row without 6 columns")));
        };
        if flt == "PASS" {
            let al: Vec<&str> = std::iter::once(r).chain(alt.split(',')).collect();
            dv.insert(
                (chrom.to_string(), cell_int(pos)?),
                called_letters(&al, gt)?,
            );
        }
        Ok(())
    })?;
    let mut called = 0u64;
    each_line(txt, |l| {
        let c: Vec<&str> = l.split('\t').collect();
        if l.starts_with('#') || c[0] == "rsid" {
            return Ok(());
        }
        match table.get(c[0]) {
            Some(t) if t[1] != "MT" => {}
            _ => return Ok(()),
        }
        let last = c[c.len() - 1];
        if !matches!(last, "--" | "0") && !c.iter().skip(3).any(|x| *x == "0") {
            called += 1;
        }
        Ok(())
    })?;
    let (mut got, mut moved, mut compared, mut nonref_n, mut bad_n) =
        (0u64, 0u64, 0u64, 0u64, 0u64);
    each_line(lifted, |l| {
        let c: Vec<&str> = l.split('\t').collect();
        let [rid, chrom, pos, gt, r, alt] = c[..] else {
            return Err(bad(format!("{lifted}: a row without 6 columns")));
        };
        let Some(t) = table.get(rid) else {
            return Ok(());
        };
        if t[1] == "MT" {
            return Ok(());
        }
        let idx: Vec<&str> = gt_parts(gt)
            .into_iter()
            .filter(|i| !matches!(*i, "." | ""))
            .collect();
        if idx.is_empty() {
            return Ok(());
        }
        got += 1;
        let p = cell_int(pos)?;
        if p != cell_int(&t[4])? {
            moved += 1;
            if moved <= 5 {
                writeln!(
                    err,
                    "moved {rid}: GRCh37 {}:{} lifts to {chrom}:{pos}, dbSNP's GRCh38 position is {}",
                    t[1], t[2], t[4]
                )?;
            }
            return Ok(());
        }
        let al: Vec<&str> = std::iter::once(r).chain(alt.split(',')).collect();
        let have = idx
            .iter()
            .map(|i| letter(&al, i))
            .collect::<Result<BTreeSet<String>, Fail>>()?;
        let reference = BTreeSet::from([t[5].clone()]);
        let want = dv
            .get(&(chrom.to_string(), p))
            .cloned()
            .unwrap_or_else(|| reference.clone());
        compared += 1;
        nonref_n += u64::from(want != reference);
        if have != want {
            bad_n += 1;
            if bad_n <= 5 {
                let list = |s: &BTreeSet<String>| repr_list(&s.iter().cloned().collect::<Vec<_>>());
                writeln!(
                    err,
                    "mismatch {rid} {chrom}:{pos} lifted {} calls {}",
                    list(&have),
                    list(&want)
                )?;
            }
        }
        Ok(())
    })?;
    writeln!(out, "{called} {got} {moved} {compared} {nonref_n} {bad_n}")?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::super::testutil::{arg, dir, gz, run};

    const HEAD: &str =
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n";

    /// A gVCF and site table covering every branch of the reading: a SNV,
    /// a filtered and a missing call, a reference block, a deletion over a
    /// site, chrX outside and inside the PAR, chrM, a minus-strand site and
    /// a site no record covers.
    fn inputs(d: &std::path::Path) {
        let mut g = String::from(HEAD);
        for l in [
            "chr1\t100\t.\tA\tG,<*>\t30\tPASS\t.\tGT:GQ\t0/1:30",
            "chr1\t101\t.\tC\t<*>\t.\t.\tEND=150\tGT:GQ\t0/0:40",
            "chr1\t160\t.\tTAC\tT,<*>\t30\tPASS\t.\tGT:GQ\t1/1:30",
            "chr1\t200\t.\tG\tA,<*>\t30\tLowQual\t.\tGT:GQ\t0/1:30",
            "chr1\t210\t.\tC\tT,<*>\t30\tPASS\t.\tGT:GQ\t./1:30",
            "chr1\t220\t.\tc\tg,<*>\t30\tPASS\t.\tGT:GQ\t1|1:30",
            "chr1\t230\t.\tA\tN,<*>\t30\tPASS\t.\tGT:GQ\t0/1:30",
            "chrX\t1000\t.\tG\tA,<*>\t30\tPASS\t.\tGT:GQ\t1/1:30",
            "chrX\t1001\t.\tG\t<*>\t.\t.\tEND=1100\tGT:GQ\t0/0:30",
            "chrX\t3000\t.\tT\tC,<*>\t30\tPASS\t.\tGT:GQ\t0/1:30",
            "chrM\t50\t.\tA\tG\t30\tPASS\t.\tGT\t1",
        ] {
            g.push_str(l);
            g.push('\n');
        }
        gz(&d.join("g.vcf.gz"), &g);
        let mut s = String::from("# rsid\tchrom\tpos37\tref37\tpos38\tref38\tstrand\n");
        for l in [
            "rs1\t1\t9\tA\t100\tA\t+",
            "rs2\t1\t9\tC\t120\tC\t+",
            "rs3\t1\t9\tC\t120\tC\t-",
            "rs4\t1\t9\tA\t161\tA\t+",
            "rs5\t1\t9\tG\t200\tG\t+",
            "rs6\t1\t9\tC\t210\tC\t+",
            "rs7\t1\t9\tC\t220\tC\t+",
            "rs8\t1\t9\tA\t230\tA\t+",
            "rs9\t1\t9\tA\t999\tA\t+",
            "rs10\tX\t9\tG\t1000\tG\t+",
            "rs11\tX\t9\tG\t1050\tG\t+",
            "rs12\tX\t9\tT\t3000\tT\t+",
            "rs13\tMT\t9\tA\t50\tA\t+",
            "rs14\t1\t9\tA\t100\tT\t+",
        ] {
            s.push_str(l);
            s.push('\n');
        }
        gz(&d.join("sites.tsv.gz"), &s);
        std::fs::write(d.join("par.bed"), "chrX\t2999\t4000\n\n").expect("par");
    }

    #[test]
    fn expected_reads_every_kind_of_record_as_the_heredoc_did() {
        let d = dir("expected");
        inputs(&d);
        let (code, out, err) = run(
            "consumer-expected",
            &[
                &arg(&d, "g.vcf.gz"),
                &arg(&d, "sites.tsv.gz"),
                &arg(&d, "par.bed"),
            ],
        );
        assert_eq!(code, 0, "{err}");
        assert_eq!(
            out,
            "rs1\tAG\nrs2\tCC\nrs3\tGG\nrs4\t--\nrs5\t--\nrs6\t--\nrs7\tGG\nrs8\t--\n\
rs9\t--\nrs10\tA\nrs11\tG\nrs12\tCT\nrs13\tG\nrs14\t--\n"
        );
    }

    #[test]
    fn expected_stops_on_a_record_without_a_gt() {
        let d = dir("expected-nogt");
        inputs(&d);
        gz(
            &d.join("g.vcf.gz"),
            &format!("{HEAD}chr1\t5\t.\tA\tG\t1\tPASS\t.\tGQ\t30\n"),
        );
        let (code, _, err) = run(
            "consumer-expected",
            &[
                &arg(&d, "g.vcf.gz"),
                &arg(&d, "sites.tsv.gz"),
                &arg(&d, "par.bed"),
            ],
        );
        assert_eq!(code, 1);
        assert!(err.contains("without a GT"), "{err}");
    }

    #[test]
    fn roundtrip_counts_called_lifted_moved_and_mismatched_sites() {
        let d = dir("roundtrip");
        gz(
            &d.join("sites.tsv.gz"),
            "# header\nrs1\t1\t10\tA\t100\tA\t+\nrs2\t1\t20\tC\t200\tC\t+\nrs3\t1\t30\tG\t300\tG\t+\n\
rs4\tMT\t40\tT\t40\tT\t+\nrs5\t1\t50\tT\t500\tT\t+\n",
        );
        std::fs::write(
            d.join("calls.tsv"),
            "chr1\t100\tA\tG\tPASS\t0/1\nchr1\t200\tC\tT\tLowQual\t1/1\nchr1\t500\tT\tA,C\tPASS\t1/2\n",
        )
        .expect("calls");
        std::fs::write(
            d.join("chip.txt"),
            "# chip\nrsid\tchromosome\tposition\tgenotype\nrs1\t1\t10\tAG\nrs2\t1\t20\tCC\n\
rs3\t1\t30\t--\nrs4\tMT\t40\tT\nrs5\t1\t50\tAC\nrs9\t1\t9\tAA\n",
        )
        .expect("chip");
        std::fs::write(
            d.join("lifted.tsv"),
            "rs1\tchr1\t100\t0/1\tA\tG\nrs2\tchr1\t200\t0/0\tC\tT\nrs3\tchr1\t301\t1/1\tG\tA\n\
rs4\tchrM\t40\t0\tT\t.\nrs5\tchr1\t500\t0/1\tT\tA\nrs6\tchr1\t600\t0/0\tA\t.\nrs1\tchr1\t100\t./.\tA\tG\n",
        )
        .expect("lifted");
        let (code, out, err) = run(
            "consumer-roundtrip",
            &[
                &arg(&d, "lifted.tsv"),
                &arg(&d, "calls.tsv"),
                &arg(&d, "sites.tsv.gz"),
                &arg(&d, "chip.txt"),
            ],
        );
        assert_eq!(code, 0, "{err}");
        assert_eq!(out, "3 4 1 3 2 1\n");
        assert_eq!(
            err,
            "moved rs3: GRCh37 1:30 lifts to chr1:301, dbSNP's GRCh38 position is 300\n\
mismatch rs5 chr1:500 lifted ['A', 'T'] calls ['A', 'C']\n"
        );
    }
}
