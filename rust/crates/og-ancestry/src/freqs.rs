//! `og ancestry-fine freqs`: each reference group's allele frequencies,
//! counted from the panel, in place of a supervised ADMIXTURE run at
//! install.
//!
//! Supervised ADMIXTURE with every reference person's group given fixes
//! each person's ancestry to that one group, so the likelihood it then
//! maximises over the frequencies splits into one binomial per group and
//! SNV: its maximum is the group's sample frequency. ADMIXTURE climbs to it
//! and stops at its own tolerance; this counts it exactly, clamped to
//! ADMIXTURE's bounds `[1e-5, 1 - 1e-5]`.
//!
//! From the LD-pruned VCF of the reference people and `samples.tsv`, it
//! writes, in `--out`:
//!
//! - `ref.bim`: one SNV per record, as `plink2 --make-bed --set-all-var-ids
//!   '@:#'` wrote it: chromosome without `chr`, `CHR:POS`, 0, the position,
//!   A1 (the ALT) and A2 (the REF);
//! - `ref.fam`: the people counted, in the VCF's order, as plink2 wrote it;
//! - `<level>.P` for `continent` and `region`, ADMIXTURE's layout: one row
//!   per SNV, one column per group, the frequency of the `.bim`'s A2 with
//!   six decimals (what `admixture -P` reads at run time);
//! - `levels.tsv`: each level, its K and its group codes in the `.P`'s
//!   column order, which is the order the groups first appear in the VCF,
//!   as ADMIXTURE numbers them;
//! - `ref.geno` (`og_admix::geno`): each person's copies of A1, for the
//!   projection's bootstrap over reference people.
//!
//! A person counts when `samples.tsv` flags them a reference, they have
//! both groups and `--exclude` does not list them; other VCF columns are
//! skipped. A genotype with a missing allele counts for neither allele. A
//! group with no called allele at a SNV gets 0.5.

use crate::{Args, Fail, Level, Pending, Result, Spec};
use og_admix::geno::RefGeno;
use std::collections::HashSet;
use std::io::Write;
use std::path::{Path, PathBuf};

pub(crate) const SPEC: Spec = Spec {
    verb: "freqs",
    usage: "usage: og ancestry-fine freqs --vcf FILE --samples samples.tsv --out DIR [--exclude FILE]\n",
    options: &["--vcf", "--samples", "--out", "--exclude"],
    required: &["--vcf", "--samples", "--out"],
};

/// ADMIXTURE's bounds on a frequency.
const LOW: f64 = 1e-5;
const HIGH: f64 = 1.0 - 1e-5;

/// The reference people of `samples.tsv` not in `exclude`: id to
/// (continental group, regional group).
fn references(
    path: &Path,
    exclude: &HashSet<String>,
) -> Result<std::collections::HashMap<String, [String; 2]>> {
    let lines = crate::read_lines(path)?;
    let header: Vec<&str> = lines
        .first()
        .map(|l| l.split('\t').collect())
        .unwrap_or_default();
    if header != crate::panel::HEADER {
        return Err(Fail(format!(
            "{} does not start with the header of samples.tsv",
            path.display()
        )));
    }
    let mut out = std::collections::HashMap::new();
    for (n, line) in lines.iter().enumerate().skip(1) {
        if line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 6 {
            return Err(Fail(format!(
                "{}: line {}: {} columns, not 6",
                path.display(),
                n + 1,
                f.len()
            )));
        }
        let groups = Level::ALL.map(|l| f[l.column()]);
        let grouped = groups.iter().all(|g| !g.is_empty() && *g != "-");
        if f[4] == "1" && grouped && !exclude.contains(f[0]) {
            out.insert(f[0].to_owned(), groups.map(str::to_owned));
        }
    }
    Ok(out)
}

/// The ids of an `--exclude` file: whitespace or commas between them,
/// `#` to the end of a line a comment.
fn excluded(path: &Path) -> Result<HashSet<String>> {
    let mut out = HashSet::new();
    for line in crate::read_lines(path)? {
        let line = line.split('#').next().unwrap_or("");
        out.extend(
            line.split(|c: char| c == ',' || c.is_whitespace())
                .filter(|w| !w.is_empty())
                .map(str::to_owned),
        );
    }
    Ok(out)
}

/// One level's groups: the code of each, in order of first appearance, and
/// the group of each counted person.
struct Groups {
    codes: Vec<String>,
    of: Vec<usize>,
}

impl Groups {
    fn new() -> Groups {
        Groups {
            codes: Vec::new(),
            of: Vec::new(),
        }
    }

    fn add(&mut self, code: &str) {
        let g = match self.codes.iter().position(|c| c == code) {
            Some(g) => g,
            None => {
                self.codes.push(code.to_owned());
                self.codes.len() - 1
            }
        };
        self.of.push(g);
    }
}

/// The copies of ALT (the `.bim`'s A1) in a diploid GT, `None` when an
/// allele is missing.
fn alt_copies(gt: &str) -> std::result::Result<Option<u8>, String> {
    let alleles: Vec<&str> = gt.split(['|', '/']).collect();
    if alleles.len() != 2 {
        return Err(format!("GT {gt} is not diploid"));
    }
    let mut copies = 0u8;
    for a in alleles {
        match a {
            "0" => {}
            "1" => copies += 1,
            "." => return Ok(None),
            _ => return Err(format!("GT {gt} has an allele other than 0 and 1")),
        }
    }
    Ok(Some(copies))
}

/// The frequency of A2 (the REF) among `called` alleles with `alt` ALTs,
/// clamped to ADMIXTURE's bounds.
pub fn a2_freq(alt: u32, called: u32) -> f64 {
    if called == 0 {
        return 0.5;
    }
    (f64::from(called - alt) / f64::from(called)).clamp(LOW, HIGH)
}

/// A base of a biallelic SNV.
fn base(s: &str) -> bool {
    matches!(s, "A" | "C" | "G" | "T")
}

pub(crate) fn main(a: &Args, out: &mut dyn Write) -> Result<()> {
    let vcf = Path::new(a.req("--vcf"));
    let exclude = match a.get("--exclude") {
        Some(p) => excluded(Path::new(p))?,
        None => HashSet::new(),
    };
    let refs = references(Path::new(a.req("--samples")), &exclude)?;
    // (VCF column, person) of each person counted, and both levels' groups.
    let mut cols: Vec<usize> = Vec::new();
    let mut fam = String::new();
    let mut levels = [Groups::new(), Groups::new()];
    let mut geno: Option<RefGeno> = None;
    let mut bim = String::new();
    let mut p_rows = [String::new(), String::new()];
    let mut n = 0u64;
    let mut row: Vec<Option<u8>> = Vec::new();
    crate::for_each_line(vcf, |line| {
        n += 1;
        if line.starts_with("##") {
            return Ok(());
        }
        let at = |what: String| Fail(format!("{}: line {n}: {what}", vcf.display()));
        if let Some(header) = line.strip_prefix('#') {
            if geno.is_some() {
                return Err(at("a second #CHROM line".to_owned()));
            }
            for (i, id) in header.split('\t').enumerate().skip(9) {
                if let Some(groups) = refs.get(id) {
                    cols.push(i);
                    fam.push_str(&format!("0\t{id}\t0\t0\t0\t-9\n"));
                    for (g, code) in levels.iter_mut().zip(groups) {
                        g.add(code);
                    }
                }
            }
            if cols.is_empty() {
                return Err(at("no reference person of the samples table".to_owned()));
            }
            geno = Some(RefGeno::new(cols.len()));
            return Ok(());
        }
        let Some(geno) = geno.as_mut() else {
            return Err(at("a record before the #CHROM line".to_owned()));
        };
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 10 {
            return Err(at(format!("{} columns, at least 10 needed", f.len())));
        }
        let (chrom, pos, ref_, alt) = (f[0], f[1], f[3], f[4]);
        if !base(ref_) || !base(alt) || ref_ == alt {
            return Err(at(format!("{ref_}>{alt} is not a biallelic SNV")));
        }
        let Some(gt_at) = f[8].split(':').position(|k| k == "GT") else {
            return Err(at("no GT in FORMAT".to_owned()));
        };
        row.clear();
        for &c in &cols {
            let gt = f
                .get(c)
                .and_then(|s| s.split(':').nth(gt_at))
                .ok_or_else(|| at(format!("no GT in column {}", c + 1)))?;
            row.push(alt_copies(gt).map_err(at)?);
        }
        geno.push_site(&row).map_err(|e| at(e.0))?;
        let chrom = chrom.strip_prefix("chr").unwrap_or(chrom);
        bim.push_str(&format!(
            "{chrom}\t{chrom}:{pos}\t0\t{pos}\t{alt}\t{ref_}\n"
        ));
        for (g, text) in levels.iter().zip(p_rows.iter_mut()) {
            let mut counts = vec![(0u32, 0u32); g.codes.len()];
            for (&grp, copies) in g.of.iter().zip(&row) {
                if let Some(c) = copies {
                    counts[grp].0 += u32::from(*c);
                    counts[grp].1 += 2;
                }
            }
            let cells: Vec<String> = counts
                .iter()
                .map(|&(alt, called)| format!("{:.6}", a2_freq(alt, called)))
                .collect();
            text.push_str(&cells.join(" "));
            text.push('\n');
        }
        Ok(())
    })?;
    let Some(geno) = geno else {
        return Err(Fail(format!("{}: no #CHROM line", vcf.display())));
    };
    if geno.sites() == 0 {
        return Err(Fail(format!("{} has no SNV", vcf.display())));
    }
    let dir = PathBuf::from(a.req("--out"));
    std::fs::create_dir_all(&dir).map_err(|e| crate::cannot_write(&dir, &e))?;
    let mut pending = Pending::new();
    pending.put(dir.join("ref.bim"), bim.as_bytes())?;
    pending.put(dir.join("ref.fam"), fam.as_bytes())?;
    let mut table = String::new();
    for ((level, g), text) in Level::ALL.iter().zip(&levels).zip(&p_rows) {
        pending.put(dir.join(format!("{}.P", level.name())), text.as_bytes())?;
        table.push_str(&format!(
            "{}\t{}\t{}\n",
            level.name(),
            g.codes.len(),
            g.codes.join(",")
        ));
    }
    let geno_path = dir.join("ref.geno");
    let shown = geno_path.clone();
    let w = pending.create(geno_path)?;
    geno.write(w).map_err(|e| crate::cannot_write(&shown, &e))?;
    pending.put(dir.join("levels.tsv"), table.as_bytes())?;
    pending.commit_all()?;
    writeln!(
        out,
        "{} SNVs of {} reference people; {}",
        geno.sites(),
        geno.people(),
        Level::ALL
            .iter()
            .zip(&levels)
            .map(|(l, g)| format!("{}: {} groups", l.name(), g.codes.len()))
            .collect::<Vec<_>>()
            .join(", ")
    )
    .map_err(|e| Fail(format!("cannot write stdout: {e}")))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testutil::{Dir, run};

    const SAMPLES: &str = "sample\tpopulation\tcontinent\tregion\treference\tpaint\n\
        a1\tYRI\tAFR\tWAF\t1\t1\n\
        e1\tCEU\tWEA\tNWE\t1\t0\n\
        a2\tYRI\tAFR\tWAF\t1\t0\n\
        s1\tIBS\tWEA\tIBE\t1\t0\n\
        x1\tASW\t-\t-\t0\t0\n\
        r1\tCEU\tWEA\tNWE\t0\t0\n\
        h1\tHan\tEAS\tEAS\t1\t0\n";

    // Columns: a1 e1 x1 a2 s1 r1 h1 (x1 has no group, r1 is no reference).
    const VCF: &str = "##fileformat=VCFv4.2\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta1\te1\tx1\ta2\ts1\tr1\th1\n\
        chr21\t100\t.\tA\tG\t.\t.\t.\tGT\t1|1\t0|1\t1|1\t1|0\t0|0\t1|1\t0|0\n\
        chr21\t200\t.\tC\tT\t.\t.\t.\tGT:DS\t0|0:0\t.|1:1\t0|0:0\t0|0:0\t1|1:2\t0|0:0\t1/1:2\n\
        chr22\t300\t.\tG\tC\t.\t.\t.\tGT\t0|1\t1|1\t0|0\t0|1\t1|0\t0|0\t0|0\n";

    fn fit(d: &Dir, extra: &[&str]) -> (u8, String, String) {
        let (vcf, samples, out) = (
            d.put("adm.vcf", VCF),
            d.put("samples.tsv", SAMPLES),
            d.path("model"),
        );
        let mut args = vec!["freqs", "--vcf", &vcf, "--samples", &samples, "--out", &out];
        args.extend_from_slice(extra);
        run(&args)
    }

    #[test]
    fn hand_counted_frequencies_of_three_groups() {
        let d = Dir::new("freqs-three");
        let (rc, out, err) = fit(&d, &[]);
        assert_eq!((rc, err.as_str()), (0, ""));
        assert_eq!(
            out,
            "3 SNVs of 5 reference people; continent: 3 groups, region: 4 groups\n"
        );
        // Groups in order of first appearance: AFR (a1 a2), WEA (e1 s1), EAS (h1).
        assert_eq!(
            d.read("model/levels.tsv"),
            "continent\t3\tAFR,WEA,EAS\nregion\t4\tWAF,NWE,IBE,EAS\n"
        );
        // REF (A2) frequencies, hand-counted:
        // 100 A>G: AFR 1/4 REF, WEA 3/4, EAS 2/2 (clamped).
        // 200 C>T: AFR 4/4 (clamped), WEA e1 missing so s1 only 0/2 (clamped), EAS 0/2.
        // 300 G>C: AFR 2/4, WEA 1/4, EAS 2/2.
        assert_eq!(
            d.read("model/continent.P"),
            "0.250000 0.750000 0.999990\n\
             0.999990 0.000010 0.000010\n\
             0.500000 0.250000 0.999990\n"
        );
        assert_eq!(
            d.read("model/region.P"),
            "0.250000 0.500000 0.999990 0.999990\n\
             0.999990 0.500000 0.000010 0.000010\n\
             0.500000 0.000010 0.500000 0.999990\n"
        );
        assert_eq!(
            d.read("model/ref.bim"),
            "21\t21:100\t0\t100\tG\tA\n21\t21:200\t0\t200\tT\tC\n22\t22:300\t0\t300\tC\tG\n"
        );
        assert_eq!(
            d.read("model/ref.fam"),
            "0\ta1\t0\t0\t0\t-9\n0\te1\t0\t0\t0\t-9\n0\ta2\t0\t0\t0\t-9\n0\ts1\t0\t0\t0\t-9\n0\th1\t0\t0\t0\t-9\n"
        );
        let geno = RefGeno::load(&d.0.join("model/ref.geno")).unwrap();
        assert_eq!((geno.sites(), geno.people()), (3, 5));
        // ref.geno counts A1 (ALT): a1 at 100 is 1|1, e1 at 200 is missing.
        assert_eq!(geno.get(0, 0), Some(2));
        assert_eq!(geno.get(1, 1), None);
        assert_eq!(geno.get(2, 3), Some(1));
        assert_eq!(geno.a1_count(0, &[0, 2]), (3, 2));
    }

    #[test]
    fn excluded_people_are_left_out() {
        let d = Dir::new("freqs-exclude");
        let ex = d.put("ex.txt", "# held out\ns1, h1\n");
        let (rc, out, err) = fit(&d, &["--exclude", &ex]);
        assert_eq!((rc, err.as_str()), (0, ""));
        assert!(
            out.starts_with("3 SNVs of 3 reference people; continent: 2 groups"),
            "{out}"
        );
        assert_eq!(
            d.read("model/levels.tsv"),
            "continent\t2\tAFR,WEA\nregion\t2\tWAF,NWE\n"
        );
        assert_eq!(
            d.read("model/continent.P").lines().next(),
            Some("0.250000 0.500000")
        );
    }

    #[test]
    fn a_bad_input_writes_nothing() {
        for (name, vcf, want) in [
            (
                "indel",
                VCF.replace("\tC\tT\t", "\tC\tCT\t"),
                "C>CT is not a biallelic SNV",
            ),
            (
                "haploid",
                VCF.replace("1|0\t0|0\t0|0\n", "1|0\t0|0\t1\n"),
                "GT 1 is not diploid",
            ),
            (
                "allele2",
                VCF.replace("1|0\t0|0\t0|0\n", "1|0\t0|0\t1|2\n"),
                "GT 1|2 has an allele other than 0 and 1",
            ),
            ("nogt", VCF.replace("GT:DS", "DS"), "no GT in FORMAT"),
            (
                "nobody",
                VCF.replace(
                    "\ta1\te1\tx1\ta2\ts1\tr1\th1\n",
                    "\tq1\tq2\tq3\tq4\tq5\tq6\tq7\n",
                ),
                "no reference person",
            ),
        ] {
            let d = Dir::new(&format!("freqs-bad-{name}"));
            d.put("adm.vcf", &vcf);
            let (vcf_p, samples, out) = (
                d.path("adm.vcf"),
                d.put("samples.tsv", SAMPLES),
                d.path("model"),
            );
            let (rc, _, err) = run(&[
                "freqs",
                "--vcf",
                &vcf_p,
                "--samples",
                &samples,
                "--out",
                &out,
            ]);
            assert_eq!(rc, 1, "{name}");
            assert!(err.contains(want), "{name}: {err}");
            assert!(!d.0.join("model/levels.tsv").exists(), "{name}");
        }
        let d = Dir::new("freqs-bad-header");
        let (vcf, samples, out) = (
            d.put("adm.vcf", VCF),
            d.put("samples.tsv", "id\tpop\n"),
            d.path("model"),
        );
        let (rc, _, err) = run(&["freqs", "--vcf", &vcf, "--samples", &samples, "--out", &out]);
        assert_eq!(rc, 1);
        assert!(
            err.contains("does not start with the header of samples.tsv"),
            "{err}"
        );
    }

    #[test]
    fn frequencies_are_clamped_as_admixture_does() {
        assert_eq!(a2_freq(0, 4), HIGH);
        assert_eq!(a2_freq(4, 4), LOW);
        assert_eq!(a2_freq(1, 4), 0.75);
        assert_eq!(a2_freq(0, 0), 0.5);
        assert_eq!(format!("{:.6} {:.6}", LOW, HIGH), "0.000010 0.999990");
    }
}
