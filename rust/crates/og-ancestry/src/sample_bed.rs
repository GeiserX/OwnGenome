//! `sample-bed`: ADMIXTURE's input for one person.
//!
//! A `.bed`/`.bim`/`.fam` of the reference `.bim`'s SNVs the person has a
//! genotype for, in the `.bim`'s order and allele coding, and each level's
//! frequencies (`.P.in`, which `admixture -P` reads) for the same SNVs.
//! ADMIXTURE stops on a SNV without a genotype, so a missing site, or a
//! genotype with an allele the `.bim` does not have, is left out. Also
//! `run_admixture.sh` ([`BOOT_SH`]), which the ADMIXTURE image runs (it has
//! sh and awk, no python): the projection at each level, then B bootstrap
//! replicates of it. ADMIXTURE's own `-B` crashes with `-P`, so a replicate
//! is a draw, with replacement, of up to 100 blocks of neighbouring SNVs,
//! written as a PLINK `.ped` (allele 2 is the `.bim`'s A1) and projected
//! the same way.
//!
//! The VCF streams: only the records at a `.bim` position are kept, so its
//! size does not matter. Every output is written under a temporary name and
//! renamed into place only after every input has been read and checked;
//! `run_admixture.sh` comes last, so a run that stops halfway leaves no
//! script for the next command to run on partial input.

use crate::{Args, Fail, Pending, Result, Spec};
use std::collections::HashMap;
use std::io::Write;
use std::path::{Path, PathBuf};

pub(crate) const SPEC: Spec = Spec {
    verb: "sample-bed",
    usage: "usage: og ancestry-fine sample-bed --bim ref.bim --vcf genotypes.vcf.gz --sample S --level LEVEL:K:P [--level ...] [--bootstrap B] [--min-sites N] --out-dir DIR\n",
    options: &[
        "--bim",
        "--vcf",
        "--sample",
        "--level",
        "--bootstrap",
        "--min-sites",
        "--out-dir",
    ],
    required: &["--bim", "--vcf", "--sample", "--level", "--out-dir"],
};

/// The bootstrap script; `@S@`, `@B@` and `@NB@` are the sample, the
/// number of replicates and the number of blocks per replicate.
pub const BOOT_SH: &str = include_str!("run_admixture.sh.in");

/// The reference `.bim`: its text in one buffer, and each SNV as offsets
/// into it, so a panel of any size costs a few dozen bytes per SNV.
pub struct Bim {
    text: String,
    rows: Vec<BimRow>,
}

/// One SNV: the line (with its ending), the chromosome without `chr`, the
/// alleles A1 and A2 as spans of [`Bim::text`], and the position.
#[derive(Debug, Clone, Copy)]
struct BimRow {
    line: (u32, u32),
    chrom: (u32, u32),
    a1: (u32, u32),
    a2: (u32, u32),
    pos: i64,
}

impl Bim {
    pub fn len(&self) -> usize {
        self.rows.len()
    }

    pub fn is_empty(&self) -> bool {
        self.rows.is_empty()
    }

    fn at(&self, (a, b): (u32, u32)) -> &str {
        &self.text[a as usize..b as usize]
    }

    /// Row `i`'s line as read, with its ending.
    pub fn line(&self, i: usize) -> &str {
        self.at(self.rows[i].line)
    }
}

/// plink2 writes chr22 as 22 in a `.bim`; the VCFs say chr22.
fn norm_chrom(c: &str) -> &str {
    c.strip_prefix("chr").unwrap_or(c)
}

pub fn read_bim(path: &Path) -> Result<Bim> {
    let mut bim = Bim {
        text: String::new(),
        rows: Vec::new(),
    };
    let mut n = 0u64;
    let too_big = || Fail(format!("{}: larger than 4 GB", path.display()));
    crate::for_each_raw_line(path, |line| {
        n += 1;
        let c: Vec<&str> = crate::split_ws(line).collect();
        if c.len() < 6 {
            return Ok(());
        }
        let pos = crate::int(c[3]).ok_or_else(|| {
            Fail(format!(
                "{}: line {n}: position is not a whole number: {}",
                path.display(),
                crate::repr(c[3])
            ))
        })?;
        let base = bim.text.len();
        let span = |s: &str| -> Result<(u32, u32)> {
            // Each field is a slice of `line`, so its offset is exact.
            let a = base + (s.as_ptr() as usize - line.as_ptr() as usize);
            let b = a + s.len();
            Ok((
                u32::try_from(a).map_err(|_| too_big())?,
                u32::try_from(b).map_err(|_| too_big())?,
            ))
        };
        let row = BimRow {
            line: span(line)?,
            chrom: span(norm_chrom(c[0]))?,
            a1: span(c[4])?,
            a2: span(c[5])?,
            pos,
        };
        bim.text.push_str(line);
        bim.rows.push(row);
        Ok(())
    })?;
    Ok(bim)
}

/// A GT of exactly two called alleles, phased or not, as written.
fn gt_alleles(gt: &str) -> Option<[&str; 2]> {
    let mut parts = gt.split(['/', '|']);
    let (a, b) = (parts.next()?, parts.next()?);
    if parts.next().is_some() || a == "." || b == "." {
        return None;
    }
    Some([a, b])
}

/// The PLINK `.bed` code of a genotype (REF and ALT alleles, the two allele
/// indices as written) for a `.bim` row with alleles A1, A2: 0 homozygous
/// A1, 2 heterozygous, 3 homozygous A2; `None` for no call (an index that
/// is no allele, or an allele the `.bim` does not have).
pub fn bed_code(alleles: &[&str], idx: [&str; 2], a1: &str, a2: &str) -> Option<u8> {
    let n = alleles.len() as i64;
    let (mut n1, mut n2) = (0, 0);
    for i in idx {
        let i = crate::int(i)?;
        // A negative index counts from the end, as Python's list indexing.
        let at = if i < 0 { n + i } else { i };
        let called = *alleles.get(usize::try_from(at).ok()?)?;
        n1 += usize::from(called == a1);
        n2 += usize::from(called == a2);
    }
    if n1 + n2 != 2 {
        return None;
    }
    Some(match n1 {
        2 => 0,
        1 => 2,
        _ => 3,
    })
}

/// The `.bed` code of each `.bim` row, streamed from the VCF: the last
/// record with a called GT at a row's position decides it. Only one byte
/// per row is kept, whatever the VCF's size.
pub fn codes_at(bim: &Bim, vcf: &Path) -> Result<Vec<Option<u8>>> {
    // Every row by (chromosome, position), sorted, so a record finds its
    // rows by binary search (a .bim may list two SNVs at one position).
    // Chromosome names are numbered in order of first appearance.
    let mut names: HashMap<&str, u32> = HashMap::new();
    let mut index: Vec<(u32, i64, u32)> = Vec::with_capacity(bim.len());
    for (i, r) in bim.rows.iter().enumerate() {
        let next = u32::try_from(names.len()).expect("fewer chromosomes than rows");
        let c = *names.entry(bim.at(r.chrom)).or_insert(next);
        let i = u32::try_from(i).map_err(|_| Fail("more than 4 billion SNVs".into()))?;
        index.push((c, r.pos, i));
    }
    index.sort_unstable();
    let mut codes: Vec<Option<u8>> = vec![None; bim.len()];
    let mut n = 0u64;
    crate::for_each_line(vcf, |line| {
        n += 1;
        if line.starts_with('#') {
            return Ok(());
        }
        let c: Vec<&str> = line.split('\t').collect();
        if c.len() < 10 {
            return Ok(());
        }
        let Some(gt_at) = c[8].split(':').position(|f| f == "GT") else {
            return Ok(());
        };
        let Some(idx) = c[9].split(':').nth(gt_at).and_then(gt_alleles) else {
            return Ok(());
        };
        let pos = crate::int(c[1]).ok_or_else(|| {
            Fail(format!(
                "{}: line {n}: POS is not a whole number: {}",
                vcf.display(),
                crate::repr(c[1])
            ))
        })?;
        let Some(&chrom) = names.get(norm_chrom(c[0])) else {
            return Ok(());
        };
        let from = index.partition_point(|&(ch, p, _)| (ch, p) < (chrom, pos));
        let rows = index[from..]
            .iter()
            .take_while(|&&(ch, p, _)| (ch, p) == (chrom, pos));
        let mut alleles = vec![c[3]];
        if c[4] != "." {
            alleles.extend(c[4].split(','));
        }
        for &(_, _, i) in rows {
            let r = &bim.rows[i as usize];
            codes[i as usize] = bed_code(&alleles, idx, bim.at(r.a1), bim.at(r.a2));
        }
        Ok(())
    })?;
    Ok(codes)
}

/// The called SNVs: (row indices, `.bed` codes).
pub fn kept(codes: &[Option<u8>]) -> (Vec<usize>, Vec<u8>) {
    codes
        .iter()
        .enumerate()
        .filter_map(|(i, c)| c.map(|c| (i, c)))
        .unzip()
}

/// The `.bed` file: magic, SNP-major, then one byte per SNV for one person.
pub fn bed_bytes(codes: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(3 + codes.len());
    out.extend_from_slice(&[0x6c, 0x1b, 0x01]);
    out.extend_from_slice(codes);
    out
}

/// The allele pair of a `.bed` code in a PLINK "12" `.ped`: allele 2 is the
/// `.bim`'s A1 (ADMIXTURE's projection of the `.ped` equals the `.bed`'s
/// only this way round).
fn ped_alleles(code: u8) -> &'static str {
    match code {
        0 => "2 2",
        2 => "1 2",
        _ => "1 1",
    }
}

/// The script for sample `s`, `b` replicates of `nb` blocks.
pub fn boot_script(s: &str, b: i64, nb: usize) -> String {
    BOOT_SH
        .replace("@S@", s)
        .replace("@B@", &b.to_string())
        .replace("@NB@", &nb.to_string())
}

pub(crate) fn main(a: &Args, out: &mut dyn Write) -> Result<()> {
    let bim_arg = a.req("--bim");
    let vcf_arg = a.req("--vcf");
    let sample = a.req("--sample");
    let bootstrap = a.int("--bootstrap", 100).unwrap_or(100);
    let min_sites = a.int("--min-sites", 100).unwrap_or(100);
    let bim = read_bim(Path::new(bim_arg))?;
    if bim.is_empty() {
        return Err(Fail(format!("{bim_arg} has no SNV")));
    }
    let (keep, codes) = kept(&codes_at(&bim, Path::new(vcf_arg))?);
    if (keep.len() as i64) < min_sites {
        return Err(Fail(format!(
            "only {} of the {} reference SNVs have a genotype in {vcf_arg} (at least {min_sites} needed)",
            keep.len(),
            bim.len()
        )));
    }
    let dir = PathBuf::from(a.req("--out-dir"));
    std::fs::create_dir_all(&dir).map_err(|e| crate::cannot_write(&dir, &e))?;
    let at = |name: String| dir.join(name);
    let bed = bed_bytes(&codes);
    let mut pending = Pending::new();
    let mut levels = String::new();
    for spec in a.all("--level") {
        let mut parts = spec.splitn(3, ':');
        let (Some(level), Some(k), Some(p_path)) = (parts.next(), parts.next(), parts.next())
        else {
            return Err(Fail(format!("--level {spec}: not LEVEL:K:P_FILE")));
        };
        let base = format!("{sample}_{level}");
        // The .P rows of the kept SNVs, streamed into the output; the row
        // count is checked before anything is renamed into place.
        let p_out = at(format!("{base}.{k}.P.in"));
        let shown = p_out.clone();
        let w = pending.create(p_out)?;
        let (mut rows, mut next) = (0usize, keep.iter().peekable());
        crate::for_each_raw_line(Path::new(p_path), |line| {
            if next.peek() == Some(&&rows) {
                next.next();
                w.write_all(line.as_bytes())
                    .map_err(|e| crate::cannot_write(&shown, &e))?;
            }
            rows += 1;
            Ok(())
        })?;
        if rows != bim.len() {
            return Err(Fail(format!(
                "{p_path} has {rows} rows but {bim_arg} {} SNVs",
                bim.len()
            )));
        }
        pending.put(at(format!("{base}.bed")), &bed)?;
        let bim_out = at(format!("{base}.bim"));
        let shown = bim_out.clone();
        let w = pending.create(bim_out)?;
        for &i in &keep {
            w.write_all(bim.line(i).as_bytes())
                .map_err(|e| crate::cannot_write(&shown, &e))?;
        }
        pending.put(
            at(format!("{base}.fam")),
            format!("{sample} {sample} 0 0 0 -9\n").as_bytes(),
        )?;
        levels.push_str(&format!("{level} {k}\n"));
    }
    let geno_out = at(format!("{sample}.geno"));
    let shown = geno_out.clone();
    let w = pending.create(geno_out)?;
    for &c in &codes {
        writeln!(w, "{}", ped_alleles(c)).map_err(|e| crate::cannot_write(&shown, &e))?;
    }
    let map_out = at(format!("{sample}.map"));
    let shown = map_out.clone();
    let w = pending.create(map_out)?;
    for &i in &keep {
        let first4: Vec<&str> = crate::split_ws(bim.line(i)).take(4).collect();
        writeln!(w, "{}", first4.join(" ")).map_err(|e| crate::cannot_write(&shown, &e))?;
    }
    pending.put(at("levels.txt".to_owned()), levels.as_bytes())?;
    let sites = format!(
        "{} of the {} SNVs ADMIXTURE uses genotyped",
        keep.len(),
        bim.len()
    );
    pending.put(
        at(format!("{sample}_sites.txt")),
        format!("{sites}\n").as_bytes(),
    )?;
    let nb = (keep.len() / 5).clamp(1, 100);
    pending.put(
        at("run_admixture.sh".to_owned()),
        boot_script(sample, bootstrap, nb).as_bytes(),
    )?;
    pending.commit_all()?;
    writeln!(out, "{sites}").map_err(|e| Fail(format!("cannot write stdout: {e}")))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testutil::{Dir, run};

    const BIM: &str = "22\t22:100\t0\t100\tG\tA\n22\t22:200\t0\t200\tT\tC\n22\t22:300\t0\t300\tA\tG\n22\t22:400\t0\t400\tC\tT\n22\t22:500\t0\t500\tC\tA\n";
    const VCF: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n\
        chr22\t100\t.\tA\tG\t.\t.\t.\tGT\t1/1\n\
        chr22\t200\t.\tC\tT\t.\t.\t.\tGT\t0|1\n\
        chr22\t300\t.\tG\tA\t.\t.\t.\tGT\t0/0\n\
        chr22\t400\t.\tT\tG\t.\t.\t.\tGT\t0/1\n";

    #[test]
    fn codes_of_each_genotype() {
        let d = Dir::new("bed-codes");
        let bim = read_bim(Path::new(&d.put("ref.bim", BIM))).unwrap();
        let (keep, codes) = kept(&codes_at(&bim, Path::new(&d.put("g.vcf", VCF))).unwrap());
        assert_eq!(
            keep,
            [0, 1, 2],
            "chr22 meets 22; another ALT and a missing site are left out"
        );
        assert_eq!(bed_bytes(&codes), b"\x6c\x1b\x01\x00\x02\x03");
    }

    #[test]
    fn a_later_called_record_wins_and_an_uncalled_one_does_not() {
        let d = Dir::new("bed-later");
        let bim = read_bim(Path::new(&d.put("ref.bim", "22\ta\t0\t100\tG\tA\n"))).unwrap();
        let v = "chr22\t100\t.\tA\tG\t.\t.\t.\tGT\t0/0\nchr22\t100\t.\tA\tG\t.\t.\t.\tGT\t1/1\nchr22\t100\t.\tA\tG\t.\t.\t.\tGT\t./.\n";
        let codes = codes_at(&bim, Path::new(&d.put("g.vcf", v))).unwrap();
        assert_eq!(codes, [Some(0)]);
    }

    #[test]
    fn python_list_indexing_of_alleles() {
        assert_eq!(bed_code(&["A", "G"], ["0", "-1"], "G", "A"), Some(2));
        assert_eq!(bed_code(&["A"], ["0", "2"], "G", "A"), None);
        assert_eq!(bed_code(&["A", "G"], ["x", "1"], "G", "A"), None);
        assert_eq!(gt_alleles("0/1/1"), None);
        assert_eq!(gt_alleles("./1"), None);
        assert_eq!(gt_alleles("1"), None);
    }

    #[test]
    fn two_snvs_at_one_position_each_get_their_code() {
        let d = Dir::new("bed-twice");
        let bim = read_bim(Path::new(&d.put(
            "ref.bim",
            "22	a	0	100	G	A
22	b	0	100	T	A
",
        )))
        .unwrap();
        let v = "chr22\t100\t.\tA\tG,T\t.\t.\t.\tGT\t0/2\n";
        assert_eq!(
            codes_at(&bim, Path::new(&d.put("g.vcf", v))).unwrap(),
            [None, Some(2)]
        );
        assert_eq!(bim.line(1), "22\tb\t0\t100\tT\tA\n");
    }

    #[test]
    fn the_command_writes_every_input_of_admixture() {
        let d = Dir::new("bed-cmd");
        let bim = d.put("ref.bim", BIM);
        let vcf = d.put("g.vcf", VCF);
        let p = d.put(
            "continent.P",
            "0.0 0.9\n0.1 0.8\n0.2 0.7\n0.3 0.6\n0.4 0.5\n",
        );
        let out = d.path("adm");
        let level = format!("continent:2:{p}");
        let (rc, stdout, err) = run(&[
            "sample-bed",
            "--bim",
            &bim,
            "--vcf",
            &vcf,
            "--sample",
            "S",
            "--level",
            &level,
            "--bootstrap",
            "7",
            "--min-sites",
            "3",
            "--out-dir",
            &out,
        ]);
        assert_eq!((rc, err.as_str()), (0, ""));
        assert_eq!(stdout, "3 of the 5 SNVs ADMIXTURE uses genotyped\n");
        assert_eq!(
            d.read("adm/S_continent.2.P.in"),
            "0.0 0.9\n0.1 0.8\n0.2 0.7\n"
        );
        assert_eq!(d.read("adm/S_continent.bim").lines().count(), 3);
        assert_eq!(d.read("adm/S.geno"), "2 2\n1 2\n1 1\n");
        assert_eq!(
            d.read("adm/S.map"),
            "22 22:100 0 100\n22 22:200 0 200\n22 22:300 0 300\n"
        );
        assert_eq!(d.read("adm/levels.txt"), "continent 2\n");
        let sh = d.read("adm/run_admixture.sh");
        assert!(sh.contains("-le 7 ]") && sh.contains("S_${level}.bed") && sh.contains("nb=1 "));
        // Too few genotyped SNVs, and a .P of the wrong length: nothing new
        // is written, and what was there stays.
        let fresh = d.path("fresh");
        let (rc, _, err) = run(&[
            "sample-bed",
            "--bim",
            &bim,
            "--vcf",
            &vcf,
            "--sample",
            "S",
            "--level",
            &level,
            "--out-dir",
            &fresh,
        ]);
        assert_eq!(rc, 1);
        assert!(
            err.contains("only 3 of the 5 reference SNVs have a genotype in"),
            "{err}"
        );
        let short = d.put("short.P", "0.1 0.9\n");
        let (rc, _, err) = run(&[
            "sample-bed",
            "--bim",
            &bim,
            "--vcf",
            &vcf,
            "--sample",
            "S",
            "--level",
            &level,
            "--level",
            &format!("region:2:{short}"),
            "--min-sites",
            "3",
            "--out-dir",
            &fresh,
        ]);
        assert_eq!(rc, 1);
        assert!(
            err.ends_with(&format!("{short} has 1 rows but {bim} 5 SNVs\n")),
            "{err}"
        );
        assert_eq!(
            std::fs::read_dir(&fresh).unwrap().count(),
            0,
            "a failure leaves no partial input"
        );
        let (rc, _, err) = run(&[
            "sample-bed",
            "--bim",
            &bim,
            "--vcf",
            &vcf,
            "--sample",
            "S",
            "--level",
            "continent",
            "--min-sites",
            "3",
            "--out-dir",
            &fresh,
        ]);
        assert_eq!(
            (rc, err.as_str()),
            (1, "ERROR: --level continent: not LEVEL:K:P_FILE\n")
        );
    }

    #[test]
    fn the_bootstrap_draws_only_blocks_that_hold_snvs() {
        // 4 SNVs in 3 blocks round up to blocks of 2, so there are 2 blocks
        // and every replicate holds 3 x 2 SNVs.
        let awk = BOOT_SH
            .split_once("-v k=\"$k\" '")
            .and_then(|(_, rest)| rest.split_once("' \"@S@.geno\""))
            .map(|(prog, _)| prog)
            .unwrap();
        let d = Dir::new("bed-boot");
        d.put("g", "1 0\n1 1\n1 2\n1 3\n");
        d.put("m", "22 s0 0 0\n22 s1 0 1\n22 s2 0 2\n22 s3 0 3\n");
        d.put("p", "0.0 0.5\n0.1 0.5\n0.2 0.5\n0.3 0.5\n");
        for r in 1..=20 {
            let ok = std::process::Command::new("awk")
                .args([
                    "-v",
                    &format!("r={r}"),
                    "-v",
                    "nb=3",
                    "-v",
                    "k=2",
                    awk,
                    "g",
                    "m",
                    "p",
                ])
                .current_dir(&d.0)
                .status()
                .unwrap();
            assert!(ok.success());
            assert_eq!(d.read("rep.map").lines().count(), 6, "seed {r}");
        }
    }
}
