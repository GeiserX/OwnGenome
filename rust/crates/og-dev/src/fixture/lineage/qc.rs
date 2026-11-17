//! The sample-QC cases' data cut to the fixture's contigs, and the
//! synthetic woman of the mosaic-loss case.

use super::{Fail, PyRandom, bad, cell_int, create, each_line, lines, want};
use flate2::Compression;
use flate2::write::GzEncoder;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Write};

/// The contig names of a `.fai`, in its order.
fn contig_order(fai: &str) -> Result<HashMap<String, usize>, Fail> {
    let mut order = HashMap::new();
    let mut i = 0;
    each_line(fai, |l| {
        order.insert(l.split('\t').next().unwrap_or("").to_string(), i);
        i += 1;
        Ok(())
    })?;
    Ok(order)
}

/// `og dev fixture contig-cut FAI SITES.vcf.gz OUT_SITES PANEL OUT_PANEL`:
/// somalier's sites and VerifyBamID2's panel (PANEL.bed, .mu and .UD, one
/// marker a line in each) cut to the contigs of FAI, since somalier stops on
/// a site whose contig the FASTA lacks. Prints the markers kept.
pub fn contig_cut(args: &[String], out: &mut dyn Write, _err: &mut dyn Write) -> Result<(), Fail> {
    want(args, 5)?;
    let contigs: HashSet<String> = contig_order(&args[0])?.into_keys().collect();
    let first = |l: &str| l.split('\t').next().unwrap_or("").to_string();
    let mut o = create(&args[2])?;
    each_line(&args[1], |l| {
        if l.starts_with('#') || contigs.contains(&first(l)) {
            writeln!(o, "{l}")?;
        }
        Ok(())
    })?;
    o.flush()?;
    let (panel, fix) = (&args[3], &args[4]);
    let bed = lines(&format!("{panel}.bed"))?;
    let keep: Vec<usize> = (0..bed.len())
        .filter(|&i| contigs.contains(&first(&bed[i])))
        .collect();
    for ext in ["bed", "mu", "UD"] {
        let all = lines(&format!("{panel}.{ext}"))?;
        if all.len() != bed.len() {
            return Err(bad(format!(
                "{panel}.{ext} has {} lines, the .bed {}",
                all.len(),
                bed.len()
            )));
        }
        let mut o = create(&format!("{fix}.{ext}"))?;
        for &i in &keep {
            writeln!(o, "{}", all[i])?;
        }
        o.flush()?;
    }
    writeln!(out, "{} panel markers on the fixture's contigs", keep.len())?;
    Ok(())
}

/// `og dev fixture add-sites SITES.vcf EXTRA.tsv FAI OUT.vcf`: the sites VCF
/// (its first 8 columns) with the rows of EXTRA added where SITES has no
/// record at that contig and position, sorted by FAI's contig order and
/// position. Prints the site count.
pub fn add_sites(args: &[String], out: &mut dyn Write, _err: &mut dyn Write) -> Result<(), Fail> {
    want(args, 4)?;
    let order = contig_order(&args[2])?;
    let (mut head, mut cols, mut recs) = (Vec::new(), None, Vec::<Vec<String>>::new());
    each_line(&args[0], |l| {
        if l.starts_with("##") {
            head.push(l.to_string());
        } else if l.starts_with('#') {
            cols = Some(l.to_string());
        } else {
            recs.push(l.split('\t').take(8).map(str::to_string).collect());
        }
        Ok(())
    })?;
    let cols = cols.ok_or_else(|| bad(format!("{}: no #CHROM line", args[0])))?;
    let seen: HashSet<(String, String)> = recs
        .iter()
        .map(|r| (r[0].clone(), r.get(1).cloned().unwrap_or_default()))
        .collect();
    each_line(&args[1], |l| {
        let r: Vec<String> = l.split('\t').map(str::to_string).collect();
        if !seen.contains(&(r[0].clone(), r.get(1).cloned().unwrap_or_default())) {
            recs.push(r);
        }
        Ok(())
    })?;
    let mut keyed = recs
        .into_iter()
        .map(|r| {
            let pos = cell_int(r.get(1).map_or("", String::as_str))?;
            Ok(((order.get(&r[0]).copied().unwrap_or(1 << 30), pos), r))
        })
        .collect::<Result<Vec<_>, Fail>>()?;
    keyed.sort_by_key(|x| x.0);
    let mut o = create(&args[3])?;
    for h in &head {
        writeln!(o, "{h}")?;
    }
    writeln!(o, "{cols}")?;
    for (_, r) in &keyed {
        writeln!(o, "{}", r.join("\t"))?;
    }
    o.flush()?;
    writeln!(out, "{} sites", keyed.len())?;
    Ok(())
}

/// A gzip file to write.
fn gz(path: &str) -> Result<GzEncoder<BufWriter<File>>, Fail> {
    Ok(GzEncoder::new(create(path)?, Compression::default()))
}

/// `og dev fixture sex-loss-woman WINDOWS OUT_DIR NAME LOSS`: a woman who
/// lost one X in LOSS of her cells (seed 11), from a man's mosdepth windows:
/// `OUT_DIR/mosdepth/NAME.regions.bed.gz` (no chrY, chrX depth doubled less
/// the loss), `OUT_DIR/vcf/NAME.vcf.gz` (2,000 heterozygous chr20 sites at
/// half and 300 chrX sites shifted by the loss) and the sample-QC table
/// `og qc sex-loss` adds to.
pub fn sex_loss_woman(
    args: &[String],
    _out: &mut dyn Write,
    _err: &mut dyn Write,
) -> Result<(), Fail> {
    want(args, 4)?;
    let (win, dir, name) = (&args[0], &args[1], &args[2]);
    let loss: f64 = args[3]
        .trim()
        .parse()
        .map_err(|_| Fail::Usage(format!("LOSS '{}' is not a number", args[3])))?;
    let d = loss / (2.0 * (2.0 - loss));
    let mut rnd = PyRandom::new(11);
    let mut o = gz(&format!("{dir}/mosdepth/{name}.regions.bed.gz"))?;
    each_line(win, |l| {
        let c: Vec<&str> = l.split('\t').collect();
        let [chrom, s, e, dp] = c[..] else {
            return Err(bad(format!("{win}: a window without 4 columns")));
        };
        match chrom {
            "chrY" => {}
            "chrX" => {
                let x: f64 = dp
                    .trim()
                    .parse()
                    .map_err(|_| bad(format!("{win}: depth '{dp}'")))?;
                let dp = pycompat::fmt::fixed(x * 2.0 * (1.0 - loss / 2.0), 2);
                writeln!(o, "{chrom}\t{s}\t{e}\t{dp}")?;
            }
            _ => writeln!(o, "{chrom}\t{s}\t{e}\t{dp}")?,
        }
        Ok(())
    })?;
    o.finish()?.flush()?;
    let mut o = gz(&format!("{dir}/vcf/{name}.vcf.gz"))?;
    write!(
        o,
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTH\n"
    )?;
    let site = |o: &mut dyn Write,
                rnd: &mut PyRandom,
                chrom: &str,
                pos: i64,
                p: f64|
     -> Result<(), Fail> {
        let n = rnd.randint(20, 40);
        let a = (0..n).filter(|_| rnd.random() < p).count() as i64;
        writeln!(
            o,
            "{chrom}\t{pos}\t.\tA\tG\t40\tPASS\t.\tGT:AD\t0/1:{},{a}",
            n - a
        )?;
        Ok(())
    };
    for i in 0..2000 {
        site(&mut o, &mut rnd, "chr20", 10_000_000 + 250 * i, 0.5)?;
    }
    for i in 0..300 {
        let p = 0.5 + if rnd.random() < 0.5 { d } else { -d };
        site(&mut o, &mut rnd, "chrX", 73_700_000 + 1000 * i, p)?;
    }
    o.finish()?.flush()?;
    let mut o = create(&format!("{dir}/qc/{name}_sample_qc.tsv"))?;
    write!(
        o,
        "key\tvalue\nsample\t{name}\ninferred_sex\tunknown\ndeclared_sex\tfemale\n"
    )?;
    o.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::super::testutil::{dir, fnv, run};
    use std::io::Read;

    fn gunzip(p: &std::path::Path) -> Vec<u8> {
        let mut v = Vec::new();
        flate2::read::GzDecoder::new(std::fs::File::open(p).expect("open"))
            .read_to_end(&mut v)
            .expect("gunzip");
        v
    }

    #[test]
    fn contig_cut_keeps_the_fai_contigs_in_step() {
        let d = dir("contig-cut");
        let w = |n: &str, s: &str| std::fs::write(d.join(n), s).expect("write");
        w("r.fai", "chr1\t10\t1\t60\t61\nchr20\t10\t20\t60\t61\n");
        w(
            "sites.vcf",
            "##h\n#CHROM\tPOS\nchr1\t5\nchr11\t6\nchr20\t7\n",
        );
        w("p.bed", "chr11\t1\t2\nchr1\t3\t4\nchr20\t5\t6\n");
        w("p.mu", "a\nb\nc\n");
        w("p.UD", "x y\nz w\nu v\n");
        let p = |n: &str| d.join(n).to_str().expect("utf8").to_string();
        let (code, out, err) = run(
            "contig-cut",
            &[&p("r.fai"), &p("sites.vcf"), &p("o.vcf"), &p("p"), &p("f")],
        );
        assert_eq!(
            (code, out.as_str()),
            (0, "2 panel markers on the fixture's contigs\n"),
            "{err}"
        );
        let r = |n: &str| std::fs::read_to_string(d.join(n)).expect("read");
        assert_eq!(r("o.vcf"), "##h\n#CHROM\tPOS\nchr1\t5\nchr20\t7\n");
        assert_eq!(
            (r("f.bed").as_str(), r("f.mu").as_str(), r("f.UD").as_str()),
            ("chr1\t3\t4\nchr20\t5\t6\n", "b\nc\n", "z w\nu v\n")
        );
        w("p.mu", "a\nb\n");
        let (code, _, _) = run(
            "contig-cut",
            &[&p("r.fai"), &p("sites.vcf"), &p("o.vcf"), &p("p"), &p("f")],
        );
        assert_eq!(code, 1);
    }

    #[test]
    fn add_sites_adds_new_positions_in_fai_order() {
        let d = dir("add-sites");
        let w = |n: &str, s: &str| std::fs::write(d.join(n), s).expect("write");
        w("r.fai", "chr2\t10\nchrX\t10\n");
        w(
            "s.vcf",
            "##a\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\nchrX\t50\t.\tA\tG\t.\tPASS\tAF=0.5\tGT\nchr2\t9\t.\tC\tT\t.\tPASS\tAF=0.5\n",
        );
        w(
            "x.tsv",
            "chrX\t50\t.\tA\tC\t.\tPASS\tAF=0.5\nchrX\t20\t.\tG\tT\t.\tPASS\tAF=0.5\nchr7\t1\t.\tG\tT\t.\tPASS\tAF=0.5\n",
        );
        let p = |n: &str| d.join(n).to_str().expect("utf8").to_string();
        let (code, out, err) = run(
            "add-sites",
            &[&p("s.vcf"), &p("x.tsv"), &p("r.fai"), &p("o.vcf")],
        );
        assert_eq!((code, out.as_str()), (0, "4 sites\n"), "{err}");
        assert_eq!(
            std::fs::read_to_string(d.join("o.vcf")).expect("read"),
            "##a\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\nchr2\t9\t.\tC\tT\t.\tPASS\tAF=0.5\nchrX\t20\t.\tG\tT\t.\tPASS\tAF=0.5\nchrX\t50\t.\tA\tG\t.\tPASS\tAF=0.5\nchr7\t1\t.\tG\tT\t.\tPASS\tAF=0.5\n"
        );
    }

    #[test]
    fn sex_loss_woman_is_the_pythons_content() {
        let d = dir("woman");
        for s in ["mosdepth", "vcf", "qc"] {
            std::fs::create_dir_all(d.join("o").join(s)).expect("dir");
        }
        std::fs::write(d.join("w.bed"), "chr20\t0\t500\t30.25\nchrX\t0\t500\t15.125\nchrY\t0\t500\t14.00\nchrX\t500\t1000\t14.375\n").expect("bed");
        let p = |n: &str| d.join(n).to_str().expect("utf8").to_string();
        let (code, _, err) = run("sex-loss-woman", &[&p("w.bed"), &p("o"), "W", "0.3"]);
        assert_eq!(code, 0, "{err}");
        let o = d.join("o");
        assert_eq!(
            String::from_utf8(gunzip(&o.join("mosdepth/W.regions.bed.gz"))).expect("utf8"),
            "chr20\t0\t500\t30.25\nchrX\t0\t500\t25.71\nchrX\t500\t1000\t24.44\n"
        );
        // FNV-1a of the VCF the Python heredoc of mosaic-loy-lox-1.sh wrote.
        let vcf = gunzip(&o.join("vcf/W.vcf.gz"));
        assert_eq!(fnv(&vcf), 0xc45b_5dc4_9d43_4229, "{:#x}", fnv(&vcf));
        assert_eq!(
            std::fs::read_to_string(o.join("qc/W_sample_qc.tsv")).expect("read"),
            "key\tvalue\nsample\tW\ninferred_sex\tunknown\ndeclared_sex\tfemale\n"
        );
    }
}
