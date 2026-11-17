//! A synthetic AADR-like panel in the AADR's own formats, for step 44.

use super::{Fail, Fasta, PyRandom, bad, cell_int, create, each_line, read_bed, read_fai, want};
use std::collections::HashMap;
use std::io::Write;

/// The made-up GRCh37 coordinates sit this far below the GRCh38 ones.
const OFF: i64 = 1000;

/// One person of the panel.
struct Person {
    id: String,
    sex: &'static str,
    group: &'static str,
    modern: bool,
    bp: u32,
    country: &'static str,
    lat: u32,
    lon: u32,
    /// Dosage of the `.snp`'s first allele (REF) per site; 9 is missing.
    dosage: Vec<u8>,
}

/// Python's `int(c[3:])` when `c[3:].isdigit()`: the number of `chr<N>`.
fn autosome(c: &str) -> Option<i64> {
    let d = c.get(3..)?;
    if !d.is_empty() && d.bytes().all(|b| b.is_ascii_digit()) {
        d.parse().ok()
    } else {
        None
    }
}

/// Pseudo-haploid dosages: 30% missing, else REF with probability
/// `1 - alt / 2`.
fn pseudo(rnd: &mut PyRandom, alts: &[f64]) -> Vec<u8> {
    alts.iter()
        .map(|&a| {
            if rnd.random() < 0.3 {
                9
            } else if rnd.random() < a / 2.0 {
                0
            } else {
                2
            }
        })
        .collect()
}

/// `og dev fixture ancient-panel REGIONS TRUTH CONFIDENT FASTA OUT_DIR CHAIN`:
/// the AADR-like panel of `tests/e2e/ancient-similarity-1.sh` (seed 44).
/// Sites are the truth set's SNVs in the autosomal slices of REGIONS, plus
/// a hom-ref site every 150 bases, inside the confident regions. Three
/// present-day groups of 30 (French and Finnish from random frequencies,
/// Sardinian from HG002's genotypes) and four ancient ones: six
/// pseudo-haploid copies each of HG002 (Planted_HG002) and of the French
/// and Finnish frequencies, and one person outside West Eurasia. Writes
/// `aadr.anno`, `.ind`, `.snp`, the transposed packed `.geno` (TGENO) and
/// `expected_sites.tsv` in OUT_DIR, and a chain moving the made-up GRCh37
/// slices back to GRCh38.
pub fn panel(args: &[String], out: &mut dyn Write, _err: &mut dyn Write) -> Result<(), Fail> {
    want(args, 6)?;
    let (regions, truth, bed, reference, dir, chain) =
        (&args[0], &args[1], &args[2], &args[3], &args[4], &args[5]);
    let fai = read_fai(&format!("{reference}.fai"))?;
    let mut fa = Fasta::open(reference)?;
    let mut slices = Vec::new();
    each_line(regions, |l| {
        let c: Vec<&str> = l.split('\t').collect();
        if c.len() < 3 {
            return Err(bad(format!("{regions}: a line with fewer than 3 columns")));
        }
        if fai.contains_key(c[0]) && autosome(c[0]).is_some_and(|n| n <= 22) {
            let s = cell_int(c[1])?;
            if s > OFF + 10 {
                slices.push((c[0].to_string(), s, cell_int(c[2])?));
            }
        }
        Ok(())
    })?;
    let conf = read_bed(bed)?;
    let confident = |c: &str, p: i64| conf.get(c).is_some_and(|b| b.has(p));
    let mut truth_gt: HashMap<(String, i64), (u8, u8, Option<u8>)> = HashMap::new();
    each_line(truth, |l| {
        if l.starts_with('#') {
            return Ok(());
        }
        let f: Vec<&str> = l.split('\t').collect();
        if f.len() < 10 {
            return Err(bad(format!("{truth}: a record with fewer than 10 columns")));
        }
        if f[3].len() == 1 && f[4].len() == 1 {
            let g = f[9].split(':').next().unwrap_or("").replace('|', "/");
            let n = match g.as_str() {
                "0/1" | "1/0" => Some(1),
                "1/1" => Some(2),
                _ => None,
            };
            truth_gt.insert(
                (f[0].to_string(), cell_int(f[1])?),
                (f[3].as_bytes()[0], f[4].as_bytes()[0], n),
            );
        }
        Ok(())
    })?;
    // (contig, GRCh38 position, ref, alt, HG002's ALT count)
    let mut sites: Vec<(String, i64, u8, u8, u8)> = Vec::new();
    for (c, s, e) in &slices {
        for p in s + 1..=*e {
            let t = truth_gt.get(&(c.clone(), p));
            if let Some(&(r, a, Some(n))) = t
                && confident(c, p)
                && fa.base(c, p as u64)? == Some(r)
            {
                sites.push((c.clone(), p, r, a, n));
            } else if p % 150 == 0 && t.is_none() && confident(c, p) {
                let alt = match fa.base(c, p as u64)? {
                    Some(b'A') => b'G',
                    Some(b'C') => b'T',
                    Some(b'G') => b'A',
                    Some(b'T') => b'C',
                    _ => continue,
                };
                let r = fa.base(c, p as u64)?.unwrap_or(b'N');
                sites.push((c.clone(), p, r, alt, 0));
            }
        }
    }
    sites.sort_by_key(|x| (autosome(&x.0).unwrap_or(0), x.1));
    let n = sites.len();
    if n < 1000 {
        return Err(bad(format!(
            "only {n} sites in the slices, fewer than 1000"
        )));
    }
    let mut rnd = PyRandom::new(44);
    let french: Vec<f64> = (0..n).map(|_| rnd.uniform(0.05, 0.95)).collect();
    let finnish: Vec<f64> = (0..n).map(|_| rnd.uniform(0.05, 0.95)).collect();
    let sardinian: Vec<f64> = sites
        .iter()
        .map(|x| [0.1, 0.5, 0.9][usize::from(x.4)])
        .collect();
    let mut people = Vec::new();
    for (g, f) in [
        ("French", &french),
        ("Finnish", &finnish),
        ("Sardinian", &sardinian),
    ] {
        for i in 0..30 {
            let dosage = f
                .iter()
                .map(|&q| {
                    let k = u8::from(rnd.random() < q) + u8::from(rnd.random() < q);
                    2 - k
                })
                .collect();
            people.push(Person {
                id: format!("{g}{i}.DG"),
                sex: "M",
                group: g,
                modern: true,
                bp: 0,
                country: "Europe",
                lat: 45,
                lon: 5,
                dosage,
            });
        }
    }
    let hg002: Vec<f64> = sites.iter().map(|x| f64::from(x.4)).collect();
    let fr2: Vec<f64> = french.iter().map(|q| 2.0 * q).collect();
    let fi2: Vec<f64> = finnish.iter().map(|q| 2.0 * q).collect();
    let ancient = |id: String, sex, group, bp, country, lat, lon, dosage| Person {
        id,
        sex,
        group,
        modern: false,
        bp,
        country,
        lat,
        lon,
        dosage,
    };
    for i in 0..6 {
        let d = pseudo(&mut rnd, &hg002);
        people.push(ancient(
            format!("PL{i}.AG"),
            "F",
            "Planted_HG002",
            5000,
            "Testland",
            48,
            10,
            d,
        ));
        let d = pseudo(&mut rnd, &fr2);
        people.push(ancient(
            format!("FR{i}.AG"),
            "M",
            "Anc_French_like",
            7000,
            "Franceland",
            46,
            2,
            d,
        ));
        let d = pseudo(&mut rnd, &fi2);
        people.push(ancient(
            format!("FI{i}.AG"),
            "M",
            "Anc_Finnish_like",
            3500,
            "Finland",
            62,
            25,
            d,
        ));
    }
    let d = pseudo(&mut rnd, &vec![1.0; n]);
    people.push(ancient(
        "OUT0.AG".into(),
        "M",
        "Far_Away",
        5000,
        "China",
        40,
        116,
        d,
    ));

    let mut o = create(&format!("{dir}/aadr.anno"))?;
    writeln!(
        o,
        "Genetic ID\tIndividual ID\tDate mean in BP in years before 1950 CE\tDate standard deviation in BP\t\
         Full Date\tGroup ID\tLocality\tPolitical Entity\tLatitude\tLongitude\tData type\t\
         SNPs hit on autosomal targets (Computed using easystats on 1240k snpset)\tASSESSMENT\tASSESSMENT WARNINGS"
    )?;
    for p in &people {
        let short = p.id.split('.').next().unwrap_or("");
        let (kind, hits) = if p.modern {
            ("Shotgun.diploid", "1100000")
        } else {
            ("1240k", "500000")
        };
        writeln!(
            o,
            "{}\t{short}\t{}\t50\tsynthetic\t{}\tSomewhere\t{}\t{}\t{}\t{kind}\t{hits}\tPASS\t",
            p.id, p.bp, p.group, p.country, p.lat, p.lon
        )?;
    }
    o.flush()?;
    let mut o = create(&format!("{dir}/aadr.ind"))?;
    for p in &people {
        writeln!(o, "{:>20} {} {}", p.id, p.sex, p.group)?;
    }
    o.flush()?;
    let mut o = create(&format!("{dir}/aadr.snp"))?;
    for (i, (c, p, r, a, _)) in sites.iter().enumerate() {
        writeln!(
            o,
            "{:>20} {:>2} 0.000000 {:>12} {} {}",
            format!("syn{i}"),
            autosome(c).unwrap_or(0),
            p - OFF,
            *r as char,
            *a as char
        )?;
    }
    o.flush()?;
    // TGENO: a 48-byte header, then a record of ceil(n / 4) bytes a person,
    // two bits a site (first site in the high bits), 3 for missing.
    let mut o = create(&format!("{dir}/aadr.geno"))?;
    let mut head = format!("TGENO {:7} {:7} 0 0", people.len(), n).into_bytes();
    head.resize(head.len().max(48), 0);
    o.write_all(&head)?;
    for p in &people {
        let mut rec = vec![0xffu8; n.div_ceil(4)];
        for (k, &x) in p.dosage.iter().enumerate() {
            let v = if x == 9 { 3 } else { x };
            let sh = (3 - k % 4) * 2;
            rec[k / 4] = (rec[k / 4] & !(3 << sh)) | (v << sh);
        }
        o.write_all(&rec)?;
    }
    o.flush()?;
    let mut o = create(chain)?;
    let mut chroms: Vec<&str> = sites.iter().map(|x| x.0.as_str()).collect();
    chroms.sort_by_key(|c| autosome(c).unwrap_or(0));
    chroms.dedup();
    let mut k = 0;
    for c in chroms {
        let sz = fai.get(c).map_or(0, |f| f[0]);
        for (_, s, e) in slices.iter().filter(|x| x.0 == c) {
            k += 1;
            write!(
                o,
                "chain 1000 {c} {sz} + {} {} {c} {sz} + {s} {e} {k}\n{}\n\n",
                s - OFF,
                e - OFF,
                e - s
            )?;
        }
    }
    o.flush()?;
    let mut o = create(&format!("{dir}/expected_sites.tsv"))?;
    for (i, (c, p, r, a, _)) in sites.iter().enumerate() {
        writeln!(o, "syn{i}\t{c}\t{p}\t{}\t{}", *r as char, *a as char)?;
    }
    o.flush()?;
    writeln!(out, "synthetic panel: {n} sites, {} people", people.len())?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::super::testutil::{dir, fasta, fnv, run};

    #[test]
    fn panel_is_the_pythons_bytes() {
        let d = dir("ancient");
        fasta(
            &d.join("ref.fa"),
            &[("chr1", 30_000), ("chrX", 5_000), ("chr2", 30_000)],
            3,
        );
        std::fs::write(
            d.join("regions.bed"),
            "chr2\t2000\t28000\tb\nchr1\t1500\t25000\nchrX\t1500\t4000\nchr1\t500\t900\n",
        )
        .expect("bed");
        std::fs::write(d.join("conf.bed"), "chr1\t0\t30000\nchr2\t3000\t30000\n").expect("bed");
        // A truth SNV whose REF is the reference base, at every 13th base.
        let mut truth = String::from("#h\n");
        let mut fa = super::Fasta::open(d.join("ref.fa").to_str().expect("utf8")).expect("fa");
        for c in ["chr1", "chr2"] {
            for p in (1600..28000).step_by(13) {
                let r = fa.base(c, p).expect("read").expect("base") as char;
                let gt = ["0/1", "1|1", "1/0", "0/0"][(p % 4) as usize];
                truth.push_str(&format!(
                    "{c}\t{p}\t.\t{r}\tN\t50\tPASS\t.\tGT:DP\t{gt}:30\n"
                ));
            }
        }
        std::fs::write(d.join("truth.vcf"), truth).expect("truth");
        let p = |n: &str| d.join(n).to_str().expect("utf8").to_string();
        std::fs::create_dir_all(d.join("out")).expect("out");
        let (code, out, err) = run(
            "ancient-panel",
            &[
                &p("regions.bed"),
                &p("truth.vcf"),
                &p("conf.bed"),
                &p("ref.fa"),
                &p("out"),
                &p("chain"),
            ],
        );
        assert_eq!(code, 0, "{err}");
        assert_eq!(out, "synthetic panel: 3091 sites, 109 people\n");
        // FNV-1a of what the Python heredoc of ancient-similarity-1.sh wrote
        // from the same inputs.
        let pins = [
            ("out/aadr.anno", 0xba31_6084_c893_a392),
            ("out/aadr.ind", 0xbdd9_323d_97f2_62d0),
            ("out/aadr.snp", 0x6c21_020e_2e5d_3d88),
            ("out/aadr.geno", 0xa2bb_b113_80a3_177b),
            ("out/expected_sites.tsv", 0x99e9_3c28_9f00_6476),
            ("chain", 0xf2f0_f4b3_b11e_2b48),
        ];
        for (f, want) in pins {
            let got = fnv(&std::fs::read(d.join(f)).expect("read"));
            assert_eq!(got, want, "{f}: {got:#x}");
        }
    }
}
