//! `og ancestry-fine simulate`: an admixed genome with known truth, built
//! from held-out people of the panel, to check the painting and the
//! projection against. The model and the files are described in
//! [`og_admix::simulate`].

use crate::{ArgError, Args, Fail, Result, Spec, repr};
use og_admix::map::GeneticMap;
use og_admix::simulate::{self, Config};
use std::io::Write;
use std::path::{Path, PathBuf};

pub(crate) const SPEC: Spec = Spec {
    verb: "simulate",
    usage: "usage: og ancestry-fine simulate --panel DIR --maps DIR --people FILE --mu GROUP=P[,GROUP=P...] --generations T --out OUT.vcf[.gz] --truth OUT.bed [--chroms chrA,chrB] [--switch-rate R] [--seed N] [--sample NAME]\n",
    options: &[
        "--panel",
        "--maps",
        "--people",
        "--mu",
        "--generations",
        "--out",
        "--truth",
        "--chroms",
        "--switch-rate",
        "--seed",
        "--sample",
    ],
    required: &[
        "--panel",
        "--maps",
        "--people",
        "--mu",
        "--generations",
        "--out",
        "--truth",
    ],
};

/// `og ancestry-fine simulate ARGS`: the exit code.
pub fn cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let a = match Args::parse(&SPEC, args).and_then(|a| {
        a.float("--generations", 0.0)?;
        a.float("--switch-rate", 0.0)?;
        a.int("--seed", 0)?;
        Ok(a)
    }) {
        Ok(a) => a,
        Err(ArgError::Help) => {
            let _ = write!(out, "{}", SPEC.usage);
            return 0;
        }
        Err(ArgError::Usage(msg)) => {
            let _ = write!(err, "{}", SPEC.usage);
            let _ = writeln!(err, "og ancestry-fine {}: error: {msg}", SPEC.verb);
            return 2;
        }
    };
    match main(&a, out) {
        Ok(()) => 0,
        Err(Fail(msg)) => {
            let _ = out.flush();
            let _ = writeln!(err, "ERROR: {msg}");
            1
        }
    }
}

/// `GROUP=P,GROUP=P`, in the order given.
fn parse_mu(s: &str) -> Result<Vec<(String, f64)>> {
    s.split(',')
        .map(|item| {
            let (g, p) = item.split_once('=').unwrap_or((item, ""));
            let p = pycompat::num::parse_float(p).ok().filter(|_| !g.is_empty());
            p.map(|p| (g.trim().to_owned(), p))
                .ok_or_else(|| Fail(format!("--mu: {} is not GROUP=PROPORTION", repr(item))))
        })
        .collect()
}

/// The held-out people: sample and group, the first two columns of each
/// line (`paint_map.tsv`'s layout); blank lines and `#` lines skipped.
fn read_people(path: &Path) -> Result<Vec<(String, String)>> {
    let mut out = Vec::new();
    let mut n = 0u64;
    crate::for_each_line(path, |line| {
        n += 1;
        if line.trim().is_empty() || line.starts_with('#') {
            return Ok(());
        }
        let c: Vec<&str> = crate::split_ws(line).collect();
        if c.len() < 2 {
            return Err(Fail(format!(
                "{}: line {n}: a sample and its group, not {}",
                path.display(),
                repr(line)
            )));
        }
        out.push((c[0].to_owned(), c[1].to_owned()));
        Ok(())
    })?;
    if out.is_empty() {
        return Err(Fail(format!("{}: no held-out people", path.display())));
    }
    Ok(out)
}

/// The chromosomes asked for, else every `<chrom>.vcf.gz` of the panel
/// that has a `<chrom>.map`, in chromosome order.
fn chroms(a: &Args, panel: &Path, maps: &Path) -> Result<Vec<String>> {
    if let Some(list) = a.get("--chroms") {
        let c: Vec<String> = list
            .split([',', ' '])
            .filter(|c| !c.is_empty())
            .map(str::to_owned)
            .collect();
        if c.is_empty() {
            return Err(Fail("--chroms names no chromosome".into()));
        }
        return Ok(c);
    }
    let entries = std::fs::read_dir(panel).map_err(|e| crate::cannot_read(panel, &e))?;
    let mut out = Vec::new();
    for e in entries {
        let e = e.map_err(|e| crate::cannot_read(panel, &e))?;
        let name = e.file_name().to_string_lossy().into_owned();
        if let Some(c) = name.strip_suffix(".vcf.gz")
            && maps.join(format!("{c}.map")).is_file()
        {
            out.push(c.to_owned());
        }
    }
    if out.is_empty() {
        return Err(Fail(format!(
            "no <chrom>.vcf.gz in {} with a <chrom>.map in {}",
            panel.display(),
            maps.display()
        )));
    }
    out.sort_by(|x, y| crate::chrom_key(x).cmp(&crate::chrom_key(y)));
    Ok(out)
}

pub(crate) fn main(a: &Args, out: &mut dyn Write) -> Result<()> {
    let seed = a.int("--seed", 1).unwrap_or(1);
    let cfg = Config {
        generations: a.float("--generations", 0.0).unwrap_or(0.0),
        switch_rate: a.float("--switch-rate", 1.0).unwrap_or(1.0),
        mu: parse_mu(a.req("--mu"))?,
        seed: u64::try_from(seed)
            .map_err(|_| Fail(format!("--seed {seed}: a seed is 0 or more")))?,
    };
    cfg.check().map_err(|e| Fail(e.0))?;
    let people = read_people(Path::new(a.req("--people")))?;
    let (panel, maps) = (Path::new(a.req("--panel")), Path::new(a.req("--maps")));
    let chroms = chroms(a, panel, maps)?;
    let sample = a.get("--sample").unwrap_or("SIM");
    if sample.is_empty() || sample.contains(['\t', '\n']) {
        return Err(Fail(format!(
            "--sample {}: not a sample name",
            repr(sample)
        )));
    }
    let vcf_path = PathBuf::from(a.req("--out"));
    let mut pending = crate::Pending::new();
    let file = pending.create(vcf_path.clone())?;
    let write_err = |e: std::io::Error| crate::cannot_write(&vcf_path, &e);
    let names: Vec<&str> = chroms.iter().map(String::as_str).collect();
    let mut done = Vec::with_capacity(chroms.len());
    let mut run = |vcf: &mut dyn Write| -> Result<()> {
        simulate::write_header(vcf, sample, &names, &cfg).map_err(|e| Fail(e.0))?;
        for c in &chroms {
            let map_path = maps.join(format!("{c}.map"));
            let map = GeneticMap::read(&map_path).map_err(|e| Fail(e.0))?;
            let panel_path = panel.join(format!("{c}.vcf.gz"));
            let chrom = simulate::simulate_chrom(&panel_path, &map, &people, &cfg, vcf)
                .map_err(|e| Fail(e.0))?;
            if chrom.name != *c {
                return Err(Fail(format!(
                    "{}: its records are on {}, not {c}",
                    panel_path.display(),
                    chrom.name
                )));
            }
            done.push(chrom);
        }
        Ok(())
    };
    if vcf_path.extension().is_some_and(|e| e == "gz") {
        let mut bgzf =
            og_io::bgzf::Writer::new(&mut *file, og_io::bgzf::workers()).map_err(write_err)?;
        run(&mut bgzf)?;
        bgzf.finish().map_err(write_err)?;
    } else {
        run(file)?;
    }
    let mut bed = Vec::new();
    simulate::write_truth(&mut bed, &done, &cfg).map_err(|e| Fail(e.0))?;
    pending.put(PathBuf::from(a.req("--truth")), &bed)?;
    pending.commit_all()?;
    let mut bp = vec![0i64; cfg.mu.len()];
    for c in &done {
        for r in &c.rows {
            bp[r.ancestry] += r.end - r.start;
        }
        writeln!(
            out,
            "{}: {} sites, {} panel records left out, {} truth rows",
            c.name,
            c.sites,
            c.skipped,
            c.rows.len()
        )
        .map_err(|e| Fail(format!("cannot write stdout: {e}")))?;
    }
    let total: i64 = bp.iter().sum();
    let shares: Vec<String> = cfg
        .mu
        .iter()
        .zip(&bp)
        .map(|((g, _), b)| format!("{g} {:.4}", *b as f64 / total.max(1) as f64))
        .collect();
    writeln!(
        out,
        "ancestry shares of the truth (bp): {}",
        shares.join(", ")
    )
    .map_err(|e| Fail(format!("cannot write stdout: {e}")))
}

#[cfg(test)]
mod tests {
    use crate::testutil::{Dir, run};

    fn setup(d: &Dir) {
        let mut a = String::from(
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tY1\tH1\tO1\n",
        );
        let mut b = a.clone();
        for i in 1..=300u32 {
            a.push_str(&format!(
                "chr21\t{}\t.\tA\tG\t.\t.\t.\tGT\t0|0\t1|1\t0|1\n",
                i * 1000
            ));
            b.push_str(&format!(
                "chr2\t{}\t.\tG\tT\t.\t.\t.\tGT\t1|1\t0|0\t0|1\n",
                i * 1000
            ));
        }
        d.put("panel/chr21.vcf.gz", &a);
        d.put("panel/chr2.vcf.gz", &b);
        d.put("panel/chr9.vcf.gz", &b);
        d.put("maps/chr21.map", "chr21 . 0.0 1\nchr21 . 30.0 300000\n");
        d.put("maps/chr2.map", "chr2 . 0.0 1\nchr2 . 30.0 300000\n");
        d.put(
            "people.tsv",
            "# sample\tgroup\nY1\tAFR\nH1\tEAS\n\nO1\tOCE\n",
        );
    }

    fn args<'a>(d: &'a Dir, extra: &[&'a str]) -> Vec<String> {
        let mut v: Vec<String> = [
            "simulate",
            "--panel",
            &d.path("panel"),
            "--maps",
            &d.path("maps"),
            "--people",
            &d.path("people.tsv"),
            "--mu",
            "AFR=0.5,EAS=0.5",
            "--generations",
            "20",
            "--truth",
            &d.path("truth.bed"),
        ]
        .iter()
        .map(|s| (*s).to_owned())
        .collect();
        v.extend(extra.iter().map(|s| (*s).to_owned()));
        v
    }

    fn run_v(v: &[String]) -> (u8, String, String) {
        let refs: Vec<&str> = v.iter().map(String::as_str).collect();
        run(&refs)
    }

    #[test]
    fn writes_the_genome_and_its_truth() {
        let d = Dir::new("simulate");
        setup(&d);
        let vcf = d.path("sim.vcf");
        let (rc, out, err) = run_v(&args(&d, &["--out", &vcf, "--seed", "3"]));
        assert_eq!((rc, err.as_str()), (0, ""));
        // chr9 has no map: by default only chr2 and chr21, in that order.
        assert!(
            out.starts_with("chr2: 300 sites, 0 panel records left out, "),
            "{out}"
        );
        assert!(out.contains("\nchr21: 300 sites"), "{out}");
        assert!(
            out.contains("ancestry shares of the truth (bp): AFR 0."),
            "{out}"
        );
        let text = d.read("sim.vcf");
        assert!(text.contains("##contig=<ID=chr2>\n##contig=<ID=chr21>\n"));
        assert!(text.contains("\tFORMAT\tSIM\n"));
        assert_eq!(text.lines().filter(|l| !l.starts_with('#')).count(), 600);
        let bed = d.read("truth.bed");
        assert!(bed.starts_with("#chrom\tstart\tend\thaplotype\tancestry\nchr2\t999\t"));
        assert!(
            bed.lines()
                .skip(1)
                .all(|l| l.ends_with("\tAFR") || l.ends_with("\tEAS"))
        );
        // Gzip output is BGZF, and the same seed gives the same genome.
        let gz = d.path("sim.vcf.gz");
        let (rc, _, err) = run_v(&args(&d, &["--out", &gz, "--seed", "3"]));
        assert_eq!((rc, err.as_str()), (0, ""));
        let raw = std::fs::read(&gz).unwrap();
        assert_eq!(
            &raw[..4],
            &[0x1f, 0x8b, 8, 4],
            "a BGZF member with its extra field"
        );
        let mut back = String::new();
        std::io::Read::read_to_string(&mut og_io::open(&gz).unwrap(), &mut back).unwrap();
        assert_eq!(back, text);
        // One chromosome asked for, another sample name.
        let (rc, out, _) = run_v(&args(
            &d,
            &["--out", &vcf, "--chroms", "chr21", "--sample", "MIX"],
        ));
        assert_eq!(rc, 0);
        assert!(out.starts_with("chr21: "), "{out}");
        assert!(d.read("sim.vcf").contains("\tFORMAT\tMIX\n"));
    }

    #[test]
    fn refuses_bad_settings_before_writing() {
        let d = Dir::new("simulate-bad");
        setup(&d);
        let vcf = d.path("sim.vcf");
        let (rc, _, err) = run_v(&args(&d, &["--out", &vcf, "--chroms", "chr9"]));
        assert_eq!(rc, 1);
        assert!(err.starts_with("ERROR: cannot read "), "{err}");
        let mut v = args(&d, &["--out", &vcf]);
        v[8] = "AFR=0.5,EAS".into();
        let (rc, _, err) = run_v(&v);
        assert_eq!(rc, 1);
        assert!(err.contains("--mu: 'EAS' is not GROUP=PROPORTION"), "{err}");
        v[8] = "AFR=0.5,EAS=0.6".into();
        assert!(run_v(&v).2.contains("sum to 1.1, not 1"));
        v[8] = "AFR=0.5,SIB=0.5".into();
        assert!(run_v(&v).2.contains("no held-out person of ancestry 'SIB'"));
        let (rc, _, err) = run_v(&args(&d, &["--out", &vcf, "--seed", "-2"]));
        assert_eq!(rc, 1);
        assert!(err.contains("--seed -2: a seed is 0 or more"), "{err}");
        let (rc, _, err) = run_v(&args(&d, &["--out", &vcf, "--generations", "x"]));
        assert_eq!(rc, 2);
        assert!(
            err.contains("argument --generations: invalid float value: 'x'"),
            "{err}"
        );
        let (rc, _, err) = run(&["simulate", "--og-scaffold-probe"]);
        assert_eq!(rc, 2);
        assert!(
            err.contains("the following arguments are required"),
            "{err}"
        );
        assert!(!std::path::Path::new(&vcf).exists(), "nothing written");
        assert!(!std::path::Path::new(&d.path("truth.bed")).exists());
        let (rc, out, _) = run(&["simulate", "--help"]);
        assert_eq!(rc, 0);
        assert!(out.starts_with("usage: og ancestry-fine simulate"));
    }
}
