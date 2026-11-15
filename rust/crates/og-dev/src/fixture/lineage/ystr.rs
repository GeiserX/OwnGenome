//! The Y cases: the private-SNV red control's seeded VCF, and the Y-STR
//! case's repeat check, read-length count and planted repeat cut.

use super::{Fail, Op, bad, cell_int, cigar, create, each_line, int, want};
use std::collections::{HashMap, HashSet};
use std::io::Write;

/// The fixture's TBL1Y stretch of chrY (scripts/ci/build-fixture.sh).
const FREE_FROM: i64 = 6_800_000;
const FREE_TO: i64 = 7_000_000;

/// Python's `s.isdigit()` for ASCII text.
fn is_digits(s: &str) -> bool {
    !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit())
}

/// `og dev fixture y-seed CALLS.vcf.gz MASK.fasta.gz PRIVATE.tsv OUT.vcf TABLE...`:
/// the VCF with two made-up haploid chrY calls (C>T) in the TBL1Y stretch at
/// positions no table and no record names: the first the callable mask
/// passes (`P`) and the first it does not, at the sample's median chrY depth
/// from PRIVATE's header (at least 10, 20 without one). Prints both
/// positions.
pub fn y_seed(args: &[String], out: &mut dyn Write, _err: &mut dyn Write) -> Result<(), Fail> {
    if args.len() < 4 {
        return Err(Fail::Usage("takes CALLS MASK PRIVATE OUT TABLE...".into()));
    }
    let (vcf, mask_path, private, dst) = (&args[0], &args[1], &args[2], &args[3]);
    let mut mask = Vec::new();
    each_line(mask_path, |l| {
        if !l.starts_with('>') {
            mask.extend_from_slice(l.trim().as_bytes());
        }
        Ok(())
    })?;
    let mut taken = HashSet::new();
    for t in &args[4..] {
        each_line(t, |l| {
            let c: Vec<&str> = l.split('\t').collect();
            if c.len() >= 4 && is_digits(c[3]) {
                taken.insert(cell_int(c[3])?);
            }
            Ok(())
        })?;
    }
    let (mut head, mut body) = (Vec::new(), Vec::new());
    each_line(vcf, |l| {
        if l.starts_with('#') {
            head.push(l.to_string());
        } else {
            let c: Vec<&str> = l.split('\t').collect();
            if (c[0] == "chrY" || c[0] == "Y") && c.len() > 1 && is_digits(c[1]) {
                taken.insert(cell_int(c[1])?);
            }
            body.push(l.to_string());
        }
        Ok(())
    })?;
    let text = std::fs::read_to_string(private).map_err(|e| bad(format!("{private}: {e}")))?;
    let depth = match text.find("median_chrY_depth=") {
        Some(i) => {
            let rest = &text[i + "median_chrY_depth=".len()..];
            let n = rest
                .find(|ch: char| !(ch.is_ascii_digit() || ch == '.'))
                .unwrap_or(rest.len());
            match &rest[..n] {
                "" => 20,
                v => {
                    let x = pycompat::num::parse_float(v)
                        .map_err(|_| bad(format!("{private}: median depth '{v}'")))?;
                    (x.trunc() as i64).max(10)
                }
            }
        }
        None => 20,
    };
    let mask_at = |p: i64| -> Result<u8, Fail> {
        usize::try_from(p - 1)
            .ok()
            .and_then(|i| mask.get(i).copied())
            .ok_or_else(|| bad(format!("{mask_path}: shorter than {p}")))
    };
    let (mut pin, mut pout) = (None, None);
    for p in (FREE_FROM..FREE_TO).filter(|p| !taken.contains(p)) {
        let passed = mask_at(p)? == b'P';
        if passed && pin.is_none() {
            pin = Some(p);
        } else if !passed && pout.is_none() {
            pout = Some(p);
        }
        if pin.is_some() && pout.is_some() {
            break;
        }
    }
    let (Some(pin), Some(pout)) = (pin, pout) else {
        return Err(bad("no free position inside and outside the mask"));
    };
    let mut o = create(dst)?;
    for l in head.iter().chain(&body) {
        writeln!(o, "{l}")?;
    }
    for p in [pin, pout] {
        writeln!(
            o,
            "chrY\t{p}\t.\tC\tT\t50\tPASS\t.\tGT:AD:DP:GQ\t1:0,{depth}:{depth}:60"
        )?;
    }
    o.flush()?;
    writeln!(out, "{pin} {pout}")?;
    Ok(())
}

/// `og dev fixture ystr-periodic REGIONS.bed REGIONS.fa`: how many regions
/// (`CHROM START END PERIOD ...`, sequences from `samtools faidx -r`) repeat
/// with their period: at least 75% of the bases equal to the base one
/// period further on.
pub fn periodic(args: &[String], out: &mut dyn Write, _err: &mut dyn Write) -> Result<(), Fail> {
    want(args, 2)?;
    let mut seqs: HashMap<String, Vec<u8>> = HashMap::new();
    let mut name: Option<String> = None;
    each_line(&args[1], |l| {
        let l = l.trim();
        if let Some(h) = l.strip_prefix('>') {
            seqs.insert(h.to_string(), Vec::new());
            name = Some(h.to_string());
        } else if let Some(n) = &name {
            seqs.get_mut(n)
                .expect("inserted with its header")
                .extend(l.bytes().map(|b| b.to_ascii_uppercase()));
        }
        Ok(())
    })?;
    let mut ok = 0u64;
    each_line(&args[0], |l| {
        let c: Vec<&str> = l.split('\t').collect();
        if c.len() < 4 {
            return Err(bad(format!(
                "{}: a region with fewer than 4 columns",
                args[0]
            )));
        }
        let s = seqs
            .get(&format!("{}:{}-{}", c[0], c[1], c[2]))
            .map_or(&[][..], Vec::as_slice);
        let p = usize::try_from(cell_int(c[3])?).map_err(|_| bad("a negative period"))?;
        let span = cell_int(c[2])? - cell_int(c[1])? + 1;
        if s.len() as i64 == span && s.len() > p {
            let n = s.len() - p;
            let same = (0..n).filter(|&i| s[i] == s[i + p]).count();
            ok += u64::from(same as f64 >= 0.75 * n as f64);
        }
        Ok(())
    })?;
    writeln!(out, "{ok}")?;
    Ok(())
}

/// The reference end of an alignment.
fn end_of(pos: i64, ops: &[Op]) -> i64 {
    pos + ops
        .iter()
        .filter(|(_, op)| b"MDN=X".contains(op))
        .map(|(n, _)| n)
        .sum::<i64>()
        - 1
}

/// One callable marker of the panel.
struct Region {
    start: i64,
    end: i64,
    period: i64,
    ref_copies: i64,
    kind: String,
}

/// `og dev fixture ystr-readcount PANEL CALLS READS.sam`: for each marker
/// step 45 called high, the count its spanning reads give (10 bases of flank
/// each side; the region's length plus the CIGAR's insertions minus its
/// deletions in the region or its flanks; the modal length), as `marker
/// step-count read-count reads kind`.
pub fn readcount(args: &[String], out: &mut dyn Write, _err: &mut dyn Write) -> Result<(), Fail> {
    want(args, 3)?;
    let mut regions: Vec<(String, Region)> = Vec::new();
    each_line(&args[0], |l| {
        let c: Vec<&str> = l.split('\t').collect();
        if l.starts_with('#')
            || c[0] == "order"
            || !c.get(3).is_some_and(|k| *k == "single" || *k == "span")
        {
            return Ok(());
        }
        let parts: Vec<Vec<&str>> = c
            .get(4)
            .ok_or_else(|| bad(format!("{}: a marker without parts", args[0])))?
            .split(';')
            .map(|x| x.splitn(4, '/').collect())
            .collect();
        let span_of = |part: &Vec<&str>| -> Result<(i64, i64), Fail> {
            let range = part[0]
                .split(':')
                .nth(1)
                .ok_or_else(|| bad(format!("part '{}'", part[0])))?;
            let (s, e) = range
                .split_once('-')
                .ok_or_else(|| bad(format!("part '{}'", part[0])))?;
            Ok((cell_int(s)?, cell_int(e)?))
        };
        let start = span_of(&parts[0])?.0;
        let end = span_of(parts.last().expect("one part at least"))?.1;
        let period = cell_int(
            parts[0]
                .get(1)
                .ok_or_else(|| bad("a part without a period"))?,
        )?;
        let mut ref_copies = 0;
        for x in &parts {
            let copies = x.get(2).ok_or_else(|| bad("a part without copies"))?;
            let v = pycompat::num::parse_float(copies)
                .map_err(|_| bad(format!("copies '{copies}'")))?;
            ref_copies += (v + 0.5).trunc() as i64;
        }
        let r = Region {
            start,
            end,
            period,
            ref_copies,
            kind: c[3].to_string(),
        };
        match regions.iter_mut().find(|(m, _)| m == c[2]) {
            Some(slot) => slot.1 = r,
            None => regions.push((c[2].to_string(), r)),
        }
        Ok(())
    })?;
    // Per marker, the length changes in the order first seen, with counts.
    let mut lens: HashMap<String, Vec<(i64, u64)>> = HashMap::new();
    each_line(&args[2], |l| {
        let f: Vec<&str> = l.split('\t').collect();
        if f.len() < 6 || f[2] != "chrY" || f[5] == "*" {
            return Ok(());
        }
        let pos = cell_int(f[3])?;
        let ops = cigar(f[5]);
        let end = end_of(pos, &ops);
        for (m, g) in &regions {
            if pos > g.start - 10 || end < g.end + 10 {
                continue;
            }
            let (mut r, mut diff) = (pos, 0i64);
            let near = |r: i64| g.start - 10 <= r && r <= g.end + 10;
            for &(n, op) in &ops {
                match op {
                    b'M' | b'=' | b'X' => r += n,
                    b'D' | b'N' => {
                        if near(r) {
                            diff -= n;
                        }
                        r += n;
                    }
                    b'I' if near(r) => diff += n,
                    _ => {}
                }
            }
            let v = lens.entry(m.clone()).or_default();
            match v.iter_mut().find(|x| x.0 == diff) {
                Some(x) => x.1 += 1,
                None => v.push((diff, 1)),
            }
        }
        Ok(())
    })?;
    each_line(&args[1], |l| {
        let c: Vec<&str> = l.split('\t').collect();
        if c.len() < 5 || c[4] != "high" {
            return Ok(());
        }
        let (Some((_, g)), Some(v)) = (regions.iter().find(|(m, _)| m == c[2]), lens.get(c[2]))
        else {
            return Ok(());
        };
        // Counter.most_common(1): the first of the highest counts.
        let Some(&(diff, n)) = v.iter().fold(None::<&(i64, u64)>, |best, x| match best {
            Some(b) if b.1 >= x.1 => Some(b),
            _ => Some(x),
        }) else {
            return Ok(());
        };
        let (whole, rest) = (diff.div_euclid(g.period), diff.rem_euclid(g.period));
        let frac = if rest != 0 {
            format!(".{rest}")
        } else {
            String::new()
        };
        writeln!(
            out,
            "{}\t{}\t{}{frac}\t{n}\t{}",
            c[2],
            c[3],
            g.ref_copies + whole,
            g.kind
        )?;
        Ok(())
    })
}

/// `og dev fixture ystr-cut START END READS.sam`: the SAM with one 4-base
/// repeat unit cut out at START: chrY reads that span START..END (10 bases
/// of flank each side) with a plain match get bases START..START+3 removed
/// (CIGAR xM4DyM, MD and NM dropped); other reads over the region are
/// dropped. The count goes to stderr.
pub fn cut(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> Result<(), Fail> {
    want(args, 3)?;
    let (s, e): (i64, i64) = (int(&args[0], "START")?, int(&args[1], "END")?);
    let mut n_cut = 0u64;
    each_line(&args[2], |l| {
        if l.starts_with('@') {
            writeln!(out, "{l}")?;
            return Ok(());
        }
        let mut f: Vec<String> = l.split('\t').map(str::to_string).collect();
        if f.len() < 11 {
            return Err(bad(format!(
                "{}: a read with fewer than 11 columns",
                args[2]
            )));
        }
        if f[2] != "chrY" || f[5] == "*" {
            writeln!(out, "{l}")?;
            return Ok(());
        }
        let pos = cell_int(&f[3])?;
        let ops = cigar(&f[5]);
        let end = end_of(pos, &ops);
        if end < s || pos > e {
            writeln!(out, "{l}")?;
            return Ok(());
        }
        let plain = f[5].strip_suffix('M').filter(|d| is_digits(d));
        let Some(m) = plain else { return Ok(()) };
        if pos > s - 10 || end < e + 10 {
            return Ok(());
        }
        let n = cell_int(m)?;
        let off = usize::try_from(s - pos).map_err(|_| bad("a negative offset"))?;
        let snip = |x: &str| -> String {
            let b = x.as_bytes();
            let head = &b[..off.min(b.len())];
            let tail = b.get(off + 4..).unwrap_or(&[]);
            String::from_utf8_lossy(&[head, tail].concat()).into_owned()
        };
        f[9] = snip(&f[9]);
        if f[10] != "*" {
            f[10] = snip(&f[10]);
        }
        f[5] = format!("{off}M4D{}M", n - off as i64 - 4);
        f.retain(|x| !x.starts_with("MD:Z:") && !x.starts_with("NM:i:"));
        writeln!(out, "{}", f.join("\t"))?;
        n_cut += 1;
        Ok(())
    })?;
    writeln!(err, "planted the cut in {n_cut} reads")?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::super::testutil::{dir, run};

    #[test]
    fn periodic_counts_repeats_of_their_period() {
        let d = dir("periodic");
        std::fs::write(d.join("r.bed"), "chrY\t10\t21\t4\tx\tDYS1\nchrY\t30\t41\t4\tx\tDYS2\nchrY\t50\t61\t4\nchrY\t70\t72\t4\n").expect("bed");
        std::fs::write(d.join("r.fa"), ">chrY:10-21\nAGATAGAT\nagat\n>chrY:30-41\nACGTTGCAACGG\n>chrY:50-61\nAGATAGATAGA\n>chrY:70-72\nAAA\n").expect("fa");
        let p = |n: &str| d.join(n).to_str().expect("utf8").to_string();
        assert_eq!(run("ystr-periodic", &[&p("r.bed"), &p("r.fa")]).1, "1\n");
        std::fs::write(d.join("one.bed"), "chrY\t10\t21\t1\n").expect("bed");
        assert_eq!(run("ystr-periodic", &[&p("one.bed"), &p("r.fa")]).1, "0\n");
    }

    #[test]
    fn readcount_takes_the_modal_length_change() {
        let d = dir("readcount");
        std::fs::write(
            d.join("panel.tsv"),
            "#c\norder\tpanel\tmarker\tkind\tparts\n1\tY-12\tDYS1\tsingle\tchrY:100-139/4/10\n2\tY-12\tDYS2\tspan\tchrY:200-219/4/5.5;chrY:230-249/4/4.4\n3\tY-12\tDYS3\tmulti\tchrY:300-310/4/3\n",
        )
        .expect("panel");
        std::fs::write(
            d.join("calls.tsv"),
            "order\tpanel\tmarker\tcount\tconfidence\n1\tY-12\tDYS1\t9\thigh\n2\tY-12\tDYS2\t12\thigh\n3\tY-12\tDYS3\t3\thigh\n",
        )
        .expect("calls");
        let mut sam = String::new();
        // DYS1: two reads lose 4 bases (one unit), one read gains 2.
        for (i, cig) in ["40M4D60M", "40M4D60M", "40M2I58M", "100M"]
            .iter()
            .enumerate()
        {
            sam.push_str(&format!("r{i}\t0\tchrY\t60\t60\t{cig}\t*\t0\t0\t*\t*\n"));
        }
        // DYS2 (190..259 with flanks): one read, no change.
        sam.push_str("s\t0\tchrY\t150\t60\t150M\t*\t0\t0\t*\t*\n");
        sam.push_str("t\t0\tchr1\t150\t60\t150M\t*\t0\t0\t*\t*\n");
        std::fs::write(d.join("reads.sam"), sam).expect("sam");
        let p = |n: &str| d.join(n).to_str().expect("utf8").to_string();
        let (code, out, err) = run(
            "ystr-readcount",
            &[&p("panel.tsv"), &p("calls.tsv"), &p("reads.sam")],
        );
        assert_eq!(code, 0, "{err}");
        assert_eq!(out, "DYS1\t9\t9\t2\tsingle\nDYS2\t12\t10\t1\tspan\n");
    }

    #[test]
    fn cut_removes_one_unit_from_plain_spanning_reads() {
        let d = dir("cut");
        let sam = "@HD\tVN:1.6\n\
                   a\t0\tchrY\t90\t60\t40M\t*\t0\t0\tAAAAAAAAAACCCCGGGGGGGGGGGGGGGGGGGGGGGGGG\tIIIIIIIIIIJJJJIIIIIIIIIIIIIIIIIIIIIIIIII\tMD:Z:40\tNM:i:0\tRG:Z:x\n\
                   b\t0\tchrY\t90\t60\t20M1I19M\t*\t0\t0\t*\t*\n\
                   c\t0\tchrY\t10\t60\t20M\t*\t0\t0\t*\t*\n\
                   d\t0\tchr1\t90\t60\t40M\t*\t0\t0\t*\t*\n";
        std::fs::write(d.join("in.sam"), sam).expect("sam");
        let (code, out, err) = run(
            "ystr-cut",
            &["100", "105", d.join("in.sam").to_str().expect("utf8")],
        );
        assert_eq!(code, 0, "{err}");
        assert_eq!(
            out,
            "@HD\tVN:1.6\na\t0\tchrY\t90\t60\t10M4D26M\t*\t0\t0\tAAAAAAAAAAGGGGGGGGGGGGGGGGGGGGGGGGGG\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\tRG:Z:x\n\
             c\t0\tchrY\t10\t60\t20M\t*\t0\t0\t*\t*\nd\t0\tchr1\t90\t60\t40M\t*\t0\t0\t*\t*\n"
        );
        assert_eq!(err, "planted the cut in 1 reads\n");
    }

    #[test]
    fn y_seed_picks_free_positions_in_and_out_of_the_mask() {
        let d = dir("y-seed");
        // The mask passes nothing before 6,800,003; 6,800,000 is in a table,
        // 6,800,001 is a record, so OUT is 6,800,002 and IN 6,800,004.
        let mut mask = vec![b'N'; 7_000_000];
        mask[6_800_002] = b'P';
        mask[6_800_003] = b'P';
        let mut fa = b">chrY\n".to_vec();
        for chunk in mask.chunks(60) {
            fa.extend_from_slice(chunk);
            fa.push(b'\n');
        }
        std::fs::write(d.join("mask.fa"), fa).expect("mask");
        std::fs::write(d.join("t1.txt"), "chrY\tA1\tR\t6800000\tC->T\tC\tT\nbad\n").expect("t1");
        std::fs::write(d.join("t2.txt"), "chrY\tA2\tR\t6800003\tC->T\tC\tT\n").expect("t2");
        std::fs::write(
            d.join("c.vcf"),
            "##fileformat=VCFv4.2\n#CHROM\tPOS\nchrY\t6800001\t.\tA\tG\n",
        )
        .expect("vcf");
        std::fs::write(
            d.join("p.tsv"),
            "# gates\tmin_depth=10 median_chrY_depth=31.9 x\n",
        )
        .expect("priv");
        let p = |n: &str| d.join(n).to_str().expect("utf8").to_string();
        let (code, out, err) = run(
            "y-seed",
            &[
                &p("c.vcf"),
                &p("mask.fa"),
                &p("p.tsv"),
                &p("o.vcf"),
                &p("t1.txt"),
                &p("t2.txt"),
            ],
        );
        assert_eq!((code, out.as_str()), (0, "6800004 6800002\n"), "{err}");
        assert_eq!(
            std::fs::read_to_string(d.join("o.vcf")).expect("vcf"),
            "##fileformat=VCFv4.2\n#CHROM\tPOS\nchrY\t6800001\t.\tA\tG\n\
             chrY\t6800004\t.\tC\tT\t50\tPASS\t.\tGT:AD:DP:GQ\t1:0,31:31:60\n\
             chrY\t6800002\t.\tC\tT\t50\tPASS\t.\tGT:AD:DP:GQ\t1:0,31:31:60\n"
        );
        std::fs::write(d.join("p.tsv"), "no depth\n").expect("priv");
        let (_, out, _) = run(
            "y-seed",
            &[&p("c.vcf"), &p("mask.fa"), &p("p.tsv"), &p("o.vcf")],
        );
        assert_eq!(out, "6800003 6800000\n");
        assert!(
            std::fs::read_to_string(d.join("o.vcf"))
                .expect("vcf")
                .contains("1:0,20:20:60")
        );
    }
}
