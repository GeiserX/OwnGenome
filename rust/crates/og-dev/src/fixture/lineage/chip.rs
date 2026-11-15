//! The CHIP case's planted reads, and the chip chrM round trip's GRCh37
//! side and its comparison.

use super::{Fail, Op, bad, cell_int, cigar, create, each_line, fai_contigs, int, want};
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::io::Write;

/// The read offset of 1-based reference `pos` in a read starting at `start`,
/// when an M, = or X operation covers it.
fn offset_at(ops: &[Op], start: i64, pos: i64) -> Option<i64> {
    let (mut r, mut q, mut hit) = (start, 0i64, None);
    for &(n, op) in ops {
        match op {
            b'M' | b'=' | b'X' => {
                if r <= pos && pos < r + n {
                    hit = Some(q + pos - r);
                }
                r += n;
                q += n;
            }
            b'I' | b'S' => q += n,
            b'D' | b'N' => r += n,
            _ => {}
        }
    }
    hit
}

/// `og dev fixture plant-snv EVERY FAI READS.sam POS`: the reads of a SAM
/// body as a SAM with the header of FAI's contigs and one read group
/// (`giab`, sample HG002). Reads on a contig FAI lacks are dropped, a mate
/// on such a contig is unset, and on chr2 every EVERY-th read with a C at
/// POS gets a T there (0 plants nothing).
pub fn plant_snv(args: &[String], out: &mut dyn Write, _err: &mut dyn Write) -> Result<(), Fail> {
    want(args, 4)?;
    let every: u64 = int(&args[0], "EVERY")?;
    let pos: i64 = int(&args[3], "POS")?;
    let contigs = fai_contigs(&args[1])?;
    let names: HashSet<&str> = contigs.iter().map(|c| c.0.as_str()).collect();
    writeln!(out, "@HD\tVN:1.6\tSO:unsorted")?;
    for (c, l) in &contigs {
        writeln!(out, "@SQ\tSN:{c}\tLN:{l}")?;
    }
    writeln!(out, "@RG\tID:giab\tSM:HG002\tPL:ILLUMINA")?;
    let mut k = 0u64;
    each_line(&args[2], |l| {
        let mut f: Vec<String> = l.split('\t').map(str::to_string).collect();
        if f.len() < 11 {
            return Err(bad(format!(
                "{}: a read with fewer than 11 columns",
                args[2]
            )));
        }
        if !names.contains(f[2].as_str()) {
            return Ok(());
        }
        f.retain(|x| !x.starts_with("RG:Z:"));
        f.push("RG:Z:giab".into());
        if f[6] != "=" && f[6] != "*" && !names.contains(f[6].as_str()) {
            f[6] = "*".into();
            f[7] = "0".into();
        }
        if every > 0 && f[2] == "chr2" && f[5] != "*" {
            let hit = offset_at(&cigar(&f[5]), cell_int(&f[3])?, pos);
            if let Some(i) = hit.and_then(|h| usize::try_from(h).ok())
                && f[9].as_bytes().get(i) == Some(&b'C')
            {
                if k.is_multiple_of(every) {
                    f[9].replace_range(i..=i, "T");
                }
                k += 1;
            }
        }
        writeln!(out, "{}", f.join("\t"))?;
        Ok(())
    })
}

/// `og dev fixture base-count POS BASE [READS.sam]`: `A N`, where N is the
/// reads of a SAM body (`samtools view`; stdin without a file) whose
/// alignment covers POS with a read base, and A those of them with BASE
/// there.
pub fn base_count(args: &[String], out: &mut dyn Write, _err: &mut dyn Write) -> Result<(), Fail> {
    if !(2..=3).contains(&args.len()) {
        return Err(Fail::Usage("takes POS BASE [READS.sam]".into()));
    }
    let pos: i64 = int(&args[0], "POS")?;
    let base = match args[1].as_bytes() {
        [b] => *b,
        _ => return Err(Fail::Usage("BASE is one letter".into())),
    };
    let (mut a, mut n) = (0u64, 0u64);
    each_line(args.get(2).map_or("-", String::as_str), |line| {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 10 || f[5] == "*" {
            return Ok(());
        }
        if let Some(i) = offset_at(&cigar(f[5]), cell_int(f[3])?, pos) {
            n += 1;
            a += u64::from(
                usize::try_from(i).ok().and_then(|i| f[9].as_bytes().get(i)) == Some(&base),
            );
        }
        Ok(())
    })?;
    writeln!(out, "{a} {n}")?;
    Ok(())
}

/// `og dev fixture sites-fasta SITES.tsv.gz OUT.fasta`: a GRCh37 FASTA (and
/// its `.fai`) with the consumer table's base at each site (`id chrom pos
/// base` rows) and N elsewhere, contigs 1 to 22, X, Y, MT in that order,
/// each one base longer than its last site.
pub fn sites_fasta(
    args: &[String],
    _out: &mut dyn Write,
    _err: &mut dyn Write,
) -> Result<(), Fail> {
    want(args, 2)?;
    let mut bases: HashMap<String, BTreeMap<usize, u8>> = HashMap::new();
    each_line(&args[0], |l| {
        if l.starts_with('#') {
            return Ok(());
        }
        let c: Vec<&str> = l.split('\t').collect();
        if c.len() < 4 || c[3].len() != 1 {
            return Err(bad(format!("{}: a row without a one-letter base", args[0])));
        }
        let p = usize::try_from(cell_int(c[2])?).map_err(|_| bad("a negative position"))?;
        bases
            .entry(c[1].to_string())
            .or_default()
            .insert(p, c[3].as_bytes()[0]);
        Ok(())
    })?;
    let order = (1..=22)
        .map(|i| i.to_string())
        .chain(["X", "Y", "MT"].map(String::from));
    let (mut fa, mut fai) = (create(&args[1])?, create(&format!("{}.fai", args[1]))?);
    let mut off = 0usize;
    for chrom in order {
        let Some(sites) = bases.get(&chrom) else {
            continue;
        };
        let n = sites.keys().next_back().map_or(0, |m| m + 1);
        let head = format!(">{chrom}\n");
        fa.write_all(head.as_bytes())?;
        off += head.len();
        writeln!(fai, "{chrom}\t{n}\t{off}\t60\t61")?;
        let mut seq = vec![b'N'; n];
        for (&p, &b) in sites {
            if p == 0 {
                return Err(bad("a site at position 0"));
            }
            seq[p - 1] = b;
        }
        for chunk in seq.chunks(60) {
            fa.write_all(chunk)?;
            fa.write_all(b"\n")?;
        }
        off += n + n.div_ceil(60);
    }
    fa.flush()?;
    fai.flush()?;
    Ok(())
}

/// The allele indexes of a GT (`/` or `|`), without missing ones.
fn gt_indexes(gt: &str) -> Vec<&str> {
    gt.split(['/', '|'])
        .filter(|i| !i.is_empty() && *i != ".")
        .collect()
}

/// The first letters of the alleles a GT names.
fn first_letters(alleles: &[&str], idx: &[&str]) -> Result<BTreeSet<String>, Fail> {
    idx.iter()
        .map(|i| {
            let a = i
                .parse::<usize>()
                .ok()
                .and_then(|i| alleles.get(i))
                .ok_or_else(|| bad(format!("allele index '{i}' out of range")))?;
            Ok(a.chars().next().map(String::from).unwrap_or_default())
        })
        .collect()
}

/// Python's `repr(sorted(set))` of plain strings: `['A', 'G']`.
fn repr(s: &BTreeSet<String>) -> String {
    let items: Vec<String> = s.iter().map(|x| format!("'{x}'")).collect();
    format!("[{}]", items.join(", "))
}

/// `og dev fixture mt-roundtrip LIFTED CALLS SITES TXT`: `called got
/// on_place bad` for a chip file read back by chip-to-vcf.sh. called: the
/// chrM sites of the table (SITES) TXT has a call for; got: the chrM
/// records of LIFTED (`ID CHROM POS REF ALT GT`) with a call; on_place:
/// those on the table's rCRS position with its reference base; bad: those
/// of them whose alleles are not the case's own call there (CALLS, `CHROM
/// POS REF ALT FILTER GT`, PASS only; no call means the reference base).
pub fn mt_roundtrip(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> Result<(), Fail> {
    want(args, 4)?;
    let (lifted, calls, sites, txt) = (&args[0], &args[1], &args[2], &args[3]);
    let mut table: HashMap<String, Vec<String>> = HashMap::new();
    each_line(sites, |l| {
        if !l.starts_with('#') {
            let c: Vec<String> = l.split('\t').map(str::to_string).collect();
            if c.get(1).is_some_and(|x| x == "MT") {
                if c.len() < 6 {
                    return Err(bad(format!("{sites}: an MT row with fewer than 6 columns")));
                }
                table.insert(c[0].clone(), c);
            }
        }
        Ok(())
    })?;
    let mut dv: HashMap<i64, BTreeSet<String>> = HashMap::new();
    each_line(calls, |l| {
        let c: Vec<&str> = l.split('\t').collect();
        let [_, pos, r, alt, flt, gt] = c[..] else {
            return Err(bad(format!("{calls}: a row without 6 columns")));
        };
        if flt == "PASS" {
            let al: Vec<&str> = std::iter::once(r).chain(alt.split(',')).collect();
            dv.insert(cell_int(pos)?, first_letters(&al, &gt_indexes(gt))?);
        }
        Ok(())
    })?;
    let mut called = 0u64;
    each_line(txt, |l| {
        if l.starts_with('#') {
            return Ok(());
        }
        let c: Vec<&str> = l.split('\t').collect();
        if table.contains_key(c[0]) {
            let g = c
                .get(3)
                .ok_or_else(|| bad(format!("{txt}: a row with fewer than 4 columns")))?;
            called += u64::from(*g != "--");
        }
        Ok(())
    })?;
    let (mut got, mut on_place, mut bad_n) = (0u64, 0u64, 0u64);
    each_line(lifted, |l| {
        let c: Vec<&str> = l.split('\t').collect();
        let [rid, _, pos, r, alt, gt] = c[..] else {
            return Err(bad(format!("{lifted}: a row without 6 columns")));
        };
        let idx = gt_indexes(gt);
        let Some(t) = table.get(rid) else {
            return Ok(());
        };
        if idx.is_empty() {
            return Ok(());
        }
        got += 1;
        if cell_int(pos)? != cell_int(&t[4])? || r != t[5] {
            if got - on_place <= 5 {
                writeln!(
                    err,
                    "off place {rid}: rCRS {} {}, read back at chrM:{pos} {r}",
                    t[4], t[5]
                )?;
            }
            return Ok(());
        }
        on_place += 1;
        let al: Vec<&str> = std::iter::once(r).chain(alt.split(',')).collect();
        let have = first_letters(&al, &idx)?;
        let want = dv
            .get(&cell_int(pos)?)
            .cloned()
            .unwrap_or_else(|| BTreeSet::from([t[5].clone()]));
        if have != want {
            bad_n += 1;
            if bad_n <= 5 {
                writeln!(
                    err,
                    "mismatch {rid} chrM:{pos} read back {}, case 21 {}",
                    repr(&have),
                    repr(&want)
                )?;
            }
        }
        Ok(())
    })?;
    writeln!(out, "{called} {got} {on_place} {bad_n}")?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::super::testutil::{dir, run};

    #[test]
    fn plant_snv_plants_every_nth_c_and_cleans_the_header() {
        let d = dir("plant");
        std::fs::write(
            d.join("ref.fa.fai"),
            "chr2\t1000\t6\t60\t61\nchr9\t500\t1030\t60\t61\n",
        )
        .expect("fai");
        let mut sam = String::new();
        // Each read covers 105 with a C there (offset 5 of the M after a 2S).
        for i in 0..5 {
            sam.push_str(&format!("r{i}\t0\tchr2\t100\t60\t2S8M\tchr7\t300\t0\tAAAAAAACAA\tIIIIIIIIII\tRG:Z:old\tNM:i:0\n"));
        }
        sam.push_str("x\t0\tchr7\t100\t60\t5M\t=\t1\t0\tAAAAA\tIIIII\n");
        sam.push_str("y\t0\tchr2\t104\t60\t1M1D3M\tchr9\t1\t0\tAAAA\tIIII\n");
        std::fs::write(d.join("in.sam"), sam).expect("sam");
        let p = |n: &str| d.join(n).to_str().expect("utf8").to_string();
        let (code, out, err) = run("plant-snv", &["2", &p("ref.fa.fai"), &p("in.sam"), "105"]);
        assert_eq!(code, 0, "{err}");
        let l: Vec<&str> = out.lines().collect();
        assert_eq!(
            &l[..4],
            [
                "@HD\tVN:1.6\tSO:unsorted",
                "@SQ\tSN:chr2\tLN:1000",
                "@SQ\tSN:chr9\tLN:500",
                "@RG\tID:giab\tSM:HG002\tPL:ILLUMINA"
            ]
        );
        assert_eq!(l.len(), 10, "{out}");
        assert_eq!(
            l[4],
            "r0\t0\tchr2\t100\t60\t2S8M\t*\t0\t0\tAAAAAAATAA\tIIIIIIIIII\tNM:i:0\tRG:Z:giab"
        );
        assert!(
            l[5].contains("\tAAAAAAACAA\t") && l[6].contains("\tAAAAAAATAA\t"),
            "{out}"
        );
        assert_eq!(
            l[9],
            "y\t0\tchr2\t104\t60\t1M1D3M\tchr9\t1\t0\tAAAA\tIIII\tRG:Z:giab"
        );
        let (_, none, _) = run("plant-snv", &["0", &p("ref.fa.fai"), &p("in.sam"), "105"]);
        assert!(!none.contains("AAAAAAATAA"), "{none}");
        std::fs::write(d.join("planted.sam"), &out).expect("sam");
        // r0, r2 and r4 carry the T; y covers 105 with a deletion, so not counted.
        let (code, count, _) = run("base-count", &["105", "T", &p("planted.sam")]);
        assert_eq!((code, count.as_str()), (0, "3 5\n"));
        assert_eq!(
            run("base-count", &["105", "C", &p("planted.sam")]).1,
            "2 5\n"
        );
        assert_eq!(run("base-count", &["105", "TT", &p("planted.sam")]).0, 2);
    }

    #[test]
    fn sites_fasta_and_mt_roundtrip() {
        let d = dir("mt");
        std::fs::write(
            d.join("sites.tsv"),
            "#id\tchrom\tpos\tbase\trcrs_pos\trcrs_ref\nrs1\t1\t3\tA\nrs2\tMT\t2\tG\t2\tG\ni3\tMT\t70\tT\t71\tT\ni4\tMT\t5\tC\t5\tC\ni5\tMT\t9\tA\t9\tA\nrs6\tX\t61\tC\n",
        )
        .expect("sites");
        let p = |n: &str| d.join(n).to_str().expect("utf8").to_string();
        let (code, _, err) = run("sites-fasta", &[&p("sites.tsv"), &p("h.fa")]);
        assert_eq!(code, 0, "{err}");
        let fa = std::fs::read_to_string(d.join("h.fa")).expect("fa");
        assert_eq!(
            fa,
            format!(
                ">1\nNNAN\n>X\n{}\nCN\n>MT\n{}\n{}\n",
                "N".repeat(60),
                "NGNNCNNNA".to_string() + &"N".repeat(51),
                "NNNNNNNNNTN"
            )
        );
        assert_eq!(
            std::fs::read_to_string(d.join("h.fa.fai")).expect("fai"),
            "1\t4\t3\t60\t61\nX\t62\t11\t60\t61\nMT\t71\t79\t60\t61\n"
        );
        std::fs::write(
            d.join("calls.tsv"),
            "chrM\t5\tC\tT,G\tPASS\t1/2\nchrM\t9\tA\tG\tLowQual\t1/1\n",
        )
        .expect("calls");
        std::fs::write(
            d.join("chip.txt"),
            "# head\nrs2\tMT\t2\tGG\ni3\tMT\t70\t--\ni4\tMT\t5\tTG\ni5\tMT\t9\tAA\nrs1\t1\t3\tAA\n",
        )
        .expect("txt");
        std::fs::write(
            d.join("lifted.tsv"),
            "rs2\tchrM\t2\tG\t.\t0/0\ni4\tchrM\t5\tC\tT\t0/1\ni5\tchrM\t9\tA\tG\t1/1\ni3\tchrM\t70\tT\t.\t./.\nzz\tchrM\t1\tA\t.\t0\n",
        )
        .expect("lifted");
        let (code, out, err) = run(
            "mt-roundtrip",
            &[
                &p("lifted.tsv"),
                &p("calls.tsv"),
                &p("sites.tsv"),
                &p("chip.txt"),
            ],
        );
        assert_eq!(code, 0, "{err}");
        assert_eq!(out, "3 3 3 2\n");
        assert_eq!(
            err,
            "mismatch i4 chrM:5 read back ['C', 'T'], case 21 ['G', 'T']\nmismatch i5 chrM:9 read back ['G'], case 21 ['A']\n"
        );
    }
}
