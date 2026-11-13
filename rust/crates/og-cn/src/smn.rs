//! The one rule that reads SMN1's copy number at exon 7 from step 35's SMN
//! table (Parascopy's `res.samples.bed`, as `chrom start end locus
//! agCN_filter agCN agCN_qual psCN_filter psCN psCN_qual
//! homologous_regions`).
//!
//! Rows in file order; each row covers its own region and its homologous
//! regions. The first row with a region that holds exon 7 gives the psCN of
//! that region, and the first row with a region that holds a copy of the
//! SMN1 g.27134T>G site gives the agCN. A copy number counts only with
//! filter PASS, quality 20 or more and a whole number. This is the rule of
//! `paraphase_report.parascopy_smn1`; the SMN1 marker table took the last
//! matching row before, so a table with two rows over exon 7 made the
//! marker and the Paraphase comparison read different rows.

use crate::Error;

/// SMN1 exon 7 (c.840, chr5:70,951,946).
pub const SMN1_EXON7: (&str, i64) = ("chr5", 70_951_946);

/// What the SMN table says at SMN1 exon 7 and at the marker site.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SmnRead {
    /// SMN1 psCN at exon 7, when resolved.
    pub pscn: Option<i64>,
    /// How the psCN reads, for the table.
    pub text: String,
    /// The aggregate SMN1 + SMN2 copies at the marker site, when resolved.
    pub agcn: Option<i64>,
}

/// One value of a `psCN`/`agCN` column (comma separated, one per region).
fn value(row: &[(&str, &str)], kind: &str, idx: usize) -> (Option<i64>, String) {
    let get = |k: &str| {
        row.iter()
            .rev()
            .find(|(h, _)| *h == k)
            .map(|(_, v)| *v)
            .unwrap_or("")
    };
    let vals: Vec<&str> = get(kind).split(',').collect();
    let quals: Vec<&str> = get(&format!("{kind}_qual")).split(',').collect();
    let filt = match get(&format!("{kind}_filter")) {
        "" => ".",
        f => f,
    };
    let v = vals.get(idx).copied().unwrap_or("?");
    let q = quals.get(idx).copied().unwrap_or("*");
    let ok = filt == "PASS"
        && pycompat::num::parse_float(q).is_ok_and(|q| q >= 20.0)
        && !v.is_empty()
        && v.bytes().all(|c| c.is_ascii_digit());
    let n = if ok { v.parse().ok() } else { None };
    let text = format!(
        "{v} (quality {q}, filter {filt}){}",
        if n.is_some() { "" } else { ", not resolved" }
    );
    (n, text)
}

/// The row's region and its homologous regions, 1-based closed.
fn regions(row: &[(&str, &str)]) -> Option<Vec<(String, i64, i64)>> {
    let get = |k: &str| row.iter().rev().find(|(h, _)| *h == k).map(|(_, v)| *v);
    let start = pycompat::num::parse_int(get("start")?).ok()?;
    let end = pycompat::num::parse_int(get("end")?).ok()?;
    let mut out = vec![(
        get("chrom")?.to_owned(),
        i64::try_from(start).ok()? + 1,
        i64::try_from(end).ok()?,
    )];
    for h in get("homologous_regions").unwrap_or("").split(',') {
        if let Some(r) = homologous(crate::strip(h)) {
            out.push(r);
        }
    }
    Some(out)
}

/// `^([^:]+):(\d+)-(\d+)` at the start of `h`.
fn homologous(h: &str) -> Option<(String, i64, i64)> {
    let (chrom, rest) = h.split_once(':')?;
    if chrom.is_empty() {
        return None;
    }
    let digits = |s: &str| s.bytes().take_while(u8::is_ascii_digit).count();
    let n = digits(rest);
    let (s, rest) = rest.split_at(n);
    let rest = rest.strip_prefix('-')?;
    let e = &rest[..digits(rest)];
    if s.is_empty() || e.is_empty() {
        return None;
    }
    Some((chrom.to_owned(), s.parse().ok()?, e.parse().ok()?))
}

/// Reads step 35's SMN table at `path` (`None`, or a path that does not
/// exist: no table). `marker_copies` are the marker's sites on `chrom`.
pub fn read(path: Option<&str>, chrom: &str, marker_copies: &[i64]) -> Result<SmnRead, Error> {
    let mut out = SmnRead {
        pscn: None,
        text: "not read (no SMN table from step 35)".to_owned(),
        agcn: None,
    };
    let Some(path) = path.filter(|p| !p.is_empty() && std::path::Path::new(p).exists()) else {
        return Ok(out);
    };
    out.text = "not read (no row of the SMN table covers SMN1 exon 7)".to_owned();
    let (mut ps_done, mut ag_done) = (false, false);
    let mut header: Option<Vec<String>> = None;
    crate::each_line(path, |line| {
        let Some(h) = &header else {
            header = Some(line.split('\t').map(str::to_owned).collect());
            return Ok(());
        };
        let row: Vec<(&str, &str)> = h.iter().map(String::as_str).zip(line.split('\t')).collect();
        let Some(regs) = regions(&row) else {
            return Ok(());
        };
        if !ps_done
            && let Some(i) = regs
                .iter()
                .position(|(c, s, e)| c == SMN1_EXON7.0 && *s <= SMN1_EXON7.1 && SMN1_EXON7.1 <= *e)
        {
            (out.pscn, out.text) = value(&row, "psCN", i);
            ps_done = true;
        }
        if !ag_done
            && regs
                .iter()
                .any(|(c, s, e)| c == chrom && marker_copies.iter().any(|p| s <= p && p <= e))
        {
            out.agcn = value(&row, "agCN", 0).0;
            ag_done = true;
        }
        Ok(())
    })?;
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    const HEAD: &str = "chrom\tstart\tend\tlocus\tagCN_filter\tagCN\tagCN_qual\tpsCN_filter\tpsCN\tpsCN_qual\thomologous_regions\n";

    /// A table file that is removed, with its folder once empty, when the
    /// test ends.
    struct Table(String);

    impl Drop for Table {
        fn drop(&mut self) {
            let p = std::path::Path::new(&self.0);
            let _ = std::fs::remove_file(p);
            if let Some(dir) = p.parent() {
                let _ = std::fs::remove_dir(dir);
            }
        }
    }

    fn table(rows: &[&str]) -> Table {
        let dir =
            std::env::temp_dir().join(format!("og-cn-smn-{}-{}", std::process::id(), rows.len()));
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join(format!("t{}.tsv", rows.join("").len()));
        std::fs::write(&p, format!("{HEAD}{}", rows.join(""))).unwrap();
        Table(p.to_string_lossy().into_owned())
    }

    const SITES: [i64; 2] = [70_952_074, 70_076_654];

    #[test]
    fn the_first_row_over_exon_7_wins() {
        // Two rows over exon 7: the first is unresolved, the second says 2.
        // The first row's reading stands, as the Paraphase comparison reads it.
        let p = table(&[
            "chr5\t70951000\t70953000\tSMN1\tPASS\t4\t30\tLowQual\t1,3\t12,40\tchr5:70075000-70078000\n",
            "chr5\t70951900\t70952000\tSMN1\tPASS\t4\t30\tPASS\t2,2\t40,40\tchr5:70076000-70076100\n",
        ]);
        let r = read(Some(&p.0), "chr5", &SITES).unwrap();
        assert_eq!(r.pscn, None);
        assert_eq!(r.text, "1 (quality 12, filter LowQual), not resolved");
        assert_eq!(r.agcn, Some(4));
    }

    #[test]
    fn a_resolved_row_and_the_placeholders() {
        let p = table(&[
            "chr5\t70925000\t70955000\tSMN1\tPASS\t4\t30\tPASS\t2\t40\tchr5:70049000-70080000\n",
        ]);
        let r = read(Some(&p.0), "chr5", &SITES).unwrap();
        assert_eq!(
            (r.pscn, r.text.as_str()),
            (Some(2), "2 (quality 40, filter PASS)")
        );
        // The homologous copy holds exon 7 too, but psCN has one value.
        let p = table(&[
            "chr5\t70049000\t70050000\tSMN2\tPASS\t4\t30\tPASS\t2\t40\tchr5:70951000-70952000\n",
        ]);
        let r = read(Some(&p.0), "chr5", &SITES).unwrap();
        assert_eq!(r.text, "? (quality *, filter PASS), not resolved");
    }

    #[test]
    fn no_table_and_no_row() {
        assert_eq!(
            read(None, "chr5", &SITES).unwrap().text,
            "not read (no SMN table from step 35)"
        );
        let p = table(&["chr1\t1\t2\tX\tPASS\t2\t30\tPASS\t2\t40\t.\n"]);
        assert_eq!(
            read(Some(&p.0), "chr5", &SITES).unwrap().text,
            "not read (no row of the SMN table covers SMN1 exon 7)"
        );
    }
}
