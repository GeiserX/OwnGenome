//! Small, local smoke fixtures and table checks. No subprocess or network.
use og_io::{AtomicFile, Lines, Utf8};
use pycompat::csv::{DictReader, Quoting, Reader};
use std::io::Write;

fn lines(path: &str) -> Result<Lines<Box<dyn std::io::BufRead>>, String> {
    Ok(Lines::new(
        og_io::open(path).map_err(|e| format!("{path}: {e}"))?,
        Utf8::Strict,
    ))
}
pub(crate) fn csv_lines(
    path: &str,
) -> Result<impl Iterator<Item = std::io::Result<String>>, String> {
    let mut input = lines(path)?;
    Ok(std::iter::from_fn(move || {
        let mut line = String::new();
        match input.read_line(&mut line) {
            Ok(true) => Some(Ok(line)),
            Ok(false) => None,
            Err(e) => Some(Err(e)),
        }
    }))
}
fn fasta(path: &str) -> Result<Vec<u8>, String> {
    let mut seq = Vec::new();
    for line in lines(path)? {
        let line = line.map_err(|e| e.to_string())?;
        if !line.starts_with('>') {
            if seq.len().saturating_add(line.trim().len()) > 1_003_000 {
                return Err("smoke reference exceeds 1003000 bases".into());
            }
            seq.extend(line.trim().bytes().map(|b| b.to_ascii_uppercase()));
        }
    }
    if seq.len() > 1_003_000 {
        return Err("smoke reference exceeds 1003000 bases".into());
    }
    if !seq.is_ascii() {
        return Err("FASTA sequence is not ASCII".into());
    }
    Ok(seq)
}
fn write_fasta(path: &str, seq: &[u8]) -> Result<(), String> {
    let mut out = AtomicFile::create(path).map_err(|e| e.to_string())?;
    out.write_all(b">chr20\n").map_err(|e| e.to_string())?;
    for chunk in seq.chunks(60) {
        out.write_all(chunk)
            .and_then(|()| out.write_all(b"\n"))
            .map_err(|e| e.to_string())?;
    }
    out.commit().map_err(|e| e.to_string())
}
fn number(s: &str) -> Result<usize, String> {
    s.parse().map_err(|_| format!("invalid integer: {s}"))
}
fn mini(args: &[String]) -> Result<(), String> {
    if args.len() != 5 {
        return Err("mini-ref SRC ORIG MINI POS LEN".into());
    }
    let mut seq = fasta(&args[0])?;
    if seq.len() != 500_000 {
        return Err(format!("mini-ref needs 500000 bases, got {}", seq.len()));
    }
    let pos = number(&args[3])?;
    let n = number(&args[4])?;
    if pos > seq.len() || n > 1_000_000 {
        return Err("plant outside bounded mini reference".into());
    }
    write_fasta(&args[1], &seq)?;
    // Stable, portable xorshift64 seed. These bases are deliberately independent
    // of CPython's RNG. They only create a deletion in reads of the original.
    let mut state = 20u64;
    let insert = (0..n).map(|_| {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        b"ACGT"[(state & 3) as usize]
    });
    seq.splice(pos..pos, insert);
    write_fasta(&args[2], &seq)
}
fn long_reads(args: &[String], stdout: &mut dyn Write) -> Result<(), String> {
    if args.len() != 3 {
        return Err("long-reads HAP1 HAP2 OUT.gz".into());
    }
    let mut out = AtomicFile::create(&args[2]).map_err(|e| e.to_string())?;
    let mut gz = flate2::GzBuilder::new()
        .mtime(0)
        .write(&mut out, flate2::Compression::default());
    let mut k = 0;
    for (h, path) in args[..2].iter().enumerate() {
        let seq = fasta(path)?;
        if seq.len() < 10000 {
            continue;
        }
        for start in (0..=seq.len() - 10000).step_by(1000) {
            let mut r = seq[start..start + 10000].to_vec();
            if k % 2 == 1 {
                r.reverse();
                for b in &mut r {
                    *b = match b {
                        b'A' => b'T',
                        b'C' => b'G',
                        b'G' => b'C',
                        b'T' => b'A',
                        _ => *b,
                    };
                }
            }
            writeln!(gz, "@hap{}_{start}", h + 1)
                .and_then(|()| gz.write_all(&r))
                .and_then(|()| gz.write_all(b"\n+\n"))
                .and_then(|()| gz.write_all(&[b'?'; 10000]))
                .and_then(|()| gz.write_all(b"\n"))
                .map_err(|e| e.to_string())?;
            k += 1;
        }
    }
    gz.finish().map_err(|e| e.to_string())?;
    out.commit().map_err(|e| e.to_string())?;
    writeln!(stdout, "{k} reads").map_err(|e| e.to_string())
}
fn catalog(args: &[String], stdout: &mut dyn Write) -> Result<(), String> {
    if args.len() != 3 {
        return Err("eh-catalog REF OUT END".into());
    }
    let seq = fasta(&args[0])?;
    let end = number(&args[2])?;
    let re = regex::bytes::Regex::new("(?:CA){12,}").expect("constant regex");
    let m = re
        .find_iter(&seq)
        .find(|m| m.start() > 20000 && m.end().checked_add(20000).is_some_and(|n| n < end))
        .ok_or("no suitable CA repeat")?;
    let data = format!(
        "[\n {{\n  \"LocusId\": \"SMOKE1\",\n  \"LocusStructure\": \"(CA)*\",\n  \"ReferenceRegion\": \"chr20:{}-{}\",\n  \"VariantType\": \"Repeat\"\n }}\n]",
        m.start(),
        m.end()
    );
    let mut out = AtomicFile::create(&args[1]).map_err(|e| e.to_string())?;
    out.write_all(data.as_bytes()).map_err(|e| e.to_string())?;
    out.commit().map_err(|e| e.to_string())?;
    writeln!(
        stdout,
        "repeat {} {} {} units",
        m.start(),
        m.end(),
        m.len() / 2
    )
    .map_err(|e| e.to_string())
}
fn csv_cell(args: &[String], out: &mut dyn Write) -> Result<(), String> {
    if args.len() != 3 {
        return Err("csv-cell FILE KEY=VALUE[,KEY=VALUE] COLUMN".into());
    }
    let want = args[1]
        .split(',')
        .map(|kv| kv.split_once('=').ok_or("expected KEY=VALUE"))
        .collect::<Result<Vec<_>, _>>()?;
    for row in DictReader::new(csv_lines(&args[0])?, ',', Quoting::Minimal) {
        let row = row.map_err(|e| e.to_string())?;
        if want.iter().all(|(k, v)| row.get(k) == Some(Some(v))) {
            let val = row
                .get(&args[2])
                .ok_or("missing CSV column")?
                .unwrap_or("None");
            writeln!(out, "{val}").map_err(|e| e.to_string())?;
            break;
        }
    }
    Ok(())
}
fn fields(line: &str) -> Vec<&str> {
    line.split_whitespace().collect()
}
fn float(s: Option<&&str>) -> Result<f64, String> {
    pycompat::num::parse_float(s.ok_or("short table row")?).map_err(|e| e.to_string())
}
#[derive(Default)]
struct PcaGroup {
    n: usize,
    first: f64,
    min: f64,
    max: f64,
    sum: f64,
    correction: f64,
}
impl PcaGroup {
    fn add(&mut self, x: f64) {
        if self.n == 0 {
            self.first = x;
            self.min = x;
            self.max = x;
        }
        // Python's min/max keep the first element when a comparison is
        // unordered (NaN), unlike f64::min/max, which skip NaN.
        if x < self.min {
            self.min = x;
        }
        if x > self.max {
            self.max = x;
        }
        let t = self.sum + x;
        self.correction += if self.sum.abs() >= x.abs() {
            (self.sum - t) + x
        } else {
            (x - t) + self.sum
        };
        self.sum = t;
        self.n += 1;
    }
    fn mean(&self) -> f64 {
        (self.sum
            + if self.correction != 0. && self.correction.is_finite() {
                self.correction
            } else {
                0.
            })
            / self.n as f64
    }
}
fn parascopy(path: &str, out: &mut dyn Write) -> Result<bool, String> {
    let mut found = false;
    let mut n = 0;
    for row in Reader::new(csv_lines(path)?, '\t', Quoting::Minimal) {
        let row = row.map_err(|e| e.to_string())?;
        if row.is_empty() || row[0].starts_with('#') {
            continue;
        }
        if row.len() < 8 {
            return Err("short Parascopy row".into());
        }
        found |= row[0] == "chr5"
            && !["", "*"].contains(&row[6].as_str())
            && !["", "*"].contains(&row[7].as_str());
        n += 1;
    }
    writeln!(
        out,
        "parascopy: {n} rows; {}",
        if found { "OK" } else { "FAIL" }
    )
    .map_err(|e| e.to_string())?;
    Ok(found)
}
fn check(kind: &str, path: &str, out: &mut dyn Write) -> Result<bool, String> {
    if kind == "parascopy" {
        return parascopy(path, out);
    }
    // The windows keep only their first covering row, matching Python's next().
    let mut n = 0;
    let mut ab = 0;
    let mut ab_len = 0.;
    let mut coverage = 0.;
    let mut outside = 0;
    let mut groups = [
        PcaGroup::default(),
        PcaGroup::default(),
        PcaGroup::default(),
    ];
    let mut codes = std::collections::BTreeMap::new();
    let mut at = [None, None];
    for (i, line) in lines(path)?.enumerate() {
        let line = line.map_err(|e| e.to_string())?;
        let r = fields(&line);
        match kind {
            "hap-ibd" => {
                if r.len() < 8 {
                    return Err("short hap-ibd row".into());
                }
                n += 1;
                if (r[0] == "A" && r[2] == "B") || (r[0] == "B" && r[2] == "A") {
                    ab += 1;
                    ab_len = float(r.get(7))?;
                }
            }
            "ibdmix" => {
                if i == 0 {
                    continue;
                }
                if r.first() == Some(&"S")
                    && float(r.get(4))? >= 4.
                    && float(r.get(3))? - float(r.get(2))? >= 50000.
                {
                    let a = float(r.get(2))?;
                    let b = float(r.get(3))?;
                    n += 1;
                    coverage += (b.min(2600000.) - a.max(2000000.)).max(0.);
                    if b <= 2000000. || a >= 2600000. {
                        outside += 1;
                    }
                }
            }
            "eigen" => {
                if line.trim_start().starts_with('#') || r.is_empty() {
                    continue;
                }
                let g = match r.last() {
                    Some(&"PopA") => 0,
                    Some(&"PopB") => 1,
                    Some(&"Proj") => 2,
                    _ => continue,
                };
                groups[g].add(float(r.get(1))?);
            }
            "rfmix" => {
                if i == 0 {
                    for pair in line
                        .split_once(':')
                        .ok_or("missing population codes")?
                        .1
                        .split_whitespace()
                    {
                        let (k, v) = pair.split_once('=').ok_or("bad population code")?;
                        codes.insert(k.to_owned(), v.to_owned());
                    }
                }
                if line.starts_with('#') {
                    continue;
                }
                let a = float(r.get(1))?;
                let b = float(r.get(2))?;
                n += 1;
                for (j, p) in [5000000., 15000000.].iter().enumerate() {
                    if at[j].is_none() && a <= *p && *p < b {
                        at[j] = Some([
                            r.get(6).ok_or("short RFMix row")?.to_string(),
                            r.get(7).ok_or("short RFMix row")?.to_string(),
                        ]);
                    }
                }
            }
            _ => return Err(format!("unknown smoke check {kind}")),
        }
    }
    let ok = match kind {
        "hap-ibd" => n == 1 && ab == 1 && (8. ..=11.).contains(&ab_len),
        "ibdmix" => coverage >= 540000. && outside == 0,
        "eigen" => {
            if groups.iter().any(|g| g.n == 0) {
                return Err("missing PCA group".into());
            }
            let a = &groups[0];
            let b = &groups[1];
            let p = groups[2].first;
            (a.max < b.min || a.min > b.max) && (p - a.mean()).abs() < (p - b.mean()).abs()
        }
        "rfmix" => ["POPA", "POPB"].iter().enumerate().all(|(i, p)| {
            codes
                .get(*p)
                .is_some_and(|c| at[i].as_ref() == Some(&[c.clone(), c.clone()]))
        }),
        _ => false,
    };
    writeln!(out, "{kind}: {n} rows; {}", if ok { "OK" } else { "FAIL" })
        .map_err(|e| e.to_string())?;
    Ok(ok)
}

/// Runs a smoke generator, CSV lookup or table assertion.
pub fn cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let result = (|| {
        let (verb, args) = args.split_first().ok_or("smoke verb required")?;
        match verb.as_str() {
            "mini-ref" => mini(args).map(|()| true),
            "long-reads" => long_reads(args, out).map(|()| true),
            "eh-catalog" => catalog(args, out).map(|()| true),
            "csv-cell" => csv_cell(args, out).map(|()| true),
            "check" if args.len() == 2 => check(&args[0], &args[1], out),
            _ => {
                Err("usage: og dev smoke mini-ref|long-reads|eh-catalog|csv-cell|check ...".into())
            }
        }
    })();
    match result {
        Ok(true) => 0,
        Ok(false) => 1,
        Err(e) => {
            let _ = writeln!(err, "ERROR: {e}");
            2
        }
    }
}
