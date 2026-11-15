//! Unmapped viral reads to plant in a BAM, for step 47.

use super::{Fail, PyRandom, bad, each_line, want};
use std::collections::HashMap;
use std::io::Write;

/// The read length.
const L: usize = 150;

/// `og dev fixture viral-spike VIRAL.fa HUMAN_DEPTH FRACTION`: unmapped SAM
/// records (150 bp, 0.3% substitutions, random strand, seed 47): HHV-6B
/// (NC_000898.1) at FRACTION of HUMAN_DEPTH, 300 CMV reads and 200 phiX.
pub fn spike(args: &[String], out: &mut dyn Write, _err: &mut dyn Write) -> Result<(), Fail> {
    want(args, 3)?;
    let num = |s: &str, what: &str| -> Result<f64, Fail> {
        pycompat::num::parse_float(s)
            .map_err(|_| Fail::Usage(format!("{what} '{s}' is not a number")))
    };
    let (human, frac) = (num(&args[1], "HUMAN_DEPTH")?, num(&args[2], "FRACTION")?);
    let mut seqs: HashMap<String, Vec<u8>> = HashMap::new();
    let mut name: Option<String> = None;
    each_line(&args[0], |l| {
        let l = l.trim();
        if let Some(h) = l.strip_prefix('>') {
            let n = h.split_whitespace().next().unwrap_or("").to_string();
            seqs.insert(n.clone(), Vec::new());
            name = Some(n);
        } else if !l.is_empty() {
            let n = name
                .as_ref()
                .ok_or_else(|| bad("sequence before the first header"))?;
            seqs.get_mut(n)
                .expect("inserted with its header")
                .extend(l.bytes().map(|b| b.to_ascii_uppercase()));
        }
        Ok(())
    })?;
    let mut rng = PyRandom::new(47);
    let mut reads = |acc: &str, n: u64, tag: &str, out: &mut dyn Write| -> Result<(), Fail> {
        let s = seqs
            .get(acc)
            .ok_or_else(|| bad(format!("{}: no {acc}", args[0])))?;
        if s.len() <= L {
            return Err(bad(format!("{acc} is shorter than a read")));
        }
        let qual = "I".repeat(L);
        for i in 0..n {
            let p = rng.below((s.len() - L) as u64) as usize;
            let mut r = s[p..p + L].to_vec();
            for b in r.iter_mut() {
                if rng.random() < 0.003 {
                    let other: Vec<u8> = b"ACGT".iter().copied().filter(|x| x != b).collect();
                    *b = other[rng.below(other.len() as u64) as usize];
                }
            }
            if rng.random() < 0.5 {
                r = r
                    .iter()
                    .rev()
                    .map(|&b| match b {
                        b'A' => b'T',
                        b'C' => b'G',
                        b'G' => b'C',
                        b'T' => b'A',
                        x => x,
                    })
                    .collect();
            }
            let r = String::from_utf8_lossy(&r);
            writeln!(out, "{tag}_{i}\t4\t*\t0\t0\t*\t*\t0\t0\t{r}\t{qual}")?;
        }
        Ok(())
    };
    let hhv6b = seqs.get("NC_000898.1").map_or(0, Vec::len) as f64;
    let n = pycompat::fmt::round_int(frac * human * hhv6b / L as f64)
        .map_err(|_| bad("the HHV-6B read count is not finite"))?;
    reads("NC_000898.1", n as u64, "spk6b", out)?;
    reads("NC_006273.2", 300, "spkcmv", out)?;
    reads("NC_001422.1", 200, "spkphix", out)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::super::testutil::{dir, fnv, run};

    #[test]
    fn spike_is_the_pythons_bytes() {
        let d = dir("viral");
        let mut rng = super::PyRandom::new(3);
        let mut fa = String::new();
        for (acc, n) in [
            ("NC_000898.1 HHV-6B", 3000),
            ("NC_006273.2", 900),
            ("NC_001664.4", 400),
            ("NC_001422.1", 700),
        ] {
            fa.push_str(&format!(">{acc}\n"));
            let seq: String = (0..n)
                .map(|_| ['A', 'c', 'G', 't', 'N'][rng.below(5) as usize])
                .collect();
            for chunk in seq.as_bytes().chunks(70) {
                fa.push_str(std::str::from_utf8(chunk).expect("ascii"));
                fa.push('\n');
            }
        }
        std::fs::write(d.join("v.fa"), fa).expect("fa");
        let (code, out, err) = run(
            "viral-spike",
            &[d.join("v.fa").to_str().expect("utf8"), "31.27", "0.5"],
        );
        assert_eq!(code, 0, "{err}");
        // round(0.5 * 31.27 * 3000 / 150) = round(312.7) = 313 HHV-6B reads.
        assert_eq!(out.lines().filter(|l| l.starts_with("spk6b_")).count(), 313);
        assert_eq!(out.lines().count(), 813);
        // FNV-1a of what the Python heredoc of viral-reads-1.sh wrote.
        assert_eq!(
            fnv(out.as_bytes()),
            0x085f_613c_1fff_c640,
            "{:#x}",
            fnv(out.as_bytes())
        );
    }
}
