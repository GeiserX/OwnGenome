//! The always-report list (`og clinvar deep-intronic`): the sample's call
//! at every site of a list of known pathogenic alleles, carrier or not.
//!
//! The list (`assets/deep_intronic_sites.tsv`) has `#` comment lines, a
//! header row and the columns `chrom pos ref alt gene hgvs clinvar_id`. The
//! calls are `bcftools query
//! '%CHROM\t%POS\t%END\t%REF\t%ALT\t[%GT]\t[%GQ]\t[%DP]\n'` of the split PASS
//! VCF and of the gVCF at those positions. One row per site comes out, with
//! two columns more:
//!
//! - `call`: the genotype in words when the sample carries the site's
//!   allele; `ref` for a reference call (a `0/0` record, or a gVCF reference
//!   block over the site), `ref (low confidence)` when its GQ is below the
//!   floor; `other allele` when the sample carries another allele there;
//!   `no call` otherwise.
//! - `detail`: the call's `GQ n, DP n` (`GQ n (gVCF block)` for a block), the
//!   floor a low-confidence reference call is below, the other allele, or
//!   why there is no call.

use std::io::{self, Write};

use crate::text::{is_digits, nth, same, split};
use crate::zygosity;

/// The floor and whether there is a gVCF to read reference blocks from.
pub struct Options {
    /// The GQ floor as given, for the words.
    pub min_gq: String,
    /// The GQ floor.
    pub min_gq_value: u64,
    pub gvcf: bool,
}

/// The number at the start of a field (after blanks), 0 when there is none.
fn leading_number(s: &[u8]) -> f64 {
    let t = s.trim_ascii_start();
    let n = t.iter().take_while(|c| c.is_ascii_digit()).count();
    std::str::from_utf8(&t[..n])
        .ok()
        .and_then(|d| d.parse().ok())
        .unwrap_or(0.0)
}

fn carries(gt: &[u8]) -> bool {
    gt.iter().any(|c| (b'1'..=b'9').contains(c))
}

fn is_ref(gt: &[u8]) -> bool {
    matches!(gt, [b'0', b'/' | b'|', b'0'])
}

impl Options {
    fn below(&self, gq: &[u8]) -> bool {
        is_digits(gq) && leading_number(gq) < self.min_gq_value as f64
    }

    fn refcall(&self, gq: &[u8]) -> &'static [u8] {
        if self.below(gq) {
            b"ref (low confidence)"
        } else {
            b"ref"
        }
    }

    fn floor_note(&self, gq: &[u8]) -> Vec<u8> {
        if self.below(gq) {
            [b", below the call floor of GQ ", self.min_gq.as_bytes()].concat()
        } else {
            Vec::new()
        }
    }
}

/// `GQ n, DP n`, each when it is a whole number.
fn quality(gq: &[u8], dp: &[u8]) -> Vec<u8> {
    let mut out = Vec::new();
    if is_digits(gq) {
        out.extend([b"GQ ", gq].concat());
    }
    if is_digits(dp) {
        if !out.is_empty() {
            out.extend_from_slice(b", ");
        }
        out.extend([b"DP ", dp].concat());
    }
    out
}

/// The lines of a text, without line endings.
fn lines(text: &[u8]) -> Vec<&[u8]> {
    let mut v: Vec<&[u8]> = text.split(|&c| c == b'\n').collect();
    if text.is_empty() || text.ends_with(b"\n") {
        v.pop();
    }
    v
}

/// The call at one site.
fn call(st: &[&[u8]], vcf: &[Vec<&[u8]>], gvcf: &[Vec<&[u8]>], o: &Options) -> (Vec<u8>, Vec<u8>) {
    let (mut call, mut detail, mut other) = (Vec::new(), Vec::new(), Vec::new());
    for vr in vcf {
        if !same(nth(vr, 1), nth(st, 1)) || !same(nth(vr, 2), nth(st, 2)) {
            continue;
        }
        let gt = nth(vr, 6);
        if same(nth(vr, 4), nth(st, 3)) && same(nth(vr, 5), nth(st, 4)) {
            if carries(gt) {
                call = zygosity(gt);
                detail = quality(nth(vr, 7), nth(vr, 8));
            } else if is_ref(gt) {
                call = o.refcall(nth(vr, 7)).to_vec();
                detail = [quality(nth(vr, 7), nth(vr, 8)), o.floor_note(nth(vr, 7))].concat();
            }
        } else if carries(gt) && other.is_empty() {
            other = [nth(vr, 4), b">", nth(vr, 5), b" ", &zygosity(gt)].concat();
        }
    }
    if call.is_empty() && !other.is_empty() {
        return (b"other allele".to_vec(), other);
    }
    if call.is_empty() {
        let pos = leading_number(nth(st, 2));
        let over = gvcf.iter().find(|gr| {
            same(nth(gr, 1), nth(st, 1))
                && leading_number(nth(gr, 2)) <= pos
                && leading_number(nth(gr, 3)) >= pos
        });
        if let Some(gr) = over {
            let (gt, gq) = (nth(gr, 6), nth(gr, 7));
            if is_ref(gt) {
                let sep: &[u8] = if is_digits(gq) { b" " } else { b"" };
                return (
                    o.refcall(gq).to_vec(),
                    [
                        quality(gq, b""),
                        sep.to_vec(),
                        b"(gVCF block)".to_vec(),
                        o.floor_note(gq),
                    ]
                    .concat(),
                );
            }
            return (
                b"no call".to_vec(),
                [
                    &b"the gVCF record over the site is not a PASS reference call (GT "[..],
                    gt,
                    b")",
                ]
                .concat(),
            );
        }
        let why: &[u8] = if o.gvcf {
            b"no record in the VCF or the gVCF covers the site"
        } else {
            b"no record in the VCF, and no gVCF to tell a reference call from no call"
        };
        return (b"no call".to_vec(), why.to_vec());
    }
    (call, detail)
}

/// Writes the header and one row per site of `sites`, in its order.
pub fn deep_intronic(
    sites: &[u8],
    vcf: &[u8],
    gvcf: &[u8],
    o: &Options,
    out: &mut impl Write,
) -> io::Result<()> {
    let table = |t| {
        lines(t)
            .into_iter()
            .map(|l| split(l, b'\t'))
            .collect::<Vec<_>>()
    };
    let (vcf, gvcf) = (table(vcf), table(gvcf));
    out.write_all(b"chrom\tpos\tref\talt\tgene\thgvs\tclinvar_id\tcall\tdetail\n")?;
    for line in lines(sites) {
        let st = split(line, b'\t');
        if line.starts_with(b"#") || nth(&st, 1) == b"chrom" || st.len() < 7 {
            continue;
        }
        let (c, d) = call(&st, &vcf, &gvcf, o);
        let cols: Vec<&[u8]> = (1..=7)
            .map(|k| nth(&st, k))
            .chain([&c[..], &d[..]])
            .collect();
        out.write_all(&cols.join(&b'\t'))?;
        out.write_all(b"\n")?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run(vcf: &str, gvcf: &str, gvcf_on: bool) -> String {
        let sites = "# test\nchrom\tpos\tref\talt\tgene\thgvs\tclinvar_id\n\
                     chr1\t100\tA\tG\tG1\tc.1A>G\t1\n";
        let o = Options {
            min_gq: "30".into(),
            min_gq_value: 30,
            gvcf: gvcf_on,
        };
        let mut out = Vec::new();
        deep_intronic(
            sites.as_bytes(),
            vcf.as_bytes(),
            gvcf.as_bytes(),
            &o,
            &mut out,
        )
        .unwrap();
        let s = String::from_utf8(out).unwrap();
        s.lines()
            .nth(1)
            .unwrap()
            .split('\t')
            .skip(7)
            .collect::<Vec<_>>()
            .join("|")
    }

    #[test]
    fn every_kind_of_call() {
        assert_eq!(
            run("chr1\t100\t100\tA\tG\t0/1\t40\t30\n", "", true),
            "het|GQ 40, DP 30"
        );
        assert_eq!(
            run("chr1\t100\t100\tA\tG\t0/0\t20\t30\n", "", true),
            "ref (low confidence)|GQ 20, DP 30, below the call floor of GQ 30"
        );
        assert_eq!(
            run("chr1\t100\t100\tA\tT\t1/1\t40\t30\n", "", true),
            "other allele|A>T hom"
        );
        assert_eq!(
            run("", "chr1\t90\t120\tA\t<*>\t0/0\t50\t.\n", true),
            "ref|GQ 50 (gVCF block)"
        );
        assert_eq!(
            run("", "chr1\t90\t120\tA\t<*>\t0/0\t.\t.\n", true),
            "ref|(gVCF block)"
        );
        assert_eq!(
            run("", "chr1\t90\t120\tA\t<*>\t./.\t.\t.\n", true),
            "no call|the gVCF record over the site is not a PASS reference call (GT ./.)"
        );
        assert_eq!(
            run("", "", true),
            "no call|no record in the VCF or the gVCF covers the site"
        );
        assert_eq!(
            run("", "", false),
            "no call|no record in the VCF, and no gVCF to tell a reference call from no call"
        );
    }
}
