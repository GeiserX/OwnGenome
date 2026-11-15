//! The always-report list (`og clinvar deep-intronic`): the sample's call
//! at every site of a list of known pathogenic alleles, carrier or not.
//!
//! The list (`assets/deep_intronic_sites.tsv`) has `#` comment lines, a
//! header row and the columns `chrom pos ref alt gene hgvs clinvar_id`. The
//! calls are `bcftools query
//! '%CHROM\t%POS\t%END\t%REF\t%ALT\t[%GT]\t[%GQ]\t[%DP]\t%FILTER\n'` of the
//! split PASS VCF and of the gVCF at those positions (an input without the
//! FILTER column reads as FILTER `.`). One row per site comes out, with two
//! columns more:
//!
//! - `call`: the genotype in words when the sample carries the site's
//!   allele (`hemizygous` for one copy on chrX or chrY outside the
//!   pseudoautosomal regions); `ref` for a reference call (a record whose
//!   alleles are all `0`, or a gVCF reference block over the site) with
//!   FILTER `PASS`, `.` or `RefCall`, a GQ at or above the call floor and a
//!   depth at or above the depth floor; `ref (low confidence)` for a
//!   reference call with another FILTER, a GQ below the call floor or no GQ;
//!   `ref (shallow)` for one whose only fault is a depth below the depth
//!   floor; `other allele` when the sample carries another allele there;
//!   `no call` otherwise.
//! - `detail`: the call's `GQ n, DP n` (`GQ n (gVCF block)` for a block,
//!   with its DP when that is below the depth floor), each reason a
//!   reference call is not plain `ref`, the other allele, or what the record
//!   at the site holds when there is no call.

use std::io::{self, Write};

use crate::text::{is_digits, nth, same, split};
use crate::zygosity;

/// The floors and whether there is a gVCF to read reference blocks from.
pub struct Options {
    /// The GQ floor as given, for the words.
    pub min_gq: String,
    /// The GQ floor.
    pub min_gq_value: u64,
    /// The depth floor as given, for the words.
    pub min_dp: String,
    /// The depth floor.
    pub min_dp_value: u64,
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

/// GRCh38's pseudoautosomal regions, 1-based and inclusive, as in
/// `assets/par_grch38.bed`.
const PAR_X: [(u64, u64); 2] = [(10_001, 2_781_479), (155_701_383, 156_030_895)];
const PAR_Y: [(u64, u64); 2] = [(10_001, 2_781_479), (56_887_903, 57_217_415)];

/// True where a caller may write one copy: chrX and chrY outside the
/// pseudoautosomal regions.
fn one_copy(chrom: &[u8], pos: u64) -> bool {
    let par = match chrom.strip_prefix(b"chr").unwrap_or(chrom) {
        b"X" => &PAR_X,
        b"Y" => &PAR_Y,
        _ => return false,
    };
    !par.iter().any(|&(a, b)| (a..=b).contains(&pos))
}

/// What a record's GT says.
#[derive(Clone, Copy, PartialEq)]
enum Gt {
    /// Every allele `0`.
    Ref,
    /// An alt allele.
    Carries,
    /// No allele called.
    NoGenotype,
    /// One allele where the genome has two.
    OneCopyWhereTwo,
    /// Anything else, such as `0/.`.
    Other,
}

impl Gt {
    fn of(gt: &[u8], one_copy_site: bool) -> Gt {
        let alleles: Vec<&[u8]> = gt
            .split(|&c| c == b'/' || c == b'|')
            .filter(|a| !a.is_empty())
            .collect();
        if alleles.iter().all(|&a| a == b".") {
            Gt::NoGenotype
        } else if alleles.len() == 1 && !one_copy_site {
            Gt::OneCopyWhereTwo
        } else if alleles.iter().all(|&a| a == b"0") {
            Gt::Ref
        } else if gt.iter().any(|c| (b'1'..=b'9').contains(c)) {
            Gt::Carries
        } else {
            Gt::Other
        }
    }
}

/// A carried GT in words: `hemizygous` for one copy, else as the hits word it.
fn words(gt: &[u8]) -> Vec<u8> {
    if gt.contains(&b'/') || gt.contains(&b'|') {
        zygosity(gt)
    } else {
        b"hemizygous".to_vec()
    }
}

/// True for a FILTER that leaves a reference call plain: `PASS`, `.`, none,
/// or DeepVariant's `RefCall` (the model reads the site as reference).
fn clean(filter: &[u8]) -> bool {
    matches!(filter, b"" | b"." | b"PASS" | b"RefCall")
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

impl Options {
    /// A reference call: its words and detail from its GQ, depth and FILTER.
    fn reference(&self, gq: &[u8], dp: &[u8], filter: &[u8], block: bool) -> (Vec<u8>, Vec<u8>) {
        let filtered = !clean(filter);
        let no_gq = !is_digits(gq);
        let low_gq = !no_gq && leading_number(gq) < self.min_gq_value as f64;
        let shallow = is_digits(dp) && leading_number(dp) < self.min_dp_value as f64;
        let call: &[u8] = if filtered || no_gq || low_gq {
            b"ref (low confidence)"
        } else if shallow {
            b"ref (shallow)"
        } else {
            b"ref"
        };
        let mut detail = quality(gq, if block && !shallow { b"" } else { dp });
        if block {
            if !detail.is_empty() {
                detail.push(b' ');
            }
            detail.extend_from_slice(b"(gVCF block)");
        }
        let mut notes: Vec<Vec<u8>> = Vec::new();
        if filtered {
            notes.push([b"FILTER ", filter].concat());
        }
        if low_gq {
            notes.push([b"below the call floor of GQ ", self.min_gq.as_bytes()].concat());
        }
        if no_gq {
            notes.push(
                [
                    b"no GQ to check against the call floor of GQ ",
                    self.min_gq.as_bytes(),
                ]
                .concat(),
            );
        }
        if shallow {
            notes.push([b"below the depth floor of DP ", self.min_dp.as_bytes()].concat());
        }
        for n in notes {
            if !detail.is_empty() {
                detail.extend_from_slice(b", ");
            }
            detail.extend(n);
        }
        (call.to_vec(), detail)
    }
}

/// Why a record (`what`) gives no call.
fn unread(what: &[u8], gt: &[u8], filter: &[u8], state: Gt) -> Vec<u8> {
    let gt: &[u8] = if gt.is_empty() { b"." } else { gt };
    match state {
        Gt::NoGenotype => [what, b" has no genotype (GT ", gt, b")"].concat(),
        Gt::OneCopyWhereTwo => [
            what,
            b" has one copy (GT ",
            gt,
            b") where the genome has two",
        ]
        .concat(),
        _ if !clean(filter) => [
            what,
            b" is not a reference call (GT ",
            gt,
            b", FILTER ",
            filter,
            b")",
        ]
        .concat(),
        _ => [what, b" is not a reference call (GT ", gt, b")"].concat(),
    }
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
    let pos = leading_number(nth(st, 2));
    let haploid = one_copy(nth(st, 1), pos as u64);
    let (mut found, mut other, mut why) = (None, Vec::new(), None);
    for vr in vcf {
        if !same(nth(vr, 1), nth(st, 1)) || !same(nth(vr, 2), nth(st, 2)) {
            continue;
        }
        let (gt, gq, dp, filter) = (nth(vr, 6), nth(vr, 7), nth(vr, 8), nth(vr, 9));
        let state = Gt::of(gt, haploid);
        if same(nth(vr, 4), nth(st, 3)) && same(nth(vr, 5), nth(st, 4)) {
            match state {
                Gt::Carries => found = Some((words(gt), quality(gq, dp))),
                Gt::Ref => found = Some(o.reference(gq, dp, filter, false)),
                _ => {
                    why.get_or_insert_with(|| {
                        unread(b"the VCF record at the site", gt, filter, state)
                    });
                }
            }
        } else if state == Gt::Carries && other.is_empty() {
            other = [nth(vr, 4), b">", nth(vr, 5), b" ", &words(gt)].concat();
        }
    }
    if let Some(f) = found {
        return f;
    }
    if !other.is_empty() {
        return (b"other allele".to_vec(), other);
    }
    let over = gvcf.iter().find(|gr| {
        same(nth(gr, 1), nth(st, 1))
            && leading_number(nth(gr, 2)) <= pos
            && leading_number(nth(gr, 3)) >= pos
    });
    if let Some(gr) = over {
        let (gt, filter) = (nth(gr, 6), nth(gr, 9));
        return match Gt::of(gt, haploid) {
            Gt::Ref => o.reference(nth(gr, 7), nth(gr, 8), filter, true),
            state => (
                b"no call".to_vec(),
                unread(b"the gVCF record over the site", gt, filter, state),
            ),
        };
    }
    if let Some(w) = why {
        return (b"no call".to_vec(), w);
    }
    let why: &[u8] = if o.gvcf {
        b"no record in the VCF or the gVCF covers the site"
    } else {
        b"no record in the VCF, and no gVCF to tell a reference call from no call"
    };
    (b"no call".to_vec(), why.to_vec())
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

    /// The call and detail at `site` (CHROM POS REF ALT), joined by `|`.
    fn run_at(site: &str, vcf: &str, gvcf: &str, gvcf_on: bool) -> String {
        let sites = format!(
            "# test\nchrom\tpos\tref\talt\tgene\thgvs\tclinvar_id\n{site}\tG1\tc.1A>G\t1\n"
        );
        let o = Options {
            min_gq: "30".into(),
            min_gq_value: 30,
            min_dp: "10".into(),
            min_dp_value: 10,
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

    fn run(vcf: &str, gvcf: &str, gvcf_on: bool) -> String {
        run_at("chr1\t100\tA\tG", vcf, gvcf, gvcf_on)
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
            run("", "", true),
            "no call|no record in the VCF or the gVCF covers the site"
        );
        assert_eq!(
            run("", "", false),
            "no call|no record in the VCF, and no gVCF to tell a reference call from no call"
        );
    }

    /// Only PASS, `.` (or no FILTER column) and RefCall read as plain
    /// reference; any other FILTER on a `0/0` is a low-confidence reference
    /// call that names it.
    #[test]
    fn a_filtered_reference_call_is_low_confidence() {
        for f in ["PASS", ".", "RefCall"] {
            assert_eq!(
                run(
                    "",
                    &format!("chr1\t90\t120\tA\t<*>\t0/0\t99\t.\t{f}\n"),
                    true
                ),
                "ref|GQ 99 (gVCF block)",
                "FILTER {f}"
            );
        }
        assert_eq!(
            run("", "chr1\t90\t120\tA\t<*>\t0/0\t99\t.\tLowQual\n", true),
            "ref (low confidence)|GQ 99 (gVCF block), FILTER LowQual"
        );
        assert_eq!(
            run("chr1\t100\t100\tA\tG\t0/0\t99\t30\tLowGQ;LowDP\n", "", true),
            "ref (low confidence)|GQ 99, DP 30, FILTER LowGQ;LowDP"
        );
    }

    /// A reference call on fewer reads than the step's depth floor is
    /// shallow, from the record's DP or the block's.
    #[test]
    fn a_shallow_reference_call_says_so() {
        assert_eq!(
            run("chr1\t100\t100\tA\tG\t0/0\t99\t1\n", "", true),
            "ref (shallow)|GQ 99, DP 1, below the depth floor of DP 10"
        );
        assert_eq!(
            run("", "chr1\t90\t120\tA\t<*>\t0/0\t99\t1\tPASS\n", true),
            "ref (shallow)|GQ 99, DP 1 (gVCF block), below the depth floor of DP 10"
        );
        assert_eq!(
            run("chr1\t100\t100\tA\tG\t0/0\t99\t10\n", "", true),
            "ref|GQ 99, DP 10"
        );
        // Below both floors: the GQ makes it low confidence, and both are named.
        assert_eq!(
            run("chr1\t100\t100\tA\tG\t0/0\t3\t1\n", "", true),
            "ref (low confidence)|GQ 3, DP 1, below the call floor of GQ 30, below the depth floor of DP 10"
        );
        // A carrier's depth is shown, not floored.
        assert_eq!(
            run("chr1\t100\t100\tA\tG\t0/1\t45\t4\n", "", true),
            "het|GQ 45, DP 4"
        );
    }

    /// A reference call with no GQ was never checked against the floor.
    #[test]
    fn a_reference_call_without_gq_is_low_confidence() {
        assert_eq!(
            run("chr1\t100\t100\tA\tG\t0/0\t.\t30\n", "", true),
            "ref (low confidence)|DP 30, no GQ to check against the call floor of GQ 30"
        );
        assert_eq!(
            run("", "chr1\t90\t120\tA\t<*>\t0/0\t.\t.\n", true),
            "ref (low confidence)|(gVCF block), no GQ to check against the call floor of GQ 30"
        );
    }

    /// GLA c.640-801G>A, on chrX outside the pseudoautosomal regions.
    const GLA: &str = "chrX\t101399747\tC\tT";

    /// One copy of chrX or chrY outside the pseudoautosomal regions: an
    /// alt `1` is hemizygous.
    #[test]
    fn a_haploid_alt_call_is_hemizygous() {
        assert_eq!(
            run_at(
                GLA,
                "chrX\t101399747\t101399747\tC\tT\t1\t99\t30\n",
                "",
                true
            ),
            "hemizygous|GQ 99, DP 30"
        );
        assert_eq!(
            run_at(
                GLA,
                "chrX\t101399747\t101399747\tC\tG\t1\t99\t30\n",
                "",
                true
            ),
            "other allele|C>G hemizygous"
        );
        assert_eq!(
            run_at(
                "chrY\t3000000\tC\tT",
                "chrY\t3000000\t3000000\tC\tT\t1\t99\t30\n",
                "",
                true
            ),
            "hemizygous|GQ 99, DP 30"
        );
    }

    /// There a `0` is a reference call, under the same floors. Elsewhere
    /// (inside PAR1, on an autosome) a one-copy GT is no call, with why.
    #[test]
    fn a_haploid_reference_call_is_ref() {
        assert_eq!(
            run_at(
                GLA,
                "chrX\t101399747\t101399747\tC\tT\t0\t99\t30\n",
                "",
                true
            ),
            "ref|GQ 99, DP 30"
        );
        assert_eq!(
            run_at(
                GLA,
                "chrX\t101399747\t101399747\tC\tT\t0\t20\t30\n",
                "",
                true
            ),
            "ref (low confidence)|GQ 20, DP 30, below the call floor of GQ 30"
        );
        assert_eq!(
            run_at(
                GLA,
                "",
                "chrX\t101399700\t101399800\tC\t<NON_REF>\t0\t50\t.\n",
                true
            ),
            "ref|GQ 50 (gVCF block)"
        );
        // Inside PAR1, and on an autosome, the genome has two copies.
        assert_eq!(
            run_at(
                "chrX\t100000\tC\tT",
                "chrX\t100000\t100000\tC\tT\t0\t99\t30\n",
                "",
                true
            ),
            "no call|the VCF record at the site has one copy (GT 0) where the genome has two"
        );
        assert_eq!(
            run("", "chr1\t90\t120\tA\t<*>\t1\t50\t.\n", true),
            "no call|the gVCF record over the site has one copy (GT 1) where the genome has two"
        );
    }

    /// A record without a genotype says so, and never that no record exists.
    #[test]
    fn no_call_says_what_was_there() {
        assert_eq!(
            run("", "chr1\t90\t120\tA\t<*>\t./.\t.\t.\n", true),
            "no call|the gVCF record over the site has no genotype (GT ./.)"
        );
        assert_eq!(
            run("chr1\t100\t100\tA\tG\t./.\t.\t.\n", "", false),
            "no call|the VCF record at the site has no genotype (GT ./.)"
        );
        assert_eq!(
            run("chr1\t100\t100\tA\tG\t./.\t.\t.\n", "", true),
            "no call|the VCF record at the site has no genotype (GT ./.)"
        );
        assert_eq!(
            run("", "chr1\t100\t100\tA\tG,<*>\t0/1\t30\t20\tLowQual\n", true),
            "no call|the gVCF record over the site is not a reference call (GT 0/1, FILTER LowQual)"
        );
        assert_eq!(
            run("", "chr1\t90\t120\tA\t<*>\t0/.\t30\t.\n", true),
            "no call|the gVCF record over the site is not a reference call (GT 0/.)"
        );
    }
}
