//! GA4GH VRS 2.0 computed identifiers (`ga4gh:VA.…`) for the alleles a fact
//! cites, worked out here from the GRCh38 reference with no network.
//!
//! The steps are those of the reference implementation, vrs-python's
//! `normalize.py` with bioutils' trim and roll, and the VRS 2.0 spec
//! (<https://vrs.ga4gh.org/en/2.0/conventions/normalization.html>):
//!
//! 1. Trim the bases the reference and the alternate share, first on the
//!    left, then on the right.
//! 2. Both still hold bases: a substitution, stated as the alternate bases
//!    (`LiteralSequenceExpression`).
//! 3. One is empty (an insertion or a deletion): roll it left and right as
//!    far as the reference allows and state the allele over that whole
//!    stretch. A deletion, or an insertion that copies the reference, is a
//!    `ReferenceLengthExpression`; any other insertion is literal.
//!
//! The identifier is the sha512t24u digest of the allele's canonical JSON,
//! whose location is replaced by its own digest. The sequence part is the
//! chromosome's refget identifier from `assets/ownhealth/grch38_refget.tsv`:
//! the no-ALT analysis set the pipeline runs masks parts of chr5, chr19,
//! chr22 and chrY, so the digest of its own sequence is not GRCh38's.

use sha2::{Digest, Sha512};

/// Reads reference bases, 0-based and half-open, upper case.
pub trait Bases {
    /// The length of the chromosome.
    fn chrom_len(&self) -> u64;
    /// The bases `[beg, end)`.
    fn get(&mut self, beg: u64, end: u64) -> Result<Vec<u8>, String>;
}

/// The state of a normalised allele.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum State {
    /// The alternate bases as written.
    Literal(String),
    /// A stretch of `length` bases copied from the reference with a repeat
    /// unit of `repeat` bases.
    RefLength { length: u64, repeat: u64 },
}

/// An allele after normalisation: interbase start and end, and its state.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Allele {
    pub start: u64,
    pub end: u64,
    pub state: State,
}

/// sha512 truncated to 24 bytes, base64url without padding.
pub fn sha512t24u(data: &[u8]) -> String {
    base64url(&Sha512::digest(data)[..24])
}

/// Base64 with the URL-safe alphabet; 24 bytes need no padding.
fn base64url(bytes: &[u8]) -> String {
    const ABC: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-_";
    let mut out = String::with_capacity(bytes.len() * 4 / 3 + 3);
    for chunk in bytes.chunks(3) {
        let n = chunk.iter().fold(0u32, |n, b| n << 8 | u32::from(*b)) << (8 * (3 - chunk.len()));
        for i in 0..=chunk.len() {
            out.push(char::from(ABC[(n >> (18 - 6 * i) & 63) as usize]));
        }
    }
    out
}

fn trim_left(a: &[u8], b: &[u8]) -> usize {
    a.iter().zip(b).take_while(|(x, y)| x == y).count()
}

fn trim_right(a: &[u8], b: &[u8]) -> usize {
    a.iter()
        .rev()
        .zip(b.iter().rev())
        .take_while(|(x, y)| x == y)
        .count()
}

/// An `N` where an indel rolls: the no-ALT analysis set masks stretches
/// that GRCh38 itself holds, so the roll would stop where GRCh38's would
/// not, and the location would differ from the one the identifier names.
fn masked(base: &[u8], pos: u64) -> Result<(), String> {
    if base.first() == Some(&b'N') {
        return Err(format!(
            "the reference has N at {pos}, where the indel rolls: masked or a gap, so GRCh38's bases are not known"
        ));
    }
    Ok(())
}

/// How far `allele` rolls left from `ref_pos` (the base before it), at
/// most to 0; `ref_pos` is `start - 1`, so `None` when the allele starts
/// the chromosome.
fn roll_left(seq: &mut dyn Bases, allele: &[u8], start: u64) -> Result<u64, String> {
    let n = allele.len() as u64;
    let mut d = 0u64;
    while d < start {
        let base = seq.get(start - 1 - d, start - d)?;
        masked(&base, start - d)?;
        // Python's a[-(d + 1) % n]
        let i = (n - 1 - d % n) as usize;
        if base.first() != Some(&allele[i]) {
            break;
        }
        d += 1;
    }
    Ok(d)
}

/// How far `allele` rolls right from `end`, at most to the last base.
fn roll_right(seq: &mut dyn Bases, allele: &[u8], end: u64) -> Result<u64, String> {
    let n = allele.len() as u64;
    let len = seq.chrom_len();
    let mut d = 0u64;
    while end + d < len {
        let base = seq.get(end + d, end + d + 1)?;
        masked(&base, end + d + 1)?;
        if base.first() != Some(&allele[(d % n) as usize]) {
            break;
        }
        d += 1;
    }
    Ok(d)
}

/// The factors of `n`, largest first.
fn factors(n: u64) -> Vec<u64> {
    let mut high = Vec::new();
    let mut low = Vec::new();
    let mut i = 1;
    while i * i <= n {
        if n.is_multiple_of(i) {
            high.push(n / i);
            if n / i != i {
                low.push(i);
            }
        }
        i += 1;
    }
    high.extend(low.into_iter().rev());
    high
}

/// Whether the bases of `target` past `template` repeat `template[start..]`.
fn valid_cycle(start: usize, template: &[u8], target: &[u8]) -> bool {
    let unit = &template[start..];
    target[template.len()..]
        .iter()
        .zip(unit.iter().cycle())
        .all(|(a, b)| a == b)
}

/// Normalises the VCF allele `pos` (1-based) `ref_`>`alt`. The reference
/// bases at `pos` must be `ref_`, or the reference is not the one the
/// call was made on.
pub fn normalize(seq: &mut dyn Bases, pos: u64, ref_: &str, alt: &str) -> Result<Allele, String> {
    let r = ref_.to_ascii_uppercase().into_bytes();
    let a = alt.to_ascii_uppercase().into_bytes();
    if pos == 0 || a.is_empty() || !a.iter().all(|b| b"ACGTN".contains(b)) {
        return Err(format!("{pos} {ref_}>{alt}: not a sequence allele"));
    }
    let start = pos - 1;
    let end = start + r.len() as u64;
    if end > seq.chrom_len() {
        return Err(format!(
            "{pos} {ref_}>{alt}: past the end of the chromosome"
        ));
    }
    let on_ref = seq.get(start, end)?;
    if on_ref != r {
        return Err(format!(
            "{pos} {ref_}>{alt}: the reference has {} there",
            String::from_utf8_lossy(&on_ref)
        ));
    }
    if r == a {
        return Err(format!(
            "{pos} {ref_}>{alt}: the alternate is the reference"
        ));
    }
    let left = trim_left(&r, &a);
    let (r1, a1) = (&r[left..], &a[left..]);
    let right = trim_right(r1, a1);
    let (tr, ta) = (&r1[..r1.len() - right], &a1[..a1.len() - right]);
    let (ts, te) = (start + left as u64, end - right as u64);
    if !tr.is_empty() && !ta.is_empty() {
        return Ok(Allele {
            start: ts,
            end: te,
            state: State::Literal(String::from_utf8_lossy(ta).into_owned()),
        });
    }
    // One side is empty: roll the other through the reference both ways.
    let seed = if tr.is_empty() { ta } else { tr };
    let seed_len = seed.len() as u64;
    let ld = roll_left(seq, seed, ts)?;
    let rd = roll_right(seq, seed, te)?;
    let (ns, ne) = (ts - ld, te + rd);
    let lseq = seq.get(ns, ts)?;
    let rseq = seq.get(te, ne)?;
    let ext_ref = seq.get(ns, ne)?;
    let mut ext_alt = lseq;
    ext_alt.extend_from_slice(ta);
    ext_alt.extend_from_slice(&rseq);
    let (len_ref, len_alt) = (ext_ref.len() as u64, ext_alt.len() as u64);
    let literal = |bases: &[u8]| State::Literal(String::from_utf8_lossy(bases).into_owned());
    let state = if ext_ref.is_empty() {
        literal(&ext_alt)
    } else if len_alt < len_ref {
        State::RefLength {
            length: len_alt,
            repeat: seed_len,
        }
    } else {
        factors(seed_len)
            .into_iter()
            .filter(|c| *c <= len_ref)
            .find(|c| valid_cycle((len_ref - c) as usize, &ext_ref, &ext_alt))
            .map_or_else(
                || literal(&ext_alt),
                |c| State::RefLength {
                    length: len_alt,
                    repeat: c,
                },
            )
    };
    Ok(Allele {
        start: ns,
        end: ne,
        state,
    })
}

/// The `ga4gh:VA.` identifier of a normalised allele on the sequence with
/// refget identifier `refget` (`SQ.…`).
pub fn allele_id(refget: &str, allele: &Allele) -> String {
    let location = format!(
        r#"{{"end":{},"sequenceReference":{{"refgetAccession":"{refget}","type":"SequenceReference"}},"start":{},"type":"SequenceLocation"}}"#,
        allele.end, allele.start
    );
    let state = match &allele.state {
        State::Literal(s) => format!(r#"{{"sequence":"{s}","type":"LiteralSequenceExpression"}}"#),
        State::RefLength { length, repeat } => format!(
            r#"{{"length":{length},"repeatSubunitLength":{repeat},"type":"ReferenceLengthExpression"}}"#
        ),
    };
    let body = format!(
        r#"{{"location":"{}","state":{state},"type":"Allele"}}"#,
        sha512t24u(location.as_bytes())
    );
    format!("ga4gh:VA.{}", sha512t24u(body.as_bytes()))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A chromosome held in memory.
    pub struct Mem(pub Vec<u8>);

    impl Bases for Mem {
        fn chrom_len(&self) -> u64 {
            self.0.len() as u64
        }
        fn get(&mut self, beg: u64, end: u64) -> Result<Vec<u8>, String> {
            Ok(self.0[beg as usize..end as usize].to_vec())
        }
    }

    /// GRCh38 chr19's refget identifier, public (vrs-python's docstring).
    const CHR19: &str = "SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul";

    #[test]
    fn digest_matches_the_vrs_example() {
        // The location of ga4gh.core's ga4gh_digest docstring (vrs-python 2.3.3).
        let loc = format!(
            r#"{{"end":44908822,"sequenceReference":{{"refgetAccession":"{CHR19}","type":"SequenceReference"}},"start":44908821,"type":"SequenceLocation"}}"#
        );
        assert_eq!(
            sha512t24u(loc.as_bytes()),
            "4t6JnYWqHwYw9WzBT_lmWBb3tLQNalkT"
        );
    }

    /// Vectors from vrs-python 2.3.3 (`normalize`, then `ga4gh_identify`) on
    /// the sequence of `SEQ` under the accession of `ACC`: pos, ref, alt,
    /// start, end, the state and the identifier.
    const SEQ: &[u8] = b"GATTACACACAGTTTACGAAAAAATTTGGCC";
    const ACC: &str = "SQ.0123456789abcdefghijklmnopqrstuv";
    #[allow(clippy::type_complexity)]
    const VECTORS: &[(u64, &str, &str, u64, u64, Option<&str>, (u64, u64), &str)] = &[
        (
            2,
            "A",
            "T",
            1,
            2,
            Some("T"),
            (0, 0),
            "uPhiAPBlXs-7Rd9ejI5aoQSedPtNSVQp",
        ),
        (
            5,
            "ACA",
            "A",
            4,
            11,
            None,
            (5, 2),
            "J8DfFu8HeQmDpzwn393XgBthEOEcSKzP",
        ),
        (
            5,
            "A",
            "ACA",
            4,
            11,
            None,
            (9, 2),
            "7xgSkv1EZeoqqyb8OGOCTOCrnjI-5WKs",
        ),
        (
            4,
            "T",
            "TAC",
            4,
            11,
            None,
            (9, 2),
            "7xgSkv1EZeoqqyb8OGOCTOCrnjI-5WKs",
        ),
        (
            12,
            "G",
            "GAGTCC",
            12,
            12,
            Some("AGTCC"),
            (0, 0),
            "mFM2OkyIyB1Y0SYpobncoO0hyMXcNPiQ",
        ),
        (
            18,
            "G",
            "GA",
            18,
            24,
            None,
            (7, 1),
            "PtkXhzcGI2GPzPY3z0K2pp_YPlLnlZbz",
        ),
        (
            18,
            "GA",
            "G",
            18,
            24,
            None,
            (5, 1),
            "2VAPcqYW5GXCNVDj7yzFUNZqaycnRomq",
        ),
        (
            1,
            "G",
            "C",
            0,
            1,
            Some("C"),
            (0, 0),
            "KNkKfJNR6bKKEE6RHv5re1dEfzcp8zPS",
        ),
        (
            3,
            "TT",
            "AG",
            2,
            4,
            Some("AG"),
            (0, 0),
            "OxIKxl7A8Vm9kQ_jiHmWlNFQMhEuhpFh",
        ),
        (
            26,
            "TT",
            "T",
            24,
            27,
            None,
            (2, 1),
            "8UJCLYd2IruR79t-fLGp69ZL02hMJj7N",
        ),
        (
            17,
            "C",
            "CTTTAC",
            12,
            17,
            None,
            (10, 5),
            "ATK4xIW0ryIU0ldhrDBG_1AkPZS3hre3",
        ),
        (
            12,
            "GTT",
            "GAT",
            12,
            13,
            Some("A"),
            (0, 0),
            "U4eI5medOMgm7cyRu0WkSLg5rQz5TGRb",
        ),
        (
            1,
            "GA",
            "G",
            1,
            2,
            None,
            (0, 1),
            "4ejeBEJ7raIOAgTtMcHG93I4ud-WFoMD",
        ),
        (
            30,
            "CC",
            "C",
            29,
            31,
            None,
            (1, 1),
            "TPETeqCf4ggh3KNr375VvmApoMvLTa16",
        ),
        (
            30,
            "C",
            "CC",
            29,
            31,
            None,
            (3, 1),
            "_JVQeRO5trNUz-yD_0zet71N9iCBEf94",
        ),
        (
            9,
            "ACAG",
            "A",
            9,
            12,
            None,
            (0, 3),
            "7_ZYTDiuP3GezrNBOUHEsDppBfa-v1ZK",
        ),
        (
            14,
            "TTA",
            "TTATTA",
            13,
            16,
            None,
            (6, 3),
            "R32FBJZmzxAyQma-YGDaXV296om6q5y-",
        ),
        (
            13,
            "T",
            "TTT",
            12,
            15,
            None,
            (5, 2),
            "-M0_cCdsTTodvoJH9-tOzcZRSZq7W003",
        ),
        (
            2,
            "ATTA",
            "A",
            1,
            5,
            None,
            (1, 3),
            "ejIMVz9TYvdN0qGwSAVApbodEzceYW0y",
        ),
        (
            2,
            "A",
            "AAA",
            1,
            2,
            None,
            (3, 1),
            "A7I0xVFLHMT8nM-_boTHtDhHY8GUuMqz",
        ),
        (
            2,
            "A",
            "AAAA",
            1,
            2,
            None,
            (4, 1),
            "GTy8T3v0RWLgrXQMFOWWJR_3-oGEC_wr",
        ),
        (
            2,
            "A",
            "ACA",
            1,
            2,
            Some("ACA"),
            (0, 0),
            "cnEzR-_IIhsMxu-3p20VUqhimd4-_h90",
        ),
    ];

    #[test]
    fn alleles_match_vrs_python() {
        for &(pos, r, a, start, end, lit, (length, repeat), id) in VECTORS {
            let got = normalize(&mut Mem(SEQ.to_vec()), pos, r, a).expect("normalises");
            let state = match lit {
                Some(s) => State::Literal(s.into()),
                None => State::RefLength { length, repeat },
            };
            assert_eq!(got, Allele { start, end, state }, "{pos} {r}>{a}");
            assert_eq!(
                allele_id(ACC, &got),
                format!("ga4gh:VA.{id}"),
                "{pos} {r}>{a}"
            );
        }
    }

    #[test]
    fn refuses_an_allele_the_reference_does_not_hold() {
        let err = normalize(&mut Mem(SEQ.to_vec()), 19, "G", "GA").expect_err("A at 19");
        assert!(err.contains("the reference has A"), "{err}");
        assert!(normalize(&mut Mem(SEQ.to_vec()), 2, "A", "A").is_err());
        assert!(normalize(&mut Mem(SEQ.to_vec()), 2, "A", "<DEL>").is_err());
        assert!(normalize(&mut Mem(SEQ.to_vec()), 31, "CC", "C").is_err());
    }

    #[test]
    fn refuses_to_roll_an_indel_into_a_masked_base() {
        // The A run of a one-base deletion reaches an N: GRCh38's own
        // bases there are not known, so no location is stated.
        let e = normalize(&mut Mem(b"GAAAANAAAG".to_vec()), 1, "GA", "G").unwrap_err();
        assert!(e.contains("N at 6"), "{e}");
        let ok = normalize(&mut Mem(b"GAAAACAAAG".to_vec()), 1, "GA", "G").unwrap();
        assert_eq!((ok.start, ok.end), (1, 5));
    }

    #[test]
    fn factors_largest_first() {
        assert_eq!(factors(12), vec![12, 6, 4, 3, 2, 1]);
        assert_eq!(factors(1), vec![1]);
        assert_eq!(factors(9), vec![9, 3, 1]);
    }
}
