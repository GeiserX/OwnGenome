use smallvec::SmallVec;

/// A genotype (FORMAT/GT).
///
/// `alleles` holds each allele index, `None` for `.`. `phased` is true when
/// every separator is `|`: `0|1` and VCF 4.4's haploid `|1` are phased; `0/1`,
/// `0|1/2` and a haploid `1` are not.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Gt {
    pub alleles: SmallVec<[Option<u32>; 2]>,
    pub phased: bool,
}

impl Gt {
    /// Reads a GT value; `None` when it is not one (an empty value, a
    /// separator with no allele after it, an allele that is not a number).
    pub fn parse(s: &[u8]) -> Option<Gt> {
        let mut alleles = SmallVec::new();
        let mut seps = 0usize;
        let mut phased_seps = 0usize;
        // VCF 4.4 allows a phase sign before the first allele.
        let mut rest = s;
        if let Some((&c, after)) = rest.split_first()
            && (c == b'|' || c == b'/')
        {
            seps += 1;
            phased_seps += usize::from(c == b'|');
            rest = after;
        }
        loop {
            let end = rest
                .iter()
                .position(|&c| c == b'/' || c == b'|')
                .unwrap_or(rest.len());
            alleles.push(allele(&rest[..end])?);
            let Some(&sep) = rest.get(end) else {
                break;
            };
            seps += 1;
            phased_seps += usize::from(sep == b'|');
            rest = &rest[end + 1..];
        }
        Some(Gt {
            alleles,
            phased: seps > 0 && phased_seps == seps,
        })
    }
}

fn allele(s: &[u8]) -> Option<Option<u32>> {
    if s == b"." {
        return Some(None);
    }
    if s.is_empty() || !s.iter().all(u8::is_ascii_digit) {
        return None;
    }
    std::str::from_utf8(s).ok()?.parse().ok().map(Some)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn gt(s: &str) -> Option<(Vec<Option<u32>>, bool)> {
        Gt::parse(s.as_bytes()).map(|g| (g.alleles.to_vec(), g.phased))
    }

    #[test]
    fn genotypes() {
        assert_eq!(gt("0/1"), Some((vec![Some(0), Some(1)], false)));
        assert_eq!(gt("0|1"), Some((vec![Some(0), Some(1)], true)));
        assert_eq!(gt("1|0/2"), Some((vec![Some(1), Some(0), Some(2)], false)));
        assert_eq!(gt("./."), Some((vec![None, None], false)));
        assert_eq!(gt(".|."), Some((vec![None, None], true)));
        assert_eq!(gt("1/."), Some((vec![Some(1), None], false)));
        assert_eq!(gt("1"), Some((vec![Some(1)], false)));
        assert_eq!(gt("."), Some((vec![None], false)));
        assert_eq!(gt("|1"), Some((vec![Some(1)], true)));
        assert_eq!(gt("/0|1"), Some((vec![Some(0), Some(1)], false)));
        assert_eq!(gt("12|300"), Some((vec![Some(12), Some(300)], true)));
        for bad in [
            "",
            "/",
            "0/",
            "0//1",
            "a/1",
            "-1/0",
            "0/+1",
            "4294967296/0",
            "0 /1",
        ] {
            assert_eq!(gt(bad), None, "{bad:?}");
        }
    }
}
