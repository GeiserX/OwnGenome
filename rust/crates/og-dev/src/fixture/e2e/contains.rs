//! `contains FILE PART`: exit 0 when FILE holds the bytes of PART anywhere,
//! 1 when it does not, 2 when either cannot be read. The paralog case asks
//! it whether the pinned Paraphase binary embeds its config byte for byte.

use super::{Fail, read};
use std::io::Write;
use std::path::Path;

/// Whether `part` occurs in `hay`. An empty part occurs everywhere.
pub fn holds(hay: &[u8], part: &[u8]) -> bool {
    let Some((&first, rest)) = part.split_first() else {
        return true;
    };
    if part.len() > hay.len() {
        return false;
    }
    let last = hay.len() - part.len();
    let mut i = 0;
    while i <= last {
        match hay[i..=last].iter().position(|&b| b == first) {
            None => return false,
            Some(k) => {
                let at = i + k;
                if &hay[at + 1..at + part.len()] == rest {
                    return true;
                }
                i = at + 1;
            }
        }
    }
    false
}

pub fn run(args: &[String], err: &mut dyn Write) -> u8 {
    if args.len() != 2 {
        let _ = writeln!(err, "og dev fixture contains: takes FILE PART");
        return 2;
    }
    let load = |p: &String| {
        read(Path::new(p)).map_err(|f| match f {
            Fail::Run(e) | Fail::Usage(e) => e,
        })
    };
    match load(&args[0]).and_then(|hay| load(&args[1]).map(|part| (hay, part))) {
        Ok((hay, part)) => u8::from(!holds(&hay, &part)),
        Err(e) => {
            let _ = writeln!(err, "og dev fixture contains: {e}");
            2
        }
    }
}

#[cfg(test)]
mod tests {
    use super::holds;
    use crate::fixture::e2e::tests::{run, scratch};
    use std::fs;

    #[test]
    fn finds_a_part_anywhere_and_only_whole() {
        assert!(holds(b"xxsmn1:\n  genes: Xyy", b"smn1:\n  genes: X"));
        assert!(holds(b"abc", b"abc"));
        assert!(holds(b"abc", b""));
        assert!(holds(b"aab", b"ab"));
        assert!(!holds(b"abc", b"abd"));
        assert!(!holds(b"ab", b"abc"));
        assert!(!holds(b"a\nb", b"a b"));
    }

    #[test]
    fn exit_codes() {
        let d = scratch("contains");
        fs::write(d.join("bin"), b"\x00\x01data/38/config.yaml\n\xff").unwrap();
        fs::write(d.join("in"), b"config.yaml\n").unwrap();
        fs::write(d.join("out"), b"other.yaml\n").unwrap();
        let p = |n: &str| d.join(n).to_str().unwrap().to_owned();
        assert_eq!(run("contains", &[&p("bin"), &p("in")]).0, Some(0));
        assert_eq!(run("contains", &[&p("bin"), &p("out")]).0, Some(1));
        assert_eq!(run("contains", &[&p("bin"), &p("none")]).0, Some(2));
        assert_eq!(run("contains", &[&p("bin")]).0, Some(2));
        fs::remove_dir_all(d).unwrap();
    }
}
