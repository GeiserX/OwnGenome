//! Version order. Never `sort -V` over raw tag lists: it ranks conda build
//! strings and suffixed tags (3.14-slim, 1.10.0-gpu, 3.15.0rc1) as
//! versions. A tag is compared only with tags of the same shape as the pin.

use crate::{Lookup, failed};
use regex::Regex;
use std::cmp::Ordering;
use std::sync::LazyLock;

/// A run of digits as a number of any size: ordered by value.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Num(String);

impl Num {
    fn new(digits: &str) -> Num {
        let d = digits.trim_start_matches('0');
        Num(if d.is_empty() {
            "0".into()
        } else {
            d.to_owned()
        })
    }
}

impl Ord for Num {
    fn cmp(&self, o: &Self) -> Ordering {
        self.0.len().cmp(&o.0.len()).then_with(|| self.0.cmp(&o.0))
    }
}

impl PartialOrd for Num {
    fn partial_cmp(&self, o: &Self) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}

/// Every run of digits in `text`, as numbers: `release_116.2` is (116, 2).
/// Compared as Python compares tuples: element by element, a prefix first.
pub fn ver(text: &str) -> Vec<Num> {
    static D: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"\d+").expect("digits"));
    D.find_iter(text).map(|m| Num::new(m.as_str())).collect()
}

/// The first element whose key is greatest (Python's `max(xs, key=...)`).
pub fn max_by<T, K: Ord>(xs: impl IntoIterator<Item = T>, key: impl Fn(&T) -> K) -> Option<T> {
    let mut best: Option<(K, T)> = None;
    for x in xs {
        let k = key(&x);
        if best.as_ref().is_none_or(|(bk, _)| k > *bk) {
            best = Some((k, x));
        }
    }
    best.map(|(_, x)| x)
}

static BIOCONDA: LazyLock<Regex> = LazyLock::new(|| {
    Regex::new(r"^(?P<v>\d+(?:\.\d+)*)--(?P<hash>[A-Za-z0-9]+)_(?P<build>\d+)$").expect("bioconda")
});

/// A tag's version, when it is a biocontainers tag (`1.3.6--h43da1c4_0`).
pub fn bioconda_version(tag: &str) -> Option<&str> {
    BIOCONDA
        .captures(tag)
        .map(|c| c.name("v").expect("v").as_str())
}

/// The sort key of a tag: its version numbers and its build.
pub type Key = (Vec<Num>, Num);

/// `tag`'s sort key when it has the same shape as the pinned tag `like`;
/// `None` when it does not.
pub fn tag_key(tag: &str, like: &str) -> Option<Key> {
    if BIOCONDA.is_match(like) {
        let c = BIOCONDA.captures(tag)?;
        return Some((ver(&c["v"]), Num::new(&c["build"])));
    }
    static D: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"\d+").expect("digits"));
    let shape = |s: &str| D.replace_all(s, "N").into_owned();
    if !like.bytes().any(|b| b.is_ascii_digit()) || shape(tag) != shape(like) {
        return None;
    }
    Some((ver(tag), Num::new("0")))
}

/// The newest tag overall and the newest on the pin's major, each newer
/// than the pin, or `None`.
pub fn newer_tags(pin: &str, tags: &[String]) -> Lookup<(Option<String>, Option<String>)> {
    let Some(pk) = tag_key(pin, pin) else {
        return failed(format!("cannot read a version from the pinned tag {pin}"));
    };
    let keyed: Vec<(Key, &String)> = tags
        .iter()
        .filter_map(|t| tag_key(t, pin).filter(|k| *k > pk).map(|k| (k, t)))
        .collect();
    let newest = keyed.iter().max().map(|(_, t)| (*t).clone());
    let same = keyed
        .iter()
        .filter(|(k, _)| k.0.first() == pk.0.first())
        .max()
        .map(|(_, t)| (*t).clone());
    Ok((newest, same))
}

/// How far `new` is from the pin: `major`, `minor`, `patch` or `build`.
pub fn change_kind(pin: &str, new: &str) -> &'static str {
    let a = tag_key(pin, pin).map(|k| k.0).unwrap_or_default();
    let b = tag_key(new, pin).map(|k| k.0).unwrap_or_default();
    let head = |v: &[Num], n: usize| v[..n.min(v.len())].to_vec();
    if head(&a, 1) != head(&b, 1) {
        "major"
    } else if head(&a, 2) != head(&b, 2) {
        "minor"
    } else if a != b {
        "patch"
    } else {
        "build"
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn v(xs: &[&str]) -> Vec<String> {
        xs.iter().map(|x| (*x).to_owned()).collect()
    }

    fn pair(a: Option<&str>, b: Option<&str>) -> (Option<String>, Option<String>) {
        (a.map(str::to_owned), b.map(str::to_owned))
    }

    #[test]
    fn the_sort_v_traps() {
        let bio = v(&[
            "1.3.6--h43da1c4_0",
            "1.3.10--h5ca1c30_0",
            "1.3.9--hdeadbe_2",
            "1.3.6--habc123_1",
            "latest",
        ]);
        assert_eq!(
            newer_tags("1.3.6--h43da1c4_0", &bio).unwrap(),
            pair(Some("1.3.10--h5ca1c30_0"), Some("1.3.10--h5ca1c30_0"))
        );
        assert_eq!(
            change_kind("1.3.6--h43da1c4_0", "1.3.6--habc123_1"),
            "build"
        );
        let py = v(&[
            "3.11",
            "3.9",
            "3.14",
            "3.14-slim",
            "3.15.0rc1",
            "3.11.9",
            "3.15-rc",
            "windowsservercore",
            "3.14.0",
        ]);
        assert_eq!(
            newer_tags("3.11", &py).unwrap(),
            pair(Some("3.14"), Some("3.14"))
        );
        let vep = v(&[
            "release_116.0",
            "release_116.2",
            "release_117.0",
            "latest",
            "116.2",
        ]);
        assert_eq!(
            newer_tags("release_116.0", &vep).unwrap(),
            pair(Some("release_117.0"), Some("release_116.2"))
        );
        assert_eq!(change_kind("release_116.0", "release_117.0"), "major");
        assert_eq!(
            newer_tags("1.21", &v(&["1.21", "1.20", "1.21-2"])).unwrap(),
            pair(None, None)
        );
        assert!(newer_tags("latest", &v(&["latest"])).is_err());
    }

    #[test]
    fn numbers_of_any_size_order_by_value() {
        assert!(ver("1.10") > ver("1.9"));
        assert!(ver("1.0.0") > ver("1.0"));
        assert!(ver("99999999999999999999999999999999999999999") > ver("9"));
        assert_eq!(ver("007"), ver("7"));
        assert_eq!(max_by(["1.0", "1.00", "0.9"], |s| ver(s)), Some("1.0"));
    }
}
