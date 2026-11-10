//! File times and file lists as Python sees them: `os.path.getmtime` as a
//! float, `os.stat().st_mtime_ns`, `glob.glob` in its own order,
//! `max(paths, key=os.path.getmtime)` with its tie rule, and `sorted()` of
//! paths. The summary picks a step's file with these (`first_existing`,
//! `newest`), so the same folder must give the same file.
//!
//! Paths are bytes, as on Linux. Where Python works on `str`, a byte that
//! is not UTF-8 stands for U+DC80 to U+DCFF (`surrogateescape`), as
//! `os.fsdecode` makes it; that decides `?`, `[a-z]` and sort order.

use std::cmp::Ordering;
use std::ffi::{OsStr, OsString};
use std::fs;
use std::io;
use std::os::unix::ffi::{OsStrExt, OsStringExt};
use std::os::unix::fs::MetadataExt;
use std::path::{Path, PathBuf};

/// The float CPython makes of a stat time: `sec + nsec * 1e-9`, rounded
/// twice, so times a few hundred nanoseconds apart can be equal.
pub fn from_parts(sec: i64, nsec: i64) -> f64 {
    sec as f64 + nsec as f64 * 1e-9
}

/// `os.path.getmtime(path)`, following symlinks.
pub fn mtime(path: impl AsRef<Path>) -> io::Result<f64> {
    let m = fs::metadata(path)?;
    Ok(from_parts(m.mtime(), m.mtime_nsec()))
}

/// `os.stat(path).st_mtime_ns`, exact.
pub fn mtime_ns(path: impl AsRef<Path>) -> io::Result<i128> {
    let m = fs::metadata(path)?;
    Ok(i128::from(m.mtime()) * 1_000_000_000 + i128::from(m.mtime_nsec()))
}

/// `max(paths, key=os.path.getmtime)`, or `None` for no paths: the first
/// of the paths whose float mtime is the largest. An error reading one
/// stops it, as the exception does in Python.
pub fn newest<P: AsRef<Path>>(paths: impl IntoIterator<Item = P>) -> io::Result<Option<P>> {
    let mut best: Option<(f64, P)> = None;
    for p in paths {
        let t = mtime(&p)?;
        if best.as_ref().is_none_or(|(b, _)| t > *b) {
            best = Some((t, p));
        }
    }
    Ok(best.map(|(_, p)| p))
}

/// Compares two paths as Python compares the `str` of each: by code point,
/// a stray byte as its surrogate. For UTF-8 this is byte order.
pub fn py_cmp(a: impl AsRef<OsStr>, b: impl AsRef<OsStr>) -> Ordering {
    code_points(a.as_ref().as_bytes()).cmp(code_points(b.as_ref().as_bytes()))
}

/// `paths.sort()` as Python sorts the `str` of each path.
pub fn sort<P: AsRef<OsStr>>(paths: &mut [P]) {
    paths.sort_by(|a, b| py_cmp(a, b));
}

/// `glob.escape(path)`: `*`, `?` and `[` each wrapped in brackets.
pub fn escape(path: &str) -> String {
    let mut out = String::with_capacity(path.len());
    for c in path.chars() {
        if matches!(c, '*' | '?' | '[') {
            out.push('[');
            out.push(c);
            out.push(']');
        } else {
            out.push(c);
        }
    }
    out
}

/// `fnmatch.fnmatchcase(name, pattern)`: `*`, `?` and `[...]` (with `!`
/// for not, and Python's rules for a `-`, a `]` and a range written
/// backwards), case-sensitive, over the whole name.
pub fn fnmatch(name: impl AsRef<OsStr>, pattern: impl AsRef<OsStr>) -> bool {
    let name: Vec<u32> = code_points(name.as_ref().as_bytes()).collect();
    Pattern::new(pattern.as_ref().as_bytes()).matches(&name)
}

/// `glob.glob(pattern, recursive=recursive)`, in the order Python returns:
/// each folder as the system lists it (no sorting), `**` depth first.
/// Names starting with `.` match only a pattern that starts with `.`, and
/// `**` never enters them. A folder that cannot be read counts as empty.
pub fn glob(pattern: impl AsRef<OsStr>, recursive: bool) -> Vec<PathBuf> {
    let p = pattern.as_ref().as_bytes();
    let mut out = Glob { recursive }.iglob(p, false);
    // glob.iglob drops the "" that "**" yields for its own folder.
    if (p.is_empty() || recursive && p.starts_with(b"**")) && out.first().is_some_and(Vec::is_empty)
    {
        out.remove(0);
    }
    out.into_iter()
        .map(|b| PathBuf::from(OsString::from_vec(b)))
        .collect()
}

/// The code points Python's `os.fsdecode` gives for these bytes.
fn code_points(b: &[u8]) -> impl Iterator<Item = u32> + '_ {
    b.utf8_chunks().flat_map(|c| {
        c.valid()
            .chars()
            .map(u32::from)
            .chain(c.invalid().iter().map(|&x| 0xDC00 + u32::from(x)))
    })
}

fn has_magic(s: &[u8]) -> bool {
    s.iter().any(|c| matches!(c, b'*' | b'?' | b'['))
}

fn is_hidden(name: &[u8]) -> bool {
    name.first() == Some(&b'.')
}

/// `os.path.split`.
fn split(p: &[u8]) -> (&[u8], &[u8]) {
    let i = p.iter().rposition(|&c| c == b'/').map_or(0, |i| i + 1);
    let (mut head, tail) = p.split_at(i);
    if head.iter().any(|&c| c != b'/') {
        while head.last() == Some(&b'/') {
            head = &head[..head.len() - 1];
        }
    }
    (head, tail)
}

/// `os.path.join(a, b)`.
fn join(a: &[u8], b: &[u8]) -> Vec<u8> {
    if b.first() == Some(&b'/') || a.is_empty() {
        b.to_vec()
    } else if a.last() == Some(&b'/') {
        [a, b].concat()
    } else {
        [a, b"/", b].concat()
    }
}

/// glob's `_join`: the other one when either is empty.
fn join_nonempty(a: &[u8], b: &[u8]) -> Vec<u8> {
    if a.is_empty() || b.is_empty() {
        [a, b].concat()
    } else {
        join(a, b)
    }
}

fn path(b: &[u8]) -> &Path {
    Path::new(OsStr::from_bytes(b))
}

fn lexists(p: &[u8]) -> bool {
    fs::symlink_metadata(path(p)).is_ok()
}

fn isdir(p: &[u8]) -> bool {
    fs::metadata(path(p)).is_ok_and(|m| m.is_dir())
}

/// glob's `_listdir`: the names in one folder, in the system's order;
/// with `dironly`, only those that are folders or symlinks to one.
fn listdir(dir: &[u8], dironly: bool) -> Vec<Vec<u8>> {
    let dir = if dir.is_empty() { b"." } else { dir };
    let Ok(entries) = fs::read_dir(path(dir)) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for e in entries {
        let Ok(e) = e else { break };
        let keep = !dironly
            || match e.file_type() {
                Ok(t) if t.is_symlink() => fs::metadata(e.path()).is_ok_and(|m| m.is_dir()),
                Ok(t) => t.is_dir(),
                Err(_) => false,
            };
        if keep {
            out.push(e.file_name().into_vec());
        }
    }
    out
}

struct Glob {
    recursive: bool,
}

impl Glob {
    fn iglob(&self, pathname: &[u8], dironly: bool) -> Vec<Vec<u8>> {
        let (dirname, basename) = split(pathname);
        if !has_magic(pathname) {
            let found = if basename.is_empty() {
                isdir(dirname)
            } else {
                lexists(pathname)
            };
            return if found {
                vec![pathname.to_vec()]
            } else {
                Vec::new()
            };
        }
        if dirname.is_empty() {
            return self.in_dir(b"", basename, dironly);
        }
        let dirs = if dirname != pathname && has_magic(dirname) {
            self.iglob(dirname, true)
        } else {
            vec![dirname.to_vec()]
        };
        let mut out = Vec::new();
        for d in dirs {
            for name in self.in_dir(&d, basename, dironly) {
                out.push(join(&d, &name));
            }
        }
        out
    }

    /// `_glob0`, `_glob1` or `_glob2`, by what the last part holds.
    fn in_dir(&self, dir: &[u8], pattern: &[u8], dironly: bool) -> Vec<Vec<u8>> {
        if self.recursive && pattern == b"**" {
            let mut out = Vec::new();
            if dir.is_empty() || isdir(dir) {
                out.push(Vec::new());
            }
            rlistdir(dir, dironly, b"", &mut out);
            out
        } else if has_magic(pattern) {
            let pat = Pattern::new(pattern);
            let hidden_ok = is_hidden(pattern);
            listdir(dir, dironly)
                .into_iter()
                .filter(|n| hidden_ok || !is_hidden(n))
                .filter(|n| {
                    let cps: Vec<u32> = code_points(n).collect();
                    pat.matches(&cps)
                })
                .collect()
        } else if pattern.is_empty() {
            if isdir(dir) {
                vec![Vec::new()]
            } else {
                Vec::new()
            }
        } else if lexists(&join_nonempty(dir, pattern)) {
            vec![pattern.to_vec()]
        } else {
            Vec::new()
        }
    }
}

/// glob's `_rlistdir`: every name under `dir`, depth first, skipping
/// hidden names, each joined to `rel`, the path from where `**` started.
fn rlistdir(dir: &[u8], dironly: bool, rel: &[u8], out: &mut Vec<Vec<u8>>) {
    for x in listdir(dir, dironly) {
        if is_hidden(&x) {
            continue;
        }
        let here = join_nonempty(rel, &x);
        out.push(here.clone());
        rlistdir(&join_nonempty(dir, &x), dironly, &here, out);
    }
}

/// One position of a compiled `fnmatch` pattern.
#[derive(Debug)]
enum Tok {
    Star,
    Lit(u32),
    Any,
    /// `[...]`: negated, then single characters and inclusive ranges.
    Set(bool, Vec<(u32, u32)>),
    /// `[]` that can match nothing (Python's `(?!)`).
    Never,
}

impl Tok {
    fn one(&self, c: u32) -> bool {
        match self {
            Tok::Lit(l) => *l == c,
            Tok::Any => true,
            Tok::Set(neg, items) => items.iter().any(|&(lo, hi)| lo <= c && c <= hi) != *neg,
            Tok::Star | Tok::Never => false,
        }
    }
}

struct Pattern(Vec<Tok>);

const BANG: u32 = '!' as u32;
const DASH: u32 = '-' as u32;
const CLOSE: u32 = ']' as u32;
const OPEN: u32 = '[' as u32;
const CARET: u32 = '^' as u32;
const BACKSLASH: u32 = '\\' as u32;

impl Pattern {
    /// `fnmatch._translate`, kept as tokens instead of a regex.
    fn new(pattern: &[u8]) -> Pattern {
        let pat: Vec<u32> = code_points(pattern).collect();
        let n = pat.len();
        let mut toks = Vec::new();
        let mut i = 0;
        while i < n {
            let c = pat[i];
            i += 1;
            match char::from_u32(c) {
                Some('*') => {
                    if !matches!(toks.last(), Some(Tok::Star)) {
                        toks.push(Tok::Star);
                    }
                }
                Some('?') => toks.push(Tok::Any),
                Some('[') => {
                    let mut j = i;
                    if j < n && pat[j] == BANG {
                        j += 1;
                    }
                    if j < n && pat[j] == CLOSE {
                        j += 1;
                    }
                    while j < n && pat[j] != CLOSE {
                        j += 1;
                    }
                    if j >= n {
                        toks.push(Tok::Lit(c));
                    } else {
                        toks.push(set(&pat, i, j));
                        i = j + 1;
                    }
                }
                _ => toks.push(Tok::Lit(c)),
            }
        }
        Pattern(toks)
    }

    /// Whole-name match. Each token but `*` takes one character, so taking
    /// the earliest place after each `*` and retrying from the last `*` on
    /// a miss finds a match whenever Python's regex does.
    fn matches(&self, s: &[u32]) -> bool {
        let t = &self.0;
        let (mut ti, mut si) = (0, 0);
        let mut star: Option<(usize, usize)> = None;
        while si < s.len() {
            match t.get(ti) {
                Some(Tok::Star) => {
                    star = Some((ti, si));
                    ti += 1;
                    continue;
                }
                Some(tok) if tok.one(s[si]) => {
                    ti += 1;
                    si += 1;
                    continue;
                }
                _ => {}
            }
            match star {
                Some((st, ss)) => {
                    ti = st + 1;
                    si = ss + 1;
                    star = Some((st, ss + 1));
                }
                None => return false,
            }
        }
        t[ti..].iter().all(|x| matches!(x, Tok::Star))
    }
}

/// The set `pat[start..end]` stands for (`end` is its `]`): the regex set
/// `fnmatch._translate` writes, then read as `re` reads it. Each character
/// carries whether the regex has it escaped.
fn set(pat: &[u32], start: usize, end: usize) -> Tok {
    let stuff = &pat[start..end];
    let mut body: Vec<(u32, bool)> = if !stuff.contains(&DASH) {
        stuff.iter().map(|&c| (c, c == BACKSLASH)).collect()
    } else {
        // Split at each '-' that can make a range: one right after the
        // start (or after a leading '!') or right after a range is literal.
        let mut chunks: Vec<Vec<u32>> = Vec::new();
        let mut i = start;
        let mut k = if pat[start] == BANG {
            start + 2
        } else {
            start + 1
        };
        while let Some(at) = (k < end)
            .then(|| pat[k..end].iter().position(|&c| c == DASH).map(|o| k + o))
            .flatten()
        {
            chunks.push(pat[i..at].to_vec());
            i = at + 1;
            k = at + 3;
        }
        if i < end {
            chunks.push(pat[i..end].to_vec());
        } else {
            chunks.last_mut().expect("a '-' was found").push(DASH);
        }
        // Drop each range written backwards, with both its ends.
        for k in (1..chunks.len()).rev() {
            if chunks[k - 1].last() > chunks[k].first() {
                let next = chunks.remove(k);
                let prev = &mut chunks[k - 1];
                prev.pop();
                prev.extend_from_slice(&next[1..]);
            }
        }
        let mut body = Vec::new();
        for (n, chunk) in chunks.iter().enumerate() {
            if n > 0 {
                body.push((DASH, false));
            }
            body.extend(chunk.iter().map(|&c| (c, c == BACKSLASH || c == DASH)));
        }
        body
    };
    match body[..] {
        [] => return Tok::Never,
        [(BANG, false)] => return Tok::Any,
        _ => {}
    }
    let neg = body[0] == (BANG, false);
    if neg {
        body.remove(0);
    } else if matches!(body[0].0, CARET | OPEN) {
        body[0].1 = true;
    }
    // re's set parser: a '-' between two characters makes a range; any
    // other character, escaped or not, stands for itself.
    let mut items = Vec::new();
    let mut p = 0;
    while p < body.len() {
        let lo = body[p].0;
        p += 1;
        if p + 1 < body.len() && body[p] == (DASH, false) {
            let hi = body[p + 1].0;
            p += 2;
            if hi < lo {
                // re.error in Python; the merge above leaves none.
                return Tok::Never;
            }
            items.push((lo, hi));
        } else {
            items.push((lo, lo));
        }
    }
    Tok::Set(neg, items)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicU32, Ordering as AtomicOrdering};

    /// A fresh folder under the system temp folder, removed on drop.
    struct Tree(PathBuf);

    impl Tree {
        fn new() -> Tree {
            static N: AtomicU32 = AtomicU32::new(0);
            let p = std::env::temp_dir().join(format!(
                "pycompat-mtime-{}-{}",
                std::process::id(),
                N.fetch_add(1, AtomicOrdering::Relaxed)
            ));
            fs::create_dir_all(&p).unwrap();
            Tree(p)
        }

        fn file(&self, rel: &str, ns: i64) -> PathBuf {
            let p = self.0.join(rel);
            fs::create_dir_all(p.parent().unwrap()).unwrap();
            let f = fs::File::create(&p).unwrap();
            let t = std::time::UNIX_EPOCH + std::time::Duration::from_nanos(ns as u64);
            f.set_modified(t).unwrap();
            p
        }

        fn rel(&self, found: &[PathBuf]) -> Vec<String> {
            let mut v: Vec<String> = found
                .iter()
                .map(|p| {
                    // Not strip_prefix, which drops a trailing '/'.
                    let root = format!("{}/", self.0.to_str().unwrap());
                    p.to_str().unwrap().strip_prefix(&root).unwrap().to_owned()
                })
                .collect();
            v.sort();
            v
        }
    }

    impl Drop for Tree {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.0);
        }
    }

    #[test]
    fn float_mtime_and_ties() {
        assert_eq!(
            from_parts(1_700_000_000, 123_456_789),
            1_700_000_000.123_456_7
        );
        assert_eq!(from_parts(-2, 500_000_000), -1.5);
        // 1 ns apart, one float: max() keeps the first.
        assert_eq!(from_parts(1_700_000_000, 1), from_parts(1_700_000_000, 0));
        let t = Tree::new();
        let base = 1_700_000_000_000_000_000;
        let a = t.file("a.txt", base + 1);
        let b = t.file("b.txt", base);
        let c = t.file("c.txt", base - 1_000_000);
        assert_eq!(mtime_ns(&a).unwrap(), i128::from(base + 1));
        assert_eq!(mtime(&a).unwrap(), mtime(&b).unwrap());
        assert_eq!(newest([&c, &b, &a]).unwrap(), Some(&b));
        assert_eq!(newest([&a, &b, &c]).unwrap(), Some(&a));
        assert_eq!(newest(Vec::<PathBuf>::new()).unwrap(), None);
        let gone = t.0.join("gone");
        assert_eq!(
            newest([&a, &gone]).unwrap_err().kind(),
            io::ErrorKind::NotFound
        );
    }

    #[test]
    fn glob_matches_python() {
        let t = Tree::new();
        for f in [
            "a.txt",
            "b.vcf",
            ".h.txt",
            "d/x.txt",
            "d/e/y.txt",
            "d/.hid/z.txt",
            "[a].txt",
        ] {
            t.file(f, 0);
        }
        let root = escape(t.0.to_str().unwrap());
        let g = |p: &str, r: bool| t.rel(&glob(format!("{root}/{p}"), r));
        assert_eq!(g("*.txt", false), ["[a].txt", "a.txt"]);
        assert_eq!(g(".*", false), [".h.txt"]);
        assert_eq!(g("[[]a].txt", false), ["[a].txt"]);
        assert_eq!(
            g("**/*.txt", true),
            ["[a].txt", "a.txt", "d/e/y.txt", "d/x.txt"]
        );
        assert_eq!(g("**/*.txt", false), ["d/x.txt"]);
        assert_eq!(
            g("**", true),
            [
                "",
                "[a].txt",
                "a.txt",
                "b.vcf",
                "d",
                "d/e",
                "d/e/y.txt",
                "d/x.txt"
            ]
        );
        assert_eq!(g("*/", false), ["d/"]);
        assert_eq!(g("d", false), ["d"]);
        assert_eq!(g("nope/*", false), Vec::<String>::new());
        // A depth-first walk: each folder right before what it holds.
        let all = glob(format!("{root}/**"), true);
        let d = all.iter().position(|p| p.ends_with("d")).unwrap();
        let de = all.iter().position(|p| p.ends_with("d/e")).unwrap();
        let y = all.iter().position(|p| p.ends_with("d/e/y.txt")).unwrap();
        assert!(d < de && de < y);
    }

    #[test]
    fn fnmatch_matches_python() {
        let yes = [
            ("abc", "a*c"),
            ("abc", "*"),
            ("", "*"),
            ("a.b.c", "*.*.c"),
            ("b", "[a-c]"),
            ("-", "[a-]"),
            ("x", "[!a-c]"),
            ("]", "[]]"),
            ("[!]", "[!]"),
            ("[", "["),
            ("a\nb", "a?b"),
            ("é", "?"),
            ("\\", "[\\]"),
            ("x", "[!z-a]"),
            ("ab", "*b*"),
        ];
        for (name, pat) in yes {
            assert!(fnmatch(name, pat), "{name:?} ~ {pat:?}");
        }
        let no = [
            ("abd", "a*c"),
            ("A", "a"),
            ("b", "[!a-c]"),
            ("a", "[z-a]"),
            ("ab", "?"),
            ("x", "[]"),
        ];
        for (name, pat) in no {
            assert!(!fnmatch(name, pat), "{name:?} !~ {pat:?}");
        }
    }

    #[test]
    fn python_sort_order() {
        // Bytes put the lone 0xF5 after U+1D400; Python's str order puts
        // its surrogate U+DCF5 before it.
        let astral = OsString::from("\u{1d400}");
        let stray = OsString::from_vec(vec![0xF5]);
        let mut v = vec![astral.clone(), stray.clone(), OsString::from("a")];
        sort(&mut v);
        assert_eq!(v, [OsString::from("a"), stray, astral]);
        assert_eq!(escape("a*b?[c]"), "a[*]b[?][[]c]");
    }
}
