//! `og pycompat-probe`: reads one case per line and writes one result per
//! line. The test at the end of this file holds the results to what CPython
//! wrote for the same cases, frozen in `tests/rust/vectors/pycompat`.
//!
//! A case is tab-separated: the operation, then its arguments. Floats travel
//! as the 16 hex digits of their bits, text as the hex of its UTF-8, so no
//! value is rounded or mangled on the way. A result is the Python result,
//! or the name of the exception Python raised, escaped by [`escape`].

use og_io::{Lines, Utf8};
use pycompat::csv::{self, DictReader, Quoting, Reader};
use pycompat::{fmt, json, mtime, num, stats, text, value};
use std::ffi::OsString;
use std::io::{self, BufRead, Write};
use std::os::unix::ffi::{OsStrExt, OsStringExt};
use std::path::PathBuf;

/// `\` `\n` `\r` `|` and NUL escaped, so every result is one line of
/// text and a list of lines can be joined with `|`.
fn escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.chars() {
        match c {
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '|' => out.push_str("\\|"),
            '\0' => out.push_str("\\0"),
            c => out.push(c),
        }
    }
    out
}

fn bad(msg: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg.into())
}

fn float(s: &str) -> io::Result<f64> {
    u64::from_str_radix(s, 16)
        .map(f64::from_bits)
        .map_err(|_| bad(format!("bad float bits {s:?}")))
}

fn bytes(s: &str) -> io::Result<Vec<u8>> {
    if !s.len().is_multiple_of(2) {
        return Err(bad("odd hex length"));
    }
    (0..s.len())
        .step_by(2)
        .map(|i| u8::from_str_radix(&s[i..i + 2], 16).map_err(|_| bad("bad hex")))
        .collect()
}

fn text(s: &str) -> io::Result<String> {
    String::from_utf8(bytes(s)?).map_err(|_| bad("hex text is not UTF-8"))
}

fn int<T: std::str::FromStr>(s: &str) -> io::Result<T> {
    s.parse().map_err(|_| bad(format!("bad integer {s:?}")))
}

fn floats(s: &str) -> io::Result<Vec<f64>> {
    if s.is_empty() {
        Ok(Vec::new())
    } else {
        s.split(',').map(float).collect()
    }
}

fn dump_opts(s: &str) -> io::Result<json::Dump<'static>> {
    let f: Vec<&str> = s.split(',').collect();
    let [indent, sep, ascii, sort, nan] = f[..] else {
        return Err(bad(format!("bad json options {s:?}")));
    };
    Ok(json::Dump {
        indent: if indent == "-" {
            None
        } else {
            Some(int(indent)?)
        },
        separators: match sep {
            "d" => None,
            "c" => Some((",", ":")),
            "p" => Some((", ", ": ")),
            "i" => Some((",", ": ")),
            _ => return Err(bad(format!("bad separators {sep:?}"))),
        },
        ensure_ascii: ascii == "1",
        sort_keys: sort == "1",
        allow_nan: nan == "1",
    })
}

/// A path or name sent as the hex of its bytes.
fn os(s: &str) -> io::Result<OsString> {
    Ok(OsString::from_vec(bytes(s)?))
}

fn hex(b: &[u8]) -> String {
    b.iter().map(|x| format!("{x:02x}")).collect()
}

fn os_list(s: &str) -> io::Result<Vec<OsString>> {
    if s.is_empty() {
        Ok(Vec::new())
    } else {
        s.split(',').map(os).collect()
    }
}

/// A path found, as hex, or Python's `None`.
fn maybe_path(p: Option<PathBuf>) -> String {
    p.map_or_else(|| "None".to_owned(), |p| hex(p.as_os_str().as_bytes()))
}

/// The OSError subclass Python raises for this error.
fn os_exception(e: &io::Error) -> &'static str {
    match e.kind() {
        io::ErrorKind::NotFound => "FileNotFoundError",
        io::ErrorKind::NotADirectory => "NotADirectoryError",
        io::ErrorKind::PermissionDenied => "PermissionError",
        _ => "OSError",
    }
}

fn py_bool(b: bool) -> String {
    if b { "True" } else { "False" }.to_owned()
}

/// The files the patterns find, each pattern's in glob order, as
/// collect_summary's first_existing and newest gather them.
fn files(patterns: &[OsString], sorted: bool) -> Vec<PathBuf> {
    let mut out = Vec::new();
    for p in patterns {
        let mut found = mtime::glob(p, true);
        if sorted {
            mtime::sort(&mut found);
        }
        out.extend(found.into_iter().filter(|f| f.is_file()));
    }
    out
}

/// The csv and text results: `%XX` for the characters that give a result
/// its structure, so no two results look alike.
fn pe(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.chars() {
        if "%[],{}=~\\|".contains(c) || u32::from(c) < 32 {
            out.push_str(&format!("%{:02X}", u32::from(c)));
        } else {
            out.push(c);
        }
    }
    out
}

/// A list as `[len:item,item]`, so `[]` and `[""]` differ.
fn render_list<S: AsRef<str>>(xs: &[S]) -> String {
    let items: Vec<String> = xs.iter().map(|x| pe(x.as_ref())).collect();
    format!("[{}:{}]", xs.len(), items.join(","))
}

/// `-` for an empty list, else its items' hex joined by `,`.
fn texts(s: &str) -> io::Result<Vec<String>> {
    if s == "-" {
        Ok(Vec::new())
    } else {
        s.split(',').map(text).collect()
    }
}

/// `-` for None, else the hex of the text.
fn opt_text(s: &str) -> io::Result<Option<String>> {
    if s == "-" {
        Ok(None)
    } else {
        text(s).map(Some)
    }
}

fn opt_int(s: &str) -> io::Result<Option<isize>> {
    if s == "-" { Ok(None) } else { int(s).map(Some) }
}

fn delimiter(s: &str) -> io::Result<char> {
    let d = text(s)?;
    let mut cs = d.chars();
    match (cs.next(), cs.next()) {
        (Some(c), None) => Ok(c),
        _ => Err(bad(format!("bad delimiter {d:?}"))),
    }
}

/// The rows a reader gives, then `!<exception>` if it stopped on one.
fn csv_result(op: &str, d: char, quoting: Quoting, lines: Vec<String>) -> String {
    let lines = lines.into_iter().map(Ok::<_, io::Error>);
    let mut out = String::new();
    let run = |out: &mut String| -> Result<(), csv::Error> {
        if op == "csvr" {
            for row in Reader::new(lines, d, quoting) {
                out.push_str(&render_list(&row?));
            }
            return Ok(());
        }
        let mut reader = DictReader::new(lines, d, quoting);
        match reader.header()? {
            None => out.push('~'),
            Some(h) => out.push_str(&render_list(h.names())),
        }
        for row in reader {
            let row = row?;
            let mut items: Vec<String> = row
                .iter()
                .map(|(k, v)| format!("{}={}", pe(k), v.map_or("~".to_owned(), pe)))
                .collect();
            if let Some(extra) = row.extra() {
                items.push(format!("~={}", render_list(extra)));
            }
            out.push_str(&format!("{{{}}}", items.join(",")));
        }
        Ok(())
    };
    if let Err(e) = run(&mut out) {
        out.push_str(&format!("!{}", e.exception()));
    }
    out
}

/// Per code point from `start`, `count` of them (surrogates skipped): the
/// isspace, isdigit, isdecimal and isalpha bits as a hex digit, then
/// `[upper,lower]`.
fn chars(start: u32, count: u32) -> String {
    let mut out = String::new();
    for c in (start..(start + count).min(0x110000)).filter_map(char::from_u32) {
        let f = [
            text::is_space,
            text::is_digit,
            text::is_decimal,
            text::is_alpha,
        ]
        .iter()
        .enumerate()
        .map(|(n, f)| u8::from(f(c)) << n)
        .sum::<u8>();
        let s = c.to_string();
        out.push_str(&format!(
            "{f:X}[{},{}]",
            pe(&text::upper(&s)),
            pe(&text::lower(&s))
        ));
    }
    out
}

/// One case's result.
fn eval(f: &[&str]) -> io::Result<String> {
    let arg = |i: usize| {
        f.get(i)
            .copied()
            .ok_or_else(|| bad(format!("{}: missing argument {i}", f[0])))
    };
    Ok(match f[0] {
        "repr" => fmt::repr(float(arg(1)?)?),
        "f" => fmt::fixed(float(arg(1)?)?, int(arg(2)?)?),
        "e" => fmt::sci(float(arg(1)?)?, int(arg(2)?)?),
        "g" => fmt::general(float(arg(1)?)?, int(arg(2)?)?),
        "pct" => fmt::percent(float(arg(1)?)?, int(arg(2)?)?),
        "grpf" => fmt::fixed_grouped(float(arg(1)?)?, int(arg(2)?)?),
        "grp" => fmt::int_grouped(int(arg(1)?)?),
        "zpad" => fmt::zero_pad(int(arg(1)?)?, int(arg(2)?)?),
        "ljust" => fmt::ljust(&text(arg(1)?)?, int(arg(2)?)?),
        "rjust" => fmt::rjust(&text(arg(1)?)?, int(arg(2)?)?),
        "round0" => match fmt::round_int(float(arg(1)?)?) {
            Ok(v) => fmt::int_str(v),
            Err(e) => e.exception().to_owned(),
        },
        "round" => match fmt::round(float(arg(1)?)?, int(arg(2)?)?) {
            Ok(v) => fmt::repr(v),
            Err(e) => e.exception().to_owned(),
        },
        "float" => match num::parse_float(&text(arg(1)?)?) {
            Ok(v) => fmt::repr(v),
            Err(e) => e.exception().to_owned(),
        },
        "int" => match num::parse_int(&text(arg(1)?)?) {
            Ok(v) => v.to_string(),
            Err(e) => e.exception().to_owned(),
        },
        "fsum" => match num::fsum(floats(arg(1)?)?) {
            Ok(v) => fmt::repr(v),
            Err(e) => e.exception().to_owned(),
        },
        "sum" => fmt::repr(num::sum(floats(arg(1)?)?)),
        "json" => {
            let opts = dump_opts(arg(2)?)?;
            match json::loads(&text(arg(1)?)?).and_then(|v| json::dumps(&v, &opts)) {
                Ok(s) => s,
                Err(e) => e.exception().to_owned(),
            }
        }
        "truthy" => match json::loads(&text(arg(1)?)?) {
            Ok(v) => py_bool(value::truthy(&v)),
            Err(e) => e.exception().to_owned(),
        },
        "vstr" | "vrepr" => {
            let v = json::loads(&text(arg(1)?)?).map_err(|e| bad(e.to_string()))?;
            let r = if f[0] == "vstr" {
                value::py_str(&v)
            } else {
                value::py_repr(&v)
            };
            r.unwrap_or_else(|e| e.exception().to_owned())
        }
        "getor" => {
            let load = |s: &str| -> io::Result<json::Value> {
                json::loads(&text(s)?).map_err(|e| bad(e.to_string()))
            };
            let (d, default) = (load(arg(1)?)?, load(arg(3)?)?);
            match value::get_or(&d, &text(arg(2)?)?, &default) {
                Ok(v) => json::dumps(v, &json::Dump::default())
                    .unwrap_or_else(|e| e.exception().to_owned()),
                Err(e) => e.exception().to_owned(),
            }
        }
        "fnmatch" => py_bool(mtime::fnmatch(os(arg(1)?)?, os(arg(2)?)?)),
        "escape" => mtime::escape(&text(arg(1)?)?),
        "sortnames" => {
            let mut names = os_list(arg(1)?)?;
            mtime::sort(&mut names);
            let hexes: Vec<String> = names.iter().map(|n| hex(n.as_bytes())).collect();
            hexes.join(",")
        }
        "glob" => {
            let found = mtime::glob(os(arg(1)?)?, arg(2)? == "1");
            let hexes: Vec<String> = found
                .iter()
                .map(|p| hex(p.as_os_str().as_bytes()))
                .collect();
            hexes.join("|")
        }
        "first" => maybe_path(files(&os_list(arg(1)?)?, true).into_iter().next()),
        "newest" => match mtime::newest(files(&os_list(arg(1)?)?, false)) {
            Ok(p) => maybe_path(p),
            Err(e) => os_exception(&e).to_owned(),
        },
        "mtime" => match mtime::mtime(os(arg(1)?)?) {
            Ok(t) => fmt::repr(t),
            Err(e) => os_exception(&e).to_owned(),
        },
        "mtimens" => match mtime::mtime_ns(os(arg(1)?)?) {
            Ok(t) => t.to_string(),
            Err(e) => os_exception(&e).to_owned(),
        },
        "lines" => {
            let lines: Vec<String> = Lines::new(io::Cursor::new(bytes(arg(1)?)?), Utf8::Replace)
                .collect::<io::Result<_>>()?;
            return Ok(lines
                .iter()
                .map(|l| escape(l))
                .collect::<Vec<_>>()
                .join("|"));
        }
        // readline() keeps each line's ending, so a missing or extra "\n" shows.
        "readline" => {
            let mut l = Lines::new(io::Cursor::new(bytes(arg(1)?)?), Utf8::Replace);
            let mut line = String::new();
            let mut all = Vec::new();
            while l.read_line(&mut line)? {
                all.push(escape(&line));
            }
            return Ok(all.join("|"));
        }
        "csvw" => {
            let mut w = csv::Writer::new(Vec::new(), delimiter(arg(1)?)?);
            w.write_row(texts(arg(2)?)?)?;
            String::from_utf8(w.into_inner()).map_err(|_| bad("csv wrote bad UTF-8"))?
        }
        "csvr" | "dictr" | "dictf" => {
            let quoting = match arg(2)? {
                "m" => Quoting::Minimal,
                "n" => Quoting::None,
                q => return Err(bad(format!("bad quoting {q:?}"))),
            };
            let lines = if f[0] == "dictf" {
                // A file read in Python's text mode, endings kept.
                let mut l = Lines::new(io::Cursor::new(bytes(arg(3)?)?), Utf8::Replace);
                let (mut line, mut all) = (String::new(), Vec::new());
                while l.read_line(&mut line)? {
                    all.push(line.clone());
                }
                all
            } else {
                texts(arg(3)?)?
            };
            csv_result(f[0], delimiter(arg(1)?)?, quoting, lines)
        }
        "split" => {
            let s = text(arg(1)?)?;
            let sep = opt_text(arg(2)?)?;
            let n: isize = int(arg(3)?)?;
            render_list(&text::split(&s, sep.as_deref(), usize::try_from(n).ok()))
        }
        "strip" | "lstrip" | "rstrip" => {
            let s = text(arg(1)?)?;
            let chars = opt_text(arg(2)?)?;
            let strip = match f[0] {
                "strip" => text::strip,
                "lstrip" => text::lstrip,
                _ => text::rstrip,
            };
            strip(&s, chars.as_deref()).to_owned()
        }
        "upper" => text::upper(&text(arg(1)?)?),
        "lower" => text::lower(&text(arg(1)?)?),
        "is" => {
            let s = text(arg(1)?)?;
            [text::isspace, text::isdigit, text::isdecimal, text::isalpha]
                .iter()
                .map(|f| if f(&s) { '1' } else { '0' })
                .collect()
        }
        "html" => text::html_escape(&text(arg(1)?)?).into_owned(),
        "chars" => chars(int(arg(1)?)?, int(arg(2)?)?),
        "slice" => text::slice(&text(arg(1)?)?, opt_int(arg(2)?)?, opt_int(arg(3)?)?).to_owned(),
        "median" | "fmean" | "variance" => {
            let stat: fn(Vec<f64>) -> Result<f64, stats::Error> = match f[0] {
                "median" => stats::median,
                "fmean" => stats::fmean,
                _ => stats::variance,
            };
            stat(floats(arg(1)?)?).map_or_else(|e| e.exception().to_owned(), fmt::repr)
        }
        "mediani" => {
            let xs = arg(1)?;
            let xs: Vec<i64> = if xs.is_empty() {
                Vec::new()
            } else {
                xs.split(',').map(int).collect::<io::Result<_>>()?
            };
            stats::median_int(xs).map_or_else(|e| e.exception().to_owned(), |m| m.to_string())
        }
        op => return Err(bad(format!("unknown operation {op:?}"))),
    })
}

pub fn run(input: impl BufRead, output: impl Write) -> io::Result<()> {
    let mut out = io::BufWriter::new(output);
    for line in input.lines() {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();
        let result = eval(&fields)?;
        // A lines result is escaped line by line already.
        if matches!(fields[0], "lines" | "readline") {
            writeln!(out, "{result}")?;
        } else {
            writeln!(out, "{}", escape(&result))?;
        }
    }
    out.flush()
}

/// `og pycompat-probe` against what CPython (`python3 -I` in PYTHON_IMAGE)
/// answered: 1000 cases per operation, drawn once by the oracle of the
/// commit that recorded them, frozen in `tests/rust/vectors/pycompat` (see
/// its README). `cases.tsv.gz` holds each case line and the answer line
/// after it; a `chars` answer is `@CHARS@`, spelled out from the table
/// `chars.tsv.gz` (`R lo hi F`: a run of code points that are their own
/// upper and lower case and need no escaping, flags F; `X cp F upper lower`
/// for the rest). The glob, first, newest, mtime and mtimens cases walked a
/// tree under `/tree`, listed in `tree.tsv.gz` (GNU find: type, path, link
/// target, mtime); the test builds it again and runs those cases on Linux
/// only (a case-insensitive file system cannot hold it). Directory order
/// belongs to the file system, so glob answers are compared as sorted
/// lists, and a `newest` answer may be any file tied with CPython's on
/// mtime.
#[cfg(test)]
mod vectors {
    use super::{files, hex, os_list, run};
    use std::collections::HashMap;
    use std::ffi::OsStr;
    use std::io::BufRead;
    use std::os::unix::ffi::OsStrExt;
    use std::path::{Path, PathBuf};
    use std::time::{Duration, SystemTime, UNIX_EPOCH};

    const TREE_OPS: &[&str] = &["glob", "first", "newest", "mtime", "mtimens"];
    /// "/tree", where the oracle built its tree, in hex.
    const TREE_HEX: &str = "2f74726565";

    fn dir() -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/rust/vectors/pycompat")
    }

    fn gz_lines(name: &str) -> Vec<Vec<u8>> {
        og_io::open(dir().join(name))
            .unwrap_or_else(|e| panic!("{name}: {e}"))
            .split(b'\n')
            .map(Result::unwrap)
            .collect()
    }

    /// The `chars` answer of each code point.
    fn chars_table() -> HashMap<u32, String> {
        let mut table = HashMap::new();
        for line in gz_lines("chars.tsv.gz") {
            let line = String::from_utf8(line).unwrap();
            let f: Vec<&str> = line.split('\t').collect();
            match f.as_slice() {
                ["R", lo, hi, flags] => {
                    for cp in lo.parse::<u32>().unwrap()..=hi.parse().unwrap() {
                        let c = char::from_u32(cp).unwrap();
                        assert!(table.insert(cp, format!("{flags}[{c},{c}]")).is_none());
                    }
                }
                ["X", cp, flags, upper, lower] => {
                    let entry = format!("{flags}[{upper},{lower}]");
                    assert!(table.insert(cp.parse().unwrap(), entry).is_none());
                }
                _ => panic!("chars.tsv: {line:?}"),
            }
        }
        table
    }

    /// `/tree` of a path in `bytes` written as `root`.
    fn at_root(bytes: &[u8], root: &[u8]) -> Vec<u8> {
        let tree = b"/tree";
        if bytes.starts_with(tree) && matches!(bytes.get(tree.len()), None | Some(b'/')) {
            [root, &bytes[tree.len()..]].concat()
        } else {
            bytes.to_vec()
        }
    }

    /// Each hex path item of `text` (items start the text or follow a tab,
    /// `,` or `|`) under `/tree` moved under the hex root `root`.
    fn swap_root(text: &str, root: &str) -> String {
        let mut out = String::with_capacity(text.len());
        let mut start = true;
        let mut rest = text;
        while !rest.is_empty() {
            if start && rest.starts_with(TREE_HEX) {
                let after = &rest[TREE_HEX.len()..];
                if after.is_empty()
                    || after.starts_with("2f")
                    || after.starts_with(['\t', ',', '|', '\\'])
                {
                    out.push_str(root);
                    rest = after;
                    start = false;
                    continue;
                }
            }
            let c = rest.chars().next().unwrap();
            out.push(c);
            start = matches!(c, '\t' | ',' | '|');
            rest = &rest[c.len_utf8()..];
        }
        out
    }

    /// `%T@` of GNU find (seconds, ten decimals, maybe negative) as a time.
    fn find_time(s: &str) -> SystemTime {
        // GNU find prints the floor second and a positive fraction, so
        // "-5.25" is -5 s + 0.25 s, not -5.25 s.
        let (sec, frac) = s.split_once('.').unwrap();
        let ns = i128::from(sec.parse::<i64>().unwrap()) * 1_000_000_000
            + i128::from(frac[..9].parse::<i64>().unwrap());
        assert!(
            frac[9..].bytes().all(|b| b == b'0'),
            "{s}: finer than a nanosecond"
        );
        let d = Duration::from_nanos(u64::try_from(ns.unsigned_abs()).unwrap());
        if ns < 0 {
            UNIX_EPOCH - d
        } else {
            UNIX_EPOCH + d
        }
    }

    /// The oracle's tree again under `root`: folders, empty files, symlinks
    /// (their targets moved under `root`), then each folder's and file's
    /// mtime, deepest first. A symlink's own mtime is never read (Python's
    /// getmtime and stat follow it), so it is left as made.
    fn build_tree(root: &Path) {
        let rows: Vec<Vec<Vec<u8>>> = gz_lines("tree.tsv.gz")
            .into_iter()
            .map(|l| l.split(|b| *b == b'\t').map(<[u8]>::to_vec).collect())
            .collect();
        let path = |rel: &[u8]| {
            if rel.is_empty() {
                root.to_path_buf()
            } else {
                root.join(OsStr::from_bytes(rel))
            }
        };
        let rb = root.as_os_str().as_bytes();
        for r in &rows {
            let [kind, rel, target, _] = r.as_slice() else {
                panic!("tree.tsv: {r:?}");
            };
            match (kind.as_slice(), rel.is_empty()) {
                (_, true) => {}
                (b"d", _) => std::fs::create_dir(path(rel)).unwrap(),
                (b"f", _) => drop(std::fs::File::create(path(rel)).unwrap()),
                (b"l", _) => {
                    std::os::unix::fs::symlink(OsStr::from_bytes(&at_root(target, rb)), path(rel))
                        .unwrap()
                }
                _ => panic!("tree.tsv: kind {kind:?}"),
            }
        }
        let mut timed: Vec<&Vec<Vec<u8>>> = rows.iter().filter(|r| r[0] != b"l").collect();
        let depth = |r: &Vec<Vec<u8>>| {
            if r[1].is_empty() {
                0
            } else {
                1 + r[1].iter().filter(|b| **b == b'/').count()
            }
        };
        timed.sort_by_key(|r| std::cmp::Reverse(depth(r)));
        for r in timed {
            let t = find_time(std::str::from_utf8(&r[3]).unwrap());
            let times = std::fs::FileTimes::new().set_accessed(t).set_modified(t);
            std::fs::File::open(path(&r[1]))
                .unwrap()
                .set_times(times)
                .unwrap();
        }
    }

    fn unhex(s: &str) -> PathBuf {
        let b: Vec<u8> = (0..s.len())
            .step_by(2)
            .map(|i| u8::from_str_radix(&s[i..i + 2], 16).unwrap())
            .collect();
        PathBuf::from(OsStr::from_bytes(&b))
    }

    /// Whether og's answer `got` stands for CPython's `want` on `case`.
    fn agrees(case: &str, want: &str, got: &str) -> bool {
        if want == got {
            return true;
        }
        let mut args = case.split('\t');
        match args.next().unwrap() {
            "glob" => {
                fn sorted(s: &str) -> Vec<&str> {
                    let mut v: Vec<&str> = s.split("\\|").collect();
                    v.sort_unstable();
                    v
                }
                sorted(want) == sorted(got)
            }
            "newest" if want != "None" && got != "None" && !got.contains('|') => {
                let t = |h: &str| pycompat::mtime::mtime(unhex(h)).ok();
                // A tie on mtime: any file the case's patterns match will do.
                let eligible = files(&os_list(args.next().unwrap_or("")).unwrap(), false);
                t(want).is_some() && t(want) == t(got) && eligible.contains(&unhex(got))
            }
            _ => false,
        }
    }

    #[test]
    fn og_pycompat_probe_matches_cpython_on_the_committed_cases() {
        let lines = gz_lines("cases.tsv.gz");
        assert!(
            lines.len() > 2 && lines.len().is_multiple_of(2),
            "cases.tsv: pairs of lines"
        );
        let chars = chars_table();
        let tree_ok = cfg!(target_os = "linux");
        let root = std::env::temp_dir()
            .canonicalize()
            .unwrap()
            .join(format!("og-pycompat-vectors-{}", std::process::id()));
        let root_s = root.to_str().unwrap();
        assert!(
            !root_s.contains(['*', '?', '[']),
            "{root_s}: a glob character"
        );
        let _ = std::fs::remove_dir_all(&root);
        std::fs::create_dir(&root).unwrap();
        if tree_ok {
            build_tree(&root);
        }
        let root_hex = hex(root.as_os_str().as_bytes());
        let (mut cases, mut want) = (Vec::new(), Vec::new());
        let mut skipped = 0;
        for pair in lines.chunks(2) {
            let case = String::from_utf8(pair[0].clone()).unwrap();
            let answer = String::from_utf8(pair[1].clone()).unwrap();
            let op = case.split('\t').next().unwrap();
            if TREE_OPS.contains(&op) {
                if !tree_ok {
                    skipped += 1;
                    continue;
                }
                cases.push(swap_root(&case, &root_hex));
                want.push(swap_root(&answer, &root_hex));
            } else if op == "chars" {
                assert_eq!(answer, "@CHARS@", "{case}");
                let f: Vec<u32> = case
                    .split('\t')
                    .skip(1)
                    .map(|x| x.parse().unwrap())
                    .collect();
                let line: String = (f[0]..(f[0] + f[1]).min(0x11_0000))
                    .filter(|c| !(0xD800..0xE000).contains(c))
                    .map(|c| chars[&c].as_str())
                    .collect();
                cases.push(case);
                want.push(line);
            } else {
                cases.push(case);
                want.push(answer);
            }
        }
        let mut input = cases.join("\n");
        input.push('\n');
        let mut out = Vec::new();
        run(input.as_bytes(), &mut out).unwrap();
        let got: Vec<&str> = std::str::from_utf8(&out).unwrap().lines().collect();
        assert_eq!(got.len(), cases.len(), "one answer per case");
        let mut bad: std::collections::BTreeMap<&str, usize> = Default::default();
        let mut report = String::new();
        for ((case, want), got) in cases.iter().zip(&want).zip(&got) {
            let op = case.split('\t').next().unwrap();
            if !agrees(case, want, got) {
                *bad.entry(op).or_default() += 1;
                if bad.values().sum::<usize>() <= 20 {
                    report.push_str(&format!(
                        "case:     {case}\npython:   {want}\nog:       {got}\n\n"
                    ));
                }
            }
        }
        let _ = std::fs::remove_dir_all(&root);
        println!(
            "pycompat: {} cases against CPython{}",
            cases.len(),
            if skipped > 0 {
                format!("; {skipped} tree cases need Linux")
            } else {
                String::new()
            }
        );
        assert!(
            bad.is_empty(),
            "{report}og pycompat-probe differs from CPython: {bad:?}"
        );
    }
}
