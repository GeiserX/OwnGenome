//! `og dev parity-compare KIND ORACLE NEW`: compares a file the Python
//! oracle wrote with the one its port wrote, by meaning.
//!
//! Exit 0 when they are the same, 1 when they differ (the first difference
//! on stdout), 2 when the kind is unknown or a file cannot be read.
//!
//! | Kind | Same when |
//! |---|---|
//! | `bytes` | the bytes are equal |
//! | `gz-content` | both are gzip and their decompressed bytes are equal |
//! | `content` | as `gz-content` when the expected file is gzip, else as `bytes`: a case whose outputs mix both |
//! | `zip-content` | the same member names with the same decompressed bytes (timestamps and deflate bytes may differ) |
//! | `png-pixels` | the same size and the same RGBA pixels |
//! | `tsv[:keys=A,B][:unordered]` | `#` lines equal and in order, the header equal, the rows equal cell by cell (csv quotes removed), in order unless `unordered`, or matched by the `keys` columns |
//! | `kv` | the same `key<TAB>value` rows in the same order |
//! | `vcf` | header lines equal and in order, records equal field by field, INFO as an ordered list of keys |
//! | `json[:paths=A,B][:mask=K,L][:tol=X]` | the same values: key order, int against float and floats compared exactly; `paths` limits it to those subtrees, `mask` drops keys of those names anywhere, `tol` lets floats differ by up to X |
//! | `mask:REGEX` | the lines are equal once every line REGEX matches is replaced by one marker |
//! | `stderr-first-line` | the first lines are equal |
//!
//! Text kinds read plain or gzip files (sniffed), with universal newlines.

use crate::zip;
use og_io::{Lines, Utf8};
use pycompat::json::{self, Value};
use std::collections::BTreeMap;
use std::fs;
use std::io::{BufRead, Read};
use std::path::{Path, PathBuf};

/// A comparison kind with its options.
#[derive(Debug)]
pub enum Kind {
    Bytes,
    GzContent,
    Content,
    ZipContent,
    PngPixels,
    Tsv {
        keys: Option<Vec<String>>,
        unordered: bool,
    },
    Kv,
    Vcf,
    Json {
        paths: Vec<String>,
        mask: Vec<String>,
        tol: Option<f64>,
    },
    Mask(regex::Regex),
    StderrFirstLine,
}

impl Kind {
    /// Reads `bytes`, `tsv:keys=a,b:unordered`, `mask:REGEX` and the rest.
    pub fn parse(spec: &str) -> Result<Kind, String> {
        if let Some(re) = spec.strip_prefix("mask:") {
            return regex::Regex::new(re)
                .map(Kind::Mask)
                .map_err(|e| format!("mask: {e}"));
        }
        let mut parts = spec.split(':');
        let name = parts.next().unwrap_or_default();
        let opts: Vec<&str> = parts.collect();
        let no_opts = |k: Kind| {
            if opts.is_empty() {
                Ok(k)
            } else {
                Err(format!("{name} takes no options"))
            }
        };
        match name {
            "bytes" => no_opts(Kind::Bytes),
            "gz-content" => no_opts(Kind::GzContent),
            "content" => no_opts(Kind::Content),
            "zip-content" => no_opts(Kind::ZipContent),
            "png-pixels" => no_opts(Kind::PngPixels),
            "kv" => no_opts(Kind::Kv),
            "vcf" => no_opts(Kind::Vcf),
            "stderr-first-line" => no_opts(Kind::StderrFirstLine),
            "tsv" => {
                let (mut keys, mut unordered) = (None, false);
                for o in opts {
                    match o.split_once('=') {
                        Some(("keys", v)) if !v.is_empty() => {
                            keys = Some(v.split(',').map(str::to_owned).collect());
                        }
                        None if o == "unordered" => unordered = true,
                        _ => return Err(format!("tsv: unknown option '{o}'")),
                    }
                }
                Ok(Kind::Tsv { keys, unordered })
            }
            "json" => {
                let (mut paths, mut mask, mut tol) = (Vec::new(), Vec::new(), None);
                for o in opts {
                    match o.split_once('=') {
                        Some(("paths", v)) if !v.is_empty() => {
                            paths = v.split(',').map(str::to_owned).collect();
                        }
                        Some(("mask", v)) if !v.is_empty() => {
                            mask = v.split(',').map(str::to_owned).collect();
                        }
                        Some(("tol", v)) => {
                            let t: f64 = v
                                .parse()
                                .ok()
                                .filter(|t: &f64| t.is_finite() && *t >= 0.0)
                                .ok_or(format!("json: tol '{v}' is not a number of 0 or more"))?;
                            tol = Some(t);
                        }
                        _ => return Err(format!("json: unknown option '{o}'")),
                    }
                }
                Ok(Kind::Json { paths, mask, tol })
            }
            _ => Err(format!("unknown kind '{spec}'")),
        }
    }
}

/// `Ok(None)` when the files are the same, `Ok(Some(why))` when they
/// differ, `Err` when one cannot be read as `kind` says.
pub type Outcome = Result<Option<String>, String>;

/// Compares `oracle` with `new` as `kind` says.
pub fn compare(kind: &Kind, oracle: &Path, new: &Path) -> Outcome {
    match kind {
        Kind::Bytes => bytes(&mut open_raw(oracle)?, &mut open_raw(new)?),
        Kind::GzContent => {
            for p in [oracle, new] {
                let mut magic = [0u8; 2];
                let n = open_raw(p)?.read(&mut magic).map_err(|e| err(p, e))?;
                if n < 2 || magic != [0x1f, 0x8b] {
                    return Ok(Some(format!("{} is not gzip", p.display())));
                }
            }
            bytes(&mut open_text(oracle)?, &mut open_text(new)?)
        }
        Kind::Content => {
            let mut magic = [0u8; 2];
            let n = open_raw(oracle)?
                .read(&mut magic)
                .map_err(|e| err(oracle, e))?;
            let kind = if n == 2 && magic == [0x1f, 0x8b] {
                Kind::GzContent
            } else {
                Kind::Bytes
            };
            compare(&kind, oracle, new)
        }
        Kind::ZipContent => zip_content(&read(oracle)?, &read(new)?),
        Kind::PngPixels => Ok(differ(&png_pixels(oracle)?, &png_pixels(new)?, "pixels")),
        Kind::Tsv { keys, unordered } => {
            tsv(&lines(oracle)?, &lines(new)?, keys.as_deref(), *unordered)
        }
        Kind::Kv => Ok(list(&kv(&lines(oracle)?), &kv(&lines(new)?), "row")),
        Kind::Vcf => Ok(vcf(&lines(oracle)?, &lines(new)?)),
        Kind::Json { paths, mask, tol } => {
            let parse = |p: &Path| -> Result<Value, String> {
                let text = String::from_utf8(read(p)?).map_err(|e| err(p, e))?;
                let mut v =
                    json::loads(&text).map_err(|e| err(p, format!("{} ({e:?})", e.exception())))?;
                strip_keys(&mut v, mask);
                Ok(v)
            };
            let (a, b) = (parse(oracle)?, parse(new)?);
            if paths.is_empty() {
                return Ok(json_diff(&a, &b, "$", *tol));
            }
            for p in paths {
                let (x, y) = (at_path(&a, p), at_path(&b, p));
                let here = format!("$.{p}");
                match (x, y) {
                    (None, None) => {}
                    (Some(x), Some(y)) => {
                        if let Some(d) = json_diff(x, y, &here, *tol) {
                            return Ok(Some(d));
                        }
                    }
                    (x, _) => {
                        let side = if x.is_none() { "oracle" } else { "new" };
                        return Ok(Some(format!("{here}: not in the {side} file")));
                    }
                }
            }
            Ok(None)
        }
        Kind::Mask(re) => {
            let mark = |ls: Vec<String>| -> Vec<String> {
                ls.into_iter()
                    .map(|l| {
                        if re.is_match(&l) {
                            "\u{0}masked".to_owned()
                        } else {
                            l
                        }
                    })
                    .collect()
            };
            Ok(list(&mark(lines(oracle)?), &mark(lines(new)?), "line"))
        }
        Kind::StderrFirstLine => {
            let first = |p: &Path| -> Result<String, String> {
                Ok(lines(p)?.into_iter().next().unwrap_or_default())
            };
            let (a, b) = (first(oracle)?, first(new)?);
            Ok((a != b).then(|| format!("first line: oracle {a:?}, new {b:?}")))
        }
    }
}

fn err(p: &Path, e: impl std::fmt::Display) -> String {
    format!("{}: {e}", p.display())
}

fn open_raw(p: &Path) -> Result<std::io::BufReader<fs::File>, String> {
    let f = fs::File::open(p).map_err(|e| err(p, e))?;
    Ok(std::io::BufReader::with_capacity(1 << 16, f))
}

fn open_text(p: &Path) -> Result<Box<dyn BufRead>, String> {
    og_io::open(p).map_err(|e| err(p, e))
}

fn read(p: &Path) -> Result<Vec<u8>, String> {
    fs::read(p).map_err(|e| err(p, e))
}

/// Every line, plain or gzip, universal newlines, strict UTF-8.
fn lines(p: &Path) -> Result<Vec<String>, String> {
    let mut r = Lines::new(open_text(p)?, Utf8::Strict);
    let mut out = Vec::new();
    let mut l = String::new();
    while r.read_line(&mut l).map_err(|e| err(p, e))? {
        out.push(l.strip_suffix('\n').unwrap_or(&l).to_owned());
    }
    Ok(out)
}

/// Streams both readers; the first differing byte, with its line.
fn bytes(a: &mut dyn BufRead, b: &mut dyn BufRead) -> Outcome {
    let (mut offset, mut line) = (0u64, 1u64);
    loop {
        let x = a.fill_buf().map_err(|e| e.to_string())?;
        let y = b.fill_buf().map_err(|e| e.to_string())?;
        if x.is_empty() || y.is_empty() {
            return Ok(match (x.is_empty(), y.is_empty()) {
                (true, true) => None,
                (true, false) => Some(format!(
                    "the oracle ends at byte {offset} (line {line}); new goes on"
                )),
                _ => Some(format!(
                    "new ends at byte {offset} (line {line}); the oracle goes on"
                )),
            });
        }
        let n = x.len().min(y.len());
        if let Some(i) = (0..n).find(|&i| x[i] != y[i]) {
            line += x[..i].iter().filter(|&&c| c == b'\n').count() as u64;
            return Ok(Some(format!(
                "byte {} (line {line}): oracle {:?}, new {:?}",
                offset + i as u64,
                char::from(x[i]),
                char::from(y[i])
            )));
        }
        line += x[..n].iter().filter(|&&c| c == b'\n').count() as u64;
        offset += n as u64;
        a.consume(n);
        b.consume(n);
    }
}

fn differ<T: PartialEq>(a: &T, b: &T, what: &str) -> Option<String> {
    (a != b).then(|| format!("the {what} differ"))
}

/// Ordered lists: the first index where they differ, or the shorter one.
fn list<T: PartialEq + std::fmt::Debug>(a: &[T], b: &[T], what: &str) -> Option<String> {
    for (i, (x, y)) in a.iter().zip(b).enumerate() {
        if x != y {
            return Some(format!("{what} {}: oracle {x:?}, new {y:?}", i + 1));
        }
    }
    match a.len().cmp(&b.len()) {
        std::cmp::Ordering::Equal => None,
        std::cmp::Ordering::Less => Some(format!(
            "new has {} {what}s more: first {:?}",
            b.len() - a.len(),
            b[a.len()]
        )),
        std::cmp::Ordering::Greater => Some(format!(
            "new misses {} {what}s: first {:?}",
            a.len() - b.len(),
            a[b.len()]
        )),
    }
}

fn zip_content(a: &[u8], b: &[u8]) -> Outcome {
    let (x, y) = (zip::members(a)?, zip::members(b)?);
    let names = |m: &BTreeMap<String, Vec<u8>>| m.keys().cloned().collect::<Vec<_>>();
    if x.len() != y.len() || x.keys().ne(y.keys()) {
        return Ok(Some(format!(
            "members: oracle {:?}, new {:?}",
            names(&x),
            names(&y)
        )));
    }
    for (name, content) in &x {
        if let Some(d) = bytes(&mut content.as_slice(), &mut y[name].as_slice())? {
            return Ok(Some(format!("member {name}: {d}")));
        }
    }
    Ok(None)
}

/// Width, height and every pixel as 16-bit RGBA (8-bit samples scaled by 257).
fn png_pixels(p: &Path) -> Result<(u32, u32, Vec<u16>), String> {
    let mut dec = png::Decoder::new(open_raw(p)?);
    dec.set_transformations(png::Transformations::EXPAND);
    let mut r = dec.read_info().map_err(|e| err(p, e))?;
    let mut buf = vec![0; r.output_buffer_size().ok_or_else(|| err(p, "too large"))?];
    let info = r.next_frame(&mut buf).map_err(|e| err(p, e))?;
    let buf = &buf[..info.buffer_size()];
    let samples: Vec<u16> = match info.bit_depth {
        png::BitDepth::Sixteen => buf
            .chunks_exact(2)
            .map(|c| u16::from_be_bytes([c[0], c[1]]))
            .collect(),
        _ => buf.iter().map(|&v| u16::from(v) * 257).collect(),
    };
    let rgba: Vec<u16> = match info.color_type {
        png::ColorType::Grayscale => samples.iter().flat_map(|&g| [g, g, g, u16::MAX]).collect(),
        png::ColorType::GrayscaleAlpha => samples
            .chunks_exact(2)
            .flat_map(|c| [c[0], c[0], c[0], c[1]])
            .collect(),
        png::ColorType::Rgb => samples
            .chunks_exact(3)
            .flat_map(|c| [c[0], c[1], c[2], u16::MAX])
            .collect(),
        png::ColorType::Rgba => samples,
        png::ColorType::Indexed => return Err(err(p, "a palette PNG that did not expand")),
    };
    Ok((info.width, info.height, rgba))
}

/// A cell as Python's csv module reads a quoted one: outer quotes removed,
/// doubled quotes made single.
fn unquote(cell: &str) -> String {
    match cell.strip_prefix('"').and_then(|c| c.strip_suffix('"')) {
        Some(inner) => inner.replace("\"\"", "\""),
        None => cell.to_owned(),
    }
}

fn tsv(a: &[String], b: &[String], keys: Option<&[String]>, unordered: bool) -> Outcome {
    // `#` lines stay lines; the first other line is the header.
    let split = |ls: &[String]| -> (Vec<String>, Option<Vec<String>>, Vec<Vec<String>>) {
        let comments: Vec<String> = ls.iter().filter(|l| l.starts_with('#')).cloned().collect();
        let mut rest = ls.iter().filter(|l| !l.starts_with('#'));
        let header = rest.next().map(|h| h.split('\t').map(unquote).collect());
        let rows = rest.map(|l| l.split('\t').map(unquote).collect()).collect();
        (comments, header, rows)
    };
    let (ca, ha, mut ra) = split(a);
    let (cb, hb, mut rb) = split(b);
    if let Some(d) = list(&ca, &cb, "# line") {
        return Ok(Some(d));
    }
    if ha != hb {
        return Ok(Some(format!("header: oracle {ha:?}, new {hb:?}")));
    }
    let header = ha.unwrap_or_default();
    if let Some(keys) = keys {
        let idx: Vec<usize> = keys
            .iter()
            .map(|k| {
                header
                    .iter()
                    .position(|h| h == k)
                    .ok_or(format!("key column {k} is not in the header"))
            })
            .collect::<Result<_, _>>()?;
        let keyed = |rows: Vec<Vec<String>>| -> Result<BTreeMap<Vec<String>, Vec<String>>, String> {
            let mut m = BTreeMap::new();
            for r in rows {
                let k: Vec<String> = idx
                    .iter()
                    .map(|&i| r.get(i).cloned().unwrap_or_default())
                    .collect();
                if m.insert(k.clone(), r).is_some() {
                    return Err(format!("key {k:?} appears twice"));
                }
            }
            Ok(m)
        };
        let (ka, kb) = (keyed(ra)?, keyed(rb)?);
        for (k, r) in &ka {
            match kb.get(k) {
                None => return Ok(Some(format!("row {k:?}: not in new"))),
                Some(s) if s != r => {
                    return Ok(Some(row_diff(&header, r, s, &format!("row {k:?}"))));
                }
                _ => {}
            }
        }
        return Ok(kb
            .keys()
            .find(|k| !ka.contains_key(*k))
            .map(|k| format!("row {k:?}: not in the oracle")));
    }
    if unordered {
        ra.sort();
        rb.sort();
    }
    for (i, (x, y)) in ra.iter().zip(&rb).enumerate() {
        if x != y {
            return Ok(Some(row_diff(&header, x, y, &format!("row {}", i + 1))));
        }
    }
    Ok(list(&ra, &rb, "row"))
}

fn row_diff(header: &[String], x: &[String], y: &[String], here: &str) -> String {
    let i = (0..x.len().max(y.len()))
        .find(|&i| x.get(i) != y.get(i))
        .unwrap_or(0);
    let col = header
        .get(i)
        .cloned()
        .unwrap_or_else(|| format!("column {}", i + 1));
    format!("{here}, {col}: oracle {:?}, new {:?}", x.get(i), y.get(i))
}

fn kv(ls: &[String]) -> Vec<(String, String)> {
    ls.iter()
        .filter(|l| !l.is_empty())
        .map(|l| match l.split_once('\t') {
            Some((k, v)) => (k.to_owned(), v.to_owned()),
            None => (l.clone(), String::new()),
        })
        .collect()
}

fn vcf(a: &[String], b: &[String]) -> Option<String> {
    let head = |ls: &[String]| {
        ls.iter()
            .take_while(|l| l.starts_with('#'))
            .cloned()
            .collect::<Vec<_>>()
    };
    let (ha, hb) = (head(a), head(b));
    if let Some(d) = list(&ha, &hb, "header line") {
        return Some(d);
    }
    const FIELDS: [&str; 8] = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"];
    let (ra, rb) = (&a[ha.len()..], &b[hb.len()..]);
    for (i, (x, y)) in ra.iter().zip(rb).enumerate() {
        if x == y {
            continue;
        }
        let (fx, fy): (Vec<&str>, Vec<&str>) = (x.split('\t').collect(), y.split('\t').collect());
        let here = format!("record {} ({}:{})", i + 1, fx[0], fx.get(1).unwrap_or(&""));
        for j in 0..fx.len().max(fy.len()) {
            let name = FIELDS
                .get(j)
                .map_or_else(|| format!("column {}", j + 1), |s| (*s).to_owned());
            let (u, v) = (fx.get(j), fy.get(j));
            if u == v {
                continue;
            }
            if j == 7
                && let (Some(u), Some(v)) = (u, v)
            {
                let info = |s: &str| s.split(';').map(str::to_owned).collect::<Vec<_>>();
                let keys = |s: &str| {
                    s.split(';')
                        .map(|e| e.split('=').next().unwrap_or("").to_owned())
                        .collect::<Vec<_>>()
                };
                if keys(u) != keys(v) {
                    return Some(format!(
                        "{here}, INFO keys: oracle {:?}, new {:?}",
                        keys(u),
                        keys(v)
                    ));
                }
                return list(&info(u), &info(v), "INFO entry").map(|d| format!("{here}, {d}"));
            }
            return Some(format!("{here}, {name}: oracle {u:?}, new {v:?}"));
        }
    }
    list(ra, rb, "record")
}

/// Removes every object key named in `mask`, at any depth.
fn strip_keys(v: &mut Value, mask: &[String]) {
    match v {
        Value::Object(m) => {
            m.retain(|k, _| !mask.contains(k));
            m.values_mut().for_each(|x| strip_keys(x, mask));
        }
        Value::Array(a) => a.iter_mut().for_each(|x| strip_keys(x, mask)),
        _ => {}
    }
}

/// The value at a dotted path (`a.b.0`): object keys, or array indexes.
fn at_path<'a>(v: &'a Value, path: &str) -> Option<&'a Value> {
    path.split('.').try_fold(v, |v, step| match v {
        Value::Object(m) => m.get(step),
        Value::Array(a) => step.parse::<usize>().ok().and_then(|i| a.get(i)),
        _ => None,
    })
}

fn kind_of(v: &Value) -> &'static str {
    match v {
        Value::Null => "null",
        Value::Bool(_) => "bool",
        Value::Int(_) => "int",
        Value::Float(_) => "float",
        Value::Str(_) => "str",
        Value::Array(_) => "array",
        Value::Object(_) => "object",
    }
}

/// The first difference between two JSON values, by path.
fn json_diff(a: &Value, b: &Value, here: &str, tol: Option<f64>) -> Option<String> {
    match (a, b) {
        (Value::Float(x), Value::Float(y)) => {
            let same = x.to_bits() == y.to_bits()
                || (x.is_nan() && y.is_nan())
                || tol.is_some_and(|t| (x - y).abs() <= t);
            (!same).then(|| {
                format!(
                    "{here}: oracle {}, new {}",
                    pycompat::fmt::repr(*x),
                    pycompat::fmt::repr(*y)
                )
            })
        }
        (Value::Object(x), Value::Object(y)) => {
            if x.keys().ne(y.keys()) {
                let k =
                    |m: &indexmap::IndexMap<String, Value>| m.keys().cloned().collect::<Vec<_>>();
                return Some(format!("{here}: keys: oracle {:?}, new {:?}", k(x), k(y)));
            }
            x.iter()
                .find_map(|(k, v)| json_diff(v, &y[k], &format!("{here}.{k}"), tol))
        }
        (Value::Array(x), Value::Array(y)) => {
            if x.len() != y.len() {
                return Some(format!(
                    "{here}: oracle has {} items, new {}",
                    x.len(),
                    y.len()
                ));
            }
            x.iter()
                .zip(y)
                .enumerate()
                .find_map(|(i, (u, v))| json_diff(u, v, &format!("{here}[{i}]"), tol))
        }
        _ if kind_of(a) != kind_of(b) => Some(format!(
            "{here}: oracle {} {a:?}, new {} {b:?}",
            kind_of(a),
            kind_of(b)
        )),
        _ => (a != b).then(|| format!("{here}: oracle {a:?}, new {b:?}")),
    }
}

/// `og dev parity-compare ARGS`.
pub fn cli(args: &[String]) -> u8 {
    match args {
        [flag] if flag == "--self-test" => crate::parity::self_test::run(),
        [kind, oracle, new] => {
            let kind = match Kind::parse(kind) {
                Ok(k) => k,
                Err(e) => {
                    eprintln!("og dev parity-compare: {e}");
                    return 2;
                }
            };
            match compare(&kind, &PathBuf::from(oracle), &PathBuf::from(new)) {
                Ok(None) => 0,
                Ok(Some(d)) => {
                    println!("{d}");
                    1
                }
                Err(e) => {
                    eprintln!("og dev parity-compare: {e}");
                    2
                }
            }
        }
        _ => {
            eprintln!(
                "usage: og dev parity-compare KIND ORACLE NEW\n       og dev parity-compare --self-test"
            );
            2
        }
    }
}

pub mod self_test;
