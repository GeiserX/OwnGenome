//! `og dev pictures`: the docs pictures of DEMO-001, taken in a headless
//! Chrome over the DevTools protocol (`og-browser`), and their checks. It
//! replaces `tests/demo/render_pictures.py` and `tests/demo/check_pictures.py`.
//!
//! ```text
//! og dev pictures take --genome-dir DIR --out OUT [--only html|cpic]
//!                      [--html-stop-before TITLES] [--cpic-head N | --cpic-full]
//!                      [--terms-file FILE]
//! og dev pictures live --site SITE --out OUT [--phone] [--terms-file FILE]
//! og dev pictures check PNG...
//! og dev pictures check --self-test [--terms-file FILE]
//! og dev pictures same-sample REPORT LIVE_REPORT
//! ```
//!
//! `take` reads DIR/DEMO-001 after steps 24, 27 and 28 ran on it and writes
//! into OUT:
//!
//! - `demo-html-report.png`: the step 24 report in a 1180-pixel-wide window,
//!   from the top down to the first card whose title starts with one of
//!   TITLES (comma-separated; by default the first table that lists genes or
//!   variants, so no gene of the invented sample shows);
//! - `demo-multiqc-report.png`: the top 1600 x 1118 pixels of the step 28
//!   report, with MultiQC's first-visit bar closed;
//! - `demo-cpic-report.png`: the top of the step 27 text report as
//!   `head -n N` prints it in a 900-pixel-wide dark terminal, N ending on the
//!   last row of the Gene Results table (`--cpic-head N` sets N,
//!   `--cpic-full` draws the whole file as `cat` prints it).
//!
//! `live` takes the live-site demo SITE (as `tests/demo/build_live_demo.py
//! --out SITE` builds it), served from 127.0.0.1, at the last moment of its
//! replay: `demo-live-site.png` is DEMO-001's person page in a 1280 x 1000
//! dark window; with `--phone`, `demo-live-phone.png` is DEMO-002's person
//! page in a 390-pixel-wide dark window as tall as the page.
//!
//! Every picture is written from its pixels alone (IHDR, IDAT and IEND), so
//! it holds no text, time or EXIF chunk. Before it is written, the text the
//! picture shows (each text node of the page inside the picture's area, or
//! the CPIC lines drawn) must name the sample and match none of the private
//! terms: home and volume paths, private IPv4 addresses, tailnet names,
//! plus one case-insensitive regular expression per line of `--terms-file`
//! (keep that file outside the repository: it lists your own names and host
//! names). A finding still writes the picture, so you can look at it, and
//! exits 1.
//!
//! `check` reads each PNG's chunks: only IHDR, PLTE, IDAT and IEND pass. It
//! prints each picture's size. `--self-test` plants one bad case for each
//! check (a tEXt chunk, a home path, a private address, a tailnet name, a
//! text without the sample, and the first line of the terms file) and exits
//! 1 unless every one is caught.
//!
//! `same-sample` reads the Total variants line of two step 24 HTML reports,
//! the one the docs picture shows and the live demo's
//! `runs/DEMO-001/report.html`, and exits 1 when they differ: the docs call
//! the live demo the same sample.
//!
//! Exit 0 clean, 1 on a finding, 2 on a usage error or when a picture
//! cannot be taken.
//!
//! Where it differs from the Python on purpose: the text check reads the
//! page's own text where the Python read the picture with macOS Vision OCR,
//! so it runs on every system and in CI, and it reads every character
//! exactly. Text drawn into a canvas or an image is not read. The
//! terminal picture is drawn by Chrome from an HTML page in a monospace
//! font (Menlo, else DejaVu Sans Mono) instead of by Pillow, so there is no
//! `--font` option. The phone picture's width is checked here, not by the
//! workflow.

use std::collections::{BTreeMap, BTreeSet};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::LazyLock;
use std::thread::sleep;
use std::time::Duration;

use og_browser::{Browser, StaticServer, find_chrome, js_str};
use pycompat::json::Value;

const USAGE: &str = "usage: og dev pictures take --genome-dir DIR --out OUT [--only html|cpic] [--html-stop-before TITLES] [--cpic-head N | --cpic-full] [--terms-file FILE]\n       og dev pictures live --site SITE --out OUT [--phone] [--terms-file FILE]\n       og dev pictures check PNG...\n       og dev pictures check --self-test [--terms-file FILE]\n       og dev pictures same-sample REPORT LIVE_REPORT\n";

const SAMPLE: &str = "DEMO-001";
/// DEMO-001's invented partner, whose page the phone picture shows.
const PARTNER: &str = "DEMO-002";
/// Cards of the step 24 report that list genes or variants, in the order
/// the report has them. The picture stops above the first one present.
const STOP_BEFORE: &str =
    "Secondary-Findings Genes,Genes With a Non-Normal,ClinVar Hits,Steps Not Run";
/// The chunks a picture may hold: pixels and their header, nothing else.
const ALLOWED: [&[u8; 4]; 4] = [b"IHDR", b"PLTE", b"IDAT", b"IEND"];
const PNG_SIGNATURE: &[u8; 8] = b"\x89PNG\r\n\x1a\n";
/// How long a page may take to show what a picture needs.
const WAIT: Duration = Duration::from_secs(60);
const PHONE_WIDTH: i128 = 390;

/// Why a verb stopped: a usage error, or a picture that could not be taken.
#[derive(Debug)]
enum Fail {
    Usage(String),
    Error(String),
}

impl From<String> for Fail {
    fn from(e: String) -> Fail {
        Fail::Error(e)
    }
}

/// `og dev pictures ARGS`.
pub fn cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let rest = args.get(1..).unwrap_or_default();
    let r = match args.first().map(String::as_str) {
        Some("take") => take(rest, out),
        Some("live") => live(rest, out),
        Some("check") => check(rest, out),
        Some("same-sample") => same_sample(rest, out),
        Some("-h" | "--help") => {
            let _ = out.write_all(USAGE.as_bytes());
            return 0;
        }
        _ => Err(Fail::Usage("give take, live, check or same-sample".into())),
    };
    match r {
        Ok(code) => code,
        Err(Fail::Usage(m)) => {
            let _ = writeln!(err, "{USAGE}og dev pictures: error: {m}");
            2
        }
        Err(Fail::Error(m)) => {
            let _ = writeln!(err, "ERROR: {m}");
            2
        }
    }
}

/// The options of one verb: `--name value` (or `--name=value`) for the
/// names in `valued`, bare switches for those in `switches`, the rest free.
#[derive(Default)]
struct Opts {
    values: BTreeMap<String, String>,
    switches: BTreeSet<String>,
    free: Vec<String>,
}

impl Opts {
    fn parse(args: &[String], valued: &[&str], switches: &[&str]) -> Result<Opts, Fail> {
        let mut o = Opts::default();
        let mut it = args.iter();
        while let Some(a) = it.next() {
            if let Some((k, v)) = a.split_once('=').filter(|(k, _)| valued.contains(k)) {
                o.values.insert(k.to_owned(), v.to_owned());
            } else if valued.contains(&a.as_str()) {
                let v = it
                    .next()
                    .ok_or_else(|| Fail::Usage(format!("{a} needs a value")))?;
                o.values.insert(a.clone(), v.clone());
            } else if switches.contains(&a.as_str()) {
                o.switches.insert(a.clone());
            } else if a.starts_with("--") {
                return Err(Fail::Usage(format!("unknown option {a}")));
            } else {
                o.free.push(a.clone());
            }
        }
        Ok(o)
    }

    fn get(&self, k: &str) -> Option<&str> {
        self.values.get(k).map(String::as_str)
    }

    fn need(&self, k: &str) -> Result<&str, Fail> {
        self.get(k)
            .ok_or_else(|| Fail::Usage(format!("{k} is required")))
    }

    fn has(&self, k: &str) -> bool {
        self.switches.contains(k)
    }

    fn no_free(&self) -> Result<(), Fail> {
        match self.free.first() {
            Some(a) => Err(Fail::Usage(format!("unexpected argument {a}"))),
            None => Ok(()),
        }
    }
}

// ---------------------------------------------------------------- terms

/// One pattern the text of a picture must not match, and what it stands
/// for (a finding names the kind, never the match: the match may be the
/// private term itself).
struct Term {
    kind: &'static str,
    rx: fancy_regex::Regex,
}

const OCTET: &str = "(?:25[0-5]|2[0-4][0-9]|1[0-9][0-9]|[1-9]?[0-9])";

/// The built-in patterns: home and volume paths, private IPv4 addresses
/// (with the shared range a tailnet hands out), tailnet host names.
static BUILT_IN: LazyLock<Vec<(&'static str, String)>> = LazyLock::new(|| {
    vec![
        (
            "home or volume path",
            r"(?i)/(?:Users|home|Volumes|mnt/user)/".to_owned(),
        ),
        (
            "private address",
            format!(
                r"(?<![0-9.])(?:10\.{o}|172\.(?:1[6-9]|2[0-9]|3[01])|192\.168|100\.(?:6[4-9]|[7-9][0-9]|1[01][0-9]|12[0-7]))\.{o}\.{o}",
                o = OCTET
            ),
        ),
        ("tailnet name", r"(?i)\.ts\.net\b".to_owned()),
    ]
});

/// The built-in patterns plus every active line of `path` (not blank, not
/// a `#` comment), each case-insensitive. A file with no active line is
/// refused: a check that reads it would pass on nothing.
fn load_terms(path: Option<&str>) -> Result<Vec<Term>, Fail> {
    let mut terms = BUILT_IN
        .iter()
        .map(|(kind, p)| Term {
            kind,
            rx: fancy_regex::Regex::new(p).expect("the built-in patterns compile"),
        })
        .collect::<Vec<_>>();
    if let Some(path) = path {
        let text =
            std::fs::read_to_string(path).map_err(|e| Fail::Error(crate::py_os_error(&e, path)))?;
        for (n, line) in text.lines().enumerate() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            // The error names the line number only: the line is private.
            let rx = fancy_regex::Regex::new(&format!("(?i){line}")).map_err(|_| {
                Fail::Error(format!(
                    "line {} of the terms file {path} is not a regular expression",
                    n + 1
                ))
            })?;
            terms.push(Term {
                kind: "private term",
                rx,
            });
        }
        if terms.len() == BUILT_IN.len() {
            return Err(Fail::Error(format!(
                "the terms file {path} has no active pattern (only blank or # lines)"
            )));
        }
    }
    Ok(terms)
}

/// What is wrong with the text a picture shows: `expect` missing, or a
/// term found. The text is read as it is and with its white space folded,
/// so a term split over two lines is found too. A pattern that cannot be
/// evaluated counts as found.
fn text_problems(text: &str, expect: &str, terms: &[Term]) -> Vec<String> {
    let flat = text.split_whitespace().collect::<Vec<_>>().join(" ");
    let mut problems = Vec::new();
    let unhyphen = |s: &str| s.replace(['-', '\u{2014}'], "");
    if !expect.is_empty() && !flat.contains(expect) && !unhyphen(&flat).contains(&unhyphen(expect))
    {
        problems.push(format!("its text does not show {expect}"));
    }
    for (i, t) in terms.iter().enumerate() {
        if t.rx.is_match(text).unwrap_or(true) || t.rx.is_match(&flat).unwrap_or(true) {
            problems.push(format!(
                "its text holds a {} (pattern {} of the list)",
                t.kind,
                i + 1
            ));
        }
    }
    problems
}

// ---------------------------------------------------------------- PNG

/// The four-letter kind of every chunk of `data`, in file order.
fn chunk_kinds(data: &[u8]) -> Result<Vec<[u8; 4]>, String> {
    if !data.starts_with(PNG_SIGNATURE) {
        return Err("not a PNG".into());
    }
    let mut kinds = Vec::new();
    let mut i = PNG_SIGNATURE.len();
    while i < data.len() {
        let head = data.get(i..i + 8).ok_or("a chunk header is cut short")?;
        let n = u32::from_be_bytes([head[0], head[1], head[2], head[3]]) as usize;
        kinds.push([head[4], head[5], head[6], head[7]]);
        i = i
            .checked_add(12 + n)
            .filter(|&end| end <= data.len())
            .ok_or("a chunk is cut short")?;
    }
    Ok(kinds)
}

/// `chunk <kind>` for every kind of chunk outside [`ALLOWED`], sorted.
fn chunk_problems(data: &[u8]) -> Result<Vec<String>, String> {
    let bad: BTreeSet<String> = chunk_kinds(data)?
        .iter()
        .filter(|k| !ALLOWED.contains(k))
        .map(|k| String::from_utf8_lossy(k).into_owned())
        .collect();
    Ok(bad.into_iter().map(|k| format!("chunk {k}")).collect())
}

/// Width and height from the IHDR of a PNG whose chunks were read.
fn png_size(data: &[u8]) -> (u32, u32) {
    let be = |i: usize| u32::from_be_bytes([data[i], data[i + 1], data[i + 2], data[i + 3]]);
    (be(16), be(20))
}

/// `data` written again as an 8-bit RGB PNG from its pixels alone (IHDR,
/// IDAT, IEND), with its width and height. Alpha is dropped, as Pillow's
/// `convert("RGB")` drops it.
fn clean_png(data: &[u8]) -> Result<(Vec<u8>, u32, u32), String> {
    let mut dec = png::Decoder::new(std::io::Cursor::new(data));
    dec.set_transformations(png::Transformations::EXPAND | png::Transformations::STRIP_16);
    let mut r = dec
        .read_info()
        .map_err(|e| format!("a PNG Chrome wrote: {e}"))?;
    let mut buf = vec![0; r.output_buffer_size().ok_or("a PNG too large to read")?];
    let info = r
        .next_frame(&mut buf)
        .map_err(|e| format!("a PNG Chrome wrote: {e}"))?;
    let px = &buf[..info.buffer_size()];
    let rgb: Vec<u8> = match info.color_type {
        png::ColorType::Rgb => px.to_vec(),
        png::ColorType::Rgba => px
            .chunks_exact(4)
            .flat_map(|p| [p[0], p[1], p[2]])
            .collect(),
        png::ColorType::Grayscale => px.iter().flat_map(|&g| [g, g, g]).collect(),
        png::ColorType::GrayscaleAlpha => px.chunks_exact(2).flat_map(|p| [p[0]; 3]).collect(),
        png::ColorType::Indexed => return Err("a palette PNG that did not expand".into()),
    };
    let mut out = Vec::new();
    {
        let mut e = png::Encoder::new(&mut out, info.width, info.height);
        e.set_color(png::ColorType::Rgb);
        e.set_depth(png::BitDepth::Eight);
        e.set_compression(png::Compression::High);
        let mut w = e.write_header().map_err(|e| e.to_string())?;
        w.write_image_data(&rgb).map_err(|e| e.to_string())?;
        w.finish().map_err(|e| e.to_string())?;
    }
    Ok((out, info.width, info.height))
}

/// Standard base64 (RFC 4648), as `Page.captureScreenshot` sends a picture.
fn base64_decode(s: &str) -> Result<Vec<u8>, String> {
    let mut out = Vec::with_capacity(s.len() / 4 * 3);
    let (mut acc, mut bits) = (0u32, 0u32);
    for c in s.bytes() {
        let v = match c {
            b'A'..=b'Z' => c - b'A',
            b'a'..=b'z' => c - b'a' + 26,
            b'0'..=b'9' => c - b'0' + 52,
            b'+' => 62,
            b'/' => 63,
            b'=' => break,
            _ => return Err(format!("base64: unexpected byte {c:#04x}")),
        };
        acc = (acc << 6) | u32::from(v);
        bits += 6;
        if bits >= 8 {
            bits -= 8;
            out.push((acc >> bits) as u8);
            acc &= (1 << bits) - 1;
        }
    }
    Ok(out)
}

// ---------------------------------------------------------------- Chrome

/// A folder of this run, removed when it is dropped.
struct Scratch(PathBuf);

impl Scratch {
    fn new() -> Result<Scratch, String> {
        let d = std::env::temp_dir().join(format!("og-dev-pictures-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).map_err(|e| format!("{}: {e}", d.display()))?;
        Ok(Scratch(d))
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

fn obj<const N: usize>(kv: [(&str, Value); N]) -> Value {
    Value::Object(kv.into_iter().map(|(k, v)| (k.to_owned(), v)).collect())
}

fn int(n: i128) -> Value {
    Value::Int(n)
}

/// A headless Chrome with one page, its profile in `scratch`.
fn open_browser(scratch: &Scratch) -> Result<Browser, String> {
    let chrome = find_chrome().ok_or(
        "no Chrome found: install Google Chrome or Chromium, or set CHROME to its executable",
    )?;
    Browser::launch(&chrome, &scratch.0.join("chrome"))
}

/// A `width` x `height` window at one device pixel per CSS pixel, dark
/// when `dark` (the page's `prefers-color-scheme`).
fn window(b: &mut Browser, width: i128, height: i128, dark: bool) -> Result<(), String> {
    b.call(
        "Emulation.setDeviceMetricsOverride",
        obj([
            ("width", int(width)),
            ("height", int(height)),
            ("deviceScaleFactor", int(1)),
            ("mobile", Value::Bool(false)),
        ]),
    )?;
    let scheme = if dark { "dark" } else { "light" };
    b.call(
        "Emulation.setEmulatedMedia",
        obj([(
            "features",
            Value::Array(vec![obj([
                ("name", Value::Str("prefers-color-scheme".into())),
                ("value", Value::Str(scheme.into())),
            ])]),
        )]),
    )?;
    Ok(())
}

/// Open `url` and wait until the page and everything it loads are in.
fn open(b: &mut Browser, url: &str) -> Result<(), String> {
    b.go(url)?;
    let ready = format!(
        "document.readyState === 'complete' && location.href === {}",
        js_str(url)
    );
    b.wait(&ready, WAIT)
        .map(|_| ())
        .ok_or_else(|| format!("{url} did not finish loading in {} s", WAIT.as_secs()))
}

/// A JavaScript value as a whole number.
fn as_int(v: &Value, what: &str) -> Result<i128, String> {
    match v {
        Value::Int(n) => Ok(*n),
        Value::Float(f) if f.is_finite() => Ok(f.floor() as i128),
        other => Err(format!("{what}: not a number: {other:?}")),
    }
}

fn eval_int(b: &mut Browser, expr: &str) -> Result<i128, String> {
    let v = b.eval(expr)?;
    as_int(&v, expr)
}

/// The top-left `width` x `height` pixels of the page, cleaned.
fn capture(b: &mut Browser, width: i128, height: i128) -> Result<(Vec<u8>, u32, u32), String> {
    let r = b.call(
        "Page.captureScreenshot",
        obj([
            ("format", Value::Str("png".into())),
            (
                "clip",
                obj([
                    ("x", int(0)),
                    ("y", int(0)),
                    ("width", int(width)),
                    ("height", int(height)),
                    ("scale", int(1)),
                ]),
            ),
            ("captureBeyondViewport", Value::Bool(true)),
        ]),
    )?;
    let Value::Object(m) = &r else {
        return Err("Page.captureScreenshot answered without an object".into());
    };
    let Some(Value::Str(data)) = m.get("data") else {
        return Err("Page.captureScreenshot answered without data".into());
    };
    clean_png(&base64_decode(data)?)
}

/// The text of every visible text node of the page that falls, at least
/// in part, inside the page's top-left `width` x `height` pixels: the text
/// the picture of that area shows.
fn shown_text(b: &mut Browser, width: i128, height: i128) -> Result<String, String> {
    let js = format!(
        r#"((x1, y1) => {{
  const out = [];
  const skip = ['SCRIPT', 'STYLE', 'NOSCRIPT', 'TEMPLATE', 'TITLE'];
  const walk = document.createTreeWalker(document.documentElement, NodeFilter.SHOW_TEXT);
  for (let n = walk.nextNode(); n; n = walk.nextNode()) {{
    const el = n.parentElement;
    if (!n.textContent.trim() || !el || skip.includes(el.tagName.toUpperCase())) continue;
    const cs = getComputedStyle(el);
    if (cs.visibility !== 'visible' || cs.display === 'none') continue;
    const r = document.createRange();
    r.selectNodeContents(n);
    for (const b of r.getClientRects()) {{
      const top = b.top + scrollY, left = b.left + scrollX;
      if (b.width > 0 && b.height > 0 && left < x1 && left + b.width > 0 && top < y1 && top + b.height > 0) {{
        out.push(n.textContent);
        break;
      }}
    }}
  }}
  return out.join('\n');
}})({width}, {height})"#
    );
    match b.eval(&js)? {
        Value::Str(s) => Ok(s),
        Value::Null => Ok(String::new()),
        other => Err(format!("the page's text: not a string: {other:?}")),
    }
}

/// One picture taken: written to `path` and its text checked. True when
/// its text passed.
fn write_picture(
    path: &Path,
    png: &[u8],
    text: &str,
    expect: &str,
    terms: &[Term],
    note: &str,
    out: &mut dyn Write,
) -> Result<bool, String> {
    std::fs::write(path, png).map_err(|e| format!("{}: {e}", path.display()))?;
    let problems = text_problems(text, expect, terms);
    if problems.is_empty() {
        let _ = writeln!(
            out,
            "  {}{note}; its text shows {expect} and no private term",
            path.display()
        );
        Ok(true)
    } else {
        let _ = writeln!(out, "FAIL {}: {}", path.display(), problems.join("; "));
        Ok(false)
    }
}

fn file_url(p: &Path) -> Result<String, String> {
    let abs = std::fs::canonicalize(p).map_err(|e| format!("{}: {e}", p.display()))?;
    Ok(format!("file://{}", abs.display()))
}

// ---------------------------------------------------------------- take

fn take(args: &[String], out: &mut dyn Write) -> Result<u8, Fail> {
    let o = Opts::parse(
        args,
        &[
            "--genome-dir",
            "--out",
            "--only",
            "--html-stop-before",
            "--cpic-head",
            "--terms-file",
        ],
        &["--cpic-full"],
    )?;
    o.no_free()?;
    let genome = PathBuf::from(o.need("--genome-dir")?);
    let dest = PathBuf::from(o.need("--out")?);
    let only = o.get("--only");
    if let Some(x) = only.filter(|x| !matches!(*x, "html" | "cpic")) {
        return Err(Fail::Usage(format!("--only: html or cpic, not {x}")));
    }
    let head =
        match o.get("--cpic-head") {
            None => None,
            Some(n) => Some(n.parse::<usize>().ok().filter(|&n| n > 0).ok_or_else(|| {
                Fail::Usage(format!("--cpic-head: a whole number above 0, not {n}"))
            })?),
        };
    let full = o.has("--cpic-full");
    if full && head.is_some() {
        return Err(Fail::Usage(
            "give --cpic-head or --cpic-full, not both".into(),
        ));
    }
    let stops: Vec<String> = o
        .get("--html-stop-before")
        .unwrap_or(STOP_BEFORE)
        .split(',')
        .map(str::trim)
        .filter(|t| !t.is_empty())
        .map(str::to_owned)
        .collect();
    let terms = load_terms(o.get("--terms-file"))?;
    std::fs::create_dir_all(&dest).map_err(|e| format!("{}: {e}", dest.display()))?;
    let scratch = Scratch::new()?;
    let mut b = open_browser(&scratch)?;
    let mut ok = true;
    if only != Some("html") {
        ok &= cpic_picture(&mut b, &scratch, &genome, &dest, head, full, &terms, out)?;
    }
    if only != Some("cpic") {
        ok &= html_pictures(&mut b, &genome.join(SAMPLE), &dest, &stops, &terms, out)?;
    }
    Ok(if ok { 0 } else { 1 })
}

/// The card tops: the page's y of the first `.card` whose `h2` starts with
/// one of `titles`, or null.
fn cut_js(titles: &[String]) -> String {
    let list = titles
        .iter()
        .map(|t| js_str(t))
        .collect::<Vec<_>>()
        .join(", ");
    format!(
        "(() => {{ const titles = [{list}]; \
         const stop = [...document.querySelectorAll('.card')].find(c => {{ const h = c.querySelector('h2'); \
         return h && titles.some(t => h.textContent.trim().startsWith(t)); }}); \
         return stop ? Math.floor(stop.getBoundingClientRect().top + window.scrollY) : null; }})()"
    )
}

fn html_pictures(
    b: &mut Browser,
    sample_dir: &Path,
    dest: &Path,
    stops: &[String],
    terms: &[Term],
    out: &mut dyn Write,
) -> Result<bool, Fail> {
    let report = sample_dir.join(format!("{SAMPLE}_report.html"));
    let multiqc = sample_dir.join("multiqc/multiqc_report.html");
    for p in [&report, &multiqc] {
        if !p.is_file() {
            return Err(Fail::Error(format!(
                "{} not found: run steps 24 and 28 first",
                p.display()
            )));
        }
    }
    let (w, h) = (1180, 900);
    window(b, w, h, false)?;
    open(b, &file_url(&report)?)?;
    // A window as tall as the page, so the whole crop is laid out as shown.
    let tall = eval_int(b, "document.documentElement.scrollHeight")?;
    window(b, w, tall.max(h), false)?;
    let cut = match b.eval(&cut_js(stops))? {
        Value::Null => {
            return Err(Fail::Error(format!(
                "no card titled {stops:?} in {}; the crop has nothing to stop above",
                report.display()
            )));
        }
        v => as_int(&v, "the card's top")?,
    };
    // The grid gap is 20 pixels: leave 12 below the last card shown.
    let ch = cut - 8;
    let (png, _, _) = capture(b, w, ch)?;
    let text = shown_text(b, w, ch)?;
    let mut ok = write_picture(
        &dest.join("demo-html-report.png"),
        &png,
        &text,
        SAMPLE,
        terms,
        &format!(": stops above the card at y={cut}"),
        out,
    )?;

    let (w, h) = (1600, 1118);
    window(b, w, h, false)?;
    open(b, &file_url(&multiqc)?)?;
    // A browser that never opened a MultiQC report gets a "Welcome!" bar
    // over the top of it: close it, as a returning reader would see it.
    // MultiQC versions without the bar leave nothing to click.
    let _ = b.wait(
        "(() => { const el = [...document.querySelectorAll('a, button')]\
         .find(e => e.textContent.includes(\"don't show again\")); \
         if (el) { el.click(); return true; } return false; })()",
        Duration::from_secs(3),
    );
    // The plots draw after load.
    sleep(Duration::from_secs(3));
    let (png, _, _) = capture(b, w, h)?;
    let text = shown_text(b, w, h)?;
    ok &= write_picture(
        &dest.join("demo-multiqc-report.png"),
        &png,
        &text,
        SAMPLE,
        terms,
        "",
        out,
    )?;
    Ok(ok)
}

/// The number of lines up to the last row of the Gene Results table.
fn gene_table_end(lines: &[&str]) -> Result<usize, String> {
    let start = lines
        .iter()
        .position(|l| l.starts_with("Gene Results:"))
        .ok_or("the CPIC report has no 'Gene Results:' section")?;
    Ok((start + 1..lines.len())
        .find(|&i| lines[i].trim().is_empty())
        .unwrap_or(lines.len()))
}

fn escape(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
}

/// The terminal page: a 30-pixel title bar with three dots, then the
/// prompt line and `lines` in 14-pixel monospace on 20-pixel rows, wrapped
/// at the last column as a terminal wraps them.
fn terminal_html(command: &str, lines: &[&str]) -> String {
    let body = lines
        .iter()
        .map(|l| escape(l))
        .collect::<Vec<_>>()
        .join("\n");
    let dots = ["#ff5f57", "#febc2e", "#28c840"]
        .iter()
        .enumerate()
        .map(|(i, c)| format!("<i style=\"left:{}px;background:{c}\"></i>", 11 + 20 * i))
        .collect::<String>();
    format!(
        "<!doctype html><meta charset=\"utf-8\"><style>\
         html,body{{margin:0;background:#1e1e28}}\
         .bar{{height:30px;background:#2c2c38;position:relative}}\
         .bar i{{position:absolute;top:9px;width:12px;height:12px;border-radius:50%}}\
         pre{{margin:0;padding:22px;font:14px/20px Menlo,\"DejaVu Sans Mono\",monospace;color:#dedee4;\
         white-space:pre-wrap;word-break:break-all}}\
         b{{font-weight:normal;color:#50fa7b}}</style>\
         <div class=\"bar\">{dots}</div><pre><b>$ </b>{}\n{body}</pre>",
        escape(command)
    )
}

#[allow(clippy::too_many_arguments)]
fn cpic_picture(
    b: &mut Browser,
    scratch: &Scratch,
    genome: &Path,
    dest: &Path,
    head: Option<usize>,
    full: bool,
    terms: &[Term],
    out: &mut dyn Write,
) -> Result<bool, Fail> {
    let rel = format!("{SAMPLE}/cpic/{SAMPLE}_cpic_recommendations.txt");
    let path = genome.join(&rel);
    let text = std::fs::read_to_string(&path)
        .map_err(|e| Fail::Error(crate::py_os_error(&e, &path.display().to_string())))?;
    let mut lines: Vec<&str> = text.lines().collect();
    let command = if full {
        format!("cat {rel}")
    } else {
        let n = match head {
            Some(n) => n,
            None => gene_table_end(&lines)?,
        };
        lines.truncate(n);
        format!("head -n {n} {rel}")
    };
    let page = scratch.0.join("terminal.html");
    std::fs::write(&page, terminal_html(&command, &lines))
        .map_err(|e| format!("{}: {e}", page.display()))?;
    let w = 900;
    window(b, w, 200, false)?;
    open(b, &file_url(&page)?)?;
    let h = eval_int(b, "document.documentElement.scrollHeight")?;
    window(b, w, h, false)?;
    let (png, _, ph) = capture(b, w, h)?;
    let rows = (i128::from(ph) - 30 - 44) / 20;
    let shown = format!("$ {command}\n{}", lines.join("\n"));
    Ok(write_picture(
        &dest.join("demo-cpic-report.png"),
        &png,
        &shown,
        SAMPLE,
        terms,
        &format!(": {rows} terminal lines"),
        out,
    )?)
}

// ---------------------------------------------------------------- live

fn live(args: &[String], out: &mut dyn Write) -> Result<u8, Fail> {
    let o = Opts::parse(args, &["--site", "--out", "--terms-file"], &["--phone"])?;
    o.no_free()?;
    let site = PathBuf::from(o.need("--site")?);
    let dest = PathBuf::from(o.need("--out")?);
    let phone = o.has("--phone");
    let terms = load_terms(o.get("--terms-file"))?;
    let person = if phone { PARTNER } else { SAMPLE };
    if !site.join("runs").join(person).join("index.html").is_file() {
        return Err(Fail::Error(format!(
            "{} is not a built live demo (no runs/{person}/index.html): build it with tests/demo/build_live_demo.py --out {}",
            site.display(),
            site.display()
        )));
    }
    std::fs::create_dir_all(&dest).map_err(|e| format!("{}: {e}", dest.display()))?;
    let scratch = Scratch::new()?;
    // Served as a static host serves the docs site.
    let srv =
        StaticServer::start(&site, &[]).map_err(|e| format!("serve {}: {e}", site.display()))?;
    let mut b = open_browser(&scratch)?;
    let (w, mut h) = if phone {
        (PHONE_WIDTH, 844)
    } else {
        (1280, 1000)
    };
    window(&mut b, w, h, true)?;
    let url = srv.url(&format!("runs/{person}/"));
    open(&mut b, &url)?;
    b.wait(
        "document.documentElement.dataset.pgpStatus === 'ready'",
        WAIT,
    )
    .ok_or_else(|| format!("{url}: the page never said ready"))?;
    // The last moment: the run ended, the charts are drawn, the re-check is in.
    b.eval(
        "(() => { const s = document.getElementById('replay-slider'); s.value = s.max; \
         s.dispatchEvent(new Event('input')); return true; })()",
    )?;
    b.wait(
        "Number(document.documentElement.dataset.pgpCharts) > 0 \
         && Number(document.documentElement.dataset.pgpResults) > 0",
        WAIT,
    )
    .ok_or_else(|| {
        format!("{url}: the charts or the results never showed at the end of the replay")
    })?;
    sleep(Duration::from_secs(1));
    if phone {
        // A window as tall as the page, so the tab bar sits at its foot, as
        // a reader scrolled to the end sees it.
        h = eval_int(&mut b, "document.documentElement.scrollHeight")?;
        window(&mut b, w, h, true)?;
        sleep(Duration::from_millis(500));
    }
    let (png, pw, _) = capture(&mut b, w, h)?;
    let text = shown_text(&mut b, w, h)?;
    let name = if phone {
        "demo-live-phone.png"
    } else {
        "demo-live-site.png"
    };
    let mut ok = write_picture(&dest.join(name), &png, &text, person, &terms, "", out)?;
    if phone && i128::from(pw) != PHONE_WIDTH {
        let _ = writeln!(out, "FAIL {name}: {pw} px wide, want {PHONE_WIDTH}");
        ok = false;
    }
    Ok(if ok { 0 } else { 1 })
}

// ---------------------------------------------------------------- check

fn check(args: &[String], out: &mut dyn Write) -> Result<u8, Fail> {
    let o = Opts::parse(args, &["--terms-file"], &["--self-test"])?;
    if o.has("--self-test") {
        o.no_free()?;
        return self_test(o.get("--terms-file"), out);
    }
    if o.get("--terms-file").is_some() {
        return Err(Fail::Usage(
            "--terms-file goes with --self-test, take or live: check reads the chunks only".into(),
        ));
    }
    if o.free.is_empty() {
        return Err(Fail::Usage("give the PNGs to check".into()));
    }
    let mut failed = 0;
    for p in &o.free {
        let data = std::fs::read(p).map_err(|e| Fail::Error(crate::py_os_error(&e, p)))?;
        match chunk_problems(&data) {
            Ok(problems) if problems.is_empty() => {
                let (w, h) = png_size(&data);
                let _ = writeln!(out, "PASS {p}: {w} x {h}, image chunks only");
            }
            Ok(problems) => {
                failed += 1;
                let _ = writeln!(out, "FAIL {p}: {}", problems.join("; "));
            }
            Err(e) => {
                failed += 1;
                let _ = writeln!(out, "FAIL {p}: {e}");
            }
        }
    }
    Ok(u8::from(failed > 0))
}

/// A small white PNG, with a tEXt chunk after its IHDR when `text` is set.
fn planted_png(text: Option<&str>) -> Vec<u8> {
    let mut data = Vec::new();
    {
        let mut e = png::Encoder::new(&mut data, 8, 4);
        e.set_color(png::ColorType::Rgb);
        e.set_depth(png::BitDepth::Eight);
        let mut w = e.write_header().expect("a PNG header");
        w.write_image_data(&[255; 8 * 4 * 3]).expect("PNG pixels");
        w.finish().expect("a PNG");
    }
    if let Some(t) = text {
        let mut body = b"Comment\0".to_vec();
        body.extend_from_slice(t.as_bytes());
        let mut kind_body = b"tEXt".to_vec();
        kind_body.extend_from_slice(&body);
        let mut chunk = (body.len() as u32).to_be_bytes().to_vec();
        chunk.extend_from_slice(&kind_body);
        chunk.extend_from_slice(&crc32(&kind_body).to_be_bytes());
        // After the signature (8) and IHDR (25).
        data.splice(33..33, chunk);
    }
    data
}

/// CRC-32 (ISO 3309), the checksum of a PNG chunk.
fn crc32(bytes: &[u8]) -> u32 {
    let mut c = !0u32;
    for &b in bytes {
        c ^= u32::from(b);
        for _ in 0..8 {
            c = if c & 1 == 1 {
                0xedb8_8320 ^ (c >> 1)
            } else {
                c >> 1
            };
        }
    }
    !c
}

/// Every check against a planted case: 0 when each behaves as it must.
fn self_test(terms_file: Option<&str>, out: &mut dyn Write) -> Result<u8, Fail> {
    let terms = load_terms(terms_file)?;
    let chunks_ok = |d: &[u8]| chunk_problems(d).is_ok_and(|p| p.is_empty());
    let text_ok = |t: &str| text_problems(t, SAMPLE, &terms).is_empty();
    let clean = planted_png(None);
    let marked = planted_png(Some("written on a laptop"));
    let cleaned = clean_png(&marked).map(|(d, _, _)| d).unwrap_or_default();
    // Built here, so this file holds no path, address or host name the
    // personal-data check would flag.
    let home = format!("/{}", ["Users", "someone", "genome"].join("/"));
    let addr = ["192", "168", "20", "31"].join(".");
    let host = ["box", "example", "ts", "net"].join(".");
    let mut cases: Vec<(String, bool, bool)> = vec![
        ("a clean picture passes".into(), chunks_ok(&clean), true),
        ("a tEXt chunk is caught".into(), chunks_ok(&marked), false),
        (
            "a cleaned picture has lost its tEXt chunk".into(),
            chunks_ok(&cleaned),
            true,
        ),
        (
            "a clean text passes".into(),
            text_ok("Sample: DEMO-001\nGenomic Analysis Report"),
            true,
        ),
        (
            "a home path is caught".into(),
            text_ok(&format!("Sample: DEMO-001\nOutput: {home}")),
            false,
        ),
        (
            "a private address is caught".into(),
            text_ok(&format!("Sample: DEMO-001\nhost {addr}")),
            false,
        ),
        (
            "a tailnet name is caught".into(),
            text_ok(&format!("Sample: DEMO-001\nhost {host}")),
            false,
        ),
        (
            "a text without DEMO-001 is caught".into(),
            text_ok("Sample: OTHER-002"),
            false,
        ),
    ];
    if terms_file.is_some()
        && let Some(first) = terms.iter().find(|t| t.kind == "private term")
    {
        let word: String = first
            .rx
            .as_str()
            .trim_start_matches("(?i)")
            .chars()
            .filter(|c| c.is_ascii_alphanumeric() || *c == ' ')
            .collect();
        let word = if word.is_empty() {
            "x".to_owned()
        } else {
            word
        };
        cases.push((
            "the first line of the terms file is caught".into(),
            text_ok(&format!("Sample: DEMO-001\nRun by {word}")),
            false,
        ));
    }
    let mut bad = 0;
    for (desc, passed, should_pass) in &cases {
        let right = passed == should_pass;
        bad += usize::from(!right);
        let _ = writeln!(
            out,
            "[{}] self-test: {desc}",
            if right { "PASS" } else { "FAIL" }
        );
    }
    if bad > 0 {
        let _ = writeln!(
            out,
            "self-test: {bad} planted case(s) were not caught as expected"
        );
        return Ok(1);
    }
    let _ = writeln!(out, "self-test: every planted case behaved as expected");
    Ok(0)
}

// ---------------------------------------------------------------- same-sample

static TOTAL: LazyLock<regex::Regex> = LazyLock::new(|| {
    regex::Regex::new(
        r#"<span class="label">Total variants</span><span class="value">([^<]*)</span>"#,
    )
    .expect("the Total variants pattern compiles")
});

/// The Total variants value of a step 24 HTML report.
fn total_variants(path: &str) -> Result<String, Fail> {
    let html =
        std::fs::read_to_string(path).map_err(|e| Fail::Error(crate::py_os_error(&e, path)))?;
    TOTAL
        .captures(&html)
        .map(|c| c[1].trim().to_owned())
        .ok_or_else(|| {
            Fail::Error(format!(
                "{path} has no Total variants line: is it a step 24 HTML report?"
            ))
        })
}

fn same_sample(args: &[String], out: &mut dyn Write) -> Result<u8, Fail> {
    let o = Opts::parse(args, &[], &[])?;
    let [report, live] = o.free.as_slice() else {
        return Err(Fail::Usage("give REPORT and LIVE_REPORT".into()));
    };
    let (a, b) = (total_variants(report)?, total_variants(live)?);
    if a == b {
        let _ = writeln!(out, "PASS: both reports say Total variants {a}");
        Ok(0)
    } else {
        let _ = writeln!(
            out,
            "FAIL: the pictured report {report} says Total variants {a}, the live demo's {live} says {b}: \
             the docs call them the same sample, so write both with one record count"
        );
        Ok(1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run(args: &[&str]) -> (u8, String, String) {
        let args: Vec<String> = args.iter().map(|s| (*s).to_owned()).collect();
        let (mut out, mut err) = (Vec::new(), Vec::new());
        let code = cli(&args, &mut out, &mut err);
        (
            code,
            String::from_utf8(out).unwrap(),
            String::from_utf8(err).unwrap(),
        )
    }

    fn tmp(name: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!(
            "og-dev-pictures-test-{name}-{}",
            std::process::id()
        ));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    #[test]
    fn base64_reads_what_chrome_sends() {
        assert_eq!(base64_decode("aGVsbG8=").unwrap(), b"hello");
        assert_eq!(base64_decode("aGVsbG8hIQ==").unwrap(), b"hello!!");
        assert_eq!(base64_decode("").unwrap(), b"");
        assert!(base64_decode("a b").is_err());
    }

    #[test]
    fn crc_matches_the_png_reference() {
        // The IEND chunk's CRC every PNG ends with.
        assert_eq!(crc32(b"IEND"), 0xae42_6082);
    }

    #[test]
    fn a_text_chunk_fails_and_cleaning_drops_it() {
        let marked = planted_png(Some("x"));
        assert_eq!(chunk_problems(&marked).unwrap(), vec!["chunk tEXt"]);
        let (clean, w, h) = clean_png(&marked).unwrap();
        assert_eq!((w, h), (8, 4));
        assert!(chunk_problems(&clean).unwrap().is_empty());
        let kinds = chunk_kinds(&clean).unwrap();
        assert_eq!(kinds.first(), Some(b"IHDR"));
        assert_eq!(kinds.last(), Some(b"IEND"));
        assert!(chunk_kinds(b"GIF89a").is_err());
        assert!(chunk_kinds(&marked[..marked.len() - 3]).is_err());
    }

    #[test]
    fn check_names_the_bad_chunk_and_the_size() {
        let d = tmp("check");
        let good = d.join("good.png");
        let bad = d.join("bad.png");
        std::fs::write(&good, planted_png(None)).unwrap();
        std::fs::write(&bad, planted_png(Some("x"))).unwrap();
        let (code, out, _) = run(&["check", good.to_str().unwrap()]);
        assert_eq!(code, 0, "{out}");
        assert!(out.contains("8 x 4, image chunks only"), "{out}");
        let (code, out, _) = run(&["check", good.to_str().unwrap(), bad.to_str().unwrap()]);
        assert_eq!(code, 1, "{out}");
        assert!(out.contains("chunk tEXt"), "{out}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn the_text_check_finds_each_kind() {
        let terms = load_terms(None).unwrap();
        assert!(text_problems("Sample: DEMO-001", SAMPLE, &terms).is_empty());
        // A hyphen the page breaks or drops still counts.
        assert!(text_problems("Sample: DEMO001", SAMPLE, &terms).is_empty());
        assert_eq!(text_problems("Sample: X", SAMPLE, &terms).len(), 1);
        let home = format!("DEMO-001 /{}/a/b", "home");
        assert!(text_problems(&home, SAMPLE, &terms)[0].contains("home or volume path"));
        let ip = format!("DEMO-001 {}", ["10", "1", "2", "3"].join("."));
        assert!(text_problems(&ip, SAMPLE, &terms)[0].contains("private address"));
        // A public address, and a version string, are not private.
        assert!(text_problems("DEMO-001 8.8.8.8 v110.0.1.2", SAMPLE, &terms).is_empty());
    }

    #[test]
    fn a_terms_file_adds_its_lines_and_an_empty_one_is_refused() {
        let d = tmp("terms");
        let f = d.join("terms.txt");
        std::fs::write(&f, "# mine\nsecretname\n").unwrap();
        let terms = load_terms(Some(f.to_str().unwrap())).unwrap();
        assert_eq!(terms.len(), BUILT_IN.len() + 1);
        let p = text_problems("DEMO-001 by SecretName", SAMPLE, &terms);
        assert_eq!(
            p,
            vec!["its text holds a private term (pattern 4 of the list)"]
        );
        std::fs::write(&f, "# only a comment\n\n").unwrap();
        assert!(load_terms(Some(f.to_str().unwrap())).is_err());
        std::fs::write(&f, "a(b\n").unwrap();
        let Err(Fail::Error(e)) = load_terms(Some(f.to_str().unwrap())) else {
            panic!("a bad pattern must be refused");
        };
        assert!(e.contains("line 1") && !e.contains("a(b"), "{e}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn the_self_test_passes_and_uses_the_terms_file() {
        let (code, out, _) = run(&["check", "--self-test"]);
        assert_eq!(code, 0, "{out}");
        assert_eq!(out.matches("[PASS]").count(), 8, "{out}");
        let d = tmp("selftest");
        let f = d.join("terms.txt");
        std::fs::write(&f, "secretname\n").unwrap();
        let (code, out, _) = run(&["check", "--self-test", "--terms-file", f.to_str().unwrap()]);
        assert_eq!(code, 0, "{out}");
        assert!(
            out.contains("[PASS] self-test: the first line of the terms file is caught"),
            "{out}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn the_gene_table_ends_at_its_first_blank_line() {
        let lines = [
            "CPIC",
            "",
            "Gene Results:",
            "  CYP2C19 *1/*2",
            "  CYP2D6 *1/*1",
            "",
            "Drugs",
        ];
        assert_eq!(gene_table_end(&lines).unwrap(), 5);
        assert_eq!(gene_table_end(&lines[..5]).unwrap(), 5);
        assert!(gene_table_end(&["no table"]).is_err());
    }

    #[test]
    fn the_terminal_page_escapes_the_report() {
        let page = terminal_html("head -n 2 x", &["a <b> & c"]);
        assert!(page.contains("a &lt;b&gt; &amp; c"));
        assert!(page.contains("<b>$ </b>head -n 2 x\n"));
    }

    fn report(total: &str) -> String {
        format!(
            "<div class=\"card\"><h2>Variant Calling</h2>\n    <div class=\"stat\"><span class=\"label\">Total variants</span><span class=\"value\">{total}</span></div></div>"
        )
    }

    #[test]
    fn same_sample_fails_on_two_record_counts() {
        let d = tmp("same");
        let (a, b, c) = (d.join("a.html"), d.join("b.html"), d.join("c.html"));
        // The pair the docs had: the picture from 4,700,000 records, the
        // live demo from 20,000.
        std::fs::write(&a, report("4700004")).unwrap();
        std::fs::write(&b, report("20000")).unwrap();
        std::fs::write(&c, "<p>no card</p>").unwrap();
        let (code, out, _) = run(&["same-sample", a.to_str().unwrap(), b.to_str().unwrap()]);
        assert_eq!(code, 1, "{out}");
        assert!(
            out.contains("says Total variants 4700004") && out.contains("says 20000"),
            "{out}"
        );
        std::fs::write(&a, report("20000")).unwrap();
        let (code, out, _) = run(&["same-sample", a.to_str().unwrap(), b.to_str().unwrap()]);
        assert_eq!(code, 0, "{out}");
        let (code, _, err) = run(&["same-sample", a.to_str().unwrap(), c.to_str().unwrap()]);
        assert_eq!(code, 2);
        assert!(err.contains("no Total variants line"), "{err}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn usage_errors_exit_2() {
        assert_eq!(run(&[]).0, 2);
        assert_eq!(run(&["take", "--out", "x"]).0, 2);
        assert_eq!(
            run(&["take", "--genome-dir", "g", "--out", "o", "--only", "live"]).0,
            2
        );
        assert_eq!(
            run(&[
                "take",
                "--genome-dir",
                "g",
                "--out",
                "o",
                "--cpic-head",
                "3",
                "--cpic-full"
            ])
            .0,
            2
        );
        assert_eq!(run(&["check"]).0, 2);
        assert_eq!(run(&["same-sample", "a"]).0, 2);
        assert_eq!(run(&["live", "--site", "s", "--out", "o", "--bogus"]).0, 2);
    }
}
