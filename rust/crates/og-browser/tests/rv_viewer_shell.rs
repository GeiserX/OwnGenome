//! The genome viewer's page (assets/viewer: index.html and viewer.js) in a headless Chrome, on a synthetic
//! sample: a random 4 kb reference with its .fai and a VCF with one variant. The page, its script, the pinned
//! igv.js release, config.json and the sample's files are served by a stub that answers the paths `og view
//! serve` answers, with that server's Content-Security-Policy, and can make config.json answer an error.
//! og-browser does not depend on og-live, so the stub keeps a copy of the policy, and
//! `the_policy_is_og_view_serves` keeps the two equal.
//!
//! - At 1440x900 and at 390x844 the page is dark: the root's colour scheme is dark and the body has an opaque
//!   dark background (main's page had none, so the browser's white showed under the header).
//! - The igv.js search input and variant track label have readable text on their backgrounds.
//! - At the live map's short iframe height, the variant track fits without vertical scrolling.
//! - On the 390 px phone, with igv.js started, the page is no wider than the screen: igv.js's own bar
//!   scrolls inside its panel (on main it widened the page past 390 px).
//! - The page has one h1, which names the sample, one main landmark holding the error box and igv.js, and a
//!   line that says search takes coordinates, not gene names.
//! - When config.json answers 500, #error has role alert, names the failure, says what to check and offers
//!   Try again; once the server answers again, Try again starts the viewer.
//!
//! The PR records the failing baseline for each check. Without Chrome the test says SKIP, and fails
//! in CI. Set `OG_BROWSER_SHOTS` to a folder to save a screenshot of each view there.

use og_browser::{Browser, Scratch, Value, chrome_or_skip, js_str, repo_root};
use pycompat::json;
use std::collections::HashMap;
use std::fs;
use std::io::{BufRead, BufReader, Write};
use std::net::{Shutdown, TcpListener, TcpStream};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex, OnceLock};
use std::thread;
use std::time::Duration;

/// `CSP` of rust/crates/og-live/src/viewer.rs, the policy `og view serve` sends.
const CSP: &str = "default-src 'none'; script-src 'self'; connect-src 'self'; img-src 'self' data: blob:; \
     style-src 'self' 'unsafe-inline'; font-src 'self' data:; worker-src 'self' blob:; \
     base-uri 'none'; form-action 'none'; frame-ancestors 'none'";
const PHONE: (u32, u32) = (390, 844);
const DESKTOP: (u32, u32) = (1440, 900);
const STATUS: &str = "document.documentElement.dataset.pgpStatus || ''";
const SETTLED: &str = "(function(){var s=document.documentElement.dataset.pgpStatus||'';return s==='ready'||s.indexOf('error')===0;})()";

/// How the stub answers config.json.
#[derive(Clone)]
enum Config {
    Json(String),
    Status(u16),
}

/// The paths `og view serve` answers, on a free port of 127.0.0.1: `/` and `/index.html`, `/viewer.js`,
/// `/igv.min.js`, `/config.json` and `/data/<i>/<name>` (byte ranges too).
struct Stub {
    port: u16,
    config: Arc<Mutex<Config>>,
    stop: Arc<AtomicBool>,
}

impl Stub {
    fn start(files: HashMap<String, PathBuf>, config: Config) -> Stub {
        let listener = TcpListener::bind("127.0.0.1:0").unwrap();
        let port = listener.local_addr().unwrap().port();
        let config = Arc::new(Mutex::new(config));
        let stop = Arc::new(AtomicBool::new(false));
        let (files, c, s) = (Arc::new(files), Arc::clone(&config), Arc::clone(&stop));
        thread::spawn(move || {
            for conn in listener.incoming() {
                if s.load(Ordering::SeqCst) {
                    break;
                }
                if let Ok(conn) = conn {
                    let (files, c) = (Arc::clone(&files), Arc::clone(&c));
                    thread::spawn(move || answer(conn, &files, &c));
                }
            }
        });
        Stub { port, config, stop }
    }

    fn set_config(&self, c: Config) {
        *self.config.lock().unwrap() = c;
    }

    fn url(&self, path: &str) -> String {
        format!("http://127.0.0.1:{}{path}", self.port)
    }
}

impl Drop for Stub {
    fn drop(&mut self) {
        self.stop.store(true, Ordering::SeqCst);
        let _ = TcpStream::connect(("127.0.0.1", self.port));
    }
}

fn answer(mut c: TcpStream, files: &HashMap<String, PathBuf>, config: &Mutex<Config>) {
    let _ = c.set_read_timeout(Some(Duration::from_secs(10)));
    let mut r = BufReader::new(&c);
    let mut first = String::new();
    if r.read_line(&mut first).is_err() {
        return;
    }
    let (mut line, mut range) = (String::new(), None);
    while r.read_line(&mut line).is_ok_and(|n| n > 2) {
        if let Some((k, v)) = line.split_once(':')
            && k.eq_ignore_ascii_case("range")
        {
            range = Some(v.trim().to_owned());
        }
        line.clear();
    }
    let target = first.split_whitespace().nth(1).unwrap_or("/");
    let path = target.split(['?', '#']).next().unwrap_or("/");
    let ctype = |p: &Path| match p.extension().and_then(|e| e.to_str()) {
        Some("html") => "text/html; charset=utf-8",
        Some("js") => "text/javascript; charset=utf-8",
        _ => "application/octet-stream",
    };
    let (status, ctype, body, extra) = if path == "/config.json" {
        match config.lock().unwrap().clone() {
            Config::Json(j) => (
                "200 OK".to_owned(),
                "application/json",
                j.into_bytes(),
                String::new(),
            ),
            Config::Status(n) => (
                format!("{n} Stub"),
                "text/plain",
                b"stub error\n".to_vec(),
                String::new(),
            ),
        }
    } else if let Some(body) = files.get(path).and_then(|f| fs::read(f).ok()) {
        let size = body.len();
        match range.as_deref().and_then(|r| byte_range(r, size)) {
            Some((a, b)) => (
                "206 Partial Content".to_owned(),
                ctype(&files[path]),
                body[a..=b].to_vec(),
                format!("Content-Range: bytes {a}-{b}/{size}\r\n"),
            ),
            None => (
                "200 OK".to_owned(),
                ctype(&files[path]),
                body,
                String::new(),
            ),
        }
    } else {
        (
            "404 Not Found".to_owned(),
            "text/plain",
            b"not found\n".to_vec(),
            String::new(),
        )
    };
    let head = format!(
        "HTTP/1.1 {status}\r\nContent-Type: {ctype}\r\nContent-Length: {}\r\nCache-Control: no-store\r\n\
         Accept-Ranges: bytes\r\nContent-Security-Policy: {CSP}\r\n{extra}Connection: close\r\n\r\n",
        body.len()
    );
    let _ = c.write_all(head.as_bytes());
    let _ = c.write_all(&body);
    let _ = c.shutdown(Shutdown::Both);
}

/// The first and last byte of `bytes=a-b` (or `a-`) within `size`, else `None` (the whole file).
fn byte_range(header: &str, size: usize) -> Option<(usize, usize)> {
    let (a, b) = header.strip_prefix("bytes=")?.split_once('-')?;
    let a: usize = a.trim().parse().ok()?;
    let b = match b.trim() {
        "" => size.checked_sub(1)?,
        b => b.parse::<usize>().ok()?.min(size.checked_sub(1)?),
    };
    (a <= b).then_some((a, b))
}

/// A random sequence of `n` bases from `seed`.
fn sequence(seed: &mut u64, n: usize) -> String {
    (0..n)
        .map(|_| {
            *seed = seed
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            b"ACGT"[(*seed >> 62) as usize] as char
        })
        .collect()
}

/// The sample in `dir` (reference/synthetic.fasta with its .fai, S.vcf) and the stub's file routes.
fn make_sample(dir: &Path) -> HashMap<String, PathBuf> {
    let mut seed = 7;
    let contigs = [
        ("chr1", sequence(&mut seed, 3000)),
        ("chr2", sequence(&mut seed, 1000)),
    ];
    let (mut fasta, mut fai, mut off) = (String::new(), String::new(), 0);
    for (n, s) in &contigs {
        off += n.len() + 2;
        fai.push_str(&format!("{n}\t{}\t{off}\t60\t61\n", s.len()));
        let lines: Vec<&str> = s
            .as_bytes()
            .chunks(60)
            .map(|c| std::str::from_utf8(c).unwrap())
            .collect();
        fasta.push_str(&format!(">{n}\n{}\n", lines.join("\n")));
        off += s.len() + lines.len();
    }
    let r = &contigs[0].1[450..451];
    let alt = "ACGT".replace(r, "")[..1].to_owned();
    let vcf = format!(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=3000>\n##contig=<ID=chr2,length=1000>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\nchr1\t451\t.\t{r}\t{alt}\t50\tPASS\t.\tGT\t0/1\n"
    );
    fs::write(dir.join("synthetic.fasta"), fasta).unwrap();
    fs::write(dir.join("synthetic.fasta.fai"), fai).unwrap();
    fs::write(dir.join("S.vcf"), vcf).unwrap();
    let assets = repo_root().join("assets/viewer");
    HashMap::from([
        ("/".to_owned(), assets.join("index.html")),
        ("/index.html".to_owned(), assets.join("index.html")),
        ("/viewer.js".to_owned(), assets.join("viewer.js")),
        ("/igv.min.js".to_owned(), pinned_igv().clone()),
        (
            "/data/0/synthetic.fasta".to_owned(),
            dir.join("synthetic.fasta"),
        ),
        (
            "/data/1/synthetic.fasta.fai".to_owned(),
            dir.join("synthetic.fasta.fai"),
        ),
        ("/data/2/S.vcf".to_owned(), dir.join("S.vcf")),
    ])
}

/// config.json as `og view serve` writes it for this sample (a VCF without an index is read whole).
const CONFIG: &str = r#"{"sample": "S", "locus": "chr1:451",
 "reference": {"id": "reference", "name": "synthetic.fasta", "fastaURL": "data/0/synthetic.fasta",
  "indexURL": "data/1/synthetic.fasta.fai"},
 "tracks": [{"name": "S small variants", "type": "variant", "format": "vcf", "url": "data/2/S.vcf",
  "indexed": false}]}"#;

/// The igv.js release versions.env pins, downloaded once per test run into the target's tmp folder and
/// checked by sha256.
fn pinned_igv() -> &'static PathBuf {
    static IGV: OnceLock<PathBuf> = OnceLock::new();
    IGV.get_or_init(|| {
        let env = fs::read_to_string(repo_root().join("versions.env")).unwrap();
        let pin = |k: &str| {
            env.lines()
                .find_map(|l| l.strip_prefix(&format!("{k}=\"")))
                .and_then(|v| v.strip_suffix('"'))
                .unwrap_or_else(|| panic!("versions.env sets no {k}"))
                .to_owned()
        };
        let (version, sha) = (pin("IGV_JS_VERSION"), pin("IGV_JS_SHA256"));
        let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
        let path = dir.join(format!("igv-{version}.min.js"));
        let sum = |p: &Path| {
            let out = Command::new("shasum")
                .args(["-a", "256"])
                .arg(p)
                .output()
                .or_else(|_| Command::new("sha256sum").arg(p).output())
                .expect("shasum or sha256sum");
            String::from_utf8_lossy(&out.stdout)
                .split_whitespace()
                .next()
                .unwrap_or_default()
                .to_owned()
        };
        if !path.is_file() || sum(&path) != sha {
            let part = dir.join(format!("igv-{version}.min.js.{}", std::process::id()));
            let url = format!("https://cdn.jsdelivr.net/npm/igv@{version}/dist/igv.min.js");
            let ok = Command::new("curl")
                .args(["-fsSL", "--retry", "3", "-o"])
                .arg(&part)
                .arg(&url)
                .status()
                .is_ok_and(|s| s.success());
            assert!(ok, "cannot download {url}");
            let got = sum(&part);
            if got != sha {
                let _ = fs::remove_file(&part);
                panic!("igv.min.js has sha256 {got}, versions.env pins {sha}");
            }
            fs::rename(&part, &path).unwrap();
        }
        path
    })
}

/// A Chrome whose page is `w` x `h` CSS pixels; a phone (`mobile`) honours the page's viewport tag.
fn browser(chrome: &Path, profile: &Path, (w, h): (u32, u32), mobile: bool) -> Browser {
    let mut b = Browser::launch(chrome, profile).unwrap();
    let params =
        format!(r#"{{"width": {w}, "height": {h}, "deviceScaleFactor": 1, "mobile": {mobile}}}"#);
    b.call(
        "Emulation.setDeviceMetricsOverride",
        json::loads(&params).unwrap(),
    )
    .unwrap();
    b.call(
        "Emulation.setEmulatedMedia",
        json::loads(r#"{"features":[{"name":"prefers-color-scheme","value":"dark"}]}"#).unwrap(),
    )
    .unwrap();
    b
}

fn text(v: Result<Value, String>) -> String {
    match v {
        Ok(Value::Str(s)) => s,
        Ok(v) => format!("{v:?}"),
        Err(e) => format!("error: {e}"),
    }
}

/// The value at `key` of the JSON object the expression `expr` returns as a string.
fn probe(b: &mut Browser, expr: &str) -> Value {
    let raw = text(b.eval(expr));
    json::loads(&raw).unwrap_or_else(|e| panic!("{expr}: {raw}: {e:?}"))
}

fn get<'a>(v: &'a Value, key: &str) -> &'a Value {
    match v {
        Value::Object(m) => m.get(key).unwrap_or(&Value::Null),
        _ => &Value::Null,
    }
}

/// The colour scheme of the root, and whether `sel`'s background is opaque and dark (luminance under 60 of
/// 255), as a JSON string.
fn dark(sel: &str) -> String {
    format!(
        "(function(){{var c=getComputedStyle(document.querySelector({})).backgroundColor;\
         var m=(c.match(/[0-9.]+/g)||[]).map(Number);var a=m.length>3?m[3]:1;\
         return JSON.stringify({{color:c,scheme:getComputedStyle(document.documentElement).colorScheme,\
         dark:m.length>=3&&a>=0.99&&0.2126*m[0]+0.7152*m[1]+0.0722*m[2]<60}});}})()",
        js_str(sel)
    )
}

/// Saves a screenshot as `<OG_BROWSER_SHOTS>/viewer-<name>.png` when that variable names a folder.
fn shot(b: &mut Browser, name: &str) {
    let Some(dir) = std::env::var_os("OG_BROWSER_SHOTS").filter(|d| !d.is_empty()) else {
        return;
    };
    let r = b
        .call(
            "Page.captureScreenshot",
            json::loads(r#"{"format": "png"}"#).unwrap(),
        )
        .unwrap();
    let Value::Str(data) = get(&r, "data") else {
        panic!("no screenshot data")
    };
    fs::write(
        Path::new(&dir).join(format!("viewer-{name}.png")),
        base64(data),
    )
    .unwrap();
}

fn base64(text: &str) -> Vec<u8> {
    let val = |c: u8| match c {
        b'A'..=b'Z' => c - b'A',
        b'a'..=b'z' => c - b'a' + 26,
        b'0'..=b'9' => c - b'0' + 52,
        b'+' => 62,
        _ => 63,
    };
    let bytes: Vec<u8> = text.bytes().filter(|c| *c != b'=').map(val).collect();
    bytes
        .chunks(4)
        .flat_map(|q| {
            let n = q
                .iter()
                .enumerate()
                .fold(0u32, |n, (i, v)| n | u32::from(*v) << (18 - 6 * i));
            n.to_be_bytes()[1..q.len()].to_vec()
        })
        .collect()
}

/// Opens the viewer at `size` on a fresh stub and waits until it reports ready or an error.
fn open(
    chrome: &Path,
    scratch: &Scratch,
    name: &str,
    size: (u32, u32),
    mobile: bool,
) -> (Stub, Browser) {
    let site = scratch.path().join(format!("site-{name}"));
    fs::create_dir_all(&site).unwrap();
    let stub = Stub::start(make_sample(&site), Config::Json(CONFIG.to_owned()));
    let mut b = browser(
        chrome,
        &scratch.path().join(format!("chrome-{name}")),
        size,
        mobile,
    );
    b.go(&stub.url("/#chr1:400-500")).unwrap();
    b.wait(SETTLED, Duration::from_secs(60));
    (stub, b)
}

#[test]
fn the_policy_is_og_view_serves() {
    let src = fs::read_to_string(repo_root().join("rust/crates/og-live/src/viewer.rs")).unwrap();
    let start = src
        .find("pub const CSP: &str = \"")
        .expect("og-live's viewer.rs sets CSP")
        + 23;
    let raw = &src[start..start + src[start..].find("\";").unwrap()];
    // A Rust string continuation: a backslash, the line break and the next line's leading blanks are dropped.
    let joined: String = raw
        .split("\\\n")
        .enumerate()
        .map(|(i, p)| if i == 0 { p } else { p.trim_start() })
        .collect();
    assert_eq!(joined, CSP, "the stub must send og view serve's policy");
}

#[test]
fn the_page_is_dark_with_one_heading_and_a_main_landmark() {
    let Some(chrome) = chrome_or_skip("load the viewer page") else {
        return;
    };
    let scratch = Scratch::new(Path::new(env!("CARGO_TARGET_TMPDIR")), "rv-viewer-desktop");
    let (_stub, mut b) = open(&chrome, &scratch, "desktop", DESKTOP, false);
    assert_eq!(
        text(b.eval(STATUS)),
        "ready",
        "igv.js did not start on the synthetic sample"
    );
    assert_eq!(
        text(b.eval("document.documentElement.dataset.pgpBlocked || ''")),
        "",
        "the policy blocked what the page needs"
    );
    let d = probe(&mut b, &dark("body"));
    assert!(
        get(&d, "dark") == &Value::Bool(true) && get(&d, "scheme") == &Value::Str("dark".into()),
        "desktop: the page must be dark: {d:?}"
    );
    let outline = probe(
        &mut b,
        "JSON.stringify({h1:document.querySelectorAll('h1').length,main:document.querySelectorAll('main').length,\
         title:(document.querySelector('h1')||{}).textContent||'',\
         inside:!!document.querySelector('main #igv')&&!!document.querySelector('main #error'),\
         about:(document.querySelector('main')||document.body).textContent.replace(/\\s+/g,' ')})",
    );
    assert!(
        get(&outline, "h1") == &Value::Int(1) && get(&outline, "main") == &Value::Int(1),
        "one h1 and one main: {outline:?}"
    );
    assert_eq!(
        get(&outline, "title"),
        &Value::Str("Reads of S".into()),
        "the h1 names the sample"
    );
    assert_eq!(
        get(&outline, "inside"),
        &Value::Bool(true),
        "main holds the error box and igv.js"
    );
    let about = text(Ok(get(&outline, "about").clone()));
    assert!(
        about.contains("coordinates")
            && about.contains("gene names do not work")
            && about.contains("a few kb"),
        "the page says what search takes and when reads show: {about}"
    );
    shot(&mut b, "desktop");
}

#[test]
fn igv_search_and_track_label_have_readable_text() {
    let Some(chrome) = chrome_or_skip("check the igv.js text contrast") else {
        return;
    };
    let scratch = Scratch::new(Path::new(env!("CARGO_TARGET_TMPDIR")), "rv-viewer-contrast");
    let (_stub, mut b) = open(&chrome, &scratch, "contrast", DESKTOP, false);
    assert_eq!(text(b.eval(STATUS)), "ready");

    let contrast = probe(
        &mut b,
        r#"JSON.stringify((() => {
            const igv = document.querySelector('#igv').shadowRoot;
            const search = igv.querySelector('.igv-search-input');
            const track = Array.from(igv.querySelectorAll('.igv-track-label'))
                .find(label => label.textContent === 'S small variants');
            // The pinned igv.js styles both surfaces with opaque RGB backgrounds.
            const luminance = rgb => {
                const channels = rgb.match(/[\d.]+/g).slice(0, 3).map(channel => {
                    const value = Number(channel) / 255;
                    return value <= 0.04045 ? value / 12.92 : ((value + 0.055) / 1.055) ** 2.4;
                });
                return 0.2126 * channels[0] + 0.7152 * channels[1] + 0.0722 * channels[2];
            };
            const colors = element => {
                const { color, backgroundColor } = getComputedStyle(element);
                const fg = luminance(color);
                const bg = luminance(backgroundColor);
                return {
                    color, backgroundColor,
                    readable: (Math.max(fg, bg) + 0.05) / (Math.min(fg, bg) + 0.05) >= 4.5
                };
            };
            return { search: colors(search), track: colors(track) };
        })())"#,
    );
    assert!(
        ["search", "track"]
            .iter()
            .all(|key| get(get(&contrast, key), "readable") == &Value::Bool(true)),
        "the search input and variant label must each have 4.5:1 text contrast: {contrast:?}"
    );
}

#[test]
fn the_tracks_fit_at_the_live_map_frame_height() {
    let Some(chrome) = chrome_or_skip("check the viewer at the map iframe height") else {
        return;
    };
    let scratch = Scratch::new(Path::new(env!("CARGO_TARGET_TMPDIR")), "rv-viewer-short");
    let layouts = [600, 900].map(|width| {
        let (_stub, mut b) = open(&chrome, &scratch, &width.to_string(), (width, 260), false);
        assert_eq!(text(b.eval(STATUS)), "ready");
        probe(
            &mut b,
            r#"JSON.stringify((() => {
                const igv = document.querySelector('#igv');
                const track = igv.shadowRoot.querySelector('.igv-viewport[data-track-type=variant]');
                const rect = track.getBoundingClientRect();
                return {
                    width: innerWidth, height: innerHeight,
                    panelTop: Math.ceil(igv.getBoundingClientRect().top),
                    trackBottom: Math.ceil(rect.bottom),
                    fits: rect.height > 0 && rect.bottom <= innerHeight
                };
            })())"#,
        )
    });
    assert!(
        layouts.iter().all(|layout| {
            matches!(get(layout, "panelTop"), Value::Int(0..=60))
                && get(layout, "fits") == &Value::Bool(true)
        }),
        "at the map iframe height, igv.js must start within 60 px and show the variant track: {layouts:?}"
    );
}

#[test]
fn on_a_phone_the_page_is_dark_and_no_wider_than_the_screen() {
    let Some(chrome) = chrome_or_skip("load the viewer page on a phone") else {
        return;
    };
    let scratch = Scratch::new(Path::new(env!("CARGO_TARGET_TMPDIR")), "rv-viewer-phone");
    let (_stub, mut b) = open(&chrome, &scratch, "phone", PHONE, true);
    assert_eq!(
        text(b.eval(STATUS)),
        "ready",
        "igv.js did not start on the synthetic sample"
    );
    // igv.js has drawn its bar and tracks.
    assert!(
        b.wait(
            "document.getElementById('igv').scrollHeight > 150",
            Duration::from_secs(20)
        )
        .is_some(),
        "igv.js drew nothing"
    );
    thread::sleep(Duration::from_millis(500));
    let d = probe(&mut b, &dark("body"));
    assert_eq!(
        get(&d, "dark"),
        &Value::Bool(true),
        "phone: the page must be dark: {d:?}"
    );
    let w = probe(
        &mut b,
        "JSON.stringify({scroll:document.scrollingElement.scrollWidth,client:document.documentElement.clientWidth})",
    );
    let wide = |k: &str| match get(&w, k) {
        Value::Int(n) => *n,
        v => panic!("{k}: {v:?}"),
    };
    assert!(
        wide("client") == i128::from(PHONE.0) && wide("scroll") <= i128::from(PHONE.0),
        "on a {} px phone the page must not be wider than the screen: {w:?}",
        PHONE.0
    );
    shot(&mut b, "phone");
}

#[test]
fn a_config_error_says_what_to_check_and_try_again_starts_the_viewer() {
    let Some(chrome) = chrome_or_skip("load the viewer page with a failing config") else {
        return;
    };
    let scratch = Scratch::new(Path::new(env!("CARGO_TARGET_TMPDIR")), "rv-viewer-error");
    let site = scratch.path().join("site");
    fs::create_dir_all(&site).unwrap();
    let stub = Stub::start(make_sample(&site), Config::Status(500));
    let mut b = browser(&chrome, &scratch.path().join("chrome"), PHONE, true);
    b.go(&stub.url("/#chr1:400-500")).unwrap();
    b.wait(SETTLED, Duration::from_secs(60));
    assert_eq!(text(b.eval(STATUS)), "error: config.json answered 500");
    let e = probe(
        &mut b,
        "(function(){var e=document.getElementById('error');var r=e.getBoundingClientRect();\
         var btn=e.querySelector('button');return JSON.stringify({role:e.getAttribute('role'),\
         shown:!e.hidden&&r.height>0,text:e.textContent,button:btn?btn.textContent:''});})()",
    );
    assert_eq!(
        get(&e, "role"),
        &Value::Str("alert".into()),
        "#error is announced: {e:?}"
    );
    assert_eq!(get(&e, "shown"), &Value::Bool(true), "#error shows: {e:?}");
    let said = text(Ok(get(&e, "text").clone()));
    assert!(
        said.contains("The viewer could not start: config.json answered 500")
            && said.contains("is still running, then try again."),
        "#error names the failure and says what to check: {said}"
    );
    assert_eq!(
        get(&e, "button"),
        &Value::Str("Try again".into()),
        "#error offers Try again: {e:?}"
    );
    let d = probe(&mut b, &dark("#error"));
    assert_eq!(
        get(&d, "dark"),
        &Value::Bool(true),
        "the error box is dark: {d:?}"
    );
    shot(&mut b, "error");
    // The server answers again: Try again starts the viewer.
    stub.set_config(Config::Json(CONFIG.to_owned()));
    b.eval("document.querySelector('#error button').click()")
        .unwrap();
    let ready = b.wait(
        "document.documentElement.dataset.pgpStatus === 'ready' && document.getElementById('error').hidden",
        Duration::from_secs(30),
    );
    assert!(
        ready.is_some(),
        "Try again did not start the viewer: {}",
        text(b.eval(STATUS))
    );
}
