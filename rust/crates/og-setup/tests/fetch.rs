//! `og setup fetch`, `plan` and `prune` against a local HTTP server, and
//! the repository checks that keep setup's downloads in one place: every
//! pin of `assets/data-sources.tsv` written once, no installer deleting a
//! raw download, no source streamed without being saved.
//!
//! The server counts requests per path and can fail the first N of them,
//! or cut a body short, so retries, resumes and reuse are observed, not
//! assumed. `curl` does the transfers, as in setup.

use og_setup::fetch::{
    self, Bad, Ctx, Done, Fetch, Kind, Knobs, OnWrong, Setting, Sources, digest, literal_pins,
    saved_checksum_path, shell_settings, source_for,
};
use std::collections::HashMap;
use std::fs;
use std::io::{BufRead, BufReader, Write};
use std::net::{TcpListener, TcpStream};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

// --- a local HTTP server ---------------------------------------------------

#[derive(Default)]
struct File {
    body: Vec<u8>,
    /// Answer 503 to this many requests first.
    fail_first: usize,
    /// Send only this many bytes of the first full answer, then close.
    cut_first: Option<usize>,
    /// Answer every request with the whole body, as a server that cannot
    /// resume does.
    ignore_range: bool,
}

#[derive(Default)]
struct State {
    files: HashMap<String, File>,
    /// Per path: the Range header of each request ("" for none).
    seen: HashMap<String, Vec<String>>,
}

struct Server {
    base: String,
    state: Arc<Mutex<State>>,
}

impl Server {
    fn start() -> Server {
        let l = TcpListener::bind("127.0.0.1:0").unwrap();
        let base = format!("http://{}", l.local_addr().unwrap());
        let state = Arc::new(Mutex::new(State::default()));
        let s = state.clone();
        std::thread::spawn(move || {
            for c in l.incoming().flatten() {
                let s = s.clone();
                std::thread::spawn(move || serve(c, &s));
            }
        });
        Server { base, state }
    }

    fn put(&self, path: &str, f: File) -> String {
        self.state.lock().unwrap().files.insert(path.to_owned(), f);
        format!("{}{path}", self.base)
    }

    fn requests(&self, path: &str) -> Vec<String> {
        self.state
            .lock()
            .unwrap()
            .seen
            .get(path)
            .cloned()
            .unwrap_or_default()
    }
}

fn serve(mut c: TcpStream, s: &Mutex<State>) {
    let mut r = BufReader::new(c.try_clone().unwrap());
    let mut line = String::new();
    if r.read_line(&mut line).is_err() {
        return;
    }
    // Files are found by path; a query is ignored.
    let target = line.split_whitespace().nth(1).unwrap_or("/");
    let path = target.split('?').next().unwrap_or(target).to_owned();
    let mut range = String::new();
    loop {
        let mut h = String::new();
        if r.read_line(&mut h).is_err() || h.trim().is_empty() {
            break;
        }
        if let Some(v) = h.to_ascii_lowercase().strip_prefix("range:") {
            range = v.trim().to_owned();
        }
    }
    let (status, body, cut) = {
        let mut st = s.lock().unwrap();
        st.seen.entry(path.clone()).or_default().push(range.clone());
        match st.files.get_mut(&path) {
            None => ("404 Not Found", vec![], None),
            Some(f) if f.fail_first > 0 => {
                f.fail_first -= 1;
                ("503 Service Unavailable", vec![], None)
            }
            Some(f) => {
                let from = if f.ignore_range {
                    0
                } else {
                    range
                        .strip_prefix("bytes=")
                        .and_then(|x| x.trim_end_matches('-').parse::<usize>().ok())
                        .unwrap_or(0)
                };
                let cut = f.cut_first.take();
                if from > 0 {
                    (
                        "206 Partial Content",
                        f.body[from.min(f.body.len())..].to_vec(),
                        cut,
                    )
                } else {
                    ("200 OK", f.body.clone(), cut)
                }
            }
        }
    };
    let mut head = format!(
        "HTTP/1.1 {status}\r\nContent-Length: {}\r\nConnection: close\r\n",
        body.len()
    );
    if status.starts_with("206") {
        let total = s.lock().unwrap().files[&path].body.len();
        head.push_str(&format!(
            "Content-Range: bytes {}-{}/{total}\r\n",
            total - body.len(),
            total - 1
        ));
    }
    head.push_str("\r\n");
    let _ = c.write_all(head.as_bytes());
    let n = cut.unwrap_or(body.len()).min(body.len());
    let _ = c.write_all(&body[..n]);
    let _ = c.flush();
}

// --- helpers ---------------------------------------------------------------

fn knobs(tries: u32) -> Knobs {
    Knobs {
        tries,
        wait: Some("0".into()),
        stall: 5,
    }
}

fn tmp(name: &str) -> PathBuf {
    // One folder per call: the tests run in parallel.
    static N: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
    let n = N.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let base = std::env::var_os("TMPDIR").map_or_else(std::env::temp_dir, PathBuf::from);
    let d = base.join(format!("og-setup-fetch-{name}-{}-{n}", std::process::id()));
    let _ = fs::remove_dir_all(&d);
    fs::create_dir_all(&d).unwrap();
    d
}

fn body(n: usize) -> Vec<u8> {
    (0..n).map(|i| (i * 7 % 251) as u8).collect()
}

fn sha256_of(b: &[u8]) -> String {
    let d = tmp("sha");
    let p = d.join("x");
    fs::write(&p, b).unwrap();
    digest(Kind::Sha256, &p).unwrap()
}

fn run(f: &Fetch, k: &Knobs) -> (Result<Done, ()>, String) {
    let mut err = vec![];
    let r = fetch::fetch(f, k, &mut err);
    (r, String::from_utf8(err).unwrap())
}

fn plain<'a>(url: &'a str, dest: &'a Path, kind: Kind, want: &'a str) -> Fetch<'a> {
    Fetch {
        url,
        dest,
        kind,
        want,
        on_wrong: OnWrong::Replace,
        source: None,
        prune: None,
    }
}

/// A URL on a port nothing listens on: a host that is down, and the two ends of the connection
/// that hold the port. The port is the local end of a connection to a listener of our own: a
/// connection to it is refused at once, and the system gives it to no other socket while it is
/// held.
fn down_url(path: &str) -> (String, [TcpStream; 2]) {
    let l = TcpListener::bind("127.0.0.1:0").unwrap();
    let c = TcpStream::connect(l.local_addr().unwrap()).unwrap();
    let url = format!("http://{}{path}", c.local_addr().unwrap());
    (url, [c, l.accept().unwrap().0])
}

// --- the fetch behaviours ---------------------------------------------------

#[test]
fn a_server_that_fails_the_first_tries_is_tried_again() {
    let s = Server::start();
    let b = body(5000);
    let url = s.put(
        "/flaky.bin",
        File {
            body: b.clone(),
            fail_first: 2,
            ..File::default()
        },
    );
    let d = tmp("flaky");
    let dest = d.join("flaky.bin");
    let sha = sha256_of(&b);
    let (r, log) = run(&plain(&url, &dest, Kind::Sha256, &sha), &knobs(3));
    assert_eq!(r, Ok(Done::Downloaded), "{log}");
    assert_eq!(fs::read(&dest).unwrap(), b);
    assert_eq!(s.requests("/flaky.bin").len(), 3);
    assert!(
        log.contains(&format!("  Download attempt 1/3 failed: {url}")),
        "{log}"
    );
    assert!(log.contains("  Download attempt 2/3 failed"), "{log}");
    // Two tries are not enough: the shell's last message, and the .part is kept.
    let url2 = s.put(
        "/down.bin",
        File {
            body: b.clone(),
            fail_first: 9,
            ..File::default()
        },
    );
    let dest2 = d.join("down.bin");
    let (r, log) = run(&plain(&url2, &dest2, Kind::Sha256, &sha), &knobs(2));
    assert_eq!(r, Err(()));
    assert!(
        log.contains(&format!(
            "ERROR: could not download {url2} (the partial file {}.part is kept",
            dest2.display()
        )),
        "{log}"
    );
    assert_eq!(s.requests("/down.bin").len(), 2);
}

#[test]
fn a_broken_transfer_is_resumed_from_its_part() {
    let s = Server::start();
    let b = body(200_000);
    let url = s.put(
        "/big.bin",
        File {
            body: b.clone(),
            cut_first: Some(70_000),
            ..File::default()
        },
    );
    let d = tmp("resume");
    let dest = d.join("big.bin");
    let (r, log) = run(&plain(&url, &dest, Kind::Sha256, &sha256_of(&b)), &knobs(3));
    assert_eq!(r, Ok(Done::Downloaded), "{log}");
    assert_eq!(fs::read(&dest).unwrap(), b);
    let seen = s.requests("/big.bin");
    assert_eq!(seen.len(), 2, "{seen:?}");
    assert_eq!(
        seen[1], "bytes=70000-",
        "the second try resumes the .part: {seen:?}"
    );
    assert!(!d.join("big.bin.part").exists());
}

#[test]
fn a_wrong_checksum_is_refused_never_renamed_and_names_the_source() {
    let s = Server::start();
    let url = s.put(
        "/bad.bin",
        File {
            body: body(3000),
            ..File::default()
        },
    );
    let d = tmp("wrong");
    let dest = d.join("bad.bin");
    let pin = "0".repeat(64);
    let f = Fetch {
        source: Some("somalier_sites"),
        on_wrong: OnWrong::Refuse,
        ..plain(&url, &dest, Kind::Sha256, &pin)
    };
    let (r, log) = run(&f, &knobs(2));
    assert_eq!(r, Err(()));
    assert!(
        !dest.exists(),
        "a download that failed its check was renamed"
    );
    assert!(
        !d.join("bad.bin.part").exists(),
        "the failed download is removed"
    );
    assert!(log.contains("sha256 checksum of bad.bin is"), "{log}");
    assert!(log.contains("Removed the download after 2 tries"), "{log}");
    assert_eq!(s.requests("/bad.bin").len(), 2, "each try downloads anew");
}

#[test]
fn a_file_already_there_with_the_right_checksum_makes_no_request() {
    let s = Server::start();
    let b = body(4000);
    let url = s.put(
        "/here.bin",
        File {
            body: b.clone(),
            ..File::default()
        },
    );
    let d = tmp("reuse");
    let dest = d.join("here.bin");
    fs::write(&dest, &b).unwrap();
    for kind in [Kind::Sha256, Kind::Md5, Kind::Sum, Kind::Crc32c] {
        let want = digest(kind, &dest).unwrap();
        let f = Fetch {
            source: Some("here"),
            on_wrong: OnWrong::Refuse,
            ..plain(&url, &dest, kind, &want)
        };
        let (r, log) = run(&f, &knobs(3));
        assert_eq!(r, Ok(Done::Reused), "{kind:?}: {log}");
        assert!(log.contains("is already here and matches its"), "{log}");
    }
    assert!(
        s.requests("/here.bin").is_empty(),
        "a matching file was requested again"
    );
    // A checksum read from a file:// list costs no request either.
    let list = d.join("md5sums.txt");
    fs::write(
        &list,
        format!("{}  here.bin\n", digest(Kind::Md5, &dest).unwrap()),
    )
    .unwrap();
    let want = format!("file://{}", list.display());
    let (r, log) = run(&plain(&url, &dest, Kind::Md5, &want), &knobs(3));
    assert_eq!(r, Ok(Done::Reused), "{log}");
    assert!(s.requests("/here.bin").is_empty());
}

#[test]
fn a_file_there_with_the_wrong_checksum_is_never_overwritten_silently() {
    let s = Server::start();
    let b = body(4000);
    let url = s.put(
        "/mine.bin",
        File {
            body: b.clone(),
            ..File::default()
        },
    );
    let d = tmp("refuse");
    let dest = d.join("mine.bin");
    fs::write(&dest, b"a file the user put here").unwrap();
    let sha = sha256_of(&b);
    // In the downloads folder: refused, left as it is, no request.
    let f = Fetch {
        source: Some("mine"),
        prune: Some("mine"),
        on_wrong: OnWrong::Refuse,
        ..plain(&url, &dest, Kind::Sha256, &sha)
    };
    let (r, log) = run(&f, &knobs(3));
    assert_eq!(r, Err(()));
    assert_eq!(fs::read(&dest).unwrap(), b"a file the user put here");
    assert!(
        log.contains("is there but is not the file of mine"),
        "{log}"
    );
    assert!(log.contains("(or run og setup prune mine)"), "{log}");
    assert!(s.requests("/mine.bin").is_empty());
    // A file of no source of one file: no prune hint (prune would refuse it).
    let f = Fetch {
        on_wrong: OnWrong::Refuse,
        ..plain(&url, &dest, Kind::Sha256, &sha)
    };
    let (r, log) = run(&f, &knobs(3));
    assert_eq!(r, Err(()));
    assert!(
        log.contains("is not the file of mine.bin") && !log.contains("og setup prune"),
        "{log}"
    );
    // A file setup installs itself: replaced, and it says so.
    let (r, log) = run(&plain(&url, &dest, Kind::Sha256, &sha), &knobs(3));
    assert_eq!(r, Ok(Done::Downloaded), "{log}");
    assert!(
        log.contains("is not the pinned file (sha256 checksum of mine.bin is"),
        "{log}"
    );
    assert_eq!(fs::read(&dest).unwrap(), b);
}

#[test]
fn without_a_checksum_a_file_there_proves_nothing_and_is_downloaded_again() {
    let s = Server::start();
    let url = s.put(
        "/plain.txt",
        File {
            body: b"new\n".to_vec(),
            ..File::default()
        },
    );
    let d = tmp("nocheck");
    let dest = d.join("plain.txt");
    fs::write(&dest, b"old\n").unwrap();
    let (r, _) = run(&plain(&url, &dest, Kind::None, ""), &knobs(1));
    assert_eq!(r, Ok(Done::Downloaded));
    assert_eq!(fs::read(&dest).unwrap(), b"new\n");
    assert_eq!(s.requests("/plain.txt").len(), 1);
}

#[test]
fn in_the_downloads_folder_a_whole_file_with_no_checksum_is_used_and_a_broken_one_refused() {
    let s = Server::start();
    let url = s.put(
        "/annotsv.tar.gz",
        File {
            body: b"new".to_vec(),
            ..File::default()
        },
    );
    let d = tmp("dlnocheck");
    let gz = d.join("annotsv.tar.gz");
    let mut e = flate2::write::GzEncoder::new(vec![], flate2::Compression::default());
    e.write_all(b"a bundle a user brought\n").unwrap();
    let whole = e.finish().unwrap();
    fs::write(&gz, &whole).unwrap();
    let refuse = |dest| Fetch {
        on_wrong: OnWrong::Refuse,
        ..plain(&url, dest, Kind::None, "")
    };
    let (r, log) = run(&refuse(&gz), &knobs(3));
    assert_eq!(r, Ok(Done::Reused), "{log}");
    assert!(
        log.contains("matches its gzip check; not downloaded"),
        "{log}"
    );
    assert!(
        s.requests("/annotsv.tar.gz").is_empty(),
        "a whole gzip there was requested again"
    );
    // A cut one is not the file: refused and left, never replaced silently.
    fs::write(&gz, &whole[..10]).unwrap();
    let (r, log) = run(&refuse(&gz), &knobs(3));
    assert_eq!(r, Err(()));
    assert!(
        log.contains("is not a complete gzip file. It was left as it is"),
        "{log}"
    );
    assert_eq!(fs::read(&gz).unwrap(), &whole[..10]);
    assert!(s.requests("/annotsv.tar.gz").is_empty());
}

#[test]
fn a_url_dest_call_inside_the_downloads_folder_refuses_outside_it_replaces() {
    let env = |dl: Option<&'static str>, g: Option<&'static str>| {
        move |k: &str| match k {
            "OG_DOWNLOADS" => dl.map(str::to_owned),
            "GENOME_DIR" => g.map(str::to_owned),
            _ => None,
        }
    };
    assert!(fetch::in_downloads(
        Path::new("/g/downloads/x.zip"),
        &env(None, Some("/g"))
    ));
    assert!(!fetch::in_downloads(
        Path::new("/g/reference/x.zip"),
        &env(None, Some("/g"))
    ));
    assert!(fetch::in_downloads(
        Path::new("/shared/dl/a/b.gz"),
        &env(Some("/shared/dl"), Some("/g"))
    ));
    assert!(
        !fetch::in_downloads(
            Path::new("/g/downloads/x.zip"),
            &env(Some("/shared/dl"), Some("/g"))
        ),
        "OG_DOWNLOADS moves the folder"
    );
    assert!(
        !fetch::in_downloads(
            Path::new("/shared/dlx/b.gz"),
            &env(Some("/shared/dl"), None)
        ),
        "a sibling with a longer name is outside"
    );
    assert!(!fetch::in_downloads(
        Path::new("/g/downloads/x"),
        &env(None, None)
    ));
}

#[test]
fn a_checksum_list_at_a_url_names_the_file() {
    let s = Server::start();
    let b = body(1234);
    let url = s.put(
        "/d/x.gz.bin",
        File {
            body: b.clone(),
            ..File::default()
        },
    );
    let d = tmp("list");
    let probe = d.join("probe");
    fs::write(&probe, &b).unwrap();
    let md5 = digest(Kind::Md5, &probe).unwrap();
    let list = s.put(
        "/d/MD5SUMS",
        File {
            body: format!("ffff  other\n{md5}  ./x.gz.bin\n").into_bytes(),
            ..File::default()
        },
    );
    let dest = d.join("x.gz.bin");
    let (r, log) = run(&plain(&url, &dest, Kind::Md5, &list), &knobs(1));
    assert_eq!(r, Ok(Done::Downloaded), "{log}");
    let none = s.put(
        "/d/EMPTY",
        File {
            body: b"a  b\nc  d\n".to_vec(),
            ..File::default()
        },
    );
    let (r, log) = run(&plain(&url, &d.join("y"), Kind::Md5, &none), &knobs(3));
    assert_eq!(r, Err(()));
    assert!(
        log.contains("ERROR: could not read the checksum of x.gz.bin from"),
        "{log}"
    );
    assert_eq!(
        s.requests("/d/EMPTY").len(),
        1,
        "a list that names no checksum is not read again"
    );
}

#[test]
fn a_server_that_cannot_resume_gets_one_fresh_download() {
    let s = Server::start();
    let b = body(6000);
    let url = s.put(
        "/norange.bin",
        File {
            body: b.clone(),
            ignore_range: true,
            ..File::default()
        },
    );
    let d = tmp("norange");
    let dest = d.join("norange.bin");
    // A part left by a cut transfer: curl asks to resume it, the server sends
    // the whole file (200, not 206), and curl stops with exit 33.
    fs::write(d.join("norange.bin.part"), &b[..100]).unwrap();
    let (r, log) = run(&plain(&url, &dest, Kind::Sha256, &sha256_of(&b)), &knobs(1));
    assert_eq!(r, Ok(Done::Downloaded), "{log}");
    assert_eq!(fs::read(&dest).unwrap(), b);
    assert_eq!(
        s.requests("/norange.bin"),
        ["bytes=100-", ""],
        "one resume, refused, then one fresh download"
    );
}

#[test]
fn a_checksum_list_that_fails_first_is_tried_again() {
    let s = Server::start();
    let b = body(1500);
    let url = s.put(
        "/r/x.gz.bin",
        File {
            body: b.clone(),
            ..File::default()
        },
    );
    let d = tmp("listretry");
    let probe = d.join("probe");
    fs::write(&probe, &b).unwrap();
    let md5 = digest(Kind::Md5, &probe).unwrap();
    let list = s.put(
        "/r/MD5SUMS",
        File {
            body: format!("{md5}  x.gz.bin\n").into_bytes(),
            fail_first: 2,
            ..File::default()
        },
    );
    let (r, log) = run(
        &plain(&url, &d.join("x.gz.bin"), Kind::Md5, &list),
        &knobs(3),
    );
    assert_eq!(r, Ok(Done::Downloaded), "{log}");
    assert_eq!(s.requests("/r/MD5SUMS").len(), 3, "{log}");
    assert!(
        log.contains("Checksum download attempt 2/3 failed"),
        "{log}"
    );
}

#[test]
fn a_kept_file_checked_against_a_list_is_used_when_the_list_cannot_be_read() {
    let s = Server::start();
    let b = body(2500);
    let url = s.put(
        "/k/x.gz.bin",
        File {
            body: b.clone(),
            ..File::default()
        },
    );
    let d = tmp("savedsum");
    let probe = d.join("probe");
    fs::write(&probe, &b).unwrap();
    let md5 = digest(Kind::Md5, &probe).unwrap();
    let list = s.put(
        "/k/MD5SUMS",
        File {
            body: format!("{md5}  x.gz.bin\n").into_bytes(),
            ..File::default()
        },
    );
    let dest = d.join("downloads/x.gz.bin");
    fn kept<'a>(url: &'a str, dest: &'a Path, want: &'a str) -> Fetch<'a> {
        Fetch {
            on_wrong: OnWrong::Refuse,
            ..plain(url, dest, Kind::Md5, want)
        }
    }
    // Downloaded and checked: the checksum it matched is saved beside it.
    let (r, log) = run(&kept(&url, &dest, &list), &knobs(1));
    assert_eq!(r, Ok(Done::Downloaded), "{log}");
    assert_eq!(
        fs::read_to_string(saved_checksum_path(&dest)).unwrap(),
        format!("md5 {md5}\n")
    );
    // The list up: one small request for the list, none for the file.
    let (r, log) = run(&kept(&url, &dest, &list), &knobs(1));
    assert_eq!(r, Ok(Done::Reused), "{log}");
    assert_eq!(s.requests("/k/MD5SUMS").len(), 2);
    assert_eq!(s.requests("/k/x.gz.bin").len(), 1);
    // The list's host down: the kept file matches what was saved, so it is used.
    let (gone, _down) = down_url("/k/MD5SUMS");
    let (r, log) = run(&kept(&url, &dest, &gone), &knobs(2));
    assert_eq!(r, Ok(Done::Reused), "{log}");
    assert!(
        log.contains("could not be read;")
            && log.contains("matches the md5 saved when it was checked"),
        "{log}"
    );
    assert_eq!(s.requests("/k/x.gz.bin").len(), 1, "{log}");
    // A kept file changed since: the saved checksum proves nothing, refused.
    fs::write(&dest, b"changed").unwrap();
    let (r, log) = run(&kept(&url, &dest, &gone), &knobs(1));
    assert_eq!(r, Err(()), "{log}");
    assert!(
        log.contains("ERROR: could not read the checksum of x.gz.bin"),
        "{log}"
    );
    assert_eq!(fs::read(&dest).unwrap(), b"changed");
    // No saved checksum (a file kept by an older setup): refused too.
    fs::write(&dest, &b).unwrap();
    fs::remove_file(saved_checksum_path(&dest)).unwrap();
    let (r, log) = run(&kept(&url, &dest, &gone), &knobs(1));
    assert_eq!(r, Err(()), "{log}");
    // The list up again: checked, and the checksum saved again.
    let (r, log) = run(&kept(&url, &dest, &list), &knobs(1));
    assert_eq!(r, Ok(Done::Reused), "{log}");
    assert!(saved_checksum_path(&dest).is_file());
    // A file setup installs itself (not kept) saves nothing.
    let other = d.join("installed/x.gz.bin");
    let (r, log) = run(&plain(&url, &other, Kind::Md5, &list), &knobs(1));
    assert_eq!(r, Ok(Done::Downloaded), "{log}");
    assert!(!saved_checksum_path(&other).exists());
}

#[test]
fn a_local_checksum_list_is_read_as_a_file_even_with_a_space() {
    let s = Server::start();
    let b = body(900);
    let url = s.put(
        "/l/x.gz.bin",
        File {
            body: b.clone(),
            ..File::default()
        },
    );
    let d = tmp("local").join("Reference Data");
    fs::create_dir_all(&d).unwrap();
    let probe = d.join("probe");
    fs::write(&probe, &b).unwrap();
    let md5 = digest(Kind::Md5, &probe).unwrap();
    let list = d.join("md5.txt");
    fs::write(&list, format!("ffff  other\n{md5}  x.gz.bin\n")).unwrap();
    let want = format!("file://{}", list.display());
    let dest = d.join("x.gz.bin");
    let (r, log) = run(&plain(&url, &dest, Kind::Md5, &want), &knobs(1));
    assert_eq!(r, Ok(Done::Downloaded), "{log}");
    let (r, log) = run(&plain(&url, &dest, Kind::Md5, &want), &knobs(1));
    assert_eq!(r, Ok(Done::Reused), "{log}");
    assert_eq!(s.requests("/l/x.gz.bin").len(), 1);
    let gone = format!("file://{}", d.join("nope.txt").display());
    let (r, log) = run(&plain(&url, &d.join("y"), Kind::Md5, &gone), &knobs(3));
    assert_eq!(r, Err(()));
    assert!(log.contains("cannot read the checksum list"), "{log}");
}

#[test]
fn a_checksum_is_named_by_the_url_path_without_its_query() {
    let s = Server::start();
    let b = body(700);
    let url = s.put(
        "/q/x.gz.bin",
        File {
            body: b.clone(),
            ..File::default()
        },
    );
    let d = tmp("query");
    let probe = d.join("probe");
    fs::write(&probe, &b).unwrap();
    let md5 = digest(Kind::Md5, &probe).unwrap();
    let list = d.join("md5.txt");
    fs::write(&list, format!("ffff  other\n{md5}  x.gz.bin\n")).unwrap();
    let want = format!("file://{}", list.display());
    let with_query = format!("{url}?download=1");
    let (r, log) = run(
        &plain(&with_query, &d.join("x.gz.bin"), Kind::Md5, &want),
        &knobs(1),
    );
    assert_eq!(
        r,
        Ok(Done::Downloaded),
        "the checksum of x.gz.bin, not of x.gz.bin?download=1: {log}"
    );
}

#[test]
fn an_empty_checksum_is_refused_before_any_request() {
    let s = Server::start();
    let url = s.put(
        "/e/x",
        File {
            body: body(10),
            ..File::default()
        },
    );
    let d = tmp("empty");
    let (r, log) = run(&plain(&url, &d.join("x"), Kind::Md5, ""), &knobs(3));
    assert_eq!(r, Err(()));
    assert!(log.contains("no md5 value given"), "{log}");
    assert!(s.requests("/e/x").is_empty(), "{log}");
}

#[test]
fn a_file_that_appears_during_a_download_is_never_replaced() {
    let d = tmp("publish");
    let (part, dest) = (d.join("x.part"), d.join("x"));
    fs::write(&part, b"checked").unwrap();
    fs::write(&dest, b"the user's").unwrap();
    let e = fetch::publish_new(&part, &dest).unwrap_err();
    assert_eq!(e.kind(), std::io::ErrorKind::AlreadyExists);
    assert_eq!(fs::read(&dest).unwrap(), b"the user's");
    assert_eq!(
        fs::read(&part).unwrap(),
        b"checked",
        "the checked file stays"
    );
    fs::remove_file(&dest).unwrap();
    fetch::publish_new(&part, &dest).unwrap();
    assert_eq!(fs::read(&dest).unwrap(), b"checked");
    assert!(!part.exists());
}

#[test]
fn two_fetches_of_one_download_take_turns() {
    let s = Server::start();
    let b = body(5000);
    let sum = sha256_of(&b);
    let url = s.put(
        "/t/x.bin",
        File {
            body: b.clone(),
            ..File::default()
        },
    );
    let d = tmp("turns");
    let dest = d.join("x.bin");
    // Another fetch holds the lock.
    let held = fs::File::create(d.join(".x.bin.lock")).unwrap();
    held.lock().unwrap();
    let (u, p, w) = (url.clone(), dest.clone(), sum.clone());
    let t = std::thread::spawn(move || {
        let mut f = plain(&u, &p, Kind::Sha256, &w);
        f.on_wrong = OnWrong::Refuse;
        let mut err = vec![];
        let r = fetch::fetch(&f, &knobs(1), &mut err);
        (r, String::from_utf8(err).unwrap())
    });
    std::thread::sleep(std::time::Duration::from_millis(500));
    assert!(!t.is_finished(), "the second fetch waits for the first");
    assert!(s.requests("/t/x.bin").is_empty());
    // The first one publishes the file and lets go.
    fs::write(&dest, &b).unwrap();
    held.unlock().unwrap();
    let (r, log) = t.join().unwrap();
    assert_eq!(r, Ok(Done::Reused), "{log}");
    assert!(log.contains("another fetch of"), "{log}");
    assert!(s.requests("/t/x.bin").is_empty(), "{log}");
}

#[test]
fn the_downloads_folder_is_found_through_dot_dot_and_symlinks() {
    let d = tmp("resolve");
    fs::create_dir_all(d.join("staging")).unwrap();
    fs::create_dir_all(d.join("downloads")).unwrap();
    std::os::unix::fs::symlink(d.join("downloads"), d.join("link")).unwrap();
    let dl = d.join("downloads").display().to_string();
    let env = move |k: &str| (k == "OG_DOWNLOADS").then(|| dl.clone());
    assert!(fetch::in_downloads(
        &d.join("staging/../downloads/new/x.gz"),
        &env
    ));
    assert!(fetch::in_downloads(&d.join("link/x.gz"), &env));
    assert!(!fetch::in_downloads(
        &d.join("downloads/../staging/x.gz"),
        &env
    ));
}

#[test]
fn checks_are_the_shell_checks() {
    let d = tmp("checks");
    let gz = d.join("a.gz");
    let mut e = flate2::write::GzEncoder::new(vec![], flate2::Compression::default());
    e.write_all(b"hello\n").unwrap();
    fs::write(&gz, e.finish().unwrap()).unwrap();
    assert_eq!(fetch::check(Kind::None, "", &gz, &gz, "a.gz"), Ok(()));
    let cut = d.join("b.gz");
    fs::write(&cut, &fs::read(&gz).unwrap()[..10]).unwrap();
    assert_eq!(
        fetch::check(Kind::None, "", &cut, &cut, "b.gz"),
        Err(Bad::Wrong("b.gz is not a complete gzip file".into()))
    );
    let empty = d.join("e.txt");
    fs::write(&empty, b"").unwrap();
    assert_eq!(
        fetch::check(Kind::None, "", &empty, &empty, "e.txt"),
        Err(Bad::Wrong("e.txt downloaded empty".into()))
    );
    assert!(matches!(
        fetch::check(Kind::Sum, "12a", &gz, &gz, "a.gz"),
        Err(Bad::Fatal(_))
    ));
    let sum = digest(Kind::Sum, &gz).unwrap();
    assert_eq!(
        fetch::check(Kind::Sum, &format!("000{sum}"), &gz, &gz, "a.gz"),
        Ok(())
    );
    let ends = format!(
        "{}:{}:{}",
        fs::metadata(&gz).unwrap().len(),
        digest(Kind::Sha256, &gz).unwrap(),
        digest(Kind::Sha256, &gz).unwrap()
    );
    assert_eq!(
        fetch::check(Kind::GzipEnds, &ends, &gz, &gz, "a.gz"),
        Ok(())
    );
}

#[test]
fn bsd_sum_is_the_sum_command() {
    // `printf 'abc' | sum` prints 00085 on GNU and BSD alike: 98+99 after a's rotations.
    let d = tmp("sum");
    let p = d.join("abc");
    fs::write(&p, b"abc").unwrap();
    let got = digest(Kind::Sum, &p).unwrap();
    let out = std::process::Command::new("sum").arg(&p).output().unwrap();
    let want: u64 = String::from_utf8_lossy(&out.stdout)
        .split_whitespace()
        .next()
        .unwrap()
        .parse()
        .unwrap();
    assert_eq!(got, want.to_string());
}

// --- the table, plan and prune ---------------------------------------------

const TABLE_HEAD: &str =
    "id\tinstallers\tsteps\tpath\tparts\turl\tcheck\tvalue\tbytes\tlicence\tinstalled\n";

fn checkout(rows: &str, versions: &str) -> PathBuf {
    let d = tmp("checkout");
    fs::create_dir_all(d.join("assets")).unwrap();
    fs::write(
        d.join("assets/data-sources.tsv"),
        format!("# a table\n{TABLE_HEAD}{rows}"),
    )
    .unwrap();
    fs::write(d.join("versions.env"), versions).unwrap();
    d
}

fn ctx(root: &Path, env: &[(&str, String)]) -> Ctx {
    let env: HashMap<String, String> = env
        .iter()
        .map(|(k, v)| ((*k).to_owned(), v.clone()))
        .collect();
    Ctx::new(root, Box::new(move |k| env.get(k).cloned())).unwrap()
}

#[test]
fn rows_expand_pins_parts_and_the_downloads_folder() {
    let root = checkout(
        "vep\tsetup,vep-cache\t13\tdownloads/vep_${VEP}.tar.gz\t\thttps://h/release-${VEP}/v.tar.gz\tsum\thttps://h/release-${VEP}/CHECKSUMS\t26000000000\tEnsembl\tvep/${VEP}/info.txt\n\
         ref\tsetup\t02\treference/ref.fna.gz\t\t${REF_URL:-https://h/ref.fna.gz}\tmd5\thttps://h/md5\t900000000\tNCBI\t-\n\
         chr\tarchaic\t43\tdownloads/arch/{part}.vcf.gz\t21..22\thttps://h/chr{part}.vcf.gz\tmd5\t@list\t\tMPI\tarch/chr{part}.done\n\
         list\tarchaic\t43\tdownloads/arch/md5.txt\t\thttps://h/md5.txt\tsha256\tabc\t10\tMPI\tarch/DONE\n",
        "VEP=\"116\"\n",
    );
    let s = Sources::read(&root).unwrap();
    let g = root.join("g");
    let c = ctx(&root, &[("GENOME_DIR", g.display().to_string())]);
    assert_eq!(
        c.path(s.get("vep").unwrap(), None).unwrap(),
        g.join("downloads/vep_116.tar.gz")
    );
    assert_eq!(
        c.url(s.get("vep").unwrap(), None).unwrap(),
        "https://h/release-116/v.tar.gz"
    );
    assert_eq!(
        c.url(s.get("ref").unwrap(), None).unwrap(),
        "https://h/ref.fna.gz"
    );
    let o = ctx(
        &root,
        &[
            ("GENOME_DIR", g.display().to_string()),
            ("REF_URL", "https://m/r.gz".into()),
            ("OG_DOWNLOADS", "/shared".into()),
        ],
    );
    assert_eq!(
        o.url(s.get("ref").unwrap(), None).unwrap(),
        "https://m/r.gz",
        "the environment override wins"
    );
    assert_eq!(
        o.path(s.get("vep").unwrap(), None).unwrap(),
        PathBuf::from("/shared/vep_116.tar.gz")
    );
    assert_eq!(
        o.path(s.get("ref").unwrap(), None).unwrap(),
        g.join("reference/ref.fna.gz"),
        "an installed file stays in GENOME_DIR"
    );
    assert_eq!(
        c.url(s.get("chr").unwrap(), Some("22")).unwrap(),
        "https://h/chr22.vcf.gz"
    );
    assert!(c.url(s.get("chr").unwrap(), Some("23")).is_err());
    assert!(c.url(s.get("chr").unwrap(), None).is_err());
    assert_eq!(
        c.value(&s, s.get("chr").unwrap(), Some("21")).unwrap(),
        format!("file://{}", g.join("downloads/arch/md5.txt").display())
    );
    let bad = checkout(
        "x\tsetup\t1\tdownloads/x\t\thttps://h/${NOPE}\t\t\t\t\tx.done\n",
        "",
    );
    let s2 = Sources::read(&bad).unwrap();
    assert!(
        ctx(&bad, &[("GENOME_DIR", "/g".into())])
            .url(s2.get("x").unwrap(), None)
            .unwrap_err()
            .contains("NOPE is not a pin")
    );
}

#[test]
fn the_table_refuses_a_repeated_id_or_path_and_a_bad_check() {
    let r = "a\tsetup\t1\tdownloads/a\t\thttps://h/a\t\t\t\t\ta.done\n";
    let again = "a\tsetup\t1\tdownloads/a2\t\thttps://h/a2\t\t\t\t\ta.done\n";
    assert!(
        Sources::parse(&format!("{TABLE_HEAD}{r}{again}"))
            .unwrap_err()
            .contains("twice")
    );
    let p = "b\tsetup\t1\tdownloads/a\t\thttps://h/b\t\t\t\t\tb.done\n";
    assert!(
        Sources::parse(&format!("{TABLE_HEAD}{r}{p}"))
            .unwrap_err()
            .contains("also a's")
    );
    let k = "c\tsetup\t1\tdownloads/c\t\thttps://h/c\tcrc64\t\t\t\tc.done\n";
    assert!(
        Sources::parse(&format!("{TABLE_HEAD}{k}"))
            .unwrap_err()
            .contains("unknown checksum type")
    );
    let up = "d\tsetup\t1\t../d\t\thttps://h/d\t\t\t\t\t-\n";
    assert!(Sources::parse(&format!("{TABLE_HEAD}{up}")).is_err());
}

#[test]
fn plan_lists_rows_sizes_state_kept_raw_size_and_free_disk() {
    let root = checkout(
        "big\tsetup\t13\tdownloads/big.tar.gz\t\thttps://h/big.tar.gz\t\t\t26000000000\tEnsembl\tbig/DONE\n\
         small\tsetup\t33\treference/small.vcf.gz\t\thttps://h/small.vcf.gz\t\t\t9800000\tMIT\t-\n\
         fam\tsetup\t43\tdownloads/f/{part}.gz\t1..4\thttps://h/{part}.gz\t\t\t4000\tMPI\tf/{part}.done\n\
         other\tarchaic\t43\tdownloads/o.gz\t\thttps://h/o.gz\t\t\t1\tMPI\to.done\n",
        "",
    );
    let g = root.join("g");
    fs::create_dir_all(g.join("reference")).unwrap();
    fs::create_dir_all(g.join("downloads/f")).unwrap();
    fs::write(g.join("reference/small.vcf.gz"), b"x").unwrap();
    fs::write(g.join("downloads/f/1.gz"), b"x").unwrap();
    let c = ctx(&root, &[("GENOME_DIR", g.display().to_string())]);
    let s = Sources::read(&root).unwrap();
    let rows = og_setup::plan::chosen(&s, &["setup".to_owned()]);
    let lines = og_setup::plan::lines(&c, &s, &rows).unwrap();
    let mut out = vec![];
    og_setup::plan::write_plan(&c, &lines, &mut out).unwrap();
    let out = String::from_utf8(out).unwrap();
    let want = format!(
        "Downloads folder: {}\nsource\tfile\tsize\tstate\n\
         big\tdownloads/big.tar.gz\t26.0 GB\tmissing\n\
         small\treference/small.vcf.gz\t9.8 MB\tpresent\n\
         fam\tdownloads/f/{{part}}.gz\t4.0 kB\t1 of 4 files present\n\
         To download: 2 sources, 26.0 GB; already on disk: 9.8 MB.\n\
         Raw downloads kept after the install: 26.0 GB (--no-keep-downloads or SETUP_KEEP_DOWNLOADS=false deletes each once installed; og setup prune deletes them later).\n",
        g.join("downloads").display()
    );
    assert!(
        out.starts_with(&want),
        "{out}\n---- wanted the start:\n{want}"
    );
    assert!(out.contains("Free disk: "), "{out}");
    assert!(
        !out.contains("other"),
        "a row of another installer is listed: {out}"
    );
    // Nothing but the plan is free disk: 26 GB is more than a temp folder has
    // on a test machine only by luck, so the warning is checked through the
    // arithmetic instead.
    assert_eq!(
        og_setup::plan::shortfall(26_000_000_000, 1_000),
        Some(25_999_999_000)
    );
    assert_eq!(og_setup::plan::shortfall(1, 1_000), None);
}

#[test]
fn the_table_refuses_a_download_with_no_install_and_a_bad_installed_file() {
    let row = |path: &str, parts: &str, installed: &str| {
        format!("{TABLE_HEAD}a\tsetup\t1\t{path}\t{parts}\thttps://h/a\t\t\t\t\t{installed}\n")
    };
    assert!(
        Sources::parse(&row("downloads/a", "", "-"))
            .unwrap_err()
            .contains("is a raw download: name the files of its install")
    );
    assert!(
        Sources::parse(&row("reference/a", "", "reference/a.done"))
            .unwrap_err()
            .contains("'installed' is for a raw download")
    );
    assert!(
        Sources::parse(&row("downloads/a", "", "a/{part}.done"))
            .unwrap_err()
            .contains("{part} in 'installed' goes with a parts list")
    );
    for bad in ["/a.done", "downloads/a.done", "x/../../a.done"] {
        assert!(
            Sources::parse(&row("downloads/a", "", bad))
                .unwrap_err()
                .contains("must be relative, in GENOME_DIR, outside downloads/"),
            "{bad}"
        );
    }
    let s = Sources::parse(&row("downloads/a", "", "a/x a/y")).unwrap();
    assert_eq!(s.get("a").unwrap().installed, ["a/x", "a/y"]);
}

/// `og setup plan ARGS` on ROOT with ENV: (exit code, stdout, stderr).
fn plan_cli(root: &Path, args: &[&str], env: &[(&str, String)]) -> (u8, String, String) {
    let env: HashMap<String, String> = env
        .iter()
        .map(|(k, v)| ((*k).to_owned(), v.clone()))
        .collect();
    let mut a: Vec<String> = args.iter().map(|x| (*x).to_owned()).collect();
    a.extend(["--root".to_owned(), root.display().to_string()]);
    let (mut out, mut err) = (vec![], vec![]);
    let rc = og_setup::plan::cli_with(
        &a,
        Box::new(move |k| env.get(k).cloned()),
        &mut out,
        &mut err,
    );
    (
        rc,
        String::from_utf8(out).unwrap(),
        String::from_utf8(err).unwrap(),
    )
}

#[test]
fn plan_check_refuses_downloads_that_do_not_fit_and_nothing_else() {
    // 10^18 bytes fit on no disk.
    let root = checkout(
        "huge\tsetup\t13\tdownloads/huge.tar.gz\t\thttps://h/huge.tar.gz\t\t\t1000000000000000000\tX\thuge/DONE\n\
         big_installed\tbig\t13\treference/huge.vcf.gz\t\thttps://h/huge.vcf.gz\t\t\t1000000000000000000\tX\t-\n",
        "",
    );
    let g = root.join("g");
    let genome = ("GENOME_DIR", g.display().to_string());
    let dl = g.join("downloads");
    // Not there: refused, the shortfall named.
    let (rc, out, err) = plan_cli(&root, &["--check"], std::slice::from_ref(&genome));
    assert_eq!(rc, 1, "{out}{err}");
    assert!(
        err.contains(
            "ERROR: the raw downloads this setup still needs, 1000000.0 TB, do not fit: they need "
        ) && err.contains(&format!("TB more than is free in {}", dl.display()))
            && err.contains("--no-keep-downloads"),
        "{err}"
    );
    assert!(
        out.contains("huge\tdownloads/huge.tar.gz\t1000000.0 TB\tmissing"),
        "{out}"
    );
    // Deleted once installed: nothing is kept, so nothing is refused.
    let (rc, _, err) = plan_cli(
        &root,
        &["--check"],
        &[genome.clone(), ("SETUP_KEEP_DOWNLOADS", "false".into())],
    );
    assert_eq!(rc, 0, "{err}");
    // An install with no download (a file setup installs, not kept): not counted.
    let (rc, out, err) = plan_cli(
        &root,
        &["--for", "big", "--check"],
        std::slice::from_ref(&genome),
    );
    assert_eq!(rc, 0, "{out}{err}");
    assert!(out.contains("To download: 1 source, 1000000.0 TB"), "{out}");
    // Its install done by an older setup, which deleted the download, and no
    // downloads folder at all: nothing to download.
    fs::create_dir_all(g.join("huge")).unwrap();
    fs::write(g.join("huge/DONE"), b"").unwrap();
    assert!(!dl.exists());
    let (rc, out, err) = plan_cli(&root, &["--check"], std::slice::from_ref(&genome));
    assert_eq!(rc, 0, "{out}{err}");
    assert!(
        out.contains("huge\tdownloads/huge.tar.gz\t1000000.0 TB\tinstalled, not kept"),
        "{out}"
    );
    assert!(out.contains("To download: 0 sources, 0 bytes"), "{out}");
    // The file there: nothing to download.
    fs::remove_file(g.join("huge/DONE")).unwrap();
    fs::create_dir_all(&dl).unwrap();
    fs::write(dl.join("huge.tar.gz"), b"x").unwrap();
    let (rc, out, err) = plan_cli(&root, &["--check"], std::slice::from_ref(&genome));
    assert_eq!(rc, 0, "{out}{err}");
    assert!(
        out.contains("huge\tdownloads/huge.tar.gz\t1000000.0 TB\tpresent"),
        "{out}"
    );
}

/// A checkout holding the real table, versions.env and the two scripts whose
/// settings it names, with EDIT applied to the file named first.
fn real_checkout(edit: Option<(&str, &str, &str)>) -> PathBuf {
    let d = tmp("real");
    for f in [
        "assets/data-sources.tsv",
        "versions.env",
        "scripts/lib/common.sh",
        "scripts/setup.sh",
    ] {
        let mut t = read(f);
        if let Some((file, from, to)) = edit
            && file == f
        {
            assert_eq!(t.matches(from).count(), 1, "{from} in {f}");
            t = t.replace(from, to);
        }
        fs::create_dir_all(d.join(f).parent().unwrap()).unwrap();
        fs::write(d.join(f), t).unwrap();
    }
    d
}

#[test]
fn a_rerun_of_an_install_built_by_the_shell_setup_needs_no_vep_download() {
    // The VEP cache row made too large for any disk; the release is versions.env's.
    let t = read("assets/data-sources.tsv");
    let row = t.lines().find(|l| l.starts_with("vep_cache\t")).unwrap();
    let mut f: Vec<&str> = row.split('\t').collect();
    f[8] = "1000000000000000000";
    let big = f.join("\t");
    let root = real_checkout(Some(("assets/data-sources.tsv", row, &big)));
    let rel = fetch::versions_pins(&read("versions.env"))["VEP_CACHE_RELEASE"].clone();
    let g = tmp("vepdone");
    let genome = ("GENOME_DIR", g.display().to_string());
    let (rc, _, err) = plan_cli(
        &root,
        &["--for", "vep-cache", "--check"],
        std::slice::from_ref(&genome),
    );
    assert_eq!(rc, 1, "with no cache the download is needed: {err}");
    // What install_vep_cache_opt_in checks, as an older setup left it: the
    // unpacked cache, no tarball and no downloads folder.
    let info = g.join(format!("vep_cache/homo_sapiens/{rel}_GRCh38/info.txt"));
    fs::create_dir_all(info.parent().unwrap()).unwrap();
    fs::write(&info, b"").unwrap();
    let (rc, out, err) = plan_cli(
        &root,
        &["--for", "vep-cache", "--check"],
        std::slice::from_ref(&genome),
    );
    assert_eq!(rc, 0, "{out}{err}");
    assert!(out.contains("\tinstalled, not kept"), "{out}");
    assert!(!g.join("downloads").exists());
}

/// The real table with row ID made too large for any disk.
fn checkout_with_huge(id: &str) -> PathBuf {
    let t = read("assets/data-sources.tsv");
    let row = t
        .lines()
        .find(|l| l.starts_with(&format!("{id}\t")))
        .unwrap();
    let mut f: Vec<&str> = row.split('\t').collect();
    f[8] = "1000000000000000000";
    let big = f.join("\t");
    real_checkout(Some(("assets/data-sources.tsv", row, &big)))
}

#[test]
fn an_archaic_folder_holding_only_the_sites_file_and_a_lock_is_still_to_download() {
    let root = checkout_with_huge("archaic_vcf");
    let g = tmp("archaic");
    let genome = ("GENOME_DIR", g.display().to_string());
    let c = ctx(&root, std::slice::from_ref(&genome));
    let x = |s: &str| c.expand(s, None).unwrap();
    let v = x("${ARCHAIC_DATA_VERSION}");
    let adl = g.join(format!("downloads/archaic-{v}"));
    fs::create_dir_all(&adl).unwrap();
    // An interrupted --archaic run: the sites file (another row, fetched
    // first), fetch's lock, a VCF's .part and a saved checksum.
    for f in [
        "1kg_sites.vcf.gz",
        ".1kg_sites.vcf.gz.lock",
        "Altai.chr1.vcf.gz.part",
        "1kg_sites.vcf.gz.checksum",
    ] {
        fs::write(adl.join(f), b"x").unwrap();
    }
    let check = || {
        plan_cli(
            &root,
            &["--for", "archaic", "--check"],
            std::slice::from_ref(&genome),
        )
    };
    let row = format!("archaic_vcf\tdownloads/archaic-{v}/\t1000000.0 TB\t");
    let (rc, out, err) = check();
    assert_eq!(rc, 1, "{out}{err}");
    assert!(out.contains(&format!("{row}missing\n")), "{out}");
    // A VCF of its own: still to download, less what it holds.
    fs::write(adl.join("Altai.chr1.vcf.gz"), b"data").unwrap();
    let (rc, out, err) = check();
    assert_eq!(rc, 1, "{out}{err}");
    assert!(
        out.contains(&format!("{row}missing, 4 bytes of it already here\n")),
        "{out}"
    );
    // The install done, stamped with today's pins: present.
    let stamp = x(
        "data=${ARCHAIC_DATA_VERSION} regions=all bcftools=${BCFTOOLS_IMAGE} picard=${PICARD_IMAGE} masks=m code=c",
    );
    let done = g.join(format!("reference/archaic/{v}/DONE"));
    fs::create_dir_all(done.parent().unwrap()).unwrap();
    fs::write(&done, format!("{stamp}\n")).unwrap();
    let (rc, out, err) = check();
    assert_eq!(rc, 0, "{out}{err}");
    assert!(out.contains(&format!("{row}present\n")), "{out}");
    // A stamp of an older bcftools: setup redoes it, so it is to download.
    fs::write(&done, stamp.replace("bcftools=", "bcftools=old-")).unwrap();
    let (rc, out, err) = check();
    assert_eq!(rc, 1, "{out}{err}");
    // Done, its VCFs pruned: nothing to download.
    fs::write(&done, format!("{stamp}\n")).unwrap();
    fs::remove_file(adl.join("Altai.chr1.vcf.gz")).unwrap();
    let (rc, out, err) = check();
    assert_eq!(rc, 0, "{out}{err}");
    assert!(
        out.contains(&format!("{row}installed, not kept\n")),
        "{out}"
    );
}

#[test]
fn a_row_of_many_installed_files_or_a_family_counts_each_file_it_lacks() {
    let root = real_checkout(None);
    let s = Sources::read(&root).unwrap();
    let g = tmp("markers");
    let c = ctx(&root, &[("GENOME_DIR", g.display().to_string())]);
    let x = |t: &str| c.expand(t, None).unwrap();
    // GENCODE: the genes table there, the exons BED not.
    let genes = g.join(x("reference/gencode.v${GENCODE_RELEASE}.basic.genes.gtf"));
    fs::create_dir_all(genes.parent().unwrap()).unwrap();
    fs::write(&genes, b"x").unwrap();
    // MoChA: 21 of 22 chromosomes cut, each stamped as setup stamps it.
    let panel = g.join(x("reference/mocha-${MOCHA_PANEL_VERSION}/panel"));
    fs::create_dir_all(&panel).unwrap();
    let stamp =
        x("manifest=${MOCHA_PANEL_MANIFEST_SHA256} regions=all bcftools=${MOCHA_IMAGE} cnp=1");
    for n in 1..=21 {
        fs::write(panel.join(format!("chr{n}.vcf.gz.tbi")), b"x").unwrap();
        fs::write(panel.join(format!("chr{n}.stamp")), format!("{stamp}\n")).unwrap();
    }
    let rows = og_setup::plan::chosen(&s, &["setup".to_owned(), "mocha-data".to_owned()]);
    let l = og_setup::plan::lines(&c, &s, &rows).unwrap();
    let get = |id: &str| l.iter().find(|x| x.id == id).unwrap().clone();
    let gc = get("gencode");
    assert_eq!((gc.present, gc.installed), (0, 0), "MP1: {gc:?}");
    let mp = get("mocha_panel");
    assert_eq!(
        (mp.files, mp.present, mp.installed),
        (22, 0, 21),
        "MP2: {mp:?}"
    );
    let mut out = vec![];
    og_setup::plan::write_plan(&c, std::slice::from_ref(&mp), &mut out).unwrap();
    let out = String::from_utf8(out).unwrap();
    let part = og_setup::plan::human(mp.bytes.unwrap() / 22);
    assert!(
        out.contains("\t0 of 22 files present, 21 installed and not kept\n")
            && out.contains(&format!("To download: 1 source, {part};")),
        "one part counted: {out}"
    );
}

#[test]
fn a_stamped_install_is_done_only_while_its_stamp_holds_todays_pins() {
    let root = real_checkout(None);
    let s = Sources::read(&root).unwrap();
    let g = tmp("stamps");
    let c = ctx(&root, &[("GENOME_DIR", g.display().to_string())]);
    let x = |t: &str| c.expand(t, None).unwrap();
    let put = |p: &str, text: &str| {
        let p = g.join(x(p));
        fs::create_dir_all(p.parent().unwrap()).unwrap();
        fs::write(p, text).unwrap();
    };
    let done = |id: &str, part: Option<&str>| c.installed(s.get(id).unwrap(), part).unwrap();
    // gnomAD-SV: AnnotSV's merge consumes the converted file; the stamp stays.
    let sv =
        "annotsv_annotations/Annotations_Human/SVincludedInFt/BenignSV/GRCh38/gnomad-sv.installed";
    assert!(!done("gnomad_sv", None));
    put(
        sv,
        &x("gnomAD-SV ${GNOMAD_SV_VERSION} sha256 ${GNOMAD_SV_SHA256}\n"),
    );
    assert!(done("gnomad_sv", None), "the stamp of today's release");
    put(sv, "gnomAD-SV 4.0 sha256 0000\n");
    assert!(
        !done("gnomad_sv", None),
        "an older release: setup adds the new one"
    );
    // MANE: the transcripts and their release stamp.
    put("annotsv_annotations/mane_select_refseq.txt", "NM_1.1\n");
    assert!(!done("mane_summary", None), "no stamp");
    put(
        "annotsv_annotations/mane_select_refseq.release",
        &x("MANE ${MANE_RELEASE} sha256 ${MANE_SUMMARY_SHA256}\n"),
    );
    assert!(done("mane_summary", None));
    // MoChA: a chromosome stamped with another image is cut again.
    let panel = "reference/mocha-${MOCHA_PANEL_VERSION}/panel";
    let stamp =
        x("manifest=${MOCHA_PANEL_MANIFEST_SHA256} regions=all bcftools=${MOCHA_IMAGE} cnp=1");
    put(&format!("{panel}/chr7.vcf.gz.tbi"), "x");
    put(&format!("{panel}/chr7.stamp"), &stamp);
    assert!(done("mocha_panel", Some("7")) && done("mocha_panel_index", Some("7")));
    put(
        &format!("{panel}/chr7.stamp"),
        &stamp.replace("bcftools=", "bcftools=x"),
    );
    assert!(!done("mocha_panel", Some("7")), "a new MOCHA_IMAGE");
    // The fine ancestry panel: the same, by its pins and bcftools image.
    let fine = "reference/ancestry_fine/${ANCESTRY_FINE_PANEL}/panel";
    let pstamp = x(
        "pins=${ANCESTRY_FINE_PANEL_SHA256} regions=all thin=1 keep=k bcftools=${BCFTOOLS_HTTPS_IMAGE}",
    );
    put(&format!("{fine}/chr3.vcf.gz.tbi"), "x");
    put(&format!("{fine}/chr3.stamp"), &pstamp);
    assert!(done("ancestry_fine_bcf", Some("3")));
    put(
        &format!("{fine}/chr3.stamp"),
        &pstamp.replace("pins=", "pins=0"),
    );
    assert!(!done("ancestry_fine_bcf", Some("3")), "new panel pins");
    // A word must be a whole word of the stamp, and a large file is no stamp.
    put(
        sv,
        &x("gnomAD-SV ${GNOMAD_SV_VERSION}x sha256 ${GNOMAD_SV_SHA256}\n"),
    );
    assert!(!done("gnomad_sv", None));
    put(
        sv,
        &format!(
            "{}{}",
            x("gnomAD-SV ${GNOMAD_SV_VERSION} sha256 ${GNOMAD_SV_SHA256}\n"),
            " ".repeat(70_000)
        ),
    );
    assert!(!done("gnomad_sv", None));
    // The table refuses an empty stamp word.
    let bad = format!("{TABLE_HEAD}a\tsetup\t1\tdownloads/a\t\thttps://h/a\t\t\t\t\ta.done#\n");
    assert!(
        Sources::parse(&bad)
            .unwrap_err()
            .contains("a stamp is FILE#WORD")
    );
}

#[test]
fn the_next_gencode_release_is_one_line_of_common_sh() {
    let now = read("scripts/lib/common.sh");
    let line = now
        .lines()
        .find(|l| l.starts_with("GENCODE_RELEASE="))
        .expect("common.sh sets GENCODE_RELEASE");
    let root = real_checkout(Some(("scripts/lib/common.sh", line, "GENCODE_RELEASE=99")));
    let s = Sources::read(&root).unwrap();
    let g = tmp("gencode99");
    let c = ctx(&root, &[("GENOME_DIR", g.display().to_string())]);
    let r = s.get("gencode").unwrap();
    assert_eq!(
        c.path(r, None).unwrap(),
        g.join("downloads/gencode.v99.basic.annotation.gtf.gz")
    );
    assert!(c.url(r, None).unwrap().contains("/release_99/gencode.v99."));
    assert!(c.value(&s, r, None).unwrap().contains("/release_99/"));
    // And og setup plan sees that file.
    fs::create_dir_all(g.join("downloads")).unwrap();
    fs::write(
        g.join("downloads/gencode.v99.basic.annotation.gtf.gz"),
        b"x",
    )
    .unwrap();
    let rows = og_setup::plan::chosen(&s, &["setup".to_owned()]);
    let l = og_setup::plan::lines(&c, &s, &rows).unwrap();
    let gl = l.iter().find(|x| x.id == "gencode").unwrap();
    assert_eq!(gl.present, 1, "{gl:?}");
}

#[test]
fn shell_settings_are_read_as_the_shell_sets_them() {
    let t = "A_B=12\nHLA=${HLA:-3.66.0}  # a comment\nQ=\"x y\"\nR=\"1.2\"\n  INDENTED=1\nlower=1\nS=$(date)\nT=${OTHER:-1}\nA_B=13\n";
    let m = shell_settings(t);
    assert_eq!(
        m.get("A_B"),
        Some(&Setting::Fixed("12".into())),
        "the first one"
    );
    assert_eq!(m.get("HLA"), Some(&Setting::EnvOr("3.66.0".into())));
    assert_eq!(m.get("R"), Some(&Setting::Fixed("1.2".into())));
    for k in ["Q", "INDENTED", "lower", "S", "T"] {
        assert!(!m.contains_key(k), "{k}");
    }
    let root = checkout("a\tsetup\t1\treference/${HLA}/${A_B}\t\t-\t\t\t\t\t-\n", "");
    fs::create_dir_all(root.join("scripts/lib")).unwrap();
    fs::write(root.join("scripts/lib/common.sh"), t).unwrap();
    let s = Sources::read(&root).unwrap();
    let r = s.get("a").unwrap();
    let c = ctx(&root, &[("GENOME_DIR", "/g".into())]);
    assert_eq!(
        c.path(r, None).unwrap(),
        PathBuf::from("/g/reference/3.66.0/12")
    );
    let o = ctx(
        &root,
        &[
            ("GENOME_DIR", "/g".into()),
            ("HLA", "3.70.0".into()),
            ("A_B", "99".into()),
        ],
    );
    assert_eq!(
        o.path(r, None).unwrap(),
        PathBuf::from("/g/reference/3.70.0/12"),
        "the environment overrides HLA=${{HLA:-...}}, never A_B=12"
    );
}

#[test]
fn a_refusal_names_the_row_of_its_file() {
    let root = checkout(
        "one\tsetup\t1\tdownloads/one.zip\t\thttps://h/one.zip\t\t\t\t\tone.done\n\
         fam\tsetup\t1\tdownloads/fam/{part}.gz\ta b\thttps://h/{part}.gz\t\t\t\t\tfam/{part}.done\n\
         dir\tsetup\t1\tdownloads/dir/\t\t-\t\t\t\t\tdir.done\n\
         inner\tsetup\t1\tdownloads/dir/inner.gz\t\t-\t\t\t\t\tdir.done\n",
        "",
    );
    let s = Sources::read(&root).unwrap();
    let g = tmp("sourcefor");
    let c = ctx(&root, &[("GENOME_DIR", g.display().to_string())]);
    let at = |p: &str| source_for(&c, &s, &g.join(p));
    assert_eq!(
        at("downloads/one.zip"),
        Some(("one".into(), Some("one".into())))
    );
    assert_eq!(at("downloads/fam/b.gz"), Some(("fam b".into(), None)));
    assert_eq!(at("downloads/dir/x/y.gz"), Some(("dir".into(), None)));
    assert_eq!(
        at("downloads/dir/inner.gz"),
        Some(("inner".into(), Some("inner".into()))),
        "a file row inside a folder row is its own"
    );
    assert_eq!(at("downloads/two.zip"), None);
    assert_eq!(at("downloads/fam/c.gz"), None);
}

#[test]
fn prune_removes_only_downloads() {
    let root = checkout(
        "raw\tsetup\t13\tdownloads/raw.tar.gz\t\thttps://h/raw.tar.gz\t\t\t\t\traw.done\n\
         fam\tsetup\t43\tdownloads/f/{part}.gz\t1..2\thttps://h/{part}.gz\t\t\t\t\tf{part}.done\n\
         inst\tsetup\t33\treference/inst.vcf.gz\t\thttps://h/inst.vcf.gz\t\t\t\t\t-\n",
        "",
    );
    let g = root.join("g");
    for f in [
        "downloads/raw.tar.gz",
        "downloads/raw.tar.gz.part",
        "downloads/f/2.gz",
        "downloads/mine.txt",
        "reference/inst.vcf.gz",
    ] {
        fs::create_dir_all(g.join(f).parent().unwrap()).unwrap();
        fs::write(g.join(f), b"data").unwrap();
    }
    let c = ctx(&root, &[("GENOME_DIR", g.display().to_string())]);
    let s = Sources::read(&root).unwrap();
    let rows = og_setup::prune::rows(&s, &[]).unwrap();
    let t = og_setup::prune::targets(&c, &rows).unwrap();
    let names: Vec<String> = t
        .iter()
        .map(|(p, _)| p.strip_prefix(&g).unwrap().display().to_string())
        .collect();
    assert_eq!(
        names,
        [
            "downloads/raw.tar.gz",
            "downloads/raw.tar.gz.part",
            "downloads/f/2.gz"
        ]
    );
    // A folder of the table moved away and linked back: its files are
    // outside the downloads folder, and nothing is removed.
    let away = tmp("away");
    fs::write(away.join("1.gz"), b"elsewhere").unwrap();
    fs::remove_dir_all(g.join("downloads/f")).unwrap();
    std::os::unix::fs::symlink(&away, g.join("downloads/f")).unwrap();
    let e = og_setup::prune::targets(&c, &rows).unwrap_err();
    assert!(e.contains("resolves outside the downloads folder"), "{e}");
    assert!(away.join("1.gz").exists());
    let named = og_setup::prune::rows(&s, &["inst".to_owned()]).unwrap_err();
    assert!(
        named.contains("prune removes only the downloads folder"),
        "{named}"
    );
}

// --- the repository ---------------------------------------------------------

fn repo() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../..")
        .canonicalize()
        .unwrap()
}

fn read(p: &str) -> String {
    fs::read_to_string(repo().join(p)).unwrap_or_else(|e| panic!("{p}: {e}"))
}

/// Pins of the table also written in one of FILES: (id, pin, file).
fn duplicated_pins(s: &Sources, files: &[(String, String)]) -> Vec<(String, String, String)> {
    let mut out = vec![];
    for (id, pin) in literal_pins(s) {
        for (name, text) in files {
            if text.contains(&pin) {
                out.push((id.clone(), pin.clone(), name.clone()));
            }
        }
    }
    out
}

/// versions.env and every script under scripts/ and scripts/lib/ (setup,
/// its library and the steps): setup.sh first, then common.sh, then the rest.
fn pin_files() -> Vec<(String, String)> {
    let mut names = vec![
        "versions.env".to_owned(),
        "scripts/setup.sh".to_owned(),
        "scripts/lib/common.sh".to_owned(),
    ];
    for dir in ["scripts", "scripts/lib"] {
        let mut more: Vec<String> = fs::read_dir(repo().join(dir))
            .unwrap()
            .flatten()
            .map(|e| format!("{dir}/{}", e.file_name().to_string_lossy()))
            .filter(|n| n.ends_with(".sh") && !names.contains(n))
            .collect();
        more.sort();
        names.extend(more);
    }
    names.into_iter().map(|f| (f.clone(), read(&f))).collect()
}

#[test]
fn every_pin_of_the_table_is_written_once() {
    let s = Sources::read(&repo()).unwrap();
    let dup = duplicated_pins(&s, &pin_files());
    assert!(
        dup.is_empty(),
        "pins the table holds, written again elsewhere: {dup:#?}"
    );
    // Every row expands: each ${NAME} is a pin of versions.env, each part fits.
    let c = ctx(&repo(), &[("GENOME_DIR", "/g".into())]);
    for r in s.rows.values() {
        let parts: Vec<Option<&str>> = if r.parts.is_empty() {
            vec![None]
        } else {
            r.parts.iter().map(|p| Some(p.as_str())).collect()
        };
        for p in parts {
            c.path(r, p).unwrap();
            c.url(r, p).unwrap();
            c.value(&s, r, p).unwrap();
        }
    }
}

#[test]
fn a_pin_planted_twice_is_caught() {
    let s = Sources::read(&repo()).unwrap();
    let (id, pin) = literal_pins(&s)
        .into_iter()
        .find(|(_, p)| !p.starts_with("http"))
        .expect("the table holds a checksum");
    let mut files = pin_files();
    files[1].1.push_str(&format!("\nPLANTED_SHA256={pin}\n"));
    let dup = duplicated_pins(&s, &files);
    assert_eq!(dup, [(id, pin, "scripts/setup.sh".to_owned())]);
}

/// Defaults of the table, `${NAME:-default}`, for a NAME that versions.env,
/// setup.sh or common.sh sets: (id, NAME, file). The default then lives twice.
fn repeated_defaults(s: &Sources, files: &[(String, String)]) -> Vec<(String, String, String)> {
    let mut out = vec![];
    for r in s.rows.values() {
        let fields = [&r.path, &r.url, &r.value]
            .into_iter()
            .chain(r.installed.iter())
            .cloned()
            .collect::<Vec<String>>()
            .join("\t");
        for (i, _) in fields.match_indices("${") {
            let inner = &fields[i + 2..];
            let Some((name, _)) = inner.split_once(":-") else {
                continue;
            };
            if name.contains('}') {
                continue;
            }
            for (f, text) in files {
                let sets = text.lines().any(|l| {
                    let l = l.trim_start();
                    let l = l
                        .strip_prefix("export ")
                        .or_else(|| l.strip_prefix("readonly "))
                        .unwrap_or(l);
                    l.strip_prefix(name).is_some_and(|x| x.starts_with('='))
                });
                if sets {
                    out.push((r.id.clone(), name.to_owned(), f.clone()));
                }
            }
        }
    }
    out
}

/// True when TEXT spells a URL whose start is PREFIX and whose end is TAIL:
/// on one line, or split over a variable, `base=PREFIX...` on one line and
/// `${base}...TAIL` (or `$base`) on another. A URL with no `${` (TAIL
/// empty) split that way: a variable set to the URL's start (a host and a
/// folder at least), then used before the rest of it.
fn spells_url(text: &str, url: &str, prefix: &str, tail: &str) -> bool {
    if text.lines().any(|l| {
        l.find(prefix)
            .is_some_and(|i| l[i + prefix.len()..].contains(tail))
    }) {
        return true;
    }
    let mut bases: Vec<(String, String)> = vec![];
    for l in text.lines() {
        for (k, v) in assignments(l) {
            let v = v.trim_matches(|c| c == '"' || c == '\'');
            if !tail.is_empty() && v.starts_with(prefix) {
                bases.push((k, tail.to_owned()));
            } else if tail.is_empty()
                && v.len() < url.len()
                && url.starts_with(v)
                && v.split_once("://")
                    .is_some_and(|(_, rest)| rest.contains('/'))
            {
                bases.push((k, url[v.len()..].to_owned()));
            }
        }
    }
    text.lines().any(|l| {
        bases.iter().any(|(k, rest)| {
            [format!("${{{k}}}"), format!("${k}")].iter().any(|use_| {
                l.match_indices(use_.as_str()).any(|(i, _)| {
                    let after = &l[i + use_.len()..];
                    let word_end = use_.ends_with('}')
                        || !after.starts_with(|c: char| c.is_ascii_alphanumeric() || c == '_');
                    word_end && after.contains(rest.as_str())
                })
            })
        })
    })
}

/// URL with a leading `${NAME}` (a pin) replaced by NAME's value in
/// versions.env (ENV); None when ENV does not set NAME.
fn literal_url(url: &str, env: &str) -> Option<String> {
    let Some(rest) = url.strip_prefix("${") else {
        return Some(url.to_owned());
    };
    let (name, after) = rest.split_once('}')?;
    let v = env.lines().find_map(|l| {
        let l = l.trim_start();
        let l = l.strip_prefix("export ").unwrap_or(l);
        l.strip_prefix(name)?.strip_prefix('=')
    })?;
    let v = match v.strip_prefix('"') {
        Some(q) => q.split('"').next()?,
        None => v.split_whitespace().next()?,
    };
    Some(format!("{v}{after}"))
}

/// URLs of the table that a script also spells (`spells_url`): the start
/// of a row's URL up to its first `${` or `{part}` (at least a host and a
/// folder) and, after it, the URL's end after its last `}` (the file name's
/// tail): (id, prefix, file). A row whose URL is a pin, `${NAME}...`, is
/// tested with NAME's value from versions.env. A row whose installer builds
/// its URL says `-`; every other script reads its URL from the table.
fn repeated_urls(s: &Sources, files: &[(String, String)]) -> Vec<(String, String, String)> {
    let env = files
        .iter()
        .find(|(f, _)| f == "versions.env")
        .map_or("", |(_, t)| t.as_str());
    let mut out = vec![];
    for r in s.rows.values() {
        if r.url == "-" {
            continue;
        }
        let Some(url) = literal_url(&r.url, env) else {
            continue;
        };
        let end = [url.find("${"), url.find("{part}")]
            .into_iter()
            .flatten()
            .min()
            .unwrap_or(url.len());
        let prefix = url[..end].trim_end_matches('/');
        let tail = if end < url.len() {
            url.rsplit('}').next().unwrap_or("")
        } else {
            ""
        };
        let host_and_folder = prefix
            .split_once("://")
            .is_some_and(|(_, rest)| rest.trim_end_matches('/').contains('/'));
        if !host_and_folder {
            continue;
        }
        for (f, text) in files {
            if *f != "versions.env" && spells_url(text, &url, prefix, tail) {
                out.push((r.id.clone(), prefix.to_owned(), f.clone()));
            }
        }
    }
    out
}

#[test]
fn no_default_or_url_of_the_table_is_written_again() {
    let s = Sources::read(&repo()).unwrap();
    let files = pin_files();
    let d = repeated_defaults(&s, &files);
    assert!(d.is_empty(), "defaults the table repeats: {d:#?}");
    let u = repeated_urls(&s, &files);
    assert!(
        u.is_empty(),
        "URLs of the table spelled again by an installer: {u:#?}"
    );
}

#[test]
fn a_planted_default_or_url_is_caught() {
    let t = read("assets/data-sources.tsv");
    let row = t.lines().find(|l| l.starts_with("gencode\t")).unwrap();
    // (a) The GENCODE release written again as the table's own default.
    let planted = t.replace(
        row,
        &row.replacen("${GENCODE_RELEASE}", "${GENCODE_RELEASE:-50}", 1),
    );
    let s = Sources::parse(&planted).unwrap();
    assert_eq!(
        repeated_defaults(&s, &pin_files()),
        [(
            "gencode".to_owned(),
            "GENCODE_RELEASE".to_owned(),
            "scripts/lib/common.sh".to_owned()
        )]
    );
    // (b) An installer that builds the VEP URL again, as main's common.sh did.
    let s = Sources::read(&repo()).unwrap();
    let mut files = pin_files();
    files[1].1.push_str(
        "\n  printf 'https://ftp.ensembl.org/pub/release-%s/variation/indexed_vep_cache/homo_sapiens_vep_%s_GRCh38.tar.gz' \"$1\" \"$1\"\n",
    );
    let u = repeated_urls(&s, &files);
    assert!(
        u.iter().any(|(id, p, f)| id == "vep_cache"
            && p == "https://ftp.ensembl.org/pub/release-"
            && f == "scripts/setup.sh"),
        "{u:#?}"
    );
    // PN2: the same URL split over a base variable, in a step script.
    let mut files = pin_files();
    let step = files
        .iter()
        .position(|(f, _)| f == "scripts/25-prs.sh")
        .expect("pin_files reads the step scripts");
    files[step].1.push_str(
        "\n  base=\"https://ftp.ensembl.org/pub/release-${rel}/variation/indexed_vep_cache\"\n  url=\"${base}/homo_sapiens_vep_${rel}_GRCh38.tar.gz\"\n",
    );
    let u = repeated_urls(&s, &files);
    assert!(
        u.iter()
            .any(|(id, _, f)| id == "vep_cache" && f == "scripts/25-prs.sh"),
        "PN2 {u:#?}"
    );
    // And a URL with no pin in it, split the same way.
    let mut files = pin_files();
    files[step].1.push_str(
        "\n  z=https://zenodo.org/records/7072866/files\n  fetch \"$z/revel-v1.3_all_chromosomes.zip\" \"$out\"\n",
    );
    let u = repeated_urls(&s, &files);
    assert!(
        u.iter()
            .any(|(id, _, f)| id == "revel" && f == "scripts/25-prs.sh"),
        "a literal URL split {u:#?}"
    );
    // U3: a row whose URL is a pin, its value from versions.env written again.
    let mito = literal_url("${GNOMAD_MITO_URL}", &read("versions.env")).unwrap();
    assert!(mito.starts_with("https://"), "{mito}");
    let mut files = pin_files();
    let cpsr = files
        .iter()
        .position(|(f, _)| f == "scripts/17-cpsr.sh")
        .expect("pin_files reads the step scripts");
    files[cpsr]
        .1
        .push_str(&format!("\n  fetch \"{mito}\" \"$out\"\n"));
    let u = repeated_urls(&s, &files);
    assert!(
        u.iter()
            .any(|(id, _, f)| id == "gnomad_mito" && f == "scripts/17-cpsr.sh"),
        "U3 {u:#?}"
    );
}

/// The shell functions of TEXT (and its top level, as one more block):
/// (name, body lines with their line numbers).
fn blocks(text: &str) -> Vec<(String, Vec<(usize, String)>)> {
    let mut out: Vec<(String, Vec<(usize, String)>)> = vec![("<top level>".into(), vec![])];
    let mut cur: Option<usize> = None;
    for (n, l) in text.lines().enumerate() {
        if cur.is_none() {
            if let Some(name) = l
                .strip_suffix("() {")
                .filter(|x| !x.starts_with(' ') && !x.contains(' '))
            {
                out.push((name.to_owned(), vec![]));
                cur = Some(out.len() - 1);
                continue;
            }
            out[0].1.push((n + 1, l.to_owned()));
        } else if l == "}" {
            cur = None;
        } else {
            out[cur.unwrap()].1.push((n + 1, l.to_owned()));
        }
    }
    out
}

/// `$x` written as `${x}`, `${x:?msg}`, `${x:-y}` and `${x?}` as `${x}` (the
/// same path when x is set), and quotes dropped, so two spellings of a path
/// compare.
fn norm(s: &str) -> String {
    let s = s.trim_matches(|c| c == '"' || c == '\'');
    let mut out = String::new();
    let b: Vec<char> = s.chars().collect();
    let name_at = |j: usize| {
        let mut k = j;
        while k < b.len() && (b[k].is_ascii_alphanumeric() || b[k] == '_') {
            k += 1;
        }
        k
    };
    // `${NAME` followed by `?`, `:?` or `:-`: where NAME ends.
    let modified = |i: usize| {
        if b[i] != '$' || b.get(i + 1) != Some(&'{') {
            return None;
        }
        let k = name_at(i + 2);
        let m = b.get(k) == Some(&'?')
            || (b.get(k) == Some(&':') && matches!(b.get(k + 1), Some('?' | '-')));
        (k > i + 2 && m).then_some(k)
    };
    let mut i = 0;
    while i < b.len() {
        if b[i] == '$' && i + 1 < b.len() && (b[i + 1].is_ascii_alphabetic() || b[i + 1] == '_') {
            let j = name_at(i + 1);
            out.push_str(&format!("${{{}}}", b[i + 1..j].iter().collect::<String>()));
            i = j;
        } else if let Some(k) = modified(i) {
            // The modifier's word runs to its closing brace.
            let mut depth = 1;
            let mut j = k;
            while j < b.len() && depth > 0 {
                match b[j] {
                    '{' => depth += 1,
                    '}' => depth -= 1,
                    _ => {}
                }
                j += 1;
            }
            out.push_str(&format!("${{{}}}", b[i + 2..k].iter().collect::<String>()));
            i = j;
        } else {
            if b[i] != '"' {
                out.push(b[i]);
            }
            i += 1;
        }
    }
    out
}

/// The shape of a normalised path: its leading `${NAME}` kept, every later
/// `${...}` or `$(...)` a `*`, so `${dir}/d/${g}.chr${c}.gz` and
/// `${dir}/d/${v%%:*}.chr${c}.gz` compare equal.
fn shape(p: &str) -> String {
    let lead = if p.starts_with("${") {
        p.find('}').map_or(0, |i| i + 1)
    } else {
        0
    };
    let (head, rest) = p.split_at(lead);
    let mut out = head.to_owned();
    let b: Vec<char> = rest.chars().collect();
    let mut i = 0;
    while i < b.len() {
        if b[i] == '$' && i + 1 < b.len() && (b[i + 1] == '{' || b[i + 1] == '(') {
            let (open, close) = if b[i + 1] == '{' {
                ('{', '}')
            } else {
                ('(', ')')
            };
            let mut depth = 0;
            let mut j = i + 1;
            while j < b.len() {
                if b[j] == open {
                    depth += 1;
                } else if b[j] == close {
                    depth -= 1;
                    if depth == 0 {
                        break;
                    }
                }
                j += 1;
            }
            out.push('*');
            i = j + 1;
        } else {
            out.push(b[i]);
            i += 1;
        }
    }
    out
}

/// The words of a shell line, quotes, `${...}` and `$(...)` kept inside a word.
fn words(l: &str) -> Vec<String> {
    let mut out = vec![];
    let mut cur = String::new();
    let (mut q, mut depth, mut prev) = (None::<char>, 0i32, ' ');
    for c in l.chars() {
        match (q, c) {
            (None, '"' | '\'') if depth == 0 => {
                q = Some(c);
                cur.push(c);
            }
            (Some(x), _) if x == c && depth == 0 => {
                q = None;
                cur.push(c);
            }
            (_, '(') => {
                depth += 1;
                cur.push(c);
            }
            (_, '{') if prev == '$' || depth > 0 => {
                depth += 1;
                cur.push(c);
            }
            (_, ')' | '}') if depth > 0 => {
                depth -= 1;
                cur.push(c);
            }
            (None, ' ' | '\t' | ';' | '|' | '&' | '{' | '}') if depth == 0 => {
                if !cur.is_empty() {
                    out.push(std::mem::take(&mut cur));
                }
            }
            _ => cur.push(c),
        }
        prev = c;
    }
    if !cur.is_empty() {
        out.push(cur);
    }
    out
}

/// The `NAME=VALUE` words of a line (after `local`, in a list, after `||`).
fn assignments(l: &str) -> Vec<(String, String)> {
    words(l)
        .into_iter()
        .filter_map(|w| {
            let (k, v) = w.split_once('=')?;
            let name = !k.is_empty()
                && k.starts_with(|c: char| c.is_ascii_alphabetic() || c == '_')
                && k.chars().all(|c| c.is_ascii_alphanumeric() || c == '_');
            name.then(|| (k.to_owned(), v.to_owned()))
        })
        .collect()
}

/// Every place in TEXT that deletes a raw download: in one function, an
/// `rm` (or a `find ... -delete`/`-exec rm`) whose argument is the
/// destination of a `fetch` there, a folder holding it, or a path in the
/// downloads folder (a variable set from `$(dl ...)`, `source_of` or
/// OG_DOWNLOADS anywhere on a line, or from such a variable). Paths compare
/// by their shape: `${dir}/x/${g}.gz` is `${dir}/x/${v%%:*}.gz`. A fetch's
/// own `.part`, which it resumes, is not a download. `done_with` is the one
/// place that may: it deletes a download only when the user asked for it
/// (`SETUP_KEEP_DOWNLOADS=false`). A line that says `# not kept:` and why
/// deletes no raw download (a staging copy, an unpacked tree, a file no
/// source of the table names).
fn deleted_downloads(name: &str, text: &str) -> Vec<String> {
    let mut out = vec![];
    for (func, lines) in blocks(text) {
        if func == "done_with" {
            continue;
        }
        let mut sets: Vec<(String, String)> = vec![];
        let mut dests: Vec<String> = vec![];
        for (_, l) in &lines {
            if l.trim_start().starts_with('#') {
                continue;
            }
            let w = words(l);
            for (i, x) in w.iter().enumerate() {
                if x == "fetch" && i + 2 < w.len() && !w[i + 1].starts_with('-') {
                    dests.push(norm(&w[i + 2]));
                }
            }
            for (k, v) in assignments(l) {
                sets.push((format!("${{{k}}}"), v));
            }
        }
        // Downloads variables, and what they lead to.
        let mut dl_vars: Vec<String> = vec![];
        for _ in 0..4 {
            for (k, v) in &sets {
                let n = norm(v);
                let from_dl = v.contains("$(dl ")
                    || v.contains("$(source_of ")
                    || v.contains("OG_DOWNLOADS")
                    || dl_vars
                        .iter()
                        .any(|d| n == *d || n.starts_with(&format!("{d}/")));
                if from_dl && !dl_vars.contains(k) {
                    dl_vars.push(k.clone());
                }
            }
            // A fetch into `$f` fetches into what f was set to.
            for (k, v) in &sets {
                if dests.contains(k) {
                    let n = norm(v);
                    if !dests.contains(&n) {
                        dests.push(n);
                    }
                }
            }
        }
        let dests: Vec<String> = dests.iter().map(|d| shape(d)).collect();
        for (n, l) in &lines {
            if l.trim_start().starts_with('#') || l.contains("# not kept:") {
                continue;
            }
            let w = words(l);
            // `find PATH... -delete` or `-exec rm`: the paths before its first
            // test or action; else the arguments of `rm`.
            let found = w
                .iter()
                .position(|x| x == "find")
                .filter(|_| w.iter().any(|x| x == "-delete") || l.contains("-exec rm"));
            let args: Vec<&String> = if let Some(f) = found {
                w[f + 1..]
                    .iter()
                    .take_while(|a| !a.starts_with(['-', '(', '!']))
                    .collect()
            } else if let Some(at) = w.iter().position(|x| x == "rm") {
                w[at + 1..].iter().filter(|a| !a.starts_with('-')).collect()
            } else {
                continue;
            };
            for a in args {
                let a = shape(&norm(a));
                if a.ends_with(".part") {
                    continue;
                }
                // A holds D: A is D, a file inside D, or a folder above D.
                let holds = |d: &str| {
                    a == d
                        || a.starts_with(&format!("{d}/"))
                        || d.strip_prefix(a.as_str())
                            .is_some_and(|r| r.starts_with('/') || r.starts_with('*'))
                };
                let hit = a.contains("OG_DOWNLOADS")
                    || a.contains("$(dl ")
                    || dl_vars.iter().any(|v| holds(v))
                    || dests.iter().any(|d| holds(d));
                if hit {
                    out.push(format!("{name}:{n} ({func}): {}", l.trim()));
                    break;
                }
            }
        }
    }
    out
}

/// The scripts the delete and stream guards read: setup.sh, common.sh and
/// every other scripts/*.sh that calls fetch or done_with.
fn guarded_scripts() -> Vec<(String, String)> {
    let mut names = vec![
        "scripts/setup.sh".to_owned(),
        "scripts/lib/common.sh".to_owned(),
    ];
    let mut more: Vec<String> = std::fs::read_dir(repo().join("scripts"))
        .unwrap()
        .filter_map(|e| {
            let p = e.ok()?.path();
            let n = p.file_name()?.to_str()?.to_owned();
            n.ends_with(".sh").then(|| format!("scripts/{n}"))
        })
        .filter(|n| !names.contains(n) && calls_fetch(&read(n)))
        .collect();
    more.sort();
    names.extend(more);
    names
        .into_iter()
        .map(|n| {
            let t = read(&n);
            (n, t)
        })
        .collect()
}

/// True when TEXT calls fetch or done_with on a path (not in a comment).
fn calls_fetch(text: &str) -> bool {
    text.lines()
        .filter(|l| !l.trim_start().starts_with('#'))
        .any(|l| {
            words(l).windows(2).any(|p| {
                (p[0] == "fetch" || p[0] == "done_with")
                    && (p[1].starts_with('"') || p[1].starts_with('$'))
            })
        })
}

fn delete_guard(files: &[(String, String)]) -> Vec<String> {
    files
        .iter()
        .flat_map(|(n, t)| deleted_downloads(n, t))
        .collect()
}

fn stream_guard(files: &[(String, String)]) -> Vec<String> {
    files
        .iter()
        .flat_map(|(n, t)| unsaved_streams(n, t))
        .collect()
}

/// FILES with LINE planted right after the pgsc_calc fetch of 25-prs.sh,
/// which the guards must read.
fn planted_in_prs(line: &str) -> Vec<(String, String)> {
    let mut files = guarded_scripts();
    for n in ["scripts/25-prs.sh", "scripts/04b-gridss.sh"] {
        assert!(
            files.iter().any(|(f, _)| f == n),
            "{n} calls fetch, so the guards must read it"
        );
    }
    let at = "  fetch \"$PGSC_URL\" \"$PGSC_TGZ\" sha256 \"$PGSC_CALC_SHA256\"\n";
    let prs = files
        .iter_mut()
        .find(|(f, _)| f == "scripts/25-prs.sh")
        .unwrap();
    assert!(prs.1.contains(at), "25-prs.sh no longer has: {at}");
    prs.1 = prs.1.replacen(at, &format!("{at}{line}\n"), 1);
    files
}

#[test]
fn no_installer_deletes_a_raw_download() {
    let bad = delete_guard(&guarded_scripts());
    assert!(
        bad.is_empty(),
        "raw downloads deleted outside done_with:\n{}",
        bad.join("\n")
    );
}

#[test]
fn a_planted_delete_of_a_download_is_caught() {
    let planted = "install_x() {\n  local tgz\n  tgz=$(dl x.tar.gz)\n  fetch \"$url\" \"$tgz\" sha256 \"$k\"\n  tar -xzf \"$tgz\" -C \"$d\"\n  rm -f \"$tgz\"\n}\n\
                   install_y() {\n  fetch \"$url\" \"${dir}/raw/${f}${e}\" md5 \"$m\"\n  rm -f \"${dir}/raw/${f}\" \"${dir}/raw/${f}.tbi\"\n}\n\
                   install_z() {\n  fetch \"$u\" \"${dir}/download/${g}.vcf.gz\" md5 \"$m\"\n  rm -rf \"${dir}/download\" \"${dir}/work\"\n}\n\
                   install_ok() {\n  fetch \"$u\" \"$dest\" md5 \"$m\" || { rm -f \"${dest}.part\"; return 1; }\n  done_with \"$dest\"\n  rm -f \"${dest}.tmp\" \"${destdir}\"\n}\n\
                   install_marked() {\n  tgz=$(dl x.tar.gz)\n  rm -rf \"$tgz\"  # not kept: a staging copy\n}\n";
    let bad = deleted_downloads("planted.sh", planted);
    assert_eq!(bad.len(), 3, "{bad:#?}");
    assert!(
        bad[0].contains("install_x")
            && bad[1].contains("install_y")
            && bad[2].contains("install_z"),
        "{bad:#?}"
    );
    // P1: the AADR files, a downloads folder set after `||`, deleted at once.
    let p1 = "install_ancient_panel() {\n  local raw f src=${ANCIENT_SOURCE_DIR:-}\n  raw=\"${base}/raw_${AADR_VERSION}\"\n  [ -n \"$src\" ] || raw=$(dl \"aadr-${AADR_VERSION}\")\n  rm -rf \"$raw\"\n}\n";
    assert_eq!(deleted_downloads("p1.sh", p1).len(), 1, "P1");
    // main's AADR: the fetch into $f, f set from $raw, then raw deleted.
    let aadr = "install_ancient_panel() {\n  raw=\"${base}/raw_${AADR_VERSION}\"\n  for e in anno ind; do\n    f=\"${raw}/aadr.${e}\"\n    fetch \"$url\" \"$f\" md5 \"$md5\"\n  done\n  rm -rf \"$raw\"\n}\n";
    assert_eq!(deleted_downloads("aadr.sh", aadr).len(), 1, "main's AADR");
    // main's archaic VCFs: the same path spelt with another variable.
    let archaic = "install_archaic_data() {\n  for v in $genomes; do\n    g=${v%%:*}\n    fetch \"${ARCHAIC_URL}/${g}/chr${c}.vcf.gz\" \"${dir}/download/${g}.chr${c}.vcf.gz\" md5 \"file://${dir}/md5/${g}.txt\"\n  done\n  for v in $genomes; do rm -f \"${dir}/download/${v%%:*}.chr${c}.vcf.gz\"; done\n}\n";
    assert_eq!(
        deleted_downloads("archaic.sh", archaic).len(),
        1,
        "main's archaic VCFs"
    );
    // G1: find names the download before -exec rm (the rm word is no path).
    let g1 = "install_pgsc_calc() {\n  url=$(source_of pgsc_calc url) && tgz=$(source_of pgsc_calc path) || return 1\n  fetch \"$url\" \"$tgz\" sha256 \"$K\" || return 1\n  find \"$tgz\" -exec rm -f {} +\n}\n";
    assert_eq!(deleted_downloads("g1.sh", g1).len(), 1, "G1 find -exec rm");
    let g1d = g1.replace("-exec rm -f {} +", "-maxdepth 0 -delete");
    assert_eq!(deleted_downloads("g1.sh", &g1d).len(), 1, "G1 find -delete");
    let unpacked = "install_x() {\n  tgz=$(dl x.tar.gz)\n  fetch \"$u\" \"$tgz\" md5 \"$m\"\n  find \"${tree}.part\" -name '*.tmp' -exec rm -f {} +\n}\n";
    assert!(deleted_downloads("tree.sh", unpacked).is_empty());
    // G1b: `${tgz:?}` is `${tgz}`, and so are `${tgz:-x}` and `${tgz?}`.
    for planted in ["\"${tgz:?}\"", "\"${tgz:-/none}\"", "\"${tgz?unset}\""] {
        let g1b = g1.replace(
            "find \"$tgz\" -exec rm -f {} +",
            &format!("rm -f {planted}"),
        );
        assert_eq!(deleted_downloads("g1b.sh", &g1b).len(), 1, "G1b {planted}");
    }
    // The staged copy, marked, deletes no download.
    let staged = "install_ancient_panel() {\n  [ -n \"$src\" ] || raw=$(dl \"aadr-${AADR_VERSION}\")\n  if [ -n \"$src\" ]; then\n    rm -rf \"$raw\"  # not kept: the staged copy\n  else\n    done_with \"$raw\"\n  fi\n}\n";
    assert!(deleted_downloads("staged.sh", staged).is_empty());
    // D8: a step script that fetches is read too: pgsc_calc deleted in 25-prs.sh.
    let bad = delete_guard(&planted_in_prs("  rm -f \"$PGSC_TGZ\""));
    assert!(
        bad.iter().any(|b| b.starts_with("scripts/25-prs.sh:")),
        "D8 {bad:#?}"
    );
}

/// Why a network read may stream: `slice` (only the regions a user chose,
/// kept), or a read of no source's data.
const NET_REASONS: [&str; 6] = [
    "slice",
    "pip",
    "git clone",
    "plugin install",
    "sample header",
    "checksum",
];

/// True when the comment lines right above line I of LINES (the block that
/// ends on the line before it) hold `# net: REASON`, REASON one of
/// NET_REASONS: each read gives its own reason.
fn says_why(lines: &[&str], i: usize) -> bool {
    lines[..i]
        .iter()
        .rev()
        .take_while(|a| a.trim_start().starts_with('#'))
        .any(|a| {
            a.trim_start().strip_prefix("# net:").is_some_and(|r| {
                let r = r.trim_start();
                NET_REASONS.iter().any(|x| {
                    r.strip_prefix(x)
                        .is_some_and(|t| t.is_empty() || t.starts_with([',', ' ', ':']))
                })
            })
        })
}

/// Where L calls a body read (`_get`, `og setup fetch --get`).
fn get_calls(l: &str) -> Vec<usize> {
    l.match_indices("_get ")
        .map(|(i, _)| i)
        .filter(|&i| {
            i == 0
                || !l[..i]
                    .chars()
                    .next_back()
                    .is_some_and(|c| c.is_ascii_alphanumeric())
        })
        .chain(l.match_indices("fetch --get ").map(|(i, _)| i))
        .collect()
}

/// True when L streams a body read into another command: to stdout (`-`,
/// `/dev/stdout`, `/dev/fd/N`) through a pipe or into an output process
/// substitution `>(...)`, or inside a process substitution `<(...)` or `>(...)`.
fn piped_get(l: &str) -> bool {
    let calls = get_calls(l);
    let substituted = l
        .match_indices("<(")
        .chain(l.match_indices(">("))
        .any(|(s, _)| {
            let end = l[s..].find(')').map_or(l.len(), |e| s + e);
            calls.iter().any(|&c| c > s && c < end)
        });
    let Some(&at) = calls.iter().min() else {
        return false;
    };
    let rest = &l[at..];
    let to_stdout = |w: &String| {
        w == "-"
            || w.trim_matches('"') == "/dev/stdout"
            || w.trim_matches('"')
                .strip_prefix("/dev/fd/")
                .is_some_and(|n| !n.is_empty() && n.bytes().all(|c| c.is_ascii_digit()))
    };
    let into = [rest.find('|'), rest.find(">(")]
        .into_iter()
        .flatten()
        .min();
    substituted || into.is_some_and(|p| words(&rest[..p]).iter().any(to_stdout))
}

/// Every network read in TEXT that saves nothing and has no `# net:` comment
/// in the comment lines right above it saying why (NET_REASONS): a container's
/// network read (`run_in --net`), a body read piped into a command, and curl
/// or wget writing a body anywhere but through fetch.
fn unsaved_streams(name: &str, text: &str) -> Vec<String> {
    let lines: Vec<&str> = text.lines().collect();
    let mut out = vec![];
    for (i, l) in lines.iter().enumerate() {
        let t = l.trim_start();
        if t.starts_with('#') {
            continue;
        }
        let w = words(l);
        let net = w.iter().any(|x| x == "run_in")
            && w.iter().any(|x| x == "--net" || x.contains("net[@]"));
        if (net || piped_get(l)) && !says_why(&lines, i) {
            out.push(format!("{name}:{}: {}", i + 1, l.trim()));
            continue;
        }
        let tool = w.iter().position(|x| x == "curl" || x == "wget");
        if let Some(p) = tool {
            let head_only = w[p + 1..]
                .iter()
                .any(|x| x.starts_with("-") && !x.starts_with("--") && x.contains('I'));
            let message = w[..p]
                .iter()
                .any(|x| x == "echo" || x == "printf" || x == "command");
            if !head_only && !message {
                out.push(format!("{name}:{}: {}", i + 1, l.trim()));
            }
        }
    }
    out
}

#[test]
fn no_source_is_streamed_without_being_saved() {
    let bad = stream_guard(&guarded_scripts());
    assert!(
        bad.is_empty(),
        "network reads that save nothing, or say nothing:\n{}",
        bad.join("\n")
    );
}

#[test]
fn a_planted_stream_is_caught() {
    let planted = "x() {\n  if ! run_in --net --rw \"$dir\" \"$IMG\" bcftools view \"$URL\" > out\n}\n\
                   y() {\n  # net: slice, the regions of LOCI only\n  run_in --net --rw \"$dir\" \"$IMG\" tabix \"$URL\" chr1\n}\n\
                   z() {\n  curl -fsSL \"$URL\" | gzip -dc > out\n  curl -fsI \"$URL\"\n  command -v curl\n}\n";
    let bad = unsaved_streams("planted.sh", planted);
    assert_eq!(bad.len(), 2, "{bad:#?}");
    assert!(bad[0].contains(":2:") && bad[1].contains(":9:"), "{bad:#?}");
    // P3: a source unpacked straight from the network.
    let p3 = "install_x() {\n  _get \"$URL\" - | tar -xzf - -C \"$d\"\n}\n";
    assert_eq!(unsaved_streams("p3.sh", p3).len(), 1, "P3 _get");
    for get in [
        "  og_setup fetch --get \"$URL\" - | tar -xzf -\n",
        "  fetch --get \"$URL\" - 2>/dev/null | gzip -dc > out\n",
        "  md5=$(_get \"${url}/md5.txt\" - | awk '{print $1}')\n",
    ] {
        assert_eq!(unsaved_streams("p3.sh", get).len(), 1, "{get}");
    }
    let saved = "_get() {\n  og_setup fetch --get \"$@\"\n}\n  _get \"$URL\" \"$out\"\n  # net: checksum, the md5 of the release\n  md5=$(_get \"${url}/md5.txt\" - | awk '{print $1}')\n";
    assert!(unsaved_streams("ok.sh", saved).is_empty());
    // P4: a whole panel read with a comment that gives no reason.
    let p4 = "y() {\n  # net: the panel\n  run_in --net \"$IMG\" bcftools view \"$URL\" > out\n}\n";
    assert_eq!(unsaved_streams("p4.sh", p4).len(), 1, "P4");
    let header = "y() {\n  # net: sample header, the sample list of the panel\n  run_in --net \"$IMG\" bcftools query -l \"$URL\"\n}\n";
    assert!(unsaved_streams("header.sh", header).is_empty());
    // G2: stdout spelt as a device.
    for out in ["/dev/stdout", "\"/dev/stdout\"", "/dev/fd/1"] {
        let g2 = format!("install_x() {{\n  _get \"$url\" {out} | tar -xzf - -C \"$dir\"\n}}\n");
        assert_eq!(unsaved_streams("g2.sh", &g2).len(), 1, "G2 {out}");
    }
    // G2b: a body read in a process substitution.
    for get in ["_get \"$url\" -", "og_setup fetch --get \"$url\" -"] {
        let g2b = format!("install_x() {{\n  tar -xzf <({get}) -C \"$dir\"\n}}\n");
        assert_eq!(unsaved_streams("g2b.sh", &g2b).len(), 1, "G2b {get}");
    }
    // G3b: one reason covers the read right below it, not the next one.
    let g3b = "y() {\n  # net: slice, the regions of LOCI only\n  run_in --net \"$IMG\" tabix \"$URL\" chr1\n  run_in --net \"$IMG\" bcftools view \"$URL\" > out\n}\n";
    let bad = unsaved_streams("g3b.sh", g3b);
    assert_eq!(bad.len(), 1, "G3b {bad:#?}");
    assert!(bad[0].contains(":4:"), "G3b {bad:#?}");
    let block = "y() {\n  # net: slice, the regions of LOCI only\n  # (htslib keeps the index here)\n  run_in --net \"$IMG\" tabix \"$URL\" chr1\n}\n";
    assert!(
        unsaved_streams("block.sh", block).is_empty(),
        "a reason in the comment block right above"
    );
    // S6: a step script that fetches is read too: pgsc_calc streamed in 25-prs.sh.
    let bad = stream_guard(&planted_in_prs(
        "  curl -fsSL \"$PGSC_URL\" | tar -xzf - -C \"${PGSC_CALC_DIR}.part\"",
    ));
    assert!(
        bad.iter().any(|b| b.starts_with("scripts/25-prs.sh:")),
        "S6 {bad:#?}"
    );
    // S3: a body read into an output process substitution.
    let s3 = "install_x() {\n  _get \"$url\" - > >(tar -xzf - -C \"$dir\")\n}\n";
    assert_eq!(unsaved_streams("s3.sh", s3).len(), 1, "S3");
    // S5: a quoted word before curl makes no message; echo, printf or command does.
    let s5 = "install_x() {\n  docker run --rm \"$IMG\" curl -fsSL \"$url\" | tar -xzf - -C \"$dir\"\n}\n";
    assert_eq!(unsaved_streams("s5.sh", s5).len(), 1, "S5");
    let said = "x() {\n  printf '%s\\n' curl wget\n  echo curl is missing\n}\n";
    assert!(unsaved_streams("said.sh", said).is_empty(), "a message");
}

#[test]
fn common_sh_and_og_agree_on_the_downloads_folder() {
    // One rule: common.sh exports OG_DOWNLOADS (OG_DOWNLOADS, or
    // GENOME_DIR/downloads), which og reads.
    let c = read("scripts/lib/common.sh");
    assert!(
        c.contains("export OG_DOWNLOADS=\"${OG_DOWNLOADS:-${GENOME_DIR}/downloads}\""),
        "common.sh no longer sets the downloads folder as og does"
    );
    let root = checkout(
        "a\tsetup\t1\tdownloads/a\t\thttps://h/a\t\t\t\t\ta.done\n",
        "",
    );
    let s = Sources::read(&root).unwrap();
    let g = ctx(&root, &[("GENOME_DIR", "/g".into())]);
    assert_eq!(
        g.path(s.get("a").unwrap(), None).unwrap(),
        PathBuf::from("/g/downloads/a")
    );
}
