//! `og setup fetch`, `plan` and `prune` against a local HTTP server, and
//! the repository checks that keep setup's downloads in one place: every
//! pin of `assets/data-sources.tsv` written once, no installer deleting a
//! raw download, no source streamed without being saved.
//!
//! The server counts requests per path and can fail the first N of them,
//! or cut a body short, so retries, resumes and reuse are observed, not
//! assumed. `curl` does the transfers, as in setup.

use og_setup::fetch::{
    self, Bad, Ctx, Done, Fetch, Kind, Knobs, OnWrong, Sources, digest, literal_pins,
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
    let path = line.split_whitespace().nth(1).unwrap_or("/").to_owned();
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
                let from = range
                    .strip_prefix("bytes=")
                    .and_then(|x| x.trim_end_matches('-').parse::<usize>().ok())
                    .unwrap_or(0);
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
    }
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
    assert!(log.contains("og setup prune mine"), "{log}");
    assert!(s.requests("/mine.bin").is_empty());
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

const TABLE_HEAD: &str = "id\tinstallers\tsteps\tpath\tparts\turl\tcheck\tvalue\tbytes\tlicence\n";

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
        "vep\tsetup,vep-cache\t13\tdownloads/vep_${VEP}.tar.gz\t\thttps://h/release-${VEP}/v.tar.gz\tsum\thttps://h/release-${VEP}/CHECKSUMS\t26000000000\tEnsembl\n\
         ref\tsetup\t02\treference/ref.fna.gz\t\t${REF_URL:-https://h/ref.fna.gz}\tmd5\thttps://h/md5\t900000000\tNCBI\n\
         chr\tarchaic\t43\tdownloads/arch/{part}.vcf.gz\t21..22\thttps://h/chr{part}.vcf.gz\tmd5\t@list\t\tMPI\n\
         list\tarchaic\t43\tdownloads/arch/md5.txt\t\thttps://h/md5.txt\tsha256\tabc\t10\tMPI\n",
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
        "x\tsetup\t1\tdownloads/x\t\thttps://h/${NOPE}\t\t\t\t\n",
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
    let r = "a\tsetup\t1\tdownloads/a\t\thttps://h/a\t\t\t\t\n";
    let again = "a\tsetup\t1\tdownloads/a2\t\thttps://h/a2\t\t\t\t\n";
    assert!(
        Sources::parse(&format!("{TABLE_HEAD}{r}{again}"))
            .unwrap_err()
            .contains("twice")
    );
    let p = "b\tsetup\t1\tdownloads/a\t\thttps://h/b\t\t\t\t\n";
    assert!(
        Sources::parse(&format!("{TABLE_HEAD}{r}{p}"))
            .unwrap_err()
            .contains("also a's")
    );
    let k = "c\tsetup\t1\tdownloads/c\t\thttps://h/c\tcrc64\t\t\t\n";
    assert!(
        Sources::parse(&format!("{TABLE_HEAD}{k}"))
            .unwrap_err()
            .contains("unknown checksum type")
    );
    let up = "d\tsetup\t1\t../d\t\thttps://h/d\t\t\t\t\n";
    assert!(Sources::parse(&format!("{TABLE_HEAD}{up}")).is_err());
}

#[test]
fn plan_lists_rows_sizes_state_kept_raw_size_and_free_disk() {
    let root = checkout(
        "big\tsetup\t13\tdownloads/big.tar.gz\t\thttps://h/big.tar.gz\t\t\t26000000000\tEnsembl\n\
         small\tsetup\t33\treference/small.vcf.gz\t\thttps://h/small.vcf.gz\t\t\t9800000\tMIT\n\
         fam\tsetup\t43\tdownloads/f/{part}.gz\t1..4\thttps://h/{part}.gz\t\t\t4000\tMPI\n\
         other\tarchaic\t43\tdownloads/o.gz\t\thttps://h/o.gz\t\t\t1\tMPI\n",
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
    let lines = og_setup::plan::lines(&c, &rows).unwrap();
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
fn prune_removes_only_downloads() {
    let root = checkout(
        "raw\tsetup\t13\tdownloads/raw.tar.gz\t\thttps://h/raw.tar.gz\t\t\t\t\n\
         fam\tsetup\t43\tdownloads/f/{part}.gz\t1..2\thttps://h/{part}.gz\t\t\t\t\n\
         inst\tsetup\t33\treference/inst.vcf.gz\t\thttps://h/inst.vcf.gz\t\t\t\t\n",
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
fn duplicated_pins(s: &Sources, files: &[(&str, String)]) -> Vec<(String, String, String)> {
    let mut out = vec![];
    for (id, pin) in literal_pins(s) {
        for (name, text) in files {
            if text.contains(&pin) {
                out.push((id.clone(), pin.clone(), (*name).to_owned()));
            }
        }
    }
    out
}

fn pin_files() -> Vec<(&'static str, String)> {
    ["versions.env", "scripts/setup.sh", "scripts/lib/common.sh"]
        .into_iter()
        .map(|f| (f, read(f)))
        .collect()
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

/// `$x` written as `${x}`, and quotes dropped, so two spellings of a path compare.
fn norm(s: &str) -> String {
    let s = s.trim_matches(|c| c == '"' || c == '\'');
    let mut out = String::new();
    let b: Vec<char> = s.chars().collect();
    let mut i = 0;
    while i < b.len() {
        if b[i] == '$' && i + 1 < b.len() && (b[i + 1].is_ascii_alphabetic() || b[i + 1] == '_') {
            let mut j = i + 1;
            while j < b.len() && (b[j].is_ascii_alphanumeric() || b[j] == '_') {
                j += 1;
            }
            out.push_str(&format!("${{{}}}", b[i + 1..j].iter().collect::<String>()));
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

/// The words of a shell line, quotes kept inside a word.
fn words(l: &str) -> Vec<String> {
    let mut out = vec![];
    let mut cur = String::new();
    let (mut q, mut depth) = (None::<char>, 0i32);
    for c in l.chars() {
        match (q, c) {
            (None, '"' | '\'') => {
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
            (_, ')') if depth > 0 => {
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
    }
    if !cur.is_empty() {
        out.push(cur);
    }
    out
}

/// Every place in TEXT that deletes a raw download: in one function, an
/// `rm` (or a `find ... -delete`/`-exec rm`) whose argument is the
/// destination of a `fetch` there, a folder holding it, or a path in the
/// downloads folder. A fetch's own `.part`, which it resumes, is not a
/// download. `done_with` is the one place that may: it deletes a download
/// only when the user asked for it (`SETUP_KEEP_DOWNLOADS=false`). A line
/// that says `# not kept:` and why deletes no raw download (a staging copy,
/// an unpacked tree, a file no source of the table names).
fn deleted_downloads(name: &str, text: &str) -> Vec<String> {
    let mut out = vec![];
    for (func, lines) in blocks(text) {
        if func == "done_with" {
            continue;
        }
        let mut dests: Vec<String> = vec![];
        let mut dl_vars: Vec<String> = vec![];
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
            if let Some((lhs, rhs)) = l.trim().trim_start_matches("local ").split_once('=')
                && !lhs.contains(' ')
                && (rhs.contains("$(dl ")
                    || rhs.contains("source_of ")
                    || rhs.contains("OG_DOWNLOADS"))
            {
                dl_vars.push(format!("${{{lhs}}}"));
            }
        }
        for (n, l) in &lines {
            if l.trim_start().starts_with('#') || l.contains("# not kept:") {
                continue;
            }
            let w = words(l);
            let del = w.iter().position(|x| x == "rm").or_else(|| {
                (w.iter().any(|x| x == "find")
                    && (w.iter().any(|x| x == "-delete") || l.contains("-exec rm")))
                .then_some(0)
            });
            let Some(at) = del else { continue };
            for a in w[at + 1..].iter().filter(|a| !a.starts_with('-')) {
                let a = norm(a);
                if a.ends_with(".part") {
                    continue;
                }
                // A holds D: A is D, a file inside D, or a folder above D.
                let holds = |d: &str| {
                    a == d
                        || a.starts_with(&format!("{d}/"))
                        || d.strip_prefix(a.as_str())
                            .is_some_and(|r| r.starts_with('/') || r.starts_with("${"))
                };
                let hit = a.contains("OG_DOWNLOADS")
                    || a.contains("$(dl ")
                    || dl_vars
                        .iter()
                        .any(|v| holds(v) || a.starts_with(&format!("{v}/")))
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

#[test]
fn no_installer_deletes_a_raw_download() {
    let mut bad = deleted_downloads("scripts/setup.sh", &read("scripts/setup.sh"));
    bad.extend(deleted_downloads(
        "scripts/lib/common.sh",
        &read("scripts/lib/common.sh"),
    ));
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
}

/// Every network read inside a container in TEXT that has no `# net:`
/// comment in the three lines above it saying why it reads the network.
/// A read of a source's data has to say `# net: slice`: it reads only the
/// regions a user chose and keeps what it read, never a whole file unsaved.
fn unsaved_streams(name: &str, text: &str) -> Vec<String> {
    let lines: Vec<&str> = text.lines().collect();
    let mut out = vec![];
    for (i, l) in lines.iter().enumerate() {
        let w = words(l);
        if !(w.iter().any(|x| x == "run_in")
            && w.iter().any(|x| x == "--net" || x.contains("net[@]")))
        {
            continue;
        }
        let above = &lines[i.saturating_sub(3)..i];
        if !above.iter().any(|a| a.trim_start().starts_with("# net:")) {
            out.push(format!("{name}:{}: {}", i + 1, l.trim()));
        }
    }
    for (i, l) in lines.iter().enumerate() {
        let t = l.trim_start();
        if t.starts_with('#') {
            continue;
        }
        // curl or wget writing a body anywhere but through fetch.
        let w = words(l);
        let tool = w.iter().position(|x| x == "curl" || x == "wget");
        if let Some(p) = tool {
            let head_only = w[p + 1..]
                .iter()
                .any(|x| x.starts_with("-") && !x.starts_with("--") && x.contains('I'));
            let message = w[..p]
                .iter()
                .any(|x| x == "echo" || x.starts_with("\"") || x == "command");
            if !head_only && !message {
                out.push(format!("{name}:{}: {}", i + 1, l.trim()));
            }
        }
    }
    out
}

#[test]
fn no_source_is_streamed_without_being_saved() {
    let mut bad = unsaved_streams("scripts/setup.sh", &read("scripts/setup.sh"));
    bad.extend(unsaved_streams(
        "scripts/lib/common.sh",
        &read("scripts/lib/common.sh"),
    ));
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
    let root = checkout("a\tsetup\t1\tdownloads/a\t\thttps://h/a\t\t\t\t\n", "");
    let s = Sources::read(&root).unwrap();
    let g = ctx(&root, &[("GENOME_DIR", "/g".into())]);
    assert_eq!(
        g.path(s.get("a").unwrap(), None).unwrap(),
        PathBuf::from("/g/downloads/a")
    );
}
