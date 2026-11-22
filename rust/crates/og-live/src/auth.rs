//! The live site's login: the password hash, the sessions, the limit on
//! failed logins, the client behind a trusted proxy, and HTTPS.
//!
//! The password is never stored. `og live serve hash-password` prints a
//! PBKDF2-SHA256 hash (`pbkdf2_sha256$<iterations>$<salt>$<hex>`, the format
//! of `bin/owngenome_live.py hash-password`, so a site keeps its login), and
//! the server reads only that hash, from `--password-file` or
//! `OWNGENOME_LIVE_PASSWORD_HASH`. A login sets a random cookie; the server
//! keeps only the SHA-256 of it, and the access log only the first 8
//! characters of that.

use crate::server::{COOKIE, Req, quote};
use crate::state::{lock, now};
use crate::static_files as sf;
use sha2::{Digest, Sha256};
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, Read, Write};
use std::net::{IpAddr, TcpStream};
use std::path::Path;
use std::process::{Command, Stdio};
use std::sync::{Arc, Mutex};
use std::time::Duration;

/// The rounds `hash-password` writes.
pub const ITERATIONS: u32 = 600_000;
/// The window failed logins are counted in.
pub const FAIL_WINDOW_SECONDS: f64 = 15.0 * 60.0;
/// How long a client may take to finish the TLS handshake.
pub const TLS_HANDSHAKE: Duration = Duration::from_secs(10);
/// The largest login form read.
const MAX_FORM: i64 = 4096;

// --- the password ---------------------------------------------------------------------------

/// The two SHA-256 states of HMAC with `key`: after the inner and after
/// the outer padded key.
fn hmac_keys(key: &[u8]) -> (Sha256, Sha256) {
    let mut k = [0u8; 64];
    if key.len() > 64 {
        k[..32].copy_from_slice(&Sha256::digest(key));
    } else {
        k[..key.len()].copy_from_slice(key);
    }
    let mut inner = Sha256::new();
    inner.update(k.map(|b| b ^ 0x36));
    let mut outer = Sha256::new();
    outer.update(k.map(|b| b ^ 0x5c));
    (inner, outer)
}

/// HMAC-SHA256 of `data` under the key of `keys`.
fn hmac(keys: &(Sha256, Sha256), data: &[u8]) -> [u8; 32] {
    let mut inner = keys.0.clone();
    inner.update(data);
    let mut outer = keys.1.clone();
    outer.update(inner.finalize());
    let mut out = [0u8; 32];
    out.copy_from_slice(&outer.finalize());
    out
}

/// PBKDF2-HMAC-SHA256 of `password` and `salt`, one 32-byte block
/// (`hashlib.pbkdf2_hmac("sha256", ...)`).
pub fn pbkdf2_sha256(password: &[u8], salt: &[u8], iterations: u32) -> [u8; 32] {
    let keys = hmac_keys(password);
    let mut first = salt.to_vec();
    first.extend_from_slice(&1u32.to_be_bytes());
    let mut u = hmac(&keys, &first);
    let mut out = u;
    for _ in 1..iterations {
        u = hmac(&keys, &u);
        for (o, x) in out.iter_mut().zip(u) {
            *o ^= x;
        }
    }
    out
}

fn hex(bytes: &[u8]) -> String {
    bytes.iter().map(|b| format!("{b:02x}")).collect()
}

/// `n` random bytes from the system.
pub fn random_bytes(n: usize) -> io::Result<Vec<u8>> {
    let mut b = vec![0u8; n];
    File::open("/dev/urandom")?.read_exact(&mut b)?;
    Ok(b)
}

/// URL-safe base64 without padding (`secrets.token_urlsafe`).
pub fn base64url(bytes: &[u8]) -> String {
    const ABC: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-_";
    let mut out = String::with_capacity(bytes.len().div_ceil(3) * 4);
    for chunk in bytes.chunks(3) {
        let n = chunk
            .iter()
            .enumerate()
            .fold(0u32, |n, (i, b)| n | (u32::from(*b) << (16 - 8 * i)));
        for i in 0..=chunk.len() {
            out.push(char::from(ABC[(n >> (18 - 6 * i)) as usize & 63]));
        }
    }
    out
}

/// A random token of 32 bytes, as 43 URL-safe characters.
pub fn new_token() -> io::Result<String> {
    random_bytes(32).map(|b| base64url(&b))
}

/// The hash line of `password` with a new random salt.
pub fn hash_password(password: &str, iterations: u32) -> io::Result<String> {
    let salt = hex(&random_bytes(16)?);
    let dk = pbkdf2_sha256(password.as_bytes(), salt.as_bytes(), iterations);
    Ok(format!("pbkdf2_sha256${iterations}${salt}${}", hex(&dk)))
}

/// (iterations, salt, hex digest) of a hash line, or why it is not one.
pub fn parse_hash(encoded: &str) -> Result<(u32, &str, &str), String> {
    let bad = || {
        "not a og live serve hash-password line (pbkdf2_sha256$<iterations>$<salt>$<hash>)"
            .to_owned()
    };
    let parts: Vec<&str> = encoded.trim().split('$').collect();
    let [kind, rounds, salt, digest] = parts[..] else {
        return Err(bad());
    };
    let rounds = (!rounds.is_empty() && rounds.bytes().all(|b| b.is_ascii_digit()))
        .then(|| rounds.parse::<u32>().ok())
        .flatten()
        .filter(|n| *n > 0);
    let digest_ok = digest.len() == 64
        && digest
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b));
    match rounds {
        Some(n) if kind == "pbkdf2_sha256" && !salt.is_empty() && digest_ok => {
            Ok((n, salt, digest))
        }
        _ => Err(bad()),
    }
}

/// Equal bytes, in time that does not depend on where they differ.
pub fn same_secret(a: &[u8], b: &[u8]) -> bool {
    a.len() == b.len() && a.iter().zip(b).fold(0u8, |acc, (x, y)| acc | (x ^ y)) == 0
}

/// Whether `password` is the one `encoded` hashes; false for a line that is
/// not a hash.
pub fn verify_password(password: &str, encoded: &str) -> bool {
    let Ok((rounds, salt, digest)) = parse_hash(encoded) else {
        return false;
    };
    let dk = pbkdf2_sha256(password.as_bytes(), salt.as_bytes(), rounds);
    same_secret(hex(&dk).as_bytes(), digest.as_bytes())
}

/// The hash from `--password-file` (else `env`, the value of
/// `OWNGENOME_LIVE_PASSWORD_HASH`): `None` for none set, an error for a file
/// that cannot be read or a value that is not a hash line.
pub fn read_password_hash(
    file: Option<&Path>,
    env: Option<&str>,
) -> Result<Option<String>, String> {
    let h = match file {
        Some(p) => std::fs::read_to_string(p)
            .map_err(|e| format!("{}: {e}", p.display()))?
            .trim()
            .to_owned(),
        None => env.unwrap_or("").trim().to_owned(),
    };
    if h.is_empty() {
        return Ok(None);
    }
    parse_hash(&h)?;
    Ok(Some(h))
}

/// Sets the terminal's echo on or off; nothing when stdin is no terminal.
fn echo(on: bool) {
    let _ = Command::new("stty")
        .arg(if on { "echo" } else { "-echo" })
        .stdin(Stdio::inherit())
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status();
}

/// One line of `input` without its line break.
fn read_line(input: &mut dyn BufRead) -> String {
    let mut line = String::new();
    let _ = input.read_line(&mut line);
    line.trim_end_matches(['\r', '\n']).to_owned()
}

/// `og live serve hash-password`: reads a password (twice, without echo,
/// from a terminal; else one line of `input`) and prints its hash.
pub fn hash_password_cli(
    args: &[String],
    input: &mut dyn BufRead,
    tty: bool,
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> u8 {
    if args.iter().any(|a| a == "-h" || a == "--help") {
        let _ = writeln!(
            out,
            "usage: og live serve hash-password\n\nPrint the hash of a password for --password-file."
        );
        return 0;
    }
    if let Some(a) = args.first() {
        let _ = writeln!(
            err,
            "og live serve hash-password: error: unrecognized arguments: {a}"
        );
        return 2;
    }
    let pw = if tty {
        let ask = |err: &mut dyn Write, input: &mut dyn BufRead, prompt: &str| {
            let _ = write!(err, "{prompt}");
            let _ = err.flush();
            echo(false);
            let line = read_line(input);
            echo(true);
            let _ = writeln!(err);
            line
        };
        let first = ask(err, input, "Password: ");
        if first != ask(err, input, "Again: ") {
            let _ = writeln!(err, "ERROR: the two passwords differ");
            return 2;
        }
        first
    } else {
        read_line(input)
    };
    if pw.chars().count() < 8 {
        let _ = writeln!(err, "ERROR: use a password of 8 characters or more");
        return 2;
    }
    match hash_password(&pw, ITERATIONS) {
        Ok(h) => {
            let _ = writeln!(out, "{h}");
            0
        }
        Err(e) => {
            let _ = writeln!(err, "ERROR: no random salt: {e}");
            1
        }
    }
}

// --- sessions and failed logins ---------------------------------------------------------------

/// The SHA-256 of a cookie, in hex: the only form of it the server keeps.
pub fn session_id(token: &str) -> String {
    hex(&Sha256::digest(token.as_bytes()))
}

#[derive(Clone, Copy)]
struct Session {
    created: f64,
    seen: f64,
}

/// The password and the logins of a running site.
pub struct Auth {
    /// The hash logins are checked against; `None`: no login.
    pub hash: Option<String>,
    /// A login lasts at most this long.
    pub session_hours: f64,
    /// A login unused this long ends.
    pub session_idle_minutes: f64,
    /// Failed logins per client in [`FAIL_WINDOW_SECONDS`].
    pub login_max_failures: usize,
    sessions: Mutex<HashMap<String, Session>>,
    failures: Mutex<HashMap<String, Vec<f64>>>,
}

impl Default for Auth {
    fn default() -> Auth {
        Auth::new(None)
    }
}

impl Auth {
    /// A site with password hash `hash` and the default limits.
    pub fn new(hash: Option<String>) -> Auth {
        Auth {
            hash,
            session_hours: 12.0,
            session_idle_minutes: 120.0,
            login_max_failures: 5,
            sessions: Mutex::new(HashMap::new()),
            failures: Mutex::new(HashMap::new()),
        }
    }

    fn expired(&self, s: Session, t: f64) -> bool {
        t - s.created > self.session_hours * 3600.0 || t - s.seen > self.session_idle_minutes * 60.0
    }

    /// A new session: its cookie.
    pub fn new_session(&self) -> io::Result<String> {
        let token = new_token()?;
        let t = now();
        lock(&self.sessions).insert(
            session_id(&token),
            Session {
                created: t,
                seen: t,
            },
        );
        Ok(token)
    }

    /// The session id of cookie `token` while it is live (this counts as a
    /// visit); an expired one is dropped.
    pub fn session(&self, token: Option<&str>) -> Option<String> {
        let sid = session_id(token.filter(|t| !t.is_empty())?);
        let t = now();
        let mut all = lock(&self.sessions);
        let s = all.get_mut(&sid)?;
        if self.expired(*s, t) {
            all.remove(&sid);
            return None;
        }
        s.seen = t;
        Some(sid)
    }

    /// Whether cookie `token`'s session is still live, without counting
    /// this as a visit: an open event stream must not keep an idle session
    /// alive.
    pub fn alive(&self, token: Option<&str>) -> bool {
        let Some(token) = token.filter(|t| !t.is_empty()) else {
            return false;
        };
        let t = now();
        lock(&self.sessions)
            .get(&session_id(token))
            .is_some_and(|s| !self.expired(*s, t))
    }

    /// Ends session `sid` (a logout).
    pub fn end_session(&self, sid: &str) {
        lock(&self.sessions).remove(sid);
    }

    /// False while `client` is blocked; else counts this attempt now, under
    /// the same lock, so parallel wrong guesses cannot all pass the check (a
    /// right password removes the count).
    pub fn reserve_attempt(&self, client: &str) -> bool {
        let t = now();
        let mut all = lock(&self.failures);
        let f = all.entry(client.to_owned()).or_default();
        f.retain(|x| t - x < FAIL_WINDOW_SECONDS);
        if f.len() >= self.login_max_failures {
            return false;
        }
        f.push(t);
        true
    }

    /// Forgets `client`'s failed logins (after a right password).
    pub fn clear_failures(&self, client: &str) {
        lock(&self.failures).remove(client);
    }
}

/// The session cookie line: `value` for `max_age` seconds, `Secure` when
/// the browser reaches the site over HTTPS.
pub fn cookie(value: &str, max_age: i64, secure: bool) -> String {
    let s = if secure { "; Secure" } else { "" };
    format!("{COOKIE}={value}; Path=/; HttpOnly; SameSite=Strict; Max-Age={max_age}{s}")
}

/// The value of the site's cookie in a `Cookie` header.
pub fn cookie_token(header: Option<&str>) -> Option<String> {
    header?.split(';').find_map(|part| {
        let (k, v) = part.trim().split_once('=').unwrap_or((part.trim(), ""));
        (k == COOKIE).then(|| v.to_owned())
    })
}

// --- the client behind a proxy --------------------------------------------------------------

/// An address block of `--trust-proxy`: an address, or `ADDR/PREFIX`
/// (host bits ignored, like `ipaddress.ip_network(strict=False)`).
#[derive(Debug, Clone, PartialEq)]
pub struct Net {
    addr: IpAddr,
    prefix: u8,
}

impl Net {
    pub fn parse(text: &str) -> Result<Net, String> {
        let bad = || format!("--trust-proxy {text}: not an address or ADDR/PREFIX");
        let (a, p) = text.trim().split_once('/').unwrap_or((text.trim(), ""));
        let addr: IpAddr = a
            .trim_matches(|c| c == '[' || c == ']')
            .parse()
            .map_err(|_| bad())?;
        let max = if addr.is_ipv4() { 32 } else { 128 };
        let prefix = if p.is_empty() {
            max
        } else {
            p.parse::<u8>().ok().filter(|n| *n <= max).ok_or_else(bad)?
        };
        Ok(Net { addr, prefix })
    }

    pub fn contains(&self, ip: IpAddr) -> bool {
        let bits = |ip: IpAddr| -> Option<(u128, u32)> {
            match ip {
                IpAddr::V4(v) => Some((u128::from(u32::from(v)), 32)),
                IpAddr::V6(v) => Some((u128::from(v), 128)),
            }
        };
        let (Some((a, w)), Some((b, w2))) = (bits(self.addr), bits(ip)) else {
            return false;
        };
        if w != w2 {
            return false;
        }
        let host = w - u32::from(self.prefix);
        host >= 128 || (a >> host) == (b >> host)
    }
}

/// Who a request comes from: the peer, or the last `X-Forwarded-For` entry
/// when the peer is a proxy of `--trust-proxy`.
pub fn client(peer: &str, forwarded_for: Option<&str>, trust: &[Net]) -> String {
    if let Ok(ip) = peer.parse::<IpAddr>()
        && trust.iter().any(|n| n.contains(ip))
        && let Some(last) = forwarded_for
            .unwrap_or("")
            .split(',')
            .map(str::trim)
            .rfind(|x| !x.is_empty())
    {
        return last.to_owned();
    }
    peer.to_owned()
}

// --- the routes ---------------------------------------------------------------------------

/// Reads the request's session: sets its id and the 8 characters the
/// access log keeps. Whether the request may see the runs: a live session,
/// or a site with no password.
pub fn begin(req: &mut Req<'_>) -> bool {
    let auth = Arc::clone(&req.live.cfg.auth);
    if auth.hash.is_none() {
        return true;
    }
    req.sid = auth.session(cookie_token(req.header("cookie")).as_deref());
    req.session8 = req.sid.as_deref().map_or("-".into(), |s| s[..8].to_owned());
    req.sid.is_some()
}

/// What a request without a session gets: a page (`/` or a path ending in
/// `/`) goes to the login form, anything else is 401.
pub fn refuse(req: &mut Req<'_>) -> io::Result<()> {
    let page = req.path == "/" || req.path.ends_with('/');
    req.log("no-session", "-", "-", if page { 303 } else { 401 });
    if page {
        let to = format!("/login?next={}", quote(&req.path));
        return req.redirect(303, &to);
    }
    req.deny(401, None)
}

/// Whether the browser reaches the site over HTTPS (TLS here, or at a
/// proxy that says so with `--cookie-secure`).
fn secure(req: &Req<'_>) -> bool {
    req.live.cfg.cookie_secure || req.live.cfg.tls.is_some()
}

fn login_page(
    req: &mut Req<'_>,
    code: u16,
    message: &str,
    extra: &[(&str, String)],
) -> io::Result<()> {
    let next = sf::query_value(&req.query, "next").unwrap_or_else(|| "/".into());
    match sf::login_page(&req.live.cfg.assets, message, &next) {
        Ok(page) => req.send(
            code,
            "text/html; charset=utf-8",
            page.as_bytes(),
            extra,
            crate::csp::CSP,
        ),
        Err(_) => req.deny(500, None),
    }
}

/// Where a login goes next: a path of this site, else `/`.
fn next_ok(next: &str) -> bool {
    next.starts_with('/')
        && !next.starts_with("//")
        && next[1..].bytes().all(|b| {
            b.is_ascii_alphanumeric() || matches!(b, b'.' | b'_' | b'~' | b'/' | b'%' | b'-')
        })
}

/// `/login`: the form; a POST checks the password.
pub fn login(req: &mut Req<'_>) -> io::Result<()> {
    if req.method != "POST" {
        return login_page(req, 200, "", &[]);
    }
    let auth = Arc::clone(&req.live.cfg.auth);
    let client = req.client.clone();
    let Some(hash) = auth.hash.clone() else {
        return req.redirect(303, "/");
    };
    if !auth.reserve_attempt(&client) {
        req.log("login-blocked", "-", "-", 429);
        return login_page(
            req,
            429,
            "Too many failed logins. Wait 15 minutes.",
            &[("Retry-After", "900".into())],
        );
    }
    let n = match req.header("content-length") {
        None => 0,
        Some(v) => v.trim().parse::<i64>().unwrap_or(-1),
    };
    if !(0..=MAX_FORM).contains(&n) {
        return req.deny(413, Some("too large"));
    }
    let body = req.read_body(n as usize)?;
    let form = String::from_utf8_lossy(&body).into_owned();
    let pw = sf::query_value(&form, "password").unwrap_or_default();
    let next = sf::query_value(&form, "next")
        .filter(|n| next_ok(n))
        .unwrap_or_else(|| "/".into());
    if !verify_password(&pw, &hash) {
        req.log("login-fail", "-", "-", 401);
        return login_page(req, 401, "Wrong password.", &[]);
    }
    auth.clear_failures(&client);
    let token = auth.new_session()?;
    req.session8 = session_id(&token)[..8].to_owned();
    req.log("login", "-", "-", 303);
    let max_age = (auth.session_hours * 3600.0) as i64;
    let set = cookie(&token, max_age, secure(req));
    req.send(
        303,
        "text/plain",
        b"",
        &[("Location", next), ("Set-Cookie", set)],
        crate::csp::CSP,
    )
}

/// `POST /logout`: ends the session and clears the cookie.
pub fn logout(req: &mut Req<'_>) -> io::Result<()> {
    if let Some(sid) = req.sid.clone() {
        req.live.cfg.auth.end_session(&sid);
    }
    req.log("logout", "-", "-", 303);
    let clear = cookie("", 0, secure(req));
    req.send(
        303,
        "text/plain",
        b"",
        &[("Location", "/login".into()), ("Set-Cookie", clear)],
        crate::csp::CSP,
    )
}

// --- HTTPS --------------------------------------------------------------------------------

/// The TLS settings of `--tls-cert` (a PEM chain) and `--tls-key`.
pub fn tls_config(cert: &Path, key: &Path) -> Result<Arc<rustls::ServerConfig>, String> {
    use rustls::pki_types::pem::PemObject;
    use rustls::pki_types::{CertificateDer, PrivateKeyDer};
    let certs = CertificateDer::pem_file_iter(cert)
        .and_then(Iterator::collect::<Result<Vec<_>, _>>)
        .map_err(|e| format!("--tls-cert {}: {e}", cert.display()))?;
    if certs.is_empty() {
        return Err(format!("--tls-cert {}: no certificate", cert.display()));
    }
    let key = PrivateKeyDer::from_pem_file(key)
        .map_err(|e| format!("--tls-key {}: {e}", key.display()))?;
    let provider = Arc::new(rustls::crypto::ring::default_provider());
    let cfg = rustls::ServerConfig::builder_with_provider(provider)
        .with_safe_default_protocol_versions()
        .map_err(|e| e.to_string())?
        .with_no_client_auth()
        .with_single_cert(certs, key)
        .map_err(|e| format!("--tls-cert and --tls-key: {e}"))?;
    Ok(Arc::new(cfg))
}

/// A connection after its TLS handshake, done in the connection's own
/// thread within [`TLS_HANDSHAKE`], so a client that connects and sends
/// nothing holds up no other; `None` when the handshake fails.
pub fn accept_tls(
    cfg: &Arc<rustls::ServerConfig>,
    mut tcp: TcpStream,
) -> Option<rustls::StreamOwned<rustls::ServerConnection, TcpStream>> {
    let mut conn = rustls::ServerConnection::new(Arc::clone(cfg)).ok()?;
    tcp.set_read_timeout(Some(TLS_HANDSHAKE)).ok()?;
    tcp.set_write_timeout(Some(TLS_HANDSHAKE)).ok()?;
    while conn.is_handshaking() {
        conn.complete_io(&mut tcp).ok()?;
    }
    Some(rustls::StreamOwned::new(conn, tcp))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pbkdf2_matches_the_rfc_7914_vector() {
        // RFC 7914 section 11: P="passwd", S="salt", c=1, dkLen=64 (the
        // first 32 bytes).
        assert_eq!(
            hex(&pbkdf2_sha256(b"passwd", b"salt", 1)),
            "55ac046e56e3089fec1691c22544b605f94185216dde0465e68b9d57c20dacbc"
        );
        // A key longer than a SHA-256 block is hashed first.
        let long = [b'k'; 100];
        assert_eq!(pbkdf2_sha256(&long, b"s", 2).len(), 32);
    }

    #[test]
    fn a_hash_verifies_its_password_and_no_other() {
        let h = hash_password("synthetic test password 7", 1000).unwrap();
        assert!(h.starts_with("pbkdf2_sha256$1000$"));
        assert!(verify_password("synthetic test password 7", &h));
        assert!(!verify_password("synthetic test password 7x", &h));
        assert!(!verify_password("x", "plain password"));
        for bad in [
            "",
            "pbkdf2_sha256$x$s$aa",
            "pbkdf2_sha256$0$s$",
            "md5$1$s$00",
        ] {
            assert!(parse_hash(bad).is_err(), "{bad}");
        }
    }

    #[test]
    fn tokens_are_43_url_safe_characters() {
        assert_eq!(base64url(b"\xfb\xff"), "-_8");
        assert_eq!(base64url(b"any carnal pleas"), "YW55IGNhcm5hbCBwbGVhcw");
        let t = new_token().unwrap();
        assert_eq!(t.len(), 43);
        assert!(
            t.bytes()
                .all(|b| b.is_ascii_alphanumeric() || b == b'-' || b == b'_')
        );
        assert_ne!(t, new_token().unwrap());
    }

    #[test]
    fn trusted_proxies_and_the_forwarded_client() {
        let trust = [
            Net::parse("127.0.0.1").unwrap(),
            Net::parse("203.0.113.3/24").unwrap(),
        ];
        assert_eq!(
            client("127.0.0.1", Some("192.0.2.1, 198.51.100.7"), &trust),
            "198.51.100.7"
        );
        assert_eq!(
            client("203.0.113.9", Some("192.0.2.1"), &trust),
            "192.0.2.1"
        );
        assert_eq!(
            client("192.0.2.5", Some("198.51.100.7"), &trust),
            "192.0.2.5"
        );
        assert_eq!(client("127.0.0.1", None, &trust), "127.0.0.1");
        assert!(Net::parse("::1").unwrap().contains("::1".parse().unwrap()));
        assert!(
            !Net::parse("::1")
                .unwrap()
                .contains("127.0.0.1".parse().unwrap())
        );
        assert!(
            Net::parse("0.0.0.0/0")
                .unwrap()
                .contains("192.0.2.1".parse().unwrap())
        );
        for bad in ["", "nope", "203.0.113.0/33", "203.0.113.0/x"] {
            assert!(Net::parse(bad).is_err(), "{bad}");
        }
    }

    #[test]
    fn the_cookie_is_read_from_its_header() {
        let h = "a=1; owngenome_live=tok; b=2";
        assert_eq!(cookie_token(Some(h)).as_deref(), Some("tok"));
        assert_eq!(cookie_token(Some("owngenome_livex=1")), None);
        assert_eq!(cookie_token(None), None);
        assert!(cookie("t", 5, true).ends_with("Max-Age=5; Secure"));
    }
}
