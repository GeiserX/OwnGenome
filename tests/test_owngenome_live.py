#!/usr/bin/env python3
"""bin/owngenome_live.py and assets/live/ on the synthetic DEMO-001 run (tests/fixtures/live).

Checks:
  1. hash-password: the hash verifies the password and no other; the server
     reads only a hash; a non-loopback --bind without one exits 2;
  2. login: no cookie gives a redirect to /login (pages) or 401 (data); a wrong
     password is 401 with a login-fail row; the right one sets an HttpOnly,
     SameSite=Strict cookie; an expired session is 401; after 5 failures the
     6th try is 429; the plain password is in no log, no stderr and no stored
     message; a negative or non-numeric Content-Length on /login is 413;
  3. Host: a foreign Host header is 403, an --allowed-host name is served;
  4. ingest: POST /ingest on the browser port is 404; a 2 MiB body is 413,
     a negative Content-Length 413, a body that is not JSON 400, a run named
     '../x' 400; metadata.parameters
     (an HPO list, a VCF path) and a task's script are not stored;
  5. an open event stream ends at logout;
  6. a plugin (tests/fixtures/live/plugins) adds a flag, a route, a run file
     and a state hook, all behind the session; the page gets one script tag
     per file in assets/live/plugins;
  7. in a headless Chrome (skipped without Chrome, except in CI): a page
     without a session goes to the login form, which leads back to the run
     page;
  8. with --tls-cert, a client that connects and sends nothing does not hold
     up another client's HTTPS request.
Checks 2, 4, 6 and 8 are also run against a planted mistake and must fail
there: a check that cannot fail is no check.

The rest of the site is checked against og live serve by cargo test -p og-live
(tests/port_live_server*.rs): the event stream's messages and Last-Event-ID,
state.json and the headlines, the policy on every answer, the run files, the
map, the re-check plugin and the step map in Chrome. The scan of assets/live/
is cargo test -p og-browser (tests/live_assets.rs), the replay mode
tests/live_replay.rs.

Run: python3 tests/test_owngenome_live.py
"""
import contextlib
import http.client
import io
import json
import os
import re
import select
import shutil
import socket
import subprocess
import sys
import tempfile
import threading
import time
import urllib.parse

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))
import owngenome_live  # noqa: E402
import run_state  # noqa: E402

FIX = os.path.join(REPO, "tests", "fixtures", "live")
ASSETS = os.path.join(REPO, "assets", "live")
RUN = "DEMO-001"
PASSWORD = "synthetic test password 7"
FAILS = []


def check(desc, ok, detail=""):
    print(("[PASS] " if ok else "[FAIL] ") + desc + ("" if ok else f" -- {str(detail)[:800]}"))
    if not ok:
        FAILS.append(desc)


def control(desc, still_passes):
    """A planted mistake: the check must fail on it."""
    check(f"control: {desc}", not still_passes, "the check passed on the planted mistake")


def request(port, method, path, body=None, headers=None):
    c = http.client.HTTPConnection("127.0.0.1", port, timeout=10)
    h = {"Host": f"127.0.0.1:{port}"}
    h.update(headers or {})
    c.request(method, path, body=body, headers=h)
    r = c.getresponse()
    data = r.read()
    hdr = {k.lower(): v for k, v in r.getheaders()}
    c.close()
    return r.status, hdr, data


def make_root(work):
    root = os.path.join(work, "genome")
    d = os.path.join(root, RUN)
    shutil.copytree(os.path.join(FIX, RUN), d)
    os.makedirs(os.path.join(d, "clinvar"))
    os.makedirs(os.path.join(d, "extra"))
    with open(os.path.join(d, "clinvar", "hits.tsv"), "w") as f:
        f.write("gene\tvariant\nGENEA\tchr1:100 A>G\n")
    with open(os.path.join(d, "clinvar", "unlinked.tsv"), "w") as f:
        f.write("not linked\n")
    with open(os.path.join(d, "extra", "note.txt"), "w") as f:
        f.write("plugin file\n")
    outside = os.path.join(work, "outside.txt")
    with open(outside, "w") as f:
        f.write("outside the root\n")
    os.symlink(outside, os.path.join(d, "clinvar", "escape.tsv"))
    # The report links CPSR's folder, which holds its HTML report and its table
    os.makedirs(os.path.join(d, "cpsr"))
    with open(os.path.join(d, "cpsr", f"{RUN}.cpsr.grch38.html"), "w") as f:
        f.write("<!doctype html><title>CPSR</title>\n")
    with open(os.path.join(d, "cpsr", f"{RUN}.cpsr.grch38.classification.tsv.gz"), "w") as f:
        f.write("not linked\n")
    os.makedirs(os.path.join(d, "twopages"))
    for n in ("a.html", "b.html"):
        with open(os.path.join(d, "twopages", n), "w") as f:
            f.write("<!doctype html>\n")
    with open(os.path.join(d, f"{RUN}_report.html"), "w") as f:
        f.write('<!doctype html><html><head><style>body{color:#eee}</style></head><body>'
                '<a href="clinvar/hits.tsv">hits</a> <a href="clinvar/escape.tsv">escape</a> '
                '<a href="cpsr/">the cpsr folder</a> <a href="twopages/">two pages</a> '
                '<a href="https://example.org/x">outside</a></body></html>\n')
    return root


def serve(argv, app=None):
    live, web, ingest = owngenome_live.start(argv, app)
    for s in (web, ingest):
        threading.Thread(target=s.serve_forever, daemon=True).start()
    return live, web, ingest


def stop(*servers):
    for live, web, ingest in servers:
        live.stopping = True
        web.shutdown()
        ingest.shutdown()
        web.server_close()
        ingest.server_close()


def login(port, password=PASSWORD):
    code, h, _ = request(port, "POST", "/login", urllib.parse.urlencode({"password": password, "next": "/"}),
                         {"Content-Type": "application/x-www-form-urlencoded"})
    m = re.match(r"owngenome_live=([^;]+)", h.get("set-cookie", ""))
    return code, h, (m.group(1) if m else None)


def post(port, run, ev):
    body = ev if isinstance(ev, (bytes, str)) else json.dumps(ev)
    return request(port, "POST", f"/ingest/{run}", body, {"Content-Type": "application/json"})


class Stream:
    """An open /events stream, read in a thread."""

    def __init__(self, port, cookie, last=None):
        self.closed = False
        self.items, self.sock = [], socket.create_connection(("127.0.0.1", port), timeout=30)
        hdr = f"GET /runs/{RUN}/events HTTP/1.1\r\nHost: 127.0.0.1:{port}\r\nCookie: owngenome_live={cookie}\r\n"
        if last is not None:
            hdr += f"Last-Event-ID: {last}\r\n"
        self.sock.sendall((hdr + "\r\n").encode())
        threading.Thread(target=self.read, daemon=True).start()

    def read(self):
        buf, cur = b"", {}
        try:
            while True:
                chunk = self.sock.recv(65536)
                if not chunk:
                    self.closed = True
                    return
                buf += chunk
                while b"\n" in buf:
                    line, buf = buf.split(b"\n", 1)
                    line = line.decode().rstrip("\r")
                    if line.startswith("id: "):
                        cur["id"] = int(line[4:])
                    elif line.startswith("data: "):
                        cur["data"] = json.loads(line[6:])
                    elif line == "" and "data" in cur:
                        self.items.append((cur["id"], cur["data"]))
                        cur = {}
        except OSError:
            return

    def wait(self, n, timeout):
        end = time.time() + timeout
        while time.time() < end and len(self.items) < n:
            time.sleep(0.02)
        return len(self.items) >= n

    def close(self):
        try:
            self.sock.close()
        except OSError:
            pass


def raw(port, head, body=b"", half_close=True, timeout=10):
    """Send a hand-made request and return the status line ('' when the server sent nothing)."""
    sk = socket.create_connection(("127.0.0.1", port), timeout=timeout)
    try:
        sk.sendall(head.encode() + body)
        if half_close:
            sk.shutdown(socket.SHUT_WR)  # a server that reads to the end returns instead of hanging
        return sk.recv(200).decode(errors="replace").split("\r\n", 1)[0]
    except OSError as e:
        return f"error {type(e).__name__}"
    finally:
        sk.close()


def tls_checks(root, work, pwfile):
    """11. One silent client must not stop the HTTPS listener. Returns (ok, detail) for a given
    wrap: the real server, and a planted one that shakes hands in accept()."""
    cert, key = os.path.join(work, "c.pem"), os.path.join(work, "k.pem")
    if not shutil.which("openssl"):
        ci = os.environ.get("CI") == "true"
        print(("[FAIL] " if ci else "[SKIP] ") + "no openssl to make a test certificate")
        if ci:
            FAILS.append("no openssl in CI")
        return
    r = subprocess.run(["openssl", "req", "-x509", "-newkey", "rsa:2048", "-nodes", "-keyout", key, "-out", cert,
                        "-days", "1", "-subj", "/CN=localhost"], capture_output=True, text=True, timeout=120)
    if r.returncode != 0:
        check("openssl makes a test certificate", False, r.stderr)
        return

    def one(planted):
        argv = ["--root", root, "--state-dir", os.path.join(work, "state-tls"), "--port", "0", "--ingest-port", "0",
                "--password-file", pwfile, "--tls-cert", cert, "--tls-key", key, "--quiet"]
        real_wrap = owngenome_live.ssl.SSLContext.wrap_socket
        if planted:
            def wrap(self, sock, **kw):
                kw.pop("do_handshake_on_connect", None)
                return real_wrap(self, sock, **kw)
            owngenome_live.ssl.SSLContext.wrap_socket = wrap
        try:
            with contextlib.redirect_stderr(io.StringIO()):
                srv = serve(argv)
        finally:
            owngenome_live.ssl.SSLContext.wrap_socket = real_wrap
        live, web, ingest = srv
        port = web.server_address[1]
        silent = socket.create_connection(("127.0.0.1", port), timeout=10)
        time.sleep(0.3)  # let accept() take the silent one first
        ctx = owngenome_live.ssl.create_default_context()
        ctx.check_hostname, ctx.verify_mode = False, owngenome_live.ssl.CERT_NONE
        try:
            c = http.client.HTTPSConnection("127.0.0.1", port, timeout=3, context=ctx)
            c.request("GET", "/login", headers={"Host": f"127.0.0.1:{port}"})
            code = c.getresponse().status
            c.close()
        except OSError as e:
            code = type(e).__name__
        silent.close()
        stop(srv)
        return code
    code = one(False)
    check("with --tls-cert, a silent client does not hold up another's HTTPS request", code == 200, code)
    control("the handshake done in accept()", one(True) == 200)


def asset_problems(texts):
    """(file, problem) for each outside URL, inline script or on* attribute in assets/live."""
    out = []
    for name, t in texts.items():
        for u in re.findall(r"(?:https?:)?//[A-Za-z0-9.-]+\.[a-z]{2,}[^\s'\"]*", t):
            if u != "http://www.w3.org/2000/svg":  # the SVG namespace name, never fetched
                out.append((name, u))
        if name.endswith(".html"):
            for tag in re.findall(r"<script\b[^>]*>", t, re.I):
                if "src=" not in tag:
                    out.append((name, tag))
            out += [(name, a) for a in re.findall(r"\son[a-z]+\s*=", t, re.I)]
        if name.endswith(".js") and re.search(r"\beval\s*\(|new\s+Function\s*\(", t):
            out.append((name, "eval"))
    return out


# --- the browser over the DevTools pipe ------------------------------------------------------

def find_chrome():
    for c in (os.environ.get("CHROME"), "google-chrome", "google-chrome-stable", "chromium", "chromium-browser",
              "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"):
        if c and (shutil.which(c) or os.path.isfile(c)):
            return shutil.which(c) or c
    return None


class Chrome:
    """Headless Chrome driven through --remote-debugging-pipe (fd 3 in, fd 4 out): no package needed."""

    def __init__(self, chrome, profile, tries=3):
        # A headless Chrome on a busy CI runner sometimes never answers its first
        # call (Target.createTarget, 30 s). Start a fresh one, up to three times.
        for attempt in range(1, tries + 1):
            try:
                return self._start(chrome, f"{profile}-{attempt}" if attempt > 1 else profile)
            except RuntimeError as e:
                self.p.kill()
                self.p.wait()
                os.close(self.to_w)
                os.close(self.from_r)
                if attempt == tries:
                    raise
                print(f"Chrome did not start ({e}); trying again ({attempt + 1} of {tries})", file=sys.stderr)

    def _start(self, chrome, profile):
        to_r, self.to_w = os.pipe()
        self.from_r, from_w = os.pipe()

        def fds():
            os.dup2(to_r, 3)
            os.dup2(from_w, 4)
        self.p = subprocess.Popen([chrome, "--headless=new", "--disable-gpu", "--no-sandbox", "--no-first-run",
                                   "--no-default-browser-check", f"--user-data-dir={profile}",
                                   "--remote-debugging-pipe", "about:blank"],
                                  stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, preexec_fn=fds, pass_fds=(3, 4))
        os.close(to_r)
        os.close(from_w)
        self.n, self.buf = 0, b""
        t = self.call("Target.createTarget", {"url": "about:blank"})["targetId"]
        self.session = self.call("Target.attachToTarget", {"targetId": t, "flatten": True})["sessionId"]

    def call(self, method, params=None, session=None, timeout=30):
        self.n += 1
        msg = {"id": self.n, "method": method, "params": params or {}}
        if session:
            msg["sessionId"] = session
        os.write(self.to_w, json.dumps(msg).encode() + b"\0")
        end = time.time() + timeout
        while time.time() < end:
            while b"\0" in self.buf:
                raw, self.buf = self.buf.split(b"\0", 1)
                m = json.loads(raw)
                if m.get("id") == self.n:
                    if "error" in m:
                        raise RuntimeError(f"{method}: {m['error']}")
                    return m.get("result", {})
            r, _, _ = select.select([self.from_r], [], [], 0.5)
            if r:
                chunk = os.read(self.from_r, 1 << 16)
                if not chunk:
                    raise RuntimeError("Chrome closed the pipe")
                self.buf += chunk
        raise RuntimeError(f"{method}: no answer in {timeout} s")

    def eval(self, expr):
        r = self.call("Runtime.evaluate", {"expression": expr, "returnByValue": True}, self.session)
        return r.get("result", {}).get("value")

    def go(self, url):
        self.call("Page.navigate", {"url": url}, self.session)

    def wait(self, expr, timeout=20):
        end = time.time() + timeout
        while time.time() < end:
            try:
                v = self.eval(expr)
            except RuntimeError:
                v = None
            if v:
                return v
            time.sleep(0.2)
        return None

    def close(self):
        try:
            self.call("Browser.close", timeout=5)
        except (RuntimeError, OSError):
            pass
        try:
            self.p.wait(timeout=10)
        except subprocess.TimeoutExpired:
            self.p.kill()


def browser_checks(port, work):
    chrome = find_chrome()
    if not chrome:
        ci = os.environ.get("CI") == "true"
        print(("[FAIL] " if ci else "[SKIP] ") + "no Chrome or Chromium to load the page in")
        if ci:
            FAILS.append("no Chrome in CI")
        return
    b = Chrome(chrome, os.path.join(work, "chrome"))
    try:
        base = f"http://127.0.0.1:{port}"
        b.go(f"{base}/runs/{RUN}/")
        check("browser: a page without a session goes to the login form",
              b.wait("location.pathname === '/login' && !!document.getElementById('password')"),
              b.eval("location.href"))
        b.eval(f"document.getElementById('password').value = {json.dumps(PASSWORD)}; "
               "document.querySelector('form').submit(); true")
        ok = b.wait(f"location.pathname === '/runs/{RUN}/' && document.documentElement.dataset.pgpStatus === 'ready'")
        check("browser: the login form leads back to the run page, which loads", ok,
              b.eval("location.href + ' ' + document.documentElement.outerHTML.slice(0, 600)"))
    finally:
        b.close()


INGEST_PORT = None


def main():
    global INGEST_PORT
    work = tempfile.mkdtemp(prefix="owngenome-live-")
    servers = []
    err = io.StringIO()
    try:
        root = make_root(work)
        state_dir = os.path.join(work, "state")
        pwfile = os.path.join(work, "password")

        # 1. the password
        r = subprocess.run([sys.executable, os.path.join(REPO, "bin", "owngenome_live.py"), "hash-password"],
                           input=PASSWORD + "\n", capture_output=True, text=True, timeout=60)
        h = r.stdout.strip()
        check("hash-password prints a pbkdf2_sha256 hash of 600000 rounds, not the password",
              r.returncode == 0 and h.startswith("pbkdf2_sha256$600000$") and PASSWORD not in r.stdout + r.stderr,
              (r.returncode, r.stderr))
        check("the hash verifies the password and no other",
              owngenome_live.verify_password(PASSWORD, h) and not owngenome_live.verify_password(PASSWORD + "x", h))
        with open(pwfile, "w") as f:
            f.write(h + "\n")
        env = {k: v for k, v in os.environ.items() if k != "OWNGENOME_LIVE_PASSWORD_HASH"}
        r = subprocess.run([sys.executable, os.path.join(REPO, "bin", "owngenome_live.py"), "--root", root, "--state-dir",
                            state_dir, "--bind", "0.0.0.0", "--port", "0", "--ingest-port", "0"],
                           capture_output=True, text=True, timeout=60, env=env)
        check("a non-loopback --bind without a password hash exits 2", r.returncode == 2 and "password" in r.stderr,
              (r.returncode, r.stderr))

        app = owngenome_live.App()
        argv = ["--root", root, "--state-dir", state_dir, "--port", "0", "--ingest-port", "0", "--password-file",
                pwfile, "--allowed-host", "live.example.org", "--plugin-dir", os.path.join(FIX, "plugins"),
                "--fixture-word", "planted-flag", "--quiet"]
        with contextlib.redirect_stderr(err):
            srv = serve(argv, app)
        servers.append(srv)
        live, web, ingest = srv
        port, INGEST_PORT = web.server_address[1], ingest.server_address[1]
        for line in open(os.path.join(FIX, "events.jsonl")):
            code, _, _ = post(INGEST_PORT, RUN, line.strip())
            if code != 200:
                check("the fixture's messages are taken", False, code)

        # 2. login
        code, hd, _ = request(port, "GET", f"/runs/{RUN}/")
        code2, _, _ = request(port, "GET", f"/runs/{RUN}/state.json")
        check("no cookie: a page redirects to /login, data is 401",
              code == 303 and hd.get("location", "").startswith("/login") and code2 == 401, (code, hd, code2))
        with contextlib.redirect_stderr(err):
            code, _, body = login(port, "wrong password")
        log = open(os.path.join(state_dir, "access.log")).read()
        check("a wrong password is 401 and writes a login-fail row", code == 401 and "\tlogin-fail\t" in log, (code, log))
        with contextlib.redirect_stderr(err):
            code, hd, cookie = login(port)
        sc = hd.get("set-cookie", "")
        check("the right password sets an HttpOnly, SameSite=Strict session cookie",
              code == 303 and cookie and "HttpOnly" in sc and "SameSite=Strict" in sc, (code, sc))
        auth = {"Cookie": f"owngenome_live={cookie}"}
        code, _, _ = request(port, "GET", f"/runs/{RUN}/state.json", headers=auth)
        check("a valid session reads state.json", code == 200, code)
        code, _, _ = request(port, "GET", "/", headers={"Cookie": "owngenome_live=made-up"})
        check("a made-up cookie is not a session", code == 303, code)
        _, _, old = login(port)
        sid = __import__("hashlib").sha256(old.encode()).hexdigest()
        live.sessions[sid]["created"] -= live.cfg.session_hours * 3600 + 1
        code, _, _ = request(port, "GET", f"/runs/{RUN}/state.json", headers={"Cookie": f"owngenome_live={old}"})
        check("an expired session is 401", code == 401, code)
        live.failures.clear()
        codes = [login(port, f"wrong {i}")[0] for i in range(6)]
        check("five wrong passwords are 401, the sixth try is 429", codes == [401] * 5 + [429], codes)
        code, _, _ = login(port)
        check("the right password is refused too while blocked", code == 429, code)

        def burst(n=12):
            live.failures.clear()
            gate = threading.Barrier(n)
            got = []

            def one(i):
                gate.wait()
                try:
                    got.append(login(port, f"parallel wrong {i}")[0])
                except OSError:  # a 429 sent before the body is read can reset the connection
                    got.append("reset")
            ts = [threading.Thread(target=one, args=(i,)) for i in range(n)]
            for x in ts:
                x.start()
            for x in ts:
                x.join(60)
            return got
        got = burst()
        check("parallel wrong passwords get at most five checks; the rest are refused",
              got.count(401) <= live.cfg.login_max_failures and len(got) == 12
              and all(x in (401, 429, "reset") for x in got), sorted(map(str, got)))
        reserve = live.reserve_attempt
        live.reserve_attempt = lambda client: True  # planted: the attempt is never counted
        try:
            planted = burst()
        finally:
            live.reserve_attempt = reserve
        control("an attempt that is not counted before the password check",
                planted.count(401) <= live.cfg.login_max_failures)
        live.failures.clear()
        bad_len = [raw(port, f"POST /login HTTP/1.1\r\nHost: 127.0.0.1:{port}\r\nContent-Type: "
                             f"application/x-www-form-urlencoded\r\nContent-Length: {n}\r\n\r\n", b"password=x")
                   for n in ("-1", "abc")]
        check("a negative or non-numeric Content-Length on /login is 413, before the body is read",
              all(" 413 " in x for x in bad_len), bad_len)
        live.failures.clear()

        def leaks():
            files = [os.path.join(state_dir, "access.log"), live.store.path(RUN)]
            return any(PASSWORD in open(p).read() for p in files if os.path.exists(p)) or PASSWORD in err.getvalue()
        check("the plain password is in no log, no stderr and no stored message", not leaks())
        app.access.write("127.0.0.1", "-", "planted", RUN, PASSWORD, 200)
        control("the password written into the access log", not leaks())
        with open(os.path.join(state_dir, "access.log"), "w") as f:
            f.write(owngenome_live.AccessLog.HEADER)

        # 3. Host
        code, _, _ = request(port, "GET", "/", headers=dict(auth, Host="evil.example.net"))
        code2, _, _ = request(port, "GET", "/", headers=dict(auth, Host="live.example.org:8443"))
        check("a foreign Host is 403, an --allowed-host name is served", code == 403 and code2 == 200, (code, code2))

        # 4. ingest
        code, _, _ = request(port, "POST", f"/ingest/{RUN}", "{}", auth)
        check("POST /ingest on the browser port is 404", code == 404, code)
        # The server answers 413 from the Content-Length, before the body: send the head and a little.
        sk = socket.create_connection(("127.0.0.1", INGEST_PORT), timeout=10)
        sk.sendall(f"POST /ingest/{RUN} HTTP/1.1\r\nHost: 127.0.0.1\r\nContent-Type: application/json\r\n"
                   f"Content-Length: {2 * 1024 * 1024}\r\n\r\n".encode() + b"x" * 1024)
        head = sk.recv(200).decode(errors="replace")
        sk.close()
        check("a 2 MiB message is 413", head.startswith("HTTP/1.0 413") or head.startswith("HTTP/1.1 413"), head)
        line = raw(INGEST_PORT, f"POST /ingest/{RUN} HTTP/1.1\r\nHost: 127.0.0.1\r\nContent-Type: application/json"
                                f"\r\nContent-Length: -1\r\n\r\n", b'{"event":"started"}')
        check("a negative Content-Length on /ingest is 413", " 413 " in line, line)
        code, _, _ = post(INGEST_PORT, RUN, '{"event": "process_started", "trace": {"attempt": NaN, "exit": Infinity}}')
        check("a message with NaN or Infinity is 400, so events.jsonl stays plain JSON", code == 400, code)
        cfgf = os.path.join(work, "live.conf")
        with open(cfgf, "w") as f:
            f.write(f"root = {root}\n")
        got = owngenome_live.configure(["--config", cfgf], owngenome_live.App())
        check("--root may come from --config", got.root == root, got.root)
        try:
            owngenome_live.configure([], owngenome_live.App())
            no_root = "accepted"
        except SystemExit as e:
            no_root = e.code
        check("no --root anywhere stops with exit 2", no_root == 2, no_root)
        numeric = owngenome_live.cut({"event": "process_started", "trace": {"name": 123, "process": 7, "tag": 9, "status": 1}})
        try:
            owngenome_live.messages([numeric], run_state.read_steps())
            run_state.event_tasks([numeric], RUN)
            survived = True
        except TypeError as e:
            survived = e
        check("a trace with numeric name, process, tag or status is stored as text and read without error",
              survived is True and numeric["trace"]["name"] == "123", (survived, numeric))

        def dropped(p, handler, t):
            """Does a client that promises a body and sends none lose its connection within 4 s?"""
            saved = handler.timeout
            handler.timeout = t
            sk = socket.create_connection(("127.0.0.1", p), timeout=4)
            try:
                sk.sendall(f"POST /login HTTP/1.1\r\nHost: 127.0.0.1:{p}\r\nContent-Length: 4096\r\n\r\n".encode())
                while sk.recv(4096):
                    pass
                return True
            except OSError:
                return False
            finally:
                sk.close()
                handler.timeout = saved
        check("both listeners bound a request's reads (BrowserHandler and IngestHandler timeout set)",
              owngenome_live.BrowserHandler.timeout and owngenome_live.IngestHandler.timeout,
              (owngenome_live.BrowserHandler.timeout, owngenome_live.IngestHandler.timeout))
        check("a client that stops sending mid-request is dropped (browser and ingest, timeout 1 s)",
              dropped(port, owngenome_live.BrowserHandler, 1) and dropped(INGEST_PORT, owngenome_live.IngestHandler, 1))
        control("no timeout on the browser listener", dropped(port, owngenome_live.BrowserHandler, None))
        code, _, _ = post(INGEST_PORT, RUN, "not json")
        check("a message that is not JSON is 400", code == 400, code)
        code, _, _ = post(INGEST_PORT, "../x", {"event": "started"})
        code2, _, _ = post(INGEST_PORT, "-x", {"event": "started"})
        check("a run named '../x' or '-x' is 400", code == 400 and code2 == 400, (code, code2))
        code, _, _ = request(INGEST_PORT, "POST", f"/ingest/{RUN}", '{"event":"started"}',
                             {"Origin": "https://evil.example.net"})
        check("a browser's post (an Origin header) is 403", code == 403, code)
        stored = open(live.store.path(RUN)).read()

        def stored_leaks(text):
            return [x for x in ("HP:0001250", "DEMO-001.vcf.gz", "bcftools roh", "commandLine") if x in text]
        check("metadata.parameters, the command line and the task script are not stored",
              not stored_leaks(stored) and '"outdir":"/data/genome"' in stored, stored_leaks(stored))
        control("the raw message stored", not stored_leaks(open(os.path.join(FIX, "events.jsonl")).read()))

        # 5. the event stream
        n0 = len(live.store.events(RUN))
        _, _, gone = login(port)
        st3 = Stream(port, gone)
        st3.wait(n0, 5)
        k = len(st3.items)
        live.end_session(__import__("hashlib").sha256(gone.encode()).hexdigest())
        post(INGEST_PORT, RUN, {"event": "pgp_step", "step": "06", "state": "ok", "utcTime": "2026-01-02T10:06:40Z"})
        end = time.time() + 3
        while time.time() < end and not st3.closed:
            time.sleep(0.02)
        check("an open stream ends when its session does, and sends nothing after",
              k == n0 and st3.closed and len(st3.items) == k, (st3.closed, len(st3.items), k, n0))
        st3.close()

        # 6. plugins
        _, _, body = request(port, "GET", f"/runs/{RUN}/state.json", headers=auth)
        s = json.loads(body)
        check("a state hook's answer is under plugins", s["plugins"].get("fixture", {}).get("steps") == len(s["steps"]),
              s.get("plugins"))
        code, _, body = request(port, "GET", f"/fixture/{RUN}", headers=auth)
        code2, _, _ = request(port, "GET", f"/fixture/{RUN}")
        check("a plugin route answers with its own flag behind the session, 401 without",
              code == 200 and json.loads(body) == {"word": "planted-flag", "run": RUN, "dir": True} and code2 == 401,
              (code, body, code2))
        code, _, body = request(port, "GET", f"/runs/{RUN}/extra/note.txt", headers=auth)
        check("a plugin's run file is served", code == 200 and body == b"plugin file\n", (code, body))
        app.routes = [(m, rx, fn, False) for m, rx, fn, _ in app.routes]
        code2, _, _ = request(port, "GET", f"/fixture/{RUN}")
        control("the plugin route registered with auth skipped", code2 == 401)
        app.routes = [(m, rx, fn, True) for m, rx, fn, _ in app.routes]
        adir = os.path.join(work, "assets")
        shutil.copytree(ASSETS, adir)
        # start from no plugins, so the plugins this repo ships do not change the count
        for f in os.listdir(os.path.join(adir, "plugins")):
            os.remove(os.path.join(adir, "plugins", f))
        for f in ("a.js", "b.js", "notes.txt"):
            open(os.path.join(adir, "plugins", f), "w").write("// plugin\n")
        live.cfg.assets = adir
        _, _, body = request(port, "GET", f"/runs/{RUN}/", headers=auth)
        tags = re.findall(r'<script src="/assets/plugins/([^"]+)"', body.decode())
        check("the run page gets one script tag per plugin file", tags == ["a.js", "b.js"], tags)
        code, _, _ = request(port, "GET", "/assets/plugins/a.js")
        check("a plugin script is served", code == 200, code)
        live.cfg.assets = ASSETS

        # 7. the browser
        browser_checks(port, work)

        # 8. TLS
        tls_checks(root, work, pwfile)
    finally:
        stop(*servers)
        shutil.rmtree(work, ignore_errors=True)
    print()
    print("RESULT: ALL PASS" if not FAILS else f"RESULT: {len(FAILS)} FAILED")
    sys.exit(1 if FAILS else 0)


if __name__ == "__main__":
    main()
