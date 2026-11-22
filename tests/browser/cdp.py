"""The headless Chrome the Python live-site tests drive, moved here from
tests/test_owngenome_live.py when that file's checks became cargo tests
(rust/crates/og-live/tests/port_live_auth_ingest*.rs). New browser checks
are Rust (the og-browser crate, tests/browser/README.md); these helpers stay
only for the Python tests that still use them.

  find_chrome()          the Chrome or Chromium to run, or None
  Chrome(chrome, dir)    a headless Chrome over --remote-debugging-pipe
  asset_problems(texts)  outside URLs, inline scripts, on* attributes, eval
"""
import json
import os
import re
import select
import shutil
import subprocess
import sys
import time


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
