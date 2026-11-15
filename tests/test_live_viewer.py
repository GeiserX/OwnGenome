#!/usr/bin/env python3
"""The genome viewer on the live site: bin/owngenome_live_viewer.py, started by
scripts/owngenome-live.sh, on the synthetic sample of tests/test_viewer.py.

Checks:
  1. owngenome-live.sh downloads the igv.js pinned in versions.env (sha256
     checked) and starts the site with the viewer: /runs/<run>/viewer/ without
     a login sends the browser to /login, and the config and the reads answer
     401;
  2. after the login the page carries the viewer's Content-Security-Policy
     and says it is served by the site; the config names the run's BAM, VCF
     and SV VCF; a byte range of the BAM is 206 with exactly those bytes, past
     the end 416; a file the config does not name (a wrong name, a path out
     of data/) is 404; a run that does not exist is 404;
  3. the report the site serves links each position to /runs/<run>/viewer/
     with the site's note, not to the local viewer on 127.0.0.1:8765;
  4. opening the page writes one access-log row, and igv.js's byte ranges
     write none;
  5. a run whose BAM is a symlink out of GENOME_DIR gets no reads track, a
     run whose VCF index is one gets no VCF track and neither file is
     served, and a reference named with a space, '#' and '?' is served from
     the encoded URL its config gives;
  6. a site started without the viewer flags answers the viewer with 404 and
     the reason, and leaves the report's links as they are;
  7. in a headless Chrome (skipped without Chrome, except in CI) the page,
     loaded from the site, starts igv.js (status ready, 3 tracks) and the
     policy blocks nothing.
  8. the map opens a same-origin reads iframe, follows the camera, hides
     closed findings and returns focus on Escape; a planted page policy
     without frame-src blocks the iframe. Other routes keep DENY.
Checks 2, 3 and 5 are also run against planted mistakes (a viewer that serves
any file of the run folder, a report filter that does nothing, a viewer that
compares the still-quoted name, one that serves an index out of GENOME_DIR)
and must fail.

Run: python3 tests/test_live_viewer.py
"""
import http.client
import base64
import gzip
import json
import importlib.util
import os
import re
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
sys.path.insert(0, os.path.join(REPO, "tests"))
import owngenome_live  # noqa: E402
import render_report  # noqa: E402
import test_viewer as tv  # noqa: E402
import test_owngenome_live as tl  # noqa: E402

PASSWORD = "viewer test password"
FAILED = []


def check(desc, ok, detail=""):
    print(("[PASS] " if ok else "[FAIL] ") + desc + ("" if ok else f" -- {str(detail)[:600]}"))
    if not ok:
        FAILED.append(desc)


def request(port, path, headers=None, method="GET", body=None):
    c = http.client.HTTPConnection("127.0.0.1", port, timeout=20)
    try:
        c.request(method, path, body=body, headers=dict({"Host": f"127.0.0.1:{port}"}, **(headers or {})))
        r = c.getresponse()
        return r.status, {k.lower(): v for k, v in r.getheaders()}, r.read()
    finally:
        c.close()


def login(port):
    code, h, _ = request(port, "/login", {"Content-Type": "application/x-www-form-urlencoded"}, "POST",
                         urllib.parse.urlencode({"password": PASSWORD, "next": "/"}))
    m = re.match(r"owngenome_live=([^;]+)", h.get("set-cookie", ""))
    return {"Cookie": f"owngenome_live={m.group(1)}"} if code == 303 and m else None


def wait_port(port, proc, seconds=120):
    for _ in range(seconds * 2):
        if proc.poll() is not None:
            return False
        try:
            with socket.create_connection(("127.0.0.1", port), timeout=1):
                return True
        except OSError:
            time.sleep(0.5)
    return False


def in_process(argv, plugin_dir=None):
    """The site in this process (plugins from plugin_dir, default bin/), on free ports."""
    argv = argv + ["--port", "0", "--ingest-port", "0", "--quiet"] + (["--plugin-dir", plugin_dir] if plugin_dir else [])
    live, web, ingest = owngenome_live.start(argv)
    for s in (web, ingest):
        threading.Thread(target=s.serve_forever, daemon=True).start()
    return live, web, ingest


def stop(live, web, ingest):
    live.stopping = True
    for s in (web, ingest):
        s.shutdown()
        s.server_close()


def viewer_checks(port, auth, s, label=""):
    """Checks 2: what the viewer serves to a signed-in browser. Returns the config."""
    import json
    code, h, body = request(port, "/runs/S/viewer/", auth)
    csp = h.get("content-security-policy", "")
    check(f"{label}page: 200 with the viewer's policy, served by the site",
          code == 200 and csp == tv.serve.CSP.replace("frame-ancestors 'none'", "frame-ancestors 'self'") and b"served by the run site" in body
          and b"served from this computer only" not in body, (code, csp, body[:300]))
    code, _, body = request(port, "/runs/S/viewer/igv.min.js", auth)
    c2, _, js = request(port, "/runs/S/viewer/viewer.js", auth)
    check(f"{label}page: igv.min.js and viewer.js are served beside it",
          code == 200 and len(body) > 100000 and c2 == 200 and b"config.json" in js, (code, c2))
    code, _, body = request(port, "/runs/S/viewer/config.json", auth)
    cfg = json.loads(body) if code == 200 else {}
    names = [t["url"].rsplit("/", 1)[1] for t in cfg.get("tracks", [])]
    check(f"{label}config: the run's BAM, VCF and SV VCF, by relative URL",
          names == ["S_sorted.bam", "S.vcf.gz", "S_sv_consensus.vcf.gz"]
          and all(t["url"].startswith("data/") for t in cfg.get("tracks", [])), (code, cfg))
    raw = open(f"{s}/aligned/S_sorted.bam", "rb").read()
    bam = f"/runs/S/viewer/{cfg['tracks'][0]['url']}" if cfg.get("tracks") else "/runs/S/viewer/data/x"
    code, h, body = request(port, bam, dict(auth, Range="bytes=10-109"))
    check(f"{label}reads: a BAM range is 206 with exactly those bytes",
          code == 206 and body == raw[10:110] and h.get("content-range") == f"bytes 10-109/{len(raw)}", (code, h))
    code, h, _ = request(port, bam, dict(auth, Range=f"bytes={len(raw) + 5}-"))
    check(f"{label}reads: a range past the end is 416", code == 416, code)
    bad = [request(port, p, auth)[0] for p in (
        "/runs/S/viewer/data/2/S_report.html", "/runs/S/viewer/data/99/S_sorted.bam",
        "/runs/S/viewer/S_report.html", "/runs/S/viewer/data/2/..%2F..%2FS_report.html",
        "/runs/S/viewer/serve.py", "/runs/NOPE/viewer/config.json")]
    check(f"{label}a file the config does not name, and a run that does not exist, are 404",
          bad == [404] * len(bad), bad)
    return cfg


def report_checks(port, auth, label=""):
    """Checks 3: the served report's links point at the site's viewer."""
    code, _, body = request(port, "/runs/S/report.html", auth)
    page = body.decode(errors="replace")
    check(f"{label}report: positions link to /runs/S/viewer/ with the site's note, none to 127.0.0.1:8765",
          code == 200 and 'href="/runs/S/viewer/#chr1:451"' in page and render_report.VIEWER_URL not in page
          and "this site's genome viewer" in page
          and render_report.VIEWER_NOTE not in page, (code, [x for x in page.splitlines() if "viewer" in x][:4]))


def outside_checks(port, auth, label=""):
    """Checks 5: files out of GENOME_DIR are left out; an odd reference name is served."""
    import json
    code, _, body = request(port, "/runs/T/viewer/config.json", auth)
    tracks = json.loads(body).get("tracks") if code == 200 else None
    check(f"{label}a BAM symlinked out of GENOME_DIR is left out of the config", code == 200 and tracks == [],
          (code, body[:300]))
    code, _, body = request(port, "/runs/U/viewer/config.json", auth)
    cfg = json.loads(body) if code == 200 else {}
    c2, _, _ = request(port, "/runs/U/viewer/data/2/U.vcf.gz", dict(auth, Range="bytes=0-9"))
    c3, _, _ = request(port, "/runs/U/viewer/data/3/U.vcf.gz.tbi", dict(auth, Range="bytes=0-9"))
    check(f"{label}a VCF whose index is symlinked out of GENOME_DIR: no track, neither file served",
          code == 200 and cfg.get("tracks") == [] and c2 == 404 and c3 == 404, (code, cfg, c2, c3))
    url = cfg.get("reference", {}).get("fastaURL", "")
    code, _, body = request(port, f"/runs/U/viewer/{url}", dict(auth, Range="bytes=0-4"))
    check(f"{label}a reference named with a space, '#' or '?' is served from its encoded URL",
          url == "data/0/my%20ref%231%3F.fasta" and code == 206 and body == b">chr1", (url, code, body))


def framing_checks(port):
    for path in ("/runs/S/viewer/", "/runs/S/viewer/viewer.js", "/runs/S/viewer/config.json",
                 "/runs/S/viewer/data/2/S_sorted.bam"):
        code, h, _ = request(port, path)
        csp = h.get("content-security-policy", "")
        check("same-origin viewer framing: " + path,
              code == 200 and "frame-ancestors 'self'" in csp
              and h.get("x-frame-options") == "SAMEORIGIN", (code, h))
    for path in ("/", "/login", "/assets/map.js", "/runs/S/", "/runs/S/state.json",
                 "/runs/S/report.html", "/not-found"):
        _, h, _ = request(port, path)
        check("non-viewer framing stays closed: " + path,
              "frame-ancestors 'none'" in h.get("content-security-policy", "")
              and h.get("x-frame-options") == "DENY", h)


def map_reads_checks(srv, work, chrome):
    """The real map and viewer, on invented map rows and test_viewer's synthetic reads.

    No product data loader is bypassed: the fixture endpoint supplies the same
    gzip document as map.json. The page supplies its existing veil callbacks.
    """
    app, port = srv[0].app, srv[1].server_address[1]
    policy = {"block": False}
    page = ('<!doctype html><html data-theme="dark"><head><meta name="viewport" content="width=device-width, initial-scale=1">'
            '<link rel="stylesheet" href="/assets/live.css">'
            '<link rel="stylesheet" href="/assets/map.css">'
            '<script src="/assets/core.js" defer></script>'
            '<script src="/assets/map.js" defer></script></head><body></body></html>')
    doc = {"v": 1, "people": [{"id": "S"}], "dict": {"veil": ["closed"]},
           "ref": {"chroms": [{"n": "chr1", "len": 60000, "masked": []}], "bands": {}},
           "per": {"S": {"findings": {
               "id": ["open", "closed"], "level": ["review", "act"], "gate": ["", ""],
               "chrom": ["chr1", "chr1"], "pos": [451, 551], "s": [450, 550], "e": [451, 551],
               "gene": ["", ""], "title": ["Synthetic finding", "Closed finding"],
               "href": ["#open", "#closed"], "phrase": ["Synthetic row", "Closed row"], "veil": [-1, 0]}}},
           "layers": {"findings": {"state": "landed"}, "roh": {"state": "not_run"}}, "text": {}}
    data = gzip.compress(json.dumps(doc).encode())

    def fixture(req, _):
        csp = owngenome_live.CSP
        if policy["block"]:
            csp = csp.replace("frame-src 'self'; ", "")
        return req.send_body(200, "text/html", page, csp=csp)

    app.route("GET", r"/map-1g/", fixture)
    app.route("GET", r"/map-1g.json", lambda req, _: req.send_body(200, "application/gzip", data))
    b = tl.Chrome(chrome, f"{work}/chrome-map-1g")
    try:
        visit = 0

        def js(expr):
            reply = b.call("Runtime.evaluate", {"expression": expr, "returnByValue": True}, b.session)
            if reply.get("exceptionDetails"):
                raise RuntimeError(str(reply["exceptionDetails"]))
            return reply.get("result", {}).get("value")

        b.call("Emulation.setEmulatedMedia", {"features": [{"name": "prefers-color-scheme", "value": "dark"},
                                                          {"name": "prefers-reduced-motion", "value": "reduce"}]}, b.session)
        b.call("Emulation.setDeviceMetricsOverride", {"width": 1440, "height": 900, "deviceScaleFactor": 1, "mobile": False}, b.session)

        def mount(demo=False):
            nonlocal visit
            visit += 1
            b.go(f"http://127.0.0.1:{port}/map-1g/?visit={visit}#chr1:1-49000")
            if not b.wait("!!(window.PGPMap && window.PGPUI)"):
                return False
            js("window.closedResult = true;document.documentElement.dataset.run='S'; "
                   "document.documentElement.dataset.mode = " + json.dumps("replay" if demo else "live") + "; "
                   "window.testMap = PGPMap.create({url:'/map-1g.json',"
                   "veil:function(id){return id==='closed' && window.closedResult ? 'closed':'open';},"
                   "known:function(){return true;},now:function(){return true;}});"
                   "document.body.appendChild(testMap.node);testMap.refresh('one');")
            return b.wait("!!document.querySelector('.gmap[data-level=reads]')")

        def loaded():
            return b.wait("(function(){var f=document.querySelector('.gm-reads iframe');"
                          "return f && f.contentDocument && f.contentDocument.documentElement.dataset.pgpStatus==='ready';})()", 30)

        check("map: opens at a 49 kb view", mount())
        check("map: under 50 kb offers the reads pill", b.eval("!!document.querySelector('.gm-reads-pill:not([hidden])')"))
        if not b.eval("!!document.querySelector('.gm-reads-pill')"):
            return
        b.eval("document.querySelector('.gm-reads-pill').click()")
        check("map: the reads panel loads the viewer in an iframe at the nearest open finding", loaded(),
              b.eval("document.querySelector('.gm-reads').textContent"))
        check("map: the iframe locus is 1-based, finding plus or minus 50",
              b.eval("document.querySelector('.gm-reads iframe').contentWindow.location.hash") == "#chr1:401-501")
        check("map: a closed result has no card or reads link",
              b.eval("!document.querySelector('.gm-marks [data-id=closed]') && !document.querySelector('.gm-card-reads[data-id=closed]')"))
        # A moved view must update the existing frame, with no page reload.
        b.eval("window.savedFrame=document.querySelector('.gm-reads iframe');location.hash='chr1:1-49001'")
        check("map: moving updates the same iframe to the map's locus",
              b.wait("savedFrame===document.querySelector('.gm-reads iframe') && savedFrame.contentWindow.location.hash==='#chr1:1-49001'"))
        check("viewer: alignments remain visible at 49 kb",
              b.wait("(function(){var h=savedFrame.contentDocument.getElementById('igv');"
                     "var v=h && h.shadowRoot && h.shadowRoot.querySelector('.igv-viewport[data-track-type=alignment]');"
                     "var n=v && v.querySelector('.igv-zoom-in-notice-container');"
                     "return !!(v && v.querySelector('canvas') && n && n.style.display!=='flex');})()"),
              b.eval("JSON.stringify(Array.from(savedFrame.contentDocument.getElementById('igv').shadowRoot.querySelectorAll('.igv-viewport')).map(function(v){"
                     "var n=v.querySelector('.igv-zoom-in-notice-container');return {type:v.dataset.trackType,"
                     "notice:n && n.style.display,canvas:v.querySelectorAll('canvas').length,text:v.textContent.slice(0,180)};}))"))
        if os.environ.get("MAP_SCREENS"):
            os.makedirs(os.environ["MAP_SCREENS"], exist_ok=True)
            for width, height in ((1440, 900), (390, 844)):
                b.call("Emulation.setDeviceMetricsOverride", {"width": width, "height": height,
                       "deviceScaleFactor": 1, "mobile": False}, b.session)
                b.eval("document.querySelector('.gm-reads').scrollIntoView({block:'center'})")
                png = b.call("Page.captureScreenshot", {"format": "png"}, b.session)["data"]
                with open(os.path.join(os.environ["MAP_SCREENS"], f"map-reads-{width}.png"), "wb") as out:
                    out.write(base64.b64decode(png))
            b.call("Emulation.setDeviceMetricsOverride", {"width": 1440, "height": 900,
                   "deviceScaleFactor": 1, "mobile": False}, b.session)
        b.eval("savedFrame.contentDocument.dispatchEvent(new KeyboardEvent('keydown',{key:'Escape',bubbles:true}))")
        check("map: Escape closes the reads and returns focus without moving the map",
              b.eval("document.querySelector('.gm-reads').hidden && document.activeElement===document.querySelector('.gm-canvas') && "
                     "document.querySelector('.gmap').dataset.view==='chr1:1-49001'"))
        b.eval("document.querySelector('.gm-canvas').dispatchEvent(new KeyboardEvent('keydown',{key:'r',bubbles:true}))")
        check("map: R reopens the panel", loaded())
        b.eval("document.querySelector('.gm-reads-close').click();document.querySelector('[data-stop=reads]').click()")
        check("map: Reads gauge opens a 2 kb view and the panel", b.wait("!document.querySelector('.gm-reads').hidden && "
                                                                 "document.querySelector('.gmap').dataset.level==='reads'"))
        js("document.querySelector('.gm-reads-close').click();location.hash='chr1:1-2000'")
        b.wait("document.querySelector('.gmap').dataset.view==='chr1:1-2000'")
        js("document.querySelector('.gm-marks [data-id=open]').focus();document.querySelector('.gm-card-reads').click()")
        check("map: the card opens reads at its finding", b.wait("document.querySelector('.gm-reads iframe').contentWindow.location.hash==='#chr1:401-501'"))
        js("window.closedResult=false;testMap.refresh('one');document.querySelector('.gm-marks [data-id=closed]').focus();"
               "document.querySelector('.gm-card-reads').click()")
        check("map: an opened finding offers reads", b.wait("document.querySelector('.gm-reads iframe').contentWindow.location.hash==='#chr1:501-601'"))
        b.eval("window.closedResult=true;testMap.refresh('one')")
        check("map: closing a result removes its open reads panel and link",
              b.eval("document.querySelector('.gm-reads').hidden && !document.querySelector('.gm-card-reads') && "
                     "!document.querySelector('.gm-reads a[href]') && "
                     "!document.querySelector('.gm-reads').textContent.includes('Closed finding')"))

        # No finding in view: use the view itself, and no viewer flags: show the server's reason.
        b.eval("location.hash='chr1:10000-59000'")
        b.wait("document.querySelector('.gmap').dataset.view==='chr1:10000-59000'")
        fasta = srv[0].cfg.viewer_fasta
        srv[0].cfg.viewer_fasta = None
        try:
            b.eval("document.querySelector('.gm-reads-pill').click()")
            check("map: a missing viewer shows the server's reason, without a broken frame",
                  b.wait("document.querySelector('.gm-reads').textContent.includes('not set up')")
                  and b.eval("!document.querySelector('.gm-reads iframe')"))
        finally:
            srv[0].cfg.viewer_fasta = fasta

        # CI proves the browser assertion detects the original policy mistake.
        policy["block"] = True
        mount()
        b.eval("window.blockedFrames=0;document.addEventListener('securitypolicyviolation',function(e){"
               "if(e.violatedDirective==='frame-src')window.blockedFrames++;});document.querySelector('.gm-reads-pill').click()")
        check("control: without frame-src the reads iframe is blocked", b.wait("window.blockedFrames>0") and b.eval(
              "(function(){var f=document.querySelector('.gm-reads iframe');return !f || !f.contentDocument || "
              "f.contentDocument.documentElement.dataset.pgpStatus!=='ready';})()"),
              b.eval("document.documentElement.dataset.pgpBlocked"))
        policy["block"] = False
        mount(demo=True)
        b.eval("document.querySelector('.gm-reads-pill').click()")
        check("map: the static demo explains why it has no viewer",
              b.wait("document.querySelector('.gm-reads').textContent.includes('The genome viewer is not in the demo')")
              and b.eval("!document.querySelector('.gm-reads iframe')"))
    finally:
        b.close()


def mutant_dir(work, name, old, new):
    """A plugin folder holding bin/owngenome_live_viewer.py with one planted mistake (beside a link
    to assets/, which the plugin finds from its own folder)."""
    top = os.path.join(work, f"mutant-{name}")
    d = os.path.join(top, "bin")
    os.makedirs(d)
    os.symlink(os.path.join(REPO, "assets"), os.path.join(top, "assets"))
    src = open(os.path.join(REPO, "bin", "owngenome_live_viewer.py")).read()
    if old not in src:
        check(f"mutant '{name}': its line is in owngenome_live_viewer.py", False, old)
    with open(os.path.join(d, "owngenome_live_viewer.py"), "w") as f:
        f.write(src.replace(old, new))
    return d


def main():
    work = tempfile.mkdtemp(prefix="pgp-live-viewer-")
    proc = None
    try:
        gd = f"{work}/genome"
        _, contigs = tv.make_sample(gd)
        # A 49 kb alignment check needs a reference longer than the original 3 kb sample.
        contigs[0] = ("chr1", contigs[0][1] + "ACGT" * 14250)
        fasta = f"{gd}/reference/synthetic.fasta"
        with open(fasta, "w") as fa, open(fasta + ".fai", "w") as fai:
            offset = 0
            for name, seq in contigs:
                offset += len(name) + 2
                fai.write(f"{name}\t{len(seq)}\t{offset}\t60\t61\n")
                fa.write(f">{name}\n" + "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)) + "\n")
                offset += len(seq) + (len(seq) + 59) // 60
        tv.write_bam(f"{gd}/S/aligned/S_sorted.bam", contigs,
                     [(0, p, contigs[0][1][p:p + 50]) for p in range(400, 700, 10)])
        s = f"{gd}/S"
        os.makedirs(f"{s}/nextflow")  # the site lists a folder with nextflow/ as a run
        html, _ = tv.findings_report(work)
        with open(f"{s}/S_report.html", "w") as f:
            f.write(html)
        state = f"{work}/state"
        os.makedirs(state)
        with open(f"{state}/password", "w") as f:
            f.write(owngenome_live.hash_password(PASSWORD) + "\n")

        # 1. owngenome-live.sh, with a password
        port, iport = tv.free_port(), tv.free_port()
        env = dict(os.environ, GENOME_DIR=gd, REF_FASTA="reference/synthetic.fasta", OWNGENOME_LIVE_PORT=str(port),
                   OWNGENOME_LIVE_INGEST_PORT=str(iport), OWNGENOME_LIVE_STATE_DIR=state)
        log = open(f"{work}/site.log", "w")
        proc = subprocess.Popen(["bash", f"{REPO}/scripts/owngenome-live.sh", "--quiet"], env=env,
                                stdout=log, stderr=subprocess.STDOUT)
        up = wait_port(port, proc)
        check("owngenome-live.sh: the site starts", up, open(f"{work}/site.log").read())
        if not up:
            raise SystemExit(1)
        pin = tv.igv_pin()
        igv = f"{gd}/tools/igv.js-{pin['IGV_JS_VERSION']}/igv.min.js"
        import hashlib
        check("owngenome-live.sh: the igv.js pinned in versions.env is downloaded, with its sha256",
              os.path.isfile(igv) and hashlib.sha256(open(igv, "rb").read()).hexdigest() == pin["IGV_JS_SHA256"],
              open(f"{work}/site.log").read())
        code, h, _ = request(port, "/runs/S/viewer/")
        c2, _, _ = request(port, "/runs/S/viewer/config.json")
        c3, _, _ = request(port, "/runs/S/viewer/data/2/S_sorted.bam", {"Range": "bytes=0-9"})
        check("no login: the page sends the browser to /login, the config and the reads are 401",
              code == 303 and h.get("location") == "/login?next=/runs/S/viewer/" and c2 == 401 and c3 == 401,
              (code, h.get("location"), c2, c3))
        auth = login(port)
        check("login: the password opens a session", auth is not None)
        auth = auth or {}

        # 2 and 3. the viewer and the served report
        acc = f"{state}/access.log"
        before = open(acc).read().count("\tviewer\t") if os.path.isfile(acc) else 0
        cfg = viewer_checks(port, auth, s)
        report_checks(port, auth)
        code, h, _ = request(port, "/runs/S/viewer", auth)
        check("/runs/S/viewer redirects to /runs/S/viewer/", code == 301 and h.get("location") == "/runs/S/viewer/",
              (code, h.get("location")))

        # 4. the access log
        time.sleep(0.5)
        rows = [r for r in open(acc).read().splitlines() if "\tviewer\t" in r]
        mine = rows[before:]
        check("access log: one row for the page, none for config.json, igv.js or the byte ranges",
              len(mine) == 2 and all("\t200" in r or "\t404" in r for r in mine), mine)
        proc.terminate()
        proc.wait(timeout=10)
        proc = None

        # 5. a BAM, and a VCF index, symlinked out of GENOME_DIR; an odd reference name
        t = f"{gd}/T"
        os.makedirs(f"{t}/aligned")
        os.makedirs(f"{t}/nextflow")
        shutil.copy(f"{s}/aligned/S_sorted.bam", f"{work}/outside.bam")
        shutil.copy(f"{s}/aligned/S_sorted.bam.bai", f"{work}/outside.bam.bai")
        os.symlink(f"{work}/outside.bam", f"{t}/aligned/T_sorted.bam")
        os.symlink(f"{work}/outside.bam.bai", f"{t}/aligned/T_sorted.bam.bai")
        u = f"{gd}/U"
        os.makedirs(f"{u}/vcf")
        os.makedirs(f"{u}/nextflow")
        shutil.copy(f"{s}/vcf/S.vcf.gz", f"{u}/vcf/U.vcf.gz")
        shutil.copy(f"{s}/vcf/S.vcf.gz.tbi", f"{work}/outside.vcf.gz.tbi")
        os.symlink(f"{work}/outside.vcf.gz.tbi", f"{u}/vcf/U.vcf.gz.tbi")
        odd = f"{gd}/reference/my ref#1?.fasta"
        shutil.copy(f"{gd}/reference/synthetic.fasta", odd)
        shutil.copy(f"{gd}/reference/synthetic.fasta.fai", odd + ".fai")
        flags = ["--root", gd, "--state-dir", f"{work}/state2", "--viewer-fasta", odd, "--viewer-igv-js", igv]
        srv = in_process(flags)
        try:
            p = srv[1].server_address[1]
            outside_checks(p, {})
            framing_checks(p)

            # 7. headless Chrome, on the site without a password (loopback)
            chrome = tv.find_chrome()
            if not chrome:
                ci = os.environ.get("CI") == "true"
                print(("[FAIL] " if ci else "[SKIP] ") + "no Chrome or Chromium to load the page in")
                if ci:
                    FAILED.append("no Chrome in CI")
            else:
                b = tl.Chrome(chrome, f"{work}/chrome")
                try:
                    b.go(f"http://127.0.0.1:{p}/runs/S/viewer/#chr1:451")
                    check("browser: the page from the site starts igv.js with the run's 3 tracks",
                          b.wait("document.documentElement.dataset.pgpStatus==='ready' && document.documentElement.dataset.pgpTracks==='3'"),
                          b.eval("document.documentElement.dataset.pgpStatus"))
                    check("browser: the policy blocked nothing the page needed", not b.eval("document.documentElement.dataset.pgpBlocked"))
                finally:
                    b.close()
                map_reads_checks(srv, work, chrome)
        finally:
            stop(*srv)

        # 6. a site without the viewer flags
        srv = in_process(["--root", gd, "--state-dir", f"{work}/state3"])
        try:
            p = srv[1].server_address[1]
            code, _, body = request(p, "/runs/S/viewer/")
            _, _, rep = request(p, "/runs/S/report.html")
            check("without the flags: the viewer is 404 with the reason, the report keeps its local links",
                  code == 404 and b"not set up" in body and render_report.VIEWER_URL.encode() in rep, (code, body))
        finally:
            stop(*srv)

        # Planted mistakes: checks 2, 3 and 5 must fail on them.
        for name, old, new in (
                ("serves any file", "        n = int(dm.group(1)) if dm else len(files)\n",
                 "        if rest and os.path.isfile(os.path.join(d, rest)):\n"
                 "            return serve.send_file(req, head_writer(req), 'application/octet-stream', "
                 "os.path.join(d, rest))\n        n = int(dm.group(1)) if dm else len(files)\n"),
                ("no report filter", "    app.report_filter(point_links)\n", ""),
                ("quoted name", "== urllib.parse.unquote(dm.group(2)):", "== dm.group(2):"),
                ("index out served", "    files = [f if n < 2 or inside(f) else None for n, f in enumerate(files)]\n",
                 "")):
            pd = mutant_dir(work, name.replace(" ", "-"), old, new)
            srv = in_process(flags[:3] + [f"{work}/state-{name.replace(' ', '-')}"] + flags[4:], pd)
            try:
                p = srv[1].server_address[1]
                quiet = []
                real, globals()["check"] = check, (lambda desc, ok, detail="": quiet.append(ok))
                try:
                    viewer_checks(p, {}, s)
                    report_checks(p, {})
                    outside_checks(p, {})
                finally:
                    globals()["check"] = real
                check(f"control: the checks fail on the planted mistake '{name}'", not all(quiet), quiet)
            finally:
                stop(*srv)
    finally:
        if proc and proc.poll() is None:
            proc.terminate()
            try:
                proc.wait(timeout=10)
            except subprocess.TimeoutExpired:
                proc.kill()
        shutil.rmtree(work, ignore_errors=True)

    print()
    print("RESULT: ALL PASS" if not FAILED else f"RESULT: {len(FAILED)} FAILED")
    sys.exit(1 if FAILED else 0)


if __name__ == "__main__":
    main()
