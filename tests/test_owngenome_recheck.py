#!/usr/bin/env python3
"""og recheck from the outside, the live-site plugin that shows its result, and the ClinVar lock
block of run-all.sh and step 06 that `og recheck --refresh-clinvar` interlocks with. No network:
the upstream answers come from tests/fixtures/recheck/offline.json.

The check itself (classes, steps, failed lookups, the locks, the rerun comparison, --launch-dir,
digests, unknown and skipped pins, steps 50 to 53) is tested by `cargo test -p og-recheck`
(rust/crates/og-recheck/src/tests.rs). Here:
  1. og recheck on the fixture sample: exit 0, ClinVar refresh-available with the steps 06 and
     downstream, state.json replaced whole and log.tsv grown by one row per item on a second check;
  2. an offline answer that failed makes og recheck exit 1, and a sample with no manifest exit 2;
  3. --rerun refuses (exit 3) while run_status.tsv has no finished_utc;
  4. --affected prints the steps of a key, and bin/owngenome_recheck.py, kept for crontabs, runs
     og recheck with the same arguments and says it is deprecated;
  5. the live-site plugin (bin/owngenome_live_recheck.py) lets the browser fetch
     recheck/state.json, and its state hook gives each step its badges ("newer ClinVar
     available", "redone <date>") and {} when there is no check;
  6. the ClinVar lock block of run-all.sh and of step 06 waits while a refresh holds the lock and
     goes on once it ends (needs flock; it fails in CI without it).
Each check is also run against a planted mistake and must fail there.

Run: python3 tests/test_owngenome_recheck.py (OG names the og binary, as tests/og_binary.py reads it)
"""
import fcntl
import gzip
import importlib.util
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "tests"))
from og_binary import og_binary  # noqa: E402

FIX = os.path.join(REPO, "tests", "fixtures", "recheck")
SHIM = os.path.join(REPO, "bin", "owngenome_recheck.py")
FAILS = 0
CLINVAR_STEPS = ["06", "13", "20", "23", "24", "27", "30", "31", "38", "39"]


def check(desc, ok, detail=""):
    global FAILS
    print(f"[{'PASS' if ok else 'FAIL'}] {desc}{'' if ok else ' -- ' + str(detail)[:800]}")
    if not ok:
        FAILS += 1


def genome(tmp, finished=True):
    g = os.path.join(tmp, "genome")
    s = os.path.join(g, "S1")
    os.makedirs(os.path.join(s, "logs"), exist_ok=True)
    shutil.copy(os.path.join(FIX, "run_manifest.tsv"), os.path.join(s, "run_manifest.tsv"))
    os.makedirs(os.path.join(g, "clinvar"), exist_ok=True)
    with gzip.open(os.path.join(g, "clinvar", "clinvar_pathogenic_chr.vcf.gz"), "wt") as f:
        f.write("##fileformat=VCFv4.1\n##fileDate=20260901\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    with open(os.path.join(s, "logs", "run_status.tsv"), "w") as f:
        f.write("meta\tstarted_epoch\t1000\nmeta\tstarted_utc\t2026-09-02 10:00:00 UTC\n"
                "meta\tcommand\techo replayed\nmeta\tswitches\tGENOME_DIR TOOLS\n")
        if finished:
            f.write("meta\tfinished_utc\t2026-09-02 12:00:00 UTC\n")
    return g


def og(cmd, *args):
    """cmd (og's path, or [interpreter, script]) with 'recheck' and args: (rc, stdout, stderr)."""
    argv = list(cmd) + (["recheck"] if len(cmd) == 1 else []) + list(args)
    p = subprocess.run(argv, capture_output=True, text=True, check=False, timeout=60)
    return p.returncode, p.stdout, p.stderr


def cli_checks(cmd):
    """Checks 1-3 on og recheck as cmd runs it: {name: (ok, detail)}."""
    out = {}
    tmp = tempfile.mkdtemp(prefix="recheck.")
    try:
        g = genome(tmp)
        off = os.path.join(FIX, "offline.json")
        base = ["--sample", "S1", "--genome-dir", g, "--offline", off]
        rc, o, e = og(cmd, *base)
        d = os.path.join(g, "S1", "recheck")
        st = json.load(open(os.path.join(d, "state.json"))) if os.path.isfile(os.path.join(d, "state.json")) else {}
        items = {i["item"]: i for i in st.get("items", [])}
        cv = items.get("clinvar_file_date", {})
        out["1. og recheck: ClinVar refresh-available, its steps 06 and downstream, exit 0"] = (
            rc == 0 and cv.get("class") == "refresh-available" and cv.get("steps") == CLINVAR_STEPS
            and st.get("affected_steps") == CLINVAR_STEPS, (rc, e, cv, st.get("affected_steps")))
        log1 = open(os.path.join(d, "log.tsv")).read() if os.path.isfile(os.path.join(d, "log.tsv")) else ""
        og(cmd, *base)
        log2 = open(os.path.join(d, "log.tsv")).read() if os.path.isfile(os.path.join(d, "log.tsv")) else ""
        left = sorted(os.listdir(d)) if os.path.isdir(d) else []
        out["1. state.json is replaced whole and log.tsv grows one row per item"] = (
            log2.startswith(log1) and len(log2.splitlines()) - len(log1.splitlines()) == len(items) > 0
            and left == ["log.tsv", "state.json"], (left, len(log1.splitlines()), len(log2.splitlines())))
        bad = os.path.join(tmp, "a500.json")
        data = json.load(open(off))
        data["PHARMCAT_IMAGE"] = {"error": "HTTP 500"}
        json.dump(data, open(bad, "w"))
        rc5, _, _ = og(cmd, "--sample", "S1", "--genome-dir", g, "--offline", bad)
        rc_none, _, e_none = og(cmd, "--sample", "S9", "--genome-dir", g, "--offline", off)
        out["2. a failed lookup exits 1, a sample with no manifest exits 2"] = (
            rc5 == 1 and rc_none == 2 and "run_manifest.tsv" in e_none, (rc5, rc_none, e_none))
        gr = genome(os.path.join(tmp, "running"), finished=False)
        rc3, _, e3 = og(cmd, "--sample", "S1", "--genome-dir", gr, "--offline", off, "--rerun")
        out["3. --rerun refuses while a run is going on"] = (rc3 == 3 and "going on" in e3, (rc3, e3))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    return out


def fake_og(tmp, body):
    """An executable that stands in for og: a planted mistake."""
    p = os.path.join(tmp, "og")
    with open(p, "w") as f:
        f.write("#!/bin/sh\n" + body + "\n")
    os.chmod(p, 0o755)
    return [p]


def lock_block(path):
    """The ClinVar lock block of a script: from 'CVLOCK=' through its closing 'fi'."""
    x = re.search(r"^CVLOCK=.*?\n  fi\nfi\n", open(path).read(), re.S | re.M)
    return x.group(0) if x else ""


def lock_waits(block):
    """block, run in bash while this process holds the lock exclusive: it must wait, then go on
    once the lock is let go."""
    import time
    tmp = tempfile.mkdtemp(prefix="recheck-cvlock.")
    try:
        g = os.path.join(tmp, "genome")
        os.makedirs(os.path.join(g, "clinvar"))
        lock = os.path.join(g, "clinvar", ".refresh.lock")
        with open(lock, "a") as held:
            fcntl.flock(held, fcntl.LOCK_EX | fcntl.LOCK_NB)
            p = subprocess.Popen(["bash", "-c", "set -euo pipefail\n" + block + "echo through\n"],
                                 env=dict(os.environ, G=g, GENOME_DIR=g), stdout=subprocess.PIPE, text=True)
            time.sleep(1.5)
            waited = p.poll() is None
        try:
            o, _ = p.communicate(timeout=20)
        except subprocess.TimeoutExpired:
            p.kill()
            o, _ = p.communicate()
        return (bool(block) and waited and p.returncode == 0 and "Waiting" in o and o.rstrip().endswith("through"),
                (waited, p.returncode, o))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def load_plugin():
    spec = importlib.util.spec_from_file_location("owngenome_live_recheck_test", os.path.join(REPO, "bin", "owngenome_live_recheck.py"))
    p = importlib.util.module_from_spec(spec)
    keep, sys.dont_write_bytecode = sys.dont_write_bytecode, True
    try:
        spec.loader.exec_module(p)
    finally:
        sys.dont_write_bytecode = keep
    return p


def plugin_ok(p):
    class App:
        def __init__(self):
            self.globs, self.hooks = [], {}

        def run_file(self, g):
            self.globs.append(g)

        def state_hook(self, name, fn):
            self.hooks[name] = fn
    app = App()
    p.register(app)
    tmp = tempfile.mkdtemp(prefix="recheck-live.")
    try:
        run = os.path.join(tmp, "S1")
        none = app.hooks["recheck"](run, {}, [])
        os.makedirs(os.path.join(run, "recheck"))
        json.dump({"checked_utc": "2026-10-06T00:00:00Z", "items": [
            {"item": "clinvar_file_date", "label": "ClinVar (step 06's file)", "class": "refresh-available", "steps": ["06", "24"]},
            {"item": "VEP_CACHE_RELEASE", "label": "Ensembl / VEP cache", "class": "no-newer-release", "steps": ["13"]}],
            "rerun": {"finished_utc": "2026-10-05T10:00:00Z", "predicted_vs_redone": {"redone": ["06"]}}},
            open(os.path.join(run, "recheck", "state.json"), "w"))
        got = app.hooks["recheck"](run, {}, [])
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    ok = (app.globs == ["recheck/state.json"] and none == {}
          and got.get("steps") == {"06": ["newer ClinVar available", "redone 2026-10-05"], "24": ["newer ClinVar available"]}
          and got.get("banner") == "Database re-check of 2026-10-06T00:00:00Z: newer ClinVar available")
    return ok, (app.globs, none, got)


def plugin_odd_rerun_ok(p):
    """A rerun block that is not an object: the banner still comes, no badge from it."""
    tmp = tempfile.mkdtemp(prefix="recheck-live-odd.")
    try:
        run = os.path.join(tmp, "S1")
        os.makedirs(os.path.join(run, "recheck"))
        out = []
        for rr in ([1], {"finished_utc": "2026-10-05T10:00:00Z", "predicted_vs_redone": ["06"]}):
            json.dump({"checked_utc": "2026-10-06T00:00:00Z", "items": [], "rerun": rr},
                      open(os.path.join(run, "recheck", "state.json"), "w"))
            try:
                out.append(p.recheck_state(run, {}, []))
            except (AttributeError, TypeError) as e:
                out.append(e)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    ok = all(isinstance(g, dict) and g.get("steps") == {} and "nothing newer" in g.get("banner", "") for g in out)
    return ok, out


def plugin_stale_ok(p):
    """A check older than the run's run_manifest.tsv describes an earlier run: no banner, no badge.
    One written after the manifest is still shown."""
    tmp = tempfile.mkdtemp(prefix="recheck-live-stale.")
    try:
        run = os.path.join(tmp, "S1")
        os.makedirs(os.path.join(run, "recheck"))
        json.dump({"checked_utc": "2026-10-06T00:00:00Z", "items": [
            {"item": "clinvar_file_date", "label": "ClinVar", "class": "refresh-available", "steps": ["06"]}]},
            open(os.path.join(run, "recheck", "state.json"), "w"))
        out = []
        for written in ("2026-10-07 09:00:00 UTC", "2026-10-05 09:00:00 UTC"):
            with open(os.path.join(run, "run_manifest.tsv"), "w") as f:
                f.write(f"run\twritten_utc\t{written}\n")
            out.append(p.recheck_state(run, {}, []))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    ok = out[0] == {} and out[1].get("steps") == {"06": ["newer ClinVar available"]}
    return ok, out


def plugin_from_source(src, name):
    p = os.path.join(REPO, "tmp") if os.path.isdir(os.path.join(REPO, "tmp")) else tempfile.gettempdir()
    path = os.path.join(tempfile.mkdtemp(prefix="recheck-plugin.", dir=p), name + ".py")
    open(path, "w").write(src)
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main():
    ogp = og_binary()
    results = cli_checks([ogp])
    for name, (ok, detail) in results.items():
        check(name, ok, detail)
    tmp = tempfile.mkdtemp(prefix="recheck-mutant.")
    try:
        # og that writes nothing and always exits 0: every check above must fail on it.
        silent = cli_checks(fake_og(tmp, "exit 0"))
        for name, (ok, _) in silent.items():
            check(f"fails on a planted mistake: {name} (an og that does nothing)", not ok, "passed on the mutant")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    rc, o, _ = og([ogp], "--affected", "clinvar_file_date")
    check("4. og recheck --affected names ClinVar's steps", rc == 0 and o.split() == CLINVAR_STEPS, (rc, o))
    # The shim runs the og beside it (bin/og). A copy of it beside an og that runs this og checks
    # it on any machine, with or without a build in bin/.og.
    tmp = tempfile.mkdtemp(prefix="recheck-shim.")
    try:
        shim = os.path.join(tmp, "owngenome_recheck.py")
        shutil.copy(SHIM, shim)
        with open(os.path.join(tmp, "og"), "w") as f:
            f.write(f"#!/bin/sh\nexec '{ogp}' \"$@\"\n")
        os.chmod(os.path.join(tmp, "og"), 0o755)
        rc, o, e = og([sys.executable, shim], "--affected", "clinvar_file_date")
        check("4. bin/owngenome_recheck.py runs og recheck and says it is deprecated",
              rc == 0 and o.split() == CLINVAR_STEPS and "deprecated" in e, (rc, o, e))
        rc, _, _ = og([sys.executable, shim], "--sample", "S1", "--genome-dir", "/nonexistent")
        check("4. bin/owngenome_recheck.py passes og recheck's exit code on", rc == 2, rc)
        with open(os.path.join(tmp, "og"), "w") as f:
            f.write("#!/bin/sh\nexit 0\n")
        rc, o, _ = og([sys.executable, shim], "--affected", "clinvar_file_date")
        check("fails on a planted mistake: bin/owngenome_recheck.py runs og recheck (an og that does nothing)",
              o.split() != CLINVAR_STEPS, "passed on the mutant")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    ok, detail = plugin_ok(load_plugin())
    check("5. the live-site plugin serves the check and gives each step its badges", ok, detail)
    bad = load_plugin()
    bad.ACT = dict(bad.ACT, **{"refresh-available": "newer {name}"})
    check("fails on a planted mistake: the live-site plugin serves the check and gives each step its badges",
          not plugin_ok(bad)[0], "passed on the mutant")
    ok, detail = plugin_odd_rerun_ok(load_plugin())
    check("5. the live-site plugin survives a rerun block that is not an object", ok, detail)
    src = open(os.path.join(REPO, "bin", "owngenome_live_recheck.py")).read()
    guard = "rr = rr if isinstance(rr, dict) else {}"
    check("planted mistake is in the copy: the live-site plugin survives a rerun block that is not an object", guard in src)
    unguarded = plugin_from_source(src.replace(guard, "rr = rr or {}"), "owngenome_live_recheck_unguarded")
    check("fails on a planted mistake: the live-site plugin survives a rerun block that is not an object",
          not plugin_odd_rerun_ok(unguarded)[0], "passed on the mutant")
    ok, detail = plugin_stale_ok(load_plugin())
    check("5. the live-site plugin leaves out a check older than the run's manifest", ok, detail)
    stale = "if utc14(manifest_written(run_dir)) > utc14(st.get(\"checked_utc\")):"
    check("planted mistake is in the copy: the live-site plugin leaves out a check older than the run's manifest",
          stale in src)
    unfiltered = plugin_from_source(src.replace(stale, "if False:"), "owngenome_live_recheck_unfiltered")
    check("fails on a planted mistake: the live-site plugin leaves out a check older than the run's manifest",
          not plugin_stale_ok(unfiltered)[0], "passed on the mutant")

    for script in ("scripts/run-all.sh", "scripts/06-clinvar-screen.sh"):
        block = lock_block(os.path.join(REPO, script))
        if not shutil.which("flock"):
            check(f"6. the ClinVar lock block of {script} waits for a refresh (flock is not installed: "
                  f"{'FAIL in CI' if os.environ.get('CI') else 'not run here'})", not os.environ.get("CI"))
            continue
        ok, detail = lock_waits(block)
        check(f"6. the ClinVar lock block of {script} waits for a refresh, then goes on", ok, detail)
        planted = re.sub(r"(?m)^    flock -s -n 7 .*$", "    :", block)
        check(f"planted mistake is in the copy: the ClinVar lock block of {script}", planted != block)
        check(f"fails on a planted mistake: the ClinVar lock block of {script} waits for a refresh (no flock call)",
              not lock_waits(planted)[0], "passed on the mutant")
    print("\nRESULT:", "ALL PASS" if FAILS == 0 else f"{FAILS} FAILED")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
