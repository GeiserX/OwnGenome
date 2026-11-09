#!/usr/bin/env python3
"""A run that stopped writing reads 'stopped', and run-all.sh prints the run's page address.

Checks (synthetic run folders, no genome data):
  1. a run with no end record whose files were last written three days ago
     reads 'stopped', with last_heard_utc the newest of those files;
  2. the same run written a minute ago reads 'running';
  3. a run quiet for longer than 6 h but shorter than twice the longest step's
     usual hours still reads 'running' (a long DeepVariant task posts nothing);
  4. a stopped run that gets a new message reads 'running' again;
  5. a run with an end record (meta finished_utc) reads 'finished', however old;
  6. the live site's /runs.json gives the old run 'stopped' with its
     last_heard_utc, and the fresh one 'running', though both write their
     trace under the same pipeline_info (another sample's run going on does
     not keep a dead one alive);
  7. run-all.sh's live block prints '<OWNGENOME_LIVE_PAGE_URL>/runs/<sample>/'
     when the setting is given, and no page line without it;
  8. collect() without stale=True (what the rerun and ClinVar refresh guards
     read) still calls the three-day-quiet run 'running'.
Checks 1, 4, 7 and 8 are then run against a planted mistake (no stale limit;
the page line removed; stale on by default) and must fail there: a check that
cannot fail is no check.

Run: python3 tests/test_run_state_stopped.py
"""
import importlib.util
import json
import os
import re
import shutil
import socket
import subprocess
import sys
import tempfile
import time
import urllib.request

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))
FAILS = 0
DAY = 86400


def check(desc, ok, detail=""):
    global FAILS
    print(f"[{'PASS' if ok else 'FAIL'}] {desc}{'' if ok else ' -- ' + str(detail)[:600]}")
    if not ok:
        FAILS += 1


def load(src=None, name="run_state"):
    spec = importlib.util.spec_from_file_location(name, src or os.path.join(REPO, "bin", "run_state.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def make_run(top, sample, age_s, finished=False):
    """A run folder under top: run_status.tsv with two planned steps, one trace row, all aged age_s."""
    d = os.path.join(top, sample)
    os.makedirs(os.path.join(d, "logs"), exist_ok=True)
    os.makedirs(os.path.join(top, "pipeline_info"), exist_ok=True)
    t = time.time() - age_s
    started = t - 3600
    with open(os.path.join(d, "logs", "run_status.tsv"), "w") as f:
        f.write(f"meta\tstarted_epoch\t{started:.0f}\nstep\t16\tplanned\nstep\t03\tplanned\n")
        if finished:
            f.write("meta\tfinished_utc\t" + time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime(t)) + "\n")
    trace = os.path.join(top, "pipeline_info", f"trace_{sample}.txt")
    ms = int((t - 600) * 1000)
    with open(trace, "w") as f:
        f.write("task_id\thash\tnative_id\tname\tstatus\texit\tsubmit\tstart\tcomplete\n")
        f.write(f"1\tab/cdef12\t1\tFASTP ({sample})\tCOMPLETED\t0\t{ms}\t{ms}\t{ms + 300000}\n")
    for p in (os.path.join(d, "logs", "run_status.tsv"), trace):
        os.utime(p, (t, t))
    return d


def state_checks(rs, work):
    out = {}
    old = make_run(os.path.join(work, "a"), "old1", 3 * DAY)
    guard = rs.collect("old1", old)
    out["collect() without stale=True (the rerun and refresh guards) still reads running"] = (
        guard["run"]["state"] == "running" and "last_heard_utc" not in guard["run"], guard["run"])
    st = rs.collect("old1", old, stale=True)
    heard = st["run"].get("last_heard_utc")
    want = time.strftime("%Y-%m-%d", time.gmtime(time.time() - 3 * DAY))
    out["a run silent for three days with no end record reads stopped"] = (
        st["run"]["state"] == "stopped" and bool(heard) and heard.startswith(want), st["run"])

    fresh = make_run(os.path.join(work, "b"), "new1", 60)
    st = rs.collect("new1", fresh, stale=True)
    out["a run written a minute ago reads running"] = (st["run"]["state"] == "running", st["run"])

    mid = make_run(os.path.join(work, "c"), "mid1", 8 * 3600)
    st = rs.collect("mid1", mid, stale=True)
    out["a run quiet for 8 h (under twice the longest step) still reads running"] = (
        st["run"]["state"] == "running", st["run"])

    ev = os.path.join(work, "old1.events.jsonl")
    with open(ev, "w") as f:
        f.write(json.dumps({"runName": "old1", "event": "process_started",
                            "utcTime": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                            "trace": {"task_id": 2, "name": "BWA_MEM2 (old1)", "status": "RUNNING",
                                      "start": int(time.time() * 1000)}}) + "\n")
    st = rs.collect("old1", old, events_path=ev, stale=True)
    out["a stopped run that gets a new message reads running again"] = (st["run"]["state"] == "running", st["run"])

    done = make_run(os.path.join(work, "d"), "end1", 30 * DAY, finished=True)
    st = rs.collect("end1", done, stale=True)
    out["a run with an end record reads finished, however old"] = (st["run"]["state"] == "finished", st["run"])
    return out


def free_port():
    s = socket.socket()
    s.bind(("127.0.0.1", 0))
    p = s.getsockname()[1]
    s.close()
    return p


def site_checks(work):
    root = os.path.join(work, "site")
    make_run(root, "old1", 3 * DAY)
    make_run(root, "new1", 60)
    web, ing = free_port(), free_port()
    p = subprocess.Popen([sys.executable, os.path.join(REPO, "bin", "owngenome_live.py"), "--root", root,
                          "--state-dir", os.path.join(work, "state"), "--port", str(web), "--ingest-port", str(ing),
                          "--quiet"], stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    data, err = None, ""
    try:
        for _ in range(100):
            try:
                with urllib.request.urlopen(f"http://127.0.0.1:{web}/runs.json", timeout=5) as r:
                    data = json.load(r)
                break
            except OSError:
                time.sleep(0.1)
    finally:
        p.terminate()
        try:
            err = p.communicate(timeout=10)[1].decode(errors="replace")
        except subprocess.TimeoutExpired:
            p.kill()
    runs = {r["name"]: r for r in (data or {}).get("runs", [])}
    o, n = runs.get("old1", {}), runs.get("new1", {})
    check("the live site's /runs.json gives the old run stopped with when it was last heard from",
          o.get("state") == "stopped" and bool(o.get("last_heard_utc")), o or err)
    check("the live site's /runs.json gives the fresh run running", n.get("state") == "running", n or err)


def live_block(src):
    m = re.search(r'^if \[ -n "\$\{OWNGENOME_LIVE_URL:-\}" \]; then\n.*?^fi$', src, re.S | re.M)
    return m.group(0) if m else ""


def page_checks(block, work):
    out = {}

    def run(env):
        s = os.path.join(work, "S")
        os.makedirs(os.path.join(s, "nextflow"), exist_ok=True)
        script = f'S={s}\nSAMPLE=demo1\nNF_WEBLOG_PLUGIN=nf-weblog@0\nNF=()\n{block}\n'
        e = {k: v for k, v in os.environ.items() if not k.startswith("OWNGENOME_LIVE")}
        e.update(env)
        return subprocess.run(["bash", "-c", script], capture_output=True, text=True, env=e, timeout=30)
    r = run({"OWNGENOME_LIVE_URL": "http://127.0.0.1:8788", "OWNGENOME_LIVE_PAGE_URL": "https://genome.example.org/"})
    out["with OWNGENOME_LIVE_PAGE_URL, run-all.sh prints the run's page"] = (
        r.returncode == 0 and "https://genome.example.org/runs/demo1/\n" in r.stdout and "org//runs" not in r.stdout,
        r.stdout + r.stderr)
    r = run({"OWNGENOME_LIVE_URL": "http://127.0.0.1:8788"})
    out["without it, run-all.sh prints no page line"] = (
        r.returncode == 0 and "/runs/" not in r.stdout and "ingest/demo1" in r.stdout, r.stdout + r.stderr)
    return out


def main():
    work = tempfile.mkdtemp(prefix="pgp-stopped-")
    try:
        rs = load()
        for name, (ok, detail) in state_checks(rs, os.path.join(work, "good")).items():
            check(name, ok, detail)
        site_checks(work)
        block = live_block(open(os.path.join(REPO, "scripts", "run-all.sh")).read())
        check("run-all.sh has its live block", bool(block))
        for name, (ok, detail) in page_checks(block, work).items():
            check(name, ok, detail)

        # Planted mistakes: the checks above must fail on each.
        bad = load(name="run_state_mutant")
        bad.stale_after = lambda steps: float("inf")  # the state before this fix: never stopped
        got = state_checks(bad, os.path.join(work, "mutant"))
        check("fails on a planted mistake: a run with no stale limit is called running",
              not got["a run silent for three days with no end record reads stopped"][0], "passed on the mutant")
        bad.stale_after = lambda steps: 0.0  # everything stopped at once
        got = state_checks(bad, os.path.join(work, "mutant0"))
        check("fails on a planted mistake: a run that gets a message stays stopped",
              not got["a stopped run that gets a new message reads running again"][0], "passed on the mutant")
        eager = load(name="run_state_mutant_eager")
        plain = eager.collect
        eager.collect = lambda *a, **k: plain(*a, **{**k, "stale": True})  # stale on for every caller
        got = state_checks(eager, os.path.join(work, "mutant-eager"))
        check("fails on a planted mistake: collect() calls a quiet run stopped for every caller",
              not got["collect() without stale=True (the rerun and refresh guards) still reads running"][0],
              "passed on the mutant")
        cut = "\n".join(line for line in block.splitlines() if "OWNGENOME_LIVE_PAGE_URL" not in line)
        got = page_checks(cut, os.path.join(work, "mutant-sh"))
        check("fails on a planted mistake: run-all.sh without the page line",
              not got["with OWNGENOME_LIVE_PAGE_URL, run-all.sh prints the run's page"][0], "passed on the mutant")
    finally:
        shutil.rmtree(work, ignore_errors=True)
    print("\nRESULT:", "ALL PASS" if FAILS == 0 else f"{FAILS} FAILED")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
