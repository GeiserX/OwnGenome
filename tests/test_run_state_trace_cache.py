#!/usr/bin/env python3
"""bin/run_state.py read_trace: a trace file is parsed once until it changes.

/runs.json builds the state of every run on the server, and each build read every
trace file of every run again. read_trace now keeps the parsed rows of each file,
keyed on its path and the file's identity (device, inode, mtime, ctime, size). Checks:
  1. the second read of an unchanged trace file opens nothing and gives the same rows;
  2. a file that grew (a task row appended) is read again and gives the new row;
  3. the rows handed out are copies: a caller that changes one does not change the cache;
  4. two samples of one trace file each get only their own rows from the one parse;
  5. the cache keeps at most TRACE_CACHE_MAX files;
  6. a trace replaced by another file of the same size and mtime (another run's
     file moved in, its times kept) is read again and gives the new rows;
  7. a duration cell that is no number ('1.2.3s', '.s') costs that row its end
     time only: the file still reads, every other row as before.
Check 1 is also run against a planted mistake (a cache that keeps nothing) and must
fail there, and check 6 against the old key (mtime and size alone): a check that
cannot fail is no check.

Run: python3 tests/test_run_state_trace_cache.py
"""
import builtins
import os
import shutil
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))
import run_state  # noqa: E402

FIX = os.path.join(REPO, "tests", "fixtures", "run_state", "trace.txt")
FAILS = []


def check(desc, ok, detail=""):
    print(("[PASS] " if ok else "[FAIL] ") + desc + ("" if ok else f" -- {str(detail)[:600]}"))
    if not ok:
        FAILS.append(desc)


def opens_of(path, fn):
    """(result of fn(), how many times fn opened path)."""
    n, real = [0], builtins.open

    def counting(file, *a, **k):
        if os.path.realpath(str(file)) == os.path.realpath(path):
            n[0] += 1
        return real(file, *a, **k)
    builtins.open = counting
    try:
        out = fn()
    finally:
        builtins.open = real
    return out, n[0]


def second_read_opens(trace):
    run_state.TRACE_CACHE.clear()
    first = run_state.read_trace([trace], "demo1")
    again, n = opens_of(trace, lambda: run_state.read_trace([trace], "demo1"))
    return first, again, n


def main():
    work = tempfile.mkdtemp(prefix="test_trace_cache.")
    try:
        trace = os.path.join(work, "trace.txt")
        shutil.copy(FIX, trace)

        # 1. an unchanged file is not opened again
        first, again, n = second_read_opens(trace)
        check("1. the second read of an unchanged trace file opens nothing and gives the same rows",
              n == 0 and first == again and len(first) > 0, (n, len(first), len(again)))

        # planted mistake: a cache that keeps nothing
        keep = run_state.TRACE_CACHE_MAX
        run_state.TRACE_CACHE_MAX = 0
        try:
            _, _, n0 = second_read_opens(trace)
        finally:
            run_state.TRACE_CACHE_MAX = keep
        check("control: a cache that keeps nothing fails check 1", n0 > 0, n0)

        # 2. a file that grew is read again
        run_state.TRACE_CACHE.clear()
        before = run_state.read_trace([trace], "demo1")
        with open(trace, "a") as f:
            f.write("99\tzz/000099\t9999\tPGX:CLINVAR_SCREEN (demo1)\tCOMPLETED\t0\t2026-01-02 11:00:00.000\t5s\t"
                    "4s\t50.0%\t1 MB\t1 MB\t1 MB\t1 MB\n")
        after, n = opens_of(trace, lambda: run_state.read_trace([trace], "demo1"))
        check("2. a trace file that grew is read again and gives the new row",
              n == 1 and len(after) == len(before) + 1 and after[-1]["id"] == "99", (n, len(before), len(after)))

        # 3. copies
        rows = run_state.read_trace([trace], "demo1")
        rows[0]["status"] = "CHANGED BY A CALLER"
        rows2 = run_state.read_trace([trace], "demo1")
        check("3. a caller changing a row does not change the next read", rows2[0]["status"] != "CHANGED BY A CALLER",
              rows2[0])

        # 4. two samples from one parse
        with open(trace, "a") as f:
            f.write("100\tzz/000100\t10000\tPGX:CLINVAR_SCREEN (demo2)\tCOMPLETED\t0\t2026-01-02 11:00:00.000\t5s\t"
                    "4s\t50.0%\t1 MB\t1 MB\t1 MB\t1 MB\n")
        run_state.TRACE_CACHE.clear()
        one = run_state.read_trace([trace], "demo1")
        two, n = opens_of(trace, lambda: run_state.read_trace([trace], "demo2"))
        check("4. a second sample of the same trace file gets its own rows from the one parse",
              n == 0 and [r["id"] for r in two] == ["100"] and "100" not in [r["id"] for r in one], (n, two))

        # 5. bounded
        run_state.TRACE_CACHE.clear()
        for i in range(run_state.TRACE_CACHE_MAX + 5):
            p = os.path.join(work, f"t{i}.txt")
            shutil.copy(FIX, p)
            run_state.read_trace([p], "demo1")
        check("5. the cache keeps at most TRACE_CACHE_MAX files", len(run_state.TRACE_CACHE) == run_state.TRACE_CACHE_MAX,
              len(run_state.TRACE_CACHE))

        # 6. a file replaced by one of the same size and mtime
        def replaced(path):
            with open(FIX) as f:
                text = f.read()
            with open(path, "w") as f:
                f.write(text)
            run_state.TRACE_CACHE.clear()
            old = run_state.read_trace([path], "demo1")
            st = os.stat(path)
            new = path + ".new"
            with open(new, "w") as f:  # the same length: one task's status COMPLETED becomes FAILED...
                f.write(text.replace("\tCOMPLETED\t0\t", "\tFAILED...\t1\t", 1))
            os.utime(new, ns=(st.st_atime_ns, st.st_mtime_ns))
            os.replace(new, path)
            st2 = os.stat(path)
            got = run_state.read_trace([path], "demo1")
            return old, got, (st.st_size, st.st_mtime_ns) == (st2.st_size, st2.st_mtime_ns)
        p6 = os.path.join(work, "replaced.txt")
        old, got, same = replaced(p6)
        check("6. a trace replaced by a file of the same size and mtime is read again",
              same and [r["status"] for r in got] != [r["status"] for r in old] and "FAILED..." in [r["status"] for r in got],
              (same, [r["status"] for r in old], [r["status"] for r in got]))
        real_stat = os.stat

        class OldKey:  # the key before: mtime and size alone (os.stat is run_state's os.stat too)
            def __init__(self, st):
                self.st = st
                self.st_dev = self.st_ino = self.st_ctime_ns = 0

            def __getattr__(self, k):
                return getattr(self.st, k)
        run_state.os.stat = lambda p, *a, **k: OldKey(real_stat(p, *a, **k))
        try:
            old, got, same = replaced(p6)
        finally:
            run_state.os.stat = real_stat
        check("control: a key of mtime and size alone fails check 6",
              same and [r["status"] for r in got] == [r["status"] for r in old], (same, [r["status"] for r in got]))

        # 7. a malformed duration
        with open(FIX) as f:
            text = f.read()
        lines = text.splitlines(True)
        bad = lines[2].split("\t")
        bad[7] = "1.2.3s"
        lines[2] = "\t".join(bad)
        lines[1] = lines[1].replace("\t2s\t", "\t.s\t", 1)
        p7 = os.path.join(work, "durations.txt")
        with open(p7, "w") as f:
            f.write("".join(lines))
        run_state.TRACE_CACHE.clear()
        want = run_state.read_trace([FIX], "demo1")
        try:
            got = run_state.read_trace([p7], "demo1")
        except ValueError as e:
            got = e
        ok = isinstance(got, list) and len(got) == len(want) and got[0]["end"] is None and got[1]["end"] is None \
            and got[2:] == want[2:] and got[1]["start"] == want[1]["start"]
        check("7. a duration that is no number costs its row the end time only", ok, got)
    finally:
        shutil.rmtree(work, ignore_errors=True)
    if FAILS:
        print(f"\n{len(FAILS)} check(s) failed")
        return 1
    print("\nAll checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
