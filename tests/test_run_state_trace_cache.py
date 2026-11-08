#!/usr/bin/env python3
"""bin/run_state.py read_trace: a trace file is parsed once until it changes.

/runs.json builds the state of every run on the server, and each build read every
trace file of every run again. read_trace now keeps the parsed rows of each file,
keyed on its path, mtime and size. Checks:
  1. the second read of an unchanged trace file opens nothing and gives the same rows;
  2. a file that grew (a task row appended) is read again and gives the new row;
  3. the rows handed out are copies: a caller that changes one does not change the cache;
  4. two samples of one trace file each get only their own rows from the one parse;
  5. the cache keeps at most TRACE_CACHE_MAX files.
Check 1 is also run against a planted mistake (a cache that keeps nothing) and must
fail there: a check that cannot fail is no check.

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
    finally:
        shutil.rmtree(work, ignore_errors=True)
    if FAILS:
        print(f"\n{len(FAILS)} check(s) failed")
        return 1
    print("\nAll checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
