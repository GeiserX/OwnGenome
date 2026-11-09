#!/usr/bin/env python3
"""scripts/ci/check-report-only.py: report-only processes carry the ignore strategy.

  1. the check passes on this checkout and names the two cDNA read checks;
  2. on a copy whose report-only list drops CLINVAR_CDNA_CHECK (the state
     before it was added) it fails and names that process; so it does when
     CPSR_CDNA_CHECK is dropped, or MITO_DELETIONS's own block is removed;
  3. on a copy where a second consumer reads CLINVAR_CDNA_CHECK's channel, the
     process feeds another step, so the check no longer asks for it;
  4. a KEEP_STOP name the scan no longer finds is refused: on a copy where
     ROH's channel gets a second consumer, the check fails naming ROH.

Run: python3 tests/test_contract_report_only.py
"""
import os
import re
import shutil
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CHECK = os.path.join(REPO, "scripts", "ci", "check-report-only.py")
FAILS = 0


def check(desc, ok, detail=""):
    global FAILS
    print("%s: %s" % ("ok" if ok else "FAIL", desc))
    if not ok:
        FAILS += 1
        if detail:
            print(detail)


def run(repo):
    p = subprocess.run([sys.executable, "-I", CHECK, "--repo", repo],
                       capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr


def copy():
    d = tempfile.mkdtemp(prefix="report-only-")
    shutil.copy(os.path.join(REPO, "main.nf"), d)
    shutil.copy(os.path.join(REPO, "nextflow.config"), d)
    shutil.copytree(os.path.join(REPO, "workflows"), os.path.join(d, "workflows"))
    return d


def edit(path, old, new, count=1):
    text = open(path).read()
    assert old in text, "planted edit not found in %s: %r" % (path, old)
    open(path, "w").write(text.replace(old, new, count))


rc, out = run(REPO)
check("checkout passes", rc == 0, out)
check("checkout names both cDNA checks",
      "ok: CLINVAR_CDNA_CHECK" in out and "ok: CPSR_CDNA_CHECK" in out, out)

for name, plant in [
    ("CLINVAR_CDNA_CHECK", lambda d: edit(os.path.join(d, "nextflow.config"),
                                          "|CLINVAR_CDNA_CHECK|", "|")),
    ("CPSR_CDNA_CHECK", lambda d: edit(os.path.join(d, "nextflow.config"),
                                       "|CPSR_CDNA_CHECK|", "|")),
    ("MITO_DELETIONS", lambda d: edit(os.path.join(d, "nextflow.config"),
                                      "withName: 'MITO_DELETIONS'", "withName: 'NO_SUCH_PROCESS'")),
]:
    d = copy()
    plant(d)
    rc, out = run(d)
    check("dropping %s from the report-only list fails and names it" % name,
          rc == 1 and re.search(r"ERROR: %s \(workflows/" % name, out) is not None, out)
    shutil.rmtree(d)

d = copy()
pgx = os.path.join(d, "workflows", "pgx.nf")
edit(os.path.join(d, "nextflow.config"), "|CLINVAR_CDNA_CHECK|", "|")
edit(pgx, "    clinvar_cdna         = ch_clinvar_cdna",
     "    clinvar_cdna_more    = ch_clinvar_cdna.map { it }\n    clinvar_cdna         = ch_clinvar_cdna")
rc, out = run(d)
check("a second consumer makes CLINVAR_CDNA_CHECK a feeder, not asked for",
      rc == 0 and "CLINVAR_CDNA_CHECK" not in out, out)
shutil.rmtree(d)

d = copy()
clin = os.path.join(d, "workflows", "clinical.nf")
text = open(clin).read()
m = re.search(r"^\s*roh_regions\s*=\s*(ch_\w+)", text, re.M)
check("ROH has its emit in workflows/clinical.nf", m is not None)
if m:
    edit(clin, m.group(0), "    roh_more = %s.map { it }\n%s" % (m.group(1), m.group(0)))
    rc, out = run(d)
    check("a KEEP_STOP name the scan no longer finds fails",
          rc == 1 and "KEEP_STOP names ROH" in out, out)
shutil.rmtree(d)

print("FAILED: %d" % FAILS if FAILS else "all report-only checks passed")
sys.exit(1 if FAILS else 0)
