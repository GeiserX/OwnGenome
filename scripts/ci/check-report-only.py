#!/usr/bin/env python3
"""Every process whose output only the report reads has the report-only error strategy.

nextflow.config gives report-only processes a withName errorStrategy that
retries a task killed for memory or time once and otherwise ignores the
failure, so one odd input cannot stop the run before DeepVariant and the
report. A process left out keeps conf/base.config's stop. This check finds
the processes whose output only the report reads and fails for each one no
such block names.

A process counts as report-only when, in its workflows/<name>.nf, every
reference to <PROCESS>.out.<channel> is either `ch_versions =
ch_versions.mix(<PROCESS>.out.versions)` or a plain `ch_x = <PROCESS>.out.<y>`,
that ch_x is used nowhere else in the file but in its own assignments and in
the emit block, and main.nf uses every such emit only in a line
`.join(<WORKFLOW>.out.<emit>... remainder: true)`, the report's joins that keep
a sample whose step is missing. Anything else (a map on the channel, a second
consumer, a join without remainder) makes the process a feeder of something
else, and it is not asked for. The check errs toward silence: it never asks for
a process that feeds another step.

KEEP_STOP names the report-only processes that keep the stop by choice. A name
there that the scan no longer finds as report-only is an error too, so the
list cannot go stale.

Usage: check-report-only.py [--repo DIR]
Exit 0 when every report-only process is covered, 1 when one is not, 2 on a
usage error. Only the standard library.
"""
import argparse
import os
import re
import sys

# Report-only by their channels, kept fatal until a decision says otherwise:
# ROH (step 17's runs of homozygosity) and CPIC_LOOKUP (the CPIC
# recommendations) have not been through the report-only review yet.
KEEP_STOP = {"ROH", "CPIC_LOOKUP"}

IGNORE_BLOCK = re.compile(
    r"withName:\s*'([^']*)'\s*\{\s*errorStrategy\s*=\s*\{[^}]*'ignore'[^}]*\}", re.S)


def strip_comments(text):
    return "\n".join(line.split("//")[0] for line in text.splitlines())


def report_only_patterns(config_text):
    return [re.compile("^(?:%s)$" % m.group(1))
            for m in IGNORE_BLOCK.finditer(strip_comments(config_text))]


def report_only_processes(repo):
    """Yield (workflow file, process, [emit, ...]) for each report-only process."""
    main = strip_comments(open(os.path.join(repo, "main.nf")).read())
    wdir = os.path.join(repo, "workflows")
    for fname in sorted(os.listdir(wdir)):
        if not fname.endswith(".nf"):
            continue
        text = strip_comments(open(os.path.join(wdir, fname)).read())
        wm = re.search(r"^workflow\s+(\w+)", text, re.M)
        if not wm or "emit:" not in text:
            continue
        wf = wm.group(1)
        emits = re.findall(r"^\s*(\w+)\s*=\s*(ch_\w+)\s*$", text.split("emit:")[-1], re.M)
        for proc in sorted(set(re.findall(r"\b([A-Z][A-Z0-9_]+)\.out\.", text))):
            used = []
            ok = True
            for line in re.findall(r"^.*\b%s\.out\.\w+.*$" % proc, text, re.M):
                if re.match(r"^\s*ch_versions\s*=\s*ch_versions\.mix\(\s*%s\.out\.versions\s*\)\s*$"
                            % proc, line):
                    continue
                a = re.match(r"^\s*(ch_\w+)\s*=\s*%s\.out\.\w+\s*$" % proc, line)
                if not a:
                    ok = False
                    break
                ch = a.group(1)
                others = [l for l in re.findall(r"^.*\b%s\b.*$" % ch, text, re.M)
                          if not re.match(r"^\s*(%s\s*=|\w+\s*=\s*%s\s*$)" % (ch, ch), l)]
                names = [e for e, c in emits if c == ch]
                if others or not names:
                    ok = False
                    break
                for emit in names:
                    lines = re.findall(r"^.*\b%s\.out\.%s\b.*$" % (wf, emit), main, re.M)
                    if not lines or not all(
                            re.search(r"\.join\(\s*%s\.out\.%s\b.*remainder:\s*true" % (wf, emit), l)
                            for l in lines):
                        ok = False
                    used.append(emit)
            if ok and used:
                yield fname, proc, used


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--repo", default=os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__)))))
    args = ap.parse_args()
    try:
        patterns = report_only_patterns(open(os.path.join(args.repo, "nextflow.config")).read())
    except OSError as e:
        print("ERROR: %s" % e, file=sys.stderr)
        return 2
    if not patterns:
        print("ERROR: no report-only withName block in nextflow.config", file=sys.stderr)
        return 1
    bad = 0
    found = set()
    for fname, proc, emits in report_only_processes(args.repo):
        found.add(proc)
        covered = any(p.match(proc) for p in patterns)
        if proc in KEEP_STOP:
            print("keep-stop: %s (workflows/%s), by choice" % (proc, fname))
        elif covered:
            print("ok: %s (workflows/%s)" % (proc, fname))
        else:
            print("ERROR: %s (workflows/%s) feeds only the report (%s, joined with remainder: true "
                  "in main.nf) but no report-only withName block in nextflow.config names it, so "
                  "one failed task stops the whole run" % (proc, fname, ", ".join(emits)),
                  file=sys.stderr)
            bad += 1
    for proc in sorted(KEEP_STOP - found):
        print("ERROR: KEEP_STOP names %s, which the scan no longer finds as report-only; "
              "drop it from the list" % proc, file=sys.stderr)
        bad += 1
    if not found:
        print("ERROR: the scan found no report-only process at all; the patterns no longer "
              "match the code", file=sys.stderr)
        bad += 1
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
