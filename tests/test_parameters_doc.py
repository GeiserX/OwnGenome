#!/usr/bin/env python3
"""docs/parameters.md lists every parameter of nextflow_schema.json.

Checks:
  1. every parameter of the schema (each group's properties, and any
     top-level ones) has a row "| `--<name>` |" in docs/parameters.md;
  2. every row of the page names a schema parameter, or a parameter main.nf
     reads that the schema does not hold yet (DOC_ONLY, each checked against
     main.nf), so a removed or renamed parameter cannot keep its row;
  3. controls: the page with one row taken out fails check 1 and names that
     parameter; the page with a planted row fails check 2.
"""
import json
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCHEMA = os.path.join(ROOT, "nextflow_schema.json")
PAGE = os.path.join(ROOT, "docs", "parameters.md")
MAIN = os.path.join(ROOT, "main.nf")
# Read by main.nf, not (yet) in the schema; the page writes their rows by hand.
DOC_ONLY = {"pedigree"}
ROW = re.compile(r"^\| `--([A-Za-z0-9_]+)` \|", re.M)

failures = 0


def check(desc, ok, detail=""):
    global failures
    print(f"[{'PASS' if ok else 'FAIL'}] {desc}" + ("" if ok else f" -- {detail}"))
    if not ok:
        failures += 1


def schema_params(schema):
    names = set(schema.get("properties", {}))
    for group in schema.get("definitions", schema.get("$defs", {})).values():
        names |= set(group.get("properties", {}))
    return names


def problems(params, page_text, main_text):
    rows = ROW.findall(page_text)
    missing = sorted(params - set(rows))
    extra = sorted(r for r in set(rows) - params
                   if r not in DOC_ONLY or not re.search(rf"\bparams\.{r}\b|'{r}'", main_text))
    return missing, extra


with open(SCHEMA, encoding="utf-8") as f:
    params = schema_params(json.load(f))
with open(PAGE, encoding="utf-8") as f:
    page = f.read()
with open(MAIN, encoding="utf-8") as f:
    main_nf = f.read()

check("the schema has parameters to compare", len(params) > 20, len(params))
missing, extra = problems(params, page, main_nf)
check("every schema parameter has a row in docs/parameters.md", not missing, f"no row for: {', '.join(missing)}")
check("every row of docs/parameters.md is a schema parameter or one main.nf reads", not extra,
      f"not in the schema: {', '.join(extra)}")
check("each parameter has one row", len(ROW.findall(page)) == len(set(ROW.findall(page))),
      sorted(r for r in set(ROW.findall(page)) if ROW.findall(page).count(r) > 1))

# Controls: the same comparison must catch a dropped row and a planted one.
victim = sorted(params)[0]
dropped = re.sub(rf"^\| `--{victim}` \|.*\n", "", page, flags=re.M)
check(f"control: the page without the --{victim} row fails, naming it",
      problems(params, dropped, main_nf)[0] == [victim], problems(params, dropped, main_nf))
planted = page + "| `--no_such_parameter` | none | planted |\n"
check("control: a planted row for a parameter the schema lacks fails",
      problems(params, planted, main_nf)[1] == ["no_such_parameter"], problems(params, planted, main_nf))

print(f"\nRESULT: {'ALL PASS' if not failures else f'{failures} FAILED'}")
sys.exit(1 if failures else 0)
